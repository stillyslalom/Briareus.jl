# The drive, optimization and clearance analyses through a DesignSession, inline (as a
# script would run them) and through a window-style delivery queue. A coarse hohlraum keeps
# the view factors to seconds.

const SESS = DesignSession(resolution = :coarse)

@testset "drive (no GL)" begin
    sess = SESS
    @test task_status(sess, :drive) == (:todo, "not computed")
    @test task_status(sess, :optimize)[1] === :todo
    @test task_status(sess, :clearance)[1] === :todo
    @test task_status(sess, :talis)[1] === :todo
    @test_throws ArgumentError task_status(sess, :results)

    build_scenario!(sess)                                # inline: spawn[] is false
    dr = sess.drive
    @test dr.scenario[] isa Scenario && !dr.building[] && !drive_stale(sess)
    @test occursin("built in", dr.status[])
    res = dr.result[]
    @test res isa DriveResult
    @test [m.l for m in res.modes] == collect(MODE_DEGREES)
    # the example cone pointing is mirrored top/bottom: no odd Legendre moments; the
    # ports' azimuths are not mirror images, so odd m ≠ 0 modes stay small
    m = Dict(r.l => r for r in res.modes)
    @test abs(m[1].legendre) < 1e-9 && abs(m[3].legendre) < 1e-9
    @test m[1].epsilon < 0.01 && m[3].epsilon < 0.01
    @test 0 < m[2].epsilon < 0.2
    @test m[2].epsilon ≈ abs(m[2].legendre) rtol = 0.05        # nearly axisymmetric
    @test 1e14 < to_Wcm2(res.mean_flux) < 1e16
    @test occursin("capsule", drive_summary(res)) && occursin("ε₂", drive_summary(res))
    @test task_status(sess, :drive)[1] === :ok && startswith(task_status(sess, :drive)[2], "Tr ")
    # power accounting sums to the delivered power, less a residual below 1e-3 of it;
    # parts below that fraction share one row, as in `report`
    P = res.evaluation.deposition.delivered
    @test sum(p.power for p in res.power) ≈ P rtol = 1e-3
    @test any(p -> p.name == "Wall", res.power)
    @test all(p -> abs(p.power) >= 1e-3 * P || occursin(" < ", p.name), res.power)
    @test count(p -> occursin(" < ", p.name), res.power) <= 1
    # element fields (W/cm²)
    t = sess.target.target[]
    lv = element_values(res, :laser)
    @test length(lv) == nelements(t)
    @test all(isnan, lv[elements(t, :leh_top)]) && maximum(filter(isfinite, lv)) > 1e13
    @test all(isfinite, element_values(res, :radiosity)) && all(isfinite, element_values(res, :incident))
    @test_throws ArgumentError set_field!(dr, :temperature)
    set_field!(dr, :radiosity)
    @test dr.field[] === :radiosity
    # the capsule θ–φ map: the observed flux on a regular grid
    mp = flux_map(dr)
    @test mp.kind === :capsule && size(mp.values) == (length(mp.x), length(mp.y))
    # (from the radiosity solve; the flux response agrees to round-off)
    @test all(isapprox.(extrema(mp.values), extrema(to_Wcm2.(res.evaluation.flux)); rtol = 1e-9))
    set_view!(dr, :map)
    @test dr.view[] === :map
    @test_throws ArgumentError set_view!(dr, :movie)
    set_view!(dr, :scene)

    # pointing edits re-evaluate without a rebuild; on/off too
    s0 = dr.scenario[]
    edit_beams!(sess.beams, 1, :z, "3.4")
    @test dr.scenario[] === s0 && dr.result[] !== res
    P = dr.result[].evaluation.deposition.delivered
    set_enabled!(sess.beams, 1, false)
    @test dr.scenario[] === s0
    @test dr.result[].evaluation.deposition.delivered ≈ P - 8 * 8.0
    set_enabled!(sess.beams, 1, true)
    edit_beams!(sess.beams, 1, :z, "3.09")

    # a target change makes the scenario stale and clears the evaluation
    edit_target!(sess.target, :length, "10.0")
    @test drive_stale(sess) && dr.result[] === nothing
    @test task_status(sess, :drive)[1] === :attention
    edit_target!(sess.target, :length, "10.1")
    @test drive_stale(sess)                                  # a new Target object
    build_scenario!(sess)
    @test !drive_stale(sess) && dr.result[] isa DriveResult
end

@testset "optimization (no GL)" begin
    sess = SESS
    st = sess.optimize
    @test_throws ArgumentError set_free!(st, :spot, true)
    set_free!(st, :power, true)
    set_free!(st, :power, false)
    @test st.fields[] == [:z, :r]
    @test edit_optimize!(st, :z, "0.5") && st.deltas[][:z] == 0.5
    @test !edit_optimize!(st, :power, "1.5") && occursin("fraction", st.error[])
    @test !edit_optimize!(st, :w2, "-1")
    @test edit_optimize!(st, :w1, "0") && st.weights[][1] == 0
    @test edit_optimize!(st, :iterations, "4") && st.iterations[] == 4
    @test_throws ArgumentError set_objective!(st, :drive)
    @test !edit_optimize!(st, :nonsense, "1")

    start_optimize!(sess; spawn = false)
    r = st.result[]
    @test r isa OptimizeResult && r.completed && !st.running[]
    @test occursin(r"^(converged|reached the iteration limit)", st.status[])
    @test length(r.names) == 8 && r.names[1] == "inner23.z"
    @test r.value <= r.initial_value
    @test all(r.lower .<= r.x .<= r.upper)
    @test r.evaluations > 0 && st.evaluations[] == r.evaluations
    @test r.iterations isa Int && 0 < r.iterations <= 4         # from optimize_pointing
    @test r.reason in (:converged, :iterations)
    # the capsule keep-out is on by default and reported before and after
    @test st.capsule_keepout[] && isfinite(r.after.capsule_clearance)
    @test r.after.capsule >= 0
    @test length(r.before.modes) == 4 && length(r.after.modes) == 4
    rows = parameter_rows(r)
    @test length(rows) == 8 && rows[1].start ≈ 3.09
    @test task_status(sess, :optimize)[1] === :attention       # not applied yet

    accept_optimum!(sess)
    @test r.applied && task_status(sess, :optimize)[1] === :ok
    rws = beam_rows(sess.beams)
    @test rws[1].z ≈ r.x[1] && rws[1].r ≈ r.x[2]
    @test sess.drive.result[].modes[2].epsilon ≈ r.after.modes[2].epsilon rtol = 1e-6
    # a beam table that no longer holds the optimized beams marks the result stale
    bs0 = sess.beams.beams[]
    sess.beams.beams[] = reverse(bs0)
    @test task_status(sess, :optimize)[1] === :attention
    sess.beams.beams[] = bs0
    @test task_status(sess, :optimize)[1] === :ok

    # 3D objective with the LEH keep-out, on a subset of beams (dropped cone)
    set_objective!(st, :ylm)
    st.keepout[] = true
    set_enabled!(sess.beams, 4, false)
    st.iterations[] = 2
    start_optimize!(sess; spawn = false)
    r2 = st.result[]
    @test r2.completed && length(r2.beams) == 48 - beam_rows(sess.beams)[4].count
    @test r2.indices == findall(sess.beams.enabled[])
    set_enabled!(sess.beams, 4, true)
    st.keepout[] = false
    set_objective!(st, :legendre)

    # the foil objective needs a foil
    set_objective!(st, :irradiance)
    start_optimize!(sess; spawn = false)
    @test occursin("needs a foil", st.status[])
    set_objective!(st, :legendre)

    # cancel only sets the flag the core callback returns
    job_st = OptimizeState()
    cancel_optimize!(job_st)
    @test job_st.cancel[] && !job_st.running[]

    # cancel through optimize_pointing's callback, inline: the progress delivered
    # after the first objective evaluation cancels, and the callback stops there
    st.iterations[] = 50
    idx = findall(sess.beams.enabled[])
    seen = Int[]
    cancel_after_first = f -> (f(); push!(seen, st.evaluations[]);
                               st.evaluations[] >= 1 && cancel_optimize!(st))
    start_optimize!(st, sess.drive.scenario[], sess.beams.beams[][idx], idx;
                    deliver = cancel_after_first, spawn = false)
    rc = st.result[]
    @test !rc.completed && rc.reason === :callback && rc.evaluations == 1
    @test first(seen) == 1                                   # progress reached `deliver`
    @test rc.value ≈ rc.initial_value                        # the start is the best so far
    @test rc.x ≈ rc.x0 && occursin("canceled; kept the best design", st.status[])
    @test task_status(sess, :optimize)[1] === :attention

    # the evaluation budget ends a run whatever the line searches do
    st.iterations[] = 50
    @test edit_optimize!(st, :max_evaluations, "3") && st.max_evaluations[] == 3
    start_optimize!(sess; spawn = false)
    rb = st.result[]
    @test rb.reason === :evaluations && rb.evaluations == 3 && rb.completed
    @test occursin("evaluation limit", st.status[])
    st.max_evaluations[] = 1000
end

@testset "capsule keep-out in the optimization (no GL)" begin
    # inner cones pulled toward the waist put their envelopes on the capsule; with the
    # keep-out the optimizer moves them off
    sess = DesignSession(resolution = :coarse)
    for (cone, z, r) in (("inner23", "2.2", "-0.5"), ("inner30", "2.4", "-0.3"))
        edit_beam_node!(sess.beams, :group, cone, :z, z)
        edit_beam_node!(sess.beams, :group, cone, :r, r)
    end
    build_scenario!(sess)
    st = sess.optimize
    st.iterations[] = 40
    st.weights[] = Dict(1 => 0.0, 2 => 1.0, 3 => 0.0, 4 => 1.0)
    start_optimize!(sess; spawn = false)
    on = st.result[]
    @test on.before.capsule_clearance < 0
    @test on.after.capsule_clearance > -0.02 && on.after.capsule < 1e-4
    st.capsule_keepout[] = false
    start_optimize!(sess; spawn = false)
    off = st.result[]
    @test off.after.capsule_clearance < on.after.capsule_clearance
end

@testset "clearance (no GL)" begin
    sess = SESS
    rv = sess.review
    @test rv.targets[] == [:wall, :capsule] && rv.hardware[] == [:cap_top, :cap_bottom]
    @test_throws ArgumentError set_role!(rv, :wall, :ignore)
    set_role!(rv, :capsule, :none)
    @test rv.targets[] == [:wall]
    set_role!(rv, :capsule, :target)
    @test !edit_review!(rv, :peak_fraction, "2")
    @test edit_review!(rv, :peak_fraction, "1e-3") && rv.peak_fraction[] == 1e-3
    @test !edit_review!(rv, :pointing, "-0.1")
    @test edit_review!(rv, :focus, "0.5") && rv.errors[][:focus] == 0.5
    kw = clearance_kwargs(rv)
    @test kw.targets == (:wall, :capsule) && kw.level == 1 && kw.focus == 0.5

    # three quads
    bt = sess.beams
    keep = [1, 9, 17]
    for i in eachindex(bt.enabled[])
        bt.enabled[][i] = i in keep
    end
    notify(bt.enabled)
    start_clearance!(sess; spawn = false)
    rep = rv.report[]
    @test rep isa ClearanceReport && length(rep.rows) == 3 * 2
    @test rv.progress[] == (3, 3) && startswith(rv.status[], "done")
    @test Set(x.beam for x in rep.rows) == Set(b.name for b in bt.beams[][keep])
    rows = clearance_rows(rv)
    rank(s) = s === :fail ? 0 : s === :review ? 1 : 2
    @test issorted([(rank(x.status), x.clearance) for x in rows])
    @test [x.component for x in worst_rows(rv)] == [:cap_top, :cap_bottom]
    @test occursin("3 beams", review_brief(rv))
    @test !rep.canceled
    # the batched check equals one clearance call without a callback
    ref = clearance(sess.target.target[], active_beams(bt); clearance_kwargs(rv)...)
    key(r) = sort([(x.beam, x.component, x.clearance) for x in r.rows])
    a, b = key(rep), key(ref)
    @test [x[1:2] for x in a] == [x[1:2] for x in b] && all(((x, y),) -> x[3] ≈ y[3], zip(a, b))
    path = joinpath(mktempdir(), "clearance.csv")
    export_clearance(rv, path)
    @test countlines(path) == 1 + length(rep.rows)
    @test !review_stale(sess)
    edit_beams!(bt, 1, :z, "3.2")
    @test review_stale(sess) && task_status(sess, :clearance)[1] === :attention

    # cancel through clearance's callback, inline: after the first batch of 2 beams
    progress = Tuple{Int,Int}[]
    cancel_after_batch = f -> (f(); push!(progress, rv.progress[]);
                               rv.progress[][1] >= 2 && cancel_clearance!(rv))
    four = active_beams(bt)[1:3]
    start_clearance!(rv, sess.target.target[], [four; bt.beams[][25]]; spawn = false, batch = 2,
                     deliver = cancel_after_batch)
    rc = rv.report[]
    @test progress[1:2] == [(0, 4), (2, 4)]                   # done = 0 first, then per batch
    @test rc.canceled && length(unique(x.beam for x in rc.rows)) == 2
    @test rv.status[] == "canceled after 2 of 4 beams" && !rv.running[]
    @test occursin("canceled, partial", review_brief(rv))
    @test task_status(sess, :clearance)[1] === :attention
    # canceled before any beam: no report
    start_clearance!(rv, sess.target.target[], four; spawn = false, batch = 2,
                     deliver = f -> (f(); cancel_clearance!(rv)))
    @test rv.report[] === nothing && rv.status[] == "canceled after 0 of 3 beams"

    rv2 = ReviewState()
    start_clearance!(rv2, sess.target.target[], Beam{Float64}[]; spawn = false)
    @test rv2.status[] == "no beams are enabled" && rv2.report[] === nothing
    @test_throws ArgumentError export_clearance(rv2, path)
    fill!(bt.enabled[], true)
    notify(bt.enabled)
end

@testset "background jobs through a delivery queue (no GL)" begin
    # As in the window: jobs run on workers and their updates wait in a queue for
    # the GUI thread; nothing changes until the queue is drained.
    q = Channel{Any}(Inf)
    sess = DesignSession(resolution = :coarse, deliver = f -> put!(q, f))
    sess.spawn[] = true
    build_scenario!(sess)
    @test sess.drive.building[]
    @test task_status(sess, :drive)[1] === :busy
    drain_until(() -> sess.drive.result[] !== nothing, q)
    @test !sess.drive.building[] && sess.drive.scenario[] isa Scenario

    # a stale build result is dropped: the later request wins
    build_scenario!(sess)
    g = sess.drive.generation[]
    build_scenario!(sess)
    @test sess.drive.generation[] == g + 1
    drain_until(() -> !sess.drive.building[] && !sess.drive.evaluating[], q)
    @test sess.drive.result[] isa DriveResult

    # optimize on a worker, cancel once it reports progress: the best design so far
    st = sess.optimize
    st.iterations[] = 200
    set_free!(st, :defocus, true)
    start_optimize!(sess; spawn = true)
    @test st.running[]
    drain_until(() -> st.evaluations[] > 0 || !st.running[], q)
    cancel_optimize!(sess)
    drain_until(() -> !st.running[], q)
    r = st.result[]
    @test r isa OptimizeResult
    @test !r.completed && occursin("canceled", st.status[]) || r.completed
    @test r.value <= r.initial_value
    @test length(r.x) == 12

    # clearance on a worker with progress, canceled after the first batch
    bt = sess.beams
    for i in eachindex(bt.enabled[])
        bt.enabled[][i] = i <= 4
    end
    notify(bt.enabled)
    rv = sess.review
    start_clearance!(sess; spawn = true, batch = 2)
    @test rv.running[] && task_status(sess, :clearance)[1] === :busy
    drain_until(() -> rv.progress[][1] > 0 || !rv.running[], q)
    cancel_clearance!(sess)
    drain_until(() -> !rv.running[], q)
    @test rv.report[] isa ClearanceReport
    n = length(unique(x.beam for x in rv.report[].rows))
    @test n == 2 && rv.report[].canceled && occursin("canceled", rv.status[]) ||
          n == 4 && !rv.report[].canceled
    @test review_batch() >= 4
    Controllers._abandon_jobs!(sess)
    @test !sess.drive.building[]
end

@testset "foil target (no GL)" begin
    sess = DesignSession(kind = :foil, resolution = :coarse)
    @test sess.beams.preset[] === :omega_backlighter
    @test sess.review.targets[] == [FOIL_NAME] && sess.review.hardware[] == [:main]
    build_scenario!(sess)
    res = sess.drive.result[]
    @test isempty(res.modes) && isnan(res.mean_flux)
    @test res.spot !== nothing && res.spot.mean > 0 && 0 < res.spot.fraction <= 1
    @test occursin("spot", drive_summary(res))
    mp = flux_map(sess.drive)
    @test mp.kind === :foil && maximum(mp.values) > 0
    sess.drive.spot_radius[] = 0.3
    @test sess.drive.result[].spot.radius == 0.3              # re-evaluated

    st = sess.optimize
    set_objective!(st, :irradiance)
    st.fields[] = [:r, :phi]
    st.by[] = :beam
    st.iterations[] = 2
    start_optimize!(sess; spawn = false)
    r = st.result[]
    @test r.completed && r.value <= r.initial_value
    @test r.before.spot !== nothing && r.after.spot !== nothing
    set_objective!(st, :legendre)
    start_optimize!(sess; spawn = false)
    @test occursin("needs a hohlraum", st.status[])

    # back to a hohlraum: the NIF quads return
    set_target_kind!(sess.target, :hohlraum)
    @test sess.beams.preset[] === :nif_quads && drive_stale(sess)
    @test sess.review.hardware[] == [:cap_top, :cap_bottom]
end

@testset "windows in the drive (no GL)" begin
    sess = DesignSession(resolution = :coarse)
    # a closed and an open window: separate components, a valid closed enclosure
    @test add_window!(sess.target)
    @test add_window!(sess.target; phi = 200.0, z = 2.0, width = 1.5, height = 1.5, open = true)
    build_scenario!(sess)
    res = sess.drive.result[]
    t = sess.target.target[]
    @test res isa DriveResult && res.evaluation.scenario.target === t
    @test sum(p.power for p in res.power) ≈ res.evaluation.deposition.delivered rtol = 1e-3
    @test any(p -> p.name == "window_1", res.power)
    @test all(isfinite, element_values(res, :radiosity))
    @test all(isnan, element_values(res, :laser)[elements(t, :window_2)])    # a hole
    # x-rays leave through the open window: its incident flux is a loss
    @test sum(res.solution.H[e] * t.elem_area[e] for e in elements(t, :window_2)) > 0
end
