# The Drive, Optimize and Review steps through a DesignSession, inline (as a script
# would run them) and through a window-style delivery queue. A coarse hohlraum keeps
# the view factors to seconds.

const SESS = DesignSession(resolution = :coarse)

@testset "Drive step (no GL)" begin
    sess = SESS
    @test step_status(sess, :target)[1] === :ok
    @test step_status(sess, :beams)[1] === :ok
    @test step_status(sess, :drive) == (:todo, "not computed")
    @test step_status(sess, :optimize)[1] === :todo
    @test step_status(sess, :review)[1] === :todo
    set_step!(sess, :drive)
    @test sess.step[] === :drive
    @test_throws ArgumentError set_step!(sess, :results)

    build_scenario!(sess)                                # inline: spawn[] is false
    dr = sess.drive
    @test dr.scenario[] isa Scenario && !dr.building[] && !drive_stale(sess)
    @test occursin("built in", dr.status[])
    res = dr.result[]
    @test res isa DriveResult
    @test [m.l for m in res.modes] == collect(MODE_DEGREES)
    # the published pointing is mirrored top/bottom: no odd Legendre moments; the
    # ports' azimuths are not mirror images, so odd m ≠ 0 modes stay small
    m = Dict(r.l => r for r in res.modes)
    @test abs(m[1].legendre) < 1e-9 && abs(m[3].legendre) < 1e-9
    @test m[1].epsilon < 0.01 && m[3].epsilon < 0.01
    @test 0 < m[2].epsilon < 0.2
    @test m[2].epsilon ≈ abs(m[2].legendre) rtol = 0.05        # nearly axisymmetric
    @test 1e14 < to_Wcm2(res.mean_flux) < 1e16
    @test occursin("capsule", drive_summary(res)) && occursin("ε₂", drive_summary(res))
    @test step_status(sess, :drive)[1] === :ok
    # power accounting sums to the delivered power
    @test sum(p.power for p in res.power) ≈ res.evaluation.deposition.delivered rtol = 1e-6
    @test any(p -> p.name == "wall", res.power)
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
    @test extrema(mp.values) == extrema(to_Wcm2.(res.evaluation.flux))
    set_view!(dr, :map)
    @test dr.view[] === :map
    @test_throws ArgumentError set_view!(dr, :image)
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
    @test step_status(sess, :drive)[1] === :attention
    edit_target!(sess.target, :length, "10.1")
    @test drive_stale(sess)                                  # a new Target object
    build_scenario!(sess)
    @test !drive_stale(sess) && dr.result[] isa DriveResult
end

@testset "Optimize step (no GL)" begin
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
    @test startswith(st.status[], "done")
    @test length(r.names) == 8 && r.names[1] == "inner23.z"
    @test r.value <= r.initial_value
    @test all(r.lower .<= r.x .<= r.upper)
    @test r.evaluations > 0 && st.evaluations[] == r.evaluations
    @test length(r.before.modes) == 4 && length(r.after.modes) == 4
    rows = parameter_rows(r)
    @test length(rows) == 8 && rows[1].start ≈ 3.09
    @test step_status(sess, :optimize)[1] === :attention       # not applied yet

    accept_optimum!(sess)
    @test r.applied && step_status(sess, :optimize)[1] === :ok
    rws = beam_rows(sess.beams)
    @test rws[1].z ≈ r.x[1] && rws[1].r ≈ r.x[2]
    @test sess.drive.result[].modes[2].epsilon ≈ r.after.modes[2].epsilon rtol = 1e-6

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

    # cancel only sets the flag the running objective checks
    job_st = OptimizeState()
    cancel_optimize!(job_st)
    @test job_st.cancel[] && !job_st.running[]
end

@testset "Review step (no GL)" begin
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

    # three quads, one per chunk boundary
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
    # the chunked check equals one clearance call
    ref = clearance(sess.target.target[], active_beams(bt); clearance_kwargs(rv)...)
    key(r) = sort([(x.beam, x.component, x.clearance) for x in r.rows])
    a, b = key(rep), key(ref)
    @test [x[1:2] for x in a] == [x[1:2] for x in b] && all(((x, y),) -> x[3] ≈ y[3], zip(a, b))
    path = joinpath(mktempdir(), "clearance.csv")
    export_clearance(rv, path)
    @test countlines(path) == 1 + length(rep.rows)
    @test !review_stale(sess)
    edit_beams!(bt, 1, :z, "3.2")
    @test review_stale(sess) && step_status(sess, :review)[1] === :attention

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
    @test step_status(sess, :drive)[1] === :busy
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

    # clearance on a worker with progress, canceled after the first chunk
    bt = sess.beams
    for i in eachindex(bt.enabled[])
        bt.enabled[][i] = i <= 2 * Controllers.REVIEW_CHUNK
    end
    notify(bt.enabled)
    rv = sess.review
    start_clearance!(sess; spawn = true)
    @test rv.running[] && step_status(sess, :review)[1] === :busy
    drain_until(() -> rv.progress[][1] > 0 || !rv.running[], q)
    cancel_clearance!(sess)
    drain_until(() -> !rv.running[], q)
    @test rv.report[] isa ClearanceReport
    n = length(unique(x.beam for x in rv.report[].rows))
    @test n == Controllers.REVIEW_CHUNK && occursin("canceled", rv.status[]) ||
          n == 2 * Controllers.REVIEW_CHUNK
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
