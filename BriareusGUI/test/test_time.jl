# The time panel's data (no GL): the pulses plotted, the drive history and the fields it
# keeps per step, the time cursor and its step, staleness, the background job, and the
# retune's power histories.

@testset "time panel: pulses and charts (no GL)" begin
    sess = DesignSession(resolution = :coarse)
    tm = sess.time
    for c in TIME_CHARTS
        ch = time_chart(sess, c)
        @test !isempty(ch.empty) && isempty(ch.left.series) && isnan(ch.cursor)
    end
    @test_throws ArgumentError time_chart(sess, :nope)
    start_history!(sess)                                       # nothing to fire yet
    @test occursin("load a pulse table", tm.status[]) && tm.history[] === nothing

    beams = sess.beams.beams[]
    pul = cone_fraction_pulses(beams, [0.0, 1.0, 2.0, 3.0, 4.0], [0.0, 100.0, 200.0, 300.0, 0.0],
                               [0.3, 0.3, 0.28, 0.25, 0.25])
    dir = mktempdir()
    table = write_pulse_table(joinpath(dir, "pulses.csv"), pul)
    load_pulses!(sess, table)
    @test tm.pulses[].by === :group && tm.chart[] === :pulses
    # the pulses chart: each group's power is the sum over its enabled beams
    ch = time_chart(sess)
    @test ch.left.label == "Laser power (TW)" && ch.xlabel == "t (ns)" && isempty(ch.empty)
    groups = unique(sess.beams.cluster[])
    @test [s.label for s in ch.left.series] == [String.(groups); "total"]
    @test ch.left.series[end].style === :dash
    @test [r.label for r in ch.legend] == [s.label for s in ch.left.series]
    cc = cluster_colors(sess.beams)
    for (g, s) in zip(groups, ch.left.series)
        n = count(b -> b.group === g, beams)
        @test s.x == [0.0, 1.0, 2.0, 3.0, 4.0]
        @test s.y ≈ [n * pul[g](τ) for τ in s.x]
        @test s.color == cc[g]                                  # the model tree's swatch
    end
    @test ch.left.series[end].y ≈ sum(s.y for s in ch.left.series[1:end-1])
    # a disabled beam leaves its group's power
    g1 = groups[1]
    i = findfirst(b -> b.group === g1, beams)
    en = copy(sess.beams.enabled[]); en[i] = false; sess.beams.enabled[] = en
    s1 = time_chart(sess).left.series[1]
    @test s1.y ≈ [(count(b -> b.group === g1, beams) - 1) * pul[g1](τ) for τ in s1.x]
    en[i] = true; sess.beams.enabled[] = en
    # a table whose columns are not the beams' groups, quads or names
    bad = write_pulse_table(joinpath(dir, "bad.csv"), Dict(:x => PulseShape([0.0, 1.0], [1.0, 1.0])))
    @test_throws ArgumentError load_pulses!(sess, bad)
end

@testset "time panel: drive history and the time cursor (no GL)" begin
    sess = DesignSession(resolution = :coarse)
    tm = sess.time
    use_pulse_shape!(sess; step = 1.0)
    @test tm.pulses[].kind === :shape && tm.pulses[].times == collect(0.0:1.0:9.0)
    @test time_step(tm) === nothing
    # no drive yet: computed first (inline here), then the history
    start_history!(sess)
    @test sess.drive.result[] !== nothing && !tm.queued[] && !tm.running[]
    h = tm.history[]
    @test h isa DriveHistory && tm.chart[] === :drive && startswith(tm.status[], "$(length(h.time)) steps")
    @test tm.progress[] == (length(h.time), length(h.time))
    shape = tm.pulses[].shape
    @test h.time == filter(τ -> shape(τ) > 0, 0.0:1.0:9.0)
    @test all(M -> size(M) == (nelements(sess.drive.scenario[].target), length(h.time)), values(h.fields))
    @test Set(keys(h.fields)) == Set(DRIVE_FIELDS)

    # independent check: gold's albedo is constant, so each step is the drive at the
    # table's powers scaled by the pulse shape there
    res = sess.drive.result[]
    for (k, τ) in enumerate(h.time)
        f = shape(τ)
        @test h.fields[:incident][:, k] ≈ f .* to_Wcm2.(res.solution.H) rtol = 1e-4
        @test h.Tr[k] ≈ radiation_temperature(f * res.mean_flux) rtol = 1e-4
        @test h.legendre[2, k] ≈ res.modes[2].legendre atol = 1e-6
    end
    lf = element_values(res, :laser)
    @test all(isnan.(h.fields[:laser][:, 1]) .== isnan.(lf))
    @test isfinite(h.dante[1]) && 0.5 < h.dante[argmax(h.Tr)] / maximum(h.Tr) < 1.5
    # the energy split: cumulative, and the laser's energy accounted for
    E = Dict(s.name => s.values for s in h.energy)
    @test Set(keys(E)) == Set([:laser, :wall, :capsule, :holes])
    @test E[:laser][end] ≈ sum((h.time[k+1] - h.time[k]) * (h.history.setup[k] + h.history.setup[k+1]) / 2
                               for k in 1:length(h.time)-1)
    @test all(v -> issorted(v), values(E))
    @test 0.9 < (E[:wall][end] + E[:capsule][end] + E[:holes][end]) / E[:laser][end] < 1.001   # (view-factor closure)
    for c in (:drive, :energy)
        ch = time_chart(sess, c)
        @test isempty(ch.empty) && all(s -> length(s.x) == length(s.y) == length(h.time), ch.left.series)
    end
    ch = time_chart(sess, :drive)
    @test [s.label for s in ch.left.series] == ["Capsule Tr", "Dante 1"]
    @test [s.label for s in ch.right.series] == ["a₂/a₀", "a₄/a₀"]
    @test ch.right.series[1].y ≈ 100 .* h.legendre[2, :]
    @test ch.left.label == "Tr (eV)" && ch.right.label == "aₗ/a₀ (%)"
    @test time_chart(sess, :energy).left.label == "Energy (kJ)"

    # the cursor: the nearest solved time, clamped to the history; colors the 3D view
    @test sess.view.color[] === :component
    set_time_cursor!(sess, 3.4)
    @test time_step(tm) == findfirst(==(3.0), h.time) && time_label(sess) == "t = 3.00 ns"
    @test sess.view.color[] === :incident                     # from the components
    @test cursor_values(sess, :incident) == Float64.(h.fields[:incident][:, time_step(tm)])
    @test cursor_range(sess, :incident) == h.range[:incident]
    @test h.range[:incident][2] ≈ maximum(h.fields[:incident])
    @test time_chart(sess).cursor == 3.0
    set_time_cursor!(sess, -5.0)
    @test time_step(tm) == 1
    set_time_cursor!(sess, 50.0)
    @test time_step(tm) == length(h.time)
    set_time_step!(sess, 2)
    @test tm.cursor[] == h.time[2]
    set_color!(sess, :radiosity)
    @test cursor_values(sess, :radiosity) == Float64.(h.fields[:radiosity][:, 2])
    set_time_cursor!(sess, NaN)                                # back to the table's powers
    @test time_step(tm) === nothing && cursor_values(sess, :radiosity) === nothing && time_label(sess) == ""
    @test sess.view.color[] === :radiosity

    # a beam edit makes the history stale: the 3D view leaves it
    set_time_step!(sess, 3)
    @test !history_stale(sess)
    edit_beams!(sess.beams, 1, :power, "1.0")
    @test history_stale(sess) && cursor_values(sess, :radiosity) === nothing && time_label(sess) == ""
    @test occursin("out of date", time_chart(sess, :drive).title)
    start_history!(sess)
    @test !history_stale(sess) && tm.history[] !== h
    @test tm.history[].time[3] == 2.0 && time_label(sess) == "t = 2.00 ns"
    # a target without a capsule: no history, and why
    set_target_kind!(sess.target, :foil)
    use_pulse_shape!(sess; step = 1.0)
    start_history!(sess)
    @test occursin("follows the capsule", tm.status[])
end

@testset "time panel: the history on a worker (no GL)" begin
    q = Channel{Any}(Inf)
    sess = DesignSession(resolution = :coarse, deliver = f -> put!(q, f))
    sess.spawn[] = true
    tm = sess.time
    use_pulse_shape!(sess; step = 0.5)
    start_history!(sess)                          # the drive first, on a worker
    @test tm.queued[] && sess.drive.building[]
    drain_until(() -> tm.running[], q)            # then the history, queued behind it
    n = tm.progress[][2]
    @test n == count(τ -> tm.pulses[].shape(τ) > 0, tm.pulses[].times)
    seen = Int[]
    on(p -> push!(seen, p[1]), tm.progress)
    drain_until(() -> !tm.running[], q)
    @test tm.history[] isa DriveHistory && length(tm.history[].time) == n
    @test issorted(seen) && last(seen) == n
    # canceled: nothing kept
    h = tm.history[]
    use_pulse_shape!(sess; step = 0.1)
    start_history!(sess)
    @test tm.running[]
    cancel_history!(sess)
    drain_until(() -> !tm.running[], q)
    @test tm.status[] == "drive history canceled" && tm.history[] === h
end

@testset "time panel: the retune's power histories (no GL)" begin
    s = DesignSession(resolution = :coarse)
    dir = mktempdir()
    beams = s.beams.beams[]
    pul = cone_fraction_pulses(beams, [0.0, 1.0, 3.0, 4.0], [0.0, 100.0, 300.0, 0.0],
                               [0.3, 0.3, 0.25, 0.25])
    table = write_pulse_table(joinpath(dir, "reference.csv"), pul)
    out = joinpath(dir, "retuned.csv")
    build_scenario!(s)
    s.beams.enabled[] = [nif_bundle(b) != "31" for b in beams]          # bundle 31 down
    retune_to_reference!(s, table, out; cap = 1.5)
    @test isfile(out)                                                   # the CSV, as before
    tm = s.time
    @test tm.chart[] === :retune && tm.retune[] !== nothing
    @test tm.pulses[].by === :group && tm.pulses[].label == "reference.csv"   # the reference fires
    ch = time_chart(s)
    groups = unique(s.beams.cluster[])
    ref, new = ch.left.series[1:length(groups)], ch.left.series[length(groups)+1:end]
    @test all(r -> r.style === :dash, ref) && all(r -> r.style === :solid, new)
    @test [r.label for r in ref] == [r.label for r in new] == String.(groups)
    back = read_pulse_table(out)
    for (g, r, n) in zip(groups, ref, new)
        @test r.y ≈ [count(b -> b.group === g, beams) * pul[g](τ) for τ in r.x]
        # the retuned powers: the enabled beams on the written pulses
        @test n.y ≈ [sum(back[Symbol(b.quad)](τ) for (b, on) in zip(beams, s.beams.enabled[]) if on && b.group === g)
                     for τ in n.x]
        @test all(n.y .<= 1.5 .* r.y .+ 1e-9)
    end
    @test ch.right.label == "ε₁ (%)" && [x.label for x in ch.right.series] == ["ε₁ before", "ε₁ retuned"]
    @test maximum(ch.right.series[2].y) < maximum(ch.right.series[1].y) / 10
    @test any(r -> r.label == "reference (dashed)", ch.legend)
end
