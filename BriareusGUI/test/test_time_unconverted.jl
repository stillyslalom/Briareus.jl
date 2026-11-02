# The time panel with unconverted (1ω, 2ω) light: the drive history and the retune take
# the drive's ratios on the beams whose light takes part, as the core's `drive_history`
# and `retune_drive` take `unconverted`; without it they are the 3ω ones, bit for bit.
# The unconverted light lands ≥ 4.8 mm off the 3ω axis (`nif_unconverted`), outside a
# 2.9 mm hohlraum: a mock hohlraum 6 mm in radius with 5.5 mm entrance holes takes some in.

# A session on the wide mock hohlraum.
function _wide_hohlraum_session()
    s = DesignSession(resolution = :coarse)
    edit_part!(s.target, 1, :radius, "6.0") && edit_part!(s.target, 1, :leh_radius, "5.5") ||
        error("the mock hohlraum did not build")
    return s
end

# Equal to 1e-12 relative: numbers, arrays, or Dicts of them with the same keys.
_near(a, b) = a isa AbstractDict ? keys(a) == keys(b) && all(isapprox(a[k], b[k]; rtol = 1e-12) for k in keys(a)) :
              isapprox(a, b; rtol = 1e-12)

# The core's drive history of a session's latest pulses on its scenario, with the GUI
# history's observed parts, wall and times, and the keywords `kw`.
function _core_history(sess, times; kw...)
    s = sess.drive.scenario[]
    src, beams, enabled = sess.time.pulses[], sess.beams.beams[], sess.beams.enabled[]
    observe = Controllers._observed_parts(s)
    names = component_names(s.target)
    walls = Tuple(n for n in wall_names(sess.target) if n in names)
    return drive_history(s.model, s.deposition, τ -> Controllers._pulse_beams(src, beams, enabled, τ), times;
                         observe, wall = isempty(walls) ? observe : walls, lmax = 4, kw...)
end

_history_fields(h) = (h.Tr, h.flux, h.moments, h.wall_Tr, h.setup, h.delivered, h.laser, h.source, h.absorbed)

@testset "time panel: drive histories with unconverted light (no GL)" begin
    sess = _wide_hohlraum_session()
    tm, bt, dr = sess.time, sess.beams, sess.drive
    use_pulse_shape!(sess; step = 1.0)
    start_history!(sess)
    h0 = tm.history[]
    @test h0 isa DriveHistory && h0.history.unconverted === nothing
    s = dr.scenario[]
    src, beams, enabled = tm.pulses[], bt.beams[], bt.enabled[]
    at(τ) = Controllers._pulse_beams(src, beams, enabled, τ)
    observe = Controllers._observed_parts(s)
    names = component_names(s.target)
    walls = Tuple(n for n in wall_names(sess.target) if n in names)
    wall = isempty(walls) ? observe : walls
    core(; kw...) = drive_history(s.model, s.deposition, at, h0.time; observe, wall, lmax = 4, kw...)
    fields(h) = (h.Tr, h.flux, h.moments, h.wall_Tr, h.setup, h.delivered, h.laser, h.source, h.absorbed)

    # without unconverted light, the 3ω history bit for bit (the core's without the keyword)
    ref = core()
    @test isequal(fields(h0.history), fields(ref))
    # a ratio with no beam's 1ω or 2ω light taking part: still the 3ω history
    set_unconverted_ratio!(dr, 0.05, 0.02)
    @test !history_stale(sess)
    start_history!(sess)
    @test isequal(fields(tm.history[].history), fields(ref)) && tm.history[].history.unconverted === nothing

    # every beam's 1ω and 2ω light: stale, then the core's history with the same keywords
    all_beams = collect(eachindex(beams))
    set_participation!(bt, all_beams, 1, true)
    @test history_stale(sess)
    set_participation!(bt, all_beams, 2, true)
    start_history!(sess)
    h1 = tm.history[]
    @test !history_stale(sess) && h1.history.unconverted !== nothing
    uc = core(; unconverted = (omega1 = 0.05, omega2 = 0.02), colors = nif_unconverted())
    near(a, b) = a isa AbstractDict ? keys(a) == keys(b) && all(isapprox(a[k], b[k]; rtol = 1e-12) for k in keys(a)) :
                  isapprox(a, b; rtol = 1e-12)
    @test all(near.(fields(h1.history), fields(uc)))
    @test h1.history.unconverted.delivered == uc.unconverted.delivered
    u = dr.result[].unconverted                                 # some of it absorbed (not vacuous)
    @test u.omega1.absorbed > 0 && u.omega2.absorbed > 0
    @test all(h1.Tr .> h0.Tr)                                   # the unconverted light adds drive
    # the energy chart: each color's delivered energy, part of the laser's
    E = Dict(x.name => x.values for x in h1.energy)
    upto(p) = [sum((h1.time[k+1] - h1.time[k]) * (p[k] + p[k+1]) / 2 for k in 1:j-1; init = 0.0)
               for j in eachindex(h1.time)]
    @test E[:omega1] ≈ upto(uc.unconverted.delivered[:omega1]) rtol = 1e-12
    @test E[:omega2] ≈ upto(uc.unconverted.delivered[:omega2]) rtol = 1e-12
    @test E[:laser][end] ≈ Dict(x.name => x.values for x in h0.energy)[:laser][end] + E[:omega1][end] + E[:omega2][end] rtol = 1e-12
    ch = time_chart(sess, :energy)
    @test count(x -> x.label in ("1ω delivered", "2ω delivered") && x.style === :dot, ch.left.series) == 2

    # a changed ratio makes it stale; back to it, current again
    set_unconverted_ratio!(dr, 0.06, 0.02)
    @test history_stale(sess) && occursin("out of date", time_chart(sess, :drive).title)
    set_unconverted_ratio!(dr, 0.05, 0.02)
    @test !history_stale(sess)
    # a color on a disabled beam only does not count
    bt.enabled[] = [i != 1 for i in all_beams]
    set_participation!(bt, [1], 1, false)
    @test Controllers._history_unconverted(sess) == (omega1 = 0.05, omega2 = 0.02)
    # 1ω on some of the enabled beams only (all but the first): a fraction per beam, the
    # core's history with the same fractions keyed by beam name
    bt.enabled[] = fill(true, length(beams))
    @test history_stale(sess)
    start_history!(sess)
    h2 = tm.history[]
    @test h2 !== h1 && !history_stale(sess) && h2.time == h0.time
    @test allunique(Symbol(b.name) for b in beams)
    keyed = Dict(Symbol(b.name) => (i == 1 ? 0.0 : 0.05) for (i, b) in enumerate(beams))
    uk = core(; unconverted = (omega1 = keyed, omega2 = 0.02), by = :beam, colors = nif_unconverted())
    @test all(near.(fields(h2.history), fields(uk)))
    @test all(near(h2.history.unconverted.delivered[c], uk.unconverted.delivered[c]) for c in (:omega1, :omega2))
    # (the first beam's 1ω light is left out: less 1ω delivered than on every beam, the same 2ω)
    @test all(0 .< uk.unconverted.delivered[:omega1] .< uc.unconverted.delivered[:omega1]) &&
          near(uk.unconverted.delivered[:omega2], uc.unconverted.delivered[:omega2])
end

@testset "time panel: a retune with unconverted light (no GL)" begin
    s = _wide_hohlraum_session()
    dir = mktempdir()
    beams = s.beams.beams[]
    pul = cone_fraction_pulses(beams, [0.0, 1.0, 3.0, 4.0], [0.0, 100.0, 300.0, 0.0],
                               [0.3, 0.3, 0.25, 0.25])
    table = write_pulse_table(joinpath(dir, "reference.csv"), pul)
    out = joinpath(dir, "retuned.csv")
    build_scenario!(s)
    all_beams = collect(eachindex(beams))
    s.beams.enabled[] = [nif_bundle(b) != "31" for b in beams]          # bundle 31 down
    keep = findall(s.beams.enabled[])
    kept = beams[keep]
    # 1ω on all but one kept beam: the retuned drive is a core history of the kept beams
    # on the new pulses with the same fractions, keyed by beam name
    set_unconverted_ratio!(s.drive, 0.05, 0.02)
    off = first(keep)
    set_participation!(s.beams, [i for i in all_beams if i != off], 1, true)
    retune_to_reference!(s, table, out; cap = 1.5)
    @test isfile(out)
    rm(out)
    r = s.time.retune[]
    sc = s.drive.scenario[]
    observe = Controllers._observed_parts(sc)
    names = component_names(sc.target)
    walls = Tuple(n for n in wall_names(s.target) if n in names)
    wall = isempty(walls) ? observe : walls
    @test r.reference.unconverted !== nothing && r.after.unconverted !== nothing
    @test allunique(Symbol(b.name) for b in beams)
    keyed = Dict(Symbol(b.name) => (i == off ? 0.0 : 0.05) for (i, b) in enumerate(beams))
    after = drive_history(sc.model, sc.deposition[keep], τ -> with_pulses(kept, r.pulses, τ; by = :quad), r.after.time;
                          observe, wall, lmax = 4, ylm = true, unconverted = (omega1 = keyed, omega2 = 0.0), by = :beam,
                          colors = nif_unconverted())
    @test isapprox(after.Tr, r.after.Tr; rtol = 1e-12) && isapprox(after.moments, r.after.moments; rtol = 1e-12)
    @test isapprox(after.unconverted.delivered[:omega1], r.after.unconverted.delivered[:omega1]; rtol = 1e-12) &&
          all(>(0), after.unconverted.delivered[:omega1])

    # every beam's 1ω and 2ω light
    set_participation!(s.beams, all_beams, 1, true)
    set_participation!(s.beams, all_beams, 2, true)
    retune_to_reference!(s, table, out; cap = 1.5)
    @test isfile(out)
    r = s.time.retune[]
    sc = s.drive.scenario[]
    @test r.reference.unconverted !== nothing && r.after.unconverted !== nothing

    # the retuned drive: a core history of the kept beams on the new pulses, with the
    # same unconverted light
    back = r.pulses
    uc = (; unconverted = (omega1 = 0.05, omega2 = 0.02), colors = nif_unconverted())
    after = drive_history(sc.model, sc.deposition[keep], τ -> with_pulses(kept, back, τ; by = :quad), r.after.time;
                          observe, wall, lmax = 4, ylm = true, uc...)
    @test isapprox(after.Tr, r.after.Tr; rtol = 1e-12) && isapprox(after.moments, r.after.moments; rtol = 1e-12)
    @test after.unconverted.delivered == r.after.unconverted.delivered

    # the retune's prediction (linear about the reference's albedos) against that history:
    # the capsule's level and shape
    v = retune_pulses(sc, beams, s.beams.enabled[], pul; by = :quad, cap = 1.5,
                      wall = Tuple(wall_names(s.target)), unconverted = uc.unconverted)
    @test v.after.Tr ≈ r.after.Tr rtol = 1e-12
    pred = v.retune
    c00 = v.after.ylm[1, :]
    @test isapprox(pred.drive, c00 ./ v.reference.ylm[1, :]; rtol = 1e-3)
    @test isapprox(pred.shape, v.after.ylm ./ c00'; atol = 1e-3)
    @test maximum(mode_amplitude(v.after.ylm[:, k], 1) for k in eachindex(v.after.time)) <
          maximum(mode_amplitude(v.before.ylm[:, k], 1) for k in eachindex(v.before.time)) / 10

    # with a 1ω prepulse over the reference's first firing time (1 ns): the reference and
    # the retuned drive take it, as the core's histories with one prepulse for every beam
    r0 = r
    P = 2e-3
    set_prepulse!(s.time, P, 0.0; start = 0.0, stop = 2.0)
    retune_to_reference!(s, table, out; cap = 1.5)
    r = s.time.retune[]
    pre = (omega1 = PulseShape([0.0, 2.0], [P, P]),)
    after = drive_history(sc.model, sc.deposition[keep], τ -> with_pulses(kept, r.pulses, τ; by = :quad), r.after.time;
                          observe, wall, lmax = 4, ylm = true, uc..., prepulse = pre)
    @test isapprox(after.Tr, r.after.Tr; rtol = 1e-12) && isapprox(after.moments, r.after.moments; rtol = 1e-12)
    @test all(_near(after.unconverted.delivered[c], r.after.unconverted.delivered[c]) for c in (:omega1, :omega2))
    ref = drive_history(sc.model, sc.deposition, τ -> with_pulses(beams, r.reference_pulses, τ; by = r.reference_by),
                        r.reference.time;
                        observe, wall, lmax = 4, ylm = true, uc..., prepulse = pre)
    @test _near(ref.Tr, r.reference.Tr) && _near(ref.unconverted.delivered[:omega1], r.reference.unconverted.delivered[:omega1])
    # (more 1ω light inside the window than without the prepulse, the same after it)
    d0, d1 = r0.reference.unconverted.delivered[:omega1], r.reference.unconverted.delivered[:omega1]
    inside = r.reference.time .<= 2.0
    @test r.reference.time == r0.reference.time && any(inside) && all(d1[inside] .> d0[inside]) &&
          _near(d1[.!inside], d0[.!inside])
    # a prepulse on some beams of a retuned key only cannot be keyed by it
    @test any(b -> b.group === beams[1].group, beams[2:end])
    pp =(omega1 = [i == 1 ? P : 0.0 for i in all_beams], omega2 = zeros(length(beams)), window = (0.0, 2.0))
    @test_throws ArgumentError Controllers._prepulse_schedules(pp, beams, :group)
    @test Controllers._prepulse_schedules(pp, beams, :beam)[:omega1][Symbol(beams[1].name)](1.0) == P
end

@testset "time panel: a prepulse and the irradiance checks (no GL)" begin
    sess = _wide_hohlraum_session()
    tm, bt, dr = sess.time, sess.beams, sess.drive
    use_pulse_shape!(sess; step = 1.0)
    beams = bt.beams[]
    all_beams = collect(eachindex(beams))
    set_unconverted_ratio!(dr, 0.05, 0.02)
    set_participation!(bt, all_beams, 1, true)
    set_participation!(bt, all_beams, 2, true)
    start_history!(sess)
    h1 = tm.history[]
    @test tm.prepulse[] == DEFAULT_PREPULSE && tm.threshold[] == DEFAULT_THRESHOLD == 1e8

    # a 1ω prepulse makes the history stale; zero power, whatever the window, does not
    set_prepulse!(tm, 0.0, 0.0; start = -4.0, stop = 2.5)
    @test !history_stale(sess)
    P = 2e-3
    set_prepulse!(tm, P, 0.0)
    @test tm.prepulse[] == (P, 0.0, -4.0, 2.5) && history_stale(sess)
    @test_throws ArgumentError set_prepulse!(tm, P, 0.0; start = 3.0, stop = 2.0)
    @test_throws ArgumentError set_prepulse!(tm, -1.0, 0.0)
    start_history!(sess)
    h2 = tm.history[]
    @test !history_stale(sess) && h2.time == h1.time

    # the drive history: the core's with the same ratios and one prepulse for every beam
    pre = PulseShape([-4.0, 2.5], [P, P])
    uc = (; unconverted = (omega1 = 0.05, omega2 = 0.02), colors = nif_unconverted())
    ref = _core_history(sess, h1.time; uc..., prepulse = (omega1 = pre,))
    @test all(_near.(_history_fields(h2.history), _history_fields(ref)))
    @test all(_near(h2.history.unconverted.delivered[c], ref.unconverted.delivered[c]) for c in (:omega1, :omega2))
    # (more 1ω light inside the window than without the prepulse, the same after it)
    d1, d2 = h1.history.unconverted.delivered[:omega1], ref.unconverted.delivered[:omega1]
    inside = h1.time .<= 2.5
    @test any(inside) && all(d2[inside] .> d1[inside]) && _near(d2[.!inside], d1[.!inside])

    # the irradiance histories and checks
    thr = DEFAULT_THRESHOLD
    start_irradiance!(sess)
    r = tm.irradiance[]
    @test r !== nothing && !irradiance_stale(sess) && tm.chart[] === :irradiance && r.foil === nothing
    src = tm.pulses[]
    ts = r.history.time
    # sampled at the pulses' times and knots, and on both sides of each end
    edges = (prevfloat(-4.0), -4.0, 2.5, nextfloat(2.5), prevfloat(src.shape.t[1]), nextfloat(src.shape.t[end]))
    @test issorted(ts) && allunique(ts) && all(in(ts), src.times) && all(in(ts), src.nodes) && all(in(ts), edges)
    # the core's history at those times: every beam on one pulse (the same table power)
    s = dr.scenario[]
    Pt = beams[1].power
    @test all(b -> b.power == Pt, beams)
    pulse = τ -> Pt * src.shape(τ)
    hc = unconverted_history(s.target, beams, pulse, ts; ratio = uc.unconverted, prepulse = (omega1 = pre,))
    @test all(_near(r.history.element[c], hc.element[c]) && _near(r.history.backside[c], hc.backside[c])
              for c in (:omega1, :omega2))
    sc = specular_check(hc; threshold = thr)
    @test [x.part for x in r.specular] == [x.part for x in sc] && [x.specular for x in r.specular] == [x.specular for x in sc]
    @test _near([x.peak for x in r.specular], [x.peak for x in sc]) &&
          all(isapprox.([x.crossing for x in r.specular], [x.crossing for x in sc]; rtol = 1e-12))
    # 2 GW of 1ω per beam is above the threshold on the hohlraum: the first crossing at the
    # prepulse's start
    c = threshold_crossing(hc, thr)
    @test r.crossing.part === c.part === :hohlraum && isapprox(r.crossing.time, c.time; rtol = 1e-12) &&
          isapprox(c.time, -4.0; atol = 1e-12)
    rows = irradiance_checks(sess)
    @test length(rows) == 1 && rows[1].check === :specular && !rows[1].pass && occursin("-4.00 ns", rows[1].text)
    # a threshold above every peak: the checks follow from the histories kept
    set_irradiance_threshold!(tm, 2 * maximum(x.peak for x in sc))
    @test tm.irradiance[].history === r.history && all(x -> x.specular, tm.irradiance[].specular) &&
          irradiance_checks(sess)[1].pass && isinf(tm.irradiance[].crossing.time)
    @test_throws ArgumentError set_irradiance_threshold!(tm, 0.0)
    set_irradiance_threshold!(tm, thr)

    # the chart: the threshold, and the hohlraum's peak at each time (the largest 1ω + 2ω
    # irradiance on its elements and back sides, from the core's history), floored for the log axis
    ch = time_chart(sess, :irradiance)
    @test ch.scale === :log && isempty(ch.empty)
    line = only(filter(x -> startswith(x.label, "threshold"), ch.left.series))
    @test line.y == [thr, thr] && line.style === :dash
    t = s.target
    E = hc.element[:omega1] .+ hc.element[:omega2]
    B = hc.backside[:omega1] .+ hc.backside[:omega2]
    on_part(ci) = t.components[ci].part === :hohlraum && t.components[ci].kind === :surface
    els = [e for e in 1:nelements(t) if on_part(t.elem_comp[e])]
    backs = [ci for ci in eachindex(t.components) if on_part(ci) && !t.components[ci].two_sided]
    peak = vec(maximum(vcat(E[els, :], B[backs, :]); dims = 1))
    hs = only(filter(x -> x.label == "hohlraum", ch.left.series))
    @test hs.x == ts && _near(hs.y, max.(peak, 1e-4 * thr)) && maximum(peak) > thr

    # a prepulse edit makes the checks stale; back to it, current again
    set_prepulse!(tm, P, 0.0; stop = 3.0)
    @test irradiance_stale(sess) && occursin("out of date", time_chart(sess, :irradiance).title)
    set_prepulse!(tm, P, 0.0; stop = 2.5)
    @test !irradiance_stale(sess)
end

@testset "time panel: the pre-plasma check on a foil, sessions (no GL)" begin
    f = DesignSession(starter = :foil, resolution = :coarse, preset = :nif_quads)
    tm = f.time
    use_pulse_shape!(f; step = 1.0)
    beams = f.beams.beams[]
    set_unconverted_ratio!(f.drive, 0.05, 0.0)
    set_participation!(f.beams, collect(eachindex(beams)), 1, true)
    P = 2e-3
    set_prepulse!(tm, P, 0.0; start = -4.0, stop = 0.5)
    start_irradiance!(f)
    r = tm.irradiance[]
    @test r.foil === active_foil(f.target) === :foil
    # the core's histories at the same times: the unconverted light, and the 3ω light on
    # the drive's deposition operator
    s = f.drive.scenario[]
    ts = r.history.time
    shape = tm.pulses[].shape
    Pt = beams[1].power
    @test all(b -> b.power == Pt, beams)
    pulse = τ -> Pt * shape(τ)
    hc = unconverted_history(s.target, beams, pulse, ts; ratio = (omega1 = 0.05, omega2 = 0.0),
                             prepulse = (omega1 = PulseShape([-4.0, 0.5], [P, P]),))
    ac = irradiance_history(s.deposition, beams, pulse, ts)
    pc = preplasma_check(hc, ac; foil = :foil)
    p = r.preplasma
    @test p.preplasma == pc.preplasma && isapprox(p.onset, pc.onset; rtol = 1e-12) &&
          isapprox(p.arrival, pc.arrival; rtol = 1e-12)
    # the prepulse reaches the threshold on the foil at its start; the 3ω light, zero
    # before the shape's first knot, not before it
    @test pc.preplasma && isapprox(pc.onset, -4.0; atol = 1e-12) && pc.arrival >= prevfloat(shape.t[1]) &&
          pc.lead ≈ pc.arrival - pc.onset
    rows = irradiance_checks(f)
    @test [x.check for x in rows] == [:specular, :preplasma] && !rows[2].pass && occursin("foil", rows[2].text)
    # the chart's 3ω line: the foil's peak 3ω irradiance at each time (the core's)
    t = s.target
    els = [e for e in 1:nelements(t) if t.components[t.elem_comp[e]].part === :foil]
    peak3 = vec(maximum(ac.element[:omega3][els, :]; dims = 1))
    line = only(filter(x -> x.label == "3ω on foil", time_chart(f, :irradiance).left.series))
    @test line.style === :dot && _near(line.y, max.(peak3, 1e-4 * tm.threshold[]))

    # sessions keep the prepulse and the threshold; older ones load the defaults
    set_irradiance_threshold!(tm, 3e8)
    st = session_state(f)
    @test st["time"]["prepulse"] == [P, 0.0, -4.0, 0.5] && st["time"]["threshold"] == 3e8
    set_prepulse!(tm, 0.0, 0.0; start = -5.0, stop = 0.0)
    set_irradiance_threshold!(tm, 1e9)
    restore!(f, st)
    @test tm.prepulse[] == (P, 0.0, -4.0, 0.5) && tm.threshold[] == 3e8
    old = deepcopy(st)
    delete!(old["time"], "prepulse")
    delete!(old["time"], "threshold")
    restore!(f, old)
    @test tm.prepulse[] == DEFAULT_PREPULSE && tm.threshold[] == DEFAULT_THRESHOLD
end
