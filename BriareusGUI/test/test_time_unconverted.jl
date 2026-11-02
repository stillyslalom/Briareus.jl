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
    # a color on a disabled beam only does not count; on some enabled beams only, no history
    bt.enabled[] = [i != 1 for i in all_beams]
    set_participation!(bt, [1], 1, false)
    @test Controllers._history_unconverted(sess) == (omega1 = 0.05, omega2 = 0.02)
    bt.enabled[] = fill(true, length(beams))
    @test history_stale(sess)
    start_history!(sess)
    @test startswith(tm.status[], "no drive history: the 1ω light takes part on some of the beams only")
    @test tm.history[] === h1
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
    set_unconverted_ratio!(s.drive, 0.05, 0.02)
    set_participation!(s.beams, all_beams[2:end], 1, true)        # 1ω on all but one beam
    retune_to_reference!(s, table, out; cap = 1.5)
    @test startswith(s.drive.retune_status[], "no retune: the 1ω light takes part on some") && !isfile(out)
    set_participation!(s.beams, all_beams, 1, true)
    set_participation!(s.beams, all_beams, 2, true)
    s.beams.enabled[] = [nif_bundle(b) != "31" for b in beams]          # bundle 31 down
    retune_to_reference!(s, table, out; cap = 1.5)
    @test isfile(out)
    r = s.time.retune[]
    @test r.reference.unconverted !== nothing && r.after.unconverted !== nothing

    # the retuned drive: a core history of the kept beams on the new pulses, with the
    # same unconverted light
    sc = s.drive.scenario[]
    keep = findall(s.beams.enabled[])
    kept = beams[keep]
    back = r.pulses
    observe = Controllers._observed_parts(sc)
    names = component_names(sc.target)
    walls = Tuple(n for n in wall_names(s.target) if n in names)
    wall = isempty(walls) ? observe : walls
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
end
