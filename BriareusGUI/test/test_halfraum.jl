# Halfraum ends and packages (assemblies phase C): a hohlraum part's open ends and
# per-end entrance holes against core `hohlraum(; top_end, bottom_end)`, the walls and
# enclosure of a halfraum, the flat observed face (its face map and flux statistics, no
# capsule modes; its drive history and chart, its retune), the two halfraum
# starters, and the hohlraum presets aimed in their entrance holes.

using LinearAlgebra: norm

@testset "halfraum: open ends and per-end holes against core (no GL)" begin
    ts = TargetSetup(resolution = :coarse)
    h = ts.mesh[].h
    i = part_index(ts, :hohlraum)
    @test edit_part!(ts, i, :bottom_end, "2")
    p = ts.parts[][i]
    @test hohlraum_pieces(p) == [:wall, :cap_top, :leh_top]
    ref = hohlraum(; length = 10.1, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.0,
                   wall = gold(albedo = 0.85, conversion = 0.8), h, leh_open_radius = 0.84 * 1.55,
                   bottom_end = :open)
    cs = part_components(p, h)
    @test [c.name for c in cs] == [:hohlraum_wall, :hohlraum_cap_top, :hohlraum_leh_top]
    @test all(_same_component(a, b) for (a, b) in zip(cs, ref.components)) && !ref.closed
    # an open end has no cap to hold the wall role, and the enclosure is not closed
    @test wall_names(ts) == [:hohlraum_wall, :hohlraum_cap_top]
    @test !closed_design(ts)
    # the bottom entrance hole of its own radius (0: the LEH radius)
    @test edit_part!(ts, i, :bottom_end, "1") && edit_part!(ts, i, :bottom_leh_radius, "1.0")
    t = ts.target[]
    @test sum(t.elem_area[elements(t, :hohlraum_leh_bottom)]) ≈ π rtol = 1e-12
    @test sum(t.elem_area[elements(t, :hohlraum_leh_top)]) ≈ π * 1.55^2 rtol = 1e-12
    @test closed_design(ts)
    @test !edit_part!(ts, i, :bottom_end, "3") && occursin("open (2)", ts.error[])
    @test !edit_part!(ts, i, :bottom_leh_radius, "3.0") && occursin("bottom_end", ts.error[])
end

@testset "halfraum starters and the flat observed face (no GL)" begin
    # the shock tube: two halfraums facing each other across the package, its top face
    # observed and flat
    ps = starter_parts(:shock_tube)
    @test [q.name for q in ps] == [:halfraum_top, :halfraum_bottom, :package, :package_top, :package_bottom]
    @test [hohlraum_pieces(q) for q in ps[1:2]] == [[:wall, :cap_top, :leh_top], [:wall, :cap_bottom, :leh_bottom]]
    sess = DesignSession(starter = :shock_tube, resolution = :coarse)
    ts = sess.target
    @test observed_name(ts) === :package_top && observed_flat(ts) && closed_design(ts)   # two cavities
    t = ts.target[]
    # the faces close each halfraum's open end: they lie in its plane, with its radius
    zs(n) = extrema(v[3] for c in t.components if c.name === n for v in c.mesh.vertices)
    @test all(isapprox.(zs(:package_top), 1.0; atol = 1e-12)) && zs(:halfraum_top_wall)[1] ≈ 1.0
    @test all(isapprox.(zs(:package_bottom), -1.0; atol = 1e-12)) && zs(:halfraum_bottom_wall)[2] ≈ -1.0
    # and meet the walls' end rings vertex for vertex, at the coarse and standard element
    # sizes (no gap at the seam: the disks take the wall's azimuthal count)
    function ring(t, n, z)
        vs = [v for c in t.components if c.name === n for v in c.mesh.vertices
              if abs(v[3] - z) < 1e-9 && abs(hypot(v[1], v[2]) - 1.5) < 1e-9]
        az = [mod(atan(v[2], v[1]), 2π) for v in vs]
        return sort(unique(round.([a > 2π - 1e-6 ? 0.0 : a for a in az]; digits = 6)))
    end
    seams(t) = [ring(t, :halfraum_top_wall, 1.0) == ring(t, :package_top, 1.0),
                ring(t, :halfraum_bottom_wall, -1.0) == ring(t, :package_bottom, -1.0)]
    @test all(seams(t))
    st = TargetSetup(starter = :shock_tube, resolution = :standard)
    @test all(seams(st.target[])) && length(ring(st.target[], :package_top, 1.0)) == 40
    build_scenario!(sess)
    res = sess.drive.result[]
    # a flat observed face: no capsule modes; its mean flux and its statistics over the face
    @test isempty(res.modes) && res.face !== nothing
    @test res.face.low <= to_Wcm2(res.mean_flux) <= res.face.peak && 0 <= res.face.rms < 1
    s = sess.drive.scenario[]
    H = res.solution.H[s.observed]
    A = s.target.elem_area[s.observed]
    @test res.mean_flux ≈ sum(H .* A) / sum(A) rtol = 1e-12
    @test occursin("observed face", drive_summary(res)) && occursin("rms", drive_summary(res; short = true))
    mp = flux_map(sess.drive)
    @test mp.kind === :face && size(mp.values) == (length(mp.x), length(mp.y)) && all(isfinite, mp.values)
    @test extrema(mp.values) == extrema(to_Wcm2.(H))
    # the capsule symmetry objective needs a sphere
    @test_throws ArgumentError Controllers._make_objective(Controllers._objective_settings(sess.optimize; observed = :package_top), s, nothing, 0.5)

    # the opacity halfraum: the sample in the open end's plane, facing into the cavity;
    # what misses it leaves (the scenario's view factors are not balanced)
    sess = DesignSession(starter = :opacity, resolution = :coarse)
    ts = sess.target
    @test [q.name for q in ts.parts[]] == [:halfraum, :sample] && observed_name(ts) === :sample && observed_flat(ts)
    t = ts.target[]
    es = elements(t, :sample)
    @test all(e -> t.elem_normal[e] ≈ Briareus.Vec3(0, 0, 1) && abs(t.elem_centroid[e][3] + 1.5) < 1e-12, es)
    build_scenario!(sess)
    res = sess.drive.result[]
    @test isempty(res.modes) && res.face !== nothing && res.face.mean > 0
    eb = energy_balance(res.solution)
    @test eb.lost > 0
    @test sum(values(eb.absorbed)) + sum(values(eb.escaped)) + eb.lost ≈ eb.source rtol = 1e-10

    # the drive history of the flat face: per step, the face's statistics and its Tr from
    # the area-mean flux, not Legendre moments
    use_pulse_shape!(sess; step = 1.0)
    start_history!(sess)
    h = sess.time.history[]
    @test h isa DriveHistory && h.face !== nothing && size(h.legendre) == (0, length(h.time))
    # independent check: the core `drive_history` on the same scenario, pulses and wall,
    # recording each step's radiosity solution's incident flux, reduced here directly
    s = sess.drive.scenario[]
    src, beams, enabled = sess.time.pulses[], sess.beams.beams[], sess.beams.enabled[]
    observe = Controllers._observed_parts(s)
    walls = Tuple(n for n in wall_names(ts) if n in component_names(s.target))
    core = drive_history(s.model, s.deposition, τ -> Controllers._pulse_beams(src, beams, enabled, τ), h.time;
                         observe, wall = isempty(walls) ? observe : walls, lmax = 4,
                         record = sol -> sol.H[s.observed])
    @test core.time == h.time
    A = s.target.elem_area[s.observed]
    ks = unique([1, cld(length(h.time), 2), length(h.time)])
    Hs = [core.record[k] for k in ks]
    μ = [sum(H .* A) / sum(A) for H in Hs]
    rms = [sqrt(sum(A .* (H .- m) .^ 2) / sum(A)) / m for (H, m) in zip(Hs, μ)]
    @test isapprox(h.face.mean[ks], 1e14 .* μ; rtol = 1e-12)
    @test isapprox(h.face.peak[ks], 1e14 .* maximum.(Hs); rtol = 1e-12)
    @test isapprox(h.face.low[ks], 1e14 .* minimum.(Hs); rtol = 1e-12)
    @test isapprox(h.face.rms[ks], rms; rtol = 1e-12)
    @test isapprox(h.Tr[ks], radiation_temperature.(μ); rtol = 1e-12)
    # (the core's Tr is from the Legendre moment P0 over the face: close, not the same)
    @test h.Tr != core.Tr && isapprox(h.Tr, core.Tr; rtol = 0.01)
    # the drive chart: the face's Tr, and its statistics on the right axis
    ch = time_chart(sess, :drive)
    @test isempty(ch.empty) && startswith(ch.title, "Observed face") && ch.left.series[1].label == "Face Tr"
    @test ch.left.series[1].y == h.Tr && ch.right.label == "Face nonuniformity (%)"
    @test [q.label for q in ch.right.series] == ["RMS", "peak", "lowest"]
    @test ch.right.series[1].y == 100 .* h.face.rms
    @test ch.right.series[2].y ≈ 100 .* (h.face.peak ./ h.face.mean .- 1) rtol = 1e-12
    @test all(ch.right.series[3].y .<= 0 .<= ch.right.series[2].y)
    # the retune of the flat face: the first quad down, the others retuned to the
    # reference's flux over the face
    dir = mktempdir()
    pul = Dict(g => PulseShape([0.0, 1.0, 2.0, 3.0], [0.0, 1.0, 2.0, 0.0]) for g in unique(b.group for b in beams))
    table = write_pulse_table(joinpath(dir, "reference.csv"), pul)
    out = joinpath(dir, "retuned.csv")
    sess.beams.enabled[] = [b.quad != beams[1].quad for b in beams]
    retune_to_reference!(sess, table, out; cap = 1.5)
    r = sess.time.retune[]
    @test r.kind === :face && isfile(out) && occursin("face worst RMS", sess.drive.retune_status[])
    # independent check: core histories of the enabled beams on the reference and on the
    # retuned pulses (as written, before the CSV's rounding), the face's mean and RMS from
    # each step's flux on its elements
    kept = findall(sess.beams.enabled[])
    stats(P, by) = begin
        hc = drive_history(s.model, s.deposition[kept], τ -> with_pulses(beams[kept], P, τ; by), r.after.time;
                           observe, wall = isempty(walls) ? observe : walls, lmax = 4,
                           record = sol -> sol.H[s.observed])
        μ = [sum(H .* A) / sum(A) for H in hc.record]
        (; mean = 1e14 .* μ, rms = [sqrt(sum(A .* (H .- m) .^ 2) / sum(A)) / m for (H, m) in zip(hc.record, μ)])
    end
    after, before = stats(r.pulses, :quad), stats(pul, :group)
    @test isapprox(read_pulse_table(out)[first(keys(r.pulses))].p, first(values(r.pulses)).p; rtol = 1e-6)
    @test isapprox(r.after.mean, after.mean; rtol = 1e-12) && isapprox(r.after.rms, after.rms; rtol = 1e-12)
    @test isapprox(r.before.mean, before.mean; rtol = 1e-12) && isapprox(r.before.rms, before.rms; rtol = 1e-12)
    # the retune restores the mean and moves the RMS toward the reference's
    Δ(f) = maximum(abs.(f.mean ./ r.reference.mean .- 1))
    @test Δ(r.after) < 1e-3 < Δ(r.before)
    @test maximum(abs.(r.after.rms .- r.reference.rms)) < maximum(abs.(r.before.rms .- r.reference.rms))
    # its chart: the face's RMS and mean, before and after
    ch = time_chart(sess, :retune)
    @test ch.right.label == "Face RMS, mean − reference (%)"
    @test [x.label for x in ch.right.series] == ["RMS before", "RMS retuned", "RMS reference", "mean before",
                                                 "mean retuned"]
    @test ch.right.series[2].y == 100 .* r.after.rms && ch.right.series[5].y ≈ 100 .* (r.after.mean ./ r.reference.mean .- 1)
    # nothing to hold (a foil asked for, and the target has none): refused
    @test Controllers._retune_kind(s, nothing, :foil) === nothing
    @test_throws ArgumentError retune_pulses(s, beams, enabled, pul; on = :foil)
end

@testset "hohlraum presets aim in the starters' entrance holes (no GL)" begin
    V = Briareus.Vec3
    top, bot = Symbol("hohlraum.leh_top"), Symbol("hohlraum.leh_bottom")
    upper(bs) = [b.port[3] > 0 for b in bs]
    # the default design: the quads in its entrance holes' frames, their chamber aims the
    # example cone pointing's own, bit for bit (`isequal` tells -0.0 from 0.0)
    cone = Controllers._cone_pointing(nif_quads(power = 8.0), Controllers.NIF_CONE_POINTING)
    s = DesignSession(resolution = :coarse)
    bt = s.beams
    up = upper(bt.beams[])
    @test all(isequal(a.aim, b.aim) for (a, b) in zip(bt.beams[], cone)) && bt.beams[] == cone
    @test bt.aim_frame[] == [u ? top : bot for u in up] && all(bt.enabled[])
    # the drive equals the drive of the same beams held at TCC (the former preset), bit for bit
    s0 = DesignSession(resolution = :coarse, beams = cone)
    @test all(==(:tcc), s0.beams.aim_frame[])
    build_scenario!(s)
    build_scenario!(s0)
    r, r0 = s.drive.result[], s0.drive.result[]
    @test r.modes == r0.modes && r.mean_flux == r0.mean_flux && r.solution.H == r0.solution.H &&
          [x.power for x in r.power] == [x.power for x in r0.power]
    # the OMEGA hohlraum beams at the entrance holes' centers, as before (±length/2 on z)
    om = preset_pointing(:omega_hohlraum; target = s.target)
    old = Controllers._cone_pointing(omega_beams(power = 0.5), Dict(g => (0.0, 10.1 / 2) for g in (:cone1, :cone2, :cone3)))
    @test om.beams == old && om.aim_frame == [u ? top : bot for u in upper(om.beams)]
    @test all(==(:tcc), preset_pointing(:nif_quads).aim_frame)   # no target: at TCC

    # the hohlraum moved: its beams by the same displacement
    loc = copy(bt.aim_local[])
    aims = [b.aim for b in bt.beams[]]
    i = part_index(s.target, :hohlraum)
    @test nudge_part!(s.target, i, 0.4; frame = :tcc, axis = :x) && nudge_part!(s.target, i, -0.25; frame = :tcc, axis = :z)
    d = V(0.4, 0, -0.25)
    @test maximum(norm(b.aim - (a + d)) for (a, b) in zip(aims, bt.beams[])) < 1e-12 && bt.aim_local[] == loc

    # the shock tube: upper beams in the upper halfraum's hole, lower in the lower one's,
    # at the default design's offsets in those frames
    st = DesignSession(starter = :shock_tube, resolution = :coarse)
    sb = st.beams
    ft, fb = Symbol("halfraum_top.leh_top"), Symbol("halfraum_bottom.leh_bottom")
    @test facing_lehs(st.target) == (ft, fb)
    @test sb.aim_frame[] == [u ? ft : fb for u in upper(sb.beams[])] && all(sb.enabled[])
    @test maximum(norm.(sb.aim_local[] .- loc)) < 1e-12
    fr = part_frames(st.target)
    @test fr[ft].origin == V(0, 0, 3.5) && fr[fb].origin == V(0, 0, -3.5)
    # the upper aims lie at the default offsets from the upper hole's plane (z = 3.5 here,
    # 5.05 in the default hohlraum), the lower ones mirrored about the lower hole's
    @test all(abs((b.aim[3] - 3.5) - (c.aim[3] - 10.1 / 2)) < 1e-12 && norm(b.aim[1:2] - c.aim[1:2]) < 1e-12
              for (b, c, u) in zip(sb.beams[], cone, upper(cone)) if u)
    @test all(abs((b.aim[3] + 3.5) - (c.aim[3] + 10.1 / 2)) < 1e-12 for (b, c, u) in zip(sb.beams[], cone, upper(cone)) if !u)
    # moving the upper halfraum carries only its beams
    a0 = [b.aim for b in sb.beams[]]
    @test nudge_part!(st.target, part_index(st.target, :halfraum_top), 0.3; frame = :tcc, axis = :z)
    moved = [norm(b.aim - a) for (a, b) in zip(a0, sb.beams[])]
    @test all(abs(m - (u ? 0.3 : 0.0)) < 1e-12 for (m, u) in zip(moved, upper(sb.beams[])))
    # another preset: the OMEGA beams at each hole's center
    load_preset!(st, :omega_hohlraum)
    fr = part_frames(st.target)
    @test all(b.aim == fr[f].origin for (b, f) in zip(sb.beams[], sb.aim_frame[])) &&
          sb.aim_frame[] == [u ? ft : fb for u in upper(sb.beams[])]

    # the opacity halfraum: one hole, facing up; the lower beams are disabled at TCC
    op = DesignSession(starter = :opacity, resolution = :coarse)
    ob = op.beams
    u = upper(ob.beams[])
    @test facing_lehs(op.target) == (Symbol("halfraum.leh_top"), nothing)
    @test ob.enabled[] == u && all(ob.aim_frame[][u] .=== Symbol("halfraum.leh_top")) && all(ob.aim_frame[][.!u] .=== :tcc)
    @test [b.aim for b in ob.beams[][.!u]] == [b.aim for b in cone[.!u]]
    @test maximum(norm.(ob.aim_local[][u] .- loc[u])) < 1e-12

    # inserted into a design whose beams drive no hohlraum: the beams at TCC go to the new
    # holes, edits kept as offsets; a second hohlraum leaves them with the first
    e = DesignSession(starter = :none, resolution = :coarse)
    eb = e.beams
    @test all(==(:tcc), eb.aim_frame[])
    edit_beam_set!(eb, [1], :aim_x, "0.2")
    @test insert_starter!(e, :shock_tube)
    @test eb.aim_frame[] == [u ? ft : fb for u in up] && all(eb.enabled[])
    want = Controllers._reaim(V(0.2, cone[1].aim[2], cone[1].aim[3]), part_frames(e.target)[eb.aim_frame[][1]],
                              Controllers._leh_reference(:nif_quads)[up[1] ? 1 : 2])
    @test eb.beams[][1].aim == want
    held = copy(eb.aim_frame[])
    @test insert_starter!(e, :hohlraum) && eb.aim_frame[] == held
    # File ▸ New target: a hohlraum preset is loaded again for holes it does not aim in
    @test new_target!(e, :hohlraum)
    @test eb.beams[] == cone && eb.aim_frame[] == [u ? top : bot for u in up]

    # sessions: the frames and local aims round trip; an older session (aims at TCC) loads
    # as saved, not aimed again
    f = joinpath(mktempdir(), "lehs.briar.toml")
    save_session(f, st)
    s2 = DesignSession(resolution = :coarse)
    load_session!(s2, f)
    @test s2.beams.aim_frame[] == sb.aim_frame[] && s2.beams.aim_local[] == sb.aim_local[] &&
          s2.beams.beams[] == sb.beams[]
    state = session_state(s0)
    @test all(r -> r["aim_frame"] == "tcc", state["beams"]["beam"])
    state["beams"]["preset"] = "nif_quads"
    s3 = DesignSession(starter = :shock_tube, resolution = :coarse)
    restore!(s3, state)
    @test all(==(:tcc), s3.beams.aim_frame[]) && s3.beams.beams[] == cone && s3.beams.preset[] === :nif_quads
end

@testset "per-cavity wall temperature (no GL)" begin
    # the default design, one cavity: its drive and drive history equal those computed
    # without cavities (the path before them), bit for bit, and the cavity's wall
    # temperature is the history's wall temperature
    sess = DesignSession(resolution = :coarse)
    ts = sess.target
    @test cavity_walls(ts) == [:hohlraum => [:hohlraum_wall, :hohlraum_cap_top, :hohlraum_cap_bottom]]
    @test observed_cavity(ts) === :hohlraum
    build_scenario!(sess)
    dr = sess.drive
    s, res = dr.scenario[], dr.result[]
    old = Controllers._evaluate_drive(s, display_beams(sess.beams), dr.lookup[], dr.spot_radius[];
                                      ratio = unconverted_ratios(sess), foil = active_foil(ts))
    @test isempty(old.cavities) && [c.name for c in res.cavities] == [:hohlraum]
    @test res.modes == old.modes && res.mean_flux == old.mean_flux && res.solution.H == old.solution.H &&
          res.solution.B == old.solution.B && [x.power for x in res.power] == [x.power for x in old.power] &&
          res.map.values == old.map.values && res.face === old.face === nothing
    @test drive_cavity(sess) === :hohlraum
    use_pulse_shape!(sess; step = 1.0)
    start_history!(sess)
    h = sess.time.history[]
    src, beams, enabled = sess.time.pulses[], sess.beams.beams[], sess.beams.enabled[]
    h0 = Controllers._history_job(s, beams, enabled, src, h.time; diags = Controllers._recorded_definitions(sess),
                                  wall = Tuple(wall_names(ts)))
    @test isempty(h0.cavities) && isempty(h0.history.cavity_Tr)
    @test all(k -> k === :cavity_Tr || isequal(h.history[k], h0.history[k]), keys(h.history))
    # (`isequal`: the laser field is NaN on apertures)
    @test isequal((h.Tr, h.legendre, h.fields, h.range, h.energy), (h0.Tr, h0.legendre, h0.fields, h0.range, h0.energy)) &&
          isequal([(d.Tr, d.intensity) for d in h.diagnostics], [(d.Tr, d.intensity) for d in h0.diagnostics])
    @test [c.name for c in h.cavities] == [:hohlraum] && h.cavities[1].Tr == h.history.wall_Tr
    # one cavity: the drive chart as before (the observed surface's Tr is the headline)
    @test !any(q -> endswith(q.label, "wall Tr"), time_chart(sess, :drive).left.series)

    # the shock tube, two cavities, half the lower beams off: each halfraum's wall temperature
    # against Σ H A / Σ A of its own wall pieces' elements (found by component name here)
    st = DesignSession(starter = :shock_tube, resolution = :coarse)
    ts = st.target
    walls = [:halfraum_top => [:halfraum_top_wall, :halfraum_top_cap_top],
             :halfraum_bottom => [:halfraum_bottom_wall, :halfraum_bottom_cap_bottom]]
    @test cavity_walls(ts) == walls
    # the observed package face closes the upper halfraum's open end
    @test observed_cavity(ts) === :halfraum_top
    st.beams.enabled[] = [b.port[3] > 0 || isodd(k) for (k, b) in enumerate(st.beams.beams[])]
    build_scenario!(st)
    r = st.drive.result[]
    t = st.drive.scenario[].target
    function direct(H, ws)
        es = [e for e in 1:nelements(t) if t.components[t.elem_comp[e]].name in ws]
        return radiation_temperature(sum(H[e] * t.elem_area[e] for e in es) / sum(t.elem_area[e] for e in es))
    end
    @test [c.name for c in r.cavities] == first.(walls)
    @test all(isapprox(c.Tr, direct(r.solution.H, last(w)); rtol = 1e-12) for (c, w) in zip(r.cavities, walls))
    @test r.cavities[1].Tr > r.cavities[2].Tr + 1     # the lit halfraum is the hotter
    @test all(c -> isapprox(c.Tr, radiation_temperature(c.flux); rtol = 1e-12), r.cavities)
    # the panel's cavity: the observed surface's, until another is chosen
    @test drive_cavity(st) === :halfraum_top
    @test set_drive_cavity!(st, :halfraum_bottom) && drive_cavity(st) === :halfraum_bottom
    @test !set_drive_cavity!(st, :package) && drive_cavity(st) === :halfraum_bottom
    # the drive history: each halfraum's wall temperature per step, against the core
    # history's recorded incident flux reduced here
    use_pulse_shape!(st; step = 1.0)
    start_history!(st)                            # Dante 1 (the default): the whole target in a 60 mm field
    hs = st.time.history[]
    s = st.drive.scenario[]
    src, beams, enabled = st.time.pulses[], st.beams.beams[], st.beams.enabled[]
    observe = Controllers._observed_parts(s)
    hole_view = dante1(field_of_view = nothing, through = :halfraum_bottom_leh_bottom)
    core = drive_history(s.model, s.deposition, τ -> Controllers._pulse_beams(src, beams, enabled, τ), hs.time;
                         observe, wall = Tuple(wall_names(ts)), lmax = 4,
                         record = sol -> (copy(sol.H), dante(hole_view, sol).intensity))
    @test [c.name for c in hs.cavities] == first.(walls)
    @test all(isapprox(c.Tr, [direct(first(r), last(w)) for r in core.record]; rtol = 1e-12)
              for (c, w) in zip(hs.cavities, walls))
    # the drive chart: a wall series per cavity
    ser = Dict(q.label => q for q in time_chart(st, :drive).left.series)
    @test ser["halfraum_top wall Tr"].y == hs.cavities[1].Tr && ser["halfraum_bottom wall Tr"].y == hs.cavities[2].Tr
    @test ser["halfraum_top wall Tr"].color != ser["halfraum_bottom wall Tr"].color
    # Dante 1 over the whole pair: from below it sees the lower halfraum through its hole
    # (the package is black, every one-sided outside dark), as the view through that hole
    whole = only(x for x in hs.diagnostics if x.id == "nif_dante1")
    @test isempty(whole.note) && whole.view.field_of_view == 60.0 && whole.view.through === :halfraum_bottom_leh_bottom
    @test all(>(0), whole.intensity) && isapprox(whole.intensity, last.(core.record); rtol = 2e-3)   # 3.5e-4
end

@testset "per-cavity balancing (no GL)" begin
    # the shock tube: each halfraum and the package face closing its open end form a
    # cavity (the package's outer face in neither), each balanced on its own
    st = DesignSession(starter = :shock_tube, resolution = :coarse)
    ts = st.target
    t = ts.target[]
    named(pre...) = sort([e for e in 1:nelements(t)
                          if any(p -> startswith(String(t.components[t.elem_comp[e]].name), p), pre)])
    cs = Controllers._closed_cavities(ts)
    @test [c.name for c in cs] == [:halfraum_top, :halfraum_bottom] && all(c -> c.balanced, cs)
    @test sort(cs[1].elements) == named("halfraum_top_", "package_top") &&
          sort(cs[2].elements) == named("halfraum_bottom_", "package_bottom")
    side = elements(t, :package)
    @test !isempty(side) && !any(e -> e in cs[1].elements || e in cs[2].elements, side)
    @test all(c -> c.leak < 1e-5, cs)        # the holes meet the caps vertex for vertex (7.5e-3, 7.8e-3 at half the count)
    build_scenario!(st)
    @test occursin("balanced per cavity: halfraum_top, halfraum_bottom", st.drive.status[])
    s = st.drive.scenario[]
    F = s.model.F
    rows(F, c) = vec(sum(F[c, c]; dims = 2))
    F0 = view_factors(s.target)                         # the raw view factors
    @test all(c -> maximum(abs.(rows(F0, c.elements) .- 1)) > 1e-5, cs)
    @test all(c -> maximum(abs.(rows(F, c.elements) .- 1)) < 1e-10, cs)
    # energy, exactly per cavity: a source in the upper halfraum's wall alone; its power
    # equals what its own elements absorb and its hole lets out, Σ (1 − α) H A over them,
    # and the lower cavity receives nothing
    P = zeros(nelements(t))
    P[elements(t, :halfraum_top_wall)] .= 1.0
    sol = solve(s.model, P)
    absorbed(es) = sum((1 - sol.albedo[e]) * sol.H[e] * t.elem_area[e] for e in es)
    @test absorbed(cs[1].elements) ≈ sum(sol.source .* t.elem_area) rtol = 1e-12
    @test absorbed(cs[2].elements) == 0

    # the default design: one cavity of every element, balanced as the whole target
    ts = TargetSetup(resolution = :coarse)
    cs = Controllers._closed_cavities(ts)
    @test length(cs) == 1 && cs[1].balanced && length(cs[1].elements) == nelements(ts.target[])
    @test Controllers._closure_kw(ts.target[], cs) === true
    # two hohlraums side by side, each with its capsule: two cavities
    @test insert_starter!(ts, :hohlraum) && edit_part!(ts, part_index(ts, :hohlraum_2), :x, "8")
    t = ts.target[]
    cs = Controllers._closed_cavities(ts)
    @test [c.name for c in cs] == [:hohlraum, :hohlraum_2] && all(c -> c.balanced, cs)
    @test sort(cs[2].elements) == sort(elements(t, :capsule_2, (Symbol(n) for n in component_names(t)
                                                                  if startswith(String(n), "hohlraum_2_"))...))
    # a piece left out leaves a hole (its cavity unbalanced, the other one kept)
    @test set_excluded!(ts, :hohlraum_2_leh_top, true)
    cs = Controllers._closed_cavities(ts)
    @test cs[1].balanced && cs[2].holed && !cs[2].balanced
    @test occursin("hohlraum_2 not balanced: a piece of its boundary", Controllers._cavities_note(cs))
    # a capsule left out leaves its cavity closed
    @test set_excluded!(ts, :hohlraum_2_leh_top, false) && set_excluded!(ts, :capsule_2, true)
    @test closed_design(ts)

    # the opacity halfraum: open where the sample does not close its end, so not balanced
    ts = TargetSetup(starter = :opacity, resolution = :coarse)
    c = only(Controllers._closed_cavities(ts))
    @test !c.balanced && c.leak > 0.05
end

@testset "piece roles: an observed cap (no GL)" begin
    sess = DesignSession(resolution = :coarse)
    ts = sess.target
    cap, walls = :hohlraum_cap_top, [:hohlraum_wall, :hohlraum_cap_bottom]
    # defaults: the wall and caps hold the wall role, the holes none
    p = ts.parts[][part_index(ts, :hohlraum)]
    @test [piece_role(p, s) for s in hohlraum_pieces(p)] == [:wall, :wall, :wall, :none, :none]
    @test piece_roles(ts, cap).roles == PIECE_ROLES && piece_roles(ts, :hohlraum_wall).roles == (:none, :wall)
    @test piece_roles(ts, :hohlraum_leh_top).roles == (:none,) && piece_roles(ts, :capsule) === nothing
    # refused: an aperture, a curved piece observed, a part's own component
    @test !set_piece_role!(ts, :hohlraum_leh_top, :wall) && occursin("aperture", ts.error[])
    @test !set_piece_role!(ts, :hohlraum_wall, :observed) && occursin("curved", ts.error[])
    @test !set_piece_role!(ts, :capsule, :wall)
    @test_throws ArgumentError set_piece_role!(ts, cap, :foil)
    target = ts.target[]
    record!(sess)
    # the top cap observed: the capsule gives the role up, the cap leaves the walls, and
    # the cavity holding it is the hohlraum
    @test set_piece_role!(ts, cap, :observed) && record!(sess)
    @test ts.target[] === target                              # the geometry is unchanged
    @test observed_name(ts) === cap && observed_flat(ts) && !ts.parts[][part_index(ts, :capsule)].observed
    @test cavity_walls(ts) == [:hohlraum => walls] && wall_names(ts) == walls
    @test observed_cavity(ts) === :hohlraum && observed_components(ts) == [cap]
    @test closed_design(ts)
    row = only(r for r in model_tree(sess) if r.id == Controllers.component_id(ts, cap))
    @test occursin("observed", row.detail)
    # the drive: the cap's face map and statistics, no modes; its mean flux and the
    # cavity's wall temperature against Σ H A / Σ A over elements found by name here
    build_scenario!(sess)
    res = sess.drive.result[]
    s = sess.drive.scenario[]
    t = s.target
    named(ns) = [e for e in 1:nelements(t) if t.components[t.elem_comp[e]].name in ns]
    mean(H, es) = sum(H[e] * t.elem_area[e] for e in es) / sum(t.elem_area[e] for e in es)
    H = res.solution.H
    @test Controllers._observed_parts(s) == (cap,) && sort(s.observed) == named((cap,))
    @test isempty(res.modes) && res.face !== nothing && flux_map(sess.drive).kind === :face
    @test res.mean_flux ≈ mean(H, named((cap,))) rtol = 1e-12
    @test [c.name for c in res.cavities] == [:hohlraum]
    @test res.cavities[1].Tr ≈ radiation_temperature(mean(H, named(walls))) rtol = 1e-12
    # (counting the cap would change it)
    @test abs(res.cavities[1].Tr / radiation_temperature(mean(H, named([walls; cap]))) - 1) > 1e-6
    # the retune holds the face; the capsule symmetry objective refuses it
    @test Controllers._retune_kind(s, nothing, :auto) === :face
    @test_throws ArgumentError Controllers._make_objective(Controllers._objective_settings(sess.optimize; observed = cap),
                                                           s, nothing, 0.5)
    # the drive history: the face's Tr per step, the cavity's from its own walls, and the
    # energy kept apart for the cap, against core histories on the same scenario
    use_pulse_shape!(sess; step = 1.0)
    start_history!(sess)
    h = sess.time.history[]
    @test h.face !== nothing && [c.name for c in h.cavities] == [:hohlraum]
    src, beams, enabled = sess.time.pulses[], sess.beams.beams[], sess.beams.enabled[]
    at = τ -> Controllers._pulse_beams(src, beams, enabled, τ)
    core = drive_history(s.model, s.deposition, at, h.time; observe = (cap,), wall = Tuple(walls), lmax = 4,
                         record = sol -> copy(sol.H))
    @test core.time == h.time
    @test isapprox(h.Tr, [radiation_temperature(mean(r, named((cap,)))) for r in core.record]; rtol = 1e-12)
    @test isapprox(h.cavities[1].Tr, [radiation_temperature(mean(r, named(walls))) for r in core.record]; rtol = 1e-12)
    eb = energy_budget(h.history)
    byp = energy_budget(core)                       # by part: the cap summed into the hohlraum
    @test haskey(eb.heating, cap) && !haskey(byp.heating, cap)
    @test sum(values(eb.heating)) ≈ sum(values(byp.heating)) rtol = 1e-12
    @test sum(v for (k, v) in eb.heating if startswith(String(k), "hohlraum_")) ≈ byp.heating[:hohlraum] rtol = 1e-12
    @test only(x for x in h.energy if x.name === :capsule).values[end] ≈ eb.heating[cap] rtol = 1e-12
    # undo takes back the loaded pulses, then gives the capsule its role back; redo the cap's
    @test undo!(sess) && sess.time.pulses[] === nothing && observed_name(ts) === cap
    @test undo!(sess) && observed_name(ts) === :capsule
    @test cavity_walls(ts) == [:hohlraum => [:hohlraum_wall, cap, :hohlraum_cap_bottom]]
    @test redo!(sess) && observed_name(ts) === cap
    # sessions: the role round trips (version 11); a version-10 session loads with the
    # pieces' default roles
    f = joinpath(mktempdir(), "cap.briar.toml")
    save_session(f, sess)
    saved = TOML.parsefile(f)
    hp = only(q for q in saved["target"]["parts"] if q["name"] == "hohlraum")
    @test saved["version"] == Controllers.SESSION_VERSION && hp["piece_roles"] == Dict("cap_top" => "observed")
    s2 = DesignSession(resolution = :coarse)
    load_session!(s2, f)
    @test observed_name(s2.target) === cap && cavity_walls(s2.target) == [:hohlraum => walls]
    old = session_state(sess)
    old["version"] = 10
    foreach(q -> delete!(q, "piece_roles"), old["target"]["parts"])
    old["target"]["parts"][findfirst(q -> q["name"] == "capsule", old["target"]["parts"])]["observed"] = true
    restore!(s2, old)
    @test observed_name(s2.target) === :capsule && all(q -> isempty(q.piece_roles), s2.target.parts[])
    # the capsule given the role again: the cap gives it up and is a wall again
    @test set_part_role!(ts, part_index(ts, :capsule), :observed)
    @test observed_name(ts) === :capsule && isempty(ts.parts[][part_index(ts, :hohlraum)].piece_roles)
    # a cap given no role leaves the walls; back to its default, the override goes
    @test set_piece_role!(ts, :hohlraum_cap_bottom, :none) && wall_names(ts) == [:hohlraum_wall, cap]
    bottom = only(r for r in model_tree(sess) if r.id == Controllers.component_id(ts, :hohlraum_cap_bottom))
    @test occursin("no role", bottom.detail)
    @test set_piece_role!(ts, :hohlraum_cap_bottom, :wall) && isempty(ts.parts[][part_index(ts, :hohlraum)].piece_roles)
end

@testset "the face objective on the opacity halfraum (no GL)" begin
    sess = DesignSession(starter = :opacity, resolution = :coarse)
    build_scenario!(sess)
    st = sess.optimize
    # the capsule objectives refuse a flat face, and say which objective serves it
    start_optimize!(sess; spawn = false)
    @test st.result[] === nothing && occursin("Observed face uniformity", st.status[])
    set_objective!(st, :face)
    @test edit_optimize!(st, :uniformity, "1e6") && st.uniformity[] == 1e6
    @test !edit_optimize!(st, :uniformity, "-1")
    start_optimize!(sess; spawn = false)
    r = st.result[]
    @test r isa OptimizeResult && isempty(r.after.modes)
    # the face's RMS before and after against a radiosity solve of the optimized beams
    s = sess.drive.scenario[]
    sub = subset(s, r.indices)
    rms(beams) = (ev = evaluate(sub, beams); H = solve(sub.model, ev.deposition).H[sub.observed];
                  rms_nonuniformity(H, sub.target.elem_area[sub.observed]))
    @test r.before.face.rms ≈ rms(sess.beams.beams[][r.indices]) rtol = 1e-9
    @test r.after.face.rms ≈ rms(r.beams) rtol = 1e-9
    @test r.after.face.rms < 0.5 * r.before.face.rms     # measured 7.3e-5 → 8.9e-6
    # the setting round trips in a session; an older one loads with the default weight
    state = session_state(sess)
    @test state["optimize"]["uniformity"] == 1e6 && state["optimize"]["objective"] == "face"
    delete!(state["optimize"], "uniformity")
    restore!(sess, state)
    @test st.uniformity[] == 1.0
end
