# Halfraum ends and packages (assemblies phase C): a hohlraum part's open ends and
# per-end entrance holes against core `hohlraum(; top_end, bottom_end)`, the walls and
# enclosure of a halfraum, the flat observed face (its face map and flux statistics, no
# capsule modes; its drive history and chart, the retune's refusal), the two halfraum
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
    @test !Controllers._closed_enclosure(ts)
    # the bottom entrance hole of its own radius (0: the LEH radius)
    @test edit_part!(ts, i, :bottom_end, "1") && edit_part!(ts, i, :bottom_leh_radius, "1.0")
    t = ts.target[]
    @test sum(t.elem_area[elements(t, :hohlraum_leh_bottom)]) ≈ π rtol = 1e-12
    @test sum(t.elem_area[elements(t, :hohlraum_leh_top)]) ≈ π * 1.55^2 rtol = 1e-12
    @test Controllers._closed_enclosure(ts)
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
    @test observed_name(ts) === :package_top && observed_flat(ts) && !Controllers._closed_enclosure(ts)
    t = ts.target[]
    # the faces close each halfraum's open end: they lie in its plane, with its radius
    zs(n) = extrema(v[3] for c in t.components if c.name === n for v in c.mesh.vertices)
    @test all(isapprox.(zs(:package_top), 1.0; atol = 1e-12)) && zs(:halfraum_top_wall)[1] ≈ 1.0
    @test all(isapprox.(zs(:package_bottom), -1.0; atol = 1e-12)) && zs(:halfraum_bottom_wall)[2] ≈ -1.0
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
    # the retune restores capsule modes: refused on a flat face, before reading the table
    retune_to_reference!(sess, joinpath(mktempdir(), "absent.csv"), joinpath(mktempdir(), "out.csv"))
    @test !sess.drive.retune_running[] && sess.drive.retune_status[] == Controllers._RETUNE_FLAT &&
          sess.time.retune[] === nothing
    @test_throws ArgumentError retune_pulses(s, beams, enabled, Dict{Symbol,PulseShape}())
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
