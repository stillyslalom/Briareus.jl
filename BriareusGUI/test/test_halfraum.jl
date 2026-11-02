# Halfraum ends and packages (assemblies phase C): a hohlraum part's open ends and
# per-end entrance holes against core `hohlraum(; top_end, bottom_end)`, the walls and
# enclosure of a halfraum, the flat observed face (its face map and flux statistics, no
# capsule modes), and the two halfraum starters.

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
end
