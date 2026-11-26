# This file tests participation and visibility: the beam matrix (1ω, 2ω and 3ω light in
# the calculations and drawn), component display modes, the exact clipping of the cutaway,
# their session state, and the unconverted-light drive (no GL).

using LinearAlgebra: norm, cross, dot, normalize

_area(a, b, c) = norm(cross(b - a, c - a)) / 2
_tris_area(ts) = sum((_area(t...) for t in ts); init = 0.0)

@testset "cutaway geometry (no GL)" begin
    v = ViewState(cut_components = (:hohlraum_wall,))
    # the square plate [-1, 1]² at z = 0, as two triangles
    P(x, y) = Briareus.Vec3(x, y, 0.0)
    sq = [(P(-1, -1), P(1, -1), P(1, 1)), (P(-1, -1), P(1, 1), P(-1, 1))]
    kept(v; name = :hohlraum_wall) = (r = cut_region(v); sum(_tris_area(cut_triangle(r, name, t...)[1]) for t in sq))
    # the default quadrant facing −x −y; another component is not cut
    @test kept(v) ≈ 3 rtol = 1e-12
    @test kept(v; name = :capsule) ≈ 4 rtol = 1e-12
    # a 60° sector about +x: the wedge |y| < x tan 30° (closed form 4 − tan 30°)
    set_cut!(v; azimuth = 0, width = 60)
    @test kept(v) ≈ 4 - tand(30) rtol = 1e-12
    # a 300° sector centered at 10° keeps 160° to 220°: ½ (tan 40° + tan 20°)
    set_cut!(v; azimuth = 10, width = 300)
    @test kept(v) ≈ (tand(40) + tand(20)) / 2 rtol = 1e-12
    # A cut bounded along its axis, a half about y, removes x > 0 within |y| < 0.5, a
    # 1 × 1 square. The azimuth runs from +z toward +x, so 90° is +x.
    set_cut!(v; axis = :y, azimuth = 90, width = 180, range = (-0.5, 0.5))
    @test kept(v) ≈ 3 rtol = 1e-12
    set_cut!(v; axis = :z, range = (-Inf, Inf))
    # The plane x + y > √2/2 removes the corner triangle of legs 2 − √2/2.
    set_cut_plane!(v; normal = (1, 1, 0), point = (0, 0, 0))
    @test v.cut_kind[] === :plane
    @test kept(v) ≈ 2 rtol = 1e-12
    set_cut_plane!(v; offset = 0.5)
    @test cut_offset(v) ≈ 0.5
    @test kept(v) ≈ 4 - (2 - 0.5 * sqrt(2))^2 / 2 rtol = 1e-12
    # the section edge of a plane cut lies on the plane, its length the chord
    r = cut_region(v)
    secs = reduce(vcat, (cut_triangle(r, :hohlraum_wall, t...)[3] for t in sq))
    @test all(s -> all(p -> abs(dot(Briareus.Vec3(1, 1, 0) / sqrt(2), p) - 0.5) < 1e-12, s), secs)
    @test sum(s -> norm(s[2] - s[1]), secs) ≈ sqrt(2) * (2 - 0.5 * sqrt(2)) rtol = 1e-12
    # every cell of any split is wholly kept or removed; the cells tile the triangle
    set_cut!(v; kind = :sector, azimuth = 37, width = 75, range = (-0.3, 0.8))
    r = cut_region(v)
    tri = (Briareus.Vec3(-1.2, 0.3, -0.7), Briareus.Vec3(1.1, 0.9, 1.3), Briareus.Vec3(0.4, -1.5, 0.2))
    cells = cut_pieces(r, tri...)
    @test length(cells) > 1
    @test sum(c -> sum(_area(c.polygon[1], c.polygon[i], c.polygon[i+1]) for i in 2:length(c.polygon)-1),
              cells) ≈ _area(tri...) rtol = 1e-12
    # (each cell's corners, nudged into it, fall on its side of the cut)
    inner(c) = (m = sum(c.polygon) / length(c.polygon); [p + 1e-6 * (m - p) for p in c.polygon])
    @test all(c -> all(p -> is_cut(r, :hohlraum_wall, p) == c.removed, inner(c)), cells)
    @test any(c -> c.removed, cells) && !all(c -> c.removed, cells)
    @test_throws ArgumentError set_cut!(v; range = (1.0, 0.0))
    @test_throws ArgumentError set_cut_plane!(v; normal = (0, 0, 0))
end

@testset "arbitrary cut planes (no GL)" begin
    V(x...) = Briareus.Vec3(x...)
    v = ViewState(cut_components = (:hohlraum_wall,))
    # A normal typed as port angles takes the polar angle from +z and the azimuth from +x
    # toward +y.
    set_cut_plane!(v; angles = (60, 135), point = (0.5, 0.0, -1.0), offset = 0.25)
    @test all(isapprox.(v.cut_normal[], (sind(60) * cosd(135), sind(60) * sind(135), cosd(60)); atol = 1e-15))
    @test all(isapprox.(cut_angles(v), (60, 135); rtol = 1e-12))
    @test cut_offset(v) == 0.25 && v.cut_point[] == (0.5, 0.0, -1.0)
    # the plane lies `offset` along the normal from the point
    @test cut_region(v).c ≈ dot(V(v.cut_normal[]...), V(0.5, 0.0, -1.0)) + 0.25 rtol = 1e-12
    set_cut_plane!(v; angles = (90, 0))
    @test all(isapprox.(v.cut_normal[], (1, 0, 0); atol = 1e-15)) && cut_offset(v) == 0.25
    @test_throws ArgumentError set_cut_plane!(v; normal = :x, angles = (0, 0))

    # a line of sight through `o` toward (θ, φ): the plane contains the line, its normal
    # is perpendicular to it, horizontal at rotation 0, and turns about it by the rotation
    o, (θ, φ) = V(0.2, -0.4, 0.7), (70.0, 40.0)
    d = port_vector(θ, φ)
    set_cut_line!(v; angles = (θ, φ), origin = o, rotation = 0)
    n0 = V(v.cut_normal[]...)
    r = cut_region(v)
    @test all(s -> abs(dot(r.n, o + s * d) - r.c) < 1e-12, (-3.0, 0.0, 5.0))
    @test abs(dot(n0, d)) < 1e-15 && abs(n0[3]) < 1e-15 && cut_offset(v) == 0
    @test all(isapprox.(v.cut_line[], Tuple(d); atol = 1e-15))
    set_cut_line!(v; rotation = 40)
    n40 = V(v.cut_normal[]...)
    @test dot(n0, n40) ≈ cosd(40) rtol = 1e-12
    @test abs(dot(n40, d)) < 1e-15 && isapprox(dot(cross(n0, n40), d), sind(40); rtol = 1e-12)
    @test abs(dot(n40, o) - cut_region(v).c) < 1e-12
    # For a line along +y at rotation 0 the normal is ẑ × ŷ = −x̂, and along z it falls
    # back to ŷ × ẑ.
    set_cut_line!(v; direction = (0, 2, 0), origin = (0, 0, 0), rotation = 0)
    @test all(isapprox.(v.cut_normal[], (-1, 0, 0); atol = 1e-15))
    set_cut_line!(v; direction = (0, 0, 1))
    @test all(isapprox.(v.cut_normal[], (1, 0, 0); atol = 1e-15))
    # flipping keeps the plane and the line, turns the normal over: rotation + 180
    set_cut_line!(v; angles = (θ, φ), origin = o, rotation = 40)
    c40 = cut_region(v).c
    flip_cut!(v)
    @test all(isapprox.(v.cut_normal[], Tuple(-n40); atol = 1e-15)) && isapprox(cut_region(v).c, -c40; rtol = 1e-12)
    @test v.cut_rotation[] == 220 && all(isapprox.(v.cut_line[], Tuple(d); atol = 1e-15))
    # a plane cut that is not on a line releases it
    set_cut_plane!(v; normal = :z)
    @test v.cut_line[] == (0.0, 0.0, 0.0)
    @test_throws ArgumentError set_cut_line!(v; rotation = 10)

    # A plane through a flat part, a disk at c tilted to (60°, 30°), takes the centroid and
    # normal of the disk. Across the plate [-1, 1]² at the height of the disk it cuts the
    # chord n₁x + n₂y = n₁c₁ + n₂c₂. With n₂/n₁ = tan 30° and the chord crossing y = ±1, the
    # chord length is 2/cos 30°.
    sess = DesignSession(resolution = :coarse)
    c = V(0.3, -0.2, 3.0)
    add_part!(sess.target, :disk; name = :tilted, x = c[1], y = c[2], z = c[3], theta = 60, phi = 30, r1 = 0.4)
    add_window!(sess.target; phi = 90, z = 0.0, width = 1.0)
    Controllers.await_target!(sess.target)
    t = sess.target.target[]
    p, n = component_plane(t, :tilted)
    @test norm(p - c) < 1e-12 && isapprox(abs(dot(n, port_vector(60, 30))), 1; rtol = 1e-12)
    v = sess.view
    cut_through_component!(v, t, :tilted)
    @test all(isapprox.(v.cut_point[], Tuple(c); atol = 1e-12)) && cut_offset(v) == 0
    plate = [(V(-1, -1, 3), V(1, -1, 3), V(1, 1, 3)), (V(-1, -1, 3), V(1, 1, 3), V(-1, 1, 3))]
    r = cut_region(v)
    secs = reduce(vcat, (cut_triangle(r, :hohlraum_wall, tr...)[3] for tr in plate))
    @test sum(s -> norm(s[2] - s[1]), secs) ≈ 2 / cosd(30) rtol = 1e-12
    @test all(s -> all(q -> abs(dot(r.n, q) - dot(r.n, c)) < 1e-12, s), secs)
    @test _tris_area(reduce(vcat, (cut_triangle(r, :hohlraum_wall, tr...)[1] for tr in plate))) < 4
    # a window on the wall at φ = 90° is flat enough: its normal is radial (±ŷ); the
    # closed wall is not, and a cut through it keeps the current normal
    pw, nw = component_plane(t, piece_name(:hohlraum, window_name(1)))
    @test nw !== nothing && abs(nw[2]) > 0.99 && abs(pw[1]) < 1e-9 && pw[2] > 0
    @test component_plane(t, :hohlraum_wall)[2] === nothing
    n_before = v.cut_normal[]
    cut_through_component!(v, t, :hohlraum_wall)
    @test v.cut_normal[] == n_before
    @test_throws ArgumentError component_plane(t, :nothing_here)
end

@testset "display modes and the beam matrix (no GL)" begin
    sess = DesignSession(resolution = :coarse)
    v, bt = sess.view, sess.beams
    # Display modes are set per component and for several at once, and opaque is the
    # default.
    @test display_mode(v, :hohlraum_wall) === :opaque
    set_display!(v, [:hohlraum_wall, :hohlraum_cap_top], :transparent)
    set_display!(v, :capsule, :wireframe)
    @test display_mode(v, :hohlraum_wall) === display_mode(v, :hohlraum_cap_top) === :transparent
    set_display!(v, :hohlraum_wall, :opaque)
    @test !haskey(v.display[], :hohlraum_wall) && display_mode(v, :capsule) === :wireframe
    @test_throws ArgumentError set_display!(v, :hohlraum_wall, :ghost)
    set_display!(v, :capsule, :invisible)
    rows = model_tree(sess)
    cap = only(r for r in rows if r.id == "part:capsule")
    @test cap.display == "invisible" && cap.on == 0 && cap.matrix == ""
    # the beam matrix: 3ω in the calculations and drawn, unconverted light neither
    beams = only(r for r in rows if r.id == "beams")
    @test beams.matrix == "001001" && beams.display == ""
    names = [b.name for b in bt.beams[][1:2]]
    set_beams_drawn!(v, names, 1, true)
    set_beams_drawn!(v, names[1:1], 3, false)
    @test beam_drawn(v, names[1], 1) && !beam_drawn(v, names[1], 3) && beam_drawn(v, names[2], 3)
    set_participation!(bt, [1], 2, true)
    @test participates(bt, 1, 2) && !participates(bt, 2, 2) && participates(bt, 1, 3)
    @test beam_matrix(sess, [1, 2]) == "021102"
    @test beam_matrix(sess, collect(eachindex(bt.beams[]))) == "021202"
    # 3ω participation is the enabled flag, nothing else
    set_participation!(bt, [2], 3, false)
    @test !bt.enabled[][2] && beam_matrix(sess, [2]) == "000101"
    set_participation!(bt, [2], 3, true)
    # a beam without an unconverted-light footprint (not a NIF quad)
    load_preset!(bt, :omega_hohlraum)
    @test !has_unconverted(bt.beams[][1]) && beam_matrix(sess, [1]) == "--1--1"
    @test isempty(unconverted_light(bt.beams[][1], 1; length = 5.0).faces)
end

@testset "unconverted footprints as deposited (no GL)" begin
    b = first(nif_quads())
    bls = unconverted_beamlets(b)
    @test length(bls) == 4
    fp = unconverted_footprint(1)
    uc = first(nif_unconverted())
    @test (fp.edge, fp.size, fp.focus) == (uc.edge, uc.footprint, uc.focal_shift)
    # the drawn square is the zero contour of the footprint distance at every axial
    # position (`talis_check`'s square with the nominal footprint)
    sq(bl, s) = unconverted_square(bl, 1, s)
    dist(x, (k, f, ô)) = Briareus.GUIAPI.omega1_square(x - f, k, ô, fp.edge + fp.size / 2, fp.size / 2,
                                                 fp.focus, fp.lens)
    @test all(abs(dist(p, bl)) < 1e-9 for bl in bls for s in (-30.0, 0.0, 12.0) for p in sq(bl, s))
    # the TaLIS modes are the contours `talis_check` measures to
    f10, fw = unconverted_footprint(1; mode = :talis_10pct), unconverted_footprint(1; mode = :talis_worst)
    @test f10.edge == 4.45 && f10.size ≈ 22.7 && fw.edge ≈ 4.175 && fw.size ≈ 23.25
    c10 = Briareus.omega1_contour()
    @test all(abs(Briareus.GUIAPI.omega1_square(p - bl[2], bl[1], bl[3], c10.c0, c10.h0, f10.focus, f10.lens)) < 1e-9
              for bl in bls for s in (-30.0, 0.0, 12.0) for p in unconverted_square(bl, 1, s; mode = :talis_10pct))
    @test_throws ArgumentError unconverted_footprint(2; mode = :talis_10pct)
    # its near edge lies `edge` from the 3ω axis at the 3ω focus
    k, f, ô = bls[1]
    @test minimum(dot(p - f, ô) for p in sq(bls[1], 0.0)) ≈ fp.edge rtol = 1e-12
    u = unconverted_light(b, 1; length = 10.0)
    @test length(u.vertices) == 32 && length(u.faces) == 32 && length(u.outlines) == 16
    # the drawn 1ω light follows the chosen contour: its near edge at the 3ω focus, and
    # the 2ω light keeps its nominal footprint
    near(u) = minimum(dot(p - f, ô) for seg in u.outlines[1:4] for p in seg)
    @test near(u) ≈ fp.edge && near(unconverted_light(b, 1; length = 10.0, mode = :talis_10pct)) ≈ c10.edge
    @test near(unconverted_light(b, 1; length = 10.0, mode = :talis_worst)) ≈ Briareus.omega1_contour(mode = :talis_worst).edge
    @test unconverted_light(b, 2; length = 10.0, mode = :talis_worst).outlines == unconverted_light(b, 2; length = 10.0).outlines
    vs = ViewState()
    @test vs.omega1_contour[] === :nominal && set_omega1_contour!(vs, :talis_10pct).omega1_contour[] === :talis_10pct
    @test_throws ArgumentError set_omega1_contour!(vs, :other)
end

@testset "1ω closest-approach geometry (no GL)" begin
    # a point 0.5 mm outside the near side of the 10% contour of one beamlet
    bl = first(nif_beams())
    k = direction(bl)
    a = omega1_closest_approach(Target(foil(:f, (0, 0, 0), (0, 0, 1); width = 1.0, n = 2)), [bl],
                                1 => Tuple(unconverted_beamlets(bl)[1][3] * 3.95))
    ô = a.beamlet[3]
    @test a.face_mask == 0x01 && a.focus_shift == 1.0
    g = closest_approach_surface(a; span = 6.0)
    # the outline is the near side of the contour at the worst focal plane, 12 mm long
    c10 = Briareus.omega1_contour()
    dist(x) = Briareus.GUIAPI.omega1_square(x - (a.beamlet[2] + a.focus_shift * k), k, ô, c10.c0, c10.h0,
                                           c10.focal_shift, a.contour.lens_distance)
    @test length(g.outline) == 4 && all(abs(dist(p)) < 1e-9 for s in g.outline for p in s)
    ss = [dot(p - a.beamlet[2] - a.focus_shift * k, k) for s in g.outline for p in s]
    @test minimum(ss) ≈ a.s_star - 6 && maximum(ss) ≈ a.s_star + 6
    # the clearance plane lies outside the contour by the errors, the tolerance and the
    # required clearance, and the hardware point is `margin` beyond it
    grow = a.budget.pointing + a.budget.target_offset + a.budget.rotation_tol + a.budget.required
    @test length(g.vertices) == 4 && g.faces == [(1, 2, 3), (1, 3, 4)]
    @test isapprox(dist((g.vertices[1] + g.vertices[2]) / 2), grow; atol = 1e-9) &&
          isapprox(dist((g.vertices[3] + g.vertices[4]) / 2), grow; atol = 1e-9)
    @test a.distance - grow ≈ a.margin atol = 1e-12
    @test g.caliper == (a.x_hw, a.x_beam) && length(g.budget) == 8
    @test startswith(g.budget[1], "nominal") && startswith(g.budget[end], "margin") && endswith(g.budget[end], "mm")
    # a corner shows two sides
    ac = omega1_closest_approach(Target(foil(:f, (0, 0, 0), (0, 0, 1); width = 1.0, n = 2)), [bl],
                                 1 => Tuple(3.9 * ô - 12.0 * cross(k, ô)))
    @test count_ones(ac.face_mask) == 2 && length(closest_approach_surface(ac).faces) == 4
    # the daylight view looks down the beamlet at the middle of the caliper; the grazing
    # view looks along the near side with the beam up the view
    dv, gv = daylight_view(a), grazing_view(a)
    @test dv.direction == k && dv.center ≈ (a.x_hw + a.x_beam) / 2 && abs(dot(dv.up, k)) < 1e-12
    @test dv.width ≈ max(4 * abs(a.distance), 1.0) && abs(dot(dv.up, a.x_hw - a.x_beam)) < 1e-9
    @test gv.up == k && abs(dot(gv.direction, k)) < 1e-12 && abs(dot(gv.direction, a.x_hw - a.x_beam)) < 1e-9
    # the cone of a reflection: its axis, its half-angle and its rim
    x = Omega1Reflection("b", :c, 1, k, Briareus.Vec3(1, 2, 3), -k, -k, bl.port, 1, 4.0, 3.0, 0.0, -20.0, 16.0, 16.0, 0.03)
    rc = reflection_cone(x; length = 10.0, n = 16)
    @test rc.half_angle == 20 + 1 + 3 && rc.incident == (x.point - 10 * k, x.point) && rc.reflected[2] ≈ x.point - 10 * k
    @test length(rc.vertices) == 17 && length(rc.faces) == 16 && rc.port == (x.point, x.point + 10 * bl.port)
    @test all(v -> isapprox(acosd(dot(normalize(v - x.point), -k)), 24; atol = 1e-9) && norm(v - x.point) ≈ 10,
              rc.vertices[2:end])
    # the ring that marks the port: a closed circle about the port direction, 3° in radius
    @test length(rc.ring) == 16 && all(rc.ring[i][2] ≈ rc.ring[mod1(i + 1, 16)][1] for i in 1:16)
    @test all(s -> isapprox(acosd(dot(normalize(s[1] - x.point), bl.port)), PORT_MARK_ANGLE; atol = 1e-6) &&
                   norm(s[1] - x.point) ≈ 10, rc.ring)
end

@testset "participation sets the drive and visibility does not (no GL)" begin
    sess = DesignSession(resolution = :coarse)
    bt, v = sess.beams, sess.view
    b = bt.beams[][1]
    # a two-sided disk centered in the first beamlet's 1ω footprint, facing its light
    k, f, ô = unconverted_beamlets(b)[1]
    c = sum(unconverted_square((k, f, ô), 1, 0.0)) / 4
    θ, φ = acosd(-k[3]), atand(-k[2], -k[1])
    add_part!(sess.target, :disk; name = :shield, two_sided = true, x = c[1], y = c[2], z = c[3],
              theta = θ, phi = φ, r1 = 3.0, conversion = 0.5)
    build_scenario!(sess)
    res0 = sess.drive.result[]
    t = res0.evaluation.scenario.target
    shield(res) = (cp = component_power(res.evaluation.deposition, t); get(cp, :shield, 0.0) + get(cp, :shield_outside, 0.0))
    @test res0.unconverted === nothing && shield(res0) < 1e-9

    # A hidden beam still deposits exactly as when it is drawn.
    set_beams_drawn!(v, [b.name], 3, false)
    @test sess.drive.result[] === res0 && res0.evaluation.deposition.per_beam[1] > 0
    set_beams_drawn!(v, [b.name], 3, true)

    # 1ω participation without a ratio changes nothing; with one, its power reaches the disk
    set_participation!(bt, [1], 1, true)
    @test unconverted_ratios(sess) === nothing && sess.drive.result[].unconverted === nothing
    set_unconverted_ratio!(sess.drive, 0.05, 0.0)
    res1 = sess.drive.result[]
    @test res1.unconverted !== nothing
    @test res1.unconverted.omega1.delivered ≈ 0.05 * b.power rtol = 1e-12
    @test shield(res1) > 0.01 * 0.05 * b.power
    @test res1.evaluation.deposition.delivered ≈ res0.evaluation.deposition.delivered + 0.05 * b.power rtol = 1e-12
    @test occursin("1ω", unconverted_summary(res1))
    # the x-ray source of the unconverted light scales with its color's `conversion`
    # (two colors sets that differ only in it): the absorbed power is the same
    dr = sess.drive
    ratio = unconverted_ratios(sess)
    half = [Briareus.UnconvertedColor(c.name; footprint = c.footprint, edge = c.edge,
                                      focal_shift = c.focal_shift, absorption = c.absorption,
                                      conversion = c.conversion / 2) for c in nif_unconverted()]
    run(colors) = Controllers._evaluate_drive(dr.scenario[], display_beams(bt), dr.lookup[],
                                              dr.spot_radius[]; ratio, colors)
    full, halved = run(nif_unconverted()), run(half)
    @test full.evaluation.deposition.element == halved.evaluation.deposition.element
    P1 = unconverted_deposition(t, display_beams(bt); ratio).omega1.element
    expected = nif_unconverted()[1].conversion / 2 .* Briareus.conversions(t) .* P1 ./ t.elem_area
    @test isapprox(full.solution.source - halved.solution.source, expected; rtol = 1e-12, atol = 1e-15)
    @test sum(full.solution.source .- halved.solution.source) > 0          # not a vacuous check
    @test full.unconverted.omega1.conversion == 2 * halved.unconverted.omega1.conversion
    # drawing the 1ω light (or not) leaves that drive untouched
    set_beams_drawn!(v, [b.name], 1, true)
    set_beams_drawn!(v, [b.name], 1, false)
    @test sess.drive.result[] === res1
    # Out of the calculations again, the drive is identical to the drive of the 3ω light
    # only.
    set_participation!(bt, [1], 1, false)
    @test sess.drive.result[].evaluation.deposition.element == res0.evaluation.deposition.element

    # the ratio and participation are design state: saved and undone
    record!(sess)
    set_participation!(bt, [1, 2], 1, true)
    record!(sess)
    @test undo!(sess) && !participates(bt, 2, 1) && sess.drive.unconverted_ratio[] == (0.05, 0.0)
end

@testset "display and cut state in session files (no GL)" begin
    s = DesignSession(resolution = :coarse)
    v = s.view
    set_display!(v, :capsule, :wireframe)
    set_display!(v, :hohlraum_cap_top, :invisible)
    set_beams_drawn!(v, [s.beams.beams[][3].name], 3, false)
    set_beams_drawn!(v, [s.beams.beams[][4].name], 2, true)
    set_participation!(s.beams, [5], 1, true)
    set_unconverted_ratio!(s.drive, 0.03, 0.01)
    set_cut!(v; axis = :x, azimuth = 30, width = 120, range = (-1.0, Inf), components = [:hohlraum_wall])
    set_cut_plane!(v; normal = :y, offset = 0.4)
    path = joinpath(mktempdir(), "display.briar.toml")
    s_los = DesignSession(resolution = :coarse)
    set_cut_line!(s_los.view; angles = (80, 200), origin = (0.1, 0.2, 0.3), rotation = 35)
    set_cut_plane!(s_los.view; offset = -0.3)
    save_session(path, s_los)
    s_los2 = DesignSession()
    load_session!(s_los2, path)
    w1, w2 = s_los.view, s_los2.view
    @test session_state(s_los2) == session_state(s_los)
    @test w2.cut_line[] == w1.cut_line[] && w2.cut_rotation[] == 35 && w2.cut_point[] == (0.1, 0.2, 0.3)
    @test cut_offset(w2) == -0.3 && w2.cut_normal[] == w1.cut_normal[] && w2.cut_kind[] === :plane
    save_session(path, s)
    s2 = DesignSession()
    load_session!(s2, path)
    v2 = s2.view
    @test session_state(s2) == session_state(s)
    @test v2.display[] == v.display[] && v2.hidden_beams[] == v.hidden_beams[] && v2.shown_2w[] == v.shown_2w[]
    @test v2.cut_kind[] === :plane && v2.cut_axis[] === :x && v2.cut_range[] == (-1.0, Inf)
    @test all(isapprox.(v2.cut_normal[], (0.0, 1.0, 0.0))) && cut_offset(v2) ≈ 0.4
    @test participates(s2.beams, 5, 1) && s2.drive.unconverted_ratio[] == (0.03, 0.01)
end
