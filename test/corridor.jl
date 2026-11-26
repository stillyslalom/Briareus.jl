@testset "corridor map" begin
    R = 100.0
    zero_tol = (; tolerance = 0.0, pointing = 0.0, target_offset = 0.0, hardware_tolerance = 0.0, focus = 0.0)

    # A cylindrical sight line along +z, 20 mm across: on the sphere the distance to its
    # side is R sinθ - 10 in front of the center plane, and behind it the distance to the
    # starting disk.
    cyl = SightLine("cyl", 0, 0; field_of_view = 20.0)
    m = corridor_map(Beam[]; sight_lines = [cyl], radius = R, ntheta = 60, nphi = 24)
    @test m.surface === :sphere && size(m.margin) == (60, 24) && m.u[1] == 1.5 && m.v[1] == 7.5
    exact(θ) = θ < 90 ? R * sind(θ) - 10 : hypot(R * cosd(θ), max(R * sind(θ) - 10, 0))
    @test all(abs(m.sight[i, j] - exact(m.u[i])) < 1e-10 for i in 1:60, j in 1:24)
    @test m.margin == m.sight && all(isinf, m.omega3) && all(isinf, m.omega1)
    @test all(==(3), m.nearest .+ 2) && m.sources == [(:sight, "cyl")]
    @test all(norm(m.points[i, j]) ≈ R && m.points[i, j][3] ≈ R * cosd(m.u[i]) for i in 1:60, j in 1:24)
    cls = corridor_classes(m)
    @test all(cls[i, j] == (m.u[i] < asind(0.1) ? 0x04 : 0x00) for i in 1:60, j in 1:24)
    # with `through`, the cylinder also covers the far pole
    mt = corridor_map(Beam[]; sight_lines = [SightLine("cyl", 0, 0; field_of_view = 20.0, through = true)],
                      radius = R, ntheta = 60, nphi = 24)
    @test all(abs(mt.sight[i, j] - (R * sind(mt.u[i]) - 10)) < 1e-10 for i in 1:60, j in 1:24)

    # A cone with its apex at the center, of half-angle α about the direction (θ0, φ0):
    # the angular margin is the angle from the axis less α, inside and outside.
    α, θ0, φ0 = 12.0, 50.0, 200.0
    cone = SightLine("cone", θ0, φ0; field_of_view = 0.0, distance = 5000.0, aperture = 2 * 5000 * tand(α))
    mc = corridor_map(Beam[]; sight_lines = [cone], radius = R)
    ax = port_vector(θ0, φ0)
    off = [acosd(clamp(dot(x, ax) / R, -1, 1)) for x in mc.points]
    am = angular_margin(mc)
    near = off .< 85
    @test maximum(abs.(am[near] .- (off[near] .- α))) < 1e-9
    @test count(<(0), am) > 0 && all(am[off .> α + 1e-9] .> 0)
    # the open share of the solid angle: all but the cap of half-angle α, sampled on 2° cells
    @test abs(corridor_open_fraction(mc) - (1 - (1 - cosd(α)) / 2)) < 2e-3
    @test corridor_open_fraction(mc; margin = R * sind(10.0)) < corridor_open_fraction(mc)
    @test_throws ArgumentError angular_margin(corridor_map(Beam[]; plane = (zeros(3), (0, 0, 1)), extent = 5.0, n = 5))

    # One 3ω beam with a round super-Gaussian spot along -z. The 1% contour of
    # exp(-(r/a)^n) lies at a ln(100)^(1/n), and the envelope grows by |s| / (2F).
    a, n, F = 0.6, 4.0, 8.0
    b = Beam("b"; theta = 0, phi = 0, power = 1.0, spot = SuperGaussianSpot(a; order = n), fnumber = F,
             aim = (0.3, -0.2, 0.1))
    mb = corridor_map([b]; radius = 30.0, ntheta = 45, nphi = 36, clearance_3w = 0.0, zero_tol...)
    ρc = a * log(100)^(1 / n)
    env(x) = hypot(x[1] - 0.3, x[2] + 0.2) - ρc - abs(x[3] - 0.1) / (2F)
    @test maximum(abs(mb.omega3[i] - env(mb.points[i])) for i in eachindex(mb.points)) < 1e-9
    @test all(isinf, mb.omega1) && mb.sources == [(:omega3, "b"), (:omega1, "b")]
    # the errors and the required clearance come off the distance: 0.075 mm of position,
    # the focus offset over 2F, the rotation about TCC and 0.5 mm
    x = Vec3(4.0, 3.0, 12.0)
    q = corridor_margin([b], x)
    @test q.omega3 ≈ env(x) - 0.075 - 1 / (2F) - tand(0.5) * norm(x) - 0.5 atol = 1e-10
    @test q.margin == q.omega3 && q.nearest == (:omega3, "b") && q.sight == Inf

    # The 1ω footprint of one NIF beam. On the 3ω axis of the beamlet, s past the focus,
    # the near side of the square is c0 (1 + s/L) - h0 |1 - s/f| away.
    bl = first(nif_beams())
    k = direction(bl)
    o = omega1_contour()
    f0 = Vec3(bl.aim) + bl.defocus * k
    on_axis(s) = corridor_margin([bl], f0 + s * k; omega1_margin = 0.0, zero_tol...).omega1
    ss = (-40.0, -3.0, 0.0, 25.0, 200.0)
    @test all(abs(on_axis(s) - (o.c0 * (1 + s / 7700) - o.h0 * abs(1 - s / o.focal_shift))) < 1e-9 for s in ss)
    # off the axis the margin is the worst-case distance of the closest approach less 0.2 mm
    t1 = Target(Component(:p, rectangle_mesh(0.1, 0.1; nu = 2, nv = 2), gold(); two_sided = true))
    xs = [f0 + s * k + 14.0 * normalize(cross(k, Vec3(0.3, 0.5, 0.8))) for s in (-20.0, 0.0, 30.0)]
    @test all(abs(corridor_margin([bl], x).omega1 - omega1_closest_approach(t1, [bl], bl.name => x).margin) < 1e-10
              for x in xs)
    qd = corridor_margin([bl], xs[1]; omega1_margin = :diagnostic)
    @test qd.omega1 ≈ corridor_margin([bl], xs[1]).omega1 + 0.2 - (10 + 0.015 * norm(xs[1]))
    @test_throws ArgumentError corridor_margin([bl], xs[1]; omega1_margin = :other)
    # both colors and a sight line in one map: every class is present and `nearest` names
    # the keep-out with the smallest margin
    mx = corridor_map([bl]; sight_lines = [cyl], radius = 60.0)
    cx = corridor_classes(mx)
    @test any(!iszero, cx .& 0x01) && any(!iszero, cx .& 0x02) && any(!iszero, cx .& 0x04) && any(iszero, cx)
    pick(i) = (mx.omega3[i], mx.omega1[i], mx.sight[i])[mx.nearest[i]]
    @test all(pick(i) == mx.margin[i] for i in eachindex(mx.margin))
    @test corridor_map([bl]; radius = 60.0, omega1 = false).omega1 == fill(Inf, 90, 180)
    fl = corridor_fills(mx)
    @test all(fl[i] == (mx.omega3[i] <= 0 ? 1 : mx.omega1[i] <= 0 ? 2 : mx.sight[i] <= 0 ? 3 : mx.margin[i] < 1 ? 4 : 5)
              for i in eachindex(fl))

    # a plane across the through-going cylinder: the margin is the radius less 10 mm
    mp = corridor_map(Beam[]; sight_lines = [SightLine("cyl", 0, 0; field_of_view = 20.0, through = true)],
                      plane = ((0, 0, -7), (0, 0, 1)), extent = 25.0, n = 51)
    @test mp.surface === :plane && isnan(mp.radius) && mp.axes == (Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1))
    @test all(abs(mp.margin[i, j] - (hypot(mp.u[i], mp.v[j]) - 10)) < 1e-10 for i in 1:51, j in 1:51)
    @test abs(corridor_open_fraction(mp) - (1 - π * 100 / 2500)) < 0.02
    tilted = corridor_map(Beam[]; plane = ((1, 2, 3), (1, 0, 1)), extent = 2.0, n = 3)
    @test all(abs(dot(x - Vec3(1, 2, 3), normalize(Vec3(1, 0, 1)))) < 1e-12 for x in tilted.points)
    @test all(isinf, tilted.margin) && corridor_open_fraction(tilted) == 1 && all(iszero, tilted.nearest)
    @test corridor_margin(Beam[], (0, 0, 0)).nearest === nothing
    @test_throws ArgumentError corridor_map(Beam[])
    @test_throws ArgumentError corridor_map(Beam[]; radius = 10.0, plane = (zeros(3), (0, 0, 1)), extent = 1.0)
    @test_throws ArgumentError corridor_map(Beam[]; plane = (zeros(3), (0, 0, 1)))

    # sight lines of the built-in definitions
    d1 = sight_line("nif_dante1")
    @test d1.axis ≈ port_vector(143.25, 274.0) && d1.radius == 30.0 && d1.end_radius == 30.0 && isinf(d1.distance)
    sx = sight_line("nif_sxi_upper")
    @test sx.radius ≈ 24e-3 * 2048 / 2 / 2 && sx.distance == 1626.7 && sx.end_radius ≈ 0.0125
    @test_throws ArgumentError sight_line("nif_inflight_radiography")
    rg = sight_line("nif_inflight_radiography"; phi = 78.0, field_of_view = 2.0, through = true)
    @test rg.axis ≈ port_vector(90.0, 78.0) && rg.through && rg.radius == 1.0
    @test_throws ArgumentError sight_line("no_such_diagnostic")
    @test_throws ArgumentError SightLine("x", 0, 0; field_of_view = 2.0, aperture = 1.0)
    # the pinhole cone: inside on the axis, and beside the pinhole the perpendicular
    # distance to the slanted side
    @test sight_distance(sx, 800.0 * sx.axis) < 0 && sight_distance(sx, -1.0 * sx.axis) ≈ 1.0
    slant = cos(atan((sx.radius - sx.end_radius) / sx.distance))
    @test sight_distance(sx, 1626.7 * sx.axis + 5.0 * normalize(cross(sx.axis, Vec3(1, 0, 0)))) ≈
          (5.0 - 0.0125) * slant

    # The parts that enter a line of sight along z of radius 1 mm, against closed forms:
    # a plate beside the cylinder 0.5 mm from its side, a plate across it (deepest on the
    # axis, 1 mm inside) and a plate behind the plane of the center (3 mm from it).
    up = SightLine("up", 0, 0; field_of_view = 2.0)
    side = foil(:side, (2.5, 0.0, 5.0), (0, 0, 1); width = 2.0, n = 4)
    across = foil(:across, (0.0, 0.0, 5.0), (0, 0, 1); width = 4.0, n = 4)
    behind = foil(:behind, (0.0, 0.0, -3.0), (0, 0, 1); width = 4.0, n = 4)
    tp = Target([side, across, behind])
    si = sight_intrusions(tp, up)
    h = 0.1                                     # the default spacing, a tenth of the radius
    @test [r.part for r in si] == [:across, :side, :behind] && [r.intrudes for r in si] == [true, false, false]
    @test all(r -> 0 <= r.margin - Dict(:across => -1.0, :side => 0.5, :behind => 3.0)[r.part] <= 2h / 3, si)
    @test [r.intrudes for r in sight_intrusions(tp, up; tolerance = 0.6)] == [true, true, false]
    fine = sight_intrusions(tp, up; spacing = 0.01)
    @test all(r -> 0 <= r.margin - Dict(:across => -1.0, :side => 0.5, :behind => 3.0)[r.part] <= 0.01, fine)
    @test_throws ArgumentError sight_intrusions(tp, up; spacing = 0.0)

    # drawings: plain text, the same for the same map, one legend entry per fill
    svg = corridor_map_svg(mx)
    @test startswith(svg, "<svg") && endswith(svg, "</svg>\n") && svg == corridor_map_svg(mx)
    @test all(occursin(c, svg) for c in Briareus._CORRIDOR_FILLS) && occursin("% open", svg)
    @test count("<path", svg) > 100 && count("<circle", svg) == 6
    @test !occursin("NaN", svg) && !occursin("Inf", svg)
    ps = corridor_map_svg(mp; size = 200)
    @test occursin("width=\"200\"", ps) && count("<rect", ps) > 51 && !occursin("<path", ps)
    @test_throws ArgumentError corridor_map_svg(mx; size = 50)
    path = joinpath(mktempdir(), "map.svg")
    @test write_corridor_map_svg(path, mx) == path && read(path, String) == svg
    @test occursin("open)", sprint(show, mx))
end
