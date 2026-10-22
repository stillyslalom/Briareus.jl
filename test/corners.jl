# Rounded wall–cap corners: `hohlraum(corner_radius = rc)`.

using Briareus, LinearAlgebra, Test

@testset "rounded hohlraum corners" begin
    kw = (; length = 9.4, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9)
    # rc = 0 is the square-cornered hohlraum, element for element
    t0, tz = hohlraum(; kw..., h = 0.5), hohlraum(; kw..., h = 0.5, corner_radius = 0)
    @test t0.elem_centroid == tz.elem_centroid && t0.elem_area == tz.elem_area

    rc = 1.0
    t = hohlraum(; kw..., h = 0.5, corner_radius = rc)
    names = [c.name for c in t.components]
    k = count(n -> startswith(String(n), "corner_top_"), names)
    @test k == 4                     # about one wall cell (0.45 mm around) of arc per band
    @test count(n -> startswith(String(n), "corner_bottom_"), names) == k
    corner = [e for n in names if startswith(String(n), "corner_") for e in elements(t, n)]
    # each band is its exact frustum: the fillet's area is the chord polygon's, just under
    # the quarter torus 2π rc (π/2 (R − rc) + rc) per corner
    A = sum(t.elem_area[corner])
    torus = 2 * 2π * rc * (π / 2 * (2.875 - rc) + rc)
    @test 0 < 1 - A / torus < (π / 4k)^2
    # wall and caps shrink by the fillets
    @test sum(t.elem_area[elements(t, :wall)]) ≈ 2π * 2.875 * (9.4 - 2rc) rtol = 1e-12
    @test sum(t.elem_area[elements(t, :cap_top)]) ≈ π * ((2.875 - rc)^2 - 1.55^2) rtol = 1e-12
    # bands face the cavity: their normals point toward the fillet's center of curvature
    for e in corner
        c, n = t.elem_centroid[e], t.elem_normal[e]
        ρ̂ = normalize(Vec3(c[1], c[2], 0))
        center = (2.875 - rc) * ρ̂ + Vec3(0, 0, sign(c[3]) * (4.7 - rc))
        @test dot(n, center - c) > 0
    end
    # the cavity is closed: raw row sums are near one, as for the square hohlraum
    F = view_factors(t)
    @test maximum(abs, closure_error(F)) < 1e-2
    @test maximum(abs, closure_error(F)[corner]) < 1e-2
    # windows must lie on the straight wall
    @test hohlraum(; kw..., h = 0.5, corner_radius = rc,
                   windows = [WallWindow(phi = 0, z = 0, width = 1.0, height = 1.0)]) isa Target
    @test_throws ArgumentError hohlraum(; kw..., h = 0.5, corner_radius = rc,
                                        windows = [WallWindow(phi = 0, z = 3.5, width = 1.0, height = 1.0)])
    @test_throws ArgumentError hohlraum(; kw..., corner_radius = 1.4)    # wider than the cap
    @test_throws ArgumentError hohlraum(; kw..., corner_radius = -0.1)

    # the D2-symmetric NIF quads put no forbidden modes on the capsule
    q = [remake(b; aim = (0.0, 0.0, sign(b.port[3]) * (b.group in (:inner23, :inner30) ? 3.3 : 4.4)))
         for b in nif_quads()]
    tq = hohlraum(; kw..., nz = 16, nphi = 24, corner_radius = rc, nr_corner = 3)
    sc = Scenario(tq, q; observe = :capsule)
    cf = flux_ylm(evaluate(sc, q); lmax = 4)
    forbidden = [(l, m) for l in 0:4 for m in -l:l if (m >= 0 ? isodd(l + m) : iseven(l + m))]
    @test maximum(abs(cf[mode_index(l, m)]) for (l, m) in forbidden) < 1e-10 * cf[1]
    # NIF cones cannot reach the fillet through the LEH (a 50° beam reaches z ≤ 3.6 here: the
    # cap shadows it); a 70° beam can, and operator and ray trace agree on the curved bands
    b = Beam("c"; theta = 70, phi = 0, power = 1.0, spot = SuperGaussianSpot(0.15; order = 4.0),
             fnumber = 20.0, aim = (-2.70, 0.0, 4.15))
    cs = [n for n in (c.name for c in tq.components) if startswith(String(n), "corner_top_")]
    d = deposit(DepositionOperator(tq, [b]), [b])
    oncorner = sum(sum(d.element[elements(tq, n)]) for n in cs)
    @test oncorner > 0.9
    r = raytrace_deposition(tq, [b]; nspot = 121)
    @test sum(sum(r.element[elements(tq, n)]) for n in cs) ≈ oncorner rtol = 0.02
end
