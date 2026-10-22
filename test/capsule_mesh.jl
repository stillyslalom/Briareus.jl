# The geodesic (compact) sphere mesh and `hohlraum(capsule_mesh = :compact)`.

using Briareus, LinearAlgebra, Test

@testset "compact capsule mesh" begin
    R = 0.9
    for n in (1, 3, 8)
        m = icosphere_mesh(R; n)
        @test npatches(m) == 20n^2
        tris = [Briareus.triangle_vertices(m, k) for k in 1:Briareus.ntriangles(m)]
        @test Briareus.ntriangles(m) == npatches(m)                 # single-triangle patches
        @test all(v -> norm(v) ≈ R, m.vertices)
        # outward, and closed: every directed edge meets its reverse exactly once
        @test all(dot(cross(b - a, c - a), a + b + c) > 0 for (a, b, c) in tris)
        edges = Dict{Tuple{Vec3,Vec3},Int}()
        for (a, b, c) in tris, e in ((a, b), (b, c), (c, a))
            edges[e] = get(edges, e, 0) + 1
        end
        @test all(==(1), values(edges)) && all(haskey(edges, (b, a)) for (a, b) in keys(edges))
        @test length(unique(m.vertices)) == 10n^2 + 2
        # the flat facets fall short of the sphere by about (h/R)²/4; on the surface, exact
        flat = sum(norm(cross(b - a, c - a)) / 2 for (a, b, c) in tris)
        hR = atan(2.0) / n
        @test 0 < 1 - flat / (4π * R^2) < 0.3hR^2
        t = Target(Component(:s, m, gold(); surface = SphereSurface(R)))
        @test sum(t.elem_area) ≈ 4π * R^2 rtol = 1e-12
        # nearly equilateral: edges within -20% / +20% of the nominal size
        ls = [norm(p - q) for (a, b, c) in tris for (p, q) in ((a, b), (b, c), (c, a))]
        n > 1 && @test 0.8 < minimum(ls) / (hR * R) && maximum(ls) / (hR * R) < 1.2
    end
    @test all(dot(cross(b - a, c - a), a + b + c) < 0
              for (a, b, c) in (Briareus.triangle_vertices(icosphere_mesh(1.0; n = 2, outward = false), k)
                                for k in 1:80))
    # element count ∝ 1/h²
    @test icosphere_frequency(R, 0.1) == 10 && icosphere_frequency(R, 0.05) == 20
    @test icosphere_frequency(R, 10.0) == 1
    @test_throws ArgumentError icosphere_frequency(R, 0.0)

    # convex faceting: from inside, every element sees all of the sphere (F_ij = A_j/4πR²)
    t = Target(Component(:s, icosphere_mesh(2.0; n = 3, outward = false), gold(); surface = SphereSurface(2.0)))
    F = view_factors(t)
    @test maximum(abs.(F .- (t.elem_area ./ sum(t.elem_area))')) < 1e-12

    # a uniform field has no l ≥ 1 moments on the compact capsule (exact solid angles)
    tc = Target(Component(:c, icosphere_mesh(R; n = 6), ablator(); surface = SphereSurface(R)))
    es = collect(1:nelements(tc))
    c = ylm_weights(tc, es; lmax = 4) * ones(length(es))
    @test c[1] ≈ sqrt(4π) rtol = 1e-12
    @test maximum(abs, c[2:end]) < 1e-12
    a = legendre_moments(ones(length(es)), tc, es; lmax = 4)
    @test maximum(abs, a[2:end]) < 1e-12 * a[1]
    # an axisymmetric P2 + P4 field, piecewise constant on the elements
    u = [normalize(x) for x in tc.elem_centroid]
    v = [1 + 0.1legendre(2, x[3]) - 0.03legendre(4, x[3]) for x in u]
    a = legendre_moments(v, tc, es; lmax = 4)
    @test a[3] / a[1] ≈ 0.1 rtol = 0.02
    @test a[5] / a[1] ≈ -0.03 rtol = 0.05

    # in a hohlraum: a uniform wall source drives the capsule as on the lat-long mesh
    kw = (; h = 0.6, capsule_h = 0.15, wall = Material(:Au; albedo = 0.85, conversion = 1.0))
    th = hohlraum(; kw..., capsule_mesh = :compact)
    tl = hohlraum(; kw..., capsule_mesh = :latlong)
    nc = icosphere_frequency(0.9, 0.15)
    @test length(elements(th, :capsule)) == 20nc^2
    @test th.components[component_index(th, :capsule)].surface isa SphereSurface
    @test sum(th.elem_area[elements(th, :capsule)]) ≈ 4π * 0.9^2 rtol = 1e-12
    @test length(elements(th, :wall)) == length(elements(tl, :wall))
    drive(t) = (P = zeros(nelements(t)); w = elements(t, :wall); P[w] .= t.elem_area[w];
                sol = solve(RadiationModel(t), P); cap = elements(t, :capsule);
                legendre_moments(sol.H, t, cap; lmax = 4))
    ah, al = drive(th), drive(tl)
    @test ah[1] ≈ al[1] rtol = 2e-3
    @test ah[3] / ah[1] ≈ al[3] / al[1] atol = 2e-3
    @test abs(ah[2] / ah[1]) < 1e-6           # the icosphere is centrally symmetric
    @test hohlraum(; capsule_mesh = :compact, capsule_ntheta = 12) isa Target   # sized from the counts
    @test_throws ArgumentError hohlraum(; capsule_mesh = :geodesic)
end
