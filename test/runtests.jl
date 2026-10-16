using Briareus
using Briareus: closest_hit, any_hit, _point_polygon_vf
using LinearAlgebra
using ForwardDiff
using Random
using Test

const R_SPHERE = 2.0

@testset "Briareus" begin

@testset "meshes" begin
    s = sphere_mesh(R_SPHERE; ntheta = 16, nphi = 24)
    t = Target(Component(:s, s, gold()))
    @test nelements(t) == 16 * 24
    @test isapprox(sum(t.elem_area), 4π * R_SPHERE^2; rtol = 0.03)
    # equal-area bands
    @test maximum(t.elem_area) / minimum(t.elem_area) < 1.1
    # outward normals
    @test all(dot(t.elem_normal[e], t.elem_centroid[e]) > 0 for e in 1:nelements(t))
    c = Target(Component(:c, cylinder_mesh(1.0, -1, 1; nz = 4, nphi = 12), gold()))
    @test all(dot(c.elem_normal[e], Vec3(c.elem_centroid[e][1], c.elem_centroid[e][2], 0)) < 0
              for e in 1:nelements(c))
    a = annulus_mesh(0.5, 1.0, 0.0; nr = 3, nphi = 64, up = false)
    ta = Target(Component(:a, a, gold()))
    @test isapprox(sum(ta.elem_area), π * (1 - 0.25); rtol = 0.01)
    @test all(n -> n[3] ≈ -1, ta.elem_normal)
    f = foil(:f, (1.0, 2.0, 3.0), (0, 1, 0); width = 2.0, n = 4)
    tf = Target(f)
    @test nelements(tf) == 2 * 16                # two-sided
    @test count(n -> n ≈ Vec3(0, 1, 0), tf.elem_normal) == 16
    @test count(n -> n ≈ Vec3(0, -1, 0), tf.elem_normal) == 16
    @test all(c -> c[2] ≈ 2.0, tf.elem_centroid)
    R = rotation_to(Vec3(1, 1, 1))
    @test R * Vec3(0, 0, 1) ≈ normalize(Vec3(1, 1, 1))
    @test det(R) ≈ 1
end

@testset "ray casting" begin
    t = hohlraum(nz = 8, nphi = 12, capsule_ntheta = 6, capsule_nphi = 8)
    rng = MersenneTwister(1)
    for _ in 1:300
        o = Vec3(randn(rng, 3)...) .* 2
        d = normalize(Vec3(randn(rng, 3)...))
        tb, idb = closest_hit(t.bvh, o, d)
        # brute force over opaque triangles
        best, bi = Inf, 0
        for k in eachindex(t.tri_area)
            t.tri_opaque[k] || continue
            tt = Briareus._intersect(o, d, t.tri_v0[k], t.tri_e1[k], t.tri_e2[k], best)
            tt < best && (best = tt; bi = k)
        end
        @test tb ≈ best
        @test any_hit(t.bvh, o, d) == (bi != 0)
    end
end

@testset "view factors" begin
    # Interior of a sphere: F_ij = A_j / (4πR²) for every pair, rows close to 1.
    s = Target(Component(:s, sphere_mesh(R_SPHERE; ntheta = 10, nphi = 16, outward = false), gold()))
    F = view_factors(s)
    A = s.elem_area
    @test maximum(abs.(closure_error(F))) < 0.01
    # two mid-latitude facets: F_ij ≈ A_j / ΣA (polar slivers deviate more, being flat)
i, j = 2 * 16 + 1, 7 * 16 + 9
@test isapprox(F[i, j], A[j] / sum(A); rtol = 0.03)
    # reciprocity
    @test maximum(abs.(A .* F .- (A .* F)')) < 1e-12

    # Coaxial parallel disks, radii r1, r2 at separation h (analytic).
    r1, r2, h = 1.0, 0.7, 0.8
    d1 = Component(:d1, disk_mesh(r1, 0.0; nr = 8, nphi = 32, up = true), gold())
    d2 = Component(:d2, disk_mesh(r2, h; nr = 8, nphi = 32, up = false), gold())
    t = Target(d1, d2)
    F = view_factors(t)
    e1, e2 = elements(t, :d1), elements(t, :d2)
    F12 = sum(t.elem_area[i] * F[i, j] for i in e1, j in e2) / sum(t.elem_area[e1])
    R1, R2 = r1 / h, r2 / h
    X = 1 + (1 + R2^2) / R1^2
    F12_exact = (X - sqrt(X^2 - 4 * (R2 / R1)^2)) / 2
    @test isapprox(F12, F12_exact; rtol = 0.01)

    # Occlusion: an opaque plate between the disks blocks the view.
    blocker = Component(:b, disk_mesh(2.0, h / 2; nr = 2, nphi = 16), gold(); two_sided = true)
    tb = Target(d1, d2, blocker)
    Fb = view_factors(tb)
    @test sum(Fb[i, j] for i in elements(tb, :d1), j in elements(tb, :d2)) == 0

    # Exact point-to-polygon factor: unit square at height 1 above a corner point,
    # against a fine midpoint quadrature of the point kernel.
    p, n = Vec3(0, 0, 0), Vec3(0, 0, 1)
    tri1 = (Vec3(0, 0, 1), Vec3(1, 0, 1), Vec3(1, 1, 1))
    tri2 = (Vec3(0, 0, 1), Vec3(1, 1, 1), Vec3(0, 1, 1))
    exact = _point_polygon_vf(p, n, tri1) + _point_polygon_vf(p, n, tri2)
    m = 400
    quad = sum(begin
        x, y = (i - 0.5) / m, (j - 0.5) / m
        r2 = x^2 + y^2 + 1
        1 / (π * r2^2) / m^2
    end for i in 1:m, j in 1:m)
    @test isapprox(exact, quad; rtol = 1e-4)
end

@testset "radiosity" begin
    # Integrating sphere with uniform source S and albedo α: B = H = S/(1-α).
    α, η = 0.7, 0.5
    s = Target(Component(:s, sphere_mesh(R_SPHERE; ntheta = 8, nphi = 12, outward = false),
                         Material(:m; albedo = α, conversion = η)))
    m = RadiationModel(s)
    P = 0.3 .* s.elem_area                      # 0.3 TW/mm² of absorbed laser
    sol = solve(m, P)
    expected = η * 0.3 / (1 - α)
    @test isapprox(sum(sol.H .* s.elem_area) / sum(s.elem_area), expected; rtol = 0.04)
    eb = energy_balance(sol)
    @test isapprox(eb.source, sum(values(eb.absorbed)) + eb.lost; rtol = 1e-10)
    @test abs(eb.lost) < 0.04 * eb.source
    # with enforced closure the integrating sphere is exact
mc = RadiationModel(s; closed = true)
@test maximum(abs.(closure_error(mc.F))) < 1e-9
@test maximum(abs.(s.elem_area .* mc.F .- (s.elem_area .* mc.F)')) < 1e-12
solc = solve(mc, P)
@test all(h -> isapprox(h, expected; rtol = 1e-6), solc.H)
# flux_response reproduces the direct solve
    obs = collect(1:5:nelements(s))
    R = flux_response(m, obs)
    @test R * P ≈ sol.H[obs]

    # Temperature-dependent albedo converges and stays in range.
    w = WallLossAlbedo()
    s2 = Target(Component(:s, sphere_mesh(R_SPHERE; ntheta = 8, nphi = 12, outward = false),
                          Material(:Au; albedo = w, conversion = 0.7)))
    sol2 = solve(RadiationModel(s2), P; time = 2.0)
    @test isapprox(sol2.albedo[1], Briareus.albedo(w, 2.0, radiation_temperature(sol2.H[1])); atol = 1e-3)
    @test 0 < sol2.albedo[1] < 0.95
    @test radiation_temperature(SIGMA_SB * 300.0^4) ≈ 300.0
end

@testset "beam profile" begin
    spot = SuperGaussianSpot(0.4, 0.25; order = 4.0, angle = 0.3)
    b = Beam("b"; theta = 30, phi = 40, power = 2.0, spot, fnumber = 8.0, aim = (0.1, -0.2, 0.3),
             defocus = 1.5)
    k = direction(b)
    e1, e2 = Briareus.spot_axes(b)
    @test abs(dot(e1, k)) < 1e-12 && abs(dot(e2, k)) < 1e-12 && abs(dot(e1, e2)) < 1e-12
    # power through planes perpendicular to the beam at several axial positions
    for s in (-3.0, 0.0, 1.5, 4.0)
        c = Vec3(b.aim) + s * k
        h = 0.02
        tot = sum(intensity(b, c + u * e1 + v * e2) for u in -2:h:2, v in -2:h:2) * h^2
        @test isapprox(tot, 2.0; rtol = 1e-3)
    end
end

@testset "deposition vs ray trace" begin
    # Normal incidence on a large foil: everything deposits, centroid at the aim point.
    tf = Target(foil(:f, (0, 0, 0), (0, 0, 1); width = 6.0, n = 40))
    b = Beam("b"; theta = 0, phi = 0, power = 3.0, spot = SuperGaussianSpot(0.5; order = 6.0),
             fnumber = 6.7, aim = (0.4, -0.3, 0.0), defocus = 2.0)
    op = DepositionOperator(tf, [b])
    dep = deposit(op, [b])
    rt = raytrace_deposition(tf, [b])
    @test isapprox(sum(dep.element), 3.0; rtol = 2e-3)
    @test isapprox(sum(rt.element), 3.0; rtol = 1e-9)
    cen(d) = sum(d.element .* tf.elem_centroid) / sum(d.element)
    @test isapprox(cen(dep), Vec3(0.4, -0.3, 0); atol = 0.01)
    # second moment of the defocused spot: analytic profile vs ray trace
    m2(d) = (c = cen(d); sum(d.element .* [sum(abs2, x - c) for x in tf.elem_centroid]) / sum(d.element))
    @test isapprox(m2(dep), m2(rt); rtol = 0.03)

    # Shadowing by an object near focus: a small disk in front of the foil.
    sh = Component(:shield, disk_mesh(0.3, 0.0; nr = 4, nphi = 24), gold(); two_sided = true)
    ts = Target(foil(:f, (0, 0, -2.0), (0, 0, 1); width = 6.0, n = 40), translate(sh, (0.2, 0, 0)))
    b2 = Beam("b"; theta = 0, phi = 0, power = 1.0, spot = SuperGaussianSpot(0.5; order = 4.0),
              fnumber = 6.7)
    op2 = DepositionOperator(ts, [b2]; ncone = 13)
    d2 = deposit(op2, [b2])
    r2 = raytrace_deposition(ts, [b2]; ncone = 13)
    pf(d) = sum(d.element[elements(ts, :f)])
    pshield(d) = sum(d.element[elements(ts, :shield)])
    @test isapprox(pshield(d2), pshield(r2); rtol = 0.03)
    @test isapprox(pf(d2), pf(r2); rtol = 0.03)
end

@testset "hohlraum + gradients" begin
    t = hohlraum(nz = 16, nphi = 24, capsule_ntheta = 12, capsule_nphi = 16)
    q = nif_quads()
    p = GroupPointing(q; fields = (:z, :r, :power))
    x0 = Float64[]
    for k in p.keys
        append!(x0, k in (:inner23, :inner30) ? [3.3, -0.8, 6.0] : [5.0, 0.0, 9.0])
    end
    sc = Scenario(t, apply(p, x0); observe = :capsule, parameterization = p)
    ev = evaluate(sc, x0)
    @test clipped_fraction(ev) < 0.01
    @test abs(missed_fraction(ev)) < 0.01
    # analytic deposition agrees with the ray trace on wall power and capsule symmetry
    rt = raytrace_deposition(t, ev.beams; nspot = 31)
    w = elements(t, :wall)
    @test isapprox(sum(ev.deposition.element[w]), sum(rt.element[w]); rtol = 0.01)
    a = flux_moments(ev; lmax = 4)
    sol = solve(sc.model, rt)
    b = legendre_moments(sol.H, t, sc.observed; lmax = 4)
    @test isapprox(a[1], b[1]; rtol = 0.01)
    @test abs(a[3] / a[1] - b[3] / b[1]) < 0.005
    # response operator equals a direct solve
    sol_a = solve(sc.model, ev.deposition)
    @test ev.flux ≈ sol_a.H[sc.observed]
    # top/bottom symmetric pointing ⇒ no odd moments
    @test abs(a[2] / a[1]) < 1e-6 && abs(a[4] / a[1]) < 1e-6
    # AD gradient matches central differences
    obj = symmetry_objective()
    f(x) = obj(evaluate(sc, x))
    g = ForwardDiff.gradient(f, x0)
    for i in (1, 2, 3, 10)
        h = 1e-5
        xp = copy(x0); xp[i] += h
        xm = copy(x0); xm[i] -= h
        @test isapprox(g[i], (f(xp) - f(xm)) / 2h; rtol = 1e-4, atol = 1e-9)
    end
    # energy accounting closes up to view-factor discretization
    eb = energy_balance(sol_a)
    @test abs(eb.lost) < 0.04 * eb.source
end

@testset "mode analysis" begin
    t = Target(Component(:s, sphere_mesh(1.0; ntheta = 40, nphi = 40), gold()))
    es = collect(1:nelements(t))
    μ = [c[3] / norm(c) for c in t.elem_centroid]
    v = 1 .+ 0.1 .* legendre.(2, μ) .- 0.03 .* legendre.(4, μ)
    a = legendre_moments(v, t, es; lmax = 6)
    @test isapprox(a[3] / a[1], 0.1; atol = 1e-3)
    @test isapprox(a[5] / a[1], -0.03; atol = 1e-3)
    @test abs(a[2]) < 1e-10
    C = ylm_moments(v, t, es; lmax = 4)
    mp = mode_power(C)
    @test isapprox(mp[3], 0.1 / sqrt(5); rtol = 0.02)
    # a tilted P2 has the same rotation-invariant power
    vt = 1 .+ 0.1 .* legendre.(2, [c[1] / norm(c) for c in t.elem_centroid])
    @test isapprox(mode_power(ylm_moments(vt, t, es; lmax = 4))[3], mp[3]; rtol = 0.03)
end

@testset "facilities" begin
    q = nif_quads()
    @test length(q) == 48
    @test length(nif_beams()) == 192
    angs = port_angles.(q)
    tops = filter(a -> a[1] < 90, angs)
    bots = filter(a -> a[1] > 90, angs)
    # every bottom port is (180-θ, 360-φ) of a top port
    @test all(any(isapprox(180 - b[1], a[1]; atol = 1e-6) && isapprox(mod(360 - b[2], 360), a[2]; atol = 1e-6)
                  for a in tops) for b in bots)
    @test count(b -> b.group === :inner23, q) == 8
    @test count(b -> b.group === :outer50, q) == 16
    ob = omega_beams(; cones = :all)
    @test length(ob) == 60
    θs = sort(unique(round.(first.(port_angles.(ob)); digits = 2)))
    @test θs ≈ [21.42, 42.02, 58.85, 81.25, 98.75, 121.15, 137.98, 158.58]
    @test length(omega_beams()) == 40
end

@testset "foil optimization + imaging" begin
    t = Target(foil(:bl, (0, 0, 0), (0, 0, 1); width = 3.0, n = 30, two_sided = false))
    beams = [Beam("b$i"; theta = 40, phi = 90i, power = 0.5, spot = omega_spot(:SG4),
                  fnumber = 6.7, group = Symbol(:b, i), aim = (0.3cos(i), 0.3sin(i), 0.0)) for i in 1:4]
    p = GroupPointing(beams; by = :group, fields = (:r, :phi, :z), mirror = false)
    sc = Scenario(t, beams; parameterization = p)
    @test sc.model === nothing
    obj = irradiance_objective(:bl; center = (0, 0, 0), radius = 0.4, spill = 0.0)
    x0 = initial_parameters(p)
    @test apply(p, x0)[1].aim ≈ beams[1].aim
    lo, hi = bounds(p, x0; r = 0.5, phi = 0.5, z = 0.0)
    res = optimize_pointing(sc, obj; x0, lower = lo, upper = hi, iterations = 50)
    @test res.value < 0.5 * res.initial_value
    dep = evaluate(sc, res.x).deposition
    img = render(t, laser_flux(dep, t), OrthoCamera((0, 0, -1); up = (0, 1, 0), width = 4.0, nx = 64))
    @test count(isfinite, img) ≈ 64^2 * (3 / 4)^2 rtol = 0.05
    @test maximum(filter(isfinite, img)) > 0
    # x-ray power toward a detector above the foil equals ΣB A cosθ / π
    B = ones(nelements(t))
    @test line_of_sight_power(t, B, (0, 0, 1)) ≈ 9.0 / π rtol = 1e-6
    @test line_of_sight_power(t, B, (0, 0, -1)) == 0
end

@testset "drive history" begin
    t = hohlraum(nz = 10, nphi = 16, capsule_ntheta = 8, capsule_nphi = 12,
                 wall = Material(:Au; albedo = WallLossAlbedo(), conversion = 0.8))
    q = nif_quads()
    q = [remake(b; aim = (0.0, 0.0, sign(b.port[3]) * 4.0)) for b in q]
    m = RadiationModel(t)
    op = DepositionOperator(t, q)
    pulse = PulseShape([0.0, 1.0, 5.0, 6.0], [0.0, 2.0, 8.0, 0.0])
    @test pulse_energy(pulse) ≈ 1 + 20 + 4
    h = drive_history(m, op, τ -> with_pulses(q, pulse, τ), [0.5, 2.0, 4.0, 5.5])
    @test issorted(h.Tr[1:3])
    @test all(0 .< h.wall_albedo .< 0.95)
    @test h.delivered ≈ 48 .* pulse.([0.5, 2.0, 4.0, 5.5])
end

end
