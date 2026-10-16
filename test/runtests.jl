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
    top, rest = split_mesh(s, c -> c[3] > 1.5)
    @test npatches(top) + npatches(rest) == npatches(s)
    @test npatches(top) == 2 * 24                 # the two uppermost bands
    tt = Target(Component(:a, top, gold()), Component(:b, rest, gold()))
    @test sum(tt.elem_area) ≈ sum(t.elem_area)
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
    # quad names from the LMF chamber-geometry figures: a bijection onto the ports,
    # every bundle one top and one bottom quad, FABS on Q31B (30°) and Q36B (130°, 186°)
    @test length(unique(b.name for b in q)) == 48
    bundles = nif_bundles(q)
    @test length(bundles) == 24
    @test all(sort([q[i].name[end] for i in v]) == ['B', 'T'] for v in values(bundles))
    @test port_angles(q[findfirst(b -> b.name == "Q36B", q)])[1] ≈ 130
    @test port_angles(q[findfirst(b -> b.name == "Q36B", q)])[2] ≈ 185.62
    @test q[findfirst(b -> b.name == "Q31B", q)].group === :inner30
    @test nif_port("Q11T") == 17 && nif_port("Q11B") == 50
    kinds = Dict{Tuple{Bool,Bool},Int}()
    for v in values(bundles)
        inner = sort([q[i].group in (:inner23, :inner30) for i in v])
        kinds[(inner[1], inner[2])] = get(kinds, (inner[1], inner[2]), 0) + 1
    end
    @test kinds == Dict((false, false) => 12, (false, true) => 8, (true, true) => 4)
    @test all(nif_bundle(b) == nif_bundle(b.quad) for b in nif_beams())
    # the LMF per-beam table: each quad's mean beam direction is its port, which also
    # confirms the quad names read from the chamber-geometry figures
    nb = nif_beams()
    @test length(nb) == 192 && nb[1].name == "B111" && nb[1].quad == "Q11T"
    for b in q
        mean_dir = normalize(sum(x.port for x in nb if x.quad == b.quad))
        @test acosd(clamp(dot(mean_dir, b.port), -1, 1)) < 0.01
    end
    @test all(x.group === q[findfirst(y -> y.quad == x.quad, q)].group for x in nb)
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

@testset "clearance" begin
    @test enclosed_radius(2.0, 1 - exp(-1)) ≈ 1
    # point–ellipse distance against brute force, inside and outside
    rng = MersenneTwister(3)
    for _ in 1:200
        a, b = 0.3 + rand(rng), 0.3 + rand(rng)
        u, v = 3 .* (rand(rng, 2) .- 0.5)
        d = Briareus._ellipse_distance(u, v, a, b)
        bf = minimum(hypot(u - a * cos(θ), v - b * sin(θ)) for θ in range(0, 2π; length = 20001))
        @test abs(d) ≈ bf atol = 1e-5
        @test (d < 0) == ((u / a)^2 + (v / b)^2 < 1)
    end
    # Beam straight down onto a foil at z = 0; a small disk upstream at z = 2 with its
    # nearest edge 1.0 mm from the axis, and one below the foil (downstream, ignored).
    spot = SuperGaussianSpot(0.3; order = 4.0)
    b = Beam("b"; theta = 0, phi = 0, power = 1.0, spot, fnumber = 8.0)
    up = translate(Component(:up, disk_mesh(0.05, 0.0; nr = 1, nphi = 32), gold()), (1.05, 0, 2.0))
    down = translate(Component(:down, disk_mesh(0.05, 0.0; nr = 1, nphi = 32), gold()), (0.2, 0, -1.0))
    t = Target(foil(:f, (0, 0, 0), (0, 0, 1); width = 4.0, n = 8), up, down)
    r = clearance(t, [b]; targets = (:f,), enclosed = 0.99, margin = 0.1, pointing = 0.05,
                  target_offset = 0.02, hardware_tolerance = 0.03, focus = 0.4, wings = [])
    row(c) = only(x for x in r.rows if x.component === c)
    nominal = 1.0 - (enclosed_radius(4.0, 0.99) * 0.3 + 2.0 / 16 + 0.1)
    @test row(:up).nominal ≈ nominal atol = 1e-9
    @test row(:up).clearance ≈ nominal - (0.05 + 0.02 + 0.03) - 0.4 / 16 atol = 1e-9
    @test row(:up).s ≈ -2.0
    @test row(:down).clearance == Inf
    @test isempty(violations(r))
    # intensity contours: the peak-fraction and absolute-intensity envelopes agree
    I0 = focal_profile(spot, 1.0, 0.0)
    ρ = envelope_radius(spot, 1.0; peak_fraction = 1e-6)
    @test ρ ≈ log(1e6)^(1 / 4) rtol = 1e-9
    @test envelope_radius(spot, 1.0; intensity = 1e-6 * I0 * 1e14) ≈ ρ rtol = 1e-9
    @test envelope_radius(spot, 1.0; intensity = 2 * I0 * 1e14) == 0
    # a halo widens the envelope far more than the core alone
    hspot = SuperGaussianSpot(0.3; order = 4.0, halo = 1e-3, halo_scale = 4.0)
    @test envelope_radius(hspot, 1.0; peak_fraction = 1e-6) > 3 * ρ
    @test enclosed_fraction(hspot, 50.0) ≈ 1
    # the halo carries its power through deposition and the ray trace
    hb = remake(b; spot = hspot)
    tf = Target(foil(:f, (0, 0, 0), (0, 0, 1); width = 12.0, n = 60))
    @test sum(raytrace_deposition(tf, [hb]).element) ≈ 1 rtol = 1e-9
    @test sum(deposit(DepositionOperator(tf, [hb]), [hb]).element) ≈ 1 rtol = 2e-3
    # pointing toward hardware raises the worst-case power and intensity on it
    near = translate(Component(:near, disk_mesh(0.2, 0.0; nr = 2, nphi = 32), gold()), (0.75, 0, 1.0))
    tn = Target(foil(:f, (0, 0, 0), (0, 0, 1); width = 4.0, n = 8), near)
    nowing = (; wings = [], target_offset = 0.0, focus = 0.0)
    r0 = only(clearance(tn, [b]; targets = (:f,), pointing = 0.0, nowing...).rows)
    rep1 = clearance(tn, [b]; targets = (:f,), pointing = 0.1, nowing...)
    r1 = only(rep1.rows)
    @test r1.power > r0.power && r1.peak_intensity > r0.peak_intensity
    @test r1.clearance ≈ r0.clearance - 0.1
    strict = (fail_intensity = Inf, fail_fraction = Inf, review_intensity = 0.5 * r1.peak_intensity,
              review_fraction = Inf)
    rs = clearance(tn, [b]; targets = (:f,), pointing = 0.1, thresholds = strict, nowing...)
    @test only(rs.rows).status === :review
    @test length(violations(rs)) == 1 && isempty(violations(rs; status = (:fail,)))
    # wings: a large foil catches all of the halo and scatter power; the far scatter
    # reaches hardware the core envelope clears by millimetres
    wl = nif_wings()
    @test sum(w.fraction for w in wl) ≈ 0.011
    tw = Target(foil(:f, (0, 0, 0), (0, 0, 1); width = 60.0, n = 30),
                translate(Component(:far, disk_mesh(0.5, 0.0; nr = 2, nphi = 16), gold()), (8.0, 0, 3.0)))
    rw = only(clearance(tw, [b]; targets = (:f,), wings = wl, focus = 0.0, intensity = 1e7).rows)
    rn = only(clearance(tw, [b]; targets = (:f,), wings = [], focus = 0.0, intensity = 1e7).rows)
    @test rn.power == 0 && rw.power > 0
    @test rw.peak_intensity ≈ 1e14 * 1e-3 / (π * 23.0^2) rtol = 1e-6   # scatter floor only
    @test rw.clearance < 0 < rn.clearance        # inside the 1e7 W/cm² scatter disk
    Lw = zeros(2)
    Briareus._raytrace_wing!(Lw, tw, b, wl[1], Briareus._wing_rh(wl[1], b))
    @test Lw[1] ≈ 0.01 rtol = 1e-9
    # A plate across the beam upstream intercepts everything and overlaps.
    plate = Component(:plate, translate(rectangle_mesh(3.0, 3.0; nu = 6, nv = 6), (0, 0, 1.0)), gold();
                      two_sided = true)
    r2 = clearance(Target(foil(:f, (0, 0, 0), (0, 0, 1); width = 4.0, n = 8), plate), [b]; targets = (:f,),
                   wings = [])
    @test only(r2.rows).clearance < 0
    @test only(r2.rows).intercepted ≈ 1 rtol = 1e-3
    @test only(r2.rows).status === :fail
    @test length(violations(r2)) == 1
    path = write_clearance_csv(tempname() * ".csv", r)
    @test countlines(path) == 1 + length(r.rows)
end

@testset "I/O" begin
    q = [remake(b; aim = (0.1, -0.2, 0.3), dlambda = 1.5, defocus = 0.7) for b in nif_quads()]
    path = write_beam_table(tempname() * ".csv", q)
    r = read_beam_table(path)
    @test length(r) == length(q)
    @test all(a.name == b.name && a.port ≈ b.port && a.aim ≈ b.aim && a.defocus ≈ b.defocus &&
              a.spot.a ≈ b.spot.a && a.spot.b ≈ b.spot.b && a.spot.order == b.spot.order &&
              a.spot.angle ≈ b.spot.angle && a.group === b.group && a.dlambda ≈ b.dlambda &&
              a.fnumber == b.fnumber for (a, b) in zip(q, r))
    # a foreign, semicolon-separated pointing table in cm updates a template beam set
    f = tempname() * ".csv"
    write(f, "# exported pointing\nBeam;X (cm);Z (cm)\nQ15T;0.01;0.42\n\"Q13T\";-0.02;0.36\n")
    fb = read_beam_table(f; columns = Dict(:name => "Beam", :aim_x => "X (cm)", :aim_z => "Z (cm)"),
                         scale = Dict(:aim_x => 10, :aim_z => 10), template = q)
    @test fb[1].aim ≈ Vec3(0.1, -0.2, 4.2)
    @test fb[2].aim ≈ Vec3(-0.2, -0.2, 3.6)
    @test fb[1].port == q[1].port && fb[1].spot.a == q[1].spot.a
    @test_throws ArgumentError read_beam_table(f)            # no name column without mapping
    t = hohlraum(nz = 4, nphi = 8, capsule_ntheta = 4, capsule_nphi = 6)
    p2 = write_element_table(tempname() * ".csv", t, :one => ones(nelements(t)))
    hdr, rows = read_csv(p2)
    @test hdr[end] == "one" && length(rows) == nelements(t)
    p3 = write_moments(tempname() * ".csv", [2.0, 0.0, -0.2])
    @test read_csv(p3)[2][3][3] == "-0.1"
end

@testset "laser-plasma models" begin
    beams = [Beam("i$k"; theta = 23.5, phi = 90k, power = 2.0, group = :inner) for k in 1:4]
    append!(beams, [Beam("o$k"; theta = 50, phi = 90k + 45, power = 3.0, group = :outer) for k in 1:4])
    bs = apply_lpi(Backscatter(:inner => 0.1, :outer => ([0.0, 10.0], [0.0, 0.2])), beams, 5.0)
    @test all(b.power ≈ 1.8 for b in bs if b.group === :inner)
    @test all(b.power ≈ 2.7 for b in bs if b.group === :outer)
    # calibrated transfer conserves power and responds to Δλ
    ct = CalibratedTransfer((:outer,), (:inner,), linear_transfer(0.05, 0.1))
    shifted = [remake(b; dlambda = b.group === :inner ? 1.0 : 0.0) for b in beams]
    eff = apply_lpi(ct, shifted)
    @test sum(b.power for b in eff) ≈ sum(b.power for b in beams)
    @test sum(b.power for b in eff if b.group === :outer) ≈ 12 * (1 - 0.15)
    chain = apply_lpi((Backscatter(:inner => 0.1), ct), shifted)
    @test sum(b.power for b in chain) ≈ 20 - 0.8
    tab = tabulated_transfer([-1.0, 0.0, 2.0], [0.0, 0.1, 0.3])
    @test tab(1.0, 0.0) ≈ 0.2
    @test ForwardDiff.derivative(x -> tab(x, 0.0), 1.0) ≈ 0.1
    # gradients flow through dlambda in a scenario
    t = hohlraum(nz = 10, nphi = 16, capsule_ntheta = 8, capsule_nphi = 12)
    q = [remake(b; aim = (0.0, 0.0, sign(b.port[3]) * (b.group in (:inner23, :inner30) ? 3.3 : 5.0)))
         for b in nif_quads()]
    p = GroupPointing(q; fields = (:dlambda,))
    lpi = CalibratedTransfer((:outer44, :outer50), (:inner23, :inner30), linear_transfer(0.0, 0.05))
    sc = Scenario(t, q; observe = :capsule, parameterization = p, lpi)
    f(x) = (a = flux_moments(evaluate(sc, x); lmax = 2); a[3] / a[1])
    x0 = zeros(nparams(p))
    x1 = [k in (:inner23, :inner30) ? 1.0 : 0.0 for k in p.keys]
    @test f(x1) < f(x0)               # more inner-cone (waist) drive lowers flux P2
    g = ForwardDiff.gradient(f, x1)
    h = 1e-6
    @test g[1] ≈ (f(x1 .+ h .* (1:4 .== 1)) - f(x1 .- h .* (1:4 .== 1))) / 2h rtol = 1e-4
    ev = evaluate(sc, x1)
    @test sum(b.power for b in ev.effective) ≈ sum(b.power for b in ev.beams)
end


@testset "CBET" begin
    @test plasma_dispersion(1.0) ≈ complex(-1.0761, 0.6520) atol = 1e-4
    @test plasma_dispersion(0.0) ≈ complex(0.0, sqrt(π))
    # independent cross-check (research notes): NIF-LEH-like CH plasma, 20° crossing,
    # 351 nm: |Im Γ| peaks ≈ 153 near Δλ(3ω) ≈ 2 Å without flow
    m = CBET(leh_z = 5.0, mach = 0.0)
    pl = cbet_plasma(m, 0.351, 0.0)
    k = 2 * (2π / 0.351e-6) * sqrt(1 - 0.06) * sind(10)
    Ω(dl) = 2π * 2.99792458e8 * dl * 1e-10 / 0.351e-6^2
    ims = [abs(imag(cbet_coupling(Ω(dl), k, 0.0; ne = pl.ne, Te = 5000.0, Ti = 1500.0, ions = CH_IONS)))
           for dl in 0:0.02:4]
    @test maximum(ims) ≈ 153 rtol = 0.02
    @test 1.8 < (0:0.02:4)[argmax(ims)] < 2.1
    # the redder beam gains; coupling is antisymmetric
    mk(φ, dl) = Beam("x"; theta = 10, phi = φ, power = 1.0, dlambda = dl, fnumber = 8.0,
                     spot = SuperGaussianSpot(0.5; order = 4.0), aim = (0, 0, 5.0))
    @test cbet_gain(m, mk(0, 1.5), mk(180, 0.0)) > 0
    @test cbet_gain(m, mk(0, 1.5), mk(180, 0.0)) ≈ -cbet_gain(m, mk(180, 0.0), mk(0, 1.5))
    @test cbet_gain(m, mk(0, 0.0), mk(180, 0.0)) == 0          # no detuning, no flow
    # hohlraum: power conserved, transfer to the redder inner cones, differentiable
    q = [remake(b; aim = (0.0, 0.0, sign(b.port[3]) * (b.group in (:inner23, :inner30) ? 3.3 : 5.0)))
         for b in nif_quads()]
    inner = (:inner23, :inner30)
    mf = CBET(leh_z = 5.05)
    shifted = [remake(b; dlambda = b.group in inner ? 1.0 : 0.0) for b in q]
    eff = apply_lpi(mf, shifted)
    @test sum(b.power for b in eff) ≈ sum(b.power for b in q) rtol = 1e-12
    @test receiver_ratio(mf, q, 2.0; receivers = inner) > receiver_ratio(mf, q, 0.5; receivers = inner)
    f(dl) = receiver_ratio(mf, q, dl; receivers = inner)
    h = 1e-5
    @test ForwardDiff.derivative(f, 1.0) ≈ (f(1.0 + h) - f(1.0 - h)) / 2h rtol = 1e-5
    @test ion_wave_amplitude(mf, eff[1], eff[end]) ≥ 0
    # calibration recovers a known multiplier from synthetic transfer data
    truth = CBET(leh_z = 5.05, multiplier = 0.4, clamp = nothing)
    data = [(dl, receiver_ratio(truth, q, dl; receivers = inner)) for dl in (0.5, 1.0, 2.0)]
    fitted, res = calibrate_cbet(CBET(leh_z = 5.05, clamp = nothing), q, data; receivers = inner)
    @test fitted.multiplier ≈ 0.4 rtol = 1e-3
end


@testset "unconverted light" begin
    b = Beam("b"; theta = 0, phi = 0, power = 1.0)          # e1 = θ̂ = x̂ at the pole
    plate = Component(:plate, translate(rectangle_mesh(2.0, 2.0; nu = 2, nv = 2), (10.0, 0, 0)), gold();
                      two_sided = true)
    near = Component(:near, translate(sphere_mesh(0.5; ntheta = 6, nphi = 8), (1.0, 0, 0)), gold())
    t = Target(plate, near)
    r = unconverted_light_check(t, [b])
    row(c, col) = only(x for x in r.rows if x.component === c && x.color === col)
    @test row(:plate, :omega1).exposure === :possible       # some offset direction reaches it
    @test row(:near, :omega1).exposure === :clear           # inside the 4.8 mm clear radius
    @test row(:near, :omega2).exposure === :clear
    @test r.flat[:plate] && !r.flat[:near]
    @test r.zone[:near] ≈ 0.5 rtol = 1e-6
    # with a known offset direction the footprint is a definite square
    rx = unconverted_light_check(t, [b]; offsets = Dict("b" => 0.0))      # toward +x
    ry = unconverted_light_check(t, [b]; offsets = Dict("b" => 90.0))     # toward +y
    @test only(x for x in rx.rows if x.component === :plate && x.color === :omega1).exposure === :exposed
    @test only(x for x in ry.rows if x.component === :plate && x.color === :omega1).exposure === :clear
end


@testset "spherical harmonics and 3D symmetry" begin
    t = Target(Component(:s, sphere_mesh(1.0; ntheta = 40, nphi = 80), gold()))
    es = collect(1:nelements(t))
    W = ylm_weights(t, es; lmax = 4)
    @test size(W) == (25, length(es))
    # x and y components land in m = +1 and m = -1 (azimuth from +x)
    u = [normalize(c) for c in t.elem_centroid]
    cx = W * [1 + 0.1 * v[1] for v in u]
    cy = W * [1 + 0.1 * v[2] for v in u]
    @test cx[mode_index(1, 1)] > 0 && abs(cx[mode_index(1, -1)]) < 1e-6
    @test cy[mode_index(1, -1)] > 0 && abs(cy[mode_index(1, 1)]) < 1e-6
    # Y_11 = √(3/4π) x, so 0.1 x has c_11 = 0.1 √(4π/3) relative to c_00 = √(4π)
    @test cx[mode_index(1, 1)] / cx[1] ≈ 0.1 / sqrt(3) rtol = 1e-3
    # axisymmetric field: relative_mode_power = (a_l/a_0)²/(2l+1)
    μ = [v[3] for v in u]
    f = 1 .+ 0.1 .* legendre.(2, μ)
    c = W * f
    @test relative_mode_power(c, 2) ≈ 0.1^2 / 5 rtol = 5e-3
    @test ylm_matrix(c, 4) ≈ ylm_moments(f, t, es; lmax = 4)

    # NIF hohlraum: the quad set is D2-symmetric, so cos-type modes with l+m odd and
    # sin-type modes with l+m even vanish; dropping a bundle breaks it
    th = hohlraum(nz = 16, nphi = 24, capsule_ntheta = 12, capsule_nphi = 24)
    q = [remake(b; aim = (0.0, 0.0, sign(b.port[3]) * (b.group in (:inner23, :inner30) ? 3.3 : 5.0)))
         for b in nif_quads()]
    p = GroupPointing(q; by = :quad, fields = (:power,))
    sc = Scenario(th, q; observe = :capsule, parameterization = p)
    cf = flux_ylm(evaluate(sc, q); lmax = 4)
    forbidden = [(l, m) for l in 0:4 for m in -l:l if (m >= 0 ? isodd(l + m) : iseven(l + m))]
    @test maximum(abs(cf[mode_index(l, m)]) for (l, m) in forbidden) < 1e-10 * cf[1]
    keep = findall(b -> nif_bundle(b) != "11", q)
    sub = subset(sc, keep)
    cd = flux_ylm(evaluate(sub, q[keep]); lmax = 4)
    @test relative_mode_power(cd, 1) > 1e-6
    # subset ≡ zero power on the dropped quads in the full scenario
    q0 = [nif_bundle(b) == "11" ? remake(b; power = 0.0) : b for b in q]
    @test flux_ylm(evaluate(sc, q0); lmax = 4) ≈ cd rtol = 1e-12
    # influence matrix: modes are linear in quad powers
    M = influence_matrix(sc, q; lmax = 4)
    P = [b.power for b in q0]
    @test M * P ≈ cd rtol = 1e-10
    # dropping also the bundle related to 11 by the 180° rotation about z restores that
    # symmetry, which removes every odd-m mode (l = 1, m = 0 is allowed and stays tiny)
    p11 = q[findfirst(x -> x.name == "Q11T", q)].port
    partner = nif_bundle(q[findfirst(b -> b.port ≈ Vec3(-p11[1], -p11[2], p11[3]), q)])
    @test partner == "34"
    c2 = flux_ylm(evaluate(sc, [nif_bundle(b) in ("11", partner) ? remake(b; power = 0.0) : b
                                for b in q]); lmax = 4)
    @test maximum(abs(c2[mode_index(l, m)]) for l in 1:4 for m in -l:l if isodd(m)) < 1e-10 * c2[1]
    @test abs(c2[mode_index(1, 0)]) < 1e-5 * c2[1]
    # objective gradient and absolute power caps
    obj = asymmetry_objective(components = ((1, 1) => 5.0,))
    x0 = initial_parameters(sub.param)
    g = ForwardDiff.gradient(x -> obj(evaluate(sub, x)), x0)
    h = 1e-6
    e1 = zeros(length(x0)); e1[3] = 1
    @test g[3] ≈ (obj(evaluate(sub, x0 .+ h .* e1)) - obj(evaluate(sub, x0 .- h .* e1))) / 2h rtol = 1e-4
    lo, hi = bounds(sub.param, x0; power = 0.5, power_max = 9.0)
    @test all(hi .== 9.0) && all(lo .== 4.0)
end

end
