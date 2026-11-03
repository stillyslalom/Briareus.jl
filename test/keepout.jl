# Keep-outs against direct geometry: a round spot through a round entrance hole and past a
# sphere, where the clearance is a distance between circles in the beam's normal plane.

@testset "keep-outs against direct geometry" begin
    # the envelope of a round super-Gaussian exp(-ρⁿ) at a fraction f of its peak:
    # ρ = (ln 1/f)^(1/n) 1/e radii
    sp = SuperGaussianSpot(0.3; order = 4.0)
    ρ = log(100)^(1 / 4)
    @test envelope_radius(sp, 1.0; peak_fraction = 0.01) ≈ ρ rtol = 1e-10
    F = 8.0
    cone = 1 / (2F)
    tr = 0.05 + 0.025                                   # pointing + target offset (defaults)

    # a beam down the axis of the top entrance hole (radius 1.55 at z = 5.05), aimed off
    # center in the hole's plane: the nearest rim point is R − |x₀| from the axis
    th = hohlraum(length = 10.1, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9,
                  nz = 6, nphi = 16, capsule_ntheta = 6, capsule_nphi = 8)
    ko = ApertureKeepOut(th)
    beam(x0, d = 0.0; theta = 0) = Beam("b"; theta, phi = 0, power = 1.0, spot = sp, fnumber = F,
                                        aim = (x0, 0.0, 5.05), defocus = d)
    x0 = [0.0, 0.2, -0.4, 0.7]
    ref = [1.55 - abs(x) - ρ * 0.3 - tr - 1.0 * cone for x in x0]
    # (the rim is a polygon on the circle: within R(1 − cos(π/96)) ≈ 1e-3 mm)
    @test aperture_clearance(ko, beam.(x0)) ≈ ref atol = 2e-3
    # defocus moves best focus off the hole's plane by d: the f-cone grows the envelope by d/2F
    @test aperture_clearance(ko, [beam(0.2, 1.5)]) ≈ [ref[2] - 1.5 * cone] atol = 2e-3
    # margin, error radii and the focus allowance subtract linearly
    km = ApertureKeepOut(th; margin = 0.1, pointing = 0.0, target_offset = 0.0, hardware_tolerance = 0.02,
                         focus = 2.0)
    @test aperture_clearance(km, [beam(0.2)]) ≈ [1.55 - 0.2 - ρ * 0.3 - 0.1 - 0.02 - 2cone] atol = 2e-3
    # a beam whose port faces neither hole has no entrance aperture
    @test aperture_clearance(ko, [beam(0.0; theta = 90)]) == [Inf]

    # a sphere of radius 0.9 at the origin, a foil behind it: a beam down −z whose axis
    # passes the center at h has clearance h − (envelope + R)
    cap = Component(:capsule, sphere_mesh(0.9; ntheta = 6, nphi = 8), gold(); surface = SphereSurface(0.9))
    ts = Target(cap, foil(:f, (0, 0, -3), (0, 0, 1); width = 8.0, n = 4))
    ks = SphereKeepOut(ts)
    @test ks.names == [:capsule]
    sbeam(h, d = 0.0) = Beam("s"; theta = 0, phi = 0, power = 1.0, spot = sp, fnumber = F,
                             aim = (h, 0.0, 0.0), defocus = d)
    h = [1.0, 1.4, 2.0]
    @test sphere_clearance(ks, sbeam.(h)) ≈ [x - ρ * 0.3 - 0.9 - tr - cone for x in h] rtol = 1e-12
    # best focus 1 mm past the center: the envelope there has grown by 1/2F
    @test sphere_clearance(ks, [sbeam(1.4, 1.0)]) ≈ [1.4 - ρ * 0.3 - cone - 0.9 - tr - cone] rtol = 1e-12
    # a foil between the beam and the sphere shadows it: no clearance to keep
    tb = Target(cap, foil(:f, (0, 0, 2), (0, 0, 1); width = 8.0, n = 4))
    @test sphere_clearance(SphereKeepOut(tb), [sbeam(1.4)]) == [Inf]
    # the penalty: weight × Σ max(0, minimum − c)²
    kp = SphereKeepOut(ts; minimum = 0.5, weight = 3.0)
    c = sphere_clearance(kp, sbeam.(h))
    @test Briareus._keepout_penalty(kp, sbeam.(h)) ≈ 3 * sum(max(0, 0.5 - x)^2 for x in c) rtol = 1e-12
end

# The aperture keep-out on holes of any shape, against direct geometric computations:
# closed-form distances from a beam axis to a hole's edges at normal incidence, and the
# minimum over densely sampled analytic rims (round spots, so the envelope distance is
# exact) for oblique, defocused beams.

@testset "aperture keep-out: hole shapes" begin
    z0 = 1.0
    floor = Component(:floor, transform(rectangle_mesh(8.0, 8.0; nu = 2, nv = 2), Matrix(1.0I, 3, 3), Vec3(0, 0, -2)), BLACK)
    hole(m; surface = nothing) = Target(floor, aperture(:hole, transform(m, Matrix(1.0I, 3, 3), Vec3(0, 0, z0)); surface))
    ko(t) = ApertureKeepOut(t; pointing = 0.0, target_offset = 0.0, focus = 0.0)
    ρ(b) = Briareus.envelope_radius(b.spot, b.power; peak_fraction = 0.01)
    # the clearance from rim points `xs` of a round-spot beam, computed directly
    function direct(b, xs)
        k = direction(b)
        f = b.aim + b.defocus * k
        r = ρ(b) * b.spot.a
        return minimum(norm((x - f) - dot(x - f, k) * k) - r - abs(dot(x - f, k)) / (2b.fnumber) for x in xs)
    end
    seg(p, q, n) = (p + (q - p) * i / n for i in 0:n)
    polyline(vs, n) = [x for i in eachindex(vs) for x in seg(vs[i], vs[mod1(i + 1, length(vs))], n)]
    # the distance from a point to a segment, in closed form
    segdist(x, p, q) = (τ = clamp(dot(x - p, q - p) / dot(q - p, q - p), 0, 1); norm(x - p - τ * (q - p)))

    # a 2 × 1 mm rectangle, with a beam down its normal and one oblique and defocused
    W, H = 2.0, 1.0
    tr = hole(rectangle_mesh(W, H; nu = 3, nv = 2))
    kr = ko(tr)
    @test length(kr.rims[1]) == 10 && all(!r.arc for r in kr.rims[1])     # the boundary edges
    down = Beam("down"; theta = 0, phi = 0, spot = SuperGaussianSpot(0.15), aim = (0.3, 0.1, z0))
    @test only(aperture_clearance(kr, [down])) ≈ min(W / 2 - 0.3, H / 2 - 0.1) - ρ(down) * 0.15 rtol = 1e-12
    # the rim at the binding edge moves with the aim: ∂c/∂y = −1, ∂c/∂x = 0
    g = ForwardDiff.gradient(a -> only(aperture_clearance(kr, [remake(down; aim = (a[1], a[2], z0))])), [0.3, 0.1])
    @test g ≈ [0.0, -1.0] atol = 1e-7       # the refined rim point, to ~1e-9 mm
    corners = [Vec3(-W / 2, -H / 2, z0), Vec3(W / 2, -H / 2, z0), Vec3(W / 2, H / 2, z0), Vec3(-W / 2, H / 2, z0)]
    tilt = Beam("tilt"; theta = 35, phi = 20, spot = SuperGaussianSpot(0.2), aim = (0.2, -0.15, z0 - 0.5),
                defocus = 0.8, fnumber = 8.0)
    @test only(aperture_clearance(kr, [tilt])) ≈ direct(tilt, polyline(corners, 50_000)) atol = 1e-9

    # an ellipse, semi-axes 1.5 and 0.8, meshed as a 64-gon (no `DiskSurface`: the mesh is
    # the geometry): the polygon's edges, and within their sag of the ellipse
    a, b, n = 1.5, 0.8, 64
    te = hole(transform(disk_mesh(1.0, 0.0; nr = 3, nphi = n), Diagonal([a, b, 1.0]), zero(Vec3)))
    ke = ko(te)
    vs = [Vec3(a * cos(2π * j / n), b * sin(2π * j / n), z0) for j in 0:n-1]
    @test length(ke.rims[1]) == n
    x = Vec3(0.2, 0.1, z0)
    de = Beam("down"; theta = 0, phi = 0, spot = SuperGaussianSpot(0.15), aim = x)
    gap = minimum(segdist(x, vs[i], vs[mod1(i + 1, n)]) for i in 1:n)
    ce = only(aperture_clearance(ke, [de]))
    @test ce ≈ gap - ρ(de) * 0.15 rtol = 1e-12
    @test 0 < direct(de, [Vec3(a * cos(θ), b * sin(θ), z0) for θ in range(0, 2π; length = 200_000)]) - ce <
          a * (1 - cos(π / n))
    @test only(aperture_clearance(ke, [tilt])) ≈ direct(tilt, polyline(vs, 2_000)) atol = 1e-9

    # a circle on its `DiskSurface` is the circle itself, whatever the mesh: one disk of
    # the full radius, and a part made of a disk inside a ring of twice its count (the
    # seam's arcs meet, their chords leave a crack: inside the hole either way)
    disk = DiskSurface(center = (0, 0, z0))
    piece(name, m) = Component(name, transform(m, Matrix(1.0I, 3, 3), Vec3(0, 0, z0)), BLACK; kind = :aperture,
                               surface = disk, part = :leh)
    cut = Target(floor, piece(:inner, disk_mesh(0.6, 0.0; nr = 2, nphi = 16)),
                 piece(:outer, annulus_mesh(0.6, 1.0, 0.0; nr = 2, nphi = 32)))
    kc = ko(cut)
    whole = ko(hole(disk_mesh(1.0, 0.0; nr = 3, nphi = 24); surface = disk))
    @test kc.names == [:leh] && length(kc.rims[1]) == 32 && all(r -> r.arc && r.ρ ≈ 1.0, kc.rims[1])
    @test all(r -> r.arc && r.ρ ≈ 1.0, whole.rims[1])
    circle = [Vec3(cos(θ), sin(θ), z0) for θ in range(0, 2π; length = 400_000)]
    cw = only(aperture_clearance(whole, [tilt]))
    @test cw ≈ direct(tilt, circle) atol = 1e-9
    @test only(aperture_clearance(kc, [tilt])) ≈ cw atol = 1e-9
end
