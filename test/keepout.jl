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
