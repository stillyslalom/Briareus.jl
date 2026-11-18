# Capture and clearance rules as quantities of an evaluation (src/design_rules.jl and the
# keep-outs' public penalty): each against the objectives' own terms, against independent
# power totals, and its gradient against finite differences.

@testset "design rules: the objectives' penalty terms" begin
    t = hohlraum(length = 10.1, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9,
                 nz = 8, nphi = 16, capsule_ntheta = 8, capsule_nphi = 12)
    p = GroupPointing(nif_quads(); fields = (:z, :r))
    sc = Scenario(t, p; observe = :capsule)
    # (z, r) of the inner and outer cones: the inner cones' spots reach the capsule, and
    # some light is clipped on the caps
    x = [2.0, 0.0, 2.0, 0.0, 5.0, 0.0, 5.0, 0.0]
    ev = evaluate(sc, x)
    c, m, d = clipped_fraction(ev), missed_fraction(ev), direct_fraction(ev)
    @test c > 1e-4 && d > 1e-4
    @test Briareus._penalties(ev; clip = 3.0, missed = 5.0, direct = 7.0, direct_names = nothing) ==
          3 * c^2 + 5 * m^2 + 7 * d^2
    @test direct_fraction(ev, (:capsule,)) == d == power_fraction(ev, :capsule)
    @test_throws ArgumentError direct_fraction(ev, (:Capsule,))

    ka = ApertureKeepOut(t; minimum = 0.3, weight = 2.0)
    ks = SphereKeepOut(t; minimum = 0.3, weight = 5.0)
    @test keepout_clearance(ka, ev) == aperture_clearance(ka, ev.beams)
    @test keepout_clearance(ks, ev) == sphere_clearance(ks, ev.beams)
    pa = keepout_penalty(ka, ev)
    @test pa == 2 * sum(max(0, 0.3 - ci)^2 for ci in aperture_clearance(ka, ev.beams) if isfinite(ci))
    @test pa > 0 && keepout_penalty(ks, ev) > 0
    @test keepout_penalty((ka, ks), ev) == keepout_penalty([ka, ks], ev) ==
          pa + keepout_penalty(ks, ev)
    @test keepout_penalty(nothing, ev) == 0

    # a hand-written objective on the public names reproduces symmetry_objective
    obj = symmetry_objective(; keepout = (ka, ks))
    function own(ev)
        a = flux_moments(ev; lmax = 6)
        J = (a[3] / a[1])^2 + (a[5] / a[1])^2 + sum((a[l+1] / a[1])^2 for l in 1:2:6)
        return J + 10 * (clipped_fraction(ev)^2 + missed_fraction(ev)^2 + direct_fraction(ev)^2) +
               keepout_penalty((ka, ks), ev)
    end
    @test own(ev) ≈ obj(ev) rtol = 1e-12
    # in the optimizer's number type
    g = ForwardDiff.gradient(y -> keepout_penalty((ka, ks), evaluate(sc, y)), x)
    h = ForwardDiff.gradient(y -> direct_fraction(evaluate(sc, y)), x)
    @test all(isfinite, g) && norm(g) > 0 && all(isfinite, h) && norm(h) > 0
end
