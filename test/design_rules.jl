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

# A mock-up hohlraum with hardware beside the top entrance hole: a one-sided shield plate
# facing away from the beams that pass it (its back side catches their light) and a
# two-sided stalk (a closed cylinder) across the hole. Shared by the testsets below.
const DR_TARGET = let
    th = hohlraum(length = 10.1, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9,
                  nz = 8, nphi = 16, capsule_ntheta = 8, capsule_nphi = 12)
    shield = foil(:shield, (2.4, 0.0, 7.0), (-1.0, 0.0, 0.0); width = 3.0, height = 2.0, n = 6,
                  two_sided = false)
    stalk = Component(:stalk, transform(cylinder_mesh(0.3, 0.0, 4.0; nz = 4, nphi = 12, inward = false),
                                        Briareus.rotation_to(Vec3(0, 1, 0)), Vec3(0.0, -4.5, 6.2)), gold())
    Target([th.components..., shield, stalk]; closed = th.closed)
end
# (z, r) per cone: the NIF quads through the entrance holes onto the wall
const DR_X = [3.09, -1.09, 3.63, -0.6, 5.2, 0.0, 5.01, 0.0]

@testset "design rules: power on named parts" begin
    t = DR_TARGET
    sc = Scenario(t, GroupPointing(nif_quads(); fields = (:z, :r)))
    ev = evaluate(sc, DR_X)
    P = ev.deposition.delivered
    on(e, n) = power_fraction(e, n) + clipped_fraction(e, n)
    @test clipped_fraction(ev, :shield) > 0.01 && power_fraction(ev, :shield) == 0
    @test power_fraction(ev, :stalk) > 0.01
    # the totals component_power reports, and every part's with the missed light closes
    cp = component_power(ev.deposition, t)
    @test clipped_fraction(ev, :shield) * P ≈ cp[:shield_outside] rtol = 1e-12
    @test [power_fraction(ev, :stalk), clipped_fraction(ev, :stalk)] .* P ≈ [cp[:stalk], cp[:stalk_outside]] rtol = 1e-12
    @test clipped_fraction(ev, (c.name for c in t.components)...) ≈ clipped_fraction(ev) rtol = 1e-12
    @test sum(n -> on(ev, n), parts(t)) + missed_fraction(ev) ≈ 1 rtol = 1e-12
    # against the reference ray trace (both sides of each component; the operator's
    # sampling of the spot edge at the hardware: 2e-3 and 8e-4 measured)
    rt = raytrace_component_power(t, ev.beams; nspot = 121)
    @test [on(ev, n) for n in (:shield, :stalk)] ≈
          [sum(rt[:, part_pieces(t, n)]) / P for n in (:shield, :stalk)] rtol = 5e-3
    # gradients against central differences (1e-9 relative measured)
    for n in (:shield, :stalk)
        f(y) = on(evaluate(sc, y), n)
        g = ForwardDiff.gradient(f, DR_X)
        h = 1e-5
        fd = [(f(DR_X + h * (eachindex(DR_X) .== i)) - f(DR_X - h * (eachindex(DR_X) .== i))) / 2h
              for i in eachindex(DR_X)]
        @test g ≈ fd rtol = 1e-6
    end
end
