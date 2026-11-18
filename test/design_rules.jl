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

@testset "design rules: hardware keep-out against clearance" begin
    t = DR_TARGET
    beams = apply(GroupPointing(nif_quads(); fields = (:z, :r)), DR_X)
    rep = clearance(t, beams; targets = (:wall, :capsule), hardware = (:shield, :stalk),
                    peak_fraction = 0.01, wings = [], nangles = 4, rtol = 0, resolution = 1e-7, nspot = 5)
    worst(part) = [only(r.clearance for r in rep.rows if r.beam == b.name && r.component === part) for b in beams]
    ks = [HardwareKeepOut(t, part; targets = (:wall, :capsule)) for part in (:shield, :stalk)]
    c = [hardware_clearance(k, beams) for k in ks]
    w = [worst(:shield), worst(:stalk)]
    @test any(<(0), c[1]) && any(>(0), c[1]) && any(<(0), c[2]) && any(>(0), c[2])
    # where both judge the same points (clearance skips those hidden from the beam, so
    # it reads higher above 1.2 mm on the shield, and deep inside the envelope its
    # point–ellipse distance loses digits next to the major axis): within its 1e-7 mm
    # resolution. Everywhere the keep-out is the lower: conservative
    near = [findall(x -> -0.5 < x < 1.0, wi) for wi in w]
    @test sum(length, near) >= 12
    @test maximum(maximum(abs.(c[j][near[j]] .- w[j][near[j]])) for j in 1:2) < 1e-6
    @test all(all(isfinite.(c[j]) .>= isfinite.(w[j])) for j in 1:2)
    @test all(all(c[j][i] <= w[j][i] + 1e-4 for i in eachindex(beams) if isfinite(w[j][i])) for j in 1:2)

    # derivatives in each beam's aim, defocus and spot axes against central differences,
    # deep inside the envelope (the medial strip) included; not where the closest point
    # lies on the landing plane's cut (one bottom beam, 10 mm from the stalk behind the
    # hohlraum), whose motion the derivatives leave out
    err = 0.0
    for k in ks, b in beams
        only(hardware_clearance(k, [b])) < 3 || continue
        mk(y) = remake(b; aim = (y[1], y[2], y[3]), defocus = y[4], spot = with_axes(b.spot, y[5], y[6]))
        y0 = [b.aim..., b.defocus, b.spot.a, b.spot.b]
        f(y) = only(hardware_clearance(k, [mk(y)]))
        g = ForwardDiff.gradient(f, y0)
        fd = [(f(y0 + 1e-6 * (1:6 .== j)) - f(y0 - 1e-6 * (1:6 .== j))) / 2e-6 for j in 1:6]
        err = max(err, maximum(abs.(g .- fd)))
    end
    @test err < 1e-6

    # the keep-out's penalty, alone and beside another keep-out, and its quantity on an
    # evaluation
    kp = HardwareKeepOut(t, (:shield, :stalk); targets = (:wall, :capsule), minimum = 0.2, weight = 3.0)
    cp = hardware_clearance(kp, beams)
    @test all(isapprox.(cp, min.(c[1], c[2]); atol = 1e-12))
    @test keepout_penalty(kp, beams) ≈ 3 * sum(max(0, 0.2 - x)^2 for x in cp if isfinite(x)) rtol = 1e-12
    ka = ApertureKeepOut(t)
    @test keepout_penalty((kp, ka), beams) ≈ keepout_penalty(kp, beams) + keepout_penalty(ka, beams) rtol = 1e-12
    ev = evaluate(Scenario(t, GroupPointing(nif_quads(); fields = (:z, :r))), DR_X)
    @test keepout_clearance(kp, ev) == cp
end

@testset "design rules: TaLIS 3ω clearance as a keep-out" begin
    # closed form: two square plates along a beam down −z, one near TCC and one far from
    # it; each plate's minimum lies at the middle of its far edge (the f-cone and the
    # tolerance both grow upstream), and the closest point jumps from the far plate to
    # the near one as the aim crosses x = 3·cone
    near = foil(:near, (1.0, 0.0, 2.0), (-1.0, 0.0, 0.0); width = 2.0, n = 2)
    far = foil(:far, (-1.0, 0.0, 8.0), (1.0, 0.0, 0.0); width = 2.0, n = 2)
    tp = Target([near, far])
    kt = HardwareKeepOut(tp, (:near, :far); targets = (), tolerance = 0.5, miss = :all)
    b0 = Beam("b"; theta = 0, phi = 0, spot = SuperGaussianSpot(0.2), aim = (0.0, 0.0, 0.0))
    R = envelope_radius(b0.spot, b0.power; peak_fraction = 0.01) * b0.spot.a
    cone, tα = 1 / (2b0.fnumber), tand(0.5)
    xs = -0.1:0.001:0.1
    c = hardware_clearance(kt, [remake(b0; aim = (x, 0.0, 0.0)) for x in xs])
    exact = [min(1 - x - 3cone - tα * sqrt(10), 1 + x - 9cone - tα * sqrt(82)) - R - 0.075 - cone for x in xs]
    @test c ≈ exact rtol = 1e-12
    # the tolerance at the closest point without it, taken after the search: a step
    post = [(x < 3cone ? 1 + x - 9cone - tα * sqrt(82) : 1 - x - 3cone - tα * sqrt(10)) - R - 0.075 - cone for x in xs]
    @test maximum(abs, diff(post)) > 0.04
    @test maximum(abs, diff(c)) <= 1.001e-3          # 1-Lipschitz in the aim
    # the landing: a beam along −x through the near plate, aimed at TCC's plane x = 0
    side = Beam("side"; theta = 90, phi = 0, spot = SuperGaussianSpot(0.2), aim = (0.0, 0.0, 2.0))
    hc(part; kw...) = only(hardware_clearance(HardwareKeepOut(tp, part; kw...), [side]))
    @test hc(:near; targets = (:near,)) == Inf                 # lands on the part itself
    @test hc(:far; targets = (:near,), miss = :all) == Inf      # downstream of the landing
    @test hc(:near; targets = (:far,)) < 0                     # crosses it before the aim plane
    @test hc(:far; targets = ()) == Inf && isfinite(hc(:far; targets = (), miss = :all))

    # against talis_check's :clearance_3w worst values on the mock-up, every row
    t = DR_TARGET
    beams = apply(GroupPointing(nif_quads(); fields = (:z, :r)), DR_X)
    rep = talis_check(t, beams; capture = (:wall, :capsule), hardware = (:shield, :stalk), omega1 = false,
                      resolution = 1e-7, rtol = 0)
    ks = [HardwareKeepOut(t, part; targets = (:wall, :capsule), tolerance = 0.5, miss = :all, minimum = 0.5)
          for part in (:shield, :stalk)]
    rows = [[only(r for r in rep.rows if r.beam == b.name && r.rule === :clearance_3w && r.component === part)
             for b in beams] for part in (:shield, :stalk)]
    c = [hardware_clearance(k, beams) for k in ks]
    w = [[r.worst for r in rs] for rs in rows]
    @test all(all(isfinite.(c[j]) .== isfinite.(w[j])) for j in 1:2)
    @test sum(j -> count(<(0.5), c[j]), 1:2) >= 20 && sum(j -> count(x -> 0.5 < x < Inf, c[j]), 1:2) >= 10
    # talis_check's worst points inside the envelope next to the major axis, where its
    # point–ellipse distance loses digits: there the keep-out is checked against the
    # closed form on the axis, d = −B √(1 − u²/(A² − B²)), at talis_check's point
    medial = [falses(length(beams)) for _ in 1:2]
    onaxis = [fill(Inf, length(beams)) for _ in 1:2]
    for j in 1:2, (i, b) in enumerate(beams)
        isfinite(w[j][i]) || continue
        env = BeamEnvelope(b; peak_fraction = 0.01)
        x = rows[j][i].point
        d = x - env.focus
        u, v, s = dot(d, env.e1), dot(d, env.e2), dot(d, env.k)
        A, B = max(env.a, env.b), min(env.a, env.b)
        u, v = env.a >= env.b ? (u, v) : (v, u)
        if abs(v) < 1e-9 && abs(u) < (A^2 - B^2) / A
            medial[j][i] = true
            onaxis[j][i] = -B * sqrt(1 - u^2 / (A^2 - B^2)) - env.cone * abs(s) - 0.075 - env.cone - tα * norm(x)
        end
    end
    # elsewhere within twice the 1e-7 mm resolution, and never below its lower bound
    rest = [isfinite.(w[j]) .& .!medial[j] for j in 1:2]
    @test sum(sum, rest) >= 40
    @test all(all(0 .<= c[j][rest[j]] .- w[j][rest[j]] .<= 2e-7) for j in 1:2)
    @test sum(sum, medial) >= 2
    @test all(all(c[j][medial[j]] .<= onaxis[j][medial[j]] .+ 1e-9) for j in 1:2)
    @test all(all(0 .<= c[j][medial[j]] .- w[j][medial[j]] .<= 1e-5) for j in 1:2)

    # derivatives in each beam's aim, defocus and spot axes against central differences,
    # the tolerance about TCC and about the beam's own aim point
    err = 0.0
    for pivot in (:tcc, :aim), part in (:shield, :stalk)
        k = HardwareKeepOut(t, part; targets = (:wall, :capsule), tolerance = 0.5, pivot, miss = :all)
        for b in beams
            only(hardware_clearance(k, [b])) < 3 || continue
            mk(y) = remake(b; aim = (y[1], y[2], y[3]), defocus = y[4], spot = with_axes(b.spot, y[5], y[6]))
            y0 = [b.aim..., b.defocus, b.spot.a, b.spot.b]
            f(y) = only(hardware_clearance(k, [mk(y)]))
            g = ForwardDiff.gradient(f, y0)
            fd = [(f(y0 + 1e-6 * (1:6 .== j)) - f(y0 - 1e-6 * (1:6 .== j))) / 2e-6 for j in 1:6]
            err = max(err, maximum(abs.(g .- fd)))
        end
    end
    @test err < 1e-6
    # the rule's penalty
    kp = HardwareKeepOut(t, (:shield, :stalk); targets = (:wall, :capsule), tolerance = 0.5, miss = :all,
                         minimum = 0.5, weight = 2.0)
    cp = hardware_clearance(kp, beams)
    @test all(isapprox.(cp, min.(c[1], c[2]); atol = 1e-12))
    @test keepout_penalty(kp, beams) ≈ 2 * sum(max(0, 0.5 - x)^2 for x in cp if isfinite(x)) rtol = 1e-12
end
