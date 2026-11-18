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
    # it reads higher above 1.2 mm on the shield), deep violations included: within its
    # 1e-7 mm resolution. Everywhere the keep-out is the lower: conservative
    near = [findall(x -> x < 1.0, wi) for wi in w]
    @test sum(length, near) >= 12 && any(j -> any(<(-0.75), w[j][near[j]]), 1:2)
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
    # talis_check's worst points inside the envelope next to the major axis (where its
    # point–ellipse distance lost digits before 2026-11-18, reading up to 6e-6 mm low):
    # there the keep-out is also checked against the closed form on the axis,
    # d = −B √(1 − u²/(A² − B²)), at talis_check's point
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
    # every row, those points included, within twice the 1e-7 mm resolution, and never
    # below its lower bound
    rest = [isfinite.(w[j]) .& .!medial[j] for j in 1:2]
    @test sum(sum, rest) >= 40
    @test sum(sum, medial) >= 2
    @test all(all(0 .<= c[j][isfinite.(w[j])] .- w[j][isfinite.(w[j])] .<= 2e-7) for j in 1:2)
    @test all(all(c[j][medial[j]] .<= onaxis[j][medial[j]] .+ 1e-9) for j in 1:2)

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

@testset "design rules: TaLIS 3ω capture as a quantity" begin
    # closed form: a square foil (half side h) at TCC normal to a round beam down −z, the
    # aim moved across its edge x = h. The margin lies at the middle of that edge (the
    # tolerance about TCC there h·tan α); past the edge the continuation less twice the
    # axis's distance from it is the same line
    h = 1.0
    sq = Target([foil(:sq, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0); width = 2h, n = 4)])
    k = CaptureKeepIn(sq, :sq)
    b0 = Beam("b"; theta = 0, phi = 0, spot = SuperGaussianSpot(0.2), aim = (0.0, 0.0, 0.0))
    R = envelope_radius(b0.spot, b0.power; peak_fraction = 0.1) * b0.spot.a
    cone, tα = 1 / (2b0.fnumber), tand(0.5)
    xs = h-0.2:0.01:h+0.2
    c = capture_margin(k, [remake(b0; aim = (x, 0.0, 0.0)) for x in xs])
    @test c ≈ [h - x - R - tα * h - 0.075 - cone for x in xs] rtol = 1e-12
    @test keepout_penalty(k, [b0]) == 0 && keepout_clearance(k, [b0]) == capture_margin(k, [b0])
    # a closed part has no substrate edge
    @test only(capture_margin(CaptureKeepIn(Target([Component(:ball, icosphere_mesh(1.0; n = 4), gold())]), :ball), [b0])) == Inf

    # against talis_check's :capture_3w rows: a 5 × 4 mm one-sided foil at TCC tilted 20°,
    # the 24 upper NIF quads aimed on a grid across and past its edges
    n = Vec3(0.0, sind(20), cosd(20))
    t = Target([foil(:foil, (0.0, 0.0, 0.0), n; width = 5.0, height = 4.0, n = 5, two_sided = false)])
    beams = [remake(b; aim = (-2.6 + 5.2 * ((i - 1) % 6) / 5, -2.1 + 4.2 * ((i - 1) ÷ 6) / 3, 0.0))
             for (i, b) in enumerate(b for b in nif_quads() if b.port[3] > 0)]
    # (the aims on the foil's plane)
    beams = [remake(b; aim = Tuple(b.aim - (dot(b.aim, n) / dot(direction(b), n)) * direction(b))) for b in beams]
    rep = talis_check(t, beams; capture = (:foil,), hardware = (), omega1 = false, resolution = 1e-7, rtol = 0)
    w = [only(r.worst for r in rep.rows if r.beam == b.name && r.rule === :capture_3w) for b in beams]
    kc = CaptureKeepIn(t, :foil; minimum = 0.5, weight = 2.0)
    c = capture_margin(kc, beams)
    landed = isfinite.(w)
    @test count(landed) >= 10 && count(!, landed) >= 8 && count(>(0.5), c) >= 1
    # where the axis lands, within the report's 1e-7 mm resolution (its lower bound)
    @test all(0 .<= c[landed] .- w[landed] .<= 2e-7)
    # where it misses (the report's −Inf), finite and negative
    @test all(isfinite, c) && all(<(0), c[.!landed])
    @test keepout_penalty(kc, beams) ≈ 2 * sum(max(0, 0.5 - x)^2 for x in c) rtol = 1e-12
    # continuity across the edge: a beam's aim moved along the foil's long axis through
    # its end: no step where the axis leaves the foil (3-Lipschitz in the aim), and falling
    # beyond it
    b = beams[findmax(c)[2]]
    e = Vec3(-1.0, 0.0, 0.0)                 # in the foil's plane, along its long side
    ys = -0.3:0.005:0.6
    m = capture_margin(kc, [remake(b; aim = Tuple(b.aim + (2.5 - dot(b.aim, e) + y) * e)) for y in ys])
    @test maximum(abs, diff(m)) <= 3 * 0.005 * 1.0001
    @test all(diff(m[ys .>= 0.05]) .< 0)
    # derivatives in each beam's aim, defocus and spot axes against central differences,
    # landed or not
    err = 0.0
    for b in beams
        mk(y) = remake(b; aim = (y[1], y[2], y[3]), defocus = y[4], spot = with_axes(b.spot, y[5], y[6]))
        y0 = [b.aim..., b.defocus, b.spot.a, b.spot.b]
        f(y) = only(capture_margin(kc, [mk(y)]))
        g = ForwardDiff.gradient(f, y0)
        fd = [(f(y0 + 1e-6 * (1:6 .== j)) - f(y0 - 1e-6 * (1:6 .== j))) / 2e-6 for j in 1:6]
        err = max(err, maximum(abs.(g .- fd)))
    end
    @test err < 1e-6
end

@testset "design rules: TaLIS 3ω incidence as a quantity" begin
    # closed form: a part of two plates, one normal to a round beam down −z (it lands
    # there) and one tilted 60° beside it, its near edge at x = 1.2 on the focal plane.
    # The tilted plate counts in full once the footprint reaches it and fades at
    # `sharpness` degrees per mm of clearance outside it
    pa = foil(:pa, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0); width = 2.0, n = 4)
    pb = foil(:pb, (1.2 + 0.5cosd(60), 0.0, 0.5sind(60)), (-sind(60), 0.0, cosd(60)); width = 1.0, height = 2.0, n = 2)
    tp = Target([Component(c.name, c.mesh, c.material; part = :p) for c in (pa, pb)])
    k = IncidenceRule(tp, :p; max_incidence = 70.0, minimum = 5.0)
    b0 = Beam("b"; theta = 0, phi = 0, spot = SuperGaussianSpot(0.2), aim = (0.0, 0.0, 0.0))
    R = envelope_radius(b0.spot, b0.power; peak_fraction = 0.1) * b0.spot.a
    cone = 1 / (2b0.fnumber)
    hc = atand(cone)
    ds = -0.1:0.01:0.7                        # the tilted plate's clearance from the footprint
    bs = [remake(b0; aim = (1.2 - R - 0.075 - cone - d, 0.0, 0.0)) for d in ds]
    θ = incidence_angle(k, bs)
    @test θ ≈ [max(0.0, 60 - 100 * max(0.0, d)) + hc for d in ds] rtol = 1e-12
    @test keepout_clearance(k, bs) == 70 .- θ
    @test keepout_penalty(k, bs) ≈ 100 * sum(max(0, 5 - (70 - x))^2 for x in θ) rtol = 1e-12
    # talis_check steps where the footprint reaches the plate; the quantity lies at or
    # above it, equal outside the fading band
    w = [only(r.worst for r in talis_check(tp, [b]; capture = (:p,), hardware = (), omega1 = false).rows
              if r.rule === :incidence_3w) for b in bs]
    @test all(θ .>= w .- 1e-9) && all(isapprox.(θ[ds .<= 0], w[ds .<= 0]; rtol = 1e-12)) &&
          all(isapprox.(θ[ds .>= 0.61], w[ds .>= 0.61]; rtol = 1e-12))
    # the derivative in the aim where the plate fades in: 100°/mm
    f(x) = only(incidence_angle(k, [remake(b0; aim = (x[1], 0.0, 0.0))]))
    @test ForwardDiff.gradient(f, [1.2 - R - 0.075 - cone - 0.3])[1] ≈ 100 rtol = 1e-9
    # no landing, no incidence
    @test only(incidence_angle(k, [remake(b0; aim = (5.0, 0.0, 0.0))])) == -Inf

    # against talis_check's :incidence_3w rows on the mock-up hohlraum: the wall and the
    # capsule with their exact normals, every quad (the largest angle on each lies
    # inside the footprint: no fading)
    t = DR_TARGET
    beams = apply(GroupPointing(nif_quads(); fields = (:z, :r)), DR_X)
    rep = talis_check(t, beams; capture = (:wall, :capsule), hardware = (), omega1 = false)
    w = [only(r.worst for r in rep.rows if r.beam == b.name && r.rule === :incidence_3w) for b in beams]
    θ = incidence_angle(IncidenceRule(t, (:wall, :capsule)), beams)
    @test θ ≈ w rtol = 1e-12
    @test count(>(70), θ) >= 4
end

@testset "design rules: TaLIS 1ω clearance as a quantity" begin
    # closed form: one beamlet's footprint (an `offsets` entry: the square displaced along
    # +x) without the taper (a prism) and without the focus error, and a plate across the
    # gap between the beam and the footprint, along the beam at x = 4: the footprint's
    # near edge lies at x = aim + 4.45, so every point of the plate is 0.45 + aim from it,
    # and the tolerance is largest at the plate's corners, √21 mm from TCC
    b0 = Beam("b"; theta = 0, phi = 0, spot = SuperGaussianSpot(0.2), aim = (0.0, 0.0, 0.0))
    pl = Target([foil(:pl, (4.0, 0.0, 1.0), (-1.0, 0.0, 0.0); width = 2.0, n = 2)])
    k = Omega1KeepOut(pl, :pl; offsets = Dict("b" => 0.0), omega1_focus = Inf, lens_distance = Inf, focus = 0.0)
    xs = -0.3:0.05:0.3
    c = omega1_clearance(k, [remake(b0; aim = (x, 0.0, 0.0)) for x in xs])
    @test c ≈ [0.45 + x - 0.075 - tand(0.5) * sqrt(21) for x in xs] rtol = 1e-9
    @test only(omega1_clearance(k, [remake(b0; name = "other")])) == Inf    # no footprint

    # against talis_check's :clearance_1w rows on the mock-up, every NIF quad
    t = DR_TARGET
    beams = apply(GroupPointing(nif_quads(); fields = (:z, :r)), DR_X)
    rep = talis_check(t, beams; capture = (:wall, :capsule), hardware = (), omega1_components = (:shield, :stalk),
                      resolution = 1e-7, rtol = 0)
    for part in (:shield, :stalk)
        w = [only(r.worst for r in rep.rows if r.beam == b.name && r.rule === :clearance_1w && r.component === part)
             for b in beams]
        c = omega1_clearance(Omega1KeepOut(t, part), beams)
        # the report's lower bound, the quantity a found value: within the 1e-7 mm resolution
        @test all(0 .<= c .- w .<= 1.1e-7)
    end
    ko = Omega1KeepOut(t, :stalk; resolution = 1e-11, minimum = 0.2, weight = 3.0)
    c = omega1_clearance(ko, beams)
    @test any(<(0.2), c) && any(>(0.2), c)
    @test keepout_penalty(ko, beams) ≈ 3 * sum(max(0, 0.2 - x)^2 for x in c) rtol = 1e-12
    # continuity: an aim scan, no step beyond the Lipschitz bound
    b = beams[argmin(c)]
    ys = -0.5:0.01:0.5
    m = omega1_clearance(ko, [remake(b; aim = Tuple(b.aim + Vec3(0.0, y, 0.0))) for y in ys])
    @test maximum(abs, diff(m)) <= 0.01 * (ko.lip + tand(0.5)) + 1e-10
    # derivatives in each beam's aim and defocus against central differences
    err = 0.0
    for b in beams
        mk(y) = remake(b; aim = (y[1], y[2], y[3]), defocus = y[4])
        y0 = [b.aim..., b.defocus]
        f(y) = only(omega1_clearance(ko, [mk(y)]))
        f(y0) < 1 || continue
        g = ForwardDiff.gradient(f, y0)
        fd = [(f(y0 + 1e-4 * (1:4 .== j)) - f(y0 - 1e-4 * (1:4 .== j))) / 2e-4 for j in 1:4]
        err = max(err, maximum(abs.(g .- fd)))
    end
    @test err < 1e-8
end
