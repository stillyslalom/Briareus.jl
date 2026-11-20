@testset "picket fluence and the bubble's picket ratio" begin
    # closed forms on a flat plate at normal incidence (collimated beams, focus on the
    # plate): one spot of energy E with best-focus peak fluence F₀ = E/(π a b Γ(1 + 2/n))
    # has the energy-weighted mean ∫F² dA/∫F dA = 2^(-2/n) F₀ (half for a Gaussian). The
    # element average lowers ∫F² by about h²/(6a²) (Gaussian: 1.1e-3 at h = 0.04 mm,
    # a = 0.5 mm), the tolerance's scale.
    plate = Target([foil(:plate, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0); width = 4.0, n = 100)])
    sg(a, b, n) = π * a * b * Briareus.gamma(1 + 2 / n)
    for spot in (SuperGaussianSpot(0.5), SuperGaussianSpot(0.6, 0.4; order = 4.5))
        b = Beam("v"; theta = 0.0, phi = 0.0, spot, fnumber = 1e6, power = 2.0)
        op = DepositionOperator(plate, [b])
        f = picket_fluence(op, [b]; groups = (b.group,), components = (:plate,), duration = 1.5)
        F0 = 3.0 / sg(spot.a, spot.b, spot.order)
        @test f.energy ≈ 3.0 rtol = 1e-4
        @test f.mean ≈ 2^(-2 / spot.order) * F0 rtol = 3e-3
        @test f.mean < 2^(-2 / spot.order) * F0                 # the element average is low
        @test f.peak ≈ F0 rtol = 1e-2
        @test f.area ≈ f.energy / f.mean rtol = 1e-14
        @test all(iszero, f.fluence[setdiff(1:nelements(plate), f.elements)])
        # the pulse form: a constant 2 TW over a 1.5 ns window is the same 3 kJ
        fp = picket_fluence(op, [b], PulseShape([-1.0, 5.0], [2.0, 2.0]), (0.0, 1.5);
                            groups = b.group, components = :plate)
        @test fp.mean ≈ f.mean rtol = 1e-14
    end

    # two Gaussian spots (1/e radius a) a distance d apart: ∫F² gains the cross term
    # 2∫F₁F₂ = 2 e^(-d²/2a²) ∫F₁², so the mean is (1 + e^(-d²/2a²)) times one spot's:
    # twice when they coincide, once when separated
    wide = Target([foil(:plate, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0); width = 6.0, n = 150)])
    g = SuperGaussianSpot(0.5)
    pair = [Beam("a"; theta = 0.0, phi = 0.0, spot = g, fnumber = 1e6, power = 1.0),
            Beam("b"; theta = 0.0, phi = 0.0, spot = g, fnumber = 1e6, power = 1.0)]
    opw = DepositionOperator(wide, pair)
    fl(bs) = picket_fluence(opw, bs; groups = (pair[1].group,), components = (:plate,)).mean
    apart(d) = [remake(pair[1]; aim = (-d / 2, 0.0, 0.0)), remake(pair[2]; aim = (d / 2, 0.0, 0.0))]
    one_spot = fl([pair[1], remake(pair[2]; power = 0.0)])
    @test fl(apart(0.0)) / one_spot ≈ 2 rtol = 1e-12
    @test all(abs(fl(apart(d)) / one_spot - (1 + exp(-d^2 / (2 * 0.5^2)))) < 2e-3 for d in (0.5, 1.0))
    @test fl(apart(3.0)) / one_spot ≈ 1 atol = 1e-4
    @test_throws ArgumentError picket_fluence(opw, pair; groups = (:outer50,), components = (:plate,))
    @test_throws KeyError picket_fluence(opw, pair; groups = (pair[1].group,), components = (:capsule,))

    # a 50° quad on the basecamp's cylinder, split tangentially (Izumi 2021's split-quad
    # direction: the same energy, a lower fluence on the wall): its two pairs of beams
    # separate along the spots' long axes, so the mean falls toward half the unsplit one
    t = hohlraum(length = 10.13, radius = 2.875, leh_radius = 1.685, capsule_radius = 0.908, h = 0.2,
                 capsule_h = 0.6)
    q = first(b.quad for b in nif_beams() if b.group === :outer50)
    quad = [remake(b; power = 0.05) for b in nif_beams() if b.quad == q]
    split_quad(s) = hybride_mock_pointing(quad; table = Dict(:outer50 => (z = 5.01, r = 0.0, s = s)))
    opq = DepositionOperator(t, split_quad(0.0))
    fq(s) = picket_fluence(opq, split_quad(s))
    # (at s = 0.75 mm the pairs are 1.5 mm apart and still clear the entrance hole)
    f0, f5, f75 = fq(0.0), fq(0.5), fq(0.75)
    @test isapprox(f5.energy, f0.energy; rtol = 1e-3) && isapprox(f75.energy, f0.energy; rtol = 1e-3)
    @test f75.mean < f5.mean < f0.mean
    @test f75.mean / f0.mean ≈ 0.5 atol = 0.03
    @test f5.peak < f0.peak

    # the ratio, and the gradient in the pointing (an axial shift and a tangential split
    # of the quad) against central differences, through the ratio into the bubble's speed
    @test picket_ratio(f0, f0) == 1
    @test picket_ratio(f5, f0; at = 1.1) ≈ 1.1 * f5.mean / f0.mean rtol = 1e-15
    eφ = [begin
              φq = atan(sum(x -> x.port[2], quad), sum(x -> x.port[1], quad))
              δ = mod(atan(b.port[2], b.port[1]) - φq + π, 2π) - π
              -sign(δ) * Briareus.Vec3(-sin(φq), cos(φq), 0)
          end for b in quad]
    base = split_quad(0.0)
    moved(x) = [remake(b; aim = b.aim + x[1] * Briareus.Vec3(0, 0, 1) + x[2] * eφ[i]) for (i, b) in enumerate(base)]
    bubble = BubbleAbsorption(radius = 2.875, turn_on = 0.3)
    speed(x) = bubble_speed(remake(bubble; picket_ratio = picket_ratio(picket_fluence(opq, moved(x)), f0)))
    x = [0.1, 0.3]
    gd = ForwardDiff.gradient(speed, x)
    fd = [(speed(x + 1e-4 * e) - speed(x - 1e-4 * e)) / 2e-4 for e in ([1.0, 0.0], [0.0, 1.0])]
    @test gd ≈ fd rtol = 1e-5
    @test gd[2] < 0                                            # splitting slows the bubble

    # the pulse window: exact for a piecewise-linear pulse, zero outside the table
    P = PulseShape([0.0, 1.0, 2.0], [0.0, 2.0, 2.0])
    @test Briareus._window_energy(P, 0.5, 1.5) ≈ 1.75 rtol = 1e-15
    @test Briareus._window_energy(P, -1.0, 3.0) == pulse_energy(P)
    @test Briareus._window_energy(P, 2.5, 3.0) == 0
    @test_throws ArgumentError picket_fluence(opq, quad, P, (1.0, 0.5))

    # the calibration platform's reference: all of the picket energy lands on the wall
    bc = basecamp_picket_fluence(h = 0.4)
    @test bc.energy ≈ 26.0 rtol = 1e-3
    @test picket_ratio(bc, bc) == 1
end

# The picket irradiance as a constraint, the bubble that follows it, and the prediction
# for a re-pointed design. References: `picket_fluence` (the map, checked against closed
# forms above) for the values; central differences for the gradients; bisection on the
# constrained quantity for the solution of a one-parameter run; `remake(bubble;
# picket_ratio)` by hand for the linked bubble; a synthetic history with a known swing.
@testset "picket ceiling and the linked bubble" begin
    geom = (; length = 10.1, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9)
    t = hohlraum(; geom..., h = 0.45)
    aims = Dict(:outer50 => (0.0, 5.01), :outer44 => (0.0, 5.20), :inner30 => (0.60, 3.63),
                :inner23 => (1.09, 3.09))
    quads = map(nif_quads(power = 8.0)) do q
        r, z = aims[q.group]
        φ = port_angles(q)[2]
        remake(q; aim = (-r * cosd(φ), -r * sind(φ), sign(q.port[3]) * z))
    end
    p = GroupPointing(quads; fields = (:z, :phi))
    s = Scenario(t, p)
    x0 = initial_parameters(p)
    ref = evaluate(s, x0)
    # a picket of 2 TW per quad in 0.3–1.3 ns on every cone, then a trough and a rise
    pulse = PulseShape([0.0, 0.3, 1.3, 2.0, 3.0, 8.0], [0.0, 2.0, 2.0, 0.5, 6.0, 6.0])
    pulses = Dict(g => pulse for g in (:outer44, :outer50, :inner23, :inner30))
    w = (0.3, 1.3)
    kw = (; pulses, window = w)

    @testset "values against the picket fluence" begin
        x = x0 .+ 0.05 .* sin.(1:length(x0))
        ev = evaluate(s, x)
        fl = picket_fluence(s.deposition, ev.beams, pulses, w)        # kJ/mm² over 1 ns
        @test picket_irradiance(ev, :wall; kw...) ≈ fl.peak rtol = 1e-12
        @test picket_irradiance(ev, :wall; kw..., measure = :weighted) ≈ fl.mean rtol = 1e-12
        # a picket power per beam, or the beam powers (the outer cones' deposition)
        f3 = picket_fluence(s.deposition, [remake(b; power = 3.0) for b in ev.beams])
        @test picket_irradiance(ev, :wall; power = b -> 3.0, measure = :weighted) ≈ f3.mean rtol = 1e-12
        outer = b -> b.group in (:outer44, :outer50) ? 1.0 : 0.0
        @test picket_irradiance(ev, :wall) == peak_irradiance(ev, :wall; weight = outer)
        # other groups count when named
        fi = picket_fluence(s.deposition, ev.beams, pulses, w; groups = (:inner23, :inner30))
        @test picket_irradiance(ev, :wall; kw..., groups = (:inner23, :inner30)) ≈ fi.peak rtol = 1e-12
        # the ceilings: the stand-in at the reference's peak, the reported value the peak
        c, cm = picket_ceiling(ref, :wall; kw...), picket_ceiling(ref, :wall; kw..., measure = :weighted)
        @test c.name === :picket_peak && cm.name === :picket_weighted && c.upper
        @test only(c.limit) == picket_irradiance(ref, :wall; kw...)
        @test only(cm.limit) == picket_irradiance(ref, :wall; kw..., measure = :weighted)
        @test c.f(ref) ≈ only(c.limit) rtol = 1e-12
        @test check_constraints(ev, c)[1].value == picket_irradiance(ev, :wall; kw...)
        @test check_constraints(ev, picket_ceiling(ref, :wall; kw..., exact = false))[1].value == c.f(ev)
        @test cm.f(ev) == check_constraints(ev, cm)[1].value
        @test cm.f(ev) ≈ fl.mean rtol = 1e-12
        @test only(picket_ceiling(ref, :wall; kw..., limit = 1.5).limit) == 1.5
        # the picket does not see the scenario's laser–plasma model
        sb = remake(s; lpi = Backscatter(:outer50 => 0.5))
        @test c.f(evaluate(sb, x)) == c.f(ev)
        @test_throws ArgumentError picket_irradiance(ev, :wall; pulses)              # no window
        @test_throws ArgumentError picket_irradiance(ev, :wall; window = w)          # no pulses
        @test_throws ArgumentError picket_irradiance(ev, :wall; kw..., power = b -> 1.0)
        @test_throws ArgumentError picket_irradiance(ev, :wall; kw..., measure = :max)
        @test_throws ArgumentError picket_ceiling(ref, :wall; kw..., measure = :max)
        @test_throws ArgumentError picket_ceiling(ref, :no_such_part; kw...)
        other = Scenario(hohlraum(; geom..., h = 0.8), quads)
        @test_throws ArgumentError c.f(evaluate(other, quads))
        @test_throws ArgumentError check_constraints(evaluate(other, quads), c)

        # gradients through the optimizer's duals against central differences
        d = sin.(1:length(x0))
        along(f) = τ -> f(evaluate(s, x .+ τ .* d))
        @test all(f -> begin
                      ad = ForwardDiff.derivative(along(f), 0.0)
                      fd = (along(f)(1e-6) - along(f)(-1e-6)) / 2e-6
                      isapprox(ad, fd; rtol = 1e-5) && ad != 0
                  end, (c.f, cm.f))
    end

    # The published quantity (absorbed outer power on the parts over a fixed area) and the
    # corridor near the inner centerlines, against sums and distances by hand
    @testset "the spatially averaged fluence and the inner-beam corridor" begin
        x = x0 .+ 0.05 .* sin.(1:length(x0))
        ev = evaluate(s, x)
        es = elements(t, :wall)
        outer2 = [remake(b; power = b.group in (:outer44, :outer50) ? 2.0 : 0.0) for b in ev.beams]  # the window's mean
        P = deposit(s.deposition, outer2).element
        @test picket_irradiance(ev, :wall; kw..., measure = :average) ≈ sum(P[es]) / 53 rtol = 1e-12
        @test picket_irradiance(ev, :wall; kw..., measure = :average, area = 10.0) ≈ sum(P[es]) / 10 rtol = 1e-12
        @test picket_irradiance(ev, :wall; kw..., measure = :average) ≈
              picket_fluence(s.deposition, ev.beams, pulses, w).energy / 53 rtol = 1e-12
        # `deposit` skips the beams at exactly zero power (their light is exact zeros): the
        # outer beams alone through the operator on their ports, and the others' totals 0
        live = findall(b -> b.group in (:outer44, :outer50), ev.beams)
        d2 = deposit(s.deposition, outer2)
        @test isapprox(d2.element, deposit(s.deposition[live], outer2[live]).element; rtol = 1e-13) &&
              all(iszero, d2.per_beam[setdiff(eachindex(outer2), live)])
        # net of backscatter: the effective beams at their own powers
        sb = remake(s; lpi = Backscatter(:outer50 => 0.5))
        Pb = deposit(s.deposition, [remake(b; power = b.group === :outer50 ? 4.0 : b.group === :outer44 ? 8.0 : 0.0)
                                    for b in ev.beams]).element
        @test picket_irradiance(evaluate(sb, x), :wall; lpi = true, measure = :average) ≈ sum(Pb[es]) / 53 rtol = 1e-12
        @test picket_irradiance(ev, :wall; lpi = true, measure = :weighted) == picket_irradiance(ev, :wall; measure = :weighted)
        @test_throws ArgumentError picket_irradiance(ev, :wall; kw..., lpi = true)
        # with lpi = true a ceiling on the reference target checks the reference's
        # model and time, as on another target: another model, none against one, or a model
        # at another time is an error; a model built alike is the same model
        refb = evaluate(sb, x0)
        cl = picket_ceiling(refb, :wall; measure = :average, lpi = true)
        evb = evaluate(sb, x)
        @test cl.f(evb) == picket_irradiance(evb, :wall; lpi = true, measure = :average)
        @test cl.f(evaluate(remake(s; lpi = Backscatter(:outer50 => 0.5)), x)) == cl.f(evb)
        @test_throws ArgumentError cl.f(ev)
        @test_throws ArgumentError cl.f(evaluate(remake(s; lpi = Backscatter(:outer50 => 0.4)), x))
        @test_throws ArgumentError cl.f(evaluate(remake(sb; time = 2.0), x))
        @test_throws ArgumentError check_constraints(ev, cl)
        @test_throws ArgumentError picket_ceiling(ref, :wall; measure = :average, lpi = true).f(evb)
        # without lpi the model does not enter, and none is checked
        @test picket_ceiling(refb, :wall; kw..., measure = :average).f(ev) == picket_irradiance(ev, :wall; kw..., measure = :average)
        @test_throws ArgumentError picket_irradiance(ev, :wall; kw..., measure = :average, area = 0.0)
        @test_throws ArgumentError picket_irradiance(ev, :wall; kw..., measure = :corridor, edge = 3.0)
        @test_throws ArgumentError picket_irradiance(ev, :wall; kw..., measure = :corridor, inner = (:none,))
        # the corridor: with no centroid on the soft edge, the hard mean within the
        # distance; the distance to each inner centerline by hand
        lines = [(b.aim, -b.port / norm(b.port)) for b in ev.beams if b.group in (:inner23, :inner30)]
        dist = [minimum(norm(cross(t.elem_centroid[e] - a, k)) for (a, k) in lines) for e in es]
        ds = sort(filter(d -> 1.0 < d < 1.6, dist))      # the widest gap between centroid distances
        g = argmax(diff(ds))
        r, e0 = (ds[g] + ds[g+1]) / 2, (ds[g+1] - ds[g]) / 2
        @test e0 > 1e-3 && !any(d -> r - e0 / 2 <= d <= r + e0 / 2, dist)
        m = dist .< r
        @test picket_irradiance(ev, :wall; kw..., measure = :corridor, corridor = r, edge = e0) ≈
              sum(P[es][m]) / sum(t.elem_area[es][m]) rtol = 1e-12
        # the soft edge: a quintic smoothstep per centerline, joined as a union
        sstep(u) = clamp(u, 0, 1)^3 * (6clamp(u, 0, 1)^2 - 15clamp(u, 0, 1) + 10)
        wt(d1, d2) = 1 - (1 - sstep((1.45 - d1) / 0.3)) * (1 - sstep((1.45 - d2) / 0.3))
        two = [(Briareus.SVector(0.0, 0.0, 0.0), Vec3(0, 0, 1)), (Briareus.SVector(3.0, 0.0, 0.0), Vec3(0, 0, 1))]
        pts = [Vec3(dx, 0.0, 0.7) for dx in (0.2, 1.2, 1.3, 1.4, 1.5, 1.7, 2.0)]
        @test all(p -> isapprox(Briareus._corridor_weight(p, two, 1.15, 1.45, 0.3), wt(p[1], 3 - p[1]); atol = 1e-14), pts)
        # the hottest elements: the peak first, sorted, with their distances
        hot = picket_hot_elements(ev, :wall; kw..., n = 5)
        @test length(hot) == 5 && hot[1].irradiance == picket_irradiance(ev, :wall; kw...)
        @test issorted([h.irradiance for h in hot]; rev = true)
        @test all(h -> isapprox(h.distance, dist[findfirst(==(h.element), es)]; rtol = 1e-12), hot)
        @test all(h -> h.z == t.elem_centroid[h.element][3] && 0 <= h.phi < 360, hot)
        # the ceilings: at the reference's value, the value the measure
        ca = picket_ceiling(ref, :wall; kw..., measure = :average)
        cc = picket_ceiling(ref, :wall; kw..., measure = :corridor)
        @test ca.name === :picket_average && cc.name === :picket_corridor
        @test only(ca.limit) == picket_irradiance(ref, :wall; kw..., measure = :average)
        @test only(cc.limit) == picket_irradiance(ref, :wall; kw..., measure = :corridor)
        @test ca.f(ev) == check_constraints(ev, ca)[1].value == picket_irradiance(ev, :wall; kw..., measure = :average)
        @test cc.f(ev) == check_constraints(ev, cc)[1].value == picket_irradiance(ev, :wall; kw..., measure = :corridor)
        # gradients against central differences (the 50° quads up toward the holes, where
        # the light starts to clip, for the average; the corridor moves with the inner aims)
        xc = copy(x)
        xc[findall(n -> n == "outer50.z", parameter_names(p))] .= 6.0
        d = sin.(1:length(x0))
        at(f, x1) = τ -> f(evaluate(s, x1 .+ τ .* d))
        @test all(((f, x1),) -> begin
                      ad = ForwardDiff.derivative(at(f, x1), 0.0)
                      fd = (at(f, x1)(1e-6) - at(f, x1)(-1e-6)) / 2e-6
                      isapprox(ad, fd; rtol = 1e-5) && ad != 0
                  end, ((ca.f, xc), (cc.f, x)))
    end

    # Inside a constrained run the evaluations keep the element powers of the group sets
    # that ceilings at the beam powers read (`_with_split`, per phase on a phased
    # scenario), and the ceilings read them instead of depositing: the sums of their own
    # deposition (bitwise where the compiler builds both loops alike; held here to
    # 1e-12 in values and partials), and the deposition unchanged. Other powers, and
    # setup beams that the laser–plasma model changed, deposit as before.
    @testset "the deposition in the evaluation serves the ceilings" begin
        x = x0 .+ 0.05 .* sin.(1:length(x0))
        xd = [ForwardDiff.Dual{:split}(x[i], (sin(i), cos(i))) for i in eachindex(x)]
        sb = remake(s; lpi = Backscatter(:outer50 => 0.5))
        cs = (picket_ceiling(ref, :wall), picket_ceiling(ref, :wall; measure = :weighted),
              picket_ceiling(ref, :wall; measure = :average), picket_ceiling(ref, :wall; measure = :corridor),
              picket_ceiling(ref, :wall; groups = (:outer50,), measure = :average),   # a set within another
              picket_ceiling(ref, :wall; kw..., measure = :average),                  # picket powers: not kept
              picket_ceiling(evaluate(sb, x0), :wall; measure = :average, lpi = true))
        req = Briareus._split_request(cs)
        @test sort(unique(last.(req))) == [[:outer44, :outer50], [:outer50]] && all(isnothing, first.(req))
        vp(a) = a isa ForwardDiff.Dual ? [ForwardDiff.value(a); collect(ForwardDiff.partials(a))] : [a]
        same(a, b) = isapprox(vp(a), vp(b); rtol = 1e-12)
        flat(v) = reduce(vcat, vp.(v))
        for (sc, n) in ((s, 6), (sb, 7)), xx in (x, xd)
            ev, evs = evaluate(sc, xx), evaluate(Briareus._with_split(sc, req), xx)
            @test isapprox(flat(ev.deposition.element), flat(evs.deposition.element); rtol = 1e-12)
            @test all(c -> same(c.f(ev), c.f(evs)) && same(something(Briareus._exact_value(c.f, ev), 0.0),
                                                            something(Briareus._exact_value(c.f, evs), 0.0)), cs[1:n])
            # read from the evaluation's deposition, except the setup beams that the
            # backscatter changed
            kept(lpi) = Briareus._split_map(evs, Briareus._OwnPower(), (:outer50, :outer44), lpi) !== nothing
            @test kept(true) && kept(false) == (sc === s)
        end
        # a phased scenario: the phase whose constraint reads the groups keeps them
        ps = PhasedScenario(p, Phase(:a, s), Phase(:b, s; power = 4.0, scaled = false))
        cb = in_phase(:b, cs[3])
        pss = Briareus._with_split(ps, Briareus._split_request((cb, in_phase(:a, cs[6]))))
        @test isempty(getfield(pss[:a], :_split)) && getfield(pss[:b], :_split) == [[:outer44, :outer50]]
        xp = initial_parameters(ps.param) .+ 0.01
        @test cb.f(evaluate(pss, xp)) ≈ cb.f(evaluate(ps, xp)) rtol = 1e-12
        @test Briareus._with_split(ps, Briareus._split_request((in_phase(:a, cs[6]),))) === ps
    end

    # One free parameter, the 50° cone's z, pulled toward the waist by the objective: the
    # picket ceiling stops it where the constrained quantity reaches the limit (a root
    # found by bisection on the quantity itself; each is monotone over the bracket)
    @testset "constrained runs end on the picket ceiling" begin
        h = HeldParameters(p, n -> n != "outer50.z")
        sh = remake(s; parameterization = h)
        pull(ev) = sum(b.aim[3]^2 for b in ev.beams if b.group === :outer50) / 25
        at(z) = evaluate(sh, [z])
        function root(g, a, b)              # bisection for g = 0 with g(a) > 0 > g(b)
            @assert g(a) > 0 > g(b)
            for _ in 1:40
                m = (a + b) / 2
                g(m) > 0 ? (a = m) : (b = m)
            end
            return (a + b) / 2
        end
        # the energy-weighted mean rises as the 50° spots move onto the 44° ones (over
        # 4.75–5.1 mm; the box keeps the run there)
        mean_at(z) = picket_irradiance(at(z), :wall; kw..., measure = :weighted)
        cm = picket_ceiling(ref, :wall; kw..., measure = :weighted, limit = mean_at(4.9))
        zm = root(z -> mean_at(z) - only(cm.limit), 4.8, 5.01)
        rm = optimize_pointing(sh, pull; x0 = [5.01], lower = [4.76], upper = [5.1], constraints = (cm,),
                               feasibility = 1e-6)
        @test rm.outcome === :interior && abs(only(rm.x) - zm) < 1e-5
        @test only(rm.constraints).status === :active && only(rm.constraints).multiplier > 0
        # the peak: its hottest element rises over 4.8–4.6 mm as the 50° spots cross the
        # 44° ones. Its stand-in (exact = false) ends active at its root; the exact
        # peak ends at or below the limit (here 0.4% below: the stand-in runs above the
        # peak near this limit, and the working limit is only lowered), no lower than
        # the peak's root
        peak_at(z) = picket_irradiance(at(z), :wall; kw...)
        lim = (peak_at(4.65) + peak_at(4.75)) / 2
        cs = picket_ceiling(ref, :wall; kw..., limit = lim, exact = false)
        zs = root(z -> cs.f(at(z)) - lim, 4.6, 4.8)
        rs = optimize_pointing(sh, pull; x0 = [4.85], lower = [4.6], upper = [5.1], constraints = (cs,),
                               feasibility = 1e-6)
        @test rs.outcome === :interior && abs(only(rs.x) - zs) < 1e-5
        @test only(rs.constraints).status === :active && only(rs.constraints).multiplier > 0
        cp = picket_ceiling(ref, :wall; kw..., limit = lim)
        zp = root(z -> peak_at(z) - lim, 4.6, 4.8)
        rp = optimize_pointing(sh, pull; x0 = [4.85], lower = [4.6], upper = [5.1], constraints = (cp,),
                               feasibility = 1e-6)
        row = only(rp.constraints)
        @test rp.outcome === :interior && row.status !== :violated && row.multiplier > 0
        @test row.value == peak_at(only(rp.x)) && 0.99 * lim <= row.value <= lim * (1 + 1e-6)
        @test zp - 1e-6 <= only(rp.x) < zp + 0.01
        # the spatially averaged fluence falls as the 50° spots clip on the holes (above
        # ~5.6 mm); the corridor's mean falls as they move away from the inner beams' paths
        # (and the average at the beam powers, which the run's evaluations serve from
        # their own deposition)
        for (meas, zl, zbox, z0, pk) in ((:average, 6.0, (5.7, 6.4), 6.3, kw), (:corridor, 5.2, (4.9, 6.0), 5.8, kw),
                                         (:average, 6.0, (5.7, 6.4), 6.3, (;)))
            val(z) = picket_irradiance(at(z), :wall; pk..., measure = meas)
            cv = picket_ceiling(ref, :wall; pk..., measure = meas, limit = val(zl))
            zr = root(z -> val(z) - only(cv.limit), zl - 0.2, zl + 0.2)
            rv = optimize_pointing(sh, pull; x0 = [z0], lower = [zbox[1]], upper = [zbox[2]], constraints = (cv,),
                                   feasibility = 1e-6)
            @test rv.outcome === :interior && abs(only(rv.x) - zr) < 1e-5
            @test only(rv.constraints).status === :active && only(rv.constraints).multiplier > 0
        end
    end

    @testset "the linked bubble" begin
        bubble = BubbleAbsorption(radius = 2.875, turn_on = 0.3, picket_ratio = 1.1)
        linked = PicketLinkedBubble(bubble, t, quads; kw...)
        @test length(linked.beams) == count(b -> b.group in (:outer44, :outer50), quads)
        @test occursin("PicketLinkedBubble", sprint(show, linked))
        # at the reference: bitwise the bubble itself
        lb = linked_bubble(linked, quads)
        @test lb.picket_ratio === bubble.picket_ratio && lb.collision === bubble.collision
        for τ in (5.0, 7.0)
            @test [b.power for b in apply_lpi(linked, quads, τ)] == [b.power for b in apply_lpi(bubble, quads, τ)]
            eff = apply_lpi(bubble, quads, τ)
            @test lpi_deposit(linked, s.deposition, eff, τ).element == lpi_deposit(bubble, s.deposition, eff, τ).element
        end
        # a re-pointed design (the outer cones down 0.3 mm toward the waist): E from its
        # picket fluence against the reference's, as `remake` by hand
        design = [b.group in (:outer44, :outer50) ? remake(b; aim = b.aim - sign(b.port[3]) * 0.3 * Vec3(0, 0, 1)) : b
                  for b in quads]
        hand(bs) = picket_fluence(DepositionOperator(t, bs), bs, pulses, w)
        E = picket_ratio(hand(design), hand(quads); at = 1.1)
        @test abs(E - 1.1) > 1e-3
        ld = linked_bubble(linked, design)
        @test ld.picket_ratio ≈ E rtol = 1e-12
        byhand = remake(bubble; picket_ratio = E)
        @test ld.collision ≈ byhand.collision rtol = 1e-12
        @test bubble_speed(ld) ≈ bubble_speed(byhand) rtol = 1e-12
        st = remake(s; lpi = linked, time = 7.0)
        @test evaluate(st, design).deposition.element ≈
              evaluate(remake(s; lpi = byhand, time = 7.0), design).deposition.element rtol = 1e-12
        @test bubble_fractions(linked, design, 6.0).absorbed == bubble_fractions(ld, design, 6.0).absorbed
        # in a chain after another model; a dropped bundle deposits nothing itself
        @test evaluate(remake(st; lpi = (Backscatter(:outer44 => 0.0), linked)), design).deposition.element ==
              evaluate(st, design).deposition.element
        bundle = nif_bundles(quads)["11"]
        keep = setdiff(eachindex(quads), bundle)
        opq = DepositionOperator(t, quads)
        at2(bs) = [remake(b; power = 2.0) for b in bs]     # the window's mean power
        Ek = picket_ratio(picket_fluence(opq, [i in bundle ? remake(b; power = 0.0) : b for (i, b) in enumerate(at2(quads))]),
                          picket_fluence(opq, at2(quads)); at = 1.1)
        @test any(b -> b.group in (:outer44, :outer50), quads[bundle]) && Ek < 1.1
        @test linked_bubble(linked, quads[keep]).picket_ratio ≈ Ek rtol = 1e-12
        # differentiable through the pointing, the bubble included (at the collision,
        # where the ridge stop rises with time)
        d = sin.(1:length(x0))
        f(τ, m) = sum(evaluate(remake(s; lpi = m, time = bubble.collision), x0 .+ τ .* d).deposition.element .*
                      (1:nelements(t)) ./ nelements(t))
        ad = ForwardDiff.derivative(τ -> f(τ, linked), 0.0)
        fd = (f(1e-6, linked) - f(-1e-6, linked)) / 2e-6
        @test ad ≈ fd rtol = 1e-5
        @test abs(ad - ForwardDiff.derivative(τ -> f(τ, bubble), 0.0)) > 1e-6 * abs(ad)   # the link counts
        # errors: a beam the map does not know, no picket, the per-beam sensitivities
        @test_throws ArgumentError linked_bubble(linked, [remake(first(b for b in quads if b.group === :outer50); name = "X")])
        @test_throws ArgumentError PicketLinkedBubble(bubble, t, quads)
        @test_throws ArgumentError PicketLinkedBubble(bubble, t, quads; kw..., groups = (:none,))
        @test_throws ArgumentError Briareus._ridge_split(linked, s.deposition, quads[1], 7.0)

        # in a drive history: bitwise the resolved bubble's (E is the same number)
        t8 = hohlraum(; geom..., h = 0.8, capsule_h = 0.4, wall = Material(:Au; albedo = 0.9, conversion = 0.8))
        model, op8 = RadiationModel(t8), DepositionOperator(t8, quads)
        times = [1.0, 6.0, 7.0]
        run(bs, m) = drive_history(model, op8, τ -> with_pulses(bs, pulses, τ), times; lpi = m)
        hl = run(design, linked)
        @test hl.moments == run(design, ld).moments
        # the prediction for the design against the reference: the picket and the bubble
        pr = bubble_prediction(linked, design, quads)
        fd_, fr_ = hand(design), hand(quads)
        @test pr.fluence_ratio ≈ fd_.mean / fr_.mean rtol = 1e-12
        @test pr.peak_ratio ≈ fd_.peak / fr_.peak rtol = 1e-12
        @test pr.picket_ratio == ld.picket_ratio && pr.reference_picket_ratio === 1.1
        @test pr.speed_ratio ≈ pr.fluence_ratio^0.2 rtol = 1e-12
        @test pr.collision == ld.collision && pr.reference_collision === bubble.collision && pr.history === nothing
        # with the history: the design's P2/P0 is the linked model's history
        pp = bubble_prediction(linked, design, quads; history = run, fixed = true, swing = (0.5, 8.0))
        @test pp.history.p2 == hl.moments[3, :] ./ hl.moments[1, :]
        @test pp.history.fixed_p2 == (r -> r.moments[3, :] ./ r.moments[1, :])(run(design, bubble))
        @test pp.history.through_bubble.minimum + pp.history.through_deposition.minimum ≈ pp.history.shift.minimum
        # the swing metrics against a synthetic history whose minimum (−0.2 at 5 ns + δ)
        # and upward zero (+0.5 ns after it) move with the bubble's speed
        ts = collect(4.0:0.05:8.0)
        κ = 0.05 / (bubble_speed(ld) - bubble_speed(bubble))       # a 0.05 ns shift for the design
        synth(bs, m) = (time = ts, moments = [ones(length(ts))'; zeros(length(ts))';
                                              [0.8 * (τ - 5 - κ * (bubble_speed(m) - bubble_speed(bubble)))^2 - 0.2
                                               for τ in ts]'])
        sp = bubble_prediction(linked, design, quads; history = synth)
        @test sp.history.swing.reference.minimum ≈ -0.2 && sp.history.swing.reference.at ≈ 5.0
        @test sp.history.swing.reference.crossing ≈ 5.5 atol = 1e-3
        @test sp.history.shift.crossing ≈ 0.05 atol = 2e-3
        @test sp.history.fixed_p2 === nothing && sp.history.through_bubble === nothing
        @test bubble_prediction(linked, design, quads; history = synth,
                                reference_history = synth(quads, bubble)).history.shift == sp.history.shift
    end
end
