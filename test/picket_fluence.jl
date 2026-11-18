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
        @test picket_irradiance(ev, :wall; kw..., measure = :mean) ≈ fl.mean rtol = 1e-12
        # a picket power per beam, or the beams' own powers (the outer cones' deposition)
        f3 = picket_fluence(s.deposition, [remake(b; power = 3.0) for b in ev.beams])
        @test picket_irradiance(ev, :wall; power = b -> 3.0, measure = :mean) ≈ f3.mean rtol = 1e-12
        outer = b -> b.group in (:outer44, :outer50) ? 1.0 : 0.0
        @test picket_irradiance(ev, :wall) == peak_irradiance(ev, :wall; weight = outer)
        # other groups count when named
        fi = picket_fluence(s.deposition, ev.beams, pulses, w; groups = (:inner23, :inner30))
        @test picket_irradiance(ev, :wall; kw..., groups = (:inner23, :inner30)) ≈ fi.peak rtol = 1e-12
        # the ceilings: the stand-in at the reference's peak, the reported value the peak
        c, cm = picket_ceiling(ref, :wall; kw...), picket_ceiling(ref, :wall; kw..., measure = :mean)
        @test c.name === :picket_peak && cm.name === :picket_mean && c.upper
        @test only(c.limit) == picket_irradiance(ref, :wall; kw...)
        @test only(cm.limit) == picket_irradiance(ref, :wall; kw..., measure = :mean)
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
        mean_at(z) = picket_irradiance(at(z), :wall; kw..., measure = :mean)
        cm = picket_ceiling(ref, :wall; kw..., measure = :mean, limit = mean_at(4.9))
        zm = root(z -> mean_at(z) - only(cm.limit), 4.8, 5.01)
        rm = optimize_pointing(sh, pull; x0 = [5.01], lower = [4.76], upper = [5.1], constraints = (cm,),
                               feasibility = 1e-6)
        @test rm.outcome === :interior && abs(only(rm.x) - zm) < 1e-5
        @test only(rm.constraints).status === :active && only(rm.constraints).multiplier > 0
        # the peak: its hottest element rises over 4.8–4.6 mm as the 50° spots cross the
        # 44° ones. Its stand-in (exact = false) ends active at its own root; the exact
        # peak ends at or below the limit (here 0.4 % below: the stand-in runs above the
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
    end

    @testset "the linked bubble" begin
        bubble = BubbleAbsorption(radius = 2.875, turn_on = 0.3, picket_ratio = 1.1)
        linked = PicketLinkedBubble(bubble, t, quads; kw...)
        @test length(linked.beams) == count(b -> b.group in (:outer44, :outer50), quads)
        @test occursin("PicketLinkedBubble", sprint(show, linked))
        # at the reference: the bubble itself, bit for bit
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
        # in a chain after another model; a dropped bundle deposits nothing of its own
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

        # in a drive history: the resolved bubble's, bit for bit (E is the same number)
        t8 = hohlraum(; geom..., h = 0.8, capsule_h = 0.4, wall = Material(:Au; albedo = 0.9, conversion = 0.8))
        model, op8 = RadiationModel(t8), DepositionOperator(t8, quads)
        times = [1.0, 6.0, 7.0]
        run(bs, m) = drive_history(model, op8, τ -> with_pulses(bs, pulses, τ), times; lpi = m)
        hl = run(design, linked)
        @test hl.moments == run(design, ld).moments
    end
end
