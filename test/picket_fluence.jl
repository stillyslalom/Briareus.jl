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
