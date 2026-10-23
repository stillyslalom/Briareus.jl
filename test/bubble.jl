@testset "inner-beam bubble absorption" begin
    pointing = Dict(:outer50 => (0.0, 5.01), :outer44 => (0.0, 5.20),
                    :inner30 => (0.60, 3.63), :inner23 => (1.09, 3.09))
    function aimed(b)
        r, z = pointing[b.group]
        φ = port_angles(b)[2]
        return remake(b; aim = (-r * cosd(φ), -r * sind(φ), sign(b.port[3]) * z))
    end
    quads = aimed.(nif_quads())
    inner(bs) = [b.group in (:inner23, :inner30) for b in bs]

    # Izumi 2021's transparency anchors: 29 % / 54 % at 5.5 ns (0.3 / 0.15 mg/cc) in a
    # 2.88 mm hohlraum, launched 0.2 ns after a 0.3 ns turn-on
    m3 = BubbleAbsorption(radius = 2.88, turn_on = 0.3)
    m15 = BubbleAbsorption(radius = 2.88, turn_on = 0.3, fill = 0.15)
    @test bubble_transmission(m3, 5.5) ≈ 0.29 rtol = 1e-10
    @test bubble_transmission(m15, 5.5) ≈ 0.54 rtol = 1e-10
    @test bubble_depth(m3, 5.5) ≈ 1.75
    @test bubble_speed(m15) / bubble_speed(m3) ≈ 2^0.2
    @test bubble_speed(BubbleAbsorption(radius = 2.88, energy = 1.5)) ≈ 0.35 * 1.5^0.2
    @test 1.7 < m3.density_exponent < 2.0                     # near IB's ρ²
    @test bubble_transmission(m3, 0.4) == 1 && bubble_depth(m3, 0.4) == 0
    @test bubble_transmission(m3, 20.0) == 0

    # the clipped fraction against a closed form: a Gaussian spot (1/e radius a) on a
    # vertical, nearly collimated beam keeps exp(-ρ²/a²) of its power outside radius ρ
    b = Beam("v"; theta = 0.0, phi = 0.0, spot = SuperGaussianSpot(0.5), fnumber = 1e6)
    for ρ in (0.2, 0.5, 0.8)
        @test iris_fraction(b, -1.0, ρ; edge = 0.01, nr = 400, nangle = 8, ncone = 1) ≈
              exp(-(ρ / 0.5)^2) atol = 3e-3
    end

    m = BubbleAbsorption(radius = 2.875)
    # before the launch nothing moves
    @test [x.power for x in apply_lpi(m, quads, 0.1)] == [x.power for x in quads]
    # power is conserved, in each hemisphere, and moves from inner to 50° beams only
    for t in (3.0, 5.0, 7.0)
        eff = apply_lpi(m, quads, t)
        for h in (1, -1)
            side(bs) = sum(x.power for x in bs if (x.port[3] > 0) == (h > 0))
            @test side(eff) ≈ side(quads) rtol = 1e-13
        end
        @test all(x.power < y.power for (x, y) in zip(eff, quads) if x.group in (:inner23, :inner30))
        @test all(x.power > y.power for (x, y) in zip(eff, quads) if x.group === :outer50)
        @test all(x.power == y.power for (x, y) in zip(eff, quads) if x.group === :outer44)
    end
    # the absorbed fraction grows with time, from zero to all of the inner power
    ts = 0.0:0.25:9.0
    A = [sum(bubble_fractions(m, quads, t).absorbed) / count(inner(quads)) for t in ts]
    @test A[1] == 0
    @test all(diff(A) .>= 0)
    @test A[end] > 0.99
    fr = bubble_fractions(m, quads, 5.0)
    @test fr.footprint_z[1] ≈ 5.01 - 2.875 / tand(50) rtol = 1e-3    # 50° axes reach the wall
    @test fr.footprint_z[2] ≈ -fr.footprint_z[1]
    @test all(fr.absorbed .<= fr.clipped)
    # a moving wall moves the footprint plane toward the entrance hole
    mw = BubbleAbsorption(radius = 2.875, wall = WallMotion(speed = 0.1))
    @test bubble_fractions(mw, quads, 5.0).footprint_z[1] > fr.footprint_z[1]

    # differentiable in the calibration (speed, opacity) and in the pointing
    function moved(v, τ, dz)
        bs = [x.group === :inner23 ? remake(x; aim = x.aim .+ (0, 0, sign(x.port[3]) * dz)) : x for x in quads]
        eff = apply_lpi(remake(m; speed = v, opacity = τ), bs, 5.0)
        return sum(x.power for x in eff if x.group === :outer50)
    end
    p = [0.35, 0.144, 0.0]
    g = ForwardDiff.gradient(q -> moved(q...), p)
    for i in 1:3
        δ = zeros(3)
        δ[i] = 1e-5
        @test g[i] ≈ (moved((p .+ δ)...) - moved((p .- δ)...)) / 2e-5 rtol = 1e-5
    end
    @test g[1] > 0 && g[2] > 0 && g[3] != 0

    # in a drive history the moved power raises the pole drive late in the pulse
    t = hohlraum(length = 10.1, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9,
                 h = 0.8, capsule_h = 0.4, wall = Material(:Au; albedo = 0.9, conversion = 0.8))
    model, op = RadiationModel(t), DepositionOperator(t, quads)
    times = [1.0, 6.0]
    h0 = drive_history(model, op, τ -> quads, times)
    h1 = drive_history(model, op, τ -> quads, times; lpi = m)
    p2(h) = h.moments[3, :] ./ h.moments[1, :]
    @test p2(h1)[1] ≈ p2(h0)[1] atol = 1e-6
    @test p2(h1)[2] > p2(h0)[2] + 0.02
end
