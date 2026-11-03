# Laser–plasma models at their limits: zero group power, fractions bounded by [0, 1].

@testset "laser-plasma models at zero power and bounded fractions" begin
    group(pin, pout; dl = 0.0) =
        [[Beam("i$k"; theta = 23.5, phi = 90k, power = pin, group = :inner, dlambda = dl) for k in 1:4];
         [Beam("o$k"; theta = 50, phi = 90k + 45, power = pout, group = :outer, dlambda = k / 4) for k in 1:4]]
    ct = CalibratedTransfer((:outer,), (:inner,), linear_transfer(0.1, 0.02))
    # donors without power give nothing, and the derivatives with respect to their power
    # are the transfer's: receivers gain f × (donor power) shared by their power
    bs = apply_lpi(ct, group(2.0, 0.0), 0.0)
    @test all(b -> isfinite(b.power), bs) && sum(b.power for b in bs) == 8.0
    @test all(b -> b.power == 2.0, bs[1:4]) && all(b -> b.power == 0.0, bs[5:8])
    # f at Δλ = 0 − mean(1/4 … 1) = −0.625 Å: 0.1 − 0.02 × 0.625 = 0.0875
    f = 0.1 - 0.02 * 0.625
    gain(x) = sum(b.power for b in apply_lpi(ct, group(2.0, x), 0.0) if b.group === :inner)
    @test ForwardDiff.derivative(gain, 0.0) ≈ 4f rtol = 1e-12
    keep(x) = sum(b.power for b in apply_lpi(ct, group(2.0, x), 0.0) if b.group === :outer)
    @test ForwardDiff.derivative(keep, 0.0) ≈ 4 * (1 - f) rtol = 1e-12
    # receivers without power take nothing: the beams come back unchanged
    bs = apply_lpi(ct, group(0.0, 3.0), 0.0)
    @test [b.power for b in bs] == [b.power for b in group(0.0, 3.0)]
    # the summary of a group without power: ratio 1 (no change)
    s = lpi_summary(group(0.0, 3.0), bs)
    @test s[:inner].ratio == 1 && s[:outer].ratio == 1
    @test s[:inner].setup == 0 && s[:inner].effective == 0

    # backscattered fractions lie in [0, 1]
    @test_throws ArgumentError Backscatter(:inner => 1.2)
    @test_throws ArgumentError Backscatter(:inner => -0.1)
    @test_throws ArgumentError Backscatter(:inner => ([0.0, 5.0], [0.0, 1.5]))
    @test_throws ArgumentError apply_lpi(Backscatter(:inner => t -> 0.2t), group(2.0, 3.0), 6.0)
    @test apply_lpi(Backscatter(:inner => t -> 0.2t), group(2.0, 3.0), 2.0)[1].power ≈ 1.2
    @test all(b -> b.power == 0, apply_lpi(Backscatter(:inner => 1.0), group(2.0, 3.0))[1:4])
    @test Backscatter(Dict(:inner => 0.1)).fractions[:inner] == 0.1
end
