# Power per beam and component, from the deposition operator and from the reference ray
# trace, against the closed form of a round Gaussian spot split by the edge between two
# coplanar foils: the fraction past the edge at distance x₀ from the axis is ½ erfc(−x₀/a).

using SpecialFunctions: erfc

@testset "beam and ray-trace component power" begin
    # two 4 × 4 mm foils meeting along x = 0 in the plane z = 0
    t = Target(foil(:left, (-2, 0, 0), (0, 0, 1); width = 4.0, n = 40),
               foil(:right, (2, 0, 0), (0, 0, 1); width = 4.0, n = 40))
    a = 0.5
    sp = SuperGaussianSpot(a; order = 2.0)              # exp(-ρ²): a 1D Gaussian of variance a²/2
    beam(x0; P = 2.0, theta = 0) = Beam("b$theta"; theta, phi = 0, power = P, spot = sp, fnumber = 8.0,
                                        aim = (x0, 0.0, 0.0))
    right(x0) = erfc(-x0 / a) / 2
    x0 = [0.0, 0.2, -0.35]
    # the ray trace: every beam's power on the foils, split by the closed form to its
    # sampling (a column of the 121 × 121 spot grid straddling the edge lands whole on one
    # foil: ~h/(2√π) = 0.025 of the power), and mirrored exactly with the aim
    R = raytrace_component_power(t, beam.(x0); nspot = 121)
    @test size(R) == (3, 2)
    @test sum(R; dims = 2) ≈ fill(2.0, 3, 1) rtol = 1e-7
    @test R[:, 2] ./ 2 ≈ right.(x0) atol = 0.02
    Rm = raytrace_component_power(t, beam.([0.2, -0.2]); nspot = 121)
    @test Rm[1, 2] ≈ Rm[2, 1] rtol = 1e-12
    # the operator (built once for the port, serving every aim), to its element quadrature
    op = DepositionOperator(t, [beam(0.0)])
    M = beam_component_power(op, beam.(x0))
    @test sum(M; dims = 2) ≈ fill(2.0, 3, 1) rtol = 2e-3
    @test M[:, 2] ./ 2 ≈ right.(x0) rtol = 5e-3
    # the rows are the beams' own: linear in power, and the per-beam totals of `deposit`
    @test beam_component_power(op, [beam(0.2; P = 6.0)]) ≈ 3 .* M[2:2, :] rtol = 1e-12
    d = deposit(op, [beam(0.2)])
    @test sum(M[2, :]) ≈ sum(d.element) + sum(d.backside) rtol = 1e-12
    # a beam from a port the operator does not hold
    @test_throws ArgumentError beam_component_power(op, [beam(0.0; theta = 10)])
end
