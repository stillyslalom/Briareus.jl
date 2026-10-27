# Named parameter vectors (`parameters`) and the scenario built from a parameterization.

using Briareus, Test

@testset "named parameters" begin
    quads = nif_quads(power = 8.0)
    p = GroupPointing(quads; fields = (:z, :r, :power))
    x = parameters(p; inner23 = (z = 3.09, r = -1.09), outer50 = (z = 5.01,))
    names = parameter_names(p)
    at(n) = x[findfirst(==(n), names)]
    x0 = initial_parameters(p)
    @test (at("inner23.z"), at("inner23.r"), at("outer50.z")) == (3.09, -1.09, 5.01)
    set = Set(["inner23.z", "inner23.r", "outer50.z"])
    @test all(x[i] == x0[i] for i in eachindex(x) if !(names[i] in set))
    @test parameters(p, :inner23 => (z = 3.09, r = -1.09), :outer50 => (z = 5.01,)) == x
    @test parameters(p) == x0
    @test_throws ArgumentError parameters(p; inner99 = (z = 1.0,))
    @test_throws ArgumentError parameters(p; inner23 = (defocus = 1.0,))
    q = GroupPointing(quads; by = :quad, fields = (:power,))
    k = first(q.keys)
    @test parameters(q, k => (power = 0.0,))[1] == 0.0

    # The scenario from a parameterization is the one from its beams.
    t = Target([foil(:foil, (0, 0, 0), (0, 0, 1); width = 4.0, n = 8)])
    a = Scenario(t, p)
    b = Scenario(t, quads; parameterization = p)
    @test a.param === p
    @test evaluate(a, x).deposition.element == evaluate(b, x).deposition.element
end
