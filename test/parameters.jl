# Named parameter vectors (`parameters`), the scenario built from a parameterization,
# the drive level (`mean_flux`), optimization results and `report`.

using Briareus, Test

# A parameterization other than GroupPointing: shifts of the parameters of `p` from `x0`.
struct _Shift
    p::GroupPointing{Symbol}
    x0::Vector{Float64}
end
Briareus.apply(s::_Shift, x::AbstractVector) = apply(s.p, s.x0 .+ x)
Briareus.initial_parameters(s::_Shift) = zeros(length(s.x0))

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
    b = Scenario(t, quads)
    @test a.param === p
    @test evaluate(a, x).deposition.element == evaluate(b, apply(p, x)).deposition.element
    # any other parameterization: the operator for its beams at the initial parameters
    c = Scenario(t, _Shift(p, x))
    @test evaluate(c, zeros(nparams(p))).deposition.element ==
          evaluate(Scenario(t, apply(p, x)), apply(p, x)).deposition.element
end

@testset "drive level, optimization result, report" begin
    t = hohlraum(nz = 8, nphi = 16, capsule_ntheta = 6, capsule_nphi = 8)
    p = GroupPointing(nif_quads(power = 8.0); fields = (:z, :r, :power))
    s = Scenario(t, p; observe = :capsule)
    x0 = parameters(p; inner23 = (z = 3.09, r = -1.09), inner30 = (z = 3.63, r = -0.60),
                    outer44 = (z = 5.20,), outer50 = (z = 5.01,))
    ev = evaluate(s, x0)
    @test mean_flux(ev) == flux_moments(ev)[1]

    # report: one line per part above `minor` of the delivered power, the rest summed
    lpi = CalibratedTransfer((:outer44, :outer50), (:inner23, :inner30), linear_transfer(0.0, 0.05))
    out = sprint(report, evaluate(Scenario(t, p; observe = :capsule, lpi), x0))
    lines = split(out, '\n')
    @test startswith(lines[1], "Setup laser power: 384 TW; after the laser–plasma model:")
    @test startswith(lines[2], "  inner23 ")
    @test count(l -> startswith(l, "  wall "), lines) == 1
    @test count(l -> occursin("parts < 0.1%", l), lines) == 1
    @test !occursin("cap_top_outside", out) && !occursin("imbalance", out)
    full = sprint((io, e) -> report(io, e; minor = 0), ev)
    @test occursin("cap_top_outside", full)
    @test occursin("missed", full) || occursin("quadrature excess", full)

    # the result names the parameters a narrow box holds at a bound
    lo, hi = bounds(p, x0; z = 1e-3, r = 1e-3, power = 1e-3)
    res = optimize_pointing(s, symmetry_objective(; drive = mean_flux(ev)); x0, lower = lo,
                            upper = hi, iterations = 5)
    @test res isa PointingResult && res.names == parameter_names(p)
    side = [x <= l ? -1 : x >= h ? 1 : 0 for (x, l, h) in zip(res.x, lo, hi)]
    @test !isempty(res.at_bounds) && res.at_bounds == res.names[side .!= 0]
    shown = sprint(show, MIME"text/plain"(), res)
    @test occursin("at a bound", shown) && all(n -> occursin(n, shown), res.at_bounds)
    @test occursin("parameters at bounds", sprint(show, res))
    # infinite and collapsed boxes, and names of a parameterization without parameter_names
    @test Briareus._bound_side.([0.0, 1.0, 5.0, 2.0], [0.0, 0.0, 0.0, 2.0], [Inf, Inf, 5.0, 2.0]) ==
          [-1, 0, 1, 0]
    @test Briareus._names(_Shift(p, x0), 2) == ["x[1]", "x[2]"]
end
