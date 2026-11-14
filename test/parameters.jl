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

# A quad as the unit of pointing (`unit = :quad`): the radial and tangential offsets in the
# quad's port azimuth, so a quad's beams that share their values share one aim point.
@testset "GroupPointing: a quad as the pointing unit" begin
    quads, beams = nif_quads(power = 8.0), nif_beams(power = 2.0)
    vals = (inner23 = (z = 3.09, r = -1.09, phi = 0.2), inner30 = (z = 3.63, r = -0.6, phi = -0.1),
            outer44 = (z = 5.2, r = 0.4, phi = 0.0), outer50 = (z = 5.01, r = 0.0, phi = 0.3))
    fields = (:z, :r, :phi)
    # quads held as one beam: the same either way, bit for bit
    pb, pq = GroupPointing(quads; fields), GroupPointing(quads; fields, unit = :quad)
    x = parameters(pb; vals...)
    @test initial_parameters(pb) == initial_parameters(pq) && pq.radial == pb.radial
    @test [b.aim for b in apply(pq, x)] == [b.aim for b in apply(pb, x)]
    # four beams per quad: each quad's beams at one point, where its quad beam aims (the mean
    # of a quad's four ports is its port to 0.004°: under a micron here)
    p4 = GroupPointing(beams; fields, unit = :quad)
    x4 = parameters(p4; vals...)
    b4 = apply(p4, x4)
    qaim = Dict(b.quad => b.aim for b in apply(pq, x))
    byquad = Dict{String,Vector{Int}}()
    foreach(i -> push!(get!(byquad, beams[i].quad, Int[]), i), eachindex(beams))
    spread(bs, idx) = maximum(norm(bs[i].aim - bs[j].aim) for i in idx, j in idx)
    @test length(byquad) == 48 && all(v -> length(v) == 4, values(byquad))
    @test all(idx -> spread(b4, idx) < 1e-12, values(byquad))
    @test maximum(norm(b4[idx[1]].aim - qaim[q]) for (q, idx) in byquad) < 1e-3
    # and the values read back: the applied beams give the parameters again
    @test initial_parameters(GroupPointing(b4; fields, unit = :quad)) ≈ x4 rtol = 1e-12 atol = 1e-12
    # per beam (the default) the same values put a quad's aims r sin Δφ apart: the 23.5°
    # cone's beams are up to 5.73° from their quad's azimuth (closed form, to the last digit
    # of the port table's angles)
    b1 = apply(GroupPointing(beams; fields), parameters(GroupPointing(beams; fields); vals...))
    inner = [idx for (q, idx) in byquad if beams[idx[1]].group === :inner23]
    φ(v) = atand(v[2], v[1])
    δ(idx) = maximum(abs(mod(φ(beams[i].port) - φ(sum(beams[j].port for j in idx)) + 180, 360) - 180) for i in idx)
    @test all(idx -> spread(b1, idx) > 0.1, inner)
    @test maximum(δ, inner) ≈ 5.73 atol = 0.01
    idx = inner[argmax(δ.(inner))]
    # (two beams at ±δ from the quad's azimuth, both offsets turned by δ: 2 |(r, phi)| sin δ)
    @test spread(b1, idx) ≈ 2 * hypot(1.09, 0.2) * sind(δ(idx)) rtol = 0.02
    # per quad and per beam rows, and a rigid shift (`relative`), keep a quad's beams together
    for by in (:quad, :beam)
        pr = GroupPointing(b4; by, fields, unit = :quad)
        xr = initial_parameters(pr)
        @test all(idx -> spread(apply(pr, xr .+ 0.25), idx) < 1e-12, values(byquad))
    end
    prel = GroupPointing(b4; fields = (:r, :phi), unit = :quad, relative = true)
    @test all(idx -> spread(apply(prel, fill(0.3, Briareus.nparams(prel))), idx) < 1e-12, values(byquad))
    # differentiable, and kept by a subset
    g = ForwardDiff.gradient(y -> sum(b -> sum(b.aim), apply(p4, y)), x4)
    @test all(isfinite, g) && any(!=(0), g)
    @test Briareus._subset_param(p4, 1:8).unit === :quad && pb.unit === :beam
    # a beam without a quad name keeps its own azimuth
    om = omega_beams()
    @test GroupPointing(om; fields, unit = :quad).radial == GroupPointing(om; fields).radial
    @test_throws ArgumentError GroupPointing(quads; unit = :cone)
end

@testset "drive level, optimization result, report" begin
    t = hohlraum(nz = 8, nphi = 16, capsule_ntheta = 6, capsule_nphi = 8)
    p = GroupPointing(nif_quads(power = 8.0); fields = (:z, :r, :power))
    s = Scenario(t, p; observe = :capsule)
    x0 = parameters(p; inner23 = (z = 3.09, r = -1.09), inner30 = (z = 3.63, r = -0.60),
                    outer44 = (z = 5.20,), outer50 = (z = 5.01,))
    ev = evaluate(s, x0)
    @test mean_flux(ev) == flux_moments(ev)[1]
    # a given model fixes the view factors: viewfactor_kw beside it is an error
    @test_throws ArgumentError Scenario(t, p; observe = :capsule, model = s.model,
                                        viewfactor_kw = (; closed = false))

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
