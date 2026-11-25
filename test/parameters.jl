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
    # Any other parameterization gets the operator for its beams at the initial
    # parameters.
    c = Scenario(t, _Shift(p, x))
    @test evaluate(c, zeros(nparams(p))).deposition.element ==
          evaluate(Scenario(t, apply(p, x)), apply(p, x)).deposition.element
end

# A quad as the unit of pointing (`unit = :quad`, the default): the radial and tangential
# offsets in the quad's port azimuth, so a quad's beams that share their values share one
# aim point.
@testset "GroupPointing: a quad as the pointing unit" begin
    quads, beams = nif_quads(power = 8.0), nif_beams(power = 2.0)
    vals = (inner23 = (z = 3.09, r = -1.09, phi = 0.2), inner30 = (z = 3.63, r = -0.6, phi = -0.1),
            outer44 = (z = 5.2, r = 0.4, phi = 0.0), outer50 = (z = 5.01, r = 0.0, phi = 0.3))
    fields = (:z, :r, :phi)
    # Quads held as one beam give identical results either way.
    pb, pq = GroupPointing(quads; fields, unit = :beam), GroupPointing(quads; fields)
    x = parameters(pb; vals...)
    @test initial_parameters(pb) == initial_parameters(pq) && pq.radial == pb.radial
    @test [b.aim for b in apply(pq, x)] == [b.aim for b in apply(pb, x)]
    # With four beams per quad, the beams of each quad aim at one point, the aim of its quad
    # beam. The mean of the four ports of a quad is its port to 0.004°, which is under a
    # micron in this geometry.
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
    # The applied beams give the parameters again.
    @test initial_parameters(GroupPointing(b4; fields, unit = :quad)) ≈ x4 rtol = 1e-12 atol = 1e-12
    # per beam (`unit = :beam`) the same values put a quad's aims r sin Δφ apart: the 23.5°
    # cone's beams are up to 5.73° from their quad's azimuth (closed form, to the last digit
    # of the port table's angles)
    p1 = GroupPointing(beams; fields, unit = :beam)
    b1 = apply(p1, parameters(p1; vals...))
    inner = [idx for (q, idx) in byquad if beams[idx[1]].group === :inner23]
    φ(v) = atand(v[2], v[1])
    δ(idx) = maximum(abs(mod(φ(beams[i].port) - φ(sum(beams[j].port for j in idx)) + 180, 360) - 180) for i in idx)
    @test all(idx -> spread(b1, idx) > 0.1, inner)
    @test maximum(δ, inner) ≈ 5.73 atol = 0.01
    idx = inner[argmax(δ.(inner))]
    # Two beams lie at ±δ from the quad azimuth, with both offsets turned by δ, so the
    # spread is 2 |(r, phi)| sin δ.
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
    @test Briareus._subset_param(p4, 1:8).unit === :quad && pb.unit === :beam && pq.unit === :quad
    @test Briareus._subset_param(p1, 1:8).unit === :beam
    # a beam without a quad name keeps its azimuth
    om = omega_beams()
    @test GroupPointing(om; fields).radial == GroupPointing(om; fields, unit = :beam).radial
    @test_throws ArgumentError GroupPointing(quads; unit = :cone)
end

# A quad's beams split about its aim point (`:split`, in the quad's frame), against beams
# split by hand: each beam's aim moved along its quad's φ̂ (the mean port's azimuth, from
# the port table) by ±split with the side of its port, and against the mock HYBRID-E
# pointing's outer-cone split (`hybride_mock_pointing`, another construction of the same
# rule: its `s` crosses the beams, `split = -s`).
@testset "GroupPointing: a quad's split" begin
    beams = nif_beams(power = 2.0)
    byquad = Dict{String,Vector{Int}}()
    foreach(i -> push!(get!(byquad, beams[i].quad, Int[]), i), eachindex(beams))
    φq = Dict(q => (m = sum(beams[i].port for i in idx); atan(m[2], m[1])) for (q, idx) in byquad)
    side(b) = sign(mod(atan(b.port[2], b.port[1]) - φq[b.quad] + π, 2π) - π)
    eφ(b) = Vec3(-sin(φq[b.quad]), cos(φq[b.quad]), 0.0)
    er(b) = Vec3(cos(φq[b.quad]), sin(φq[b.quad]), 0.0)
    @test all(b -> side(b) != 0, beams)
    # The absolute parameterization has one (z, r, phi, split) per quad.
    p = GroupPointing(beams; by = :quad, fields = (:z, :r, :phi, :split))
    x = [(4.5, -0.5, 0.2, 0.35)[mod1(i, 4)] + 0.01 * sin(i) for i in 1:nparams(p)]
    at(q, f) = x[findfirst(==(string(q, '.', f)), parameter_names(p))]
    hand = [at(b.quad, "r") * er(b) + (at(b.quad, "phi") + side(b) * at(b.quad, "split")) * eφ(b) +
            Vec3(0, 0, sign(b.port[3]) * at(b.quad, "z")) for b in beams]
    bs = apply(p, x)
    @test maximum(norm(b.aim - h) for (b, h) in zip(bs, hand)) < 1e-12
    # the parameters read back from the split beams, and the template's split as the start
    @test maximum(abs.(parameters(p, bs) .- x)) < 1e-12
    @test maximum(abs.(initial_parameters(GroupPointing(bs; by = :quad, fields = (:z, :r, :phi, :split))) .- x)) < 1e-12
    # the mock HYBRID-E pointing per cone: z, r and split = -s
    hyb = hybride_mock_pointing(beams)
    pc = GroupPointing(beams; fields = (:z, :r, :split))
    xc = parameters(pc; (c => (z = v.z, r = v.r, split = -v.s) for (c, v) in Briareus.HYBRIDE_MOCK_POINTING)...)
    @test maximum(norm(a.aim - b.aim) for (a, b) in zip(apply(pc, xc), hyb)) < 1e-12
    ph = GroupPointing(hyb; fields = (:z, :r, :split))
    @test maximum(abs.(initial_parameters(ph) .- xc)) < 1e-12
    # with the split held at the template's, the other fields keep it: the hybrid's quads
    # moved in z keep their ±s
    xz = initial_parameters(ph) .+ [f == "z" ? 0.1 : 0.0 for f in last.(split.(parameter_names(ph), '.'))]
    @test maximum(norm(a.aim - b.aim - Vec3(0, 0, 0.1 * sign(b.port[3]))) for (a, b) in zip(apply(ph, xz), hyb)) < 1e-12
    # The relative parameterization shifts the split of each beam and keeps the
    # differences between the members.
    pr = GroupPointing(hyb; by = :quad, fields = (:phi, :split), relative = true)
    δ = [0.02 * cos(i) for i in 1:nparams(pr)]
    atr(q, f) = δ[findfirst(==(string(q, '.', f)), parameter_names(pr))]
    @test initial_parameters(pr) == zeros(nparams(pr))
    @test maximum(norm(a.aim - (b.aim + (atr(b.quad, "phi") + side(b) * atr(b.quad, "split")) * eφ(b)))
                  for (a, b) in zip(apply(pr, δ), hyb)) < 1e-12
    @test maximum(abs.(parameters(pr, apply(pr, δ)) .- δ)) < 1e-12
    # The aim is differentiable, and d aim / d split is exactly the side of the beam along
    # the φ̂ of its quad.
    @test parameter_names(p)[4] == beams[1].quad * ".split"
    J = ForwardDiff.jacobian(s -> apply(p, [i == 4 ? s[1] : oftype(s[1], x[i]) for i in eachindex(x)])[1].aim, [x[4]])
    @test J[:, 1] ≈ side(beams[1]) * eφ(beams[1]) rtol = 1e-14
    # A quad held as one beam is not split. A split per beam, without quads or in the
    # beam frame, is refused, and spherical fields exclude it.
    q = nif_quads()
    pq = GroupPointing(q; by = :quad, fields = (:z, :split))
    @test [b.aim for b in apply(pq, fill(0.3, nparams(pq)))] == [b.aim for b in apply(pq, repeat([0.3, 0.0], 48))]
    @test_throws ArgumentError GroupPointing(beams; fields = (:split,), unit = :beam)
    @test_throws ArgumentError GroupPointing(omega_beams(); fields = (:split,))
    @test_throws ArgumentError GroupPointing(beams; fields = (:theta_s, :split), radius = 1.0)
    @test bounds(p, x; split = 0.2)[2][4] == x[4] + 0.2
end

# Parameters read from beams (`parameters(p, beams)`): a quad re-pointing handed to a
# per-beam parameterization on the same template, and back to quads as group means.
@testset "parameters from beams" begin
    tmpl = hybride_mock_pointing(nif_beams())          # outer quads' beams split ±0.5 mm
    kw = (fields = (:z, :phi), mirror = false, relative = true, unit = :quad)
    quad, per = GroupPointing(tmpl; by = :quad, kw...), GroupPointing(tmpl; by = :beam, kw...)
    xq = [0.3 * sin(1.7i) for i in 1:nparams(quad)]
    bq = apply(quad, xq)
    xb = parameters(per, bq)
    # each beam takes its quad's shift (by name), and the beams come back
    at(p, x, n) = x[findfirst(==(n), parameter_names(p))]
    want = [at(quad, xq, string(tmpl[findfirst(b -> b.name == split(n, '.')[1], tmpl)].quad, '.', split(n, '.')[2]))
            for n in parameter_names(per)]
    @test maximum(abs.(xb .- want)) < 1e-12
    @test maximum(norm(a.aim - b.aim) for (a, b) in zip(apply(per, xb), bq)) < 1e-12
    @test parameters(per, reverse(bq)) == xb
    # Per-beam shifts taken by quad give the mean shift of each quad.
    δ = [0.05 * cos(2.3i) for i in eachindex(xb)]
    xm = parameters(quad, apply(per, xb .+ δ))
    mean_by_quad = [sum(at(per, xb .+ δ, string(tmpl[i].name, '.', split(n, '.')[2]))
                        for i in findall(b -> b.quad == split(n, '.')[1], tmpl)) / 4
                    for n in parameter_names(quad)]
    @test maximum(abs.(xm .- mean_by_quad)) < 1e-12
    # With absolute fields the template gives the initial parameters, and an applied
    # vector gives itself.
    p = GroupPointing(tmpl; fields = (:z, :r, :defocus, :power, :spot))
    x0 = initial_parameters(p)
    @test maximum(abs.(parameters(p, tmpl) .- x0)) < 1e-12
    x1 = x0 .+ [0.1 * sin(i) for i in eachindex(x0)]
    @test maximum(abs.(parameters(p, apply(p, x1)) .- x1)) < 1e-12
    @test_throws ArgumentError parameters(per, bq[2:end])
    @test_throws ArgumentError parameters(GroupPointing(nif_quads(); fields = (:theta_s,), radius = 1.1),
                                          nif_quads())
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

    # The report has one line per part above `minor` of the delivered power and sums the
    # rest.
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
    @test Briareus._parameter_names(_Shift(p, x0), 2) == ["x[1]", "x[2]"]
end
