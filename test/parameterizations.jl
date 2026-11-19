# Composed parameterizations: parameters held at given values (`HeldParameters`) and a
# power scale per group of beams (`PowerScales`), checked against the beams composed by
# hand (the inner parameterization applied, then powers remade) and against closed-form
# gradients.

using Briareus, ForwardDiff, LinearAlgebra, Test

@testset "composed parameterizations" begin
    quads = nif_quads(power = 8.0)
    isinner(b) = b.group in (:inner23, :inner30)
    cone2(b) = isinner(b) ? :inner : :outer
    gp = GroupPointing(quads; by = :quad, fields = (:z, :phi), mirror = false, relative = true)
    x0 = initial_parameters(gp)
    names = parameter_names(gp)
    inner = Set(b.quad for b in quads if isinner(b))
    byinner(n) = first(split(n, '.')) in inner
    held = findall(byinner, names)
    @test length(held) == 2 * 16                      # 16 inner quads, two fields each

    @testset "HeldParameters" begin
        values = x0 .+ 0.01 .* (1:length(x0))          # held values away from the start
        h = HeldParameters(gp, byinner; values)
        free = setdiff(1:length(x0), held)
        @test nparams(h) == length(free) && parameter_names(h) == names[free]
        @test initial_parameters(h) == values[free]
        # indices, a mask, names and a set of names select the same parameters
        @test HeldParameters(gp, held; values).free == h.free
        @test HeldParameters(gp, [n in held for n in eachindex(names)]; values).free == h.free
        @test HeldParameters(gp, names[held]; values).free == h.free
        @test HeldParameters(gp, Set(names[held]); values).free == h.free
        # apply: the inner parameterization at the held values and the free parameters
        y = initial_parameters(h) .+ 0.05 .* cos.(1:nparams(h))
        x = copy(values)
        x[free] = y
        @test apply(h, y) == apply(gp, x)
        # bounds: the inner box about the same point, at the free parameters
        lo, hi = bounds(gp, x; z = 0.3, phi = 0.2)
        @test bounds(h, y; z = 0.3, phi = 0.2) == (lo[free], hi[free])
        # nothing held, everything held
        @test apply(HeldParameters(gp, Int[]), x0) == apply(gp, x0)
        hall = HeldParameters(gp, n -> true; values)
        @test nparams(hall) == 0 && apply(hall, Float64[]) == apply(gp, values)
        @test_throws ArgumentError HeldParameters(gp, ["Q99T.z"])
        @test_throws DimensionMismatch HeldParameters(gp, [true, false])
        @test_throws DimensionMismatch HeldParameters(gp, held; values = [0.0])
        @test_throws DimensionMismatch apply(h, zeros(3))
    end

    @testset "PowerScales" begin
        p = GroupPointing(quads; fields = (:z, :r))
        xp = initial_parameters(p)
        ps = PowerScales(p; by = cone2)
        @test Set(ps.keys) == Set([:inner, :outer])
        @test nparams(ps) == length(xp) + 2
        @test parameter_names(ps) == [parameter_names(p); string.(ps.keys, ".power_scale")]
        @test initial_parameters(ps) == [xp; 1.0; 1.0]
        @test apply(ps, initial_parameters(ps)) == apply(p, xp)
        s = Dict(:inner => 0.8, :outer => 1.1)
        xs = xp .+ 0.1
        y = [xs; [s[k] for k in ps.keys]]
        byhand = [remake(b; power = b.power * s[cone2(b)]) for b in apply(p, xs)]
        @test apply(ps, y) == byhand
        # the cone fraction follows the scales (closed form: 16 inner and 32 outer quads)
        P = [b.power for b in apply(ps, y)]
        @test sum(P[isinner.(quads)]) / sum(P) ≈ 16 * 0.8 / (16 * 0.8 + 32 * 1.1) rtol = 1e-14
        # bounds: the inner box, then each scale in y0 (1 ± scale) within [min, max]
        lo, hi = bounds(ps, y; scale = 0.5, scale_max = 1.0, z = 0.4)
        plo, phi = bounds(p, xs; z = 0.4)
        sv = [s[k] for k in ps.keys]
        @test lo == [plo; 0.5 .* sv] && hi == [phi; min.(1.5 .* sv, 1.0)]
        @test_throws ArgumentError bounds(ps, y; scale_min = 2.0)
        # by cone (group), by quad
        @test PowerScales(p).keys == unique(b.group for b in quads)
        @test length(PowerScales(p; by = :quad).keys) == 48
        @test_throws ArgumentError PowerScales(GroupPointing(omega_beams()); by = :quad)
        @test_throws DimensionMismatch apply(ps, xp)
    end

    @testset "compositions, gradients, scenarios" begin
        # held inner aims, free outer aims and two cone scales
        values = x0 .+ 0.02
        c = HeldParameters(PowerScales(gp; by = cone2), byinner; values = [values; 1.0; 1.0])
        free = setdiff(1:length(x0), held)
        @test parameter_names(c)[end-1:end] == string.(c.inner.keys, ".power_scale")
        y = [values[free] .- 0.03; 0.9; 1.2]
        s = Dict(zip(c.inner.keys, y[end-1:end]))
        x = copy(values)
        x[free] = y[1:end-2]
        byhand = [remake(b; power = b.power * s[cone2(b)]) for b in apply(gp, x)]
        @test apply(c, y) == byhand
        # the other order of composition: the same beams
        c2 = PowerScales(HeldParameters(gp, byinner; values); by = cone2)
        @test c2.keys == c.inner.keys && apply(c2, y) == byhand
        # gradients (closed forms): Σ aim_z changes by 1 per free z shift (relative, not
        # mirrored, one beam per quad) and not at all with a phi shift or a scale; the
        # total power by the power each scale multiplies
        g = ForwardDiff.gradient(v -> sum(b -> b.aim[3], apply(c, v)), y)
        nf = parameter_names(c)
        @test g == [endswith(n, ".z") ? 1.0 : 0.0 for n in nf]
        gP = ForwardDiff.gradient(v -> sum(b -> b.power, apply(c, v)), y)
        @test all(iszero, gP[1:end-2])
        @test gP[end-1:end] ≈ [8.0 * count(b -> cone2(b) === k, quads) for k in c.inner.keys] rtol = 1e-14
        # a FramePointing inside: the frames' aims, scaled powers
        fp = FramePointing(quads, CHAMBER; fields = (:u, :v, :w))
        fps = HeldParameters(PowerScales(fp; by = cone2), n -> endswith(n, ".w"))
        yf = initial_parameters(fps) .+ 0.01
        xf = yf[1:end-2]
        full = initial_parameters(fp)
        full[setdiff(1:length(full), findall(n -> endswith(n, ".w"), parameter_names(fp)))] = xf
        sf = Dict(zip(fps.inner.keys, yf[end-1:end]))
        @test apply(fps, yf) == [remake(b; power = b.power * sf[cone2(b)]) for b in apply(fp, full)]

        # in a scenario and its subset: evaluations of the composed beams
        t = Target([foil(:foil, (0, 0, 0), (0, 0, 1); width = 8.0, n = 8)])
        sc = Scenario(t, quads)
        keep = findall(b -> nif_bundle(b) != "31", quads)
        kept = quads[keep]
        ck = HeldParameters(PowerScales(GroupPointing(kept; by = :quad, fields = (:z, :phi), mirror = false,
                                                      relative = true); by = cone2),
                            n -> first(split(n, '.')) in inner)
        sk = subset(sc, keep; parameterization = ck)
        yk = initial_parameters(ck) .+ 0.05
        @test evaluate(sk, yk).deposition.element == evaluate(subset(sc, keep), apply(ck, yk)).deposition.element
        # and the optimizer on its default box, which `bounds` gives the composition
        objective(ev) = sum(abs2, ev.deposition.element)
        r = optimize_pointing(sk, objective; iterations = 2)
        @test r.names == parameter_names(ck) && length(r.x) == nparams(ck)
        @test r.value <= r.initial_value
        @test occursin("$(nparams(ck)) parameters", sprint(show, sk))
    end

    @testset "subsets and named values" begin
        # four beams per quad by quad, held inner quads and a power scale per quad; kept:
        # all but bundle 31 (its quads and their scales dropped) and two beams of an
        # inner quad (its parameters kept, held at their values)
        beams = nif_beams(power = 2.0)
        gq = GroupPointing(beams; by = :quad, fields = (:z, :phi), mirror = false, relative = true)
        nq = nparams(gq)
        values = [initial_parameters(gq) .+ 0.01 .* (1:nq); 1.0 .+ 0.001 .* (1:48)]
        c = HeldParameters(PowerScales(gq; by = :quad), byinner; values)
        part = beams[findfirst(isinner, beams)].quad
        @test nif_bundle(part) != "31"
        keep = [i for (i, b) in enumerate(beams)
                if nif_bundle(b) != "31" && !(i in findall(x -> x.quad == part, beams)[1:2])]
        t = Target([foil(:foil, (0, 0, 0), (0, 0, 1); width = 8.0, n = 2)])
        sk = subset(Scenario(t, c), keep)
        sub = sk.param
        # by hand: the kept beams by quad, a scale per kept quad, the old values by name
        byname = Dict(zip(parameter_names(c.inner), values))
        ps = PowerScales(GroupPointing(beams[keep]; by = :quad, fields = (:z, :phi), mirror = false,
                                       relative = true); by = :quad)
        hand = HeldParameters(ps, byinner; values = [byname[n] for n in parameter_names(ps)])
        @test sub isa HeldParameters && sub.inner isa PowerScales
        @test parameter_names(sub) == parameter_names(hand) && sub.free == hand.free
        @test initial_parameters(sub) == initial_parameters(hand)
        y = initial_parameters(hand) .+ 0.02 .* sin.(1:nparams(hand))
        @test apply(sub, y) == apply(hand, y)
        @test evaluate(sk, y).beams == apply(hand, y)
        # the partly dropped quad keeps its parameters, held at their values; bundle 31's
        # quads and their scales are gone
        inner_names = parameter_names(sub.inner)
        @test all(n -> n in inner_names && !(n in parameter_names(sub)) &&
                       sub.values[findfirst(==(n), inner_names)] == byname[n],
                  ("$part.z", "$part.phi", "$part.power_scale"))
        @test !any(n -> nif_bundle(first(split(n, '.'))) == "31", inner_names)
        @test length(sub.inner.keys) == 48 - count(q -> nif_bundle(q) == "31", unique(b.quad for b in beams))
        # a power scale per cone alone: as built by hand on the kept beams
        sp = subset(Scenario(t, PowerScales(gq; by = cone2)), keep).param
        hp = PowerScales(GroupPointing(beams[keep]; by = :quad, fields = (:z, :phi), mirror = false,
                                       relative = true); by = cone2)
        @test parameter_names(sp) == parameter_names(hp) && sp.group == hp.group
        yp = initial_parameters(hp) .+ 0.03
        @test apply(sp, yp) == apply(hp, yp)
        # a group (cone) that loses every beam loses its scale
        ss = Briareus._subset_param(PowerScales(gq), findall(b -> b.group !== :outer44, beams))
        @test ss.keys == [g for g in unique(b.group for b in beams) if g !== :outer44]
        # a GroupPointing's own subset keeps its groups too: the full pointing's names, as
        # built by hand on the kept beams by quad
        sg = subset(Scenario(t, gq), keep).param
        hg = GroupPointing(beams[keep]; by = :quad, fields = (:z, :phi), mirror = false, relative = true)
        @test parameter_names(sg) == parameter_names(hg) && all(in(parameter_names(gq)), parameter_names(sg))
        @test apply(sg, initial_parameters(hg) .+ 0.04) == apply(hg, initial_parameters(hg) .+ 0.04)

        # named values: a quad's fields and its scale by the quad's name, passed through
        # to the inner parameterization; a held parameter is refused
        oq = beams[findfirst(!isinner, beams)].quad
        free_names = parameter_names(c)
        yn = parameters(c, oq => (z = 0.3, power_scale = 1.2))
        yh = copy(initial_parameters(c))
        yh[findfirst(==("$oq.z"), free_names)] = 0.3
        yh[findfirst(==("$oq.power_scale"), free_names)] = 1.2
        @test yn == yh
        @test parameters(c) == initial_parameters(c)
        @test_throws ArgumentError parameters(c, part => (z = 0.3,))                # held
        @test_throws ArgumentError parameters(c, oq => (r = 0.3,))                  # not a field
        @test_throws ArgumentError parameters(c, "Q99T" => (power_scale = 1.0,))    # no such scale
        # scales keyed by cone beside groups keyed by cone: one keyword sets both
        pc = PowerScales(GroupPointing(quads; fields = (:z, :r)))
        xn = parameters(pc; inner23 = (z = 3.0, power_scale = 0.9), outer50 = (power_scale = 1.1,))
        xh = initial_parameters(pc)
        pn = parameter_names(pc)
        xh[findfirst(==("inner23.z"), pn)] = 3.0
        xh[findfirst(==("inner23.power_scale"), pn)] = 0.9
        xh[findfirst(==("outer50.power_scale"), pn)] = 1.1
        @test xn == xh
        # the other order: scales outside, the hold inside
        c2 = PowerScales(HeldParameters(gq, byinner; values = values[1:nq]); by = cone2)
        y2 = parameters(c2, oq => (phi = -0.2,); outer = (power_scale = 0.8,))
        n2 = parameter_names(c2)
        @test y2[findfirst(==("$oq.phi"), n2)] == -0.2 && y2[findfirst(==("outer.power_scale"), n2)] == 0.8 &&
              y2[findfirst(==("inner.power_scale"), n2)] == 1.0
    end
end
