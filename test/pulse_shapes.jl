# A pulse-shape budget: grouping per-key multipliers over phases into a few shapes, and
# the shapes written out as pulses.
#
# References: exhaustive enumeration of every assignment of the keys to the groups (k^n of
# them, written here independently of the package's enumeration of partitions) for the
# least squared error in the log multipliers, and of every split of the shapes between
# the classes; the geometric mean per phase in closed form; for the pulses, the base
# pulse times the interpolated multiplier at the table's points in closed form.

using Briareus, Test

# The least squared error of grouping the rows of X into exactly k groups, over every
# assignment (groups left empty skipped).
function _brute_sse(X, k)
    n = size(X, 1)
    best = Inf
    for code in 0:(k^n - 1)
        a = k == 1 ? ones(Int, n) : digits(code; base = k, pad = n) .+ 1
        length(unique(a)) == k || continue
        sse = 0.0
        for g in 1:k
            m = findall(==(g), a)
            c = sum(X[m, :]; dims = 1) ./ length(m)
            sse += sum(abs2, X[m, :] .- c)
        end
        best = min(best, sse)
    end
    return best
end

# deterministic multipliers: three loose clusters in four phases, and noise
function _mock_multipliers(n, d; seed = 0)
    centers = [0.8 0.85 0.9 0.95; 1.0 1.02 1.05 1.08; 1.3 1.25 1.2 1.1]
    return [centers[1 + (i % 3), j] * (1 + 0.05 * sin(7.1i + 3.3j + seed)) for i in 1:n, j in 1:d]
end

@testset "pulse shapes" begin
    @testset "grouping against exhaustive enumeration" begin
        sse(g, M) = length(g.keys) * size(M, 2) * g.rms^2
        cases = [(_mock_multipliers(7, 4), 3), (_mock_multipliers(8, 2; seed = 1), 2),
                 (_mock_multipliers(9, 3; seed = 2), 3), (_mock_multipliers(8, 1; seed = 3), 4)]
        ok_exact, ok_restarts = Bool[], Bool[]
        for (M, k) in cases
            keys = ["Q$i" for i in 1:size(M, 1)]
            ref = _brute_sse(log.(M), k)
            ge = pulse_shape_groups(M, keys; shapes = k)                     # enumeration (or 1-D DP)
            gr = pulse_shape_groups(M, keys; shapes = k, exact_limit = 0)    # restarts (1-D: DP)
            push!(ok_exact, ge.exact && isapprox(sse(ge, M), ref; rtol = 1e-12, atol = 1e-15))
            push!(ok_restarts, isapprox(sse(gr, M), ref; rtol = 1e-12, atol = 1e-15))
        end
        @test all(ok_exact) && all(ok_restarts)
        # restarts are flagged as such, except in one dimension (exact by dynamic programming)
        M = _mock_multipliers(9, 3; seed = 2)
        @test !pulse_shape_groups(M, 1:9; shapes = 3, exact_limit = 0).exact
        @test pulse_shape_groups(M[:, 1:1], 1:9; shapes = 3, exact_limit = 0).exact
        # the seed fixes the result; a single restart may land in a worse local optimum,
        # never a better one than the enumeration's
        a = pulse_shape_groups(M, 1:9; shapes = 3, exact_limit = 0, seed = 7)
        b = pulse_shape_groups(M, 1:9; shapes = 3, exact_limit = 0, seed = 7)
        @test a.assignment == b.assignment && a.multipliers == b.multipliers
        one = pulse_shape_groups(M, 1:9; shapes = 3, exact_limit = 0, restarts = 1)
        @test one.rms >= pulse_shape_groups(M, 1:9; shapes = 3).rms * (1 - 1e-12)
    end

    @testset "classes, split and geometric means" begin
        M = _mock_multipliers(10, 3; seed = 4)
        keys = ["Q$i" for i in 1:10]
        cls = [i <= 4 ? :inner : :outer for i in 1:10]
        g = pulse_shape_groups(M, keys; shapes = 4, class = cls)
        # the split by exhaustive search: each class's best grouping per count
        X = log.(M)
        inn, out = X[1:4, :], X[5:10, :]
        splits = [(ki, 4 - ki) for ki in 1:3]
        tot = [_brute_sse(inn, ki) + _brute_sse(out, ko) for (ki, ko) in splits]
        j = argmin(tot)
        @test g.split == [:inner => splits[j][1], :outer => splits[j][2]]
        @test g.rms ≈ sqrt(tot[j] / 30) rtol = 1e-12
        @test [s for (s, _) in g.splits] == [[1, 3], [2, 2], [3, 1]]
        @test all(isapprox.([r for (_, r) in g.splits], sqrt.(tot ./ 30); rtol = 1e-12))
        # a class never shares a shape; shapes named per class; members and assignment agree
        @test all(c -> all(k -> cls[findfirst(==(k), keys)] === g.class[c], g.members[c]), eachindex(g.names))
        @test all(i -> g.assignment[keys[i]] in g.names[g.class .=== cls[i]], 1:10)
        @test sort(reduce(vcat, g.members)) == sort(keys)
        @test all(startswith.(string.(g.names), string.(g.class)))
        # each shape's multipliers: the geometric mean of its members' per phase
        gm = [prod(M[findfirst(==(k), keys), p] for k in g.members[c])^(1 / length(g.members[c]))
              for c in eachindex(g.names), p in 1:3]
        @test g.multipliers ≈ gm rtol = 1e-12
        # in a class, shapes numbered by increasing mean log multiplier
        @test all(c -> issorted(vec(sum(log.(g.multipliers[g.class .=== c, :]); dims = 2))), (:inner, :outer))
        # a class as a function or a dictionary; the dictionary forms of the multipliers
        @test pulse_shape_groups(M, keys; shapes = 4, class = k -> cls[parse(Int, k[2:end])]).assignment == g.assignment
        @test pulse_shape_groups(M, keys; shapes = 4, class = Dict(zip(keys, cls))).assignment == g.assignment
        byvec = Dict(k => M[i, :] for (i, k) in enumerate(keys))
        pairs = Dict((k, p) => M[i, j] for (i, k) in enumerate(keys) for (j, p) in enumerate((:picket, :foot, :rise)))
        gd = pulse_shape_groups(byvec; shapes = 4, class = Dict(zip(keys, cls)))
        gp = pulse_shape_groups(pairs; phases = (:picket, :foot, :rise), shapes = 4, class = Dict(zip(keys, cls)))
        @test gd.assignment == gp.assignment == g.assignment && isapprox(gd.rms, g.rms; rtol = 1e-12)
        # no class: one class, shapes named :shape1, ...; a shape per key fits exactly
        @test pulse_shape_groups(M, keys; shapes = 2).names == [:shape1, :shape2]
        @test pulse_shape_groups(M, keys; shapes = 10).rms == 0
        @test_throws ArgumentError pulse_shape_groups(M, keys; shapes = 1, class = cls)
        @test_throws ArgumentError pulse_shape_groups(M, keys; shapes = 11)
        @test_throws ArgumentError pulse_shape_groups(-M, keys; shapes = 2)
        @test_throws DimensionMismatch pulse_shape_groups(M, keys[1:9]; shapes = 2)
        @test occursin("rms", sprint(show, MIME"text/plain"(), g)) && occursin("inner1", sprint(show, MIME"text/plain"(), g))

        # a tolerance: the fewest shapes up to `shapes` whose best split fits within it,
        # against the best split per count by exhaustive search
        best(m) = sqrt(minimum(_brute_sse(inn, ki) + _brute_sse(out, m - ki) for ki in 1:min(4, m - 1) if m - ki <= 6) / 30)
        r2, r3, r4 = best(2), best(3), best(4)
        @test r2 > r3 > r4 > 0 && g.counts == [4 => g.rms]
        g3 = pulse_shape_groups(M, keys; shapes = 4, class = cls, tolerance = (r2 + r3) / 2)
        @test length(g3.names) == 3 && first.(g3.counts) == [2, 3]
        @test all(isapprox.(last.(g3.counts), [r2, r3]; rtol = 1e-12)) && isapprox(g3.rms, r3; rtol = 1e-12)
        @test length(pulse_shape_groups(M, keys; shapes = 4, class = cls, tolerance = r2 * (1 + 1e-9)).names) == 2
        # none within it: `shapes`, the same grouping as without a tolerance
        g0 = pulse_shape_groups(M, keys; shapes = 4, class = cls, tolerance = 0)
        @test g0.assignment == g.assignment && g0.rms == g.rms && first.(g0.counts) == [2, 3, 4]
        @test_throws ArgumentError pulse_shape_groups(M, keys; shapes = 4, tolerance = -0.01)
    end

    @testset "shapes as pulses" begin
        base = PulseShape([0.0, 1.0, 2.0, 4.0, 6.0, 8.0], [0.0, 0.4, 0.2, 0.6, 2.0, 2.0])
        t = [1.0, 2.0, 4.0]
        m = [0.5, 1.5, 1.2]
        p = scaled_pulse(base, t, m)
        # at the base's points: base times the multiplier, interpolated inside the phases
        # and held outside them
        @test p.t == base.t
        @test p.p ≈ base.p .* [0.5, 0.5, 1.5, 1.2, 1.2, 1.2] rtol = 1e-14
        # a phase time between the base's points is a point of the table
        q = scaled_pulse(base, [1.5, 4.0], [2.0, 1.0])
        @test q.t == [0.0, 1.0, 1.5, 2.0, 4.0, 6.0, 8.0]
        @test isapprox(q(1.5), base(1.5) * 2.0; rtol = 1e-14) && isapprox(q(2.0), base(2.0) * (2.0 - 0.5 / 2.5); rtol = 1e-14)
        # the peak unscaled: the multiplier ramps from the last phase's to 1 by `peak`
        r = scaled_pulse(base, t, m; peak = 5.0)
        @test r.t == [0.0, 1.0, 2.0, 4.0, 5.0, 6.0, 8.0]
        @test isapprox(r(5.0), base(5.0); rtol = 1e-14) && r(6.0) == base(6.0) && r(8.0) == base(8.0)
        @test r(4.0) ≈ base(4.0) * 1.2 rtol = 1e-14
        @test_throws ArgumentError scaled_pulse(base, t, m; peak = 3.0)
        @test_throws ArgumentError scaled_pulse(base, [2.0, 1.0], [1.0, 1.0])
        @test_throws DimensionMismatch scaled_pulse(base, t, m[1:2])

        # every shape of a grouping, from its class's base pulse, through the pulse table
        M = _mock_multipliers(10, 3; seed = 4)
        keys = ["Q$i" for i in 1:10]
        g = pulse_shape_groups(M, keys; shapes = 4, class = [i <= 4 ? :inner : :outer for i in 1:10])
        inner = PulseShape([0.0, 1.0, 2.0, 4.0, 6.0, 8.0], [0.0, 0.06, 0.05, 0.3, 2.0, 2.0])
        ps = shape_pulses(g, Dict(:inner => inner, :outer => base), t; peak = 5.0)
        @test sort(collect(Base.keys(ps))) == sort(g.names)
        @test all(i -> ps[g.names[i]] == scaled_pulse(g.class[i] === :inner ? inner : base, t,
                                                      g.multipliers[i, :]; peak = 5.0), eachindex(g.names))
        path = write_pulse_table(joinpath(mktempdir(), "shapes.csv"), ps)
        back = read_pulse_table(path)
        # (the table carries ten significant digits)
        @test all(n -> all(τ -> isapprox(back[n](τ), ps[n](τ); rtol = 1e-9, atol = 1e-15), 0:0.5:8), g.names)
        @test_throws DimensionMismatch shape_pulses(g, base, t[1:2])
        # the laser's limit checked on the table, a shape's limit per class
        @test_throws ArgumentError scaled_pulse(base, t, m; power_max = 2.0)
        @test scaled_pulse(base, t, m; peak = 5.0, power_max = 2.0) == r
        @test_throws ArgumentError shape_pulses(g, Dict(:inner => inner, :outer => base), t; peak = 5.0,
                                                power_max = Dict(:inner => 2.0, :outer => 1.0))
    end

    # The laser's limit between the phase samples. `pulse_ceiling` holds every point of the
    # written pulse (`scaled_pulse` of the base with the scales as multipliers, ramped to 1
    # at `peak`) within the limit; checked against the written pulse sampled every 1e-4 ns
    # (every point of the tables below lies on that grid), on random designs: a design
    # meets the constraint if and only if its sampled pulse stays within the limit. The box
    # of `bounds(...; pulses, times, peak)` holds the conditions on one scale alone and
    # excludes no design within the limit.
    @testset "power cap between the phases" begin
        quads = nif_quads(power = 1.0)[1:4]
        tph = [1.0, 2.25, 3.75, 5.0]
        times = (a = 1.0, b = 2.25, c = 3.75, d = 5.0, peak = 6.0)
        names = collect(keys(times))
        pk, L = 5.5, 2.0
        dense(q) = maximum(q, 0.0:1e-4:8.0)
        phased(base, pw = [fill(base(times[k]), 4) for k in names]) =
            PhasedPointing(GroupPointing(quads; by = :quad, fields = (:z,)), names, pw; scaled = names[1:4], by = b -> :all)
        # a steep rise between the last two phases to a peak at the limit inside the ramp (as
        # on a multi-phase study); a hump between two phases away from both; a pulse of 141
        # random points, a measured table's density
        rise = PulseShape([0.0, 1.0, 2.25, 3.75, 4.5, 5.0, 5.15, 5.5, 8.0], [0.1, 0.2, 0.2, 0.6, 1.2, 1.94, 2.0, 1.98, 1.98])
        hump = PulseShape([0.0, 1.0, 1.6, 2.25, 3.0, 3.75, 5.0, 8.0], [0.1, 0.3, 1.9, 0.5, 1.7, 0.6, 1.5, 1.5])
        tr = collect(0.0:0.05:7.0)
        rough = PulseShape(tr, ifelse.(tr .< pk, 0.2 .+ 1.7 .* rand(Random.Xoshiro(3), length(tr)), 1.9))
        rng = Random.Xoshiro(11)
        for (base, smax) in ((rise, 1.12), (hump, 1.3), (rough, 1.3))
            pp = phased(base)
            c = pulse_ceiling(pp; pulses = base, times, power_max = L, peak = pk)
            lo, hi = bounds(pp; scale = Inf, scale_min = 0.5, power_max = L, pulses = base, times, peak = pk)
            S = [0.5 .+ (smax - 0.5) .* rand(rng, 4) for _ in 1:300]
            meets = [all(c.f(apply(pp, [zeros(pp.n); s])) .<= c.limit .* (1 + 1e-12)) for s in S]
            within = [dense(scaled_pulse(base, tph, s; peak = pk)) <= L * (1 + 1e-12) for s in S]
            @test meets == within && 30 < count(within) < 270
            @test all(k -> !within[k] || all(lo[pp.n+1:end] .<= S[k] .<= hi[pp.n+1:end]), eachindex(S))
            # the hull keeps fewer points than the tables carry
            @test length(c.limit) < length(scaled_pulse(base, tph, ones(4); peak = pk).t)
        end
        # the study's case: the earlier scales at 1.06, the last below 1; the former rule (each
        # scale under the largest base on the interval it governs) capped the step's at
        # 2/1.94, while the written pulse stays within the limit
        pr = phased(rise)
        c = pulse_ceiling(pr; pulses = rise, times, power_max = L, peak = pk)
        _, hi = bounds(pr; scale = Inf, power_max = L, pulses = rise, times, peak = pk)
        s = [1.06, 1.06, 1.06, 0.99]
        @test L / maximum(rise, 2.25:1e-4:5.0) < 1.04 && all(s .<= hi[pr.n+1:end])
        @test all(c.f(apply(pr, [zeros(pr.n); s])) .<= c.limit) && dense(scaled_pulse(rise, tph, s; peak = pk)) <= L
        @test scaled_pulse(rise, tph, s; peak = pk, power_max = L) isa PulseShape
        # the last scale at the peak inside the ramp: at most 1, exactly
        @test isapprox(hi[end], 1.0; rtol = 1e-12)
        @test c.name === :pulse && all(==(L), c.limit) && "all, 5.15 ns" in c.labels
        # a design over the limit at the hump between the first two phases, inside the box:
        # the constraint reports it at 1.6 ns, and the table check refuses the pulse
        ph = phased(hump)
        c = pulse_ceiling(ph; pulses = hump, times, power_max = L, peak = pk)
        lo, hi = bounds(ph; scale = Inf, scale_min = 0.5, power_max = L, pulses = hump, times, peak = pk)
        s = [1.1, 1.1, 1.0, 1.0]
        v = c.f(apply(ph, [zeros(ph.n); s]))
        @test all(lo[ph.n+1:end] .<= s .<= hi[ph.n+1:end]) && c.labels[argmax(v .- c.limit)] == "all, 1.6 ns"
        @test maximum(v) ≈ 1.9 * 1.1 rtol = 1e-14
        @test_throws ArgumentError scaled_pulse(hump, tph, s; peak = pk, power_max = L)
        # generic in the number type: the components' gradient, the constant coefficients
        x = [zeros(ph.n); s]
        J = ForwardDiff.jacobian(y -> c.f(apply(ph, y)), x)
        k = findfirst(==("all, 1.6 ns"), c.labels)
        w = (1.6 - 1.0) / 1.25
        @test isapprox(J[k, ph.n+1:end], [1.9 * (1 - w), 1.9 * w, 0, 0]; rtol = 1e-12) && all(iszero, J[:, 1:ph.n])
        # a vector of times in the phases' order, and a function of the beam, alike
        c2 = pulse_ceiling(ph; pulses = b -> hump, times = collect(values(times)), power_max = b -> L, peak = pk)
        @test c2.limit == c.limit && c2.labels == c.labels && c2.f(apply(ph, x)) == v
        @test bounds(ph; scale = Inf, power_max = L, pulses = b -> hump, times = collect(values(times)), peak = pk) ==
              bounds(ph; scale = Inf, power_max = L, pulses = hump, times, peak = pk)
        # the phases' powers must be the base pulse's at their times
        pw = [fill(hump(times[k]), 4) for k in names]
        bad = phased(hump, [[pw[1] .* 1.1]; pw[2:end]])
        @test_throws ArgumentError bounds(bad; power_max = L, pulses = hump, times, peak = pk)
        @test_throws ArgumentError pulse_ceiling(bad; pulses = hump, times, power_max = L, peak = pk)
        @test_throws ArgumentError bounds(ph; power_max = L, pulses = hump)
        @test_throws ArgumentError pulse_ceiling(ph; pulses = hump, times, power_max = L, peak = 4.0)
    end
end
