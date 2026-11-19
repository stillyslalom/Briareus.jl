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
    end
end
