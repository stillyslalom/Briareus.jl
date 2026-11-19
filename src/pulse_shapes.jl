# A pulse-shape budget: per-key power multipliers over a pulse's phases (one per quad and
# phase, from a re-pointing run with power scales) grouped into a few distinct shapes, each
# of which needs facility validation, and the shapes written out as pulses.

"""
    PulseShapeGroups

Multipliers per key and phase grouped into shapes, from [`pulse_shape_groups`](@ref).
Fields:

- `keys`: the keys (e.g. quad names), in the order of the input;
- `names`: the shapes' names, `Symbol(class, i)` per class (`:shape1`, … without
  classes), numbered within a class by increasing mean log multiplier;
- `class`: each shape's class;
- `multipliers`: shapes × phases, each shape's multiplier per phase, the geometric mean
  of its members';
- `members`: each shape's keys, in input order;
- `assignment`: a `Dict` from each key to its shape's name;
- `rms`: the fit's root-mean-square residual in the natural log of the multipliers, over
  every key and phase;
- `split`: the number of shapes per class, as `class => n` pairs in order of first
  appearance; `splits`: every admissible split of the chosen number of shapes (counts
  per class in that order) and its rms;
- `counts`: each number of shapes tried and the rms of its best split, in increasing
  order (`shapes` alone without a `tolerance`);
- `exact`: `true` when every class's grouping is the exact least-squares optimum (by
  enumeration or, with one phase, by dynamic programming), `false` where seeded restarts
  found it.
"""
struct PulseShapeGroups{K,C}
    keys::Vector{K}
    names::Vector{Symbol}
    class::Vector{C}
    multipliers::Matrix{Float64}
    members::Vector{Vector{K}}
    assignment::Dict{K,Symbol}
    rms::Float64
    split::Vector{Pair{C,Int}}
    splits::Vector{Pair{Vector{Int},Float64}}
    counts::Vector{Pair{Int,Float64}}
    exact::Bool
end

"""
    pulse_shape_groups(multipliers::AbstractMatrix, keys; shapes, tolerance = nothing,
                       class = nothing, seed = 1, restarts = 100,
                       exact_limit = 10^6) -> PulseShapeGroups
    pulse_shape_groups(multipliers::AbstractDict; phases = nothing, kwargs...)

Group per-key power multipliers over a pulse's phases into `shapes` distinct pulse shapes:
k-means on the natural log of the multipliers (a key's point has one coordinate per
phase), so that a shape's multiplier per phase is the geometric mean of its members' and
the fit is least squares in the relative error. Rows of the matrix are `keys` (e.g. quad
names) and columns phases; a dictionary maps each key to its vector of multipliers per
phase, or, with `phases` (the phases in order), each `(key, phase)` pair to a multiplier,
the layout of power scales keyed per quad and phase. Keys of a dictionary are sorted (by
their string) for a deterministic order. Every multiplier must be positive.

`class` (a function of the key, a dictionary or a vector aligned with the keys; `nothing`:
one class) separates keys that never share a shape, such as the inner and the outer cones'
quads, which follow different base pulses. Each class is grouped on its own, with at least
one shape and at most one per key, and the `shapes` are split between the classes by the
smallest total squared error over every admissible split (ties to the first in order).

`tolerance` (an rms in the natural log of the multipliers; `nothing`, the default, for a
fixed count) makes `shapes` an upper bound: the result is the smallest number of shapes,
from one per class up to `shapes`, whose best split fits with an rms at or below
`tolerance`, or `shapes` when none does (`counts` lists the rms of each number tried).
Shapes that differ by less than the laser delivers distinctly need not be separate
pulses. For scale, not a default: NIF's laser performance model predicts total energies
"within 2–3 percent" and gives "beam-to-beam energy balance better than 5 percent for
full system shots" (NIF User Guide, section 4.8, Laser Performance Operations Model,
<https://nifuserguide.llnl.gov/home/4-laser-system/48-laser-performance-operations-model>;
the page states no statistic for the balance); a relative error ε is about ε in the log,
so 0.02–0.05 spans those figures. The choice of tolerance is the study's.

Method, per class of `n` keys into `k` groups: exact where the size allows, by
enumerating every partition into `k` groups when there are at most `exact_limit` of them
(the Stirling number S(n, k)), or by dynamic programming over the sorted values when there
is a single phase (an optimal one-dimensional grouping is contiguous); otherwise the best
of `restarts` k-means++ seedings refined by single-key transfers until no transfer
lowers the squared error (Hartigan's rule, a local optimum at least as good as Lloyd's).
The seedings draw from a SplitMix64 stream started at `seed`, so the result is the same on
every Julia version and machine. Restarts are a heuristic: the result is the global
optimum only with high probability (it matched exhaustive enumeration on every small case
tested).

[`shape_pulses`](@ref) writes the shapes out as pulses, and the `assignment` ties the
power scales of a re-optimization per shape ([`PowerScales`](@ref) with
`by = b -> groups.assignment[b.quad]`).

```julia
inner = Set(b.quad for b in beams if b.group in (:inner23, :inner30))
g = pulse_shape_groups(m, quad_names; shapes = 6, class = q -> q in inner ? :inner : :outer)
g.split, g.rms, g.assignment["Q15T"]
# the fewest shapes, up to 8, that fit within 3 % rms
g3 = pulse_shape_groups(m, quad_names; shapes = 8, tolerance = 0.03,
                        class = q -> q in inner ? :inner : :outer)
length(g3.names), g3.counts
```
"""
function pulse_shape_groups(M::AbstractMatrix{<:Real}, keys::AbstractVector; shapes::Integer,
                            class = nothing, seed::Integer = 1, restarts::Integer = 100,
                            exact_limit::Real = 10^6, tolerance::Union{Nothing,Real} = nothing)
    n, d = size(M)
    tolerance === nothing || tolerance >= 0 || throw(ArgumentError("tolerance must not be negative"))
    length(keys) == n || throw(DimensionMismatch("$(length(keys)) keys for $n rows of multipliers"))
    allunique(keys) || throw(ArgumentError("the keys must be distinct"))
    n > 0 && d > 0 || throw(ArgumentError("no multipliers"))
    all(x -> isfinite(x) && x > 0, M) || throw(ArgumentError("multipliers must be positive and finite"))
    restarts >= 1 || throw(ArgumentError("restarts must be at least 1"))
    # the class of each key, resolved here: the grouping below takes no callable
    cl = class === nothing ? fill(nothing, n) :
         class isa AbstractVector ? collect(class) :
         class isa AbstractDict ? [class[k] for k in keys] : [class(k) for k in keys]
    length(cl) == n || throw(DimensionMismatch("$(length(cl)) classes for $n keys"))
    return _pulse_shape_groups(log.(Float64.(M)), collect(keys), cl, Int(shapes), UInt64(seed),
                               Int(restarts), Float64(exact_limit),
                               tolerance === nothing ? nothing : Float64(tolerance))
end

function _pulse_shape_groups(X::Matrix{Float64}, keys::Vector{K}, cl::Vector{C}, shapes::Int,
                             seed::UInt64, restarts::Int, exact_limit::Float64,
                             tolerance::Union{Nothing,Float64}) where {K,C}
    n, d = size(X)
    classes = unique(cl)
    rows = [findall(isequal(c), cl) for c in classes]
    sizes = length.(rows)
    nc = length(classes)
    nc <= shapes <= n ||
        throw(ArgumentError("$shapes shapes for $n keys in $nc classes: at least one shape per class, at most one per key"))
    # each class's best grouping into every number of groups it can take
    kmax = [min(sizes[j], shapes - nc + 1) for j in 1:nc]
    fits = [[_kmeans_log(X[rows[j], :], k, seed, restarts, exact_limit) for k in 1:kmax[j]] for j in 1:nc]
    # every split of a count of shapes between the classes, and its rms: the count
    # `shapes`, or with a tolerance the smallest count up to it whose best split meets it
    local splits, best, bestsse
    counts = Pair{Int,Float64}[]
    for m in (tolerance === nothing ? shapes : nc):shapes
        splits = Pair{Vector{Int},Float64}[]
        best, bestsse = Int[], Inf
        for s in _compositions(m, [min(sizes[j], m - nc + 1) for j in 1:nc])
            sse = sum(fits[j][s[j]].sse for j in 1:nc)
            push!(splits, s => sqrt(sse / (n * d)))
            sse < bestsse && ((best, bestsse) = (s, sse))
        end
        push!(counts, m => sqrt(bestsse / (n * d)))
        tolerance === nothing || last(counts[end]) <= tolerance && break
    end
    names, shclass, mult, members = Symbol[], C[], Vector{Float64}[], Vector{K}[]
    assignment = Dict{K,Symbol}()
    exact = true
    for j in 1:nc
        f = fits[j][best[j]]
        exact &= f.exact
        # shapes in a class by increasing mean log multiplier
        order = sortperm([sum(f.centers[g, :]) for g in 1:best[j]])
        for (i, g) in enumerate(order)
            name = classes[j] === nothing ? Symbol(:shape, i) : Symbol(classes[j], i)
            mem = rows[j][f.assignment .== g]
            push!(names, name)
            push!(shclass, classes[j])
            push!(mult, exp.(vec(sum(X[mem, :]; dims = 1)) ./ length(mem)))
            push!(members, collect(K, keys[mem]))
            foreach(m -> assignment[keys[m]] = name, mem)
        end
    end
    allunique(names) || throw(ArgumentError("class names give clashing shape names: $(names)"))
    return PulseShapeGroups{K,C}(keys, names, shclass, permutedims(reduce(hcat, mult)),
                                 members, assignment, sqrt(bestsse / (n * d)),
                                 [classes[j] => best[j] for j in 1:nc], splits, counts, exact)
end

function pulse_shape_groups(d::AbstractDict; phases = nothing, kwargs...)
    if phases === nothing
        ks = sort!(collect(keys(d)); by = string)
        vs = [collect(Float64, d[k]) for k in ks]
        allequal(length.(vs)) || throw(DimensionMismatch("keys carry different numbers of phases"))
        return pulse_shape_groups(permutedims(reduce(hcat, vs)), ks; kwargs...)
    end
    ph = collect(phases)
    ks = sort!(unique(first(k) for k in keys(d)); by = string)
    M = [Float64(d[(k, p)]) for k in ks, p in ph]
    length(d) == length(M) || throw(ArgumentError("the dictionary has pairs outside the keys and phases given"))
    return pulse_shape_groups(M, ks; kwargs...)
end

function Base.show(io::IO, ::MIME"text/plain", g::PulseShapeGroups)
    @printf(io, "PulseShapeGroups: %d keys in %d shapes (%s), rms %.4f in log multiplier%s",
            length(g.keys), length(g.names), join(("$(c === nothing ? "all" : c) $n" for (c, n) in g.split), ", "),
            g.rms, g.exact ? "" : " (seeded restarts)")
    for (i, name) in enumerate(g.names)
        @printf(io, "\n  %-10s %3d keys, multipliers %s", name, length(g.members[i]),
                join((@sprintf("%.3f", m) for m in g.multipliers[i, :]), " "))
    end
end

Base.show(io::IO, g::PulseShapeGroups) =
    print(io, "PulseShapeGroups(", length(g.keys), " keys, ", length(g.names), " shapes, rms ",
          @sprintf("%.4f", g.rms), ")")

# Every vector s with 1 ≤ s[j] ≤ kmax[j] and sum(s) == total, in lexicographic order.
function _compositions(total::Integer, kmax::Vector{Int})
    out = Vector{Int}[]
    s = zeros(Int, length(kmax))
    function rec(j, left)
        if j == length(kmax)
            1 <= left <= kmax[j] && (s[j] = left; push!(out, copy(s)))
            return
        end
        for v in 1:min(kmax[j], left - (length(kmax) - j))
            s[j] = v
            rec(j + 1, left - v)
        end
    end
    rec(1, total)
    return out
end

# The squared error of a grouping: Σ over groups of Σ |x - mean|².
function _grouping_sse(X::Matrix{Float64}, a::Vector{Int}, k::Int)
    n, d = size(X)
    S, cnt = zeros(k, d), zeros(Int, k)
    for i in 1:n
        cnt[a[i]] += 1
        for j in 1:d
            S[a[i], j] += X[i, j]
        end
    end
    sse = 0.0
    for i in 1:n, j in 1:d
        sse += (X[i, j] - S[a[i], j] / cnt[a[i]])^2
    end
    return sse
end

_centers(X, a, k) = [sum(X[i, j] for i in eachindex(a) if a[i] == g) / count(==(g), a) for g in 1:k, j in 1:size(X, 2)]

# Stirling number of the second kind S(n, k), as Float64 (to compare with a limit).
function _stirling2(n::Int, k::Int)
    S = zeros(Float64, k + 1)
    S[1] = 1.0                      # S(0, 0)
    for m in 1:n
        for j in min(m, k):-1:1
            S[j+1] = j * S[j+1] + S[j]
        end
        S[1] = 0.0
    end
    return S[k+1]
end

# k-means of the rows of X into k groups: (; assignment, centers, sse, exact).
function _kmeans_log(X::Matrix{Float64}, k::Int, seed::UInt64, restarts::Integer, exact_limit::Real)
    n = size(X, 1)
    a = if k == 1
        ones(Int, n)
    elseif k == n
        collect(1:n)
    elseif _stirling2(n, k) <= exact_limit
        _kmeans_enumerate(X, k)
    elseif size(X, 2) == 1
        _kmeans_1d(vec(X), k)
    else
        nothing
    end
    exact = a !== nothing
    a === nothing && (a = _kmeans_restarts(X, k, seed, restarts))
    return (; assignment = a, centers = _centers(X, a, k), sse = _grouping_sse(X, a, k), exact)
end

# Every partition of 1:n into exactly k groups (restricted growth strings), the least
# squared error kept (the first of equal ones).
function _kmeans_enumerate(X::Matrix{Float64}, k::Int)
    n = size(X, 1)
    a = ones(Int, n)
    best, bestsse = copy(a), Inf
    # a[i] ≤ max(a[1:i-1]) + 1, with exactly k groups used
    function rec(i, m)
        if i > n
            m == k || return
            sse = _grouping_sse(X, a, k)
            sse < bestsse && (bestsse = sse; copyto!(best, a))
            return
        end
        # enough keys left to open the groups still missing
        for g in 1:min(m + 1, k)
            m + (g > m) + (n - i) >= k || continue
            a[i] = g
            rec(i + 1, max(m, g))
        end
    end
    a[1] = 1
    rec(2, 1)
    return best
end

# One dimension: an optimal grouping is contiguous in the sorted values; dynamic
# programming over the split points (O(k n²)).
function _kmeans_1d(x::Vector{Float64}, k::Int)
    n = length(x)
    p = sortperm(x)
    y = x[p]
    c1, c2 = [0.0; cumsum(y)], [0.0; cumsum(y .^ 2)]
    cost(i, j) = (s = c1[j+1] - c1[i]; c2[j+1] - c2[i] - s^2 / (j - i + 1))   # y[i:j]
    D = fill(Inf, k, n)
    arg = zeros(Int, k, n)
    for j in 1:n
        D[1, j] = cost(1, j)
    end
    for g in 2:k, j in g:n
        for i in g:j                 # the last group is y[i:j]
            v = D[g-1, i-1] + cost(i, j)
            v < D[g, j] && (D[g, j] = v; arg[g, j] = i)
        end
    end
    a = zeros(Int, n)
    j = n
    for g in k:-1:1
        i = g == 1 ? 1 : arg[g, j]
        a[p[i:j]] .= g
        j = i - 1
    end
    return a
end

# SplitMix64 (Steele, Lea and Flood, OOPSLA 2014): a stream fixed by its seed on every
# Julia version, unlike `Random`'s default generators.
mutable struct _SplitMix
    s::UInt64
end
function _next!(r::_SplitMix)
    r.s += 0x9e3779b97f4a7c15
    z = r.s
    z = (z ⊻ (z >> 30)) * 0xbf58476d1ce4e5b9
    z = (z ⊻ (z >> 27)) * 0x94d049bb133111eb
    return z ⊻ (z >> 31)
end
_uniform!(r::_SplitMix) = (_next!(r) >> 11) * 2.0^-53          # in [0, 1)

# The best of `restarts` k-means++ seedings (Arthur and Vassilvitskii, SODA 2007), each
# assigned to its nearest seed and refined by Hartigan's single-key transfers.
function _kmeans_restarts(X::Matrix{Float64}, k::Int, seed::UInt64, restarts::Integer)
    n, d = size(X)
    rng = _SplitMix(seed)
    best, bestsse = zeros(Int, n), Inf
    dist2(i, c) = sum(abs2, view(X, i, :) .- view(c, :))
    for _ in 1:restarts
        seeds = Int[1 + Int(_next!(rng) % UInt64(n))]
        D = [dist2(i, X[seeds[1], :]) for i in 1:n]
        while length(seeds) < k
            tot = sum(D)
            if tot <= 0             # fewer distinct points than groups: any key not yet a seed
                push!(seeds, first(setdiff(1:n, seeds)))
            else
                u, acc, pick = _uniform!(rng) * tot, 0.0, 0
                for i in 1:n
                    acc += D[i]
                    acc > u && (pick = i; break)
                end
                pick == 0 && (pick = findlast(>(0), D))
                push!(seeds, pick)
            end
            D .= min.(D, [dist2(i, X[seeds[end], :]) for i in 1:n])
        end
        a = [argmin([dist2(i, X[s, :]) for s in seeds]) for i in 1:n]
        for (g, s) in enumerate(seeds)
            a[s] = g                 # each seed in its own group: none empty
        end
        _hartigan!(a, X, k)
        sse = _grouping_sse(X, a, k)
        sse < bestsse && (bestsse = sse; copyto!(best, a))
    end
    return best
end

# Single-key transfers: move a key from group a to group b when it lowers the squared
# error, n_b/(n_b+1) |x - c_b|² < n_a/(n_a-1) |x - c_a|², until no move does.
function _hartigan!(a::Vector{Int}, X::Matrix{Float64}, k::Int)
    n, d = size(X)
    cnt = [count(==(g), a) for g in 1:k]
    S = zeros(k, d)
    for i in 1:n, j in 1:d
        S[a[i], j] += X[i, j]
    end
    moved, sweeps = true, 0
    while moved && sweeps < 1000
        moved, sweeps = false, sweeps + 1
        for i in 1:n
            g = a[i]
            cnt[g] > 1 || continue
            dg = sum(j -> (X[i, j] - S[g, j] / cnt[g])^2, 1:d) * cnt[g] / (cnt[g] - 1)
            bestg, bestd = g, dg
            for h in 1:k
                h == g && continue
                dh = sum(j -> (X[i, j] - S[h, j] / cnt[h])^2, 1:d) * cnt[h] / (cnt[h] + 1)
                dh < bestd * (1 - 1e-12) && ((bestg, bestd) = (h, dh))
            end
            bestg == g && continue
            for j in 1:d
                S[g, j] -= X[i, j]
                S[bestg, j] += X[i, j]
            end
            cnt[g] -= 1
            cnt[bestg] += 1
            a[i] = bestg
            moved = true
        end
    end
    return a
end

"""
    scaled_pulse(base::PulseShape, times, multipliers; peak = nothing) -> PulseShape

The pulse `base` with its power multiplied by a multiplier given per phase: `multipliers[i]`
at `times[i]` (ns, increasing), interpolated linearly between the phase times and held at
the first and last values outside them. The result is tabulated at the union of `base`'s
times and the phase times (inside `base`'s span), each power the base's times the
multiplier there; between those points it is linear, as every `PulseShape` is (the exact
product of two linear pieces is quadratic: the table is the multiplier applied at the
points, not that product).

`peak`, a time after the last phase, leaves the pulse from there on unscaled: the
multiplier then runs linearly from the last phase's value to 1 between the last phase time
and `peak`, and is 1 after it. This is the treatment for a peak at the laser's limit,
which carries no multiplier: the phases sample the pulse at their times only, and the ramp
between the last phase and the peak is an assumption of this function, not a fitted
quantity.
"""
function scaled_pulse(base::PulseShape, times, multipliers; peak = nothing)
    t = collect(Float64, times)
    m = collect(Float64, multipliers)
    length(t) == length(m) || throw(DimensionMismatch("$(length(t)) phase times for $(length(m)) multipliers"))
    isempty(t) && throw(ArgumentError("no phases"))
    issorted(t; lt = <=) || throw(ArgumentError("phase times must increase"))
    if peak !== nothing
        peak > t[end] || throw(ArgumentError("peak ($peak ns) must follow the last phase time ($(t[end]) ns)"))
        push!(t, Float64(peak))
        push!(m, 1.0)
    end
    lo, hi = base.t[1], base.t[end]
    ts = sort!(unique!([base.t; filter(τ -> lo < τ < hi, t)]))
    return PulseShape(ts, [base(τ) * _interp(t, m, τ) for τ in ts])
end

"""
    shape_pulses(groups::PulseShapeGroups, base, times; peak = nothing) -> Dict{Symbol,PulseShape}

The pulse of each shape of `groups` ([`pulse_shape_groups`](@ref)), keyed by the shape's
name: [`scaled_pulse`](@ref) of its class's base pulse with the shape's multipliers at the
phase `times` (and `peak`). `base` is one [`PulseShape`](@ref) for every class, or a
dictionary or function from class to pulse (e.g. the inner and the outer cones' pulses).
[`write_pulse_table`](@ref) writes the result as a table, one column per shape, and
[`with_pulses`](@ref) applies it to beams keyed by `groups.assignment`.
"""
function shape_pulses(g::PulseShapeGroups, base, times; peak = nothing)
    size(g.multipliers, 2) == length(times) ||
        throw(DimensionMismatch("$(length(times)) phase times for $(size(g.multipliers, 2)) phases"))
    pick(c) = base isa PulseShape ? base : base isa AbstractDict ? base[c] : base(c)
    return Dict{Symbol,PulseShape}(name => scaled_pulse(pick(g.class[i]), times, g.multipliers[i, :]; peak)
                                   for (i, name) in enumerate(g.names))
end
