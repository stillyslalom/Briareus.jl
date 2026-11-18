# Parameterizations built from others: some parameters held at given values
# (`HeldParameters`), and a power scale per group of beams beside the inner
# parameterization's parameters (`PowerScales`). Both wrap any parameterization with
# `apply` and `initial_parameters` (a `GroupPointing`, a `FramePointing`, or another
# wrapper), so they compose.

"""
    HeldParameters(p, held; values = initial_parameters(p))

The parameterization `p` with the parameters `held` fixed at `values` and the rest free:
its parameter vector holds the free parameters of `p`, in their order, and
`apply(h, y)` is `apply(p, x)` with `x = values` except at the free indices, where
`x = y`. `held` selects parameters of `p` by

- indices (a vector of integers) or a mask (a `Bool` per parameter);
- names (a vector of `parameter_names(p)` entries, e.g. `"Q31B.z"`);
- a predicate on the name, e.g. `n -> startswith(n, "inner")` for the inner cones of a
  `GroupPointing` grouped by cone, or `n -> first(split(n, '.')) in inner_quads` for a
  set of quad names when grouped by quad.

An unknown name is an `ArgumentError`. `p` may itself be a composition, such as a
[`PowerScales`](@ref); `p` needs a [`parameter_names`](@ref) method to select by name.

The full parameterization interface follows `p`: [`initial_parameters`](@ref) are
`values` at the free indices, [`parameter_names`](@ref) and [`nparams`](@ref) those of
the free parameters, and [`bounds`](@ref)`(h, y0; kwargs...)` the box of `p` about the
held values and `y0` (the keywords of `p`'s `bounds`), restricted to the free
parameters. `apply` is generic in the number type: the held values enter as constants
of the element type of `y`, so ForwardDiff duals flow through the free parameters only.

```julia
pointing = GroupPointing(quads; by = :quad, fields = (:z, :phi), relative = true)
inner = Set(b.quad for b in quads if b.group in (:inner23, :inner30))
outers_only = HeldParameters(pointing, n -> first(split(n, '.')) in inner)
```
"""
struct HeldParameters{P}
    inner::P
    values::Vector{Float64}
    free::Vector{Int}
end

function HeldParameters(p, held; values = initial_parameters(p))
    v = collect(Float64, values)
    n = length(v)
    hasmethod(nparams, Tuple{typeof(p)}) && nparams(p) != n &&
        throw(DimensionMismatch("$(length(v)) values for $(nparams(p)) parameters"))
    mask = _held_mask(p, held, n)
    return HeldParameters(p, v, findall(!, mask))
end

# A `Bool` per parameter of `p` (`n` of them): whether `held` selects it.
function _held_mask(p, held::AbstractVector{Bool}, n)
    length(held) == n || throw(DimensionMismatch("a mask of $(length(held)) for $n parameters"))
    return collect(Bool, held)
end
function _held_mask(p, held::AbstractVector{<:Integer}, n)
    all(i -> 1 <= i <= n, held) || throw(ArgumentError("held indices outside 1:$n"))
    m = falses(n)
    m[held] .= true
    return collect(Bool, m)
end
function _held_mask(p, held::AbstractVector{<:AbstractString}, n)
    names = parameter_names(p)
    for h in held
        h in names || throw(ArgumentError("no parameter \"$h\"; the parameters are " * join(names, ", ")))
    end
    return [nm in held for nm in names]
end
_held_mask(p, held::AbstractVector, n) = isempty(held) ? fill(false, n) :
    throw(ArgumentError("held parameters are indices, a mask, names or a predicate on the name"))
_held_mask(p, held::AbstractString, n) = _held_mask(p, [held], n)
_held_mask(p, held::AbstractSet, n) = _held_mask(p, collect(String, held), n)
_held_mask(p, held, n) = Bool[held(nm)::Bool for nm in parameter_names(p)]

"Number of free parameters of a `HeldParameters`."
nparams(h::HeldParameters) = length(h.free)

"Names of the free parameters of a `HeldParameters` (those of its inner parameterization)."
parameter_names(h::HeldParameters) = _inner_names(h.inner, length(h.values))[h.free]

"""
    initial_parameters(h::HeldParameters) -> Vector{Float64}

The held values at the free indices.
"""
initial_parameters(h::HeldParameters) = h.values[h.free]

# The inner parameter vector: the held values, `y` at the free indices, in eltype(y).
function _full_parameters(h::HeldParameters, y::AbstractVector{T}) where {T}
    length(y) == length(h.free) ||
        throw(DimensionMismatch("expected $(length(h.free)) parameters"))
    x = Vector{T}(undef, length(h.values))
    for i in eachindex(h.values)
        x[i] = h.values[i]
    end
    for (j, i) in enumerate(h.free)
        x[i] = y[j]
    end
    return x
end

"""
    apply(h::HeldParameters, y) -> Vector{Beam}

Beams of the inner parameterization at the held values, with `y` at the free
parameters (element type follows `y`).
"""
apply(h::HeldParameters, y::AbstractVector) = apply(h.inner, _full_parameters(h, y))

"""
    bounds(h::HeldParameters, y0 = initial_parameters(h); kwargs...) -> (lower, upper)

The inner parameterization's box (its `bounds` with `kwargs`) about the held values and
`y0`, at the free parameters.
"""
function bounds(h::HeldParameters, y0 = initial_parameters(h); kwargs...)
    lo, hi = bounds(h.inner, _full_parameters(h, collect(Float64, y0)); kwargs...)
    return lo[h.free], hi[h.free]
end

_nparams_shown(h::HeldParameters) = nparams(h)

"""
    PowerScales(p; by = :group)

The parameterization `p` with a power scale per group of beams after its own
parameters: the parameter vector is `[x; s]`, and `apply(ps, [x; s])` is `apply(p, x)`
with each beam's power multiplied by the scale `s` of its group. Groups are keyed by
`by` on the beams of `p` at its initial parameters: `:group` (e.g. the cone), `:quad`,
`:beam`, or a function `beam -> key`, such as two keys for the inner and the outer cones,

```julia
cones = PowerScales(pointing; by = b -> b.group in (:inner23, :inner30) ? :inner : :outer)
```

which varies the cone fraction (and the total power) beside the pointing. The scales
start at 1 and are named `"<key>.power_scale"`. They multiply whatever power `p` gives,
a `:power` field of `p` included, so the two are redundant together; holding one of
them ([`HeldParameters`](@ref)) removes the redundancy.

[`bounds`](@ref)`(ps, y0; scale = 0.5, scale_min = 0.0, scale_max = Inf, kwargs...)`
puts each scale within `y0 (1 ± scale)`, clipped to `[scale_min, scale_max]`
(engineering defaults, not facility limits; `scale_max = 1` holds each group at or below
the power of `p`), and the other parameters in `p`'s box with `kwargs`. `apply` is
generic in the number type.
"""
struct PowerScales{P,K}
    inner::P
    n::Int                 # parameters of `inner`
    keys::Vector{K}        # group keys, in order of first appearance
    group::Vector{Int}     # per beam: its group's index in `keys`
end

function PowerScales(p; by = :group)
    x0 = initial_parameters(p)
    beams = apply(p, x0)
    by === :quad && any(b -> isempty(b.quad), beams) &&
        throw(ArgumentError("by = :quad needs quad names on every beam (OMEGA beams have none)"))
    key = by === :group ? (b -> b.group) : by === :quad ? (b -> b.quad) :
          by === :beam ? (b -> b.name) : by
    ks = map(key, beams)
    uk = unique(ks)
    return PowerScales(p, length(x0), uk, [findfirst(isequal(k), uk) for k in ks])
end

"Number of parameters of a `PowerScales`: the inner parameterization's plus one per group."
nparams(ps::PowerScales) = ps.n + length(ps.keys)

"""
Names of a `PowerScales`' parameters: the inner parameterization's (`"x[i]"` for one
without names), then `"<key>.power_scale"`.
"""
parameter_names(ps::PowerScales) =
    vcat(_inner_names(ps.inner, ps.n), [string(k, ".power_scale") for k in ps.keys])

# The parameter names of `p` (`n` parameters), "x[i]" without a `parameter_names` method.
_inner_names(p, n) = hasmethod(parameter_names, Tuple{typeof(p)}) ?
                     collect(String, parameter_names(p)) : ["x[$i]" for i in 1:n]

"""
    initial_parameters(ps::PowerScales) -> Vector{Float64}

The inner parameterization's initial parameters, then a scale of 1 per group.
"""
initial_parameters(ps::PowerScales) = vcat(initial_parameters(ps.inner), ones(length(ps.keys)))

"""
    apply(ps::PowerScales, y) -> Vector{Beam}

Beams of the inner parameterization at `y[1:n]`, each beam's power multiplied by its
group's scale in `y[n+1:end]` (element type follows `y`).
"""
function apply(ps::PowerScales, y::AbstractVector{T}) where {T}
    length(y) == nparams(ps) || throw(DimensionMismatch("expected $(nparams(ps)) parameters"))
    beams = apply(ps.inner, y[1:ps.n])
    length(beams) == length(ps.group) ||
        throw(DimensionMismatch("the inner parameterization gives $(length(beams)) beams, not $(length(ps.group))"))
    return _scaled_powers(beams, ps.group, y, ps.n)
end

# `beams` with each power times its group's scale `y[n + group[i]]` (a function barrier:
# the inner parameterization's beam type is known here).
_scaled_powers(beams::AbstractVector{<:Beam}, group::Vector{Int}, y::AbstractVector, n::Int) =
    [remake(b; power = b.power * y[n+group[i]]) for (i, b) in enumerate(beams)]

"""
    bounds(ps::PowerScales, y0 = initial_parameters(ps); scale = 0.5, scale_min = 0.0,
           scale_max = Inf, kwargs...) -> (lower, upper)

The inner parameterization's box (its `bounds` with `kwargs`) for the first parameters,
and each power scale within `y0 (1 ± scale)`, clipped to `[scale_min, scale_max]`.
"""
function bounds(ps::PowerScales, y0 = initial_parameters(ps); scale = 0.5, scale_min = 0.0,
                scale_max = Inf, kwargs...)
    length(y0) == nparams(ps) || throw(DimensionMismatch("expected $(nparams(ps)) parameters"))
    lo, hi = bounds(ps.inner, collect(Float64, y0[1:ps.n]); kwargs...)
    s = collect(Float64, y0[ps.n+1:end])
    slo = max.(s .* (1 - scale), scale_min)
    shi = min.(s .* (1 + scale), scale_max)
    all(slo .<= shi) || throw(ArgumentError("empty power-scale bounds: scale_min > scale_max or above the scales"))
    return vcat(lo, slo), vcat(hi, shi)
end

_nparams_shown(ps::PowerScales) = nparams(ps)
