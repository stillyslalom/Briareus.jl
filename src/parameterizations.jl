# This file builds parameterizations from others. `HeldParameters` holds some parameters
# at given values, and `PowerScales` adds a power scale per group of beams to the
# parameters of the inner parameterization. Both wrap any parameterization with `apply`
# and `initial_parameters`, such as `GroupPointing`, `FramePointing` and other wrappers,
# so they compose.

"""
    HeldParameters(p, held; values = initial_parameters(p))

The parameterization `p` with the parameters `held` fixed at `values` and the rest free. Its
parameter vector contains the free parameters of `p` in their order. `apply(h, y)` is
`apply(p, x)` with `x = values` except at the free indices, where `x = y`. `held` selects
parameters of `p` by

- indices (a vector of integers) or a mask (a `Bool` per parameter);
- names (a vector of `parameter_names(p)` entries, e.g. `"Q31B.z"`);
- a predicate on the name, e.g. `n -> startswith(n, "inner")` for the inner cones of a
  `GroupPointing` grouped by cone, or `n -> first(split(n, '.')) in inner_quads` for a
  set of quad names when grouped by quad.

An unknown name throws an `ArgumentError`. `p` may itself be a composition, such as a
[`PowerScales`](@ref); `p` needs a [`parameter_names`](@ref) method to select by name.

The full parameterization interface follows `p`:

- [`initial_parameters`](@ref) are `values` at the free indices;
- [`parameter_names`](@ref) and [`nparams`](@ref) are those of the free parameters;
- [`bounds`](@ref)`(h, y0; kwargs...)` is the box of `p` about the held values and `y0`,
  restricted to the free parameters, with `kwargs` the keywords of the `bounds` of `p`.

`apply` is generic in the number type. The held values enter as constants of the element
type of `y`, so ForwardDiff duals pass only through the free parameters. Named values pass
through to `p` ([`parameters`](@ref)), and a [`subset`](@ref) of a scenario contains the
names, values and holds of the kept parameters.

```julia
pointing = GroupPointing(quads; by = :quad, fields = (:z, :phi), relative = true)
inner = Set(b.quad for b in quads if b.group in (:inner23, :inner30))
outers_only = HeldParameters(pointing, n -> first(split(n, '.')) in inner)
```
"""
struct HeldParameters{P} <: AbstractParameterization
    inner::P
    values::Vector{Float64}
    free::Vector{Int}
end

function HeldParameters(p, held; values = initial_parameters(p))
    v = collect(Float64, values)
    n = length(v)
    p isa AbstractParameterization && nparams(p) != n &&
        throw(DimensionMismatch("$(length(v)) values for $(nparams(p)) parameters"))
    mask = _held_mask(p, held, n)
    return HeldParameters(p, v, findall(!, mask))
end

# Returns a `Bool` per parameter of `p` (`n` of them), true where `held` selects it.
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
parameter_names(h::HeldParameters) = _parameter_names(h.inner, length(h.values))[h.free]

"""
    initial_parameters(h::HeldParameters) -> Vector{Float64}

Returns the held values at the free indices.
"""
initial_parameters(h::HeldParameters) = h.values[h.free]

# Builds the inner parameter vector in eltype(y) from the held values and `y` at the
# free indices.
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

Returns the beams of the inner parameterization at the held values, with `y` at the free
parameters. The element type follows `y`.
"""
apply(h::HeldParameters, y::AbstractVector) = apply(h.inner, _full_parameters(h, y))

"""
    bounds(h::HeldParameters, y0 = initial_parameters(h); kwargs...) -> (lower, upper)

Returns the box of the inner parameterization, its `bounds` with `kwargs`, about the held
values and `y0`, at the free parameters.
"""
function bounds(h::HeldParameters, y0 = initial_parameters(h); kwargs...)
    lo, hi = bounds(h.inner, _full_parameters(h, collect(Float64, y0)); kwargs...)
    return lo[h.free], hi[h.free]
end

_nparams_shown(h::HeldParameters) = nparams(h)

"""
    PowerScales(p; by = :group)

The parameterization `p` with a power scale per group of beams after its parameters. The
parameter vector is `[x; s]`, and `apply(ps, [x; s])` is `apply(p, x)` with the power of
each beam multiplied by the scale `s` of its group. Groups are keyed by `by` on the beams
of `p` at its initial parameters. `by` is `:group`, which keys cones, or `:quad`, `:beam`
or a function `beam -> key`, such as one with two keys for the inner and outer cones,

```julia
cones = PowerScales(pointing; by = b -> b.group in (:inner23, :inner30) ? :inner : :outer)
```

With this grouping, the cone fraction and the total power vary beside the pointing.
Beams that share a key share a scale, so a grouping ties them. With the shapes of
[`pulse_shape_groups`](@ref), keyed by quad,
`PowerScales(pointing; by = b -> g.assignment[b.quad])` gives one scale per pulse shape
for a re-optimization. The scales start at 1 and are named `"<key>.power_scale"`. They
multiply whatever power `p` gives, including a `:power` field of `p`, so the two are
redundant together. Holding one of them ([`HeldParameters`](@ref)) removes the
redundancy.

[`bounds`](@ref)`(ps, y0; scale = 0.5, scale_min = 0.0, scale_max = Inf, kwargs...)` puts
each scale within `y0 (1 ± scale)`, clipped to `[scale_min, scale_max]`, and the other
parameters in the box of `p` with `kwargs`. The defaults are engineering defaults, and no
facility limit sets them; with `scale_max = 1`, each group stays at or below the power of
`p`. `apply` is generic in the number type. A scale is set by name as the field
`power_scale` of its key ([`parameters`](@ref)), and a [`subset`](@ref) of a scenario
contains the scales of the groups that still have a beam.
"""
struct PowerScales{P,K} <: AbstractParameterization
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

"Number of parameters of a `PowerScales`: those of the inner parameterization, plus one per
group."
nparams(ps::PowerScales) = ps.n + length(ps.keys)

"""
Returns the names of the parameters of a `PowerScales`. These are the names of the inner
parameterization, or `"x[i]"` for one without names, then `"<key>.power_scale"`.
"""
parameter_names(ps::PowerScales) =
    vcat(_parameter_names(ps.inner, ps.n), [string(k, ".power_scale") for k in ps.keys])


"""
    initial_parameters(ps::PowerScales) -> Vector{Float64}

Returns the initial parameters of the inner parameterization, then a scale of 1 per group.
"""
initial_parameters(ps::PowerScales) = vcat(initial_parameters(ps.inner), ones(length(ps.keys)))

"""
    apply(ps::PowerScales, y) -> Vector{Beam}

Returns the beams of the inner parameterization at `y[1:n]`, with the power of each beam
multiplied by the scale of its group in `y[n+1:end]`. The element type follows `y`.
"""
function apply(ps::PowerScales, y::AbstractVector{T}) where {T}
    length(y) == nparams(ps) || throw(DimensionMismatch("expected $(nparams(ps)) parameters"))
    beams = apply(ps.inner, y[1:ps.n])
    length(beams) == length(ps.group) ||
        throw(DimensionMismatch("the inner parameterization gives $(length(beams)) beams, not $(length(ps.group))"))
    return _scaled_powers(beams, ps.group, y, ps.n)
end

# Returns `beams` with each power times the scale `y[n + group[i]]` of its group. This
# function is a barrier, so the beam type of the inner parameterization is known inside it.
_scaled_powers(beams::AbstractVector{<:Beam}, group::Vector{Int}, y::AbstractVector, n::Int) =
    [remake(b; power = b.power * y[n+group[i]]) for (i, b) in enumerate(beams)]

"""
    bounds(ps::PowerScales, y0 = initial_parameters(ps); scale = 0.5, scale_min = 0.0,
           scale_max = Inf, kwargs...) -> (lower, upper)

Returns the box of the inner parameterization, its `bounds` with `kwargs`, for the first
parameters, and each power scale within `y0 (1 ± scale)`, clipped to
`[scale_min, scale_max]`.
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

# ---------------------------------------------------------------------------------------
# Named values and subsets through the wrappers
# ---------------------------------------------------------------------------------------

"""
    parameters(h::HeldParameters, values::Pair...; group = (field = value, ...), ...)
    parameters(ps::PowerScales, values::Pair...; group = (field = value, ...), ...)
        -> Vector{Float64}

Returns the parameter vector of a composed parameterization with named values, given as
for a [`GroupPointing`](@ref). The values are `key => (; field = value, ...)` pairs, or
keywords for `Symbol` keys. The fields pass through to the wrapped `GroupPointing` or
[`FramePointing`](@ref). The scale of a [`PowerScales`](@ref) is the field `power_scale`
of its key, so one pair can set the aim of a quad and its scale:

```julia
p = HeldParameters(PowerScales(GroupPointing(beams; by = :quad, fields = (:z, :phi)); by = :quad),
                   n -> startswith(n, "Q31"))
y0 = parameters(p, "Q15T" => (z = 5.2, power_scale = 0.9))
```

Values not named keep [`initial_parameters`](@ref), which are the values of a
`HeldParameters` at its free parameters and scales of 1. A held parameter, a scale key
with no scale, an unknown group or a field that is not a parameter throws an
`ArgumentError`.
"""
parameters(h::HeldParameters, values::Pair...; kwargs...) = _with_named(h, (values..., pairs(kwargs)...))
parameters(ps::PowerScales, values::Pair...; kwargs...) = _with_named(ps, (values..., pairs(kwargs)...))

# Returns the indices into the parameters of `p` and the values set by the named values
# `kv`, which are `key => (; field = value)` pairs.
function _named_values(p::Union{GroupPointing,FramePointing}, kv)
    out = Pair{Int,Any}[]
    for (key, fv) in kv
        g = findfirst(==(key), p.keys)
        g === nothing && throw(ArgumentError("no group $(repr(key)); the groups are " *
                                             join(repr.(p.keys), ", ")))
        for (f, v) in pairs(fv)
            fi = findfirst(==(f), p.fields)
            fi === nothing && throw(ArgumentError("$f is not a parameter; the fields are " *
                                                  join(p.fields, ", ")))
            push!(out, _pidx(p, g, fi) => v)
        end
    end
    return out
end

# Sets a scale per key from `power_scale` and passes the other fields to the inner
# parameterization.
function _named_values(ps::PowerScales, kv)
    out = Pair{Int,Any}[]
    rest = Pair[]
    for (key, fv) in kv
        others = Pair{Symbol,Any}[]
        for (f, v) in pairs(fv)
            if f === :power_scale
                g = findfirst(isequal(key), ps.keys)
                g === nothing && throw(ArgumentError("no power scale for $(repr(key)); the scales' keys are " *
                                                     join(repr.(ps.keys), ", ")))
                push!(out, ps.n + g => v)
            else
                push!(others, f => v)
            end
        end
        isempty(others) || push!(rest, key => (; others...))
    end
    return append!(out, _named_values(ps.inner, rest))
end

# Maps the named values of the inner parameterization to the free parameters and refuses
# a held one.
function _named_values(h::HeldParameters, kv)
    pos = zeros(Int, length(h.values))
    pos[h.free] = eachindex(h.free)
    out = Pair{Int,Any}[]
    for (i, v) in _named_values(h.inner, kv)
        pos[i] == 0 && throw(ArgumentError("parameter $(_parameter_names(h.inner, length(h.values))[i]) is held " *
                                           "(at $(h.values[i]))"))
        push!(out, pos[i] => v)
    end
    return out
end

_named_values(p, kv) = isempty(kv) ? Pair{Int,Any}[] :
    throw(ArgumentError("named values need a GroupPointing or FramePointing inside, not a $(typeof(p))"))

# Implements the default of `subset` for the wrappers (`subset` in optimize.jl). The
# subset of the inner parameterization, which keeps its groups, is wrapped again. The
# kept parameters, matched by name, keep their values and stay held or free as they were.
function _subset_param(h::HeldParameters, keep)
    inner = _subset_param(h.inner, keep)
    old = _parameter_names(h.inner, length(h.values))
    new = _parameter_names(inner, length(initial_parameters(inner)))
    idx = [something(findfirst(==(n), old), 0) for n in new]
    any(iszero, idx) && throw(ArgumentError("the subset's parameters are not among the parameterization's"))
    isfree = falses(length(h.values))
    isfree[h.free] .= true
    return HeldParameters(inner, h.values[idx], findall(i -> isfree[i], idx))
end

# Keeps the scales of the kept groups, in order of first appearance among the kept beams,
# as the constructor orders them.
function _subset_param(ps::PowerScales, keep)
    inner = _subset_param(ps.inner, keep)
    g = ps.group[keep]
    used = unique(g)
    return PowerScales(inner, length(initial_parameters(inner)), ps.keys[used], [findfirst(==(i), used) for i in g])
end
