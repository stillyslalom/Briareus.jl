# Laser–plasma interaction corrections applied between the laser setup and deposition.
#
# A model maps the *setup* beams (what the laser delivers to the chamber: pointing,
# power, wavelength) at a given time to *effective* beams whose powers are what reaches
# and is absorbed in the target. Models must stay generic in the number type so
# gradients flow through power and `dlambda`.

"""
    LaserPlasmaModel

Supertype of laser–plasma corrections. Implement `apply_lpi(model, beams, time)`
returning beams with modified powers (same order, ports and pointing).
"""
abstract type LaserPlasmaModel end

"""
    apply_lpi(model, beams, time = 0.0) -> Vector{Beam}

Effective beams after the laser–plasma `model` (`nothing` = no correction). A tuple
or vector of models is applied in order.
"""
apply_lpi(::Nothing, beams, time = 0.0) = beams
apply_lpi(ms::Union{Tuple,AbstractVector}, beams, time = 0.0) =
    foldl((bs, m) -> apply_lpi(m, bs, time), ms; init = beams)

"""
    lpi_deposit(model, op::DepositionOperator, beams, time = 0.0) -> Deposition

Deposition of the effective `beams` (from [`apply_lpi`](@ref)) under the laser–plasma
`model`: `deposit(op, beams)`, unless the model deposits part of the light where no
beam lands (the ridge stop of [`BubbleAbsorption`](@ref)). Of a tuple or vector of
models, at most one may do so.
"""
lpi_deposit(m, op, beams, time = 0.0) = deposit(op, beams)
function lpi_deposit(ms::Union{Tuple,AbstractVector}, op, beams, time = 0.0)
    ds = [m for m in ms if _deposits(m)]
    length(ds) > 1 && throw(ArgumentError("more than one laser–plasma model deposits light"))
    return isempty(ds) ? deposit(op, beams) : lpi_deposit(only(ds), op, beams, time)
end
_deposits(m) = false

# The same, beam by beam (for the beam-error sensitivities, which deposit each beam
# alone): `nothing`, or `(; g, es, w)`: the operator deposits 1 − g of beam `b`, and
# g of it lands on elements `es` with weights `w` (summing to 1).
_ridge_split(m, op, b, time) = nothing
function _ridge_split(ms::Union{Tuple,AbstractVector}, op, b, time)
    ds = [m for m in ms if _deposits(m)]
    length(ds) > 1 && throw(ArgumentError("more than one laser–plasma model deposits light"))
    return isempty(ds) ? nothing : _ridge_split(only(ds), op, b, time)
end

# A per-group quantity: a number, a function of time, or a (times, values) table.
_tvalue(f::Real, t) = f
_tvalue(f::Function, t) = f(t)
_tvalue(f::Tuple{<:AbstractVector,<:AbstractVector}, t) = _interp(f[1], f[2], t)
_tvalue(f::PulseShape, t) = f(t)

"""
    Backscatter(fractions)

Per-group backscattered power fraction (SRS + SBS), removed before deposition.
`fractions` maps a beam group to a number, a function `t -> fraction`, or a
`(times, values)` table (ns); groups not listed lose nothing. Typically set from
measured FABS/NBI backscatter per cone.

    Backscatter(:inner23 => 0.10, :inner30 => ([0, 5, 8], [0.0, 0.05, 0.12]))
"""
struct Backscatter <: LaserPlasmaModel
    fractions::Dict{Symbol,Any}
end
Backscatter(ps::Pair...) = Backscatter(Dict{Symbol,Any}(Symbol(k) => v for (k, v) in ps))

function apply_lpi(m::Backscatter, beams::AbstractVector{<:Beam}, time = 0.0)
    return [haskey(m.fractions, b.group) ?
            remake(b; power = b.power * (1 - _tvalue(m.fractions[b.group], time))) : b
            for b in beams]
end

"""
    CalibratedTransfer(donors, receivers, fraction)

Empirical cross-beam transfer: a fraction of the total power of the `donors` groups
moves to the `receivers` groups (shared in proportion to each receiving beam's
power). `fraction(Δλ, t)` is a calibration — e.g. fit to shot data — where
`Δλ = ⟨dlambda⟩_receivers − ⟨dlambda⟩_donors` (Å, power-weighted means), so a
positive `Δλ` means the receivers are red-shifted relative to the donors, the
direction that drives transfer toward them. Total power is conserved.

Build `fraction` with `linear_transfer` or `tabulated_transfer`, or pass any
differentiable function.
"""
struct CalibratedTransfer{F} <: LaserPlasmaModel
    donors::Vector{Symbol}
    receivers::Vector{Symbol}
    fraction::F
end
CalibratedTransfer(donors, receivers, f) =
    CalibratedTransfer(Symbol[donors...], Symbol[receivers...], f)

"""
    linear_transfer(f0, slope; range = (0, 1))

`fraction(Δλ, t) = clamp(f0 + slope·Δλ, range...)` (slope per Å).
"""
linear_transfer(f0, slope; range = (0.0, 1.0)) =
    (Δλ, t) -> clamp(f0 + slope * Δλ, range...)

"""
    tabulated_transfer(Δλs, fractions)

Piecewise-linear `fraction(Δλ, t)` through the table, constant beyond its ends.
"""
tabulated_transfer(xs, ys) = (Δλ, t) -> _interp_ad(collect(Float64, xs), collect(Float64, ys), Δλ)

# Linear interpolation that keeps dual numbers in the abscissa.
function _interp_ad(xs, ys, x)
    x <= xs[1] && return ys[1] + 0 * x
    x >= xs[end] && return ys[end] + 0 * x
    k = searchsortedlast(xs, ForwardDiff.value(x))
    w = (x - xs[k]) / (xs[k+1] - xs[k])
    return (1 - w) * ys[k] + w * ys[k+1]
end

function apply_lpi(m::CalibratedTransfer, beams::AbstractVector{<:Beam}, time = 0.0)
    isd = [b.group in m.donors for b in beams]
    isr = [b.group in m.receivers for b in beams]
    any(isd) && any(isr) || return beams
    Pd = sum(b.power for (b, d) in zip(beams, isd) if d)
    Pr = sum(b.power for (b, r) in zip(beams, isr) if r)
    λd = sum(b.power * b.dlambda for (b, d) in zip(beams, isd) if d) / Pd
    λr = sum(b.power * b.dlambda for (b, r) in zip(beams, isr) if r) / Pr
    f = m.fraction(λr - λd, time)
    moved = f * Pd
    return [isd[i] ? remake(b; power = b.power * (1 - f)) :
            isr[i] ? remake(b; power = b.power * (1 + moved / Pr)) : b
            for (i, b) in enumerate(beams)]
end

"""
    lpi_summary(setup, effective) -> Dict{Symbol,NamedTuple}

Per-group setup and effective power (TW) and their ratio.
"""
function lpi_summary(setup::AbstractVector{<:Beam}, eff::AbstractVector{<:Beam})
    out = Dict{Symbol,NamedTuple}()
    for g in unique(b.group for b in setup)
        ps = sum(ForwardDiff.value(b.power) for b in setup if b.group === g)
        pe = sum(ForwardDiff.value(b.power) for b in eff if b.group === g)
        out[g] = (setup = ps, effective = pe, ratio = pe / ps)
    end
    return out
end
