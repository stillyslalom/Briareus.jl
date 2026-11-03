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
lpi_deposit(m, op, beams, time = 0.0) = _lpi_deposit(m, op, beams, time, nothing)

# `lpi_deposit`, with `split` as in `_deposit` (the light a model deposits off the
# operator added to its beam's group too).
_lpi_deposit(m, op, beams, time, split) = _deposit(op, beams, split)
function _lpi_deposit(ms::Union{Tuple,AbstractVector}, op, beams, time, split)
    ds = [m for m in ms if _deposits(m)]
    length(ds) > 1 && throw(ArgumentError("more than one laser–plasma model deposits light"))
    return isempty(ds) ? _deposit(op, beams, split) : _lpi_deposit(only(ds), op, beams, time, split)
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
measured FABS/NBI backscatter per cone. A fraction lies in [0, 1]: numbers and tables
outside it throw an `ArgumentError` at construction, a function's value outside it
when applied.

    Backscatter(:inner23 => 0.10, :inner30 => ([0, 5, 8], [0.0, 0.05, 0.12]))
"""
struct Backscatter <: LaserPlasmaModel
    fractions::Dict{Symbol,Any}
    function Backscatter(fractions::Dict{Symbol,Any})
        for (g, f) in fractions
            ok = f isa Real ? 0 <= f <= 1 :
                 f isa Tuple{<:AbstractVector,<:AbstractVector} ? all(x -> 0 <= x <= 1, f[2]) : true
            ok || throw(ArgumentError("backscattered fraction of $g outside [0, 1]: $f"))
        end
        return new(fractions)
    end
end
Backscatter(ps::Pair...) = Backscatter(Dict{Symbol,Any}(Symbol(k) => v for (k, v) in ps))
Backscatter(d::AbstractDict) = Backscatter(Dict{Symbol,Any}(Symbol(k) => v for (k, v) in d))

function apply_lpi(m::Backscatter, beams::AbstractVector{<:Beam}, time = 0.0)
    return [haskey(m.fractions, b.group) ?
            remake(b; power = b.power * (1 - _backscattered(m.fractions[b.group], b.group, time))) : b
            for b in beams]
end

function _backscattered(f, g, time)
    x = _tvalue(f, time)
    0 <= x <= 1 || throw(ArgumentError("backscattered fraction of $g at $time ns outside [0, 1]: $x"))
    return x
end

"""
    CalibratedTransfer(donors, receivers, fraction)

Empirical cross-beam transfer: a fraction of the total power of the `donors` groups
moves to the `receivers` groups (shared in proportion to each receiving beam's
power). `fraction(Δλ, t)` is a calibration — e.g. fit to shot data — where
`Δλ = ⟨dlambda⟩_receivers − ⟨dlambda⟩_donors` (Å, power-weighted means), so a
positive `Δλ` means the receivers are red-shifted relative to the donors, the
direction that drives transfer toward them. Total power is conserved.

At zero power the transfer has its limit: donors without power give nothing (`f · 0`),
their `⟨dlambda⟩` then the unweighted mean, so that the derivatives with respect to
their powers are finite; receivers without power take nothing (the light is shared in
proportion to their powers, so there is nothing to share it by), and the beams are
returned unchanged.

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
    ForwardDiff.value(Pr) > 0 || return beams           # no receiving power to share by
    λd = _mean_dlambda(beams, isd, Pd)
    λr = _mean_dlambda(beams, isr, Pr)
    f = m.fraction(λr - λd, time)
    moved = f * Pd
    return [isd[i] ? remake(b; power = b.power * (1 - f)) :
            isr[i] ? remake(b; power = b.power * (1 + moved / Pr)) : b
            for (i, b) in enumerate(beams)]
end

# The power-weighted mean detuning of the beams `sel` (total power `P`); the unweighted
# mean when they have no power (the weights stay in the power's number type).
function _mean_dlambda(beams, sel, P)
    on = ForwardDiff.value(P) > 0
    w(b) = on ? b.power : one(b.power)
    return sum(w(b) * b.dlambda for (b, s) in zip(beams, sel) if s) /
           sum(w(b) for (b, s) in zip(beams, sel) if s)
end

"""
    lpi_summary(setup, effective) -> Dict{Symbol,NamedTuple}

Per-group setup and effective power (TW) and their ratio: 1 for a group without power
in either (no change), `Inf` for one that gains power from none.
"""
function lpi_summary(setup::AbstractVector{<:Beam}, eff::AbstractVector{<:Beam})
    out = Dict{Symbol,NamedTuple}()
    for g in unique(b.group for b in setup)
        ps = sum(ForwardDiff.value(b.power) for b in setup if b.group === g)
        pe = sum((ForwardDiff.value(b.power) for b in eff if b.group === g); init = 0.0)
        out[g] = (setup = ps, effective = pe, ratio = ps > 0 ? pe / ps : pe == 0 ? 1.0 : Inf)
    end
    return out
end

"""
    kemp2026_sbs() -> NamedTuple

Measured inner-cone SBS backscatter of N210808 per quad against time: Kemp et al., PoP 33,
082705 (2026) Fig. 2, digitized (`data/kemp2026_fig2.csv`, provenance in its header).
Columns `t` (ns, 4.4–7.1 every 0.05), `inner23` and `inner30`: all backscattered SBS
power of one 23° or 30° quad (TW; the in-lens signal × 3.4), cone-averaged over both
hemispheres. The paper's energies are 2.8 and 0.9 kJ per quad (eight quads per cone);
single-shot errors reach 100 %, N210808 repeats vary by 30–50 %. SBS starts as the inner
power reaches its peak (the pre-CBET quad power is at half its 9 TW peak at 4.74 ns on
this axis), peaks at 5.3 (23°) and 5.75 ns (30°) and ends by 6.6–6.9 ns, the 30° cone
first (refraction and absorption by an on-axis density spike, the paper's §III C). Used
by `examples/n221204.jl`'s `backscatter=history`. Research notes, "Inner-cone SBS time
histories and the late-time refraction (2026-10-25)".
"""
kemp2026_sbs() = _data_table("kemp2026_fig2.csv")
