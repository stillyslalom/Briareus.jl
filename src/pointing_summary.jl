# A one-line account of a pointing against reference beams: aim shifts, power and cone
# fraction, keep-out clearances and the fractions of laser power clipped and missed.

"""
    PointingSummary

A pointing (e.g. the beams of a re-pointing result) against reference beams, from
[`pointing_summary`](@ref). Fields:

- `nbeams`; `shift_mean`, `shift_rms` and `shift_max` (mm): the mean, root-mean-square
  and largest aim shift `|aim - reference aim|` over the beams, matched by name (the rms
  is the square root of [`shift_objective`](@ref) with unit weights);
- `power` (TW): the total setup power, and `cone_fraction`: the inner groups' share of it;
- `clearance` and `reference_clearance` (mm): per keep-out, the least clearance of the
  beams and of the reference beams, leaving out `+Inf` (a beam the keep-out does not
  apply to; `Inf` when every beam is such) but keeping `-Inf` (a beam the rule cannot
  place, e.g. a [`CaptureKeepIn`](@ref) beam whose axis misses every capture part) and
  `NaN`, which then are the least and show as such; and `margin` (mm): the clearance
  less the keep-out's `minimum` (negative: inside the keep-out);
- `clipped` and `missed`: [`clipped_fraction`](@ref) and [`missed_fraction`](@ref) of
  the beams in the scenario;
- with a `drive_reference`: `drive_ratio`, the capsule drive's mean over the
  reference's (`c₀₀/r₀₀` of the [`flux_ylm`](@ref) coefficients), and `departure`, the
  [`mode_departure`](@ref) per degree `l = 1:lmax` from the reference's shape, split by
  order into `departure_zero` (`m = 0`, `mode_departure(...; m = :zero)`) and
  `departure_nonzero` (`m ≠ 0`), whose squares add to `departure`'s; `NaN` and empty
  without one;
- with `gates`: `gate`, the gate per degree, and `gate_departure`, the departure it
  bounds (a [`mode_gate`](@ref)'s own value, from its reference's shape, shifted if so;
  otherwise `departure`); with `floor`: `floor`, the pointing floor per degree, against
  which `departure` is shown. Empty without them.

Its display is one line, suited to a log of a study's runs, with the departure over
the gate and over the floor per degree where given (1.00: at the gate).
"""
struct PointingSummary
    nbeams::Int
    shift_mean::Float64
    shift_rms::Float64
    shift_max::Float64
    power::Float64
    cone_fraction::Float64
    clearance::Vector{Float64}
    reference_clearance::Vector{Float64}
    margin::Vector{Float64}
    clipped::Float64
    missed::Float64
    drive_ratio::Float64
    departure::Vector{Float64}
    gate::Vector{Float64}
    gate_departure::Vector{Float64}
    floor::Vector{Float64}
    departure_zero::Vector{Float64}
    departure_nonzero::Vector{Float64}
end

PointingSummary(nbeams, shift_mean, shift_rms, shift_max, power, cone_fraction, clearance,
                reference_clearance, margin, clipped, missed, drive_ratio, departure) =
    PointingSummary(nbeams, shift_mean, shift_rms, shift_max, power, cone_fraction, clearance,
                    reference_clearance, margin, clipped, missed, drive_ratio, departure,
                    Float64[], Float64[], Float64[], Float64[], Float64[])

"""
    pointing_summary(scenario, beams; reference = beams of the scenario's start,
                     keepout = nothing, inner = (:inner23, :inner30),
                     drive_reference = nothing, lmax = 4) -> PointingSummary
    pointing_summary(scenario, result::PointingResult; kwargs...)

Summary of `beams` (or a result's `beams`) in `scenario` against the `reference` beams,
by default the scenario parameterization's beams at its initial parameters (the start
of `optimize_pointing`). Each beam is matched to the reference beam of its name (an
`ArgumentError` when there is none), so a subset's beams compare with the full set's.
`keepout` is an [`ApertureKeepOut`](@ref), a [`SphereKeepOut`](@ref), a
[`HardwareKeepOut`](@ref) (any keep-out with a [`keepout_clearance`](@ref) method), a
tuple or vector of them, or `nothing`; `inner` names the groups whose power is the cone fraction (as for
[`cone_fraction_pulses`](@ref)). The clipped and missed fractions come from
`evaluate(scenario, beams)`, after the scenario's laser–plasma model.

`drive_reference`, an [`Evaluation`](@ref) of a reference design (e.g. the pointing
before a bundle was dropped, in the full scenario) or its [`flux_ylm`](@ref)
coefficients through `l = lmax`, adds the capsule drive against it: the drive ratio
and the departure per degree, the quantities of [`drive_objective`](@ref) and
[`mode_gate`](@ref). It needs a scenario that observes the capsule.

`gates` adds the departure against gates per degree: a vector of gates (in the units of
[`mode_departure`](@ref), against the departure from `drive_reference`), a
[`mode_uncertainty`](@ref) result (its floor, `√(2l+1) × rms_mode`), or a
[`mode_gate`](@ref) constraint, whose own value (the departure from its reference's
shape, shifted if so, with its check of the laser–plasma model and time) and limits are
shown. `floor` (a vector, or a `mode_uncertainty` result) adds the departure from
`drive_reference` against the pointing floor. Both need degrees `1:lmax`.

```julia
r = optimize_pointing(scenario, objective)
pointing_summary(scenario, r; keepout = ApertureKeepOut(target), drive_reference = ref)
```
"""
function pointing_summary(s::Scenario, beams::AbstractVector{<:Beam};
                          reference = apply(s.param, initial_parameters(s.param)),
                          keepout = nothing, inner = (:inner23, :inner30),
                          drive_reference = nothing, lmax::Integer = 4, gates = nothing,
                          floor = nothing)
    ref = Dict(b.name => b for b in reference)
    for b in beams
        haskey(ref, b.name) || throw(ArgumentError("beam $(b.name) has no reference beam"))
    end
    bs = [_value_beam(b) for b in beams]
    shifts = [norm(b.aim - _value_beam(ref[b.name]).aim) for b in bs]
    P = sum(b -> b.power, bs; init = 0.0)
    Pin = sum(b -> b.power, (b for b in bs if b.group in inner); init = 0.0)
    ks = keepout === nothing ? Any[] : keepout isa Union{Tuple,AbstractVector} ? collect(Any, keepout) : Any[keepout]
    # +Inf (a beam the keep-out does not apply to) left out; -Inf (a capture failure) and
    # NaN kept, so that they show (`minimum` propagates NaN)
    least(k, b) = minimum(filter(!=(Inf), keepout_clearance(k, b)); init = Inf)
    refb = [_value_beam(ref[b.name]) for b in bs]
    c = Float64[least(k, bs) for k in ks]
    cref = Float64[least(k, refb) for k in ks]
    margin = Float64[ci - k.minimum for (ci, k) in zip(c, ks)]
    ev = evaluate(s, bs)
    ratio, dep, dz, dn = NaN, Float64[], Float64[], Float64[]
    if drive_reference !== nothing
        cr = drive_reference isa Evaluation ? flux_ylm(drive_reference; lmax) :
             collect(Float64, drive_reference)
        cd = flux_ylm(ev; lmax)
        ratio = cd[1] / cr[1]
        dep = mode_departure(cd, cr; lmax)
        dz, dn = mode_departure(cd, cr; lmax, m = :zero), mode_departure(cd, cr; lmax, m = :nonzero)
    end
    gate, gdep = _summary_gate(gates, ev, dep, lmax)
    fl = floor === nothing ? Float64[] : _summary_floor(floor, dep, lmax)
    return PointingSummary(length(bs), isempty(shifts) ? 0.0 : sum(shifts) / length(shifts),
                           isempty(shifts) ? 0.0 : sqrt(sum(abs2, shifts) / length(shifts)),
                           maximum(shifts; init = 0.0), P, P > 0 ? Pin / P : NaN, c, cref, margin,
                           Float64(clipped_fraction(ev)), Float64(missed_fraction(ev)), ratio, dep,
                           gate, gdep, fl, dz, dn)
end

# The summary's gates per degree and the departure they bound: a mode gate's own value
# and limits, or gates (a vector, a `mode_uncertainty` floor) against `dep`.
_summary_gate(::Nothing, ev, dep, lmax) = (Float64[], Float64[])
function _summary_gate(c::Constraint, ev::Evaluation, dep, lmax)
    v = _gate_value(c.f, ev)
    v === nothing && throw(ArgumentError("gates: constraint $(c.name) is not a mode_gate"))
    g = Float64.(ForwardDiff.value.(v))
    lim = length(c.limit) == length(g) ? copy(c.limit) : fill(c.limit[1], length(g))
    return lim, g
end
function _summary_gate(gates, ev, dep, lmax)
    isempty(dep) && throw(ArgumentError("gates as departures need a `drive_reference` (or give a mode_gate)"))
    g = collect(Float64, _gate_values(gates))
    length(g) == lmax || throw(DimensionMismatch("$(length(g)) gates for degrees 1:$lmax"))
    return g, dep
end

# A mode gate's value at `ev` (`nothing` for another constraint's function).
_gate_value(f, ev) = nothing
_gate_value(g::_GateValue, ev::Evaluation) = g(ev)

function _summary_floor(floor, dep, lmax)
    isempty(dep) && throw(ArgumentError("a floor needs a `drive_reference`"))
    f = collect(Float64, _gate_values(floor))
    length(f) == lmax || throw(DimensionMismatch("a floor of $(length(f)) degrees for 1:$lmax"))
    return f
end

pointing_summary(s::Scenario, r::PointingResult; kwargs...) = pointing_summary(s, r.beams; kwargs...)

"""
    shift_objective(scenario; reference = beams of the scenario's start, weight = b -> 1.0)
    shift_objective(reference::AbstractVector{<:Beam}; weight = b -> 1.0)

An objective for [`optimize_pointing`](@ref): the weighted mean squared aim shift (mm²)
of an evaluation's beams from the `reference` beams, matched by name,
`Σ w_b |aim_b - aim_ref(b)|² / Σ w_b`, by default from the scenario parameterization's
beams at its initial parameters (the start of `optimize_pointing`). With gates as
`constraints` it gives the least re-pointing that meets them, e.g. to cost each
entrance-hole size against a fixed reference drive:

```julia
r = optimize_pointing(scenario, shift_objective(scenario); constraints = (gates, keepout))
sqrt(r.value)            # the rms aim shift (mm), as pointing_summary(scenario, r).shift_rms
```

`weight`, a function of a reference beam (its group, quad or name), is evaluated once on
the reference beams, e.g. `b -> b.group in (:inner23, :inner30) ? 4.0 : 1.0` makes the
inner cones' shifts dearer; the weights must be finite and non-negative with a positive
sum over the beams evaluated. A beam without a reference beam of its name is an
`ArgumentError` at evaluation; reference beams outside the evaluation are ignored (so a
subset's beams compare with the full set's). Only the aim enters, not the defocus. The
objective is generic in the number type and thread-safe (it reads only its tables).
"""
shift_objective(s::Scenario; reference = apply(s.param, initial_parameters(s.param)), weight = b -> 1.0) =
    shift_objective(reference; weight)

function shift_objective(reference::AbstractVector{<:Beam}; weight = b -> 1.0)
    ref = Dict{String,Tuple{SVector{3,Float64},Float64}}()
    for b in reference
        haskey(ref, b.name) && throw(ArgumentError("two reference beams are named $(b.name)"))
        w = Float64(weight(b))
        isfinite(w) && w >= 0 || throw(ArgumentError("beam $(b.name)'s weight $w is not finite and non-negative"))
        ref[b.name] = (SVector{3,Float64}(_value_beam(b).aim), w)
    end
    return _ShiftObjective(ref)
end

struct _ShiftObjective
    ref::Dict{String,Tuple{SVector{3,Float64},Float64}}
end

(o::_ShiftObjective)(ev::Evaluation) = _mean_shift2(o, ev.beams)

function _mean_shift2(o::_ShiftObjective, beams::AbstractVector{Beam{T}}) where {T}
    s, W = zero(T), 0.0
    for b in beams
        r = get(o.ref, b.name, nothing)
        r === nothing && throw(ArgumentError("beam $(b.name) has no reference beam"))
        a, w = r
        w == 0 && continue
        s += w * sum(abs2, b.aim - a)
        W += w
    end
    W > 0 || throw(ArgumentError("the weights of the beams evaluated sum to zero"))
    return s / W
end

function Base.show(io::IO, p::PointingSummary)
    f(x) = @sprintf("%.3f", x)
    print(io, "PointingSummary(", p.nbeams, " beams: aim shifts mean ", f(p.shift_mean),
          " mm, rms ", f(p.shift_rms), " mm, max ", f(p.shift_max), " mm; ", @sprintf("%.1f", p.power), " TW, cone fraction ",
          f(p.cone_fraction))
    for i in eachindex(p.clearance)
        print(io, "; keep-out ", length(p.clearance) > 1 ? "$i " : "", "clearance ", f(p.clearance[i]),
              " mm (reference ", f(p.reference_clearance[i]), ", margin ", f(p.margin[i]), ")")
    end
    print(io, @sprintf("; clipped %.1e, missed %.1e", p.clipped, p.missed))
    isnan(p.drive_ratio) ||
        print(io, @sprintf("; drive %.5f of the reference, departure (l = 1–%d) ", p.drive_ratio,
                           length(p.departure)),
              join((@sprintf("%.4f", 100d) for d in p.departure), " "), " %")
    isempty(p.gate) || print(io, "; over the gates ", join((@sprintf("%.2f", d / g) for (d, g) in zip(p.gate_departure, p.gate)), " "))
    isempty(p.floor) || print(io, "; over the floor ", join((@sprintf("%.2f", d / f) for (d, f) in zip(p.departure, p.floor)), " "))
    print(io, ")")
end
