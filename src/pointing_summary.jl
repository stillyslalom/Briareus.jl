# A one-line account of a pointing against reference beams: aim shifts, power and cone
# fraction, keep-out clearances and the fractions of laser power clipped and missed.

"""
    PointingSummary

A pointing (e.g. the beams of a re-pointing result) against reference beams, from
[`pointing_summary`](@ref). Fields:

- `nbeams`; `shift_mean` and `shift_max` (mm): the mean and largest aim shift
  `|aim - reference aim|` over the beams, matched by name;
- `power` (TW): the total setup power, and `cone_fraction`: the inner groups' share of it;
- `clearance` and `reference_clearance` (mm): per keep-out, the least finite clearance
  of the beams and of the reference beams (`Inf` with none), and `margin` (mm): the
  clearance less the keep-out's `minimum` (negative: inside the keep-out);
- `clipped` and `missed`: [`clipped_fraction`](@ref) and [`missed_fraction`](@ref) of
  the beams in the scenario.

Its display is one line, suited to a log of a study's runs.
"""
struct PointingSummary
    nbeams::Int
    shift_mean::Float64
    shift_max::Float64
    power::Float64
    cone_fraction::Float64
    clearance::Vector{Float64}
    reference_clearance::Vector{Float64}
    margin::Vector{Float64}
    clipped::Float64
    missed::Float64
end

"""
    pointing_summary(scenario, beams; reference = beams of the scenario's start,
                     keepout = nothing, inner = (:inner23, :inner30)) -> PointingSummary
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

```julia
r = optimize_pointing(scenario, objective)
pointing_summary(scenario, r; keepout = ApertureKeepOut(target))
```
"""
function pointing_summary(s::Scenario, beams::AbstractVector{<:Beam};
                          reference = apply(s.param, initial_parameters(s.param)),
                          keepout = nothing, inner = (:inner23, :inner30))
    ref = Dict(b.name => b for b in reference)
    for b in beams
        haskey(ref, b.name) || throw(ArgumentError("beam $(b.name) has no reference beam"))
    end
    bs = [_value_beam(b) for b in beams]
    shifts = [norm(b.aim - _value_beam(ref[b.name]).aim) for b in bs]
    P = sum(b -> b.power, bs; init = 0.0)
    Pin = sum(b -> b.power, (b for b in bs if b.group in inner); init = 0.0)
    ks = keepout === nothing ? Any[] : keepout isa Union{Tuple,AbstractVector} ? collect(Any, keepout) : Any[keepout]
    least(k, b) = minimum(filter(isfinite, keepout_clearance(k, b)); init = Inf)
    refb = [_value_beam(ref[b.name]) for b in bs]
    c = Float64[least(k, bs) for k in ks]
    cref = Float64[least(k, refb) for k in ks]
    margin = Float64[ci - k.minimum for (ci, k) in zip(c, ks)]
    ev = evaluate(s, bs)
    return PointingSummary(length(bs), isempty(shifts) ? 0.0 : sum(shifts) / length(shifts),
                           maximum(shifts; init = 0.0), P, P > 0 ? Pin / P : NaN, c, cref, margin,
                           Float64(clipped_fraction(ev)), Float64(missed_fraction(ev)))
end

pointing_summary(s::Scenario, r::PointingResult; kwargs...) = pointing_summary(s, r.beams; kwargs...)

function Base.show(io::IO, p::PointingSummary)
    f(x) = @sprintf("%.3f", x)
    print(io, "PointingSummary(", p.nbeams, " beams: aim shifts mean ", f(p.shift_mean),
          " mm, max ", f(p.shift_max), " mm; ", @sprintf("%.1f", p.power), " TW, cone fraction ",
          f(p.cone_fraction))
    for i in eachindex(p.clearance)
        print(io, "; keep-out ", length(p.clearance) > 1 ? "$i " : "", "clearance ", f(p.clearance[i]),
              " mm (reference ", f(p.reference_clearance[i]), ", margin ", f(p.margin[i]), ")")
    end
    print(io, @sprintf("; clipped %.1e, missed %.1e)", p.clipped, p.missed))
end
