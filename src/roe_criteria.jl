# Built-in criteria of the rules-of-engagement scan (roe.jl). Each constructor wraps an
# existing quantity of an `Evaluation`. The default thresholds are engineering defaults
# (research notes, "Beam rules of engagement"), to be replaced by the limits of the
# campaign.

"""
    roe_drive(; green = 0.02, amber = 0.05, name = :drive, case = :default) -> ROECriterion

Fractional loss of the mean x-ray flux on the observed elements ([`mean_flux`](@ref),
the solid-angle mean about the mode center), for a capsule. The thresholds are
engineering defaults and not facility criteria.
"""
roe_drive(; green::Real = 0.02, amber::Real = 0.05, name = :drive, case::Symbol = :default) =
    ROECriterion(name, mean_flux; green, amber, measure = :loss, case, label = "fraction")

"""
    roe_mode(l; green = 0.01, amber = 0.03, name = Symbol(:P, l), case = :default) -> ROECriterion

Absolute change of the Legendre ratio `aₗ/a₀` of the observed flux
([`flux_moments`](@ref)), such as `P2/P0` and `P4/P0` on a capsule. The thresholds are
engineering defaults and not facility criteria.
"""
function roe_mode(l::Integer; green::Real = 0.01, amber::Real = 0.03, name = Symbol(:P, l),
                  case::Symbol = :default)
    1 <= l <= _CACHED_LMAX || throw(ArgumentError("the mode must lie in 1:$_CACHED_LMAX, got $l"))
    f = ev -> (a = flux_moments(ev; lmax = l); a[l+1] / a[1])
    return ROECriterion(name, f; green, amber, measure = :change, case, label = "of P0")
end

"""
    roe_hotspot_velocity(; green = 30.0, amber = 65.0, measure = :change,
                         name = :hotspot_velocity, case = :default, kwargs...) -> ROECriterion

Hot-spot velocity (km/s) that mode 1 of the capsule flux drives
([`hotspot_velocity`](@ref), which takes `kwargs`). With `measure = :change` the score is
the magnitude of the velocity vector that the outage adds to the velocity of the full
shot, so a nominal mode 1 is not charged to every outage. With `:value` the score is the
speed with the outage, the nominal mode 1 included. The default thresholds are the speeds
at which the piston model loses 2% and 10% of the yield at an implosion velocity of
380 km/s ([`mode1_tolerance`](@ref) gives 29.7 and 67.4 km/s), rounded. They are
engineering defaults and not facility criteria.
"""
function roe_hotspot_velocity(; green::Real = 30.0, amber::Real = 65.0, measure::Symbol = :change,
                              name = :hotspot_velocity, case::Symbol = :default, kwargs...)
    f = ev -> hotspot_velocity(ev; kwargs...)
    return ROECriterion(name, f; green, amber, measure, case, label = "km/s")
end

"""
    roe_yield(; green = 0.02, amber = 0.10, name = :yield, case = :default, kwargs...) -> ROECriterion

Fractional loss of the piston-model yield factor of [`implosion_metrics`](@ref) (which
takes `kwargs`, such as `vimp`), from the hot-spot velocity of mode 1. The thresholds
are engineering defaults and not facility criteria.
"""
function roe_yield(; green::Real = 0.02, amber::Real = 0.10, name = :yield, case::Symbol = :default,
                   kwargs...)
    f = ev -> implosion_metrics(ev; kwargs...).total
    return ROECriterion(name, f; green, amber, measure = :loss, case, label = "fraction")
end

"""
    roe_face_flux(; green = 0.02, amber = 0.05, name = :face_flux, case = :default) -> ROECriterion

Fractional loss of the area-mean x-ray flux on the observed face
([`face_mean_flux`](@ref)), such as a package face or a sample. The thresholds are
engineering defaults and not facility criteria.
"""
roe_face_flux(; green::Real = 0.02, amber::Real = 0.05, name = :face_flux, case::Symbol = :default) =
    ROECriterion(name, face_mean_flux; green, amber, measure = :loss, case, label = "fraction")

"""
    roe_face_nonuniformity(; green = 0.005, amber = 0.02, name = :face_rms,
                           case = :default) -> ROECriterion

Rise of the RMS nonuniformity of the x-ray flux on the observed face
([`face_nonuniformity`](@ref), a fraction of the mean). The thresholds are engineering
defaults and not facility criteria.
"""
roe_face_nonuniformity(; green::Real = 0.005, amber::Real = 0.02, name = :face_rms,
                       case::Symbol = :default) =
    ROECriterion(name, face_nonuniformity; green, amber, measure = :rise, case, label = "of mean")

"""
    roe_temperature(; green = 2.0, amber = 5.0, name = :Tr, case = :default) -> ROECriterion

Change (eV) of the radiation temperature of the area-mean flux on the observed elements,
`radiation_temperature(face_mean_flux(ev))`. With the wall of one cavity observed, it is
the wall temperature of that cavity as [`wall_temperature`](@ref) defines it, at the
albedos of the scenario. A target with several cavities takes one case per cavity. The
thresholds are engineering defaults and not facility criteria.
"""
roe_temperature(; green::Real = 2.0, amber::Real = 5.0, name = :Tr, case::Symbol = :default) =
    ROECriterion(name, ev -> radiation_temperature(face_mean_flux(ev)); green, amber,
                 measure = :change, case, label = "eV")

"""
    roe_flux_balance(a, b; green = 0.02, amber = 0.05, name = :balance, case = :default) -> ROECriterion

Change of the drive balance `(H̄_a − H̄_b) / (H̄_a + H̄_b)` between two sets of observed
elements, with `H̄` the area-mean x-ray flux on the parts or components `a` and `b` (a
name or a collection of names, all observed by the scenario of the case). It scores the
upper against the lower halfraum of a shock tube, or the two sides of a package. The
thresholds are engineering defaults and not facility criteria.
"""
function roe_flux_balance(a, b; green::Real = 0.02, amber::Real = 0.05, name = :balance,
                          case::Symbol = :default)
    na, nb = _roe_names(a), _roe_names(b)
    function f(ev)
        s = ev.scenario
        t = s.target
        H = face_flux(ev)
        pos = Dict(e => i for (i, e) in enumerate(s.observed))
        mean_on(names) = begin
            es = elements(t, names...)
            all(e -> haskey(pos, e), es) && !isempty(es) ||
                throw(ArgumentError("the scenario does not observe every element of $(join(names, ", "))"))
            sum(t.elem_area[e] * H[pos[e]] for e in es) / sum(t.elem_area[e] for e in es)
        end
        Ha, Hb = mean_on(na), mean_on(nb)
        return (Ha - Hb) / (Ha + Hb)
    end
    return ROECriterion(name, f; green, amber, measure = :change, case, label = "balance")
end

_roe_names(n::Symbol) = (n,)
_roe_names(n) = Tuple(Symbol(x) for x in n)

# The elements of a laser-lit region: side `side` of the parts `names`, within `radius`
# of `center` when a center is given.
function _roe_region(t::Target, names, center, radius, side)
    es = [e for e in elements(t, names...) if t.elem_side[e] == side &&
          (center === nothing || norm(t.elem_centroid[e] - center) <= radius)]
    isempty(es) && throw(ArgumentError("no elements of $(join(names, ", ")) in the region"))
    return es
end

_roe_center(c) = c === nothing ? nothing : Vec3(c...)

"""
    roe_foil_irradiance(names; center = nothing, radius = Inf, side = 1, green = 0.05,
                        amber = 0.15, name = :irradiance, case = :default) -> ROECriterion

Fractional loss of the mean absorbed laser intensity on a foil: the power on side `side`
(1 front, -1 back) of the parts or components `names`, within `radius` of `center` when
given, over the area of that region. It applies to a backlighter foil, and to any other
laser-lit surface. The thresholds are engineering defaults and not facility criteria.
"""
function roe_foil_irradiance(names; center = nothing, radius::Real = Inf, side::Integer = 1,
                             green::Real = 0.05, amber::Real = 0.15, name = :irradiance,
                             case::Symbol = :default)
    ns, c0 = _roe_names(names), _roe_center(center)
    function f(ev)
        t = ev.scenario.target
        es = _roe_region(t, ns, c0, radius, side)
        return sum(ev.deposition.element[e] for e in es) / sum(t.elem_area[e] for e in es)
    end
    return ROECriterion(name, f; green, amber, measure = :loss, case, label = "fraction")
end

"""
    roe_foil_nonuniformity(names; center = nothing, radius = Inf, side = 1, green = 0.05,
                           amber = 0.15, name = :irradiance_rms, case = :default) -> ROECriterion

Rise of the area-weighted RMS nonuniformity of the absorbed laser intensity on the foil
region of [`roe_foil_irradiance`](@ref), as a fraction of the mean. The thresholds are
engineering defaults and not facility criteria.
"""
function roe_foil_nonuniformity(names; center = nothing, radius::Real = Inf, side::Integer = 1,
                                green::Real = 0.05, amber::Real = 0.15, name = :irradiance_rms,
                                case::Symbol = :default)
    ns, c0 = _roe_names(names), _roe_center(center)
    function f(ev)
        t = ev.scenario.target
        es = _roe_region(t, ns, c0, radius, side)
        A = t.elem_area[es]
        return rms_nonuniformity(ev.deposition.element[es] ./ A, A)
    end
    return ROECriterion(name, f; green, amber, measure = :rise, case, label = "of mean")
end

"""
    roe_foil_centroid(names; center = nothing, radius = Inf, side = 1, green = 0.02,
                      amber = 0.05, name = :centroid, case = :default) -> ROECriterion

Shift (mm) of the intensity-weighted centroid of the absorbed laser power on the foil
region of [`roe_foil_irradiance`](@ref): the displacement of the x-ray source that a
crystal imager or a pinhole is aligned to. The thresholds are engineering defaults and
not facility criteria.
"""
function roe_foil_centroid(names; center = nothing, radius::Real = Inf, side::Integer = 1,
                           green::Real = 0.02, amber::Real = 0.05, name = :centroid,
                           case::Symbol = :default)
    ns, c0 = _roe_names(names), _roe_center(center)
    function f(ev)
        t = ev.scenario.target
        es = _roe_region(t, ns, c0, radius, side)
        P = sum(ev.deposition.element[e] for e in es)
        return sum(ev.deposition.element[e] * t.elem_centroid[e] for e in es) / P
    end
    return ROECriterion(name, f; green, amber, measure = :change, case, label = "mm")
end

"""
    roe_backlighter_signal(los; line, pulse, distance, green = 0.10, amber = 0.25,
                           name = :backlighter, case = :default, kwargs...) -> ROECriterion

Fractional loss of the He-α photons at the detector, `backlighter_signal(target,
deposition, los; line, pulse, distance, kwargs...).photons` ([`backlighter_signal`](@ref)).
The conversion efficiency depends on the irradiance, so the signal is evaluated in full
for every outage. The thresholds are engineering defaults and not facility criteria.
"""
function roe_backlighter_signal(los; line, pulse::Real, distance::Real, green::Real = 0.10,
                                amber::Real = 0.25, name = :backlighter, case::Symbol = :default,
                                kwargs...)
    f = ev -> backlighter_signal(ev.scenario.target, ev.deposition, los; line, pulse, distance,
                                 kwargs...).photons
    return ROECriterion(name, f; green, amber, measure = :loss, case, label = "fraction")
end

"""
    roe_illumination_rms(target, names = (:capsule,); lmax = 30, green = 0.02, amber = 0.05,
                         measure = :rise, name = :illumination_rms, case = :default) -> ROECriterion

RMS nonuniformity `mode_rms` (1 ≤ l ≤ `lmax`) of the absorbed laser intensity on the
parts `names` of `target`, for hard-sphere direct-drive illumination ([`mode_rms`](@ref)).
The harmonic weights of the elements are built once, at construction. With `measure =
:rise` the score is the increase of the nonuniformity over the full shot, and with
`:value` the nonuniformity with the outage. The thresholds are engineering defaults and
not facility criteria.
"""
function roe_illumination_rms(t::Target, names = (:capsule,); lmax::Integer = 30, green::Real = 0.02,
                              amber::Real = 0.05, measure::Symbol = :rise,
                              name = :illumination_rms, case::Symbol = :default)
    es = elements(t, _roe_names(names)...)
    isempty(es) && throw(ArgumentError("no elements of $(join(_roe_names(names), ", "))"))
    W = ylm_weights(t, es; lmax)
    A = t.elem_area[es]
    f = ev -> mode_rms(W * (ev.deposition.element[es] ./ A), lmax)
    return ROECriterion(name, f; green, amber, measure, case, label = "of mean")
end

"""
    roe_corona_rms(op::CoronaOperator; lmax = op.lmax, green = 0.02, amber = 0.05,
                   measure = :rise, name = :corona_rms, case = :default) -> ROECriterion

RMS nonuniformity `mode_rms` (1 ≤ l ≤ `lmax`) of the laser power absorbed in the corona
of `op` ([`CoronaOperator`](@ref)), for direct drive. The corona deposition is computed
again for every outage from the effective beams. With `measure = :rise` the score is the
increase over the full shot, and with `:value` the nonuniformity with the outage. The
thresholds are engineering defaults and not facility criteria.
"""
function roe_corona_rms(op::CoronaOperator; lmax::Integer = op.lmax, green::Real = 0.02,
                        amber::Real = 0.05, measure::Symbol = :rise, name = :corona_rms,
                        case::Symbol = :default)
    f = ev -> mode_rms(deposit(op, ev.effective).ylm, lmax)
    return ROECriterion(name, f; green, amber, measure, case, label = "of mean")
end
