# This file defines the diagnostics of an evaluation (fractions, moments, face flux) and the
# objective factories built on them.

"""
    AbstractKeepOut

Supertype of the smooth geometric rules on the beams that an optimization holds:
[`ApertureKeepOut`](@ref), [`SphereKeepOut`](@ref), [`HardwareKeepOut`](@ref),
[`CaptureKeepIn`](@ref), [`IncidenceRule`](@ref) and [`Omega1KeepOut`](@ref). A subtype
has a [`keepout_clearance`](@ref) method, the clearance per beam, and the fields
`minimum`, the clearance to keep, and `weight`, the weight of its penalty
([`keepout_penalty`](@ref)). [`clearance_floor`](@ref) makes a constraint of any of them.
"""
abstract type AbstractKeepOut end

"Fraction of delivered laser power striking the outside (back side) of any component."
clipped_fraction(ev::Evaluation) = sum(ev.deposition.backside) / ev.deposition.delivered
"Fraction of delivered power that strikes no surface at all."
missed_fraction(ev::Evaluation) = 1 - sum(ev.deposition.per_beam) / ev.deposition.delivered
"Fraction of delivered power absorbed by the named parts or components (see `elements`)."
power_fraction(ev::Evaluation, names::Symbol...) =
    sum(ev.deposition.element[e] for e in elements(ev.scenario.target, names...); init = zero(ev.deposition.delivered)) /
    ev.deposition.delivered

"""
    flux_moments(ev; lmax = 8, axis = ẑ, center = mode_center(target, observed)) -> Vector

Legendre moments `aₗ` of the observed incident flux (index `l+1`), about the center of
the observed capsule by default (`mode_center`).
"""
function flux_moments(ev::Evaluation; lmax::Integer = 8, axis = Vec3(0, 0, 1),
                      center = nothing)
    s = ev.scenario
    _radiative(s) || throw(ArgumentError("scenario has no radiation response"))
    if lmax <= _CACHED_LMAX && Vec3(axis...) == Vec3(0, 0, 1) &&
       (center === nothing || Vec3(center...) == mode_center(s.target, s.observed))
        return _cached_product(getfield(ev, :_moments), s.moment_response, ev.deposition.element)[1:lmax+1]
    end
    return legendre_moments(_full(ev), s.target, s.observed; lmax, axis, center)
end

"""
    mean_flux(ev) -> Real

Mean incident x-ray flux on the observed elements (TW/mm²): the solid-angle mean seen
from the mode center, `flux_moments(ev)[1]`. It is the drive level that the `drive`
keyword of [`symmetry_objective`](@ref) and [`asymmetry_objective`](@ref) holds, so
With `drive = mean_flux(evaluate(scenario, x0))`, the drive of the starting design is
held.
"""
mean_flux(ev::Evaluation) = flux_moments(ev; lmax = 0)[1]

function _full(ev::Evaluation)
    ev.flux === nothing && throw(ArgumentError("scenario has no radiation response"))
    v = zeros(eltype(ev.flux), nelements(ev.scenario.target))
    v[ev.scenario.observed] .= ev.flux
    return v
end

"""
    symmetry_objective(; modes = (2 => 1.0, 4 => 1.0), odd = 1.0, lmax = 6,
                       clip = 10.0, missed = 10.0, direct = 10.0, direct_names = nothing,
                       drive = nothing, drive_weight = 1.0, center = nothing, keepout = nothing,
                       uncertainty = nothing)

Objective for indirect-drive symmetry. The objective is the sum of these terms:

- `Σ wₗ (aₗ/a₀)²` over `modes`;
- `odd ×` the sum of squared odd moments, which measures the top/bottom imbalance;
- `clip`, `missed` and `direct` times the squared fractions of laser power clipped on
  the outside of the target, missing it, or hitting the observed elements directly.
  The observed elements are those named by the scenario `observe`, the capsule whose
  drive is measured. With `direct_names`, the `direct` term counts the named parts or
  components instead, and a name that is not in the target is an error;
- with `drive` (TW/mm²) given, `drive_weight × (a₀/drive - 1)²`, with `a₀` the
  [`mean_flux`](@ref). With `drive = mean_flux(evaluate(scenario, x0))`, the objective holds
  the drive of the design `x0`;
- with a `keepout`, its penalty on the worst-case clearance of the beams. A keep-out is
  an `ApertureKeepOut` for the entrance-hole rims, a `SphereKeepOut` for the capsule, or
  a tuple of them.

The `direct` penalty counts only power already on the capsule, so keep a margin with a
`SphereKeepOut`.

With `uncertainty = (; pointing, power)`, the mode terms are their expected values over
those errors, linearized. The errors are the rms beam-to-beam pointing in mm and the
relative power, as in `mode_uncertainty`. Each `(aₗ/a₀)²` gains the variance due to the
errors, so the optimum also makes the drive shape insensitive to them (robust design).
This costs about three gradients per evaluation. The shared errors of
`mode_uncertainty` (`group_pointing`, `group_power`, `groups`, `positioning`) are
accepted too. They add a fourth per-beam sensitivity, for about four gradients per
evaluation.
"""
function symmetry_objective(; modes = (2 => 1.0, 4 => 1.0), odd = 1.0, lmax::Integer = 6,
                            clip = 10.0, missed = 10.0, direct = 10.0,
                            direct_names = nothing, drive = nothing, drive_weight = 1.0,
                            center = nothing, keepout = nothing, uncertainty = nothing)
    errs = _uncertainty_spec(uncertainty)
    errs === nothing || (center === nothing && lmax <= _CACHED_LMAX) ||
        throw(ArgumentError("uncertainty needs the default center and lmax ≤ $_CACHED_LMAX"))
    return function (ev::Evaluation)
        a = flux_moments(ev; lmax, center)
        J = sum(w * (a[l+1] / a[1])^2 for (l, w) in modes)
        odd > 0 && (J += odd * sum((a[l+1] / a[1])^2 for l in 1:2:lmax))
        if errs !== nothing
            v = _error_variance(ev.scenario, ev.effective, ev.scenario.moment_response, lmax + 1, errs)
            J += sum(w * v[l+1] for (l, w) in modes)
            odd > 0 && (J += odd * sum(v[l+1] for l in 1:2:lmax))
        end
        J += _penalties(ev; clip, missed, direct, direct_names)
        J += keepout_penalty(keepout, ev.beams)
        drive === nothing || (J += drive_weight * (a[1] / drive - 1)^2)
        return J
    end
end

# Returns the errors of the `uncertainty` keyword of an objective, which takes the
# keywords of `mode_uncertainty`, or `nothing`.
_uncertainty_spec(::Nothing) = nothing
function _uncertainty_spec(u)
    e = (; pointing = Float64(get(u, :pointing, 0.0)), power = Float64(get(u, :power, 0.0)),
         group_pointing = Float64(get(u, :group_pointing, 0.0)),
         group_power = Float64(get(u, :group_power, 0.0)), groups = get(u, :groups, :quad),
         positioning = _positioning_rms(get(u, :positioning, 0.0)))
    min(e.pointing, e.power, e.group_pointing, e.group_power) >= 0 ||
        throw(ArgumentError("uncertainty must be nonnegative rms errors"))
    return e
end

function _penalties(ev::Evaluation; clip, missed, direct, direct_names)
    J = zero(ev.deposition.delivered)
    clip > 0 && (J += clip * clipped_fraction(ev)^2)
    missed > 0 && (J += missed * missed_fraction(ev)^2)
    direct > 0 && (J += direct * _direct_fraction(ev, direct_names)^2)
    return J
end

# Returns the fraction of the delivered laser power absorbed by the observed elements, or
# by the named parts or components. A fitted CAD capsule is the part of its pieces.
_direct_fraction(ev::Evaluation, ::Nothing) =
    sum(ev.deposition.element[e] for e in ev.scenario.observed; init = zero(ev.deposition.delivered)) /
    ev.deposition.delivered
function _direct_fraction(ev::Evaluation, names)
    cs = ev.scenario.target.components
    for n in names
        any(c -> c.name === n || c.part === n, cs) ||
            throw(ArgumentError("direct_names: the target has no part or component $(repr(n)) " *
                                "(parts: $(join(repr.(parts(ev.scenario.target)), ", ")))"))
    end
    return power_fraction(ev, names...)
end

"""
    flux_ylm(ev; lmax = 6) -> Vector

Real spherical-harmonic coefficients of the observed incident flux in the flat
`mode_index` order, from the cached `ylm_weights` of the scenario. The polar axis is ẑ,
and the azimuth is measured from +x about the center of the observed capsule
(`mode_center`, which is TCC unless the capsule is offset). Use `ylm_matrix(c, lmax)`
for the `C[l+1, lmax+1+m]` layout.
"""
function flux_ylm(ev::Evaluation; lmax::Integer = _CACHED_YLM_LMAX)
    _radiative(ev.scenario) || throw(ArgumentError("scenario has no radiation response"))
    lmax <= _CACHED_YLM_LMAX || throw(ArgumentError("cached Y_lm weights reach lmax = $_CACHED_YLM_LMAX"))
    return _cached_product(getfield(ev, :_ylm), ev.scenario.ylm_response, ev.deposition.element)[1:(lmax+1)^2]
end

"""
    asymmetry_objective(; modes = (1 => 1.0, 2 => 1.0, 3 => 1.0, 4 => 1.0), components = (),
                        clip = 10.0, missed = 10.0, direct = 10.0, direct_names = nothing,
                        drive = nothing, drive_weight = 1.0, keepout = nothing,
                        reference = nothing, uncertainty = nothing)

Objective for three-dimensional drive symmetry. For each `l => w` in `modes` it adds
`w × relative_mode_power(c, l)` — the rotation-invariant power of all `2l+1` real
harmonics of degree `l` over the mean (`(aₗ/a₀)²/(2l+1)` for an axisymmetric field) — and
for each `(l, m) => w` in `components` adds `w × (c_lm/c_00)²`, so a particular
azimuthal mode (e.g. the `l = 1, m = ±1` drift of a dropped bundle) can be weighted
separately. Penalties for clipping, missing and direct capsule illumination and the
`drive` target and `keepout` are as in `symmetry_objective`, with `drive` compared to
the mean flux `c_00/√(4π)`.

With `reference` (a coefficient vector from `flux_ylm`, e.g. of a tuned design before a
change to the target or the beams), the modes and components measure the departure of
`c_lm/c_00` from `c_lm/c_00` of the reference instead of from zero, so the objective
restores a known-good drive shape.

With `uncertainty = (; pointing, power)` the mode and component terms are their
expected values over independent beam errors of those rms sizes (linearized; see
`mode_uncertainty` and `symmetry_objective`): a robust design. Shared errors
(`group_pointing`, `group_power`, `groups`, `positioning`) are accepted as in
`mode_uncertainty`.
"""
function asymmetry_objective(; modes = (1 => 1.0, 2 => 1.0, 3 => 1.0, 4 => 1.0),
                             components = (), clip = 10.0, missed = 10.0, direct = 10.0,
                             direct_names = nothing, drive = nothing, drive_weight = 1.0,
                             keepout = nothing, reference = nothing, uncertainty = nothing)
    errs = _uncertainty_spec(uncertainty)
    lmax = maximum(vcat([l for (l, _) in modes], [lm[1] for (lm, _) in components], 0))
    nref = (max(lmax, 1) + 1)^2
    reference === nothing || length(reference) >= nref ||
        throw(ArgumentError("reference needs coefficients through l = $(max(lmax, 1))"))
    shape = reference === nothing ? zeros(nref) : reference[1:nref] ./ reference[1]
    return function (ev::Evaluation)
        c = flux_ylm(ev; lmax = max(lmax, 1))
        d = c .- c[1] .* shape                 # departure from the reference shape
        J = zero(c[1])
        for (l, w) in modes
            J += w * sum(abs2, view(d, l^2+1:(l+1)^2)) / c[1]^2   # relative_mode_power
        end
        for ((l, m), w) in components
            J += w * (d[mode_index(l, m)] / c[1])^2
        end
        if errs !== nothing
            v = _error_variance(ev.scenario, ev.effective, ev.scenario.ylm_response, nref, errs)
            for (l, w) in modes
                J += w * sum(@view v[l^2+1:(l+1)^2])
            end
            for ((l, m), w) in components
                J += w * v[mode_index(l, m)]
            end
        end
        J += _penalties(ev; clip, missed, direct, direct_names)
        J += keepout_penalty(keepout, ev.beams)
        drive === nothing || (J += drive_weight * (c[1] / sqrt(4π) / drive - 1)^2)
        return J
    end
end

"""
    influence_matrix(scenario, beams; lmax = 6) -> Matrix

Mode response of the observed flux to each beam: column `b` is the real-harmonic
coefficient vector (`mode_index` order) produced by beam `b` by itself at 1 TW, with the
pointing, spot and focus of the beam. Without a laser–plasma model the flux is linear in
beam powers, so the coefficients of any power assignment `P` are `M * P` — the basis of
power-rebalancing solves and of tables of what losing each beam costs. `beams` are the
effective beams (after `apply_lpi`); a laser–plasma model that deposits light itself
(the ridge stop of `BubbleAbsorption`, linear in the power of each beam) is included
through `lpi_deposit` at the scenario time.
"""
function influence_matrix(s::Scenario, beams::AbstractVector{<:Beam};
                          lmax::Integer = _CACHED_YLM_LMAX)
    _radiative(s) || throw(ArgumentError("scenario has no radiation response"))
    nb = length(beams)
    M = zeros((lmax + 1)^2, nb)
    lmax <= _CACHED_YLM_LMAX || throw(ArgumentError("cached Y_lm weights reach lmax = $_CACHED_YLM_LMAX"))
    Wy = view(s.ylm_response, 1:(lmax+1)^2, :)
    for i in 1:nb
        unit = [remake(b; power = j == i ? 1.0 : 0.0) for (j, b) in enumerate(beams)]
        dep = lpi_deposit(s.lpi, s.deposition, unit, s.time)
        M[:, i] = Wy * dep.element
    end
    return M
end

"""
    irradiance_objective(name; center, radius, intensity = nothing, uniformity = 1.0,
                         spill = 1.0, intensity_weight = 1.0, side = 1)

Objective for a laser-irradiated foil such as a backlighter. The region is the set of
elements of component `name` within `radius` of `center`, on side `side` (1 front, -1
back). The objective is the RMS nonuniformity of the absorbed laser intensity over the
region, plus `spill ×` the squared fraction of delivered power landing outside it. With
`intensity` (TW/mm²) given, the objective adds `intensity_weight × (mean/intensity - 1)²`.
"""
function irradiance_objective(name::Symbol; center, radius, intensity = nothing,
                              uniformity = 1.0, spill = 1.0, intensity_weight = 1.0, side = 1)
    c0 = Vec3(center...)
    return function (ev::Evaluation)
        t = ev.scenario.target
        es = [e for e in elements(t, name) if t.elem_side[e] == side &&
              norm(t.elem_centroid[e] - c0) <= radius]
        isempty(es) && throw(ArgumentError("no elements of $name within radius"))
        P = ev.deposition.element[es]
        A = t.elem_area[es]
        J = uniformity * rms_nonuniformity(P ./ A, A)^2
        spill > 0 && (J += spill * (1 - sum(P) / ev.deposition.delivered)^2)
        if intensity !== nothing
            J += intensity_weight * (sum(P) / sum(A) / intensity - 1)^2
        end
        return J
    end
end

"""
    face_flux(ev) -> Vector

Incident x-ray flux `H_e` (TW/mm²) on each observed element of the scenario: the
full linear flux response (`scenario.response`, computed when first read) times the
absorbed laser power, generic in the number type (the duals of the optimizer). On a
dense model, the same numbers as the `flux` of the evaluation; on a large (iterative) model
it builds the response (a solve per observed element) instead of one iterative solve, so
meant for an observed face of moderate size.
"""
function face_flux(ev::Evaluation)
    s = ev.scenario
    _radiative(s) || throw(ArgumentError("scenario has no radiation response"))
    return _mul(_response(s)::Matrix{Float64}, ev.deposition.element)
end

"""
    face_mean_flux(ev) -> Real

Area mean `H̄ = Σ A_e H_e / Σ A_e` (TW/mm²) of the incident x-ray flux on the observed
elements of the scenario ([`face_flux`](@ref)): the level of [`FaceFlux`](@ref) and the
level that the `flux` keyword of [`face_flux_objective`](@ref) holds, so
With `flux = face_mean_flux(evaluate(scenario, x0))`, the drive of the design `x0` is
held.
"""
face_mean_flux(ev::Evaluation) = _face_terms(ev)[1]

"""
    face_nonuniformity(ev) -> Real

Area-weighted RMS nonuniformity `σ = rms_nonuniformity(H, A)` (relative to the mean) of
the incident x-ray flux on the observed elements of the scenario ([`face_flux`](@ref)):
the shape term of [`face_flux_objective`](@ref), and the scale of its `rms_scale`
keyword, so `rms_scale = face_nonuniformity(evaluate(scenario, x0))` weighs the
nonuniformity at the start as 1.
"""
face_nonuniformity(ev::Evaluation) = _face_terms(ev)[2]

# The area mean and the RMS nonuniformity of the observed face's flux.
function _face_terms(ev::Evaluation)
    H = face_flux(ev)
    A = ev.scenario.target.elem_area[ev.scenario.observed]
    return sum(A .* H) / sum(A), rms_nonuniformity(H, A)
end

"""
    face_flux_objective(; uniformity = 1.0, rms_scale = nothing, flux = nothing,
                        flux_weight = 1.0, clip = 10.0, missed = 10.0, direct = 10.0,
                        direct_names = nothing, keepout = nothing)

Objective for the x-ray drive of a flat observed face, such as a package face in a
halfraum or an opacity sample, which has no capsule modes. The shape term is
`uniformity × (σ/σ₀)²`. Here `σ = rms_nonuniformity(H, A)` is the area-weighted RMS
nonuniformity of the incident flux on the observed elements of the scenario
([`face_nonuniformity`](@ref)). It is computed through the linear flux response of the
scenario, so it is differentiable, and it is the shape term of [`FaceFlux`](@ref) for a
uniform target. `σ₀ = rms_scale`, and `nothing` sets it to 1, which gives the absolute
`σ²`. With `flux` (TW/mm²) given, the objective adds `flux_weight × (H̄/flux - 1)²`,
with `H̄` the [`face_mean_flux`](@ref). The `clip`, `missed`, `direct` and `keepout`
penalties are as in [`symmetry_objective`](@ref). Here `direct` counts laser light on
the face itself, or with `direct_names` on the named parts or components.

`rms_scale = face_nonuniformity(evaluate(scenario, x0))` makes the term equal to
`uniformity` at the design `x0`, so `uniformity = 1` balances the shape against the
penalties on any face. Without it, the term is the absolute `σ²`, which is as small as
the face is uniform. On a uniform 1 mm sample, σ ~ 1e-4 and σ² ~ 1e-8, so the penalties
(10 × a squared light fraction, ~1e-2) set the optimum. The modes of a capsule, by
contrast, have `(Pₗ/P₀)²` ~ 1e-2 before tuning, the same size as the penalties. The
scale is a number fixed at construction, so the objective remains a function of only the
beams, thread-safe and differentiable.
"""
function face_flux_objective(; uniformity = 1.0, rms_scale = nothing, flux = nothing,
                             flux_weight = 1.0, clip = 10.0, missed = 10.0, direct = 10.0,
                             direct_names = nothing, keepout = nothing)
    uniformity >= 0 || throw(ArgumentError("uniformity must not be negative"))
    rms_scale === nothing || (isfinite(rms_scale) && rms_scale > 0) ||
        throw(ArgumentError("rms_scale must be a positive RMS nonuniformity"))
    flux === nothing || flux > 0 || throw(ArgumentError("the flux level must be positive"))
    w = rms_scale === nothing ? Float64(uniformity) : uniformity / Float64(rms_scale)^2
    return function (ev::Evaluation)
        H̄, σ = _face_terms(ev)
        J = w * σ^2
        J += _penalties(ev; clip, missed, direct, direct_names)
        J += keepout_penalty(keepout, ev.beams)
        flux === nothing || (J += flux_weight * (H̄ / flux - 1)^2)
        return J
    end
end
