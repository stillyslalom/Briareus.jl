# This file defines the x-ray albedo models, the heat front and transmission of thin foils,
# the closing gaps and lips, the capsule radius schedules and the laser conversion.

"""
    AlbedoModel

Supertype of x-ray albedo models: `ConstantAlbedo`, `TabulatedAlbedo`,
`WallLossAlbedo`, `BurnthroughAlbedo`, `GapClosure`, `LipClosure`. Evaluate with `albedo`.
"""
abstract type AlbedoModel end

"""
    ConstantAlbedo(α)

Time-independent x-ray albedo (re-emitted / incident flux).
"""
struct ConstantAlbedo <: AlbedoModel
    value::Float64
end

"""
    TabulatedAlbedo(times, values)

Albedo interpolated linearly in time (ns) and held constant outside the table.
"""
struct TabulatedAlbedo <: AlbedoModel
    t::Vector{Float64}
    value::Vector{Float64}
    function TabulatedAlbedo(t, v)
        length(t) == length(v) || throw(ArgumentError("times and values differ in length"))
        issorted(t) || throw(ArgumentError("times must be sorted"))
        new(collect(Float64, t), collect(Float64, v))
    end
end

"""
    albedo(model, t, T = nothing)
    albedo(material, t = 0.0)

X-ray albedo at time `t` (ns); temperature-dependent models use the radiation
temperature `T` (eV) when given.
"""
albedo(m::ConstantAlbedo, t, T = nothing) = m.value
albedo(m::TabulatedAlbedo, t, T = nothing) = _interp(m.t, m.value, t)

"""
    TemperatureAlbedo(temperatures, values; reference_temperature = 250.0)

Albedo interpolated linearly in the incident-flux radiation temperature (eV) and held
constant outside the table; `reference_temperature` is used where no flux is known.
`solve` iterates it to self-consistency, like `WallLossAlbedo`.
"""
struct TemperatureAlbedo <: AlbedoModel
    T::Vector{Float64}
    value::Vector{Float64}
    reference_temperature::Float64
    function TemperatureAlbedo(T, v; reference_temperature = 250.0)
        length(T) == length(v) || throw(ArgumentError("temperatures and values differ in length"))
        issorted(T) || throw(ArgumentError("temperatures must be sorted"))
        new(collect(Float64, T), collect(Float64, v), Float64(reference_temperature))
    end
end
albedo(m::TemperatureAlbedo, t, T = nothing) =
    _interp(m.T, m.value, (T === nothing || T <= 0) ? m.reference_temperature : T)

"""
    HeatFront(; coefficient = 9.51e-4, temperature_exponent = 1.914, time_exponent = 0.5156)

Areal mass (g/cm²) behind the ablative (subsonic) Marshak heat front that radiation of
temperature `T` drives into a gold wall: at constant `T`,

    m_F = c (T/100 eV)^a (t/1 ns)^b,

with the defaults from the self-similar solution of Hammer & Rosen, Phys. Plasmas 10,
1829 (2003), eq. 112 for a constant drive (q₀ = 0), built on their power-law fits to
the opacity and equation of state of gold. For a drive `T(t)` the front is integrated as

    m_F(t) = c [∫₀ᵗ (T(t′)/100 eV)^(a/b) dt′]^b,

which is exact at constant `T`, depends only on the drive history, and for a ramp
`T ∝ tᵏ` gives `m_F = c (1 + k a/b)^(-b) T(t)^a t^b`: within 2% of Hammer & Rosen's
exact coefficient for `k = 0.182` and within 1% of Shussman & Heizler's (Phys. Plasmas
22, 082109 (2015)) for `k = 0.123` (see `heat_front` and the diagnostic-windows
tutorial). Lindl et al. (Phys. Plasmas 11, 339 (2004), eq. 4-24) give
`c = 1.20e-3, a = 1.86, b = 0.54`, the scaling the `WallLossAlbedo` defaults come
from; it heats about 20% more mass at 300 eV, 1 ns. `T` is the incident-flux
radiation temperature (Lindl's convention; Hammer & Rosen drive with the surface
temperature, a few percent lower for gold late in a pulse).
"""
Base.@kwdef struct HeatFront
    coefficient::Float64 = 9.51e-4
    temperature_exponent::Float64 = 1.914
    time_exponent::Float64 = 0.5156
end

# `_front_rate` is the integrand (T/100 eV)^(a/b) of the drive integral of the front, and
# `_front_mass` is the areal mass for a value of the integral.
_front_rate(f::HeatFront, T) = (max(T, zero(T)) / 100)^(f.temperature_exponent / f.time_exponent)
_front_mass(f::HeatFront, I) = f.coefficient * max(I, zero(I))^f.time_exponent
_front_integral(f::HeatFront, m) = (m / f.coefficient)^(1 / f.time_exponent)

"""
    WallLossAlbedo(; coefficient = 0.32, temperature_exponent = -0.7, time_exponent = -0.38,
                   reference_temperature = 250.0, t0 = 0.1, range = (0.0, 0.95),
                   factor = 1.0, factor_exponent = 0.0, factor_temperature = 300.0,
                   liner = 0.0, front = HeatFront(), release = 1.0)

Marshak-wave wall-loss albedo, `1 - α = c (T/100 eV)^p (t/1 ns)^q`, with defaults for
gold from Lindl's scaling (T in heV, t in ns). `t` is clamped below at `t0`, and α to
`range`. The radiation temperature is the incident-flux temperature of each element, so
`solve` iterates it to self-consistency; `reference_temperature` (eV) is used where
no flux is known. The form is Lindl's for a wall held at `T` since `t = 0`;
[`drive_history`](@ref) evaluates it instead at the effective time of each element, from
the temperature history of the element (`wall_history`).

**Wall material.** `factor` f and `factor_exponent` δ scale the gold loss to that of
another wall material:

    1 - α = r(T) c (T/100 eV)^p (t/1 ns)^q,    r(T) = f (T / factor_temperature)^δ.

The factor `r(T)` is the ratio of the Marshak wall losses of the two materials. Rosen
(2005) gives `E/A ∝ T^3.1 t^0.57` for U₀.₆Dy₀.₂Au₀.₂ and `T^3.3 t^0.6` for Au, a ratio
∝ T^-0.22 t^-0.05. The factor applies to the gold loss after its `range` clamp, so the
albedo of a DU wall may exceed the gold cap. `du_wall_loss()` is the preset for an
Au-lined depleted-uranium wall. A loss above 1 leaves albedo 0, and the wall then absorbs
all it receives; such a loss arises for a lossier material (`r > 1`) where the gold loss is
already large. A `range` reaching below 0, as for the unclamped power law of a closure
study, lowers that floor to the lower end of the range. The constructor throws an
`ArgumentError` for a negative `coefficient`, `factor`, `liner` or `release`, a nonpositive
`t0`, `reference_temperature` or `factor_temperature`, and a `range` with `lo > hi` or
`hi > 1`.

**Liner.** `liner` is the areal mass (g/cm², 0 for none) of a gold liner on the inside of
such a wall. Until the heat front has crossed the liner, the wall loses like gold. The
areal mass of the front is the constant-temperature similarity form of `front` (a
[`HeatFront`](@ref)) at the `T` and `t` of the element, `m_F = c_F (T/100 eV)^a t^b`.
This form does not depend on the temperature history. It is exact at constant `T` and
within a few percent of the drive integral for slow ramps. The ratio becomes

    r_eff = w + (1 - w) r(T),    w = min(1, liner / m_F),

where `w` is the share of the heated mass in the gold liner. This energy-content blend is
an *engineering* choice. It is continuous in `T` and `t`, so the albedo iteration stays
smooth; a sharp front would switch to the loss rate of the substrate as soon as `m_F`
passes the liner. With the defaults (`factor = 1`, no liner) the albedo is identical to
that of gold.

**Re-emission.** `release` scales the stored energy that the wall gives back under a falling
drive in a `drive_history` (`wall_release`). The released energy is the linear (Duhamel)
response of the self-similar loss of the wall, times `release`. The default 1.0 is
*calibrated* against 1-D Lagrangian radiation hydrodynamics of an ablating gold wall in
`examples/late_time_wall.jl`. That model uses the Hammer & Rosen 2003 power laws,
`P = 0.25 e/V`, flux-limited diffusion and a Marshak boundary, and the closure is given the
constant-T scaling of that wall. Drops of 10, 30 and 50% in the drive temperature after a
NIF-like peak need 1.07, 0.99 and 0.91 for the energy the wall takes in over the next
1.5 ns. On the pulses of three NIF shots, 0-D hohlraums need 0.92–1.06 (range
0.9–1.07). That wall reproduces Hammer & Rosen's ablative
self-similar solution (Rosen, OSTI 885380 (2005), Lecture 1) to 1–3% in stored energy
and 3–7% in heated mass at constant and ramped temperature, within the first-order
accuracy of that solution. A static (non-ablating) wall needs 0.77 (the value before
2026-10-25): its cooled surface throttles the release more, while Lindl's loss, which
this model uses, is the ablative one. Research notes, "Wall re-emission and the capsule
radius (2026-10-24)" and "Capsule absorption, conversion and the ablating-wall release
(2026-10-25)".
"""
Base.@kwdef struct WallLossAlbedo <: AlbedoModel
    coefficient::Float64 = 0.32
    temperature_exponent::Float64 = -0.7
    time_exponent::Float64 = -0.38
    reference_temperature::Float64 = 250.0
    t0::Float64 = 0.1
    range::Tuple{Float64,Float64} = (0.0, 0.95)
    factor::Float64 = 1.0
    factor_exponent::Float64 = 0.0
    factor_temperature::Float64 = 300.0
    liner::Float64 = 0.0
    front::HeatFront = HeatFront()
    release::Float64 = 1.0
    function WallLossAlbedo(coefficient, temperature_exponent, time_exponent,
                            reference_temperature, t0, range, factor, factor_exponent,
                            factor_temperature, liner, front, release)
        bad(what) = throw(ArgumentError("WallLossAlbedo: $what"))
        coefficient >= 0 || bad("coefficient must not be negative")
        factor >= 0 || bad("factor must not be negative")
        liner >= 0 || bad("liner must not be negative")
        release >= 0 || bad("release must not be negative")
        t0 > 0 || bad("t0 must be positive")
        reference_temperature > 0 || bad("reference_temperature must be positive")
        factor_temperature > 0 || bad("factor_temperature must be positive")
        range[1] <= range[2] || bad("range must satisfy lo ≤ hi, got $range")
        range[2] <= 1 || bad("an albedo cannot exceed 1: range $range")
        return new(coefficient, temperature_exponent, time_exponent, reference_temperature, t0,
                   range, factor, factor_exponent, factor_temperature, liner, front, release)
    end
end

function albedo(m::WallLossAlbedo, t, T = nothing)
    Te = (T === nothing || T <= 0) ? m.reference_temperature : T
    tc = max(t, m.t0)
    α = clamp(1 - _gold_loss(m, Te, tc), m.range...)
    m.factor == 1 && m.factor_exponent == 0 && return α
    a = 1 - _loss_ratio(m, Te, tc) * (1 - α)
    return max(a, min(m.range[1], zero(a)))     # with a loss above 1 the wall absorbs all
end

# `_gold_loss` returns the loss of gold in Lindl's form, unclamped, and `_loss_ratio`
# returns the ratio `r_eff` of the loss of the material to it. Both take the element
# temperature `Te` and the clamped time `tc`.
_gold_loss(m::WallLossAlbedo, Te, tc) =
    m.coefficient * (Te / 100)^m.temperature_exponent * tc^m.time_exponent
function _loss_ratio(m::WallLossAlbedo, Te, tc)
    r = m.factor * (Te / m.factor_temperature)^m.factor_exponent
    if m.liner > 0
        f = m.front
        w = min(1.0, m.liner / (f.coefficient * (Te / 100)^f.temperature_exponent * tc^f.time_exponent))
        r = w + (1 - w) * r
    end
    return r
end

"""
    wall_loss_factor(model::WallLossAlbedo, t, T)

The ratio `r_eff` of the wall loss of the model to that of gold at time `t` (ns) and
radiation temperature `T` (eV): 1 for gold, `f (T/T_f)^δ` for another material, blended
toward 1 while the heat front is inside a gold `liner`.
"""
function wall_loss_factor(m::WallLossAlbedo, t, T)
    m.factor == 1 && m.factor_exponent == 0 && return 1.0
    Te = (T === nothing || T <= 0) ? m.reference_temperature : T
    return _loss_ratio(m, Te, max(t, m.t0))
end

# The depleted-uranium loss factor at 300 eV. It is calibrated so that an Au-lined DU wall
# is worth 6.5% of peak laser power on a NIF high-foot hohlraum, with the wall loss at its
# effective time. The research notes section "Energetics defaults and the capsule share
# (2026-10-25)" records the calibration.
const DU_LOSS_FACTOR = 0.836

"Areal mass (g/cm²) of the 0.7 µm gold liner of Döppner et al., PRL 115, 055001 (2015)."
const DU_LINER = 19.3 * 0.7e-4

"""
    du_wall_loss(; factor = 0.836, factor_exponent = -0.2, liner = 19.3 × 0.7e-4, kwargs...)

`WallLossAlbedo` for an Au-lined depleted-uranium hohlraum wall: the gold loss times
`r(T) = factor (T/300 eV)^factor_exponent`, gold until the heat front crosses the
`liner` (g/cm²). `kwargs` pass through (the Lindl coefficients for gold by default).

- `factor_exponent` −0.2: the temperature scaling of the cocktail/Au loss ratio, T^-0.22
  from Rosen's Hammer–Rosen fits (M. D. Rosen, "Fundamentals of ICF hohlraums", UCRL-JRNL-
  211678 / OSTI 885380 (2005), §VI) and T^-0.19 between LASNEX's 97.6 / 90.8 / 88.6% at
  185 / 275 / 305 eV (Schein et al., PRL 98, 175003 (2007)).
- `factor` 0.836 at 300 eV is *calibrated* in this model. With it the DU wall is worth 6.5%
  of laser power against gold at the same capsule radiation temperature, the measured
  equivalence of the first DU shots (Döppner 2015, citing Kline 2013 and Callahan 2012).
  Callahan et al., PoP 27, 072704 (2020), report a drive benefit of about 6.5% in both
  simulations and experiments. The calibration case is the NIF high-foot hohlraum of Döppner
  et al. 2015 at peak, with 5.75 × 9.43 mm, 3.1 mm LEHs, a 1.125 mm capsule, 390 TW × 0.86
  coupling, conversion 0.9 and the wall source in inner and outer bands. The wall loss is
  taken at the effective time 2.5 ns, which the wall history of `drive_history` gives a wall
  about 3 ns into a NIF peak. Over effective times 2–4 ns the factor ranges over
  0.827–0.838. At equal power it raises the capsule Tr by 4.8 eV against the measured
  +5 ± 6 eV of Döppner's T-1 pair (N140311 Au 317 eV, N140520 DU 322 eV, 391/393 TW). It
  agrees with the ratio for the U:Dy:Au cocktail, Rosen's "about 0.83" and Schein's 17.5%
  (0.825). Before 2026-10-24 the factor was calibrated at 13 ns since the pulse start with
  the constant-temperature evaluation, and it was 0.785. That wall was too reflective, and
  the loss reduction was worth too little.
- `liner` is the 0.7 µm gold of Döppner 2015's DU hohlraums (7 µm DU, 22.3 µm gold
  backing); HYBRID-E's liner thickness is unpublished. The front crosses it in about
  0.4 ns at 150 eV, so it matters only in the picket.
"""
du_wall_loss(; factor = DU_LOSS_FACTOR, factor_exponent = -0.2, liner = DU_LINER, kwargs...) =
    WallLossAlbedo(; factor, factor_exponent, liner, kwargs...)

"""
    heat_front(front::HeatFront, times, Tr) -> Vector

Areal mass (g/cm²) behind the heat front at each of `times` (ns) for the radiation
temperature history `Tr` (eV) sampled at those times, piecewise linear between them
(trapezoidal drive integral), from zero before `times[1]`.
"""
function heat_front(f::HeatFront, times, Tr)
    length(times) == length(Tr) || throw(ArgumentError("times and temperatures differ in length"))
    I = 0.0
    m = zeros(length(times))
    for k in 2:length(times)
        I += (times[k] - times[k-1]) * (_front_rate(f, Tr[k-1]) + _front_rate(f, Tr[k])) / 2
        m[k] = _front_mass(f, I)
    end
    return m
end

"""
    eddington_transmission(τ)

Fraction of the incident flux that a gray, planar foil of optical depth `τ` transmits
in the Eddington approximation with a Marshak boundary condition,

    F_trans / F_inc = (4/3) / (τ + 2/3) = 4 / (3τ + 2),

eq. 2 of R. S. Lester et al., Phys. Rev. E 114, L023201 (2026), doi:10.1103/k8dz-pjdx.
It is a diffusion result for optically thick foils: it exceeds 1 below `τ = 2/3`,
where it is capped at 1. `BurnthroughAlbedo` takes its default transmission from it
(see `optical_depth`). Generic in the number type.
"""
function eddington_transmission(τ)
    τ >= 0 || throw(DomainError(τ, "the optical depth must not be negative"))
    T = 4 / (3τ + 2)
    return min(one(T), T)
end

# The default foil optical depth is scaled from τ ∼ 5 for a gold foil 2.5 µm thick
# (Lester et al. 2026, after eq. 1), which is 19.3 g/cm³ × 2.5 µm of gold.
const _FOIL_REFERENCE_DEPTH = 5.0
const _FOIL_REFERENCE_AREAL_MASS = 19.3 * 2.5e-4      # g/cm²

_foil_optical_depth(areal_mass) = _FOIL_REFERENCE_DEPTH * areal_mass / _FOIL_REFERENCE_AREAL_MASS

"""
    RosselandTransmission(; reference_temperature = 325.0, exponent = 1.5)

Post-burnthrough transmission of a gold foil that follows the radiation temperature
`T` (eV) of the flux on it, for `BurnthroughAlbedo(; transmission =
RosselandTransmission())`. The optical depth of the foil scales with the Rosseland
opacity of gold,

    τ(T) = τ_ref (T / reference_temperature)^(-exponent),

and the transmission is `eddington_transmission(τ(T))`, where `τ_ref` is the `optical_depth`
of the foil (5 at 2.5 µm). The exponent is that of Hammer & Rosen's power-law fit to the
Rosseland opacity of gold, `1/K = g T^1.5 ρ^-0.2` (Phys. Plasmas 10, 1829 (2003), eq. 31,
fitted over 100–200 eV), so the mean free path grows as `T^1.5`. The reference
temperature is the middle of the 300–350 eV at which Lester et al. (Phys. Rev. E 114,
L023201 (2026)) evaluate `τ ∼ 5`: at the reference the transmission equals the constant
default. Hammer & Rosen's fit also reproduces that point: with the as-built areal mass
of 2.5 µm of gold at Lester's expanded density of 1 g/cm³ it gives `τ = 5.9` at 325 eV
(*derived*, extrapolating the fit). The density is held at that expansion; since
`K ∝ ρ^0.2`, solid gold would double `τ`. Where no flux is known (`T = nothing`) the
reference temperature is used.
"""
Base.@kwdef struct RosselandTransmission
    reference_temperature::Float64 = 325.0
    exponent::Float64 = 1.5
end

# Return the transmitted fraction of a burned-through foil of optical depth τ_ref at
# temperature T.
_transmission(x::Real, τref, Tr) = x
function _transmission(m::RosselandTransmission, τref, Tr)
    Te = (Tr === nothing || Tr <= 0) ? m.reference_temperature : Tr
    return eddington_transmission(τref * (Te / m.reference_temperature)^(-m.exponent))
end

"""
    DiffusionTransmission(; reference_temperature = 325.0, exponent = 1.5)

Post-burnthrough transmission of a gold foil whose Rosseland opacity follows its
temperature through its thickness, for `BurnthroughAlbedo(; transmission =
DiffusionTransmission())` (opt-in; the default is the constant gray estimate). The
optical depth at the temperature `T` (eV) of the flux on the foil is that of
[`RosselandTransmission`](@ref), `τ = τ_ref (T / reference_temperature)^(-exponent)`;
the transmitted fraction `t` is the steady planar diffusion flux through a foil whose
opacity varies as `T^(-α)` (`α = exponent`) across it, its front at `T` and its rear
radiating into vacuum (`T_rear⁴ = t T⁴/2`, the boundaries of Lester et al.'s eq. 2):

    (3τ/4) t = (4/(4 + α)) [1 - (t/2)^((4 + α)/4)].

At `α = 0` this is [`eddington_transmission`](@ref) exactly; for `τ ≫ 1` it tends to
`4/(4 + α)` of it, 8/11 for Hammer & Rosen's `α = 1.5` (Phys. Plasmas 10, 1829 (2003),
eq. 31): the rear of the foil, colder than its front, is more opaque. It is *derived*,
with no parameter fitted to a measurement, and capped at 1 like the gray estimate.

The arrival of the front and the `width` of the transition are unchanged by this choice.
`examples/tamped_foil.jl` models 1-D gold foils with Hammer & Rosen's power laws,
flux-limited diffusion and hydrodynamics. On a public ignition-scale drive, N221204's Dante
temperature at ×1 and ×0.95, those foils transmit over 4–9 ns 0.76–0.79 of the energy this
model gives for 1.5–3 µm of gold. They transmit 0.49–0.54 of the energy the constant
default gives and 0.59–0.62 of that of `RosselandTransmission`. The remainder of the gap
is the later onset of the 1-D foil. The research notes section "A diffusion transmission
for burned-through windows (2026-11-18)" gives the comparison. Generic in the number type.
"""
Base.@kwdef struct DiffusionTransmission
    reference_temperature::Float64 = 325.0
    exponent::Float64 = 1.5
end

function _transmission(m::DiffusionTransmission, τref, Tr)
    Te = (Tr === nothing || Tr <= 0) ? m.reference_temperature : Tr
    return diffusion_transmission(τref * (Te / m.reference_temperature)^(-m.exponent), m.exponent)
end

"""
    diffusion_transmission(τ, α = 1.5)

Steady transmitted fraction of a planar foil of optical depth `τ` at its front, whose
opacity varies as `T^(-α)` through it (front at the drive temperature, rear radiating
into vacuum; see [`DiffusionTransmission`](@ref)): the root of
`(3τ/4) t = (4/(4 + α)) [1 - (t/2)^((4 + α)/4)]`, capped at 1. `α = 0` gives
`eddington_transmission(τ)`. Newton's method from the gray value, which bounds the root
from above (the root is at most the gray value for `α ≥ 0`); generic in the number type.
"""
function diffusion_transmission(τ, α = 1.5)
    τ >= 0 || throw(DomainError(τ, "the optical depth must not be negative"))
    α >= 0 || throw(DomainError(α, "the opacity's temperature exponent must not be negative"))
    k, p = 4 / (4 + α), (4 + α) / 4
    t = 4 / (3τ + 2)                          # gray root, with g(t) ≥ 0 and g convex
    for _ in 1:100
        u = t / 2
        g = 0.75τ * t - k * (1 - u^p)
        tn = t - g / (0.75τ + k * p * u^(p - 1) / 2)
        done = abs(tn - t) <= 4eps(Float64) * abs(t)
        t = tn
        done && break
    end
    return min(one(t), t)
end

"""
    BurnthroughAlbedo(; thickness, density = 19.3, before = WallLossAlbedo(),
                      transmission = eddington_transmission(τ), after = nothing,
                      width = 0.4, burnthrough = nothing, front = HeatFront())

Albedo of a thin wall layer, such as a thin-foil diagnostic window, that re-emits like
the wall until the Marshak heat front (`front`, a `HeatFront`) has heated its whole areal
mass `ρ × thickness`, and then transmits part of the flux that reaches it. The layer
is `thickness` µm of a material of `density` g/cm³, gold by default.

Once burned through, a fraction `f = transmission` of the incident flux passes through
the foil and leaves the hohlraum (a number, `RosselandTransmission()` for a fraction
that follows the radiation temperature, or `DiffusionTransmission()` for one that also
follows the temperature profile of the foil). The rest meets the foil as before and is
re-emitted with the albedo of `before`. Across the transition

    α(t) = α_before (1 - s f),    transmitted fraction s f,

with `s = s((t - t_b)/width)` the smooth step `3x² - 2x³` on [0, 1] and `t_b` the
burnthrough time, so the burned-through albedo is `(1 - f) α_before`. Of the flux
`H A` on the foil, `s f H A` is transmitted and `(1 - s f)(1 - α_before) H A` absorbed
(`energy_balance`, `drive_history`). `before` is an albedo model or a number, by
default the Marshak loss of the wall (`WallLossAlbedo`).

The default `f` is the gray-foil estimate of R. S. Lester et al., Phys. Rev. E 114,
L023201 (2026), doi:10.1103/k8dz-pjdx: `f = 4/(3τ + 2)` (their eq. 2,
`eddington_transmission`). The optical depth is their eq. 1, `τ = L/λ_Au(T_h)` at a
fixed Rosseland mean free path ``λ_{Au}``, scaled linearly in areal mass (`optical_depth`)
from their point `τ ∼ 5` at `L = 2.5 µm`:

    τ = 5 × ρ thickness / (19.3 g/cm³ × 2.5 µm),

so `f ≈ 0.33`, `0.235` and `0.17` for 1.7, 2.5 and 3.6 µm of gold. Limitations:
- ``λ_{Au}`` is held at their evaluation, for a hohlraum radiation temperature `T_h` of
  typically 300–350 eV and an assumed gold density of about 1 g/cm³ (the expanded
  foil), so `f` is constant once burned through.
  With `transmission = RosselandTransmission()`, ``λ_{Au}`` follows the temperature of
  the flux on the foil instead; its density dependence is not modeled.
- Pairing that density with the 2.5 µm as-built thickness is theirs. Taken literally, eq. 1
  at 1 g/cm³ describes a layer of 0.25 mg/cm², against 4.8 mg/cm² in the foil. The scaling
  with thickness (areal mass) is derived in Briareus from eq. 1.
- The estimate is gray, planar and quasi-static, and for gold: pass `transmission`
  for another material. Below `τ = 2/3` (0.33 µm of gold) it is capped at 1.

`transmission = 1` makes a burned-through foil an open hole whose whole incident flux
is transmitted.

An explicit `after` (an albedo model or number; exclusive with `transmission`) sets
the burned-through albedo directly instead,

    α(t) = α_before + s (α_after - α_before),

and counts as transmitted only the loss beyond that of the intact foil,
`(α_before - α) H A`; the `(1 - α_before) H A` of the intact foil stays absorbed.
`after = 0` is the open-hole bound on the loss with this split (the default before
`transmission` existed).

The burnthrough time `t_b` is, in order:
- the per-element time passed to `solve`/`albedos` (`burnthrough = ...`), which
  `drive_history` integrates from the incident-flux temperature
  history of each element (`heat_front`);
- the explicit `burnthrough` (ns), when given;
- otherwise the foil is intact (`t_b = ∞`): a single solve has no history to estimate
  it from. `burnthrough_time` gives the estimate for a temperature history.

`before` is evaluated at `t`; in a standalone `solve` or `albedos`, `t` is the time since
the pulse start, as for a wall without history. In a [`drive_history`](@ref) with
`wall_history` on, a `WallLossAlbedo` `before` instead follows the temperature
history of the element, as the wall elements do: an intact foil equals the wall around
it, and across and after the transition `α_before` follows that history (with the
re-emission of stored energy from the wall scaled by `α/α_before`).

Before `t_b` the foil loses exactly what a thick wall loses and transmits nothing: in
the diffusion picture the cold material ahead of a nonlinear heat front does not
change the flux into it, so no early loss is modeled. Measured burnthrough of 0.8–3 µm
gold foils on a 260 eV hohlraum wall (Nova; Lindl et al., Phys. Plasmas 11, 339
(2004), Fig. 4-6) reaches half signal 0.1–0.25 ns after the computed front arrival,
0.19 ns on average, which sets the default `width` (half transmission at
`t_b + width/2`; see the diagnostic-windows tutorial).

**Caution on an ignition-scale drive.** The onset, the `width` and the constant gray `f`
come from bare Nova foils and the Letter's estimate. `examples/tamped_foil.jl` models a
1-D gold foil with Hammer & Rosen's power laws, flux-limited diffusion and hydrodynamics.
On an ignition-scale drive (N221204's Dante temperature) that foil starts transmitting
about 0.6 ns after the arrival of the front, and its transmission rises slowly. Over
4–9 ns it passes 0.49–0.54 of the energy this default gives it (1.5–3 µm). On Nova's
drive the modeled onset comes 0.2–0.45 ns after the arrival, and the measured half signal
0.09–0.23 ns after it. The one ignition-scale measurement, the window flux
Lester et al. extracted from Dante (its solid-angle calibration pending), corresponds to
about 0.18 of the flux on a 2.5 µm window near 8.5 ns, between 0.16 for the 1-D foil and
0.235 for this default. A study should bracket the transmission and start time of a window
accordingly: `transmission = DiffusionTransmission()` (lower, temperature-following)
against the default, and the arrival of the front against a start up to 0.6 ns later.
Research notes, "Backed windows: what delays the emergent flux, and eq. 3 (2026-11-18)"
and "A diffusion transmission for burned-through windows (2026-11-18)".
"""
struct BurnthroughAlbedo <: AlbedoModel
    thickness::Float64
    density::Float64
    before::AlbedoModel
    after::Union{Nothing,AlbedoModel}       # an explicit burned-through albedo, or
    transmission::Union{Nothing,Float64,RosselandTransmission,DiffusionTransmission}  # the transmitted fraction (exactly one is set)
    width::Float64
    burnthrough::Union{Nothing,Float64}
    front::HeatFront
    function BurnthroughAlbedo(thickness, density, before, after, transmission, width,
                               burnthrough, front)
        thickness > 0 && density > 0 || throw(ArgumentError("thickness and density must be positive"))
        width >= 0 || throw(ArgumentError("the transition width must not be negative"))
        after === nothing || transmission === nothing ||
            throw(ArgumentError("give `after` or `transmission`, not both"))
        if after === nothing
            transmission === nothing &&
                (transmission = eddington_transmission(_foil_optical_depth(density * thickness * 1e-4)))
            transmission isa Union{RosselandTransmission,DiffusionTransmission} || 0 <= transmission <= 1 ||
                throw(ArgumentError("transmission must lie in [0, 1]"))
        end
        new(thickness, density, _albedo_model(before),
            after === nothing ? nothing : _albedo_model(after),
            transmission isa Real ? Float64(transmission) : transmission, width,
            burnthrough === nothing ? nothing : Float64(burnthrough), front)
    end
end
_albedo_model(x) = x isa AlbedoModel ? x : ConstantAlbedo(x)

BurnthroughAlbedo(; thickness, density = 19.3, before = WallLossAlbedo(), after = nothing,
                  transmission = nothing, width = 0.4, burnthrough = nothing,
                  front = HeatFront()) =
    BurnthroughAlbedo(thickness, density, before, after, transmission, width, burnthrough, front)

"Areal mass (g/cm²) of a burnthrough layer."
areal_mass(m::BurnthroughAlbedo) = m.density * m.thickness * 1e-4

"""
    optical_depth(model::BurnthroughAlbedo)

Gray optical depth of a burned-through gold foil, scaled linearly in areal mass from
`τ ∼ 5` at 2.5 µm (Lester et al., Phys. Rev. E 114, L023201 (2026), eq. 1 at a fixed
Rosseland mean free path): `5 × areal_mass(model) / (19.3 g/cm³ × 2.5 µm)`. The
default `transmission` is `eddington_transmission(optical_depth(model))`. An explicit
`transmission` or `after` does not change this estimate.
"""
optical_depth(m::BurnthroughAlbedo) = _foil_optical_depth(areal_mass(m))

"""
    foil_transmission(model::BurnthroughAlbedo, T = nothing)

Fraction of the incident flux that a burned-through foil transmits when the flux on it has
radiation temperature `T` (eV): the `transmission` of the model, evaluated at `T` for a
`RosselandTransmission` or `DiffusionTransmission` (at its reference temperature when
`T` is `nothing`). `NaN`
for a model with an explicit `after`.
"""
foil_transmission(m::BurnthroughAlbedo, T = nothing) =
    m.after === nothing ? _transmission(m.transmission, optical_depth(m), T) : NaN

# Return the transmitting fraction of the transition at time t for burnthrough time tb.
function _burnthrough_step(m::BurnthroughAlbedo, t, tb)
    isfinite(tb) || return 0.0
    m.width == 0 && return t >= tb ? 1.0 : 0.0
    x = clamp((t - tb) / m.width, 0.0, 1.0)
    return x^2 * (3 - 2x)
end

# Return the albedo at time t, the intact (`before`) albedo, and the fraction of the
# incident flux transmitted out of the enclosure.
function _burnthrough_split(m::BurnthroughAlbedo, t, T, burnthrough)
    tb = burnthrough !== nothing && !isnan(burnthrough) ? burnthrough :
         m.burnthrough !== nothing ? m.burnthrough : Inf
    αb = albedo(m.before, t, T)
    s = _burnthrough_step(m, t, tb)
    s == 0 && return αb, αb, zero(αb)
    if m.after === nothing
        τ = s * _transmission(m.transmission, optical_depth(m), T)
        return αb * (1 - τ), αb, τ + zero(αb)
    end
    α = αb + s * (albedo(m.after, t, T) - αb)
    return α, αb, αb - α
end

albedo(m::BurnthroughAlbedo, t, T = nothing; burnthrough = nothing) =
    first(_burnthrough_split(m, t, T, burnthrough))

"""
    GapClosure(; closure, width = 0.4, after = WallLossAlbedo())

Albedo of an opening that closes during the pulse, such as the assembly gap around a
diagnostic window, `WallWindow(; gap, gap_material = Material(:Au; albedo =
GapClosure(closure = 3.0)))`. Until `closure` (ns) the gap is open: all the flux that
reaches it leaves the hohlraum, as through an aperture. Across the next `width` ns it
closes, and then it re-emits like the wall with the albedo of `after` (a model or a
number):

    α(t) = s α_after,    transmitted fraction 1 - s,

with `s = s((t - closure)/width)` the smooth step `3x² - 2x³` on [0, 1], as for
`BurnthroughAlbedo`. `energy_balance` books `(1 - s) H A` as transmitted, and
`drive_history` counts it in `window_loss`. The closure time is an input: a gap tens of
micrometers wide closes as plasma blown off its walls fills it, which this model does
not compute. Unlike an aperture, the gap is opaque to laser light and to view-factor
visibility rays, which matters only for gaps a beam would reach.
"""
struct GapClosure <: AlbedoModel
    closure::Float64
    width::Float64
    after::AlbedoModel
    function GapClosure(closure, width, after)
        width >= 0 || throw(ArgumentError("the closure width must not be negative"))
        new(Float64(closure), Float64(width), _albedo_model(after))
    end
end
GapClosure(; closure, width = 0.4, after = WallLossAlbedo()) = GapClosure(closure, width, after)

# Return the closed fraction of a gap at time t, 0 when open and 1 when closed.
function _closure_step(m::GapClosure, t)
    m.width == 0 && return t >= m.closure ? 1.0 : 0.0
    x = clamp((t - m.closure) / m.width, 0.0, 1.0)
    return x^2 * (3 - 2x)
end

# Return the albedo at time t, the closed albedo, and the transmitted fraction.
function _closure_split(m::GapClosure, t, T)
    αa = albedo(m.after, t, T)
    s = _closure_step(m, t)
    return s * αa, αa, (1 - s) + zero(αa)
end

albedo(m::GapClosure, t, T = nothing) = first(_closure_split(m, t, T))

"""
    LipClosure(; speed = 0.022, exponent = 0.5, depth = 0.24, hold_time = 0.4,
               hold = (6.0, 17.0), after = WallLossAlbedo())

An entrance-hole lip whose closure is computed from the local conditions
(`hohlraum(; leh_closure = :computed)`). The gold blown off the lip fills the hole from
its rim inward; the plasma edge (the radius at which the blow-off turns opaque) sits a
depth `d(t)` inside the initial radius `r₀` and moves as

    dd/dt = (1 − h) v(T) + h (depth − d)/hold_time,    v(T) = speed (T/100 eV)^exponent.

Here `T` is the radiation temperature incident on the lip; the sound speed of x-ray-ablated
gold, `c_s ∝ T^½`, sets the free expansion. The laser hold `h ∈ [0, 1]` is the smooth step
`3x² − 2x³` of the laser intensity through the hole between `hold[1]` and `hold[2]`. That
intensity is the power of the beams entering the hole over its initial area (TW/mm²). While
the beams pass at high power, the blow-off stagnates against the laser-heated plasma in a
dense ridge `depth` inside the rim. Chen et al., PoP 27, 022702 (2020), observe the ridge
forming on the rise to peak where the outer beams cross the lip blow-off, after which it
persists. They attribute the open LEH to the laser heating of the plasma inside the LEH
region. A lip that closed further in a long foot is pushed back to the ridge; Schneider
et al. 2015 observe closure in the low-power part and reopening at peak. Below `hold[1]`, in
the foot and after the laser turns off, the edge moves in freely at `v(T)`.

`drive_history` integrates `d(t)` per lip from `times[1]` (open there) and maps it onto
the lip elements: a lip meshed in rings between its inner radius and `r₀` (the floor
of the closure) has each ring the closed fraction `s` of its area outside the edge, so
the open area is exact at any radial resolution. A closed part re-emits with the albedo
of `after`, `α = s α_after`, and `1 − s` escapes (as for `GapClosure`). A single
`solve` treats the lip as open unless given its closed fraction per element as
`burnthrough`. `lip_depth` integrates the closure for given histories.

The defaults are *calibrated*, and the research notes section "Computed LEH closure
(2026-10-24)" records the fit. `speed`, `depth` and `hold_time` are fitted to the
measured effective LEH radius of five HDC basecamp shots (Chen 2020 Fig. 27, 3.37 mm
LEHs, `chen2020_leh()`), with the lip temperature that the model computes for the
basecamp history. The hold intensities bracket the rise, with free closure at 110 TW and
the hold from about 300 TW of 400. The fit has an rms of 0.033 mm over 31 points, about
the shot-to-shot scatter. `depth` 0.24 mm matches the closure at peak of Young et al.
2020 (SXI, from 1.83 to 1.57 mm, 0.26) and Schneider et al. 2015 (0.84 of 1.55 mm, 0.25). On
3.1–3.64 mm LEHs a common depth and a common fraction of r₀ are indistinguishable, and
the depth is used because the ridge stands where the beams meet the blow-off. After
the laser turns off the measured edge moves in ≤ 0.02 mm/ns over 1.7 ns (LASNEX
0.07–0.11), as the fitted free speed gives at the falling temperature.
"""
struct LipClosure <: AlbedoModel
    speed::Float64
    exponent::Float64
    depth::Float64
    hold_time::Float64
    hold::NTuple{2,Float64}
    after::AlbedoModel
    function LipClosure(speed, exponent, depth, hold_time, hold, after)
        speed >= 0 && depth >= 0 || throw(ArgumentError("speed and depth must not be negative"))
        hold_time > 0 || throw(ArgumentError("hold_time must be positive"))
        0 <= hold[1] < hold[2] || throw(ArgumentError("need 0 ≤ hold[1] < hold[2]"))
        new(Float64(speed), Float64(exponent), Float64(depth), Float64(hold_time),
            (Float64(hold[1]), Float64(hold[2])), _albedo_model(after))
    end
end
LipClosure(; speed = 0.022, exponent = 0.5, depth = 0.24, hold_time = 0.4, hold = (6.0, 17.0),
           after = WallLossAlbedo()) = LipClosure(speed, exponent, depth, hold_time, hold, after)

# Return the laser hold on the lip at intensity I (TW/mm² through the hole), 0 when the
# lip is free and 1 when it is held.
function _lip_hold(m::LipClosure, I)
    x = clamp((I - m.hold[1]) / (m.hold[2] - m.hold[1]), 0.0, 1.0)
    return x^2 * (3 - 2x)
end

# Return the depth after a step Δ from d0, with the mean free speed v and hold h of the
# step. The update is exact for constant coefficients, since the relaxation of the hold is
# stiff on 0.25 ns steps.
function _lip_advance(m::LipClosure, d0, Δ, v, h)
    h <= 0 && return d0 + v * Δ
    b = h / m.hold_time
    deq = ((1 - h) * v + b * m.depth) / b
    return deq + (d0 - deq) * exp(-b * Δ)
end
_lip_speed(m::LipClosure, T) = m.speed * (max(T, zero(T)) / 100)^m.exponent

"""
    lip_depth(model::LipClosure, times, Tr, intensity; depth = 0.0, floor = Inf) -> Vector

Closure depth `d` (mm inside the initial LEH radius) of a `LipClosure` lip at `times`
(ns), given the radiation temperature on the lip `Tr` (eV) and the laser intensity
through the hole `intensity` (TW/mm² over its initial area) at those times: the
trapezoidal mean of `v(T)` and of the hold over each step, integrated exactly for the
relaxation of the hold (as `drive_history` does). The integration starts from `depth` at
`times[1]`, and the depth is capped at `floor` (the radial extent of the lip).
"""
function lip_depth(m::LipClosure, times, Tr, intensity; depth::Real = 0.0, floor::Real = Inf)
    n = length(times)
    length(Tr) == n && length(intensity) == n ||
        throw(ArgumentError("times, temperatures and intensities differ in length"))
    d = zeros(n)
    d[1] = min(depth, floor)
    for k in 2:n
        v = (_lip_speed(m, Tr[k-1]) + _lip_speed(m, Tr[k])) / 2
        h = (_lip_hold(m, intensity[k-1]) + _lip_hold(m, intensity[k])) / 2
        d[k] = clamp(_lip_advance(m, d[k-1], times[k] - times[k-1], v, h), 0.0, floor)
    end
    return d
end

# Return the closed fraction of a lip ring [ra, rb], with radii from the hole center, when
# the plasma edge is at radius r. The fraction is the share of the ring area outside r.
_ring_closed(ra, rb, r) = r <= ra ? 1.0 : r >= rb ? 0.0 : (rb^2 - r^2) / (rb^2 - ra^2)

function _lip_split(m::LipClosure, t, T, s)
    αa = albedo(m.after, t, T)
    s = s isa Real && !isnan(s) ? clamp(s, 0.0, 1.0) : 0.0
    return s * αa, αa, (1 - s) + zero(αa)
end
albedo(m::LipClosure, t, T = nothing) = first(_lip_split(m, t, T, 0.0))

"""
    RadiusSchedule

A capsule radius history, `R(t)/R₀` against time (ns), for `CapsuleAlbedo`'s `radius`
and `moving_hohlraum`'s `capsule_shrink`; [`radius_fraction`](@ref)`(s, t)` evaluates it.
Concrete schedules: [`TabulatedRadius`](@ref) (a table, linear, held at its ends),
[`ScaledRadius`](@ref) (another schedule with its time and radius axes scaled) and
[`TrajectoryRadius`](@ref) (a kinematic trajectory through published anchors). A new kind
(a rocket model, say) subtypes this and defines `radius_fraction`. Both
keywords also take a `(times, R/R₀)` table, which becomes a `TabulatedRadius`, or a
function `t -> R/R₀`.

Which radius a published history follows (ablation front, ablator–fuel interface, fuel
center of mass) and which shot it belongs to matters: research notes, "Capsule radius
schedules: published histories and the API (2026-10-26)".
"""
abstract type RadiusSchedule end

"""
    TabulatedRadius(t, r; source = "")

A [`RadiusSchedule`](@ref) from a table: `r` = R/R₀ at the times `t` (ns, strictly
increasing), linear between them and held at the ends. `source` records where the table
comes from (paper and figure, and whether it is published or derived).
[`kritcher2024_radius`](@ref Briareus.Benchmarks.kritcher2024_radius) gives the published N221204 and N210808 tables.
"""
struct TabulatedRadius <: RadiusSchedule
    t::Vector{Float64}
    r::Vector{Float64}
    source::String
end
function TabulatedRadius(t, r; source::AbstractString = "")
    length(t) == length(r) && !isempty(t) ||
        throw(ArgumentError("a radius table needs as many times as radii, at least one"))
    all(>(0), diff(t)) || throw(ArgumentError("radius table times must increase strictly"))
    all(x -> isfinite(x) && x >= 0, r) || throw(ArgumentError("R/R₀ must be finite and ≥ 0"))
    return TabulatedRadius(collect(Float64, t), collect(Float64, r), String(source))
end

"""
    ScaledRadius(base; time_scale = 1.0, radius_scale = 1.0, source = "")

The [`RadiusSchedule`](@ref) `base` stretched in time and radius:

    R(t)/R₀ = radius_scale · base(t / time_scale),

so a feature of `base` at time t₁ (its minimum, say) falls at `time_scale · t₁`. For a
[`TabulatedRadius`](@ref) base the table times are scaled once (`base.t .* time_scale`,
interpolated as a separate table). A scaled schedule is *derived*: it stands in for a
shot whose own history is unpublished (e.g. Zylstra's LEH pair on Kritcher 2024's
N210808 capsule, its minimum at 9.18 ns moved to the bang time of each shot). `source`
defaults to the `source` of the base with the factors appended.
"""
struct ScaledRadius{S} <: RadiusSchedule
    base::S
    time_scale::Float64
    radius_scale::Float64
    source::String
    times::Vector{Float64}   # a tabulated base's scaled times (empty otherwise)
end
function ScaledRadius(base; time_scale::Real = 1.0, radius_scale::Real = 1.0,
                      source::AbstractString = "")
    time_scale > 0 && radius_scale > 0 ||
        throw(ArgumentError("time_scale and radius_scale must be positive"))
    b = _radius_schedule(base)
    b isa Union{RadiusSchedule,Function} ||
        throw(ArgumentError("ScaledRadius needs a RadiusSchedule, a (times, R/R₀) table or a function"))
    src = isempty(source) ?
          string(b isa RadiusSchedule ? b.source : "function", " (derived: time × ", time_scale,
                 ", R × ", radius_scale, ")") : String(source)
    times = b isa TabulatedRadius ? b.t .* Float64(time_scale) : Float64[]
    return ScaledRadius(b, Float64(time_scale), Float64(radius_scale), src, times)
end

"""
    TrajectoryRadius(; start, peak_time, peak_radius, velocity, bang, minimum, radius,
                     source = "")

A kinematic [`RadiusSchedule`](@ref) through the published anchors of an implosion, for a
shot without a published radius history. The surface stays at its initial radius `radius`
(mm) until `start` (ns), when the first shock reaches it, as at the shock-1 breakout of a
keyhole from the ablator. It then falls along a cubic Hermite segment to `peak_radius` (mm)
at `peak_time`, at rest at `start` and with the inward speed `velocity` (km/s) at
`peak_time`. A second segment takes it to `minimum` (R/R₀) at `bang`, at rest there, and
the radius stays at that value after bang; the rebound is not modeled.
`R/R₀` = R/`radius`. The anchors are a keyhole breakout, the bang time and the values that
HYDRA tables give: the peak implosion velocity, its time and the limb-minimum radius at
peak velocity. Kritcher et al., PRE 106, 025201 (2022), list these values in Table I.
Both segments must stay monotone, with the speed at most 3× the segment mean; otherwise
the constructor throws an `ArgumentError`.

The schedule is *derived*, an interpolant that does not model ablation physics. It is
checked against Kritcher 2024 Fig. 1, which was not used to build it. N221204's anchors
are a keyhole start at 3.864 ns, bang at 9.53 ns and 383 km/s, with peak velocity
0.44 ns before bang at 197 µm; these two values come from the N210207/N210307/N210808
HYDRA rows. With these anchors, R/R₀ lies within 0.019 of the published history (0.029
in R²) from 3 ns to peak velocity. Research notes, "Capsule radius schedules:
published histories and the API (2026-10-26)".
"""
struct TrajectoryRadius <: RadiusSchedule
    start::Float64
    peak_time::Float64
    peak_fraction::Float64
    peak_speed::Float64      # R₀ per ns
    bang::Float64
    minimum::Float64
    source::String
end
function TrajectoryRadius(; start::Real, peak_time::Real, peak_radius::Real, velocity::Real,
                          bang::Real, minimum::Real, radius::Real, source::AbstractString = "")
    start < peak_time < bang || throw(ArgumentError("need start < peak_time < bang"))
    radius > 0 && velocity > 0 || throw(ArgumentError("radius and velocity must be positive"))
    rv, v = peak_radius / radius, velocity / 1000 / radius
    1 > rv > minimum >= 0 || throw(ArgumentError("need radius > peak_radius > minimum · radius ≥ 0"))
    # A cubic Hermite segment with end slopes 0 and β times the mean slope stays monotone
    # for β ≤ 3 (Fritsch & Carlson 1980).
    v <= 3 * (1 - rv) / (peak_time - start) && v <= 3 * (rv - minimum) / (bang - peak_time) ||
        throw(ArgumentError("velocity too high for a monotone trajectory through these anchors"))
    return TrajectoryRadius(Float64(start), Float64(peak_time), Float64(rv), Float64(v),
                            Float64(bang), Float64(minimum), String(source))
end
function radius_fraction(s::TrajectoryRadius, t)
    t <= s.start && return one(t) * 1.0
    t >= s.bang && return one(t) * s.minimum
    t <= s.peak_time && return _hermite(t, s.start, 1.0, 0.0, s.peak_time, s.peak_fraction, -s.peak_speed)
    return _hermite(t, s.peak_time, s.peak_fraction, -s.peak_speed, s.bang, s.minimum, 0.0)
end

"""
    radius_fraction(s, t) -> R(t)/R₀

The capsule radius over its initial radius at time `t` (ns) on the schedule `s`: a
[`RadiusSchedule`](@ref), a `(times, R/R₀)` table, a function `t -> R/R₀`, or `nothing`
(R = R₀, 1).
"""
radius_fraction(s::TabulatedRadius, t) = _interp(s.t, s.r, t)
radius_fraction(s::ScaledRadius{TabulatedRadius}, t) = s.radius_scale * _interp(s.times, s.base.r, t)
radius_fraction(s::ScaledRadius, t) = s.radius_scale * radius_fraction(s.base, t / s.time_scale)
radius_fraction(s, t) = _radius_fraction(_radius_schedule(s), t)

# Convert a `(times, R/R₀)` table to a `TabulatedRadius`, and return anything else as given.
_radius_schedule(r::Tuple) = TabulatedRadius(r[1], r[2])
_radius_schedule(r) = r

"""
    CapsuleAlbedo(; albedo = 0.2, radius = nothing, outside_absorptance = 0.0)

Albedo of a capsule that absorbs as a sphere of effective radius `R(t)` inside its
fixed mesh of radius `R₀`, plus a share `outside_absorptance` of what crosses the mesh
outside that sphere:

    α_eff(t) = 1 - (1 - α_c) s² - a_out (1 - s²),    s = R(t)/R₀,

so that the capsule mesh absorbs `(1 - α_c) H 4πR(t)²`, the absorption of the smaller
sphere at the same incident flux, plus `a_out H 4π(R₀² - R(t)²)`. The mesh re-emits the
rest, which approximates light crossing the surroundings of the shrunken capsule.
`albedo` ``α_c`` is a number or an `AlbedoModel`, such as `TabulatedAlbedo` in time or
`TemperatureAlbedo` in the incident radiation temperature. `radius` is `nothing`
(R = R₀), a [`RadiusSchedule`](@ref) such as [`TabulatedRadius`](@ref) or
[`ScaledRadius`](@ref), a `(times, R/R₀)` table, or a function `t -> R/R₀`. A table is in
ns, linear and held at its ends, and it is stored as a `TabulatedRadius`.
`outside_absorptance` ``a_{out}`` (``0 ≤ a_{out} ≤ 1``, default 0) is the absorptance of
the blow-off between the ablation front and `R₀` for the flux that misses the front. This
explicit term for absorption outside the front, as by W-doped ablator plasma, is
unpublished, so it is 0 unless a study sets it. With the default the model is only the
absorption of the front, and without a `radius` the term vanishes. `energy_balance` then
books the absorbed capsule power at the effective radius plus the outside term, and
`energy_budget` books its integral ``E_{cap}``.

The default ``α_c`` = 0.2 is the 0-D value ``f_{cap} = 1 - α_{cap} ≈ 0.8`` of Swadling et
al., PoP 32, 052707 (2025). For a NIF HDC capsule in a drive history, Young et al.'s
quasi-1D HYDRA ablator albedo (PoP 27, 082702 (2020) Fig. 3a, `young2020_albedos()`) is
the published ``α_c(t)``; `examples/n221204.jl` uses it on N221204's time axis.
"""
struct CapsuleAlbedo{F} <: AlbedoModel
    albedo::AlbedoModel
    radius::F
    outside_absorptance::Float64
end
function CapsuleAlbedo(; albedo = 0.2, radius = nothing, outside_absorptance::Real = 0.0)
    0 <= outside_absorptance <= 1 ||
        throw(ArgumentError("outside_absorptance must lie in [0, 1]"))
    return CapsuleAlbedo(_albedo_model(albedo),
                         _radius_schedule(radius), Float64(outside_absorptance))
end

"""
    radius_fraction(m::CapsuleAlbedo, t) -> R(t)/R₀

The effective radius over the mesh radius of a [`CapsuleAlbedo`](@ref) at time `t` (ns):
that of its `radius` schedule, 1 without one.
"""
radius_fraction(m::CapsuleAlbedo, t) = _radius_fraction(m.radius, t)
_radius_fraction(::Nothing, t) = 1.0
_radius_fraction(r::RadiusSchedule, t) = radius_fraction(r, t)
_radius_fraction(r, t) = r(t)

function albedo(m::CapsuleAlbedo, t, T = nothing)
    αc = albedo(m.albedo, t, T)
    m.radius === nothing && return αc
    s2 = radius_fraction(m, t)^2
    m.outside_absorptance == 0 && return 1 - (1 - αc) * s2
    return 1 - (1 - αc) * s2 - m.outside_absorptance * (1 - s2)
end

# Return the albedo of an albedo model, the albedo it would have intact (closed), and the
# fraction of the incident flux it transmits out of the enclosure.
_albedo_split(m::AlbedoModel, t, T, burnthrough) = (α = albedo(m, t, T); (α, α, zero(α)))
_albedo_split(m::BurnthroughAlbedo, t, T, burnthrough) = _burnthrough_split(m, t, T, burnthrough)
_albedo_split(m::GapClosure, t, T, burnthrough) = _closure_split(m, t, T)
_albedo_split(m::LipClosure, t, T, closed) = _lip_split(m, t, T, closed)

# `_WallHistory` is the wall-loss history of a `WallLossAlbedo` element through a
# `drive_history` step.
# Lindl's wall energy at constant T is E_w ∝ T^a t^b, with a = 4 + p = 3.3 and
# b = 1 + q = 0.62. For a varying T(t) it is carried over as E_w ∝ [∫ r(T) dt]^b with
# r(T) = (T/100 eV)^(a/b), the form of `HeatFront`'s drive integral. Its rate dE_w/dt is
# the constant-T loss at the current T and the effective time t_eff = ∫ r(T) dt / r(T),
# the time the wall would have needed at today's T to hold today's energy. Over a step,
# the trapezoidal integral gives t_eff = A / r(T) + B with A = I₀ + Δ r₀/2 and B = Δ/2.
# At the first step A = 0 and B = max(t, 0), which is constant T since t = 0. No history
# precedes a first step at or before t = 0, such as the step of a prepulse.
#
# Re-emission under a falling drive (`drive_history(...; wall_release = true)`) works as
# follows. The loss above is the heating response. It is evaluated at the highest
# temperature of the element up to the current time, T_h = max(T_max, T). A drop below
# T_max is the deficit D(t) = φ(T_max) − φ(T), with φ = (T/100 eV)^a. The wall responds
# to the deficit with the same self-similar kernel, superposed (Duhamel), so the net flux
# into the wall is
#
#     F = L(T_h, t_eff) − f_r b C ∫ (t − t′)^(b−1) dD(t′),    C = L t_eff^(1−b) / (b φ(T_h)).
#
# Here C is the coefficient of the loss; for gold's unclamped Lindl loss it is the constant
# c σ (100 eV)⁴ / b. A constant-T wall then loses C b φ t^(b−1), and a step down of D
# releases f_r C D (t − t_d)^b by time t, with f_r the model's `release`. D is piecewise
# linear between the steps, and the integral is exact for it. The state contains `G`, the
# part of the earlier steps at the time of this step; `D`, the deficit of the previous
# step; and `Δ`, the length of this step. With no deficit (T never below its running
# maximum) the albedo equals that of the heating loss, with no re-emission term.
struct _WallHistory
    A::Float64
    B::Float64
    Tmax::Float64          # highest temperature before this step (0 if none)
    G::Float64
    D::Float64
    Δ::Float64
end
_WallHistory(A, B) = _WallHistory(A, B, 0.0, 0.0, 0.0, 0.0)

_loss_rate(m::WallLossAlbedo, T) =
    (max(T, zero(T)) / 100)^((4 + m.temperature_exponent) / (1 + m.time_exponent))
# Return the potential φ = (T/100 eV)^a of the wall energy. When φ drops, the wall
# re-emits the released energy.
_wall_potential(m::WallLossAlbedo, T) = (max(T, zero(T)) / 100)^(4 + m.temperature_exponent)

function _effective_time(m::WallLossAlbedo, h::_WallHistory, T)
    Te = (T === nothing || T <= 0) ? m.reference_temperature : T
    h.A == 0 && return h.B
    return h.A / _loss_rate(m, Te) + h.B
end

# The state of a wall-loss element with history `h` at temperature T. The function returns
# the albedo α, an emission E (TW/mm², added to its x-ray source) and the flux Fc that the
# element re-emits from stored energy. The net flux into the wall is
# (1 − α) H − E = L(T_h) − Fc. While T is at or above its running maximum, the albedo is
# that of the heating loss and E = Fc. Below the maximum, the implicit part of Fc, the
# term C Δ^(b−1) φ(T) from the deficit of this step, goes into the albedo, and the rest,
# independent of T, goes into E. The hot surface layer responds stiffly at short times, so
# a T-dependent emission would diverge in `solve`'s albedo iteration. The albedo
# 1 − C Δ^(b−1) φ(T)/H ∝ T^-0.7 instead iterates like Lindl's. The radiosity αH + E is
# kept ≥ 0.
function _wall_split(m::WallLossAlbedo, h::_WallHistory, T)
    # A heated wall without flux, on a step where no beam fires, starts below its running
    # maximum, where its stored energy re-emits; `solve`'s iteration then finds its T.
    (T === nothing || T <= 0) && h.Tmax > 0 && (T = h.Tmax / 2)
    Te = (T === nothing || T <= 0) ? m.reference_temperature : T
    if Te >= h.Tmax && h.G == 0 && h.D == 0
        α = albedo(m, _effective_time(m, h, T), T)
        return α, zero(α), zero(α)
    end
    Th = max(Te, h.Tmax)
    b = 1 + m.time_exponent
    K = h.Δ^(b - 1)
    tc = max(_effective_time(m, h, Th), m.t0)
    αh = albedo(m, tc, Th)
    Lh = (1 - αh) * SIGMA_SB * Th^4
    C = m.release * Lh * tc^(1 - b) / (b * _wall_potential(m, Th))
    φ = _wall_potential(m, Te)
    Fc = C * (h.G + (_wall_potential(m, Th) - φ - h.D) * K)
    α, E = Te >= h.Tmax ? (αh, Fc) : (1 - C * K * φ / (SIGMA_SB * Te^4), C * K * φ - (Lh - Fc))
    # A wall that cooled and is heated again takes its surface layer back (Fc < 0), but
    # no more than the flux that falls on it, so the radiosity αH + E stays ≥ 0. Where it
    # would go negative, α = 0 and E = 0, and the re-emitted flux makes the net flux into
    # the wall equal H.
    H = SIGMA_SB * Te^4
    α * H + E >= 0 && return α, E, Fc
    return zero(α), zero(α), Lh - H
end

_albedo_split(m::WallLossAlbedo, t, T, h::_WallHistory) =
    (α = first(_wall_split(m, h, T)); (α, α, zero(α)))

# The state of a thin-foil element in a `drive_history` whose `before` is a
# `WallLossAlbedo` that follows the temperature history of the element. The state contains
# the burnthrough time, NaN when none is known, and that wall history.
struct _FoilHistory
    tb::Float64
    wall::_WallHistory
end

# The per-element state that `drive_history` passes to `solve` in place of the public
# `burnthrough` vector. `values` is that vector: the burnthrough time of each thin-foil
# element (NaN when none is known) and the closed fraction of each closing-lip element.
# `wall` contains the wall-loss history of the elements whose loss follows their
# temperature history, and `nothing` elsewhere.
struct _StepState
    values::Vector{Float64}
    wall::Vector{Union{Nothing,_WallHistory}}
end

# Returns the burnthrough time or closed fraction of element `e`, and its wall history,
# from the `burnthrough` argument of `solve`.
_element_value(::Nothing, e) = NaN
_element_value(b::AbstractVector, e) = b[e]
_element_value(s::_StepState, e) = s.values[e]
_element_wall(b, e) = nothing
_element_wall(s::_StepState, e) = s.wall[e]

# Return the albedo, intact albedo, transmitted fraction and emission of such a foil, and
# the factor on the wall emission. While the foil is intact, these equal the values of the
# wall, from `_wall_split` of `before`. Across the transition they blend
# `_burnthrough_split` with that `before`, and the re-emitted stored energy is scaled with
# the albedo by α/α_before. That factor is 1 - s f, or the blend of an explicit `after`.
function _foil_history_split(m::BurnthroughAlbedo, t, T, st::_FoilHistory)
    αb, E, _ = _wall_split(m.before::WallLossAlbedo, st.wall, T)
    tb = !isnan(st.tb) ? st.tb : m.burnthrough !== nothing ? m.burnthrough : Inf
    s = _burnthrough_step(m, t, tb)
    s == 0 && return αb, αb, zero(αb), E, one(αb)
    if m.after === nothing
        τ = s * _transmission(m.transmission, optical_depth(m), T)
        return αb * (1 - τ), αb, τ + zero(αb), E * (1 - τ), 1 - τ
    end
    α = αb + s * (albedo(m.after, t, T) - αb)
    k = αb > 0 ? α / αb : 1 - s
    return α, αb, αb - α, E * k, k
end
# Test whether a model transmits part of its incident flux out of the enclosure.
_transmits(::AlbedoModel) = false
_transmits(::Union{BurnthroughAlbedo,GapClosure,LipClosure}) = true

# Return the progress of the transition of a model at time t for burnthrough time tb, 0
# before the transition and 1 after.
_transition_step(m::BurnthroughAlbedo, t, tb) = _burnthrough_step(m, t, tb)
_transition_step(m::GapClosure, t, tb) = _closure_step(m, t)

"""
    burnthrough_time(model::BurnthroughAlbedo, times, Tr) -> Float64

Time (ns) at which the heat front driven by the radiation temperature history `Tr`
(eV, at `times`, piecewise linear; `heat_front`) reaches the back of the layer, or `Inf`
if it does not within `times`. The function returns only this estimate and does not
use an explicit `burnthrough` of the model.
"""
function burnthrough_time(m::BurnthroughAlbedo, times, Tr)
    length(times) == length(Tr) || throw(ArgumentError("times and temperatures differ in length"))
    f = m.front
    need = _front_integral(f, areal_mass(m))
    I = 0.0
    for k in 2:length(times)
        r0, r1 = _front_rate(f, Tr[k-1]), _front_rate(f, Tr[k])
        Δ = times[k] - times[k-1]
        dI = Δ * (r0 + r1) / 2
        if I + dI >= need
            # The integrand is linear across the step, so solve the quadratic for the
            # crossing.
            return times[k-1] + _linear_crossing(r0, r1, Δ, need - I)
        end
        I += dI
    end
    return Inf
end

# Return x in [0, Δ] with ∫₀ˣ (r0 + (r1 - r0) s/Δ) ds = J, for 0 ≤ J ≤ Δ(r0 + r1)/2.
function _linear_crossing(r0, r1, Δ, J)
    a = (r1 - r0) / 2Δ
    abs(a) * Δ^2 <= 1e-12 * max(J, eps()) && return r0 > 0 ? J / r0 : Δ
    # Solve a x² + r0 x - J = 0 for the root in [0, Δ], in a form that avoids cancellation.
    return 2J / (r0 + sqrt(max(r0^2 + 4a * J, 0.0)))
end

# Interpolate linearly through (ts, vs), the interpolation used across the package. The
# times ts are sorted, and a time may repeat once for a step, which takes the value after
# it. Outside the table the end values are held. The function is generic in `t` for
# ForwardDiff. The end values come back in the number type of t, because a Float64 beside
# duals made callers infer a Union, and the interval is found on the value of t. Albedo,
# radius and conversion tables, pulses, LPI transfer tables, response kernels and
# drive-history resampling all use it. `TabulatedCurve` uses a separate form,
# v₁ + w (v₂ − v₁), which rounds differently.
function _interp(ts, vs, t)
    t <= ts[1] && return vs[1] * one(t)
    t >= ts[end] && return vs[end] * one(t)
    k = searchsortedlast(ts, ForwardDiff.value(t))
    w = (t - ts[k]) / (ts[k+1] - ts[k])
    return (1 - w) * vs[k] + w * vs[k+1]
end

"""
    Conversion(; efficiency = 0.9, time_exponent = 0.0, reference_time = 1.0, start = 0.0,
               t0 = 0.1, intensities = nothing, efficiencies = nothing, range = (0.0, 1.0))

Laser-to-x-ray conversion efficiency that varies in time and with the absorbed laser
intensity, for `Material(...; conversion = Conversion(...))`:

    η(t, I) = clamp(η₀(I) ((t - start)/reference_time)^s, range...),

with `s` = `time_exponent`, `t - start` clamped below at `t0` (ns), and `η₀(I)` either
the constant `efficiency` or, with `intensities` (W/cm²) and `efficiencies`, the table
interpolated linearly in log₁₀ I and held at its ends. The default `efficiency` 0.9 is the
conversion of [`gold`](@ref), the input of a published static VisRad model of a NIF hohlraum
(Young et al., PoP 27, 082702 (2020)); 0.75 suits 1 ns Nova and Omega hohlraums. Research
notes, "Material defaults and their sources (2026-11-18)". `I` is the absorbed laser flux
`P/A` of each element (×1e14 W/cm²): a property of the element mesh too (coarse elements
average a spot over more area). Evaluate with `conversion_efficiency`.

Published forms (see `nif_conversion` for the preset):
- in time, η ∝ t^0.12–0.2: `0.7 t^0.12` for Nova's 1 ns hohlraums, which makes the
  measured `T ∝ t^0.15` (Lindl et al., PoP 11, 339 (2004), Fig. 4-4 caption);
  `0.7 t^0.2` (Rosen 2005, Nova LASNEX); `0.85 t^0.2` (Schein 2007, Omega);
- in intensity, the backscatter-corrected plateau on Au, U and cocktail spheres:
  0.95 at 10¹⁴ and 0.80 at 10¹⁵ W/cm² (Dewald et al., PoP 15, 072706 (2008),
  abstract), `Conversion(intensities = (1e14, 1e15), efficiencies = (0.95, 0.80))`.
  The η of a hohlraum exceeds that of a disk at the same intensity (Lindl 2004, Fig. 4-10).

Generic in the number type of the laser power (the intensity table is differentiable);
for a material with a plain number, that number is used exactly.
"""
struct Conversion
    efficiency::Float64
    time_exponent::Float64
    reference_time::Float64
    start::Float64
    t0::Float64
    log_intensity::Vector{Float64}     # log₁₀(W/cm²); empty without intensity dependence
    table::Vector{Float64}
    range::Tuple{Float64,Float64}
end

function Conversion(; efficiency::Real = 0.9, time_exponent::Real = 0.0, reference_time::Real = 1.0,
                    start::Real = 0.0, t0::Real = 0.1, intensities = nothing, efficiencies = nothing,
                    range = (0.0, 1.0))
    (intensities === nothing) == (efficiencies === nothing) ||
        throw(ArgumentError("give both `intensities` and `efficiencies`, or neither"))
    li = intensities === nothing ? Float64[] : log10.(collect(Float64, intensities))
    tb = efficiencies === nothing ? Float64[] : collect(Float64, efficiencies)
    length(li) == length(tb) || throw(ArgumentError("intensities and efficiencies differ in length"))
    issorted(li) || throw(ArgumentError("intensities must be sorted"))
    reference_time > 0 && t0 > 0 || throw(ArgumentError("reference_time and t0 must be positive"))
    return Conversion(efficiency, time_exponent, reference_time, start, t0, li, tb,
                      (Float64(range[1]), Float64(range[2])))
end

"""
    conversion_efficiency(model, t = 0.0, I = nothing)

Laser-to-x-ray conversion efficiency of `model` (a `Conversion` or a number) at time `t`
(ns) and absorbed laser intensity `I` (W/cm²; required by an intensity table).
"""
conversion_efficiency(η::Real, t = 0.0, I = nothing) = η
function conversion_efficiency(m::Conversion, t = 0.0, I = nothing)
    η0 = if isempty(m.table)
        m.efficiency
    else
        I === nothing && throw(ArgumentError("this conversion depends on the laser intensity: pass it"))
        _interp(m.log_intensity, m.table, log10(max(I, 1.0)))
    end
    m.time_exponent == 0 && return clamp(η0, m.range...)
    τ = max(t - m.start, m.t0) / m.reference_time
    return clamp(η0 * τ^m.time_exponent, m.range...)
end

_intensity_dependent(::Real) = false
_intensity_dependent(m::Conversion) = !isempty(m.table)

"""
    nif_conversion(; efficiency = 0.7, time_exponent = 0.12, range = (0.0, 0.95), kwargs...)

The conversion preset for NIF hohlraums: Lindl's Nova `0.7 t^0.12` (t from the start
of the pulse, ns) capped at 0.95. Extrapolated to a NIF pulse it gives 0.80 at 3 ns,
0.85 at 5 ns and 0.90 at 8 ns, the 85–90% inferred from NIF vacuum hohlraums at
235–345 eV (Olson et al., PoP 19, 053301 (2012), abstract); the cap is the 10¹⁴ W/cm²
plateau of Dewald 2008 and the value of the only published VisRad energetics
calibration (η 0.95 with albedo 0.77, Au, 400 TW; Swadling et al., PoP 32, 052707
(2025)). The preset is *derived*, and the research notes section "Wall material, conversion
and capsule energetics (2026-10-24)" gives the derivation. It is the conversion of the
drive-history walls [`nif_gold`](@ref) and [`nif_uranium`](@ref); [`gold`](@ref) and
[`uranium`](@ref) keep the constant 0.9.
"""
nif_conversion(; efficiency = 0.7, time_exponent = 0.12, range = (0.0, 0.95), kwargs...) =
    Conversion(; efficiency, time_exponent, range, kwargs...)
