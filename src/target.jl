# Materials, target components, and the assembled `Target` (elements + ray-casting BVH).

"""
    AlbedoModel

Supertype of x-ray albedo models: `ConstantAlbedo`, `TabulatedAlbedo`,
`WallLossAlbedo`, `BurnthroughAlbedo`, `GapClosure`. Evaluate with `albedo`.
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
gold's opacity and equation of state. For a drive `T(t)` the front is integrated as

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

# The integrand (T/100 eV)^(a/b) of the front's drive integral, and the areal mass for
# a value of the integral.
_front_rate(f::HeatFront, T) = (max(T, zero(T)) / 100)^(f.temperature_exponent / f.time_exponent)
_front_mass(f::HeatFront, I) = f.coefficient * max(I, zero(I))^f.time_exponent
_front_integral(f::HeatFront, m) = (m / f.coefficient)^(1 / f.time_exponent)

"""
    WallLossAlbedo(; coefficient = 0.32, temperature_exponent = -0.7, time_exponent = -0.38,
                   reference_temperature = 250.0, t0 = 0.1, range = (0.0, 0.95),
                   factor = 1.0, factor_exponent = 0.0, factor_temperature = 300.0,
                   liner = 0.0, front = HeatFront())

Marshak-wave wall-loss albedo, `1 - α = c (T/100 eV)^p (t/1 ns)^q`, with defaults for
gold from Lindl's scaling (T in heV, t in ns). `t` is clamped below at `t0`, and α to
`range`. The radiation temperature is each element's incident-flux temperature, so
`solve` iterates it to self-consistency; `reference_temperature` (eV) is used where
no flux is known. The form is Lindl's for a wall held at `T` since `t = 0`;
[`drive_history`](@ref) evaluates it at each element's effective time instead, from
its own temperature history (`wall_history`).

**Wall material.** `factor` f and `factor_exponent` δ scale gold's loss to another wall
material's,

    1 - α = r(T) c (T/100 eV)^p (t/1 ns)^q,    r(T) = f (T / factor_temperature)^δ,

applied to gold's loss after its `range` clamp (so a DU wall's albedo may exceed gold's
cap), the ratio of the two materials' Marshak wall losses (Rosen 2005: for U₀.₆Dy₀.₂Au₀.₂
against Au, `E/A ∝ T^3.1 t^0.57` against `T^3.3 t^0.6`, a ratio ∝ T^-0.22 t^-0.05).
`du_wall_loss()` is the preset for an Au-lined depleted-uranium wall.

**Liner.** `liner` (g/cm², 0 for none) of gold on the inside of such a wall: until the
heat front has crossed it the wall loses like gold. The front's areal mass is the
constant-temperature similarity form of `front` (a [`HeatFront`](@ref)) at the
element's own `T` and `t`, `m_F = c_F (T/100 eV)^a t^b` (no history: exact at constant
`T`, within a few percent of the drive integral for slow ramps), and the ratio becomes

    r_eff = w + (1 - w) r(T),    w = min(1, liner / m_F),

the gold liner's share of the heated mass (an energy-content blend, *engineering*:
continuous in `T` and `t`, so the albedo iteration stays smooth; a sharp front would
switch the loss rate to the substrate's as soon as `m_F` passes the liner). With the
defaults (`factor = 1`, no liner) the albedo is gold's, bit for bit.
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
end

function albedo(m::WallLossAlbedo, t, T = nothing)
    Te = (T === nothing || T <= 0) ? m.reference_temperature : T
    tc = max(t, m.t0)
    loss = m.coefficient * (Te / 100)^m.temperature_exponent * tc^m.time_exponent
    α = clamp(1 - loss, m.range...)
    m.factor == 1 && m.factor_exponent == 0 && return α
    r = m.factor * (Te / m.factor_temperature)^m.factor_exponent
    if m.liner > 0
        f = m.front
        w = min(1.0, m.liner / (f.coefficient * (Te / 100)^f.temperature_exponent * tc^f.time_exponent))
        r = w + (1 - w) * r
    end
    return 1 - r * (1 - α)
end

"""
    wall_loss_factor(model::WallLossAlbedo, t, T)

The ratio `r_eff` of the model's wall loss to gold's at time `t` (ns) and radiation
temperature `T` (eV): 1 for gold, `f (T/T_f)^δ` for another material, blended toward 1
while the heat front is inside a gold `liner`.
"""
function wall_loss_factor(m::WallLossAlbedo, t, T)
    g = WallLossAlbedo(; m.coefficient, m.temperature_exponent, m.time_exponent,
                       m.reference_temperature, m.t0, range = (-Inf, Inf))
    d = WallLossAlbedo(; m.coefficient, m.temperature_exponent, m.time_exponent,
                       m.reference_temperature, m.t0, range = (-Inf, Inf), m.factor,
                       m.factor_exponent, m.factor_temperature, m.liner, m.front)
    return (1 - albedo(d, t, T)) / (1 - albedo(g, t, T))
end

# The depleted-uranium loss factor at 300 eV: calibrated so that an Au-lined DU wall is
# worth 6.5 % of peak laser power on a NIF high-foot hohlraum (research notes, "Wall
# material, conversion and capsule energetics (2026-10-24)").
const DU_LOSS_FACTOR = 0.785

"Areal mass (g/cm²) of the 0.7 µm gold liner of Döppner et al., PRL 115, 055001 (2015)."
const DU_LINER = 19.3 * 0.7e-4

"""
    du_wall_loss(; factor = 0.785, factor_exponent = -0.2, liner = 19.3 × 0.7e-4, kwargs...)

`WallLossAlbedo` for an Au-lined depleted-uranium hohlraum wall: gold's loss times
`r(T) = factor (T/300 eV)^factor_exponent`, gold until the heat front crosses the
`liner` (g/cm²). `kwargs` pass through (gold's Lindl coefficients by default).

- `factor_exponent` −0.2: the cocktail/Au loss ratio's temperature scaling, T^-0.22 from
  Rosen's Hammer–Rosen fits (M. D. Rosen, "Fundamentals of ICF hohlraums", UCRL-JRNL-
  211678 / OSTI 885380 (2005), §VI) and T^-0.19 between LASNEX's 97.6 / 90.8 / 88.6 % at
  185 / 275 / 305 eV (Schein et al., PRL 98, 175003 (2007)).
- `factor` 0.785 at 300 eV is *calibrated* in this model: on the NIF high-foot hohlraum
  of Döppner et al. 2015 at peak (5.75 × 9.43 mm, 3.1 mm LEHs, 1.125 mm capsule,
  390 TW × 0.86 coupling at 13 ns, conversion 0.9, wall source in inner/outer bands) the
  DU wall is then worth 6.5 % of laser power against gold at the same capsule radiation
  temperature, the measured equivalence of the first DU shots (Döppner 2015, citing
  Kline 2013 and Callahan 2012; Callahan et al., PoP 27, 072704 (2020): "in both
  simulations and experiments ... about 6.5 % in drive"). Range 0.78–0.79 over 12–14 ns
  and conversion 0.75–0.95. At equal power it raises the capsule Tr by 4.9 eV against
  the measured +5 ± 6 eV of Döppner's T-1 pair (N140311 Au 317 eV, N140520 DU 322 eV,
  391/393 TW). Rosen's "about 0.83" and Schein's 17.5 % (0.825), for the U:Dy:Au
  cocktail, give 5.1 % here: gold's Lindl loss leaves the wall only ~30 % of the x-ray
  source at a NIF peak, so the measured 6.5 % needs a larger loss reduction than a
  wall-dominated accounting (Schein's 17.5 % → 10 %) would.
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

# The point the default foil optical depth is scaled from: τ ∼ 5 for a gold foil
# 2.5 µm thick (Lester et al. 2026, after eq. 1), i.e. 19.3 g/cm³ × 2.5 µm of gold.
const _FOIL_REFERENCE_DEPTH = 5.0
const _FOIL_REFERENCE_AREAL_MASS = 19.3 * 2.5e-4      # g/cm²

_foil_optical_depth(areal_mass) = _FOIL_REFERENCE_DEPTH * areal_mass / _FOIL_REFERENCE_AREAL_MASS

"""
    RosselandTransmission(; reference_temperature = 325.0, exponent = 1.5)

Post-burnthrough transmission of a gold foil that follows the radiation temperature
`T` (eV) of the flux on it, for `BurnthroughAlbedo(; transmission =
RosselandTransmission())`. The foil's optical depth scales with gold's Rosseland
opacity,

    τ(T) = τ_ref (T / reference_temperature)^(-exponent),

and the transmission is `eddington_transmission(τ(T))`, where `τ_ref` is the foil's
`optical_depth` (5 at 2.5 µm). The exponent is that of Hammer & Rosen's power-law fit
to gold's Rosseland opacity, `1/K = g T^1.5 ρ^-0.2` (Phys. Plasmas 10, 1829 (2003), eq.
31, fitted over 100–200 eV), so the mean free path grows as `T^1.5`. The reference
temperature is the middle of the 300–350 eV at which Lester et al. (Phys. Rev. E 114,
L023201 (2026)) evaluate `τ ∼ 5`: at the reference the transmission equals the constant
default. Hammer & Rosen's fit also reproduces that point: with the as-built areal mass
of 2.5 µm of gold at Lester's expanded density of 1 g/cm³ it gives `τ = 5.9` at 325 eV
(*derived*, extrapolating the fit). The density is held at that expansion (`K ∝ ρ^0.2`:
solid gold would double `τ`). Where no flux is known (`T = nothing`) the reference
temperature is used.
"""
Base.@kwdef struct RosselandTransmission
    reference_temperature::Float64 = 325.0
    exponent::Float64 = 1.5
end

# Transmitted fraction of a burned-through foil of optical depth τ_ref at temperature T.
_transmission(x::Real, τref, Tr) = x
function _transmission(m::RosselandTransmission, τref, Tr)
    Te = (Tr === nothing || Tr <= 0) ? m.reference_temperature : Tr
    return eddington_transmission(τref * (Te / m.reference_temperature)^(-m.exponent))
end

"""
    BurnthroughAlbedo(; thickness, density = 19.3, before = WallLossAlbedo(),
                      transmission = eddington_transmission(τ), after = nothing,
                      width = 0.4, burnthrough = nothing, front = HeatFront())

Albedo of a thin wall layer (a thin-foil diagnostic window: `thickness` µm of a
material of `density` g/cm³, gold by default) that re-emits like the wall until the
Marshak heat front (`front`, a `HeatFront`) has heated its whole areal mass
`ρ × thickness`, and then lets part of the flux that reaches it through.

Once burned through, a fraction `T = transmission` of the incident flux passes through
the foil and leaves the hohlraum (a number, or `RosselandTransmission()` for a fraction
that follows the radiation temperature). The rest meets the foil as before and is re-emitted
with the albedo of `before`. Across the transition

    α(t) = α_before (1 - s T),    transmitted fraction s T,

with `s = s((t - t_b)/width)` the smooth step `3x² - 2x³` on [0, 1] and `t_b` the
burnthrough time, so the burned-through albedo is `(1 - T) α_before`. Of the flux
`H A` on the foil, `s T H A` is transmitted and `(1 - s T)(1 - α_before) H A` absorbed
(`energy_balance`, `drive_history`). `before` is an albedo model or a number, by
default the wall's Marshak loss (`WallLossAlbedo`).

The default `T` is the gray-foil estimate of R. S. Lester et al., Phys. Rev. E 114,
L023201 (2026), doi:10.1103/k8dz-pjdx: `T = 4/(3τ + 2)` (their eq. 2,
`eddington_transmission`). The optical depth is their eq. 1, `τ = L/λ_Au(T_h)` at a
fixed Rosseland mean free path λ_Au, scaled linearly in areal mass (`optical_depth`)
from their point `τ ∼ 5` at `L = 2.5 µm`:

    τ = 5 × ρ thickness / (19.3 g/cm³ × 2.5 µm),

so `T ≈ 0.33`, `0.235` and `0.17` for 1.7, 2.5 and 3.6 µm of gold. Limitations:
- λ_Au is held at their evaluation, for a hohlraum radiation temperature `T_h` of
  typically 300–350 eV and an assumed gold density of about 1 g/cm³ (the expanded
  foil), so `T` is constant once burned through. `transmission =
  RosselandTransmission()` lets λ_Au follow the temperature of the flux on the foil
  instead; its density dependence is not modeled.
- Pairing that density with the 2.5 µm as-built thickness is theirs. Taken literally,
  eq. 1 at 1 g/cm³ describes a layer of 0.25 mg/cm², not the foil's 4.8. The scaling
  with thickness (areal mass) is ours, from eq. 1.
- The estimate is gray, planar and quasi-static, and for gold: pass `transmission`
  for another material. Below `τ = 2/3` (0.33 µm of gold) it is capped at 1.

`transmission = 1` makes a burned-through foil an open hole whose whole incident flux
is transmitted.

An explicit `after` (an albedo model or number; exclusive with `transmission`) sets
the burned-through albedo directly instead,

    α(t) = α_before + s (α_after - α_before),

and counts as transmitted only the loss beyond the intact foil's, `(α_before - α) H A`;
the intact foil's own `(1 - α_before) H A` stays absorbed. `after = 0` is the open-hole
bound on the loss with this split (the default before `transmission` existed).

The burnthrough time `t_b` is, in order:
- the per-element time passed to `solve`/`albedos` (`burnthrough = ...`), which
  `drive_history` integrates from each element's own incident-flux temperature
  history (`heat_front`);
- the explicit `burnthrough` (ns), when given;
- otherwise the foil is intact (`t_b = ∞`): a single solve has no history to estimate
  it from. `burnthrough_time` gives the estimate for a temperature history.

Before `t_b` the foil loses exactly what a thick wall loses and transmits nothing: in
the diffusion picture the cold material ahead of a nonlinear heat front does not
change the flux into it, so no early loss is modeled. Measured burnthrough of 0.8–3 µm
gold foils on a 260 eV hohlraum wall (Nova; Lindl et al., Phys. Plasmas 11, 339
(2004), Fig. 4-6) reaches half signal 0.1–0.25 ns after the computed front arrival,
0.19 ns on average, which sets the default `width` (half transmission at
`t_b + width/2`; see the diagnostic-windows tutorial).
"""
struct BurnthroughAlbedo <: AlbedoModel
    thickness::Float64
    density::Float64
    before::AlbedoModel
    after::Union{Nothing,AlbedoModel}       # an explicit burned-through albedo, or
    transmission::Union{Nothing,Float64,RosselandTransmission}  # the transmitted fraction (exactly one is set)
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
            transmission isa RosselandTransmission || 0 <= transmission <= 1 ||
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

Fraction of the incident flux a burned-through foil transmits when the flux on it has
radiation temperature `T` (eV): the model's `transmission`, evaluated at `T` for a
`RosselandTransmission` (at its reference temperature when `T` is `nothing`). `NaN`
for a model with an explicit `after`.
"""
foil_transmission(m::BurnthroughAlbedo, T = nothing) =
    m.after === nothing ? _transmission(m.transmission, optical_depth(m), T) : NaN

# Transmitting fraction of the transition at time t for burnthrough time tb.
function _burnthrough_step(m::BurnthroughAlbedo, t, tb)
    isfinite(tb) || return 0.0
    m.width == 0 && return t >= tb ? 1.0 : 0.0
    x = clamp((t - tb) / m.width, 0.0, 1.0)
    return x^2 * (3 - 2x)
end

# The albedo at time t, the intact (`before`) albedo, and the fraction of the incident
# flux transmitted out of the enclosure.
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
micrometres wide closes as plasma blown off its walls fills it, which this model does
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

# Closed fraction of a gap at time t (0 open, 1 closed).
function _closure_step(m::GapClosure, t)
    m.width == 0 && return t >= m.closure ? 1.0 : 0.0
    x = clamp((t - m.closure) / m.width, 0.0, 1.0)
    return x^2 * (3 - 2x)
end

# The albedo at time t, the closed albedo, and the transmitted fraction.
function _closure_split(m::GapClosure, t, T)
    αa = albedo(m.after, t, T)
    s = _closure_step(m, t)
    return s * αa, αa, (1 - s) + zero(αa)
end

albedo(m::GapClosure, t, T = nothing) = first(_closure_split(m, t, T))

"""
    CapsuleAlbedo(; albedo = 0.2, radius = nothing)

Albedo of a capsule that absorbs as a sphere of effective radius `R(t)` inside its
fixed mesh of radius `R₀`:

    α_eff(t) = 1 - (1 - α_c) (R(t)/R₀)²,

so that the capsule mesh absorbs `(1 - α_c) H 4πR(t)²`, the absorption of the smaller
sphere at the same incident flux, and passes the rest on (re-emitted from the mesh, an
approximation of light crossing the shrunken capsule's surroundings). `albedo` α_c is
a number or an `AlbedoModel` (`TabulatedAlbedo` in time, `TemperatureAlbedo` in the
incident radiation temperature); `radius` is `nothing` (R = R₀), a `(times, R/R₀)`
table (ns; linear, held at its ends) or a function `t -> R/R₀`. `energy_balance` then
books the capsule's absorbed power at the effective radius, and `energy_budget` its
integral E_cap.

The default α_c = 0.2 is the 0-D value f_cap = 1 − α_cap ≈ 0.8 of Swadling et al., PoP
32, 052707 (2025); no measured α_cap(t) for W-doped HDC is published.
"""
struct CapsuleAlbedo{F} <: AlbedoModel
    albedo::AlbedoModel
    radius::F
end
CapsuleAlbedo(; albedo = 0.2, radius = nothing) =
    CapsuleAlbedo(_albedo_model(albedo),
                  radius isa Tuple ? (collect(Float64, radius[1]), collect(Float64, radius[2])) : radius)

"Effective radius over the mesh radius, `R(t)/R₀`, of a `CapsuleAlbedo`."
capsule_radius_fraction(m::CapsuleAlbedo, t) = _radius_fraction(m.radius, t)
_radius_fraction(::Nothing, t) = 1.0
_radius_fraction(r::Tuple, t) = _interp(r[1], r[2], t)
_radius_fraction(r, t) = r(t)

function albedo(m::CapsuleAlbedo, t, T = nothing)
    αc = albedo(m.albedo, t, T)
    m.radius === nothing && return αc
    return 1 - (1 - αc) * capsule_radius_fraction(m, t)^2
end

# An albedo model's albedo, the albedo it would have intact (closed), and the fraction of
# the incident flux it transmits out of the enclosure.
_albedo_split(m::AlbedoModel, t, T, burnthrough) = (α = albedo(m, t, T); (α, α, zero(α)))
_albedo_split(m::BurnthroughAlbedo, t, T, burnthrough) = _burnthrough_split(m, t, T, burnthrough)
_albedo_split(m::GapClosure, t, T, burnthrough) = _closure_split(m, t, T)

# The wall-loss history of a `WallLossAlbedo` element through a `drive_history` step.
# Lindl's wall energy at constant T, E_w ∝ T^a t^b (a = 4 + p = 3.3, b = 1 + q = 0.62),
# is carried to a varying T(t) as E_w ∝ [∫ r(T) dt]^b, r(T) = (T/100 eV)^(a/b) (the
# form of `HeatFront`'s drive integral); its rate dE_w/dt is the constant-T loss at the
# current T and the effective time t_eff = ∫ r(T) dt / r(T): the time the wall would
# have needed at today's T to hold today's energy. Over a step, with the trapezoidal
# integral, t_eff = A / r(T) + B (A = I₀ + Δ r₀/2, B = Δ/2; at the first step A = 0,
# B = t, i.e. constant T since t = 0).
struct _WallHistory
    A::Float64
    B::Float64
end

_loss_rate(m::WallLossAlbedo, T) =
    (max(T, zero(T)) / 100)^((4 + m.temperature_exponent) / (1 + m.time_exponent))

function _effective_time(m::WallLossAlbedo, h::_WallHistory, T)
    Te = (T === nothing || T <= 0) ? m.reference_temperature : T
    h.A == 0 && return h.B
    return h.A / _loss_rate(m, Te) + h.B
end

_albedo_split(m::WallLossAlbedo, t, T, h::_WallHistory) =
    (α = albedo(m, _effective_time(m, h, T), T); (α, α, zero(α)))

# Whether a model lets part of the flux on it leave the enclosure.
_transmits(::AlbedoModel) = false
_transmits(::Union{BurnthroughAlbedo,GapClosure}) = true

# Progress of a model's transition at time t (burnthrough time tb): 0 before, 1 after.
_transition_step(m::BurnthroughAlbedo, t, tb) = _burnthrough_step(m, t, tb)
_transition_step(m::GapClosure, t, tb) = _closure_step(m, t)

"""
    burnthrough_time(model::BurnthroughAlbedo, times, Tr) -> Float64

Time (ns) at which the heat front driven by the radiation temperature history `Tr`
(eV, at `times`, piecewise linear; `heat_front`) reaches the back of the layer, or `Inf`
if it does not within `times`. This is the estimate alone; an explicit `burnthrough`
of the model is not consulted.
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
            # the integrand is linear across the step: solve the quadratic for the crossing
            return times[k-1] + _linear_crossing(r0, r1, Δ, need - I)
        end
        I += dI
    end
    return Inf
end

# x in [0, Δ] with ∫₀ˣ (r0 + (r1 - r0) s/Δ) ds = J (0 ≤ J ≤ Δ(r0 + r1)/2).
function _linear_crossing(r0, r1, Δ, J)
    a = (r1 - r0) / 2Δ
    abs(a) * Δ^2 <= 1e-12 * max(J, eps()) && return r0 > 0 ? J / r0 : Δ
    # a x² + r0 x - J = 0, root in [0, Δ], written to avoid cancellation
    return 2J / (r0 + sqrt(max(r0^2 + 4a * J, 0.0)))
end

function _interp(ts, vs, t)
    t <= ts[1] && return vs[1]
    t >= ts[end] && return vs[end]
    k = searchsortedlast(ts, t)
    w = (t - ts[k]) / (ts[k+1] - ts[k])
    return (1 - w) * vs[k] + w * vs[k+1]
end

"""
    Conversion(; efficiency = 0.75, time_exponent = 0.0, reference_time = 1.0, start = 0.0,
               t0 = 0.1, intensities = nothing, efficiencies = nothing, range = (0.0, 1.0))

Laser → x-ray conversion efficiency that varies in time and with the absorbed laser
intensity, for `Material(...; conversion = Conversion(...))`:

    η(t, I) = clamp(η₀(I) ((t - start)/reference_time)^s, range...),

with `s` = `time_exponent`, `t - start` clamped below at `t0` (ns), and `η₀(I)` either
the constant `efficiency` or, with `intensities` (W/cm²) and `efficiencies`, the table
interpolated linearly in log₁₀ I and held at its ends. `I` is each element's absorbed
laser flux `P/A` (×1e14 W/cm²): a property of the element mesh too (coarse elements
average a spot over more area). Evaluate with `conversion_efficiency`.

Published forms (see `nif_conversion` for the preset):
- in time, η ∝ t^0.12–0.2: `0.7 t^0.12` for Nova's 1 ns hohlraums, which makes the
  measured `T ∝ t^0.15` (Lindl et al., PoP 11, 339 (2004), Fig. 4-4 caption);
  `0.7 t^0.2` (Rosen 2005, Nova LASNEX); `0.85 t^0.2` (Schein 2007, Omega);
- in intensity, the backscatter-corrected plateau on Au, U and cocktail spheres:
  0.95 at 10¹⁴ and 0.80 at 10¹⁵ W/cm² (Dewald et al., PoP 15, 072706 (2008),
  abstract), `Conversion(intensities = (1e14, 1e15), efficiencies = (0.95, 0.80))`.
  A hohlraum's η exceeds a disk's at the same intensity (Lindl 2004, Fig. 4-10).

Generic in the number type of the laser power (the intensity table is differentiable);
a material with a plain number keeps that number, bit for bit.
"""
struct Conversion
    efficiency::Float64
    time_exponent::Float64
    reference_time::Float64
    start::Float64
    t0::Float64
    log_intensity::Vector{Float64}     # log₁₀(W/cm²); empty: no intensity dependence
    table::Vector{Float64}
    range::Tuple{Float64,Float64}
end

function Conversion(; efficiency::Real = 0.75, time_exponent::Real = 0.0, reference_time::Real = 1.0,
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

Laser → x-ray conversion efficiency of `model` (a `Conversion` or a number) at time `t`
(ns) and absorbed laser intensity `I` (W/cm²; required by an intensity table).
"""
conversion_efficiency(η::Real, t = 0.0, I = nothing) = η
function conversion_efficiency(m::Conversion, t = 0.0, I = nothing)
    η0 = if isempty(m.table)
        m.efficiency
    else
        I === nothing && throw(ArgumentError("this conversion depends on the laser intensity: pass it"))
        _interp_ad(m.log_intensity, m.table, log10(max(I, 1.0)))
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
0.85 at 5 ns and 0.90 at 8 ns, the 85–90 % inferred from NIF vacuum hohlraums at
235–345 eV (Olson et al., PoP 19, 053301 (2012), abstract); the cap is the 10¹⁴ W/cm²
plateau of Dewald 2008 and the value of the only published VisRad energetics
calibration (η 0.95 with albedo 0.77, Au, 400 TW; Swadling et al., PoP 32, 052707
(2025)). *Derived* (research notes, "Wall material, conversion and capsule energetics
(2026-10-24)"); the materials' default stays the constant 0.75.
"""
nif_conversion(; efficiency = 0.7, time_exponent = 0.12, range = (0.0, 0.95), kwargs...) =
    Conversion(; efficiency, time_exponent, range, kwargs...)

"""
    Material(name; albedo = 0.0, conversion = 0.0)

Surface properties for the radiation model. `albedo` (a number or `AlbedoModel`) is the
x-ray re-emission fraction; `conversion` (a number or a [`Conversion`](@ref)) is the
fraction of absorbed laser power re-emitted as x-rays.
"""
struct Material
    name::Symbol
    albedo::AlbedoModel
    conversion::Union{Float64,Conversion}
    Material(name, albedo, conversion) =
        new(name, albedo, conversion isa Conversion ? conversion : Float64(conversion))
end

Material(name::Symbol; albedo = 0.0, conversion = 0.0) =
    Material(name, albedo isa AlbedoModel ? albedo : ConstantAlbedo(albedo), conversion)

albedo(m::Material, t = 0.0) = albedo(m.albedo, t)

"Gold hohlraum wall (defaults: albedo 0.8, laser→x-ray conversion 0.75)."
gold(; albedo = 0.8, conversion = 0.75) = Material(:Au; albedo, conversion)
"Au-lined depleted-uranium hohlraum wall (defaults: `du_wall_loss()`, conversion 0.75, as gold's)."
uranium(; albedo = du_wall_loss(), conversion = 0.75) = Material(:DU; albedo, conversion)
"Capsule ablator (defaults: albedo 0.3, conversion 0.3)."
ablator(; albedo = 0.3, conversion = 0.3) = Material(:ablator; albedo, conversion)
"Perfect absorber that emits nothing."
const BLACK = Material(:black; albedo = 0.0, conversion = 0.0)

"""
    Component(name, mesh, material; kind = :surface, two_sided = false, surface = nothing)

A named piece of a target. `kind = :surface` components are opaque: they block laser
rays and lines of sight, absorb laser light, and radiate from their front side (and
back side when `two_sided`). `kind = :aperture` components (e.g. laser entrance holes)
are transparent diagnostic surfaces that record the x-ray flux crossing them and
re-emit nothing, or, with a `GapClosure` albedo, re-emit what a closing lip of plasma
would (see `aperture`). `surface` is the exact `CurvedSurface` the mesh approximates, if any;
the radiation model then uses the curved geometry (see `view_factors`).
"""
struct Component
    name::Symbol
    mesh::Mesh
    material::Material
    kind::Symbol
    two_sided::Bool
    surface::Union{Nothing,CurvedSurface}
    function Component(name, mesh, material; kind::Symbol = :surface, two_sided::Bool = false,
                       surface::Union{Nothing,CurvedSurface} = nothing)
        kind in (:surface, :aperture) || throw(ArgumentError("kind must be :surface or :aperture"))
        kind === :aperture && two_sided && throw(ArgumentError("apertures are one-sided"))
        new(Symbol(name), mesh, material, kind, two_sided, surface)
    end
end

"""
    aperture(name, mesh; surface = nothing, closure = nothing)

Transparent diagnostic surface (e.g. a laser entrance hole): rays pass through it and
its incident x-ray flux is the flux escaping through it.

With `closure`, a `GapClosure`, it is an opening that plasma closes to x-rays during the
pulse while laser light and lines of sight still pass (the lip of a closing laser entrance
hole, `hohlraum(; leh_closure)`): it re-emits the fraction α(t) = s α_after of its
incident flux and lets `1 - s` escape, booked as escaped by `energy_balance`. It
converts no laser light. The Monte Carlo reference (`Briareus.montecarlo_incident`) lets
paths escape through open apertures and re-emits them diffusely from a lip with the
lip's albedo (given `albedo = albedos(target, time)`).
"""
aperture(name, mesh; surface = nothing, closure::Union{Nothing,GapClosure} = nothing) =
    Component(name, mesh, closure === nothing ? BLACK : Material(:plasma; albedo = closure, conversion = 0.0);
              kind = :aperture, surface)

"Component with its mesh rigidly moved (`x ↦ R*x + t`)."
transform(c::Component, R::AbstractMatrix, t = zero(Vec3)) =
    Component(c.name, transform(c.mesh, R, Vec3(t...)), c.material; kind = c.kind,
              two_sided = c.two_sided, surface = _transform_surface(c.surface, R, Vec3(t...)))
translate(c::Component, t) = transform(c, SMatrix{3,3,Float64}(I), t)

"""
    Target(components; closed = false, visibility = :centroid)

Assembled target geometry. Each patch side that radiates becomes an *element* (the
unknowns of the radiation model): one element per patch for one-sided surfaces and
apertures, two for two-sided surfaces. The back side of a one-sided surface absorbs
laser light as a loss but is not modeled radiatively.

`closed = true` declares that the elements (apertures included) form a closed enclosure,
so `RadiationModel` balances its view factors by default (`enforce_closure!`).
`hohlraum` sets it; a target assembled from its components does not, since added parts
(a foil, a sombrero) may lie outside the enclosure.

`visibility` is the `view_factors` visibility rule `RadiationModel` uses by default:
`:centroid` (one ray per pair) or `:sampled` (the visible fraction of pairs partly
shadowed by a rim). Targets with rims set `:sampled`: `cad_target` and `visrad_target`
(bored entrance holes, gaps, CAD walls), where one centroid ray leaves row sums up to
1% off on walls and 6% on LEH apertures; `hohlraum` keeps `:centroid` (flush apertures,
row-sum errors ~0.1% that balancing removes; research notes, "Partial visibility").
"""
struct Target
    components::Vector{Component}
    # triangles
    tri_v0::Vector{Vec3}
    tri_e1::Vector{Vec3}
    tri_e2::Vector{Vec3}
    tri_normal::Vector{Vec3}
    tri_area::Vector{Float64}
    tri_comp::Vector{Int32}
    tri_front::Vector{Int32}       # element on the front side (0 = none)
    tri_back::Vector{Int32}        # element on the back side (0 = none)
    tri_opaque::BitVector
    # elements
    elem_centroid::Vector{Vec3}
    elem_normal::Vector{Vec3}      # unit normal of the radiating side
    elem_area::Vector{Float64}
    elem_comp::Vector{Int32}
    elem_side::Vector{Int8}        # +1 front, -1 back
    elem_tris::Vector{UnitRange{Int}}
    elem_curved::Vector{Int8}      # radiating side relative to the curved surface's outward normal (0: flat)
    comp_elems::Vector{Vector{Int}}
    bvh::BVH
    scale::Float64                 # bounding-sphere radius about the origin (mm)
    closed::Bool                   # the elements form a closed enclosure
    visibility::Symbol             # default `view_factors` visibility rule
end

function Target(components::AbstractVector{Component}; closed::Bool = false,
                visibility::Symbol = :centroid)
    visibility in (:centroid, :sampled) ||
        throw(ArgumentError("visibility must be :centroid or :sampled, got :$visibility"))
    names = [c.name for c in components]
    allunique(names) || throw(ArgumentError("component names must be unique: $names"))
    tv0, te1, te2, tn, ta = Vec3[], Vec3[], Vec3[], Vec3[], Float64[]
    tcomp, tfront, tback = Int32[], Int32[], Int32[]
    topaque = BitVector()
    ec, en, ea, ecomp, eside, ecurv = Vec3[], Vec3[], Float64[], Int32[], Int8[], Int8[]
    etris = UnitRange{Int}[]
    comp_elems = Vector{Int}[]
    for (ci, c) in enumerate(components)
        m = c.mesh
        t_offset = length(tv0)
        for (a, b, d) in m.triangles
            v0, v1, v2 = m.vertices[a], m.vertices[b], m.vertices[d]
            cr = cross(v1 - v0, v2 - v0)
            push!(tv0, v0); push!(te1, v1 - v0); push!(te2, v2 - v0)
            ar = norm(cr) / 2
            push!(ta, ar); push!(tn, ar > 0 ? cr / (2ar) : zero(Vec3))
            push!(tcomp, ci); push!(tfront, 0); push!(tback, 0)
            push!(topaque, c.kind === :surface)
        end
        mine = Int[]
        for side in (c.two_sided ? (1, -1) : (1,))
            for r in m.patches
                gr = r .+ t_offset
                A = sum(ta[gr])
                A > 0 || continue
                cen = sum(ta[k] * (tv0[k] + (te1[k] + te2[k]) / 3) for k in gr) / A
                nrm = normalize(sum(ta[k] * tn[k] for k in gr))
                cen = _onto_patch(tv0, te1, te2, gr, cen, nrm)
                curv = Int8(0)
                if c.surface !== nothing && !_at_apex(c.surface, tv0, te1, te2, gr)
                    # radiation uses the curved patch: its area, and normals on the side
                    # the flat element radiates from
                    A = c.surface isa DiskSurface ?
                        _sector_area(_sector(c.surface, _patch_vertices(tv0, te1, te2, gr))) :
                        sum(_curved_area(c.surface, tv0[k], tv0[k] + te1[k], tv0[k] + te2[k]) for k in gr)
                    curv = Int8(sign(dot(side * nrm, last(surface_point(c.surface, cen)))))
                end
                push!(ec, cen); push!(en, side * nrm); push!(ea, A); push!(ecurv, curv)
                push!(ecomp, ci); push!(eside, side); push!(etris, gr)
                e = length(ec)
                push!(mine, e)
                for k in gr
                    side == 1 ? (tfront[k] = e) : (tback[k] = e)
                end
            end
        end
        push!(comp_elems, mine)
    end
    opaque_idx = findall(topaque)
    bvh = _sah_bvh(BVH(tv0[opaque_idx], te1[opaque_idx], te2[opaque_idx], Int32.(opaque_idx)))
    scale = isempty(tv0) ? 1.0 :
            maximum(max(norm(v), norm(v + a), norm(v + b)) for (v, a, b) in zip(tv0, te1, te2))
    return Target(collect(components), tv0, te1, te2, tn, ta, tcomp, tfront, tback, topaque,
                  ec, en, ea, ecomp, eside, etris, ecurv, comp_elems, bvh, scale, closed, visibility)
end

Target(cs::Component...; kwargs...) = Target(collect(cs); kwargs...)

# Whether an element touches the apex of the cone it lies on: the true surface is
# singular there (the kernel between points on either side of the apex grows as 1/s²
# toward it), which the curved elements' point rules miss by a fixed fraction however
# fine the mesh (the ring around the apex of a pointed cone saw 3.4% too little); such an
# element stays flat, where near pairs use the exact polygon kernel with adaptive
# emitters.
_at_apex(s, tv0, te1, te2, gr) = false
function _at_apex(s::ConeSurface, tv0, te1, te2, gr)
    dmin, dmax = Inf, 0.0
    for k in gr, v in (tv0[k], tv0[k] + te1[k], tv0[k] + te2[k])
        d = norm(v - s.apex)
        dmin, dmax = min(dmin, d), max(dmax, d)
    end
    return dmin <= 1e-3 * dmax
end

# A folded patch (a remeshed CAD element spans facets within a few degrees of one
# plane) has its area centroid off its facets: behind them on a convex surface, where
# rays from the centroid (visibility, imaging) start inside the solid and are blocked by
# the patch itself. Such a centroid is moved along the mean normal onto the facet the
# line through it meets (`remesh_mesh` makes sure it meets one); a planar patch keeps
# its centroid unchanged.
function _onto_patch(tv0, te1, te2, gr, c, n)
    length(gr) > 1 || return c
    h = 0.0
    for k in gr, v in (tv0[k], tv0[k] + te1[k], tv0[k] + te2[k])
        h = max(h, norm(v - c))
    end
    best = Inf
    for k in gr
        nt = cross(te1[k], te2[k])
        den = dot(n, nt)
        den == 0 && continue
        a = tv0[k]
        s = dot(a - c, nt) / den
        p = c + s * n
        b, d = a + te1[k], a + te2[k]
        tol = -1e-9 * dot(nt, nt)
        dot(_tri_cross(a, b, p), nt) >= tol && dot(_tri_cross(b, d, p), nt) >= tol &&
            dot(_tri_cross(d, a, p), nt) >= tol && abs(s) < abs(best) && (best = s)
    end
    return isfinite(best) && abs(best) > 1e-9 * h ? c + best * n : c
end

"Number of radiating elements (patch sides) of a target."
nelements(t::Target) = length(t.elem_area)
ntriangles(t::Target) = length(t.tri_area)

"""
    component_index(target, name) -> Int
"""
function component_index(t::Target, name::Symbol)
    i = findfirst(c -> c.name === name, t.components)
    i === nothing && throw(KeyError(name))
    return i
end

"""
    elements(target, names...) -> Vector{Int}

Element indices of the named components (both sides of two-sided components).
"""
elements(t::Target, names::Symbol...) =
    reduce(vcat, (t.comp_elems[component_index(t, n)] for n in names); init = Int[])

"""
    element_material(target, e)
"""
element_material(t::Target, e::Integer) = t.components[t.elem_comp[e]].material
is_aperture(t::Target, e::Integer) = t.components[t.elem_comp[e]].kind === :aperture

function Base.show(io::IO, ::MIME"text/plain", t::Target)
    println(io, "Target with ", length(t.components), " components, ", nelements(t),
            " elements, ", ntriangles(t), " triangles")
    for (i, c) in enumerate(t.components)
        A = sum(t.elem_area[e] for e in t.comp_elems[i] if t.elem_side[e] == 1; init = 0.0)
        @printf(io, "  %-14s %-8s %-9s %6d elements  area %9.4f mm²%s\n", c.name, c.kind,
                c.material.name, length(t.comp_elems[i]), A, c.two_sided ? "  (two-sided)" : "")
    end
end
Base.show(io::IO, t::Target) = print(io, "Target(", length(t.components), " components, ",
                                     nelements(t), " elements)")

# Small offset used to start rays just off a surface (mm).
_ray_eps(t::Target) = 1e-7 * max(t.scale, 1.0)

"""
    visible(target, a, b) -> Bool

Whether the open segment between points `a` and `b` is free of opaque surfaces.
"""
visible(t::Target, a::Vec3, b::Vec3) = !any_hit(t.bvh, a, b - a, 1 - 1e-9)

"""
    element_quadrature(target, e, level) -> (points, weights, normals)

Quadrature over element `e` from `4^level` sub-triangles per triangle: sub-triangle
centroids weighted by area, with the element normal. On a component with a
`CurvedSurface`, the points are projected onto the surface and carry its normals
(on the element's radiating side) and curved-patch areas; an element of more than two
facets (imported CAD) instead gets as many points as a two-triangle element, from its
facets' pieces lumped into cells of the element's tangent plane. Weights sum to the
element area.
"""
function element_quadrature(t::Target, e::Integer, level::Integer)
    pts = Vec3[]
    ws = Float64[]
    for k in t.elem_tris[e]
        _subtri!(pts, ws, t.tri_v0[k], t.tri_v0[k] + t.tri_e1[k], t.tri_v0[k] + t.tri_e2[k], level)
    end
    s = t.components[t.elem_comp[e]].surface
    if s === nothing || t.elem_curved[e] == 0
        return pts, ws, fill(t.elem_normal[e], length(pts))
    end
    if s isa DiskSurface
        r0, r1, a, b = element_sector(t, e)
        nr, nφ = 2^level, 2^(level + 1)
        rs = sqrt.(range(r0^2, r1^2; length = nr + 1))        # equal-area rings
        φs = range(a, b; length = nφ + 1)
        cpts, cws = Vec3[], Float64[]
        for i in 1:nr, k in 1:nφ
            push!(cpts, _sector_point(s, sqrt((rs[i]^2 + rs[i+1]^2) / 2), (φs[k] + φs[k+1]) / 2))
            push!(cws, (rs[i+1]^2 - rs[i]^2) / 2 * (φs[k+1] - φs[k]))
        end
        return cpts, cws, fill(t.elem_curved[e] * s.normal, length(cpts))
    end
    length(t.elem_tris[e]) > 2 && return _lumped_quadrature(t, e, level, s)
    cpts, cws, cns = Vec3[], Float64[], Vec3[]
    for k in t.elem_tris[e]
        _curved_subtri!(cpts, cws, cns, s, t.elem_curved[e], t.tri_v0[k], t.tri_v0[k] + t.tri_e1[k],
                        t.tri_v0[k] + t.tri_e2[k], level)
    end
    return cpts, cws, cns
end

# Quadrature of a curved element made of more than two facets (CAD faceting, e.g. a
# remeshed element across many strips of a cylinder), with as many points as a
# two-triangle element would have, wherever its facets lie: each facet is bisected along
# its longest edge until that edge is under half a cell, the pieces (projected onto
# the surface, exact curved areas) are lumped into square cells of side
# size / (2^level √2) in the element's tangent plane (`size` = √2 × the largest
# centroid-to-vertex distance, as in `view_factors`; a square element then has 2·4^level
# cells), and each cell becomes one point: the projection of its pieces' weighted mean
# point, with their total area. Weights sum to the element's exact curved area.
function _lumped_quadrature(t::Target, e::Integer, level::Integer, s::CurvedSurface)
    C, n, sgn = t.elem_centroid[e], t.elem_normal[e], t.elem_curved[e]
    r = 0.0
    for k in t.elem_tris[e], v in (t.tri_v0[k], t.tri_v0[k] + t.tri_e1[k], t.tri_v0[k] + t.tri_e2[k])
        r = max(r, norm(v - C))
    end
    hq = r / 2^level                                   # √2 r / (2^level √2)
    e1 = normalize(cross(n, abs(n[1]) < 0.9 ? Vec3(1, 0, 0) : Vec3(0, 1, 0)))
    e2 = cross(n, e1)
    cells = Dict{Tuple{Int,Int},Tuple{Float64,Vec3}}()
    for k in t.elem_tris[e]
        a = t.tri_v0[k]
        _lump_pieces!(cells, s, C, e1, e2, hq, a, a + t.tri_e1[k], a + t.tri_e2[k])
    end
    pts, ws, ns = Vec3[], Float64[], Vec3[]
    for key in sort!(collect(keys(cells)))
        w, wp = cells[key]
        w > 0 || continue
        p, nn = surface_point(s, wp / w)
        push!(pts, p); push!(ws, w); push!(ns, sgn * nn)
    end
    return pts, ws, ns
end

function _lump_pieces!(cells, s, C, e1, e2, hq, a, b, c)
    lab, lbc, lca = sum(abs2, b - a), sum(abs2, c - b), sum(abs2, a - c)
    if max(lab, lbc, lca) > (hq / 2)^2
        if lab >= lbc && lab >= lca
            m = (a + b) / 2
            _lump_pieces!(cells, s, C, e1, e2, hq, a, m, c)
            _lump_pieces!(cells, s, C, e1, e2, hq, m, b, c)
        elseif lbc >= lca
            m = (b + c) / 2
            _lump_pieces!(cells, s, C, e1, e2, hq, a, b, m)
            _lump_pieces!(cells, s, C, e1, e2, hq, a, m, c)
        else
            m = (c + a) / 2
            _lump_pieces!(cells, s, C, e1, e2, hq, a, b, m)
            _lump_pieces!(cells, s, C, e1, e2, hq, m, b, c)
        end
        return
    end
    g = (a + b + c) / 3
    key = (floor(Int, dot(g - C, e1) / hq), floor(Int, dot(g - C, e2) / hq))
    w = _curved_area(s, a, b, c)
    p = first(surface_point(s, g))
    w0, p0 = get(cells, key, (0.0, zero(Vec3)))
    cells[key] = (w0 + w, p0 + w * p)
    return
end

_patch_vertices(tv0, te1, te2, gr) =
    unique(vcat([[tv0[k], tv0[k] + te1[k], tv0[k] + te2[k]] for k in gr]...))

"""
    element_sector(target, e) -> (r0, r1, φa, φb)

Polar sector of an element of a component with a `DiskSurface`.
"""
element_sector(t::Target, e::Integer) =
    _sector(t.components[t.elem_comp[e]].surface,
            _patch_vertices(t.tri_v0, t.tri_e1, t.tri_e2, t.elem_tris[e]))

function _curved_subtri!(pts, ws, ns, s, sgn, a, b, c, level)
    if level == 0
        p, n = surface_point(s, (a + b + c) / 3)
        push!(pts, p); push!(ns, sgn * n); push!(ws, _curved_area(s, a, b, c))
        return
    end
    ab, bc, ca = (a + b) / 2, (b + c) / 2, (c + a) / 2
    _curved_subtri!(pts, ws, ns, s, sgn, a, ab, ca, level - 1)
    _curved_subtri!(pts, ws, ns, s, sgn, ab, b, bc, level - 1)
    _curved_subtri!(pts, ws, ns, s, sgn, ca, bc, c, level - 1)
    _curved_subtri!(pts, ws, ns, s, sgn, ab, bc, ca, level - 1)
end

function _subtri!(pts, ws, a, b, c, level)
    if level == 0
        push!(pts, (a + b + c) / 3)
        push!(ws, norm(cross(b - a, c - a)) / 2)
        return
    end
    ab, bc, ca = (a + b) / 2, (b + c) / 2, (c + a) / 2
    _subtri!(pts, ws, a, ab, ca, level - 1)
    _subtri!(pts, ws, ab, b, bc, level - 1)
    _subtri!(pts, ws, ca, bc, c, level - 1)
    _subtri!(pts, ws, ab, bc, ca, level - 1)
end
