# This file models inner-beam absorption in the wall bubble. The laser–plasma model moves
# part of the power of the inner cones to the footprint of the outer cones as the Au
# bubble grows into their path.
#
# The model follows Izumi et al., Phys. Plasmas 28, 022706 (2021), Eqs. 3 and 6, with the
# clipping computed from the Briareus beams instead of an annulus. See
# `reference/research-notes.md`, "Inner-beam bubble absorption".

# Calibration (research notes, "Hohlraum dynamics calibrated to published data
# (2026-10-24)"). The tip speed and launch are fitted to Izumi 2021 Fig. 8a (28 measured
# points, seven shots, `izumi2021_bubble()`). The opacity is fitted to the HYDRA absorbed
# fraction of Fig. 10a (`izumi2021_absorption()`) on the basecamp geometry with the
# bubble at z = ±3 mm.
const _BUBBLE_SPEED = 0.36      # mm/ns at 0.3 mg/cc and the basecamp picket (0.347–0.372)
const _BUBBLE_DELAY = 0.25      # ns after the outer turn-on (0.07–0.42)
const _BUBBLE_OPACITY = 0.111   # τ_b (0.100–0.124 statistical; 0.085–0.158 with the plane and launch)
const _BUBBLE_Q = 1.8           # ρ_bub² ∝ (S ρ₀)², S ∝ ρ₀^-0.1 (Izumi 2021 Eq. 7)

# Basecamp geometry of Izumi 2018/2021: hohlraum and capsule outer radius (mm), and the
# ablator–wall collision 5.5 ns after the outer turn-on (inner spots separate at 5.8 ns).
const _BASECAMP = (radius = 2.875, capsule = 0.908, collision = 5.5)

"""
    BubbleAbsorption(; radius, turn_on = 0.0, fill = 0.3, picket_ratio = 1.0, kwargs...)

Inner-beam absorption in the wall bubble (a [`LaserPlasmaModel`](@ref)). The picket of
the outer cones launches a bubble of wall plasma from the wall under their footprint.
Its tip has moved

    x_bal(t) = v_b · max(0, t − t_L),   t_L = turn_on + delay,
    v_b = speed · E^0.2 · (fill / 0.3)^-0.2

in from the initial wall `radius` R (mm, ns) while it is ballistic. E = `picket_ratio`
is the picket fluence of the outer cones over that of the calibration platform, the HDC
basecamp of Izumi 2021, and is dimensionless. `fill` is the gas fill in mg/cc. E is set
by hand by default. [`picket_ratio`](@ref) computes it from the picket fluence of a
pointing on the wall ([`picket_fluence`](@ref)), for `remake(m; picket_ratio)`, and
[`PicketLinkedBubble`](@ref) computes it for the beams of every evaluation; both are
opt-in. The tip stagnates smoothly at the depth `stagnation` ``x_{max}``
([`bubble_depth`](@ref)). It is ballistic up to ``x_{max} - w`` (w = `stagnation_width`),
then at ``x_{max} - w e^{-(x_{bal} - x_{max} + w)/w}``. In the plane `z_b` of the outer
footprint, at the axial position where the beam axes of the footprint cones reach the
wall, the bubble fills the annulus from R in to radius R − ``x_{bub}``. Each inner beam
loses

    absorbed = f · (1 − T_b),
    T_b = exp(−τ_b (fill/0.3)^q (x_bub/x_ref)^p (R/(R − x_bub))²),

where `f` is the fraction of its power crossing `z_b` inside the bubble and `T_b` the
transmission of the bubble. `f` is computed from the beam itself, with its best-focus
spot carried along its `ncone` f-cone directions, as in the deposition. `T_b` uses the
`opacity` ``τ_b``, the exponent q = `density_exponent` of the fill density and the
exponent p = `depth_exponent` of the depth. The absorbed power is conducted to the wall
under the outer footprint (Ralph et al. 2018). It is shared equally by the `footprint`
beams in the hemisphere of the inner beam, those whose port is on the same side, so it
is deposited where those beams land. Beam powers are conserved, since the model moves
power between beams and loses none. The fraction of the footprint beams that misses the
wall, a few 1e-4 through the entrance holes, applies to the moved power as well.

**Ridge stop.** After the ablator–wall collision at `collision` ``t_c``, the inner power
that leaves the bubble stops at the stagnation ridge, at radius `ridge`·R about the axis
(Izumi 2023), and does not reach the wall. A fraction g(t) of it
([`ridge_fraction`](@ref), rising over ±`ridge_rise` ns about ``t_c``) is deposited by
[`lpi_deposit`](@ref) on an axisymmetric band of the `ridge_components` (the hohlraum
wall). The band is centered on the height where the beam axis crosses that radius,
Gaussian in z with rms `ridge_width` (mm), and area-weighted. The heated ridge radiates
at the conversion of the wall, from the wall just outside it. The rest of the beam power
is deposited as usual. Power is conserved, and the ridge power is booked to the
`per_beam` of the beam. The ridge stop acts only where the deposition goes through
`lpi_deposit`: in [`drive_history`](@ref), in `Scenario` evaluations, in the beam-error
sensitivities ([`mode_sensitivity`](@ref) and the rest of `uncertainty.jl`) and in
[`influence_matrix`](@ref).

The edge of the bubble is smoothed over `edge` (mm, a logistic step), which keeps `f`
smooth in time, in the bubble depth and in the pointing of the beams. Everything is
generic in the number type. `speed`, `opacity`, `stagnation`, `collision` and
`picket_ratio`, and the beams, can carry ForwardDiff duals ([`remake`](@ref)), so the
model can be used in [`drive_history`](@ref)`(...; lpi)`, `Scenario(...; lpi)` and
calibrations. The model composes with CBET as `lpi = (cbet, bubble)`. The transfer acts
first, and the bubble then takes its share of the power in the inner beams.

**Keywords**, by role (defaults in parentheses):
- Shot and geometry: `radius` R (required; the initial wall radius, mm), `turn_on`
  (0.0; the turn-on of the outer cones, ns), `fill` (0.3 mg/cc), `picket_ratio` E (1.0),
  `capsule_radius` ``r_{cap}`` (0.908 mm), `inner` (`(:inner23, :inner30)`) and
  `footprint` (`(:outer50,)`), the beam groups that lose and receive the power;
  `wall` (`nothing`), a [`WallMotion`](@ref) that moves the footprint plane with the
  radiating wall. The axes of the outer beams reach the moved wall nearer the entrance
  hole. The bubble tip, the density factor R/(R − ``x_{bub}``) and the ridge stay on the
  initial wall radius in either case: every radius of the model is a laboratory radius
  about the hohlraum axis. The tip positions of Izumi 2021 and the ridge radius of Izumi
  2023 were measured from the initial wall of hohlraums whose walls moved, so the fitted
  `speed`, `opacity` and `ridge` already contain that motion, and a tip measured from
  the moved wall would count it twice. `footprint_z` (`nothing`) fixes |``z_b``| instead.
- Calibration: `speed` (0.36 mm/ns), `delay` (0.25 ns), `opacity` ``τ_b`` (0.111),
  `density_exponent` q (1.8), `depth_exponent` p (0.5) and its reference depth `x_ref`
  (1.0 mm).
- Tip stagnation: `stagnation` (0.80 R) and `stagnation_width` (0.08 R), both defaulting
  to fractions of `radius`; with `stagnation = Inf` the tip stays ballistic.
- Ridge sub-model: `collision` ``t_c`` (derived from the other keywords, below;
  `collision = Inf` turns the ridge stop off), `ridge` (0.73), `ridge_width` (0.5 mm),
  `ridge_rise` (0.2 ns), `ridge_components` (`(:wall,)`).
- Numerics: `edge` (0.05 mm); spot points, `nr` (12) equal-power radii by `nangle` (24)
  angles per f-cone direction, over `ncone` (13) directions.

Defaults (research notes, "Hohlraum dynamics calibrated to published data (2026-10-24)";
the data are [`izumi2021_bubble`](@ref Briareus.Benchmarks.izumi2021_bubble) and [`izumi2021_absorption`](@ref Briareus.Benchmarks.izumi2021_absorption)):
- `speed` 0.36 mm/ns at 0.3 mg/cc and the basecamp picket, ∝ E^0.2 ρ₀^-0.2, and
  `delay` 0.25 ns after the outer turn-on are *fitted* to the 28 measured tip positions
  of Izumi et al., PoP 28, 022706 (2021) Fig. 8a. The data cover seven HDC SymCap shots
  at 0.15–0.3 mg/cc with pickets of 0.73–1.5×, and the fit has χ²/dof 0.89 with the
  ±40 µm of the paper. It gives 0.360 (0.347–0.372) mm/ns and a delay of 0.25
  (0.07–0.42) ns, with the two correlated. Meezan et al. 2020 measured 321–461 µm/ns on
  the same shots.
- `density_exponent` q = 1.8: inverse-bremsstrahlung absorption goes as ``ρ_{bub}^2``,
  and Eq. 7 of Izumi has ``ρ_{bub}`` = S ρ₀ R/(R − x) with the shock compression
  S ∝ ρ₀^-0.1, so ``ρ_{bub}^2`` ∝ ρ₀^1.8 (R/(R − x))². `depth_exponent` 0.5 is the chord
  (``x_{bub}``/1 mm)^0.5 of Eq. 6.
- `opacity` ``τ_b`` 0.111 is *fitted* to the inner-power fraction absorbed in the bubble
  in the 2D HYDRA simulation by Izumi of the 0.3 mg/cc basecamp (Fig. 10a: 0.02, 0.23,
  0.55 at 4, 5, 6 ns). The fit uses the clipping that this model computes for NIF beams
  pointed for a 5.75 mm hohlraum, with the bubble plane at z = ±3 mm as in Izumi
  (`footprint_z = 3.0`). The fitted model absorbs 0.03, 0.20, 0.36 and 0.56 at 4, 5, 5.5
  and 6 ns, with χ² 1.1 at ±0.03. The statistical range is 0.100–0.124. The plane
  (2.8–3.2 mm) and the launch move ``τ_b`` over 0.085–0.158 but leave the absorbed
  fraction, the calibrated quantity, within 0.01. ``T_b`` at 5.5 ns is then 0.36 at
  0.3 mg/cc and 0.58 at 0.15 mg/cc; the model of Izumi gives 0.29 and 0.54, with an
  annulus for f.
- `stagnation` 0.80 R and `stagnation_width` 0.08 R are *engineering* values. The tip is
  exactly ballistic to 0.72 R (2.07 mm in a 2.88 mm hohlraum), past the measured depths
  at 0.3 mg/cc (2.07 mm at 6.2 ns) and near those at 0.15 mg/cc (up to 2.28 mm,
  0.79 R). In Izumi 2023 the tip slows on back-pressure near the axis. A depth of
  R − 0.73 R, the depth of the ridge at the waist, would stop the tip before it reaches
  the inner beams at ``z_b``, 1.3 mm from the wall.
- `collision` ``t_c`` = `turn_on` + 5.5 ns · (R − ``r_{cap}``)/(2.875 − 0.908 mm) · E^-0.16
  is *fitted* to Izumi et al., RSI 89, 10K111, 2018. There the inner spots separate at
  5.8 ns for the 0.3 ns turn-on of the basecamp, at 5.45 ns with a 1.5× picket, and about
  0.4 ns earlier with the 1.18× capsule. The gap between the outer radius of the capsule
  `capsule_radius` ``r_{cap}`` and the wall gives 0.44 ns for that capsule, so a constant
  closing speed of ablator and wall plasma, 0.36 mm/ns, accounts for the capsule shot.
  `ridge` 0.73 corresponds to the ridge of Izumi 2023 at r ≈ 2.1 mm in a 2.88 mm
  hohlraum. `ridge_width` 0.5 mm and `ridge_rise` 0.2 ns are engineering choices and
  were not measured.
"""
struct BubbleAbsorption{S<:Real} <: LaserPlasmaModel
    radius::Float64
    turn_on::Float64
    launch::Float64
    speed::S
    opacity::S
    fill::Float64
    picket_ratio::S
    density_exponent::Float64
    depth_exponent::Float64
    x_ref::Float64
    stagnation::S
    stagnation_width::Float64
    collision::S
    ridge::Float64
    ridge_width::Float64
    ridge_rise::Float64
    inner::Vector{Symbol}
    footprint::Vector{Symbol}
    ridge_components::Vector{Symbol}
    wall::Union{Nothing,WallMotion}
    footprint_z::Union{Nothing,Float64}
    edge::Float64
    nr::Int
    nangle::Int
    ncone::Int
end

# Returns the ablator–wall collision time (Izumi 2018), 5.5 ns after the outer turn-on for
# the basecamp, scaled with the capsule–wall gap and the picket fluence ratio.
_collision_time(turn_on, radius, capsule_radius, picket_ratio) =
    turn_on + _BASECAMP.collision * (radius - capsule_radius) /
              (_BASECAMP.radius - _BASECAMP.capsule) * picket_ratio^-0.16

function BubbleAbsorption(; radius::Real, turn_on::Real = 0.0, delay::Real = _BUBBLE_DELAY,
                          speed::Real = _BUBBLE_SPEED, opacity::Real = _BUBBLE_OPACITY,
                          fill::Real = 0.3, picket_ratio::Real = 1.0,
                          density_exponent::Real = _BUBBLE_Q,
                          depth_exponent::Real = 0.5, x_ref::Real = 1.0,
                          stagnation::Real = 0.80 * radius,
                          stagnation_width::Real = 0.08 * radius,
                          capsule_radius::Real = _BASECAMP.capsule,
                          collision::Real = _collision_time(turn_on, radius, capsule_radius, picket_ratio),
                          ridge::Real = 0.73, ridge_width::Real = 0.5, ridge_rise::Real = 0.2,
                          inner = (:inner23, :inner30), footprint = (:outer50,),
                          ridge_components = (:wall,),
                          wall::Union{Nothing,WallMotion} = nothing, footprint_z = nothing,
                          edge::Real = 0.05, nr::Integer = 12, nangle::Integer = 24,
                          ncone::Integer = 13)
    radius > 0 || throw(ArgumentError("radius must be positive"))
    fill > 0 && picket_ratio > 0 || throw(ArgumentError("fill and picket_ratio must be positive"))
    edge > 0 || throw(ArgumentError("edge must be positive"))
    stagnation_width > 0 && stagnation > stagnation_width ||
        throw(ArgumentError("stagnation must exceed stagnation_width > 0"))
    0 < ridge < 1 && ridge_width > 0 && ridge_rise > 0 ||
        throw(ArgumentError("need 0 < ridge < 1, ridge_width > 0, ridge_rise > 0"))
    nr >= 1 && nangle >= 1 && 1 <= ncone <= 16 || throw(ArgumentError("bad sampling counts"))
    s, τ, xs, tc, E = promote(speed, opacity, stagnation, collision, picket_ratio)
    return BubbleAbsorption(Float64(radius), Float64(turn_on), Float64(turn_on + delay), s, τ,
                            Float64(fill), E, Float64(density_exponent), Float64(depth_exponent),
                            Float64(x_ref), xs, Float64(stagnation_width), tc, Float64(ridge),
                            Float64(ridge_width), Float64(ridge_rise), Symbol[inner...],
                            Symbol[footprint...], Symbol[ridge_components...], wall,
                            footprint_z === nothing ? nothing : Float64(abs(footprint_z)),
                            Float64(edge), Int(nr), Int(nangle), Int(ncone))
end

"""
    remake(m::BubbleAbsorption; speed = m.speed, opacity = m.opacity,
           stagnation = m.stagnation, collision, picket_ratio = m.picket_ratio)

The model with new calibration parameters (e.g. ForwardDiff duals in a fit), or a new
picket fluence ratio E (`picket_ratio`, e.g. from [`picket_ratio`](@ref) for a
re-pointed design). With a new E, the picket dependence of the model enters the ridge
stop as in the constructor, with ``t_c`` − `turn_on` ∝ E^-0.16. `collision` defaults
to `turn_on + (m.collision − turn_on)·(E/m.picket_ratio)^-0.16`. That is the stored
value when E is unchanged, whether the constructor derived it or it was given, and
`Inf` stays `Inf`. An explicit `collision` sets the value. Every argument may be a
ForwardDiff dual.
"""
function remake(m::BubbleAbsorption; speed = m.speed, opacity = m.opacity,
                stagnation = m.stagnation, picket_ratio = m.picket_ratio,
                collision = _rescaled_collision(m, picket_ratio))
    s, τ, xs, tc, E = promote(speed, opacity, stagnation, collision, picket_ratio)
    return BubbleAbsorption(m.radius, m.turn_on, m.launch, s, τ, m.fill, E, m.density_exponent,
                            m.depth_exponent, m.x_ref, xs, m.stagnation_width, tc, m.ridge,
                            m.ridge_width, m.ridge_rise, m.inner, m.footprint, m.ridge_components,
                            m.wall, m.footprint_z, m.edge, m.nr, m.nangle, m.ncone)
end

# Returns the collision time at picket ratio E, scaled from the stored time
# (t_c − turn_on ∝ E^-0.16). It is identical to the stored value when E is the stored
# ratio or the ridge stop is off.
function _rescaled_collision(m::BubbleAbsorption, E)
    E > 0 || throw(ArgumentError("picket_ratio must be positive"))
    (E === m.picket_ratio || isinf(m.collision)) && return m.collision + zero(E)
    return m.turn_on + (m.collision - m.turn_on) * (E / m.picket_ratio)^-0.16
end

function Base.show(io::IO, m::BubbleAbsorption)
    @printf(io, "BubbleAbsorption(R = %.3g mm, v_b = %.3g mm/ns from %.3g ns, τ_b = %.3g, x_max = %.3g mm, ridge from %.3g ns, %s → %s)",
            m.radius, ForwardDiff.value(bubble_speed(m)), m.launch, ForwardDiff.value(m.opacity),
            ForwardDiff.value(m.stagnation), ForwardDiff.value(m.collision),
            join(m.inner, "/"), join(m.footprint, "/"))
end

"Bubble tip speed (mm/ns): `speed · picket_ratio^0.2 · (fill/0.3)^-0.2`."
bubble_speed(m::BubbleAbsorption) = m.speed * m.picket_ratio^0.2 * (m.fill / 0.3)^-0.2

"""
    bubble_depth(m::BubbleAbsorption, t)

Bubble tip displacement from the initial wall (mm) at time `t` (ns). It is the ballistic
depth ``v_b (t - t_L)`` up to the knee `stagnation − stagnation_width`, and beyond the
knee `stagnation − w·exp(−(x − knee)/w)`, which stops at `stagnation` (C¹ at the knee).
"""
function bubble_depth(m::BubbleAbsorption, t)
    x = bubble_speed(m) * max(zero(t), t - m.launch)
    return _stagnate(x, m.stagnation, m.stagnation_width)
end

"""
    ridge_fraction(m::BubbleAbsorption, t) -> g

The fraction of the remaining power of the inner beams that stops at the ablator–wall
ridge at time `t` (ns). It is a C¹ smoothstep from 0 at `collision − ridge_rise` to 1 at
`collision + ridge_rise`, and 0 throughout with `collision = Inf`.
"""
function ridge_fraction(m::BubbleAbsorption, t)
    s = (t - m.collision + m.ridge_rise) / (2 * m.ridge_rise)
    s <= 0 && return zero(s)
    s >= 1 && return one(s)
    return s^2 * (3 - 2s)
end

"""
    bubble_transmission(m::BubbleAbsorption, t) -> T_b

Transmission of the bubble for the inner-beam light it clips, at time `t` (ns); see
[`BubbleAbsorption`](@ref). 1 before the launch, 0 once the tip reaches the axis.
"""
function bubble_transmission(m::BubbleAbsorption, t)
    x = bubble_depth(m, t)
    x <= 0 && return one(x)
    R = m.radius
    x >= R && return zero(x)
    return exp(-m.opacity * (m.fill / 0.3)^m.density_exponent * (x / m.x_ref)^m.depth_exponent *
               (R / (R - x))^2)
end

# Hemisphere of a beam: +1 for ports above the equator.
_hemisphere(b::Beam) = b.port[3] >= 0 ? 1 : -1

# Returns the height at which the axis of a beam reaches the cylinder of radius R, at the
# far crossing.
function _axis_wall_z(b::Beam, R)
    k = direction(b)
    a = b.aim
    kk = k[1]^2 + k[2]^2
    ak = a[1] * k[1] + a[2] * k[2]
    aa = a[1]^2 + a[2]^2
    s = (-ak + sqrt(max(ak^2 - kk * (aa - R^2), zero(ak)))) / kk
    return a[3] + s * k[3]
end

# Returns the footprint plane height in hemisphere h (signed), or nothing without
# footprint beams.
function _footprint_z(m::BubbleAbsorption, beams, h, t)
    fz, wall = m.footprint_z, m.wall          # locals, since Julia 1.10 narrows these and not fields
    fz !== nothing && return h * fz
    R = wall === nothing ? m.radius : wall_radius(wall, m.radius, t)
    zs = [_axis_wall_z(b, R) for b in beams if b.group in m.footprint && _hemisphere(b) == h]
    isempty(zs) && return nothing
    return sum(zs) / length(zs)
end

# Returns the equal-power radii (core 1/e units) of a super-Gaussian of order n, the
# midpoints in the enclosed fraction P(2/n, ρⁿ).
_equal_power_radii(n, nr) =
    [gamma_inc_inv(2 / n, (i - 0.5) / nr, 1 - (i - 0.5) / nr)^(1 / n) for i in 1:nr]

@inline _logistic(x) = x > 0 ? 1 / (1 + exp(-x)) : exp(x) / (1 + exp(x))

"""
    iris_fraction(beam, z, ρ; edge = 0.05, nr = 12, nangle = 24, ncone = 13)

Fraction of the power of `beam` that crosses the plane at height `z` (mm) outside the
cylinder of radius `ρ` about the hohlraum axis. The best-focus spot of the beam (core
and halo) is carried along `ncone` f-cone directions, and the edge of the cylinder is
smoothed over `edge` mm. Generic in the number type of the beam, `z` and `ρ`.
"""
function iris_fraction(b::Beam, z, ρ; edge = 0.05, nr::Integer = 12, nangle::Integer = 24,
                       ncone::Integer = 13)
    k = direction(b)
    e1, e2 = spot_axes(b)
    focus = b.aim + b.defocus * k
    ds, wd = cone_directions(b, ncone)
    sp = b.spot
    cs = [cos(2π * (j - 0.5) / nangle) for j in 1:nangle]
    sn = [sin(2π * (j - 0.5) / nangle) for j in 1:nangle]
    acc = _iris_component(focus, e1, e2, ds, wd, sp.a, sp.b, _equal_power_radii(sp.order, nr),
                          cs, sn, z, ρ, edge) * (1 - sp.halo)
    if sp.halo > 0
        sc = sp.halo_scale
        acc += sp.halo * _iris_component(focus, e1, e2, ds, wd, sc * sp.a, sc * sp.b,
                                         _equal_power_radii(sp.halo_order, nr), cs, sn, z, ρ, edge)
    end
    return acc
end

function _iris_component(focus, e1, e2, ds, wd, a, b, radii, cs, sn, z, ρ, edge)
    acc = zero(promote_type(eltype(focus), typeof(a), typeof(z), typeof(ρ)))
    w = 1 / (length(radii) * length(cs))
    for (d, wj) in zip(ds, wd)
        sub = zero(acc)
        for r in radii, j in eachindex(cs)
            p = focus + (a * r * cs[j]) * e1 + (b * r * sn[j]) * e2
            s = (z - p[3]) / d[3]
            qx, qy = p[1] + s * d[1], p[2] + s * d[2]
            sub += _logistic((sqrt(qx^2 + qy^2) - ρ) / edge)
        end
        acc += wj * w * sub
    end
    return acc
end

"""
    bubble_fractions(m::BubbleAbsorption, beams, t) -> NamedTuple

The bubble at time `t` (ns) and its effect on `beams`. The fields are `depth` ``x_{bub}``
(mm), `transmission` ``T_b``, `footprint_z` (``z_b`` of the upper and lower hemispheres,
`nothing` where there are no footprint beams), `ridge` (the [`ridge_fraction`](@ref) g),
and per beam `clipped` (f, zero for beams that are not inner) and `absorbed`
(``f (1 - T_b)``). Of the power of an inner beam, `absorbed` goes to the footprint beams,
and g(1 − absorbed) stops at the ridge.
"""
function bubble_fractions(m::BubbleAbsorption, beams::AbstractVector{Beam{S}}, t) where {S<:Real}
    x = bubble_depth(m, t)
    T = bubble_transmission(m, t)
    zb = (_footprint_z(m, beams, 1, t), _footprint_z(m, beams, -1, t))
    ρ = m.radius - x
    Tf = promote_type(typeof(x), typeof(T), S)
    clipped = zeros(Tf, length(beams))
    if x > 0
        for (i, b) in enumerate(beams)
            b.group in m.inner || continue
            z = zb[_hemisphere(b) == 1 ? 1 : 2]
            z === nothing && continue
            clipped[i] = iris_fraction(b, z, ρ; edge = m.edge, nr = m.nr, nangle = m.nangle,
                                       ncone = m.ncone)
        end
    end
    return (depth = x, transmission = T, footprint_z = zb, ridge = ridge_fraction(m, t),
            clipped = clipped, absorbed = clipped .* (1 - T))
end
# Beams of mixed number types are first converted to one type, which the quantities of
# the bubble then carry.
bubble_fractions(m::BubbleAbsorption, beams::AbstractVector{<:Beam}, t) =
    bubble_fractions(m, _common_type(beams), t)

function apply_lpi(m::BubbleAbsorption, beams::AbstractVector{Beam{S}}, time = 0.0) where {S<:Real}
    # (The named tuple of the fractions is not concrete, since `footprint_z` may be
    # `nothing`, so the closures below capture only its concrete field.)
    absorbed = bubble_fractions(m, beams, time).absorbed
    P = [b.power * one(eltype(absorbed)) for b in beams]
    for h in (1, -1)
        rec = [i for (i, b) in enumerate(beams) if b.group in m.footprint && _hemisphere(b) == h]
        isempty(rec) && continue
        moved = zero(eltype(P))
        for (i, b) in enumerate(beams)
            if b.group in m.inner && _hemisphere(b) == h
                δ = b.power * absorbed[i]
                P[i] -= δ
                moved += δ
            end
        end
        for i in rec
            P[i] += moved / length(rec)
        end
    end
    return [remake(b; power = P[i]) for (i, b) in enumerate(beams)]
end
apply_lpi(m::BubbleAbsorption, beams::AbstractVector{<:Beam}, time = 0.0) =
    apply_lpi(m, _common_type(beams), time)

_deposits(::BubbleAbsorption) = true

# The ridge stop. The power of the inner beams times g(t) goes to a band of the ridge
# components at the height where the axis of each beam crosses the ridge radius. The rest
# of every beam is deposited by the operator. Power is conserved and booked per beam.
function _lpi_deposit(m::BubbleAbsorption, op::DepositionOperator,
                      beams::AbstractVector{Beam{S}}, time, split) where {S<:Real}
    g = ridge_fraction(m, time)
    g == 0 && return _deposit(op, beams, split)
    isinner = [b.group in m.inner for b in beams]
    dep = _deposit(op, [remake(b; power = b.power * (isinner[i] ? 1 - g : one(g)))
                        for (i, b) in enumerate(beams)], split)
    es = _ridge_elements(m, op.target)
    T = promote_type(eltype(dep.element), typeof(g), S)
    E, L, PB = Vector{T}(dep.element), Vector{T}(dep.backside), Vector{T}(dep.per_beam)
    added = zero(T)
    for (i, b) in enumerate(beams)
        isinner[i] || continue
        P = b.power * g
        δ = P .* _ridge_weights(m, op.target, es, b)
        E[es] .+= δ
        split === nothing || split[1][i] == 0 || (split[2][es, split[1][i]] .+= ForwardDiff.value.(δ))
        PB[i] += P
        added += P
    end
    return Deposition(E, L, PB, convert(T, dep.delivered + added))
end
_lpi_deposit(m::BubbleAbsorption, op::DepositionOperator, beams::AbstractVector{<:Beam}, time, split) =
    _lpi_deposit(m, op, _common_type(beams), time, split)

function _ridge_elements(m::BubbleAbsorption, t)
    es = reduce(vcat, (t.comp_elems[ci] for (ci, c) in enumerate(t.components)
                       if c.name in m.ridge_components || c.part in m.ridge_components); init = Int[])
    isempty(es) && throw(ArgumentError("no ridge components $(m.ridge_components) in the target"))
    return es
end

# Returns the weights (sum 1) of the ridge band of beam b over the ridge elements `es`.
function _ridge_weights(m::BubbleAbsorption, t, es, b::Beam)
    zr = _axis_wall_z(b, m.ridge * m.radius)
    w = [t.elem_area[e] * exp(-((t.elem_centroid[e][3] - zr) / m.ridge_width)^2 / 2) for e in es]
    W = sum(w)
    W > 0 || throw(ArgumentError("ridge band of $(b.name) misses the ridge components"))
    return w ./ W
end

function _ridge_split(m::BubbleAbsorption, op::DepositionOperator, b::Beam, time)
    g = ridge_fraction(m, time)
    (g == 0 || !(b.group in m.inner)) && return nothing
    es = _ridge_elements(m, op.target)
    return (; g, es, w = _ridge_weights(m, op.target, es, b))
end
