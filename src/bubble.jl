# Inner-beam absorption in the wall bubble: a laser–plasma model that moves part of the
# inner cones' power to the outer cones' footprint as the Au bubble grows into their path.
#
# Izumi et al., Phys. Plasmas 28, 022706 (2021), Eqs. 3 and 6, with the clipping taken
# from Briareus's own beams instead of an annulus. See `reference/research-notes.md`,
# "Inner-beam bubble absorption".

# Calibration (research notes, "Hohlraum dynamics calibrated to published data
# (2026-10-24)"): tip speed and launch fitted to Izumi 2021 Fig. 8a (28 measured points,
# seven shots, `izumi2021_bubble()`), the opacity to the HYDRA absorbed fraction of Fig. 10a
# (`izumi2021_absorption()`) on the basecamp geometry with the bubble at z = ±3 mm.
const _BUBBLE_SPEED = 0.36      # mm/ns at 0.3 mg/cc and the basecamp picket (0.347–0.372)
const _BUBBLE_DELAY = 0.25      # ns after the outer turn-on (0.07–0.42)
const _BUBBLE_OPACITY = 0.111   # τ_b (0.100–0.124 statistical; 0.085–0.158 with the plane and launch)
const _BUBBLE_Q = 1.8           # ρ_bub² ∝ (S ρ₀)², S ∝ ρ₀^-0.1 (Izumi 2021 Eq. 7)

# Basecamp geometry of Izumi 2018/2021: hohlraum and capsule outer radius (mm), and the
# ablator–wall collision 5.5 ns after the outer turn-on (inner spots separate at 5.8 ns).
const _BASECAMP = (radius = 2.875, capsule = 0.908, collision = 5.5)

"""
    BubbleAbsorption(; radius, turn_on = 0.0, delay = 0.25, speed = 0.36, fill = 0.3,
                     energy = 1.0, opacity = 0.111, density_exponent = 1.8,
                     exponent = 0.5, x_ref = 1.0, stagnation = 0.80 radius,
                     stagnation_width = 0.08 radius, collision = <t_c>, capsule_radius = 0.908,
                     ridge = 0.73, ridge_width = 0.5, ridge_rise = 0.2,
                     inner = (:inner23, :inner30), footprint = (:outer50,),
                     ridge_components = (:wall,), wall = nothing, footprint_z = nothing,
                     edge = 0.05, nr = 12, nangle = 24, ncone = 13)

Inner-beam absorption in the wall bubble (a [`LaserPlasmaModel`](@ref)). The outer
cones' picket launches a bubble of wall plasma from the wall under their footprint. Its
tip has moved

    x_bal(t) = v_b · max(0, t − t_L),   t_L = turn_on + delay,
    v_b = speed · energy^0.2 · (fill / 0.3)^-0.2

in from the initial wall `radius` R (mm, ns; `energy` is the outer picket fluence over
the reference design's, `fill` the gas fill in mg/cc) while it is ballistic. It
stagnates smoothly at the depth `stagnation` ``x_{max}`` ([`bubble_depth`](@ref)):
ballistic up to ``x_{max} - w`` (w = `stagnation_width`), then
``x_{max} - w e^{-(x_{bal} - x_{max} + w)/w}``. In the plane `z_b` of the outer
footprint, at the axial position where the footprint cones' beam axes reach the wall,
the bubble fills the annulus from R in to radius R − ``x_{bub}``. Each inner beam loses

    absorbed = f · (1 − T_b),
    T_b = exp(−τ_b (fill/0.3)^q (x_bub/x_ref)^p (R/(R − x_bub))²),

where `f` is the fraction of its power crossing `z_b` inside the bubble, computed from
the beam itself (its best-focus spot carried along its `ncone` f-cone directions, as
deposition does), and `T_b` the bubble's transmission (`opacity` ``τ_b``,
`density_exponent` q, `exponent` p). The absorbed power is conducted to the wall
under the outer footprint (Ralph et al. 2018): it is shared equally by the `footprint`
beams of the inner beam's hemisphere (those whose port is on the same side), so it is
deposited where those beams land. Beam powers are conserved: the model moves power
between beams and loses none. (Whatever fraction of the footprint beams themselves
misses the wall, a few 1e-4 through the entrance holes, applies to the moved power too.)

**Ridge stop.** After the ablator–wall collision at `collision` ``t_c``, the inner power
the bubble leaves stops at the stagnation ridge, radius `ridge`·R about the axis
(Izumi 2023), instead of reaching the wall. A fraction g(t) of it ([`ridge_fraction`](@ref),
rising over ±`ridge_rise` ns about ``t_c``) is deposited, by [`lpi_deposit`](@ref), on an
axisymmetric band of the `ridge_components` (the hohlraum wall) centered on the height
where the beam's axis crosses that radius, Gaussian in z with rms `ridge_width` (mm),
area-weighted: the heated ridge radiates at the wall's conversion, from the wall just
outside it. The rest of the beam's power is deposited as usual. Power is conserved
(the ridge power is booked to the beam's `per_beam`). The ridge stop acts only where
the deposition goes through `lpi_deposit`: [`drive_history`](@ref), `Scenario`
evaluations, the beam-error sensitivities ([`mode_sensitivity`](@ref) and the rest of
`uncertainty.jl`) and [`influence_matrix`](@ref).

The edge of the bubble is smoothed over `edge` (mm, a logistic step), which keeps `f`
smooth in time, in the bubble depth and in the beams' pointing. Spot points: `nr`
equal-power radii by `nangle` angles per f-cone direction. Everything is generic in the
number type: `speed`, `opacity`, `stagnation` and `collision` (and the beams) can carry
ForwardDiff duals ([`remake`](@ref)), so the model enters
[`drive_history`](@ref)`(...; lpi)`, `Scenario(...; lpi)` and calibrations. Compose with
CBET as `lpi = (cbet, bubble)` (transfer first, then the bubble takes its share of what
the inner beams carry). `wall`, a [`WallMotion`](@ref), moves the footprint plane with
the radiating wall (the outer beams' axes reach the moved wall nearer the entrance
hole); the bubble tip and the ridge are measured from the initial wall either way.
`footprint_z` fixes |``z_b``| instead. `stagnation = Inf` keeps the tip ballistic;
`collision = Inf` turns the ridge stop off.

Defaults (research notes, "Hohlraum dynamics calibrated to published data (2026-10-24)";
the data are [`izumi2021_bubble`](@ref) and [`izumi2021_absorption`](@ref)):
- `speed` 0.36 mm/ns at 0.3 mg/cc and the basecamp picket, ∝ E^0.2 ρ₀^-0.2, and
  `delay` 0.25 ns after the outer turn-on: *fitted* to the 28 measured tip positions of
  Izumi et al., PoP 28, 022706 (2021) Fig. 8a (seven HDC SymCap shots, 0.15–0.3 mg/cc,
  picket 0.73–1.5×; χ²/dof 0.89 with the paper's ±40 µm): 0.360 (0.347–0.372) mm/ns,
  delay 0.25 (0.07–0.42) ns, the two correlated. Meezan et al. 2020 measured
  321–461 µm/ns on the same shots.
- `density_exponent` q = 1.8: inverse-bremsstrahlung absorption goes as ``ρ_{bub}^2``, and
  Izumi's Eq. 7 has ``ρ_{bub}`` = S ρ₀ R/(R − x) with the shock compression S ∝ ρ₀^-0.1,
  so ``ρ_{bub}^2`` ∝ ρ₀^1.8 (R/(R − x))². `exponent` 0.5 is Eq. 6's chord
  (``x_{bub}``/1 mm)^0.5.
- `opacity` ``τ_b`` 0.111 is *fitted* to the inner-power fraction absorbed in the bubble
  in Izumi's 2D HYDRA simulation of the 0.3 mg/cc basecamp (Fig. 10a: 0.02, 0.23, 0.55 at
  4, 5, 6 ns), with this model's own clipping of NIF beams pointed for a 5.75 mm hohlraum
  and the bubble plane at Izumi's z = ±3 mm (`footprint_z = 3.0`): absorbed 0.03, 0.20,
  0.36, 0.56 at 4, 5, 5.5, 6 ns, χ² 1.1 at ±0.03. Statistical range 0.100–0.124; the
  plane (2.8–3.2 mm) and the launch move ``τ_b`` over 0.085–0.158 but leave the absorbed
  fraction within 0.01, which is the calibrated quantity. ``T_b`` at 5.5 ns is then 0.36 at
  0.3 mg/cc and 0.58 at 0.15 (Izumi's own model: 0.29 and 0.54, with an annulus for f).
- `stagnation` 0.80 R, `stagnation_width` 0.08 R (*engineering*): the tip is exactly
  ballistic to 0.72 R (2.07 mm in a 2.88 mm hohlraum), past the measured depths at
  0.3 mg/cc (2.07 mm at 6.2 ns) and near those at 0.15 (up to 2.28 mm, 0.79 R). Izumi
  2023: the tip slows on back-pressure near the axis. (R − 0.73 R, the ridge's depth at
  the waist, would stop the tip before it reaches the inner beams at ``z_b``, 1.3 mm from
  the wall.)
- `collision` ``t_c`` = `turn_on` + 5.5 ns · (R − ``r_{cap}``)/(2.875 − 0.908 mm) ·
  energy^-0.16
  (*fitted*, Izumi et al., RSI 89, 10K111, 2018: the inner spots separate at 5.8 ns for
  the basecamp's 0.3 ns turn-on, 5.45 ns with a 1.5× picket, and about 0.4 ns earlier
  with the 1.18× capsule; the gap between the capsule's outer radius `capsule_radius`
  ``r_{cap}`` and the wall gives 0.44 ns for that capsule, so a constant closing speed of
  ablator and wall plasma, 0.36 mm/ns, explains the capsule shot). `ridge` 0.73: Izumi
  2023's ridge at r ≈ 2.1 mm in a 2.88 mm hohlraum. `ridge_width` 0.5 mm and
  `ridge_rise` 0.2 ns are engineering choices (not measured).
"""
struct BubbleAbsorption{S<:Real} <: LaserPlasmaModel
    radius::Float64
    launch::Float64
    speed::S
    opacity::S
    fill::Float64
    energy::Float64
    density_exponent::Float64
    exponent::Float64
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

# Ablator–wall collision time (Izumi 2018): the basecamp's 5.5 ns after the outer turn-on,
# scaled with the capsule–wall gap and the picket energy.
_collision_time(turn_on, radius, capsule_radius, energy) =
    turn_on + _BASECAMP.collision * (radius - capsule_radius) /
              (_BASECAMP.radius - _BASECAMP.capsule) * energy^-0.16

function BubbleAbsorption(; radius::Real, turn_on::Real = 0.0, delay::Real = _BUBBLE_DELAY,
                          speed::Real = _BUBBLE_SPEED, opacity::Real = _BUBBLE_OPACITY,
                          fill::Real = 0.3, energy::Real = 1.0,
                          density_exponent::Real = _BUBBLE_Q,
                          exponent::Real = 0.5, x_ref::Real = 1.0,
                          stagnation::Real = 0.80 * radius,
                          stagnation_width::Real = 0.08 * radius,
                          capsule_radius::Real = _BASECAMP.capsule,
                          collision::Real = _collision_time(turn_on, radius, capsule_radius, energy),
                          ridge::Real = 0.73, ridge_width::Real = 0.5, ridge_rise::Real = 0.2,
                          inner = (:inner23, :inner30), footprint = (:outer50,),
                          ridge_components = (:wall,),
                          wall::Union{Nothing,WallMotion} = nothing, footprint_z = nothing,
                          edge::Real = 0.05, nr::Integer = 12, nangle::Integer = 24,
                          ncone::Integer = 13)
    radius > 0 || throw(ArgumentError("radius must be positive"))
    fill > 0 && energy > 0 || throw(ArgumentError("fill and energy must be positive"))
    edge > 0 || throw(ArgumentError("edge must be positive"))
    stagnation_width > 0 && stagnation > stagnation_width ||
        throw(ArgumentError("stagnation must exceed stagnation_width > 0"))
    0 < ridge < 1 && ridge_width > 0 && ridge_rise > 0 ||
        throw(ArgumentError("need 0 < ridge < 1, ridge_width > 0, ridge_rise > 0"))
    nr >= 1 && nangle >= 1 && 1 <= ncone <= 16 || throw(ArgumentError("bad sampling counts"))
    s, τ, xs, tc = promote(speed, opacity, stagnation, collision)
    return BubbleAbsorption(Float64(radius), Float64(turn_on + delay), s, τ, Float64(fill),
                            Float64(energy), Float64(density_exponent), Float64(exponent),
                            Float64(x_ref), xs, Float64(stagnation_width), tc, Float64(ridge),
                            Float64(ridge_width), Float64(ridge_rise), Symbol[inner...],
                            Symbol[footprint...], Symbol[ridge_components...], wall,
                            footprint_z === nothing ? nothing : Float64(abs(footprint_z)),
                            Float64(edge), Int(nr), Int(nangle), Int(ncone))
end

"""
    remake(m::BubbleAbsorption; speed = m.speed, opacity = m.opacity,
           stagnation = m.stagnation, collision = m.collision)

The model with new calibration parameters (e.g. ForwardDiff duals in a fit).
"""
function remake(m::BubbleAbsorption; speed = m.speed, opacity = m.opacity,
                stagnation = m.stagnation, collision = m.collision)
    s, τ, xs, tc = promote(speed, opacity, stagnation, collision)
    return BubbleAbsorption(m.radius, m.launch, s, τ, m.fill, m.energy, m.density_exponent,
                            m.exponent, m.x_ref, xs, m.stagnation_width, tc, m.ridge,
                            m.ridge_width, m.ridge_rise, m.inner, m.footprint, m.ridge_components,
                            m.wall, m.footprint_z, m.edge, m.nr, m.nangle, m.ncone)
end

function Base.show(io::IO, m::BubbleAbsorption)
    @printf(io, "BubbleAbsorption(R = %.3g mm, v_b = %.3g mm/ns from %.3g ns, τ_b = %.3g, x_max = %.3g mm, ridge from %.3g ns, %s → %s)",
            m.radius, ForwardDiff.value(bubble_speed(m)), m.launch, ForwardDiff.value(m.opacity),
            ForwardDiff.value(m.stagnation), ForwardDiff.value(m.collision),
            join(m.inner, "/"), join(m.footprint, "/"))
end

"""
    izumi2021_bubble() -> NamedTuple

Au bubble tip displacement from the initial wall against time on the HDC SymCap series of
Izumi et al., PoP 28, 022706 (2021), Fig. 8 (5.75 mm Au hohlraum, outer turn-on 0.3 ns),
extracted from the figure's vector graphics (`data/izumi2021_fig8.csv`, provenance in
its header). Columns `shot`, `source` (`"data"`: measured, Fig. 8a, error bars ±40 µm;
`"hydra"`: 2D HYDRA, Fig. 8b), `energy` (outer picket fluence over the basecamp's),
`fill` (mg/cc), `capsule` (radius over the basecamp's), `t` (ns), `x_bub` (mm) and
`u` (mm). The calibration target of [`BubbleAbsorption`](@ref)'s `speed` and `delay`.
"""
izumi2021_bubble() = _data_table("izumi2021_fig8.csv")

"""
    izumi2021_absorption() -> NamedTuple

The fraction of the inner-cone power absorbed in the bubble against time for the
0.3 mg/cc basecamp of Izumi et al., PoP 28, 022706 (2021), Fig. 10a, from the figure's
vector graphics (`data/izumi2021_fig10a.csv`): `source` `"hydra"` (the 2D simulation:
0, 0.02, 0.23, 0.55 at 3, 4, 5, 6 ns) or `"model"` (the paper's Eq. 6, every 0.1 ns),
`t` (ns), `absorbed` and `u`. The calibration target of [`BubbleAbsorption`](@ref)'s
`opacity`.
"""
izumi2021_absorption() = _data_table("izumi2021_fig10a.csv")

"Bubble tip speed (mm/ns): `speed · energy^0.2 · (fill/0.3)^-0.2`."
bubble_speed(m::BubbleAbsorption) = m.speed * m.energy^0.2 * (m.fill / 0.3)^-0.2

"""
    bubble_depth(m::BubbleAbsorption, t)

Bubble tip displacement from the initial wall (mm) at time `t` (ns): the ballistic
depth ``v_b (t - t_L)`` up to the knee `stagnation − stagnation_width`, beyond it
`stagnation − w·exp(−(x − knee)/w)`, which stops at `stagnation` (C¹ at the knee).
"""
function bubble_depth(m::BubbleAbsorption, t)
    x = bubble_speed(m) * max(zero(t), t - m.launch)
    return _stagnate(x, m.stagnation, m.stagnation_width)
end

function _stagnate(x, xmax, w)
    isinf(xmax) && return x + zero(xmax)
    k = xmax - w
    return x <= k ? x + zero(k) : xmax - w * exp(-(x - k) / w)
end

"""
    ridge_fraction(m::BubbleAbsorption, t) -> g

The fraction of the inner beams' remaining power that stops at the ablator–wall ridge at
time `t` (ns): a C¹ smoothstep from 0 at `collision − ridge_rise` to 1 at
`collision + ridge_rise` (0 throughout with `collision = Inf`).
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
    return exp(-m.opacity * (m.fill / 0.3)^m.density_exponent * (x / m.x_ref)^m.exponent *
               (R / (R - x))^2)
end

# Hemisphere of a beam: +1 for ports above the equator.
_hemisphere(b::Beam) = b.port[3] >= 0 ? 1 : -1

# Height at which a beam's axis reaches the cylinder of radius R (the far crossing).
function _axis_wall_z(b::Beam, R)
    k = direction(b)
    a = b.aim
    kk = k[1]^2 + k[2]^2
    ak = a[1] * k[1] + a[2] * k[2]
    aa = a[1]^2 + a[2]^2
    s = (-ak + sqrt(max(ak^2 - kk * (aa - R^2), zero(ak)))) / kk
    return a[3] + s * k[3]
end

# Footprint plane height in hemisphere h (signed), or nothing without footprint beams.
function _footprint_z(m::BubbleAbsorption, beams, h, t)
    m.footprint_z !== nothing && return h * m.footprint_z
    R = m.wall === nothing ? m.radius : wall_radius(m.wall, m.radius, t)
    zs = [_axis_wall_z(b, R) for b in beams if b.group in m.footprint && _hemisphere(b) == h]
    isempty(zs) && return nothing
    return sum(zs) / length(zs)
end

# Equal-power radii (core 1/e units) of a super-Gaussian of order n: midpoints in the
# enclosed fraction P(2/n, ρⁿ).
_equal_power_radii(n, nr) =
    [gamma_inc_inv(2 / n, (i - 0.5) / nr, 1 - (i - 0.5) / nr)^(1 / n) for i in 1:nr]

@inline _logistic(x) = x > 0 ? 1 / (1 + exp(-x)) : exp(x) / (1 + exp(x))

"""
    iris_fraction(beam, z, ρ; edge = 0.05, nr = 12, nangle = 24, ncone = 13)

Fraction of `beam`'s power that crosses the plane at height `z` (mm) outside the
cylinder of radius `ρ` about the hohlraum axis: its best-focus spot (core and halo)
carried along `ncone` f-cone directions, with the cylinder's edge smoothed over `edge`
mm. Generic in the number type of the beam, `z` and `ρ`.
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

The bubble at time `t` (ns) and what it does to `beams`: `depth` ``x_{bub}`` (mm),
`transmission` ``T_b``, `footprint_z` (``z_b`` of the upper and lower hemispheres, `nothing`
where there are no footprint beams), `ridge` (the [`ridge_fraction`](@ref) g), and per
beam `clipped` (f, zero for beams that are not inner) and `absorbed` (``f (1 - T_b)``).
Of an inner beam's power, `absorbed` goes to the footprint beams, and g(1 − absorbed)
stops at the ridge.
"""
function bubble_fractions(m::BubbleAbsorption, beams::AbstractVector{<:Beam}, t)
    x = bubble_depth(m, t)
    T = bubble_transmission(m, t)
    zb = (_footprint_z(m, beams, 1, t), _footprint_z(m, beams, -1, t))
    ρ = m.radius - x
    Tf = promote_type(typeof(x), typeof(T), map(b -> typeof(b.power), beams)...,
                      map(b -> eltype(b.aim), beams)...)
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

function apply_lpi(m::BubbleAbsorption, beams::AbstractVector{<:Beam}, time = 0.0)
    fr = bubble_fractions(m, beams, time)
    P = [b.power * one(eltype(fr.absorbed)) for b in beams]
    for h in (1, -1)
        rec = [i for (i, b) in enumerate(beams) if b.group in m.footprint && _hemisphere(b) == h]
        isempty(rec) && continue
        moved = zero(eltype(P))
        for (i, b) in enumerate(beams)
            if b.group in m.inner && _hemisphere(b) == h
                δ = b.power * fr.absorbed[i]
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

_deposits(::BubbleAbsorption) = true

# The ridge stop: the inner beams' power times g(t) goes to a band of the ridge
# components at the height where each beam's axis crosses the ridge radius; the rest
# of every beam is deposited by the operator. Power is conserved and booked per beam.
function _lpi_deposit(m::BubbleAbsorption, op::DepositionOperator,
                      beams::AbstractVector{<:Beam}, time, split)
    g = ridge_fraction(m, time)
    g == 0 && return _deposit(op, beams, split)
    isinner = [b.group in m.inner for b in beams]
    dep = _deposit(op, [remake(b; power = b.power * (isinner[i] ? 1 - g : one(g)))
                        for (i, b) in enumerate(beams)], split)
    es = _ridge_elements(m, op.target)
    T = promote_type(eltype(dep.element), typeof(g), map(b -> eltype(b.aim), beams)...)
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

function _ridge_elements(m::BubbleAbsorption, t)
    es = reduce(vcat, (t.comp_elems[ci] for (ci, c) in enumerate(t.components)
                       if c.name in m.ridge_components || c.part in m.ridge_components); init = Int[])
    isempty(es) && throw(ArgumentError("no ridge components $(m.ridge_components) in the target"))
    return es
end

# Weights (sum 1) of beam b's ridge band over the ridge elements `es`.
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
