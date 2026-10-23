# Inner-beam absorption in the wall bubble: a laser–plasma model that moves part of the
# inner cones' power to the outer cones' footprint as the Au bubble grows into their path.
#
# Izumi et al., Phys. Plasmas 28, 022706 (2021), Eqs. 3 and 6, with the clipping taken
# from Briareus's own beams instead of an annulus. See `reference/research-notes.md`,
# "Inner-beam bubble absorption".

# Izumi 2021's transparency anchors (bubble transparency at 5.5 ns, computed there for a
# 5.75 mm hohlraum, R = 2.88 mm): 29 % at 0.3 mg/cc, 54 % at 0.15 mg/cc.
const _BUBBLE_ANCHORS = (radius = 2.88, time = 5.5, launch = 0.5, fills = (0.3, 0.15),
                         transmission = (0.29, 0.54))

# Opacity τ_b and density exponent q that put T_b through both anchors, with the
# default speed (0.35 mm/ns at 0.3 mg/cc, ∝ ρ₀^-0.2), x_ref = 1 mm and p = 0.5.
function _bubble_anchor_fit(; speed = 0.35, exponent = 0.5, x_ref = 1.0, fill_ref = 0.3,
                            a = _BUBBLE_ANCHORS)
    R = a.radius
    x = [speed * (ρ / fill_ref)^-0.2 * (a.time - a.launch) for ρ in a.fills]
    g = [(xi / x_ref)^exponent * (R / (R - xi))^2 for xi in x]
    τ = -log.(collect(a.transmission))
    q = log(τ[1] / τ[2] * g[2] / g[1]) / log(a.fills[1] / a.fills[2])
    τb = τ[1] / ((a.fills[1] / fill_ref)^q * g[1])
    return τb, q
end

const _BUBBLE_FIT = _bubble_anchor_fit()

"""
    BubbleAbsorption(; radius, turn_on = 0.0, delay = 0.2, speed = 0.35, fill = 0.3,
                     energy = 1.0, opacity = 0.144, density_exponent = 1.86,
                     exponent = 0.5, x_ref = 1.0, inner = (:inner23, :inner30),
                     footprint = (:outer50,), wall = nothing, footprint_z = nothing,
                     edge = 0.05, nr = 12, nangle = 24, ncone = 13)

Inner-beam absorption in the wall bubble (a [`LaserPlasmaModel`](@ref)). The outer
cones' picket launches a bubble of wall plasma from the wall under their footprint. Its
tip has moved

    x_bub(t) = v_b · max(0, t − t_L),   t_L = turn_on + delay,
    v_b = speed · energy^0.2 · (fill / 0.3)^-0.2

in from the initial wall `radius` R (mm, ns; `energy` is the outer picket fluence over
the reference design's, `fill` the gas fill in mg/cc). In the plane `z_b` of the outer
footprint, at the axial position where the footprint cones' beam axes reach the wall,
the bubble fills the annulus from R in to radius R − x_bub. Each inner beam loses

    absorbed = f · (1 − T_b),
    T_b = exp(−τ_b (fill/0.3)^q (x_bub/x_ref)^p (R/(R − x_bub))²),

where `f` is the fraction of its power crossing `z_b` inside the bubble, computed from
the beam itself (its best-focus spot carried along its `ncone` f-cone directions, as
deposition does), and `T_b` the bubble's transmission (`opacity` τ_b,
`density_exponent` q, `exponent` p). The absorbed power is conducted to the wall
under the outer footprint (Ralph et al. 2018): it is shared equally by the `footprint`
beams of the inner beam's hemisphere (those whose port is on the same side), so it is
deposited where those beams land. Beam powers are conserved: the model moves power
between beams and loses none. (Whatever fraction of the footprint beams themselves
misses the wall, a few 1e-4 through the entrance holes, applies to the moved power too.)

The edge of the bubble is smoothed over `edge` (mm, a logistic step), which keeps `f`
smooth in time, in the bubble depth and in the beams' pointing. Spot points: `nr`
equal-power radii by `nangle` angles per f-cone direction. Everything is generic in the
number type: `speed` and `opacity` (and the beams) can carry ForwardDiff duals, so the
model enters [`drive_history`](@ref)`(...; lpi)`, `Scenario(...; lpi)` and
calibrations. Compose with CBET as `lpi = (cbet, bubble)` (transfer first, then the
bubble takes its share of what the inner beams carry). `wall`, a [`WallMotion`](@ref),
moves the footprint plane with the radiating wall (the outer beams' axes reach the
moved wall nearer the entrance hole); the bubble tip is measured from the initial wall
either way. `footprint_z` fixes |z_b| instead.

Engineering defaults (research notes, "Inner-beam bubble absorption"):
- `speed` 0.35 mm/ns at 0.3 mg/cc and the reference picket, ∝ E^0.2 ρ₀^-0.2: Izumi
  et al., PoP 28, 022706 (2021) (345 µm/ns baseline); Meezan et al. 2020 measured
  321–461 µm/ns. `delay` 0.2 ns: Izumi's launch at the outer turn-on, fit intercept
  0.5 ± 0.1 ns for a 0.3 ns turn-on.
- `opacity` and `density_exponent` are *fitted* so that T_b passes through Izumi
  2021's computed transparencies at 5.5 ns, 29 % (0.3 mg/cc) and 54 % (0.15), in a
  2.88 mm-radius hohlraum (`_bubble_anchor_fit`). The fitted q ≈ 1.86 is near the
  ρ² of inverse-bremsstrahlung absorption. `exponent` 0.5 is Izumi's Eq. 6
  (x_bub/1 mm)^0.5; the (R/(R − x))² factor stands for the compressed bubble's density
  (Izumi's ρ_bub ∝ (R_hohl/R_bub)^k with k = 1; the paper's exponent is not legible).
- The tip moves ballistically. Izumi 2023 finds it stagnates in long low-trough pulses
  and slows near the axis; late in a long peak the model clips the inner beams
  entirely.
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
    inner::Vector{Symbol}
    footprint::Vector{Symbol}
    wall::Union{Nothing,WallMotion}
    footprint_z::Union{Nothing,Float64}
    edge::Float64
    nr::Int
    nangle::Int
    ncone::Int
end

function BubbleAbsorption(; radius::Real, turn_on::Real = 0.0, delay::Real = 0.2,
                          speed::Real = 0.35, opacity::Real = _BUBBLE_FIT[1], fill::Real = 0.3,
                          energy::Real = 1.0, density_exponent::Real = _BUBBLE_FIT[2],
                          exponent::Real = 0.5, x_ref::Real = 1.0,
                          inner = (:inner23, :inner30), footprint = (:outer50,),
                          wall::Union{Nothing,WallMotion} = nothing, footprint_z = nothing,
                          edge::Real = 0.05, nr::Integer = 12, nangle::Integer = 24,
                          ncone::Integer = 13)
    radius > 0 || throw(ArgumentError("radius must be positive"))
    fill > 0 && energy > 0 || throw(ArgumentError("fill and energy must be positive"))
    edge > 0 || throw(ArgumentError("edge must be positive"))
    nr >= 1 && nangle >= 1 && 1 <= ncone <= 16 || throw(ArgumentError("bad sampling counts"))
    s, τ = promote(speed, opacity)
    return BubbleAbsorption(Float64(radius), Float64(turn_on + delay), s, τ, Float64(fill),
                            Float64(energy), Float64(density_exponent), Float64(exponent),
                            Float64(x_ref), Symbol[inner...], Symbol[footprint...], wall,
                            footprint_z === nothing ? nothing : Float64(abs(footprint_z)),
                            Float64(edge), Int(nr), Int(nangle), Int(ncone))
end

"""
    remake(m::BubbleAbsorption; speed = m.speed, opacity = m.opacity)

The model with new calibration parameters (e.g. ForwardDiff duals in a fit).
"""
function remake(m::BubbleAbsorption; speed = m.speed, opacity = m.opacity)
    s, τ = promote(speed, opacity)
    return BubbleAbsorption(m.radius, m.launch, s, τ, m.fill, m.energy, m.density_exponent,
                            m.exponent, m.x_ref, m.inner, m.footprint, m.wall, m.footprint_z,
                            m.edge, m.nr, m.nangle, m.ncone)
end

function Base.show(io::IO, m::BubbleAbsorption)
    @printf(io, "BubbleAbsorption(R = %.3g mm, v_b = %.3g mm/ns from %.3g ns, τ_b = %.3g, %s → %s)",
            m.radius, ForwardDiff.value(bubble_speed(m)), m.launch, ForwardDiff.value(m.opacity),
            join(m.inner, "/"), join(m.footprint, "/"))
end

"Bubble tip speed (mm/ns): `speed · energy^0.2 · (fill/0.3)^-0.2`."
bubble_speed(m::BubbleAbsorption) = m.speed * m.energy^0.2 * (m.fill / 0.3)^-0.2

"Bubble tip displacement from the initial wall (mm) at time `t` (ns)."
bubble_depth(m::BubbleAbsorption, t) = bubble_speed(m) * max(zero(t), t - m.launch)

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

The bubble at time `t` (ns) and what it does to `beams`: `depth` x_bub (mm),
`transmission` T_b, `footprint_z` (z_b of the upper and lower hemispheres, `nothing`
where there are no footprint beams), and per beam `clipped` (f, zero for beams that
are not inner) and `absorbed` (f·(1 − T_b)).
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
    return (depth = x, transmission = T, footprint_z = zb, clipped = clipped,
            absorbed = clipped .* (1 - T))
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
