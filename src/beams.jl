# Laser beams: port geometry, focal-spot profiles, pointing, and defocus.
#
# Conventions
# * A beam enters from its port, a unit vector `port` from target-chamber center (TCC)
#   toward the final optics; it propagates along `k = -port`. Pointing translates the
#   beam axis to pass through `aim`; the direction stays fixed (the final optics are
#   meters away, so mm-scale repointing changes the direction by ≲1e-4 rad).
# * The best-focus plane sits at `aim + defocus * k`: positive defocus places focus
#   beyond the aim point, so the beam is still converging there.
# * Spot axes: `a` lies along `cos(ψ) θ̂ + sin(ψ) φ̂` of the port, `b` perpendicular,
#   where θ̂ and φ̂ are the spherical unit vectors at the port.

"""
    SuperGaussianSpot(a, b = a; order = 2.0, angle = 0.0, halo = 0.0, halo_scale = 3.0,
                      halo_order = 2.0)

Elliptical super-Gaussian best-focus intensity `I ∝ exp(-[(u/a)² + (v/b)²]^(order/2))`
with 1/e semi-axes `a`, `b` (mm). `angle` (radians) rotates the `a` axis from the
port's θ̂ toward φ̂.

Measured phase-plate spots have low-level wings beyond the super-Gaussian core fit.
A fraction `halo` of the power can be put in a second, broader component with semi-axes
`halo_scale·(a, b)` and order `halo_order` — the halo matters for clearance and
clipping checks, not for drive symmetry.
"""
struct SuperGaussianSpot{T<:Real}
    a::T
    b::T
    order::Float64
    angle::Float64
    halo::Float64
    halo_scale::Float64
    halo_order::Float64
end

function SuperGaussianSpot(a::Real, b::Real = a; order = 2.0, angle = 0.0, halo = 0.0,
                           halo_scale = 3.0, halo_order = 2.0)
    0 <= halo < 1 || throw(ArgumentError("halo fraction must be in [0, 1)"))
    halo > 0 && halo_scale <= 0 && throw(ArgumentError("halo_scale must be positive"))
    return SuperGaussianSpot(promote(a, b)..., Float64(order), Float64(angle), Float64(halo),
                             Float64(halo_scale), Float64(halo_order))
end

"Same spot profile with semi-axes `a`, `b` (number type follows them)."
with_axes(s::SuperGaussianSpot, a, b) =
    SuperGaussianSpot(promote(a, b)..., s.order, s.angle, s.halo, s.halo_scale, s.halo_order)

_retype(s::SuperGaussianSpot, ::Type{T}) where {T} =
    SuperGaussianSpot{T}(s.a, s.b, s.order, s.angle, s.halo, s.halo_scale, s.halo_order)

# ∫ exp(-ρ^n) dA over the unit-semi-axis ellipse: 2π Γ(2/n) / n
_sg_norm(n) = 2π * gamma(2 / n) / n
# <u²>/a² for the super-Gaussian: Γ(4/n) / (2 Γ(2/n))
_sg_moment(n) = gamma(4 / n) / (2 * gamma(2 / n))

"""
    Beam(name; theta, phi, power = 1.0, spot, fnumber = 20.0, aim = (0,0,0),
         defocus = 0.0, group = :none, quad = "")

One laser beam entering from the port at polar angle `theta` and azimuth `phi`
(degrees, target-chamber frame, z up). `power` in TW. `group` labels beams that share
pointing parameters (e.g. a cone). `wavelength` is the nominal vacuum wavelength (µm)
and `dlambda` the detuning from it (Å), used by laser–plasma models.
"""
struct Beam{T<:Real}
    name::String
    port::Vec3
    aim::SVector{3,T}
    power::T
    defocus::T
    spot::SuperGaussianSpot{T}
    fnumber::Float64
    group::Symbol
    quad::String
    wavelength::Float64       # nominal vacuum wavelength (µm)
    dlambda::T                # wavelength offset from nominal (Å)
end

function Beam(name::AbstractString; theta, phi, power = 1.0, spot = SuperGaussianSpot(0.5),
              fnumber = 20.0, aim = (0.0, 0.0, 0.0), defocus = 0.0, group = :none, quad = "",
              wavelength = 0.351, dlambda = 0.0)
    p = port_vector(theta, phi)
    T = promote_type(eltype(aim), typeof(power), typeof(defocus), typeof(spot.a), typeof(dlambda))
    sp = _retype(spot, T)
    return Beam{T}(String(name), p, SVector{3,T}(aim...), T(power), T(defocus), sp,
                   Float64(fnumber), Symbol(group), String(quad), Float64(wavelength), T(dlambda))
end

"""
    port_vector(theta, phi)

Unit vector from TCC to a port at polar angle `theta`, azimuth `phi` (degrees).
"""
port_vector(θ, φ) = Vec3(sind(θ) * cosd(φ), sind(θ) * sind(φ), cosd(θ))

"Polar and azimuthal angles (degrees) of a beam's port."
port_angles(b::Beam) = (acosd(clamp(b.port[3], -1, 1)), mod(atand(b.port[2], b.port[1]), 360))

"Unit propagation direction of a beam, `-port`."
direction(b::Beam) = -b.port

"""
    port_basis(port) -> (θ̂, φ̂)

Spherical unit vectors at the port; the transverse frame for spot orientation.
"""
function port_basis(p::Vec3)
    θ = acos(clamp(p[3], -1, 1))
    φ = atan(p[2], p[1])
    eθ = Vec3(cos(θ) * cos(φ), cos(θ) * sin(φ), -sin(θ))
    eφ = Vec3(-sin(φ), cos(φ), 0)
    return eθ, eφ
end

"""
    spot_axes(beam) -> (a, b)

Unit vectors of the spot's `a` and `b` axes, transverse to the beam: the port's θ̂ and
φ̂ turned by the spot's `angle` toward φ̂.
"""
function spot_axes(b::Beam)
    eθ, eφ = port_basis(b.port)
    ψ = b.spot.angle
    return cos(ψ) * eθ + sin(ψ) * eφ, -sin(ψ) * eθ + cos(ψ) * eφ
end

"""
    remake(beam; aim, power, defocus, spot, dlambda, port, fnumber, group, name)

Copy of `beam` with the given fields replaced; the element type widens as needed (so
ForwardDiff duals can flow through pointing parameters).
"""
function remake(b::Beam; aim = b.aim, power = b.power, defocus = b.defocus, spot = b.spot,
                dlambda = b.dlambda, port = b.port, fnumber = b.fnumber, group = b.group,
                name = b.name)
    T = promote_type(eltype(aim), typeof(power), typeof(defocus), typeof(spot.a), typeof(spot.b),
                     typeof(dlambda))
    sp = _retype(spot, T)
    return Beam{T}(String(name), Vec3(port), SVector{3,T}(aim...), T(power), T(defocus), sp,
                   Float64(fnumber), Symbol(group), b.quad, b.wavelength, T(dlambda))
end

"Wavelength (µm) including the detuning `dlambda`."
wavelength(b::Beam) = b.wavelength + 1e-4 * b.dlambda

function Base.show(io::IO, b::Beam)
    θ, φ = port_angles(b)
    @printf(io, "Beam(%s, θ=%.2f°, φ=%.2f°, %s TW, aim=(%s), defocus=%s mm)", b.name, θ, φ,
            _fmt(b.power), join(_fmt.(b.aim), ", "), _fmt(b.defocus))
end
_fmt(x::AbstractFloat) = @sprintf("%.4g", x)
_fmt(x) = string(ForwardDiff.value(x))

# ---------------------------------------------------------------------------------------
# Intensity evaluation
# ---------------------------------------------------------------------------------------

"""
    BeamFrame(beam)

Precomputed quantities for evaluating a beam's intensity at many points.
"""
# One super-Gaussian component of a spot in a beam frame.
struct _SGComponent{T}
    a2::T
    b2::T
    halfn::Float64
    defocus_c::Float64        # (1/(2F))² / (4 m_n): adds s²·c to a², b²
    power_norm::T             # component power / (2π Γ(2/n)/n)
    qcut::Float64             # (u/a)² + (v/b)² beyond which the profile is < e⁻⁴⁰ of peak
    peak::T                   # power_norm / (a b): the best-focus peak intensity
end

function _SGComponent(a2, b2, halfn, defocus_c, power_norm)
    T = promote_type(typeof(a2), typeof(b2), typeof(power_norm))
    return _SGComponent{T}(a2, b2, halfn, defocus_c, power_norm, 40.0^(1 / halfn),
                           power_norm / sqrt(a2 * b2))
end

# Best-focus profile of one component at focal-plane coordinates (u, v); zero where it
# is below e⁻⁴⁰ of its peak, which skips the transcendental functions off the spot.
@inline function _sg_focal(c::_SGComponent, u, v)
    q = u^2 / c.a2 + v^2 / c.b2
    q > c.qcut && return zero(c.power_norm * q)
    q > 0 || return c.peak * one(q)
    return c.peak * exp(-_pow(q, c.halfn))
end

# `q^h` for q > 0 as exp(h log q): half the cost of `^` (whose extra accuracy, ~1e-15
# relative, the profile does not need), and a quarter on duals, where `^` takes a
# second power for the derivative.
@inline _pow(q, h) = exp(h * log(q))

# Whether (u, v) is inside the cutoff of the profile (`_sg_focal` is nonzero there).
@inline _sg_inside(c::_SGComponent, u, v) = u^2 / c.a2 + v^2 / c.b2 <= c.qcut

@inline function _sg_eval(c::_SGComponent, u, v, s, cscale = 1.0)
    a2 = c.a2 + cscale * c.defocus_c * s^2
    b2 = c.b2 + cscale * c.defocus_c * s^2
    return c.power_norm / sqrt(a2 * b2) * exp(-(u^2 / a2 + v^2 / b2)^c.halfn)
end

struct BeamFrame{T}
    k::Vec3
    e1::Vec3
    e2::Vec3
    focus::SVector{3,T}
    core::_SGComponent{T}
    halo::_SGComponent{T}
    has_halo::Bool
    fnumber::Float64
end

function BeamFrame(b::Beam)
    k = direction(b)
    e1, e2 = spot_axes(b)
    sp = b.spot
    focus = b.aim + b.defocus * k
    cone2 = (1 / (2b.fnumber))^2
    n, nh, h, sc = sp.order, sp.halo_order, sp.halo, sp.halo_scale
    core = _SGComponent(sp.a^2, sp.b^2, n / 2, cone2 / (4 * _sg_moment(n)),
                        (1 - h) * b.power / _sg_norm(n))
    halo = _SGComponent((sc * sp.a)^2, (sc * sp.b)^2, nh / 2, cone2 / (4 * _sg_moment(nh)),
                        h * b.power / _sg_norm(nh))
    return BeamFrame(k, e1, e2, focus, core, halo, h > 0, b.fnumber)
end

"""
    intensity(frame_or_beam, x)

Beam intensity (TW/mm², per unit area normal to the beam) at point `x`. Away from best
focus the spot is the focal profile convolved with the geometric f-cone disk; this is
approximated by widening the semi-axes so second moments add exactly
(`a(s)² = a² + s²/(16 F² m_n)` with `m_n = Γ(4/n)/(2Γ(2/n))`) at fixed order `n`, for
the core and halo components separately. This smooth profile serves envelope and
clearance estimates; far from focus its shape departs from the convolution, so
deposition (`DepositionOperator`, `raytrace_deposition`) sums over the cone directions
instead.
"""
@inline function intensity(f::BeamFrame, x)
    d = x - f.focus
    s = dot(d, f.k)
    u = dot(d, f.e1)
    v = dot(d, f.e2)
    I = _sg_eval(f.core, u, v, s)
    f.has_halo && (I += _sg_eval(f.halo, u, v, s))
    return I
end

"""
    focal_profile(spot, power, ρ)

Best-focus intensity (TW/mm²) at elliptical radius `ρ` (in units of the core 1/e
semi-axes) for a beam of `power` TW.
"""
function focal_profile(sp::SuperGaussianSpot, P, ρ)
    I = (1 - sp.halo) * P / (_sg_norm(sp.order) * sp.a * sp.b) * exp(-ρ^sp.order)
    if sp.halo > 0
        sc = sp.halo_scale
        I += sp.halo * P / (_sg_norm(sp.halo_order) * sc^2 * sp.a * sp.b) *
             exp(-(ρ / sc)^sp.halo_order)
    end
    return I
end

"""
    enclosed_fraction(spot, ρ)

Fraction of the beam power inside elliptical radius `ρ` (core 1/e units) at best focus.
"""
function enclosed_fraction(sp::SuperGaussianSpot, ρ)
    f = (1 - sp.halo) * gamma_inc(2 / sp.order, ρ^sp.order)[1]
    sp.halo > 0 && (f += sp.halo * gamma_inc(2 / sp.halo_order, (ρ / sp.halo_scale)^sp.halo_order)[1])
    return f
end

intensity(b::Beam, x) = intensity(BeamFrame(b), x)

# ---------------------------------------------------------------------------------------
# VisRad's beam profile
# ---------------------------------------------------------------------------------------

"""
    VisRadProfile(; envelope = 1.53645, focal_length = 7695.0)

VisRad's analytic beam model (VR/power_sources/set_laser_profile_props.html), for
deposition that reproduces VisRad's (`DepositionOperator(...; profile)`). The beam is
not the focal spot convolved over the f-cone, as in Briareus's default model. Instead
the super-Gaussian keeps its shape along the beam with semi-axes

    a(z) = a₀ √((z/z_A)² + 1),   z_A = FL [(FL / (2a₀ f#))² − 1]^(-1/2),

and likewise b(z) with z_B, for distance `z` from best focus and final-optic distance
`FL` (`focal_length`, mm; VisRad's port radius, 769.5 cm at NIF). Then a(FL) is the
f-cone's radius at the lens. The profile is cut off at the elliptical radius
`envelope` (in 1/e semi-axes), and its peak is renormalized so the envelope holds the
beam's power in every transverse plane. A surface absorbs the profile times the
incidence factor |n·k| of the beam axis, where the axis reaches it unobstructed.

What is not modeled, for lack of VisRad's constants (see `reference/research-notes.md`,
"The faint wings"):
- the exponent's change along the beam, n(z) = n₀[1 + 0.2 z/z₀], with z₀ unstated;
- the blend toward a square beam away from focus, whose weight
  f_sq = [z/√(z² + ab)]^Γ is not small at sub-mm defocus unless Γ is large;
- occultation by VisRad's ray grid. Here the beam axis through each point decides it.

Against a VisRad foil calculation the elements above 5% of peak agree within 1.4%.
The wings at sub-mm defocus differ by up to ~20%. VisRad's effective envelope there lay
near 1.64, against the file's ratio of 1.53645; the ratio's exact definition is
unresolved. `visrad_profile` takes the parameters from a workspace.
"""
struct VisRadProfile
    envelope::Float64
    focal_length::Float64
    function VisRadProfile(envelope, focal_length)
        envelope > 0 || throw(ArgumentError("envelope must be positive"))
        focal_length > 0 || throw(ArgumentError("focal_length must be positive"))
        return new(Float64(envelope), Float64(focal_length))
    end
end

VisRadProfile(; envelope = 1.53645, focal_length = 7695.0) = VisRadProfile(envelope, focal_length)

# One component of a beam under VisRad's profile: focal semi-axes squared, their
# Rayleigh-like ranges z_A², z_B², and the power normalization within the envelope.
struct _VRComponent{T}
    a2::T
    b2::T
    za2::T
    zb2::T
    halfn::Float64
    power_norm::T             # component power / (2π Γ(2/n)/n × fraction inside the envelope)
    q2::Float64               # envelope², in (u/a)² + (v/b)²
end

# z_R² = FL² / ((FL / (2 r f#))² − 1) for a 1/e semi-axis r
function _vr_range2(r2, F, FL)
    x = FL^2 / (4 * r2 * F^2)
    x > 1 || throw(ArgumentError("spot semi-axis exceeds the lens radius FL/(2f#)"))
    return FL^2 / (x - 1)
end

function _VRComponent(c::_SGComponent, F, prof::VisRadProfile)
    n = 2 * c.halfn
    inside = gamma_inc(2 / n, prof.envelope^n)[1]
    return _VRComponent(c.a2, c.b2, _vr_range2(c.a2, F, prof.focal_length),
                        _vr_range2(c.b2, F, prof.focal_length), c.halfn, c.power_norm / inside,
                        prof.envelope^2)
end

# Intensity (per unit area normal to the beam) at focal-frame coordinates (u, v), a
# distance s from best focus; 0 outside the envelope.
@inline function _vr_eval(c::_VRComponent, u, v, s)
    a2 = c.a2 * (1 + s^2 / c.za2)
    b2 = c.b2 * (1 + s^2 / c.zb2)
    q = u^2 / a2 + v^2 / b2
    q > c.q2 && return zero(c.power_norm * q)
    q > 0 || return c.power_norm / sqrt(a2 * b2) * one(q)
    return c.power_norm / sqrt(a2 * b2) * exp(-_pow(q, c.halfn))
end

@inline _vr_outside(c::_VRComponent, u, v, s) =
    u^2 / (c.a2 * (1 + s^2 / c.za2)) + v^2 / (c.b2 * (1 + s^2 / c.zb2)) > c.q2 * (1 + 1e-9)

struct _VRFrame{T}
    k::Vec3
    e1::Vec3
    e2::Vec3
    focus::SVector{3,T}
    core::_VRComponent{T}
    halo::_VRComponent{T}
    has_halo::Bool
end

_VRFrame(f::BeamFrame, prof::VisRadProfile) =
    _VRFrame(f.k, f.e1, f.e2, f.focus, _VRComponent(f.core, f.fnumber, prof),
             _VRComponent(f.halo, f.fnumber, prof), f.has_halo)

@inline function _vr_coords(f::_VRFrame, x)
    d = x - f.focus
    return dot(d, f.e1), dot(d, f.e2), dot(d, f.k)
end

@inline function _vr_intensity(f::_VRFrame, x)
    u, v, s = _vr_coords(f, x)
    I = _vr_eval(f.core, u, v, s)
    f.has_halo && (I += _vr_eval(f.halo, u, v, s))
    return I
end

@inline function _vr_beyond(f::_VRFrame, x)
    u, v, s = _vr_coords(f, x)
    return _vr_outside(f.core, u, v, s) && (!f.has_halo || _vr_outside(f.halo, u, v, s))
end

"""
    intensity(profile::VisRadProfile, beam, x)

Intensity (TW/mm², per unit area normal to the beam) of `beam` at point `x` under
VisRad's beam model (`VisRadProfile`).
"""
intensity(prof::VisRadProfile, b::Beam, x) = _vr_intensity(_VRFrame(BeamFrame(b), prof), x)

"""
    cone_directions(beam, n) -> (dirs, weights)

`n`-point quadrature of the beam's converging f-cone directions (center plus rings,
equal-area weights). Directions are unit propagation vectors.
"""
function cone_directions(b::Beam, n::Integer)
    k = direction(b)
    n <= 1 && return [k], [1.0]
    e1, e2 = port_basis(b.port)
    θmax = atan(1 / (2b.fnumber))
    pts = _disk_quadrature(n)
    dirs = [normalize(k + tan(θmax) * (r * cos(ϕ) * e1 + r * sin(ϕ) * e2)) for (r, ϕ) in pts]
    return dirs, fill(1.0 / length(dirs), length(dirs))
end

# Equal-weight points on the unit disk: the center plus rings of 6, 12, 18, … points,
# each ring at the rms radius of its equal-area annulus. The radii are then scaled so
# the set reproduces the disk's second moment ⟨r²⟩ = 1/2 exactly, so that the
# defocused spot of the deposition operator and the ray trace, a sum over these
# directions, grows at the exact geometric-optics rate.
function _disk_quadrature(n)
    n <= 1 && return [(0.0, 0.0)]
    counts = Int[]
    while sum(counts; init = 0) < n - 1
        push!(counts, min(6 * (length(counts) + 1), n - 1 - sum(counts; init = 0)))
    end
    pts = [(0.0, 0.0)]
    lo = sqrt(1 / n)
    for (j, c) in enumerate(counts)
        hi = sqrt((1 + sum(counts[1:j])) / n)
        r = sqrt((hi^2 + lo^2) / 2)
        for q in 0:c-1
            push!(pts, (r, 2π * (q + 0.5 * (j - 1)) / c))
        end
        lo = hi
    end
    scale = sqrt(0.5 / (sum(p[1]^2 for p in pts) / n))
    return [(scale * r, ϕ) for (r, ϕ) in pts]
end
