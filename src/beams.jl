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

The axis keeps the port's direction and passes through `aim`, so direct-drive pointings
are aims: a point `s` on a target sphere about `c` is `aim = s` (the axis first meets
the sphere there when `s` faces the port, `port · (s - c) > 0`; for a point on the
sphere at polar angle θ and azimuth φ, `s = c + R * port_vector(θ, φ)`), a shift `Δ`
normal to the beam is `aim = c + Δ * u` with `u ⟂ port` (toward the equator for a
northern port: `u = θ̂` of `port_basis`), and an axis crossing the polar axis at height
`z` is `aim = (0, 0, z)`. Best focus is at `aim + defocus * k` in each case.
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
    VisRadProfile(; envelope = nothing, focal_length = 7695.0, growth = 0.67, measured = true)

VisRad's analytic beam model (VR/power_sources/set_laser_profile_props.html), for
deposition that reproduces VisRad 21's (`DepositionOperator(...; profile)`). The beam is
not the focal spot convolved over the f-cone, as in Briareus's default model, but a
super-Gaussian whose semi-axes a(z), b(z) and exponent n(z) change with the distance `z`
from best focus. It is cut off at the elliptical radius ρ (the envelope ratio, in 1/e
semi-axes a(z), b(z)) and renormalized so the envelope holds the beam's power in every
transverse plane. A surface absorbs the profile times the incidence factor |n·k| of the
beam axis, where the axis reaches it unobstructed.

**Measured plates** (`measured = true`, the default). For the four NIF phase plates
VisRad's profile along the beam is measured (`_VR_MEASURED`: one beam per quad deposited
on a plate normal to it at |z| = 0–15 mm, VisRad 21; `reference/research-notes.md`,
"VisRad deposition parity"): a(z)/a₀, b(z)/b₀ and n(z)/n₀ are tabulated in |z|
(symmetric about focus) and interpolated linearly, the last segment extended beyond
the table. n grows by about 1.1%/mm (inner plates) and 1.6%/mm (outer), the manual's
n₀(1 + 0.2|z|/z₀) with z₀ ≈ 18.5 and 12.7 mm; the semi-axes follow no hyperbola, the
minor growing faster than the major. A beam takes a plate's tables when its focal
exponent matches to 0.01 and its 1/e semi-axes to 1% (as `visrad_profile` reads them
from a workspace). With VisRad's quad-relative spot rotation (`read_visrad_workspace`)
this reproduces VisRad's hohlraum wall deposition to 0.9% rms per element.

**Other spots** (and `measured = false`) take the manual's hyperbolic semi-axes at fixed
n,

    a(z) = a₀ √((z/z_A)² + 1),   z_A = (γ/ρ) FL [(FL / (2a₀ f#))² − 1]^(-1/2),

and likewise b(z) with z_B, for final-optic distance `FL` (`focal_length`, mm; VisRad's
port radius, 769.5 cm at NIF) and `growth` γ. The manual gives z_A without the factor
γ/ρ (`growth = nothing` uses it so). γ = 0.67 was fitted to VisRad's hohlraum
deposition before the plates were measured (2.8% element rms there).

`envelope` is the ratio ρ for every beam (a number), per beam name (a `Dict`; beams
not in it take the default), or `nothing` for VisRad's default (ln 1000)^(1/n₀), the
radius where the focal profile falls to 10⁻³ of its peak (it reproduces the NIF
plates' ratios 1.53645, 1.49576, 1.27326 and 1.28955). `visrad_profile` takes the
envelopes and focal length from a workspace.

Not modeled: the blend toward a square beam (f_sq, zero for beams circular at the
final optics, as NIF's behave), the envelope's own dependence on n(z) (the cut is at
ρ in 1/e semi-axes; the profile there is ≲ 10⁻³ of peak), and occultation by VisRad's
ray grid (here the beam axis through each point decides it).
"""
struct VisRadProfile
    envelope::Union{Nothing,Float64,Dict{String,Float64}}
    focal_length::Float64
    growth::Union{Nothing,Float64}
    measured::Bool
    function VisRadProfile(envelope, focal_length, growth, measured = true)
        env = envelope isa AbstractDict ? Dict{String,Float64}(String(k) => v for (k, v) in envelope) :
              envelope === nothing ? nothing : Float64(envelope)
        all(>(0), env isa Dict ? values(env) : something(env, 1.0)) ||
            throw(ArgumentError("envelope must be positive"))
        focal_length > 0 || throw(ArgumentError("focal_length must be positive"))
        growth === nothing || growth > 0 || throw(ArgumentError("growth must be positive"))
        return new(env, Float64(focal_length), growth === nothing ? nothing : Float64(growth), Bool(measured))
    end
end

VisRadProfile(; envelope = nothing, focal_length = 7695.0, growth = 0.67, measured = true) =
    VisRadProfile(envelope, focal_length, growth, measured)

"""
    visrad_envelope(order) -> Float64

VisRad's default envelope ratio for a super-Gaussian of exponent `order`: the radius,
in 1/e units, where the focal profile falls to 10⁻³ of its peak, (ln 1000)^(1/n).
Inferred: it reproduces every NIF plate's ratio in VisRad 21 workspaces.
"""
visrad_envelope(order) = log(1000.0)^(1 / order)

# The envelope ratio `prof` gives beam `b`.
function _vr_envelope(prof::VisRadProfile, b::Beam)
    e = prof.envelope
    e isa Float64 && return e
    e isa Dict && haskey(e, b.name) && return e[b.name]
    return visrad_envelope(b.spot.order)
end

# Factor on z_A² (and z_B²): (γ/ρ)², or 1 for the manual's z_A.
_vr_zfactor(prof::VisRadProfile, env) = prof.growth === nothing ? 1.0 : (prof.growth / env)^2

# VisRad 21's profile along the beam for the four NIF phase plates, measured with one beam
# per quad focused onto a plate normal to it at |z| from best focus (48 plates; thor
# visrad/profile_probe, `reference/research-notes.md`, "VisRad deposition parity"). Each
# plate was fitted by a rotated super-Gaussian (0–2.5% rms): the semi-axes and exponent
# over their focal values, the same before and after focus, normalized at z = 0. Matched
# to a beam by its focal exponent and 1/e semi-axes (mm).
const _VR_MEASURED = [
    (plate = "Rev3-Inner-2w (23°)", order = 7.6, a = 0.934, b = 0.6926,
     z = Float64[0, 2, 3, 4, 6, 9, 12],
     a_ratio = [1.00000, 1.00311, 1.01015, 1.01974, 1.04499, 1.09513, 1.15740],
     b_ratio = [1.00000, 1.01861, 1.03125, 1.04633, 1.08392, 1.15838, 1.25306],
     n_ratio = [1.00000, 1.02632, 1.03718, 1.04715, 1.06562, 1.09150, 1.11713]),
    (plate = "1.0 Rev1a-Inner-2w (30°)", order = 8.0, a = 0.9003, b = 0.6386,
     z = Float64[0, 2, 3, 4, 6, 9, 12],
     a_ratio = [1.00000, 1.00216, 1.00940, 1.01941, 1.04603, 1.09888, 1.16440],
     b_ratio = [1.00000, 1.02001, 1.03389, 1.05058, 1.09247, 1.17617, 1.28271],
     n_ratio = [1.00000, 1.02631, 1.03740, 1.04773, 1.06667, 1.09346, 1.12053]),
    (plate = "Rev3-Outer-2w (44°)", order = 4.8, a = 0.6908, b = 0.3994,
     z = Float64[0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15],
     a_ratio = [1.00000, 1.01046, 1.02604, 1.04585, 1.06938, 1.09625, 1.12622, 1.15903, 1.19449, 1.27238, 1.35796, 1.49642],
     b_ratio = [1.00000, 1.01897, 1.04567, 1.08037, 1.12310, 1.17373, 1.23194, 1.29742, 1.36980, 1.53353, 1.71984, 2.03510],
     n_ratio = [1.00000, 1.01959, 1.03594, 1.05063, 1.06445, 1.07798, 1.09152, 1.10537, 1.11990, 1.15129, 1.18648, 1.24972]),
    (plate = "Rev1-Outer-2w (50°)", order = 4.5, a = 0.6726, b = 0.3905,
     z = Float64[0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15],
     a_ratio = [1.00000, 1.01234, 1.02980, 1.05158, 1.07723, 1.10640, 1.13884, 1.17431, 1.21256, 1.29637, 1.38807, 1.53547],
     b_ratio = [1.00000, 1.02042, 1.04904, 1.08621, 1.13193, 1.18603, 1.24822, 1.31815, 1.39545, 1.57030, 1.76936, 2.10609],
     n_ratio = [1.00000, 1.01938, 1.03587, 1.05078, 1.06499, 1.07901, 1.09319, 1.10784, 1.12328, 1.15688, 1.19543, 1.26526]),
]

# The measured plate of a focal spot (exponent `n`, 1/e semi-axes `a`, `b` in mm), if any:
# its index in `_VR_MEASURED` (0 for none; an index, so frame construction allocates nothing)
function _vr_measured_index(n, a, b)
    for (i, m) in enumerate(_VR_MEASURED)
        abs(n - m.order) <= 0.01 && abs(a / m.a - 1) <= 0.01 && abs(b / m.b - 1) <= 0.01 && return i
    end
    return 0
end
_vr_measured(n, a, b) = (i = _vr_measured_index(n, a, b); i == 0 ? nothing : _VR_MEASURED[i])
_vr_measured(prof::VisRadProfile, sp::SuperGaussianSpot) =
    prof.measured ? _vr_measured(sp.order, ForwardDiff.value(max(sp.a, sp.b)), ForwardDiff.value(min(sp.a, sp.b))) :
    nothing

# One component of a beam under VisRad's profile: focal semi-axes squared, their
# Rayleigh-like ranges z_A², z_B², and the power normalization within the envelope; for a
# measured plate (`_VR_MEASURED`), its tables in |z| (empty otherwise): semi-axis and
# exponent ratios and the normalization relative to focus (the exponent's change moves
# the power inside the envelope).
struct _VRComponent{T}
    a2::T
    b2::T
    za2::T
    zb2::T
    halfn::Float64
    power_norm::T             # component power / (2π Γ(2/n)/n × fraction inside the envelope)
    q2::Float64               # envelope², in (u/a)² + (v/b)²
    zs::Vector{Float64}
    ra::Vector{Float64}       # a(z)/a₀ (a is the component's first axis, `spot_axes`)
    rb::Vector{Float64}
    rn::Vector{Float64}
    rnorm::Vector{Float64}
end

# z_R² = FL² / ((FL / (2 r f#))² − 1) for a 1/e semi-axis r
function _vr_range2(r2, F, FL)
    x = FL^2 / (4 * r2 * F^2)
    x > 1 || throw(ArgumentError("spot semi-axis exceeds the lens radius FL/(2f#)"))
    return FL^2 / (x - 1)
end

function _VRComponent(c::_SGComponent, F, prof::VisRadProfile, env)
    n = 2 * c.halfn
    inside = gamma_inc(2 / n, env^n)[1]
    g = _vr_zfactor(prof, env)
    a, b = sqrt(ForwardDiff.value(c.a2)), sqrt(ForwardDiff.value(c.b2))
    i = prof.measured ? _vr_measured_index(n, max(a, b), min(a, b)) : 0
    zs, ra, rb, rn, rnorm = if i == 0
        _VR_NONE, _VR_NONE, _VR_NONE, _VR_NONE, _VR_NONE
    else
        m = _VR_MEASURED[i]
        # a table's a/b follow the major/minor axis, whichever of the component's axes it is
        ra, rb = a >= b ? (m.a_ratio, m.b_ratio) : (m.b_ratio, m.a_ratio)
        m.z, ra, rb, m.n_ratio, _vr_rnorm(i, n, env)
    end
    return _VRComponent(c.a2, c.b2, g * _vr_range2(c.a2, F, prof.focal_length),
                        g * _vr_range2(c.b2, F, prof.focal_length), c.halfn, c.power_norm / inside,
                        env^2, zs, ra, rb, rn, rnorm)
end

const _VR_NONE = Float64[]           # no table (shared; never mutated)

# A measured plate's normalization relative to focus at each table node, for exponent `n`
# and envelope `env`: built once per (plate, n, env), since each node takes a `gamma_inc`.
# Frames are built per beam and work unit, so reads take no lock and allocate nothing: a
# miss copies the cache, adds the entry and swaps it in under the lock.
const _VRNormKey = Tuple{Int,Float64,Float64}       # (plate index, n, env): isbits
const _VR_RNORM = Ref(Dict{_VRNormKey,Vector{Float64}}())
const _VR_RNORM_LOCK = ReentrantLock()
function _vr_rnorm(i, n, env)
    key = (i, Float64(n), Float64(env))
    v = get(_VR_RNORM[], key, nothing)
    v === nothing || return v
    lock(_VR_RNORM_LOCK)
    try
        d = copy(_VR_RNORM[])
        v = get!(d, key) do
            norm0 = _sg_norm(n) * gamma_inc(2 / n, env^n)[1]
            [norm0 / (_sg_norm(n * r) * gamma_inc(2 / (n * r), env^(n * r))[1]) for r in _VR_MEASURED[i].n_ratio]
        end
        _VR_RNORM[] = d
        return v
    finally
        unlock(_VR_RNORM_LOCK)
    end
end

# Squared semi-axes, half exponent and normalization factor a distance s from best focus;
# tables are interpolated linearly in |s|, the last segment extended beyond them
@inline function _vr_axes(c::_VRComponent, s)
    if isempty(c.zs)
        return c.a2 * (1 + s^2 / c.za2), c.b2 * (1 + s^2 / c.zb2), c.halfn * one(s), one(s)
    end
    zs = c.zs
    z = abs(s)
    zv = ForwardDiff.value(z)
    k = 1
    @inbounds while k < length(zs) - 1 && zs[k+1] < zv
        k += 1
    end
    @inbounds w = (z - zs[k]) / (zs[k+1] - zs[k])
    ra, rb = _vr_lerp(c.ra, k, w), _vr_lerp(c.rb, k, w)
    return c.a2 * ra^2, c.b2 * rb^2, c.halfn * _vr_lerp(c.rn, k, w), _vr_lerp(c.rnorm, k, w)
end

# (no closure over the loop-updated `k`: it would be boxed)
@inline _vr_lerp(y, k, w) = @inbounds y[k] + (y[k+1] - y[k]) * w

# Intensity (per unit area normal to the beam) at focal-frame coordinates (u, v), a
# distance s from best focus; 0 outside the envelope.
@inline function _vr_eval(c::_VRComponent, u, v, s)
    a2, b2, hn, g = _vr_axes(c, s)
    q = u^2 / a2 + v^2 / b2
    q > c.q2 && return zero(c.power_norm * q * g)
    I0 = g * c.power_norm / sqrt(a2 * b2)
    q > 0 || return I0 * one(q)
    return I0 * exp(-_pow(q, hn))
end

@inline function _vr_outside(c::_VRComponent, u, v, s)
    a2, b2, _, _ = _vr_axes(c, s)
    return u^2 / a2 + v^2 / b2 > c.q2 * (1 + 1e-9)
end

struct _VRFrame{T}
    k::Vec3
    e1::Vec3
    e2::Vec3
    focus::SVector{3,T}
    core::_VRComponent{T}
    halo::_VRComponent{T}
    has_halo::Bool
end

# `env`: the beam's envelope ratio (`_vr_envelope`)
_VRFrame(f::BeamFrame, prof::VisRadProfile, env) =
    _VRFrame(f.k, f.e1, f.e2, f.focus, _VRComponent(f.core, f.fnumber, prof, env),
             _VRComponent(f.halo, f.fnumber, prof, env), f.has_halo)

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
intensity(prof::VisRadProfile, b::Beam, x) =
    _vr_intensity(_VRFrame(BeamFrame(b), prof, _vr_envelope(prof, b)), x)

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
