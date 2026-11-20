# This file defines laser beams, with port geometry, focal-spot profiles, pointing and
# defocus.
#
# Conventions
# * A beam enters from its port, a unit vector `port` from target-chamber center (TCC)
#   toward the final optics; it propagates along `k = -port`. Pointing translates the
#   beam axis to pass through `aim`, and the direction stays fixed. The final optics are
#   meters away, so mm-scale repointing changes the direction by ≲1e-4 rad.
# * The best-focus plane lies at `aim + defocus * k`. Positive defocus places focus
#   beyond the aim point, so the beam is still converging there.
# * Spot axes: `a` lies along `cos(ψ) θ̂ + sin(ψ) φ̂` of the port, `b` perpendicular,
#   where θ̂ and φ̂ are the spherical unit vectors at the port.

"""
    SuperGaussianSpot(a, b = a; order = 2.0, angle = 0.0, halo = 0.0, halo_scale = 3.0,
                      halo_order = 2.0)

Elliptical super-Gaussian best-focus intensity `I ∝ exp(-[(u/a)² + (v/b)²]^(order/2))`
with 1/e semi-axes `a`, `b` (mm). `angle` (radians) rotates the `a` axis from θ̂ of the
port toward φ̂.

Measured phase-plate spots have low-level wings beyond the super-Gaussian core fit.
A fraction `halo` of the power can be put in a second, broader component with semi-axes
`halo_scale·(a, b)` and order `halo_order`. The halo matters for clearance and clipping
checks; drive symmetry is insensitive to it.
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

The axis has the direction of the port and passes through `aim`, so direct-drive
pointings are aims:

- a point `s` on a target sphere about `c` is `aim = s`. The axis first meets the sphere
  there when `s` faces the port, `port · (s - c) > 0`. For a point on the sphere at polar
  angle θ and azimuth φ, `s = c + R * port_vector(θ, φ)`;
- a shift `Δ` normal to the beam is `aim = c + Δ * u` with `u ⟂ port`, for example
  toward the equator for a northern port with `u = θ̂` of `port_basis`;
- an axis crossing the polar axis at height `z` is `aim = (0, 0, z)`.

Best focus is at `aim + defocus * k` in each case.
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

Spherical unit vectors at the port, which are the transverse frame for spot orientation.
The spot `a` axis of a beam at `angle` is `cos(angle) θ̂ + sin(angle) φ̂`
([`spot_axes`](@ref)). An aim displaced along θ̂ of the port by `d` mm is:

```julia
eθ, eφ = port_basis(b.port)
moved = remake(b; aim = Tuple(b.aim + d * eθ))
```
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

Unit vectors of the `a` and `b` axes of the spot, transverse to the beam. They are θ̂
and φ̂ of the port turned by the `angle` of the spot toward φ̂.
"""
function spot_axes(b::Beam)
    eθ, eφ = port_basis(b.port)
    ψ = b.spot.angle
    return cos(ψ) * eθ + sin(ψ) * eφ, -sin(ψ) * eθ + cos(ψ) * eφ
end

"""
    remake(beam; aim, power, defocus, spot, dlambda, port, fnumber, group, name)

Copy of `beam` with the given fields replaced. The element type widens as needed, so
that ForwardDiff duals can flow through pointing parameters.
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

# One super-Gaussian component of a spot in a beam frame.
struct _SGComponent{T}
    a2::T
    b2::T
    halfn::Float64
    defocus_c::Float64        # (1/(2F))² / (4 m_n): adds s²·c to a², b²
    power_norm::T             # component power / (2π Γ(2/n)/n)
    qcut::Float64             # (u/a)² + (v/b)² beyond which the profile is < e⁻⁴⁰ of peak
    peak::T                   # power_norm / (a b), the best-focus peak intensity
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

# Computes `q^h` for q > 0 as exp(h log q). This costs half as much as `^`, whose extra
# accuracy (about 1e-15 relative) the profile does not need, and a quarter as much on
# duals, where `^` takes a second power for the derivative.
@inline _pow(q, h) = exp(h * log(q))

# Whether (u, v) is inside the cutoff of the profile (`_sg_focal` is nonzero there).
@inline _sg_inside(c::_SGComponent, u, v) = u^2 / c.a2 + v^2 / c.b2 <= c.qcut

@inline function _sg_eval(c::_SGComponent, u, v, s, cscale = 1.0)
    a2 = c.a2 + cscale * c.defocus_c * s^2
    b2 = c.b2 + cscale * c.defocus_c * s^2
    return c.power_norm / sqrt(a2 * b2) * exp(-(u^2 / a2 + v^2 / b2)^c.halfn)
end

"""
    BeamFrame(beam)

Precomputed quantities for evaluating the intensity of a beam at many points.
"""
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
the core and halo components separately. This smooth profile is used for envelope and
clearance estimates. Far from focus its shape departs from the convolution, so
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
# The VisRad beam profile
# ---------------------------------------------------------------------------------------

"""
    VisRadProfile(; envelope = nothing, taper = 0)

The beam model of VisRad 21, for deposition that reproduces the VisRad deposition
(`DepositionOperator(...; profile)`). The beam is the analytic profile of the manual
(`VR/power_sources/set_laser_profile_props.html`), a super-Gaussian blended toward a
square away from best focus. The default model of Briareus instead convolves the focal
spot over the f-cone. At a distance `z` from best focus, with (u, v) along the spot
axes,

    I = I₀(z) exp(−α),   α = f [(|u|/ΔL)^M + (|v|/ΔL)^M] + (1 − f) [(u/a)² + (v/b)²]^(n/2).

For focal 1/e semi-axes a₀ ≥ b₀, exponent n₀, f-number F and envelope ratio ρ, with
``z_A = \\sqrt{2} a_0 F`` and ``z_B = \\sqrt{2} b_0 F``:

- a(z)² = a₀² + z²/(2F²), and likewise b(z);
- ``n(z) = n_0 (1 + 0.2 |z|/z_A)``;
- the weight of the square, ``f(z) = |z| / \\sqrt{z^2 + z_A z_B}``;
- its half-width ΔL(z) = ρ a₀/√2 + 0.99605 |z|/(2F);
- its exponent ``M(z) = 10 (1 + |z|/\\sqrt{z_A z_B})``.

The profile is cut where α reaches ln 1000 (10⁻³ of its local peak), and I₀(z) puts the
power of the beam in every transverse plane. A surface absorbs I |n·k| of the beam axis,
where the axis reaches it unobstructed.

The forms of a, b, n, f and ΔL are those of the VisRad parameters, which VisRad exports
with the beam–target intersection points. These were exported at 299 beam hits
(|z| ≤ 15 mm), and the forms reproduce them to 5 × 10⁻⁵ in a and b, 1.4 × 10⁻⁴ in n,
5 × 10⁻⁴ in f, and 1.3 × 10⁻³ in ΔL. Only M and the cut were fitted, to the VisRad
deposition on 48 plates normal to single NIF beams (0.06–0.22% rms per element). On
32 tilted plates and on normal plates up to 15 mm past focus, the model is then within
0.12–0.27% rms of VisRad, against 1.3–4.6% for the tabulated profile that it replaces
(`reference/research-notes.md`, "VisRad's beam profile law").

`envelope` sets ρ for every beam (a number), per beam name (a `Dict`, where beams not in
it take the default), or `nothing` for the VisRad default (ln 1000)^(1/n₀). The default
puts the focal profile at 10⁻³ of its peak, and it reproduces the ratios 1.53645,
1.49576, 1.27326 and 1.28955 of the NIF plates. VisRad 21 replaced larger ratios (2 and
3) with the default, so ρ is capped there. ρ enters only ΔL: with a smaller ratio, 1.2,
the VisRad deposition does not show a cut-off at q = ρ. `visrad_profile` takes the envelopes
from a workspace.

`taper` is an option for optimization and is not part of the VisRad model. The default,
0, is the hard cut of VisRad. At the cut the profile drops from 10⁻³ of its local peak
to zero. As an aim moves, the deposition, and every objective or constraint computed
from it, therefore jumps wherever an operator point crosses the cut, and finite
differences and line searches encounter steps that derivatives miss. A positive `taper`
multiplies exp(−α) by the quintic smoothstep S(t) = 10t³ − 15t⁴ + 6t⁵ of
t = (ln 1000 − α)/taper (S = 1 for t ≥ 1), and I₀(z) is renormalized so that each plane
still contains the power of the beam. The profile is then twice continuously
differentiable across the cut, and its support, which the sampling and screens of the
operator use, is unchanged. The band lies where the law is least constrained, since the
VisRad deposition fits a cut at α = ln 1000 or none equally well. `taper = 1` (one
e-fold, intensities between 10⁻³ and e × 10⁻³ of the local peak) changes the element rms
against VisRad on 59 probe plates by at most 0.011 percentage points, and NIF wall
deposition by at most 7 × 10⁻⁵ of its peak (`reference/research-notes.md`, "VisRad's
beam profile law").

Some dependences were not identified:
- any dependence on the final-optic distance (every run had 769.5 cm) or on the
  f-number beyond F itself (every run had f/20.755, where the ΔL slope was measured);
- the orientation of the square for a beam whose spot VisRad rotates, which is taken
  on the spot axes;
- occultation by the ray grid of VisRad, which here is set by the beam axis through each
  point.

VisRad also renormalizes the deposited total of each beam, which per-plane
normalization misses on steep surfaces (+0.2% at 40° incidence, +2.2% at 78°).
"""
struct VisRadProfile
    envelope::Union{Nothing,Float64,Dict{String,Float64}}
    taper::Float64
    function VisRadProfile(envelope, taper)
        env = envelope isa AbstractDict ? Dict{String,Float64}(String(k) => v for (k, v) in envelope) :
              envelope === nothing ? nothing : Float64(envelope)
        all(>(0), env isa Dict ? values(env) : something(env, 1.0)) ||
            throw(ArgumentError("envelope must be positive"))
        0 <= taper < _VR_CUT || throw(ArgumentError("taper must be in [0, ln 1000)"))
        return new(env, Float64(taper))
    end
end

VisRadProfile(; envelope = nothing, taper::Real = 0.0) = VisRadProfile(envelope, taper)

Base.:(==)(p::VisRadProfile, q::VisRadProfile) = p.envelope == q.envelope && p.taper == q.taper
Base.hash(p::VisRadProfile, h::UInt) = hash(p.envelope, hash(p.taper, hash(:VisRadProfile, h)))

"""
    visrad_envelope(order) -> Float64

The default envelope ratio of VisRad for a super-Gaussian of exponent `order`. It is the
radius, in 1/e units, where the focal profile falls to 10⁻³ of its peak,
(ln 1000)^(1/n). The value is inferred, and it reproduces the ratio of every NIF plate
in VisRad 21 workspaces.
"""
visrad_envelope(order) = log(1000.0)^(1 / order)

# Returns the envelope ratio that `prof` gives beam `b`, capped at the default as VisRad
# caps it.
function _vr_envelope(prof::VisRadProfile, b::Beam)
    e = prof.envelope
    ρd = visrad_envelope(b.spot.order)
    e isa Float64 && return min(e, ρd)
    e isa Dict && haskey(e, b.name) && return min(e[b.name], ρd)
    return ρd
end

# Constants of the law (`VisRadProfile`): the exponent of the square at focus, the slope
# of its half-width in units of 1/(2F) (measured at f/20.755), and the cut-off α
const _VR_M0 = 10.0
const _VR_DL_SLOPE = 0.99605
const _VR_CUT = log(1000.0)

# One component of a beam under the VisRad profile. a2, b2 are the squared focal 1/e
# semi-axes of the component, along its axes, and the first may be the minor one. The
# lengths of the law are taken from the major and minor axes. ntab is the normalization
# ∫∫ exp(−α) du dv (times the taper, when there is one) over a(z) b(z) at
# |z| = 0, h, 2h, …, held constant beyond the table.
struct _VRComponent{T}
    a2::T
    b2::T
    g2::Float64               # 1/(2F²): a(z)² = a₀² + z² g2
    halfn::Float64            # n₀/2
    inv_za::T                 # 1/z_A
    zab::T                    # z_A z_B
    dl0::T                    # ρ a₀/√2
    dls::Float64              # dΔL/d|z|
    power::T                  # power of the component (TW)
    ntab::Vector{Float64}
    h::Float64
    taper::Float64            # width in α of the taper below the cut (0 for the hard cut)
    qfar::Float64             # bound of `_vr_far` on the q of the Gaussian term
    gk::Int                   # n₀/2 ≥ gk + gj/4, the power of q in the bound of `_vr_far`
    gj::Int
end

function _VRComponent(c::_SGComponent, F, env, taper, active::Bool)
    n = 2 * c.halfn
    a, b = sqrt(c.a2), sqrt(c.b2)
    amaj = max(a, b)
    za, zb = sqrt(2) * amaj * F, sqrt(2) * min(a, b) * F
    dl0 = env * amaj / sqrt(2)
    dls = _VR_DL_SLOPE / (2F)
    ntab = active ? _vr_ntab(ForwardDiff.value(a), ForwardDiff.value(b), n, Float64(F), Float64(env),
                             Float64(taper)) : _VR_NONE
    return _VRComponent(c.a2, c.b2, 1 / (2F^2), c.halfn, 1 / za, za * zb, dl0, dls,
                        c.power_norm * _sg_norm(n), ntab, _VR_NTAB_H, Float64(taper),
                        max(1.0, _VR_FAR_CUT^(1 / c.halfn)), _vr_far_powers(c.halfn)...)
end

# Returns the integer and quarter parts (k, j) of a Gaussian exponent n₀/2 ≥ k + j/4, for
# the bound q^(n₀/2) ≥ q^k (q^(1/4))^j of `_vr_far` at q ≥ 1.
function _vr_far_powers(hn)
    k = floor(Int, hn)
    return k, floor(Int, 4 * (hn - k))
end

const _VR_NONE = Float64[]           # empty table for an inactive halo (shared, never mutated)

# Returns the α of the law at (u, v) a distance s from best focus, and a(z)², b(z)², |z|.
@inline function _vr_alpha(c::_VRComponent, u, v, s)
    z = abs(s)
    a2, b2 = c.a2 + z^2 * c.g2, c.b2 + z^2 * c.g2
    q = u^2 / a2 + v^2 / b2
    hn = c.halfn * (1 + 0.2 * z * c.inv_za)
    rz = sqrt(c.zab)
    f = z / sqrt(z^2 + c.zab)
    dl = c.dl0 + c.dls * z
    M = _VR_M0 * (1 + z / rz)
    sq = (u == 0 ? zero(dl) : _pow(abs(u) / dl, M)) + (v == 0 ? zero(dl) : _pow(abs(v) / dl, M))
    gauss = q > 0 ? _pow(q, hn) : zero(q)
    return f * sq + (1 - f) * gauss, a2, b2, z
end

# Returns the normalization over a(z) b(z) at |z|, interpolated linearly in the table.
@inline function _vr_norm(c::_VRComponent, z)
    t = c.ntab
    x = z / c.h
    xv = ForwardDiff.value(x)
    xv >= length(t) - 1 && return t[end] * one(x)
    k = floor(Int, xv) + 1
    w = x - (k - 1)
    return @inbounds t[k] + (t[k+1] - t[k]) * w
end

# Returns the intensity (per unit area normal to the beam) at focal-frame coordinates
# (u, v), a distance s from best focus, and 0 beyond the cut-off. With a taper, exp(−α)
# is multiplied by the smoothstep of (ln 1000 − α)/taper, which reaches zero at the cut
# with its first two derivatives, so the support is unchanged.
@inline function _vr_eval(c::_VRComponent, u, v, s)
    α, a2, b2, z = _vr_alpha(c, u, v, s)
    α >= _VR_CUT && return zero(c.power * α)
    I = c.power / (sqrt(a2 * b2) * _vr_norm(c, z)) * exp(-α)
    return c.taper > 0 ? I * _vr_smoothstep((_VR_CUT - α) / c.taper) : I
end

# The quintic smoothstep: 0 at t = 0 and 1 from t = 1 on, with continuous first and
# second derivatives (C²).
@inline _vr_smoothstep(t) = t >= 1 ? one(t) : t^3 * (10 - t * (15 - 6t))

@inline _vr_outside(c::_VRComponent, u, v, s) = _vr_alpha(c, u, v, s)[1] > _VR_CUT * (1 + 1e-9)

# Whether (u, v, s) lies beyond the cut-off, by a bound that needs no logarithm or
# exponential. With α = f sq + (1 − f) gauss and 0 ≤ f < 1:
# - sq ≥ r^M with r = max(|u|, |v|)/ΔL and M ≥ M₀ = 10, so sq ≥ r^10 for r ≥ 1;
# - gauss = q^hn with hn ≥ n₀/2 ≥ k + j/4, so gauss ≥ q^k (q^(1/4))^j for q ≥ 1.
# The test first requires both terms to be at least c, with r ≥ c^(1/10) and
# q ≥ max(1, c^(2/n₀)). It then tests their weighted sum with f = z/√(z² + z_A z_B)
# against c = ln 1000 (1 + 1e-6), far above the allowance of 1e-9 of the screens and the
# rounding of α. A point that passes lies beyond the cut of `_vr_outside`, and
# `_vr_eval` would give an exact zero there.
const _VR_FAR_CUT = _VR_CUT * (1 + 1e-6)
const _VR_FAR_T = _VR_FAR_CUT^(1 / _VR_M0)
@inline _vr_far(c::_VRComponent, u, v, s) = (z = abs(s); _vr_far(c, abs(u), abs(v), z, z))

# The same bound over a set of points with |u| ≥ au, |v| ≥ av and zlo ≤ |s| ≤ zhi. ΔL, a(z)
# and b(z) grow with |z|, so r and q are least at zhi. f also grows with |z|, and it is
# taken at whichever end gives the smaller sum.
@inline function _vr_far(c::_VRComponent, au, av, zlo, zhi)
    r = max(au, av) / (c.dl0 + c.dls * zhi)
    q = au^2 / (c.a2 + zhi^2 * c.g2) + av^2 / (c.b2 + zhi^2 * c.g2)
    r >= _VR_FAR_T && q >= c.qfar && return true
    sq = r >= 1 ? (r2 = r * r; r4 = r2 * r2; r4 * r4 * r2) : zero(r)
    g = zero(q)
    if q >= 1
        g = one(q)
        for _ in 1:c.gk
            g *= q
        end
        if c.gj > 0
            q4 = sqrt(sqrt(q))
            for _ in 1:c.gj
                g *= q4
            end
        end
    end
    flo = zlo / sqrt(zlo^2 + c.zab)
    fhi = zhi == zlo ? flo : zhi / sqrt(zhi^2 + c.zab)
    return min(flo * sq + (1 - flo) * g, fhi * sq + (1 - fhi) * g) >= _VR_FAR_CUT
end

# The normalization table of a component with focal semi-axes a, b (along its axes),
# exponent n, f-number F, envelope env and taper, on |z| = 0:h:zmax. Each node integrates
# exp(−α) over the quarter plane in scaled polar coordinates (u = a(z) r cos θ,
# v = b(z) r sin θ), out to the cut-off radius along each θ, since α grows with r. The
# integral uses composite Gauss–Legendre rules in θ and r, and the share of a taper is
# subtracted, integrated on the taper band. The table is built once per
# (a, b, n, F, env, taper) and cached. Frames are built per beam and work unit, so reads
# take no lock and allocate nothing. On a miss, the cache is copied, the entry is added
# and the copy is swapped in under the lock.
const _VR_NTAB_H = 0.1
const _VR_NTAB_ZMAX = 60.0
const _VRNormKey = NTuple{6,Float64}
const _VR_NTABS = Ref(Dict{_VRNormKey,Vector{Float64}}())
const _VR_NTABS_LOCK = ReentrantLock()

function _vr_ntab(a, b, n, F, env, taper)
    key = (a, b, n, F, env, taper)
    v = get(_VR_NTABS[], key, nothing)
    v === nothing || return v
    lock(_VR_NTABS_LOCK)
    try
        d = copy(_VR_NTABS[])
        v = get!(() -> _vr_ntab_build(a, b, n, F, env, taper), d, key)
        _VR_NTABS[] = d
        return v
    finally
        unlock(_VR_NTABS_LOCK)
    end
end

function _vr_ntab_build(a, b, n, F, env, taper)
    sp = _SGComponent(a^2, b^2, n / 2, 0.0, 1 / _sg_norm(n))
    c = _VRComponent(sp, F, env, taper, false)
    xg, wg = _gauss_legendre(16)                     # corona.jl
    zs = 0:_VR_NTAB_H:_VR_NTAB_ZMAX
    t = zeros(length(zs))
    Threads.@threads for i in eachindex(zs)
        t[i] = _vr_norm_at(c, zs[i], xg, wg)
    end
    return t
end

function _vr_norm_at(c::_VRComponent, z, xg, wg)
    az, bz = sqrt(c.a2 + z^2 * c.g2), sqrt(c.b2 + z^2 * c.g2)
    α(r, θ) = _vr_alpha(c, az * r * cos(θ), bz * r * sin(θ), z)[1]
    total = 0.0
    nθ, nr = 3, 3                                    # panels
    for pθ in 1:nθ, (xθ, wθ) in zip(xg, wg)
        θ = (π / 2) * (pθ - 1 + (xθ + 1) / 2) / nθ
        # the cut-off radius along θ, where α(r) = ln 1000, is bracketed and bisected
        hi = 1.0
        while α(hi, θ) < _VR_CUT
            hi *= 2
        end
        lo = 0.0
        for _ in 1:44
            mid = (lo + hi) / 2
            α(mid, θ) < _VR_CUT ? (lo = mid) : (hi = mid)
        end
        R = (lo + hi) / 2
        inner = 0.0
        for pr in 1:nr, (xr, wr) in zip(xg, wg)
            r = R * (pr - 1 + (xr + 1) / 2) / nr
            inner += wr * r * exp(-α(r, θ))
        end
        if c.taper > 0
            total += wθ * (inner * R / (2nr) - _vr_taper_loss(r -> α(r, θ), R, c.taper, xg, wg))
        else
            total += wθ * inner * R / (2nr)
        end
    end
    return 4 * total * (π / 2) / (2nθ)
end

# Returns ∫ r exp(−α) (1 − smoothstep) dr over the band of the taper along one direction,
# from the radius where α = ln 1000 − taper (bisected) to the cut-off radius R.
function _vr_taper_loss(α, R, taper, xg, wg)
    lo, hi = 0.0, R
    for _ in 1:44
        mid = (lo + hi) / 2
        α(mid) < _VR_CUT - taper ? (lo = mid) : (hi = mid)
    end
    R1 = (lo + hi) / 2
    loss = 0.0
    np = 2                                           # panels
    for p in 1:np, (x, w) in zip(xg, wg)
        r = R1 + (R - R1) * (p - 1 + (x + 1) / 2) / np
        a = α(r)
        loss += w * r * exp(-a) * (1 - _vr_smoothstep((_VR_CUT - a) / taper))
    end
    return loss * (R - R1) / (2np)
end

struct _VRFrame{T}
    k::Vec3
    e1::Vec3
    e2::Vec3
    focus::SVector{3,T}
    core::_VRComponent{T}
    halo::_VRComponent{T}
    has_halo::Bool
    c0::SVector{3,T}          # coordinates of the origin, whose partials every point shares
end

# `env` is the envelope ratio of the beam (`_vr_envelope`), and `taper` that of the profile.
function _VRFrame(f::BeamFrame, env, taper)
    d0 = zero(Vec3) - f.focus
    return _VRFrame(f.k, f.e1, f.e2, f.focus, _VRComponent(f.core, f.fnumber, env, taper, true),
                    _VRComponent(f.halo, f.fnumber, env, taper, f.has_halo), f.has_halo,
                    SVector(dot(d0, f.e1), dot(d0, f.e2), dot(d0, f.k)))
end

@inline function _vr_coords(f::_VRFrame, x)
    d = x - f.focus
    return dot(d, f.e1), dot(d, f.e2), dot(d, f.k)
end

# Returns the coordinates of a point in the frame `f` from those in its Float64 frame
# (`uvs`, of the values of the beam). On duals these are the values with the partials of
# the coordinates of the origin, which every point shares, since the axes of the frame
# carry none. This replaces three dual dot products per point.
@inline _vr_coords(f::_VRFrame{Float64}, x, uvs) = uvs
@inline function _vr_coords(f::_VRFrame{ForwardDiff.Dual{G,Float64,N}}, x, (u, v, s)) where {G,N}
    D = ForwardDiff.Dual{G,Float64,N}
    return D(u, ForwardDiff.partials(f.c0[1])), D(v, ForwardDiff.partials(f.c0[2])),
           D(s, ForwardDiff.partials(f.c0[3]))
end
@inline _vr_coords(f::_VRFrame, x, uvs) = _vr_coords(f, x)

@inline _vr_intensity(f::_VRFrame, x) = _vr_intensity(f, _vr_coords(f, x)...)
@inline function _vr_intensity(f::_VRFrame, u, v, s)
    I = _vr_eval(f.core, u, v, s)
    f.has_halo && (I += _vr_eval(f.halo, u, v, s))
    return I
end

# Whether a point lies beyond the cut-off by the law itself, which is the reference
# against which the bound of `_vr_far` is checked.
@inline _vr_beyond(f::_VRFrame, x) = _vr_beyond(f, _vr_coords(f, x)...)
@inline _vr_beyond(f::_VRFrame, u, v, s) =
    _vr_outside(f.core, u, v, s) && (!f.has_halo || _vr_outside(f.halo, u, v, s))

# `_vr_beyond` by the bound of `_vr_far`, true only where `_vr_beyond` is (a Float64 frame).
@inline _vr_far(f::_VRFrame{Float64}, x) = _vr_far(f, _vr_coords(f, x)...)
@inline _vr_far(f::_VRFrame{Float64}, u, v, s) =
    _vr_far(f.core, u, v, s) && (!f.has_halo || _vr_far(f.halo, u, v, s))

# `_vr_far` for every point of a ball (center and radius). The coordinates of the frame
# are orthonormal, so the coordinates of each point lie within the radius of those of the
# center.
@inline function _vr_far_ball(f::_VRFrame{Float64}, (x, y, z, R)::NTuple{4,Float64})
    u, v, s = _vr_coords(f, Vec3(x, y, z))
    au, av = max(abs(u) - R, 0.0), max(abs(v) - R, 0.0)
    zlo, zhi = max(abs(s) - R, 0.0), abs(s) + R
    return _vr_far(f.core, au, av, zlo, zhi) && (!f.has_halo || _vr_far(f.halo, au, av, zlo, zhi))
end

"""
    intensity(profile::VisRadProfile, beam, x)

Intensity (TW/mm², per unit area normal to the beam) of `beam` at point `x` under the
VisRad beam model (`VisRadProfile`).
"""
intensity(prof::VisRadProfile, b::Beam, x) =
    _vr_intensity(_VRFrame(BeamFrame(b), _vr_envelope(prof, b), prof.taper), x)

"""
    cone_directions(beam, n) -> (dirs, weights)

`n`-point quadrature of the converging f-cone directions of the beam (center plus
rings, equal-area weights). Directions are unit propagation vectors.
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

# Equal-weight points on the unit disk, the center plus rings of 6, 12, 18, … points,
# with each ring at the rms radius of its equal-area annulus. The radii are then scaled
# so that the set reproduces the second moment of the disk, ⟨r²⟩ = 1/2, exactly. The
# defocused spot of the deposition operator and the ray trace, a sum over these
# directions, then grows at the exact geometric-optics rate.
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
