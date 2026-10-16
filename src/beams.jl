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
    SuperGaussianSpot(a, b = a; order = 2.0, angle = 0.0)

Elliptical super-Gaussian best-focus intensity `I ∝ exp(-[(u/a)² + (v/b)²]^(order/2))`
with 1/e semi-axes `a`, `b` (mm). `angle` (radians) rotates the `a` axis from the
port's θ̂ toward φ̂.
"""
struct SuperGaussianSpot{T<:Real}
    a::T
    b::T
    order::Float64
    angle::Float64
end

SuperGaussianSpot(a::Real, b::Real = a; order = 2.0, angle = 0.0) =
    SuperGaussianSpot(promote(a, b)..., Float64(order), Float64(angle))

# ∫ exp(-ρ^n) dA over the unit-semi-axis ellipse: 2π Γ(2/n) / n
_sg_norm(n) = 2π * gamma(2 / n) / n
# <u²>/a² for the super-Gaussian: Γ(4/n) / (2 Γ(2/n))
_sg_moment(n) = gamma(4 / n) / (2 * gamma(2 / n))

"""
    Beam(name; theta, phi, power = 1.0, spot, fnumber = 20.0, aim = (0,0,0),
         defocus = 0.0, group = :none, quad = "")

One laser beam entering from the port at polar angle `theta` and azimuth `phi`
(degrees, target-chamber frame, z up). `power` in TW. `group` labels beams that share
pointing parameters (e.g. a cone).
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
end

function Beam(name::AbstractString; theta, phi, power = 1.0, spot = SuperGaussianSpot(0.5),
              fnumber = 20.0, aim = (0.0, 0.0, 0.0), defocus = 0.0, group = :none, quad = "")
    p = port_vector(theta, phi)
    T = promote_type(eltype(aim), typeof(power), typeof(defocus), typeof(spot.a))
    sp = SuperGaussianSpot{T}(spot.a, spot.b, spot.order, spot.angle)
    return Beam{T}(String(name), p, SVector{3,T}(aim...), T(power), T(defocus), sp,
                   Float64(fnumber), Symbol(group), String(quad))
end

"""
    port_vector(theta, phi)

Unit vector from TCC to a port at polar angle `theta`, azimuth `phi` (degrees).
"""
port_vector(θ, φ) = Vec3(sind(θ) * cosd(φ), sind(θ) * sind(φ), cosd(θ))

"Polar and azimuthal angles (degrees) of a beam's port."
port_angles(b::Beam) = (acosd(clamp(b.port[3], -1, 1)), mod(atand(b.port[2], b.port[1]), 360))

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

function spot_axes(b::Beam)
    eθ, eφ = port_basis(b.port)
    ψ = b.spot.angle
    return cos(ψ) * eθ + sin(ψ) * eφ, -sin(ψ) * eθ + cos(ψ) * eφ
end

"""
    remake(beam; aim, power, defocus, spot)

Copy of `beam` with the given fields replaced; the element type widens as needed (so
ForwardDiff duals can flow through pointing parameters).
"""
function remake(b::Beam; aim = b.aim, power = b.power, defocus = b.defocus, spot = b.spot)
    T = promote_type(eltype(aim), typeof(power), typeof(defocus), typeof(spot.a), typeof(spot.b))
    sp = SuperGaussianSpot{T}(spot.a, spot.b, spot.order, spot.angle)
    return Beam{T}(b.name, b.port, SVector{3,T}(aim...), T(power), T(defocus), sp, b.fnumber,
                   b.group, b.quad)
end

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
struct BeamFrame{T}
    k::Vec3
    e1::Vec3
    e2::Vec3
    focus::SVector{3,T}
    a2::T
    b2::T
    halfn::Float64
    defocus_c::Float64        # (1/(2F))² / (4 m_n): adds s²·c to a², b²
    power_norm::T             # power / (2π Γ(2/n)/n)
end

function BeamFrame(b::Beam)
    k = direction(b)
    e1, e2 = spot_axes(b)
    n = b.spot.order
    focus = b.aim + b.defocus * k
    c = (1 / (2b.fnumber))^2 / (4 * _sg_moment(n))
    return BeamFrame(k, e1, e2, focus, b.spot.a^2, b.spot.b^2, n / 2, c, b.power / _sg_norm(n))
end

"""
    intensity(frame_or_beam, x)

Beam intensity (TW/mm², per unit area normal to the beam) at point `x`. Away from best
focus the spot is the focal profile convolved with the geometric f-cone disk; this is
approximated by widening the semi-axes so second moments add exactly
(`a(s)² = a² + s²/(16 F² m_n)` with `m_n = Γ(4/n)/(2Γ(2/n))`) at fixed order `n`.
"""
@inline function intensity(f::BeamFrame, x)
    d = x - f.focus
    s = dot(d, f.k)
    u = dot(d, f.e1)
    v = dot(d, f.e2)
    a2 = f.a2 + f.defocus_c * s^2
    b2 = f.b2 + f.defocus_c * s^2
    ρ2 = u^2 / a2 + v^2 / b2
    return f.power_norm / sqrt(a2 * b2) * exp(-ρ2^f.halfn)
end

intensity(b::Beam, x) = intensity(BeamFrame(b), x)

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

# Equal-weight points on the unit disk: center + rings (each ring at its area centroid).
function _disk_quadrature(n)
    n <= 1 && return [(0.0, 0.0)]
    pts = [(0.0, 0.0)]
    rest = n - 1
    nrings = rest <= 8 ? 1 : 2
    counts = nrings == 1 ? [rest] : [rest ÷ 3, rest - rest ÷ 3]
    edges = [sqrt(1 / n)]
    for (j, c) in enumerate(counts)
        lo = edges[end]
        hi = sqrt((1 + sum(counts[1:j])) / n)
        r = (2 / 3) * (hi^3 - lo^3) / (hi^2 - lo^2)
        for q in 0:c-1
            push!(pts, (r, 2π * (q + 0.5 * (j - 1)) / c))
        end
        push!(edges, hi)
    end
    return pts
end
