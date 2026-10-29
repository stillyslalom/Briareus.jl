# Rigid frames: a rotation and an origin in the chamber, the frames of ports, foils and
# target parts, coordinates in a frame, and a named graph of frames placed in each other.

"""
    Frame(R, origin)
    Frame(origin; z, x = θ̂ of z)

A right-handed rigid frame in the chamber: `R` holds the frame's unit axes as columns
and `origin` is its origin (mm), so a point with local coordinates `x` lies at
`R * x + origin` in the chamber ([`to_chamber`](@ref)), and a chamber point `p` has
local coordinates `R' * (p - origin)` ([`to_frame`](@ref)). [`CHAMBER`](@ref) is the
identity frame at target-chamber center (TCC).

The keyword form builds the frame from its `z` axis (any length) and a vector `x` toward
which its `x` axis points: `x` minus its component along `z`, normalized; by default the
polar unit vector θ̂ of `z` (`port_basis`; along chamber `x` when `z` is the pole).
`y = z × x` completes the right-handed triad. An `x` parallel to `z` is an
`ArgumentError`.

A frame is Float64 and constant; the points mapped through it may carry ForwardDiff
duals (the aims of [`FramePointing`](@ref)).
"""
struct Frame
    R::SMatrix{3,3,Float64,9}
    origin::Vec3
end

Frame(R::AbstractMatrix, origin) = Frame(SMatrix{3,3,Float64}(R), Vec3(origin...))

function Frame(origin; z, x = nothing)
    ez = normalize(Vec3(z...))
    all(isfinite, ez) || throw(ArgumentError("the frame's z axis must be a nonzero vector"))
    xr = x === nothing ? port_basis(ez)[1] : Vec3(x...)
    xp = xr - dot(xr, ez) * ez
    norm(xp) > 1e-9 * max(norm(xr), 1e-300) ||
        throw(ArgumentError("the frame's x direction is parallel to its z axis"))
    ex = normalize(xp)
    ey = cross(ez, ex)
    return Frame(hcat(ex, ey, ez), Vec3(origin...))
end

"""
    CHAMBER

The chamber frame: identity axes at target-chamber center.
"""
const CHAMBER = Frame(SMatrix{3,3,Float64}(I), zero(Vec3))

Base.:(==)(a::Frame, b::Frame) = a.R == b.R && a.origin == b.origin
Base.hash(f::Frame, h::UInt) = hash(f.origin, hash(f.R, hash(:Frame, h)))
Base.isapprox(a::Frame, b::Frame; atol::Real = 1e-12, rtol::Real = 0) =
    isapprox(a.R, b.R; atol, rtol) && isapprox(a.origin, b.origin; atol, rtol)

function Base.show(io::IO, f::Frame)
    v(x) = join((@sprintf("%.4g", c) for c in x), ", ")
    print(io, "Frame(origin (", v(f.origin), "), x (", v(f.R[:, 1]), "), y (", v(f.R[:, 2]),
          "), z (", v(f.R[:, 3]), "))")
end

"""
    to_chamber(frame, x; direction = false) -> SVector{3}

Chamber coordinates of the local point `x` of `frame`, `R * x + origin`; with
`direction = true`, of the local vector `x`, `R * x` (no translation). The element type
follows `x`, so duals pass through.
"""
to_chamber(f::Frame, x::AbstractVector; direction::Bool = false) =
    direction ? f.R * SVector{3}(x[1], x[2], x[3]) : f.R * SVector{3}(x[1], x[2], x[3]) + f.origin
to_chamber(f::Frame, x::Tuple; direction::Bool = false) = to_chamber(f, SVector(x...); direction)

"""
    to_frame(frame, p; direction = false) -> SVector{3}

Local coordinates in `frame` of the chamber point `p`, `R' * (p - origin)`; with
`direction = true`, of the chamber vector `p`, `R' * p`. The inverse of
[`to_chamber`](@ref).
"""
to_frame(f::Frame, p::AbstractVector; direction::Bool = false) =
    direction ? f.R' * SVector{3}(p[1], p[2], p[3]) :
    f.R' * (SVector{3}(p[1], p[2], p[3]) - f.origin)
to_frame(f::Frame, p::Tuple; direction::Bool = false) = to_frame(f, SVector(p...); direction)

"""
    compose(parent, child) -> Frame

The chamber frame of `child`, a frame given in `parent`'s local coordinates:
`to_chamber(compose(parent, child), x) == to_chamber(parent, to_chamber(child, x))` (to
round-off).
"""
compose(a::Frame, b::Frame) = Frame(a.R * b.R, a.R * b.origin + a.origin)

"""
    inv(frame) -> Frame

The inverse rigid map: `to_chamber(inv(f), p) == to_frame(f, p)` (to round-off), and
`compose(f, inv(f))` is the identity.
"""
Base.inv(f::Frame) = Frame(f.R', -(f.R' * f.origin))

# Rotations about z and y by `a` degrees (right-handed).
_Rz(a) = SMatrix{3,3,Float64}(cosd(a), sind(a), 0, -sind(a), cosd(a), 0, 0, 0, 1)
_Ry(a) = SMatrix{3,3,Float64}(cosd(a), 0, -sind(a), 0, 1, 0, sind(a), 0, cosd(a))

"""
    euler_rotation(polar, azimuth, rotation = 0) -> SMatrix{3,3}

The rotation `Rz(azimuth) · Ry(polar) · Rz(rotation)` (degrees, right-handed): turn by
`rotation` about `z`, tilt by `polar` about `y`, then turn by `azimuth` about `z`. Its
`z` column is `port_vector(polar, azimuth)`, and at `rotation = 0` its `x` and `y`
columns are the polar and azimuthal unit vectors θ̂ and φ̂ there (`port_basis`), so a
frame with this rotation has its `z` axis toward the direction (polar, azimuth) and
`rotation` turns its `x` axis from θ̂ toward φ̂. This is the order of VisRad's
orientation angles (Polar, Azimuthal, Rotation; the GUI's part placement uses the same).
"""
euler_rotation(polar, azimuth, rotation = 0.0) = _Rz(azimuth) * _Ry(polar) * _Rz(rotation)

"""
    euler_frame(origin, polar, azimuth, rotation = 0) -> Frame

The frame at `origin` with axes [`euler_rotation`](@ref)`(polar, azimuth, rotation)`.
"""
euler_frame(origin, polar, azimuth, rotation = 0.0) =
    Frame(euler_rotation(polar, azimuth, rotation), Vec3(origin...))

"""
    port_frame(theta, phi; roll = 0, origin = (0, 0, 0)) -> Frame

The frame of the port at polar angle `theta` and azimuth `phi` (degrees): `z` toward the
port (`port_vector(theta, phi)`, so a beam from it travels along `-z`), `x` the port's
θ̂ turned by `roll` (degrees) toward φ̂, at `origin` (TCC by default). It equals
[`euler_frame`](@ref)`(origin, theta, phi, roll)`. The `x` axis along θ̂ is an
assumption at low confidence (research notes, "Frames and placement (2026-10-28)"):
VisRad's port frames have not been seen in a sample workspace.
"""
port_frame(θ, φ; roll = 0.0, origin = (0.0, 0.0, 0.0)) = euler_frame(origin, θ, φ, roll)

"""
    foil_frame(los; standoff, tilt = 0.0, tilt_toward = 0.0, rotation = 0.0,
               center = (0, 0, 0)) -> Frame

The frame of a backlighter foil on the line of sight `los` = (θ, φ) in degrees, as
[`backlighter_foil`](@ref) places it: origin at the foil's center,
`center − standoff · port_vector(θ, φ)`; `x` along the foil's `width` edge, `y` along
its `height` edge and `z` along its front normal (toward the detector). Before the
tilt, `x` is the line of sight's θ̂ turned by `rotation` toward φ̂; `tilt` turns the
whole frame off the line of sight toward `cos(ψ) θ̂ + sin(ψ) φ̂`, `ψ = tilt_toward`
(degrees). Local `(u, v, 0)` is a point on the foil's face, `u` across its width.
"""
function foil_frame(los; standoff::Real, tilt::Real = 0.0, tilt_toward::Real = 0.0,
                    rotation::Real = 0.0, center = (0.0, 0.0, 0.0))
    d = port_vector(los...)
    eθ, eφ = port_basis(d)
    u = cosd(rotation) * eθ + sind(rotation) * eφ
    v = -sind(rotation) * eθ + cosd(rotation) * eφ
    e = cosd(tilt_toward) * eθ + sind(tilt_toward) * eφ
    β = deg2rad(tilt)
    u, v, nrm = _tilted(u, d, e, β), _tilted(v, d, e, β), _tilted(d, d, e, β)
    R = SMatrix{3,3,Float64}(u[1], u[2], u[3], v[1], v[2], v[3], nrm[1], nrm[2], nrm[3])
    return Frame(R, Vec3(center...) - standoff * d)
end

# Unit vector `los` rotated by `tilt` (rad) toward `e` (a unit vector ⟂ los), and the
# rotation that does it, applied to `v`.
function _tilted(v::Vec3, los::Vec3, e::Vec3, tilt)
    k = cross(los, e)                       # rotation axis (unit: los ⟂ e)
    return cos(tilt) * v + sin(tilt) * cross(k, v) + (1 - cos(tilt)) * dot(k, v) * k
end

"""
    part_centroid(target, part) -> Vec3
    part_centroid(target) -> Dict{Symbol,Vec3}

The area-weighted centroid of the triangles of `part` (all its pieces, see
[`part_pieces`](@ref)): the center of a foil, a disk or an entrance hole, a reference
point to aim beams at. The second form gives it for every part of the target.
"""
function part_centroid(t::Target, name::Symbol)
    cis = Set(part_pieces(t, name))
    c, a = zero(Vec3), 0.0
    for k in eachindex(t.tri_area)
        t.tri_comp[k] in cis || continue
        A = t.tri_area[k]
        c += A * (t.tri_v0[k] + (t.tri_e1[k] + t.tri_e2[k]) / 3)
        a += A
    end
    a > 0 || throw(ArgumentError("part $name has no area"))
    return c / a
end
part_centroid(t::Target) = Dict(p => part_centroid(t, p) for p in parts(t))

_surface_axis(s::CylinderSurface) = s.axis
_surface_axis(s::ConeSurface) = s.axis
_surface_axis(s::TorusSurface) = s.axis
_surface_axis(s::DiskSurface) = s.normal
_surface_axis(::Any) = nothing

"""
    component_frame(target, part) -> Frame

A frame for a part of a target built in a script (the GUI's parts have their own
placement frames): origin at its [`part_centroid`](@ref), and `z` axis, in order of
preference,
1. the axis of its exact surface: a `CylinderSurface`, `ConeSurface` or `TorusSurface`
   axis or a `DiskSurface` normal, when every piece with such a surface has a parallel
   one (the first piece's, with its sign);
2. else the area-weighted mean front normal of its triangles, when that mean is at
   least 0.1 in length (a flat foil: its front normal; a hemispherical cup: 0.5);
3. else the principal axis of its triangles' area moments whose moment differs most
   from the other two (the symmetry axis of a body of revolution), signed toward
   chamber `+z` (else `+x`, else `+y`); when the three moments agree within 1e-6 (a
   sphere, a cube) chamber `z`.
The `x` axis is θ̂ of `z` (see [`Frame`](@ref)). The choice is ambiguous for shapes
with no axis (a sphere or a closed box) and for the sign of a principal axis; the 0.1
threshold of rule 2 is an engineering default. Give the frame explicitly when it matters.
"""
function component_frame(t::Target, name::Symbol)
    cis = part_pieces(t, name)
    o = part_centroid(t, name)
    axes = [a for a in (_surface_axis(t.components[ci].surface) for ci in cis) if a !== nothing]
    if !isempty(axes) && all(a -> abs(dot(a, axes[1])) > 1 - 1e-12, axes)
        return Frame(o; z = axes[1])
    end
    cs = Set(cis)
    tris = [k for k in eachindex(t.tri_area) if t.tri_comp[k] in cs]
    area = sum(k -> t.tri_area[k], tris)
    nsum = sum(k -> t.tri_area[k] * t.tri_normal[k], tris)
    norm(nsum) >= 0.1 * area && return Frame(o; z = nsum)
    # area moments about the centroid, each triangle by its vertices' and centroid's
    # exact second moment: A/12 (Σ vᵢvᵢᵀ + 9 c cᵀ) about the origin o
    M = zeros(SMatrix{3,3,Float64,9})
    for k in tris
        v0 = t.tri_v0[k] - o
        v1, v2 = v0 + t.tri_e1[k], v0 + t.tri_e2[k]
        c = (v0 + v1 + v2) / 3
        M += t.tri_area[k] / 12 * (v0 * v0' + v1 * v1' + v2 * v2' + 9 * c * c')
    end
    λ, V = eigen(Symmetric(Matrix(M)))
    span = λ[3] - λ[1]
    span <= 1e-6 * max(abs(λ[3]), 1e-300) && return Frame(o; z = Vec3(0, 0, 1))
    # the eigenvalue farthest from the mean of the other two
    j = argmax(i -> abs(λ[i] - (sum(λ) - λ[i]) / 2), 1:3)
    z = Vec3(V[:, j]...)
    s = abs(z[3]) > 1e-12 ? z[3] : abs(z[1]) > 1e-12 ? z[1] : z[2]
    return Frame(o; z = s < 0 ? -z : z)
end

"""
    coords_of(p, system = :cartesian; frame = CHAMBER) -> NTuple{3,Float64}

Coordinates of the chamber point `p` in `frame`: `:cartesian` (x, y, z),
`:cylindrical` (r, φ, z: distance from the frame's `z` axis, azimuth in degrees in
[0, 360) from `x` toward `y`, height) or `:spherical` (r, θ, φ: distance from the
origin, polar angle from `z` and azimuth, degrees). At r = 0 the angles are 0.
[`point_of`](@ref) is the inverse.
"""
function coords_of(p, system::Symbol = :cartesian; frame::Frame = CHAMBER)
    x, y, z = to_frame(frame, Vec3(p...))
    system === :cartesian && return (x, y, z)
    φ = mod(atand(y, x), 360)
    system === :cylindrical && return (hypot(x, y), φ, z)
    system === :spherical || throw(ArgumentError("system must be :cartesian, :cylindrical or :spherical"))
    r = sqrt(x^2 + y^2 + z^2)
    return (r, r > 0 ? acosd(clamp(z / r, -1, 1)) : 0.0, φ)
end

"""
    point_of(c, system = :cartesian; frame = CHAMBER) -> Vec3

The chamber point with coordinates `c` in `frame`, in the order of
[`coords_of`](@ref): `:cartesian` (x, y, z), `:cylindrical` (r, φ, z) or `:spherical`
(r, θ, φ), angles in degrees. A negative radius is an `ArgumentError`.
"""
function point_of(c, system::Symbol = :cartesian; frame::Frame = CHAMBER)
    a, b, d = c
    system === :cartesian && return Vec3(to_chamber(frame, Vec3(a, b, d)))
    a >= 0 || throw(ArgumentError("the radius must not be negative"))
    system === :cylindrical && return Vec3(to_chamber(frame, Vec3(a * cosd(b), a * sind(b), d)))
    system === :spherical ||
        throw(ArgumentError("system must be :cartesian, :cylindrical or :spherical"))
    return Vec3(to_chamber(frame, Vec3(a * sind(b) * cosd(d), a * sind(b) * sind(d), a * cosd(b))))
end

abstract type PlacementStep end

"""
    FrameDef(parent = :tcc; offset = (0, 0, 0), rotation = (0, 0, 0), steps = [])

A frame placed in another: its origin at `offset` (mm) in the frame named `parent`, and
its axes turned by `rotation` = (polar, azimuth, rotation) in degrees in the parent
([`euler_rotation`](@ref)), so its chamber frame is
`compose(parent's frame, euler_frame(offset, rotation...))`; then the
[`PlacementStep`](@ref)s `steps` apply in order (relations to other frames, relative
moves). Degrees of freedom no step fixes keep the values of `offset` and `rotation`. A
port's frame is `FrameDef(:tcc; rotation = (θ, φ, roll))` ([`port_frame`](@ref)).
Resolve a set of named definitions with [`resolve_frames`](@ref).
"""
struct FrameDef
    parent::Symbol
    offset::Vec3
    rotation::NTuple{3,Float64}
    steps::Vector{PlacementStep}
end

FrameDef(parent::Symbol = :tcc; offset = (0.0, 0.0, 0.0), rotation = (0.0, 0.0, 0.0),
         steps = PlacementStep[]) =
    FrameDef(parent, Vec3(offset...), Tuple(Float64.(rotation))::NTuple{3,Float64},
             collect(PlacementStep, steps))

"""
    resolve_frames(defs, base = Dict(:tcc => CHAMBER)) -> Dict{Symbol,Frame}

Chamber frames of the named definitions `defs` (a dictionary or a collection of
`name => FrameDef` pairs), each placed in its parent, which is another definition or one
of the `base` frames (analytic frames: TCC, a template's entrance holes, a foil; `:tcc`
is added as [`CHAMBER`](@ref) when absent), then moved by its placement steps. The
result holds the base frames too. Definitions resolve depth first, the frames a
definition's steps name before it, so a frame follows every ancestor and every frame its
steps refer to exactly. A definition that reuses a base frame's name, a name given
twice, an unknown parent or step frame, frames that depend on each other in a loop (the
error names the loop, e.g. `a → b → a`), a step that cannot apply (a facing point at the
origin) and a step that undoes an earlier one (see [`placement_report`](@ref)) are each
an `ArgumentError`.
"""
resolve_frames(defs, base::AbstractDict = Dict(:tcc => CHAMBER)) = first(_resolve_frames(defs, base))

function _resolve_frames(defs, base::AbstractDict)
    out = Dict{Symbol,Frame}(Symbol(k) => v for (k, v) in base)
    haskey(out, :tcc) || (out[:tcc] = CHAMBER)
    pending = Dict{Symbol,FrameDef}()
    for (k, d) in defs
        name = Symbol(k)
        haskey(out, name) && throw(ArgumentError("frame $name is a base frame and cannot be redefined"))
        haskey(pending, name) && throw(ArgumentError("frame $name is defined twice"))
        pending[name] = d
    end
    reports = Dict{Symbol,PlacementReport}()
    path = Symbol[]
    function visit(name::Symbol)
        haskey(out, name) && return out[name]
        haskey(pending, name) || throw(ArgumentError("no frame named $name"))
        if name in path
            loop = path[findfirst(==(name), path):end]
            throw(ArgumentError("frames placed in each other in a loop: " *
                                join(string.([loop; name]), " → ")))
        end
        push!(path, name)
        d = pending[name]
        parent = visit(d.parent)
        for s in d.steps, r in _step_refs(s)
            visit(r)
        end
        pop!(path)
        f, reports[name] = _apply_steps(name, compose(parent, euler_frame(d.offset, d.rotation...)),
                                        d.steps, out)
        out[name] = f
        return f
    end
    for name in keys(pending)
        visit(name)
    end
    return out, reports
end
