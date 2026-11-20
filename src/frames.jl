# Rigid frames are defined in this file, each a rotation and an origin in the chamber. The
# file also contains the frames of ports, foils and target parts, coordinates in a frame,
# and a named graph of frames placed in each other.

"""
    Frame(R, origin)
    Frame(origin; z, x = θ̂ of z)

A right-handed rigid frame in the chamber. The columns of `R` are the unit axes of the
frame, and `origin` is its origin (mm). A point with local coordinates `x` therefore lies at
`R * x + origin` in the chamber ([`to_chamber`](@ref)), and a chamber point `p` has local
coordinates `R' * (p - origin)` ([`to_frame`](@ref)). [`CHAMBER`](@ref) is the identity
frame at target-chamber center (TCC).

The keyword form builds the frame from its `z` axis, of any length, and a vector `x` that
sets the direction of its `x` axis. The `x` axis is `x` minus its component along `z`,
normalized. By default `x` is the polar unit vector θ̂ of `z` (`port_basis`), which lies
along chamber `x` when `z` is the pole. `y = z × x` completes the right-handed triad. An
`x` parallel to `z` throws an `ArgumentError`.

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

Returns the chamber coordinates `R * x + origin` of the local point `x` of `frame`. With
`direction = true`, it returns those of the local vector `x`, `R * x`, without
translation. The element type follows `x`, so duals pass through.
"""
to_chamber(f::Frame, x::AbstractVector; direction::Bool = false) =
    direction ? f.R * SVector{3}(x[1], x[2], x[3]) : f.R * SVector{3}(x[1], x[2], x[3]) + f.origin
to_chamber(f::Frame, x::Tuple; direction::Bool = false) = to_chamber(f, SVector(x...); direction)

"""
    to_frame(frame, p; direction = false) -> SVector{3}

Returns the local coordinates `R' * (p - origin)` in `frame` of the chamber point `p`.
With `direction = true`, it returns those of the chamber vector `p`, `R' * p`. It is the
inverse of [`to_chamber`](@ref).
"""
to_frame(f::Frame, p::AbstractVector; direction::Bool = false) =
    direction ? f.R' * SVector{3}(p[1], p[2], p[3]) :
    f.R' * (SVector{3}(p[1], p[2], p[3]) - f.origin)
to_frame(f::Frame, p::Tuple; direction::Bool = false) = to_frame(f, SVector(p...); direction)

"""
    compose(parent, child) -> Frame

Returns the chamber frame of `child`, a frame given in the local coordinates of
`parent`, so that
`to_chamber(compose(parent, child), x) == to_chamber(parent, to_chamber(child, x))`
within round-off.
"""
compose(a::Frame, b::Frame) = Frame(a.R * b.R, a.R * b.origin + a.origin)

"""
    inv(frame) -> Frame

Returns the inverse rigid map: `to_chamber(inv(f), p) == to_frame(f, p)` (within round-off),
and `compose(f, inv(f))` is the identity.
"""
Base.inv(f::Frame) = Frame(f.R', -(f.R' * f.origin))

# Return the right-handed rotations about z and y by `a` degrees.
_Rz(a) = SMatrix{3,3,Float64}(cosd(a), sind(a), 0, -sind(a), cosd(a), 0, 0, 0, 1)
_Ry(a) = SMatrix{3,3,Float64}(cosd(a), 0, -sind(a), 0, 1, 0, sind(a), 0, cosd(a))

"""
    euler_rotation(polar, azimuth, rotation = 0) -> SMatrix{3,3}

Returns the rotation `Rz(azimuth) · Ry(polar) · Rz(rotation)` (degrees, right-handed). The
rotation is a turn by `rotation` about `z`, a tilt by `polar` about `y`, and then a turn by
`azimuth` about `z`. Its `z` column is `port_vector(polar, azimuth)`. At `rotation = 0` its
`x` and `y` columns are the polar and azimuthal unit vectors θ̂ and φ̂ at that direction
(`port_basis`). A frame with this rotation therefore has its `z` axis along the direction
(polar, azimuth), and `rotation` turns its `x` axis from θ̂ toward φ̂. This is the order of
the VisRad orientation angles Polar, Azimuthal and Rotation, which part placement in the GUI
also uses.
"""
euler_rotation(polar, azimuth, rotation = 0.0) = _Rz(azimuth) * _Ry(polar) * _Rz(rotation)

"""
    euler_frame(origin, polar, azimuth, rotation = 0) -> Frame

Returns the frame at `origin` with axes
[`euler_rotation`](@ref)`(polar, azimuth, rotation)`.
"""
euler_frame(origin, polar, azimuth, rotation = 0.0) =
    Frame(euler_rotation(polar, azimuth, rotation), Vec3(origin...))

"""
    port_frame(theta, phi; roll = 0, origin = (0, 0, 0)) -> Frame

    port_frame(beam::Beam; roll = 0, origin = (0, 0, 0)) -> Frame
    port_frame(facility::Symbol, port; roll = 0, origin = (0, 0, 0)) -> Frame

Returns the frame of the port at polar angle `theta` and azimuth `phi` (degrees), with
its origin at `origin` (TCC by default):

- `z` lies along `port_vector(theta, phi)`, toward the port. A beam from the port
  therefore travels along `-z`, as the beam convention `k = -port` requires.
- `x` is the port θ̂ turned by `roll` (degrees) toward φ̂.

The frame equals [`euler_frame`](@ref)`(origin, theta, phi, roll)`. The second form takes
the port of a beam ([`port_angles`](@ref)), and the third takes a facility port by name
([`port_angles`](@ref)`(facility, port)`). According to the documentation of the VisRad
coordinate transform tool, its diagnostic-port frames share the origin and `z`. The `x`
axis along θ̂ is an assumption at low confidence, because VisRad port frames have not been
seen in a sample workspace (research notes, "Frames and placement (2026-10-28)").
"""
port_frame(θ, φ; roll = 0.0, origin = (0.0, 0.0, 0.0)) = euler_frame(origin, θ, φ, roll)
port_frame(b::Beam; roll = 0.0, origin = (0.0, 0.0, 0.0)) = port_frame(port_angles(b)...; roll, origin)
port_frame(facility::Symbol, port; roll = 0.0, origin = (0.0, 0.0, 0.0)) =
    port_frame(port_angles(facility, port)...; roll, origin)

"""
    port_angles(facility::Symbol, port) -> (theta, phi)

Returns the polar angle and azimuth (degrees) of a facility port from the facility table:
- `:nif`: an indirect-drive quad port, by its number (`NIF_PORTS`, NIF User Guide
  Table 4-1) or its quad name such as `"Q31B"` ([`nif_port`](@ref)), or a diagnostic
  port by name such as `"P143-274"` ([`NIF_DIAGNOSTIC_PORTS`](@ref), Table 5-1);
- `:omega`: a beam port by name, `"B01"`–`"B60"`, in the hohlraum frame of
  [`omega_ports`](@ref).
An unknown port throws a `KeyError`, and another facility throws an `ArgumentError`.
"""
function port_angles(facility::Symbol, port)
    if facility === :nif
        d = port isa Integer ? nothing : _nif_diagnostic_port(String(port))
        d === nothing || return (d[2], d[3])
        n = port isa Integer ? Int(port) : nif_port(String(port))
        k = findfirst(r -> r[1] == n, NIF_PORTS)
        k === nothing && throw(KeyError(port))
        return (Float64(NIF_PORTS[k][2]), Float64(NIF_PORTS[k][3]))
    elseif facility === :omega
        k = findfirst(r -> r.name == String(port), omega_ports())
        k === nothing && throw(KeyError(port))
        r = omega_ports()[k]
        return (r.theta, r.phi)
    end
    throw(ArgumentError("facility must be :nif or :omega"))
end

"""
    transform_point(from::Frame, to::Frame, c; input = :cartesian, output = :cartesian,
                    direction = false) -> NTuple{3,Float64}
    transform_point(frames, from::Symbol, to::Symbol, c; …)

Returns the coordinates in frame `to` of the point with coordinates `c` in frame `from`,
after the VisRad coordinate transform tool. `c` is given in the system `input` and the
result in `output`, each `:cartesian`, `:cylindrical` or `:spherical`
([`coords_of`](@ref), [`point_of`](@ref)), in mm and degrees. With `direction` the three
numbers are a vector, which is turned by the frame rotations and not moved by the frame
origins. Its cylindrical or spherical coordinates are those of the vector tip from the
origin. The second form looks the frames up by name in `frames`, a dictionary such as
[`resolve_frames`](@ref) returns, and a missing name throws an `ArgumentError`.
"""
function transform_point(from::Frame, to::Frame, c; input::Symbol = :cartesian,
                         output::Symbol = :cartesian, direction::Bool = false)
    direction || return coords_of(point_of(c, input; frame = from), output; frame = to)
    v = to_chamber(from, point_of(c, input); direction = true)
    return coords_of(to_frame(to, v; direction = true), output)
end

function transform_point(frames::AbstractDict, from::Symbol, to::Symbol, c; kwargs...)
    F, T = get(frames, from, nothing), get(frames, to, nothing)
    F === nothing && throw(ArgumentError("no frame named $from"))
    T === nothing && throw(ArgumentError("no frame named $to"))
    return transform_point(F, T, c; kwargs...)
end

"""
    foil_frame(los; standoff, tilt = 0.0, tilt_toward = 0.0, rotation = 0.0,
               center = (0, 0, 0)) -> Frame

Returns the frame of a backlighter foil on the line of sight `los` = (θ, φ) in degrees,
as [`backlighter_foil`](@ref) places it:

- The origin is at the foil center, `center − standoff · port_vector(θ, φ)`.
- `x` lies along the `width` edge of the foil, `y` along its `height` edge, and `z` along
  its front normal, toward the detector.

Before the tilt, `x` is θ̂ of the line of sight turned by `rotation` toward φ̂. `tilt`
then turns the whole frame off the line of sight toward `cos(ψ) θ̂ + sin(ψ) φ̂`, with
`ψ = tilt_toward` in degrees. Local `(u, v, 0)` is a point on the foil face, with `u`
across its width.
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

# Rotates `v` by `tilt` (rad) about the axis los × e. This rotation turns the unit vector
# `los` toward `e`, a unit vector ⟂ los.
function _tilted(v::Vec3, los::Vec3, e::Vec3, tilt)
    k = cross(los, e)                       # rotation axis (unit: los ⟂ e)
    return cos(tilt) * v + sin(tilt) * cross(k, v) + (1 - cos(tilt)) * dot(k, v) * k
end

"""
    part_centroid(target, part) -> Vec3
    part_centroid(target) -> Dict{Symbol,Vec3}

Returns the area-weighted centroid of the triangles of `part`, over all its pieces (see
[`part_pieces`](@ref)). For parts such as foils, disks and entrance holes the centroid
is the center, a reference point to aim beams at. The second form gives the centroid of
every part of the target.
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

Returns a frame for a part of a target built in a script. Parts made in the GUI have
separate placement frames. The origin is at the [`part_centroid`](@ref) of the part, and
the `z` axis is the first of these that applies:

1. The axis of its exact surface, such as the axis of a `CylinderSurface`,
   `ConeSurface` or `TorusSurface` or the normal of a `DiskSurface`. This rule applies
   when every piece with such a surface has a parallel one, and the axis of the first
   piece sets the sign.
2. The area-weighted mean front normal of its triangles, when that mean is at least 0.1
   in length. For a flat foil this is its front normal; for a hemispherical cup the
   mean is 0.5 long.
3. The principal axis of its triangle area moments whose moment differs most from the
   other two, which is the symmetry axis of a body of revolution. The axis is signed
   toward chamber `+z`, else `+x`, else `+y`. When the three moments agree within 1e-6,
   as for spheres and cubes, the axis is chamber `z`.

The `x` axis is θ̂ of `z` (see [`Frame`](@ref)). The choice is ambiguous for shapes with
no axis, such as spheres and closed boxes, and for the sign of a principal axis. The 0.1
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
    # Sum the area moments about the centroid. The exact second moment of each triangle
    # about the origin o, from its vertices and centroid, is A/12 (Σ vᵢvᵢᵀ + 9 c cᵀ).
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
    # Pick the eigenvalue farthest from the mean of the other two.
    j = argmax(i -> abs(λ[i] - (sum(λ) - λ[i]) / 2), 1:3)
    z = Vec3(V[:, j]...)
    s = abs(z[3]) > 1e-12 ? z[3] : abs(z[1]) > 1e-12 ? z[1] : z[2]
    return Frame(o; z = s < 0 ? -z : z)
end

"""
    coords_of(p, system = :cartesian; frame = CHAMBER) -> NTuple{3,Float64}

Returns the coordinates of the chamber point `p` in `frame`, in one of three systems:

- `:cartesian` (the default) gives (x, y, z).
- `:cylindrical` gives (r, φ, z): the distance from the frame `z` axis, the azimuth in
  degrees in [0, 360) from `x` toward `y`, and the height.
- `:spherical` gives (r, θ, φ): the distance from the origin, the polar angle from `z`
  and the azimuth, in degrees.

At r = 0 the angles are 0. [`point_of`](@ref) is the inverse.
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

Returns the chamber point with coordinates `c` in `frame`, in the order of
[`coords_of`](@ref): `:cartesian` (x, y, z), `:cylindrical` (r, φ, z) or `:spherical`
(r, θ, φ), with angles in degrees. A negative radius throws an `ArgumentError`.
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

A frame placed in another frame. Its origin lies at `offset` (mm) in the frame named
`parent`, and its axes are turned by `rotation` = (polar, azimuth, rotation) in degrees in
the parent ([`euler_rotation`](@ref)). Its chamber frame is therefore
`compose(parent's frame, euler_frame(offset, rotation...))`. The [`PlacementStep`](@ref)s
`steps` then apply in order, as relations to other frames and relative moves. Degrees of
freedom that no step fixes keep the values of `offset` and `rotation`. The frame of a port
is `FrameDef(:tcc; rotation = (θ, φ, roll))` ([`port_frame`](@ref)). Resolve a set of named
definitions with [`resolve_frames`](@ref).
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

Returns the chamber frames of the named definitions `defs`, given as a dictionary or a
collection of `name => FrameDef` pairs. Each frame is placed in its parent and then moved by
its placement steps. The parent is another definition or one of the `base` frames, which are
analytic frames such as TCC, the entrance holes of a template and foils. `:tcc` is added as
[`CHAMBER`](@ref) when absent, and the result also contains the base frames too. Definitions
are resolved depth first, and the frames named by the steps of a definition are resolved
before it. A frame therefore moves exactly with every ancestor and with every frame that its
steps refer to. Each of these is an `ArgumentError`:

- a definition that reuses the name of a base frame, or a name given twice;
- an unknown parent or step frame;
- frames that depend on each other in a loop, where the error names the loop, e.g.
  `a → b → a`;
- a step that cannot apply, such as a facing point at the origin;
- a step that undoes an earlier one (see [`placement_report`](@ref)).
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
