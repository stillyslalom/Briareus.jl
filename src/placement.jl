# Placement steps: relations to other frames (coincident, on an axis, concentric, at a
# distance, facing, parallel, rolled) and relative moves, applied in order to a frame
# after its offset and rotation, with the degrees of freedom each step fixes and
# conflicts between steps.

"""
    PlacementStep

A step of a frame's placement ([`FrameDef`](@ref)'s `steps`), applied in order after the
frame's offset and rotation in its parent. Relations fix degrees of freedom:
[`Coincident`](@ref), [`OnAxis`](@ref), [`Concentric`](@ref), [`Distance`](@ref),
[`Facing`](@ref), [`Parallel`](@ref) and [`Roll`](@ref) with a reference direction.
Relative moves fix none: [`Offset`](@ref), [`Rotate`](@ref) and [`Roll`](@ref) without
one. Each step has a closed form; a rotation turns the frame about its own origin.
Frames named by steps (`ref`) are resolved first, so a frame follows them as it follows
its parent. [`placement_report`](@ref) gives what each step fixes.
"""
PlacementStep

# Unit vector of an axis given as :x, :y, :z or a 3-vector (local coordinates).
function _step_axis(axis)
    axis === :x && return Vec3(1, 0, 0)
    axis === :y && return Vec3(0, 1, 0)
    axis === :z && return Vec3(0, 0, 1)
    axis isa Symbol && throw(ArgumentError("axis must be :x, :y, :z or a 3-vector, not :$axis"))
    v = Vec3(axis...)
    n = norm(v)
    n > 0 && isfinite(n) || throw(ArgumentError("an axis must be a nonzero vector"))
    return v / n
end

"""
    Coincident(ref; point = (0, 0, 0))

Move the frame's origin onto the local `point` of the frame `ref` (mm). Fixes the three
translations.
"""
struct Coincident <: PlacementStep
    ref::Symbol
    point::Vec3
end
Coincident(ref::Symbol; point = (0.0, 0.0, 0.0)) = Coincident(ref, Vec3(point...))

"""
    OnAxis(ref; axis = :z, distance = nothing)

Move the frame's origin onto the line through `ref`'s origin along its `axis` (`:x`,
`:y`, `:z` or a vector in `ref`'s coordinates): with `distance` (mm), to that signed
distance from `ref`'s origin along the axis, fixing the three translations; without,
to the nearest point of the line, fixing the two translations across it.
"""
struct OnAxis <: PlacementStep
    ref::Symbol
    axis::Vec3
    distance::Union{Nothing,Float64}
end
OnAxis(ref::Symbol; axis = :z, distance = nothing) =
    OnAxis(ref, _step_axis(axis), distance === nothing ? nothing : Float64(distance))

"""
    Concentric(ref; axis = :z, flip = false)

Turn the frame so its `z` axis lies along `ref`'s `axis` (against it with `flip`), by the
smallest rotation, and move its origin onto that axis line at the nearest point. Fixes
the two translations and the two rotations across the axis; the position along it and
the roll about it stay free.
"""
struct Concentric <: PlacementStep
    ref::Symbol
    axis::Vec3
    flip::Bool
end
Concentric(ref::Symbol; axis = :z, flip::Bool = false) = Concentric(ref, _step_axis(axis), flip)

"""
    Distance(ref, distance; point = (0, 0, 0))

Move the frame's origin along the line from `ref`'s local `point` through it, to
`distance` (mm, ≥ 0) from the point. Fixes the translation along that line. An origin
at the point is an `ArgumentError` (no line).
"""
struct Distance <: PlacementStep
    ref::Symbol
    distance::Float64
    point::Vec3
end
function Distance(ref::Symbol, distance::Real; point = (0.0, 0.0, 0.0))
    distance >= 0 || throw(ArgumentError("a distance must not be negative"))
    return Distance(ref, Float64(distance), Vec3(point...))
end

"""
    Facing(ref; point = (0, 0, 0), side = :front)

Turn the frame by the smallest rotation so its `z` axis points at `ref`'s local `point`
(`side = :front`) or away from it (`side = :back`). Fixes the two rotations across that
direction; the roll about it stays free. A point at the frame's origin is an
`ArgumentError`.
"""
struct Facing <: PlacementStep
    ref::Symbol
    point::Vec3
    side::Symbol
end
function Facing(ref::Symbol; point = (0.0, 0.0, 0.0), side::Symbol = :front)
    side in (:front, :back) || throw(ArgumentError("side must be :front or :back"))
    return Facing(ref, Vec3(point...), side)
end

"""
    Parallel(ref; axis = :z, flip = false)

Turn the frame by the smallest rotation so its `z` axis is parallel to `ref`'s `axis`
(antiparallel with `flip`). Fixes the two rotations across the axis.
"""
struct Parallel <: PlacementStep
    ref::Symbol
    axis::Vec3
    flip::Bool
end
Parallel(ref::Symbol; axis = :z, flip::Bool = false) = Parallel(ref, _step_axis(axis), flip)

"""
    Roll(angle; from = nothing)

Turn the frame about its own `z` axis. With `from` (a frame name, its `x` axis, or
`(name, axis)`), set the frame's `x` axis to `angle` degrees from that direction's
projection across `z`, turning toward `y`; this fixes the roll. A direction along `z`
has no projection and is an `ArgumentError`. Without `from`, turn by `angle` from the
current roll, a relative move.
"""
struct Roll <: PlacementStep
    angle::Float64
    ref::Union{Nothing,Symbol}
    axis::Vec3
end
function Roll(angle::Real; from = nothing)
    from === nothing && return Roll(Float64(angle), nothing, Vec3(1, 0, 0))
    from isa Symbol && return Roll(Float64(angle), from, Vec3(1, 0, 0))
    return Roll(Float64(angle), Symbol(from[1]), _step_axis(from[2]))
end

"""
    Offset(dx, dy, dz; frame = nothing)

Move the frame by `(dx, dy, dz)` (mm) along its own axes, or along the axes of the frame
named `frame`. A relative move: it fixes nothing.
"""
struct Offset <: PlacementStep
    delta::Vec3
    ref::Union{Nothing,Symbol}
end
Offset(dx::Real, dy::Real, dz::Real; frame::Union{Nothing,Symbol} = nothing) =
    Offset(Vec3(dx, dy, dz), frame)

"""
    Rotate(axis, angle)

Turn the frame by `angle` degrees (right-handed) about its own `axis` (`:x`, `:y`, `:z`
or a vector in its coordinates) through its origin. A relative move: it fixes nothing.
"""
struct Rotate <: PlacementStep
    axis::Vec3
    angle::Float64
end
Rotate(axis, angle::Real) = Rotate(_step_axis(axis), Float64(angle))

_step_refs(s::PlacementStep) = s.ref === nothing ? () : (s.ref,)
_step_refs(::Rotate) = ()
_relative(::PlacementStep) = false
_relative(::Union{Offset,Rotate}) = true
_relative(s::Roll) = s.ref === nothing

function _step_label(k, s::PlacementStep)
    refs = _step_refs(s)
    return "step $k ($(nameof(typeof(s)))" * (isempty(refs) ? "" : " $(only(refs))") * ")"
end

# Angle between two vectors (degrees).
_angle_between(a, b) = atand(norm(cross(a, b)), dot(a, b))

# Rotation by `angle` (rad) about the unit vector `k` (Rodrigues).
function _axis_rotation(k::Vec3, angle)
    K = SMatrix{3,3,Float64}(0, k[3], -k[2], -k[3], 0, k[1], k[2], -k[1], 0)
    return cos(angle) * SMatrix{3,3,Float64}(I) + sin(angle) * K + (1 - cos(angle)) * (k * k')
end

# The frame turned about its origin by the smallest rotation taking its z axis to the
# unit vector `t`; z opposite to `t` turns by 180° about the frame's x axis.
function _turn_z(f::Frame, t::Vec3)
    z = f.R[:, 3]
    c = cross(z, t)
    s, d = norm(c), dot(z, t)
    s <= 1e-15 && d > 0 && return f
    Q = s <= 1e-15 ? _axis_rotation(f.R[:, 1], π) : _axis_rotation(c / s, atan(s, d))
    return Frame(Q * f.R, f.origin)
end

_ref_axis(fr, ref::Symbol, axis::Vec3) = Vec3(to_chamber(fr[ref], axis; direction = true))
_ref_point(fr, ref::Symbol, p::Vec3) = Vec3(to_chamber(fr[ref], p))

const _AXES3 = (Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1))
_across(a::Vec3) = port_basis(a)

# Each step: the new frame, the chamber directions of the translations it fixes and the
# axes of the rotations it fixes.
function _place(s::Coincident, f::Frame, fr)
    return Frame(f.R, _ref_point(fr, s.ref, s.point)), _AXES3, ()
end

function _place(s::OnAxis, f::Frame, fr)
    p0, a = fr[s.ref].origin, _ref_axis(fr, s.ref, s.axis)
    s.distance === nothing || return Frame(f.R, p0 + s.distance * a), _AXES3, ()
    return Frame(f.R, p0 + dot(f.origin - p0, a) * a), _across(a), ()
end

function _place(s::Concentric, f::Frame, fr)
    p0, a = fr[s.ref].origin, _ref_axis(fr, s.ref, s.axis)
    g = _turn_z(f, s.flip ? -a : a)
    return Frame(g.R, p0 + dot(g.origin - p0, a) * a), _across(a), _across(a)
end

function _place(s::Distance, f::Frame, fr)
    q = _ref_point(fr, s.ref, s.point)
    r = f.origin - q
    n = norm(r)
    n > 1e-12 * max(1.0, norm(q)) ||
        throw(ArgumentError("the origin is at the point, so the line to it is undefined"))
    u = r / n
    return Frame(f.R, q + s.distance * u), (u,), ()
end

function _place(s::Facing, f::Frame, fr)
    t = _ref_point(fr, s.ref, s.point) - f.origin
    n = norm(t)
    n > 1e-12 * max(1.0, norm(f.origin)) ||
        throw(ArgumentError("the point is at the frame's origin, so it has no direction"))
    t = s.side === :front ? t / n : -t / n
    return _turn_z(f, t), (), _across(t)
end

function _place(s::Parallel, f::Frame, fr)
    a = _ref_axis(fr, s.ref, s.axis)
    a = s.flip ? -a : a
    return _turn_z(f, a), (), _across(a)
end

function _place(s::Roll, f::Frame, fr)
    z = f.R[:, 3]
    s.ref === nothing &&
        return Frame(_axis_rotation(z, deg2rad(s.angle)) * f.R, f.origin), (), ()
    e = _roll_x(s, z, fr)
    x = cosd(s.angle) * e + sind(s.angle) * cross(z, e)
    return Frame(hcat(x, cross(z, x), z), f.origin), (), (z,)
end

function _roll_x(s::Roll, z::Vec3, fr)
    r = _ref_axis(fr, s.ref, s.axis)
    p = r - dot(r, z) * z
    norm(p) > 1e-9 || throw(ArgumentError("the reference direction lies along the frame's z axis"))
    return normalize(p)
end

function _place(s::Offset, f::Frame, fr)
    R = s.ref === nothing ? f.R : fr[s.ref].R
    return Frame(f.R, f.origin + R * s.delta), (), ()
end

_place(s::Rotate, f::Frame, fr) =
    Frame(_axis_rotation(f.R * s.axis, deg2rad(s.angle)) * f.R, f.origin), (), ()

# How far the frame is from satisfying a step: (mm, degrees). Relative moves: zero.
_residual(::PlacementStep, f::Frame, fr) = (0.0, 0.0)
_residual(s::Coincident, f::Frame, fr) = (norm(f.origin - _ref_point(fr, s.ref, s.point)), 0.0)
function _residual(s::OnAxis, f::Frame, fr)
    p0, a = fr[s.ref].origin, _ref_axis(fr, s.ref, s.axis)
    s.distance === nothing && return (norm(cross(f.origin - p0, a)), 0.0)
    return (norm(f.origin - (p0 + s.distance * a)), 0.0)
end
function _residual(s::Concentric, f::Frame, fr)
    p0, a = fr[s.ref].origin, _ref_axis(fr, s.ref, s.axis)
    return (norm(cross(f.origin - p0, a)), _angle_between(f.R[:, 3], s.flip ? -a : a))
end
_residual(s::Distance, f::Frame, fr) =
    (abs(norm(f.origin - _ref_point(fr, s.ref, s.point)) - s.distance), 0.0)
function _residual(s::Facing, f::Frame, fr)
    t = _ref_point(fr, s.ref, s.point) - f.origin
    return (0.0, _angle_between(f.R[:, 3], s.side === :front ? t : -t))
end
function _residual(s::Parallel, f::Frame, fr)
    a = _ref_axis(fr, s.ref, s.axis)
    return (0.0, _angle_between(f.R[:, 3], s.flip ? -a : a))
end
function _residual(s::Roll, f::Frame, fr)
    s.ref === nothing && return (0.0, 0.0)
    z = f.R[:, 3]
    e = _roll_x(s, z, fr)
    return (0.0, _angle_between(f.R[:, 1], cosd(s.angle) * e + sind(s.angle) * cross(z, e)))
end

# Tolerances on a step's residual (mm, degrees): satisfied, and disturbed by a later step.
const _PLACE_TOL = (1e-9, 1e-7)

# Add the directions `vs` to the orthonormal basis `B`; the number added.
function _extend!(B::Vector{Vec3}, vs)
    n = 0
    for v in vs
        w = v - sum((dot(v, b) * b for b in B); init = zero(Vec3))
        norm(w) > 1e-9 || continue
        push!(B, normalize(w))
        n += 1
    end
    return n
end

# An orthonormal basis of the directions not in the orthonormal basis `B`.
function _complement(B::Vector{Vec3})
    C = copy(B)
    _extend!(C, _AXES3)
    return C[length(B)+1:end]
end

"""
    PlacementReport

What a frame's placement steps fix ([`placement_report`](@ref)). Fields: `frame` (its
name), `steps`, `fixed` (per step, the translations and rotations it fixes that no
earlier step fixed, as `(translations, rotations)`), `redundant` (indices of relations
that fix nothing new), `free_translation` and `free_rotation` (orthonormal chamber
directions and rotation axes no step fixes; they keep the values of the frame's offset
and rotation) and `residuals` (per step, how far the final frame is from it,
`(mm, degrees)`: zero for a relation no later relative move disturbs, and for every
relative move).

The rotations a step fixes are counted at the frame's orientation after the step (the
axes across a `z` axis set by [`Facing`](@ref), [`Parallel`](@ref) or
[`Concentric`](@ref), and that `z` axis for [`Roll`](@ref)).
"""
struct PlacementReport
    frame::Symbol
    steps::Vector{PlacementStep}
    fixed::Vector{NTuple{2,Int}}
    redundant::Vector{Int}
    free_translation::Vector{Vec3}
    free_rotation::Vector{Vec3}
    residuals::Vector{NTuple{2,Float64}}
end

function Base.show(io::IO, r::PlacementReport)
    nt, nr = length(r.free_translation), length(r.free_rotation)
    print(io, "PlacementReport(", r.frame, ": ", length(r.steps), " steps, free: ", nt,
          " translation", nt == 1 ? "" : "s", ", ", nr, " rotation", nr == 1 ? "" : "s")
    isempty(r.redundant) || print(io, "; redundant: ", join(r.redundant, ", "))
    print(io, ")")
end

# Apply `steps` in order to the frame `f` of the frame `name`; the placed frame and its
# report. A relation that disturbs an earlier one (raises its residual past the
# tolerance) is an `ArgumentError` naming both steps and the size of the change.
function _apply_steps(name::Symbol, f::Frame, steps::AbstractVector{<:PlacementStep}, fr)
    T, Rb = Vec3[], Vec3[]
    fixed, redundant = NTuple{2,Int}[], Int[]
    res = NTuple{2,Float64}[]
    for (k, s) in enumerate(steps)
        g, tf, rf = try
            _place(s, f, fr)
        catch err
            err isa ArgumentError || rethrow()
            throw(ArgumentError("frame $name, $(_step_label(k, s)): $(err.msg)"))
        end
        now = [_residual(steps[j], g, fr) for j in 1:k-1]
        if !_relative(s)
            for j in 1:k-1
                _relative(steps[j]) && continue
                dmm, ddeg = now[j] .- res[j]
                (dmm > _PLACE_TOL[1] || ddeg > _PLACE_TOL[2]) || continue
                what = dmm > _PLACE_TOL[1] ? @sprintf("moves it by %.4g mm", dmm) :
                       @sprintf("turns it by %.4g°", ddeg)
                throw(ArgumentError("frame $name: $(_step_label(k, s)) conflicts with " *
                                    "$(_step_label(j, steps[j])), which it undoes ($what); " *
                                    "steps apply in order"))
            end
        end
        nt, nr = _extend!(T, tf), _extend!(Rb, rf)
        push!(fixed, (nt, nr))
        !_relative(s) && nt + nr == 0 && push!(redundant, k)
        res = push!(now, _residual(s, g, fr))
        f = g
    end
    return f, PlacementReport(name, collect(PlacementStep, steps), fixed, redundant,
                              _complement(T), _complement(Rb), res)
end

"""
    placement_report(defs, base = Dict(:tcc => CHAMBER)) -> Dict{Symbol,PlacementReport}

A [`PlacementReport`](@ref) for every definition in `defs` (as in
[`resolve_frames`](@ref), which raises the same errors): the degrees of freedom each
placement step fixes, redundant steps, what stays free, and the residuals.
"""
placement_report(defs, base::AbstractDict = Dict(:tcc => CHAMBER)) = last(_resolve_frames(defs, base))
