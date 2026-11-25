# This file covers the frames of the parts, with user and port frames, the frame a part is
# placed in, the placement steps, and moves and turns along the axes of a frame.

# ---------------------------------------------------------------- user frames

"Whether part `p` is a port frame: a `:frame` part tied to a port (its `object`)."
is_port_frame(p::PartSpec) = p.kind === :frame && !isempty(p.object)

"""
    user_frames(ts::TargetSetup) -> Vector{PartSpec}

The user frames of a design (its `:frame` parts, port frames included), in order.
"""
user_frames(ts::TargetSetup) = PartSpec[p for p in ts.parts[] if p.kind === :frame]

"""
    frame_users(ts::TargetSetup, name::Symbol) -> Vector{Symbol}

The parts and frames placed in the frame `name` (or a piece frame of part `name`), or
with a placement step that refers to it. When the frame is removed they keep their
place, in the nearest frame above that stays ([`remove_parts!`](@ref)); beams aimed in
a removed frame keep their aims, held in `:tcc`.
"""
function frame_users(ts::TargetSetup, name::Symbol)
    own(f) = _own_frame_name(name, f)
    return Symbol[q.name for q in ts.parts[] if q.name !== name &&
                  (own(q.frame) || any(st -> any(own, Briareus.GUIAPI.step_references(st)), q.steps))]
end

"""
    add_frame!(ts::TargetSetup; name = nothing, frame = :tcc, values...) -> Bool

Add a user frame (a `:frame` part, see [`PartSpec`](@ref)) named `name` (default
`frame_1`, …) placed in the frame `frame` by the placement fields given as keywords
(`x`, `y`, `z` in mm, `theta`, `phi`, `rotation` in degrees, as a part's). Parts can be
placed in it and beams aimed in it like any part's frame. A rejected frame sets `error`
and returns `false`.
"""
add_frame!(ts::TargetSetup; name = nothing, frame::Symbol = :tcc, values...) =
    add_part!(ts, :frame; name = name === nothing ? _free_name(ts, :frame) : name, frame, values...)

"""
    add_port_frame!(ts::TargetSetup, port, theta, phi; roll = 0, name = nothing) -> Bool

Add the frame of the port named `port` at polar angle `theta` and azimuth `phi`
(degrees): `Briareus.port_frame(theta, phi; roll)`, at TCC with `z` toward the port (a beam from it
travels along `-z`). Its name defaults to `port_<port>`. Only its roll can be edited
afterward; another frame placed in it carries offsets along the port axis.
"""
function add_port_frame!(ts::TargetSetup, port::AbstractString, θ::Real, φ::Real; roll::Real = 0.0,
                         name = nothing)
    nm = name === nothing ? _free_name(ts, Symbol(:port_, replace(String(port), r"[^A-Za-z0-9_]" => "_"));
                                       numbered = false) : name
    part = try
        isempty(port) && throw(ArgumentError("a port frame needs its port's name"))
        PartSpec(:frame; name = nm, object = port,
                 values = Dict(:theta => Float64(θ), :phi => Float64(φ), :rotation => Float64(roll)))
    catch err
        err isa ArgumentError || rethrow()
        ts.error[] = _errmsg(err)
        return false
    end
    return _try_build!(ts; parts = [ts.parts[]; part])
end

"""
    set_part_frame!(ts::TargetSetup, i, frame::Symbol; keep = true) -> Bool

Place part `i` in another frame (`:tcc`, another part's or a hohlraum part's piece
frame: see [`part_frames`](@ref)). With `keep` the part stays where it is in the
chamber: the part frame is re-expressed in the new one (its center, and the angles of
its rotation there; the angles stay as they are when the two frames have the same
axes); otherwise its placement values are kept and it moves
with the frame. A frame that leads back to the part is a loop: refused, with `error`
set.
"""
function set_part_frame!(ts::TargetSetup, i::Integer, frame::Symbol; keep::Bool = true)
    _check_part_index(ts, i)
    p = ts.parts[][i]
    p.frame === frame && return true
    try
        _own_frame_name(p.name, frame) && throw(ArgumentError("a part cannot be placed in the frame named after it"))
        trial = copy(ts.parts[])
        trial[i] = _with(p; frame)
        part_frames(trial)                     # a loop, a missing frame
        ps = copy(ts.parts[])
        ps[i] = keep ? _kept_in(p, frame, part_frames(ts)) : _with(p; frame)
        return _try_build!(ts; parts = ps)
    catch err
        err isa ArgumentError || rethrow()
        ts.error[] = _errmsg(err)
        return false
    end
end

# Part `p` placed in the frame named `frame`, where it is in the chamber (`fr`: the
# design's chamber frames, `part_frames`): its center and the angles of its rotation
# re-expressed there (the angles kept when the two frames have the same axes).
function _kept_in(p::PartSpec, frame::Symbol, fr::AbstractDict)
    vals = copy(p.values)
    P0, P1 = fr[p.frame], fr[frame]
    vals[:x], vals[:y], vals[:z] = to_frame(P1, to_chamber(P0, part_center(p)))
    if P0.R != P1.R
        local_R = P1.R' * (P0.R * euler_rotation(p.values[:theta], p.values[:phi], p.values[:rotation]))
        vals[:theta], vals[:phi], vals[:rotation] = _euler_angles(local_R)
    end
    return _with(p; frame, values = vals)
end

"""
Kinds of placement step a part can have ([`add_part_step!`](@ref)), in menu order:
relations that fix degrees of freedom (`:coincident`, `:on_axis`, `:concentric`,
`:distance`, `:facing`, `:parallel`, `:roll`) and relative moves (`:offset`, `:rotate`).
"""
const PART_STEP_KINDS = (:coincident, :on_axis, :concentric, :distance, :facing, :parallel,
                         :roll, :offset, :rotate)

const _STEP_TYPES = Dict(:coincident => Coincident, :on_axis => OnAxis, :concentric => Concentric,
                         :distance => Distance, :facing => Facing, :parallel => Parallel,
                         :roll => Roll, :offset => Offset, :rotate => Rotate)

"Menu labels of [`PART_STEP_KINDS`](@ref)."
const PART_STEP_LABELS = Dict(:coincident => "Coincident", :on_axis => "On axis",
                              :concentric => "Concentric", :distance => "Distance",
                              :facing => "Facing", :parallel => "Parallel", :roll => "Roll",
                              :offset => "Offset", :rotate => "Rotate")

"""
    step_kind(step) -> Symbol

The [`PART_STEP_KINDS`](@ref) entry of a placement step.
"""
step_kind(st::PlacementStep) = only(k for (k, T) in _STEP_TYPES if st isa T)

# A step as a dictionary of its fields (the session file; `nothing` fields left out), and back.
function _step_state(st::PlacementStep)
    d = Dict{String,Any}("kind" => String(step_kind(st)))
    for f in fieldnames(typeof(st))
        v = getfield(st, f)
        v === nothing && continue
        d[String(f)] = v isa Symbol ? String(v) : v isa Vec3 ? collect(v) : v
    end
    return d
end

function _state_step(d::AbstractDict)
    kind = Symbol(d["kind"])
    haskey(_STEP_TYPES, kind) || throw(ArgumentError("unknown placement step \"$kind\""))
    T = _STEP_TYPES[kind]
    args = map(fieldnames(T)) do f
        v, FT = get(d, String(f), nothing), fieldtype(T, f)
        if v === nothing
            Nothing <: FT || throw(ArgumentError("a $kind step needs its $f"))
            return nothing
        end
        FT === Vec3 ? Vec3(Float64.(v)...) :
        FT <: Union{Nothing,Symbol} ? Symbol(v) :
        FT === Bool ? Bool(v) : Float64(v)
    end
    return _checked_step(T(args...))
end

# The checks of the steps' keyword constructors, for steps built from their fields.
function _checked_step(st::PlacementStep)
    st isa Distance && st.distance < 0 && throw(ArgumentError("a distance must not be negative"))
    st isa Facing && !(st.side in (:front, :back)) && throw(ArgumentError("side must be front or back"))
    for f in (:axis, :point, :delta)
        hasfield(typeof(st), f) || continue
        all(isfinite, getfield(st, f)) || throw(ArgumentError("the $f must be finite"))
    end
    return st
end

# A new step of `kind` on the frame `ref`, with values that leave the part's frame `F`
# where it is where the step allows (`fr`: the design's frames).
function _step_default(kind::Symbol, ref::Symbol, F::Frame, fr)
    R = fr[ref]
    z = F.R[:, 3]
    kind === :coincident && return Coincident(ref; point = to_frame(R, F.origin))
    kind === :on_axis && return OnAxis(ref)
    kind === :concentric && return Concentric(ref; flip = dot(z, R.R[:, 3]) < 0)
    kind === :distance && return Distance(ref, norm(F.origin - R.origin))
    kind === :facing && return Facing(ref; side = dot(z, R.origin - F.origin) < 0 ? :back : :front)
    kind === :parallel && return Parallel(ref; flip = dot(z, R.R[:, 3]) < 0)
    if kind === :roll
        for a in (:x, :y)
            r = to_chamber(R, Briareus.GUIAPI.step_axis(a); direction = true)
            e = r - dot(r, z) * z
            norm(e) > 1e-9 || continue
            e = normalize(e)
            return Roll(atand(dot(F.R[:, 1], cross(z, e)), dot(F.R[:, 1], e)); from = (ref, a))
        end
    end
    kind === :offset && return Offset(0.0, 0.0, 0.0)
    kind === :rotate && return Rotate(:z, 0.0)
    throw(ArgumentError("placement step kind must be one of $PART_STEP_KINDS"))
end

function _set_steps!(ts::TargetSetup, i::Integer, steps)
    ps = copy(ts.parts[])
    ps[i] = _with(ps[i]; steps = collect(PlacementStep, steps))
    return _try_build!(ts; parts = ps)
end

function _step_edit(f, ts::TargetSetup)
    try
        return f()
    catch err
        err isa ArgumentError || rethrow()
        ts.error[] = _errmsg(err)
        return false
    end
end

"""
    add_part_step!(ts::TargetSetup, i, kind; ref = the part's frame) -> Bool

Append a placement step of `kind` ([`PART_STEP_KINDS`](@ref)) to part `i`, referring to
the frame `ref`, with values that leave the part where it is where the step allows
(coincident at its origin's current place, the current distance, roll and flips); a
relative move starts at zero. A step that cannot apply, or that undoes an earlier one,
is refused with `error` set.
"""
function add_part_step!(ts::TargetSetup, i::Integer, kind::Symbol; ref::Union{Nothing,Symbol} = nothing)
    _check_part_index(ts, i)
    p = ts.parts[][i]
    return _step_edit(ts) do
        fr = part_frames(ts)
        r = something(ref, p.frame)
        r === p.name && throw(ArgumentError("a placement step cannot refer to the frame of its part"))
        haskey(fr, r) || throw(ArgumentError("no frame named $r"))
        _set_steps!(ts, i, [p.steps; _step_default(kind, r, fr[p.name], fr)])
    end
end

"""
    remove_part_step!(ts::TargetSetup, i, k) -> Bool
    move_part_step!(ts::TargetSetup, i, k, by) -> Bool

Remove placement step `k` of part `i`, or move it `by` places (−1 earlier, +1 later;
the order matters, since steps apply in sequence).
"""
function remove_part_step!(ts::TargetSetup, i::Integer, k::Integer)
    _check_part_index(ts, i)
    steps = copy(ts.parts[][i].steps)
    checkbounds(steps, k)
    deleteat!(steps, k)
    return _step_edit(() -> _set_steps!(ts, i, steps), ts)
end

function move_part_step!(ts::TargetSetup, i::Integer, k::Integer, by::Integer)
    _check_part_index(ts, i)
    steps = copy(ts.parts[][i].steps)
    checkbounds(steps, k)
    j = clamp(k + by, 1, length(steps))
    j == k && return true
    st = popat!(steps, k)
    insert!(steps, j, st)
    return _step_edit(() -> _set_steps!(ts, i, steps), ts)
end

"""
    edit_part_step!(ts::TargetSetup, i, k, key, value) -> Bool

Set a field of placement step `k` of part `i`: `:ref` (a frame name; `""` for none: the
part axes for an offset, a relative roll), `:px`/`:py`/`:pz` (the point in `ref`,
mm), `:dx`/`:dy`/`:dz` (an offset, mm), `:distance` (mm; `""` for none on an on-axis
step), `:angle` (degrees), `:axis` (`"x"`, `"y"` or `"z"`), `:flip` (Bool) or `:side`
(`"front"` or `"back"`); numbers as text. A rejected value sets `error`.
"""
function edit_part_step!(ts::TargetSetup, i::Integer, k::Integer, key::Symbol, value)
    _check_part_index(ts, i)
    p = ts.parts[][i]
    checkbounds(p.steps, k)
    st = p.steps[k]
    return _step_edit(ts) do
        d = _step_state(st)
        kind = d["kind"]
        coord = _STEP_COORDS
        if key === :ref
            s = strip(string(value))
            hasfield(typeof(st), :ref) || throw(ArgumentError("a $kind step has no frame"))
            isempty(s) && !(st isa Union{Roll,Offset}) && throw(ArgumentError("a $kind step needs a frame"))
            Symbol(s) === p.name && throw(ArgumentError("a placement step cannot refer to the frame of its part"))
            isempty(s) ? delete!(d, "ref") : (d["ref"] = s)
        elseif haskey(coord, key)
            f, c = coord[key]
            haskey(d, f) || throw(ArgumentError("a $kind step has no $f"))
            d[f][c] = _parse_number(value, String(key))
        elseif key === :distance && st isa OnAxis && isempty(strip(string(value)))
            delete!(d, "distance")
        elseif key in (:distance, :angle)
            hasfield(typeof(st), key) || throw(ArgumentError("a $kind step has no $key"))
            d[String(key)] = _parse_number(value, String(key))
        elseif key === :axis
            haskey(d, "axis") || throw(ArgumentError("a $kind step has no axis"))
            d["axis"] = collect(Briareus.GUIAPI.step_axis(Symbol(value)))
        elseif key === :flip
            haskey(d, "flip") || throw(ArgumentError("a $kind step has no flip"))
            d["flip"] = Bool(value)
        elseif key === :side
            haskey(d, "side") || throw(ArgumentError("a $kind step has no side"))
            d["side"] = String(value)
        else
            throw(ArgumentError("unknown placement step field :$key"))
        end
        new = _state_step(d)
        _step_state(new) == _step_state(st) && (ts.error[] = ""; return true)
        steps = copy(p.steps)
        steps[k] = new
        _set_steps!(ts, i, steps)
    end
end

"""
    part_placement(ts::TargetSetup, i) -> PlacementReport

What part `i`'s placement steps fix (`placement_report`): the degrees of freedom per
step, redundant steps, what stays free, and the residuals left by relative moves.
"""
function part_placement(ts::TargetSetup, i::Integer)
    _check_part_index(ts, i)
    ps = ts.parts[]
    defs = [d for p in ps for d in (p.name => part_frame_def(p), piece_frame_defs(p)...)]
    return placement_report(defs)[ps[i].name]
end

# The fields of a step's point and offset coordinates: (field, index).
const _STEP_COORDS = Dict(:px => ("point", 1), :py => ("point", 2), :pz => ("point", 3),
                          :dx => ("delta", 1), :dy => ("delta", 2), :dz => ("delta", 3))

"""
    step_numbers(step) -> Vector{Tuple{Symbol,String,String}}

The numbers of a placement step, as `(key, label, unit)`, keys as in
[`edit_part_step!`](@ref): a point's or an offset's coordinates (mm), a distance (mm;
an on-axis step's only when set) and an angle (degrees). A sweep can vary each
([`set_sweep_step!`](@ref)).
"""
function step_numbers(st::PlacementStep)
    out = Tuple{Symbol,String,String}[]
    hasfield(typeof(st), :point) && append!(out, [(:px, "point x", "mm"), (:py, "point y", "mm"), (:pz, "point z", "mm")])
    hasfield(typeof(st), :delta) && append!(out, [(:dx, "offset x", "mm"), (:dy, "offset y", "mm"), (:dz, "offset z", "mm")])
    (st isa Distance || st isa OnAxis && st.distance !== nothing) && push!(out, (:distance, "distance", "mm"))
    hasfield(typeof(st), :angle) && push!(out, (:angle, "angle", "°"))
    return out
end

"""
    step_number(step, key) -> Float64

The value of the number `key` of a placement step ([`step_numbers`](@ref)).
"""
function step_number(st::PlacementStep, key::Symbol)
    haskey(_STEP_COORDS, key) && return Float64(getfield(st, Symbol(_STEP_COORDS[key][1]))[_STEP_COORDS[key][2]])
    v = hasfield(typeof(st), key) ? getfield(st, key) : nothing
    v isa Real || throw(ArgumentError("a $(step_kind(st)) step has no number $key"))
    return Float64(v)
end

# The step with its number `key` set to `v` (the checks of `_state_step`).
function _step_with(st::PlacementStep, key::Symbol, v::Real)
    step_number(st, key)                              # (a number the step has)
    d = _step_state(st)
    haskey(_STEP_COORDS, key) ? (d[_STEP_COORDS[key][1]][_STEP_COORDS[key][2]] = v) : (d[String(key)] = v)
    return _state_step(d)
end

"""
    _euler_angles(R) -> (polar, azimuth, rotation)

Angles (degrees) of a rotation matrix in `euler_rotation`'s order, so that
`euler_rotation(_euler_angles(R)...) ≈ R`: polar in [0, 180], azimuth and rotation in
[0, 360). With the axis along `±z` (polar 0 or 180, to 1e-12 in the sine) the azimuth is
0 and the whole turn is the rotation.
"""
function _euler_angles(R::AbstractMatrix)
    s = hypot(R[1, 3], R[2, 3])
    θ = atand(s, R[3, 3])
    if s < 1e-12
        ψ = R[3, 3] > 0 ? atand(R[2, 1], R[1, 1]) : atand(R[2, 1], R[2, 2])
        return θ, 0.0, mod(ψ, 360)
    end
    return θ, mod(atand(R[2, 3], R[1, 3]), 360), mod(atand(R[3, 2], -R[3, 1]), 360)
end

"""
    nudge_part!(ts::TargetSetup, i, Δ; frame = the part's frame, axis = :z) -> Bool

Move part `i` by `Δ` mm along `axis` (`:x`, `:y`, `:z`, or a vector, normalized) of the
frame named `frame` ([`part_frames`](@ref); by default the frame the part is placed in),
as one edit: the part's center moves by `Δ` times that axis in the chamber, its
orientation and the frame it is placed in stay. Parts placed in its frame, and beams
aimed in it, move with it. A rejected move sets `error` and returns `false`.
"""
function nudge_part!(ts::TargetSetup, i::Integer, Δ::Real; frame::Union{Nothing,Symbol} = nothing,
                     axis = :z)
    _check_part_index(ts, i)
    try
        return _try_build!(ts; parts = nudged_parts(ts.parts[], i, Δ; frame, axis))
    catch err
        err isa ArgumentError || rethrow()
        ts.error[] = _errmsg(err)
        return false
    end
end

"""
    nudged_parts(parts, i, Δ; frame = the part's frame, axis = :z) -> Vector{PartSpec}
    turned_parts(parts, i, θ; frame = the part's frame, axis = :z) -> Vector{PartSpec}

The parts after the edit [`nudge_part!`](@ref) or [`turn_part!`](@ref) makes, without
making it (the placement triad's preview draws them). Throws `ArgumentError` for a move
the edit refuses.
"""
function nudged_parts(parts::AbstractVector{PartSpec}, i::Integer, Δ::Real;
                      frame::Union{Nothing,Symbol} = nothing, axis = :z)
    p = parts[i]
    frame = something(frame, p.frame)
    isfinite(Δ) || throw(ArgumentError("the move must be a finite length"))
    fr = part_frames(parts)
    haskey(fr, frame) || throw(ArgumentError("no frame named $frame"))
    d = Δ * to_chamber(fr[frame], _axis_vector3(axis); direction = true)
    P = fr[p.frame]
    c = part_center(p) + to_frame(P, d; direction = true)
    vals = copy(p.values)
    vals[:x], vals[:y], vals[:z] = c
    ps = copy(parts)
    ps[i] = _with(p; values = vals)
    # placement steps that fix the position along the move leave nothing to move
    if !isempty(p.steps) && Δ != 0 &&
       norm(part_frames(ps)[p.name].origin - fr[p.name].origin) < 1e-12 * abs(Δ)
        throw(ArgumentError("its placement steps fix its position along this axis: edit or sweep a step instead"))
    end
    return ps
end

"""
    turn_part!(ts::TargetSetup, i, θ; frame = the part's frame, axis = :z) -> Bool

Turn part `i` by `θ` degrees (right-handed) about `axis` (`:x`, `:y`, `:z`, or a vector,
normalized) of the frame named `frame` ([`part_frames`](@ref); by default the frame the
part is placed in), through the part origin, as one edit: the angles of its
rotation (`theta`, `phi`, `rotation`) are re-expressed in the frame it is placed in, its
origin and that frame stay. A part without placement steps ends where a `Rotate` step
about the same chamber axis would put it (within round-off in the angles). Parts placed in
its frame, and beams aimed in it, turn with it. The rotational counterpart of
[`nudge_part!`](@ref) (the placement triad's rings make it). A port frame (its axis is
its port's) and a turn its placement steps undo are refused with `error` set, returning
`false`.
"""
function turn_part!(ts::TargetSetup, i::Integer, θ::Real; frame::Union{Nothing,Symbol} = nothing,
                    axis = :z)
    _check_part_index(ts, i)
    try
        return _try_build!(ts; parts = turned_parts(ts.parts[], i, θ; frame, axis))
    catch err
        err isa ArgumentError || rethrow()
        ts.error[] = _errmsg(err)
        return false
    end
end

function turned_parts(parts::AbstractVector{PartSpec}, i::Integer, θ::Real;
                      frame::Union{Nothing,Symbol} = nothing, axis = :z)
    p = parts[i]
    frame = something(frame, p.frame)
    isfinite(θ) || throw(ArgumentError("the turn must be a finite angle"))
    is_port_frame(p) && throw(ArgumentError("a port frame's axis is its port's: choose another port instead"))
    fr = part_frames(parts)
    haskey(fr, frame) || throw(ArgumentError("no frame named $frame"))
    e = to_chamber(fr[frame], _axis_vector3(axis); direction = true)
    Q = Briareus.GUIAPI.axis_rotation(Vec3(e), deg2rad(θ))
    P = fr[p.frame]
    F0 = _own_frame(p, P)                          # before the steps
    O = fr[p.name].origin                          # the pivot is the part origin
    vals = copy(p.values)
    # (without steps the pivot is its center, kept as typed)
    isempty(p.steps) || ((vals[:x], vals[:y], vals[:z]) = to_frame(P, O + Q * (F0.origin - O)))
    vals[:theta], vals[:phi], vals[:rotation] = _euler_angles(P.R' * (Q * F0.R))
    ps = copy(parts)
    ps[i] = _with(p; values = vals)
    if !isempty(p.steps) && θ != 0
        G = part_frames(ps)[p.name]
        norm(G.R - fr[p.name].R) < 1e-12 && norm(G.origin - O) < 1e-12 * max(1.0, norm(O)) &&
            throw(ArgumentError("its placement steps fix its orientation about this axis: edit or sweep a step instead"))
    end
    return ps
end

# A frame axis: `:x`, `:y`, `:z` or a vector (normalized).
function _axis_vector3(axis)
    axis === :x && return Vec3(1, 0, 0)
    axis === :y && return Vec3(0, 1, 0)
    axis === :z && return Vec3(0, 0, 1)
    axis isa Symbol && throw(ArgumentError("the axis is :x, :y, :z or a vector"))
    v = Vec3(axis...)
    norm(v) > 0 || throw(ArgumentError("the axis must not be zero"))
    return normalize(v)
end
