# The placement triad of the 3D view (assemblies design, decision 6): arrows that move
# the selected part along, and rings that turn it about, the axes of the frame it is
# placed in. A drag previews its move and, on release, makes one edit through
# `nudge_part!` or `turn_part!` (the Move-along row's commit), so it is undone as one and
# the placement stays parametric. The canvas (canvas/placement_triad.jl) maps the mouse
# to a handle and an amount; this file holds what is framework-free.

"""
The placement triad's handles, `(kind, k)`: `(:move, k)` the arrow along axis `k` (1 x,
2 y, 3 z) of the frame the part is placed in, `(:turn, k)` the ring about it.
"""
const TRIAD_HANDLES = ((:move, 1), (:move, 2), (:move, 3), (:turn, 1), (:turn, 2), (:turn, 3))

"""
Snap increment (degrees) of a placement-triad ring; an arrow snaps to the Move-along
row's step. Ctrl held at the press divides either by [`TRIAD_FINE`](@ref). An
engineering default (15° is the arrow keys' orbit step).
"""
const TRIAD_TURN_STEP = 15.0

"Divisor of the placement triad's snap increments with Ctrl held at the press."
const TRIAD_FINE = 10

"""
    triad_subject(sess) -> Union{Nothing,Int}

The part (its index in `sess.target.parts[]`) the placement triad moves: with one node
selected, the part, assembly or user frame it names, or the part a piece, a window or
a part's component belongs to (pieces are not separately movable). `nothing` with
nothing or several nodes selected, a beam or an analysis selected, or a port frame (its
axis is its port's).
"""
function triad_subject(sess)
    ids = selected_ids(sess.view)
    length(ids) == 1 || return nothing
    ts = sess.target
    kind, key = selection_kind(only(ids))
    name = Symbol(key)
    i = if kind in (:part, :assembly, :frame)
        part_index(ts, name)
    elseif kind === :window
        w = window_of(ts, name)
        w === nothing ? nothing : w[1]
    elseif kind === :comp
        _component_part(ts, name)
    else
        nothing
    end
    (i === nothing || is_port_frame(ts.parts[][i])) && return nothing
    return i
end

# The part a target component (drawn, or a ghost left out of the calculations) belongs to.
function _component_part(ts, name::Symbol)
    for t in (ts.target[], ts.ghost[])
        t === nothing && continue
        owner = get(_component_owners(t, ts.parts[]), name, nothing)
        owner === nothing || return part_index(ts, owner)
    end
    return nothing
end

"""
    triad_frame(parts, i) -> Frame
    triad_frame(ts::TargetSetup, i) -> Frame

Where the placement triad of part `i` is drawn: at the part's own origin, with the axes
of the frame it is placed in ("Placed in"), along which its arrows move the part and
about which its rings turn it.
"""
function triad_frame(parts::AbstractVector{PartSpec}, i::Integer)
    fr = part_frames(parts)
    p = parts[i]
    return Frame(fr[p.frame].R, fr[p.name].origin)
end
triad_frame(ts::TargetSetup, i::Integer) = triad_frame(ts.parts[], i)

"""
    triad_parts(parts, i, handle, amount) -> Vector{PartSpec}

The parts after a placement-triad drag of `amount` on `handle` ([`TRIAD_HANDLES`](@ref))
of part `i`: mm along an arrow ([`nudged_parts`](@ref)), degrees about a ring
([`turned_parts`](@ref)), in the axes of the frame the part is placed in. Throws
`ArgumentError` for a move the edit refuses.
"""
function triad_parts(parts::AbstractVector{PartSpec}, i::Integer, handle::Tuple{Symbol,Int}, amount::Real)
    kind, k = handle
    axis = (:x, :y, :z)[k]
    frame = parts[i].frame
    kind === :move && return nudged_parts(parts, i, amount; frame, axis)
    kind === :turn && return turned_parts(parts, i, amount; frame, axis)
    throw(ArgumentError("a triad handle is :move or :turn"))
end

"""
    commit_triad!(ts::TargetSetup, i, handle, amount) -> Bool

Make a placement-triad drag as one edit: [`nudge_part!`](@ref) by `amount` mm along an
arrow's axis, or [`turn_part!`](@ref) by `amount` degrees about a ring's, in the frame
part `i` is placed in (the edit the Move-along row's buttons make). A zero amount makes
no edit and returns `true`; a refused one sets `error` and returns `false`.
"""
function commit_triad!(ts::TargetSetup, i::Integer, handle::Tuple{Symbol,Int}, amount::Real)
    amount == 0 && return true
    kind, k = handle
    axis = (:x, :y, :z)[k]
    frame = ts.parts[][i].frame
    kind === :move && return nudge_part!(ts, i, amount; frame, axis)
    kind === :turn && return turn_part!(ts, i, amount; frame, axis)
    throw(ArgumentError("a triad handle is :move or :turn"))
end

"""
    triad_snap(x, step) -> Float64

`x` rounded to a whole number of `step`s, as the decimal it prints as (12 significant
digits): three steps of 0.1 mm are 0.3 mm, the move the Move-along row makes with a step
of 0.3.
"""
function triad_snap(x::Real, step::Real)
    step > 0 || throw(ArgumentError("the snap step must be positive"))
    n = round(x / step)
    return n == 0 ? 0.0 : round(n * step; sigdigits = 12)
end

"""
    triad_preview(sess, parts) -> (; components, aims)

The design as `parts` ([`triad_parts`](@ref)) would place it, without rebuilding it:
the components of each part whose frame changes, moved rigidly by that frame's map
(drawn and excluded components alike), and `aims`, `beam index => aim` for the beams
aimed in a moved frame (the aims a rebuilt design gives them, to round-off).
"""
function triad_preview(sess, parts::AbstractVector{PartSpec})
    ts, bt = sess.target, sess.beams
    old, new = part_frames(ts), part_frames(parts)
    maps = Dict{Symbol,Frame}(n => compose(F, inv(old[n])) for (n, F) in new if haskey(old, n) && F != old[n])
    comps = Component[]
    for t in (ts.target[], ts.ghost[])
        t === nothing && continue
        owner = _component_owners(t, ts.parts[])
        for c in t.components
            M = get(maps, get(owner, c.name, :none), nothing)
            M === nothing || push!(comps, transform(c, M.R, M.origin))
        end
    end
    aims = Pair{Int,Vec3}[i => Vec3(to_chamber(maps[f], b.aim))
                          for (i, (b, f)) in enumerate(zip(bt.beams[], bt.aim_frame[])) if haskey(maps, f)]
    return (; components = comps, aims)
end
