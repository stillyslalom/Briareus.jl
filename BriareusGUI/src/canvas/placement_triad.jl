# The placement triad on the selected part (controllers/placement_triad.jl): three
# arrows and three rings in the axes of the frame the part is placed in, a constant size
# on screen, drawn over the geometry in thin lines. A left drag from a handle previews
# the move (the moved components, the beams aimed in the moved frames and the amount)
# and commits it on release as one edit; a press released within the click slop is a
# click, so the beams and parts behind a handle stay selectable.
#
# Snapping: an arrow moves by whole Move-along steps (the row's "by" field, `triad_step`),
# a ring by whole `TRIAD_TURN_STEP`s; Ctrl held at the press snaps to a tenth of either.

const _TRIAD_ARROW_PX = 70.0       # arrow length on screen (pixels)
const _TRIAD_RING = 0.62           # ring radius / arrow length
const _TRIAD_GRAB_PX = 6.0         # how near (pixels) a press must fall to a handle
const _TRIAD_RING_SEGMENTS = 48
const _TRIAD_ACTIVE = RGBAf(1.0, 0.85, 0.1, 1.0)
const _TRIAD_PREVIEW = RGBAf(0.1, 0.6, 1.0, 0.35)

"""
    TriadDrag

A placement-triad drag in progress: the part (`part`, its index, and `name`), the
`handle` (`TRIAD_HANDLES`), the triad's frame `F` and arrow length `L` (mm) at the press,
the `step` it snaps to (mm or degrees), the pointer's position along the arrow at the
press (`s0`, mm) or its accumulated screen angle about the triad's origin (`angle`,
radians, with `last` the previous pointer angle), the snapped `amount`, the previewed
`parts` (`nothing` when the move is refused, with the reason in `error`).
"""
mutable struct TriadDrag
    part::Int
    name::Symbol
    handle::Tuple{Symbol,Int}
    F::Frame
    L::Float64
    step::Float64
    s0::Float64
    angle::Float64
    last::Float64
    amount::Float64
    parts::Any
    error::String
end

# Length (mm) of a triad arrow at `O`: `_TRIAD_ARROW_PX` pixels at its depth.
function _triad_length(c::CadCamera, O)
    f = camera_frame(c)[1]
    z = c.ortho ? norm(c.lookat - c.eye) : max(dot(Makie.Vec3d(O...) - c.eye, f), 1e-9)
    return _TRIAD_ARROW_PX * 2 * _tanhalf(c) * z / max(c.viewport[4], 1.0)
end

# The polyline (chamber points) of each handle of a triad at frame `F`, arrows `L` long.
function _triad_polyline(F::Frame, L, handle)
    kind, k = handle
    O = Makie.Vec3d(F.origin...)
    e(j) = Makie.Vec3d(F.R[:, mod1(j, 3)]...)
    kind === :move && return [O + 0.18L * e(k), O + L * e(k)]
    ρ = _TRIAD_RING * L
    return [O + ρ * (cos(t) * e(k + 1) + sin(t) * e(k + 2))
            for t in range(0, 2π; length = _TRIAD_RING_SEGMENTS + 1)]
end

# Whether a handle can be dragged as the camera sees it: an arrow not along the line of
# sight (at least 12 pixels long on screen), a ring not edge-on (its axis at most 81°
# from the line of sight).
function _triad_usable(c::CadCamera, F::Frame, L, handle)
    kind, k = handle
    e = Makie.Vec3d(F.R[:, k]...)
    if kind === :move
        a = project_point(c, Makie.Vec3d(F.origin...))
        b = project_point(c, Makie.Vec3d(F.origin...) + L * e)
        return hypot(b[1] - a[1], b[2] - a[2]) >= 12
    end
    return abs(dot(e, camera_frame(c)[1])) >= 0.15
end

"""
    triad_handle_at(c::CadCamera, F, L, px) -> Union{Nothing,Tuple{Symbol,Int}}

The placement-triad handle (`TRIAD_HANDLES`) whose drawn line passes nearest the pixel
`px`, within `_TRIAD_GRAB_PX`, of a triad at frame `F` with arrows `L` mm long; handles
that cannot be dragged from this view (an arrow along the line of sight, a ring
edge-on) are skipped.
"""
function triad_handle_at(c::CadCamera, F::Frame, L, px)
    f = camera_frame(c)[1]
    best, hit = _TRIAD_GRAB_PX, nothing
    for h in TRIAD_HANDLES
        _triad_usable(c, F, L, h) || continue
        pts = _triad_polyline(F, L, h)
        for j in 1:length(pts)-1
            (dot(pts[j] - c.eye, f) > 0 && dot(pts[j+1] - c.eye, f) > 0) || continue
            d = _segment_distance(px, project_point(c, pts[j]), project_point(c, pts[j+1]))
            d < best && ((best, hit) = (d, h))
        end
    end
    return hit
end

# Position (mm) along the axis `e` through `O` nearest the pixel's ray.
function _axis_position(c::CadCamera, O, e, px)
    o, d = pixel_ray(c, px)
    w = Makie.Vec3d(O...) - o
    b = dot(e, d)
    den = 1 - b^2
    den < 1e-9 && return 0.0
    return (b * dot(d, w) - dot(e, w)) / den
end

# Screen angle (radians) of the pixel about the projected origin of `F`.
function _screen_angle(c::CadCamera, F::Frame, px)
    o = project_point(c, Makie.Vec3d(F.origin...))
    return atan(px[2] - o[2], px[1] - o[1])
end

"""
    triad_amount(c::CadCamera, drag::TriadDrag, px) -> Float64

The unsnapped amount of a placement-triad drag with the pointer at `px`: for an arrow,
the move (mm) along its axis from the press, where the pointer's ray passes nearest the
axis; for a ring, the turn (degrees, right-handed about its axis) the pointer has swept
about the triad's origin on screen since the press (counterclockwise is positive with
the axis toward the viewer). Updates the ring's accumulated angle.
"""
function triad_amount(c::CadCamera, g::TriadDrag, px)
    kind, k = g.handle
    e = Makie.Vec3d(g.F.R[:, k]...)
    kind === :move && return _axis_position(c, g.F.origin, e, px) - g.s0
    a = _screen_angle(c, g.F, px)
    g.angle += rem2pi(a - g.last, RoundNearest)
    g.last = a
    s = dot(e, camera_frame(c)[1]) < 0 ? 1.0 : -1.0      # the axis toward the viewer: counterclockwise
    return rad2deg(s * g.angle)
end

# Press: grab a handle of the shown triad (`c.triad_at`).
function _triad_press!(c, sess, px, ctrl::Bool)
    at = c.triad_at[]
    at === nothing && return false
    i, name, F = at
    L = _triad_length(c.camera, F.origin)
    h = triad_handle_at(c.camera, F, L, px)
    h === nothing && return false
    base = h[1] === :move ? c.triad_step[] : TRIAD_TURN_STEP
    step = ctrl ? base / TRIAD_FINE : base
    s0 = h[1] === :move ? _axis_position(c.camera, F.origin, Makie.Vec3d(F.R[:, h[2]]...), px) : 0.0
    a0 = _screen_angle(c.camera, F, px)
    c.placing[] = TriadDrag(i, name, h, F, L, step, s0, 0.0, a0, 0.0, nothing, "")
    return true
end

# Drag: the snapped amount, and its preview when it changes.
function _triad_move!(c, sess, px)
    g = c.placing[]
    g === nothing && return
    amount = triad_snap(triad_amount(c.camera, g, px), g.step)
    amount == g.amount && return
    g.amount = amount
    ts = sess.target
    i = part_index(ts, g.name)
    if i === nothing || amount == 0
        g.parts, g.error = nothing, ""
    else
        try
            g.parts, g.error = triad_parts(ts.parts[], i, g.handle, amount), ""
        catch err
            err isa ArgumentError || rethrow()
            g.parts, g.error = nothing, Controllers._errmsg(err)
        end
    end
    _draw_placement_triad!(c, sess)
    return
end

# Release: one edit (none for a zero amount), then the triad where the part now is.
function _triad_release!(c, sess, px)
    _triad_move!(c, sess, px)
    g = c.placing[]
    c.placing[] = nothing
    g === nothing && return
    ts = sess.target
    i = part_index(ts, g.name)
    if i !== nothing && g.amount != 0
        ok = commit_triad!(ts, i, g.handle, g.amount)
        kind, k = g.handle
        frame = ts.parts[][i].frame
        sess.status[] = !ok ? "not moved: " * ts.error[] :
                        kind === :move ? @sprintf("moved %s by %+g mm along %s of %s", g.name, g.amount,
                                                  ("x", "y", "z")[k], frame) :
                        @sprintf("turned %s by %+g° about %s of %s", g.name, g.amount, ("x", "y", "z")[k], frame)
    end
    _place_triad!(c, sess)
    return
end

_triad_cancel!(c, sess) = (c.placing[] = nothing; _draw_placement_triad!(c, sess); nothing)

# Where the triad is: the selection's part and its triad frame, or nothing.
function _place_triad!(c, sess)
    i = triad_subject(sess)
    c.triad_at[] = i === nothing ? nothing :
                   (try
                        (i, sess.target.parts[][i].name, triad_frame(sess.target, i))
                    catch err
                        err isa ArgumentError || rethrow()
                        nothing
                    end)
    _draw_placement_triad!(c, sess)
    return c
end

# The triad (at its previewed place during a drag, the dragged handle highlighted), the
# preview's components and beams, and the amount.
function _draw_placement_triad!(c, sess)
    at, g = c.triad_at[], c.placing[]
    if at === nothing
        _update!(c.place_arrows, _NOSEG3; color = [NAN_COLOR, NAN_COLOR], visible = false)
        _update!(c.place_rings, _NOSEG3; color = [NAN_COLOR, NAN_COLOR], visible = false)
        _update!(c.place_ghost, _NOTRI...; visible = false)
        _update!(c.place_beams, _NOSEG3; visible = false)
        _update!(c.place_label, _NOPOINT3; text = [""], visible = false)
        c.dirty[] = true
        return c
    end
    F = at[3]
    if g !== nothing && g.parts !== nothing
        F = try
            triad_frame(g.parts, g.part)
        catch err
            err isa ArgumentError || rethrow()
            F
        end
    end
    L = _triad_length(c.camera, F.origin)
    arrows, acols, rings, rcols = Point3f[], RGBAf[], Point3f[], RGBAf[]
    for h in TRIAD_HANDLES
        col = g !== nothing && g.handle == h ? _TRIAD_ACTIVE : TRIAD_COLORS[h[2]]
        pts = _triad_polyline(F, L, h)
        if h[1] === :move
            a, b = pts[1], pts[2]
            # an arrowhead in the plane facing the camera
            e = (b - a) / norm(b - a)
            side = cross(e, camera_frame(c.camera)[1])
            side = norm(side) < 1e-9 ? Makie.Vec3d(0, 0, 0) : side / norm(side)
            base = b - 0.16L * e
            for q in (a, b, b, base + 0.07L * side, b, base - 0.07L * side)
                push!(arrows, Point3f(q...))
                push!(acols, col)
            end
        else
            for j in 1:length(pts)-1
                push!(rings, Point3f(pts[j]...), Point3f(pts[j+1]...))
                push!(rcols, col, col)
            end
        end
    end
    _update!(c.place_arrows, arrows; color = acols, visible = true)
    _update!(c.place_rings, rings; color = rcols, visible = true)
    if g === nothing
        _update!(c.place_ghost, _NOTRI...; visible = false)
        _update!(c.place_beams, _NOSEG3; visible = false)
        _update!(c.place_label, _NOPOINT3; text = [""], visible = false)
    else
        ghost, segs = Point3f[], Point3f[]
        if g.parts !== nothing
            pv = triad_preview(sess, g.parts)
            for comp in pv.components, (i, j, k) in comp.mesh.triangles
                vs = comp.mesh.vertices
                push!(ghost, Point3f(vs[i]...), Point3f(vs[j]...), Point3f(vs[k]...))
            end
            paths = c.paths[]
            for (b, aim) in pv.aims
                b <= length(paths) || continue
                p = paths[b]
                push!(segs, Point3f((p.start + (aim - p.aim))...), Point3f(aim...))
            end
        end
        isempty(ghost) ? _update!(c.place_ghost, _NOTRI...; visible = false) :
                         _update!(c.place_ghost, ghost, _faces(length(ghost) ÷ 3); visible = true)
        _update!(c.place_beams, isempty(segs) ? _NOSEG3 : segs; visible = !isempty(segs))
        unit = g.handle[1] === :move ? " mm" : "°"
        text = isempty(g.error) ? @sprintf("%+g%s", g.amount, unit) : g.error
        _update!(c.place_label, [Point3f((Makie.Vec3d(F.origin...) + 1.05L * Makie.Vec3d(F.R[:, g.handle[2]]...))...)];
                 text = [text], visible = true)
    end
    c.dirty[] = true
    return c
end
