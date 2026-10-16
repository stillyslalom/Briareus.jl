# The 3D viewer: the target mesh colored by component or by a drive field, and the
# beam axes from their ports to their aim points. Under Qt a canvas has a current GL
# context only while Qt renders it, and GLMakie creates or frees GPU objects as soon as
# a plot is added to or removed from a displayed figure. So every plot is created
# before the figure is first shown, and afterwards only its inputs change, atomically
# through `update!` (applied lazily at render time).

# Atomic update of a plot's positional arguments and attributes (positional data go in
# as `arg1…argN`; a lone positional would hit the Dict method).
_update!(plot, args...; attrs...) =
    Makie.update!(plot; (Symbol(:arg, i) => a for (i, a) in enumerate(args))..., attrs...)

const GLTriangleFace = Makie.GeometryBasics.GLTriangleFace
const COMPONENT_COLORS = Makie.wong_colors()
const NAN_COLOR = RGBAf(0.55, 0.55, 0.55, 1.0)
const STATUS_COLORS = Dict(:fail => RGBAf(0.85, 0.15, 0.15, 1), :review => RGBAf(0.95, 0.6, 0.1, 1),
                           :pass => RGBAf(0.2, 0.7, 0.3, 1))
const FIELD_COLORMAP = :inferno
# Components cut away (the quadrant facing the default camera) to show the capsule.
const CUTAWAY_COMPONENTS = (:wall, :cap_top, :cap_bottom)
const _NOPOINT3 = [Point3f(NaN, NaN, NaN)]

"""
    SceneCanvas

The 3D viewer of the Briareus window: `fig`/`ax` (an `Axis3`), the target `body` mesh,
beam axis `segments` and `aims`, the diagnostic-window `outlines`, the `colorbar`,
`cutaway` (cut the quadrant of the hohlraum shell facing the camera; display only),
and a `dirty` flag the shell reads to request a redraw.
"""
struct SceneCanvas
    fig::Figure
    ax::Axis3
    body::Any
    segments::Any
    aims::Any
    outlines::Any
    colorbar::Colorbar
    cutaway::Base.RefValue{Bool}
    dirty::Base.RefValue{Bool}
    ntriangles::Base.RefValue{Int}         # triangles drawn (tests)
end

"""
    scene_canvas(sess::DesignSession) -> SceneCanvas

Build the 3D viewer for `sess` and keep it in sync with the session.
"""
function scene_canvas(sess::DesignSession)
    fig = Figure(size = (900, 800), figure_padding = 6)
    ax = Axis3(fig[1, 1]; aspect = :data, viewmode = :fit, xlabel = "x (mm)",
               ylabel = "y (mm)", zlabel = "z (mm)", protrusions = 30)
    tri = [Point3f(0, 0, 0), Point3f(1, 0, 0), Point3f(0, 1, 0)]          # replaced at once
    body = mesh!(ax, tri, [GLTriangleFace(1, 2, 3)]; color = fill(NAN_COLOR, 3))
    segments = linesegments!(ax, [Point3f(NaN, NaN, NaN), Point3f(NaN, NaN, NaN)];
                             color = [NAN_COLOR, NAN_COLOR], linewidth = 1.5)
    aims = scatter!(ax, _NOPOINT3; color = [NAN_COLOR], markersize = 7)
    outlines = linesegments!(ax, [Point3f(NaN, NaN, NaN), Point3f(NaN, NaN, NaN)];
                             color = [NAN_COLOR, NAN_COLOR], linewidth = 2.5)
    cb = Colorbar(fig[1, 2]; colormap = [NAN_COLOR, NAN_COLOR], limits = (0.0, 1.0), label = " ",
                  height = Relative(0.6))
    c = SceneCanvas(fig, ax, body, segments, aims, outlines, cb, Ref(true), Ref(true), Ref(0))
    dr = sess.drive
    onany((_...) -> _draw_body!(c, sess), sess.target.target, sess.step, dr.result, dr.field)
    onany((_...) -> _draw_beams!(c, sess), sess.target.target, sess.beams.beams,
          sess.beams.enabled, sess.step, sess.review.report, sess.review.view,
          sess.talis.report)
    on(_ -> _draw_windows!(c, sess), sess.target.target)
    _draw_body!(c, sess)
    _draw_beams!(c, sess)
    _draw_windows!(c, sess)
    return c
end

"""
    set_cutaway!(c::SceneCanvas, sess, on::Bool)

Cut away (or restore) the quadrant of the hohlraum shell facing the default view.
"""
function set_cutaway!(c::SceneCanvas, sess::DesignSession, on::Bool)
    c.cutaway[] == on && return c
    c.cutaway[] = on
    _draw_body!(c, sess)
    return c
end

"""
    reset_view!(c::SceneCanvas)

The default viewing angles and limits around the whole scene.
"""
function reset_view!(c::SceneCanvas)
    c.ax.azimuth[] = 1.275π
    c.ax.elevation[] = π / 8
    reset_limits!(c.ax)
    c.dirty[] = true
    return c
end

# What the body is colored by on the current step: a drive field (Drive, Optimize)
# or the components.
function _body_coloring(sess::DesignSession)
    res = sess.drive.result[]
    res === nothing && return nothing
    sess.step[] in (:drive, :optimize) || return nothing
    res.evaluation.scenario.target === sess.target.target[] || return nothing
    return element_values(res, sess.drive.field[])
end

function _draw_body!(c::SceneCanvas, sess::DesignSession)
    t = sess.target.target[]
    vals = _body_coloring(sess)
    cut = c.cutaway[]
    cutcomp = Set(i for (i, comp) in enumerate(t.components) if comp.name in CUTAWAY_COMPONENTS)
    # legend entries: the opaque components (apertures are drawn as holes)
    shown = [i for (i, comp) in enumerate(t.components) if comp.kind === :surface]
    legend = Dict(ci => k for (k, ci) in enumerate(shown))
    pts = Point3f[]
    tvals = Float64[]
    for k in eachindex(t.tri_area)
        t.tri_opaque[k] || continue                 # apertures are holes
        v0, e1, e2 = t.tri_v0[k], t.tri_e1[k], t.tri_e2[k]
        if cut && t.tri_comp[k] in cutcomp
            m = v0 + (e1 + e2) / 3
            (m[1] < 0 && m[2] < 0) && continue
        end
        push!(pts, Point3f(v0...), Point3f((v0 + e1)...), Point3f((v0 + e2)...))
        e = t.tri_front[k]
        push!(tvals, vals === nothing ? legend[t.tri_comp[k]] : e > 0 ? vals[e] : NaN)
    end
    n = length(tvals)
    faces = [GLTriangleFace(3i - 2, 3i - 1, 3i) for i in 1:n]
    cb = c.colorbar
    if vals === nothing
        names = [String(t.components[i].name) for i in shown]
        nc = length(names)
        cmap = [RGBAf(COMPONENT_COLORS[mod1(i, length(COMPONENT_COLORS))]) for i in 1:nc]
        cols = [cmap[Int(v)] for v in tvals]
        # bands of a continuous colormap (a categorical cgrad cannot replace a
        # continuous colormap on a live Colorbar)
        cb.colormap[] = repeat(cmap; inner = 64)
        cb.limits[] = (0.5, nc + 0.5)
        cb.ticks[] = (collect(1.0:nc), names)
        cb.label[] = "component"
    else
        fin = filter(isfinite, tvals)
        lo, hi = isempty(fin) ? (0.0, 1.0) : extrema(fin)
        sess.drive.field[] === :laser && (lo = 0.0)
        hi > lo || (hi = lo + 1)
        cmap = Makie.to_colormap(FIELD_COLORMAP)
        cols = [isfinite(v) ? _cmap_color(cmap, (v - lo) / (hi - lo)) : NAN_COLOR for v in tvals]
        cb.colormap[] = FIELD_COLORMAP
        cb.limits[] = (lo, hi)
        cb.ticks[] = Makie.automatic
        cb.label[] = DRIVE_FIELD_LABELS[sess.drive.field[]]
    end
    vcols = repeat(cols; inner = 3)
    _update!(c.body, pts, faces; color = vcols)
    c.ntriangles[] = n
    c.dirty[] = true
    return c
end

_cmap_color(cmap, f) = RGBAf(cmap[clamp(round(Int, f * (length(cmap) - 1)) + 1, 1, length(cmap))])

function _draw_beams!(c::SceneCanvas, sess::DesignSession)
    t = sess.target.target[]
    bs = sess.beams.beams[]
    en = sess.beams.enabled[]
    if isempty(bs)
        _update!(c.segments, [Point3f(NaN, NaN, NaN), Point3f(NaN, NaN, NaN)];
                 color = [NAN_COLOR, NAN_COLOR])
        _update!(c.aims, _NOPOINT3; color = [NAN_COLOR])
        c.dirty[] = true
        return c
    end
    L = 1.6 * t.scale
    groups = unique(b.group for b in bs)
    status = _beam_status(sess)
    segs = Point3f[]
    aims = Point3f[]
    cols = RGBAf[]
    for (b, on) in zip(bs, en)
        aim = Point3f(b.aim...)
        push!(segs, Point3f((b.aim + L * b.port)...), aim)
        push!(aims, aim)
        col = status !== nothing ? get(STATUS_COLORS, get(status, b.name, :none), NAN_COLOR) :
              RGBAf(COMPONENT_COLORS[mod1(findfirst(==(b.group), groups), length(COMPONENT_COLORS))])
        on || (col = RGBAf(0.5, 0.5, 0.5, 0.25))
        push!(cols, col)
    end
    _update!(c.segments, segs; color = repeat(cols; inner = 2))
    _update!(c.aims, aims; color = cols)
    c.dirty[] = true
    return c
end

# On the Review step, each beam's worst status in the shown report (clearance or
# TaLIS; the TAS stay-in rows are per hemisphere, not per beam).
function _beam_status(sess::DesignSession)
    sess.step[] === :review || return nothing
    r = sess.review.view[] === :talis ? sess.talis.report[] : sess.review.report[]
    r === nothing && return nothing
    rank = Dict(:fail => 0, :review => 1, :pass => 2)
    out = Dict{String,Symbol}()
    for row in r.rows
        row isa TaLISRow && row.rule === :tas_stay_in && continue
        s = get(out, row.beam, :pass)
        out[row.beam] = rank[row.status] < rank[s] ? row.status : s
    end
    return out
end

const WINDOW_COLORS = (closed = RGBAf(0.08, 0.08, 0.08, 1), open = RGBAf(0.9, 0.3, 0.1, 1))

# Outline each diagnostic window of a hohlraum on the wall, just inside and just
# outside it (an open window is a hole; a closed one is a patch colored as its own
# component): arcs along φ at its z edges and axial lines at its φ edges.
function _draw_windows!(c::SceneCanvas, sess::DesignSession)
    ts = sess.target
    ws = ts.kind[] === :hohlraum ? ts.windows[] : WindowSpec[]
    if isempty(ws)
        _update!(c.outlines, [Point3f(NaN, NaN, NaN), Point3f(NaN, NaN, NaN)];
                 color = [NAN_COLOR, NAN_COLOR])
        c.dirty[] = true
        return c
    end
    R = ts.params[][:radius]
    segs = Point3f[]
    cols = RGBAf[]
    narc = 12
    for w in ws, r in (R - 0.01, R + 0.01)
        φa = deg2rad(w.phi) - w.width / 2R
        Δφ = w.width / R
        za, zb = w.z - w.height / 2, w.z + w.height / 2
        p(φ, z) = Point3f(r * cos(φ), r * sin(φ), z)
        col = w.open ? WINDOW_COLORS.open : WINDOW_COLORS.closed
        for z in (za, zb), k in 1:narc
            push!(segs, p(φa + (k - 1) * Δφ / narc, z), p(φa + k * Δφ / narc, z))
            push!(cols, col, col)
        end
        for φ in (φa, φa + Δφ)
            push!(segs, p(φ, za), p(φ, zb))
            push!(cols, col, col)
        end
    end
    _update!(c.outlines, segs; color = cols)
    c.dirty[] = true
    return c
end
