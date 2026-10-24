# The 3D viewer: the target mesh colored by component or by a drive field, the beams
# as rays (axis from a start point to where it lands) or volumes (envelope cut where
# it meets the target), the selection highlighted, an orientation triad, and the CAD
# camera (cad_camera.jl) with picking by ray casting into the target.
#
# Under Qt a canvas has a current GL context only while Qt renders it, and GLMakie
# creates or frees GPU objects as soon as a plot is added to or removed from a
# displayed figure. So every plot is created before the figure is first shown, and
# afterwards only its inputs change, atomically through `update!` (applied lazily at
# render time).

# Atomic update of a plot's positional arguments and attributes (positional data go in
# as `arg1…argN`; a lone positional would hit the Dict method).
_update!(plot, args...; attrs...) =
    Makie.update!(plot; (Symbol(:arg, i) => a for (i, a) in enumerate(args))..., attrs...)

const GLTriangleFace = Makie.GeometryBasics.GLTriangleFace
_hexcolor(h::AbstractString) = RGBAf(Makie.to_color(h))
const NAN_COLOR = RGBAf(0.55, 0.55, 0.55, 1.0)
const STATUS_COLORS = Dict(:fail => RGBAf(0.85, 0.15, 0.15, 1), :review => RGBAf(0.95, 0.6, 0.1, 1),
                           :pass => RGBAf(0.2, 0.7, 0.3, 1))
const FIELD_COLORMAP = :inferno
const HIGHLIGHT = RGBAf(0.1, 0.6, 1.0, 0.8)
const GHOST = RGBAf(0.6, 0.6, 0.65, 0.18)       # components left out of the calculations
const BEAM_SELECTED = RGBAf(1.0, 0.95, 0.2, 1.0)
const BEAM_OFF = RGBAf(0.5, 0.5, 0.5, 0.25)
const VOLUME_ALPHA = 0.22f0
const _NOPOINT3 = [Point3f(NaN, NaN, NaN)]
const _NOSEG3 = [Point3f(NaN, NaN, NaN), Point3f(NaN, NaN, NaN)]
const _NOTRI = ([Point3f(0, 0, 0), Point3f(0, 0, 0), Point3f(0, 0, 0)], [GLTriangleFace(1, 2, 3)])
const _PICK_PIXELS = 7.0           # how near (pixels) a click must fall to a beam's axis

"""
    SceneCanvas

The 3D viewer of the Briareus window: `fig`, the `lscene`, the `camera` (a
[`CadCamera`](@ref)), plots (`body`, `highlight`, `rays`, `selected_rays`, `aims`,
`hits`, `volumes`, `outlines`), the `colorbar`, the triad scene, the beam `paths` last
drawn, which triangles are drawn (`drawn`, for picking), and a `dirty` flag the shell
reads to request a redraw.
"""
struct SceneCanvas
    fig::Figure
    lscene::LScene
    camera::CadCamera
    body::Any
    highlight::Any
    rays::Any
    selected_rays::Any
    aims::Any
    hits::Any
    volumes::Any
    outlines::Any
    colorbar::Colorbar
    triad::Scene
    paths::Base.RefValue{Vector{BeamPath}}
    drawn::Base.RefValue{BitVector}
    fitted::Base.RefValue{Any}                  # the target the camera was last fitted to
    dirty::Base.RefValue{Bool}
    ntriangles::Base.RefValue{Int}               # triangles drawn (tests)
    ghost::Any                                   # components left out of the calculations
end

"""
    scene_canvas(sess::DesignSession) -> SceneCanvas

Build the 3D viewer for `sess` and keep it in sync with the session.
"""
function scene_canvas(sess::DesignSession)
    fig = Figure(size = (900, 800), figure_padding = 0)
    ls = LScene(fig[1, 1]; show_axis = false, scenekw = (camera = s -> nothing,))
    sc = ls.scene
    tri = _NOTRI
    body = mesh!(sc, tri...; color = fill(NAN_COLOR, 3))
    highlight = mesh!(sc, tri...; color = HIGHLIGHT, transparency = true, depth_shift = -1f-4,
                      visible = false)
    volumes = mesh!(sc, tri...; color = fill(RGBAf(0, 0, 0, 0), 3), transparency = true,
                    visible = false)
    ghost = mesh!(sc, tri...; color = GHOST, transparency = true, visible = false)
    rays = linesegments!(sc, _NOSEG3; color = [NAN_COLOR, NAN_COLOR], linewidth = 1.5)
    selected_rays = linesegments!(sc, _NOSEG3; color = BEAM_SELECTED, linewidth = 4)
    aims = scatter!(sc, _NOPOINT3; color = [NAN_COLOR], markersize = 6, strokewidth = 0.5,
                    strokecolor = :white)
    hits = scatter!(sc, _NOPOINT3; color = [NAN_COLOR], markersize = 8, marker = :diamond)
    outlines = linesegments!(sc, _NOSEG3; color = [NAN_COLOR, NAN_COLOR], linewidth = 2.5)
    cb = Colorbar(fig[1, 2]; colormap = [NAN_COLOR, NAN_COLOR], limits = (0.0, 1.0), label = " ",
                  height = Relative(0.6), tellheight = false)
    triad = _triad_scene(sc)
    cam = CadCamera()
    c = SceneCanvas(fig, ls, cam, body, highlight, rays, selected_rays, aims, hits, volumes,
                    outlines, cb, triad, Ref(BeamPath[]), Ref(falses(0)), Ref{Any}(nothing),
                    Ref(true), Ref(0), ghost)
    cam.pick = (o, d) -> min(_pick_surface(c, sess, o, d)[1], _pick_ghost(sess, o, d)[1])
    cam.pivot = () -> selection_center(sess)
    cam.on_click = (o, d, px) -> _click!(c, sess, o, d, px)
    cam.on_double = (o, d, px) -> _double_click!(c, sess, o, d, px)
    cam.on_change = () -> (_update_triad!(c); c.dirty[] = true)
    attach_cad_camera!(sc, cam)
    dr, v = sess.drive, sess.view
    onany((_...) -> _draw_body!(c, sess), sess.target.target, dr.result, v.color, dr.field,
          sess.time.cursor, sess.time.history, v.cutaway, v.hidden, v.cut_azimuth, v.cut_width, v.cut_components)
    onany((_...) -> _draw_beams!(c, sess), sess.target.target, sess.beams.beams,
          sess.beams.enabled, sess.beams.cluster, v.beams, v.beam_color, v.beam_length, sess.review.report,
          sess.talis.report, v.hidden_beams, v.show_off)
    onany((_...) -> _draw_selection!(c, sess), v.selection, v.extra, sess.target.target,
          v.cutaway, v.cut_azimuth, v.cut_width, v.cut_components, sess.beams.beams,
          sess.beams.cluster)
    on(_ -> _draw_windows!(c, sess), sess.target.target)
    onany((_...) -> _draw_ghost!(c, sess), sess.target.ghost, v.hidden)
    _draw_ghost!(c, sess)
    _draw_body!(c, sess)
    _draw_beams!(c, sess)
    _draw_windows!(c, sess)
    _draw_selection!(c, sess)
    return c
end

"""
    reset_view!(c::SceneCanvas) / fit_all!(c) / standard_view!(c, name)

The default isometric view fitted to the target; fit the target from the current
direction; look along a standard direction (`STANDARD_VIEWS`).
"""
reset_view!(c::SceneCanvas) = (orient_camera!(c.camera, :iso); c)
fit_all!(c::SceneCanvas) = (fit_camera!(c.camera); c)
standard_view!(c::SceneCanvas, name::Symbol) = (orient_camera!(c.camera, name); c)

"""
    fit_selection!(c::SceneCanvas, sess)

Fit the view to the selected component or beams (the whole target without a selection).
"""
function fit_selection!(c::SceneCanvas, sess::DesignSession)
    b = selection_bounds(sess)
    b === nothing ? fit_camera!(c.camera) : fit_camera!(c.camera, b...)
    return c
end

"""
    look_along_beam!(c::SceneCanvas, sess)

Look down the first selected beam from its port toward its aim (the view of the
target the beam's lens has), fitted to the target.
"""
function look_along_beam!(c::SceneCanvas, sess::DesignSession)
    idx = selected_beams(sess)
    isempty(idx) && return c
    b = sess.beams.beams[][idx[1]]
    k = direction(b)
    up = abs(k[3]) < 0.95 ? Makie.Vec3d(0, 0, 1) : Makie.Vec3d(1, 0, 0)
    cam = c.camera
    cam.lookat = Makie.Vec3d(b.aim...)
    orient_camera!(cam, Makie.Vec3d(k...), up)
    fit_camera!(cam, Makie.Vec3d(b.aim...), sess.target.target[].scale)
    return c
end

"""
    look_from_port!(c::SceneCanvas, sess, theta, phi)

Look at the target center from the port at polar angle `theta` and azimuth `phi`
(degrees; a beam or diagnostic port), fitted to the target.
"""
function look_from_port!(c::SceneCanvas, sess::DesignSession, θ::Real, φ::Real)
    k = -port_vector(θ, φ)
    up = abs(k[3]) < 0.95 ? Makie.Vec3d(0, 0, 1) : Makie.Vec3d(1, 0, 0)
    cam = c.camera
    cam.lookat = Makie.Vec3d(0, 0, 0)
    orient_camera!(cam, Makie.Vec3d(k...), up)
    fit_camera!(cam, Makie.Vec3d(0, 0, 0), sess.target.target[].scale)
    return c
end

"""
    set_projection!(c::SceneCanvas, ortho::Bool) / set_turntable!(c, on::Bool)
"""
set_projection!(c::SceneCanvas, ortho::Bool) = (c.camera.ortho = ortho; c.camera.on_change(); c)
set_turntable!(c::SceneCanvas, on::Bool) = (c.camera.turntable = on; c)
release_input!(c::SceneCanvas) = (release_input!(c.lscene.scene, c.camera); c)

# ---------------------------------------------------------------- picking

# Triangles drawn (visible, not cut away): picking passes through the others.
function _drawn_triangles(t::Target, v::ViewState)
    hid = v.hidden[]
    out = falses(length(t.tri_area))
    for k in eachindex(t.tri_area)
        t.tri_opaque[k] || continue
        name = t.components[t.tri_comp[k]].name
        name in hid && continue
        # whole patches (the element's centroid), so the cut edge follows the grid
        e = t.tri_front[k]
        m = e > 0 ? t.elem_centroid[e] : t.tri_v0[k] + (t.tri_e1[k] + t.tri_e2[k]) / 3
        is_cut(v, name, m) && continue
        out[k] = true
    end
    return out
end

# First drawn surface along a ray: (distance, triangle), (Inf, 0) when none.
function _pick_surface(c::SceneCanvas, sess::DesignSession, o, d)
    t = sess.target.target[]
    drawn = c.drawn[]
    length(drawn) == length(t.tri_area) || return (Inf, 0)
    o, d = Briareus.Vec3(o...), Briareus.Vec3(d...)
    traveled = 0.0
    for _ in 1:64
        τ, tri = Briareus.closest_hit(t.bvh, o, d)
        tri == 0 && return (Inf, 0)
        traveled += τ
        drawn[tri] && return (traveled, Int(tri))
        step = 1e-6 * max(t.scale, 1.0)
        o = o + (τ + step) * d
        traveled += step
    end
    return (Inf, 0)
end

# First shown ghost (a component left out of the calculations) along a ray: (distance,
# component name), (Inf, :none) when none. Ghosts are drawn whole (no cutaway).
function _pick_ghost(sess::DesignSession, o, d)
    g = sess.target.ghost[]
    g === nothing && return (Inf, :none)
    hid = sess.view.hidden[]
    o, d = Briareus.Vec3(o...), Briareus.Vec3(d...)
    traveled = 0.0
    for _ in 1:64
        τ, tri = Briareus.closest_hit(g.bvh, o, d)
        tri == 0 && return (Inf, :none)
        traveled += τ
        name = g.components[g.tri_comp[tri]].name
        name in hid || return (traveled, name)
        step = 1e-6 * max(g.scale, 1.0)
        o = o + (τ + step) * d
        traveled += step
    end
    return (Inf, :none)
end

# The beam whose drawn axis passes nearest the pixel, within `_PICK_PIXELS`.
function _pick_beam(c::SceneCanvas, sess::DesignSession, px)
    sess.view.beams[] === :off && return nothing
    best, bi = _PICK_PIXELS, nothing
    hid = sess.view.hidden_beams[]
    bs = sess.beams.beams[]
    en, off = sess.beams.enabled[], sess.view.show_off[]
    for (i, p) in enumerate(c.paths[])
        i <= length(bs) && (bs[i].name in hid || !(off || en[i])) && continue
        a, b = project_point(c.camera, Makie.Vec3d(p.start...)), project_point(c.camera, Makie.Vec3d(p.stop...))
        d = _segment_distance(px, a, b)
        d < best && ((best, bi) = (d, i))
    end
    return bi
end

function _segment_distance(p, a, b)
    ab = (b[1] - a[1], b[2] - a[2])
    L2 = ab[1]^2 + ab[2]^2
    s = L2 > 0 ? clamp(((p[1] - a[1]) * ab[1] + (p[2] - a[2]) * ab[2]) / L2, 0.0, 1.0) : 0.0
    return hypot(p[1] - a[1] - s * ab[1], p[2] - a[2] - s * ab[2])
end

function _click!(c::SceneCanvas, sess::DesignSession, o, d, px)
    bi = _pick_beam(c, sess, px)
    τ, tri = _pick_surface(c, sess, o, d)
    τg, ghost = _pick_ghost(sess, o, d)
    add = c.camera.click_ctrl
    if bi !== nothing
        b = sess.beams.beams[][bi]
        add ? select_add!(sess.view, beam_id(b)) : select!(sess.view, beam_id(b))
    elseif τg < τ                       # a part left out of the calculations, in front
        id = component_id(sess.target, ghost)
        add ? select_add!(sess.view, id) : select!(sess.view, id)
    elseif tri != 0
        t = sess.target.target[]
        id = component_id(sess.target, t.components[t.tri_comp[tri]].name)
        add ? select_add!(sess.view, id) : select!(sess.view, id)
    else
        add || select!(sess.view, "")
    end
    c.dirty[] = true
    return nothing
end

# Double click: fit the clicked component (or beam), or the whole target on empty space.
function _double_click!(c::SceneCanvas, sess::DesignSession, o, d, px)
    _click!(c, sess, o, d, px)
    fit_selection!(c, sess)
    return nothing
end

"""
    selection_bounds(sess) -> Union{Nothing,Tuple{Vec3d,Float64}}

Bounding sphere of the selected component or beams (their aim points and landing
points), `nothing` without a component or beam selection.
"""
function selection_bounds(sess::DesignSession)
    t = sess.target.target[]
    names = with_lips(selected_components(sess))      # an entrance hole with its lip
    pts = Briareus.Vec3[]
    if !isempty(names)
        # the selected components, in the calculations or excluded (drawn as ghosts)
        for tt in (t, sess.target.ghost[])
            tt === nothing && continue
            for k in eachindex(tt.tri_area)
                tt.components[tt.tri_comp[k]].name in names || continue
                push!(pts, tt.tri_v0[k], tt.tri_v0[k] + tt.tri_e1[k], tt.tri_v0[k] + tt.tri_e2[k])
            end
        end
    else
        idx = selected_beams(sess)
        isempty(idx) && return nothing
        for i in idx
            p = beam_path(t, sess.beams.beams[][i])
            push!(pts, p.aim, p.stop)
        end
    end
    isempty(pts) && return nothing
    lo = reduce((a, b) -> min.(a, b), pts)
    hi = reduce((a, b) -> max.(a, b), pts)
    ctr = (lo + hi) / 2
    return Makie.Vec3d(ctr...), max(norm(hi - lo) / 2, 0.05 * t.scale)
end

"Center of the selection (the rotation pivot), or `nothing`."
selection_center(sess::DesignSession) = (b = selection_bounds(sess); b === nothing ? nothing : b[1])

# ---------------------------------------------------------------- target

# What the body is colored by: a drive field when chosen and computed for this target,
# at the time cursor's step of the drive history when it shows one (`cursor_values`).
function _body_coloring(sess::DesignSession)
    field = sess.view.color[]
    field === :component && return nothing
    t = sess.target.target[]
    dr = sess.drive
    vals = cursor_values(sess, field)
    if vals === nothing
        res = dr.result[]
        res === nothing && return nothing
        s = res.evaluation.scenario
        s.target === t && return element_values(res, field)
        s === dr.scenario[] || return nothing
        vals = element_values(res, field)
    else
        dr.scenario[].target === t && return vals
    end
    # a drive on the cavity: its elements are some of the target's
    idx = dr.cavity_index[]
    (idx !== nothing && dr.key[] !== nothing && dr.key[][1] === t) || return nothing
    full = fill(NaN, nelements(t))
    full[idx] = vals
    return full
end

# Components left out of the calculations, drawn faint (not pickable: rays and the
# picker see only the calculation target).
function _draw_ghost!(c::SceneCanvas, sess::DesignSession)
    g = sess.target.ghost[]
    hid = sess.view.hidden[]
    pts = Point3f[]
    if g !== nothing
        for k in eachindex(g.tri_area)
            g.components[g.tri_comp[k]].name in hid && continue
            v0, e1, e2 = g.tri_v0[k], g.tri_e1[k], g.tri_e2[k]
            push!(pts, Point3f(v0...), Point3f((v0 + e1)...), Point3f((v0 + e2)...))
        end
    end
    if isempty(pts)
        _update!(c.ghost, _NOTRI...; visible = false)
    else
        _update!(c.ghost, pts, [GLTriangleFace(3i - 2, 3i - 1, 3i) for i in 1:length(pts)÷3];
                 visible = true)
    end
    c.dirty[] = true
    return c
end

function _draw_body!(c::SceneCanvas, sess::DesignSession)
    t = sess.target.target[]
    g = sess.target.ghost[]
    # with the beams drawn, and the parts left out of the calculations (CAD mounting
    # hardware can reach 10² mm): the clip planes hold them all
    c.camera.bounds = (Makie.Vec3d(0, 0, 0), max(2.2 * t.scale, g === nothing ? 0.0 : 1.05 * g.scale))
    c.camera.home = (Makie.Vec3d(0, 0, 0), t.scale)
    if c.fitted[] === nothing
        c.camera.lookat = Makie.Vec3d(0, 0, 0)
        orient_camera!(c.camera, :iso)
        c.fitted[] = t
    end
    vals = _body_coloring(sess)
    drawn = _drawn_triangles(t, sess.view)
    c.drawn[] = drawn
    shown = [i for (i, comp) in enumerate(t.components) if comp.kind === :surface]
    legend = Dict(ci => k for (k, ci) in enumerate(shown))
    pts = Point3f[]
    tvals = Float64[]
    for k in eachindex(t.tri_area)
        drawn[k] || continue
        v0, e1, e2 = t.tri_v0[k], t.tri_e1[k], t.tri_e2[k]
        push!(pts, Point3f(v0...), Point3f((v0 + e1)...), Point3f((v0 + e2)...))
        e = t.tri_front[k]
        push!(tvals, vals === nothing ? legend[t.tri_comp[k]] : e > 0 ? vals[e] : NaN)
    end
    n = length(tvals)
    cb = c.colorbar
    if n == 0
        _update!(c.body, _NOTRI...; color = fill(RGBAf(0, 0, 0, 0), 3))
    elseif vals === nothing
        # component colors are keyed in the model tree (its swatches), not on the canvas
        pal = component_colors(t)
        cmap = [_hexcolor(pal[t.components[i].name]) for i in shown]
        cols = [cmap[Int(v)] for v in tvals]
        _show_colorbar!(c, false)
        _update!(c.body, pts, [GLTriangleFace(3i - 2, 3i - 1, 3i) for i in 1:n];
                 color = repeat(cols; inner = 3))
    else
        fin = filter(isfinite, tvals)
        lo, hi = isempty(fin) ? (0.0, 1.0) : extrema(fin)
        sess.view.color[] === :laser && (lo = 0.0)
        hi > lo || (hi = lo + 1)
        # a drive history's step: one scale for the whole history
        rng = cursor_range(sess, sess.view.color[])
        rng === nothing || ((lo, hi) = rng)
        at = time_label(sess)
        cmap = Makie.to_colormap(FIELD_COLORMAP)
        cols = [isfinite(v) ? _cmap_color(cmap, (v - lo) / (hi - lo)) : NAN_COLOR for v in tvals]
        cb.colormap[] = FIELD_COLORMAP
        cb.limits[] = (lo, hi)
        cb.ticks[] = Makie.automatic
        cb.label[] = DRIVE_FIELD_LABELS[sess.view.color[]] * (isempty(at) ? "" : " at " * at[5:end])
        _show_colorbar!(c, true)
        _update!(c.body, pts, [GLTriangleFace(3i - 2, 3i - 1, 3i) for i in 1:n];
                 color = repeat(cols; inner = 3))
    end
    c.ntriangles[] = n
    c.dirty[] = true
    return c
end

_cmap_color(cmap, f) = RGBAf(cmap[clamp(round(Int, f * (length(cmap) - 1)) + 1, 1, length(cmap))])

# The colorbar takes room only for a drive field: hidden, it has no width, ticks or label.
function _show_colorbar!(c::SceneCanvas, show::Bool)
    cb = c.colorbar
    cb.blockscene.visible[] == show && return
    cb.blockscene.visible[] = show
    for a in (cb.ticklabelsvisible, cb.ticksvisible, cb.labelvisible)
        a[] = show
    end
    cb.width[] = show ? Makie.Auto() : 0
    Makie.colgap!(c.fig.layout, show ? 12 : 0)
    return
end

function _draw_selection!(c::SceneCanvas, sess::DesignSession)
    t = sess.target.target[]
    name = selected_component(sess)
    drawn = c.drawn[]
    pts = Point3f[]
    if name !== nothing && length(drawn) == length(t.tri_area)
        ci = findfirst(x -> x.name === name, t.components)
        # a part left out of the calculations is highlighted on its ghost
        g = sess.target.ghost[]
        tt = ci === nothing && g !== nothing ? g : t
        # an entrance hole with its closing lip
        cis = findall(x -> x.name in with_lips((name,)), tt.components)
        for k in eachindex(tt.tri_area)
            tt.tri_comp[k] in cis || continue
            # apertures and cut-away parts are not drawn; highlight them anyway
            v0, e1, e2 = tt.tri_v0[k], tt.tri_e1[k], tt.tri_e2[k]
            push!(pts, Point3f(v0...), Point3f((v0 + e1)...), Point3f((v0 + e2)...))
        end
    end
    if isempty(pts)
        _update!(c.highlight, _NOTRI...; visible = false)
    else
        _update!(c.highlight, pts, [GLTriangleFace(3i - 2, 3i - 1, 3i) for i in 1:length(pts)÷3];
                 visible = true)
    end
    _draw_selected_beams!(c, sess)
    c.dirty[] = true
    return c
end

# ---------------------------------------------------------------- beams

_beam_length(sess::DesignSession) =
    (L = sess.view.beam_length[]; L > 0 ? L : 2 * sess.target.target[].scale)

function _beam_colors(sess::DesignSession)
    bs = sess.beams.beams[]
    en = sess.beams.enabled[]
    pal = cluster_colors(sess.beams)
    status = _beam_status(sess)
    hid = sess.view.hidden_beams[]
    return map(zip(bs, en, sess.beams.cluster[])) do (b, on, g)
        b.name in hid && return RGBAf(0, 0, 0, 0)               # not drawn
        on || return sess.view.show_off[] ? BEAM_OFF : RGBAf(0, 0, 0, 0)
        status !== nothing && return get(STATUS_COLORS, get(status, b.name, :none), NAN_COLOR)
        _hexcolor(pal[g])
    end
end

function _draw_beams!(c::SceneCanvas, sess::DesignSession)
    t = sess.target.target[]
    bs = sess.beams.beams[]
    style = sess.view.beams[]
    L = _beam_length(sess)
    paths = isempty(bs) ? BeamPath[] : beam_paths(t, bs; length = L)
    c.paths[] = paths
    cols = _beam_colors(sess)
    if style === :off || isempty(bs)
        _update!(c.rays, _NOSEG3; color = [NAN_COLOR, NAN_COLOR], visible = false)
        _update!(c.aims, _NOPOINT3; color = [NAN_COLOR], visible = false)
        _update!(c.hits, _NOPOINT3; color = [NAN_COLOR], visible = false)
    else
        segs = Point3f[]
        for p in paths
            push!(segs, Point3f(p.start...), Point3f(p.stop...))
        end
        _update!(c.rays, segs; color = repeat(cols; inner = 2), visible = style === :rays,
                 linewidth = 1.5)
        _update!(c.aims, [Point3f(p.aim...) for p in paths]; color = cols, visible = true)
        landed = [p.component !== :missed for p in paths]
        # landing markers just in front of the surface, so they do not show through it
        δ = 2e-3 * t.scale
        hp = [Point3f((p.stop - δ * normalize(p.stop - p.start))...) for (p, l) in zip(paths, landed) if l]
        _update!(c.hits, isempty(hp) ? _NOPOINT3 : hp;
                 color = isempty(hp) ? [NAN_COLOR] : cols[landed], visible = true)
    end
    if style === :volumes && !isempty(bs)
        verts = Point3f[]
        faces = GLTriangleFace[]
        vcols = RGBAf[]
        for (b, col) in zip(bs, cols)
            b.name in sess.view.hidden_beams[] && continue
            col.alpha == 0 && continue                          # an off beam not drawn
            v, f = beam_volume(t, b; length = L)
            o = length(verts)
            append!(verts, Point3f.(v))
            append!(faces, (GLTriangleFace(x[1] + o, x[2] + o, x[3] + o) for x in f))
            append!(vcols, fill(RGBAf(col.r, col.g, col.b, col.alpha * VOLUME_ALPHA), length(v)))
        end
        _update!(c.volumes, verts, faces; color = vcols, visible = true)
    else
        _update!(c.volumes, _NOTRI...; color = fill(RGBAf(0, 0, 0, 0), 3), visible = false)
    end
    _draw_selected_beams!(c, sess)
    c.dirty[] = true
    return c
end

function _draw_selected_beams!(c::SceneCanvas, sess::DesignSession)
    idx = selected_beams(sess)
    paths = c.paths[]
    segs = Point3f[]
    for i in idx
        i <= length(paths) || continue
        push!(segs, Point3f(paths[i].start...), Point3f(paths[i].stop...))
    end
    if isempty(segs) || sess.view.beams[] === :off
        _update!(c.selected_rays, _NOSEG3; visible = false)
    else
        _update!(c.selected_rays, segs; visible = true)
    end
    return c
end

# Each beam's worst status in the report the beams are colored by (the TaLIS stay-in
# rows are per hemisphere, not per beam).
function _beam_status(sess::DesignSession)
    mode = sess.view.beam_color[]
    mode === :group && return nothing
    r = mode === :talis ? sess.talis.report[] : sess.review.report[]
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
        _update!(c.outlines, _NOSEG3; color = [NAN_COLOR, NAN_COLOR])
        c.dirty[] = true
        return c
    end
    R = ts.params[][:radius]
    segs = Point3f[]
    cols = RGBAf[]
    narc = 12
    for w in ws, r in (R - 0.01, R + 0.01)
        p(φ, z) = Point3f(r * cos(φ), r * sin(φ), z)
        col = w.open ? WINDOW_COLORS.open : WINDOW_COLORS.closed
        if w.round
            # a circle on the unrolled wall: arc length R φ by z
            φ0, a = deg2rad(w.phi), w.width / 2
            q(t) = p(φ0 + a * cos(t) / R, w.z + a * sin(t))
            tq = range(0, 2π; length = 4narc + 1)
            for k in 1:4narc
                push!(segs, q(tq[k]), q(tq[k+1]))
                push!(cols, col, col)
            end
            continue
        end
        φa = deg2rad(w.phi) - w.width / 2R
        Δφ = w.width / R
        za, zb = w.z - w.height / 2, w.z + w.height / 2
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

# ---------------------------------------------------------------- triad

const _TRIAD_PX = 84

# A small scene in the viewer's lower left showing the x (red), y (green), z (blue)
# axes as the camera sees them.
function _triad_scene(parent::Scene)
    vp = lift(v -> Rect2i(Makie.origin(v) .+ 6, (_TRIAD_PX, _TRIAD_PX)), parent, Makie.viewport(parent))
    s = Scene(parent; viewport = vp, camera = s -> nothing, clear = false)
    cols = [RGBAf(0.85, 0.2, 0.2, 1), RGBAf(0.2, 0.65, 0.25, 1), RGBAf(0.2, 0.4, 0.9, 1)]
    o = Point3f(0, 0, 0)
    axes = [Point3f(1, 0, 0), Point3f(0, 1, 0), Point3f(0, 0, 1)]
    linesegments!(s, [p for a in axes for p in (o, a)]; color = repeat(cols; inner = 2), linewidth = 3)
    text!(s, [1.3f0 * a for a in axes]; text = ["x", "y", "z"], color = cols, fontsize = 14,
          align = (:center, :center))
    return s
end

function _update_triad!(c::SceneCanvas)
    f, _, u = camera_frame(c.camera)
    V = Makie.lookat(Makie.Vec3d(-3f), Makie.Vec3d(0, 0, 0), u)
    P = Makie.orthographicprojection(-1.6, 1.6, -1.6, 1.6, 0.1, 10.0)
    Makie.set_proj_view!(Makie.camera(c.triad), Makie.Mat4d(P), Makie.Mat4d(V))
    return c
end
