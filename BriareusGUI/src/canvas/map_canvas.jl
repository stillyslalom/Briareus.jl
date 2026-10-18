# The flux map: the capsule's incident x-ray flux over polar angle and azimuth (a
# hohlraum), or the laser irradiance over a foil's face. Plots are created up front;
# updates change their inputs only (see scene_canvas.jl).

"""
    MapCanvas

The flux-map viewer of the Briareus window: `fig`/`ax`, the `heatmap`, the
`colorbar`, and a `dirty` flag the shell reads to request a redraw.
"""
struct MapCanvas
    fig::Figure
    ax::Axis
    heatmap::Any
    colorbar::Colorbar
    dirty::Base.RefValue{Bool}
end

"""
    map_canvas(sess::DesignSession) -> MapCanvas

Build the flux-map viewer for `sess` and keep it in sync with the Drive step.
"""
function map_canvas(sess::DesignSession)
    fig = Figure(size = (900, 600), figure_padding = 8)
    ax = Axis(fig[1, 1]; xlabel = "φ (°)", ylabel = "θ (°)", yreversed = true,
              title = "Compute the drive to see the flux map")
    hm = heatmap!(ax, [0.0, 1.0], [0.0, 1.0], fill(NaN32, 2, 2); colormap = FIELD_COLORMAP,
                  colorrange = (0.0, 1.0), nan_color = :transparent)
    cb = Colorbar(fig[1, 2], hm; label = " ")
    c = MapCanvas(fig, ax, hm, cb, Ref(true))
    onany((_...) -> _draw_map!(c, sess), sess.drive.result, sess.drive.image, sess.drive.view)
    _draw_map!(c, sess)
    return c
end

function _draw_map!(c::MapCanvas, sess::DesignSession)
    sess.drive.view[] === :image && return _draw_image!(c, sess)
    m = flux_map(sess.drive)
    ax = c.ax
    if m === nothing
        _update!(c.heatmap, [0.0, 1.0], [0.0, 1.0], fill(NaN32, 2, 2); colorrange = (0.0, 1.0))
        ax.title[] = sess.drive.result[] === nothing ? "Compute the drive to see the flux map" :
                     "This target has no capsule or foil to map"
        c.colorbar.label[] = " "
    else
        fin = filter(isfinite, m.values)
        lo, hi = isempty(fin) ? (0.0, 1.0) : extrema(fin)
        hi > lo || (hi = lo + max(abs(lo), 1.0) * 1e-6)
        _update!(c.heatmap, m.x, m.y, Float32.(m.values); colorrange = (lo, hi))
        ax.xlabel[] = m.xlabel
        ax.ylabel[] = m.ylabel
        ax.yreversed[] = m.kind === :capsule
        ax.aspect[] = m.kind === :capsule ? AxisAspect(2) : DataAspect()     # equal degrees / mm
        ax.title[] = m.kind === :capsule ? "Capsule incident x-ray flux" : "Foil laser irradiance"
        c.colorbar.label[] = "W/cm²"
        reset_limits!(ax)
    end
    c.dirty[] = true
    return c
end

# The x-ray image (see `xray_image!`): radiance over the image plane, mm by mm.
function _draw_image!(c::MapCanvas, sess::DesignSession)
    im = sess.drive.image[]
    ax = c.ax
    if im === nothing
        _update!(c.heatmap, [0.0, 1.0], [0.0, 1.0], fill(NaN32, 2, 2); colorrange = (0.0, 1.0))
        ax.title[] = "Computing the drive for the x-ray image…"
        c.colorbar.label[] = " "
    else
        fin = filter(isfinite, im.values)
        lo, hi = isempty(fin) ? (0.0, 1.0) : (0.0, maximum(fin))
        hi > lo || (hi = lo + 1.0)
        _update!(c.heatmap, collect(im.x), collect(im.y), Float32.(im.values); colorrange = (lo, hi))
        ax.xlabel[] = "image x (mm)"
        ax.ylabel[] = "image y (mm)"
        ax.yreversed[] = false
        ax.aspect[] = DataAspect()
        ax.title[] = im.title
        c.colorbar.label[] = "W/cm²/sr"
        reset_limits!(ax)
    end
    c.dirty[] = true
    return c
end
