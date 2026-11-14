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
                     "This target has no capsule, observed face or foil to map"
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
        ax.title[] = m.kind === :capsule ? "Capsule incident x-ray flux" :
                     m.kind === :face ? "Observed face incident x-ray flux" : "Foil laser irradiance"
        c.colorbar.label[] = "W/cm²"
        ax.limits[] = (nothing, nothing)                  # not an image's window
        reset_limits!(ax)
    end
    c.dirty[] = true
    return c
end

"""
    image_canvas() -> MapCanvas

A canvas for one image, the Diagnostics tab's: [`show_image!`](@ref) draws a recorded
imager's image on it.
"""
function image_canvas()
    fig = Figure(size = (700, 300), figure_padding = 6)
    ax = Axis(fig[1, 1]; xlabel = "image x (mm)", ylabel = "image y (mm)", aspect = DataAspect(), title = " ")
    hm = heatmap!(ax, [0.0, 1.0], [0.0, 1.0], fill(NaN32, 2, 2); colormap = FIELD_COLORMAP,
                  colorrange = (0.0, 1.0), nan_color = :transparent)
    cb = Colorbar(fig[1, 2], hm; label = " ")
    return MapCanvas(fig, ax, hm, cb, Ref(true))
end

# The x-ray image (see `xray_image!`): radiance over the image plane, mm by mm.
_draw_image!(c::MapCanvas, sess::DesignSession) =
    show_image!(c, sess.drive.image[]; empty = "Computing the drive for the x-ray image…")

"""
    show_image!(c::MapCanvas, im; empty = " ")

Draw the image `im` (`(; x, y, values, title, unit, window)`, as
`Controllers.diagnostic_image` gives it) on `c`, opened on its region with signal;
`nothing` clears the canvas and titles it `empty`.
"""
function show_image!(c::MapCanvas, im; empty::AbstractString = " ")
    ax = c.ax
    if im === nothing
        _update!(c.heatmap, [0.0, 1.0], [0.0, 1.0], fill(NaN32, 2, 2); colorrange = (0.0, 1.0))
        ax.title[] = empty
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
        c.colorbar.label[] = get(im, :unit, "W/cm²/sr")
        # an imager's image opens on its region with signal: Ctrl + click returns to it,
        # Ctrl + Shift + click (Makie's autolimits) shows the whole field
        w = get(im, :window, nothing)
        ax.limits[] = w === nothing ? (nothing, nothing) : ((w[1], w[2]), (w[3], w[4]))
        reset_limits!(ax)
    end
    c.dirty[] = true
    return c
end
