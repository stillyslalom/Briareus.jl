# The time panel's chart: what `time_chart` describes (pulses, the drive history, the
# energy split, the retune), on a left and an optional right y axis over time, and the
# time cursor. A click or drag on the chart moves the cursor, which colors the 3D view
# by that step of the history. The sweep panel's chart (`sweep_chart`: the least
# clearance and the blocked line of sight against the shift, the selected step as the
# cursor) is the same canvas on another description. Lines are created up front, a pool
# per axis; updates set their inputs and hide the unused ones (see scene_canvas.jl). The
# key is drawn by the panel (QML), as the model tree keys the 3D view's colors.

# Series per axis and line style, at most (a line's style is fixed when it is created)
const _TIME_LINES = Dict(:solid => 24, :dash => 12, :dot => 4)           # left axis
const _TIME_RIGHT = Dict(:solid => 4, :dash => 8, :dot => 2)             # right axis

"""
    TimeCanvas

A chart of the Briareus window (the time chart, or the sweep's): `fig`, the axes `ax`
(left) and `ax2` (right), the line pools `left` and `right` (by line style), the `cursor`
line, a `dirty` flag the shell reads to request a redraw, and `chart` (the function
returning what to plot: `time_chart`'s description, or `sweep_chart`'s).
"""
struct TimeCanvas
    fig::Figure
    ax::Axis
    ax2::Axis
    left::Dict{Symbol,Vector{Any}}
    right::Dict{Symbol,Vector{Any}}
    cursor::Any
    dirty::Base.RefValue{Bool}
    chart::Any
end

const _LINESTYLES = Dict(:solid => nothing, :dash => :dash, :dot => :dot)

# Ticks on a y axis `h` pixels tall: about one per 45 px, two to five. Makie's
# `LinearTicks` takes the count as given (PlotUtils' Wilkinson search with a low maximum
# warns "No strict ticks found" on some spans).
_yticks_for(h) = LinearTicks(clamp(round(Int, h / 45), 2, 5))

"""
    time_canvas(sess::DesignSession) -> TimeCanvas

Build the time chart for `sess` and keep it in sync with the time panel's state.
"""
function time_canvas(sess::DesignSession)
    tm = sess.time
    c = _chart_canvas(() -> time_chart(sess))
    onany((_...) -> _draw_time!(c, sess), tm.chart, tm.history, tm.pulses, tm.retune, sess.beams.beams,
          sess.beams.enabled, sess.beams.cluster)
    on(_ -> _draw_cursor!(c, sess), tm.cursor)
    # a click on the chart, or a drag along it, puts the cursor at that time
    _on_chart_click(x -> sess.time.history[] === nothing || set_time_cursor!(sess, x), c)
    _draw_time!(c, sess)
    return c
end

"""
    sweep_canvas(sess::DesignSession) -> TimeCanvas

Build the sweep panel's chart (`sweep_chart`) for `sess`: a click on it selects the
nearest step.
"""
function sweep_canvas(sess::DesignSession)
    sw = sess.sweep
    c = _chart_canvas(() -> sweep_chart(sess))
    onany((_...) -> _draw_time!(c, sess), sw.result, sw.running, sess.target.target, sess.beams.beams)
    on(_ -> _draw_cursor!(c, sess), sw.selected)
    _on_chart_click(c) do x
        r = sw.result[]
        r === nothing || isempty(r.rows) ||
            select_sweep_step!(sw, argmin(k -> abs(r.rows[k].value - x), eachindex(r.rows)))
    end
    _draw_time!(c, sess)
    return c
end

# A chart's figure, axes and line pools; `chart()` describes what to draw.
function _chart_canvas(chart)
    fig = Figure(size = (900, 260), figure_padding = (6, 6, 2, 6))
    ticks = _yticks_for(220)
    ax = Axis(fig[1, 1]; xlabel = "t (ns)", ylabel = " ", xgridvisible = false, yticks = ticks)
    ax2 = Axis(fig[1, 1]; yaxisposition = :right, ylabel = " ", xgridvisible = false, ygridvisible = false,
               yticks = ticks)
    hidespines!(ax2)
    hidexdecorations!(ax2)
    linkxaxes!(ax, ax2)
    # the panel's height sets the tick count, so that a short chart's labels do not overlap
    on(ax.scene.viewport) do vp
        t = _yticks_for(vp.widths[2])
        t == ax.yticks[] && return
        ax.yticks[] = t
        ax2.yticks[] = t
    end
    # the chart is for reading and pointing at times: no zoom or pan gestures
    for a in (ax, ax2), name in collect(keys(Makie.interactions(a)))
        deregister_interaction!(a, name)
    end
    pt = [Point2f(0, 0)]
    pool(a, n) = Dict(st => Any[lines!(a, pt; color = :black, linestyle = _LINESTYLES[st],
                                       linewidth = st === :solid ? 2.0 : 1.6, visible = false) for _ in 1:n[st]]
                      for st in keys(n))
    left, right = pool(ax, _TIME_LINES), pool(ax2, _TIME_RIGHT)
    cursor = vlines!(ax, [0.0]; color = RGBAf(0.1, 0.1, 0.1, 0.55), linewidth = 1.5, visible = false)
    return TimeCanvas(fig, ax, ax2, left, right, cursor, Ref(true), chart)
end

# A click on the chart, or a drag along it, calls `f(x)` at the pointer's abscissa.
function _on_chart_click(f, c::TimeCanvas)
    pressed = Ref(false)
    on(events(c.fig).mousebutton) do e
        e.button == Mouse.left || return Consume(false)
        if e.action == Mouse.press && is_mouseinside(c.ax.scene)
            pressed[] = true
            f(mouseposition(c.ax.scene)[1])
            return Consume(true)
        end
        e.action == Mouse.release && (pressed[] = false)
        return Consume(false)
    end
    on(events(c.fig).mouseposition) do _
        pressed[] && f(mouseposition(c.ax.scene)[1])
        return Consume(false)
    end
    return c
end

function _set_line!(p, s)
    pts = [Point2f(x, y) for (x, y) in zip(s.x, s.y)]
    isempty(pts) && (pts = [Point2f(0, 0)])
    _update!(p, pts; color = _hexcolor(s.color), visible = true)
    return
end

# Draw `series` with a pool's lines of their style; the lines left over are hidden, on a
# point of the drawn data (out of the autolimits' way).
function _draw_series!(pools, series, at)
    used = Dict(st => 0 for st in keys(pools))
    for s in series
        pool = pools[s.style]
        used[s.style] < length(pool) || continue
        _set_line!(pool[used[s.style] += 1], s)
    end
    for (st, pool) in pools, k in used[st]+1:length(pool)
        _update!(pool[k], [at]; visible = false)
    end
end

# Unicode subscripts (aₗ/a₀), which Makie's font lacks, as rich-text subscripts.
const _SUBSCRIPTS = Dict(zip("₀₁₂₃₄₅₆₇₈₉ₗ", "0123456789l"))
function _axis_label(s::AbstractString)
    any(in(keys(_SUBSCRIPTS)), s) || return s
    parts = Any[]
    for (sub, run) in _runs(c -> haskey(_SUBSCRIPTS, c), s)
        push!(parts, sub ? subscript(map(c -> _SUBSCRIPTS[c], run)) : run)
    end
    return rich(parts...)
end

# Maximal runs of `s` over which `f` is constant, as `(f, run)` pairs.
function _runs(f, s::AbstractString)
    out = Tuple{Bool,String}[]
    for c in s
        v = f(c)
        if !isempty(out) && out[end][1] == v
            out[end] = (v, out[end][2] * c)
        else
            push!(out, (v, string(c)))
        end
    end
    return out
end

function _draw_time!(c::TimeCanvas, sess::DesignSession)
    ch = c.chart()
    sl = ch.left.series
    sr = ch.right === nothing ? NamedTuple[] : ch.right.series
    _draw_series!(c.left, sl, isempty(sl) ? Point2f(0, 0) : Point2f(sl[1].x[1], sl[1].y[1]))
    _draw_series!(c.right, sr, isempty(sr) ? Point2f(0, 0) : Point2f(sr[1].x[1], sr[1].y[1]))
    ax, ax2 = c.ax, c.ax2
    ax.xlabel[] = ch.xlabel
    ax.ylabel[] = isempty(ch.left.label) ? " " : ch.left.label
    two = !isempty(sr)
    ax2.ylabel[] = two ? _axis_label(ch.right.label) : " "
    for a in (ax2.yticklabelsvisible, ax2.yticksvisible, ax2.ylabelvisible)
        a[] = two
    end
    # powers and energies from zero; temperatures, modes and clearances where they lie
    ax.limits[] = ch.chart in (:drive, :sweep) ? (nothing, nothing) : (nothing, (0, nothing))
    reset_limits!(ax)
    two && reset_limits!(ax2)
    _draw_cursor!(c, sess)
    return c
end

# The cursor at the chart's `cursor` (the time chart's step, the sweep's selected step;
# `NaN`: none).
function _draw_cursor!(c::TimeCanvas, sess::DesignSession)
    x = c.chart().cursor
    if isnan(x)
        _update!(c.cursor; visible = false)
    else
        _update!(c.cursor, [x]; visible = true)
    end
    c.dirty[] = true
    return c
end
