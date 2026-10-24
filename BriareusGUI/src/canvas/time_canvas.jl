# The time panel's chart: what `time_chart` describes (pulses, the drive history, the
# energy split, the retune), on a left and an optional right y axis over time, and the
# time cursor. A click or drag on the chart moves the cursor, which colors the 3D view
# by that step of the history. Lines are created up front, a pool per axis; updates set
# their inputs and hide the unused ones (see scene_canvas.jl). The key is drawn by the
# panel (QML), as the model tree keys the 3D view's colors.

# Series per axis and line style, at most (a line's style is fixed when it is created)
const _TIME_LINES = Dict(:solid => 24, :dash => 12, :dot => 4)           # left axis
const _TIME_RIGHT = Dict(:solid => 4, :dash => 2, :dot => 2)             # right axis

"""
    TimeCanvas

The time chart of the Briareus window: `fig`, the axes `ax` (left) and `ax2` (right),
the line pools `left` and `right` (by line style), the `cursor` line, and a `dirty` flag the shell reads
to request a redraw.
"""
struct TimeCanvas
    fig::Figure
    ax::Axis
    ax2::Axis
    left::Dict{Symbol,Vector{Any}}
    right::Dict{Symbol,Vector{Any}}
    cursor::Any
    dirty::Base.RefValue{Bool}
end

const _LINESTYLES = Dict(:solid => nothing, :dash => :dash, :dot => :dot)

"""
    time_canvas(sess::DesignSession) -> TimeCanvas

Build the time chart for `sess` and keep it in sync with the time panel's state.
"""
function time_canvas(sess::DesignSession)
    fig = Figure(size = (900, 260), figure_padding = (6, 6, 2, 6))
    ticks = WilkinsonTicks(4; k_min = 3)            # a short panel: few ticks
    ax = Axis(fig[1, 1]; xlabel = "t (ns)", ylabel = " ", xgridvisible = false, yticks = ticks)
    ax2 = Axis(fig[1, 1]; yaxisposition = :right, ylabel = " ", xgridvisible = false, ygridvisible = false,
               yticks = ticks)
    hidespines!(ax2)
    hidexdecorations!(ax2)
    linkxaxes!(ax, ax2)
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
    c = TimeCanvas(fig, ax, ax2, left, right, cursor, Ref(true))
    tm = sess.time
    onany((_...) -> _draw_time!(c, sess), tm.chart, tm.history, tm.pulses, tm.retune, sess.beams.beams,
          sess.beams.enabled, sess.beams.cluster)
    on(_ -> _draw_cursor!(c, sess), tm.cursor)
    # a click on the chart, or a drag along it, puts the cursor at that time
    pressed = Ref(false)
    on(events(fig).mousebutton) do e
        e.button == Mouse.left || return Consume(false)
        if e.action == Mouse.press && is_mouseinside(ax.scene)
            pressed[] = true
            _point_cursor!(c, sess)
            return Consume(true)
        end
        e.action == Mouse.release && (pressed[] = false)
        return Consume(false)
    end
    on(events(fig).mouseposition) do _
        pressed[] && _point_cursor!(c, sess)
        return Consume(false)
    end
    _draw_time!(c, sess)
    return c
end

function _point_cursor!(c::TimeCanvas, sess::DesignSession)
    sess.time.history[] === nothing && return
    set_time_cursor!(sess, mouseposition(c.ax.scene)[1])
    return
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

function _draw_time!(c::TimeCanvas, sess::DesignSession)
    ch = time_chart(sess)
    sl = ch.left.series
    sr = ch.right === nothing ? NamedTuple[] : ch.right.series
    _draw_series!(c.left, sl, isempty(sl) ? Point2f(0, 0) : Point2f(sl[1].x[1], sl[1].y[1]))
    _draw_series!(c.right, sr, isempty(sr) ? Point2f(0, 0) : Point2f(sr[1].x[1], sr[1].y[1]))
    ax, ax2 = c.ax, c.ax2
    ax.xlabel[] = ch.xlabel
    ax.ylabel[] = isempty(ch.left.label) ? " " : ch.left.label
    two = !isempty(sr)
    ax2.ylabel[] = two ? ch.right.label : " "
    for a in (ax2.yticklabelsvisible, ax2.yticksvisible, ax2.ylabelvisible)
        a[] = two
    end
    # powers and energies from zero; temperatures and modes where they lie
    ax.limits[] = ch.chart === :drive ? (nothing, nothing) : (nothing, (0, nothing))
    reset_limits!(ax)
    two && reset_limits!(ax2)
    _draw_cursor!(c, sess)
    return c
end

function _draw_cursor!(c::TimeCanvas, sess::DesignSession)
    k = time_step(sess.time)
    if k === nothing
        _update!(c.cursor; visible = false)
    else
        _update!(c.cursor, [sess.time.history[].time[k]]; visible = true)
    end
    c.dirty[] = true
    return c
end
