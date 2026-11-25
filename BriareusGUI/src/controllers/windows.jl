# This file holds the edits of the diagnostic windows of a hohlraum part.

# ---------------------------------------------------------------- windows

# The index of the hohlraum part whose windows an edit takes: the part named `part`,
# else the design's first hohlraum part. Throws `ArgumentError` when there is none.
function _window_part(ts::TargetSetup, part)
    ps = ts.parts[]
    i = part === nothing ? findfirst(p -> p.kind === :hohlraum, ps) : part_index(ts, Symbol(part))
    i === nothing && throw(ArgumentError(part === nothing ?
                                         "diagnostic windows are cut into a hohlraum wall: insert a hohlraum" :
                                         "no part named $part"))
    ps[i].kind === :hohlraum ||
        throw(ArgumentError("diagnostic windows are cut into a hohlraum wall: $(ps[i].name) is a $(ps[i].kind)"))
    return i
end

"""
    hohlraum_windows(ts::TargetSetup; part = nothing) -> Vector{WindowSpec}

The diagnostic windows of hohlraum part `part` (by default the design's first hohlraum
part; none without one).
"""
function hohlraum_windows(ts::TargetSetup; part = nothing)
    i = part === nothing ? findfirst(p -> p.kind === :hohlraum, ts.parts[]) : part_index(ts, Symbol(part))
    return i === nothing ? WindowSpec[] : ts.parts[][i].windows
end

# Commit a new window list of part `i` if the hohlraum accepts it; otherwise keep the
# target and report why. The report leaves out the part name, since the edit names the
# window.
function _set_windows!(ts::TargetSetup, i::Integer, ws::Vector{WindowSpec})
    ps = copy(ts.parts[])
    name = String(ps[i].name)
    ps[i] = _with(ps[i]; windows = ws)
    ok = _try_build!(ts; parts = ps)
    if !ok
        e = replace(ts.error[], "ArgumentError: " => "")
        ts.error[] = startswith(e, name * ": ") ? e[ncodeunits(name)+3:end] : e
    end
    return ok
end

function _check_window_index(ts::TargetSetup, i::Integer, k::Integer)
    1 <= k <= length(ts.parts[][i].windows) || throw(ArgumentError("no window $k"))
    return k
end

# `f(i)` for the window part of `part`, or `false` with `error` set when there is none.
function _on_window_part(f, ts::TargetSetup, part)
    i = try
        _window_part(ts, part)
    catch err
        err isa ArgumentError || rethrow()
        ts.error[] = _errmsg(err)
        return false
    end
    return f(i)
end

"""
    add_window!(ts::TargetSetup; part = nothing, kwargs...) -> Bool

Add a diagnostic window to the wall of hohlraum part `part` (by default the design's
first hohlraum part): a [`WindowSpec`](@ref) from `kwargs`, or by default a closed round
window 1 mm across at the hohlraum's midplane with the wall's albedo and conversion, at
the first of φ = 90°, 270°, 0°, 180°, … where it does not overlap another window. A
window the geometry rejects, or a design without a hohlraum, sets `error` and returns
`false`.
"""
function add_window!(ts::TargetSetup; part = nothing, kwargs...)
    _on_window_part(ts, part) do i
        p = ts.parts[][i]
        base = (; z = 0.0, width = 1.0, height = 1.0, open = false, albedo = p.values[:albedo],
                conversion = p.values[:conversion])
        haskey(kwargs, :phi) && return _set_windows!(ts, i, [p.windows; WindowSpec(; base..., kwargs...)])
        for φ in (90.0, 270.0, 0.0, 180.0, 45.0, 135.0, 225.0, 315.0)
            ws = [p.windows; WindowSpec(; base..., phi = φ, kwargs...)]
            try
                # The windows are checked as the wall is meshed, on the hohlraum by itself.
                part_components(_with(p; windows = ws), ts.mesh[].h)
            catch err
                err isa ArgumentError || rethrow()
                continue
            end
            return _set_windows!(ts, i, ws)
        end
        ts.error[] = "no free place for another window at the default positions: add one with an explicit φ"
        return false
    end
end

"""
    remove_window!(ts::TargetSetup, k; part = nothing) -> Bool

Remove window `k` of hohlraum part `part` (by default the design's first hohlraum part;
later windows move up one: window `k` is component `<part>_window_k`).
"""
function remove_window!(ts::TargetSetup, k::Integer; part = nothing)
    _on_window_part(ts, part) do i
        _check_window_index(ts, i, k)
        ws = copy(ts.parts[][i].windows)
        deleteat!(ws, k)
        return _set_windows!(ts, i, ws)
    end
end

"""
    edit_window!(ts::TargetSetup, k, key, text; part = nothing) -> Bool

Set field `key` (one of [`WINDOW_FIELDS`](@ref)) of window `k` of hohlraum part `part`
(by default the design's first hohlraum part) from text-field text. A value the geometry
rejects (overlapping windows, a window past the end of the wall, a fraction outside
[0, 1]) leaves the target unchanged, sets `error`, and returns `false`.
"""
function edit_window!(ts::TargetSetup, k::Integer, key::Symbol, text; part = nothing)
    any(f -> f[1] === key, WINDOW_FIELDS) || throw(ArgumentError("unknown window field :$key"))
    _on_window_part(ts, part) do i
        _check_window_index(ts, i, k)
        try
            label = only(f[2] for f in WINDOW_FIELDS if f[1] === key)
            v = _parse_number(text, "the window $(lowercase(label))")
            w = ts.parts[][i].windows[k]
            getfield(w, key) == v && (ts.error[] = ""; return true)
            ws = copy(ts.parts[][i].windows)
            ws[k] = WindowSpec((f === key ? v : getfield(w, f) for f in fieldnames(WindowSpec))...)
            ok = _set_windows!(ts, i, ws)
            ts.error[] = replace(ts.error[], "ArgumentError: " => "")
            # say which window, unless the geometry's message names them
            ok || occursin(r"windows? \d", ts.error[]) || (ts.error[] = "window $k: " * ts.error[])
            return ok
        catch err
            err isa ArgumentError || rethrow()
            ts.error[] = _errmsg(err)
            return false
        end
    end
end

"""
    set_window_open!(ts::TargetSetup, k, open::Bool; part = nothing) -> Bool

Make window `k` of hohlraum part `part` (by default the first) an open hole (`true`) or
a closed wall patch with a separate albedo and conversion (`false`).
"""
function set_window_open!(ts::TargetSetup, k::Integer, open::Bool; part = nothing)
    _on_window_part(ts, part) do i
        _check_window_index(ts, i, k)
        w = ts.parts[][i].windows[k]
        w.open == open && return true
        ws = copy(ts.parts[][i].windows)
        ws[k] = WindowSpec(w.phi, w.z, w.width, w.height, open, w.albedo, w.conversion, w.round)
        return _set_windows!(ts, i, ws)
    end
end

"""
    set_window_round!(ts::TargetSetup, k, round::Bool; part = nothing) -> Bool

Make window `k` of hohlraum part `part` (by default the first) a disk (diameter = its
width) or a rectangle.
"""
function set_window_round!(ts::TargetSetup, k::Integer, round::Bool; part = nothing)
    _on_window_part(ts, part) do i
        _check_window_index(ts, i, k)
        w = ts.parts[][i].windows[k]
        w.round == round && return true
        ws = copy(ts.parts[][i].windows)
        ws[k] = WindowSpec(w.phi, w.z, w.width, w.height, w.open, w.albedo, w.conversion, round)
        return _set_windows!(ts, i, ws)
    end
end

"The name of window `k` among a hohlraum part's pieces (its component is `<part>_window_k`)."
window_name(k::Integer) = Symbol(:window_, k)
