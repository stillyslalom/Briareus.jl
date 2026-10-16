# The Target step: a hohlraum or a backlighter foil, its dimensions and wall
# properties, and the mesh resolution. Every edit rebuilds the `Target` (meshing takes
# milliseconds; the view factors come later, on the Drive step).

"""
Target kinds of the Target step: `:hohlraum` (cylindrical hohlraum with LEHs and an
optional capsule, `Briareus.hohlraum`) and `:foil` (a backlighter foil facing TCC with
an optional black sphere standing in for the main target).
"""
const TARGET_KINDS = (:hohlraum, :foil)

"""
Editable dimensions per target kind: `(key, label, unit, default)`. Lengths in mm;
albedo and conversion are fractions.
"""
const TARGET_FIELDS = (
    hohlraum = ((:length, "Length", "mm", 10.1),
                (:radius, "Radius", "mm", 2.875),
                (:leh_radius, "LEH radius", "mm", 1.55),
                (:capsule_radius, "Capsule radius", "mm", 0.9),
                (:wall_albedo, "Wall albedo", "", 0.85),
                (:wall_conversion, "Wall laser → x-ray conversion", "", 0.8)),
    foil = ((:distance, "Distance from TCC", "mm", 6.0),
            (:width, "Foil width", "mm", 2.0),
            (:conversion, "Laser → x-ray conversion", "", 0.02),
            (:main_radius, "Main-target sphere radius", "mm", 1.0)))

"""
Mesh resolutions of the Target step. `:standard` is `hohlraum`'s default mesh (about
3,300 elements); `:coarse` (about 550) keeps the view factors to a few seconds.
"""
const MESH_RESOLUTIONS = (:coarse, :standard, :fine)

const _HOHLRAUM_MESH = Dict(
    :coarse => (; nz = 12, nphi = 16, nr_cap = 3, nr_leh = 2, capsule_ntheta = 8, capsule_nphi = 12),
    :standard => (; nz = 36, nphi = 48, nr_cap = 6, nr_leh = 4, capsule_ntheta = 24, capsule_nphi = 32),
    :fine => (; nz = 48, nphi = 64, nr_cap = 8, nr_leh = 4, capsule_ntheta = 32, capsule_nphi = 48))
const _FOIL_MESH = Dict(:coarse => (16, 8, 12), :standard => (40, 16, 24), :fine => (60, 24, 32))

"Name of the foil component of a `:foil` target; the main-target sphere is `:main`."
const FOIL_NAME = :backlighter

"""
    WindowSpec(; phi = 0.0, z = 0.0, width = 1.0, height = 1.0, open = false,
               albedo = 0.85, conversion = 0.8)

A diagnostic window cut into the hohlraum wall (a core `WallWindow`): centered at
azimuth `phi` (degrees from +x) and axial position `z` (mm), `width` mm of arc by
`height` mm along the axis. A closed window (`open = false`) is an opaque wall patch
with its own x-ray `albedo` and laser → x-ray `conversion`; an open one is a hole
(an aperture: laser light and x-rays leave through it). Window `k` of the table is the
target component `:window_k`.
"""
struct WindowSpec
    phi::Float64
    z::Float64
    width::Float64
    height::Float64
    open::Bool
    albedo::Float64
    conversion::Float64
end
WindowSpec(; phi = 0.0, z = 0.0, width = 1.0, height = 1.0, open::Bool = false,
           albedo = 0.85, conversion = 0.8) =
    WindowSpec(phi, z, width, height, open, albedo, conversion)

"""
Editable window fields: `(key, label, unit)`. Albedo and conversion apply to closed
windows only.
"""
const WINDOW_FIELDS = ((:phi, "φ", "°"), (:z, "z", "mm"), (:width, "Width", "mm of arc"),
                       (:height, "Height", "mm"), (:albedo, "Albedo", ""),
                       (:conversion, "Conversion", ""))

# The core window of a spec (validates the fractions; `WallWindow` checks the size).
function _wall_window(w::WindowSpec, k::Integer)
    if !w.open
        0 <= w.albedo <= 1 || throw(ArgumentError("window $k: the albedo must be between 0 and 1"))
        0 <= w.conversion <= 1 ||
            throw(ArgumentError("window $k: the conversion must be between 0 and 1"))
    end
    w.width > 0 && w.height > 0 ||
        throw(ArgumentError("window $k: the width and height must be positive"))
    mat = w.open ? nothing : Material(:window; albedo = w.albedo, conversion = w.conversion)
    return WallWindow(; phi = w.phi, z = w.z, width = w.width, height = w.height,
                      material = mat, open = w.open)
end

"""
    TargetSetup(; kind = :hohlraum, resolution = :standard, windows = WindowSpec[])

State of the Target step: `kind` (one of [`TARGET_KINDS`](@ref)), `params` (every
[`TARGET_FIELDS`](@ref) value of both kinds, so switching kinds keeps the other's
edits), `resolution` (one of [`MESH_RESOLUTIONS`](@ref)), `windows` (the hohlraum's
diagnostic windows, [`WindowSpec`](@ref)s; a foil ignores them), the built `target`,
and `error` (the last rejected edit). Edit with [`set_target_kind!`](@ref),
[`edit_target!`](@ref), [`set_resolution!`](@ref), [`add_window!`](@ref),
[`edit_window!`](@ref), [`set_window_open!`](@ref) and [`remove_window!`](@ref).
"""
struct TargetSetup
    kind::Observable{Symbol}
    params::Observable{Dict{Symbol,Float64}}
    resolution::Observable{Symbol}
    windows::Observable{Vector{WindowSpec}}
    target::Observable{Target}
    error::Observable{String}
end

function TargetSetup(; kind::Symbol = :hohlraum, resolution::Symbol = :standard,
                     windows::AbstractVector{WindowSpec} = WindowSpec[])
    kind in TARGET_KINDS || throw(ArgumentError("target kind must be one of $TARGET_KINDS"))
    resolution in MESH_RESOLUTIONS ||
        throw(ArgumentError("resolution must be one of $MESH_RESOLUTIONS"))
    params = Dict{Symbol,Float64}(k => v for fs in TARGET_FIELDS for (k, _, _, v) in fs)
    ws = collect(WindowSpec, windows)
    t = build_target(kind, params, resolution, ws)
    return TargetSetup(Observable(kind), Observable(params), Observable(resolution),
                       Observable(ws), Observable(t), Observable(""))
end

"""
    build_target(kind, params, resolution, windows = WindowSpec[]) -> Target
    build_target(ts::TargetSetup) -> Target

The `Target` for a kind, its dimensions (`params`, as in [`TARGET_FIELDS`](@ref)), a
mesh resolution and, for a hohlraum, its diagnostic windows. Throws `ArgumentError`
for dimensions or windows the geometry rejects (overlapping windows, a window past
the end of the wall, …).
"""
function build_target(kind::Symbol, p::AbstractDict, resolution::Symbol,
                      windows::AbstractVector{WindowSpec} = WindowSpec[])
    frac(k, what) = 0 <= p[k] <= 1 || throw(ArgumentError("$what must be between 0 and 1"))
    if kind === :hohlraum
        frac(:wall_albedo, "the wall albedo")
        frac(:wall_conversion, "the wall conversion")
        p[:length] > 0 || throw(ArgumentError("the length must be positive"))
        p[:capsule_radius] >= 0 || throw(ArgumentError("the capsule radius must not be negative"))
        return hohlraum(; length = p[:length], radius = p[:radius], leh_radius = p[:leh_radius],
                        capsule_radius = p[:capsule_radius],
                        wall = gold(albedo = p[:wall_albedo], conversion = p[:wall_conversion]),
                        windows = [_wall_window(w, k) for (k, w) in enumerate(windows)],
                        _HOHLRAUM_MESH[resolution]...)
    elseif kind === :foil
        frac(:conversion, "the conversion")
        p[:width] > 0 || throw(ArgumentError("the foil width must be positive"))
        p[:main_radius] >= 0 || throw(ArgumentError("the sphere radius must not be negative"))
        p[:distance] > p[:main_radius] + p[:width] / 2 ||
            throw(ArgumentError("the foil must sit outside the main-target sphere"))
        n, nt, np = _FOIL_MESH[resolution]
        comps = Component[foil(FOIL_NAME, foil_center(p), (-1, 0, 0); width = p[:width], n,
                               material = Material(:Ti; conversion = p[:conversion]),
                               two_sided = false)]
        p[:main_radius] > 0 &&
            push!(comps, Component(:main, sphere_mesh(p[:main_radius]; ntheta = nt, nphi = np), BLACK))
        return Target(comps)
    end
    throw(ArgumentError("unknown target kind :$kind"))
end
build_target(ts::TargetSetup) = build_target(ts.kind[], ts.params[], ts.resolution[], ts.windows[])

"""
    foil_center(params) -> Vec3

Center of the foil of a `:foil` target: `distance` mm from TCC along +x, facing TCC.
"""
foil_center(p::AbstractDict) = Vec3(p[:distance], 0.0, 0.0)
foil_center(ts::TargetSetup) = foil_center(ts.params[])

# Rebuild after a change; a rejected change is undone by the caller.
function _rebuild!(ts::TargetSetup)
    ts.target[] = build_target(ts)
    ts.error[] = ""
    return ts
end

"""
    set_target_kind!(ts::TargetSetup, kind::Symbol)

Switch between a hohlraum and a foil (one of [`TARGET_KINDS`](@ref)).
"""
function set_target_kind!(ts::TargetSetup, kind::Symbol)
    kind in TARGET_KINDS || throw(ArgumentError("target kind must be one of $TARGET_KINDS"))
    ts.kind[] == kind && return ts
    ts.kind[] = kind
    return _rebuild!(ts)
end

"""
    set_resolution!(ts::TargetSetup, resolution::Symbol)

Mesh resolution, one of [`MESH_RESOLUTIONS`](@ref).
"""
function set_resolution!(ts::TargetSetup, resolution::Symbol)
    resolution in MESH_RESOLUTIONS ||
        throw(ArgumentError("resolution must be one of $MESH_RESOLUTIONS"))
    ts.resolution[] == resolution && return ts
    ts.resolution[] = resolution
    return _rebuild!(ts)
end

"""
    edit_target!(ts::TargetSetup, key::Symbol, text) -> Bool

Set dimension `key` (see [`TARGET_FIELDS`](@ref)) from text-field text. A value the
geometry rejects leaves the target unchanged, sets `error`, and returns `false`.
"""
function edit_target!(ts::TargetSetup, key::Symbol, text)
    haskey(ts.params[], key) || throw(ArgumentError("unknown target field :$key"))
    label = _target_label(key)
    old = ts.params[][key]
    try
        v = _parse_number(text, label)
        v == old && (ts.error[] = ""; return true)
        p = copy(ts.params[])
        p[key] = v
        t = build_target(ts.kind[], p, ts.resolution[], ts.windows[])   # validate before committing
        ts.params[] = p
        ts.target[] = t
        ts.error[] = ""
        return true
    catch err
        err isa ArgumentError || rethrow()
        ts.error[] = _errmsg(err)
        return false
    end
end

function _target_label(key::Symbol)
    for fs in TARGET_FIELDS, (k, label, _, _) in fs
        k === key && return label
    end
    return String(key)
end

# ---------------------------------------------------------------- windows

# Commit a new window list if the hohlraum accepts it; otherwise keep the target and
# report why. Windows only exist on a hohlraum.
function _set_windows!(ts::TargetSetup, ws::Vector{WindowSpec})
    try
        ts.kind[] === :hohlraum ||
            throw(ArgumentError("diagnostic windows are cut into a hohlraum wall"))
        t = build_target(ts.kind[], ts.params[], ts.resolution[], ws)
        ts.windows[] = ws
        ts.target[] = t
        ts.error[] = ""
        return true
    catch err
        err isa ArgumentError || rethrow()
        ts.error[] = _errmsg(err)
        return false
    end
end

function _check_window_index(ts::TargetSetup, i::Integer)
    1 <= i <= length(ts.windows[]) || throw(ArgumentError("no window $i"))
    return i
end

"""
    add_window!(ts::TargetSetup; kwargs...) -> Bool

Add a diagnostic window to the hohlraum wall: a [`WindowSpec`](@ref) from `kwargs`,
or by default a closed 1 × 1 mm window at the hohlraum's midplane with the wall's
albedo and conversion, at the first of φ = 90°, 270°, 0°, 180°, … where it does not
overlap another window. A window the geometry rejects sets `error` and returns
`false`.
"""
function add_window!(ts::TargetSetup; kwargs...)
    p = ts.params[]
    base = (; z = 0.0, width = 1.0, height = 1.0, open = false, albedo = p[:wall_albedo],
            conversion = p[:wall_conversion])
    if haskey(kwargs, :phi)
        return _set_windows!(ts, [ts.windows[]; WindowSpec(; base..., kwargs...)])
    end
    for φ in (90.0, 270.0, 0.0, 180.0, 45.0, 135.0, 225.0, 315.0)
        ws = [ts.windows[]; WindowSpec(; base..., phi = φ, kwargs...)]
        try
            ts.kind[] === :hohlraum && build_target(ts.kind[], p, ts.resolution[], ws)
        catch err
            err isa ArgumentError || rethrow()
            continue
        end
        return _set_windows!(ts, ws)
    end
    ts.error[] = "no free place for another window at the default positions: add one with an explicit φ"
    return false
end

"""
    remove_window!(ts::TargetSetup, i)

Remove window `i` (later windows move up one: window `k` is component `:window_k`).
"""
function remove_window!(ts::TargetSetup, i::Integer)
    _check_window_index(ts, i)
    ws = copy(ts.windows[])
    deleteat!(ws, i)
    ts.kind[] === :hohlraum || (ts.windows[] = ws; return ts)
    _set_windows!(ts, ws)
    return ts
end

"""
    edit_window!(ts::TargetSetup, i, key, text) -> Bool

Set field `key` (one of [`WINDOW_FIELDS`](@ref)) of window `i` from text-field text. A
value the geometry rejects (overlapping windows, a window past the end of the wall,
a fraction outside [0, 1]) leaves the target unchanged, sets `error`, and returns
`false`.
"""
function edit_window!(ts::TargetSetup, i::Integer, key::Symbol, text)
    any(f -> f[1] === key, WINDOW_FIELDS) || throw(ArgumentError("unknown window field :$key"))
    _check_window_index(ts, i)
    try
        label = only(f[2] for f in WINDOW_FIELDS if f[1] === key)
        v = _parse_number(text, "the window $(lowercase(label))")
        w = ts.windows[][i]
        getfield(w, key) == v && (ts.error[] = ""; return true)
        ws = copy(ts.windows[])
        ws[i] = WindowSpec((f === key ? v : getfield(w, f) for f in fieldnames(WindowSpec))...)
        ok = _set_windows!(ts, ws)
        ts.error[] = replace(ts.error[], "ArgumentError: " => "")
        # say which window, unless the geometry's message names them
        ok || occursin(r"windows? \d", ts.error[]) || (ts.error[] = "window $i: " * ts.error[])
        return ok
    catch err
        err isa ArgumentError || rethrow()
        ts.error[] = _errmsg(err)
        return false
    end
end

"""
    set_window_open!(ts::TargetSetup, i, open::Bool) -> Bool

Make window `i` an open hole (`true`) or a closed wall patch with its own albedo and
conversion (`false`).
"""
function set_window_open!(ts::TargetSetup, i::Integer, open::Bool)
    _check_window_index(ts, i)
    w = ts.windows[][i]
    w.open == open && return true
    ws = copy(ts.windows[])
    ws[i] = WindowSpec(w.phi, w.z, w.width, w.height, open, w.albedo, w.conversion)
    return _set_windows!(ts, ws)
end

"Target component name of window `i`."
window_name(i::Integer) = Symbol(:window_, i)

"""
    target_summary(ts::TargetSetup) -> String

Element count and component list, for the step rail and page.
"""
function target_summary(ts::TargetSetup)
    t = ts.target[]
    kind = ts.kind[] === :hohlraum ? "hohlraum" : "foil"
    nw = ts.kind[] === :hohlraum ? length(ts.windows[]) : 0
    nw > 0 && (kind *= nw == 1 ? " with 1 window" : " with $nw windows")
    return "$kind · $(nelements(t)) elements · " *
           join((String(c.name) for c in t.components), ", ")
end

"Names of a target's components, in order."
component_names(t::Target) = [c.name for c in t.components]

"Names of a target's opaque (non-aperture) components."
surface_names(t::Target) = [c.name for c in t.components if c.kind === :surface]
