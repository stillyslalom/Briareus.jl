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
    TargetSetup(; kind = :hohlraum, resolution = :standard)

State of the Target step: `kind` (one of [`TARGET_KINDS`](@ref)), `params` (every
[`TARGET_FIELDS`](@ref) value of both kinds, so switching kinds keeps the other's
edits), `resolution` (one of [`MESH_RESOLUTIONS`](@ref)), the built `target`, and
`error` (the last rejected edit). Edit with [`set_target_kind!`](@ref),
[`edit_target!`](@ref) and [`set_resolution!`](@ref).
"""
struct TargetSetup
    kind::Observable{Symbol}
    params::Observable{Dict{Symbol,Float64}}
    resolution::Observable{Symbol}
    target::Observable{Target}
    error::Observable{String}
end

function TargetSetup(; kind::Symbol = :hohlraum, resolution::Symbol = :standard)
    kind in TARGET_KINDS || throw(ArgumentError("target kind must be one of $TARGET_KINDS"))
    resolution in MESH_RESOLUTIONS ||
        throw(ArgumentError("resolution must be one of $MESH_RESOLUTIONS"))
    params = Dict{Symbol,Float64}(k => v for fs in TARGET_FIELDS for (k, _, _, v) in fs)
    t = build_target(kind, params, resolution)
    return TargetSetup(Observable(kind), Observable(params), Observable(resolution),
                       Observable(t), Observable(""))
end

"""
    build_target(kind, params, resolution) -> Target
    build_target(ts::TargetSetup) -> Target

The `Target` for a kind, its dimensions (`params`, as in [`TARGET_FIELDS`](@ref)) and
a mesh resolution. Throws `ArgumentError` for dimensions the geometry rejects.
"""
function build_target(kind::Symbol, p::AbstractDict, resolution::Symbol)
    frac(k, what) = 0 <= p[k] <= 1 || throw(ArgumentError("$what must be between 0 and 1"))
    if kind === :hohlraum
        frac(:wall_albedo, "the wall albedo")
        frac(:wall_conversion, "the wall conversion")
        p[:length] > 0 || throw(ArgumentError("the length must be positive"))
        p[:capsule_radius] >= 0 || throw(ArgumentError("the capsule radius must not be negative"))
        return hohlraum(; length = p[:length], radius = p[:radius], leh_radius = p[:leh_radius],
                        capsule_radius = p[:capsule_radius],
                        wall = gold(albedo = p[:wall_albedo], conversion = p[:wall_conversion]),
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
build_target(ts::TargetSetup) = build_target(ts.kind[], ts.params[], ts.resolution[])

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
        t = build_target(ts.kind[], p, ts.resolution[])     # validate before committing
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

"""
    target_summary(ts::TargetSetup) -> String

Element count and component list, for the step rail and page.
"""
function target_summary(ts::TargetSetup)
    t = ts.target[]
    kind = ts.kind[] === :hohlraum ? "hohlraum" : "foil"
    return "$kind · $(nelements(t)) elements · " *
           join((String(c.name) for c in t.components), ", ")
end

"Names of a target's components, in order."
component_names(t::Target) = [c.name for c in t.components]

"Names of a target's opaque (non-aperture) components."
surface_names(t::Target) = [c.name for c in t.components if c.kind === :surface]
