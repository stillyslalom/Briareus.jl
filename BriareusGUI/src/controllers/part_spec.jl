# This file defines the windows of a hohlraum part, the part kinds and their fields, the profile
# and pieces of a hohlraum part, the frames of the parts, and the components a part makes.

# ---------------------------------------------------------------- windows

"""
    WindowSpec(; phi = 0.0, z = 0.0, width = 1.0, height = 1.0, open = false,
               albedo = 0.8, conversion = 0.9, round = true)

A diagnostic window cut into a hohlraum part's wall (a core `WallWindow`): centered at
azimuth `phi` (degrees from the part's +x) and axial position `z` (mm, in the part's
frame); `round` (the default, as diagnostic windows usually are) a disk of diameter
`width` mm on the unrolled wall, otherwise a rectangle `width` mm of arc by `height` mm
along the axis. A closed window (`open = false`) is an opaque wall patch with a separate
x-ray `albedo` and laser-to-x-ray `conversion` (by default the hohlraum wall's, core
`gold()`); an open one is a hole (an aperture: laser light and x-rays leave through it).
Window `k` of part `name` is the component `<name>_window_k`.
"""
struct WindowSpec
    phi::Float64
    z::Float64
    width::Float64
    height::Float64
    open::Bool
    albedo::Float64
    conversion::Float64
    round::Bool
end
WindowSpec(; phi = 0.0, z = 0.0, width = 1.0, height = 1.0, open::Bool = false,
           albedo = 0.8, conversion = 0.9, round::Bool = true) =
    WindowSpec(phi, z, width, height, open, albedo, conversion, round)

"""
Editable window fields: `(key, label, unit)`. Albedo and conversion apply to closed
windows only; a round window's `width` is its diameter and it has no `height`.
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
    w.round && return WallWindow(; phi = w.phi, z = w.z, diameter = w.width, material = mat,
                                 open = w.open)
    return WallWindow(; phi = w.phi, z = w.z, width = w.width, height = w.height,
                      material = mat, open = w.open)
end

# ---------------------------------------------------------------- parts

"""
Kinds of free parts: `:sphere`, `:disk`, `:annulus` (a washer), `:rectangle`,
`:cylinder` (open tube), `:cone` (open frustum), `:mesh` (a triangle mesh from an STL
file, such as a CAD export), `:visrad` (an object of a VisRad workspace) and `:hohlraum`
(a cylindrical hohlraum, `Briareus.hohlraum`, in pieces: see [`hohlraum_pieces`](@ref)).
`:frame` is a user frame: a named frame with no geometry, placed like a part
([`add_frame!`](@ref), [`add_port_frame!`](@ref)); the model tree lists it under Frames.
`:assembly` is an assembly of the model tree: a named frame with no geometry,
placed like a part, whose members are the parts and assemblies that name it as their
`assembly` ([`add_assembly!`](@ref)).
"""
const PART_KINDS = (:sphere, :disk, :annulus, :rectangle, :cylinder, :cone, :mesh, :visrad, :hohlraum, :frame,
                    :assembly)

const PART_KIND_LABELS = Dict(:sphere => "Sphere", :disk => "Disk", :annulus => "Washer (annulus)",
                              :rectangle => "Rectangle", :cylinder => "Cylinder (tube)",
                              :cone => "Cone (frustum)", :mesh => "Mesh (STL file)",
                              :visrad => "VisRad object", :hohlraum => "Hohlraum", :frame => "Frame",
                              :assembly => "Assembly")

# Whether a part kind, or part `p`, has no geometry: a user frame or an assembly.
_no_geometry(kind::Symbol) = kind === :frame || kind === :assembly
_no_geometry(p) = _no_geometry(p.kind)

"""
Size fields of each part kind: `(key, label, unit, default)`. The part's local axis is
its `z`: a disk, washer or rectangle faces `+z` (its front), a cylinder or cone runs
along `z` centered on the part's center (a cone from `r1` at `−length/2` to `r2` at
`+length/2`).
"""
const PART_SIZE_FIELDS = Dict(
    :sphere => ((:r1, "Radius", "mm", 0.5),),
    :disk => ((:r1, "Radius", "mm", 0.5),),
    :annulus => ((:r1, "Inner radius", "mm", 0.3), (:r2, "Outer radius", "mm", 0.8)),
    :rectangle => ((:r1, "Width", "mm", 1.0), (:r2, "Height", "mm", 1.0)),
    :cylinder => ((:r1, "Radius", "mm", 0.5), (:length, "Length", "mm", 2.0)),
    :cone => ((:r1, "Radius at −L/2", "mm", 0.2), (:r2, "Radius at +L/2", "mm", 0.6),
              (:length, "Length", "mm", 1.0)),
    :mesh => ((:scale, "Scale", "mm per file unit", 1.0),),
    :visrad => (),
    :frame => (),
    :assembly => (),
    :hohlraum => ((:length, "Length", "mm", 10.1), (:radius, "Radius", "mm", 2.875),
                  (:leh_radius, "LEH radius", "mm", 1.55)))

"""
Values of a `:hohlraum` part besides its sizes (`key => default`): its wall's albedo and
laser-to-x-ray conversion (the placement fields `albedo` and `conversion`, which default
to 0.8 and 0.9 for a hohlraum: core `gold()`, the static NIF wall of Young et al., PoP 27,
082702 (2020); sessions keep the values they saved), each end's opening (`top_end`,
`bottom_end`: 1 an entrance hole, 0 closed: the cap covers the end, 2 open: neither cap
nor hole, the cavity open to a package placed there or to nothing, as a halfraum's; see
[`HOHLRAUM_ENDS`](@ref)), `bottom_leh_radius` (mm: the bottom entrance hole's radius; 0,
the default, `leh_radius`, which is the top's), the closing entrance holes:
`leh_closure` (1 on, 0 off) with the [`LEH_CLOSURE_FIELDS`](@ref), and the wall's profile
(`profile`, see [`HOHLRAUM_PROFILES`](@ref)): for a frustraum `end_radius` (mm, the top
end's radius; `radius` is the waist's) and `bottom_end_radius` (mm; 0, the default:
`end_radius`).
"""
const HOHLRAUM_VALUES = (:albedo => albedo(gold()), :conversion => gold().conversion, :top_end => 1.0, :bottom_end => 1.0,
                         :bottom_leh_radius => 0.0, :leh_closure => 0.0,
                         (k => d for (k, _, _, d) in LEH_CLOSURE_FIELDS)...,
                         :profile => 0.0, :end_radius => 0.0, :bottom_end_radius => 0.0)

"""
The wall profiles of a hohlraum part, by value of its `profile`
([`HOHLRAUM_VALUES`](@ref)): 0 `:cylinder` (one radius, `radius`), 1 `:frustraum` (two
cones meeting at the waist: `radius` at the waist, `end_radius` and `bottom_end_radius`
at the ends, core `frustraum_profile`) and 2 `:sections` (the part's `sections`, a radius
per break, core `HohlraumProfile`). See [`hohlraum_profile`](@ref).
"""
const HOHLRAUM_PROFILES = (:cylinder, :frustraum, :sections)

# The profile kind of hohlraum values `v`.
function _profile_kind(v)
    x = v[:profile]
    x in (0.0, 1.0, 2.0) || throw(ArgumentError("a profile is a cylinder (0), a frustraum (1) or sections (2)"))
    return HOHLRAUM_PROFILES[Int(x)+1]
end

"""
    hohlraum_profile(p::PartSpec) -> Union{Nothing,HohlraumProfile}

The wall profile of hohlraum part `p` for core `hohlraum(; profile)`: `nothing` for a
cylinder, `frustraum_profile` of its length, waist radius (`radius`) and end radii for a
frustraum, and for sections a `HohlraumProfile` of its `sections` with the first and
last breaks moved to `∓length/2`. Throws `ArgumentError` for a profile the core refuses
(fewer than two breaks, unsorted heights, a radius not positive).
"""
function hohlraum_profile(p)                    # p::PartSpec
    v = p.values
    k = _profile_kind(v)
    k === :cylinder && return nothing
    L = v[:length]
    if k === :frustraum
        re = v[:end_radius]
        rb = v[:bottom_end_radius] > 0 ? v[:bottom_end_radius] : re
        re > 0 || throw(ArgumentError("a frustraum needs its end radius"))
        return frustraum_profile(; length = L, waist_radius = v[:radius], end_radius = re, bottom_end_radius = rb)
    end
    s = p.sections
    length(s) >= 2 || throw(ArgumentError("a section profile needs at least two breaks"))
    z = [first(b) for b in s]
    z[1], z[end] = -L / 2, L / 2
    return HohlraumProfile(z, [last(b) for b in s])
end

# The length of hohlraum part `p` (its ends at ±length/2 in its frame) and the wall radius
# at height `z` (the narrower side of a step).
_hohlraum_length(p) = p.values[:length]
function _wall_radius_at(p, z)
    prof = try
        hohlraum_profile(p)
    catch err
        err isa ArgumentError || rethrow()
        nothing                                 # a profile the geometry refuses
    end
    prof === nothing && return p.values[:radius]
    r = Briareus.GUIAPI.profile_radius(prof, z)
    return isfinite(r) ? r : p.values[:radius]
end

# The profile a hohlraum part's sections start from when switched to a section list:
# the breaks of its current profile, a cylinder split in four equal sections.
function _default_sections(p)
    prof = hohlraum_profile(_with(p; values = merge(p.values, Dict(:profile => min(p.values[:profile], 1.0)))))
    L, R = p.values[:length], p.values[:radius]
    prof === nothing && return [(z, R) for z in (-L / 2, -L / 4, 0.0, L / 4, L / 2)]
    return collect(zip(prof.z, prof.r))
end

"""
Values of a `:sphere` part besides its radius (`key => default`): `capsule_mesh` (1:
meshed as a capsule, the geodesic mesh of `Briareus.hohlraum`'s capsule with elements of
the target's capsule size, `MeshSpec.capsule_h`, unless the part sets a separate `h`; 0,
the default: the latitude–longitude `sphere_mesh` with elements of `h`). The starter
hohlraum's capsule sets it, so it meshes as the former hohlraum template's capsule did.
"""
const SPHERE_VALUES = (:capsule_mesh => 0.0,)

"""
The openings of a hohlraum part's end, by value of its `top_end` or `bottom_end`
([`HOHLRAUM_VALUES`](@ref)): 0 `:closed`, 1 `:leh` (an entrance hole), 2 `:open` (core
`hohlraum(; top_end = :open)`: no cap, no hole). A hohlraum part with an open end is a
halfraum.
"""
const HOHLRAUM_ENDS = (:closed, :leh, :open)

# The opening of end `e` (`:top_end` or `:bottom_end`) of hohlraum values `v`.
function _end_kind(v, e::Symbol)
    x = v[e]
    x in (0.0, 1.0, 2.0) || throw(ArgumentError("an end is closed (0), an entrance hole (1) or open (2)"))
    return HOHLRAUM_ENDS[Int(x)+1]
end

# The values a kind has besides the placement and size fields.
_kind_values(kind::Symbol) = kind === :hohlraum ? HOHLRAUM_VALUES : kind === :sphere ? SPHERE_VALUES : ()

"""
    hohlraum_pieces(p::PartSpec) -> Vector{Symbol}

The piece names of a `:hohlraum` part `p` without its prefix, in the order of its
components (each component is named `<part>_<piece>` and carries `part = <part>`):
`:wall` (with a profile of several sections, `:wall_1`, `:wall_2`, … from the bottom),
then `:cap_top` and `:cap_bottom` for each end that is not open, then for each end with
an entrance hole `:leh_top` or `:leh_bottom` (and its closing lip `:leh_top_lip`, … with
`leh_closure`), then the windows `:window_1`, ….
"""
function hohlraum_pieces(p)                     # p::PartSpec (defined below)
    p.kind === :hohlraum || throw(ArgumentError("$(p.name) is not a hohlraum"))
    v = p.values
    out = _wall_pieces(p)
    for (e, cap) in ((:top_end, :cap_top), (:bottom_end, :cap_bottom))
        _end_kind(v, e) === :open || push!(out, cap)
    end
    for (e, leh) in ((:top_end, :leh_top), (:bottom_end, :leh_bottom))
        _end_kind(v, e) === :leh || continue
        push!(out, leh)
        v[:leh_closure] != 0 && push!(out, Symbol(leh, :_lip))
    end
    append!(out, (window_name(k) for k in eachindex(p.windows)))
    return out
end

# The wall pieces of hohlraum part `p`: `:wall`, or one per section of its profile (a
# profile of one cylindrical section is the cylinder's `:wall`, as in core `hohlraum`).
function _wall_pieces(p)
    prof = try
        hohlraum_profile(p)
    catch err
        err isa ArgumentError || rethrow()
        nothing                                 # a profile the geometry refuses
    end
    n = prof === nothing ? 1 : length(prof.z) - 1
    n == 1 && (prof === nothing || prof.r[1] == prof.r[2]) && return [:wall]
    return [Symbol(:wall_, k) for k in 1:n]
end

"The component name of piece `piece` of hohlraum part `name`: `<name>_<piece>`."
piece_name(name::Symbol, piece::Symbol) = Symbol(name, :_, piece)

"""
Placement and material fields of every part: `(key, label, unit, default)`. The part
is rotated by `rotation` about its axis, then tilted to polar angle `theta` and
azimuth `phi` (VisRad's Polar/Azimuthal/Rotation order: the axis then points along
`(θ, φ)`), then moved to `(x, y, z)`. `h` (mm) overrides the mesh element size (0:
the target's; for a VisRad primitive, the VisRad grid, see [`part_size_unit`](@ref)).
A user frame and an assembly have the placement fields only: `albedo`, `conversion` and
`h` belong to a surface.
"""
const PART_FIELDS = ((:x, "x", "mm", 0.0), (:y, "y", "mm", 0.0), (:z, "z", "mm", 0.0),
                     (:theta, "Axis polar θ", "°", 0.0), (:phi, "Axis azimuth φ", "°", 0.0),
                     (:rotation, "Rotation about axis", "°", 0.0),
                     (:albedo, "Albedo", "fraction re-emitted", 0.0), (:conversion, "Conversion", "laser → x-ray", 0.0),
                     (:h, "Element size", "mm (0: target's)", 0.0))

"""
Labels of the named frames that are not a part's: `:tcc`, the target chamber (VisRad's
reference frame), the base every part is placed in by default and beams aim in. Every
part also has a frame named after the part, and a hohlraum part its pieces'
([`part_frames`](@ref)): a part placed in one moves and turns with it.
"""
const FRAME_LABELS = Dict(:tcc => "Target chamber (TCC)")

# The fields of `PART_FIELDS` that describe a surface, which a frame or assembly lacks.
const _SURFACE_FIELDS = (:albedo, :conversion, :h)

"""
    PartGeometry(; two_sided = false, inward = false, aperture = false, file = "",
                 foil = false, observed = false, wall = false, windows = [],
                 sections = [], piece_roles = Dict())

What a part with geometry has beyond its placement, the `geometry` of a
[`PartSpec`](@ref): the flags of its surface, its file, its role flags and, for a
hohlraum, its windows, sections and piece roles. The fields are described with
[`PartSpec`](@ref), which also gives them as properties of the part.
"""
struct PartGeometry
    two_sided::Bool
    inward::Bool
    aperture::Bool
    file::String
    foil::Bool
    observed::Bool
    wall::Bool
    windows::Vector{WindowSpec}
    sections::Vector{NTuple{2,Float64}}
    piece_roles::Dict{Symbol,Symbol}
end
PartGeometry(; two_sided::Bool = false, inward::Bool = false, aperture::Bool = false,
             file::AbstractString = "", foil::Bool = false, observed::Bool = false, wall::Bool = false,
             windows = WindowSpec[], sections = NTuple{2,Float64}[], piece_roles = Dict{Symbol,Symbol}()) =
    PartGeometry(two_sided, inward, aperture, String(file), foil, observed, wall, collect(WindowSpec, windows),
                 NTuple{2,Float64}[(Float64(z), Float64(r)) for (z, r) in sections],
                 Dict{Symbol,Symbol}(Symbol(k) => Symbol(r) for (k, r) in piece_roles))

# What a frame or an assembly answers for a geometry property: no flag, file, window,
# section or piece role.
const _NO_GEOMETRY = PartGeometry()

_with(g::PartGeometry; two_sided = g.two_sided, inward = g.inward, aperture = g.aperture, file = g.file,
      foil = g.foil, observed = g.observed, wall = g.wall, windows = g.windows, sections = g.sections,
      piece_roles = g.piece_roles) =
    PartGeometry(two_sided, inward, aperture, file, foil, observed, wall, windows, sections, piece_roles)

"""
    PartSpec(kind; name, values, two_sided = false, inward = false, aperture = false,
             file = "", frame = :tcc, coords = :cartesian, object = "", steps = [],
             foil = false, observed = false, wall = false, windows = [],
             assembly = TOP_LEVEL, sections = [], piece_roles = Dict())

A free target part: `kind` (one of [`PART_KINDS`](@ref)), its component `name`,
`values` (every [`PART_FIELDS`](@ref) key and its kind's [`PART_SIZE_FIELDS`](@ref)),
`two_sided` (both sides radiate and absorb), `inward` (a sphere, cylinder or cone
radiates from its inside; a flat part from its back; a mesh from the other side) and
`aperture` (a transparent diagnostic surface: rays pass, its incident flux is what
escapes through it). The role flags mark the part for the analyses
([`set_part_role!`](@ref)): `foil` a backlighter foil ([`active_foil`](@ref)),
`observed` the surface whose x-ray drive is analyzed ([`observed_name`](@ref)), `wall`
a re-emitting cavity wall ([`wall_names`](@ref)); the name does not carry a role. A `:mesh`
part reads its triangles from the STL `file`. The
placement fields are in the frame named `frame` (`:tcc`, another part's name or a piece
frame of a hohlraum part: see [`part_frames`](@ref)); `coords` ([`PART_COORDINATES`](@ref)) is
how its position is shown and edited (the values keep `x`, `y`, `z`). The placement
`steps` (`PlacementStep`s: relations to other frames, relative moves) apply in order
after the placement fields; the degrees of freedom they leave free keep the fields'
values ([`part_placement`](@ref)). A `:hohlraum` part has its [`HOHLRAUM_VALUES`](@ref)
too, `windows`, its diagnostic windows ([`WindowSpec`](@ref)s; other kinds have none),
and `sections`, the `(z, r)` breaks of its wall when its profile is a section list
([`hohlraum_profile`](@ref); the first and last `z` follow its length), and
`piece_roles`, the roles its pieces hold in place of their default ones (a piece's name
without the part's prefix => one of [`PIECE_ROLES`](@ref); see [`piece_role`](@ref)).
A `:frame` part is a user frame: its placement fields and steps define the frame and it
has no components. With `object` set it is a port frame ([`is_port_frame`](@ref)): the
frame of the port named `object`, at TCC with `z` toward the port (`theta`, `phi` the
port's angles, `rotation` its roll), placed in `:tcc` without steps.

The struct holds what every node of the model tree has (`kind`, `name`, `values`,
`frame`, `coords`, `object`, `steps`, `assembly`) and `geometry`, a
[`PartGeometry`](@ref) with the flags, the file, the roles, the windows, the sections
and the piece roles. A `:frame` or `:assembly` has no geometry (`geometry` is
`nothing`), its `values` are the placement fields only, and the constructor refuses a
geometry keyword for it. The geometry fields are also properties of the part
(`p.observed`, `p.windows`, …); a frame or an assembly answers `false` or empty.

`assembly` is the part's place in the model tree: the name of the assembly it is a member
of (an `:assembly` part), or [`TOP_LEVEL`](@ref). Membership is separate from placement:
a member is placed in its assembly's frame by default ([`add_assembly!`](@ref),
[`move_into_assembly!`](@ref)) but may be placed in any frame. An `:assembly` part's
`assembly` is its parent. User frames are in no assembly.
"""
struct PartSpec
    kind::Symbol
    name::Symbol
    values::Dict{Symbol,Float64}
    frame::Symbol
    coords::Symbol
    object::String
    steps::Vector{PlacementStep}
    assembly::Symbol
    geometry::Union{Nothing,PartGeometry}
    function PartSpec(kind::Symbol, name::Symbol, values::Dict{Symbol,Float64}, frame::Symbol, coords::Symbol,
                      object::String, steps::Vector{PlacementStep}, assembly::Symbol,
                      geometry::Union{Nothing,PartGeometry})
        _no_geometry(kind) == (geometry === nothing) ||
            throw(ArgumentError(geometry === nothing ? "a $kind part needs its geometry" : "a $kind has no geometry"))
        return new(kind, name, values, frame, coords, object, steps, assembly, geometry)
    end
end

function Base.getproperty(p::PartSpec, s::Symbol)
    hasfield(PartSpec, s) && return getfield(p, s)
    g = getfield(p, :geometry)
    return getfield(g === nothing ? _NO_GEOMETRY : g, s)
end
Base.propertynames(::PartSpec) = (fieldnames(PartSpec)..., fieldnames(PartGeometry)...)

"""
The top level of the model tree: the `assembly` ([`PartSpec`](@ref)) of a part or
assembly that is a member of none, `Symbol("")`.
"""
const TOP_LEVEL = Symbol("")

"""
Coordinates a part's position is shown and edited in: `:cartesian` (x, y, z),
`:cylindrical` (r, φ, z about the frame's z axis) and `:spherical` (r, θ, φ).
"""
const PART_COORDINATES = (:cartesian, :cylindrical, :spherical)

# Position fields of each coordinate system: (key, label, unit).
const PART_POSITION_FIELDS = Dict(
    :cartesian => ((:x, "x", "mm"), (:y, "y", "mm"), (:z, "z", "mm")),
    :cylindrical => ((:pos_r, "r (from the axis)", "mm"), (:pos_phi, "φ (position)", "°"), (:z, "z", "mm")),
    :spherical => ((:pos_r, "r (from the origin)", "mm"), (:pos_theta, "θ (position)", "°"),
                   (:pos_phi, "φ (position)", "°")))

"""
    part_position(p::PartSpec) -> Vector{Pair{Symbol,Float64}}

Position of a part in its coordinates (`coords_of` its center in its frame):
`key => value` for its [`PART_POSITION_FIELDS`](@ref).
"""
part_position(p::PartSpec) =
    [k => v for ((k, _, _), v) in zip(PART_POSITION_FIELDS[p.coords], coords_of(part_center(p), p.coords))]

# x, y, z from a part's position in its coordinates with `key` set to `v` (`point_of`).
function _set_position(p::PartSpec, key::Symbol, v::Real)
    pos = part_position(p)
    k = findfirst(kv -> first(kv) === key, pos)
    k === nothing && throw(ArgumentError("$(p.coords) coordinates have no field :$key"))
    c = Float64[last(kv) for kv in pos]
    c[k] = v
    return Tuple(point_of(c, p.coords))
end

function PartSpec(kind::Symbol; name = kind, values = Dict{Symbol,Float64}(),
                  two_sided::Bool = false, inward::Bool = false, aperture::Bool = false,
                  file::AbstractString = "", frame::Symbol = :tcc, coords::Symbol = :cartesian,
                  object::AbstractString = "", steps = PlacementStep[], foil::Bool = false,
                  observed::Bool = false, wall::Bool = false, windows = WindowSpec[],
                  assembly::Symbol = TOP_LEVEL, sections = NTuple{2,Float64}[],
                  piece_roles = Dict{Symbol,Symbol}())
    kind in PART_KINDS || throw(ArgumentError("part kind must be one of $PART_KINDS"))
    _own_frame_name(Symbol(name), frame) && throw(ArgumentError("a part cannot be placed in the frame named after it"))
    coords in PART_COORDINATES || throw(ArgumentError("coordinates must be one of $PART_COORDINATES"))
    kind === :mesh && isempty(file) && throw(ArgumentError("a mesh part needs an STL file"))
    kind === :visrad && (isempty(file) || isempty(object)) &&
        throw(ArgumentError("a VisRad part needs a workspace file and an object name"))
    bare = _no_geometry(kind)
    bare && (two_sided || inward || aperture || foil || observed || wall || !isempty(file)) &&
        throw(ArgumentError("a $kind has no surface: it takes no flag, file or role"))
    v = Dict{Symbol,Float64}(k => d for (k, _, _, d) in (PART_FIELDS..., PART_SIZE_FIELDS[kind]...)
                             if !(bare && k in _SURFACE_FIELDS))
    for (k, d) in _kind_values(kind)
        v[k] = d
    end
    kind === :hohlraum || isempty(windows) || throw(ArgumentError("only a hohlraum has windows"))
    kind === :hohlraum || isempty(sections) || throw(ArgumentError("only a hohlraum has sections"))
    kind === :hohlraum || isempty(piece_roles) ||
        throw(ArgumentError("piece roles are allowed only on a hohlraum part"))
    all(r -> Symbol(r) in PIECE_ROLES, Base.values(piece_roles)) ||
        throw(ArgumentError("a piece's role must be one of $PIECE_ROLES"))
    for (k, x) in values
        haskey(v, k) || throw(ArgumentError("a $kind has no field :$k"))
        v[k] = x
    end
    any(st -> Symbol(name) in Briareus.GUIAPI.step_references(st), steps) &&
        throw(ArgumentError("a placement step cannot refer to the frame of its part"))
    geometry = bare ? nothing :
               PartGeometry(; two_sided, inward, aperture, file, foil, observed, wall, windows, sections, piece_roles)
    return _port_checked(PartSpec(kind, Symbol(name), v, frame, coords, String(object),
                                  collect(PlacementStep, steps), assembly, geometry))
end

# A port frame stays at TCC, in `:tcc`, without steps (`_with` checks every edit).
function _port_checked(p::PartSpec)
    p.kind === :frame && !isempty(p.object) || return p
    v = p.values
    (p.frame === :tcc && isempty(p.steps) && v[:x] == v[:y] == v[:z] == 0) ||
        throw(ArgumentError("a port frame's origin is TCC and its axis its port's: change its roll, " *
                            "or place another frame in it"))
    return p
end

# Whether `frame` is the frame of part `name` or one of its pieces' (`name.leh_top`, …).
_own_frame_name(name::Symbol, frame::Symbol) = frame === name || startswith(String(frame), String(name) * ".")

# Part `p` with some fields replaced. The other keywords are fields of its geometry
# (`PartGeometry`), which a frame or an assembly lacks.
function _with(p::PartSpec; name = p.name, values = p.values, frame = p.frame, coords = p.coords,
               object = p.object, steps = p.steps, assembly = p.assembly, geometry...)
    g = p.geometry
    if !isempty(geometry)
        g === nothing && throw(ArgumentError("a $(p.kind) has no $(join(keys(geometry), ", "))"))
        g = _with(g; geometry...)
    end
    return _port_checked(PartSpec(p.kind, name, values, frame, coords, object, steps, assembly, g))
end

# A part with the frame named `old` renamed `new`: the frame it is placed in, the frames
# its placement steps refer to, and the assembly it is a member of.
_renamed(p::PartSpec, old::Symbol, new::Symbol) =
    _with(p; frame = p.frame === old ? new : p.frame,
          steps = PlacementStep[_rename_step(st, old, new) for st in p.steps],
          assembly = p.assembly === old ? new : p.assembly)

function _rename_step(st::PlacementStep, old::Symbol, new::Symbol)
    T = typeof(st)
    old in Briareus.GUIAPI.step_references(st) || return st
    return T((f === :ref ? new : getfield(st, f) for f in fieldnames(T))...)
end

"Axis (unit vector) of a part in its frame: polar `theta` and azimuth `phi` in degrees."
part_axis(p::PartSpec) = Vec3(sind(p.values[:theta]) * cosd(p.values[:phi]),
                              sind(p.values[:theta]) * sind(p.values[:phi]), cosd(p.values[:theta]))

"Center of a part in its frame."
part_center(p::PartSpec) = Vec3(p.values[:x], p.values[:y], p.values[:z])

# The part frame in the chamber before its placement steps, placed in its parent's
# chamber frame `frame`: the rotation `euler_rotation(theta, phi, rotation)` (VisRad's
# order) at its center.
_own_frame(p::PartSpec, frame::Frame = CHAMBER) =
    compose(frame, euler_frame(part_center(p), p.values[:theta], p.values[:phi], p.values[:rotation]))

"""
    part_frame_def(p::PartSpec) -> FrameDef

The frame of a part as placed in its parent `p.frame`: origin at its center
([`part_center`](@ref)), axes turned by its `theta`, `phi` and `rotation`
(`euler_rotation`), so its `z` axis is the part's axis; then moved by its placement
`steps`.
"""
part_frame_def(p::PartSpec) =
    FrameDef(p.frame; offset = part_center(p), rotation = (p.values[:theta], p.values[:phi], p.values[:rotation]),
             steps = p.steps)

"""
    part_frames(parts) -> Dict{Symbol,Frame}
    part_frames(ts::TargetSetup) -> Dict{Symbol,Frame}

Every named frame of a design in the chamber: `:tcc`, each part frame
([`part_frame_def`](@ref), named after the part) and its pieces'
([`piece_frame_defs`](@ref)), resolved with `resolve_frames`, so a part placed in
another part's frame follows every ancestor exactly. Throws `ArgumentError` for a part
placed in a frame that does not exist, or parts placed in each other in a loop (the
error names the loop), and for a model tree that is not one ([`check_tree`](@ref)).
"""
function part_frames(parts::AbstractVector)
    check_tree(parts)
    defs = _frame_defs(parts)
    names = Set{Symbol}(first.(defs))
    for p in parts
        (p.frame === :tcc || p.frame in names) && continue
        throw(ArgumentError("no frame named $(p.frame) to place $(p.name) in"))
    end
    return resolve_frames(defs)
end

# The frame definitions of the parts and their pieces.
_frame_defs(parts) = Pair{Symbol,FrameDef}[d for p in parts for d in (p.name => part_frame_def(p), piece_frame_defs(p)...)]

"""
    piece_frame_defs(p::PartSpec) -> Vector{Pair{Symbol,FrameDef}}

The child frames of a part's pieces, placed in the part frame: for a `:hohlraum`
part `name` of length `L`, the entrance holes' `name.leh_top` (origin at the top end's
center `(0, 0, L/2)`, the part's axes, so `z` points out of the hohlraum) and
`name.leh_bottom` (at `(0, 0, −L/2)`, turned 180° about the part's `x`: `z` out of the
hohlraum, `y` reversed, so one placement or aim mirrors between the two ends, as the
former template's entrance-hole frames did), and each window's `name.window_k` (origin at the
window's center on the wall, the part's axes). The frames of the ends are also used for
their caps, since the centroid of a cap is the center of its end, and the wall aims in
the part frame, since its centroid is the center of the part. Other kinds have none.
"""
function piece_frame_defs(p::PartSpec)
    out = Pair{Symbol,FrameDef}[]
    p.kind === :hohlraum || return out
    L = p.values[:length]
    push!(out, Symbol(p.name, ".leh_top") => FrameDef(p.name; offset = (0.0, 0.0, L / 2)),
          Symbol(p.name, ".leh_bottom") => FrameDef(p.name; offset = (0.0, 0.0, -L / 2),
                                                    rotation = (180.0, 0.0, 180.0)))
    for (k, w) in enumerate(p.windows)
        R = _wall_radius_at(p, w.z)
        push!(out, Symbol(p.name, ".", window_name(k)) =>
                   FrameDef(p.name; offset = (R * cosd(w.phi), R * sind(w.phi), w.z)))
    end
    return out
end

"""
    piece_frame(name::Symbol, part::Symbol) -> Symbol

The frame component `name` of part `part` aims in as a piece of it: a hohlraum part's
wall in the part's frame, its caps and entrance holes (and their lips) in its ends'
frames (`part.leh_top`, `part.leh_bottom`), a window in a separate frame (`part.window_k`), see
[`piece_frame_defs`](@ref); any other component in its part's frame.
"""
function piece_frame(name::Symbol, part::Symbol)
    s, pre = String(name), String(part) * "_"
    startswith(s, pre) || return part
    piece = s[ncodeunits(pre)+1:end]
    occursin(r"^(cap|leh)_top", piece) && return Symbol(part, ".leh_top")
    occursin(r"^(cap|leh)_bottom", piece) && return Symbol(part, ".leh_bottom")
    occursin(r"^window_\d+$", piece) && return Symbol(part, ".", piece)
    return part
end

_count(len, h; min = 1, multiple = 1) = multiple * max(cld(min, multiple), ceil(Int, len / h / multiple))

# The azimuthal patch count of a disk, washer, cylinder or cone of radius `r` with
# elements of `h`: core `hohlraum_counts`'s rule (a multiple of 8, at least 8), so a part
# of a hohlraum part's radius meshed with its element size meets the wall's end ring
# at shared vertices (a package disk closing a halfraum's open end leaves no gap at the
# seam; with multiples of 4 the standard size left a 36-gon inside the wall's 40-gon).
_nphi(r, h) = _count(2π * r, h; min = 8, multiple = 8)

"""
    part_component(p::PartSpec, h; frame) -> Component

The core component of a part, meshed with elements of size `h` (mm) unless the part
sets a separate size, carrying its exact curved surface (sphere, disk, washer, cylinder) for
the view factors, and placed in the chamber by `frame`, the chamber `Frame` of the part
(its entry in [`part_frames`](@ref), placement steps included). By default the part is
placed by its placement fields in the chamber frame, without its steps.
"""
function part_component(p::PartSpec, h::Real; frame::Frame = _own_frame(p))
    p.kind === :visrad && throw(ArgumentError("a VisRad object can make several components: use part_components"))
    _no_geometry(p) && throw(ArgumentError("a $(p.kind) has no component"))
    return only(part_components(p, h; frame))
end

"""
    part_components(p::PartSpec, h; frame, raw = false) -> Vector{Component}

The core components of a part, placed by its chamber `frame` as in
[`part_component`](@ref): one for every kind but `:visrad`, whose object can be a
primitive hohlraum (wall, caps and entrance holes). A VisRad object keeps its place
in the chamber (VisRad's frames), moved further by the part's placement; its windings
are repaired and its CAD surface remeshed into elements of the part's size. A VisRad
primitive (rectangle, sphere, cylindrical hohlraum) keeps VisRad's grid, which
reproduces VisRad's element table, unless the part sets a separate element size. With `raw`
a VisRad CAD part keeps its file's facets as they are (no winding repair, no remeshing:
for drawing a part left out of the calculations, at no meshing cost).

A sphere with `capsule_mesh` set (value 1, see [`SPHERE_VALUES`](@ref)) is meshed as a
capsule: the geodesic mesh of `Briareus.hohlraum`'s capsule (`icosphere_mesh`,
`capsule_mesh = :compact`) with elements of `capsule_h` (the target's capsule size,
[`MeshSpec`](@ref)) unless the part sets a separate size. A `:hohlraum` part is core
`hohlraum(...)` without a capsule, meshed with elements of `h`, its components named and
grouped as in [`hohlraum_pieces`](@ref). At the chamber frame it is not moved at all, so
it is identical to the hohlraum of the former template.
"""
function part_components(p::PartSpec, h::Real; frame::Frame = _own_frame(p), raw::Bool = false,
                         capsule_h::Real = h)
    # the same components (identical objects) for the same part, size and frame, so the
    # targets built from them are found again (`_target_of`): an edit elsewhere, an undo
    lock(() -> length(_PART_MEMO) > 512 && empty!(_PART_MEMO), _STL_LOCK)
    hc = _capsule_meshed(p) ? Float64(capsule_h) : 0.0
    return _cached(() -> _part_components(p, h, frame, raw, hc), _PART_MEMO, _part_key(p, h, frame, raw, hc))
end
