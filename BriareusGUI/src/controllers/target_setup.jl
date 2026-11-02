# The target: parts placed anywhere (hohlraums, spheres, disks, washers, rectangles,
# cylinders, cones, as VisRad's target components; CAD meshes and VisRad objects), the
# starters that insert assemblies of them, and the mesh element sizes. Every edit
# rebuilds the `Target` (meshing takes milliseconds; the view factors are computed on
# request).

"""
    MeshSpec(h, capsule_h)

Element sizes (mm): `h` for hohlraum walls, end caps and LEHs (`hohlraum_counts`) and
the parts; `capsule_h` for the surface whose drive is studied: a sphere meshed as a
capsule (`capsule_counts`, see [`SPHERE_VALUES`](@ref)) and a starter's backlighter foil
when it is inserted. A part may override `h` with its own.
"""
struct MeshSpec
    h::Float64
    capsule_h::Float64
end

"""
Mesh presets (`MeshSpec`s), from the work-precision studies in
`examples/mesh_convergence.jl` and `examples/capsule_mesh.jl` (reference/research-notes.md,
"Hohlraum mesh work-precision" and "Compact capsule mesh") on the NIF hohlraum, whose
capsule is the geodesic `capsule_mesh = :compact`. The wall sizes set capsule P2/P0
errors against the converged meshes of about 0.5 percentage points (`:coarse`, view
factors in under a second), 0.1 (`:standard`, a few seconds) and 0.02 (`:fine`, tens of
seconds: converged for design work). The capsule sizes keep its share within those:
P2/P0 within about 0.17, 0.05 and 0.02 points and P4/P0 within 0.09, 0.03 and 0.01
(320, 980 and 2,880 capsule elements; P2 falls as 1/N on any capsule mesh, P4 ten times
faster on this one than on the lat-long mesh, which needed 0.2, 0.09 and 0.05 mm).
"""
const MESH_PRESETS = (coarse = MeshSpec(0.6, 0.3), standard = MeshSpec(0.28, 0.16),
                      fine = MeshSpec(0.17, 0.09))

"Mesh presets in order; `:custom` marks sizes edited by hand."
const MESH_RESOLUTIONS = (:coarse, :standard, :fine)

_preset_of(m::MeshSpec) = something(findfirst(==(m), MESH_PRESETS), :custom)

"""
Closing entrance holes of a hohlraum part (`Briareus.hohlraum(; leh_closure,
leh_open_radius)`), among its values ([`HOHLRAUM_VALUES`](@ref)): `:leh_closure` (1 on,
0 off) and these fields `(key, label, unit, default)`, edited with
[`edit_part!`](@ref). On, each entrance hole is open out to `leh_open_fraction ×` the
LEH radius and its lip beyond closes to x-rays at `leh_closure_time` over
`leh_closure_width` (a `GapClosure` whose closed albedo is the wall's). 0.84 is the
measured LEH radius at peak (research notes, "Closing entrance holes"); the time and
width are an engineering example (closure at the rise to peak), to be set per design.
"""
const LEH_CLOSURE_FIELDS = ((:leh_closure_time, "Closure time", "ns", 5.0),
                            (:leh_closure_width, "Closure duration", "ns", 0.5),
                            (:leh_open_fraction, "Open radius", "× LEH radius", 0.84))

# The closing lip of an entrance hole (`:hohlraum_leh_top_lip`, …) is a split of the
# hole for the calculations: the GUI lists, selects, highlights and leaves it out with
# its hole.
const _LIP_NAME = r"^((?:[A-Za-z0-9_]+_)?leh_(?:top|bottom))_lip$"
"Whether `c` is a hohlraum part's closing entrance-hole lip."
is_lip(c::Component) =
    c.kind === :aperture && c.material.albedo isa GapClosure && occursin(_LIP_NAME, String(c.name))
"The component a name is shown as: an entrance hole for its lip, else itself."
shown_component(name::Symbol) = (m = match(_LIP_NAME, String(name)); m === nothing ? name : Symbol(m[1]))
"Components `names` with the lips of the entrance holes among them."
with_lips(names) = Symbol[n for x in names
                          for n in (occursin(r"(^|_)leh_(top|bottom)$", String(x)) ? (x, Symbol(x, :_lip)) : (x,))]

# ---------------------------------------------------------------- windows

"""
    WindowSpec(; phi = 0.0, z = 0.0, width = 1.0, height = 1.0, open = false,
               albedo = 0.85, conversion = 0.8, round = true)

A diagnostic window cut into a hohlraum part's wall (a core `WallWindow`): centered at
azimuth `phi` (degrees from the part's +x) and axial position `z` (mm, in the part's
frame); `round` (the default, as diagnostic windows usually are) a disk of diameter
`width` mm on the unrolled wall, otherwise a rectangle `width` mm of arc by `height` mm
along the axis. A closed window (`open = false`) is an opaque wall patch with its own
x-ray `albedo` and laser → x-ray `conversion`; an open one is a hole (an aperture: laser
light and x-rays leave through it). Window `k` of part `name` is the component
`<name>_window_k`.
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
           albedo = 0.85, conversion = 0.8, round::Bool = true) =
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
"""
const PART_KINDS = (:sphere, :disk, :annulus, :rectangle, :cylinder, :cone, :mesh, :visrad, :hohlraum, :frame)

const PART_KIND_LABELS = Dict(:sphere => "Sphere", :disk => "Disk", :annulus => "Washer (annulus)",
                              :rectangle => "Rectangle", :cylinder => "Cylinder (tube)",
                              :cone => "Cone (frustum)", :mesh => "Mesh (STL file)",
                              :visrad => "VisRad object", :hohlraum => "Hohlraum", :frame => "Frame")

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
    :hohlraum => ((:length, "Length", "mm", 10.1), (:radius, "Radius", "mm", 2.875),
                  (:leh_radius, "LEH radius", "mm", 1.55)))

"""
Values of a `:hohlraum` part besides its sizes (`key => default`): its wall's albedo and
laser → x-ray conversion (the placement fields `albedo` and `conversion`, which default
to 0.85 and 0.8 for a hohlraum, as the former hohlraum template had), each end's opening (`top_end`,
`bottom_end`: 1 an entrance hole of radius `leh_radius`, 0 closed: the cap covers the
end), and the closing entrance holes: `leh_closure` (1 on, 0 off) with the
[`LEH_CLOSURE_FIELDS`](@ref).
"""
const HOHLRAUM_VALUES = (:albedo => 0.85, :conversion => 0.8, :top_end => 1.0, :bottom_end => 1.0,
                         :leh_closure => 0.0,
                         (k => d for (k, _, _, d) in LEH_CLOSURE_FIELDS)...)

"""
Values of a `:sphere` part besides its radius (`key => default`): `capsule_mesh` (1:
meshed as a capsule, the geodesic mesh of `Briareus.hohlraum`'s capsule with elements of
the target's capsule size, `MeshSpec.capsule_h`, unless the part sets its own `h`; 0,
the default: the latitude–longitude `sphere_mesh` with elements of `h`). The starter
hohlraum's capsule sets it, so it meshes as the former hohlraum template's capsule did.
"""
const SPHERE_VALUES = (:capsule_mesh => 0.0,)

# The values a kind has besides the placement and size fields.
_kind_values(kind::Symbol) = kind === :hohlraum ? HOHLRAUM_VALUES : kind === :sphere ? SPHERE_VALUES : ()

"""
    hohlraum_pieces(p::PartSpec) -> Vector{Symbol}

The piece names of a `:hohlraum` part `p` without its prefix, in the order of its
components (each component is named `<part>_<piece>` and carries `part = <part>`):
`:wall`, `:cap_top`, `:cap_bottom`, then for each end with an entrance hole `:leh_top`
or `:leh_bottom` (and its closing lip `:leh_top_lip`, … with `leh_closure`), then the
windows `:window_1`, ….
"""
function hohlraum_pieces(p)                     # p::PartSpec (defined below)
    p.kind === :hohlraum || throw(ArgumentError("$(p.name) is not a hohlraum"))
    v = p.values
    out = [:wall, :cap_top, :cap_bottom]
    for (e, leh) in ((:top_end, :leh_top), (:bottom_end, :leh_bottom))
        v[e] != 0 || continue
        push!(out, leh)
        v[:leh_closure] != 0 && push!(out, Symbol(leh, :_lip))
    end
    append!(out, (window_name(k) for k in eachindex(p.windows)))
    return out
end

"The component name of piece `piece` of hohlraum part `name`: `<name>_<piece>`."
piece_name(name::Symbol, piece::Symbol) = Symbol(name, :_, piece)

"""
Placement and material fields of every part: `(key, label, unit, default)`. The part
is rotated by `rotation` about its own axis, then tilted to polar angle `theta` and
azimuth `phi` (VisRad's Polar/Azimuthal/Rotation order: the axis then points along
`(θ, φ)`), then moved to `(x, y, z)`. `h` (mm) overrides the mesh element size (0:
the target's; for a VisRad primitive, VisRad's own grid, see [`part_size_unit`](@ref)).
"""
const PART_FIELDS = ((:x, "x", "mm", 0.0), (:y, "y", "mm", 0.0), (:z, "z", "mm", 0.0),
                     (:theta, "Axis polar θ", "°", 0.0), (:phi, "Axis azimuth φ", "°", 0.0),
                     (:rotation, "Rotation about axis", "°", 0.0),
                     (:albedo, "Albedo", "", 0.0), (:conversion, "Conversion", "laser → x-ray", 0.0),
                     (:h, "Element size", "mm (0: target's)", 0.0))

"""
Labels of the named frames that are not a part's: `:tcc`, the target chamber (VisRad's
reference frame), the base every part is placed in by default and beams aim in. Every
part has its own frame too, named after the part, and a hohlraum part its pieces'
([`part_frames`](@ref)): a part placed in one moves and turns with it.
"""
const FRAME_LABELS = Dict(:tcc => "Target chamber (TCC)")

"""
    PartSpec(kind; name, values, two_sided = false, inward = false, aperture = false,
             file = "", frame = :tcc, coords = :cartesian, object = "", steps = [],
             foil = false, observed = false, wall = false, windows = [])

A free target part: `kind` (one of [`PART_KINDS`](@ref)), its component `name`,
`values` (every [`PART_FIELDS`](@ref) key and its kind's [`PART_SIZE_FIELDS`](@ref)),
`two_sided` (both sides radiate and absorb), `inward` (a sphere, cylinder or cone
radiates from its inside; a flat part from its back; a mesh from the other side) and
`aperture` (a transparent diagnostic surface: rays pass, its incident flux is what
escapes through it). The role flags mark the part for the analyses
([`set_part_role!`](@ref)): `foil` a backlighter foil ([`active_foil`](@ref)),
`observed` the surface whose x-ray drive is analyzed ([`observed_name`](@ref)), `wall`
a re-emitting cavity wall ([`wall_names`](@ref)); the name carries no role. A `:mesh`
part reads its triangles from the STL `file`. The
placement fields are in the frame named `frame` (`:tcc`, another part's name or a piece
frame of a hohlraum part: see [`part_frames`](@ref)); `coords` ([`PART_COORDINATES`](@ref)) is
how its position is shown and edited (the values keep `x`, `y`, `z`). The placement
`steps` (`PlacementStep`s: relations to other frames, relative moves) apply in order
after the placement fields; the degrees of freedom they leave free keep the fields'
values ([`part_placement`](@ref)). A `:hohlraum` part has its [`HOHLRAUM_VALUES`](@ref)
too and `windows`, its diagnostic windows ([`WindowSpec`](@ref)s; other kinds have none).
A `:frame` part is a user frame: its placement fields and steps define the frame and it
has no components. With `object` set it is a port frame ([`is_port_frame`](@ref)): the
frame of the port named `object`, at TCC with `z` toward the port (`theta`, `phi` the
port's angles, `rotation` its roll), placed in `:tcc` without steps.
"""
struct PartSpec
    kind::Symbol
    name::Symbol
    values::Dict{Symbol,Float64}
    two_sided::Bool
    inward::Bool
    aperture::Bool
    file::String
    frame::Symbol
    coords::Symbol
    object::String
    steps::Vector{PlacementStep}
    foil::Bool
    observed::Bool
    wall::Bool
    windows::Vector{WindowSpec}
end

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
                  observed::Bool = false, wall::Bool = false, windows = WindowSpec[])
    kind in PART_KINDS || throw(ArgumentError("part kind must be one of $PART_KINDS"))
    _own_frame_name(Symbol(name), frame) && throw(ArgumentError("a part cannot be placed in its own frame"))
    coords in PART_COORDINATES || throw(ArgumentError("coordinates must be one of $PART_COORDINATES"))
    kind === :mesh && isempty(file) && throw(ArgumentError("a mesh part needs an STL file"))
    kind === :visrad && (isempty(file) || isempty(object)) &&
        throw(ArgumentError("a VisRad part needs a workspace file and an object name"))
    v = Dict{Symbol,Float64}(k => d for (k, _, _, d) in (PART_FIELDS..., PART_SIZE_FIELDS[kind]...))
    for (k, d) in _kind_values(kind)
        v[k] = d
    end
    kind === :hohlraum || isempty(windows) || throw(ArgumentError("only a hohlraum has windows"))
    for (k, x) in values
        haskey(v, k) || throw(ArgumentError("a $kind has no field :$k"))
        v[k] = x
    end
    any(st -> Symbol(name) in Briareus._step_refs(st), steps) &&
        throw(ArgumentError("a part's placement step cannot refer to its own frame"))
    return _port_checked(PartSpec(kind, Symbol(name), v, two_sided, inward, aperture, String(file), frame,
                                  coords, String(object), collect(PlacementStep, steps), foil, observed, wall,
                                  collect(WindowSpec, windows)))
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

_with(p::PartSpec; kind = p.kind, name = p.name, values = p.values, two_sided = p.two_sided,
      inward = p.inward, aperture = p.aperture, file = p.file, frame = p.frame,
      coords = p.coords, object = p.object, steps = p.steps, foil = p.foil,
      observed = p.observed, wall = p.wall, windows = p.windows) =
    _port_checked(PartSpec(kind, name, values, two_sided, inward, aperture, file, frame, coords, object, steps,
                           foil, observed, wall, windows))

# A part with the frame named `old` renamed `new`: the frame it is placed in and the
# frames its placement steps refer to.
_renamed(p::PartSpec, old::Symbol, new::Symbol) =
    _with(p; frame = p.frame === old ? new : p.frame,
          steps = PlacementStep[_rename_step(st, old, new) for st in p.steps])

function _rename_step(st::PlacementStep, old::Symbol, new::Symbol)
    T = typeof(st)
    old in Briareus._step_refs(st) || return st
    return T((f === :ref ? new : getfield(st, f) for f in fieldnames(T))...)
end

"Axis (unit vector) of a part in its frame: polar `theta` and azimuth `phi` in degrees."
part_axis(p::PartSpec) = Vec3(sind(p.values[:theta]) * cosd(p.values[:phi]),
                              sind(p.values[:theta]) * sind(p.values[:phi]), cosd(p.values[:theta]))

"Center of a part in its frame."
part_center(p::PartSpec) = Vec3(p.values[:x], p.values[:y], p.values[:z])

# The part's own frame in the chamber before its placement steps, placed in its parent's
# chamber frame `frame`: the rotation `euler_rotation(theta, phi, rotation)` (VisRad's
# order) at its center.
_own_frame(p::PartSpec, frame::Frame = CHAMBER) =
    compose(frame, euler_frame(part_center(p), p.values[:theta], p.values[:phi], p.values[:rotation]))

"""
    part_frame_def(p::PartSpec) -> FrameDef

A part's own frame as placed in its parent `p.frame`: origin at its center
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

Every named frame of a design in the chamber: `:tcc`, each part's own
([`part_frame_def`](@ref), named after the part) and its pieces'
([`piece_frame_defs`](@ref)), resolved with `resolve_frames`, so a part placed in
another part's frame follows every ancestor exactly. Throws `ArgumentError` for a part
placed in a frame that does not exist, or parts placed in each other in a loop (the
error names the loop).
"""
function part_frames(parts::AbstractVector)
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

The child frames of a part's pieces, placed in the part's own frame: for a `:hohlraum`
part `name` of length `L`, the entrance holes' `name.leh_top` (origin at the top end's
center `(0, 0, L/2)`, the part's axes, so `z` points out of the hohlraum) and
`name.leh_bottom` (at `(0, 0, −L/2)`, turned 180° about the part's `x`: `z` out of the
hohlraum, `y` reversed, so one placement or aim mirrors between the two ends, as the
former template's entrance-hole frames did), and each window's `name.window_k` (origin at the
window's center on the wall, the part's axes). The ends' frames serve their caps too (a
cap's centroid is its end's center), and the wall aims in the part's own frame (its
centroid is the part's center). Other kinds have none.
"""
function piece_frame_defs(p::PartSpec)
    out = Pair{Symbol,FrameDef}[]
    p.kind === :hohlraum || return out
    L, R = p.values[:length], p.values[:radius]
    push!(out, Symbol(p.name, ".leh_top") => FrameDef(p.name; offset = (0.0, 0.0, L / 2)),
          Symbol(p.name, ".leh_bottom") => FrameDef(p.name; offset = (0.0, 0.0, -L / 2),
                                                    rotation = (180.0, 0.0, 180.0)))
    for (k, w) in enumerate(p.windows)
        push!(out, Symbol(p.name, ".", window_name(k)) =>
                   FrameDef(p.name; offset = (R * cosd(w.phi), R * sind(w.phi), w.z)))
    end
    return out
end

"""
    piece_frame(name::Symbol, part::Symbol) -> Symbol

The frame component `name` of part `part` aims in as a piece of it: a hohlraum part's
wall in the part's frame, its caps and entrance holes (and their lips) in its ends'
frames (`part.leh_top`, `part.leh_bottom`), a window in its own (`part.window_k`), see
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

"""
    part_component(p::PartSpec, h; frame) -> Component

The core component of a part, meshed with elements of size `h` (mm) unless the part
gives its own, carrying its exact curved surface (sphere, disk, washer, cylinder) for
the view factors, and placed in the chamber by `frame`, the part's own chamber `Frame`
(its entry in [`part_frames`](@ref), placement steps included). By default the part is
placed by its placement fields in the chamber frame, without its steps.
"""
function part_component(p::PartSpec, h::Real; frame::Frame = _own_frame(p))
    p.kind === :visrad && throw(ArgumentError("a VisRad object can make several components: use part_components"))
    p.kind === :frame && throw(ArgumentError("a frame has no component"))
    return only(part_components(p, h; frame))
end

"""
    part_components(p::PartSpec, h; frame, raw = false) -> Vector{Component}

The core components of a part, placed by its own chamber `frame` as in
[`part_component`](@ref): one for every kind but `:visrad`, whose object can be a
primitive hohlraum (wall, caps and entrance holes). A VisRad object keeps its own place
in the chamber (VisRad's frames), moved further by the part's placement; its windings
are repaired and its CAD surface remeshed into elements of the part's size. A VisRad
primitive (rectangle, sphere, cylindrical hohlraum) keeps VisRad's grid, which
reproduces VisRad's element table, unless the part sets its own element size. With `raw`
a VisRad CAD part keeps its file's facets as they are (no winding repair, no remeshing:
for drawing a part left out of the calculations, at no meshing cost).

A sphere with `capsule_mesh` set (value 1, see [`SPHERE_VALUES`](@ref)) is meshed as a
capsule: the geodesic mesh of `Briareus.hohlraum`'s capsule (`icosphere_mesh`,
`capsule_mesh = :compact`) with elements of `capsule_h` (the target's capsule size,
[`MeshSpec`](@ref)) unless the part gives its own. A `:hohlraum` part is core `hohlraum(...)` without a capsule, meshed with
elements of `h`, its components named and grouped as in [`hohlraum_pieces`](@ref); at
the chamber frame it is not moved at all, so it equals the former template's hohlraum bit for
bit.
"""
function part_components(p::PartSpec, h::Real; frame::Frame = _own_frame(p), raw::Bool = false,
                         capsule_h::Real = h)
    # the same components (identical objects) for the same part, size and frame, so the
    # targets built from them are found again (`_target_of`): an edit elsewhere, an undo
    lock(() -> length(_PART_MEMO) > 512 && empty!(_PART_MEMO), _STL_LOCK)
    hc = _capsule_meshed(p) ? Float64(capsule_h) : 0.0
    return _cached(() -> _part_components(p, h, frame, raw, hc), _PART_MEMO, _part_key(p, h, frame, raw, hc))
end

const _PART_MEMO = Dict{Any,Vector{Component}}()

# Whether a part is meshed as the capsule (`capsule_h`, geodesic).
_capsule_meshed(p::PartSpec) = p.kind === :sphere && p.values[:capsule_mesh] != 0

function _part_key(p::PartSpec, h::Real, frame, raw::Bool, hc::Real = 0.0)
    stamp = isfile(p.file) ? (mtime(p.file), filesize(p.file)) : nothing
    return (p.kind, p.name, sort!(collect(p.values)), p.two_sided, p.inward, p.aperture, p.file,
            p.object, stamp, Float64(h), frame.R, frame.origin, raw, Float64(hc), p.windows)
end

function _part_components(p::PartSpec, h::Real, frame::Frame, raw::Bool, hc::Real = 0.0)
    p.kind === :frame && return Component[]               # a user frame: no geometry
    R, center = frame.R, frame.origin
    v = p.values
    hp = v[:h] > 0 ? v[:h] : Float64(h)
    hp > 0 || throw(ArgumentError("$(p.name): the element size must be positive"))
    for (k, label, _, _) in PART_SIZE_FIELDS[p.kind]
        v[k] > 0 || throw(ArgumentError("$(p.name): the $(lowercase(label)) must be positive"))
    end
    0 <= v[:albedo] <= 1 || throw(ArgumentError("$(p.name): the albedo must be between 0 and 1"))
    0 <= v[:conversion] <= 1 ||
        throw(ArgumentError("$(p.name): the conversion must be between 0 and 1"))
    p.kind === :hohlraum && return _hohlraum_components(p, hp, frame)
    out = !p.inward
    k = p.kind
    if k === :sphere && hc > 0
        r = v[:r1]
        m = icosphere_mesh(r; n = icosphere_frequency(r, v[:h] > 0 ? v[:h] : Float64(hc)), outward = out)
        surf = SphereSurface(r)
    elseif k === :sphere
        r = v[:r1]
        (; capsule_ntheta, capsule_nphi) = capsule_counts(r, hp)
        m = sphere_mesh(r; ntheta = capsule_ntheta, nphi = capsule_nphi, outward = out)
        surf = SphereSurface(r)
    elseif k === :disk
        r = v[:r1]
        m = disk_mesh(r, 0.0; nr = _count(r, hp), nphi = _count(2π * r, hp; min = 8, multiple = 4), up = out)
        surf = DiskSurface()
    elseif k === :annulus
        ri, ro = v[:r1], v[:r2]
        ro > ri || throw(ArgumentError("$(p.name): the outer radius must exceed the inner"))
        m = annulus_mesh(ri, ro, 0.0; nr = _count(ro - ri, hp), nphi = _count(2π * ro, hp; min = 8, multiple = 4),
                         up = out)
        surf = DiskSurface()
    elseif k === :rectangle
        w, ht = v[:r1], v[:r2]
        m = rectangle_mesh(w, ht; nu = _count(w, hp), nv = _count(ht, hp))
        out || (m = transform(m, [1 0 0; 0 -1 0; 0 0 -1]))          # face −z
        surf = nothing
    elseif k === :visrad
        cs = _visrad_components(p, hp, transpose(R) * (zero(Vec3) - center); raw)
        return [transform(c, R, center) for c in cs]
    elseif k === :mesh
        v[:scale] > 0 || throw(ArgumentError("$(p.name): the scale must be positive"))
        # facing the target center, remeshed into elements of the part's size
        toward = transpose(R) * (zero(Vec3) - center)
        m = _stl_mesh(p, hp, toward; raw)
        # a part that is one cylinder or sphere (a tube, a shell) gets it for the view factors
        surf = raw ? nothing : _whole_surface(m, _stl_mesh(p, hp, toward; raw = true))
    elseif k === :cylinder
        r, L = v[:r1], v[:length]
        m = cylinder_mesh(r, -L / 2, L / 2; nz = _count(L, hp), nphi = _count(2π * r, hp; min = 8, multiple = 4),
                          inward = !out)
        surf = CylinderSurface(r)
    else
        r1, r2, L = v[:r1], v[:r2], v[:length]
        m = frustum_mesh(r1, -L / 2, r2, L / 2; nz = _count(L, hp),
                         nphi = _count(2π * max(r1, r2), hp; min = 8, multiple = 4), inward = !out)
        surf = nothing
    end
    c = p.aperture ? aperture(p.name, m; surface = surf) :
        Component(p.name, m, Material(p.name; albedo = v[:albedo], conversion = v[:conversion]);
                  two_sided = p.two_sided, surface = surf)
    return [transform(c, R, center)]
end

# The components of a `:hohlraum` part: core `hohlraum` (gold wall, no capsule) with
# elements of `h`, each component renamed `<part>_<name>` with `part = <part>`; a closed
# end's entrance hole (and lip) merged into its cap, of the wall's material. Placed by
# `frame`; at the chamber frame the meshes are not moved.
function _hohlraum_components(p::PartSpec, h::Real, frame::Frame)
    v = p.values
    for e in (:top_end, :bottom_end)
        v[e] in (0.0, 1.0) || throw(ArgumentError("$(p.name): an end is an entrance hole (1) or closed (0)"))
    end
    wall = gold(albedo = v[:albedo], conversion = v[:conversion])
    closing = v[:leh_closure] != 0
    lip = closing ? _leh_closure(v, wall) : nothing
    fopen = closing ? v[:leh_open_fraction] : 0.84
    t = try
        hohlraum(; length = v[:length], radius = v[:radius], leh_radius = v[:leh_radius],
                 capsule_radius = 0.0, wall,
                 windows = [_wall_window(w, k) for (k, w) in enumerate(p.windows)], h,
                 leh_closure = lip, leh_open_radius = fopen * v[:leh_radius])
    catch err
        err isa ArgumentError || rethrow()
        throw(ArgumentError("$(p.name): " * err.msg))
    end
    comps = copy(t.components)
    for (e, cap, leh) in ((:top_end, :cap_top, :leh_top), (:bottom_end, :cap_bottom, :leh_bottom))
        v[e] == 0 || continue
        hole = filter(c -> c.name in (leh, Symbol(leh, :_lip)), comps)
        k = findfirst(c -> c.name === cap, comps)
        c = comps[k]
        comps[k] = Component(cap, merge_meshes(c.mesh, (d.mesh for d in hole)...), c.material;
                             surface = c.surface)
        filter!(c -> !(c in hole), comps)
    end
    out = Component[Component(piece_name(p.name, c.name), c.mesh, c.material; kind = c.kind,
                              two_sided = c.two_sided, surface = c.surface, part = p.name)
                    for c in comps]
    frame == CHAMBER && return out
    return Component[transform(c, frame.R, frame.origin) for c in out]
end

# A VisRad object's components in the chamber frame (windings repaired toward `toward`,
# CAD surfaces remeshed into elements of `h`; with `raw` the file's facets as they are;
# primitives on VisRad's grid unless the part sets its own element size), kept per
# workspace version and settings: every edit rebuilds the target, and remeshing a CAD
# part of 10⁵ facets takes tens of seconds.
const _VRW_PARTS = Dict{Any,Vector{Component}}()

# `cache[key]`, computed by `f()` outside the lock when missing (parts mesh in parallel,
# and a long meshing never holds up the GUI thread's lookups of other parts). A key
# being computed is computed once: other callers wait for it (an edit made while a
# part remeshes in the background asks for the same part again).
const _INFLIGHT = Dict{Any,Base.Event}()

function _cached(f, cache::AbstractDict, key)
    ev, mine = lock(_STL_LOCK) do
        haskey(cache, key) && return (nothing, false)
        e = get(_INFLIGHT, key, nothing)
        e === nothing || return (e, false)
        e = _INFLIGHT[key] = Base.Event()
        return (e, true)
    end
    ev === nothing && return lock(() -> cache[key], _STL_LOCK)
    if !mine
        wait(ev)
        v = lock(() -> get(cache, key, nothing), _STL_LOCK)
        return v === nothing ? _cached(f, cache, key) : v      # it failed there: try here
    end
    try
        v = f()
        lock(() -> (cache[key] = v), _STL_LOCK)
        return v
    finally
        lock(() -> delete!(_INFLIGHT, key), _STL_LOCK)
        notify(ev)
    end
end

_iscached(cache::AbstractDict, key) = lock(() -> haskey(cache, key), _STL_LOCK)

# The cache key of a VisRad part's components, the element size a primitive takes
# (`nothing`: VisRad's grid, which reproduces VisRad's element table) and whether the
# part is drawn from its file's facets (`raw`, CAD only).
function _visrad_key(p::PartSpec, o::VisRadObject, h, toward, raw::Bool)
    v = p.values
    cad = o.kind === :mesh
    ph = cad || v[:h] <= 0 ? nothing : v[:h]
    raw &= cad
    key = (abspath(p.file), mtime(p.file), p.object, raw ? 0.0 : cad ? Float64(h) : something(ph, 0.0),
           raw || !cad ? nothing : round.(Tuple(toward); digits = 9), p.two_sided, v[:albedo],
           v[:conversion], p.name, raw)
    return key, ph, raw
end

function _visrad_components(p::PartSpec, h, toward; raw::Bool = false)
    v = p.values
    o = _visrad_object(p.file, p.object)
    cad = o.kind === :mesh
    key, ph, raw = _visrad_key(p, o, h, toward, raw)
    _cached(_VRW_PARTS, key) do
        mat = Material(p.name; albedo = v[:albedo], conversion = v[:conversion])
        # one component per CAD part (`fit = false`): the window shows and keys the parts;
        # the drive fits their curved faces in its own target (`_fitted_target`)
        cs = raw ? visrad_component(o; material = mat, two_sided = p.two_sided) :
             visrad_component(o; material = mat, two_sided = p.two_sided, orient = !p.two_sided,
                              toward, remesh = cad ? h : nothing, primitive_h = ph, fit = false)
        isempty(cs) && throw(ArgumentError("$(p.name): VisRad object type $(o.type) is not supported"))
        # a CAD part that is one cylinder or sphere (a tube, a shell) gets it for the view
        # factors, fitted on its file's facets
        if cad && !raw
            ref = visrad_component(o; material = mat, two_sided = p.two_sided)
            cs = [length(ref) == 1 ? _with_surface(c, _whole_surface(c.mesh, only(ref).mesh)) : c for c in cs]
        end
        # components carry the part's name (renamed, or given a role), not the object's
        String(p.name) == o.name ? cs :
        [Component(String(p.name) * String(c.name)[ncodeunits(o.name)+1:end], c.mesh, c.material;
                   kind = c.kind, two_sided = c.two_sided, surface = c.surface) for c in cs]
    end
end

"""
    part_size_unit(p::PartSpec) -> String

The unit text of a part's element-size field, with what 0 means for it (the target's
size, or for a VisRad primitive VisRad's own grid).
"""
function part_size_unit(p::PartSpec)
    p.kind === :visrad || return "mm (0: target's)"
    cad = try
        _visrad_object(p.file, p.object).kind === :mesh
    catch err
        err isa ArgumentError || rethrow()
        true
    end
    return cad ? "mm (0: target's)" : "mm (0: VisRad's grid)"
end

"Names of the components a VisRad object makes (without meshing its CAD)."
visrad_component_names(o::VisRadObject) =
    o.kind === :mesh ? [Symbol(o.name)] : [c.name for c in visrad_component(o)]

# A VisRad workspace's object, the workspace read once per file version.
const _VRW_CACHE = Dict{Any,VisRadWorkspace}()

function _visrad_workspace(file::AbstractString)
    isfile(file) || throw(ArgumentError("VisRad workspace not found: $file"))
    key = (abspath(file), mtime(file), filesize(file))
    return _cached(() -> read_visrad_workspace(file), _VRW_CACHE, key)
end

function _visrad_object(file::AbstractString, object::AbstractString)
    ws = _visrad_workspace(file)
    i = findfirst(o -> o.name == object, ws.objects)
    i === nothing && throw(ArgumentError("no object $object in $(basename(file))"))
    return ws.objects[i]
end

# An STL part's mesh: its file's facets (a whole file, or solid `object` of a
# multi-solid file) scaled to mm, oriented (`orient_mesh`: closed solids face out of
# their material, open sheets toward `toward`, flipped for `inward`) and remeshed into
# elements of `h`; with `raw` the facets as they are (drawing only). Every edit rebuilds
# the target, so the result is kept per file version and settings.
const _STL_CACHE = Dict{Any,Mesh}()
const _STL_LOCK = ReentrantLock()

function _stl_key(p::PartSpec, h, toward, raw::Bool)
    isfile(p.file) || throw(ArgumentError("STL file not found: $(p.file)"))
    return (abspath(p.file), mtime(p.file), filesize(p.file), p.object, p.values[:scale],
            raw ? 0.0 : Float64(h), raw ? nothing : round.(Tuple(toward); digits = 9),
            !raw && p.inward)
end

function _stl_mesh(p::PartSpec, h, toward; raw::Bool = false)
    scale = p.values[:scale]
    _cached(_STL_CACHE, _stl_key(p, h, toward, raw)) do
        m0 = if isempty(p.object)
            read_stl(p.file; scale)
        else
            solids = _stl_solids(p.file)
            k = parse(Int, p.object)
            1 <= k <= length(solids) || throw(ArgumentError("$(basename(p.file)) has no solid $k"))
            m = last(solids[k])
            Mesh([scale * v for v in m.vertices], m.triangles, m.patches)
        end
        raw && return m0
        m = remesh_mesh(orient_mesh(m0; toward); size = h)
        p.inward ? Mesh(m.vertices, [(a, c, b) for (a, b, c) in m.triangles], m.patches) : m
    end
end

"""
    _whole_surface(mesh, reference) -> Union{Nothing,CylinderSurface,SphereSurface}

The cylinder or sphere a remeshed CAD part lies on when the whole part is one
(`fit_surfaces` on its file's facets, `reference`, whose vertices lie on the true
surface; every element of `mesh` in that one region), else `nothing`. The part then
keeps its name and its elements and is integrated on the true surface: a tube or shell
tessellated in slivers otherwise carries hundreds of facets per element, and its view
factors cost ~20× more on the flat path (research notes, "Finely faceted curved CAD on
the flat path"). A part of several faces (a hohlraum half: its bore, walls, ends) stays
one flat component in the window, which shows and keys parts; the drive splits it into
its fitted pieces in its own target (`_fitted_target`).
"""
function _whole_surface(m::Mesh, ref::Mesh)
    npatches(m) >= 6 && length(ref.triangles) >= 6 || return nothing
    # `fit_surfaces` grows every face of a part (minutes for a 4 × 10⁴-facet CAD capsule
    # with its fill tube): only a part whose vertices are all on one sphere or cylinder
    # is worth it
    _on_one_surface(ref) || return nothing
    regs = try
        fit_surfaces(m; reference = ref, kinds = (:cylinder, :sphere))
    catch err
        err isa ArgumentError || rethrow()
        return nothing
    end
    length(regs) == 1 && length(regs[1].patches) == npatches(m) || return nothing
    s = regs[1].surface
    # a nearly flat part fits a sphere or cylinder of any large radius; visibility ignores
    # a curved part's facets between points inside it, so only a true tube or shell counts
    lo, hi = extrema(v -> v[1], ref.vertices), extrema(v -> v[2], ref.vertices)
    zl = extrema(v -> v[3], ref.vertices)
    extent = hypot(lo[2] - lo[1], hi[2] - hi[1], zl[2] - zl[1])
    return s !== nothing && s.radius <= extent ? s : nothing
end

# Whether every vertex of `m` lies on one sphere or one cylinder, within twice
# `fit_surfaces`' default tolerance (max(0.1 µm, 10⁻³ R)): algebraic least-squares fits to
# all vertices, the cylinder's axis the direction the facet normals are all
# perpendicular to. Milliseconds for 10⁵ vertices.
function _on_one_surface(m::Mesh)
    P = unique(m.vertices)
    length(P) >= 6 || return false
    tol(r) = 2 * max(1e-4, 1e-3 * r)
    # sphere: |p|² = 2 c·p + (r² − |c|²)
    A = [2p[k] for p in P, k in 1:3]
    x = [A ones(length(P))] \ [dot(p, p) for p in P]
    c = Vec3(x[1], x[2], x[3])
    r = sqrt(max(x[4] + dot(c, c), 0.0))
    r > 0 && maximum(p -> abs(norm(p - c) - r), P) <= tol(r) && return true
    # cylinder: the facet normals span the plane across the axis
    M = zeros(3, 3)
    for (a, b, d) in m.triangles
        n = cross(m.vertices[b] - m.vertices[a], m.vertices[d] - m.vertices[a])   # area-weighted
        M .+= n * transpose(n) ./ max(norm(n), eps())
    end
    F = eigen(Symmetric(M))
    F.values[1] <= 1e-3 * F.values[3] || return false
    ax = Vec3(F.vectors[:, 1]...)
    e1 = Vec3(F.vectors[:, 2]...)
    e2 = Vec3(F.vectors[:, 3]...)
    U = [dot(p, e1) for p in P]
    V = [dot(p, e2) for p in P]
    y = [2U 2V ones(length(P))] \ (U .^ 2 .+ V .^ 2)
    rc = sqrt(max(y[3] + y[1]^2 + y[2]^2, 0.0))
    return rc > 0 && maximum(k -> abs(hypot(U[k] - y[1], V[k] - y[2]) - rc), eachindex(P)) <= tol(rc)
end

_with_surface(c::Component, surf) =
    surf === nothing ? c : Component(c.name, c.mesh, c.material; kind = c.kind, two_sided = c.two_sided,
                                     surface = surf)

# The solids of a multi-solid (ASCII) STL file, read once per version.
const _STL_SOLIDS = Dict{Any,Vector{Pair{String,Mesh}}}()

function _stl_solids(file::AbstractString)
    isfile(file) || throw(ArgumentError("STL file not found: $file"))
    key = (abspath(file), mtime(file), filesize(file))
    return _cached(() -> read_stl_solids(file), _STL_SOLIDS, key)
end

# Whether part `p` (element size `h`, placed by its own chamber `frame`) is meshed
# already: CAD parts (VisRad CAD, STL) remesh for seconds to minutes and are cached; the
# other kinds mesh at once. A part the geometry rejects counts as ready (the build says why).
function _part_ready(p::PartSpec, h::Real, frame::Frame)
    p.kind in (:visrad, :mesh) || return true
    _iscached(_PART_MEMO, _part_key(p, h, frame, false)) && return true
    try
        R, center = frame.R, frame.origin
        hp = p.values[:h] > 0 ? p.values[:h] : Float64(h)
        toward = transpose(R) * (zero(Vec3) - center)
        if p.kind === :visrad
            o = _visrad_object(p.file, p.object)
            o.kind === :mesh || return true
            return _iscached(_VRW_PARTS, first(_visrad_key(p, o, hp, toward, false)))
        end
        return _iscached(_STL_CACHE, _stl_key(p, hp, toward, false))
    catch err
        err isa ArgumentError || rethrow()
        return true
    end
end

# Names of the CAD parts of a design that are in the calculations and not meshed yet.
function _unmeshed_parts(mesh::MeshSpec, parts, excluded)
    out = Symbol[]
    any(p -> p.kind in (:visrad, :mesh) && !(p.name in excluded), parts) || return out
    frames = try
        part_frames(parts)
    catch err
        err isa ArgumentError || rethrow()
        return out                                    # (the build says why)
    end
    for p in parts
        (p.kind in (:visrad, :mesh) && !(p.name in excluded)) || continue
        _part_ready(p, mesh.h, frames[p.name]) || push!(out, p.name)
    end
    return out
end

# Mesh `parts` (of the design's parts `all`) on worker threads, each once (the part
# caches then serve the target build, its commit and later edits); `raw` parts (names)
# keep their file's facets. Progress goes to `report`.
function _mesh_parts(parts, mesh, raw, report, label; all = parts)
    n = length(parts)
    finished = Threads.Atomic{Int}(0)
    report("$label: meshing $n parts…")
    frames = part_frames(all)
    tasks = map(parts) do p
        Threads.@spawn begin
            part_components(p, mesh.h; frame = frames[p.name], raw = p.name in raw,
                            capsule_h = mesh.capsule_h)
            k = Threads.atomic_add!(finished, 1) + 1
            k < n && report("$label: meshed $k of $n parts…")
        end
    end
    for t in tasks
        try
            fetch(t)
        catch err
            err isa TaskFailedException ? throw(err.task.exception) : rethrow()
        end
    end
    return nothing
end

"""
    add_stl_parts!(ts::TargetSetup, file; values...) -> Bool

Add the parts of an STL file: one `:mesh` part for a file of one solid (as
[`add_part!`](@ref)), or one per named solid of an ASCII file of several (a multi-body CAD
export, `read_stl_solids`), each named after its solid (letters, digits and `_`, made
unique) and placed and sized by `values` alike. The parts are added together: one the
geometry rejects (an empty solid, …) adds none, sets `error` and returns `false`.
"""
function add_stl_parts!(ts, file::AbstractString; values...)      # ts::TargetSetup (defined below)
    solids = try
        _stl_solids(file)
    catch err
        err isa ArgumentError || rethrow()
        ts.error[] = _errmsg(err)
        return false
    end
    length(solids) <= 1 && return add_part!(ts, :mesh; file, values...)
    taken = union(Set(all_component_names(ts)), (p.name for p in ts.parts[]))
    parts = copy(ts.parts[])
    vals = Dict{Symbol,Float64}(k => Float64(v) for (k, v) in values)
    for (k, (sname, _)) in enumerate(solids)
        base = replace(strip(sname), r"[^A-Za-z0-9_]" => "_")
        (isempty(base) || !occursin(r"^[A-Za-z_]", base)) && (base = "solid_" * base)
        name, j = Symbol(base), 1
        while name in taken
            j += 1
            name = Symbol(base, :_, j)
        end
        push!(taken, name)
        push!(parts, PartSpec(:mesh; name, file, object = string(k), values = vals))
    end
    return _try_build!(ts; parts)
end

# ---------------------------------------------------------------- the setup

"""
    TargetSetup(; starter = :hohlraum, resolution = :standard, mesh = nothing,
                parts = PartSpec[])

The target being designed: `parts` (its parts, [`PartSpec`](@ref)s; when none are given,
the parts of `starter`, one of [`STARTERS`](@ref) placed at TCC
([`starter_parts`](@ref)), or with `starter = :none` a sphere at TCC: a target is never
empty), `mesh` (the element sizes, a [`MeshSpec`](@ref)), `resolution` (the preset the
sizes match, one of [`MESH_RESOLUTIONS`](@ref), or `:custom`), `frames` (every named
frame of the design in the chamber, [`part_frames`](@ref): updated with the parts,
before the target is built), `renames` (while an edit updates `frames`, the parts it
renamed, old name → new, so that beams aimed in their frames follow; empty otherwise),
the built `target`, `error` (the last rejected edit), `excluded` (component names left
out of the calculations, see [`set_excluded!`](@ref)), `ghost` (those components, for
drawing, or `nothing`) and `active_foil` (the foil chosen with
[`set_active_foil!`](@ref), or `nothing`; it follows a rename, and
[`active_foil`](@ref) gives the foil the analyses take).

CAD parts (VisRad CAD, STL) are remeshed when their element size, the target's, or
their placement changes, and when an undo or a session brings back a size not meshed
before: seconds to minutes for 10⁵ facets. With `spawn[]` set (a window shares its
session's) that runs on worker threads: the edit is committed at once with those parts
drawn from their file's facets, `pending` names them, and the target is built again
from the current design when they are ready (through `deliver[]`, on the GUI thread).
Without it the edit waits for the meshing.

Building the `Target` itself (its BVH and element tables) takes a second for 10⁵
triangles, so with `spawn[]` a target of more than [`BACKGROUND_TRIANGLES`](@ref)
triangles that has not been built before is built on a worker too: the edit's inputs
(`parts`, `mesh`, …) are committed at once, `target` and `ghost` keep the previous
design's until the new ones are ready, and `building` is set meanwhile. A CAD part
given the role of the capsule has its sphere fitted on a worker as well (`fitting`
names it; the role is committed with the fit). The drive and the checks wait for both
([`target_ready`](@ref)).

Edit with [`insert_starter!`](@ref), [`new_target!`](@ref), [`set_resolution!`](@ref),
[`edit_mesh!`](@ref), the window functions ([`add_window!`](@ref), …) and the part
functions ([`add_part!`](@ref), …).
"""
struct TargetSetup
    mesh::Observable{MeshSpec}
    resolution::Observable{Symbol}
    parts::Observable{Vector{PartSpec}}
    target::Observable{Target}
    error::Observable{String}
    excluded::Observable{Set{Symbol}}
    ghost::Observable{Union{Nothing,Target}}
    pending::Observable{Set{Symbol}}       # CAD parts being remeshed (drawn from their facets)
    spawn::Base.RefValue{Bool}             # remesh on workers (a window's session sets it)
    deliver::Base.RefValue{Any}            # hands a worker's result to the GUI thread
    mesh_job::Base.RefValue{Int}           # the latest remeshing or build job (earlier ones are dropped)
    building::Observable{Bool}             # the target of the committed inputs is being built
    fitting::Observable{Set{Symbol}}       # CAD parts whose capsule sphere is being fitted
    fit_jobs::Dict{Symbol,Int}             # their jobs (a fit whose job is gone is dropped)
    after_build::Vector{Any}               # run once the target being built is committed
    frames::Observable{Dict{Symbol,Frame}} # the named frames of the committed design (`part_frames`)
    renames::Base.RefValue{Dict{Symbol,Symbol}} # old => new part names of the edit committing `frames`
    active_foil::Observable{Union{Nothing,Symbol}} # the chosen foil (`active_foil` resolves it)
end

function TargetSetup(; starter::Symbol = :hohlraum, resolution::Symbol = :standard,
                     mesh::Union{Nothing,MeshSpec} = nothing,
                     parts::AbstractVector{PartSpec} = PartSpec[])
    if mesh === nothing
        resolution in MESH_RESOLUTIONS ||
            throw(ArgumentError("resolution must be one of $MESH_RESOLUTIONS"))
        mesh = MESH_PRESETS[resolution]
    end
    ps = isempty(parts) ? _starter_or_sphere(starter, mesh) : collect(PartSpec, parts)
    t = build_target(mesh, ps)
    return TargetSetup(Observable(mesh), Observable(_preset_of(mesh)), Observable(ps),
                       Observable(t), Observable(""), Observable(Set{Symbol}()),
                       Observable{Union{Nothing,Target}}(nothing), Observable(Set{Symbol}()),
                       Ref(false), Ref{Any}(f -> f()), Ref(0), Observable(false),
                       Observable(Set{Symbol}()), Dict{Symbol,Int}(), Any[],
                       Observable(part_frames(ps)), Ref(Dict{Symbol,Symbol}()),
                       Observable{Union{Nothing,Symbol}}(nothing))
end

# The parts of a new target: a starter's, or (`:none`) a sphere at TCC.
function _starter_or_sphere(starter::Symbol, mesh::MeshSpec)
    starter === :none && return [PartSpec(:sphere; name = :sphere_1)]
    starter in STARTERS || throw(ArgumentError("a starter is one of $STARTERS, or :none"))
    return starter_parts(starter; capsule_h = mesh.capsule_h)
end

part_frames(ts::TargetSetup) = ts.frames[]

"""
Targets of more triangles than this (≈ 30 ms to build) are built on a worker in a window
(see [`TargetSetup`](@ref)); smaller ones, such as the starters' meshes, at once. Set it
to change the threshold for a session (tests set 0).
"""
const BACKGROUND_TRIANGLES = Ref(30_000)

"""
    target_ready(ts::TargetSetup) -> Bool

Whether `target` is the target of the committed design: no CAD part is being remeshed
(`pending`), no target is being built (`building`) and no capsule sphere is being
fitted (`fitting`). The drive and the checks wait for it in a window.
"""
target_ready(ts::TargetSetup) = isempty(ts.pending[]) && !ts.building[] && isempty(ts.fitting[])

# What the target is waiting for, for status lines ("" when ready).
function target_waiting(ts::TargetSetup)
    out = String[]
    isempty(ts.pending[]) || push!(out, "remeshing " * join(sort!(String.(collect(ts.pending[]))), ", "))
    isempty(ts.fitting[]) ||
        push!(out, "fitting the capsule to " * join(sort!(String.(collect(ts.fitting[]))), ", "))
    ts.building[] && isempty(ts.pending[]) && push!(out, "building the target")
    return join(out, " · ")
end

# Run `f()` once the target of the committed design is in `target` (at once if it is).
function _after_build!(f, ts::TargetSetup)
    ts.building[] ? push!(ts.after_build, f) : f()
    return nothing
end

"""
    await_target!(ts::TargetSetup) -> Bool

Make `target` the committed design's now if a worker is still building it (the build
in flight is waited for, not repeated), as the checks that read it need. CAD parts
still remeshing stay `pending` (their facets are drawn) unless the session no longer
spawns jobs (a closed window), when they are meshed here. Returns whether the target
is current.
"""
function await_target!(ts::TargetSetup)
    (ts.building[] || (!ts.spawn[] && !isempty(ts.pending[]))) || return true
    return _try_build!(ts; sync = true)
end

"""
    build_target(mesh, parts) -> Target
    build_target(ts::TargetSetup) -> Target

The `Target` of `parts` ([`PartSpec`](@ref)s) meshed with the element sizes `mesh` (a
[`MeshSpec`](@ref) or a preset name). Throws `ArgumentError` for parts the geometry
rejects (overlapping windows, a window past the end of a hohlraum's wall, a duplicate
component name, …) and for no parts at all.
"""
build_target(resolution::Symbol, parts) = build_target(MESH_PRESETS[resolution], parts)

build_target(mesh::MeshSpec, parts::AbstractVector{PartSpec}; raw = ()) =
    _target_of(_target_components(mesh, parts; raw))
build_target(ts::TargetSetup) = build_target(ts.mesh[], ts.parts[])

# The components of a design (`raw`: names of CAD parts drawn from their file's facets).
function _target_components(mesh::MeshSpec, parts::AbstractVector{PartSpec}; raw = ())
    mesh.h > 0 && mesh.capsule_h > 0 || throw(ArgumentError("element sizes must be positive"))
    comps = Component[]
    frames = isempty(parts) ? Dict{Symbol,Frame}() : part_frames(parts)
    for part in parts
        new = part_components(part, mesh.h; frame = frames[part.name], raw = part.name in raw,
                              capsule_h = mesh.capsule_h)
        for c in new
            any(d -> d.name === c.name, comps) && throw(ArgumentError("two components are named $(c.name)"))
        end
        append!(comps, new)
    end
    isempty(comps) && throw(ArgumentError("the target is empty: add a part"))
    return comps
end

# The lips' `GapClosure`: closed, they re-emit with the wall's albedo.
function _leh_closure(p::AbstractDict, wall::Material)
    p[:leh_closure_width] >= 0 || throw(ArgumentError("the closure duration must not be negative"))
    0 < p[:leh_open_fraction] < 1 ||
        throw(ArgumentError("the open radius must be between 0 and 1 of the LEH radius"))
    return GapClosure(closure = p[:leh_closure_time], width = p[:leh_closure_width], after = wall.albedo)
end

# The `Target` of `comps`, built once for the same components (identical objects): the
# few latest are kept, so an edit to a part left out of the calculations, or an undo,
# finds its targets again (a NIF CAD workspace's take a second to build: BVHs of
# 2 × 10⁵ facets). A target being built on a worker is waited for, not built twice.
# Keyed by the components' object ids: the entry holds the components, so no other
# object takes their ids while it is kept.
const _TARGETS = Dict{Vector{UInt},Pair{Vector{Component},Target}}()
const _TARGETS_ORDER = Vector{UInt}[]          # least recently used first
const _TARGETS_KEPT = 6                        # calculation targets and ghosts of 3 designs

_target_key(comps::Vector{Component}) = UInt[objectid(c) for c in comps]
_target_built(comps::Vector{Component}) = _iscached(_TARGETS, _target_key(comps))

function _target_of(comps::Vector{Component})
    key = _target_key(comps)
    e = _cached(() -> comps => Target(comps), _TARGETS, key)
    lock(_STL_LOCK) do
        filter!(!=(key), _TARGETS_ORDER)
        push!(_TARGETS_ORDER, key)
        while length(_TARGETS_ORDER) > _TARGETS_KEPT
            delete!(_TARGETS, popfirst!(_TARGETS_ORDER))
        end
    end
    return e.second
end

# A design's components split into the target of those in the calculations and that of
# those left out (`excluded`, drawn as a ghost: `nothing` when none is).
function _split_components(comps::Vector{Component}, excluded)
    keep, gh = _split_lists(comps, excluded)
    return _targets_of(keep, gh)
end

# The components in the calculations and those left out (`nothing` when none is).
function _split_lists(comps::Vector{Component}, excluded)
    # a lip with its hole; a piece with its part
    out(c) = (is_lip(c) ? shown_component(c.name) : c.name) in excluded || c.part in excluded
    keep = [c for c in comps if !out(c)]
    isempty(keep) && throw(ArgumentError("every component is left out of the calculations"))
    length(keep) == length(comps) && return keep, nothing
    return keep, [c for c in comps if out(c)]
end

_targets_of(keep, gh) = _target_of(keep), gh === nothing ? nothing : _target_of(gh)

# Whether building the targets of these lists is worth a worker: not built yet, and
# more than `BACKGROUND_TRIANGLES` triangles to build.
function _build_later(keep, gh)
    n = 0
    for cs in (keep, gh)
        (cs === nothing || _target_built(cs)) && continue
        n += sum(c -> length(c.mesh.triangles), cs; init = 0)
    end
    return n > 0 && n > BACKGROUND_TRIANGLES[]
end

# Rebuild with changed inputs, committing them only if the geometry accepts them.
# `prebuilt` is the `(target, ghost)` of exactly these inputs when a worker built it
# already ([`import_visrad_workspace!`](@ref)). With `spawn[]`, CAD parts not meshed yet
# are drawn from their file's facets and meshed on workers (`pending`), and the target
# is built again when they are ready; a large target not built before is built on a
# worker (`building`) unless `sync`, the inputs committed at once.
function _try_build!(ts::TargetSetup; mesh = ts.mesh[], parts = ts.parts[], excluded = ts.excluded[],
                     prebuilt = nothing, sync::Bool = !ts.spawn[], renames = Dict{Symbol,Symbol}())
    try
        frames = part_frames(parts)
        waiting = prebuilt === nothing && ts.spawn[] ? _unmeshed_parts(mesh, parts, excluded) : Symbol[]
        # components left out of the calculations (those still in the target)
        built, lists = if prebuilt === nothing
            keep, gh = _split_lists(_target_components(mesh, parts; raw = union(Set{Symbol}(excluded), waiting)),
                                    excluded)
            !sync && _build_later(keep, gh) ? (nothing, (keep, gh)) : (_targets_of(keep, gh), nothing)
        else
            prebuilt, nothing
        end
        _assign!(ts.excluded, Set{Symbol}(excluded))
        built === nothing || (ts.ghost[] === built[2] || (ts.ghost[] = built[2]))
        # (first: what follows a frame sees the design's frames;
        # beams aimed in a renamed part's frame follow its new name)
        ts.renames[] = Dict{Symbol,Symbol}(renames)
        try
            _assign!(ts.frames, frames)
        finally
            ts.renames[] = Dict{Symbol,Symbol}()
        end
        if mesh != ts.mesh[]
            ts.mesh[] = mesh
            _assign!(ts.resolution, _preset_of(mesh))
        end
        # the chosen foil keeps its role through a rename
        a = ts.active_foil[]
        a !== nothing && haskey(renames, a) && (ts.active_foil[] = renames[a])
        parts === ts.parts[] || (ts.parts[] = parts)
        g = (ts.mesh_job[] += 1)                  # a meshing or build job for earlier inputs is dropped
        ts.error[] = ""
        if built === nothing
            # the target and ghost of the previous inputs stay until these are built
            _assign!(ts.building, true)
            _assign!(ts.pending, Set{Symbol}(waiting))
            deliver = ts.deliver[]
            job = () -> _targets_of(lists...)
            _run_job(() -> _deliver_compiled(deliver, out -> _built!(ts, g, out), _try_job(job)), true)
        else
            built[1] === ts.target[] || (ts.target[] = built[1])
            _assign!(ts.pending, Set{Symbol}(waiting))
            _built_now!(ts)
        end
        isempty(waiting) || _mesh_in_background!(ts, g, waiting, mesh, parts)
        return true
    catch err
        err isa ArgumentError || rethrow()
        ts.error[] = _errmsg(err)
        return false
    end
end

# A worker built the target of job `g`: commit it unless a later edit superseded it.
function _built!(ts::TargetSetup, g::Integer, out)
    g == ts.mesh_job[] || return false
    if out.err !== nothing                       # (a bug: the inputs were checked)
        ts.error[] = "building the target failed: " * _errmsg(out.err)
        empty!(ts.after_build)
        _assign!(ts.building, false)
        return false
    end
    t, ghost = out.value
    ts.ghost[] === ghost || (ts.ghost[] = ghost)
    t === ts.target[] || (ts.target[] = t)
    _built_now!(ts)
    return true
end

# The target of the committed design is in place: run what waited for it.
function _built_now!(ts::TargetSetup)
    _assign!(ts.building, false)
    isempty(ts.after_build) && return
    fs = copy(ts.after_build)
    empty!(ts.after_build)
    foreach(f -> f(), fs)
    return
end

# Mesh the parts `names` on workers, then build the target again from the design as it
# is then (job `g`: a later edit's job supersedes it; a part both ask for is meshed once).
function _mesh_in_background!(ts::TargetSetup, g::Integer, names, mesh, parts)
    deliver = ts.deliver[]
    todo = [p for p in parts if p.name in names]
    job = () -> _mesh_parts(todo, mesh, Set{Symbol}(), _ -> nothing, ""; all = parts)
    _run_job(() -> _deliver_compiled(deliver, out -> _meshed!(ts, g, out), _try_job(job)), true)
    return nothing
end

function _meshed!(ts::TargetSetup, g::Integer, out)
    g == ts.mesh_job[] || return false
    if out.err !== nothing
        ts.error[] = "remeshing " * join(sort!(String.(collect(ts.pending[]))), ", ") * " failed: " *
                     _errmsg(out.err)
        return false
    end
    return _try_build!(ts)
end

"""
    set_excluded!(ts::TargetSetup, name::Symbol, excluded::Bool) -> Bool

Leave component `name` out of the calculations (`excluded = true`) or put it back.
An excluded component is not part of `target` (no laser, view factors, clearance or
checks see it) but stays listed and drawn as a ghost (`ghost`), so it can be put
back; a hohlraum with an excluded component is no longer treated as closed.
"""
function set_excluded!(ts::TargetSetup, name::Symbol, excluded::Bool)
    ex = copy(ts.excluded[])
    excluded ? push!(ex, name) : delete!(ex, name)
    ex == ts.excluded[] && return true
    return _try_build!(ts; excluded = ex)
end

"Every component name of the design: in the calculations, then excluded ones."
all_component_names(ts::TargetSetup) =
    [component_names(ts.target[]); ts.ghost[] === nothing ? Symbol[] : component_names(ts.ghost[])]

"""
    set_leh_closure!(ts::TargetSetup, on::Bool; part = nothing) -> Bool

Close hohlraum part `part`'s entrance holes (by default the design's first hohlraum
part's) during the pulse (`on`), with the settings of [`LEH_CLOSURE_FIELDS`](@ref), or
keep them open (its value `leh_closure`). Each hole is then split into the hole and its
closing lip (`<part>_leh_top_lip`, …), shown, selected and left out with the hole.
"""
function set_leh_closure!(ts::TargetSetup, on::Bool; part = nothing)
    try
        i = _window_part(ts, part)
        p = ts.parts[][i]
        (p.values[:leh_closure] != 0) == on && return true
        return edit_part!(ts, i, :leh_closure, on ? "1" : "0")
    catch err
        err isa ArgumentError || rethrow()
        ts.error[] = _errmsg(err)
        return false
    end
end

"""
    new_target!(ts::TargetSetup, starter = :hohlraum) -> Bool

Start a new target (File ▸ New target): the parts of `starter` alone (one of
[`STARTERS`](@ref), placed at TCC, [`starter_parts`](@ref)), or with `starter = :none`
a sphere at TCC; nothing is left out of the calculations, and the element sizes are
kept.
"""
function new_target!(ts::TargetSetup, starter::Symbol = :hohlraum)
    parts = try
        _starter_or_sphere(starter, ts.mesh[])
    catch err
        err isa ArgumentError || rethrow()
        ts.error[] = _errmsg(err)
        return false
    end
    _cancel_fits!(ts)
    _try_build!(ts; parts, excluded = Set{Symbol}()) || return false
    ts.active_foil[] = nothing
    return true
end

"""
    set_resolution!(ts::TargetSetup, resolution::Symbol) -> Bool

Mesh element sizes from a preset, one of [`MESH_RESOLUTIONS`](@ref) (see
[`MESH_PRESETS`](@ref)).
"""
function set_resolution!(ts::TargetSetup, resolution::Symbol)
    resolution in MESH_RESOLUTIONS ||
        throw(ArgumentError("resolution must be one of $MESH_RESOLUTIONS"))
    ts.mesh[] == MESH_PRESETS[resolution] && return true
    return _try_build!(ts; mesh = MESH_PRESETS[resolution])
end

"""
    edit_mesh!(ts::TargetSetup, key::Symbol, text) -> Bool

Set an element size (`:h` or `:capsule_h`, mm) from text; the resolution becomes the
preset the sizes match, or `:custom`.
"""
function edit_mesh!(ts::TargetSetup, key::Symbol, text)
    key in (:h, :capsule_h) || throw(ArgumentError("unknown mesh field :$key"))
    try
        v = _parse_number(text, key === :h ? "the element size" : "the capsule element size";
                          positive = true)
        m = ts.mesh[]
        m = key === :h ? MeshSpec(v, m.capsule_h) : MeshSpec(m.h, v)
        return _try_build!(ts; mesh = m)
    catch err
        err isa ArgumentError || rethrow()
        ts.error[] = _errmsg(err)
        return false
    end
end

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
# target and report why (without the part's name: the edit names the window).
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
                # the windows are checked as the wall is meshed (the hohlraum alone)
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
a closed wall patch with its own albedo and conversion (`false`).
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

# ---------------------------------------------------------------- parts

function _check_part_index(ts::TargetSetup, i::Integer)
    1 <= i <= length(ts.parts[]) || throw(ArgumentError("no part $i"))
    return i
end

"""
    add_part!(ts::TargetSetup, kind::Symbol; name, file = "", frame = :tcc, values...) -> Bool

Add a free part (one of [`PART_KINDS`](@ref)) named `name` (default: the kind and the
first free number) with [`PART_FIELDS`](@ref) and size values given as keywords
(e.g. `add_part!(ts, :disk; z = 6.0, r1 = 1.0, theta = 180)`), placed in the frame named
`frame` ([`part_frames`](@ref)). A `:mesh` part reads the
STL `file` (scaled by `scale`, mm per file unit), placed like any part, with its
windings repaired (closed solids facing out of their material, open sheets toward the
target center) and its surface remeshed into elements of the part's size (`h`, or the
target's; `remesh_mesh`). A part the geometry rejects sets `error` and returns `false`.
"""
function add_part!(ts::TargetSetup, kind::Symbol; name = nothing, two_sided::Bool = false,
                   inward::Bool = false, aperture::Bool = false, file::AbstractString = "",
                   frame::Symbol = :tcc, values...)
    # (the parts too: the target may still be the previous design's, being built)
    taken = union(Set(component_names(ts.target[])), (p.name for p in ts.parts[]))
    if name === nothing
        k = 1
        while Symbol(kind, :_, k) in taken
            k += 1
        end
        name = Symbol(kind, :_, k)
    end
    part = try
        PartSpec(kind; name, values = Dict{Symbol,Float64}(k => Float64(v) for (k, v) in values),
                 two_sided, inward, aperture, file, frame)
    catch err
        err isa ArgumentError || rethrow()
        ts.error[] = _errmsg(err)
        return false
    end
    return _try_build!(ts; parts = [ts.parts[]; part])
end

"""
    remove_part!(ts::TargetSetup, i) -> Bool
"""
function remove_part!(ts::TargetSetup, i::Integer)
    _check_part_index(ts, i)
    p = ts.parts[][i]
    if p.kind === :frame
        users = frame_users(ts, p.name)
        isempty(users) ||
            (ts.error[] = "frame $(p.name) is in use: $(join(String.(users), ", ")) placed in it or related to it";
             return false)
    end
    ps = copy(ts.parts[])
    deleteat!(ps, i)
    return _try_build!(ts; parts = ps)
end

# ---------------------------------------------------------------- user frames

"Whether part `p` is a port frame: a `:frame` part tied to a port (its `object`)."
is_port_frame(p::PartSpec) = p.kind === :frame && !isempty(p.object)

"""
    user_frames(ts::TargetSetup) -> Vector{PartSpec}

The user frames of a design (its `:frame` parts, port frames included), in order.
"""
user_frames(ts::TargetSetup) = PartSpec[p for p in ts.parts[] if p.kind === :frame]

"""
    frame_users(ts::TargetSetup, name::Symbol) -> Vector{Symbol}

The parts and frames placed in the frame `name` (or a piece frame of part `name`), or
with a placement step that refers to it. A user frame in use is not removed
([`remove_part!`](@ref)); beams aimed in a removed frame keep their aims, held in `:tcc`.
"""
function frame_users(ts::TargetSetup, name::Symbol)
    own(f) = _own_frame_name(name, f)
    return Symbol[q.name for q in ts.parts[] if q.name !== name &&
                  (own(q.frame) || any(st -> any(own, Briareus._step_refs(st)), q.steps))]
end

"""
    add_frame!(ts::TargetSetup; name = nothing, frame = :tcc, values...) -> Bool

Add a user frame (a `:frame` part, see [`PartSpec`](@ref)) named `name` (default
`frame_1`, …) placed in the frame `frame` by the placement fields given as keywords
(`x`, `y`, `z` in mm, `theta`, `phi`, `rotation` in degrees, as a part's). Parts can be
placed in it and beams aimed in it like any part's frame. A rejected frame sets `error`
and returns `false`.
"""
add_frame!(ts::TargetSetup; name = nothing, frame::Symbol = :tcc, values...) =
    add_part!(ts, :frame; name = name === nothing ? _free_name(ts, :frame) : name, frame, values...)

"""
    add_port_frame!(ts::TargetSetup, port, theta, phi; roll = 0, name = nothing) -> Bool

Add the frame of the port named `port` at polar angle `theta` and azimuth `phi`
(degrees): `Briareus.port_frame(theta, phi; roll)`, at TCC with `z` toward the port (a beam from it
travels along `-z`). Its name defaults to `port_<port>`. Only its roll can be edited
afterward; another frame placed in it carries offsets along the port axis.
"""
function add_port_frame!(ts::TargetSetup, port::AbstractString, θ::Real, φ::Real; roll::Real = 0.0,
                         name = nothing)
    nm = name === nothing ? _free_name(ts, Symbol(:port_, replace(String(port), r"[^A-Za-z0-9_]" => "_"));
                                       numbered = false) : name
    part = try
        isempty(port) && throw(ArgumentError("a port frame needs its port's name"))
        PartSpec(:frame; name = nm, object = port,
                 values = Dict(:theta => Float64(θ), :phi => Float64(φ), :rotation => Float64(roll)))
    catch err
        err isa ArgumentError || rethrow()
        ts.error[] = _errmsg(err)
        return false
    end
    return _try_build!(ts; parts = [ts.parts[]; part])
end

# A name not taken by a component, part or frame: `base_1`, `base_2`, … (or `base`
# itself first, unless `numbered`).
function _free_name(ts::TargetSetup, base::Symbol; numbered::Bool = true)
    taken = union(Set(component_names(ts.target[])), (p.name for p in ts.parts[]), keys(part_frames(ts)))
    !numbered && !(base in taken) && return base
    k = 1
    while Symbol(base, :_, k) in taken
        k += 1
    end
    return Symbol(base, :_, k)
end

"""
    edit_part!(ts::TargetSetup, i, key, text) -> Bool

Set a field of part `i` from text: a [`PART_FIELDS`](@ref) or size key, or `:name`.
"""
function edit_part!(ts::TargetSetup, i::Integer, key::Symbol, text)
    _check_part_index(ts, i)
    p = ts.parts[][i]
    try
        p.kind === :frame && key in (:albedo, :conversion, :h) && throw(ArgumentError("a frame has no field :$key"))
        is_port_frame(p) && key in (:theta, :phi) &&
            throw(ArgumentError("a port frame's axis is its port's: change its roll, or place another frame in it"))
        if key === :name
            s = strip(String(text))
            occursin(r"^[A-Za-z_][A-Za-z0-9_]*$", s) ||
                throw(ArgumentError("a part name is letters, digits and _, starting with a letter"))
            Symbol(s) == p.name && return true
            new = _with(p; name = Symbol(s))
            # parts placed in this one's frame or its pieces', or with steps that refer to
            # them, follow the rename, and so do the beams aimed in them
            renames = Dict(p.name => new.name)
            for ((a, _), (b, _)) in zip(piece_frame_defs(p), piece_frame_defs(new))
                renames[a] = b
            end
            ps = [foldl((q, (a, b)) -> _renamed(q, a, b), renames; init = q) for q in ts.parts[]]
            ps[i] = new
            return _try_build!(ts; parts = ps, renames)
        elseif key in (:pos_r, :pos_phi, :pos_theta)
            x, y, z = _set_position(p, key, _parse_number(text, String(key)))
            vals = copy(p.values)
            vals[:x], vals[:y], vals[:z] = x, y, z
            new = _with(p; values = vals)
        else
            haskey(p.values, key) || throw(ArgumentError("a $(p.kind) has no field :$key"))
            v = _parse_number(text, String(key))
            v == p.values[key] && (ts.error[] = ""; return true)
            vals = copy(p.values)
            vals[key] = v
            new = _with(p; values = vals)
        end
        ps = copy(ts.parts[])
        ps[i] = new
        return _try_build!(ts; parts = ps)
    catch err
        err isa ArgumentError || rethrow()
        ts.error[] = _errmsg(err)
        return false
    end
end

"""
    set_part_coords!(ts::TargetSetup, i, coords::Symbol) -> Bool

Show and edit part `i`'s position in other coordinates ([`PART_COORDINATES`](@ref)); the
part does not move.
"""
function set_part_coords!(ts::TargetSetup, i::Integer, coords::Symbol)
    _check_part_index(ts, i)
    coords in PART_COORDINATES || throw(ArgumentError("coordinates must be one of $PART_COORDINATES"))
    p = ts.parts[][i]
    p.coords === coords && return true
    ps = copy(ts.parts[])
    ps[i] = _with(p; coords)
    ts.parts[] = ps                     # the geometry is unchanged: no rebuild
    return true
end

"""
    set_part_flag!(ts::TargetSetup, i, flag, on::Bool) -> Bool

Set `:two_sided`, `:inward` or `:aperture` of part `i`.
"""
function set_part_flag!(ts::TargetSetup, i::Integer, flag::Symbol, on::Bool)
    _check_part_index(ts, i)
    flag in (:two_sided, :inward, :aperture) || throw(ArgumentError("unknown part flag :$flag"))
    p = ts.parts[][i]
    getfield(p, flag) == on && return true
    ps = copy(ts.parts[])
    ps[i] = _with(p; (flag => on,)...)
    return _try_build!(ts; parts = ps)
end

"""
    set_part_frame!(ts::TargetSetup, i, frame::Symbol; keep = true) -> Bool

Place part `i` in another frame (`:tcc`, another part's or a hohlraum part's piece
frame: see [`part_frames`](@ref)). With `keep` the part stays where it is in the
chamber: its own frame is re-expressed in the new one (its center, and the angles of
its rotation there; the angles stay as they are when the two frames have the same
axes); otherwise its placement values are kept and it moves
with the frame. A frame that leads back to the part is a loop: refused, with `error`
set.
"""
function set_part_frame!(ts::TargetSetup, i::Integer, frame::Symbol; keep::Bool = true)
    _check_part_index(ts, i)
    p = ts.parts[][i]
    p.frame === frame && return true
    try
        _own_frame_name(p.name, frame) && throw(ArgumentError("a part cannot be placed in its own frame"))
        trial = copy(ts.parts[])
        trial[i] = _with(p; frame)
        part_frames(trial)                     # a loop, a missing frame
        vals = copy(p.values)
        if keep
            fr = part_frames(ts)
            P0, P1 = fr[p.frame], fr[frame]
            vals[:x], vals[:y], vals[:z] = to_frame(P1, to_chamber(P0, part_center(p)))
            if P0.R != P1.R
                local_R = P1.R' * (P0.R * euler_rotation(p.values[:theta], p.values[:phi], p.values[:rotation]))
                vals[:theta], vals[:phi], vals[:rotation] = _euler_angles(local_R)
            end
        end
        ps = copy(ts.parts[])
        ps[i] = _with(p; frame, values = vals)
        return _try_build!(ts; parts = ps)
    catch err
        err isa ArgumentError || rethrow()
        ts.error[] = _errmsg(err)
        return false
    end
end

"""
Kinds of placement step a part can have ([`add_part_step!`](@ref)), in menu order:
relations that fix degrees of freedom (`:coincident`, `:on_axis`, `:concentric`,
`:distance`, `:facing`, `:parallel`, `:roll`) and relative moves (`:offset`, `:rotate`).
"""
const PART_STEP_KINDS = (:coincident, :on_axis, :concentric, :distance, :facing, :parallel,
                         :roll, :offset, :rotate)

const _STEP_TYPES = Dict(:coincident => Coincident, :on_axis => OnAxis, :concentric => Concentric,
                         :distance => Distance, :facing => Facing, :parallel => Parallel,
                         :roll => Roll, :offset => Offset, :rotate => Rotate)

"Menu labels of [`PART_STEP_KINDS`](@ref)."
const PART_STEP_LABELS = Dict(:coincident => "Coincident", :on_axis => "On axis",
                              :concentric => "Concentric", :distance => "Distance",
                              :facing => "Facing", :parallel => "Parallel", :roll => "Roll",
                              :offset => "Offset", :rotate => "Rotate")

"""
    step_kind(step) -> Symbol

The [`PART_STEP_KINDS`](@ref) entry of a placement step.
"""
step_kind(st::PlacementStep) = only(k for (k, T) in _STEP_TYPES if st isa T)

# A step as a dictionary of its fields (the session file; `nothing` fields left out), and back.
function _step_state(st::PlacementStep)
    d = Dict{String,Any}("kind" => String(step_kind(st)))
    for f in fieldnames(typeof(st))
        v = getfield(st, f)
        v === nothing && continue
        d[String(f)] = v isa Symbol ? String(v) : v isa Vec3 ? collect(v) : v
    end
    return d
end

function _state_step(d::AbstractDict)
    kind = Symbol(d["kind"])
    haskey(_STEP_TYPES, kind) || throw(ArgumentError("unknown placement step \"$kind\""))
    T = _STEP_TYPES[kind]
    args = map(fieldnames(T)) do f
        v, FT = get(d, String(f), nothing), fieldtype(T, f)
        if v === nothing
            Nothing <: FT || throw(ArgumentError("a $kind step needs its $f"))
            return nothing
        end
        FT === Vec3 ? Vec3(Float64.(v)...) :
        FT <: Union{Nothing,Symbol} ? Symbol(v) :
        FT === Bool ? Bool(v) : Float64(v)
    end
    return _checked_step(T(args...))
end

# The checks of the steps' keyword constructors, for steps built from their fields.
function _checked_step(st::PlacementStep)
    st isa Distance && st.distance < 0 && throw(ArgumentError("a distance must not be negative"))
    st isa Facing && !(st.side in (:front, :back)) && throw(ArgumentError("side must be front or back"))
    for f in (:axis, :point, :delta)
        hasfield(typeof(st), f) || continue
        all(isfinite, getfield(st, f)) || throw(ArgumentError("the $f must be finite"))
    end
    return st
end

# A new step of `kind` on the frame `ref`, with values that leave the part's frame `F`
# where it is where the step allows (`fr`: the design's frames).
function _step_default(kind::Symbol, ref::Symbol, F::Frame, fr)
    R = fr[ref]
    z = F.R[:, 3]
    kind === :coincident && return Coincident(ref; point = to_frame(R, F.origin))
    kind === :on_axis && return OnAxis(ref)
    kind === :concentric && return Concentric(ref; flip = dot(z, R.R[:, 3]) < 0)
    kind === :distance && return Distance(ref, norm(F.origin - R.origin))
    kind === :facing && return Facing(ref; side = dot(z, R.origin - F.origin) < 0 ? :back : :front)
    kind === :parallel && return Parallel(ref; flip = dot(z, R.R[:, 3]) < 0)
    if kind === :roll
        for a in (:x, :y)
            r = to_chamber(R, Briareus._step_axis(a); direction = true)
            e = r - dot(r, z) * z
            norm(e) > 1e-9 || continue
            e = normalize(e)
            return Roll(atand(dot(F.R[:, 1], cross(z, e)), dot(F.R[:, 1], e)); from = (ref, a))
        end
    end
    kind === :offset && return Offset(0.0, 0.0, 0.0)
    kind === :rotate && return Rotate(:z, 0.0)
    throw(ArgumentError("placement step kind must be one of $PART_STEP_KINDS"))
end

function _set_steps!(ts::TargetSetup, i::Integer, steps)
    ps = copy(ts.parts[])
    ps[i] = _with(ps[i]; steps = collect(PlacementStep, steps))
    return _try_build!(ts; parts = ps)
end

function _step_edit(f, ts::TargetSetup)
    try
        return f()
    catch err
        err isa ArgumentError || rethrow()
        ts.error[] = _errmsg(err)
        return false
    end
end

"""
    add_part_step!(ts::TargetSetup, i, kind; ref = the part's frame) -> Bool

Append a placement step of `kind` ([`PART_STEP_KINDS`](@ref)) to part `i`, referring to
the frame `ref`, with values that leave the part where it is where the step allows
(coincident at its origin's current place, the current distance, roll and flips); a
relative move starts at zero. A step that cannot apply, or that undoes an earlier one,
is refused with `error` set.
"""
function add_part_step!(ts::TargetSetup, i::Integer, kind::Symbol; ref::Union{Nothing,Symbol} = nothing)
    _check_part_index(ts, i)
    p = ts.parts[][i]
    return _step_edit(ts) do
        fr = part_frames(ts)
        r = something(ref, p.frame)
        r === p.name && throw(ArgumentError("a part's placement step cannot refer to its own frame"))
        haskey(fr, r) || throw(ArgumentError("no frame named $r"))
        _set_steps!(ts, i, [p.steps; _step_default(kind, r, fr[p.name], fr)])
    end
end

"""
    remove_part_step!(ts::TargetSetup, i, k) -> Bool
    move_part_step!(ts::TargetSetup, i, k, by) -> Bool

Remove placement step `k` of part `i`, or move it `by` places (−1 earlier, +1 later;
the order matters: steps apply in sequence).
"""
function remove_part_step!(ts::TargetSetup, i::Integer, k::Integer)
    _check_part_index(ts, i)
    steps = copy(ts.parts[][i].steps)
    checkbounds(steps, k)
    deleteat!(steps, k)
    return _step_edit(() -> _set_steps!(ts, i, steps), ts)
end

function move_part_step!(ts::TargetSetup, i::Integer, k::Integer, by::Integer)
    _check_part_index(ts, i)
    steps = copy(ts.parts[][i].steps)
    checkbounds(steps, k)
    j = clamp(k + by, 1, length(steps))
    j == k && return true
    st = popat!(steps, k)
    insert!(steps, j, st)
    return _step_edit(() -> _set_steps!(ts, i, steps), ts)
end

"""
    edit_part_step!(ts::TargetSetup, i, k, key, value) -> Bool

Set a field of placement step `k` of part `i`: `:ref` (a frame name; `""` for none: the
part's own axes for an offset, a relative roll), `:px`/`:py`/`:pz` (the point in `ref`,
mm), `:dx`/`:dy`/`:dz` (an offset, mm), `:distance` (mm; `""` for none on an on-axis
step), `:angle` (degrees), `:axis` (`"x"`, `"y"` or `"z"`), `:flip` (Bool) or `:side`
(`"front"` or `"back"`); numbers as text. A rejected value sets `error`.
"""
function edit_part_step!(ts::TargetSetup, i::Integer, k::Integer, key::Symbol, value)
    _check_part_index(ts, i)
    p = ts.parts[][i]
    checkbounds(p.steps, k)
    st = p.steps[k]
    return _step_edit(ts) do
        d = _step_state(st)
        kind = d["kind"]
        coord = _STEP_COORDS
        if key === :ref
            s = strip(string(value))
            hasfield(typeof(st), :ref) || throw(ArgumentError("a $kind step has no frame"))
            isempty(s) && !(st isa Union{Roll,Offset}) && throw(ArgumentError("a $kind step needs a frame"))
            Symbol(s) === p.name && throw(ArgumentError("a part's placement step cannot refer to its own frame"))
            isempty(s) ? delete!(d, "ref") : (d["ref"] = s)
        elseif haskey(coord, key)
            f, c = coord[key]
            haskey(d, f) || throw(ArgumentError("a $kind step has no $f"))
            d[f][c] = _parse_number(value, String(key))
        elseif key === :distance && st isa OnAxis && isempty(strip(string(value)))
            delete!(d, "distance")
        elseif key in (:distance, :angle)
            hasfield(typeof(st), key) || throw(ArgumentError("a $kind step has no $key"))
            d[String(key)] = _parse_number(value, String(key))
        elseif key === :axis
            haskey(d, "axis") || throw(ArgumentError("a $kind step has no axis"))
            d["axis"] = collect(Briareus._step_axis(Symbol(value)))
        elseif key === :flip
            haskey(d, "flip") || throw(ArgumentError("a $kind step has no flip"))
            d["flip"] = Bool(value)
        elseif key === :side
            haskey(d, "side") || throw(ArgumentError("a $kind step has no side"))
            d["side"] = String(value)
        else
            throw(ArgumentError("unknown placement step field :$key"))
        end
        new = _state_step(d)
        _step_state(new) == _step_state(st) && (ts.error[] = ""; return true)
        steps = copy(p.steps)
        steps[k] = new
        _set_steps!(ts, i, steps)
    end
end

"""
    part_placement(ts::TargetSetup, i) -> PlacementReport

What part `i`'s placement steps fix (`placement_report`): the degrees of freedom per
step, redundant steps, what stays free, and the residuals left by relative moves.
"""
function part_placement(ts::TargetSetup, i::Integer)
    _check_part_index(ts, i)
    ps = ts.parts[]
    defs = [d for p in ps for d in (p.name => part_frame_def(p), piece_frame_defs(p)...)]
    return placement_report(defs)[ps[i].name]
end

# The fields of a step's point and offset coordinates: (field, index).
const _STEP_COORDS = Dict(:px => ("point", 1), :py => ("point", 2), :pz => ("point", 3),
                          :dx => ("delta", 1), :dy => ("delta", 2), :dz => ("delta", 3))

"""
    step_numbers(step) -> Vector{Tuple{Symbol,String,String}}

The numbers of a placement step, as `(key, label, unit)`, keys as in
[`edit_part_step!`](@ref): a point's or an offset's coordinates (mm), a distance (mm;
an on-axis step's only when set) and an angle (degrees). A sweep can vary each
([`set_sweep_step!`](@ref)).
"""
function step_numbers(st::PlacementStep)
    out = Tuple{Symbol,String,String}[]
    hasfield(typeof(st), :point) && append!(out, [(:px, "point x", "mm"), (:py, "point y", "mm"), (:pz, "point z", "mm")])
    hasfield(typeof(st), :delta) && append!(out, [(:dx, "offset x", "mm"), (:dy, "offset y", "mm"), (:dz, "offset z", "mm")])
    (st isa Distance || st isa OnAxis && st.distance !== nothing) && push!(out, (:distance, "distance", "mm"))
    hasfield(typeof(st), :angle) && push!(out, (:angle, "angle", "°"))
    return out
end

"""
    step_number(step, key) -> Float64

The value of the number `key` of a placement step ([`step_numbers`](@ref)).
"""
function step_number(st::PlacementStep, key::Symbol)
    haskey(_STEP_COORDS, key) && return Float64(getfield(st, Symbol(_STEP_COORDS[key][1]))[_STEP_COORDS[key][2]])
    v = hasfield(typeof(st), key) ? getfield(st, key) : nothing
    v isa Real || throw(ArgumentError("a $(step_kind(st)) step has no number $key"))
    return Float64(v)
end

# The step with its number `key` set to `v` (the checks of `_state_step`).
function _step_with(st::PlacementStep, key::Symbol, v::Real)
    step_number(st, key)                              # (a number the step has)
    d = _step_state(st)
    haskey(_STEP_COORDS, key) ? (d[_STEP_COORDS[key][1]][_STEP_COORDS[key][2]] = v) : (d[String(key)] = v)
    return _state_step(d)
end

"""
    _euler_angles(R) -> (polar, azimuth, rotation)

Angles (degrees) of a rotation matrix in `euler_rotation`'s order, so that
`euler_rotation(_euler_angles(R)...) ≈ R`: polar in [0, 180], azimuth and rotation in
[0, 360). With the axis along `±z` (polar 0 or 180, to 1e-12 in the sine) the azimuth is
0 and the whole turn is the rotation.
"""
function _euler_angles(R::AbstractMatrix)
    s = hypot(R[1, 3], R[2, 3])
    θ = atand(s, R[3, 3])
    if s < 1e-12
        ψ = R[3, 3] > 0 ? atand(R[2, 1], R[1, 1]) : atand(R[2, 1], R[2, 2])
        return θ, 0.0, mod(ψ, 360)
    end
    return θ, mod(atand(R[2, 3], R[1, 3]), 360), mod(atand(R[3, 2], -R[3, 1]), 360)
end

"""
    nudge_part!(ts::TargetSetup, i, Δ; frame = the part's frame, axis = :z) -> Bool

Move part `i` by `Δ` mm along `axis` (`:x`, `:y`, `:z`, or a vector, normalized) of the
frame named `frame` ([`part_frames`](@ref); by default the frame the part is placed in),
as one edit: the part's center moves by `Δ` times that axis in the chamber, its
orientation and the frame it is placed in stay. Parts placed in its frame, and beams
aimed in it, move with it. A rejected move sets `error` and returns `false`.
"""
function nudge_part!(ts::TargetSetup, i::Integer, Δ::Real; frame::Union{Nothing,Symbol} = nothing,
                     axis = :z)
    _check_part_index(ts, i)
    p = ts.parts[][i]
    frame = something(frame, p.frame)
    try
        isfinite(Δ) || throw(ArgumentError("the move must be a finite length"))
        fr = part_frames(ts)
        haskey(fr, frame) || throw(ArgumentError("no frame named $frame"))
        d = Δ * to_chamber(fr[frame], _axis_vector3(axis); direction = true)
        P = fr[p.frame]
        c = part_center(p) + to_frame(P, d; direction = true)
        vals = copy(p.values)
        vals[:x], vals[:y], vals[:z] = c
        ps = copy(ts.parts[])
        ps[i] = _with(p; values = vals)
        # placement steps that fix the position along the move leave nothing to move
        if !isempty(p.steps) && Δ != 0 &&
           norm(part_frames(ps)[p.name].origin - fr[p.name].origin) < 1e-12 * abs(Δ)
            throw(ArgumentError("its placement steps fix its position along this axis: edit or sweep a step instead"))
        end
        return _try_build!(ts; parts = ps)
    catch err
        err isa ArgumentError || rethrow()
        ts.error[] = _errmsg(err)
        return false
    end
end

# A frame axis: `:x`, `:y`, `:z` or a vector (normalized).
function _axis_vector3(axis)
    axis === :x && return Vec3(1, 0, 0)
    axis === :y && return Vec3(0, 1, 0)
    axis === :z && return Vec3(0, 0, 1)
    axis isa Symbol && throw(ArgumentError("the axis is :x, :y, :z or a vector"))
    v = Vec3(axis...)
    norm(v) > 0 || throw(ArgumentError("the axis must not be zero"))
    return normalize(v)
end

"""
Roles a part can take in the analyses, held by the part's flags (its name carries none):
`:none`; `:capsule`, the observed role held by a sphere (flag `observed`: its x-ray
drive's modes, the θ–φ map and the capsule keep-out; one part at a time, see
[`observed_name`](@ref)); `:foil` (a backlighter foil, a flat part or a mesh: any number
of parts can be foils, and the active one ([`active_foil`](@ref)) is the foil whose spot
statistics, map, objective, line of sight and captured power the analyses take); and
`:wall` (a re-emitting cavity surface: the drive history's wall temperature, and a
surface the clearance and TaLIS checks land the beams on; any number of parts, see
[`wall_names`](@ref)).
"""
const PART_ROLES = (:none, :capsule, :foil, :wall)

const PART_ROLE_LABELS = Dict(:none => "None", :capsule => "Capsule (observed x-ray drive)",
                              :foil => "Backlighter foil (spot uniformity)",
                              :wall => "Wall (re-emitting cavity surface)")

"Role of a part (one of [`PART_ROLES`](@ref)), from its flags."
part_role(p::PartSpec) = p.observed ? :capsule : p.foil ? :foil : p.wall ? :wall : :none

"""
    observed_name(ts::TargetSetup) -> Union{Symbol,Nothing}

The surface whose x-ray drive the analyses observe (the drive's modes and θ–φ map, the
capsule keep-out and summaries, the drive history, the retune, the results export): the
part with the observed role (`PartSpec.observed`, given by [`set_part_role!`](@ref); a
starter's capsule holds it) while it is in the calculations, else `nothing`. One surface
at a time.
"""
function observed_name(ts::TargetSetup)
    ex = ts.excluded[]
    for p in ts.parts[]
        p.observed && !(p.name in ex) && return p.name
    end
    return nothing
end

"""
    wall_names(ts::TargetSetup) -> Vector{Symbol}

The re-emitting cavity walls in the calculations (the drive history's wall temperature;
surfaces the clearance and TaLIS checks land the beams on), in the parts' order: each
part with the role `:wall`, and each hohlraum part's wall and caps (its pieces
`<part>_wall`, `<part>_cap_top`, `<part>_cap_bottom`: they hold the role by being a
hohlraum's; the checks keep the caps as hardware, the entrance holes' edges).
"""
function wall_names(ts::TargetSetup)
    ex = ts.excluded[]
    out = Symbol[]
    for p in ts.parts[]
        p.name in ex && continue
        p.wall && push!(out, p.name)
        p.kind === :hohlraum && append!(out, (piece_name(p.name, s) for s in _HOHLRAUM_WALL))
    end
    return out
end

# The pieces of a hohlraum part that hold the wall role.
const _HOHLRAUM_WALL = (:wall, :cap_top, :cap_bottom)

"""
    foil_names(ts::TargetSetup) -> Vector{Symbol}

The backlighter foils of the design in the calculations: each part with the role
`:foil`, in the parts' order. Foils left out of the calculations are not listed.
"""
foil_names(ts::TargetSetup) = Symbol[p.name for p in ts.parts[] if p.foil && !(p.name in ts.excluded[])]

"""
    active_foil(ts::TargetSetup) -> Union{Symbol,Nothing}

The foil the foil analyses take (the drive's spot statistics and foil map, the foil
objective, the sweep's line of sight and captured power): the one chosen with
[`set_active_foil!`](@ref) while it is still a foil in the calculations, else the first
of [`foil_names`](@ref), else `nothing`. A part given the foil role when the design has
no foil becomes the chosen one.
"""
function active_foil(ts::TargetSetup)
    fs = foil_names(ts)
    isempty(fs) && return nothing
    a = ts.active_foil[]
    return a in fs ? a : first(fs)
end

"""
    set_active_foil!(ts::TargetSetup, name::Symbol) -> Bool

Make foil `name` (one of [`foil_names`](@ref)) the one the foil analyses take
([`active_foil`](@ref)). A name that is not a foil in the calculations sets `error` and
returns `false`. Saved with the session, and undone like any edit.
"""
function set_active_foil!(ts::TargetSetup, name::Symbol)
    if !(name in foil_names(ts))
        ts.error[] = "$name is not a backlighter foil in the calculations"
        return false
    end
    ts.error[] = ""
    ts.active_foil[] === name || (ts.active_foil[] = name)
    return true
end

"""
    foil_components(ts::TargetSetup, t = ts.target[]) -> Vector{Symbol}

The components of `t` that belong to a backlighter foil, active or not: the components
of each part with the foil role (a VisRad part's carry its name as a prefix). The
clearance and TaLIS checks take them as surfaces the beams land on.
"""
foil_components(ts::TargetSetup, t::Target = ts.target[]) = _held_components(ts, t, p -> p.foil, ())

"""
    observed_components(ts::TargetSetup, t = ts.target[]) -> Vector{Symbol}
    wall_components(ts::TargetSetup, t = ts.target[]) -> Vector{Symbol}

The components of `t` that hold the observed role ([`observed_name`](@ref)) or the wall
role ([`wall_names`](@ref)): the components of the parts with the role (a VisRad part's
carry its name as a prefix) and, for the wall, each hohlraum part's wall and caps.
"""
function observed_components(ts::TargetSetup, t::Target = ts.target[])
    o = observed_name(ts)
    o === nothing && return Symbol[]
    return _held_components(ts, t, p -> p.name === o, ())
end

wall_components(ts::TargetSetup, t::Target = ts.target[]) =
    _held_components(ts, t, p -> p.wall,
                     [piece_name(p.name, s) for p in ts.parts[] if p.kind === :hohlraum for s in _HOHLRAUM_WALL])

# The components of `t` of the parts `holds(p)` selects (a VisRad part's carry its name as
# a prefix), and the components named `named` (a hohlraum part's walls).
function _held_components(ts::TargetSetup, t::Target, holds, named)
    fs = Set{Symbol}(p.name for p in ts.parts[] if holds(p))
    vr = Set{Symbol}(p.name for p in ts.parts[] if holds(p) && p.kind === :visrad)
    out = Symbol[]
    for c in t.components
        n = c.name
        if n in named || n in fs || c.part in fs || any(v -> startswith(String(n), String(v)), vr)
            push!(out, n)
        end
    end
    return out
end

# The surfaces the clearance check lands the beams on by default (`lands_on`: walls, the
# observed surface, the foils) and those the TaLIS check leaves out (`observed`), as part
# and component names of `t` (the checks list parts: the parts of the wall and observed
# components join them; the foils by component).
function _landing(ts::TargetSetup, t::Target = ts.target[])
    obs = _with_parts(t, observed_components(ts, t))
    # (a hohlraum part's caps hold the wall role for the drive, but the checks keep them
    # as hardware: they are the entrance holes' edges, as the former template's were)
    walls = filter(n -> !occursin(r"_cap_(top|bottom)$", String(n)) || !_piece_of_hohlraum(ts, n),
                   wall_components(ts, t))
    capture = unique!([_with_parts(t, walls); foil_components(ts, t)])
    return (; lands_on = unique!([capture; obs]), capture, observed = obs)
end

# Components `ns` of `t` and their parts (a piece of a part in pieces stands alone).
_with_parts(t::Target, ns) = unique!([ns; [c.part for c in t.components if c.name in ns && !_is_piece(c)]])

# Whether component `name` is a piece of one of the design's hohlraum parts.
_piece_of_hohlraum(ts::TargetSetup, name::Symbol) =
    any(p -> p.kind === :hohlraum && startswith(String(name), String(p.name) * "_"), ts.parts[])

# The walls' and the observed surface's components and parts, as sets.
_held_sets(ts::TargetSetup, t::Target = ts.target[]) =
    (; walls = Set(setdiff(_landing(ts, t).capture, foil_components(ts, t))),
     observed = Set(_with_parts(t, observed_components(ts, t))))

# The object kind a part stands for (a VisRad object's own kind).
function _part_shape(p::PartSpec)
    p.kind === :visrad || return p.kind
    k = _visrad_object(p.file, p.object).kind
    return k === :mesh ? :mesh : k
end

"""
    set_part_role!(ts::TargetSetup, i, role::Symbol; done = () -> nothing) -> Bool

Give part `i` a role in the analyses (one of [`PART_ROLES`](@ref)); the part keeps its
name, and its components (and so the target) are unchanged. `:foil` (a disk, rectangle
or mesh, or a VisRad rectangle or CAD part: a backlighter whose spot the foil objective
makes uniform) marks the part as a foil; the part becomes the active foil
([`active_foil`](@ref)) when the design had none. `:capsule` makes the part the observed
surface ([`observed_name`](@ref)); another part that held it gives it up. A hohlraum
part takes no role (its wall and caps hold the wall role). The capsule is a sphere (the mode analysis and the
keep-out take its center and radius): a sphere part or a VisRad sphere takes the role
itself; a CAD part (a capsule with its fill tube, say) is replaced in the calculations by
a sphere part fitted to it ([`fit_sphere`](@ref): named `capsule`, or `capsule_2`, … when
that is taken; its element size the target's capsule size, its albedo and conversion the
part's), which takes the role, and left out of the calculations itself. `:wall` (any
surface but an aperture) marks a re-emitting cavity wall ([`wall_names`](@ref)). `:none`
clears the part's role. A role the part cannot take sets `error` and returns `false`. In
a window (`spawn[]`) a large CAD capsule is fitted on a worker (see
[`TargetSetup`](@ref)): `true` then means the fit is under way, and `done()` runs once the
role is given (at once otherwise).
"""
function set_part_role!(ts::TargetSetup, i::Integer, role::Symbol; done = () -> nothing)
    _check_part_index(ts, i)
    role in PART_ROLES || throw(ArgumentError("a part's role must be one of $PART_ROLES"))
    ps = copy(ts.parts[])
    p = ps[i]
    part_role(p) === role && (done(); return true)
    try
        p.kind === :frame && throw(ArgumentError("a frame has no surface to take a role"))
        p.kind === :hohlraum &&
            throw(ArgumentError("a hohlraum's wall and caps hold the wall role already; its holes are apertures"))
        shape = _part_shape(p)
        if role === :capsule && shape !== :sphere
            shape === :mesh ||
                throw(ArgumentError("the capsule is a sphere, or a CAD part a sphere is fitted to"))
            return _capsule_from_part!(ts, i; done)
        end
        role === :foil && !(shape in (:rectangle, :disk, :mesh)) &&
            throw(ArgumentError("the backlighter foil is a flat part (a rectangle or disk) or a mesh"))
        role === :wall && p.aperture &&
            throw(ArgumentError("an aperture is transparent: a wall is a surface"))
    catch err
        err isa ArgumentError || rethrow()
        ts.error[] = _errmsg(err)
        return false
    end
    role === :capsule && _drop_observed!(ps)            # one observed surface at a time
    ps[i] = _with(ps[i]; foil = role === :foil, observed = role === :capsule, wall = role === :wall)
    first_foil = role === :foil && active_foil(ts) === nothing
    # the geometry is unchanged: the parts are committed without building the target
    ts.error[] = ""
    ts.parts[] = ps
    first_foil && (ts.active_foil[] = ps[i].name)
    done()
    return true
end

# Clear the observed role of every part of `ps` (in place).
function _drop_observed!(ps::Vector{PartSpec})
    for k in eachindex(ps)
        ps[k].observed && (ps[k] = _with(ps[k]; observed = false))
    end
    return ps
end

"""
    fit_sphere(points, weights = ones; inlier = 0.1) -> (center, radius)

A sphere through most of `points` (mm), weighted by `weights` (e.g. the surface area
each point stands for, so a finely tessellated fill tube does not outweigh the capsule
it is attached to): a sphere fitted to every point (weighted least squares, algebraic),
then refitted to the points within 0.3, then four times within `inlier`, × its radius.
Robust to a minority of the weight off the sphere. A shell modeled with both faces
fits between them; the radius returned is the outer face's (the weighted 75th
percentile of the distances within 20 % of the fit).
"""
function fit_sphere(points::AbstractVector, weights::AbstractVector = ones(length(points));
                    inlier::Real = 0.1)
    P = [Vec3(p...) for p in points]
    w = Float64.(weights)
    length(P) >= 4 && length(w) == length(P) && all(>=(0), w) && sum(w) > 0 ||
        throw(ArgumentError("a sphere needs at least four points with nonnegative weights"))
    function wmedian(x, w)
        o = sortperm(x)
        cw = cumsum(w[o])
        return x[o[searchsortedfirst(cw, cw[end] / 2)]]
    end
    c = sum(w .* P) / sum(w)
    r = wmedian([norm(p - c) for p in P], w)
    for it in 1:6
        # every point first, then those near the sphere
        band = (it == 1 ? Inf : it == 2 ? 0.3 : Float64(inlier)) * r
        sel = [k for k in eachindex(P) if abs(norm(P[k] - c) - r) <= band && w[k] > 0]
        length(sel) >= 4 || break
        # |p|² = 2 c·p + (r² − |c|²): linear in (c, k); rows weighted by √w
        sw = sqrt.(w[sel])
        A = reduce(vcat, [sw[j] * [2P[k][1] 2P[k][2] 2P[k][3] 1.0] for (j, k) in enumerate(sel)])
        b = [sw[j] * dot(P[k], P[k]) for (j, k) in enumerate(sel)]
        x = A \ b
        c = Vec3(x[1], x[2], x[3])
        r = sqrt(max(x[4] + dot(c, c), 0.0))
    end
    # a shell's inner and outer faces fit in between: take the outer one (the radiating
    # surface), the weighted 75th percentile of the distances near the sphere
    near = [k for k in eachindex(P) if abs(norm(P[k] - c) - r) <= 0.2r]
    if length(near) >= 4
        d = [norm(P[k] - c) for k in near]
        o = sortperm(d)
        cw = cumsum(w[near][o])
        r = d[o[searchsortedfirst(cw, 0.75 * cw[end])]]
    end
    return c, r
end

# Role `:capsule` for CAD part `i`: a sphere part with the role, fitted to it, in its place in
# the calculations. With `spawn[]` a part of more than `BACKGROUND_TRIANGLES` facets is
# fitted on a worker (`fitting`), and the role is given when the fit is done, to the part
# as it is then (fitted again if it moved meanwhile; dropped if it went away).
function _capsule_from_part!(ts::TargetSetup, i::Integer; done = () -> nothing)
    p = ts.parts[][i]
    cs = part_components(p, ts.mesh[].h; frame = part_frames(ts)[p.name], raw = true)
    if ts.spawn[] && sum(c -> length(c.mesh.triangles), cs; init = 0) > BACKGROUND_TRIANGLES[]
        job = ts.fit_jobs[p.name] = get(ts.fit_jobs, p.name, 0) + 1
        _assign!(ts.fitting, union(ts.fitting[], (p.name,)))
        ts.error[] = ""
        deliver = ts.deliver[]
        _run_job(() -> (out = _try_job(() -> _fit_capsule(cs, p.name));
                        deliver(() -> _capsule_fitted!(ts, p, job, out, done))), true)
        return true
    end
    ok = _apply_capsule!(ts, i, _fit_capsule(cs, p.name)...)
    ok && done()
    return ok
end

# The sphere (center, radius) of a CAD capsule's components (raw facets in the chamber).
function _fit_capsule(cs, name)
    # each triangle's vertices weighted by a third of its area
    pts, wts = Vec3[], Float64[]
    for comp in cs, (a, b, d) in comp.mesh.triangles
        va, vb, vd = comp.mesh.vertices[a], comp.mesh.vertices[b], comp.mesh.vertices[d]
        A = norm(cross(vb - va, vd - va)) / 6
        push!(pts, va, vb, vd)
        append!(wts, (A, A, A))
    end
    c, r = fit_sphere(pts, wts)
    r > 0 || throw(ArgumentError("$name: no sphere fits it"))
    return c, r
end

# A worker fitted the capsule to part `p` (job `job`): give it the role now.
function _capsule_fitted!(ts::TargetSetup, p::PartSpec, job::Integer, out, done = () -> nothing)
    get(ts.fit_jobs, p.name, 0) == job || return false       # asked again, or canceled
    delete!(ts.fit_jobs, p.name)
    ok = try
        i = findfirst(q -> q === p, ts.parts[])
        j = part_index(ts, p.name)
        if out.err !== nothing
            ts.error[] = _errmsg(out.err)
            false
        elseif i !== nothing
            _apply_capsule!(ts, i, out.value...) && (done(); true)
        elseif j !== nothing                                  # edited meanwhile: fit it again
            set_part_role!(ts, j, :capsule; done)
        else
            ts.error[] = "$(p.name) went away before its capsule sphere was fitted"
            false
        end
    catch err
        err isa ArgumentError || rethrow()
        ts.error[] = _errmsg(err)
        false
    end
    # the role is committed (its target building, or built) before the fit stops
    # holding the drive back, unless the part is being fitted again
    haskey(ts.fit_jobs, p.name) || _assign!(ts.fitting, setdiff(ts.fitting[], (p.name,)))
    return ok
end

# Drop the capsule fits in flight (an undo, a session loaded, a window closed).
function _cancel_fits!(ts::TargetSetup)
    empty!(ts.fit_jobs)
    _assign!(ts.fitting, Set{Symbol}())
    return ts
end

function _apply_capsule!(ts::TargetSetup, i::Integer, c, r)
    ps = _drop_observed!(copy(ts.parts[]))
    p = ps[i]
    ex = copy(ts.excluded[])
    taken = union(Set(all_component_names(ts)), (q.name for q in ps))
    name, k = :capsule, 1
    while name in taken
        k += 1
        name = Symbol(:capsule_, k)
    end
    push!(ps, PartSpec(:sphere; name, observed = true,
                       values = Dict(:x => c[1], :y => c[2], :z => c[3], :r1 => r,
                                     :h => ts.mesh[].capsule_h, :albedo => p.values[:albedo],
                                     :conversion => p.values[:conversion])))
    push!(ex, p.name)
    return _try_build!(ts; parts = ps, excluded = ex)
end

"""
    part_index(ts::TargetSetup, name::Symbol) -> Union{Nothing,Int}
"""
part_index(ts::TargetSetup, name::Symbol) = findfirst(p -> p.name === name, ts.parts[])

# ---------------------------------------------------------------- starters

"""
Starter assemblies, inserted with [`insert_starter!`](@ref): `:hohlraum` (*Hohlraum with
capsule*: a hohlraum part and a capsule sphere with the observed role) and `:foil`
(*Backlighter foil*: a foil rectangle facing TCC and an absorbing sphere standing in for
the main target). Their dimensions are illustrative (the NIF-scale defaults of the
former templates), not a facility's design.
"""
const STARTERS = (:hohlraum, :foil)

const STARTER_LABELS = Dict(:hohlraum => "Hohlraum with capsule", :foil => "Backlighter foil")

"""
    starter_parts(kind; frame = :tcc, taken = (), capsule_h = 0.3) -> Vector{PartSpec}

The parts of starter `kind` (one of [`STARTERS`](@ref)), placed in the frame named
`frame`, their names made unique against `taken` (`hohlraum`, then `hohlraum_2`, …).

- `:hohlraum`: a `:hohlraum` part `hohlraum` (length 10.1 mm, radius 2.875 mm, entrance
  holes of 1.55 mm radius, wall albedo 0.85 and conversion 0.8) and, placed in its frame
  at its center, the sphere `capsule` (radius 0.9 mm, albedo and conversion 0.3: the
  ablator's, meshed as a capsule, `capsule_mesh`) with the observed role. At TCC it
  meshes as the former hohlraum template did, element for element.
- `:foil`: the rectangle `foil` (2 mm square, conversion 0.02, its element size the
  target's capsule size `capsule_h` at the insert), 6 mm from the frame's origin
  along its `x` and facing it (its edges along the frame's `z` and `y`), with the foil
  role, and the absorbing sphere `main` (radius 1 mm, albedo and conversion 0) at the
  origin. At TCC both mesh as the former foil template's did.
"""
function starter_parts(kind::Symbol; frame::Symbol = :tcc, taken = (), capsule_h::Real = 0.3)
    kind in STARTERS || throw(ArgumentError("a starter is one of $STARTERS"))
    used = Set{Symbol}(taken)
    function unique_name(base::Symbol)
        name, k = base, 1
        while name in used || any(n -> startswith(String(n), String(name) * "."), used)
            k += 1
            name = Symbol(base, :_, k)
        end
        push!(used, name)
        return name
    end
    if kind === :hohlraum
        h = unique_name(:hohlraum)
        c = unique_name(:capsule)
        return [PartSpec(:hohlraum; name = h, frame),
                PartSpec(:sphere; name = c, frame = h, observed = true,
                         values = Dict(:r1 => 0.9, :albedo => 0.3, :conversion => 0.3, :capsule_mesh => 1.0))]
    end
    f = unique_name(:foil)
    m = unique_name(:main)
    return [PartSpec(:rectangle; name = f, frame, foil = true,
                     values = Dict(:x => 6.0, :theta => 90.0, :phi => 180.0, :rotation => 180.0,
                                   :r1 => 2.0, :r2 => 2.0,
                                   :conversion => 0.02, :h => Float64(capsule_h))),
            PartSpec(:sphere; name = m, frame, values = Dict(:r1 => 1.0))]
end

"""
    insert_starter!(ts::TargetSetup, kind; frame = :tcc) -> Bool

Add the parts of starter `kind` ([`starter_parts`](@ref)) to the design, placed in the
frame named `frame` (a part's, a piece's or TCC), their names made unique against the
design's parts and components (`hohlraum_2`, …). The first foil a starter brings
becomes the active foil when the design had none; a starter's observed surface takes
the role from the part that held it. A design the geometry rejects sets `error` and
returns `false`.
"""
function insert_starter!(ts::TargetSetup, kind::Symbol; frame::Symbol = :tcc)
    try
        haskey(part_frames(ts), frame) || throw(ArgumentError("no frame named $frame"))
        taken = union(Set(all_component_names(ts)), (p.name for p in ts.parts[]), keys(part_frames(ts)))
        new = starter_parts(kind; frame, taken, capsule_h = ts.mesh[].capsule_h)
        ps = copy(ts.parts[])
        any(p -> p.observed, new) && _drop_observed!(ps)
        first_foil = active_foil(ts) === nothing ? findfirst(p -> p.foil, new) : nothing
        _try_build!(ts; parts = [ps; new]) || return false
        first_foil === nothing || (ts.active_foil[] = new[first_foil].name)
        return true
    catch err
        err isa ArgumentError || rethrow()
        ts.error[] = _errmsg(err)
        return false
    end
end

# ---------------------------------------------------------------- summary

"""
    target_summary(ts::TargetSetup) -> String

The hohlraum parts and their windows, the other parts, and the element count.
"""
function target_summary(ts::TargetSetup)
    t = ts.target[]
    ps = ts.parts[]
    nh = count(p -> p.kind === :hohlraum, ps)
    nw = sum(p -> length(p.windows), ps; init = 0)
    no = count(p -> !(p.kind in (:hohlraum, :frame)), ps)
    kind = nh == 0 ? "target" : nh == 1 ? "hohlraum" : "$nh hohlraums"
    nw > 0 && (kind *= nw == 1 ? " with 1 window" : " with $nw windows")
    no > 0 && (kind *= nh == 0 ? (no == 1 ? " of 1 part" : " of $no parts") : no == 1 ? " + 1 part" : " + $no parts")
    w = target_waiting(ts)
    return "$kind · $(nelements(t)) elements" * (isempty(w) ? "" : " · " * w)
end

"""
Names of a target's components, in order: its parts (the drive's fitted pieces of a CAD
part are the part), then the pieces of the parts made in pieces (a hohlraum part's
`<part>_wall`, …), which are addressed on their own too.
"""
component_names(t::Target) = unique!([Briareus.parts(t); [c.name for c in t.components if _is_piece(c)]])

"""
Names of a target's opaque (non-aperture) components, as the checks list them: each part
(the drive's fitted pieces of a CAD part are the part), but each piece of a part made in
pieces (a hohlraum part's wall and caps).
"""
surface_names(t::Target) = unique!([_is_piece(c) ? c.name : c.part for c in t.components if c.kind === :surface])

# Whether `c` is a piece of a part made in pieces (named `<part>_<piece>`; a hohlraum
# part's), as opposed to a whole part or a fitted piece of a CAD part (`<part>/cyl1`, …).
_is_piece(c::Component) = c.part !== c.name && startswith(String(c.name), String(c.part) * "_")

# The component of `names` that a name in a diagnostic definition (written for a hohlraum's
# pieces: `leh_top`) refers to: itself when the target has it, else the first hohlraum
# part's piece of that name (`hohlraum_leh_top`), else itself.
function _piece_named(names, n::Symbol)
    n in names && return n
    k = findfirst(m -> endswith(String(m), "_" * String(n)), names)
    return k === nothing ? n : names[k]
end

# A diagnostic's view (`Briareus.diagnostic`) on a target with components `names`: the
# components it refers to (`through`, `center`, `area`) resolved with `_piece_named`, with
# a cache of its own; the view itself when it refers to none the target lacks.
function _view_on(v, names)
    T = typeof(v)
    refs = (:through, :center, :area)
    fs = fieldnames(T)
    any(f -> f in refs && getfield(v, f) isa Symbol && !(getfield(v, f) in names), fs) || return v
    vals = map(fs) do f
        x = getfield(v, f)
        f in refs && x isa Symbol ? _piece_named(names, x) :
        x isa IdDict ? typeof(x)() : x isa ReentrantLock ? ReentrantLock() : x
    end
    return T(vals...)
end
