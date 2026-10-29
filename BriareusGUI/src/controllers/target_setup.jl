# The target: a template (a hohlraum, a backlighter foil, or none) with its dimensions
# and wall properties, diagnostic windows in a hohlraum wall, free parts placed
# anywhere (spheres, disks, washers, rectangles, cylinders, cones, as VisRad's target
# components), and the mesh element sizes. Every edit rebuilds the `Target` (meshing
# takes milliseconds; the view factors are computed on request).

"""
Target templates: `:hohlraum` (cylindrical hohlraum with LEHs and an optional capsule,
`Briareus.hohlraum`), `:foil` (a backlighter foil facing TCC with an optional absorbing
sphere, `BLACK`, standing in for the main target) and `:custom` (parts only).
"""
const TARGET_KINDS = (:hohlraum, :foil, :custom)

const TARGET_KIND_LABELS = Dict(:hohlraum => "Hohlraum", :foil => "Backlighter foil",
                                :custom => "Custom (parts only)")

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
                (:wall_conversion, "Wall conversion", "laser → x-ray", 0.8)),
    foil = ((:distance, "Distance from TCC", "mm", 6.0),
            (:width, "Foil width", "mm", 2.0),
            (:conversion, "Conversion", "laser → x-ray", 0.02),
            (:main_radius, "TCC sphere radius", "mm (0: none)", 1.0)),
    custom = ())

"""
    MeshSpec(h, capsule_h)

Element sizes (mm): `h` for the hohlraum wall, end caps and LEHs (`hohlraum_counts`),
a foil target's main-target sphere and the parts; `capsule_h` for the surface whose
drive is studied: the capsule (`capsule_counts`) or the backlighter foil. A part may
override `h` with its own.
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
Closing entrance holes of a hohlraum (`Briareus.hohlraum(; leh_closure, leh_open_radius)`),
kept in the target's `params` with the dimensions: `:leh_closure` (1 on, 0 off; set
with [`set_leh_closure!`](@ref)) and these fields `(key, label, unit, default)`, edited
with [`edit_target!`](@ref). On, each entrance hole is open out to `leh_open_fraction ×`
the LEH radius and its lip beyond closes to x-rays at `leh_closure_time` over
`leh_closure_width` (a `GapClosure` whose closed albedo is the wall's). 0.84 is the
measured LEH radius at peak (research notes, "Closing entrance holes"); the time and
width are an engineering example (closure at the rise to peak), to be set per design.
"""
const LEH_CLOSURE_FIELDS = ((:leh_closure_time, "Closure time", "ns", 5.0),
                            (:leh_closure_width, "Closure duration", "ns", 0.5),
                            (:leh_open_fraction, "Open radius", "× LEH radius", 0.84))

# The closing lip of an entrance hole (`:leh_top_lip`, `:leh_bottom_lip`) is a split of
# the hole for the calculations: the GUI lists, selects, highlights and leaves it out
# with its hole.
const _LIP_NAME = r"^(leh_(?:top|bottom))_lip$"
"Whether `c` is a hohlraum's closing entrance-hole lip."
is_lip(c::Component) =
    c.kind === :aperture && c.material.albedo isa GapClosure && occursin(_LIP_NAME, String(c.name))
"The component a name is shown as: an entrance hole for its lip, else itself."
shown_component(name::Symbol) = (m = match(_LIP_NAME, String(name)); m === nothing ? name : Symbol(m[1]))
"Components `names` with the lips of the entrance holes among them."
with_lips(names) = Symbol[n for x in names for n in (x in (:leh_top, :leh_bottom) ? (x, Symbol(x, :_lip)) : (x,))]

"Name of the foil component of a `:foil` target; the main-target sphere is `:main`."
const FOIL_NAME = :backlighter

# ---------------------------------------------------------------- windows

"""
    WindowSpec(; phi = 0.0, z = 0.0, width = 1.0, height = 1.0, open = false,
               albedo = 0.85, conversion = 0.8, round = true)

A diagnostic window cut into the hohlraum wall (a core `WallWindow`): centered at
azimuth `phi` (degrees from +x) and axial position `z` (mm); `round` (the default, as
diagnostic windows usually are) a disk of diameter `width` mm on the unrolled wall,
otherwise a rectangle `width` mm of arc by `height` mm along the axis. A closed window
(`open = false`) is an opaque wall patch
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
`:cylinder` (open tube), `:cone` (open frustum) and `:mesh` (a triangle mesh from an STL
file, such as a CAD export).
"""
const PART_KINDS = (:sphere, :disk, :annulus, :rectangle, :cylinder, :cone, :mesh, :visrad)

const PART_KIND_LABELS = Dict(:sphere => "Sphere", :disk => "Disk", :annulus => "Washer (annulus)",
                              :rectangle => "Rectangle", :cylinder => "Cylinder (tube)",
                              :cone => "Cone (frustum)", :mesh => "Mesh (STL file)",
                              :visrad => "VisRad object")

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
    :visrad => ())

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
Frames of the target templates (VisRad's reference frames), which parts are placed in
and beams aim in: `:tcc` (the target chamber), `:leh_top` and `:leh_bottom` (a
hohlraum's entrance holes: origin at the hole's center, `z` pointing out of the
hohlraum, so one placement mirrors between the two ends), and `:backlighter` (a
backlighter target's foil, [`backlighter_frame`](@ref): origin at its center, `z` along
its front normal, toward TCC). A part in a template's frame moves with it when the
target's dimensions change. Every part has its own frame too, named after the part
([`part_frames`](@ref)): a part placed in it moves and turns with that part.
"""
const TEMPLATE_FRAMES = (:tcc, :leh_top, :leh_bottom, :backlighter)

const FRAME_LABELS = Dict(:tcc => "Target chamber (TCC)", :leh_top => "Top entrance hole",
                          :leh_bottom => "Bottom entrance hole", :backlighter => "Backlighter foil")

"""
    backlighter_frame(params) -> Frame

The frame of a backlighter target's foil (dimensions `params`, see
[`TARGET_FIELDS`](@ref)): origin at the foil's center ([`foil_center`](@ref)), `z` along
its front normal (toward TCC, chamber `-x`), `x` along chamber `z` and `y` along chamber
`y`, as the foil's mesh is laid out (its `width` edges along `x` and `y`). Local
`(u, v, 0)` is a point of the foil's face; `w < 0` lies behind it.
"""
backlighter_frame(p::AbstractDict) = Frame(rotation_to(Vec3(-1, 0, 0)), foil_center(p))

"""
    template_frames(kind, params) -> Dict{Symbol,Frame}

The [`TEMPLATE_FRAMES`](@ref) a target of template `kind` with dimensions `params` has:
`:tcc` always, the entrance holes' for a hohlraum, `:backlighter` for a foil target.
"""
function template_frames(kind::Symbol, p::AbstractDict)
    fr = Dict{Symbol,Frame}(:tcc => CHAMBER)
    if kind === :hohlraum
        L = p[:length]
        fr[:leh_top] = Frame([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0], Vec3(0, 0, L / 2))
        fr[:leh_bottom] = Frame([1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 -1.0], Vec3(0, 0, -L / 2))
    elseif kind === :foil
        fr[:backlighter] = backlighter_frame(p)
    end
    return fr
end

"""
    PartSpec(kind; name, values, two_sided = false, inward = false, aperture = false,
             file = "", frame = :tcc, coords = :cartesian, object = "", steps = [],
             foil = false)

A free target part: `kind` (one of [`PART_KINDS`](@ref)), its component `name`,
`values` (every [`PART_FIELDS`](@ref) key and its kind's [`PART_SIZE_FIELDS`](@ref)),
`two_sided` (both sides radiate and absorb), `inward` (a sphere, cylinder or cone
radiates from its inside; a flat part from its back; a mesh from the other side) and
`aperture` (a transparent diagnostic surface: rays pass, its incident flux is what
escapes through it). `foil` marks a part with the role of a backlighter foil
([`set_part_role!`](@ref), [`active_foil`](@ref)). A `:mesh` part reads its triangles from the STL `file`. The
placement fields are in the frame named `frame` (one of [`TEMPLATE_FRAMES`](@ref), or
another part's name: see [`part_frames`](@ref)); `coords` ([`PART_COORDINATES`](@ref)) is
how its position is shown and edited (the values keep `x`, `y`, `z`). The placement
`steps` (`PlacementStep`s: relations to other frames, relative moves) apply in order
after the placement fields; the degrees of freedom they leave free keep the fields'
values ([`part_placement`](@ref)).
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
                  object::AbstractString = "", steps = PlacementStep[], foil::Bool = false)
    kind in PART_KINDS || throw(ArgumentError("part kind must be one of $PART_KINDS"))
    frame === Symbol(name) && throw(ArgumentError("a part cannot be placed in its own frame"))
    coords in PART_COORDINATES || throw(ArgumentError("coordinates must be one of $PART_COORDINATES"))
    kind === :mesh && isempty(file) && throw(ArgumentError("a mesh part needs an STL file"))
    kind === :visrad && (isempty(file) || isempty(object)) &&
        throw(ArgumentError("a VisRad part needs a workspace file and an object name"))
    v = Dict{Symbol,Float64}(k => d for (k, _, _, d) in (PART_FIELDS..., PART_SIZE_FIELDS[kind]...))
    for (k, x) in values
        haskey(v, k) || throw(ArgumentError("a $kind has no field :$k"))
        v[k] = x
    end
    any(st -> Symbol(name) in Briareus._step_refs(st), steps) &&
        throw(ArgumentError("a part's placement step cannot refer to its own frame"))
    return PartSpec(kind, Symbol(name), v, two_sided, inward, aperture, String(file), frame,
                    coords, String(object), collect(PlacementStep, steps), foil)
end

_with(p::PartSpec; kind = p.kind, name = p.name, values = p.values, two_sided = p.two_sided,
      inward = p.inward, aperture = p.aperture, file = p.file, frame = p.frame,
      coords = p.coords, object = p.object, steps = p.steps, foil = p.foil) =
    PartSpec(kind, name, values, two_sided, inward, aperture, file, frame, coords, object, steps, foil)

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
    part_frames(kind, params, parts) -> Dict{Symbol,Frame}
    part_frames(ts::TargetSetup) -> Dict{Symbol,Frame}

Every named frame of a design in the chamber: the template's ([`template_frames`](@ref))
and each part's own ([`part_frame_def`](@ref), named after the part), resolved with
`resolve_frames`, so a part placed in another part's frame follows every ancestor
exactly. Throws `ArgumentError` for a template frame the template lacks, a part placed in
a frame that does not exist, or parts placed in each other in a loop (the error names
the loop).
"""
function part_frames(kind::Symbol, params::AbstractDict, parts::AbstractVector)
    base = template_frames(kind, params)
    names = Set{Symbol}(p.name for p in parts)
    for p in parts
        (haskey(base, p.frame) || p.frame in names) && continue
        p.frame in (:leh_top, :leh_bottom) &&
            throw(ArgumentError("only a hohlraum has entrance holes to place a part at"))
        p.frame === :backlighter &&
            throw(ArgumentError("only a backlighter target has a foil frame to place a part in"))
        throw(ArgumentError("no frame named $(p.frame) to place $(p.name) in"))
    end
    return resolve_frames([p.name => part_frame_def(p) for p in parts], base)
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
"""
function part_components(p::PartSpec, h::Real; frame::Frame = _own_frame(p), raw::Bool = false)
    # the same components (identical objects) for the same part, size and frame, so the
    # targets built from them are found again (`_target_of`): an edit elsewhere, an undo
    lock(() -> length(_PART_MEMO) > 512 && empty!(_PART_MEMO), _STL_LOCK)
    return _cached(() -> _part_components(p, h, frame, raw), _PART_MEMO, _part_key(p, h, frame, raw))
end

const _PART_MEMO = Dict{Any,Vector{Component}}()

function _part_key(p::PartSpec, h::Real, frame, raw::Bool)
    stamp = isfile(p.file) ? (mtime(p.file), filesize(p.file)) : nothing
    return (p.kind, p.name, sort!(collect(p.values)), p.two_sided, p.inward, p.aperture, p.file,
            p.object, stamp, Float64(h), frame.R, frame.origin, raw)
end

function _part_components(p::PartSpec, h::Real, frame::Frame, raw::Bool)
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
    out = !p.inward
    k = p.kind
    if k === :sphere
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
function _unmeshed_parts(kind, params, mesh::MeshSpec, parts, excluded)
    out = Symbol[]
    any(p -> p.kind in (:visrad, :mesh) && !(p.name in excluded), parts) || return out
    frames = try
        part_frames(kind, params, parts)
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
function _mesh_parts(parts, kind, params, mesh, raw, report, label; all = parts)
    n = length(parts)
    finished = Threads.Atomic{Int}(0)
    report("$label: meshing $n parts…")
    frames = part_frames(kind, params, all)
    tasks = map(parts) do p
        Threads.@spawn begin
            part_components(p, mesh.h; frame = frames[p.name], raw = p.name in raw)
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
    TargetSetup(; kind = :hohlraum, resolution = :standard, mesh = nothing,
                windows = WindowSpec[], parts = PartSpec[])

The target being designed: `kind` (one of [`TARGET_KINDS`](@ref)), `params` (every
[`TARGET_FIELDS`](@ref) value of every kind, so switching kinds keeps the others'
edits), `mesh` (the element sizes, a [`MeshSpec`](@ref)), `resolution` (the preset the
sizes match, one of [`MESH_RESOLUTIONS`](@ref), or `:custom`), `windows` (the
hohlraum's diagnostic windows, [`WindowSpec`](@ref)s; other templates ignore them),
`parts` (free parts, [`PartSpec`](@ref)s), `frames` (every named frame of the design
in the chamber, [`part_frames`](@ref): updated with the parts, before the target is
built), `renames` (while an edit updates `frames`, the parts it renamed, old name →
new, so that beams aimed in their frames follow; empty otherwise), the built `target`, `error` (the last
rejected edit), `excluded` (component names left out of the calculations, see
[`set_excluded!`](@ref)), `ghost` (those components, for drawing, or `nothing`) and
`active_foil` (the foil chosen with [`set_active_foil!`](@ref), or `nothing`; it follows
a rename, and [`active_foil`](@ref) gives the foil the analyses take).

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
(`parts`, `params`, …) are committed at once, `target` and `ghost` keep the previous
design's until the new ones are ready, and `building` is set meanwhile. A CAD part
given the role of the capsule has its sphere fitted on a worker as well (`fitting`
names it; the role is committed with the fit). The drive and the checks wait for both
([`target_ready`](@ref)).

Edit with [`set_target_kind!`](@ref), [`edit_target!`](@ref),
[`set_resolution!`](@ref), [`edit_mesh!`](@ref), the window functions
([`add_window!`](@ref), …) and the part functions ([`add_part!`](@ref), …).
"""
struct TargetSetup
    kind::Observable{Symbol}
    params::Observable{Dict{Symbol,Float64}}
    mesh::Observable{MeshSpec}
    resolution::Observable{Symbol}
    windows::Observable{Vector{WindowSpec}}
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

function TargetSetup(; kind::Symbol = :hohlraum, resolution::Symbol = :standard,
                     mesh::Union{Nothing,MeshSpec} = nothing,
                     windows::AbstractVector{WindowSpec} = WindowSpec[],
                     parts::AbstractVector{PartSpec} = PartSpec[])
    kind in TARGET_KINDS || throw(ArgumentError("target kind must be one of $TARGET_KINDS"))
    if mesh === nothing
        resolution in MESH_RESOLUTIONS ||
            throw(ArgumentError("resolution must be one of $MESH_RESOLUTIONS"))
        mesh = MESH_PRESETS[resolution]
    end
    params = Dict{Symbol,Float64}(k => v for fs in TARGET_FIELDS for (k, _, _, v) in fs)
    params[:leh_closure] = 0.0
    for (k, _, _, v) in LEH_CLOSURE_FIELDS
        params[k] = v
    end
    ws = collect(WindowSpec, windows)
    ps = collect(PartSpec, parts)
    kind === :custom && isempty(ps) && push!(ps, PartSpec(:sphere; name = :sphere_1))
    t = build_target(kind, params, mesh, ws, ps)
    return TargetSetup(Observable(kind), Observable(params), Observable(mesh),
                       Observable(_preset_of(mesh)), Observable(ws), Observable(ps),
                       Observable(t), Observable(""), Observable(Set{Symbol}()),
                       Observable{Union{Nothing,Target}}(nothing), Observable(Set{Symbol}()),
                       Ref(false), Ref{Any}(f -> f()), Ref(0), Observable(false),
                       Observable(Set{Symbol}()), Dict{Symbol,Int}(), Any[],
                       Observable(part_frames(kind, params, ps)), Ref(Dict{Symbol,Symbol}()),
                       Observable{Union{Nothing,Symbol}}(nothing))
end

part_frames(ts::TargetSetup) = ts.frames[]

"""
Targets of more triangles than this (≈ 30 ms to build) are built on a worker in a window
(see [`TargetSetup`](@ref)); smaller ones, such as the templates' meshes, at once. Set it
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
    build_target(kind, params, mesh, windows = WindowSpec[], parts = PartSpec[]) -> Target
    build_target(ts::TargetSetup) -> Target

The `Target` of a template kind, its dimensions (`params`, as in
[`TARGET_FIELDS`](@ref)), the element sizes (a [`MeshSpec`](@ref) or a preset name), a
hohlraum's diagnostic windows, and free parts. Throws `ArgumentError` for dimensions,
windows or parts the geometry rejects (overlapping windows, a window past the end of
the wall, a duplicate component name, …).
"""
build_target(kind::Symbol, p::AbstractDict, resolution::Symbol, windows = WindowSpec[],
             parts = PartSpec[]) =
    build_target(kind, p, MESH_PRESETS[resolution], windows, parts)

function build_target(kind::Symbol, p::AbstractDict, mesh::MeshSpec,
                      windows::AbstractVector{WindowSpec} = WindowSpec[],
                      parts::AbstractVector{PartSpec} = PartSpec[]; raw = ())
    return _target_of(_target_components(kind, p, mesh, windows, parts; raw))
end
build_target(ts::TargetSetup) = build_target(ts.kind[], ts.params[], ts.mesh[], ts.windows[], ts.parts[])

# The components of a design (`raw`: names of CAD parts drawn from their file's facets).
function _target_components(kind::Symbol, p::AbstractDict, mesh::MeshSpec,
                            windows::AbstractVector{WindowSpec} = WindowSpec[],
                            parts::AbstractVector{PartSpec} = PartSpec[]; raw = ())
    frac(k, what) = 0 <= p[k] <= 1 || throw(ArgumentError("$what must be between 0 and 1"))
    mesh.h > 0 && mesh.capsule_h > 0 || throw(ArgumentError("element sizes must be positive"))
    comps = Component[]
    if kind === :hohlraum
        frac(:wall_albedo, "the wall albedo")
        frac(:wall_conversion, "the wall conversion")
        p[:length] > 0 || throw(ArgumentError("the length must be positive"))
        p[:capsule_radius] >= 0 || throw(ArgumentError("the capsule radius must not be negative"))
        wall = gold(albedo = p[:wall_albedo], conversion = p[:wall_conversion])
        closing = get(p, :leh_closure, 0.0) != 0
        lip = closing ? _leh_closure(p, wall) : nothing
        fopen = closing ? p[:leh_open_fraction] : 0.84
        t = hohlraum(; length = p[:length], radius = p[:radius], leh_radius = p[:leh_radius],
                     capsule_radius = p[:capsule_radius], wall,
                     windows = [_wall_window(w, k) for (k, w) in enumerate(windows)],
                     h = mesh.h, capsule_h = mesh.capsule_h,
                     leh_closure = lip, leh_open_radius = fopen * p[:leh_radius])
        append!(comps, t.components)
    elseif kind === :foil
        frac(:conversion, "the conversion")
        p[:width] > 0 || throw(ArgumentError("the foil width must be positive"))
        p[:main_radius] >= 0 || throw(ArgumentError("the sphere radius must not be negative"))
        p[:distance] > p[:main_radius] + p[:width] / 2 ||
            throw(ArgumentError("the foil must sit outside the main-target sphere"))
        n = max(4, ceil(Int, p[:width] / mesh.capsule_h))
        # the foil's square laid out in its frame (`backlighter_frame`: front along local z)
        F = backlighter_frame(p)
        push!(comps, Component(FOIL_NAME, transform(rectangle_mesh(p[:width], p[:width]; nu = n, nv = n),
                                                    F.R, F.origin),
                               Material(:Ti; conversion = p[:conversion]); two_sided = false))
        if p[:main_radius] > 0
            (; capsule_ntheta, capsule_nphi) = capsule_counts(p[:main_radius], mesh.h)
            push!(comps, Component(:main, sphere_mesh(p[:main_radius]; ntheta = capsule_ntheta,
                                                     nphi = capsule_nphi), BLACK;
                                   surface = SphereSurface(p[:main_radius])))
        end
    elseif kind !== :custom
        throw(ArgumentError("unknown target kind :$kind"))
    end
    frames = isempty(parts) ? Dict{Symbol,Frame}() : part_frames(kind, p, parts)
    for part in parts
        new = part_components(part, mesh.h; frame = frames[part.name], raw = part.name in raw)
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

"""
    foil_center(params) -> Vec3

Center of the foil of a `:foil` target: `distance` mm from TCC along +x, facing TCC.
"""
foil_center(p::AbstractDict) = Vec3(p[:distance], 0.0, 0.0)
foil_center(ts::TargetSetup) = foil_center(ts.params[])

# A design's components split into the target of those in the calculations and that of
# those left out (`excluded`, drawn as a ghost: `nothing` when none is).
function _split_components(comps::Vector{Component}, excluded)
    keep, gh = _split_lists(comps, excluded)
    return _targets_of(keep, gh)
end

# The components in the calculations and those left out (`nothing` when none is).
function _split_lists(comps::Vector{Component}, excluded)
    out(c) = (is_lip(c) ? shown_component(c.name) : c.name) in excluded   # a lip with its hole
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
function _try_build!(ts::TargetSetup; kind = ts.kind[], params = ts.params[], mesh = ts.mesh[],
                     windows = ts.windows[], parts = ts.parts[], excluded = ts.excluded[],
                     prebuilt = nothing, sync::Bool = !ts.spawn[], renames = Dict{Symbol,Symbol}())
    try
        frames = part_frames(kind, params, parts)
        waiting = prebuilt === nothing && ts.spawn[] ?
                  _unmeshed_parts(kind, params, mesh, parts, excluded) : Symbol[]
        # components left out of the calculations (those still in the target)
        built, lists = if prebuilt === nothing
            keep, gh = _split_lists(_target_components(kind, params, mesh, windows, parts;
                                                       raw = union(Set{Symbol}(excluded), waiting)),
                                    excluded)
            !sync && _build_later(keep, gh) ? (nothing, (keep, gh)) : (_targets_of(keep, gh), nothing)
        else
            prebuilt, nothing
        end
        _assign!(ts.excluded, Set{Symbol}(excluded))
        built === nothing || (ts.ghost[] === built[2] || (ts.ghost[] = built[2]))
        # (first: what follows a frame, or a new kind's beams, sees the design's frames;
        # beams aimed in a renamed part's frame follow its new name)
        ts.renames[] = Dict{Symbol,Symbol}(renames)
        try
            _assign!(ts.frames, frames)
        finally
            ts.renames[] = Dict{Symbol,Symbol}()
        end
        kind == ts.kind[] || (ts.kind[] = kind)
        params === ts.params[] || (ts.params[] = params)
        if mesh != ts.mesh[]
            ts.mesh[] = mesh
            _assign!(ts.resolution, _preset_of(mesh))
        end
        windows === ts.windows[] || (ts.windows[] = windows)
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
        isempty(waiting) || _mesh_in_background!(ts, g, waiting, kind, params, mesh, parts)
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
function _mesh_in_background!(ts::TargetSetup, g::Integer, names, kind, params, mesh, parts)
    deliver = ts.deliver[]
    todo = [p for p in parts if p.name in names]
    job = () -> _mesh_parts(todo, kind, params, mesh, Set{Symbol}(), _ -> nothing, ""; all = parts)
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
    set_leh_closure!(ts::TargetSetup, on::Bool) -> Bool

Close a hohlraum's entrance holes during the pulse (`on`), with the settings of
[`LEH_CLOSURE_FIELDS`](@ref), or keep them open. Each hole is then split into the hole
and its closing lip (`:leh_top_lip`, …), shown, selected and left out with the hole.
"""
function set_leh_closure!(ts::TargetSetup, on::Bool)
    (get(ts.params[], :leh_closure, 0.0) != 0) == on && return true
    p = copy(ts.params[])
    p[:leh_closure] = on ? 1.0 : 0.0
    return _try_build!(ts; params = p)
end

"""
    set_target_kind!(ts::TargetSetup, kind::Symbol) -> Bool

Switch the template (one of [`TARGET_KINDS`](@ref)). Parts are kept; a custom target
without parts gets a sphere at TCC.
"""
function set_target_kind!(ts::TargetSetup, kind::Symbol)
    kind in TARGET_KINDS || throw(ArgumentError("target kind must be one of $TARGET_KINDS"))
    ts.kind[] == kind && return true
    parts = ts.parts[]
    kind === :custom && isempty(parts) && (parts = [PartSpec(:sphere; name = :sphere_1)])
    return _try_build!(ts; kind, parts)
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
        return _try_build!(ts; params = p)
    catch err
        err isa ArgumentError || rethrow()
        ts.error[] = _errmsg(err)
        return false
    end
end

function _target_label(key::Symbol)
    for fs in (TARGET_FIELDS..., LEH_CLOSURE_FIELDS), (k, label, _, _) in fs
        k === key && return label
    end
    return String(key)
end

# ---------------------------------------------------------------- windows

# Commit a new window list if the hohlraum accepts it; otherwise keep the target and
# report why. Windows only exist on a hohlraum.
function _set_windows!(ts::TargetSetup, ws::Vector{WindowSpec})
    if ts.kind[] !== :hohlraum
        ts.error[] = "diagnostic windows are cut into a hohlraum wall"
        return false
    end
    return _try_build!(ts; windows = ws)
end

function _check_window_index(ts::TargetSetup, i::Integer)
    1 <= i <= length(ts.windows[]) || throw(ArgumentError("no window $i"))
    return i
end

"""
    add_window!(ts::TargetSetup; kwargs...) -> Bool

Add a diagnostic window to the hohlraum wall: a [`WindowSpec`](@ref) from `kwargs`,
or by default a closed round window 1 mm across at the hohlraum's midplane with the wall's
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
            # the windows are checked as the wall is meshed (no target to build, no CAD
            # part to remesh: those are drawn from their facets here)
            ts.kind[] === :hohlraum && _target_components(ts.kind[], p, ts.mesh[], ws, ts.parts[];
                                                          raw = Set(q.name for q in ts.parts[]))
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
    ws[i] = WindowSpec(w.phi, w.z, w.width, w.height, open, w.albedo, w.conversion, w.round)
    return _set_windows!(ts, ws)
end

"""
    set_window_round!(ts::TargetSetup, i, round::Bool) -> Bool

Make window `i` a disk (diameter = its width) or a rectangle.
"""
function set_window_round!(ts::TargetSetup, i::Integer, round::Bool)
    _check_window_index(ts, i)
    w = ts.windows[][i]
    w.round == round && return true
    ws = copy(ts.windows[])
    ws[i] = WindowSpec(w.phi, w.z, w.width, w.height, w.open, w.albedo, w.conversion, round)
    return _set_windows!(ts, ws)
end

"Target component name of window `i`."
window_name(i::Integer) = Symbol(:window_, i)

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
    ps = copy(ts.parts[])
    deleteat!(ps, i)
    return _try_build!(ts; parts = ps)
end

"""
    edit_part!(ts::TargetSetup, i, key, text) -> Bool

Set a field of part `i` from text: a [`PART_FIELDS`](@ref) or size key, or `:name`.
"""
function edit_part!(ts::TargetSetup, i::Integer, key::Symbol, text)
    _check_part_index(ts, i)
    p = ts.parts[][i]
    try
        if key === :name
            s = strip(String(text))
            occursin(r"^[A-Za-z_][A-Za-z0-9_]*$", s) ||
                throw(ArgumentError("a part name is letters, digits and _, starting with a letter"))
            Symbol(s) == p.name && return true
            new = _with(p; name = Symbol(s))
            # parts placed in this one's frame, or with steps that refer to it, follow the rename
            ps = [_renamed(q, p.name, Symbol(s)) for q in ts.parts[]]
            ps[i] = new
            return _try_build!(ts; parts = ps, renames = Dict(p.name => new.name))
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

Place part `i` in another frame (one of [`TEMPLATE_FRAMES`](@ref), or another part's
name: see [`part_frames`](@ref)). With `keep` the part stays where it is in the
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
        frame === p.name && throw(ArgumentError("a part cannot be placed in its own frame"))
        kind, params, parts = ts.kind[], ts.params[], ts.parts[]
        trial = copy(parts)
        trial[i] = _with(p; frame)
        part_frames(kind, params, trial)                     # a loop, a missing frame
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
    return placement_report([p.name => part_frame_def(p) for p in ps],
                            template_frames(ts.kind[], ts.params[]))[ps[i].name]
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
           norm(part_frames(ts.kind[], ts.params[], ps)[p.name].origin - fr[p.name].origin) < 1e-12 * abs(Δ)
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
Roles a part can take in the analyses: `:none`, `:capsule` (the observed surface: its
x-ray drive's modes and the capsule keep-out; a sphere, named `:capsule` as the drive and
the optimizer look for it) and `:foil` (a backlighter foil, a flat part or a mesh that
keeps its name: any number of parts can be foils, and the active one
([`active_foil`](@ref)) is the foil whose spot statistics, map, objective, line of sight
and captured power the analyses take).
"""
const PART_ROLES = (:none, :capsule, :foil)

const PART_ROLE_LABELS = Dict(:none => "None", :capsule => "Capsule (observed x-ray drive)",
                              :foil => "Backlighter foil (spot uniformity)")

"Role of a part (one of [`PART_ROLES`](@ref)): the capsule by its name, a foil by its flag."
part_role(p::PartSpec) = p.name === :capsule ? :capsule : p.foil ? :foil : :none

"""
    foil_names(ts::TargetSetup) -> Vector{Symbol}

The backlighter foils of the design in the calculations: a foil target's foil
(`FOIL_NAME`), then each part with the role `:foil` in the parts' order. Foils left out of
the calculations are not listed.
"""
function foil_names(ts::TargetSetup)
    out = Symbol[]
    ex = ts.excluded[]
    ts.kind[] === :foil && !(FOIL_NAME in ex) && push!(out, FOIL_NAME)
    for p in ts.parts[]
        p.foil && !(p.name in ex) && push!(out, p.name)
    end
    return out
end

"""
    active_foil(ts::TargetSetup) -> Union{Symbol,Nothing}

The foil the foil analyses take (the drive's spot statistics and foil map, the foil
objective, the sweep's line of sight and captured power): the one chosen with
[`set_active_foil!`](@ref) while it is still a foil in the calculations, else the first
of [`foil_names`](@ref) (a foil target's foil, then the first foil part), else `nothing`.
A part given the foil role when the design has no foil becomes the chosen one.
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

The components of `t` that belong to a backlighter foil, active or not: a foil target's
foil and the components of each part with the foil role (a VisRad part's carry its name
as a prefix). The clearance and TaLIS checks take them as surfaces the beams land on.
"""
function foil_components(ts::TargetSetup, t::Target = ts.target[])
    fs = Set{Symbol}(p.name for p in ts.parts[] if p.foil)
    vr = Set{Symbol}(p.name for p in ts.parts[] if p.foil && p.kind === :visrad)
    out = Symbol[]
    for c in t.components
        n = c.name
        if (ts.kind[] === :foil && n === FOIL_NAME) || n in fs || c.part in fs ||
           any(v -> startswith(String(n), String(v)), vr)
            push!(out, n)
        end
    end
    return out
end

# The object kind a part stands for (a VisRad object's own kind).
function _part_shape(p::PartSpec)
    p.kind === :visrad || return p.kind
    k = _visrad_object(p.file, p.object).kind
    return k === :mesh ? :mesh : k
end

# A part's name without a role: a VisRad object's own name, else its kind and a free number.
function _plain_name(p::PartSpec, taken)
    p.kind === :visrad && return Symbol(p.object)
    k = 1
    while Symbol(p.kind, :_, k) in taken
        k += 1
    end
    return Symbol(p.kind, :_, k)
end

"""
    set_part_role!(ts::TargetSetup, i, role::Symbol; done = () -> nothing) -> Bool

Give part `i` a role in the analyses (one of [`PART_ROLES`](@ref)). `:foil` (a disk,
rectangle or mesh, or a VisRad rectangle or CAD part: a backlighter whose spot the foil
objective makes uniform) marks the part as a foil and keeps its name; other foils keep
theirs, and the part becomes the active foil ([`active_foil`](@ref)) when the design had
none. `:capsule` names the part as the drive and the optimizer expect. The capsule is a
sphere (the mode analysis and the keep-out take its center and radius): a sphere part or
a VisRad sphere is renamed; a CAD part (a capsule with its fill tube, say) is replaced
in the calculations by a sphere part `:capsule` fitted to it ([`fit_sphere`](@ref): its
element size the target's capsule size, its albedo and conversion the part's), and left
out of the calculations itself. A part that held the capsule role gives it up (back to
its own name). `:none` clears the foil role, or renames the capsule to its own name (a
VisRad object's, or its kind and a number). Parts placed in a renamed part's frame
follow. A role the part cannot take, or the capsule role a hohlraum's capsule already
holds, sets `error` and returns `false`. In a window
(`spawn[]`) a large CAD capsule is fitted on a worker (see [`TargetSetup`](@ref)):
`true` then means the fit is under way, and `done()` runs once the role is given (at
once otherwise).
"""
function set_part_role!(ts::TargetSetup, i::Integer, role::Symbol; done = () -> nothing)
    _check_part_index(ts, i)
    role in PART_ROLES || throw(ArgumentError("a part's role must be one of $PART_ROLES"))
    ps = copy(ts.parts[])
    p = ps[i]
    part_role(p) === role && (done(); return true)
    ex = copy(ts.excluded[])
    renames = Dict{Symbol,Symbol}()
    try
        shape = _part_shape(p)
        if role === :capsule && shape !== :sphere
            shape === :mesh ||
                throw(ArgumentError("the capsule is a sphere, or a CAD part a sphere is fitted to"))
            return _capsule_from_part!(ts, i; done)
        end
        role === :foil && !(shape in (:rectangle, :disk, :mesh)) &&
            throw(ArgumentError("the backlighter foil is a flat part (a rectangle or disk) or a mesh"))
        taken = Set{Symbol}(q.name for q in ps)
        function rename(k, new)
            old = ps[k].name
            renames[old] = new
            ps[k] = _with(ps[k]; name = new)
            for j in eachindex(ps)
                ps[j] = _renamed(ps[j], old, new)
            end
            old in ex && (delete!(ex, old); push!(ex, new))    # stays left out
        end
        # the foil role is a flag; the capsule's is the name `:capsule`
        ps[i] = _with(ps[i]; foil = role === :foil)
        if role === :capsule
            j = findfirst(q -> q.name === :capsule, ps)
            j === nothing || rename(j, _plain_name(ps[j], taken))
            rename(i, :capsule)
        elseif p.name === :capsule
            rename(i, _plain_name(p, taken))
        end
    catch err
        err isa ArgumentError || rethrow()
        ts.error[] = _errmsg(err)
        return false
    end
    first_foil = role === :foil && active_foil(ts) === nothing
    ok = _try_build!(ts; parts = ps, excluded = ex, renames)
    ok && first_foil && (ts.active_foil[] = ps[i].name)
    ok && done()
    return ok
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

# Role `:capsule` for CAD part `i`: a sphere part `:capsule` fitted to it in its place in
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
    ps = copy(ts.parts[])
    p = ps[i]
    ex = copy(ts.excluded[])
    taken = Set{Symbol}(q.name for q in ps)
    j = findfirst(q -> q.name === :capsule, ps)
    renames = Dict{Symbol,Symbol}()
    if j !== nothing
        old = ps[j].name
        ps[j] = _with(ps[j]; name = _plain_name(ps[j], taken))
        renames[old] = ps[j].name
        for k in eachindex(ps)
            ps[k].frame === old && (ps[k] = _with(ps[k]; frame = ps[j].name))
        end
        old in ex && (delete!(ex, old); push!(ex, ps[j].name))
    end
    push!(ps, PartSpec(:sphere; name = :capsule,
                       values = Dict(:x => c[1], :y => c[2], :z => c[3], :r1 => r,
                                     :h => ts.mesh[].capsule_h, :albedo => p.values[:albedo],
                                     :conversion => p.values[:conversion])))
    push!(ex, p.name)
    return _try_build!(ts; parts = ps, excluded = ex, renames)
end

"""
    part_index(ts::TargetSetup, name::Symbol) -> Union{Nothing,Int}
"""
part_index(ts::TargetSetup, name::Symbol) = findfirst(p -> p.name === name, ts.parts[])

# ---------------------------------------------------------------- summary

"""
    target_summary(ts::TargetSetup) -> String

Template, windows and parts, and the element count.
"""
function target_summary(ts::TargetSetup)
    t = ts.target[]
    kind = ts.kind[] === :hohlraum ? "hohlraum" : ts.kind[] === :foil ? "foil" : "custom target"
    nw = ts.kind[] === :hohlraum ? length(ts.windows[]) : 0
    nw > 0 && (kind *= nw == 1 ? " with 1 window" : " with $nw windows")
    np = length(ts.parts[])
    np > 0 && ts.kind[] !== :custom && (kind *= np == 1 ? " + 1 part" : " + $np parts")
    ts.kind[] === :custom && (kind *= np == 1 ? " (1 part)" : " ($np parts)")
    w = target_waiting(ts)
    return "$kind · $(nelements(t)) elements" * (isempty(w) ? "" : " · " * w)
end

"Names of a target's components, in order: its parts (the drive's fitted pieces of a CAD part are the part)."
component_names(t::Target) = Briareus.parts(t)

"Names of a target's opaque (non-aperture) components, as parts (see `component_names`)."
surface_names(t::Target) = unique!([c.part for c in t.components if c.kind === :surface])
