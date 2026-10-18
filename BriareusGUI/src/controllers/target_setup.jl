# The target: a template (a hohlraum, a backlighter foil, or none) with its dimensions
# and wall properties, diagnostic windows in a hohlraum wall, free parts placed
# anywhere (spheres, disks, washers, rectangles, cylinders, cones, as VisRad's target
# components), and the mesh element sizes. Every edit rebuilds the `Target` (meshing
# takes milliseconds; the view factors are computed on request).

"""
Target templates: `:hohlraum` (cylindrical hohlraum with LEHs and an optional capsule,
`Briareus.hohlraum`), `:foil` (a backlighter foil facing TCC with an optional black
sphere standing in for the main target) and `:custom` (parts only).
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
                (:wall_conversion, "Wall laser → x-ray conversion", "", 0.8)),
    foil = ((:distance, "Distance from TCC", "mm", 6.0),
            (:width, "Foil width", "mm", 2.0),
            (:conversion, "Laser → x-ray conversion", "", 0.02),
            (:main_radius, "Main-target sphere radius", "mm", 1.0)),
    custom = ())

"""
    MeshSpec(h, capsule_h)

Element sizes (mm): `h` for the hohlraum wall, end caps and LEHs (`hohlraum_counts`),
a foil target's black main sphere and the parts; `capsule_h` for the surface whose
drive is studied: the capsule (`capsule_counts`) or the backlighter foil. A part may
override `h` with its own.
"""
struct MeshSpec
    h::Float64
    capsule_h::Float64
end

"""
Mesh presets (`MeshSpec`s), from the work-precision study in
`examples/mesh_convergence.jl` (reference/research-notes.md, "Hohlraum mesh
work-precision") on the NIF hohlraum: errors in capsule P2/P0 and P4/P0 against the
converged meshes of about 0.5 percentage points (`:coarse`, ≈ 1,300 elements, view
factors in under a second), 0.1 (`:standard`, ≈ 5,900 elements, a few seconds) and
0.02 (`:fine`, ≈ 16,000 elements, tens of seconds: converged for design work). The
capsule mesh sets P4: its error falls as `capsule_h²`.
"""
const MESH_PRESETS = (coarse = MeshSpec(0.6, 0.2), standard = MeshSpec(0.28, 0.09),
                      fine = MeshSpec(0.17, 0.05))

"Mesh presets in order; `:custom` marks sizes edited by hand."
const MESH_RESOLUTIONS = (:coarse, :standard, :fine)

_preset_of(m::MeshSpec) = something(findfirst(==(m), MESH_PRESETS), :custom)

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
    :mesh => ((:scale, "Scale (mm per file unit)", "", 1.0),),
    :visrad => ())

"""
Placement and material fields of every part: `(key, label, unit, default)`. The part
is rotated by `rotation` about its own axis, then tilted to polar angle `theta` and
azimuth `phi` (VisRad's Polar/Azimuthal/Rotation order: the axis then points along
`(θ, φ)`), then moved to `(x, y, z)`. `h` (mm) overrides the mesh element size (0:
the target's).
"""
const PART_FIELDS = ((:x, "x", "mm", 0.0), (:y, "y", "mm", 0.0), (:z, "z", "mm", 0.0),
                     (:theta, "Axis polar θ", "°", 0.0), (:phi, "Axis azimuth φ", "°", 0.0),
                     (:rotation, "Rotation about axis", "°", 0.0),
                     (:albedo, "Albedo", "", 0.0), (:conversion, "Laser → x-ray conversion", "", 0.0),
                     (:h, "Element size (0: target's)", "mm", 0.0))

"""
Frames a part can be placed in (VisRad's reference frames): `:tcc` (the target chamber),
`:leh_top` and `:leh_bottom` (a hohlraum's entrance holes: origin at the hole's center,
`z` pointing out of the hohlraum, so one placement mirrors between the two ends), and
`:foil` (a backlighter foil's center, `z` along its back normal, away from TCC). A part
in a target's frame moves with it when its dimensions change. A part can also be
placed in another part's own frame (its anchor is that part's name): it then moves and
turns with that part.
"""
const PART_ANCHORS = (:tcc, :leh_top, :leh_bottom, :foil)

const PART_ANCHOR_LABELS = Dict(:tcc => "Target chamber (TCC)", :leh_top => "Top entrance hole",
                                :leh_bottom => "Bottom entrance hole", :foil => "Backlighter foil")

"""
    anchor_frame(anchor, kind, params, parts = PartSpec[]) -> (R, origin)

Rotation and origin of a part's frame for a target of template `kind` with dimensions
`params`: one of [`PART_ANCHORS`](@ref), or the name of one of `parts` (that part's own
frame: its axis is the frame's `z`). A point `x` of the frame is `R * x + origin` in the
chamber. Throws `ArgumentError` for a frame the template lacks, an unknown part, or
parts placed relative to each other in a loop.
"""
function anchor_frame(anchor::Symbol, kind::Symbol, params::AbstractDict,
                      parts::AbstractVector = PartSpec[]; _visiting = Set{Symbol}())
    I3 = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    if !(anchor in PART_ANCHORS)
        i = findfirst(q -> q.name === anchor, parts)
        i === nothing && throw(ArgumentError("no part named $anchor to place a part in"))
        anchor in _visiting && throw(ArgumentError("parts placed relative to each other in a loop ($anchor)"))
        push!(_visiting, anchor)
        q = parts[i]
        Rq, tq = anchor_frame(q.anchor, kind, params, parts; _visiting)
        delete!(_visiting, anchor)
        return Rq * _part_rotation(q), Rq * part_center(q) + tq
    end
    anchor === :tcc && return I3, zero(Vec3)
    if anchor === :leh_top || anchor === :leh_bottom
        kind === :hohlraum || throw(ArgumentError("only a hohlraum has entrance holes to place a part at"))
        L = params[:length]
        return anchor === :leh_top ? (I3, Vec3(0, 0, L / 2)) :
               ([1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 -1.0], Vec3(0, 0, -L / 2))
    end
    anchor === :foil || throw(ArgumentError("a part's frame must be one of $PART_ANCHORS"))
    kind === :foil || throw(ArgumentError("only a backlighter target has a foil to place a part at"))
    return [0.0 0.0 1.0; 0.0 1.0 0.0; -1.0 0.0 0.0], foil_center(params)
end

"""
    PartSpec(kind; name, values, two_sided = false, inward = false, aperture = false,
             file = "", anchor = :tcc)

A free target part: `kind` (one of [`PART_KINDS`](@ref)), its component `name`,
`values` (every [`PART_FIELDS`](@ref) key and its kind's [`PART_SIZE_FIELDS`](@ref)),
`two_sided` (both sides radiate and absorb), `inward` (a sphere, cylinder or cone
radiates from its inside; a flat part from its back; a mesh from the other side) and
`aperture` (a transparent diagnostic surface: rays pass, its incident flux is what
escapes through it). A `:mesh` part reads its triangles from the STL `file`. The
placement fields are in the frame `anchor` ([`PART_ANCHORS`](@ref), or another part's
name); `coords` ([`PART_COORDINATES`](@ref)) is how its position is shown and edited
(the values keep `x`, `y`, `z`).
"""
struct PartSpec
    kind::Symbol
    name::Symbol
    values::Dict{Symbol,Float64}
    two_sided::Bool
    inward::Bool
    aperture::Bool
    file::String
    anchor::Symbol
    coords::Symbol
    object::String
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

"Position of a part in its coordinates: `key => value` for its position fields."
function part_position(p::PartSpec)
    x, y, z = p.values[:x], p.values[:y], p.values[:z]
    p.coords === :cartesian && return [:x => x, :y => y, :z => z]
    φ = mod(atand(y, x), 360)
    p.coords === :cylindrical && return [:pos_r => hypot(x, y), :pos_phi => φ, :z => z]
    r = sqrt(x^2 + y^2 + z^2)
    return [:pos_r => r, :pos_theta => r > 0 ? acosd(clamp(z / r, -1, 1)) : 0.0, :pos_phi => φ]
end

# x, y, z from a part's position in its coordinates with `key` set to `v`.
function _set_position(p::PartSpec, key::Symbol, v::Real)
    pos = Dict(part_position(p))
    haskey(pos, key) || throw(ArgumentError("$(p.coords) coordinates have no field :$key"))
    pos[key] = v
    p.coords === :cartesian && return pos[:x], pos[:y], pos[:z]
    if p.coords === :cylindrical
        pos[:pos_r] >= 0 || throw(ArgumentError("the radius must not be negative"))
        return pos[:pos_r] * cosd(pos[:pos_phi]), pos[:pos_r] * sind(pos[:pos_phi]), pos[:z]
    end
    r, θ, φ = pos[:pos_r], pos[:pos_theta], pos[:pos_phi]
    r >= 0 || throw(ArgumentError("the radius must not be negative"))
    return r * sind(θ) * cosd(φ), r * sind(θ) * sind(φ), r * cosd(θ)
end

function PartSpec(kind::Symbol; name = kind, values = Dict{Symbol,Float64}(),
                  two_sided::Bool = false, inward::Bool = false, aperture::Bool = false,
                  file::AbstractString = "", anchor::Symbol = :tcc, coords::Symbol = :cartesian,
                  object::AbstractString = "")
    kind in PART_KINDS || throw(ArgumentError("part kind must be one of $PART_KINDS"))
    anchor === Symbol(name) && throw(ArgumentError("a part cannot be placed in its own frame"))
    coords in PART_COORDINATES || throw(ArgumentError("coordinates must be one of $PART_COORDINATES"))
    kind === :mesh && isempty(file) && throw(ArgumentError("a mesh part needs an STL file"))
    kind === :visrad && (isempty(file) || isempty(object)) &&
        throw(ArgumentError("a VisRad part needs a workspace file and an object name"))
    v = Dict{Symbol,Float64}(k => d for (k, _, _, d) in (PART_FIELDS..., PART_SIZE_FIELDS[kind]...))
    for (k, x) in values
        haskey(v, k) || throw(ArgumentError("a $kind has no field :$k"))
        v[k] = x
    end
    return PartSpec(kind, Symbol(name), v, two_sided, inward, aperture, String(file), anchor,
                    coords, String(object))
end

_with(p::PartSpec; kind = p.kind, name = p.name, values = p.values, two_sided = p.two_sided,
      inward = p.inward, aperture = p.aperture, file = p.file, anchor = p.anchor,
      coords = p.coords, object = p.object) =
    PartSpec(kind, name, values, two_sided, inward, aperture, file, anchor, coords, object)

"Axis (unit vector) of a part in its frame: polar `theta` and azimuth `phi` in degrees."
part_axis(p::PartSpec) = Vec3(sind(p.values[:theta]) * cosd(p.values[:phi]),
                              sind(p.values[:theta]) * sind(p.values[:phi]), cosd(p.values[:theta]))

"Center of a part in its frame."
part_center(p::PartSpec) = Vec3(p.values[:x], p.values[:y], p.values[:z])

# Rz(φ) Ry(θ) Rz(ψ): VisRad's order (rotation about the own axis, polar tilt, azimuth).
function _part_rotation(p::PartSpec)
    θ, φ, ψ = deg2rad.((p.values[:theta], p.values[:phi], p.values[:rotation]))
    Rz(a) = [cos(a) -sin(a) 0; sin(a) cos(a) 0; 0 0 1]
    Ry(a) = [cos(a) 0 sin(a); 0 1 0; -sin(a) 0 cos(a)]
    return Rz(φ) * Ry(θ) * Rz(ψ)
end

_count(len, h; min = 1, multiple = 1) = multiple * max(cld(min, multiple), ceil(Int, len / h / multiple))

"""
    part_component(p::PartSpec, h; frame = (I, 0)) -> Component

The core component of a part, meshed with elements of size `h` (mm) unless the part
gives its own, carrying its exact curved surface (sphere, disk, washer, cylinder) for
the view factors, and placed in the chamber through its `frame` (`(R, origin)`, see
[`anchor_frame`](@ref)).
"""
function part_component(p::PartSpec, h::Real;
                        frame = ([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0], zero(Vec3)))
    p.kind === :visrad && throw(ArgumentError("a VisRad object can make several components: use part_components"))
    return only(part_components(p, h; frame))
end

"""
    part_components(p::PartSpec, h; frame = (I, 0)) -> Vector{Component}

The core components of a part: one for every kind but `:visrad`, whose object can be a
primitive hohlraum (wall, caps and entrance holes). A VisRad object keeps its own place
in the chamber (VisRad's frames), moved further by the part's placement; its windings
are repaired and its CAD surface remeshed into elements of the part's size.
"""
function part_components(p::PartSpec, h::Real;
                         frame = ([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0], zero(Vec3)))
    Rf, tf = frame
    R = Rf * _part_rotation(p)
    center = Rf * part_center(p) + tf
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
        cs = _visrad_components(p, hp, transpose(R) * (zero(Vec3) - center))
        return [transform(c, R, center) for c in cs]
    elseif k === :mesh
        v[:scale] > 0 || throw(ArgumentError("$(p.name): the scale must be positive"))
        # facing the target center, remeshed into elements of the part's size
        toward = transpose(R) * (zero(Vec3) - center)
        m = _stl_part_mesh(p.file, v[:scale], hp, toward, p.inward)
        surf = nothing
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
# CAD surfaces remeshed into elements of `h`), kept per workspace version and settings: every edit
# rebuilds the target, and orienting a CAD part of 10⁵ facets takes seconds.
const _VRW_PARTS = Dict{Any,Vector{Component}}()

function _visrad_components(p::PartSpec, h, toward)
    v = p.values
    key = (abspath(p.file), mtime(p.file), p.object, Float64(h), round.(Tuple(toward); digits = 9),
           p.two_sided, v[:albedo], v[:conversion], p.name)
    lock(_STL_LOCK) do
        get!(_VRW_PARTS, key) do
            o = _visrad_object(p.file, p.object)
            mat = Material(p.name; albedo = v[:albedo], conversion = v[:conversion])
            cs = visrad_component(o; material = mat, two_sided = p.two_sided, orient = !p.two_sided,
                                  toward, remesh = o.kind === :mesh ? h : nothing)
            isempty(cs) && throw(ArgumentError("$(p.name): VisRad object type $(o.type) is not supported"))
            cs
        end
    end
end

"Names of the components a VisRad object makes (without meshing its CAD)."
visrad_component_names(o::VisRadObject) =
    o.kind === :mesh ? [Symbol(o.name)] : [c.name for c in visrad_component(o)]

# A VisRad workspace's object, the workspace read once per file version.
const _VRW_CACHE = Dict{Any,VisRadWorkspace}()

function _visrad_object(file::AbstractString, object::AbstractString)
    isfile(file) || throw(ArgumentError("VisRad workspace not found: $file"))
    key = (abspath(file), mtime(file), filesize(file))
    ws = lock(_STL_LOCK) do
        get!(() -> read_visrad_workspace(file), _VRW_CACHE, key)
    end
    i = findfirst(o -> o.name == object, ws.objects)
    i === nothing && throw(ArgumentError("no object $object in $(basename(file))"))
    return ws.objects[i]
end

# An STL file's mesh, oriented (`orient_mesh`: closed solids face out of their material,
# open sheets toward `toward`, flipped for `inward`) and remeshed into elements of `h`.
# Every edit rebuilds the target, so the result is kept per file version and settings.
const _STL_CACHE = Dict{Any,Mesh}()
const _STL_LOCK = ReentrantLock()

function _stl_part_mesh(file::AbstractString, scale, h, toward, inward::Bool)
    isfile(file) || throw(ArgumentError("STL file not found: $file"))
    key = (abspath(file), mtime(file), filesize(file), Float64(scale), Float64(h),
           round.(Tuple(toward); digits = 9), inward)
    lock(_STL_LOCK) do
        get!(_STL_CACHE, key) do
            m = remesh_mesh(orient_mesh(read_stl(file; scale); toward); size = h)
            inward ? Mesh(m.vertices, [(a, c, b) for (a, b, c) in m.triangles], m.patches) : m
        end
    end
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
`parts` (free parts, [`PartSpec`](@ref)s), the built `target`, `error` (the last
rejected edit), `excluded` (component names left out of the calculations, see
[`set_excluded!`](@ref)) and `ghost` (those components, for drawing, or `nothing`).
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
    ws = collect(WindowSpec, windows)
    ps = collect(PartSpec, parts)
    kind === :custom && isempty(ps) && push!(ps, PartSpec(:sphere; name = :sphere_1))
    t = build_target(kind, params, mesh, ws, ps)
    return TargetSetup(Observable(kind), Observable(params), Observable(mesh),
                       Observable(_preset_of(mesh)), Observable(ws), Observable(ps),
                       Observable(t), Observable(""), Observable(Set{Symbol}()),
                       Observable{Union{Nothing,Target}}(nothing))
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
                      parts::AbstractVector{PartSpec} = PartSpec[])
    frac(k, what) = 0 <= p[k] <= 1 || throw(ArgumentError("$what must be between 0 and 1"))
    mesh.h > 0 && mesh.capsule_h > 0 || throw(ArgumentError("element sizes must be positive"))
    comps = Component[]
    if kind === :hohlraum
        frac(:wall_albedo, "the wall albedo")
        frac(:wall_conversion, "the wall conversion")
        p[:length] > 0 || throw(ArgumentError("the length must be positive"))
        p[:capsule_radius] >= 0 || throw(ArgumentError("the capsule radius must not be negative"))
        t = hohlraum(; length = p[:length], radius = p[:radius], leh_radius = p[:leh_radius],
                     capsule_radius = p[:capsule_radius],
                     wall = gold(albedo = p[:wall_albedo], conversion = p[:wall_conversion]),
                     windows = [_wall_window(w, k) for (k, w) in enumerate(windows)],
                     h = mesh.h, capsule_h = mesh.capsule_h)
        append!(comps, t.components)
    elseif kind === :foil
        frac(:conversion, "the conversion")
        p[:width] > 0 || throw(ArgumentError("the foil width must be positive"))
        p[:main_radius] >= 0 || throw(ArgumentError("the sphere radius must not be negative"))
        p[:distance] > p[:main_radius] + p[:width] / 2 ||
            throw(ArgumentError("the foil must sit outside the main-target sphere"))
        n = max(4, ceil(Int, p[:width] / mesh.capsule_h))
        push!(comps, foil(FOIL_NAME, foil_center(p), (-1, 0, 0); width = p[:width], n,
                          material = Material(:Ti; conversion = p[:conversion]), two_sided = false))
        if p[:main_radius] > 0
            (; capsule_ntheta, capsule_nphi) = capsule_counts(p[:main_radius], mesh.h)
            push!(comps, Component(:main, sphere_mesh(p[:main_radius]; ntheta = capsule_ntheta,
                                                     nphi = capsule_nphi), BLACK;
                                   surface = SphereSurface(p[:main_radius])))
        end
    elseif kind !== :custom
        throw(ArgumentError("unknown target kind :$kind"))
    end
    for part in parts
        new = part_components(part, mesh.h; frame = anchor_frame(part.anchor, kind, p, parts))
        for c in new
            any(d -> d.name === c.name, comps) && throw(ArgumentError("two components are named $(c.name)"))
        end
        append!(comps, new)
    end
    isempty(comps) && throw(ArgumentError("the target is empty: add a part"))
    return Target(comps)
end
build_target(ts::TargetSetup) = build_target(ts.kind[], ts.params[], ts.mesh[], ts.windows[], ts.parts[])

"""
    foil_center(params) -> Vec3

Center of the foil of a `:foil` target: `distance` mm from TCC along +x, facing TCC.
"""
foil_center(p::AbstractDict) = Vec3(p[:distance], 0.0, 0.0)
foil_center(ts::TargetSetup) = foil_center(ts.params[])

# Rebuild with changed inputs, committing them only if the geometry accepts them.
function _try_build!(ts::TargetSetup; kind = ts.kind[], params = ts.params[], mesh = ts.mesh[],
                     windows = ts.windows[], parts = ts.parts[], excluded = ts.excluded[])
    try
        t = build_target(kind, params, mesh, windows, parts)
        # components left out of the calculations (those still in the target)
        out = Set(c.name for c in t.components if c.name in excluded)
        ghost = nothing
        if !isempty(out)
            keep = [c for c in t.components if !(c.name in out)]
            isempty(keep) && throw(ArgumentError("every component is left out of the calculations"))
            ghost = Target([c for c in t.components if c.name in out])
            t = Target(keep)
        end
        _assign!(ts.excluded, Set{Symbol}(excluded))
        ts.ghost[] = ghost
        kind == ts.kind[] || (ts.kind[] = kind)
        params === ts.params[] || (ts.params[] = params)
        if mesh != ts.mesh[]
            ts.mesh[] = mesh
            _assign!(ts.resolution, _preset_of(mesh))
        end
        windows === ts.windows[] || (ts.windows[] = windows)
        parts === ts.parts[] || (ts.parts[] = parts)
        ts.target[] = t
        ts.error[] = ""
        return true
    catch err
        err isa ArgumentError || rethrow()
        ts.error[] = _errmsg(err)
        return false
    end
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
    for fs in TARGET_FIELDS, (k, label, _, _) in fs
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
            ts.kind[] === :hohlraum && build_target(ts.kind[], p, ts.mesh[], ws, ts.parts[])
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
    add_part!(ts::TargetSetup, kind::Symbol; name, file = "", values...) -> Bool

Add a free part (one of [`PART_KINDS`](@ref)) named `name` (default: the kind and the
first free number) with [`PART_FIELDS`](@ref) and size values given as keywords
(e.g. `add_part!(ts, :disk; z = 6.0, r1 = 1.0, theta = 180)`). A `:mesh` part reads the
STL `file` (scaled by `scale`, mm per file unit), placed like any part, with its
windings repaired (closed solids facing out of their material, open sheets toward the
target center) and its surface remeshed into elements of the part's size (`h`, or the
target's; `remesh_mesh`). A part the geometry rejects sets `error` and returns `false`.
"""
function add_part!(ts::TargetSetup, kind::Symbol; name = nothing, two_sided::Bool = false,
                   inward::Bool = false, aperture::Bool = false, file::AbstractString = "",
                   anchor::Symbol = :tcc, values...)
    taken = Set(component_names(ts.target[]))
    if name === nothing
        k = 1
        while Symbol(kind, :_, k) in taken
            k += 1
        end
        name = Symbol(kind, :_, k)
    end
    part = try
        PartSpec(kind; name, values = Dict{Symbol,Float64}(k => Float64(v) for (k, v) in values),
                 two_sided, inward, aperture, file, anchor)
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
            # parts placed in this one's frame follow the rename
            ps = [q.anchor === p.name ? _with(q; anchor = Symbol(s)) : q for q in ts.parts[]]
            ps[i] = new
            return _try_build!(ts; parts = ps)
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
    set_part_anchor!(ts::TargetSetup, i, anchor::Symbol; keep = true) -> Bool

Place part `i` in another frame ([`PART_ANCHORS`](@ref), or another part's name). With
`keep` the part stays
where it is in the chamber: its center and axis are converted to the new frame (and
its rotation about the axis kept, exact when the frames differ by a translation or a
half turn); otherwise its placement values are kept and it moves with the frame.
"""
function set_part_anchor!(ts::TargetSetup, i::Integer, anchor::Symbol; keep::Bool = true)
    _check_part_index(ts, i)
    p = ts.parts[][i]
    p.anchor === anchor && return true
    try
        anchor === p.name && throw(ArgumentError("a part cannot be placed in its own frame"))
        kind, params, parts = ts.kind[], ts.params[], ts.parts[]
        # (a frame that leads back to this part is a loop: caught by anchor_frame)
        trial = copy(parts)
        trial[i] = _with(p; anchor)
        anchor_frame(anchor, kind, params, trial)
        vals = copy(p.values)
        if keep
            R0, t0 = anchor_frame(p.anchor, kind, params, parts)
            R1, t1 = anchor_frame(anchor, kind, params, parts)
            l = transpose(R1) * (R0 * part_center(p) + t0 - t1)
            vals[:x], vals[:y], vals[:z] = l
            a = transpose(R1) * (R0 * part_axis(p))             # the same chamber direction
            vals[:theta] = acosd(clamp(a[3], -1, 1))
            vals[:phi] = mod(atand(a[2], a[1]), 360)
        end
        ps = copy(ts.parts[])
        ps[i] = _with(p; anchor, values = vals)
        return _try_build!(ts; parts = ps)
    catch err
        err isa ArgumentError || rethrow()
        ts.error[] = _errmsg(err)
        return false
    end
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
    return "$kind · $(nelements(t)) elements"
end

"Names of a target's components, in order."
component_names(t::Target) = [c.name for c in t.components]

"Names of a target's opaque (non-aperture) components."
surface_names(t::Target) = [c.name for c in t.components if c.kind === :surface]
