# VisRad workspace (`.vrw`) reader: target components (primitives and CAD meshes) and
# laser beams with their pointing and spots.
#
# The workspace layout is undocumented (VR/files/workspace_files.html says only that it
# is a text file). Everything here was read off two workspaces (Workspace Format ID 11,
# target-component format 29): a NIF hohlraum with bundled CAD (VisRad 20.0.0, laser-beam
# format 22) and a primitive hohlraum + foil (VisRad 21.0.0, laser-beam format 23). Details those files
# do not fix are labeled *assumption*; see `reference/research-notes.md`, "VisRad
# workspace files".

"""
    VisRadObject

One target component of a VisRad workspace, as read by `read_visrad_workspace`.

Fields: `name`; VisRad's object `type` code and the `kind` Briareus recognizes
(`:mesh` for CAD/generalized meshes, `:sphere`, `:rectangle`, `:cylindrical_hohlraum`,
`:other`); the `folder` path in VisRad's component tree (`""` at the top); `params`,
every `key = value` line of the object verbatim; `rotation` and `position` (mm), the
object's frame in the target-chamber frame (`x_chamber = rotation * x_object +
position`); `mesh`, the surface in the chamber frame (mm) for meshes (one patch per
VisRad surface element) and `nothing` otherwise; `albedo`, `conversion` (x-ray
conversion efficiency) and `reflectivity` (laser); the flags `radiosity`, `laser`
(direct laser deposition) and `visible`; and `color` (RGB in [0, 1]).
"""
struct VisRadObject
    name::String
    type::Int
    kind::Symbol
    folder::String
    params::Dict{String,String}
    rotation::SMatrix{3,3,Float64,9}
    position::Vec3
    mesh::Union{Nothing,Mesh}
    albedo::Float64
    conversion::Float64
    reflectivity::Float64
    radiosity::Bool
    laser::Bool
    visible::Bool
    color::NTuple{3,Float64}
end

function Base.show(io::IO, o::VisRadObject)
    print(io, "VisRadObject(\"", o.name, "\", :", o.kind)
    o.mesh === nothing || print(io, ", ", ntriangles(o.mesh), " triangles")
    print(io, ")")
end

"""
    VisRadWorkspace

Contents of a VisRad workspace file (`read_visrad_workspace`): `path`, VisRad `version`
string, `objects` (`VisRadObject`s in file order), `beams` (`Beam`s), `clusters` (beam
cluster names in file order), and `rows`, the raw beam records (one `NamedTuple` per
beam with every column of the beam table, the reference frame name and any power
tables).
"""
struct VisRadWorkspace
    path::String
    version::String
    objects::Vector{VisRadObject}
    beams::Vector{Beam}
    clusters::Vector{String}
    rows::Vector{NamedTuple}
end

function Base.show(io::IO, ws::VisRadWorkspace)
    print(io, "VisRadWorkspace(\"", basename(ws.path), "\", VisRad ", ws.version, ", ",
          length(ws.objects), " objects, ", length(ws.beams), " beams, ",
          count(b -> _vr_v(b.power) > 0, ws.beams), " on)")
end

# VisRad object type codes seen in workspaces. Meanings of 1, 4 and 5 are inferred from
# their parameters and element counts (*assumption* for the code itself); 13 holds an
# explicit node/connectivity mesh. 11 is the small invisible "reference cone" VisRad
# creates for a multi-shape STEP import, which the shapes use as their reference frame
# (VR/grid/import_export_cad_files.html); it is not meshed.
const _VRW_KINDS = Dict(1 => :sphere, 4 => :rectangle, 5 => :cylindrical_hohlraum, 11 => :cone,
                        13 => :mesh)

# ---------------------------------------------------------------------------------------
# Text scanning
# ---------------------------------------------------------------------------------------

function _vrw_kv(line::AbstractString)
    i = findfirst('=', line)
    i === nothing && return nothing
    return strip(line[1:prevind(line, i)]), strip(line[nextind(line, i):end])
end

_vrw_float(d, k, default = 0.0) = (v = get(d, k, ""); x = tryparse(Float64, v); x === nothing ? Float64(default) : x)
_vrw_int(d, k, default = 0) = round(Int, _vrw_float(d, k, default))

# "... # table times = N ..." -> N
function _vrw_table_count(line)
    m = match(r"#\s*table\s+times\s*=\s*(\d+)", line)
    return m === nothing ? 0 : parse(Int, m[1])
end

function _vrw_find(lines, pat, from = 1, to = length(lines))
    for i in from:to
        startswith(lines[i], pat) && return i
    end
    return nothing
end

# ---------------------------------------------------------------------------------------
# Target components
# ---------------------------------------------------------------------------------------

# Rotation of VisRad's orientation angles (VR/grid/specifying_orient_angles.html):
# rotate about z by the Rotation angle, then about the original y by the Polar angle,
# then about the original z by the Azimuthal angle. The workspace stores "Tilt angle",
# "Rotation angle" and "Spin angle"; reading them as Polar, Azimuthal and Rotation, and
# the right-hand sense, are an *assumption* (the workspaces seen have only a 90° tilt).
_vrw_Rz(a) = SMatrix{3,3,Float64}(cosd(a), sind(a), 0, -sind(a), cosd(a), 0, 0, 0, 1)
_vrw_Ry(a) = SMatrix{3,3,Float64}(cosd(a), 0, -sind(a), 0, 1, 0, sind(a), 0, cosd(a))
_vrw_rotation(tilt, rot, spin) = _vrw_Rz(rot) * _vrw_Ry(tilt) * _vrw_Rz(spin)

# Position (cm, degrees) in Cartesian (0), cylindrical (1: r, z, φ, as in the summary
# table; *assumption*) or spherical (2: r, θ, φ) coordinates, to mm.
function _vrw_point(geom::Integer, p1, p2, p3)
    v = geom == 0 ? Vec3(p1, p2, p3) :
        geom == 1 ? Vec3(p1 * cosd(p3), p1 * sind(p3), p2) :
        geom == 2 ? p1 * port_vector(p2, p3) :
        throw(ArgumentError("unknown VisRad coordinate geometry $geom"))
    return 10 * v
end

# Raw object block -> (params, nodes, faces)
function _vrw_object_block(lines, i, stop)
    params = Dict{String,String}()
    nodes = Vec3[]
    faces = Vector{Vector{Int}}()
    nv = ns = 0
    while i <= stop
        l = lines[i]
        if startswith(l, "End Parameters for Object")
            return params, nodes, faces, i
        elseif startswith(l, "Time-Dependent")
            n = _vrw_table_count(l)
            n > 0 && (params[strip(split(l, " for ")[1])] = join(lines[i+1:i+n], '\n'))
            i += n
        elseif startswith(l, "Node positions:")
            sizehint!(nodes, nv)
            for k in 1:nv
                f = split(lines[i+k])
                push!(nodes, Vec3(parse(Float64, f[2]), parse(Float64, f[3]), parse(Float64, f[4])))
            end
            i += nv
        elseif startswith(l, "Surface element connectivity:")
            sizehint!(faces, ns)
            for k in 1:ns
                f = split(lines[i+k])
                vs = Int[]
                for s in @view f[2:end]
                    v = parse(Int, s)
                    v < 0 && break
                    push!(vs, v + 1)                 # 0-based in the file
                end
                push!(faces, vs)
            end
            i += ns
        else
            kv = _vrw_kv(l)
            if kv !== nothing
                k, v = kv
                params[k] = v
                k == "Number of vertices" && (nv = parse(Int, v))
                k == "Number of surfaces" && (ns = parse(Int, v))
            end
        end
        i += 1
    end
    throw(ArgumentError("unterminated object block (truncated file?)"))
end

# Directory tree: "Number of target components = n", n names, "Number of child folders
# = m", then m × (folder name, subtree).
function _vrw_folders!(out, lines, i, path)
    n = parse(Int, last(_vrw_kv(lines[i])))
    for k in 1:n
        out[strip(lines[i+k])] = path
    end
    i += n + 1
    m = parse(Int, last(_vrw_kv(lines[i])))
    i += 1
    for _ in 1:m
        name = strip(lines[i])
        i = _vrw_folders!(out, lines, i + 1, isempty(path) ? name : path * "/" * name)
    end
    return i
end

function _vrw_flip(m::Mesh)
    Mesh(m.vertices, [(a, c, b) for (a, b, c) in m.triangles], m.patches)
end

# Primitive surfaces in the object's own frame (mm), as Briareus components.
# Grid conventions inferred from VisRad's own element counts (power-source table of the
# primitive workspace): a rectangle with 10 × 10 mesh points has 81 elements, a
# cylindrical hohlraum with mesh points (11, 12, 4) has 216 — 10 × 12 wall, 2 × 3 × 12
# end-cap and 2 × 12 laser-entrance-hole elements. Periodic azimuths count elements,
# other directions count nodes. With `h` (mm) the grid is Briareus's for elements of
# about that size instead (`hohlraum_counts`, `capsule_counts`; VisRad's counts ignored).
function _vrw_primitive(kind, name, p, mat; h = nothing)
    L(k) = 10 * _vrw_float(p, "Length $k value")
    N(k) = _vrw_int(p, "Number of mesh pts $k")
    h === nothing || h > 0 || throw(ArgumentError("VisRad object $name: the element size must be positive"))
    cells(len) = max(1, ceil(Int, len / h))
    flip = _vrw_int(p, "Normals in/out") != 0
    _vrw_int(p, "Use non-zero wall thickness") != 0 &&
        @warn "VisRad object $name has a wall thickness; Briareus models its outer surface only"
    if kind === :rectangle
        # Length 1 × Length 2 (x × y) centered in the object's x-y plane, normal +z
        m = h === nothing ? rectangle_mesh(L(1), L(2); nu = max(1, N(1) - 1), nv = max(1, N(2) - 1)) :
            rectangle_mesh(L(1), L(2); nu = cells(L(1)), nv = cells(L(2)))
        return [Component(name, flip ? _vrw_flip(m) : m, mat)]
    elseif kind === :sphere
        # Length 1 = radius; Lengths 2–5 = angular limits (full sphere only)
        Ls = Tuple(_vrw_float(p, "Length $k value") for k in 2:5)
        Ls == (0.0, 360.0, 0.0, 180.0) ||
            throw(ArgumentError("VisRad sphere $name: partial spheres are not supported"))
        nt, np = h === nothing ? (max(2, N(1) - 1), max(3, N(2))) : values(capsule_counts(L(1), h))
        m = sphere_mesh(L(1); ntheta = nt, nphi = np, outward = !flip)
        return [Component(name, m, mat; surface = SphereSurface(L(1)))]
    elseif kind === :cylindrical_hohlraum
        # Length 1 = radius, Length 2 = length, Length 3 = LEH radius; mesh points 1 =
        # axial nodes, 2 = azimuthal elements, 3 = radial nodes on each end cap. VisRad
        # hohlraum normals point in by default. Each LEH is one ring of elements.
        R, H, r = L(1), L(2), L(3)
        # VisRad's grid: one ring of LEH elements; Briareus's: `hohlraum`'s (LEH rings of
        # about h, half the wall's azimuthal count)
        nz, nphi, nr, nr_leh, nphi_leh = if h === nothing
            max(1, N(1) - 1), max(3, N(2)), max(1, N(3) - 1), 1, max(3, N(2))
        else
            c = hohlraum_counts(; length = H, radius = R, leh_radius = min(r, R), h)
            c.nz, c.nphi, c.nr_cap, c.nr_leh, c.nphi ÷ 2
        end
        z = H / 2
        cyl = CylinderSurface(R)
        nm(s) = Symbol(name, "_", s)
        comps = Component[Component(nm(:wall), cylinder_mesh(R, -z, z; nz, nphi, inward = !flip), mat;
                                    surface = cyl)]
        if r < R
            push!(comps, Component(nm(:cap_top), annulus_mesh(r, R, z; nr, nphi, up = flip), mat;
                                   surface = DiskSurface(center = (0, 0, z))))
            push!(comps, Component(nm(:cap_bottom), annulus_mesh(r, R, -z; nr, nphi, up = !flip), mat;
                                   surface = DiskSurface(center = (0, 0, -z))))
        end
        if r > 0
            push!(comps, aperture(nm(:leh_top), disk_mesh(min(r, R), z; nr = nr_leh, nphi = nphi_leh, up = false);
                                  surface = DiskSurface(center = (0, 0, z))))
            push!(comps, aperture(nm(:leh_bottom), disk_mesh(min(r, R), -z; nr = nr_leh, nphi = nphi_leh, up = true);
                                  surface = DiskSurface(center = (0, 0, -z))))
        end
        return comps
    end
    return Component[]
end

function _vrw_objects(lines, a, b)
    blocks = []
    i = a
    while (j = _vrw_find(lines, "Begin Parameters for Object", i, b)) !== nothing
        name = String(strip(last(_vrw_kv(lines[j]))))
        params, nodes, faces, i = _vrw_object_block(lines, j + 1, b)
        push!(blocks, (; name, params, nodes, faces))
    end
    folders = Dict{String,String}()
    d = findnext(l -> occursin("Begin Target Component Directory Structure", l), lines, a)
    if d !== nothing
        k = findnext(l -> occursin("Number of target components", l), lines, d)
        k === nothing || _vrw_folders!(folders, lines, k, "")
    end
    byname = Dict(blk.name => blk for blk in blocks)
    frames = Dict{String,Tuple{SMatrix{3,3,Float64,9},Vec3}}()
    function frame(name, depth = 0)
        haskey(frames, name) && return frames[name]
        depth > length(blocks) && throw(ArgumentError("cyclic VisRad reference frames at $name"))
        p = byname[name].params
        _vrw_int(p, "Set position of ID") == 0 ||
            throw(ArgumentError("VisRad object $name is positioned by a key point other than its centroid (not supported)"))
        ref = get(p, "Reference coord. system name", "Target Chamber")
        R0, t0 = _vrw_frame_of(ref, n -> haskey(byname, n) ? frame(n, depth + 1) : nothing)
        x = _vrw_point(_vrw_int(p, "Coord. sys. geometry type"), _vrw_float(p, "Position 1"),
                       _vrw_float(p, "Position 2"), _vrw_float(p, "Position 3"))
        R = _vrw_rotation(_vrw_float(p, "Tilt angle"), _vrw_float(p, "Rotation angle"),
                          _vrw_float(p, "Spin angle"))
        return frames[name] = (R0 * R, R0 * x + t0)
    end
    objs = VisRadObject[]
    for blk in blocks
        p = blk.params
        type = _vrw_int(p, "Object type", -1)
        kind = get(_VRW_KINDS, type, :other)
        R, t = frame(blk.name)
        mesh = nothing
        if kind === :mesh
            # nodes are in the object's own frame (cm), centered on its centroid; the
            # "Gen. mesh object origin" (= -Position) maps them back to the CAD file's
            # frame and is not needed here
            vs = [R * (10 * v) + t for v in blk.nodes]
            tris = NTuple{3,Int}[]
            patches = UnitRange{Int}[]
            for f in blk.faces
                length(f) >= 3 || continue
                t0 = length(tris) + 1
                for k in 2:length(f)-1
                    push!(tris, (f[1], f[k], f[k+1]))
                end
                push!(patches, t0:length(tris))
            end
            mesh = Mesh(vs, tris, patches)
            _vrw_int(p, "Normals in/out") != 0 && (mesh = _vrw_flip(mesh))
        end
        for (key, what) in (("Albedo model ID", "albedo"), ("Xray conv. effic. model ID", "x-ray conversion"),
                            ("Laser reflectivity model ID", "laser reflectivity"))
            _vrw_int(p, key) != 0 && @warn "VisRad object $(blk.name): time-dependent $what is not imported; using its constant value"
        end
        push!(objs, VisRadObject(blk.name, type, kind, get(folders, blk.name, ""), p, R, t, mesh,
                                 _vrw_float(p, "Albedo value"), _vrw_float(p, "Xray conv. effic. value"),
                                 _vrw_float(p, "Laser reflectivity value"),
                                 _vrw_int(p, "Include in radiosity calc.") != 0,
                                 _vrw_int(p, "Include in dir laser dep calc.") != 0,
                                 _vrw_int(p, "Object visibility") != 0,
                                 (_vrw_float(p, "Red value"), _vrw_float(p, "Green value"),
                                  _vrw_float(p, "Blue value"))))
    end
    return objs
end

function _vrw_frame_of(ref::AbstractString, lookup)
    r = strip(ref)
    (isempty(r) || lowercase(r) in ("target chamber", "tcc")) && return (SMatrix{3,3,Float64}(I), zero(Vec3))
    f = lookup(r)
    f === nothing && throw(ArgumentError("unknown VisRad reference coordinate system \"$r\" " *
                                         "(target mounts and diagnostic ports are not supported)"))
    return f
end

# ---------------------------------------------------------------------------------------
# Laser beams
# ---------------------------------------------------------------------------------------

# Columns of the per-beam row (laser-beam data formats 22 and 23), named from the
# file's own header. "rot_id" is blank in both files, leaving 52 tokens of 53.
const _VRW_BEAM_COLUMNS = (:port, :theta, :phi, :port_radius, :custom, :cone_min, :cone_max,
    :fnumber, :red, :green, :blue, :transparency, :method, :port_distance, :geometry,
    :p1, :p2, :p3, :focus, :rg_model, :on_blue, :on_green, :on_red, :vis_blue, :vis_green,
    :vis_red, :mode_blue, :mode_green, :mode_red, :power_blue, :power_green, :power_red,
    :random_mult, :energy_norm, :energy, :use_mult_blue, :use_mult_green, :use_mult_red,
    :mult_blue, :mult_green, :mult_red, :profile, :size_model, :major, :elliptical, :minor,
    :rotation, :rot_model, :rot_id, :order, :envelope_ratio, :lens_ratio, :dpr)

function _vrw_beam_row(line)
    f = split(line)
    n = length(_VRW_BEAM_COLUMNS)
    length(f) == n - 1 && insert!(f, findfirst(==(:rot_id), _VRW_BEAM_COLUMNS), "")
    length(f) == n || throw(ArgumentError("VisRad beam row has $(length(f)) fields, expected $n: $line"))
    vals = [k === :rot_id ? String(s) : parse(Float64, s) for (k, s) in zip(_VRW_BEAM_COLUMNS, f)]
    return NamedTuple{_VRW_BEAM_COLUMNS}(Tuple(vals))
end

# Built-in phase plate with these radii (µm) and exponent, if any: the workspace stores
# a profile-model number with the plate's tabulated numbers copied into the row, so
# the numbers identify the plate (and with it the intensity level the radii refer to).
function _vrw_plate_name(major, minor, order)
    for (k, v) in VISRAD_PHASE_PLATES
        isapprox(v.major, major; rtol = 1e-6) && isapprox(v.minor, minor; rtol = 1e-6) &&
            isapprox(v.order, order; rtol = 1e-6) && return k
    end
    return ""
end

function _vrw_beams(lines, a, b)
    clusters = String[]
    c = _vrw_find(lines, "Cluster Names:", a, b)
    if c !== nothing
        k = c + 1
        while k <= b && !isempty(strip(lines[k]))
            push!(clusters, strip(lines[k]))
            k += 1
        end
    end
    useran = false
    for l in @view lines[a:something(c, b)]
        kv = _vrw_kv(l)
        kv !== nothing && kv[1] == "use random beam pwr mult" && (useran = kv[2] != "0")
    end
    rows = NamedTuple[]
    i = a
    while (j = _vrw_find(lines, "Laser Beam Name", i, b)) !== nothing
        name = String(strip(last(_vrw_kv(lines[j]))))
        cluster = String(strip(last(_vrw_kv(lines[j+1]))))
        row = _vrw_beam_row(lines[j+2])
        k = j + 3
        info = Dict{String,String}()
        tables = Dict{String,Vector{String}}()
        while k <= b && !startswith(lines[k], "Laser Beam Name") && !startswith(lines[k], "[End")
            l = lines[k]
            m = match(r"^(Ref\. Coord\. Sys\.|D toward Port ID|Pulse Shape Name):\s*(.*)$", l)
            m === nothing || (info[m[1]] = strip(m[2]))
            if occursin("table values:", l)
                n = _vrw_table_count(l)
                tables[strip(split(l, ":")[1])] = lines[k+1:k+n]
                k += n
            end
            k += 1
        end
        push!(rows, merge((; name, cluster, ref = get(info, "Ref. Coord. Sys.", "Target Chamber"),
                           port_id = get(info, "D toward Port ID", ""),
                           pulse = get(info, "Pulse Shape Name", ""), tables, random_used = useran), row))
        i = k
    end
    return clusters, rows
end

# Fraction of the power column that is in the blue (3ω) cone. VisRad's red/green power
# models (VR/power_sources/set_laser_beam_energetics.html), numbered in the dialog's
# order: 0 None, 1 Absolute Power (TW), 2 Fraction of Total Beam Power, 3 Relative to
# Blue Cone Power. Only under model 2 is the power column the total, with the blue cone
# carrying 1 - f_G - f_R of it; the fractions count whether or not those cones are on.
# The numbering and the on/off rule are read off the foil workspace, where four 1 TW
# beams with f_G = f_R = 0.05 (green off) deposit 3.60 TW.
_vrw_blue_fraction(r) = round(Int, r.rg_model) == 2 ? 1 - r.power_green - r.power_red : 1.0

# Beam record (as `_vr_apply` takes it) from a raw row.
function _vrw_record(r, frame_of)
    method = round(Int, r.method)
    aim = if method == 0
        zero(Vec3)
    elseif method == 2
        R, t = _vrw_frame_of(r.ref, frame_of)
        R * _vrw_point(round(Int, r.geometry), r.p1, r.p2, r.p3) + t
    else
        throw(ArgumentError("beam $(r.name): pointing toward a diagnostic port (method $method) is not supported"))
    end
    on = r.on_blue != 0
    (r.on_green != 0 || r.on_red != 0) &&
        @warn "beam $(r.name): unconverted (2ω/1ω) light is on in VisRad; only the 3ω power is imported (see `nif_unconverted`)"
    r.mode_blue != 0 && @warn "beam $(r.name): time-dependent power is not imported; using the constant power column"
    r.energy_norm != 0 && @warn "beam $(r.name): power is normalized to an energy in VisRad; using the power column"
    r.dpr != 0 && @warn "DPR is on for beam $(r.name); Briareus does not model DPR"
    P = r.power_blue * _vrw_blue_fraction(r) * (r.use_mult_blue != 0 ? r.mult_blue : 1.0) *
        (r.random_used ? r.random_mult : 1.0)
    minor = r.elliptical != 0 ? r.minor : r.major
    plate = _vrw_plate_name(r.major, minor, r.order)
    spot = (model = plate, major = r.major, minor = minor, rot = r.rotation, order = r.order, level = nothing)
    return (name = r.name, on = on, theta = r.theta, phi = r.phi, aim = aim, defocus = 10 * r.focus,
            power = on ? P : 0.0, spot = spot)
end

# ---------------------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------------------

"""
    read_visrad_workspace(path; template = nothing, missing = :off, spot_definition = :auto,
                          rotation_sense = :beam, angle_tol = 0.05) -> VisRadWorkspace

Read a VisRad workspace (`.vrw`, the text format of VisRad 20–21, "Workspace Format
ID = 11"): its target components, including CAD imported as generalized meshes, and its
laser beams. Lengths are converted from cm to mm; every object and beam frame is
resolved to the target-chamber frame. Build a Briareus target from the components with
`visrad_target`.

**Components** become `VisRadObject`s. CAD meshes keep VisRad's surface elements as
patches (one per triangle). Spheres, rectangles and cylindrical hohlraums are meshed on
VisRad's grid when converted by `visrad_target`. Objects are placed by VisRad's
orientation angles and position (VR/grid/specifying_orient_angles.html), relative to the
chamber or to another component.

**Beams.** Port angles and f-number come from the file. Pointing at TCC or at a point
given in the chamber's or a component's frame is supported; the point is the aim and
the focus offset (cm) the defocus (mm), as in `read_visrad_summary`. Power is the 3ω
(blue cone) power in TW for beams with 3ω on, and 0 otherwise: the power column times
the beam's power multiplier when enabled and, under VisRad's "Fraction of Total Beam
Power" red/green model (where that column is the total), times 1 − f_G − f_R. Spots are VisRad's super-Gaussian radii (µm), exponent and minor-axis
rotation, converted as in `read_visrad_summary`: radii that match a built-in phase
plate (`VISRAD_PHASE_PLATES`) are read at that plate's level (HWHM for NIF CPPs), others
at `spot_definition` (`:auto` = 1/e). Unconverted light, time-dependent powers, energy
normalization and DPR are warned about and not imported.

With a `template` (e.g. `nif_beams()`), the file's pointing, power and spots are
applied to the template's beams by name and its other properties (halo, wavelength,
quad, group) are kept; port angles are checked to `angle_tol` degrees and `missing`
handles template beams absent from the file, as in `read_visrad_summary`. Without one,
beams are built from the file, grouped by VisRad cluster.
"""
function read_visrad_workspace(path::AbstractString; template = nothing, missing::Symbol = :off,
                               spot_definition = :auto, rotation_sense::Symbol = :beam,
                               angle_tol = 0.05)
    _vr_sense(rotation_sense)
    lines = _vr_lines(path)
    occursin("Workspace Format ID", first(lines)) ||
        throw(ArgumentError("$path is not a VisRad text workspace (first line: \"$(first(lines))\")"))
    vi = findfirst(l -> occursin("VisRad Version", l), lines)
    version = vi === nothing ? "" : String(strip(replace(lines[vi], r".*VisRad Version_?" => "")))
    section(name) = (a = _vrw_find(lines, "[" * name * "]"); a === nothing ? nothing :
                     (a, something(_vrw_find(lines, "[End " * name * "]", a), length(lines))))
    ts = section("Target Component Data")
    objects = ts === nothing ? VisRadObject[] : _vrw_objects(lines, ts...)
    frames = Dict(o.name => (o.rotation, o.position) for o in objects)
    frame_of = n -> get(frames, n, nothing)
    ls = section("Laser Beam Data")
    clusters, rows = ls === nothing ? (String[], NamedTuple[]) : _vrw_beams(lines, ls...)
    recs = Dict{String,Any}()
    order = String[]
    for r in rows
        k = _vr_key(r.name)
        recs[k] = _vrw_record(r, frame_of)
        push!(order, k)
    end
    beams = if template === nothing
        out = Beam[]
        for (k, r) in zip(order, rows)
            rec = recs[k]
            sp = _vr_spot(rec.spot, nothing, spot_definition, rotation_sense)
            push!(out, Beam(r.name; theta = r.theta, phi = r.phi, power = rec.power, aim = rec.aim,
                            defocus = rec.defocus, spot = sp, fnumber = r.fnumber,
                            group = Symbol(r.cluster)))
        end
        identity.(out)
    else
        _vr_apply(template, recs, missing; angle_tol, spot_definition, sense = rotation_sense)
    end
    return VisRadWorkspace(String(path), version, objects, beams, clusters, rows)
end

"""
    visrad_profile(workspace; growth = 0.67) -> VisRadProfile

VisRad's beam model with the workspace's parameters: each beam's envelope ratio (by
name; VisRad's default ratio differs by phase plate, 1.27–1.54 for the NIF plates) and
the final-optic distance (port radius, cm → mm) of the beams that are on (of all beams
if none is on; beams that differ in it are warned about, and the largest is used).
`growth` is `VisRadProfile`'s. Pass it to `DepositionOperator(target, beams; profile)`,
or to a `Scenario` through `deposition_kw = (; profile)`.
"""
function visrad_profile(ws::VisRadWorkspace; growth = 0.67)
    rows = filter(r -> r.on_blue != 0, ws.rows)
    isempty(rows) && (rows = ws.rows)
    isempty(rows) && return VisRadProfile(; growth)
    fl = unique(10 * r.port_radius for r in rows)
    length(fl) > 1 && @warn "beams differ in final-optic distance; using the largest" distances_mm = fl
    return VisRadProfile(; envelope = Dict(r.name => r.envelope_ratio for r in ws.rows),
                         focal_length = maximum(fl), growth)
end

"""
    visrad_component(object; material = nothing, two_sided = false, orient = false,
                     toward = (0, 0, 0), coarsen = nothing, remesh = nothing,
                     fit = false, primitive_h = nothing) -> Vector{Component}

Briareus components for one `VisRadObject`, in the chamber frame (mm). The material
defaults to the object's albedo and x-ray conversion efficiency. CAD meshes give one
component radiating from the front of VisRad's surface elements (counterclockwise
winding, *assumption*), or from both sides when `two_sided`. VisRad does not let users
edit the normals of imported CAD (VR/grid/import_export_cad_files.html), and tessellated
CAD often faces the wrong way: in the NIF workspace this reader was built from, the
inner hohlraum wall faces away from the cavity. `orient = true` repairs the winding
with `orient_mesh`: consistent across shared edges, closed solids facing out of their
material, open sheets facing the point `toward` (mm, chamber frame; by default the
target chamber center). `two_sided = true` instead radiates from both sides (thin
shells are then modeled as two radiating sides). CAD needs fewer, compact elements
before view factors: `remesh` (mm) divides the surface into near-planar elements of
about that size whatever the facets' shapes (`remesh_mesh`; long slivers are cut),
while `coarsen` (mm) only merges whole near-coplanar facets up to that size
(`coarsen_mesh`; slivers longer than the size stay). `fit = true` (or a named tuple of
`curved_components` and `fit_surfaces` keywords) then splits a CAD part into its
cylinders, cones and spheres, each with its exact surface fitted to the CAD facets, and
the rest (`curved_components`), so the radiation model integrates on the curved walls
as it does for a primitive hohlraum. Rectangles give one component (normal +z before orientation), spheres one with
a `SphereSurface`, and cylindrical hohlraums five — `<name>_wall`, `_cap_top`,
`_cap_bottom` and the apertures `_leh_top`, `_leh_bottom` — with the curved surfaces
`hohlraum` uses, all on VisRad's grid (`Number of mesh pts`, which reproduces VisRad's
element table), or with `primitive_h` (mm) on Briareus's grid for elements of about
that size (`hohlraum_counts`, `capsule_counts`), e.g. to resolve a spot on a foil that
VisRad meshes coarsely. Other object types return no components (with a warning).
"""
function visrad_component(o::VisRadObject; material = nothing, two_sided::Bool = false,
                          orient::Bool = false, toward = (0.0, 0.0, 0.0), coarsen = nothing,
                          remesh = nothing, fit = false, primitive_h = nothing)
    coarsen === nothing || remesh === nothing ||
        throw(ArgumentError("give one of `coarsen` and `remesh`, not both"))
    mat = material === nothing ? Material(Symbol(o.name); albedo = o.albedo, conversion = o.conversion) :
          material
    o.reflectivity > 0 && @warn "VisRad object $(o.name) reflects laser light ($(o.reflectivity)); Briareus absorbs it all"
    if o.kind === :mesh
        facets = orient ? orient_mesh(o.mesh; toward) : o.mesh
        mesh = coarsen === nothing ? facets : coarsen_mesh(facets; size = coarsen)
        remesh === nothing || (mesh = remesh_mesh(mesh; size = remesh))
        c = Component(o.name, mesh, mat; two_sided)
        fit === false && return [c]
        # fitted on the CAD facets, whose vertices lie on the surfaces
        return curved_components(c; reference = mesh === facets ? nothing : facets,
                                 (fit === true ? (;) : fit)...)
    end
    cs = _vrw_primitive(o.kind, o.name, o.params, mat; h = primitive_h)
    cs = [c.kind === :surface && two_sided ? Component(c.name, c.mesh, c.material; two_sided, surface = c.surface) : c
          for c in cs]
    isempty(cs) && (@warn "VisRad object $(o.name) (type $(o.type)) is not supported; skipped"; return cs)
    return [transform(c, o.rotation, o.position) for c in cs]
end

"""
    visrad_target(workspace; select = o -> o.radiosity || o.laser, materials = Dict(),
                  two_sided = false, orient = false, toward = (0, 0, 0),
                  coarsen = nothing, remesh = nothing, fit = false,
                  inside = nothing) -> Target

Target built from the components of a VisRad workspace (`read_visrad_workspace`) for
which `select(object)` is true — by default those VisRad includes in its radiosity or
direct-laser calculation. `materials` maps object names to a `Material` overriding the
file's albedo and conversion efficiency. `two_sided` (a `Bool`, or a function of the
object) makes surfaces radiate from both sides; see `visrad_component` for why CAD
meshes often need it, and for `orient`, `toward`, `remesh`, `coarsen` and `fit` (which
apply to CAD meshes: repaired windings, elements of about `remesh` or `coarsen` mm, and
exact cylinders, cones and spheres for the radiation model).

CAD meshes can hold hundreds of thousands of triangles, and each becomes a radiating
element. That is fine for ray tracing (deposition, clearance, imaging) but too many for
view factors; select the parts that matter for the radiation problem, e.g.
`select = o -> o.folder == "Hohlraum"` or by name, and remesh them. `inside`, a point
(mm, chamber frame) in the cavity, keeps only the elements that bound the free space
around it (`enclosure` with its defaults; call it on the full target for the options
and a per-component account of what was dropped). It needs consistent windings
(`orient = true`). The selection runs on the CAD facets as oriented, `remesh` or
`coarsen` then divide only the facets kept (much cheaper when most of the assembly lies
outside the cavity), and a second selection on the resulting elements drops those
that average too little flux. `fit` then applies to the elements kept, fitted on the
kept facets. The CAD meshes go through these steps in `cad_target`, which takes the
other components as they are (into the selection, but not remeshed or fitted); it
does the same for CAD read from STL files.
"""
function visrad_target(ws::VisRadWorkspace; select = o -> o.radiosity || o.laser,
                       materials = Dict{String,Material}(), two_sided = false,
                       orient::Bool = false, toward = (0.0, 0.0, 0.0), coarsen = nothing,
                       remesh = nothing, fit = false, inside = nothing)
    coarsen === nothing || remesh === nothing ||
        throw(ArgumentError("give one of `coarsen` and `remesh`, not both"))
    parts = Union{Component,Pair{Symbol,Mesh}}[]
    mats = Dict{Symbol,Material}()        # of the CAD parts, for cad_target
    sides = Dict{Symbol,Bool}()
    for o in ws.objects
        select(o) || continue
        (o.radiosity && o.laser) || @warn "VisRad object $(o.name) is included in only one of VisRad's radiosity and laser calculations; Briareus uses it in both"
        ts = two_sided isa Bool ? two_sided : Bool(two_sided(o))
        material = get(materials, o.name, nothing)
        if o.kind === :mesh
            # the CAD part as VisRad has it (material, sidedness); cad_target does the rest
            c = only(visrad_component(o; material, two_sided = ts))
            push!(parts, c.name => c.mesh)
            mats[c.name], sides[c.name] = c.material, c.two_sided
        else
            append!(parts, visrad_component(o; material, two_sided = ts, orient, toward, coarsen, remesh, fit))
        end
    end
    isempty(parts) && throw(ArgumentError("no VisRad components selected"))
    return cad_target(parts; orient, toward, inside, remesh, coarsen, fit, material = mats, two_sided = sides)
end

"""
    write_visrad_workspace(path, template, beams; angle_tol = 0.05) -> path

A VisRad workspace with the pointing, focus, power and on/off of `beams` written into a
copy of the workspace `template` (a `.vrw` file, e.g. the one a design was imported
from with `read_visrad_workspace`). Each beam row of the template whose name matches a
beam gets that beam's aim, as a Cartesian position in the row's own reference frame
(cm), its focus offset `defocus` (cm), its power as the 3ω power column (divided by
the row's power multipliers and, under the "Fraction of Total Beam Power" red/green
model, by its blue fraction, so VisRad's blue cone carries the beam's power) and on (power > 0)
or off. Everything else is copied byte for byte: the target components and their
frames, the spots, pulse shapes and VisRad's settings stay VisRad's own. Numbers are
written in place at the template's column widths and decimals (0.1 µm for positions).

Rows without a matching beam are left as they were, and beams without a row are
reported in a warning: the template sets which beams exist. A beam whose port differs
from its row's by more than `angle_tol` degrees is refused (ports cannot be changed).
Time-dependent power tables and energy normalization in the template are left as they
are, so a row that uses them does not take the constant power (warned).
"""
function write_visrad_workspace(path::AbstractString, template::AbstractString,
                                beams::AbstractVector{<:Beam}; angle_tol = 0.05)
    ws = read_visrad_workspace(template)
    frames = Dict(o.name => (o.rotation, o.position) for o in ws.objects)
    byname = Dict(b.name => b for b in beams)
    rows = Dict(r.name => r for r in ws.rows)
    text = read(template, String)
    lines = split(text, '\n')                       # line endings (\r) stay on the lines
    used = Set{String}()
    for (j, l) in enumerate(lines)
        startswith(l, "Laser Beam Name") || continue
        name = String(strip(last(_vrw_kv(rstrip(l, '\r')))))
        haskey(byname, name) && haskey(rows, name) || continue
        b, r = byname[name], rows[name]
        θ, φ = port_angles(b)
        (abs(θ - r.theta) <= angle_tol && abs(mod(φ - r.phi + 180, 360) - 180) <= angle_tol) ||
            throw(ArgumentError("beam $name: port (θ $θ, φ $φ) differs from the workspace's ($(r.theta), $(r.phi))"))
        (r.mode_blue != 0 || r.energy_norm != 0) &&
            @warn "beam $name: VisRad uses a power table or energy normalization for it; its constant power is written but not used"
        R, t = _vrw_frame_of(r.ref, n -> get(frames, n, nothing))
        p = transpose(R) * (Vec3(_vr_v.(b.aim)) - t) / 10
        P = _vr_v(b.power)
        mult = _vrw_blue_fraction(r) * (r.use_mult_blue != 0 ? r.mult_blue : 1.0) *
               (r.random_used ? r.random_mult : 1.0)
        mult > 0 || throw(ArgumentError("beam $name: the workspace's power multiplier (or blue fraction) is zero"))
        new = Dict(:method => 2, :geometry => 0, :p1 => p[1], :p2 => p[2], :p3 => p[3],
                   :focus => _vr_v(b.defocus) / 10, :on_blue => P > 0 ? 1 : 0,
                   :power_blue => P > 0 ? P / mult : r.power_blue)
        lines[j+2] = _vrw_patch_row(lines[j+2], new)
        push!(used, name)
    end
    missing_rows = [b.name for b in beams if !(b.name in used)]
    isempty(missing_rows) ||
        @warn "$(length(missing_rows)) beams have no row in the template workspace and were not written" beams = missing_rows
    write(path, join(lines, '\n'))
    return path
end


const _VRW_DECIMALS = Dict(:p1 => 6, :p2 => 6, :p3 => 6, :focus => 6, :power_blue => 4)

# A beam row with the given columns rewritten in place: each value right-aligned where
# the old token ended, with the old token's decimals (integers for integer tokens), at
# least `_VRW_DECIMALS` of them (VisRad's own: 0.1 µm for positions in cm).
function _vrw_patch_row(line::AbstractString, new::AbstractDict)
    cr = endswith(line, '\r')
    s = cr ? chop(line) : String(line)
    spans = [m.offset:m.offset+ncodeunits(m.match)-1 for m in eachmatch(r"\S+", s)]
    n = length(_VRW_BEAM_COLUMNS)
    blank = length(spans) == n - 1            # the "rot_id" column is blank
    length(spans) == n || blank || throw(ArgumentError("unexpected VisRad beam row: $line"))
    ri = findfirst(==(:rot_id), _VRW_BEAM_COLUMNS)
    tok(k) = (i = findfirst(==(k), _VRW_BEAM_COLUMNS); blank && i > ri ? i - 1 : i)
    edits = sort([(tok(k), k, v) for (k, v) in new]; by = first, rev = true)
    for (i, k, v) in edits
        sp = spans[i]
        old = s[sp]
        d = occursin('.', old) ? max(length(old) - findfirst('.', old), get(_VRW_DECIMALS, k, 0)) : 0
        txt = d == 0 ? string(round(Int, v)) : @sprintf("%.*f", d, v)
        from = i == 1 ? 1 : last(spans[i-1]) + 1       # the field: its leading spaces and token
        w = last(sp) - from + 1
        s = s[1:from-1] * lpad(txt, max(w, length(txt) + 1)) * s[last(sp)+1:end]
    end
    return cr ? s * "\r" : s
end
