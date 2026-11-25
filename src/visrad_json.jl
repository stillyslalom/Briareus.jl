# This file reads the JSON workspace of VisRad (`.vrw`, "Workspace Format ID": 1001),
# which VisRad 21 writes when it saves a workspace that it has computed. The JSON format
# holds the content of the text format under readable keys. This file maps it onto the
# structures of the text reader, so that `read_visrad_workspace` and
# `visrad_power_sources` take either format. The layout was taken from VisRad 21.0.0
# files (research notes, "VisRad JSON workspaces").

function _vrj_parse(path::AbstractString)
    d = JSON.parsefile(path)
    get(d, "Workspace Format ID", nothing) isa Integer ||
        throw(ArgumentError("$path is JSON but not a VisRad workspace (no \"Workspace Format ID\")"))
    return d
end

# A JSON value as text: VisRad writes numbers in arrays as strings and scalars as numbers.
_vrj_str(x::AbstractString) = String(x)
_vrj_str(x::Integer) = string(x)
_vrj_str(x::Real) = isinteger(x) && abs(x) < 1e15 ? string(Int(x)) : string(x)
_vrj_str(x) = string(x)
_vrj_num(x::Real) = Float64(x)
_vrj_num(x::AbstractString) = parse(Float64, x)

# The elements "<prefix>[1]", "<prefix>[2]", … of a JSON object, in order.
function _vrj_elements(d, prefix)
    out = []
    k = 1
    while haskey(d, "$prefix[$k]")
        push!(out, d["$prefix[$k]"])
        k += 1
    end
    return out
end

function _vrj_section(d, key, path)
    haskey(d, key) || throw(ArgumentError("$path: no \"$key\" in the JSON workspace"))
    return d[key]
end

# Target components, as the object block of the text format that `_vrw_build_objects`
# takes. Keys are those of the text format, the parameters that `_vrw_build_objects` and
# `_vrw_primitive` use, except "Object color" (RGB) and "Gen. mesh object origin", which
# are arrays, and the node array and connectivity columns.
function _vrj_block(o, path)
    name = String(strip(_vrj_str(o["Object name"])))
    params = Dict{String,String}()
    for (k, v) in o
        v isa AbstractVector || v isa AbstractDict || (params[k] = _vrj_str(v))
    end
    if haskey(o, "Object color")
        rgb = _vrj_str.(o["Object color"])
        length(rgb) == 3 || throw(ArgumentError("$path: object $name: \"Object color\" is not RGB"))
        params["Red value"], params["Green value"], params["Blue value"] = rgb
    end
    nodes = Vec3[]
    faces = Vector{Vector{Int}}()
    if haskey(o, "Object node array-Num elements")
        nv = Int(o["Object node array-Num elements"])
        sizehint!(nodes, nv)
        for k in 1:nv
            v = o["Object node array[$k]"]
            push!(nodes, Vec3(_vrj_num(v[1]), _vrj_num(v[2]), _vrj_num(v[3])))
        end
        ns = Int(o["Object connectivity array-Num rows"])
        nc = Int(o["Object connectivity array-Num cols"])
        cols = [o["Object connectivity array-Values Col $c"] for c in 1:nc]
        all(c -> length(c) == ns, cols) ||
            throw(ArgumentError("$path: object $name: connectivity columns are not $ns long"))
        sizehint!(faces, ns)
        for r in 1:ns
            vs = Int[]
            for c in 1:nc
                v = parse(Int, _vrj_str(cols[c][r]))
                v < 0 && break
                push!(vs, v + 1)                         # 0-based in the file
            end
            push!(faces, vs)
        end
        params["Number of vertices"], params["Number of surfaces"] = string(nv), string(ns)
    end
    return (; name, params, nodes, faces)
end

# Folders from "Target comp list structure", with each component as "Target
# Components/<folder path>/<name>". The folder of the text reader is the path between the
# two, and "" at the top.
function _vrj_folders(d)
    out = Dict{String,String}()
    s = get(d, "Target comp list structure", nothing)
    s === nothing && return out
    for p in get(s, "Target comp name", [])
        parts = split(String(p), '/')
        length(parts) >= 2 || continue
        out[parts[end]] = join(parts[2:end-1], '/')
    end
    return out
end

# The beam row of the text format (`_VRW_BEAM_COLUMNS`) from a JSON beam, as
# `key => column`. "Beam cone F numuber" is the spelling of VisRad. The cone half-angles
# (`cone_min`, `cone_max`, not stored) are those of VisRad, 0 and atan(1/2F), in degrees.
const _VRJ_BEAM_KEYS = (
    :port => "Port index", :theta => "Port polar angle", :phi => "Port azimuthal angle",
    :port_radius => "Port position radius", :custom => "Beam is a custom beam",
    :fnumber => "Beam cone F numuber", :transparency => "Beam cone transparency",
    :method => "Pointing method id", :port_distance => "Distance to port - distance",
    :geometry => "Coord sys type for pointing", :p1 => "LPRP position 1", :p2 => "LPRP position 2",
    :p3 => "LPRP position 3", :focus => "Offset from LPRP", :rg_model => "Red-green power spec id",
    :on_blue => "Blue cone power is on", :on_green => "Green cone power is on", :on_red => "Red cone power is on",
    :vis_blue => "Blue cone is visible", :vis_green => "Green cone is visible", :vis_red => "Red cone is visible",
    :mode_blue => "Power selection mode", :mode_green => "Green power selection mode",
    :mode_red => "Red power selection mode", :power_blue => "Blue beam power param",
    :power_green => "Green beam power param", :power_red => "Red beam power param",
    :random_mult => "Beam power random multiplier", :energy_norm => "Normalize beam energy",
    :energy => "Beam energy", :use_mult_blue => "Use laser power multiplier",
    :use_mult_green => "Use green laser power multiplier", :use_mult_red => "Use red laser power multiplier",
    :mult_blue => "Laser power multiplier", :mult_green => "Green laser power multiplier",
    :mult_red => "Red laser power multiplier", :profile => "Spatial power profile model",
    :size_model => "Spot size define method", :major => "Spot size major radius",
    :elliptical => "Beam profile is elliptical", :minor => "Spot size minor radius",
    :rotation => "Beam profile rotation angle", :rot_model => "Beam rotation angle method",
    :order => "Supergauss exponent", :envelope_ratio => "Env rad over spot size at best focus",
    :lens_ratio => "Lens rad over spot size at final optic", :dpr => "DPR beam split is on")

const _VRJ_TABLES = ("Laser power table" => "Laser power table values",
                     "Green cone power table" => "Green laser power table values",
                     "Red cone power table" => "Red laser power table values")

function _vrj_beam(b, useran, path)
    name = String(strip(_vrj_str(b["Laser beam name"])))
    get1(k) = haskey(b, k) ? b[k] : throw(ArgumentError("$path: beam $name has no \"$k\""))
    vals = Dict{Symbol,Any}(c => _vrj_num(get1(k)) for (c, k) in _VRJ_BEAM_KEYS)
    rgb = _vrj_num.(get1("Beam cone color"))
    vals[:red], vals[:green], vals[:blue] = rgb
    vals[:cone_min] = 0.0
    vals[:cone_max] = atand(1 / (2 * vals[:fnumber]))
    vals[:rot_id] = String(strip(_vrj_str(get(b, "Beam rotation angle orientation id", ""))))
    row = NamedTuple{_VRW_BEAM_COLUMNS}(Tuple(vals[c] for c in _VRW_BEAM_COLUMNS))
    # Power tables are stored under the labels of the text reader, with one string per
    # row and its columns joined by spaces, as the text format lists them. The sample
    # files have no power tables, so the column keys ("<table>-Values Col c", as for the
    # mesh connectivity) are an *assumption*. The text reader keeps a label for each table
    # line of its file, empty or not. The JSON file lists all three tables, so only
    # non-empty tables are kept, and the two readers then agree on whether a beam has a
    # power table.
    tables = Dict{String,Vector{String}}()
    for (k, label) in _VRJ_TABLES
        n = Int(_vrj_num(get(b, "$k-Num rows", 0)))
        n == 0 && continue
        nc = Int(_vrj_num(get(b, "$k-Num cols", 0)))
        cols = [get(b, "$k-Values Col $c", nothing) for c in 1:nc]
        all(c -> c isa AbstractVector && length(c) == n, cols) ||
            throw(ArgumentError("$path: beam $name: $k has $n rows but no \"$k-Values Col 1..$nc\" columns"))
        tables[label] = [join((_vrj_str(c[r]) for c in cols), " ") for r in 1:n]
    end
    return merge((; name, cluster = String(_vrj_str(get1("Cluster name"))),
                  ref = String(strip(_vrj_str(get(b, "Reference coord sys name", "Target Chamber")))),
                  port_id = String(strip(_vrj_str(get(b, "Dist to port - port id", "")))),
                  pulse = String(strip(_vrj_str(get(b, "Pulse shape name", "")))), tables, random_used = useran), row)
end

# Returns (version, objects, clusters, rows) of a JSON workspace, as the text reader
# builds them.
function _vrj_workspace(d, path)
    version = String(_vrj_str(get(get(d, "Header data", Dict()), "VisRad version", "")))
    tc = get(d, "Target component data", nothing)
    blocks = tc === nothing ? [] : [_vrj_block(o, path) for o in _vrj_elements(tc, "Target component element")]
    if tc !== nothing
        n = get(tc, "Number of target components", length(blocks))
        n == length(blocks) || throw(ArgumentError("$path: $n target components declared, $(length(blocks)) found"))
    end
    objects = isempty(blocks) ? VisRadObject[] : _vrw_build_objects(blocks, _vrj_folders(d))
    lb = get(d, "Laser beam data", nothing)
    lb === nothing && return version, objects, String[], NamedTuple[]
    clusters = String[String(_vrj_str(c)) for c in get(lb, "Cluster names list", [])]
    useran = _vrj_num(get(lb, "Use random beam pwr mult", 0)) != 0
    rows = NamedTuple[_vrj_beam(b, useran, path) for b in _vrj_elements(lb, "Laser beam element")]
    n = get(lb, "Total number of laser beams", length(rows))
    n == length(rows) || throw(ArgumentError("$path: $n laser beams declared, $(length(rows)) found"))
    return version, objects, clusters, rows
end

# The power-source columns of the JSON format ("<label>-Values") in the order of the text
# table, with the names and units of the text reader.
const _VRJ_PS_COLUMNS = (
    "Direct laser pwr src" => (:direct_laser_pwr_source, "TW/cm2"), "Uniform pwr src" => (:uniform_source, "TW/cm2"),
    "Albedo" => (:albedo, ""), "Xray conv eff" => (:x_ray_conv_efficiency, ""),
    "Laser reflectivity" => (:laser_reflectivity, ""), "Pre-reflect laser pwr" => (:pre_reflect_laser_pwr, "TW/cm2"),
    "Num red cones hit surf" => (:n_red_cones_hitting_surf, ""), "Num blue cones hit surf" => (:n_blue_cones_hitting_surf, ""),
    "Min angle overlap 1w" => (:min_angle_overlap_1w, "deg"), "Min angle overlap 3w" => (:min_angle_overlap_3w, "deg"))
const _VRJ_PS_INT = (:n_red_cones_hitting_surf, :n_blue_cones_hitting_surf)

# `visrad_power_sources` of a JSON workspace. Its "Power source data" lists one
# "Power source element[k]" per target component, in the order of the components, since
# the objects have no names there.
function _vrj_power_sources(d, path, stale)
    ps = _vrj_section(d, "Power source data", path)
    if _vrj_num(get(ps, "Num pwr src up to date", 0)) == 0
        stale === :nothing && return nothing
        throw(ArgumentError("$path: VisRad's power sources are not up to date (compute them in VisRad and save the workspace)"))
    end
    times = _vrj_elements(ps, "Sim time element")
    isempty(times) && throw(ArgumentError("$path: no simulation time in \"Power source data\""))
    length(times) == 1 || @warn "$path: power sources for $(length(times)) simulation times; reading the first"
    st = times[1]
    hd = get(d, "Header data", Dict())
    ht = _vrj_elements(hd, "Sim time element")
    time = isempty(ht) ? 0.0 : _vrj_num(get(ht[1], "Simulation time", 0))
    names = [String(strip(_vrj_str(o["Object name"]))) for o in
             _vrj_elements(_vrj_section(d, "Target component data", path), "Target component element")]
    elems = _vrj_elements(st, "Power source element")
    length(elems) == length(names) ||
        throw(ArgumentError("$path: $(length(elems)) power-source objects for $(length(names)) target components"))
    tables = Dict{String,NamedTuple}()
    units = Dict{Symbol,String}()
    for (name, e) in zip(names, elems)
        n = Int(e["Num surface elements"])
        cols = Symbol[:surface_index]
        vals = Any[collect(0:n-1)]
        for (label, (col, unit)) in _VRJ_PS_COLUMNS
            haskey(e, "$label-Values") || continue
            v = _vrj_num.(e["$label-Values"])
            length(v) == n || throw(ArgumentError("$path: $name: $(length(v)) $label values for $n surfaces"))
            push!(cols, col)
            push!(vals, col in _VRJ_PS_INT ? round.(Int, v) : v)
            units[col] = unit
        end
        units[:surface_index] = ""
        tables[name] = NamedTuple{Tuple(cols)}(Tuple(vals))
    end
    return (; time, units, tables)
end
