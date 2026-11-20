# A synthetic VisRad workspace for the tests (core and GUI): object blocks and beam rows
# in the layout VisRad 20/21 writes, copied from real workspaces with made-up names.
#
#   synthetic_workspace(; meshes = ()) -> String (the file's text)

# The CAD object is a centered right tetrahedron in cm, with outward counterclockwise
# faces and 0-based node indices.
const TETRA_NODES = [Tuple(v .- 0.025) for v in [(0.0, 0.0, 0.0), (0.1, 0.0, 0.0), (0.0, 0.1, 0.0), (0.0, 0.0, 0.1)]]
const TETRA_FACES = [(0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3)]

function synthetic_workspace(; meshes = ())
    # A workspace in the layout VisRad 20/21 writes (Workspace Format ID 11), with
    # object blocks and beam rows copied from real workspaces: the primitive hohlraum
    # + foil of a VisRad 21.0.0 file, a CAD mesh placed through a reference cone, and
    # NIF beam rows. Lengths are cm in the file. `meshes`: further CAD objects after
    # these, `(name, nodes, faces, position)` (cm; nodes relative to the position,
    # 0-based faces), placed in the chamber frame.
    function vobj(io, name, type; pos = (0, 0, 0), geom = 0, angles = (0, 0, 0),
                  lengths = zeros(6), pts = (0, 0, 0), ref = "Target Chamber", albedo = 0.7,
                  include = 1, mesh = nothing)
        println(io, "Begin Parameters for Object = ", name)
        println(io, "Object type                 = ", type, "\n")
        for k in 1:6
            println(io, "   Length $k model ID              = 0")
            println(io, "   Length $k value                 = ", lengths[k])
            println(io, "   Length $k string                = ")
        end
        println(io, "   Use non-zero wall thickness    = 0")
        for k in 1:3
            println(io, "   Number of mesh pts $k           = ", pts[k])
        end
        println(io, "   Normals in/out                 = 0")
        println(io, "   Red value                      = 1\n   Green value                    = 0.5\n   Blue value                     = 0")
        println(io, "   Object visibility              = 1")
        println(io, "   Reference coord. system type   = ", ref == "Target Chamber" ? 0 : 1)
        println(io, "   Reference coord. system name   = ", ref)
        println(io, "   Set position of ID             = 0")
        println(io, "   Coord. sys. geometry type      = ", geom)
        for k in 1:3
            println(io, "   Position $k                     = ", pos[k])
        end
        println(io, "   Tilt angle                     = ", angles[1])
        println(io, "   Rotation angle                 = ", angles[2])
        println(io, "   Spin angle                     = ", angles[3])
        println(io, "   Albedo model ID                = 0\n   Albedo value                   = ", albedo)
        println(io, "   Xray conv. effic. model ID     = 0\n   Xray conv. effic. value        = 0.6")
        println(io, "   Laser reflectivity model ID    = 0\n   Laser reflectivity value       = 0")
        println(io, "   Include in dir laser dep calc. = ", include)
        println(io, "   Include in ref laser dep calc. = ", include)
        println(io, "   Include in radiosity calc.     = ", include)
        println(io, "Time-Dependent Albedo Values        for <$name>:  # table times = 0    [format=1]\n")
        println(io, "   Number of Bounding Volumes = 0\n")
        if mesh !== nothing
            nodes, faces = mesh
            println(io, "   Gen. mesh object origin: x     = ", -pos[1])
            println(io, "   Gen. mesh object origin: y     = ", -pos[2])
            println(io, "   Gen. mesh object origin: z     = ", -pos[3])
            println(io, "Number of vertices = ", length(nodes))
            println(io, "Number of surfaces = ", length(faces), "\n\nNode positions:")
            for (i, v) in enumerate(nodes)
                println(io, "  ", i - 1, "  ", join(v, "  "))
            end
            println(io, "Surface element connectivity:")
            for (i, f) in enumerate(faces)
                println(io, "  ", i - 1, "  ", join(f, "  "), "  -1")
            end
        end
        println(io, "\nEnd Parameters for Object: ", name, "\n")
    end
    tn, tf = TETRA_NODES, TETRA_FACES
    # beam rows as VisRad writes them (format 22/23 columns; "rotat. angle ID" blank)
    row(id, θ, φ, method, p, focus, on, model, major, minor, n; rotid = "", rg = 2) =
        "    $id     $θ     $φ    769.500      0      0.000      1.380     20.755     0.137     0.137     0.455     0.500      $method      0.000000      0   $(p[1])    $(p[2])    $(p[3])    $focus    $rg      $on      0      0      0      0      0      0      0      0      1.000     0.050     0.050    1.0000      0    1.000      0      0      0    1.00000    1.00000    1.00000     $model      0  $major      1  $minor    0.00000      0     $rotid         $n    1.53645    0.40000      0"
    io = IOBuffer()
    println(io, "Workspace Format ID  = 11\n[Header Data]:\n   Format ID = 4\n   VisRad Version_21.0.0\n[End Header Data]\n")
    println(io, "[Target Component Data]:\n   Format ID = 2\n#\n# File format ID:  29\n#\nNumber of View Factor Objects =  5\n")
    vobj(io, "Rectangle_0", 4; pos = (0.1, 0, 0), angles = (90, 0, 0), lengths = [0.4, 0.4, 0, 360, 180, 0], pts = (10, 10, 0))
    vobj(io, "Cyl_Hohlraum_0", 5; pos = (-0.4, 0, 0), lengths = [0.255, 0.935, 0.152, 0, 180, 0], pts = (11, 12, 4))
    vobj(io, "CAD_ref", 11; pos = (0, 0, 0.2), angles = (0, 90, 0), lengths = [0.01, 0, 0.02, 360, 180, 0], include = 0)
    vobj(io, "CAD_0", 13; pos = (0.1, 0, 0), ref = "CAD_ref", albedo = 0.5, mesh = (tn, tf))
    vobj(io, "Marker", 1; pos = (0.32, 90, 150), geom = 2, lengths = [0.1, 0, 360, 0, 180, 0], pts = (10, 10, 0), include = 0)
    for (name, nodes, faces, pos) in meshes
        vobj(io, name, 13; pos, mesh = (nodes, faces))
    end
    println(io, """
      Begin Target Component Directory Structure:

        Target Components
          Number of target components = 2
            Cyl_Hohlraum_0
            Marker
          Number of child folders = 2
            Foil
              Number of target components = 1
                Rectangle_0
              Number of child folders = 0
            CAD
              Number of target components = 2
                CAD_ref
                CAD_0
              Number of child folders = 0

      End Target Component Directory Structure

    [End Target Component Data]

    [Laser Beam Data]:
    #
    #VF LaserParmFile: Format ID = 23
    #
    Laser System ID = 1
    use random beam pwr mult      =  0

    Cluster Names:
    Ring_050
    Ring_045
    """)
    for (name, cl, r, ref) in (
            ("B113", "Ring_050", row(1, 52.382, 37.056, 2, (-0.024500, 0.037750, 0.590500), "-0.010000", 1, 52, "620.00000", "360.00000", "4.50000"), "Target Chamber"),
            ("B111", "Ring_050", row(2, 52.382, 42.284, 2, (0.0, 0.0, 0.0), "0.000000", 1, 0, "300.00000", "300.00000", "2.00000"; rotid = "R1", rg = 0), "CAD_0"),
            ("B121", "Ring_045", row(5, 46.888, 59.624, 0, (0.1, 0.2, 0.3), "0.000000", 0, 51, "640.00000", "370.00000", "4.80000"), "Target Chamber"))
        println(io, "Laser Beam Name   =  ", name, "\nCluster Name      = ", cl, "\n", r)
        println(io, "Ref. Coord. Sys.:  ", ref, "\nD toward Port ID:  P   0 -   0\nPulse Shape Name:  ")
        println(io, "  GraphicsItems format ID    = 8\n   Rotation angle       = 0.00000")
        println(io, "Laser power table values:  # table times = 0    [format=1]")
    end
    println(io, "[End Laser Beam Data]\n\n[Diagnostic Port Data]:\n[End Diagnostic Port Data]")
    return String(take!(io))
end

# The [Power Source Data] section VisRad writes after computing the power sources, in
# the layout of real VisRad 20/21 workspaces (header lines copied verbatim): per object,
# `name => direct laser source per element (TW/cm²)`; the other columns are filler.
function power_source_section(objects; up_to_date = true)
    io = IOBuffer()
    println(io, "[Power Source Data]:")
    if !up_to_date
        println(io, "   Power sources are not up to date for time = 0.000000e+00\n[End Power Source Data]\n")
        return String(take!(io))
    end
    println(io, "   Power sources are up to date for time = 0.000000e+00")
    println(io, "#\n#**************** RESULTS FOR POWER SOURCE TERMS **************#")
    println(io, "#    Number of Objects =      ", length(objects), "\n#\n#")
    for (name, q) in objects
        println(io, "#  --- Power Source Summary for ", name, "\n#")
        println(io, "#      Format ID          = 6\n#      Number of surfaces = ", length(q), "\n#\n#")
        println(io, "# Surface     Direct Laser   Uniform      Albedo        X-ray Conv.   Laser         Pre-reflect   # Red Cones   # Blue Cones  Min Angle     Min Angle")
        println(io, "#  Index      Pwr Source     Source                     Efficiency    Reflectivity  Laser Pwr     Hitting Surf  Hitting Surf  Overlap-1w    Overlap-3w")
        println(io, "#             (TW/cm2)       (TW/cm2)                                               (TW/cm2)                                  (deg)         (deg)")
        println(io, "#")
        for (i, x) in enumerate(q)
            println(io, lpad(i - 1, 7), lpad(string(x), 14), "       0.00000       0.70000       0.60000",
                    "       0.00000       0.00000             0", lpad(x > 0 ? 2 : 0, 14), "     360.00000     360.00000")
        end
    end
    println(io, "[End Power Source Data]\n")
    return String(take!(io))
end

# The JSON workspace VisRad 21 saves after computing ("Workspace Format ID": 1001), built
# from a parsed text workspace in the layout of real VisRad 21.0.0 JSON files. Object
# parameters are under the keys of the text format, and "Object color" and "Gen. mesh object
# origin" are arrays. The node array is "Object node array[k]", the connectivity is given as
# columns, and numbers in arrays are strings. Beams are under readable keys at the precision
# `ports` gives; VisRad writes them unrounded, and the text table rounds to 3 decimals. The
# power sources are per object in the order of the components. `sources` gives, per object
# name, the direct laser source per element, or `nothing` for "not up to date".
function synthetic_json_workspace(ws; sources = nothing, ports = Dict{String,NTuple{3,Float64}}())
    q(x) = string(x)
    obj = Dict{String,Any}[]
    for o in ws.objects
        d = Dict{String,Any}("Object name" => o.name)
        for (k, v) in o.params
            k in ("Red value", "Green value", "Blue value", "Number of vertices", "Number of surfaces") && continue
            startswith(k, "Gen. mesh object origin") && continue
            x = tryparse(Float64, v)
            d[k] = x === nothing ? v : isinteger(x) ? Int(x) : x
        end
        d["Object color"] = [o.params[k] for k in ("Red value", "Green value", "Blue value")]
        if o.mesh !== nothing
            Rt, t = transpose(o.rotation), o.position
            nodes = [Rt * (v - t) / 10 for v in o.mesh.vertices]          # cm, object frame
            d["Gen. mesh object origin"] = q.(-t / 10)
            d["Object node array-Num elements"] = length(nodes)
            for (k, v) in enumerate(nodes)
                d["Object node array[$k]"] = q.(round.(collect(v); sigdigits = 12))
            end
            d["Object connectivity array-Num rows"] = length(o.mesh.triangles)
            d["Object connectivity array-Num cols"] = 4
            for c in 1:4
                d["Object connectivity array-Values Col $c"] = [c == 4 ? "-1" : q(t[c] - 1) for t in o.mesh.triangles]
            end
        end
        push!(obj, d)
    end
    tc = Dict{String,Any}("Exploded view reference object" => "Target Chamber", "Number of target components" => length(obj))
    for (k, d) in enumerate(obj)
        tc["Target component element[$k]"] = d
    end
    path(o) = join(filter(!isempty, ["Target Components", o.folder, o.name]), "/")
    beams = Dict{String,Any}("Number of fixed laser beams" => length(ws.rows), "Total number of laser beams" => length(ws.rows),
                             "Use random beam pwr mult" => 0, "Number of laser beam clusters" => length(ws.clusters),
                             "Cluster names list" => ws.clusters)
    for (k, r) in enumerate(ws.rows)
        θ, φ, F = get(ports, r.name, (r.theta, r.phi, r.fnumber))
        b = Dict{String,Any}("Laser beam name" => " " * r.name, "Cluster name" => r.cluster,
                             "Beam cone color" => q.([r.red, r.green, r.blue]),
                             "Beam rotation angle orientation id" => r.rot_id, "Reference coord sys name" => r.ref,
                             "Dist to port - port id" => r.port_id, "Pulse shape name" => r.pulse,
                             "Laser power table-Num rows" => 0, "Laser power table-Num cols" => 2,
                             "Blue cone" => Dict("Main window BSO text box" => Dict("Visibility" => 1)))
        for (c, key) in Briareus._VRJ_BEAM_KEYS
            b[key] = c === :theta ? θ : c === :phi ? φ : c === :fnumber ? F : getfield(r, c)
        end
        beams["Laser beam element[$k]"] = b
    end
    ps = Dict{String,Any}("Num pwr src up to date" => sources === nothing ? 0 : 1)
    if sources !== nothing
        st = Dict{String,Any}("Number of objects" => length(obj))
        for (k, o) in enumerate(ws.objects)
            v = get(sources, o.name, Float64[])
            e = Dict{String,Any}("Num surface elements" => length(v))
            for (label, (col, _)) in Briareus._VRJ_PS_COLUMNS
                vals = col === :direct_laser_pwr_source ? v : col === :albedo ? fill(0.7, length(v)) :
                       col === :n_blue_cones_hitting_surf ? [x > 0 ? 2.0 : 0.0 for x in v] :
                       col in (:min_angle_overlap_1w, :min_angle_overlap_3w) ? fill(360.0, length(v)) : zeros(length(v))
                e["$label-Num elems"] = length(v)
                e["$label-Values"] = [Briareus.@sprintf("%.6e", x) for x in vals]
            end
            st["Power source element[$k]"] = e
        end
        ps["Sim time element[1]"] = st
    end
    return Briareus.JSON.json(Dict{String,Any}(
        "Workspace Format ID" => 1001,
        "Header data" => Dict("VisRad version" => ws.version, "Number of simulation times" => 1,
                              "Sim time element[1]" => Dict("Simulation time" => 0, "Selected for run" => 1)),
        "Target component data" => tc,
        "Target comp list structure" => Dict("Folder name" => ["Target Components"], "Target comp name" => path.(ws.objects)),
        "Laser beam data" => beams,
        "Power source data" => ps,
        "Radiosity data" => Dict("Num radiosity up to date" => 0)); pretty = 3)
end
