# A synthetic VisRad workspace for the tests (core and GUI): object blocks and beam rows
# in the layout VisRad 20/21 writes, copied from real workspaces with made-up names.
#
#   synthetic_workspace() -> String (the file's text)

# the CAD object: a right tetrahedron (cm), centered, outward counterclockwise faces
# (0-based node indices)
const TETRA_NODES = [Tuple(v .- 0.025) for v in [(0.0, 0.0, 0.0), (0.1, 0.0, 0.0), (0.0, 0.1, 0.0), (0.0, 0.0, 0.1)]]
const TETRA_FACES = [(0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3)]

function synthetic_workspace()
    # A workspace in the layout VisRad 20/21 writes (Workspace Format ID 11), with
    # object blocks and beam rows copied from real workspaces: the primitive hohlraum
    # + foil of a VisRad 21.0.0 file, a CAD mesh placed through a reference cone, and
    # NIF beam rows. Lengths are cm in the file.
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
