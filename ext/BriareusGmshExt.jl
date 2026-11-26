# The STEP reader of Briareus (`read_step`), loaded with Gmsh. OpenCASCADE, inside Gmsh,
# reads the file with its units, Gmsh meshes the faces, and the type and parameters of
# each face are taken from the B-rep. Everything else is in src/step.jl.
module BriareusGmshExt

using Briareus
using Briareus: Vec3, StepSolid, _step_solid
using Gmsh: Gmsh, gmsh

function _read_step(path::AbstractString; size::Real = 0.5, minsize::Real = size / 10, angle::Real = 30.0,
                    toward = (0.0, 0.0, 0.0), tol::Real = 1e-9)
    size > 0 || throw(ArgumentError("size must be positive"))
    0 < minsize <= size || throw(ArgumentError("minsize must lie between 0 and size"))
    0 < angle < 90 || throw(ArgumentError("angle must be between 0 and 90 degrees"))
    started = Gmsh.initialize(["-v", "0"]; finalize_atexit = false)
    try
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.option.setNumber("General.Verbosity", 0)
        gmsh.model.add("briareus_step")
        gmsh.option.setString("Geometry.OCCTargetUnit", "MM")
        gmsh.model.occ.importShapes(String(path), false)   # sheets too
        gmsh.model.occ.synchronize()
        gmsh.option.setNumber("Mesh.MeshSizeMax", Float64(size))
        gmsh.option.setNumber("Mesh.MeshSizeMin", Float64(minsize))
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", ceil(360 / angle))
        # sizes from the edges of a small fillet would otherwise spread over the faces it joins
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.ElementOrder", 1)
        gmsh.model.mesh.generate(2)
        ntags, xyz, _ = gmsh.model.mesh.getNodes()
        index = Dict{Int,Int}(Int(t) => i for (i, t) in enumerate(ntags))
        point(t) = (i = index[Int(t)]; Vec3(xyz[3i-2], xyz[3i-1], xyz[3i]))
        # each solid, then each face outside every solid
        groups = Tuple{String,Vector{Int}}[]
        bounding = Set{Int}()
        for (_, v) in gmsh.model.getEntities(3)
            fs = [abs(Int(t)) for (_, t) in gmsh.model.getBoundary([(3, v)], false, false, false)]
            union!(bounding, fs)
            push!(groups, (_label(gmsh.model.getEntityName(3, v), "solid$(length(groups) + 1)"), fs))
        end
        nsheet = 0
        for (_, s) in gmsh.model.getEntities(2)
            Int(s) in bounding && continue
            nsheet += 1
            push!(groups, (_label(gmsh.model.getEntityName(2, s), "sheet$nsheet"), [Int(s)]))
        end
        seen = Dict{String,Int}()
        solids = StepSolid[]
        for (name, fs) in groups
            n = seen[name] = get(seen, name, 0) + 1
            P = Vec3[]
            local_id = Dict{Int,Int}()
            vid(t) = get!(() -> (push!(P, point(t)); length(P)), local_id, Int(t))
            tris = NTuple{3,Int}[]
            ranges = UnitRange{Int}[]
            info = Tuple{Int,String,Vector{Float64}}[]
            for s in fs
                k0 = length(tris) + 1
                types, _, nodes = gmsh.model.mesh.getElements(2, s)
                for (et, nd) in zip(types, nodes)
                    et == 2 || throw(ArgumentError("face $s of $path was meshed with elements of type $et, not triangles"))
                    for j in 1:3:length(nd)
                        push!(tris, (vid(nd[j]), vid(nd[j+1]), vid(nd[j+2])))
                    end
                end
                push!(ranges, k0:length(tris))
                _, reals = gmsh.model.getEntityProperties(2, s)
                push!(info, (s, String(gmsh.model.getType(2, s)), Float64.(reals)))
            end
            push!(solids, _step_solid(n == 1 ? name : "$name#$n", P, tris, ranges, info; toward, tol))
        end
        return solids
    finally
        gmsh.model.remove()
        started && gmsh.finalize()
    end
end

# The last segment of an OpenCASCADE label such as "Shapes/assembly/part", or `default`
# for an empty label.
function _label(name::AbstractString, default::AbstractString)
    s = strip(last(split(name, '/')))
    return isempty(s) ? String(default) : String(s)
end

__init__() = (Briareus._STEP_READER[] = _read_step)

end
