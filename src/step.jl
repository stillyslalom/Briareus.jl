# STEP import. The reader is the package extension BriareusGmshExt, which loads with
# Gmsh. Gmsh is a weak dependency under the GPL, so the default install of Briareus holds
# none of it. This file holds the solids that the reader returns, their exact surfaces
# and the components built from them.

"""
    StepFace

One face of the boundary representation of a solid read by [`read_step`](@ref):

- `tag`: the number of the face in the STEP model (as Gmsh numbers it)
- `kind`: `:plane`, `:cylinder`, `:cone`, `:sphere`, `:torus`, or `:other` for the
  free-form faces (B-spline, Bézier, revolved, extruded and offset surfaces)
- `type`: the surface type as OpenCASCADE names it, such as `"Cylinder"` or
  `"BSpline surface"`
- `surface`: the exact `CylinderSurface`, `ConeSurface`, `SphereSurface` or
  `TorusSurface` of a curved face, from the parameters stored in the file (`nothing` for
  planes and free-form faces)
- `patches`: the indices of its patches in the mesh of the solid
"""
struct StepFace
    tag::Int
    kind::Symbol
    type::String
    surface::Union{Nothing,CylinderSurface,ConeSurface,SphereSurface,TorusSurface}
    patches::Vector{Int}
end

"""
    StepSolid

A solid, or a sheet outside every solid, read from a STEP file by [`read_step`](@ref):

- `name`: the label of the solid in the file (made unique by `#2`, `#3`, … where
  labels repeat)
- `mesh`: its triangles, one patch per triangle, with the windings made consistent by
  [`orient_mesh`](@ref), so a solid faces out of its material
- `faces`: its B-rep faces ([`StepFace`](@ref))
- `regions`: the faces that lie on one plane or analytic surface, merged where several
  faces share it (such as the two halves of a split cylinder), as the
  [`SurfaceRegion`](@ref)s that [`fit_surfaces`](@ref) would give, by decreasing area.
  Free-form faces are in no region.

[`step_components`](@ref) gives the components of the radiation model, and
[`cad_target`](@ref) takes solids as parts.
"""
struct StepSolid
    name::String
    mesh::Mesh
    faces::Vector{StepFace}
    regions::Vector{SurfaceRegion}
end

function Base.show(io::IO, s::StepSolid)
    kinds = Dict{Symbol,Int}()
    for f in s.faces
        kinds[f.kind] = get(kinds, f.kind, 0) + 1
    end
    order = (:plane, :cylinder, :cone, :sphere, :torus, :other)
    desc = join(("$(kinds[k]) $k" for k in order if haskey(kinds, k)), ", ")
    print(io, "StepSolid(\"", s.name, "\", ", npatches(s.mesh), " patches, ", length(s.faces), " faces: ", desc, ")")
end

# Set by the extension BriareusGmshExt when Gmsh is loaded.
const _STEP_READER = Ref{Any}(nothing)

"""
    read_step(path; size = 0.5, minsize = size / 10, angle = 30.0, toward = (0, 0, 0),
              tol = 1e-9) -> Vector{StepSolid}

The solids of a STEP file (AP203, AP214 or AP242), each meshed into triangles with the
exact analytic surface of every curved face. The file is read by OpenCASCADE through the
package Gmsh, loaded with `using Gmsh`. The reader is a package extension of Briareus,
and Gmsh is not installed with Briareus, since it is under the GPL. Without Gmsh the
function throws an `ArgumentError`.

Coordinates are in mm whatever the units of the file, since STEP stores its units and
OpenCASCADE converts them. Every face is meshed by Gmsh with triangles of at most `size`
mm, and finer on curved faces, with one triangle per `angle` degrees of turn about the
smaller radius of curvature, but none under `minsize` mm. The size on each face is set
by that face alone, so the fine triangles of a small fillet do not spread over the walls
that it joins. Fine meshing is only for the facets, which ray casting uses, since the
radiation model integrates over the exact surfaces whatever the size of the elements. A
thick hohlraum with filleted inner corners (radii 2.9 and 3.2 mm, 10 mm long, 0.25 mm
fillets, 1.6 mm entrance holes) gives 7,580 elements at the defaults, 2,094 of them on
the fillets. The triangles of neighboring faces share their edges. Each triangle is one patch
of the mesh, and [`coarsen_mesh`](@ref) merges them where fewer elements are wanted.
Windings are made consistent by [`orient_mesh`](@ref) with `toward` and `tol`, so a
closed solid faces out of its material and an open sheet faces `toward`.

The type and parameters of each face come from the file: plane coefficients; the axis
and radius of a cylinder; the axis, reference radius and half-angle of a cone; the
center and radius of a sphere; and the axis and radii of a torus. The curved faces then
carry their exact surfaces ([`StepFace`](@ref)), and the components built from them
([`step_components`](@ref), [`cad_target`](@ref)) need no [`fit_surfaces`](@ref). The
1ω curvature classification ([`talis_check`](@ref)) and the radiation model use the
exact surfaces. Free-form faces, such as B-spline fillets, stay as flat triangles with
the curvature estimated from the mesh.

Each solid of the file gives one [`StepSolid`](@ref), named by its label. The faces outside
every solid (sheets, such as a foil drawn as a surface) give one more each, named by
their label or `sheet<k>`. The defaults of `size`, `minsize` and `angle` are engineering
defaults (research notes, "STEP import").
"""
function read_step(path::AbstractString; kwargs...)
    reader = _STEP_READER[]
    reader === nothing &&
        throw(ArgumentError("read_step needs the package Gmsh: install it (`using Pkg; Pkg.add(\"Gmsh\")`) " *
                            "and load it with `using Gmsh`, which loads the STEP reader of Briareus"))
    isfile(path) || throw(ArgumentError("no file $path"))
    return reader(path; kwargs...)::Vector{StepSolid}
end

"""
    step_components(solid, material; two_sided = false, kinds = (:cylinder, :cone, :sphere, :torus),
                    minfraction = 0.0, maxradius = Inf) -> Vector{Component}

The components of a solid read by [`read_step`](@ref): one per analytic surface that its
faces lie on, with the exact surface, and the rest, as [`curved_components`](@ref) splits
a fitted CAD part. They are named `<name>/cyl1`, …, `<name>/rest`, and their `part` is
the name of the solid. Only the surfaces of `kinds` are given to components, and
`minfraction` and `maxradius` leave small or gently curved faces in the rest, as in
`curved_components`. `minfraction` defaults to 0 here, so the exact surface is given to every
curved face, since no fit is needed. `material` and `two_sided` are those of every
component.
"""
function step_components(s::StepSolid, material::Material; two_sided::Bool = false,
                         kinds = (:cylinder, :cone, :sphere, :torus), minfraction::Real = 0.0,
                         maxradius::Real = Inf)
    c = Component(Symbol(s.name), s.mesh, material; two_sided)
    return curved_components(c; regions = s.regions, kinds, minfraction, maxradius)
end

# Returns the exact surface model of a face from the type and the real parameters that
# Gmsh gives for it (`getType`, `getEntityProperties`), or `nothing` for a free-form face.
function _step_model(type::AbstractString, r::AbstractVector{<:Real})
    v(i) = Vec3(r[i], r[i+1], r[i+2])
    if type == "Plane" && length(r) >= 4
        n = v(1)
        L = norm(n)
        L > 0 || return nothing
        return _SurfaceModel(:plane, -r[4] * n / L^2, n / L, Inf)
    elseif type == "Cylinder" && length(r) >= 7
        return _SurfaceModel(:cylinder, v(1), normalize(v(4)), Float64(r[7]))
    elseif type == "Cone" && length(r) >= 8
        # the radius is r[7] at the point v(1) and grows along the axis by tan(half-angle)
        return _SurfaceModel(:cone, v(1), normalize(v(4)), Float64(r[7]), tan(r[8]))
    elseif type == "Sphere" && length(r) >= 4
        return _SurfaceModel(:sphere, v(1), zero(Vec3), Float64(r[4]))
    elseif type == "Torus" && length(r) >= 8
        return _SurfaceModel(:torus, v(1), normalize(v(4)), Float64(r[7]), Float64(r[8]))
    end
    return nothing
end

# Builds a StepSolid from the mesh that the reader gives: the vertices P, the triangles,
# the triangles of each face (`facetris`, consecutive ranges) and, per face, its tag, its
# type and the real parameters of its surface.
function _step_solid(name::AbstractString, P::Vector{Vec3}, tris::Vector{NTuple{3,Int}},
                     facetris::Vector{UnitRange{Int}}, info; toward = (0.0, 0.0, 0.0), tol = 1e-9)
    mesh = Mesh(P, tris, [k:k for k in eachindex(tris)])
    mesh = orient_mesh(mesh; toward, tol)
    models = [_step_model(type, r) for (_, type, r) in info]
    faces = [StepFace(tag, m === nothing ? :other : m.kind, String(type),
                      m === nothing ? nothing : _model_surface(m), collect(facetris[i]))
             for (i, ((tag, type, _), m)) in enumerate(zip(info, models))]
    isempty(P) && return StepSolid(String(name), mesh, faces, SurfaceRegion[])
    L = max(norm(maximum(P) - minimum(P)), eps())
    # merge the faces that lie on one surface: same kind and radius, and the vertices of
    # one within the tolerance of the surface of the other
    t = sqrt(tol) * L
    group = collect(eachindex(faces))
    for j in eachindex(faces), i in 1:j-1
        group[i] == i || continue
        (mi, mj) = (models[i], models[j])
        (mi === nothing || mj === nothing || mi.kind !== mj.kind) && continue
        isfinite(mi.R) && abs(_model_radius(mi) - _model_radius(mj)) > t && continue
        mi.kind === :plane && abs(abs(dot(mi.u, mj.u)) - 1) > sqrt(tol) && continue
        if all(abs(_model_distance(mi, P[v])) <= t for k in facetris[j] for v in tris[k])
            group[j] = i
            break
        end
    end
    D = _fit_data(mesh, tol)
    scratch = falses(length(D.P))
    regs = SurfaceRegion[]
    for i in eachindex(faces)
        (group[i] == i && models[i] !== nothing) || continue
        members = reduce(vcat, (faces[j].patches for j in eachindex(faces) if group[j] == i); init = Int[])
        isempty(members) && continue
        push!(regs, _region(models[i].kind, members, models[i], D, scratch))
    end
    sort!(regs; by = r -> -r.area)
    return StepSolid(String(name), mesh, faces, regs)
end
