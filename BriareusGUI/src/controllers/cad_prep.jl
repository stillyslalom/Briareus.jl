# This file covers the guided preparation of a target imported from CAD: its parts split
# into their connected pieces, the roles that the geometry suggests, and the curved
# surfaces fitted to the parts that take a role.

"""
Steps of the guided preparation of a CAD target ([`cad_prep_rows`](@ref),
[`split_cad_part!`](@ref), [`suggest_cad_roles`](@ref), [`apply_cad_roles!`](@ref),
[`curvature_preview`](@ref)): `:split`, orient the STL parts and split each into its
connected pieces; `:detect`, detect cavities, capsules, foils and hardware; `:assign`,
give the roles and preview the curved surfaces fitted to the parts that take one.
"""
const CAD_PREP_STEPS = (:split, :detect, :assign)

const CAD_PREP_LABELS = Dict(:split => "1. Orient and split", :detect => "2. Detect cavities, foils and hardware",
                             :assign => "3. Assign roles and preview curvature fits")

"""
    cad_prep_rows(ts::TargetSetup) -> Vector{NamedTuple}

The STL parts of the design (kind `:mesh`), one row `(; index, name, file, facets,
pieces)` each: the part index, its name, the file it comes from, its facets and the
number of its connected pieces (`Briareus.connected_pieces`, vertices welded). The
facets of each part face out of a closed solid, or toward the target center for an open
sheet, as the part is meshed (`Briareus.orient_mesh`). A part that cannot be read gives
`facets = pieces = 0`.
"""
function cad_prep_rows(ts::TargetSetup)
    out = NamedTuple{(:index, :name, :file, :facets, :pieces),Tuple{Int,Symbol,String,Int,Int}}[]
    for (i, p) in enumerate(ts.parts[])
        p.kind === :mesh || continue
        n, k = try
            m = _stl_facets(p.file, p.object, p.values[:scale])
            occursin('#', p.object) ? (npatches(m), 1) : (npatches(m), length(_stl_pieces(p.file, p.object, m)))
        catch err
            err isa ArgumentError || err isa SystemError || rethrow()
            (0, 0)
        end
        push!(out, (; index = i, name = p.name, file = basename(p.file), facets = n, pieces = k))
    end
    return out
end

"""
    split_cad_part!(ts::TargetSetup, i) -> Bool

Replace STL part `i` by one part per connected piece of its facets, as a multi-body CAD
export would give them: each named `<part>_<k>` (made unique), with the placement, size,
scale and flags of the part, as members of a new assembly named after the part in the
assembly of the part. A piece keeps no role, so roles are given afterward
([`apply_cad_roles!`](@ref)). A part of one piece, a part other parts are placed in and
anything the geometry rejects leave the design as it was, set `error` and return
`false`.
"""
function split_cad_part!(ts::TargetSetup, i::Integer)
    _check_part_index(ts, i)
    ps = ts.parts[]
    p = ps[i]
    try
        p.kind === :mesh || throw(ArgumentError("only an STL part is split into pieces"))
        occursin('#', p.object) && throw(ArgumentError("$(p.name) is one piece already"))
        ownedby(f) = Symbol(first(split(String(f), '.'; limit = 2))) === p.name     # a piece frame's part
        any(q -> q !== p && (ownedby(q.frame) || any(st -> any(ownedby, Briareus.GUIAPI.step_references(st)), q.steps)), ps) &&
            throw(ArgumentError("other parts are placed in $(p.name): place them elsewhere first"))
        m = _stl_facets(p.file, p.object, p.values[:scale])
        n = length(_stl_pieces(p.file, p.object, m))
        n > 1 || throw(ArgumentError("$(p.name) is one connected piece"))
        taken = setdiff(union(Set(all_component_names(ts)), (q.name for q in ps), keys(part_frames(ts))),
                        (p.name,))
        asm = p.name
        push!(taken, asm)
        pieces = [_with(p; name = _free_identifier(string(p.name, "_", k), taken, "piece"),
                        object = string(p.object, "#", k), assembly = asm,
                        foil = false, observed = false, wall = false) for k in 1:n]
        new = vcat(ps[1:i-1], pieces, ps[i+1:end],
                   [PartSpec(:assembly; name = asm, frame = assembly_frame(p.assembly), assembly = p.assembly)])
        return _try_build!(ts; parts = new)
    catch err
        err isa ArgumentError || rethrow()
        ts.error[] = _errmsg(err)
        return false
    end
end

"""
    suggest_cad_roles(sess::DesignSession; kwargs...) -> (; rows, cavities, text)

The roles that the geometry of the target suggests for the parts of the design
(`Briareus.cad_roles`, with its keywords; the cavities are looked for at target chamber
center and at the center of each part that encloses it). One row per part that takes a
role in the window, `(; index, name, role, kind, cavity, current)`: the suggested role
(one of [`PART_ROLES`](@ref)), the kind of shape (`:capsule`, `:wall`, `:sheet`,
`:solid` or `:hardware`), the cavity it bounds (0 for none) and the role it has now.
Hohlraum parts are left out, since their pieces hold their roles. `text` lists the
cavities and the rows.
"""
function suggest_cad_roles(sess::DesignSession; kwargs...)
    ts = sess.target
    t = ts.target[]
    r = cad_roles(t; kwargs...)
    rows = NamedTuple{(:index, :name, :role, :kind, :cavity, :current),Tuple{Int,Symbol,Symbol,Symbol,Int,Symbol}}[]
    for x in r.parts
        i = part_index(ts, x.part)
        i === nothing && continue
        p = ts.parts[][i]
        p.kind === :hohlraum && continue
        push!(rows, (; index = i, name = x.part, role = x.role, kind = x.kind, cavity = x.cavity, current = part_role(p)))
    end
    io = IOBuffer()
    if isempty(r.cavities)
        println(io, "no cavity found")
    else
        for (k, c) in enumerate(r.cavities)
            @printf(io, "cavity %d at (%.3g, %.3g, %.3g) mm, %.0f%% enclosed: %s\n", k, c.center..., 100 * c.enclosed,
                    join(c.walls, ", "))
        end
    end
    for x in rows
        @printf(io, "%s: %s, suggested %s%s\n", x.name, x.kind, PART_ROLE_LABELS[x.role],
                x.current === x.role ? "" : " (now " * PART_ROLE_LABELS[x.current] * ")")
    end
    return (; rows, r.cavities, text = String(take!(io)))
end

"""
    apply_cad_roles!(sess::DesignSession, rows = suggest_cad_roles(sess).rows) -> Vector{String}

Give each part of `rows` its suggested role ([`set_part_role!`](@ref)), the observed role
last since one part holds it at a time. A CAD part suggested as the observed capsule gets
a sphere fitted to it ([`fit_sphere`](@ref)). Returns one line per part whose role could
not be given, with the reason; the other parts keep their new roles.
"""
function apply_cad_roles!(sess::DesignSession, rows = suggest_cad_roles(sess).rows)
    ts = sess.target
    failed = String[]
    for x in sort(collect(rows); by = x -> x.role === :observed)
        i = part_index(ts, x.name)
        i === nothing && continue
        part_role(ts.parts[][i]) === x.role && continue
        set_part_role!(ts, i, x.role) || push!(failed, string(x.name, ": ", ts.error[]))
    end
    return failed
end

"""
    curvature_preview(sess::DesignSession; roles = (:observed, :foil, :wall)) -> Vector{NamedTuple}

The curved surfaces of the parts of the target with one of `roles`, as the drive fits
them ([`FIT_CAD`](@ref): `Briareus.curved_components` on the file facets of each CAD part,
whose vertices lie on the true surfaces): one row `(; part, cylinders, cones, spheres,
tori, flat)` per part, the count of each kind of curved piece (a sphere or cylinder part
is one piece of its kind) and the share of the area of the part left flat.
"""
function curvature_preview(sess::DesignSession; roles = (:observed, :foil, :wall))
    ts = sess.target
    t = ts.target[]
    ps = [p for p in ts.parts[] if part_role(p) in roles]
    named = Set(p.name for p in ps)
    refs = filter(kv -> first(kv) in named, _fit_references(t, ts.mesh[], ts.parts[]))
    tf, _ = _fitted_target(t, refs)
    out = NamedTuple{(:part, :cylinders, :cones, :spheres, :tori, :flat),Tuple{Symbol,Int,Int,Int,Int,Float64}}[]
    for n in unique(c.part for c in tf.components if c.part in named && c.kind !== :aperture)
        cs = [c for c in tf.components if c.part === n && c.kind !== :aperture]
        kinds = [c.surface for c in cs]
        area = [sum(_mesh_areas(c.mesh)) for c in cs]
        flat = sum(area[k] for k in eachindex(cs) if kinds[k] === nothing; init = 0.0)
        push!(out, (; part = n, cylinders = count(s -> s isa CylinderSurface, kinds),
                    cones = count(s -> s isa ConeSurface, kinds), spheres = count(s -> s isa SphereSurface, kinds),
                    tori = count(s -> s isa TorusSurface, kinds), flat = flat / max(sum(area), eps())))
    end
    return out
end

# The area of each triangle of a mesh.
_mesh_areas(m::Mesh) = [norm(cross(m.vertices[b] - m.vertices[a], m.vertices[c] - m.vertices[a])) / 2
                        for (a, b, c) in m.triangles]
