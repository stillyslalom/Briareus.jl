# The model tree of parts and assemblies: an assembly is a part of kind `:assembly` (a
# named frame without geometry, placed and swept like a part), and each part's and
# assembly's `assembly` names the assembly it is a member of (`TOP_LEVEL`: none). The
# tree's editing commands keep the members where they are in the chamber. The analyses
# see the flattened components (`_target_components`) and never the tree.

"Whether part `p` is an assembly (an `:assembly` part, see [`PART_KINDS`](@ref))."
is_assembly(p) = p.kind === :assembly

"""
    assemblies(ts::TargetSetup) -> Vector{PartSpec}

The assemblies of a design (its `:assembly` parts), in the parts' order.
"""
assemblies(ts::TargetSetup) = PartSpec[p for p in ts.parts[] if is_assembly(p)]

"""
    check_tree(parts)

Throw `ArgumentError` unless the parts' memberships ([`PartSpec`](@ref)'s `assembly`)
form a tree: each names [`TOP_LEVEL`](@ref) or an assembly of `parts`, no assembly is
inside itself (at any depth), and no user frame is a member (frames stay in the Frames
folder).
"""
function check_tree(parts)
    kind = Dict{Symbol,Symbol}(p.name => p.kind for p in parts)
    parent = Dict{Symbol,Symbol}(p.name => p.assembly for p in parts)
    for p in parts
        a = p.assembly
        a === TOP_LEVEL && continue
        p.kind === :frame && throw(ArgumentError("$(p.name) is a user frame: frames are in no assembly"))
        get(kind, a, nothing) === :assembly || throw(ArgumentError("no assembly named $a to put $(p.name) in"))
        seen = Set{Symbol}((p.name,))
        while a !== TOP_LEVEL
            a in seen && throw(ArgumentError("assembly $a is inside itself (through $(p.name))"))
            push!(seen, a)
            a = parent[a]
        end
    end
    return nothing
end

# Each part's place among its siblings in the tree: its index in the parts, an assembly's
# the least of its own and its members' (at any depth), so an assembly is listed where
# its first member was.
function _tree_ranks(parts)
    rank = Dict{Symbol,Int}(p.name => i for (i, p) in enumerate(parts))
    parent = Dict{Symbol,Symbol}(p.name => p.assembly for p in parts)
    for (i, p) in enumerate(parts)
        a, n = p.assembly, 0
        while a !== TOP_LEVEL && haskey(rank, a) && n <= length(parts)
            rank[a] = min(rank[a], i)
            a, n = parent[a], n + 1
        end
    end
    return rank
end

"""
    assembly_children(parts, name = TOP_LEVEL) -> Vector{PartSpec}

The members of assembly `name` (with [`TOP_LEVEL`](@ref): the parts and assemblies in
none), user frames left out, in the tree's order: an assembly where its first member is
among the parts, the others in the parts' order.
"""
function assembly_children(parts, name::Symbol = TOP_LEVEL)
    rank = _tree_ranks(parts)
    kids = PartSpec[p for p in parts if p.assembly === name && p.kind !== :frame]
    return sort!(kids; by = p -> rank[p.name])
end

"""
    assembly_members(parts, name) -> Vector{Symbol}

Every part and assembly inside assembly `name`, at any depth, depth first in the tree's
order ([`assembly_children`](@ref)).
"""
function assembly_members(parts, name::Symbol)
    out = Symbol[]
    for p in assembly_children(parts, name)
        push!(out, p.name)
        is_assembly(p) && append!(out, assembly_members(parts, p.name))
    end
    return out
end

"""
    assembly_parts(parts, name) -> Vector{Symbol}

The parts with geometry inside assembly `name`, at any depth ([`assembly_members`](@ref)
without the assemblies), in the tree's order.
"""
assembly_parts(parts, name::Symbol) =
    Symbol[n for n in assembly_members(parts, name) if (k = findfirst(p -> p.name === n, parts); !is_assembly(parts[k]))]

"""
    assembly_frame(assembly::Symbol) -> Symbol

The frame a member of `assembly` is placed in by default: the assembly's own, or `:tcc`
at the top level ([`TOP_LEVEL`](@ref)).
"""
assembly_frame(a::Symbol) = a === TOP_LEVEL ? :tcc : a

# The assemblies `a` is inside, itself first, then its parent, … (none for TOP_LEVEL).
function _ancestors(parts, a::Symbol)
    out = Symbol[]
    while a !== TOP_LEVEL && !(a in out)
        push!(out, a)
        k = findfirst(p -> p.name === a, parts)
        k === nothing && break
        a = parts[k].assembly
    end
    return out
end

# `parts` with the parts and assemblies `names` moved into assembly `a` (TOP_LEVEL: out
# of every assembly). A member placed in its old assembly's frame (or in `:tcc`, at the
# top level) is placed in the new one's, where it is in the chamber (`fr`: the chamber
# frames, `part_frames`); one placed in any other frame keeps its placement.
function _moved_into(parts::AbstractVector{PartSpec}, names, a::Symbol, fr::AbstractDict)
    out = copy(parts)
    if a !== TOP_LEVEL
        k = findfirst(p -> p.name === a, out)
        (k !== nothing && is_assembly(out[k])) || throw(ArgumentError("no assembly named $a"))
    end
    up = _ancestors(out, a)
    for n in names
        n = Symbol(n)
        k = findfirst(p -> p.name === n, out)
        k === nothing && throw(ArgumentError("no part or assembly named $n"))
        p = out[k]
        p.kind === :frame && throw(ArgumentError("$n is a user frame: frames are in no assembly"))
        n in up && throw(ArgumentError(n === a ? "an assembly cannot be put in itself" :
                                       "$a is inside $n: an assembly cannot be put in itself"))
        p.assembly === a && continue
        q = p.frame === assembly_frame(p.assembly) ? _kept_in(p, assembly_frame(a), fr) : p
        out[k] = _with(q; assembly = a)
    end
    return out
end

# Run `f()`, an edit that returns whether it was committed; a rejected one (an
# `ArgumentError`) sets `error` and returns `false`.
function _tree_edit(f, ts::TargetSetup)
    try
        return f()
    catch err
        err isa ArgumentError || rethrow()
        ts.error[] = _errmsg(err)
        return false
    end
end

"""
    add_assembly!(ts::TargetSetup; name = nothing, members = [], parent = nothing) -> Bool

Add an assembly named `name` (default `assembly_1`, …) and put the parts and assemblies
`members` (names) in it, where they are in the chamber ([`move_into_assembly!`](@ref)).
It is a member of `parent` (an assembly's name, or [`TOP_LEVEL`](@ref)); by default of
the assembly its members share, else of none. Its frame is its parent's (at its origin,
with its axes: placed in [`assembly_frame`](@ref)`(parent)` by zero offsets), so placing
the members in it moves nothing; it is edited, placed and swept like any part's
afterward ([`edit_part!`](@ref), [`set_part_frame!`](@ref), the placement steps). One
edit: a rejected assembly (a name taken, a member that would be inside itself) adds
nothing, sets `error` and returns `false`.
"""
function add_assembly!(ts::TargetSetup; name = nothing, members = Symbol[], parent = nothing)
    _tree_edit(ts) do
        ps = ts.parts[]
        ms = Symbol[Symbol(m) for m in members]
        if parent === nothing
            owners = unique(p.assembly for p in ps if p.name in ms)
            par = length(owners) == 1 ? only(owners) : TOP_LEVEL
        else
            par = Symbol(parent)
        end
        nm = name === nothing ? _free_name(ts, :assembly) : Symbol(name)
        occursin(r"^[A-Za-z_][A-Za-z0-9_]*$", String(nm)) ||
            throw(ArgumentError("an assembly name is letters, digits and _, starting with a letter"))
        new = [ps; PartSpec(:assembly; name = nm, frame = assembly_frame(par), assembly = par)]
        return _try_build!(ts; parts = _moved_into(new, ms, nm, part_frames(new)))
    end
end

"""
    move_into_assembly!(ts::TargetSetup, names, assembly::Symbol) -> Bool

Make the parts and assemblies `names` members of `assembly` (an assembly's name, or
[`TOP_LEVEL`](@ref) for none). They stay where they are in the chamber, as
[`set_part_frame!`](@ref) keeps a part: a member placed in its old assembly's frame (or
in TCC, at the top level) is placed in the new one's; one placed in another frame (a
capsule in its hohlraum's) keeps its placement. User frames are in no assembly, and an
assembly cannot go inside itself: either sets `error` and returns `false`, and so does
any edit the geometry rejects.
"""
move_into_assembly!(ts::TargetSetup, names, a::Symbol) =
    _tree_edit(() -> _try_build!(ts; parts = _moved_into(ts.parts[], names, a, part_frames(ts))), ts)

"""
    ungroup_assembly!(ts::TargetSetup, name::Symbol) -> Bool

Remove assembly `name` and keep its members: they become members of its parent, and
everything placed in its frame (its members by default) is placed in the parent's
([`assembly_frame`](@ref)), where it is in the chamber. Refused (with `error` set) while a
placement step refers to the assembly's frame; beams aimed in it keep their aims, held
in `:tcc` ([`follow_frames!`](@ref)).
"""
function ungroup_assembly!(ts::TargetSetup, name::Symbol)
    _tree_edit(ts) do
        ps = ts.parts[]
        k = part_index(ts, name)
        (k !== nothing && is_assembly(ps[k])) || throw(ArgumentError("no assembly named $name"))
        a = ps[k]
        refs = [q.name for q in ps if any(st -> name in Briareus._step_refs(st), q.steps)]
        isempty(refs) || throw(ArgumentError("a placement step of $(join(String.(refs), ", ")) refers to $name"))
        fr, into = part_frames(ts), assembly_frame(a.assembly)
        out = PartSpec[]
        for q in ps
            q === a && continue
            q.frame === name && (q = _kept_in(q, into, fr))
            q.assembly === name && (q = _with(q; assembly = a.assembly))
            push!(out, q)
        end
        return _try_build!(ts; parts = out)
    end
end
