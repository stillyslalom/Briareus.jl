# This file holds the edits of the parts. It covers adding, removing and renaming them, their
# fields and flags, and the breaks of the section list of a hohlraum part.

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
    # The parts count too, since the target may still be that of the previous design,
    # being built.
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

Remove part `i` ([`remove_parts!`](@ref) of its name).
"""
remove_part!(ts::TargetSetup, i::Integer) = remove_parts!(ts, [ts.parts[][_check_part_index(ts, i)].name])

"""
    remove_parts!(ts::TargetSetup, names) -> Bool

Remove the parts, assemblies and user frames `names` in one edit. An assembly goes with
its subtree: every part and assembly in it, at any depth ([`assembly_members`](@ref)).
What stays keeps its place in the chamber: a part or frame that was placed in a removed
frame (a capsule in its hohlraum's, a part in a user frame or an assembly's) is placed
in the nearest frame above it that stays, where it is; a placement step that refers to
a removed frame is dropped, the part keeping its place when its other steps are
relations. [`replaced_parts`](@ref) lists them beforehand. Beams aimed in a removed
frame keep their aims, held in `:tcc`. Removing every part leaves an empty target.
"""
function remove_parts!(ts::TargetSetup, names)
    try
        return _try_build!(ts; parts = _without_parts(ts, names).parts)
    catch err
        err isa ArgumentError || rethrow()
        ts.error[] = _errmsg(err)
        return false
    end
end

"""
    replaced_parts(ts::TargetSetup, names) -> Vector{Symbol}

The parts and frames that [`remove_parts!`](@ref) of `names` would place in another
frame or take a placement step from: those that stay and are placed in a removed frame
or refer to one.
"""
replaced_parts(ts::TargetSetup, names) = _without_parts(ts, names).moved

# The design's parts without `names` (assemblies with their subtrees): `(; parts, gone,
# moved)`, `moved` the names of the parts that stay and were placed in a removed frame or
# referred to one, each kept where it is (`remove_parts!`).
function _without_parts(ts::TargetSetup, names)
    all = ts.parts[]
    gone = Set{Symbol}()
    for n in names
        i = findfirst(p -> p.name === n, all)
        i === nothing && throw(ArgumentError("no part named $n"))
        push!(gone, n)
        all[i].kind === :assembly && union!(gone, assembly_members(all, n))
    end
    owner(f) = Symbol(first(split(String(f), '.'; limit = 2)))      # a piece frame's part
    removed(f) = f !== :tcc && owner(f) in gone
    byname = Dict{Symbol,PartSpec}(p.name => p for p in all)
    # the nearest frame above `f` that stays
    function kept_frame(f)
        while removed(f)
            f = byname[owner(f)].frame
        end
        return f
    end
    fr = part_frames(ts)
    out, moved = PartSpec[], Symbol[]
    for p in all
        p.name in gone && continue
        refs = any(st -> any(removed, Briareus._step_refs(st)), p.steps)
        if !removed(p.frame) && !refs
            push!(out, p.assembly in gone ? _with(p; assembly = TOP_LEVEL) : p)
            continue
        end
        push!(moved, p.name)
        steps = PlacementStep[st for st in p.steps if !any(removed, Briareus._step_refs(st))]
        frame = kept_frame(p.frame)
        # where it is now, as its placement fields in the frame it is placed in from here
        # on; relative steps that stay apply on top of the fields, which are then kept
        vals = copy(p.values)
        if !is_port_frame(p) && !any(Briareus._relative, steps)
            F, P1 = fr[p.name], fr[frame]
            vals[:x], vals[:y], vals[:z] = to_frame(P1, F.origin)
            vals[:theta], vals[:phi], vals[:rotation] = _euler_angles(P1.R' * F.R)
        elseif frame !== p.frame
            vals = _kept_in(p, frame, fr).values
        end
        q = _with(p; frame, values = vals, steps)
        push!(out, q.assembly in gone ? _with(q; assembly = TOP_LEVEL) : q)
    end
    return (; parts = out, gone, moved)
end

# A name not taken by a component, part or frame: `base_1`, `base_2`, … (or `base`
# itself first, unless `numbered`).
function _free_name(ts::TargetSetup, base::Symbol; numbered::Bool = true)
    taken = union(Set(component_names(ts.target[])), (p.name for p in ts.parts[]), keys(part_frames(ts)))
    !numbered && !(base in taken) && return base
    k = 1
    while Symbol(base, :_, k) in taken
        k += 1
    end
    return Symbol(base, :_, k)
end

"""
    rename_map(ts::TargetSetup, p::PartSpec, new::Symbol) -> Dict{Symbol,Symbol}

Every name that changes when part `p` is renamed `new`, old name → new: the part's (and
its frame's), its piece frames' (`<part>.leh_top`, …), and its components', in the
calculations or left out (`<part>_wall`, a CAD part's `<part>_<object piece>`: the part's
name is their prefix). A rename commits it as the edit's `renames`, and what is keyed by
these names (left-out components, display modes, cutaway components, roles in the
checks, films, the sweep's subject) follows it.
"""
function rename_map(ts::TargetSetup, p::PartSpec, new::Symbol)
    old = p.name
    rn = Dict{Symbol,Symbol}(old => new)
    for ((a, _), (b, _)) in zip(piece_frame_defs(p), piece_frame_defs(_with(p; name = new)))
        rn[a] = b
    end
    pre = String(old) * "_"
    for t in (ts.target[], ts.ghost[])
        t === nothing && continue
        for c in t.components
            s = String(c.name)
            (c.part === old || c.name === old) && startswith(s, pre) &&
                (rn[c.name] = Symbol(new, s[ncodeunits(String(old))+1:end]))
        end
    end
    return rn
end

"""
    edit_part!(ts::TargetSetup, i, key, text) -> Bool

Set a field of part `i` from text: a [`PART_FIELDS`](@ref) or size key, or `:name`.
"""
function edit_part!(ts::TargetSetup, i::Integer, key::Symbol, text)
    _check_part_index(ts, i)
    p = ts.parts[][i]
    try
        _no_geometry(p) && key in (:albedo, :conversion, :h) && throw(ArgumentError("a $(p.kind) has no field :$key"))
        is_port_frame(p) && key in (:theta, :phi) &&
            throw(ArgumentError("a port frame's axis is its port's: change its roll, or place another frame in it"))
        if key === :name
            s = strip(String(text))
            occursin(r"^[A-Za-z_][A-Za-z0-9_]*$", s) ||
                throw(ArgumentError("a part name is letters, digits and _, starting with a letter"))
            Symbol(s) == p.name && return true
            new = _with(p; name = Symbol(s))
            renames = rename_map(ts, p, new.name)
            # parts placed in this one's frame or its pieces', or with steps that refer to
            # them, follow the rename, and so do the beams aimed in them; the components
            # left out stay out (the session's state follows `renames` too)
            ps = [foldl((q, (a, b)) -> _renamed(q, a, b), renames; init = q) for q in ts.parts[]]
            ps[i] = new
            excluded = Set{Symbol}(get(renames, n, n) for n in ts.excluded[])
            return _try_build!(ts; parts = ps, excluded, renames)
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
            new = key === :profile ? _switched_profile(p, vals) : _with(p; values = vals)
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

# Hohlraum part `p` with values `vals`, whose `profile` changed: a frustraum without an end
# radius gets one (60% of the waist radius, at least 1.2 times the LEH radius, to the
# micrometer), and a section list without breaks starts from the current profile's
# (`_default_sections`: a cylinder in four equal sections).
function _switched_profile(p, vals)
    p.kind === :hohlraum || return _with(p; values = vals)
    _profile_kind(vals)                         # a valid profile
    if vals[:profile] == 1 && vals[:end_radius] == 0
        vals[:end_radius] = round(max(0.6 * vals[:radius], 1.2 * vals[:leh_radius]); digits = 3)
    end
    sections = vals[:profile] == 2 && isempty(p.sections) ? _default_sections(p) : p.sections
    return _with(p; values = vals, sections)
end

# Commit hohlraum part `i` with the section breaks `s` if the geometry accepts them.
function _set_sections!(ts::TargetSetup, i::Integer, s)
    ps = copy(ts.parts[])
    ps[i] = _with(ps[i]; sections = s)
    return _try_build!(ts; parts = ps)
end

function _section_part(ts::TargetSetup, i::Integer)
    _check_part_index(ts, i)
    p = ts.parts[][i]
    p.kind === :hohlraum && _profile_kind(p.values) === :sections ||
        throw(ArgumentError("$(p.name) has no section list: set its profile to sections"))
    return p
end

"""
    edit_hohlraum_break!(ts::TargetSetup, i, k, field, text) -> Bool

Set the height (`field = :z`, mm) or the wall radius (`:r`, mm) of break `k` of hohlraum
part `i`'s section list ([`hohlraum_profile`](@ref)) from text-field text. The first and
last breaks' heights follow the part's length (edit the length instead). A value the
geometry rejects (breaks out of order, a radius not positive, an entrance hole wider
than its end) sets `error` and returns `false`, the target unchanged.
"""
function edit_hohlraum_break!(ts::TargetSetup, i::Integer, k::Integer, field::Symbol, text)
    try
        p = _section_part(ts, i)
        s = copy(p.sections)
        1 <= k <= length(s) || throw(ArgumentError("no break $k"))
        field in (:z, :r) || throw(ArgumentError("a break has a height z and a radius r"))
        v = _parse_number(text, String(field))
        field === :z && k in (1, length(s)) &&
            throw(ArgumentError("the end breaks are at ±length/2: edit the length"))
        s[k] = field === :z ? (v, s[k][2]) : (s[k][1], v)
        return _set_sections!(ts, i, s)
    catch err
        err isa ArgumentError || rethrow()
        ts.error[] = _errmsg(err)
        return false
    end
end

"""
    add_hohlraum_break!(ts::TargetSetup, i) -> Bool

Split the longest section of hohlraum part `i`'s section list at its middle (a break on
the wall, so the shape is unchanged and the section becomes two pieces).
"""
function add_hohlraum_break!(ts::TargetSetup, i::Integer)
    try
        p = _section_part(ts, i)
        prof = hohlraum_profile(p)
        z, r = prof.z, prof.r
        k = argmax([hypot(z[j+1] - z[j], r[j+1] - r[j]) for j in 1:length(z)-1])
        s = copy(p.sections)
        insert!(s, k + 1, ((z[k] + z[k+1]) / 2, (r[k] + r[k+1]) / 2))
        return _set_sections!(ts, i, s)
    catch err
        err isa ArgumentError || rethrow()
        ts.error[] = _errmsg(err)
        return false
    end
end

"""
    remove_hohlraum_break!(ts::TargetSetup, i, k) -> Bool

Remove break `k` of hohlraum part `i`'s section list, joining the two sections beside it;
the end breaks stay.
"""
function remove_hohlraum_break!(ts::TargetSetup, i::Integer, k::Integer)
    try
        p = _section_part(ts, i)
        1 < k < length(p.sections) || throw(ArgumentError("the end breaks stay: remove an inner break"))
        s = copy(p.sections)
        deleteat!(s, k)
        return _set_sections!(ts, i, s)
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
    ts.parts[] = ps                     # the geometry is unchanged, so no rebuild
    return true
end

"""
    set_part_flag!(ts::TargetSetup, i, flag, on::Bool) -> Bool

Set `:two_sided`, `:inward` or `:aperture` of part `i`. A rectangle or disk observed as a
flat face stays one-sided ([`observed_flat`](@ref)): making it two-sided sets `error` and
returns `false`.
"""
function set_part_flag!(ts::TargetSetup, i::Integer, flag::Symbol, on::Bool)
    _check_part_index(ts, i)
    flag in (:two_sided, :inward, :aperture) || throw(ArgumentError("unknown part flag :$flag"))
    p = ts.parts[][i]
    getfield(p, flag) == on && return true
    if flag === :two_sided && on && p.observed && _part_shape(p) in _FLAT_OBSERVED
        ts.error[] = "the observed face is one-sided: give the role to another part first"
        return false
    end
    ps = copy(ts.parts[])
    ps[i] = _with(p; (flag => on,)...)
    return _try_build!(ts; parts = ps)
end

"""
    part_index(ts::TargetSetup, name::Symbol) -> Union{Nothing,Int}
"""
part_index(ts::TargetSetup, name::Symbol) = findfirst(p -> p.name === name, ts.parts[])
