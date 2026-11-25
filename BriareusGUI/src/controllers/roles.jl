# This file covers the roles of parts and of the pieces of a hohlraum part (observed, foil and
# wall), the components that hold them, and the capsule fitted to a CAD part.

"""
Roles a part can take in the analyses, held by the part's flags (its name does not carry
a role):
`:none`; `:observed`, the surface whose x-ray drive is analyzed (flag `observed`; one
part at a time, see [`observed_name`](@ref)): a sphere (a capsule: its drive's modes,
the θ–φ map and the capsule keep-out) or a flat face (a rectangle or disk: a package's
face, an opacity sample; the face map of its incident x-ray flux and that flux's
statistics over the face, see [`observed_flat`](@ref)); `:foil` (a backlighter foil, a flat part or a mesh: any number
of parts can be foils, and the active one ([`active_foil`](@ref)) is the foil whose spot
statistics, map, objective, line of sight and captured power the analyses take); and
`:wall` (a re-emitting cavity surface: the drive history's wall temperature, and a
surface the clearance and TaLIS checks land the beams on; any number of parts, see
[`wall_names`](@ref)). The pieces of a hohlraum part hold separate roles
([`PIECE_ROLES`](@ref), [`set_piece_role!`](@ref)).
"""
const PART_ROLES = (:none, :observed, :foil, :wall)

const PART_ROLE_LABELS = Dict(:none => "None", :observed => "Observed (x-ray drive: a capsule or a face)",
                              :foil => "Backlighter foil (spot uniformity)",
                              :wall => "Wall (re-emitting cavity surface)")

"Role of a part (one of [`PART_ROLES`](@ref)), from its flags."
part_role(p::PartSpec) = p.observed ? :observed : p.foil ? :foil : p.wall ? :wall : :none

"""
Roles a piece of a hohlraum part can hold ([`set_piece_role!`](@ref)): `:none`,
`:observed` (a cap, flat: the observed surface's face map and flux statistics, as a
package's face) and `:wall` (a re-emitting cavity wall: the drive history's and its
cavity's wall temperature). The foil role stays on parts.
"""
const PIECE_ROLES = (:none, :observed, :wall)

const PIECE_ROLE_LABELS = Dict(:none => "None", :observed => "Observed (x-ray drive: a flat face)",
                               :wall => "Wall (re-emitting cavity surface)")

"""
    piece_role(p::PartSpec, piece::Symbol) -> Symbol

The role (one of [`PIECE_ROLES`](@ref)) piece `piece` (without the prefix, one of
[`hohlraum_pieces`](@ref)) of hohlraum part `p` holds: the one its row in the model tree
gave it (`p.piece_roles`), else its default, `:wall` for the wall (each section of it)
and the caps and `:none` for the entrance holes, their lips and the windows.
"""
piece_role(p::PartSpec, piece::Symbol) = get(p.piece_roles, piece, _default_piece_role(p, piece))

_default_piece_role(p::PartSpec, piece::Symbol) =
    piece in (:cap_top, :cap_bottom) || piece in _wall_pieces(p) ? :wall : :none

# The pieces of hohlraum part `p` that are flat (the caps, which can be observed).
const _FLAT_PIECES = (:cap_top, :cap_bottom)

# The piece of hohlraum part `p` that holds the observed role (`nothing`: none). A role
# given to a piece the part no longer has (a cap of an end since opened) is held by
# nothing until the piece is back.
function _observed_piece(p::PartSpec)
    p.kind === :hohlraum || return nothing
    ps = hohlraum_pieces(p)
    for s in ps
        piece_role(p, s) === :observed && return s
    end
    return nothing
end

"""
    piece_roles(ts::TargetSetup, name::Symbol) -> Union{Nothing,NamedTuple}

The role of component `name` when it is a piece of a hohlraum part (`<part>_<piece>`):
`(; part, piece, role, roles)`, `role` the piece's ([`piece_role`](@ref)) and `roles`
those it can take ([`PIECE_ROLES`](@ref): an entrance hole, its lip or an open window is
an aperture and takes none; only a cap, flat, can be observed). `nothing` for any other
component.
"""
function piece_roles(ts::TargetSetup, name::Symbol)
    w = _hohlraum_piece(ts, name)
    w === nothing && return nothing
    p, s = w
    return (; part = p.name, piece = s, role = piece_role(p, s), roles = _piece_roles_allowed(p, s))
end

# The hohlraum part and piece (without the prefix) of component `name`, or `nothing`.
function _hohlraum_piece(ts::TargetSetup, name::Symbol)
    n = String(name)
    for p in ts.parts[]
        p.kind === :hohlraum || continue
        pre = String(p.name) * "_"
        startswith(n, pre) || continue
        s = Symbol(n[ncodeunits(pre)+1:end])
        s in hohlraum_pieces(p) && return (p, s)
    end
    return nothing
end

# The roles piece `s` of hohlraum part `p` can take.
function _piece_roles_allowed(p::PartSpec, s::Symbol)
    _piece_aperture(p, s) && return (:none,)
    return s in _FLAT_PIECES ? PIECE_ROLES : (:none, :wall)
end

# Whether piece `s` of hohlraum part `p` is an aperture: an entrance hole, its lip, an
# open window.
function _piece_aperture(p::PartSpec, s::Symbol)
    str = String(s)
    startswith(str, "leh_") && return true
    k = findfirst(i -> window_name(i) === s, eachindex(p.windows))
    return k !== nothing && p.windows[k].open
end

"""
    set_piece_role!(ts::TargetSetup, name::Symbol, role::Symbol) -> Bool

Give piece `name` of a hohlraum part (its component, `<part>_<piece>`) a separate role
(one of [`PIECE_ROLES`](@ref)) in place of its default ([`piece_role`](@ref)): a cap
`:observed` becomes the observed surface ([`observed_name`](@ref): a flat face, its face
map and flux statistics) and leaves the wall role, so its cavity's wall temperature
([`cavity_walls`](@ref)) no longer counts it; a wall section or a cap given `:none`
leaves the walls too; a closed window given `:wall` joins them. One observed surface at
a time: a part or piece that held it gives it up. The components (and so the target)
are unchanged. Saved with the session, and undone like any edit. A role the piece cannot
take (an aperture takes none; only a cap can be observed), or a component that is not a
hohlraum part's piece, sets `error` and returns `false`.
"""
function set_piece_role!(ts::TargetSetup, name::Symbol, role::Symbol)
    role in PIECE_ROLES || throw(ArgumentError("a piece's role must be one of $PIECE_ROLES"))
    w = _hohlraum_piece(ts, name)
    if w === nothing
        ts.error[] = "$name is not a piece of a hohlraum part: a part's role is set on the part"
        return false
    end
    p, s = w
    piece_role(p, s) === role && (ts.error[] = ""; return true)
    if !(role in _piece_roles_allowed(p, s))
        ts.error[] = _piece_aperture(p, s) ? "$name is an aperture (transparent): it takes no role" :
                     "only a cap, flat, can be the observed surface: $name is curved"
        return false
    end
    ps = copy(ts.parts[])
    role === :observed && _drop_observed!(ps)            # one observed surface at a time
    i = findfirst(q -> q.name === p.name, ps)
    q = ps[i]
    roles = copy(q.piece_roles)
    role === _default_piece_role(q, s) ? delete!(roles, s) : (roles[s] = role)
    ps[i] = _with(q; piece_roles = roles)
    ts.error[] = ""
    ts.parts[] = ps                                      # the geometry is unchanged
    return true
end

"""
    observed_flat(ts::TargetSetup) -> Bool

Whether the observed surface ([`observed_name`](@ref)) is a flat face (a rectangle or a
disk, a VisRad rectangle, a hohlraum part's cap) instead of a sphere: the drive then
maps its face and gives its flux's statistics over the face, and the capsule-only views
(modes, keep-out, the θ–φ map) have no capsule. Decided on the built target's elements of
the observed surface, by the rule the drive applies to the elements its scenario
observes (every element facing one way), so the two agree; `false` when the observed
surface is not in the calculations.
"""
function observed_flat(ts::TargetSetup)
    o = observed_name(ts)
    o === nothing && return false
    t = ts.target[]
    o in component_names(t) || return false
    return _flat_face(t, elements(t, o))
end

# The shapes a part given the observed role takes it as a flat face (`set_part_role!`).
const _FLAT_OBSERVED = (:rectangle, :disk)

"""
    observed_name(ts::TargetSetup) -> Union{Symbol,Nothing}

The surface whose x-ray drive the analyses observe (the drive's modes and θ–φ map, the
capsule keep-out and summaries, the drive history, the retune, the results export): the
part with the observed role (`PartSpec.observed`, given by [`set_part_role!`](@ref); a
starter's capsule holds it), or the piece of a hohlraum part given it by
[`set_piece_role!`](@ref) (its component name, `<part>_cap_top`, say), while it is in the
calculations, else `nothing`. One surface at a time.
"""
function observed_name(ts::TargetSetup)
    ex = ts.excluded[]
    for p in ts.parts[]
        p.name in ex && continue
        p.observed && return p.name
        s = _observed_piece(p)
        s === nothing && continue
        n = piece_name(p.name, s)
        n in ex || return n
    end
    return nothing
end

"""
    wall_names(ts::TargetSetup) -> Vector{Symbol}

The re-emitting cavity walls in the calculations (the drive history's wall temperature;
surfaces the clearance and TaLIS checks land the beams on), in the parts' order: each
part with the role `:wall`, and each hohlraum part's pieces holding it ([`piece_role`](@ref):
by default its wall and caps, `<part>_wall`, `<part>_cap_top`, `<part>_cap_bottom`, an
open end having no cap; a piece given another role in the model tree leaves them; the
checks keep the caps as hardware, the entrance holes' edges).
"""
function wall_names(ts::TargetSetup)
    ex = ts.excluded[]
    out = Symbol[]
    for p in ts.parts[]
        p.name in ex && continue
        p.wall && push!(out, p.name)
        p.kind === :hohlraum && append!(out, _hohlraum_walls(p))
    end
    return out
end

# The pieces of a hohlraum part that hold the wall role (`piece_role`; by default its wall,
# each section of it, and the caps it has), in the order of its pieces.
_hohlraum_walls(p) = [piece_name(p.name, s) for s in hohlraum_pieces(p) if piece_role(p, s) === :wall]

"""
    default_cut(ts::TargetSetup) -> Set{Symbol}

The components the cutaway cuts by default: the shell of each hohlraum part, its wall
(each section of a profile's) and the caps of its ends that are not open, whatever role
a piece holds; the entrance holes, lips and windows are not cut. A session's cutaway
starts as this and follows it through edits while the cut components are this default.
"""
function default_cut(ts::TargetSetup)
    out = Set{Symbol}()
    for p in ts.parts[]
        p.kind === :hohlraum || continue
        for s in hohlraum_pieces(p)
            (s in _wall_pieces(p) || s in (:cap_top, :cap_bottom)) && push!(out, piece_name(p.name, s))
        end
    end
    return out
end

"""
    cavity_walls(ts::TargetSetup) -> Vector{Pair{Symbol,Vector{Symbol}}}

The cavities of the design for the drive's per-cavity wall temperature: one per
hohlraum part in the calculations, in the parts' order, `part => its wall pieces` (the
pieces holding the wall role, [`piece_role`](@ref): by default its wall, each section of
it, and the caps it has, so a cap made the observed surface is not among them;
[`wall_names`](@ref) lists them all together). Parts that hold the wall role without being a hohlraum part (a shield, a
VisRad hohlraum) count in the drive history's wall temperature but do not form a
separate cavity. A cavity's walls are named piece by piece: every piece of a hohlraum part,
its entrance holes included, carries the part's name as `Component.part`.
"""
cavity_walls(ts::TargetSetup) =
    Pair{Symbol,Vector{Symbol}}[p.name => _hohlraum_walls(p) for p in ts.parts[]
                                if p.kind === :hohlraum && !(p.name in ts.excluded[])]

"""
    observed_cavity(ts::TargetSetup) -> Union{Symbol,Nothing}

The hohlraum part ([`cavity_walls`](@ref)) that holds the observed surface
([`observed_name`](@ref)): the first whose closed volume (its profile about its axis, ends
included, to 10⁻⁶ of its length) contains the observed part's area-weighted centroid in
the built target. A capsule at a hohlraum's center is in it; a shock-tube package's face
closing a halfraum's open end lies in that end's plane, so it is in that halfraum, and so
does an observed cap. `nothing` without an observed surface or when no cavity holds it.
"""
function observed_cavity(ts::TargetSetup)
    o = observed_name(ts)
    t = ts.target[]
    (o === nothing || t === nothing || !_has_part(t, o)) && return nothing
    c = Briareus.part_centroid(t, o)
    frames = part_frames(ts)
    for p in ts.parts[]
        (p.kind === :hohlraum && !(p.name in ts.excluded[]) && haskey(frames, p.name)) || continue
        q = to_frame(frames[p.name], c)
        L = _hohlraum_length(p)
        tol = 1e-6 * L
        abs(q[3]) <= L / 2 + tol && hypot(q[1], q[2]) <= _wall_radius_at(p, clamp(q[3], -L / 2, L / 2)) + tol &&
            return p.name
    end
    return nothing
end

"""
    foil_names(ts::TargetSetup) -> Vector{Symbol}

The backlighter foils of the design in the calculations: each part with the role
`:foil`, in the parts' order. Foils left out of the calculations are not listed.
"""
foil_names(ts::TargetSetup) = Symbol[p.name for p in ts.parts[] if p.foil && !(p.name in ts.excluded[])]

"""
    active_foil(ts::TargetSetup) -> Union{Symbol,Nothing}

The foil the foil analyses take (the drive's spot statistics and foil map, the foil
objective, the sweep's line of sight and captured power): the one chosen with
[`set_active_foil!`](@ref) while it is still a foil in the calculations, else the first
of [`foil_names`](@ref), else `nothing`. A part given the foil role when the design has
no foil becomes the chosen one.
"""
function active_foil(ts::TargetSetup)
    fs = foil_names(ts)
    isempty(fs) && return nothing
    a = ts.active_foil[]
    return a in fs ? a : first(fs)
end

"""
    set_active_foil!(ts::TargetSetup, name::Symbol) -> Bool

Make foil `name` (one of [`foil_names`](@ref)) the one the foil analyses take
([`active_foil`](@ref)). A name that is not a foil in the calculations sets `error` and
returns `false`. Saved with the session, and undone like any edit.
"""
function set_active_foil!(ts::TargetSetup, name::Symbol)
    if !(name in foil_names(ts))
        ts.error[] = "$name is not a backlighter foil in the calculations"
        return false
    end
    ts.error[] = ""
    ts.active_foil[] === name || (ts.active_foil[] = name)
    return true
end

"""
    foil_components(ts::TargetSetup, t = ts.target[]) -> Vector{Symbol}

The components of `t` that belong to a backlighter foil, active or not: the components
of each part with the foil role (a VisRad part's carry its name as a prefix). The
clearance and TaLIS checks take them as surfaces the beams land on.
"""
foil_components(ts::TargetSetup, t::Target = ts.target[]) = _held_components(ts, t, p -> p.foil, ())

"""
    observed_components(ts::TargetSetup, t = ts.target[]) -> Vector{Symbol}
    wall_components(ts::TargetSetup, t = ts.target[]) -> Vector{Symbol}

The components of `t` that hold the observed role ([`observed_name`](@ref)) or the wall
role ([`wall_names`](@ref)): the components of the parts with the role (a VisRad part's
carry its name as a prefix) and, for the wall, each hohlraum part's wall and caps.
"""
function observed_components(ts::TargetSetup, t::Target = ts.target[])
    o = observed_name(ts)
    o === nothing && return Symbol[]
    part_index(ts, o) === nothing && return Symbol[c.name for c in t.components if c.name === o]   # a piece
    return _held_components(ts, t, p -> p.name === o, ())
end

wall_components(ts::TargetSetup, t::Target = ts.target[]) =
    _held_components(ts, t, p -> p.wall,
                     [n for p in ts.parts[] if p.kind === :hohlraum for n in _hohlraum_walls(p)])

# The components of `t` of the parts `holds(p)` selects (a VisRad part's carry its name as
# a prefix), and the components named `named` (a hohlraum part's walls).
function _held_components(ts::TargetSetup, t::Target, holds, named)
    fs = Set{Symbol}(p.name for p in ts.parts[] if holds(p))
    vr = Set{Symbol}(p.name for p in ts.parts[] if holds(p) && p.kind === :visrad)
    out = Symbol[]
    for c in t.components
        n = c.name
        if n in named || n in fs || c.part in fs || any(v -> startswith(String(n), String(v)), vr)
            push!(out, n)
        end
    end
    return out
end

# The surfaces the clearance check lands the beams on by default (`lands_on`: walls, the
# observed surface, the foils) and those the TaLIS check leaves out (`observed`), as part
# and component names of `t` (the checks list parts: the parts of the wall and observed
# components join them; the foils by component).
function _landing(ts::TargetSetup, t::Target = ts.target[])
    obs = _with_parts(t, observed_components(ts, t))
    # The caps of a hohlraum part hold the wall role for the drive, but the checks keep
    # them as hardware, since they are the edges of the entrance holes, as the caps of
    # the former template were.
    walls = filter(n -> !occursin(r"_cap_(top|bottom)$", String(n)) || !_piece_of_hohlraum(ts, n),
                   wall_components(ts, t))
    capture = unique!([_with_parts(t, walls); foil_components(ts, t)])
    return (; lands_on = unique!([capture; obs]), capture, observed = obs)
end

# Components `ns` of `t` and their parts; a piece of a part in pieces is listed by itself.
_with_parts(t::Target, ns) = unique!([ns; [c.part for c in t.components if c.name in ns && !_is_piece(c)]])

# Whether component `name` is a piece of one of the design's hohlraum parts.
_piece_of_hohlraum(ts::TargetSetup, name::Symbol) =
    any(p -> p.kind === :hohlraum && startswith(String(name), String(p.name) * "_"), ts.parts[])

# The walls' and the observed surface's components and parts, as sets.
_held_sets(ts::TargetSetup, t::Target = ts.target[]) =
    (; walls = Set(setdiff(_landing(ts, t).capture, foil_components(ts, t))),
     observed = Set(_with_parts(t, observed_components(ts, t))))

# The object kind a part stands for (the kind of a VisRad object).
function _part_shape(p::PartSpec)
    p.kind === :visrad || return p.kind
    k = _visrad_object(p.file, p.object).kind
    return k === :mesh ? :mesh : k
end

"""
    set_part_role!(ts::TargetSetup, i, role::Symbol; done = () -> nothing) -> Bool

Give part `i` a role in the analyses (one of [`PART_ROLES`](@ref)); the part keeps its name,
and its components (and so the target) are unchanged. `:foil` (a disk, rectangle or mesh, or
a VisRad rectangle or CAD part: a backlighter whose spot the foil objective makes uniform)
marks the part as a foil; the part becomes the active foil ([`active_foil`](@ref)) when the
design had none. `:observed` makes the part the observed surface ([`observed_name`](@ref));
another part or piece that held it gives it up. A hohlraum part does not take a role; its
pieces hold theirs (the wall and caps the wall role by default, see
[`set_piece_role!`](@ref)). The observed surface is a sphere (a capsule: the mode analysis
and the keep-out take its center and radius) or a flat face (a one-sided rectangle or disk,
or a VisRad rectangle: [`observed_flat`](@ref)): such a part takes the role itself; a CAD
part (a capsule with its fill tube, say) is replaced in the calculations by a sphere part
fitted to it ([`fit_sphere`](@ref): named `capsule`, or `capsule_2`, … when that is taken;
its element size the target's capsule size, its albedo and conversion the part's), which
takes the role, and left out of the calculations itself. `:wall` (any surface but an
aperture) marks a re-emitting cavity wall ([`wall_names`](@ref)). `:none` clears the part's
role. A role the part cannot take sets `error` and returns `false`. In a window (`spawn[]`)
a large CAD capsule is fitted on a worker (see [`TargetSetup`](@ref)): `true` then means the
fit is under way, and `done()` runs once the role is given (at once otherwise).
"""
function set_part_role!(ts::TargetSetup, i::Integer, role::Symbol; done = () -> nothing)
    _check_part_index(ts, i)
    role in PART_ROLES || throw(ArgumentError("a part's role must be one of $PART_ROLES"))
    ps = copy(ts.parts[])
    p = ps[i]
    part_role(p) === role && (done(); return true)
    try
        _no_geometry(p) && throw(ArgumentError("a $(p.kind) has no surface to take a role"))
        p.kind === :hohlraum &&
            throw(ArgumentError("a hohlraum part's pieces hold the roles (its wall and caps the wall role): " *
                                "select a piece to change its role"))
        shape = _part_shape(p)
        if role === :observed && shape !== :sphere && !(shape in _FLAT_OBSERVED)
            shape === :mesh ||
                throw(ArgumentError("the observed surface is a sphere, a flat face (a rectangle or disk), " *
                                    "or a CAD part a sphere is fitted to"))
            return _capsule_from_part!(ts, i; done)
        end
        # (a two-sided face faces both ways: not a flat face, `observed_flat`)
        role === :observed && shape in _FLAT_OBSERVED && p.two_sided &&
            throw(ArgumentError("the observed face is one-sided: make the part one-sided first"))
        role === :foil && !(shape in (:rectangle, :disk, :mesh)) &&
            throw(ArgumentError("the backlighter foil is a flat part (a rectangle or disk) or a mesh"))
        role === :wall && p.aperture &&
            throw(ArgumentError("an aperture is transparent: a wall is a surface"))
    catch err
        err isa ArgumentError || rethrow()
        ts.error[] = _errmsg(err)
        return false
    end
    role === :observed && _drop_observed!(ps)            # one observed surface at a time
    ps[i] = _with(ps[i]; foil = role === :foil, observed = role === :observed, wall = role === :wall)
    first_foil = role === :foil && active_foil(ts) === nothing
    # The geometry is unchanged, so the parts are committed without building the target.
    ts.error[] = ""
    ts.parts[] = ps
    first_foil && (ts.active_foil[] = ps[i].name)
    done()
    return true
end

# Clear the observed role of every part of `ps` and of their pieces (in place; a piece
# returns to its default role).
function _drop_observed!(ps::Vector{PartSpec})
    for k in eachindex(ps)
        ps[k].observed && (ps[k] = _with(ps[k]; observed = false))
        r = ps[k].piece_roles
        any(==(:observed), values(r)) &&
            (ps[k] = _with(ps[k]; piece_roles = filter(kv -> last(kv) !== :observed, r)))
    end
    return ps
end

"""
    fit_sphere(points, weights = ones; inlier = 0.1) -> (center, radius)

A sphere through most of `points` (mm), weighted by `weights` (e.g. the surface area
each point stands for, so a finely tessellated fill tube does not outweigh the capsule
it is attached to): a sphere fitted to every point (weighted least squares, algebraic),
then refitted to the points within 0.3, then four times within `inlier`, × its radius.
Robust to a minority of the weight off the sphere. A shell modeled with both faces
fits between them; the radius returned is the outer face's (the weighted 75th
percentile of the distances within 20% of the fit).
"""
function fit_sphere(points::AbstractVector, weights::AbstractVector = ones(length(points));
                    inlier::Real = 0.1)
    P = [Vec3(p...) for p in points]
    w = Float64.(weights)
    length(P) >= 4 && length(w) == length(P) && all(>=(0), w) && sum(w) > 0 ||
        throw(ArgumentError("a sphere needs at least four points with nonnegative weights"))
    function wmedian(x, w)
        o = sortperm(x)
        cw = cumsum(w[o])
        return x[o[searchsortedfirst(cw, cw[end] / 2)]]
    end
    c = sum(w .* P) / sum(w)
    r = wmedian([norm(p - c) for p in P], w)
    for it in 1:6
        # every point first, then those near the sphere
        band = (it == 1 ? Inf : it == 2 ? 0.3 : Float64(inlier)) * r
        sel = [k for k in eachindex(P) if abs(norm(P[k] - c) - r) <= band && w[k] > 0]
        length(sel) >= 4 || break
        # |p|² = 2 c·p + (r² − |c|²): linear in (c, k); rows weighted by √w
        sw = sqrt.(w[sel])
        A = reduce(vcat, [sw[j] * [2P[k][1] 2P[k][2] 2P[k][3] 1.0] for (j, k) in enumerate(sel)])
        b = [sw[j] * dot(P[k], P[k]) for (j, k) in enumerate(sel)]
        x = A \ b
        c = Vec3(x[1], x[2], x[3])
        r = sqrt(max(x[4] + dot(c, c), 0.0))
    end
    # a shell's inner and outer faces fit in between: take the outer one (the radiating
    # surface), the weighted 75th percentile of the distances near the sphere
    near = [k for k in eachindex(P) if abs(norm(P[k] - c) - r) <= 0.2r]
    if length(near) >= 4
        d = [norm(P[k] - c) for k in near]
        o = sortperm(d)
        cw = cumsum(w[near][o])
        r = d[o[searchsortedfirst(cw, 0.75 * cw[end])]]
    end
    return c, r
end

# Role `:observed` for CAD part `i`: a sphere part with the role, fitted to it, in its place in
# the calculations. With `spawn[]` a part of more than `BACKGROUND_TRIANGLES` facets is
# fitted on a worker (`fitting`), and the role is given when the fit is done, to the part
# as it is then (fitted again if it moved meanwhile; dropped if it went away).
function _capsule_from_part!(ts::TargetSetup, i::Integer; done = () -> nothing)
    p = ts.parts[][i]
    cs = part_components(p, ts.mesh[].h; frame = part_frames(ts)[p.name], raw = true)
    if ts.spawn[] && sum(c -> length(c.mesh.triangles), cs; init = 0) > BACKGROUND_TRIANGLES[]
        job = ts.fit_jobs[p.name] = get(ts.fit_jobs, p.name, 0) + 1
        _assign!(ts.fitting, union(ts.fitting[], (p.name,)))
        ts.error[] = ""
        deliver = ts.deliver[]
        _run_job(() -> (out = _try_job(() -> _fit_capsule(cs, p.name));
                        deliver(() -> _capsule_fitted!(ts, p, job, out, done))), true)
        return true
    end
    ok = _apply_capsule!(ts, i, _fit_capsule(cs, p.name)...)
    ok && done()
    return ok
end

# The sphere (center, radius) of a CAD capsule's components (raw facets in the chamber).
function _fit_capsule(cs, name)
    # each triangle's vertices weighted by a third of its area
    pts, wts = Vec3[], Float64[]
    for comp in cs, (a, b, d) in comp.mesh.triangles
        va, vb, vd = comp.mesh.vertices[a], comp.mesh.vertices[b], comp.mesh.vertices[d]
        A = norm(cross(vb - va, vd - va)) / 6
        push!(pts, va, vb, vd)
        append!(wts, (A, A, A))
    end
    c, r = fit_sphere(pts, wts)
    r > 0 || throw(ArgumentError("$name: no sphere fits it"))
    return c, r
end

# A worker fitted the capsule to part `p` (job `job`): give it the role now.
function _capsule_fitted!(ts::TargetSetup, p::PartSpec, job::Integer, out, done = () -> nothing)
    get(ts.fit_jobs, p.name, 0) == job || return false       # asked again, or canceled
    delete!(ts.fit_jobs, p.name)
    ok = try
        i = findfirst(q -> q === p, ts.parts[])
        j = part_index(ts, p.name)
        if out.err !== nothing
            ts.error[] = _errmsg(out.err)
            false
        elseif i !== nothing
            _apply_capsule!(ts, i, out.value...) && (done(); true)
        elseif j !== nothing                                  # edited meanwhile, so fit it again
            set_part_role!(ts, j, :observed; done)
        else
            ts.error[] = "$(p.name) went away before its capsule sphere was fitted"
            false
        end
    catch err
        err isa ArgumentError || rethrow()
        ts.error[] = _errmsg(err)
        false
    end
    # the role is committed (its target building, or built) before the fit stops
    # holding the drive back, unless the part is being fitted again
    haskey(ts.fit_jobs, p.name) || _assign!(ts.fitting, setdiff(ts.fitting[], (p.name,)))
    return ok
end

# Drop the capsule fits in flight, as after an undo, a loaded session or a closed window.
function _cancel_fits!(ts::TargetSetup)
    empty!(ts.fit_jobs)
    _assign!(ts.fitting, Set{Symbol}())
    return ts
end

function _apply_capsule!(ts::TargetSetup, i::Integer, c, r)
    ps = _drop_observed!(copy(ts.parts[]))
    p = ps[i]
    ex = copy(ts.excluded[])
    taken = union(Set(all_component_names(ts)), (q.name for q in ps))
    name, k = :capsule, 1
    while name in taken
        k += 1
        name = Symbol(:capsule_, k)
    end
    push!(ps, PartSpec(:sphere; name, observed = true,
                       values = Dict(:x => c[1], :y => c[2], :z => c[3], :r1 => r,
                                     :h => ts.mesh[].capsule_h, :albedo => p.values[:albedo],
                                     :conversion => p.values[:conversion])))
    push!(ex, p.name)
    return _try_build!(ts; parts = ps, excluded = ex)
end
