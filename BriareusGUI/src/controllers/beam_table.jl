# The beam table: a facility preset or an imported beam table, its beams gathered in
# groups (VisRad's clusters: the cone by default, or any named set such as
# "backlighter"), edited per group, quad, beam or any set of beams. Edits go through
# `GroupPointing` with only the edited field, so a group edit sets that field on every
# member and leaves the others as they were. NIF quads can be split into their four
# beams and merged back.

"""
Beam presets:

- `:nif_quads` — the 48 NIF indirect-drive quads (`nif_quads`, 8 TW each), pointed at
  the published cone (r, z) aims of `examples/nif_hohlraum.jl`;
- `:nif_beams` — all 192 NIF beams (`nif_beams`, 2 TW each), same cone pointing;
- `:omega_hohlraum` — the 40 OMEGA hohlraum beams (cones 1–3, `omega_beams`, 0.5 TW),
  aimed at the LEH centers of a hohlraum target (TCC otherwise);
- `:omega_backlighter` — the 60-port OMEGA set reduced to the beams that reach the
  driven face of a foil target from at least 20° off the foil axis, aimed at the foil
  center (0.4 TW, SG4 spots), as in `examples/backlighter.jl`.
"""
const BEAM_PRESETS = (:nif_quads, :nif_beams, :omega_hohlraum, :omega_backlighter)

const PRESET_LABELS = Dict(:nif_quads => "NIF quads (48)", :nif_beams => "NIF beams (192)",
                           :omega_hohlraum => "OMEGA hohlraum (40)",
                           :omega_backlighter => "OMEGA backlighter", :custom => "imported table")

"""
Facilities: `:nif`, `:omega`, `:omega_ep` and `:custom` (another facility, or none).
The facility decides the presets offered ([`FACILITY_PRESETS`](@ref)) and the checks
that apply: the TaLIS target-configuration rules are NIF's. OMEGA EP's port geometry
is not built in: its beams come from a beam table or a VisRad file.
"""
const FACILITIES = (:nif, :omega, :omega_ep, :custom)

const FACILITY_LABELS = Dict(:nif => "NIF", :omega => "OMEGA", :omega_ep => "OMEGA EP",
                             :custom => "Other (imported beams)")

"Beam presets of each facility (the first is the default)."
const FACILITY_PRESETS = Dict(:nif => (:nif_quads, :nif_beams),
                              :omega => (:omega_hohlraum, :omega_backlighter),
                              :omega_ep => (), :custom => ())

"The facility of a preset (`nothing` for `:custom`)."
function facility_of(preset::Symbol)
    for (f, ps) in FACILITY_PRESETS
        preset in ps && return f
    end
    return nothing
end

"Row groupings of the beam table: per group (`:group`, see `BeamTable`'s `cluster`), per `:quad`, per `:beam`."
const BEAM_GROUPINGS = (:group, :quad, :beam)

"Editable beam-table columns (`GroupPointing` fields): aim z, r, tangential offset, defocus (mm), power (TW)."
const BEAM_FIELDS = (:z, :r, :phi, :defocus, :power)

# Published NIF cone pointing (r toward the port azimuth, z; mm) used by
# examples/nif_hohlraum.jl for a 10.1 mm hohlraum.
const NIF_CONE_POINTING = Dict(:inner23 => (-1.09, 3.09), :inner30 => (-0.60, 3.63),
                               :outer44 => (0.0, 5.20), :outer50 => (0.0, 5.01))

"""
    BeamTable(beams = Beam{Float64}[]; preset = :custom)

The beam table: `beams` (every beam of the table, in table order), `enabled` (one flag
per beam; disabled beams are kept but deliver no power), `cluster` (each beam's group
name, as VisRad's beam clusters: the beam's cone `group` by default; a group edit,
`:group` rows and the optimizer's per-group parameters follow it, while the beam's own
`group` stays its cone for the laser–plasma models), `preset` (one of
[`BEAM_PRESETS`](@ref), or `:custom` for an imported or scripted table), `by` (row
grouping, one of [`BEAM_GROUPINGS`](@ref)), `path` (the last imported or exported
file), `error` (the last rejected edit), `warnings` (what the last VisRad import warned
about: unmatched beams, port-angle mismatches, …), `facility` (one of
[`FACILITIES`](@ref); set by a preset, or with [`set_facility!`](@ref)), `aim_ref` (per
beam, the target component its aim follows, or `:none`; see
[`set_aim_reference!`](@ref)) and `ref_point` (where that component's reference point
was when the aim last followed it).
"""
struct BeamTable
    beams::Observable{Vector{Beam{Float64}}}
    enabled::Observable{Vector{Bool}}
    cluster::Observable{Vector{Symbol}}
    preset::Observable{Symbol}
    by::Observable{Symbol}
    path::Observable{String}
    error::Observable{String}
    warnings::Observable{Vector{String}}
    facility::Observable{Symbol}
    aim_ref::Observable{Vector{Symbol}}
    ref_point::Base.RefValue{Vector{Vec3}}
end

function BeamTable(beams::AbstractVector = Beam{Float64}[]; preset::Symbol = :custom)
    bs = Beam{Float64}[_float_beam(b) for b in beams]
    return BeamTable(Observable(bs), Observable(fill(true, length(bs))),
                     Observable(Symbol[b.group for b in bs]), Observable(preset),
                     Observable(:group), Observable(""), Observable(""), Observable(String[]),
                     Observable(something(facility_of(preset), :custom)),
                     Observable(fill(:none, length(bs))), Ref(fill(zero(Vec3), length(bs))))
end

_float_beam(b::Beam{Float64}) = b
_float_beam(b::Beam) = remake(b; aim = Float64.(b.aim), power = Float64(b.power),
                              defocus = Float64(b.defocus),
                              spot = with_axes(b.spot, Float64(b.spot.a), Float64(b.spot.b)),
                              dlambda = Float64(b.dlambda))

"""
    preset_beams(preset; target_kind = :hohlraum, params = TargetSetup().params[]) -> Vector{Beam}

The beams of a [`BEAM_PRESETS`](@ref) entry, pointed for a target of `target_kind`
with dimensions `params` (see [`TARGET_FIELDS`](@ref)).
"""
function preset_beams(preset::Symbol; target_kind::Symbol = :hohlraum,
                      params::AbstractDict = Dict(k => v for fs in TARGET_FIELDS for (k, _, _, v) in fs))
    if preset === :nif_quads || preset === :nif_beams
        bs = preset === :nif_quads ? nif_quads(power = 8.0) : nif_beams(power = 2.0)
        return _cone_pointing(bs, NIF_CONE_POINTING)
    elseif preset === :omega_hohlraum
        bs = omega_beams(power = 0.5)
        target_kind === :hohlraum || return bs
        h = params[:length] / 2
        return _cone_pointing(bs, Dict(g => (0.0, h) for g in (:cone1, :cone2, :cone3)))
    elseif preset === :omega_backlighter
        c = target_kind === :foil ? foil_center(params) : Vec3(0, 0, 0)
        return [remake(b; aim = c, power = 0.4)
                for b in omega_beams(cones = :all, spot = omega_spot(:SG4))
                if b.port[1] < -cosd(70)]
    end
    throw(ArgumentError("unknown beam preset :$preset; choose from $BEAM_PRESETS"))
end

# Point each cone at (r, z) with top/bottom mirroring.
function _cone_pointing(bs, rz)
    p = GroupPointing(bs; by = :group, fields = (:r, :z))
    x = initial_parameters(p)
    for (g, k) in enumerate(p.keys)
        haskey(rz, k) || continue
        x[2g-1], x[2g] = rz[k]
    end
    return apply(p, x)
end

"""
    set_facility!(bt::BeamTable, facility; target_kind = :hohlraum, params) -> bt

Choose the facility ([`FACILITIES`](@ref)). If the table's preset belongs to another
facility, the new facility's first preset replaces it (for a target of `target_kind`
with dimensions `params`); a facility without presets keeps the beams.
"""
function set_facility!(bt::BeamTable, facility::Symbol; kwargs...)
    facility in FACILITIES || throw(ArgumentError("facility must be one of $FACILITIES"))
    _assign!(bt.facility, facility)
    ps = FACILITY_PRESETS[facility]
    if !isempty(ps) && facility_of(bt.preset[]) !== facility
        load_preset!(bt, first(ps); kwargs...)
    end
    return bt
end

"""
    set_aim_reference!(bt::BeamTable, idx, name::Symbol, point = nothing) -> bt

Aim beams `idx` at `point` (the reference point of target component `name`, such as
its centroid: a backlighter foil, an entrance hole) and keep them on it: when the
target changes, [`follow_references!`](@ref) moves their aims with the component.
Pointing edits afterwards set offsets from it, which the aims keep as they follow.
`name = :none` releases the beams (their aims stay where they are).
"""
function set_aim_reference!(bt::BeamTable, idx, name::Symbol, point = nothing)
    name === :none || point !== nothing || throw(ArgumentError("give the point to aim at"))
    refs, pts, bs = copy(bt.aim_ref[]), copy(bt.ref_point[]), copy(bt.beams[])
    for i in idx
        refs[i] = name
        name === :none && continue
        pts[i] = Vec3(point...)
        bs[i] = remake(bs[i]; aim = Vec3(point...))
    end
    bt.ref_point[] = pts
    bt.aim_ref[] = refs
    bs == bt.beams[] || (bt.beams[] = bs)
    return bt
end

"""
    follow_references!(bt::BeamTable, points::AbstractDict) -> Bool

Move every beam that follows a component ([`set_aim_reference!`](@ref)) by how far that
component's reference point moved since the aim last followed it (`points`: component
name → its point now). A beam whose component is gone is released. Returns whether
anything changed.
"""
function follow_references!(bt::BeamTable, points::AbstractDict)
    refs, pts, bs = copy(bt.aim_ref[]), copy(bt.ref_point[]), copy(bt.beams[])
    changed = false
    for i in eachindex(refs)
        r = refs[i]
        r === :none && continue
        if !haskey(points, r)
            refs[i] = :none
            changed = true
            continue
        end
        p = Vec3(points[r]...)
        δ = p - pts[i]
        iszero(δ) && continue
        bs[i] = remake(bs[i]; aim = Vec3(bs[i].aim) + δ)
        pts[i] = p
        changed = true
    end
    changed || return false
    bt.ref_point[] = pts
    refs == bt.aim_ref[] || (bt.aim_ref[] = refs)
    bs == bt.beams[] || (bt.beams[] = bs)
    return true
end

"""
    load_preset!(bt::BeamTable, preset; target_kind, params)

Replace the table with a preset's beams (see [`preset_beams`](@ref)), all enabled.
"""
function load_preset!(bt::BeamTable, preset::Symbol; kwargs...)
    set_beams!(bt, preset_beams(preset; kwargs...); preset)
    return bt
end

"""
    set_beams!(bt::BeamTable, beams; preset = :custom, enabled = all, cluster = cones)

Replace the table's beams (all enabled unless `enabled` is given, grouped by cone
unless `cluster` is given).
"""
function set_beams!(bt::BeamTable, beams::AbstractVector; preset::Symbol = :custom,
                    enabled::AbstractVector{Bool} = fill(true, length(beams)),
                    cluster::AbstractVector{Symbol} = Symbol[b.group for b in beams],
                    aim_ref::AbstractVector{Symbol} = fill(:none, length(beams)),
                    ref_point::AbstractVector = fill(zero(Vec3), length(beams)))
    length(enabled) == length(beams) || throw(DimensionMismatch("one enabled flag per beam"))
    length(cluster) == length(beams) || throw(DimensionMismatch("one group per beam"))
    length(aim_ref) == length(ref_point) == length(beams) ||
        throw(DimensionMismatch("one aim reference per beam"))
    bt.enabled.val = collect(Bool, enabled)          # observers see consistent lengths
    bt.cluster.val = collect(Symbol, cluster)
    bt.aim_ref.val = collect(Symbol, aim_ref)
    bt.ref_point[] = Vec3[Vec3(p...) for p in ref_point]
    bt.preset[] = preset
    f = facility_of(preset)
    f === nothing || _assign!(bt.facility, f)
    bt.by[] === :quad && any(b -> isempty(b.quad), beams) && (bt.by[] = :group)
    bt.beams[] = Beam{Float64}[_float_beam(b) for b in beams]
    notify(bt.enabled)
    notify(bt.cluster)
    notify(bt.aim_ref)
    bt.error[] = ""
    isempty(bt.warnings[]) || (bt.warnings[] = String[])
    return bt
end

"""
    set_grouping!(bt::BeamTable, by::Symbol)

Edit the table per cone (`:group`), `:quad`, or `:beam`.
"""
function set_grouping!(bt::BeamTable, by::Symbol)
    _check_grouping(bt.beams[], by)
    _assign!(bt.by, by)
    return bt
end

function _check_grouping(beams, by::Symbol)
    by in BEAM_GROUPINGS || throw(ArgumentError("grouping must be one of $BEAM_GROUPINGS"))
    by === :quad && any(b -> isempty(b.quad), beams) &&
        throw(ArgumentError("these beams have no quad names: group per cone or per beam"))
    return nothing
end

# Rows group by cluster with z mirrored top/bottom (one value describes a cone); per
# quad or beam, z is the aim's own coordinate.
_mirror(by::Symbol) = by === :group

# The `GroupPointing` key of a grouping: `:group` rows follow the table's clusters.
function _keyfun(bt::BeamTable, by::Symbol)
    by === :group || return by
    d = Dict(b.name => c for (b, c) in zip(bt.beams[], bt.cluster[]))
    return b -> d[b.name]
end

_pointing(bt::BeamTable, fields; by::Symbol = bt.by[]) =
    GroupPointing(bt.beams[]; by = _keyfun(bt, by), fields, mirror = _mirror(by))

"""
    beam_rows(bt::BeamTable; by = bt.by[]) -> Vector{NamedTuple}

One row per group of the grouping `by` (the table's by default): `key` (cone, quad, or beam name), `count`
(beams), `members` (beam indices), `theta`/`phi` (port angles, degrees: per cone the
mean polar angle from the nearer pole and `NaN` azimuth; per quad or beam the mean
port angles), the mean
`z`, `r`, `phi_offset`, `defocus` (mm) and `power` (TW per beam), and `enabled`
(`true`, `false`, or `missing` when only some members are enabled).
"""
function beam_rows(bt::BeamTable; by::Symbol = bt.by[])
    bs = bt.beams[]
    isempty(bs) && return NamedTuple[]
    _check_grouping(bs, by)
    p = _pointing(bt, BEAM_FIELDS; by)
    x = initial_parameters(p)
    nf = length(BEAM_FIELDS)
    en = bt.enabled[]
    return map(enumerate(p.keys)) do (g, k)
        mem = p.members[g]
        angles = port_angles.(bs[mem])
        on = count(i -> en[i], mem)
        cone = _mirror(by)
        θ = sum(a -> cone ? min(a[1], 180 - a[1]) : a[1], angles) / length(mem)
        φ = cone ? NaN : sum(last, angles) / length(mem)
        (; key = string(k), count = length(mem), members = mem, theta = θ, phi = φ,
         z = x[(g-1)*nf+1], r = x[(g-1)*nf+2], phi_offset = x[(g-1)*nf+3],
         defocus = x[(g-1)*nf+4], power = x[(g-1)*nf+5],
         enabled = on == length(mem) ? true : on == 0 ? false : missing)
    end
end

"""
    edit_beams!(bt::BeamTable, row, field, text; by = bt.by[]) -> Bool

Set `field` (one of [`BEAM_FIELDS`](@ref)) of every beam in `row` of
[`beam_rows`](@ref) (with grouping `by`) from text-field text. Only that field
changes. A rejected value sets `error` and returns `false`.
"""
function edit_beams!(bt::BeamTable, row::Integer, field::Symbol, text; by::Symbol = bt.by[])
    field in BEAM_FIELDS || throw(ArgumentError("beam field must be one of $BEAM_FIELDS"))
    try
        v = _parse_number(text, field === :power ? "power" : "$field"; nonnegative = field === :power)
        p = _pointing(bt, (field,); by)
        1 <= row <= length(p.keys) || throw(ArgumentError("no beam-table row $row"))
        x = initial_parameters(p)
        x[row] == v && (bt.error[] = ""; return true)
        x[row] = v
        # only the row's beams change (the others are not even re-rounded)
        new = apply(p, x)
        bs = copy(bt.beams[])
        for i in p.members[row]
            bs[i] = _float_beam(new[i])
        end
        bt.beams[] = bs
        bt.error[] = ""
        return true
    catch err
        err isa ArgumentError || rethrow()
        bt.error[] = _errmsg(err)
        return false
    end
end

"""
    set_enabled!(bt::BeamTable, row, on::Bool; by = bt.by[])

Enable or disable every beam in `row` of [`beam_rows`](@ref) (with grouping `by`).
Disabled beams stay in the table and deliver no power (e.g. a dropped quad or bundle).
"""
function set_enabled!(bt::BeamTable, row::Integer, on::Bool; by::Symbol = bt.by[])
    rows = beam_rows(bt; by)
    1 <= row <= length(rows) || throw(ArgumentError("no beam-table row $row"))
    en = copy(bt.enabled[])
    en[rows[row].members] .= on
    _assign!(bt.enabled, en)
    return bt
end

"""
    active_beams(bt::BeamTable) -> Vector{Beam}
    display_beams(bt::BeamTable) -> Vector{Beam}

The enabled beams; every beam with the disabled ones at zero power (same ports and
order as the table, which the drive's deposition operator is built for).
"""
active_beams(bt::BeamTable) = bt.beams[][bt.enabled[]]
display_beams(bt::BeamTable) =
    [on ? b : remake(b; power = 0.0) for (b, on) in zip(bt.beams[], bt.enabled[])]

"""
    import_beams!(bt::BeamTable, path)

Read a beam table CSV with `read_beam_table` (canonical columns, see
`write_beam_table`); the table becomes `:custom`.
"""
function import_beams!(bt::BeamTable, path::AbstractString)
    bs = read_beam_table(path)
    isempty(bs) && throw(ArgumentError("$(basename(path)) holds no beams"))
    set_beams!(bt, bs; preset = :custom)
    bt.path[] = String(path)
    return bt
end

"""
    export_beams(bt::BeamTable, path; enabled_only = false) -> path

Write the table (or only the enabled beams) with `write_beam_table`.
"""
function export_beams(bt::BeamTable, path::AbstractString; enabled_only::Bool = false)
    write_beam_table(path, enabled_only ? active_beams(bt) : bt.beams[])
    bt.path[] = String(path)
    return path
end

# ---------------------------------------------------------------- VisRad files

# Collects the warnings logged while reading a file (the VisRad readers report
# unmatched beams, port-angle mismatches and the like with @warn).
struct _WarningCollector <: Base.CoreLogging.AbstractLogger
    messages::Vector{String}
end
Base.CoreLogging.min_enabled_level(::_WarningCollector) = Base.CoreLogging.Warn
Base.CoreLogging.shouldlog(::_WarningCollector, args...) = true
Base.CoreLogging.catch_exceptions(::_WarningCollector) = false
function Base.CoreLogging.handle_message(c::_WarningCollector, level, msg, _module, group, id,
                                         file, line; kwargs...)
    level < Base.CoreLogging.Warn && return nothing
    s = string(msg)
    kw = [(k, v) for (k, v) in kwargs if k !== :maxlog]
    isempty(kw) ||
        (s *= " (" * join(("$k = $(v isa AbstractVector ? join(first(v, 12), ", ") * (length(v) > 12 ? ", …" : "") : v)"
                           for (k, v) in kw), "; ") * ")")
    push!(c.messages, length(s) > 400 ? s[1:prevind(s, 400)] * "…" : s)
    return nothing
end

function _with_warnings(f)
    c = _WarningCollector(String[])
    v = Base.CoreLogging.with_logger(f, c)
    return v, unique(c.messages)
end

# Apply beams read from a VisRad file with the table as template. A beam the file
# turns off (absent from it, or "Off") becomes disabled at its previous power, so it
# can be switched back on; the others are enabled, as VisRad turns on the beams the
# file has data for.
function _apply_visrad!(bt::BeamTable, read, path::AbstractString)
    old = bt.beams[]
    isempty(old) && throw(ArgumentError("load a beam preset first: the VisRad file is applied to the table"))
    new, warns = _with_warnings(() -> read(old))
    length(new) == length(old) || throw(ArgumentError("$(basename(path)) did not match the beam table"))
    en = copy(bt.enabled[])
    bs = Beam{Float64}[]
    noff = 0
    for (i, (a, b)) in enumerate(zip(old, new))
        b = _float_beam(b)
        if b.power == 0 && a.power > 0
            en[i] = false
            noff += 1
            b = remake(b; power = a.power)
        elseif b.power > 0
            en[i] = true
        end
        push!(bs, b)
    end
    set_beams!(bt, bs; preset = bt.preset[], enabled = en, cluster = bt.cluster[])
    noff > 0 && pushfirst!(warns, (noff == 1 ? "1 beam of the table is" : "$noff beams of the table are") *
                                  " off or absent in $(basename(path)): disabled (power kept)")
    bt.warnings[] = warns
    bt.path[] = String(path)
    return bt
end

"""
    import_visrad_sst!(bt::BeamTable, path)

Apply a NIF SST pointing CSV (`read_visrad_sst`, as VisRad exports it) to the table:
each beam in the file gets its aim point and defocus (mm; quads take the mean of their
four beams); powers and spots stay. As in VisRad, beams absent from the file are
turned off: here they are disabled at their power. Warnings land in `warnings`.
"""
import_visrad_sst!(bt::BeamTable, path::AbstractString) =
    _apply_visrad!(bt, old -> read_visrad_sst(path; template = old, missing = :off), path)

"""
    import_visrad_summary!(bt::BeamTable, path)

Apply a VisRad Laser Beam Summary CSV (`read_visrad_summary`) to the table: pointing,
defocus, power and spot of each beam in the file (port angles are checked against the
table's, with a warning). Beams that are "Off" or absent are disabled at their power.
Warnings land in `warnings`.
"""
import_visrad_summary!(bt::BeamTable, path::AbstractString) =
    _apply_visrad!(bt, old -> read_visrad_summary(path; template = old, missing = :off), path)

"""
    export_visrad_sst(bt::BeamTable, path) -> path

Write the enabled beams' pointing as a NIF SST CSV (`write_visrad_sst`; quads as their
four beams, disabled beams left out so VisRad turns them off).
"""
function export_visrad_sst(bt::BeamTable, path::AbstractString)
    write_visrad_sst(path, display_beams(bt))
    bt.path[] = String(path)
    return path
end

"""
    export_visrad_summary(bt::BeamTable, path) -> path

Write the table as a VisRad Laser Beam Summary CSV (`write_visrad_summary`; disabled
beams "Off", each beam's cluster its group).
"""
function export_visrad_summary(bt::BeamTable, path::AbstractString)
    # quads are written as their four beams: those take the quad's group
    d = Dict(b.name => c for (b, c) in zip(bt.beams[], bt.cluster[]))
    dq = Dict(b.quad => c for (b, c) in zip(bt.beams[], bt.cluster[]) if !isempty(b.quad))
    write_visrad_summary(path, display_beams(bt);
                         cluster = b -> string(get(d, b.name, get(dq, b.quad, b.group))))
    bt.path[] = String(path)
    return path
end

"""
    beams_summary(bt::BeamTable) -> String

Beam count, enabled count, and total enabled power.
"""
function beams_summary(bt::BeamTable)
    n = length(bt.beams[])
    n == 0 && return "no beams"
    act = active_beams(bt)
    P = sum(b.power for b in act; init = 0.0)
    src = get(PRESET_LABELS, bt.preset[], String(bt.preset[]))
    return @sprintf("%d of %d beams on · %.4g TW · %s", length(act), n, P, src)
end

"Ports of a beam list (what a deposition operator is tied to)."
_ports(bs) = [b.port for b in bs]

# ---------------------------------------------------------------- by name

"""
    beam_node(bt::BeamTable, by, key) -> Union{Nothing,NamedTuple}

The [`beam_rows`](@ref) row named `key` (a cone, quad or beam name) under grouping
`by`, with its index as `row`.
"""
function beam_node(bt::BeamTable, by::Symbol, key::AbstractString)
    by === :quad && any(b -> isempty(b.quad), bt.beams[]) && return nothing
    rows = beam_rows(bt; by)
    i = findfirst(r -> r.key == key, rows)
    return i === nothing ? nothing : merge(rows[i], (; row = i))
end

"""
    edit_beam_node!(bt, by, key, field, text) -> Bool
    set_node_enabled!(bt, by, key, on::Bool)

[`edit_beams!`](@ref) and [`set_enabled!`](@ref) on the cone, quad or beam named `key`.
"""
function edit_beam_node!(bt::BeamTable, by::Symbol, key::AbstractString, field::Symbol, text)
    n = beam_node(bt, by, key)
    n === nothing && (bt.error[] = "no $by named $key"; return false)
    return edit_beams!(bt, n.row, field, text; by)
end

function set_node_enabled!(bt::BeamTable, by::Symbol, key::AbstractString, on::Bool)
    n = beam_node(bt, by, key)
    n === nothing && throw(ArgumentError("no $by named $key"))
    return set_enabled!(bt, n.row, on; by)
end

# ---------------------------------------------------------------- sets of beams

"Editable fields of a set of beams: the `BEAM_FIELDS` and the Cartesian aim point."
const SET_FIELDS = (BEAM_FIELDS..., :aim_x, :aim_y, :aim_z)

function _check_set(bt::BeamTable, idx)
    isempty(idx) && throw(ArgumentError("no beams selected"))
    all(i -> 1 <= i <= length(bt.beams[]), idx) || throw(ArgumentError("no such beam"))
    return collect(Int, idx)
end

# A set spanning both hemispheres describes z as a cone does (mirrored); a set in one
# hemisphere uses the aim's own z.
_set_mirror(bs, idx) = any(i -> bs[i].port[3] > 0, idx) && any(i -> bs[i].port[3] < 0, idx)

function _set_pointing(bt::BeamTable, idx, fields)
    S = Set(bt.beams[][i].name for i in idx)
    p = GroupPointing(bt.beams[]; by = b -> b.name in S, fields,
                      mirror = _set_mirror(bt.beams[], idx))
    return p, findfirst(==(true), p.keys)
end

"""
    beam_set_values(bt::BeamTable, idx) -> NamedTuple

Mean `z`, `r`, `phi` (tangential offset), `defocus` (mm) and `power` (TW per beam) of
the beams `idx`, as [`edit_beam_set!`](@ref) reads them (z mirrored when the set spans
both hemispheres: `mirrored`), their mean Cartesian `aim`, `count`, and `enabled`
(`true`, `false`, or `missing` for some).
"""
function beam_set_values(bt::BeamTable, idx)
    idx = _check_set(bt, idx)
    p, g = _set_pointing(bt, idx, BEAM_FIELDS)
    x = initial_parameters(p)
    nf = length(BEAM_FIELDS)
    v(f) = x[(g-1)*nf+findfirst(==(f), BEAM_FIELDS)]
    bs = bt.beams[]
    aim = sum(bs[i].aim for i in idx) / length(idx)
    on = count(i -> bt.enabled[][i], idx)
    return (; z = v(:z), r = v(:r), phi = v(:phi), defocus = v(:defocus), power = v(:power), aim,
            count = length(idx), mirrored = _set_mirror(bs, idx),
            enabled = on == length(idx) ? true : on == 0 ? false : missing)
end

"""
    edit_beam_set!(bt::BeamTable, idx, field, text) -> Bool

Set `field` (one of [`SET_FIELDS`](@ref)) of every beam in `idx` from text: a
[`BEAM_FIELDS`](@ref) value as a group edit does (only that field changes), or one
coordinate of the aim point (`:aim_x`, `:aim_y`, `:aim_z`, mm), which every beam of the
set then shares, as VisRad points a cluster at a reference point. A rejected value sets
`error` and returns `false`.
"""
function edit_beam_set!(bt::BeamTable, idx, field::Symbol, text)
    field in SET_FIELDS || throw(ArgumentError("beam field must be one of $SET_FIELDS"))
    try
        idx = _check_set(bt, idx)
        v = _parse_number(text, field === :power ? "power" : "$field"; nonnegative = field === :power)
        bs = copy(bt.beams[])
        if field in BEAM_FIELDS
            p, g = _set_pointing(bt, idx, (field,))
            x = initial_parameters(p)
            x[g] == v && (bt.error[] = ""; return true)
            x[g] = v
            new = apply(p, x)
            for i in idx
                bs[i] = _float_beam(new[i])
            end
        else
            k = field === :aim_x ? 1 : field === :aim_y ? 2 : 3
            for i in idx
                a = collect(bs[i].aim)
                a[k] = v
                bs[i] = remake(bs[i]; aim = Vec3(a...))
            end
        end
        bt.beams[] = bs
        bt.error[] = ""
        return true
    catch err
        err isa ArgumentError || rethrow()
        bt.error[] = _errmsg(err)
        return false
    end
end

"""
    set_set_enabled!(bt::BeamTable, idx, on::Bool)
"""
function set_set_enabled!(bt::BeamTable, idx, on::Bool)
    idx = _check_set(bt, idx)
    en = copy(bt.enabled[])
    en[idx] .= on
    _assign!(bt.enabled, en)
    return bt
end

"Group (cluster) names of the table, in order of first appearance."
clusters(bt::BeamTable) = unique(bt.cluster[])

"""
    set_cluster!(bt::BeamTable, idx, name) -> Bool

Move the beams `idx` into the group `name` (a new group if no beam is in it yet), as
VisRad's "Move to Cluster". A group left without beams disappears.
"""
function set_cluster!(bt::BeamTable, idx, name)
    try
        idx = _check_set(bt, idx)
        s = strip(String(name))
        isempty(s) && throw(ArgumentError("give the group a name"))
        c = copy(bt.cluster[])
        c[idx] .= Symbol(s)
        _assign!(bt.cluster, c)
        bt.error[] = ""
        return true
    catch err
        err isa ArgumentError || rethrow()
        bt.error[] = _errmsg(err)
        return false
    end
end

"""
    rename_cluster!(bt::BeamTable, old, new) -> Bool
"""
rename_cluster!(bt::BeamTable, old, new) =
    set_cluster!(bt, findall(==(Symbol(old)), bt.cluster[]), new)

"""
    split_quad!(bt::BeamTable, quad) -> Bool

Replace a NIF quad that the table holds as one beam (`nif_quads`) by its four beams
(`nif_beams`: their own ports and f/20.7 cones), each with the quad's aim, defocus,
detuning, spot and group, and a quarter of its power, so its beams can be pointed
apart (split quads). The quad's enabled flag carries over. Ports change, so the drive
must be computed again.
"""
function split_quad!(bt::BeamTable, quad::AbstractString)
    try
        bs = bt.beams[]
        i = findfirst(b -> b.quad == quad, bs)
        i === nothing && throw(ArgumentError("no quad $quad in the table"))
        count(b -> b.quad == quad, bs) == 1 ||
            throw(ArgumentError("quad $quad is already split into its beams"))
        q = bs[i]
        parts = [remake(b; aim = q.aim, defocus = q.defocus, dlambda = q.dlambda, spot = q.spot,
                        power = q.power / 4, group = q.group)
                 for b in nif_beams() if b.quad == quad]
        length(parts) == 4 || throw(ArgumentError("$quad is not a NIF quad"))
        nb = vcat(bs[1:i-1], parts, bs[i+1:end])
        en = vcat(bt.enabled[][1:i-1], fill(bt.enabled[][i], 4), bt.enabled[][i+1:end])
        cl = vcat(bt.cluster[][1:i-1], fill(bt.cluster[][i], 4), bt.cluster[][i+1:end])
        rf = vcat(bt.aim_ref[][1:i-1], fill(bt.aim_ref[][i], 4), bt.aim_ref[][i+1:end])
        rp = vcat(bt.ref_point[][1:i-1], fill(bt.ref_point[][i], 4), bt.ref_point[][i+1:end])
        set_beams!(bt, nb; preset = bt.preset[], enabled = en, cluster = cl, aim_ref = rf,
                   ref_point = rp)
        return true
    catch err
        err isa ArgumentError || rethrow()
        bt.error[] = _errmsg(err)
        return false
    end
end

"""
    merge_quad!(bt::BeamTable, quad) -> Bool

Replace the four beams of a NIF quad by the quad as one beam (`nif_quads`): their mean
aim and defocus, their total power, the first beam's spot, group and enabled flag.
"""
function merge_quad!(bt::BeamTable, quad::AbstractString)
    try
        bs = bt.beams[]
        idx = findall(b -> b.quad == quad, bs)
        length(idx) == 4 || throw(ArgumentError("quad $quad is not split into four beams"))
        tmpl = findfirst(b -> b.quad == quad, nif_quads())
        tmpl === nothing && throw(ArgumentError("$quad is not a NIF quad"))
        b1 = bs[idx[1]]
        q = remake(nif_quads()[tmpl]; aim = sum(bs[i].aim for i in idx) / 4,
                   defocus = sum(bs[i].defocus for i in idx) / 4, dlambda = b1.dlambda,
                   power = sum(bs[i].power for i in idx), spot = b1.spot, group = b1.group)
        keep = setdiff(eachindex(bs), idx[2:end])
        nb = copy(bs)
        nb[idx[1]] = q
        set_beams!(bt, nb[keep]; preset = bt.preset[], enabled = bt.enabled[][keep],
                   cluster = bt.cluster[][keep], aim_ref = bt.aim_ref[][keep],
                   ref_point = bt.ref_point[][keep])
        return true
    catch err
        err isa ArgumentError || rethrow()
        bt.error[] = _errmsg(err)
        return false
    end
end
