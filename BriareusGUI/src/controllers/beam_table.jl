# The Beams step: a facility preset or an imported beam table, edited per cone
# (group), quad, or beam. Edits go through `GroupPointing` with only the edited field,
# so a group edit sets that field on every member and leaves the others as they were.

"""
Beam presets of the Beams step:

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

"Row groupings of the beam table: per cone (`:group`), per `:quad`, per `:beam`."
const BEAM_GROUPINGS = (:group, :quad, :beam)

"Editable beam-table columns (`GroupPointing` fields): aim z, r, tangential offset, defocus (mm), power (TW)."
const BEAM_FIELDS = (:z, :r, :phi, :defocus, :power)

# Published NIF cone pointing (r toward the port azimuth, z; mm) used by
# examples/nif_hohlraum.jl for a 10.1 mm hohlraum.
const NIF_CONE_POINTING = Dict(:inner23 => (-1.09, 3.09), :inner30 => (-0.60, 3.63),
                               :outer44 => (0.0, 5.20), :outer50 => (0.0, 5.01))

"""
    BeamTable(beams = Beam{Float64}[]; preset = :custom)

State of the Beams step: `beams` (every beam of the table, in table order), `enabled`
(one flag per beam; disabled beams are kept but deliver no power), `preset` (one of
[`BEAM_PRESETS`](@ref), or `:custom` for an imported or scripted table), `by` (row
grouping, one of [`BEAM_GROUPINGS`](@ref)), `path` (the last imported or exported
file), and `error` (the last rejected edit).
"""
struct BeamTable
    beams::Observable{Vector{Beam{Float64}}}
    enabled::Observable{Vector{Bool}}
    preset::Observable{Symbol}
    by::Observable{Symbol}
    path::Observable{String}
    error::Observable{String}
end

function BeamTable(beams::AbstractVector = Beam{Float64}[]; preset::Symbol = :custom)
    bs = Beam{Float64}[_float_beam(b) for b in beams]
    return BeamTable(Observable(bs), Observable(fill(true, length(bs))), Observable(preset),
                     Observable(:group), Observable(""), Observable(""))
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
    load_preset!(bt::BeamTable, preset; target_kind, params)

Replace the table with a preset's beams (see [`preset_beams`](@ref)), all enabled.
"""
function load_preset!(bt::BeamTable, preset::Symbol; kwargs...)
    set_beams!(bt, preset_beams(preset; kwargs...); preset)
    return bt
end

"""
    set_beams!(bt::BeamTable, beams; preset = :custom, enabled = all)

Replace the table's beams (all enabled unless `enabled` is given).
"""
function set_beams!(bt::BeamTable, beams::AbstractVector; preset::Symbol = :custom,
                    enabled::AbstractVector{Bool} = fill(true, length(beams)))
    length(enabled) == length(beams) || throw(DimensionMismatch("one enabled flag per beam"))
    bt.enabled.val = collect(Bool, enabled)          # observers see consistent lengths
    bt.preset[] = preset
    bt.by[] === :quad && any(b -> isempty(b.quad), beams) && (bt.by[] = :group)
    bt.beams[] = Beam{Float64}[_float_beam(b) for b in beams]
    notify(bt.enabled)
    bt.error[] = ""
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

# Rows group by cone with z mirrored top/bottom (one value describes a cone); per quad
# or beam, z is the aim's own coordinate.
_mirror(by::Symbol) = by === :group
_pointing(bt::BeamTable, fields) =
    GroupPointing(bt.beams[]; by = bt.by[], fields, mirror = _mirror(bt.by[]))

"""
    beam_rows(bt::BeamTable) -> Vector{NamedTuple}

One row per group of the current grouping: `key` (cone, quad, or beam name), `count`
(beams), `members` (beam indices), `theta`/`phi` (port angles, degrees: per cone the
mean polar angle from the nearer pole and `NaN` azimuth; per quad or beam the mean
port angles), the mean
`z`, `r`, `phi_offset`, `defocus` (mm) and `power` (TW per beam), and `enabled`
(`true`, `false`, or `missing` when only some members are enabled).
"""
function beam_rows(bt::BeamTable)
    bs = bt.beams[]
    isempty(bs) && return NamedTuple[]
    p = _pointing(bt, BEAM_FIELDS)
    x = initial_parameters(p)
    nf = length(BEAM_FIELDS)
    en = bt.enabled[]
    return map(enumerate(p.keys)) do (g, k)
        mem = p.members[g]
        angles = port_angles.(bs[mem])
        on = count(i -> en[i], mem)
        cone = _mirror(bt.by[])
        θ = sum(a -> cone ? min(a[1], 180 - a[1]) : a[1], angles) / length(mem)
        φ = cone ? NaN : sum(last, angles) / length(mem)
        (; key = string(k), count = length(mem), members = mem, theta = θ, phi = φ,
         z = x[(g-1)*nf+1], r = x[(g-1)*nf+2], phi_offset = x[(g-1)*nf+3],
         defocus = x[(g-1)*nf+4], power = x[(g-1)*nf+5],
         enabled = on == length(mem) ? true : on == 0 ? false : missing)
    end
end

"""
    edit_beams!(bt::BeamTable, row, field, text) -> Bool

Set `field` (one of [`BEAM_FIELDS`](@ref)) of every beam in `row` of
[`beam_rows`](@ref) from text-field text. Only that field changes. A rejected value
sets `error` and returns `false`.
"""
function edit_beams!(bt::BeamTable, row::Integer, field::Symbol, text)
    field in BEAM_FIELDS || throw(ArgumentError("beam field must be one of $BEAM_FIELDS"))
    try
        v = _parse_number(text, field === :power ? "power" : "$field"; nonnegative = field === :power)
        p = _pointing(bt, (field,))
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
    set_enabled!(bt::BeamTable, row, on::Bool)

Enable or disable every beam in `row` of [`beam_rows`](@ref). Disabled beams stay in
the table and deliver no power (e.g. a dropped quad or bundle).
"""
function set_enabled!(bt::BeamTable, row::Integer, on::Bool)
    rows = beam_rows(bt)
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
order as the table, which the Drive step's deposition operator is built for).
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
