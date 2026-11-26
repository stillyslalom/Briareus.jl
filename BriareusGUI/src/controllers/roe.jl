# This file covers the beam rules of engagement in the window: the scan of the drive
# scenario against the default criteria of the design, the tier of each beam for the 3D
# view, the line about a selected beam, and the port map of the Outages tab with the
# hypothetical outage of the units picked on it.

"""
    ROEState

The beam rules of engagement of a session (`Briareus.roe_scan`). `report` is the report
of the last scan (a core `ROEReport`, or `nothing`), `notes` the lines about the criteria
that were used, `status` a line about the scan and `inputs` the target, beams, drive
scenario and `stands_for` ([`roe_stands_for`](@ref)) of that scan. `picked` lists the units (beams or quads, by name) of a
hypothetical outage, and `outage` its result (`Briareus.roe_outage`, or `nothing`).
"""
struct ROEState
    report::Observable{Union{Nothing,ROEReport}}
    notes::Observable{Vector{String}}
    status::Observable{String}
    inputs::Base.RefValue{Any}
    picked::Observable{Vector{String}}
    outage::Observable{Union{Nothing,ROEOutage}}
end

ROEState() = ROEState(Observable{Union{Nothing,ROEReport}}(nothing), Observable(String[]), Observable(""), Ref{Any}(nothing),
                      Observable(String[]), Observable{Union{Nothing,ROEOutage}}(nothing))

"""
    roe_stands_for(bt::BeamTable) -> Dict{String,Int}

The number of facility beams that each beam of the table stands for, for `stands_for` of
`Briareus.roe_scan`: 4 for each NIF quad held as one beam ([`quads_to_split`](@ref)).
The scan of such a table then scores the loss of each quad and of one of its four beams,
as a quarter of the power of the quad. Beams that stand for themselves are left out.
"""
function roe_stands_for(bt::BeamTable)
    bs = bt.beams[]
    return Dict{String,Int}(bs[i].name => 4 for i in quads_to_split(bt))
end

# The units of a report that hold beam `name`: the beam unit and the quad unit, where present.
function _roe_units(r::ROEReport, name::AbstractString)
    k = findfirst(==(name), r.beams)
    k === nothing && return Int[]
    return [i for i in eachindex(r.units) if k in r.units[i].beams]
end

"""
    roe_tiers(rs::ROEState) -> Dict{String,Symbol}

The tier (`:green`, `:amber` or `:red`) of each beam of the last scan, for the coloring
of the 3D view: the tier of the beam as a unit, or of its quad when the scan has no unit
for single beams. Empty before a scan.
"""
function roe_tiers(rs::ROEState)
    r = rs.report[]
    out = Dict{String,Symbol}()
    r === nothing && return out
    for name in r.beams
        us = _roe_units(r, name)
        isempty(us) && continue
        i = something(findfirst(i -> r.units[i].kind === :beam, us), 1)
        out[name] = r.tiers[us[i]]
    end
    return out
end

"""
    roe_beam_text(rs::ROEState, name) -> String

A line about beam `name` in the last scan: for the beam and for its quad, the tier, the
criterion that binds it with its score and thresholds, and how many pairs are worse
together. Empty before a scan or for a beam that the scan left out.
"""
function roe_beam_text(rs::ROEState, name::AbstractString)
    r = rs.report[]
    r === nothing && return ""
    parts = String[]
    for i in _roe_units(r, name)
        u, o = r.units[i], r.singles[i]
        c = r.criteria[o.binding]
        s = @sprintf("%s %s: %s, bound by %s at %.3g%s (green ≤ %.3g, amber ≤ %.3g)", u.kind, u.name, r.tiers[i],
                     c.name, o.scores[o.binding], isempty(c.label) ? "" : " " * c.label, c.green, c.amber)
        n = length(r.caveats[i])
        n > 0 && (s *= @sprintf("; %d pair%s worse together", n, n == 1 ? "" : "s"))
        push!(parts, s)
    end
    return join(parts, ". ")
end

"""
    roe_marks(rs::ROEState) -> Vector{NamedTuple}

The marks of the port map of the Outages tab (`Briareus.port_map` of the beams of the
last scan, by quad when the scan has quad units), each with its `tier` (`:green`,
`:amber`, `:red`, or `:none` for a mark that the scan left out), the worst tier `beam` of
the single beams of the mark (the tier of a beam mark, or of the beam units within a
quad, `:none` without them) and whether it or one of those beams is `picked` for the
hypothetical outage. Empty before a scan.
"""
function roe_marks(rs::ROEState)
    r = rs.report[]
    r === nothing && return NamedTuple[]
    beams = rs.inputs[][2]
    by = any(u -> u.kind === :quad, r.units) ? :quad : :beam
    unit(name) = findfirst(u -> u.name == name, r.units)
    tier(name) = (i = unit(name); i === nothing ? :none : r.tiers[i])
    # the beam units within the unit of a mark
    within(i) = r.units[i].kind === :beam ? [i] :
                [j for j in eachindex(r.units) if r.units[j].kind === :beam && r.units[j].beams ⊆ r.units[i].beams]
    function beam_tier(name)
        i = unit(name)
        i === nothing && return :none
        ts = r.tiers[within(i)]
        return isempty(ts) ? :none : ts[argmax(map(t -> t === :green ? 1 : t === :amber ? 2 : 3, ts))]
    end
    function picked(name)
        name in rs.picked[] && return true
        i = unit(name)
        return i !== nothing && any(j -> r.units[j].name in rs.picked[], within(i))
    end
    return [(; m..., tier = tier(m.name), beam = beam_tier(m.name), picked = picked(m.name)) for m in port_map(beams; by)]
end

"""
    roe_port_drawing(rs::ROEState) -> String

The port map of the Outages tab as the JSON text that its canvas draws: one entry per
mark ([`roe_marks`](@ref)) with its `name`, its position `x`, `y` on the unit disk of its
hemisphere, `upper` (whether it lies in the upper hemisphere), its `tier`, the tier
`beam` of its single beams and whether it is `picked`. Empty before a scan.
"""
function roe_port_drawing(rs::ROEState)
    ms = roe_marks(rs)
    isempty(ms) && return ""
    item(m) = @sprintf("{\"name\":\"%s\",\"x\":%.5f,\"y\":%.5f,\"upper\":%s,\"tier\":\"%s\",\"beam\":\"%s\",\"picked\":%s}",
                       escape_string(m.name), m.x, m.y, m.hemisphere === :upper, m.tier, m.beam, m.picked)
    return "[" * join((item(m) for m in ms), ",") * "]"
end

"""
    roe_outage_rows(rs::ROEState) -> Vector{NamedTuple}

The criteria of the hypothetical outage of the Outages tab, one row each: `criterion`,
`score` with its `label`, the `green` and `amber` thresholds and the `tier`, and whether
the criterion `binding` sets the tier of the outage. Empty without an outage.
"""
function roe_outage_rows(rs::ROEState)
    r, o = rs.report[], rs.outage[]
    (r === nothing || o === nothing) && return NamedTuple[]
    return [(; criterion = String(c.name), score = o.scores[k], label = c.label, green = c.green, amber = c.amber,
             tier = o.tiers[k], binding = k == o.binding) for (k, c) in enumerate(r.criteria)]
end

"""
    roe_outage_text(rs::ROEState) -> String

A line about the hypothetical outage of the Outages tab: the picked units, the number of
facility beams lost, the tier and the criterion that binds it. Empty without picked
units.
"""
function roe_outage_text(rs::ROEState)
    isempty(rs.picked[]) && return ""
    o, r = rs.outage[], rs.report[]
    head = "outage of " * join(rs.picked[], ", ")
    (o === nothing || r === nothing) && return head
    c = r.criteria[o.binding]
    held = length(rs.inputs[]) > 3 ? rs.inputs[][4] : Dict{String,Int}()
    n = round(Int, sum(f * get(held, r.beams[b], 1) for (b, f) in zip(o.beams, o.shares)))
    return @sprintf("%s (%d beam%s): %s, bound by %s at %.3g%s", head, n, n == 1 ? "" : "s", o.tier, c.name,
                    o.scores[o.binding], isempty(c.label) ? "" : " " * c.label)
end
