# This file covers the beam rules of engagement in the window: the scan of the drive
# scenario against the default criteria of the design, the tier of each beam for the 3D
# view, the line about a selected beam, and the port map of the Outages tab with the
# hypothetical outage of the units picked on it.

"""
    ROEState

The beam rules of engagement of a session (`Briareus.roe_scan`). `report` is the report
of the last scan (a core `ROEReport`, or `nothing`), `notes` the lines about the criteria
that were used, `status` a line about the scan and `inputs` the target, beams and drive
scenario of that scan. `picked` lists the units (beams or quads, by name) of a
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
`:amber`, `:red`, or `:none` for a mark that the scan left out) and whether it is
`picked` for the hypothetical outage. Empty before a scan.
"""
function roe_marks(rs::ROEState)
    r = rs.report[]
    r === nothing && return NamedTuple[]
    beams = rs.inputs[][2]
    by = any(u -> u.kind === :quad, r.units) ? :quad : :beam
    tier(name) = (i = findfirst(u -> u.name == name, r.units); i === nothing ? :none : r.tiers[i])
    return [(; m..., tier = tier(m.name), picked = m.name in rs.picked[]) for m in port_map(beams; by)]
end

"""
    roe_port_drawing(rs::ROEState) -> String

The port map of the Outages tab as the JSON text that its canvas draws: one entry per
mark ([`roe_marks`](@ref)) with its `name`, its position `x`, `y` on the unit disk of its
hemisphere, `upper` (whether it lies in the upper hemisphere), its `tier` and whether it
is `picked`. Empty before a scan.
"""
function roe_port_drawing(rs::ROEState)
    ms = roe_marks(rs)
    isempty(ms) && return ""
    item(m) = @sprintf("{\"name\":\"%s\",\"x\":%.5f,\"y\":%.5f,\"upper\":%s,\"tier\":\"%s\",\"picked\":%s}",
                       escape_string(m.name), m.x, m.y, m.hemisphere === :upper, m.tier, m.picked)
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

A line about the hypothetical outage of the Outages tab: the picked units, the tier and
the criterion that binds it. Empty without picked units.
"""
function roe_outage_text(rs::ROEState)
    isempty(rs.picked[]) && return ""
    o, r = rs.outage[], rs.report[]
    head = "outage of " * join(rs.picked[], ", ")
    (o === nothing || r === nothing) && return head
    c = r.criteria[o.binding]
    return @sprintf("%s (%d beams): %s, bound by %s at %.3g%s", head, length(o.beams), o.tier, c.name,
                    o.scores[o.binding], isempty(c.label) ? "" : " " * c.label)
end
