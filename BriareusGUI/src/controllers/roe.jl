# This file covers the beam rules of engagement in the window: the scan of the drive
# scenario against the default criteria of the design, the tier of each beam for the 3D
# view, and the line about a selected beam.

"""
    ROEState

The beam rules of engagement of a session (`Briareus.roe_scan`). `report` is the report
of the last scan (a core `ROEReport`, or `nothing`), `notes` the lines about the criteria
that were used, `status` a line about the scan and `inputs` the beams and target of that
scan.
"""
struct ROEState
    report::Observable{Union{Nothing,ROEReport}}
    notes::Observable{Vector{String}}
    status::Observable{String}
    inputs::Base.RefValue{Any}
end

ROEState() = ROEState(Observable{Union{Nothing,ROEReport}}(nothing), Observable(String[]), Observable(""), Ref{Any}(nothing))

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
