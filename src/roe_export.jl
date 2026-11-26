# Tables of a rules-of-engagement report (roe.jl): rows, CSV and Markdown.

const _ROE_DECISION = Dict(:red => "Critical", :amber => "Ask RI", :green => "Auto-Drop")

_roe_unit_names(r::ROEReport, o::ROEOutage) = [r.units[i].name for i in o.units]

# The pair caveats of unit `i` as "partner: tier (criterion)".
function _roe_caveat_text(r::ROEReport, i::Int)
    out = String[]
    for k in r.caveats[i]
        p = r.pairs[k]
        partner = r.units[p.units[1] == i ? p.units[2] : p.units[1]].name
        push!(out, string(partner, ": ", p.tier, " (", r.criteria[p.binding].name, ")"))
    end
    return out
end

"""
    roe_table(report; kind = nothing) -> Vector{NamedTuple}

The rows of a rules-of-engagement report ([`ROEReport`](@ref)), one per unit, sorted by
role, then from red to green, then by name. `kind = :beam` or `:quad` keeps one kind.
Each row has:

- `role`, `unit`, `kind` and `beams` (the beam names);
- `tier`, and `decision`: `"Critical"` (red, hold the shot), `"Ask RI"` (amber) or
  `"Auto-Drop"` (green);
- `single`, the tier of the unit alone, which differs from `tier` when a pair raised it;
- `binding`, the criterion that sets the tier of the single outage, with `score`,
  `green` and `amber`, its score and thresholds;
- `scores`, a named tuple of the scores on every criterion;
- `companion` and `companion_tier`, the pre-planned second drop and the tier of the two
  together, or `nothing`;
- `caveats`, the pairs that are worse than both of their units alone, as text.
"""
function roe_table(r::ROEReport; kind::Union{Nothing,Symbol} = nothing)
    names = Tuple(c.name for c in r.criteria)
    rows = NamedTuple[]
    for (i, u) in enumerate(r.units)
        kind === nothing || u.kind === kind || continue
        o = r.singles[i]
        c = r.criteria[o.binding]
        comp = r.companions[i]
        partner = comp === nothing ? nothing :
                  r.units[comp.units[1] == i ? comp.units[2] : comp.units[1]].name
        push!(rows, (; role = u.role, unit = u.name, kind = u.kind, beams = r.beams[u.beams],
                     tier = r.tiers[i], decision = _ROE_DECISION[r.tiers[i]], single = o.tier,
                     binding = c.name, score = o.scores[o.binding], green = c.green, amber = c.amber,
                     scores = NamedTuple{names}(Tuple(o.scores)), companion = partner,
                     companion_tier = comp === nothing ? nothing : comp.tier,
                     caveats = _roe_caveat_text(r, i)))
    end
    roles = unique(u.role for u in r.units)
    order(row) = (findfirst(==(row.role), roles), row.kind === :quad ? 0 : 1, -_roe_rank(row.tier), row.unit)
    return sort!(rows; by = order)
end

"""
    write_roe_csv(path, report; kind = nothing) -> path

Writes the rows of [`roe_table`](@ref) as CSV, one line per unit: `role`, `unit`, `kind`,
`beams`, `tier`, `decision`, `single`, `binding`, `score`, `green`, `amber`, one column
per criterion with its score, `companion`, `companion_tier` and `caveats`. Lists within
a field are separated by semicolons.
"""
function write_roe_csv(path::AbstractString, r::ROEReport; kind::Union{Nothing,Symbol} = nothing)
    open(path, "w") do io
        head = ["role", "unit", "kind", "beams", "tier", "decision", "single", "binding", "score",
                "green", "amber", (string(c.name) for c in r.criteria)..., "companion",
                "companion_tier", "caveats"]
        println(io, join(head, ","))
        for row in roe_table(r; kind)
            f = Any[row.role, row.unit, row.kind, join(row.beams, ";"), row.tier, row.decision,
                    row.single, row.binding, row.score, row.green, row.amber, values(row.scores)...,
                    something(row.companion, ""), something(row.companion_tier, ""),
                    join(row.caveats, "; ")]
            println(io, join((_csvfield(x) for x in f), ","))
        end
    end
    return path
end

_roe_number(x::Real) = isfinite(x) ? @sprintf("%.4g", x) : string(x)
_roe_number(x) = "(" * join((_roe_number(y) for y in x), ", ") * ")"

"""
    roe_markdown(report; title = "Beam rules of engagement", kind = nothing) -> String

The rules-of-engagement sheet of a report as Markdown: the criteria with their thresholds
and nominal values, then for each role the budget and a table of the units from critical
to auto-drop, with the binding criterion, its score against the thresholds, the scores on
every criterion, the companion drop and the pair caveats. [`write_roe_markdown`](@ref)
writes it to a file.
"""
function roe_markdown(r::ROEReport; title::AbstractString = "Beam rules of engagement",
                      kind::Union{Nothing,Symbol} = nothing)
    io = IOBuffer()
    println(io, "# ", title, "\n")
    println(io, "An outage is a beam or quad at zero power with every other beam as planned. ",
            "Critical units hold the shot, Ask RI units need approval before the drop, and ",
            "Auto-Drop units may drop without a request.\n")
    println(io, "## Criteria\n")
    println(io, "| Criterion | Case | Measure | Green ≤ | Amber ≤ | Unit | Nominal | Path |")
    println(io, "|---|---|---|---|---|---|---|---|")
    for (c, v) in zip(r.criteria, r.baseline)
        j = c.case === :default ? 1 : findfirst(==(c.case), r.cases)
        println(io, "| ", c.name, " | ", r.cases[j], " | ", c.measure, " | ", _roe_number(c.green), " | ",
                _roe_number(c.amber), " | ", c.label, " | ", _roe_number(v), " | ",
                r.linear[j] ? "linear deposition" : "full evaluation", " |")
    end
    rows = roe_table(r; kind)
    for role in unique(row.role for row in rows)
        println(io, "\n## Role: ", role, "\n")
        for b in r.budgets
            b.role === role && (kind === nothing || b.kind === kind) || continue
            println(io, "- Budget (", b.kind, "s): ", b.text, ".")
            isempty(b.exceptions) ||
                println(io, "  Pairs that leave green: ", join((a * " + " * c for (a, c) in b.exceptions), ", "), ".")
            isempty(b.triple_exceptions) ||
                println(io, "  Triples that leave green: ", join((join(t, " + ") for t in b.triple_exceptions), ", "), ".")
        end
        println(io)
        println(io, "| Unit | Kind | Decision | Binding | Score | Green ≤ | Amber ≤ | ",
                join((string(c.name) for c in r.criteria), " | "), " | Companion | Pair caveats |")
        println(io, "|---|---|---|---|---|---|---|", "---|"^length(r.criteria), "---|---|")
        for row in rows
            row.role === role || continue
            comp = row.companion === nothing ? "" : string(row.companion, ": ", row.companion_tier)
            raised = row.tier === row.single ? "" : " (alone: $(row.single))"
            println(io, "| ", row.unit, " | ", row.kind, " | ", row.decision, raised, " | ", row.binding,
                    " | ", _roe_number(row.score), " | ", _roe_number(row.green), " | ",
                    _roe_number(row.amber), " | ", join((_roe_number(s) for s in values(row.scores)), " | "),
                    " | ", comp, " | ", join(row.caveats, "; "), " |")
        end
    end
    return String(take!(io))
end

"""
    write_roe_markdown(path, report; kwargs...) -> path

Writes [`roe_markdown`](@ref) of the report to `path`.
"""
function write_roe_markdown(path::AbstractString, r::ROEReport; kwargs...)
    write(path, roe_markdown(r; kwargs...))
    return path
end
