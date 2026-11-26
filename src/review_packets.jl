# Review packets: the 1ω review of a target and the rules-of-engagement sheet of a shot,
# each as one self-contained HTML page with its tables and vector drawings, and as a
# folder of files (HTML, CSV, SVG, JSON) for a reviewer without the package.

const _OMEGA1_RULES = (:clearance_1w, :management_1w, :diagnostic_1w)
const _DISPOSITION_ORDER = (:fail, :review, :specular_clear, :curved_safe, :tight, :clear)

"""
    Omega1Review

The 1ω review of a target, as [`omega1_review`](@ref) returns it:

- `report`, the [`TaLISReport`](@ref) of the check;
- `disposition`, the [`Omega1Disposition`](@ref) of the hardware, or `nothing`;
- `components`, one row per checked component of the disposition: `component`,
  `triangles`, `worst` (the most severe status), `margin` (the smallest 1ω margin, mm)
  and the area (mm²) of its triangles under each status (`fail`, `review`,
  `specular_clear`, `curved_safe`, `tight`, `clear`);
- `threshold`, the margin (mm) below which a passing pair counts as tight;
- `items`, one entry per failing or tight 1ω row, from the most severe: the `row`, its
  [`Omega1Approach`](@ref) as `approach`, and the `slice` across the beamlet
  ([`omega1_slice`](@ref)).
"""
struct Omega1Review
    report::TaLISReport
    disposition::Union{Nothing,Omega1Disposition}
    components::Vector{NamedTuple}
    threshold::Float64
    items::Vector{NamedTuple}
end

function Base.show(io::IO, r::Omega1Review)
    n(st) = count(x -> x.rule in _OMEGA1_RULES && x.status === st, r.report.rows)
    print(io, "Omega1Review(", n(:fail), " failing, ", n(:review), " for review, ", length(r.items),
          " approaches below ", r.threshold, " mm or not passing)")
end

# The keywords of `kw` that `f`-like calls share, by name.
_kw_subset(kw, names) = (; (k => v for (k, v) in pairs(kw) if k in names)...)

"""
    omega1_review(target, beams; threshold = 1.0, disposition = true, window = nothing,
                  kwargs...) -> Omega1Review

Run the TaLIS check of `beams` on `target` ([`talis_check`](@ref), which takes `kwargs`,
`capture` among them) and collect what a reviewer of the unconverted light looks at: the
report, the 1ω disposition of the hardware ([`omega1_disposition`](@ref), with the same
components, tolerance and error set, and `threshold` as its tight threshold) summarized
per component, and for every 1ω row that fails, is left for review or passes with less
than `threshold` (mm) of margin, the closest approach
([`omega1_closest_approach`](@ref)) with its slice ([`omega1_slice`](@ref), of half-width
`window`). `disposition = false` skips the disposition. A target of more than 20,000
triangles gets one verdict per triangle in the disposition, as its facets are fine already.

[`omega1_review_html`](@ref) and [`write_omega1_review`](@ref) write the packet.
"""
function omega1_review(t::Target, beams::AbstractVector{<:Beam}; threshold::Real = 1.0, disposition::Bool = true,
                       window = nothing, kwargs...)
    threshold >= 0 || throw(ArgumentError("threshold must be ≥ 0"))
    report = talis_check(t, beams; kwargs...)
    disp = nothing
    if disposition && get(kwargs, :omega1, true)
        shared = _kw_subset(kwargs, (:diagnostics, :tolerance, :pivot, :pointing, :target_offset, :hardware_tolerance,
                                     :focus, :omega1_margin, :omega1_edge, :omega1_size, :omega1_focus,
                                     :lens_distance, :offsets, :open_ports, :reflection_angle, :cylinder_radius,
                                     :sphere_radius, :curvature, :resolution))
        level = length(t.tri_area) > 20_000 ? 0 : get(kwargs, :level, 2)
        disp = omega1_disposition(t, beams; components = get(kwargs, :omega1_components, nothing),
                                  tight_threshold = threshold, margin_range = max(5.0, threshold), shadows = true,
                                  level, shared...)
    end
    rows = [x for x in report.rows if x.rule in _OMEGA1_RULES &&
            (x.status !== :pass || (x.unit === :mm && x.margin < threshold))]
    sort!(rows; by = x -> (_talis_rank(x), x.margin, x.beam, String(x.component)))
    items = NamedTuple[]
    for row in rows
        a = omega1_closest_approach(t, beams, report, row; offsets = get(kwargs, :offsets, nothing))
        push!(items, (; row, approach = a, slice = omega1_slice(t, a; window)))
    end
    return Omega1Review(report, disp, _disposition_components(t, disp), Float64(threshold), items)
end

function _disposition_components(t::Target, d)
    out = NamedTuple[]
    d === nothing && return out
    for (ci, c) in enumerate(t.components)
        js = [j for j in eachindex(d.status) if t.tri_comp[j] == ci && d.status[j] !== :unchecked]
        isempty(js) && continue
        area(st) = sum((t.tri_area[j] for j in js if d.status[j] === st); init = 0.0)
        worst = argmax(st -> _DISPOSITION_RANK[st], (d.status[j] for j in js))
        push!(out, (; component = c.name, triangles = length(js), worst, margin = minimum(d.margin[j] for j in js),
                    (st => area(st) for st in _DISPOSITION_ORDER)...))
    end
    return out
end

_packet_mm(x::Real) = isfinite(x) ? @sprintf("%.3f", x) : string(x)
_packet_vec(v) = "(" * join((@sprintf("%.3f", x) for x in v), ", ") * ")"
_packet_setting(v) = v isa AbstractVector || v isa Tuple ? join(v, " ") : string(v)

"""
    omega1_review_html(review; title = "1ω review", generated = nothing, images = Dict()) -> String

The 1ω review ([`Omega1Review`](@ref)) as one self-contained HTML page: the counts of the
check by rule and status, each failing or tight approach with the drawing of its slice
([`omega1_slice_svg`](@ref)) and its margin terms, the reflection records, the
disposition of the hardware per component, every row of the check in a table that sorts
on a click of a column head, and the settings and the checks that were not made. The
page refers to no other file and to nothing on the network.

`images` maps the index of an approach in `review.items` to the `src` of a picture shown
beside its slice, such as a rendered view as a `data:` URI. Such a picture needs a
graphics context, so an application supplies it and the package leaves the place empty. `generated` is a
line for the foot of the page, such as a date. Without it the page of the same review is
the same text.
"""
function omega1_review_html(r::Omega1Review; title::AbstractString = "1ω review", generated = nothing,
                            images::AbstractDict = Dict{Int,String}())
    rep = r.report
    s = rep.settings
    io = IOBuffer()
    o1 = [x for x in rep.rows if x.rule in _OMEGA1_RULES]
    lead = string(length(rep.rows), " checks, rules from ", TALIS_SOURCE, ". ", count(x -> x.status === :fail, o1), " 1ω checks fail, ",
                  count(x -> x.status === :review, o1), " need review and ",
                  count(x -> x.status === :pass && x.unit === :mm && x.margin < r.threshold, o1),
                  " pass with less than ", r.threshold, " mm of margin.")
    println(io, "<h2>Counts</h2>")
    rules = [rule for rule in _TALIS_RULE_ORDER if any(x -> x.rule === rule, rep.rows)]
    n(rule, st) = count(x -> x.rule === rule && x.status === st, rep.rows)
    _html_table(io, ["Rule", "Requirement", "Fail", "Review", "Pass"],
                [[rule, talis_requirement(rule, s), _html_num(n(rule, :fail), string(n(rule, :fail))),
                  _html_num(n(rule, :review), string(n(rule, :review))),
                  _html_num(n(rule, :pass), string(n(rule, :pass)))] for rule in rules]; sortable = false)

    println(io, "<h2>Closest approaches</h2>")
    if isempty(r.items)
        println(io, "<p>No 1ω check fails, needs review or passes with less than ", r.threshold, " mm of margin.</p>")
    else
        println(io, "<p class=\"note\">Each drawing is the plane across the 1ω beamlet through the hardware point: ",
                "the 10% contour of the footprint (solid), the worst-case square that hardware must stay outside ",
                "(dashed), the cross-section of the hardware (black) and the closest approach with its length.</p>")
    end
    anchor = Dict{Tuple{String,Symbol,Symbol},Int}()
    for (i, it) in enumerate(r.items)
        row, a = it.row, it.approach
        anchor[(row.beam, row.rule, row.component)] = i
        b = a.budget
        println(io, "<figure id=\"approach-", i, "\"><figcaption>", _xml_escape(row.beam), " / ",
                _xml_escape(row.component), ": ", row.rule, " <span class=\"chip ", row.status, "\">", row.status,
                "</span></figcaption>")
        println(io, _html_svg(omega1_slice_svg(it.slice)).s)
        haskey(images, i) && println(io, "<img alt=\"view of approach ", i, "\" width=\"360\" src=\"",
                                     _xml_escape(images[i]), "\"/>")
        _html_table(io, ["Quantity", "Value"],
                    [["margin of the row (" * string(row.unit) * ")", _html_num(row.margin, _packet_mm(row.margin))],
                     ["transverse distance, nominal (mm)", _html_num(a.nominal, _packet_mm(a.nominal))],
                     ["worst case (mm)", _html_num(a.worst, _packet_mm(a.worst))],
                     ["required (mm)", _html_num(a.required, _packet_mm(a.required))],
                     ["1ω margin (mm)", _html_num(a.margin, _packet_mm(a.margin))],
                     ["pointing (mm)", _html_num(b.pointing, _packet_mm(b.pointing))],
                     ["target offset (mm)", _html_num(b.target_offset, _packet_mm(b.target_offset))],
                     ["hardware tolerance (mm)", _html_num(b.hardware_tolerance, _packet_mm(b.hardware_tolerance))],
                     ["rule tolerance at the point (mm)", _html_num(b.rotation_tol, _packet_mm(b.rotation_tol))],
                     ["focus error, effect (mm)", _html_num(b.focus_delta, _packet_mm(b.focus_delta))],
                     ["beamlet, focal plane (mm)", string(a.beamlet_index, ", ", @sprintf("%+g", a.focus_shift))],
                     ["hardware point (mm)", _packet_vec(a.x_hw)], ["point of the footprint (mm)", _packet_vec(a.x_beam)],
                     ["note", row.note]]; sortable = false)
        println(io, "</figure>")
    end

    if !isempty(rep.reflections)
        println(io, "<h2>Specular reflections</h2>")
        println(io, "<p class=\"note\">One record per 1ω management row decided on reflection: the reflection that ",
                "comes closest to an open port, the area of the hardware in the worst-case footprint, and the ",
                "fraction of the 1ω power of the beam reflected within the required angle of an open port.</p>")
        _html_table(io, ["Beam", "Component", "Angle to port (°)", "Spread (°)", "Worst (°)", "Margin (°)",
                         "Exposed area (mm²)", "Offending area (mm²)", "Reflected fraction", "Point (mm)"],
                    [[x.beam, x.component, _html_num(x.angle; digits = 2), _html_num(x.spread; digits = 2),
                      _html_num(x.worst; digits = 2), _html_num(x.margin; digits = 2),
                      _html_num(x.exposed_area, @sprintf("%.4g", x.exposed_area)),
                      _html_num(x.offending_area, @sprintf("%.4g", x.offending_area)),
                      _html_num(x.reflected_fraction, @sprintf("%.3g", x.reflected_fraction)), _packet_vec(x.point)]
                     for x in rep.reflections])
    end

    if r.disposition !== nothing
        println(io, "<h2>Disposition of the hardware</h2>")
        println(io, "<p class=\"note\">The area (mm²) of each checked component under each 1ω verdict, with its ",
                "most severe verdict and its smallest worst-case margin.</p>")
        _html_table(io, ["Component", "Worst", "Margin (mm)", "Triangles",
                         (replace(string(st), "_" => " ") for st in _DISPOSITION_ORDER)...],
                    [[c.component, _html_chip(c.worst, replace(string(c.worst), "_" => " "); rank = _DISPOSITION_RANK[c.worst]),
                      _html_num(c.margin, _packet_mm(c.margin)), _html_num(c.triangles, string(c.triangles)),
                      (_html_num(c[st], @sprintf("%.4g", c[st])) for st in _DISPOSITION_ORDER)...]
                     for c in r.components])
    end

    println(io, "<h2>All checks</h2>")
    function beam_cell(x)
        i = get(anchor, (x.beam, x.rule, x.component), 0)
        return i == 0 ? x.beam : _HTML(string("<a href=\"#approach-", i, "\">", _xml_escape(x.beam), "</a>"))
    end
    _html_table(io, ["Beam", "Rule", "Component", "Status", "Margin", "Unit", "Worst", "Nominal", "Required", "Note"],
                [[beam_cell(x), x.rule, x.component, _html_chip(x.status; rank = _talis_rank(x)),
                  _html_num(x.margin, _packet_mm(x.margin)), x.unit, _html_num(x.worst, _packet_mm(x.worst)),
                  _html_num(x.nominal, _packet_mm(x.nominal)), _html_num(x.required, _packet_mm(x.required)), x.note]
                 for x in sort(rep.rows; by = x -> (_talis_rank(x), x.margin))]; id = "checks")

    println(io, "<h2>Settings and scope</h2>")
    println(io, "<details><summary>Settings of the check</summary>")
    _html_table(io, ["Setting", "Value"], [[k, _packet_setting(v)] for (k, v) in pairs(s)]; sortable = false)
    println(io, "</details>")
    println(io, "<p>This is a partial assessment. Not evaluated:</p><ul>")
    for x in rep.notes
        println(io, "<li>", _xml_escape(x), "</li>")
    end
    for x in TALIS_OUT_OF_SCOPE
        println(io, "<li>out of scope: ", _xml_escape(x), "</li>")
    end
    println(io, "</ul>")
    return _html_document(title, String(take!(io)); lead, generated)
end

_json_number(x::Real) = isfinite(x) ? x : string(x)
_json_vec(v) = [_json_number(x) for x in v]

function _talis_data(rep::TaLISReport)
    rows = [Dict("beam" => x.beam, "rule" => x.rule, "component" => x.component, "status" => x.status,
                 "margin" => _json_number(x.margin), "unit" => x.unit, "worst" => _json_number(x.worst),
                 "nominal" => _json_number(x.nominal), "required" => _json_number(x.required),
                 "point_mm" => _json_vec(x.point), "note" => x.note) for x in rep.rows]
    settings = Dict(string(k) => (v isa AbstractVector || v isa Tuple ? collect(Any, map(y -> y isa Real ? _json_number(y) : string(y), v)) :
                                  v isa Real ? _json_number(v) : string(v)) for (k, v) in pairs(rep.settings))
    counts = Dict(string(st) => count(x -> x.status === st, rep.rows) for st in (:fail, :review, :pass))
    return Dict("source" => TALIS_SOURCE, "settings" => settings, "counts" => counts, "rows" => rows,
                "not_evaluated" => rep.notes, "out_of_scope" => collect(TALIS_OUT_OF_SCOPE),
                "reflections" => [Dict("beam" => x.beam, "component" => x.component, "angle_deg" => _json_number(x.angle),
                                       "spread_deg" => _json_number(x.spread), "worst_deg" => _json_number(x.worst),
                                       "margin_deg" => _json_number(x.margin), "exposed_area_mm2" => _json_number(x.exposed_area),
                                       "offending_area_mm2" => _json_number(x.offending_area),
                                       "reflected_fraction" => _json_number(x.reflected_fraction),
                                       "point_mm" => _json_vec(x.point), "port" => _json_vec(x.port))
                                  for x in rep.reflections])
end

_write_json(path::AbstractString, data) = (open(io -> JSON.print(io, data, 2), path, "w"); path)

"""
    write_talis_json(path, report) -> path

Write a TaLIS report as JSON: the source of the rules, the `settings`, the counts by
status, every row (`beam`, `rule`, `component`, `status`, `margin`, `unit`, `worst`,
`nominal`, `required`, `point_mm`, `note`), the reflection records, and the checks that
were not evaluated or are out of scope. A value that is not finite is written as text.
"""
write_talis_json(path::AbstractString, rep::TaLISReport) = _write_json(path, _talis_data(rep))

function _review_data(r::Omega1Review)
    d = _talis_data(r.report)
    d["threshold_mm"] = r.threshold
    d["approaches"] = [begin
                           a = it.approach
                           Dict("beam" => a.beam, "component" => a.component, "rule" => it.row.rule,
                                "status" => it.row.status, "margin_mm" => _json_number(a.margin),
                                "nominal_mm" => _json_number(a.nominal), "worst_mm" => _json_number(a.worst),
                                "required_mm" => _json_number(a.required), "beamlet" => a.beamlet_index,
                                "focus_shift_mm" => a.focus_shift, "s_star_mm" => _json_number(a.s_star),
                                "hardware_point_mm" => _json_vec(a.x_hw), "footprint_point_mm" => _json_vec(a.x_beam),
                                "budget_mm" => Dict(string(k) => _json_number(v) for (k, v) in pairs(a.budget)),
                                "slice" => @sprintf("slice-%03d.svg", i))
                       end for (i, it) in enumerate(r.items)]
    d["disposition"] = [Dict("component" => c.component, "triangles" => c.triangles, "worst" => c.worst,
                             "margin_mm" => _json_number(c.margin),
                             "area_mm2" => Dict(string(st) => c[st] for st in _DISPOSITION_ORDER)) for c in r.components]
    return d
end

"""
    write_omega1_review(dir, review; title = "1ω review", generated = nothing, json = false,
                        images = Dict()) -> Vector{String}

Write the 1ω review packet into the folder `dir` (made if missing) and return the paths:
`review-1w.html` ([`omega1_review_html`](@ref)), `talis.csv` ([`write_talis_csv`](@ref)),
one `slice-NNN.svg` per approach of the review in the order of `review.items`, and with
`json = true` the file `review-1w.json` with the rows, the approaches with their margin
terms, and the disposition per component.
"""
function write_omega1_review(dir::AbstractString, r::Omega1Review; title::AbstractString = "1ω review",
                             generated = nothing, json::Bool = false, images::AbstractDict = Dict{Int,String}())
    mkpath(dir)
    paths = String[joinpath(dir, "review-1w.html")]
    write(paths[1], omega1_review_html(r; title, generated, images))
    push!(paths, write_talis_csv(joinpath(dir, "talis.csv"), r.report))
    for (i, it) in enumerate(r.items)
        p = joinpath(dir, @sprintf("slice-%03d.svg", i))
        write(p, omega1_slice_svg(it.slice))
        push!(paths, p)
    end
    json && push!(paths, _write_json(joinpath(dir, "review-1w.json"), _review_data(r)))
    return paths
end

# ---------------------------------------------------------------- rules of engagement

_roe_tier_rank(t::Symbol) = t === :red ? 0 : t === :amber ? 1 : 2

"""
    roe_html(report, beams = nothing; title = "Beam rules of engagement", notes = String[],
             generated = nothing) -> String

The rules-of-engagement sheet of a report ([`ROEReport`](@ref)) as one self-contained
HTML page: the criteria with their thresholds and nominal values, the chamber port map
colored by tier ([`port_map_svg`](@ref), drawn when the `beams` of the scan are given),
and for each role the multiple-outage budget and a table of the units from critical to
auto-drop with the binding criterion, its score against the thresholds, the scores on
every criterion, the companion drop and the pair caveats. The tables sort on a click of a
column head. `notes` are lines shown under the criteria, such as how the criteria were
chosen. `generated` is a line for the foot of the page. Without it the page of the same
report is the same text. The page refers to no other file and to nothing on the network.
"""
function roe_html(r::ROEReport, beams::Union{Nothing,AbstractVector{<:Beam}} = nothing;
                  title::AbstractString = "Beam rules of engagement", notes::AbstractVector{<:AbstractString} = String[],
                  generated = nothing)
    io = IOBuffer()
    nt(t) = count(==(t), r.tiers)
    lead = string(length(r.units), " units of ", length(r.beams), " beams against ", length(r.criteria),
                  " criteria: ", nt(:red), " critical, ", nt(:amber), " ask RI, ", nt(:green), " auto-drop. ",
                  "An outage is a beam or quad at zero power with every other beam as planned.")
    println(io, "<h2>Criteria</h2>")
    _html_table(io, ["Criterion", "Case", "Measure", "Green ≤", "Amber ≤", "Unit", "Nominal", "Path"],
                [begin
                     j = c.case === :default ? 1 : findfirst(==(c.case), r.cases)
                     [c.name, r.cases[j], c.measure, _html_num(c.green, _roe_number(c.green)),
                      _html_num(c.amber, _roe_number(c.amber)), c.label, _roe_number(v),
                      r.linear[j] ? "linear deposition" : "full evaluation"]
                 end for (c, v) in zip(r.criteria, r.baseline)]; sortable = false)
    println(io, "<p class=\"note\">Critical units hold the shot, Ask RI units need approval before the drop, and ",
            "Auto-Drop units may drop without a request. A unit that is green alone is raised to Ask RI when a pair ",
            "with another green unit is critical.</p>")
    for x in notes
        println(io, "<p class=\"note\">", _xml_escape(x), "</p>")
    end
    if beams !== nothing
        println(io, "<h2>Chamber port map</h2>")
        for kind in (:quad, :beam)
            any(u -> u.kind === kind, r.units) || continue
            println(io, "<figure id=\"ports-", kind, "\"><figcaption>", kind === :quad ? "Quads" : "Beams",
                    " by tier</figcaption>")
            println(io, _html_svg(port_map_svg(r, beams; kind)).s)
            println(io, "</figure>")
        end
    end
    rows = roe_table(r)
    for role in unique(row.role for row in rows)
        println(io, "<h2>Role: ", _xml_escape(role), "</h2>")
        bs = [b for b in r.budgets if b.role === role]
        if !isempty(bs)
            println(io, "<ul>")
            for b in bs
                print(io, "<li>Budget (", b.kind, "s): ", _xml_escape(b.text), ".")
                isempty(b.exceptions) || print(io, " Pairs that leave green: ",
                                               _xml_escape(join((a * " + " * c for (a, c) in b.exceptions), ", ")), ".")
                println(io, "</li>")
            end
            println(io, "</ul>")
        end
        _html_table(io, ["Unit", "Kind", "Decision", "Binding", "Score", "Green ≤", "Amber ≤",
                         (string(c.name) for c in r.criteria)..., "Companion", "Pair caveats"],
                    [[row.unit, row.kind,
                      _html_chip(row.tier, row.decision * (row.tier === row.single ? "" : " (alone: $(row.single))");
                                 rank = _roe_tier_rank(row.tier)),
                      row.binding, _html_num(row.score, _roe_number(row.score)),
                      _html_num(row.green, _roe_number(row.green)), _html_num(row.amber, _roe_number(row.amber)),
                      (_html_num(x, _roe_number(x)) for x in values(row.scores))...,
                      row.companion === nothing ? "" : string(row.companion, ": ", row.companion_tier),
                      join(row.caveats, "; ")] for row in rows if row.role === role])
    end
    return _html_document(title, String(take!(io)); lead, generated)
end

function _roe_data(r::ROEReport, notes)
    crit = [begin
                j = c.case === :default ? 1 : findfirst(==(c.case), r.cases)
                Dict("name" => c.name, "case" => r.cases[j], "measure" => c.measure, "green" => _json_number(c.green),
                     "amber" => _json_number(c.amber), "unit" => c.label,
                     "nominal" => v isa Real ? _json_number(v) : _json_vec(v), "linear" => r.linear[j])
            end for (c, v) in zip(r.criteria, r.baseline)]
    units = [Dict("role" => row.role, "unit" => row.unit, "kind" => row.kind, "beams" => row.beams, "tier" => row.tier,
                  "decision" => row.decision, "single" => row.single, "binding" => row.binding,
                  "score" => _json_number(row.score), "green" => _json_number(row.green),
                  "amber" => _json_number(row.amber),
                  "scores" => Dict(string(k) => _json_number(v) for (k, v) in pairs(row.scores)),
                  "companion" => row.companion, "companion_tier" => row.companion_tier, "caveats" => row.caveats)
             for row in roe_table(r)]
    budgets = [Dict("role" => b.role, "kind" => b.kind, "green" => b.green, "pairs" => b.pairs, "safe" => b.safe,
                    "limit" => b.limit, "exceptions" => [[a, c] for (a, c) in b.exceptions], "text" => b.text)
               for b in r.budgets]
    counts = Dict(string(t) => count(==(t), r.tiers) for t in _ROE_TIERS)
    return Dict("criteria" => crit, "notes" => collect(notes), "counts" => counts, "units" => units, "budgets" => budgets)
end

"""
    write_roe_json(path, report; notes = String[]) -> path

Write a rules-of-engagement report as JSON: the criteria with thresholds and nominal
values, the counts by tier, the rows of [`roe_table`](@ref) as `units`, the budgets, and
`notes`.
"""
write_roe_json(path::AbstractString, r::ROEReport; notes::AbstractVector{<:AbstractString} = String[]) =
    _write_json(path, _roe_data(r, notes))

"""
    write_roe_packet(dir, report, beams = nothing; title = "Beam rules of engagement",
                     notes = String[], generated = nothing, json = false) -> Vector{String}

Write the rules-of-engagement packet into the folder `dir` (made if missing) and return
the paths: `roe.html` ([`roe_html`](@ref)), `roe.csv` ([`write_roe_csv`](@ref)), `roe.md`
([`roe_markdown`](@ref)), with the `beams` of the scan `ports.svg`
([`port_map_svg`](@ref)), and with `json = true` the file `roe.json`
([`write_roe_json`](@ref)).
"""
function write_roe_packet(dir::AbstractString, r::ROEReport, beams::Union{Nothing,AbstractVector{<:Beam}} = nothing;
                          title::AbstractString = "Beam rules of engagement",
                          notes::AbstractVector{<:AbstractString} = String[], generated = nothing, json::Bool = false)
    mkpath(dir)
    paths = String[joinpath(dir, "roe.html")]
    write(paths[1], roe_html(r, beams; title, notes, generated))
    push!(paths, write_roe_csv(joinpath(dir, "roe.csv"), r))
    push!(paths, write_roe_markdown(joinpath(dir, "roe.md"), r; title))
    beams === nothing || push!(paths, write_port_map_svg(joinpath(dir, "ports.svg"), r, beams))
    json && push!(paths, write_roe_json(joinpath(dir, "roe.json"), r; notes))
    return paths
end

"""
    roe_default_criteria(scenario, beams; foils = (), share = 0.05) -> (; criteria, notes)

A starting set of criteria for [`roe_scan`](@ref), chosen from what the design contains,
with one line of `notes` per choice:

- observed elements that close on themselves, such as a capsule: [`roe_drive`](@ref),
  [`roe_mode`](@ref) for P2 and P4, and [`roe_hotspot_velocity`](@ref);
- any other observed elements, such as a package face: [`roe_face_flux`](@ref),
  [`roe_face_nonuniformity`](@ref) and [`roe_temperature`](@ref);
- each part named in `foils`: [`roe_foil_irradiance`](@ref),
  [`roe_foil_nonuniformity`](@ref) and [`roe_foil_centroid`](@ref), with the part in the
  criterion names;
- with nothing observed and no foil, each part that absorbs at least `share` of the
  deposited laser power, up to four: [`roe_foil_irradiance`](@ref).

The criteria of a design whose beams land on nothing are empty. Observed elements close on
themselves when their area-weighted normals sum to less than 5% of their area. The set
and every threshold are engineering defaults and not facility criteria: a shot sheet
needs the criteria of the experiment.
"""
function roe_default_criteria(s::Scenario, beams::AbstractVector{<:Beam}; foils = (), share::Real = 0.05)
    t = s.target
    criteria, notes = ROECriterion[], String[]
    obs = s.observed
    if !isempty(obs) && s.model !== nothing
        A = sum(t.elem_area[e] for e in obs)
        closed = norm(sum(t.elem_area[e] * t.elem_normal[e] for e in obs)) < 0.05 * A
        names = join(unique(t.components[t.elem_comp[e]].part for e in obs), ", ")
        if closed
            append!(criteria, [roe_drive(), roe_mode(2), roe_mode(4), roe_hotspot_velocity()])
            push!(notes, "The observed elements ($names) close on themselves: drive, P2, P4 and hot-spot velocity.")
        else
            append!(criteria, [roe_face_flux(), roe_face_nonuniformity(), roe_temperature()])
            push!(notes, "The observed elements ($names) form an open face: mean flux, nonuniformity and temperature.")
        end
    end
    fs = foils isa Symbol || foils isa AbstractString ? (Symbol(foils),) : Tuple(Symbol(f) for f in foils)
    for f in fs
        append!(criteria, [roe_foil_irradiance(f; name = Symbol(f, :_irradiance)),
                           roe_foil_nonuniformity(f; name = Symbol(f, :_irradiance_rms)),
                           roe_foil_centroid(f; name = Symbol(f, :_centroid))])
        push!(notes, "Foil $f: absorbed intensity, its nonuniformity and the shift of its centroid.")
    end
    if isempty(criteria)
        dep = evaluate(s, beams).deposition.element
        total = sum(dep)
        lit = Tuple{Float64,Symbol}[]
        for p in parts(t)
            es = [e for e in elements(t, p) if t.elem_side[e] == 1]
            isempty(es) && continue
            P = sum(dep[e] for e in es)
            total > 0 && P >= share * total && push!(lit, (P / total, p))
        end
        sort!(lit; rev = true)
        for (f, p) in first(lit, 4)
            push!(criteria, roe_foil_irradiance(p; name = Symbol(p, :_irradiance)))
            push!(notes, @sprintf("Part %s absorbs %.0f%% of the deposited power: absorbed intensity.", p, 100f))
        end
        isempty(criteria) && push!(notes, "The beams deposit power on no part: no criteria.")
    end
    push!(notes, "The criteria and thresholds are engineering defaults chosen from the design, not facility criteria.")
    return (; criteria, notes)
end
