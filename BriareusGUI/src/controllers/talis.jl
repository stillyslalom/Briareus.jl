# The Review step's TaLIS checks: the NIF target-configuration rules of `talis_check`
# (3ω capture, incidence and clearance, 1ω clearance and management, the TAS stay-in
# zone) on chosen capture and hardware components, with the clearance check's error stack.
# One `talis_check` call runs on a worker (it takes well under a second on the GUI's
# meshes); canceling drops its result.

"""
Rule tolerance pivots of the TaLIS checks (`talis_check`'s `pivot`): the 0.5° rule
tolerance becomes a lateral allowance `r·tan(tolerance)` with `r` the distance of the
hardware point from TCC (`:tcc`), from the beam's aim point (`:aim`), the larger of
the two (`:max`, most conservative), or from a given point (`:point`, e.g. the mount
of a component).
"""
const TALIS_PIVOTS = (:tcc, :aim, :max, :point)

"Short names of the TaLIS rules (`TaLISRow.rule`), in report order."
const TALIS_RULE_LABELS = (capture_3w = "3ω capture", incidence_3w = "3ω incidence",
                           clearance_3w = "3ω clearance", clearance_1w = "1ω clearance",
                           management_1w = "1ω management", diagnostic_1w = "1ω diagnostics",
                           tas_stay_in = "TAS stay-in zone")

"""
What the rule tolerance stands for, carried by every TaLIS report: the requirements'
"up to 0.5 degrees" of metrology, alignment and positioning tolerance is the target
component rotated off nominal (insertion or assembly error), in the least favorable
direction. It applies to the clearances and the capture margin; the 70° incidence limit
does not include it.
"""
const TALIS_TOLERANCE_NOTE =
    "The tolerance is the component rotated off nominal (target inserted rotated, or " *
    "assembled off-square), in the least favorable direction; it applies to the " *
    "clearances and the capture margin, not to the 70° incidence limit. For a " *
    "component assembled off-square, " *
    "set the pivot to its mount."

"""
    TaLISState()

State of the TaLIS checks. Settings: `capture` (components beams are
meant to land on: backlighter foils, capture surfaces, a hohlraum wall), `hardware`
(other components checked for 3ω and 1ω clearance; capture components are checked
too, for beams landing elsewhere), `diagnostics` (diagnostic hardware: checked for 3ω
clearance and kept 10–25 mm from 1ω light, `talis_check`'s `diagnostics`),
`tolerance` (rule tolerance, degrees), `pivot` (one of [`TALIS_PIVOTS`](@ref)) and
`pivot_point` (mm, for `:point`), and `omega1` (run the 1ω rules). The error stack and
the hardware sampling come from the [`ReviewState`](@ref). Job: `running`, `started`,
`status`, `error`, `report` (a core `TaLISReport`) and `inputs` (the target, beams and
`talis_check` keywords it checked). With the 1ω rules on, the check also leaves
`disposition`, the 1ω verdict of every mesh triangle (a core `Omega1Disposition`), which
the 3D view colors by ([`DISPOSITION_MODES`](@ref)); `tight` (mm) is the margin up to
which hardware that clears is drawn as tight, and `hide_curved` draws hardware that
scatters by its curvature as clear.
"""
struct TaLISState
    capture::Observable{Vector{Symbol}}
    hardware::Observable{Vector{Symbol}}
    diagnostics::Observable{Vector{Symbol}}
    tolerance::Observable{Float64}
    pivot::Observable{Symbol}
    pivot_point::Observable{NTuple{3,Float64}}
    omega1::Observable{Bool}
    running::Observable{Bool}
    started::Observable{Float64}
    status::Observable{String}
    error::Observable{String}
    report::Observable{Union{Nothing,TaLISReport}}
    disposition::Observable{Union{Nothing,Omega1Disposition}}
    tight::Observable{Float64}
    hide_curved::Observable{Bool}
    inputs::Base.RefValue{Any}
    known::Base.RefValue{Vector{Symbol}}
    generation::Base.RefValue{Int}
end

TaLISState() = TaLISState(Observable(Symbol[]), Observable(Symbol[]), Observable(Symbol[]),
                          Observable(0.5),
                          Observable(:tcc), Observable((0.0, 0.0, 0.0)), Observable(true),
                          Observable(false), Observable(0.0), Observable(""), Observable(""),
                          Observable{Union{Nothing,TaLISReport}}(nothing),
                          Observable{Union{Nothing,Omega1Disposition}}(nothing), Observable(1.0),
                          Observable(false), Ref{Any}(nothing),
                          Ref(Symbol[]), Ref(0))

"""
    default_selection!(tl::TaLISState, t::Target; capture = (), observed = ())

Choose the capture components and hardware for a target's components. The surfaces `capture`
capture the beams; they are names, and a session passes the holders of the wall and foil
roles, such as the wall of a hohlraum and every backlighter foil. The observed surface
(`observed`, such as a capsule) is left out as an engineering default: light on it is drive
physics, handled by keep-outs and the drive, and is not a TaLIS hardware question.
Everything else (LEH caps and lips, closed windows, the main target) is checked. Checked LEH
caps and lips go to `talis_check` as its `leh`, so their 3ω clearance is an LEH-clearance
review; it does not count as a failure. No component is a diagnostic by default. Components
the selection was already made for keep their roles, and a selection made afresh clears the
diagnostics too.
"""
function default_selection!(tl::TaLISState, t::Target; capture = (), observed = ())
    names = surface_names(t)
    fresh = !any(in(names), tl.capture[])
    _update_roles!(tl.capture, tl.hardware, tl.known, names,
                   n -> n in capture ? 1 : n in observed ? 0 : 2)
    _assign!(tl.diagnostics,
             fresh ? Symbol[] : filter(n -> n in names && !(n in tl.capture[]) && !(n in tl.hardware[]),
                                       tl.diagnostics[]))
    return tl
end

# Components that form an entrance-hole edge: a hohlraum part's caps and lips
# (`<part>_cap_top`, …).
_leh_edge(n::Symbol) = occursin(r"_(cap_top|cap_bottom|leh_top_lip|leh_bottom_lip)$", String(n))

"""
    set_talis_role!(tl::TaLISState, component::Symbol, role::Symbol)

Make a component a `:capture` surface, checked `:hardware`, `:diagnostic` hardware
(the 1ω diagnostic clearance instead of the 0.2 mm one), or none of them (`:none`).
"""
function set_talis_role!(tl::TaLISState, name::Symbol, role::Symbol)
    role in (:capture, :hardware, :diagnostic, :none) ||
        throw(ArgumentError("role must be :capture, :hardware, :diagnostic or :none"))
    cp = filter(!=(name), tl.capture[])
    hw = filter(!=(name), tl.hardware[])
    dg = filter(!=(name), tl.diagnostics[])
    role === :capture && push!(cp, name)
    role === :hardware && push!(hw, name)
    role === :diagnostic && push!(dg, name)
    _assign!(tl.capture, cp)
    _assign!(tl.hardware, hw)
    _assign!(tl.diagnostics, dg)
    return tl
end

"""
    talis_role(tl::TaLISState, component::Symbol) -> Symbol

`:capture`, `:hardware`, `:diagnostic` or `:none` (see [`set_talis_role!`](@ref)).
"""
talis_role(tl::TaLISState, name::Symbol) =
    name in tl.capture[] ? :capture : name in tl.hardware[] ? :hardware :
    name in tl.diagnostics[] ? :diagnostic : :none

"""
    set_talis_pivot!(tl::TaLISState, pivot::Symbol)

One of [`TALIS_PIVOTS`](@ref).
"""
function set_talis_pivot!(tl::TaLISState, pivot::Symbol)
    pivot in TALIS_PIVOTS || throw(ArgumentError("pivot must be one of $TALIS_PIVOTS"))
    _assign!(tl.pivot, pivot)
    return tl
end

"""
    edit_talis!(tl::TaLISState, key::Symbol, text) -> Bool

Edit `:tolerance` (degrees, ≥ 0) or `:pivot_point` (`"x, y, z"` in mm) from text. A
rejected value sets `error` and returns `false`.
"""
function edit_talis!(tl::TaLISState, key::Symbol, text)
    try
        if key === :tolerance
            v = _parse_number(text, "the rule tolerance"; nonnegative = true)
            v < 90 || throw(ArgumentError("the rule tolerance is an angle below 90°"))
            _assign!(tl.tolerance, v)
        elseif key === :pivot_point
            parts = split(strip(string(text)), r"[\s,;]+"; keepempty = false)
            length(parts) == 3 ||
                throw(ArgumentError("the pivot point is three coordinates x, y, z (mm), got \"$text\""))
            _assign!(tl.pivot_point, Tuple(_parse_number(p, "a pivot coordinate") for p in parts))
        else
            throw(ArgumentError("unknown TaLIS setting :$key"))
        end
        tl.error[] = ""
        return true
    catch err
        err isa ArgumentError || rethrow()
        tl.error[] = _errmsg(err)
        return false
    end
end

"""
    talis_kwargs(tl::TaLISState, rv::ReviewState) -> NamedTuple

The keyword arguments of the `talis_check` call: capture and hardware components, the
checked LEH caps and lips (`leh`), the diagnostics (3ω-checked as hardware, and kept
from 1ω light by the diagnostic rule instead of the 0.2 mm 1ω clearance), the rule
tolerance and its pivot, the clearance check's error stack (pointing, positioning,
hardware tolerance, focus) and hardware sampling level.
"""
function talis_kwargs(tl::TaLISState, rv::ReviewState)
    e = rv.errors[]
    checked = Tuple(unique([tl.capture[]; tl.hardware[]]))
    pivot = tl.pivot[] === :point ? tl.pivot_point[] : tl.pivot[]
    return (; capture = Tuple(tl.capture[]), hardware = Tuple(unique([checked...; tl.diagnostics[]])),
            leh = Tuple(n for n in tl.hardware[] if _leh_edge(n)), omega1_components = checked,
            diagnostics = Tuple(tl.diagnostics[]),
            omega1 = tl.omega1[], tolerance = tl.tolerance[], pivot,
            pointing = e[:pointing], target_offset = e[:target_offset],
            hardware_tolerance = e[:hardware_tolerance], focus = e[:focus],
            level = REVIEW_QUALITY[rv.quality[]].level)
end

"""
    start_talis!(tl::TaLISState, rv::ReviewState, target, beams; deliver = f -> f(),
                 spawn = true)

Run `talis_check` on `beams` (on a worker with `spawn`), with the settings of `tl`
and the error stack of `rv` (see [`talis_kwargs`](@ref)). The result arrives through
`deliver`. [`cancel_talis!`](@ref) drops it.
"""
function start_talis!(tl::TaLISState, rv::ReviewState, t::Target, beams::AbstractVector;
                      deliver = f -> f(), spawn::Bool = true)
    tl.running[] && return tl
    msg = isempty(beams) ? "no beams are enabled" :
          isempty(tl.capture[]) ? "choose the capture components the beams land on" : nothing
    msg === nothing || (tl.status[] = msg; return tl)
    kw = talis_kwargs(tl, rv)
    bs = collect(beams)
    g = (tl.generation[] += 1)
    tl.inputs[] = (t, bs, kw)
    tl.started[] = time()
    tl.status[] = "checking…"
    tl.running[] = true
    job = function ()
        out = _try_job(() -> (talis_check(t, bs; kw...),
                              kw.omega1 ? omega1_disposition(t, bs; disposition_kwargs(t, kw)...) : nothing))
        deliver(() -> _finish_talis!(tl, out, g))
    end
    _run_job(job, spawn)
    return tl
end

function _finish_talis!(tl::TaLISState, out, g::Integer)
    g == tl.generation[] || return tl                  # canceled or superseded
    if out.err === nothing
        tl.disposition[] = last(out.value)
        tl.report[] = first(out.value)
        tl.status[] = @sprintf("done in %s: ", _seconds(time() - tl.started[])) * talis_brief(tl)
    else
        tl.status[] = "failed: " * _errmsg(out.err)
    end
    tl.running[] = false
    return tl
end

"""
What the 3D view colors hardware by from the 1ω disposition of the TaLIS check:
`:omega1_disposition` (its verdict, in the colors of [`DISPOSITION_COLORS`](@ref)) and
`:omega1_margin` (its worst-case 1ω clearance margin, mm).
"""
const DISPOSITION_MODES = (:omega1_disposition, :omega1_margin)

"""
Colors (hex) of the 1ω verdicts: gray for hardware that clears, bluish green for hardware
in the footprint that scatters by its curvature, amber for hardware that is tight, clears
only by its reflection or needs review, vermilion for hardware that fails, and light gray
for hardware the check did not cover.
"""
const DISPOSITION_COLORS = Dict(:unchecked => "#D9D9D9", :clear => "#9E9E9E", :tight => "#E69F00",
                                :curved_safe => "#009E73", :specular_clear => "#E69F00",
                                :review => "#E69F00", :fail => "#D55E00")

"Largest margin (mm) that the disposition follows, and so the largest tight threshold."
const DISPOSITION_RANGE = 5.0

"""
    disposition_kwargs(target, kw) -> NamedTuple

The keyword arguments of the `omega1_disposition` call that goes with a `talis_check`
call with the keywords `kw` ([`talis_kwargs`](@ref)): the same components, diagnostics,
tolerance and error stack, margins followed to [`DISPOSITION_RANGE`](@ref), the shadow
flag, and one verdict per facet of a target of more than 20,000 triangles (the facets
of CAD are fine already).
"""
disposition_kwargs(t::Target, kw) =
    (; components = kw.omega1_components, diagnostics = kw.diagnostics, tolerance = kw.tolerance,
     pivot = kw.pivot, pointing = kw.pointing, target_offset = kw.target_offset,
     hardware_tolerance = kw.hardware_tolerance, focus = kw.focus, margin_range = DISPOSITION_RANGE,
     shadows = true, level = length(t.tri_area) > 20_000 ? 0 : kw.level)

"""
    set_disposition_filter!(tl::TaLISState; tight = tl.tight[], hide_curved = tl.hide_curved[]) -> tl

Set the margin (mm, 0 to [`DISPOSITION_RANGE`](@ref)) up to which hardware that clears
the 1ω footprints is drawn as tight, and whether hardware that scatters by its curvature
is drawn as clear, so that a review shows only tight and failing hardware. Neither
changes the check.
"""
function set_disposition_filter!(tl::TaLISState; tight::Real = tl.tight[], hide_curved::Bool = tl.hide_curved[])
    0 <= tight <= DISPOSITION_RANGE ||
        throw(ArgumentError("the tight threshold is between 0 and $DISPOSITION_RANGE mm"))
    _assign!(tl.tight, Float64(tight))
    _assign!(tl.hide_curved, hide_curved)
    return tl
end

"""
    disposition_of(tl::TaLISState, target) -> Union{Nothing,Omega1Disposition}

The 1ω disposition of the last check if it was computed for `target`, else `nothing`
(no check yet, the 1ω rules off, or the target changed since).
"""
function disposition_of(tl::TaLISState, t::Target)
    d = tl.disposition[]
    (d === nothing || tl.inputs[] === nothing || first(tl.inputs[]) !== t) && return nothing
    return d
end

"""
    disposition_status(tl::TaLISState, target) -> Union{Nothing,Vector{Symbol}}

The 1ω verdict of each mesh triangle of `target` as the view draws it: the verdicts of
the check with the tight threshold and the curved-safe filter of `tl` applied
([`set_disposition_filter!`](@ref)); `nothing` without a disposition for `target`.
"""
function disposition_status(tl::TaLISState, t::Target)
    d = disposition_of(tl, t)
    d === nothing && return nothing
    thr, hide = tl.tight[], tl.hide_curved[]
    return map(d.status, d.margin) do st, m
        st in (:clear, :tight) ? (m <= thr ? :tight : :clear) :
        (st === :curved_safe && hide) ? :clear : st
    end
end

"""
    disposition_info(tl::TaLISState, target, beams, triangle) -> String

What the 1ω check found at mesh triangle `triangle` of `target`: its verdict, the
worst-case margin with the beam and beamlet that set it, the radii of curvature it was
judged on, the open port its reflection comes nearest with the angle, and whether other
hardware shadows it. Empty without a disposition for `target`.
"""
function disposition_info(tl::TaLISState, t::Target, beams::AbstractVector, k::Integer)
    d = disposition_of(tl, t)
    (d === nothing || !checkbounds(Bool, d.status, k)) && return ""
    st = d.status[k]
    st === :unchecked && return "not covered by the 1ω check"
    io = IOBuffer()
    print(io, replace(String(st), "_" => " "))
    if isfinite(d.margin[k])
        b = d.beam[k]
        @printf(io, "; 1ω margin %+.3f mm (%s, beamlet %d)", d.margin[k],
                checkbounds(Bool, beams, b) ? beams[b].name : "beam $b", d.beamlet[k])
    else
        @printf(io, "; 1ω margin above %g mm", d.settings.margin_range)
    end
    rad(r) = isfinite(r) ? @sprintf("%.3g mm", r) : "flat"
    isnan(d.R1[k]) || print(io, "; radii of curvature ", rad(d.R1[k]), ", ", rad(d.R2[k]))
    d.port[k] == 0 ||
        @printf(io, "; reflection %.1f° from open port %d (%.0f° required)", d.port_angle[k], d.port[k],
                d.settings.reflection_angle)
    d.shadowed[k] && print(io, "; shadowed by hardware upstream")
    return String(take!(io))
end

"""
    cancel_talis!(tl::TaLISState)

Drop the check in progress. `talis_check` does not take a progress callback, so the
worker runs to the end and the result is discarded. The previous report stays.
"""
function cancel_talis!(tl::TaLISState)
    tl.running[] || return tl
    tl.generation[] += 1
    tl.running[] = false
    tl.status[] = "canceled"
    return tl
end

"""
    talis_rows(tl::TaLISState) -> Vector{TaLISRow}

Every row of the report: failures first, then review, then pass, each by margin.
"""
talis_rows(tl::TaLISState) =
    tl.report[] === nothing ? TaLISRow[] : violations(tl.report[]; status = (:fail, :review, :pass))

"""
    talis_rule_rows(tl::TaLISState) -> Vector{NamedTuple}

One row per rule in the report: `(; rule, label, requirement, pass, review, fail,
worst)` with `worst` the row of smallest margin.
"""
function talis_rule_rows(tl::TaLISState)
    r = tl.report[]
    r === nothing && return NamedTuple[]
    out = NamedTuple[]
    for rule in keys(TALIS_RULE_LABELS)
        rs = [x for x in r.rows if x.rule === rule]
        isempty(rs) && continue
        push!(out, (; rule, label = TALIS_RULE_LABELS[rule],
                    requirement = talis_requirement(rule, r.settings),
                    pass = count(x -> x.status === :pass, rs),
                    review = count(x -> x.status === :review, rs),
                    fail = count(x -> x.status === :fail, rs),
                    worst = rs[argmin([x.margin for x in rs])]))
    end
    return out
end

"""
The requirements the TaLIS checks leave out, in a line for the report's notes
(`Briareus.TALIS_OUT_OF_SCOPE` has each one in full, as the CSV export writes them).
"""
const TALIS_SCOPE_NOTE = "Out of scope: positioner clearance, TAS clearance and withdrawal, " *
                         "shroud swept volume, ARC beam footprints."

"""
    talis_notes(tl::TaLISState) -> Vector{NamedTuple}

What the report does not settle, as `(; kind, text)` rows, each a short line: the
rule tolerance as applied (`:tolerance`, always first, also before a check), the error
set (`:settings`), the checks the report did not evaluate (`:skipped`), the
requirements out of scope in one line (`:scope`, [`TALIS_SCOPE_NOTE`](@ref)), and the
source (`:source`). What the tolerance means is `TALIS_TOLERANCE_NOTE`, in the window's
Help.
"""
function talis_notes(tl::TaLISState)
    r = tl.report[]
    notes = NamedTuple[]
    if r === nothing
        piv = tl.pivot[]
        txt = piv === :tcc ? "r = distance from TCC" : piv === :aim ? "r = distance from the beam's aim point" :
              piv === :max ? "r = larger of the distances from TCC and from the aim point" :
              @sprintf("r = distance from (%g, %g, %g) mm", tl.pivot_point[]...)
        push!(notes, (; kind = :tolerance,
                      text = @sprintf("Rule tolerance %g° as a lateral allowance r·tan(%g°), %s.",
                                      tl.tolerance[], tl.tolerance[], txt)))
    else
        s = r.settings
        push!(notes, (; kind = :tolerance,
                      text = @sprintf("Rule tolerance %g° as a lateral allowance r·tan(%g°), %s.",
                                      s.tolerance, s.tolerance, s.pivot_text)))
        push!(notes, (; kind = :settings,
                      text = @sprintf("Errors, added linearly: transverse %g mm (pointing %g + positioning %g + hardware %g), focus ±%g mm.",
                                      s.transverse, s.pointing, s.target_offset,
                                      s.hardware_tolerance, s.focus)))
        for n in r.notes
            push!(notes, (; kind = :skipped, text = "Not evaluated: " * n))
        end
    end
    push!(notes, (; kind = :scope, text = TALIS_SCOPE_NOTE))
    push!(notes, (; kind = :source, text = "Rules from " * TALIS_SOURCE))
    return notes
end

"""
    talis_brief(tl::TaLISState) -> String

Fail/review counts and the number of checks, in a line.
"""
function talis_brief(tl::TaLISState)
    r = tl.report[]
    r === nothing && return "not run"
    isempty(r.rows) && return "no checks"
    nf = count(x -> x.status === :fail, r.rows)
    nr = count(x -> x.status === :review, r.rows)
    nb = length(unique(x.beam for x in r.rows if x.rule !== :tas_stay_in))
    return @sprintf("%d beams, %d checks: %d fail, %d review", nb, length(r.rows), nf, nr)
end

"""
    export_talis(tl::TaLISState, path) -> path

Write the report with `write_talis_csv`.
"""
function export_talis(tl::TaLISState, path::AbstractString)
    tl.report[] === nothing && throw(ArgumentError("run the TaLIS checks first"))
    write_talis_csv(path, tl.report[])
    tl.status[] = "wrote $(basename(path))"
    return path
end

"""
    talis_stale(tl::TaLISState, rv::ReviewState, target, beams) -> Bool

Whether anything the last check depends on changed since: the target, the enabled
beams, or any setting (the `talis_kwargs` of `tl` and `rv`: component roles, rule
tolerance and pivot, 1ω rules, error stack, sampling).
"""
function talis_stale(tl::TaLISState, rv::ReviewState, t::Target, beams::AbstractVector)
    inp = tl.inputs[]
    inp === nothing && return false
    return inp[1] !== t || inp[2] != beams || inp[3] != talis_kwargs(tl, rv)
end
