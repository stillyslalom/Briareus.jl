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
The caveat on the rule tolerance that every TaLIS report carries: the requirements
add "up to 0.5 degrees" of metrology, alignment and positioning tolerance without
saying how an angle becomes a distance, so its interpretation (`pivot`) is an
assumption to confirm.
"""
const TALIS_TOLERANCE_CAVEAT =
    "ASSUMED: the requirements do not say how the angular tolerance becomes a distance; " *
    "it is applied here as a lateral allowance r·tan(tolerance). Confirm this " *
    "interpretation with the facility before relying on the margins."

"""
    TaLISState()

State of the TaLIS checks. Settings: `capture` (components beams are
meant to land on: backlighter foils, capture surfaces, a hohlraum wall), `hardware`
(other components checked for 3ω and 1ω clearance; capture components are checked
too, for beams landing elsewhere), `tolerance` (rule tolerance, degrees), `pivot` (one
of [`TALIS_PIVOTS`](@ref)) and `pivot_point` (mm, for `:point`), and `omega1` (run the
1ω rules). The error stack and the hardware sampling come from the
[`ReviewState`](@ref). Job: `running`, `started`, `status`, `error`, `report` (a core
`TaLISReport`) and `inputs` (the target and beams it checked).
"""
struct TaLISState
    capture::Observable{Vector{Symbol}}
    hardware::Observable{Vector{Symbol}}
    tolerance::Observable{Float64}
    pivot::Observable{Symbol}
    pivot_point::Observable{NTuple{3,Float64}}
    omega1::Observable{Bool}
    running::Observable{Bool}
    started::Observable{Float64}
    status::Observable{String}
    error::Observable{String}
    report::Observable{Union{Nothing,TaLISReport}}
    inputs::Base.RefValue{Any}
    known::Base.RefValue{Vector{Symbol}}
    generation::Base.RefValue{Int}
end

TaLISState() = TaLISState(Observable(Symbol[]), Observable(Symbol[]), Observable(0.5),
                          Observable(:tcc), Observable((0.0, 0.0, 0.0)), Observable(true),
                          Observable(false), Observable(0.0), Observable(""), Observable(""),
                          Observable{Union{Nothing,TaLISReport}}(nothing), Ref{Any}(nothing),
                          Ref(Symbol[]), Ref(0))

"""
    default_selection!(tl::TaLISState, t::Target)

Capture components and hardware for a target's components: a hohlraum's wall or a
backlighter foil captures the beams; the capsule is left out (light on it is drive
physics, handled by keep-outs and the drive, not a TaLIS hardware question; an
engineering default); everything else (LEH caps and lips, closed windows, the main
target) is checked. Checked LEH caps and lips go to `talis_check` as its `leh`, so
their 3ω clearance is an LEH-clearance review rather than a failure. Components the
selection was already made for keep their roles.
"""
function default_selection!(tl::TaLISState, t::Target)
    _update_roles!(tl.capture, tl.hardware, tl.known, surface_names(t),
                   n -> n in (:wall, FOIL_NAME) ? 1 : n === :capsule ? 0 : 2)
    return tl
end

# Components of the hohlraum template that form an entrance-hole edge.
const _LEH_EDGES = (:cap_top, :cap_bottom, :leh_top_lip, :leh_bottom_lip)

"""
    set_talis_role!(tl::TaLISState, component::Symbol, role::Symbol)

Make a component a `:capture` surface, checked `:hardware`, or neither (`:none`).
"""
function set_talis_role!(tl::TaLISState, name::Symbol, role::Symbol)
    role in (:capture, :hardware, :none) ||
        throw(ArgumentError("role must be :capture, :hardware or :none"))
    cp = filter(!=(name), tl.capture[])
    hw = filter(!=(name), tl.hardware[])
    role === :capture && push!(cp, name)
    role === :hardware && push!(hw, name)
    _assign!(tl.capture, cp)
    _assign!(tl.hardware, hw)
    return tl
end

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
checked LEH caps and lips (`leh`), the rule tolerance and its pivot, the clearance
check's error stack (pointing, positioning, hardware tolerance, focus) and hardware
sampling level.
"""
function talis_kwargs(tl::TaLISState, rv::ReviewState)
    e = rv.errors[]
    checked = Tuple(unique([tl.capture[]; tl.hardware[]]))
    pivot = tl.pivot[] === :point ? tl.pivot_point[] : tl.pivot[]
    return (; capture = Tuple(tl.capture[]), hardware = checked,
            leh = Tuple(n for n in tl.hardware[] if n in _LEH_EDGES), omega1_components = checked,
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
    tl.inputs[] = (t, bs)
    tl.started[] = time()
    tl.status[] = "checking…"
    tl.running[] = true
    job = function ()
        out = _try_job(() -> talis_check(t, bs; kw...))
        deliver(() -> _finish_talis!(tl, out, g))
    end
    _run_job(job, spawn)
    return tl
end

function _finish_talis!(tl::TaLISState, out, g::Integer)
    g == tl.generation[] || return tl                  # canceled or superseded
    if out.err === nothing
        tl.report[] = out.value
        tl.status[] = @sprintf("done in %s: ", _seconds(time() - tl.started[])) * talis_brief(tl)
    else
        tl.status[] = "failed: " * _errmsg(out.err)
    end
    tl.running[] = false
    return tl
end

"""
    cancel_talis!(tl::TaLISState)

Drop the check in progress (`talis_check` has no progress callback: the worker runs to
its end and its result is discarded). The previous report stays.
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
    talis_notes(tl::TaLISState) -> Vector{NamedTuple}

What the report does not settle, as `(; kind, text)` rows: the rule-tolerance
interpretation (`:assumed`, always first, also before a check), the error set
(`:settings`), the checks the report did not evaluate (`:skipped`), the requirements
out of scope (`:scope`, `TALIS_OUT_OF_SCOPE`), and the source (`:source`).
"""
function talis_notes(tl::TaLISState)
    r = tl.report[]
    notes = NamedTuple[]
    if r === nothing
        piv = tl.pivot[]
        txt = piv === :tcc ? "r = distance from TCC" : piv === :aim ? "r = distance from the beam's aim point" :
              piv === :max ? "r = larger of the distances from TCC and from the aim point" :
              @sprintf("r = distance from (%g, %g, %g) mm", tl.pivot_point[]...)
        push!(notes, (; kind = :assumed,
                      text = @sprintf("Rule tolerance %g° as a lateral allowance r·tan(%g°), %s. ",
                                      tl.tolerance[], tl.tolerance[], txt) * TALIS_TOLERANCE_CAVEAT))
    else
        s = r.settings
        push!(notes, (; kind = :assumed,
                      text = @sprintf("Rule tolerance %g° as a lateral allowance r·tan(%g°), %s. ",
                                      s.tolerance, s.tolerance, s.pivot_text) * TALIS_TOLERANCE_CAVEAT))
        push!(notes, (; kind = :settings,
                      text = @sprintf("Error set added linearly: transverse %g mm (pointing %g + positioning %g + hardware %g), focus ±%g mm.",
                                      s.transverse, s.pointing, s.target_offset,
                                      s.hardware_tolerance, s.focus)))
        for n in r.notes
            push!(notes, (; kind = :skipped, text = "Not evaluated: " * n))
        end
    end
    for n in TALIS_OUT_OF_SCOPE
        push!(notes, (; kind = :scope, text = "Out of scope: " * n))
    end
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
    talis_stale(tl::TaLISState, target, beams) -> Bool

Whether the target or the enabled beams changed since the last check.
"""
function talis_stale(tl::TaLISState, t::Target, beams::AbstractVector)
    inp = tl.inputs[]
    inp === nothing && return false
    return inp[1] !== t || inp[2] != beams
end
