# The Review step: worst-case beam-to-hardware clearance (`clearance`) of the enabled
# beams on chosen hardware components. Each beam's rows are independent, so the check
# runs in chunks of beams on a worker: progress between chunks, cancel keeps the rows
# already computed.

"""
Clearance-check quality presets: `:quick` (hardware sampled at triangle vertices and
4 sub-triangle points, 8 pointing-error angles, a 31-point spot grid for power) and
`:standard` (`clearance`'s defaults: 16 sub-points, 16 angles, 61-point grid).
"""
const REVIEW_QUALITY = Dict(:quick => (; level = 1, nangles = 8, nspot = 31),
                            :standard => (; level = 2, nangles = 16, nspot = 61))

# Beams per clearance call: progress and cancel granularity.
const REVIEW_CHUNK = 4

"""
    ReviewState()

State of the Review step. Settings: `targets` (components the beams should land on;
hardware beyond a beam's landing is not checked), `hardware` (components checked),
`peak_fraction` (envelope contour, fraction of each spot's peak), `wings` (include
the clearance-only `nif_wings` far field), the error stack `pointing`,
`target_offset`, `hardware_tolerance`, `focus` (mm), and `quality` (a
[`REVIEW_QUALITY`](@ref) key). Job: `running`, `progress` `(beams done, total)`,
`started`, `status`, `report` (a core `ClearanceReport`, partial after a cancel),
and `inputs` (the target and beams it checked).
"""
struct ReviewState
    targets::Observable{Vector{Symbol}}
    hardware::Observable{Vector{Symbol}}
    peak_fraction::Observable{Float64}
    wings::Observable{Bool}
    errors::Observable{Dict{Symbol,Float64}}
    quality::Observable{Symbol}
    running::Observable{Bool}
    progress::Observable{Tuple{Int,Int}}
    started::Observable{Float64}
    status::Observable{String}
    error::Observable{String}
    report::Observable{Union{Nothing,ClearanceReport}}
    inputs::Base.RefValue{Any}
    cancel::Threads.Atomic{Bool}
end

ReviewState() = ReviewState(Observable(Symbol[]), Observable(Symbol[]), Observable(1e-4),
                            Observable(true),
                            Observable(Dict(:pointing => 0.05, :target_offset => 0.025,
                                            :hardware_tolerance => 0.0, :focus => 1.0)),
                            Observable(:quick), Observable(false), Observable((0, 0)),
                            Observable(0.0), Observable(""), Observable(""),
                            Observable{Union{Nothing,ClearanceReport}}(nothing), Ref{Any}(nothing),
                            Threads.Atomic{Bool}(false))

"""
    default_selection!(rv::ReviewState, t::Target)

Choose the landing targets and hardware for a target: a hohlraum's beams land on the
wall (and capsule) and the LEH caps are checked; a foil's beams land on the foil and
the main target is checked. Keeps a selection that still names existing components.
"""
function default_selection!(rv::ReviewState, t::Target)
    names = surface_names(t)
    keep(v) = !isempty(v) && all(in(names), v)
    if !(keep(rv.targets[]) && all(in(names), rv.hardware[]))
        if :wall in names
            _assign!(rv.targets, [n for n in (:wall, :capsule) if n in names])
            _assign!(rv.hardware, [n for n in (:cap_top, :cap_bottom) if n in names])
        elseif FOIL_NAME in names
            _assign!(rv.targets, [FOIL_NAME])
            _assign!(rv.hardware, filter(!=(FOIL_NAME), names))
        else
            _assign!(rv.targets, names[1:min(1, end)])
            _assign!(rv.hardware, names[2:end])
        end
    end
    return rv
end

"""
    set_role!(rv::ReviewState, component::Symbol, role::Symbol)

Make a component a landing `:target`, checked `:hardware`, or neither (`:none`).
"""
function set_role!(rv::ReviewState, name::Symbol, role::Symbol)
    role in (:target, :hardware, :none) || throw(ArgumentError("role must be :target, :hardware or :none"))
    tg = filter(!=(name), rv.targets[])
    hw = filter(!=(name), rv.hardware[])
    role === :target && push!(tg, name)
    role === :hardware && push!(hw, name)
    _assign!(rv.targets, tg)
    _assign!(rv.hardware, hw)
    return rv
end

"""
    edit_review!(rv::ReviewState, key::Symbol, text) -> Bool

Edit `:peak_fraction` or an error-stack entry (`:pointing`, `:target_offset`,
`:hardware_tolerance`, `:focus`, mm) from text. A rejected value sets `error` and
returns `false`.
"""
function edit_review!(rv::ReviewState, key::Symbol, text)
    try
        if key === :peak_fraction
            v = _parse_number(text, "the contour level"; positive = true)
            v < 1 || throw(ArgumentError("the contour level is a fraction of the peak, below 1"))
            _assign!(rv.peak_fraction, v)
        elseif haskey(rv.errors[], key)
            v = _parse_number(text, "the $(replace(String(key), '_' => ' ')) error"; nonnegative = true)
            e = copy(rv.errors[]); e[key] = v; _assign!(rv.errors, e)
        else
            throw(ArgumentError("unknown review setting :$key"))
        end
        rv.error[] = ""
        return true
    catch err
        err isa ArgumentError || rethrow()
        rv.error[] = _errmsg(err)
        return false
    end
end

"""
    clearance_kwargs(rv::ReviewState) -> NamedTuple

The keyword arguments of the `clearance` call the Review step makes.
"""
function clearance_kwargs(rv::ReviewState)
    e = rv.errors[]
    return (; targets = Tuple(rv.targets[]), hardware = Tuple(rv.hardware[]),
            peak_fraction = rv.peak_fraction[], wings = rv.wings[] ? nif_wings() : SpotWing[],
            pointing = e[:pointing], target_offset = e[:target_offset],
            hardware_tolerance = e[:hardware_tolerance], focus = e[:focus],
            REVIEW_QUALITY[rv.quality[]]...)
end

"""
    start_clearance!(rv::ReviewState, target, beams; deliver = f -> f(), spawn = true)

Check `beams` against the selected hardware of `target` with `clearance`, in chunks
of [`REVIEW_CHUNK`](@ref) beams (on a worker with `spawn`). [`cancel_clearance!`](@ref)
stops after the chunk in flight and keeps the rows already computed.
"""
function start_clearance!(rv::ReviewState, t::Target, beams::AbstractVector;
                          deliver = f -> f(), spawn::Bool = true)
    rv.running[] && return rv
    msg = isempty(beams) ? "no beams are enabled" :
          isempty(rv.targets[]) ? "choose the components the beams land on" :
          isempty(rv.hardware[]) ? "choose hardware to check" : nothing
    msg === nothing || (rv.status[] = msg; return rv)
    kw = clearance_kwargs(rv)
    bs = collect(beams)
    rv.cancel[] = false
    rv.inputs[] = (t, bs)
    rv.progress[] = (0, length(bs))
    rv.started[] = time()
    rv.status[] = "checking…"
    rv.running[] = true
    job = function ()
        rows = ClearanceRow[]
        head = nothing
        outcome = try
            for lo in 1:REVIEW_CHUNK:length(bs)
                rv.cancel[] && throw(JobCancelled())
                hi = min(lo + REVIEW_CHUNK - 1, length(bs))
                rep = clearance(t, bs[lo:hi]; kw...)
                head === nothing && (head = rep)
                append!(rows, rep.rows)
                deliver(() -> (rv.progress[] = (hi, length(bs))))
            end
            :done
        catch err
            err isa JobCancelled ? :canceled : err
        end
        report = head === nothing ? nothing :
                 ClearanceReport(rows, head.contour, head.transverse, head.errors, head.margin,
                                 head.targets, head.thresholds)
        deliver(() -> _finish_clearance!(rv, report, outcome))
    end
    _run_job(job, spawn)
    return rv
end

function _finish_clearance!(rv::ReviewState, report, outcome)
    rv.report[] = report
    done, total = rv.progress[]
    rv.status[] = outcome === :done ? "done: " * review_brief(rv) :
                  outcome === :canceled ? "canceled after $done of $total beams" :
                  "failed: " * _errmsg(outcome)
    rv.running[] = false
    return rv
end

"""
    cancel_clearance!(rv::ReviewState)

Stop after the chunk of beams in flight; rows already computed are kept.
"""
cancel_clearance!(rv::ReviewState) = (rv.cancel[] = true; rv)

"""
    clearance_rows(rv::ReviewState) -> Vector{ClearanceRow}

Every row of the report, failures first, then review, then pass, each by worst-case
clearance.
"""
clearance_rows(rv::ReviewState) =
    rv.report[] === nothing ? ClearanceRow[] : violations(rv.report[]; status = (:fail, :review, :pass))

"""
    worst_rows(rv::ReviewState) -> Vector{ClearanceRow}

The closest beam for each checked component (`min_clearance`), in component order.
"""
function worst_rows(rv::ReviewState)
    rv.report[] === nothing && return ClearanceRow[]
    m = min_clearance(rv.report[])
    return [m[k] for k in rv.hardware[] if haskey(m, k)]
end

"""
    review_brief(rv::ReviewState) -> String

Fail/review counts and the worst-case clearance, in a line.
"""
function review_brief(rv::ReviewState)
    r = rv.report[]
    r === nothing && return "not run"
    isempty(r.rows) && return "no rows"
    nf = count(x -> x.status === :fail, r.rows)
    nr = count(x -> x.status === :review, r.rows)
    worst = minimum(x -> x.clearance, r.rows)
    nb = length(unique(x.beam for x in r.rows))
    return @sprintf("%d beams: %d fail, %d review · worst clearance %.3f mm", nb, nf, nr, worst)
end

"""
    export_clearance(rv::ReviewState, path) -> path

Write the report with `write_clearance_csv`.
"""
function export_clearance(rv::ReviewState, path::AbstractString)
    rv.report[] === nothing && throw(ArgumentError("run the clearance check first"))
    write_clearance_csv(path, rv.report[])
    rv.status[] = "wrote $(basename(path))"
    return path
end

"""
    review_stale(rv::ReviewState, target, beams) -> Bool

Whether the target or the enabled beams changed since the last check.
"""
function review_stale(rv::ReviewState, t::Target, beams::AbstractVector)
    inp = rv.inputs[]
    inp === nothing && return false
    return inp[1] !== t || inp[2] != beams
end
