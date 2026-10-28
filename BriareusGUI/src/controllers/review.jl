# The Review step: worst-case beam-to-hardware clearance (`clearance`) of the enabled
# beams on chosen hardware components, on a worker. `clearance` checks the beams in
# batches and reports progress through its callback between them; canceling returns
# the rows of the beams already checked. The step's TaLIS checks live in talis.jl.

"""
Clearance-check quality presets: `:quick` (hardware sampled at triangle vertices and
4 sub-triangle points, 8 pointing-error angles, a 31-point spot grid for power) and
`:standard` (`clearance`'s defaults: 16 sub-points, 16 angles, 61-point grid).
"""
const REVIEW_QUALITY = Dict(:quick => (; level = 1, nangles = 8, nspot = 31),
                            :standard => (; level = 2, nangles = 16, nspot = 61))

"""
    review_batch() -> Int

Beams per `clearance` batch (its `batch` keyword): the progress and cancel
granularity of the clearance check. One beam per thread, at least 4.
"""
review_batch() = max(4, Threads.nthreads())

"Views of the Review step: worst-case clearance, or the TaLIS rule checks."
const REVIEW_VIEWS = (:clearance, :talis)

"""
    ReviewState()

State of the clearance check. Settings: `targets` (components the beams should land on;
hardware beyond a beam's landing is not checked), `hardware` (components checked),
`peak_fraction` (envelope contour, fraction of each spot's peak), `wings` (include
the clearance-only `nif_wings` far field), the error stack `pointing`,
`target_offset`, `hardware_tolerance`, `focus` (mm; the TaLIS checks use it too), and
`quality` (a [`REVIEW_QUALITY`](@ref) key). `view` is the shown part of the step (one
of [`REVIEW_VIEWS`](@ref)). Job: `running`, `progress` `(beams done, total)`,
`started`, `status`, `report` (a core `ClearanceReport`, partial with `canceled` set
after a cancel), and `inputs` (the target, beams and `clearance` keywords it checked).
"""
struct ReviewState
    targets::Observable{Vector{Symbol}}
    hardware::Observable{Vector{Symbol}}
    peak_fraction::Observable{Float64}
    wings::Observable{Bool}
    errors::Observable{Dict{Symbol,Float64}}
    quality::Observable{Symbol}
    view::Observable{Symbol}
    running::Observable{Bool}
    progress::Observable{Tuple{Int,Int}}
    started::Observable{Float64}
    status::Observable{String}
    error::Observable{String}
    report::Observable{Union{Nothing,ClearanceReport}}
    inputs::Base.RefValue{Any}
    known::Base.RefValue{Vector{Symbol}}         # components the roles were chosen for
    cancel::Threads.Atomic{Bool}
end

ReviewState() = ReviewState(Observable(Symbol[]), Observable(Symbol[]), Observable(1e-4),
                            Observable(true),
                            Observable(Dict(:pointing => 0.05, :target_offset => 0.025,
                                            :hardware_tolerance => 0.0, :focus => 1.0)),
                            Observable(:quick), Observable(:clearance), Observable(false),
                            Observable((0, 0)), Observable(0.0), Observable(""), Observable(""),
                            Observable{Union{Nothing,ClearanceReport}}(nothing), Ref{Any}(nothing),
                            Ref(Symbol[]), Threads.Atomic{Bool}(false))

# Whether a component is where beams are meant to land by default (a hohlraum wall or
# capsule, a backlighter foil); everything else is hardware to check.
_lands_on(name::Symbol) = name in (:wall, :capsule, FOIL_NAME)

# Keep two role lists (`first`, `second`) in step with a target's components: roles of
# known components stay; new components get `role(name)` (`1`, `2`, or `0` for none);
# if no component keeps the first role, every role is chosen afresh.
function _update_roles!(first::Observable, second::Observable, known::Base.RefValue,
                        names::Vector{Symbol}, role)
    a = filter(in(names), first[])
    b = filter(in(names), second[])
    fresh = isempty(a)
    for n in names
        (fresh || !(n in known[])) || continue
        fresh && (filter!(!=(n), a); filter!(!=(n), b))
        r = role(n)
        r == 1 ? push!(a, n) : r == 2 && push!(b, n)
    end
    if isempty(a) && !isempty(names)              # no default landing component
        a = names[1:1]
        b = names[2:end]
    end
    known[] = copy(names)
    _assign!(first, a)
    _assign!(second, b)
    return nothing
end

"""
    default_selection!(rv::ReviewState, t::Target)

Choose the landing targets and hardware for a target's components: a hohlraum's beams
land on the wall and capsule and everything else (LEH caps, closed diagnostic
windows) is checked; a foil's beams land on the foil and the main target is checked.
Components the selection was already made for keep their roles; new components get
the default role.
"""
function default_selection!(rv::ReviewState, t::Target)
    _update_roles!(rv.targets, rv.hardware, rv.known, surface_names(t),
                   n -> _lands_on(n) ? 1 : 2)
    return rv
end

"""
    set_review_view!(rv::ReviewState, view::Symbol)

Show the clearance check or the TaLIS checks (one of [`REVIEW_VIEWS`](@ref)); the
viewer colors the beams by the shown report.
"""
function set_review_view!(rv::ReviewState, view::Symbol)
    view in REVIEW_VIEWS || throw(ArgumentError("view must be one of $REVIEW_VIEWS"))
    _assign!(rv.view, view)
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

The keyword arguments of the `clearance` call the clearance check makes.
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
    start_clearance!(rv::ReviewState, target, beams; deliver = f -> f(), spawn = true,
                     batch = review_batch())

Check `beams` against the selected hardware of `target` with `clearance` (on a worker
with `spawn`). `clearance` checks `batch` beams at a time and reports `(; done,
total)` through its callback after each batch, which updates `progress` through
`deliver`; [`cancel_clearance!`](@ref) makes the callback stop the check after the
batch in flight, and the report keeps the rows already computed (`canceled` set).
"""
function start_clearance!(rv::ReviewState, t::Target, beams::AbstractVector;
                          deliver = f -> f(), spawn::Bool = true, batch::Integer = review_batch())
    rv.running[] && return rv
    msg = isempty(beams) ? "no beams are enabled" :
          isempty(rv.targets[]) ? "choose the components the beams land on" :
          isempty(rv.hardware[]) ? "choose hardware to check" : nothing
    msg === nothing || (rv.status[] = msg; return rv)
    kw = clearance_kwargs(rv)
    bs = collect(beams)
    rv.cancel[] = false
    rv.inputs[] = (t, bs, kw)
    rv.progress[] = (0, length(bs))
    rv.started[] = time()
    rv.status[] = "checking…"
    rv.running[] = true
    job = function ()
        # runs on the job's task between batches: report progress, poll the cancel flag
        progress = function (p)
            done, total = p.done, p.total
            deliver(() -> (rv.progress[] = (done, total)))
            return rv.cancel[]
        end
        out = _try_job(() -> clearance(t, bs; kw..., callback = progress, batch))
        deliver(() -> _finish_clearance!(rv, out))
    end
    _run_job(job, spawn)
    return rv
end

function _finish_clearance!(rv::ReviewState, out)
    rep = out.value
    if out.err !== nothing
        rv.status[] = "failed: " * _errmsg(out.err)
    elseif rep.canceled
        # nothing checked yet: keep no report
        rv.report[] = isempty(rep.rows) ? nothing : rep
        done, total = rv.progress[]
        rv.status[] = "canceled after $done of $total beams"
    else
        rv.report[] = rep
        rv.status[] = "done: " * review_brief(rv)
    end
    rv.running[] = false
    return rv
end

"""
    cancel_clearance!(rv::ReviewState)

Stop after the batch of beams in flight (`clearance`'s callback returns `true`); the
rows already computed are kept.
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
    return @sprintf("%d beams%s: %d fail, %d review · worst clearance %.3f mm", nb,
                    r.canceled ? " (canceled, partial)" : "", nf, nr, worst)
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

Whether anything the last check depends on changed since: the target, the enabled
beams, or any setting (`clearance_kwargs`: landing targets and hardware, contour,
wings, error stack, sampling).
"""
function review_stale(rv::ReviewState, t::Target, beams::AbstractVector)
    inp = rv.inputs[]
    inp === nothing && return false
    return inp[1] !== t || inp[2] != beams || inp[3] != clearance_kwargs(rv)
end
