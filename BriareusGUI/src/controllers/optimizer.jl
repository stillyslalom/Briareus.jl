# The Optimize step: a GroupPointing parameterization of the enabled beams (fields,
# grouping, box bounds), an objective (capsule symmetry, axisymmetric or 3D, or foil
# uniformity) and `optimize_pointing` on a worker. Its callback reports each objective
# evaluation (progress) and stops the optimization on cancel, which then returns the
# best design found so far.

"""
Objectives of the pointing optimization: `:legendre` (`symmetry_objective`: weighted
`(Pₗ/P₀)²` of the capsule flux, for axisymmetric designs), `:ylm` (`asymmetry_objective`:
weighted rotation-invariant mode power over all `m`, for 3D designs such as a dropped
quad), `:irradiance` (`irradiance_objective`: foil-spot uniformity and spill).
"""
const OBJECTIVES = (:legendre, :ylm, :irradiance)

const OBJECTIVE_LABELS = Dict(:legendre => "Capsule symmetry, axisymmetric (Pₗ/P₀)",
                              :ylm => "Capsule symmetry, 3D (all m)",
                              :irradiance => "Foil spot uniformity")

# Default bound half-widths: mm for z, r, phi, defocus; a fraction of the start for power.
const DEFAULT_DELTAS = Dict(:z => 1.0, :r => 1.0, :phi => 1.0, :defocus => 2.0, :power => 0.4)

"""
    OptimizeResult

A finished (or canceled) optimization: `names` (parameter names), `x0`, `x`, `lower`,
`upper`, `initial_value`, `value`, `beams` (the optimized enabled beams), `indices`
(their positions in the beam table), `before`/`after` (quick summaries: `modes`,
`mean_flux`, `spot`, `clipped`, `capsule` (power fraction on the capsule),
`capsule_clearance` (worst-case, mm)), `evaluations` (objective evaluations),
`iterations` (L-BFGS-B iterations completed), `reason` (why it ended: `:converged`,
`:iterations`, `:evaluations` or `:callback`), `seconds`, `completed` (`false` when
canceled: the best design found so far), and `applied` (set by
[`accept_optimum!`](@ref)).
"""
mutable struct OptimizeResult
    names::Vector{String}
    x0::Vector{Float64}
    x::Vector{Float64}
    lower::Vector{Float64}
    upper::Vector{Float64}
    initial_value::Float64
    value::Float64
    beams::Vector{Beam{Float64}}
    indices::Vector{Int}
    before::NamedTuple
    after::NamedTuple
    evaluations::Int
    iterations::Int
    reason::Symbol
    seconds::Float64
    completed::Bool
    applied::Bool
end

"""
    OptimizeState()

State of the pointing optimization. Settings: `by` (grouping, one of
[`BEAM_GROUPINGS`](@ref)), `fields` (free `GroupPointing` fields, a subset of
[`BEAM_FIELDS`](@ref)), `deltas` (bound half-widths per field: mm, or for `:power` a
fraction of the start), `objective` (one of [`OBJECTIVES`](@ref)), `weights` (mode
weight per degree `l`), `hold_drive` (keep the capsule drive at its start value),
`keepout`/`keepout_min` (an `ApertureKeepOut` on the LEH rims with this minimum
worst-case clearance, mm), `capsule_keepout`/`capsule_min` (a `SphereKeepOut` on the
capsule, on by default: without it the optimizer may move the inner cones onto the
capsule, which the deposition penalty barely sees), `spill` (foil objective weight on
power outside the spot), `iterations` (L-BFGS-B iterations) and `max_evaluations`
(objective evaluations; a line search can take many per iteration). Job: `running`,
`evaluations`, `iteration`, `current`/`best` objective values, `started`, `status`,
and `result` (an [`OptimizeResult`](@ref)).
"""
struct OptimizeState
    by::Observable{Symbol}
    fields::Observable{Vector{Symbol}}
    deltas::Observable{Dict{Symbol,Float64}}
    objective::Observable{Symbol}
    weights::Observable{Dict{Int,Float64}}
    hold_drive::Observable{Bool}
    keepout::Observable{Bool}
    keepout_min::Observable{Float64}
    capsule_keepout::Observable{Bool}
    capsule_min::Observable{Float64}
    spill::Observable{Float64}
    iterations::Observable{Int}
    max_evaluations::Observable{Int}
    running::Observable{Bool}
    evaluations::Observable{Int}
    iteration::Observable{Int}
    current::Observable{Float64}
    best::Observable{Float64}
    started::Observable{Float64}
    status::Observable{String}
    error::Observable{String}
    result::Observable{Union{Nothing,OptimizeResult}}
    cancel::Threads.Atomic{Bool}
end

OptimizeState() =
    OptimizeState(Observable(:group), Observable([:z, :r]), Observable(copy(DEFAULT_DELTAS)),
                  Observable(:legendre), Observable(Dict(1 => 1.0, 2 => 1.0, 3 => 1.0, 4 => 1.0)),
                  Observable(true), Observable(false), Observable(0.0), Observable(true),
                  Observable(0.0), Observable(0.2), Observable(100), Observable(1000),
                  Observable(false), Observable(0), Observable(0), Observable(NaN),
                  Observable(NaN), Observable(0.0), Observable(""), Observable(""),
                  Observable{Union{Nothing,OptimizeResult}}(nothing), Threads.Atomic{Bool}(false))

"""
    set_free!(st::OptimizeState, field::Symbol, on::Bool)

Free (`true`) or fix a `GroupPointing` field (one of [`BEAM_FIELDS`](@ref)).
"""
function set_free!(st::OptimizeState, field::Symbol, on::Bool)
    field in BEAM_FIELDS || throw(ArgumentError("field must be one of $BEAM_FIELDS"))
    fs = filter(!=(field), st.fields[])
    on && push!(fs, field)
    _assign!(st.fields, [f for f in BEAM_FIELDS if f in fs])     # canonical order
    return st
end

"""
    edit_optimize!(st::OptimizeState, key::Symbol, text) -> Bool

Edit a setting from text: a bound half-width (`:z`, `:r`, `:phi`, `:defocus`,
`:power`), a mode weight (`:w1`…`:w4`), `:keepout_min`, `:capsule_min`, `:spill`,
`:iterations` or `:max_evaluations`.
A rejected value sets `error` and returns `false`.
"""
function edit_optimize!(st::OptimizeState, key::Symbol, text)
    try
        if key in BEAM_FIELDS
            v = _parse_number(text, "the $key bound"; positive = true)
            key === :power && v >= 1 && throw(ArgumentError("the power bound is a fraction below 1"))
            d = copy(st.deltas[]); d[key] = v; _assign!(st.deltas, d)
        elseif key in (:w1, :w2, :w3, :w4)
            l = parse(Int, String(key)[2:end])
            w = copy(st.weights[]); w[l] = _parse_number(text, "the weight"; nonnegative = true)
            _assign!(st.weights, w)
        elseif key === :keepout_min
            _assign!(st.keepout_min, _parse_number(text, "the keep-out clearance"))
        elseif key === :capsule_min
            _assign!(st.capsule_min, _parse_number(text, "the capsule clearance"))
        elseif key === :spill
            _assign!(st.spill, _parse_number(text, "the spill weight"; nonnegative = true))
        elseif key === :iterations
            n = round(Int, _parse_number(text, "the iteration limit"; positive = true))
            _assign!(st.iterations, n)
        elseif key === :max_evaluations
            n = round(Int, _parse_number(text, "the evaluation limit"; positive = true))
            _assign!(st.max_evaluations, n)
        else
            throw(ArgumentError("unknown optimization setting :$key"))
        end
        st.error[] = ""
        return true
    catch err
        err isa ArgumentError || rethrow()
        st.error[] = _errmsg(err)
        return false
    end
end

"""
    set_objective!(st::OptimizeState, objective::Symbol)

One of [`OBJECTIVES`](@ref).
"""
function set_objective!(st::OptimizeState, objective::Symbol)
    objective in OBJECTIVES || throw(ArgumentError("objective must be one of $OBJECTIVES"))
    _assign!(st.objective, objective)
    return st
end

# The objective settings, captured on the GUI thread when a job starts.
_objective_settings(st::OptimizeState) =
    (; kind = st.objective[], weights = copy(st.weights[]), hold_drive = st.hold_drive[],
     keepout = st.keepout[], keepout_min = st.keepout_min[],
     capsule_keepout = st.capsule_keepout[], capsule_min = st.capsule_min[], spill = st.spill[])

# The objective for a scenario of the enabled beams, given its start evaluation `ev0`.
function _make_objective(cfg::NamedTuple, s::Scenario, ev0::Evaluation, spot_radius::Real)
    kind = cfg.kind
    if kind === :irradiance
        FOIL_NAME in component_names(s.target) ||
            throw(ArgumentError("the foil objective needs a foil target"))
        c = sum(s.target.elem_centroid[e] for e in elements(s.target, FOIL_NAME)) /
            length(elements(s.target, FOIL_NAME))
        return irradiance_objective(FOIL_NAME; center = c, radius = spot_radius, spill = cfg.spill)
    end
    ev0.flux === nothing && throw(ArgumentError("capsule symmetry needs a hohlraum with a capsule"))
    modes = Tuple(l => w for (l, w) in sort(collect(cfg.weights)) if w > 0)
    isempty(modes) && throw(ArgumentError("give at least one mode a positive weight"))
    keepouts = Any[]
    cfg.keepout && push!(keepouts, ApertureKeepOut(s.target; minimum = cfg.keepout_min))
    cfg.capsule_keepout && :capsule in component_names(s.target) &&
        push!(keepouts, SphereKeepOut(s.target; spheres = (:capsule,), minimum = cfg.capsule_min))
    keepout = isempty(keepouts) ? nothing : Tuple(keepouts)
    if kind === :legendre
        drive = cfg.hold_drive ? flux_moments(ev0)[1] : nothing
        return symmetry_objective(; modes, odd = 0.0, lmax = max(maximum(first, modes), 2),
                                  drive, keepout)
    end
    drive = cfg.hold_drive ? flux_ylm(ev0)[1] / sqrt(4π) : nothing
    return asymmetry_objective(; modes, drive, keepout)
end

# Capsule modes / foil spot and clipping of an evaluation, for before/after tables.
function _quick_summary(ev::Evaluation, spot_radius::Real)
    modes = NamedTuple[]
    mean = NaN
    if ev.flux !== nothing
        c = flux_ylm(ev; lmax = maximum(MODE_DEGREES))
        a = flux_moments(ev; lmax = maximum(MODE_DEGREES))
        mean = a[1]
        modes = [(; l, epsilon = mode_amplitude(c, l), legendre = a[l+1] / a[1]) for l in MODE_DEGREES]
    end
    t = ev.scenario.target
    has = :capsule in component_names(t)
    cap = has ? minimum(sphere_clearance(SphereKeepOut(t; spheres = (:capsule,)), ev.beams);
                        init = Inf) : Inf
    return (; modes, mean_flux = mean, spot = _spot_stats(ev, spot_radius),
            clipped = clipped_fraction(ev), capsule = has ? power_fraction(ev, :capsule) : 0.0,
            capsule_clearance = cap)
end

"""
    start_optimize!(st::OptimizeState, scenario, beams, indices; spot_radius = 0.5,
                    deliver = f -> f(), spawn = true)

Optimize the enabled `beams` (positions `indices` in the beam table; `scenario` is
built for the whole table) with the current settings. Runs `optimize_pointing` on a
worker with `spawn`; per-group parameters follow `clusters` (the beams' group names in
the beam table) when given, otherwise the beams' cones. Its callback hands progress (`evaluations`, `current`, `best`,
at most every 0.1 s) to `deliver`, and so does the result. [`cancel_optimize!`](@ref)
makes the callback stop the optimization after the next objective evaluation; the
result is then the best design found so far.
"""
function start_optimize!(st::OptimizeState, s::Scenario, beams::AbstractVector,
                         indices::AbstractVector{<:Integer}; spot_radius::Real = 0.5,
                         clusters = nothing, deliver = f -> f(), spawn::Bool = true)
    st.running[] && return st
    try
        isempty(beams) && throw(ArgumentError("no beams are enabled"))
        isempty(st.fields[]) && throw(ArgumentError("free at least one field"))
        _check_grouping(beams, st.by[])
    catch err
        st.status[] = _errmsg(err)
        return st
    end
    by, fields = st.by[], Tuple(st.fields[])
    iterations, evaluations = st.iterations[], st.max_evaluations[]
    st.cancel[] = false
    st.evaluations[] = 0
    st.iteration[] = 0
    st.current[] = NaN
    st.best[] = NaN
    st.started[] = time()
    st.status[] = "starting…"
    st.running[] = true
    settings = (; deltas = copy(st.deltas[]))
    cfg = _objective_settings(st)
    job = function ()
        t0 = time()
        key = by === :group && clusters !== nothing ?
              (d = Dict(b.name => c for (b, c) in zip(beams, clusters)); b -> d[b.name]) : by
        p = GroupPointing(beams; by = key, fields, mirror = _mirror(by))
        sub = subset(s, collect(indices); parameterization = p)
        x0 = initial_parameters(p)
        d = settings.deltas
        lower, upper = bounds(p, x0; z = d[:z], r = d[:r], phi = d[:phi],
                              defocus = d[:defocus], power = d[:power])
        ev0 = evaluate(sub, x0)
        before = _quick_summary(ev0, spot_radius)
        obj = _make_objective(cfg, sub, ev0, spot_radius)
        # called on this task after every objective evaluation; `true` stops
        reported = Ref(0.0)
        progress = function (info)
            if time() - reported[] > 0.1
                reported[] = time()
                n, it, J, b = info.evaluations, info.iteration, info.value, info.best_value
                deliver(() -> (st.evaluations[] = n; st.iteration[] = it; st.current[] = J;
                               st.best[] = b))
            end
            return st.cancel[]
        end
        res = optimize_pointing(sub, obj; x0, lower, upper, iterations, evaluations,
                                callback = progress)
        after = _quick_summary(evaluate(sub, res.x), spot_radius)
        return OptimizeResult(parameter_names(p), x0, collect(Float64, res.x), lower, upper,
                              res.initial_value, res.value,
                              Beam{Float64}[_float_beam(b) for b in res.beams],
                              collect(Int, indices), before, after, res.evaluations,
                              res.iterations, res.reason, time() - t0, res.reason !== :callback,
                              false)
    end
    run = function ()
        out = _try_job(job)
        deliver(() -> _finish_optimize!(st, out))
    end
    _run_job(run, spawn)
    return st
end

function _finish_optimize!(st::OptimizeState, out)
    if out.err === nothing
        v = out.value
        st.result[] = v
        st.evaluations[] = v.evaluations
        st.best[] = v.value
        how = v.reason === :converged ? "converged" :
              v.reason === :iterations ? "reached the iteration limit" :
              v.reason === :evaluations ? "reached the evaluation limit; kept the best design" :
              "canceled; kept the best design so far"
        st.status[] = how * @sprintf(": objective %.4g → %.4g in %d iterations, %d evaluations, %s",
                                     v.initial_value, v.value, v.iterations, v.evaluations,
                                     _seconds(v.seconds))
    else
        st.status[] = "failed: " * _errmsg(out.err)
    end
    st.running[] = false
    return st
end

"""
    cancel_optimize!(st::OptimizeState)

Stop the optimization after the objective evaluation in progress (the
`optimize_pointing` callback returns `true`); the best design found so far becomes
the result.
"""
cancel_optimize!(st::OptimizeState) = (st.cancel[] = true; st)

"""
    optimize_progress(st::OptimizeState) -> String

Evaluations, the latest and best objective values, and elapsed time while running.
"""
function optimize_progress(st::OptimizeState)
    st.running[] || return st.status[]
    n = st.evaluations[]
    el = _seconds(time() - st.started[])
    n == 0 && return "setting up… · $el"
    return @sprintf("iteration %d of %d · %d of %d evaluations · objective %.4g (best %.4g) · %s",
                    st.iteration[], st.iterations[], n, st.max_evaluations[], st.current[],
                    st.best[], el)
end

"""
    parameter_rows(res::OptimizeResult) -> Vector{NamedTuple}

Rows `(; name, start, optimized, lower, upper, at_bound)` of an optimization's
parameters (mm, or TW per beam for power); `at_bound` flags a parameter within 0.1 %
of the bound range from a bound.
"""
function parameter_rows(res::OptimizeResult)
    return map(eachindex(res.names)) do i
        lo, hi, x = res.lower[i], res.upper[i], res.x[i]
        tol = 1e-3 * max(hi - lo, eps())
        (; name = res.names[i], start = res.x0[i], optimized = x, lower = lo, upper = hi,
         at_bound = x - lo < tol || hi - x < tol)
    end
end
