# This file runs the pointing optimization: `optimize_pointing`, the L-BFGS-B driver with
# its threaded ForwardDiff gradients, and the `PointingResult` with its report.

"""
    optimize_pointing(scenario, objective; from = nothing,
                      x0 = initial_parameters(scenario.param),
                      lower = bounds(scenario.param, x0)[1], upper = bounds(scenario.param, x0)[2],
                      iterations = 200, evaluations = 20 * iterations, show_trace = false,
                      g_tol = 1e-8, callback = nothing, constraints = (),
                      feasibility = 1e-3, outer = 30, time_limit = Inf, phase1 = true,
                      relaxation = false, inner_reltol = 1e-10,
                      inner_evaluations = nothing, smooth_cut = true) -> PointingResult

Minimize `objective(evaluate(scenario, x))` within the box `lower ≤ x ≤ upper` by
L-BFGS-B with forward-mode automatic-differentiation gradients, from `x0` (clamped to
the box). Returns a [`PointingResult`](@ref), whose display lists the parameters that
end at a bound.

The default box is `bounds(scenario.param, x0)` with its default widths. These are aims
within ±1 mm of `x0` (`z`, `r`, `phi`), defocus ±2 mm, power and spot scale ±50%,
Δλ ±2 Å and the spherical angles ±5–10°. They are engineering defaults and are not
facility limits. A parameter at its bound usually means that the box set its value and
the objective without the box would have moved it further. To widen the box, pass keywords
of [`bounds`](@ref), e.g.
`lower, upper = bounds(pointing, x0; z = 1.5, power = 0.4, power_max = 12.0)`, or pass
any `lower` and `upper` vectors. A parameterization other than `GroupPointing` and
[`FramePointing`](@ref) has no `bounds` method and needs `lower` and `upper`.

`scenario` may be a [`PhasedScenario`](@ref), which contains one pointing over several
phases of a pulse with power scales per (group, phase). The objective and every
constraint then receive the [`PerPhase`](@ref) of the evaluations of the phases at each
parameter vector, evaluated once. The `beams` of the result are a `PerPhase` of beam
vectors. Everything else is the same, including constraints, `from`, `phase1`,
`relaxation` and the time limit.

The run ends after `iterations` L-BFGS-B iterations or `evaluations` objective
evaluations, whichever comes first; each gradient counts as one evaluation. A line
search may take many evaluations in one iteration, and dozens on a kinked objective
such as keep-out hinges or a coarse mesh, so the iteration limit alone does not bound
the run time. `time_limit`, in seconds of wall time (default none), ends the run at the
first evaluation after it has passed. `reason` is one of:

- `:converged`: the projected gradient fell below `g_tol`, or the run made no further
  progress;
- `:iterations`, `:evaluations`, `:time_limit` or `:callback`: that limit or the
  callback ended the run.

On `:evaluations`, `:time_limit` and `:callback`, the result is the best design
evaluated before the stop, with `stopped = true` and `result = nothing`. `iterations` is
the number of iterations completed in every case.

`callback`, if given, is called after every (non-derivative) objective evaluation with
`(; evaluations, iteration, value, best_value, best_x)`; returning `true` stops the
optimization. Use it for progress reports and cancellation.

**Constraints.** `constraints` is a tuple or vector of [`Constraint`](@ref)s: ceilings
`f(ev) ≤ limit` and floors `f(ev) ≥ limit`, such as [`irradiance_ceiling`](@ref) and
[`mode_gate`](@ref). They are handled as constraints by an augmented Lagrangian; they are
not penalty terms of the objective. Its inner problems are box-constrained L-BFGS-B
runs on `objective + Σ μ/2 [max(0, gᵢ + λᵢ/μ)² - (λᵢ/μ)²]`, where `gᵢ` are the
normalized violations `±(value - limit)/scale`. Each inner run starts from the last, and
the multipliers `λ` and the penalty `μ` are updated between them by the rules of
ALGENCAN (Andreani et al., SIAM J. Optim. 18, 1286, 2008).

The run converges when three conditions hold:

- every violation is at most `feasibility`, in the scaled units of the constraints
  (1e-3, or 0.1% of a limit, by default);
- the multipliers are complementary to the same tolerance;
- the last inner run has converged at the final tolerance.

The tolerance of the inner runs on the projected gradient tightens tenfold per outer
iteration, from 1e-3 to `g_tol`. The run ends only once an inner run at `g_tol` has
converged. An inner run converges when its projected gradient falls below the tolerance,
or when its value changes by less than `inner_reltol` relative per iteration, since a
drive objective is flat along the limits. A start already within the limits, as on a
restart or a hand-off to a finer parameterization, is therefore still optimized. A
constraint that cannot be met ends the run at the verdict of phase 1 (below). With
`phase1 = false`, it ends the run once the penalty has grown by 10⁸ without progress.

`outer` caps the outer iterations, `iterations` caps each inner run, and `evaluations`
is the budget of the whole run. The `iterations` and `evaluations` of the result count
over all inner runs. The multipliers change only between inner runs, so the cap sets
the number of iterations spent on each set of multipliers. An inner run ended by its
cap has not solved its subproblem, and the penalty test (has the
infeasibility–complementarity measure halved?) requires a solved one. After a capped
run that ends above the limits, the multipliers are updated, while the penalty, the
inner tolerance and the reference of the test are held. A cap therefore never grows the
penalty toward the infeasible verdict. The usual test is applied to a capped run that
ends within the limits. `capped` in the result counts both kinds (research notes, "Least
violation, capped inner runs and phase 1's budget"). A run that ends above its limits,
infeasible or stopped by a budget, the time limit or the callback, returns its iterate
of least violation instead of its last (the `least_violation` and `last_feasible` fields
of [`PointingResult`](@ref)).

The public cases of up to about a hundred parameters converge with a large cap in about
a thousand evaluations; the dropped-bundle tutorial passes 1,000. Small early caps were
measured as an alternative: they saved evaluations on some of these cases and cost drive
on others. A problem with hundreds of parameters, whose evaluations take seconds and
whose start is far from the limits, gains from a cap of about 100, which updates the
multipliers sooner. A 276-parameter study completed one outer iteration in an hour with
`iterations = 1000`, and six in 40 minutes with 100 (research notes, "The inner
budget").

The cap on iterations does not bound the evaluations of an inner run, which line
searches consume. An inner run can spend hundreds of them on 1e-4 of the objective
while limits stay violated. `inner_evaluations` (default none) caps the evaluations of
each inner run. An inner run that reaches it ends at its best design, as one ended by
`iterations` does, and the outer loop updates the multipliers and the penalty from
there. Phase 1 is not capped; the runs of the relaxations are. A capped inner run never
counts as converged, so the run still ends only after an inner run converges within
the cap at the final tolerance, and a small cap needs a larger `outer`. A coarse version
of the run within limits in the dropped-bundle tutorial has 92 parameters, and its
default run spent 795 of its 1,044 evaluations in one inner run. With a cap of 100, it
converged in 10 outer iterations and 766 evaluations at 0.06% less drive. With a cap of
30, it ended at `outer = 30` unconverged (research notes, "Inner-run options, a
wall-loss scenario warning and laser delivery figures").

`inner_reltol` (default 1e-10; 0 turns it off) is the relative value tolerance of the
inner runs above, also used by phase 1. The default ends the inner runs of a drive
objective where the objective is flat along the limits, short of stationarity; the
distance is recorded in the `stationarity` of the result. A tighter tolerance costs
evaluations. On the run within limits in the dropped-bundle tutorial, a tolerance of 1e-12
gives a stationarity of 1.9e-5 against 1.2e-3 for 2.1 times the evaluations, with the drive
the same within the 1e-4 spread of the flat family (research notes, "Multipliers from the
converged active set"). When the multipliers or bound sensitivities of a result matter,
continue it tighter through `from` with a smaller `inner_reltol`.

The `constraints` and `outcome` of the result report each constraint against its limit,
with the multipliers that best satisfy stationarity at the final design. They also
report whether the box, the constraints or the budget ended the run
([`PointingResult`](@ref)). The named tuple of the callback then also contains `outer`,
and its `value` is that of the inner problem, the objective with the constraint terms.
Without constraints, the run is the plain L-BFGS-B above.

**Infeasible constraints.** With `phase1 = true` (the default), a start that violates a
constraint first goes through phase 1. Phase 1 minimizes the infeasibility
`½ Σ max(0, gᵢ)²` over the normalized violation of every component, without the
objective, by the same L-BFGS-B. This infeasibility is the limit of the augmented
Lagrangian as its penalty grows. Phase 1 stops at the first design where every
violation is within `feasibility`, and the main run starts at that design. Phase 1 can
instead converge or stall with a violation above `feasibility`; it stalls when its least
infeasibility falls by less than 1% over 100 evaluations, an engineering default. The
run then ends there with `outcome` `:infeasible`, or `:infeasible_box` when parameters
are at a bound, and no outer iteration. That design is the least infeasible one the
search reached from this start, a local verdict, and its `constraints` rows show the
violation of each constraint at that design.

`relaxation = true` implies `phase1`. After such a verdict, it adds one constrained run
per constraint violated or active at that design. Each run computes the smallest largest
violation
of the components of its constraint with every other constraint met, which is how far
that one constraint must be relaxed ([`Relaxation`](@ref), in
`PointingResult.relaxations`). These runs share the budget and the time limit.

`phase1` is on by default and `relaxation` off. A run whose start violates a constraint
starts its main run from the design of phase 1, and it can end at another design of a
flat family than a run started at `x0`. `phase1 = false` starts the main run at `x0`. A
start at which every constraint is met, and a restart with `from`, skip phase 1. The
`outer` of the callback is 0 during phase 1. Phase 1's L-BFGS-B run is not capped by
`iterations`, which caps the subproblems of the augmented Lagrangian. It ends on its
convergence, its stall rule, the `evaluations` of the run, the time limit or the
callback. A stop before a verdict is named in the result (`phase1.reason`) and its
display.

On the dropped-bundle test case made infeasible, phase 1 reached the verdict in 577
evaluations, where the augmented Lagrangian had none after 2,195. On feasible starts,
phase 1 cost 8% fewer to 29% more evaluations (research notes, "A verdict on
conflicting constraints"). On the constrained cases of the test suite that start
outside their limits, it cost −23% to +54%, with drive designs 3e-4 apart on a flat
family (research notes, "Defaults changed (2026-11-19)").

**VisRad's hard cut.** A scenario can deposit with a [`VisRadProfile`](@ref) at its
hard cut (`taper = 0`, its default). Every deposition-derived quantity then jumps where
an operator point crosses the cut, and line searches stall on the jumps. With
`smooth_cut = true` (the default), the run optimizes the same scenario with the profile
tapered over one e-fold below the cut (`taper = 1`). The points and screens of the
operator and the responses of the scenario are shared, and nothing is rebuilt. The run
then reports the design under the hard cut of the scenario. The `value`,
`initial_value`, `constraints` and `outcome` of the result are those of the hard cut at
the design and start of the run, and `tapered_value` contains the objective the tapered
run reached.

The taper moves wall elements by at most 7e-5 of the wall peak on NIF hohlraums, well
inside the default `feasibility` of 1e-3. A design within its limits under the taper is
therefore normally within them under the cut. Where it is not, its row shows `:violated`
and `outcome` is `:infeasible`. The multipliers (fitted under the taper), stationarity,
bound sensitivities, phase 1, relaxations and `state` are those of the tapered run, and
a restart with `from` continues it. `smooth_cut = false` optimizes the scenario profile
as it is. Other scenarios, such as the f-cone model or a profile already tapered, are
optimized as they are. Each phase of a [`PhasedScenario`](@ref) is tapered separately.

**Restarts.** `from`, the [`PointingResult`](@ref) of a previous run, continues that
run. `x0`, `lower` and `upper` default to its last iterate and box. The last iterate is
`state.x` of a constrained run, and `x` otherwise or whenever the run returned its last
iterate. A constrained run starts from the multipliers, penalty parameter, inner
tolerance and working limits of the previous run (`result.state`) instead of starting
them again from zero. A run cut short by `evaluations` or by the callback and continued
with `from` resumes the interrupted inner problem. The counts of the new result
(`iterations`, `evaluations`, `outer`) cover only that run, and its `initial_value` is
the objective at its start. The constraints must have the same components as those of
the previous run, or a `DimensionMismatch` is thrown. An unconstrained `from` result
only sets the start and the box.

Compared with penalty terms at fixed weights, the weights of the constraints are carried
by the multipliers from one inner run to the next. The limits are met without stiff weights,
which stall the line search, and without stages of loosened limits set by hand. The
design also ends at its limits. With penalty terms, it ends where the terms balance the
objective. On the dropped-bundle tutorial, this gave 0.35% more drive in a comparable
number of evaluations, 997 against 1,072 for the penalty form (research notes,
"Constrained re-pointing").
"""
Base.@nospecializeinfer function optimize_pointing(
                           s::_AbstractScenario, @nospecialize(objective);
                           @nospecialize(from = nothing),
                           x0 = from === nothing ? initial_parameters(s.param) : _restart_x(from),
                           lower = from === nothing ? bounds(s.param, x0)[1] : from.lower,
                           upper = from === nothing ? bounds(s.param, x0)[2] : from.upper,
                           iterations::Integer = 200, evaluations::Integer = 20 * iterations,
                           show_trace::Bool = false, g_tol = 1e-8, @nospecialize(callback = nothing),
                           @nospecialize(constraints = ()), feasibility::Real = 1e-3,
                           outer::Integer = 30, time_limit::Real = Inf, phase1::Bool = true,
                           relaxation::Bool = false, inner_reltol::Real = 1e-10,
                           @nospecialize(inner_evaluations = nothing), smooth_cut::Bool = true)
    # This body is inferred and compiled once for every objective and callback. Each
    # objective closure type used to recompile it and the L-BFGS-B of Optim, at about 1 s
    # each. The objective is reached by dynamic dispatch once per evaluation.
    deadline = time() + Float64(time_limit)
    time_limit > 0 || throw(ArgumentError("time_limit must be positive (seconds)"))
    # Under the VisRad hard cut, the run uses the profile tapered below the cut, and the
    # result is reported under the cut of the scenario.
    st = smooth_cut ? _tapered(s) : nothing
    if st !== nothing
        r = optimize_pointing(st, objective; from, x0, lower, upper, iterations, evaluations,
                              show_trace, g_tol, callback, constraints, feasibility, outer,
                              time_limit, phase1, relaxation, inner_reltol, inner_evaluations,
                              smooth_cut = false)
        return _hard_cut_result(s, r, objective, constraints, clamp.(x0, lower, upper),
                                Float64(feasibility))
    end
    from === nothing || from isa PointingResult ||
        throw(ArgumentError("`from` takes a PointingResult, got $(typeof(from))"))
    xs = clamp.(x0, lower, upper)
    inner_evaluations === nothing || inner_evaluations isa Integer ||
        throw(ArgumentError("inner_evaluations takes an integer or `nothing`"))
    inner_evaluations === nothing || !isempty(constraints) ||
        throw(ArgumentError("inner_evaluations caps the inner runs of a constrained run; without " *
                            "constraints `evaluations` is the run's budget"))
    isempty(constraints) ||
        return _constrained_pointing(s, objective, Tuple(constraints), x0, xs, lower, upper;
                                     iterations, evaluations, show_trace, g_tol, callback,
                                     feasibility = Float64(feasibility), outer = Int(outer),
                                     state = from === nothing ? nothing : from.state, deadline,
                                     phase1, relaxation, inner_reltol = Float64(inner_reltol),
                                     inner_evaluations = inner_evaluations === nothing ? typemax(Int) :
                                                         Int(inner_evaluations))
    r = _lbfgsb(s, objective, x0, xs, lower, upper; iterations, evaluations, show_trace, g_tol,
                callback, deadline)
    return _pointing_result(s, r.x, r.value, r.initial_value, r.iterations, r.evaluations,
                            r.stopped, r.reason, r.result, lower, upper; objective)
end

# Returns the start of a restart. It is the iterate that the multipliers in the `state` of
# the result belong to (`state.x`), or else the `x` of the result. The two differ only
# when the run ended above its limits and returned its iterate of least violation.
# The return type is declared `Vector{Float64}` because `from` is not specialized on. An
# inferred `Any` would compile the body of `optimize_pointing` a second time.
function _restart_x(r)::Vector{Float64}
    st = r.state
    return st isa NamedTuple && haskey(st, :x) ? st.x : r.x
end

# Returns the gradient of `F(evaluate(s, x))` at `x`, through the dual tag and chunks of
# the optimizer. It is identical to the chunk-mode gradient of ForwardDiff (see
# `_pointing_gradient!`).
Base.@nospecializeinfer function _scenario_gradient(s::_AbstractScenario, @nospecialize(F), x::Vector{Float64})
    J(x) = F(evaluate(s, x))
    cfg = ForwardDiff.GradientConfig(J, x, ForwardDiff.Chunk(min(length(x), 12)),
                                     ForwardDiff.Tag(_PointingTag(), Float64))
    out = ForwardDiff.DiffResults.GradientResult(similar(x))
    return copy(ForwardDiff.DiffResults.gradient(_pointing_gradient!(out, J, x, cfg)))
end

# Returns the Jacobian of the vector `F(evaluate(s, x))` at `x`, with one row per entry,
# through the dual tag and chunks of the optimizer. The dual pass of each chunk computes
# every entry, so the Jacobian costs one gradient whatever its number of rows. The passes
# run concurrently as in `_pointing_gradient!`. Each row is identical to the chunk-mode
# gradient of that entry from ForwardDiff.
Base.@nospecializeinfer function _scenario_jacobian(s::_AbstractScenario, @nospecialize(F), x::Vector{Float64})
    J(x) = F(evaluate(s, x))
    cfg = ForwardDiff.GradientConfig(J, x, ForwardDiff.Chunk(min(length(x), 12)),
                                     ForwardDiff.Tag(_PointingTag(), Float64))
    return _pointing_jacobian(J, x, cfg)
end

function _pointing_jacobian(f, x::Vector{Float64}, cfg::ForwardDiff.GradientConfig{T,Float64,N}) where {T,N}
    nc = cld(length(x), N)
    blocks = if nc <= 1 || Threads.nthreads() == 1 || !_PARALLEL_CHUNKS[]
        [_jacobian_chunk(f, x, cfg.seeds, T, c) for c in 1:nc]
    else
        tasks = [Threads.@spawn _jacobian_chunk(f, x, cfg.seeds, T, c) for c in 1:nc]
        map(tasks) do t
            try
                fetch(t)
            catch err
                err isa TaskFailedException ? throw(err.task.exception) : rethrow()
            end
        end
    end
    return reduce(hcat, blocks)
end

# Runs dual pass `c` of `_pointing_jacobian`, which computes the columns of its
# parameters.
function _jacobian_chunk(f, x, seeds::NTuple{N,ForwardDiff.Partials{N,Float64}}, ::Type{T}, c) where {N,T}
    lo = (c - 1) * N
    zp = zero(ForwardDiff.Partials{N,Float64})
    xd = [ForwardDiff.Dual{T}(x[i], lo < i <= lo + N ? seeds[i-lo] : zp) for i in eachindex(x)]
    y = f(xd)
    k = min(N, length(x) - lo)
    return Float64[y[r] isa ForwardDiff.Dual ? ForwardDiff.partials(y[r], j) : 0.0
                   for r in eachindex(y), j in 1:k]
end

# Runs L-BFGS-B once on `objective(evaluate(s, x))` in the box, from `xs`; `x0` sizes the
# gradient configuration. The run returns the parameters, objective value, iterations,
# evaluations and stop of `optimize_pointing`, and the result of Optim, which is
# `nothing` when stopped. `f_reltol`, an Optim option that is off at its default of 0,
# ends the run on a small relative change of the value. The inner runs of the
# constrained path use it.
Base.@nospecializeinfer function _lbfgsb(s::_AbstractScenario, @nospecialize(objective), x0, xs, lower,
                                         upper; iterations::Integer, evaluations::Integer,
                                         show_trace::Bool, g_tol, @nospecialize(callback),
                                         f_reltol::Float64 = 0.0, deadline::Float64 = Inf,
                                         cap::Int = typemax(Int))
    best_x, best_value, count, iter = copy(xs), Inf, Ref(0), Ref(0)
    J(x) = objective(evaluate(s, x))
    # Every objective value passed to the optimizer goes through `seen` for the
    # bookkeeping, the callback and the budget.
    function seen(x, v)
        count[] += 1
        v < best_value && ((best_value, best_x) = (v, copy(x)))
        callback === nothing ||
            callback((; evaluations = count[], iteration = iter[], value = v, best_value,
                      best_x)) === true && throw(_OptimizationStop(:callback))
        count[] >= evaluations && throw(_OptimizationStop(:evaluations))
        time() > deadline && throw(_OptimizationStop(:time_limit))
        count[] >= cap && throw(_OptimizationStop(:cap))
        return v
    end
    f(x) = seen(x, J(x))
    # One dual tag is used for every objective, so the compiled deposition and objective
    # code, including the precompiled copy in the package, works for any objective. The
    # tag check of ForwardDiff (`Val(false)` in `_pointing_gradient!`) compares the
    # function of the tag with `J` and would refuse this shared tag. That check guards
    # against a config reused for another function; nesting is a separate matter. An
    # objective that differentiates internally makes a separate tag. The ordering of
    # `_PointingTag` (`≺`) nests that tag inside this one, whatever the tag counts
    # (test/nested_gradient.jl).
    cfg = ForwardDiff.GradientConfig(J, x0, ForwardDiff.Chunk(min(length(x0), 12)),
                                     ForwardDiff.Tag(_PointingTag(), Float64))
    # The dual passes of the chunks run concurrently (`_pointing_gradient!`), and the
    # result is identical to the serial chunk-mode gradient of ForwardDiff.
    out = ForwardDiff.DiffResults.GradientResult(similar(xs))
    g!(G, x) = copyto!(G, ForwardDiff.DiffResults.gradient(_pointing_gradient!(out, J, x, cfg)))
    # L-BFGS-B evaluates the value and gradient together, and the dual pass computes both.
    function fg!(G, x)
        _pointing_gradient!(out, J, x, cfg)
        copyto!(G, ForwardDiff.DiffResults.gradient(out))
        return seen(x, ForwardDiff.DiffResults.value(out))
    end
    initial_value = J(xs)
    track(st) = (iter[] = st.iteration; false)
    res = try
        Optim.optimize(Optim.OnceDifferentiable(f, g!, fg!, xs, 0.0), lower, upper, xs,
                       Optim.LBFGSB(),
                       Optim.Options(; iterations, show_trace, g_abstol = g_tol, f_reltol,
                                     callback = track))
    catch err
        err isa _OptimizationStop || rethrow()
        err
    end
    if res isa _OptimizationStop
        # The cap on the evaluations of one run ends that run as the iteration limit does.
        # It does not stop the whole optimization.
        capped = res.reason === :cap
        return (; x = best_x, value = best_value, initial_value, iterations = iter[],
                evaluations = count[], stopped = !capped, reason = capped ? :iterations : res.reason,
                result = nothing)
    end
    n = Optim.iterations(res)
    return (; x = Optim.minimizer(res), value = Optim.minimum(res), initial_value,
            iterations = n, evaluations = count[], stopped = false,
            reason = Optim.converged(res) || n < iterations ? :converged : :iterations,
            result = res)
end

"""
    ConstraintRow

One component of a [`Constraint`](@ref) at a design, from [`check_constraints`](@ref) or
in `PointingResult.constraints`. Fields:

- `name` — the name of the constraint, with the label of the component for a vector
  constraint, e.g. `"mode_gate[l = 2]"`.
- `value`, `limit` — in the units of the constraint.
- `upper` — `true` for a ceiling `value ≤ limit`, `false` for a floor.
- `violation` — the normalized excess `±(value - limit)/scale`, positive when violated
  and negative inside.
- `multiplier` — the Lagrange multiplier of the component at the design, in the
  normalized units. In the result of a constrained run, it is the nonnegative
  least-squares solution of the stationarity condition on the parameters inside the box
  ([`PointingResult`](@ref)). The fit covers the components within `feasibility` of
  their limits or beyond them, and the other components get zero. It is 0 when not
  computed, as in [`check_constraints`](@ref).
- `status` — `:violated` (`violation > feasibility`), `:active` (within `feasibility`
  of the limit) or `:satisfied`.
"""
struct ConstraintRow
    name::String
    value::Float64
    limit::Float64
    upper::Bool
    violation::Float64
    multiplier::Float64
    status::Symbol
end

"""
    Relaxation

How far a single constraint must be relaxed for the others to be met. A relaxation is
an entry of `PointingResult.relaxations` after a run that [`optimize_pointing`](@ref)
found infeasible with `relaxation = true`. Fields:

- `name` — the name of the constraint.
- `violation` — the largest normalized violation of its components at a design where
  every other constraint is met. This is the relaxation in units of its `scale`; where
  the scale is the limit (the default), the relaxed limit is `1 + violation` times the
  limit.
- `limit`, `relaxed_limit` — one per component, in the units of the constraint. The
  relaxed limit is a ceiling raised, or a floor lowered, by `violation × scale`.
- `upper` — `true` for a ceiling.
- `status` — `:relaxed`, `:others_conflict` or `:stopped` (below).
- `x` — that design.
- `constraints` — the [`ConstraintRow`](@ref) of every constraint at `x`.
- `evaluations` — the objective evaluations of the search.

`status` is `:relaxed` when the others are met at `x`. It is `:others_conflict` when the
other constraints cannot all be met either, so relaxing only this one does not
suffice; `violation` is then `NaN`. It is `:stopped` when the budget, the time limit or
the callback ended the search first; `violation` is then that of the best design, and
the others are possibly unmet.

The relaxation is achievable, since the others are met at the design `x` with this one at
`relaxed_limit`. It is the least one only locally. The run that finds it minimizes a
smooth stand-in for the largest violation, a 64-norm of the positive violations, from
the infeasible design, which is a local search.
"""
struct Relaxation
    name::String
    violation::Float64
    limit::Vector{Float64}
    relaxed_limit::Vector{Float64}
    upper::Bool
    status::Symbol
    x::Vector{Float64}
    constraints::Vector{ConstraintRow}
    evaluations::Int
end

"""
    PointingResult

Result of [`optimize_pointing`](@ref). Fields:

- `x` — the best parameters.
- `beams` — `apply(scenario.param, x)`, a [`PerPhase`](@ref) of beam vectors for a
  phased scenario.
- `value`, `initial_value` — the objective at `x` and at the clamped `x0`.
- `iterations` — the L-BFGS-B iterations completed.
- `evaluations` — the objective evaluations, each gradient counting as one.
- `reason` — `:converged`, `:iterations`, `:evaluations`, `:time_limit` or `:callback`.
- `stopped` — `true` for the last three reasons, which return the best design evaluated
  before the stop.
- `result` — the result of Optim, `nothing` when `stopped`.
- `names` — the `parameter_names` of the parameterization, or `x[i]` for one without a
  method.
- `lower`, `upper` — the box.
- `at_bounds` — the names of the parameters at a bound, within 1e-8 of the box width.
  Parameters fixed by `lower == upper` are excluded.

The display lists the parameters at a bound with their values and bounds.

What ended the run is recorded in two more fields. `constraints` holds one
[`ConstraintRow`](@ref) per constraint component at `x`, and it is empty without
`constraints`. The display shows a constraint of more than eight components, such as a
keep-out per beam, as its worst row and a count by status. `outcome` is one of:

- `:interior`: the run converged with every constraint met and no parameter at a bound;
- `:box`: the run converged with every constraint met and parameters at a bound. The box
  is active, and a wider box may lower the objective;
- `:infeasible_box`: a constraint is still violated and parameters are at a bound, so the
  constraint cannot be met within the box;
- `:infeasible`: a constraint is still violated and no parameter is at a bound. The
  design is a local point where the constraints cannot be met, and another start may
  help;
- `:iterations`, `:evaluations`, `:time_limit` or `:callback`: the budget, the time
  limit or the callback ended the run first. Its constraints, if any, are in
  `constraints`.

`reason` is `:converged` for the first four. `outer` is the number of outer
(multiplier) iterations of a constrained run, and 0 without constraints. `iterations`
and `evaluations` are then totals over them, and `value` and `initial_value` are the
values of the objective without the constraint terms.

`state` is what `optimize_pointing(...; from = result)` continues from. It is `nothing`
without constraints. For a constrained run, it is a named tuple of:

- the multipliers `λ`, in normalized units, one per constraint component in the order
  of `constraints`;
- the penalty parameter `μ` and its ceiling `μmax`;
- the tolerance `tol` of the inner runs;
- the last infeasibility–complementarity measure `Vold`;
- the working limits `limits`, which contain the lowered stand-in limit of an exact
  ceiling (see [`irradiance_ceiling`](@ref));
- the iterate `x` they belong to. It is the `x` of the result unless the run returned
  its iterate of least violation (below).

For a run ended by its budget, time limit or callback, the state of the interrupted
inner problem is kept, so that the restart resumes it.

**Multipliers.** The multipliers in `constraints` are those that best satisfy
stationarity at `x`, which differ from the iterated `state.λ`. They are the nonnegative
least-squares solution of `∂f/∂xᵢ + Σⱼ λⱼ ∂gⱼ/∂xᵢ = 0` over the parameters inside the
box. Here `gⱼ` are the normalized violations, with the stand-in of an exact ceiling at
its working limit. The fit covers the components within `feasibility` of their limits
or beyond them, and the other components get zero (complementarity). The iterated
estimate `λ ← max(0, λ + μ g)` after each inner run is exact only when the inner run
ended at a stationary point. The inner runs of a drive objective end on their value
tolerance instead. At a large penalty `μ`, the update then zeroes the multipliers of
active limits at which the design ends a few 1e-5 inside (research notes, "Multipliers
from the converged active set").

The fit costs one Jacobian of the objective and the constraints after the run, the
cost of a gradient, which is not counted in `evaluations`. Several sets of multipliers
can fit equally, as for components with the same gradient, such as the clearances of
beams related by the port symmetry. The fit then returns the one of least norm.
`stationarity` is the residual of the fit: the largest `|∂L/∂xᵢ|` over the parameters
inside the box, with `L = f + Σⱼ λⱼ gⱼ`. It is in objective units per parameter unit,
the measure that the inner runs compare with `g_tol`, and it is `NaN` without
constraints or when phase 1 ended the run. Near zero, the design is a
Karush–Kuhn–Tucker point with these multipliers. A drive run whose inner runs ended on
their value tolerance can end with a stationarity that is a sizable fraction of the
objective gradient. The multipliers and bound sensitivities are then estimates at a
design short of stationarity. On the case of the dropped-bundle tutorial, the
stationarity is 1.2e-3, against 3.0e-3 per mm for the largest derivative of the drive
objective.

An exact ceiling, such as [`irradiance_ceiling`](@ref) or [`picket_ceiling`](@ref) with
`exact = true`, reports the exact value in its row. It is included in the fit through
its smooth stand-in at the working limit. Its row can therefore show `:satisfied`, with
the exact value inside the limit by more than `feasibility`, together with a positive
multiplier from the active stand-in. The display notes such rows; the multiplier is
that of the stand-in.

**Bound sensitivities.** For each parameter in `at_bounds`, in the same order,
`bound_gradient` contains the derivative `∂f/∂xᵢ` of the objective at `x`.
`bound_sensitivity` contains that of the Lagrangian, `∂L/∂xᵢ` with `L = f + Σⱼ λⱼ gⱼ`,
using the multipliers of `constraints` above. Without constraints the two are equal. By
the envelope theorem, `bound_sensitivity` is the derivative of the optimal objective
with respect to the bound value, with the active set held. When the run has converged,
it is positive at a lower bound and negative at an upper bound. Widening the box by δ in
that parameter, by lowering a lower bound or raising an upper one, then lowers the
optimum by about `|bound_sensitivity| δ`, to first order. The unit is objective units
per parameter unit. `bound_gradient` is the gain of the objective without the constraint
terms, before the constraints take their share.

Without constraints, the sensitivities cost one gradient after the run. With
constraints, they come from the Jacobian of the multiplier fit. Neither cost is counted
in `evaluations`. The sensitivities are meaningful at convergence. The active
constraints and bounds can leave a flat family of designs, as for a drive objective
under several active gates, and the run can end short of stationarity when its inner
runs end on their value tolerance. In both cases the derivative can differ from the
gain of a wider box by its whole size. A sign that would move the parameter inward
marks a bound that does not hold it (research notes, "Multipliers from the converged
active set"). A restart with `from` in the wider box measures the gain.

**VisRad's hard cut.** When `optimize_pointing` optimized a scenario with the hard cut
on a profile tapered below its cut (`smooth_cut`, its default), `tapered_value` is the
objective that the tapered run reached. `value`, `initial_value`, `constraints` and
`outcome` are then those of the scenario under its hard cut. Otherwise `tapered_value`
is `NaN`.

**Phase 1 and relaxations** (`phase1 = true`, `relaxation = true`). `phase1` is
`nothing` when phase 1 did not run. Otherwise it is a named tuple of:

- `evaluations` and `iterations`, which are included in the totals;
- `feasible`, `true` when phase 1 reached a design within the limits, where the main run
  started;
- `verdict`, `true` when phase 1 converged or stalled above `feasibility`. The run then
  ended at that design, with `outer` 0;
- its design `x`;
- `reason`, the cause of the end of phase 1. It is `:feasible`, `:converged` or
  `:stalled` for a verdict or a design within the limits. Otherwise it is the stop that
  came first, `:evaluations`, `:time_limit` or `:callback`. No verdict is reached on such
  a stop, and the run ends with its design of least violation.

With `relaxation = true`, `relaxations` contains one [`Relaxation`](@ref) per constraint
violated or active at the design of a verdict; otherwise it is empty. Their evaluations
count in `evaluations`. The display lists them under the constraints.

**Least violation and the last design within the limits.** These fields apply to
constrained runs and are `nothing` without constraints. The run records the largest
normalized violation of each iterate over the reported rows, `max(0, maxⱼ violation)`
of [`ConstraintRow`](@ref). For an exact ceiling, the row contains the exact value
against its limit, and the violation is zero when the constraint is met for every
component. The record is taken
at the start, at the design of phase 1 and after every inner run. `least_violation` is
the iterate of least violation, with violations within `feasibility` counted as equal
and the lower objective then preferred. `last_feasible` is the latest iterate with every
violation within `feasibility`, or `nothing` when none had. Each is a named tuple
`(; x, violation, value, outer)` that contains the parameters, the largest violation, the
objective and the outer iteration of that iterate. The outer iteration is 0 for the
start and for the design of phase 1.

A run can end with a violation above `feasibility`: with `outcome` `:infeasible` or
`:infeasible_box`, or after a budget, time-limit or callback stop above the limits. It
then returns `least_violation` as `x`, with its beams, value and rows. The last iterate
is not returned, since it can be worse than the start on the limits that the start
met. `state.x` is then the last iterate, to which the multipliers of `state` belong and
from which a restart with `from` continues. When an earlier iterate met the limits, the
limits are met at the returned design, and a run that ended at the penalty ceiling has
`reason` and `outcome` `:iterations`. Under VisRad's hard cut (`tapered_value` finite), the
records are those of the tapered run.

`capped` counts the outer iterations whose inner run ended on its cap (`iterations` or
`inner_evaluations`) instead of converging. After each such run that ended above the
limits, the multipliers were updated, and the penalty, the inner tolerance and the
reference of the penalty test were held ([`optimize_pointing`](@ref)). The count is
shown in the display.
"""
struct PointingResult{B}
    x::Vector{Float64}
    beams::B
    value::Float64
    initial_value::Float64
    iterations::Int
    evaluations::Int
    stopped::Bool
    reason::Symbol
    result::Any
    names::Vector{String}
    lower::Vector{Float64}
    upper::Vector{Float64}
    at_bounds::Vector{String}
    constraints::Vector{ConstraintRow}
    outcome::Symbol
    outer::Int
    state::Union{Nothing,NamedTuple}
    bound_gradient::Vector{Float64}
    bound_sensitivity::Vector{Float64}
    phase1::Union{Nothing,NamedTuple}
    relaxations::Vector{Relaxation}
    stationarity::Float64
    tapered_value::Float64
    least_violation::Union{Nothing,NamedTuple}
    last_feasible::Union{Nothing,NamedTuple}
    capped::Int
end


# Returns 0 inside the box, -1 at the lower bound and 1 at the upper bound, within 1e-8
# of the width. A parameter fixed by lower == upper counts as inside.
function _bound_side(x, lo, hi)
    lo < hi || return 0
    tol = isfinite(hi - lo) ? 1e-8 * (hi - lo) : 0.0
    return x - lo <= tol ? -1 : hi - x <= tol ? 1 : 0
end

# `x`, the values and the Optim `result` are left unspecialized. The argument types of the
# calls from the several exits of the constrained run are inferred differently, and the
# calls share one body.
Base.@nospecializeinfer function _pointing_result(
        s::_AbstractScenario, @nospecialize(x), @nospecialize(value), @nospecialize(initial_value),
        iterations, evaluations, stopped, reason, @nospecialize(result),
        lower, upper; constraints::Vector{ConstraintRow} = ConstraintRow[], outer::Int = 0,
        @nospecialize(state = nothing), @nospecialize(objective = nothing),
        gradients::Union{Nothing,NTuple{2,Vector{Float64}}} = nothing, stationarity::Float64 = NaN,
        @nospecialize(phase1 = nothing), relaxations::Vector{Relaxation} = Relaxation[],
        @nospecialize(least_violation = nothing), @nospecialize(last_feasible = nothing), capped::Int = 0)
    x = collect(Float64, x)
    names = _parameter_names(s.param, length(x))
    side = map(_bound_side, x, lower, upper)
    atb = names[side .!= 0]
    outcome = _outcome(reason, constraints, atb)
    # The derivatives at the bound parameters are computed only when a parameter is at a
    # bound. A constrained run supplies its gradients of the objective and the Lagrangian
    # at `x`; otherwise one gradient of the objective is computed.
    gf, gL = Float64[], Float64[]
    at = findall(!=(0), side)
    if gradients !== nothing
        gf, gL = gradients[1][at], gradients[2][at]
    elseif !isempty(atb) && objective !== nothing
        gf = _scenario_gradient(s, objective, x)[at]
        gL = copy(gf)
    end
    return PointingResult(x, apply(s.param, x), Float64(value), Float64(initial_value),
                          Int(iterations), Int(evaluations), stopped, reason, result, names,
                          collect(Float64, lower), collect(Float64, upper), atb,
                          constraints, outcome, outer, state, gf, gL, phase1, relaxations, stationarity,
                          NaN, least_violation, last_feasible, capped)
end

# Returns what ended a run. This is the reason of the stop or, at convergence, whether the
# constraints are met and whether parameters are at a bound (`atb` holds their names).
function _outcome(reason::Symbol, rows::Vector{ConstraintRow}, atb)
    reason !== :converged && return reason
    met = all(r -> r.status !== :violated, rows)
    return met ? (isempty(atb) ? :interior : :box) : (isempty(atb) ? :infeasible : :infeasible_box)
end

# Returns the scenario that the optimizer runs on for a scenario with the VisRad profile
# at its hard cut. The operator is the same, with its points, normals, visibility, work
# plan and envelope ratios, but the cut is tapered over `_OPT_TAPER` e-folds. The
# responses and caches are shared. For any other scenario, the result is `nothing`. The
# taper does not change the support of the profile, so the operator points and the
# screens that skip points outside the support are valid for both.
const _OPT_TAPER = 1.0

function _tapered(op::DepositionOperator)
    prof = op.profile
    (prof === nothing || prof.taper > 0) && return nothing
    return DepositionOperator(op.target, op.ports, op.dirs, op.pts, op.w, op.normal, op.mask, op.dest,
                              op.plan, VisRadProfile(prof.envelope, _OPT_TAPER), op.envelope)
end

_tapered(::_AbstractScenario) = nothing

function _tapered(s::Scenario)
    op = _tapered(s.deposition)
    op === nothing && return nothing
    return Scenario(s.target, s.param, op, s.model, s.observed, s.moment_weights, s.ylm_weights,
                    s.moment_response, s.ylm_response, s.lpi, s.time, getfield(s, :_response), s.albedo,
                    getfield(s, :_split))
end

# Some constraints read the element powers of sets of beam groups (`_split_sets`), such as
# the outer cones of a picket ceiling. For such constraints, the optimizer runs on the
# same scenario with those sets recorded. Its evaluations then keep the powers of the sets
# from the deposition of the evaluation (`evaluate`), and the constraints read them
# instead of depositing the groups again. Each entry of `req` pairs a phase with a group
# set, and the phase `nothing` stands for the scenario itself. `_with_split` returns `s`
# unchanged when no set applies.
const _SplitRequest = Vector{Pair{Union{Nothing,Symbol},Vector{Symbol}}}

_split_sets(f) = _SplitRequest()

function _split_request(cs::Tuple)
    req = _SplitRequest()
    for c in cs
        append!(req, _split_sets(c.f))
    end
    return req
end

_with_split(s::_AbstractScenario, req) = s

function _with_split(s::Scenario, req)
    old = getfield(s, :_split)
    sets = copy(old)
    for (k, g) in req
        k === nothing && !(g in sets) && push!(sets, g)
    end
    length(sets) == length(old) && return s
    return Scenario(s.target, s.param, s.deposition, s.model, s.observed, s.moment_weights, s.ylm_weights,
                    s.moment_response, s.ylm_response, s.lpi, s.time, getfield(s, :_response), s.albedo, sets)
end

# Reports a run on the tapered profile (`_tapered`) under the hard cut of the scenario `s`.
# The objective is evaluated at the design of the run and at its start `xs`, with the
# constraint rows there and the outcome. The rows contain the multipliers of the run. The
# value of the tapered run is kept as `tapered_value`, and every other field is taken
# from the run.
Base.@nospecializeinfer function _hard_cut_result(s::_AbstractScenario, r::PointingResult,
                                                  @nospecialize(objective), @nospecialize(constraints),
                                                  xs::AbstractVector, feasibility::Float64)
    ev = evaluate(s, r.x)
    value = Float64(ForwardDiff.value(objective(ev)))
    initial = Float64(ForwardDiff.value(objective(evaluate(s, collect(Float64, xs)))))
    rows = isempty(r.constraints) ? r.constraints :
           _constraint_rows(ev, Tuple(constraints), [c.multiplier for c in r.constraints], feasibility)
    return PointingResult(r.x, r.beams, value, initial, r.iterations, r.evaluations, r.stopped,
                          r.reason, r.result, r.names, r.lower, r.upper, r.at_bounds, rows,
                          _outcome(r.reason, rows, r.at_bounds), r.outer, r.state, r.bound_gradient,
                          r.bound_sensitivity, r.phase1, r.relaxations, r.stationarity, r.value,
                          r.least_violation, r.last_feasible, r.capped)
end

function Base.show(io::IO, r::PointingResult)
    print(io, "PointingResult(", r.reason, ", objective ", round(r.initial_value; sigdigits = 4),
          " → ", round(r.value; sigdigits = 4), ", ", length(r.at_bounds), " of ",
          length(r.x), " parameters at bounds")
    isempty(r.constraints) ||
        print(io, ", ", count(c -> c.status === :violated, r.constraints), " of ",
              length(r.constraints), " constraints violated")
    print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", r::PointingResult)
    how = r.outcome in (:infeasible, :infeasible_box) ? "ended with constraints violated" :
          r.reason === :converged ? "converged" :
          r.reason === :iterations ? "reached the iteration limit" :
          r.reason === :evaluations ? "reached the evaluation limit (best design kept)" :
          r.reason === :time_limit ? "reached the time limit (best design kept)" :
          "stopped by the callback (best design kept)"
    @printf(io, "Pointing optimization %s: objective %.4g → %.4g\n", how, r.initial_value,
            r.value)
    isnan(r.tapered_value) ||
        @printf(io, "  optimized on VisRad's profile tapered below its cut (objective %.4g there); reported under the scenario's hard cut\n",
                r.tapered_value)
    @printf(io, "  %d iterations, %d evaluations, %d parameters", r.iterations,
            r.evaluations, length(r.x))
    if isempty(r.at_bounds)
        print(io, ", none at a bound")
    else
        print(io, "; at a bound (widen the box if the objective should set them):")
        w = maximum(length, r.at_bounds)
        for (j, n) in enumerate(r.at_bounds)
            i = findfirst(==(n), r.names)
            lo = _bound_side(r.x[i], r.lower[i], r.upper[i]) < 0
            @printf(io, "\n    %s = %.6g  (%s bound; box %.6g to %.6g)", rpad(n, w), r.x[i],
                    lo ? "lower" : "upper", r.lower[i], r.upper[i])
            j <= length(r.bound_sensitivity) &&
                @printf(io, "  d objective/d bound %.3g", r.bound_sensitivity[j])
        end
    end
    isempty(r.constraints) && return nothing
    what = r.outcome === :interior ? "constraints met inside the box" :
           r.outcome === :box ? "constraints met, parameters at a bound" :
           r.outcome === :infeasible_box ? "constraints violated with parameters at a bound (the box blocks them)" :
           r.outcome === :infeasible ? "constraints violated inside the box" :
           "ended by the budget, the time limit or the callback"
    p1 = r.phase1
    if p1 !== nothing && p1.verdict
        @printf(io, "\n  phase 1 (%d evaluations, %s): no design meets every constraint from this start; the least infeasibility found:",
                p1.evaluations, get(p1, :reason, :converged) === :stalled ? "stalled" : "converged")
    else
        if p1 !== nothing
            p1.feasible ? @printf(io, "\n  phase 1: a design within the limits in %d evaluations", p1.evaluations) :
                          @printf(io, "\n  phase 1 (%d evaluations): no design within the limits and no verdict; %s ended it first",
                                  p1.evaluations, _phase1_stop_text(get(p1, :reason, :converged)))
        end
        @printf(io, "\n  %d outer iterations", r.outer)
        r.capped > 0 && @printf(io, " (%d ended on the inner cap)", r.capped)
        print(io, ": ", what)
        isnan(r.stationarity) ||
            @printf(io, "\n  stationarity: largest |∂L/∂x| inside the box %.2e (objective units per parameter unit)",
                    r.stationarity)
    end
    lv, lf = r.least_violation, r.last_feasible
    if lv !== nothing && r.state !== nothing && haskey(r.state, :x) && r.x == lv.x && r.state.x != r.x
        @printf(io, "\n  returned the iterate of least violation (%s, largest violation %.3g), not the last (`state.x`)",
                lv.outer == 0 ? "the start or phase 1's design" : "outer iteration $(lv.outer)", lv.violation)
    end
    if lf === nothing
        r.outcome in (:interior, :box) || print(io, "\n  no iterate met every limit")
    elseif r.x != lf.x
        @printf(io, "\n  the last iterate within the limits: %s, objective %.4g (`last_feasible`)",
                lf.outer == 0 ? "the start or phase 1's design" : "outer iteration $(lf.outer)", lf.value)
    end
    _show_constraint_rows(io, r.constraints)
    # The row of an exact ceiling shows the exact value, inside its limit, while its
    # multiplier is that of the stand-in, which is active at its working limit.
    stand = unique(first(split(x.name, '[')) for x in r.constraints if x.status === :satisfied && x.multiplier > 0)
    isempty(stand) ||
        print(io, "\n  ", join(stand, ", "), ": :satisfied with a positive multiplier, the smooth stand-in's at its ",
              "working limit, not the reported value's (an exact ceiling; PointingResult, \"Multipliers\")")
    isempty(r.relaxations) && return nothing
    print(io, "\n  each violated constraint relaxed alone, the others met:")
    for x in r.relaxations
        print(io, "\n    ", rpad(x.name, 12), " ", _relaxation_text(x))
    end
    return nothing
end

# Describes in words what stopped a phase 1 short of a verdict.
_phase1_stop_text(why::Symbol) = why === :evaluations ? "the evaluation budget" :
                                 why === :time_limit ? "the time limit" : why === :callback ? "the callback" :
                                 string(why)

# Describes a relaxation in words. The relaxed limits are given in the units of the
# constraint, and as a multiple of the limit where the scale is the magnitude of the
# limit.
function _relaxation_text(x::Relaxation)
    x.status === :others_conflict && return "not enough alone: the others conflict too"
    v = x.violation
    s = x.status === :stopped ? " (stopped first: the others possibly unmet)" : ""
    # The relaxation is relative when the scale is the magnitude of the limit, so each
    # limit moves by the same multiple.
    m = x.relaxed_limit ./ x.limit
    rel = all(l -> l != 0, x.limit) &&
          all(i -> isapprox(abs(x.relaxed_limit[i] - x.limit[i]), v * abs(x.limit[i]); rtol = 1e-9),
              eachindex(x.limit))
    lims = length(x.limit) == 1 ? @sprintf("limit %.5g → %.5g", x.limit[1], x.relaxed_limit[1]) :
           "limits " * join((@sprintf("%.4g", l) for l in x.relaxed_limit), ", ")
    return @sprintf("by %.4g × its scale%s: %s%s (%d evaluations)", v,
                    rel ? @sprintf(" (the limit × %.4g)", m[1]) : "", lims, s, x.evaluations)
end

# A display lists all the rows (components) of a constraint up to this many. Beyond it,
# the display shows the worst row, the one of largest violation, and a count by status.
# A keep-out per beam has hundreds of rows, and all of them stay in the vector.
const _ROWS_SHOWN = 8

function _show_constraint_rows(io::IO, rows::AbstractVector{ConstraintRow}; indent = "\n    ")
    # Consecutive rows of one constraint share the name before the label of the
    # component.
    stem(r) = first(split(r.name, '['))
    groups = Vector{ConstraintRow}[]
    for r in rows
        isempty(groups) || stem(groups[end][1]) != stem(r) ? push!(groups, [r]) : push!(groups[end], r)
    end
    shown = [length(g) <= _ROWS_SHOWN ? g : [g[argmax([r.violation for r in g])]] for g in groups]
    w = maximum(r -> length(r.name), Iterators.flatten(shown); init = 10)
    @printf(io, "%s%s %12s %12s %10s %10s  %s", indent, rpad("constraint", w), "value", "limit",
            "violation", "multiplier", "status")
    for (g, rs) in zip(groups, shown), r in rs
        @printf(io, "\n    %s %12.5g %s%11.5g %10.2e %10.3g  %s", rpad(r.name, w), r.value,
                r.upper ? "≤" : "≥", r.limit, r.violation, r.multiplier, r.status)
        length(g) > _ROWS_SHOWN &&
            print(io, " (the worst of ", length(g), ": ", count(x -> x.status === :violated, g), " violated, ",
                  count(x -> x.status === :active, g), " active)")
    end
end

Base.show(io::IO, ::MIME"text/plain", rows::AbstractVector{ConstraintRow}) =
    _show_constraint_rows(io, rows; indent = "    ")

struct _OptimizationStop <: Exception
    reason::Symbol
end

struct _PointingTag end

# The optimizer's dual tag precedes every other tag. ForwardDiff orders the tags of
# nested derivatives by `tagcount`, a counter read when a tag type is first compiled,
# so that a derivative taken inside the objective (its tag `Tag{typeof(g), <:Dual}`)
# nests inside the optimizer's. The precompile workload compiles this tag, so its count
# is the one from the precompile session, while the counter restarts on loading: a tag
# made at run time can then get an equal or smaller count, and its duals would not nest
# in the optimizer's ("Cannot determine ordering of Dual tags"). The gradient of
# `optimize_pointing` is always the outermost derivative (its parameters and scenario
# are Float64, and the tag is not public), so stating the order is exact.
const _PTAG = ForwardDiff.Tag{_PointingTag,Float64}
ForwardDiff.:≺(::Type{_PTAG}, ::Type{ForwardDiff.Tag{F,V}}) where {F,V} = true
ForwardDiff.:≺(::Type{ForwardDiff.Tag{F,V}}, ::Type{_PTAG}) where {F,V} = false
ForwardDiff.:≺(::Type{_PTAG}, ::Type{_PTAG}) = false

# Computes the chunk-mode `gradient!(out, f, x, cfg)` of ForwardDiff with the dual passes
# of the chunks run concurrently, one task each. Their inner threaded loops nest on the
# same pool. Each pass is seeded exactly as ForwardDiff seeds it: the last chunk has N
# partials, with its trailing seeds zero. Every pass is deterministic, so the gradient and
# value are bitwise those of the serial `ForwardDiff.gradient!`. The value is taken from
# the last chunk, as ForwardDiff takes it. The passes run serially on one thread, for a
# single chunk, or with `_PARALLEL_CHUNKS[] = false`, the reference path.
const _PARALLEL_CHUNKS = Ref(true)

function _pointing_gradient!(out, f, x::AbstractVector{Float64},
                             cfg::ForwardDiff.GradientConfig{T,Float64,N}) where {T,N}
    nc = cld(length(x), N)
    (nc <= 1 || Threads.nthreads() == 1 || !_PARALLEL_CHUNKS[]) &&
        return ForwardDiff.gradient!(out, f, x, cfg, Val(false))
    G = ForwardDiff.DiffResults.gradient(out)
    vals = Vector{Float64}(undef, nc)
    tasks = [Threads.@spawn _gradient_chunk!(G, vals, f, x, cfg.seeds, T, c) for c in 1:nc]
    for t in tasks
        try
            wait(t)
        catch err
            # rethrow the exception of the chunk (e.g. a stop request) unwrapped
            err isa TaskFailedException ? throw(err.task.exception) : rethrow()
        end
    end
    return ForwardDiff.DiffResults.value!(out, vals[nc])
end

# Runs dual pass `c` of `_pointing_gradient!`. The pass writes only its gradient slots
# and the value.
function _gradient_chunk!(G, vals, f, x, seeds::NTuple{N,ForwardDiff.Partials{N,Float64}},
                          ::Type{T}, c) where {N,T}
    lo = (c - 1) * N
    zp = zero(ForwardDiff.Partials{N,Float64})
    xd = [ForwardDiff.Dual{T}(x[i], lo < i <= lo + N ? seeds[i-lo] : zp) for i in eachindex(x)]
    y = f(xd)
    for k in 1:min(N, length(x) - lo)
        G[lo+k] = y isa ForwardDiff.Dual ? ForwardDiff.partials(y, k) : 0.0
    end
    vals[c] = ForwardDiff.value(y)
    return nothing
end

"""
    report([io = stdout,] ev; lmax = 6, minor = 1e-3) -> nothing

Print a power accounting and flux-moment summary of an evaluation. The summary lists:

- the setup and effective power per group, when the scenario has a laser–plasma model;
- the delivered laser power and its split over the parts (`component_power`). The parts
  below `minor` × the delivered power are summed on one line;
- the residual, unless it is below that fraction;
- the mean of the observed flux ([`mean_flux`](@ref)), its radiation temperature and
  its Legendre moments `Pₗ/P₀` up to `lmax`.

The residual is the delivered power minus the power absorbed and clipped. It is the
light that does not strike a surface plus the error of the deposition quadrature;
`DepositionOperator` conserves the power of each beam to 1e-4 when nothing clips. The
residual prints as `missed` when positive and as `quadrature excess` when negative,
where the quadrature counts more power than the beams deliver.
"""
function report(io::IO, ev::Evaluation; lmax::Integer = 6, minor::Real = 1e-3)
    d = ev.deposition
    P = ForwardDiff.value(d.delivered)
    if ev.scenario.lpi !== nothing
        Ps = sum(ForwardDiff.value(b.power) for b in ev.beams)
        @printf(io, "Setup laser power: %.4g TW; after the laser–plasma model:\n", Ps)
        for (g, v) in sort(collect(lpi_summary(ev.beams, ev.effective)); by = first)
            @printf(io, "  %-20s %9.4g → %9.4g TW  (×%.3f)\n", g, v.setup, v.effective, v.ratio)
        end
    end
    @printf(io, "Delivered laser power: %.4g TW\n", P)
    cp = component_power(d, ev.scenario.target)
    missed = pop!(cp, :missed)
    rows = sort(collect(cp); by = first)
    small = filter(kv -> abs(kv[2]) < minor * P, rows)
    for (k, v) in rows
        abs(v) < minor * P || @printf(io, "  %-20s %9.4g TW  (%6.2f%%)\n", k, v, 100v / P)
    end
    if !isempty(small)
        v = sum(last, small)
        @printf(io, "  %-20s %9.4g TW  (%6.2f%%)\n",
                string(length(small), length(small) == 1 ? " part" : " parts", " < ",
                       100minor, "%"), v, 100v / P)
    end
    abs(missed) < minor * P ||
        @printf(io, "  %-20s %9.4g TW  (%6.2f%%)\n", missed > 0 ? "missed" : "quadrature excess",
                missed, 100missed / P)
    if _radiative(ev.scenario)
        a = ForwardDiff.value.(flux_moments(ev; lmax))
        @printf(io, "Observed flux: mean %.4g TW/mm² (%.4g W/cm², Tr = %.1f eV)\n", a[1],
                to_Wcm2(a[1]), radiation_temperature(a[1]))
        for l in 1:lmax
            @printf(io, "  P%d/P0 = %+7.3f%%\n", l, 100 * a[l+1] / a[1])
        end
    end
    return nothing
end
report(ev::Evaluation; kwargs...) = report(stdout, ev; kwargs...)
