# This file solves the constrained pointing problem: the augmented Lagrangian loop, the
# multiplier fit at the solution, phase 1 on the infeasibility and the relaxations.

# Concatenate the normalized violations of all constraints at `ev`, as Float64.
_all_violations(cs::Tuple, ev) =
    reduce(vcat, map(c -> Float64.(ForwardDiff.value.(_violations(c, ev))), cs); init = Float64[])

# The augmented Lagrangian (Powell–Hestenes–Rockafellar, for inequalities) of an
# objective with multipliers `λ` and penalty `μ` is
# f + μ/2 Σᵢ [max(0, gᵢ + λᵢ/μ)² - (λᵢ/μ)²], with gᵢ the normalized violations. `λ` contains
# one vector per constraint. The fields are fixed during an inner run, so the value depends
# only on the evaluation, and the function is thread-safe.
struct _Augmented{O,C<:Tuple,L<:Tuple}
    objective::O
    constraints::C
    λ::L
    μ::Float64
end

function (a::_Augmented)(ev)
    J = a.objective(ev)
    terms = map((c, λ) -> _al_term(_violations(c, ev), λ, a.μ), a.constraints, a.λ)
    return J + sum(terms)
end

function _al_term(g::AbstractVector{T}, λ::Vector{Float64}, μ::Float64) where {T}
    length(g) == length(λ) ||
        throw(DimensionMismatch("a constraint changed its number of components"))
    s = zero(T)
    for i in eachindex(g, λ)
        s += max(zero(T), g[i] + λ[i] / μ)^2 - (λ[i] / μ)^2
    end
    return μ / 2 * s
end

# The constrained run of `optimize_pointing` is an augmented Lagrangian. `_lbfgsb` solves
# its inner problems, the box-constrained minimizations of `_Augmented`, each warm-started
# from the last. Multipliers and penalty are updated as in ALGENCAN (Andreani, Birgin,
# Martínez and Schuverdt, SIAM J. Optim. 18, 1286, 2008; Birgin and Martínez, "Practical
# Augmented Lagrangian Methods for Constrained Optimization", SIAM 2014):
# - The multipliers are updated as λ ← max(0, λ + μ g).
# - μ grows tenfold when the infeasibility–complementarity measure
#   V = maxᵢ |max(gᵢ, -λᵢ/μ)| has not halved.
# - The first μ balances the objective against the initial violations.
# - The inner tolerance on the projected gradient tightens tenfold per outer iteration,
#   down to `g_tol`.
# - An inner run also ends when an iteration changes its value by less than
#   `inner_reltol` relative, 1e-10 by default. The drive objectives are flat along the
#   constraints, and without this test an inner run creeps for hundreds of iterations
#   at a 1e-6 gain.
#
# An inner run ended by its cap (`iterations` or `inner_evaluations`) has not solved its
# subproblem, which is the premise of the penalty test (ALGENCAN's step 1, and the inner
# tolerance ω_k of Conn, Gould and Toint, SIAM J. Numer. Anal. 28, 545, 1991). After a
# capped run that ends above the limits, the multipliers are updated as after any inner
# run, but the penalty, the inner tolerance and the reference V of the penalty test are
# held, so that a cap never grows the penalty toward the infeasible verdict. A capped run
# that ends within the limits takes the usual test, and `capped` counts both kinds.
#
# The largest normalized violation of the reported rows is recorded at the start, at
# phase 1's design and after every inner run. A run that ends with a violation above
# `feasibility` returns the iterate of least violation (`least_violation`), whether the
# penalty ceiling, `outer`, the budget, the time limit or the callback stopped it. The
# latest iterate within the limits is kept as `last_feasible`.
Base.@nospecializeinfer function _constrained_pointing(
        s::_AbstractScenario, @nospecialize(objective), @nospecialize(cs::Tuple), x0, xs, lower, upper;
        iterations::Integer, evaluations::Integer, show_trace::Bool, g_tol,
        @nospecialize(callback), feasibility::Float64, outer::Int, @nospecialize(state = nothing),
        deadline::Float64 = Inf, phase1::Bool = false, relaxation::Bool = false,
        inner_reltol::Float64 = 1e-10, inner_evaluations::Int = typemax(Int))
    for c in cs
        c isa Constraint || throw(ArgumentError("constraints must be `Constraint`s, got $(typeof(c))"))
    end
    feasibility > 0 || throw(ArgumentError("feasibility must be positive"))
    inner_reltol >= 0 || throw(ArgumentError("inner_reltol must be nonnegative"))
    inner_evaluations >= 1 || throw(ArgumentError("inner_evaluations must be at least 1"))
    phase1 |= relaxation
    # The picket maps of the constraints are read from the deposition of the evaluations.
    s = _with_split(s, _split_request(cs))
    ev = evaluate(s, xs)
    f0 = Float64(ForwardDiff.value(objective(ev)))
    g = _all_violations(cs, ev)
    lens = map(c -> length(_violations(c, ev)), cs)
    x, ne, ni, k = copy(xs), 0, 0, 0
    p1 = nothing
    # The iterates of least violation and the latest within the limits; outer iteration 0
    # is the start or phase 1's design.
    least = rec0 = _iterate_record(cs, ev, xs, f0, 0, feasibility)
    feas = least.violation <= feasibility ? least : nothing
    rec1, ev1 = rec0, ev           # phase 1's record and evaluation, when it runs
    # Phase 1 runs when a fresh run starts at a design that violates a constraint. It
    # minimizes the infeasibility from the start without the objective. When the
    # infeasibility stays above `feasibility`, phase 1 returns a verdict; otherwise the main
    # run starts from its first design within the limits.
    if phase1 && state === nothing && maximum(g; init = -Inf) > feasibility
        r1 = _phase1(s, cs, x0, xs, lower, upper; evaluations, show_trace, g_tol,
                     callback, feasibility, deadline, inner_reltol)
        ne, ni = r1.evaluations, r1.iterations
        p1 = (; evaluations = r1.evaluations, iterations = r1.iterations, feasible = r1.feasible,
              x = copy(r1.x), verdict = r1.verdict, reason = r1.why)
        ev1 = evaluate(s, r1.x)
        rec1 = _iterate_record(cs, ev1, r1.x, Float64(ForwardDiff.value(objective(ev1))), 0, feasibility)
        _less_violating(rec1, least, feasibility) && (least = rec1)
        rec1.violation <= feasibility && (feas = rec1)
        if r1.verdict || (r1.stopped && !r1.feasible)
            # The problem is infeasible, because phase 1 converged above `feasibility`, or
            # the budget, the time limit or the callback stopped phase 1 before it reached a
            # design within the limits. The result is the design of least violation, which
            # is the phase-1 design unless the start has less violation.
            xr = least.x
            ev1 = xr === r1.x ? ev1 : evaluate(s, xr)
            rows = _constraint_rows(ev1, cs, nothing, feasibility)
            rel = Relaxation[]
            if r1.verdict && relaxation
                rel = _relaxations(s, cs, rows, lens, collect(Float64, x0), collect(Float64, xr),
                                   lower, upper; iterations,
                                   evaluations = evaluations - ne, show_trace, g_tol,
                                   feasibility, outer, deadline, inner_reltol, inner_evaluations)
                ne += sum(r -> r.evaluations, rel; init = 0)
            end
            return _pointing_result(s, xr, Float64(ForwardDiff.value(objective(ev1))), f0, ni, ne,
                                    r1.stopped, r1.verdict ? :converged : r1.reason, r1.result,
                                    lower, upper; constraints = rows, phase1 = p1, relaxations = rel,
                                    least_violation = least, last_feasible = feas)
        end
    end
    # The main run starts from the design of phase 1 when that design is within the
    # limits, and from the start otherwise. The objective and violations of that design
    # set the first penalty.
    fμ = f0
    if p1 !== nothing && p1.feasible
        x = copy(p1.x)
        fμ = rec1.value
        g = _all_violations(cs, ev1)
    end
    # The inner problems work on these constraints. For an exact ceiling they contain the
    # stand-in at a working limit, which is lowered while the exact value ends above the
    # limit.
    work = cs
    if state === nothing
        λ = zeros(length(g))
        μ = clamp(10 * max(1.0, abs(fμ)) / max(1.0, sum(abs2, max.(g, 0.0)) / 2), 1e-6, 1e6)
        μmax = 1e8 * μ
        tol = max(Float64(g_tol), 1e-3)
        Vold = Inf
    else
        # A warm restart (`from`) takes the multipliers, penalty, inner tolerance and
        # working limits from the end of the previous run.
        length(state.λ) == length(g) && length(state.limits) == length(cs) &&
            all(i -> length(state.limits[i]) == length(cs[i].limit), eachindex(state.limits)) ||
            throw(DimensionMismatch("`from` is a run with other constraints"))
        λ, μ, μmax, tol, Vold = copy(state.λ), state.μ, state.μmax, state.tol, state.Vold
        work = map((c, L) -> Constraint(c.name, c.f, copy(L), c.scale, c.upper, c.labels), cs,
                   Tuple(state.limits))
    end
    last = nothing
    reason = :converged
    evx = ev
    capped = 0
    rec = least
    limits(w) = [copy(c.limit) for c in w]
    while true
        k += 1
        λs = _split_multipliers(λ, lens)
        aug = _Augmented(objective, work, λs, μ)
        cb = callback === nothing ? nothing : let ne = ne, ni = ni, k = k
            st -> callback(merge(st, (; evaluations = st.evaluations + ne,
                                      iteration = st.iteration + ni, outer = k)))
        end
        r = _lbfgsb(s, aug, x0, x, lower, upper; iterations,
                    evaluations = evaluations - ne, show_trace, g_tol = tol, callback = cb,
                    f_reltol = inner_reltol, deadline, cap = inner_evaluations)
        x, ne, ni, last = copy(r.x), ne + r.evaluations, ni + r.iterations, r
        evx = evaluate(s, x)
        g = _all_violations(work, evx)
        V = maximum(abs, max.(g, -λ ./ μ); init = 0.0)
        λinner = λ
        λ = max.(0.0, λ .+ μ .* g)
        rec = _iterate_record(cs, evx, x, Float64(ForwardDiff.value(objective(evx))), k, feasibility)
        _less_violating(rec, least, feasibility) && (least = rec)
        rec.violation <= feasibility && (feas = rec)
        if r.stopped
            reason = r.reason
            # A restart resumes the interrupted inner problem.
            state = (; λ = λinner, μ, μmax, tol, Vold, limits = limits(work), x = copy(x))
            break
        end
        # An inner run ended by its cap (`iterations`, `inner_evaluations`) leaves its
        # subproblem unsolved. Above the limits, the multipliers are updated, and the
        # penalty, the tolerance and the reference of the penalty test are held. A V that
        # has not halved after such a run measures the cap, and a growing penalty would
        # lead to the infeasible verdict. Within the limits the usual test applies, since V
        # is then the complementarity and no verdict can follow.
        capped_run = r.reason === :iterations
        capped += capped_run
        if capped_run && maximum(g; init = -Inf) > feasibility
            state = (; λ, μ, μmax, tol, Vold, limits = limits(work), x = copy(x))
            if k >= outer
                reason = :iterations
                break
            end
            continue
        end
        feasible = maximum(g; init = -Inf) <= feasibility
        # The run has converged when it is feasible and complementary and the inner problem
        # is solved at the final tolerance. Without the tolerance condition, a start that
        # is already feasible and active, as on a restart, would pass on the first, loose
        # inner run, with zero multipliers wherever the gradient of the objective is below
        # 1e-3.
        if feasible && V <= feasibility && r.reason === :converged && tol <= g_tol
            # For an exact ceiling above its limit, the working limit of its stand-in is
            # scaled by limit/exact, and the run continues from this design with the same
            # multipliers.
            lowered = _lower_exact_limits(cs, work, evx, feasibility)
            if lowered === nothing
                state = (; λ, μ, μmax, tol, Vold = V, limits = limits(work), x = copy(x))
                break
            end
            work = lowered
            # The violation of the new limit is not a stall, so the penalty is kept.
            V = Inf
        end
        # Set the penalty and tolerance of the next outer iteration, which a restart also
        # continues with.
        μnext = V > Vold / 2 && V < Inf ? 10μ : μ
        tolnext = max(Float64(g_tol), tol / 10)
        state = (; λ, μ = μnext, μmax, tol = tolnext, Vold = V, limits = limits(work), x = copy(x))
        if k >= outer
            reason = :iterations
            break
        end
        μ = μnext
        if !feasible && μ > μmax
            # The run is infeasible at the largest penalty. After an iterate within the
            # limits, the penalty ceiling ends the run in the same way as the outer budget.
            feas === nothing || (reason = :iterations)
            break
        end
        Vold = V
        tol = tolnext
    end
    # A run that ends with a violation above `feasibility` returns its iterate of least
    # violation, which is the start, the phase-1 design or an outer iterate.
    if rec.violation > feasibility && least !== rec
        x = copy(least.x)
        evx = evaluate(s, x)
    end
    # The reported multipliers are those that best satisfy stationarity at the design. A
    # restart continues from the iteration estimates in `state`.
    kkt = _active_set_multipliers(s, objective, work, x, lower, upper,
                                  _all_violations(work, evx), feasibility)
    rows = _constraint_rows(evx, cs, kkt.λ, feasibility)
    stopped = last.stopped
    return _pointing_result(s, x, Float64(ForwardDiff.value(objective(evx))), f0, ni, ne,
                            stopped, reason, last.result, lower, upper; constraints = rows,
                            outer = k, state, objective, phase1 = p1, gradients = kkt.gradients,
                            stationarity = kkt.stationarity, least_violation = least,
                            last_feasible = feas, capped)
end

# The record of an iterate for `least_violation` and `last_feasible`. It contains the
# parameters, the largest normalized violation over the reported rows, the objective and
# the outer iteration. The violation takes the exact value of an exact ceiling and is zero
# when the constraint is met for every component. Outer iteration 0 is the start or the
# phase-1 design.
function _iterate_record(cs::Tuple, ev, x, value::Float64, outer::Int, feasibility::Float64)
    rows = _constraint_rows(ev, cs, nothing, feasibility)
    v = max(0.0, maximum(r -> r.violation, rows; init = 0.0))
    return (; x = collect(Float64, x), violation = v, value, outer)
end

# Whether iterate `a` has less violation than `b`: violations within `feasibility` count as
# equal, and then the lower objective is less.
function _less_violating(a, b, feasibility::Float64)
    va, vb = max(a.violation, feasibility), max(b.violation, feasibility)
    return va < vb || va == vb && a.value < b.value
end

# Compute the multipliers reported at the end of a constrained run. They are the
# nonnegative least-squares solution of the stationarity condition
# ∂f/∂xᵢ + Σⱼ λⱼ ∂gⱼ/∂xᵢ = 0, set up as follows:
# - The rows are the parameters free to move inside the box. A parameter at a bound is
#   left out of the fit, because its row is the multiplier of the bound, ∂L/∂xᵢ, of
#   either sign. A sign that would move the parameter inward shows that the bound does
#   not hold it.
# - The unknowns are the multipliers of the components within `feasibility` of their
#   working limits or beyond them (`gⱼ ≥ -feasibility`). By complementarity the others
#   are zero.
# - Where several solutions fit equally, the least-norm one is taken.
# The first-order update λ ← max(0, λ + μ g) is exact only after an inner run solved to
# stationarity. An inner run ended by its value tolerance leaves an error of μ times its
# constraint error, which zeroes the multipliers of limits a few 1e-5 inside at a large μ
# (research notes, "Multipliers from the converged active set"). The function also
# returns the gradients of the objective and of the Lagrangian f + Σⱼ λⱼ gⱼ at `x`, for
# the bound sensitivities, and the largest |∂L/∂xᵢ| over the free parameters. That value
# is the stationarity left with these multipliers, on the scale of the `g_tol` of the
# inner runs. The cost is one Jacobian, which costs as much as a gradient.
Base.@nospecializeinfer function _active_set_multipliers(s::_AbstractScenario, @nospecialize(objective),
                                                         @nospecialize(work::Tuple), x::Vector{Float64},
                                                         lower, upper, g::Vector{Float64},
                                                         feasibility::Float64)
    λ = zeros(length(g))
    side = map(_bound_side, x, lower, upper)
    act = findall(>=(-feasibility), g)
    Jm = _scenario_jacobian(s, _KKTValues(objective, work), x)     # f first, then each gᵢ
    gf = Jm[1, :]
    # The stationarity rows are the parameters free to move, inside the box and not fixed.
    free = findall(i -> side[i] == 0 && lower[i] < upper[i], eachindex(x))
    if !isempty(act)
        A = Jm[1 .+ act, free]'
        # Where the multipliers are not unique, as for components with the same gradient,
        # take the least-norm solution. A weight of 1e-6 times the largest column norm on
        # ‖λ‖ selects it, and it moves the multipliers of well-conditioned columns by about
        # 1e-12 relative.
        m = length(act)
        ρ = 1e-6 * maximum(norm, eachcol(A); init = 0.0)
        λ[act] .= _nnls([A; ρ * Matrix(1.0I, m, m)], [-gf[free]; zeros(m)])
    end
    gL = gf .+ view(Jm, 2:size(Jm, 1), :)' * λ
    return (; λ, gradients = (gf, gL), stationarity = maximum(i -> abs(gL[i]), free; init = 0.0))
end

# Stack the objective and the normalized violations of the constraints at `ev`. These are
# the rows of the Jacobian of `_active_set_multipliers`.
struct _KKTValues{O,C<:Tuple}
    objective::O
    constraints::C
end

function (a::_KKTValues)(ev)
    J = a.objective(ev)
    gs = map(c -> _violations(c, ev), a.constraints)
    out = promote_type(typeof(J), map(eltype, gs)...)[J]
    return _append_all!(out, gs)
end

_append_all!(out, ::Tuple{}) = out
_append_all!(out, gs::Tuple) = _append_all!(append!(out, first(gs)), Base.tail(gs))

# Phase 1 of a constrained run minimizes the infeasibility ψ = ½ Σᵢ max(0, gᵢ)² over the
# box, with gᵢ the normalized violation of each component. ψ is the limit of the augmented
# Lagrangian as μ → ∞, and it is zero exactly on the feasible set. Phase 1 leaves out the
# objective and uses the same L-BFGS-B with its inner-run value tolerance. It stops at the
# first design with ψ ≤ ½ feasibility², where every violation is within `feasibility`
# (`feasible`). When it converges above that level with a violation over `feasibility`,
# the design is a local minimum of the infeasibility. That result is the verdict that the
# constraints cannot all be met from this start (`verdict`).
struct _Infeasibility{C<:Tuple}
    constraints::C
end
(a::_Infeasibility)(ev) = sum(map(c -> _hinge2(_violations(c, ev)), a.constraints))

function _hinge2(g::AbstractVector{T}) where {T}
    s = zero(T)
    for x in g
        s += max(zero(T), x)^2
    end
    return s / 2
end

Base.@nospecializeinfer function _phase1(s::_AbstractScenario, @nospecialize(cs::Tuple), x0, xs, lower, upper;
                                         evaluations::Integer, show_trace::Bool,
                                         g_tol, @nospecialize(callback), feasibility::Float64,
                                         deadline::Float64, inner_reltol::Float64)
    thr = feasibility^2 / 2
    found, stalled = Ref(false), Ref(false)
    best = Float64[]                       # least ψ after each evaluation
    W = _PHASE1_WINDOW
    # Phase 1 passes `outer = 0` to the callback.
    stop = st -> begin
        st.value <= thr && (found[] = true; return true)
        push!(best, st.best_value)
        n = length(best)
        n > W && best[n] > (1 - _PHASE1_STALL) * best[n-W] && (stalled[] = true; return true)
        return callback === nothing ? false : callback(merge(st, (; outer = 0))) === true
    end
    # Phase 1 has a separate budget with no cap on its iterations. The `iterations` cap of
    # the inner runs applies to the subproblems of the augmented Lagrangian, and phase 1 is
    # not one of them. Phase 1 therefore ends on convergence, on its stall rule, on the
    # evaluation budget of the run, at the time limit or by the callback.
    r = _lbfgsb(s, _Infeasibility(cs), x0, xs, lower, upper; iterations = typemax(Int), evaluations,
                show_trace, g_tol, callback = stop, f_reltol = inner_reltol, deadline)
    # A converged or stalled phase 1 ends at a local minimum of ψ, or as near one as the
    # stall rule can resolve.
    ended = !r.stopped && r.reason === :converged || r.stopped && r.reason === :callback && stalled[]
    feasible = found[] || r.value <= thr
    g = feasible ? Float64[] : _all_violations(cs, evaluate(s, r.x))
    verdict = !feasible && ended && maximum(g; init = -Inf) > feasibility
    # An ended phase 1 whose violations are all within `feasibility` counts as feasible.
    # This case arises when many components each lie just inside it.
    feasible |= ended && !verdict
    stopped = r.stopped && !found[] && !stalled[]
    # Record why phase 1 ended. `:feasible`, `:converged` and `:stalled` end it with a
    # verdict or with a design within the limits. Otherwise the reason of the stop is
    # recorded (`:evaluations`, `:time_limit` or `:callback`).
    why = found[] || r.value <= thr ? :feasible : stalled[] ? :stalled : stopped ? r.reason : :converged
    return (; r.x, r.value, r.iterations, r.evaluations, stopped, r.reason, r.result, feasible, verdict, why)
end

# The stall rule of phase 1 uses engineering defaults. Phase 1 ends when its least
# infeasibility has fallen by less than `_PHASE1_STALL` (relative) over the last
# `_PHASE1_WINDOW` evaluations. On a drive problem whose limits can be reached, phase 1
# reaches them in tens of evaluations. When they cannot be reached, phase 1 creeps toward
# its local minimum of ψ for thousands (research notes, "A verdict on conflicting
# constraints").
const _PHASE1_WINDOW = 100
const _PHASE1_STALL = 0.01

# The objective of a relaxation is a smooth largest positive normalized violation of the
# components of one constraint, `m (Σₖ max(0, gₖ/m)ᵖ)^(1/p)`. Here `m` is the largest
# violation, held constant, which is valid because the norm is homogeneous. The value
# exceeds the largest by at most a factor n^(1/p) for n equal components, and it is zero
# once the constraint is met for every component.
struct _Relaxed{C}
    constraint::C
    p::Float64
end
(r::_Relaxed)(ev) = _positive_norm(_violations(r.constraint, ev), r.p)

function _positive_norm(g::AbstractVector{T}, p::Float64) where {T}
    m = maximum(x -> Float64(ForwardDiff.value(x)), g; init = 0.0)
    m > 0 || return zero(T)
    s = zero(T)
    for x in g
        x > 0 && (s += (x / m)^p)
    end
    return m * s^(1 / p)
end

# `_RELAX_P` is the p of the smooth largest violation of a relaxation, an engineering
# default. With it, 4 equal components, as for the gate multiples of a study, give at most
# 2.2% above their largest.
const _RELAX_P = 64.0

# Compute the relaxation of each constraint that is violated or active at the infeasible
# design `x1`, whose rows are `rows`. A constraint of the conflict can end at its limit at
# `x1`, with the violation carried by the others. Each relaxation is the smallest largest
# violation of the components of the constraint, with the others as constraints. It is a
# constrained run from `x1`, with phase 1 on the others first, within the remaining budget
# and time.
Base.@nospecializeinfer function _relaxations(s::_AbstractScenario, @nospecialize(cs::Tuple), rows,
                                              @nospecialize(lens::Tuple), x0::Vector{Float64},
                                              x1::Vector{Float64}, lower, upper; iterations::Integer,
                                              evaluations::Integer, show_trace::Bool, g_tol,
                                              feasibility::Float64, outer::Int, deadline::Float64,
                                              inner_reltol::Float64, inner_evaluations::Int)
    offs = cumsum((0, lens...))
    out = Relaxation[]
    used = 0
    for j in eachindex(cs)
        any(i -> rows[i].status !== :satisfied, offs[j]+1:offs[j+1]) || continue
        c = cs[j]
        others = cs[[i for i in eachindex(cs) if i != j]]
        budget = evaluations - used
        if budget <= 0 || time() > deadline
            push!(out, Relaxation(string(c.name), NaN, copy(c.limit), fill(NaN, length(c.limit)),
                                  c.upper, :stopped, copy(x1), ConstraintRow[], 0))
            continue
        end
        obj = _Relaxed(c, _RELAX_P)
        r = if isempty(others)
            # With no other constraint, the run is box-constrained, as `optimize_pointing`
            # runs without constraints.
            q = _lbfgsb(s, obj, x0, x1, lower, upper; iterations, evaluations = budget, show_trace, g_tol,
                        callback = nothing, deadline)
            _pointing_result(s, q.x, q.value, q.initial_value, q.iterations, q.evaluations, q.stopped,
                             q.reason, q.result, lower, upper)
        else
            _constrained_pointing(s, obj, Tuple(others), x0, x1, lower, upper; iterations,
                                  evaluations = budget, show_trace, g_tol, callback = nothing,
                                  feasibility, outer, deadline, phase1 = true, inner_reltol,
                                  inner_evaluations)
        end
        used += r.evaluations
        evr = evaluate(s, r.x)
        allrows = _constraint_rows(evr, cs, nothing, feasibility)
        status = r.outcome in (:interior, :box) || (isempty(others) && r.reason === :converged) ? :relaxed :
                 r.outcome in (:infeasible, :infeasible_box) ? :others_conflict : :stopped
        v = status === :others_conflict ? NaN :
            max(0.0, maximum(i -> allrows[i].violation, offs[j]+1:offs[j+1]))
        sg = c.upper ? 1 : -1
        n = max(length(c.limit), length(c.scale))
        lim = [c.limit[length(c.limit) == 1 ? 1 : i] for i in 1:n]
        sc = [c.scale[length(c.scale) == 1 ? 1 : i] for i in 1:n]
        push!(out, Relaxation(string(c.name), v, lim, lim .+ sg .* v .* sc, c.upper, status,
                              copy(r.x), allrows, r.evaluations))
    end
    return out
end

# Lower the working limits of exact ceilings. For each exact ceiling whose exact value at
# `ev` exceeds its limit by more than `feasibility` (normalized), the working limit of its
# stand-in is multiplied by limit/exact. The function returns the new working constraints,
# or `nothing` when every exact value is within the tolerance. Near the limit the stand-in,
# a p-norm, is nearly proportional to the exact peak, so the scaled limit puts the exact
# value near the limit, and a few rounds remove the remainder.
function _lower_exact_limits(cs::Tuple, work::Tuple, ev, feasibility::Float64)
    lowered = false
    new = map(cs, work) do c, w
        e = _exact_value(c.f, ev)
        (e === nothing || !c.upper) && return w
        L, S = c.limit[1], c.scale[1]
        (e - L) / S > feasibility || return w
        lowered = true
        return Constraint(w.name, w.f, w.limit .* (L / e), w.scale, w.upper, w.labels)
    end
    return lowered ? new : nothing
end

_split_multipliers(λ::Vector{Float64}, lens::Tuple) =
    (offs = cumsum((0, lens...)); ntuple(i -> λ[offs[i]+1:offs[i+1]], length(lens)))
