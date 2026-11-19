# Inequality constraints for `optimize_pointing` (an augmented Lagrangian around its
# L-BFGS-B), and the ready-made ones: a laser-irradiance ceiling on named parts and a
# per-degree gate on the capsule's mode shape; the drive objective they go with.

"""
    Constraint(name, f; upper = nothing, lower = nothing, scale = nothing, labels = nothing)

An inequality constraint for [`optimize_pointing`](@ref): `f(ev) ≤ upper` (a ceiling) or
`f(ev) ≥ lower` (a floor), with `f` any function of an [`Evaluation`](@ref) (its beams
`ev.beams` included), generic in the number type so that the optimizer's
automatic-differentiation gradients pass through it. Give exactly one of `upper` and
`lower`.

`f` may return a number or a vector; a vector constraint has one component per entry,
each with its own multiplier and row in the report, against a `limit` that is a number
(shared) or a vector of the same length. `labels` name the components in reports
(default `"1"`, `"2"`, …): a collection of names, or a function of the evaluation that
returns them, for a constraint whose components follow the evaluation (one per beam,
as in [`clearance_floor`](@ref)). Within one run the number of components is fixed; a
constraint with a shared limit and labels from the evaluation serves scenarios with
different numbers of beams (a stage on quads, then one on their beams).

`scale` (a number or one per component, positive) sets the units of the violation: the
optimizer works on `±(f(ev) - limit)/scale` and its `feasibility` tolerance is in those
units. By default `scale = |limit|`, so a violation is relative (1e-3: 0.1 % of the
limit); a zero limit needs a `scale`.

Ready-made constraints: [`irradiance_ceiling`](@ref), [`mode_gate`](@ref),
[`clearance_floor`](@ref) (a keep-out's clearance per beam), [`light_ceiling`](@ref)
and [`light_floor`](@ref) (the laser light on named parts). Any other quantity of the
beams alone is `Constraint(:name, ev -> quantity(ev.beams); lower = minimum)`.
"""
struct Constraint{F,L}
    name::Symbol
    f::F
    limit::Vector{Float64}
    scale::Vector{Float64}
    upper::Bool
    labels::L              # Vector{String}, or a function of the evaluation
end

function Constraint(name, f; upper = nothing, lower = nothing, scale = nothing, labels = nothing)
    (upper === nothing) == (lower === nothing) &&
        throw(ArgumentError("give exactly one of `upper` and `lower`"))
    lim = upper === nothing ? lower : upper
    L = collect(Float64, lim isa Real ? (lim,) : lim)
    S = scale === nothing ? abs.(L) : collect(Float64, scale isa Real ? (scale,) : scale)
    !isempty(L) && all(isfinite, L) || throw(ArgumentError("the limit must be finite"))
    length(S) == 1 || length(S) == length(L) || length(L) == 1 ||
        throw(DimensionMismatch("give one scale or one per component of the limit"))
    all(x -> isfinite(x) && x > 0, S) ||
        throw(ArgumentError("the scale must be positive and finite (a zero limit needs a `scale`)"))
    lab = labels === nothing ? String[] : labels isa Function ? labels : collect(String, string.(labels))
    return Constraint(Symbol(name), f, L, S, upper !== nothing, lab)
end

# The labels of a constraint's `n` components at `ev` (empty: numbered).
_labels(c::Constraint{<:Any,Vector{String}}, ev, n) = c.labels
function _labels(c::Constraint, ev, n)
    lab = collect(String, string.(c.labels(ev)))
    length(lab) == n ||
        throw(DimensionMismatch("constraint $(c.name): $(length(lab)) labels for $n components"))
    return lab
end

Base.show(io::IO, c::Constraint) =
    print(io, "Constraint(:", c.name, c.upper ? " ≤ " : " ≥ ",
          length(c.limit) == 1 ? c.limit[1] : c.limit, ")")

# The constraint's components normalized as violations: `±(value - limit)/scale`
# (positive when violated), in the evaluation's number type.
function _violations(c::Constraint, ev)
    v = c.f(ev)
    s = c.upper ? 1 : -1
    return v isa Number ? [s * (v - c.limit[1]) / c.scale[1]] : s .* (v .- c.limit) ./ c.scale
end

"""
    check_constraints(ev, constraints; feasibility = 1e-3) -> Vector{ConstraintRow}

The value of each component of `constraints` (a [`Constraint`](@ref), or a tuple or
vector of them) at the evaluation `ev`, against its limit, as [`ConstraintRow`](@ref)s:
the report of `PointingResult.constraints`, for any design (a reference, a staged
run's result). `feasibility` is the tolerance on the normalized violation that
separates `:violated`, `:active` and `:satisfied`.
"""
check_constraints(ev::Evaluation, c::Constraint; kwargs...) = check_constraints(ev, (c,); kwargs...)
check_constraints(ev::Evaluation, cs; feasibility::Real = 1e-3) =
    _constraint_rows(ev, Tuple(cs), nothing, Float64(feasibility))

function _constraint_rows(ev, cs::Tuple, multipliers, feasibility::Float64)
    rows = ConstraintRow[]
    k = 0
    for c in cs
        # the reported value: the exact quantity where the constraint's function is a
        # smooth stand-in for one (`_exact_value`), else the function's own
        v = something(_exact_value(c.f, ev), c.f(ev))
        vals = v isa Number ? [Float64(ForwardDiff.value(v))] : Float64.(ForwardDiff.value.(v))
        sg = c.upper ? 1 : -1
        n = length(vals)
        labels = _labels(c, ev, n)
        for i in 1:n
            k += 1
            lim = c.limit[length(c.limit) == 1 ? 1 : i]
            gi = sg * (vals[i] - lim) / c.scale[length(c.scale) == 1 ? 1 : i]
            label = n == 1 && isempty(labels) ? string(c.name) :
                    string(c.name, "[", isempty(labels) ? string(i) : labels[i], "]")
            status = gi > feasibility ? :violated : gi >= -feasibility ? :active : :satisfied
            push!(rows, ConstraintRow(label, vals[i], lim, c.upper, gi,
                                      multipliers === nothing ? 0.0 : multipliers[k], status))
        end
    end
    return rows
end

# The exact quantity that a constraint's function stands in for smoothly (a scalar, at a
# Float64 evaluation), or `nothing` where the function is itself the quantity: the
# irradiance ceiling's peak when it holds the peak (`exact = true`, its default).
_exact_value(f, ev) = nothing

# All constraints' normalized violations at `ev`, concatenated (Float64).
_all_violations(cs::Tuple, ev) =
    reduce(vcat, map(c -> Float64.(ForwardDiff.value.(_violations(c, ev))), cs); init = Float64[])

# The augmented Lagrangian (Powell–Hestenes–Rockafellar, for inequalities) of an
# objective with multipliers `λ` (one vector per constraint) and penalty `μ`:
# f + μ/2 Σᵢ [max(0, gᵢ + λᵢ/μ)² - (λᵢ/μ)²], with gᵢ the normalized violations. A function
# of the evaluation only (its fields are fixed during an inner run): thread-safe.
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

# The constrained run of `optimize_pointing`: an augmented Lagrangian whose inner
# problems (the box-constrained minimizations of `_Augmented`) are solved by `_lbfgsb`,
# each warm-started from the last. Multipliers and penalty are updated as in ALGENCAN
# (Andreani, Birgin, Martínez and Schuverdt, SIAM J. Optim. 18, 1286, 2008; Birgin and
# Martínez, "Practical Augmented Lagrangian Methods for Constrained Optimization",
# SIAM 2014): λ ← max(0, λ + μ g); μ grows tenfold when the infeasibility–
# complementarity measure V = maxᵢ |max(gᵢ, -λᵢ/μ)| has not halved; the first μ balances
# the objective against the initial violations; the inner tolerance on the projected
# gradient tightens tenfold per outer iteration down to `g_tol`, and an inner run also
# ends when an iteration changes its value by less than `_AL_F_RELTOL` relative (the
# drive objectives are flat along the constraints: without it an inner run creeps for
# hundreds of iterations at a 1e-6 gain).
const _AL_F_RELTOL = Ref(1e-10)

Base.@nospecializeinfer function _constrained_pointing(
        s::_AbstractScenario, @nospecialize(objective), @nospecialize(cs::Tuple), x0, xs, lower, upper;
        iterations::Integer, evaluations::Integer, show_trace::Bool, g_tol,
        @nospecialize(callback), feasibility::Float64, outer::Int, @nospecialize(state = nothing),
        deadline::Float64 = Inf, phase1::Bool = false, relaxation::Bool = false)
    for c in cs
        c isa Constraint || throw(ArgumentError("constraints must be `Constraint`s, got $(typeof(c))"))
    end
    feasibility > 0 || throw(ArgumentError("feasibility must be positive"))
    phase1 |= relaxation
    ev = evaluate(s, xs)
    f0 = Float64(ForwardDiff.value(objective(ev)))
    g = _all_violations(cs, ev)
    lens = map(c -> length(_violations(c, ev)), cs)
    x, ne, ni, k = copy(xs), 0, 0, 0
    p1 = nothing
    # phase 1 (a fresh run whose start violates a constraint): the smallest infeasibility
    # from the start, ignoring the objective; a verdict when it stays above `feasibility`,
    # else the main run starts from its first design within the limits
    if phase1 && state === nothing && maximum(g; init = -Inf) > feasibility
        r1 = _phase1(s, cs, x0, xs, lower, upper; iterations, evaluations, show_trace, g_tol,
                     callback, feasibility, deadline)
        ne, ni = r1.evaluations, r1.iterations
        p1 = (; evaluations = r1.evaluations, iterations = r1.iterations, feasible = r1.feasible,
              x = copy(r1.x), verdict = r1.verdict)
        if r1.verdict || (r1.stopped && !r1.feasible)
            # infeasible (phase 1 converged above `feasibility`), or stopped by the budget,
            # the time limit or the callback before a design within the limits
            ev1 = evaluate(s, r1.x)
            rows = _constraint_rows(ev1, cs, nothing, feasibility)
            rel = Relaxation[]
            if r1.verdict && relaxation
                rel = _relaxations(s, cs, rows, lens, collect(Float64, x0), collect(Float64, r1.x),
                                   lower, upper; iterations,
                                   evaluations = evaluations - ne, show_trace, g_tol,
                                   feasibility, outer, deadline)
                ne += sum(r -> r.evaluations, rel; init = 0)
            end
            return _pointing_result(s, r1.x, Float64(ForwardDiff.value(objective(ev1))), f0, ni, ne,
                                    r1.stopped, r1.verdict ? :converged : r1.reason, r1.result,
                                    lower, upper; constraints = rows, phase1 = p1, relaxations = rel)
        end
    end
    # the main run's start: phase 1's design within the limits, if any (its objective and
    # violations set the first penalty), else the start
    fμ = f0
    if p1 !== nothing && p1.feasible
        x = copy(p1.x)
        evs = evaluate(s, x)
        fμ = Float64(ForwardDiff.value(objective(evs)))
        g = _all_violations(cs, evs)
    end
    # the constraints as the inner problems see them: an exact ceiling's stand-in at a
    # working limit, lowered while the exact value ends above the limit
    work = cs
    if state === nothing
        λ = zeros(length(g))
        μ = clamp(10 * max(1.0, abs(fμ)) / max(1.0, sum(abs2, max.(g, 0.0)) / 2), 1e-6, 1e6)
        μmax = 1e8 * μ
        tol = max(Float64(g_tol), 1e-3)
        Vold = Inf
    else
        # a warm restart (`from`): the multipliers, penalty, inner tolerance and working
        # limits where the previous run left them
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
                    f_reltol = _AL_F_RELTOL[], deadline)
        x, ne, ni, last = copy(r.x), ne + r.evaluations, ni + r.iterations, r
        evx = evaluate(s, x)
        g = _all_violations(work, evx)
        V = maximum(abs, max.(g, -λ ./ μ); init = 0.0)
        λinner = λ
        λ = max.(0.0, λ .+ μ .* g)
        if r.stopped
            reason = r.reason
            # a restart resumes the interrupted inner problem
            state = (; λ = λinner, μ, μmax, tol, Vold, limits = limits(work))
            break
        end
        feasible = maximum(g; init = -Inf) <= feasibility
        # converged: feasible and complementary, with the inner problem solved at the final
        # tolerance (a start already feasible and active, as on a restart, would otherwise
        # pass on the first, loose inner run with zero multipliers wherever the objective's
        # gradient is below 1e-3)
        if feasible && V <= feasibility && r.reason === :converged && tol <= g_tol
            # an exact ceiling above its limit: its stand-in's working limit scaled by
            # limit/exact, and the run continued from here with the same multipliers
            lowered = _lower_exact_limits(cs, work, evx, feasibility)
            if lowered === nothing
                state = (; λ, μ, μmax, tol, Vold = V, limits = limits(work))
                break
            end
            work = lowered
            # the new limit's violation is not a stall: the penalty stays
            V = Inf
        end
        # the next outer iteration's penalty and tolerance (what a restart continues with)
        μnext = V > Vold / 2 && V < Inf ? 10μ : μ
        tolnext = max(Float64(g_tol), tol / 10)
        state = (; λ, μ = μnext, μmax, tol = tolnext, Vold = V, limits = limits(work))
        if k >= outer
            reason = :iterations
            break
        end
        μ = μnext
        !feasible && μ > μmax && break                    # infeasible at the largest penalty
        Vold = V
        tol = tolnext
    end
    # the reported multipliers: those that best satisfy stationarity at the design (the
    # iteration's own, in `state`, are what a restart continues from)
    kkt = _active_set_multipliers(s, objective, work, x, lower, upper,
                                  _all_violations(work, evx), feasibility)
    rows = _constraint_rows(evx, cs, kkt.λ, feasibility)
    stopped = last.stopped
    return _pointing_result(s, x, Float64(ForwardDiff.value(objective(evx))), f0, ni, ne,
                            stopped, reason, last.result, lower, upper; constraints = rows,
                            outer = k, state, objective, phase1 = p1, gradients = kkt.gradients,
                            stationarity = kkt.stationarity)
end

# The multipliers reported at the end of a constrained run: the nonnegative least-squares
# solution of the stationarity condition ∂f/∂xᵢ + Σⱼ λⱼ ∂gⱼ/∂xᵢ = 0 on the parameters free
# to move (inside the box), over the components within `feasibility` of their (working)
# limits or beyond them (`gⱼ ≥ -feasibility`; zero for the others: complementary), the
# least-norm one where several fit equally. A parameter at a bound is left out of the
# fit: its row is the bound's multiplier, ∂L/∂xᵢ, of either sign (a sign that would move
# the parameter inward shows that the bound does not hold it). The first-order update
# λ ← max(0, λ + μ g) is exact only after an inner run solved to stationarity; an inner
# run ended by its value tolerance leaves an error of μ times its constraint error, which
# zeroes the multipliers of limits a few 1e-5 inside at a large μ (research notes,
# "Multipliers from the converged active set"). Also the gradients of the objective and
# of the Lagrangian f + Σⱼ λⱼ gⱼ at `x`, for the bound sensitivities, and the largest
# |∂L/∂xᵢ| over the free parameters (the stationarity left with these multipliers, on
# the scale of the inner runs' `g_tol`). One Jacobian (a gradient's cost).
Base.@nospecializeinfer function _active_set_multipliers(s::_AbstractScenario, @nospecialize(objective),
                                                         @nospecialize(work::Tuple), x::Vector{Float64},
                                                         lower, upper, g::Vector{Float64},
                                                         feasibility::Float64)
    λ = zeros(length(g))
    side = map(_bound_side, x, lower, upper)
    act = findall(>=(-feasibility), g)
    Jm = _scenario_jacobian(s, _KKTValues(objective, work), x)     # rows: f, then each gᵢ
    gf = Jm[1, :]
    # the stationarity rows: the parameters free to move (inside the box, not fixed)
    free = findall(i -> side[i] == 0 && lower[i] < upper[i], eachindex(x))
    if !isempty(act)
        A = Jm[1 .+ act, free]'
        # where the multipliers are not unique (components with the same gradient) the
        # least-norm solution: a weight 1e-6 × the largest column norm on ‖λ‖, which moves
        # the multipliers of well-conditioned columns by about 1e-12 relative
        m = length(act)
        ρ = 1e-6 * maximum(norm, eachcol(A); init = 0.0)
        λ[act] .= _nnls([A; ρ * Matrix(1.0I, m, m)], [-gf[free]; zeros(m)])
    end
    gL = gf .+ view(Jm, 2:size(Jm, 1), :)' * λ
    return (; λ, gradients = (gf, gL), stationarity = maximum(i -> abs(gL[i]), free; init = 0.0))
end

# The objective and the constraints' normalized violations at `ev`, stacked (the rows of
# the Jacobian of `_active_set_multipliers`).
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

# Nonnegative least squares, argmin ‖A z - b‖ over z ≥ 0, by the active-set method of
# Lawson and Hanson ("Solving Least Squares Problems", Prentice-Hall 1974, ch. 23).
# Rank-deficient passive sets are solved by pivoted QR (a basic solution).
function _nnls(A::Matrix{Float64}, b::Vector{Float64})
    m, n = size(A)
    z = zeros(n)
    n == 0 && return z
    P = falses(n)
    tol = 10 * eps() * max(m, n) * max(opnorm(A, 1), 1e-300) * max(norm(b, Inf), 1e-300)
    w = A' * b
    # columns whose entry left the solution unchanged (its least-squares entry came out
    # nonpositive at once, a rounding-level gain): skipped until the solution moves
    blocked = falses(n)
    for _ in 1:3n
        j = 0
        for i in 1:n
            !P[i] && !blocked[i] && w[i] > tol && (j == 0 || w[i] > w[j]) && (j = i)
        end
        j == 0 && break
        P[j] = true
        entered = true
        while true
            idx = findall(P)
            t = zeros(n)
            t[idx] .= A[:, idx] \ b
            if all(>(0), view(t, idx))
                z = t
                break
            end
            if entered && t[j] <= 0
                P[j] = false
                blocked[j] = true
                break
            end
            entered = false
            α = minimum(z[i] / (z[i] - t[i]) for i in idx if t[i] <= 0)
            z .+= α .* (t .- z)
            for i in idx
                z[i] <= tol && (P[i] = false; z[i] = 0.0)
            end
        end
        blocked[j] || (blocked .= false)
        w = A' * (b - A * z)
    end
    return z
end

# Phase 1 of a constrained run: the box-constrained minimization of the infeasibility
# ψ = ½ Σᵢ max(0, gᵢ)² over every component's normalized violation (the augmented
# Lagrangian's limit as μ → ∞, zero exactly on the feasible set), ignoring the
# objective, by the same L-BFGS-B and its inner-run value tolerance. It stops at the
# first design with ψ ≤ ½ feasibility², where every violation is within `feasibility`
# (`feasible`). Converged above it with a violation over `feasibility`, the design is a
# local minimum of the infeasibility: the verdict that the constraints cannot all be met
# from this start (`verdict`).
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
                                         iterations::Integer, evaluations::Integer, show_trace::Bool,
                                         g_tol, @nospecialize(callback), feasibility::Float64,
                                         deadline::Float64)
    thr = feasibility^2 / 2
    found, stalled = Ref(false), Ref(false)
    best = Float64[]                       # the least ψ after each evaluation
    W = _PHASE1_WINDOW
    # the callback's view: phase 1 reports `outer = 0`
    stop = st -> begin
        st.value <= thr && (found[] = true; return true)
        push!(best, st.best_value)
        n = length(best)
        n > W && best[n] > (1 - _PHASE1_STALL) * best[n-W] && (stalled[] = true; return true)
        return callback === nothing ? false : callback(merge(st, (; outer = 0))) === true
    end
    r = _lbfgsb(s, _Infeasibility(cs), x0, xs, lower, upper; iterations, evaluations, show_trace,
                g_tol, callback = stop, f_reltol = _AL_F_RELTOL[], deadline)
    # converged, or stalled: a local minimum of ψ, or as near one as the stall rule tells
    ended = !r.stopped && r.reason === :converged || r.stopped && r.reason === :callback && stalled[]
    feasible = found[] || r.value <= thr
    g = feasible ? Float64[] : _all_violations(cs, evaluate(s, r.x))
    verdict = !feasible && ended && maximum(g; init = -Inf) > feasibility
    # an ended phase 1 whose violations are all within `feasibility` (many components each
    # just inside it) counts as feasible
    feasible |= ended && !verdict
    return (; r.x, r.value, r.iterations, r.evaluations,
            stopped = r.stopped && !found[] && !stalled[], r.reason, r.result, feasible, verdict)
end

# Phase 1's stall rule (engineering defaults): it ends when its least infeasibility has
# fallen by less than `_PHASE1_STALL` (relative) over the last `_PHASE1_WINDOW`
# evaluations. A drive problem's phase 1 that can reach the limits does so in tens of
# evaluations; one that cannot creeps toward its local minimum of ψ for thousands
# (research notes, "A verdict on conflicting constraints").
const _PHASE1_WINDOW = 100
const _PHASE1_STALL = 0.01

# The objective of a relaxation: a smooth largest positive normalized violation of one
# constraint's components, `m (Σₖ max(0, gₖ/m)ᵖ)^(1/p)` with `m` the largest (a constant:
# the norm is homogeneous), above the largest by at most a factor n^(1/p) for n equal
# components; zero once every component holds.
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

# The p of a relaxation's smooth largest violation (engineering default: as the gate
# multiples of a study, 4 components equal at most 2.2 % above their largest)
const _RELAX_P = 64.0

# The relaxation of each constraint violated or active at the infeasible design `x1`
# (rows `rows`; a constraint of the conflict can end at its limit there, the violation
# carried by the others): the smallest largest violation of its components with the
# others as constraints, each a constrained run (phase 1 on the others first) from `x1`,
# within the remaining budget and time.
Base.@nospecializeinfer function _relaxations(s::_AbstractScenario, @nospecialize(cs::Tuple), rows,
                                              @nospecialize(lens::Tuple), x0::Vector{Float64},
                                              x1::Vector{Float64}, lower, upper; iterations::Integer,
                                              evaluations::Integer, show_trace::Bool, g_tol,
                                              feasibility::Float64, outer::Int, deadline::Float64)
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
            # the constraint alone: a box-constrained run, as `optimize_pointing` without constraints
            q = _lbfgsb(s, obj, x0, x1, lower, upper; iterations, evaluations = budget, show_trace, g_tol,
                        callback = nothing, deadline)
            _pointing_result(s, q.x, q.value, q.initial_value, q.iterations, q.evaluations, q.stopped,
                             q.reason, q.result, lower, upper)
        else
            _constrained_pointing(s, obj, Tuple(others), x0, x1, lower, upper; iterations,
                                  evaluations = budget, show_trace, g_tol, callback = nothing,
                                  feasibility, outer, deadline, phase1 = true)
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

# The working constraints with each exact ceiling whose exact value at `ev` exceeds its
# limit by more than `feasibility` (normalized) given its stand-in's working limit times
# limit/exact, or `nothing` when every exact value is within it. The stand-in tracks the
# exact value nearly in proportion near the limit (a p-norm against the peak), so the
# scaled limit lands the exact value near the limit; a few rounds settle the remainder.
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

# ---------------------------------------------------------------------------------------
# Laser irradiance on named parts
# ---------------------------------------------------------------------------------------

# The elements of the named parts or components (an error for a name not in the target).
_named_elements(t::Target, names) = (_check_names(t, names); elements(t, names...))

function _check_names(t::Target, names)
    isempty(names) && throw(ArgumentError("name the parts (e.g. :wall)"))
    for n in names
        any(c -> c.name === n || c.part === n, t.components) ||
            throw(ArgumentError("the target has no part or component $(repr(n)) " *
                                "(parts: $(join(repr.(parts(t)), ", ")))"))
    end
    return nothing
end

"""
    peak_irradiance(ev, names...; weight = nothing) -> Real

The largest absorbed laser irradiance (TW/mm²; `to_Wcm2` converts) over the elements of
the named parts or components: element power over element area, as in
[`laser_flux`](@ref). Exact over the elements, so it has a kink wherever the peak moves
to another element: the reported quantity of [`irradiance_ceiling`](@ref), whose
optimizer works on the smooth [`irradiance_norm`](@ref). The element average resolves
a peak only as finely as the mesh: a coarser mesh averages the hottest spot over more
area and reports a lower peak (on the dropped-bundle hohlraum, 8 % lower at 0.45 mm
than at 0.2 mm).

`weight`, a function of a beam, counts each beam's deposition times `weight(beam)`
instead of all beams at their own powers: `0` leaves a beam out, so
`weight = b -> b.group in (:outer44, :outer50) ? 1.0 : 0.0` gives the outer cones'
irradiance, and a picket energy over the beam's power gives a fluence map (the
deposition is linear in each beam's power). The weighted deposition is one more
deposit per evaluation, of the evaluation's effective beams through the scenario's
laser–plasma model. The function sees the evaluation's setup beams (their names,
groups, quads and ports).
"""
function peak_irradiance(ev::Evaluation, names::Symbol...; weight = nothing)
    t = ev.scenario.target
    return _peak(_counted_power(ev, weight), t, _named_elements(t, names))
end

_peak(P, t::Target, es) = maximum(P[e] / t.elem_area[e] for e in es)

# The element powers whose irradiance counts: the evaluation's deposition, or that of its
# effective beams with each beam's power times `weight(setup beam)`.
_counted_power(ev::Evaluation, ::Nothing) = ev.deposition.element
function _counted_power(ev::Evaluation, weight)
    s = ev.scenario
    eff = ev.effective
    beams = [remake(b; power = weight(ev.beams[i]) * b.power) for (i, b) in enumerate(eff)]
    return lpi_deposit(s.lpi, s.deposition, beams, s.time).element
end

"""
    irradiance_norm(ev, names...; p = 64, weight = nothing) -> Real

The area-weighted `p`-norm of the absorbed laser irradiance `I` (TW/mm²) over the
elements of the named parts or components, `(Σ Aₑ Iₑᵖ / Σ Aₑ)^(1/p)`: a smooth stand-in for
the [`peak_irradiance`](@ref), below it and approaching it as `p` grows (the power-mean
inequality), with a gradient from every element near the peak instead of from the
hottest one alone. As an area mean it is a property of the irradiance field rather
than of the elements, so it depends less on the mesh than the peak does. `weight` is
as in [`peak_irradiance`](@ref). Generic in the number type.
"""
function irradiance_norm(ev::Evaluation, names::Symbol...; p::Real = 64, weight = nothing)
    t = ev.scenario.target
    return _irradiance_norm(_counted_power(ev, weight), t, _named_elements(t, names), Float64(p))
end

function _irradiance_norm(P::AbstractVector{T}, t::Target, es, p::Float64) where {T}
    p >= 1 || throw(ArgumentError("p must be at least 1"))
    # scaled by the largest value (a constant: the norm is homogeneous) against overflow
    m = maximum(e -> ForwardDiff.value(P[e]) / t.elem_area[e], es; init = 0.0)
    m > 0 || return zero(T)
    s = zero(T)
    A = 0.0
    for e in es
        a = t.elem_area[e]
        s += a * (P[e] / (a * m))^p
        A += a
    end
    return m * (s / A)^(1 / p)
end

"""
    irradiance_ceiling(reference, names...; limit = peak_irradiance(reference, names...),
                       p = 64, weight = nothing, scale = limit, exact = true) -> Constraint

A ceiling on the absorbed laser irradiance over the elements of the named parts or
components: a [`Constraint`](@ref) named `:irradiance` that holds the
[`peak_irradiance`](@ref) at or below `limit` (by default; `exact`, below). Its
function is the smooth stand-in the optimizer works on, the [`irradiance_norm`](@ref)
scaled to equal the peak at the `reference` evaluation, `‖I‖ₚ × peak_ref/‖I‖ₚ,ref`.
With the default `limit` (the reference's peak) the ceiling holds the irradiance at the
reference design's, the
reference being, for example, the pointing before a bundle was dropped; a number in
TW/mm² sets another, measured on the same scale. The reference evaluation must be on
the same target as the optimized scenario (a `subset` of it shares the target), with
the same laser–plasma model (`lpi`, compared field by field, so that a model built
alike passes) at the same time; otherwise evaluating the ceiling throws an
`ArgumentError`. A reference evaluated without the scenario's CBET model gives the peak
and the stand-in's scale of another quantity, and the ceiling would hold the design to
a limit it does not share, so the mismatch is an error rather than a warning:
evaluate the reference in a scenario with the same model, e.g.
`evaluate(remake(reference_scenario; lpi = scenario.lpi), beams)`.
`weight` selects and weighs the beams that count, as in [`peak_irradiance`](@ref) (the
outer cones only, or picket energies for a fluence ceiling), for the reference and the
design alike.

The ceiling is on an absolute irradiance: where the power can fall (power scales,
[`PowerScales`](@ref), or powers among the parameters) the optimizer can meet it by
cutting power, at a loss of drive (on a re-pointing study with a power scale per cone
at most 1, a large one). A [`drive_floor`](@ref) beside it holds the drive.

A ceiling on the wall's peak irradiance at the reference's value is an engineering
limit, not a published criterion: the published scalings of the gold bubble's growth
use the outer cones' picket fluence on the wall (research notes), which `weight`
expresses on the same machinery.

With `exact = true` (the default) the ceiling holds the peak itself. The optimizer
works on the norm, but once the constrained run has converged it compares the
[`peak_irradiance`](@ref) with the limit; while the peak exceeds it by more than the
run's `feasibility` (in units of `scale`), the norm's working limit is scaled by
limit/peak and the run continues from the same design, multipliers and penalty. The
reported value ([`check_constraints`](@ref), `PointingResult.constraints`) is the peak
against the limit, and a converged run ends with the peak at or below the limit within
`feasibility`. The cost is the further outer iterations, when the peak ends above the
limit (research notes, "Constrained re-pointing").

`exact = false` bounds the norm alone, and the norm is the reported value: at the limit
the field's hot spots are as hot, in the `p`-norm mean, as the reference's, while the
peak over the elements may exceed the limit, by a margin that depends on where the run
ends: on the tutorial's mesh 0.1 % at `p = 64` and 3 % at `p = 16`, and on the
dropped-bundle test case (0.45 mm wall elements) 3 % at `p = 64` for one of two designs
that meet the same gates, under another floating-point rounding of the same run
(research notes, "Constrained re-pointing").

Across meshes the ceiling is a ratio to the reference on the same
mesh, which changes less with the element size than the peak's ratio does (rms
difference 1.0 % against 1.3 % between 0.2 and 0.45 mm elements near the limit), while
the peak itself falls on a coarser mesh.
"""
function irradiance_ceiling(reference::Evaluation, names::Symbol...;
                            limit::Union{Nothing,Real} = nothing, p::Real = 64,
                            weight = nothing, scale::Union{Nothing,Real} = nothing,
                            exact::Bool = true)
    t = reference.scenario.target
    es = _named_elements(t, names)
    P = Float64.(ForwardDiff.value.(_counted_power(reference, weight)))
    pk = _peak(P, t, es)
    pk > 0 || throw(ArgumentError("the reference deposits no laser power on $(join(names, ", "))"))
    n_ref = _irradiance_norm(P, t, es, Float64(p))
    lim = limit === nothing ? pk : Float64(limit)
    lim > 0 || throw(ArgumentError("the irradiance limit must be positive"))
    rs = reference.scenario
    value = _CeilingValue(t, es, Float64(p), pk / n_ref, weight, exact, rs.lpi, rs.time,
                          Ref{Any}(rs.lpi))
    return Constraint(:irradiance, value; upper = lim, scale = something(scale, lim))
end

# The ceiling's value: the irradiance norm on the reference's peak scale; with `exact`,
# the reported value is the peak itself. `lpi` and `time` are the reference scenario's
# laser–plasma model and time, which every evaluation must share; `checked` holds the
# last model found equal to it (an identity check from then on).
struct _CeilingValue{W}
    target::Target
    elements::Vector{Int}
    p::Float64
    factor::Float64
    weight::W
    exact::Bool
    lpi::Any
    time::Float64
    checked::Base.RefValue{Any}
end

function (c::_CeilingValue)(ev::Evaluation)
    _check_ceiling(c, ev)
    return c.factor * _irradiance_norm(_counted_power(ev, c.weight), c.target, c.elements, c.p)
end

function _check_ceiling(c::_CeilingValue, ev::Evaluation)
    s = ev.scenario
    s.target === c.target ||
        throw(ArgumentError("irradiance_ceiling: the evaluation is on another target than the reference"))
    _check_reference_model("irradiance_ceiling", c.lpi, c.time, c.checked, s)
    return nothing
end

# A constraint's reference evaluation and the optimized scenario share the laser–plasma
# model (compared field by field; `checked` holds the last model found equal, an identity
# check from then on) and, with a model, its time; otherwise an `ArgumentError` from `who`.
function _check_reference_model(who::String, ref_lpi, ref_time::Float64, checked::Base.RefValue{Any},
                                s::Scenario)
    lpi = s.lpi
    if lpi !== checked[]
        _same_model(lpi, ref_lpi) ||
            throw(ArgumentError("$who: the reference was evaluated with another laser–plasma " *
                                "model than the scenario's ($(_model_name(ref_lpi)) against " *
                                "$(_model_name(lpi))), so its limit bounds another quantity: evaluate the " *
                                "reference in a scenario with the same `lpi`, e.g. " *
                                "`evaluate(remake(reference_scenario; lpi = scenario.lpi), beams)`"))
        checked[] = lpi
    end
    lpi === nothing || s.time == ref_time ||
        throw(ArgumentError("$who: the reference was evaluated at time $(ref_time) ns and the " *
                            "scenario's laser–plasma model at $(s.time) ns: evaluate the reference at the " *
                            "scenario's time (`remake(reference_scenario; time)`)"))
    return nothing
end

# The laser–plasma model and time of a reference evaluation, for a constraint that
# refuses evaluations with another (`_check_reference_model`).
struct _ReferenceModel
    lpi::Any
    time::Float64
    checked::Base.RefValue{Any}
end
_ReferenceModel(s::Scenario) = _ReferenceModel(s.lpi, s.time, Ref{Any}(s.lpi))

_check_model(::Nothing, who, ev::Evaluation) = nothing
_check_model(m::_ReferenceModel, who, ev::Evaluation) =
    _check_reference_model(who, m.lpi, m.time, m.checked, ev.scenario)

_model_name(::Nothing) = "none"
_model_name(m) = string(nameof(typeof(m)))
_model_name(m::Union{Tuple,AbstractVector}) = "(" * join(map(_model_name, m), ", ") * ")"

# Structural equality of two laser–plasma models (two built alike are the same model):
# identical, or of one type with equal fields, arrays and dictionaries by their entries.
function _same_model(a, b)
    a === b && return true
    typeof(a) === typeof(b) || return false
    a isa Union{AbstractArray,Tuple} &&
        return length(a) == length(b) && all(i -> _same_model(a[i], b[i]), eachindex(a, b))
    a isa Union{AbstractDict,AbstractString,Number,Symbol} && return isequal(a, b)
    n = nfields(a)
    n == 0 && return isequal(a, b)
    return all(i -> _same_model(getfield(a, i), getfield(b, i)), 1:n)
end

function _exact_value(c::_CeilingValue, ev::Evaluation)
    c.exact || return nothing
    _check_ceiling(c, ev)
    return _peak(Float64.(ForwardDiff.value.(_counted_power(ev, c.weight))), c.target, c.elements)
end

# ---------------------------------------------------------------------------------------
# Keep-outs and the light on named parts
# ---------------------------------------------------------------------------------------

"""
    clearance_floor(keepout; minimum = keepout.minimum, scale = 1.0, name = :clearance) -> Constraint

A floor on each beam's worst-case clearance (mm) from a keep-out: the
[`Constraint`](@ref) `keepout_clearance(keepout, ev) .≥ minimum`, one component per beam
of the evaluation, labeled by the beam's name. The keep-out is any with a
[`keepout_clearance`](@ref) method and a `minimum` field: an [`ApertureKeepOut`](@ref),
a [`SphereKeepOut`](@ref), a [`HardwareKeepOut`](@ref). A beam the keep-out does not
apply to (clearance `Inf`) gives a component that always holds; a clearance of `-Inf` or
`NaN` throws an `ArgumentError` naming the keep-out and the beam, as in
[`keepout_penalty`](@ref). The keep-out's `weight`
is not used: the constraint replaces its penalty ([`keepout_penalty`](@ref)), so pass
the keep-out to the objective's `keepout` keyword or as a constraint, not both.

`scale` (mm) is the unit of the violation: with the default 1 mm and `feasibility =
1e-3`, a clearance within 1 µm of `minimum` is `:active` (an engineering default, finer
than the keep-out's own accuracy against [`clearance`](@ref), a few micrometers). Give
each of several keep-outs its own `name`.

The floor is per beam, not on the smallest clearance over the beams: in the augmented
Lagrangian each beam's component has its own multiplier, which stays zero for a beam
clear of its minimum (its term then vanishes from the inner problem), and each term is
as smooth as that beam's clearance, while the smallest clearance over the beams has a
kink wherever the closest beam changes (research notes, "Constrained re-pointing").
The components follow the evaluation's beams, so the same constraint serves a stage on
quads and a stage on their beams.

On a [`PhasedScenario`](@ref) the constraint takes the [`PerPhase`](@ref) evaluation
directly, without [`in_phase`](@ref): each beam's component is its least clearance over
the phases, and its rows carry no phase in their names (`leh[B211]`). The phases share
the aims, so a keep-out that does not read the beams' power (an `ApertureKeepOut` at a
`peak_fraction`) gives every phase the same clearance, and the value is any one phase's
bit for bit; a keep-out at an absolute `intensity` widens with the power, and the
brightest phase sets the clearance.
"""
function clearance_floor(keepout; minimum::Union{Nothing,Real} = nothing, scale::Real = 1.0,
                         name = :clearance)
    hasmethod(keepout_clearance, Tuple{typeof(keepout),Vector{Beam{Float64}}}) ||
        throw(ArgumentError("$(typeof(keepout)) has no keepout_clearance method"))
    lim = minimum === nothing ? keepout.minimum : minimum
    return Constraint(name, _ClearanceValue(keepout); lower = lim, scale, labels = _beam_names)
end

# The constraint's value: the keep-out's clearance per beam of the evaluation, `-Inf` and
# `NaN` refused as by `keepout_penalty` (a violation of +Inf would make the augmented
# Lagrangian infinite without a message).
struct _ClearanceValue{K}
    keepout::K
end
function (c::_ClearanceValue)(ev::Evaluation)
    v = keepout_clearance(c.keepout, ev.beams)
    _check_clearances(c.keepout, ev.beams, v)
    return v
end

_beam_names(ev::Evaluation) = [b.name for b in ev.beams]

"""
    light_fraction(ev, names...) -> Real

All the laser light the named parts or components receive, from either side, as a
fraction of the delivered (effective, post-LPI) power:
[`power_fraction`](@ref)`(ev, names...) + `[`clipped_fraction`](@ref)`(ev, names...)`,
the absorbed power on their elements and the power striking their outside (back side).
The value of [`light_ceiling`](@ref) and [`light_floor`](@ref); differentiable, in the
evaluation's number type. A name not in the target is an `ArgumentError`.
"""
function light_fraction(ev::Evaluation, names::Symbol...)
    _check_names(ev.scenario.target, names)
    return power_fraction(ev, names...) + clipped_fraction(ev, names...)
end

"""
    light_ceiling(reference, names...; limit = light_fraction(reference, names...),
                  scale = max(limit, 1e-3), name = :light) -> Constraint
    light_ceiling(limit::Real, names...; scale = max(limit, 1e-3), name = :light) -> Constraint

A ceiling on the laser light the named parts or components receive: the
[`Constraint`](@ref) `light_fraction(ev, names...) ≤ limit` ([`light_fraction`](@ref):
absorbed plus striking their outside, a fraction of the delivered power), labeled by
the names. With a `reference` evaluation the limit is its light on the parts (no more
light on a window or a shield than an accepted design gives it); a number sets the
fraction directly.

`scale` sets the units of the violation: relative to the limit (the default) when the
limit is above 1e-3 of the delivered power, otherwise absolute in units of 1e-3, so
that a ceiling at or near zero (a part meant to stay dark) is met to 1e-6 of the
delivered power at `feasibility = 1e-3` (an engineering default). The fraction comes
from the deposition operator, which samples the spots' edges: compare a ray trace
([`raytrace_deposition`](@ref)) where the light on hardware must be known to better than
the operator's sampling (research notes, "Capture and clearance rules as optimization
quantities").
"""
light_ceiling(reference::Evaluation, names::Symbol...; limit = nothing, kwargs...) =
    _light_constraint(true, _light_limit(reference, names, limit), names; kwargs...)
light_ceiling(limit::Real, names::Symbol...; kwargs...) =
    _light_constraint(true, Float64(limit), names; kwargs...)

"""
    light_floor(reference, names...; limit = light_fraction(reference, names...),
                scale = max(limit, 1e-3), name = :light) -> Constraint
    light_floor(limit::Real, names...; scale = max(limit, 1e-3), name = :light) -> Constraint

A floor on the laser light the named parts or components receive,
`light_fraction(ev, names...) ≥ limit`: the capture by a part meant to be hit (a
backlighter foil, a landing surface), at least the `reference` design's or a given
fraction. Otherwise as [`light_ceiling`](@ref).
"""
light_floor(reference::Evaluation, names::Symbol...; limit = nothing, kwargs...) =
    _light_constraint(false, _light_limit(reference, names, limit), names; kwargs...)
light_floor(limit::Real, names::Symbol...; kwargs...) =
    _light_constraint(false, Float64(limit), names; kwargs...)

_light_limit(reference::Evaluation, names, limit) =
    limit === nothing ? Float64(ForwardDiff.value(light_fraction(reference, names...))) : Float64(limit)

function _light_constraint(upper::Bool, limit::Float64, names::Tuple; scale = nothing, name = :light)
    isempty(names) && throw(ArgumentError("name the parts (e.g. :cap_top)"))
    0 <= limit <= 1 || throw(ArgumentError("the limit is a fraction of the delivered power, got $limit"))
    s = scale === nothing ? max(limit, 1e-3) : scale
    f = _LightValue(names)
    lab = [join(string.(names), "+")]
    return upper ? Constraint(name, f; upper = limit, scale = s, labels = lab) :
                   Constraint(name, f; lower = limit, scale = s, labels = lab)
end

# The constraint's value: the light on the named parts.
struct _LightValue{N}
    names::NTuple{N,Symbol}
end
(l::_LightValue)(ev::Evaluation) = light_fraction(ev, l.names...)

# ---------------------------------------------------------------------------------------
# Mode shape against a reference
# ---------------------------------------------------------------------------------------

"""
    mode_departure(c, reference; lmax = 4) -> Vector
    mode_departure(ev, reference; lmax = 4) -> Vector

The departure of a drive shape from a reference shape per degree `l = 1:lmax`:
`δₗ = √((2l+1) Σₘ (c_lm/c_00 - r_lm/r_00)²)`, the [`mode_amplitude`](@ref) of the
difference of the two normalized shapes (`c` and `reference` are coefficient vectors
from [`flux_ylm`](@ref); with an evaluation, its own). For an axisymmetric difference
this is `|aₗ/a₀ - aₗ'/a₀'|`. Generic in the number type: the value of
[`mode_gate`](@ref).
"""
function mode_departure(c::AbstractVector, reference::AbstractVector; lmax::Integer = 4)
    n = (lmax + 1)^2
    length(c) >= n && length(reference) >= n ||
        throw(ArgumentError("the coefficient vectors need degrees through l = $lmax"))
    # a zero of the result type rather than the type itself: the comprehension's closure
    # stores a type as a `DataType` field, and Julia 1.10 then dispatches at run time
    z = zero(c[1] / c[1] - reference[1] / reference[1])
    return [_departure(c, reference, l, z) for l in 1:lmax]
end
mode_departure(ev::Evaluation, reference::AbstractVector; lmax::Integer = 4) =
    mode_departure(flux_ylm(ev; lmax), reference; lmax)

function _departure(c, r, l, z)
    s = z
    for m in -l:l
        i = mode_index(l, m)
        s += (c[i] / c[1] - r[i] / r[1])^2
    end
    d2 = (2l + 1) * s
    # √ at an exact match is a cone tip: its subgradient 0 instead of 0/0
    return d2 > 0 ? sqrt(d2) : zero(d2)
end

"""
    mode_gate(reference::Evaluation, gates; fraction = 1, floor = nothing, shift = nothing)
        -> Constraint
    mode_gate(reference::AbstractVector, gates; fraction = 1, floor = nothing, shift = nothing)
        -> Constraint

A gate per degree on the departure of the capsule's drive shape from a `reference`
shape: the [`Constraint`](@ref) `mode_departure(ev, reference; lmax) ≤ gates`, named
`:mode_gate`, one component per degree `l = 1:lmax` (labeled `l = 1`, …), each in units
of its own gate.

`gates` is a vector of departures (`gates[l]` for degree `l`, in the units of
[`mode_departure`](@ref): for an axisymmetric difference `|Δaₗ/a₀|`, so `0.005` at
`l = 2` is 0.5 % in P2/P0), or the result of [`mode_uncertainty`](@ref), whose
`rms_mode` gives the gate `√(2l+1) × rms_mode[l]`: the departure that random errors of
that size (e.g. NIF's 50 µm rms pointing) give each degree. That gate is a lower bound,
what pointing error alone produces: laser delivery and the hohlraum plasma move the
shape from shot to shot by more, in ways a view-factor model does not carry, so a
tolerance per degree in physical units, a choice of the design, is the usual gate. At
the pointing floor the gates can cost re-pointing and drive that a tolerance within
that shot-to-shot spread does not (research notes, "Mode gates in physical units").
With [`drive_objective`](@ref), the largest drive within the gates; with
[`shift_objective`](@ref), the least re-pointing that meets them.

With a reference [`Evaluation`](@ref) (e.g. the pointing before a bundle was dropped,
in a scenario that observes the capsule), the shape is its [`flux_ylm`](@ref)
coefficients, and the gate records the reference scenario's laser–plasma model (`lpi`)
and time: evaluating the gate on a scenario with another model, or the same model at
another time, throws an `ArgumentError`, as for [`irradiance_ceiling`](@ref) (the shape
of a drive without the scenario's CBET model is another shape, and the gates would hold
the design to it). The reference may be on another target (a hohlraum without the
design's windows). A coefficient vector records no model, so that method cannot check
it: the vector must come from an evaluation with the scenario's `lpi` and time.

**Keywords.**

- `fraction` (positive, default 1) scales a pointing floor: with `gates` a
  [`mode_uncertainty`](@ref) result the gate is `fraction × √(2l+1) × rms_mode[l]`
  (`fraction = 0.8`: 0.8 of the floor).
- `floor` (a [`mode_uncertainty`](@ref) result, or a vector of departures per degree)
  makes the floor a lower bound of physical gates: the gate per degree is the larger of
  `gates[l]` and `fraction × floor[l]`, so that a tolerance in physical units is never
  tighter than the departure pointing error alone produces (research notes, "Mode
  gates in physical units"). `gates` and `floor` cover the same degrees.
- `shift` moves the target shape away from the reference's: the departure is measured
  from the reference's coefficients with each listed mode shifted, the reference's
  laser–plasma model and time still checked. A shift is a pair `l => Δ` (the axisymmetric
  mode `(l, 0)`) or `(l, m) => Δ`, or a collection of them, with `Δ` in the units of
  [`mode_departure`](@ref): for `m = 0`, `Δ` is the change of `aₗ/a₀` in the Legendre
  expansion of the capsule flux `Σ aₗ Pₗ(cos θ)` (`2 => 0.024` sets the target's P2/P0
  0.024 above the reference's; with the orthonormal coefficients of [`flux_ylm`](@ref),
  `c_l0 → c_l0 + Δ c_00/√(2l+1)`), and for `m ≠ 0` the same change of the real
  coefficient `c_lm`, so that the shifted shape departs from the reference's by `|Δ|` in
  degree `l` ([`shift_modes`](@ref)). The drive `c_00` is unchanged.
"""
mode_gate(reference::Evaluation, gates; fraction::Real = 1, floor = nothing, shift = nothing) =
    _mode_gate_of(reference, _gate_vector(gates, fraction, floor), shift)
mode_gate(reference::AbstractVector, gates; fraction::Real = 1, floor = nothing, shift = nothing) =
    _mode_gate_of(reference, _gate_vector(gates, fraction, floor), shift)

function _mode_gate_of(reference::Evaluation, g, shift)
    c = flux_ylm(reference; lmax = max(length(g), 1))
    return _mode_gate(shift === nothing ? c : shift_modes(c, shift), g, _ReferenceModel(reference.scenario))
end
_mode_gate_of(reference::AbstractVector, g, shift) =
    _mode_gate(shift === nothing ? reference : shift_modes(reference, shift), g, nothing)

_gate_values(gates::AbstractVector{<:Real}) = gates
_gate_values(u::NamedTuple) = [sqrt(2l + 1) * r for (l, r) in enumerate(u.rms_mode)]
_gate_lmax(gates) = max(length(_gate_values(gates)), 1)

# The gates per degree from `mode_gate`'s arguments: the gates themselves, `fraction` of
# a floor, or the larger of physical gates and `fraction` of a floor.
function _gate_vector(gates, fraction::Real, floor)
    isfinite(fraction) && fraction > 0 || throw(ArgumentError("the fraction must be positive and finite, got $fraction"))
    if floor === nothing
        gates isa NamedTuple && return fraction == 1 ? _gate_values(gates) : fraction .* _gate_values(gates)
        fraction == 1 ||
            throw(ArgumentError("`fraction` scales a pointing floor: give the gates as a mode_uncertainty " *
                                "result, or the floor as `floor`"))
        return _gate_values(gates)
    end
    gates isa NamedTuple &&
        throw(ArgumentError("give the floor once: physical gates (a vector) with `floor`, or the floor as the gates"))
    g, f = _gate_values(gates), _gate_values(floor)
    length(g) == length(f) ||
        throw(DimensionMismatch("$(length(g)) gates for a floor of $(length(f)) degrees"))
    return max.(collect(Float64, g), fraction .* f)
end

"""
    shift_modes(c, shift) -> Vector{Float64}

The coefficient vector `c` ([`flux_ylm`](@ref): orthonormal real harmonics in
[`mode_index`](@ref) order) with each mode of `shift` moved, the target shape of
[`mode_gate`](@ref)`(...; shift)`: a pair `l => Δ` (the mode `(l, 0)`) or `(l, m) => Δ`,
or a collection of them, with `Δ` in the units of [`mode_departure`](@ref):
`c_lm → c_lm + Δ c_00/√(2l+1)`. For `m = 0`, `Δ` is the change of `aₗ/a₀` in the Legendre expansion
`Σ aₗ Pₗ(cos θ)` (since `aₗ/a₀ = √(2l+1) c_l0/c_00`); a single shift departs from `c` by
`|Δ|` in its degree. `c_00` (the drive) is unchanged, and the vector records no
laser–plasma model: `mode_gate(reference::Evaluation, gates; shift)` keeps the
reference's check.
"""
function shift_modes(c::AbstractVector{<:Real}, shift)
    out = collect(Float64, c)
    for (lm, Δ) in _shift_pairs(shift)
        l, m = lm isa Integer ? (Int(lm), 0) : (Int(lm[1]), Int(lm[2]))
        l >= 1 && -l <= m <= l || throw(ArgumentError("a shift needs a mode with l ≥ 1 and |m| ≤ l, got ($l, $m)"))
        isfinite(Δ) || throw(ArgumentError("the shift of mode ($l, $m) must be finite, got $Δ"))
        i = mode_index(l, m)
        i <= length(out) || throw(ArgumentError("the coefficients stop before mode ($l, $m)"))
        out[i] += Δ * out[1] / sqrt(2l + 1)
    end
    return out
end

_shift_pairs(p::Pair) = (p,)
_shift_pairs(ps) = ps

function _mode_gate(reference::AbstractVector, gates::AbstractVector{<:Real}, model)
    lmax = length(gates)
    lmax >= 1 || throw(ArgumentError("give a gate for at least l = 1"))
    all(g -> isfinite(g) && g > 0, gates) || throw(ArgumentError("gates must be positive and finite"))
    length(reference) >= (lmax + 1)^2 ||
        throw(ArgumentError("the reference needs coefficients through l = $lmax"))
    ref = collect(Float64, reference[1:(lmax+1)^2])
    return Constraint(:mode_gate, _GateValue(ref, lmax, model); upper = collect(Float64, gates),
                      labels = ["l = $l" for l in 1:lmax])
end

# The gate's value; `model` the reference's laser–plasma model and time (`nothing` for a
# coefficient vector: unchecked).
struct _GateValue{M<:Union{Nothing,_ReferenceModel}}
    reference::Vector{Float64}
    lmax::Int
    model::M
end
function (g::_GateValue)(ev::Evaluation)
    _check_model(g.model, "mode_gate", ev)
    return mode_departure(ev, g.reference; lmax = g.lmax)
end

"""
    drive_objective(; drive, clip = 10.0, missed = 10.0, direct = 10.0,
                    direct_names = nothing, keepout = nothing)

Objective for the largest drive: `-a₀/drive`, with `a₀` the [`mean_flux`](@ref) on the
capsule and `drive` (TW/mm²) its scale (e.g. `mean_flux(evaluate(scenario, x0))`, or a
reference design's), plus the `clip`, `missed`, `direct` and `keepout` penalties of
[`symmetry_objective`](@ref). Meant to run with constraints that hold the shape, such
as [`mode_gate`](@ref): alone it trades every symmetry for drive.
"""
function drive_objective(; drive::Real, clip = 10.0, missed = 10.0, direct = 10.0,
                         direct_names = nothing, keepout = nothing)
    drive > 0 || throw(ArgumentError("the drive scale must be positive"))
    return function (ev::Evaluation)
        J = -mean_flux(ev) / drive
        J += _penalties(ev; clip, missed, direct, direct_names)
        J += keepout_penalty(keepout, ev.beams)
        return J
    end
end

"""
    drive_floor(reference; fraction, scale = limit, name = :drive) -> Constraint

A floor on the capsule's drive: the [`Constraint`](@ref) `mean_flux(ev) ≥ fraction ×
mean_flux(reference)` ([`mean_flux`](@ref): the solid-angle mean of the incident flux on
the observed capsule, TW/mm²), with `reference` an [`Evaluation`](@ref) of a reference
design in a scenario that observes the capsule (the pointing before a bundle was
dropped, or the start). `fraction` (positive, no default) sets how much of the
reference's drive the design must keep: no value holds for every study, as a reference
with more beams than the design delivers more drive.

The constraint for a design whose power can fall: with power scales
([`PowerScales`](@ref)) or powers among the parameters, a ceiling on an absolute
quantity (an [`irradiance_ceiling`](@ref), a [`light_ceiling`](@ref)) is met by cutting
power, and the drive falls with it. The default `scale` is the limit, so the default
`feasibility = 1e-3` meets the floor to 0.1 % of the limit. As for
[`mode_gate`](@ref), the floor records the reference scenario's laser–plasma model and
time, and evaluating it with another throws an `ArgumentError`; the reference may be on
another target. [`drive_ceiling`](@ref) is the upper side, for a band.
"""
function drive_floor(reference::Evaluation; fraction::Real, scale::Union{Nothing,Real} = nothing,
                     name = :drive)
    lim = _drive_limit(reference, fraction)
    return Constraint(name, _DriveValue(_ReferenceModel(reference.scenario), false); lower = lim,
                      scale = something(scale, lim))
end

"""
    drive_ceiling(reference; fraction, scale = limit, name = :drive_ceiling) -> Constraint

A ceiling on the capsule's drive: the [`Constraint`](@ref) `mean_flux(ev) ≤ fraction ×
mean_flux(reference)`, the upper side of [`drive_floor`](@ref), with the same reference,
`scale` and check of the reference's laser–plasma model and time (an `ArgumentError` on
an evaluation with another; the reference may be on another target). With a floor it
holds the drive in a band about the reference's: `drive_floor(reference; fraction =
1.0)` and `drive_ceiling(reference; fraction = 1.02)` keep it within 0–2 % above, so
that power scales among the parameters cannot buy shape with drive beyond the band.
`fraction` (positive, no default) is a choice of the study, not a facility limit.
"""
function drive_ceiling(reference::Evaluation; fraction::Real, scale::Union{Nothing,Real} = nothing,
                       name = :drive_ceiling)
    lim = _drive_limit(reference, fraction)
    return Constraint(name, _DriveValue(_ReferenceModel(reference.scenario), true); upper = lim,
                      scale = something(scale, lim))
end

# The limit `fraction × mean_flux(reference)` of a drive floor or ceiling.
function _drive_limit(reference::Evaluation, fraction::Real)
    isfinite(fraction) && fraction > 0 || throw(ArgumentError("the fraction must be positive and finite"))
    a0 = Float64(ForwardDiff.value(mean_flux(reference)))
    a0 > 0 || throw(ArgumentError("the reference has no drive on the capsule"))
    return fraction * a0
end

# The floor's or ceiling's value: the capsule's mean flux, on the reference's laser–plasma
# model.
struct _DriveValue
    model::_ReferenceModel
    ceiling::Bool
end
function (d::_DriveValue)(ev::Evaluation)
    _check_model(d.model, d.ceiling ? "drive_ceiling" : "drive_floor", ev)
    return mean_flux(ev)
end
