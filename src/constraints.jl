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

function _constraint_rows(ev::Evaluation, cs::Tuple, multipliers, feasibility::Float64)
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
# irradiance ceiling's peak when it holds the peak (`exact = true`).
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

function (a::_Augmented)(ev::Evaluation)
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
        s::Scenario, @nospecialize(objective), @nospecialize(cs::Tuple), x0, xs, lower, upper;
        iterations::Integer, evaluations::Integer, show_trace::Bool, g_tol,
        @nospecialize(callback), feasibility::Float64, outer::Int, @nospecialize(state = nothing),
        deadline::Float64 = Inf)
    for c in cs
        c isa Constraint || throw(ArgumentError("constraints must be `Constraint`s, got $(typeof(c))"))
    end
    feasibility > 0 || throw(ArgumentError("feasibility must be positive"))
    ev = evaluate(s, xs)
    f0 = Float64(ForwardDiff.value(objective(ev)))
    g = _all_violations(cs, ev)
    lens = map(c -> length(_violations(c, ev)), cs)
    # the constraints as the inner problems see them: an exact ceiling's stand-in at a
    # working limit, lowered while the exact value ends above the limit
    work = cs
    if state === nothing
        λ = zeros(length(g))
        μ = clamp(10 * max(1.0, abs(f0)) / max(1.0, sum(abs2, max.(g, 0.0)) / 2), 1e-6, 1e6)
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
    x, ne, ni, k = copy(xs), 0, 0, 0
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
    rows = _constraint_rows(evx, cs, λ, feasibility)
    stopped = last.stopped
    return _pointing_result(s, x, Float64(ForwardDiff.value(objective(evx))), f0, ni, ne,
                            stopped, reason, last.result, lower, upper; constraints = rows,
                            outer = k, state, objective,
                            lagrangian = _Lagrangian(objective, work, _split_multipliers(λ, lens)))
end

# The Lagrangian f + Σᵢ λᵢ gᵢ of an objective and the constraints' normalized violations,
# for the bound sensitivities of a result.
struct _Lagrangian{O,C<:Tuple,L<:Tuple}
    objective::O
    constraints::C
    λ::L
end

function (a::_Lagrangian)(ev::Evaluation)
    J = a.objective(ev)
    terms = map((c, λ) -> _lagrangian_term(_violations(c, ev), λ), a.constraints, a.λ)
    return J + sum(terms)
end

function _lagrangian_term(g::AbstractVector{T}, λ::Vector{Float64}) where {T}
    s = zero(T)
    for i in eachindex(g, λ)
        s += λ[i] * g[i]
    end
    return s
end

# The working constraints with each exact ceiling whose exact value at `ev` exceeds its
# limit by more than `feasibility` (normalized) given its stand-in's working limit times
# limit/exact, or `nothing` when every exact value is within it. The stand-in tracks the
# exact value nearly in proportion near the limit (a p-norm against the peak), so the
# scaled limit lands the exact value near the limit; a few rounds settle the remainder.
function _lower_exact_limits(cs::Tuple, work::Tuple, ev::Evaluation, feasibility::Float64)
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
                       p = 64, weight = nothing, scale = limit, exact = false) -> Constraint

A ceiling on the absorbed laser irradiance over the elements of the named parts or
components: a [`Constraint`](@ref) named `:irradiance` whose value is the smooth
[`irradiance_norm`](@ref) scaled to equal the [`peak_irradiance`](@ref) at the
`reference` evaluation, `‖I‖ₚ × peak_ref/‖I‖ₚ,ref ≤ limit`. With the default `limit`
(the reference's peak) the ceiling holds the irradiance at the reference design's, the
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

A ceiling on the wall's peak irradiance at the reference's value is an engineering
limit, not a published criterion: the published scalings of the gold bubble's growth
use the outer cones' picket fluence on the wall (research notes), which `weight`
expresses on the same machinery.

The constraint bounds the norm, not the peak itself: at the limit the field's hot spots
are as hot, in the `p`-norm mean, as the reference's, while the peak over the elements
may exceed the limit, by a margin that depends on where the run ends: on the tutorial's
mesh 0.1 % at `p = 64` and 3 % at `p = 16`, and on the dropped-bundle test case (0.45 mm
wall elements) 3 % at `p = 64` for one of two designs that meet the same gates, under
another floating-point rounding of the same run (research notes, "Constrained
re-pointing").

`exact = true` holds the peak itself. The optimizer still works on the norm, but once
the constrained run has converged it compares the [`peak_irradiance`](@ref) with the
limit; while the peak exceeds it by more than the run's `feasibility` (in units of
`scale`), the norm's working limit is scaled by limit/peak and the run continues from
the same design, multipliers and penalty. The reported value
([`check_constraints`](@ref), `PointingResult.constraints`) is then the peak against
the limit, and a converged run ends with the peak at or below the limit within
`feasibility`. The cost is the further outer iterations (research notes, "Constrained
re-pointing").

Across meshes the ceiling is a ratio to the reference on the same
mesh, which changes less with the element size than the peak's ratio does (rms
difference 1.0 % against 1.3 % between 0.2 and 0.45 mm elements near the limit), while
the peak itself falls on a coarser mesh.
"""
function irradiance_ceiling(reference::Evaluation, names::Symbol...;
                            limit::Union{Nothing,Real} = nothing, p::Real = 64,
                            weight = nothing, scale::Union{Nothing,Real} = nothing,
                            exact::Bool = false)
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
    lpi = s.lpi
    if lpi !== c.checked[]
        _same_model(lpi, c.lpi) ||
            throw(ArgumentError("irradiance_ceiling: the reference was evaluated with another laser–plasma " *
                                "model than the scenario's ($(_model_name(c.lpi)) against " *
                                "$(_model_name(lpi))), so its limit bounds another quantity: evaluate the " *
                                "reference in a scenario with the same `lpi`, e.g. " *
                                "`evaluate(remake(reference_scenario; lpi = scenario.lpi), beams)`"))
        c.checked[] = lpi
    end
    lpi === nothing || s.time == c.time ||
        throw(ArgumentError("irradiance_ceiling: the reference was evaluated at time $(c.time) ns and the " *
                            "scenario's laser–plasma model at $(s.time) ns: evaluate the reference at the " *
                            "scenario's time (`remake(reference_scenario; time)`)"))
    return nothing
end

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
apply to (clearance `Inf`) gives a component that always holds. The keep-out's `weight`
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
"""
function clearance_floor(keepout; minimum::Union{Nothing,Real} = nothing, scale::Real = 1.0,
                         name = :clearance)
    hasmethod(keepout_clearance, Tuple{typeof(keepout),Vector{Beam{Float64}}}) ||
        throw(ArgumentError("$(typeof(keepout)) has no keepout_clearance method"))
    lim = minimum === nothing ? keepout.minimum : minimum
    return Constraint(name, _ClearanceValue(keepout); lower = lim, scale, labels = _beam_names)
end

# The constraint's value: the keep-out's clearance per beam of the evaluation.
struct _ClearanceValue{K}
    keepout::K
end
(c::_ClearanceValue)(ev::Evaluation) = keepout_clearance(c.keepout, ev.beams)

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
    mode_gate(reference, gates) -> Constraint

A gate per degree on the departure of the capsule's drive shape from a `reference`
shape (a coefficient vector from [`flux_ylm`](@ref)): the [`Constraint`](@ref)
`mode_departure(ev, reference; lmax) ≤ gates`, named `:mode_gate`, one component per
degree `l = 1:lmax` (labeled `l = 1`, …), each in units of its own gate.

`gates` is a vector of departures (`gates[l]` for degree `l`), or the result of
[`mode_uncertainty`](@ref), whose `rms_mode` gives the gate `√(2l+1) × rms_mode[l]`:
the departure that random errors of that size (e.g. NIF's 50 µm rms pointing) give
each degree, so that a design within the gates differs from the reference by no more
than the shot-to-shot spread. With [`drive_objective`](@ref), the largest drive within
the gates.
"""
function mode_gate(reference::AbstractVector, gates::AbstractVector{<:Real})
    lmax = length(gates)
    lmax >= 1 || throw(ArgumentError("give a gate for at least l = 1"))
    all(g -> isfinite(g) && g > 0, gates) || throw(ArgumentError("gates must be positive and finite"))
    length(reference) >= (lmax + 1)^2 ||
        throw(ArgumentError("the reference needs coefficients through l = $lmax"))
    ref = collect(Float64, reference[1:(lmax+1)^2])
    return Constraint(:mode_gate, _GateValue(ref, lmax); upper = collect(Float64, gates),
                      labels = ["l = $l" for l in 1:lmax])
end
mode_gate(reference::AbstractVector, u::NamedTuple) =
    mode_gate(reference, [sqrt(2l + 1) * r for (l, r) in enumerate(u.rms_mode)])

struct _GateValue
    reference::Vector{Float64}
    lmax::Int
end
(g::_GateValue)(ev::Evaluation) = mode_departure(ev, g.reference; lmax = g.lmax)

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
