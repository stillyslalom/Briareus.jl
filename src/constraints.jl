# This file defines inequality constraints for `optimize_pointing`, which wraps its
# L-BFGS-B in an augmented Lagrangian. It also defines ready-made constraints, such as a
# laser-irradiance ceiling on named parts and a per-degree gate on the capsule mode shape,
# and the drive objective used with them.

"""
    Constraint(name, f; upper = nothing, lower = nothing, scale = nothing, labels = nothing)

An inequality constraint for [`optimize_pointing`](@ref): `f(ev) ≤ upper` (a ceiling) or
`f(ev) ≥ lower` (a floor), with `f` any function of an [`Evaluation`](@ref) (its beams
`ev.beams` included), generic in the number type so that the automatic-differentiation
gradients of the optimizer pass through it. Give exactly one of `upper` and
`lower`.

`f` may return a number or a vector; a vector constraint has one component per entry,
each with a separate multiplier and row in the report, against a `limit` that is a number
(shared) or a vector of the same length. `labels` name the components in reports
(default `"1"`, `"2"`, …): a collection of names, or a function of the evaluation that
returns them, for a constraint whose components follow the evaluation (one per beam,
as in [`clearance_floor`](@ref)). Within one run the number of components is fixed; a
constraint with a shared limit and labels from the evaluation can be used for scenarios
with different numbers of beams, such as a stage on quads and then one on their beams.

`scale` (a number or one per component, positive) sets the units of the violation: the
optimizer works on `±(f(ev) - limit)/scale` and its `feasibility` tolerance is in those
units. By default `scale = |limit|`, so a violation is relative (1e-3: 0.1% of the
limit); a zero limit needs a `scale`.

Ready-made constraints: [`irradiance_ceiling`](@ref), [`mode_gate`](@ref),
[`clearance_floor`](@ref) (the clearance from a keep-out per beam), [`light_ceiling`](@ref)
and [`light_floor`](@ref) (the laser light on named parts), [`pulse_ceiling`](@ref)
(the limit of the laser on the pulses of a phased design). Any other quantity of the
beams only is `Constraint(:name, ev -> quantity(ev.beams); lower = minimum)`.
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

# Return the labels of the `n` components of a constraint at `ev`. Empty labels mean that
# the components are numbered.
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

# Return the components of the constraint normalized as violations,
# `±(value - limit)/scale`, in the number type of the evaluation. A violated component is
# positive.
function _violations(c::Constraint, ev)
    v = c.f(ev)
    s = c.upper ? 1 : -1
    return v isa Number ? [s * (v - c.limit[1]) / c.scale[1]] : s .* (v .- c.limit) ./ c.scale
end

"""
    check_constraints(ev, constraints; feasibility = 1e-3) -> Vector{ConstraintRow}

The value of each component of `constraints` (a [`Constraint`](@ref), or a tuple or
vector of them) at the evaluation `ev`, against its limit, as [`ConstraintRow`](@ref)s:
the report of `PointingResult.constraints`, for any design (a reference, the result of a
staged run). `feasibility` is the tolerance on the normalized violation that
separates `:violated`, `:active` and `:satisfied`.
"""
check_constraints(ev::Evaluation, c::Constraint; kwargs...) = check_constraints(ev, (c,); kwargs...)
check_constraints(ev::Evaluation, cs; feasibility::Real = 1e-3) =
    _constraint_rows(ev, Tuple(cs), nothing, Float64(feasibility))

function _constraint_rows(ev, cs::Tuple, multipliers, feasibility::Float64)
    rows = ConstraintRow[]
    k = 0
    for c in cs
        # Report the exact quantity where the function of the constraint is a smooth
        # stand-in for one (`_exact_value`), and the function value otherwise.
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

# Return the exact quantity for which the function of a constraint is a smooth stand-in,
# as a scalar at a Float64 evaluation. Return `nothing` where the function is itself the
# quantity. One example is the peak of an irradiance ceiling with `exact = true`, its
# default.
_exact_value(f, ev) = nothing

# ---------------------------------------------------------------------------------------
# Laser irradiance on named parts
# ---------------------------------------------------------------------------------------

# Return the elements of the named parts or components. A name not in the target is an
# error.
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

The largest absorbed laser irradiance (TW/mm²; `to_Wcm2` converts) over the elements of the
named parts or components. The irradiance of an element is its power over its area, as in
[`laser_flux`](@ref). The peak is exact over the elements, so it has a kink wherever it
moves to another element. It is the reported quantity of [`irradiance_ceiling`](@ref), whose
optimizer works on the smooth [`irradiance_norm`](@ref). The element average resolves a peak
only as finely as the mesh: a coarser mesh averages the hottest spot over more area and
reports a lower peak (on the dropped-bundle hohlraum, 8% lower at 0.45 mm than at 0.2 mm).

`weight`, a function of a beam, counts the deposition of each beam times `weight(beam)`
instead of all beams at their own powers: `0` leaves a beam out, so
`weight = b -> b.group in (:outer44, :outer50) ? 1.0 : 0.0` gives the irradiance of the
outer cones, and a picket energy over the beam power gives a fluence map (the deposition
is linear in the power of each beam). The weighted deposition is one more deposit per
evaluation, of the effective beams of the evaluation through the laser–plasma model of
the scenario. The function is called with the setup beams of the evaluation (their names,
groups, quads and ports).
"""
function peak_irradiance(ev::Evaluation, names::Symbol...; weight = nothing)
    t = ev.scenario.target
    return _peak(_counted_power(ev, weight), t, _named_elements(t, names))
end

_peak(P, t::Target, es) = maximum(P[e] / t.elem_area[e] for e in es)

# Return the element powers whose irradiance counts. Without a weight they are the
# deposition of the evaluation. With a weight they are the deposition of its effective
# beams, with the power of each beam multiplied by `weight` of its setup beam.
_counted_power(ev::Evaluation, ::Nothing) = ev.deposition.element
function _counted_power(ev::Evaluation, weight)
    s = ev.scenario
    eff = ev.effective
    beams = [remake(b; power = weight(ev.beams[i]) * b.power) for (i, b) in enumerate(eff)]
    return lpi_deposit(s.lpi, s.deposition, beams, s.time).element
end

"""
    irradiance_norm(ev, names...; p = 64, weight = nothing) -> Real

The area-weighted `p`-norm of the absorbed laser irradiance `I` (TW/mm²) over the elements
of the named parts or components, `(Σ Aₑ Iₑᵖ / Σ Aₑ)^(1/p)`. It is a smooth stand-in for the
[`peak_irradiance`](@ref). It lies below the peak and approaches it as `p` grows (the
power-mean inequality), and its gradient comes from every element near the peak instead of
only the hottest one. As an area mean it is a property of the irradiance field and not of
the elements, so it depends less on the mesh than the peak does. `weight` is as in
[`peak_irradiance`](@ref). Generic in the number type.
"""
function irradiance_norm(ev::Evaluation, names::Symbol...; p::Real = 64, weight = nothing)
    t = ev.scenario.target
    return _irradiance_norm(_counted_power(ev, weight), t, _named_elements(t, names), Float64(p))
end

function _irradiance_norm(P::AbstractVector{T}, t::Target, es, p::Float64) where {T}
    p >= 1 || throw(ArgumentError("p must be at least 1"))
    # Scale by the largest value against overflow; the scale is held constant, which is
    # valid because the norm is homogeneous.
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
components. The [`Constraint`](@ref) is named `:irradiance`, and by default it bounds the
[`peak_irradiance`](@ref) by `limit` (see `exact` below). Its function is the smooth
stand-in the optimizer works on, the [`irradiance_norm`](@ref) scaled to equal the peak at
the `reference` evaluation, `‖I‖ₚ × peak_ref/‖I‖ₚ,ref`. With the default `limit`, the peak
of the reference, the irradiance is held at that of the reference design, for example the
pointing before a bundle was dropped; a
number in TW/mm² sets another, measured on the same scale. The reference evaluation must be
on the same target as the optimized scenario (a `subset` of it shares the target), with the
same laser–plasma model (`lpi`, compared on every field, so that a model built alike passes)
at the same time; otherwise evaluating the ceiling throws an `ArgumentError`. A reference
evaluated without the CBET model of the scenario gives the peak and the stand-in scale of
another quantity, and the ceiling would hold the design to a limit on that quantity, so the
mismatch throws an error. Evaluate the reference in a scenario with the
same model, e.g. `evaluate(remake(reference_scenario; lpi = scenario.lpi), beams)`. `weight`
selects and weighs the beams that count, as in [`peak_irradiance`](@ref) (the outer cones
only, or picket energies for a fluence ceiling), for the reference and the design alike.

The ceiling is on an absolute irradiance: where the power can fall (power scales,
[`PowerScales`](@ref), or powers among the parameters) the optimizer can meet it by
cutting power, at a loss of drive (on a re-pointing study with a power scale per cone
at most 1, a large one). With a [`drive_floor`](@ref) beside it, the drive is held.

A ceiling on the peak wall irradiance at the value of the reference is an engineering
limit and is not a published criterion. The published scalings of the growth of the gold
bubble are written in terms of the picket fluence of the outer cones on the wall (research
notes), and a `weight` computes that fluence on the same machinery.

With `exact = true` (the default) the peak itself is held at the limit. The optimizer
works on the norm, but once the constrained run has converged it compares the
[`peak_irradiance`](@ref) with the limit; while the peak exceeds it by more than the
`feasibility` of the run (in units of `scale`), the working limit of the norm is scaled by
limit/peak and the run continues from the same design, multipliers and penalty. The
reported value ([`check_constraints`](@ref), `PointingResult.constraints`) is the peak
against the limit, and a converged run ends with the peak at or below the limit within
`feasibility`. The cost is the further outer iterations, when the peak ends above the
limit (research notes, "Constrained re-pointing").

`exact = false` bounds only the norm, and the norm is the reported value. At the limit the
hot spots of the field are as hot, in the `p`-norm mean, as those of the reference, but
the peak over the elements may exceed the limit by a margin that depends on where the run
ends. On the tutorial mesh the margin is 0.1% at `p = 64` and 3% at `p = 16`. On the
dropped-bundle test case (0.45 mm wall elements) it is 3% at `p = 64` for one of two
designs that meet the same gates, under another floating-point rounding of the same run
(research notes, "Constrained re-pointing").

Across meshes the ceiling is a ratio to the reference on the same mesh, which changes
less with the element size than the ratio of the peaks does (rms
difference 1.0% against 1.3% between 0.2 and 0.45 mm elements near the limit), while
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
    value = _CeilingValue(t, es, Float64(p), _peak_factor(pk, P, t, es, p), weight, exact,
                          _ReferenceModel(reference.scenario))
    return _ceiling_constraint(:irradiance, value, limit, pk, scale, "irradiance")
end

# The parts that the ceilings on an irradiance share (`irradiance_ceiling`,
# `picket_ceiling`). `_peak_factor` is the factor that puts the smooth p-norm of the
# irradiance `P/A` over the elements `es` on the scale of its reference value `ref`, so
# that the norm equals `ref` at the reference. `_ceiling_constraint` makes the constraint
# `value ≤ limit`, with the reference value as the default limit and the limit as the
# default scale.
_peak_factor(ref, P, t::Target, es, p) = ref / _irradiance_norm(P, t, es, Float64(p))

function _ceiling_constraint(name, value, limit, ref, scale, what::String)
    lim = limit === nothing ? ref : Float64(limit)
    lim > 0 || throw(ArgumentError("the $what limit must be positive"))
    return Constraint(name, value; upper = lim, scale = something(scale, lim))
end

# `_ReferenceModel` contains the laser–plasma model and time of a reference evaluation. A
# constraint with one throws on evaluations with another model (`_check_model`). The
# model type is a parameter, so the comparison with the model of a scenario, which is
# untyped there, is inferred after one type test.
struct _ReferenceModel{L}
    lpi::L
    time::Float64
end
_ReferenceModel(s::Scenario) = _ReferenceModel(s.lpi, s.time)

# The value of the ceiling is the irradiance norm on the peak scale of the reference. With
# `exact`, the reported value is the peak itself. `model` contains the laser–plasma model
# and time of the reference scenario, which every evaluation must share.
struct _CeilingValue{W,M<:_ReferenceModel}
    target::Target
    elements::Vector{Int}
    p::Float64
    factor::Float64
    weight::W
    exact::Bool
    model::M
end

function (c::_CeilingValue)(ev::Evaluation)
    _check_ceiling(c, ev)
    return c.factor * _irradiance_norm(_counted_power(ev, c.weight), c.target, c.elements, c.p)
end

function _check_ceiling(c::_CeilingValue, ev::Evaluation)
    ev.scenario.target === c.target ||
        throw(ArgumentError("irradiance_ceiling: the evaluation is on another target than the reference"))
    _check_model(c.model, "irradiance_ceiling", ev)
    return nothing
end

# Check that the reference evaluation of a constraint and the optimized scenario share the
# laser–plasma model and, with a model, its time. Otherwise throw an `ArgumentError` that
# names `who`. The models are compared on every field (`_same_model`) at each evaluation,
# and identical models are equal at once. The check keeps no state, so the evaluations of
# concurrent gradient chunks share nothing through it.
_check_model(::Nothing, who, ev::Evaluation) = nothing
function _check_model(m::_ReferenceModel{L}, who::String, ev::Evaluation) where {L}
    s = ev.scenario
    lpi = s.lpi
    lpi isa L && _same_model(lpi, m.lpi) ||
        throw(ArgumentError("$who: the reference was evaluated with another laser–plasma " *
                            "model than the scenario's ($(_model_name(m.lpi)) against " *
                            "$(_model_name(lpi))), so its limit bounds another quantity: evaluate the " *
                            "reference in a scenario with the same `lpi`, e.g. " *
                            "`evaluate(remake(reference_scenario; lpi = scenario.lpi), beams)`"))
    lpi === nothing || s.time == m.time ||
        throw(ArgumentError("$who: the reference was evaluated at time $(m.time) ns and the " *
                            "scenario's laser–plasma model at $(s.time) ns: evaluate the reference at the " *
                            "scenario's time (`remake(reference_scenario; time)`)"))
    return nothing
end

# The name of a model in the messages: its type, which for a tuple of models lists theirs.
# One method on an untyped argument, so the messages add no dispatch to the checks.
Base.@nospecializeinfer _model_name(@nospecialize(m)) = m === nothing ? "none" : string(typeof(m)::DataType)

# Test two laser–plasma models for structural equality, so that two models built alike
# are the same model. They are equal when they are identical, or when they have one type
# and equal fields, with arrays compared by their entries and dictionaries, strings,
# numbers and symbols by `isequal`. The field comparison is generated per type, so models
# of one concrete type with concretely typed fields are compared without dispatch.
function _same_model(a, b)
    a === b && return true
    typeof(a) === typeof(b) || return false
    return _same_fields(a, b)
end

@generated function _same_fields(a::T, b::T) where {T}
    T <: AbstractArray && return :(_same_entries(a, b))
    (T <: Union{AbstractDict,AbstractString,Number,Symbol} || fieldcount(T) == 0) && return :(isequal(a, b))
    ex = :(true)
    for i in fieldcount(T):-1:1
        ex = :(_same_model(getfield(a, $i), getfield(b, $i)) && $ex)
    end
    return ex
end

_same_entries(a::AbstractArray, b::AbstractArray) =
    length(a) == length(b) && all(i -> _same_model(a[i], b[i]), eachindex(a, b))

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

A floor on the worst-case clearance (mm) of each beam from a keep-out: the
[`Constraint`](@ref) `keepout_clearance(keepout, ev) .≥ minimum`, one component per beam
of the evaluation, labeled by the beam name. The keep-out is any
[`AbstractKeepOut`](@ref), such as an [`ApertureKeepOut`](@ref), a
[`SphereKeepOut`](@ref) or a [`HardwareKeepOut`](@ref). A beam the keep-out does not
apply to (clearance `Inf`) gives a component that is always met; a clearance of `-Inf` or
`NaN` throws an `ArgumentError` naming the keep-out and the beam, as in
[`keepout_penalty`](@ref). The `weight` of the keep-out is not used: the constraint
replaces its penalty ([`keepout_penalty`](@ref)), so pass the keep-out to the `keepout`
keyword of the objective or as a constraint, not both.

`scale` (mm) is the unit of the violation: with the default 1 mm and `feasibility =
1e-3`, a clearance within 1 µm of `minimum` is `:active` (an engineering default, finer
than the keep-out accuracy against [`clearance`](@ref), a few micrometers). Give
each of several keep-outs a distinct `name`.

Where the limits conflict, `scale` is also the weight of the keep-out against the other
constraints, because phase 1, the augmented Lagrangian and the least-violation design of
[`PointingResult`](@ref) all trade normalized violations. At 1 mm, a beam 0.1 mm inside
the keep-out weighs as much as a mode gate exceeded by 10% of itself or a ceiling 10% over
its limit, with one component per beam. A smaller scale, e.g. a pointing tolerance of
0.05 mm, makes the same shortfall weigh twenty times more. On a re-pointing study whose
gates conflicted with the keep-out and a ceiling, that change moved the violation from the
keep-out onto the gates and the ceiling, and each of them ended several times further
over (research notes, "Least violation, capped inner runs and phase 1's budget", § 6). The
scale should be chosen for the trade the design accepts.

The floor applies to each beam separately. In the augmented Lagrangian each beam
component has a separate multiplier, which stays zero for a beam clear of its minimum, and
the term of that beam then vanishes from the inner problem. Each term is as smooth as the
clearance of its beam. The smallest clearance over the beams would instead have a kink
wherever the closest beam changes (research notes, "Constrained re-pointing"). The
components follow the beams of the evaluation, so the same constraint can be used for a
stage on quads and a stage on their beams.

On a [`PhasedScenario`](@ref) the constraint takes the [`PerPhase`](@ref) evaluation
directly, without [`in_phase`](@ref): the component of each beam is its least clearance over
the phases, and its rows carry no phase in their names (`leh[B211]`). The phases share the
aims, so a keep-out that does not read the beam power (an `ApertureKeepOut` at a
`peak_fraction`) gives every phase the same clearance, identical to that of any one phase; a
keep-out at an absolute `intensity` widens with the power, and the brightest phase sets the
clearance.
"""
function clearance_floor(keepout::AbstractKeepOut; minimum::Union{Nothing,Real} = nothing, scale::Real = 1.0,
                         name = :clearance)
    lim = minimum === nothing ? keepout.minimum : minimum
    return Constraint(name, _ClearanceValue(keepout); lower = lim, scale, labels = _beam_names)
end

# The value of the constraint is the clearance from the keep-out per beam of the
# evaluation. `-Inf` and `NaN` throw, as in `keepout_penalty`, because a violation of +Inf
# would make the augmented Lagrangian infinite without a message.
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
number type of the evaluation. A name not in the target is an `ArgumentError`.
"""
function light_fraction(ev::Evaluation, names::Symbol...)
    _check_names(ev.scenario.target, names)
    return power_fraction(ev, names...) + clipped_fraction(ev, names...)
end

"""
    light_ceiling(reference, names...; limit = light_fraction(reference, names...),
                  scale = max(limit, 1e-3), name = :light) -> Constraint
    light_ceiling(limit::Real, names...; scale = max(limit, 1e-3), name = :light) -> Constraint

A ceiling on the laser light the named parts or components receive, the
[`Constraint`](@ref) `light_fraction(ev, names...) ≤ limit`, labeled by the names.
[`light_fraction`](@ref) is the light absorbed on the parts plus the light striking their
outside, as a fraction of the delivered power. With a `reference` evaluation the limit is
the light of the reference on the parts, so that a window or a shield receives no more
light than in an accepted design. A number sets the fraction directly.

`scale` sets the units of the violation: relative to the limit (the default) when the
limit is above 1e-3 of the delivered power, otherwise absolute in units of 1e-3, so
that a ceiling at or near zero (a part meant to stay dark) is met to 1e-6 of the
delivered power at `feasibility = 1e-3` (an engineering default). The fraction comes
from the deposition operator, which samples the spot edges: compare a ray trace
([`raytrace_deposition`](@ref)) where the light on hardware must be known to better than
the sampling of the operator (research notes, "Capture and clearance rules as optimization
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
`light_fraction(ev, names...) ≥ limit`: the capture by a part meant to be hit, such as a
backlighter foil or a landing surface, at least that of the `reference` design or a given
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

# The value of the constraint is the light on the named parts.
struct _LightValue{N}
    names::NTuple{N,Symbol}
end
(l::_LightValue)(ev::Evaluation) = light_fraction(ev, l.names...)

# ---------------------------------------------------------------------------------------
# Mode shape against a reference
# ---------------------------------------------------------------------------------------

"""
    mode_departure(c, reference; lmax = 4, m = :all) -> Vector
    mode_departure(ev, reference; lmax = 4, m = :all) -> Vector

The departure of a drive shape from a reference shape per degree `l = 1:lmax`:
`δₗ = √((2l+1) Σₘ (c_lm/c_00 - r_lm/r_00)²)`, the [`mode_amplitude`](@ref) of the
difference of the two normalized shapes (`c` and `reference` are coefficient vectors
from [`flux_ylm`](@ref); with an evaluation, those of the evaluation). For an axisymmetric difference
this is `|aₗ/a₀ - aₗ'/a₀'|`. Generic in the number type: the value of
[`mode_gate`](@ref).

`m = :zero` sums only the zonal order `m = 0`, `√(2l+1) |Δ(c_l0/c_00)|`, which is the
change of `aₗ/a₀` about the polar axis of the coefficients. `m = :nonzero` sums the other
orders, so that `δₗ² = δₗ(:zero)² + δₗ(:nonzero)²`, as in the split of
[`mode_rms`](@ref)'s `m`. A gate per degree bounds `δₗ`, and the split shows whether the
departure in a degree is axisymmetric. For example, `l = 1` with `m = ±1` is a sideways
shift of the drive from a dropped bundle.
"""
function mode_departure(c::AbstractVector, reference::AbstractVector; lmax::Integer = 4, m::Symbol = :all)
    n = (lmax + 1)^2
    length(c) >= n && length(reference) >= n ||
        throw(ArgumentError("the coefficient vectors need degrees through l = $lmax"))
    m in (:all, :zero, :nonzero) || throw(ArgumentError("m must be :all, :zero or :nonzero, got :$m"))
    # Pass a zero of the result type, not the type itself. The closure of the comprehension
    # would store a type as a `DataType` field, and Julia 1.10 would then dispatch at run
    # time.
    z = zero(c[1] / c[1] - reference[1] / reference[1])
    return [_departure(c, reference, l, z, m) for l in 1:lmax]
end
mode_departure(ev::Evaluation, reference::AbstractVector; lmax::Integer = 4, m::Symbol = :all) =
    mode_departure(flux_ylm(ev; lmax), reference; lmax, m)

# Return the departure of degree l over the orders that `m` selects (:all, :zero or
# :nonzero).
function _departure(c, r, l, z, m::Symbol)
    s = z
    for k in (m === :zero ? (0:0) : -l:l)
        m === :nonzero && k == 0 && continue
        i = mode_index(l, k)
        s += (c[i] / c[1] - r[i] / r[1])^2
    end
    d2 = (2l + 1) * s
    # At an exact match √ has a cone tip, so the subgradient 0 is returned instead of 0/0.
    return d2 > 0 ? sqrt(d2) : zero(d2)
end

"""
    mode_gate(reference::Evaluation, gates; fraction = 1, floor = nothing, shift = nothing)
        -> Constraint
    mode_gate(reference::AbstractVector, gates; fraction = 1, floor = nothing, shift = nothing)
        -> Constraint

A gate per degree on the departure of the capsule drive shape from a `reference` shape.
The [`Constraint`](@ref) `mode_departure(ev, reference; lmax) ≤ gates` is named
`:mode_gate` and has one component per degree `l = 1:lmax`, labeled `l = 1`, …, each in
units of its gate.

`gates` is either a vector of departures or the result of [`mode_uncertainty`](@ref). In
a vector, `gates[l]` is the gate for degree `l` in the units of [`mode_departure`](@ref),
which for an axisymmetric difference are `|Δaₗ/a₀|`; so `0.005` at `l = 2` is 0.5% in
P2/P0. For a [`mode_uncertainty`](@ref) result the gate is its `floor`,
`√(2l+1) × rms_mode[l]`, the departure that random errors of that size (e.g. NIF's 50 µm
rms pointing) give each degree. That gate is a lower bound, because it is what pointing
error alone produces. Laser delivery and the hohlraum plasma move the shape from shot to
shot by more, in ways a view-factor model does not carry. A tolerance per degree in
physical units, chosen for the design, is therefore the usual gate. At
the pointing floor the gates can cost re-pointing and drive that a tolerance within
that shot-to-shot spread does not (research notes, "Mode gates in physical units").
The share of laser delivery can be included in the floor through `mode_uncertainty`'s
`power`, an rms fraction per pointing unit. On the case of the dropped-bundle tutorial,
3% per quad raises the floor by 40% at `l = 1` and 21% at `l = 2`, and 5% raises it by 92%
and 51%. NIF's published beam-to-beam balance is better than 5% (research notes,
"Inner-run options, a wall-loss scenario warning and laser delivery figures", § 7).
With [`drive_objective`](@ref) the optimizer finds the largest drive within the gates;
with [`shift_objective`](@ref) it finds the least re-pointing that meets them.

With a reference [`Evaluation`](@ref) (e.g. the pointing before a bundle was dropped, in a
scenario that observes the capsule), the shape is its [`flux_ylm`](@ref) coefficients, and
the gate records the laser–plasma model (`lpi`) and time of the reference scenario:
evaluating the gate on a scenario with another model, or the same model at another time,
throws an `ArgumentError`, as for [`irradiance_ceiling`](@ref) (the shape of a drive without
the CBET model of the scenario is another shape, and the gates would hold the design to it).
The reference may be on another target (a hohlraum without the windows of the design). A
coefficient vector does not record a model, so that method cannot check it: the vector must
come from an evaluation with the `lpi` and time of the scenario.

**Keywords.**

- `fraction` (positive, default 1) scales a pointing floor. With `gates` a
  [`mode_uncertainty`](@ref) result, the gate is `fraction × √(2l+1) × rms_mode[l]`, that
  is `fraction × u.floor[l]`; `fraction = 0.8` gives 0.8 of the floor. A vector as long as
  the floor gives one multiple per degree; `fraction = [1, 2, 2, 2]` gives the floor at
  `l = 1` and twice the floor at `l = 2–4`.
- `floor` (a [`mode_uncertainty`](@ref) result, or a vector of departures per degree)
  makes the floor a lower bound of physical gates: the gate per degree is the larger of
  `gates[l]` and `fraction × floor[l]`, so that a tolerance in physical units is never
  tighter than the departure pointing error alone produces (research notes, "Mode
  gates in physical units"). `gates` and `floor` cover the same degrees.
- `shift` moves the target shape away from the reference shape: the departure is measured
  from the reference coefficients with each listed mode shifted, the laser–plasma model and
  time of the reference still checked. A shift is a pair `l => Δ` (the axisymmetric mode
  `(l, 0)`) or `(l, m) => Δ`, or a collection of them, with `Δ` in the units of
  [`mode_departure`](@ref): for `m = 0`, `Δ` is the change of `aₗ/a₀` in the Legendre
  expansion of the capsule flux `Σ aₗ Pₗ(cos θ)` (`2 => 0.024` sets the target P2/P0 0.024
  above that of the reference; with the orthonormal coefficients of [`flux_ylm`](@ref),
  `c_l0 → c_l0 + Δ c_00/√(2l+1)`), and for `m ≠ 0` the same change of the real coefficient
  `c_lm`, so that the shifted shape departs from the reference shape by `|Δ|` in degree `l`
  ([`shift_modes`](@ref)). The drive `c_00` is unchanged.
"""
mode_gate(reference::Evaluation, gates; fraction::Union{Real,AbstractVector{<:Real}} = 1, floor = nothing,
          shift = nothing) =
    _mode_gate_of(reference, _gate_vector(gates, fraction, floor), shift)
mode_gate(reference::AbstractVector, gates; fraction::Union{Real,AbstractVector{<:Real}} = 1, floor = nothing,
          shift = nothing) =
    _mode_gate_of(reference, _gate_vector(gates, fraction, floor), shift)

function _mode_gate_of(reference::Evaluation, g, shift)
    c = flux_ylm(reference; lmax = max(length(g), 1))
    return _mode_gate(shift === nothing ? c : shift_modes(c, shift), g, _ReferenceModel(reference.scenario))
end
_mode_gate_of(reference::AbstractVector, g, shift) =
    _mode_gate(shift === nothing ? reference : shift_modes(reference, shift), g, nothing)

_gate_values(gates::AbstractVector{<:Real}) = gates
_gate_values(u::NamedTuple) = _mode_floor(u.rms_mode)

# Return the gates per degree from the arguments of `mode_gate`. They are the gates
# themselves, `fraction` of a floor, or the larger of physical gates and `fraction` of a
# floor. `fraction` is one number or one per degree of the floor.
function _gate_vector(gates, fraction::Union{Real,AbstractVector{<:Real}}, floor)
    all(f -> isfinite(f) && f > 0, fraction) ||
        throw(ArgumentError("the fraction must be positive and finite, got $fraction"))
    if floor === nothing
        gates isa NamedTuple && return _scaled_floor(fraction, _gate_values(gates))
        all(==(1), fraction) ||
            throw(ArgumentError("`fraction` scales a pointing floor: give the gates as a mode_uncertainty " *
                                "result, or the floor as `floor`"))
        fraction isa Real || length(fraction) == length(_gate_values(gates)) ||
            throw(DimensionMismatch("$(length(fraction)) fractions for $(length(_gate_values(gates))) gates"))
        return _gate_values(gates)
    end
    gates isa NamedTuple &&
        throw(ArgumentError("give the floor once: physical gates (a vector) with `floor`, or the floor as the gates"))
    g, f = _gate_values(gates), _gate_values(floor)
    length(g) == length(f) ||
        throw(DimensionMismatch("$(length(g)) gates for a floor of $(length(f)) degrees"))
    return max.(collect(Float64, g), _scaled_floor(fraction, f))
end

# Multiply the floor `f` by `fraction`, which is one number or one per degree.
_scaled_floor(fraction::Real, f) = fraction == 1 ? f : fraction .* f
function _scaled_floor(fraction::AbstractVector{<:Real}, f)
    length(fraction) == length(f) ||
        throw(DimensionMismatch("$(length(fraction)) fractions for a floor of $(length(f)) degrees"))
    return Float64.(fraction) .* f
end

"""
    shift_modes(c, shift) -> Vector{Float64}

The coefficient vector `c` with each mode of `shift` moved, the target shape of
[`mode_gate`](@ref)`(...; shift)`. The coefficients are those of [`flux_ylm`](@ref),
orthonormal real harmonics in [`mode_index`](@ref) order. A shift is a pair `l => Δ` for
the mode `(l, 0)`, a pair `(l, m) => Δ`, or a collection of them, with `Δ` in the units of
[`mode_departure`](@ref):
`c_lm → c_lm + Δ c_00/√(2l+1)`. For `m = 0`, `Δ` is the change of `aₗ/a₀` in the Legendre
expansion `Σ aₗ Pₗ(cos θ)` (since `aₗ/a₀ = √(2l+1) c_l0/c_00`); a single shift departs from
`c` by `|Δ|` in its degree. `c_00` (the drive) is unchanged, and the vector does not record
a laser–plasma model: `mode_gate(reference::Evaluation, gates; shift)` still checks the
reference model.
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

# The value of a gate. `model` contains the laser–plasma model and time of the reference, or
# `nothing` for a coefficient vector, which is not checked.
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
that of a reference design), plus the `clip`, `missed`, `direct` and `keepout` penalties of
[`symmetry_objective`](@ref). The objective is meant to run with constraints that hold
the shape, such as [`mode_gate`](@ref). Without them it rewards drive at any cost to
symmetry.
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

A floor on the capsule drive, the [`Constraint`](@ref)
`mean_flux(ev) ≥ fraction × mean_flux(reference)`. [`mean_flux`](@ref) is the solid-angle
mean of the incident flux on the observed capsule (TW/mm²). `reference` is an
[`Evaluation`](@ref) of a reference design, such as the pointing before a bundle was
dropped or the start, in a scenario that observes the capsule. `fraction` (positive, no
default) sets how much of the reference drive the design must keep: no value suits every
study, as a reference with more beams than the design delivers more drive.

This floor is the constraint for a design whose power can fall. With power scales
([`PowerScales`](@ref)) or powers among the parameters, a ceiling on an absolute
quantity, such as an [`irradiance_ceiling`](@ref) or a [`light_ceiling`](@ref), is met by
cutting power, and the drive falls with it. The default `scale` is the limit, so the default
`feasibility = 1e-3` meets the floor to 0.1% of the limit. As for
[`mode_gate`](@ref), the floor records the laser–plasma model and time of the reference
scenario, and evaluating it with another throws an `ArgumentError`; the reference may be on
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

A ceiling on the capsule drive: the [`Constraint`](@ref) `mean_flux(ev) ≤ fraction ×
mean_flux(reference)`, the upper side of [`drive_floor`](@ref), with the same reference,
`scale` and check of the laser–plasma model and time of the reference (an `ArgumentError`
on an evaluation with another; the reference may be on another target). With a floor the
drive is held in a band about the reference drive: `drive_floor(reference; fraction =
1.0)` and `drive_ceiling(reference; fraction = 1.02)` keep it within 0–2% above, so
that power scales among the parameters cannot buy shape with drive beyond the band.
`fraction` (positive, no default) is a choice of the study and is not a facility limit.
"""
function drive_ceiling(reference::Evaluation; fraction::Real, scale::Union{Nothing,Real} = nothing,
                       name = :drive_ceiling)
    lim = _drive_limit(reference, fraction)
    return Constraint(name, _DriveValue(_ReferenceModel(reference.scenario), true); upper = lim,
                      scale = something(scale, lim))
end

# Return the limit `fraction × mean_flux(reference)` of a drive floor or ceiling.
function _drive_limit(reference::Evaluation, fraction::Real)
    isfinite(fraction) && fraction > 0 || throw(ArgumentError("the fraction must be positive and finite"))
    a0 = Float64(ForwardDiff.value(mean_flux(reference)))
    a0 > 0 || throw(ArgumentError("the reference has no drive on the capsule"))
    return fraction * a0
end

# The value of the floor or ceiling is the mean flux on the capsule, evaluated with the
# laser–plasma model of the reference.
struct _DriveValue{M<:_ReferenceModel}
    model::M
    ceiling::Bool
end
function (d::_DriveValue)(ev::Evaluation)
    _check_model(d.model, d.ceiling ? "drive_ceiling" : "drive_floor", ev)
    return mean_flux(ev)
end
