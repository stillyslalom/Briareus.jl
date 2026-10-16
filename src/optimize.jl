# Pointing parameterizations, objectives, and gradient-based optimization.

const _CACHED_LMAX = 12

const _FIELDS = (:z, :r, :phi, :defocus, :power, :spot)

"""
    GroupPointing(beams; by = :group, fields = (:z, :r, :defocus), mirror = true)

Maps a parameter vector to beams. Beams are grouped by `by` (`:group` — e.g. cone —,
`:quad`, `:beam`, or a function `beam -> key`); each group gets one value per field:

- `:z` — axial aim coordinate (mm). With `mirror`, beams from the lower hemisphere
  (port `z < 0`) aim at `-z`, so one value describes a top/bottom-symmetric cone.
- `:r` — radial aim offset (mm) toward the beam's own port azimuth.
- `:phi` — tangential aim offset (mm), counterclockwise about `z`.
- `:defocus` — focal-plane offset along the beam (mm, see `Beam`).
- `:power` — power per beam (TW).
- `:spot` — multiplier on the best-focus spot semi-axes.

Fields not listed keep each beam's template value.
"""
struct GroupPointing{K}
    template::Vector{Beam{Float64}}
    keys::Vector{K}
    members::Vector{Vector{Int}}
    fields::Vector{Symbol}
    mirror::Bool
end

function GroupPointing(beams::AbstractVector{<:Beam}; by = :group,
                       fields = (:z, :r, :defocus), mirror::Bool = true)
    for f in fields
        f in _FIELDS || throw(ArgumentError("unknown field $f; choose from $_FIELDS"))
    end
    key = by === :group ? (b -> b.group) : by === :quad ? (b -> b.quad) :
          by === :beam ? (b -> b.name) : by
    ks = map(key, beams)
    uk = unique(ks)
    members = [findall(==(k), ks) for k in uk]
    tmpl = [remake(b; aim = Float64.(ForwardDiff.value.(b.aim)),
                   power = Float64(ForwardDiff.value(b.power)),
                   defocus = Float64(ForwardDiff.value(b.defocus)),
                   spot = SuperGaussianSpot(Float64(ForwardDiff.value(b.spot.a)),
                                            Float64(ForwardDiff.value(b.spot.b));
                                            order = b.spot.order, angle = b.spot.angle))
            for b in beams]
    return GroupPointing(tmpl, uk, members, collect(Symbol, fields), mirror)
end

nparams(p::GroupPointing) = length(p.keys) * length(p.fields)

"Names `\"<group>.<field>\"` for each parameter."
parameter_names(p::GroupPointing) = [string(k, ".", f) for k in p.keys for f in p.fields]

_pidx(p::GroupPointing, g, fi) = (g - 1) * length(p.fields) + fi

_zsign(p::GroupPointing, b::Beam) = (p.mirror && b.port[3] < 0) ? -1.0 : 1.0
_radial(b::Beam) = (φ = atan(b.port[2], b.port[1]); (Vec3(cos(φ), sin(φ), 0), Vec3(-sin(φ), cos(φ), 0)))

"""
    initial_parameters(p::GroupPointing) -> Vector{Float64}

Parameters reproducing the template beams (group means where members differ).
"""
function initial_parameters(p::GroupPointing)
    x = zeros(nparams(p))
    for (g, mem) in enumerate(p.members), (fi, f) in enumerate(p.fields)
        vals = map(mem) do bi
            b = p.template[bi]
            er, eφ = _radial(b)
            f === :z ? _zsign(p, b) * b.aim[3] :
            f === :r ? dot(b.aim, er) :
            f === :phi ? dot(b.aim, eφ) :
            f === :defocus ? b.defocus :
            f === :power ? b.power : 1.0
        end
        x[_pidx(p, g, fi)] = sum(vals) / length(vals)
    end
    return x
end

"""
    apply(p::GroupPointing, x) -> Vector{Beam}

Beams for parameter vector `x` (element type follows `x`, so this is differentiable).
"""
function apply(p::GroupPointing, x::AbstractVector{T}) where {T}
    length(x) == nparams(p) || throw(DimensionMismatch("expected $(nparams(p)) parameters"))
    out = Vector{Beam{T}}(undef, length(p.template))
    for (g, mem) in enumerate(p.members)
        val(f) = (fi = findfirst(==(f), p.fields); fi === nothing ? nothing : x[_pidx(p, g, fi)])
        z, r, ϕ, dfc, pw, sp = val(:z), val(:r), val(:phi), val(:defocus), val(:power), val(:spot)
        for bi in mem
            b = p.template[bi]
            er, eφ = _radial(b)
            aim0 = b.aim
            rr = r === nothing ? dot(aim0, er) : r
            tt = ϕ === nothing ? dot(aim0, eφ) : ϕ
            zz = z === nothing ? aim0[3] : _zsign(p, b) * z
            aim = rr * er + tt * eφ + SVector{3,T}(0, 0, zz)
            spot = sp === nothing ? SuperGaussianSpot{T}(b.spot.a, b.spot.b, b.spot.order, b.spot.angle) :
                   SuperGaussianSpot{T}(sp * b.spot.a, sp * b.spot.b, b.spot.order, b.spot.angle)
            out[bi] = remake(b; aim, spot,
                             power = pw === nothing ? T(b.power) : pw,
                             defocus = dfc === nothing ? T(b.defocus) : dfc)
        end
    end
    return out
end

"""
    bounds(p::GroupPointing, x0; z = 1.0, r = 1.0, phi = 1.0, defocus = 2.0,
           power = 0.5, spot = 0.5) -> (lower, upper)

Box bounds around `x0`: `±z`, `±r`, `±phi`, `±defocus` mm, and relative `±power`,
`±spot` fractions.
"""
function bounds(p::GroupPointing, x0 = initial_parameters(p); z = 1.0, r = 1.0, phi = 1.0,
                defocus = 2.0, power = 0.5, spot = 0.5)
    lo, hi = copy(x0), copy(x0)
    for g in eachindex(p.members), (fi, f) in enumerate(p.fields)
        i = _pidx(p, g, fi)
        if f in (:power, :spot)
            frac = f === :power ? power : spot
            lo[i] = x0[i] * (1 - frac); hi[i] = x0[i] * (1 + frac)
        else
            δ = f === :z ? z : f === :r ? r : f === :phi ? phi : defocus
            lo[i] = x0[i] - δ; hi[i] = x0[i] + δ
        end
    end
    return lo, hi
end

# ---------------------------------------------------------------------------------------
# Scenario: precomputed operators + evaluation
# ---------------------------------------------------------------------------------------

"""
    Scenario(target, beams; observe = nothing, radiation = true, time = 0.0,
             parameterization = GroupPointing(beams), deposition_kw = (;), viewfactor_kw = (;))

Everything precomputed for fast repeated evaluation of a pointing design:
the `DepositionOperator`, and — when `radiation` and `observe` (component names whose
incident x-ray flux is the figure of merit, e.g. `:capsule`) are given — the view
factors and the linear flux response of the observed elements.
"""
struct Scenario{P}
    target::Target
    param::P
    deposition::DepositionOperator
    model::Union{Nothing,RadiationModel}
    observed::Vector{Int}
    response::Union{Nothing,Matrix{Float64}}
    moment_weights::Matrix{Float64}   # legendre_weights of the observed elements (lmax = 12)
end

function Scenario(t::Target, beams::AbstractVector{<:Beam}; observe = nothing,
                  radiation::Bool = observe !== nothing, time = 0.0,
                  parameterization = GroupPointing(beams), deposition_kw = (;),
                  viewfactor_kw = (;))
    op = DepositionOperator(t, beams; deposition_kw...)
    obs = observe === nothing ? Int[] : elements(t, (observe isa Symbol ? (observe,) : observe)...)
    model = radiation ? RadiationModel(t; viewfactor_kw...) : nothing
    R = (model !== nothing && !isempty(obs)) ? flux_response(model, obs; time) : nothing
    W = isempty(obs) ? zeros(0, 0) : legendre_weights(t, obs; lmax = _CACHED_LMAX)
    return Scenario(t, parameterization, op, model, obs, R, W)
end

function Base.show(io::IO, s::Scenario)
    print(io, "Scenario(", nelements(s.target), " elements, ", length(s.deposition.ports),
          " beams, ", s.param isa GroupPointing ? nparams(s.param) : "?", " parameters, ",
          length(s.observed), " observed elements)")
end

"""
    Evaluation

Result of `evaluate`: the `beams`, their `deposition`, and `flux` — incident x-ray flux
on the scenario's observed elements (TW/mm²; `nothing` without radiation).
"""
struct Evaluation{T,B}
    scenario::Scenario
    beams::Vector{B}
    deposition::Deposition{T}
    flux::Union{Nothing,Vector{T}}
end

"""
    evaluate(scenario, x) / evaluate(scenario, beams) -> Evaluation
"""
evaluate(s::Scenario, x::AbstractVector{<:Real}) = evaluate(s, apply(s.param, x))

function evaluate(s::Scenario, beams::AbstractVector{<:Beam})
    dep = deposit(s.deposition, beams)
    H = s.response === nothing ? nothing : s.response * dep.element
    return Evaluation(s, beams, dep, H)
end

"Fraction of delivered laser power striking the outside (back side) of any component."
clipped_fraction(ev::Evaluation) = sum(ev.deposition.backside) / ev.deposition.delivered
"Fraction of delivered power that strikes no surface at all."
missed_fraction(ev::Evaluation) = 1 - sum(ev.deposition.per_beam) / ev.deposition.delivered
"Fraction of delivered power absorbed by the named components."
power_fraction(ev::Evaluation, names::Symbol...) =
    sum(ev.deposition.element[e] for e in elements(ev.scenario.target, names...); init = zero(ev.deposition.delivered)) /
    ev.deposition.delivered

"""
    flux_moments(ev; lmax = 8, axis = ẑ, center = 0) -> Vector

Legendre moments `aₗ` of the observed incident flux (index `l+1`).
"""
function flux_moments(ev::Evaluation; lmax::Integer = 8, axis = Vec3(0, 0, 1),
                      center = zero(Vec3))
    s = ev.scenario
    ev.flux === nothing && throw(ArgumentError("scenario has no radiation response"))
    if lmax <= _CACHED_LMAX && Vec3(axis...) == Vec3(0, 0, 1) && Vec3(center...) == zero(Vec3)
        return view(s.moment_weights, 1:lmax+1, :) * ev.flux
    end
    return legendre_moments(_full(ev), s.target, s.observed; lmax, axis, center)
end

function _full(ev::Evaluation)
    ev.flux === nothing && throw(ArgumentError("scenario has no radiation response"))
    v = zeros(eltype(ev.flux), nelements(ev.scenario.target))
    v[ev.scenario.observed] .= ev.flux
    return v
end

"""
    symmetry_objective(; modes = (2 => 1.0, 4 => 1.0), odd = 1.0, lmax = 6,
                       clip = 10.0, missed = 10.0, direct = 10.0, direct_names = (:capsule,),
                       drive = nothing, drive_weight = 1.0, center = 0)

Objective for indirect-drive symmetry: `Σ wₗ (aₗ/a₀)²` over `modes`, plus `odd ×` the
sum of squared odd moments (top/bottom imbalance), plus `clip`, `missed` and `direct`
times the squared fractions of laser power clipped on the outside of the target,
missing it, or hitting `direct_names` (e.g. the capsule) directly. When `drive`
(TW/mm²) is given, adds `drive_weight × (a₀/drive - 1)²`.
"""
function symmetry_objective(; modes = (2 => 1.0, 4 => 1.0), odd = 1.0, lmax::Integer = 6,
                            clip = 10.0, missed = 10.0, direct = 10.0,
                            direct_names = (:capsule,), drive = nothing, drive_weight = 1.0,
                            center = zero(Vec3))
    return function (ev::Evaluation)
        a = flux_moments(ev; lmax, center)
        J = sum(w * (a[l+1] / a[1])^2 for (l, w) in modes)
        odd > 0 && (J += odd * sum((a[l+1] / a[1])^2 for l in 1:2:lmax))
        clip > 0 && (J += clip * clipped_fraction(ev)^2)
        missed > 0 && (J += missed * missed_fraction(ev)^2)
        if direct > 0 && all(n -> any(c -> c.name === n, ev.scenario.target.components), direct_names)
            J += direct * power_fraction(ev, direct_names...)^2
        end
        drive === nothing || (J += drive_weight * (a[1] / drive - 1)^2)
        return J
    end
end

"""
    irradiance_objective(name; center, radius, intensity = nothing, uniformity = 1.0,
                         spill = 1.0, intensity_weight = 1.0, side = 1)

Objective for a laser-irradiated foil (e.g. a backlighter): RMS nonuniformity of the
absorbed laser intensity over the elements of component `name` within `radius` of
`center` (on side `side`: 1 front, -1 back), plus `spill ×` the squared fraction of
delivered power landing outside that region, plus — if `intensity` (TW/mm²) is given —
`intensity_weight × (mean/intensity - 1)²`.
"""
function irradiance_objective(name::Symbol; center, radius, intensity = nothing,
                              uniformity = 1.0, spill = 1.0, intensity_weight = 1.0, side = 1)
    c0 = Vec3(center...)
    return function (ev::Evaluation)
        t = ev.scenario.target
        es = [e for e in elements(t, name) if t.elem_side[e] == side &&
              norm(t.elem_centroid[e] - c0) <= radius]
        isempty(es) && throw(ArgumentError("no elements of $name within radius"))
        P = ev.deposition.element[es]
        A = t.elem_area[es]
        J = uniformity * rms_nonuniformity(P ./ A, A)^2
        spill > 0 && (J += spill * (1 - sum(P) / ev.deposition.delivered)^2)
        if intensity !== nothing
            J += intensity_weight * (sum(P) / sum(A) / intensity - 1)^2
        end
        return J
    end
end

"""
    optimize_pointing(scenario, objective; x0, lower, upper, iterations = 200,
                      show_trace = false) -> NamedTuple

Minimize `objective(evaluate(scenario, x))` within box bounds by L-BFGS-B with
forward-mode automatic-differentiation gradients. Returns `(; x, beams, value,
initial_value, result)`.
"""
function optimize_pointing(s::Scenario, objective; x0 = initial_parameters(s.param),
                           lower = bounds(s.param, x0)[1], upper = bounds(s.param, x0)[2],
                           iterations::Integer = 200, show_trace::Bool = false,
                           g_tol = 1e-8)
    f(x) = objective(evaluate(s, x))
    cfg = ForwardDiff.GradientConfig(f, x0, ForwardDiff.Chunk(min(length(x0), 12)))
    g!(G, x) = ForwardDiff.gradient!(G, f, x, cfg)
    xs = clamp.(x0, lower, upper)
    res = Optim.optimize(f, g!, lower, upper, xs, Optim.LBFGSB(),
                         Optim.Options(; iterations, show_trace, g_abstol = g_tol))
    x = Optim.minimizer(res)
    return (; x, beams = apply(s.param, x), value = Optim.minimum(res),
            initial_value = f(xs), result = res)
end

"""
    report(ev; lmax = 6) -> nothing

Print a power accounting and flux-moment summary of an evaluation.
"""
function report(io::IO, ev::Evaluation; lmax::Integer = 6)
    d = ev.deposition
    P = ForwardDiff.value(d.delivered)
    @printf(io, "Delivered laser power: %.4g TW\n", P)
    cp = component_power(d, ev.scenario.target)
    missed = pop!(cp, :missed)
    for (k, v) in sort(collect(cp); by = first)
        v > 1e-6 * P && @printf(io, "  %-20s %9.4g TW  (%6.2f%%)\n", k, v, 100v / P)
    end
    # negative values are quadrature imbalance of the analytic deposition
    abs(missed) > 1e-6 * P &&
        @printf(io, "  %-20s %9.4g TW  (%6.2f%%)\n", missed > 0 ? "missed" : "imbalance", missed, 100missed / P)
    if ev.flux !== nothing
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
