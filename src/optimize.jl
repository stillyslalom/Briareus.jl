# Pointing parameterizations, objectives, and gradient-based optimization.

const _CACHED_LMAX = 12

const _FIELDS = (:z, :r, :phi, :defocus, :power, :spot, :dlambda)

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
- `:dlambda` — wavelength detuning (Å), used by laser–plasma models.

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
    by === :quad && any(b -> isempty(b.quad), beams) &&
        throw(ArgumentError("by = :quad needs quad names on every beam (OMEGA beams have none)"))
    key = by === :group ? (b -> b.group) : by === :quad ? (b -> b.quad) :
          by === :beam ? (b -> b.name) : by
    ks = map(key, beams)
    uk = unique(ks)
    members = [findall(==(k), ks) for k in uk]
    tmpl = [remake(b; aim = Float64.(ForwardDiff.value.(b.aim)),
                   power = Float64(ForwardDiff.value(b.power)),
                   defocus = Float64(ForwardDiff.value(b.defocus)),
                   spot = with_axes(b.spot, Float64(ForwardDiff.value(b.spot.a)),
                                    Float64(ForwardDiff.value(b.spot.b))),
                   dlambda = Float64(ForwardDiff.value(b.dlambda)))
            for b in beams]
    return GroupPointing(tmpl, uk, members, collect(Symbol, fields), mirror)
end

"Number of parameters of a `GroupPointing` (groups × fields)."
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
            f === :power ? b.power :
            f === :dlambda ? b.dlambda : 1.0
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
        dl = val(:dlambda)
        for bi in mem
            b = p.template[bi]
            er, eφ = _radial(b)
            aim0 = b.aim
            rr = r === nothing ? dot(aim0, er) : r
            tt = ϕ === nothing ? dot(aim0, eφ) : ϕ
            zz = z === nothing ? aim0[3] : _zsign(p, b) * z
            aim = rr * er + tt * eφ + SVector{3,T}(0, 0, zz)
            spot = sp === nothing ? _retype(b.spot, T) :
                   _retype(with_axes(b.spot, sp * b.spot.a, sp * b.spot.b), T)
            out[bi] = remake(b; aim, spot,
                             power = pw === nothing ? T(b.power) : pw,
                             defocus = dfc === nothing ? T(b.defocus) : dfc,
                             dlambda = dl === nothing ? T(b.dlambda) : dl)
        end
    end
    return out
end

"""
    bounds(p::GroupPointing, x0; z = 1.0, r = 1.0, phi = 1.0, defocus = 2.0,
           power = 0.5, spot = 0.5, dlambda = 2.0, power_min = 0.0,
           power_max = Inf) -> (lower, upper)

Box bounds around `x0`: `±z`, `±r`, `±phi`, `±defocus` mm, relative `±power`, `±spot`
fractions, and `±dlambda` Å. Power bounds are further clipped to the absolute limits
`power_min ≤ P ≤ power_max` (TW per beam, e.g. a facility's per-quad capability).
"""
function bounds(p::GroupPointing, x0 = initial_parameters(p); z = 1.0, r = 1.0, phi = 1.0,
                defocus = 2.0, power = 0.5, spot = 0.5, dlambda = 2.0, power_min = 0.0,
                power_max = Inf)
    lo, hi = copy(x0), copy(x0)
    for g in eachindex(p.members), (fi, f) in enumerate(p.fields)
        i = _pidx(p, g, fi)
        if f === :power
            lo[i] = max(x0[i] * (1 - power), power_min)
            hi[i] = min(x0[i] * (1 + power), power_max)
            lo[i] <= hi[i] || throw(ArgumentError("empty power bounds for group $(p.keys[g])"))
        elseif f === :spot
            lo[i] = x0[i] * (1 - spot); hi[i] = x0[i] * (1 + spot)
        else
            δ = f === :z ? z : f === :r ? r : f === :phi ? phi : f === :dlambda ? dlambda : defocus
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
             parameterization = GroupPointing(beams), lpi = nothing,
             deposition_kw = (;), viewfactor_kw = (;))

Everything precomputed for fast repeated evaluation of a pointing design:
the `DepositionOperator`, and — when `radiation` and `observe` (component names whose
incident x-ray flux is the figure of merit, e.g. `:capsule`) are given — the view
factors and the linear flux response of the observed elements. `lpi` is a laser–plasma
model (`Backscatter`, `CalibratedTransfer`, `CBET`, or a tuple of them) applied at `time`
between the setup beams and deposition.
"""
struct Scenario{P}
    target::Target
    param::P
    deposition::DepositionOperator
    model::Union{Nothing,RadiationModel}
    observed::Vector{Int}
    response::Union{Nothing,Matrix{Float64}}
    moment_weights::Matrix{Float64}   # legendre_weights of the observed elements (lmax = 12)
    ylm_weights::Matrix{Float64}      # ylm_weights of the observed elements (lmax = 6)
    lpi::Any                          # laser–plasma model (see apply_lpi)
    time::Float64
end

const _CACHED_YLM_LMAX = 6

function Scenario(t::Target, beams::AbstractVector{<:Beam}; observe = nothing,
                  radiation::Bool = observe !== nothing, time = 0.0,
                  parameterization = GroupPointing(beams), deposition_kw = (;),
                  viewfactor_kw = (;), lpi = nothing)
    op = DepositionOperator(t, beams; deposition_kw...)
    obs = observe === nothing ? Int[] : elements(t, (observe isa Symbol ? (observe,) : observe)...)
    model = radiation ? RadiationModel(t; viewfactor_kw...) : nothing
    R = (model !== nothing && !isempty(obs)) ? flux_response(model, obs; time) : nothing
    W = isempty(obs) ? zeros(0, 0) : legendre_weights(t, obs; lmax = _CACHED_LMAX)
    Y = isempty(obs) ? zeros(0, 0) : ylm_weights(t, obs; lmax = _CACHED_YLM_LMAX)
    return Scenario(t, parameterization, op, model, obs, R, W, Y, lpi, Float64(time))
end

"""
    subset(scenario, keep; parameterization = GroupPointing(kept beams; same fields))

Scenario for the beams `keep` (indices into the scenario's beams) — e.g. after
dropping a bundle — reusing the view factors, flux response and mode weights; only the
deposition operator is restricted. The kept beams' template values (pointing, power)
come from the scenario's parameterization; the default new parameterization keeps its
fields and mirroring and groups by `group` (pass one to group by quad or beam).
"""
function subset(s::Scenario, keep::AbstractVector{<:Integer}; parameterization = nothing)
    s.param isa GroupPointing || parameterization !== nothing ||
        throw(ArgumentError("give a parameterization for a non-GroupPointing scenario"))
    beams = s.param.template[keep]
    param = parameterization === nothing ?
            GroupPointing(beams; fields = Tuple(s.param.fields), mirror = s.param.mirror) :
            parameterization
    return Scenario(s.target, param, s.deposition[keep], s.model, s.observed, s.response,
                    s.moment_weights, s.ylm_weights, s.lpi, s.time)
end

function Base.show(io::IO, s::Scenario)
    print(io, "Scenario(", nelements(s.target), " elements, ", length(s.deposition.ports),
          " beams, ", s.param isa GroupPointing ? nparams(s.param) : "?", " parameters, ",
          length(s.observed), " observed elements)")
end

"""
    Evaluation

Result of `evaluate`: the setup `beams`, the `effective` beams after the scenario's
laser–plasma model, their `deposition`, and `flux` — incident x-ray flux on the
scenario's observed elements (TW/mm²; `nothing` without radiation). Deposition
fractions are relative to effective (post-LPI) power.
"""
struct Evaluation{T,B,E}
    scenario::Scenario
    beams::Vector{B}
    effective::Vector{E}
    deposition::Deposition{T}
    flux::Union{Nothing,Vector{T}}
end

"""
    evaluate(scenario, x) / evaluate(scenario, beams) -> Evaluation
"""
evaluate(s::Scenario, x::AbstractVector{<:Real}) = evaluate(s, apply(s.param, x))

function evaluate(s::Scenario, beams::AbstractVector{<:Beam})
    eff = apply_lpi(s.lpi, beams, s.time)
    dep = deposit(s.deposition, eff)
    H = s.response === nothing ? nothing : s.response * dep.element
    return Evaluation(s, beams, eff, dep, H)
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
                       drive = nothing, drive_weight = 1.0, center = 0, keepout = nothing)

Objective for indirect-drive symmetry: `Σ wₗ (aₗ/a₀)²` over `modes`, plus `odd ×` the
sum of squared odd moments (top/bottom imbalance), plus `clip`, `missed` and `direct`
times the squared fractions of laser power clipped on the outside of the target,
missing it, or hitting `direct_names` (e.g. the capsule) directly. When `drive`
(TW/mm²) is given, adds `drive_weight × (a₀/drive - 1)²`. A `keepout` (an
`ApertureKeepOut`) adds its penalty on the beams' worst-case clearance from the
entrance-hole rims.
"""
function symmetry_objective(; modes = (2 => 1.0, 4 => 1.0), odd = 1.0, lmax::Integer = 6,
                            clip = 10.0, missed = 10.0, direct = 10.0,
                            direct_names = (:capsule,), drive = nothing, drive_weight = 1.0,
                            center = zero(Vec3), keepout = nothing)
    return function (ev::Evaluation)
        a = flux_moments(ev; lmax, center)
        J = sum(w * (a[l+1] / a[1])^2 for (l, w) in modes)
        odd > 0 && (J += odd * sum((a[l+1] / a[1])^2 for l in 1:2:lmax))
        J += _penalties(ev; clip, missed, direct, direct_names)
        J += _keepout_penalty(keepout, ev.beams)
        drive === nothing || (J += drive_weight * (a[1] / drive - 1)^2)
        return J
    end
end

function _penalties(ev::Evaluation; clip, missed, direct, direct_names)
    J = zero(ev.deposition.delivered)
    clip > 0 && (J += clip * clipped_fraction(ev)^2)
    missed > 0 && (J += missed * missed_fraction(ev)^2)
    if direct > 0 && all(n -> any(c -> c.name === n, ev.scenario.target.components), direct_names)
        J += direct * power_fraction(ev, direct_names...)^2
    end
    return J
end

"""
    flux_ylm(ev; lmax = 6) -> Vector

Real spherical-harmonic coefficients of the observed incident flux in the flat
`mode_index` order (polar axis ẑ, azimuth from +x about TCC), from the scenario's cached
`ylm_weights`. Use `ylm_matrix(c, lmax)` for the `C[l+1, lmax+1+m]` layout.
"""
function flux_ylm(ev::Evaluation; lmax::Integer = _CACHED_YLM_LMAX)
    ev.flux === nothing && throw(ArgumentError("scenario has no radiation response"))
    lmax <= _CACHED_YLM_LMAX || throw(ArgumentError("cached Y_lm weights reach lmax = $_CACHED_YLM_LMAX"))
    return view(ev.scenario.ylm_weights, 1:(lmax+1)^2, :) * ev.flux
end

"""
    asymmetry_objective(; modes = (1 => 1.0, 2 => 1.0, 3 => 1.0, 4 => 1.0), components = (),
                        clip = 10.0, missed = 10.0, direct = 10.0, direct_names = (:capsule,),
                        drive = nothing, drive_weight = 1.0, keepout = nothing)

Objective for three-dimensional drive symmetry. For each `l => w` in `modes` it adds
`w × relative_mode_power(c, l)` — the rotation-invariant power of all `2l+1` real
harmonics of degree `l` over the mean (`(aₗ/a₀)²/(2l+1)` for an axisymmetric field) — and
for each `(l, m) => w` in `components` adds `w × (c_lm/c_00)²`, so a particular
azimuthal mode (e.g. the `l = 1, m = ±1` drift of a dropped bundle) can be weighted
separately. Penalties for clipping, missing and direct capsule illumination and the
`drive` target and `keepout` are as in `symmetry_objective`, with `drive` compared to the mean flux
`c_00/√(4π)`.
"""
function asymmetry_objective(; modes = (1 => 1.0, 2 => 1.0, 3 => 1.0, 4 => 1.0),
                             components = (), clip = 10.0, missed = 10.0, direct = 10.0,
                             direct_names = (:capsule,), drive = nothing, drive_weight = 1.0,
                             keepout = nothing)
    lmax = maximum(vcat([l for (l, _) in modes], [lm[1] for (lm, _) in components], 0))
    return function (ev::Evaluation)
        c = flux_ylm(ev; lmax = max(lmax, 1))
        J = zero(c[1])
        for (l, w) in modes
            J += w * relative_mode_power(c, l)
        end
        for ((l, m), w) in components
            J += w * (c[mode_index(l, m)] / c[1])^2
        end
        J += _penalties(ev; clip, missed, direct, direct_names)
        J += _keepout_penalty(keepout, ev.beams)
        drive === nothing || (J += drive_weight * (c[1] / sqrt(4π) / drive - 1)^2)
        return J
    end
end

"""
    influence_matrix(scenario, beams; lmax = 6) -> Matrix

Mode response of the observed flux to each beam: column `b` is the real-harmonic
coefficient vector (`mode_index` order) produced by beam `b` alone at 1 TW, with the
beam's pointing, spot and focus. Without a laser–plasma model the flux is linear in
beam powers, so the coefficients of any power assignment `P` are `M * P` — the basis of
power-rebalancing solves and of tables of what losing each beam costs.
"""
function influence_matrix(s::Scenario, beams::AbstractVector{<:Beam};
                          lmax::Integer = _CACHED_YLM_LMAX)
    s.response === nothing && throw(ArgumentError("scenario has no radiation response"))
    nb = length(beams)
    M = zeros((lmax + 1)^2, nb)
    Wy = view(s.ylm_weights, 1:(lmax+1)^2, :)
    for i in 1:nb
        unit = [remake(b; power = j == i ? 1.0 : 0.0) for (j, b) in enumerate(beams)]
        dep = deposit(s.deposition, unit)
        M[:, i] = Wy * (s.response * dep.element)
    end
    return M
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
                      show_trace = false, g_tol = 1e-8, callback = nothing) -> NamedTuple

Minimize `objective(evaluate(scenario, x))` within box bounds by L-BFGS-B with
forward-mode automatic-differentiation gradients. Returns `(; x, beams, value,
initial_value, iterations, evaluations, stopped, result)`.

`callback`, if given, is called after every (non-derivative) objective evaluation with
`(; evaluations, value, best_value, best_x)`; returning `true` stops the optimization,
which then returns the best design found so far with `stopped = true` and `result =
nothing`. Use it for progress reports and cancellation.
"""
function optimize_pointing(s::Scenario, objective; x0 = initial_parameters(s.param),
                           lower = bounds(s.param, x0)[1], upper = bounds(s.param, x0)[2],
                           iterations::Integer = 200, show_trace::Bool = false,
                           g_tol = 1e-8, callback = nothing)
    xs = clamp.(x0, lower, upper)
    best_x, best_value, count = copy(xs), Inf, Ref(0)
    function f(x)
        J = objective(evaluate(s, x))
        if eltype(x) <: AbstractFloat
            count[] += 1
            J < best_value && ((best_value, best_x) = (J, copy(x)))
            callback === nothing ||
                callback((; evaluations = count[], value = J, best_value, best_x)) === true &&
                throw(_StopOptimization())
        end
        return J
    end
    cfg = ForwardDiff.GradientConfig(f, x0, ForwardDiff.Chunk(min(length(x0), 12)))
    g!(G, x) = ForwardDiff.gradient!(G, f, x, cfg)
    initial_value = objective(evaluate(s, xs))
    res = try
        Optim.optimize(f, g!, lower, upper, xs, Optim.LBFGSB(),
                       Optim.Options(; iterations, show_trace, g_abstol = g_tol))
    catch err
        err isa _StopOptimization || rethrow()
        nothing
    end
    if res === nothing
        return (; x = best_x, beams = apply(s.param, best_x), value = best_value,
                initial_value, iterations = nothing, evaluations = count[], stopped = true,
                result = nothing)
    end
    x = Optim.minimizer(res)
    return (; x, beams = apply(s.param, x), value = Optim.minimum(res), initial_value,
            iterations = Optim.iterations(res), evaluations = count[], stopped = false,
            result = res)
end

struct _StopOptimization <: Exception end

"""
    report(ev; lmax = 6) -> nothing

Print a power accounting and flux-moment summary of an evaluation.
"""
function report(io::IO, ev::Evaluation; lmax::Integer = 6)
    d = ev.deposition
    P = ForwardDiff.value(d.delivered)
    if ev.scenario.lpi !== nothing
        Ps = sum(ForwardDiff.value(b.power) for b in ev.beams)
        @printf(io, "Setup laser power: %.4g TW; after laser–plasma model:
", Ps)
        for (g, v) in sort(collect(lpi_summary(ev.beams, ev.effective)); by = first)
            @printf(io, "  %-20s %9.4g → %9.4g TW  (×%.3f)
", g, v.setup, v.effective, v.ratio)
        end
    end
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
