# This file defines the `Scenario`, which holds the precomputed deposition operator and flux
# response of a target and a parameterization, and the `Evaluation` of beams in it.

# ---------------------------------------------------------------------------------------
# Scenario: precomputed operators + evaluation
# ---------------------------------------------------------------------------------------

# The full flux response of a scenario, computed when first read.
mutable struct _LazyResponse
    const lock::ReentrantLock
    R::Union{Nothing,Matrix{Float64}}
end
_LazyResponse() = _LazyResponse(ReentrantLock(), nothing)

# What `optimize_pointing` runs on: a `Scenario`, or a `PhasedScenario` (phased.jl).
abstract type _AbstractScenario end

"""
    Scenario(target, pointing::GroupPointing; observe = nothing,
             radiation = observe !== nothing, time = 0.0, albedo = nothing, lpi = nothing,
             deposition_kw = (;), viewfactor_kw = (;), model = nothing)
    Scenario(target, beams; kwargs...)
    Scenario(target, parameterization; kwargs...)

A scenario contains everything precomputed for fast repeated evaluation of a pointing
design. Its parameterization (`scenario.param`) maps the parameter vectors of `evaluate`
and `optimize_pointing` to beams. The three forms differ in the parameterization:

- The first form builds the deposition operator for the template beams of `pointing`,
  which may also be a [`FramePointing`](@ref).
- The second form is `Scenario(target, GroupPointing(beams))`, with fields `:z`, `:r`
  and `:defocus` per group.
- The third form takes any other parameterization with methods of
  `apply(parameterization, x)` and `initial_parameters(parameterization)`, and builds
  the operator for the beams at its initial parameters.

The operator resolves the beam footprints where those beams aim (`reach`,
`DepositionOperator`). Besides the `DepositionOperator`, the scenario contains the view
factors and the linear flux response of the observed elements when `radiation` and
`observe` are given. `observe` names the components whose incident x-ray flux is the
figure of merit, e.g. `:capsule`. `lpi` is a laser–plasma model (`Backscatter`,
`CalibratedTransfer`, `CBET`, or a tuple of them), applied at `time` between the setup
beams and deposition.

The flux response is linear in the laser source at fixed albedos. By default these are
the albedos of the materials at `time`, with temperature-dependent ones at their
reference temperatures and thin foils intact. An explicit per-element `albedo` replaces
them, such as the converged `sol.albedo` of a solution at the design point (e.g. the
peak of a drive history, with windows burned through). For a time partway through a
shaped pulse (a phase), the material albedos at `time` differ from those of the wall.
The albedo of a `WallLossAlbedo` depends on the temperature history of each element,
and neither the reference temperature nor a standalone `solve` at `time` includes that
history. Pass the history albedos as `st.time, st.albedo`, with
`st = history_state(target, history, τ)` ([`history_state`](@ref)).

A scenario with a response (`observe` and a model) warns at construction when it is
built without `albedo` on a target whose materials include a [`WallLossAlbedo`](@ref),
or a [`BurnthroughAlbedo`](@ref) whose `before` is one. The loss is then taken at the
reference temperature (250 eV by default). Mid-pulse, this overstated the drive against
a drive history by up to 1.64 times on a public mock-up (research notes, "Drive
histories: sliver steps and the albedo of a phase"). `albedo = :reference` marks the
reference temperature as intended and builds the same scenario without the warning.

Objectives need only the moments of the observed flux, so the response is stored
contracted with the Legendre (`lmax = 12`) and real-harmonic (`lmax = 6`) weights. Each
evaluation then costs a 13 × n and a 49 × n product, and the setup costs a 62-column solve.
The full observed × n response `scenario.response`, and with it the `flux` of an evaluation,
is computed when first read. On a large (iterative) model, the `flux` of an evaluation is
one solve instead, until the response has been read.

Options of the operators pass through two named tuples:

- `deposition_kw` is passed to [`DepositionOperator`](@ref). Its options are:
  - `profile`, a [`VisRadProfile`](@ref) that selects the VisRad beam model instead of
    the f-cone superposition;
  - `ncone`, the number of f-cone directions (13);
  - `reach`, the aim motion in mm over which the operator stays accurate (2.0);
  - `resolution`, `level`, `maxlevel` and `shadow`, the quadrature refinement.
- `viewfactor_kw` is passed to [`RadiationModel`](@ref) and through it to
  [`view_factors`](@ref). Its options are:
  - `closed`, the closure correction, by default the declaration of the target that
    `hohlraum` makes;
  - `F`, given view factors, e.g. those of a target that differs only in its materials;
  - `storage`, one of `Float64`, `Float32` or `:hmatrix`;
  - `solver`, `visibility` (`:centroid` or `:sampled`), `visibility_grid`, `level` and
    `near`.

`model`, a `RadiationModel` of the target, is used instead of building one. It excludes
`viewfactor_kw`, and giving both throws an `ArgumentError`.

[`remake`](@ref)`(scenario; ...)` replaces the laser–plasma model, the parameterization,
the deposition options or the observed elements, time and albedos of the response, and
rebuilds only what depends on them.
"""
struct Scenario{P} <: _AbstractScenario
    target::Target
    param::P
    deposition::DepositionOperator
    model::Union{Nothing,RadiationModel}
    observed::Vector{Int}
    moment_weights::Matrix{Float64}   # legendre_weights of the observed elements (lmax = 12)
    ylm_weights::Matrix{Float64}      # ylm_weights of the observed elements (lmax = 6)
    moment_response::Matrix{Float64}  # moment_weights * response: moments per element power
    ylm_response::Matrix{Float64}     # ylm_weights * response
    lpi::Any                          # laser–plasma model (see apply_lpi)
    time::Float64
    _response::_LazyResponse          # the full response, once read (shared by subsets)
    albedo::Union{Nothing,Vector{Float64}}   # explicit albedos of the response
    _split::Vector{Vector{Symbol}}    # group sets whose element powers `evaluate` records
end

# Builds a scenario whose evaluations record the powers of no group sets. `_with_split`
# adds the sets.
Scenario(t::Target, param, op::DepositionOperator, model, obs, W, Y, WR, YR, lpi, time, response::_LazyResponse,
         albedo) = Scenario(t, param, op, model, obs, W, Y, WR, YR, lpi, time, response, albedo, Vector{Symbol}[])

const _CACHED_YLM_LMAX = 6

Base.getproperty(s::Scenario, name::Symbol) =
    name === :response ? _response(s) : getfield(s, name)
Base.propertynames(::Scenario) =
    (:target, :param, :deposition, :model, :observed, :response, :moment_weights,
     :ylm_weights, :moment_response, :ylm_response, :lpi, :time, :albedo)

# Whether the scenario has a radiation response (observed elements and a model).
_radiative(s::Scenario) = s.model !== nothing && !isempty(s.observed)

function _response(s::Scenario)
    _radiative(s) || return nothing
    c, m = getfield(s, :_response), s.model::RadiationModel
    lock(c.lock) do
        c.R === nothing && (c.R = flux_response(m, s.observed; time = s.time, s.albedo))
        c.R
    end
end

Scenario(t::Target, p::GroupPointing; kwargs...) = _scenario(t, p.template, p; kwargs...)
Scenario(t::Target, beams::AbstractVector{<:Beam}; kwargs...) =
    _scenario(t, beams, GroupPointing(beams); kwargs...)
Scenario(t::Target, p; kwargs...) = _scenario(t, apply(p, initial_parameters(p)), p; kwargs...)

function _scenario(t::Target, beams::AbstractVector{<:Beam}, parameterization;
                   observe = nothing, radiation::Bool = observe !== nothing, time = 0.0,
                   albedo = nothing, deposition_kw = (;), viewfactor_kw = (;), lpi = nothing,
                   model = nothing)
    model === nothing || model.target === t ||
        throw(ArgumentError("model is a radiation model of another target"))
    model === nothing || isempty(viewfactor_kw) ||
        throw(ArgumentError("give a model or viewfactor_kw, not both (the model's view factors are fixed)"))
    # `:reference` selects the material albedos at `time` and marks them as intended.
    meant = albedo === :reference
    albedo isa Symbol && !meant &&
        throw(ArgumentError("albedo takes albedos per element, `nothing` or `:reference`, got :$albedo"))
    meant && (albedo = nothing)
    op = DepositionOperator(t, beams; deposition_kw...)
    obs = _observed(t, observe)
    model = !radiation ? nothing : model === nothing ? RadiationModel(t; viewfactor_kw...) : model
    meant || albedo !== nothing || model === nothing || isempty(obs) || _warn_reference_albedo(t, time)
    W, Y, WR, YR = _moment_responses(t, model, obs, time, albedo)
    return Scenario(t, parameterization, op, model, obs, W, Y, WR, YR, lpi, Float64(time),
                    _LazyResponse(), albedo === nothing ? nothing : collect(Float64, albedo))
end

# Without albedos, the loss in a scenario response on a target whose wall loss depends on
# a temperature history is taken at the reference temperature of the model. Such a loss is a
# `WallLossAlbedo`, or the `BurnthroughAlbedo` of a thin foil over one. Mid-pulse, the
# reference loss overstates the drive against the history (research notes, "Drive
# histories: sliver steps and the albedo of a phase"). `_history_albedo` flags these
# albedos, and the scenario warns once per construction.
_history_albedo(a) = a isa WallLossAlbedo || a isa BurnthroughAlbedo && a.before isa WallLossAlbedo

function _warn_reference_albedo(t::Target, time)
    i = findfirst(c -> c.kind === :surface && _history_albedo(c.material.albedo), t.components)
    i === nothing && return nothing
    c = t.components[i]
    w = c.material.albedo isa BurnthroughAlbedo ? c.material.albedo.before : c.material.albedo
    @warn "Scenario: $(c.part)'s WallLossAlbedo is evaluated at its reference temperature " *
          "($(w.reference_temperature) eV) at time = $(Float64(time)) ns, as a wall held at it since the " *
          "pulse start; mid-pulse that can overstate the drive against a drive_history (by up to 1.64 times " *
          "on a public mock-up). For a phase of a pulse pass the history's albedos, " *
          "`st = history_state(target, history, time)` with `time = st.time, albedo = st.albedo` " *
          "(or `remake(scenario; history, time)`); `albedo = :reference` marks the reference " *
          "temperature is intended."
    return nothing
end

# Returns the observed elements of `observe`, a set of component or part names. `nothing`
# gives none.
_observed(t::Target, observe) =
    observe === nothing ? Int[] : elements(t, (observe isa Symbol ? (observe,) : observe)...)

# Returns the mode weights of the observed elements and, with a model, their responses.
# One flux-response solve covers all weights.
function _moment_responses(t::Target, model, obs, time, albedo)
    W = isempty(obs) ? zeros(0, 0) : legendre_weights(t, obs; lmax = _CACHED_LMAX)
    Y = isempty(obs) ? zeros(0, 0) : ylm_weights(t, obs; lmax = _CACHED_YLM_LMAX)
    WR, YR = zeros(0, 0), zeros(0, 0)
    if model !== nothing && !isempty(obs)
        WY = flux_response(model, obs; time, albedo, weights = vcat(W, Y))
        WR, YR = WY[1:size(W, 1), :], WY[size(W, 1)+1:end, :]
    end
    return W, Y, WR, YR
end

# Returns the beams that the deposition operator of a scenario is built for, from its
# parameterization.
_scenario_beams(p::GroupPointing) = p.template
_scenario_beams(p) = apply(p, initial_parameters(p))

# Marks a keyword of `remake(::Scenario)` that was not given.
struct _Unchanged end
const _UNCHANGED = _Unchanged()

"""
    remake(scenario; parameterization, lpi, deposition_kw, observe, time, albedo, model,
           history) -> Scenario

Copy of `scenario` with the given parts replaced. What depends on a replaced part is
rebuilt, and everything else is shared:

- `lpi` (a laser–plasma model or `nothing`): nothing is rebuilt.
- `parameterization`: the deposition operator is kept when the beams of the new
  parameterization enter through the ports of the operator in the same order with the
  same f-cone directions, and, under a [`VisRadProfile`](@ref), with the same envelope
  ratios. These beams are the template of a `GroupPointing` or `FramePointing`, or the
  beams of any other parameterization at its initial parameters. Otherwise the change is
  refused with an `ArgumentError`, unless `deposition_kw` is given too. A kept operator
  resolves the footprints where the original beams aimed (its `reach`), as for
  [`subset`](@ref).
- `deposition_kw`: the operator is rebuilt for the beams of the parameterization with
  these options. Options not given take their defaults, because a scenario does not
  record the options its operator was built with. Pass all options that differ from the
  defaults, such as `profile`, `reach` and `resolution`.
- `observe`, `time`, `albedo` and `model` (a `RadiationModel` of the same target, or
  `nothing` for none): the mode weights and moment responses are recomputed with one
  flux-response solve, and the full `response` is computed again when read.
  With `albedo = nothing` or `:reference`, the material albedos at `time` are used
  without a warning; any warning was issued when the scenario was constructed.

The target, `radiation` and `viewfactor_kw` do not change through `remake`, and giving
them throws an `ArgumentError`. They need a new `Scenario`, or a `model` built for the
purpose.

`history`, a [`drive_history`](@ref) of the scenario target run with
`keep_albedo = true`, sets with `time` the albedos and time of the history state at
`time`. `remake(scenario; history = h, time = τ)` is exactly
`st = `[`history_state`](@ref)`(scenario.target, h, τ)` followed by
`remake(scenario; st.time, st.albedo)`. Remade so, a phase of a design and its reference
at the same time have the same `time`, which the reference checks of
[`mode_gate`](@ref) and [`drive_floor`](@ref) compare.
"""
function remake(s::Scenario; parameterization = _UNCHANGED, lpi = _UNCHANGED,
                deposition_kw = _UNCHANGED, observe = _UNCHANGED, time = _UNCHANGED,
                albedo = _UNCHANGED, model = _UNCHANGED, history = nothing, kwargs...)
    isempty(kwargs) ||
        throw(ArgumentError("remake(::Scenario) cannot change $(join(keys(kwargs), ", ")): build a new Scenario " *
                            "(or pass a `model`); it changes parameterization, lpi, deposition_kw, observe, " *
                            "time, albedo and model"))
    if history !== nothing
        time === _UNCHANGED && throw(ArgumentError("a `history` needs the `time` (ns) of its state"))
        albedo === _UNCHANGED ||
            throw(ArgumentError("give the albedos by `history` or by `albedo`, not both"))
        st = history_state(s.target, history, time)
        time, albedo = st.time, st.albedo
    end
    t = s.target
    param = parameterization === _UNCHANGED ? s.param : parameterization
    op = s.deposition
    if deposition_kw !== _UNCHANGED
        op = DepositionOperator(t, _scenario_beams(param); deposition_kw...)
    elseif parameterization !== _UNCHANGED
        _check_operator(op, _scenario_beams(param))
    end
    lp = lpi === _UNCHANGED ? s.lpi : lpi
    if all(v -> v === _UNCHANGED, (observe, time, albedo, model))
        return Scenario(t, param, op, s.model, s.observed, s.moment_weights, s.ylm_weights,
                        s.moment_response, s.ylm_response, lp, s.time, getfield(s, :_response),
                        s.albedo)
    end
    m = model === _UNCHANGED ? s.model : model
    m === nothing || m.target === t ||
        throw(ArgumentError("model is a radiation model of another target"))
    obs = observe === _UNCHANGED ? s.observed : _observed(t, observe)
    tm = time === _UNCHANGED ? s.time : Float64(time)
    alb = albedo === _UNCHANGED ? s.albedo : albedo === nothing || albedo === :reference ? nothing :
          collect(Float64, albedo)
    W, Y, WR, YR = _moment_responses(t, m, obs, tm, alb)
    return Scenario(t, param, op, m, obs, W, Y, WR, YR, lp, tm, _LazyResponse(), alb)
end

# Checks that the operator `op` is valid for `beams`. Their ports and f-cone directions must
# match in order, and so must their envelope ratios under the VisRad profile. A mismatch
# throws an `ArgumentError`.
function _check_operator(op::DepositionOperator, beams)
    fail(why) = throw(ArgumentError("the parameterization's beams $why the deposition operator's: " *
                                    "give deposition_kw to rebuild it"))
    length(beams) == length(op.ports) || fail("are $(length(beams)), not the $(length(op.ports)) of")
    all(i -> beams[i].port == op.ports[i], eachindex(beams)) || fail("enter through other ports than")
    all(i -> cone_directions(beams[i], length(op.dirs[i]))[1] == op.dirs[i], eachindex(beams)) ||
        fail("have other f-cone directions (f-numbers) than")
    op.profile === nothing || [_vr_envelope(op.profile, b) for b in beams] == op.envelope ||
        fail("have other VisRad envelope ratios than")
    return nothing
end

"""
    subset(scenario, keep; parameterization = restricted to the kept beams)

Scenario for the beams `keep` (indices into the scenario beams), e.g. after dropping a
bundle. The subset reuses the view factors, flux response and mode weights, and only
the deposition operator is restricted. The template values (pointing, power) of the
kept beams come from the scenario parameterization.

By default, the new parameterization has the groups of the scenario. For a
[`GroupPointing`](@ref), the keys of the groups that still contain a beam are kept, so the
parameter names of the subset are those of the full pointing. The fields, mirroring,
`relative`, `unit` and sphere carry over. With `unit = :quad`, the default, the azimuth and
split of a quad are those of its kept beams. Pass a parameterization to group the kept beams
otherwise. For a [`FramePointing`](@ref) scenario, the frames, groups, fields and
`relative` of the kept beams are kept by default.

For a [`HeldParameters`](@ref) or [`PowerScales`](@ref) scenario, the default wraps the
subset of the inner parameterization, so the kept parameters keep their names. The
result equals the wrappers built by hand on the kept beams:

- For a parameter whose group still contains some of its beams, the value and the held
  or free state are kept. The value applies to the kept beams as it applied to all, as an
  absolute aim or a shift.
- A parameter whose group contains none of the kept beams is dropped, held or free.
- The power scale of a group with no kept beam is dropped too.

Any other parameterization needs an explicit `parameterization`.
"""
function subset(s::Scenario, keep::AbstractVector{<:Integer}; parameterization = nothing)
    param = parameterization === nothing ? _subset_param(s.param, keep) : parameterization
    return Scenario(s.target, param, s.deposition[keep], s.model, s.observed,
                    s.moment_weights, s.ylm_weights, s.moment_response, s.ylm_response,
                    s.lpi, s.time, getfield(s, :_response), s.albedo)
end

# Restricts `p` to the beams `keep`. The keys of the groups that still contain a beam are
# kept, and so are their parameter names. The fields, mirroring, `relative`, `unit` and
# sphere of `p` carry over.
function _subset_param(p::GroupPointing, keep)
    key = similar(p.keys, length(p.template))
    for (g, mem) in enumerate(p.members), bi in mem
        key[bi] = p.keys[g]
    end
    return _group_pointing(p.template[keep], key[keep]; fields = Tuple(p.fields), mirror = p.mirror,
                           relative = p.relative, unit = p.unit, radius = _spherical(p) ? p.radius : nothing,
                           center = p.center)
end
_subset_param(p, keep) =
    throw(ArgumentError("give a parameterization for a scenario that is not on a GroupPointing or FramePointing"))

_nparams_shown(p::GroupPointing) = nparams(p)
_nparams_shown(p) = "?"

function Base.show(io::IO, s::Scenario)
    print(io, "Scenario(", nelements(s.target), " elements, ", length(s.deposition.ports),
          " beams, ", _nparams_shown(s.param), " parameters, ",
          length(s.observed), " observed elements)")
end

"""
    Evaluation

Result of `evaluate`: the setup `beams`, the `effective` beams after the laser–plasma
model of the scenario, their `deposition`, and `flux` — incident x-ray flux on the
observed elements of the scenario (TW/mm²; `nothing` without radiation). Deposition
fractions are relative to effective (post-LPI) power. The flux is computed when first
read; `flux_moments` and `flux_ylm` do not need it.
"""
struct Evaluation{T,B,E}
    scenario::Scenario
    beams::Vector{B}
    effective::Vector{E}
    deposition::Deposition{T}
    _flux::Base.RefValue{Union{Nothing,Vector{T}}}   # observed flux, once computed
    _split::Union{Nothing,_GroupSplit{T}}            # powers of the scenario group sets
    _moments::Base.RefValue{Union{Nothing,Vector{T}}}  # cached Legendre moments, once computed
    _ylm::Base.RefValue{Union{Nothing,Vector{T}}}      # cached Y_lm coefficients, once computed
end

Evaluation(s::Scenario, beams::Vector, effective::Vector, d::Deposition{T},
           flux::Base.RefValue{Union{Nothing,Vector{T}}}, split = nothing) where {T} =
    Evaluation(s, beams, effective, d, flux, split, Ref{Union{Nothing,Vector{T}}}(nothing),
               Ref{Union{Nothing,Vector{T}}}(nothing))

# Returns the cached product of the moment (`moment_response`) or Y_lm (`ylm_response`)
# response of the scenario with the deposition. Each product is computed once per
# evaluation, and a gate, a drive floor, a ceiling and the objective on one phase all
# read the same copy.
function _cached_product(c::Base.RefValue{Union{Nothing,Vector{T}}}, R::Matrix{Float64},
                         P::Vector{T}) where {T}
    v = c[]
    v === nothing || return v
    w = _mul(R, P)
    c[] = w
    return w
end

Base.getproperty(ev::Evaluation, name::Symbol) =
    name === :flux ? _flux(ev) : getfield(ev, name)
Base.propertynames(::Evaluation) = (:scenario, :beams, :effective, :deposition, :flux)

function _flux(ev::Evaluation)
    s = getfield(ev, :scenario)
    _radiative(s) || return nothing
    c = getfield(ev, :_flux)
    c[] === nothing && (c[] = _observed_flux(s, getfield(ev, :deposition).element))
    return c[]
end

# Returns the observed flux for absorbed laser power P. The flux is the full response
# times P. On a large model whose response has not been read, it is one iterative solve
# instead, since the response would take a solve per observed element.
function _observed_flux(s::Scenario, P)
    m = s.model
    (m.iterative && getfield(s, :_response).R === nothing) || return _mul(s.response, P)
    t = m.target
    α = s.albedo === nothing ? albedos(t, s.time) : s.albedo
    η = _constant_conversion(t) ? conversions(t) : conversions(t, s.time)
    B = _iterative_radiosity(m, α, η .* P ./ t.elem_area, nothing)
    return _mul(m.F, B)[s.observed]
end

"""
    evaluate(scenario, x) / evaluate(scenario, beams) -> Evaluation
"""
evaluate(s::Scenario, x::AbstractVector{<:Real}) = evaluate(s, apply(s.param, x))

function evaluate(s::Scenario, beams::AbstractVector{<:Beam})
    lpi = _resolve_lpi(s.lpi, beams)
    return _evaluate(s, beams, apply_lpi(lpi, beams, s.time), lpi)
end

# Computes `evaluate` for the setup `beams` and their effective beams `eff` under `lpi`,
# the laser–plasma model of the scenario resolved for the setup beams (`_resolve_lpi`).
function _evaluate(s::Scenario, beams::AbstractVector{<:Beam}, eff, lpi)
    sets = getfield(s, :_split)
    (isempty(sets) || !_plain_deposit(lpi)) && return _evaluation(s, beams, eff, lpi_deposit(lpi, s.deposition, eff, s.time), nothing)
    # The element powers of the group sets that the constraints of the run read
    # (`_with_split`) are computed in the same pass.
    T = mapreduce(b -> typeof(b.power), promote_type, eff)
    split = _GroupSplit(sets, eff, nelements(s.target), T)
    return _evaluation(s, beams, eff, _deposit(s.deposition, eff, split), split)
end

_evaluation(s, beams, eff, dep::Deposition{T}, split) where {T} =
    Evaluation(s, beams, eff, dep, Ref{Union{Nothing,Vector{T}}}(nothing), split)

# Returns whether the effective beams of the laser–plasma model are deposited as
# `deposit` deposits them. This holds when no model deposits light off the operator
# (`_deposits`).
_plain_deposit(m) = !_deposits(m)
_plain_deposit(ms::Union{Tuple,AbstractVector}) = !any(_deposits, ms)

# Computes `A * x` for a Float64 matrix and a vector of duals through BLAS. The values and
# each partial are the columns of one matrix product; a generic dual mat-vec runs about
# 10 times slower. Other element types use the generic product.
_mul(A::AbstractMatrix, x::AbstractVector) = A * x
function _mul(A::Matrix{Float64}, x::Vector{ForwardDiff.Dual{Tag,Float64,N}}) where {Tag,N}
    X = Matrix{Float64}(undef, length(x), N + 1)
    for (i, d) in enumerate(x)
        X[i, 1] = ForwardDiff.value(d)
        p = ForwardDiff.partials(d)
        for j in 1:N
            X[i, j+1] = p[j]
        end
    end
    Y = A * X
    return [ForwardDiff.Dual{Tag}(Y[i, 1], ntuple(j -> Y[i, j+1], Val(N))) for i in axes(Y, 1)]
end
