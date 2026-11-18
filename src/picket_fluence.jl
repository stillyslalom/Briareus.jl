# The outer cones' picket fluence on the wall (energy per area), the quantity the
# published gold-bubble scalings are written in, and `BubbleAbsorption`'s picket ratio
# from it, and the picket irradiance as a constraint. Research notes, "The swing and the
# bubble without new data (2026-11-18)", §5 "The bubble's speed from the pointing's
# picket fluence"; literature in `reference/litreview/bubble_drivers.md`.

"""
    PicketFluence{T}

The picket fluence of a set of beams on a target's wall, from [`picket_fluence`](@ref):
- `fluence`: per element of the target (kJ/mm²; ×1e5 → J/cm²), zero outside `elements`;
- `elements`: the elements counted (those of the named components);
- `energy`: the picket energy deposited on them, ∫F dA (kJ);
- `mean`: the energy-weighted mean fluence ∫F² dA / ∫F dA (kJ/mm²), the summary that
  [`picket_ratio`](@ref) reads;
- `peak`: the largest element fluence (kJ/mm²), for information;
- `area`: the effective footprint area `energy / mean` = (∫F dA)² / ∫F² dA (mm²), over
  all counted elements (both hemispheres of a hohlraum: a symmetric pointing has the
  `mean` of either side and twice its `area`).
"""
struct PicketFluence{T<:Real}
    fluence::Vector{T}
    elements::Vector{Int}
    energy::T
    mean::T
    peak::T
    area::T
end

function Base.show(io::IO, f::PicketFluence)
    @printf(io, "PicketFluence(%.4g kJ on %d elements: mean %.4g kJ/mm² over %.4g mm², peak %.4g kJ/mm²)",
            ForwardDiff.value(f.energy), length(f.elements), ForwardDiff.value(f.mean),
            ForwardDiff.value(f.area), ForwardDiff.value(f.peak))
end

"""
    picket_fluence(op, beams; groups = (:outer44, :outer50), components = (:wall,),
                   duration = 1.0) -> PicketFluence
    picket_fluence(op, beams, pulses, window; by = :group, groups, components)

The picket fluence F (energy per area, kJ/mm²) of the beams in `groups` (the outer
cones) on the elements of `components` (part or component names, as
[`elements`](@ref)), and its summaries ([`PicketFluence`](@ref)). The first form takes
each beam's `power` as its picket power (TW) over a picket of `duration` ns, or, with
`duration = 1`, as its picket energy (kJ). The second integrates each beam's pulse
(`pulses` keyed as for [`with_pulses`](@ref), with `by`) over `window = (t0, t1)` ns,
e.g. the first ns after the outer turn-on (the convention of the research notes for the
basecamp's and N221204's pickets; use the same in a design and its reference). The
other beams are deposited at zero power. Deposition is linear in power, so F is the
deposition of those beams at their picket energies over each element's area: the
fluence on the wall, with the incidence and the overlap of the spots. Backscatter or
transfer in the picket can be taken off the beams first (Ralph et al. 2018 use the
picket energy net of outer backscatter).

**The summary.** Every published bubble scaling is written in the outer cones' picket
fluence as a nominal energy over a nominal area: Ralph et al., PoP 25, 082701 (2018),
the outer picket energy over the outer phase-plate spot area (speed ∝ I^0.36); Callahan
et al., PoP 25, 056305 (2018), the outer picket energy over the one-side outer spot area,
53 mm², the footprint of one side's outer beams together (size ∝ (E/A)^0.5); Izumi et
al., PoP 28, 022706 (2021), the picket fluence relative to the basecamp's (speed
∝ E^0.2), 27 % lower when the quads and cones are split at the same energy, i.e. with
the overlap of the spots on the wall counted. `mean` = ∫F² dA / ∫F dA is that quantity
for a computed map: the fluence averaged over the deposited energy, the energy over the
effective area (∫F dA)²/∫F² dA. One isolated spot gives its energy over its effective
area (a flat top's E/A exactly; for a super-Gaussian of order n, 2^(-2/n) of its peak
fluence, half for a Gaussian); k coincident spots give k times one, separated spots the
same as one; and as an integral it is smooth in the pointing and nearly independent of
the mesh. `peak` follows one element and the mesh, and a mean over the area above half
the peak inherits the peak's threshold; neither is what the papers' energy over area
measures.

**Mesh.** An element average lowers ∫F² dA by about h²/(6a²) for a Gaussian spot of 1/e
radius a on elements of size h, so `mean` converges as O(h²) and is low on a coarse
mesh: with the NIF outer spots (0.67 × 0.39 mm 1/e semi-axes) it is about 1 % low at
h = 0.2 mm and 0.3 % at 0.1 mm on NIF hohlraums, and its ratio between two hohlraums
0.35 % high at 0.2 mm (research notes). Compare a design with its reference on the same
mesh size, where most of the error cancels in the ratio; a
radiation mesh of 0.8 mm is too coarse for the map, and a separate target for the
deposition alone is cheap.

Generic in the number type of the beams (ForwardDiff duals through the pointing).
*Derived*: a link built on the published definitions and exponents, not a facility
criterion.
"""
function picket_fluence(op::DepositionOperator, beams::AbstractVector{<:Beam};
                        groups = (:outer44, :outer50), components = (:wall,),
                        duration::Real = 1.0)
    gs, es = _symbols(groups), _picket_elements(op.target, components)
    any(b -> b.group in gs, beams) || throw(ArgumentError("no beams in the groups $gs"))
    picket = [remake(b; power = b.group in gs ? b.power * duration : zero(b.power)) for b in beams]
    return _picket_fluence(deposit(op, picket).element, op.target, es)
end

function picket_fluence(op::DepositionOperator, beams::AbstractVector{<:Beam}, pulses, window;
                        by = :group, groups = (:outer44, :outer50), components = (:wall,))
    t0, t1 = window
    t0 < t1 || throw(ArgumentError("the picket window must run forward, got $window"))
    key, gs = _beam_key(by), _symbols(groups)
    picket = [remake(b; power = b.group in gs ? _window_energy(_beam_pulse(pulses, key, b), t0, t1) * one(b.power) :
                                zero(b.power)) for b in beams]
    return picket_fluence(op, picket; groups = gs, components)
end

_symbols(s::Symbol) = (s,)
_symbols(s::Tuple{Vararg{Symbol}}) = s
_symbols(s) = Tuple(Symbol.(collect(s)))

function _picket_elements(t::Target, components)
    es = elements(t, _symbols(components)...)
    isempty(es) && throw(ArgumentError("no elements in the components $components"))
    return es
end

# The map and its summaries from element energies E (kJ).
function _picket_fluence(E::AbstractVector{T}, t::Target, es::Vector{Int}) where {T}
    F = zeros(T, length(E))
    s1, s2, pk = zero(T), zero(T), zero(T)
    for e in es
        f = E[e] / t.elem_area[e]
        F[e] = f
        s1 += E[e]
        s2 += E[e] * f
        pk = max(pk, f)
    end
    s1 > 0 || throw(ArgumentError("the picket deposits no energy on the components"))
    mean = s2 / s1
    return PicketFluence{T}(F, es, s1, mean, pk, s1 / mean)
end

# ∫ s dt over [t0, t1] (kJ per beam): exact for the piecewise-linear pulse, zero outside
# its table.
function _window_energy(s::PulseShape, t0, t1)
    a, b = max(t0, s.t[1]), min(t1, s.t[end])
    a < b || return 0.0
    ks = findall(τ -> a < τ < b, s.t)
    return _trapz([a; s.t[ks]; b], [s(a); s.p[ks]; s(b)])
end

"""
    picket_ratio(design::PicketFluence, reference::PicketFluence; at = 1.0) -> E

[`BubbleAbsorption`](@ref)'s `picket_ratio` E for `design`: `at` times the ratio of the
two designs' energy-weighted mean picket fluences (`mean`, [`picket_fluence`](@ref)), so
that E follows the pointing (overlap, incidence, spot sizes) and the picket energy.

**Relative to what.** The bubble's tip speed (0.36 mm/ns, ∝ E^0.2, Izumi 2021) was
fitted at E = 1: the HDC basecamp of Izumi et al., PoP 28, 022706 (2021), a 5.75 mm Au
hohlraum with unsplit outer quads (research notes, "Hohlraum dynamics calibrated to
published data (2026-10-24)"). With `reference =` [`basecamp_picket_fluence`](@ref)`()`
and `at = 1`, E is relative to that calibration; its pointing is a proxy (the facility's
is unpublished) and its picket energy is read by eye, so the absolute E is uncertain
by more than its changes between designs. A re-pointing study is better referred to a
design whose E is already set: `picket_ratio(new, old; at = E_old)` moves E by the
pointing alone (`E_old` is 1.1 in `examples/n221204.jl`, from the picket energy and the
split quads). Compute both maps on the same mesh size.

`remake(bubble; picket_ratio = E)` then carries the design into the tip speed (∝ E^0.2)
and the collision time (``t_c`` − `turn_on` ∝ E^-0.16). Opt-in: no default reads it, and
the exponents and calibrated constants are the model's. Generic in the number type.
*Derived*, not a facility criterion: the exponents are published (Izumi 2021: 0.2 ± 0.03
on seven shots), their application to a re-pointing at fixed energy is not (no
published re-pointing study; Izumi 2021's split-quad shot is the nearest, research
notes).
"""
picket_ratio(design::PicketFluence, reference::PicketFluence; at::Real = 1.0) =
    at * design.mean / reference.mean

# Izumi 2021's HDC SymCap basecamp (`data/izumi2021_fig8.csv`'s header): length and LEH
# radius (mm), with `_BASECAMP`'s wall and capsule radii; the outer picket energy (kJ, all
# outer beams, first ns after the outer turn-on: Izumi 2021 Fig. 3a read by eye, research
# notes "The P2 swing: capsule radius, backscatter and the collision"); the outer aims
# (mm on the axis, upper hemisphere): the package's proxy pointing for a 5.75 mm hohlraum.
const _BASECAMP_PICKET = (length = 10.13, leh_radius = 1.685, energy = 26.0,
                          aim = (outer44 = 5.20, outer50 = 5.01))

"""
    basecamp_picket_fluence(; energy = 26.0, h = 0.2, capsule_h = 0.4) -> PicketFluence

The outer picket fluence of [`BubbleAbsorption`](@ref)'s calibration platform, its
`picket_ratio` E = 1: the HDC SymCap basecamp of Izumi et al., PoP 28, 022706 (2021), a
5.75 × 10.13 mm Au hohlraum with 3.37 mm entrance holes and a 0.908 mm capsule, mesh
size `h` (mm; `capsule_h` the capsule's). Its 128 outer beams ([`nif_beams`](@ref), the
NIF outer phase-plate spots) carry equal shares of `energy`, the outer picket energy
(kJ, all outer beams; ~26 kJ in the first ns after the outer turn-on, Izumi 2021 Fig. 3a
read by eye, *L*), with the quads unsplit and aimed on the axis at z = ±5.20 mm (44.5°)
and ±5.01 mm (50°): the package's proxy pointing for this hohlraum, the one the bubble
opacity was calibrated with (the facility's basecamp pointing is unpublished, *L*). The
reference of [`picket_ratio`](@ref) for an absolute E; compute the design's map on the
same `h`.
"""
function basecamp_picket_fluence(; energy::Real = _BASECAMP_PICKET.energy, h::Real = 0.2,
                                 capsule_h::Real = 0.4)
    t = hohlraum(; length = _BASECAMP_PICKET.length, radius = _BASECAMP.radius,
                 leh_radius = _BASECAMP_PICKET.leh_radius, capsule_radius = _BASECAMP.capsule,
                 h, capsule_h)
    outer = [b for b in nif_beams() if b.group in (:outer44, :outer50)]
    beams = [remake(b; aim = (0.0, 0.0, sign(b.port[3]) * getfield(_BASECAMP_PICKET.aim, b.group)),
                    power = energy / length(outer)) for b in outer]
    return picket_fluence(DepositionOperator(t, beams), beams)
end

# ---------------------------------------------------------------------------------------
# Picket powers per beam
# ---------------------------------------------------------------------------------------

# A beam's picket power (TW): the mean of its pulse over the window [t0, t1] (the
# window's energy over its length), the pulse keyed as for `with_pulses` (`by` a symbol
# or a function of the beam: a concrete type either way, unlike `_beam_key`'s closures).
struct _WindowPower{P,K}
    pulses::P
    by::K
    t0::Float64
    t1::Float64
end
(w::_WindowPower)(b::Beam) =
    _window_energy(_beam_pulse(w.pulses, Base.Fix1(_pulse_key, w.by), b), w.t0, w.t1) / (w.t1 - w.t0)

_pulse_key(by::Symbol, b::Beam) = by === :group ? b.group : by === :quad ? Symbol(b.quad) :
                                  by === :beam ? Symbol(b.name) :
                                  throw(ArgumentError("by is :group, :quad, :beam or a function, got $(repr(by))"))
_pulse_key(by, b::Beam) = by(b)

# A beam's own power as its picket power.
struct _OwnPower end
(::_OwnPower)(b::Beam) = b.power

# The picket power per beam from the keywords: `pulses` over `window`, or `power` (a
# function of the beam), or (`own`) each beam's own power.
function _picket_power(pulses, window, by, power, own::Bool)
    if pulses !== nothing
        power === nothing || throw(ArgumentError("give the picket by `pulses` and `window` or by `power`, not both"))
        window === nothing && throw(ArgumentError("`pulses` need a picket `window = (t0, t1)` (ns)"))
        t0, t1 = window
        t0 < t1 || throw(ArgumentError("the picket window must run forward, got $window"))
        return _WindowPower(pulses, by, Float64(t0), Float64(t1))
    end
    window === nothing || throw(ArgumentError("a picket `window` needs `pulses`"))
    power === nothing || return power
    own && return _OwnPower()
    throw(ArgumentError("give the picket: `pulses` and `window`, or `power` (a function of the beam)"))
end

# The beams at their picket powers (TW): those in the groups `gs`, the others at zero.
_picket_beams(power, gs, beams) =
    [remake(b; power = b.group in gs ? oftype(b.power, power(b)) : zero(b.power)) for b in beams]

# The energy-weighted mean ∫I² dA / ∫I dA of element powers P over the elements `es`
# (the arithmetic of `_picket_fluence`).
function _energy_weighted(P::AbstractVector{T}, t::Target, es) where {T}
    s1, s2 = zero(T), zero(T)
    for e in es
        s1 += P[e]
        s2 += P[e] * (P[e] / t.elem_area[e])
    end
    s1 > 0 || throw(ArgumentError("the picket deposits no power on the parts"))
    return s2 / s1
end

# ---------------------------------------------------------------------------------------
# The picket irradiance on named parts, and a ceiling on it
# ---------------------------------------------------------------------------------------

"""
    picket_irradiance(ev, names...; pulses = nothing, window = nothing, by = :group,
                      power = nothing, groups = (:outer44, :outer50), measure = :peak) -> Real

The laser irradiance (TW/mm²) during the picket on the elements of the named parts or
components, from the beams of an [`Evaluation`](@ref) at their picket powers: the
largest over the elements (`measure = :peak`, the maximum irradiance at the picket), or
the energy-weighted mean ∫I² dA / ∫I dA (`measure = :mean`, the summary that
[`picket_ratio`](@ref) and the gold bubble's published scalings read; see
[`picket_fluence`](@ref)). The picket fluence is this irradiance times the picket's
length, so its ratios between designs are the fluence's.

**The picket powers.** With `pulses` (keyed as for [`with_pulses`](@ref), with `by`)
and `window = (t0, t1)` (ns), each beam's mean power over the window; with `power`, a
function of the beam that returns its picket power (TW), e.g. by cone from a picket
table; with neither, each beam's own power (the scenario's beams carry their picket
powers, or `power = b -> r * b.power` scales them by a picket-to-peak ratio `r`). The
beams counted are those in `groups`, the outer cones by default: every published
bubble scaling is written in the outer cones' picket fluence (Ralph et al., PoP 25,
082701 (2018); Callahan et al., PoP 25, 056305 (2018); Izumi et al., PoP 28, 022706
(2021); research notes, "The swing and the bubble without new data (2026-11-18)" §4);
the other beams are deposited at zero power.

The picket is deposited by the scenario's deposition operator with the evaluation's
setup beams (`ev.beams`: their pointing, spots and focus), without the scenario's
laser–plasma model, which is set for its own `time`: take backscatter or transfer in
the picket off the picket powers. One deposition per call. Generic in the number type
(ForwardDiff duals through the pointing). The peak follows one element and the mesh
(on a coarser mesh the hottest spot is averaged over more area); the mean is an
integral, nearly independent of the mesh ([`picket_fluence`](@ref)).
"""
function picket_irradiance(ev::Evaluation, names::Symbol...; pulses = nothing, window = nothing,
                           by = :group, power = nothing, groups = (:outer44, :outer50),
                           measure::Symbol = :peak)
    t = ev.scenario.target
    es = _named_elements(t, names)
    P = _picket_map(ev, _picket_power(pulses, window, by, power, true), _symbols(groups))
    return _picket_measure(P, t, es, measure)
end

# The element powers (TW) of the evaluation's setup beams at their picket powers.
_picket_map(ev::Evaluation, power, gs) = deposit(ev.scenario.deposition, _picket_beams(power, gs, ev.beams)).element

function _picket_measure(P, t::Target, es, measure::Symbol)
    measure === :peak && return _peak(P, t, es)
    measure === :mean && return _energy_weighted(P, t, es)
    throw(ArgumentError("measure is :peak or :mean, got $(repr(measure))"))
end

"""
    picket_ceiling(reference, names...; pulses = nothing, window = nothing, by = :group,
                   power = nothing, groups = (:outer44, :outer50), measure = :peak,
                   limit = picket_irradiance(reference, names...; ...), p = 64,
                   scale = limit, exact = true) -> Constraint

A ceiling on the laser irradiance during the picket over the elements of the named
parts or components ([`picket_irradiance`](@ref)), at or below the `reference`
evaluation's by default (a number in TW/mm² sets another limit): a [`Constraint`](@ref)
for [`optimize_pointing`](@ref). The picket powers (`pulses` and `window`, or `power`)
and the beams that count (`groups`) are those of [`picket_irradiance`](@ref), the same
for the reference and the design.

- `measure = :peak` (named `:picket_peak`): the maximum irradiance over the elements
  during the picket. The optimizer works on the smooth area-weighted `p`-norm of the
  picket irradiance scaled to the reference's peak (as [`irradiance_ceiling`](@ref)
  does at the scenario's time); with `exact = true` (the default) the run holds the
  peak itself, lowering the stand-in's working limit while the peak exceeds the limit,
  and reports the peak against the limit ([`check_constraints`](@ref),
  `PointingResult.constraints`).
- `measure = :mean` (named `:picket_mean`): the energy-weighted mean irradiance, the
  quantity that [`picket_ratio`](@ref) carries into [`BubbleAbsorption`](@ref)'s tip
  speed (∝ ratio^0.2): smooth, held directly.

The reference evaluation must be on the same target as the optimized scenario (a
[`subset`](@ref) of it shares the target), or evaluating the constraint throws an
`ArgumentError`; its laser–plasma model and time do not enter (the picket is deposited
without them). One deposition per evaluation, in the evaluation's number type.

*Derived*, not a facility criterion: the picket fluence is the measure of the published
bubble scalings, and holding its peak at the reference's is the re-pointing study's
choice (user, 2026-11-18), an engineering limit on the picket's hottest element; the
model's prediction of the bubble a design grows is the mean's ratio to the 0.2 power
([`picket_ratio`](@ref)), confidence low (research
notes, "The swing and the bubble without new data (2026-11-18)" §5).
"""
function picket_ceiling(reference::Evaluation, names::Symbol...; pulses = nothing, window = nothing,
                        by = :group, power = nothing, groups = (:outer44, :outer50),
                        measure::Symbol = :peak, limit::Union{Nothing,Real} = nothing, p::Real = 64,
                        scale::Union{Nothing,Real} = nothing, exact::Bool = true)
    measure in (:peak, :mean) || throw(ArgumentError("measure is :peak or :mean, got $(repr(measure))"))
    t = reference.scenario.target
    es = _named_elements(t, names)
    pw = _picket_power(pulses, window, by, power, true)
    gs = collect(Symbol, _symbols(groups))
    P = Float64.(ForwardDiff.value.(_picket_map(reference, pw, gs)))
    sum(P[e] for e in es) > 0 ||
        throw(ArgumentError("the reference's picket deposits no power on $(join(names, ", "))"))
    ref = _picket_measure(P, t, es, measure)
    factor = measure === :peak ? ref / _irradiance_norm(P, t, es, Float64(p)) : 1.0
    lim = limit === nothing ? ref : Float64(limit)
    lim > 0 || throw(ArgumentError("the picket irradiance limit must be positive"))
    value = _PicketValue(t, es, gs, pw, measure, Float64(p), factor, exact)
    return Constraint(measure === :peak ? :picket_peak : :picket_mean, value; upper = lim,
                      scale = something(scale, lim))
end

# The picket ceiling's value: the picket irradiance's p-norm on the reference's peak
# scale (`:peak`), or its energy-weighted mean (`:mean`); with `exact`, the reported
# value of `:peak` is the peak itself (`_exact_value`).
struct _PicketValue{F}
    target::Target
    elements::Vector{Int}
    groups::Vector{Symbol}
    power::F
    measure::Symbol
    p::Float64
    factor::Float64
    exact::Bool
end

function _check_picket(c::_PicketValue, ev::Evaluation)
    ev.scenario.target === c.target ||
        throw(ArgumentError("picket_ceiling: the evaluation is on another target than the reference"))
    return nothing
end

function (c::_PicketValue)(ev::Evaluation)
    _check_picket(c, ev)
    P = _picket_map(ev, c.power, c.groups)
    return c.measure === :mean ? _energy_weighted(P, c.target, c.elements) :
           c.factor * _irradiance_norm(P, c.target, c.elements, c.p)
end

function _exact_value(c::_PicketValue, ev::Evaluation)
    c.exact && c.measure === :peak || return nothing
    _check_picket(c, ev)
    return _peak(Float64.(ForwardDiff.value.(_picket_map(ev, c.power, c.groups))), c.target, c.elements)
end
