# This file computes the picket fluence of the outer cones on the wall (energy per area),
# the quantity in which the published gold-bubble scalings are written, and the picket
# ratio of `BubbleAbsorption` from it. It also provides the picket irradiance as a
# constraint, the bubble that follows the picket of the beams, and the prediction for a
# re-pointed design. See the research notes, "The swing and the bubble without new data
# (2026-11-18)", §5 "The bubble's speed from the pointing's picket fluence"; the
# literature is in `reference/litreview/bubble_drivers.md`.

"""
    PicketFluence{T}

The picket fluence of a set of beams on the wall of a target, from
[`picket_fluence`](@ref):
- `fluence`: per element of the target (kJ/mm²; ×1e5 → J/cm²), zero outside `elements`;
- `elements`: the elements counted (those of the named components);
- `energy`: the picket energy deposited on them, ∫F dA (kJ). Over a fixed spot area it
  is the spatially averaged fluence of the published bubble scalings, whose ratio
  between designs is the energy ratio (the `:average` of [`picket_irradiance`](@ref));
- `mean`: the energy-weighted mean fluence ∫F² dA / ∫F dA (kJ/mm²), the summary that
  [`picket_ratio`](@ref) uses;
- `peak`: the largest element fluence (kJ/mm²), for information;
- `area`: the effective footprint area `energy / mean` = (∫F dA)² / ∫F² dA (mm²), over
  all counted elements. These include both hemispheres of a hohlraum, so a symmetric
  pointing has the `mean` of either side and twice its `area`.
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
cones) on the elements of `components`, and its summaries ([`PicketFluence`](@ref)).
`components` are part or component names, as for [`elements`](@ref). The first form
takes the `power` of each beam as its picket power (TW) over a picket of `duration` ns,
or, with `duration = 1`, as its picket energy (kJ). The second form integrates the pulse
of each beam (`pulses` keyed as for [`with_pulses`](@ref), with `by`) over
`window = (t0, t1)` ns. An example window is the first ns after the outer turn-on, the
convention of the research notes for the pickets of the basecamp and N221204; the same
window should be used in a design and its reference. The other beams are deposited at
zero power. Deposition is linear in power, so F is the deposition of those beams at
their picket energies over the area of each element. It is the fluence on the wall,
with the incidence and the overlap of the spots. Backscatter or transfer in the picket
can be taken off the beams first; Ralph et al. 2018 use the picket energy net of outer
backscatter.

**The summaries.** Every published bubble scaling is written in the picket fluence of the
outer cones, as a nominal energy over a nominal area:

- Ralph et al., PoP 25, 082701 (2018) use the outer picket energy net of backscatter
  over the outer phase-plate spot area (speed ∝ I^0.36).
- Callahan et al., PoP 25, 056305 (2018) use the outer picket energy over the one-side
  outer spot area, 53 mm² (size ∝ (E/A)^0.5).
- Izumi et al., PoP 28, 022706 (2021) use the spatially averaged outer picket fluence
  over the first ns (tip position ∝ E^0.2). It is 27% lower when the quads and cones are
  split at the same energy.

Over a fixed area that quantity is `energy`, whose ratios between designs are the energy
ratios (the `:average` of [`picket_irradiance`](@ref)). `mean` = ∫F² dA / ∫F dA is the
energy over the effective area (∫F dA)²/∫F² dA, and it also counts the overlap of the
spots. One isolated spot gives its energy over its effective area. For a flat top that
is exactly E/A; for a super-Gaussian of order n it is 2^(-2/n) of the peak fluence, and
half the peak for a Gaussian. k coincident spots give k times the value of one, and
separated spots the same value as one. `mean` therefore includes the split-quad change
of Izumi 2021, which a fixed area does not. As an integral, it is smooth in the pointing
and nearly independent of the mesh. The two measures differ when the picket energies of
the beams are uneven or their spots move onto one another (research notes, "The
published picket fluence, the trough and the inner-beam corridor (2026-11-18)").
`peak` follows one element and the mesh, and it measures neither quantity.

**Mesh.** An element average lowers ∫F² dA by about h²/(6a²) for a Gaussian spot of 1/e
radius a on elements of size h. `mean` therefore converges as O(h²) and is low on a
coarse mesh. With the NIF outer spots (0.67 × 0.39 mm 1/e semi-axes) on NIF hohlraums,
it is about 1% low at h = 0.2 mm and 0.3% low at 0.1 mm, and its ratio between two
hohlraums is 0.35% high at 0.2 mm (research notes). A design and its reference should
be compared on the same mesh size, where most of the error cancels in the ratio. A
radiation mesh of 0.8 mm is too coarse for the map, and a separate target used only for
the deposition is cheap.

Generic in the number type of the beams (ForwardDiff duals through the pointing).
*Derived*: a link built on the published definitions and exponents. No facility
criterion sets it.
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

# Returns ∫ s dt over [t0, t1] (kJ per beam). The integral is exact for the
# piecewise-linear pulse and zero outside its table.
function _window_energy(s::PulseShape, t0, t1)
    a, b = max(t0, s.t[1]), min(t1, s.t[end])
    a < b || return 0.0
    ks = findall(τ -> a < τ < b, s.t)
    return _trapz([a; s.t[ks]; b], [s(a); s.p[ks]; s(b)])
end

"""
    picket_ratio(design::PicketFluence, reference::PicketFluence; at = 1.0) -> E

The `picket_ratio` E of [`BubbleAbsorption`](@ref) for `design`. E is `at` times the
ratio of the energy-weighted mean picket fluences of the two designs (`mean`,
[`picket_fluence`](@ref)), so that E follows the pointing (overlap, incidence, spot
sizes) and the picket energy.

**Relative to what.** The tip speed of the bubble (0.36 mm/ns, ∝ E^0.2, Izumi 2021) was
fitted at E = 1, on the HDC basecamp of Izumi et al., PoP 28, 022706 (2021). That
platform is a 5.75 mm Au hohlraum with unsplit outer quads (research notes, "Hohlraum
dynamics calibrated to published data (2026-10-24)"). With
`reference =` [`basecamp_picket_fluence`](@ref)`()` and `at = 1`, E is relative to that
calibration. Its pointing is a proxy, since the facility pointing is unpublished, and
its picket energy was read from a figure by eye. The absolute E is therefore uncertain
by more than its changes between designs. A re-pointing study is better referred to a
design whose E is already set: `picket_ratio(new, old; at = E_old)` changes E through
the pointing only. `E_old` is 1.1 in `examples/n221204.jl`, from the picket energy and
the split quads. Both maps should be computed on the same mesh size.

With `remake(bubble; picket_ratio = E)`, the design then enters the tip speed (∝ E^0.2)
and the collision time (``t_c`` − `turn_on` ∝ E^-0.16). The link is opt-in, and no
default uses it; the exponents and calibrated constants are those of the model. Generic
in the number type. *Derived*; no facility criterion sets it. The exponents are
published (Izumi 2021: 0.2 ± 0.03 on seven shots), but their application to a
re-pointing at fixed energy is not. No re-pointing study is published, and the nearest
test is the split-quad shot of Izumi 2021 (research notes).
"""
picket_ratio(design::PicketFluence, reference::PicketFluence; at::Real = 1.0) =
    at * (design.mean / reference.mean)     # `at` itself, exactly, at the reference

# The HDC SymCap basecamp of Izumi 2021, from the header of `data/izumi2021_fig8.csv`.
# The fields are:
# - the length and LEH radius (mm), with the wall and capsule radii of `_BASECAMP`;
# - the outer picket energy (kJ, all outer beams, first ns after the outer turn-on), read
#   by eye from Izumi 2021 Fig. 3a (research notes, "The P2 swing: capsule radius,
#   backscatter and the collision");
# - the outer aims (mm on the axis, upper hemisphere), the proxy pointing of the package
#   for a 5.75 mm hohlraum.
const _BASECAMP_PICKET = (length = 10.13, leh_radius = 1.685, energy = 26.0,
                          aim = (outer44 = 5.20, outer50 = 5.01))

"""
    basecamp_picket_fluence(; energy = 26.0, h = 0.2, capsule_h = 0.4) -> PicketFluence

The outer picket fluence of the calibration platform of [`BubbleAbsorption`](@ref),
where the `picket_ratio` is E = 1. The platform is the HDC SymCap basecamp of Izumi et
al., PoP 28, 022706 (2021), a 5.75 × 10.13 mm Au hohlraum with 3.37 mm entrance holes and
a 0.908 mm capsule. The mesh size is `h` (mm), and `capsule_h` is that of the capsule.
The 128 outer beams ([`nif_beams`](@ref), with the NIF outer phase-plate spots) carry
equal shares of `energy`, the outer picket energy in kJ over all outer beams. The value
is about 26 kJ in the first ns after the outer turn-on, read by eye from Izumi 2021
Fig. 3a (*L*). The quads are unsplit and aimed on the axis at z = ±5.20 mm (44.5°) and
±5.01 mm (50°). This is the proxy pointing of the package for this hohlraum, the one
with which the bubble opacity was calibrated; the basecamp pointing of the facility is
unpublished (*L*). The result is the reference of [`picket_ratio`](@ref) for an absolute
E, and the map of the design should be computed on the same `h`.
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

# The picket power of a beam (TW), the mean of its pulse over the window [t0, t1], which
# is the energy in the window over its length. The pulse is keyed as for `with_pulses`.
# `by` is a symbol or a function of the beam, and either way it has a concrete type,
# unlike the closures of `_beam_key`.
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

# The power of each beam as its picket power.
struct _OwnPower end
(::_OwnPower)(b::Beam) = b.power

# Returns the picket power per beam from the keywords. It is `pulses` over `window`, or
# `power` (a function of the beam), or with `own` the power of each beam.
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

# Returns the beams at their picket powers (TW) for the groups `gs`, and the other beams
# at zero power.
_picket_beams(power, gs, beams) =
    [remake(b; power = b.group in gs ? oftype(b.power, power(b)) : zero(b.power)) for b in beams]

# Returns the energy-weighted mean ∫I² dA / ∫I dA of element powers P over the elements
# `es`, with the arithmetic of `_picket_fluence`.
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

# The measures of a picket map and their options, the keywords of `picket_irradiance`.
struct _PicketMeasure
    measure::Symbol
    area::Float64           # `:average`: the fixed spot area (mm²)
    corridor::Float64       # `:corridor`: the distance from an inner centerline (mm) …
    edge::Float64           # … and the width of the soft edge of the weight (mm)
    inner::Vector{Symbol}   # the groups whose centerlines the corridor follows
end

function _picket_measure_options(measure::Symbol, area::Real, corridor::Real, edge::Real, inner)
    measure in (:peak, :weighted, :average, :corridor) ||
        throw(ArgumentError("measure is :peak, :weighted, :average or :corridor, got $(repr(measure))"))
    area > 0 || throw(ArgumentError("the spot area must be positive, got $area"))
    0 < edge < 2corridor ||
        throw(ArgumentError("the corridor's edge must be positive and less than twice its distance, got $edge and $corridor"))
    return _PicketMeasure(measure, Float64(area), Float64(corridor), Float64(edge), collect(Symbol, _symbols(inner)))
end

# With `lpi = true`, the effective beams of the evaluation are used at their present
# powers.
_check_lpi(lpi::Bool, pulses, power) =
    lpi && (pulses !== nothing || power !== nothing) &&
    throw(ArgumentError("lpi = true deposits the evaluation's effective beams at their own powers: " *
                        "give no `pulses` or `power`"))

"""
    picket_irradiance(ev, names...; pulses = nothing, window = nothing, by = :group,
                      power = nothing, groups = (:outer44, :outer50), lpi = false,
                      measure = :peak, area = 53.0, corridor = 1.3, edge = 0.3,
                      inner = (:inner23, :inner30)) -> Real

The laser irradiance (TW/mm²) of the beams of an [`Evaluation`](@ref) in `groups` (the
outer cones by default) on the elements of the named parts or components. `measure`
selects the summary:

- `:peak`: the largest over the elements;
- `:weighted`: the energy-weighted mean ∫I² dA / ∫I dA, the summary that
  [`picket_ratio`](@ref) and [`PicketLinkedBubble`](@ref) use ([`picket_fluence`](@ref));
- `:average`: the spatially averaged irradiance over a fixed spot area, which is the
  power absorbed on the elements over `area` (mm²) and the quantity of the published
  bubble scalings (below);
- `:corridor`: the mean irradiance on the elements near the centerlines of the beams in
  `inner`, ∫w I dA / ∫w dA with a weight w that is 1 within `corridor − edge/2` of a
  centerline and 0 beyond `corridor + edge/2` (below).

The irradiance times the length of the phase is its fluence, so the ratios between
designs are the fluence ratios. The powers are those of whatever phase the evaluation
represents: the picket, the trough or the foot.

**The powers.** With `pulses` (keyed as for [`with_pulses`](@ref), with `by`) and
`window = (t0, t1)` (ns), the power of each beam is its mean power over the window. With
`power`, a function of the beam that returns its power (TW), the power can be set, for
example, by cone from a picket table. With neither, the power of each beam is used, as
in an evaluation run at the powers of a phase, or with `power = b -> r * b.power` to
scale them. The other beams are deposited at zero power, by the deposition operator of
the scenario with the setup beams of the evaluation (`ev.beams`, which set the pointing,
spots and focus). With `lpi = true` the beams are instead the effective beams of the
evaluation (`ev.effective`), after the laser–plasma model of the scenario at its `time`,
for example net of a [`Backscatter`](@ref) and after a transfer. They are used at their
present powers, and `pulses` and `power` are then an error. Otherwise the laser–plasma
model is not used, and backscatter or transfer is taken off the powers given. Each call
performs one deposition. On an evaluation inside [`optimize_pointing`](@ref) that
already contains the element powers of the groups at their present powers
([`picket_ceiling`](@ref)), no deposition is needed. Generic in the number type
(ForwardDiff duals through the pointing and powers).

**The published quantity** (`:average`). The published scalings use these definitions:

- Callahan et al., PoP 25, 056305 (2018) divide the picket energy of the outer cones by
  the one-side outer spot area of the NIF phase plates, 53 mm². That area is the default
  `area`, so that the irradiance times the picket length in ns is their E/A in kJ/mm².
- Ralph et al., PoP 25, 082701 (2018) take the outer picket energy net of backscatter
  over the outer spot area.
- Izumi et al., PoP 28, 022706 (2021) take the spatially averaged outer picket fluence
  over the first ns.

In Briareus terms, the quantity is the power *absorbed on the named parts*, over a fixed
area. For the wall, light clipped on the outside of the entrance holes, light landing on
the capsule and light leaving the target do not grow the bubble on the wall (*L*). The
power is *net of backscatter* when the powers are, either with `lpi = true` and a
backscatter model or with powers net of backscatter. A ratio to a reference does not
depend on the area, and an absolute value uses the convention of Callahan. A fixed area
does not count the overlap of the spots, which `:weighted` does. The two measures
differ for uneven powers between beams and for spots moved onto one another (research
notes, "The published picket fluence, the trough and the inner-beam corridor
(2026-11-18)").

**The inner-beam corridor** (`:corridor`). Izumi et al., PoP 30, 032706 (2023) found
that outer bubble spikes in line with the 23.5° beams clip those beams and make the core
more oblate. The corridor measures the outer light where the inner beams pass. A
centerline is the line through the aim of a beam along its direction, for each beam of
`inner` with nonzero power in the evaluation. The distance is that from the element
centroid, in 3-D, to the nearest centerline. The weight is one minus the product over
the centerlines of (1 − s), where s rises from 0 to 1 across the edge as a quintic
smoothstep. The weight is therefore smooth (C²) in the aims of the inner beams as they
move, and exactly 0 or 1 away from the edge; elements beyond every edge are skipped.
`corridor = 1.3` mm and `edge = 0.3` mm are engineering defaults. The value 1.3 mm was
chosen for a re-pointing study (user, 2026-11-18). The edge is wider than an element of a
0.2 mm deposition map, so that the mean changes smoothly as centroids cross it. Neither
value is a published criterion.

The peak follows one element and the mesh, since on a coarser mesh the hottest spot is
averaged over more area. The other three measures are integrals, nearly independent of
the mesh ([`picket_fluence`](@ref)). [`picket_hot_elements`](@ref) lists the hottest
elements with their distance to the nearest inner centerline.
"""
function picket_irradiance(ev::Evaluation, names::Symbol...; pulses = nothing, window = nothing,
                           by = :group, power = nothing, groups = (:outer44, :outer50), lpi::Bool = false,
                           measure::Symbol = :peak, area::Real = 53.0, corridor::Real = 1.3, edge::Real = 0.3,
                           inner = (:inner23, :inner30))
    t = ev.scenario.target
    es = _named_elements(t, names)
    mo = _picket_measure_options(measure, area, corridor, edge, inner)
    _check_lpi(lpi, pulses, power)
    P = _picket_map(ev, _picket_power(pulses, window, by, power, true), _symbols(groups), lpi)
    return _picket_measure(P, t, es, mo, ev.beams)
end

# Returns the element powers (TW) of the setup beams of the evaluation, or with `lpi` of
# its effective beams, at their picket powers. At their present powers, the deposition of
# the groups kept in the evaluation is used where available (`_split_map`); otherwise the
# beams are deposited.
function _picket_map(ev::Evaluation, power, gs, lpi::Bool)
    P = _split_map(ev, power, gs, lpi)
    P === nothing || return P
    return lpi ? _picket_map(ev.scenario.deposition, power, gs, ev.effective) :
                 _picket_map(ev.scenario.deposition, power, gs, ev.beams)
end
_picket_map(op::DepositionOperator, power, gs, beams) = deposit(op, _picket_beams(power, gs, beams)).element

# Returns the element powers of the groups `gs` from the deposition in the evaluation, or
# `nothing`. The evaluation keeps them through the split of the `Scenario`, which the
# optimizer sets for the constraints that use them. They are valid for the beams at their
# present powers when the evaluation deposited those beams as `deposit` does. Those are
# its effective beams (`lpi`), or its setup beams when the laser–plasma model left them
# unchanged. The sums are identical to the deposition in `_picket_map` (`_GroupSplit`),
# and equal within rounding when phases of a `PhasedScenario` were deposited in one pass
# (`_deposit_sets`).
function _split_map(ev::Evaluation, power, gs, lpi::Bool)
    sp = getfield(ev, :_split)
    (sp === nothing || !(power isa _OwnPower)) && return nothing
    lpi || ev.effective === ev.beams || return nothing
    j = findfirst(==(_group_set(gs)), sp.sets)
    return j === nothing ? nothing : sp.element[j]
end

function _picket_measure(P, t::Target, es, m::_PicketMeasure, beams)
    m.measure === :peak && return _peak(P, t, es)
    m.measure === :weighted && return _energy_weighted(P, t, es)
    m.measure === :average && return _absorbed(P, es) / m.area
    return _corridor_mean(P, t, es, _centerlines(beams, m.inner), m)
end

function _absorbed(P::AbstractVector{T}, es) where {T}
    s = zero(T)
    for e in es
        s += P[e]
    end
    s > 0 || throw(ArgumentError("the picket deposits no power on the parts"))
    return s
end

# The centerlines (aim, unit direction) of the beams in `inner` with nonzero power.
function _centerlines(beams::AbstractVector{Beam{S}}, inner) where {S<:Real}
    lines = Tuple{SVector{3,S},Vec3}[]
    for b in beams
        (b.group in inner && ForwardDiff.value(b.power) != 0) || continue
        push!(lines, (b.aim, normalize(direction(b))))
    end
    return lines
end
_centerlines(beams::AbstractVector{<:Beam}, inner) = _centerlines(_common_type(beams), inner)

# Returns the weight of the corridor at x, 1 − Π (1 − s(u)) over the centerlines, with s
# the quintic smoothstep of u = (hi − d)/edge for a distance d between lo = hi − edge and
# hi. Each distance is first classified by its value. Beyond hi the factor is 1, and
# within lo the weight is 1. Both are exact because s is flat there, so only the edge
# carries duals.
function _corridor_weight(x::Vec3, lines::Vector{Tuple{SVector{3,S},Vec3}}, lo::Float64, hi::Float64,
                          edge::Float64) where {S<:Real}
    q = one(S)
    for (a, k) in lines
        v = x - ForwardDiff.value.(a)
        d2 = dot(v, v) - dot(v, k)^2
        d2 >= hi^2 && continue
        d2 <= lo^2 && return one(S)
        w = x - a
        d = sqrt(max(dot(w, w) - dot(w, k)^2, zero(S)))
        u = (hi - d) / edge
        q *= 1 - u^3 * (10 - u * (15 - 6u))
    end
    return 1 - q
end

function _corridor_mean(P::AbstractVector{T}, t::Target, es, lines::Vector{Tuple{SVector{3,S},Vec3}},
                        m::_PicketMeasure) where {T,S<:Real}
    isempty(lines) && throw(ArgumentError("no beams with power in the inner groups $(Tuple(m.inner)) to center the corridor on"))
    lo, hi = m.corridor - m.edge / 2, m.corridor + m.edge / 2
    R = promote_type(T, S)
    num, den = zero(R), zero(R)
    for e in es
        w = _corridor_weight(t.elem_centroid[e], lines, lo, hi, m.edge)
        iszero(w) && continue      # a zero value has zero partials (s is flat at the edge)
        num += w * P[e]
        den += w * t.elem_area[e]
    end
    den > 0 || throw(ArgumentError("no element of the parts lies within $(m.corridor) mm of an inner centerline"))
    return num / den
end

"""
    picket_hot_elements(ev, names...; n = 8, pulses = nothing, window = nothing, by = :group,
                        power = nothing, groups = (:outer44, :outer50), lpi = false,
                        inner = (:inner23, :inner30)) -> Vector{NamedTuple}

The `n` elements of the named parts or components with the largest irradiance from the
beams in `groups`, hottest first. The powers and beams are as in
[`picket_irradiance`](@ref). Each entry gives the index of the element (`element`), its
`irradiance` (TW/mm²), centroid `z` (mm), azimuth `phi` (degrees, 0–360) and
`distance` (mm). The distance runs from the centroid to the nearest centerline of a beam
in `inner` with power (`Inf` without one), and it is the distance weighted by the
`:corridor` measure of [`picket_irradiance`](@ref). The report is in Float64, with the
values of any duals.
"""
function picket_hot_elements(ev::Evaluation, names::Symbol...; n::Integer = 8, pulses = nothing,
                             window = nothing, by = :group, power = nothing, groups = (:outer44, :outer50),
                             lpi::Bool = false, inner = (:inner23, :inner30))
    n >= 1 || throw(ArgumentError("n must be at least 1, got $n"))
    t = ev.scenario.target
    es = _named_elements(t, names)
    _check_lpi(lpi, pulses, power)
    P = Float64.(ForwardDiff.value.(_picket_map(ev, _picket_power(pulses, window, by, power, true), _symbols(groups), lpi)))
    lines = [(Vec3(ForwardDiff.value.(a)), k) for (a, k) in _centerlines(ev.beams, collect(Symbol, _symbols(inner)))]
    I = [P[e] / t.elem_area[e] for e in es]
    hot = partialsortperm(I, 1:min(n, length(es)); rev = true)
    return map(hot) do i
        c = t.elem_centroid[es[i]]
        d = minimum((v = c - a; sqrt(max(dot(v, v) - dot(v, k)^2, 0.0))) for (a, k) in lines; init = Inf)
        (element = es[i], irradiance = I[i], z = c[3], phi = mod(atand(c[2], c[1]), 360.0), distance = d)
    end
end

"""
    picket_ceiling(reference, names...; pulses = nothing, window = nothing, by = :group,
                   power = nothing, groups = (:outer44, :outer50), lpi = false,
                   measure = :peak, area = 53.0, corridor = 1.3, edge = 0.3,
                   inner = (:inner23, :inner30),
                   limit = picket_irradiance(reference, reference_names...; ...), p = 64,
                   scale = limit, exact = true, reference_names = names) -> Constraint

A ceiling on the laser irradiance of the beams in `groups` over the elements of the
named parts or components, as a [`Constraint`](@ref) for [`optimize_pointing`](@ref)
named `:picket_<measure>`. The irradiance is that of [`picket_irradiance`](@ref), whose
keywords set the powers, the beams and the measure, the same for the reference and the
design. By default the limit is the value of the `reference` evaluation; a number in
TW/mm² sets another limit.

- `measure = :peak`: the maximum irradiance over the elements. The optimizer works on
  the smooth area-weighted `p`-norm of the irradiance, scaled to the peak of the
  reference, as [`irradiance_ceiling`](@ref) does at the time of the scenario. With
  `exact = true` (the default) the run enforces the peak itself. The working limit of
  the stand-in is lowered while the peak exceeds the limit, and the peak is reported
  against the limit ([`check_constraints`](@ref), `PointingResult.constraints`).
- `measure = :weighted`: the energy-weighted mean irradiance, the quantity that
  [`picket_ratio`](@ref) carries into the tip speed of [`BubbleAbsorption`](@ref)
  (∝ ratio^0.2). It is smooth and enforced directly.
- `measure = :average`: the spatially averaged irradiance over a fixed area, the
  quantity of the published bubble scalings. Against the reference it is the absorbed
  power on the parts over that of the reference. It is linear in the powers and enforced
  directly.
- `measure = :corridor`: the mean irradiance near the centerlines of the inner beams of
  the design, which move with its aims. It is smooth and enforced directly.

**A reference on another target.** For `:peak` the reference evaluation must be on the
same target as the optimized scenario; a [`subset`](@ref) of the scenario shares the
target. Otherwise evaluating the constraint throws an `ArgumentError`. The reason is
that the peak follows one element and the mesh, and the scale of the stand-in is the
ratio of peak to norm on the reference mesh. The three integral measures are valid
across meshes and accept a reference on another target, such as a hohlraum without the
windows of the design or a separate mesh for the reference. They are:

- `:average`, the absorbed power over a fixed area, which is a quantity of the field and
  does not depend on the mesh;
- `:weighted`, an integral that converges as O(h²) with the element size h. Targets
  should be meshed at the same `h`, where the bias of about 1% at 0.2 mm cancels in the
  ratio;
- `:corridor`, an area mean.

On an evaluation of another target, the constraint measures the parts `names` of that
target, and the limit is the value of the reference over `reference_names` (by default
`names`). A design with windows in its mesh can then hold the light on its wall and
windows at the light of the reference on its wall, with `picket_ceiling(reference,
:wall, :window; reference_names = (:wall,), measure = :average)`.

The laser–plasma model is used only with `lpi = true`. The reference, on the same target
or another, must then have the model and time of the scenario, as for
[`irradiance_ceiling`](@ref) and [`mode_gate`](@ref). Evaluating the constraint with
another model, or at another time with a model, throws an `ArgumentError`. The
constraint performs one deposition per evaluation, in the number type of the
evaluation. Inside [`optimize_pointing`](@ref), with the beams at their present powers (no
`pulses` or `power`), no deposition is needed. The element powers of the groups are then
kept from the same pass of the deposition in the evaluation, with identical sums. The
exceptions are a laser–plasma model of the scenario that deposits light itself, and one that
changes the beams with `lpi = false`. For a phase other than that of the optimized scenario,
such as a trough or a foot, the powers of the phase are given by `pulses` and `window` or by
`power`, or the constraint is built on an evaluation of that phase.

*Derived*; no facility criterion sets the limit. The picket fluence is the measure of
the published bubble scalings, and holding it at the value of a reference was chosen for
a re-pointing study (user, 2026-11-18). That study limits first its peak, as an
engineering limit on the hottest element of the picket, then the published `:average` in
the picket, the trough and the second step, with the `:corridor` reported. The model
predicts the bubble of a design from the ratio of the means to the 0.2 power
([`PicketLinkedBubble`](@ref), [`bubble_prediction`](@ref)), with low confidence
(research notes, "The swing and the bubble without new data (2026-11-18)" §5, and "The
published picket fluence, the trough and the inner-beam corridor (2026-11-18)").
"""
function picket_ceiling(reference::Evaluation, names::Symbol...; pulses = nothing, window = nothing,
                        by = :group, power = nothing, groups = (:outer44, :outer50), lpi::Bool = false,
                        measure::Symbol = :peak, area::Real = 53.0, corridor::Real = 1.3, edge::Real = 0.3,
                        inner = (:inner23, :inner30), limit::Union{Nothing,Real} = nothing, p::Real = 64,
                        scale::Union{Nothing,Real} = nothing, exact::Bool = true,
                        reference_names = names)
    mo = _picket_measure_options(measure, area, corridor, edge, inner)
    _check_lpi(lpi, pulses, power)
    t = reference.scenario.target
    rn = Tuple(Symbol.(collect(reference_names)))
    same = rn == names
    isempty(names) && throw(ArgumentError("name the parts (e.g. :wall)"))
    res = _named_elements(t, rn)
    pw = _picket_power(pulses, window, by, power, true)
    gs = collect(Symbol, _symbols(groups))
    P = Float64.(ForwardDiff.value.(_picket_map(reference, pw, gs, lpi)))
    sum(P[e] for e in res) > 0 ||
        throw(ArgumentError("the reference's picket deposits no power on $(join(rn, ", "))"))
    ref = Float64(ForwardDiff.value(_picket_measure(P, t, res, mo, reference.beams)))
    factor = measure === :peak ? ref / _irradiance_norm(P, t, res, Float64(p)) : 1.0
    lim = limit === nothing ? ref : Float64(limit)
    lim > 0 || throw(ArgumentError("the picket irradiance limit must be positive"))
    # parts of the design on the target of the reference, resolved at evaluation when the
    # targets differ
    es = same ? res : Int[]
    value = _PicketValue(t, es, collect(Symbol, names), gs, pw, lpi, mo, Float64(p), factor, exact,
                         lpi ? _ReferenceModel(reference.scenario) : nothing)
    return Constraint(Symbol(:picket_, measure), value; upper = lim, scale = something(scale, lim))
end

# The value of the picket ceiling. For `:peak` it is the p-norm of the picket irradiance
# on the peak scale of the reference; otherwise it is the measure itself. With `exact`,
# the reported value of `:peak` is the peak itself (`_exact_value`). `elements` are those
# of `names` on the target of the reference, or empty when they are resolved at
# evaluation. With `lpi`, `model` is the laser–plasma model and time of the reference,
# checked on an evaluation of another target; without `lpi` it is `nothing`.
struct _PicketValue{F,M<:Union{Nothing,_ReferenceModel}}
    target::Target
    elements::Vector{Int}
    names::Vector{Symbol}
    groups::Vector{Symbol}
    power::F
    lpi::Bool
    options::_PicketMeasure
    p::Float64
    factor::Float64
    exact::Bool
    model::M
end

# Returns the target and elements that the ceiling measures at `ev`. They are the target
# of the reference and its elements, or, for a measure valid across meshes, the target of
# the evaluation and the elements of the named parts on it. With `lpi`, the laser–plasma
# model and time of the evaluation must be those of the reference on either target.
function _picket_parts(c::_PicketValue, ev::Evaluation)
    t = ev.scenario.target
    t === c.target || c.options.measure !== :peak ||
        throw(ArgumentError("picket_ceiling: the evaluation is on another target than the reference: a ceiling " *
                            "on the peak holds on the reference's mesh only (the measures :average, :weighted " *
                            "and :corridor take a reference on another target)"))
    _check_model(c.model, "picket_ceiling", ev)
    t === c.target && return t, isempty(c.elements) ? _part_elements(t, c.names) : c.elements
    return t, _part_elements(t, c.names)
end

# `_named_elements` for a vector of names (`elements(t, names...)`, in its order, without
# the splat).
function _part_elements(t::Target, names::Vector{Symbol})
    _check_names(t, names)
    return reduce(vcat, (t.comp_elems[ci] for n in names for ci in part_pieces(t, n)); init = Int[])
end

function (c::_PicketValue)(ev::Evaluation)
    t, es = _picket_parts(c, ev)
    P = _picket_map(ev, c.power, c.groups, c.lpi)
    return c.options.measure === :peak ? c.factor * _irradiance_norm(P, t, es, c.p) :
           _picket_measure(P, t, es, c.options, ev.beams)
end

function _exact_value(c::_PicketValue, ev::Evaluation)
    c.exact && c.options.measure === :peak || return nothing
    t, es = _picket_parts(c, ev)
    return _peak(Float64.(ForwardDiff.value.(_picket_map(ev, c.power, c.groups, c.lpi))), t, es)
end

# Returns the group set that a ceiling uses, for the split of the optimizer
# (`_with_split`). The set applies only to beams at their present powers, since other
# powers are not those of the evaluation.
_split_sets(c::_PicketValue) =
    c.power isa _OwnPower ? _SplitRequest([nothing => _group_set(c.groups)]) : _SplitRequest()

# ---------------------------------------------------------------------------------------
# The bubble that follows the picket fluence of the pointing
# ---------------------------------------------------------------------------------------

"""
    PicketLinkedBubble(bubble, target, reference; pulses = nothing, window = nothing,
                       by = :group, power = nothing, groups = (:outer44, :outer50),
                       components = (:wall,), deposition_kw = (;))

A [`BubbleAbsorption`](@ref) whose picket fluence ratio E follows the beams it is
applied to. The model is an opt-in [`LaserPlasmaModel`](@ref). For beams `bs`,

    E(bs) = bubble.picket_ratio × F̄(bs) / F̄(reference),

with F̄ the energy-weighted mean outer picket fluence on the `components` of `target`
([`picket_fluence`](@ref)). The bubble is `remake(bubble; picket_ratio = E(bs))`, with
tip speed ∝ E^0.2 and collision time − `turn_on` ∝ E^-0.16. `bubble.picket_ratio` is
the E of the `reference` pointing, set by hand or from [`picket_ratio`](@ref) against
the calibration platform. At the reference the model is therefore identical to
`bubble`. A re-pointed design is given the bubble that its picket grows, inside
`evaluate(scenario, beams)`, the constraints and objectives built on it
([`mode_gate`](@ref)) and [`drive_history`](@ref)`(...; lpi)`. The model composes as
`bubble` does, for example `lpi = (cbet, linked)`.

The picket map is computed on `target`. It is the hohlraum on a mesh fine enough for the map
(h ≤ 0.2 mm), a separate target used only for the deposition; the radiation mesh of a
drive history is usually too coarse (research notes). A deposition operator is built on
it for the `reference` beams in `groups` (with `deposition_kw`). The beams that the
model is applied to are matched to those beams by name, and they must enter through the
same ports. A beam of `groups` not among them is an error, and a missing one, as in a
dropped bundle, deposits nothing. The picket powers are those of
[`picket_irradiance`](@ref): `pulses` over a `window`, for example the mean power over
the first ns after the outer turn-on, or `power`, a function of the beam. One of the two
is required. It must not use the `power` of the beam, which at the `time` of the model
is not the picket power; in a chain after CBET it is not even the setup power. E depends
only on the pointing, spots and focus of the beams, so it is constant through a drive
history.

Cost: two picket depositions per evaluation, one in [`apply_lpi`](@ref) and one in
[`lpi_deposit`](@ref). Generic in the number type. E and the bubble carry the
ForwardDiff duals of the beams, so a pointing gradient includes the response of the
bubble to the picket fluence. Thread-safe, since the model has no state. The per-beam
sensitivities of the uncertainty analysis ([`mode_sensitivity`](@ref),
[`mode_uncertainty`](@ref)) deposit beams one at a time and cannot resolve E, so they
throw an `ArgumentError`. The resolved bubble [`linked_bubble`](@ref)`(model, beams)` is
passed to them instead.

*Derived* (confidence low). The exponents are published for picket energy and quad
splitting on one platform: E^0.2 from Izumi et al., PoP 28, 022706 (2021), on seven
shots, and E^-0.16 from Izumi et al., RSI 89, 10K111 (2018). Their use for a fluence
changed by re-pointing at fixed energy has no published test. F̄ counts the overlap of
the spots, as the split-quad shot of Izumi 2021 does. The energy over a fixed area used
in the papers (the `:average` of [`picket_irradiance`](@ref)) does not, and it differs
from F̄ on uneven powers. See the research notes, "The swing and the bubble without new
data (2026-11-18)" §5, and "The published picket fluence, the trough and the inner-beam
corridor (2026-11-18)" §5, which measures the link through the fixed area instead. No
facility criterion sets the model, and no default uses it.
"""
struct PicketLinkedBubble{S<:Real,F} <: LaserPlasmaModel
    bubble::BubbleAbsorption{S}
    operator::DepositionOperator
    beams::Vector{Beam{Float64}}        # beams of the operator, those of the reference in `groups`
    index::Dict{String,Int}             # operator index of each beam name
    power::F                            # picket power per beam (TW)
    groups::Vector{Symbol}
    elements::Vector{Int}
    reference::PicketFluence{Float64}   # map and summaries of the reference
end

function PicketLinkedBubble(bubble::BubbleAbsorption, target::Target, reference::AbstractVector{<:Beam};
                            pulses = nothing, window = nothing, by = :group, power = nothing,
                            groups = (:outer44, :outer50), components = (:wall,), deposition_kw = (;))
    pw = _picket_power(pulses, window, by, power, false)
    gs = collect(Symbol, _symbols(groups))
    rb = Beam{Float64}[b for b in reference if b.group in gs]
    isempty(rb) && throw(ArgumentError("no reference beams in the groups $(Tuple(gs))"))
    names = [b.name for b in rb]
    allunique(names) || throw(ArgumentError("the reference beams' names must be unique"))
    op = DepositionOperator(target, rb; deposition_kw...)
    es = _picket_elements(target, components)
    index = Dict(n => i for (i, n) in enumerate(names))
    ref = _linked_fluence(op, rb, index, pw, gs, es, reference)
    return PicketLinkedBubble(bubble, op, rb, index, pw, gs, es, ref)
end

function Base.show(io::IO, m::PicketLinkedBubble)
    print(io, "PicketLinkedBubble(E = ", ForwardDiff.value(m.bubble.picket_ratio), " at the reference, ",
          length(m.beams), " ", join(m.groups, "/"), " beams, ", m.reference, ")")
end

# Returns the picket fluence of `beams` on the map of the model. Each beam of the groups
# of the model is placed at its picket power on the operator beam of the same name, and
# the other beams of the operator are at zero.
_linked_fluence(m::PicketLinkedBubble, beams) =
    _linked_fluence(m.operator, m.beams, m.index, m.power, m.groups, m.elements, beams)

function _linked_fluence(op, rb, index, power, gs, es, beams::AbstractVector{Beam{T}}) where {T<:Real}
    picket = Vector{Beam{T}}(undef, length(rb))
    found = falses(length(rb))
    for b in beams
        b.group in gs || continue
        i = get(index, b.name, 0)
        i == 0 && throw(ArgumentError("beam $(b.name) is not among the linked bubble's reference beams"))
        picket[i] = remake(b; power = oftype(b.power, power(b)))
        found[i] = true
    end
    for i in eachindex(rb)
        found[i] || (picket[i] = remake(rb[i]; power = zero(T)))
    end
    return _picket_fluence(deposit(op, picket).element, op.target, es)
end
_linked_fluence(op, rb, index, power, gs, es, beams::AbstractVector{<:Beam}) =
    _linked_fluence(op, rb, index, power, gs, es, _common_type(beams))

# Converts a vector of beams of mixed number types to one type, the promotion of their
# number types, which is abstract if one of them is.
_number_type(::Beam{T}) where {T} = T
function _common_type(beams)
    T = mapreduce(_number_type, promote_type, beams)
    return Beam{T}[Beam{T}(b.name, b.port, SVector{3,T}(b.aim), T(b.power), T(b.defocus), _retype(b.spot, T),
                           b.fnumber, b.group, b.quad, b.wavelength, T(b.dlambda)) for b in beams]
end

"""
    linked_bubble(m::PicketLinkedBubble, beams) -> BubbleAbsorption

The bubble of a [`PicketLinkedBubble`](@ref) for `beams`, which is `remake(m.bubble;
picket_ratio = E)` with E from the picket fluence of the beams relative to that of the
reference. At the reference it is `m.bubble` itself. The result is in the number type of
the beams.
"""
linked_bubble(m::PicketLinkedBubble, beams::AbstractVector{<:Beam}) =
    remake(m.bubble; picket_ratio = picket_ratio(_linked_fluence(m, beams), m.reference;
                                                 at = m.bubble.picket_ratio))

apply_lpi(m::PicketLinkedBubble, beams::AbstractVector{<:Beam}, time = 0.0) =
    apply_lpi(linked_bubble(m, beams), beams, time)

_deposits(::PicketLinkedBubble) = true

_lpi_deposit(m::PicketLinkedBubble, op::DepositionOperator, beams::AbstractVector{<:Beam}, time, split) =
    _lpi_deposit(linked_bubble(m, beams), op, beams, time, split)

_ridge_split(m::PicketLinkedBubble, op::DepositionOperator, b::Beam, time) =
    throw(ArgumentError("a PicketLinkedBubble resolves its picket ratio from all beams together, and the " *
                        "per-beam sensitivities deposit one beam at a time: use the resolved bubble, " *
                        "`linked_bubble(model, beams)`, in the scenario's laser–plasma model"))

bubble_fractions(m::PicketLinkedBubble, beams::AbstractVector{<:Beam}, t) =
    bubble_fractions(linked_bubble(m, beams), beams, t)

# ---------------------------------------------------------------------------------------
# What a re-pointing does to the bubble, and through it to the swing
# ---------------------------------------------------------------------------------------

"""
    bubble_prediction(model::PicketLinkedBubble, design, reference; history = nothing,
                      reference_history = nothing, fixed = false, swing = (4.5, 9.0))
        -> NamedTuple

The prediction of the gold bubble model for a re-pointed `design` (beams) against a
`reference` pointing. Both are mapped by `model`, a [`PicketLinkedBubble`](@ref) with its
picket powers, target and bubble. Usually `reference` is the pointing on which the model
was built. The fields are:

- `peak_ratio`: the maximum picket fluence (equivalently irradiance) over the elements of
  the model, design over reference (the `peak` of [`picket_fluence`](@ref));
- `fluence_ratio`: the energy-weighted mean picket fluence, design over reference;
- `picket_ratio`, `reference_picket_ratio`: E of each ([`linked_bubble`](@ref));
- `speed_ratio`: the tip speed of the bubble, design over reference, `fluence_ratio^0.2`;
- `collision`, `reference_collision`: the ablator–wall collision times (ns);
- `bubble`, `reference_bubble`: the two [`BubbleAbsorption`](@ref)s; `fluence`,
  `reference_fluence`: the two [`PicketFluence`](@ref) maps;
- `history`: `nothing`, or, with `history`, the P2/P0 of the capsule drive through time
  with the bubble linked.

`history(beams, bubble)` runs a drive history of `beams` (the result of
[`drive_history`](@ref), with `lmax ≥ 2`) with the laser–plasma chain built on `bubble`,
for example
`(bs, bub) -> drive_history(model, τ -> with_pulses(bs, pulses, τ), times; lpi = (cbet, bub))`.
It is run for the design with its bubble and for the reference with a separate bubble.
`reference_history`, a result of the same call for the reference, replaces that second
run. With `fixed = true` it is also run for the design with the bubble of the reference,
which splits the effect of the pointing into the parts through the deposition and
through the bubble. `history` then contains `time`, `p2`, `reference_p2` and `fixed_p2`
(or `nothing`). For each run (`design`, `reference`, `fixed`) it gives the `minimum` of
P2/P0 in the swing over the `swing` window (ns), the time at which it is reached (`at`),
and the upward zero `crossing` after it (linear between steps; `NaN` if none). It also
gives `shift`, the value of the design minus that of the reference, and with `fixed`,
`through_bubble` (design − fixed) and `through_deposition` (fixed − reference).

*Derived* (confidence low). The response of the bubble to the pointing rests on the
published exponents (tip speed ∝ E^0.2, Izumi 2021; collision ∝ E^-0.16, Izumi 2018)
applied to a fluence changed by re-pointing, which has no published test
([`PicketLinkedBubble`](@ref)). Ratios between pointings on one mesh converge to < 1%
(research notes, "The swing and the bubble without new data (2026-11-18)" §5).
"""
Base.@nospecializeinfer function bubble_prediction(m::PicketLinkedBubble, design::AbstractVector{<:Beam},
                                                   reference::AbstractVector{<:Beam};
                                                   @nospecialize(history = nothing),
                                                   @nospecialize(reference_history = nothing),
                                                   fixed::Bool = false, swing = (4.5, 9.0))
    fd, fr = _linked_fluence(m, design), _linked_fluence(m, reference)
    at = m.bubble.picket_ratio
    Ed, Er = picket_ratio(fd, m.reference; at), picket_ratio(fr, m.reference; at)
    bd, br = remake(m.bubble; picket_ratio = Ed), remake(m.bubble; picket_ratio = Er)
    h = history === nothing ? nothing :
        _swing_history(history, design, reference, bd, br, reference_history, fixed, swing)
    return (; peak_ratio = fd.peak / fr.peak, fluence_ratio = fd.mean / fr.mean, picket_ratio = Ed,
            reference_picket_ratio = Er, speed_ratio = bubble_speed(bd) / bubble_speed(br),
            collision = bd.collision, reference_collision = br.collision, bubble = bd,
            reference_bubble = br, fluence = fd, reference_fluence = fr, history = h)
end

Base.@nospecializeinfer function _swing_history(@nospecialize(history), design, reference, bd, br,
                                                @nospecialize(reference_history), fixed::Bool, window)
    p2(r) = Float64.(ForwardDiff.value.(r.moments[3, :] ./ r.moments[1, :]))
    rd = history(design, bd)
    rr = reference_history === nothing ? history(reference, br) : reference_history
    rf = fixed ? history(design, br) : nothing
    t = collect(Float64, rd.time)
    collect(Float64, rr.time) == t && (rf === nothing || collect(Float64, rf.time) == t) ||
        throw(ArgumentError("the histories must share their times"))
    pd, pr = p2(rd), p2(rr)
    pf = rf === nothing ? nothing : p2(rf)
    sd, sr = _swing(t, pd, window), _swing(t, pr, window)
    sf = pf === nothing ? nothing : _swing(t, pf, window)
    diff(a, b) = (minimum = a.minimum - b.minimum, at = a.at - b.at, crossing = a.crossing - b.crossing)
    return (; time = t, p2 = pd, reference_p2 = pr, fixed_p2 = pf,
            swing = (; design = sd, reference = sr, fixed = sf), shift = diff(sd, sr),
            through_bubble = sf === nothing ? nothing : diff(sd, sf),
            through_deposition = sf === nothing ? nothing : diff(sf, sr))
end

# Returns the swing of a P2/P0 history: its minimum over the window, the time of the
# minimum, and the first upward zero crossing after it (linear between the steps; NaN if
# none).
function _swing(t::Vector{Float64}, p::Vector{Float64}, window)
    ks = [k for k in eachindex(t) if window[1] <= t[k] <= window[2] && isfinite(p[k])]
    isempty(ks) && return (minimum = NaN, at = NaN, crossing = NaN)    # no step in the window
    km = ks[argmin(p[ks])]
    kz = findfirst(k -> k > km && p[k-1] < 0 <= p[k], eachindex(p))
    tz = kz === nothing ? NaN : t[kz-1] - p[kz-1] * (t[kz] - t[kz-1]) / (p[kz] - p[kz-1])
    return (minimum = p[km], at = t[km], crossing = tz)
end
