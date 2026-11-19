# The outer cones' picket fluence on the wall (energy per area), the quantity the
# published gold-bubble scalings are written in, and `BubbleAbsorption`'s picket ratio
# from it; the picket irradiance as a constraint, the bubble that follows the beams'
# picket, and the prediction for a re-pointed design. Research notes, "The swing and the
# bubble without new data (2026-11-18)", §5 "The bubble's speed from the pointing's
# picket fluence"; literature in `reference/litreview/bubble_drivers.md`.

"""
    PicketFluence{T}

The picket fluence of a set of beams on a target's wall, from [`picket_fluence`](@ref):
- `fluence`: per element of the target (kJ/mm²; ×1e5 → J/cm²), zero outside `elements`;
- `elements`: the elements counted (those of the named components);
- `energy`: the picket energy deposited on them, ∫F dA (kJ); over a fixed spot area it
  is the spatially averaged fluence of the published bubble scalings (its ratio
  between designs is the energy's: [`picket_irradiance`](@ref)'s `:average`);
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

**The summaries.** Every published bubble scaling is written in the outer cones' picket
fluence as a nominal energy over a nominal area: Ralph et al., PoP 25, 082701 (2018),
the outer picket energy net of backscatter over the outer phase-plate spot area (speed
∝ I^0.36); Callahan et al., PoP 25, 056305 (2018), the outer picket energy over the
one-side outer spot area, 53 mm² (size ∝ (E/A)^0.5); Izumi et al., PoP 28, 022706
(2021), the spatially averaged outer picket fluence over the first ns (tip position
∝ E^0.2), 27 % lower when the quads and cones are split at the same energy. Over a fixed
area that quantity is `energy` (its ratios between designs are the energy's;
[`picket_irradiance`](@ref)'s `:average`). `mean` = ∫F² dA / ∫F dA, the energy over the
effective area (∫F dA)²/∫F² dA, counts the overlap of the spots as well: one isolated
spot gives its energy over its effective area (a flat top's E/A exactly; for a
super-Gaussian of order n, 2^(-2/n) of its peak fluence, half for a Gaussian), k
coincident spots k times one, separated spots the same as one, so it carries Izumi
2021's split-quad change, which a fixed area does not; as an integral it is smooth in
the pointing and nearly independent of the mesh. The two part when the beams' picket
energies are uneven or their spots move onto one another (research notes, "The
published picket fluence, the trough and the inner-beam corridor (2026-11-18)").
`peak` follows one element and the mesh; it measures neither.

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
    at * (design.mean / reference.mean)     # `at` itself, exactly, at the reference

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

# The measures of a picket map and their options (`picket_irradiance`'s keywords).
struct _PicketMeasure
    measure::Symbol
    area::Float64           # `:average`: the fixed spot area (mm²)
    corridor::Float64       # `:corridor`: the distance from an inner centerline (mm) …
    edge::Float64           # … and the width of the weight's soft edge (mm)
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

# `lpi = true` reads the evaluation's effective beams, at their own powers.
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
outer cones by default) on the elements of the named parts or components, summarized by
`measure`:

- `:peak`: the largest over the elements;
- `:weighted`: the energy-weighted mean ∫I² dA / ∫I dA, the summary that
  [`picket_ratio`](@ref) and [`PicketLinkedBubble`](@ref) read ([`picket_fluence`](@ref));
- `:average`: the spatially averaged irradiance over a fixed spot area, the power
  absorbed on the elements over `area` (mm²): the quantity of the published bubble
  scalings (below);
- `:corridor`: the mean irradiance on the elements near the centerlines of the beams in
  `inner`, ∫w I dA / ∫w dA with a weight w that is 1 within `corridor − edge/2` of a
  centerline and 0 beyond `corridor + edge/2` (below).

The irradiance times the phase's length is its fluence, so the ratios between designs
are the fluence's. The powers are those of whatever phase the evaluation stands for:
the picket, the trough or the foot.

**The powers.** With `pulses` (keyed as for [`with_pulses`](@ref), with `by`) and
`window = (t0, t1)` (ns), each beam's mean power over the window; with `power`, a
function of the beam that returns its power (TW), e.g. by cone from a picket table; with
neither, each beam's own power (an evaluation run at a phase's powers, or `power = b ->
r * b.power` to scale them). The other beams are deposited at zero power, by the
scenario's deposition operator with the evaluation's setup beams (`ev.beams`: their
pointing, spots and focus). With `lpi = true` the beams are instead the evaluation's
effective beams (`ev.effective`, after the scenario's laser–plasma model at its own
`time`, e.g. net of a [`Backscatter`](@ref) and after a transfer), at their own powers
(`pulses` and `power` are then an error); otherwise the laser–plasma model does not
enter, and backscatter or transfer is taken off the powers given. One deposition per
call. Generic in the number type (ForwardDiff duals through the pointing and powers).

**The published quantity** (`:average`). Callahan et al., PoP 25, 056305 (2018) divide
the outer cones' picket energy by the one-side outer spot area of NIF's phase plates,
53 mm² (the default `area`, so that the irradiance times the picket's length in ns is
his E/A in kJ/mm²); Ralph et al., PoP 25, 082701 (2018) take the outer picket energy net
of backscatter over the outer spot area; Izumi et al., PoP 28, 022706 (2021) the
spatially averaged outer picket fluence over the first ns. In Briareus terms: the power
*absorbed on the named parts* (the wall: light clipped on the entrance holes' outside,
landing on the capsule or leaving the target does not grow the wall's bubble, *L*),
*net of backscatter* when the powers are (`lpi = true` with a backscatter model, or
powers net of it), over a fixed area. A ratio to a reference does not depend on the
area; an absolute value takes Callahan's convention. A fixed area does not count the
overlap of the spots, which `:weighted` does: the two part on uneven powers between beams
and on spots moved onto one another (research notes, "The published picket fluence,
the trough and the inner-beam corridor (2026-11-18)").

**The inner-beam corridor** (`:corridor`). Izumi et al., PoP 30, 032706 (2023) found
outer bubble spikes in line with the 23.5° beams to clip them and turn the core more
oblate; the corridor measures the outer light where the inner beams pass. A centerline
is the line through a beam's aim along its direction, for each beam of `inner` with
nonzero power in the evaluation; the distance is the element centroid's, in 3-D, to
the nearest one. The weight is one minus the product over the centerlines of
(1 − s), s rising from 0 to 1 across the edge as a quintic smoothstep: smooth (C²) in
the inner beams' aims as they move, exactly 0 or 1 away from the edge (elements beyond
every edge are skipped). `corridor = 1.3` mm and `edge = 0.3` mm are engineering
defaults: 1.3 mm is a re-pointing study's choice (user, 2026-11-18), and the edge is
wider than an element of a 0.2 mm deposition map, so that the mean moves smoothly as
centroids cross it; neither is a published criterion.

The peak follows one element and the mesh (on a coarser mesh the hottest spot is
averaged over more area); the other three are integrals, nearly independent of the mesh
([`picket_fluence`](@ref)). [`picket_hot_elements`](@ref) lists the hottest elements
with their distance to the nearest inner centerline.
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

# The element powers (TW) of the evaluation's setup beams (or with `lpi` its effective
# beams) at their picket powers.
_picket_map(ev::Evaluation, power, gs, lpi::Bool) =
    lpi ? _picket_map(ev.scenario.deposition, power, gs, ev.effective) :
          _picket_map(ev.scenario.deposition, power, gs, ev.beams)
_picket_map(op::DepositionOperator, power, gs, beams) = deposit(op, _picket_beams(power, gs, beams)).element

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

# The corridor's weight at x: 1 − Π (1 − s(u)) over the centerlines, s the quintic
# smoothstep of u = (hi − d)/edge for a distance d between lo = hi − edge and hi. Each
# distance is classified on its value first: beyond hi a factor 1, within lo the weight 1
# (both exact: s is flat there), so only the edge carries duals.
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
beams in `groups` (the powers and beams as in [`picket_irradiance`](@ref)), hottest
first: per element its index (`element`), `irradiance` (TW/mm²), centroid `z` (mm) and
azimuth `phi` (degrees, 0–360), and `distance` (mm), from its centroid to the nearest
centerline of a beam in `inner` with power (`Inf` without one), the distance that
[`picket_irradiance`](@ref)'s `:corridor` weighs. A report in Float64 (the values of
duals).
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
                   limit = picket_irradiance(reference, names...; ...), p = 64,
                   scale = limit, exact = true) -> Constraint

A ceiling on the laser irradiance of the beams in `groups` over the elements of the
named parts or components ([`picket_irradiance`](@ref), whose keywords set the powers,
the beams and the measure, the same for the reference and the design), at or below the
`reference` evaluation's by default (a number in TW/mm² sets another limit): a
[`Constraint`](@ref) for [`optimize_pointing`](@ref), named `:picket_<measure>`.

- `measure = :peak`: the maximum irradiance over the elements. The optimizer works on
  the smooth area-weighted `p`-norm of the irradiance scaled to the reference's peak
  (as [`irradiance_ceiling`](@ref) does at the scenario's time); with `exact = true`
  (the default) the run holds the peak itself, lowering the stand-in's working limit
  while the peak exceeds the limit, and reports the peak against the limit
  ([`check_constraints`](@ref), `PointingResult.constraints`).
- `measure = :weighted`: the energy-weighted mean irradiance, the quantity that
  [`picket_ratio`](@ref) carries into [`BubbleAbsorption`](@ref)'s tip speed
  (∝ ratio^0.2): smooth, held directly.
- `measure = :average`: the spatially averaged irradiance over a fixed area, the
  quantity of the published bubble scalings (against the reference: the absorbed power
  on the parts over the reference's); linear in the powers, held directly.
- `measure = :corridor`: the mean irradiance near the inner beams' centerlines, the
  design's own (they move with its aims): smooth, held directly.

The reference evaluation must be on the same target as the optimized scenario (a
[`subset`](@ref) of it shares the target), or evaluating the constraint throws an
`ArgumentError`; its laser–plasma model and time enter only with `lpi = true`. One
deposition per evaluation, in the evaluation's number type. For a phase other than the
optimized scenario's (a trough, a foot), give the phase's powers by `pulses` and
`window` or `power`, or build the constraint on an evaluation of that phase.

*Derived*, not a facility criterion: the picket fluence is the measure of the published
bubble scalings, and holding it at a reference's is a re-pointing study's choice (user,
2026-11-18): first its peak, as an engineering limit on the picket's hottest element,
then the published `:average` in the picket, the trough and the second step, with the
`:corridor` reported. The model's prediction of the bubble a design grows is the mean's
ratio to the 0.2 power ([`PicketLinkedBubble`](@ref), [`bubble_prediction`](@ref)),
confidence low (research notes, "The swing and the bubble without new data
(2026-11-18)" §5, and "The published picket fluence, the trough and the inner-beam
corridor (2026-11-18)").
"""
function picket_ceiling(reference::Evaluation, names::Symbol...; pulses = nothing, window = nothing,
                        by = :group, power = nothing, groups = (:outer44, :outer50), lpi::Bool = false,
                        measure::Symbol = :peak, area::Real = 53.0, corridor::Real = 1.3, edge::Real = 0.3,
                        inner = (:inner23, :inner30), limit::Union{Nothing,Real} = nothing, p::Real = 64,
                        scale::Union{Nothing,Real} = nothing, exact::Bool = true)
    mo = _picket_measure_options(measure, area, corridor, edge, inner)
    _check_lpi(lpi, pulses, power)
    t = reference.scenario.target
    es = _named_elements(t, names)
    pw = _picket_power(pulses, window, by, power, true)
    gs = collect(Symbol, _symbols(groups))
    P = Float64.(ForwardDiff.value.(_picket_map(reference, pw, gs, lpi)))
    sum(P[e] for e in es) > 0 ||
        throw(ArgumentError("the reference's picket deposits no power on $(join(names, ", "))"))
    ref = Float64(ForwardDiff.value(_picket_measure(P, t, es, mo, reference.beams)))
    factor = measure === :peak ? ref / _irradiance_norm(P, t, es, Float64(p)) : 1.0
    lim = limit === nothing ? ref : Float64(limit)
    lim > 0 || throw(ArgumentError("the picket irradiance limit must be positive"))
    value = _PicketValue(t, es, gs, pw, lpi, mo, Float64(p), factor, exact)
    return Constraint(Symbol(:picket_, measure), value; upper = lim, scale = something(scale, lim))
end

# The picket ceiling's value: the picket irradiance's p-norm on the reference's peak
# scale (`:peak`), or the measure itself; with `exact`, the reported value of `:peak` is
# the peak itself (`_exact_value`).
struct _PicketValue{F}
    target::Target
    elements::Vector{Int}
    groups::Vector{Symbol}
    power::F
    lpi::Bool
    options::_PicketMeasure
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
    P = _picket_map(ev, c.power, c.groups, c.lpi)
    return c.options.measure === :peak ? c.factor * _irradiance_norm(P, c.target, c.elements, c.p) :
           _picket_measure(P, c.target, c.elements, c.options, ev.beams)
end

function _exact_value(c::_PicketValue, ev::Evaluation)
    c.exact && c.options.measure === :peak || return nothing
    _check_picket(c, ev)
    return _peak(Float64.(ForwardDiff.value.(_picket_map(ev, c.power, c.groups, c.lpi))), c.target, c.elements)
end

# ---------------------------------------------------------------------------------------
# The bubble that follows the pointing's picket fluence
# ---------------------------------------------------------------------------------------

"""
    PicketLinkedBubble(bubble, target, reference; pulses = nothing, window = nothing,
                       by = :group, power = nothing, groups = (:outer44, :outer50),
                       components = (:wall,), deposition_kw = (;))

A [`BubbleAbsorption`](@ref) whose picket fluence ratio E follows the beams it is
applied to (a [`LaserPlasmaModel`](@ref), opt-in): for beams `bs`,

    E(bs) = bubble.picket_ratio × F̄(bs) / F̄(reference),

with F̄ the energy-weighted mean outer picket fluence on the `components` of `target`
([`picket_fluence`](@ref)), and the bubble `remake(bubble; picket_ratio = E(bs))`
(tip speed ∝ E^0.2, collision time − `turn_on` ∝ E^-0.16). `bubble.picket_ratio` is
the E of the `reference` pointing (set by hand, or from [`picket_ratio`](@ref)
against the calibration platform), so at the reference the model is `bubble` itself,
bit for bit, and a re-pointed design gets the bubble its own picket grows, inside
`evaluate(scenario, beams)`, the constraints and objectives built on it
([`mode_gate`](@ref)) and [`drive_history`](@ref)`(...; lpi)`. Compose as `bubble` is
composed, e.g. `lpi = (cbet, linked)`.

`target` carries the picket map: the hohlraum on a mesh fine enough for it (h ≤ 0.2
mm; the radiation mesh of a drive history is usually too coarse; research notes), a
separate target for the deposition alone. A deposition operator is built on it for the
`reference` beams in `groups` (with `deposition_kw`); the beams the model is applied to
are matched to them by name (a beam of `groups` not among them is an error; one of them
missing, as a dropped bundle's, deposits nothing), and must enter through the same
ports. The picket powers are those of [`picket_irradiance`](@ref): `pulses` over a
`window` (the mean power over e.g. the first ns after the outer turn-on), or `power`, a
function of the beam; one of the two is required, and it must not read the beam's
`power`, which at the model's `time` is not the picket's (in a chain after CBET it is
not even the setup's). E depends on the beams' pointing, spots and focus alone, so it
is constant through a drive history.

Cost: two picket depositions per evaluation (one in [`apply_lpi`](@ref), one in
[`lpi_deposit`](@ref)). Generic in the number type: E and the bubble carry the beams'
ForwardDiff duals, so a pointing gradient includes the bubble's response to the picket
fluence. Thread-safe (no state). The per-beam sensitivities of the uncertainty
analysis ([`mode_sensitivity`](@ref), [`mode_uncertainty`](@ref)) deposit beams one by
one and cannot resolve E: they throw an `ArgumentError`; pass them the resolved bubble
[`linked_bubble`](@ref)`(model, beams)` instead.

*Derived* (confidence low): the exponents (E^0.2, Izumi et al., PoP 28, 022706 (2021),
seven shots; E^-0.16, Izumi et al., RSI 89, 10K111 (2018)) are published for picket
energy and quad splitting on one platform; their use for a fluence changed by
re-pointing at fixed energy has no published test. F̄ counts the overlap of the spots,
as Izumi 2021's split-quad shot does; the papers' energy over a fixed area
([`picket_irradiance`](@ref)'s `:average`) does not, and parts from F̄ on uneven powers
(research notes, "The swing and the bubble without new data (2026-11-18)" §5, and "The
published picket fluence, the trough and the inner-beam corridor (2026-11-18)" §5,
which measures the link read through the fixed area instead). Not a facility
criterion; no default reads it.
"""
struct PicketLinkedBubble{S<:Real,F} <: LaserPlasmaModel
    bubble::BubbleAbsorption{S}
    operator::DepositionOperator
    beams::Vector{Beam{Float64}}        # the operator's beams: the reference's in `groups`
    index::Dict{String,Int}             # their names → their operator index
    power::F                            # picket power per beam (TW)
    groups::Vector{Symbol}
    elements::Vector{Int}
    reference::PicketFluence{Float64}   # the reference's map and summaries
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

# The picket fluence of `beams` on the model's map: each beam of the model's groups at its
# picket power on the operator's beam of its name, the operator's other beams at zero.
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

# A vector of beams of mixed number types as one type (the promotion of the beams' own,
# abstract if one of them is).
_number_type(::Beam{T}) where {T} = T
function _common_type(beams)
    T = mapreduce(_number_type, promote_type, beams)
    return Beam{T}[Beam{T}(b.name, b.port, SVector{3,T}(b.aim), T(b.power), T(b.defocus), _retype(b.spot, T),
                           b.fnumber, b.group, b.quad, b.wavelength, T(b.dlambda)) for b in beams]
end

"""
    linked_bubble(m::PicketLinkedBubble, beams) -> BubbleAbsorption

The bubble of a [`PicketLinkedBubble`](@ref) for `beams`: `remake(m.bubble;
picket_ratio = E)` with E from the beams' picket fluence relative to the reference's
(`m.bubble` itself at the reference). In the beams' number type.
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

What the gold bubble model predicts for a re-pointed `design` (beams) against a
`reference` pointing, both mapped by `model` ([`PicketLinkedBubble`](@ref): its picket
powers, target and bubble; usually `reference` is the pointing the model was built on):

- `peak_ratio`: the maximum picket fluence (equivalently irradiance) over the model's
  elements, design over reference ([`picket_fluence`](@ref)'s `peak`);
- `fluence_ratio`: the energy-weighted mean picket fluence, design over reference;
- `picket_ratio`, `reference_picket_ratio`: E of each ([`linked_bubble`](@ref));
- `speed_ratio`: the bubble's tip speed, design over reference, `fluence_ratio^0.2`;
- `collision`, `reference_collision`: the ablator–wall collision times (ns);
- `bubble`, `reference_bubble`: the two [`BubbleAbsorption`](@ref)s; `fluence`,
  `reference_fluence`: the two [`PicketFluence`](@ref) maps;
- `history`: `nothing`, or, with `history`, the capsule drive's P2/P0 through time with
  the bubble linked.

`history(beams, bubble)` runs a drive history (the result of [`drive_history`](@ref),
with `lmax ≥ 2`) of `beams` with the laser–plasma chain built on `bubble`, e.g.
`(bs, bub) -> drive_history(model, τ -> with_pulses(bs, pulses, τ), times; lpi = (cbet, bub))`.
It runs for the design with its own bubble and for the reference with its own
(`reference_history`, a result of the same call for the reference, skips that run),
and with `fixed = true` also for the design with the reference's bubble, which splits
the change into the pointing's effect through the deposition and through the bubble.
`history` then holds `time`, `p2`, `reference_p2` and `fixed_p2` (or `nothing`); per
run (`design`, `reference`, `fixed`) the swing's `minimum` of P2/P0 over the `swing`
window (ns), the time it is reached (`at`) and the upward zero `crossing` after it
(linear between steps; `NaN` if none); `shift`, the design's minus the reference's;
and with `fixed`, `through_bubble` (design − fixed) and `through_deposition` (fixed −
reference).

*Derived* (confidence low): the bubble's response to the pointing rests on the
published exponents (tip speed ∝ E^0.2, Izumi 2021; collision ∝ E^-0.16, Izumi 2018)
applied to a fluence changed by re-pointing, which has no published test
([`PicketLinkedBubble`](@ref)); ratios between pointings on one mesh converge to < 1 %
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

# The swing of a P2/P0 history: its minimum over the window, when, and the first upward
# zero crossing after it (linear between the steps; NaN if none).
function _swing(t::Vector{Float64}, p::Vector{Float64}, window)
    ks = [k for k in eachindex(t) if window[1] <= t[k] <= window[2] && isfinite(p[k])]
    isempty(ks) && return (minimum = NaN, at = NaN, crossing = NaN)    # no step in the window
    km = ks[argmin(p[ks])]
    kz = findfirst(k -> k > km && p[k-1] < 0 <= p[k], eachindex(p))
    tz = kz === nothing ? NaN : t[kz-1] - p[kz-1] * (t[kz] - t[kz-1]) / (p[kz] - p[kz-1])
    return (minimum = p[km], at = t[km], crossing = tz)
end
