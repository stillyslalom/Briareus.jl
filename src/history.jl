# This file models the time-dependent drive with pulse shapes and quasi-static radiosity
# histories.

"""
    PulseShape(times, powers)

Piecewise-linear power history (ns, TW per beam), zero outside the table. Callable:
`pulse(t)`. Two pulse shapes are equal (`==`, `isequal`, `hash`) when their times and
powers are.
"""
struct PulseShape
    t::Vector{Float64}
    p::Vector{Float64}
    function PulseShape(t, p)
        length(t) == length(p) || throw(ArgumentError("times and powers differ in length"))
        issorted(t) || throw(ArgumentError("times must be sorted"))
        new(collect(Float64, t), collect(Float64, p))
    end
end

Base.:(==)(a::PulseShape, b::PulseShape) = a.t == b.t && a.p == b.p
Base.isequal(a::PulseShape, b::PulseShape) = isequal(a.t, b.t) && isequal(a.p, b.p)
Base.hash(s::PulseShape, h::UInt) = hash(s.p, hash(s.t, hash(:PulseShape, h)))

(s::PulseShape)(t) = (t < s.t[1] || t > s.t[end]) ? 0.0 : _interp(s.t, s.p, t)

# The trapezoid rule for ∫ y dt through the samples (t, y), of which there are at least two.
_trapz(t, y) = sum((t[k+1] - t[k]) * (y[k] + y[k+1]) / 2 for k in 1:length(t)-1)

"Laser energy (kJ = TW·ns) of one beam following the pulse."
pulse_energy(s::PulseShape) = length(s.t) < 2 ? 0.0 : _trapz(s.t, s.p)

"""
    with_pulses(beams, pulses, t; by = :group) -> Vector{Beam}

Beams with power set to `pulses[key](t)`. `key` is the `group` of each beam (pulses
given as a `Dict` or NamedTuple keyed by group), its `quad` as a `Symbol` (`by = :quad`),
its `name` as a `Symbol` (`by = :beam`) or `by(beam)` for a function. `pulses` may also
be a single `PulseShape` for all beams. Every pulse reader and writer keys by `Symbol`:
[`read_pulse_table`](@ref) by column header, [`cone_fraction_pulses`](@ref) by group and
[`read_visrad_powers`](@ref) by beam name (`by = :beam`).
"""
function with_pulses(beams, pulses, t; by = :group)
    key = _beam_key(by)
    return [remake(b; power = _beam_pulse(pulses, key, b)(t)) for b in beams]
end

_beam_key(by) = by === :group ? (b -> b.group) : by === :quad ? (b -> Symbol(b.quad)) :
                by === :beam ? (b -> Symbol(b.name)) : by

# The pulse of beam `b`: one `PulseShape` for all, or `pulses[key(b)]`.
_beam_pulse(p::PulseShape, key, b) = p
function _beam_pulse(p, key, b)
    k = key(b)
    haskey(p, k) || throw(KeyError(k))
    return p[k]
end

"""
    drive_history(model, op, beams_at, times; observe = :capsule, lmax = 6,
                  wall = (:wall,), cavities = (), lpi = nothing, holes = nothing, ylm = false,
                  observed_flux = false, transition = nothing, flux_change = nothing, min_step = 1e-3,
                  tol = 1e-4, keep_albedo = false, record = nothing,
                  wall_history = true, wall_release = true, mband = MBandSplit(),
                  spectral = nothing, unconverted = nothing, colors = nif_unconverted(),
                  unconverted_kw = (;), prepulse = nothing, by = :group,
                  per_component = false) -> NamedTuple

Quasi-static drive history. At each time, `beams_at(t)` is deposited with `op` after the
laser–plasma model `lpi`, evaluated at that time. The radiosity is solved at that time,
with time- and temperature-dependent albedos, and the following records are kept:

- `time` — the times solved (ns): `times`, plus any steps inserted by `transition` or
  `flux_change`
- `delivered` — laser power (TW)
- `flux`, `Tr` — mean incident flux on `observe` (TW/mm²) and its radiation
  temperature (eV)
- `moments` — `(lmax+1) × nt` Legendre moments of the observed flux (`aₗ/a₀` rows
  are `moments[l+1, :] ./ moments[1, :]`)
- `ylm` — with `ylm = true`, the `(lmax+1)² × nt` real spherical-harmonic
  coefficients of the observed flux (`ylm_weights`, in `mode_index` order; columns
  work with `mode_amplitude`), otherwise `nothing`
- `observed_flux` — with `observed_flux = true`, the `n × nt` incident flux (TW/mm²) on
  each of the `n` observed elements (in the order of `elements(target, observe...)`),
  which is the state that the retune of a flat face holds ([`FaceFlux`](@ref));
  otherwise `nothing`
- `wall_Tr`, `wall_albedo` — area-weighted mean incident-flux temperature and albedo
  of the `wall` components
- `cavity_Tr` — per cavity of `cavities` (pairs `name => wall names`, e.g.
  `[:upper => (:upper_wall, :upper_cap_top), :lower => (:lower_wall, :lower_cap_bottom)]`
  for two halfraums), the temperature of the area-weighted mean incident flux on the
  wall components of that cavity ([`wall_temperature`](@ref) of the solution of each
  step), keyed by cavity name; empty without `cavities` (the default). The walls of a
  cavity are named one by one, since grouping by part would also include the entrance
  holes of a hohlraum.
- `escaped` — power (TW) through each aperture, and `transmitted` through each
  thin-foil (`BurnthroughAlbedo`) and closing-gap (`GapClosure`) component, as in
  `energy_balance` (by default the fraction `s T` of the flux on each burned-through
  element, and `1 - s` on a gap a fraction `s` closed)
- `window_loss` — power lost through diagnostic windows: `transmitted` plus `escaped`
  through every aperture except the entrance `holes`. By default the holes are the
  apertures not on a `CylinderSurface`, which are the LEHs of a hohlraum.
- `front` — per thin-foil component, the area-mean heat-front areal mass as a fraction
  of the foil areal mass (`≥ 1` once burned through), from its flux history even when
  the model has an explicit `burnthrough`
- `burnthrough` — per thin-foil component, the burnthrough time of each element (ns;
  `Inf` if it has not burned through)
- `albedo` — with `keep_albedo = true`, the converged albedo per element at each step;
  otherwise `nothing`. This is the state about which a flux response is linearized, for
  example for `retune_drive`. Where a wall re-emits stored energy, it is the effective
  albedo that includes the emission, so that `B = ηP/A + αH`. Alongside it,
  `transmission` is the fraction of the incident flux per element transmitted out of the
  enclosure at each step (from the `RadiationSolution`), otherwise `nothing`.
  [`history_state`](@ref) interpolates both to a time of the pulse.
- `setup` — the power of `beams_at(t)` (TW, before `lpi`); `laser` — absorbed laser
  power per component (TW, back side included); `source` — x-ray source power per
  component (Σ ηP); `absorbed` — x-ray power absorbed per surface component (and per
  closing lip), as in `energy_balance`. `absorbed` is net of `released`, the power that
  wall components re-emit from stored energy (see below). `lost` — x-ray power lost to
  open space. `energy_budget` integrates these records. As in `energy_balance`,
  `escaped`, `transmitted`, `laser`, `source` and `absorbed` are keyed by part, and the
  pieces of a part split by `curved_components` are summed into it. `front` and
  `burnthrough` stay per component. With `per_component = true`, all of these records,
  `released` and the laser power of `mband` stay per component. This separates one piece
  of a part from the others, for example the cap of a hohlraum observed while its wall
  heats.
- `lip_radius` — per computed entrance-hole lip (`LipClosure`, `hohlraum(; leh_closure =
  :computed)`), the radius of its plasma edge at each step (mm)
- `record` — with `record = f`, the vector of `f(sol)` for the `RadiationSolution` of
  each step (e.g. `sol -> dante(dante1(), sol)` for a synthetic Dante, or a
  `spectral_solve`), otherwise `nothing`.
- `mband_moments`, `mband` — the `(lmax+1) × nt` Legendre moments of the observed flux
  at photon energies above the lower M-band edge, and the laser-produced M-band power
  per component (TW, keyed by part like `laser`). The moments are those of
  `band_flux(s, mband.emission.band[1])` for the spectral solution `s` of the step, which
  includes laser-produced and recirculating M-band and the thermal tail. Both are
  `nothing` with `mband = nothing`, when the target has no surface of the materials of
  the model, or, for a per-group model, when no beam of `beams_at(times[1])` is in its
  groups.
- `unconverted` — with `unconverted`, `(; delivered, skipped)`: per color, the
  unconverted power delivered at each step (TW, before absorption), and the beams whose
  unconverted light cannot be placed (see [`unconverted_deposition`](@ref)); otherwise
  `nothing`.

**Unconverted light.** `unconverted` gives the unconverted power per color of `colors`
([`UnconvertedColor`](@ref)s, the NIF 1ω and 2ω by default) as a fraction of the 3ω
power of each beam. The fraction is a number or a function of time (ns), keyed by color
name, for example `(omega1 = 0.05, omega2 = 0.01)` or
`(omega1 = PulseShape(ts, fs), omega2 = 0.0)`, as for the `ratio` of
[`unconverted_history`](@ref). There is no default, since the fraction is an input. The
fraction of a color applies to every beam, or each beam has a separate fraction. A
separate fraction is a vector with one entry per beam of `beams_at(τ)`, in its order,
which the pointing fixes, or it is keyed by `by` as `prepulse` is. A beam without an
entry is a `KeyError`, as in `with_pulses`. An example for 1ω on the outer cones only is
`(omega1 = Dict(:outer44 => 0.05, :outer50 => 0.05, :inner23 => 0.0, :inner30 => 0.0),
omega2 = 0.01)`. At each step the absorbed powers of the colors are those of
[`unconverted_deposition`](@ref) for the beams of `beams_at(τ)`, at the fractions of each
beam at that step. The beam powers are those before `lpi`, since unconverted light does
not take part in the 3ω laser–plasma model. The step is solved with these powers as the
`unconverted` argument of `solve`, and each color is converted at its absorbed flux.

`prepulse` adds light not proportional to the 3ω pulse, as the `prepulse` of
`unconverted_history` does. Per color it is optional, and it is a power per beam (TW) as
a function of time: one function for all beams, or keyed by `by` (`:group`, `:quad`,
`:beam` or a function, as in [`with_pulses`](@ref)). The power of a beam in color `c` at
step τ is then `r_c(τ) · P₃ω(τ) + prepulse[c](τ)`, with `r_c` its fraction in color `c`
and `P₃ω` the power of the beam in `beams_at(τ)`. A prepulse without unconverted
fractions uses zero fractions (`unconverted = (omega1 = 0.0, omega2 = 0.0)`). The
history is solved only at `times`, so a prepulse reaches the drive through the steps in
its window. A step where only the prepulse fires (every 3ω power zero, also before the
main pulse at negative times) is solved with only the prepulse light. Its x-rays heat
the wall, and the memory of the wall (below) carries them into the main pulse. Between a
prepulse and the main pulse, a step where nothing fires lets the wall cool and re-emit
the energy it stored. Without such a step, the wall-energy integral runs across the gap
at the interpolated temperature. `unconverted_kw` passes the `offsets`,
`lens_distance`, `coarse` and `resolution` of `unconverted_deposition`.

The light is linear in the power of each beam, so the history deposits it once per set
of beams whose powers are equal at every step and whose fractions and prepulses are the
same (equal numbers or the same function). For beams on shared pulses this is once per
(pulse, fraction, prepulse) schedule, as in `unconverted_history`. A set whose light is
zero is never deposited, and the deposited patterns are scaled by the powers in time.
The pointing must not change in time. `laser` then includes the absorbed power of each
color, and `setup` and `delivered` its delivered power, so `energy_budget` balances. The
M-band is still that of the 3ω light. Without `unconverted` the history is identical to
the 3ω history.

**M-band.** `mband` is an [`MBandSplit`](@ref), which contains the M-band emission model
and the spectral split with which it is solved. The default is `MBandSplit()`, with
`nif_mband()` (the per-cone conversions of Dewald 2020 on gold and Au-lined uranium) in
the calibrated two-group split. The laser-produced M-band of each step then comes from
the effective beams of that step after `lpi` (CBET, backscatter, the bubble), as
`mband_power(mband.emission, op, beams; lpi, time)` gives it. It is the share of each
beam group (cone) of every element, split off in the deposition pass of that step, so
the M-band follows the cone balance in time. An intensity-form
`MBandSplit(MBandEmission())` converts the deposition of the step directly. The gray
solution of the step is then split by `spectral_solve` with the `edges`, `thermalize`
and group albedos of the split, and `spectral` contains the other keywords of
`spectral_solve` (`spot`, `tol`, ...). The M-band is a diagnostic, and the gray history
(every other record) is identical with or without it. It costs about 6% of a small
history (600 elements), and less than the run-to-run spread on N221204 (research notes,
"Per-cone M-band in the drive history (2026-10-27)"). `mband = nothing` turns it off.

The temperature-dependent albedos of each step are iterated to `tol` (see `solve`),
starting from the temperatures of the previous step.

**Wall-loss history.** With `wall_history = true` (the default), each element whose
albedo is a [`WallLossAlbedo`](@ref) loses energy as its temperature history requires.
The loss is not that of a wall held at the present temperature since `t = 0`. The
constant-T wall energy of Lindl, `E_w ∝ T^a t^b` (a = 4 + p = 3.3, b = 1 + q = 0.62 for
the defaults), is extended to a varying `T(t)` as `E_w ∝ [∫ (T/100 eV)^(a/b) dt]^b`.
This is the form of the drive integral of [`HeatFront`](@ref), exact at constant T. Its
rate is the constant-T loss at the effective time

    t_eff = ∫₀ᵗ (T(t′)/T(t))^(a/b) dt′,

integrated trapezoidally between the steps from the incident-flux temperature of each
element. The temperature is taken as constant from `t = 0` to `times[1]`, and no such
interval is added before a first step at or before `t = 0`, such as that of a prepulse.
A wall heated in a low foot holds little energy at the peak temperature. Its loss early
in the peak is therefore several times the constant-T value at `t` (`t_eff` ≈ the time
since the rise), and the drive rises through a flat peak.

A thin-foil element ([`BurnthroughAlbedo`](@ref)) whose `before` is a `WallLossAlbedo`
follows the same history for as long as `before` enters its albedo. While intact it is a
wall element with that loss, so an intact window equals the wall around it. Across its
transition and after it, `α_before` comes from the history, and the stored energy that
it re-emits (below) is that of the wall element scaled by `α/α_before` (`1 − s f`). Its
heat front does not depend on this. `wall_history = false` evaluates the albedo at `t`,
for walls and foils alike. This was the behavior before 2026-10-24, and for foils before
2026-11-19. See the research notes, "N221204 late-time drive (2026-10-24)" and "Intact
windows follow the wall's history (2026-11-19)".

**Wall re-emission.** The self-similar loss never gives energy back: when `T` falls,
`t_eff` grows and the loss only shrinks. With `wall_release = true` (the default, with
`wall_history`), a falling drive makes the wall re-emit stored energy. The loss above is
taken as the heating response to the highest temperature of the element so far,
`T_h = max T`. A drop below it gives the deficit `D = φ(T_h) − φ(T)` in the wall-energy
potential `φ = (T/100 eV)^a`, and the response to the deficit is the same self-similar
kernel, superposed. This follows the principle of Duhamel for the diffusion in the wall,
which is linear in `φ` (*engineering*):

    F = L(T_h, t_eff) − f_r b C ∫ (t − t′)^(b−1) dD(t′),    C = L t_eff^(1−b) / (b φ(T_h)),

is the net flux into the wall. `C`, the coefficient of the loss, is constant for the
Lindl loss of gold. `f_r` is the `release` of `WallLossAlbedo`, 1.0, calibrated against
1-D radiation hydrodynamics of an ablating gold wall (`examples/late_time_wall.jl`). A
constant-T wall loses `C b φ t^(b−1)`. A drop `D` at `t_d` gives back
`f_r C D (t − t_d)^b` by `t`, quickly at first from the hot surface layer and then as
the `t^b` of the heat wave. The effective albedo `1 − F/H` exceeds 1 where the wall gives
back more than it takes, and `solve` includes the release as an emission added to the
source of the element. While no element temperature falls below its running maximum,
the history is identical to that of the heating loss. At a step where nothing fires,
the flux of a heated wall is its re-emission, and the albedo iteration of `solve` starts
it at half its running maximum. An element without flux would otherwise take
`reference_temperature`. A wall heated again after cooling takes back the energy of its
surface layer, but never more than the flux on it: its radiosity is kept ≥ 0 (albedo 0,
no emission, *engineering*). `released` records the re-emitted power per component, and
`energy_budget` integrates it; `absorbed` is net of it. See the research notes, "Wall
re-emission and the capsule radius (2026-10-24)" for the closure and "Capsule
absorption, conversion and the ablating-wall release (2026-10-25)" for `release` 1.0.

**Computed lip closure.** Each entrance-hole lip with a [`LipClosure`](@ref) albedo has
a closure depth `d`, open at `times[1]`. The depth is advanced over each step from the
trapezoidal mean of the free speed `v(T)` of the lip and of the laser hold. `T` is the
area-mean incident-flux temperature of the lip at the start of the step and, corrected,
at its end. The hold is the power of the beams in `beams_at`, after `lpi`, whose ports
lie on the side of the hole, over the initial area of the hole. The update is exact for
the relaxation of the hold. The step is solved with each lip ring closed outside the
edge `r₀ − d`, and solved again until the corrected depth changes by less than 1e-5 mm.

Each thin-foil element without an explicit `burnthrough` time tracks a separate heat
front. The front is the drive integral of `heat_front` over the incident-flux
temperature history of that element, trapezoidal between the steps from `times[1]`, so
its burnthrough is consistent with the flux that the windowed hohlraum puts on it. A
step is solved with the front extended from the previous step at constant temperature.
It is solved again until the temperatures of that step no longer move a burnthrough
that affects its albedo.

A burnthrough transition takes the `width` of the foil (0.4 ns by default), and a step
much longer than that jumps across it. Two optional criteria insert steps between
`times`:

- `transition` — the transition of no thin-foil element (the smooth step `s` of
  `BurnthroughAlbedo`, 0 intact to 1 burned through) and of no closing gap
  (`GapClosure`, 0 open to 1 closed) advances by more than this in one step.
  Steps are planned from the known or predicted burnthrough time of each element. A
  step is halved if the temperatures of that step move a burnthrough so that the advance
  exceeds 1.5 `transition`.
- `flux_change` — the observed flux relative to its mean, `H/⟨H⟩` on each element of
  `observe`, changes by no more than this in one step. A step that exceeds it is halved
  and solved again.

Neither criterion shortens a step below `min_step` (ns). A planned step that would leave
less than `min_step` before the next of `times` is extended to that time instead.
Without the criteria (the default) the steps are exactly `times`. With `wall_release`,
the albedo of a cooling wall element in a step of length Δ falls as `1 − c Δ^(b−1)`,
the implicit part of the release kernel. A step of 1e-5 ns gave albedos near −27, which
the unpreconditioned GMRES of a large model does not solve, so steps much shorter than
`min_step` are best left out of `times` as well.
"""
drive_history(m::RadiationModel, op::DepositionOperator, beams_at, times; kwargs...) =
    _drive_history(m.target, τ -> m, τ -> op, beams_at, times; kwargs...)

# Computes the history on the models `model_at(τ)` and deposition operators
# `operator_at(τ)`, whose elements are those of `t` (a fixed target, or the keyframes of a
# moving one). The function is compiled and inferred once for every closure, laser–plasma
# model and recorder. Before, the `beams_at` closure of each call site recompiled the
# whole march, which took about 2 s. Their values reach the solves by dynamic dispatch, at
# no measurable run-time cost, and the march is in the precompile workload of the package
# (research notes, "Compile latency").
# The options of a history that are `nothing` or a value, in one concrete type, so that the
# march is compiled once for every combination of them.
struct _HistoryOptions
    transition::Union{Nothing,Float64}
    flux_change::Union{Nothing,Float64}
    emission::Union{Nothing,MBandEmission}
    ulight::Union{Nothing,_UnconvertedLight}
    min_step::Float64
    tol::Float64
    lmax::Int
    ylm::Bool
    observed_flux::Bool
    keep_albedo::Bool
    wall_release::Bool
    per_component::Bool
    has_record::Bool
end

# The state of the march at the last step taken: the drive integral `I`, burnthrough time
# `tb` and integrand `rate` of each thin-foil element, the wall-energy integral `Iw` and
# integrand `ratew` of each wall-loss element, the highest temperature `Tmax`, deficit
# `Dw` and `cooled` flag of the wall release with the deficits and times of the steps
# taken, and the depth, temperature and hole intensity of each computed lip.
mutable struct _MarchState
    I::Vector{Float64}
    tb::Vector{Float64}
    rate::Vector{Float64}
    Iw::Vector{Float64}
    ratew::Vector{Float64}
    Tmax::Vector{Float64}
    Dw::Vector{Float64}
    deficits::Vector{Vector{Float64}}
    dtimes::Vector{Float64}
    cooled::BitVector
    ldepth::Vector{Float64}
    lT::Vector{Float64}
    lI::Vector{Float64}
end

Base.@nospecializeinfer function _drive_history(
                        t::Target, @nospecialize(model_at), @nospecialize(operator_at),
                        @nospecialize(beams_at), @nospecialize(times);
                        observe = :capsule, lmax::Integer = 6, wall = (:wall,),
                        @nospecialize(cavities = ()), @nospecialize(lpi = nothing), holes = nothing, ylm::Bool = false,
                        observed_flux::Bool = false, transition = nothing, flux_change = nothing, min_step::Real = 1e-3,
                        tol = 1e-4, keep_albedo::Bool = false, @nospecialize(record = nothing),
                        wall_history::Bool = true, wall_release::Bool = true,
                        mband = MBandSplit(), @nospecialize(spectral = nothing),
                        @nospecialize(unconverted = nothing), colors = nif_unconverted(),
                        @nospecialize(unconverted_kw = (;)), @nospecialize(prepulse = nothing),
                        @nospecialize(by = :group), per_component::Bool = false)
    cx, st, rec, spectral = _init_history_state(
        t, beams_at, times; observe, lmax, wall, cavities, holes, ylm, observed_flux, transition, flux_change,
        min_step, tol, keep_albedo, record, wall_history, wall_release, mband, spectral, unconverted, colors,
        unconverted_kw, prepulse, by, per_component)
    _march_history!(cx, st, rec, model_at, operator_at, beams_at, times, lpi, record, spectral)
    return _assemble_history_output(cx, st, rec)
end

# Checks the options of a history on target `t` and returns what the march works with: the
# context `cx` (the element sets and constants of the target with the options), the
# state `st` of the march before the first step, the empty records `rec` with one
# container per quantity, and the `spectral_solve` keywords of the M-band.
Base.@nospecializeinfer function _init_history_state(
                        t::Target, @nospecialize(beams_at), @nospecialize(times);
                        observe, lmax::Integer, wall, @nospecialize(cavities), holes, ylm::Bool,
                        observed_flux::Bool, transition, flux_change, min_step::Real, tol, keep_albedo::Bool,
                        @nospecialize(record), wall_history::Bool, wall_release::Bool, mband,
                        @nospecialize(spectral), @nospecialize(unconverted), colors,
                        @nospecialize(unconverted_kw), @nospecialize(prepulse), @nospecialize(by),
                        per_component::Bool)
    # fractions and depositions of the unconverted light (`nothing` without it)
    ulight = unconverted === nothing ? nothing :
             _UnconvertedLight(unconverted, prepulse, by, colors, unconverted_kw)
    unconverted === nothing && !isempty(unconverted_kw) &&
        throw(ArgumentError("unconverted_kw sets the unconverted light's deposition; without `unconverted` it has no effect"))
    unconverted === nothing && prepulse !== nothing &&
        throw(ArgumentError("a prepulse needs `unconverted`, the colors' fractions of the 3ω power (0 for none)"))
    # the M-band model of each step (`nothing` when off) and its `spectral_solve` keywords
    emission = nothing
    if mband === nothing
        spectral === nothing ||
            throw(ArgumentError("`spectral` sets the M-band's spectral_solve; with mband = nothing it has no effect"))
    else
        mband isa MBandSplit ||
            throw(ArgumentError("mband must be an MBandSplit (e.g. MBandSplit(nif_mband())) or nothing"))
        spectral = spectral === nothing ? (;) : NamedTuple(spectral)
        split_keys = filter(in((:edges, :thermalize, :albedo, :mband, :mband_power)), collect(keys(spectral)))
        isempty(split_keys) || throw(ArgumentError("`spectral` takes spectral_solve's other keywords; " *
                                                   join(split_keys, ", ") * " are the MBandSplit's (`mband`)"))
        albedo = mband.albedo isa Float64 ?
                 mband_albedo(t, mband.emission, mband.edges; albedo = mband.albedo) : mband.albedo
        spectral = (; spectral..., edges = mband.edges, thermalize = mband.thermalize, albedo)
        # The model is on only where there is M-band to make, on a target with an emitting
        # surface lit by a beam group that the model converts (any beam, for the
        # intensity form).
        em = mband.emission
        emits = any(e -> _emits_mband(em, t, e), 1:nelements(t))
        lit = em.groups === nothing || isempty(times) ||
              any(b -> haskey(em.groups, b.group), beams_at(first(times)))
        emits && lit && (emission = em)
    end
    transition === nothing || transition > 0 || throw(ArgumentError("transition must be positive"))
    flux_change === nothing || flux_change > 0 || throw(ArgumentError("flux_change must be positive"))
    min_step > 0 || throw(ArgumentError("min_step must be positive"))
    issorted(times) || throw(ArgumentError("times must be sorted"))
    obs = elements(t, (observe isa Symbol ? (observe,) : observe)...)
    we = elements(t, wall...)
    # the wall elements of each cavity (`cavities` maps names to wall names)
    cavity_elems = _cavity_elements(t, cavities)
    weights = _TargetCache(tg -> ylm_weights(tg, obs; lmax))
    N = nelements(t)
    # thin-foil components and the elements whose burnthrough is tracked
    foils = [ci for (ci, c) in enumerate(t.components)
             if c.kind === :surface && c.material.albedo isa BurnthroughAlbedo]
    foil_elems = [e for ci in foils for e in t.comp_elems[ci]]
    # components that transmit part of their flux out: thin foils and closing gaps
    transmitters = [ci for (ci, c) in enumerate(t.components)
                    if c.kind === :surface && _transmits(c.material.albedo)]
    model(e) = element_material(t, e).albedo
    tracked = [e for e in foil_elems if model(e).burnthrough === nothing]
    # transitions a step can resolve (closing entrance-hole lips included)
    smooth = [e for (ci, c) in enumerate(t.components)
              if _transmits(c.material.albedo) && !(c.material.albedo isa LipClosure)
              for e in t.comp_elems[ci] if model(e).width > 0]
    # entrance-hole lips whose closure is computed (`LipClosure`): their geometry, and
    # the closure depth, lip temperature and hole intensity at the last step taken
    lips = _lip_geometry(t)
    ldepth, lT, lI = zeros(length(lips)), zeros(length(lips)), zeros(length(lips))
    need = zeros(N)
    for e in foil_elems
        need[e] = _front_integral(model(e).front, areal_mass(model(e)))
    end
    tb = fill(NaN, N)              # burnthrough times (NaN if not burned)
    for e in foil_elems
        model(e).burnthrough === nothing || (tb[e] = model(e).burnthrough)
    end
    I = zeros(N)                   # drive integral of each thin-foil element
    rate = zeros(N)                # its integrand at the previous time
    # Wall-loss elements whose loss follows the temperature history of the element are
    # walls and thin foils whose `before` is a wall loss. `wl` returns the
    # `WallLossAlbedo` followed.
    wl(e) = (m = model(e); m isa BurnthroughAlbedo ? m.before : m)
    lossy = wall_history ?
            [e for e in 1:N if !is_aperture(t, e) && wl(e) isa WallLossAlbedo] : Int[]
    Iw, ratew = zeros(N), zeros(N) # their wall-energy integral and integrand
    # For re-emission (`wall_release`), the highest temperature of each element, the
    # deficit φ(T_max) − φ(T) at every step taken, and the times of those steps
    Tmax, Dw = zeros(N), zeros(N)
    # (the deficits of the `lossy` elements only, in their order: `lossy_index`)
    deficits, dtimes = Vector{Float64}[], Float64[]
    lossy_index = zeros(Int, N)
    lossy_index[lossy] .= 1:length(lossy)
    cooled = falses(N)             # elements whose deficit has ever been nonzero
    holes = holes === nothing ?
            [c.name for c in t.components if c.kind === :aperture && !(c.surface isa CylinderSurface)] :
            collect(Symbol, holes)
    apertures = [c.name for c in t.components if c.kind === :aperture]
    partof = Dict(c.name => c.part for c in t.components)
    # records, one entry per step taken
    time = Float64[]
    delivered, flux, Tr = Float64[], Float64[], Float64[]
    wall_Tr, wall_albedo = Float64[], Float64[]
    cavity_Tr = Dict{Symbol,Vector{Float64}}(first(c) => Float64[] for c in cavity_elems)
    moments = Vector{Float64}[]
    clm = Vector{Float64}[]
    Hobs = Vector{Float64}[]
    escaped = Dict{Symbol,Vector{Float64}}(n => Float64[] for n in apertures)
    transmitted = Dict{Symbol,Vector{Float64}}(t.components[ci].name => Float64[] for ci in transmitters)
    front = Dict{Symbol,Vector{Float64}}(t.components[ci].name => Float64[] for ci in foils)
    window_loss = Float64[]
    albedo_steps = Vector{Float64}[]
    transmission_steps = Vector{Float64}[]
    recorded = Any[]
    names = [c.name for c in t.components]
    setup, lost = Float64[], Float64[]
    laser = Dict{Symbol,Vector{Float64}}(n => Float64[] for n in names)
    source = Dict{Symbol,Vector{Float64}}(n => Float64[] for n in names)
    released = Dict{Symbol,Vector{Float64}}(n => Float64[] for n in names)
    absorbed = Dict{Symbol,Vector{Float64}}(c.name => Float64[] for c in t.components
                    if c.kind !== :aperture || _transmits(c.material.albedo))
    lip_radius = Dict{Symbol,Vector{Float64}}(L.name => Float64[] for L in lips)
    mband_moments = Vector{Float64}[]
    mband_laser = Dict{Symbol,Vector{Float64}}(n => Float64[] for n in names)
    uc_delivered = Dict{Symbol,Vector{Float64}}(c.name => Float64[] for c in (ulight === nothing ? UnconvertedColor[] : ulight.colors))
    options = _HistoryOptions(transition === nothing ? nothing : Float64(transition),
                              flux_change === nothing ? nothing : Float64(flux_change), emission, ulight,
                              Float64(min_step), Float64(tol), Int(lmax), ylm, observed_flux, keep_albedo,
                              wall_release, per_component, record !== nothing)
    cx = (; t, obs, we, cavity_elems, weights, N, foils, foil_elems, transmitters, tracked, smooth,
     lips, need, lossy, lossy_index, holes, apertures, partof, names, options)
    st = _MarchState(I, tb, rate, Iw, ratew, Tmax, Dw, deficits, dtimes, cooled, ldepth, lT, lI)
    rec = (; time, delivered, flux, Tr, wall_Tr, wall_albedo, cavity_Tr, moments, clm, Hobs,
     escaped, transmitted, front, window_loss, albedo_steps, transmission_steps, recorded, setup,
     lost, laser, source, released, absorbed, lip_radius, mband_moments, mband_laser, uc_delivered)
    return cx, st, rec, spectral
end

# Takes the step from τ0 (`nothing` at the first time) to τ and returns its solution
# and the thin-foil state at τ. The state at τ0 is left untouched, so a step can be
# retried. `Il` is the laser intensity through the hole of each computed lip at τ,
# and `b0` the beams at τ before the laser–plasma model, which give the powers of the
# unconverted light.
Base.@nospecializeinfer function _solve_history_step!(cx, st::_MarchState, @nospecialize(model_at), τ0, τ, dep, guess,
                                                      Il, b0)
    (; t, obs, we, cavity_elems, weights, N, foils, foil_elems, transmitters, tracked, smooth,
     lips, need, lossy, lossy_index, holes, apertures, partof, names) = cx
    (; transition, flux_change, emission, ulight, min_step, tol, lmax, ylm, observed_flux,
     keep_albedo, wall_release, per_component) = cx.options
    model(e) = element_material(t, e).albedo
    wl(e) = (m = model(e); m isa BurnthroughAlbedo ? m.before : m)
    (; I, tb, rate, Iw, ratew, Tmax, Dw, deficits, dtimes, cooled, ldepth, lT, lI) = st
    Δ = τ0 === nothing ? 0.0 : τ - τ0
    # The predictor advances each front through the step at its last temperature.
    pred = copy(tb)
    for e in tracked
        isnan(tb[e]) && rate[e] > 0 && I[e] + Δ * rate[e] >= need[e] &&
            (pred[e] = τ0 + (need[e] - I[e]) / rate[e])
    end
    mτ = model_at(τ)
    uc = ulight === nothing ? nothing : _unconverted_at!(ulight, mτ.target, b0, τ)
    ucs = uc === nothing ? nothing : _with_colors(uc, ulight.colors)   # as `solve` takes it
    # The per-element state used by the albedos holds the burnthrough times and the
    # wall-loss history of the `lossy` elements (`_WallHistory`).
    hist = Dict(e => (τ0 === nothing ? _WallHistory(0.0, max(τ, 0.0)) :
                      !wall_release ? _WallHistory(Iw[e] + Δ * ratew[e] / 2, Δ / 2) :
                      _WallHistory(Iw[e] + Δ * ratew[e] / 2, Δ / 2, Tmax[e],
                                   cooled[e] ? _deficit_sum(wl(e), deficits, dtimes, lossy_index[e], τ) : 0.0,
                                   Dw[e], Δ)) for e in lossy)
    # closure depths of the lips at τ, from the lip temperatures `T1` at τ
    depths(T1) = [τ0 === nothing ? 0.0 :
                  clamp(_lip_advance(L.model, ldepth[k], Δ,
                                     (_lip_speed(L.model, lT[k]) + _lip_speed(L.model, T1[k])) / 2,
                                     (_lip_hold(L.model, lI[k]) + _lip_hold(L.model, Il[k])) / 2),
                        0.0, L.r0 - L.rin) for (k, L) in enumerate(lips)]
    ds = depths(lT)            # the predictor uses the lip temperatures of the last step
    # the wall histories per element, as `solve` takes them (`_StepState`)
    walls = Vector{Union{Nothing,_WallHistory}}(nothing, N)
    for e in lossy
        walls[e] = hist[e]
    end
    function state(bt)
        isempty(lossy) && isempty(lips) && return bt
        s = copy(bt)
        for (k, L) in enumerate(lips), (j, e) in enumerate(L.elems)
            s[e] = _ring_closed(L.ra[j], L.rb[j], L.r0 - ds[k])
        end
        return _StepState(s, walls)
    end
    sol = solve(mτ, dep; time = τ, burnthrough = state(pred), guess, tol, unconverted = ucs)
    # The corrector sets the lip depths from the lip temperatures of this step.
    for _ in 1:5
        (isempty(lips) || τ0 === nothing) && break
        dn = depths(_lip_temperatures(lips, sol))
        maximum(abs.(dn .- ds)) <= 1e-5 && break
        ds = dn
        sol = solve(mτ, dep; time = τ, burnthrough = state(pred), guess = sol, tol, unconverted = ucs)
    end
    # The corrector computes the trapezoidal integral with the temperatures of this
    # step, from the state at τ0 in `I0, tb0, rate0`. The `I, tb, rate` of the march
    # are reassigned after each step, so if captured here they would be boxed, and
    # possibly undefined to the report_call of JET.
    I0, tb0, rate0 = I, tb, rate
    function advance(sol)
        # With `local`, the per-step `Ik, tk, rk` below are not captured. A capture
        # would box them and could leave them undefined to the report_call of JET.
        local Ik, tk, rk = copy(I0), copy(tb0), copy(rate0)
        for e in foil_elems
            rk[e] = _front_rate(model(e).front, radiation_temperature(ForwardDiff.value(sol.H[e])))
            τ0 === nothing && continue
            dI = Δ * (rate0[e] + rk[e]) / 2
            if isnan(tk[e]) && Ik[e] + dI >= need[e]   # (explicit times are never NaN)
                tk[e] = τ0 + _linear_crossing(rate0[e], rk[e], Δ, need[e] - Ik[e])
            end
            Ik[e] += dI
        end
        return Ik, tk, rk
    end
    Ik, tk, rk = advance(sol)
    # The step is solved again while its temperatures move a burnthrough that affects
    # its albedo. This is a fixed point, and the transition is smooth, so it converges
    # in a few passes.
    used = pred
    for _ in 1:20
        moved(e) = !isequal(tk[e], used[e]) && min(_nan_inf(tk[e]), _nan_inf(used[e])) <= τ &&
                   !(abs(tk[e] - used[e]) <= 1e-12 * max(abs(τ), 1.0))
        any(moved, tracked) || break
        used = tk
        sol = solve(mτ, dep; time = τ, burnthrough = state(used), guess = sol, tol, unconverted = ucs)
        Ik, tk, rk = advance(sol)
    end
    # the wall-energy integrals at τ with the temperatures of this step (the heating
    # response at the highest temperature so far when the release is on)
    Iwk, rwk = copy(Iw), copy(ratew)
    for e in lossy
        T = radiation_temperature(ForwardDiff.value(sol.H[e]))
        wall_release && (T = max(T, Tmax[e]))
        rwk[e] = _loss_rate(wl(e), T)
        Iwk[e] = τ0 === nothing ? max(τ, 0.0) * rwk[e] : hist[e].A + Δ * rwk[e] / 2
    end
    return sol, Ik, tk, rk, Iwk, rwk, hist, ds, uc
end

# Marches the history through `times`: plans each step, solves it (`_solve_history_step!`),
# halves it while a transition or the flux shape moves too far, and records it
# (`_record_history_step!`).
Base.@nospecializeinfer function _march_history!(cx, st::_MarchState, rec, @nospecialize(model_at),
                                                 @nospecialize(operator_at), @nospecialize(beams_at),
                                                 @nospecialize(times), @nospecialize(lpi),
                                                 @nospecialize(record), @nospecialize(spectral))
    (; t, obs, we, cavity_elems, weights, N, foils, foil_elems, transmitters, tracked, smooth,
     lips, need, lossy, lossy_index, holes, apertures, partof, names) = cx
    (; transition, flux_change, emission, ulight, min_step, tol, lmax, ylm, observed_flux,
     keep_albedo, wall_release, per_component) = cx.options
    model(e) = element_material(t, e).albedo
    wl(e) = (m = model(e); m isa BurnthroughAlbedo ? m.before : m)
    # the largest advance of a smooth transition from τ0 to τ, given burnthrough times
    advance_of(τ0, τ, tbs) =
        maximum((_transition_step(model(e), τ, _nan_inf(tbs[e])) -
                 _transition_step(model(e), τ0, _nan_inf(tbs[e])) for e in smooth); init = 0.0)
    shape(sol, a0) = ForwardDiff.value.(sol.H[obs]) ./ ForwardDiff.value(a0)

    τ0 = nothing                   # time of the last step taken
    last = nothing                 # its solution, which starts the iteration of the next step
    shape0 = nothing               # its observed flux shape
    for τgoal in times
        while true
            (; I, tb, rate) = st       # the state at the last step taken
            τ = τgoal
            if τ0 !== nothing && transition !== nothing && !isempty(smooth)
                # plan the step from the burnthrough times known or predicted at τ0
                tc = copy(tb)
                for e in tracked
                    isnan(tb[e]) && rate[e] > 0 && (tc[e] = τ0 + (need[e] - I[e]) / rate[e])
                end
                if advance_of(τ0, τ, tc) > transition
                    lo, hi = τ0 + min_step, τ       # the longest step within the tolerance
                    if advance_of(τ0, lo, tc) <= transition
                        while hi - lo > 1e-3 * min_step
                            mid = (lo + hi) / 2
                            advance_of(τ0, mid, tc) > transition ? (hi = mid) : (lo = mid)
                        end
                    end
                    τ = min(lo, τgoal)
                    # A remainder shorter than `min_step` is added to this step. A sliver
                    # left before the goal is pure discretization, and the short-time
                    # kernel (Δ^(b-1)) of the wall release makes its albedo strongly
                    # negative.
                    τgoal - τ < min_step && (τ = τgoal)
                end
            end
            local sol, Ik, tk, rk, Iwk, rwk, hist, ds, Il, a, dep, P0, bl, msplit, uc
            while true
                b0 = beams_at(τ)
                P0 = sum(b -> ForwardDiff.value(b.power), b0; init = 0.0)
                lr = _resolve_lpi(lpi, b0)      # resolved once for the step
                bl = apply_lpi(lr, b0, τ)
                # for the per-group M-band, the cones are split in the same deposition pass
                msplit = emission === nothing || emission.groups === nothing ? nothing :
                         _mband_split(emission, t, bl)
                dep = _lpi_deposit(lr, operator_at(τ), bl, τ, msplit)
                Il = _lip_intensities(lips, bl)
                sol, Ik, tk, rk, Iwk, rwk, hist, ds, uc =
                    _solve_history_step!(cx, st, model_at, τ0, τ, dep, last, Il, b0)
                a = legendre_moments(sol.H, sol.model.target, obs; lmax)
                τ0 === nothing && break
                h = τ - τ0
                h / 2 >= min_step || break
                too_far = (transition !== nothing && advance_of(τ0, τ, tk) > 1.5 * transition) ||
                          (flux_change !== nothing &&
                           maximum(abs.(shape(sol, a[1]) .- shape0)) > flux_change)
                too_far || break
                τ = τ0 + h / 2
            end
            _record_history_step!(cx, st, rec, record, spectral, τ, sol, dep, hist, ds, Il, uc, P0, msplit, a,
                                  Ik, tk, rk, Iwk, rwk)
            τ0, last = τ, sol
            shape0 = shape(sol, a[1])
            τ == τgoal && break
        end
    end
    return nothing
end

# Commits the step solved at `τ` to the state of the march and appends its records: the
# thin-foil and wall states `Ik, tk, rk, Iwk, rwk`, the lip state (`ds`, `Il`), the wall
# release, and one entry per recorded quantity from the solution `sol`, the deposition
# `dep`, the unconverted light `uc`, the setup power `P0`, the M-band split `msplit` and
# the moments `a`.
Base.@nospecializeinfer function _record_history_step!(cx, st::_MarchState, rec, @nospecialize(record),
                                                       @nospecialize(spectral), τ, sol, dep, hist, ds, Il, uc, P0,
                                                       msplit, a, Ik, tk, rk, Iwk, rwk)
    (; t, obs, we, cavity_elems, weights, N, foils, foil_elems, transmitters, tracked, smooth,
     lips, need, lossy, lossy_index, holes, apertures, partof, names) = cx
    (; transition, flux_change, emission, ulight, min_step, tol, lmax, ylm, observed_flux,
     keep_albedo, wall_release, per_component) = cx.options
    model(e) = element_material(t, e).albedo
    wl(e) = (m = model(e); m isa BurnthroughAlbedo ? m.before : m)
    (; time, delivered, flux, Tr, wall_Tr, wall_albedo, cavity_Tr, moments, clm, Hobs, escaped,
     transmitted, front, window_loss, albedo_steps, transmission_steps, recorded, setup, lost,
     laser, source, released, absorbed, lip_radius, mband_moments, mband_laser, uc_delivered) = rec
    st.I, st.tb, st.rate = Ik, tk, rk
    st.Iw, st.ratew = Iwk, rwk
    (; I, tb, rate, Iw, ratew, Tmax, Dw, deficits, dtimes, cooled, ldepth, lT, lI) = st
    if !isempty(lips)
        ldepth .= ds
        lT .= _lip_temperatures(lips, sol)
        lI .= Il
        for (k, L) in enumerate(lips)
            push!(lip_radius[L.name], L.r0 - ds[k])
        end
    end
    # The release is the power re-emitted from stored energy per component. The
    # running maximum and deficit of each element at τ follow.
    rel = zeros(N)
    if wall_release
        for e in lossy
            T = radiation_temperature(ForwardDiff.value(sol.H[e]))
            rel[e] = _wall_split(wl(e), hist[e], T)[3] * sol.model.target.elem_area[e]
            # A thin foil re-emits the stored wall energy scaled with its albedo.
            model(e) isa BurnthroughAlbedo &&
                (rel[e] *= _foil_history_split(model(e), τ, T, _FoilHistory(tb[e], hist[e]))[5])
            Th = max(T, Tmax[e])
            Dw[e] = _wall_potential(wl(e), Th) - _wall_potential(wl(e), T)
            Tmax[e] = Th
            Dw[e] != 0 && (cooled[e] = true)
        end
        push!(deficits, Dw[lossy])
        push!(dtimes, τ)
    end
    push!(time, τ)
    # the unconverted light delivered and absorbed, per color (TW)
    udel = 0.0
    if uc !== nothing
        for c in ulight.colors
            push!(uc_delivered[c.name], uc[c.name].delivered)
            udel += uc[c.name].delivered
        end
    end
    push!(delivered, uc === nothing ? dep.delivered : dep.delivered + udel)
    push!(moments, a)
    ylm && push!(clm, weights(sol.model.target) * sol.H[obs])
    observed_flux && push!(Hobs, sol.H[obs])
    Aw = sol.model.target.elem_area[we]
    push!(flux, a[1])
    push!(Tr, radiation_temperature(a[1]))
    Hw = sum(sol.H[we] .* Aw) / sum(Aw)
    push!(wall_Tr, radiation_temperature(Hw))
    for (name, es) in cavity_elems
        push!(cavity_Tr[name], _wall_temperature(sol, es))
    end
    # The emission that a re-emitting wall adds to the source of the solution
    # (TW/mm²). The wall albedo is the effective one, reflected over incident flux.
    xe = wall_release && !isempty(lossy) ?
         ForwardDiff.value.(sol.source) .- ForwardDiff.value.(sol.converted) :
         nothing
    αw = xe === nothing || all(iszero, xe) ? sol.albedo[we] :
         sol.albedo[we] .+ xe[we] ./ ForwardDiff.value.(sol.H[we])
    push!(wall_albedo, sum(αw .* Aw) / sum(Aw))
    keep_albedo && push!(albedo_steps, xe === nothing || all(iszero, xe) ? copy(sol.albedo) :
                           sol.albedo .+ xe ./ max.(ForwardDiff.value.(sol.H), eps()))
    keep_albedo && push!(transmission_steps, copy(sol.transmission))
    record === nothing || push!(recorded, record(sol))
    if emission !== nothing
        # The effective beams of this step give the share of each cone in the
        # deposition, from only the beams of that cone, which equals
        # `mband_power(emission, op, bl; lpi, time)`. The intensity form converts
        # the deposition of the step directly.
        tm = sol.model.target
        PM = msplit === nothing ? mband_power(emission, tm, dep.element) :
             _mband_from_split(emission, tm, msplit)
        sm = spectral_solve(sol; spectral..., mband = emission, mband_power = PM)
        push!(mband_moments, legendre_moments(band_flux(sm, emission.band[1]), tm, obs; lmax))
        for (ci, c) in enumerate(tm.components)
            push!(mband_laser[c.name], sum(sm.mband[e] * tm.elem_area[e] for e in tm.comp_elems[ci]; init = 0.0))
        end
    end
    bal = energy_balance(sol; by = :component)
    loss = 0.0
    for n in apertures
        push!(escaped[n], bal.escaped[n])
        n in holes || partof[n] in holes || (loss += bal.escaped[n])
    end
    for ci in transmitters
        c = t.components[ci]
        push!(transmitted[c.name], bal.transmitted[c.name])
        loss += bal.transmitted[c.name]
    end
    # the target of this step, whose element areas differ from those of `t` when
    # the target moves
    tg = sol.model.target
    for ci in foils
        c = t.components[ci]
        es = t.comp_elems[ci]
        push!(front[c.name], sum(_front_mass(model(e).front, I[e]) / areal_mass(model(e)) * tg.elem_area[e]
                                 for e in es) / sum(tg.elem_area[es]))
    end
    push!(window_loss, loss)
    push!(setup, uc === nothing ? P0 : P0 + udel)
    push!(lost, bal.lost)
    v = ForwardDiff.value
    # The re-emission of the wall is part of the source of the solution, so the
    # x-ray source is booked as converted laser (ηP) and the absorbed power net of
    # the emission.
    emitted = xe === nothing ? nothing : xe .* tg.elem_area
    for (ci, c) in enumerate(tg.components)
        es = tg.comp_elems[ci]
        Ec = emitted === nothing ? 0.0 : sum(emitted[e] for e in es; init = 0.0)
        Pl = sum(v(sol.laser[e]) for e in es; init = 0.0) + v(dep.backside[ci])
        uc === nothing || (Pl += sum(uc[col.name].backside[ci] for col in ulight.colors; init = 0.0))
        push!(laser[c.name], Pl)
        push!(source[c.name], sum(v(sol.source[e]) * tg.elem_area[e] for e in es; init = 0.0) - Ec)
        push!(released[c.name], sum(rel[e] for e in es; init = 0.0))
        haskey(absorbed, c.name) && push!(absorbed[c.name], bal.absorbed[c.name] - Ec)
    end
    return nothing
end

# Returns the records of a history as its result, with the burnthrough times of the foils
# and the pieces of split parts summed into their parts.
function _assemble_history_output(cx, st::_MarchState, rec)
    (; t, obs, we, cavity_elems, weights, N, foils, foil_elems, transmitters, tracked, smooth,
     lips, need, lossy, lossy_index, holes, apertures, partof, names) = cx
    (; transition, flux_change, emission, ulight, min_step, tol, lmax, ylm, observed_flux,
     keep_albedo, wall_release, per_component) = cx.options
    model(e) = element_material(t, e).albedo
    wl(e) = (m = model(e); m isa BurnthroughAlbedo ? m.before : m)
    (; time, delivered, flux, Tr, wall_Tr, wall_albedo, cavity_Tr, moments, clm, Hobs, escaped,
     transmitted, front, window_loss, albedo_steps, transmission_steps, recorded, setup, lost,
     laser, source, released, absorbed, lip_radius, mband_moments, mband_laser, uc_delivered) = rec
    (; tb) = st
    record = cx.options.has_record ? recorded : nothing
    burnthrough = Dict(t.components[ci].name => [_nan_inf(tb[e]) for e in t.comp_elems[ci]] for ci in foils)
    if !per_component && any(c -> c.part !== c.name, t.components)   # pieces of split parts summed into them
        escaped, transmitted, laser, source, absorbed, released =
            _by_part(t, escaped), _by_part(t, transmitted), _by_part(t, laser), _by_part(t, source),
            _by_part(t, absorbed), _by_part(t, released)
        emission === nothing || (mband_laser = _by_part(t, mband_laser))
    end
    return (; time, delivered, flux, Tr, moments = _columns(moments, lmax + 1),
            ylm = ylm ? _columns(clm, (lmax + 1)^2) : nothing,
            observed_flux = observed_flux ? _columns(Hobs, length(obs)) : nothing, wall_Tr,
            wall_albedo, cavity_Tr, escaped, transmitted, window_loss, front, burnthrough,
            albedo = keep_albedo ? albedo_steps : nothing,
            transmission = keep_albedo ? transmission_steps : nothing, setup, laser, source, absorbed, lost,
            released, lip_radius, record = record === nothing ? nothing : [r for r in recorded],
            mband_moments = emission === nothing ? nothing : _columns(mband_moments, lmax + 1),
            mband = emission === nothing ? nothing : mband_laser,
            unconverted = ulight === nothing ? nothing :
                          (; delivered = uc_delivered, skipped = copy(ulight.skipped)))
end

"""
    history_state(target, history, time; mean = ()) -> NamedTuple

The radiation state at `time` (ns) of a [`drive_history`](@ref) of `target`, run with
`keep_albedo = true`. A linear scenario at that time of the pulse, such as a phase of a
shaped pulse, needs this state in place of the albedos of the materials.

- `time` — the time;
- `albedo` — the albedo of each element, the converged albedo of the history
  (`history.albedo`), interpolated linearly in time between the two steps around `time`.
  At the time of a step it is the value of that step. It follows the temperature history
  of each wall element. It includes the emission of a re-emitting wall (the effective
  albedo `B = ηP/A + αH`) and the state of a thin foil (`(1 − s𝒯) α_before` while it
  burns through);
- `transmission` — the fraction of the incident flux on each element transmitted out of
  the enclosure (the `transmission` of the `RadiationSolution`, from burned-through foils
  and open gaps), interpolated the same way;
- `windows` — per thin-foil (`BurnthroughAlbedo`) component, `(; transition,
  transmission, burned)`. These are the area mean of the transition `s` of its elements
  at `time` (0 intact, 1 burned through, from the `burnthrough` times of the history and
  the `width` of the foil), the area mean of their transmitted fraction, and whether
  every element has finished its transition.

`mean` (a component or part name, or a collection of them) replaces the albedo of their
elements by its area mean, one constant over all of them. `mean = :wall` gives the
`wall_albedo` of the history on the wall.

`time` and `albedo` are the keywords of [`Scenario`](@ref) and `remake(::Scenario)`.
With `st = history_state(target, h, τ)`, `Scenario(target, beams; observe = :capsule,
st.time, st.albedo)` or `remake(scenario; st.time, st.albedo)` is the flux response about
the state of the history at `τ`, linear in the laser source as the optimizer requires.
At the time of a step and with the history beams, it gives the observed flux of the
history. A standalone `solve` at a time partway through a pulse, or the default albedos
of a scenario, evaluate a [`WallLossAlbedo`](@ref) as if the wall had been at its
present (or reference) temperature since `t = 0`, without the history or the
re-emission. On a NIF-like mock-up the capsule drive then differs from that of the
history by −15 to +64%, and P2/P0 by up to 3.6 percentage points (research notes, "Drive
histories: sliver steps and the albedo of a phase (2026-11-18)").
"""
function history_state(t::Target, h, τ::Real; mean = ())
    alb = get(h, :albedo, nothing)
    alb === nothing && throw(ArgumentError("the history needs keep_albedo = true"))
    ts = h.time
    isempty(ts) && throw(ArgumentError("the history has no steps"))
    ts[1] <= τ <= ts[end] ||
        throw(ArgumentError("time $τ ns lies outside the history's steps, $(ts[1])–$(ts[end]) ns"))
    N = nelements(t)
    length(alb[1]) == N ||
        throw(DimensionMismatch("the history has $(length(alb[1])) elements, the target $N"))
    # the steps around τ and the weight of the later one (0 at the time of a step)
    k = searchsortedlast(ts, τ)
    w = k == length(ts) ? 0.0 : (τ - ts[k]) / (ts[k+1] - ts[k])
    at(v) = w == 0 ? copy(v[k]) : (1 - w) .* v[k] .+ w .* v[k+1]
    α, tr = at(alb), at(h.transmission)
    names = mean isa Symbol ? (mean,) : Tuple(mean)
    if !isempty(names)
        es = elements(t, names...)
        A = t.elem_area[es]
        α[es] .= sum(α[es] .* A) / sum(A)
    end
    windows = Dict{Symbol,NamedTuple{(:transition, :transmission, :burned),Tuple{Float64,Float64,Bool}}}()
    for (ci, c) in enumerate(t.components)
        c.kind === :surface && c.material.albedo isa BurnthroughAlbedo || continue
        es = t.comp_elems[ci]
        A = t.elem_area[es]
        tb = h.burnthrough[c.name]
        s = [_burnthrough_step(c.material.albedo, τ, x) for x in tb]
        windows[c.name] = (; transition = sum(s .* A) / sum(A), transmission = sum(tr[es] .* A) / sum(A),
                           burned = all(==(1), s))
    end
    return (; time = Float64(τ), albedo = α, transmission = tr, windows)
end

# Returns the computed entrance-hole lips of a target (`LipClosure` apertures). For each
# lip it gives the name, model, hole center and normal, initial (outer) radius `r0`, inner
# radius `rin` (the floor of the closure), and elements. It also gives the radial span
# [ra, rb] of each element from the center, the extreme radii of its vertices in the
# plane of the hole.
function _lip_geometry(t::Target)
    lips = []
    for (ci, c) in enumerate(t.components)
        c.kind === :aperture && c.material.albedo isa LipClosure || continue
        s = c.surface
        s isa DiskSurface || throw(ArgumentError("a LipClosure lip needs a DiskSurface (component $(c.name))"))
        es = t.comp_elems[ci]
        radius(v) = (d = v - s.center; norm(d - dot(d, s.normal) * s.normal))
        span(e) = extrema(radius(v) for k in t.elem_tris[e]
                          for v in (t.tri_v0[k], t.tri_v0[k] + t.tri_e1[k], t.tri_v0[k] + t.tri_e2[k]))
        sp = span.(es)
        ra, rb = first.(sp), last.(sp)
        push!(lips, (; name = c.name, model = c.material.albedo, center = s.center, normal = s.normal,
                     r0 = maximum(rb), rin = minimum(ra), elems = es, ra, rb,
                     area = t.elem_area[es]))
    end
    return lips
end

# Returns the incident-flux radiation temperature of each lip, from the area mean of H
# over its elements.
_lip_temperatures(lips, sol) =
    [radiation_temperature(sum(ForwardDiff.value(sol.H[e]) * A for (e, A) in zip(L.elems, L.area)) / sum(L.area))
     for L in lips]

# Returns the laser intensity through the hole of each lip (TW/mm² over its initial
# area). It is the power of the beams entering through the hole, whose ports are on the
# side of the hole from the chamber center.
_lip_intensities(lips, beams) =
    [sum((ForwardDiff.value(b.power) for b in beams if dot(b.port, L.center) > 0); init = 0.0) / (π * L.r0^2)
     for L in lips]

# Returns the wall-release Duhamel sum at time τ over the steps taken, for the element at
# position e of the deficit vectors (those of the wall-loss elements of the history). It is
# the deficit D, piecewise linear between the steps, convolved with the integral of the
# kernel,
# Σₖ (Dₖ − Dₖ₋₁)/(tₖ − tₖ₋₁) [(τ − tₖ₋₁)^b − (τ − tₖ)^b].
function _deficit_sum(m::WallLossAlbedo, D, ts, e, τ)
    b = 1 + m.time_exponent
    G = 0.0
    for k in 2:length(ts)
        dD = D[k][e] - D[k-1][e]
        dD == 0 && continue
        G += dD / (ts[k] - ts[k-1]) * ((τ - ts[k-1])^b - (τ - ts[k])^b)
    end
    return G
end

# Sums per-step series keyed by component into the parts of the components.
function _by_part(t::Target, d::AbstractDict)
    out = Dict{Symbol,Vector{Float64}}()
    for c in t.components
        haskey(d, c.name) || continue
        v = d[c.name]
        out[c.part] = haskey(out, c.part) ? out[c.part] .+ v : copy(v)
    end
    return out
end

"""
    energy_budget(history; from = nothing, to = nothing) -> NamedTuple

Energies (kJ = TW·ns) of a `drive_history`, the trapezoidal integrals of its power
records over the solved times. With `[from, to]`, the integrals run over that interval,
and outside the solved times the powers count only through the trapezoids that they
bound. The energies are:

- `setup` — the laser power that `beams_at(t)` delivers (before the laser–plasma
  model); `coupled` — after it (`delivered`); `backscattered` — their difference, which
  is the power that `Backscatter` and other lossy models remove (CBET conserves power);
- `laser` — absorbed laser energy per component (front elements plus back side),
  `unabsorbed` — `coupled − Σ laser`, the light that leaves through the holes;
- `source` — x-ray energy converted per component (Σ η P); `absorbed` — x-ray energy
  absorbed per surface component (`(1 − α) H A`, for a `CapsuleAlbedo` capsule at its
  effective radius, ``E_{cap}``); `escaped` per aperture; `transmitted` per thin-foil or
  closing-gap component; `lost` to open space;
- `heating` — per surface component, the energy that stays in it, which is unconverted
  laser plus absorbed x-rays, `laser − source + absorbed`;
- `released` — per component, the stored wall energy re-emitted under a falling drive
  (the `wall_release` of `drive_history`). `absorbed` is net of it: a wall component took
  in `absorbed + released` and holds `absorbed`.

`setup ≈ backscattered + unabsorbed + Σ heating + Σ escaped + Σ transmitted + lost`, to
the closure of the view factors.
"""
function energy_budget(h; from = nothing, to = nothing)
    t = h.time
    keep = [(from === nothing || τ >= from) && (to === nothing || τ <= to) for τ in t]
    tt = t[keep]
    ∫(p) = (q = p[keep]; sum((tt[k+1] - tt[k]) * (q[k] + q[k+1]) / 2 for k in 1:length(tt)-1; init = 0.0))
    d(x) = Dict(n => ∫(v) for (n, v) in x)
    setup, coupled = ∫(h.setup), ∫(h.delivered)
    laser, source, absorbed = d(h.laser), d(h.source), d(h.absorbed)
    heating = Dict(n => laser[n] - source[n] + absorbed[n] for n in keys(absorbed))
    return (; setup, coupled, backscattered = setup - coupled, laser,
            unabsorbed = coupled - sum(values(laser); init = 0.0), source, absorbed,
            escaped = d(h.escaped), transmitted = d(h.transmitted), lost = ∫(h.lost), heating,
            released = haskey(h, :released) ? d(h.released) : Dict{Symbol,Float64}())
end

_columns(v, n) = isempty(v) ? zeros(n, 0) : reduce(hcat, v)
_nan_inf(x) = isnan(x) ? Inf : x

# A value computed from a target and kept while the target stays the same object.
mutable struct _TargetCache{F}
    f::F
    target::Union{Nothing,Target}
    value::Any
end
_TargetCache(f) = _TargetCache(f, nothing, nothing)
function (c::_TargetCache)(t::Target)
    c.target === t || (c.value = c.f(t); c.target = t)
    return c.value
end
