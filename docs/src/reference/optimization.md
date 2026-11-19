# Optimization

[`GroupPointing`](@ref) maps parameters to beams, per cone, quad or beam: axial, radial
and tangential aim (in the quad's port azimuth by default) and a quad's split about its
aim point, defocus, power, spot scale and Δλ, with top–bottom mirroring;
[`parameters`](@ref) builds its parameter vector from values named by group and field,
and [`Scenario`](@ref) precomputes everything an evaluation needs.
[`symmetry_objective`](@ref) and [`asymmetry_objective`](@ref) target capsule drive
symmetry (axisymmetric and three-dimensional), [`face_flux_objective`](@ref) the
uniformity of a flat observed face's x-ray drive (a halfraum's package, an opacity
sample), and [`irradiance_objective`](@ref) foil uniformity. [`optimize_pointing`](@ref) runs bounded L-BFGS-B on forward-mode
automatic-differentiation gradients, with a callback for progress and cancellation.
[`influence_matrix`](@ref) gives each beam's mode response.
[`HeldParameters`](@ref) holds some parameters of a parameterization at given values, and
[`PowerScales`](@ref) adds a power scale per group of beams (e.g. per cone) beside its
parameters; both compose with `GroupPointing`, `FramePointing` and each other.

Limits that a design must meet go to [`optimize_pointing`](@ref) as
[`Constraint`](@ref)s, which an augmented Lagrangian around the same L-BFGS-B handles;
[`check_constraints`](@ref) reports them on any evaluation, and the
[`PointingResult`](@ref) records each one's value against its limit, what ended the
run and, for each parameter held by the box, the optimum's derivative with respect to
that bound. `from = result` continues a run with its multipliers and penalty, and
`time_limit` bounds its wall time. `phase1 = true` first looks for a design within the
limits and ends with `:infeasible` when the search finds none, and `relaxation = true`
then reports how far each conflicting constraint must be relaxed alone
([`Relaxation`](@ref)). Ready-made constraints: [`irradiance_ceiling`](@ref), a
ceiling on the laser irradiance over named parts (its exact value
[`peak_irradiance`](@ref), held itself by default, its smooth stand-in
[`irradiance_norm`](@ref), held alone with `exact = false`), and [`mode_gate`](@ref), a gate per degree on the
[`mode_departure`](@ref) from a reference drive shape, which goes with
[`drive_objective`](@ref) for the largest drive within the gates, and
[`drive_floor`](@ref), the capsule's drive at least a fraction of a reference's;
[`clearance_floor`](@ref), each beam's worst-case clearance from a keep-out at or above
its minimum; [`light_ceiling`](@ref) and [`light_floor`](@ref), the laser light on
named parts ([`light_fraction`](@ref)) at most or at least a limit. The outer cones'
irradiance in a phase (the picket, the trough, the foot), the quantity of the gold
bubble's published scalings, has its own ceiling, [`picket_ceiling`](@ref) (Laser–plasma
models).

[`mode_sensitivity`](@ref) gives the exact sensitivity of the drive shape to each beam's
pointing and power errors, and [`mode_uncertainty`](@ref) the asymmetry independent
errors leave by linear propagation ([`sample_modes`](@ref) is the Monte Carlo check).
Errors can be shared by a group of beams or be the target's positioning. The
objectives' `uncertainty` keyword optimizes the expected symmetry over them, and
[`infer_errors`](@ref) finds the most probable errors behind a measured asymmetry.

```@autodocs
Modules = [Briareus]
Pages = ["src/optimize.jl", "src/constraints.jl", "src/uncertainty.jl", "src/parameterizations.jl"]
Private = false
```

[`pointing_summary`](@ref) reports a pointing against reference beams (a re-pointing
result against its start): aim shifts, power and cone fraction, keep-out clearances,
the clipped and missed fractions and, against a reference design's drive, the drive
ratio and the departure per degree. [`shift_objective`](@ref), the mean squared aim shift
from reference beams, gives with gates as constraints the least re-pointing that meets
them.

```@autodocs
Modules = [Briareus]
Pages = ["src/pointing_summary.jl"]
Private = false
```

## Several phases under one pointing

A [`PhasedScenario`](@ref) holds one pointing over several phases of a pulse (a
picket, a foot, the peak before and after windows burn through): each [`Phase`](@ref)
has its own target, radiation model and options and its own beam powers, given per
beam or from pulses at a time, and the parameter vector ([`PhasedPointing`](@ref)) is
the shared pointing's, which sets the aims of every phase, followed by a power scale per
(group, phase) for the phases not fixed. [`evaluate`](@ref) gives one evaluation per
phase in a [`PerPhase`](@ref), which [`optimize_pointing`](@ref) passes to the
objective and to every constraint; [`in_phase`](@ref) applies any constraint or
objective of one evaluation to one phase, the phase in its reported name.
`bounds(scenario.param; power_max)` caps every beam in every phase at the laser's
limit, [`phase_beams`](@ref) puts reference beams at each phase's powers, and
[`pointing_summary`](@ref) reports each phase.

```@autodocs
Modules = [Briareus]
Pages = ["src/phased.jl"]
Private = false
```

## Pulse shapes

Each distinct pulse shape needs facility validation. [`pulse_shape_groups`](@ref) groups
per-key power multipliers over a pulse's phases (one per quad and phase, from power
scales) into a given number of shapes, by k-means on the log multipliers per class (e.g.
the inner and the outer cones), the shapes split between the classes by fit;
[`shape_pulses`](@ref) and [`scaled_pulse`](@ref) write the shapes out as pulses for
[`write_pulse_table`](@ref), and the grouping's `assignment` ties the power scales of a
re-optimization per shape ([`PowerScales`](@ref)).

```@autodocs
Modules = [Briareus]
Pages = ["src/pulse_shapes.jl"]
Private = false
```
