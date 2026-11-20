# Optimization

[`GroupPointing`](@ref) maps parameters to beams, per cone, quad or beam. The parameters
are the axial, radial and tangential aim (in the port azimuth of the quad by default),
the split of a quad about its aim point, defocus, power, spot scale and Δλ, with
top–bottom mirroring. [`parameters`](@ref) builds its parameter vector from values named
by group and field, and [`Scenario`](@ref) precomputes everything an evaluation needs.
The objectives are:

- [`symmetry_objective`](@ref) and [`asymmetry_objective`](@ref), for capsule drive
  symmetry (axisymmetric and three-dimensional);
- [`face_flux_objective`](@ref), for the uniformity of the x-ray drive on a flat
  observed face, such as a halfraum package or an opacity sample;
- [`irradiance_objective`](@ref), for foil uniformity.

[`optimize_pointing`](@ref) runs bounded L-BFGS-B on forward-mode
automatic-differentiation gradients, with a callback for progress and cancellation.
[`influence_matrix`](@ref) gives the mode response of each beam.
[`HeldParameters`](@ref) holds some parameters of a parameterization at given values,
and [`PowerScales`](@ref) adds a power scale per group of beams (for example per cone)
beside its parameters; both compose with `GroupPointing`, `FramePointing` and each
other.

Limits that a design must meet go to [`optimize_pointing`](@ref) as
[`Constraint`](@ref)s, which an augmented Lagrangian around the same L-BFGS-B handles.
[`check_constraints`](@ref) reports them on any evaluation. The [`PointingResult`](@ref)
records the value of each constraint against its limit and what ended the run, and, for
each parameter held by the box, the derivative of the optimum with respect to that
bound. `from = result` continues a run with its multipliers and penalty, and
`time_limit` bounds its wall time. On [`VisRadProfile`](@ref) with its hard cut, a run
optimizes the profile tapered below the cut and reports the design under the hard cut
(`smooth_cut`). Phase 1 (`phase1`, on by default) first searches for a design within the
limits and ends with `:infeasible` when the search finds none; `relaxation = true` then
reports how far each conflicting constraint must be relaxed by itself
([`Relaxation`](@ref)).

The ready-made constraints are:

- [`irradiance_ceiling`](@ref), a ceiling on the laser irradiance over named parts. By
  default the exact value, [`peak_irradiance`](@ref), is constrained; with
  `exact = false` only the smooth stand-in [`irradiance_norm`](@ref) is.
- [`mode_gate`](@ref), a gate per degree on the [`mode_departure`](@ref) from a
  reference drive shape, which goes with [`drive_objective`](@ref) for the largest drive
  within the gates. Keywords of `mode_gate` set the gates at a fraction of the pointing
  floor, set physical gates with the floor as their lower bound, or shift the target
  shape from that of the reference by [`shift_modes`](@ref).
- [`drive_floor`](@ref), a lower bound on the capsule drive as a fraction of that of a
  reference, and [`drive_ceiling`](@ref), the matching upper bound, for a band.
- [`clearance_floor`](@ref), a lower bound on the worst-case clearance of each beam from a
  keep-out at its minimum.
- [`light_ceiling`](@ref) and [`light_floor`](@ref), an upper and a lower bound on the
  laser light on named parts ([`light_fraction`](@ref)).
- [`picket_ceiling`](@ref) (Laser–plasma models), a dedicated ceiling on the irradiance
  of the outer cones in a phase such as the picket, trough or foot, the quantity of the
  published gold-bubble scalings.

[`mode_sensitivity`](@ref) gives the exact sensitivity of the drive shape to the
pointing and power errors of each beam. [`mode_uncertainty`](@ref) gives the asymmetry
that independent errors leave by linear propagation ([`sample_modes`](@ref) is the Monte
Carlo check), with its `floor`, the spread per degree in the units of
[`mode_departure`](@ref) and the lower bound of a gate. Errors can be shared by a group
of beams or be the positioning error of the target. The `uncertainty` keyword of the
objectives optimizes the expected symmetry over them, and [`infer_errors`](@ref) finds
the most probable errors behind a measured asymmetry.

```@autodocs
Modules = [Briareus]
Pages = ["src/optimize.jl", "src/constraints.jl", "src/uncertainty.jl", "src/parameterizations.jl"]
Private = false
```

[`pointing_summary`](@ref) reports a pointing against reference beams, such as a
re-pointing result against its start. It lists the aim shifts, power and cone fraction,
keep-out clearances, and the clipped and missed fractions, and, against the drive of a
reference design, the drive ratio and the departure per degree.
[`shift_objective`](@ref), the mean squared aim shift from reference beams, gives with
gates as constraints the least re-pointing that meets them.

```@autodocs
Modules = [Briareus]
Pages = ["src/pointing_summary.jl"]
Private = false
```

## Several phases under one pointing

A [`PhasedScenario`](@ref) holds one pointing over several phases of a pulse, such as
the picket, the foot, and the peak before and after windows burn through. Each
[`Phase`](@ref) has a separate target, radiation model, options and beam powers, given
per beam or from pulses at a time. The parameter vector ([`PhasedPointing`](@ref)) is
that of the shared pointing, which sets the aims of every phase, followed by a power
scale per (group, phase) for the phases not fixed. [`evaluate`](@ref) gives one
evaluation per phase in a [`PerPhase`](@ref), which [`optimize_pointing`](@ref) passes
to the objective and to every constraint. [`in_phase`](@ref) applies any constraint or
objective of one evaluation to one phase, with the phase in its reported name, while the
[`clearance_floor`](@ref) of a keep-out takes every phase at once. The albedos and time
of a phase come from a drive history with `Phase(name, scenario; history, time)`, and a
reference at the same time from `remake(scenario; history, time)`.

Further tools for phased designs are:

- `bounds(scenario.param; power_max)` caps every beam in every phase at the limit of the
  laser, and [`pulse_ceiling`](@ref) constrains the pulses written for a design
  ([`scaled_pulse`](@ref)) to stay within it between the phase samples.
- [`phase_beams`](@ref) puts reference beams at the powers of each phase.
- [`power_scales`](@ref) returns the scales inside an objective (differentiable),
  [`power_scale_objective`](@ref) penalizes their departure from a reference drive (each
  scale weighted by the beams it governs), and [`power_scale_table`](@ref) tabulates
  them per key and phase.
- `remake(scenario; deposition_kw)` rebuilds the deposition operator of every phase, for
  example to check a finished design on a finer one.
- [`pointing_summary`](@ref) reports each phase against its gates and pointing floor,
  with the departure per degree split into `m = 0` and `m ≠ 0` and the quads moved past
  given shifts.

```@autodocs
Modules = [Briareus]
Pages = ["src/phased.jl"]
Private = false
```

## Pulse shapes

Each distinct pulse shape needs facility validation. [`pulse_shape_groups`](@ref) groups
per-key power multipliers over the phases of a pulse (one per quad and phase, from power
scales) into a given number of shapes. It runs k-means on the log multipliers per class
(for example the inner and the outer cones) and splits the shapes between the classes by
fit. [`shape_pulses`](@ref) and [`scaled_pulse`](@ref) write the shapes out as pulses
for [`write_pulse_table`](@ref), and the `assignment` of the grouping ties the power
scales of a re-optimization per shape ([`PowerScales`](@ref)).

```@autodocs
Modules = [Briareus]
Pages = ["src/pulse_shapes.jl"]
Private = false
```
