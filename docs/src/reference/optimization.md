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
`time_limit` bounds its wall time. Ready-made constraints: [`irradiance_ceiling`](@ref), a
ceiling on the laser irradiance over named parts (its exact value
[`peak_irradiance`](@ref), held itself by default, its smooth stand-in
[`irradiance_norm`](@ref), held alone with `exact = false`), and [`mode_gate`](@ref), a gate per degree on the
[`mode_departure`](@ref) from a reference drive shape, which goes with
[`drive_objective`](@ref) for the largest drive within the gates;
[`clearance_floor`](@ref), each beam's worst-case clearance from a keep-out at or above
its minimum; [`light_ceiling`](@ref) and [`light_floor`](@ref), the laser light on
named parts ([`light_fraction`](@ref)) at most or at least a limit. The laser
irradiance during the picket, the quantity of the gold bubble's published scalings, has
its own ceiling, [`picket_ceiling`](@ref) (Laser–plasma models).

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
ratio and the departure per degree.

```@autodocs
Modules = [Briareus]
Pages = ["src/pointing_summary.jl"]
Private = false
```
