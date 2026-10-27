# Optimization

[`GroupPointing`](@ref) maps parameters to beams, per cone, quad or beam: axial, radial
and tangential aim, defocus, power, spot scale and Δλ, with top–bottom mirroring;
[`parameters`](@ref) builds its parameter vector from values named by group and field,
and [`Scenario`](@ref) precomputes everything an evaluation needs.
[`symmetry_objective`](@ref) and [`asymmetry_objective`](@ref) target capsule drive
symmetry (axisymmetric and three-dimensional), and [`irradiance_objective`](@ref) foil
uniformity. [`optimize_pointing`](@ref) runs bounded L-BFGS-B on forward-mode
automatic-differentiation gradients, with a callback for progress and cancellation.
[`influence_matrix`](@ref) gives each beam's mode response.

[`mode_sensitivity`](@ref) gives the exact sensitivity of the drive shape to each beam's
pointing and power errors, and [`mode_uncertainty`](@ref) the asymmetry independent
errors leave by linear propagation ([`sample_modes`](@ref) is the Monte Carlo check).
Errors can be shared by a group of beams or be the target's positioning. The
objectives' `uncertainty` keyword optimizes the expected symmetry over them, and
[`infer_errors`](@ref) finds the most probable errors behind a measured asymmetry.

```@autodocs
Modules = [Briareus]
Pages = ["src/optimize.jl", "src/uncertainty.jl"]
Private = false
```
