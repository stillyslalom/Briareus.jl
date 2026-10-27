# Radiation

[`view_factors`](@ref) computes occlusion-aware view factors (the exact point-to-polygon
kernel near the emitter, adaptive toward shared edges; about 0.4 s for a 3,000-element
hohlraum), and [`solve`](@ref) the gray diffuse radiosity. For closed enclosures,
`closed = true` removes discretization leakage while keeping reciprocity. Wall albedo
can be constant, tabulated in time, or the Marshak wall-loss scaling
[`WallLossAlbedo`](@ref), iterated to self-consistency with the local radiation
temperature. [`flux_response`](@ref) precomputes the linear map from deposition to
capsule flux, so optimizers never re-solve. Large models are never densified: above 10⁴
elements the view factors are a [`Float32Matrix`](@ref) with iterative solves, above
1.5 × 10⁴ an [`HMatrix`](@ref). [`spectral_solve`](@ref) adds multigroup radiosity about
a gray solution, with the laser-produced gold M-band.

[`drive_history`](@ref) follows the drive in time from a [`PulseShape`](@ref) per cone,
quad or beam. Each wall element loses energy by its own temperature history and
re-emits stored energy when the drive falls; the wall moves ([`WallMotion`](@ref)), the
entrance holes close as prescribed ([`GapClosure`](@ref)) or computed
([`LipClosure`](@ref)), and the capsule absorbs on a shrinking radius
([`CapsuleAlbedo`](@ref)). [`dante1`](@ref) and [`dante_unfold`](@ref) read the drive as
the Dante instrument would. [`retune_drive`](@ref) finds the power histories that
reproduce a reference shot's capsule drive under per-key power limits; with NIF bundle
31 down, retuning the other 46 quads takes the mode-1 drive asymmetry from 2 % to
0.02 % at the peak.

```@autodocs
Modules = [Briareus]
Pages = ["src/storage.jl", "src/viewfactors.jl", "src/hmatrix.jl", "src/enclosure.jl", "src/radiation.jl", "src/history.jl", "src/spectral.jl", "src/xray_optics.jl", "src/dante.jl", "src/dante_response.jl", "src/dynamics.jl", "src/retune.jl", "src/montecarlo.jl"]
Private = false
```
