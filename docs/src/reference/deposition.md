# Laser deposition

Two methods deposit laser power on a target. [`DepositionOperator`](@ref) with
[`deposit`](@ref) is a precomputed quadrature over the f-cone directions with
per-direction visibility, refined where each beam can land: each surface point receives,
from each cone direction that reaches it unobstructed, the focal profile where that ray
crosses the best-focus plane, weighted by the ray's own incidence angle. The result is
smooth in pointing, defocus, spot size and power, so ForwardDiff gradients flow through
it; 48 quads on a 3,000-element hohlraum take about 45 ms on one thread, and each beam's
power is conserved to 1e-4. [`raytrace_deposition`](@ref), an energy-conserving
geometric ray trace, is the reference.

Power that clips the outside of a target, misses it or strikes the capsule is accounted
for separately. The operator samples a target's outside only for its total, with Gauss
rules graded toward the entrance-hole rims, where clipped light concentrates; clipped
power agrees with converged ray traces to 1 % from 1e-4 to 5 % of a beam.
[`render`](@ref) draws orthographic images of any per-element field from any direction
or port, and [`line_of_sight_power`](@ref) gives the x-ray power toward a distant
detector.

```@autodocs
Modules = [Briareus]
Pages = ["src/deposition.jl", "src/corona.jl", "src/imaging.jl"]
Private = false
```
