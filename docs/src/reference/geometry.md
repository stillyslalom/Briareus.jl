# Geometry and targets

Targets are assembled from meshed primitives (cylinder, annulus or disk, sphere,
rectangle, frustum) with per-component materials; [`hohlraum`](@ref) and [`foil`](@ref)
build the common ones. Surfaces are one- or two-sided, and laser entrance holes are
transparent apertures that measure the escaping x-ray flux. A component can carry the
exact surface its mesh approximates ([`SphereSurface`](@ref), [`CylinderSurface`](@ref),
[`DiskSurface`](@ref)), and the radiation model then integrates over the true curved
surface; hohlraums do this by default.

[`WallWindow`](@ref) cuts rectangular or round diagnostic windows, with optional assembly
gaps, into a hohlraum wall at exact areas, and the mesh conforms to their edges. Closed
windows re-emit with their own material; open ones are apertures. A thin-foil window
burns through when the Marshak heat front reaches its back ([`BurnthroughAlbedo`](@ref))
and then transmits the gray Eddington fraction ``4/(3\tau + 2)`` of the flux on it
(Lester et al., Phys. Rev. E 114, L023201 (2026)), optionally with the optical depth
following the radiation temperature ([`RosselandTransmission`](@ref)). Gaps can close
during the pulse ([`GapClosure`](@ref)). [`Sombrero`](@ref) adds LEH shields.

CAD parts (for example the solids of a multi-solid STL file) become a radiation-ready
target in one call, [`cad_target`](@ref): windings repaired into closed solids facing out
of their material ([`orient_mesh`](@ref)), the cavity around a point kept
([`enclosure`](@ref)), the surface remeshed into compact, near-planar elements of a
chosen size ([`remesh_mesh`](@ref), which also cuts the long slivers CAD tessellates
curved walls into; [`coarsen_mesh`](@ref) only merges whole facets), and cylinders,
spheres, cones and tori fitted ([`fit_surfaces`](@ref)). Primitive components can be
passed with the parts. [Radiation in an imported CAD hohlraum](@ref) walks through each
step.

```@docs
Briareus
```

```@autodocs
Modules = [Briareus]
Pages = ["src/geometry.jl", "src/meshops.jl", "src/target.jl", "src/targets.jl", "src/cadfit.jl"]
Private = false
```
