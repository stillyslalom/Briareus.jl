# Briareus.jl

Briareus computes laser deposition, view-factor radiation transport and drive symmetry
for laser-driven high-energy-density experiments: indirect-drive hohlraums,
backlighter foils and other laser-irradiated targets. It covers the calculations a
view-factor code such as VisRad is used for, and it adds gradient-based optimization of
pointing, defocus, power balance and wavelength detuning, reduced laser–plasma
corrections (backscatter, cross-beam energy transfer), and shot-review checks of beam
clearance and unconverted light.

## The model in brief

- **Targets** are assembled from meshed surfaces: cylinders, annuli, spheres,
  rectangles, frusta. Each patch side that radiates is an element of the radiation
  model. Laser entrance holes are transparent apertures that record the flux escaping
  through them.
- **Beams** enter from facility ports with elliptical super-Gaussian phase-plate spots,
  an f-number cone, pointing and defocus. A precomputed quadrature makes deposition a
  smooth function of those parameters, so derivatives come from automatic
  differentiation. An energy-conserving ray trace is the reference.
- **Radiation** is gray and diffuse. View factors account for occlusion, and a
  radiosity solve gives the incident flux on every element. Wall albedo can depend on
  time and on the local radiation temperature. Because the system is linear in the
  laser source, the response of the capsule flux to deposition is computed once and
  reused by optimizers.
- **Symmetry** is measured by Legendre and real spherical-harmonic moments, integrated
  exactly over each facet.

Units are mm, ns, TW (TW·ns = kJ), TW/mm² for flux (`to_Wcm2` multiplies by 10¹⁴)
and eV for temperature. Port angles are in degrees: polar angle from +z, azimuth
counterclockwise from +x, which on NIF points to Switchyard 1.

## Choose a path

- **Learn by running a calculation.** [Laser spot on a foil](@ref) checks deposition
  and defocus against exact results. [Integrating sphere](@ref) checks the radiation
  model against the sphere-multiplier formula. [Freeing a bundle for radiography](@ref)
  chooses a NIF bundle to drive a backlighter and restores the drive symmetry it costs,
  within worst-case beam clearance.
- **Understand what the checks establish.** [Verification and validation](@ref)
  separates the activities the tutorials and tests carry out.
- **Look up exact behavior.** The reference pages document every exported name.

Facility data and physics defaults are sourced in `reference/research-notes.md` in the
repository, with a confidence level for each number.
