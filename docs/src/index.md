# Briareus.jl

Briareus computes laser deposition, view-factor radiation transport and drive symmetry
for laser-driven high-energy-density experiments: indirect-drive hohlraums, backlighter
foils and other laser-irradiated targets. It covers the calculations that a view-factor
code such as VisRad is used for and exchanges files with VisRad. It adds gradient-based
optimization of pointing, defocus, power balance and wavelength detuning; reduced
laser–plasma corrections (backscatter, cross-beam energy transfer, inner-beam
absorption); a time-dependent drive with a moving wall, closing entrance holes and a
synthetic Dante; and shot-review checks of worst-case beam clearance, NIF's
target-configuration rules and unconverted light.

```julia
pkg> add Briareus
```

A first calculation evaluates a NIF hohlraum at an illustrative cone pointing (aim
height `z` and radial offset `r` per cone, mm; an example, not a shot's pointing), then
moves the inner and outer rings of cones along the axis and re-balances their power to
null the capsule flux's P2 and P4 at the same mean flux, with the beams kept clear of
the entrance-hole rims and the capsule. Once the package has loaded, the scenario
builds in 2.5 s on 16 threads and the optimization takes 5 s on its first call (mostly
compilation) and 1 s after:

```julia
using Briareus

target = hohlraum(length = 10.1, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9)
cones = GroupPointing(nif_quads(power = 8.0); fields = (:z, :r))
quads = apply(cones, parameters(cones; inner23 = (z = 3.09, r = -1.09),
                                inner30 = (z = 3.63, r = -0.60),
                                outer44 = (z = 5.20,), outer50 = (z = 5.01,)))

ring(b) = b.group in (:inner23, :inner30) ? :inner : :outer
pointing = GroupPointing(quads; by = ring, fields = (:z, :power), relative = true)
scenario = Scenario(target, pointing; observe = :capsule)

x0 = initial_parameters(pointing)    # no aim shift, 8 TW per quad
report(evaluate(scenario, x0))       # 384 TW, Tr = 262.7 eV, P2/P0 = -15.737 %

drive = mean_flux(evaluate(scenario, x0))    # mean capsule flux, held fixed
keepout = (ApertureKeepOut(target), SphereKeepOut(target))
objective = symmetry_objective(modes = (2 => 1.0, 4 => 1.0); drive, keepout)
lower, upper = bounds(pointing, x0; power = 0.6)     # ±1 mm of aim, ±60 % of power
result = optimize_pointing(scenario, objective; x0, lower, upper)  # 20 iterations
report(evaluate(scenario, result.x)) # 422.1 TW, P2/P0 = -0.000 %, P4/P0 = +0.000 %
```

The starting P2/P0 of −15.7 % is expected: equal power per quad puts 1/3 of it on the
inner cones, near HYBRID-E's peak cone fraction, and a static solve without cross-beam
energy transfer or the hohlraum's dynamics is not balanced there. The N221204 benchmark
has a static P2/P0 of −26 % at 5.5 ns; the shot's balance comes from wavelength-detuned
transfer, wall motion, the gold bubble and closing entrance holes, which
[`drive_history`](@ref) models. The optimum lowers the inner-cone fraction to 0.15 and
stops where the 30° beams reach the entrance-hole keep-out.

The same workflow runs in a desktop window: see [Desktop application](gui.md).

## Start here

**Evaluating Briareus against VisRad or a design workflow.**
[For VisRad users](@ref) maps VisRad's concepts and files to Briareus's, states where
the two codes differ, and walks through porting a VisRad study.
[Verification and validation](@ref) lists every check by kind (code verification,
solution verification, validation, calibration), including the comparison of laser
deposition with VisRad 21 and the N221204 benchmark against published data.
[Polar direct drive: a published NIF pointing](@ref) reproduces a published pointing
against an independent quadrature. [Retuning a drive to a reference shot](@ref) restores
a reference shot's capsule drive through time after a bundle is lost.
[Radiation in an imported CAD hohlraum](@ref) takes tessellated CAD solids to a
radiation model. [VisRad interchange](@ref) documents the workspace, pointing and power
files read and written.

**Learning the physics with the code.** [Physics background](@ref) introduces indirect
drive, laser spots, view factors, drive symmetry and the hohlraum's time dependence,
with pointers into the code and further reading. The tutorials build up in order, and
each ends with "What this checks", the exact references it is held to.

1. [Laser spot on a foil](@ref): a phase-plate spot, defocus and oblique incidence.
2. [Integrating sphere](@ref): view factors and radiosity in a cavity with an exact
   solution.
3. [Freeing a bundle for radiography](@ref): capsule drive symmetry in a NIF hohlraum,
   the port symmetry, and power rebalancing within beam clearance.
4. One topic each:
   [Entrance-hole size and the power balance](@ref),
   [Diagnostic windows in the hohlraum wall](@ref),
   [Pointing and power errors](@ref),
   [From drive asymmetry to implosion metrics](@ref).

The reference pages document every exported name.

## The model in brief

- **Targets** are assembled from meshed surfaces (cylinders, annuli, spheres,
  rectangles, frusta) or imported CAD. Each patch side that radiates is an element of
  the radiation model. Curved components carry the exact surface their mesh
  approximates, and the radiation model integrates over it. Laser entrance holes are
  transparent apertures that record the flux escaping through them.
- **Beams** enter from facility ports with elliptical super-Gaussian phase-plate spots,
  an f-number cone, pointing and defocus. A precomputed quadrature over the cone
  directions makes deposition a smooth function of those parameters, so derivatives
  come from automatic differentiation. An energy-conserving ray trace is the reference.
- **Radiation** is diffuse, gray or multigroup. View factors account for occlusion, and
  a radiosity solve gives the incident flux on every element. Wall albedo can depend on
  time and on the local radiation temperature. Because the system is linear in the
  laser source, the response of the capsule flux to deposition is computed once and
  reused by optimizers.
- **Symmetry** is measured by Legendre and real spherical-harmonic moments, integrated
  exactly over each facet's solid angle.
- **Time dependence** is quasi-static: at each time, the radiosity is solved for the
  deposition, the wall state and the geometry of that time.

## Glossary

- **View factor** ``F_{ij}``: the fraction of the diffuse radiation leaving element
  ``i`` that arrives at element ``j``, set by geometry and occlusion alone.
- **Radiosity**: the power per unit area leaving a surface, emission plus reflection
  (re-emission, for a hohlraum wall). With the view factors it closes a linear system
  for the flux on every element; [Integrating sphere](@ref) solves the case with an
  exact answer.
- **Albedo**: the fraction of the incident x-ray flux a wall re-emits. Gold hohlraum
  walls reach 0.8–0.9 at peak drive.
- **LEH**: the laser entrance hole at each end of a hohlraum. It admits the beams and
  loses the x rays that reach it.
- **P2/P0** (and Pₗ/P0): the ratio of the capsule flux's Legendre moment of degree
  ``l`` to its mean. P2 > 0 is pole-hot, P2 < 0 waist-hot. Real spherical harmonics
  ``Y_{lm}`` extend this to three-dimensional asymmetry, such as mode 1 from a missing
  bundle.
- **Cone fraction**: the inner cones' share of the total laser power on NIF. It is the
  main control of P2 through the pulse.
- **CBET**: cross-beam energy transfer, the exchange of power between crossing beams
  through ion-acoustic waves in the LEH plasma. A wavelength detuning Δλ between cones
  sets its direction and size. See [Laser–plasma models](@ref).
- **Dante**: the absolutely calibrated, filtered x-ray diode array that measures the
  radiation temperature seen through an LEH.

Lindl et al., Phys. Plasmas 11, 339 (2004) review indirect-drive hohlraum physics; the
tutorials cite the primary sources for each topic.

## Units and conventions

Units are mm, ns, TW (TW·ns = kJ), TW/mm² for flux (`to_Wcm2` multiplies by 10¹⁴)
and eV for temperature. Port angles are in degrees: polar angle from +z, azimuth
counterclockwise from +x, which on NIF points to Switchyard 1.

A beam's direction is fixed by its port (``\mathbf{k} = -\hat{\mathbf{p}}``), and
pointing translates its axis through `aim`. Best focus is at `aim + defocus·k`: positive
defocus moves the focus along the beam, as in VisRad and the NIF CMT. Spot axes are given
in the port's (``\hat\theta``, ``\hat\phi``) frame.

Away from focus, the spot is the focal profile convolved with the f-cone. The deposition
operator samples that convolution at the cone directions: each surface point receives,
from each cone direction that reaches it unobstructed, the focal profile where that ray
crosses the best-focus plane, weighted by the ray's own incidence angle. Power is
conserved at any distance from focus and incidence, and shadows near focus are
correlated with the part of the spot they block.

## Provenance

Facility data and physics defaults are sourced in `reference/research-notes.md` in the
repository, with a confidence level for each number and the data each calibration was
fitted to. Non-public facility data (shot pointing exports, CAD, VisRad workspaces) is
not in the repository; the HYBRID-E pointing in the examples is a public mock-up
([`hybride_mock_pointing`](@ref)).
