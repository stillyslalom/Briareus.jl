# Physics background

This page collects the physics a Briareus calculation rests on, at the level needed to
set one up and to judge its results. It is written for a reader who knows electromagnetism
and thermal radiation but not inertial-confinement fusion. Each section ends with what
Briareus computes for that piece of physics and where the tutorial and the reference
documentation for it are. The tutorials derive and check the formulas used here; the
reference texts at the end treat each subject in full.

## Indirect drive and the hohlraum

In indirect drive, a fuel capsule is imploded by x rays rather than by laser light
directly. The capsule sits at the center of a hohlraum, a high-Z cylinder (gold, or
depleted uranium lined with gold) a few millimeters across. Laser beams enter through a
laser entrance hole (LEH) at each end and strike the inner wall. The wall's laser-heated
plasma converts the absorbed light into soft x rays, with a conversion efficiency of
roughly 0.7–0.95 depending on intensity and time (research notes, "Literature review against Kritcher et al. 2024"). The
x rays fill the cavity, are partly absorbed by the wall and partly re-emitted, and
bathe the capsule. The capsule's ablator absorbs them, ablates outward, and the reaction
drives the shell and fuel inward. NIF's ignition-era hohlraums reach radiation
temperatures near 300 eV: Dante measured about 307–313 eV on the HYBRID-E shots
(Abu-Shawareb et al., Phys. Rev. Lett. 132, 065102 (2024), Table I).

The cavity is what makes the drive uniform. The capsule sees the wall from a distance
comparable to the hohlraum radius, so each point of the capsule receives flux from a
wide patch of wall, and the fine structure of the laser spots averages out. What
remains are the low-order variations set by where on the wall the laser power lands,
which the beam layout controls. The price is efficiency: most of the x-ray energy is
absorbed by the wall or escapes through the LEHs rather than reaching the capsule
(Lindl 1995; Lindl et al. 2004, Sec. IV).

NIF's 192 beams are grouped into 48 quads of four. The indirect-drive ports lie on four
cones in each hemisphere, at polar angles of 23.5° and 30° (the *inner* cones) and
44.5° and 50° (the *outer* cones). Outer-cone beams strike the wall near the LEHs;
inner-cone beams cross the hohlraum at a shallower angle and strike it near the waist.
The share of the laser power carried by the inner cones, the **cone fraction**, is the
main control of the pole-to-waist balance of the drive. The HYBRID-E designs run a
peak cone fraction near 33 % (research notes, "Wavelength detuning (NIF)").

**In Briareus.** A [`Target`](@ref) is a list of [`Component`](@ref)s: meshed surfaces
with a [`Material`](@ref) (albedo and conversion efficiency) and, for curved walls, the
exact surface the mesh approximates. LEHs are [`aperture`](@ref)s, transparent elements
that record the flux escaping through them. [`hohlraum`](@ref) builds the standard
cylinder with LEHs and capsule, and [`nif_quads`](@ref) and [`nif_beams`](@ref) give the
NIF beams with their cones as `group`. The tutorial [Freeing a bundle for
radiography](@ref) works through a NIF hohlraum; the reference pages are [Geometry and
targets](@ref) and [Beams and facilities](@ref).

## Laser spots

A high-power laser beam focused by a lens alone has a spot shaped by the aberrations
of the beam. A phase plate before the lens breaks the beam into many small, mutually
phased elements, each of which illuminates the whole focal region. Averaged over the
resulting speckle, the spot is smooth, with a shape set by the plate: elliptical and
flat-topped for NIF's continuous phase plates (CPPs), round for most of OMEGA's
distributed phase plates (DPPs). A super-Gaussian of order ``n`` with semi-axes ``a``
and ``b`` describes these spots well,

```math
I(u, v) \propto \exp\!\left[-\left(\frac{u^2}{a^2} + \frac{v^2}{b^2}\right)^{n/2}\right],
```

where ``u`` and ``v`` are transverse coordinates. Briareus quotes ``a`` and ``b`` at
1/e of the peak. Facility tables quote radii at other levels: VisRad's NIF CPP table
gives half widths at half maximum, and a radius ``r_f`` at the fraction ``f`` of the
peak converts as ``a = r_f/(\ln 1/f)^{1/n}``.

The **f-number** ``F`` is the ratio of the focal length to the lens diameter; light
converges onto the focus from a cone of half-angle ``\arctan(1/2F)``. A NIF beam has
``F \approx 20``, and a quad of four beams acts as a beam of ``F \approx 8``. Away from
best focus each point receives light from a disk of focal points, of radius
``|s|/2F`` at a distance ``s`` along the beam, so in geometric optics the defocused spot
is the focal profile convolved with that disk. **Defocus** moves best focus along the
beam on purpose, enlarging the spot on the wall to lower its intensity. In Briareus and
in VisRad a positive defocus moves the focus along the propagation direction.

A beam striking a surface at incidence angle ``\iota`` spreads its power over a
footprint stretched by ``1/\cos\iota``, so the absorbed power per unit area is the beam
intensity times ``\cos\iota``. In a hohlraum, inner-cone beams strike the wall at
grazing incidence, and their footprints are long and faint.

**In Briareus.** A [`Beam`](@ref) has a port (which fixes its direction), an aim, a
defocus, a power, an f-number and a [`SuperGaussianSpot`](@ref). A
[`DepositionOperator`](@ref) carries out the f-cone convolution by summing the focal
profile over a quadrature of cone directions at points fixed on the target; it conserves
each beam's power and is differentiable in aim, defocus and power. By default each lit surface
absorbs all the light that reaches it; no reflection is followed, and backscatter is a
laser–plasma correction (below). [`raytrace_deposition`](@ref) is the reference. The
tutorial [Laser spot on a foil](@ref) checks the spot's power, its growth through focus
and its footprint at oblique incidence against exact results; the reference page is
[Laser deposition](@ref).

## View factors and radiosity

Hohlraum walls emit diffusely, with a radiance independent of direction (Lambertian
emission), and Briareus treats the radiation as gray, a single energy group, unless a
multigroup solve is requested. The exchange of diffuse radiation between surfaces is
then purely geometric. The **view factor** ``F_{ij}`` is the fraction of the radiation
leaving surface element ``i`` that arrives at element ``j``:

```math
F_{ij} = \frac{1}{A_i}\int_{A_i}\!\int_{A_j} \frac{\cos\theta_i\,\cos\theta_j}{\pi r^2}\,
V(\mathbf{x}_i, \mathbf{x}_j)\; dA_j\, dA_i ,
```

where ``\theta_i`` and ``\theta_j`` are the angles between the line joining the two
points and the surface normals, ``r`` its length, and ``V`` is 1 where the two points see
each other and 0 where another surface blocks the line. Two identities follow from the
definition. **Reciprocity**, ``A_i F_{ij} = A_j F_{ji}``, holds because the integrand
is symmetric. **Closure**, ``\sum_j F_{ij} = 1``, holds in a closed enclosure, where all
radiation leaving an element lands somewhere; a hohlraum is closed once its LEHs are
counted as elements. A discretized view-factor matrix satisfies closure only to its
quadrature error, and that error is a useful measure of its accuracy.

The **radiosity** ``B_i`` is the power per unit area leaving element ``i``, and the
incident flux ``H_i`` the power per unit area arriving there. An element that absorbs
laser power ``P_i`` over area ``A_i``, converts a fraction ``\eta_i`` of it to x rays,
and re-emits a fraction ``\alpha_i`` (the **albedo**) of the x rays it receives obeys

```math
B_i = \frac{\eta_i P_i}{A_i} + \alpha_i H_i, \qquad H_i = \sum_j F_{ij} B_j .
```

The two equations close into the linear system ``(I - \operatorname{diag}(\alpha)F)B =
S`` for the source ``S_i = \eta_i P_i / A_i``. An LEH has ``\alpha = \eta = 0``: the
flux it receives is the flux that escapes. The radiation temperature of a flux is
``T_r = (H/\sigma)^{1/4}``.

The albedo of a hohlraum wall is not a surface property in the ordinary sense. X rays
absorbed by the wall drive a radiative heat wave (a Marshak wave) into it, and the
energy that wave stores is the wall loss; what the wall does not store, it re-emits.
The loss falls as the wall heats and as time goes on, so the albedo rises through the
pulse, to about 0.8–0.9 at peak drive for gold. Briareus's default for gold is Lindl's
power-law fit, ``1 - \alpha = 0.32\,(T_r/100\ \mathrm{eV})^{-0.7}(t/1\ \mathrm{ns})^{-0.38}``
(Lindl et al. 2004, Sec. IV.A; Hammer and Rosen, Phys. Plasmas 10, 1829 (2003), for the
self-similar wave).

A sphere is the one cavity with a closed-form answer. Every point of a sphere's
interior sees every other with the same kernel, so the incident flux is uniform
wherever the source lies, and with a port taking a fraction ``f`` of its area,

```math
H = \frac{\eta P}{A}\,\frac{1}{1 - \alpha(1 - f)} .
```

The factor ``1/[1 - \alpha(1 - f)]`` is the cavity's multiplier: about 9 at
``\alpha = 0.95`` and ``f = 6\,\%``. It also multiplies any radiation a discretization
loses or creates, so a closure error of 10⁻³ per row matters at high albedo. At
``\alpha = 0.925``, a 5 × 10⁻⁴ shift of the wall's mean row sum moves the capsule drive
of a NIF hohlraum by 0.3 % and its P2 by about 0.05 % ([`RadiationModel`](@ref)).

The flux on the capsule is linear in the laser deposition, as long as the albedos are
fixed: ``H_{\text{capsule}} = R\,P`` for a response matrix ``R`` that depends on the
geometry and the albedos alone. A pointing study that changes only where the laser
power lands therefore needs ``R`` once. A temperature-dependent albedo makes the
problem nonlinear, and the solve then iterates the albedos to self-consistency with
the flux; ``R`` is the linearization about one albedo state.

**In Briareus.** [`view_factors`](@ref) computes occluded view factors with the exact
point-to-polygon kernel near the emitter, on the exact curved surface where a
component carries one: on a sphere, flat facets misstate the exchange between
neighboring elements by an amount of order one, and the curved calculation is exact to
round-off. `closed = true` removes the closure error while keeping reciprocity
([`enforce_closure!`](@ref); [`closure_error`](@ref) reports it). [`solve`](@ref) solves
the radiosity system, [`WallLossAlbedo`](@ref) is the Marshak wall loss, and
[`flux_response`](@ref) is the matrix ``R``. The tutorials [Integrating sphere](@ref) and
[Entrance-hole size and the power balance](@ref) check these against the sphere's
closed form, the uniform-flux power balance and a Monte Carlo photon walk; the
reference page is [Radiation](@ref).

## Drive symmetry

A capsule converges to a small fraction of its initial radius, and a small variation
of the drive around it grows into a large distortion of the compressed fuel. The
tolerance is of order 1 % in the time-integrated flux for the lowest modes (Lindl 1995;
Lindl et al. 2004). Variations of the flux over the capsule are measured by its
expansion about the capsule center. For a hohlraum with azimuthal symmetry, a Legendre
expansion in the polar angle ``\theta`` from the hohlraum axis suffices,

```math
H(\theta) = \sum_l a_l P_l(\cos\theta), \qquad
a_l = \frac{2l + 1}{4\pi}\int H\,P_l(\cos\theta)\,d\Omega ,
```

and the asymmetry of mode ``l`` is the ratio ``a_l/a_0``, written P2/P0 for ``l = 2``.
Because ``P_2(\pm 1) = 1`` and ``P_2(0) = -1/2``, a positive P2/P0 means the poles are
hotter than the waist (pole-hot, an oblate implosion) and a negative one that the waist
is hotter (waist-hot, a prolate implosion). Real spherical harmonics ``Y_{lm}`` extend
the expansion to flux that varies in azimuth as well, such as the mode 1 left by a
missing bundle or a power imbalance between hemispheres.

The geometry filters the modes. Each capsule point integrates the wall's emission over a
wide solid angle, which suppresses high orders; what reaches the capsule is dominated
by ``l = 2`` and ``l = 4``, with P6 smaller still. A drive that is symmetric between
the two hemispheres has no odd Legendre modes. The low modes have different
consequences. Mode 1 pushes the hot spot sideways: about 100 km/s of hot-spot velocity
per percent of mode-1 flux over the rise and peak (Casey et al., Phys. Rev. Lett. 126,
025002 (2021); MacGowan et al., High Energy Density Phys. 40, 100944 (2021); Young et
al., Phys. Plasmas 27, 082702 (2020)), a velocity that carries kinetic energy the hot
spot no longer converts to heat. P2 and P4 distort the hot spot and shell, with the
same loss of stagnation pressure and yield.

Pointing moves P2 because each cone's power lands at a latitude of the wall as seen
from the capsule. Raising the cone fraction puts more power near the waist and drives
P2 negative; moving a ring's aim along the axis moves its spots in latitude. With two
rings per hemisphere, the cone fraction and the two rings' axial aims are enough to
null P2 and P4 together, as the first calculation on the home page does.
On NIF, about 2.7 µm of hot-spot P2 follows from each percent of cone fraction
(Callahan et al., Phys. Plasmas 25, 056305 (2018)).

**In Briareus.** [`legendre_moments`](@ref) and [`ylm_moments`](@ref) integrate the
harmonics over each facet's exact solid angle as seen from the capsule center (a
centroid rule leaks P0 into P2 at the 0.5 % level); [`mode_amplitude`](@ref) gives the
rotation-invariant size of each degree, and [`flux_moments`](@ref) the moments of a
scenario's capsule flux. [`symmetry_objective`](@ref) and [`optimize_pointing`](@ref)
null chosen modes at a held drive, and [`hotspot_velocity`](@ref) converts mode 1 to a
velocity. The tutorials are [Freeing a bundle for radiography](@ref) (the port symmetry,
mode 1, rebalancing), [Pointing and power errors](@ref) and [From drive asymmetry to
implosion metrics](@ref); the reference pages are [Mode analysis](@ref),
[Optimization](@ref) and [Implosion metrics](@ref).

## Time dependence

The laser pulse of a NIF hohlraum shot lasts several nanoseconds (5.5–14.1 ns over the
39 shots of Callahan et al. 2018). A low-power picket launches the
first shock, a trough holds the drive, and a rise to peak power sets the final
implosion velocity. Each cone has its own pulse, so the cone fraction varies through
the shot. Light crosses a hohlraum in tens of picoseconds, so the radiation field
follows the drive at each instant; the slow parts of the problem are the wall's stored
energy and the hohlraum's plasma.

The hohlraum changes during the pulse in ways that move the drive symmetry:

- **Wall motion.** The ablated wall plasma moves the x-ray emission surface inward,
  at about 50–100 µm/ns (research notes, "Literature review against Kritcher et al.
  2024").
- **LEH closure.** Gold blown off the LEH rims shrinks the holes, lowering the loss
  through them and moving the outer-cone footprints.
- **The gold bubble.** The outer cones' picket blows off a bubble of high-Z plasma under
  their footprint. Its tip moves inward at 300–460 µm/ns, and it absorbs inner-cone
  light before that light reaches the waist, so the inner-cone deposition moves toward
  the LEHs and the drive turns pole-hot (Callahan et al. 2018). On N221204 the capsule
  flux's P2/P0 swings from waist-hot to pole-hot during peak power (Kritcher et al.,
  Phys. Rev. E 109, 025204 (2024)).
- **Cross-beam energy transfer (CBET).** Where beams cross in the LEH plasma, the beat
  of two beams drives an ion-acoustic wave that scatters power from one into the
  other: from the higher frequency to the lower, in the frame of the flowing plasma. A
  wavelength detuning Δλ between cones sets the direction and size of the transfer. With
  Δλ = λ(inner) − λ(outer) positive, power moves from the outer cones to the inner ones;
  on NIF it is quoted at 1ω, three times its 3ω value. N221204 ran 2.75 Å at 1ω (research
  notes, "Wavelength detuning (NIF)"; Michel et al., Phys. Rev. Lett. 102, 025004
  (2009)).
- **Backscatter.** Stimulated Brillouin and Raman scattering return part of a beam's
  power out of the LEH before it reaches the wall (Kruer 1988).

**What Briareus models, and what it does not.** The drive is quasi-static: at each time
the radiosity is solved for that time's deposition, wall albedos and geometry. Each wall
element's loss follows its own temperature history. Wall motion, LEH closure, the bubble,
CBET and backscatter are reduced models, calibrated to published simulations and
measurements that the research notes record. Briareus solves no hydrodynamics, traces
no refraction, and stores no energy in the hohlraum plasma; after the laser turns off
the synthetic Dante therefore reads low, by 37–54 eV at 9 ns on three shots
([Verification and validation](@ref)). The CBET model assumes a uniform LEH plasma and
must be calibrated before its magnitudes are trusted. The implosion itself is outside
the code: [`capsule_drive`](@ref) hands the capsule's flux history to a hydrodynamics
code.

**In Briareus.** [`PulseShape`](@ref) and [`with_pulses`](@ref) set powers in time,
and [`drive_history`](@ref) runs the quasi-static history. [`WallMotion`](@ref),
[`LipClosure`](@ref) and [`GapClosure`](@ref), [`BubbleAbsorption`](@ref),
[`CBET`](@ref) and [`TwoZoneCBET`](@ref), and [`Backscatter`](@ref) are the reduced
models; `Beam.dlambda` is the detuning at the beam's own (3ω) wavelength.
[`dante1`](@ref) reads the drive as the Dante instrument does. The tutorial
[Retuning a drive to a reference shot](@ref) runs a history; the reference pages are
[Radiation](@ref) and [Laser–plasma models](@ref), and `examples/n221204.jl` is the
benchmark against N221204.

## Further reading

- J. D. Lindl, "Development of the indirect-drive approach to inertial confinement
  fusion and the target physics basis for ignition and gain", Phys. Plasmas 2, 3933
  (1995). The physics basis of indirect drive: hohlraum energetics, symmetry and
  capsule stability.
- J. D. Lindl et al., "The physics basis for ignition using indirect-drive targets on
  the National Ignition Facility", Phys. Plasmas 11, 339 (2004). The NIF-era update;
  Sec. IV covers hohlraum energetics and wall loss.
- S. Atzeni and J. Meyer-ter-Vehn, *The Physics of Inertial Fusion* (Oxford University
  Press, 2004). A textbook treatment of implosions, ignition and indirect drive.
- M. F. Modest, *Radiative Heat Transfer*, 3rd ed. (Academic Press, 2013), and J. R.
  Howell, M. P. Mengüç and R. Siegel, *Thermal Radiation Heat Transfer*, 6th ed. (CRC
  Press, 2015). View factors, their identities and the radiosity method.
- W. L. Kruer, *The Physics of Laser Plasma Interactions* (Addison-Wesley, 1988).
  Laser absorption, stimulated scattering and the ion-acoustic coupling behind CBET.
- J. J. MacFarlane, "VISRAD—A 3-D view factor code and design tool for high-intensity
  laser and particle beam target experiments", J. Quant. Spectrosc. Radiat. Transfer
  81, 287 (2003). The view-factor method as VisRad implements it; see also
  [For VisRad users](@ref).
- A. L. Kritcher et al., Phys. Rev. E 109, 025204 (2024). The modeling chain behind
  N221204's symmetry tuning, and the benchmark data Briareus checks against.

`reference/research-notes.md` in the repository gives the source and confidence of
every number Briareus uses, and `reference/litreview/` reviews the literature on
hohlraum dynamics, CBET calibration, energetics and asymmetry transfer.
