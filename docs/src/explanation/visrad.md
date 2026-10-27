# For VisRad users

VisRad (Prism Computational Sciences; MacFarlane, J. Quant. Spectrosc. Radiat.
Transfer 81, 287 (2003)) is the view-factor code most hohlraum and backlighter pointing
work is done in. Briareus covers the same core calculation (laser deposition, 3D
view-factor radiosity and capsule drive symmetry), reads and writes VisRad's files, and
has been compared with VisRad 21 on NIF hohlraums. This page maps VisRad's concepts to
Briareus's, lists the files each reader and writer handles, summarizes how the two
codes' results compare, and states where they differ, including what VisRad does that
Briareus does not. It ends with a walkthrough of porting a VisRad study. The
[Physics background](@ref) page covers the physics both codes compute, and
[VisRad interchange](@ref) documents every reader and writer.

## Concepts

| VisRad | Briareus | Notes |
|:--|:--|:--|
| Workspace (`.vrw`) | [`VisRadWorkspace`](@ref) as read; then a [`Target`](@ref) and a vector of [`Beam`](@ref)s | A Briareus study holds no workspace object: the target, the beams and the precomputed operators ([`DepositionOperator`](@ref), [`RadiationModel`](@ref), or both in a [`Scenario`](@ref)) are separate values. |
| Target component (object), component tree | [`VisRadObject`](@ref); [`Component`](@ref) | One component per object; a cylindrical hohlraum becomes five (wall, two end caps, two LEH apertures). The tree's folder is `object.folder`, for selecting parts. |
| Surface elements, object grid ("Number of mesh pts") | Elements (radiating patch sides) | [`visrad_component`](@ref) meshes primitives on VisRad's grid, or on Briareus's with `primitive_h`. One-sided surfaces radiate from the front; `two_sided` radiates from both. |
| Albedo, x-ray conversion efficiency | [`Material`](@ref): `albedo`, `conversion` | Constant, tabulated in time, or models: [`WallLossAlbedo`](@ref) (Marshak wall loss), [`Conversion`](@ref). |
| LEH (elements of albedo 0) | [`aperture`](@ref) | Transparent; its incident flux is the flux that escapes. |
| Laser beam, beam cluster | [`Beam`](@ref), its `group` | NIF beams `B111`–`B468` ([`nif_beams`](@ref)) or quads ([`nif_quads`](@ref)); `group` is the cone. |
| Pointing Reference Point, focus offset | `aim`, `defocus` | Same sign: positive moves the focus along the propagation direction. |
| Spatial profile (super-Gaussian), built-in phase plates | [`SuperGaussianSpot`](@ref) (1/e semi-axes); [`VISRAD_PHASE_PLATES`](@ref), [`visrad_spot`](@ref) | VisRad's NIF CPP radii are half widths at half maximum; the readers convert them. |
| Beam profile away from focus | f-cone sum (default); [`VisRadProfile`](@ref) | See "The beam away from focus" below. |
| Time-dependent beam powers, pulse shapes | [`PulseShape`](@ref), [`with_pulses`](@ref), [`drive_history`](@ref) | A pulse is piecewise linear in ns and TW. |
| Direct laser power source | `deposit(op, beams).element` ([`deposit`](@ref)) | TW absorbed per element, not TW/cm²; [`visrad_deposition`](@ref) converts VisRad's. |
| Radiosity results (Flux In, Flux Out) | [`RadiationSolution`](@ref): `H`, `B` | TW/mm²; [`to_Wcm2`](@ref) converts. |
| Surface element statistics (Legendre and ``Y_{LM}`` fits) | [`legendre_moments`](@ref), [`ylm_moments`](@ref) | See "Modes" below; [`visrad_legendre_fit`](@ref) and [`visrad_ylm_fit`](@ref) reproduce VisRad's fits. |
| Capture-clearance report | [`clearance`](@ref), [`ApertureKeepOut`](@ref) | Worst case over pointing and spot tolerances. |

VisRad's workspaces and tables are in cm, TW/cm², µm for spot radii and seconds for
pulse-shape files; Briareus is in mm, TW/mm², ns and TW. The readers and writers convert.

## Files

| File | Read | Written | Notes |
|:--|:--|:--|:--|
| NIF SST pointing CSV | [`read_visrad_sst`](@ref) | [`write_visrad_sst`](@ref) | Aim (without defocus) and defocus per beam; beams absent from the file are turned off, as VisRad does. The format to set pointing in VisRad. |
| Laser Beam Summary CSV | [`read_visrad_summary`](@ref) | [`write_visrad_summary`](@ref) | Ports, pointing, focus, power, spots. The headers after "Power" are not shown in the manual, so the written ones are an assumption. |
| Time-dependent beam powers CSV | [`read_visrad_powers`](@ref) | [`write_visrad_powers`](@ref) | One column per beam, TW against ns; for NIF, VisRad assigns the powers to the 3ω cones. |
| Pulse-shape library file | [`read_visrad_pulse_shape`](@ref) | — | W against s; the column order is an assumption. |
| Workspace `.vrw` (keyed text, VisRad 20–21) | [`read_visrad_workspace`](@ref), [`visrad_target`](@ref) | [`write_visrad_workspace`](@ref) | Components (CAD as generalized meshes; spheres, rectangles, cylindrical hohlraums) and beams. The writer copies a template workspace and rewrites the beams' pointing, focus, power and on/off; everything else stays VisRad's. JSON workspaces are refused. |
| Computed power sources in a `.vrw` | [`visrad_power_sources`](@ref), [`visrad_deposition`](@ref) | — | VisRad's per-element deposition, mapped onto a target built from the same workspace. |
| Surface Element Properties (`.sed`), Power Summary (`.psd`) | [`read_visrad_table`](@ref) | — | The file layout is not documented; the reader is tolerant. |
| Benchmark kit | [`compare_visrad`](@ref) | [`write_visrad_kit`](@ref) | Pointing, beam summary, powers and target description for rebuilding a case in VisRad, with Briareus's results named as VisRad's Power Summary columns. |
| NIF CMT XML, `.bti` intersection points, beam-energy files | — | — | Not implemented: the CMT schema is undocumented and no public sample was found. |

## How the results compare

VisRad 21 computed the laser deposition on two NIF hohlraum workspaces built from
Briareus geometry, and Briareus read each one, built the same target and beams, and
compared its deposition element by element with VisRad's. With
[`VisRadProfile`](@ref), wall deposition differs from VisRad's by 0.9 % rms per
element, peaks agree within 0.5 %, and the capsule's P2 and P4, computed from the two
depositions through the same Briareus radiosity, agree to 0.014 and 0.004 percentage
points. [Verification and validation](@ref) gives the setup and all the numbers. On a
foil workspace with four beams, the totals agree to the table's precision (3.6000 TW),
elements above 20 TW/cm² within 1.4 %, and the faint wings near focus (2–30 TW/cm²)
differ by 8–15 % (research notes, "VisRad workspace files (2026-10-17)").

[`write_visrad_kit`](@ref) and [`compare_visrad`](@ref) set up the same comparison for
any case. The comparison covers deposition only: VisRad's own
radiosity, absolute drive levels and mode amplitudes have not yet been compared with
Briareus's (ROADMAP, "Benchmark against VisRad").

**Modes.** VisRad's Surface Element Statistics fit Legendre polynomials and real
spherical harmonics by least squares to the element values; the weighting is not
documented. Briareus integrates the harmonics over each element's exact solid angle.
The Legendre coefficients mean the same in both (``A_l`` and ``a_l``, with
``A_2/A_0 = a_2/a_0`` = P2/P0), but the two estimates differ on a coarse mesh; apply
[`visrad_legendre_fit`](@ref) to both codes' element values for a like-for-like
comparison. VisRad's spherical harmonics are normalized so that ``Y_{00} = 1``: they
are Briareus's orthonormal [`ylm_moments`](@ref) times ``\sqrt{4\pi}``. Whether VisRad
includes the Condon–Shortley phase is not stated, so the signs of odd-``M``
coefficients may differ; per-``L`` amplitudes do not.

## Where the codes differ

### The beam away from focus

Near best focus both codes deposit the same super-Gaussian. Away from focus they
differ. Briareus's default model is the geometric-optics one: the best-focus profile
summed over a quadrature of the f-cone's directions, which is the focal profile
convolved with the disk of focal points. It conserves power at any distance and
matches an exact convolution ([Laser spot on a foil](@ref)). VisRad's model is
analytic: a super-Gaussian whose semi-axes ``a(z)``, ``b(z)`` and order ``n(z)`` change
with the distance ``z`` from focus, cut off at an envelope ratio (by default the radius
where the profile falls to 10⁻³ of its peak) and renormalized in each transverse plane.
The two shapes differ in kind: on flat plates up to 15 mm from focus, VisRad's order
rises with ``|z|``, and the f-cone sum lowers it, to about half its focal value at
15 mm.

On a hohlraum wall, 3–8 mm past focus, the f-cone sum gave 10 % rms per element against
VisRad's deposition, with peaks 10 % higher. [`VisRadProfile`](@ref) reproduces VisRad's
model: for the four NIF phase plates it uses the semi-axes and order measured from
VisRad itself on those plates, and for other spots the manual's hyperbolic semi-axes
with a growth factor fitted to one workspace. [`visrad_profile`](@ref) reads its
parameters from a workspace. Neither model is a measurement of the real beams: the
default is the geometric-optics convolution of the focal profile, and
[`VisRadProfile`](@ref) is VisRad's model. Use [`VisRadProfile`](@ref) to reproduce a
VisRad calculation; the measured tables cover only the four NIF plates at VisRad's
default envelope and lens ratio.

Smaller differences of the beam model:

- **Spot rotation.** VisRad's rotation 0 is the quad's horizontal, which is turned from
  the beam's own port frame by 1.6–5.3° depending on the cone; the readers and writers
  convert. The sense of a nonzero rotation (viewed along the beam or from the target) is
  not stated in the manual; the keyword `rotation_sense` selects it.
- **The beam axis.** VisRad's axis runs from the port on the chamber wall through the
  Reference Point; Briareus's is parallel to the port direction. The difference at focus
  is under 1 µm for millimeter-scale aims.
- **Occultation.** VisRad tests occultation with a ray grid (25 × 24 per cone by
  default); Briareus tests each f-cone direction at each point, or, with
  [`VisRadProfile`](@ref), the beam axis through each point. VisRad puts up to
  7 × 10⁻⁴ TW/cm² on LEH washers that Briareus leaves dark.
- **Tilted surfaces.** VisRad renormalizes each beam's total on a tilted surface;
  [`VisRadProfile`](@ref) normalizes each transverse plane, which puts 0.15 % too much
  power on a foil tilted 20°.
- **Laser reflectivity.** VisRad can reflect laser light from a surface; Briareus
  absorbs all of it by default and warns when a workspace object reflects.
- **Unconverted light.** Under VisRad's "Fraction of Total Beam Power" model, the power
  column is the total of the 3ω, 2ω and 1ω cones; the readers take the 3ω share,
  ``P(1 - f_G - f_R)``. Briareus's unconverted-light check
  ([`unconverted_light_check`](@ref)) treats the 2ω and 1ω light separately.

### Radiation on curved surfaces

Briareus's radiation model integrates over the exact surface a component carries: the
cylinder and spheres of [`hohlraum`](@ref), and the cylinders, cones and spheres that
[`visrad_target`](@ref) fits to remeshed CAD (`fit`). On flat facets, the exchange
between neighboring elements is misstated by an amount of order one, which shows in
the flux next to a laser spot; on the curved surface an integrating sphere is uniform
to round-off ([Integrating sphere](@ref)). The research notes do not record VisRad's
view-factor quadrature, so how much of this difference appears in a VisRad calculation
is not known. When rebuilding a VisRad primitive by hand, note that VisRad's grid size
of a sphere, cylinder or cylindrical hohlraum is that of its outer wall, with the inner
wall inward by the wall thickness.

### Closure and visibility

A discretized view-factor matrix closes only to its quadrature error, and at high
albedo the cavity multiplies that error. For a target declared closed (apertures over
every opening), Briareus rescales the view factors so that every row sums to one while
reciprocity holds ([`enforce_closure!`](@ref)). [`hohlraum`](@ref) declares its targets
closed; [`visrad_target`](@ref) does not, and takes `closed = true` once the apertures
close the cavity ([`open_fraction`](@ref) checks for leaks). Whether VisRad applies a
comparable correction is not recorded in the research notes.

Visibility between two elements is determined from their centroids by default
(`visibility = :centroid`); `:sampled` traces a grid of rays where the centroid rays
and probe rays disagree, and is the default for CAD and VisRad targets, whose LEH
washers and rims partly shadow elements.

### Element values

VisRad's direct laser power source is an element average: centroid values differ from
it by up to 13 % next to a spot's peak. Briareus deposits the integral over each
element. [`visrad_deposition`](@ref) multiplies VisRad's values by the element areas, so
the two codes compare in TW per element.

### What VisRad has that Briareus does not

From the import warnings and ROADMAP ("VisRad import and CAD: what is skipped", "GUI
requests"):

- STEP CAD import (Briareus reads STL and the CAD bundled in workspaces), and the NIF
  CMT XML.
- Most of VisRad's primitives: cylinders, cones, disks, boxes, polygons, halfraums,
  rugby and slotted hohlraums, surfaces of revolution, multi-element objects, wall
  thickness and partial spheres. Briareus has its own mesh builders for most of these
  shapes, but [`read_visrad_workspace`](@ref) skips them with a warning.
- From workspaces: unconverted 2ω and 1ω beams, time-dependent beam powers and
  energy-normalized beams, DPR (split OMEGA beams), pointing toward a diagnostic port or
  in target-mount frames, and the time-dependent tables of albedo, conversion, lengths
  and power sources. Briareus has its own time-dependent models (above), but does not
  import VisRad's.
- Laser reflection, and power self-emission.
- The blend toward a square beam far from focus (zero for NIF's beams).

What Briareus adds (gradient-based optimization, laser–plasma corrections, the moving
wall and closing holes, a synthetic Dante, worst-case clearance and NIF's
target-configuration rules) is summarized on the home page.

## Porting a VisRad study

The code below is a sketch of the usual sequence, not run by the documentation build:
it needs a workspace. File names and component names are placeholders.

**1. Read the workspace.** With a template beam set, the workspace's pointing, power and
spots are applied to the NIF beams by name, which keeps their cones as `group`:

```julia
using Briareus

ws = read_visrad_workspace("study.vrw"; template = nif_beams())
ws.objects                  # components: name, kind, folder, albedo, conversion, flags
beams = ws.beams            # 3ω power, aim, defocus and spot per beam
profile = visrad_profile(ws)   # VisRad's beam model, with this workspace's envelopes
```

**2. Check the deposition against VisRad's own.** On a workspace that VisRad has
computed, a target on VisRad's own elements (no remeshing) takes VisRad's per-element
deposition. [`visrad_deposition`](@ref) maps only onto components on those elements, so
this target serves the comparison only; step 3 builds the radiation target:

```julia
t0 = visrad_target(ws)
dep = deposit(DepositionOperator(t0, beams; profile), beams)
P = visrad_deposition(ws, t0)          # TW per element, from VisRad's table
wall = elements(t0, :Hohlraum_wall)    # the wall object's name in the workspace
sum(dep.element[wall]), sum(P[wall])
```

**3. Build the radiation target.** CAD can hold far more triangles than a view-factor
matrix should. Select the parts that bound the cavity, repair their normals, remesh
them and keep what bounds the free space around a point in the gas:

```julia
t = visrad_target(ws; select = o -> o.folder == "Hohlraum",
                  orient = true, remesh = 0.5, inside = (0.0, 0.0, 2.0),
                  materials = Dict("Hohlraum_wall" => nif_gold()))
open_fraction(t).fraction    # pass closed = true above only if apertures close the cavity
model = RadiationModel(t)
```

**4. Solve and analyze.**

```julia
op = DepositionOperator(t, beams; profile)    # or the default f-cone model
sol = solve(model, deposit(op, beams))
capsule = elements(t, :Capsule)
a = legendre_moments(sol.H, t, capsule)
a[3] / a[1]                                   # P2/P0
```

**5. Time-dependent powers.** VisRad's power file becomes a pulse per beam, keyed by
beam name, and the drive history evaluates the beams at each time (every beam needs a
column in the file):

```julia
pulses = read_visrad_powers("powers.csv"; template = beams)
beams_at(time) = with_pulses(beams, pulses, time; by = :beam)
history = drive_history(model, op, beams_at, 0.0:0.25:9.0; observe = :Capsule)
history.moments[3, :] ./ history.moments[1, :]   # P2/P0 against time
```

**6. Optimize and write back.** Re-point and re-balance by cone, then write the design
into a copy of the original workspace, whose components, spots and settings stay
VisRad's, and as an SST file:

```julia
pointing = GroupPointing(beams; fields = (:z, :power), relative = true)
scenario = Scenario(t, pointing; observe = :Capsule, deposition_kw = (; profile))
x0 = initial_parameters(pointing)
drive = mean_flux(evaluate(scenario, x0))
objective = symmetry_objective(; drive)   # the direct penalty is on `observe`
lower, upper = bounds(pointing, x0; power = 0.3)
result = optimize_pointing(scenario, objective; x0, lower, upper)
design = apply(pointing, result.x)
write_visrad_workspace("study_repointed.vrw", "study.vrw", design)
write_visrad_sst("study_repointed_sst.csv", design)
```

**7. Benchmark a case in both codes.** For a target built in Briareus (`target`, here a
[`hohlraum`](@ref) with a component named `capsule`), the kit holds
what VisRad needs to rebuild the case and Briareus's results; VisRad's element and power
summary exports are then compared quantity by quantity:

```julia
kit = write_visrad_kit("kit", target, design; observe = :capsule)
# ... build and run the case in VisRad from kit/README.txt, export its tables ...
compare_visrad("kit"; elements = "capsule_elements.txt",
               power_summary = "power_summary.psd",
               names = Dict("hohlraum_barrel" => "wall"), io = stdout)
```

The research notes ("VisRad interchange formats" and the sections after it) record how
each format detail was established, its confidence, and the items still to confirm on
real files ("To verify with real files at work").
