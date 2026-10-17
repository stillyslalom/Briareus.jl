# Roadmap

## Done (v0.1, 2026-10-15)
Core pipeline: geometry/BVH, NIF + OMEGA beam sets, differentiable deposition operator +
reference ray trace, occluded view factors with exact near-field kernel, gray radiosity
with constant/tabulated/wall-loss albedo, flux response operator, Legendre/Y_lm
analysis with exact facet integration, orthographic rendering, quasi-static drive
history, grouped pointing parameterization + L-BFGS-B optimization, examples, 676 tests.

## Done (2026-10-16)
Shot-review layer:
- worst-case clearance over a linear error stack, with intensity-contour envelopes, spot
  wings and pass/review/fail tiers;
- the unconverted-light check;
- beam-table I/O;
- backscatter, calibrated transfer, and a reduced kinetic CBET model with calibration.

## Done (2026-10-16, overnight)
- **Facility data:** the LMF per-beam table, NIF quad names and bundles, and the
  measured CPP fits.
- **Radiation on true curved surfaces:** sphere, cylinder and disk sectors, with
  adaptive near-pair quadrature. View factors are 19× faster and use 165× less memory,
  with a SAH ray tree for every target.
- **Real spherical-harmonic analysis:** `mode_amplitude`, the influence matrix and
  `asymmetry_objective`.
- **Deposition operator rebuilt:**
  - It sums the best-focus profile over the f-cone directions, with per-direction
    incidence factors; the moment-matched profile had the wrong shape far from focus.
  - It refines triangles where each beam can land.
  - Per-beam power is conserved to 1e-4 unclipped and 2e-3 when clipped, on fine and
    coarse meshes.
- **Shot review:**
  - `ApertureKeepOut`, a differentiable worst-case LEH clearance for optimization;
  - `talis_check`, NIF's published target-configuration rules.
- **VisRad interchange:** SST pointing, Laser Beam Summary and power files, plus the
  benchmark kit.
- **Docs:** a Documenter site with tutorials that double as validation (laser spot,
  integrating sphere, freeing a bundle for radiography).
- **`BriareusGUI`:** a desktop window built on GLMakie + QML.
- **Hohlraum diagnostic windows** (`WallWindow`), exact cylinder-patch areas, cached
  radiosity factorizations (repeated solves 45× faster), and progress/cancel hooks on
  `clearance` and `optimize_pointing`.

## Done (2026-10-16, windows)
- **Circular windows and gaps:** `WallWindow(; diameter, gap)`. Disks are drawn on
  the unrolled wall and re-meshed inside a box stitched to the conforming grid. An
  equal-area edge polygon gives exact window and gap areas.
- **View-factor visibility on curved surfaces:** facets of a cylinder or sphere no
  longer occlude rays between points inside it. The centroid facing test is relaxed
  by the spread of true normals across a facet. A tilted triangle on the wall had been
  culling real exchanges, and capsule rows now close 10× better.
- **Thin-foil burnthrough:** `BurnthroughAlbedo` and `HeatFront` use Hammer & Rosen's
  subsonic gold heat front, integrated over each element's own T_r(t) in
  `drive_history`. `energy_balance` reports `transmitted`, and `drive_history`
  reports `window_loss`.
- **Tutorial "Diagnostic windows in the hohlraum wall":** the front is checked against
  the exact ramp solutions (2.2% and 0.7%) and the Nova foil-burnthrough data. It runs
  1.7, 2.5 and 3.6 µm windows with and without a 50 µm gap.
- **Adaptive drive histories:** `drive_history(...; transition, flux_change)` inserts
  steps through burnthrough transitions (planned from predicted burnthrough times) or
  wherever the observed flux shape moves too far. Each step's albedo iteration starts
  from the previous step's temperatures (`solve(...; guess)`), about a third fewer
  factorizations.
- **Post-burnthrough transmission (default):** a burned-through foil transmits
  `T = 4/(3τ + 2)` (gray planar Eddington with a Marshak boundary; Lester et al., Phys.
  Rev. E 114, L023201 (2026)), τ = 5 at 2.5 µm of gold and linear in areal mass, and
  re-emits the rest with the wall's albedo: T = 0.33, 0.235, 0.17 at 1.7, 2.5, 3.6 µm.
  `eddington_transmission`, `optical_depth`, `sol.transmission`; `after = 0` keeps the
  open-hole bound. In the tutorial the windows' flux-weighted ΔP2 drops from +4.1–5.3%
  (open hole) to +0.7–1.7%, and the capsule-energy cost from 12–15% to 2.4–5.5%.
- **Threaded deposition:** `deposit` hands out work units (element blocks, back-side
  pieces) dynamically; each unit owns its outputs, so results are bitwise independent
  of the thread count. Operator construction and the ray trace are balanced the same
  way.
- **Back-side sampling:** the outside of a target (clipped power) is sampled for its
  total, with Gauss rules graded toward open edges (LEH rims). Clipped power at 1e-4 of
  a beam had been 10–60% low with spot-resolving sub-triangles; it now matches
  converged ray traces to 1%, with 44% fewer operator points on the NIF hohlraum.
  Element values are unchanged bit for bit.

## Done (2026-10-17, GUI as a design tool)
- **BriareusGUI as a CAD-style tool:** the step rail is gone. Model tree (components,
  beams by cone/quad/beam, analyses with their state), Properties for the selection,
  result tabs, menus and a toolbar; any operation at any time, stale results flagged;
  coloring by a drive field computes the drive on request (as VisRad does).
- **Viewer:** own camera (`CadCamera`): zoom toward the pointer, orbit about the
  selection or the picked point, pan the grabbed point, double-click fit, standard and
  along-beam views, perspective/orthographic, near/far planes from the bounding
  sphere (no clipping), ray-cast picking, an orientation triad. Beams as rays to where
  they land, or as envelope volumes cut at the target.
- **Free parts** (sphere, disk, washer, rectangle, cylinder, cone) placed with
  VisRad's polar/azimuth/rotation angles, with their own material, sidedness,
  aperture flag and element size. Round diagnostic windows by default.
- **Beam groups** (VisRad's clusters, kept apart from the cone `group` the LPI models
  use): any set of beams, Ctrl-selected in the tree or the view, moved to a named
  group and pointed together, relative to the ports or at a common aim point; the
  optimizer's per-group parameters follow them. NIF quads split into their four
  beams and merge back. The tree carries the component and group colors (no color
  key on the canvas).
- **Capsule keep-out** (`SphereKeepOut`): the optimizer had moved inner cones onto the
  capsule (0.75% of their power) because the `direct` penalty is quadratic in a tiny
  fraction; the keep-out holds the worst-case envelope off it (on by default in the
  GUI). `optimize_pointing` also caps objective evaluations (`evaluations`), which a
  GUI run had overrun with 100 iterations of long line searches.
- **Mesh work-precision study** (`examples/mesh_convergence.jl`): element-size meshing
  (`hohlraum(; h, capsule_h)`), GUI presets coarse/standard/fine at ~0.5/0.1/0.02
  percentage points in P2 and P4.

## Done (2026-10-17, VisRad workspaces)
- **`read_visrad_workspace`** reads VisRad 20/21 `.vrw` workspaces, which are keyed
  text, not the JSON the manual describes. It reads the target components, including
  bundled CAD as generalized meshes (cm to mm, placed in the chamber frame through
  reference frames), and the beams with their pointing, 3ω power, f-number and spots.
  NIF CPP radii are recognized as HWHM by their numbers. `visrad_target` builds a
  `Target`, meshing spheres, rectangles and cylindrical hohlraums on VisRad's own grid
  (element counts checked against VisRad's power-source table). Checked on a real NIF
  CAD workspace: ports match `nif_beams` to 0.006°, and a ray trace of its 184 beams
  through all 210k CAD elements lands them in the hohlraum. CAD normals there are
  inconsistent (the inner wall faces out), hence `two_sided`.

## Next
1. **Sombrero and backlighter geometry.** Model the LEH shield as a component (a
   dimpled annulus), and the backlighter foil's placement relative to the diagnostic
   line of sight. The line-of-sight and clearance checks in the radiography tutorial
   then decide which bundles can drive the foil. This needs real dimensions.
2. **Diagnostic windows: design studies.** Circular windows, gaps, burnthrough and the
   gray post-burnthrough transmission are in place (see the tutorial). Two 2 mm windows
   of 1.7–3.6 µm gold give a flux-weighted ΔP2 of +0.7 to +1.7% and cost 2.4–5.5% of
   capsule energy (open-hole bound, `after = 0`: +4 to +5% and 12–15%).
   - Refine the transmission: λ_Au is held at a 300–350 eV evaluation. A
     temperature- and density-dependent Rosseland mean free path (or a tabulated
     T(t)) would let it follow the drive, and a backing could absorb and re-emit part
     of what the gold transmits.
   - Model gap closure in time, e.g. a `TabulatedAlbedo` or a closing-gap model.
   - Re-tune pointing around the window pair: m = 2 and P2 at the waist, versus time.
   - Earlier look (two open 2 × 2 mm windows, NIF quads): power rebalance is
     infeasible, and re-pointing only trims the m = 2 departure.
3. **Confirm against facility practice:**
   - the TaLIS 0.5° tolerance (pivot, and whether it replaces or adds to the error set);
   - the 1ω footprint side and shape;
   - which quads count as open for the 20° rule;
   - the VisRad items listed in the research notes (summary-CSV columns after Power,
     spot-rotation sense, SST preamble, element-export layout; in workspaces, the
     Rotation/Spin angle mapping, the cylindrical position order, primitive grid
     spacing and the spot-size definition of custom profiles: one small workspace
     with a rotated cylinder and a rotated CAD part would settle most of them);
   - the NIF CMT XML schema.
4. **Benchmark against VisRad.** Run `write_visrad_kit` cases at work and bring back the
   outputs for `compare_visrad`. This validates absolute drive levels and mode
   amplitudes.
5. **VisRad workspace import: what `read_visrad_workspace` skips.** Each is warned
   about or refused today. Most need a sample workspace that uses the feature, since
   the layout of the extra data is not seen yet.
   - Unconverted light: beams with 2ω/1ω on (the green/red power columns) as
     `nif_unconverted` colors or extra beams.
   - Time-dependent beam powers (power mode ≠ 0, `laser power table values` with
     N > 0; VisRad's time unit is seconds), energy-normalized beams, and named pulse
     shapes, as `PulseShape`s.
   - DPR beams (split into two cones, as VisRad does for OMEGA).
   - Pointing method 1 (distance toward a diagnostic port, from `[Diagnostic Port
     Data]`), and target-mount or diagnostic-port reference frames for objects and beams.
   - Objects positioned by a key point other than the centroid (cylinder top/bottom,
     box and rectangle corners).
   - Time-dependent albedo, x-ray conversion, lengths and power sources (the ten
     `Time-Dependent …` tables), laser reflectivity, and power self-emission.
   - Wall thickness (inner and side walls), partial spheres, and the other primitive
     types: cylinder, cone, disk, box, polygon, halfraum, rugby and slotted hohlraums,
     surfaces of revolution, multi-element objects.
   - CAD meshes: consistent orientation (repair inward/outward faces, e.g. by ray
     casting from inside closed parts) instead of `two_sided`, and coplanar merging or
     decimation so that selected CAD parts are small enough for view factors.
6. **GUI requests:** session save/load (a workspace file); undo/redo; x-ray images
   in the window; reference frames for parts (attach to a component or a port axis, as
   VisRad does); clipping volumes and STEP/STL import (CAD already arrives through VisRad
   workspaces, `read_visrad_workspace`; a GUI loader and decimation for view factors
   would make it usable there); a second simultaneous view
   (the capsule map beside the 3D view).
   Adaptive mesh refinement of the radiation mesh where deposition varies fastest
   (the uniform fine mesh is converged but costs ~16k elements).
7. **Time-dependent view factor.** Wall, LEH and spot motion from rad-hydro or a
   parametric model; P2(t) swing; optimization of cone fraction versus time.
8. **Shot-data calibration.** Fit effective parameters (cone multipliers, CBET
   multiplier/clamp, albedo) to measured P2/P4 and shock symmetry, with uncertainty.
   Reproduce P2-vs-Δλ curves.
9. **Robust optimization and UQ.** Monte Carlo or Jacobian propagation of 50 µm
   pointing and power imbalance; optimize the mean plus kσ; sensitivity matrices;
   mode-1 and m = 4 budgets and inverse solves (MacGowan 2021).
10. **Backlighter physics.** CE(Z, I, τ) models with uncertainty, and the photometric
   chain to the detector.
11. **Spectral radiosity.** Multigroup emission, detector responses along Dante and XRD
    lines of sight, and brightness temperatures.
12. **Scale.** Multi-ray partial visibility for view factors, hierarchical or iterative
    radiosity, Float32 storage for more than 10⁴ elements, and GPU kernels.
13. **Polar direct drive mode**, benchmarked against published GA optimizers.
