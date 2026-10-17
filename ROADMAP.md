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
     spot-rotation sense, SST preamble, element-export layout);
   - the NIF CMT XML schema.
4. **Benchmark against VisRad.** Run `write_visrad_kit` cases at work and bring back the
   outputs for `compare_visrad`. This validates absolute drive levels and mode
   amplitudes.
5. **GUI requests:** foil center and normal stored on the component; session
   save/load; x-ray images in the window. (Cached radiosity factorization, clearance and
   optimizer progress/cancel hooks are done.)
6. **Time-dependent view factor.** Wall, LEH and spot motion from rad-hydro or a
   parametric model; P2(t) swing; optimization of cone fraction versus time.
7. **Shot-data calibration.** Fit effective parameters (cone multipliers, CBET
   multiplier/clamp, albedo) to measured P2/P4 and shock symmetry, with uncertainty.
   Reproduce P2-vs-Δλ curves.
8. **Robust optimization and UQ.** Monte Carlo or Jacobian propagation of 50 µm
   pointing and power imbalance; optimize the mean plus kσ; sensitivity matrices;
   mode-1 and m = 4 budgets and inverse solves (MacGowan 2021).
9. **Backlighter physics.** CE(Z, I, τ) models with uncertainty, and the photometric
   chain to the detector.
10. **Spectral radiosity.** Multigroup emission, detector responses along Dante and XRD
    lines of sight, and brightness temperatures.
11. **Scale.** Multi-ray partial visibility for view factors, hierarchical or iterative
    radiosity, Float32 storage for more than 10⁴ elements, and GPU kernels.
12. **Polar direct drive mode**, benchmarked against published GA optimizers.
