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

## Next
1. **Sombrero and backlighter geometry.** Model the LEH shield as a component (a
   dimpled annulus), and the backlighter foil's placement relative to the diagnostic
   line of sight. The line-of-sight and clearance checks in the radiography tutorial
   then decide which bundles can drive the foil. This needs real dimensions.
2. **Diagnostic windows: design studies.** `hohlraum(...; windows)` is in place (closed or
   open windows, conforming mesh, exact curved areas). Next: a tutorial on re-tuning
   pointing around a window set (mode 1, m = 2 and m = 4 at the waist, P2), and a
   time-dependent burn-through model (closed → open). First look (two open 2 × 2 mm
   windows at φ = 78°/258°, NIF quads): −18% drive, c(2,2)/c₀₀ = +0.061, ε₄ 0.6 → 3.1%.
   Restoring the no-window shape by quad power needs −22 to 57 TW (infeasible);
   re-pointing (±1.5 mm, ±20% power, keep-out) only reduces the m = 2 departure from
   0.061 to ~0.04–0.05. Window size looks like the main lever; needs the user's window
   dimensions and goals before a tutorial.
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
