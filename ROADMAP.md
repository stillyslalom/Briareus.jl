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

## Next (from the SOTA survey, reference/research-notes.md)
1. **VisRad bridge.** There is no public format spec. Collect sanitized sample files
   (workspace JSON, pointing tables, capsule-flux and clearance-report exports) and
   write readers and writers. Add a "benchmark kit" that exports standard cases to run
   in VisRad at work and compares the outputs back here.
2. **Time-dependent view factor.** Wall, LEH and spot motion from rad-hydro or a
   parametric model; P2(t) swing; cone fraction vs. time optimization.
3. **Shot-data calibration.** Fit effective parameters (cone multipliers, CBET
   multiplier/clamp, albedo) to measured P2/P4 and shock symmetry, with uncertainty.
   Reproduce P2-vs-Δλ curves.
4. **Robust optimization and UQ.** Monte Carlo or Jacobian propagation of 50 µm
   pointing and power imbalance; optimize mean + kσ; sensitivity matrices. Also mode-1
   and m=4 budgets and inverse solves (MacGowan 2021).
5. **Independent 23.5°/30° Δλ and per-quad azimuthal knobs.**
6. **Polar direct drive mode**, benchmarked against published GA optimizers.
7. **Backlighter module.** CE(Z, I, τ) models with uncertainty; photometric chain to
   the detector.
8. **Visualization extension (Makie weakdep).**
9. **Docs** (Documenter, Diátaxis).

## Earlier backlog
1. **Benchmark against VisRad.** Reproduce published VisRad cases (MacFarlane 2003
   OMEGA hohlraum; NIF capsule P2/P4 for a known pointing) and user-supplied VisRad runs.
2. **Visualization extension (Makie weakdep).** 3D target + beam view, unwrapped wall
   and Mollweide capsule flux maps, deposition footprints, time histories.
3. **Beam data I/O.** Read real port tables (NIF per-beam, OMEGA chamber frame with
   P6/P7 orientation, OMEGA EP), VisRad-style setup files; save/load scenarios.
4. **Robustness.** Monte Carlo pointing/power-imbalance errors (e.g. ±50 µm rms),
   sensitivity reports (∂P_l/∂parameter from the AD gradients), multi-start optimization.
5. **Time-dependent geometry.** Wall motion and LEH closure (effective radii vs time),
   capsule convergence; time-dependent pointing/pulse objectives (P2(t) swing).
6. **Effective laser–plasma corrections.** Per-cone backscatter, CBET-like power transfer
   between cones as a parameterized model, 1ω/2ω unconverted light.
7. **Spectral radiosity.** Multigroup Planckian emission, detector responses (Dante /
   XRD lines of sight from facility diagnostic ports), brightness temperatures.
8. **Scale.** Multi-ray partial visibility for view factors, hierarchical/iterative
   radiosity and Float32 storage for >10⁴ elements, GPU kernels.
9. **Docs.** Documenter site (Diátaxis: tutorials from the examples, how-tos, explanation
   of the deposition/view-factor models, API reference).
