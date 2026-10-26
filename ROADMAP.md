# Roadmap

"Done" is a terse history, oldest first; the numbers behind each entry are in
`reference/research-notes.md`, under the section named. "Next" is the open work: each
item says where it stands, what remains and what data would settle it.

## Done

- **v0.1 (2026-10-15).** Geometry/BVH, NIF and OMEGA beams, differentiable deposition
  operator and reference ray trace, occluded view factors, gray radiosity, flux response,
  Legendre/Y_lm analysis, rendering, quasi-static drive history, pointing optimization.
- **Shot review (2026-10-16).** Worst-case `clearance` with envelopes, wings and tiers;
  the unconverted-light check; beam-table I/O; backscatter, calibrated transfer and the
  reduced kinetic CBET model ("Beam clearance methodology", "Cross-beam energy transfer").
- **Facility data, curved surfaces, deposition (2026-10-16, overnight).** LMF beam table,
  quads, CPP fits; radiation on true spheres, cylinders and disks; `mode_amplitude`;
  the deposition operator rebuilt on the f-cone sum; `ApertureKeepOut`, `talis_check`;
  VisRad SST/summary/power files and the kit; the docs site; `BriareusGUI`; `WallWindow`
  ("Later search (2026-10-16)", "VisRad interchange formats").
- **Diagnostic windows (2026-10-16).** Round windows and gaps, visibility on curved
  surfaces, Hammer & Rosen burnthrough, Lester's post-burnthrough transmission, adaptive
  histories, threaded deposition, graded back-side sampling ("Gold wall loss and
  thin-foil burnthrough").
- **GUI as a design tool (2026-10-17).** Model tree, CAD-style viewer, free parts, beam
  groups, `SphereKeepOut`, an evaluation cap; mesh presets ("Hohlraum mesh
  work-precision").
- **VisRad workspaces (2026-10-17).** `read_visrad_workspace`, `visrad_target` ("VisRad
  workspace files").
- **Performance (2026-10-17; 2026-10-23).** First objective call 14 s → 4 ms, gradient
  31 → 8 ms, `clearance` 6.0 → 1.8 s; concurrent gradient chunks (~2×, bitwise equal to
  serial).
- **Windows and CAD repair (2026-10-17).** `RosselandTransmission`, `GapClosure`,
  `Scenario(...; albedo)`, restoring the drive around windows; `orient_mesh`,
  `coarsen_mesh` ("Restoring the drive around a window pair").
- **GUI sessions and CAD parts (2026-10-17).** TOML sessions, undo/redo, STL parts, part
  frames, x-ray images in the window.
- **Beam errors (2026-10-17).** `mode_sensitivity`, `mode_uncertainty`, `sample_modes`,
  robust objectives ("Beam-error budgets").
- **Design workflow (2026-10-17).** `retune_drive`, pulse tables, `write_results`, VisRad
  import/export in the window, frames ("Retuning to a reference drive").
- **VisRad deposition (2026-10-17).** The red/green power model; a VisRad deposition table
  matched; `VisRadProfile` ("VisRad workspace files": "The faint wings").
- **CAD and STL (2026-10-18).** `write_stl`/`read_stl_solids`, `remesh_mesh`,
  `enclosure`, straight cuts, `fit_surfaces`/`curved_components` with cones,
  `cad_target`, the GUI on CAD workspaces, the CAD tutorial ("CAD import, remeshing and
  surface fits").
- **VisRad deposition parity (2026-10-21).** Defocus growth γ/ρ (γ = 0.67), measured NIF
  plates, spot rotation, per-beam envelopes: wall 0.9 % rms against VisRad ("VisRad
  deposition parity: hohlraum runs").
- **2026-10-22.** Shared beam errors (`group_pointing`, `group_power`, `positioning`);
  `infer_errors`; the beam model in the GUI; VisRad results read back
  (`visrad_power_sources`, `visrad_deposition`); `icosphere_mesh`; `:sampled` visibility;
  `montecarlo_incident`; rounded hohlraum corners ("Compact capsule mesh", "Partial
  visibility", "Monte Carlo radiosity reference for partial visibility").
- **Hohlraum dynamics, implosion metrics, N221204 (2026-10-23).** `src/implosion.jl`;
  `WallMotion`/`Keyframes`; closing LEH lips; `BubbleAbsorption` and `lpi_deposit`; the
  mock HYBRID-E pointing; `TwoZoneCBET` calibrated to PRE 109 Fig. 9; `examples/n221204.jl`;
  GMRES albedo iterations (history steps 32 → 3.6 s); VisRad kits built ("Implosion
  metrics", "Moving wall on keyframes", "Closing entrance holes", "Inner-beam bubble
  absorption", "Mock HYBRID-E pointing", "Two-zone CBET calibrated to PRE 109 Fig. 9",
  "N221204 time-dependent benchmark").
- **2026-10-24, targets and backlighters.** `Sombrero` LEH shields; `backlighter_foil`,
  `line_of_sight`; He-α CE, attenuation and `backlighter_signal` ("LEH shield
  (sombrero)", "Backlighter foil on the line of sight", "Backlighter conversion
  efficiency and the photometric chain").
- **2026-10-24, dynamics and CBET calibrated to published data.** Wall motion, bubble and
  collision refitted; `LipClosure`; joint multi-shot CBET (`hybride_cbet`); Young 2020
  and Nora 2023 checks ("Hohlraum dynamics calibrated to published data", "Computed LEH
  closure", "Joint multi-shot CBET calibration", "Young 2020 validation").
- **2026-10-24, energetics.** `du_wall_loss`/`uranium`, `Conversion`/`nif_conversion`,
  `CapsuleAlbedo`, `energy_budget`; wall loss from each element's temperature history
  (`wall_history`) and the measured `n221204_pulse()`; wall re-emission
  (`wall_release`); `capsule_shrink` ("Wall material, conversion and capsule
  energetics", "N221204 late-time drive", "Wall re-emission and the capsule radius";
  values since superseded in part, see 2026-10-25).
- **2026-10-24, spectra and Dante.** `spectral_solve`, the laser M-band, a synthetic
  Dante, Dante 1 channels from their hardware, the unfold, Sentinel views, thermal laser
  spots ("Spectral radiosity, M-band and a synthetic Dante", "Dante channel responses
  and the unfold", "Laser-spot spectrum").
- **2026-10-24, scale.** `Float32Matrix` and factorization-free solves above 10⁴
  elements; `HMatrix` above 1.5 × 10⁴ ("Large radiation models", "Hierarchical view
  factors").
- **2026-10-24, CAD.** Fitting on by default with `remesh`, parts of fitted CAD
  (`Component.part`, `parts`), `minfraction = 1e-3`; tori; covers for faces of none of
  the kinds; `:sampled` visibility on folded CAD elements fixed ("Fitting CAD by
  default", "Parts of fitted CAD", "Tori", "Sampled visibility on folded elements").
- **2026-10-24, GUI time.** Pulses, drive histories, the retune and a time cursor in the
  window (`BriareusGUI/README.md`).
- **2026-10-24, CompactLES stage 0 started.** `CapsuleDrive`/`capsule_drive` with every
  m ("Capsule drive hand-over and CompactLES stage 0").
- **2026-10-25, CBET null.** NC24's 1–1.55 Å is the round-hot-spot Δλ, not the null;
  `hybride_cbet_iraum()` ("Joint multi-shot CBET calibration: the null").
- **2026-10-25, scale.** H-matrix clusters by surface, exact ball masks, `tol` 1e-6, the
  GUI budget; one ray packet per triangle in the deposition operator; slivers cut across
  their length ("The H-matrix's N^1.5 remainder", "Deposition-operator rays and sliver
  view factors").
- **2026-10-25, energetics defaults (current).** Zylstra 2022 LEH pair; published
  energetics as the defaults, DU factor 0.836; late-time drive literature and an
  ablating-wall reference; the P2 swing's published inputs (K24 capsule radius, uniform
  backscatter cut); then the final set: `gold()`/`uranium()` at 0.75 and the
  drive-history presets `nif_gold()`/`nif_uranium()`, `release` 1.0, Young's capsule
  albedo on the K24 schedule, `CapsuleAlbedo(; outside)` ("Zylstra 2022 LEH pair,
  time-dependent", "Energetics defaults and the capsule share", "Late-time drive: the
  wall, the plasma and Dante", "The P2 swing: capsule radius, backscatter and the
  collision", "Capsule absorption, conversion and the ablating-wall release").
- **2026-10-25, backscatter histories.** N210808's per-cone SBS (Kemp 2026) as
  `backscatter=history`, and the late-time refraction as a bound ("Inner-cone SBS time
  histories and the late-time refraction").
- **2026-10-25, compile latency and tests.** `drive_history`/`optimize_pointing` compiled
  once for every closure, a wider workload ("Compile latency"); shared test fixtures,
  three CI shards (`BRIAREUS_TEST_SHARD`), `BRIAREUS_LONG_TESTS`. The workload then
  cut back to the common scripted session (CAD, TALIS, spectra/Dante, CBET and the VisRad
  kit compile at first use): package image 138 → 96 s.
- **2026-10-25, closure on rims, windows and CAD.** Balancing checked against continuum
  Monte Carlos on bored LEHs, a waist gap, windows and a remeshed CAD cavity;
  `open_fraction`; `cad_target`/`visrad_target(...; closed)` ("Closure on bored,
  windowed and CAD cavities"). The compact capsule's raw rows at capsule_h = h/3
  explained ("Compact capsule rows at capsule_h = h/3").

## Decisions from the CAD work (2026-10-18; the leans confirmed by the user 2026-10-22)

Parts may be split for computation, but the GUI shows and picks the original parts.
1. `fit = true` with `remesh`: done 2026-10-24 (above).
2. `minfraction = 1e-3`: done 2026-10-24; `maxradius` stays off.
3. Cones in `fit_surfaces`' default kinds: kept.
4. CAD capsule tents and LEH windows seal the capsule in a pocket: the GUI opens it
   (tent as a thin film); the library stays explicit. Open: a transmitting thin-film
   option (Next 6).
5. `enclosure` threshold 0.05 (0.01 costs 1.4×, return estimate 0.8 → 0.12 mm² of
   14.5 mm² open): kept.
6. GUI: far CAD parts drawn only on import; per-beam parameters for the foil objective;
   no cavity for non-CAD designs; the 30k-triangle background build; the colorbar flash
   at window open; the whole-part fit guard (`_on_one_surface`): all kept. Partial
   curved fits of CAD parts in the GUI wait for a workspace that benefits (Next 7).
7. `coarsen` and `remesh` both kept; the docs recommend `remesh`.
8. Remeshed planar strips as 2–4-facet rectangles on the primitive path: accepted.
9. `hohlraum()` uses the compact capsule mesh (GUI presets capsule_h 0.3/0.16/0.09);
   `Target` carries a visibility rule, `:sampled` for `cad_target`/`visrad_target`,
   `:centroid` for `hohlraum()` and hand-built targets ("Default capsule mesh and
   visibility"). Balancing's "own bias" on flat elements (+0.06 to +0.48 pp in a₂/a₀) was
   cracked facets (item 11).
10. `orient_mesh` faces a hollow solid's inner surface into the material: kept,
    documented.
11. `closed = true` stays the `hohlraum()` default (curved path, against a continuum Monte
    Carlo at every preset: "Balancing on the curved path ..."). Done 2026-10-25: with
    windows, gaps, lips and sombreros it stays; bored LEHs, gaps and CAD cavities closed
    by apertures balance to the reference loss within 0.2 % (raw up to 0.8 %, 2–3 % with
    centroid rays), a₂/a₀ within the raw factors' ± 0.08 pp with `:sampled`. Closure is
    declared, not inferred: hand-built targets pass `closed = true`, CAD ones
    `cad_target(...; closed = true)` (new keyword, default false); `open_fraction` shows
    where facets leak, but cracks (24-gon disks in 20-gon bores) and small openings give
    the same fraction. Flat facets must be watertight to balance ("Closure on bored,
    windowed and CAD cavities"). The GUI's CAD cavity passes `closed` once the user's
    apertures close it (`open_fraction` ≤ 1e-2, `CAVITY_MAX_LEAK`): done 2026-10-26.
12. Compact capsule at capsule_h = h/3: done 2026-10-25, documented. Raw capsule rows
    −3.5e-4 at h 0.25 (−9e-5 at 0.4) are the far-pair one-point rule on the capsule–wall
    pairs (the wall's size sets the tier; `mid = 24` removes it), not the geodesic
    triangles (lat-long −2.3e-4); the source-weighted loss is the centroid rule's capsule
    shadow, 0.4–0.9 % too large on any capsule mesh. Balanced tutorial results are within
    2.2σ (1.3e-4 of the source) of the exact walk on both meshes ("Compact capsule rows
    at capsule_h = h/3").
13. `:sampled` visibility on remeshed CAD: fixed 2026-10-24 (above).

## Next

1. **VisRad deposition parity away from focus.** *State:* `VisRadProfile` with
   `visrad_profile(ws)` matches VisRad near focus and on the NIF hohlraum runs (wall
   0.9 % rms); on N230729 at 0.1 mm it is 1.16 % rms per wall element, almost all on
   the 44.5°/50° shoulders (−3 to −16 % at r/ρ 0.6–0.9), and per-cone refits to that run
   degrade the others. Kits for six cases are built (thor `visrad/kits_2026-10-23`).
   *Remains:* the law behind the tabulated a(z), b(z) (plates other than the four NIF
   ones); δ2 0.012 pp high on bu34_repointed (the wall); totals (wall 0.07 % low, tilted
   foil 0.37 % high); the foil wings; the envelope ratio's definition (effective cut-off
   ~1.64 1/e radii against the file's 1.53645); whether VisRad renormalizes on tilted
   surfaces (per-plane normalization gives 1.1e-4 too much at ~40°); occultation
   (VisRad's 25 × 24 ray grid against Briareus's beam axis per point). *Data:* VisRad
   runs of the kits, each settling one of these:
   - the intersection-point export (File ▸ Export, with "write beam spatial profile
     parameters" in Preferences ▸ Lasers): a(z), b(z), n(z), f_SQ(z), ΔL_W/H(z), M(z) at
     each hit, on a hohlraum (the growth: a(z) ~6 % above a₀ at 4.5 mm for 50°); one cone
     with a non-default envelope or lens ratio says which sets the growth;
   - a foil normal to the beams (focal profile apart from the defocus law); a focus sweep
     (±1, ±3, ±10 mm); envelope 2.0 on a 41 × 41 foil; the tilted-plate probe (cos i
     0.6–0.8 at 4–5 mm, thor `visrad/build_profile_probe.jl`); a ~70° foil; shadowing (a
     washer or sphere before the foil, a hohlraum with LEH clipping);
   - end to end: a primitive hohlraum optimized here, exported, run and saved computed in
     VisRad, read back with `visrad_deposition` (its element order measured first, or
     the hohlraum given as CAD).
2. **Backlighter geometry (optional).** Shield and foil are done. Left: the shield's
   mounting sector (the sheet runs past 9.2 mm there); backlighter beams aimed off the
   foil center or foils off the equator (could bring the 13.2 mm standoff in). *Data:* a
   shot's real backlighter dimensions.
3. **Diagnostic windows.** The design studies are in the tutorial. Left: the expanded
   foil's density (held at Lester's 1 g/cm³), a backing that absorbs and re-emits part
   of what the gold transmits, a computed gap-closure time (plasma filling the gap), and
   the windows' l = 4, which neither quad powers nor pointing removed.
4. **Confirm against facility practice.** The TaLIS 0.5° tolerance (pivot; replaces or
   adds to the error set); the 1ω footprint side and shape; which quads count as open
   for the 20° rule; the NIF CMT XML schema; the VisRad items in the research notes
   ("To verify with real files at work": summary-CSV columns after Power, spot-rotation
   sense, SST preamble, element-export layout; in workspaces the Rotation/Spin mapping,
   the cylindrical position order, primitive grid spacing, custom profiles' spot-size
   definition; one small workspace with a rotated cylinder and a rotated CAD part would
   settle most).
5. **Benchmark against VisRad.** Run `write_visrad_kit` cases at work and bring the
   outputs back for `compare_visrad`: absolute drive levels and mode amplitudes.
6. **VisRad import and CAD: what is skipped.** Each is warned about or refused today;
   most need a sample workspace using the feature.
   - Unconverted light (2ω/1ω beams on) as `nif_unconverted` colors or extra beams;
     time-dependent beam powers (power mode ≠ 0; VisRad's time unit is seconds),
     energy-normalized beams and named pulse shapes as `PulseShape`s; DPR beams (split
     into two cones); pointing method 1 (toward a diagnostic port) and target-mount or
     diagnostic-port frames; objects placed by a key point other than the centroid.
   - Time-dependent albedo, conversion, lengths and power sources (the ten
     `Time-Dependent …` tables), laser reflectivity, power self-emission; wall thickness,
     partial spheres and the other primitives (cylinder, cone, disk, box, polygon,
     halfraum, rugby and slotted hohlraums, surfaces of revolution, multi-element
     objects).
   - Remeshed CAD: receivers through an element's outline instead of its ~10–14
     triangles; an adaptive `angle` for thin rods and fillets; closing cracks where two
     faces tessellate one edge differently.
   - A transmitting thin-film option for CAD tents and LEH windows (or thin-solid
     detection in the library), so the cavity can include the tent region.
   - CAD fits: rims on the outer half of a torus (a fillet meeting a cone) stay second
     order (a facet-aware (φ, ψ) map would fix it); one surface shared across parts (the
     halves' inner walls); a quadrature at a cone apex; TaLIS curvature kinds for cones
     and tori.
   - Adjacent flat CAD pairs: near-kernel evaluations scale as shared-edge length over
     r/256 (~700 per pair on a sliver tube); an exact or coarser rule along straight
     shared edges would cut them; circular plane families would extend banding to the
     NIF halves; `_earclip` falls back to raw triangles on pinched outlines.
   - Laser deposition on sliver CAD is slow in `any_hit` (~5 s per beam on the full
     157k-element NIF CAD target; moot on the cavity).
7. **GUI on CAD.** Partial curved fits of CAD parts (1.5× faster view factors and 0.6 %
   on the halves' flux on the NIF workspace; wait for a workspace that benefits); the
   0.75 s first-import tick and 0.16 s first field-coloring frame (compilation; not
   remeasured since the 2026-10-25 compile work); `include_components!` 0.3 s on the UI
   thread; a "fit the hohlraum template to CAD parts" action; the foil objective with one
   shared parameter set for beams aimed at one point does not move. `fit_surfaces` cost:
   a variable-radius fillet is quadratic while its tolerance-stopped growths
   stay under the 512-patch covers (16–50k facets: 16–107 s); the 133k-facet NIF mount
   now 29 s (was 1100 s; research notes, "Cost of fitting large CAD parts").
8. **GUI requests.** Frames on a port axis; clipping volumes; STEP import (needs a CAD
   kernel; STL and VisRad workspaces cover CAD meanwhile); an x-ray image that follows
   the camera; adaptive refinement of the radiation mesh where deposition varies fastest
   (the uniform fine mesh, ~16k elements, now runs without factorizations, item 14).

Items 9–13 and 16–17 follow the 2026-10-22 literature review against Kritcher et al. 2024
(research notes, "Literature review against Kritcher et al. 2024"; `reference/litreview/`).

9. **The P2 swing (N221204, PRE 109 Fig. 10).** *State:* moving wall, closing lips
   (prescribed or computed) and the bubble are calibrated to published data; nothing is
   fitted to Fig. 10. On the current defaults (coarse; "Capsule absorption, conversion
   and the ablating-wall release" §2) the swing reaches −0.181 at 5.25 ns and turns
   positive at 6.19 ns with the measured backscatter, −0.195 and 6.24 ns with the
   uniform cut Fig. 10 is compared with (rms 0.044); Fig. 10: −0.170 at 5.6 ns, 6.67 ns.
   Every published input tried (K24 capsule radius, HYDRA's backscatter bookkeeping and
   post-CBET cone fraction, glint, Kemp 2026's SBS history, the late refraction as a
   bound) moves the crossing by −0.10 to +0.16 ns ("The P2 swing", "Inner-cone SBS time
   histories"); ~0.3 ns remains, which only the least constrained bubble inputs can carry
   (picket fluence `energy` 0.8: +0.33 ns; collision ±0.3 ns: +0.17 / −0.11). The
   minimum's timing moves with none of them. *Remains:* the HYBRID-E bubble speed and
   collision time; Izumi 2023's low-contrast stall and relaunch; the computed lips'
   thermal-band radius against the 3–5 keV edge they are fitted to; the refraction at the
   hydro's own timing (crashed on a libuv assertion; rerun `julia --project=. -t 4
   examples/n221204.jl 0.8 0.4 full backscatter=measured refraction=0.5`). Pending the
   user: `backscatter=measured` stays the default and `history` an option (the
   coordinator's choice). Validation left: Young's dynamic albedo lowers mode 1 ~20 %
   where he reports 4.6 % (needs N170827's delivered per-beam powers); Ralph 2018's
   18.5 % cone fraction for zero flux P2 (needs the high-foot pointing and wall radius);
   the N210808 proxy (its pulse and foot cone fraction are approximated); Callahan 2018's
   39-shot curve (per-shot inputs not tabulated; only its line is digitized). On the
   2026-10-24 baseline the rise sat 0.06 above Fig. 10 at 4.5 ns (DN24's input cone
   fraction there is below Fig. 9b's) and the foot ~0.03 below (mock pointing, foot cone
   fraction read by eye); not rechecked on the current defaults. *Data:* a
   HYBRID-E thin-wall or G-LEH measurement of the bubble tip or inner spots against time;
   a HYBRID-E LEH-radius history (tests the lips' hold depth); the N221204-design 2DConA
   in-flight shell P2; the HYDRA multipliers behind Fig. 10 and where its flux is
   evaluated; N221204's shell thickness at 5–7 ns; N221204's per-cone backscatter history.
10. **Shot-data calibration (CBET and effective parameters).** *State:* `n221204_cbet()`
    (N221204 Fig. 9, rms 0.016; the N221204 example's default), `hybride_cbet()` (joint
    fit with Te ≤ 5.5 keV and the null in 1–1.55 Å: rms 0.027 / 0.016 / 0.047 on N210808 /
    N220919 / N221204) and `hybride_cbet_iraum()` (null in 0.25–0.5 Å: rms 0.024 / 0.016 /
    0.028, null 0.24 ± 0.20 Å, dCF/dΔλ 0.137 ± 0.011 per Å; the calibration example's
    default) ("Joint multi-shot CBET calibration" and its "the null"). HYDRA's flat
    cross-shot response contradicts the measured P2 sensitivity, so no flattening
    mechanism was added. *Remains:* fit effective parameters (cone multipliers, CBET
    multiplier and clamp, bubble transparency, wall speed, albedo) to measured in-flight
    and hot-spot P2/P4 and keyhole shock symmetry, with uncertainty; reproduce the P2-vs-Δλ
    sequence N210808 → N220919 → N221204 → N230729 (1.8 / 2.5 / 2.75 / 2.9 Å at 1ω:
    hot-spot P2 −2.3 / −10.8 / +0.3 / +6 µm) and Ruof 2026's rise cone-fraction response
    (0.21 → 0.24 raises shell P2 9.1 → 27.5 µm; 0.21 → 0.16 leaves it); the fitted Te's
    dip to 2.5 keV at 5.5 ns and the peak clamp (1.4×10⁻⁴) below HYDRA's decks (test a
    drive-following Te(t) prior, B16's shape); plasma-state tables as inputs, outer-cone
    depletion, an inner-quad SBS proxy (N220919's limit, the 2.9 Å overshoot), the
    23.5°/30° split (the model favours 30°, HYDRA splits evenly). HYDRA's own knobs for
    comparison: clamp 0.002 (foot) / 0.0075 (peak), 2DConA fits 0.005–0.007, laser
    degraded 9–13 % or peak multipliers 0.825–0.86 (`reference/litreview/cbet_calibration.md`).
    *Data:* N230729 transfer curves (none published); per-shot plasma states.
11. **Implosion metrics and UQ.** Done: the metrics, Young 2020's dropped-quad and offset
    checks, Nora's high-yield branch and shape term ("Implosion metrics", "Young 2020
    validation"). Left: time-dependent error budgets and multi-shot trending.
12. **Backlighter physics.** Done: He-α CE(hν, I) with its uncertainty, filters, the
    photometric chain. Left: pulse-length and spot-size dependence (confounded with the
    picket in the published data), a weighted, larger low-Z compilation (Ti, V),
    continuum and higher lines in the filter pass band, detector response, the object's
    transmission.
13. **Spectral radiosity.** Done: multigroup radiosity, the laser M-band, a synthetic
    Dante with hardware channels, the unfold, Sentinel views, thermal laser spots. Left:
    measured Ni XRD and thin-filter transmissions, Dante 2's and Omega's channel
    hardware; a measured soft spectrum of spots in a cavity (the local-Tr spot temperature
    is an assumption); the per-cone M-band inside `drive_history` (the N221204 example
    weights cones by their pre-CBET deposition); the M-band conversion's angle dependence
    (sin⁴θ, Dewald 2020) and its emission from the outer bubble; emission of the closing
    LEH plasma (Dante's initial- and open-area readings differ by 20 eV at N221204's
    peak); a time-resolved M-band P2 against Dewald's Fig. 10 on a moving hohlraum; the
    M-band split in the `CapsuleDrive` handover (item 17).
14. **Scale.** Done: Float32 and iterative solves, `HMatrix` (49.8k elements in 57 s and
    2.2 GB), packet visibility in the deposition operator, sliver cuts ("The H-matrix's
    N^1.5 remainder"), zero near-field blocks dropped and facing-cut rows refined to
    low rank (−17 % storage, −15 % F x; "H-matrix near field"), a vectorized sector
    contour (−24 to −28 % assembly, round-off; a contiguous quadrature layout gained
    nothing and was dropped). Left, in order: the
    near-field quadrature itself; rim certificates for non-convex occluders and `:sampled` visibility
    (still a ray per pair of an occupied leaf); the bored cavity's gap pairs (~60 % of its
    view-factor time; refinement toward the shared edge onto ~34-edge sector receivers);
    GPU mat-vecs (≲ 25 % end to end on the RTX 2060 while assembly stays on the CPU; not
    worth CUDA.jl yet).
15. **Polar direct drive mode**, benchmarked against published GA optimizers. *State:* literature review and staged plan in `reference/litreview/polar_direct_drive.md` (2026-10-27): Barlow 2024 is the only published GA PDD optimizer; analytic cases (Murakami & Nishi 2017, Temporal 2014) are the reproducible stage-0 checks; open decisions there.
16. **Absolute energetics on N221204** (`reference/litreview/energetics.md`,
    `n221204_benchmark.md`, `late_time_drive.md`). *State* (current defaults, coarse;
    "Capsule absorption, conversion and the ablating-wall release" §2–3): capsule share
    13.8 % and LEH 24.7 % (HYDRA 16.2 / 25.2); capsule Tr 306 eV at 8 ns; Dante 1 at
    8 ns 290.9 / 314.2 eV over the initial / open LEH area (measured 310.1); Dante at bang
    −40.5 / −46.5 eV low on N191118 / N210601. The pointing stays the public mock-up
    (`hybride_mock_pointing()`); CPPs, window geometry, corners and the liner thickness are
    unpublished and stay at repository defaults. *Remains, each with its data:*
    - **The drive after the laser falls.** The model falls far faster than Dante on all
      three shots; neither a static nor an ablating gold wall holds the drive in an
      independent 0-D hohlraum, so the missing ~250–450 kJ is in the plasma (Kline 2011's
      corona, Town 2011), with an unpublished return rate; nothing is fitted. A 0.15 ns
      source lag (`lag=0.15`) is a sensitivity, not a model. *Data:* HYDRA's energy
      inventory through the coast, a post-shot synthetic Dante for N191118 or N210601,
      Dante spectra between laser-off and bang, Moody 2017's cooling constant, gated
      LEH/bubble images after laser-off; for a fast reservoir, a conversion response time
      or Dante's rise delay on a step pulse.
    - **The capsule share** (2.4–4.2 points below HYDRA's 16.2 % between the radius
      bounds). *Decided 2026-10-26:* `outside` stays 0 until a blow-off opacity or
      absorbing radius is published; or Young's N170827 albedo is too high for N221204's
      thicker W-doped ablator. Schedules are `RadiusSchedule`s (`TabulatedRadius`,
      `ScaledRadius` for the pair's accepted, derived K24-shaped one, `TrajectoryRadius`
      from published anchors; no rocket model; "Capsule radius schedules", bit for bit).
      K24's radius is the ablator–fuel interface or fuel centre of mass, so the front at
      6 ns lies in 0.84–0.906 R₀. *Data:* the blow-off opacity or HYDRA's absorbing
      radius; the shell thickness at 5–7 ns; the pair's own radius histories or v_imp/RPV.
    - **The early peak**: over the open area the model sits 13–20 eV above Dante at
      5.0–5.5 ns (standard mesh, "Energetics defaults and the capsule share"; its rise
      ~0.15 ns ahead of the measured one).
    - **Other sinks**: the He gas (4.5 % in HYDRA), LEH windows, the M-band in the drive
      history (item 13).
    - **The Zylstra pair**: the model's LEH difference is 2–4 eV (~3–4 % power) too
      large; N191118's pointing and cone fraction are unpublished.
    - **Hot-spot velocity**: 0 to round-off on the mock pointing (no mode 1); the
      measured 36 km/s needs per-quad delivered powers, window or capsule mode-1 inputs,
      none public.
17. **Briareus → CompactLES** (`reference/litreview/hydro_coupling.md`; research notes,
    "Capsule drive hand-over and CompactLES stage 0"). LLNL's two-step method, Briareus
    as the integrated step and CompactLES as the capsule-only step, coupled in memory
    through a package extension (`CompactLESBriareusExt`, in CompactLES on branch
    `briareus-drive`, unmerged; HDF5 output optional). The handover: T_r(t),
    capsule-averaged flux, a_l/a_0 for l ≤ 8 (P2 > 0 pole-hot) and every m on a ≤ 50 ps
    grid, the evaluation radius, provenance; later a radius table and the Planck + M-band
    split (item 13). *State:* stage 0 (an absorbed-flux ablation source, P_abl ∝ T^3.3
    t^−0.09) runs in 1-D, r–z, 3-D and on refined levels. The full-model N221204 drive
    (306 eV peak) at ṁ = 2.85e5 T³ gives 340 km/s, 5.0 % remaining, 24 kJ, stagnation
    10.2 ns (bang 9.53; targets v_imp 380–400 km/s, 5–6.6 % remaining, 19–33 kJ). A is an
    effective rate fitted to the remaining mass; Olson's planar HDC rates (4.15–5e5)
    over-ablate and the grid study fails off nx 1400, so a measured rate is CompactLES
    ROADMAP A22. *Next:* the stagnation lag (rocket model or drive); then r–θ with the
    a_l (Masse's P2 sign flip near 4.5 ns), r–z with AMR; a 1 µm 3-D shell needs the band
    (set by CompactLES's covered-node timestep and tile fragility) and the band following
    the shell over a regrid. Stage 1: implicit gray flux-limited diffusion with a Marshak
    boundary carrying the angle-resolved flux, conduction at 0.15, QEOS/FEOS and TOPS
    opacities, validated on Shestakov 1998's Be/D capsule at 160 eV and then Kritcher 2024
    Table III; stage 2 multigroup with the M-band. Bang time is an input in LLNL's
    capsule-only runs and the first prediction here. Pitfalls: apply the asymmetry near
    the ablation front, not at a distant fixed radius; the M-band has its own larger P2.
