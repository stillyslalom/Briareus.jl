# Roadmap

Open work only. Each item says where it stands, what remains and what data would settle
it. What has been done is in the git history and in `reference/research-notes.md`, under
the sections named here. Items marked *(survey)* come from the code survey of 2026-11-02;
those under "Correctness" were read in the code, the rest are leads to confirm first.

## Release

1. **Releases 0.2.** Briareus 0.1.0 and BriareusGUI 0.1.0 are in General, both at the
   commit of tag `v0.1.0` (6d032b1; General PR 174205 for the GUI). Main has moved on:
   the GUI on main calls core API added since 0.1.0 (`FaceFlux`, `wall_temperature`,
   `HohlraumProfile`, `port_frame(facility, port)`, `drive_history`'s `per_component`,
   among others) while its `Project.toml` still allows `Briareus = "0.1"`. The next
   release registers core 0.2 (APIs broke since 0.1) and the GUI 0.2 with its compat
   raised to `Briareus = "0.2"`, core first. Check `pkg> app add BriareusGUI` and
   `briareus --shortcut` from the registry once the 0.1 GUI is installable.
2. **Desktop launch.**
   - The macOS app bundle and splash are untested (the macOS `qml` tool needs each Qt
     framework folder on dyld's fallback path).
   - Whether the Windows shortcut groups with its window on the taskbar (the window sets
     its own AppUserModelID; the `.lnk` does not carry it).
   - The first frame takes 21 s on Windows: 15 s loading package images and 3 s compiling
     in QML's `__init__` (upstream; a system image would remove both). A bare `briareus`
     in a terminal shows nothing until the window.
   - The precompile workload changes the default drive at round-off (P0 5e-17 against
     7e-17 built without it); harmless, unexplained.

## Correctness

3. **Bugs and weak spots found by the survey** *(survey; GUI items read in the code)*.
   - GUI: an undo of a part rename brings the old name back, but the view state keyed to
     the new one (display modes, the cut, triads) is not in the history and is dropped; a
     diagnostic definition's `center_part` is the user's file and keeps the old name.
   - Core: the `line_of_sight_power` occlusion is per centroid (all or nothing), unlike
     Dante's pixel rays; `ApertureKeepOut` takes every hole as a circle; primitive targets
     default to `visibility = :centroid` (row sums 1–6 % off at LEH rims unless
     `:sampled`).
4. **Engineering closures to validate.** The re-heated wall's radiosity clamp (≥ 0,
   albedo 0, no emission; `history.jl`, eae01d4); the conversion of overlapping colors at
   each color's own flux; the liner blend in `WallLossAlbedo`. *Data:* a radiation-
   hydrodynamics run of a wall heated, cooled and heated again (a prepulse and gap).
5. **Two cavities joined by a package.** View factors are balanced for a single closed
   hohlraum only (`_closed_enclosure`); halfraum pairs, the shock tube and other multi-
   cavity targets run unbalanced. Needs per-cavity enclosures.

## VisRad parity and facility practice

6. **VisRad deposition: shadowing and beam totals.** *State:* `VisRadProfile` is VisRad's
   own profile law: 0.1–0.3 % rms per element on plates, 0.45–0.7 % on hohlraum walls
   ("VisRad probe kits (2026-10-27): results"). *Remains:*
   - VisRad rescales each beam's total after shadowing: to its power when every ray hits
     something (per-plane normalization is +0.2–2.2 % high on steep foils), and to a share
     within ±2.5 % of the geometric one when part of the beam hits nothing (grid not
     identified).
   - VisRad's shadow rays lean outward in proportion to spot position (a virtual source
     12–19 mm upstream); axis shadowing misses by ~30 % rms behind occluders. An
     aim-dependent mask conflicts with an operator built once per geometry.
   - The LEH-washer leak (≤ 7 × 10⁻⁴ TW/cm²).

   *Data:* occluders ~2 and 10 mm upstream and a foil off focus; the Laser Beam Summary
   export of the spot-rotation workspace (`read_visrad_summary` end to end).
7. **Benchmark against VisRad.** Run `write_visrad_kit` cases at work and compare with
   `compare_visrad`: VisRad's radiosity, absolute drive levels and mode amplitudes
   (deposition alone is compared so far).
8. **VisRad import and export: what is skipped.** Each case is warned about or refused;
   most need a sample workspace that uses the feature.
   - Beams: unconverted light (read VisRad's green/red powers into the ratios of
     `unconverted_deposition`; the import still takes 3ω only); time-dependent powers
     (power mode ≠ 0; VisRad's time unit is seconds), energy-normalized beams, named
     pulse shapes; DPR beams; pointing method 1 and target-mount or diagnostic-port
     frames; objects placed by a key point other than the centroid; quad beams pointed
     differently are averaged to one aim (warning only).
   - Objects: the `Time-Dependent …` tables, laser reflectivity, power self-emission;
     wall thickness, partial spheres and the other primitives (cylinder, cone, disk, box,
     polygon, halfraum, rugby and slotted hohlraums, surfaces of revolution,
     multi-element objects); the element order of primitive spheres and hohlraums;
     `remesh`/`coarsen` do not reach primitive objects (they stay on VisRad's grid).
   - Workspace details to check on one small workspace with a rotated cylinder and a
     rotated CAD part: the Rotation/Spin mapping (tested on a 90° tilt only), the
     cylindrical position order, primitive grid spacing, custom profiles' spot-size
     definition; the export's column spellings, table order, fit weighting and odd-M sign.
   - Export: spot halos are dropped (VisRad has none); the GUI exports only an imported
     workspace.
   - Frames phase 4 (`reference/frames-design.md`): `Ref. Coord. Sys.` export, object
     hierarchies as anchored parts, the port frame's x̂ (confidence L).
9. **Confirm against facility practice.** The TaLIS 0.5° pivot per component (TCC or the
   mount); the 1ω footprint's side and shape; which quads count as open for the 20°
   reflection rule; the NIF CMT XML schema; the SST preamble and the element-export layout
   ("To verify with real files at work"); the NIF cone pointing of the GUI presets and the
   7700 mm lens distance (unsourced in `beam_table.jl`, `beam_geometry.jl`).

## Physics against shot data

Items 10–15 follow the 2026-10-22 literature review against Kritcher et al. 2024
("Literature review against Kritcher et al. 2024"; `reference/litreview/`).

10. **The P2 swing (N221204, PRE 109 Fig. 10).** *State:* nothing fitted to Fig. 10; the
    swing reaches −0.181 at 5.25 ns and turns positive at 6.19 ns (Fig. 10: −0.170 at
    5.6 ns, 6.67 ns). Every published input tried moves the crossing by −0.10 to +0.16 ns;
    ~0.3 ns remains, which only the least constrained bubble inputs can carry ("The P2
    swing", "Inner-cone SBS time histories"). *Remains:* the HYBRID-E bubble speed and
    collision time; Izumi 2023's low-contrast stall and relaunch; the computed lips'
    thermal-band radius against the 3–5 keV edge they are fitted to; the rise and foot on
    the current defaults (0.06 above and ~0.03 below Fig. 10 on the 2026-10-24 baseline).
    Validation left: Young's dynamic albedo (needs N170827's delivered per-beam powers);
    Ralph 2018's 18.5 % cone fraction (needs the high-foot pointing and wall radius); the
    N210808 proxy; Callahan 2018's 39-shot curve. *Data:* bubble-tip or inner-spot
    histories (thin-wall or G-LEH), a HYBRID-E LEH-radius history, the 2DConA in-flight
    shell P2, the HYDRA multipliers behind Fig. 10, N221204's shell thickness at 5–7 ns and
    per-cone backscatter history.
11. **Shot-data calibration (CBET and effective parameters).** *State:* `n221204_cbet()`,
    `hybride_cbet()` and `hybride_cbet_iraum()` ("Joint multi-shot CBET calibration").
    *Remains:* fit effective parameters (cone multipliers, CBET multiplier and clamp, bubble
    transparency, wall speed, albedo) to measured in-flight and hot-spot P2/P4 and keyhole
    symmetry, with uncertainty; reproduce the P2-vs-Δλ sequence N210808 → N230729 and Ruof
    2026's rise cone-fraction response (the calibrated models' dCF/dΔλ is 2–3× steeper
    than HYDRA's flat cross-shot response); the fitted Te's dip at 5.5 ns and the peak
    clamp below HYDRA's, Te at the 5.5 keV fit bound at 8 ns (try a drive-following Te(t)
    prior); plasma-state tables as inputs, outer-cone depletion, an inner-quad SBS proxy,
    the 23.5°/30° split; the two HYBRID-E models share most code and a hard-coded
    `leh_z`. *Data:* N230729 transfer curves; per-shot plasma states.
12. **Absolute energetics on N221204** (`reference/litreview/energetics.md`,
    `late_time_drive.md`). *State:* capsule share 13.8 % and LEH 24.7 % (HYDRA 16.2 /
    25.2); Dante 1 at 8 ns 290.9 / 314.2 eV over the initial / open area (measured 310.1),
    computed through the hole before Dante 1 became the whole-field view (7e-5 apart on
    the GUI default; not rerun on N221204). *Remains:*
    - The drive after the laser falls: 37–54 eV low at bang on three shots; the missing
      ~250–450 kJ is in the plasma, with an unpublished return rate. *Data:* HYDRA's energy
      inventory through the coast, a post-shot synthetic Dante, Dante spectra between
      laser-off and bang, gated LEH images after laser-off.
    - The capsule share (2.4–4.2 points below HYDRA). *Data:* the blow-off opacity or
      HYDRA's absorbing radius; the shell thickness at 5–7 ns; the radius histories.
    - The early peak: 13–20 eV above Dante at 5.0–5.5 ns over the open area.
    - Other sinks: the He gas (4.5 % in HYDRA), LEH windows; laser light reflected from a
      surface is a loss, not re-deposited (no specular or diffuse laser reflection).
    - The Zylstra pair's LEH difference is 2–4 eV too large (N191118's pointing and cone
      fraction are unpublished); the hot-spot velocity needs per-quad delivered powers.
    - The research notes' "N221204 late-time drive ▸ What remains" predates the later
      sections (its capsule share 11.6/12.3 % is superseded by 13.8 %): mark it.
13. **Spectral radiosity and diagnostics** (`reference/diagnostics-design.md`).
    - Dante: measured Ni XRD and thin-filter transmissions; the SiO₂ mirror's angle scale,
      roughness and surface layer (assumed or fitted; extrapolated below 2.2 keV); the
      within-group spectral shape (finer groups); Dante 2's hardware and fields and
      Omega's (Dante 2 and Omega's Dante look through the hole: fields not public;
      `dante2()` is not definition-backed); port names from the facility table (moves
      Dante 1 by 0.25°); a field of view per channel (the multilayer channel's 12 mm); time
      response and per-channel uncertainty, and an uncertainty band from `dante_unfold`;
      recording modes other than `history`; field-of-view checks against Chen 2024, the
      BAH 27 % ring occlusion and May 2015's cosine ratio. Sessions saved before 46d8a76
      keep their through-hole Dante 1.
    - Diagnostic definitions phase 5: the imagers' SXI soft channel, per-pixel gate times
      and near-pinhole perspective (and whether SXI and G-LEH 3–5 keV images need a
      non-thermal spot model); backlit radiography in the GUI (definitions are listed and
      can be ticked but record nothing), with ARC, opacity and streaked radiographs and a
      continuum, satellites, photon noise and scatter in the radiograph; streaked VISAR
      and SOP; only three diagnostic kinds are registered.
    - Radiosity: a measured soft spectrum of spots in a cavity; the M-band's angle
      dependence (Dewald 2020) and its emission from the outer bubble; the closing LEH
      plasma's emission; a time-resolved M-band P2 against Dewald's Fig. 10; the M-band
      split in the `CapsuleDrive` handover; `_ALPHA_M` rests on one observable (Dewald 2020
      F_M); the spectral solve is Float64 only (no gradients for M-band or Dante
      objectives); block-GMRES and thermalized solves error out on non-convergence.
14. **Implosion metrics and UQ.** Time-dependent error budgets and multi-shot trending;
    focus, spot-size and timing errors beside pointing, power and positioning (linear
    propagation, Monte Carlo as the check); groups over both hemispheres are refused;
    `burn_amplification` refuses Y_amp 40–80 (a gap in the published fits);
    `flat_kernel` stands in for Young's mode-1 Green's function; check the Yₗₘ → Young
    Cartesian sign convention before comparing vectors
    (`litreview/asymmetry_transfer.md` §7).
15. **Backlighter physics and retune.** The low-Z end of He-α CE (one Ti point; no V/Sc
    data at 351 nm); pulse length and spot size as variables (CE ignores both); continuum
    and higher lines in the filter pass band, detector response, the object's
    transmission. Retune (core `retune_drive` with `FoilIrradiance`/`FaceFlux`, level
    held; "Retuning faces and foils; the held level"): the GUI offers the reference target
    only (no uniform target at a set level), lets the user choose neither the quantity
    (`on = :auto`) nor load the retuned pulses as the Time panel's; its foil retune leaves
    out unconverted light and the prepulse; peak-to-mean is reported, not optimized; the
    retune holds the reference's albedos (linearized) and retunes powers, not pointing.
16. **Unconverted light over time.** *State:* `unconverted_history` and
    `irradiance_history`, `specular_check` and `preplasma_check` against the User Guide's
    10⁸ W/cm², the GUI's Time tab *Irradiance*, the prepulse in the drive history's Tr
    ("Unconverted light over time", "The prepulse in the drive history"). *Remains:*
    - pointing or geometry that changes in time (a deposition per configuration);
    - a fluence-based damage criterion beside the irradiance threshold, and where the
      specular 1ω reflection goes (the check tests the threshold only);
    - the footprint's phase-plate edge blur (a flat-top square now).

    *Data:* the unconverted fraction vs. pulse shape and the prepulse level; plasma,
    damage and specularity thresholds of Au, Al and Si at 1ω and 2ω.
17. **Per-color sources and a materials library.** *State:* per-color sources in `solve`,
    `drive_history`, `flux_response` and `retune_drive`; the GUI's drive, histories and
    capsule retune take the Drive panel's ratios and prepulse ("Unconverted light:
    absorption and conversion vs. wavelength"). *Remains:*
    - the M-band of unconverted light (the M-band stays the 3ω light's);
    - data on the conversion of overlapping colors (each color at its own flux is an
      engineering default) and on whether the 3ω intensity table applies to 1ω and 2ω
      (Mead's Iλ² scaling would shift it);
    - a library of material defaults (albedo, conversion vs. intensity and wavelength,
      absorption, cold reflectance), each with its data source and uncertainty, in
      place of the per-part constants (the GUI's hohlraum 0.85/0.8, capsule 0.3/0.3, foil
      conversion 0.02, free parts 0/0 carry no source).

## Coupling and new modes

18. **Briareus → CompactLES** (`reference/litreview/hydro_coupling.md`; "Capsule drive
    hand-over and CompactLES stage 0"). *State:* `CompactLESBriareusExt` is on CompactLES
    main; stage 0 runs in 1-D, r–z and 3-D, with the r–θ P2 history and r–z/3-D on refined
    levels; the N221204 drive gives 340 km/s and stagnation 10.2 ns (bang 9.53); the late
    lag is the rocket closure (CompactLES A22, branch `briareus-drive-2`, two commits
    unmerged). *Next:* a 1 µm 3-D shell (blocked by the tiled levels and the timestep);
    stage 1, gray flux-limited diffusion with a Marshak boundary, validated on Shestakov
    1998 and then Kritcher 2024 Table III; stage 2, multigroup with the M-band.
19. **Polar direct drive** (`reference/litreview/polar_direct_drive.md`). *State:* stage 0
    (hard-sphere illumination, `mode_rms`, spherical pointing, the LMJ fixture), the corona
    operator (stage 1) with its absorbed power per steradian and radius (`shells`) and the
    Ba24 P_abl and Ye21 proxies, the Ba24 σ₁ comparison and optimization. *Remains:* the
    CompactLES snapshot example; the Ye21 m ≠ 0 benchmark (needs a comparable profile);
    OMEGA fixtures (`omega_ports` is an idealized icosahedron, not chamber coordinates or
    port IDs; OMEGA EP, LMJ beyond the fixture and Z are absent); stages 2–3 (a corona
    time series, direct-drive CBET; the corona is static, inverse bremsstrahlung only);
    open decisions 2, 10 and 13–15 in the review. The corona has no `Scenario` field (an
    objective captures the `CoronaOperator`, at the cost of the scenario's unused surface
    deposition per evaluation), and the proxies are not wired into `evaluate`.
20. **Photometric forward modeling of x-ray sources.** From a source's spectrum and
    history (a backlighter's He-α chain, the hohlraum's re-emission, `spectral.jl`)
    through cold opacities of filters, shields and the object to a detector's response,
    to aid design and the choice of diagnostics:
    - the best scintillator, and the best gate time for time-integrated detectors;
    - the expected signal on film and image plates;
    - the shielding a radiograph needs next to an igniting hohlraum.

    It builds on `backlighter_signal` and "Backlighter conversion efficiency and the
    photometric chain". *Data:* public cold opacities (Henke, NIST XCOM), detector
    responses (image plate, film, scintillators), the x-ray, neutron and γ backgrounds of
    igniting shots.
21. **Beam and pointing models.** Beams carry an elliptical super-Gaussian spot (and an
    optional halo) only: no tabulated far fields, polarization or SSD/bandwidth fields
    (SSD is a clearance smear); at most 16 f-cone directions (`UInt16` mask);
    `VisRadProfile` is measured at f/20.755 only; the off-focus `intensity` widens moments
    at a fixed order (used by clearance and keep-outs). The optimizer is local L-BFGS-B
    with penalties (no multistart, no robust optimization under pointing error);
    `placement_sweep` is 1-D translation only.

## Targets, CAD and the GUI

22. **Assemblies: phase E and leftovers** (`reference/assemblies-design.md`). Phases A–D
    and C2, the placement triad, piece roles and the package seams are built. *Remains:*
    - Phase E, the double backlighter with ARC: the starter, ARC beams as a facility beam
      set (public description only; mock-up timing), per-foil results if one active foil
      is not enough.
    - Profiles: windows on cones and across breaks, fillets on a profile, a frustraum
      pointing preset; the tentless membranes.
    - Roles: several observed surfaces (decision 3, "later"); roles on VisRad and CAD
      pieces (only hohlraum-part pieces hold their own); an observed washer, cylinder,
      cone or flat CAD face; per-cavity albedo and energy series; a Dante field centered
      on a part from the GUI; VisRad primitive pieces have no frames.
    - Triad: a hover highlight, planar (two-axis) drags, several selected parts.
    - Tree: drag-and-drop into assemblies, duplicate, a bare hohlraum part in Add part,
      expanded nodes saved; Remove deletes the primary selection only.
    - The design docs keep built phases against their own rule (remove each phase once it
      lands; `assemblies-design.md` "Today", `frames-design.md` phases 1–3 and phase 5,
      done by the foil starter).
23. **GUI: workflow gaps** *(survey)*.
    - No optimizer objective for a flat observed face (halfraum and opacity designs
      cannot be optimized for face uniformity); modes stop at l = 4 in the tables,
      weights and retune.
    - Beams: spot shape, f-number, wavelength and port are saved but not editable; aims
      are Cartesian only; VisRad profiles and focal lengths only from an import; no OMEGA
      EP preset; applying SST or a beam summary needs a preset first; no aim-frame column
      or sorting in the beam table.
    - Time: no pulse editing; the loaded pulses are not saved with the session; no Time
      node in the Analyses tree (its prepulse and threshold settings sit in Drive).
    - Sessions and undo: view edits (display modes, cutaway, hidden beams, beam style)
      and loading pulses are not undoable; camera pose and chart choices are not saved;
      ticking a diagnostic makes design undo steps.
    - Views: Look from port lists beam ports only; *Save view as image* grabs the main
      canvas only; a beam within the pick distance wins over a nearer surface.
24. **GUI: Help, tooltips and layout** *(survey; the rule: explanations in tooltips, ⓘ and
    Help, never page paragraphs)*.
    - Controls outside the survey's list still without tooltips: Target, Component,
      Frames and placement-step properties, the model tree's dialogs, Look from port.
25. **GUI density and layout (user review, 2026-10-29).** Aim toward Office-ribbon
    density, not all the way; the 3D view first. *Left:* the next pass on what still takes
    room once the user has used it. Principle: frequent commands near at hand,
    explanations in tooltips and Help, the 3D view and its inspection given the most room.
26. **GUI responsiveness** *(survey)*. `can_undo` serializes and deep-compares the whole
    design on every refresh (4×/s while jobs run) and `record!` serializes again per
    edit; the model tree and every panel are rebuilt each refresh; undo rebuilds the
    target (CAD designs depend on the 6-entry target cache); targets under 30k triangles
    build on the GUI thread; each sweep step rebuilds the target.
27. **Backlighter geometry (optional).** The shield's mounting sector (the sheet runs past
    9.2 mm there); backlighter beams aimed off the foil center, or foils off the equator
    (both expressible with `backlighter_foil(los)` and `FramePointing`; neither explored).
    *Data:* a shot's real backlighter dimensions.
28. **Diagnostic windows.** The expanded foil's density (held at 1 g/cm³), a backing that
    absorbs and re-emits part of what the gold transmits, a computed gap-closure time, and
    the windows' l = 4, which neither quad powers nor pointing removed.
29. **CAD.**
    - A transmitting thin-film option for CAD tents and LEH windows, so the cavity can
      include the tent region.
    - Remeshed CAD: receivers through an element's outline; an adaptive `angle` for thin
      rods and fillets; closing cracks where two faces tessellate one edge differently.
    - Fits: rims on the outer half of a torus stay second order; one surface shared
      across parts; a quadrature at a cone apex; TaLIS curvature kinds for cones and tori;
      free-form faces stay facets; the cover and rms constants are tuned on NIF parts
      only.
    - View factors: mid-tier flat pairs by the double contour; circular plane families for
      the NIF halves; `_earclip`'s fallback on pinched outlines ("Near flat pairs by the
      double contour").
    - Cost: laser deposition on sliver CAD in `any_hit` (~5 s per beam on the 157k-element
      NIF CAD target); `fit_surfaces` on a variable-radius fillet is quadratic (16–50k
      facets: 16–107 s) and takes 1278 s on the 133k-facet mount.
    - Input: STL only (no STEP, OBJ or PLY; STEP needs a CAD kernel).
30. **GUI on CAD.** Partial curved fits of CAD parts (wait for a workspace that benefits);
    the first import's commit tick (1.5 s) and first field coloring and selection (0.2 s
    each); the CAD capsule's remesh (1.0–2.0 s, now clustering); a "fit the hohlraum
    template to CAD parts" action, which would give a CAD hohlraum piece roles and a
    per-cavity Tr. Confirm in the window that the foil objective moves with one shared
    `FramePointing` (frames phase 1 should have fixed the stuck case), then drop it.
31. **GUI requests.** Clipping volumes; an x-ray image that follows the camera; adaptive
    refinement of the radiation mesh where deposition varies fastest.
32. **Scale.** Segments between the capsule's balls (~21M rays left in the 11.4k bored
    cavity's `:sampled` build, ≲ 7 % of it); GPU mat-vecs (≲ 25 % end to end; not worth
    CUDA.jl yet) ("H-matrix visibility certificates", "Prism certificates for the waist
    gap"); a cheaper predictor for ACA give-ups (~20 % of assembly) and near-field
    quadrature reuse; an `HMatrix` has no direct solve; `storage = Float64` densifies an
    `HMatrix` without warning; `Keyframes` build a deposition operator at every time.
33. **Code health** *(survey)*. Thin or indirect tests: `lpi.jl`, `keepout.jl`,
    `xray_optics.jl`, `bvh.jl`, `storage.jl`, `facilities.jl` have no test file of their
    own; a dozen exported functions have no test reference (`beams_from_table`,
    `beam_component_power`, `raytrace_component_power`, `with_axes`, `talis_requirement`,
    `ylm_frame`, the Young 2020 helpers); GUI radiography, `ungroup_assembly!`,
    `add_hohlraum_break!` and `rename_cluster` are tested once or not at all. Shared
    helpers: the several linear interpolators, `_trapz` (in `implosion.jl`, used by
    `imager.jl`), `_data_table` (belongs in `io.jl`), the triplicated coarsen-xor-remesh
    check; benchmark toggles (`_ITERATIVE_SOLVE`, `_SOLVE_STATS`) in production code; a
    `dante_unfold` docstring names a `solver` keyword that does not exist.
