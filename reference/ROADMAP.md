# Roadmap

Open work only. Each item states where it stands, what remains and what data would settle
it. What has been done is in the git history and in the research notes
(`reference/notes/`; `reference/research-notes.md` lists every section by file), under
the sections named here. Items marked *(survey)* come from the code survey of 2026-11-02;
those under "Correctness" were read in the code, the rest are leads to confirm first.

Items are grouped by theme, and an item keeps its number for as long as it is open: code
comments, data headers and the literature reviews cite items by number. A closed item
leaves a gap (1, 3 and 5 are closed), and a new item takes the next free number wherever
its theme puts it (34–36 under "Re-pointing studies", 37 under "Correctness").

**Priorities** (user review, 2026-11-18). First, the hohlraum physics against shot data:
the gold bubble (item 10), CBET (item 11), then the energetics and the wall and spectral
closures (items 12, 4, 13, 17). Second, re-pointing studies: the public API they need
(item 34), capture and clearance rules inside the optimization (item 35) and the
constrained runs that stall at their bounds and force staged tolerances (item 36).

## Release

2. **Desktop launch.**
   - The macOS app bundle and splash are untested (the macOS `qml` tool needs each Qt
     framework folder on dyld's fallback path).
   - It is untested whether the Windows shortcut groups with its window on the taskbar.
     The window sets the AppUserModelID itself, and the `.lnk` does not carry it.
   - The first frame takes 21 s on Windows: 15 s loading package images and 3 s compiling
     in QML's `__init__` (upstream; a system image would remove both). A bare `briareus`
     in a terminal shows nothing until the window.
   - The precompile workload changes the default drive at the round-off level (P0 is
     5e-17 with it against 7e-17 built without it). The change is harmless and
     unexplained.

## Correctness

4. **Engineering closures to validate.** Three closures need validation: the radiosity
   clamp of the re-heated wall (≥ 0, albedo 0, no emission; `history.jl`, eae01d4), the
   conversion of overlapping colors at the flux of each color, and the liner blend in
   `WallLossAlbedo`. *Data:* a radiation-hydrodynamics run of a wall heated, cooled and
   heated again (a prepulse and gap).
   The wall albedo of a fixed-time scenario also needs a check (scratch study,
   2026-11-18). A standalone `solve` evaluates `WallLossAlbedo` at the time since the
   pulse start (0.754 mean over the wall at 1 ns in an ignition-scale picket), while
   `drive_history` follows the temperature history of each element (0.607 at the same
   time). A phase scenario uses the per-element albedos of the history through
   `history_state(target, history, time)`, with the wall mean as an option. `Scenario`
   and `remake` take its `time` and `albedo`, so one target is used for every phase, with
   the windows intact or burned through. On a public mock-up the drive of the standalone
   solve is 0.85–1.44 of the history's at four times of a pulse, and its P2/P0 is off by
   1.1–2.8 points; the wall mean of the history, used as one constant, is within 0.7%
   and 0.5 points ("Drive histories: sliver steps and
   the albedo of a phase"). A `Scenario` with a response on a `WallLossAlbedo` target
   (or a `BurnthroughAlbedo` over one) without `albedo` takes the 250 eV reference
   temperature (drive up to 1.64 of the history's mid-pulse) and now warns at
   construction, naming `history_state`; `albedo = :reference` marks the reference
   temperature as intended ("Inner-run options, a wall-loss scenario warning and laser
   delivery figures").
37. **`drive_history` radiosity convergence** (scratch study, 2026-11-18). A quasi-static
    history of a hohlraum with two thin-foil windows (`BurnthroughAlbedo`, a
    `WallLossAlbedo` wall, `wall_history` on, `transition = 0.05`) stopped at 0.1 ns steps
    with "iterative radiosity solve (radiosity) did not converge in 400 iterations" for
    190 beams; the same run at 0.2 ns steps, and at 0.1 ns with all 192 beams, converged.
    The cause was that the `transition` planner ended a step 8e-6 ns short of a goal
    time and took the remainder as a separate step. On a cooling wall the implicit part
    of the release scales as Δ^(b−1), which put the albedo of the wall elements near
    −27, and GMRES does not solve that system at any iteration cap. The planner now runs a
    step to the goal when less than `min_step` would be left; a public mock-up test
    checks it, and the run of the study completes. *Remains:* steps much shorter than
    `min_step` given by the caller on a cooling wall still stiffen the release
    (documented); a bound on the short-step factor in the wall-release kernel would be
    a physics change (user).

## VisRad parity and facility practice

*Data on hand (2026-11-18):* many VisRad runs are saved outside the repository: hohlraum
workspaces with VisRad's computed results at several mesh sizes, the probe kits, SST
pointing exports and a TaLIS report. They are non-public. Models may be tuned to and
checked against them; the files and their contents stay out of the repository, and the
notes describe them generically. Check them before listing a run as missing.

6. **VisRad deposition: shadowing and beam totals.** *State:* `VisRadProfile` is the
   VisRad profile law. It agrees with VisRad to 0.1–0.3% rms per element on plates and
   0.45–0.7% on hohlraum walls ("VisRad probe kits (2026-10-27): results"). *Remains:*
   - VisRad rescales each beam's total after shadowing: to its power when every ray hits
     something (per-plane normalization is +0.2–2.2% high on steep foils), and to a share
     within ±2.5% of the geometric one when part of the beam hits nothing (grid not
     identified).
   - VisRad's shadow rays lean outward in proportion to spot position (a virtual source
     12–19 mm upstream); axis shadowing misses by ~30% rms behind occluders. An
     aim-dependent mask conflicts with an operator built once per geometry.
   - The leak through the LEH washer (≤ 7 × 10⁻⁴ TW/cm²) is not reproduced.

   *Data:* occluders ~2 and 10 mm upstream and a foil off focus; the Laser Beam Summary
   export of the spot-rotation workspace (`read_visrad_summary` end to end).
7. **Benchmark against VisRad.** *State:* the saved runs contain only VisRad's deposition.
   Their radiosity and view-factor sections are empty in every file ("VisRad benchmark
   beyond deposition: the saved runs"). With VisRad's deposition as the source, Briareus on
   VisRad's profile law differs after the same radiosity solve by 0.13–0.17% rms on the
   gold, −0.05 to −0.08% in the capsule's flux and ≤ 0.03 pp in P2, and all of the
   difference comes from the deposition (the drive level follows VisRad's per-beam totals,
   item 6). The VisRad view factors, visibility and closure are not yet compared.
   `compare_visrad(workspace, solution; elements)` compares radiosity, temperatures, drive
   and modes element by element once VisRad exports them (aa19750). *Remains:* in VisRad,
   compute the radiosities of one saved hohlraum run at ~0.6 mm and one at 0.1 mm, save, and
   export the Surface Element Properties of the gold objects and the capsule with the Power
   Summary; then run `compare_visrad`, confirm the export's index base and write a reader
   for the populated radiosity section. `write_visrad_kit` cases are needed only for what
   the saved runs lack: primitives with curved surfaces, a closed entrance hole, a
   temperature-dependent albedo. A 0.1 mm run's `:sampled` view factors (60k elements) take
   more than 5 minutes on 16 threads; the comparison there used `:centroid` visibility.
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
9. **Confirm against facility practice.** The items to confirm are the TaLIS 0.5° pivot
   per component (TCC or the mount), the side and shape of the 1ω footprint, which quads
   count as open for the 20° reflection rule, the NIF CMT XML schema, the SST preamble and
   the element-export layout ("To verify with real files at work"), and the NIF cone
   pointing of the GUI presets and the 7700 mm lens distance (unsourced in
   `beam_table.jl`, `beam_geometry.jl`).

## Physics against shot data

Items 10–15 follow the 2026-10-22 literature review against Kritcher et al. 2024
("Literature review against Kritcher et al. 2024"; `reference/litreview/`).

10. **The P2 swing (N221204, PRE 109 Fig. 10).** *State:* nothing is fitted to Fig. 10. The
    swing reaches −0.181 at 5.25 ns and turns positive at 6.19 ns (Fig. 10: −0.170 at 5.6
    ns, 6.67 ns). Every published input tried moves the crossing by −0.10 to +0.16 ns; ~0.3
    ns remains, which only the least constrained bubble inputs can carry ("The P2 swing",
    "Inner-cone SBS time histories"). A remeasurement on 874dcca without new data ("The
    swing and the bubble without new data"; `litreview/bubble_drivers.md`) found the swing
    unchanged since 2026-10-25. Lips held at a thermal-band depth (0.30–0.36 mm, an
    estimate) move the crossing by ≤ 0.02 ns, the wrong way. A tip slowed through the trough
    (a 0.5 ns delay) puts the crossing at 6.58 ns but deepens the minimum by 0.04, and
    N221204's trough contrast (16× total, 12× outer) lies at Izumi 2023's 15× ballistic
    criterion, so the slowed tip is not adopted. *Remains:* a HYBRID-E bubble speed and
    collision time (none found published); a thermal-band LEH radius or the warm lip's
    opacity; the rise (0.063 above Fig. 10 at 4.5 ns and 0.075 below at 5.0 ns, since the
    model falls into the minimum later and more steeply) and the foot at 3–4 ns (0.015–0.054
    above). The benchmark keeps its hand-set picket ratio, and re-pointing studies use the
    link relative to their reference (user, 2026-11-18). The link (46d88e6, §5 of the same
    section) is opt-in: `picket_fluence`, `picket_ratio` and `remake(bubble; picket_ratio)`
    take the bubble speed from the energy-weighted picket fluence of the outer cones on the
    wall. On N221204 a ±0.5 mm re-pointing of the outer cones along z changes that fluence
    by ≤ 1% and the crossing through the bubble by ≤ 0.013 ns; unsplit outer quads raise it
    19% (−0.19 ns). The benchmark's pointing gives a ratio of 1.147 against the calibration
    platform (1.1 when set by hand), on a proxy pointing and a picket energy read by eye;
    ratios between pointings converge to < 1%. A published bubble speed at fixed picket
    energy and a different spot overlap would test it. `bubble_prediction` (74e7ae3)
    reproduces that table; the picket's peak does not bound its mean. With the outer cones
    moved +0.5 mm in z, the peak changes by ×1.000 and the mean by ×1.010. Validation
    remains against Young's dynamic albedo (which needs N170827's delivered per-beam
    powers), Ralph 2018's 18.5% cone fraction (which needs the high-foot pointing and wall
    radius), the N210808 proxy and Callahan 2018's 39-shot curve. *Data:* bubble-tip or
    inner-spot histories (thin-wall or G-LEH), a HYBRID-E LEH-radius history, the 2DConA
    in-flight shell P2, the HYDRA multipliers behind Fig. 10, N221204's shell thickness at
    5–7 ns and per-cone backscatter history; a tip position during the trough.
11. **Shot-data calibration (CBET and effective parameters).** *State:* `n221204_cbet()`,
    `hybride_cbet()` and `hybride_cbet_iraum()` share one builder since d72ed90 ("Joint
    multi-shot CBET calibration"; "CBET without new data: the Δλ slope and the Te
    history"). The models' dCF/dΔλ, 2–3× steeper than HYDRA's cross-shot response, is
    not shown to be a model defect. On the one step with the pulse held (from N220919 to
    N221204) the models give a hot-spot P2 step of +10.7 to +13.6 µm, the measurement
    +11.1 µm and HYDRA's runs +3.8 µm, and no model term (the clamp, a local depletion,
    the plasma state, the outside zone) reaches HYDRA's slope. A Te(t) that follows the
    Dante drive removes the fitted Te's dip at 5.5 ns and lifts the peak clamp into
    HYDRA's range, but it doubles the loss and worsens the blind Fig. 10 check on three
    shots, so it is not adopted. *Remains:*
    - Fit effective parameters (cone multipliers, CBET multiplier and clamp, bubble
      transparency, wall speed, albedo) to measured in-flight and hot-spot P2/P4 and
      keyhole symmetry, with uncertainty, through a flux-P2 to hot-spot-P2 transfer.
      HYDRA's flat cross-shot dCF/dΔλ is no longer a target (user, 2026-11-18).
    - Ruof 2026's asymmetric rise cone-fraction response (the reduced CBET amplifies
      cuts and raises alike, 1.2–1.5×) remains open.
    - Te is at the 5.5 keV fit bound at 8 ns; a laser-power Te prior is not run.
    - Plasma-state tables as inputs, an inner-quad SBS proxy and the 23.5°/30° split
      remain. *Data:* N230729
    transfer curves; per-shot plasma states. Haines et al., Phys. Plasmas 33, 112712
    (2026), Fig. 9, is such a curve for the neutron image against a flux P2 fixed in
    time. It is flat for applied P2/P0 of −0.03 to +0.035 on the ignition cliff, so a
    linear transfer needs in-flight shapes, and three nominal repeats of one design span 0.046
    in the inferred applied P2/P0, which bounds what a same-pulse P2 step can calibrate
    ("Drive asymmetry inferred from post-shot capsule simulations: Haines 2026").
    *Remains from it:* a fixed-in-time equivalent of a drive history's P2 (a weighted
    time average of `drive_history`'s capsule modes; the weight to be chosen), to set
    beside the three inferred values as a check.
12. **Absolute energetics on N221204** (`reference/litreview/energetics.md`,
    `late_time_drive.md`; "Absolute energetics without new data: the early peak and the
    missing sinks"). *State* (coarse mesh, 90e28c4): the capsule share is 13.9% and the LEH
    share 24.2% (HYDRA 16.2 / 25.2); Dante 1 at 8 ns gives 289.7 / 312.9 eV over the initial
    / open area (310.1 measured); the whole-field view equals the view through the hole on
    N221204 within 0.1 eV. *Remains:*
    - The rise through the peak: Dante rises 72 eV from 5 to 8 ns, the model 47 eV over
      the initial area and 57 eV over the open area. The "early peak" (13–20 eV above
      Dante at 5.0–5.5 ns) is the open-area face of this one error; on the Dante
      convention the model is within 6 eV there and 0.07 ns early. No existing input
      within its documented range adds more than 4 eV to the rise. *Data:* a
      thermal-band LEH radius, or HYDRA's loss history through the peak. A second
      public history to check the rise on: N230729's laser power and Dante temperature
      (Haines et al. 2026, Fig. 1(a)); its pointing is not public, so the check would
      use a mock-up pointing.
    - The drive after the laser falls: 30 / 38 / 47 eV low at 8.5 / 9.0 / 9.35 ns (37–54
      eV at bang on three shots); the missing ~250–450 kJ is in the plasma, with an
      unpublished return rate. *Data:* HYDRA's energy inventory through the coast, a
      post-shot synthetic Dante, Dante spectra between laser-off and bang, gated LEH
      images after laser-off.
    - The capsule share is 2.3–4.2 points below HYDRA. It and the LEH share move against
      each other under every input except a lower wall loss, which the DU calibration
      rules out. *Data:* the blow-off opacity or HYDRA's absorbing radius; the shell
      thickness at 5–7 ns; the radius histories.
    - Some sinks are sized from published numbers and not modeled: the He gas (4.5% in
      HYDRA, or 80 kJ and 3–4 eV of drive at every time), LEH windows, diffusely scattered
      laser light (5–9% of the light not absorbed at the first hit; one bounce sends 72% of
      it back to the wall) and specular glint, each under 2 eV.
    - The LEH difference of the Zylstra pair is 2–4 eV too large (N191118's pointing and
      cone fraction are unpublished), and the hot-spot velocity needs per-quad delivered
      powers.
13. **Spectral radiosity and diagnostics** (`reference/diagnostics-design.md`).
    - Dante needs measured Ni XRD and thin-filter transmissions; the angle scale,
      roughness and surface layer of the SiO₂ mirror (assumed or fitted; extrapolated
      below 2.2 keV); the within-group spectral shape (finer groups); the hardware and
      fields of Dante 2 and of Omega's Dante, which look through the hole because their
      fields are not public (`dante2()` assumes Dante 1's channels); time response and
      per-channel uncertainty, and an uncertainty band from `dante_unfold`; recording
      modes other than `history`; and field-of-view checks against Chen 2024, the BAH 27%
      ring occlusion and May 2015's cosine ratio. Sessions saved before 46d8a76 keep
      their through-hole Dante 1.
    - Diagnostic definitions phase 5 covers the SXI soft channel of the imagers,
      per-pixel gate times and near-pinhole perspective, and whether SXI and G-LEH 3–5 keV
      images need a non-thermal spot model. It also covers backlit radiography in the GUI
      (definitions are listed and can be ticked but record nothing), with ARC, opacity
      and streaked radiographs and a continuum, satellites, photon noise and scatter in
      the radiograph, and streaked VISAR and SOP. Only three diagnostic kinds are
      registered.
    - The radiosity needs a measured soft spectrum of spots in a cavity; the angle
      dependence of the M-band (Dewald 2020) and its emission from the outer bubble; the
      emission of the closing LEH plasma; a time-resolved M-band P2 against Dewald's
      Fig. 10; and the M-band split in the `CapsuleDrive` handover. `_ALPHA_M` rests on one
      observable (Dewald 2020 F_M). The spectral solve is Float64 only, so M-band and Dante
      objectives have no gradients. Block-GMRES and thermalized solves raise an error on
      non-convergence.
14. **Implosion metrics and UQ.** Time-dependent error budgets and multi-shot trending
    are missing, as are focus, spot-size and timing errors beside pointing, power and
    positioning (linear propagation, Monte Carlo as the check). Groups over both
    hemispheres are refused. `burn_amplification` refuses Y_amp 40–80 (a gap in the
    published fits). `flat_kernel` stands in for Young's mode-1 Green's function. The
    sign convention from Yₗₘ to Young's Cartesian vector needs a check before vectors are
    compared (`litreview/asymmetry_transfer.md` §7).
15. **Backlighter physics and retune.** Open items are the low-Z end of He-α CE (one Ti
    point; no V/Sc data at 351 nm); pulse length and spot size as variables (CE ignores
    both); and continuum and higher lines in the filter pass band, detector response and
    the transmission of the object. The retune (core `retune_drive` with
    `FoilIrradiance`/`FaceFlux`, level held; "Retuning faces and foils; the held level")
    reports peak-to-mean but does not optimize it. It holds the albedos of the reference
    (linearized) and retunes powers; the pointing is not retuned.
16. **Unconverted light over time.** *State:* `unconverted_history` and
    `irradiance_history` exist, with `specular_check` and `preplasma_check` against the
    User Guide's 10⁸ W/cm², the *Irradiance* chart of the GUI's Time tab, and the
    prepulse in the Tr of the drive history ("Unconverted light over time", "The prepulse
    in the drive history"). *Remains:*
    - pointing or geometry that changes in time (a deposition per configuration);
    - a fluence-based damage criterion beside the irradiance threshold, and where the
      specular 1ω reflection goes (the check tests the threshold only);
    - the footprint's phase-plate edge blur (a flat-top square now).

    *Data:* the unconverted fraction vs. pulse shape and the prepulse level; plasma,
    damage and specularity thresholds of Au, Al and Si at 1ω and 2ω.
17. **Per-color sources and a materials library.** *State:* `solve`, `drive_history`,
    `flux_response` and `retune_drive` take per-color sources, and the drive, histories and
    capsule retune of the GUI take the ratios and prepulse of the Drive panel ("Unconverted
    light: absorption and conversion vs. wavelength"). *Remains:*
    - the M-band of unconverted light (the M-band stays the 3ω light's);
    - data on the conversion of overlapping colors (each color at its flux is an
      engineering default) and on whether the 3ω intensity table applies to 1ω and 2ω
      (Mead's Iλ² scaling would shift it);
    - materials: `material_library()` holds the defaults that have a source, each quantity
      with its source and uncertainty (154a4c5; "Material defaults and their sources").
      Since aa7135f (user decisions, 2026-11-18) `gold()` and `uranium()` convert at 0.9
      (Young 2020's static NIF setting; the capsule Tr of the default hohlraum rose 4.66%,
      with P2 unchanged), the GUI's hohlraum is `gold()`, and foils take `foil_material()`:
      gold at 0.8/0.8 (Dewald 2008 at 10¹⁵ W/cm², Lindl's wall loss), the mid-Z backlighter
      foils on a suggested range (conversion 0.1–0.6, albedo 0.3–0.8, which are not
      measurements). Left: a measured 3ω broadband conversion for the mid-Z foils and the
      ablators (CH, HDC, Be), and albedos specific to CH and Be (none found public); the
      capsule conversion 0.3 (user); a Help entry for the foil ranges if the GUI lists foil
      materials. Since 2026-11-19 the default efficiency of `Conversion()` is 0.9, as for
      `gold()` ("Defaults changed (2026-11-19)").

## Re-pointing studies

Items 34–36 come from the user's review of 2026-11-18, after a study that drops a
bundle and re-points the remaining quads to hold a reference drive shape (the public
analog is the dropped-bundle tutorial, `docs/literate/dropped_bundle.jl`).

34. **Re-pointing studies on the public API.** *State:* the public API has `subset`,
    `GroupPointing` (`relative = true` for a rigid shift per quad),
    `asymmetry_objective(; reference)`, `ApertureKeepOut` and `mode_uncertainty`. Since
    0594191 it also has held parameters and power scales per group (`HeldParameters`,
    `PowerScales`), `remake(::Scenario)`, NIF beams on VisRad's plates in the quad frame
    (`visrad_nif_spot`, `spot_frame = :quad`), quads to beams and back (`nif_beams(quads)`,
    `nif_quads(beams)`, an SST table through `read_visrad_sst(; template)`) and
    `pointing_summary` ("Re-pointing studies on the public API"). The gated objective and
    the penalty terms belong to items 36 and 35 (`drive_objective` with `mode_gate`;
    `keepout_penalty`, `light_ceiling`). *Remains:*
    - Overlap and reach estimates need geometry queries. `port_basis` and `visible` are
      exported (2acd1e9). Left: `closest_hit` as a method on the target (it takes the
      internal `BVH`), and the beam section at a plane as a `beam_section(beam, point,
      normal)` (`_beam_at_plane` is shaped for CBET).
    - `subset` keeps the parameter names, values and holds of a parameterization,
      including the groups of a `GroupPointing` since 2026-11-19 ("Defaults changed
      (2026-11-19)"), and the wrappers take named values (59b366d). Left: reading a
      wrapper's values off beams.
      `GroupPointing` points by quad by default since 33cf879, with a quad's split as
      a parameter (`:split`) and `unit = :beam` or `by = :beam` for beams pointed apart.
    - `PhasedScenario` holds one pointing over several phases ("One pointing over
      several phases of a pulse"). Each `Phase` has a separate target or scenario (for
      example a `remake` of one shared scenario with the albedos of `history_state`) and
      separate beam powers, per beam or from pulses at a time. The parameters are those
      of the shared pointing, followed by a power scale per (key, phase) for the phases
      not fixed. `evaluate` gives one evaluation per phase (`PerPhase`), which the
      objective and every constraint receive, and `in_phase` binds any constraint or
      objective of one evaluation to a phase. `bounds(...; power_max)` caps every beam in
      every phase at the laser's limit, and `optimize_pointing` runs on it with
      constraints, `from`, `time_limit`, `phase1` and `relaxation`. It replaces the four
      workarounds of a multi-phase study (beam copies under one name, a call-counting
      `by`, a carrier scenario, a cache behind a lock). It reproduces that construction
      to rounding at its saved parameters (552 of them) and takes 0.45 of its time
      per gradient. *Left:* a scale shared across phases; a runnable example in a
      tutorial (over a minute of build time); the stages of the study rerun on it.
    - The quad stage of the multi-phase study was rerun on `PhasedScenario` (scratch
      runs, 2026-11-19). The start reproduces that of the study to every printed digit,
      and the run reaches the same design region in 64% of its evaluations. What the
      study still wrote by hand is now API ("The multi-phase study's hand-written
      pieces as API"): `drive_ceiling` beside `drive_floor` for a band;
      `mode_gate(...; fraction, floor, shift)` (a gate at a fraction of the noise
      floor, the floor as a lower bound under gates in physical units, a target
      shape shifted in a_l/a₀ with the model and time check kept) and `shift_modes`;
      `clearance_floor` over every phase without naming one; `power_scales`,
      `phase_powers` and `power_scale_table`; `remake(scenario; history, time)` and
      `Phase(...; history, time)` for a phase straight from a drive history;
      `pointing_summary(...; gates, floor)`; `picket_ceiling(...; reference_names)`
      with a reference on another target for the integral measures (`:peak` and
      `irradiance_ceiling` still refuse one: they follow one element). The study's
      start rewritten on these reproduces the hand-built objective and its 212
      constraint rows to 1e-16. `picket_ceiling` with `lpi = true` now checks the
      laser–plasma model and time of the reference on the reference target too
      ("Inner-run options, a wall-loss scenario warning and laser delivery figures").
      *Left:* `shift_objective` on a `PerPhase`. No run of the study or its variants converged within 40 minutes
      each: an inner run can spend hundreds of evaluations on 1e-4 of the objective
      while limits stay violated (the iteration cap does not bound evaluations). An
      opt-in cap on one inner run's evaluations, `inner_evaluations`, now ends such a
      run and updates the multipliers ("Inner-run options, a wall-loss scenario warning
      and laser delivery figures"). *Left:* the study rerun under it.
    - A design ladder was run on the same study (scratch, 2026-11-20). What its script
      still wrote by hand is now API ("A ladder study's hand-written pieces as API
      (2026-11-20)"), and the script, rewritten on it, checks each piece against the
      original:
      - the pointing-noise floor as `mode_uncertainty(...).floor`, and gates at a
        multiple of it per degree, `mode_gate(reference, u; fraction = [1, 2, 2, 2])`;
      - `power_scale_objective`, the beam-weighted mean of `(s/s0 − 1)²` over the power
        scales of a `PhasedScenario`;
      - `remake(::PhasedScenario; deposition_kw)`, every phase's operator rebuilt (one
        per mesh) to check a finished design on a finer one;
      - `mode_departure(...; m = :zero | :nonzero)` and the summary's
        `departure_zero`/`departure_nonzero`;
      - `pointing_summary(...; moved, moved_by)`, the quads moved past given shifts.
      *Left:* uniform start scales at a fraction of a reference's drive (exact only for
      a drive linear in power; a few lines of a study).
    - A pulse-shape budget is built. `pulse_shape_groups` groups per-key multipliers
      over the phases into a given number of shapes by k-means on the log multipliers
      per class, and it chooses the split between the classes by fit. The fit is exact
      by enumeration or 1-D dynamic programming where the size allows, with seeded
      restarts otherwise, and it reproduces the grouping of the study. `shape_pulses`
      and `scaled_pulse` write the shapes out as `PulseShape`s for `write_pulse_table`
      (with `peak`, a peak at the laser's limit is left unscaled), and the `assignment`
      of the grouping is used as the `by` of `PowerScales`. `tolerance` makes the count
      an upper bound: the fewest shapes whose fit is within an rms in log multiplier.
      NIF's public figures (energies within 2–3% and beam-to-beam balance better than
      5%) give the scale of that rms but do not set a default ("Inner-run options, a
      wall-loss scenario warning and laser delivery figures"). Since 2026-11-20
      `scaled_pulse` and `shape_pulses` check a table against the laser's limit
      (`power_max`). A violation raises an error that names the time and the excess,
      and nothing is clipped. `pulse_ceiling(scenario; pulses, times, power_max, peak)`
      holds the written pulses within the limit between the phase samples. It imposes
      the exact condition, which is linear in two neighboring scales, as one component
      per vertex of the lower convex hull of the points between two phases.
      `bounds(...; pulses, times, peak)` holds the conditions on a single scale ("Least
      violation, capped inner runs and phase 1's budget", § 4a).
      *Left:* the tied re-optimization with per-phase keys on a `PhasedScenario` as
      one call; phase weights in the grouping (user).
    - The smallest re-pointing that meets given gates: `shift_objective`,
      the weighted mean squared aim shift from reference beams, and
      `PointingSummary.shift_rms`.
    - Evaluation cost across versions ("Evaluation cost, measured (2026-11-19)").
      VisRad's profile law (5.6× per deposition) and the refined, then smoothed,
      aperture keep-out (about 10× v0.1.0's gradient) had made an evaluation of a
      study-sized scenario 5× that of the 2026-10-22 tree on one thread (2.7–2.9× at 8
      threads) and 1.2× that of v0.1.0 (1.3–1.6×). The constraints of a multi-phase study
      cost as much as its objective. The fixes, in order:
      - The keep-out is computed once for phases whose beams have the same aims and
        spots, its dual search runs on values, and `deposit` skips zero-power beams and
        points beyond VisRad's cut, with identical results (an evaluation at 0.52–0.62
        of 0.3.0's, 0.7–0.8 of v0.1.0's).
      - The keep-out takes exact interval minima with chain-rule partials ("The aperture
        keep-out's exact interval minima (2026-11-20)"). Its gradient costs 0.1 of
        before and 0.53 of v0.1.0's, with clearances within 1.2e-7 µm.
      - The far screen is inlined, twice as fast again per deposition at a 1.2e-13
        relative change of Float64 VisRad depositions (an evaluation at 0.37–0.43 of
        0.3.0's).
      - Inside a constrained run, the picket ceilings take the powers of the beams from
        the deposition of the evaluation, and the group sets they read are kept in the
        same pass.
      - VisRad's screens are tightened (a bound on the weighted sum of the law, runs of
        points by bounding ball, no value screen on duals), and the coordinates of a dual
        beam share their partials. A VisRad deposition then costs 0.5 of that with the
        inlined screen.
      - The moment and Y_lm products of an evaluation are formed once. With every fix, a
        study-sized evaluation costs 0.18–0.21 of 0.3.0's.
      - The phases of a phased scenario on one mesh are deposited in one pass, with the
        profile evaluated once per mesh and scaled by the powers of each phase. The
        gradients of a multi-phase study take 0.37–0.47 of the time before on one and 8
        threads, with results at rounding.
      - The same coordinate reuse applies to the f-cone model (case A's gradient −14%,
        results at rounding). Since "The dropped-bundle tutorial's long runs under
        rounding (2026-11-20)", the dropped-bundle tutorial quotes the ranges of its long
        runs.

      *Left:* the rim search of the keep-out on values, repeated in every chunk of a
      gradient, is now its whole cost (1.8 ms for 46 quads on one thread). A cache across
      chunks would cut that to about a fifth, at the price of state in the keep-out. The
      envelope-radius test behind fix 1 costs about 3 ms a dual pass on 184 beams over
      five phases.

    The tutorial's section "A re-pointing study, end to end" uses only exported names for
    these steps (a96792f): quads held, a power scale per cone, the gates, the ceiling, a
    hardware keep-out and the light on a part as constraints, and a per-beam polish as a
    second call on the beams of the first result (`remake(scenario; parameterization)`).
    A schedule type was tried; it saved no lines over the two calls and hid the check
    between them.
35. **Capture and clearance rules in the optimization.** *State:* `clearance` and
    `talis_check` check a finished design on the worst case over pointing, positioning and
    focus errors (reports, not differentiable), and `placement_sweep` runs the capture and
    clearance checks at each step. For the optimizer, since 2c1825c ("Capture and clearance
    rules as optimization quantities"), these exist: the keep-out terms as public quantities
    (`keepout_clearance`, `keepout_penalty`, `direct_fraction`); the light on named parts
    for floors and ceilings (`power_fraction` with `clipped_fraction(ev, names...)`);
    `HardwareKeepOut`, the exact worst-case distance from a beam's envelope to named
    hardware, smooth in the beams and within 1e-6 mm of `clearance` (lower by up to 0.4 mm
    where `clearance` skips points hidden from the beam, as `talis_check` does not); and the
    TaLIS rules that carry a margin, each within 1e-7 mm of `talis_check` on a mock-up
    (every row, since `_ellipse_distance` became exact next to the major axis, 98efc9a): the
    3ω clearance (`HardwareKeepOut` with `tolerance`, the 0.5° term inside the
    minimization), 3ω capture (`CaptureKeepIn`), incidence (`IncidenceRule`) and the 1ω
    clearance (`Omega1KeepOut`). Each has a `keepout_clearance` method, so each can be used
    as a penalty or as a constraint (item 36). *Remains:*
    - a loop that runs `clearance` and `talis_check` on the result and tightens the
      minimums until the review passes (item 21);
    - capture fractions per beam (today's are totals over the beams);
    - the 1ω clearance by an exact search (branch and bound now: ~3.5 ms per beam, three
      depositions' worth for 48 quads); an incidence with a gradient on curved parts
      (the angle is constant per facet, so the gradient is zero away from facet edges);
    - 1ω management, the diagnostics-to-1ω rule, the TAS stay-in rule and the LEH review
      tier have no quantity;
    - `ApertureKeepOut`'s clearance is a smooth minimum over the rim since 2026-11-19
      (`softness`, 1e-4 mm: at most 0.18 µm below the exact minimum on NIF quads; the
      largest slope change over 144 scans 0.002 instead of a kink of 0.21; "Defaults
      changed (2026-11-19)"). Left: the residual slope change where one rim interval
      holds two local minima of the distance.

    *Check:* each surrogate against `clearance` or `talis_check` on the worst case, as
    the tutorial does for the aperture keep-out (agreement to a few micrometers); a
    re-pointed design that passes `talis_check` without repair by hand.
36. **Constrained re-pointing stalls at its bounds.** The study's objective matches the
    capsule modes of a reference drive while holding the peak laser irradiance over the
    wall's elements at or below the reference pointing's (high irradiance drives the gold
    bubble, item 10); started on the final tolerances, its penalty runs ended with
    parameters at their bounds and the objective unmet, and stages set by hand were the way
    around it. *State* (7045105, "Constrained re-pointing: an irradiance ceiling and mode
    gates"): `optimize_pointing(...; constraints)` solves inequality `Constraint`s by an
    augmented Lagrangian around the same L-BFGS-B, and the result reports each constraint
    against its limit and what stopped the run (`PointingResult.constraints`, `outcome`);
    `irradiance_ceiling` (an area-weighted 64-norm as the smooth stand-in, `peak_irradiance`
    as the reported quantity, a `weight` per beam for outer-cone or picket-fluence
    variants), `mode_gate` and `mode_departure`, `drive_objective`. On the dropped-bundle
    tutorial's case one call meets the four gates under the ceiling from three starts
    (488–765 evaluations before 320335b, which made constrained runs converge at the final
    inner tolerance only, and 1.5–6 times as many since), at 0.4% more drive than the
    penalty form reaches, staged by hand or not. On the study itself (scratch runs,
    2026-11-18) one call from the as-shot start meets every gate and the ceiling inside the
    box, every constraint active, at the drive of the earlier design without a ceiling; the
    study's penalty form, run the same way, stalls with the peak 7% over the limit and a
    parameter at its bound (the reward for symmetry inside the gates balances the ceiling's
    penalty), so the stall came from the penalty form and not from infeasibility. A per-beam
    polish (`parameters(p, beams)` hands a quad result to per-beam shifts) moves beams by
    micrometers and the drive not at all. Since bf79dc2 a run continues from a result with
    its multipliers (`from`), stops at a `time_limit` with its best design, reports the
    drive ratio and departures in `pointing_summary`, and a ceiling refuses a reference
    evaluated with another laser–plasma model. *Remains:*
    - On the study's variants (scratch runs, 2026-11-18) one call converges where the
      constraints are compatible (the base case, one CBET variant) and runs to its
      time limit where they are not: with the other CBET variant and with windows the
      gates conflict with the entrance-hole keep-out at the reference's clearance, not
      with the ceiling, and with only the outer quads re-pointed the gates conflict with
      the ceiling. Those runs left these items:
      - `VisRadProfile`'s hard cut (the profile dropped below 1e-3 of its peak) makes
        every deposition-derived quantity jump as a sample point crosses it, and the
        line search parks on a jump. `VisRadProfile(; taper = 1)` (e77cb5c, opt-in)
        tapers the last e-fold smoothly on the same support: no jumps, derivatives
        of the mode gates right to 1e-6 where they were off by 1e-3–1e-2, at a cost
        of at most 0.011 points of rms against VisRad's plates and 7e-5 of a wall
        peak; the stalled case then runs 8 outer iterations instead of 1. Since
        2026-11-19 `optimize_pointing` optimizes a hard-cut `VisRadProfile` scenario
        on the tapered profile and reports under the hard cut (`smooth_cut`;
        `tapered_value`; "Defaults changed (2026-11-19)");
      - the entrance-hole keep-out's penalty is not smooth at such a design (the
        objective's derivative is off by 2–5% against differences while the drive,
        ceiling and gates agree to 1e-6); this is being diagnosed;
      - a verdict on conflicting constraints: since 83c1f55 `phase1 = true` ends an
        infeasible run with `:infeasible` from a first phase that minimizes only the
        violations (on the study's two infeasible cases in 15–34 s instead of
        a run to its limit), and `relaxation = true` names how far each conflicting
        constraint must be relaxed, one at a time, for the rest to hold (`Relaxation`; on
        the study within 2% of the gate multiples found by hand). `phase1` is on by
        default since 2026-11-19 ("Defaults changed (2026-11-19)"; `relaxation`
        stays off). Since 2026-11-20 phase 1 has a separate budget (not the inner runs'
        `iterations`, which had cut it short of a verdict on a study variant) and
        names what ended it ("Least violation, capped inner runs and phase 1's
        budget"). Left: relaxations cost one solve per constraint and can outlast
        the time limit; the stall rule's defaults (1% over 100 evaluations);
      - a keep-out given as a penalty gives way silently (by 0.07–0.37 mm there) once
        the multipliers grow, while `clearance_floor` as a constraint held;
      - an absolute irradiance ceiling with `PowerScales` is met by cutting power:
        `drive_floor` (83c1f55) holds the drive beside it;
      - a design on its gates exceeds them under another beam model (by 7–12% of a
        gate with the f-cone model against the VisRad profile; within 0.5% across
        mesh and operator resolution): a margin on the gates is the designer's choice.
        For scale (Haines et al. 2026, read from their Fig. 9): a flux P2/P0 fixed in
        time does not cost yield within −0.01 to +0.015 on one ignition design, and nominal
        repeats of it span 0.046, against gates of tenths of a percent from the
        pointing-noise floor. The floor is an aggressive gate (user, 2026-11-18):
        variation in laser delivery moves the shape considerably more but cannot be
        readily quantified, and the tolerances of a re-pointing study can be relaxed
        without a loss of performance within that uncertainty. The `mode_gate`
        docstring and the tutorial now present the floor as a lower bound and a
        tolerance per degree in physical units as the usual gate; the tutorial adds a
        least-re-pointing run at 0.5% per degree. What a relaxed gate buys is
        measured on the case of the tutorial ("Mode gates in physical units"). The
        least re-pointing falls from 93 µm rms at the floor to 63 µm at 0.5%, and to
        nothing once the gates reach the departures of the start. The largest drive
        rises from 0.962 of the reference at the floor to 0.967 at 0.5% and 0.975 at
        five times the floor, with l = 2 active throughout and no end to the returns
        in that range. *Remains:* the
        tolerance per degree (0.5% in P2/P0 is the modelers' repeat uncertainty;
        l = 1, 3 and 4 have no published number) is the user's. NIF's public
        beam-to-beam balance (better than 5%) sizes laser delivery's share through
        `mode_uncertainty(...; power)`: 3% per quad raises the tutorial's floor by
        40% at l = 1 and 21% at l = 2 ("Inner-run options, a wall-loss scenario
        warning and laser delivery figures");
      - a drive run under gates reported zero gate multipliers at convergence. The cause
        was the first-order update after inner runs that ended on their value tolerance,
        at penalty weights of 1e2–1e6. Now the reported multipliers are fitted
        to stationarity over the near-active set (non-negative least squares),
        `bound_sensitivity` follows from them and `PointingResult.stationarity`
        shows what they leave; designs are unchanged. The display notes an exact
        ceiling whose row shows `:satisfied` with a positive multiplier (its
        stand-in's, at the working limit). On the tutorial's run they
        name the ceiling and l = 2 but leave a third of the drive's gradient: the
        design ends short of stationarity. The inner value tolerance is a keyword,
        `inner_reltol` (default 1e-10; at 1e-12 the run reaches a stationarity of
        1.9e-5 for 2.1 times the evaluations and the same drive within the flat
        family's 1e-4), and a result continues tighter through `from` ("Inner-run
        options, a wall-loss scenario warning and laser delivery figures").
    - a constrained run stopped by its budget reports the exact peak as it stands
      (`exact = true` is the default since cd6f084, and its rounds run after
      convergence): a capped polish needs convergence or `from` to meet the ceiling.
      The tutorial's per-beam polish stays capped at 120 evaluations (user,
      2026-11-19; "Defaults changed (2026-11-19)"), and its text states what the cap
      leaves and that a design to hand off runs to convergence;
    - bound sensitivities overstated the gain of a wider box on the study (about fiftyfold).
      The cause is the same: a design short of stationarity on a flat family of designs,
      which `stationarity` now shows. They match the closed form and restarts with a bound
      moved on a closed-form problem; on the tutorial, restarts with a bound widened gain
      nothing, while narrowing costs 1.4–6.4e-4 per mm against a reported 1.1–1.9e-4 (at an
      inner tolerance of 1e-12). *Remains:* a direct estimate (a restart with the bound
      moved, through `from`) built in, if wanted;
    - The ceiling the study means is the maximum irradiance at the picket, and the aim
      is prediction of the bubble, not a limit standing in for it (user, 2026-11-18).
      Since 74e7ae3 there are `picket_ceiling` (the outer cones' peak picket irradiance on the
      wall, held exactly, or its energy-weighted mean), `PicketLinkedBubble` (the
      bubble's speed follows the pointing inside an evaluation) and
      `bubble_prediction` (picket ratios, tip speed, collision and the P2 swing of a
      design against its reference). On the study, with picket powers assumed from a
      public pulse, the peak-power and picket-peak ceilings conflict by 1–4% (each is
      exceeded by the design that holds the other), neither bounds the mean fluence,
      and the predicted bubble change of either design is under 0.2% in tip speed;
      on a proxy model the re-pointed deposition itself moves the swing's crossing
      by +0.2 ns, ten times the bubble's part. Left: the study pulse for the
      picket powers; which measure to hold, or the time-dependent swing itself in
      the objective (user); the prediction on the study geometry in a
      time-dependent model.
    - The published bubble scalings use the outer cones' picket fluence on the wall (item
      10). The quantity of the scalings is the outer picket energy, net of backscatter, over
      a fixed spot area (Callahan et al., Phys. Plasmas 25, 056305 (2018); Ralph et al.,
      Phys. Plasmas 25, 082701 (2018); Izumi et al., Phys. Plasmas 28, 022706 (2021):
      "spatially averaged outer picket fluence over the first 1 ns"). It is not the
      energy-weighted ∫F²dA/∫FdA of `picket_fluence`. The two quantities differ on uneven
      power scales. On a dropped-bundle study, holding ∫F²dA/∫FdA at the value of the
      reference cost 16% of the picket drive, while a design at the published quantity of
      the reference gives 1.10× on it (scratch runs, 2026-11-18). Now ("The published picket
      fluence, the trough and the inner-beam corridor") `picket_irradiance` and
      `picket_ceiling` take two more measures. `measure = :average` is the outer power
      absorbed on the named parts over a fixed area, in any phase such as the picket, trough
      or foot (`lpi = true` for net of backscatter; Callahan's 53 mm² as the default area).
      `measure = :corridor` is the outer fluence under a smooth weight around the
      centerlines of the inner beams (1.3 mm with a 0.3 mm edge, both engineering defaults).
      `picket_hot_elements` lists the hottest outer-cone elements with their distance to the
      nearest inner centerline. On public cases uneven per-quad picket powers of 0.2 rms
      raise ∫F²dA/∫FdA by 5–7% at the published quantity unchanged, quad crowding by 19%,
      and z shifts move the corridor by 20% per 0.5 mm. The exponents are in the notes (core
      P2 as fluence^0.5, bubble position as fluence^0.2; +10% of fluence costs 5–7 µm of
      core P2). *Remains:* which quantity `PicketLinkedBubble` reads is settled for now
      (user, 2026-11-18): the energy-weighted one, which reproduces Izumi 2021's split-quad
      shot where the fixed-area one cannot see it (the measure is `:weighted`, renamed from
      `:mean`); a published bubble trajectory at fixed energy with uneven quad powers or a
      re-pointing would test the choice; whether the corridor's mean or its peak stands for
      Izumi 2023's spikes;
    - For the inner budget, `iterations` caps each inner run, and when set high it leaves
      the multipliers unchanged for most of a run. On a 276-parameter multi-phase study,
      `iterations = 1000` completed one outer iteration in an hour (the peak 18% over its
      ceiling at the time limit); `iterations = 100` completed six in 40 min and met every
      limit but one, within 1%. The measurements ("Re-pointing tools: the inner budget, the
      smallest shift, a pulse-shape budget") show that a cap that starts at 10 and doubles
      per outer iteration saves 18–56% of the evaluations on three public cases but costs
      0.16% of drive and 20% more evaluations on the tutorial's first stage and breaks its
      capped polish; a cap proportional to the parameters and a stop on a stalled violation
      are worse. The default stays, and the docstring gives the rule (about 100 for hundreds
      of slow parameters started far from the limits). Since then an opt-in cap on the
      evaluations of each inner run exists, `inner_evaluations` (default none). On the
      coarse dropped-bundle test case 100 converged in 766 evaluations against 1,044 (one
      default inner run took 795) at 6e-4 less drive, and 30 needed more than 30 outer
      iterations ("Inner-run options, a wall-loss scenario warning and laser delivery
      figures"). *Left:* a converged run of the study under the doubling schedule or a cap
      before any default changes (user). Since 2026-11-20 a capped inner run that ends above
      the limits no longer raises the penalty (on a study variant every capped run had,
      ending the run at the penalty's ceiling after nine, and on a feasible closed form with
      a false verdict), and a run that ends above its limits returns its iterate of least
      violation, with `least_violation`, `last_feasible` and `capped` in the result ("Least
      violation, capped inner runs and phase 1's budget");
    - the engineering defaults `feasibility = 1e-3`, `p = 64`, `outer = 30` (user review);
    - the scales of the ready-made constraints (`clearance_floor`: 1 mm; `light_ceiling`
      and `light_floor`: the limit, at least 1e-3 of the delivered power) are engineering
      defaults (user review).

    *Check:* the study's objective reaches its gates under the exact peak in one call.

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
    - Profiles: windows on cones and across breaks, fillets on a profile; the tentless
      membranes; NIF presets re-aimed when a hohlraum part's profile is switched in
      place (the frustraum pointing applies on load, New target and Insert only).
    - Roles: several observed surfaces (decision 3, "later"); roles on VisRad and CAD
      pieces (only hohlraum-part pieces hold their own); an observed washer, cylinder,
      cone or flat CAD face; per-cavity albedo and energy series; a Dante field centered
      on a part from the GUI; VisRad primitive pieces have no frames.
    - Triad: several selected parts.
    - Tree: drag-and-drop into assemblies, duplicate, expanded nodes saved.
    - Removing a part keeps the survivors' places when their remaining placement steps
      are relations; one left with relative steps (Offset, Rotate) keeps its placement
      fields and may move (`Controllers.remove_parts!`).
23. **GUI: workflow gaps** *(survey)*.
    - Modes stop at l = 4 in the tables, weights and retune.
    - Beams: spot shape, f-number, wavelength and port are saved but not editable; aims
      are Cartesian only; VisRad profiles and focal lengths only from an import; no OMEGA
      EP preset (choosing OMEGA EP empties the table for an import); applying SST or a
      beam summary needs a preset first.
    - Starters by facility (user review, 2026-11-13): File ▸ New target lists them under
      NIF and OMEGA and loads the facility's beams, but OMEGA has only the backlighter
      foil. *Data:* published dimensions for an OMEGA-scale hohlraum and capsule (none
      found in two searches; the research notes hold Nova's scale-1 cylinder only); a
      direct-drive starter waits for a direct-drive analysis in the window (item 19).
    - Quads as units (user review, 2026-11-13): a split quad keeps one row, one aim
      point (aim r and the tangential offset in the quad's port azimuth, `GroupPointing`'s
      `unit = :quad`: the rows, the 192-beam preset, the optimizer) and a per-quad power
      (the rows; a pulse table's cone and quad columns, `Controllers.PULSE_UNITS`; the
      retune's result), and a click on its beams selects it while its row is closed.
      Left: in the core, `with_pulses`,
      `cone_fraction_pulses` and the pulse tables are per beam (the N221204 benchmark
      fires 192 beams on them): the per-quad reading is the window's;
      a quad's four beams pointed apart ("split quads" in the pointing sense) have no
      quad-level display of their spread (the row shows the mean).
    - Time: no pulse editing; no Time node in the Analyses tree (its prepulse and
      threshold settings sit in Drive).
    - Views: Look from port lists beam ports only; *Save view as image* grabs the main
      canvas only; a beam within the pick distance wins over a nearer surface.
24. **GUI: Help, tooltips and layout** *(survey, and the user review of 2026-11-13; the
    rule: a tooltip is a one-line hint, an explanation belongs in Help, written as short
    structured entries (`qml/common/help.js`), never page paragraphs)*.
    - Controls outside the survey's list still without tooltips: Target, Component,
      Frames and placement-step properties, the model tree's dialogs, Look from port.
25. **GUI density and layout (user review, 2026-10-29).** The aim is Office-ribbon
    density, but not all the way, with the 3D view first. *Left:* the next pass on what
    still takes room once the user has used it. The principle is to keep frequent
    commands near at hand and explanations in tooltips and Help, and to give the 3D view
    and its inspection the most room.
26. **GUI responsiveness** *(survey; measured 2026-11-03 through the window's shell
    offscreen, on the default design, 5.2k elements, and a CAD mock-up of 40 solids, 10.8k
    elements)*. *Done:* `record!` and `can_undo` answer from a count of the session's
    observable changes between edits. Before, a snapshot took 0.3–0.6 ms and its comparison
    0.25–0.5 ms at every event, and twice per refresh with no undo step; a refresh took
    1.8–2.0 ms then and takes 0.4–0.9 ms now. An undo step whose target section is the
    design's keeps the target, and restoring the view no longer makes a sector cut a plane
    cut on the way (two redraws of the 3D view each time); on the mock-up a view undo took
    88 ms and now takes 6 ms, a beam undo 126 ms and now 25 ms, and a part undo 188 ms and
    now 82 ms, each with its refresh. *Left:* the model tree and every panel are rebuilt
    each refresh (the tree in 0.1–0.3 ms and the panels together in under 0.3 ms, not worth
    caching yet); a `restore!` of an unchanged CAD design (a session reopened) still takes
    13 ms in `_try_build!` (each part's cache key reads its file's time stamp) and 17 ms in
    `set_beams!` (the beams redrawn); a part undo rebuilds the target (CAD designs depend on
    the 6-entry target cache); targets under 30k triangles build on the GUI thread (~30 ms);
    each sweep step rebuilds the target (a background job).
27. **Backlighter geometry (optional).** The shield's mounting sector (the sheet runs past
    9.2 mm there); backlighter beams aimed off the foil center, or foils off the equator
    (both expressible with `backlighter_foil(los)` and `FramePointing`; neither explored).
    *Data:* a shot's real backlighter dimensions.
28. **Diagnostic windows.** The expanded foil's density (held at 1 g/cm³), a backing that
    absorbs and re-emits part of what the gold transmits, a computed gap-closure time
    (the one public datum: a gap between window and wall of at most 400 µm as built,
    expected to close to 200 µm, on a shot that matched its nominal prediction; Haines
    et al., Phys. Plasmas 33, 112712 (2026), which also gives the windows as "≈ 3 µm"
    of gold against the Letter's 2.5 µm), the capsule asymmetry of an open gap of given
    width as an uncertainty term, and
    the windows' l = 4, which neither quad powers nor pointing removed. From a scratch
    study (2026-11-18):
    - The burnthrough time of a backed foil. The bare heat front (`HeatFront`, Hammer &
      Rosen) puts a 2.5 µm gold window backed by 35 µm of epoxy through at 4.8 ns on the
      ignition-scale pulse whose window-emergent flux Lester et al. (Phys. Rev. E 114,
      L023201 (2026)) measured leading at 7.5 ns (8 ns predicted): 2.7 ns early. The
      tamper holding the gold dense is the likely cause; no term models it. Calibrated to
      that edge, a fitted heat-front coefficient (×0.43) and a constant offset (+2.7 ns)
      put 2 µm of gold through 0.45 ns apart (6.6 and 7.05 ns), mid-peak either way.
      "Backed windows" (`examples/tamped_foil.jl`, `litreview/window_burnthrough.md`)
      accounts for the 2.7 ns with nothing fitted to that edge. The backing does not hold
      the gold back: a backed foil starts transmitting 0.03–0.07 ns earlier than a bare
      one in a 1-D model that reproduces Nova's bare foils to 0.04–0.23 ns. The measured
      edge is the flux emerging through gold and backing, and the backing burns through
      1.8–2.4 ns (35 µm of epoxy, at the measured ablation rate of polyimide) or 2.7–3.3 ns
      (50 µm of Be) after the gold. The prediction on a public drive is 7.0 ns (6.7–7.5
      over the uncertainty of the rate; 7.4 ns for ≈ 3 µm of gold). For the hohlraum the
      loss through the window starts when the gold burns through, about 2 ns before the
      flux emerges, so both calibrations to the edge start the loss that much too late.
      No backing option was added (the albedo does not depend on the backing). The
      docstring now carries the caution for ignition-scale drives, and
      `DiffusionTransmission` is an opt-in transmitted fraction (steady diffusion through
      a foil whose opacity falls as T^−1.5 across it; no fitted parameter; Lester's eq. 2
      at a constant opacity). It gives 0.17–0.175 for 2.5 µm, against the Letter's
      extracted 0.18 (a blind check) and the default's 0.235. On a public mock-up the
      peak loss of the windows falls from 10.4 to 7.0 TW and their peak δ₂ from 2.17 to
      1.45%. The default is unchanged (user, 2026-11-18). An intact foil whose `before` is
      a `WallLossAlbedo` now follows the wall history of its element under
      `wall_history`, since an intact window has no reason to be brighter than the wall
      (user, 2026-11-19). Before burnthrough it equals the wall around it (the δ₂ of the
      windows there was up to 0.96% on a mock-up, now zero), the transition blends from
      the albedo of the history, and the stored energy that a burned-through foil
      re-emits is scaled as its albedo is (a choice left to the user). *Remains:* the
      onset under an ignition-scale drive (the 1-D foil lags the front by 0.6 ns and
      passes 0.49–0.54 of the default's energy over 4–9 ns, 0.76–0.79 of the
      option's; but an onset that follows the 1-D foil is 0.1–0.3 ns late on Nova's
      measured times, so the Nova-calibrated `width` stays); what the backing
      absorbs against what leaves, in `drive_history`'s accounts. *Data:* emergence
      through the same gold with another backing (Be emerges 0.9–1.0 ns after epoxy),
      then another gold thickness; the Letter's window flux with its solid angle
      calibrated.
    - The burned-through windows' symmetry against Lester et al.'s estimate. On the
      study's pointing, two 2 mm windows of 2 µm gold burned through change the capsule
      drive by δ2 = 2.5% (c₂₀/c₀₀ +0.58%, c₂,±₂/c₀₀ +0.84 and −0.49%) and δ4 = 0.7%,
      against their eq. 3 estimate of a −0.9% drive P2 for 2.5 µm. Reconciled
      ("Backed windows"): eq. 3 is the magnitude of the axial P2 of the windows' flux deficit;
      about the hohlraum axis the drive's P2 is positive, and a public mock-up of that
      scale gives √5 c₂₀/c₀₀ of +0.92 to +1.04% for wall albedos of 0.80–0.88 against
      the 0.9%, with a closed form on a sphere held in `test/windows.jl` (c₂₀ to 1%).
      The m = ±2 part carries three quarters of the l = 2 power. *Remains:* the
      smoothing factor of a cylinder in eq. 3's form (a sphere's overstates it by 13%);
      the windows' l = 4. *Data:* their N250622
      against N230729 (core P2 −1% against 9% prolate; M2/M0 12% along the window
      axis) as a check, not a fit target.
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
    refinement of the radiation mesh where deposition varies fastest. *An imager along
    the view* (user review, 2026-11-13): Compute ▸ X-ray image along the view as a
    recordable diagnostic of the Diagnostics tab, so that a drive history records it per
    step and the export writes it. It needs a diagnostic kind without an optic or a
    detector response (a line of sight, a field and a sampling: the unfiltered radiance
    image), a definition made from the camera and saved with the session (the tab
    already keeps session definitions), and its place among the x-ray imagers; the
    image view's current image is exported meanwhile (`image.csv`).
32. **Scale.** Segments between the capsule's balls (~21M rays left in the 11.4k bored
    cavity's `:sampled` build, ≲ 7% of it); GPU mat-vecs (≲ 25% end to end; not worth
    CUDA.jl yet) ("H-matrix visibility certificates", "Prism certificates for the waist
    gap"); a cheaper predictor for ACA give-ups (~20% of assembly) and near-field
    quadrature reuse; an `HMatrix` has no direct solve; `storage = Float64` densifies an
    `HMatrix` without warning; `Keyframes` build a deposition operator at every time.
33. **Code health** *(survey; cleanup of 2026-11-25)*. GUI radiography is untested (it
    records nothing yet). `TabulatedCurve` interpolates as v₁ + w (v₂ − v₁), not through
    `_interp`: switching changes its values at the last bit, so it waits for a change that
    remeasures the Dante responses. Left open by the cleanup of 2026-11-25:
    - **`_drive_history`** is still one 400-line function with closures over about forty
      locals. Its per-element state is typed (`_StepState`), and the split into
      `_init_history_state`, `_solve_history_step!`, `_record_history_step!` and
      `_assemble_history_output` waits for a state struct. `test/inference.jl` (one
      compiled instance for every closure, JET `report_call`) must stay green through it.
    - **Include order.** `src/Briareus.jl` still has forward references between files
      (functions only; types and constants are defined before use). The workspace and
      JSON readers of VisRad call each other. The VisRad and diagnostics names are
      candidates for submodules like `Briareus.Benchmarks`.
    - **Lines of sight through closing lips.** The Dante footprint scans the lip
      triangles linearly per ray; a small BVH over them would remove the scan.
    - **Declined, with the reasons** (2026-11-25; the items of the cleanup plan that are
      closed without the change).
      - *A type parameter for `Scenario.lpi`.* The plan wanted `Scenario{P,M,L}` so that
        the reference-model check of a constraint could be inferred. `_ReferenceModel{L}`
        gives the check the type on the constraint side, and it passes JET `report_opt`
        with `lpi::Any`. A parameter on the scenario would compile the optimizer and the
        gradient drivers once per model type, and every `remake(s; lpi)` and every CBET
        schedule given as a closure is a new type. "Compile latency (2026-10-25)"
        measured that cost on the history march, 52 s of inference in the suite.
      - *Type parameters for `Scenario.model` and `Evaluation.scenario`.* The fields
        that an objective takes through `ev.scenario` (the target, the contracted
        responses) are concretely typed, so they are inferred through the abstract
        field. JET `report_opt` finds no dispatch in `asymmetry_objective` on Float64 or
        duals, and one call of `symmetry_objective` on a cached evaluation takes 0.17 µs
        and 112 bytes. The model is used only when the full flux is built, once per
        scenario. The parameters would compile every objective once per scenario type
        for no measurable gain.
      - *The exact ellipse distance in `ApertureKeepOut` and `SphereKeepOut`.* Measured:
        the first-order distance never overstates a clearance and falls short by at
        most 6 µm for clearances up to 0.2 mm and 29 µm up to 0.5 mm ("The ellipse
        distance of the keep-outs (2026-11-25)").
      - *Merging the line parser of the VisRad readers into `read_csv`.* `read_csv`
        scans the whole file by character, since a quoted field may hold line breaks,
        takes the first line that is not a comment as the header and keeps quoted text
        as written. The VisRad readers work by line, find the header after a preamble
        of free text, strip quoted fields and accept columns separated by blanks in the
        result tables. One function for both would hold two parsers. The digitized
        tables of data/ do share one loader (`_data_table`).
