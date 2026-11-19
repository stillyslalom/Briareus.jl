# Roadmap

Open work only. Each item says where it stands, what remains and what data would settle
it. What has been done is in the git history and in `reference/research-notes.md`, under
the sections named here. Items marked *(survey)* come from the code survey of 2026-11-02;
those under "Correctness" were read in the code, the rest are leads to confirm first.

**Priorities** (user review, 2026-11-18). First, the hohlraum physics against shot data:
the gold bubble (item 10), CBET (item 11), then the energetics and the wall and spectral
closures (items 12, 4, 13, 17). Second, re-pointing studies: the public API they need
(item 34), capture and clearance rules inside the optimization (item 35) and the
constrained runs that stall at their bounds and force staged tolerances (item 36).

## Release

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

4. **Engineering closures to validate.** The re-heated wall's radiosity clamp (≥ 0,
   albedo 0, no emission; `history.jl`, eae01d4); the conversion of overlapping colors at
   each color's own flux; the liner blend in `WallLossAlbedo`. *Data:* a radiation-
   hydrodynamics run of a wall heated, cooled and heated again (a prepulse and gap).
   Also the wall albedo of a fixed-time scenario (scratch study, 2026-11-18): a standalone
   `solve` evaluates `WallLossAlbedo` at the time since the pulse start (0.754 mean over
   the wall at 1 ns in an ignition-scale picket), while `drive_history` follows each
   element's own temperature history (0.607 at the same time). Linearized phase scenarios
   need one of the two, chosen and documented, and a helper that takes a scenario's
   albedos from a history at a given time.
37. **`drive_history` radiosity convergence** (scratch study, 2026-11-18). A quasi-static
    history of a hohlraum with two thin-foil windows (`BurnthroughAlbedo`, a
    `WallLossAlbedo` wall, `wall_history` on, `transition = 0.05`) stopped at 0.1 ns steps
    with "iterative radiosity solve (radiosity) did not converge in 400 iterations" for
    190 beams; the same run at 0.2 ns steps, and at 0.1 ns with all 192 beams, converged.
    *Remains:* a reproduction on a public mock-up; then the iteration cap, or a fallback
    to the direct solve.

## VisRad parity and facility practice

*Data on hand (2026-11-18):* many VisRad runs are saved outside the repository: hohlraum
workspaces with VisRad's computed results at several mesh sizes, the probe kits, SST
pointing exports and a TaLIS report. They are non-public. Models may be tuned to and
checked against them; the files and their contents stay out of the repository, and the
notes describe them generically. Check them before listing a run as missing.

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
7. **Benchmark against VisRad.** *State:* the saved runs hold VisRad's deposition only:
   their radiosity and view-factor sections are empty in every file ("VisRad benchmark
   beyond deposition: the saved runs"). With VisRad's deposition as the source, Briareus
   on VisRad's profile law differs after the same radiosity solve by 0.13–0.17 % rms on
   the gold, −0.05 to −0.08 % in the capsule's flux and ≤ 0.03 pp in P2, all of it the
   deposition's (the drive level follows VisRad's per-beam totals, item 6). VisRad's own
   view factors, visibility and closure are not yet compared. `compare_visrad(workspace,
   solution; elements)` compares radiosity, temperatures, drive and modes element by
   element once VisRad exports them (aa19750). *Remains:* in VisRad, compute the
   radiosities of one saved hohlraum run at ~0.6 mm and one at 0.1 mm, save, and export
   the Surface Element Properties of the gold objects and the capsule with the Power
   Summary; then run `compare_visrad`, confirm the export's index base and write a
   reader for the populated radiosity section. `write_visrad_kit` cases only for what
   the saved runs lack: primitives with curved surfaces, a closed entrance hole, a
   temperature-dependent albedo. A 0.1 mm run's `:sampled` view factors (60k elements)
   take more than 5 minutes on 16 threads; the comparison there used `:centroid`
   visibility.
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
    swing", "Inner-cone SBS time histories"). Remeasured on 874dcca without new data
    ("The swing and the bubble without new data"; `litreview/bubble_drivers.md`): the
    swing is unchanged since 2026-10-25; lips held at a thermal-band depth (0.30–0.36 mm,
    an estimate) move the crossing by ≤ 0.02 ns, the wrong way; a tip slowed through the
    trough (a 0.5 ns delay) puts the crossing at 6.58 ns but deepens the minimum by 0.04,
    and N221204's trough contrast (16× total, 12× outer) sits at Izumi 2023's 15×
    ballistic criterion: not adopted. *Remains:* a HYBRID-E bubble speed and collision
    time (none found published); a thermal-band LEH radius or the warm lip's opacity; the
    rise (0.063 above Fig. 10 at 4.5 ns, 0.075 below at 5.0 ns: the model falls into the
    minimum later and more steeply) and the foot at 3–4 ns (0.015–0.054 above). The
    benchmark keeps its hand-set picket ratio, and re-pointing studies use the link
    relative to their reference (user, 2026-11-18). The link (46d88e6,
    §5 of the same section; opt-in): `picket_fluence`, `picket_ratio` and `remake(bubble;
    picket_ratio)` take the bubble's speed from the outer cones' energy-weighted picket
    fluence on the wall. On N221204 a ±0.5 mm re-pointing of the outer cones along z
    changes that fluence by ≤ 1 % and the crossing through the bubble by ≤ 0.013 ns;
    unsplit outer quads raise it 19 % (−0.19 ns). The benchmark's pointing gives a ratio
    of 1.147 against the calibration platform (hand-set: 1.1), on a proxy pointing and a
    picket energy read by eye; ratios between pointings converge to < 1 %. A published
    bubble speed at fixed picket energy and a different spot overlap would test it.
    Validation left: Young's dynamic albedo (needs N170827's delivered per-beam powers);
    Ralph 2018's 18.5 % cone fraction (needs the high-foot pointing and wall radius); the
    N210808 proxy; Callahan 2018's 39-shot curve. *Data:* bubble-tip or inner-spot
    histories (thin-wall or G-LEH), a HYBRID-E LEH-radius history, the 2DConA in-flight
    shell P2, the HYDRA multipliers behind Fig. 10, N221204's shell thickness at 5–7 ns and
    per-cone backscatter history; a tip position during the trough.
11. **Shot-data calibration (CBET and effective parameters).** *State:* `n221204_cbet()`,
    `hybride_cbet()` and `hybride_cbet_iraum()`, on one builder since d72ed90 ("Joint
    multi-shot CBET calibration"; "CBET without new data: the Δλ slope and the Te
    history"). The models' dCF/dΔλ, 2–3× steeper than HYDRA's cross-shot response, is
    not shown to be a model defect: on the one step with the pulse held (N220919 →
    N221204) the models give a hot-spot P2 step of +10.7 to +13.6 µm, the measurement
    +11.1 µm and HYDRA's runs +3.8 µm, and no model term (the clamp, a local depletion,
    the plasma state, the outside zone) reaches HYDRA's slope. A Te(t) that follows the
    Dante drive removes the fitted Te's dip at 5.5 ns and lifts the peak clamp into
    HYDRA's range, but doubles the loss and worsens the blind Fig. 10 check on three
    shots: not adopted. *Remains:* fit effective parameters (cone multipliers, CBET
    multiplier and clamp, bubble transparency, wall speed, albedo) to measured in-flight
    and hot-spot P2/P4 and keyhole symmetry, with uncertainty, through a flux-P2 to
    hot-spot-P2 transfer (HYDRA's flat cross-shot dCF/dΔλ is no longer a target: user,
    2026-11-18); Ruof 2026's asymmetric rise
    cone-fraction response (the reduced CBET amplifies cuts and raises alike, 1.2–1.5×);
    Te at the 5.5 keV fit bound at 8 ns; a laser-power Te prior (not run); plasma-state
    tables as inputs, an inner-quad SBS proxy, the 23.5°/30° split. *Data:* N230729
    transfer curves; per-shot plasma states.
12. **Absolute energetics on N221204** (`reference/litreview/energetics.md`,
    `late_time_drive.md`; "Absolute energetics without new data: the early peak and the
    missing sinks"). *State* (coarse mesh, 90e28c4): capsule share 13.9 % and LEH 24.2 %
    (HYDRA 16.2 / 25.2); Dante 1 at 8 ns 289.7 / 312.9 eV over the initial / open area
    (measured 310.1); the whole-field view equals the view through the hole on N221204
    within 0.1 eV. *Remains:*
    - The rise through the peak: Dante rises 72 eV from 5 to 8 ns, the model 47 eV over
      the initial area and 57 eV over the open area. The "early peak" (13–20 eV above
      Dante at 5.0–5.5 ns) is the open-area face of this one error; on Dante's own
      convention the model is within 6 eV there and 0.07 ns early. No existing input
      within its documented range adds more than 4 eV to the rise. *Data:* a
      thermal-band LEH radius, or HYDRA's loss history through the peak.
    - The drive after the laser falls: 30 / 38 / 47 eV low at 8.5 / 9.0 / 9.35 ns (37–54
      eV at bang on three shots); the missing ~250–450 kJ is in the plasma, with an
      unpublished return rate. *Data:* HYDRA's energy inventory through the coast, a
      post-shot synthetic Dante, Dante spectra between laser-off and bang, gated LEH
      images after laser-off.
    - The capsule share (2.3–4.2 points below HYDRA): it and the LEH share move against
      each other under every input except a lower wall loss, which the DU calibration
      rules out. *Data:* the blow-off opacity or HYDRA's absorbing radius; the shell
      thickness at 5–7 ns; the radius histories.
    - Sinks, sized from published numbers and not modeled: the He gas (4.5 % in HYDRA:
      80 kJ, 3–4 eV of drive at every time); LEH windows, diffusely scattered laser
      light (5–9 % of the light not absorbed at the first hit; one bounce sends 72 % of
      it back to the wall) and specular glint, each under 2 eV.
    - The Zylstra pair's LEH difference is 2–4 eV too large (N191118's pointing and cone
      fraction are unpublished); the hot-spot velocity needs per-quad delivered powers.
13. **Spectral radiosity and diagnostics** (`reference/diagnostics-design.md`).
    - Dante: measured Ni XRD and thin-filter transmissions; the SiO₂ mirror's angle scale,
      roughness and surface layer (assumed or fitted; extrapolated below 2.2 keV); the
      within-group spectral shape (finer groups); Dante 2's hardware and fields and
      Omega's (Dante 2 and Omega's Dante look through the hole: fields not public;
      `dante2()` assumes Dante 1's channels); time
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
    held; "Retuning faces and foils; the held level"): peak-to-mean is reported, not
    optimized; the retune holds the reference's albedos (linearized) and retunes powers,
    not pointing.
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
    - materials: `material_library()` holds the defaults that have a source, each
      quantity with its source and uncertainty (154a4c5; "Material defaults and their
      sources"). Left: the defaults without one (user): the GUI's hohlraum 0.85/0.8
      (15 eV above the core's `gold()` on the default hohlraum; Young 2020's static
      NIF setting is 0.8/0.9); `gold()`'s conversion 0.75, a 1 ns value below every
      published NIF value (0.85–0.95); the foil conversion 0.05 (core) and 0.02 (GUI),
      10–48× below the published values for gold; the capsule conversion 0.3. No
      public 3ω broadband conversion was found for the ablators (CH, HDC, Be) or the
      mid-Z foils (Ti–Ge, Ag), nor albedos specific to CH and Be.

## Re-pointing studies

Items 34–36 come from the user's review of 2026-11-18, after a study that drops a
bundle and re-points the remaining quads to hold a reference drive shape (the public
analog is the dropped-bundle tutorial, `docs/literate/dropped_bundle.jl`).

34. **Re-pointing studies on the public API.** *State:* `subset`, `GroupPointing`
    (`relative = true` for a rigid shift per quad), `asymmetry_objective(; reference)`,
    `ApertureKeepOut`, `mode_uncertainty`; since 0594191 also held parameters and power
    scales per group (`HeldParameters`, `PowerScales`), `remake(::Scenario)`, NIF beams on
    VisRad's plates in the quad frame (`visrad_nif_spot`, `spot_frame = :quad`), quads to
    beams and back (`nif_beams(quads)`, `nif_quads(beams)`, an SST table through
    `read_visrad_sst(; template)`) and `pointing_summary` ("Re-pointing studies on the
    public API"); the gated objective and the penalty terms are items 36 and 35's
    (`drive_objective` with `mode_gate`; `keepout_penalty`, `light_ceiling`). *Remains:*
    - Geometry queries for overlap and reach estimates, to settle (user): `port_basis`
      and `visible` are documented and stable (export them); `closest_hit` takes the
      internal `BVH` (a method on the target instead); the beam's section at a plane
      (`_beam_at_plane`) is shaped for CBET (a `beam_section(beam, point, normal)`);
      `BeamFrame`, `_vr_intensity` and `_vr_envelope` stay internal (`intensity` is the
      public path; an envelope method on the profile would cover the last).
    - `GroupPointing`'s default `unit = :beam` against the quads-as-units rule (item 23);
      a default `subset` and `parameters(...)` named values for the composed
      parameterizations (`subset`'s default regroups by `:group`, which renames
      parameters).
    - One pointing over several phases (scratch study, 2026-11-18). A study holding the
      drive shape at six times (picket, trough, second step, rise, and the peak before
      and after its windows burn through) with one set of aims and a power scale per
      quad and phase was built from workarounds: a copy of each beam per phase under the
      same name (so `GroupPointing(by = :beam)` gives the copies one aim); a
      `PowerScales` `by` closure that counts its calls (`by` receives no index); a main
      scenario without radiation that only carries the parameters; and a cache keyed on
      the `Evaluation`'s internal flux slot, so that the objective and each constraint
      share one evaluation per phase (thread-safe, as gradient chunks run
      concurrently). A phased scenario would replace them: shared aims, per-phase powers
      and radiation models, an evaluation per phase, and per-phase constraints and boxes
      (a laser-limit cap on each beam in each phase).
    - A pulse-shape budget (same study): each distinct pulse shape needs facility
      validation, so the per-quad multipliers were grouped afterwards into a few shapes
      (k-means on the log multipliers, the split between the inner and outer cones
      chosen by fit) and re-optimized with tied `PowerScales`. A helper for the grouping,
      and for writing the shapes out as `PulseShape` tables.
    - A ready-made smallest-shift objective (mean squared aim shift from the starting
      beams), for the minimum re-pointing that meets given gates: a per-LEH study used
      it to cost each entrance-hole size against a fixed reference drive.

    The tutorial's section "A re-pointing study, end to end" walks these on exported
    names only (a96792f): quads held, a power scale per cone, the gates, the ceiling, a
    hardware keep-out and the light on a part as constraints, and a per-beam polish as a
    second call on the first result's beams (`remake(scenario; parameterization)`; a
    schedule type saved no line over the two calls and hid the check between them).
35. **Capture and clearance rules in the optimization.** *State:* `clearance` and
    `talis_check` judge a finished design on the worst case over pointing, positioning
    and focus errors (reports, not differentiable), and `placement_sweep` runs the
    capture and clearance checks at each step. For the optimizer, since 2c1825c ("Capture
    and clearance rules as optimization quantities"): the keep-out terms as public
    quantities (`keepout_clearance`, `keepout_penalty`, `direct_fraction`); the light on
    named parts for floors and ceilings (`power_fraction` with `clipped_fraction(ev,
    names...)`); `HardwareKeepOut`, the exact worst-case distance from a beam's envelope
    to named hardware, smooth in the beams and within 1e-6 mm of `clearance` (lower by up
    to 0.4 mm where `clearance` skips points hidden from the beam, as `talis_check` does
    not); and the TaLIS rules that carry a margin, each within 1e-7 mm of `talis_check`
    on a mock-up (6e-6 mm on deep violations, below): the 3ω clearance (`HardwareKeepOut`
    with `tolerance`, the 0.5° term inside the minimization), 3ω capture
    (`CaptureKeepIn`), incidence (`IncidenceRule`) and the 1ω clearance
    (`Omega1KeepOut`). Each has a `keepout_clearance` method, so each serves as a penalty
    or as a constraint (item 36). *Remains:*
    - a loop that runs `clearance` and `talis_check` on the result and tightens the
      minimums until the review passes (item 21);
    - capture fractions per beam (today's are totals over the beams);
    - the 1ω clearance by an exact search (branch and bound now: ~3.5 ms per beam, three
      depositions' worth for 48 quads); an incidence with a gradient on curved parts
      (the angle is constant per facet, so the gradient is zero away from facet edges);
    - 1ω management, the diagnostics-to-1ω rule, the TAS stay-in rule and the LEH review
      tier have no quantity;
    - `keepout_penalty` skips every non-finite clearance, so a −Inf (a beam that misses
      every capture part with no free edge) is dropped silently: skip +Inf only (user);
    - `_ellipse_distance` (`clearance.jl`, also behind `talis_check`) loses digits for
      points inside the ellipse within ~1e-9 mm of its major axis: deep violations
      (below about −0.75 mm) read up to 6e-5 mm too low. Treating points within 1e-6 of
      the semi-minor axis as on the axis changes the distance at second order only; it
      changes report values (user).

    *Check:* each surrogate against `clearance` or `talis_check` on the worst case, as
    the tutorial does for the aperture keep-out (agreement to a few micrometers); a
    re-pointed design that passes `talis_check` without repair by hand.
36. **Constrained re-pointing stalls at its bounds.** The study's objective matches the
    capsule modes of a reference drive while holding the peak laser irradiance over the
    wall's elements at or below the reference pointing's (high irradiance drives the
    gold bubble, item 10); started on the final tolerances, its penalty runs ended with
    parameters at their bounds and the objective unmet, and stages set by hand were the
    way around it. *State* (7045105, "Constrained re-pointing: an irradiance ceiling and
    mode gates"): `optimize_pointing(...; constraints)` solves inequality `Constraint`s
    by an augmented Lagrangian around the same L-BFGS-B, and the result reports each
    constraint against its limit and what stopped the run (`PointingResult.constraints`,
    `outcome`); `irradiance_ceiling` (an area-weighted 64-norm as the smooth stand-in,
    `peak_irradiance` as the reported quantity, a `weight` per beam for outer-cone or
    picket-fluence variants), `mode_gate` and `mode_departure`, `drive_objective`. On the
    dropped-bundle tutorial's case one call meets the four gates under the ceiling from
    three starts (488–765 evaluations before 320335b, which made constrained runs
    converge at the final inner tolerance only: 1.5–6 times the evaluations since), at
    0.4 % more drive than the penalty form reaches, staged by hand or not. On the study
    itself (scratch runs, 2026-11-18) one call from the as-shot start meets every gate
    and the ceiling inside the box, every constraint active, at the drive of the earlier
    design without a ceiling; the study's penalty form, run the same way, stalls with
    the peak 7 % over the limit and a parameter at its bound (the reward for symmetry
    inside the gates balances the ceiling's penalty), so the stall was the penalty form,
    not infeasibility. A per-beam polish (`parameters(p, beams)` hands a quad result to
    per-beam shifts) moves beams by micrometers and the drive not at all. Since bf79dc2
    a run continues from a result with its multipliers (`from`), stops at a `time_limit`
    with its best design, reports the drive ratio and departures in `pointing_summary`,
    and a ceiling refuses a reference evaluated with another laser–plasma model.
    *Remains:*
    - the study's CBET and window variants, its other remedies, and a ray-traced
      check of the designs;
    - `irradiance_ceiling`'s default `exact = false` (the 64-norm alone left the peak
      0.4–0.5 % over the limit on the study and up to 3 % on a 0.45 mm mesh, depending on
      the path the run takes; `exact = true` holds the peak itself, at 0.9997 of the
      reference's on the study in one call, and lands a little below the limit with one
      or two gates up to 3 % inside theirs) (user);
    - `mode_gate` takes coefficient vectors and cannot check that its reference was
      evaluated with the scenario's laser–plasma model, as `irradiance_ceiling` now does:
      a method on an `Evaluation` that records it;
    - bound sensitivities (`PointingResult.bound_sensitivity`) match the re-optimized
      objective on a closed-form problem but overstate the box's cost about fiftyfold
      on the study, where the gates leave a flat family of designs and the multipliers
      are not unique: a direct estimate (a restart with the bound moved, through
      `from`) could be built in;
    - which quantity the ceiling should hold: the published bubble scalings use the
      outer cones' picket fluence on the wall, not the peak irradiance, and on N221204
      a z re-pointing of the outer cones changes that fluence by ≤ 1 % while spot
      overlap changes it by up to 19 % (item 10; `picket_fluence`; the peak ceiling is
      an engineering limit) (user). The scalings' quantity is not `picket_fluence`'s
      energy-weighted ∫F²dA/∫FdA but the outer picket energy, net of backscatter, over a
      fixed spot area (Callahan et al., Phys. Plasmas 25, 056305 (2018); Ralph et al.,
      Phys. Plasmas 25, 082701 (2018); Izumi et al., Phys. Plasmas 28, 022706 (2021):
      "spatially averaged outer picket fluence over the first 1 ns"). The two part on
      uneven power scales: on a dropped-bundle study, holding ∫F²dA/∫FdA at the
      reference's cost 16 % of the picket drive, while a design at the reference's
      published quantity reads 1.10× on it (scratch runs, 2026-11-18). *Remains:* the
      published quantity as a quantity and a constraint, with the trough and foot as
      well as the picket (Izumi et al., Phys. Plasmas 30, 032706 (2023): a low trough
      stalls the bubble); the outer fluence near the inner beams' centerlines (the same
      paper: bubble spikes in line with the 23.5° beams cost P2), which the study checks
      as the mean fluence on wall elements within 1.3 mm of an inner centerline; the
      scalings' exponents in the research notes (core P2 ∝ fluence^0.5, bubble position
      ∝ fluence^0.2), with the −4 to −7 µm of P2 per +10 % of fluence marked as read
      from their figures;
    - the inner budget: `iterations` caps each inner run, and set high it leaves the
      multipliers unchanged for most of a run. On a 276-parameter multi-phase study,
      `iterations = 1000` completed one outer iteration in an hour (the peak 18 % over
      its ceiling at the time limit); `iterations = 100` completed six in 40 min and met
      every limit but one, within 1 %. A default tied to the outer loop, an inner budget
      that scales with the parameters, or a documented rule (user review);
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
    - The design docs keep built phases against their own rule (remove each phase once it
      lands; `assemblies-design.md` "Today", `frames-design.md` phases 1–3 and phase 5,
      done by the foil starter).
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
      Left: in the core `GroupPointing`'s default stays `unit = :beam` (the window opts
      in; `hybride_mock_pointing` uses the quad's azimuth already), and `with_pulses`,
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
    - The documentation's screenshots (`docs/src/assets/gui`, `dev/gui_screenshots.jl`)
      predate the Diagnostics tab, the Help window and the menus of 2026-11-13.
25. **GUI density and layout (user review, 2026-10-29).** Aim toward Office-ribbon
    density, not all the way; the 3D view first. *Left:* the next pass on what still takes
    room once the user has used it. Principle: frequent commands near at hand,
    explanations in tooltips and Help, the 3D view and its inspection given the most room.
26. **GUI responsiveness** *(survey; measured 2026-11-03 through the window's shell
    offscreen, on the default design, 5.2k elements, and a CAD mock-up of 40 solids,
    10.8k elements)*. *Done:* `record!` and `can_undo` answer from a count of the
    session's observable changes between edits (before, a snapshot took 0.3–0.6 ms and
    its comparison 0.25–0.5 ms at every event, and twice per refresh with no undo step: a
    refresh took 1.8–2.0 ms then, 0.4–0.9 ms now). An undo step whose target section is
    the design's keeps the target, and restoring the view no longer makes a sector cut a
    plane cut on the way (two redraws of the 3D view each time); on the mock-up a view
    undo took 88 ms, now 6 ms, a beam undo 126 → 25 ms, a part undo 188 → 82 ms, each
    with its refresh. *Left:* the model tree and every panel are rebuilt each refresh (the
    tree 0.1–0.3 ms, the panels together under 0.3 ms: not worth caching yet); a
    `restore!` of an unchanged CAD design (a session reopened) still takes 13 ms in
    `_try_build!` (each part's cache key reads its file's time stamp) and 17 ms in
    `set_beams!` (the beams redrawn); a part undo rebuilds the target (CAD designs depend
    on the 6-entry target cache); targets under 30k triangles build on the GUI thread
    (~30 ms); each sweep step rebuilds the target (a background job).
27. **Backlighter geometry (optional).** The shield's mounting sector (the sheet runs past
    9.2 mm there); backlighter beams aimed off the foil center, or foils off the equator
    (both expressible with `backlighter_foil(los)` and `FramePointing`; neither explored).
    *Data:* a shot's real backlighter dimensions.
28. **Diagnostic windows.** The expanded foil's density (held at 1 g/cm³), a backing that
    absorbs and re-emits part of what the gold transmits, a computed gap-closure time, and
    the windows' l = 4, which neither quad powers nor pointing removed. From a scratch
    study (2026-11-18):
    - The burnthrough time of a backed foil. The bare heat front (`HeatFront`, Hammer &
      Rosen) puts a 2.5 µm gold window backed by 35 µm of epoxy through at 4.8 ns on the
      ignition-scale pulse whose window-emergent flux Lester et al. (Phys. Rev. E 114,
      L023201 (2026)) measured leading at 7.5 ns (8 ns predicted): 2.7 ns early. The
      tamper holding the gold dense is the likely cause; no term models it. Calibrated to
      that edge, a fitted heat-front coefficient (×0.43) and a constant offset (+2.7 ns)
      put 2 µm of gold through 0.45 ns apart (6.6 and 7.05 ns), mid-peak either way.
      *Remains:* a tamped-foil burnthrough model, with the backing's areal mass and
      material (epoxy, Be). *Data:* Lester et al.'s leading edge is the one point; a
      window of another thickness or backing would separate the two calibrations.
    - The burned-through windows' symmetry against Lester et al.'s estimate. On the
      study's pointing, two 2 mm windows of 2 µm gold burned through change the capsule
      drive by δ2 = 2.5 % (c₂₀/c₀₀ +0.58 %, c₂,±₂/c₀₀ +0.84 and −0.49 %) and δ4 = 0.7 %,
      against their eq. 3 estimate of a −0.9 % drive P2 for 2.5 µm. *Remains:* reconcile
      the m = 0 part with eq. 3 (sign convention and amplitude). *Data:* their N250622
      against N230729 (core P2 −1 % against 9 % prolate; M2/M0 12 % along the window
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
    cavity's `:sampled` build, ≲ 7 % of it); GPU mat-vecs (≲ 25 % end to end; not worth
    CUDA.jl yet) ("H-matrix visibility certificates", "Prism certificates for the waist
    gap"); a cheaper predictor for ACA give-ups (~20 % of assembly) and near-field
    quadrature reuse; an `HMatrix` has no direct solve; `storage = Float64` densifies an
    `HMatrix` without warning; `Keyframes` build a deposition operator at every time.
33. **Code health** *(survey)*. GUI radiography is untested (it records nothing yet).
    `TabulatedCurve` interpolates as v₁ + w (v₂ − v₁), not through `_interp`: switching
    changes its values at the last bit, so it waits for a change that remeasures the
    Dante responses.
