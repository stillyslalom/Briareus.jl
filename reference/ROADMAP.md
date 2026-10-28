# Roadmap

Open work only. Each item says where it stands, what remains and what data would settle
it. What has been done is in the git history and in `reference/research-notes.md`, under
the sections named here.

## Release

1. **Registration.** Register Briareus in General, then BriareusGUI
   (`subdir=BriareusGUI`; its `Project.toml` must carry no `[sources]`). Risk: AutoMerge's
   load test on the GUI's precompile workload. Then `pkg> app add BriareusGUI` and
   `briareus --shortcut`.
2. **Desktop launch.**
   - The macOS app bundle and splash are untested (the macOS `qml` tool needs each Qt
     framework folder on dyld's fallback path).
   - Whether the Windows shortcut groups with its window on the taskbar (the window now
     sets its own AppUserModelID; the `.lnk` does not carry it).
   - The first frame takes 21 s on Windows: 15 s loading package images and 3 s compiling
     in QML's `__init__` (upstream; a system image would remove both). A bare `briareus`
     in a terminal shows nothing until the window.

## VisRad parity and facility practice

3. **VisRad deposition: shadowing and beam totals.** *State:* `VisRadProfile` is VisRad's
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
4. **Benchmark against VisRad.** Run `write_visrad_kit` cases at work and compare with
   `compare_visrad`: VisRad's radiosity, absolute drive levels and mode amplitudes
   (deposition alone is compared so far).
5. **VisRad import: what is skipped.** Each case is warned about or refused; most need a
   sample workspace that uses the feature.
   - Beams: unconverted light (2ω/1ω on): read VisRad's green/red powers into the ratios
     of `unconverted_deposition` (the model exists; the import still takes 3ω only);
     time-dependent powers (power mode ≠ 0; VisRad's time unit is seconds),
     energy-normalized beams, named pulse shapes; DPR beams; pointing method 1 and
     target-mount or diagnostic-port frames; objects placed by a key point other than
     the centroid.
   - Objects: the `Time-Dependent …` tables, laser reflectivity, power self-emission;
     wall thickness, partial spheres and the other primitives (cylinder, cone, disk, box,
     polygon, halfraum, rugby and slotted hohlraums, surfaces of revolution,
     multi-element objects); the element order of primitive spheres and hohlraums.
   - Workspace details to check on one small workspace with a rotated cylinder and a
     rotated CAD part: the Rotation/Spin mapping, the cylindrical position order, primitive
     grid spacing, custom profiles' spot-size definition.
6. **Confirm against facility practice.** The TaLIS 0.5° pivot per component (TCC or the
   mount); the 1ω footprint's side and shape; which quads count as open for the 20°
   reflection rule; the NIF CMT XML schema; the SST preamble and the element-export layout
   ("To verify with real files at work").

## Physics against shot data

Items 7–12 follow the 2026-10-22 literature review against Kritcher et al. 2024
("Literature review against Kritcher et al. 2024"; `reference/litreview/`).

7. **The P2 swing (N221204, PRE 109 Fig. 10).** *State:* nothing fitted to Fig. 10; the
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
8. **Shot-data calibration (CBET and effective parameters).** *State:* `n221204_cbet()`,
   `hybride_cbet()` and `hybride_cbet_iraum()` ("Joint multi-shot CBET calibration").
   *Remains:* fit effective parameters (cone multipliers, CBET multiplier and clamp, bubble
   transparency, wall speed, albedo) to measured in-flight and hot-spot P2/P4 and keyhole
   symmetry, with uncertainty; reproduce the P2-vs-Δλ sequence N210808 → N230729 and Ruof
   2026's rise cone-fraction response; the fitted Te's dip at 5.5 ns and the peak clamp
   below HYDRA's (try a drive-following Te(t) prior); plasma-state tables as inputs,
   outer-cone depletion, an inner-quad SBS proxy, the 23.5°/30° split. *Data:* N230729
   transfer curves; per-shot plasma states.
9. **Absolute energetics on N221204** (`reference/litreview/energetics.md`,
   `late_time_drive.md`). *State:* capsule share 13.8 % and LEH 24.7 % (HYDRA 16.2 /
   25.2); Dante 1 at 8 ns 290.9 / 314.2 eV over the initial / open area (measured 310.1).
   *Remains:*
   - The drive after the laser falls: 37–54 eV low at bang on three shots; the missing
     ~250–450 kJ is in the plasma, with an unpublished return rate. *Data:* HYDRA's energy
     inventory through the coast, a post-shot synthetic Dante, Dante spectra between
     laser-off and bang, gated LEH images after laser-off.
   - The capsule share (2.4–4.2 points below HYDRA). *Data:* the blow-off opacity or
     HYDRA's absorbing radius; the shell thickness at 5–7 ns; the radius histories.
   - The early peak: 13–20 eV above Dante at 5.0–5.5 ns over the open area.
   - Other sinks: the He gas (4.5 % in HYDRA), LEH windows.
   - The Zylstra pair's LEH difference is 2–4 eV too large (N191118's pointing and cone
     fraction are unpublished); the hot-spot velocity needs per-quad delivered powers.
10. **Spectral radiosity.** Left: measured Ni XRD and thin-filter transmissions, Dante 2's
    and Omega's channel hardware; a measured soft spectrum of spots in a cavity; the
    M-band's angle dependence (Dewald 2020) and its emission from the outer bubble; the
    closing LEH plasma's emission; a time-resolved M-band P2 against Dewald's Fig. 10; the
    M-band split in the `CapsuleDrive` handover.
11. **Implosion metrics and UQ.** Time-dependent error budgets and multi-shot trending.
12. **Backlighter physics.** The low-Z end of He-α CE (one Ti point; no V/Sc data at
    351 nm); pulse length and spot size as variables; continuum and higher lines in the
    filter pass band, detector response, the object's transmission.

13. **Unconverted light over time.** *State:* `unconverted_deposition` is a steady
    source, and the GUI drive adds it. 1ω and 2ω light does not enter the hohlraum and is
    neglected there. *Remains:* irradiance histories on clipped components and
    backlighter foils, from the pulse shape and an unconverted fraction in time:
    - whether a surface stays below the irradiance at which it stops being specular;
    - whether a foil forms a pre-plasma before its backlighter beams arrive (NIF User
      Guide §4.4: plasma above ~10⁸ W/cm² on metals; the 1ω prepulse passes that
      3–5 ns before the pulse).

    *Data:* the unconverted fraction vs. pulse shape; plasma and damage thresholds of
    Au, Al and Si at 1ω and 2ω.
14. **Per-color x-ray sources and a materials library.** *State:* the GUI drive weights
    each color's x-ray source as an equivalent 3ω power, exact only for a conversion
    independent of intensity. *Remains:*
    - a `solve` that takes each color's source with its own intensity-dependent
      conversion;
    - a library of material defaults (albedo, conversion vs. intensity and wavelength,
      absorption, cold reflectance), each with its data source and uncertainty, in
      place of the per-template constants.

    *Data:* the sources under "Unconverted light: absorption and conversion vs.
    wavelength", with the gaps listed there.

## Coupling and new modes

15. **Briareus → CompactLES** (`reference/litreview/hydro_coupling.md`; "Capsule drive
    hand-over and CompactLES stage 0"). *State:* the extension `CompactLESBriareusExt`
    (CompactLES branch `briareus-drive`, unmerged); stage 0 runs in 1-D, r–z and 3-D; the
    N221204 drive gives 340 km/s and stagnation 10.2 ns (bang 9.53); the late lag is the
    rocket closure (CompactLES A22). *Next:* r–θ with the a_l, r–z with AMR, a 1 µm 3-D
    shell; stage 1, gray flux-limited diffusion with a Marshak boundary, validated on
    Shestakov 1998 and then Kritcher 2024 Table III; stage 2, multigroup with the M-band.
16. **Polar direct drive** (`reference/litreview/polar_direct_drive.md`). *State:* stage 0
    (hard-sphere illumination, `mode_rms`, spherical pointing, the LMJ fixture), the corona
    operator (stage 1), the Ba24 σ₁ comparison and optimization. *Remains:* stage 1's
    absorbed power per radius, the P_abl/Ye21 proxies, `Scenario`/`evaluate` integration,
    the CompactLES snapshot example; OMEGA fixtures; open decisions 13–15 in the review.

17. **Photometric forward modeling of x-ray sources.** From a source's spectrum and
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

## Targets, CAD and the GUI

18. **Backlighter geometry (optional).** The shield's mounting sector (the sheet runs past
    9.2 mm there); backlighter beams aimed off the foil center, or foils off the equator.
    *Data:* a shot's real backlighter dimensions.
19. **Diagnostic windows.** The expanded foil's density (held at 1 g/cm³), a backing that
    absorbs and re-emits part of what the gold transmits, a computed gap-closure time, and
    the windows' l = 4, which neither quad powers nor pointing removed.
20. **CAD.**
    - A transmitting thin-film option for CAD tents and LEH windows, so the cavity can
      include the tent region.
    - Remeshed CAD: receivers through an element's outline; an adaptive `angle` for thin
      rods and fillets; closing cracks where two faces tessellate one edge differently.
    - Fits: rims on the outer half of a torus stay second order; one surface shared
      across parts; a quadrature at a cone apex; TaLIS curvature kinds for cones and tori.
    - View factors: mid-tier flat pairs by the double contour; circular plane families for
      the NIF halves; `_earclip`'s fallback on pinched outlines ("Near flat pairs by the
      double contour").
    - Laser deposition on sliver CAD is slow in `any_hit` (~5 s per beam on the full
      157k-element NIF CAD target).
    - `fit_surfaces` cost: a variable-radius fillet is quadratic (16–50k facets: 16–107 s).
21. **GUI on CAD.** Partial curved fits of CAD parts (wait for a workspace that benefits);
    the first import's commit tick (1.5 s) and first field coloring and selection (0.2 s
    each); the CAD capsule's remesh (1.0–2.0 s, now clustering); a "fit the hohlraum
    template to CAD parts" action; the foil objective with one shared parameter set for
    beams aimed at one point does not move ("GUI on CAD: import, parts and first frames").
22. **GUI requests.** Frames on a port axis; clipping volumes; STEP import (needs a CAD
    kernel); an x-ray image that follows the camera; adaptive refinement of the radiation
    mesh where deposition varies fastest.
23. **Scale.** Segments between the capsule's balls (~21M rays left in the 11.4k bored
    cavity's `:sampled` build, ≲ 7 % of it); GPU mat-vecs (≲ 25 % end to end; not worth
    CUDA.jl yet) ("H-matrix visibility certificates", "Prism certificates for the waist
    gap").
