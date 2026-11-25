# Diagnostic definitions: design (2026-10-29)

Aim (user, 2026-10-29): synthetic diagnostics as good as public information allows,
and a framework that a user group extends with (non-public) details for the
configurations of that group, without editing the package. First kind: the soft x-ray power
diagnostics (NIF Dante 1 and 2, Sentinel, Omega Dante). Remove each phase once it
lands; the docstrings and research notes then hold what was built.

Phases 1–4 and the first parts of phase 5 have landed and are removed: the definitions
(TOML files, registered kinds, `extends`, CSV curves, the search path,
`check_diagnostics`), the field of view and reference area, the public definitions in
`data/diagnostics`, the GUI's Diagnostics tab and sessions that embed their definitions,
the x-ray pinhole imagers in core and GUI, backlit radiography in core, named ports and a
field of view per channel. The docs page "Diagnostic definitions" gives the schema; the
research notes give the public values ("Dante-family hardware and geometry: what is
public (2026-10-29)", "Public diagnostic definitions and the time-dependent Tr area",
"The x-ray pinhole imager model and its public definitions", "Backlit radiography: the
model and its public definitions", "Named ports and per-channel fields").

## Today

- Three kinds are registered: `xray_power` (recording mode `history`), `xray_imager` and
  `radiography` (`integrated` or `gated`).
- The GUI lists radiography definitions and does not record them.
- Not yet: time response, uncertainty, the streaked recording mode, channel voltages and
  the unfold in the GUI (its histories are gray: they need a spectral history).
  `sentinel(dim)` constructs a view in code with Dante 1's channels; `dante2()` builds
  `nif_dante2` (2026-11-03) with Dante 1's channels, assumed, the definition holding none.

## Model (open parts)

9. **Recording is part of the definition.** A `[recording]` table names the mode:
   `history` (a time series per channel, Dante), `integrated` (an image or a number over
   the shot), `gated` (frames with times, gate width, interframe time and strips),
   `streaked` (one spatial axis along a slit × time: slit position, width and
   orientation on the target, sweep window, temporal and spatial resolution). The
   streaked mode is wanted because VISAR and SOP are of value only streaked (user,
   2026-10-29); a kind declares which modes it supports.
11. **GUI.** Each ticked definition is recorded in the drive history. Still to plot: its
    channel voltages and unfolded Tr.

## What the research sets, still open (research notes, "Dante-family hardware and geometry: what is public (2026-10-29)")

- **Time response per channel**: a Gaussian FWHM (≈ 400 ps for NIF) or a CSV impulse
  response, applied to the voltage histories before the unfold.
- **Uncertainty per channel** (fraction, 1σ), for the UQ path (ROADMAP item 14).
- **Not public, so absent from the built-ins**: Dante 2's channels and solid angles,
  Omega Dante's standoff and solid angles, the fielded NIF collimators. The definitions
  say so in their comments; `check_diagnostics` names the keys to supply.
- **Checks of the field-of-view path** (never fit targets), still to run on modeled
  targets: Chen 2024's LEH / open-end pair (≈ 15 / 35 TW/sr), the BAH paper's 27%
  occlusion of Dante 2 by a ring, May 2015's cos 37° / cos 64° ratio.

## Phases

5. Next kinds, in the user's order (2026-10-29): x-ray imagers (core and GUI
   done; left: the SXI soft channel, which needs the W/B₄C multilayer's reflectivity;
   per-pixel gate times across a strip, Chen 2018; the perspective of near pinholes, Izumi's 0.38 m camera); backlit radiography
   (core done; left: its GUI, ARC radiographs, which need attenuation above 40 keV and
   a bremsstrahlung source, opacity radiographs, which need a spectrometer, image plates
   with Rosenberg's response, and the streaked mode); streaked VISAR and SOP. SOP needs optical emission and VISAR a moving reflector, neither of
   which the radiation model has: they need a hydro coupling (CompactLES, ROADMAP item
   18) or a prescribed shock, so those kinds come with the coupling.

   Every kind is a chain from the line of sight through the optic and the spectral
    response to the recorder (research notes,
   "X-ray imagers, radiography and streaked diagnostics: what is public"), so these
   are shared tables in the schema: `optic` (pinholes with positions, sizes and filters;
   a multilayer mirror; a slit and relay), `filters`, `detector` (CCD, image plate with
   Rosenberg's response model, MCP strips, hCMOS, streak tube), `recording`. A gated
   pinhole array is several cameras (each pinhole with a separate position and time, t₀
   varying across the sensor, Chen 2018). Open for imagers: whether the 3–5 keV SXI and
   G-LEH images need a non-thermal spot-emission model; to be measured against published
   images, not assumed.
