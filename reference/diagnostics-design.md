# Diagnostic definitions: design (2026-10-29)

Aim (user, 2026-10-29): synthetic diagnostics as good as public information allows,
and a framework that a user group extends with its own (non-public) details for its own
configurations, without editing the package. First kind: the soft x-ray power
diagnostics (NIF Dante 1 and 2, Sentinel, Omega Dante). Decisions 10 and 11 of
`reference/assemblies-design.md` (whole-assembly view, collimator field of view,
reference area, shields as parts) are this file's geometry. Remove each phase once it
lands; the docstrings and research notes then hold what was built.

## Today

- Definitions (phase 1, 2026-10-29): `src/diagnostics.jl`, the docs page "Diagnostic
  definitions" (schema, search path, extends, CSV curves, kinds); the kind
  `xray_power`; built-ins `nif_dante1_2021` and `nif_dante1` in `data/diagnostics`,
  from which `dante1()` and `dante1_responses()` build (bit for bit with the table they
  replaced).
- Field of view and reference area (phase 2, 2026-10-29): `DanteView(...;
  field_of_view, center, area)`, `[view]` keys `field_of_view_mm`, `center_part`,
  `center_mm`, `area_part`; the view through a hole is unchanged (bit for bit).
- Not yet: port names, an area history (Tr over A(t), Glenzer 2010), a field of view
  per channel (the multilayer channel's 12 mm; Rubery 2022 finds it covers every NIF
  hohlraum), time response, uncertainty, the recording modes other than `history`,
  session embedding. `dante2()` and `sentinel(dim)` construct views in code with Dante
  1's channels.
- XRD cathodes `:Al`, `:Cr`, `:Ni` (a Henke-yield model fitted to Rubery's curves),
  SiO₂ mirrors (Fresnel, angle and roughness scaled to Rubery Fig. 3a), the Cr/Sc
  multilayer as Rubery's digitized curve; filters by Henke's tables.
- The GUI records "Dante 1" only, and only when the target has `:leh_bottom`.

## Model

1. **One TOML file per diagnostic.** TOML: comments for provenance, already the session
   format, a standard library. A file names its schema (`schema = 1`), an `id`
   (unique, `nif_dante1`), a display `name`, a `kind`, and per value optional `source`
   and `confidence` (shown in the GUI's tooltip). Unknown keys are errors (a typo fails
   loudly, with the file and key).
2. **Kinds are registered.** `kind = "xray_power"` (the Dante family) first; a kind is a
   parser and a constructor registered in code (`register_diagnostic_kind`), so imagers
   and other kinds plug in later without a schema change.
3. **Built-ins are files.** The public definitions live in `data/diagnostics/*.toml`
   and the code constructors (`dante1()`) read them; Dante 1 reproduces today's
   channels bit for bit. Missing public values are absent, not guessed; a definition
   that lacks a value a calculation needs fails with the key to supply.
4. **Overrides by inheritance.** `extends = "nif_dante1"` takes a definition and
   replaces what the file gives: channels by `number` (a campaign's filter set), the
   field of view, the reference area. A user group keeps a directory of small files on
   top of the public ones.
5. **Measured curves by file.** Any modeled component (a cathode's sensitivity, a
   mirror's reflectivity, a filter's transmission, a whole channel's response) can be a
   two-column CSV (energy eV, value) given by a path relative to the TOML; a measured
   response replaces the model, each one stated in the result's provenance.
6. **Search path.** Built-ins, then each user directory in order, then the session's
   own; a later definition of an `id` shadows an earlier one, with a warning naming
   both files. Directories come from a user configuration file (platform config dir:
   `%APPDATA%\Briareus\config.toml`, `~/Library/Application Support/Briareus`,
   `$XDG_CONFIG_HOME/briareus`; key `diagnostics.paths`), from `BRIAREUS_DIAGNOSTICS`
   (a path list, for batch jobs and clusters), and from calls
   (`add_diagnostic_path(dir; save = true)`). The GUI's Diagnostics ▸ Folders… edits
   the saved list; a shared network directory serves a whole group.
7. **Checking.** `check_diagnostics(dir)` reads every file, resolves `extends`, builds
   each channel and reports errors and missing values by file and key; the GUI shows
   the same report.
8. **Geometry from the target.** The line of sight is a port (θ, φ) or a named port;
   the field of view a collimator (see the research below) centered on a frame (TCC by
   default, any part or piece frame); the reference area an aperture component, a part
   or piece, or a number (mm²). The view rays over the field of view, occluded by all
   included geometry (assemblies decision 10).
9. **Recording is part of the definition from the start.** A `[recording]` table names
   the mode: `history` (a time series per channel, Dante), `integrated` (an image or a
   number over the shot), `gated` (frames: times, gate width, interframe time, strips),
   `streaked` (one spatial axis along a slit × time: slit position, width and
   orientation on the target, sweep window, temporal and spatial resolution). The
   streaked mode exists from schema 1 because VISAR and SOP are of value only streaked
   (user, 2026-10-29); a kind declares which modes it supports.
10. **Sessions embed what they use.** A saved session (`.briar`) carries the definitions
    its diagnostics use, with their CSVs, as it carries CAD files (user, 2026-10-29): a
    session reproduces on any machine. On opening, an embedded definition that differs
    from the one of the same `id` on the search path is reported (file, keys) and the
    embedded one is used unless the user switches. Sessions then contain the group's
    detail wherever they are sent; the GUI's Save dialog states it when a session
    embeds a definition from outside the built-ins.
11. **GUI.** Time ▸ Diagnostics lists the loaded definitions (built-in, user, session);
   each ticked one is recorded in the drive history, its power, channel voltages and
   unfolded Tr plotted; a tooltip shows the file and each value's source.

## What the research sets (research notes, "Dante-family hardware and geometry: what is public (2026-10-29)")

- **Configurations are dated.** Dante 1's channels changed in 2009, 2012–15, 2017–18 (shot
  by shot) and October 2021; each published set becomes its own built-in file
  (`nif_dante1_2021`, `nif_dante1_N170924`, …) with `valid_from` / `shots`, and `nif_dante1`
  extends the latest. A session records which one it used.
- **Field of view per channel.** The instrument's field of view (NIF: ±30 mm, the 2004
  design) is the default and a channel may narrow it (the multilayer channel: 12 mm).
  With a 60 mm field the view is limited by occlusion and the target's extent.
- **Ports from the facility table.** The line of sight is a port name resolved in the
  facility's port table (P143-274 is θ 143.25° there; Rubery's 37° is θ 143°): adopting
  the table moves Dante 1 by 0.25°, a change measured and stated when it lands.
- **Tr's area is a setting**: `initial` (an aperture's initial area, the default and the
  convention of Rubery 2024 and the BAH paper), `projected` (a part's projected area, as
  Fournier 2010), `history` (an area per time), or a number.
- **Time response per channel**: a Gaussian FWHM (≈ 400 ps for NIF) or a CSV impulse
  response, applied to the voltage histories before the unfold.
- **Uncertainty per channel** (fraction, 1σ), for the UQ path (ROADMAP item 11).
- **Not public, so absent from the built-ins**: Dante 2's channels and solid angles,
  Omega Dante's standoff and solid angles, the fielded NIF collimators. The definitions
  say so in their comments; `check_diagnostics` names the keys to supply.
- **Checks of the field-of-view path** (never fit targets), still to run on modeled
  targets: Chen 2024's LEH / open-end pair (≈ 15 / 35 TW/sr), the BAH paper's 27 %
  occlusion of Dante 2 by a ring, May 2015's cos 37° / cos 64° ratio.

## Phases

3. Public definitions from the research notes: Dante 1 by date (2009, 2012–14, the three
   2017–18 shots), Sentinel (N170924), Omega Dante (current; no solid angles, so
   `:missing` until a user supplies them), Dante 2 (line of sight only).
4. GUI: the folders preference, the diagnostics list, several diagnostics per history.
5. Next kinds, in the user's order (2026-10-29): x-ray imagers (gated and integrated,
   from the existing `OrthoCamera` rendering through filters); backlit radiography (the
   backlighter-signal path as a definition; ARC and opacity radiographs); streaked
   VISAR and SOP. SOP needs optical emission and VISAR a moving reflector, neither of
   which the radiation model has: they need a hydro coupling (CompactLES, ROADMAP item
   15) or a prescribed shock, so the streaked mode is in the schema at once and those
   kinds come with the coupling.

   Every kind is line of sight → optic → spectral response → recorder (research notes,
   "X-ray imagers, radiography and streaked diagnostics: what is public"), so these
   are shared tables in the schema: `optic` (pinholes with positions, sizes and filters;
   a multilayer mirror; a slit and relay), `filters`, `detector` (CCD, image plate with
   Rosenberg's response model, MCP strips, hCMOS, streak tube), `recording`. A gated
   pinhole array is several cameras (each pinhole its own position and time, t₀ varying
   across the sensor, Chen 2018). Open for imagers: whether the 3–5 keV SXI and G-LEH
   images need a non-thermal spot-emission model; to be measured against published
   images, not assumed.
