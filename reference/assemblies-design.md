# Assemblies from parts: design (2026-10-29)

The user's review of the GUI (2026-10-29) asked whether the component framework
extends to the compound targets fielded on NIF, and reported that beams cannot aim
relative to the hohlraum template's components: the built-in hohlraum is a special
assembly. Cases named:
- two halfraums driving a shock-tube physics package between them (Rayleigh–Taylor and
  Richtmyer–Meshkov experiments);
- opacity halfraums (a sample heated by the cavity, viewed by a backlighter);
- the double-backlighter ARC target: two foils and ARC beams separated in time and
  space for two time-resolved radiographs
  (https://lmf.llnl.gov/about/keys-to-success/targets);
- from the same page: four-section (tentless) hohlraums, frustraums, planar packages.

This file records the design and the phased plan. Remove each phase from it once that
phase lands; the docstrings and research notes then hold what was built.

Phases A, B, B2, C, C2 and D, the placement triad, piece roles, the face and foil
retune, the package seams and per-cavity balancing have landed and are removed. Their
leftovers are in the roadmap ("Assemblies: phase E and leftovers"). Phase E is open.

## Today (before phase E)

- **Parts only.** A target is its parts (`PartSpec`, `controllers/part_spec.jl`).
  Starters (`starter_parts`, `insert_starter!`) insert parts in a frame: *Hohlraum with
  capsule*, *Backlighter foil*, *Halfraum pair with shock tube*, *Opacity halfraum* and
  the frustraum (`STARTERS`). Starter geometry is illustrative (mock-ups, public
  dimensions only).
- **A model tree.** Parts and assemblies (`:assembly` parts) nest to any depth through
  `PartSpec.assembly` (`controllers/assemblies.jl`); the analyses see only the parts'
  components.
- **Frames.** Every part has a frame, so beams aim relative to any part or piece, and
  moving a part carries the beams aimed in its frame.
- **Roles.** The observed, foil and wall roles are flags on parts (`PART_ROLES`). Any
  number of parts may hold the foil role; one foil is active (`active_foil`), and the
  analyses take that one. One `observed` surface is analyzed at a time.
- **Beams.** The beam table holds the facility presets of NIF and OMEGA. ARC is not
  modeled: `talis_check` lists the rule on ARC footprints among those it does not check.

## Model (phase E)

4. **Templates become starters.** The starter still to build is *Double backlighter*
   (two foils, the main target sphere).
5. **Beams per sub-assembly.** Nothing new in the beam model: a quad aims in a part's
   frame, so moving the part carries its beams. The presets aim in the frames of the
   starter's parts.

## Phases

E. **Double backlighter with ARC:** the starter, ARC beams as a facility beam set (public
   description only; mock-up timing), per-foil results if one active foil is not enough.

## Decisions (user, 2026-10-29) that bear on phase E

3. One `observed` surface at a time to start (as the active foil); several later if a
   case needs it.
5. Starters are insertable assemblies (Insert ▸ Hohlraum, Halfraum, Shock tube, Foil, …),
   placed in a selected frame (the selected part's, piece's or assembly's; TCC with
   nothing selected); File ▸ New is an empty target plus one insert. Inserted names are
   made unique (`hohlraum_2`).
9. `foil` and `observed` stay two roles. A `foil` is laser-driven, with spot uniformity,
   the foil objective and a line of sight. An `observed` surface is x-ray-driven, such as
   capsules, package faces and opacity samples. A flat `observed` face reuses the face view
   of the foil map.
