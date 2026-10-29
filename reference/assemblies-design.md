# Assemblies from parts: design (2026-10-29)

The user's review of the GUI (ROADMAP item 24) asked whether the component framework
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

## Today

- **Target kinds.** `TargetSetup.kind` is `:hohlraum`, `:foil` or `:custom`
  (`controllers/target_setup.jl`). A hohlraum builds core `hohlraum(...)` in the chamber
  frame (wall, caps, LEHs, lips, capsule, windows); a foil target builds its foil and the
  main-target sphere. Parts (`PartSpec`) are added on top.
- **Frames.** Every part has a frame (`part_frames`), so beams aim relative to any part.
  The template has only `TEMPLATE_FRAMES` (TCC, the two LEHs, the backlighter foil); its
  wall, caps and capsule have none. This is the reported bug.
- **Names carry meaning.** The GUI finds what to analyze by component name: `:capsule`
  (37 uses: modes, keep-out, Dante, TaLIS and clearance roles), `:wall` (9), the LEHs and
  lips, `FOIL_NAME`. A part takes the capsule role by being renamed `:capsule`.
- **Core is already generic.** Core analyses take the names as keywords (`observe`,
  `wall`, `foil`, `object`, `targets`), and the VisRad import already builds a
  cylindrical hohlraum as one part with prefixed pieces (`<name>_wall`, …,
  `src/visrad_workspace.jl`). `Component.part` and `parts`/`part_pieces` already group a
  part's pieces. The radiation model is global over components, so two cavities in one
  target need nothing new.

## Model

1. **Roles, not names (GUI).** What an analysis uses is a role held by a part or a piece:
   - `observed` (one surface whose x-ray drive is analyzed): a sphere gives modes and the
     capsule keep-out; a flat face (a foil, a package's face) gives its map and spot
     statistics;
   - `foil` (any number; one active, see the active-foil change of 2026-10-29);
   - `wall` (re-emitting cavity surfaces: drive history `wall`, LEH closure, Dante);
   - `lands_on` for clearance and TaLIS (wall, observed, foils), today `_lands_on`.

   Names become free. Sessions convert: a component named `:capsule` takes `observed`.
2. **A hohlraum part kind.** `PART_KINDS` gains `:hohlraum`: core `hohlraum(...)` placed
   by the part's frame, with its pieces named `<part>_wall`, `<part>_cap_top`, … (the
   VisRad primitive's scheme) and `Component.part = <part>`. Its fields are today's
   template fields (length, radius, LEH radius, wall albedo and conversion, LEH closure),
   its windows, and per end: an LEH (radius), closed, or open to a package. A halfraum is
   a hohlraum part with one end open. Its LEH frames are child frames
   `<part>.leh_top` and `<part>.leh_bottom`, so beams can aim in either halfraum's LEH.
3. **Every piece is an aim point.** The aim-frame list shows each part and, under a
   part with several pieces, each piece ("hohlraum ▸ wall"). A piece aims in its part's
   frame, its center the piece's centroid (`part_centroid` on the piece). This fixes the
   reported bug for the wall, caps and capsule.
4. **Templates become starters.** File ▸ New target creates parts and leaves no template
   behind: *Hohlraum with capsule* (a hohlraum part and a sphere with `observed`),
   *Backlighter foil* (a rectangle with `foil`, ROADMAP frames phase 5), *Halfraum pair
   with shock tube* (two halfraum parts facing each other, a cylinder package between,
   its end face `observed`), *Double backlighter* (two foils, the main target sphere),
   *Opacity halfraum* (a halfraum, a sample foil at its open end). `TargetSetup.kind` and
   `TEMPLATE_FRAMES` other than `:tcc` go away; the Target page keeps the mesh and the
   parts list. Starter geometry is illustrative (mock-ups, public dimensions only).
5. **Beams per sub-assembly.** Nothing new in the beam model: a quad aims in a
   halfraum's LEH frame, so moving the halfraum carries its beams (phase 1 frames). The
   presets aim in the frames of the starter's parts.
6. **Analyses with several cavities.** Drive and Dante read the `wall` and `observed`
   roles. Capsule-only views (modes, keep-out, the θ–φ map) need a sphere `observed`; a
   flat `observed` face gets the foil map's face view and spot statistics. One `observed`
   surface at a time, like the active foil.

7. **Nested assemblies (a component tree).** Components are managed as a tree, as the
   beams are (groups, quads, beams):
   - An `Assembly` has a name, a frame (a `FrameDef`, so it can be placed and swept like
     a part) and children: parts or other assemblies. Parts and assemblies may both sit
     at the top level, as in VisRad.
   - Membership is a field (`PartSpec.assembly`, an assembly's `parent`), separate from
     placement: a member's placement defaults to its assembly's frame but may refer to
     any frame through its steps. Moving an assembly moves its members and the beams
     aimed in their frames.
   - Show/hide (drawing) and include/exclude (calculations) act on a subtree, with the
     tree's mixed state for a partly shown assembly, as for beam groups. Roles stay on
     parts and pieces; an excluded subtree gives up its roles until included again (the
     active foil and the `observed` surface fall back as for a removed part).
   - A part with several pieces (a hohlraum part, a CAD part split into pieces) shows
     its pieces as children, like a quad's beams; pieces are not separately movable.
   - Tree editing: new assembly, move into an assembly (keeping the members where they
     are, as `set_part_frame!` does), ungroup, rename, remove a subtree. Sessions store
     the tree; the analyses see the flattened components and never the tree.

   One workflow this serves: a CAD workspace imported as one part or assembly, proposed
   add-ons mocked up as separate assemblies, and configurations compared by including
   or excluding them. The design does not depend on that workflow: the tree is general
   (any nesting, any mix of parts and assemblies, any number at the top).

## Bit for bit

The starter *Hohlraum with capsule* calls the same core `hohlraum(...)` at the identity
frame, so its meshes, view factors and drive equal today's template's exactly; only the
names change. A test builds both and compares element tables and the drive on Float64
(names mapped). The foil starter changes mesh counts (frames phase 5): measure and state
it.

## Phases

A. **Roles in the GUI** (no visible change, bit for bit): `observed`, `wall`,
   `lands_on` read from roles; the capsule role stops renaming; sessions version up.
B. **Hohlraum part, piece aim points, starters:** the reported bug fixed; target kinds
   removed; sessions with a template convert to parts.
B2. **Assembly tree:** `Assembly`, membership, subtree show/hide and include/exclude,
   the tree's editing commands, sessions; a CAD import as one assembly of its parts.
C. **Halfraum ends and packages:** open ends, a package `observed` face, the shock-tube
   and opacity starters.
D. **Profiles:** frustraums and multi-section hohlraums (a radius per section; four
   sections for the tentless design).
E. **Double backlighter with ARC:** the starter, ARC beams as a facility beam set (public
   description only; mock-up timing), per-foil results if one active foil is not enough.

## Phase A starting point (2026-10-29)

The several-foils change (5bef323) is phase A's pattern for one role: a flag on
`PartSpec` (`foil`), a session-level choice (`active_foil`), the analyses given the name
(`flux_lookup(s; foil)`, the optimizer's and sweep's `foil`), every holder counted as a
landing surface (`foil_components`), a session watcher for role changes that do not
rebuild the target. `observed` (today `:capsule` by name) and `wall` follow it. Uses of
the fixed names (`:capsule`, `:wall`, the LEHs, caps, lips, `FOIL_NAME`) in the GUI, by
file: `target_setup.jl` 35 (mostly the template build, which phase B replaces),
`view.jl` 9, `timeline.jl` 8, `session.jl` 6, `optimizer.jl` 5, `talis.jl` 4,
`shell.jl` 4, `sweep.jl` 3, `review.jl` 3, `drive.jl` 3, `workspace.jl` 2. Core takes
these as keywords already (`observe`, `wall`, `foil`, `object`, `targets`).

## Decisions (user, 2026-10-29)

1. Piece names `<part>_wall`, VisRad's scheme (shown in the tree as "hohlraum ▸ wall").
2. `TargetSetup.kind` disappears entirely; the starters leave only parts and assemblies.
3. One `observed` surface at a time to start (as the active foil); several later if a
   case needs it.
4. Components form a tree of nested assemblies (item 7), general enough for any
   workflow of top-level parts and assemblies, not built around one.

## Risks

- The 37 `:capsule` uses span drive, Dante, TaLIS, clearance and the optimizer; phase A
  must be a pure refactor, checked by the existing GUI suite unchanged.
- Sessions from every version so far convert through one path; keep its tests.
- A halfraum's open end meeting a package needs no gap between components; the
  enclosure checks (`closed = true`) must allow two cavities joined by a package.
