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

## Today (after phase B)

- **Parts only.** A target is its parts (`PartSpec`, `controllers/target_setup.jl`); the
  hohlraum is a `:hohlraum` part (core `hohlraum(...)` placed by its frame, pieces
  `<part>_wall`, … with `Component.part = <part>`, its windows on the part). Starters
  (`starter_parts`, `insert_starter!`: *Hohlraum with capsule*, *Backlighter foil*)
  insert parts in a frame; File ▸ New target is one starter.
- **Frames.** Every part has a frame, and a hohlraum part its pieces' (`<part>.leh_top`,
  `<part>.leh_bottom`, `<part>.window_k`), so beams aim relative to the wall, caps,
  entrance holes, windows and capsule.
- **Roles, not names.** The observed, foil and wall roles are flags on parts; a hohlraum
  part's wall and caps hold the wall role.
- **Core is already generic.** Core analyses take the names as keywords (`observe`,
  `wall`, `foil`, `object`, `targets`), and the VisRad import builds a cylindrical
  hohlraum as one part with prefixed pieces (`<name>_wall`, …,
  `src/visrad_workspace.jl`). The radiation model is global over components, so two
  cavities in one target need nothing new.

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
frame, so its meshes, view factors and drive equal the former template's exactly; only
the names change (checked: see "Phase B as built").

## Phases

A. Landed (roles in the GUI, bit for bit; session version 7): commit 3a07a50.
B. Landed (the hohlraum part kind, piece aim points, target kinds removed, sessions with a
   template converted to parts, Insert ▸ starters; session version 8): commit cd46de8 and
   the commit after it.
B2. **Assembly tree:** `Assembly`, membership, subtree show/hide and include/exclude,
   the tree's editing commands, sessions; a CAD import as one assembly of its parts; the
   placement triad (decision 6).
C. **Halfraum ends and packages:** open ends, a package `observed` face, the shock-tube
   and opacity starters.
D. **Profiles:** frustraums and multi-section hohlraums (a radius per section; four
   sections for the tentless design).
C2. **Per-cavity drive and whole-assembly Dante** (decisions 8, 10, 11): Tr per hohlraum
   part, a cavity selector in Drive; `DanteView` with a collimator FOV and a reference
   area; Dante FOVs sourced first.
E. **Double backlighter with ARC:** the starter, ARC beams as a facility beam set (public
   description only; mock-up timing), per-foil results if one active foil is not enough.

## Phase A as built (2026-11-01): what departs from the design

- The role is the flag `PartSpec.observed`, but its key in `PART_ROLES` stays `:capsule`
  (the existing API and tests use it; phase A's observed surface is a sphere). Rename it
  `:observed` when flat observed faces arrive (phase C).
- A VisRad object named `capsule` imports with the observed role.

## Phase B as built (2026-11-02): what departs from the design

- One LEH radius for both ends (core `hohlraum` takes one); each end is an entrance
  hole or closed (`top_end`, `bottom_end`: a closed end's hole merges into its cap, of
  the wall's material). Per-end radii and open ends come with phase C.
- Piece frames are analytic, so they exist before meshing and follow every edit exactly:
  `<part>.leh_top` and `<part>.leh_bottom` (the template's entrance-hole axes) serve
  the entrance holes and the caps (a cap's centroid is its end's center), each window has
  `<part>.window_k` (origin on the wall), and the wall aims in the part's own frame (its
  center). A VisRad primitive's pieces have no frames of their own yet.
- The aim-frame and placement lists are flat, each part followed by its piece frames
  (`hohlraum ▸ top LEH and cap`), not a nested menu.
- The caps hold the wall role for the drive history (decision 7; converted templates
  too: the wall Tr rises by at most 0.11 eV, 0.04 %, the capsule Tr is unchanged), but
  the clearance and TaLIS checks keep them as hardware, the entrance holes' edges (the
  template's checks): a hohlraum part's caps landing the beams would hide LEH clipping.
- A sphere meshes as the template's capsule only with `capsule_mesh` set (the starter's
  capsule sets it); the observed role does not change a sphere's mesh, so existing
  sessions keep theirs.
- The model tree lists a hohlraum part's row with its pieces under it (a step toward B2's
  tree); the part has no display switch of its own (its pieces do).
- The foil starter meshes as the foil template did (its rectangle's element size is the
  capsule size at the insert, so the counts are unchanged); afterwards the foil keeps
  that size instead of following the capsule size.
- File ▸ New target ▸ Empty starts with a sphere at TCC (a target is never empty), as the
  custom template did.
- Starters insert in the selected part's or piece's frame; the capsule is placed in its
  hohlraum's frame, so moving the hohlraum carries it.
## Decisions (user, 2026-10-29)

1. Piece names `<part>_wall`, VisRad's scheme (shown in the tree as "hohlraum ▸ wall").
2. `TargetSetup.kind` disappears entirely; the starters leave only parts and assemblies.
3. One `observed` surface at a time to start (as the active foil); several later if a
   case needs it.
4. Components form a tree of nested assemblies (item 7), general enough for any
   workflow of top-level parts and assemblies, not built around one.
5. Starters are insertable assemblies (Insert ▸ Hohlraum, Halfraum, Shock tube, Foil, …),
   placed in a selected frame (the selected part's, piece's or assembly's; TCC with
   nothing selected); File ▸ New is an empty target plus one insert. Inserted names are
   made unique (`hohlraum_2`).
6. Placement in the 3D view by a triad (SolidWorks' "move with triad"): arrows translate
   along, rings rotate about, the axes of the frame chosen in "Placed in"; each drag
   commits one placement step (one undo), as the Move-along row does, so placement stays
   parametric and the beams aimed in the moved frames follow.
7. A hohlraum part gives its wall and caps the `wall` role automatically; a piece's row in
   the tree can change a piece's role (a cap `observed`, say). Roles live on parts and
   pieces.
8. Drive is per cavity: each hohlraum part reports Tr from the area-weighted flux incident
   on its own `wall` pieces; the `observed` surface keeps its own Tr beside them, and for
   a single hohlraum with a capsule the capsule's Tr stays the headline (today's value).
9. `foil` (laser-driven: spot uniformity, the foil objective, line of sight) and
   `observed` (x-ray-driven: a capsule, a package face, an opacity sample) stay two roles;
   a flat `observed` face reuses the foil map's face view.
10. Synthetic Dante views the whole assembly, not one LEH: rays across a collimator
    field of view (a disk of given diameter at a center frame, TCC by default) sum every
    visible element's radiance, occluded by any geometry. Shock-tube walls and thinned
    windows count. The primary output is power per steradian (and the channel voltages);
    Tr is quoted over a chosen reference area (an LEH, a window or a typed area; default
    the nominal LEH, Dante's convention). Each Dante's FOV diameter goes in the research
    notes with a source before it becomes a default.
11. Dante shields (dimpled shielding that hides the LEH so other regions can be measured)
    are ordinary parts: they absorb, re-emit and occlude in every calculation.

Taken without a question (reverse any):
- Frames outlive exclusion: a beam aimed in an excluded part's frame still aims there
  (an excluded add-on does not move the beams of the configuration compared against).
- The aim-frame list is a nested menu (part ▸ piece); CAD parts with many pieces list
  them the same way, collapsed.
- Old template sessions convert to top-level parts (a hohlraum part and the capsule
  sphere `observed`), no wrapping assembly, so the bit-for-bit test compares like with
  like.
- The cropped-through-a-hole `DanteView` stays as the special case of a FOV equal to the
  hole; the FOV path is checked against it on a closed hohlraum (the wall's outside is
  one-sided and dark, so the two agree to pixel convergence, not bit for bit).

## Risks

- The 37 `:capsule` uses span drive, Dante, TaLIS, clearance and the optimizer; phase A
  must be a pure refactor, checked by the existing GUI suite unchanged.
- Sessions from every version so far convert through one path; keep its tests.
- A halfraum's open end meeting a package needs no gap between components; the
  enclosure checks (`closed = true`) must allow two cavities joined by a package.
