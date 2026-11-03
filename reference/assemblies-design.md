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

## Today (after phase C)

- **Parts only.** A target is its parts (`PartSpec`, `controllers/target_setup.jl`); the
  hohlraum is a `:hohlraum` part (core `hohlraum(...)` placed by its frame, each end an
  entrance hole, closed or open (a halfraum), pieces
  `<part>_wall`, … with `Component.part = <part>`, its windows on the part). Starters
  (`starter_parts`, `insert_starter!`: *Hohlraum with capsule*, *Backlighter foil*,
  *Halfraum pair with shock tube*, *Opacity halfraum*)
  insert parts in a frame; File ▸ New target is one starter.
- **A model tree.** Parts and assemblies (`:assembly` parts) nest to any depth through
  `PartSpec.assembly` (`controllers/assemblies.jl`); the analyses see only the parts'
  components.
- **Frames.** Every part has a frame, and a hohlraum part its pieces' (`<part>.leh_top`,
  `<part>.leh_bottom`, `<part>.window_k`), so beams aim relative to the wall, caps,
  entrance holes, windows and capsule.
- **Roles, not names.** The observed (a sphere or a flat face), foil and wall roles are
  flags on parts; a hohlraum
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
B2. Landed (the assembly tree: assemblies, membership, subtree show/hide and
   include/exclude, the tree's editing commands, sessions, a CAD or VisRad import as one
   assembly; session version 10). The placement triad (decision 6) landed after it, with
   no session version change: see "Placement triad as built".
C. Landed (halfraum ends and packages: per-end entrance holes and open ends in core
   `hohlraum`, the role key `:observed`, flat observed faces, the shock-tube and opacity
   starters; no session version change): see "Phase C as built".
D. Landed (profiles: frustraums and multi-section hohlraums, core `hohlraum(; profile)`
   and the GUI hohlraum part's *Profile*, the frustraum starter; no session version
   change): see "Phase D as built".
C2. Landed (per-cavity wall Tr in core `drive_history(; cavities)` and `wall_temperature`,
   the Drive panel's cavity selector and the history's per-cavity series; the built-in
   `nif_dante1` views the whole target; no session version change): see "Phase C2 as
   built".
Piece roles (decision 7's second half: a piece's own role, a cap observed; core
   `drive_history(; per_component)`; session version 11): see "Piece roles as built".
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

## Phase B2 as built (2026-11-02): what departs from the design

- An assembly is a part of kind `:assembly` (a `PartSpec`, like a user frame: no
  geometry), not a separate `Assembly` type. Its frame is the part frame
  (`part_frame_def`: placement fields and steps), so placement, Move along, steps,
  renames, sweeps and undo serve it unchanged. `PartSpec.assembly` is the membership (an
  assembly's is its parent); `TOP_LEVEL` (`Symbol("")`) is the top. `check_tree` (called
  by `part_frames`) refuses a member of a missing assembly, an assembly inside itself and
  a user frame in an assembly (frames stay in the Frames folder).
- Show/hide and include/exclude act on the members, as for beam groups: an assembly
  holds no state of its own, and its box and its line ("partly in calculations") are its
  components' together. Excluding an assembly puts its parts' names in `excluded`, so the
  roles fall back through the existing paths. Putting a part back (alone or with its
  assembly) now also puts back its pieces left out on their own: a part comes back whole.
- A hohlraum part's row now has a box too (its pieces' state, mixed included), changing
  phase B's "no display switch of its own".
- New assembly: placed at its parent's origin with its axes (zero offsets in
  `assembly_frame(parent)`), so grouping moves nothing; the parent is the one its members
  share, else the top. Moving into an assembly re-places a member only when it is placed
  in its old assembly's frame (TCC at the top), keeping it where it is (`_kept_in`,
  factored out of `set_part_frame!`); a member placed elsewhere (a capsule in its
  hohlraum) keeps its placement. Ungroup re-places everything placed in the assembly's
  frame in its parent's, and is refused while a placement step refers to it; beams
  aimed in it keep their aims in `:tcc`. Removing an assembly removes its subtree, refused
  while a part left behind is placed in it.
- Tree order: an assembly is listed where its first member is among the parts; a new
  assembly is appended to the parts, so no part's index changes.
- A VisRad workspace imports as one assembly named after the file (appended last, its
  objects placed in it, open in the tree); so does a multi-body STL file. Insert with an
  assembly selected puts the starter in it (window only).
- The tree's expanded nodes are not saved (as before): assemblies load collapsed.
- Bit for bit: grouping the starter leaves its frames, its components (the same objects)
  and its target unchanged, and the drive equal (test). The default design's drive is
  equal before and after this phase built with the same precompile settings; built with
  and without the GUI's precompile workload it differs at round-off (P0/P4 in the last
  digits), unrelated to this phase.
- Not built in this phase: the placement triad (decision 6); see "Placement triad as
  built".

## Placement triad as built (2026-11-02)

- One placement model: an arrow drag commits `nudge_part!` in the part's *Placed in*
  frame (the Move-along row's edit, frame values equal exactly); a ring drag commits the
  new `turn_part!`, its rotational counterpart, which re-expresses the part's angles
  (`theta`, `phi`, `rotation`) in that frame and keeps its origin. No placement step is
  appended per drag (the steps list would grow with every turn); a part without steps
  ends where a `Rotate` step about the same chamber axis puts it, to round-off. A turn
  its steps undo is refused as a fixed move is. `nudged_parts`/`turned_parts` give the
  edit's parts without making it (the preview).
- Shown for one selected part, assembly, user frame, or a piece, window or component
  (its part's triad); hidden with nothing, several things, a beam or a port frame
  selected. Drawn at the part's origin, 70 px long on screen, thin lines over the
  geometry; rings 0.62 of the arrow length.
- Snapping: an arrow by whole steps of the Move-along row's *by* field, a ring by 15°
  (`TRIAD_TURN_STEP`, the arrow keys' orbit step); Ctrl at the press snaps to a tenth
  (`TRIAD_FINE`). The amount is the decimal it prints as (`triad_snap`: three 0.1 mm
  steps are 0.3 mm).
- An arrow's amount is where the pointer's ray passes nearest the axis; a ring's is the
  pointer's screen angle swept about the projected origin (signed by whether the axis
  faces the viewer). An arrow along the line of sight (under 12 px) and a ring within 9°
  of edge-on are not grabbed.
- The preview edits nothing (the undo history records only the release): the moved
  parts' components drawn faint, the beams aimed in moved frames from their shifted
  start to their new aim, and the amount by the triad. A press released within the click
  slop is a click, so the beams and parts behind a handle stay selectable; focus lost
  mid-drag cancels it.
- Not built: hover highlighting of a handle; a planar (two-axis) drag; dragging several
  selected parts at once (the Move-along row moves them all).

## Phase C as built (2026-11-02): what departs from the design

- Core: `hohlraum(; top_end, bottom_end)`, each `nothing` (an entrance hole of
  `leh_radius`), a radius, `:closed` or `:open`. The closed end moved from the GUI into
  core unchanged (the cap and the hole's and lip's meshes merged). An open end has no cap
  and no hole; the target is then not declared `closed`, and an open end refuses corner
  fillets and sombreros (a fillet at an open end would be a lip; a sombrero shields a
  hole). With a closing lip, an end of its own radius scales `leh_open_radius` by the
  ratio of the radii.
- Bit for bit: the components, element centroids and areas of nine `hohlraum`
  configurations, and the view factors, deposition, B, H and moments of a coarse NIF
  hohlraum under `nif_quads`, serialized before and after: identical. The GUI hohlraum
  part's components (default, closed ends with and without lips, a placed part, windows)
  likewise. Validation: the coaxial-disk closed form and an aperture-closed energy balance
  (`test/halfraum.jl`; research notes, "Halfraum ends and starters").
- The GUI hohlraum part's ends take 2 (open) besides 0 and 1 (`HOHLRAUM_ENDS`), and
  `bottom_leh_radius` (0: `leh_radius`) gives the bottom hole its own radius; the top's
  is `leh_radius`. Sessions store values by key and fill missing keys with defaults, and
  the role is the flag `observed`, so no session version change and no conversion: older
  sessions load unchanged.
- The role key is `:observed` (`PART_ROLES`), labeled *Observed*. A flat observed face is
  a rectangle or disk (or a VisRad rectangle) (`observed_flat`); the drive's scenario
  tests flatness from the observed elements' normals (`_observed_flat`). A flat face gets
  a face map of its incident x-ray flux (`flux_lookup` kind `:face`, in-plane axes `u`,
  `v`), `DriveResult.face` (mean, peak, lowest, RMS nonuniformity over the whole face,
  not the foil's spot radius), `mean_flux` its area mean, and no modes; the symmetry
  objective and the capsule keep-out refuse it.
- Starters `:shock_tube` and `:opacity` (`starter_parts`): plain parts, not wrapped in an
  assembly, as the earlier starters. The package's end disks have the halfraums' radius
  and close their open ends; the shock-tube pair is two cavities, so `_closed_enclosure`
  leaves it unbalanced (as any target with more than one hohlraum part), and a single
  hohlraum part with an open end too. The beam presets did not aim in the starters'
  frames at first; see "Item 5 as built".
- The drive history of a flat observed face (built after phase C): `DriveHistory.face`
  holds the face's mean, peak, lowest and RMS nonuniformity per step (`_face_stats` on
  each step's solution, as the steady drive), `Tr` comes from the area-mean flux (the
  core history's `Tr` is from P0, within 0.2 % on the opacity starter), `legendre` has
  no rows, and the drive chart's right axis shows the RMS and the peak and lowest flux
  relative to the mean. The retune first restored capsule modes only and refused a flat
  face; since 2026-11-02 it restores a flat face's flux over its elements (core
  `FaceFlux`) and, without a radiation model, the active foil's spot irradiance (core
  `FoilIrradiance`); see "Face and foil retune as built". Sphere-observed histories are unchanged:
  the default session's history fields and drive chart serialized before and after are
  equal.
- A piece's role (decision 7, a cap observed): built, see "Piece roles as built".

## Face and foil retune as built (2026-11-02)

- Core: `retune_drive(quantity, op, beams, reference; target, level, ...)` with
  `CapsuleModes`, `FaceFlux` or `FoilIrradiance`; `target = :reference` restores the
  reference's whole shape (modes, or the flux over the elements), `:uniform` asks for a
  uniform drive at `level`. The level is held (secant on a shifted penalty target), not
  a soft penalty. `drive_history(...; observed_flux = true)` keeps the face's flux.
- GUI: `retune_pulses` picks what it holds (`kind`): an observed sphere (capsule), an
  observed flat face (`_observed_flat`), else the active foil within the drive's spot
  radius (laser only; the drive needs no radiation model). The face's and the foil's
  `reference`, `before` and `after` are `(; time, mean, peak, low, rms)`, the face's from
  `_face_stats` on each history step, the foil's from depositions over `_spot_elements`
  (shared with `_spot_stats`). The retune chart's right axis then shows the RMS (before,
  after, reference) and the mean's departure from the reference's (before, after).
- The GUI's retune holds the reference's shape; the uniform target at a prescribed level
  is core only.
- The foil's GUI retune is the 3ω light's: it takes no unconverted light (`retune_pulses`
  throws if given), though core `retune_drive` adds each color's absorbed power to a
  foil's irradiance when the reference carries it.

## Item 5 as built (2026-11-02): presets in the starters' frames

- The hohlraum presets (`HOHLRAUM_PRESETS`: `nif_quads`, `nif_beams`, `omega_hohlraum`)
  aim in the entrance holes' frames (`preset_pointing`): beams with port θ < 90° in the
  hole that faces up (`<part>.leh_top` or `<part>.leh_bottom` whose `z` points up), the
  others in the hole facing down (`facing_lehs`). Of several holes facing one way the
  highest (lowest) wins, a hole under another hohlraum being shadowed by it; a tie goes
  to the first part. Closed and open ends are not holes. The shock-tube pair thus sends
  the upper beams to `halfraum_top.leh_top` and the lower to `halfraum_bottom.leh_bottom`,
  and moving a hohlraum part carries its beams.
- Offsets: the NIF cone pointing keeps its offsets from the default hohlraum part's holes
  at TCC (the pointing is for the 10.1 mm default), so a hohlraum of another length gets
  the pointing at the same offsets from its holes (before, at fixed chamber points). The
  OMEGA beams aim at the holes' centers (offset zero), as before for a hohlraum centered
  on TCC. When the hole's frame equals the reference frame (the default hohlraum at TCC),
  the aim is the reference aim itself, not a round trip through the frame.
- Bit for bit: the default design's chamber aims equal the cone pointing's (`isequal`,
  signs of zero included) and its drive (modes, mean flux, H, beam powers) equals the
  drive of the same beams held at TCC; the OMEGA hohlraum aims equal the former formula's
  (test, `test_halfraum.jl`).
- The opacity halfraum (decision taken): its one hole faces up, so the upper hemisphere
  drives it and the lower beams are disabled, aimed as without a hohlraum (at TCC). The
  lower beams would reach the sample's back through the open end, not the cavity.
  Without any hole (the foil, an empty target) every beam aims in `:tcc` as before.
- When it applies: loading a preset (`load_preset!`, the session's default table), File ▸
  New target when the beams of a hohlraum preset do not aim in the new target's holes as
  the preset would (a halfraum pair after a hohlraum: the preset is loaded again; a new
  hohlraum after a hohlraum keeps the beams and their edits), and Insert of a starter
  into a design whose beams drive no hole (`insert_starter!` on the session, then
  `aim_in_lehs!`): the beams held at TCC move to the new holes keeping their offsets from
  the reference holes (edits made at TCC carry over), enabled flags untouched. Inserting
  next to a hohlraum the beams already drive leaves them with it.
- Sessions: aims are saved with their frames and local coordinates and loaded through
  `set_beams!`, never re-aimed; a session saved before this change (aims at TCC) loads
  with its aims at TCC (test).

## Phase D as built (2026-11-02): what departs from the design

- Core: `HohlraumProfile(z, r)`, breaks from `−L/2` to `L/2` with a wall radius each;
  between breaks a cylinder (equal radii), a cone (unequal) or a flat ring (two breaks at
  one height: a step). `frustraum_profile(; length, waist_radius, end_radius,
  bottom_end_radius)`. `hohlraum(; profile)` takes no `length` or `radius` with it.
- One wall piece per section (decision taken: a `Component` carries one exact surface):
  `wall_1`, `wall_2`, … from the bottom, each with `part = :wall`, so `:wall` still names
  the whole wall in `elements`, `drive_history(; wall)`, clearance and the other
  name-taking analyses. A profile of one cylindrical section is the cylindrical hohlraum
  exactly (`:wall`). In the GUI the pieces are `<part>_wall_1`, …, all holding the wall
  role (`_hohlraum_walls`).
- Each end's cap runs from its hole to that end's radius; ends, lips, closed and open
  ends, the capsule (checked against every section in its meridian half-plane) and
  sombreros (the bore must clear the widest radius) work as on the cylinder. Windows are
  cut into cylindrical sections only, each within one; a window on a cone or across a
  break is refused, and so are corner fillets.
- Bit for bit: the cylinder path was refactored into shared helpers; the phase C
  snapshots (core: nine `hohlraum` configurations and a coarse NIF drive; GUI: five
  hohlraum parts) serialized before and after are identical. Validation: a closed
  frustraum against coaxial-disk closed forms, and LEH loss and capsule drive against
  Monte Carlo (`test/profiles.jl`; research notes, "Hohlraum profiles").
- GUI: `HOHLRAUM_VALUES` gains `profile` (0 cylinder, 1 frustraum, 2 sections,
  `HOHLRAUM_PROFILES`), `end_radius` and `bottom_end_radius` (a frustraum's; its waist is
  `radius`), and `PartSpec.sections` holds the breaks `(z, r)`, the end heights following
  the length (`hohlraum_profile`). Switching fills what is missing: a frustraum's end
  radius (60 % of the waist, at least 1.2 LEH radii), a section list from the current
  profile's breaks (a cylinder in four equal sections: the tentless four-section
  hohlraum's sections are a fabrication split with no public radii). The properties page
  has a *Profile* choice and, for sections, a table of breaks (z, r, remove; *Add break*
  splits the longest section), explained in tooltips and Help ▸ Hohlraum profiles.
- Starter `:frustraum` (*Frustraum with capsule*, Insert and File ▸ New target): the
  stated size of Baker et al. 2023's NIF frustraum (12.24 mm long, waist radius 4.88 mm,
  holes of 1.55 mm), end radius 2.282 mm derived from its 23° wall angle, capsule radius
  1.281 mm; a mock-up of the shape (its wall area is not the paper's quoted 270 mm²). The
  NIF presets aim in its entrance holes with the cylinder's offsets (item 5), not a
  frustraum pointing.
- Sessions store `sections` beside `windows` and the new values by key; missing keys take
  their defaults, so older sessions load unchanged (test) and the version stays.
- Window frames and the canvas's window outlines use the wall radius at the window's
  height (`_wall_radius_at`).
- Not built: windows on cones; fillets on a profile; the cutaway's default components
  (`DEFAULT_CUT_COMPONENTS`, the starter's `hohlraum_wall` and caps) do not name section
  walls, so a frustraum is not cut away by default; a frustraum pointing preset; the
  tentless membranes (not geometry here); GUI screenshots not regenerated.

## Phase C2 as built (2026-11-02): what departs from the design

- Core: `wall_temperature(sol, names...)`, the radiation temperature of the area-weighted
  mean incident flux on the named components; `drive_history(; cavities)` takes pairs
  `name => wall names` and returns `cavity_Tr`, one series per cavity (empty by
  default). A cavity's walls are named piece by piece, never by `Component.part`: a
  hohlraum part's entrance holes carry the part's name too (the trap of d5e6606).
- A cavity is a hohlraum part (decision taken): `cavity_walls(ts)` gives each included
  hohlraum part with its wall pieces (each section) and caps. A part holding the wall
  role without being a hohlraum part (a shield, a VisRad hohlraum) still counts in the
  history's combined `wall_Tr`, not as a cavity of its own.
- The Drive panel's result row stacks *Hohlraum wall* under the observed surface's drive:
  the chosen cavity's mean incident flux and Tr, a choice box with two or more cavities
  (tooltip; Help ▸ Hohlraum wall). The default (`drive_cavity`) is the cavity holding
  the observed surface (`observed_cavity`: the observed part's centroid inside the
  hohlraum part's profile, ends included, in the part's frame; a package face closing a
  halfraum's open end lies in it), else the first. The choice is a display choice in
  `DriveState.wall_cavity` (decision taken): not saved with the session and not undone;
  the drive and the history compute every cavity.
- The drive chart adds a dotted wall-Tr series per cavity when there are two or more
  (colored as each cavity's wall); with one cavity it is unchanged, the observed surface
  the headline (decision 8). `DriveResult.cavities` and `DriveHistory.cavities` hold the
  rows. The panel's heading reads *Observed face* without the empty mode table for a flat
  observed face, and the common choice box sizes its list for lists from Julia.
- Dante: `DanteView`'s collimated field of view, its center and the reference area
  already existed (diagnostics phase 2, 2026-10-29), as did the sourced field: ±30 mm at
  TCC, the 2004 design (Dewald et al. 2004, research notes "Dante-family hardware and
  geometry", H); the fielded collimators are not public. New: the built-in `nif_dante1`
  is Dante 1's channels and line of sight over a 60 mm field at TCC, Tr over the lower
  entrance hole: what the instrument sees (2026-11-02; C2 first kept `nif_dante1` as the
  view through the hole, beside a separate `nif_dante1_field`, so that existing Dante
  results would not change; the user rejected that, and the one definition replaced
  both). The dated Dante 1 configurations take the same field (46d8a76); Dante 2 and
  Omega's Dante look through the hole (their fields are not public).
- Bit for bit: the default design's drive (modes, mean flux, H, B, powers, map) and drive
  history (every core field, the GUI's Tr, modes, fields, ranges, energy, diagnostics)
  equal those computed without cavities, the path before this phase (test); its one
  cavity's wall Tr equals the history's `wall_Tr` exactly. Core: a history without
  `cavities` equals one with them in every other field.
- Validation (`test/cavities.jl`; GUI `test_halfraum.jl`): each cavity's Tr against Σ H A
  / Σ A over its wall elements found by component name in the test, per step of the
  core history's recorded flux, to 1e-12 (two hohlraums side by side driven unequally;
  the shock-tube starter with half its lower beams off). Dante: a field of 8 mm about
  one hohlraum receives exactly nothing from a blackbody in the other; the field reads
  the view through that hohlraum's hole to 1.2e-4, a 60 mm field both holes' sum to
  1.4e-4; on the shock tube `nif_dante1` reads the view through the lower hole to
  3.5e-4.
- Not built: per-cavity albedo or energy
  series; a Dante field centered on a part from the GUI (a definition's `center_part`
  serves).

## Piece roles as built (2026-11-02): decision 7's second half

- `PartSpec.piece_roles` (a hohlraum part's only): piece name without the prefix => role,
  holding only the roles that differ from the default (`piece_role`: `:wall` for the wall,
  each section, and the caps; `:none` for the holes, lips and windows). Keyed by the
  piece, so a rename keeps them; a role on a piece the part no longer has (a cap of an end
  since opened) is held by nothing until the piece is back.
- `PIECE_ROLES = (:none, :observed, :wall)` (decision taken): the foil role stays on parts.
  Only a cap, flat, can be observed (a wall section or a window is curved, and the
  capsule-only views need a sphere part); the entrance holes, lips and open windows are
  apertures and take none; a closed window can be a wall. `set_piece_role!(ts, component,
  role)` takes the component name (`hohlraum_cap_top`); one observed surface at a time
  (`_drop_observed!` clears parts' flags and pieces' observed roles alike).
- `observed_name` is then the piece's component name, never its part's (`part_pieces` on
  the part's name would take every piece, the trap of d5e6606); the scenario's
  `_observed_parts` gives core `observe` a part only when all its components are observed,
  else the components by name. `observed_flat` is true for an observed piece; `wall_names`,
  `cavity_walls` and `wall_components` take the pieces holding `:wall`, so an observed cap
  leaves its cavity's walls. `observed_cavity` and `observed_components` take a piece's
  name.
- Energy: core `drive_history(; per_component = true)` keeps `laser`, `source`,
  `absorbed`, `escaped`, `transmitted` and `released` per component instead of summing a
  part's pieces; the GUI's history passes it only when the observed surface is a piece, so
  the energy chart's observed-surface heat is the cap's and the walls' the rest. Default
  `false`: every other history is unchanged.
- `_closed_enclosure` now also counts parts lying wholly inside the single hohlraum (every
  vertex in its closed volume): the capsule that gave the role to a cap keeps the default
  design closed. Before, any free part without the observed role opened it (the test of a
  sphere without a role inside the hohlraum changed accordingly). The ghost check (a piece
  of the hohlraum or the observed part left out leaves a hole) is as before.
- GUI: a piece's properties have a *Role* box, and its tree row's right-click menu the
  three roles (those it cannot take disabled); the row's detail reads `observed`, `wall` or
  `no role` when the piece's role is not its default; Help ▸ Pieces ▸ A piece's role.
  Undo comes with the session state. Sessions: version 11 stores `piece_roles` per part;
  older sessions load with the defaults (test).
- Consumers checked: the drive's scenario and face map, the face statistics and retune
  kind (`_observed_flat` from the elements' normals: a cap is flat), the capsule keep-out
  and symmetry objective (refuse a flat face), the history (Tr from the face, the cavities
  from their walls, the energy), the clearance and TaLIS roles (`_held_sets`: an observed
  cap is a landing surface the TaLIS check leaves out, as any observed surface), the
  optimizer and sweep (`observed_name`). Dante's field and `ApertureKeepOut` read no role
  (geometry and apertures only).
- Bit for bit: the default design's and the shock tube's drive (modes, mean flux, H, B,
  powers, map, face, cavities, summary) and drive history (every core field, Tr, Legendre
  rows, face, fields, ranges, energy, diagnostics, cavities, the drive and energy charts),
  serialized before (fe18c7f) and after: identical.
- Validation (GUI `test_halfraum.jl`, "piece roles"): the default design with its top cap
  observed: the cap's mean flux against Σ H A / Σ A over the elements of the component
  named `hohlraum_cap_top`, found in the test, and the cavity's wall Tr against the same
  over `hohlraum_wall` and `hohlraum_cap_bottom`, to 1e-12, in the drive and per step of a
  core history; the cap's heat plus the hohlraum's other pieces' equals the by-part
  budget's hohlraum heat to 1e-12.

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
- A halfraum's open end meeting a package needs no gap between components: checked and
  fixed (see "Package seams as built"). The enclosure checks (`closed = true`) still
  leave two cavities joined by a package unbalanced (open: balancing needs the per-cavity
  enclosures, not one hohlraum part).

## Package seams as built (2026-11-02)

- The shock-tube starter's package disks have the halfraums' radius (1.5 mm) and lie in
  their open ends' planes exactly (z = ±1.0 mm): no overlap and no axial gap. The opacity
  starter's sample is a 1 mm square inside the open end, not a closure (what misses it
  leaves, by design), so it has no seam.
- The gap was azimuthal: the GUI meshed disks, washers, cylinders and cones with counts
  in multiples of 4, core `hohlraum` its wall in multiples of 8, so at the standard
  element size (0.28 mm, the default) the disk's rim was a 36-gon inside the wall's
  40-gon: 6.8e-3 mm², 9.6e-4 of the end (the inscribed polygons' areas). The coarse and
  fine sizes happened to agree (16 and 56).
- Leak, the upper halfraum closed by its disk alone (unbalanced view factors, wall
  driven, energy balance exact to 5e-15): Σ A(1 − ΣF) out of the enclosure 3.8e-4 mm²
  larger with the gap, the power lost 2.3e-5 of the source larger (−6.437e-4 → −6.672e-4:
  the view factors' own closure error, a quadrature excess, is the larger term).
- Fixed at the source (decision taken): the GUI's circular parts take core's rule
  (`_nphi`, multiples of 8), so a part of a hohlraum part's radius and element size meets
  the wall's ring vertex for vertex; after it the starter's closed configuration equals
  the one with the disk remeshed at the wall's count, exactly. Results change for any
  disk, washer, cylinder or cone meshed at a size where the two rules differ (the
  shock-tube starter at the standard size among them: 1760 → 1840 elements, the observed
  face's mean flux +1.0e-4 relative, its RMS 4.705 % → 4.710 %, each halfraum's wall Tr
  311.2789 → 311.2875 eV); the default design has none.
  Test: the seams' vertex azimuths equal at the coarse and standard sizes
  (`test_halfraum.jl`). Research notes, "Halfraum ends and starters".
