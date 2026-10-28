# Frames, placement and frame-relative pointing: design (2026-10-28)

The user asked for component frames and relative positioning that hold up in practice:
- optimize a backlighter's beam pointing in the foil's own coordinates, then step the
  foil along an axis and check for clipping;
- place a foil relative to other components (facing a radiography snout, parallel to
  another component's surface).

It draws on VisRad's coordinate transform tool
(https://prism-cs.com/Manuals/VisRad/additional_topics/coord_transform_tool.html) and on
Solidworks constraints. This file records the design and the phased plan. Remove each
phase from it once that phase lands; the docstrings and research notes then hold what
was built.

## Today (before phase 1)

- **Part frames.** Part frames exist only in the GUI (`target_setup.jl`).
  - `PartSpec.anchor` names a frame: one of `PART_ANCHORS = (:tcc, :leh_top,
    :leh_bottom, :foil)` or another part. `anchor_frame` returns `(R, origin)`, and
    `_visiting` catches cycles.
  - Placement is `x, y, z, theta, phi, rotation`, rotated Rz(φ)·Ry(θ)·Rz(ψ)
    (`_part_rotation`, the same order as core `_vrw_rotation`).
- **Beam aims.** `aim_ref` and `ref_point` (`beam_table.jl`) follow a component's
  centroid by translation only.
- **Foils.** `backlighter_foil` builds its frame (u, v, n) internally and discards it.
- **Optimizer.** It optimizes chamber-cylindrical aims (`GroupPointing`: `:z, :r,
  :phi`). A shared parameter set means different moves for beams from different ports,
  so the foil objective with shared parameters does not move ("GUI on CAD").
- **VisRad import.** Object frames are flattened into the chamber frame. Target-mount
  and diagnostic-port frames are refused, and so is pointing method 1.

## Model

1. **`Frame` (core, new `src/frames.jl`).**
   - `R::SMatrix{3,3}` holds the axes as columns and `origin::Vec3` the origin, so
     `x_chamber = R x_local + origin`.
   - Operations: `CHAMBER`, `to_chamber`, `to_frame` (points and directions), `compose`,
     `inv`.
   - Constructors:
     - `euler_rotation(polar, azimuth, rotation)` and `euler_frame`. These replace
       `_part_rotation` and `_vrw_rotation`.
     - `Frame(origin; z, x = θ̂ of z)`.
     - `port_frame(θ, φ; roll)`: z toward the port, x = θ̂. That x̂ is assumed, at
       confidence L, until a sample workspace shows it.
     - `foil_frame(los; …)`, onto which `backlighter_foil` is refactored, bit for bit.
     - `component_frame(t, part)` for script targets: the centroid, then the axis of the
       analytic surface, else the area-weighted normal, else the principal axes, with
       the ambiguity documented.
   - `coords_of` and `point_of` (Cartesian, cylindrical, spherical) move into core.
   - `to_chamber` is generic in the number type: the frame is Float64 and constant, and
     points may carry duals.
2. **`FrameDef`.**
   - Fields: `parent`, `offset`, `rotation` (polar, azimuth, rotation in the parent) and
     `steps`.
   - A part's frame and a user frame share one namespace. A port frame is a `FrameDef`
     on `:tcc`.
   - `resolve_frames(defs, base)` resolves depth first, with a loop error that names the
     cycle. It re-resolves everything on every edit, which costs microseconds.
   - Base frames come only from analytic sources: TCC, template frames, VisRad object
     frames. Frames never depend on beams or on the built target; beams depend on
     frames.
3. **Placement steps (phase 2).** Steps are applied in order, each in closed form:
   - `Coincident`, `OnAxis(ref; distance)`, `Concentric`, `Distance`;
   - `Facing(ref; point, side)` and `Parallel(ref; axis, flip)`, each a minimal rotation
     of the current frame;
   - `Roll(angle; from)`;
   - the relative moves `Offset` and `Rotate`.

   The steps run in sequence rather than through a general solver. Every requested
   relation has a closed form, the result is deterministic with explicit branches, and
   each error names one step. `PlacementReport` records:
   - the degrees of freedom each step fixes;
   - conflicts, as an error naming both steps and the size of the change;
   - redundant steps;
   - the free degrees of freedom (under-defined is allowed, and the rest keep their base
     values);
   - residuals after relative moves.
4. **Beam pointing in a frame.**
   - In the GUI, `aim_frame` and `aim_local` replace `aim_ref` and `ref_point`, with the
     invariant `aim == to_chamber(frame, aim_local)`. Every aim writer goes through
     `_set_aims!`. `follow_frames!` follows rigidly, rotation included. Defocus stays
     scalar along −port.
   - In core, `FramePointing(beams, frames; by, fields = (:u, :v, :w, :defocus, :power,
     :spot, :dlambda), relative)`:
     - `apply` computes `aim = R*SVector{3,T}(u,v,w) + origin`;
     - a shared local aim in `:leh_top` and `:leh_bottom` gives mirrored pointing
       without a flag;
     - `GroupPointing`'s non-aim tail moves into a shared helper, bit for bit;
     - `bounds`, `initial_parameters` and `parameter_names` as in `GroupPointing`;
     - inference and JET cases in `test/inference.jl`.
   - The GUI optimize fields gain `:u, :v, :w` (exclusive with the chamber aim fields).
5. **Sweep.**
   - Core `placement_sweep(design, values; clearance_kw, los, capture, optimize,
     callback)` returns a `SweepResult` with a row per step:
     - minimum clearance, the worst beam and component, fails and reviews;
     - line of sight: clear, blocked fraction, covered;
     - the minimum captured fraction;
     - the objective before and after.

     Re-optimizing is optional and warm-started. A script method moves parts along a
     frame axis, and beams marked `follow` move with them. `write_sweep_csv` writes the
     rows.
   - In the GUI, `SweepState` provides the settings (part or frame, axis, translate or
     rotate, range, which checks, quality), the job state, and
     `apply_sweep_step!`, which can be undone. `nudge_part!` steps by hand. The surface
     is a `:sweep` analysis node, `SweepProps.qml`, `SweepPanel.qml` (a table plus
     clearance and blocked fraction against position), and a ghost of the selected
     step.
6. **VisRad.**
   - `transform_point(frames, from, to, p; input, output)` and a Tools ▸ Coordinate
     transform dialog.
   - Pointing methods map: 0 → `:tcc`; 2 → that system's frame; 1 → the port frame with
     local (0, 0, distance).
   - Import parses diagnostic ports, keeps object hierarchies as anchored parts, and
     exports `Ref. Coord. Sys.`.
7. **Migration (no compatibility layer; old sessions convert on load).**
   - `anchor` becomes `frame`.
   - The `:foil` anchor becomes the `:backlighter` frame, with z along the front normal.
     Version ≤ 2 sessions are re-expressed with the old matrix inline, the only legacy
     code.
   - `aim_ref` becomes `aim_frame` (or `:tcc`). The saved chamber aims stay
     authoritative.
   - `component_points` moves to core as `part_centroid`.
   - `SESSION_VERSION` becomes 3.

## Phases

1. **Frames, frame pointing, sweep: the backlighter-foil case end to end.**
   - Core: **done (2026-10-28, uncommitted).** `src/frames.jl`,
     `src/frame_pointing.jl`, `src/sweep.jl`, the reference page "Frames and placement",
     the tutorial `backlighter_sweep.jl`, `test/frames.jl` (shard 3). Deviations from
     the model above:
     - `FrameDef` has no `steps` field yet; phase 2 adds it.
     - `FramePointing` lives in its own file after `optimize.jl`. `frames` is one
       `Frame`, a vector, or `beam -> Frame`; default fields `(:u, :v, :defocus)`.
       `subset` keeps the kept beams' frames and groups.
     - `placement_sweep`: the design callback returns `(; target, beams, pointing)`;
       the line of sight and capture take `los` + `los_kw` and `capture` +
       `capture_kw`; the convenience method translates parts only (no rotation), and
       `follow` is a Bool, indices or a predicate. The re-optimization box is `bounds`
       about each step's own initial parameters. Rows are `SweepRow`s.
     - `part_centroid(t, part)` is keyed by part (the GUI's `component_points` by
       component).
     - `_vrw_rotation` is gone: the workspace reader calls `euler_rotation`, the same
       matrices bit for bit (it was `cosd` already). The GUI's `_part_rotation` (radians
       through `cos`) differs from `euler_rotation` by ≤ 1e-15 absolute: the GUI half
       must accept that change in the last bits or keep its own.
     - `Frame` does not clash with CairoMakie, Documenter or Literate exports.
   - GUI: **done (2026-10-28, uncommitted).** Part frames through `resolve_frames`
     (`part_frames`, the `TargetSetup.frames` observable), the `:backlighter` template
     frame, `nudge_part!`, `aim_frame`/`aim_local` with `_set_aims!` and
     `follow_frames!`, the u, v, w fields, `SweepState` (`controllers/sweep.jl`,
     `SweepProps.qml`, `SweepPanel.qml`), session version 3, `test/test_frames.jl`.
     Deviations from the model above:
     - The invariant holds to round-off, `aim ≈ to_chamber(frame, aim_local)`: writers
       set chamber aims, which stay authoritative, and their local coordinates are
       derived; `follow_frames!` sets aims from local coordinates exactly.
     - `set_aim_frame!(sess, idx, frame; keep)`: `keep = false` aims at `frame_point`
       (a part's `part_centroid`, a template frame's origin), not at the frame's origin,
       as a CAD part's origin may lie anywhere.
     - The sweep's subject is a part or the `:backlighter` template frame; it moves the
       parts placed in its frames and the beams aimed in them. It runs the core's
       `design` form (template components move too; a `GroupPointing` is rebuilt per
       step, warm-started only for frame fields). Applying a step moves a part's center,
       or the template foil's `distance` (chamber x only). Re-optimization makes the
       objective per step's target (foil center, keep-outs, held drive from the step's
       start). The line of sight defaults to the foil → object direction. No ghost of
       the selected step yet.
     - A beam whose frame goes away (a part removed, or renamed by a role) keeps its aim
       and is held in `:tcc`, as `aim_ref` released it.
     - The part placement now uses `euler_rotation` (`cosd`) and `compose` (static
       matrices): part vertices move by ≤ 1.1e-15 mm against the old `_part_rotation`;
       version-2 parts in the old foil frame convert to within 9e-16 mm.
2. **Placement steps:** `PlacementReport`, a step list in the part inspector, and a
   sweep that varies a step's number.
3. **User frames:** port-axis frames, a Frames folder in the tree, the coordinate
   transform, and frame triads drawn in the view.
4. **VisRad frames,** import and export. (`_vrw_rotation` was unified in phase 1, bit
   for bit.)
5. **(Optional)** The foil template becomes a rectangle part with the foil role. That
   changes mesh counts: measure it and say so.

## Tests (independent references)

- **Frame algebra:**
  - `euler_rotation(θ,φ,0)ẑ == port_vector(θ,φ)`;
  - `compose(F, inv F) ≈ I`;
  - equality with the old `_part_rotation`;
  - `backlighter_foil` vertices unchanged.
- **Frame graph:** a loop error that names its path; a grandchild that follows its
  grandparent exactly.
- **Steps:** closed forms (Facing, Parallel, Concentric, Distance), degree-of-freedom
  counts, and conflict and redundancy reports.
- **FramePointing:**
  - in an identity frame, the aims equal `GroupPointing`'s;
  - in a rotated frame, R·u + o exactly;
  - the ForwardDiff Jacobian of the aims equals R's columns;
  - mirrored LEH frames give mirrored aims;
  - inference and JET cases.
- **Sweep:**
  - aims that follow move exactly Δ·axis;
  - a plate edge swept across a beam has clearance slope 1;
  - a foil receding behind an annular aperture starts blocking at the
    similar-triangles Δ*;
  - a canceled sweep returns partial rows.

## Risks

- Every aim writer must keep the `aim_local` invariant.
- The `:backlighter` frame's orientation differs from the old `:foil` anchor's, so the
  conversion must be exact.
- A sweep's cost: one `DepositionOperator` and one clearance per step. Default to the
  quick quality and keep optimization opt-in.
- VisRad conventions not yet seen: the port frame's x̂, the cylindrical order, target
  mounts.
