# Frames, placement and frame-relative pointing: design (2026-10-28)

The user asked for component frames and relative positioning that hold up in practice:
- optimize a backlighter's beam pointing in foil coordinates, then step the
  foil along an axis and check for clipping;
- place a foil relative to other components (facing a radiography snout, parallel to
  another component's surface).

It draws on VisRad's coordinate transform tool
(https://prism-cs.com/Manuals/VisRad/additional_topics/coord_transform_tool.html) and on
Solidworks constraints. This file records the design and the phased plan. Remove each
phase from it once that phase lands; the docstrings and research notes then hold what
was built.

Phases 1–3 and 5 have landed and are removed. Frames, frame pointing and the sweep are in
`src/frames.jl`, `src/frame_pointing.jl` and `src/sweep.jl`; placement steps are in
`src/placement.jl`; user frames, port frames and the coordinate transform are in the GUI
(a user frame is a `PartSpec` of kind `:frame`); the foil starter is a rectangle part
with the foil role. The research notes, "Frames and placement", record the conventions.
Phase 4 is open.

## Today (VisRad frames, before phase 4)

- **VisRad import.** Object frames are flattened into the chamber frame. A beam pointed
  by method 0 aims at TCC and one pointed by method 2 aims at its point in the frame
  named by `Ref. Coord. Sys.`, resolved at import (`src/visrad_workspace.jl`).
  Target-mount and diagnostic-port frames are refused, and so is pointing method 1.
- **VisRad export.** `write_visrad_workspace` does not write `Ref. Coord. Sys.`.
- **Rotations.** VisRad's rotation order is `euler_rotation` (`_vrw_rotation` was unified
  with it in phase 1, with results unchanged).
- **Port frames.** `port_frame(θ, φ; roll)` has z toward the port and x = θ̂. That x̂ is
  assumed, at confidence L, until a sample workspace shows it.

## Model (phase 4)

6. **VisRad.**
   - Pointing methods map: 0 → `:tcc`; 2 → that system's frame; 1 → the port frame with
     local (0, 0, distance).
   - Import parses diagnostic ports, keeps object hierarchies as anchored parts, and
     exports `Ref. Coord. Sys.`.

## Phases

4. **VisRad frames,** import and export.

## Risks

- VisRad conventions not yet seen: the port frame's x̂, the cylindrical order, target
  mounts.
