# Frames and placement

A [`Frame`](@ref) is a rigid frame in the chamber: its axes as the columns of a rotation
and its origin, so that local coordinates map to the chamber by
[`to_chamber`](@ref) and back by [`to_frame`](@ref). [`compose`](@ref) places one frame
in another, and `inv` inverts one. Frames come from closed-form constructors:
[`euler_rotation`](@ref) and [`euler_frame`](@ref) (VisRad's polar, azimuthal and
rotation angles), [`port_frame`](@ref) (a port's axis), [`foil_frame`](@ref) (the frame
in which [`backlighter_foil`](@ref) places its foil) and [`component_frame`](@ref) (a
part of a scripted target, about its [`part_centroid`](@ref)). [`coords_of`](@ref) and
[`point_of`](@ref) convert between chamber points and Cartesian, cylindrical or
spherical coordinates in a frame.

[`FrameDef`](@ref) places a named frame in a parent frame by an offset and three angles,
and [`resolve_frames`](@ref) resolves a set of them, depth first, from base frames such
as target-chamber center. A definition's [`PlacementStep`](@ref)s then move it relative
to other frames, in order, each in closed form: relations ([`Coincident`](@ref),
[`OnAxis`](@ref), [`Concentric`](@ref), [`Distance`](@ref), [`Facing`](@ref),
[`Parallel`](@ref), [`Roll`](@ref)) fix degrees of freedom, and relative moves
([`Offset`](@ref), [`Rotate`](@ref)) fix none. A relation that undoes an earlier one is
an error naming both; [`placement_report`](@ref) gives the degrees of freedom each step
fixes, redundant steps, what stays free, and the residuals left by relative moves. The
steps run in sequence rather than through a simultaneous solver: every relation has a
closed form, the result is deterministic, and an error names one step.

[`FramePointing`](@ref) parameterizes beam aims as local coordinates `(u, v, w)` of a
frame per beam, with the non-aim fields of [`GroupPointing`](@ref): a backlighter's
beams optimized across the foil's face, or a hohlraum's cones in the frames of their
entrance holes, one local aim giving mirrored chamber aims. It serves
[`Scenario`](@ref), [`subset`](@ref) and [`optimize_pointing`](@ref).

[`placement_sweep`](@ref) steps a design through a sequence of placements, such as a
foil moved along its line of sight with its beams following, and checks each step for
beam clearance, an unobstructed line of sight and beam capture, optionally
re-optimizing the pointing with a warm start. [`SweepResult`](@ref) holds one
[`SweepRow`](@ref) per step, and [`write_sweep_csv`](@ref) writes them. The
[Backlighter placement sweep](@ref) tutorial runs the whole sequence.

```@autodocs
Modules = [Briareus]
Pages = ["src/frames.jl", "src/placement.jl", "src/frame_pointing.jl", "src/sweep.jl"]
Private = false
```
