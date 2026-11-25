# The target consists of parts placed anywhere, the starters that insert assemblies of them,
# and the mesh element sizes. The parts are hohlraums, spheres, disks, washers, rectangles,
# cylinders and cones, as in the VisRad target components, and CAD meshes and VisRad
# objects. Every edit rebuilds the `Target`; meshing takes milliseconds, and the view
# factors are computed on request.

"""
    MeshSpec(h, capsule_h)

Element sizes (mm): `h` for hohlraum walls, end caps and LEHs (`hohlraum_counts`) and
the parts; `capsule_h` for the surface whose drive is studied: a sphere meshed as a
capsule (`capsule_counts`, see [`SPHERE_VALUES`](@ref)) and a starter's backlighter foil
when it is inserted. A part may override `h` with a separate value.
"""
struct MeshSpec
    h::Float64
    capsule_h::Float64
end

"""
Mesh presets (`MeshSpec`s), from the work-precision studies in
`examples/mesh_convergence.jl` and `examples/capsule_mesh.jl` (reference/research-notes.md,
"Hohlraum mesh work-precision" and "Compact capsule mesh") on the NIF hohlraum, whose
capsule is the geodesic `capsule_mesh = :compact`. The wall sizes set capsule P2/P0
errors against the converged meshes of about 0.5 percentage points (`:coarse`, view
factors in under a second), 0.1 (`:standard`, a few seconds) and 0.02 (`:fine`, tens of
seconds: converged for design work). The capsule sizes keep its share within those:
P2/P0 within about 0.17, 0.05 and 0.02 points and P4/P0 within 0.09, 0.03 and 0.01
(320, 980 and 2,880 capsule elements; P2 falls as 1/N on any capsule mesh, P4 ten times
faster on this one than on the lat-long mesh, which needed 0.2, 0.09 and 0.05 mm).
"""
const MESH_PRESETS = (coarse = MeshSpec(0.6, 0.3), standard = MeshSpec(0.28, 0.16),
                      fine = MeshSpec(0.17, 0.09))

"Mesh presets in order; `:custom` marks sizes edited by hand."
const MESH_RESOLUTIONS = (:coarse, :standard, :fine)

_preset_of(m::MeshSpec) = something(findfirst(==(m), MESH_PRESETS), :custom)

"""
Closing entrance holes of a hohlraum part (`Briareus.hohlraum(; leh_closure,
leh_open_radius)`), among its values ([`HOHLRAUM_VALUES`](@ref)): `:leh_closure` (1 on,
0 off) and these fields `(key, label, unit, default)`, edited with
[`edit_part!`](@ref). On, each entrance hole is open out to `leh_open_fraction ×` the
LEH radius and its lip beyond closes to x-rays at `leh_closure_time` over
`leh_closure_width` (a `GapClosure` whose closed albedo is the wall's). 0.84 is the
measured LEH radius at peak (research notes, "Closing entrance holes"); the time and
width are an engineering example (closure at the rise to peak), to be set per design.
"""
const LEH_CLOSURE_FIELDS = ((:leh_closure_time, "Closure time", "ns", 5.0),
                            (:leh_closure_width, "Closure duration", "ns", 0.5),
                            (:leh_open_fraction, "Open radius", "× LEH radius", 0.84))

# The closing lip of an entrance hole (`:hohlraum_leh_top_lip`, …) is a split of the
# hole for the calculations: the GUI lists, selects, highlights and leaves it out with
# its hole.
const _LIP_NAME = r"^((?:[A-Za-z0-9_]+_)?leh_(?:top|bottom))_lip$"
"Whether `c` is a hohlraum part's closing entrance-hole lip."
is_lip(c::Component) =
    c.kind === :aperture && c.material.albedo isa GapClosure && occursin(_LIP_NAME, String(c.name))
"The component a name is shown as: an entrance hole for its lip, else itself."
shown_component(name::Symbol) = (m = match(_LIP_NAME, String(name)); m === nothing ? name : Symbol(m[1]))
"Components `names` with the lips of the entrance holes among them."
with_lips(names) = Symbol[n for x in names
                          for n in (occursin(r"(^|_)leh_(top|bottom)$", String(x)) ? (x, Symbol(x, :_lip)) : (x,))]
