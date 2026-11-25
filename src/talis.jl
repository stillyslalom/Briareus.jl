# NIF TaLIS target-configuration checks. They cover 3ω capture and clearance, 1ω clearance
# and unconverted-light management, and the geometric hardware rules Briareus can
# represent. Each check is judged on the worst case over pointing, positioning and focus
# errors.

"Source of the rules implemented by `talis_check`."
const TALIS_SOURCE = "NIF \"TCC Configuration Requirements\" (NIF User Forum, 27 Oct 2018)"

"""
    TALIS_OUT_OF_SCOPE

Requirements of the TCC Configuration Requirements document that `talis_check` does
not evaluate, because Briareus has no model of the hardware involved.
"""
const TALIS_OUT_OF_SCOPE = [
    "positioner clearance (25 mm inside the CIVS NFOV, i.e. < 120 mm from TCC for TANDM-124 and < 150 mm otherwise; 35 mm outside): needs positioner and diagnostic geometry",
    "targets to TAS (5 mm above/below, 10 mm radially including insertion/retraction): needs TAS geometry",
    "TAS withdrawal envelope (cylinder of radius 60 mm for TARPOS 90-239 and TANDM 90-348, 29 mm for C-TARPOS 90-015; target height 38 mm): cylinder axis not specified",
    "diagnostics to the target shroud swept volume (10 mm): needs shroud geometry",
    "ARC beams must not intersect a NIF beam footprint on a target feature: ARC not modeled",
]

const _TALIS_BLOWBY = "TaLIS: configuration must be assessed for blow-by by the backscatter working group"
const _TALIS_LEH = "LEH clearance (an entrance-hole edge, not hardware): " *
                   "to be assessed for blow-by by the backscatter working group"

"""
    TaLISRow

One check of `talis_check`, for a (beam, rule, component) triple.

- `beam` — beam name (`"upper beams"`/`"lower beams"` for `:tas_stay_in`)
- `rule` — `:capture_3w`, `:incidence_3w`, `:clearance_3w`, `:clearance_1w`,
  `:management_1w`, `:diagnostic_1w` or `:tas_stay_in` (see `talis_check`)
- `component` — target component (`:none` when a beam misses every capture surface)
- `status` — `:pass`; `:fail` (the criterion is violated in the worst case); or
  `:review` (a TaLIS "evaluate further" outcome, or not decidable automatically)
- `margin` — worst-case value less the requirement, in `unit` (`:mm`, or `:deg` for
  angle rules); positive passes
- `worst`, `nominal` — the measured quantity in the worst case (errors and the
  tolerance applied) and at nominal (neither applied)
- `required` — the TaLIS threshold, in `unit`
- `point` — worst location (mm)
- `note` — consequence or explanation
"""
struct TaLISRow
    beam::String
    rule::Symbol
    component::Symbol
    status::Symbol
    margin::Float64
    unit::Symbol
    worst::Float64
    nominal::Float64
    required::Float64
    point::Vec3
    note::String
end

"""
    Omega1Reflection

The specular 1ω reflection behind a `:management_1w` row of [`talis_check`](@ref) that
was decided on reflection, as [`omega1_reflection`](@ref) returns it.

- `beam`, `component` — as in the row
- `beamlet_index`, `incident` — the beamlet whose reflection comes closest to an open
  port, and its direction
- `point`, `normal` — the centroid of the sub-triangle that reflects it, and the normal
  there
- `reflected` — the direction of the specular reflection
- `port`, `port_index` — the direction to the open port nearest that reflection, and
  its place among the open ports of the check
- `angle` — the angle (degrees) between `reflected` and `port`
- `spread` — the allowance (degrees) for the turn of the normal over the sub-triangle,
  twice the spread of the normal (0 on a flat mesh)
- `worst`, `margin` — `angle - spread` less twice the rule tolerance, and that less the
  required angle, as in the row
- `exposed_area` — the area (mm², projected along the beam) of the sub-triangles of the
  item that lie in the worst-case 1ω footprint and are judged on reflection; a
  sub-triangle that reaches the footprint anywhere counts whole
- `offending_area` — the part of `exposed_area` on sub-triangles whose reflection comes
  within the required angle of an open port

The areas separate a reflector from a sliver: the hazard is the reflected 1ω power that
enters another beamline, which scales with `offending_area` over the area of the
footprint.
"""
struct Omega1Reflection
    beam::String
    component::Symbol
    beamlet_index::Int
    incident::Vec3
    point::Vec3
    normal::Vec3
    reflected::Vec3
    port::Vec3
    port_index::Int
    angle::Float64
    spread::Float64
    worst::Float64
    margin::Float64
    exposed_area::Float64
    offending_area::Float64
end

"""
    TaLISReport

Result of `talis_check`, with `rows` (one `TaLISRow` per beam, rule and component), the
`settings` used (thresholds, error set, tolerance interpretation), `notes` on checks
that were not evaluated, and `reflections` (one [`Omega1Reflection`](@ref) per
`:management_1w` row decided on reflection). See `violations`, `write_talis_csv`.
"""
struct TaLISReport
    rows::Vector{TaLISRow}
    settings::NamedTuple
    notes::Vector{String}
    reflections::Vector{Omega1Reflection}
end

# The rule tolerance is a lateral allowance tan(α)·r, with r the distance of the point
# from the pivot.
struct _TaLISTol
    tanα::Float64
    pivot::Symbol            # :tcc, :aim, :max or :point
    point::Vec3
end

function _talis_tol(tolerance, pivot)
    tolerance >= 0 || throw(ArgumentError("tolerance must be ≥ 0 degrees"))
    if pivot isa Symbol
        pivot in (:tcc, :aim, :max) ||
            throw(ArgumentError("pivot must be :tcc, :aim, :max or a point"))
        return _TaLISTol(tand(tolerance), pivot, zero(Vec3))
    end
    return _TaLISTol(tand(tolerance), :point, Vec3(pivot...))
end

(τ::_TaLISTol)(x, aim) =
    τ.tanα * (τ.pivot === :tcc ? norm(x) : τ.pivot === :aim ? norm(x - aim) :
              τ.pivot === :max ? max(norm(x), norm(x - aim)) : norm(x - τ.point))

_pivot_text(τ::_TaLISTol) =
    τ.pivot === :tcc ? "r = distance from TCC" :
    τ.pivot === :aim ? "r = distance from the beam's aim point" :
    τ.pivot === :max ? "r = larger of the distances from TCC and from the aim point" :
    @sprintf("r = distance from (%.3g, %.3g, %.3g) mm", τ.point...)

"""
    talis_check(target, beams; capture, hardware = all surfaces, leh = (),
                omega1_components = all surfaces, diagnostics = (), tolerance = 0.5, pivot = :tcc,
                pointing = 0.05, target_offset = 0.025, hardware_tolerance = 0.0,
                focus = 1.0, ssd = 0.0, offsets = nothing, open_ports = :nif,
                curvature = Dict(), level = 2, resolution = 1e-4, rtol = 0.01, ...) -> TaLISReport

NIF TaLIS target-configuration checks from the "TCC Configuration Requirements"
(NIF User Forum, 27 Oct 2018, LMF Shot RI Resources). Every distance is measured
orthogonal to the beam, as the document states, and judged on the worst case over the
error set (and, except for incidence, the rule tolerance), never at nominal.

`capture` lists the components that beams are meant to land on, such as backlighter foils,
capture surfaces and hohlraum walls. The *landing component* of a beam is the first of
them its nominal axis meets. `leh` lists the components that form the edge of a laser
entrance hole the beams pass through (the end caps or LEH lips of a hohlraum, e.g.
`(:cap_top, :cap_bottom)` for [`hohlraum`](@ref)); they are checked for 3ω clearance like
`hardware`, whether or not `hardware` names them. Here and in `hardware`, `leh`,
`omega1_components`, `diagnostics` and `curvature`, a component may be named by its part:
a part split for computation (`curved_components`) is checked as one, and its rows are
listed under the part name.

# Rules

- `:capture_3w` — the `capture_contour` (10%)-of-peak contour must lie ≥ `capture_margin`
  (0.5 mm) inside the substrate edge of the landing component: the minimum over the edge
  of the transverse distance to the beam envelope (best-focus contour grown by the f-cone
  at the axial position of each point, see `BeamEnvelope`). Substrate edges are the free
  edges of the opaque mesh (edges with a single triangle; edges shared with another
  component are not edges). A beam whose axis misses every capture component fails with
  component `:none`.
- `:incidence_3w` — incidence < `max_incidence` (70°) to the target normal, the maximum
  over the surface of the landing component inside the worst-case 10% footprint (visible
  from the beam; exact normals on `SphereSurface`/`CylinderSurface` components). The
  angle is that of the steepest ray of the f-cone, the angle of the beam axis plus the
  cone half-angle atan(1/(2f#)) (1.38° at NIF), as VisRad's capture report measures it
  (VisRad 21 reports tilt + 1.384–1.390° on plates tilted 68–69.7° toward either spot
  axis). The incidence check includes the error set but not the rule tolerance; the
  capture margin and the clearances include both.
- `:clearance_3w` — the `clearance_contour` (1%)-of-peak contour must clear each
  `hardware` component other than the landing component by > `clearance_margin` (0.5 mm);
  hardware downstream of the landing point is skipped, and everything along the beam is
  checked if the beam misses every capture component. Capture and clearance failures have
  the TaLIS consequence, an assessment of blow-by by the backscatter working group. On an
  `leh` component the rule is an LEH-clearance assessment: not clearing it is a `:review`
  (on other hardware, a `:fail`; see "Entrance-hole edges").
- `:clearance_1w` — the 10% contour of the 1ω footprint must clear each
  `omega1_components` component by > `omega1_margin` (0.2 mm). Not clearing is a
  `:review` (the TaLIS outcome is to evaluate unconverted-light management) and is
  resolved by the next
  rule.
- `:management_1w` — for each surface the 1ω footprint does not clear: pass if its radius
  of curvature is < `cylinder_radius` (2 mm, cylindrical) or < `sphere_radius` (5 mm,
  spherical); otherwise the specular 1ω reflection from every uncleared facet must stay >
  `reflection_angle` (20°) from every open quad (`open_ports`), less twice the tolerance
  angle (a tilt by α turns a reflection by 2α). On a curved surface the angle of each
  sub-triangle is also reduced by twice the spread of the surface normal over it (the
  largest angle between the normals at its centroid and at its vertices), so a coarse mesh
  does not clear a port that lies between the reflections of its sampled normals. Both
  limits are on the unsigned radius: a concave surface diverges its reflection in the far
  field as a convex one does. Curvature
  comes from `curvature[name]` when given, and otherwise from the `CylinderSurface`,
  `SphereSurface`, `ConeSurface` or `TorusSurface` of the component or a planar mesh. On a
  cone and a torus the radii vary over the surface, and each sub-triangle in the footprint
  is judged separately on its largest radii: it scatters if its smaller principal radius is
  < `cylinder_radius` or its larger < `sphere_radius`, and it is given the reflection check
  otherwise. The radius of a cone is its distance from the axis over the cosine of the
  half-angle, so a tapered stalk may scatter near its tip and reflect near its base. The
  radii of a torus are its minor radius and `|R/cos ψ + r|` at the tube angle ψ. For a
  mesh with none of these, such as a CAD solid of a stalk and a flat flange, the
  radii of each facet are estimated from the turn of the vertex normals along its edges,
  where a vertex normal averages the facets around the vertex within 40° of the facet. A
  facet scatters if its estimate with a quarter added is below a limit, the fit residual
  is at most 0.1 and its normals spread by at most 25°; every other facet in the footprint
  is given the reflection check, with the spread of its vertex normals as on a curved
  surface, and a failing check is a `:review`. The estimate is exact on regular meshes of
  cylinders and within about 20% on spheres (the research notes). A `curvature` entry is
  `(:cylinder, r)`, `(:sphere, r)`, `:flat`, or `(:dimpled, r)` / `(:dimpled, r, half_side)` for a surface covered with spherical
  dimples of radius of curvature `r`, e.g. the dimpled 1ω "sombrero" shield around the LEH
  of an ignition hohlraum, modeled as a flat `annulus_mesh`. A dimpled surface is judged
  like a sphere against `sphere_radius`. The NIF User Guide (§5.6.2) also requires a
  dimple half-side/radius ≥ 0.25, so a dimpled surface passes only with `half_side` given
  and that ratio met, and is a `:review` otherwise. A flat component without a
  `curvature` entry is given the reflection check, which a horizontal plate usually fails
  (it mirrors each beam into the conjugate port of the same hemisphere).
- `:diagnostic_1w` — the `diagnostics` components (diagnostic hardware modeled in the
  target) must stay ≥ 10 mm from 1ω light at TCC, growing linearly to 25 mm at 1 m
  (`10 + 0.015 r` mm, r the distance from TCC; interpretation of "a cone of 10 mm at TCC
  extruded to 25 mm at 1 m"). With no `diagnostics` the rule is listed in the `notes` of
  the report as not evaluated. Diagnostic hardware given only as `hardware` is checked
  for the 3ω and 0.2 mm 1ω clearances but not for this rule.
- `:tas_stay_in` — (NIF beams) no component may extend more than `stay_in` (19 mm)
  above the lowest aim point of an upper beam or below the highest aim point of a lower
  beam (TAS alignment stay-in zone; violating it may require a different alignment
  sequence or beam selection). `hardware_tolerance` is added to the extent.

The remaining requirements of the document are listed in `TALIS_OUT_OF_SCOPE`.

# 1ω footprint

The 1ω rules apply to beams whose `quad` is a NIF quad name (or that have an `offsets`
entry). A NIF quad-level beam (port within 0.1° of its quad port, as from `nif_quads`)
is expanded into its four beamlets (LMF beam table) sharing its aim and defocus. The
1ω footprint of each beamlet at its 3ω best-focus plane is a square of side `omega1_size`
whose nearest edge is `omega1_edge` = 4.45 mm (the TaLIS 10% contour) from the 3ω
axis, displaced in φ̂ toward the beamlet side of its quad (the side of the quad
axis its port lies on; NIF User Guide Fig. 4-12 — *derived*, confidence medium). The
default side is the User Guide's 22 mm footprint grown by 4.8 − 4.45 = 0.35 mm per
side, so its near edge is the TaLIS 10% contour (22.7 mm). Away from that plane the
footprint is the cone converging to the 1ω focus `omega1_focus` (default 460 mm, from
`nif_unconverted`) beyond the 3ω focus, with its center on the chief ray from the lens
center `lens_distance` (7.7 m) upstream: at axial position `s` the near edge is
`4.45 + 0.027 s` mm, i.e. closer to the 3ω axis upstream of focus (*derived*). Pass
`omega1_focus = Inf, lens_distance = Inf` for a prism (the footprint fixed in the
transverse plane). `offsets` maps beam names to a footprint direction (degrees from the
port θ̂ toward φ̂, or a vector) and overrides the NIF convention for that beam. The
focus error moves the footprint with the 3ω focal plane; shadowing is ignored
(conservative).

# Entrance-hole edges

The TaLIS document states 3ω clearance "to target features" without distinguishing the
entrance holes of a target from hardware, and gives every 3ω violation the same
consequence: the configuration "must be assessed for blow-by by the backscatter working
group". The consequence is an assessment and does not prohibit the configuration.
Reporting LEH edges as `:review` while hardware proper, such as shields, stalks,
diagnostics and closed windows, stays `:fail` **is an assumption of Briareus; the TaLIS
document does not state it**. It is based on the fact that no NIF hohlraum with 3.1
mm-diameter LEHs meets the rule on its LEH edge. The best-focus 1% contour of an inner
(23.5°) quad has a semi-axis of 1.145 mm along φ̂, which is horizontal and orthogonal to
the beam. With the axis through the LEH center, 1.55 − 1.145 = 0.405 mm remains to the lip
before the f-cone, any error or the tolerance is counted; the rule requires 0.5 mm. Light
that reaches the lip is LEH clipping, which the deposition computes as such.

# Tolerance and error set

The document adds "up to 0.5 degrees metrology/alignment/positioning tolerance" to the 3ω
and 1ω clearances; Briareus also applies it to the 3ω capture margin. The angle is an
insertion or assembly error of the target component: the beams must still clear with the
component rotated by up to 0.5° off nominal in the least favorable direction, whether the
whole target is inserted rotated or the component is assembled off-square (user, from
facility practice, 2026-10-27; the document does not give a formula). It enters as a
lateral allowance `r·tan(tolerance)`, with `r` the distance of the hardware point from the
pivot of the rotation: the largest displacement of the point that such a rotation causes,
taken in full toward the beam (conservative). The default `pivot = :tcc` is the target
inserted rotated about TCC; a component assembled off-square rotates about its mount,
given as a point (`pivot = (0, 0, 12.0)`). `:aim` (rotation about the beam aim point) and
`:max` (the larger of the TCC and aim distances) are alternatives; with `tolerance = 0`
the allowance is dropped.

It is combined with the error set of `clearance` by linear addition (no RSS):

    worst = nominal − (pointing + target_offset + hardware_tolerance) − focus/(2F) − r·tan(tolerance)

evaluated at each point and minimized over the component (for the 1ω footprint the
focus term is the footprint shift over `±focus` instead of the f-cone). The
`target_offset` and `hardware_tolerance` partly overlap what the tolerance is meant to
cover, so this double counts conservatively; set them to 0 to rely only on the
tolerance. `ssd` (mm) widens the 3ω envelopes as in `clearance`.

Distances are minimized over the whole of each triangle (and each substrate edge), not at
sample points. Every distance in these rules is Lipschitz in the point, so a
branch-and-bound subdivision bounds the minimum from below. The Lipschitz constant is √(1
+ cone²) for the beam envelope, 1 plus its taper for the 1ω footprint and tan(tolerance)
for the tolerance allowance. The reported worst and nominal values are those lower bounds,
never above the exact minimum and within `resolution` (mm) or `rtol` times the distance of
the value from the rule threshold, whichever is larger: a bound never crosses the
threshold, so the status is exact (up to `resolution`), and only large margins are
reported coarser (1% by default). A plate crossed by a beam between its vertices fails.
The incidence and specular-reflection rules include each of the `4^level` sub-triangles of
a component (with its normal) that has any point inside the footprint, and whether a
sub-triangle reaches the footprint is determined in the same way. Clearance-only far-field
wings (`nif_wings`) are not used: the TaLIS rules refer to contours of the beam spot.
"""
function talis_check(t::Target, beams::AbstractVector{<:Beam}; capture,
                     hardware = nothing, leh = (), omega1_components = nothing, diagnostics = (),
                     tolerance = 0.5, pivot = :tcc,
                     pointing = 0.05, target_offset = 0.025, hardware_tolerance = 0.0,
                     focus = 1.0, ssd = 0.0,
                     capture_contour = 0.1, capture_margin = 0.5, max_incidence = 70.0,
                     clearance_contour = 0.01, clearance_margin = 0.5,
                     omega1::Bool = true, omega1_margin = 0.2, omega1_edge = 4.45,
                     omega1_size = nothing, omega1_focus = nothing, lens_distance = NIF_LENS_DISTANCE,
                     offsets = nothing, open_ports = :nif, reflection_angle = 20.0,
                     cylinder_radius = 2.0, sphere_radius = 5.0, curvature = Dict{Symbol,Any}(),
                     stay_in = 19.0, level::Integer = 2, resolution = 1e-4, rtol = 0.01)
    resolution > 0 || throw(ArgumentError("resolution must be > 0"))
    0 <= rtol < 1 || throw(ArgumentError("rtol must be in [0, 1)"))
    tol = _talis_tol(tolerance, pivot)
    R = pointing + target_offset + hardware_tolerance
    # Each checked item pairs a name with its components, the pieces of a part or one
    # component.
    surf = Pair{Symbol,Vector{Int}}[]
    for p in parts(t)
        cis = [ci for ci in part_pieces(t, p) if t.components[ci].kind === :surface]
        isempty(cis) || push!(surf, p => cis)
    end
    idx(names) = Pair{Symbol,Vector{Int}}[Symbol(n) => part_pieces(t, Symbol(n)) for n in names]
    cap_idx = idx(capture)
    hw_idx = hardware === nothing ? surf : idx(hardware)
    leh_idx = idx(leh)
    hw_idx = vcat(hw_idx, [g for g in leh_idx if !any(h -> first(h) === first(g), hw_idx)])
    o1_idx = omega1 ? (omega1_components === nothing ? surf : idx(omega1_components)) : empty(surf)
    dg_idx = omega1 ? idx(diagnostics) : empty(surf)
    geo = Dict(ci => _talis_samples(t, ci, level)
               for ci in unique(Int[ci for g in (cap_idx, hw_idx, o1_idx, dg_idx) for p in g for ci in last(p)]))
    pe = _free_edges(t, unique(vcat(Int[], last.(cap_idx)...)))
    edges = [vcat(NTuple{2,Vec3}[], (pe[ci] for ci in last(g))...) for g in cap_idx]
    cap_geo = [length(last(g)) == 1 ? geo[only(last(g))] : _merge_samples([geo[ci] for ci in last(g)])
               for g in cap_idx]
    o1 = omega1_contour(; edge = omega1_edge, size = omega1_size, focus = omega1_focus)
    f1 = o1.focal_shift
    vb = [_value_beam(b) for b in beams]
    table = any(_is_nif, vb) ? _nif_beam_table() : nothing
    bls = [omega1 ? _omega1_beamlets(b, offsets, table) : nothing for b in vb]
    ports = _talis_open_ports(open_ports, vb)
    cx = (; tol, R, focus = Float64(focus), ssd = Float64(ssd), cap_idx, hw_idx, o1_idx, dg_idx,
          leh = Symbol[first(g) for g in leh_idx],
          geo, edges, cap_geo, capture_contour, capture_margin, max_incidence = Float64(max_incidence),
          clearance_contour, clearance_margin, tolerance = Float64(tolerance),
          omega1_margin = Float64(omega1_margin), o1.c0, o1.h0, f1,
          L = Float64(lens_distance), ports, reflection_angle = Float64(reflection_angle),
          cylinder_radius = Float64(cylinder_radius), sphere_radius = Float64(sphere_radius),
          curvature, curv = Dict{Int,Vector{NTuple{4,Float64}}}(), curv_lock = ReentrantLock(),
          ε = Float64(resolution), rtol = Float64(rtol),
          lip1 = _omega1_lipschitz(o1.c0, o1.h0, f1, lens_distance))
    out = [TaLISRow[] for _ in vb]
    refl = [Omega1Reflection[] for _ in vb]
    Threads.@threads for bi in eachindex(vb)
        _talis_beam!(out[bi], refl[bi], t, vb[bi], bls[bi], cx)
    end
    rows = reduce(vcat, out; init = TaLISRow[])
    append!(rows, _talis_stay_in(t, vb, surf, Float64(stay_in), Float64(hardware_tolerance)))
    notes = String[]
    skipped = [b.name for (b, bl) in zip(vb, bls) if bl === nothing]
    omega1 && !isempty(skipped) && !isempty(o1_idx) &&
        push!(notes, "1ω rules for beams without a NIF quad name or an `offsets` entry: " *
                     join(first(skipped, 12), ", ") * (length(skipped) > 12 ? ", …" : ""))
    omega1 || push!(notes, "1ω rules (omega1 = false)")
    omega1 && isempty(dg_idx) && any(!isnothing, bls) &&
        push!(notes, "diagnostics to 1ω light (10 mm at TCC to 25 mm at 1 m): no `diagnostics` given")
    any(_is_nif, vb) || push!(notes, "TAS stay-in zone (no NIF beams)")
    settings = (; tolerance = Float64(tolerance), pivot = tol.pivot, pivot_point = tol.point,
                pivot_text = _pivot_text(tol), transverse = Float64(R),
                pointing = Float64(pointing), target_offset = Float64(target_offset),
                hardware_tolerance = Float64(hardware_tolerance), focus = Float64(focus),
                ssd = Float64(ssd), capture_contour = Float64(capture_contour),
                capture_margin = Float64(capture_margin), max_incidence = Float64(max_incidence),
                clearance_contour = Float64(clearance_contour),
                clearance_margin = Float64(clearance_margin), omega1_margin = Float64(omega1_margin),
                omega1_edge = Float64(omega1_edge), omega1_size = 2 * o1.h0, omega1_focus = f1,
                lens_distance = Float64(lens_distance), reflection_angle = Float64(reflection_angle),
                cylinder_radius = Float64(cylinder_radius), sphere_radius = Float64(sphere_radius),
                open_ports = length(ports), stay_in = Float64(stay_in),
                level = Int(level), resolution = Float64(resolution), rtol = Float64(rtol),
                capture = Symbol[first(g) for g in cap_idx], leh = cx.leh,
                diagnostics = Symbol[first(g) for g in dg_idx])
    return TaLISReport(rows, settings, notes, reduce(vcat, refl; init = Omega1Reflection[]))
end

# Per-beam checks (3ω capture, incidence and clearance; 1ω clearance, management and
# diagnostics). `bl` contains the 1ω beamlets of the beam, or is `nothing`.
function _talis_beam!(rows, refl, t::Target, b::Beam, bl, cx)
    aim = Vec3(b.aim)
    tol = cx.tol
    k = direction(b)
    env1 = BeamEnvelope(b; peak_fraction = cx.clearance_contour, margin = cx.ssd)
    env10 = BeamEnvelope(b; peak_fraction = cx.capture_contour, margin = cx.ssd)
    shift = cx.R + cx.focus * env1.cone          # exact worst-case transverse loss
    L3 = sqrt(1 + env1.cone^2)                   # Lipschitz constant of envelope_distance
    land, s_stop = 0, Inf                        # index of the landing item in cap_idx
    for (j, g) in enumerate(cx.cap_idx)
        s = _landing(t, env1, last(g); miss = Inf)
        s < s_stop && ((land, s_stop) = (j, s))
    end
    landed = land == 0 ? Int[] : last(cx.cap_idx[land])

    # 3ω capture (with the rule tolerance) and incidence (error set only)
    if !isempty(cx.cap_idx)
        if land == 0
            push!(rows, TaLISRow(b.name, :capture_3w, :none, :fail, -Inf, :mm, -Inf, -Inf,
                                 cx.capture_margin, aim,
                                 "beam axis misses every capture component; " * _TALIS_BLOWBY))
        else
            name = first(cx.cap_idx[land])
            dist10 = x -> envelope_distance(env10, x)[1]
            worst, wp = _lipschitz_minimum(x -> dist10(x) - shift - tol(x, aim), cx.edges[land],
                                       L3 + tol.tanα, (cx.ε, cx.rtol, cx.capture_margin), aim)
            nom, _ = _lipschitz_minimum(dist10, cx.edges[land], L3, (cx.ε, cx.rtol, cx.capture_margin), aim)
            m = worst - cx.capture_margin
            note = !isfinite(worst) ? "no free substrate edge" : m >= 0 ? "" : _TALIS_BLOWBY
            push!(rows, TaLISRow(b.name, :capture_3w, name, m >= 0 ? :pass : :fail, m, :mm, worst,
                                 nom, cx.capture_margin, wp, note))
            # The steepest ray of the f-cone is at the axis angle plus the cone half-angle.
            hc = atand(1 / (2b.fnumber))
            worst, θn, ip = _talis_incidence(t, cx.cap_geo[land], env10, k, aim, shift, L3, cx.ε)
            worst, θn = min(worst + hc, 90.0), min(θn + hc, 90.0)
            m = cx.max_incidence - worst
            push!(rows, TaLISRow(b.name, :incidence_3w, name, m > 0 ? :pass : :fail, m, :deg, worst,
                                 θn, cx.max_incidence, ip, m > 0 ? "" : _TALIS_BLOWBY))
        end
    end

    # 3ω clearance
    dist1 = x -> envelope_distance(env1, x)[1]
    axial = x -> dot(x - env1.focus, env1.k)     # s of envelope_distance
    for (name, cis) in cx.hw_idx
        any(in(landed), cis) && continue
        tris = _talis_tris(cx, cis)
        worst, wp = _lipschitz_minimum(x -> dist1(x) - shift - tol(x, aim), tris, L3 + tol.tanα,
                                   (cx.ε, cx.rtol, cx.clearance_margin), aim, axial, s_stop)
        nom, _ = _lipschitz_minimum(dist1, tris, L3, (cx.ε, cx.rtol, cx.clearance_margin), aim, axial, s_stop)
        m = worst - cx.clearance_margin
        lip = name in cx.leh
        note = !isfinite(worst) ? "downstream of the landing point" : m > 0 ? "" :
               lip ? _TALIS_LEH : _TALIS_BLOWBY
        push!(rows, TaLISRow(b.name, :clearance_3w, name, m > 0 ? :pass : lip ? :review : :fail, m,
                             :mm, worst, nom, cx.clearance_margin, wp, note))
    end

    bl === nothing && return rows
    # 1ω clearance and management
    fw = x -> _omega1_distance(x, bl, cx)[2] - cx.R - tol(x, aim)
    fn = x -> _omega1_distance(x, bl, cx)[1]
    for (name, cis) in cx.o1_idx
        tris = _talis_tris(cx, cis)
        worst, wp = _lipschitz_minimum(fw, tris, cx.lip1 + tol.tanα, (cx.ε, cx.rtol, cx.omega1_margin), aim)
        nom, _ = _lipschitz_minimum(fn, tris, cx.lip1, (cx.ε, cx.rtol, cx.omega1_margin), aim)
        m = worst - cx.omega1_margin
        if m > 0
            push!(rows, TaLISRow(b.name, :clearance_1w, name, :pass, m, :mm, worst, nom,
                                 cx.omega1_margin, wp, ""))
        else
            push!(rows, TaLISRow(b.name, :clearance_1w, name, :review, m, :mm, worst, nom,
                                 cx.omega1_margin, wp,
                                 "1ω footprint does not clear: unconverted-light management required (see management_1w)"))
            push!(rows, _talis_management!(refl, t, name, cis, b, bl, cx))
        end
    end
    # diagnostics to 1ω light
    for (name, cis) in cx.dg_idx
        tris = _talis_tris(cx, cis)
        best, wp = _lipschitz_minimum(x -> fw(x) - (10 + 0.015 * norm(x)), tris,
                                  cx.lip1 + tol.tanα + 0.015, (cx.ε, cx.rtol, 0.0), aim)
        nom, _ = _lipschitz_minimum(fn, tris, cx.lip1, (cx.ε, cx.rtol, 10.0), aim)
        req = 10 + 0.015 * norm(wp)
        worst = best + req
        push!(rows, TaLISRow(b.name, :diagnostic_1w, name, best >= 0 ? :pass : :fail,
                             best, :mm, worst, nom, req, wp,
                             "required clearance 10 mm at TCC → 25 mm at 1 m (linear in distance from TCC)"))
    end
    return rows
end

# Largest incidence angle (degrees) on a landing component inside the worst-case and
# nominal footprints, over sub-triangles visible from upstream that reach into them
# anywhere. Reaching is tested with `_talis_reaches`, and `lip` is the Lipschitz constant
# of `envelope_distance`. The sub-triangle nearest the beam axis is always included.
function _talis_incidence(t::Target, g, env, k, aim, shift, lip, res)
    far = 4 * (t.scale + norm(aim)) + 1
    ε = 1e3 * _ray_eps(t)
    θw, θn, ip = -Inf, -Inf, aim
    bestd, θc, pc = Inf, NaN, aim
    dist = x -> envelope_distance(env, x)[1]
    for (leaf, n) in zip(g.leaves, g.normals)
        d = minimum(dist, leaf)
        tri = (leaf[1], leaf[2], leaf[3])
        inw = d <= shift || _talis_reaches(dist, tri, lip, shift, res)
        (inw || d < bestd) || continue
        m = leaf[4]
        visible(t, m - far * k, m - ε * k) || continue
        θ = acosd(clamp(abs(dot(n, k)), 0.0, 1.0))
        d < bestd && ((bestd, θc, pc) = (d, θ, m))
        inw && (d <= 0 || _talis_reaches(dist, tri, lip, 0.0, res)) && (θn = max(θn, θ))
        inw && θ > θw && ((θw, ip) = (θ, m))
    end
    if isfinite(θc)
        θn = max(θn, θc)
        θc > θw && ((θw, ip) = (θc, pc))
    end
    return θw, θn, ip
end

# 1ω management of an uncleared item `name` made of the components `cis`. The result is the
# worst row of its components (status first, then margin), under the name of the item. If
# that row was decided on reflection, its reflection is pushed to `refl`, with the areas of
# every component of the item.
function _talis_management!(refl, t::Target, name::Symbol, cis, b::Beam, bl, cx)
    rank(r) = r.status === :fail ? 0 : r.status === :review ? 1 : 2
    best, bx = nothing, nothing
    xs = Omega1Reflection[]
    for ci in cis
        r, x = _talis_management(t, ci, b, bl, cx)
        x === nothing || push!(xs, x)
        (best === nothing || (rank(r), r.margin) < (rank(best), best.margin)) && ((best, bx) = (r, x))
    end
    r = best
    bx === nothing ||
        push!(refl, Omega1Reflection(bx.beam, name, bx.beamlet_index, bx.incident, bx.point, bx.normal,
                                     bx.reflected, bx.port, bx.port_index, bx.angle, bx.spread, bx.worst,
                                     bx.margin, sum(x -> x.exposed_area, xs), sum(x -> x.offending_area, xs)))
    return r.component === name ? r :
           TaLISRow(r.beam, r.rule, name, r.status, r.margin, r.unit, r.worst, r.nominal, r.required,
                    r.point, r.note)
end

# 1ω management of an uncleared component, by curvature or else by specular reflection:
# its row and, when the row was decided on reflection, the reflection.
function _talis_management(t::Target, ci::Integer, b::Beam, bl, cx)
    c = t.components[ci]
    aim = Vec3(b.aim)
    kind, rc, half = _talis_curvature(t, ci, cx.curvature)
    lim = kind === :cylinder ? cx.cylinder_radius :
          (kind === :sphere || kind === :dimpled) ? cx.sphere_radius : NaN
    if kind !== :flat && kind !== :unknown && rc < lim
        what = kind === :dimpled ? "spherical dimple" : string(kind)
        note = @sprintf("%s radius of curvature %.3g mm < %.3g mm scatters 1ω", what, rc, lim)
        # The NIF User Guide (§5.6.2) requires a dimple half-side / radius of at least 0.25.
        if kind === :dimpled && !isfinite(half)
            return TaLISRow(b.name, :management_1w, c.name, :review, lim - rc, :mm, rc, rc, lim, aim,
                            note * "; dimple spacing unknown: give (:dimpled, r, half_side) to check half-side/radius ≥ 0.25 (User Guide §5.6.2)"), nothing
        elseif kind === :dimpled && half / rc < 0.25
            return TaLISRow(b.name, :management_1w, c.name, :review, lim - rc, :mm, rc, rc, lim, aim,
                            note * @sprintf("; dimple half-side/radius %.3g < 0.25 (User Guide §5.6.2: periodicity too short)",
                                            half / rc)), nothing
        end
        return TaLISRow(b.name, :management_1w, c.name, :pass, lim - rc, :mm, rc, rc, lim, aim, note), nothing
    end
    g = cx.geo[ci]
    αw, αn, wp = Inf, Inf, aim
    fw = x -> _omega1_distance(x, bl, cx)[2] - cx.R - cx.tol(x, aim)
    fn = x -> _omega1_distance(x, bl, cx)[1]
    # On a cone or a torus the radii of curvature vary over the surface, and each
    # sub-triangle in the footprint is judged on the largest radii over it.
    varying = kind === :cone || kind === :torus
    # A mesh without a surface has its radii estimated from the turn of its normals. An
    # estimate counts where the fit is good and the facets are fine enough, with a quarter
    # added to the radius for its error.
    mc = kind === :unknown ? _component_curvature(t, ci, cx) : nothing
    nl = length(g.leaves) ÷ max(length(g.tris), 1)
    scat = (Inf, NaN, NaN, aim)                  # (limit − radius, radius, limit, point), the least
    gentle = 0.0
    kq = direction(b)
    exposed, offending = 0.0, 0.0
    wit = (0, kq, kq, kq, 0, NaN, 0.0)           # beamlet, k, normal, reflection, port, angle, spread
    for (i, (leaf, n, spread)) in enumerate(zip(g.leaves, g.normals, g.spreads))
        tri = (leaf[1], leaf[2], leaf[3])
        _talis_reaches(fw, tri, cx.lip1 + cx.tol.tanα, cx.omega1_margin, cx.ε) || continue
        if varying || mc !== nothing
            class, R1, R2, spread = _leaf_curvature(c.surface, kind, rc, half, mc, nl, i, tri, spread, cx)
            if class === :scatter
                v = R1 < cx.cylinder_radius ? (cx.cylinder_radius - R1, R1, cx.cylinder_radius, leaf[4]) :
                    (cx.sphere_radius - R2, R2, cx.sphere_radius, leaf[4])
                v[1] < scat[1] && (scat = v)
                continue
            end
            varying && (gentle = max(gentle, R1))
        end
        hitn = _talis_reaches(fn, tri, cx.lip1, 0.0, cx.ε)
        αl = Inf
        for (j, (kb, _, _)) in enumerate(bl)
            r = kb - 2 * dot(kb, n) * n
            # a normal turned by the spread turns the reflection by twice that
            raw, pj = findmin(p -> acosd(clamp(dot(r, p), -1.0, 1.0)), cx.ports)
            α = raw - 2 * spread
            αl = min(αl, α)
            if α < αw
                αw, wp = α, leaf[4]
                wit = (j, kb, n, r, pj, raw, 2 * spread)
            end
            hitn && (αn = min(αn, α))
        end
        ap = abs(dot(cross(leaf[2] - leaf[1], leaf[3] - leaf[1]), kq)) / 2
        exposed += ap
        αl - 2 * cx.tolerance <= cx.reflection_angle && (offending += ap)
    end
    if !isfinite(αw) && isfinite(scat[1])
        note = @sprintf("%s radius of curvature ≤ %.3g mm < %.3g mm over the surface in the footprint scatters 1ω",
                        varying ? string(kind) : "mesh estimate (a quarter added) of the", scat[2], scat[3])
        return TaLISRow(b.name, :management_1w, c.name, :pass, scat[1], :mm, scat[2], scat[2], scat[3],
                        scat[4], note), nothing
    end
    worst = αw - 2 * cx.tolerance
    m = worst - cx.reflection_angle
    desc = kind === :unknown ? "free-form mesh, the facets curved below the limits by the mesh estimate left out (give `curvature` to set it)" :
           kind === :flat ? "flat" :
           varying ? @sprintf("%s radius up to %.3g mm ≥ %.3g mm in the footprint", kind, gentle, cx.cylinder_radius) :
           @sprintf("%s radius %.3g mm ≥ %.3g mm", kind, rc, lim)
    status = m > 0 ? :pass : kind === :unknown ? :review : :fail
    note = m > 0 ? "$desc; specular 1ω reflection clears open quads" :
           "$desc; specular 1ω reflection within $(cx.reflection_angle)° of an open quad" *
           @sprintf(" from %.3g mm² of %.3g mm² in the footprint (projected along the beam)", offending, exposed) *
           ": TaLIS additional review, may not be allowed"
    row = TaLISRow(b.name, :management_1w, c.name, status, m, :deg, worst, αn,
                   cx.reflection_angle, wp, note)
    isfinite(αw) || return row, nothing
    j, kb, n, r, pj, raw, sp = wit
    return row, Omega1Reflection(b.name, c.name, j, kb, wp, n, r, cx.ports[pj], pj, raw, sp, worst, m,
                                 exposed, offending)
end

"""
    omega1_reflection(report::TaLISReport, row::TaLISRow) -> Union{Omega1Reflection,Nothing}

The specular reflection behind a `:management_1w` row of `report`: the reflecting point,
the incident and reflected directions, the open port that the reflection threatens and
the projected areas in the footprint. `nothing` for a row decided on curvature and for a
row of another rule.
"""
function omega1_reflection(r::TaLISReport, row::TaLISRow)
    row.rule === :management_1w || return nothing
    i = findfirst(x -> x.beam == row.beam && x.component === row.component, r.reflections)
    return i === nothing ? nothing : r.reflections[i]
end

# The unsigned principal radii of curvature (smaller, larger) of a cone or torus, each the
# largest over the triangle `tri`. On a cone they are z sin α / cos² α, the distance from
# the axis over cos α, at height z from the apex, and Inf along the ruling. On a torus
# they are the minor radius and |R / cos ψ + r| at the tube angle ψ from the outer
# equator. Both vary monotonically between the vertices, except that the second radius of
# a torus is unbounded where the triangle crosses a pole of the tube.
function _surface_radii(s::ConeSurface, tri)
    sα, cα = sincosd(s.half_angle)
    z = maximum(p -> dot(p - s.apex, s.axis), tri)
    return max(z, 0.0) * sα / cα^2, Inf
end

function _surface_radii(s::TorusSurface, tri)
    cs = map(tri) do p
        n = last(surface_point(s, p))
        dot(n, _torus_coordinates(s, p)[3])
    end
    lo, hi = extrema(cs)
    R2 = (lo <= 0 <= hi) ? Inf : maximum(c -> abs(s.major_radius / c + s.minor_radius), cs)
    return minmax(s.minor_radius, R2)
end

# Unsigned principal radii of curvature of the triangles of a mesh without a `CurvedSurface`,
# from the turn of the vertex normals along the edges of each triangle (the per-face
# shape operator of Rusinkiewicz, 3DPVT 2004). The normal at a vertex of a triangle is the
# angle-weighted mean of the facet normals around the vertex that lie within `crease`
# degrees of the normal of that triangle, so a normal does not average across a sharp
# edge. On a regular mesh of a cylinder or a sphere the vertex normals are the surface
# normals and the chords of the normals and of the points have the ratio 1/R, so the radii
# are exact.
#
# Returns, per triangle, (R1, R2, spread, quality): the radii (smaller, larger; Inf for no
# curvature), the largest angle (degrees) between the facet normal and its vertex normals,
# and the residual of the fit relative to the turn of the normals (0 for a flat triangle).
function _mesh_curvature(tris::AbstractVector{NTuple{3,Vec3}}, normals::AbstractVector{Vec3};
                         crease::Real = 40.0, scale::Real = 1.0)
    q = 1e-7 * max(scale, 1.0)
    key(v) = (round(Int, v[1] / q), round(Int, v[2] / q), round(Int, v[3] / q))
    fan = Dict{NTuple{3,Int},Vector{Tuple{Int,Float64}}}()      # vertex → (triangle, corner angle)
    for (k, tri) in enumerate(tris), i in 1:3
        a, b, c = tri[i], tri[mod1(i + 1, 3)], tri[mod1(i + 2, 3)]
        u, v = b - a, c - a
        lu, lv = norm(u), norm(v)
        (lu > 0 && lv > 0) || continue
        push!(get!(() -> Tuple{Int,Float64}[], fan, key(a)), (k, acos(clamp(dot(u, v) / (lu * lv), -1.0, 1.0))))
    end
    cosc = cosd(crease)
    out = Vector{NTuple{4,Float64}}(undef, length(tris))
    for (k, tri) in enumerate(tris)
        n = normals[k]
        vn = map(tri) do p
            m = zero(Vec3)
            for (j, w) in get(fan, key(p), Tuple{Int,Float64}[])
                dot(normals[j], n) >= cosc && (m += w * normals[j])
            end
            nm = norm(m)
            nm > 0 ? m / nm : n
        end
        spread = maximum(v -> acosd(clamp(dot(v, n), -1.0, 1.0)), vn)
        e1, e2 = _ortho_basis(n)
        # Least squares for the symmetric S = [a b; b c] with Δn = S Δp over the three edges,
        # both in the tangent plane of the facet.
        A = zeros(6, 3)
        r = zeros(6)
        for (row, (i, j)) in enumerate(((1, 2), (2, 3), (3, 1)))
            dp, dn = tri[j] - tri[i], vn[j] - vn[i]
            x, y = dot(dp, e1), dot(dp, e2)
            A[2row-1, :] .= (x, y, 0.0)
            A[2row, :] .= (0.0, x, y)
            r[2row-1], r[2row] = dot(dn, e1), dot(dn, e2)
        end
        turn = norm(r)
        if turn < 1e-9
            out[k] = (Inf, Inf, spread, 0.0)
            continue
        end
        s = A \ r
        res = norm(A * s - r) / turn
        tr, det = s[1] + s[3], s[1] * s[3] - s[2]^2
        disc = sqrt(max(tr^2 / 4 - det, 0.0))
        κ1, κ2 = abs(tr / 2 + disc), abs(tr / 2 - disc)
        R1, R2 = minmax(1 / κ1, 1 / κ2)
        out[k] = (R1, R2, spread, res)
    end
    return out
end

# The mesh curvature of the triangles of component `ci` (in the order of `cx.geo[ci].tris`),
# computed once per check.
function _component_curvature(t::Target, ci, cx)
    return lock(cx.curv_lock) do
        get!(cx.curv, ci) do
            ks = [k for k in eachindex(t.tri_area) if t.tri_comp[k] == ci && t.tri_area[k] > 0]
            _mesh_curvature(cx.geo[ci].tris, t.tri_normal[ks]; scale = t.scale)
        end
    end
end

# How sub-triangle `i` of a component is judged for 1ω management: `:scatter` (a radius
# below its limit), `:reflect` (the specular check) or `:review` (a dimpled surface whose
# spacing is unknown or too fine), with the radii used and the spread of the normals. The
# radii are those of the component for a sphere, cylinder or dimpled surface, the largest
# over `tri` on a cone or torus, and the mesh estimate with a quarter added on a mesh
# without a surface, where an estimate counts only if its fit is good and its facets fine.
function _leaf_curvature(surface, kind, rc, half, mc, nl, i, tri, spread, cx)
    if kind === :flat
        return :reflect, Inf, Inf, spread
    elseif kind === :cylinder
        return (rc < cx.cylinder_radius ? :scatter : :reflect), rc, Inf, spread
    elseif kind === :sphere
        return (rc < cx.sphere_radius ? :scatter : :reflect), rc, rc, spread
    elseif kind === :dimpled
        rc < cx.sphere_radius || return :reflect, rc, rc, spread
        return ((isfinite(half) && half / rc >= 0.25) ? :scatter : :review), rc, rc, spread
    elseif kind === :cone || kind === :torus
        R1, R2 = _surface_radii(surface, tri)
        return ((R1 < cx.cylinder_radius || R2 < cx.sphere_radius) ? :scatter : :reflect), R1, R2, spread
    end
    R1, R2, sp, res = mc[(i-1)÷nl+1]
    R1, R2 = 1.25 * R1, 1.25 * R2
    trusted = res <= _MESH_CURVATURE_RESIDUAL && sp <= _MESH_CURVATURE_SPREAD
    return ((trusted && (R1 < cx.cylinder_radius || R2 < cx.sphere_radius)) ? :scatter : :reflect), R1, R2, sp
end

# Limits on the fit residual and on the spread of the normals (degrees) below which a
# mesh estimate of the radii counts.
const _MESH_CURVATURE_RESIDUAL = 0.1
const _MESH_CURVATURE_SPREAD = 25.0

# Returns (kind, radius, dimple half-side), with kind :cylinder, :sphere, :dimpled, :flat,
# :unknown (radii from `_mesh_curvature`), or :cone or :torus (radii from `_surface_radii`).
function _talis_curvature(t::Target, ci, overrides)
    c = t.components[ci]
    if haskey(overrides, c.name) || haskey(overrides, c.part)
        v = haskey(overrides, c.name) ? overrides[c.name] : overrides[c.part]
        v === :flat && return (:flat, Inf, NaN)
        kind = Symbol(v[1])
        kind in (:cylinder, :sphere, :dimpled) ||
            throw(ArgumentError("curvature kind must be :cylinder, :sphere, :dimpled or :flat"))
        return (kind, Float64(v[2]), length(v) >= 3 ? Float64(v[3]) : NaN)
    end
    c.surface isa SphereSurface && return (:sphere, c.surface.radius, NaN)
    c.surface isa CylinderSurface && return (:cylinder, c.surface.radius, NaN)
    c.surface isa DiskSurface && return (:flat, Inf, NaN)
    c.surface isa ConeSurface && return (:cone, NaN, NaN)
    c.surface isa TorusSurface && return (:torus, NaN, NaN)
    ns = [t.tri_normal[k] for k in eachindex(t.tri_area) if t.tri_comp[k] == ci && t.tri_area[k] > 0]
    all(n -> abs(dot(n, ns[1])) > 1 - 1e-9, ns) && return (:flat, Inf, NaN)
    return (:unknown, NaN, NaN)
end

# Transverse signed distance (negative inside) from `x` to the union of the 1ω footprints
# of a beam, at nominal focus and as the minimum over focus offsets 0 and ±focus.
function _omega1_distance(x, bl, cx)
    dn, dw = Inf, Inf
    for (kb, fb, ô) in bl
        d0 = _omega1_square(x - fb, kb, ô, cx.c0, cx.h0, cx.f1, cx.L)
        dn = min(dn, d0)
        dw = min(dw, d0)
        if cx.focus > 0
            dw = min(dw, _omega1_square(x - fb + cx.focus * kb, kb, ô, cx.c0, cx.h0, cx.f1, cx.L),
                     _omega1_square(x - fb - cx.focus * kb, kb, ô, cx.c0, cx.h0, cx.f1, cx.L))
        end
    end
    return dn, dw
end

# The Lipschitz constant of the 1ω footprint distance, 1 plus the taper of the center
# (c0/L) and half-side (h0/f1) of the footprint along the beam.
_omega1_lipschitz(c0, h0, f1, L) = 1 + hypot(c0 / L + h0 / f1, h0 / f1)

# Footprint square at axial position s = d·k from the 3ω focus. Its center is c0(1 + s/L)
# along ô, the chief ray from the lens center, and its half-side is h0|1 − s/f1|, from the
# cone to the 1ω focus.
function _omega1_square(d, k, ô, c0, h0, f1, L)
    s = dot(d, k)
    c = c0 * (1 + s / L)
    h = h0 * abs(1 - s / f1)
    p = dot(d, ô) - c
    q = dot(d, cross(k, ô))
    ox, oy = abs(p) - h, abs(q) - h
    return (ox > 0 || oy > 0) ? hypot(max(ox, 0.0), max(oy, 0.0)) : max(ox, oy)
end

# `_omega1_square` with its witness: the distance, the axial station s, the offset
# (Δp, Δq) along ô and k × ô from the point to the closest point of the square boundary,
# and the sides of the square that hold that point. The bits of the mask are 0x01 for the
# side toward the 3ω axis (−ô), 0x02 for the far side (+ô), 0x04 for −k × ô and 0x08 for
# +k × ô. A corner sets two bits, as does a point inside the square on a diagonal.
function _omega1_witness(d, k, ô, c0, h0, f1, L)
    s = dot(d, k)
    c = c0 * (1 + s / L)
    h = h0 * abs(1 - s / f1)
    p = dot(d, ô) - c
    q = dot(d, cross(k, ô))
    ox, oy = abs(p) - h, abs(q) - h
    bp = p >= 0 ? 0x02 : 0x01
    bq = q >= 0 ? 0x08 : 0x04
    if ox > 0 || oy > 0
        dist = hypot(max(ox, 0.0), max(oy, 0.0))
        mask = (ox >= 0 ? bp : 0x00) | (oy >= 0 ? bq : 0x00)
        return dist, s, clamp(p, -h, h) - p, clamp(q, -h, h) - q, mask
    end
    mask = (ox >= oy ? bp : 0x00) | (oy >= ox ? bq : 0x00)
    return ox >= oy ? (ox, s, copysign(h, p) - p, 0.0, mask) : (oy, s, 0.0, copysign(h, q) - q, mask)
end

"""
    Omega1Approach

The closest approach of the 1ω footprint of a beam to a hardware point, as
[`omega1_closest_approach`](@ref) returns it: the geometry behind one number of the 1ω
rules of [`talis_check`](@ref).

- `beam`, `component` — the beam name and the part or component of the hardware point
- `beamlet_index`, `beamlet` — the beamlet whose footprint is closest, as its number
  among the beamlets of the beam and as `(k, f3w, ô)`: its direction, its nominal 3ω
  best focus and the direction in which its footprint is displaced
- `focus_shift` — the displacement (mm along `k`) of the focal plane, out of 0 and
  `±focus`, at which the footprint is closest
- `x_hw` — the hardware point
- `x_beam` — the closest point of the boundary of the footprint square, in the plane
  through `x_hw` orthogonal to `k` (the rules measure in that plane)
- `s_star` — the axial station of that plane, in mm past the shifted 3ω focus
- `face_mask` — the sides of the square pyramid that hold `x_beam`: bit `0x01` for the
  side toward the 3ω axis, `0x02` for the far side (`+ô`), `0x04` for the side along
  `-k × ô` and `0x08` for the side along `+k × ô`. A corner sets two bits.
- `distance` — the signed transverse distance from `x_hw` to `x_beam` (mm, negative
  inside the footprint) at `focus_shift`
- `nominal` — the same distance from the footprint of that beamlet with the focal
  plane at its nominal position
- `worst` — `distance` less the transverse errors and the rule tolerance of `budget`
- `required`, `margin` — the required clearance and `worst - required`
- `budget` — the terms between `nominal` and `margin`, each in mm: `pointing`,
  `target_offset`, `hardware_tolerance`, `rotation_tol` (the rule tolerance at `x_hw`),
  `focus_delta` (`nominal - distance`) and `required`, so that
  `margin = nominal - pointing - target_offset - hardware_tolerance - rotation_tol -
  focus_delta - required`
- `contour` — the footprint as `(; c0, h0, focal_shift, lens_distance)` (see
  [`omega1_contour`](@ref)), from which the square at any station follows
"""
struct Omega1Approach
    beam::String
    component::Symbol
    beamlet_index::Int
    beamlet::NTuple{3,Vec3}
    focus_shift::Float64
    x_hw::Vec3
    x_beam::Vec3
    s_star::Float64
    face_mask::UInt8
    distance::Float64
    nominal::Float64
    worst::Float64
    required::Float64
    margin::Float64
    budget::@NamedTuple{pointing::Float64, target_offset::Float64, hardware_tolerance::Float64,
                        rotation_tol::Float64, focus_delta::Float64, required::Float64}
    contour::@NamedTuple{c0::Float64, h0::Float64, focal_shift::Float64, lens_distance::Float64}
end

# The settings of the 1ω clearance that `omega1_closest_approach` shares with
# `talis_check`.
function _omega1_setup(tolerance, pivot, pointing, target_offset, hardware_tolerance, focus,
                       omega1_edge, omega1_size, omega1_focus, lens_distance)
    o = omega1_contour(; edge = omega1_edge, size = omega1_size, focus = omega1_focus)
    L = Float64(lens_distance)
    return (; tol = _talis_tol(tolerance, pivot), pointing = Float64(pointing),
            target_offset = Float64(target_offset), hardware_tolerance = Float64(hardware_tolerance),
            R = pointing + target_offset + hardware_tolerance, focus = Float64(focus),
            o.c0, o.h0, f1 = o.focal_shift, L, lip1 = _omega1_lipschitz(o.c0, o.h0, o.focal_shift, L))
end

function _omega1_named_beam(beams, id)
    i = id isa Integer ? Int(id) : findfirst(b -> b.name == String(id), beams)
    (i === nothing || !checkbounds(Bool, beams, i)) && throw(ArgumentError("no beam $(repr(id))"))
    return _value_beam(beams[i])
end

function _omega1_named_beamlets(b::Beam, offsets)
    bl = _omega1_beamlets(b, offsets, _is_nif(b) ? _nif_beam_table() : nothing)
    bl === nothing &&
        throw(ArgumentError("beam $(b.name) has no 1ω footprint (no NIF quad name or `offsets` entry)"))
    return bl
end

# The approach at the hardware point `x`, over the beamlets `bl` and the focus offsets.
function _omega1_approach(b::Beam, bl, component::Symbol, x::Vec3, cx, required::Float64)
    best, jb, δb = (Inf, 0.0, 0.0, 0.0, 0x00), 0, 0.0
    for (j, (kb, fb, ô)) in enumerate(bl), δ in (cx.focus > 0 ? (0.0, -cx.focus, cx.focus) : (0.0,))
        w = _omega1_witness(x - fb - δ * kb, kb, ô, cx.c0, cx.h0, cx.f1, cx.L)
        w[1] < best[1] && ((best, jb, δb) = (w, j, δ))
    end
    kb, fb, ô = bl[jb]
    dist, s, Δp, Δq, mask = best
    nominal = _omega1_square(x - fb, kb, ô, cx.c0, cx.h0, cx.f1, cx.L)
    rot = cx.tol(x, Vec3(b.aim))
    worst = dist - cx.R - rot
    budget = (; cx.pointing, cx.target_offset, cx.hardware_tolerance, rotation_tol = rot,
              focus_delta = nominal - dist, required)
    return Omega1Approach(b.name, component, jb, (kb, fb, ô), δb, x, x + Δp * ô + Δq * cross(kb, ô), s,
                          mask, dist, nominal, worst, required, worst - required, budget,
                          (; cx.c0, cx.h0, focal_shift = cx.f1, lens_distance = cx.L))
end

"""
    omega1_closest_approach(target, beams, beam => part; tolerance = 0.5, pivot = :tcc,
                            pointing = 0.05, target_offset = 0.025, hardware_tolerance = 0.0,
                            focus = 1.0, omega1_margin = 0.2, omega1_edge = 4.45,
                            omega1_size = nothing, omega1_focus = nothing,
                            lens_distance = NIF_LENS_DISTANCE, offsets = nothing,
                            resolution = 1e-4, rtol = 0.01) -> Omega1Approach
    omega1_closest_approach(target, beams, beam => point; ...)
    omega1_closest_approach(target, beams, row::TaLISRow; ...)
    omega1_closest_approach(target, beams, report::TaLISReport, row; offsets = nothing)

The closest approach of the 1ω footprint of a beam to hardware, with its witness geometry
(see [`Omega1Approach`](@ref)). `beam` is a beam name or an index into `beams`. The
keywords and the footprint are those of the 1ω clearance of [`talis_check`](@ref).

With a part or component name, the hardware point is the point of the part with the
smallest worst-case clearance, found by the branch and bound of `talis_check` to the same
`resolution` and `rtol`. With a point (three coordinates, mm), the approach is evaluated
at that point, and its `component` is `:point`. With a `:clearance_1w`, `:management_1w`
or `:diagnostic_1w` row of a report, the hardware point is the `point` of the row. The
form that takes the report uses the settings of that report; `offsets` is not among them
and is passed again.

The approach is exact at its hardware point: `distance` is the transverse distance
between `x_hw` and `x_beam`. The `worst` and `nominal` of a row are lower bounds of minima
over the component, so the `worst` of the approach at the point of a row lies at or above
the `worst` of the row and within its accuracy of it, and its `nominal` (one beamlet, at
one point) lies at or above the `nominal` of the row.
"""
function omega1_closest_approach(t::Target, beams::AbstractVector{<:Beam}, pair::Pair;
                                 tolerance = 0.5, pivot = :tcc, pointing = 0.05, target_offset = 0.025,
                                 hardware_tolerance = 0.0, focus = 1.0, omega1_margin = 0.2,
                                 omega1_edge = 4.45, omega1_size = nothing, omega1_focus = nothing,
                                 lens_distance = NIF_LENS_DISTANCE, offsets = nothing,
                                 resolution = 1e-4, rtol = 0.01, required = nothing)
    cx = _omega1_setup(tolerance, pivot, pointing, target_offset, hardware_tolerance, focus,
                       omega1_edge, omega1_size, omega1_focus, lens_distance)
    b = _omega1_named_beam(beams, first(pair))
    bl = _omega1_named_beamlets(b, offsets)
    req = Float64(required === nothing ? omega1_margin : required)
    what = last(pair)
    what isa Union{Symbol,AbstractString} ||
        return _omega1_approach(b, bl, :point, Vec3(what...), cx, req)
    name = Symbol(what)
    tris = NTuple{3,Vec3}[]
    for ci in part_pieces(t, name)
        append!(tris, _talis_samples(t, ci, 0).tris)
    end
    isempty(tris) && throw(ArgumentError("part $(repr(name)) has no triangles"))
    aim = Vec3(b.aim)
    fw = x -> _omega1_distance(x, bl, cx)[2] - cx.R - cx.tol(x, aim)
    _, x = _lipschitz_minimum(fw, tris, cx.lip1 + cx.tol.tanα,
                              (Float64(resolution), Float64(rtol), req), aim)
    return _omega1_approach(b, bl, name, x, cx, req)
end

function omega1_closest_approach(t::Target, beams::AbstractVector{<:Beam}, row::TaLISRow; kw...)
    row.rule in (:clearance_1w, :management_1w, :diagnostic_1w) ||
        throw(ArgumentError("rule $(repr(row.rule)) is not a 1ω rule"))
    b = _omega1_named_beam(beams, row.beam)
    a = row.rule === :diagnostic_1w ?
        omega1_closest_approach(t, beams, b.name => row.point; kw..., required = row.required) :
        omega1_closest_approach(t, beams, b.name => row.point; kw...)
    return Omega1Approach(a.beam, row.component, a.beamlet_index, a.beamlet, a.focus_shift, a.x_hw,
                          a.x_beam, a.s_star, a.face_mask, a.distance, a.nominal, a.worst, a.required,
                          a.margin, a.budget, a.contour)
end

function omega1_closest_approach(t::Target, beams::AbstractVector{<:Beam}, r::TaLISReport, row::TaLISRow;
                                 offsets = nothing)
    s = r.settings
    return omega1_closest_approach(t, beams, row; s.tolerance,
                                   pivot = s.pivot === :point ? s.pivot_point : s.pivot, s.pointing,
                                   s.target_offset, s.hardware_tolerance, s.focus, s.omega1_margin,
                                   s.omega1_edge, s.omega1_size, s.omega1_focus, s.lens_distance,
                                   offsets)
end

"""
    Omega1Disposition

The disposition of hardware toward the 1ω light of a set of beams, per mesh triangle of
the target, as [`omega1_disposition`](@ref) returns it. Each vector has one entry per
triangle of the target (`target.tri_comp` gives its component).

- `status` — `:unchecked` (a triangle outside the checked components), or, from least to
  most severe:
  - `:clear` — the worst-case clearance exceeds the required one by more than the tight
    threshold
  - `:tight` — clears, by the tight threshold or less
  - `:curved_safe` — inside the worst-case footprint, and scatters 1ω by its curvature
  - `:specular_clear` — inside the footprint, flat or gently curved, and its specular
    reflection stays the required angle from every open port
  - `:review` — inside the footprint and not decidable: a reflection toward a port from a
    mesh of unknown curvature, or a dimpled surface without a valid spacing
  - `:fail` — inside the footprint with a reflection within the required angle of an open
    port, or a diagnostic component closer to the 1ω light than its rule allows
- `margin` — the smallest worst-case clearance margin (mm) over the beams, sampled at the
  vertices of the triangle and the centroids of its sub-triangles; `Inf` beyond
  `margin_range`
- `beam`, `beamlet` — the beam (index into the beams) and its beamlet with that margin; 0
  where the margin is `Inf`
- `R1`, `R2` — the unsigned principal radii (mm) on which a triangle in the footprint was
  judged (`Inf` for none); `NaN` outside the footprint
- `port`, `port_angle` — for a triangle judged on reflection, the open port (index into
  `ports`) nearest its reflection and the worst-case angle to it (degrees, less the
  normal spread and twice the rule tolerance); 0 and `NaN` otherwise
- `shadowed` — with `shadows = true`, whether other hardware lies upstream of the
  centroid of a triangle in the footprint along the beam (the rules ignore shadowing; the
  flag separates lit hardware from shadowed hardware in a review)

The directions to the open ports are in `ports`, and the thresholds and error set used
in `settings`.
"""
struct Omega1Disposition
    status::Vector{Symbol}
    margin::Vector{Float64}
    beam::Vector{Int}
    beamlet::Vector{Int}
    R1::Vector{Float64}
    R2::Vector{Float64}
    port::Vector{Int}
    port_angle::Vector{Float64}
    shadowed::Vector{Bool}
    ports::Vector{Vec3}
    settings::NamedTuple
end

const _DISPOSITION_RANK = (clear = 0, tight = 1, curved_safe = 2, specular_clear = 3, review = 4, fail = 5)

function Base.show(io::IO, d::Omega1Disposition)
    print(io, "Omega1Disposition(", count(!=(:unchecked), d.status), " triangles checked")
    for st in (:fail, :review, :specular_clear, :curved_safe, :tight)
        n = count(==(st), d.status)
        n > 0 && print(io, ", ", n, " ", st)
    end
    print(io, ")")
end

"""
    omega1_disposition(target, beams; components = all surfaces, diagnostics = (),
                       tight_threshold = 1.0, margin_range = 5.0, shadows = false,
                       tolerance = 0.5, pivot = :tcc, pointing = 0.05, target_offset = 0.025,
                       hardware_tolerance = 0.0, focus = 1.0, omega1_margin = 0.2,
                       omega1_edge = 4.45, omega1_size = nothing, omega1_focus = nothing,
                       lens_distance = NIF_LENS_DISTANCE, offsets = nothing, open_ports = :nif,
                       reflection_angle = 20.0, cylinder_radius = 2.0, sphere_radius = 5.0,
                       curvature = Dict(), level = 2, resolution = 1e-4) -> Omega1Disposition

The 1ω rules of [`talis_check`](@ref) as a field over the hardware: for every mesh
triangle of `components` (parts or components; all surfaces by default), how it stands
toward the 1ω footprints of all `beams` together (see [`Omega1Disposition`](@ref)). The
footprint, the error set, the curvature limits and the reflection check, with their
keywords, are those of `talis_check`, whose `:clearance_1w` and `:management_1w` rows
give the same verdicts per beam and component. `diagnostics` names the components under
the diagnostic rule (10 mm from 1ω light at TCC, growing to 25 mm at 1 m), which fail
where they are closer.

Whether a triangle clears, is tight (within `tight_threshold` mm of the required
clearance) or lies in the footprint is decided on the whole triangle by the bound of
`talis_check`, to `resolution`. A triangle in the footprint takes the most severe
verdict of its `4^level` sub-triangles over the beams. `margin` is sampled and is
followed up to `margin_range` mm. With `shadows = true` the triangles in the footprint
are tested for hardware upstream. The cost grows with the number of sub-triangles in the
footprints; `level = 0` suits CAD meshes whose facets are already fine.
"""
function omega1_disposition(t::Target, beams::AbstractVector{<:Beam}; components = nothing, diagnostics = (),
                            tight_threshold = 1.0, margin_range = 5.0, shadows::Bool = false,
                            tolerance = 0.5, pivot = :tcc, pointing = 0.05, target_offset = 0.025,
                            hardware_tolerance = 0.0, focus = 1.0, omega1_margin = 0.2,
                            omega1_edge = 4.45, omega1_size = nothing, omega1_focus = nothing,
                            lens_distance = NIF_LENS_DISTANCE, offsets = nothing, open_ports = :nif,
                            reflection_angle = 20.0, cylinder_radius = 2.0, sphere_radius = 5.0,
                            curvature = Dict{Symbol,Any}(), level::Integer = 2, resolution = 1e-4)
    resolution > 0 || throw(ArgumentError("resolution must be > 0"))
    (tight_threshold >= 0 && margin_range >= tight_threshold) ||
        throw(ArgumentError("need 0 ≤ tight_threshold ≤ margin_range"))
    vb = [_value_beam(b) for b in beams]
    table = any(_is_nif, vb) ? _nif_beam_table() : nothing
    bls = [_omega1_beamlets(b, offsets, table) for b in vb]
    pieces(names) = Int[ci for n in names for ci in part_pieces(t, Symbol(n))]
    dg = pieces(diagnostics)
    cis = components === nothing ?
          Int[ci for (ci, c) in enumerate(t.components) if c.kind === :surface] : pieces(components)
    cis = unique(vcat(cis, dg))
    geo = Dict(ci => _talis_samples(t, ci, level) for ci in cis)
    tri_of = Dict(ci => [k for k in eachindex(t.tri_area) if t.tri_comp[k] == ci && t.tri_area[k] > 0] for ci in cis)
    ports = _talis_open_ports(open_ports, vb)
    cx = merge(_omega1_setup(tolerance, pivot, pointing, target_offset, hardware_tolerance, focus,
                             omega1_edge, omega1_size, omega1_focus, lens_distance),
               (; geo, ports, tolerance = Float64(tolerance), omega1_margin = Float64(omega1_margin),
                reflection_angle = Float64(reflection_angle), cylinder_radius = Float64(cylinder_radius),
                sphere_radius = Float64(sphere_radius), curvature,
                curv = Dict{Int,Vector{NTuple{4,Float64}}}(), curv_lock = ReentrantLock(),
                ε = Float64(resolution), tight = Float64(tight_threshold), range = Float64(margin_range),
                shadows, diagnostics = Set(dg), tri_of))
    nt = length(t.tri_area)
    d = Omega1Disposition(fill(:unchecked, nt), fill(Inf, nt), zeros(Int, nt), zeros(Int, nt),
                          fill(NaN, nt), fill(NaN, nt), zeros(Int, nt), fill(NaN, nt), fill(false, nt), ports,
                          (; tolerance = Float64(tolerance), pivot = cx.tol.pivot, pivot_point = cx.tol.point,
                           pointing = Float64(pointing), target_offset = Float64(target_offset),
                           hardware_tolerance = Float64(hardware_tolerance), focus = Float64(focus),
                           omega1_margin = Float64(omega1_margin), omega1_edge = Float64(omega1_edge),
                           omega1_size = 2 * cx.h0, omega1_focus = cx.f1, lens_distance = cx.L,
                           reflection_angle = Float64(reflection_angle),
                           cylinder_radius = Float64(cylinder_radius), sphere_radius = Float64(sphere_radius),
                           tight_threshold = Float64(tight_threshold), margin_range = Float64(margin_range),
                           level = Int(level), resolution = Float64(resolution)))
    for ci in cis
        d.status[tri_of[ci]] .= :clear
    end
    found = [Tuple{Int,Symbol,Float64,Int,Float64,Float64,Int,Float64,Bool}[] for _ in vb]
    Threads.@threads for bi in eachindex(vb)
        bls[bi] === nothing || _disposition_beam!(found[bi], t, vb[bi], bls[bi], cis, cx)
    end
    for (bi, fs) in enumerate(found), (k, st, m, jb, R1, R2, pj, α, sh) in fs
        if m < d.margin[k]
            d.margin[k], d.beam[k], d.beamlet[k] = m, bi, jb
        end
        r0, r1 = _DISPOSITION_RANK[d.status[k]], _DISPOSITION_RANK[st]
        if r1 > r0 || (r1 == r0 && pj != 0 && !(α >= d.port_angle[k]))
            d.status[k] = st
            r1 >= 2 && ((d.R1[k], d.R2[k], d.shadowed[k]) = (R1, R2, sh))
            d.port[k], d.port_angle[k] = pj, pj == 0 ? NaN : α
        end
    end
    return d
end

# The triangles that one beam brings within `margin_range` of its 1ω footprint, each as
# (triangle, status, margin, beamlet, R1, R2, port, port angle, shadowed).
function _disposition_beam!(out, t::Target, b::Beam, bl, cis, cx)
    aim = Vec3(b.aim)
    kq = direction(b)
    lip = cx.lip1 + cx.tol.tanα
    fw = x -> _omega1_distance(x, bl, cx)[2] - cx.R - cx.tol(x, aim)
    far = 4 * (t.scale + norm(aim)) + 1
    εr = 1e3 * _ray_eps(t)
    for ci in cis
        g = cx.geo[ci]
        isempty(g.tris) && continue
        c = t.components[ci]
        diag = ci in cx.diagnostics
        # the requirement of a diagnostic grows with the distance from TCC
        req = diag ? (x -> 10 + 0.015 * norm(x)) : (x -> cx.omega1_margin)
        f = x -> fw(x) - req(x)
        lf = lip + (diag ? 0.015 : 0.0)
        nl = length(g.leaves) ÷ length(g.tris)
        kind, rc, half = diag ? (:flat, Inf, NaN) : _talis_curvature(t, ci, cx.curvature)
        mc = nothing
        for (j, tri) in enumerate(g.tris)
            m0 = (tri[1] + tri[2] + tri[3]) / 3
            ρ = maximum(x -> norm(x - m0), tri)
            f(m0) - lf * ρ > cx.range && continue
            m = min(f(tri[1]), f(tri[2]), f(tri[3]))
            ls = (j-1)*nl+1:j*nl
            for i in ls
                m = min(m, f(g.leaves[i][4]))
            end
            m > cx.range && !_talis_reaches(f, tri, lf, cx.tight, cx.ε) && continue
            k = cx.tri_of[ci][j]
            jb = _nearest_beamlet(m0, bl, cx)
            if !_talis_reaches(f, tri, lf, 0.0, cx.ε)
                st = _talis_reaches(f, tri, lf, cx.tight, cx.ε) ? :tight : :clear
                push!(out, (k, st, max(m, 0.0), jb, NaN, NaN, 0, NaN, false))
                continue
            end
            sh = cx.shadows && !visible(t, m0 - far * kq, m0 - εr * kq)
            if diag
                push!(out, (k, :fail, min(m, 0.0), jb, NaN, NaN, 0, NaN, sh))
                continue
            end
            kind === :unknown && mc === nothing && (mc = _component_curvature(t, ci, cx))
            st, R1, R2, pj, α = :clear, NaN, NaN, 0, NaN
            for i in ls
                leaf = g.leaves[i]
                lt = (leaf[1], leaf[2], leaf[3])
                _talis_reaches(f, lt, lf, 0.0, cx.ε) || continue
                class, r1, r2, spread = _leaf_curvature(c.surface, kind, rc, half, mc, nl, i, lt, g.spreads[i], cx)
                s1, p1, α1 = class === :scatter ? :curved_safe : :review, 0, NaN
                if class === :reflect
                    n = g.normals[i]
                    α1 = Inf
                    for (kb, _, _) in bl
                        r = kb - 2 * dot(kb, n) * n
                        raw, pp = findmin(p -> acosd(clamp(dot(r, p), -1.0, 1.0)), cx.ports)
                        raw - 2 * spread < α1 && ((α1, p1) = (raw - 2 * spread, pp))
                    end
                    α1 -= 2 * cx.tolerance
                    s1 = α1 > cx.reflection_angle ? :specular_clear : kind === :unknown ? :review : :fail
                end
                a, a1 = _DISPOSITION_RANK[st], _DISPOSITION_RANK[s1]
                if a1 > a || (a1 == a && p1 != 0 && !(α1 >= α))
                    st, R1, R2, pj, α = s1, r1, r2, p1, α1
                end
            end
            # the bound may reach a triangle none of whose sub-triangles it reaches
            st === :clear && (st = :tight)
            push!(out, (k, st, min(m, st === :tight ? cx.tight : 0.0), jb, R1, R2, pj, α, sh))
        end
    end
    return out
end

# The beamlet of `bl` whose worst-case 1ω footprint is nearest the point `x`.
function _nearest_beamlet(x, bl, cx)
    best, jb = Inf, 0
    for (j, (kb, fb, ô)) in enumerate(bl), δ in (cx.focus > 0 ? (0.0, -cx.focus, cx.focus) : (0.0,))
        v = _omega1_square(x - fb - δ * kb, kb, ô, cx.c0, cx.h0, cx.f1, cx.L)
        v < best && ((best, jb) = (v, j))
    end
    return jb
end

function _talis_open_ports(spec, beams)
    spec === :nif && return [port_vector(θ, φ) for (_, θ, φ) in NIF_PORTS]
    spec === :beams && return unique([_is_nif(b) ? _nif_quad_port(b.quad) : b.port for b in beams])
    spec isa Symbol && throw(ArgumentError("open_ports must be :nif, :beams or a list of ports"))
    return [x isa Beam ? x.port : normalize(Vec3(x...)) for x in spec]
end

# TAS stay-in zone (NIF beams), as the extent of each component above the lowest
# upper-beam aim point and below the highest lower-beam aim point.
function _talis_stay_in(t::Target, beams, surf, limit, htol)
    rows = TaLISRow[]
    nif = [b for b in beams if _is_nif(b)]
    isempty(nif) && return rows
    for (grp, sgn) in (("upper beams", 1), ("lower beams", -1))
        bs = [b for b in nif if sgn * b.port[3] > 0]
        isempty(bs) && continue
        zref = sgn > 0 ? minimum(b.aim[3] for b in bs) : maximum(b.aim[3] for b in bs)
        for (name, cis) in surf
            ks = [k for k in eachindex(t.tri_area) if t.tri_comp[k] in cis]
            isempty(ks) && continue
            pts = [v for k in ks for v in (t.tri_v0[k], t.tri_v0[k] + t.tri_e1[k], t.tri_v0[k] + t.tri_e2[k])]
            e, i = findmax(sgn * (p[3] - zref) for p in pts)
            m = limit - (e + htol)
            push!(rows, TaLISRow(grp, :tas_stay_in, name, m >= 0 ? :pass : :fail, m,
                                 :mm, e + htol, e, limit, pts[i],
                                 m >= 0 ? "" : "outside the TAS stay-in zone: may need a different alignment sequence or beam selection"))
        end
    end
    return rows
end

# Mesh triangles, sub-triangles (a, b, c, centroid) and their normals (exact on sphere,
# cylinder, cone and torus surfaces) of a component. `spreads` holds, per sub-triangle of
# a curved surface, the largest angle (degrees) between the surface normal at its centroid
# and at its vertices, which bounds the normals of the surface over it; it is 0 on a flat
# mesh.
function _talis_samples(t::Target, ci, level)
    c = t.components[ci]
    curved = c.surface isa Union{SphereSurface,CylinderSurface,ConeSurface,TorusSurface}
    tris, leaves, normals, buf = NTuple{3,Vec3}[], NTuple{4,Vec3}[], Vec3[], Vec3[]
    spreads = Float64[]
    for k in eachindex(t.tri_area)
        (t.tri_comp[k] == ci && t.tri_area[k] > 0) || continue
        v0, e1, e2 = t.tri_v0[k], t.tri_e1[k], t.tri_e2[k]
        push!(tris, (v0, v0 + e1, v0 + e2))
        empty!(buf)
        _subtri_vertices!(buf, v0, v0 + e1, v0 + e2, level)
        for j in 1:3:length(buf)
            a, bb, d = buf[j], buf[j+1], buf[j+2]
            m = (a + bb + d) / 3
            push!(leaves, (a, bb, d, m))
            if curved
                n = last(surface_point(c.surface, m))
                push!(normals, n)
                push!(spreads, maximum(x -> acosd(clamp(dot(n, last(surface_point(c.surface, x))), -1.0, 1.0)),
                                       (a, bb, d)))
            else
                push!(normals, t.tri_normal[k])
                push!(spreads, 0.0)
            end
        end
    end
    return (; tris, leaves, normals, spreads)
end

# The mesh triangles of the components `cis`.
_talis_tris(cx, cis) = length(cis) == 1 ? cx.geo[only(cis)].tris :
                       reduce(vcat, (cx.geo[ci].tris for ci in cis))

# Whether `f` (Lipschitz constant `lip`) is ≤ `thr` somewhere on the simplex `s`. A
# simplex still undecided at lip·ρ ≤ `ε` counts as reached (conservative).
function _talis_reaches(f, s, lip, thr, ε)
    m = sum(s) / length(s)
    fm = f(m)
    fm <= thr && return true
    ρ = maximum(x -> norm(x - m), s)
    fm - lip * ρ > thr && return false
    lip * ρ <= ε && return true
    return any(c -> _talis_reaches(f, c, lip, thr, ε), _simplex_split(s))
end

# The samples of several components (the pieces of a part) as one.
_merge_samples(gs) = (; tris = reduce(vcat, (g.tris for g in gs)),
                      leaves = reduce(vcat, (g.leaves for g in gs)),
                      normals = reduce(vcat, (g.normals for g in gs)),
                      spreads = reduce(vcat, (g.spreads for g in gs)))

# The free edges (one opaque triangle) of the given components, as segments; vertices
# are merged by position across components.
function _free_edges(t::Target, comps)
    q = 1e-7 * max(t.scale, 1.0)
    key(v) = (round(Int, v[1] / q), round(Int, v[2] / q), round(Int, v[3] / q))
    seen = Dict{NTuple{6,Int},Tuple{Int,Int,Vec3,Vec3}}()
    for k in eachindex(t.tri_area)
        (t.tri_opaque[k] && t.tri_area[k] > 0) || continue
        vs = (t.tri_v0[k], t.tri_v0[k] + t.tri_e1[k], t.tri_v0[k] + t.tri_e2[k])
        for (i, j) in ((1, 2), (2, 3), (3, 1))
            a, b = key(vs[i]), key(vs[j])
            e = a < b ? (a..., b...) : (b..., a...)
            n, ci, va, vb = get(seen, e, (0, Int(t.tri_comp[k]), vs[i], vs[j]))
            seen[e] = (n + 1, ci, va, vb)
        end
    end
    segs = Dict(ci => NTuple{2,Vec3}[] for ci in comps)
    for (n, ci, a, b) in values(seen)
        (n == 1 && haskey(segs, ci)) || continue
        push!(segs[ci], (a, b))
    end
    return segs
end

const _TALIS_RULE_ORDER = (:capture_3w, :incidence_3w, :clearance_3w, :clearance_1w,
                           :management_1w, :diagnostic_1w, :tas_stay_in)

"""
    talis_requirement(rule, settings) -> String

Human-readable requirement of a `talis_check` rule (`:capture_3w`, `:incidence_3w`,
`:clearance_3w`, `:clearance_1w`, `:management_1w`, `:diagnostic_1w`, `:tas_stay_in`)
under the thresholds in `settings` (a `TaLISReport`'s `settings`).
"""
function talis_requirement(rule, s)
    rule === :capture_3w &&
        return @sprintf("%g%% contour ≥ %g mm inside the substrate edge", 100s.capture_contour,
                        s.capture_margin)
    rule === :incidence_3w && return @sprintf("incidence < %g° to the surface normal", s.max_incidence)
    rule === :clearance_3w &&
        return @sprintf("%g%% contour > %g mm from the component", 100s.clearance_contour,
                        s.clearance_margin)
    rule === :clearance_1w &&
        return @sprintf("1ω 10%% footprint (%g mm off the 3ω axis) > %g mm from the component",
                        s.omega1_edge, s.omega1_margin)
    rule === :management_1w &&
        return @sprintf("curvature radius < %g mm (cylinder) or < %g mm (sphere), or specular 1ω > %g° from open ports",
                        s.cylinder_radius, s.sphere_radius, s.reflection_angle)
    rule === :diagnostic_1w && return "diagnostics ≥ 10 mm (at TCC) to 25 mm (at 1 m) from 1ω light"
    rule === :tas_stay_in && return @sprintf("≤ %g mm beyond the extreme aim points (TAS)", s.stay_in)
    return string(rule)
end

_talis_rank(x::TaLISRow) = x.status === :fail ? 0 : x.status === :review ? 1 : 2

"""
    violations(report::TaLISReport; status = (:fail, :review)) -> Vector{TaLISRow}

Rows of a TaLIS report whose status is in `status`, sorted with failures first and then
by margin.
"""
violations(r::TaLISReport; status = (:fail, :review)) =
    sort(filter(x -> x.status in status, r.rows); by = x -> (_talis_rank(x), x.margin))

function Base.show(io::IO, ::MIME"text/plain", r::TaLISReport)
    s = r.settings
    nb = length(unique(x.beam for x in r.rows if x.rule !== :tas_stay_in))
    println(io, "TaLISReport: ", nb, " beams, ", length(r.rows), " checks against ", TALIS_SOURCE)
    @printf(io, "  rule tolerance %.3g° applied as a lateral allowance r·tan(%.3g°), %s\n",
            s.tolerance, s.tolerance, s.pivot_text)
    println(io, "    (component rotated off nominal by insertion or assembly; the 70° incidence limit excluded)")
    @printf(io, "  error set (added linearly): transverse %.3g mm (pointing %.3g + target %.3g + hardware %.3g), focus ±%.3g mm, SSD %.3g mm\n",
            s.transverse, s.pointing, s.target_offset, s.hardware_tolerance, s.focus, s.ssd)
    println(io, "  rule           requirement                                           pass review fail   worst margin  beam / component")
    for rule in _TALIS_RULE_ORDER
        rs = [x for x in r.rows if x.rule === rule]
        isempty(rs) && continue
        w = rs[argmin([x.margin for x in rs])]
        @printf(io, "  %-14s %-53s %4d %6d %4d  %+9.3f %-3s  %s / %s\n", rule,
                talis_requirement(rule, s), count(x -> x.status === :pass, rs),
                count(x -> x.status === :review, rs), count(x -> x.status === :fail, rs),
                w.margin, w.unit, w.beam, w.component)
    end
    bad = violations(r)
    if !isempty(bad)
        println(io, "  fail and review rows, worst first:")
        for x in first(bad, 20)
            @printf(io, "    %-6s %-12s %-14s %-14s %+9.3f %-3s at (%.3f, %.3f, %.3f)  %s\n", x.status,
                    x.beam, x.rule, x.component, x.margin, x.unit, x.point..., x.note)
        end
        length(bad) > 20 && println(io, "    … and ", length(bad) - 20, " more (see `violations`)")
    end
    for n in r.notes
        println(io, "  not evaluated: ", n)
    end
    println(io, "  out of scope (see TALIS_OUT_OF_SCOPE): positioner clearance, targets to TAS, ",
            "TAS withdrawal envelope, diagnostics to shroud, ARC footprints, blow-by assessment")
end

"""
    write_talis_csv(path, report)

Write every row of a TaLIS report as CSV, with margins in the `unit` of each row. A
header of `#` comment lines makes the file stand alone. The header contains the source of
the rules, every setting in `settings`, the requirement of each rule, the checks not
evaluated (`notes`) and the requirements out of scope (`TALIS_OUT_OF_SCOPE`). The settings
are the thresholds, the error set, the rule tolerance and pivot, the capture, LEH and
diagnostic components, and the resolution. CSV readers skip the header with a comment
option, e.g. `CSV.File(path; comment = "#")`.
"""
function write_talis_csv(path::AbstractString, r::TaLISReport)
    open(path, "w") do io
        println(io, "# TaLIS checks (Briareus talis_check), rules from ", TALIS_SOURCE)
        println(io, "# partial assessment: the requirements out of scope below are not checked")
        for (k, v) in pairs(r.settings)
            println(io, "# setting ", k, " = ", v isa AbstractVector ? join(v, " ") : v)
        end
        for rule in _TALIS_RULE_ORDER
            any(x -> x.rule === rule, r.rows) &&
                println(io, "# rule ", rule, ": ", talis_requirement(rule, r.settings))
        end
        for n in r.notes
            println(io, "# not evaluated: ", n)
        end
        for n in TALIS_OUT_OF_SCOPE
            println(io, "# out of scope: ", n)
        end
        println(io, "beam,rule,component,status,margin,unit,worst,nominal,required,x_mm,y_mm,z_mm,note")
        for x in r.rows
            println(io, join((x.beam, x.rule, x.component, x.status,
                              (@sprintf("%.6g", v) for v in (x.margin,))..., x.unit,
                              (@sprintf("%.6g", v) for v in (x.worst, x.nominal, x.required, x.point...))...,
                              "\"" * replace(x.note, "\"" => "\"\"") * "\""), ","))
        end
    end
    return path
end
