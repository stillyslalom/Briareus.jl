# NIF TaLIS target-configuration checks: 3ω capture and clearance, 1ω clearance and
# unconverted-light management, and the geometric hardware rules Briareus can
# represent, judged on the worst case over pointing, positioning and focus errors.

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
    "blow-by assessment by the backscatter working group, and positioner/alignment-sequence reviews: facility processes",
]

const _TALIS_BLOWBY = "TaLIS: configuration must be assessed for blow-by by the backscatter working group"

"""
    TaLISRow

One check of `talis_check`: a (beam, rule, component) triple.

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
    TaLISReport

Result of `talis_check`: `rows` (one `TaLISRow` per beam, rule and component), the
`settings` used (thresholds, error set, tolerance interpretation) and `notes` on checks
that were not evaluated. See `violations`, `write_talis_csv`.
"""
struct TaLISReport
    rows::Vector{TaLISRow}
    settings::NamedTuple
    notes::Vector{String}
end

# The rule tolerance: a lateral allowance tan(α)·r, r the distance of the point from
# the pivot.
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
    talis_check(target, beams; capture, hardware = all surfaces, omega1_components = all
                surfaces, diagnostics = (), tolerance = 0.5, pivot = :tcc,
                pointing = 0.05, target_offset = 0.025, hardware_tolerance = 0.0,
                focus = 1.0, ssd = 0.0, offsets = nothing, open_ports = :nif,
                curvature = Dict(), level = 2, edge_step = 0.02, ...) -> TaLISReport

NIF TaLIS target-configuration checks from the "TCC Configuration Requirements"
(NIF User Forum, 27 Oct 2018, LMF Shot RI Resources). Every distance is measured
orthogonal to the beam, as the document states, and judged on the worst case over the
error set and the rule tolerance, never at nominal.

`capture` lists the components that beams are meant to land on (backlighter foils, capture
surfaces, a hohlraum wall). Each beam's *landing component* is the first of them its
nominal axis meets. Here and in `hardware`, `omega1_components`, `diagnostics` and
`curvature`, a component may be named by its part: a part split for computation
(`curved_components`) is checked as one, its rows under the part's name.

# Rules

- `:capture_3w` — the `capture_contour` (10 %)-of-peak contour must lie ≥
  `capture_margin` (0.5 mm) inside the landing component's substrate edge: the minimum
  over the edge of the transverse distance to the beam envelope (best-focus contour grown
  by the f-cone at each point's axial position, see `BeamEnvelope`). Substrate edges are
  the free edges of the opaque mesh (edges with a single triangle; edges shared with
  another component are not edges), sampled every `edge_step` mm. A beam whose axis
  misses every capture component fails with component `:none`.
- `:incidence_3w` — incidence < `max_incidence` (70°) to the target normal, the maximum
  over the landing component's surface inside the worst-case 10 % footprint (visible
  from the beam; exact normals on `SphereSurface`/`CylinderSurface` components), plus
  the tolerance angle (a component tilted by the tolerance changes incidence by up to
  as much).
- `:clearance_3w` — the `clearance_contour` (1 %)-of-peak contour must clear each
  `hardware` component other than the landing component by > `clearance_margin`
  (0.5 mm); hardware downstream of the landing point is skipped, and everything along
  the beam is checked if the beam lands on no capture component. Capture and clearance
  failures carry the TaLIS consequence: blow-by assessment by the backscatter working
  group.
- `:clearance_1w` — the 10 % contour of the 1ω footprint must clear each
  `omega1_components` component by > `omega1_margin` (0.2 mm). Not clearing is a
  `:review` (TaLIS: evaluate unconverted-light management), decided by the next rule.
- `:management_1w` — for each surface the 1ω footprint does not clear: pass if its
  radius of curvature is < `cylinder_radius` (2 mm, cylindrical) or < `sphere_radius`
  (5 mm, spherical); otherwise the specular 1ω reflection from every uncleared facet
  must stay > `reflection_angle` (20°) from every open quad (`open_ports`), less twice
  the tolerance angle (a tilt by α turns a reflection by 2α). Curvature comes from
  `curvature[name]` — `(:cylinder, r)`, `(:sphere, r)`, `:flat`, or `(:dimpled, r)` /
  `(:dimpled, r, half_side)` for a surface covered with spherical dimples of radius of
  curvature `r` (e.g. the dimpled 1ω "sombrero" shield around an ignition hohlraum's
  LEH, modeled as a flat `annulus_mesh`), judged like a sphere against
  `sphere_radius`; with `half_side`, a dimple half-side/radius below 0.25 (NIF User
  Guide §5.6.2) makes it a `:review` — else from the component's
  `CylinderSurface`/`SphereSurface`, or a planar mesh. For other meshes the curvature is
  unknown and a failing reflection check is a `:review`. A flat component without a
  `curvature` entry gets the reflection check, which a horizontal plate usually fails
  (it mirrors each beam into the conjugate port of the same hemisphere).
- `:diagnostic_1w` — the `diagnostics` components (diagnostic hardware modeled in the
  target) must stay ≥ 10 mm from 1ω light at TCC, growing linearly to 25 mm at 1 m
  (`10 + 0.015 r` mm, r the distance from TCC; interpretation of "a cone of 10 mm at TCC
  extruded to 25 mm at 1 m").
- `:tas_stay_in` — (NIF beams) no component may extend more than `stay_in` (19 mm)
  above the lowest aim point of an upper beam or below the highest aim point of a lower
  beam (TAS alignment stay-in zone; violating it may require a different alignment
  sequence or beam selection). `hardware_tolerance` is added to the extent.

The remaining requirements of the document are listed in `TALIS_OUT_OF_SCOPE`.

# 1ω footprint

The 1ω rules apply to beams whose `quad` is a NIF quad name (or that have an `offsets`
entry). A NIF quad-level beam (port within 0.1° of its quad port, as from `nif_quads`)
is expanded into its four beamlets (LMF beam table) sharing its aim and defocus. Each
beamlet's 1ω footprint at its 3ω best-focus plane is a square of side `omega1_size`
whose nearest edge is `omega1_edge` = 4.45 mm (TaLIS: the 10 % contour) from the 3ω
axis, displaced along the beamlet's own side of its quad in φ̂ (the side of the quad
axis its port lies on; NIF User Guide Fig. 4-12 — *derived*, confidence medium). The
default side is the User Guide's 22 mm footprint grown by 4.8 − 4.45 = 0.35 mm per
side, so its near edge is the TaLIS 10 % contour (22.7 mm). Away from that plane the
footprint is the cone converging to the 1ω focus `omega1_focus` (default 460 mm, from
`nif_unconverted`) beyond the 3ω focus, with its centre on the chief ray from the lens
centre `lens_distance` (7.7 m) upstream: at axial position `s` the near edge is
`4.45 + 0.027 s` mm, i.e. closer to the 3ω axis upstream of focus (*derived*). Pass
`omega1_focus = Inf, lens_distance = Inf` for a prism (the footprint fixed in the
transverse plane). `offsets` maps beam names to a footprint direction (degrees from the
port's θ̂ toward φ̂, or a vector) and overrides the NIF convention for that beam. The
focus error moves the footprint with the 3ω focal plane; shadowing is ignored
(conservative).

# Tolerance and error set

The document adds "up to 0.5 degrees metrology/alignment/positioning tolerance" to the
3ω and 1ω clearances without saying how an angle becomes a distance. **This
interpretation is an assumption that the user must confirm with TaLIS.** By default
(`pivot = :tcc`) it is a lateral allowance `r·tan(tolerance)`, with `r` the distance of
the hardware point from TCC: the transverse displacement of the point when the target
assembly rotates by the tolerance about TCC. Alternatives: `pivot = :aim` (rotation
about the beam's aim point, so the allowance vanishes where the beam lands), `:max`
(the larger of the two, most conservative), or a point (e.g. the mount of the
component, `pivot = (0, 0, 12.0)`); `tolerance = 0` drops it.

It is combined with the error set of `clearance` by linear addition (no RSS):

    worst = nominal − (pointing + target_offset + hardware_tolerance) − focus/(2F) − r·tan(tolerance)

evaluated point by point and minimised over the component (for the 1ω footprint the
focus term is the footprint shift over `±focus` instead of the f-cone). The
`target_offset` and `hardware_tolerance` partly overlap what the tolerance is meant to
cover, so this double counts conservatively; set them to 0 to rely on the tolerance
alone. `ssd` (mm) widens the 3ω envelopes as in `clearance`.

Hardware is sampled at triangle vertices and `4^level` sub-triangles; substrate edges
every `edge_step` mm. Clearance-only far-field wings (`nif_wings`) are not used: the
TaLIS rules refer to contours of the beam spot.
"""
function talis_check(t::Target, beams::AbstractVector{<:Beam}; capture,
                     hardware = nothing, omega1_components = nothing, diagnostics = (),
                     tolerance = 0.5, pivot = :tcc,
                     pointing = 0.05, target_offset = 0.025, hardware_tolerance = 0.0,
                     focus = 1.0, ssd = 0.0,
                     capture_contour = 0.1, capture_margin = 0.5, max_incidence = 70.0,
                     clearance_contour = 0.01, clearance_margin = 0.5,
                     omega1::Bool = true, omega1_margin = 0.2, omega1_edge = 4.45,
                     omega1_size = nothing, omega1_focus = nothing, lens_distance = 7700.0,
                     offsets = nothing, open_ports = :nif, reflection_angle = 20.0,
                     cylinder_radius = 2.0, sphere_radius = 5.0, curvature = Dict{Symbol,Any}(),
                     stay_in = 19.0, level::Integer = 2, edge_step = 0.02)
    tol = _talis_tol(tolerance, pivot)
    R = pointing + target_offset + hardware_tolerance
    # checked items: name => its components (a part's pieces, or one component)
    surf = Pair{Symbol,Vector{Int}}[]
    for p in parts(t)
        cis = [ci for ci in part_pieces(t, p) if t.components[ci].kind === :surface]
        isempty(cis) || push!(surf, p => cis)
    end
    idx(names) = Pair{Symbol,Vector{Int}}[Symbol(n) => part_pieces(t, Symbol(n)) for n in names]
    cap_idx = idx(capture)
    hw_idx = hardware === nothing ? surf : idx(hardware)
    o1_idx = omega1 ? (omega1_components === nothing ? surf : idx(omega1_components)) : empty(surf)
    dg_idx = omega1 ? idx(diagnostics) : empty(surf)
    geo = Dict(ci => _talis_samples(t, ci, level)
               for ci in unique(Int[ci for g in (cap_idx, hw_idx, o1_idx, dg_idx) for p in g for ci in last(p)]))
    pe = _free_edge_points(t, unique(vcat(Int[], last.(cap_idx)...)), edge_step)
    edges = [unique(vcat(Vec3[], (pe[ci] for ci in last(g))...)) for g in cap_idx]
    cap_geo = [length(last(g)) == 1 ? geo[only(last(g))] : _merge_samples([geo[ci] for ci in last(g)])
               for g in cap_idx]
    uc = first(c for c in nif_unconverted() if c.name === :omega1)
    W0 = omega1_size === nothing ? uc.footprint + 2 * (uc.edge - omega1_edge) : Float64(omega1_size)
    f1 = omega1_focus === nothing ? uc.focal_shift : Float64(omega1_focus)
    vb = [_value_beam(b) for b in beams]
    table = any(_is_nif, vb) ? _nif_beam_table() : nothing
    bls = [omega1 ? _omega1_beamlets(b, offsets, table) : nothing for b in vb]
    ports = _talis_open_ports(open_ports, vb)
    cx = (; tol, R, focus = Float64(focus), ssd = Float64(ssd), cap_idx, hw_idx, o1_idx, dg_idx,
          geo, edges, cap_geo, capture_contour, capture_margin, max_incidence = Float64(max_incidence),
          clearance_contour, clearance_margin, tolerance = Float64(tolerance),
          omega1_margin = Float64(omega1_margin), c0 = omega1_edge + W0 / 2, h0 = W0 / 2, f1,
          L = Float64(lens_distance), ports, reflection_angle = Float64(reflection_angle),
          cylinder_radius = Float64(cylinder_radius), sphere_radius = Float64(sphere_radius),
          curvature)
    out = [TaLISRow[] for _ in vb]
    Threads.@threads for bi in eachindex(vb)
        _talis_beam!(out[bi], t, vb[bi], bls[bi], cx)
    end
    rows = reduce(vcat, out; init = TaLISRow[])
    append!(rows, _talis_stay_in(t, vb, surf, Float64(stay_in), Float64(hardware_tolerance)))
    notes = String[]
    skipped = [b.name for (b, bl) in zip(vb, bls) if bl === nothing]
    omega1 && !isempty(skipped) && !isempty(o1_idx) &&
        push!(notes, "1ω rules for beams without a NIF quad name or an `offsets` entry: " *
                     join(first(skipped, 12), ", ") * (length(skipped) > 12 ? ", …" : ""))
    omega1 || push!(notes, "1ω rules (omega1 = false)")
    any(_is_nif, vb) || push!(notes, "TAS stay-in zone (no NIF beams)")
    settings = (; tolerance = Float64(tolerance), pivot = tol.pivot, pivot_point = tol.point,
                pivot_text = _pivot_text(tol), transverse = Float64(R),
                pointing = Float64(pointing), target_offset = Float64(target_offset),
                hardware_tolerance = Float64(hardware_tolerance), focus = Float64(focus),
                ssd = Float64(ssd), capture_contour = Float64(capture_contour),
                capture_margin = Float64(capture_margin), max_incidence = Float64(max_incidence),
                clearance_contour = Float64(clearance_contour),
                clearance_margin = Float64(clearance_margin), omega1_margin = Float64(omega1_margin),
                omega1_edge = Float64(omega1_edge), omega1_size = W0, omega1_focus = f1,
                lens_distance = Float64(lens_distance), reflection_angle = Float64(reflection_angle),
                cylinder_radius = Float64(cylinder_radius), sphere_radius = Float64(sphere_radius),
                open_ports = length(ports), stay_in = Float64(stay_in),
                capture = Symbol[first(g) for g in cap_idx])
    return TaLISReport(rows, settings, notes)
end

# Per-beam checks (3ω capture, incidence and clearance; 1ω clearance, management and
# diagnostics). `bl` holds the beam's 1ω beamlets, or `nothing`.
function _talis_beam!(rows, t::Target, b::Beam, bl, cx)
    aim = Vec3(b.aim)
    tol = cx.tol
    k = direction(b)
    env1 = BeamEnvelope(b; peak_fraction = cx.clearance_contour, margin = cx.ssd)
    env10 = BeamEnvelope(b; peak_fraction = cx.capture_contour, margin = cx.ssd)
    shift = cx.R + cx.focus * env1.cone          # exact worst-case transverse loss
    land, s_stop = 0, Inf                        # index of the landing item in cap_idx
    for (j, g) in enumerate(cx.cap_idx)
        s = _landing(t, env1, last(g); miss = Inf)
        s < s_stop && ((land, s_stop) = (j, s))
    end
    landed = land == 0 ? Int[] : last(cx.cap_idx[land])

    # 3ω capture and incidence
    if !isempty(cx.cap_idx)
        if land == 0
            push!(rows, TaLISRow(b.name, :capture_3w, :none, :fail, -Inf, :mm, -Inf, -Inf,
                                 cx.capture_margin, aim,
                                 "beam axis misses every capture component; " * _TALIS_BLOWBY))
        else
            name = first(cx.cap_idx[land])
            worst, nom, wp = Inf, Inf, aim
            for x in cx.edges[land]
                d, _ = envelope_distance(env10, x)
                nom = min(nom, d)
                w = d - shift - tol(x, aim)
                w < worst && ((worst, wp) = (w, x))
            end
            m = worst - cx.capture_margin
            note = !isfinite(worst) ? "no free substrate edge" : m >= 0 ? "" : _TALIS_BLOWBY
            push!(rows, TaLISRow(b.name, :capture_3w, name, m >= 0 ? :pass : :fail, m, :mm, worst,
                                 nom, cx.capture_margin, wp, note))
            θw, θn, ip = _talis_incidence(t, cx.cap_geo[land], env10, k, aim, shift, tol)
            worst = θw + cx.tolerance
            m = cx.max_incidence - worst
            push!(rows, TaLISRow(b.name, :incidence_3w, name, m > 0 ? :pass : :fail, m, :deg, worst,
                                 θn, cx.max_incidence, ip, m > 0 ? "" : _TALIS_BLOWBY))
        end
    end

    # 3ω clearance
    for (name, cis) in cx.hw_idx
        any(in(landed), cis) && continue
        worst, nom, wp = Inf, Inf, aim
        for ci in cis, x in cx.geo[ci].verts
            d, s = envelope_distance(env1, x)
            s <= s_stop || continue
            nom = min(nom, d)
            w = d - shift - tol(x, aim)
            w < worst && ((worst, wp) = (w, x))
        end
        m = worst - cx.clearance_margin
        note = !isfinite(worst) ? "downstream of the landing point" : m > 0 ? "" : _TALIS_BLOWBY
        push!(rows, TaLISRow(b.name, :clearance_3w, name, m > 0 ? :pass : :fail, m,
                             :mm, worst, nom, cx.clearance_margin, wp, note))
    end

    bl === nothing && return rows
    # 1ω clearance and management
    for (name, cis) in cx.o1_idx
        worst, nom, wp = Inf, Inf, aim
        for ci in cis, x in cx.geo[ci].verts
            dn, dw = _omega1_distance(x, bl, cx)
            nom = min(nom, dn)
            w = dw - cx.R - tol(x, aim)
            w < worst && ((worst, wp) = (w, x))
        end
        m = worst - cx.omega1_margin
        if m > 0
            push!(rows, TaLISRow(b.name, :clearance_1w, name, :pass, m, :mm, worst, nom,
                                 cx.omega1_margin, wp, ""))
        else
            push!(rows, TaLISRow(b.name, :clearance_1w, name, :review, m, :mm, worst, nom,
                                 cx.omega1_margin, wp,
                                 "1ω footprint does not clear: unconverted-light management required (see management_1w)"))
            push!(rows, _talis_management(t, name, cis, b, bl, cx))
        end
    end
    # diagnostics to 1ω light
    for (name, cis) in cx.dg_idx
        worst, nom, req, wp = Inf, Inf, NaN, aim
        best = Inf
        for ci in cis, x in cx.geo[ci].verts
            dn, dw = _omega1_distance(x, bl, cx)
            w = dw - cx.R - tol(x, aim)
            r = 10 + 0.015 * norm(x)
            nom = min(nom, dn)
            w - r < best && ((best, worst, req, wp) = (w - r, w, r, x))
        end
        push!(rows, TaLISRow(b.name, :diagnostic_1w, name, best >= 0 ? :pass : :fail,
                             best, :mm, worst, nom, req, wp,
                             "required clearance 10 mm at TCC → 25 mm at 1 m (linear in distance from TCC)"))
    end
    return rows
end

# Largest incidence angle (degrees) on a landing component inside the worst-case and
# nominal footprints, over sub-triangles visible from upstream. The sub-triangle nearest
# the beam axis is always included.
function _talis_incidence(t::Target, g, env, k, aim, shift, tol)
    far = 4 * (t.scale + norm(aim)) + 1
    ε = 1e3 * _ray_eps(t)
    θw, θn, ip = -Inf, -Inf, aim
    bestd, θc, pc = Inf, NaN, aim
    for (leaf, n) in zip(g.leaves, g.normals)
        dn = minimum(envelope_distance(env, p)[1] for p in leaf)
        dw = minimum(envelope_distance(env, p)[1] - tol(p, aim) for p in leaf)
        (dw <= shift || dn < bestd) || continue
        m = leaf[4]
        visible(t, m - far * k, m - ε * k) || continue
        θ = acosd(clamp(abs(dot(n, k)), 0.0, 1.0))
        dn < bestd && ((bestd, θc, pc) = (dn, θ, m))
        dn <= 0 && (θn = max(θn, θ))
        dw <= shift && θ > θw && ((θw, ip) = (θ, m))
    end
    if isfinite(θc)
        θn = max(θn, θc)
        θc > θw && ((θw, ip) = (θc, pc))
    end
    return θw, θn, ip
end

# 1ω management of an uncleared item `name` made of the components `cis`: the worst of
# its components' rows (status first, then margin), under the item's name.
function _talis_management(t::Target, name::Symbol, cis, b::Beam, bl, cx)
    rank(r) = r.status === :fail ? 0 : r.status === :review ? 1 : 2
    best = nothing
    for ci in cis
        r = _talis_management(t, ci, b, bl, cx)
        (best === nothing || (rank(r), r.margin) < (rank(best), best.margin)) && (best = r)
    end
    r = best
    return r.component === name ? r :
           TaLISRow(r.beam, r.rule, name, r.status, r.margin, r.unit, r.worst, r.nominal, r.required,
                    r.point, r.note)
end

# 1ω management of an uncleared component: curvature, else specular reflection.
function _talis_management(t::Target, ci::Integer, b::Beam, bl, cx)
    c = t.components[ci]
    aim = Vec3(b.aim)
    kind, rc, half = _talis_curvature(t, ci, cx.curvature)
    lim = kind === :cylinder ? cx.cylinder_radius :
          (kind === :sphere || kind === :dimpled) ? cx.sphere_radius : NaN
    if kind !== :flat && kind !== :unknown && rc < lim
        what = kind === :dimpled ? "spherical dimple" : string(kind)
        note = @sprintf("%s radius of curvature %.3g mm < %.3g mm scatters 1ω", what, rc, lim)
        # NIF User Guide §5.6.2: dimple half-side / radius must not be below 0.25
        if kind === :dimpled && isfinite(half) && half / rc < 0.25
            return TaLISRow(b.name, :management_1w, c.name, :review, lim - rc, :mm, rc, rc, lim, aim,
                            note * @sprintf("; dimple half-side/radius %.3g < 0.25 (User Guide §5.6.2: periodicity too short)",
                                            half / rc))
        end
        return TaLISRow(b.name, :management_1w, c.name, :pass, lim - rc, :mm, rc, rc, lim, aim, note)
    end
    g = cx.geo[ci]
    αw, αn, wp = Inf, Inf, aim
    for (leaf, n) in zip(g.leaves, g.normals)
        hit, hitn = false, false
        for p in leaf
            dn, dw = _omega1_distance(p, bl, cx)
            dw - cx.R - cx.tol(p, aim) <= cx.omega1_margin && (hit = true)
            dn <= 0 && (hitn = true)
        end
        hit || continue
        for (kb, _, _) in bl
            r = kb - 2 * dot(kb, n) * n
            α = minimum(acosd(clamp(dot(r, p), -1.0, 1.0)) for p in cx.ports)
            α < αw && ((αw, wp) = (α, leaf[4]))
            hitn && (αn = min(αn, α))
        end
    end
    worst = αw - 2 * cx.tolerance
    m = worst - cx.reflection_angle
    desc = kind === :unknown ? "curvature unknown (give `curvature`)" :
           kind === :flat ? "flat" : @sprintf("%s radius %.3g mm ≥ %.3g mm", kind, rc, lim)
    status = m > 0 ? :pass : kind === :unknown ? :review : :fail
    note = m > 0 ? "$desc; specular 1ω reflection clears open quads" :
           "$desc; specular 1ω reflection within $(cx.reflection_angle)° of an open quad: TaLIS additional review, may not be allowed"
    return TaLISRow(b.name, :management_1w, c.name, status, m, :deg, worst, αn,
                    cx.reflection_angle, wp, note)
end

# (kind, radius, dimple half-side): kind :cylinder, :sphere, :dimpled, :flat or :unknown.
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
    ns = [t.tri_normal[k] for k in eachindex(t.tri_area) if t.tri_comp[k] == ci && t.tri_area[k] > 0]
    all(n -> abs(dot(n, ns[1])) > 1 - 1e-9, ns) && return (:flat, Inf, NaN)
    return (:unknown, NaN, NaN)
end

# Transverse signed distance (negative inside) from `x` to the union of a beam's 1ω
# footprints: at nominal focus, and the minimum over focus offsets 0, ±focus.
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

# Footprint square at axial position s = d·k from the 3ω focus: centre c0(1 + s/L) along
# ô (chief ray from the lens centre), half-side h0|1 − s/f1| (cone to the 1ω focus).
function _omega1_square(d, k, ô, c0, h0, f1, L)
    s = dot(d, k)
    c = c0 * (1 + s / L)
    h = h0 * abs(1 - s / f1)
    p = dot(d, ô) - c
    q = dot(d, cross(k, ô))
    ox, oy = abs(p) - h, abs(q) - h
    return (ox > 0 || oy > 0) ? hypot(max(ox, 0.0), max(oy, 0.0)) : max(ox, oy)
end

_is_nif(b::Beam) = b.quad in values(NIF_QUAD_NAMES)

function _nif_quad_port(quad)
    p = nif_port(quad)
    for (n, θ, φ) in NIF_PORTS
        n == p && return port_vector(θ, φ)
    end
end

# 1ω beamlets of a beam as (k, 3ω focus, footprint direction ô), or `nothing`.
function _omega1_beamlets(b::Beam, offsets, table)
    aim = Vec3(b.aim)
    mk(port, ô) = (-port, aim - b.defocus * port, ô)
    if offsets !== nothing && haskey(offsets, b.name)
        dir = offsets[b.name]
        k = direction(b)
        e1, e2 = port_basis(b.port)
        ô = dir isa Real ? cosd(dir) * e1 + sind(dir) * e2 :
            (v = Vec3(dir...); normalize(v - dot(v, k) * k))
        return [mk(b.port, ô)]
    end
    _is_nif(b) || return nothing
    pq = _nif_quad_port(b.quad)
    ports = acosd(clamp(dot(b.port, pq), -1.0, 1.0)) < 0.1 ?
            [port_vector(r.theta, r.phi) for r in table if "Q" * r.bundle * r.hemisphere == b.quad] :
            [b.port]
    return [mk(p, _own_side(p, pq)) for p in ports]
end

# φ̂ of the beamlet port, signed toward the side of the quad axis the port lies on.
function _own_side(p, pq)
    _, eφ = port_basis(p)
    return dot(p - pq, eφ) >= 0 ? eφ : -eφ
end

function _talis_open_ports(spec, beams)
    spec === :nif && return [port_vector(θ, φ) for (_, θ, φ) in NIF_PORTS]
    spec === :beams && return unique([_is_nif(b) ? _nif_quad_port(b.quad) : b.port for b in beams])
    spec isa Symbol && throw(ArgumentError("open_ports must be :nif, :beams or a list of ports"))
    return [x isa Beam ? x.port : normalize(Vec3(x...)) for x in spec]
end

# TAS stay-in zone (NIF beams): extent of each component above the lowest upper-beam aim
# point and below the highest lower-beam aim point.
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

# Vertices, sub-triangles (a, b, c, centroid) and their normals (exact on sphere,
# cylinder and cone surfaces) of a component.
function _talis_samples(t::Target, ci, level)
    c = t.components[ci]
    curved = c.surface isa Union{SphereSurface,CylinderSurface,ConeSurface,TorusSurface}
    verts, leaves, normals, buf = Vec3[], NTuple{4,Vec3}[], Vec3[], Vec3[]
    for k in eachindex(t.tri_area)
        (t.tri_comp[k] == ci && t.tri_area[k] > 0) || continue
        v0, e1, e2 = t.tri_v0[k], t.tri_e1[k], t.tri_e2[k]
        empty!(buf)
        _subtri_vertices!(buf, v0, v0 + e1, v0 + e2, level)
        append!(verts, buf)
        for j in 1:3:length(buf)
            a, bb, d = buf[j], buf[j+1], buf[j+2]
            m = (a + bb + d) / 3
            push!(leaves, (a, bb, d, m))
            push!(normals, curved ? last(surface_point(c.surface, m)) : t.tri_normal[k])
        end
    end
    return (; verts = unique(verts), leaves, normals)
end

# The samples of several components (the pieces of a part) as one.
_merge_samples(gs) = (; verts = unique(reduce(vcat, (g.verts for g in gs))),
                      leaves = reduce(vcat, (g.leaves for g in gs)),
                      normals = reduce(vcat, (g.normals for g in gs)))

# Points every ≤ `step` mm along the free edges (one opaque triangle) of the given
# components; vertices are merged by position across components.
function _free_edge_points(t::Target, comps, step)
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
    pts = Dict(ci => Vec3[] for ci in comps)
    for (n, ci, a, b) in values(seen)
        (n == 1 && haskey(pts, ci)) || continue
        m = max(1, ceil(Int, norm(b - a) / step))
        for i in 0:m
            push!(pts[ci], a + (i / m) * (b - a))
        end
    end
    return Dict(ci => unique(p) for (ci, p) in pts)
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

Rows of a TaLIS report whose status is in `status`; failures first, then by margin.
"""
violations(r::TaLISReport; status = (:fail, :review)) =
    sort(filter(x -> x.status in status, r.rows); by = x -> (_talis_rank(x), x.margin))

function Base.show(io::IO, ::MIME"text/plain", r::TaLISReport)
    s = r.settings
    nb = length(unique(x.beam for x in r.rows if x.rule !== :tas_stay_in))
    println(io, "TaLISReport: ", nb, " beams, ", length(r.rows), " checks against ", TALIS_SOURCE)
    @printf(io, "  rule tolerance %.3g° applied as a lateral allowance r·tan(%.3g°), %s\n",
            s.tolerance, s.tolerance, s.pivot_text)
    println(io, "    ASSUMED: the source does not say how the 0.5° is applied; confirm the interpretation with TaLIS")
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

Write every row of a TaLIS report as CSV (margins in the row's `unit`).
"""
function write_talis_csv(path::AbstractString, r::TaLISReport)
    open(path, "w") do io
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
