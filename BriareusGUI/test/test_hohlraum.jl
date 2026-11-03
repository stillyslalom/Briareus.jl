# The hohlraum part and its pieces: the starter "Hohlraum with capsule" against core
# `hohlraum` (the former template's geometry, names mapped), sessions with a template
# converted to parts against the template-era GUI's own drives, pieces and ends, the
# entrance-hole and window frames, beams aimed in piece frames, roles, unique names on
# insert.

using LinearAlgebra: norm

# The former template's component name of a starter component (`hohlraum_wall` → `wall`).
_template_name(n::Symbol) = Symbol(replace(String(n), r"^hohlraum_" => ""))

# Whether every beam aims where its frame and local aim put it.
_in_frames(bt) = all(i -> isapprox(bt.beams[][i].aim, to_chamber(bt.frames[][bt.aim_frame[][i]], bt.aim_local[][i]); atol = 1e-12),
                     eachindex(bt.beams[]))

# Whether two components are the same up to their names and parts.
_same_component(a, b) =
    a.mesh.vertices == b.mesh.vertices && a.mesh.triangles == b.mesh.triangles &&
    a.mesh.patches == b.mesh.patches && a.kind == b.kind && a.two_sided == b.two_sided &&
    a.surface == b.surface && albedo(a.material) == albedo(b.material) &&
    a.material.conversion == b.material.conversion

@testset "hohlraum part: the starter meshes as the former template (no GL)" begin
    ts = TargetSetup(resolution = :coarse)
    m = ts.mesh[]
    # the former template's call: core `hohlraum` at TCC, its capsule the geodesic mesh
    ref = hohlraum(; length = 10.1, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9,
                   wall = gold(albedo = 0.85, conversion = 0.8), h = m.h, capsule_h = m.capsule_h,
                   leh_open_radius = 0.84 * 1.55)
    t = ts.target[]
    @test _template_name.(c.name for c in t.components) == [c.name for c in ref.components]
    @test all(_same_component(a, b) for (a, b) in zip(t.components, ref.components))
    @test all(c -> c.part === :hohlraum, t.components[1:end-1]) && t.components[end].part === :capsule
    tr = Target(ref.components)
    @test t.elem_area == tr.elem_area && t.elem_centroid == tr.elem_centroid && t.elem_normal == tr.elem_normal
    @test closed_design(ts) && observed_name(ts) === :capsule
    # with a window and closing entrance holes: the same meshes, the capsule after the window
    h = ts.parts[][1]
    v = copy(h.values)
    v[:leh_closure] = 1.0
    w = WindowSpec(; phi = 90.0, albedo = 0.85, conversion = 0.8)
    tw = TargetSetup(resolution = :coarse, parts = [Controllers._with(h; values = v, windows = [w]), ts.parts[][2]])
    refw = hohlraum(; length = 10.1, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9,
                    wall = gold(albedo = 0.85, conversion = 0.8), h = m.h, capsule_h = m.capsule_h,
                    windows = [Controllers._wall_window(w, 1)],
                    leh_closure = GapClosure(closure = 5.0, width = 0.5, after = 0.85), leh_open_radius = 0.84 * 1.55)
    byname = Dict(_template_name(c.name) => c for c in tw.target[].components)
    @test Set(keys(byname)) == Set(c.name for c in refw.components)
    @test all(_same_component(c, byname[c.name]) for c in refw.components)
end

@testset "hohlraum part: sessions with a template convert to parts (no GL)" begin
    data = joinpath(@__DIR__, "data")
    old = TOML.parsefile(joinpath(data, "session_v7_expected.toml"))
    # a hohlraum template with a window, closing holes, a shield at the top LEH, beams
    # aimed in the bottom LEH's frame and display modes
    s = DesignSession(resolution = :coarse)
    load_session!(s, joinpath(data, "session_v7_hohlraum.briar.toml"))
    ts, bt = s.target, s.beams
    e = old["hohlraum"]
    @test [p.name for p in ts.parts[]] == [:hohlraum, :capsule, :shield]
    @test ts.parts[][1].kind === :hohlraum && ts.parts[][1].values[:leh_closure] == 1 &&
          only(ts.parts[][1].windows).phi == 60 && ts.parts[][2].observed && ts.parts[][2].frame === :hohlraum
    @test ts.parts[][3].frame === Symbol("hohlraum.leh_top")
    @test bt.aim_frame[] == [f == "tcc" ? :tcc : Symbol("hohlraum.", f) for f in e["aim_frames"]]
    @test [collect(b.aim) for b in bt.beams[]] == e["aims"] && _in_frames(bt)
    @test display_mode(s.view, :hohlraum_wall) === :wireframe && display_mode(s.view, :hohlraum_leh_top) === :transparent
    t = ts.target[]
    @test sort(_template_name.(c.name for c in t.components)) == sort(Symbol.(e["components"]))
    @test sort(_template_name.(s.review.targets[])) == sort(Symbol.(e["review_targets"]))
    @test sort(_template_name.(s.review.hardware[])) == sort(Symbol.(e["review_hardware"]))
    @test _template_name.(s.talis.capture[]) == Symbol.(e["talis_capture"])
    # the drive of the template-era GUI on the same session (its capsule now comes after the
    # window: to round-off), which left this design's view factors unbalanced (the shield
    # outside the hohlraum): the scenario without its cavities
    build_scenario!(s)
    rb = s.drive.result[]
    out = Controllers._build_scenario(t, bt.beams[]; observed = observed_name(ts),
                                      profile = Controllers._deposition_profile(s.drive))
    r = Controllers._evaluate_drive(out.s, display_beams(bt), out.lookup, s.drive.spot_radius[];
                                    ratio = unconverted_ratios(s), foil = active_foil(ts))
    # (the drive balances the hohlraum's cavity, the shield in none: +0.33 % mean flux)
    c = only(Controllers._closed_cavities(ts))
    @test c.balanced && isempty(intersect(c.elements, elements(t, :shield)))
    @test isapprox(rb.mean_flux, e["mean_flux"]; rtol = 5e-3) && !isapprox(rb.mean_flux, e["mean_flux"]; rtol = 1e-3)
    @test isapprox(r.mean_flux, e["mean_flux"]; rtol = 1e-10)
    @test all(isapprox(m.legendre, l; rtol = 1e-8, atol = 1e-12) for (m, l) in zip(r.modes, e["legendre"]))
    inc = element_values(r, :incident)
    sums = Dict(_template_name(c.name) => sum(inc[t.comp_elems[ci]]) for (ci, c) in enumerate(t.components))
    # the holes and lips have twice their template-era elements (the caps' azimuthal count
    # since 2026-11-03): their sums of element fluxes double, to 1.5e-3 of the element means
    hole(k) = startswith(k, "leh")
    @test all(isapprox(sums[Symbol(k)], x; rtol = 1e-10) for (k, x) in e["incident_sum"] if !hole(k))
    @test all(isapprox(sums[Symbol(k)] / 2, x; rtol = 3e-3) for (k, x) in e["incident_sum"] if hole(k))
    # saved again: the current version, read back the same
    f = joinpath(mktempdir(), "v10.briar.toml")
    save_session(f, s)
    @test TOML.parsefile(f)["version"] == Controllers.SESSION_VERSION
    s8 = DesignSession(resolution = :coarse)
    load_session!(s8, f)
    @test session_state(s8) == session_state(s)

    # a foil template with a substrate in the foil's frame and the OMEGA backlighter beams
    b = DesignSession(starter = :foil, resolution = :coarse)
    load_session!(b, joinpath(data, "session_v7_foil.briar.toml"))
    e = old["foil"]
    @test [p.name for p in b.target.parts[]] == [:foil, :main, :substrate]
    @test b.target.parts[][1].foil && active_foil(b.target) === :foil && b.target.parts[][3].frame === :foil
    @test b.beams.aim_frame[] == [f == "backlighter" ? :foil : Symbol(f) for f in e["aim_frames"]]
    @test [collect(x.aim) for x in b.beams.beams[]] == e["aims"]
    t = b.target.target[]
    @test nelements(t) == e["elements"]
    build_scenario!(b)
    r = b.drive.result[]
    las = element_values(r, :laser)
    names(n) = n === :foil ? :backlighter : n
    sums = Dict(names(c.name) => sum(filter(isfinite, las[t.comp_elems[ci]]); init = 0.0) for (ci, c) in enumerate(t.components))
    @test all(isapprox(sums[Symbol(k)], x; rtol = 1e-12) for (k, x) in e["laser_sum"])
    @test isapprox(r.spot.mean, e["spot_mean"]; rtol = 1e-12) && isapprox(r.spot.peak, e["spot_peak"]; rtol = 1e-12)
end

@testset "hohlraum part: pieces, ends, frames, aims (no GL)" begin
    p = PartSpec(:hohlraum; name = :hr, windows = [WindowSpec(; phi = 30.0, z = 1.0)],
                 values = Dict(:leh_closure => 1.0))
    @test hohlraum_pieces(p) == [:wall, :cap_top, :cap_bottom, :leh_top, :leh_top_lip, :leh_bottom,
                                 :leh_bottom_lip, :window_1]
    cs = part_components(p, 0.6)
    @test [c.name for c in cs] == [piece_name(:hr, s) for s in hohlraum_pieces(p)]
    @test all(c -> c.part === :hr, cs)
    @test_throws ArgumentError PartSpec(:sphere; windows = [WindowSpec()])
    @test_throws ArgumentError PartSpec(:hohlraum; name = :hr, frame = Symbol("hr.leh_top"))

    # a closed end: its hole merged into its cap, of the wall's material
    q = PartSpec(:hohlraum; name = :hr, values = Dict(:top_end => 0.0))
    @test hohlraum_pieces(q) == [:wall, :cap_top, :cap_bottom, :leh_bottom]
    open, closed = part_components(PartSpec(:hohlraum; name = :hr), 0.6), part_components(q, 0.6)
    area(c) = sum(norm(Briareus.cross(c.mesh.vertices[b] - c.mesh.vertices[a], c.mesh.vertices[d] - c.mesh.vertices[a])) / 2
                  for (a, b, d) in c.mesh.triangles)
    cap(cs) = cs[findfirst(c -> c.name === :hr_cap_top, cs)]
    leh = open[findfirst(c -> c.name === :hr_leh_top, open)]
    @test area(cap(closed)) ≈ area(cap(open)) + area(leh) rtol = 1e-12
    @test cap(closed).kind === :surface && albedo(cap(closed).material) == 0.85
    @test_throws ArgumentError part_components(PartSpec(:hohlraum; name = :hr, values = Dict(:top_end => 3.0)), 0.6)

    # the frames: the ends' (z out of the hohlraum at either end) and the window's, which
    # follow the part when it moves and turns
    sess = DesignSession(starter = :none, resolution = :coarse)
    ts = sess.target
    @test insert_starter!(ts, :hohlraum) && [p.name for p in ts.parts[]] == [:sphere_1, :hohlraum, :capsule]
    i = part_index(ts, :hohlraum)
    @test edit_part!(ts, i, :z, "0.4") && edit_part!(ts, i, :theta, "90") && edit_part!(ts, i, :length, "8")
    fr = part_frames(ts)
    F = fr[:hohlraum]
    top, bot = fr[Symbol("hohlraum.leh_top")], fr[Symbol("hohlraum.leh_bottom")]
    @test top.origin ≈ to_chamber(F, (0.0, 0.0, 4.0)) atol = 1e-12
    @test bot.origin ≈ to_chamber(F, (0.0, 0.0, -4.0)) atol = 1e-12
    @test isapprox(top.R[:, 3], F.R[:, 3]; atol = 1e-12) && isapprox(bot.R[:, 3], -F.R[:, 3]; atol = 1e-12)
    @test bot.R[:, 1] ≈ F.R[:, 1] atol = 1e-12
    # the capsule placed in the hohlraum's frame moved with it
    @test ts.parts[][3].frame === :hohlraum && isapprox(fr[:capsule].origin, F.origin; atol = 1e-12)
    @test piece_frame(:hohlraum_cap_bottom, :hohlraum) === Symbol("hohlraum.leh_bottom")
    @test piece_frame(:hohlraum_leh_top_lip, :hohlraum) === Symbol("hohlraum.leh_top")
    @test piece_frame(:hohlraum_wall, :hohlraum) === :hohlraum

    # beams aimed at a piece frame's center follow the hohlraum
    set_aim_frame!(sess, 1:2, Symbol("hohlraum.leh_top"); keep = false)
    @test all(b -> isapprox(b.aim, top.origin; atol = 1e-12), sess.beams.beams[][1:2])
    @test edit_part!(ts, i, :x, "0.3")
    top2 = part_frames(ts)[Symbol("hohlraum.leh_top")]
    @test all(b -> isapprox(b.aim, top2.origin; atol = 1e-12), sess.beams.beams[][1:2])
    @test frame_point(sess, Symbol("hohlraum.leh_bottom")) == part_frames(ts)[Symbol("hohlraum.leh_bottom")].origin
    # a part placed at an entrance hole of the hohlraum part
    @test add_part!(ts, :disk; name = :shield, frame = Symbol("hohlraum.leh_top"), z = 0.5, r1 = 2.0)
    sh = part_frames(ts)[:shield]
    @test sh.origin ≈ to_chamber(top2, (0.0, 0.0, 0.5)) atol = 1e-12
    @test part_placement(ts, part_index(ts, :shield)) isa Briareus.PlacementReport
end

@testset "hohlraum part: roles, exclusion, unique names on insert (no GL)" begin
    sess = DesignSession(resolution = :coarse)
    ts = sess.target
    @test wall_names(ts) == [:hohlraum_wall, :hohlraum_cap_top, :hohlraum_cap_bottom]
    @test Set(wall_components(ts)) == Set(wall_names(ts))
    # the checks land the beams on the wall and the capsule, the caps are LEH edges
    @test :hohlraum_wall in sess.review.targets[] && :capsule in sess.review.targets[]
    @test :hohlraum_cap_top in sess.review.hardware[] && :hohlraum_cap_top in sess.talis.hardware[]
    @test !set_part_role!(ts, part_index(ts, :hohlraum), :wall) && occursin("wall role", ts.error[])
    # pieces are named on their own, the part stands for all of them
    names = component_names(ts.target[])
    @test :hohlraum in names && :hohlraum_wall in names && :hohlraum_wall in surface_names(ts.target[])
    @test !(:hohlraum in surface_names(ts.target[]))
    # the tree: the part's row, its pieces under it
    sess.view.expanded[] = union(sess.view.expanded[], ("target",))
    rows = model_tree(sess)
    k = findfirst(r -> r.id == "part:hohlraum", rows)
    @test k !== nothing && rows[k+1].label == "Wall" && rows[k+1].depth == 2
    # leaving the part out leaves out every piece (and the enclosure is open)
    set_excluded!(ts, :hohlraum, true)
    @test all(c -> c.part !== :hohlraum, ts.target[].components)
    @test length(ts.ghost[].components) == length(hohlraum_pieces(ts.parts[][1]))
    @test isempty(wall_names(ts)) && !closed_design(ts)
    set_excluded!(ts, :hohlraum, false)

    # a second insert: names made unique; the new capsule takes the observed role
    @test insert_starter!(ts, :hohlraum; frame = Symbol("hohlraum.leh_top"))
    @test [p.name for p in ts.parts[]] == [:hohlraum, :capsule, :hohlraum_2, :capsule_2]
    @test ts.parts[][3].frame === Symbol("hohlraum.leh_top") && ts.parts[][4].frame === :hohlraum_2
    @test observed_name(ts) === :capsule_2 && !ts.parts[][2].observed
    @test !closed_design(ts)                  # two cavities: not one enclosure
    @test :hohlraum_2_wall in wall_names(ts)
    @test !insert_starter!(ts, :hohlraum; frame = :nowhere) && occursin("no frame", ts.error[])
end
