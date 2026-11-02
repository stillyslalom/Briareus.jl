# The hohlraum part and its pieces: the starter "Hohlraum with capsule" against the
# former hohlraum template (meshes and drive, names mapped), pieces and ends, the
# entrance-hole and window frames, beams aimed in piece frames, roles, unique names on
# insert.

using LinearAlgebra: norm

# The template's component name of a starter component (`hohlraum_wall` → `wall`).
_template_name(n::Symbol) = Symbol(replace(String(n), r"^hohlraum_" => ""))

# A coarse session of the starter "Hohlraum with capsule" alone (the template's design).
function _starter_session(; windows = WindowSpec[], closure = false)
    sess = DesignSession(resolution = :coarse)
    ps = starter_parts(:hohlraum)
    v = copy(ps[1].values)
    v[:leh_closure] = closure ? 1.0 : 0.0
    ps[1] = Controllers._with(ps[1]; values = v, windows)
    @test Controllers._try_build!(sess.target; kind = :custom, parts = ps)
    return sess
end

# Whether two components are the same up to their names and parts.
_same_component(a, b) =
    a.mesh.vertices == b.mesh.vertices && a.mesh.triangles == b.mesh.triangles &&
    a.mesh.patches == b.mesh.patches && a.kind == b.kind && a.two_sided == b.two_sided &&
    a.surface == b.surface && albedo(a.material) == albedo(b.material) &&
    a.material.conversion == b.material.conversion

@testset "hohlraum part: the starter equals the template (no GL)" begin
    tmpl = DesignSession(resolution = :coarse)
    st = _starter_session()
    ta, tb = tmpl.target.target[], st.target.target[]
    # element for element, in the same order, names mapped
    @test _template_name.(c.name for c in tb.components) == [c.name for c in ta.components]
    @test all(_same_component(a, b) for (a, b) in zip(ta.components, tb.components))
    @test all(c -> c.part === :hohlraum, tb.components[1:end-1]) && tb.components[end].part === :capsule
    @test ta.elem_area == tb.elem_area && ta.elem_centroid == tb.elem_centroid && ta.elem_normal == tb.elem_normal
    @test Controllers._closed_enclosure(st.target) && Controllers._closed_enclosure(tmpl.target)
    @test observed_name(st.target) === :capsule
    # the drive on Float64: bit for bit
    build_scenario!(tmpl)
    build_scenario!(st)
    ra, rb = tmpl.drive.result[], st.drive.result[]
    @test [(m.l, m.legendre, m.epsilon) for m in ra.modes] == [(m.l, m.legendre, m.epsilon) for m in rb.modes]
    @test ra.mean_flux == rb.mean_flux
    @test all(element_values(ra, f) == element_values(rb, f) for f in (:incident, :radiosity))
    @test isequal(element_values(ra, :laser), element_values(rb, :laser))      # (NaN on apertures)
    @test flux_map(tmpl.drive).values == flux_map(st.drive).values
    # the power accounting: a row per piece, named as the template's components
    @test [(p.name, p.power) for p in ra.power] == [(p.name, p.power) for p in rb.power]
    # the checks' default roles: the wall lands the beams, the caps stay hardware (LEH edges)
    @test sort(_template_name.(st.review.targets[])) == sort(tmpl.review.targets[])
    @test sort(_template_name.(st.talis.capture[])) == sort(tmpl.talis.capture[])
    @test sort(_template_name.(st.review.hardware[])) == sort(tmpl.review.hardware[])

    # a window and closing entrance holes: the same meshes; the capsule comes after the
    # window now (a part of its own), so the drive agrees to round-off
    w = WindowSpec(; phi = 90.0, albedo = 0.85, conversion = 0.8)
    add_window!(tmpl.target)
    set_leh_closure!(tmpl.target, true)
    stw = _starter_session(; windows = [w], closure = true)
    ta, tb = tmpl.target.target[], stw.target.target[]
    @test tmpl.target.windows[] == [w]
    byname = Dict(_template_name(c.name) => c for c in tb.components)
    @test Set(keys(byname)) == Set(c.name for c in ta.components)
    @test all(_same_component(c, byname[c.name]) for c in ta.components)
    build_scenario!(tmpl)
    build_scenario!(stw)
    ra, rb = tmpl.drive.result[], stw.drive.result[]
    @test isapprox(ra.mean_flux, rb.mean_flux; rtol = 1e-10)
    @test all(isapprox(a.legendre, b.legendre; rtol = 1e-8, atol = 1e-12) for (a, b) in zip(ra.modes, rb.modes))
    inc(t, r, n) = element_values(r, :incident)[elements(t, n)]
    @test all(isapprox(inc(ta, ra, c.name), inc(tb, rb, byname[c.name].name); rtol = 1e-9) for c in ta.components)
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
    @test_throws ArgumentError part_components(PartSpec(:hohlraum; name = :hr, values = Dict(:top_end => 2.0)), 0.6)

    # the frames: the ends' (z out of the hohlraum at either end) and the window's, which
    # follow the part when it moves and turns
    sess = DesignSession(resolution = :coarse)
    ts = sess.target
    @test insert_starter!(ts, :hohlraum)
    i = part_index(ts, :hohlraum)
    @test edit_part!(ts, i, :z, "0.4") && edit_part!(ts, i, :theta, "90") && edit_part!(ts, i, :length, "8")
    fr = part_frames(ts)
    F = fr[:hohlraum]
    top, bot = fr[Symbol("hohlraum.leh_top")], fr[Symbol("hohlraum.leh_bottom")]
    @test top.origin ≈ to_chamber(F, (0.0, 0.0, 4.0)) atol = 1e-12
    @test bot.origin ≈ to_chamber(F, (0.0, 0.0, -4.0)) atol = 1e-12
    @test isapprox(top.R[:, 3], F.R[:, 3]; atol = 1e-12) && isapprox(bot.R[:, 3], -F.R[:, 3]; atol = 1e-12)
    @test bot.R[:, 1] ≈ F.R[:, 1] atol = 1e-12
    # the capsule placed in the hohlraum's frame moved with it (`capsule_2`: the
    # template's capsule has the name)
    @test ts.parts[][2].name === :capsule_2 && isapprox(fr[:capsule_2].origin, F.origin; atol = 1e-12)
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
    @test Controllers._try_build!(ts; kind = :custom, parts = starter_parts(:hohlraum))
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
    @test isempty(wall_names(ts)) && !Controllers._closed_enclosure(ts)
    set_excluded!(ts, :hohlraum, false)

    # a second insert: names made unique; the new capsule takes the observed role
    @test insert_starter!(ts, :hohlraum; frame = Symbol("hohlraum.leh_top"))
    @test [p.name for p in ts.parts[]] == [:hohlraum, :capsule, :hohlraum_2, :capsule_2]
    @test ts.parts[][3].frame === Symbol("hohlraum.leh_top") && ts.parts[][4].frame === :hohlraum_2
    @test observed_name(ts) === :capsule_2 && !ts.parts[][2].observed
    @test !Controllers._closed_enclosure(ts)                  # two cavities: not one enclosure
    @test :hohlraum_2_wall in wall_names(ts)
    @test !insert_starter!(ts, :hohlraum; frame = :nowhere) && occursin("no frame", ts.error[])
end
