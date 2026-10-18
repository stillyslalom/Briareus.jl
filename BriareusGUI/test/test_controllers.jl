# The target, beam table and view controllers: no scenario, no GL.

@testset "TargetSetup (no GL)" begin
    ts = TargetSetup(resolution = :coarse)
    t0 = ts.target[]
    @test component_names(t0) == [:wall, :cap_top, :cap_bottom, :leh_top, :leh_bottom, :capsule]
    @test surface_names(t0) == [:wall, :cap_top, :cap_bottom, :capsule]
    @test occursin("hohlraum", target_summary(ts)) && occursin("$(nelements(t0)) elements", target_summary(ts))

    # an accepted edit rebuilds; a rejected one keeps the target and reports why
    @test edit_target!(ts, :length, "9.4")
    @test ts.target[] !== t0 && ts.params[][:length] == 9.4 && ts.error[] == ""
    t1 = ts.target[]
    @test !edit_target!(ts, :leh_radius, "3.5")               # wider than the hohlraum
    @test ts.target[] === t1 && ts.params[][:leh_radius] == 1.55 && occursin("leh_radius", ts.error[])
    @test !edit_target!(ts, :wall_albedo, "1.2")
    @test occursin("between 0 and 1", ts.error[])
    @test !edit_target!(ts, :radius, "abc")
    @test occursin("must be a number", ts.error[])
    @test edit_target!(ts, :radius, "2.875") && ts.error[] == ""      # unchanged value
    @test_throws ArgumentError edit_target!(ts, :height, "1")

    # no capsule
    @test edit_target!(ts, :capsule_radius, "0")
    @test !(:capsule in component_names(ts.target[]))
    edit_target!(ts, :capsule_radius, "0.9")

    # resolution and kind
    n = nelements(ts.target[])
    set_resolution!(ts, :standard)
    @test nelements(ts.target[]) > 4n
    set_resolution!(ts, :coarse)
    @test_throws ArgumentError set_resolution!(ts, :ultra)
    set_target_kind!(ts, :foil)
    @test component_names(ts.target[]) == [FOIL_NAME, :main]
    @test foil_center(ts) == Briareus.Vec3(6.0, 0, 0)
    @test !edit_target!(ts, :distance, "1.2")                 # inside the sphere
    @test edit_target!(ts, :main_radius, "0")
    @test component_names(ts.target[]) == [FOIL_NAME]
    @test_throws ArgumentError set_target_kind!(ts, :sphere)
    set_target_kind!(ts, :hohlraum)
    @test :capsule in component_names(ts.target[])            # the hohlraum edits were kept
    @test ts.params[][:length] == 9.4
end

@testset "BeamTable (no GL)" begin
    bt = BeamTable(preset_beams(:nif_quads); preset = :nif_quads)
    @test length(bt.beams[]) == 48 && all(bt.enabled[])
    rows = beam_rows(bt)
    @test [r.key for r in rows] == ["inner23", "inner30", "outer44", "outer50"]
    @test all(r -> r.count == 8 || r.count == 16, rows) && sum(r.count for r in rows) == 48
    r23 = rows[1]
    @test r23.z ≈ 3.09 && r23.r ≈ -1.09 && r23.power ≈ 8.0 && r23.enabled === true
    # published pointing, mirrored for the lower hemisphere
    @test all(b -> abs(b.aim[3]) ≈ 3.09, bt.beams[][r23.members])
    @test any(b -> b.aim[3] < 0, bt.beams[][r23.members])
    @test occursin("48 of 48 beams on", beams_summary(bt)) && occursin("384 TW", beams_summary(bt))

    # a group edit sets one field on every member and nothing else
    before = copy(bt.beams[])
    @test edit_beams!(bt, 1, :z, "3.3")
    after = bt.beams[]
    m = r23.members
    @test all(b -> abs(b.aim[3]) ≈ 3.3, after[m])
    @test [b.power for b in after] == [b.power for b in before]
    others = setdiff(1:48, m)
    @test after[others] == before[others]
    @test beam_rows(bt)[1].r ≈ -1.09 atol = 1e-12                   # r kept
    @test !edit_beams!(bt, 1, :power, "-2")
    @test occursin("must not be negative", bt.error[]) && bt.beams[] === after
    @test !edit_beams!(bt, 1, :z, "far")
    @test edit_beams!(bt, 1, :z, "3.3") && bt.error[] == ""
    @test_throws ArgumentError edit_beams!(bt, 1, :spot, "1")

    # per quad and per beam
    set_grouping!(bt, :quad)
    @test length(beam_rows(bt)) == 48
    q = beam_rows(bt)[5]
    @test edit_beams!(bt, 5, :defocus, "1.5")
    @test bt.beams[][only(q.members)].defocus == 1.5
    set_grouping!(bt, :beam)
    @test beam_rows(bt)[5].z ≈ bt.beams[][5].aim[3]                  # not mirrored per beam
    @test_throws ArgumentError set_grouping!(bt, :bundle)
    set_grouping!(bt, :group)

    # on/off: disabled beams keep their place and deliver no power
    set_enabled!(bt, 2, false)
    @test beam_rows(bt)[2].enabled === false
    @test length(active_beams(bt)) == 48 - beam_rows(bt)[2].count
    @test length(display_beams(bt)) == 48
    @test all(b -> b.power == 0, display_beams(bt)[beam_rows(bt)[2].members])
    m2 = beam_rows(bt)[2].members
    set_grouping!(bt, :quad)                                          # one quad back on
    set_enabled!(bt, findfirst(r -> r.members == [m2[1]], beam_rows(bt)), true)
    set_grouping!(bt, :group)
    @test beam_rows(bt)[2].enabled === missing
    set_enabled!(bt, 2, false)

    # CSV round trip
    path = joinpath(mktempdir(), "beams.csv")
    export_beams(bt, path)
    @test bt.path[] == path
    bt2 = BeamTable()
    @test beams_summary(bt2) == "no beams"
    import_beams!(bt2, path)
    @test bt2.preset[] === :custom && length(bt2.beams[]) == 48
    @test all(((a, b),) -> a.port ≈ b.port && a.aim ≈ b.aim && a.power ≈ b.power &&
                           a.defocus ≈ b.defocus && a.name == b.name,
              zip(bt.beams[], bt2.beams[]))
    export_beams(bt, path; enabled_only = true)
    @test length(read_beam_table(path)) == length(active_beams(bt))

    # presets
    for p in BEAM_PRESETS
        @test !isempty(preset_beams(p))
    end
    @test length(preset_beams(:nif_beams)) == 192
    om = preset_beams(:omega_hohlraum; params = TargetSetup(resolution = :coarse).params[])
    @test length(om) == 40 && all(b -> abs(b.aim[3]) ≈ 10.1 / 2, om)
    foilp = TargetSetup(kind = :foil, resolution = :coarse).params[]
    bl = preset_beams(:omega_backlighter; target_kind = :foil, params = foilp)
    @test all(b -> b.aim ≈ foil_center(foilp) && b.port[1] < -cosd(70), bl)
    @test_throws ArgumentError preset_beams(:lmj)
    set_grouping!(bt, :quad)
    load_preset!(bt, :omega_hohlraum)
    @test bt.preset[] === :omega_hohlraum && length(bt.enabled[]) == 40 && all(bt.enabled[])
    @test bt.by[] === :group                         # OMEGA beams have no quad names
    @test_throws ArgumentError set_grouping!(bt, :quad)
    @test_throws ArgumentError Controllers._check_grouping(bt.beams[], :quad)   # the optimizer's check
end

@testset "diagnostic windows (no GL)" begin
    ts = TargetSetup(resolution = :coarse)
    t0 = ts.target[]
    n0 = nelements(t0)
    @test add_window!(ts)                                   # closed, 1 × 1 mm at φ = 90°
    w1 = only(ts.windows[])
    @test w1.phi == 90 && w1.z == 0 && !w1.open && w1.albedo == ts.params[][:wall_albedo]
    @test :window_1 in surface_names(ts.target[]) && window_name(1) === :window_1
    @test occursin("with 1 window", target_summary(ts))
    @test add_window!(ts) && ts.windows[][2].phi == 270     # the next free default place
    t = ts.target[]

    # the geometry's rejections surface as messages and keep the target
    @test !edit_window!(ts, 2, :phi, "90.2")
    @test occursin("overlap", ts.error[]) && ts.target[] === t && ts.windows[][2].phi == 270
    @test !edit_window!(ts, 1, :z, "5.0")                   # past the end of the 10.1 mm wall
    @test occursin("past the end", ts.error[]) && ts.windows[][1].z == 0
    @test startswith(ts.error[], "window 1")                  # which window, for the user
    @test !edit_window!(ts, 1, :width, "0") && occursin("positive", ts.error[])
    @test !edit_window!(ts, 1, :albedo, "1.5") && occursin("between 0 and 1", ts.error[])
    @test !edit_window!(ts, 1, :height, "tall") && occursin("must be a number", ts.error[])
    @test !add_window!(ts; phi = 92.0) && occursin("overlap", ts.error[])
    @test length(ts.windows[]) == 2 && ts.target[] === t
    @test_throws ArgumentError edit_window!(ts, 1, :depth, "1")
    @test_throws ArgumentError edit_window!(ts, 3, :z, "1")
    @test_throws ArgumentError remove_window!(ts, 3)

    # accepted edits rebuild; a closed window re-emits with its own material
    @test edit_window!(ts, 1, :z, "1.5") && ts.error[] == "" && ts.target[] !== t
    @test edit_window!(ts, 1, :albedo, "0.3") && edit_window!(ts, 1, :conversion, "0.1")
    tw = ts.target[]
    c1 = tw.components[component_index(tw, :window_1)]
    @test c1.kind === :surface && albedo(c1.material) ≈ 0.3 && c1.material.conversion ≈ 0.1
    zs = [tw.elem_centroid[e][3] for e in elements(tw, :window_1)]
    @test all(z -> 1.0 < z < 2.0, zs)                       # the mesh conforms to the window
    # an open window is a hole (an aperture)
    @test set_window_open!(ts, 2, true)
    @test ts.target[].components[component_index(ts.target[], :window_2)].kind === :aperture
    @test !(:window_2 in surface_names(ts.target[]))
    # a dimension edit that leaves a window off the wall is rejected
    @test !edit_target!(ts, :length, "3.0") && occursin("past the end", ts.error[])
    @test ts.params[][:length] == 10.1
    # resolution and kind keep the windows; a foil has none
    set_resolution!(ts, :standard)
    @test :window_1 in component_names(ts.target[]) && nelements(ts.target[]) > 4n0
    set_resolution!(ts, :coarse)
    set_target_kind!(ts, :foil)
    @test !any(n -> startswith(String(n), "window"), component_names(ts.target[]))
    @test !occursin("window", target_summary(ts))
    @test !add_window!(ts) && occursin("hohlraum", ts.error[])
    set_target_kind!(ts, :hohlraum)
    @test component_names(ts.target[])[end-1:end] == [:window_1, :window_2]
    # removing renumbers: the open window becomes window_1
    remove_window!(ts, 1)
    @test length(ts.windows[]) == 1 && ts.windows[][1].open
    @test ts.target[].components[component_index(ts.target[], :window_1)].kind === :aperture
    remove_window!(ts, 1)
    @test isempty(ts.windows[]) && nelements(ts.target[]) == n0
    # windows given up front
    ts2 = TargetSetup(resolution = :coarse, windows = [WindowSpec(phi = 0.0, height = 2.0, open = true)])
    @test :window_1 in component_names(ts2.target[])
end

@testset "VisRad interchange through the beam table (no GL)" begin
    bt = BeamTable(preset_beams(:nif_quads); preset = :nif_quads)
    edit_beams!(bt, 1, :z, "3.3")
    edit_beams!(bt, 2, :defocus, "0.5")
    set_enabled!(bt, 4, false)                                 # the 50° cone off
    ref, en = copy(bt.beams[]), copy(bt.enabled[])
    dir = mktempdir()

    # SST pointing: one row per NIF beam of the enabled quads
    sst = export_visrad_sst(bt, joinpath(dir, "sst.csv"))
    @test bt.path[] == sst
    @test countlines(sst) == 1 + 4 * count(en)
    fresh() = BeamTable(preset_beams(:nif_quads); preset = :nif_quads)
    b2 = fresh()
    import_visrad_sst!(b2, sst)
    @test all(((a, b),) -> a.aim ≈ b.aim && a.defocus ≈ b.defocus,
              zip(bt.beams[][en], b2.beams[][en]))
    @test beam_rows(b2)[1].z ≈ 3.3 && beam_rows(b2)[2].defocus ≈ 0.5
    # the quads absent from the file come back disabled, at their power
    @test b2.enabled[] == en && all(b -> b.power == 8.0, b2.beams[])
    @test b2.preset[] === :nif_quads && b2.path[] == sst
    @test occursin("disabled", b2.warnings[][1])
    # the reader's warnings are collected: a beam the table does not have
    open(io -> println(io, "B999,1.0,2.0,3.0,0.0"), sst, "a")
    b2b = fresh()
    import_visrad_sst!(b2b, sst)
    @test any(w -> occursin("ignored", w) && occursin("B999", w), b2b.warnings[])

    # Laser Beam Summary: pointing, power and spot round trip; off beams disabled
    summ = export_visrad_summary(bt, joinpath(dir, "summary.csv"))
    b3 = fresh()
    edit_beams!(b3, 1, :power, "5")                       # the file restores it
    import_visrad_summary!(b3, summ)
    @test b3.enabled[] == en
    @test all(((a, b),) -> a.aim ≈ b.aim && a.defocus ≈ b.defocus && a.power ≈ b.power &&
                           a.spot.a ≈ b.spot.a && a.spot.b ≈ b.spot.b,
              zip(ref, b3.beams[]))
    # a file that is not VisRad's leaves the table unchanged
    bad = joinpath(dir, "bad.csv")
    write(bad, "Beam,X\nB115,1\n")
    b4 = fresh()
    before = b4.beams[]
    @test_throws ArgumentError import_visrad_sst!(b4, bad)
    @test b4.beams[] === before
    @test_throws ArgumentError import_visrad_sst!(BeamTable(), sst)       # nothing to apply to
    # loading a preset clears the warnings
    load_preset!(b2, :nif_quads)
    @test isempty(b2.warnings[])
end

@testset "mesh element sizes (no GL)" begin
    ts = TargetSetup(resolution = :coarse)
    @test ts.mesh[] == MESH_PRESETS.coarse && ts.resolution[] === :coarse
    n = nelements(ts.target[])
    @test edit_mesh!(ts, :h, "0.4")
    @test ts.mesh[].h == 0.4 && ts.resolution[] === :custom && nelements(ts.target[]) > n
    @test !edit_mesh!(ts, :capsule_h, "-1") && occursin("positive", ts.error[])
    @test edit_mesh!(ts, :h, "0.6") && ts.resolution[] === :coarse       # back on the preset
    @test_throws ArgumentError edit_mesh!(ts, :n, "1")
    # the hohlraum counts follow the element size
    t = ts.target[]
    wall = elements(t, :wall)
    hmax = maximum(sqrt(t.elem_area[e]) for e in wall)
    @test hmax <= 0.6 * 1.01
    @test set_resolution!(ts, :fine) && ts.mesh[] == MESH_PRESETS.fine
end

@testset "free parts (no GL)" begin
    ts = TargetSetup(resolution = :coarse)
    for (k, kw) in ((:sphere, (; x = 4.0, r1 = 0.5)), (:disk, (; z = 7.0, r1 = 1.0, theta = 180)),
                    (:annulus, (; z = -7.0, r1 = 0.5, r2 = 1.5)),
                    (:rectangle, (; x = -6.0, r1 = 1.0, r2 = 2.0, theta = 90, phi = 0)),
                    (:cylinder, (; y = 5.0, r1 = 0.3, length = 2.0, theta = 90, phi = 90)),
                    (:cone, (; y = -5.0, theta = 90, phi = -90)))
        @test add_part!(ts, k; kw...)
        @test ts.error[] == ""
    end
    t = ts.target[]
    @test [p.name for p in ts.parts[]] == [:sphere_1, :disk_1, :annulus_1, :rectangle_1, :cylinder_1, :cone_1]
    @test all(n -> n in component_names(t), (:sphere_1, :cone_1))
    # placement: a disk tilted to θ = 180° faces −z, at its center
    es = elements(t, :disk_1)
    @test all(e -> t.elem_centroid[e][3] ≈ 7.0, es)
    @test all(e -> t.elem_normal[e][3] < -0.99, es)
    # a rectangle along +x faces +x; its area is width × height
    er = elements(t, :rectangle_1)
    @test sum(t.elem_area[er]) ≈ 2.0 && all(e -> t.elem_normal[e][1] > 0.99, er)
    # the sphere carries its exact surface at its center
    c = t.components[component_index(t, :sphere_1)]
    @test c.surface isa SphereSurface && c.surface.center ≈ Briareus.Vec3(4, 0, 0)
    # edits: accepted, rejected, renamed, flags
    i = part_index(ts, :disk_1)
    @test edit_part!(ts, i, :r1, "1.5") && ts.parts[][i].values[:r1] == 1.5
    @test !edit_part!(ts, i, :r1, "0") && occursin("positive", ts.error[])
    @test !edit_part!(ts, i, :albedo, "2") && occursin("between 0 and 1", ts.error[])
    @test !edit_part!(ts, i, :name, "wall") && occursin("two components", ts.error[])
    @test !edit_part!(ts, i, :name, "1bad") && occursin("letters", ts.error[])
    @test edit_part!(ts, i, :name, "shield") && :shield in component_names(ts.target[])
    @test set_part_flag!(ts, i, :aperture, true)
    @test ts.target[].components[component_index(ts.target[], :shield)].kind === :aperture
    @test set_part_flag!(ts, i, :aperture, false)
    n1 = length(elements(ts.target[], :shield))
    @test set_part_flag!(ts, i, :two_sided, true)
    @test length(elements(ts.target[], :shield)) == 2 * n1
    @test !add_part!(ts, :annulus; r1 = 2.0, r2 = 1.0) && occursin("outer radius", ts.error[])
    @test !add_part!(ts, :torus) && occursin("part kind", ts.error[])
    @test remove_part!(ts, i) && !(:shield in component_names(ts.target[]))
    # parts survive a change of template; a custom target starts with a sphere
    np = length(ts.parts[])
    @test set_target_kind!(ts, :foil) && length(ts.parts[]) == np
    @test set_target_kind!(ts, :custom)
    @test component_names(ts.target[]) == [p.name for p in ts.parts[]]
    ts2 = TargetSetup(kind = :custom, resolution = :coarse)
    @test component_names(ts2.target[]) == [:sphere_1]
    @test !remove_part!(ts2, 1) && occursin("empty", ts2.error[])
end

@testset "model tree and selection (no GL)" begin
    sess = DesignSession(resolution = :coarse)
    rows = model_tree(sess)
    ids = [r.id for r in rows]
    @test ids[1] == "target" && "comp:wall" in ids && "comp:capsule" in ids && "beams" in ids
    @test "group:inner23" in ids && !("quad:Q15T" in ids)          # cones collapsed
    @test count(r -> r.kind == "task", rows) == length(TASKS)
    @test all(r -> r.status == "todo", filter(r -> r.kind == "task", rows))
    toggle_expanded!(sess.view, "group:inner23")
    rows = model_tree(sess)
    kids = [r for r in rows if r.depth == 2]
    @test length(kids) == 8 && all(r -> r.kind == "beam", kids)      # NIF quads: one beam each
    # selection by id: beams under it, components
    select!(sess.view, "group:inner23")
    @test length(selected_beams(sess)) == 8 && selected_component(sess) === nothing
    @test only(r for r in model_tree(sess) if r.id == "group:inner23").selected
    select!(sess.view, kids[1].id)
    @test length(selected_beams(sess)) == 1
    select!(sess.view, "comp:wall")
    @test selected_component(sess) === :wall && isempty(selected_beams(sess))
    add_window!(sess.target)
    @test component_id(sess.target, :window_1) == "window:1"
    select!(sess.view, "window:1")
    @test selected_component(sess) === :window_1
    remove_window!(sess.target, 1)
    @test sess.view.selection[] == "target"                         # the selection went away
    @test selection_kind("part:shield") == (:part, "shield") && selection_kind("") == (:none, "")
    # editing and switching a node by name, whatever the table's grouping
    bt = sess.beams
    set_grouping!(bt, :beam)
    @test edit_beam_node!(bt, :group, "inner30", :z, "3.5")
    @test all(b -> abs(b.aim[3]) ≈ 3.5, filter(b -> b.group === :inner30, bt.beams[]))
    @test bt.by[] === :beam                                          # the table kept its grouping
    q = bt.beams[][1].name
    set_node_enabled!(bt, :beam, q, false)
    @test !bt.enabled[][1] && count(bt.enabled[]) == 47
    @test model_tree(sess)[findfirst(r -> r.id == "beams", model_tree(sess))].on == 2
    @test !edit_beam_node!(bt, :quad, "nope", :z, "1") && occursin("no quad", bt.error[])
    # hidden components
    set_hidden!(sess.view, :wall, true)
    toggle_expanded!(sess.view, "target")
    toggle_expanded!(sess.view, "target")
    @test only(r for r in model_tree(sess) if r.id == "comp:wall").on == 0
end

@testset "beam geometry (no GL)" begin
    sess = DesignSession(resolution = :coarse)
    t = sess.target.target[]
    bs = sess.beams.beams[]
    ps = beam_paths(t, bs)
    # every published NIF quad lands on the wall, past its aim in the LEH plane
    @test all(p -> p.component === :wall, ps)
    @test all(p -> abs(p.stop[3]) < abs(p.aim[3]), ps)
    @test all(p -> abs(hypot(p.stop[1], p.stop[2]) - 2.875) < 0.05, ps)
    # the start is `length` before the aim, against the beam direction
    b = bs[1]
    p = beam_path(t, b; length = 10.0)
    @test p.start ≈ Briareus.Vec3(b.aim) - 10.0 * direction(b)
    # a beam aimed away from the target misses it
    away = remake(b; aim = Briareus.Vec3(b.aim) + 20 * BriareusGUI.normalize(BriareusGUI.cross(direction(b), Briareus.Vec3(0, 0, 1))))
    @test beam_path(t, away).component === :missed
    # volumes: closed rings around the axis, cut at the target
    v, f = beam_volume(t, b; nring = 16, nseg = 4)
    @test length(f) == 2 * 16 * (length(v) ÷ 16 - 1)
    @test all(x -> 1 <= x <= length(v), Iterators.flatten(f))
    @test maximum(x -> hypot(x[1], x[2]), v[2:end]) < 2 * t.scale + 3     # nothing far past the wall
    # pick distance from a ray through the axis is zero; off it, the offset
    @test beam_hit_distance(p, p.start - direction(b), direction(b)) < 1e-9
end

@testset "beam groups, sets and quads (no GL)" begin
    bt = BeamTable(preset_beams(:nif_quads); preset = :nif_quads)
    @test clusters(bt) == [:inner23, :inner30, :outer44, :outer50]
    # a backlighter group: two outer quads aimed at a point away from TCC
    idx = [findfirst(b -> b.group === :outer44, bt.beams[]), findfirst(b -> b.group === :outer50, bt.beams[])]
    @test set_cluster!(bt, idx, "backlighter")
    @test :backlighter in clusters(bt) && bt.cluster[][idx] == [:backlighter, :backlighter]
    @test bt.beams[][idx[1]].group === :outer44                   # the cone stays for LPI
    @test only(r for r in beam_rows(bt) if r.key == "backlighter").count == 2
    for (f, v) in ((:aim_x, "0.5"), (:aim_y, "-1.0"), (:aim_z, "6.0"))
        @test edit_beam_set!(bt, idx, f, v)
    end
    @test all(i -> bt.beams[][i].aim ≈ Briareus.Vec3(0.5, -1.0, 6.0), idx)
    v = beam_set_values(bt, idx)
    @test v.aim ≈ Briareus.Vec3(0.5, -1.0, 6.0) && v.count == 2 && v.enabled === true
    # a relative edit on a set in one hemisphere uses the aim's own z
    up = findall(b -> b.port[3] > 0 && b.group === :inner23, bt.beams[])
    @test !beam_set_values(bt, up).mirrored
    @test edit_beam_set!(bt, up, :z, "2.5") && all(i -> bt.beams[][i].aim[3] ≈ 2.5, up)
    both = findall(==(:inner30), bt.cluster[])
    @test beam_set_values(bt, both).mirrored
    @test edit_beam_set!(bt, both, :z, "3.3")
    @test all(i -> bt.beams[][i].aim[3] ≈ sign(bt.beams[][i].port[3]) * 3.3, both)
    @test !edit_beam_set!(bt, both, :power, "-1") && occursin("negative", bt.error[])
    @test !edit_beam_set!(bt, Int[], :z, "1") && occursin("no beams", bt.error[])
    set_set_enabled!(bt, idx, false)
    @test beam_set_values(bt, idx).enabled === false && count(bt.enabled[]) == 46
    @test !set_cluster!(bt, idx, "  ") && occursin("name", bt.error[])
    @test rename_cluster!(bt, :backlighter, "bl") && :bl in clusters(bt) && !(:backlighter in clusters(bt))
    # the optimizer's per-group parameters follow the groups
    p = Controllers._pointing(bt, (:z,); by = :group)
    @test :bl in p.keys && length(p.keys) == 5

    # split a quad into its four beams, edit one, merge back
    q = bt.beams[][1].quad
    g = bt.cluster[][1]
    aim, P = bt.beams[][1].aim, bt.beams[][1].power
    @test split_quad!(bt, q)
    four = findall(b -> b.quad == q, bt.beams[])
    @test length(four) == 4 && length(bt.beams[]) == 51 && all(==(g), bt.cluster[][four])
    @test all(i -> bt.beams[][i].aim ≈ aim && bt.beams[][i].power ≈ P / 4, four)
    @test all(i -> bt.beams[][i].port != bt.beams[][four[1]].port || i == four[1], four)   # own ports
    @test !split_quad!(bt, q) && occursin("already split", bt.error[])
    @test edit_beam_set!(bt, four[1:2], :z, "2.0")                 # an upper half pointed apart
    @test merge_quad!(bt, q)
    m = findfirst(b -> b.quad == q, bt.beams[])
    @test length(bt.beams[]) == 48 && bt.beams[][m].power ≈ P && bt.beams[][m].fnumber ≈ 7.9
    @test bt.beams[][m].aim[3] ≈ (2 * 2.0 + 2 * aim[3]) / 4
    @test !merge_quad!(bt, q) && occursin("not split", bt.error[])
    @test !split_quad!(bt, "Q99T") && occursin("no quad", bt.error[])
end

@testset "multi-selection and round windows (no GL)" begin
    sess = DesignSession(resolution = :coarse)
    v = sess.view
    select!(v, "group:inner23")
    select_add!(v, "group:inner30")
    @test selected_ids(v) == ["group:inner23", "group:inner30"] && length(selected_beams(sess)) == 16
    rows = model_tree(sess)
    @test count(r -> r.selected, rows) == 2
    @test only(r for r in rows if r.id == "group:inner23").color == cluster_colors(sess.beams)[:inner23]
    @test only(r for r in rows if r.id == "comp:capsule").color == component_colors(sess.target.target[])[:capsule]
    @test only(r for r in rows if r.id == "comp:leh_top").color == ""          # a hole
    select_add!(v, "group:inner23")                                            # toggles off
    @test selected_ids(v) == ["group:inner30"]
    select_add!(v, "comp:wall")                                                # not a beam: replaces
    @test selected_ids(v) == ["comp:wall"]
    # moving the selection to a new group keeps it selectable
    select!(v, "group:outer44")
    set_cluster!(sess.beams, selected_beams(sess), "drive")
    @test v.selection[] == "beams"                                             # outer44 is gone
    @test length(node_beams(sess.beams, "group:drive")) == 16
    # split quads show as quads in their group
    q = sess.beams.beams[][1].quad
    split_quad!(sess.beams, q)
    toggle_expanded!(v, "group:inner23")
    g = String(sess.beams.cluster[][1])
    @test any(r -> r.id == "quad:" * q * "@" * g, model_tree(sess))
    @test length(node_beams(sess.beams, "quad:" * q * "@" * g)) == 4
    # windows are round by default; rectangles on request
    ts = sess.target
    @test add_window!(ts) && ts.windows[][1].round
    t = ts.target[]
    A = sum(t.elem_area[e] for e in elements(t, :window_1))
    @test A ≈ π * 0.5^2 rtol = 1e-6                                           # exact disk area
    @test set_window_round!(ts, 1, false) && !ts.windows[][1].round
    t = ts.target[]
    @test sum(t.elem_area[e] for e in elements(t, :window_1)) ≈ 1.0 rtol = 1e-6
end

@testset "model rows inserted where a tree changes (no GL)" begin
    sess = DesignSession(resolution = :coarse)
    old = model_tree(sess)
    g = only(r.id for r in old if startswith(r.id, "group:inner23"))
    toggle_expanded!(sess.view, g)                       # its quads appear below it, or go
    new = model_tree(sess)
    length(new) < length(old) && ((old, new) = (new, old))
    p = BriareusGUI._common_lead(old, new)
    k = length(new) - length(old)
    @test k > 0 && old[p].id == g && [r.id for r in old[p+1:end]] == [r.id for r in new[p+k+1:end]]
    @test BriareusGUI._common_lead(new, old) == p
end
