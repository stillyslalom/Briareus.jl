# The 3D and map canvases (offscreen GLMakie). Plots are created up front and only
# their inputs change: the plot count never changes. The CAD camera's gestures are
# driven through their handlers, as the scene's mouse events would drive them.

const B = BriareusGUI

@testset "scene and map canvases (offscreen)" begin
    sess = DesignSession(resolution = :coarse)
    sc = B.scene_canvas(sess)
    mc = B.map_canvas(sess)
    nplots = length(sc.lscene.scene.plots)
    finite(pts) = count(p -> all(isfinite, p), pts)
    t = sess.target.target[]
    opaque = count(t.tri_opaque)

    # every opaque triangle but the cut quadrant, colored by component as the tree's
    # swatches are; no colorbar for components
    @test 0 < sc.ntriangles[] < opaque
    @test !sc.colorbar.blockscene.visible[]
    @test length(sc.body.color[]) == 3 * sc.ntriangles[]
    pal = component_colors(t)
    @test collect(keys(pal)) ⊆ [:wall, :cap_top, :cap_bottom, :capsule] && length(pal) == 4
    @test B._hexcolor(pal[:capsule]) in sc.body.color[]
    sess.view.cutaway[] = false
    @test sc.ntriangles[] == opaque
    sess.view.cutaway[] = true
    # rays: one segment per beam from its start to where it lands on the wall
    @test finite(sc.rays[1][]) == 2 * 48 && sc.rays.visible[]
    @test finite(sc.aims[1][]) == 48 && finite(sc.hits[1][]) == 48
    @test all(p -> abs(hypot(p[1], p[2]) - 2.875) < 0.1, sc.hits[1][])
    @test !sc.volumes.visible[]
    sess.view.beams[] = :volumes
    @test sc.volumes.visible[] && !sc.rays.visible[] && length(sc.volumes[1][]) > 48 * 24
    sess.view.beams[] = :off
    @test !sc.volumes.visible[] && !sc.rays.visible[] && !sc.aims.visible[]
    sess.view.beams[] = :rays
    @test !isempty(colorbuffer(sc.fig; px_per_unit = 1))

    # drive fields: chosen in the view, computed on request
    @test mc.ax.title[] == "Compute the drive to see the flux map"
    set_color!(sess, :laser)                          # builds the scenario (inline here)
    @test sess.drive.result[] !== nothing
    @test sc.colorbar.label[] == Controllers.DRIVE_FIELD_LABELS[:laser]
    @test sc.colorbar.blockscene.visible[]
    @test sc.colorbar.limits[][1] == 0 && sc.colorbar.limits[][2] > 1e13
    set_color!(sess, :incident)
    @test sc.colorbar.label[] == Controllers.DRIVE_FIELD_LABELS[:incident]
    m = flux_map(sess.drive)
    @test size(mc.heatmap[3][]) == size(m.values)
    @test mc.ax.title[] == "Capsule incident x-ray flux" && mc.ax.yreversed[]
    @test !isempty(colorbuffer(mc.fig; px_per_unit = 1))
    @test !isempty(colorbuffer(sc.fig; px_per_unit = 1))
    set_color!(sess, :component)
    @test !sc.colorbar.blockscene.visible[]

    # beam edits move the drawn aim; disabled beams stay drawn, dimmed
    edit_beams!(sess.beams, 3, :z, "4.0")
    mem = beam_rows(sess.beams)[3].members
    @test sc.aims[1][][mem[1]][3] ≈ sign(sess.beams.beams[][mem[1]].aim[3]) * 4.0
    set_enabled!(sess.beams, 3, false)
    @test all(c -> c.alpha < 0.5, sc.aims.color[][mem])
    set_enabled!(sess.beams, 3, true)

    # beams take their group's color; a new group gets its own
    pal = cluster_colors(sess.beams)
    @test sc.aims.color[][1] == B._hexcolor(pal[sess.beams.cluster[][1]])
    set_cluster!(sess.beams, [1, 2], "backlighter")
    @test sc.aims.color[][1] == B._hexcolor(cluster_colors(sess.beams)[:backlighter])
    @test sc.aims.color[][1] != sc.aims.color[][3]
    set_cluster!(sess.beams, [1, 2], String(sess.beams.beams[][1].group))

    # selection: a component is highlighted, a beam drawn bold
    @test !sc.highlight.visible[]
    select!(sess.view, "comp:capsule")
    @test sc.highlight.visible[] && length(sc.highlight[1][]) == count(==(6), t.tri_comp)   # faces
    # an entrance hole is highlighted with its closing lip
    set_leh_closure!(sess.target, true)
    tl = sess.target.target[]
    select!(sess.view, "comp:leh_top")
    nleh = count(k -> tl.components[k].name in (:leh_top, :leh_top_lip), tl.tri_comp)
    @test sc.highlight.visible[] && length(sc.highlight[1][]) == nleh &&
          nleh > count(k -> tl.components[k].name === :leh_top, tl.tri_comp)
    set_leh_closure!(sess.target, false)
    select!(sess.view, "group:inner23")
    @test !sc.highlight.visible[] && sc.selected_rays.visible[] && finite(sc.selected_rays[1][]) == 2 * 8
    select_add!(sess.view, "group:inner30")
    @test finite(sc.selected_rays[1][]) == 2 * 16
    c0 = B.selection_center(sess)
    @test c0 !== nothing && abs(c0[3]) < 3.5

    # beams colored by their worst clearance status
    bt = sess.beams
    for i in eachindex(bt.enabled[])
        bt.enabled[][i] = i in (1, 2)
    end
    notify(bt.enabled)
    start_clearance!(sess; spawn = false)
    sess.view.beam_color[] = :clearance
    worst = B._beam_status(sess)
    @test Set(keys(worst)) == Set(b.name for b in active_beams(bt))
    @test sc.aims.color[][1] == B.STATUS_COLORS[worst[bt.beams[][1].name]]
    sess.view.beam_color[] = :group
    @test B._beam_status(sess) === nothing

    # a target change: new mesh, map cleared until the drive is computed again
    edit_target!(sess.target, :radius, "3.0")
    @test sess.drive.result[] === nothing
    @test mc.ax.title[] == "Compute the drive to see the flux map"
    set_target_kind!(sess.target, :foil)
    @test Set(keys(component_colors(sess.target.target[]))) == Set([:backlighter, :main])
    build_scenario!(sess)
    @test mc.ax.title[] == "Foil laser irradiance" && !mc.ax.yreversed[]
    # hidden components are not drawn and not picked
    n = sc.ntriangles[]
    set_hidden!(sess.view, :main, true)
    @test sc.ntriangles[] < n
    B.reset_view!(sc)
    @test !isempty(colorbuffer(sc.fig; px_per_unit = 1))
    @test length(sc.lscene.scene.plots) == nplots
end

@testset "CAD camera (offscreen)" begin
    sess = DesignSession(resolution = :coarse)
    sc = B.scene_canvas(sess)
    c = sc.camera
    t = sess.target.target[]
    # the view holds the whole target; near and far planes bracket it
    near, far = B.clip_planes(c)
    C, R = c.bounds
    @test 0 < near < B.norm(c.eye - C) - t.scale && far > B.norm(c.eye - C) + t.scale
    # projection and pixel rays are inverses
    p = B.Makie.Vec3d(0.5, -0.3, 1.2)
    px = B.project_point(c, p)
    o, d = B.pixel_ray(c, px)
    @test B.norm(B.cross(p - o, d)) < 1e-9
    # zooming toward a pixel keeps the point under it fixed, and moves the eye in
    B.orient_camera!(c, :front)
    mid = (c.viewport[3] / 2, c.viewport[4] / 2)
    off = (mid[1] + 60, mid[2] + 40)
    q, hit = B._point_under(c, off)
    @test hit                                          # the wall is under that pixel
    d0 = B.norm(c.eye - q)
    B.camera_scroll!(c, off, 2)
    @test all(isapprox.(B.project_point(c, q), off; atol = 1e-6))
    @test B.norm(c.eye - q) ≈ 0.85^2 * d0
    # deep inside the hohlraum the near plane follows, nothing is clipped
    for _ in 1:20
        B.camera_scroll!(c, off, 1)
    end
    near, far = B.clip_planes(c)
    @test near < 0.2 * B.norm(c.eye - q)
    # orbiting about a center keeps that center under the same pixel
    B.fit_all!(sc)
    ctr = B.Makie.Vec3d(0, 0, 2.0)
    before = B.project_point(c, ctr)
    B.orbit_camera!(c, 0.4, 0.2; center = ctr)
    @test all(isapprox.(B.project_point(c, ctr), before; atol = 1e-6))
    # turntable keeps z up on screen
    @test B.camera_frame(c)[3][3] > 0
    # panning moves a point at the grab depth with the pointer
    g = B.Makie.Vec3d(0, 0, 0)
    pg = B.project_point(c, g)
    depth = B.dot(g - c.eye, B.camera_frame(c)[1])
    B.pan_camera!(c, 25, -10; depth)
    @test all(isapprox.(B.project_point(c, g), (pg[1] + 25, pg[2] - 10); atol = 1e-6))
    # orthographic: same gestures
    B.set_projection!(sc, true)
    q2, _ = B._point_under(c, off)
    B.camera_scroll!(c, off, -1)
    @test all(isapprox.(B.project_point(c, q2), off; atol = 1e-6))
    B.set_projection!(sc, false)

    # clicks select what is under the pointer: the wall in front, a beam near its axis
    B.orient_camera!(c, :front)
    B.select!(sess.view, "")
    B.camera_press!(c, mid, :left)
    B.camera_release!(c, mid)
    @test startswith(sess.view.selection[], "comp:")
    i = 5
    path = sc.paths[][i]
    bpx = B.project_point(c, (path.start + path.aim) / 2)
    c.last_click = -Inf
    B.camera_press!(c, bpx, :left)
    B.camera_release!(c, bpx)
    @test sess.view.selection[] == beam_id(sess.beams.beams[][i])
    j = 9
    path2 = sc.paths[][j]
    bpx2 = B.project_point(c, (path2.start + path2.aim) / 2)
    c.last_click = -Inf
    B.camera_press!(c, bpx2, :left; ctrl = true)
    B.camera_release!(c, bpx2)
    @test selected_beams(sess) == sort([i, j])
    select!(sess.view, beam_id(sess.beams.beams[][i]))
    # a drag is not a click; it orbits about the selection
    sel = sess.view.selection[]
    e0 = c.eye
    B.camera_press!(c, mid, :left)
    B.camera_move!(c, (mid[1] + 50, mid[2]))
    B.camera_release!(c, (mid[1] + 50, mid[2]))
    @test sess.view.selection[] == sel && c.eye != e0
    # double click fits the clicked component
    B.select!(sess.view, "comp:capsule")
    B.fit_selection!(sc, sess)
    @test B.norm(c.lookat) < 0.1                       # the capsule is at TCC
    @test B.norm(c.eye - c.lookat) < 0.5 * t.scale / sind(c.fov / 2)
    # standard views and the view down a beam
    B.standard_view!(sc, :top)
    @test B.camera_frame(c)[1] ≈ B.Makie.Vec3d(0, 0, -1)
    select!(sess.view, "group:outer44")
    B.look_along_beam!(sc, sess)
    b = sess.beams.beams[][first(selected_beams(sess))]
    @test B.camera_frame(c)[1] ≈ B.Makie.Vec3d(direction(b)...)
    # the triad follows the rotation
    @test !isempty(colorbuffer(sc.fig; px_per_unit = 1))
end

@testset "window outlines, parts and TaLIS beam colors (offscreen)" begin
    sess = DesignSession(resolution = :coarse)
    sc = B.scene_canvas(sess)
    nplots = length(sc.lscene.scene.plots)
    finite(pts) = count(p -> all(isfinite, p), pts)
    @test finite(sc.outlines[1][]) == 0
    add_window!(sess.target)
    add_window!(sess.target; phi = 200.0, open = true)
    # round windows: inside and outside circles of 48 pieces each
    @test finite(sc.outlines[1][]) == 2 * 2 * 48 * 2
    R = sess.target.params[][:radius]
    @test all(p -> abs(hypot(p[1], p[2]) - R) < 0.02, filter(p -> all(isfinite, p), sc.outlines[1][]))
    @test sc.outlines.color[][end] == B.WINDOW_COLORS.open
    # the closed window has a color; the open one is a hole
    pal = component_colors(sess.target.target[])
    @test haskey(pal, :window_1) && !haskey(pal, :window_2)
    # a rectangular window: two arcs and two axial edges per loop
    set_window_round!(sess.target, 1, false)
    @test finite(sc.outlines[1][]) == 2 * 48 * 2 + 2 * (2 * 12 + 2) * 2
    set_window_round!(sess.target, 1, true)
    # a part is drawn and selectable
    add_part!(sess.target, :disk; z = 7.0, r1 = 1.0, theta = 180)
    @test haskey(component_colors(sess.target.target[]), :disk_1)
    select!(sess.view, "part:disk_1")
    @test sc.highlight.visible[]
    @test !isempty(colorbuffer(sc.fig; px_per_unit = 1))
    # beams colored by the TaLIS report
    bt = sess.beams
    for i in eachindex(bt.enabled[])
        bt.enabled[][i] = i in (1, 2)
    end
    notify(bt.enabled)
    remove_part!(sess.target, 1)
    start_talis!(sess; spawn = false)
    @test B._beam_status(sess) === nothing                  # colored by cone
    sess.view.beam_color[] = :talis
    worst = B._beam_status(sess)
    @test Set(keys(worst)) == Set(b.name for b in active_beams(bt))
    @test sc.aims.color[][1] == B.STATUS_COLORS[worst[bt.beams[][1].name]]
    remove_window!(sess.target, 2)
    remove_window!(sess.target, 1)
    @test finite(sc.outlines[1][]) == 0
    @test length(sc.lscene.scene.plots) == nplots
end
