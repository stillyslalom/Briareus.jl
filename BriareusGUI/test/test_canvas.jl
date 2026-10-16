# The 3D and map canvases (offscreen GLMakie). Plots are created up front and only
# their inputs change: the plot count never changes.

@testset "scene and map canvases (offscreen)" begin
    sess = DesignSession(resolution = :coarse)
    sc = BriareusGUI.scene_canvas(sess)
    mc = BriareusGUI.map_canvas(sess)
    nplots = length(sc.ax.scene.plots)
    finite(pts) = count(p -> all(isfinite, p), pts)
    t = sess.target.target[]
    opaque = count(t.tri_opaque)

    # Target step: every opaque triangle but the cut quadrant, colored by component
    @test 0 < sc.ntriangles[] < opaque
    @test sc.colorbar.label[] == "component"
    @test sc.colorbar.ticks[][2] == ["wall", "cap_top", "cap_bottom", "capsule"]
    @test length(sc.body.color[]) == 3 * sc.ntriangles[]
    BriareusGUI.set_cutaway!(sc, sess, false)
    @test sc.ntriangles[] == opaque
    BriareusGUI.set_cutaway!(sc, sess, true)
    # one axis segment per beam, from the port side to the aim point
    @test finite(sc.segments[1][]) == 2 * 48
    @test finite(sc.aims[1][]) == 48
    @test !isempty(colorbuffer(sc.fig; px_per_unit = 1))

    # Drive: colored by the field; the map shows the capsule flux
    @test mc.ax.title[] == "Compute the drive to see the flux map"
    build_scenario!(sess)
    set_step!(sess, :drive)
    @test sc.colorbar.label[] == Controllers.DRIVE_FIELD_LABELS[:laser]
    @test sc.colorbar.limits[][1] == 0 && sc.colorbar.limits[][2] > 1e13
    set_field!(sess.drive, :incident)
    @test sc.colorbar.label[] == Controllers.DRIVE_FIELD_LABELS[:incident]
    m = flux_map(sess.drive)
    @test size(mc.heatmap[3][]) == size(m.values)
    @test mc.ax.title[] == "Capsule incident x-ray flux" && mc.ax.yreversed[]
    @test !isempty(colorbuffer(mc.fig; px_per_unit = 1))
    @test !isempty(colorbuffer(sc.fig; px_per_unit = 1))

    # beam edits move the drawn aim; disabled beams stay drawn, dimmed
    set_step!(sess, :beams)
    @test sc.colorbar.label[] == "component"
    edit_beams!(sess.beams, 3, :z, "4.0")
    mem = beam_rows(sess.beams)[3].members
    @test sc.aims[1][][mem[1]][3] ≈ sign(sess.beams.beams[][mem[1]].aim[3]) * 4.0
    set_enabled!(sess.beams, 3, false)
    @test all(c -> c.alpha < 0.5, sc.aims.color[][mem])

    # Review: beams colored by their worst clearance status
    set_enabled!(sess.beams, 3, true)
    bt = sess.beams
    for i in eachindex(bt.enabled[])
        bt.enabled[][i] = i in (1, 2)
    end
    notify(bt.enabled)
    start_clearance!(sess; spawn = false)
    set_step!(sess, :review)
    worst = BriareusGUI._beam_status(sess)
    @test Set(keys(worst)) == Set(b.name for b in active_beams(bt))
    @test sc.aims.color[][1] == BriareusGUI.STATUS_COLORS[worst[bt.beams[][1].name]]

    # a target change: new mesh, map cleared until the drive is computed again
    edit_target!(sess.target, :radius, "3.0")
    @test sess.drive.result[] === nothing
    @test mc.ax.title[] == "Compute the drive to see the flux map"
    set_target_kind!(sess.target, :foil)
    @test sc.colorbar.ticks[][2] == ["backlighter", "main"]
    build_scenario!(sess)
    @test mc.ax.title[] == "Foil laser irradiance" && !mc.ax.yreversed[]
    BriareusGUI.reset_view!(sc)
    @test !isempty(colorbuffer(sc.fig; px_per_unit = 1))
    @test length(sc.ax.scene.plots) == nplots
end
