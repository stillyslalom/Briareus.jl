# Drives the real Briareus window through its QML callbacks — target, parts and
# windows, the model tree and inspector, beams and VisRad files, the 3D view's camera
# and picking, the drive on request, the optimization with the capsule keep-out,
# clearance and TaLIS — and grabs window images with `request_grab` (run by
# runtests.jl in its own process when BRIAREUSGUI_QT_TESTS=true; also runnable alone).
# Exits nonzero on failure.
#
#   julia --project=BriareusGUI -t 4 BriareusGUI/test/qt_window.jl
using Test
using BriareusGUI
using BriareusGUI.Controllers
using BriareusGUI.Briareus

@testset "briareus window (Qt)" begin
    sess = DesignSession(resolution = :coarse)
    stage = Ref(0)
    t_stage = Ref(time())
    t_start = time()
    seen = Symbol[]
    shots = get(ENV, "BRIAREUSGUI_SHOTS", mktempdir())     # keep the images: set the variable
    mkpath(shots)
    files = mktempdir()
    png(name) = joinpath(shots, name * ".png")
    # a grabbed image lands asynchronously: wait for the file, then a moment
    grabbed(path) = isfile(path) && filesize(path) > 0 && time() - t_stage[] > 1
    next!() = (stage[] += 1; t_stage[] = time())
    sst = joinpath(files, "pointing.csv")
    kit = joinpath(files, "kit")
    talis_csv = joinpath(files, "talis.csv")
    session_file = joinpath(files, "design.briareus.toml")
    reference_csv = joinpath(files, "reference.csv")
    retuned_csv = joinpath(files, "retuned.csv")
    target_stl = joinpath(files, "target.stl")
    B = BriareusGUI
    B._TICK_HOOK[] = function (sh)
        s = sh.sess
        time() - t_start > 600 && return B.request_close()
        k = stage[]
        if k == 0 && time() - t_stage[] > 1
            B.br_select("target")
            B.br_edit_target("length", "10.0")                     # through the bridge
            B.br_add_window()                                      # closed, at φ = 90°
            B.br_add_window()                                      # at φ = 270°
            B.br_edit_window(1, "phi", "60")
            B.br_window_open(2, true)
            B.br_edit_window(2, "z", "9")                          # off the wall: rejected
            B.request_grab(png("target")); next!()
        elseif k == 1 && grabbed(png("target"))
            s.target.params[][:length] == 10.0 && push!(seen, :target)
            ws = s.target.windows[]
            length(sh.rows[:windowModel]) == 2 && ws[1].phi == 60 && ws[2].z == 0 && ws[2].open &&
                occursin("past the end", s.target.error[]) && push!(seen, :windows)
            # a shield disk above the top LEH, edited from the inspector
            B.br_add_part("disk")
            name = String(s.target.parts[][end].name)
            B.br_edit_part("z", "7.5")
            B.br_edit_part("r1", "2.5")
            B.br_edit_part("r1", "-1")                              # rejected
            s.view.selection[] == "part:" * name && s.target.parts[][end].values[:z] == 7.5 &&
                s.target.parts[][end].values[:r1] == 2.5 && !isempty(s.target.error[]) &&
                any(r -> r.id == "part:" * name, model_tree(s)) && push!(seen, :part)
            B.br_edit_mesh("h", "0.5")
            s.target.resolution[] === :custom && s.target.mesh[].h == 0.5 && push!(seen, :mesh)
            B.br_set_resolution("coarse")
            B.request_grab(png("part")); next!()
        elseif k == 2 && grabbed(png("part"))
            B.br_remove_selected()
            isempty(s.target.parts[]) && s.view.selection[] == "target" && push!(seen, :remove)
            # beams: a quad from the tree, the table per quad, VisRad round trip
            B.br_toggle_expand("group:inner23")
            B.br_select("group:inner23")
            B.br_edit_node("z", "3.0")
            q = first(b.name for b in s.beams.beams[] if b.group === :inner23)
            B.br_tree_on("beam:" * q, false)
            s.view.selection[] == "group:inner23" &&
                abs(s.beams.beams[][findfirst(b -> b.group === :inner23, s.beams.beams[])].aim[3]) ≈ 3.0 &&
                count(!, s.beams.enabled[]) == 1 && push!(seen, :tree_beams)
            B.br_set_grouping("group")
            B.br_export_sst(sst)
            row = findfirst(r -> r.enabled === true, beam_rows(s.beams))
            z = beam_rows(s.beams)[row].z
            B.br_edit_beam(row, "z", "4.0")
            B.br_import_sst(sst)
            isfile(sst) && beam_rows(s.beams)[row].z ≈ z && count(!, s.beams.enabled[]) == 1 &&
                !isempty(s.beams.warnings[]) && push!(seen, :visrad)
            B.br_tree_on("beams", true)
            B.br_beam_style("volumes")
            B.br_select("beams")
            B.request_grab(png("beams")); next!()
        elseif k == 3 && grabbed(png("beams"))
            all(s.beams.enabled[]) && sh.scene.volumes.visible[] && push!(seen, :volumes)
            # the camera: zoom toward a pixel, orbit, pan, a click that selects the wall
            c = sh.scene.camera
            d0 = B.norm(c.eye - c.lookat)
            mid = (c.viewport[3] / 2, c.viewport[4] / 2)
            B.camera_scroll!(c, mid, 3)
            B.norm(c.eye - c.lookat) < d0 && push!(seen, :zoom)
            B.br_beam_style("rays")
            B.br_select("")
            B.br_standard_view("front")
            o, dvec = B.pixel_ray(c, mid)
            τ, tri = B._pick_surface(sh.scene, s, o, dvec)
            B.camera_press!(c, mid, :left)
            B.camera_release!(c, mid)
            tri != 0 && startswith(s.view.selection[], "comp:") && push!(seen, :pick)
            B.camera_press!(c, mid, :left)
            B.camera_move!(c, (mid[1] + 80, mid[2] + 30))
            B.camera_release!(c, (mid[1] + 80, mid[2] + 30))
            B.camera_press!(c, mid, :right)
            B.camera_move!(c, (mid[1] - 40, mid[2]))
            B.camera_release!(c, (mid[1] - 40, mid[2]))
            B.br_fit_all()
            B.br_set_cutaway(true)
            # color by the laser: computes the drive on request
            B.br_set_color("laser"); next!()
        elseif k == 4 && s.drive.result[] !== nothing && !s.drive.evaluating[]
            s.drive.scenario[] !== nothing && push!(seen, :drive)
            B.br_select("task:drive")
            B.br_write_kit(kit)                                       # on a worker
            B.request_grab(png("drive_scene")); next!()
        elseif k == 5 && grabbed(png("drive_scene")) && !s.drive.kit_running[]
            isfile(joinpath(kit, "README.txt")) && push!(seen, :kit)
            B.br_set_view("map")
            next!()
        elseif k == 6 && time() - t_stage[] > 1
            sh.figure === sh.map.fig && push!(seen, :map)
            B.request_grab(png("drive_map")); next!()
        elseif k == 7 && grabbed(png("drive_map"))
            B.br_set_view("scene")
            B.br_set_color("component")
            B.br_select("task:optimize")
            B.br_opt_edit("iterations", "6")
            B.br_opt_edit("max_evaluations", "40")
            B.br_opt_start(); next!()
        elseif k == 8 && !s.optimize.running[] && s.optimize.result[] !== nothing
            r = s.optimize.result[]
            r.iterations <= 6 && r.evaluations <= 40 && push!(seen, :optimize)
            r.after.capsule_clearance > -0.05 && push!(seen, :capsule)
            B.request_grab(png("optimize")); next!()
        elseif k == 9 && grabbed(png("optimize"))
            B.br_opt_accept()
            s.optimize.result[].applied && push!(seen, :applied)
            # check two quads only (the check takes about a second per beam)
            B.br_set_grouping("quad")
            for i in 1:48
                B.br_set_beam_on(i, i <= 2)
            end
            B.br_select("task:clearance")
            B.br_rv_start(); next!()
        elseif k == 10 && !s.review.running[] && s.review.report[] !== nothing
            length(sh.rows[:clearModel]) == 2 * length(s.review.hardware[]) &&
                s.view.beam_color[] === :clearance && push!(seen, :clearance)
            B.request_grab(png("clearance")); next!()
        elseif k == 11 && grabbed(png("clearance"))
            B.br_select("task:talis")
            B.br_tl_edit("tolerance", "0.4")
            B.br_tl_start(); next!()
        elseif k == 12 && !s.talis.running[] && s.talis.report[] !== nothing
            !isempty(sh.rows[:talisModel]) && !isempty(sh.rows[:ruleModel]) &&
                sh.rows[:noteModel][1].kind == "assumed" && s.talis.report[].settings.tolerance == 0.4 &&
                push!(seen, :talis)
            B.request_grab(png("talis")); next!()
        elseif k == 13 && grabbed(png("talis"))
            B.br_tl_export(talis_csv)
            isfile(talis_csv) && push!(seen, :talis_csv)
            B.br_select("group:outer44")
            B.br_look_along_beam()
            B.br_projection(true)
            B.request_grab(png("along_beam")); next!()
        elseif k == 14 && grabbed(png("along_beam"))
            # a backlighter group from two Ctrl-selected quads, pointed at one point
            B.br_projection(false)
            B.br_tree_on("beams", true)
            q1, q2 = (b.name for b in s.beams.beams[] if b.group === :outer50)
            B.br_select("beam:" * q1)
            B.br_select_add("beam:" * q2)
            B.br_set_cluster("backlighter")
            ok = s.view.selection[] == "group:backlighter" &&
                 length(node_beams(s.beams, "group:backlighter")) == 2
            B.br_edit_node("aim_x", "6.0")
            B.br_edit_node("aim_z", "0.0")
            ok && all(i -> s.beams.beams[][i].aim ≈ Briareus.Vec3(6.0, s.beams.beams[][i].aim[2], 0.0),
                      node_beams(s.beams, "group:backlighter")) && push!(seen, :group)
            B.br_fit_all()
            B.request_grab(png("group")); next!()
        elseif k == 15 && grabbed(png("group"))
            # split a quad (its inspector state comes with the refresh after selecting it)
            q = first(b.name for b in s.beams.beams[] if b.group === :inner23)
            B.br_select("beam:" * q)
            next!()
        elseif k == 16 && time() - t_stage[] > 0.2
            q = first(b.name for b in s.beams.beams[] if b.group === :inner23)
            B.br_split_quad()
            count(b -> b.quad == q, s.beams.beams[]) == 4 && startswith(s.view.selection[], "quad:") &&
                push!(seen, :split)
            B.br_edit_node("defocus", "0.5")
            B.request_grab(png("split")); next!()
        elseif k == 17 && grabbed(png("split"))
            # undo the split (one step: the tick recorded it with the defocus edit),
            # after saving the session with the quad split
            B.br_save_session(session_file)
            B.br_undo(); next!()
        elseif k == 18 && time() - t_stage[] > 0.2
            q = first(b.quad for b in s.beams.beams[] if b.group === :inner23)
            count(b -> b.quad == q, s.beams.beams[]) == 1 && sh.app["canRedo"] && push!(seen, :undo)
            B.br_redo(); next!()
        elseif k == 19 && time() - t_stage[] > 0.2
            q = first(b.quad for b in s.beams.beams[] if b.group === :inner23)
            count(b -> b.quad == q, s.beams.beams[]) == 4 && push!(seen, :redo)
            # a change, then the saved session back
            B.br_edit_target("length", "11.0")
            B.br_open_session(session_file)
            s.target.params[][:length] == 10.0 && count(b -> b.quad == q, s.beams.beams[]) == 4 &&
                sh.session_path == session_file && push!(seen, :session)
            # an x-ray image down the axis, through the top entrance hole (the drive is
            # stale after the load: computed first)
            B.br_standard_view("top")
            B.br_xray_image(); next!()
        elseif k == 20 && s.drive.image[] !== nothing && !s.drive.building[] && !s.drive.evaluating[]
            count(isfinite, s.drive.image[].values) > 100 && sh.figure === sh.map.fig &&
                push!(seen, :xray)
            # the image beside the 3D view: two canvases at once (camera actions keep it)
            B.br_side_view(true)
            B.br_standard_view("iso")
            B.request_grab(png("xray")); next!()
        elseif k == 21 && grabbed(png("xray"))
            sh.figure === sh.scene.fig && sh.side_figure === sh.map.fig && sh.app["sideView"] &&
                s.drive.view[] === :image && push!(seen, :side)
            # back to the scene, a cut half the shell wide, and the view from a port
            B.br_side_view(false)
            B.br_close_image()
            B.br_cut_width(180)
            B.br_look_from_port("80", "225")         # into the cut
            s.view.cut_width[] == 180 && s.view.cutaway[] && sh.figure === sh.scene.fig &&
                startswith(s.status[], "looking from the port") && push!(seen, :port)
            # retune to a reference drive with bundle 31 down, on a worker
            bs = s.beams.beams[]
            pul = cone_fraction_pulses(bs, [0.0, 1.0, 3.0, 4.0], [0.0, 100.0, 300.0, 0.0],
                                       [0.3, 0.3, 0.25, 0.25])
            write_pulse_table(reference_csv, pul)
            s.beams.enabled[] = [nif_bundle(b) != "31" for b in bs]
            B.br_retune(reference_csv, retuned_csv, "quad", "1.5", "0")
            s.drive.retune_running[] && push!(seen, :retune_started)
            B.request_grab(png("cut_port")); next!()
        elseif k == 22 && grabbed(png("cut_port")) && !s.drive.retune_running[]
            isfile(retuned_csv) && occursin("quad powers", s.drive.retune_status[]) &&
                sh.app["retuneStatus"] == s.drive.retune_status[] && push!(seen, :retune)
            B.br_set_target_kind("custom"); next!()                # no template fields
        elseif k == 23 && time() - t_stage[] > 0.5
            s.target.kind[] === :custom && sh.app["targetKind"] == "custom" &&
                isempty(sh.rows[:targetFieldModel]) && push!(seen, :custom)
            # File ▸ Export target as STL; View ▸ Beams ▸ Show beams that are off
            B.br_export_stl(target_stl)
            B.br_show_off(false)
            isfile(target_stl) && !s.view.show_off[] && push!(seen, :stl)
            B.request_close(); next!()
        end
    end
    t0 = time()
    try
        briareus(sess)
    finally
        B._TICK_HOOK[] = nothing
    end
    @info "window session" seconds = round(time() - t0; digits = 1) shots
    @test stage[] == 24
    @test seen == [:target, :windows, :part, :mesh, :remove, :tree_beams, :visrad, :volumes, :zoom,
                   :pick, :drive, :kit, :map, :optimize, :capsule, :applied, :clearance, :talis,
                   :talis_csv, :group, :split, :undo, :redo, :session, :xray, :side, :port,
                   :retune_started, :retune, :custom, :stl]
    @test !sess.spawn[]                              # restored when the window closed
    @test all(isfile ∘ png, ("target", "part", "beams", "drive_scene", "drive_map", "optimize",
                             "clearance", "talis", "along_beam", "group", "split", "cut_port"))
end
