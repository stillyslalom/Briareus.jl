# Drives the real Briareus window through every step through its QML callbacks and
# grabs window images with `request_grab` (run by runtests.jl in its own process when
# BRIAREUSGUI_QT_TESTS=true; also runnable alone). Exits nonzero on failure.
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
    png(name) = joinpath(shots, name * ".png")
    # a grabbed image lands asynchronously: wait for the file, then a moment
    grabbed(path) = isfile(path) && filesize(path) > 0 && time() - t_stage[] > 1
    next!() = (stage[] += 1; t_stage[] = time())
    B = BriareusGUI
    B._TICK_HOOK[] = function (sh)
        s = sh.sess
        time() - t_start > 600 && return B.request_close()
        k = stage[]
        if k == 0 && time() - t_stage[] > 1
            B.br_edit_target("length", "10.0")                     # through the bridge
            B.br_set_step("target")
            B.request_grab(png("target")); next!()
        elseif k == 1 && grabbed(png("target"))
            s.target.params[][:length] == 10.0 && push!(seen, :target)
            B.br_set_step("beams")
            B.br_set_grouping("quad")
            B.br_edit_beam(1, "z", "3.2")
            B.br_set_beam_on(2, false)
            B.br_set_grouping("group")
            B.request_grab(png("beams")); next!()
        elseif k == 2 && grabbed(png("beams"))
            length(sh.rows[:beamModel]) == 4 && push!(seen, :beams)
            B.br_set_beam_on(1, true)                                # row 1 mixed → all on
            B.br_set_step("drive")
            B.br_build(); next!()
        elseif k == 3 && s.drive.result[] !== nothing && !s.drive.evaluating[]
            push!(seen, :drive)
            B.request_grab(png("drive_scene")); next!()
        elseif k == 4 && grabbed(png("drive_scene"))
            B.br_set_view("map")
            next!()
        elseif k == 5 && time() - t_stage[] > 1
            sh.figure === sh.map.fig && push!(seen, :map)
            B.request_grab(png("drive_map")); next!()
        elseif k == 6 && grabbed(png("drive_map"))
            B.br_set_view("scene")
            B.br_set_step("optimize")
            B.br_opt_edit("iterations", "6")
            B.br_opt_start(); next!()
        elseif k == 7 && !s.optimize.running[] && s.optimize.result[] !== nothing
            push!(seen, :optimize)
            B.request_grab(png("optimize")); next!()
        elseif k == 8 && grabbed(png("optimize"))
            B.br_opt_accept()
            s.optimize.result[].applied && push!(seen, :applied)
            # review two quads only (the check takes about a second per beam)
            B.br_set_grouping("quad")
            for i in 1:48
                B.br_set_beam_on(i, i <= 2)
            end
            B.br_set_step("review")
            B.br_rv_start(); next!()
        elseif k == 9 && !s.review.running[] && s.review.report[] !== nothing
            length(sh.rows[:clearModel]) == 4 && push!(seen, :review)
            B.request_grab(png("review")); next!()
        elseif k == 10 && grabbed(png("review"))
            B.br_reset_view()
            B.br_set_cutaway(false)
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
    @test stage[] == 11
    @test seen == [:target, :beams, :drive, :map, :optimize, :applied, :review]
    @test !sess.spawn[]                              # restored when the window closed
    @test all(isfile ∘ png, ("target", "beams", "drive_scene", "drive_map", "optimize", "review"))
end
