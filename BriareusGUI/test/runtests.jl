using Test
using BriareusGUI
using BriareusGUI.Controllers
using BriareusGUI.GLMakie
using BriareusGUI.Briareus
using SHA: sha256
using ZipArchives: ZipReader, zip_names, zip_readentry

# Drain a window-style delivery queue until `done()` (background-job tests).
function drain_until(done, queue; timeout = 300)
    t0 = time()
    while !done()
        time() - t0 > timeout && error("timed out waiting for a background job")
        while isready(queue)
            take!(queue)()
        end
        sleep(0.02)
    end
    while isready(queue)
        take!(queue)()
    end
    return true
end

@testset "BriareusGUI.jl" begin
    @testset "Controllers are framework-free" begin
        for name in (:Figure, :Axis, :Axis3, :mesh!, :GLMakie, :Makie, :QML, :QMLMakie)
            @test !isdefined(BriareusGUI.Controllers, name)
        end
    end

    include("test_controllers.jl")      # Target and Beams steps (no GL)
    include("test_session.jl")          # Drive, Optimize, Review, background jobs (no GL)
    include("test_talis.jl")            # Review step: TaLIS checks (no GL)
    include("test_workspace.jl")        # session files, undo and redo (no GL)
    include("test_frames.jl")           # frames, aims in frames, the placement sweep (no GL)
    include("test_display.jl")          # participation, visibility, cutaway (no GL)
    include("test_time.jl")             # pulses, drive histories, the time cursor (no GL)
    include("test_app.jl")              # the `briareus` command: shortcuts, script, splash (no GL)
    include("test_canvas.jl")           # 3D and map canvases (offscreen GL)
end

if get(ENV, "BRIAREUSGUI_QT_TESTS", "") == "true"
    @testset "briareus window (Qt, separate process)" begin
        script = joinpath(@__DIR__, "qt_window.jl")
        cmd = `$(Base.julia_cmd()) --project=$(Base.active_project()) -t 4 $script`
        @test success(pipeline(cmd; stdout, stderr))
    end
else
    @info "Skipping the Qt window test; set BRIAREUSGUI_QT_TESTS=true on a machine with a display to run it."
end
