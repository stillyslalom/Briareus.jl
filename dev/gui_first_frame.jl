# Records the signatures the Briareus window compiles for its first frames and writes them
# to BriareusGUI/src/qt/first_frame.jl, which the GUI's precompile workload replays.
#
#     julia dev/gui_first_frame.jl
#
# The workload has no GL context, so it cannot render the canvases on the window's screen
# type (`GLMakie.Screen{QMLMakie.QMLWindow}`); without these signatures the first frame
# takes several seconds to compile. A window opens for a few seconds (it needs a display).
# Rerun after updating GLMakie, QMLMakie or Makie, or after changing what the canvases plot:
# signatures that no longer resolve are skipped, but new ones are compiled at run time.
#
# Each run replaces the list. Signatures the precompile cache already holds do not show in
# the trace, so the list is emptied before the recording: the cache is then built again
# without the replay, and the trace holds everything the first frames compile. The old list
# is put back if the recording fails. Signatures with names that change from one build to
# the next are left out (dev/gui_first_frame_filter.jl).
root = dirname(@__DIR__)
gui = joinpath(root, "BriareusGUI")
out = joinpath(gui, "src", "qt", "first_frame.jl")
include(joinpath(@__DIR__, "gui_first_frame_filter.jl"))
trace = tempname() * ".jl"
script = """
using BriareusGUI
const B = BriareusGUI
n = Ref(0)
B._TICK_HOOK[] = sh -> ((n[] += 1) == 180 && B.request_close())
B.briareus()
"""
old = isfile(out) ? read(out, String) : nothing
write_signatures(out, String[])
try
    run(`$(Base.julia_cmd()) --startup-file=no --threads=auto --project=$gui --trace-compile=$trace -e $script`)
catch
    old === nothing ? rm(out; force = true) : write(out, old)
    rethrow()
end
# Signatures of the rendering stack and the GUI; those compiled while loading packages
# (CxxWrap's enums, QML's `__init__`) run before the workload's code is loaded, and those
# naming the recording script's closures do not resolve.
stack = r"GLMakie|QMLMakie|Makie\.|ComputePipeline|GeometryBasics|BriareusGUI"
recorded = filter(l -> startswith(l, "precompile(") && occursin(stack, l) && !occursin("Main.", l),
                  replace.(readlines(trace), r"\s*# recompile$" => ""))
lines = filter_signatures(recorded)
write_signatures(out, lines)
println("wrote $(length(lines)) signatures to ", relpath(out, root),
        " ($(length(unique(recorded)) - length(lines)) with unstable names left out)")
