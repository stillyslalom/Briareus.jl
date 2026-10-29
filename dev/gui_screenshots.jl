# Captures the screenshots of the Briareus window shown on the documentation's
# "Desktop application" page (docs/src/gui.md) and writes them to docs/src/assets/gui/.
#
#     julia dev/gui_screenshots.jl                  # all images
#     julia dev/gui_screenshots.jl overview talis   # only the named ones
#
# The script opens the real window (it needs a display; on a headless Linux machine run it
# under `xvfb-run`), drives it through its QML callbacks from `BriareusGUI._TICK_HOOK[]`
# as `BriareusGUI/test/qt_window.jl` does, and grabs the window's contents below the
# toolbar with `request_grab` (QML `grabToImage`, written asynchronously); the Cut
# settings popup and the Help screen are grabbed alone (`request_cut_settings`,
# `request_popup`). Nothing outside the window is captured. Every target and pointing is a
# built-in public one: the NIF hohlraum template with the NIF quads at illustrative example
# cone aims, and the OMEGA backlighter set on a foil; the 1ω and 2ω power ratios of the
# beam-matrix image are illustrative. The run takes a few minutes (the drive on the
# standard mesh, a pointing optimization, a drive history, the TaLIS rules on every quad,
# a placement sweep of the foil); stages wait for their jobs, so a slower machine gives
# the same images. Qt's high-DPI scaling is switched off so that the window has its
# default size (1640 × 1000) in screen pixels on any display scale; on a screen whose free area is
# smaller than about 1820 × 1110 pixels the window opens at 90 % of that area instead.
# Each grab is quantized to a 256-color palette by Pillow through uv (`uv run --python
# 3.12 --with pillow`), which keeps each PNG near 200 kB; without uv the grabs are copied
# as they are. Rerun after a change to the window's layout or canvases.
root = dirname(@__DIR__)
gui = joinpath(root, "BriareusGUI")
out = joinpath(root, "docs", "src", "assets", "gui")
raw = mktempdir()
const ALL = ["optimize", "overview", "xray", "time_history", "talis", "beam_matrix", "display_cut",
             "cut_plane", "cut_settings", "backlighter", "sweep", "help"]
names = isempty(ARGS) ? ALL : ARGS
all(in(ALL), names) || error("images are named among ", join(ALL, ", "))
mkpath(out)

script = """
using BriareusGUI
const B = BriareusGUI
raw, wanted = ARGS[1], split(ARGS[2], ",")
png(name) = joinpath(raw, name * ".png")
stage = Ref(0)
t_stage = Ref(time())
t_start = time()
next!() = (stage[] += 1; t_stage[] = time())
after(dt) = time() - t_stage[] > dt
# a grab lands asynchronously: wait for the file, then a moment for it to be complete
grabbed(name) = !(name in wanted) || (isfile(png(name)) && filesize(png(name)) > 0 && after(1.5))
function grab(sh, name)
    name in wanted || return
    println("grab ", name, " at ", round(Int, time() - t_start), " s: ", sh.sess.status[])
    B.request_grab(png(name))
end
idle(s) = s.drive.result[] !== nothing && !s.drive.building[] && !s.drive.evaluating[] &&
          B.Controllers.target_ready(s.target)
B._TICK_HOOK[] = function (sh)
    s = sh.sess
    time() - t_start > 900 && (println("timed out at stage ", stage[]); return B.request_close())
    k = stage[]
    c = sh.scene.camera
    mid = (c.viewport[3] / 2, c.viewport[4] / 2)
    if k == 0 && after(2)
        # color the cut-away target by the incident x-ray flux, which computes the drive
        B.br_set_cutaway(true)
        B.br_set_color("incident")
        B.br_select("task:drive")
        B.request_tab("drive"); next!()
    elseif k == 1 && idle(s) && after(1)
        # optimize the cone aims (z within ± 2 mm, r) for capsule symmetry, holding the
        # drive level, with the beams' envelopes kept off the capsule and the LEH rims
        B.br_select("task:optimize")
        B.request_tab("optimize")
        B.br_opt_edit("z", "2")
        B.br_opt_keepout(true)
        B.br_opt_edit("iterations", "40")
        B.br_opt_edit("max_evaluations", "300")
        "optimize" in wanted && B.br_opt_start(); next!()
    elseif k == 2 && !s.optimize.running[] && after(1)
        # applied to the beams: the capsule's θ–φ flux map of the optimized drive
        s.optimize.result[] === nothing || B.br_opt_accept()
        B.br_set_view("map"); next!()
    elseif k == 3 && idle(s) && after(1)
        next!()
    elseif k == 4 && after(1.5)
        grab(sh, "optimize"); next!()
    elseif k == 5 && grabbed("optimize")
        # the hero shot: the optimized design, its drive and modes
        B.br_set_view("scene")
        B.br_select("task:drive")
        B.request_tab("drive")
        B.br_standard_view("iso"); next!()
    elseif k == 6 && idle(s) && after(2)
        grab(sh, "overview"); next!()
    elseif k == 7 && grabbed("overview")
        # an x-ray image through the top entrance hole, 20° off the axis
        B.br_set_cutaway(false)
        B.br_look_from_port("20", "0")
        B.camera_scroll!(c, mid, 4)
        B.br_xray_image(); next!()
    elseif k == 8 && s.drive.image[] !== nothing && idle(s)
        # beside the 3D view, the beams drawn as volumes
        B.br_side_view(true)
        B.br_set_cutaway(true)
        B.br_beam_style("volumes")
        B.br_standard_view("iso"); next!()
    elseif k == 9 && after(1.5)
        grab(sh, "xray"); next!()
    elseif k == 10 && grabbed("xray")
        # the drive in time: N221204's pulse shape, the history, the cursor at peak Tr
        B.br_side_view(false)
        B.br_close_image()
        B.br_beam_style("rays")
        B.br_standard_view("iso")
        if "time_history" in wanted
            B.br_time_shape()
            B.br_time_run()
        end
        B.request_tab("time"); next!()
    elseif k == 11 && !s.time.running[] && after(1)
        h = s.time.history[]
        h === nothing || B.br_time_step(argmax(h.Tr))
        B.br_time_chart("drive"); next!()
    elseif k == 12 && after(1.5)
        grab(sh, "time_history"); next!()
    elseif k == 13 && grabbed("time_history")
        # the TaLIS rules on every quad, the beams colored by their worst status
        B.br_time_now()
        B.br_set_color("component")
        B.br_select("task:talis")
        B.request_tab("talis")
        "talis" in wanted && B.br_tl_start(); next!()
    elseif k == 14 && !s.talis.running[] && after(1.5)
        grab(sh, "talis"); next!()
    elseif k == 15 && grabbed("talis")
        # the beam matrix: illustrative 1ω and 2ω power ratios of 5 % each (as in the
        # VisRad workspaces seen so far; the NIF User Guide gives no fraction), and the
        # unconverted light of one inner quad drawn, its properties beside the tree's dots
        B.br_beam_color("group")
        B.br_set_color("component")
        B.br_unconverted_ratio(1, "0.05")
        B.br_unconverted_ratio(2, "0.05")
        B.br_toggle_expand("group:inner23")
        B.br_select("beam:Q15T")
        B.br_beam_dot("beam:Q15T", "draw", 1, true)
        B.br_beam_dot("beam:Q15T", "draw", 2, true)
        B.br_set_cutaway(true)
        B.br_look_from_port("60", "200")
        B.camera_scroll!(c, mid, -6); next!()
    elseif k == 16 && idle(s) && after(2)
        grab(sh, "beam_matrix"); next!()
    elseif k == 17 && grabbed("beam_matrix")
        # display modes and an exact sector cut: the top cap wireframe, the bottom cap
        # semi-transparent, and a 120° slice of the wall bounded to |z| ≤ 3 mm, toward the view
        B.br_beam_dot("beam:Q15T", "draw", 1, false)
        B.br_beam_dot("beam:Q15T", "draw", 2, false)
        B.br_unconverted_ratio(1, "")
        B.br_unconverted_ratio(2, "")
        B.br_toggle_expand("group:inner23")
        B.br_set_color("incident")
        B.br_set_display("comp:cap_bottom", "transparent")
        B.br_select("comp:cap_top")
        B.br_set_display("comp:cap_top", "wireframe")
        B.br_standard_view("iso")
        B.br_cut_kind("sector")
        B.br_cut_width(120)
        B.br_cut_toward_view()
        B.br_cut_range("-3", "3")
        B.request_tab("drive"); next!()
    elseif k == 18 && idle(s) && after(2)
        grab(sh, "display_cut"); next!()
    elseif k == 19 && grabbed("display_cut")
        # a 2 × 2 mm diagnostic window in the wall at φ = 90°, cut by the plane that
        # contains its line of sight (from the window's centroid along +y)
        B.br_set_display("comp:cap_top", "opaque")
        B.br_set_display("comp:cap_bottom", "opaque")
        B.br_set_color("component")
        B.br_add_window()
        B.br_edit_window(1, "width", "2")
        B.br_edit_window(1, "height", "2"); next!()
    elseif k == 20 && B.Controllers.target_ready(s.target) && after(1)
        # the drive and the TaLIS rules on the target with its window, so nothing is out of date
        B.br_build()
        "cut_plane" in wanted && B.br_tl_start()
        B.br_select("window:1")
        B.br_cut_line("90", "90", true, 0)
        B.br_look_from_port("70", "135"); next!()
    elseif k == 21 && idle(s) && !s.talis.running[] && after(2)
        B.br_beam_color("group"); next!()
    elseif k == 22 && after(1.5)
        grab(sh, "cut_plane"); next!()
    elseif k == 23 && grabbed("cut_plane")
        # the Cut settings popup of that plane
        "cut_settings" in wanted && B.request_cut_settings(png("cut_settings")); next!()
    elseif k == 24 && grabbed("cut_settings")
        # a backlighter foil with the OMEGA backlighter beams, colored by the absorbed
        # laser intensity: the 3D view toward the foil's front, the foil's face beside it.
        # The hohlraum's optimization does not apply to the foil: it is cleared, as in a
        # window opened on the foil.
        B.br_remove_window(1)
        B.br_set_cutaway(false)
        if "backlighter" in wanted || "sweep" in wanted
            s.optimize.result[] = nothing
            s.optimize.status[] = ""
            B.br_set_target_kind("foil")
            B.br_set_color("laser")
        end
        B.br_select("target")
        B.request_tab("drive"); next!()
    elseif k == 25 && (!("backlighter" in wanted || "sweep" in wanted) || idle(s)) && after(2)
        B.br_look_from_port("70", "160")
        B.br_set_view("map")
        B.br_side_view(true); next!()
    elseif k == 26 && after(1.5)
        grab(sh, "backlighter"); next!()
    elseif k == 27 && grabbed("backlighter")
        # a placement sweep of the foil along its front normal, from 2 mm away from the
        # target chamber center to 2 mm toward it: clearance, the line of sight through
        # the main target and the power the foil captures at each step
        B.br_set_view("scene")
        B.br_side_view(false)
        B.br_standard_view("iso")
        B.br_select("task:sweep")
        B.br_sw_check("clearance", true)
        B.br_sw_check("los", true)
        B.br_sw_check("capture", true)
        B.br_sw_edit("start", "-2")
        B.br_sw_edit("stop", "2")
        B.br_sw_edit("step", "1")
        "sweep" in wanted && B.br_sw_start()
        B.request_tab("sweep"); next!()
    elseif k == 28 && !s.sweep.running[] && after(1)
        s.sweep.result[] === nothing || B.br_sw_select(4); next!()
    elseif k == 29 && after(1.5)
        grab(sh, "sweep"); next!()
    elseif k == 30 && grabbed("sweep")
        # the Help screen, at the placement sweep's page
        "help" in wanted && B.request_popup("topics:sweep.intro", png("help")); next!()
    elseif k == 31 && grabbed("help")
        B.request_popup("")
        B.request_close(); next!()
    end
end
briareus(resolution = :standard)
"""
cmd = `$(Base.julia_cmd()) --startup-file=no --threads=auto --project=$gui -e $script $raw $(join(names, ","))`
run(addenv(cmd, "QT_ENABLE_HIGHDPI_SCALING" => "0"))

grabs = filter(isfile, [joinpath(raw, n * ".png") for n in names])
length(grabs) == length(names) || @warn "missing grabs" missing = setdiff(names, first.(splitext.(basename.(grabs))))
shrink = """
import sys
from PIL import Image
files = sys.argv[1:]
for i in range(0, len(files), 2):
    im = Image.open(files[i]).convert("RGB")
    im.quantize(256, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.NONE).save(files[i + 1], optimize=True)
"""
pairs = collect(Iterators.flatten((g, joinpath(out, basename(g))) for g in grabs))
uv = Sys.which("uv")
if uv === nothing
    @warn "uv not found: copying the grabs without quantizing them"
    foreach(g -> cp(g, joinpath(out, basename(g)); force = true), grabs)
else
    run(`$uv run --quiet --python 3.12 --with pillow python -c $shrink $pairs`)
end
for g in grabs
    f = joinpath(out, basename(g))
    println(rpad(relpath(f, root), 40), round(Int, filesize(f) / 1024), " kB")
end
