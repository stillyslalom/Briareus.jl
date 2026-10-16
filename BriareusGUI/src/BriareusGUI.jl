"""
    BriareusGUI

The Briareus window (`briareus`): a target, its beams, the laser and x-ray drive, a
pointing optimization, and a beam-clearance review, in one Qt Quick window with a
3D GLMakie viewer. Use `Controllers` to run the same steps from a script without a
window. Controller state is exposed through `Observables`.
"""
module BriareusGUI

using Briareus
using GLMakie
using QML: QML, JuliaPropertyMap, JuliaItemModel, loadqml, exec, @qmlfunction
import QMLMakie                       # registers the MakieArea QML type
import Qt6Declarative_jll
import Libdl
using Printf

# Framework-free controller layer: the submodule boundary keeps Makie names out of
# scope, so controller code cannot grow GL dependencies by accident.
module Controllers

using Briareus
using Observables
using Printf
using LinearAlgebra: dot, norm, normalize

include("controllers/shared.jl")
include("controllers/target_setup.jl")
include("controllers/beam_table.jl")
include("controllers/drive.jl")
include("controllers/optimizer.jl")
include("controllers/review.jl")
include("controllers/session.jl")             # links the steps

export JobCancelled, display_number
export TargetSetup, TARGET_KINDS, TARGET_FIELDS, MESH_RESOLUTIONS, FOIL_NAME, build_target,
       foil_center, set_target_kind!, set_resolution!, edit_target!, target_summary,
       component_names, surface_names
export BeamTable, BEAM_PRESETS, BEAM_GROUPINGS, BEAM_FIELDS, PRESET_LABELS, preset_beams,
       load_preset!, set_beams!, set_grouping!, beam_rows, edit_beams!, set_enabled!,
       active_beams, display_beams, import_beams!, export_beams, beams_summary
export DriveState, DriveResult, DRIVE_FIELDS, DRIVE_FIELD_LABELS, MODE_DEGREES, set_field!,
       set_view!, flux_lookup, flux_map, power_rows, element_values, drive_summary
export OptimizeState, OptimizeResult, OBJECTIVES, OBJECTIVE_LABELS, set_free!,
       edit_optimize!, set_objective!, start_optimize!, cancel_optimize!, optimize_progress,
       parameter_rows
export ReviewState, REVIEW_QUALITY, REVIEW_CHUNK, default_selection!, set_role!,
       edit_review!, clearance_kwargs, start_clearance!, cancel_clearance!, clearance_rows,
       worst_rows, review_brief, export_clearance, review_stale
export DesignSession, DESIGN_STEPS, set_step!, step_label, step_status, drive_stale,
       build_scenario!, accept_optimum!

end # module Controllers

using .Controllers

export DesignSession, briareus, set_step!, step_status, build_scenario!, start_optimize!,
       cancel_optimize!, accept_optimum!, start_clearance!, cancel_clearance!
export TargetSetup, set_target_kind!, set_resolution!, edit_target!, BeamTable,
       load_preset!, set_grouping!, beam_rows, edit_beams!, set_enabled!, import_beams!,
       export_beams, set_field!, set_view!, export_clearance

include("canvas/scene_canvas.jl")
include("canvas/map_canvas.jl")
include("qt/shell.jl")

using PrecompileTools: @setup_workload, @compile_workload

# Time-to-first-window workload: the controllers on a coarse hohlraum and the
# canvases' Figure construction (no GL context at precompile time, so no rendering).
@setup_workload begin
    # BRIAREUSGUI_NO_WORKLOAD=true skips the workload (faster rebuilds while developing)
    get(ENV, "BRIAREUSGUI_NO_WORKLOAD", "") == "true" || @compile_workload begin
        sess = DesignSession(resolution = :coarse)
        edit_target!(sess.target, :length, "10.0")
        set_grouping!(sess.beams, :quad)
        beam_rows(sess.beams)
        edit_beams!(sess.beams, 1, :z, "3.2")
        set_grouping!(sess.beams, :group)
        for st in DESIGN_STEPS
            set_step!(sess, st)
            step_status(sess, st)
        end
        sc = scene_canvas(sess)
        mc = map_canvas(sess)
    end
end

end # module BriareusGUI
