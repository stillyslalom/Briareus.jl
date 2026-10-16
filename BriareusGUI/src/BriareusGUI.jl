"""
    BriareusGUI

The Briareus window (`briareus`): a target (with diagnostic windows), its beams (with
VisRad pointing and beam-summary interchange), the laser and x-ray drive (and a VisRad
benchmark kit), a pointing optimization, and a review of beam clearance and the TaLIS
target-configuration rules, in one Qt Quick window with a 3D GLMakie viewer. Use `Controllers` to run the same steps from a script without a
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
include("controllers/talis.jl")
include("controllers/session.jl")             # links the steps

export display_number
export TargetSetup, TARGET_KINDS, TARGET_FIELDS, MESH_RESOLUTIONS, FOIL_NAME, build_target,
       foil_center, set_target_kind!, set_resolution!, edit_target!, target_summary,
       component_names, surface_names
export WindowSpec, WINDOW_FIELDS, add_window!, remove_window!, edit_window!, set_window_open!,
       window_name
export BeamTable, BEAM_PRESETS, BEAM_GROUPINGS, BEAM_FIELDS, PRESET_LABELS, preset_beams,
       load_preset!, set_beams!, set_grouping!, beam_rows, edit_beams!, set_enabled!,
       active_beams, display_beams, import_beams!, export_beams, beams_summary
export import_visrad_sst!, import_visrad_summary!, export_visrad_sst, export_visrad_summary
export DriveState, DriveResult, DRIVE_FIELDS, DRIVE_FIELD_LABELS, MODE_DEGREES, set_field!,
       set_view!, flux_lookup, flux_map, power_rows, element_values, drive_summary
export OptimizeState, OptimizeResult, OBJECTIVES, OBJECTIVE_LABELS, set_free!,
       edit_optimize!, set_objective!, start_optimize!, cancel_optimize!, optimize_progress,
       parameter_rows
export ReviewState, REVIEW_QUALITY, REVIEW_VIEWS, review_batch, default_selection!, set_role!,
       set_review_view!, edit_review!, clearance_kwargs, start_clearance!, cancel_clearance!,
       clearance_rows, worst_rows, review_brief, export_clearance, review_stale
export TaLISState, TALIS_PIVOTS, TALIS_RULE_LABELS, TALIS_TOLERANCE_CAVEAT, set_talis_role!,
       set_talis_pivot!, edit_talis!, talis_kwargs, start_talis!, cancel_talis!, talis_rows,
       talis_requirement, talis_rule_rows, talis_notes, talis_brief, export_talis, talis_stale
export DesignSession, DESIGN_STEPS, set_step!, step_label, step_status, drive_stale,
       build_scenario!, accept_optimum!, write_visrad_kit!

end # module Controllers

using .Controllers

export DesignSession, briareus, set_step!, step_status, build_scenario!, start_optimize!,
       cancel_optimize!, accept_optimum!, start_clearance!, cancel_clearance!, start_talis!,
       cancel_talis!, write_visrad_kit!
export TargetSetup, set_target_kind!, set_resolution!, edit_target!, add_window!,
       edit_window!, set_window_open!, remove_window!, BeamTable, load_preset!,
       set_grouping!, beam_rows, edit_beams!, set_enabled!, import_beams!, export_beams,
       import_visrad_sst!, import_visrad_summary!, export_visrad_sst, export_visrad_summary,
       set_field!, set_view!, export_clearance, export_talis

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
        add_window!(sess.target)
        edit_window!(sess.target, 1, :z, "0.5")
        remove_window!(sess.target, 1)
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
