"""
    BriareusGUI

The Briareus window (`briareus`): a laser-target design tool in the manner of VisRad. A
target (a hohlraum with diagnostic windows, a backlighter foil, or free parts), its
beams (with VisRad pointing and beam-summary interchange), and analyses run on them in
any order: the laser and x-ray drive (and a VisRad benchmark kit), a pointing
optimization, beam clearance and the TaLIS target-configuration rules. One Qt Quick
window: model tree, properties, a 3D viewer with CAD-style navigation, and result
panels. Use `Controllers` to do the same from a script without a window. Controller
state is exposed through `Observables`.
"""
module BriareusGUI

using Briareus
using GLMakie
using QML: QML, JuliaPropertyMap, JuliaItemModel, loadqml, exec, @qmlfunction
import QMLMakie                       # registers the MakieArea QML type
import Qt6Declarative_jll
import Libdl
using Printf
using LinearAlgebra: norm, dot, cross, normalize

# Framework-free controller layer: the submodule boundary keeps Makie names out of
# scope, so controller code cannot grow GL dependencies by accident.
module Controllers

using Briareus
using Observables
using Printf
using LinearAlgebra: dot, norm, normalize, cross
import TOML

include("controllers/shared.jl")
include("controllers/target_setup.jl")
include("controllers/beam_table.jl")
include("controllers/drive.jl")
include("controllers/optimizer.jl")
include("controllers/review.jl")
include("controllers/talis.jl")
include("controllers/beam_geometry.jl")
include("controllers/view.jl")
include("controllers/session.jl")             # links the parts
include("controllers/workspace.jl")           # session files, undo and redo

export display_number
export import_visrad_workspace!, export_visrad_workspace, export_results, retune_pulses,
       retune_to_reference!, session_state, restore!, save_session, load_session!, EditHistory, record!, undo!,
       redo!, can_undo, can_redo, reset_history!
export TargetSetup, TARGET_KINDS, TARGET_KIND_LABELS, TARGET_FIELDS, MeshSpec, MESH_PRESETS,
       MESH_RESOLUTIONS, FOIL_NAME, build_target, foil_center, set_target_kind!,
       set_resolution!, edit_mesh!, edit_target!, target_summary, component_names, surface_names
export PartSpec, PART_KINDS, PART_KIND_LABELS, PART_FIELDS, PART_SIZE_FIELDS, part_component,
       part_axis, part_center, add_part!, remove_part!, edit_part!, set_part_flag!, part_index,
       PART_ANCHORS, PART_ANCHOR_LABELS, anchor_frame, set_part_anchor!, set_excluded!,
       PART_COORDINATES, PART_POSITION_FIELDS, part_position, set_part_coords!, part_components,
       all_component_names
export WindowSpec, WINDOW_FIELDS, add_window!, remove_window!, edit_window!, set_window_open!,
       set_window_round!,
       window_name
export FACILITIES, FACILITY_LABELS, FACILITY_PRESETS, facility_of, set_facility!, session_tasks
export set_aim_reference!, follow_references!, component_points, aim_at_component!
export BeamTable, BEAM_PRESETS, BEAM_GROUPINGS, BEAM_FIELDS, PRESET_LABELS, preset_beams,
       load_preset!, set_beams!, set_grouping!, beam_rows, edit_beams!, set_enabled!,
       beam_node, edit_beam_node!, set_node_enabled!, SET_FIELDS, beam_set_values,
       edit_beam_set!, set_set_enabled!, clusters, set_cluster!, rename_cluster!, split_quad!,
       merge_quad!,
       active_beams, display_beams, import_beams!, export_beams, beams_summary
export import_visrad_sst!, import_visrad_summary!, export_visrad_sst, export_visrad_summary
export DriveState, DriveResult, DRIVE_FIELDS, DRIVE_FIELD_LABELS, MODE_DEGREES, set_field!, xray_image!,
       set_side_view!,
       set_view!, flux_lookup, flux_map, power_rows, element_values, drive_summary
export OptimizeState, OptimizeResult, OBJECTIVES, OBJECTIVE_LABELS, set_free!,
       edit_optimize!, set_objective!, start_optimize!, cancel_optimize!, optimize_progress,
       parameter_rows
export ReviewState, REVIEW_QUALITY, REVIEW_VIEWS, review_batch, default_selection!, set_role!,
       set_review_view!, edit_review!, clearance_kwargs, start_clearance!, cancel_clearance!,
       clearance_rows, worst_rows, review_brief, export_clearance, review_stale
export TaLISState, TALIS_PIVOTS, TALIS_RULE_LABELS, TALIS_TOLERANCE_CAVEAT, set_talis_role!,
       set_talis_pivot!, edit_talis!, talis_kwargs, start_talis!, cancel_talis!, talis_rows,
       talis_rule_rows, talis_notes, talis_brief, export_talis, talis_stale
export BEAM_STYLES, BeamPath, beam_path, beam_paths, beam_volume, beam_hit_distance
export ViewState, COLOR_MODES, COLOR_LABELS, BEAM_COLORINGS, PALETTE, component_colors,
       DEFAULT_CUT_COMPONENTS, set_cut!, is_cut,
       cluster_colors, TreeRow, select!, select_add!, select_range!, selected_ids, toggle_expanded!, set_hidden!,
       set_beams_hidden!, selected_components,
       selection_kind, component_id, beam_id, model_tree, node_beams, selected_beams,
       selected_component
export DesignSession, TASKS, TASK_LABELS, task_status, set_color!, drive_stale,
       build_scenario!, accept_optimum!, write_visrad_kit!

end # module Controllers

using .Controllers

export DesignSession, briareus, task_status, set_color!, build_scenario!, start_optimize!,
       cancel_optimize!, accept_optimum!, start_clearance!, cancel_clearance!, start_talis!,
       cancel_talis!, write_visrad_kit!
export TargetSetup, set_target_kind!, set_resolution!, edit_mesh!, edit_target!, add_window!,
       add_part!, edit_part!, remove_part!,
       edit_window!, set_window_open!, remove_window!, BeamTable, load_preset!,
       set_grouping!, beam_rows, edit_beams!, set_enabled!, import_beams!, export_beams,
       import_visrad_sst!, import_visrad_summary!, export_visrad_sst, export_visrad_summary,
       set_field!, set_view!, export_clearance, export_talis

include("canvas/cad_camera.jl")
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
        add_part!(sess.target, :disk; z = 7.0)
        model_tree(sess)
        for task in TASKS
            task_status(sess, task)
        end
        sc = scene_canvas(sess)
        mc = map_canvas(sess)
    end
end

end # module BriareusGUI
