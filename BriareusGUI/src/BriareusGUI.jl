"""
    BriareusGUI

The Briareus window (`briareus`): a laser-target design tool in the manner of VisRad. A
target (a hohlraum with diagnostic windows, a backlighter foil, or free parts), its
beams (with VisRad pointing and beam-summary interchange), and analyses run on them in
any order: the laser and x-ray drive, a pointing optimization, beam clearance and the TaLIS target-configuration rules, and the drive in
time (pulse tables, drive histories and a time cursor on the 3D view). One Qt Quick
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
using LinearAlgebra: dot, norm, normalize, cross, eigen, Symmetric
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
include("controllers/timeline.jl")            # pulses, drive histories, the time cursor
include("controllers/workspace.jl")           # session files, undo and redo

export display_number
export import_visrad_workspace!, export_visrad_workspace, export_results, retune_pulses,
       retune_to_reference!, session_state, restore!, save_session, load_session!, EditHistory, record!, undo!,
       redo!, can_undo, can_redo, reset_history!
export TargetSetup, TARGET_KINDS, TARGET_KIND_LABELS, TARGET_FIELDS, MeshSpec, MESH_PRESETS,
       MESH_RESOLUTIONS, FOIL_NAME, build_target, foil_center, set_target_kind!,
       set_resolution!, edit_mesh!, edit_target!, target_summary, component_names, surface_names,
       target_ready, await_target!, BACKGROUND_TRIANGLES, LEH_CLOSURE_FIELDS, set_leh_closure!,
       is_lip, shown_component, with_lips
export PartSpec, PART_KINDS, PART_KIND_LABELS, PART_FIELDS, PART_SIZE_FIELDS, part_component,
       part_axis, part_center, add_part!, remove_part!, edit_part!, set_part_flag!, part_index,
       PART_ANCHORS, PART_ANCHOR_LABELS, anchor_frame, set_part_anchor!, set_excluded!,
       PART_COORDINATES, PART_POSITION_FIELDS, part_position, set_part_coords!, part_components,
       all_component_names, PART_ROLES, PART_ROLE_LABELS, part_role, set_part_role!,
       export_target_stl, add_stl_parts!, fit_sphere, part_size_unit
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
       set_side_view!, set_cavity!, cavity_rows, cavity_summary, cavity_sealed,
       set_view!, flux_lookup, flux_map, power_rows, element_values, drive_summary,
       BEAM_MODELS, BEAM_MODEL_LABELS, set_beam_model!
export TimeState, DriveHistory, TIME_CHARTS, TIME_CHART_LABELS, set_chart!, load_pulses!,
       use_pulse_shape!, start_history!, cancel_history!, history_stale, set_time_cursor!,
       set_time_step!, time_step, cursor_values, cursor_range, time_label, time_chart
export OptimizeState, OptimizeResult, OBJECTIVES, OBJECTIVE_LABELS, set_free!,
       edit_optimize!, set_objective!, start_optimize!, cancel_optimize!, optimize_progress,
       parameter_rows
export ReviewState, REVIEW_QUALITY, REVIEW_VIEWS, review_batch, default_selection!, set_role!,
       set_review_view!, edit_review!, clearance_kwargs, start_clearance!, cancel_clearance!,
       clearance_rows, worst_rows, review_brief, export_clearance, review_stale
export TaLISState, TALIS_PIVOTS, TALIS_RULE_LABELS, TALIS_TOLERANCE_NOTE, set_talis_role!,
       set_talis_pivot!, edit_talis!, talis_kwargs, start_talis!, cancel_talis!, talis_rows,
       talis_rule_rows, talis_notes, talis_brief, export_talis, talis_stale
export BEAM_STYLES, BeamPath, beam_path, beam_paths, beam_volume, beam_hit_distance
export ViewState, COLOR_MODES, COLOR_LABELS, BEAM_COLORINGS, PALETTE, component_colors,
       DEFAULT_CUT_COMPONENTS, set_cut!, is_cut,
       cluster_colors, TreeRow, select!, select_add!, select_range!, selected_ids, toggle_expanded!, set_hidden!,
       set_beams_hidden!, selected_components,
       selection_kind, component_id, beam_id, model_tree, node_beams, selected_beams,
       selected_component, component_label
export DesignSession, TASKS, TASK_LABELS, task_status, set_color!, drive_stale,
       build_scenario!, accept_optimum!, radiation_size, include_components!,
       RADIATION_MAX_ELEMENTS

end # module Controllers

using .Controllers

export DesignSession, briareus, task_status, set_color!, build_scenario!, start_optimize!,
       cancel_optimize!, accept_optimum!, start_clearance!, cancel_clearance!, start_talis!,
       cancel_talis!
export TargetSetup, set_target_kind!, set_resolution!, edit_mesh!, edit_target!, add_window!,
       add_part!, edit_part!, remove_part!,
       edit_window!, set_window_open!, remove_window!, BeamTable, load_preset!,
       set_grouping!, beam_rows, edit_beams!, set_enabled!, import_beams!, export_beams,
       import_visrad_sst!, import_visrad_summary!, export_visrad_sst, export_visrad_summary,
       set_field!, set_view!, export_clearance, export_talis

include("canvas/cad_camera.jl")
include("canvas/scene_canvas.jl")
include("canvas/map_canvas.jl")
include("canvas/time_canvas.jl")
include("qt/shell.jl")

using PrecompileTools: @setup_workload, @compile_workload

# The window's first frames on its screen type, which the workload cannot render (no GL
# context): the `precompile` signatures dev/gui_first_frame.jl records, each resolved
# against the loaded packages by name. Signatures that no longer resolve (a dependency's
# update renamed a closure) are skipped. Returns how many compiled.
module _PrecompileNames end                      # the loaded packages, by name
function _replay_precompiles(path)
    isfile(path) || return 0
    include_dependency(path)                     # a new recording rebuilds the cache
    names = _PrecompileNames
    for mod in values(Base.loaded_modules)
        isdefined(names, nameof(mod)) || Core.eval(names, :(const $(nameof(mod)) = $mod))
    end
    n = 0
    for line in eachline(path)
        startswith(line, "precompile(") || continue
        n += try
            Core.eval(names, Meta.parse(line)) === true
        catch
            false
        end
    end
    return n
end

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
        # the capsule role's sphere fit (0.65 s to compile at the first CAD capsule)
        Controllers._fit_capsule(part_components(PartSpec(:sphere; name = :s), 0.3), :s)
        # the first undo and redo (0.3 s of compilation otherwise)
        record!(sess)
        edit_target!(sess.target, :length, "10.2")
        record!(sess)
        undo!(sess)
        redo!(sess)
        model_tree(sess)
        for task in TASKS
            task_status(sess, task)
        end
        # a drive on the cavity, its result and summaries (the window's first drive
        # compiled these on the GUI thread)
        set_cavity!(sess.drive, true)
        build_scenario!(sess)
        drive_summary(sess.drive.result[])
        cavity_summary(sess.drive)
        element_values(sess.drive.result[], :incident)
        for task in TASKS
            task_status(sess, task)
        end
        # the time panel: pulses, a short drive history and its charts
        use_pulse_shape!(sess; step = 3.0)
        start_history!(sess)
        for c in TIME_CHARTS
            time_chart(sess, c)
        end
        set_time_cursor!(sess, 6.0)
        sc = scene_canvas(sess)
        mc = map_canvas(sess)
        tc = time_canvas(sess)
        # as in a window (jobs on workers, their results queued for the GUI thread):
        # leaving a part out of the calculations and putting it back (0.3 s of compilation
        # on the GUI thread at the first), and a workspace import's call (its body: one
        # compiled copy for every `done`; this file is missing, so it stops at the read)
        queue = Channel{Any}(Inf)
        sess.deliver[] = f -> put!(queue, f)
        sess.spawn[] = true
        function drain()
            t0 = time()
            while !(isempty(sess.meshing[]) && target_ready(sess.target)) || isready(queue)
                time() - t0 > 60 && break
                isready(queue) ? take!(queue)() : sleep(0.01)
            end
        end
        part = sess.target.parts[][end].name
        include_components!(sess, [part], false)
        drain()
        include_components!(sess, [part], true)
        drain()
        import_visrad_workspace!(sess, joinpath(@__DIR__, "missing.vrw"); done = () -> nothing)
        drain()
        # the window's Julia side, short of GL: its shell, property map and item models,
        # their first refresh, and the QML callbacks (QML's `__init__` has created the Qt
        # application here as well); 14 s of compilation at the first window otherwise
        sh = BriareusShell(DesignSession(resolution = :coarse))
        _refresh!(sh)
        _register_qml_functions()
    end
    precompile(br_tick, ())
    _replay_precompiles(joinpath(@__DIR__, "qt", "first_frame.jl"))
    # an import's beams (0.17 s of its commit on the GUI thread)
    kw = NamedTuple{(:preset, :enabled, :cluster),Tuple{Symbol,Vector{Bool},Vector{Symbol}}}
    precompile(Core.kwcall, (kw, typeof(set_beams!), BeamTable, Vector{Beam}))
end

end # module BriareusGUI
