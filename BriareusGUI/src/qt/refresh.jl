# This file defines the shell of the window and the refresh that sets what QML shows from the
# session, one function per panel.

# ---------------------------------------------------------------- shell

const _MODEL_ROWS = (treeModel = TreeItem, partFieldModel = FieldRow,
                     stepModel = StepRow,
                     windowModel = WindowRow, beamModel = BeamRow, modeModel = ModeRow,
                     powerModel = PowerRow, optModeModel = ModeRow, paramModel = ParamRow,
                     compModel = CompRow, clearModel = ClearRow, worstModel = ClearRow,
                     talisCompModel = CompRow, ruleModel = RuleRow, talisModel = TalisRow,
                     noteModel = NoteRow, timeLegendModel = LegendRow, sweepModel = SweepRowItem,
                     sweepLegendModel = LegendRow, diagLegendModel = LegendRow, diagModel = DiagRow,
                     folderModel = FolderRow,
                     issueModel = IssueRow, diagImageModel = ImageRow)

"""
State shared between a `DesignSession` and its QML window: the property map `app`,
the item models (by their QML names) and their rows, the 3D and map canvases, and
the queue of background updates that `br_tick` drains.
"""
mutable struct BriareusShell
    sess::DesignSession
    scene::SceneCanvas
    map::MapCanvas
    time::TimeCanvas
    sweep::TimeCanvas
    diag::TimeCanvas                               # the Diagnostics tab's chart
    diag_image::MapCanvas                          # and its image, on the tab's one canvas in turn
    app::JuliaPropertyMap
    rows::Dict{Symbol,Vector}
    models::Dict{Symbol,JuliaItemModel}
    queue::Channel{Any}
    dirty::Base.RefValue{UInt16}                   # the domains to refresh (dirty.jl)
    shown::Dict{String,Any}
    figure::Any                                    # the figure on the canvas
    session_path::String                           # the session file, once saved or opened
    side_figure::Any                               # the figure on the side canvas
    blank::Any                                     # what the side canvas shows when hidden
end

function BriareusShell(sess::DesignSession; queue::Channel{Any} = Channel{Any}(Inf))
    rows = Dict{Symbol,Vector}(k => R[] for (k, R) in pairs(_MODEL_ROWS))
    models = Dict{Symbol,JuliaItemModel}(k => JuliaItemModel(rows[k]) for k in keys(_MODEL_ROWS))
    scene = scene_canvas(sess)
    sh = BriareusShell(sess, scene, map_canvas(sess), time_canvas(sess), sweep_canvas(sess),
                       diagnostics_canvas(sess), image_canvas(),
                       JuliaPropertyMap(), rows, models, queue, Ref(_DIRTY_ALL), Dict{String,Any}(),
                       scene.fig, "", nothing, Figure(size = (40, 40)))
    sh.app["canvas"] = scene.fig
    sh.app["timeCanvas"] = sh.time.fig
    sh.app["sweepCanvas"] = sh.sweep.fig
    # the Diagnostics tab's one canvas shows the chart or the image, by category
    sh.app["diagCanvas"] = sh.diag.fig
    sh.shown["_diagFigure"] = sh.diag.fig
    _set!(sh, "diagCategory", "xray_power")
    _set!(sh, "pulseShapeKeys", [String(k) for k in PULSE_SHAPES])
    _set!(sh, "pulseShapeLabels", [PULSE_SHAPE_LABELS[k] for k in PULSE_SHAPES])
    _set!(sh, "facilityKeys", [String(f) for f in FACILITIES])
    _set!(sh, "facilityLabels", [FACILITY_LABELS[f] for f in FACILITIES])
    # a figure is shown on one canvas at a time: the side canvas holds a blank one
    # unless the map sits beside the 3D view
    sh.side_figure = sh.blank
    sh.app["sideCanvas"] = sh.blank
    _set!(sh, "sideView", false)
    _set!(sh, "grabPath", "")
    _set!(sh, "resultsRequest", "")
    _set!(sh, "sizeRequest", "")
    _set!(sh, "popupRequest", "")
    _set!(sh, "popupGrabPath", "")
    _set!(sh, "helpRequest", "")
    _set!(sh, "sessionPath", "")
    _set!(sh, "canUndo", false)
    _set!(sh, "canRedo", false)
    _set!(sh, "visradSource", "")
    _set!(sh, "cadPrompt", "")
    for k in ("compInfo", "selTitle", "partKind", "partName", "partFile", "bnBy", "bn_z", "bn_r", "bn_phi",
              "bn_defocus", "bn_power", "bn_aim_x", "bn_aim_y", "bn_aim_z", "bnInfo",
              "bnCluster", "bnQuad")
        _set!(sh, k, "")
    end
    _set!(sh, "bnMirrored", false)
    _set!(sh, "bnSplit", 0)
    _set!(sh, "bnMerge", 0)
    _set!(sh, "bnPowerUnit", "beam")
    _set!(sh, "bnCount", 0)
    _set!(sh, "clusterNames", String[])
    _set!(sh, "aimChoices", ["Target chamber (TCC)"])
    _set!(sh, "bnAimFrame", 0)
    _set!(sh, "bnAimFrameLabel", "")
    _set!(sh, "beamNames", String[])
    _set!(sh, "cutSettingsRequest", false)
    _set!(sh, "cutGrabPath", "")
    _set!(sh, "facility", "nif")
    for k in ("partTwoSided", "partInward", "partAperture", "partCopy", "partSourceExists", "partLehClosure", "partCapsuleMesh")
        _set!(sh, k, false)
    end
    _set!(sh, "bnOn", 0)
    _set!(sh, "bnMatrix", "")
    _set!(sh, "compDisplay", 0)
    _set!(sh, "displayLabels", [DISPLAY_LABELS[m] for m in DISPLAY_MODES])
    _set!(sh, "partNumber", 0)
    _set!(sh, "partTopEnd", 1)
    _set!(sh, "partBottomEnd", 1)
    _set!(sh, "partProfile", 0)
    _set!(sh, "partSectionZ", String[])
    _set!(sh, "partSectionR", String[])
    _set!(sh, "partFrame", 0)
    _set!(sh, "moveFrames", String[])
    _set!(sh, "moveFrame", 0)
    _set!(sh, "moveAxis", 2)
    _set!(sh, "moveStep", "0.1")
    sh.scene.triad_step[] = 0.1
    _set!(sh, "partRole", 0)
    _set!(sh, "partFilm", false)
    _set!(sh, "cavitySummary", "")
    _set!(sh, "compIncluded", true)
    _set!(sh, "roleLabels", [PART_ROLE_LABELS[r] for r in PART_ROLES])
    _set!(sh, "pieceRoleLabels", [PIECE_ROLE_LABELS[r] for r in PIECE_ROLES])
    _set!(sh, "pieceRoleable", false)
    _set!(sh, "pieceRole", 0)
    _set!(sh, "pieceRoleAllowed", fill(false, length(PIECE_ROLES)))
    _set!(sh, "partCoords", 0)
    _set!(sh, "frameLabels", [FRAME_LABELS[:tcc]])
    _set!(sh, "stepKindLabels", [PART_STEP_LABELS[k] for k in PART_STEP_KINDS])
    _set!(sh, "placementNote", "")
    _set!(sh, "windowNumber", 0)
    _set!(sh, "hasHohlraum", false)
    _set!(sh, "framePort", "")
    _set!(sh, "frameTriad", false)
    _set!(sh, "frameInfo", "")
    _set!(sh, "assemblyInfo", "")
    _set!(sh, "assemblyChoices", ["Top level"])
    _set!(sh, "portChoices", String[])
    _set!(sh, "frameCount", 0)
    _set!(sh, "tfFrames", [FRAME_LABELS[:tcc]])
    for k in ("tfFrom", "tfTo", "tfInput", "tfOutput")
        _set!(sh, k, 0)
    end
    _set!(sh, "tfDirection", false)
    _set!(sh, "tfIn", ["0", "0", "0"])
    _set!(sh, "tfOut", ["0", "0", "0"])
    _set!(sh, "tfError", "")
    _set!(sh, "ortho", false)
    _set!(sh, "turntable", true)
    for (obs, mask) in _dirty_sources(sess)
        on(_ -> _mark!(sh.dirty, mask), obs)
    end
    _refresh!(sh)
    return sh
end

# The observables of a session that the window shows, each with the mask of the domains a
# change of it marks. Most mark every domain: only an observable whose readers are known
# marks fewer, and a new one belongs in the first list until its readers are checked
# (the test "dirty domains" compares a refresh of the marked domains with a full one).
function _dirty_sources(sess::DesignSession)
    ts, bt, dr, st, rv, tl, v, tm, sw = sess.target, sess.beams, sess.drive, sess.optimize, sess.review,
                                        sess.talis, sess.view, sess.time, sess.sweep
    cm = sess.corridor
    all = Any[sess.roe.report, sess.roe.status, sess.roe.picked, sess.roe.outage, sess.status, sess.meshing, sess.cad_changed, ts.mesh, ts.resolution,
              ts.parts, ts.target, ts.error, ts.pending, ts.building, ts.fitting, ts.active_foil,
              bt.beams, bt.enabled, bt.preset, bt.by, bt.path, bt.error, bt.warnings,
              dr.scenario, dr.building, dr.result, dr.evaluating, dr.field, dr.view,
              dr.cavity, dr.cavity_point, dr.films, dr.beam_model, dr.profile,
              dr.spot_radius, dr.status, dr.retune_running, dr.unconverted_ratio,
              dr.retune_status, dr.wall_cavity,
              st.by, st.fields, st.deltas, st.objective, st.weights, st.hold_drive, st.keepout,
              st.keepout_min, st.capsule_keepout, st.capsule_min, st.spill, st.uniformity, st.iterations,
              st.max_evaluations, st.running, st.status,
              st.error, st.result, bt.cluster, bt.aim_frame, ts.frames, v.extra,
              sw.subject, sw.frame, sw.axis, sw.start, sw.stop,
              sw.step, sw.clearance, sw.los, sw.los_theta,
              sw.los_phi, sw.los_object, sw.capture, sw.optimize,
              sw.quality, sw.running, sw.status,
              sw.error, sw.result, sw.selected,
              rv.targets, rv.hardware, rv.peak_fraction, rv.wings, rv.errors, rv.quality,
              rv.running, rv.status, rv.error, rv.report, cm.radius, cm.shown, cm.map, cm.status,
              cm.surface, cm.point, cm.normal, cm.extent, cm.error,
              tl.capture, tl.hardware, tl.diagnostics, tl.tolerance, tl.pivot, tl.pivot_point, tl.omega1,
              tl.running, tl.status, tl.error, tl.report, tl.tight, tl.hide_curved, tl.focus, tl.pinch,
              v.display, v.selection, v.hidden_beams, v.shown_1w, v.shown_2w, bt.unconverted,
              v.triads, tm.pulses, tm.history, tm.retune, tm.diagnostics,
              tm.catalog, tm.issues, tm.running,
              tm.status, tm.prepulse, tm.threshold, tm.irradiance, tm.checking]
    # the viewer's options, which the View menu and the cut settings show
    view = Any[v.beams, v.beam_color, v.show_off, v.omega1_contour, v.cut_components, v.beam_length, v.cutaway,
               v.cut_kind, v.cut_axis, v.cut_range, v.cut_point, v.cut_normal, v.cut_shift,
               v.cut_line, v.cut_rotation, v.cut_azimuth, v.cut_width, v.camera]
    # the progress of a job: its panel, its row of the model tree and the list of running jobs
    progress(d) = _dirty_mask(d, :tree, :view)
    return Pair{Any,UInt16}[[o => _DIRTY_ALL for o in all];
                            [o => _dirty_mask(:view) for o in view];
                            # the coloring: the viewer's options, and the mode of the Time tab
                            v.color => _dirty_mask(:view, :time, :drive);
                            tm.cursor => _dirty_mask(:time);
                            tm.chart => _dirty_mask(:time);
                            v.beam_sort => _dirty_mask(:beams);
                            v.expanded => _dirty_mask(:tree);
                            st.evaluations => progress(:optimize);
                            st.iteration => progress(:optimize);
                            rv.progress => progress(:review);
                            sw.progress => progress(:sweep);
                            tm.progress => progress(:time)]
end

# Assign a QML property only when its value changed.
function _set!(sh::BriareusShell, key::String, val)
    haskey(sh.shown, key) && isequal(sh.shown[key], val) && return
    sh.shown[key] = val
    sh.app[key] = val
    return
end

_num(x::Real) = isfinite(x) ? Controllers.display_number(x) : ""
_fix(x::Real, d::Integer = 3) = isfinite(x) ? @sprintf("%.*f", d, x) : "–"
# A fraction as a percentage; one that rounds to zero prints without a sign ("-0.00%").
function _pct(x::Real, d::Integer = 2)
    isfinite(x) || return "–"
    t = @sprintf("%.*f%%", d, 100x)
    return startswith(t, "-") && all(c -> c in "0.", t[2:end-2]) ? t[2:end] : t
end
# A laser power (TW) and its fraction of the delivered power: four significant digits
# and two decimals of a percent; a power below 1 GW as "< 0.001" and a fraction below
# 0.005% as "< 0.01%" (when the beams reach only the wall, the other parts sum to
# round-off, ~10⁻¹³ TW).
_power_text(p::Real) = abs(p) < 1e-3 ? "< 0.001" : @sprintf("%.4g", p)
_power_pct(x::Real) = isfinite(x) && abs(x) < 5e-5 ? "< 0.01%" : _pct(x)

"""
    _refresh!(sh::BriareusShell, mask = _DIRTY_ALL)

Set what QML shows from the session: the title, the status and the undo state, and the
panels of the domains `mask` holds ([`DIRTY_DOMAINS`](@ref)); every panel by default.
"""
function _refresh!(sh::BriareusShell, mask::UInt16 = _DIRTY_ALL)
    sess = sh.sess
    dirty(d) = _is_dirty(mask, d)
    file = isempty(sh.session_path) ? "" : basename(sh.session_path) * " | "
    _set!(sh, "title", "Briareus | " * file * target_summary(sess.target))
    _set!(sh, "sessionPath", sh.session_path)
    _set!(sh, "canUndo", can_undo(sess))
    _set!(sh, "canRedo", can_redo(sess))
    _set!(sh, "visradSource", basename(cad_source(sess, sess.visrad_source[])))
    changed = sess.cad_changed[]
    _set!(sh, "cadPrompt", isempty(changed) ? "" :
                           "These CAD files changed on disk since the session was saved:\n\n" *
                           join(cad_source.(Ref(sess), changed), "\n") *
                           "\n\nThe design reads the copies saved in the session. Read the changed files instead?")
    _set!(sh, "status", sess.status[])
    _set!(sh, "ortho", sh.scene.camera.ortho)             # (a session's camera sets them too)
    _set!(sh, "turntable", sh.scene.camera.turntable)
    dirty(:view) && _refresh_view!(sh)
    dirty(:tree) && _refresh_target!(sh)
    dirty(:beams) && _refresh_beams!(sh)
    dirty(:selection) && _refresh_selection!(sh)
    dirty(:tree) && _refresh_frames!(sh)
    dirty(:drive) && _refresh_drive!(sh)
    dirty(:time) && _refresh_time!(sh)
    dirty(:optimize) && _refresh_optimize!(sh)
    dirty(:review) && _refresh_review!(sh)
    dirty(:talis) && _refresh_talis!(sh)
    dirty(:sweep) && _refresh_sweep!(sh)
    if dirty(:tree)
        tree = [TreeItem(r...) for r in model_tree(sess)]
        _set!(sh, "selectionCount", length(selected_ids(sess.view)))
        _sync_rows!(sh.rows[:treeModel], sh.models[:treeModel], tree)
    end
    if dirty(:drive)
        # the flux map (or the x-ray image) when chosen, in place of the 3D view or beside it
        side = sess.drive.side[]
        fig = !side && sess.drive.view[] in (:map, :image) ? sh.map.fig : sh.scene.fig
        sidefig = side ? sh.map.fig : sh.blank
        if fig !== sh.figure
            sh.figure = fig
            sh.app["canvas"] = fig
        end
        if sidefig !== sh.side_figure
            sh.side_figure = sidefig
            sh.app["sideCanvas"] = sidefig
        end
        _set!(sh, "sideView", side)
    end
    return
end

function _refresh_view!(sh::BriareusShell)
    v = sh.sess.view
    _set!(sh, "colorBy", String(v.color[]))
    _set!(sh, "beamStyle", String(v.beams[]))
    _set!(sh, "beamColor", String(v.beam_color[]))
    rs = sh.sess.roe
    sel = selected_beams(sh.sess)
    bs = sh.sess.beams.beams[]
    _set!(sh, "roeText", rs.status[] * (roe_stale(sh.sess) ? " (out of date)" : ""))
    _set!(sh, "roeBeamText", isempty(sel) || rs.report[] === nothing ? "" : roe_beam_text(rs, bs[sel[1]].name))
    _set!(sh, "roePorts", roe_port_drawing(rs))
    _set!(sh, "roeOutageText", roe_outage_text(rs))
    _set!(sh, "roeOutageRows", _roe_rows_json(roe_outage_rows(rs)))
    _set!(sh, "beamLength", v.beam_length[] > 0 ? _num(v.beam_length[]) : "")
    _set!(sh, "cutaway", v.cutaway[])
    _set!(sh, "showOff", v.show_off[])
    _set!(sh, "omega1Contour", String(v.omega1_contour[]))
    _set!(sh, "cutWidth", v.cut_width[])
    _set!(sh, "cutAzimuth", v.cut_azimuth[])
    _set!(sh, "cutKind", String(v.cut_kind[]))
    _set!(sh, "cutAxis", String(v.cut_axis[]))
    lo, hi = v.cut_range[]
    _set!(sh, "cutLo", isfinite(lo) ? _num(lo) : "")
    _set!(sh, "cutHi", isfinite(hi) ? _num(hi) : "")
    n = v.cut_normal[]
    # Components below the digits shown print as 0 instead of -0 or 1e-17.
    _set!(sh, "cutNormal", join((@sprintf("%.3g", abs(x) < 5e-4 ? 0.0 : x) for x in n), ", "))
    _set!(sh, "cutOffset", cut_offset(v))
    θ, φ = cut_angles(v)
    _set!(sh, "cutTheta", @sprintf("%.2f", θ))
    _set!(sh, "cutPhi", @sprintf("%.2f", φ))
    _set!(sh, "cutPoint", [_num(round(x; digits = 4) + 0.0) for x in v.cut_point[]])
    line = v.cut_line[]
    _set!(sh, "cutLine", any(!=(0), line))
    lθ, lφ = acosd(clamp(line[3], -1.0, 1.0)), mod(atand(line[2], line[1]), 360)
    _set!(sh, "cutLineText", any(!=(0), line) ? @sprintf("θ = %.1f°, φ = %.1f°", lθ, lφ) : "")
    # the Direction fields start at the line's angles (or 90°, 0°), its origin at TCC or not
    _set!(sh, "cutLineTheta", any(!=(0), line) ? @sprintf("%.2f", lθ) : "90")
    _set!(sh, "cutLinePhi", any(!=(0), line) ? @sprintf("%.2f", lφ) : "0")
    _set!(sh, "cutAtTcc", norm(Briareus.Vec3(v.cut_point[]...)) < 1e-9)
    _set!(sh, "cutRotation", v.cut_rotation[])
    sel = selected_component(sh.sess)
    _set!(sh, "cutSelected", sel === nothing ? "" : component_label(sel))
    _set!(sh, "cutScale", 1.1 * sh.sess.target.target[].scale)
    t = sh.sess.target.target[]
    _set!(sh, "cutChoices", [String(c.name) for c in t.components if c.kind === :surface])
    _set!(sh, "cutLabels", [component_label(c.name) for c in t.components if c.kind === :surface])
    _set!(sh, "cutSet", sort([String(n) for n in v.cut_components[]]))
    busy = String[]
    dr, st, rv, tl = sh.sess.drive, sh.sess.optimize, sh.sess.review, sh.sess.talis
    isempty(sh.sess.meshing[]) || push!(busy, sh.sess.meshing[])
    waiting = Controllers.target_waiting(sh.sess.target)
    isempty(waiting) || push!(busy, waiting)
    dr.building[] && push!(busy, "drive " * Controllers._seconds(time() - dr.started[]))
    st.running[] && push!(busy, "optimizing")
    rv.running[] && push!(busy, "clearance $(rv.progress[][1])/$(rv.progress[][2])")
    sw = sh.sess.sweep
    sw.running[] && push!(busy, "sweep $(sw.progress[][1])/$(sw.progress[][2])")
    tl.running[] && push!(busy, "TaLIS")
    tm = sh.sess.time
    tm.running[] && push!(busy, "drive history $(tm.progress[][1])/$(tm.progress[][2])")
    _set!(sh, "busy", join(busy, " · "))
    return
end

function _refresh_target!(sh::BriareusShell)
    ts = sh.sess.target
    _set!(sh, "targetResolution", String(ts.resolution[]))
    _set!(sh, "mesh_h", _num(ts.mesh[].h))
    _set!(sh, "mesh_capsule_h", _num(ts.mesh[].capsule_h))
    _set!(sh, "targetSummary", target_summary(ts))
    _set!(sh, "targetError", ts.error[])
    # the windows of the selected window's hohlraum (or the selected or first hohlraum part)
    i = _window_part_index(sh.sess)
    ws = i === nothing ? WindowSpec[] : ts.parts[][i].windows
    rows = [WindowRow(k, String(window_name(k)), _num(w.phi), _num(w.z), _num(w.width),
                      _num(w.height), w.open, _num(w.albedo), _num(w.conversion), w.round)
            for (k, w) in enumerate(ws)]
    _sync_rows!(sh.rows[:windowModel], sh.models[:windowModel], rows)
    _set!(sh, "windowCount", length(ws))
    _set!(sh, "hasHohlraum", any(p -> p.kind === :hohlraum, ts.parts[]))
    _set!(sh, "partCount", length(ts.parts[]))
    # where Move into assembly can put the selection: the top level, or an assembly
    _set!(sh, "assemblyChoices", ["Top level"; [String(a.name) for a in assemblies(ts)]])
    names = component_names(ts.target[])
    obs = observed_name(ts)
    _set!(sh, "hasCapsule", obs !== nothing && obs in names)
    # the observed surface in the calculations: a capsule or a flat face (the map, the keep-out)
    _set!(sh, "observedKind", obs !== nothing && obs in names ? (observed_flat(ts) ? "face" : "capsule") : "")
    foil = active_foil(ts)
    _set!(sh, "hasFoil", foil !== nothing && Controllers._has_part(ts.target[], foil))
    # the foils to choose the active one from (shown with two or more)
    fs = foil_names(ts)
    _set!(sh, "foilNames", String.(fs))
    _set!(sh, "activeFoil", something(findfirst(==(foil), fs), 1) - 1)
    return
end

# The hohlraum part whose windows the window page and Add window take: the selected
# window's, else the selected hohlraum part (or a piece of one), else the first.
function _window_part_index(sess)
    ts = sess.target
    kind, key = selection_kind(sess.view.selection[])
    if kind === :window
        w = window_of(ts, Symbol(key))
        w === nothing || return w[1]
    elseif kind in (:part, :comp)
        n = Symbol(key)
        for (i, p) in enumerate(ts.parts[])
            p.kind === :hohlraum && (p.name === n || startswith(key, String(p.name) * "_")) && return i
        end
    end
    return findfirst(p -> p.kind === :hohlraum, ts.parts[])
end

function _refresh_beams!(sh::BriareusShell)
    bt = sh.sess.beams
    _set!(sh, "preset", String(bt.preset[]))
    _set!(sh, "presetLabel", Controllers.beams_label(bt))
    _set!(sh, "beamCount", length(bt.beams[]))
    _set!(sh, "beamsBy", String(bt.by[]))
    _set!(sh, "beamsSummary", beams_summary(bt))
    _set!(sh, "beamsError", bt.error[])
    _set!(sh, "beamsPath", bt.path[])
    _set!(sh, "beamsWarnings", join(bt.warnings[], "\n"))
    _set!(sh, "hasQuads", !isempty(bt.beams[]) && all(b -> !isempty(b.quad), bt.beams[]))
    _set!(sh, "beamNames", [b.name for b in bt.beams[]])
    _set!(sh, "facility", String(bt.facility[]))
    # (shown in the display order; each row keeps its number, which edits refer to)
    br = beam_rows(bt)
    # the power column's unit: per quad when every row is whole quads (a quad is the unit
    # whether it is held as one beam or four), per beam when every row is beams
    units = unique(r.power_unit for r in br)
    _set!(sh, "beamPowerUnit", length(units) == 1 ? String(units[1]) : "mixed")
    col, desc = sh.sess.view.beam_sort[]
    _set!(sh, "beamSort", String(col))
    _set!(sh, "beamSortDesc", desc)
    rows = [(r = br[i]; f = r.aim_frame;
             BeamRow(i, r.key, r.count, @sprintf("%.2f", r.theta), isnan(r.phi) ? "–" : @sprintf("%.2f", r.phi),
                     _fix(r.z), _fix(r.r), _fix(r.phi_offset), _fix(r.defocus), _num(r.power),
                     r.enabled === missing ? 2 : r.enabled ? 1 : 0,
                     f === nothing ? "mixed" : f === :tcc ? "chamber" : String(f),
                     f === nothing ? "The row's beams aim in different frames" :
                     "Aims held in " * (f === :tcc ? "the chamber (TCC)" : _frame_label(f))))
            for i in beam_row_order(br, col, desc)]
    _sync_rows!(sh.rows[:beamModel], sh.models[:beamModel], rows)
    return
end

const _BY_OF_KIND = Dict(:group => :group, :quad => :quad, :beam => :beam)

# The inspector shows what is selected and its editable values.
function _refresh_selection!(sh::BriareusShell)
    sess = sh.sess
    id = sess.view.selection[]
    kind, key = selection_kind(id)
    _set!(sh, "selection", id)
    _set!(sh, "selKind", String(kind))
    _set!(sh, "selKey", key)
    kind in (:comp, :part, :window) && _refresh_component_props!(sh)
    if kind in (:part, :frame, :assembly)
        _refresh_part_props!(sh, key)
    elseif kind === :window
        _refresh_window_props!(sh, key)
    elseif kind in (:group, :quad, :beam)
        _refresh_beam_props!(sh, kind, key)
    end
    return
end

# The selected component, part or window as a component: what it is, its elements and area,
# its display mode, and the role of a piece of a hohlraum part.
function _refresh_component_props!(sh::BriareusShell)
    sess = sh.sess
    t = sess.target.target[]
    name = selected_component(sess)
    # a component left out of the calculations is described from the ghost
    g = sess.target.ghost[]
    tc = any(c -> c.name === name || c.part === name, t.components) || g === nothing ? t : g
    ci = findfirst(c -> c.name === name, tc.components)
    _set!(sh, "compIncluded", !(name in sess.target.excluded[]))
    if ci === nothing && any(c -> c.part === name, tc.components)
        # a part made in pieces (a hohlraum part): its pieces together
        es = elements(tc, name)
        np = count(c -> c.part === name, tc.components)
        _set!(sh, "compInfo", @sprintf("%d pieces · %d elements · %.4g mm²%s", np, length(es),
                                       sum(tc.elem_area[e] for e in es; init = 0.0),
                                       tc === t ? "" : " · not in the calculations"))
        _set!(sh, "compDisplay", findfirst(==(display_mode(sess.view, name)), DISPLAY_MODES) - 1)
        _set!(sh, "selTitle", String(name))
    elseif ci !== nothing
        c = tc.components[ci]
        # an entrance hole is shown with its closing lip (a split for the calculations)
        li = findfirst(l -> is_lip(l) && shown_component(l.name) === name, tc.components)
        es = li === nothing ? elements(tc, name) : [elements(tc, name); elements(tc, tc.components[li].name)]
        A = sum(tc.elem_area[e] for e in es; init = 0.0)
        what = li !== nothing ?
               @sprintf("closing entrance hole (rays pass; its outer lip lets x-rays escape until %.3g ns, then re-emits them)",
                        tc.components[li].material.albedo.closure) :
               c.kind === :aperture ? "aperture (transparent: rays pass, x-rays escape)" :
               c.two_sided ? "two-sided surface" : "one-sided surface"
        mat = c.kind === :aperture ? "" :
              @sprintf(" · albedo %.3g · conversion %.3g", albedo(c.material), c.material.conversion)
        surf = c.surface === nothing ? "faceted" : "exact " * replace(string(nameof(typeof(c.surface))), "Surface" => "") * " surface"
        _set!(sh, "compInfo", @sprintf("%s · %d elements · %.4g mm² · %s%s%s", what, length(es), A,
                                       surf, mat, tc === t ? "" : " · not in the calculations"))
        _set!(sh, "compDisplay", findfirst(==(display_mode(sess.view, name)), DISPLAY_MODES) - 1)
        _set!(sh, "selTitle", String(name))
    end
    # a hohlraum part's piece: its role, and those it can take (indices into PIECE_ROLES)
    pr = name === nothing ? nothing : piece_roles(sess.target, name)
    _set!(sh, "pieceRoleable", pr !== nothing && length(pr.roles) > 1)
    _set!(sh, "pieceRole", pr === nothing ? 0 : findfirst(==(pr.role), PIECE_ROLES) - 1)
    _set!(sh, "pieceRoleAllowed", [pr !== nothing && r in pr.roles for r in PIECE_ROLES])
    return
end

# The selected part, frame or assembly `key`: its kind, name and file, its frame and the
# frames of the move row, its flags and role, its fields, and its placement steps.
function _refresh_part_props!(sh::BriareusShell, key)
    sess = sh.sess
    i = part_index(sess.target, Symbol(key))
    i === nothing && return
    p = sess.target.parts[][i]
    _set!(sh, "partNumber", i)
    _set!(sh, "partKind", is_port_frame(p) ? "Port frame" : PART_KIND_LABELS[p.kind])
    if p.kind === :frame
        _set!(sh, "selTitle", String(p.name))
        _set!(sh, "framePort", p.object)
        _set!(sh, "frameTriad", p.name in sess.view.triads[])
        _set!(sh, "frameInfo", _frame_info(sess, p))
    elseif p.kind === :assembly
        _set!(sh, "selTitle", String(p.name))
        _set!(sh, "assemblyInfo", _assembly_info(sess, p))
    end
    _set!(sh, "partName", String(p.name))
    src = cad_source(sess, p.file)
    _set!(sh, "partFile", src)
    _set!(sh, "partCopy", src != p.file)            # uses the copy saved in the session
    _set!(sh, "partSourceExists", src != p.file && isfile(src))
    # Its frame is that of the target or of another part.
    frames = _part_frame_choices(sess, p)
    _set!(sh, "frameLabels", [_frame_label(f) for f in frames])
    _set!(sh, "partFrame", something(findfirst(==(p.frame), frames), 0) - 1)
    # The move row moves along an axis of any frame, the part frame by default.
    mf = _move_frame_choices(sess)
    _set!(sh, "moveFrames", [_frame_label(f) for f in mf])
    # (the part frame when another part is selected; kept by name as frames come and go)
    get(sh.shown, "_movePart", nothing) === p.name || (sh.shown["_movePart"] = p.name;
                                                        sh.shown["_moveFrameName"] = p.name)
    _set!(sh, "moveFrame", something(findfirst(==(_move_frame(sh)), mf), 1) - 1)
    _set!(sh, "partCoords", findfirst(==(p.coords), PART_COORDINATES) - 1)
    _set!(sh, "partTwoSided", p.two_sided)
    _set!(sh, "partInward", p.inward)
    _set!(sh, "partAperture", p.aperture)
    _set!(sh, "partRole", findfirst(==(part_role(p)), PART_ROLES) - 1)
    _set!(sh, "partFilm", p.name in sess.drive.films[])
    pos = Dict(part_position(p))
    fields = is_port_frame(p) ?
             [FieldRow("rotation", "Roll about the axis", "°", _num(p.values[:rotation]))] :
             [[FieldRow(String(k), label, unit, _num(p.values[k]))
               for (k, label, unit, _) in PART_SIZE_FIELDS[p.kind]];
              [FieldRow(String(k), label, unit, _num(pos[k]))
               for (k, label, unit) in PART_POSITION_FIELDS[p.coords]];
              [FieldRow(String(k), label, k === :h ? part_size_unit(p) : unit, _num(p.values[k]))
               for (k, label, unit, _) in PART_FIELDS
               if !(k in (:x, :y, :z)) && haskey(p.values, k)]]
    # a hohlraum's ends and closing entrance holes
    hoh = p.kind === :hohlraum
    _set!(sh, "partTopEnd", hoh ? Int(p.values[:top_end]) : 1)
    _set!(sh, "partBottomEnd", hoh ? Int(p.values[:bottom_end]) : 1)
    _set!(sh, "partLehClosure", hoh && p.values[:leh_closure] != 0)
    _set!(sh, "partCapsuleMesh", p.kind === :sphere && p.values[:capsule_mesh] != 0)
    hoh && p.values[:bottom_end] == 1 &&
        insert!(fields, length(PART_SIZE_FIELDS[:hohlraum]) + 1,
                FieldRow("bottom_leh_radius", "Bottom LEH radius", "mm (0: the LEH radius)",
                         _num(p.values[:bottom_leh_radius])))
    hoh && p.values[:leh_closure] != 0 &&
        append!(fields, [FieldRow(String(k), label, unit, _num(p.values[k])) for (k, label, unit, _) in LEH_CLOSURE_FIELDS])
    # The profile of a hohlraum gives the waist and end radii of a frustraum, and
    # the breaks of a section list form a separate table (its radius field unused).
    prof = hoh ? Int(p.values[:profile]) : 0
    _set!(sh, "partProfile", prof)
    _set!(sh, "partSectionZ", hoh && prof == 2 ? [_num(z) for (z, _) in p.sections] : String[])
    _set!(sh, "partSectionR", hoh && prof == 2 ? [_num(r) for (_, r) in p.sections] : String[])
    if hoh && prof != 0
        k = findfirst(f -> f.key == "radius", fields)
        if prof == 1
            fields[k] = FieldRow("radius", "Waist radius", "mm", fields[k].value)
            insert!(fields, k + 1, FieldRow("end_radius", "End radius", "mm (top)", _num(p.values[:end_radius])))
            insert!(fields, k + 2, FieldRow("bottom_end_radius", "Bottom end radius", "mm (0: the end radius)",
                                            _num(p.values[:bottom_end_radius])))
        else
            deleteat!(fields, k)
        end
    end
    _sync_rows!(sh.rows[:partFieldModel], sh.models[:partFieldModel], fields)
    rows, note = _step_rows(sess, i)
    _sync_rows!(sh.rows[:stepModel], sh.models[:stepModel], rows)
    _set!(sh, "placementNote", note)
    return
end

# The selected window `key`: its number among the windows of its hohlraum part.
function _refresh_window_props!(sh::BriareusShell, key)
    sess = sh.sess
    w = window_of(sess.target, Symbol(key))
    _set!(sh, "windowNumber", w === nothing ? 0 : w[2])
    return
end

# The selected beams (a group, a quad, a beam, or several of them): their shared values,
# their aim frame, the quads to split or merge, and where the first beam lands.
function _refresh_beam_props!(sh::BriareusShell, kind::Symbol, key)
    sess = sh.sess
    t = sess.target.target[]
    bt = sess.beams
    idx = selected_beams(sess)
    isempty(idx) && return
    n = beam_set_values(bt, idx)
    bs = bt.beams[]
    nsel = length(selected_ids(sess.view))
    gname = split(key, '@')
    _set!(sh, "selTitle", nsel > 1 ? "$(nsel) selections · $(n.count) beams" :
                          kind === :group ? "Group " * key :
                          kind === :quad ? "Quad " * gname[1] * (length(gname) > 1 ? " (" * gname[2] * ")" : "") :
                          "Beam " * key)
    _set!(sh, "bnBy", nsel > 1 ? "set" : String(_BY_OF_KIND[kind]))
    for (f, v) in ((:z, n.z), (:r, n.r), (:phi, n.phi), (:defocus, n.defocus),
                   (:aim_x, n.aim_local[1]), (:aim_y, n.aim_local[2]), (:aim_z, n.aim_local[3]))
        _set!(sh, "bn_" * String(f), _fix(v))
    end
    _set!(sh, "bn_power", _num(n.power))
    _set!(sh, "bnPowerUnit", String(n.power_unit))
    _set!(sh, "bnMirrored", n.mirrored)
    _set!(sh, "bnOn", n.enabled === missing ? 2 : n.enabled ? 1 : 0)
    _set!(sh, "bnMatrix", beam_matrix(sess, idx))
    _set!(sh, "bnCount", n.count)
    cls = unique(bt.cluster[][idx])
    _set!(sh, "bnCluster", length(cls) == 1 ? String(cls[1]) : "")
    _set!(sh, "clusterNames", String.(clusters(bt)))
    # the frame the beams aim in (its axes and origin for the aim point's fields)
    choices = _aim_choices(sess)
    _set!(sh, "aimChoices", [_frame_label(c) for c in choices])
    f = n.aim_frame
    _set!(sh, "bnAimFrame", f === nothing ? -1 : something(findfirst(==(f), choices), 1) - 1)
    _set!(sh, "bnAimFrameLabel", f === nothing ? "the chamber (the beams aim in different frames)" :
                                 f === :tcc ? "the chamber" : _frame_label(f))
    # NIF quads held as one beam can be split, and those held as four merged.
    qs = unique(bs[i].quad for i in idx)
    _set!(sh, "bnQuad", length(qs) == 1 ? qs[1] : "")
    _set!(sh, "bnSplit", length(Controllers.quads_to_split(bt, idx)))
    _set!(sh, "bnMerge", length(Controllers.quads_to_merge(bt, idx)))
    b = bs[idx[1]]
    p = beam_path(t, b)
    land = n.count > 1 ? "" : p.component === :missed ? " · misses the target" :
           " · axis lands on " * _in_text(p.component)
    spot = @sprintf("spot %.3g × %.3g mm (order %.3g), f/%.3g", b.spot.a, b.spot.b,
                    b.spot.order, b.fnumber)
    θ, φ = port_angles(b)
    angles = n.count == 1 ? @sprintf("port θ %.2f°, φ %.2f°", θ, φ) : ""
    _set!(sh, "bnInfo", join(filter(!isempty, [@sprintf("%d beam%s", n.count, n.count == 1 ? "" : "s"),
                                                angles, spot]), " · ") * land)
    return
end

function _refresh_drive!(sh::BriareusShell)
    sess = sh.sess
    dr = sess.drive
    res = dr.result[]
    _set!(sh, "driveBuilding", dr.building[])
    _set!(sh, "driveElapsed", dr.building[] ? Controllers._seconds(time() - dr.started[]) : "")
    _set!(sh, "driveEvaluating", dr.evaluating[])
    _set!(sh, "driveBuilt", dr.scenario[] !== nothing)
    _set!(sh, "driveStale", drive_stale(sess))
    _set!(sh, "driveCurrent", Controllers._drive_current(sess))
    _set!(sh, "driveStatus", dr.status[])
    _set!(sh, "driveHasResult", res !== nothing)
    _set!(sh, "driveField", String(dr.field[]))
    _set!(sh, "driveView", String(dr.view[]))
    _set!(sh, "spotRadius", _num(dr.spot_radius[]))
    r1, r2 = dr.unconverted_ratio[]
    _set!(sh, "ucRatio1", r1 > 0 ? _num(r1) : "")
    _set!(sh, "ucRatio2", r2 > 0 ? _num(r2) : "")
    _set!(sh, "ucLive", [r1 > 0, r2 > 0])
    p1, p2, a, b = sess.time.prepulse[]
    _set!(sh, "ppOmega1", p1 > 0 ? _num(p1) : "")
    _set!(sh, "ppOmega2", p2 > 0 ? _num(p2) : "")
    _set!(sh, "ppStart", _num(a))
    _set!(sh, "ppStop", _num(b))
    _set!(sh, "irrThreshold", _num(sess.time.threshold[]))
    _set!(sh, "unconvertedText", res === nothing ? "" : unconverted_summary(res))
    _set!(sh, "driveSummary", res === nothing ? "" : drive_summary(res))
    _set!(sh, "driveMean", res === nothing || isnan(res.mean_flux) ? "" :
                           @sprintf("%.3g W/cm² · Tr %.1f eV", to_Wcm2(res.mean_flux),
                                    radiation_temperature(res.mean_flux)))
    _set!(sh, "driveFace", res !== nothing && res.face !== nothing)
    _set!(sh, "visradRows", _visrad_rows_json(visrad_parity_rows(sess)))
    # each cavity's wall temperature, the one shown chosen in the Drive panel
    cavs = res === nothing ? NamedTuple[] : res.cavities
    shown = drive_cavity(sess)
    k = something(findfirst(r -> r.name === shown, cavs), 1)
    _set!(sh, "cavityNames", String[String(r.name) for r in cavs])
    _set!(sh, "cavityIndex", k - 1)
    _set!(sh, "cavityTr", isempty(cavs) ? "" :
                          @sprintf("%.3g W/cm² · Tr %.1f eV", to_Wcm2(cavs[k].flux), cavs[k].Tr))
    s = res === nothing ? nothing : res.spot
    _set!(sh, "spotText", s === nothing ? "" :
                          @sprintf("in r ≤ %.3g mm: mean %.2e W/cm², peak %.2e W/cm², rms nonuniformity %.1f%%, %.1f%% of the delivered power",
                                   s.radius, s.mean, s.peak, 100 * s.rms, 100 * s.fraction))
    modes = res === nothing ? ModeRow[] :
            [ModeRow(m.l, _pct(m.epsilon), _pct(m.legendre), "", "") for m in res.modes]
    _sync_rows!(sh.rows[:modeModel], sh.models[:modeModel], modes)
    power = res === nothing ? PowerRow[] :
            [PowerRow(p.name, _power_text(p.power), _power_pct(p.fraction))
             for p in res.power]
    _sync_rows!(sh.rows[:powerModel], sh.models[:powerModel], power)
    _set!(sh, "delivered", res === nothing ? "" :
                           @sprintf("%.4g TW delivered", res.evaluation.deposition.delivered))
    _set!(sh, "driveCavity", dr.cavity[])
    _set!(sh, "beamModel", String(dr.beam_model[]))
    _set!(sh, "cavityPoint", join((_num(x) for x in dr.cavity_point[]), ", "))
    # the cavity's account changes with the drive's build only
    if get(sh.shown, "_enclosure", missing) !== dr.enclosure[]
        sh.shown["_enclosure"] = dr.enclosure[]
        _set!(sh, "cavitySummary", cavity_summary(dr))
    end
    _set!(sh, "retuneRunning", dr.retune_running[])
    _set!(sh, "retuneStatus", dr.retune_status[])
    # what a retune can hold (the dialog's choice when there are several, and the level's units)
    qs = Controllers._retune_quantities(sess)
    _set!(sh, "retuneOn", String[String(q.on) for q in qs])
    _set!(sh, "retuneKinds", String[String(q.kind) for q in qs])
    _set!(sh, "retuneLabels", String[q.label for q in qs])
    return
end

# A diagnostic's tooltip: its definition and where it comes from, and what keeps it from use.
function _diag_tip(r)
    origin = r.origin === :builtin ? "built-in" : r.origin === :user ? "from a diagnostic directory" :
             "saved with the session"
    lines = [r.id * " · " * origin * (r.public ? "" : " · not public")]
    r.usable || push!(lines, "Cannot be built: " * r.problem)
    isempty(r.changed) || push!(lines, "Differs from the directories' definition in: " * join(r.changed, ", "))
    return join(lines, "\n")
end

# The Diagnostics tab's rows: one per instrument, its recorded configuration (or the one
# a tick would record) and its configurations.
function _diag_rows(sess)
    rows = DiagRow[]
    for r in diagnostic_instruments(sess)
        cs = r.configs
        ids, labels = join((x.id for x in cs), "|"), join((x.config for x in cs), "|")
        ks = isempty(r.recorded) ? [something(findfirst(c -> c.usable, cs), 1)] :
             [findfirst(c -> c.id == id, cs) for id in r.recorded]
        for (j, k) in enumerate(ks)
            c = cs[k]
            # (the configuration is named when there is a choice, or when it says something)
            push!(rows, DiagRow(c.id, r.instrument, String(r.kind), !isempty(r.recorded), r.usable,
                                !isempty(c.changed), _diag_tip(c),
                                length(cs) > 1 || c.config != "default" ? c.config : "", ids, labels, !r.here, j == 1))
        end
    end
    return rows
end

# The image the Diagnostics tab shows: the chosen one of the history's (`_diagImage`, an
# `(id, frame)`), else its first; drawn again when the history, the choice or (for the
# image at the cursor) the cursor's step changes.
function _refresh_diag_image!(sh::BriareusShell, images)
    sess = sh.sess
    sel = get(sh.shown, "_diagImage", nothing)
    k = sel === nothing ? nothing : findfirst(r -> (r.id, r.frame) == sel, images)
    k === nothing && !isempty(images) && (k = 1)
    _set!(sh, "diagImageIndex", something(k, 0) - 1)
    pick = k === nothing ? nothing : (images[k].id, images[k].frame)
    key = (sess.time.history[], pick, pick !== nothing && pick[2] == 0 ? time_step(sess.time) : nothing)
    isequal(get(sh.shown, "_diagImageKey", missing), key) && return
    sh.shown["_diagImageKey"] = key
    im = pick === nothing ? nothing : Controllers.diagnostic_image(sess, pick...)
    show_image!(sh.diag_image, im; empty = "No image yet")
    return
end

function _refresh_time!(sh::BriareusShell)
    sess = sh.sess
    tm = sess.time
    # The diagnostics chart belongs to the Diagnostics tab, so a session left on it shows
    # the drive.
    tm.chart[] === :diagnostics && set_chart!(tm, :drive)
    ch = time_chart(sess)
    h = tm.history[]
    _set!(sh, "timeChart", String(ch.chart))
    _set!(sh, "timeTitle", ch.title)
    _set!(sh, "timeEmpty", ch.empty)
    _set!(sh, "timeRunning", tm.running[])
    k, n = tm.progress[]
    _set!(sh, "timeProgress", tm.running[] ? @sprintf("step %d of %d · %s", k, n,
                                                       Controllers._seconds(time() - tm.started[])) : "")
    _set!(sh, "timeFraction", n > 0 ? k / n : 0.0)
    _set!(sh, "timeStatus", tm.status[])
    _set!(sh, "timeHasPulses", tm.pulses[] !== nothing)
    _set!(sh, "timeChecking", tm.checking[])
    checks = irradiance_checks(sess)
    _set!(sh, "irrTexts", String[c.text for c in checks])
    _set!(sh, "irrPass", Bool[c.pass for c in checks])
    _set!(sh, "irrTips", String[c.tip for c in checks])
    _set!(sh, "irrStale", irradiance_stale(sess))
    _set!(sh, "timeHasHistory", h !== nothing)
    _set!(sh, "timeStale", history_stale(sess))
    _set!(sh, "timeSteps", h === nothing ? 0 : length(h.time))
    step = time_step(tm)
    _set!(sh, "timeStep", step === nothing ? 0 : step)
    _set!(sh, "timeLabel", time_label(sess))
    # The 3D view shows the drive at the beam powers or a step of the history.
    _set!(sh, "timeMode", time_label(sess) == "" ? "beams" : "history")
    step === nothing || (sh.shown["_lastStep"] = step)
    src = tm.pulses[]
    _set!(sh, "pulsesLabel", src === nothing ? "" : src.label)
    # What the cone or quad columns of a pulse table are given per, or "" where there is
    # nothing to choose: a table per beam, a shape, or beams without quads.
    _set!(sh, "pulseUnit", src !== nothing && src.kind === :table && src.by !== :beam &&
                           any(b -> !isempty(b.quad), sess.beams.beams[]) ? String(src.unit) : "")
    images = diagnostic_images(sess)
    _sync_rows!(sh.rows[:diagImageModel], sh.models[:diagImageModel],
                ImageRow[ImageRow(r.id, r.label, r.frame) for r in images])
    _refresh_diag_image!(sh, images)
    # the tab's canvas: the image on the imagers' category, the chart otherwise (a figure
    # is shown on one canvas at a time)
    dfig = get(sh.shown, "diagCategory", "") == "xray_imager" ? sh.diag_image.fig : sh.diag.fig
    if dfig !== sh.shown["_diagFigure"]
        sh.shown["_diagFigure"] = dfig
        sh.app["diagCanvas"] = dfig
        sh.diag.dirty[] = true
    end
    rows = diagnostic_catalog(sess)
    drows = _diag_rows(sess)
    _sync_rows!(sh.rows[:diagModel], sh.models[:diagModel], drows)
    # the instruments recorded, per category (the tab's toggle shows the counts)
    _set!(sh, "diagCounts", [count(r -> r.ticked && r.kind == String(k), drows) for (k, _) in Controllers._DIAGNOSTIC_GROUPS])
    _set!(sh, "diagKinds", [String(k) for (k, _) in Controllers._DIAGNOSTIC_GROUPS])
    _set!(sh, "diagKindLabels", [l for (_, l) in Controllers._DIAGNOSTIC_GROUPS])
    dch = time_chart(sess, :diagnostics)
    _set!(sh, "diagEmpty", dch.empty)
    _set!(sh, "diagTitle", dch.title)
    _sync_rows!(sh.rows[:diagLegendModel], sh.models[:diagLegendModel],
                LegendRow[LegendRow(r.label, r.color, String(r.style)) for r in dch.legend])
    if get(sh.shown, "_catalog", nothing) !== rows             # the folders and problems, once per reading
        sh.shown["_catalog"] = rows
        _sync_rows!(sh.rows[:folderModel], sh.models[:folderModel],
                    FolderRow[FolderRow(f.path, f.saved) for f in diagnostic_folders()])
        _sync_rows!(sh.rows[:issueModel], sh.models[:issueModel],
                    IssueRow[IssueRow(String(i.level), basename(i.file) * (isempty(i.key) ? "" : " [" * i.key * "]"),
                                      i.message, i.file) for i in tm.issues[]])
        _set!(sh, "diagCount", length(rows))
        _set!(sh, "diagFolders", length(sh.rows[:folderModel]))
        _set!(sh, "diagIssues", count(i -> i.level !== :warning, tm.issues[]))
    end
    _sync_rows!(sh.rows[:timeLegendModel], sh.models[:timeLegendModel],
                LegendRow[LegendRow(r.label, r.color, String(r.style)) for r in ch.legend])
    return
end

function _refresh_optimize!(sh::BriareusShell)
    st = sh.sess.optimize
    _set!(sh, "optBy", String(st.by[]))
    for f in OPTIMIZE_FIELDS
        _set!(sh, "optFree_" * String(f), f in st.fields[])
        _set!(sh, "optDelta_" * String(f), _num(st.deltas[][f]))
    end
    _set!(sh, "optObjective", String(st.objective[]))
    for l in 1:4
        _set!(sh, "optW$l", _num(get(st.weights[], l, 0.0)))
    end
    _set!(sh, "optHold", st.hold_drive[])
    _set!(sh, "optKeepout", st.keepout[])
    _set!(sh, "optKeepoutMin", _num(st.keepout_min[]))
    _set!(sh, "optCapsule", st.capsule_keepout[])
    _set!(sh, "optCapsuleMin", _num(st.capsule_min[]))
    _set!(sh, "optSpill", _num(st.spill[]))
    _set!(sh, "optUniformity", _num(st.uniformity[]))
    _set!(sh, "optIterations", st.iterations[])
    _set!(sh, "optMaxEvaluations", st.max_evaluations[])
    _set!(sh, "optRunning", st.running[])
    _set!(sh, "optProgress", optimize_progress(st))
    _set!(sh, "optStatus", st.status[])
    _set!(sh, "optError", st.error[])
    r = st.result[]
    _set!(sh, "optHasResult", r !== nothing)
    _set!(sh, "optApplied", r !== nothing && r.applied)
    _set!(sh, "optStale", optimize_stale(sh.sess))
    _set!(sh, "optValues", r === nothing ? "" : st.status[])
    _set!(sh, "optClipped", r === nothing ? "" :
                            "clipped laser power " * _pct(r.before.clipped) * " → " * _pct(r.after.clipped) *
                            (isfinite(r.after.capsule_clearance) ?
                             " · on the capsule " * _pct(r.before.capsule, 3) * " → " * _pct(r.after.capsule, 3) *
                             @sprintf(" · capsule worst-case clearance %.3f → %.3f mm",
                                      r.before.capsule_clearance, r.after.capsule_clearance) : ""))
    _set!(sh, "optSpot", r === nothing ? "" :
                         r.after.spot !== nothing ?
                         @sprintf("spot rms nonuniformity %.1f%% → %.1f%%", 100 * r.before.spot.rms,
                                  100 * r.after.spot.rms) :
                         get(r.after, :face, nothing) !== nothing ?
                         @sprintf("face rms nonuniformity %.2f%% → %.2f%% · mean flux %.3g → %.3g W/cm²",
                                  100 * r.before.face.rms, 100 * r.after.face.rms,
                                  to_Wcm2(r.before.face.mean), to_Wcm2(r.after.face.mean)) : "")
    modes = r === nothing || isempty(r.after.modes) ? ModeRow[] :
            [ModeRow(a.l, _pct(a.epsilon), _pct(a.legendre), _pct(b.epsilon), _pct(a.epsilon))
             for (b, a) in zip(r.before.modes, r.after.modes)]
    _sync_rows!(sh.rows[:optModeModel], sh.models[:optModeModel], modes)
    params = r === nothing ? ParamRow[] :
             [ParamRow(p.name, _fix(p.start), _fix(p.optimized), _fix(p.lower), _fix(p.upper), p.at_bound)
              for p in parameter_rows(r)]
    _sync_rows!(sh.rows[:paramModel], sh.models[:paramModel], params)
    return
end

# A component named within a sentence: "the top cap"-style labels in lower case, a
# user's part name as it is.
_in_text(n::Symbol) = (l = component_label(n); l == String(n) ? l : lowercasefirst(l))

function _clear_row(x)
    return ClearRow(x.beam, component_label(x.component), String(x.status), _fix(x.clearance),
                    _fix(x.nominal), @sprintf("%.3g", x.power), @sprintf("%.3g", x.peak_intensity))
end

function _refresh_review!(sh::BriareusShell)
    sess = sh.sess
    rv = sess.review
    e = rv.errors[]
    _set!(sh, "rvPeak", @sprintf("%g", rv.peak_fraction[]))
    _set!(sh, "rvWings", rv.wings[])
    for k in (:pointing, :target_offset, :hardware_tolerance, :focus)
        _set!(sh, "rv_" * String(k), _num(e[k]))
    end
    _set!(sh, "rvQuality", String(rv.quality[]))
    _set!(sh, "rvRunning", rv.running[])
    _set!(sh, "rvDone", rv.progress[][1])
    _set!(sh, "rvTotal", rv.progress[][2])
    _set!(sh, "rvStatus", rv.status[])
    _set!(sh, "rvError", rv.error[])
    _set!(sh, "rvHasReport", rv.report[] !== nothing)
    _set!(sh, "rvBrief", review_brief(rv))
    _set!(sh, "rvStale", review_stale(rv, sess.target.target[], active_beams(sess.beams)))
    cm = sess.corridor
    _set!(sh, "cmShown", cm.shown[])
    _set!(sh, "cmRadius", _num(cm.radius[]))
    _set!(sh, "cmText", cm.status[])
    _set!(sh, "cmSurface", String(cm.surface[]))
    _set!(sh, "cmPoint", join((_num(v) for v in cm.point[]), ", "))
    _set!(sh, "cmNormal", join((_num(v) for v in cm.normal[]), ", "))
    _set!(sh, "cmExtent", _num(cm.extent[]))
    _set!(sh, "cmError", cm.error[])
    _set!(sh, "cmMap", cm.map[] === nothing ? "" : corridor_drawing(cm.map[]; tight = sess.talis.tight[]))
    _set!(sh, "cmStale", cm.map[] !== nothing &&
                         corridor_stale(cm, rv, sess.talis, active_beams(sess.beams), corridor_sight_lines(sess);
                                        target = sess.target.target[]))
    _set!(sh, "cmShadows", cm.shadows[])
    comps = [CompRow(String(n), component_label(n), n in rv.targets[] ? "target" : n in rv.hardware[] ? "hardware" : "none")
             for n in surface_names(sess.target.target[])]
    _sync_rows!(sh.rows[:compModel], sh.models[:compModel], comps)
    _sync_rows!(sh.rows[:clearModel], sh.models[:clearModel], ClearRow[_clear_row(x) for x in clearance_rows(rv)])
    _sync_rows!(sh.rows[:worstModel], sh.models[:worstModel], ClearRow[_clear_row(x) for x in worst_rows(rv)])
    return
end

function _refresh_sweep!(sh::BriareusShell)
    sess = sh.sess
    sw = sess.sweep
    subs = sweep_subjects(sess)
    _set!(sh, "swSubjects", String.(subs))
    _set!(sh, "swSubject", something(findfirst(==(sweep_subject(sess)), subs), 0) - 1)
    mf = _move_frame_choices(sess)
    _set!(sh, "swFrames", [_frame_label(f) for f in mf])
    _set!(sh, "swFrame", something(findfirst(==(sweep_frame(sess)), mf), 1) - 1)
    _set!(sh, "swAxis", something(findfirst(==(sw.axis[]), SWEEP_AXES), 3) - 1)
    # a placement step's number: the subject part's steps, and the numbers of the chosen one
    _set!(sh, "swMode", sw.mode[] === :step ? 1 : 0)
    i = part_index(sess.target, sweep_subject(sess))
    psteps = i === nothing ? PlacementStep[] : sess.target.parts[][i].steps
    _set!(sh, "swStepChoices", ["$k · " * PART_STEP_LABELS[step_kind(st)] for (k, st) in enumerate(psteps)])
    k = sw.step_index[]
    _set!(sh, "swStepIndex", 1 <= k <= length(psteps) ? k - 1 : -1)
    nums = 1 <= k <= length(psteps) ? step_numbers(psteps[k]) : Tuple{Symbol,String,String}[]
    _set!(sh, "swNumbers", [l for (_, l, _) in nums])
    j = findfirst(n -> n[1] === sw.step_key[], nums)
    _set!(sh, "swNumber", something(j, 0) - 1)
    _set!(sh, "swUnit", j === nothing ? "mm" : nums[j][3])
    for k in (:start, :stop, :step, :los_theta, :los_phi)
        _set!(sh, "sw_" * String(k), _num(getfield(sw, k)[]))
    end
    for k in (:clearance, :los, :capture, :optimize)
        _set!(sh, "sw_" * String(k), getfield(sw, k)[])
    end
    objs = [:none; surface_names(sess.target.target[])]
    _set!(sh, "swObjects", [o === :none ? "(none)" : component_label(o) for o in objs])
    _set!(sh, "swObject", something(findfirst(==(sweep_los_object(sess)), objs), 1) - 1)
    _set!(sh, "swQuality", String(sw.quality[]))
    _set!(sh, "swRunning", sw.running[])
    _set!(sh, "swDone", sw.progress[][1])
    _set!(sh, "swTotal", sw.progress[][2])
    _set!(sh, "swStatus", sw.status[])
    _set!(sh, "swError", sw.error[])
    _set!(sh, "swHasResult", sw.result[] !== nothing)
    _set!(sh, "swBrief", sweep_brief(sw))
    _set!(sh, "swStale", sweep_stale(sess))
    _set!(sh, "swSelected", sw.selected[])
    steps = try
        length(sweep_values(sw))
    catch
        0
    end
    _set!(sh, "swSteps", steps)
    if get(sh.shown, "_sweepResult", missing) !== sw.result[]
        sh.shown["_sweepResult"] = sw.result[]
        rows = [SweepRowItem(k, @sprintf("%+.4g", w.value), _fix(w.clearance), w.worst, w.fails, w.reviews, w.los,
                             _pct(w.blocked, 1), _pct(w.covered, 1), _pct(w.captured, 2),
                             isnan(w.objective_after) ? "–" : @sprintf("%.3g → %.3g", w.objective_before, w.objective_after))
                for (k, w) in enumerate(sweep_rows(sw))]
        _sync_rows!(sh.rows[:sweepModel], sh.models[:sweepModel], rows)
    end
    ch = sweep_chart(sess)
    _set!(sh, "swEmpty", ch.empty)
    _set!(sh, "swTitle", ch.title)
    _sync_rows!(sh.rows[:sweepLegendModel], sh.models[:sweepLegendModel],
                LegendRow[LegendRow(r.label, r.color, String(r.style)) for r in ch.legend])
    return
end

_unit(u::Symbol) = u === :deg ? "°" : String(u)
# A margin or distance; ±Inf means the rule found nothing to measure (no free edge, all
# hardware downstream of the landing point).
_val(x::Real, u::Symbol) = isfinite(x) ? @sprintf("%+.3f %s", x, _unit(u)) : x > 0 ? "n/a" : "–"
_val0(x::Real, u::Symbol) = isfinite(x) ? @sprintf("%.3f %s", x, _unit(u)) : x > 0 ? "n/a" : "–"

function _talis_row(x::TaLISRow)
    return TalisRow(x.beam, get(TALIS_RULE_LABELS, x.rule, String(x.rule)), component_label(x.component),
                    String(x.status), _val(x.margin, x.unit), _val0(x.worst, x.unit),
                    _val0(x.nominal, x.unit), _val0(x.required, x.unit), x.note)
end

function _refresh_talis!(sh::BriareusShell)
    sess = sh.sess
    tl = sess.talis
    _set!(sh, "tlTolerance", _num(tl.tolerance[]))
    _set!(sh, "tlPivot", String(tl.pivot[]))
    _set!(sh, "tlPivotPoint", join((_num(v) for v in tl.pivot_point[]), ", "))
    _set!(sh, "tlOmega1", tl.omega1[])
    _set!(sh, "tlTight", _num(tl.tight[]))
    _set!(sh, "tlHideCurved", tl.hide_curved[])
    _set!(sh, "tlFocus", tl.focus[] === nothing ? -1 : tl.focus[].index - 1)
    _set!(sh, "tlFocusApproach", tl.focus[] !== nothing && tl.focus[].approach !== nothing)
    _set!(sh, "tlFocusText", talis_focus_text(tl))
    _set!(sh, "tlPinch", tl.pinch[])
    sl = tl.focus[] === nothing ? nothing : tl.focus[].slice
    _set!(sh, "tlSlice", sl === nothing ? "" : slice_drawing(sl))
    _set!(sh, "tlSliceText", sl === nothing ? "" : slice_text(sl))
    _set!(sh, "tlRunning", tl.running[])
    _set!(sh, "tlElapsed", tl.running[] ? Controllers._seconds(time() - tl.started[]) : "")
    _set!(sh, "tlStatus", tl.status[])
    _set!(sh, "tlError", tl.error[])
    _set!(sh, "tlHasReport", tl.report[] !== nothing)
    _set!(sh, "tlBrief", talis_brief(tl))
    _set!(sh, "tlStale", talis_stale(tl, sess.review, sess.target.target[], active_beams(sess.beams)))
    comps = [CompRow(String(n), component_label(n), String(talis_role(tl, n)))
             for n in surface_names(sess.target.target[])]
    _sync_rows!(sh.rows[:talisCompModel], sh.models[:talisCompModel], comps)
    # The rows, often hundreds, change only with the report, so they are not re-sent on
    # every refresh.
    if get(sh.shown, "_talisReport", missing) !== tl.report[]
        sh.shown["_talisReport"] = tl.report[]
        rules = [RuleRow(r.label, r.requirement, r.pass, r.review, r.fail,
                         _val(r.worst.margin, r.worst.unit), r.worst.beam * " / " * component_label(r.worst.component),
                         r.fail > 0 ? "fail" : r.review > 0 ? "review" : "pass")
                 for r in talis_rule_rows(tl)]
        _sync_rows!(sh.rows[:ruleModel], sh.models[:ruleModel], rules)
        _sync_rows!(sh.rows[:talisModel], sh.models[:talisModel],
                    TalisRow[_talis_row(x) for x in talis_rows(tl)])
    end
    _sync_rows!(sh.rows[:noteModel], sh.models[:noteModel],
                NoteRow[NoteRow(String(n.kind), n.text) for n in talis_notes(tl)])
    return
end

# The rows of the hypothetical outage of the Outages tab, as the JSON its table reads.
function _roe_rows_json(rows)
    item(r) = @sprintf("{\"criterion\":\"%s\",\"score\":\"%.3g%s\",\"green\":\"%.3g\",\"amber\":\"%.3g\",\"tier\":\"%s\",\"binding\":%s}",
                       escape_string(r.criterion), r.score, isempty(r.label) ? "" : " " * escape_string(r.label),
                       r.green, r.amber, r.tier, r.binding)
    return "[" * join((item(r) for r in rows), ",") * "]"
end

# The summary of the comparison with VisRad in the Drive tab, as the JSON its table reads.
function _visrad_rows_json(rows)
    isempty(rows) && return ""
    pct(x) = isfinite(x) ? @sprintf("%.2f%%", 100x) : "–"
    item(o) = @sprintf("{\"object\":\"%s\",\"briareus\":\"%.4g\",\"visrad\":\"%.4g\",\"relative\":\"%s\",\"rms\":\"%s\",\"max\":\"%s\",\"correlation\":\"%s\"}",
                       escape_string(o.object), o.power_briareus, o.power_visrad,
                       o.power_briareus == 0 ? "–" : pct((o.power_visrad - o.power_briareus) / abs(o.power_briareus)),
                       pct(o.rms), pct(o.max), isfinite(o.correlation) ? @sprintf("%.4f", o.correlation) : "–")
    return "[" * join((item(o) for o in rows), ",") * "]"
end
