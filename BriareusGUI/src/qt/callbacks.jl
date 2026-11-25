# This file defines the functions QML calls: the tick of the window and the `br_*` callbacks,
# with their registration.

# ---------------------------------------------------------------- QML callbacks

const _SHELL = Ref{Union{Nothing,BriareusShell}}(nothing)
# `_TICK_HOOK[](shell)` runs on every tick (scripted smoke tests drive the window
# through it); `request_close()` closes the window from Julia.
const _TICK_HOOK = Ref{Any}(nothing)
const _ON_SHOWN = Ref{Any}(nothing)            # `briareus(; shown)`, until the first tick
const _CLOSE_REQUESTED = Ref(false)
const _LAST_REFRESH = Ref(0.0)

"""
    request_close()

Close the open Briareus window (from a callback or a tick hook).
"""
request_close() = (_CLOSE_REQUESTED[] = true; nothing)

"""
    request_grab(path)

Save an image of the open window's body (panels and viewer) to `path` (PNG) on the
window's next tick, through QML `grabToImage`. The file appears asynchronously.
"""
request_grab(path::AbstractString) = (sh = _SHELL[]; sh === nothing || _set!(sh, "grabPath", String(path)); nothing)
br_grab_started() = (sh = _SHELL[]; sh === nothing || _set!(sh, "grabPath", ""); nothing)

"""
    request_tab(key)

Show a results tab of the open window (`"beams"`, `"drive"`, `"time"`, `"optimize"`,
`"clearance"`, `"sweep"`, `"talis"`), as its menus do.
"""
request_tab(key::AbstractString) = (sh = _SHELL[]; sh === nothing ||
                                    (_set!(sh, "resultsRequest", ""); _set!(sh, "resultsRequest", String(key))); nothing)

"""
    request_window_size(width, height)

Resize the open window to `width × height` (in Qt's device-independent pixels), for
scripted screenshots of the layout at several window sizes.
"""
request_window_size(width::Integer, height::Integer) =
    (sh = _SHELL[]; sh === nothing || _set!(sh, "sizeRequest", string(width, "x", height)); nothing)

"""
    br_help(topic)

Open the window's Help screen at `topic` (a key of `qml/common/help.js`, e.g.
`"talis.tolerance"`; `""` for the first page), as an ⓘ does when clicked.
"""
br_help(topic) = (sh = _SHELL[]; sh === nothing ||
                  (_set!(sh, "helpRequest", ""); _set!(sh, "helpRequest", "topic:" * String(topic))); nothing)

"The documentation's Desktop application page online (Help ▸ Documentation)."
const DOCS_URL = "https://stillyslalom.github.io/Briareus.jl/dev/gui/"

"""
    br_docs_url() -> String

The Desktop application page of the documentation: a local build of a development
checkout (`docs/build`, from `docs/make.jl`) when there is one, otherwise `DOCS_URL`.
"""
function br_docs_url()
    build = normpath(joinpath(@__DIR__, "..", "..", "..", "docs", "build"))
    for page in (joinpath(build, "gui.html"), joinpath(build, "gui", "index.html"))
        isfile(page) && return "file:///" * replace(page, '\\' => '/')
    end
    return DOCS_URL
end

"""
    request_popup(key, path = "")

Open a dialog or menu of the open window: `"retune"`, `"port"` or `"group"` (the
dialogs), `"diagnostics"` (the Diagnostics tab's list of instruments) and `"images"` (the
Diagnostics tab on its imagers), `"topics"` or `"topics:<topic>"` (the Help window, see
`br_help`; `"help"` opens it at the mouse and keys), `"rowmenu:<node id>"` (the model
tree's row menu for a beam row, a frame or an analysis), a menu bar menu by its title
without the `&` (`"File"`, `"View"`, …), or `"window"` for none.
With a `path`, save an image of the popup there (for `"window"`, of the toolbar; PNG,
written asynchronously, as by [`request_grab`](@ref)); `request_popup("")` closes it.
"""
function request_popup(key::AbstractString, path::AbstractString = "")
    sh = _SHELL[]
    sh === nothing && return nothing
    _set!(sh, "popupGrabPath", String(path))
    _set!(sh, "popupRequest", "")
    _set!(sh, "popupRequest", String(key))
    return nothing
end

"""
    request_cut_settings(path = "")

Open the cutaway's Cut settings popup of the open window, as the Cutaway button's
context menu does, and with a `path` save an image of the popup there (PNG, written
asynchronously, as by [`request_grab`](@ref)).
"""
function request_cut_settings(path::AbstractString = "")
    sh = _SHELL[]
    sh === nothing && return nothing
    _set!(sh, "cutGrabPath", String(path))
    _set!(sh, "cutSettingsRequest", false)
    _set!(sh, "cutSettingsRequest", true)
    return nothing
end

# Run a callback body on the current shell, reporting errors in the status line
# instead of letting them escape into Qt.
function _with_shell(f, dirty::UInt16 = _DIRTY_ALL)
    sh = _SHELL[]
    sh === nothing && return nothing
    try
        f(sh)
    catch err
        sh.sess.status[] = "error: " * Controllers._errmsg(err)
        @error "BriareusGUI callback failed" exception = (err, catch_backtrace())
    end
    _mark!(sh.dirty, dirty)
    return nothing
end

# The mask of a callback that only sets observables of the session: each observable marks
# the domains that show it (`_dirty_sources`), and the callback marks the edit.
const _DIRTY_OBSERVED = _dirty_mask()

# Returns 2 (close the window), or the canvases to redraw as bits: 1 the 3D view or the
# map, 4 the time chart, 8 the sweep chart, 16 the Diagnostics tab's chart and image.
function br_tick()
    sh = _SHELL[]
    sh === nothing && return 0
    # Only thread 1 runs libuv's event loop, and it sits in Qt's exec(): pump it here.
    Base.process_events()
    if _ICON_TRIES[] > 0
        _ICON_TRIES[] = _set_window_icon() ? 0 : _ICON_TRIES[] - 1
    end
    if _ON_SHOWN[] !== nothing
        f, _ON_SHOWN[] = _ON_SHOWN[], nothing
        try
            f()
        catch err
            @error "shown callback failed" exception = (err, catch_backtrace())
        end
    end
    if _TICK_HOOK[] !== nothing
        try
            _TICK_HOOK[](sh)
        catch err
            @error "tick hook failed" exception = (err, catch_backtrace())
            request_close()
        end
    end
    _CLOSE_REQUESTED[] && (_CLOSE_REQUESTED[] = false; return 2)
    while isready(sh.queue)
        f = take!(sh.queue)
        try
            f()
        catch err
            @error "BriareusGUI background update failed" exception = (err, catch_backtrace())
        end
    end
    sess = sh.sess
    # After any action, note the design for undo and then show the new state.
    mask = _take_dirty!(sh.dirty)
    if mask != _DIRTY_NONE
        record!(sess)
        _refresh!(sh, mask)
        _LAST_REFRESH[] = time()
    elseif (sess.drive.building[] || sess.optimize.running[] || sess.review.running[] ||
            sess.talis.running[] || sess.time.running[] || sess.sweep.running[] ||
            !target_ready(sess.target)) &&
           time() - _LAST_REFRESH[] > 0.25
        # Running jobs report elapsed time four times a second. This is not a design edit,
        # since a snapshot for undo every tick took a third of the GUI thread during a drive.
        _refresh!(sh)
        _LAST_REFRESH[] = time()
    end
    redraw = sh.scene.dirty[] || sh.map.dirty[]
    chart = sh.time.dirty[]
    swchart = sh.sweep.dirty[]
    dchart = sh.diag.dirty[] || sh.diag_image.dirty[]
    sh.scene.dirty[] = false
    sh.map.dirty[] = false
    sh.time.dirty[] = false
    sh.sweep.dirty[] = false
    sh.diag.dirty[] = false
    sh.diag_image.dirty[] = false
    return (redraw ? 1 : 0) | (chart ? 4 : 0) | (swchart ? 8 : 0) | (dchart ? 16 : 0)
end

_sym(x) = Symbol(String(x))
_parse(text, what) = Controllers._parse_number(string(text), what)

# Model tree and selection
br_select(id) = _with_shell(sh -> select!(sh.sess.view, String(id)))
br_select_add(id) = _with_shell(sh -> select_add!(sh.sess.view, String(id)))
br_select_range(id) = _with_shell(sh -> select_range!(sh.sess.view, [r.id for r in model_tree(sh.sess)],
                                                      String(id)))

# The nodes a tree action on `id` applies to: the whole selection when `id` is in it.
function _acted_nodes(sh, id)
    ids = selected_ids(sh.sess.view)
    return String(id) in ids ? ids : [String(id)]
end
# Component name of a component node (`nothing` for other nodes).
function _node_component(sess, id)
    kind, key = selection_kind(id)
    kind in (:comp, :part) && return Symbol(key)
    kind === :window && return Symbol(key)
    return nothing
end
_node_beam_names(sess, id) = [sess.beams.beams[][i].name for i in node_beams(sess.beams, id)]

# The NIF quads the tree's Split and Merge commands act on for the nodes acted on by a
# menu on `id`: how many are held as one beam, and how many as their four beams.
function br_node_quads(id)
    sh = _SHELL[]
    sh === nothing && return [0, 0]
    bt = sh.sess.beams
    idx = unique!(reduce(vcat, (node_beams(bt, n) for n in _acted_nodes(sh, id)); init = Int[]))
    isempty(idx) && return [0, 0]
    return [length(Controllers.quads_to_split(bt, idx)), length(Controllers.quads_to_merge(bt, idx))]
end
# What the tree's row menu needs for the rows a menu on `id` acts on, as integers: whether
# they are shown and included in the calculations (`br_node_flags`), the quads to split
# and to merge (`br_node_quads`), whether the row is a hohlraum piece that can take a
# role, its role (an index into `PIECE_ROLES`, from 0) and whether each role is allowed.
function br_node_menu(id)
    sh = _SHELL[]
    sh === nothing && return fill(0, 6 + length(PIECE_ROLES))
    f, q = br_node_flags(id), br_node_quads(id)
    kind, key = selection_kind(String(id))
    pr = kind in (:comp, :window) ? piece_roles(sh.sess.target, Symbol(key)) : nothing
    return Int[f[1], f[2], q[1], q[2], pr !== nothing && length(pr.roles) > 1,
               pr === nothing ? 0 : findfirst(==(pr.role), PIECE_ROLES) - 1,
               (pr !== nothing && r in pr.roles for r in PIECE_ROLES)...]
end
# Whether a node is shown in the 3D view and takes part in the calculations, for the
# tree's right-click menu: [shown, included].
function br_node_flags(id)
    sh = _SHELL[]
    sh === nothing && return [true, true]
    sess = sh.sess
    # For a user frame the switch shows whether its axes are drawn, in place of whether
    # it is shown.
    kind, key = selection_kind(String(id))
    kind === :frame && return [Symbol(key) in sess.view.triads[], true]
    # A part or an assembly is shown when all of its components are, and included when
    # none is left out.
    if kind in (:part, :assembly)
        cs = node_components(sess, String(id))
        g = sess.target.ghost[]
        out = g === nothing ? Set{Symbol}() : Set(c.name for c in g.components)
        return [all(c -> display_mode(sess.view, c) !== :invisible, cs), !any(in(out), cs)]
    end
    c = _node_component(sess, String(id))
    if c !== nothing
        return [display_mode(sess.view, c) !== :invisible, !(c in sess.target.excluded[])]
    end
    idx = node_beams(sess.beams, String(id))
    isempty(idx) && return [true, true]
    bs = sess.beams.beams[]
    return [all(i -> beam_drawn(sess.view, bs[i].name, 3), idx), all(i -> sess.beams.enabled[][i], idx)]
end
# Show or hide nodes: a component becomes opaque or invisible, and the 3ω light of a beam
# is drawn or not.
function br_node_show(id, on)
    _with_shell() do sh
        sess = sh.sess
        for n in _acted_nodes(sh, id)
            cs = node_components(sess, n)
            if !isempty(cs)
                foreach(c -> _show_component!(sess.view, c, Bool(on)), cs)
            else
                set_beams_drawn!(sess.view, _node_beam_names(sess, n), 3, Bool(on))
            end
        end
    end
end
# Showing an invisible component makes it opaque; another mode is kept.
_show_component!(v, c, on) = on ? (display_mode(v, c) === :invisible && set_display!(v, c, :opaque)) :
                                  set_display!(v, c, :invisible)
# Display mode of the components among the nodes acted on (the tree's or Properties'):
# a part's or an assembly's every component.
function br_set_display(id, mode)
    _with_shell() do sh
        m = mode isa Real ? DISPLAY_MODES[round(Int, mode) + 1] : _sym(mode)
        cs = unique!(reduce(vcat, (node_components(sh.sess, n) for n in _acted_nodes(sh, id)); init = Symbol[]))
        isempty(cs) || set_display!(sh.sess.view, cs, m)
    end
end
# A dot of the beam matrix: `row` "use" (takes part in the calculations) or "draw", of
# `harmonic` 1, 2 or 3, for the beams under the nodes acted on (every beam for "beams").
function br_beam_dot(id, row, h, on)
    _with_shell() do sh
        sess = sh.sess
        idx = sort!(unique!(reduce(vcat, (node_beams(sess.beams, n) for n in _acted_nodes(sh, id));
                                   init = Int[])))
        isempty(idx) && return
        h, on = Int(h), Bool(on)
        if String(row) == "use"
            set_participation!(sess.beams, idx, h, on)
        else
            set_beams_drawn!(sess.view, [sess.beams.beams[][i].name for i in idx], h, on)
        end
    end
end
function br_node_include(id, on)
    _with_shell() do sh
        sess = sh.sess
        comps = Symbol[]
        for n in _acted_nodes(sh, id)
            # An assembly stands for the parts of its subtree.
            c = _node_component(sess, n)
            kind, key = selection_kind(n)
            if kind === :assembly
                push!(comps, Symbol(key))
            elseif c !== nothing
                push!(comps, c)
            else
                idx = node_beams(sess.beams, n)
                isempty(idx) || set_set_enabled!(sess.beams, idx, Bool(on))
            end
        end
        # CAD parts put back are remeshed on workers (the target follows when ready)
        isempty(comps) || include_components!(sess, comps, Bool(on))
    end
end
br_toggle_expand(id) = _with_shell(sh -> toggle_expanded!(sh.sess.view, String(id)), _DIRTY_OBSERVED)
function br_tree_on(id, on)
    _with_shell() do sh
        sess = sh.sess
        kind, key = selection_kind(String(id))
        on = Bool(on)
        if kind in (:group, :quad, :beam)
            idx = String(id) == sess.view.selection[] ? selected_beams(sess) : node_beams(sess.beams, String(id))
            isempty(idx) || set_set_enabled!(sess.beams, idx, on)
        elseif kind === :beams
            Controllers._assign!(sess.beams.enabled, fill(on, length(sess.beams.beams[])))
        elseif kind in (:comp, :part, :window, :assembly)
            # every selected component when the box is on a selected one (a part's or an
            # assembly's every component)
            for n in _acted_nodes(sh, id), c in node_components(sess, n)
                _show_component!(sess.view, c, on)
            end
        end
    end
end
# Target
# File ▸ New target: only a starter, with the beams of its facility ("empty" for no parts).
function br_new_target(starter)
    _with_shell() do sh
        s = _sym(starter)
        new_target!(sh.sess, s === :empty ? :none : s) || return
        select!(sh.sess.view, "target")
        sh.sess.status[] = s === :empty ? "new target: empty (Insert adds parts)" :
                           "new target: " * STARTER_LABELS[s] * " · " * Controllers.beams_label(sh.sess.beams)
    end
end
# The frame of the selection, where Insert places a starter: a part's, a piece's (its
# part's piece frame), a window's, an assembly's; TCC otherwise.
function _selected_frame(sess)
    ts = sess.target
    kind, key = selection_kind(sess.view.selection[])
    kind in (:part, :comp, :window, :frame, :assembly) || return :tcc
    n = Symbol(key)
    part_index(ts, n) === nothing || return n
    for c in ts.target[].components
        c.name === n && c.part !== n && return piece_frame(n, c.part)
    end
    return :tcc
end
br_set_resolution(r) = _with_shell(sh -> set_resolution!(sh.sess.target, _sym(r)))
br_edit_mesh(key, text) = _with_shell(sh -> edit_mesh!(sh.sess.target, _sym(key), string(text)))
# Add a window to the selected (or first) hohlraum part's wall, and select it.
function br_add_window()
    _with_shell() do sh
        ts = sh.sess.target
        i = _window_part_index(sh.sess)
        part = i === nothing ? nothing : ts.parts[][i].name
        add_window!(ts; part) || return
        p = ts.parts[][_window_part_index(sh.sess)]
        select!(sh.sess.view, "window:" * String(piece_name(p.name, window_name(length(p.windows)))))
    end
end
# The hohlraum part of the window page (`_window_part_index`).
_window_part(sh) = (i = _window_part_index(sh.sess); i === nothing ? nothing : sh.sess.target.parts[][i].name)
function br_edit_window(row, key, text)
    _with_shell() do sh
        k = round(Int, row)
        # a rejected edit keeps the window: re-send its row so the field shows it again
        edit_window!(sh.sess.target, k, _sym(key), string(text); part = _window_part(sh)) ||
            (k <= length(sh.rows[:windowModel]) && (sh.models[:windowModel][k] = sh.rows[:windowModel][k]))
    end
end
function br_window_round(row, on)
    _with_shell() do sh
        k = round(Int, row)
        set_window_round!(sh.sess.target, k, Bool(on); part = _window_part(sh)) ||
            (k <= length(sh.rows[:windowModel]) && (sh.models[:windowModel][k] = sh.rows[:windowModel][k]))
    end
end
function br_window_open(row, on)
    _with_shell() do sh
        k = round(Int, row)
        set_window_open!(sh.sess.target, k, Bool(on); part = _window_part(sh)) ||
            (k <= length(sh.rows[:windowModel]) && (sh.models[:windowModel][k] = sh.rows[:windowModel][k]))
    end
end
function br_add_part(kind)
    _with_shell() do sh
        ts = sh.sess.target
        if add_part!(ts, _sym(kind))
            "target" in sh.sess.view.expanded[] || toggle_expanded!(sh.sess.view, "target")
            select!(sh.sess.view, "part:" * String(ts.parts[][end].name))
            sh.sess.status[] = "inserted " * String(ts.parts[][end].name) * " at TCC"
        else
            sh.sess.status[] = "not inserted: " * ts.error[]
        end
    end
end
function br_add_stl(url)
    _with_shell() do sh
        ts = sh.sess.target
        path = _url_to_path(String(url))
        n0 = length(ts.parts[])
        # one part, or one per named solid of a multi-body export, in one assembly
        if add_stl_parts!(ts, path)
            new = ts.parts[][n0+1:end]
            n = count(!is_assembly, new)
            if n == 1
                name = String(new[1].name)
                select!(sh.sess.view, "part:" * name)
                sh.sess.status[] = "added $name from $(basename(path)): set its scale (mm per file unit) and place in Properties"
            else
                a = new[end]
                id = part_node_id(a)
                id in sh.sess.view.expanded[] || toggle_expanded!(sh.sess.view, id)
                select!(sh.sess.view, id)
                sh.sess.status[] = "added $n parts from $(basename(path)), one per solid, in assembly $(a.name): " *
                                   "set their scale (mm per file unit) and place them in Properties"
            end
        end
    end
end
function _resend_part_fields!(sh)
    rows = sh.rows[:partFieldModel]
    for i in eachindex(rows)
        sh.models[:partFieldModel][i] = rows[i]
    end
end
# The selected parts (indices), the one shown in Properties first: edits apply to all.
function _selected_parts(sh)
    ts = sh.sess.target
    idx = Int[]
    for id in selected_ids(sh.sess.view)
        kind, key = selection_kind(id)
        kind in (:part, :frame, :assembly) || continue
        i = part_index(ts, Symbol(key))
        i === nothing || push!(idx, i)
    end
    return idx
end
function br_edit_part(key, text)
    _with_shell() do sh
        ts = sh.sess.target
        idx = _selected_parts(sh)
        isempty(idx) && return
        k = _sym(key)
        if k === :name                                  # names are per part
            i = idx[1]
            edit_part!(ts, i, k, string(text)) ?
                select!(sh.sess.view, part_node_id(ts.parts[][i])) : _resend_part_fields!(sh)
            return
        end
        ok = true
        for i in idx
            haskey(ts.parts[][i].values, k) || continue  # a field of other kinds
            ok &= edit_part!(ts, i, k, string(text))
        end
        ok || _resend_part_fields!(sh)
    end
end
# A hohlraum part's section list: break k's height or radius, a break added in the longest
# section, an inner break removed (the first selected part; k from 0 in QML).
function _section_edit(f)
    _with_shell() do sh
        idx = _selected_parts(sh)
        isempty(idx) && return
        f(sh.sess.target, idx[1]) && return
        for k in ("partSectionZ", "partSectionR")       # a rejected edit shows the table again
            haskey(sh.shown, k) && (sh.app[k] = copy(sh.shown[k]))
        end
    end
end
br_section_edit(k, field, text) =
    _section_edit((ts, i) -> edit_hohlraum_break!(ts, i, Int(k) + 1, _sym(field), string(text)))
br_section_add() = _section_edit((ts, i) -> add_hohlraum_break!(ts, i))
br_section_remove(k) = _section_edit((ts, i) -> remove_hohlraum_break!(ts, i, Int(k) + 1))
function br_part_flag(flag, on)
    _with_shell() do sh
        ts = sh.sess.target
        for i in _selected_parts(sh)
            set_part_flag!(ts, i, _sym(flag), Bool(on))
        end
    end
end
# A frame's name in the menus: TCC's label, a user frame's, a part's name, or a piece
# frame's (part ▸ piece). `ts`: the design the names belong to (the window's by default).
function _frame_label(f::Symbol, ts = (sh = _SHELL[]; sh === nothing ? nothing : sh.sess.target))
    haskey(FRAME_LABELS, f) && return FRAME_LABELS[f]
    s = String(f)
    i = findfirst('.', s)
    if i === nothing
        k = ts === nothing ? nothing : part_index(ts, f)
        p = k === nothing ? nothing : ts.parts[][k]
        p !== nothing && is_port_frame(p) && return "Port frame: " * s * " (" * p.object * ")"
        p !== nothing && p.kind === :frame && return "Frame: " * s
        p !== nothing && p.kind === :assembly && return "Assembly: " * s
        return "Part: " * s
    end
    piece = s[i+1:end]
    return s[1:i-1] * " ▸ " * get(_PIECE_FRAME_LABELS, piece, replace(piece, "_" => " "))
end
const _PIECE_FRAME_LABELS = Dict("leh_top" => "top LEH and cap", "leh_bottom" => "bottom LEH and cap")
# A part's frame and its pieces' (`piece_frame_defs`), in order.
_with_piece_frames(ps) = Symbol[n for q in ps for n in (q.name, first.(piece_frame_defs(q))...)]
# The parts in menu order: the user frames first, then the parts.
_frames_first(ps) = [filter(q -> q.kind === :frame, ps); filter(q -> q.kind !== :frame, ps)]
# Frames a part can be placed in: TCC, the user frames, then the other parts and their pieces.
function _part_frame_choices(sess, p)
    return [:tcc; _with_piece_frames(q for q in _frames_first(sess.target.parts[]) if q.name !== p.name)]
end
# Frames to move along or aim in: TCC, the user frames, then every part and its pieces.
function _move_frame_choices(sess)
    return [:tcc; _with_piece_frames(_frames_first(sess.target.parts[]))]
end
# A user frame's line under its name: where it is in the chamber, and what uses it.
function _frame_info(sess, p)
    F = part_frames(sess.target)[p.name]
    _, θ, φ = coords_of(F.R[:, 3], :spherical)
    nb = count(==(p.name), sess.beams.aim_frame[])
    nu = length(frame_users(sess.target, p.name))
    use = join(filter(!isempty, [nb == 0 ? "" : nb == 1 ? "1 beam aims in it" : "$nb beams aim in it",
                                 nu == 0 ? "" : nu == 1 ? "1 part or frame placed in it" : "$nu parts or frames placed in it"]), " · ")
    return @sprintf("origin (%.4g, %.4g, %.4g) mm · z at θ %.4g°, φ %.4g°", F.origin..., θ, φ) *
           (isempty(use) ? "" : " · " * use)
end
# An assembly's line under its name: its members, where its frame is, the beams aimed in it.
function _assembly_info(sess, p)
    ps = sess.target.parts[]
    nk, na = length(assembly_parts(ps, p.name)), count(is_assembly, assembly_children(ps, p.name))
    F = part_frames(sess.target)[p.name]
    nb = count(==(p.name), sess.beams.aim_frame[])
    what = join(filter(!isempty, [nk == 1 ? "1 part" : "$nk parts",
                                  na == 0 ? "" : na == 1 ? "1 assembly" : "$na assemblies",
                                  nb == 0 ? "" : nb == 1 ? "1 beam aims in it" : "$nb beams aim in it"]), " · ")
    return what * @sprintf(" · origin (%.4g, %.4g, %.4g) mm", F.origin...)
end
# The model tree's assembly commands. New assembly: of the selected parts and assemblies
# (an empty one with none), selected and open.
function br_new_assembly()
    _with_shell() do sh
        sess = sh.sess
        ts = sess.target
        members = [ts.parts[][i].name for i in _selected_parts(sh) if ts.parts[][i].kind !== :frame]
        add_assembly!(ts; members) || return
        a = last(assemblies(ts))
        id = part_node_id(a)
        for e in ("target", id)
            e in sess.view.expanded[] || toggle_expanded!(sess.view, e)
        end
        select!(sess.view, id)
        sess.status[] = "added assembly $(a.name)" * (isempty(members) ? "" : " of " * join(String.(members), ", "))
    end
end
# Move the selected parts and assemblies into the assembly of `assemblyChoices`, where 0 is
# the top level.
function br_move_into_assembly(index)
    _with_shell() do sh
        sess = sh.sess
        ts = sess.target
        k = Int(index)
        as = assemblies(ts)
        0 <= k <= length(as) || return
        a = k == 0 ? TOP_LEVEL : as[k].name
        names = [ts.parts[][i].name for i in _selected_parts(sh)]
        isempty(names) && return
        move_into_assembly!(ts, names, a) || (sess.status[] = ts.error[]; return)
        a === TOP_LEVEL || (id = "assembly:" * String(a); id in sess.view.expanded[] || toggle_expanded!(sess.view, id))
        sess.status[] = "moved " * join(String.(names), ", ") * (a === TOP_LEVEL ? " to the top level" : " into $a")
    end
end
# Ungroup the selected assembly; its members stay in its parent.
function br_ungroup_assembly()
    _with_shell() do sh
        sess = sh.sess
        kind, key = selection_kind(sess.view.selection[])
        kind === :assembly || return
        ungroup_assembly!(sess.target, Symbol(key)) || (sess.status[] = sess.target.error[]; return)
        sess.status[] = "ungrouped $key"
    end
end
# Rename the selected part, assembly or frame (what refers to it follows).
function br_rename_node(text)
    _with_shell() do sh
        sess = sh.sess
        kind, key = selection_kind(sess.view.selection[])
        kind in (:part, :assembly, :frame) || return
        i = part_index(sess.target, Symbol(key))
        i === nothing && return
        if edit_part!(sess.target, i, :name, string(text))
            select!(sess.view, part_node_id(sess.target.parts[][i]))
        else
            sess.status[] = "not renamed: " * sess.target.error[]
        end
    end
end
# Move the selected part to another frame, keeping it where it is (index into
# `frameLabels`).
function br_part_frame(index)
    _with_shell() do sh
        ts = sh.sess.target
        idx = _selected_parts(sh)
        isempty(idx) && return
        frame = _part_frame_choices(sh.sess, ts.parts[][idx[1]])[Int(index)+1]
        ok = true
        for i in idx
            ok &= set_part_frame!(ts, i, frame)
        end
        ok || _set!(sh, "partFrame", -1)
    end
end

# The frames a placement step of part `p` can refer to (every frame but the part frame),
# with "" first for none where the step allows it, and their menu labels.
function _step_ref_options(sess, p, st)
    frames = [f for f in _move_frame_choices(sess) if f !== p.name]
    labels = _frame_label.(frames)
    st isa Roll && return [Symbol(""); frames], ["(relative: from the current roll)"; labels]
    st isa Offset && return [Symbol(""); frames], ["(the part axes)"; labels]
    return frames, labels
end

const _STEP_AXES = ("x", "y", "z")

# What step `k` fixes, from the placement report `r`, and how far the part ends from it.
function _step_note(r, k, st)
    t, q = r.fixed[k]
    plural(n, w) = "$n $w" * (n == 1 ? "" : "s")
    note = Briareus.GUIAPI.is_relative_step(st) ? "relative move" :
           k in r.redundant ? "redundant: fixes nothing new" :
           "fixes " * join(filter(!isempty, [t > 0 ? plural(t, "translation") : "",
                                            q > 0 ? plural(q, "rotation") : ""]), ", ")
    mm, deg = r.residuals[k]
    off = filter(!isempty, [mm > 1e-9 ? _num(mm) * " mm" : "", deg > 1e-7 ? _num(deg) * "°" : ""])
    isempty(off) || (note *= "; off by " * join(off, ", ") * " after later moves")
    return note
end

# The inspector rows of part `i`'s placement steps, and a note on what stays free.
function _step_rows(sess, i)
    ts = sess.target
    p = ts.parts[][i]
    rows = StepRow[]
    isempty(p.steps) && return rows, ""
    r = try
        part_placement(ts, i)
    catch err
        err isa ArgumentError || rethrow()
        nothing
    end
    for (k, st) in enumerate(p.steps)
        kind = step_kind(st)
        push!(rows, StepRow(k, "head", "", "$k · " * PART_STEP_LABELS[kind], "",
                            r === nothing ? "" : _step_note(r, k, st), "", 0))
        if hasfield(typeof(st), :ref)
            opts, labels = _step_ref_options(sess, p, st)
            cur = something(st.ref, Symbol(""))
            label = st isa Offset ? "Along the axes of" : st isa Roll ? "From the axes of" : "Frame"
            push!(rows, StepRow(k, "choice", "ref", label, "", "", join(labels, "|"),
                                something(findfirst(==(cur), opts), 1) - 1))
        end
        if hasfield(typeof(st), :point)
            for (c, key) in enumerate(("px", "py", "pz"))
                push!(rows, StepRow(k, "number", key, "Point " * _STEP_AXES[c] * " (in the frame)", "mm",
                                    _num(st.point[c]), "", 0))
            end
        end
        if hasfield(typeof(st), :delta)
            for (c, key) in enumerate(("dx", "dy", "dz"))
                push!(rows, StepRow(k, "number", key, "Offset " * _STEP_AXES[c], "mm", _num(st.delta[c]), "", 0))
            end
        end
        if hasfield(typeof(st), :axis) && !(st isa Roll && st.ref === nothing)
            j = findfirst(a -> st.axis == Briareus.GUIAPI.step_axis(Symbol(a)), _STEP_AXES)
            choices = j === nothing ? [_STEP_AXES...; "(" * join(_num.(st.axis), ", ") * ")"] : [_STEP_AXES...]
            label = st isa Rotate ? "About its axis" : st isa Roll ? "Reference axis" : "Axis of the frame"
            push!(rows, StepRow(k, "choice", "axis", label, "", "", join(choices, "|"),
                                something(j, length(choices)) - 1))
        end
        st isa Union{OnAxis,Distance} &&
            push!(rows, StepRow(k, "number", "distance", st isa OnAxis ? "Distance along the axis" : "Distance",
                                st isa OnAxis ? "mm; blank: the nearest point of the axis" : "mm",
                                st.distance === nothing ? "" : _num(st.distance), "", 0))
        hasfield(typeof(st), :angle) &&
            push!(rows, StepRow(k, "number", "angle", "Angle", "°", _num(st.angle), "", 0))
        st isa Facing && push!(rows, StepRow(k, "choice", "side", "z axis", "", "",
                                             "toward the point|away from the point", st.side === :front ? 0 : 1))
        hasfield(typeof(st), :flip) &&
            push!(rows, StepRow(k, "check", "flip", "Opposite direction", "", "", "", st.flip ? 1 : 0))
    end
    r === nothing && return rows, ""
    nt, nr = length(r.free_translation), length(r.free_rotation)
    note = nt + nr == 0 ? "The steps fix the part's place and orientation." :
           "Free: $nt translation" * (nt == 1 ? "" : "s") * ", $nr rotation" * (nr == 1 ? "" : "s") *
           " (the placement fields above set them)."
    return rows, note
end

# Edits of the shown part's placement steps (`k` from 1); a rejected one shows the rows again.
function _step_action(f, sh)
    ts = sh.sess.target
    idx = _selected_parts(sh)
    isempty(idx) && return
    f(ts, idx[1]) === false || return
    rows = sh.rows[:stepModel]
    for i in eachindex(rows)
        sh.models[:stepModel][i] = rows[i]
    end
end
br_step_add(index) = _with_shell(sh -> _step_action((ts, i) -> add_part_step!(ts, i, PART_STEP_KINDS[Int(index)+1]), sh))
br_step_remove(k) = _with_shell(sh -> _step_action((ts, i) -> remove_part_step!(ts, i, Int(k)), sh))
br_step_move(k, by) = _with_shell(sh -> _step_action((ts, i) -> move_part_step!(ts, i, Int(k), Int(by)), sh))
br_step_edit(k, key, text) =
    _with_shell(sh -> _step_action((ts, i) -> edit_part_step!(ts, i, Int(k), _sym(key), string(text)), sh))
function br_step_choice(k, key, index)
    _with_shell() do sh
        _step_action(sh) do ts, i
            st = ts.parts[][i].steps[Int(k)]
            j, key = Int(index) + 1, _sym(key)
            value = key === :ref ? String(first(_step_ref_options(sh.sess, ts.parts[][i], st))[j]) :
                    key === :axis ? (j <= 3 ? _STEP_AXES[j] : return true) :
                    key === :side ? ("front", "back")[j] :
                    key === :flip ? j == 2 : throw(ArgumentError("unknown step choice $key"))
            edit_part_step!(ts, i, Int(k), key, value)
        end
    end
end

# The move row: the frame (index into `moveFrames`) and axis (0 x, 1 y, 2 z) it moves
# along, and the step as text.
br_move_frame(index) = _with_shell(sh -> (sh.shown["_moveFrameName"] = _move_frame_choices(sh.sess)[Int(index)+1]))
# The frame the move row moves along: the chosen one while the design has it, else the chamber's.
function _move_frame(sh)
    f = get(sh.shown, "_moveFrameName", :tcc)
    return f in _move_frame_choices(sh.sess) ? f : :tcc
end
br_move_axis(index) = _with_shell(sh -> _set!(sh, "moveAxis", Int(index)))
# (the placement triad's arrows snap to the same step)
br_move_step(text) = _with_shell() do sh
    sh.scene.triad_step[] = Controllers._parse_number(string(text), "the step"; positive = true)
    _set!(sh, "moveStep", String(strip(string(text))))
end
# Move the selected part by ± the step along the chosen axis (`sign` +1 or -1): one edit.
function br_nudge_part(sign)
    _with_shell() do sh
        ts = sh.sess.target
        idx = _selected_parts(sh)
        isempty(idx) && return
        Δ = Float64(sign) * Controllers._parse_number(sh.shown["moveStep"], "the step"; positive = true)
        frame = _move_frame(sh)
        axis = SWEEP_AXES[sh.shown["moveAxis"]+1]
        for i in idx
            nudge_part!(ts, i, Δ; frame, axis)
        end
        sh.sess.status[] = @sprintf("moved %s by %+g mm along %s of %s", ts.parts[][idx[1]].name, Δ, axis,
                                    _frame_label(frame))
    end
end
# Set up a sweep of the selected part along the move row's axis, and show it.
function br_sweep_part()
    _with_shell() do sh
        idx = _selected_parts(sh)
        isempty(idx) && return
        p = sh.sess.target.parts[][idx[1]]
        set_sweep!(sh.sess.sweep; subject = p.name, mode = :translate,
                   frame = _move_frame(sh),
                   axis = SWEEP_AXES[sh.shown["moveAxis"]+1])
        select!(sh.sess.view, "task:sweep")
        request_tab("sweep")
    end
end
# The selected part's role in the analyses (index into `roleLabels`); a refused role
# shows the current role of the part again.
function br_part_role(index)
    _with_shell() do sh
        idx = _selected_parts(sh)
        isempty(idx) && return
        ts = sh.sess.target
        role = PART_ROLES[Int(index)+1]
        name = ts.parts[][idx[1]].name
        # once given (a CAD capsule's sphere may be fitted on a worker first): select the
        # part holding the role, the fitted sphere for a CAD capsule
        done = function ()
            o = observed_name(ts)
            i = role === :observed && o !== nothing ? part_index(ts, o) : nothing
            p = i === nothing ? ts.parts[][something(part_index(ts, name), idx[1])] : ts.parts[][i]
            select!(sh.sess.view, "part:" * String(p.name))
            sh.sess.status[] = "$(p.name): " * lowercase(PART_ROLE_LABELS[part_role(p)])
            _mark!(sh.dirty)
        end
        if set_part_role!(ts, idx[1], role; done)
            name in ts.fitting[] && (sh.sess.status[] = "$name: fitting the capsule's sphere…")
        else
            delete!(sh.shown, "partRole")
        end
    end
end
# The selected piece's role (index into `pieceRoleLabels`); a refused role shows the
# current role of the piece again.
function br_piece_role(index)
    _with_shell() do sh
        name = selected_component(sh.sess)
        name === nothing && return
        role = PIECE_ROLES[Int(index)+1]
        if set_piece_role!(sh.sess.target, name, role)
            sh.sess.status[] = "$name: " * lowercase(PIECE_ROLE_LABELS[role])
            _mark!(sh.dirty)
        else
            sh.sess.status[] = sh.sess.target.error[]
        end
        delete!(sh.shown, "pieceRole")
    end
end
# The foil the foil analyses take (index into `foilNames`).
function br_active_foil(index)
    _with_shell() do sh
        fs = foil_names(sh.sess.target)
        k = Int(index) + 1
        1 <= k <= length(fs) || return
        set_active_foil!(sh.sess.target, fs[k]) &&
            (sh.sess.status[] = "the foil analyses take $(fs[k])")
        delete!(sh.shown, "activeFoil")
    end
end
# The cavity whose wall temperature the Drive panel shows (index into `cavityNames`).
function br_drive_cavity(index)
    _with_shell() do sh
        res = sh.sess.drive.result[]
        k = Int(index) + 1
        (res === nothing || !(1 <= k <= length(res.cavities))) && return
        set_drive_cavity!(sh.sess, res.cavities[k].name)
        delete!(sh.shown, "cavityIndex")
    end
end
function br_part_coords(index)
    _with_shell() do sh
        for i in _selected_parts(sh)
            set_part_coords!(sh.sess.target, i, PART_COORDINATES[Int(index)+1])
        end
    end
end
# The Frames folder's page and the transform dialog: the ports a port frame is made from
# (once per beam table) and the transform of the current frames.
function _refresh_frames!(sh)
    bs = sh.sess.beams.beams[]
    if get(sh.shown, "_portBeams", nothing) !== bs
        sh.shown["_portBeams"] = bs
        _set!(sh, "portChoices", [@sprintf("%s · θ %.4g° φ %.4g°", q.name, q.theta, q.phi)
                                  for q in port_choices(sh.sess.beams)])
    end
    _set!(sh, "frameCount", length(user_frames(sh.sess.target)))
    haskey(sh.shown, "_tf") && _refresh_transform!(sh)
    return
end
# Frames: a user frame in the selected frame (TCC with none), a port frame of a port of
# the beam table (index into `portChoices`), a frame's axes drawn or not.
function br_add_frame()
    _with_shell() do sh
        sess = sh.sess
        frame = _selected_frame(sess)
        add_frame!(sess.target; frame) || return
        p = sess.target.parts[][end]
        "frames" in sess.view.expanded[] || toggle_expanded!(sess.view, "frames")
        select!(sess.view, part_node_id(p))
        sess.status[] = "added frame $(p.name) " * (frame === :tcc ? "at TCC" : "in " * _frame_label(frame)) *
                        ": place it in Properties"
    end
end
function br_add_port_frame(index)
    _with_shell() do sh
        sess = sh.sess
        ports = port_choices(sess.beams)
        k = Int(index) + 1
        1 <= k <= length(ports) || return
        q = ports[k]
        add_port_frame!(sess.target, q.name, q.theta, q.phi) || return
        p = sess.target.parts[][end]
        "frames" in sess.view.expanded[] || toggle_expanded!(sess.view, "frames")
        select!(sess.view, part_node_id(p))
        sess.status[] = @sprintf("added port frame %s: z toward port %s (θ %.4g°, φ %.4g°), origin at TCC",
                                 p.name, q.name, q.theta, q.phi)
    end
end
# Draw the axes of the frame of node `id` (a user frame; "" for the selection's).
function br_set_triad(id, on)
    _with_shell() do sh
        kind, key = selection_kind(isempty(String(id)) ? sh.sess.view.selection[] : String(id))
        kind === :frame && set_triad!(sh.sess.view, Symbol(key), Bool(on))
    end
end
# The coordinate transform dialog: its state in `shown["_tf"]`, the result in `tf*`.
function _tf_state(sh)
    get!(sh.shown, "_tf") do
        Dict{String,Any}("from" => :tcc, "to" => :tcc, "input" => :cartesian, "output" => :cartesian,
                         "direction" => false, "c" => [0.0, 0.0, 0.0])
    end
end
function _refresh_transform!(sh)
    sess, st = sh.sess, _tf_state(sh)
    choices = _move_frame_choices(sess)
    for k in ("from", "to")
        st[k] in choices || (st[k] = :tcc)
    end
    _set!(sh, "tfFrames", [_frame_label(f) for f in choices])
    _set!(sh, "tfFrom", findfirst(==(st["from"]), choices) - 1)
    _set!(sh, "tfTo", findfirst(==(st["to"]), choices) - 1)
    _set!(sh, "tfInput", findfirst(==(st["input"]), PART_COORDINATES) - 1)
    _set!(sh, "tfOutput", findfirst(==(st["output"]), PART_COORDINATES) - 1)
    _set!(sh, "tfDirection", st["direction"])
    _set!(sh, "tfIn", [_num(x) for x in st["c"]])
    try
        out = transform_coordinates(sess.target, st["from"], st["to"], st["c"]; input = st["input"],
                                    output = st["output"], direction = st["direction"])
        _set!(sh, "tfOut", [@sprintf("%.6g", x) for x in out])
        _set!(sh, "tfError", "")
    catch err
        err isa ArgumentError || rethrow()
        _set!(sh, "tfOut", ["–", "–", "–"])
        _set!(sh, "tfError", Controllers._errmsg(err))
    end
    return
end
# Set a field of the transform: "from", "to" (indices into `tfFrames`), "input",
# "output" (indices into the coordinate systems), "direction" (Bool), "c1"–"c3" (text).
function br_tf_set(field, value)
    _with_shell() do sh
        st, f = _tf_state(sh), String(field)
        if f in ("from", "to")
            st[f] = _move_frame_choices(sh.sess)[Int(value)+1]
        elseif f in ("input", "output")
            st[f] = PART_COORDINATES[Int(value)+1]
        elseif f == "direction"
            st[f] = Bool(value)
        else
            k = parse(Int, f[2:end])
            x = tryparse(Float64, strip(string(value)))
            x === nothing ? _set!(sh, "tfError", "not a number: $(value)") : (st["c"][k] = x)
        end
        _refresh_transform!(sh)
    end
end
# Open the transform from the selected frame (TCC with none) to TCC.
function br_tf_open()
    _with_shell() do sh
        st = _tf_state(sh)
        st["from"] = _selected_frame(sh.sess)
        _refresh_transform!(sh)
    end
end
# Open the transform dialog from the selected frame.
br_open_transform() = request_popup("transform")
# A user frame made from the transform (`add_frame_from_transform!`), selected.
function br_tf_create()
    _with_shell() do sh
        sess, st = sh.sess, _tf_state(sh)
        add_frame_from_transform!(sess.target, st["from"], st["to"], st["c"]; input = st["input"],
                                  direction = st["direction"]) || (_set!(sh, "tfError", sess.target.error[]); return)
        p = sess.target.parts[][end]
        "frames" in sess.view.expanded[] || toggle_expanded!(sess.view, "frames")
        select!(sess.view, part_node_id(p))
        sess.status[] = "added frame $(p.name) in " * _frame_label(st["to"])
        _refresh_transform!(sh)
    end
end
# Remove every selected part, assembly, frame and window, as one edit (`remove_nodes!`).
function br_remove_selected()
    _with_shell() do sh
        sess = sh.sess
        r = Controllers.remove_nodes!(sess)
        if isempty(r.removed)
            sess.status[] = r.ok ? "nothing to remove: select parts, assemblies, frames or windows (a piece can be hidden)" :
                            "not removed: " * sess.target.error[]
            return
        end
        n = length(r.removed)
        what = n <= 4 ? join(r.removed, ", ") : "$n items"
        kept = isempty(r.replaced) ? "" :
               " · " * join(String.(r.replaced), ", ") * (length(r.replaced) == 1 ? " keeps its place" : " keep their places")
        sess.status[] = "removed " * what * kept * (r.ok ? "" : " · " * sess.target.error[])
    end
end
# Beams
br_load_preset(name) = _with_shell(sh -> load_preset!(sh.sess, _sym(name)))
br_set_grouping(by) = _with_shell(sh -> set_grouping!(sh.sess.beams, _sym(by)))
function br_edit_beam(row, field, text)
    _with_shell() do sh
        k = round(Int, row)
        # a rejected edit leaves the beams unchanged: re-send the row so the field
        # shows the kept value again
        edit_beams!(sh.sess.beams, k, _sym(field), string(text)) ||
            (k <= length(sh.rows[:beamModel]) && (sh.models[:beamModel][k] = sh.rows[:beamModel][k]))
    end
end
br_set_beam_on(row, on) = _with_shell(sh -> set_enabled!(sh.sess.beams, round(Int, row), Bool(on)))
# Edit every selected beam (a group, quad, beam or Ctrl-selected set) at once.
function br_edit_node(field, text)
    _with_shell() do sh
        idx = selected_beams(sh.sess)
        isempty(idx) && return
        edit_beam_set!(sh.sess.beams, idx, _sym(field), string(text)) ||
            (delete!(sh.shown, "bn_" * String(field)))       # re-send the kept value
    end
end
function br_set_cluster(name)
    _with_shell() do sh
        sess = sh.sess
        idx = selected_beams(sess)
        isempty(idx) && return
        if set_cluster!(sess.beams, idx, string(name))
            g = strip(string(name))
            "group:" * g in sess.view.expanded[] || toggle_expanded!(sess.view, "group:" * g)
            select!(sess.view, "group:" * g)
            sess.status[] = "moved $(length(idx)) beam$(length(idx) == 1 ? "" : "s") to group $g"
        end
    end
end
function br_rename_cluster(new)
    _with_shell() do sh
        kind, key = selection_kind(sh.sess.view.selection[])
        kind === :group && rename_cluster!(sh.sess.beams, key, string(new)) &&
            select!(sh.sess.view, "group:" * strip(string(new)))
    end
end
_aim_choices(sess) = _move_frame_choices(sess)
# Hold the selected beams' aims in a frame, where they are (index into `aimChoices`).
function br_aim_frame(index)
    _with_shell() do sh
        idx = selected_beams(sh.sess)
        isempty(idx) && return
        f = _aim_choices(sh.sess)[Int(index)+1]
        set_aim_frame!(sh.sess, idx, f)
        sh.sess.status[] = f === :tcc ? "the beams aim in the chamber frame" :
                           "the beams aim in the frame of " * _frame_label(f) * " and follow it"
    end
end
# Aim the selected beams at their frame's point (a part's centroid, a template frame's origin).
function br_aim_at_frame()
    _with_shell() do sh
        idx = selected_beams(sh.sess)
        isempty(idx) && return
        fs = unique(sh.sess.beams.aim_frame[][idx])
        length(fs) == 1 || throw(ArgumentError("the beams aim in different frames: choose one first"))
        set_aim_frame!(sh.sess, idx, fs[1]; keep = false)
        sh.sess.status[] = "the beams aim at the center of " * _frame_label(fs[1])
    end
end
function br_set_facility(f)
    _with_shell() do sh
        s = sh.sess
        f = _sym(f)
        f === s.beams.facility[] && return
        set_facility!(s.beams, f; target = s.target)
        s.status[] = "facility: " * FACILITY_LABELS[f] * " · " *
                     (isempty(s.beams.beams[]) ? "no built-in beams: import a beam table (File)" : beams_summary(s.beams))
    end
end
# The beams under the nodes a tree command on `id` acts on ("": the selected beams).
function _command_beams(sh, id)
    isempty(String(id)) && return selected_beams(sh.sess)
    return sort!(unique!(reduce(vcat, (node_beams(sh.sess.beams, n) for n in _acted_nodes(sh, id)); init = Int[])))
end
# The selection after quads were split (`split`) or merged: the quad node in place
# of its beam's (or its beams'); other nodes as they were.
function _reselect_quads!(sh, ids, quads, split::Bool)
    sess = sh.sess
    bt = sess.beams
    group(q) = (i = findfirst(b -> b.quad == q, bt.beams[]); i === nothing ? "" : String(bt.cluster[][i]))
    out = String[]
    for id in ids
        kind, key = selection_kind(id)
        q = kind === :beam && split && key in quads ? key :
            kind === :quad && !split && first(split_at(key)) in quads ? first(split_at(key)) : nothing
        new = q === nothing ? id : split ? "quad:" * q * "@" * group(q) : "beam:" * q
        (new in out || isempty(node_beams(bt, new))) || push!(out, new)
    end
    isempty(out) && return select!(sess.view, "beams")
    select!(sess.view, out[1])
    foreach(x -> select_add!(sess.view, x), out[2:end])
    return
end
split_at(key) = split(key, '@'; limit = 2)
# Split the NIF quads held as one beam under the nodes acted on (`id`; "": the selected
# beams) into their four beams, or merge those held as four, in one edit.
function br_split_quads(id)
    _with_shell() do sh
        sess = sh.sess
        bt = sess.beams
        idx = _command_beams(sh, id)
        quads = String[bt.beams[][i].quad for i in Controllers.quads_to_split(bt, idx)]
        ids = selected_ids(sess.view)
        n = Controllers.split_quads!(bt, idx)
        n == 0 && (sess.status[] = bt.error[]; return)
        _reselect_quads!(sh, ids, quads, true)
        sess.status[] = (n == 1 ? "split $(quads[1]) into its four beams" : "split $n quads into their beams") *
                        ": the quad stays the unit of pointing and power · compute the drive again"
    end
end
function br_merge_quads(id)
    _with_shell() do sh
        sess = sh.sess
        bt = sess.beams
        idx = _command_beams(sh, id)
        quads = Controllers.quads_to_merge(bt, idx)
        ids = selected_ids(sess.view)
        n = Controllers.merge_quads!(bt, idx)
        n == 0 && (sess.status[] = bt.error[]; return)
        _reselect_quads!(sh, ids, quads, false)
        sess.status[] = (n == 1 ? "merged $(quads[1]) into one beam" : "merged $n quads into one beam each") *
                        " · compute the drive again"
    end
end
# Session files and undo
br_undo() = _with_shell(sh -> undo!(sh.sess) || (sh.sess.status[] = "nothing to undo"))
br_redo() = _with_shell(sh -> redo!(sh.sess) || (sh.sess.status[] = "nothing to redo"))
function br_open_session(url)
    _with_shell() do sh
        path = _url_to_path(String(url))
        load_session!(sh.sess, path)
        sh.session_path = path
        select!(sh.sess.view, "target")
    end
end
# An empty url saves to the current session file; `refs` saves the CAD files by path only.
function br_save_session(url, refs)
    _with_shell() do sh
        path = isempty(String(url)) ? sh.session_path : _url_to_path(String(url))
        isempty(path) && throw(ArgumentError("choose a file to save the session to"))
        save_session(path, sh.sess; embed = !Bool(refs) && !endswith(path, ".toml"))
        sh.session_path = path
    end
end
# For CAD copies saved with a session, read the changed sources or keep the copies. The
# prompt is answered either way, and a file the target rejects stays a saved copy.
br_reload_cad_sources() = _with_shell(sh -> (reload_cad_sources!(sh.sess); sh.sess.cad_changed[] = String[]))
br_keep_cad_copies() = _with_shell(sh -> (sh.sess.cad_changed[] = String[]))
# The selected part's CAD file read from `newfile` (every part reading that file moves).
function _relink_part(sh, newfile)
    i = part_index(sh.sess.target, Symbol(selection_kind(sh.sess.view.selection[])[2]))
    i === nothing && return
    file = sh.sess.target.parts[][i].file
    relink_cad!(sh.sess, file, newfile) &&
        (sh.sess.status[] = "read the parts of $(basename(cad_source(sh.sess, file))) from $newfile")
end
br_relink_part(url) = _with_shell(sh -> _relink_part(sh, _url_to_path(String(url))))
function br_part_from_source()
    _with_shell() do sh
        i = part_index(sh.sess.target, Symbol(selection_kind(sh.sess.view.selection[])[2]))
        i === nothing || _relink_part(sh, cad_source(sh.sess, sh.sess.target.parts[][i].file))
    end
end
br_import_beams(url) = _with_shell(sh -> (import_beams!(sh.sess.beams, _url_to_path(String(url)));
                                          sh.sess.status[] = "imported " * basename(sh.sess.beams.path[])))
br_export_beams(url) = _with_shell(sh -> (export_beams(sh.sess.beams, _url_to_path(String(url)));
                                          sh.sess.status[] = "exported " * basename(sh.sess.beams.path[])))
function _visrad_io(f, url, verb)
    _with_shell() do sh
        path = _url_to_path(String(url))
        f(sh.sess.beams, path)
        sh.sess.status[] = verb * " " * basename(path)
    end
end
br_import_sst(url) = _visrad_io(import_visrad_sst!, url, "applied VisRad pointing from")
br_import_summary(url) = _visrad_io(import_visrad_summary!, url, "applied the VisRad beam summary")
br_export_sst(url) = _visrad_io(export_visrad_sst, url, "wrote VisRad pointing to")
br_export_summary(url) = _visrad_io(export_visrad_summary, url, "wrote the VisRad beam summary to")
# Start an analysis (a `TASKS` id: the model tree's Run, the Compute commands); the sweep
# shows its panel.
br_run_task(id) = _with_shell(sh -> (k = _sym(id); run_task!(sh.sess, k); k === :sweep && request_tab("sweep")))
# Drive
br_build() = br_run_task("drive")
br_set_view(v) = _with_shell(sh -> set_view!(sh.sess.drive, _sym(v)))
# Leave the x-ray image for the view shown before it.
br_close_image() = _with_shell(sh -> set_view!(sh.sess.drive, get(sh.shown, "_viewBeforeImage", :scene)))
br_side_view(on) = _with_shell(sh -> set_side_view!(sh.sess.drive, Bool(on)))
# An x-ray image seen along the 3D view's direction, framing what the view frames.
function br_xray_image()
    _with_shell() do sh
        c = sh.scene.camera
        f, _, u = camera_frame(c)
        width = 2 * norm(c.lookat - c.eye) * tand(c.fov / 2)
        current = Controllers._drive_current(sh.sess)
        v = sh.sess.drive.view[]
        v === :image || (sh.shown["_viewBeforeImage"] = v)
        xray_image!(sh.sess; direction = Tuple(f), up = Tuple(u), center = Tuple(c.lookat), width)
        current || (sh.sess.status[] = "computing the drive for the x-ray image…")
        current && sh.sess.drive.image[] !== nothing && (sh.sess.status[] = "x-ray image along the view")
    end
end
br_cavity(on) = _with_shell(sh -> set_cavity!(sh.sess.drive, Bool(on)))
function br_cavity_point(text)
    _with_shell() do sh
        xs = split(strip(string(text)), r"[,\s]+")
        length(xs) == 3 || throw(ArgumentError("give the cavity point as x, y, z (mm)"))
        set_cavity!(sh.sess.drive, sh.sess.drive.cavity[];
                    point = [Controllers._parse_number(x, "a coordinate") for x in xs])
        delete!(sh.shown, "cavityPoint")
    end
end
# The selected parts as thin films (left out of the drive's cavity) or not.
function br_part_film(on)
    _with_shell() do sh
        dr = sh.sess.drive
        f = copy(dr.films[])
        for i in _selected_parts(sh)
            n = sh.sess.target.parts[][i].name
            Bool(on) ? push!(f, n) : delete!(f, n)
        end
        set_cavity!(dr, dr.cavity[]; films = f)
    end
end
br_beam_model(model) = _with_shell(sh -> set_beam_model!(sh.sess.drive, _sym(model)))
# The 1ω (`h` = 1) or 2ω (2) power as a fraction of the 3ω power, as text; empty when off.
function br_unconverted_ratio(h, text)
    _with_shell() do sh
        dr = sh.sess.drive
        x = isempty(strip(string(text))) ? 0.0 :
            Controllers._parse_number(string(text), "the power ratio"; nonnegative = true)
        r = collect(dr.unconverted_ratio[])
        r[round(Int, h)] = x
        set_unconverted_ratio!(dr, r...)
    end
end
# The 1ω (`h` = 1) or 2ω (2) power per beam of the prepulse in TW, empty when off, or its
# start (3) or stop (4) time in ns, as text.
function br_prepulse(h, text)
    _with_shell() do sh
        tm = sh.sess.time
        k = round(Int, h)
        s = strip(string(text))
        p = collect(tm.prepulse[])
        p[k] = k <= 2 && isempty(s) ? 0.0 :
               Controllers._parse_number(String(s), k <= 2 ? "the prepulse power" : "the prepulse time";
                                         nonnegative = k <= 2)
        set_prepulse!(tm, p[1], p[2]; start = p[3], stop = p[4])
    end
end
br_irr_threshold(text) = _with_shell(sh -> set_irradiance_threshold!(sh.sess.time,
                                                                     Controllers._parse_number(string(text), "the threshold";
                                                                                               positive = true)))
br_set_spot_radius(text) = _with_shell(sh -> (sh.sess.drive.spot_radius[] =
                                                  Controllers._parse_number(string(text), "the spot radius"; positive = true)))
function br_import_visrad_workspace(url)
    _with_shell() do sh
        # read and meshed on a worker (seconds for CAD); the design changes when it is ready
        import_visrad_workspace!(sh.sess, _url_to_path(String(url));
                                 done = () -> (select!(sh.sess.view, "target"); fit_all!(sh.scene)))
    end
end
br_export_visrad_workspace(url) = _with_shell(sh -> export_visrad_workspace(sh.sess, _url_to_path(String(url))))
br_export_results(url) = _with_shell(sh -> export_results(sh.sess, _url_to_path(String(url))))
br_export_stl(url) = _with_shell(sh -> export_target_stl(sh.sess, _url_to_path(String(url))))
br_export_script(url) = _with_shell(sh -> write_session_script(_url_to_path(String(url)), sh.sess))
# The retune's errors show in its dialog, where it was asked for.
function br_retune(table, out, by, cap, stay, on, target, level, load)
    _with_shell() do sh
        try
            # what the retune holds (`auto`: the first the design offers), and the level in
            # its units (a capsule's Tᵣ, else TW/mm²)
            qs = Controllers._retune_quantities(sh.sess)
            which = _sym(on)
            k = which === :auto ? (isempty(qs) ? nothing : 1) : findfirst(q -> q.on === which, qs)
            lv = k === nothing ? nothing : Controllers._retune_level(qs[k].kind, String(level))
            # (the reference table's cone and quad columns as the loaded pulses' are read)
            retune_to_reference!(sh.sess, _url_to_path(String(table)), _url_to_path(String(out));
                                 by = _sym(by), unit = Controllers._pulse_unit(sh.sess),
                                 cap = Controllers._parse_number(cap, "the power limit"; positive = true),
                                 stay = Controllers._parse_number(stay, "the weight on the reference powers";
                                                                  nonnegative = true),
                                 target = _sym(target), level = _sym(target) === :uniform ? lv : nothing, on = which,
                                 load = Bool(load))
        catch err
            sh.sess.drive.retune_status[] = "retune failed: " * Controllers._errmsg(err)
        end
    end
end
# Time panel
br_time_chart(key) = _with_shell(sh -> set_chart!(sh.sess.time, _sym(key)), _DIRTY_OBSERVED)
# A click on a beam-table header sorts the display of the rows by that column.
br_sort_beams(key) = _with_shell(sh -> sort_beam_rows!(sh.sess.view, _sym(key)), _DIRTY_OBSERVED)
br_time_pulses(url) = _with_shell(sh -> load_pulses!(sh.sess, _url_to_path(String(url))))
# A pulse table's cone or quad columns as powers per "quad" or per "beam".
br_time_unit(unit) = _with_shell(sh -> set_pulse_unit!(sh.sess, _sym(unit)))
# Fire the beams with a published pulse shape (a `PULSE_SHAPES` key; N221204's without one).
br_time_shape(key = "n221204") = _with_shell(sh -> use_pulse_shape!(sh.sess, _sym(key)))
# What the 3D view shows: "history" (the drive history at the step shown last, or its
# last) or "beams" (the drive at the powers set on the beams).
function br_time_mode(mode)
    _with_shell() do sh
        sess = sh.sess
        if String(mode) == "history"
            h = sess.time.history[]
            h === nothing && return
            set_time_step!(sess, clamp(get(sh.shown, "_lastStep", length(h.time)), 1, length(h.time)))
        else
            set_time_cursor!(sess, NaN)
        end
    end
end
br_time_run() = _with_shell(sh -> start_history!(sh.sess))
br_time_checks() = _with_shell(sh -> start_irradiance!(sh.sess))
br_time_cancel() = _with_shell(sh -> cancel_history!(sh.sess))
br_time_step(k) = _with_shell(sh -> set_time_step!(sh.sess, round(Int, k)), _DIRTY_OBSERVED)
br_time_now() = _with_shell(sh -> set_time_cursor!(sh.sess, NaN), _DIRTY_OBSERVED)
br_time_diagnostic(id, on) = _with_shell(sh -> set_diagnostic!(sh.sess, String(id), Bool(on)))
# Record an instrument of the Diagnostics tab's list (its default configuration), or stop.
br_diag_record(instrument, on) = _with_shell(sh -> Controllers.set_instrument!(sh.sess, String(instrument), Bool(on)))
# Record definition `new` in place of the recorded `id`: a chip's configuration.
br_diag_config(id, new) = _with_shell(sh -> Controllers.set_configuration!(sh.sess, String(id), String(new)))
# The category the Diagnostics tab shows (a diagnostic kind: "xray_power", …).
br_diag_category(kind) = _with_shell(sh -> _set!(sh, "diagCategory", String(kind)))
br_diag_use_path(id) = _with_shell(sh -> use_path_diagnostic!(sh.sess, String(id)))
# The image the Diagnostics tab shows (an imager's id and frame, as `diagImageModel` lists them).
br_diag_pick_image(id, frame) = _with_shell(sh -> (sh.shown["_diagImage"] = (String(id), round(Int, frame))))
# Open an imager's image in the main image view, where it zooms and pans.
function br_diag_image(id, frame)
    _with_shell() do sh
        v = sh.sess.drive.view[]
        v === :image || (sh.shown["_viewBeforeImage"] = v)
        sh.shown["_diagImage"] = (String(id), round(Int, frame))
        show_diagnostic_image!(sh.sess, String(id), round(Int, frame))
    end
end
br_diag_reload() = _with_shell(sh -> refresh_diagnostics!(sh.sess))
br_diag_add_folder(url) = _with_shell(sh -> add_diagnostic_folder!(sh.sess, _url_to_path(String(url))))
br_diag_remove_folder(path) = _with_shell(sh -> remove_diagnostic_folder!(sh.sess, String(path)))
# Optimize
br_opt_set_by(by) = _with_shell(sh -> (by = _sym(by); by in BEAM_GROUPINGS &&
                                                      Controllers._assign!(sh.sess.optimize.by, by)))
br_opt_free(field, on) = _with_shell(sh -> set_free!(sh.sess.optimize, _sym(field), Bool(on)))
br_opt_edit(key, text) = _with_shell(sh -> edit_optimize!(sh.sess.optimize, _sym(key), string(text)))
br_opt_objective(k) = _with_shell(sh -> set_objective!(sh.sess.optimize, _sym(k)))
br_opt_hold(on) = _with_shell(sh -> (sh.sess.optimize.hold_drive[] = Bool(on)))
br_opt_keepout(on) = _with_shell(sh -> (sh.sess.optimize.keepout[] = Bool(on)))
br_opt_capsule(on) = _with_shell(sh -> (sh.sess.optimize.capsule_keepout[] = Bool(on)))
br_opt_start() = br_run_task("optimize")
br_opt_cancel() = _with_shell(sh -> cancel_optimize!(sh.sess))
br_opt_accept() = _with_shell(sh -> accept_optimum!(sh.sess))
# Clearance
br_rv_role(name, role) = _with_shell(sh -> set_role!(sh.sess.review, _sym(name), _sym(role)))
br_rv_edit(key, text) = _with_shell(sh -> edit_review!(sh.sess.review, _sym(key), string(text)))
br_roe_start() = _with_shell(sh -> start_roe!(sh.sess))
br_cm_show(on) = _with_shell(sh -> show_corridor!(sh.sess, Bool(on)))
function br_cm_radius(text)
    _with_shell() do sh
        r = tryparse(Float64, strip(string(text)))
        r === nothing || r <= 0 || !isfinite(r) ? (sh.sess.status[] = "the radius of the keep-out map is a positive length in mm") :
                                                  set_corridor_radius!(sh.sess, r)
    end
end
br_rv_wings(on) = _with_shell(sh -> (sh.sess.review.wings[] = Bool(on)))
br_rv_quality(q) = _with_shell(sh -> (q = _sym(q); haskey(REVIEW_QUALITY, q) &&
                                                   Controllers._assign!(sh.sess.review.quality, q)))
br_rv_start() = br_run_task("clearance")
br_rv_cancel() = _with_shell(sh -> cancel_clearance!(sh.sess))
br_rv_export(url) = _with_shell(sh -> export_clearance(sh.sess.review, _url_to_path(String(url))))
# Sweep
br_sw_subject(i) = _with_shell(sh -> set_sweep!(sh.sess.sweep; subject = sweep_subjects(sh.sess)[Int(i)+1]))
br_sw_frame(i) = _with_shell(sh -> set_sweep!(sh.sess.sweep; frame = _move_frame_choices(sh.sess)[Int(i)+1]))
br_sw_axis(i) = _with_shell(sh -> set_sweep!(sh.sess.sweep; axis = SWEEP_AXES[Int(i)+1]))
br_sw_mode(i) = _with_shell(sh -> set_sweep!(sh.sess.sweep; mode = SWEEP_MODES[Int(i)+1]))
# A `:step` sweep of placement step `k` of the subject part (its first number), or of
# the chosen step's `j`-th number (indices from 0).
function _sweep_step_number!(sess, part::Symbol, k::Integer, j::Integer = 0)
    i = part_index(sess.target, part)
    i === nothing && return (sess.sweep.error[] = "choose a part with placement steps"; false)
    steps = sess.target.parts[][i].steps
    1 <= k <= length(steps) || return false
    nums = step_numbers(steps[k])
    isempty(nums) && return (sess.sweep.error[] = "placement step $k has no number to sweep"; false)
    return set_sweep_step!(sess, part, k, nums[clamp(j + 1, 1, length(nums))][1])
end
br_sw_step(i) = _with_shell(sh -> _sweep_step_number!(sh.sess, sweep_subject(sh.sess), Int(i) + 1))
br_sw_number(j) = _with_shell(sh -> _sweep_step_number!(sh.sess, sweep_subject(sh.sess),
                                                        sh.sess.sweep.step_index[], Int(j)))
# Set up a sweep of the shown part's placement step `k` and show it.
function br_step_sweep(k)
    _with_shell() do sh
        idx = _selected_parts(sh)
        isempty(idx) && return
        _sweep_step_number!(sh.sess, sh.sess.target.parts[][idx[1]].name, Int(k)) || return
        select!(sh.sess.view, "task:sweep")
        request_tab("sweep")
    end
end
br_sw_edit(key, text) = _with_shell(sh -> edit_sweep!(sh.sess.sweep, _sym(key), string(text)))
br_sw_check(key, on) = _with_shell(sh -> set_sweep!(sh.sess.sweep; (_sym(key) => Bool(on),)...))
br_sw_object(i) = _with_shell(sh -> set_sweep!(sh.sess.sweep;
                                              los_object = [:none; surface_names(sh.sess.target.target[])][Int(i)+1]))
br_sw_quality(q) = _with_shell(sh -> set_sweep!(sh.sess.sweep; quality = _sym(q)))
br_sw_start() = br_run_task("sweep")
br_sw_cancel() = _with_shell(sh -> cancel_sweep!(sh.sess))
br_sw_select(k) = _with_shell(sh -> select_sweep_step!(sh.sess.sweep, round(Int, k)))
br_sw_apply(k) = _with_shell(sh -> apply_sweep_step!(sh.sess, round(Int, k)))
# TaLIS
br_tl_role(name, role) = _with_shell(sh -> set_talis_role!(sh.sess.talis, _sym(name), _sym(role)))
br_tl_edit(key, text) = _with_shell(sh -> edit_talis!(sh.sess.talis, _sym(key), string(text)))
br_tl_pivot(p) = _with_shell(sh -> set_talis_pivot!(sh.sess.talis, _sym(p)))
br_tl_omega1(on) = _with_shell(sh -> Controllers._assign!(sh.sess.talis.omega1, Bool(on)))
function br_tl_tight(text)
    _with_shell() do sh
        tl = sh.sess.talis
        try
            set_disposition_filter!(tl; tight = Controllers._parse_number(string(text), "the tight threshold";
                                                                          nonnegative = true))
            tl.error[] = ""
        catch err
            err isa ArgumentError || rethrow()
            tl.error[] = Controllers._errmsg(err)
        end
    end
end
br_tl_hide_curved(on) = _with_shell(sh -> set_disposition_filter!(sh.sess.talis; hide_curved = Bool(on)))
# Put row `i` (from 0) of the TaLIS table under inspection, or take it out when it is
# already there; turn the 3D view to it.
br_tl_focus(i) = _with_shell(sh -> (tl = sh.sess.talis; f = tl.focus[];
                                    f !== nothing && f.index == Int(i) + 1 ? clear_talis_focus!(tl) :
                                    focus_talis_row!(tl, Int(i) + 1)))
br_tl_pinch(on) = _with_shell(sh -> set_talis_pinch!(sh.sess.talis, Bool(on)))
br_tl_view(kind) = _with_shell(sh -> talis_view!(sh.sess, _sym(kind)))
br_tl_start() = br_run_task("talis")
br_tl_cancel() = _with_shell(sh -> cancel_talis!(sh.sess))
br_tl_export(url) = _with_shell(sh -> export_talis(sh.sess.talis, _url_to_path(String(url))))
function br_tl_review(url)
    _with_shell() do sh
        dir = _url_to_path(String(url))
        r = review_1w_command(sh.sess, dir; threshold = sh.sess.talis.tight[], title = "1ω review")
        sh.sess.status[] = string("wrote the 1ω review (", length(r.review.items), " approaches) to ", r.paths[1])
    end
end
br_cancel_all() = _with_shell(sh -> (cancel_optimize!(sh.sess); cancel_clearance!(sh.sess);
                                     cancel_talis!(sh.sess); cancel_sweep!(sh.sess)))
# Viewer
br_set_color(mode) = _with_shell(sh -> set_color!(sh.sess, _sym(mode)), _DIRTY_OBSERVED)
br_beam_style(s) = _with_shell(sh -> (s = _sym(s); s in BEAM_STYLES && Controllers._assign!(sh.sess.view.beams, s)), _DIRTY_OBSERVED)
br_omega1_contour(s) = _with_shell(sh -> set_omega1_contour!(sh.sess.view, _sym(s)), _DIRTY_OBSERVED)
br_beam_color(s) = _with_shell(sh -> (s = _sym(s); s in BEAM_COLORINGS &&
                                                   Controllers._assign!(sh.sess.view.beam_color, s)))
br_beam_length(text) = _with_shell(sh -> Controllers._assign!(sh.sess.view.beam_length,
                                                              isempty(strip(string(text))) ? 0.0 :
                                                              Controllers._parse_number(string(text), "the beam length";
                                                                                        nonnegative = true)))
br_show_off(on) = _with_shell(sh -> Controllers._assign!(sh.sess.view.show_off, Bool(on)), _DIRTY_OBSERVED)
br_set_cutaway(on) = _with_shell(sh -> Controllers._assign!(sh.sess.view.cutaway, Bool(on)), _DIRTY_OBSERVED)
# Cut parameters (the Cutaway button's right-click menu); any change turns the cut on.
br_cut_width(w) = _with_shell(sh -> (set_cut!(sh.sess.view; width = Float64(w), kind = :sector);
                                     Controllers._assign!(sh.sess.view.cutaway, true)), _DIRTY_OBSERVED)
br_cut_azimuth(a) = _with_shell(sh -> (set_cut!(sh.sess.view; azimuth = Float64(a), kind = :sector);
                                       Controllers._assign!(sh.sess.view.cutaway, true)), _DIRTY_OBSERVED)
# Center the cut on the side the 3D view looks from.
function br_cut_toward_view()
    _with_shell() do sh
        c = sh.scene.camera
        d = c.eye - c.lookat
        _, e1, e2 = Controllers._cut_basis(sh.sess.view.cut_axis[])
        x, y = dot(d, e1), dot(d, e2)
        hypot(x, y) > 0 || return
        set_cut!(sh.sess.view; azimuth = atand(y, x), kind = :sector)
        Controllers._assign!(sh.sess.view.cutaway, true)
    end
end
# The cut's kind ("sector" or "plane") and the sector's axis ("x", "y", "z").
br_cut_kind(k) = _with_shell(sh -> (set_cut!(sh.sess.view; kind = _sym(k));
                                    Controllers._assign!(sh.sess.view.cutaway, true)), _DIRTY_OBSERVED)
br_cut_axis(a) = _with_shell(sh -> (set_cut!(sh.sess.view; axis = _sym(a), kind = :sector);
                                    Controllers._assign!(sh.sess.view.cutaway, true)), _DIRTY_OBSERVED)
# The bounds of the sector along its axis, as text; empty when unbounded.
function br_cut_range(lo, hi)
    _with_shell() do sh
        val(x, d, what) = isempty(strip(string(x))) ? d : _parse(x, what)
        set_cut!(sh.sess.view; range = (val(lo, -Inf, "the lower bound"), val(hi, Inf, "the upper bound")))
        Controllers._assign!(sh.sess.view.cutaway, true)
    end
end
# The cut plane: normal along an axis ("x", "y", "z", "-x", …) or toward the view
# ("view": the side facing the eye is removed), through the plane's point at offset 0;
# or the other side kept ("flip": the same plane, its normal turned over).
function br_cut_plane(which)
    _with_shell() do sh
        v = sh.sess.view
        w = String(which)
        if w == "flip"
            flip_cut!(v)
        else
            n = if w == "view"
                c = sh.scene.camera
                Briareus.Vec3((c.eye - c.lookat)...)
            else
                (startswith(w, "-") ? -1 : 1) * Controllers._axis_vector(Symbol(lstrip(w, '-')))
            end
            set_cut_plane!(v; normal = n, offset = 0.0)
        end
        Controllers._assign!(v.cutaway, true)
    end
end
br_cut_offset(x) = _with_shell(sh -> (set_cut_plane!(sh.sess.view; offset = Float64(x));
                                      Controllers._assign!(sh.sess.view.cutaway, true)), _DIRTY_OBSERVED)
# The plane's normal typed as port angles θ, φ (degrees) or as a vector "x, y, z".
br_cut_angles(theta, phi) =
    _with_shell(sh -> (set_cut_plane!(sh.sess.view; angles = (_parse(theta, "the polar angle θ"),
                                                               _parse(phi, "the azimuth φ")));
                       Controllers._assign!(sh.sess.view.cutaway, true)))
br_cut_normal(text) = _with_shell(sh -> (set_cut_plane!(sh.sess.view; normal = _parse_vector(text, "the normal"));
                                         Controllers._assign!(sh.sess.view.cutaway, true)))
# The plane's point, typed (mm), or the target center (TCC).
br_cut_point(x, y, z) =
    _with_shell(sh -> (set_cut_plane!(sh.sess.view; point = (_parse(x, "x"), _parse(y, "y"), _parse(z, "z")));
                       Controllers._assign!(sh.sess.view.cutaway, true)))
function _parse_vector(text, what)
    xs = split(strip(string(text), ['(', ')', '[', ']', ' ']), r"[,\s]+"; keepempty = false)
    length(xs) == 3 || throw(ArgumentError("$what is three numbers, x, y, z"))
    return Tuple(_parse(x, what) for x in xs)
end
# The plane through the selected component's centroid, with its normal if it is flat.
function br_cut_through_selected()
    _with_shell() do sh
        name = selected_component(sh.sess)
        name === nothing && throw(ArgumentError("select a component, window or part to cut through"))
        Controllers.await_target!(sh.sess.target)
        v = sh.sess.view
        cut_through_component!(v, sh.sess.target.target[], name)
        Controllers._assign!(v.cutaway, true)
        sh.sess.status[] = "cut plane through $(component_label(name))" *
                           (v.cut_line[] == (0.0, 0.0, 0.0) ? "" : ", containing the line of sight")
    end
end
# The plane containing the line of sight at port angles θ, φ from the target center or
# (`from_selected`) the selected component's centroid, turned `rotation` degrees about it.
function br_cut_line(theta, phi, from_selected, rotation)
    _with_shell() do sh
        v = sh.sess.view
        origin = (0.0, 0.0, 0.0)
        if Bool(from_selected)
            name = selected_component(sh.sess)
            name === nothing && throw(ArgumentError("select a component, window or part for the line's origin"))
            Controllers.await_target!(sh.sess.target)
            origin = Tuple(component_plane(sh.sess.target.target[], name)[1])
        end
        set_cut_line!(v; angles = (_parse(theta, "the polar angle θ"), _parse(phi, "the azimuth φ")),
                      origin, rotation = Float64(rotation))
        Controllers._assign!(v.cutaway, true)
    end
end
br_cut_rotation(a) = _with_shell(sh -> (set_cut_line!(sh.sess.view; rotation = Float64(a));
                                        Controllers._assign!(sh.sess.view.cutaway, true)), _DIRTY_OBSERVED)
function br_cut_component(name, on)
    _with_shell() do sh
        v = sh.sess.view
        s = copy(v.cut_components[])
        Bool(on) ? push!(s, Symbol(String(name))) : delete!(s, Symbol(String(name)))
        set_cut!(v; components = s)
        Controllers._assign!(v.cutaway, true)
    end
end
# Camera actions bring the 3D view back to the canvas (unless the map sits beside it).
_show_scene!(sh) = sh.sess.drive.side[] || set_view!(sh.sess.drive, :scene)
br_standard_view(name) = _with_shell(sh -> (_show_scene!(sh); standard_view!(sh.scene, _sym(name))))
br_fit_all() = _with_shell(sh -> sh.figure === sh.map.fig ? (reset_limits!(sh.map.ax); sh.map.dirty[] = true) :
                                 fit_all!(sh.scene))
br_fit_selection() = _with_shell(sh -> fit_selection!(sh.scene, sh.sess))
function br_look_from_port(theta, phi)
    _with_shell() do sh
        θ, φ = _parse(theta, "the polar angle θ"), _parse(phi, "the azimuth φ")
        _show_scene!(sh)
        look_from_port!(sh.scene, sh.sess, θ, φ)
        sh.sess.status[] = @sprintf("looking from the port at θ = %g°, φ = %g°", θ, φ)
    end
end
# Port angles of beam `i` of the table (0-based, as the dialog lists them), as text.
function br_port_angles(i)
    sh = _SHELL[]
    sh === nothing && return ["", ""]
    bs = sh.sess.beams.beams[]
    0 <= i < length(bs) || return ["", ""]
    θ, φ = port_angles(bs[i+1])
    return [@sprintf("%.3f", θ), @sprintf("%.3f", mod(φ, 360))]
end
br_look_along_beam() = _with_shell(sh -> (_show_scene!(sh); look_along_beam!(sh.scene, sh.sess)))
br_projection(ortho) = _with_shell(sh -> (set_projection!(sh.scene, Bool(ortho)); _set!(sh, "ortho", Bool(ortho))))
br_turntable(on) = _with_shell(sh -> (set_turntable!(sh.scene, Bool(on)); _set!(sh, "turntable", Bool(on))))
# The canvas lost the keyboard (or the window lost focus): forget held keys and drags.
br_canvas_focus(on) = (sh = _SHELL[]; sh === nothing || Bool(on) || release_input!(sh.scene); nothing)

function _register_qml_functions()
    @qmlfunction br_tick br_grab_started br_select br_select_add br_toggle_expand br_tree_on
    @qmlfunction br_select_range br_node_menu br_node_show br_node_include
    @qmlfunction br_undo br_redo br_open_session br_save_session br_reload_cad_sources br_keep_cad_copies
    @qmlfunction br_relink_part br_part_from_source
    @qmlfunction br_set_cluster br_rename_cluster br_split_quads br_merge_quads br_window_round
    @qmlfunction br_set_facility br_aim_frame br_aim_at_frame
    @qmlfunction br_new_target br_set_resolution br_edit_mesh
    @qmlfunction br_add_window br_edit_window br_window_open
    @qmlfunction br_section_edit br_section_add br_section_remove
    @qmlfunction br_add_part br_add_stl br_edit_part br_part_flag br_part_frame br_part_coords
    @qmlfunction br_move_frame br_move_axis br_move_step br_nudge_part br_sweep_part
    @qmlfunction br_step_add br_step_remove br_step_move br_step_edit br_step_choice
    @qmlfunction br_sw_subject br_sw_frame br_sw_axis br_sw_edit br_sw_check br_sw_object
    @qmlfunction br_sw_quality br_sw_start br_sw_cancel br_sw_select br_sw_apply
    @qmlfunction br_sw_mode br_sw_step br_sw_number br_step_sweep
    @qmlfunction br_part_role br_piece_role br_active_foil br_drive_cavity br_export_stl br_export_script br_cavity br_cavity_point br_part_film
    @qmlfunction br_remove_selected br_help br_docs_url
    @qmlfunction br_new_assembly br_move_into_assembly br_ungroup_assembly br_rename_node
    @qmlfunction br_add_frame br_add_port_frame br_set_triad br_tf_set br_tf_open br_tf_create br_open_transform
    @qmlfunction br_sort_beams
    @qmlfunction br_load_preset br_set_grouping br_edit_beam br_set_beam_on br_edit_node
    @qmlfunction br_import_beams br_export_beams br_import_sst br_import_summary br_export_sst
    @qmlfunction br_export_summary br_run_task br_build br_set_view br_side_view br_xray_image br_close_image br_set_spot_radius br_unconverted_ratio br_beam_model br_export_results
    @qmlfunction br_import_visrad_workspace br_export_visrad_workspace br_retune
    @qmlfunction br_time_chart br_time_pulses br_time_shape br_time_run br_time_cancel br_time_step br_time_mode
    @qmlfunction br_time_unit
    @qmlfunction br_time_checks br_prepulse br_irr_threshold
    @qmlfunction br_time_now br_time_diagnostic br_diag_use_path br_diag_reload br_diag_add_folder
    @qmlfunction br_diag_remove_folder br_diag_image br_diag_record br_diag_config br_diag_category br_diag_pick_image
    @qmlfunction br_opt_set_by br_opt_free br_opt_edit br_opt_objective br_opt_hold
    @qmlfunction br_opt_keepout br_opt_capsule br_opt_start br_opt_cancel br_opt_accept
    @qmlfunction br_rv_role br_rv_edit br_rv_wings br_rv_quality br_rv_start br_rv_cancel br_cm_show br_cm_radius
    @qmlfunction br_rv_export br_tl_role br_tl_edit br_tl_pivot br_tl_omega1 br_tl_start
    @qmlfunction br_tl_cancel br_tl_export br_tl_review br_cancel_all br_tl_tight br_tl_hide_curved br_tl_focus br_tl_view br_tl_pinch
    @qmlfunction br_set_color br_beam_style br_beam_color br_roe_start br_beam_length br_set_cutaway br_show_off br_omega1_contour
    @qmlfunction br_cut_width br_cut_azimuth br_cut_toward_view br_cut_component
    @qmlfunction br_cut_kind br_cut_axis br_cut_range br_cut_plane br_cut_offset
    @qmlfunction br_cut_angles br_cut_normal br_cut_point br_cut_through_selected br_cut_line
    @qmlfunction br_cut_rotation
    @qmlfunction br_set_display br_beam_dot
    @qmlfunction br_standard_view br_fit_all br_fit_selection br_look_along_beam
    @qmlfunction br_look_from_port br_port_angles
    @qmlfunction br_projection br_turntable br_canvas_focus
    return
end
