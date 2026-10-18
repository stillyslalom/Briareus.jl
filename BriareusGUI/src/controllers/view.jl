# What the window shows: the 3D view's coloring, beam drawing, cutaway and hidden
# components, the selection, and the model tree (target components, beams by cone and
# quad, analyses) with its expanded nodes. No GL: the canvas reads this state.

"""
What the 3D view colors the target by: `:component`, or a drive field (see
[`DRIVE_FIELDS`](@ref): absorbed laser, x-ray radiosity, incident x-ray flux).
"""
const COLOR_MODES = (:component, :laser, :radiosity, :incident)

const COLOR_LABELS = Dict(:component => "Components", :laser => "Absorbed laser",
                          :radiosity => "X-ray radiosity", :incident => "Incident x-ray flux")

"Beam colorings: by group (`:group`), or by worst status in the clearance or TaLIS report."
const BEAM_COLORINGS = (:group, :clearance, :talis)

"""
Categorical colors (hex) of components and beam groups, shared by the 3D view and the
model tree's swatches: Wong's seven, then the Tableau 20 without near repeats, so
targets with many parts keep distinct colors (they cycle after 25).
"""
const PALETTE = ("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#F0E442",
                 "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B", "#E377C2",
                 "#7F7F7F", "#BCBD22", "#17BECF", "#AEC7E8", "#FFBB78", "#98DF8A", "#FF9896",
                 "#C5B0D5", "#C49C94", "#DBDB8D", "#9EDAE5")

palette_color(k::Integer) = PALETTE[mod1(k, length(PALETTE))]

"""
    component_colors(t::Target) -> Dict{Symbol,String}

Palette color of each opaque component, in component order (apertures are drawn as
holes and have none).
"""
function component_colors(t::Target)
    names = [c.name for c in t.components if c.kind === :surface]
    return Dict(n => palette_color(k) for (k, n) in enumerate(names))
end

"""
    cluster_colors(bt::BeamTable) -> Dict{Symbol,String}

Palette color of each beam group, in order of first appearance.
"""
cluster_colors(bt::BeamTable) = Dict(c => palette_color(k) for (k, c) in enumerate(clusters(bt)))

"""
    ViewState()

State of the 3D view and the model tree: `color` (one of [`COLOR_MODES`](@ref)),
`beams` (one of [`BEAM_STYLES`](@ref)), `beam_color` (one of
[`BEAM_COLORINGS`](@ref)), `beam_length` (mm drawn before the aim point; 0: twice the
target's size), `cutaway` (cut away part of the target to show its inside: the
components `cut_components` within `cut_width` degrees of azimuth centered on
`cut_azimuth`, by default the quadrant of a hohlraum shell facing −x −y; see
[`set_cut!`](@ref)), `hidden` (component names not drawn), `hidden_beams` (beam names not drawn; whether a
beam takes part in the calculations is the beam table's `enabled`), `selection` (a
tree node id, `""` for none: what Properties shows), `extra` (more nodes of the same
family, beams or components, selected with Ctrl or Shift: edited together with the
selection) and `expanded` (expanded tree node ids).

Tree node ids: `"target"`, `"comp:<name>"`, `"window:<k>"`, `"part:<name>"`,
`"beams"`, `"group:<group>"` (a beam group, see `BeamTable`'s `cluster`),
`"quad:<quad>@<group>"` (a quad's beams in a group), `"beam:<name>"`, `"analyses"`,
`"task:drive"`, `"task:optimize"`, `"task:clearance"`, `"task:talis"`.
"""
struct ViewState
    color::Observable{Symbol}
    beams::Observable{Symbol}
    beam_color::Observable{Symbol}
    beam_length::Observable{Float64}
    cutaway::Observable{Bool}
    hidden::Observable{Set{Symbol}}
    selection::Observable{String}
    extra::Observable{Vector{String}}
    expanded::Observable{Set{String}}
    cut_azimuth::Observable{Float64}
    cut_width::Observable{Float64}
    cut_components::Observable{Set{Symbol}}
    hidden_beams::Observable{Set{String}}
end

"Components the cutaway cuts by default: a hohlraum's shell."
const DEFAULT_CUT_COMPONENTS = (:wall, :cap_top, :cap_bottom)

ViewState() = ViewState(Observable(:component), Observable(:rays), Observable(:group),
                        Observable(0.0), Observable(true), Observable(Set{Symbol}()),
                        Observable(""), Observable(String[]),
                        Observable(Set(["target", "beams", "analyses"])), Observable(225.0),
                        Observable(90.0), Observable(Set{Symbol}(DEFAULT_CUT_COMPONENTS)),
                        Observable(Set{String}()))

"""
    set_cut!(v::ViewState; azimuth, width, components) -> v

Where the cutaway cuts: the sector of `width` degrees of azimuth (90 a quadrant, 180 a
half) centered on `azimuth` (degrees from +x; 225 faces −x −y), through the named
`components` (whole elements, so the cut edge follows the mesh). Unspecified settings
are kept.
"""
function set_cut!(v::ViewState; azimuth = v.cut_azimuth[], width = v.cut_width[],
                  components = v.cut_components[])
    0 < width <= 360 || throw(ArgumentError("the cut width must be in (0, 360] degrees"))
    _assign!(v.cut_azimuth, mod(Float64(azimuth), 360))
    _assign!(v.cut_width, Float64(width))
    _assign!(v.cut_components, Set{Symbol}(components))
    return v
end

"""
    is_cut(v::ViewState, name::Symbol, point) -> Bool

Whether the cutaway removes component `name` at `point` (an element's centroid).
"""
function is_cut(v::ViewState, name::Symbol, p)
    v.cutaway[] && name in v.cut_components[] || return false
    hypot(p[1], p[2]) > 0 || return false
    δ = mod(atand(p[2], p[1]) - v.cut_azimuth[] + 180, 360) - 180
    return abs(δ) < v.cut_width[] / 2
end

"""
    select!(v::ViewState, id) / select_add!(v, id) / toggle_expanded!(v, id) /
    set_hidden!(v, name, hidden)

`select!` selects one node; `select_add!` adds a node to the selection, or removes it
if already selected (Ctrl + click), when it is of the selection's family (beam nodes,
or component nodes: components, parts and windows); other nodes are selected alone.
`select_range!` (Shift + click) selects the nodes of the family from the selection to
`id` in the tree's order.
"""
function select!(v::ViewState, id::AbstractString)
    isempty(v.extra[]) || (v.extra[] = String[])
    _assign!(v.selection, String(id))
    return v
end

_is_beam_node(id) = first(selection_kind(id)) in (:group, :quad, :beam, :beams)

# Nodes that can be selected together: beam nodes, or component nodes.
function _family(id)
    k = first(selection_kind(id))
    k in (:group, :quad, :beam, :beams) && return :beams
    k in (:comp, :part, :window) && return :components
    return nothing
end

function select_add!(v::ViewState, id::AbstractString)
    id = String(id)
    cur = v.selection[]
    f = _family(id)
    (f !== nothing && f === _family(cur)) || return select!(v, id)
    ex = copy(v.extra[])
    if id == cur
        isempty(ex) && return select!(v, "")
        nxt = popfirst!(ex)
        v.extra[] = ex
        _assign!(v.selection, nxt)
    elseif id in ex
        v.extra[] = filter(!=(id), ex)
    else
        push!(ex, id)
        v.extra[] = ex
    end
    return v
end

function select_range!(v::ViewState, order::AbstractVector{<:AbstractString}, id::AbstractString)
    id, cur = String(id), v.selection[]
    f = _family(id)
    i, j = findfirst(==(cur), order), findfirst(==(id), order)
    (f !== nothing && f === _family(cur) && i !== nothing && j !== nothing) || return select!(v, id)
    span = [String(x) for x in order[min(i, j):max(i, j)] if _family(x) === f && x != cur]
    v.extra[] = span
    return v
end

"Every selected node id: the selection and the extra ones."
selected_ids(v::ViewState) = isempty(v.selection[]) ? copy(v.extra[]) : [v.selection[]; v.extra[]]

function toggle_expanded!(v::ViewState, id::AbstractString)
    e = copy(v.expanded[])
    id in e ? delete!(e, id) : push!(e, String(id))
    v.expanded[] = e
    return v
end

function set_beams_hidden!(v::ViewState, names, hidden::Bool)
    h = copy(v.hidden_beams[])
    hidden ? union!(h, String.(names)) : setdiff!(h, String.(names))
    h == v.hidden_beams[] || (v.hidden_beams[] = h)
    return v
end

function set_hidden!(v::ViewState, name::Symbol, hidden::Bool)
    h = copy(v.hidden[])
    hidden ? push!(h, name) : delete!(h, name)
    h == v.hidden[] || (v.hidden[] = h)
    return v
end

"""
    selection_kind(id) -> (kind::Symbol, key::String)

Split a tree node id: `("comp", "wall")` → `(:comp, "wall")`; `"target"` → `(:target, "")`.
"""
function selection_kind(id::AbstractString)
    isempty(id) && return (:none, "")
    i = findfirst(':', id)
    i === nothing && return (Symbol(id), "")
    return (Symbol(id[1:i-1]), String(id[i+1:end]))
end

"Tree node id of a target component (window and part components get their own kinds)."
function component_id(ts, name::Symbol)
    m = match(r"^window_(\d+)$", String(name))
    m !== nothing && ts.kind[] === :hohlraum && parse(Int, m[1]) <= length(ts.windows[]) &&
        return "window:" * m[1]
    any(p -> p.name === name, ts.parts[]) && return "part:" * String(name)
    return "comp:" * String(name)
end

"Tree node id of a beam."
beam_id(b::Beam) = "beam:" * b.name

"""
    TreeRow

One visible row of the model tree: `id`, `depth`, `label`, `detail`, `kind` (the id's
kind), `expandable`, `expanded`, `on` (−1: no switch, 0 off, 1 on, 2 some on: a
beam's or a beam set's enabled flag, a component's visibility), `status` (`""`, or
`"ok"`, `"todo"`, `"attention"`, `"busy"` for analyses), `selected`, `color` (the
swatch: a component's or beam group's color in the 3D view, `""` for none).
"""
const TreeRow = NamedTuple{(:id, :depth, :label, :detail, :kind, :expandable, :expanded, :on,
                            :status, :selected, :color),
                           Tuple{String,Int,String,String,String,Bool,Bool,Int,String,Bool,String}}

_tri(on::Int, n::Int) = on == n ? 1 : on == 0 ? 0 : 2
_row(id, depth, label, detail, expandable, expanded, on, status, sel; color = "") =
    TreeRow((id, depth, label, detail, String(first(selection_kind(id))), expandable,
             expandable && expanded, on, status, id in sel, color))

const _COMPONENT_LABELS = Dict(:wall => "Wall", :cap_top => "Top cap", :cap_bottom => "Bottom cap",
                               :leh_top => "Top LEH", :leh_bottom => "Bottom LEH",
                               :capsule => "Capsule", :backlighter => "Foil", :main => "Main target")

"""
    model_tree(sess) -> Vector{TreeRow}

The visible rows of the model tree, depth first.
"""
function model_tree(sess)
    v = sess.view
    ex, sel = v.expanded[], Set(selected_ids(v))
    rows = TreeRow[]
    ts = sess.target
    t = ts.target[]
    push!(rows, _row("target", 0, "Target: " * TARGET_KIND_LABELS[ts.kind[]],
                     @sprintf("%d elements", nelements(t)), true, "target" in ex, -1, "", sel))
    if "target" in ex
        hid = v.hidden[]
        cols = component_colors(t)
        for c in t.components
            id = component_id(ts, c.name)
            lab = get(_COMPONENT_LABELS, c.name, String(c.name))
            kindtxt = c.kind === :aperture ? "aperture" : c.two_sided ? "two-sided" : ""
            push!(rows, _row(id, 1, lab, kindtxt, false, false, c.name in hid ? 0 : 1, "", sel;
                             color = get(cols, c.name, "")))
        end
        # left out of the calculations: listed (and drawn as a ghost) to be put back
        if ts.ghost[] !== nothing
            for c in ts.ghost[].components
                lab = get(_COMPONENT_LABELS, c.name, String(c.name))
                push!(rows, _row(component_id(ts, c.name), 1, lab, "not in calculations", false,
                                 false, c.name in hid ? 0 : 1, "", sel))
            end
        end
    end
    bt = sess.beams
    bs, en = bt.beams[], bt.enabled[]
    push!(rows, _row("beams", 0, "Beams: " * get(PRESET_LABELS, bt.preset[], String(bt.preset[])),
                     @sprintf("%d of %d on", count(en), length(bs)), !isempty(bs), "beams" in ex,
                     isempty(bs) ? -1 : _tri(count(en), length(bs)), "", sel))
    if "beams" in ex && !isempty(bs)
        cl = bt.cluster[]
        ccol = cluster_colors(bt)
        hasquads = all(b -> !isempty(b.quad), bs)
        for g in clusters(bt)
            gi = findall(==(g), cl)
            gid = "group:" * String(g)
            P = sum(bs[i].power for i in gi if en[i]; init = 0.0)
            push!(rows, _row(gid, 1, String(g), @sprintf("%d beams · %.4g TW", length(gi), P), true,
                             gid in ex, _tri(count(en[gi]), length(gi)), "", sel; color = ccol[g]))
            gid in ex || continue
            if hasquads
                for q in unique(bs[i].quad for i in gi)
                    qi = [i for i in gi if bs[i].quad == q]
                    single = length(qi) == 1
                    qid = single ? beam_id(bs[qi[1]]) : "quad:" * q * "@" * String(g)
                    θ, φ = port_angles(bs[qi[1]])
                    push!(rows, _row(qid, 2, single ? bs[qi[1]].name : q,
                                     single ? @sprintf("θ %.1f° φ %.1f°", θ, φ) :
                                              @sprintf("%d beams · θ %.1f° φ %.1f°", length(qi), θ, φ),
                                     !single, qid in ex, _tri(count(en[qi]), length(qi)), "", sel))
                    (!single && qid in ex) || continue
                    for i in qi
                        push!(rows, _row(beam_id(bs[i]), 3, bs[i].name, @sprintf("%.4g TW", bs[i].power),
                                         false, false, en[i] ? 1 : 0, "", sel))
                    end
                end
            else
                for i in gi
                    θ, φ = port_angles(bs[i])
                    push!(rows, _row(beam_id(bs[i]), 2, bs[i].name, @sprintf("θ %.1f° φ %.1f°", θ, φ),
                                     false, false, en[i] ? 1 : 0, "", sel))
                end
            end
        end
    end
    push!(rows, _row("analyses", 0, "Analyses", "", true, "analyses" in ex, -1, "", sel))
    if "analyses" in ex
        for task in session_tasks(sess)
            st, txt = task_status(sess, task)
            push!(rows, _row("task:" * String(task), 1, TASK_LABELS[task], txt, false, false, -1,
                             String(st), sel))
        end
    end
    return rows
end

"""
    node_beams(bt::BeamTable, id) -> Vector{Int}

Table indices of the beams under a tree node: a beam, a quad in a group (`quad:Q@g`,
or every beam of the quad without `@g`), a group, or every beam (`"beams"`); empty for
other nodes.
"""
function node_beams(bt::BeamTable, id::AbstractString)
    kind, key = selection_kind(id)
    bs, cl = bt.beams[], bt.cluster[]
    kind === :beam && return findall(b -> b.name == key, bs)
    if kind === :quad
        q, g = occursin('@', key) ? split(key, '@'; limit = 2) : (key, nothing)
        return [i for (i, b) in enumerate(bs) if b.quad == q && (g === nothing || String(cl[i]) == g)]
    end
    kind === :group && return findall(c -> String(c) == key, cl)
    kind === :beams && return collect(eachindex(bs))
    return Int[]
end

"""
    selected_beams(sess) -> Vector{Int}

Table indices of the beams under every selected tree node (beams, quads, groups, or
every beam), in table order; empty for other selections.
"""
selected_beams(sess) =
    sort!(unique!(reduce(vcat, (node_beams(sess.beams, id) for id in selected_ids(sess.view));
                         init = Int[])))

"Names of the selected components (components, parts and windows)."
function selected_components(sess)
    out = Symbol[]
    for id in selected_ids(sess.view)
        kind, key = selection_kind(id)
        kind in (:comp, :part) && push!(out, Symbol(key))
        kind === :window && push!(out, window_name(parse(Int, key)))
    end
    return out
end

"""
    selected_component(sess) -> Union{Nothing,Symbol}

The target component of the selected tree node, if it is one.
"""
function selected_component(sess)
    kind, key = selection_kind(sess.view.selection[])
    kind in (:comp, :part) && return Symbol(key)
    kind === :window && return window_name(parse(Int, key))
    return nothing
end
