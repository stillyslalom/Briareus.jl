# What the window shows: the 3D view's coloring, beam drawing (which beams' 1ω, 2ω and
# 3ω light is drawn), each component's display mode, the cutaway, the selection, and
# the model tree (target components, beams by cone and quad, analyses) with its
# expanded nodes. No GL: the canvas reads this state.

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
Display modes of a target component in the 3D view: `:opaque`, `:transparent`
(semi-transparent), `:wireframe` (its mesh edges) and `:invisible`. A component's mode
changes the drawing only: whether it takes part in the calculations is separate
(`include_components!`).
"""
const DISPLAY_MODES = (:opaque, :transparent, :wireframe, :invisible)

const DISPLAY_LABELS = Dict(:opaque => "Opaque", :transparent => "Semi-transparent",
                            :wireframe => "Wireframe", :invisible => "Invisible")

"""
Kinds of cutaway: `:sector` (a pie slice about a coordinate axis, optionally bounded
along it) and `:plane` (everything on the side a cut plane's normal points to).
"""
const CUT_KINDS = (:sector, :plane)

"""
    ViewState()

State of the 3D view and the model tree: `color` (one of [`COLOR_MODES`](@ref)),
`beams` (one of [`BEAM_STYLES`](@ref)), `beam_color` (one of
[`BEAM_COLORINGS`](@ref)), `beam_length` (mm drawn before the aim point; 0: twice the
target's size), `display` (component name => display mode, one of
[`DISPLAY_MODES`](@ref); a component not listed is opaque; see [`set_display!`](@ref)),
`hidden_beams` (beams whose 3ω light is not drawn), `shown_1w` and `shown_2w` (beams
whose unconverted 1ω or 2ω light is drawn, none by default; see
[`set_beams_drawn!`](@ref)), `show_off` (draw the beams switched off, faint; an
imported VisRad workspace whose beams are mostly off turns it off), `selection` (a
tree node id, `""` for none: what Properties shows), `extra` (more nodes of the same
family, beams or components, selected with Ctrl or Shift: edited together with the
selection) and `expanded` (expanded tree node ids). Drawing never changes a
calculation: whether a beam's light takes part is the beam table's
([`set_participation!`](@ref)).

The cutaway (`cutaway` on or off) cuts the components `cut_components` and is of kind
`cut_kind` (one of [`CUT_KINDS`](@ref)): a sector of `cut_width` degrees of azimuth
centered on `cut_azimuth` about the coordinate axis `cut_axis` (`:z` by default; by
default the quadrant of a hohlraum shell facing −x −y), within `cut_range` (lo, hi)
along that axis, or the half-space in front of the plane with normal `cut_normal`
through the point `cut_shift` (mm) along it from `cut_point`. A plane set to contain a
line of sight also keeps the line's direction `cut_line` (zero otherwise) and the
plane's rotation `cut_rotation` (degrees) about it. See [`set_cut!`](@ref),
[`set_cut_plane!`](@ref) and [`set_cut_line!`](@ref).

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
    display::Observable{Dict{Symbol,Symbol}}
    selection::Observable{String}
    extra::Observable{Vector{String}}
    expanded::Observable{Set{String}}
    cut_azimuth::Observable{Float64}
    cut_width::Observable{Float64}
    cut_components::Observable{Set{Symbol}}
    hidden_beams::Observable{Set{String}}
    show_off::Observable{Bool}
    shown_1w::Observable{Set{String}}
    shown_2w::Observable{Set{String}}
    cut_kind::Observable{Symbol}
    cut_axis::Observable{Symbol}
    cut_range::Observable{NTuple{2,Float64}}
    cut_point::Observable{NTuple{3,Float64}}
    cut_normal::Observable{NTuple{3,Float64}}
    cut_shift::Observable{Float64}
    cut_line::Observable{NTuple{3,Float64}}
    cut_rotation::Observable{Float64}
end

"Components the cutaway cuts by default: a hohlraum's shell."
const DEFAULT_CUT_COMPONENTS = (:wall, :cap_top, :cap_bottom)

ViewState() = ViewState(Observable(:component), Observable(:rays), Observable(:group),
                        Observable(0.0), Observable(true), Observable(Dict{Symbol,Symbol}()),
                        Observable(""), Observable(String[]),
                        Observable(Set(["target", "beams", "analyses"])), Observable(225.0),
                        Observable(90.0), Observable(Set{Symbol}(DEFAULT_CUT_COMPONENTS)),
                        Observable(Set{String}()), Observable(true), Observable(Set{String}()),
                        Observable(Set{String}()), Observable(:sector), Observable(:z),
                        Observable((-Inf, Inf)), Observable((0.0, 0.0, 0.0)),
                        Observable((0.0, -1.0, 0.0)), Observable(0.0),
                        Observable((0.0, 0.0, 0.0)), Observable(0.0))

"""
    set_cut!(v::ViewState; kind, azimuth, width, axis, range, components) -> v

Where the cutaway cuts: a sector (`kind = :sector`) of `width` degrees of azimuth (90 a
quadrant, 180 a half) centered on `azimuth` about the coordinate `axis` (`:x`, `:y` or
`:z`; azimuth from +y toward +z about x, from +z toward +x about y, from +x toward +y
about z, so 225 about z faces −x −y), bounded along the axis by `range` = (lo, hi) mm
(`(-Inf, Inf)`: unbounded); or a plane (`kind = :plane`, see [`set_cut_plane!`](@ref)).
The cut goes through the named `components`, clipped exactly along the cut planes
([`cut_pieces`](@ref)). Unspecified settings are kept.
"""
function set_cut!(v::ViewState; kind = v.cut_kind[], azimuth = v.cut_azimuth[], width = v.cut_width[],
                  axis = v.cut_axis[], range = v.cut_range[], components = v.cut_components[])
    kind in CUT_KINDS || throw(ArgumentError("the cut kind must be one of $CUT_KINDS"))
    0 < width <= 360 || throw(ArgumentError("the cut width must be in (0, 360] degrees"))
    axis in (:x, :y, :z) || throw(ArgumentError("the cut axis must be :x, :y or :z"))
    lo, hi = Float64.(range)
    lo < hi || throw(ArgumentError("the cut range must have lo < hi"))
    _assign!(v.cut_kind, Symbol(kind))
    _assign!(v.cut_azimuth, mod(Float64(azimuth), 360))
    _assign!(v.cut_width, Float64(width))
    _assign!(v.cut_axis, Symbol(axis))
    _assign!(v.cut_range, (lo, hi))
    _assign!(v.cut_components, Set{Symbol}(components))
    return v
end

"""
    set_cut_plane!(v::ViewState; normal, angles, point, offset) -> v

Cut with a plane (`v.cut_kind` becomes `:plane`): the cut removes what lies on the side
the normal points to, of the plane through `point` (mm) moved `offset` (mm) along the
normal, as a section slider does. The normal is `normal`, a 3-vector or `:x`, `:y`, `:z`
(the positive axis), or `angles` = (θ, φ) in degrees, the polar angle from +z and the
azimuth from +x toward +y, as for a port (`port_vector`). A new normal releases the
plane from a line of sight ([`set_cut_line!`](@ref)). Unspecified settings are kept.
"""
function set_cut_plane!(v::ViewState; normal = nothing, angles = nothing, point = v.cut_point[],
                        offset = v.cut_shift[])
    normal === nothing || angles === nothing ||
        throw(ArgumentError("the cut plane's normal is a vector or angles, not both"))
    n = angles !== nothing ? port_vector(Float64.(Tuple(angles))...) :
        normal === nothing ? Vec3(v.cut_normal[]...) :
        normal isa Symbol ? _axis_vector(normal) : Vec3(Float64.(Tuple(normal))...)
    all(isfinite, n) && norm(n) > 0 || throw(ArgumentError("the cut plane needs a nonzero normal"))
    p = Vec3(Float64.(Tuple(point))...)
    all(isfinite, p) && isfinite(offset) ||
        throw(ArgumentError("the cut plane's point and offset must be finite"))
    (normal !== nothing || angles !== nothing) && _assign!(v.cut_line, (0.0, 0.0, 0.0))
    _assign!(v.cut_normal, Tuple(normalize(n)))
    _assign!(v.cut_point, Tuple(p))
    _assign!(v.cut_shift, Float64(offset))
    _assign!(v.cut_kind, :plane)
    return v
end

function _axis_vector(a::Symbol)
    a === :x && return Vec3(1, 0, 0)
    a === :y && return Vec3(0, 1, 0)
    a === :z && return Vec3(0, 0, 1)
    throw(ArgumentError("an axis is :x, :y or :z"))
end

"The cut plane's offset (mm) along its normal from its point `cut_point`."
cut_offset(v::ViewState) = v.cut_shift[]

"""
    cut_angles(v::ViewState) -> (θ, φ)

The cut plane's normal as port angles in degrees: the polar angle θ from +z and the
azimuth φ from +x toward +y, in [0, 360).
"""
function cut_angles(v::ViewState)
    n = v.cut_normal[]
    return acosd(clamp(n[3], -1.0, 1.0)), mod(atand(n[2], n[1]), 360)
end

"""
    set_cut_line!(v::ViewState; direction, angles, origin, rotation) -> v

Cut with a plane that contains a line of sight: the line through `origin` (mm) along
`direction` (a 3-vector) or along `angles` = (θ, φ) in degrees (as for a port), with the
plane's normal perpendicular to the line and turned `rotation` degrees about it. At
rotation 0 the normal is ẑ × d / |ẑ × d| for the unit direction d (ŷ × d for a line
along z), so the plane is vertical and contains the line and the z direction; at 90 it
is d × that normal. The plane's offset becomes 0 (it passes through the line). A
rotation alone turns the plane of the current line. Unspecified settings are kept.
"""
function set_cut_line!(v::ViewState; direction = nothing, angles = nothing, origin = v.cut_point[],
                       rotation = v.cut_rotation[])
    direction === nothing || angles === nothing ||
        throw(ArgumentError("the line of sight is a vector or angles, not both"))
    d = angles !== nothing ? port_vector(Float64.(Tuple(angles))...) :
        direction !== nothing ? Vec3(Float64.(Tuple(direction))...) : Vec3(v.cut_line[]...)
    all(isfinite, d) && norm(d) > 0 ||
        throw(ArgumentError("the cut plane needs a line of sight (a nonzero direction)"))
    d = normalize(d)
    n = _line_normal(d, Float64(rotation))
    set_cut_plane!(v; normal = n, point = origin, offset = 0.0)
    _assign!(v.cut_line, Tuple(d))
    _assign!(v.cut_rotation, mod(Float64(rotation), 360))
    return v
end

# The unit normal perpendicular to the unit direction `d`, turned `α` degrees about it
# from the horizontal reference ẑ × d (ŷ × d for d along z).
function _line_normal(d::Vec3, α::Float64)
    e1 = cross(Vec3(0, 0, 1), d)
    norm(e1) < 1e-9 && (e1 = cross(Vec3(0, 1, 0), d))
    e1 = normalize(e1)
    return cosd(α) * e1 + sind(α) * cross(d, e1)
end

"""
    flip_cut!(v::ViewState) -> v

Keep the other side of the cut plane: the plane stays where it is and its normal turns
over (a plane on a line of sight turns 180° about the line).
"""
function flip_cut!(v::ViewState)
    line, rot = v.cut_line[], v.cut_rotation[]
    set_cut_plane!(v; normal = .-v.cut_normal[], offset = -v.cut_shift[])
    if any(!=(0), line)
        _assign!(v.cut_line, line)
        _assign!(v.cut_rotation, mod(rot + 180, 360))
    end
    return v
end

"""
    is_cut(v::ViewState, name::Symbol, point) -> Bool

Whether the cutaway removes component `name` at `point`.
"""
is_cut(v::ViewState, name::Symbol, p) = is_cut(cut_region(v), name, p)

"""
    display_mode(v::ViewState, name::Symbol) -> Symbol

How component `name` is drawn (one of [`DISPLAY_MODES`](@ref)).
"""
display_mode(v::ViewState, name::Symbol) = get(v.display[], name, :opaque)

"""
    set_display!(v::ViewState, names, mode::Symbol) -> v

Draw the components `names` (a name or a collection) in `mode` (one of
[`DISPLAY_MODES`](@ref)).
"""
function set_display!(v::ViewState, names, mode::Symbol)
    mode in DISPLAY_MODES || throw(ArgumentError("the display mode must be one of $DISPLAY_MODES"))
    d = copy(v.display[])
    for n in (names isa Symbol ? (names,) : names)
        mode === :opaque ? delete!(d, Symbol(n)) : (d[Symbol(n)] = mode)
    end
    d == v.display[] || (v.display[] = d)
    return v
end

"""
    beam_drawn(v::ViewState, name, harmonic) -> Bool

Whether the 3D view draws beam `name`'s light of `harmonic` 1 (unconverted 1ω), 2
(unconverted 2ω) or 3 (the 3ω drive).
"""
function beam_drawn(v::ViewState, name::AbstractString, h::Integer)
    h == 3 && return !(name in v.hidden_beams[])
    h == 1 && return name in v.shown_1w[]
    h == 2 && return name in v.shown_2w[]
    throw(ArgumentError("the harmonic is 1, 2 or 3"))
end

"""
    set_beams_drawn!(v::ViewState, names, harmonic, on::Bool) -> v

Draw (or not) the light of `harmonic` (1, 2 or 3; see [`beam_drawn`](@ref)) of the
beams `names`. A view setting: a beam not drawn still takes part in every calculation
it takes part in.
"""
function set_beams_drawn!(v::ViewState, names, h::Integer, on::Bool)
    obs = h == 3 ? v.hidden_beams : h == 1 ? v.shown_1w : h == 2 ? v.shown_2w :
          throw(ArgumentError("the harmonic is 1, 2 or 3"))
    s = copy(obs[])
    (h == 3 ? !on : on) ? union!(s, String.(names)) : setdiff!(s, String.(names))
    s == obs[] || (obs[] = s)
    return v
end

"""
    select!(v::ViewState, id) / select_add!(v, id) / toggle_expanded!(v, id)

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

"""Tree node id of a target component (window and part components get their own kinds;
an entrance hole's closing lip is its hole's)."""
function component_id(ts, name::Symbol)
    m = match(r"^window_(\d+)$", String(name))
    m !== nothing && ts.kind[] === :hohlraum && parse(Int, m[1]) <= length(ts.windows[]) &&
        return "window:" * m[1]
    any(p -> p.name === name, ts.parts[]) && return "part:" * String(name)
    ts.kind[] === :hohlraum && return "comp:" * String(shown_component(name))
    return "comp:" * String(name)
end

"Tree node id of a beam."
beam_id(b::Beam) = "beam:" * b.name

"""
    TreeRow

One visible row of the model tree: `id`, `depth`, `label`, `detail`, `kind` (the id's
kind), `expandable`, `expanded`, `on` (−1: no switch, 0 off, 1 on, 2 some on: a
beam's or a beam set's 3ω participation, a component's visibility), `status` (`""`, or
`"ok"`, `"todo"`, `"attention"`, `"busy"` for analyses), `selected`, `color` (the
swatch: a component's or beam group's color in the 3D view, `""` for none), `display`
(a component's display mode, one of [`DISPLAY_MODES`](@ref); `""` for other rows) and
`matrix` (a beam row's [`beam_matrix`](@ref); `""` for other rows).
"""
const TreeRow = NamedTuple{(:id, :depth, :label, :detail, :kind, :expandable, :expanded, :on,
                            :status, :selected, :color, :display, :matrix),
                           Tuple{String,Int,String,String,String,Bool,Bool,Int,String,Bool,String,
                                 String,String}}

_tri(on::Int, n::Int) = on == n ? 1 : on == 0 ? 0 : 2
_row(id, depth, label, detail, expandable, expanded, on, status, sel; color = "", display = "",
     matrix = "") =
    TreeRow((id, depth, label, detail, String(first(selection_kind(id))), expandable,
             expandable && expanded, on, status, id in sel, color, display, matrix))

"""
    beam_matrix(sess, idx) -> String

The participation and drawing of the beams `idx` (table indices) as VisRad's 2 × 3
matrix of dots, in six characters: whether their 1ω, 2ω and 3ω light takes part in the
calculations ([`participates`](@ref)), then whether it is drawn
([`beam_drawn`](@ref)). Each is `'1'` (every beam), `'0'` (none) or `'2'` (some); a
1ω or 2ω entry is `'-'` when no beam of the set has an unconverted-light footprint
([`has_unconverted`](@ref)).
"""
function beam_matrix(sess, idx)
    bt, v = sess.beams, sess.view
    bs = bt.beams[]
    n = length(idx)
    c(k) = k == n ? '1' : k == 0 ? '0' : '2'
    foot = any(i -> has_unconverted(bs[i]), idx)
    use = map(h -> h < 3 && !foot ? '-' : c(count(i -> participates(bt, i, h), idx)), 1:3)
    drawn = map(h -> h < 3 && !foot ? '-' : c(count(i -> beam_drawn(v, bs[i].name, h), idx)), 1:3)
    return String([use..., drawn...])
end

const _COMPONENT_LABELS = Dict(:wall => "Wall", :cap_top => "Top cap", :cap_bottom => "Bottom cap",
                               :leh_top => "Top LEH", :leh_bottom => "Bottom LEH",
                               :capsule => "Capsule", :backlighter => "Foil", :main => "Main target")

"""
    component_label(name::Symbol) -> String

A component as the model tree names it ("Top cap" for `:cap_top`, "Top LEH lip" for
`:leh_top_lip`); a name without a label (a CAD part, a window) as it is.
"""
function component_label(name::Symbol)
    haskey(_COMPONENT_LABELS, name) && return _COMPONENT_LABELS[name]
    hole = shown_component(name)
    hole !== name && haskey(_COMPONENT_LABELS, hole) && return _COMPONENT_LABELS[hole] * " lip"
    return String(name)
end

_nbeams(n::Integer) = n == 1 ? "1 beam" : "$n beams"

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
        cols = component_colors(t)
        lips = Set(shown_component(c.name) for c in t.components if is_lip(c))
        crow(id, lab, detail, name; color = "") =
            (m = display_mode(v, name);
             _row(id, 1, lab, detail, false, false, m === :invisible ? 0 : 1, "", sel; color,
                  display = String(m)))
        for c in t.components
            is_lip(c) && continue                   # shown as part of its entrance hole
            kindtxt = c.kind === :aperture ? (c.name in lips ? "aperture · closing" : "aperture") :
                      c.two_sided ? "two-sided" : ""
            push!(rows, crow(component_id(ts, c.name), component_label(c.name), kindtxt, c.name;
                             color = get(cols, c.name, "")))
        end
        # left out of the calculations: listed (and drawn as a ghost) to be put back
        if ts.ghost[] !== nothing
            for c in ts.ghost[].components
                is_lip(c) && continue
                push!(rows, crow(component_id(ts, c.name), component_label(c.name),
                                 "not in calculations", c.name))
            end
        end
    end
    bt = sess.beams
    bs, en = bt.beams[], bt.enabled[]
    brow(id, depth, label, detail, expandable, idx; color = "") =
        _row(id, depth, label, detail, expandable, id in ex, _tri(count(en[idx]), length(idx)), "", sel;
             color, matrix = beam_matrix(sess, idx))
    push!(rows, isempty(bs) ?
                _row("beams", 0, "Beams: " * get(PRESET_LABELS, bt.preset[], String(bt.preset[])),
                     "0 of 0 on", false, false, -1, "", sel) :
                brow("beams", 0, "Beams: " * get(PRESET_LABELS, bt.preset[], String(bt.preset[])),
                     all(en) ? "all on" : @sprintf("%d of %d on", count(en), length(bs)), true,
                     collect(eachindex(bs))))
    if "beams" in ex && !isempty(bs)
        cl = bt.cluster[]
        ccol = cluster_colors(bt)
        hasquads = all(b -> !isempty(b.quad), bs)
        for g in clusters(bt)
            gi = findall(==(g), cl)
            gid = "group:" * String(g)
            P = sum(bs[i].power for i in gi if en[i]; init = 0.0)
            push!(rows, brow(gid, 1, String(g), _nbeams(length(gi)) * @sprintf(" · %.4g TW", P), true, gi;
                             color = ccol[g]))
            gid in ex || continue
            if hasquads
                for q in unique(bs[i].quad for i in gi)
                    qi = [i for i in gi if bs[i].quad == q]
                    single = length(qi) == 1
                    qid = single ? beam_id(bs[qi[1]]) : "quad:" * q * "@" * String(g)
                    θ, φ = port_angles(bs[qi[1]])
                    push!(rows, brow(qid, 2, single ? bs[qi[1]].name : q,
                                     single ? @sprintf("θ %.1f° φ %.1f°", θ, φ) :
                                              _nbeams(length(qi)) * @sprintf(" · θ %.1f° φ %.1f°", θ, φ),
                                     !single, qi))
                    (!single && qid in ex) || continue
                    for i in qi
                        push!(rows, brow(beam_id(bs[i]), 3, bs[i].name, @sprintf("%.4g TW", bs[i].power),
                                         false, [i]))
                    end
                end
            else
                for i in gi
                    θ, φ = port_angles(bs[i])
                    push!(rows, brow(beam_id(bs[i]), 2, bs[i].name, @sprintf("θ %.1f° φ %.1f°", θ, φ),
                                     false, [i]))
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
