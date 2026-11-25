# This file defines the rows the item models of the window hold, and the helpers that bring a
# model to a new list of rows.

# ---------------------------------------------------------------- bridge rows

# One visible row of the model tree (see `Controllers.model_tree`).
mutable struct TreeItem
    id::String
    depth::Int
    label::String
    detail::String
    kind::String
    expandable::Bool
    expanded::Bool
    on::Int                                        # -1 none, 0 off, 1 on, 2 mixed
    status::String
    selected::Bool
    color::String                                  # swatch (hex), "" for none
    drawmode::String                               # a component's display mode, "" for none (not
                                                   # `display`: Qt's DisplayRole name)
    matrix::String                                 # a beam row's 2 × 3 dots, "" for none
end

# An editable field of the inspector: key, label, unit, value as text, and its tooltip
# (`_FIELD_TIPS`, "" for none).
mutable struct FieldRow
    key::String
    label::String
    unit::String
    value::String
    tip::String
end
FieldRow(key, label, unit, value) = FieldRow(key, label, unit, value, get(_FIELD_TIPS, String(key), ""))

# The tooltips of the part, assembly and frame fields, by key.
const _FIELD_TIPS = Dict(
    "x" => "In the frame it is placed in", "y" => "In the frame it is placed in",
    "z" => "In the frame it is placed in",
    "pos_r" => "Distance from the placement frame's z axis (cylindrical) or origin (spherical)",
    "pos_phi" => "Azimuth about the placement frame's z, from its +x",
    "pos_theta" => "Polar angle from the placement frame's +z",
    "theta" => "The polar angle the axis (local z) is tilted to, after the rotation about it",
    "phi" => "The azimuth the axis (local z) is tilted toward",
    "rotation" => "The turn about the axis (local z), applied first",
    "albedo" => "The x-ray albedo: the fraction of the incident x-ray flux re-emitted",
    "conversion" => "The fraction of the absorbed laser power converted to x-rays",
    "h" => "The mesh element size of this part; 0: the target's (Target ▸ Mesh)",
    "scale" => "The mesh file's length unit in mm (25.4 for a file in inches)",
    "leh_radius" => "The radius of the entrance holes (the top's, and the bottom's unless set)",
    "bottom_leh_radius" => "The bottom entrance hole's radius; 0: the LEH radius",
    "end_radius" => "The radius of the frustraum's top end (and of the bottom's unless set)",
    "bottom_end_radius" => "The radius of the frustraum's bottom end; 0: the end radius",
    "leh_closure_time" => "When the lips start to close to x-rays",
    "leh_closure_width" => "How long the lips take to close, from the closure time",
    "leh_open_fraction" => "The radius each entrance hole stays open to, as a fraction of its radius")

# A placement step's heading or field: `kind` head, number, choice or check; `choices`
# joined by "|"; `pick` the chosen index, or 1 for a checked box.
mutable struct StepRow
    step::Int
    kind::String
    key::String
    label::String
    unit::String
    value::String
    choices::String
    pick::Int
end

# One beam-table row; numbers travel as text-field text. `on`: 0 off, 1 on, 2 mixed.
mutable struct BeamRow
    number::Int                                    # "index" is reserved in QML delegates
    key::String
    count::Int
    theta::String
    phi::String
    z::String
    r::String
    phiOffset::String
    defocus::String
    power::String
    on::Int
    frame::String                                  # the aim frame's name, its label the tooltip
    frameTip::String
end

mutable struct ModeRow
    l::Int
    epsilon::String
    legendre::String
    before::String                                 # optimization: ε_l before → after
    after::String
end

mutable struct PowerRow
    name::String
    power::String
    percent::String
end

mutable struct ParamRow
    name::String
    start::String
    optimized::String
    lower::String
    upper::String
    atBound::Bool
end

mutable struct CompRow
    name::String                                   # the component's name, for the callbacks
    label::String                                  # as the model tree names it
    role::String                                   # "target"/"capture", "hardware", "diagnostic", "none"
end

mutable struct ClearRow
    beam::String
    component::String
    status::String
    worst::String
    nominal::String
    power::String
    intensity::String
end

# One diagnostic window; numbers as text-field text.
mutable struct WindowRow
    number::Int
    name::String
    phi::String
    z::String
    width::String
    height::String
    open::Bool
    albedo::String
    conversion::String
    round::Bool
end

# TaLIS rows: the summary of one rule, one check and one note.
mutable struct RuleRow
    label::String
    requirement::String
    pass::Int
    review::Int
    fail::Int
    worst::String                                  # worst margin with unit
    where_::String                                 # its beam / component
    status::String                                 # worst status
end

mutable struct TalisRow
    beam::String
    rule::String
    component::String
    status::String
    margin::String
    worst::String
    nominal::String
    required::String
    note::String
end

mutable struct NoteRow
    kind::String                                   # tolerance, settings, skipped, scope, source
    text::String
end

# One step of the placement sweep; numbers as text.
mutable struct SweepRowItem
    number::Int
    value::String
    clearance::String
    worst::String
    fails::Int
    reviews::Int
    los::String
    blocked::String
    covered::String
    captured::String
    objective::String
end

# One instrument of the Diagnostics tab's list (`diagnostic_instruments`): one row per
# instrument, or one per recorded configuration of an instrument recorded in several (a
# script's, an earlier session's): the list shows the `first` rows, the chips the ticked.
mutable struct DiagRow
    id::String                                     # the recorded definition, else the one a tick records
    name::String                                   # the instrument
    kind::String                                   # the diagnostic kind (the tab's category)
    ticked::Bool                                   # recorded in the next drive history
    usable::Bool
    changed::Bool                                  # saved with the session, differs from the search path's
    tip::String
    config::String                                 # the recorded configuration ("" with one only)
    configIds::String                              # the instrument's configurations, joined by "|":
    configLabels::String                           # their definitions' ids and their labels
    foreign::Bool                                  # of another facility than the session's
    first::Bool                                    # the instrument's first row
end

# One folder of the diagnostic search path; one problem found in its definitions.
mutable struct FolderRow
    path::String
    saved::Bool
end

# One image a drive history holds of a recorded imager.
mutable struct ImageRow
    id::String
    label::String
    frame::Int
end

mutable struct IssueRow
    level::String
    where_::String                                 # file name [key]
    message::String
    file::String
end

# One entry of the time chart's key.
mutable struct LegendRow
    label::String
    color::String
    style::String                                  # solid, dash, dot
end

_same(a, b) = all(f -> getfield(a, f) == getfield(b, f), fieldnames(typeof(a)))

# How many leading rows of `old` and `new` are the same nodes (equal first fields): rows
# are inserted or removed after them.
function _common_lead(old::AbstractVector, new::AbstractVector)
    p = 0
    while p < min(length(old), length(new)) && getfield(old[p+1], 1) == getfield(new[p+1], 1)
        p += 1
    end
    return p
end

# Replace the rows of an item model. Only changed rows are re-sent, and a change in the
# count inserts or removes rows where the old and new rows part, aligned on their first
# field (the id of a tree node), so the view keeps the delegates of the others. A full reset
# rebuilt every delegate of the model tree, a noticeable delay on a CAD workspace each time
# a node was expanded or a part added. `new` may be an empty comprehension whose eltype
# inference gave up on (`Vector{Union{}}`).
function _sync_rows!(rows::Vector{R}, model::JuliaItemModel, new::AbstractVector) where {R}
    new = convert(Vector{R}, new)
    n, m = length(rows), length(new)
    if n != m && n > 0 && m > 0
        p = _common_lead(rows, new)
        if m > n
            QML.begin_insert_rows(model, p + 1, p + m - n)
            for k in 1:m-n
                insert!(rows, p + k, new[p+k])
            end
            QML.@deferredcall QML.end_insert_rows(model)
        else
            QML.remove_rows!(model, p + 1, n - m)
        end
    end
    if length(new) == length(rows)
        for (i, r) in enumerate(new)
            _same(rows[i], r) && continue
            rows[i] = r
            model[i] = r
        end
    else
        empty!(rows)
        append!(rows, new)
        QML.force_model_update(model)
    end
    return
end
