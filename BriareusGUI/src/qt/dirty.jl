# This file defines the dirty domains of the window. An edit marks the panels it can change,
# and the next tick refreshes those panels only (`_refresh!`).

"""
The domains of what the window shows, each refreshed by one function of refresh.jl:
`:tree` (the model tree, the target page and the frames), `:selection` (the inspector),
`:beams` (the beam table), `:drive`, `:time` (the Time and Diagnostics tabs), `:optimize`,
`:review` (the clearance check), `:talis`, `:sweep` and `:view` (the viewer's options and
the list of running jobs). The title, the status line and the undo state belong to no
domain: they are refreshed after every edit.
"""
const DIRTY_DOMAINS = (:tree, :selection, :beams, :drive, :time, :optimize, :review, :talis, :sweep, :view)

const _DIRTY_BITS = NamedTuple{DIRTY_DOMAINS}(ntuple(i -> UInt16(1) << (i - 1), length(DIRTY_DOMAINS)))
const _DIRTY_NONE = UInt16(0)
# An edit was made: the undo history notes it and the title and status are refreshed.
const _DIRTY_EDIT = UInt16(1) << length(DIRTY_DOMAINS)
const _DIRTY_ALL = reduce(|, values(_DIRTY_BITS)) | _DIRTY_EDIT

"""
    _dirty_mask(domains::Symbol...) -> UInt16

The mask of an edit that can change `domains` ([`DIRTY_DOMAINS`](@ref)); with no domain,
an edit that changes none of the panels. Throws `ArgumentError` for an unknown domain.
"""
function _dirty_mask(domains::Symbol...)
    mask = _DIRTY_EDIT
    for d in domains
        haskey(_DIRTY_BITS, d) || throw(ArgumentError("unknown dirty domain :$d (one of $DIRTY_DOMAINS)"))
        mask |= _DIRTY_BITS[d]
    end
    return mask
end

"Whether `mask` holds domain `d`."
_is_dirty(mask::UInt16, d::Symbol) = mask & _DIRTY_BITS[d] != 0

"The domains `mask` holds, in the order of [`DIRTY_DOMAINS`](@ref)."
_dirty_domains(mask::UInt16) = Tuple(d for d in DIRTY_DOMAINS if _is_dirty(mask, d))

# Add `mask` to the dirty state `ref`, and take the state (leaving it clean).
_mark!(ref::Base.RefValue{UInt16}, mask::UInt16 = _DIRTY_ALL) = (ref[] |= mask; nothing)
_take_dirty!(ref::Base.RefValue{UInt16}) = (mask = ref[]; ref[] = _DIRTY_NONE; mask)
