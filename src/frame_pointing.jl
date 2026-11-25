# This file parameterizes pointing in frames. Each aim is given in the local coordinates
# (u, v, w) of a frame per beam, such as the frame of a foil, an entrance hole or a port.
# The non-aim fields are those of `GroupPointing`.

const _FRAME_FIELDS = (:u, :v, :w, :defocus, :power, :spot, :dlambda)
const _FRAME_AIM = (:u, :v, :w)

"""
    FramePointing(beams, frames; by = :group, fields = (:u, :v, :defocus),
                  relative = false)

Maps a parameter vector to beams whose aims are local coordinates of a [`Frame`](@ref).
Beam `i` aims at `to_chamber(frame_i, (u, v, w))`, which is `R * SVector(u, v, w) + origin`.
`frames` is one `Frame` for every beam, a vector of frames (one per beam), or a function
`beam -> Frame`. Beams are grouped by `by` (`:group`, `:quad`, `:beam`, or a function
`beam -> key`), and each group gets one value per field:

- `:u`, `:v`, `:w` — the local coordinates of the aim (mm). On the
  [`foil_frame`](@ref) of a foil, `u` and `v` run across its face and `w` along its
  front normal.
- `:defocus`, `:power`, `:spot`, `:dlambda` — as for [`GroupPointing`](@ref).

A group shares local coordinates. Beams of one group in different frames aim at the
same local point of each frame, so their chamber aims differ. With two entrance-hole
frames, the bottom one the top one mirrored through the equator, a mirrored pointing is
given without a `mirror` flag. For aim coordinates that are not parameters, the value of
each beam is kept, which is its template aim in its frame. With no aim field at all, the
template aims are kept exactly.

`relative = true` makes the aim fields shifts added to the local coordinates
of each beam, so a group moves rigidly in its frame; its initial aim parameters are then
zero. The other fields stay absolute.

`apply` is generic in the number type (frames are Float64, parameters may be
ForwardDiff duals), so the parameterization can be used with [`Scenario`](@ref),
[`optimize_pointing`](@ref) and [`subset`](@ref) in the same way as a `GroupPointing`.
"""
struct FramePointing{K} <: AbstractParameterization
    template::Vector{Beam{Float64}}
    keys::Vector{K}
    members::Vector{Vector{Int}}
    fields::Vector{Symbol}
    relative::Bool
    frames::Vector{Frame}
    local_aims::Vector{Vec3}          # template aims in their frames
    beam_keys::Vector{K}              # group key of each beam
end

function FramePointing(beams::AbstractVector{<:Beam}, frames; by = :group,
                       fields = (:u, :v, :defocus), relative::Bool = false)
    by === :quad && any(b -> isempty(b.quad), beams) &&
        throw(ArgumentError("by = :quad needs quad names on every beam (OMEGA beams have none)"))
    key = by === :group ? (b -> b.group) : by === :quad ? (b -> b.quad) :
          by === :beam ? (b -> b.name) : by
    fr = frames isa Frame ? fill(frames, length(beams)) :
         frames isa AbstractVector ? collect(Frame, frames) : Frame[frames(b) for b in beams]
    length(fr) == length(beams) ||
        throw(DimensionMismatch("$(length(fr)) frames for $(length(beams)) beams"))
    return _frame_pointing(collect(beams), fr, map(key, beams), fields, relative)
end

function _frame_pointing(beams, frames::Vector{Frame}, ks::AbstractVector, fields, relative)
    for f in fields
        f in _FRAME_FIELDS || throw(ArgumentError("unknown field $f; choose from $_FRAME_FIELDS"))
    end
    allunique(fields) || throw(ArgumentError("fields repeat: $fields"))
    uk = unique(ks)
    members = [findall(==(k), ks) for k in uk]
    tmpl = Beam{Float64}[_value_beam(b) for b in beams]
    loc = Vec3[Vec3(to_frame(f, b.aim)) for (f, b) in zip(frames, tmpl)]
    return FramePointing(tmpl, uk, members, collect(Symbol, fields), relative, frames, loc,
                         collect(ks))
end

"Number of parameters of a `FramePointing` (groups × fields)."
nparams(p::FramePointing) = length(p.keys) * length(p.fields)

"Names `\"<group>.<field>\"` for each parameter of a `FramePointing`."
parameter_names(p::FramePointing) = [string(k, ".", f) for k in p.keys for f in p.fields]

_pidx(p::FramePointing, g, fi) = (g - 1) * length(p.fields) + fi

"""
    initial_parameters(p::FramePointing) -> Vector{Float64}

Parameters that reproduce the template beams. They are the group means of the local aim
coordinates and of the other fields. The aim fields of a `relative` parameterization are
zero shifts.
"""
function initial_parameters(p::FramePointing)
    x = zeros(nparams(p))
    for (g, mem) in enumerate(p.members), (fi, f) in enumerate(p.fields)
        p.relative && f in _FRAME_AIM && continue
        vals = map(mem) do bi
            b = p.template[bi]
            f === :u ? p.local_aims[bi][1] : f === :v ? p.local_aims[bi][2] :
            f === :w ? p.local_aims[bi][3] : f === :defocus ? b.defocus :
            f === :power ? b.power : f === :dlambda ? b.dlambda : 1.0
        end
        x[_pidx(p, g, fi)] = sum(vals) / length(vals)
    end
    return x
end

"""
    parameters(p::FramePointing, values::Pair...; group = (field = value, ...), ...)
        -> Vector{Float64}

The parameter vector of `p` with named values, as for a `GroupPointing`. Fields and
groups not named are set to their [`initial_parameters`](@ref).
"""
function parameters(p::FramePointing, values::Pair...; kwargs...)
    x = initial_parameters(p)
    for (key, fv) in (values..., pairs(kwargs)...)
        g = findfirst(==(key), p.keys)
        g === nothing && throw(ArgumentError("no group $(repr(key)); the groups are " *
                                             join(repr.(p.keys), ", ")))
        for (f, v) in pairs(fv)
            fi = findfirst(==(f), p.fields)
            fi === nothing && throw(ArgumentError("$f is not a parameter; the fields are " *
                                                  join(p.fields, ", ")))
            x[_pidx(p, g, fi)] = v
        end
    end
    return x
end

"""
    apply(p::FramePointing, x) -> Vector{Beam}

Beams for parameter vector `x`. Each aim is `R * SVector(u, v, w) + origin` in the frame
of its beam. The element type is that of `x`, so the map is differentiable.
"""
function apply(p::FramePointing, x::AbstractVector{T}) where {T}
    length(x) == nparams(p) || throw(DimensionMismatch("expected $(nparams(p)) parameters"))
    out = Vector{Beam{T}}(undef, length(p.template))
    anyaim = any(in(_FRAME_AIM), p.fields)
    for (g, mem) in enumerate(p.members)
        val(f) = (fi = findfirst(==(f), p.fields); fi === nothing ? nothing : x[_pidx(p, g, fi)])
        u, v, w = val(:u), val(:v), val(:w)
        dfc, pw, sp, dl = val(:defocus), val(:power), val(:spot), val(:dlambda)
        for bi in mem
            b = p.template[bi]
            aim = anyaim ? _frame_aim(p.frames[bi], p.local_aims[bi], u, v, w, p.relative, T) :
                  SVector{3,T}(b.aim)
            out[bi] = _pointed(b, aim, T, sp, pw, dfc, dl)
        end
    end
    return out
end

# Returns the chamber aim of local coordinates. Each of `u`, `v` and `w` is a value, a
# shift of the template coordinate in `l` (with `relative`), or `nothing` for the
# template coordinate.
@inline function _frame_aim(f::Frame, l::Vec3, u, v, w, relative::Bool, ::Type{T}) where {T}
    c(x, l0) = x === nothing ? T(l0) : relative ? l0 + x : x
    return f.R * SVector{3,T}(c(u, l[1]), c(v, l[2]), c(w, l[3])) + f.origin
end

"""
    bounds(p::FramePointing, x0; u = 1.0, v = 1.0, w = 1.0, defocus = 2.0,
           power = 0.5, spot = 0.5, dlambda = 2.0, power_min = 0.0, power_max = Inf)
           -> (lower, upper)

Box bounds around `x0`: `±u`, `±v`, `±w`, `±defocus` mm, relative `±power` and `±spot`
fractions and `±dlambda` Å (engineering defaults that are not facility limits), with power
clipped to `power_min ≤ P ≤ power_max` (TW per beam). A group at zero power gets
`[0, power · P̄]`, `P̄` the mean nonzero power, as for a `GroupPointing`.
"""
bounds(p::FramePointing, x0 = initial_parameters(p); u = 1.0, v = 1.0, w = 1.0,
       defocus = 2.0, power = 0.5, spot = 0.5, dlambda = 2.0, power_min = 0.0,
       power_max = Inf) =
    _box_bounds(p, x0, (; u, v, w, defocus, dlambda), power, spot, power_min, power_max)

Scenario(t::Target, p::FramePointing; kwargs...) = _scenario(t, p.template, p; kwargs...)
_scenario_beams(p::FramePointing) = p.template

_subset_param(p::FramePointing, keep) =
    _frame_pointing(p.template[keep], p.frames[keep], p.beam_keys[keep], p.fields, p.relative)

_nparams_shown(p::FramePointing) = nparams(p)

# Returns `p` with the frames and template aims of the beams `moved` (a Bool per beam)
# translated by `Δ`. The local aims are unchanged, so the parameters have the same
# meaning; a placement sweep uses this for a warm start.
function _translated(p::FramePointing, Δ::Vec3, moved::AbstractVector{Bool})
    fr = [m ? Frame(f.R, f.origin + Δ) : f for (f, m) in zip(p.frames, moved)]
    tm = [m ? remake(b; aim = b.aim + Δ) : b for (b, m) in zip(p.template, moved)]
    return FramePointing(tm, p.keys, p.members, p.fields, p.relative, fr, p.local_aims,
                         p.beam_keys)
end
