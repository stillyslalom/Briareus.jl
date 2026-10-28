# Pointing parameterized in frames: aims as local coordinates (u, v, w) of a frame per
# beam (a foil's, an entrance hole's, a port's), with the non-aim fields of
# `GroupPointing`.

const _FRAME_FIELDS = (:u, :v, :w, :defocus, :power, :spot, :dlambda)
const _FRAME_AIM = (:u, :v, :w)

"""
    FramePointing(beams, frames; by = :group, fields = (:u, :v, :defocus),
                  relative = false)

Maps a parameter vector to beams whose aims are local coordinates of a [`Frame`](@ref):
beam `i` aims at `to_chamber(frame_i, (u, v, w))`, `R * SVector(u, v, w) + origin`.
`frames` is one `Frame` for every beam, a vector of frames (one per beam), or a function
`beam -> Frame`. Beams are grouped by `by` (`:group`, `:quad`, `:beam`, or a function
`beam -> key`), and each group gets one value per field:

- `:u`, `:v`, `:w` — the aim's local coordinates (mm): on a foil's
  [`foil_frame`](@ref), `u` and `v` across its face and `w` along its front normal;
- `:defocus`, `:power`, `:spot`, `:dlambda` — as for [`GroupPointing`](@ref).

A group shares local coordinates, not chamber aims: beams of one group in different
frames aim at the same local point of each. Two entrance holes' frames, the bottom one
the top one mirrored through the equator, so give a mirrored pointing with no `mirror`
flag. Aim coordinates that are not parameters keep each beam's own (its template aim in
its frame); with no aim field at all the template aims are kept exactly.

`relative = true` makes the aim fields shifts added to each beam's own local
coordinates, so a group moves rigidly in its frame; its initial aim parameters are then
zero. The other fields stay absolute.

`apply` is generic in the number type (frames are Float64, parameters may be
ForwardDiff duals), so the parameterization serves [`Scenario`](@ref),
[`optimize_pointing`](@ref) and [`subset`](@ref) as a `GroupPointing` does.
"""
struct FramePointing{K}
    template::Vector{Beam{Float64}}
    keys::Vector{K}
    members::Vector{Vector{Int}}
    fields::Vector{Symbol}
    relative::Bool
    frames::Vector{Frame}
    local_aims::Vector{Vec3}          # template aims in their frames
    beam_keys::Vector{K}              # each beam's group key
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

Parameters reproducing the template beams: group means of the local aim coordinates
and of the other fields (zero shifts for the aim fields of a `relative`
parameterization).
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

The parameter vector of `p` with named values, as for a `GroupPointing`: fields and
groups not named keep [`initial_parameters`](@ref).
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

Beams for parameter vector `x`: each aim `R * SVector(u, v, w) + origin` in its beam's
frame (element type follows `x`, so this is differentiable).
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

# The chamber aim of local coordinates: each of `u`, `v`, `w` a value, a shift of the
# template's `l` (`relative`), or `nothing` (the template's).
@inline function _frame_aim(f::Frame, l::Vec3, u, v, w, relative::Bool, ::Type{T}) where {T}
    c(x, l0) = x === nothing ? T(l0) : relative ? l0 + x : x
    return f.R * SVector{3,T}(c(u, l[1]), c(v, l[2]), c(w, l[3])) + f.origin
end

"""
    bounds(p::FramePointing, x0; u = 1.0, v = 1.0, w = 1.0, defocus = 2.0,
           power = 0.5, spot = 0.5, dlambda = 2.0, power_min = 0.0, power_max = Inf)
           -> (lower, upper)

Box bounds around `x0`: `±u`, `±v`, `±w`, `±defocus` mm, relative `±power` and `±spot`
fractions and `±dlambda` Å (engineering defaults, not facility limits), with power
clipped to `power_min ≤ P ≤ power_max` (TW per beam).
"""
function bounds(p::FramePointing, x0 = initial_parameters(p); u = 1.0, v = 1.0, w = 1.0,
                defocus = 2.0, power = 0.5, spot = 0.5, dlambda = 2.0, power_min = 0.0,
                power_max = Inf)
    lo, hi = copy(x0), copy(x0)
    for g in eachindex(p.members), (fi, f) in enumerate(p.fields)
        i = _pidx(p, g, fi)
        if f === :power
            lo[i] = max(x0[i] * (1 - power), power_min)
            hi[i] = min(x0[i] * (1 + power), power_max)
            lo[i] <= hi[i] || throw(ArgumentError("empty power bounds for group $(p.keys[g])"))
        elseif f === :spot
            lo[i] = x0[i] * (1 - spot); hi[i] = x0[i] * (1 + spot)
        else
            δ = f === :u ? u : f === :v ? v : f === :w ? w : f === :dlambda ? dlambda : defocus
            lo[i] = x0[i] - δ; hi[i] = x0[i] + δ
        end
    end
    return lo, hi
end

Scenario(t::Target, p::FramePointing; kwargs...) = _scenario(t, p.template, p; kwargs...)

_subset_param(p::FramePointing, keep) =
    _frame_pointing(p.template[keep], p.frames[keep], p.beam_keys[keep], p.fields, p.relative)

_nparams_shown(p::FramePointing) = nparams(p)

# `p` with the frames and template aims of the beams `moved` (a Bool per beam)
# translated by `Δ`; the local aims stay as they are, so the parameters keep their
# meaning (a warm start of a placement sweep).
function _translated(p::FramePointing, Δ::Vec3, moved::AbstractVector{Bool})
    fr = [m ? Frame(f.R, f.origin + Δ) : f for (f, m) in zip(p.frames, moved)]
    tm = [m ? remake(b; aim = b.aim + Δ) : b for (b, m) in zip(p.template, moved)]
    return FramePointing(tm, p.keys, p.members, p.fields, p.relative, fr, p.local_aims,
                         p.beam_keys)
end
