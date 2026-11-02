# User frames beyond `add_frame!`: the ports a port frame is chosen from, and the
# coordinate transform between any two named frames (after VisRad's coordinate transform
# tool), with a frame made from a transform.

"""
    port_choices(bt::BeamTable) -> Vector{NamedTuple}

The ports of the beam table's beams, one per port, sorted by name: `(; name, theta, phi)`
(degrees). A beam's port is named by its quad (NIF quads), else by the beam. A NIF
quad's angles are the facility table's (`Briareus.port_angles(:nif, quad)`, the NIF
User Guide's port table); any other port's are its beam's (`port_angles`).
"""
function port_choices(bt::BeamTable)
    seen = Dict{String,NTuple{2,Float64}}()
    nif = bt.facility[] === :nif
    for b in bt.beams[]
        name = isempty(b.quad) ? b.name : b.quad
        haskey(seen, name) && continue
        θφ = nif && haskey(_NIF_QUAD_PORTS, name) ? port_angles(:nif, name) : port_angles(b)
        seen[name] = (Float64(θφ[1]), Float64(θφ[2]))
    end
    return [(; name, theta = seen[name][1], phi = seen[name][2]) for name in sort!(collect(keys(seen)))]
end

const _NIF_QUAD_PORTS = Dict(q => p for (p, q) in NIF_QUAD_NAMES)

"""
    transform_coordinates(ts::TargetSetup, from, to, c; input = :cartesian,
                          output = :cartesian, direction = false) -> NTuple{3,Float64}

Coordinates in the design's frame `to` of the point (or, with `direction`, the vector)
with coordinates `c` in frame `from`: `Briareus.transform_point` on the design's frames
([`part_frames`](@ref)); `input` and `output` are one of [`PART_COORDINATES`](@ref).
An unknown frame is an `ArgumentError`.
"""
transform_coordinates(ts::TargetSetup, from::Symbol, to::Symbol, c; kwargs...) =
    transform_point(part_frames(ts), from, to, Tuple(Float64.(c)); kwargs...)

"""
    add_frame_from_transform!(ts::TargetSetup, from, to, c; input = :cartesian,
                              direction = false, name = nothing) -> Bool

Add a user frame ([`add_frame!`](@ref)) made from a transform, placed in frame `to`:
for a point `c` (in `from`, coordinates `input`) its origin is the point and its axes
are those of `from`; for a direction its origin is `to`'s and its `z` axis the
direction (its `x` the direction's θ̂ in `to`). A rejected frame sets `error` and
returns `false`.
"""
function add_frame_from_transform!(ts::TargetSetup, from::Symbol, to::Symbol, c; input::Symbol = :cartesian,
                                   direction::Bool = false, name = nothing)
    try
        fr = part_frames(ts)
        for f in (from, to)
            haskey(fr, f) || throw(ArgumentError("no frame named $f"))
        end
        if direction
            r, θ, φ = transform_point(fr[from], fr[to], Tuple(Float64.(c)); input, output = :spherical, direction = true)
            r > 0 || throw(ArgumentError("a frame's axis needs a nonzero direction"))
            return add_frame!(ts; name, frame = to, theta = θ, phi = φ)
        end
        x, y, z = transform_point(fr[from], fr[to], Tuple(Float64.(c)); input)
        θ, φ, ψ = _euler_angles(fr[to].R' * fr[from].R)
        return add_frame!(ts; name, frame = to, x, y, z, theta = θ, phi = φ, rotation = ψ)
    catch err
        err isa ArgumentError || rethrow()
        ts.error[] = _errmsg(err)
        return false
    end
end
