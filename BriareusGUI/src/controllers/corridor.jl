# This file covers the keep-out map of the 3D view and the Clearance tab: a sphere about TCC
# or a plane through the mounting zone, colored by the light and the lines of sight that
# cross it, with the open region where new hardware fits.

"""
    CorridorState

The keep-out map of a session (`Briareus.corridor_map`). `surface` is `:sphere` or
`:plane`. `radius` (mm) is the radius of the sphere about TCC. The plane passes through
`point` (mm) with the normal `normal` and spans `extent` mm on each side of the point.
`shadows` is whether the hardware of the target shadows the keep-outs (`target` of
`Briareus.corridor_map`). `shown` is whether the map is drawn, `map` the map of the last update (a core
`CorridorMap`, or `nothing`), `status` a line about it, `error` the message of a rejected
setting and `inputs` the beams, sight lines, surface and settings of that update.
"""
struct CorridorState
    surface::Observable{Symbol}
    radius::Observable{Float64}
    point::Observable{NTuple{3,Float64}}
    normal::Observable{NTuple{3,Float64}}
    extent::Observable{Float64}
    shadows::Observable{Bool}
    shown::Observable{Bool}
    map::Observable{Union{Nothing,CorridorMap}}
    status::Observable{String}
    error::Observable{String}
    inputs::Base.RefValue{Any}
end

CorridorState() = CorridorState(Observable(:sphere), Observable(50.0), Observable((0.0, 0.0, 0.0)),
                                Observable((0.0, 0.0, 1.0)), Observable(20.0), Observable(true), Observable(false),
                                Observable{Union{Nothing,CorridorMap}}(nothing), Observable(""), Observable(""),
                                Ref{Any}(nothing))

# The surface of the map as `corridor_map` takes it, and the inputs it is compared by.
_corridor_where(cs::CorridorState) =
    cs.surface[] === :plane ? (; plane = (cs.point[], cs.normal[]), extent = cs.extent[]) : (; radius = cs.radius[])

"""
Colors (hex) of the keep-out map by drawing class (`Briareus.corridor_fills`): 3ω light,
1ω light, a line of sight, open with a tight margin, open, and keep-outs that the hardware
of the target shadows.
"""
const CORRIDOR_COLORS = ("#0072B2", "#D55E00", "#CC79A7", "#F0E442", "#009E73", "#999999")

"Opacity of each drawing class of the keep-out map: the open region stays faint."
const CORRIDOR_OPACITY = (0.55, 0.55, 0.55, 0.4, 0.12, 0.3)

"""
    corridor_kwargs(rv::ReviewState, tl::TaLISState) -> NamedTuple

The keyword arguments of the `corridor_map` call of a session: the error stack of the
clearance check, and the rule tolerance, its pivot and the 1ω switch of the TaLIS checks.
"""
function corridor_kwargs(rv::ReviewState, tl::TaLISState)
    e = rv.errors[]
    return (; tolerance = tl.tolerance[], pivot = tl.pivot[] === :point ? tl.pivot_point[] : tl.pivot[],
            pointing = e[:pointing], target_offset = e[:target_offset],
            hardware_tolerance = e[:hardware_tolerance], focus = e[:focus], omega1 = tl.omega1[])
end

"""
    update_corridor!(cs::CorridorState, rv, tl, beams, sight_lines; omitted = String[], target = nothing) -> cs

Compute the keep-out map of `beams` and `sight_lines` on the sphere of `cs.radius` about
TCC, or on the plane of `cs.point`, `cs.normal` and `cs.extent` when `cs.surface` is
`:plane`, with the settings of [`corridor_kwargs`](@ref), and leave it in `cs.map`. With a
`target` and `cs.shadows`, its hardware shadows the keep-outs. The map of 192 beams
takes a fraction of a second, so it is computed in place. The status line ends with the
names in `omitted`, the diagnostics that gave no line of sight.
"""
function update_corridor!(cs::CorridorState, rv::ReviewState, tl::TaLISState, beams::AbstractVector,
                          sight_lines::AbstractVector; omitted::AbstractVector{<:AbstractString} = String[],
                          target::Union{Nothing,Target} = nothing)
    kw = corridor_kwargs(rv, tl)
    bs = collect(beams)
    sl = collect(SightLine, sight_lines)
    where = _corridor_where(cs)
    t = cs.shadows[] ? target : nothing
    m = corridor_map(bs; sight_lines = sl, target = t, where..., kw...)
    cs.inputs[] = (bs, sl, where, kw, t)
    text = corridor_text(m; tight = tl.tight[])
    isempty(omitted) || (text *= "; left out, no line of sight: " * join(omitted, ", "))
    cs.status[] = text
    cs.map[] = m
    return cs
end

"""
    clear_corridor!(cs::CorridorState) -> cs

Drop the keep-out map and its drawing.
"""
function clear_corridor!(cs::CorridorState)
    cs.inputs[] = nothing
    _assign!(cs.status, "")
    cs.map[] === nothing || (cs.map[] = nothing)
    return cs
end

"""
    corridor_stale(cs::CorridorState, rv, tl, beams, sight_lines) -> Bool

Whether anything the keep-out map depends on changed since its last update: the enabled
beams, the sight lines, the surface, a setting or the `target` that shadows the keep-outs.
"""
function corridor_stale(cs::CorridorState, rv::ReviewState, tl::TaLISState, beams::AbstractVector,
                        sight_lines::AbstractVector; target::Union{Nothing,Target} = nothing)
    inp = cs.inputs[]
    inp === nothing && return false
    return inp[1] != beams || inp[2] != sight_lines || inp[3] != _corridor_where(cs) || inp[4] != corridor_kwargs(rv, tl) ||
           inp[5] !== (cs.shadows[] ? target : nothing)
end

"""
    corridor_text(map; tight = 1.0) -> String

One line about a keep-out map: the open share of the sphere or plane, the share with more
than the `tight` margin, the widest margin and where it lies, and the number of beams and
lines of sight.
"""
function corridor_text(m::CorridorMap; tight::Real = 1.0)
    n(class) = count(s -> first(s) === class, m.sources)
    io = IOBuffer()
    vec(v) = join((@sprintf("%g", x) for x in v), ", ")
    m.surface === :sphere ? @printf(io, "sphere of %g mm about TCC", m.radius) :
    @printf(io, "plane through (%s) mm, normal (%s), ±%g mm", vec(m.origin), vec(m.axes[3]), m.u[end])
    @printf(io, ": %.1f%% open, %.1f%% with more than %g mm of margin",
            100 * corridor_open_fraction(m), 100 * corridor_open_fraction(m; margin = tight), tight)
    if any(!iszero, m.shadowed)
        # the share of the surface where the hardware shadows a keep-out, and where that
        # leaves the sample open
        w = m.surface === :sphere ? [sind(θ) for θ in m.u, _ in m.v] : ones(size(m.margin))
        part = sum(w[m.shadowed .!= 0]) / sum(w)
        @printf(io, "; the target shadows keep-outs over %.1f%%, leaving %.1f%% more open", 100 * part,
                100 * (corridor_open_fraction(m; shadowed = true) - corridor_open_fraction(m)))
    end
    best, at = findmax(m.margin)
    if isfinite(best) && best > 0 && m.surface === :sphere
        @printf(io, "; widest margin %.1f mm (%.1f°) at θ %.0f°, φ %.0f°", best, angular_margin(m)[at], m.u[at[1]],
                m.v[at[2]])
    elseif isfinite(best) && best > 0
        @printf(io, "; widest margin %.1f mm at (%.1f, %.1f) mm in the plane", best, m.u[at[1]], m.v[at[2]])
    end
    @printf(io, "; %d beams%s, %d lines of sight", n(:omega3), all(isinf, m.omega1) ? "" : " with 1ω footprints",
            n(:sight))
    return String(take!(io))
end

"""
    corridor_surface(map; tight = 1.0) -> (; vertices, faces, fills)

What the 3D view draws of a keep-out map: the corners of the cells of the sphere or plane
as `vertices`, two triangles per cell as `faces` (indices into the vertices) and the
drawing class of each triangle as `fills` (`Briareus.corridor_fills`, 1 to 5, colored by
[`CORRIDOR_COLORS`](@ref)). Each sample of a plane is the center of a square cell.
"""
function corridor_surface(m::CorridorMap; tight::Real = 1.0)
    m.surface === :plane && return _plane_surface(m, corridor_fills(m; tight))
    nθ, nφ = length(m.u), length(m.v)
    e1, e2, e3 = m.axes
    corner(i, j) = (θ = 180 * i / nθ; φ = 360 * j / nφ;
                    m.origin + m.radius * (sind(θ) * (cosd(φ) * e1 + sind(φ) * e2) + cosd(θ) * e3))
    vertices = [corner(i, j) for i in 0:nθ for j in 0:nφ-1]
    id(i, j) = i * nφ + mod(j, nφ) + 1
    cls = corridor_fills(m; tight)
    faces, fills = NTuple{3,Int}[], Int[]
    for i in 1:nθ, j in 1:nφ
        a, b, c, d = id(i - 1, j - 1), id(i, j - 1), id(i, j), id(i - 1, j)
        i < nθ && (push!(faces, (a, b, c)); push!(fills, cls[i, j]))     # no sliver at a pole
        i > 1 && (push!(faces, (a, c, d)); push!(fills, cls[i, j]))
    end
    return (; vertices, faces, fills)
end

function _plane_surface(m::CorridorMap, cls)
    n = length(m.u)
    w = (m.u[end] - m.u[1]) / (n - 1)
    e1, e2, _ = m.axes
    at(x) = m.u[1] - w / 2 + w * x
    vertices = [m.origin + at(i) * e1 + at(j) * e2 for i in 0:n for j in 0:n]
    id(i, j) = i * (n + 1) + j + 1
    faces, fills = NTuple{3,Int}[], Int[]
    for i in 1:n, j in 1:n
        a, b, c, d = id(i - 1, j - 1), id(i, j - 1), id(i, j), id(i - 1, j)
        push!(faces, (a, b, c), (a, c, d))
        push!(fills, cls[i, j], cls[i, j])
    end
    return (; vertices, faces, fills)
end

"""
    corridor_extent(map) -> Float64

The distance (mm) from TCC within which a keep-out map lies: the radius of a sphere about
TCC, or the far corner of a plane. The 3D view fits the camera to it.
"""
corridor_extent(m::CorridorMap) =
    m.surface === :sphere ? norm(m.origin) + m.radius : norm(m.origin) + sqrt(2) * maximum(abs, m.u)

"""
    corridor_drawing(map; tight = 1.0) -> String

A keep-out map as the JSON text that the Clearance tab draws: the `surface` (`"sphere"` or
`"plane"`), the numbers of samples `nu` and `nv` (polar angle and azimuth on a sphere, the
two axes of a plane), the half-width `extent` of a plane in mm, and the drawing class of
each sample (`Briareus.corridor_fills`, 1 to 5) as a string of digits `fills`, the first
index fastest. A sphere is drawn as the two hemispheres about its axis, each with the
polar angle as the distance from the center of a disk, as `Briareus.corridor_map_svg`
draws it.
"""
function corridor_drawing(m::CorridorMap; tight::Real = 1.0)
    cls = corridor_fills(m; tight)
    e = m.surface === :plane ? @sprintf("%.4f", m.u[end]) : "0"
    return "{\"surface\":\"" * String(m.surface) * "\",\"nu\":" * string(length(m.u)) * ",\"nv\":" *
           string(length(m.v)) * ",\"extent\":" * e * ",\"fills\":\"" * join(cls) * "\"}"
end
