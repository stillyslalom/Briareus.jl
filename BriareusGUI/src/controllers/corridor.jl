# This file covers the keep-out map of the 3D view: the sphere about TCC colored by the
# light and the lines of sight that cross it, with the open region where new hardware fits.

"""
    CorridorState

The keep-out map of a session (`Briareus.corridor_map`). `radius` (mm) is the radius of
the sphere about TCC, `shown` whether the 3D view draws the map, `map` the map of the last
update (a core `CorridorMap`, or `nothing`), `status` a line about it and `inputs` the
beams, sight lines, radius and settings of that update.
"""
struct CorridorState
    radius::Observable{Float64}
    shown::Observable{Bool}
    map::Observable{Union{Nothing,CorridorMap}}
    status::Observable{String}
    inputs::Base.RefValue{Any}
end

CorridorState() = CorridorState(Observable(50.0), Observable(false), Observable{Union{Nothing,CorridorMap}}(nothing),
                                Observable(""), Ref{Any}(nothing))

"""
Colors (hex) of the keep-out map by drawing class (`Briareus.corridor_fills`): 3ω light,
1ω light, a line of sight, open with a tight margin, and open.
"""
const CORRIDOR_COLORS = ("#0072B2", "#D55E00", "#CC79A7", "#F0E442", "#009E73")

"Opacity of each drawing class of the keep-out map: the open region stays faint."
const CORRIDOR_OPACITY = (0.55, 0.55, 0.55, 0.4, 0.12)

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
    update_corridor!(cs::CorridorState, rv, tl, beams, sight_lines) -> cs

Compute the keep-out map of `beams` and `sight_lines` on the sphere of `cs.radius` about
TCC, with the settings of [`corridor_kwargs`](@ref), and leave it in `cs.map`. The map of
192 beams takes a fraction of a second, so it is computed in place.
"""
function update_corridor!(cs::CorridorState, rv::ReviewState, tl::TaLISState, beams::AbstractVector,
                          sight_lines::AbstractVector)
    kw = corridor_kwargs(rv, tl)
    bs = collect(beams)
    sl = collect(SightLine, sight_lines)
    m = corridor_map(bs; sight_lines = sl, radius = cs.radius[], kw...)
    cs.inputs[] = (bs, sl, cs.radius[], kw)
    cs.status[] = corridor_text(m; tight = tl.tight[])
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
beams, the sight lines, the radius or a setting.
"""
function corridor_stale(cs::CorridorState, rv::ReviewState, tl::TaLISState, beams::AbstractVector,
                        sight_lines::AbstractVector)
    inp = cs.inputs[]
    inp === nothing && return false
    return inp[1] != beams || inp[2] != sight_lines || inp[3] != cs.radius[] || inp[4] != corridor_kwargs(rv, tl)
end

"""
    corridor_text(map; tight = 1.0) -> String

One line about a keep-out map: the open share of the sphere, the share with more than the
`tight` margin, the widest margin and where it lies, and the number of beams and lines of
sight.
"""
function corridor_text(m::CorridorMap; tight::Real = 1.0)
    n(class) = count(s -> first(s) === class, m.sources)
    io = IOBuffer()
    @printf(io, "sphere of %g mm about TCC: %.1f%% open, %.1f%% with more than %g mm of margin", m.radius,
            100 * corridor_open_fraction(m), 100 * corridor_open_fraction(m; margin = tight), tight)
    best, at = findmax(m.margin)
    if isfinite(best) && best > 0
        @printf(io, "; widest margin %.1f mm (%.1f°) at θ %.0f°, φ %.0f°", best, angular_margin(m)[at], m.u[at[1]],
                m.v[at[2]])
    end
    @printf(io, "; %d beams%s, %d lines of sight", n(:omega3), all(isinf, m.omega1) ? "" : " with 1ω footprints",
            n(:sight))
    return String(take!(io))
end

"""
    corridor_surface(map; tight = 1.0) -> (; vertices, faces, fills)

What the 3D view draws of a spherical keep-out map: the corners of the cells of the sphere
as `vertices`, two triangles per cell as `faces` (indices into the vertices) and the
drawing class of each triangle as `fills` (`Briareus.corridor_fills`, 1 to 5, colored by
[`CORRIDOR_COLORS`](@ref)).
"""
function corridor_surface(m::CorridorMap; tight::Real = 1.0)
    m.surface === :sphere || throw(ArgumentError("the 3D view draws a spherical map"))
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
