# This file maps the keep-outs of a shot over a sphere about a point or over a plane: the
# 3ω beam envelopes, the 1ω footprints and the lines of sight of the diagnostics. The open
# region of the map is where new hardware may be placed. The beam distances are those of
# the TaLIS check (`talis.jl`), and no target is needed.

"""
    SightLine(name, theta, phi; field_of_view, center = (0, 0, 0), distance = Inf,
              aperture = field_of_view, through = false)

The volume that a diagnostic at polar angles `theta`, `phi` (degrees) views: the solid of
revolution about the line from `center` toward the instrument whose radius is
`field_of_view / 2` at `center` and `aperture / 2` at `distance` (all in mm). With the
default `distance = Inf` the volume is a cylinder, and with a pinhole at a finite distance
it is a cone that narrows toward the instrument. The volume starts at the plane through
`center` normal to the line and ends at `distance`. With `through = true` it also extends
behind `center` as a cylinder, as for a backlit radiograph whose source lies behind the
object. [`sight_line`](@ref) builds one from a diagnostic definition, and
[`sight_distance`](@ref) measures from it.
"""
struct SightLine
    name::String
    axis::Vec3
    center::Vec3
    radius::Float64
    end_radius::Float64
    distance::Float64
    through::Bool
end

function SightLine(name::AbstractString, theta::Real, phi::Real; field_of_view::Real, center = zero(Vec3),
                   distance::Real = Inf, aperture::Real = field_of_view, through::Bool = false)
    field_of_view >= 0 || throw(ArgumentError("field_of_view must be ≥ 0 mm"))
    aperture >= 0 || throw(ArgumentError("aperture must be ≥ 0 mm"))
    distance > 0 || throw(ArgumentError("distance must be > 0 mm"))
    isfinite(distance) || aperture == field_of_view ||
        throw(ArgumentError("an aperture different from the field of view needs a finite distance"))
    return SightLine(String(name), port_vector(theta, phi), Vec3(center...), field_of_view / 2, aperture / 2,
                     Float64(distance), through)
end

"""
    sight_line(id; paths = diagnostic_paths(), kwargs...)
    sight_line(def::DiagnosticDefinition; name = def.name, theta, phi, field_of_view, center,
               distance, aperture, through = false)

The [`SightLine`](@ref) of a diagnostic definition. The direction comes from the
`line_of_sight` table of the definition, the field of view from `view.field_of_view_mm` or
else from the detector size over the magnification, and the center from `view.center_mm`.
A definition with an `optic.distance_mm` gives a cone that narrows to the pinhole diameter
at that distance, and any other definition gives a cylinder. A keyword replaces the value
of the definition, and a value that the definition lacks, such as the azimuth of a line of
sight whose port is not published, must be given as a keyword.
"""
function sight_line(id::AbstractString; paths = diagnostic_paths(), kwargs...)
    issues = DiagnosticIssue[]
    raws = _read_all(paths, issues)
    defs = _checked(raws, issues)
    haskey(defs, id) && return sight_line(defs[id]; kwargs...)
    # a definition that lacks values of its instrument may still give a line of sight
    def = haskey(raws, id) ? _resolve(id, raws, issues) : nothing
    def === nothing && throw(ArgumentError("no diagnostic \"$id\" on the search path"))
    return sight_line(def; kwargs...)
end

function sight_line(def::DiagnosticDefinition; name = def.name, theta = nothing, phi = nothing,
                    field_of_view = nothing, center = nothing, distance = nothing, aperture = nothing,
                    through::Bool = false)
    d = def.data
    need(what) = throw(ArgumentError("diagnostic \"$(def.id)\" does not give $what: pass it as a keyword"))
    if theta === nothing || phi === nothing
        los = get(d, "line_of_sight", Dict{String,Any}())
        if haskey(los, "port")
            θ, φ = _line_of_sight(d)
        else
            θ, φ = get(los, "theta_deg", nothing), get(los, "phi_deg", nothing)
        end
        theta = something(theta, θ, Some(nothing))
        phi = something(phi, φ, Some(nothing))
        theta === nothing && need("`theta`")
        phi === nothing && need("`phi`")
    end
    v = get(d, "view", Dict{String,Any}())
    o = get(d, "optic", Dict{String,Any}())
    det = get(d, "detector", Dict{String,Any}())
    if field_of_view === nothing
        if haskey(v, "field_of_view_mm")
            field_of_view = Float64(v["field_of_view_mm"])
        elseif haskey(det, "pixel_um") && haskey(det, "pixels") && haskey(o, "magnification")
            field_of_view = 1e-3 * det["pixel_um"] * maximum(det["pixels"]) / o["magnification"]
        else
            need("`field_of_view`")
        end
    end
    if center === nothing
        c = _view_center(v)
        c isa Symbol && need("`center` as a point (the view is centered on part :$c)")
        center = c
    end
    if distance === nothing
        distance = haskey(o, "distance_mm") ? Float64(o["distance_mm"]) : Inf
    end
    if aperture === nothing
        aperture = isfinite(distance) ? (haskey(o, "diameter_um") ? 1e-3 * o["diameter_um"] : 0.0) : field_of_view
    end
    return SightLine(String(name), theta, phi; field_of_view, center, distance, aperture, through)
end

# distance from (s, ρ) to the segment from (s0, ρ0) to (s1, ρ1)
function _segment_distance(s, ρ, s0, ρ0, s1, ρ1)
    ds, dρ = s1 - s0, ρ1 - ρ0
    l2 = ds^2 + dρ^2
    τ = l2 > 0 ? clamp(((s - s0) * ds + (ρ - ρ0) * dρ) / l2, 0.0, 1.0) : 0.0
    return hypot(s - s0 - τ * ds, ρ - ρ0 - τ * dρ)
end

const _SIGHT_FAR = 1e7        # mm, the length that stands for an unbounded sight line

"""
    sight_distance(line::SightLine, x) -> Float64

The signed distance (mm) from point `x` to the surface of the volume of a
[`SightLine`](@ref), negative inside. The distance is exact in three dimensions: it is the
distance to the outline of the volume in the half-plane of axial position and radius.
"""
function sight_distance(l::SightLine, x)
    d = Vec3(x...) - l.center
    s = dot(d, l.axis)
    ρ = norm(d - s * l.axis)
    D = isfinite(l.distance) ? l.distance : _SIGHT_FAR
    dist = _segment_distance(s, ρ, 0.0, l.radius, D, l.end_radius)
    isfinite(l.distance) && (dist = min(dist, _segment_distance(s, ρ, D, 0.0, D, l.end_radius)))
    r = l.radius + (l.end_radius - l.radius) * clamp(s / D, 0.0, 1.0)
    if l.through
        dist = min(dist, _segment_distance(s, ρ, -_SIGHT_FAR, l.radius, 0.0, l.radius))
        inside = s <= D && ρ < r
    else
        dist = min(dist, _segment_distance(s, ρ, 0.0, 0.0, 0.0, l.radius))
        inside = 0 <= s <= D && ρ < r
    end
    return inside ? -dist : dist
end

# The keep-outs of a set of beams and sight lines, ready for point queries.
struct _CorridorContext{C}
    envs::Vector{BeamEnvelope}
    aims::Vector{Vec3}
    shifts::Vector{Float64}                       # R + focus × cone of each beam
    beamlets::Vector{Union{Nothing,Vector{NTuple{3,Vec3}}}}
    sights::Vector{SightLine}
    cx::C
    clearance_3w::Float64
    omega1_margin::Float64
    diagnostic::Bool                              # 1ω requirement of a diagnostic, 10 mm + 0.015 r
    sight_clearance::Float64
    sources::Vector{Tuple{Symbol,String}}
end

function _corridor_context(beams, sight_lines; tolerance = 0.5, pivot = :tcc, pointing = 0.05,
                           target_offset = 0.025, hardware_tolerance = 0.0, focus = 1.0, ssd = 0.0,
                           contour_3w = 0.01, clearance_3w = 0.5, omega1::Bool = true, omega1_margin = 0.2,
                           omega1_edge = 4.45, omega1_size = nothing, omega1_focus = nothing,
                           lens_distance = NIF_LENS_DISTANCE, offsets = nothing, sight_clearance = 0.0)
    omega1_margin === :diagnostic || omega1_margin isa Real ||
        throw(ArgumentError("omega1_margin must be a distance in mm or :diagnostic"))
    cx = _omega1_setup(tolerance, pivot, pointing, target_offset, hardware_tolerance, focus, omega1_edge,
                       omega1_size, omega1_focus, lens_distance)
    vb = [_value_beam(b) for b in beams]
    table = any(_is_nif, vb) ? _nif_beam_table() : nothing
    envs = BeamEnvelope[BeamEnvelope(b; peak_fraction = contour_3w, margin = ssd) for b in vb]
    bls = Union{Nothing,Vector{NTuple{3,Vec3}}}[omega1 ? _omega1_beamlets(b, offsets, table) : nothing for b in vb]
    sights = SightLine[l for l in sight_lines]
    sources = Tuple{Symbol,String}[]
    for class in (:omega3, :omega1), b in vb
        push!(sources, (class, String(b.name)))
    end
    append!(sources, [(:sight, l.name) for l in sights])
    return _CorridorContext(envs, Vec3[Vec3(b.aim) for b in vb], [cx.R + cx.focus * e.cone for e in envs], bls,
                            sights, cx, Float64(clearance_3w),
                            omega1_margin === :diagnostic ? 0.0 : Float64(omega1_margin),
                            omega1_margin === :diagnostic, Float64(sight_clearance), sources)
end

# The margins of point `x` to the three classes of keep-out and the index, in the
# sources, of the keep-out with the smallest margin (0 with no keep-out).
function _corridor_point(c::_CorridorContext, x::Vec3)
    nb = length(c.envs)
    m3, m1, ms, best, ib = Inf, Inf, Inf, Inf, 0
    req1 = c.diagnostic ? 10 + 0.015 * norm(x) : c.omega1_margin
    for i in 1:nb
        rot = c.cx.tol(x, c.aims[i])
        d = envelope_distance(c.envs[i], x)[1] - c.shifts[i] - rot - c.clearance_3w
        m3 = min(m3, d)
        d < best && ((best, ib) = (d, i))
        bl = c.beamlets[i]
        bl === nothing && continue
        d = _omega1_distance(x, bl, c.cx)[2] - c.cx.R - rot - req1
        m1 = min(m1, d)
        d < best && ((best, ib) = (d, nb + i))
    end
    for (j, l) in enumerate(c.sights)
        d = sight_distance(l, x) - c.sight_clearance
        ms = min(ms, d)
        d < best && ((best, ib) = (d, 2nb + j))
    end
    return m3, m1, ms, ib
end

"""
    corridor_margin(beams, x; sight_lines = SightLine[], kwargs...)
        -> (; margin, omega3, omega1, sight, nearest)

The margins (mm) of point `x` to the keep-outs of a shot, as [`corridor_map`](@ref)
evaluates them at each sample and with the same keywords: to the 3ω envelopes, to the 1ω
footprints and to the lines of sight, the smallest of the three as `margin`, and the
keep-out that sets it as `nearest`, a pair of its class (`:omega3`, `:omega1` or `:sight`)
and the name of the beam or sight line (`nothing` with no keep-out). A margin is negative
inside a keep-out and `Inf` for a class with no member.
"""
function corridor_margin(beams::AbstractVector{<:Beam}, x; sight_lines = SightLine[], kwargs...)
    c = _corridor_context(beams, sight_lines; kwargs...)
    m3, m1, ms, ib = _corridor_point(c, Vec3(x...))
    return (; margin = min(m3, m1, ms), omega3 = m3, omega1 = m1, sight = ms,
            nearest = ib == 0 ? nothing : c.sources[ib])
end

"""
    CorridorMap

The keep-outs of a shot sampled over a sphere or a plane, as [`corridor_map`](@ref)
returns them.

- `surface` — `:sphere` or `:plane`
- `origin`, `radius` — the center and radius (mm) of the sphere, or a point of the plane
  and `NaN`
- `axes` — three orthonormal vectors: the polar axis of the sphere last, or the two
  in-plane axes and the normal of the plane
- `u`, `v` — the sample coordinates: polar angle θ and azimuth φ (degrees) about the polar
  axis, or the in-plane coordinates `p`, `q` (mm)
- `points` — the sample points, `points[i, j]` at `(u[i], v[j])`
- `omega3`, `omega1`, `sight` — the margin (mm) of each sample to the 3ω envelopes, the 1ω
  footprints and the lines of sight, negative inside and `Inf` for a class with no member
- `margin` — the smallest of the three
- `nearest` — the index in `sources` of the keep-out that sets `margin` (0 with none)
- `sources` — the keep-outs as pairs of class (`:omega3`, `:omega1`, `:sight`) and name

[`corridor_classes`](@ref) gives the classes that cover each sample,
[`angular_margin`](@ref) the margin on a sphere as an angle, [`corridor_open_fraction`](@ref) the
open share of the surface and [`corridor_map_svg`](@ref) a drawing.
"""
struct CorridorMap
    surface::Symbol
    origin::Vec3
    radius::Float64
    axes::NTuple{3,Vec3}
    u::Vector{Float64}
    v::Vector{Float64}
    points::Matrix{Vec3}
    omega3::Matrix{Float64}
    omega1::Matrix{Float64}
    sight::Matrix{Float64}
    margin::Matrix{Float64}
    nearest::Matrix{Int32}
    sources::Vector{Tuple{Symbol,String}}
end

function Base.show(io::IO, m::CorridorMap)
    print(io, "CorridorMap(", m.surface, ", ", length(m.u), "×", length(m.v), " samples, ",
          @sprintf("%.1f%%", 100 * corridor_open_fraction(m)), " open)")
end

# two unit vectors that complete `n` to a right-handed frame (e1, e2, n)
function _corridor_axes(n::Vec3)
    n = normalize(n)
    z = Vec3(0.0, 0.0, 1.0)
    abs(dot(n, z)) > 1 - 1e-12 && return (Vec3(1.0, 0.0, 0.0), Vec3(0.0, sign(n[3]), 0.0), n)
    e1 = normalize(cross(z, n))
    return (e1, cross(n, e1), n)
end

"""
    corridor_map(beams; radius, center = (0, 0, 0), axis = (0, 0, 1), ntheta = 90, nphi = 180, kwargs...)
    corridor_map(beams; plane = (point, normal), extent, n = 121, kwargs...) -> CorridorMap

Map where new hardware, such as stalks, shields, foils and diagnostic snouts, may be
placed: the margin of each sample of a surface to the keep-outs of the shot. The surface
is the sphere of `radius` (mm) about `center`, sampled at the centers of `ntheta` by
`nphi` cells in the polar angle and azimuth about `axis`, or the square of half-width
`extent` (mm) about `point` in the plane with the given `normal`, sampled on `n` by `n`
points. The in-plane axes of a plane are `ẑ × normal` and its cross product with the
normal (`x̂` and `±ŷ` for a normal along `ẑ`).

The keep-outs are those of the TaLIS check ([`talis_check`](@ref)), measured the same way
and with the same keywords and defaults:

- 3ω light — the envelope of each beam at `contour_3w = 0.01` of the peak
  ([`BeamEnvelope`](@ref)), grown by `ssd`, the position errors (`pointing`,
  `target_offset`, `hardware_tolerance`), the `focus` offset and the rotation `tolerance`
  about `pivot`, with the required `clearance_3w = 0.5` mm
- 1ω light — the footprint squares of each beam ([`omega1_contour`](@ref), with
  `omega1_edge`, `omega1_size`, `omega1_focus`, `lens_distance` and `offsets`) at focus
  offsets 0 and `±focus`, grown by the same position errors, with the required
  `omega1_margin = 0.2` mm. `omega1_margin = :diagnostic` applies the requirement of a
  diagnostic, 10 mm at TCC growing to 25 mm at 1 m. `omega1 = false` leaves the class out.
- lines of sight — the volumes of `sight_lines` ([`SightLine`](@ref)), with the required
  `sight_clearance = 0` mm

A beam distance is measured in the plane through the sample normal to the beam, as the
rules measure it, and a sight-line distance in three dimensions. A margin is the distance
less the errors and the required clearance, so a sample is open where its margin is
positive. Beams upstream and downstream of best focus both count, and hardware that would
shadow a keep-out is not considered. The map is computed without a target and applies to
any kind of experiment. `beams` may be empty.
"""
function corridor_map(beams::AbstractVector{<:Beam}; sight_lines = SightLine[], radius = nothing,
                      center = zero(Vec3), axis = Vec3(0.0, 0.0, 1.0), ntheta::Integer = 90, nphi::Integer = 180,
                      plane = nothing, extent = nothing, n::Integer = 121, kwargs...)
    c = _corridor_context(beams, sight_lines; kwargs...)
    if plane === nothing
        radius === nothing && throw(ArgumentError("give `radius` for a sphere or `plane` and `extent` for a plane"))
        radius > 0 || throw(ArgumentError("radius must be > 0 mm"))
        ntheta >= 2 && nphi >= 4 || throw(ArgumentError("ntheta must be ≥ 2 and nphi ≥ 4"))
        o = Vec3(center...)
        ax = _corridor_axes(Vec3(axis...))
        u = [(i - 0.5) * 180 / ntheta for i in 1:ntheta]
        v = [(j - 0.5) * 360 / nphi for j in 1:nphi]
        pts = [o + radius * (sind(θ) * (cosd(φ) * ax[1] + sind(φ) * ax[2]) + cosd(θ) * ax[3]) for θ in u, φ in v]
        return _corridor_fill(c, :sphere, o, Float64(radius), ax, u, v, pts)
    end
    radius === nothing || throw(ArgumentError("give one of `radius` and `plane`"))
    extent !== nothing && extent > 0 || throw(ArgumentError("a plane needs `extent` > 0 mm"))
    n >= 2 || throw(ArgumentError("n must be ≥ 2"))
    o = Vec3(plane[1]...)
    ax = _corridor_axes(Vec3(plane[2]...))
    u = collect(range(-Float64(extent), Float64(extent); length = n))
    pts = [o + p * ax[1] + q * ax[2] for p in u, q in u]
    return _corridor_fill(c, :plane, o, NaN, ax, u, copy(u), pts)
end

function _corridor_fill(c::_CorridorContext, surface, o, R, ax, u, v, pts::Matrix{Vec3})
    m3, m1, ms = similar(pts, Float64), similar(pts, Float64), similar(pts, Float64)
    near = similar(pts, Int32)
    Threads.@threads for i in eachindex(pts)
        m3[i], m1[i], ms[i], near[i] = _corridor_point(c, pts[i])
    end
    return CorridorMap(surface, o, R, ax, u, v, pts, m3, m1, ms, min.(m3, m1, ms), near, c.sources)
end

"""
    corridor_classes(map) -> Matrix{UInt8}

The classes of keep-out that cover each sample of a [`CorridorMap`](@ref), as bits: `0x01`
for 3ω light, `0x02` for 1ω light and `0x04` for a line of sight. An open sample is `0x00`.
"""
corridor_classes(m::CorridorMap) =
    UInt8.(m.omega3 .<= 0) .| (UInt8.(m.omega1 .<= 0) .<< 1) .| (UInt8.(m.sight .<= 0) .<< 2)

"""
    angular_margin(map) -> Matrix{Float64}

The margins of a spherical [`CorridorMap`](@ref) as angles (degrees) at the center of the
sphere: `asind(margin / radius)`, the angle through which a sample moves along the sphere
to cover the margin. The angle is exact for a keep-out that is a cone with its apex at the
center, and margins beyond the radius give ±90°.
"""
function angular_margin(m::CorridorMap)
    m.surface === :sphere || throw(ArgumentError("angular margins are defined on a sphere"))
    return asind.(clamp.(m.margin ./ m.radius, -1.0, 1.0))
end

# the weight of each sample: the solid angle of its cell on a sphere, uniform on a plane
_corridor_weights(m::CorridorMap) =
    m.surface === :sphere ? [sind(θ) for θ in m.u, _ in m.v] : ones(size(m.margin))

"""
    corridor_open_fraction(map; margin = 0.0) -> Float64

The share of the surface of a [`CorridorMap`](@ref) whose margin exceeds `margin` (mm):
of the solid angle on a sphere, or of the area on a plane.
"""
function corridor_open_fraction(m::CorridorMap; margin::Real = 0.0)
    w = _corridor_weights(m)
    return sum(w[m.margin .> margin]) / sum(w)
end

# the fill of a sample: by the class that covers it (3ω first, then 1ω, then a line of
# sight), or open, with a second shade below the `tight` margin
function _corridor_fill_index(m3, m1, ms, tight)
    m3 <= 0 && return 1
    m1 <= 0 && return 2
    ms <= 0 && return 3
    return min(m3, m1, ms) < tight ? 4 : 5
end

const _CORRIDOR_FILLS = ("#0072B2", "#D55E00", "#CC79A7", "#F0E442", "#B7E3D4")
const _CORRIDOR_LABELS = ("3ω light", "1ω light", "line of sight", "open, tight", "open")

"""
    corridor_fills(map; tight = 1.0) -> Matrix{Int}

The drawing class of each sample of a [`CorridorMap`](@ref): 1 inside 3ω light, 2 inside
1ω light, 3 inside a line of sight (the first that applies), 4 open with a margin below
`tight` (mm) and 5 open. [`corridor_map_svg`](@ref) fills the samples by this class.
"""
corridor_fills(m::CorridorMap; tight::Real = 1.0) =
    [_corridor_fill_index(m.omega3[i], m.omega1[i], m.sight[i], tight) for i in CartesianIndices(m.margin)]

"""
    corridor_map_svg(map; size = 480, tight = 1.0) -> String

A vector drawing (SVG) of a [`CorridorMap`](@ref). A sphere is drawn as two polar disks
`size` pixels across, the hemisphere about the polar axis on the left and the opposite one
on the right. The radius on a disk is proportional to the angle from the pole of that
hemisphere, and the azimuth runs counterclockwise from the right on both, so the right
disk shows its hemisphere as seen through the sphere from the side of the polar axis. A
plane is drawn as one square with `p` to the right and `q` up. The samples are filled by
[`corridor_fills`](@ref): 3ω light, 1ω light, a line of sight, and the open region in two
shades, the darker one where the margin is below `tight` (mm). Rings mark every 30° of
polar angle. The text of the drawing is the same for the same map: coordinates are written
with two decimals.
"""
function corridor_map_svg(m::CorridorMap; size::Integer = 480, tight::Real = 1.0)
    size >= 100 || throw(ArgumentError("size must be at least 100 pixels"))
    fills = corridor_fills(m; tight)
    io = IOBuffer()
    f(x) = @sprintf("%.2f", x)
    sphere = m.surface === :sphere
    W, H = sphere ? 2size + 16 : size, size + 44
    println(io, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"", W, "\" height=\"", H, "\" viewBox=\"0 0 ", W,
            " ", H, "\" font-family=\"sans-serif\" font-size=\"12\">")
    println(io, "<rect width=\"", W, "\" height=\"", H, "\" fill=\"#ffffff\"/>")
    if sphere
        _corridor_disks!(io, m, fills, size, f)
    else
        _corridor_square!(io, m, fills, size, f)
    end
    x = 8.0
    for (c, label) in zip(_CORRIDOR_FILLS, _CORRIDOR_LABELS)
        println(io, "<rect x=\"", f(x), "\" y=\"", f(size + 26.0), "\" width=\"12\" height=\"12\" fill=\"", c,
                "\" stroke=\"#000000\" stroke-width=\"0.5\"/>")
        println(io, "<text x=\"", f(x + 16), "\" y=\"", f(size + 37.0), "\">", label, "</text>")
        x += 28 + 7 * length(label)
    end
    title = sphere ? @sprintf("sphere of %g mm: %.1f%% open", m.radius, 100 * corridor_open_fraction(m)) :
            @sprintf("plane, ±%g mm: %.1f%% open", m.u[end], 100 * corridor_open_fraction(m))
    println(io, "<text x=\"8.00\" y=\"", f(size + 18.0), "\">", title, "</text>")
    print(io, "</svg>\n")
    return String(take!(io))
end

# runs of equal fill along a row, as (first, last, fill)
function _fill_runs(row, maxlen)
    runs = Tuple{Int,Int,Int}[]
    a = 1
    for j in 2:length(row)+1
        if j > length(row) || row[j] != row[a] || j - a >= maxlen
            push!(runs, (a, j - 1, row[a]))
            a = j
        end
    end
    return runs
end

function _corridor_disks!(io, m::CorridorMap, fills, size, f)
    nθ, nφ = length(m.u), length(m.v)
    dθ, dφ = 180 / nθ, 360 / nφ
    Rpx = size / 2 - 4
    for south in (false, true)
        cx = (south ? size + 16 : 0) + size / 2
        cy = size / 2
        xy(r, φ) = (cx + r * cosd(φ), cy - r * sind(φ))
        for i in 1:nθ
            θ0, θ1 = (i - 1) * dθ, i * dθ
            south && ((θ0, θ1) = (180 - θ1, 180 - θ0))
            θ0 >= 90 - 1e-9 && continue
            r0, r1 = Rpx * θ0 / 90, Rpx * min(θ1, 90.0) / 90
            # runs of at most a third of a turn, so that each arc is a minor arc
            for (a, b, c) in _fill_runs(view(fills, i, :), max(1, nφ ÷ 3))
                φ0, φ1 = (a - 1) * dφ, b * dφ
                p1, p2, p3, p4 = xy(r1, φ0), xy(r1, φ1), xy(r0, φ1), xy(r0, φ0)
                println(io, "<path d=\"M", f(p1[1]), " ", f(p1[2]), "A", f(r1), " ", f(r1), " 0 0 0 ", f(p2[1]), " ",
                        f(p2[2]), "L", f(p3[1]), " ", f(p3[2]), "A", f(r0), " ", f(r0), " 0 0 1 ", f(p4[1]), " ",
                        f(p4[2]), "Z\" fill=\"", _CORRIDOR_FILLS[c], "\" stroke=\"", _CORRIDOR_FILLS[c],
                        "\" stroke-width=\"0.5\"/>")
            end
        end
        for θ in (30, 60, 90)
            println(io, "<circle cx=\"", f(cx), "\" cy=\"", f(cy), "\" r=\"", f(Rpx * θ / 90),
                    "\" fill=\"none\" stroke=\"#000000\" stroke-width=\"0.5\" stroke-opacity=\"0.5\"/>")
        end
        println(io, "<line x1=\"", f(cx), "\" y1=\"", f(cy), "\" x2=\"", f(cx + Rpx), "\" y2=\"", f(cy),
                "\" stroke=\"#000000\" stroke-width=\"0.5\" stroke-opacity=\"0.5\"/>")
        println(io, "<text x=\"", f(cx - Rpx), "\" y=\"14.00\">", south ? "θ from 180° (center) to 90°" :
                "θ from 0° (center) to 90°", "</text>")
        println(io, "<text x=\"", f(cx + Rpx - 34), "\" y=\"", f(cy - 4), "\">φ = 0</text>")
    end
    return io
end

function _corridor_square!(io, m::CorridorMap, fills, size, f)
    n = length(m.u)
    w = size / n
    for j in 1:n, (a, b, c) in _fill_runs(view(fills, :, j), n)
        println(io, "<rect x=\"", f((a - 1) * w), "\" y=\"", f(size - j * w), "\" width=\"", f((b - a + 1) * w),
                "\" height=\"", f(w), "\" fill=\"", _CORRIDOR_FILLS[c], "\" stroke=\"", _CORRIDOR_FILLS[c],
                "\" stroke-width=\"0.5\"/>")
    end
    println(io, "<rect x=\"0.00\" y=\"0.00\" width=\"", f(Float64(size)), "\" height=\"", f(Float64(size)),
            "\" fill=\"none\" stroke=\"#000000\" stroke-width=\"0.5\"/>")
    return io
end

"""
    write_corridor_map_svg(path, map; size = 480, tight = 1.0) -> path

Write the drawing of a [`CorridorMap`](@ref) ([`corridor_map_svg`](@ref)) to `path`. The
drawing is plain text and is made without a graphics context.
"""
function write_corridor_map_svg(path::AbstractString, m::CorridorMap; size::Integer = 480, tight::Real = 1.0)
    write(path, corridor_map_svg(m; size, tight))
    return path
end
