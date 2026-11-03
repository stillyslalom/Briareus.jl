# Ready-made target geometries.

"""
    WallWindow(; phi, z = 0.0, width, height, material = nothing, open = false, name = nothing,
               gap = 0.0, gap_material = nothing)
    WallWindow(; phi, z = 0.0, diameter, segments = 0, ...)

Diagnostic window cut into a hohlraum wall (see `hohlraum`), centered at azimuth `phi`
(degrees from `+x`) and axial position `z` (mm). Its shape is drawn on the unrolled
wall (arc length along the wall by `z`, an isometry of the cylinder, so areas are those
of the flat shape):

- a rectangle `width` mm of arc by `height` mm along the axis, or
- a disk of `diameter` mm, area π d²/4. Its edge is a regular polygon of `segments`
  sides (default: about one per wall cell, a multiple of 8, at least 32) enlarged to
  the disk's area, so the window's area is exact while its edge stays within
  `π² d/6 segments²` of the circle (3 µm for a 2 mm window with 32 sides).

- A closed window (`open = false`) is an opaque patch of wall that re-emits with its
  own `material` (albedo and laser conversion, e.g. a thin gold foil whose albedo
  follows `BurnthroughAlbedo`); `material = nothing` keeps the wall's.
- An open window (`open = true`) is a hole: an aperture that laser light passes
  through and whose incident x-ray flux escapes.

`gap > 0` (mm) surrounds the window with an assembly gap of that width (a ring around
a disk, a frame around a rectangle), its own component `Symbol(name, :_gap)`: open by
default, or a closed patch of `gap_material`. `name` defaults to `:window_<k>` for the
`k`-th window.
"""
struct WallWindow
    phi::Float64
    z::Float64
    width::Float64                 # a disk's diameter
    height::Float64                # a disk's diameter
    circular::Bool
    segments::Int                  # sides of a disk's polygon (0: automatic)
    material::Union{Nothing,Material}
    open::Bool
    name::Union{Nothing,Symbol}
    gap::Float64
    gap_material::Union{Nothing,Material}
    function WallWindow(phi, z, width, height, circular, segments, material, open, name, gap,
                        gap_material)
        width > 0 && height > 0 || throw(ArgumentError("window width and height must be positive"))
        circular && width != height && throw(ArgumentError("a circular window has a single diameter"))
        open && material !== nothing && throw(ArgumentError("an open window has no material"))
        gap >= 0 || throw(ArgumentError("the window gap must not be negative"))
        segments == 0 || segments >= 8 || throw(ArgumentError("a circular window needs at least 8 segments"))
        new(phi, z, width, height, circular, segments, material, open,
            name === nothing ? nothing : Symbol(name), gap, gap_material)
    end
end
function WallWindow(; phi, z = 0.0, width = nothing, height = nothing, diameter = nothing,
                    material = nothing, open::Bool = false, name = nothing, gap = 0.0,
                    gap_material = nothing, segments::Integer = 0)
    if diameter === nothing
        (width === nothing || height === nothing) &&
            throw(ArgumentError("give a window's `diameter`, or its `width` and `height`"))
        return WallWindow(phi, z, width, height, false, 0, material, open, name, gap, gap_material)
    end
    width === nothing && height === nothing ||
        throw(ArgumentError("give a window's `diameter` or its `width` and `height`, not both"))
    return WallWindow(phi, z, diameter, diameter, true, segments, material, open, name, gap,
                      gap_material)
end

# Azimuthal (radians: start, angular width) and axial (mm) extent of the rectangle of
# half-sizes `a` (mm of arc) by `b` (mm) centered on window `w`, on a wall of radius R.
_rect_extent(w::WallWindow, R, a, b) = (deg2rad(w.phi) - a / R, 2a / R, w.z - b, w.z + b)

# Rectangles whose edges the wall grid must contain: a rectangular window and the outer
# edge of its gap frame; for a disk, the square box (the disk, its gap, and half a wall
# cell `hcell` of margin) that `_circular_window` re-meshes. The last is the footprint.
function _window_extents(w::WallWindow, R, hcell)
    if w.circular
        B = w.width / 2 + w.gap + hcell / 2
        return [_rect_extent(w, R, B, B)]
    end
    inner = _rect_extent(w, R, w.width / 2, w.height / 2)
    w.gap > 0 || return [inner]
    return [inner, _rect_extent(w, R, w.width / 2 + w.gap, w.height / 2 + w.gap)]
end

# Whether the point `c` (a patch centroid) lies inside the extent.
function _in_window((φa, Δφ, za, zb), c)
    za < c[3] < zb || return false
    return mod(atan(c[2], c[1]) - φa, 2π) < Δφ
end

function _check_windows(windows, R, h, hcell)
    exts = [_window_extents(w, R, hcell) for w in windows]
    tol = 1e-9 * max(R, h)
    for (i, e) in enumerate(exts)
        _, Δφ, za, zb = last(e)
        Δφ <= 2π + 1e-12 || throw(ArgumentError("window $i is wider than the wall circumference"))
        -h - tol <= za && zb <= h + tol ||
            throw(ArgumentError("window $i extends past the end of the wall"))
    end
    for i in eachindex(exts), j in i+1:lastindex(exts)
        (ai, wi, zai, zbi), (aj, wj, zaj, zbj) = last(exts[i]), last(exts[j])
        zai < zbj - tol && zaj < zbi - tol || continue
        (mod(aj - ai, 2π) < wi - 1e-12 || mod(ai - aj, 2π) < wj - 1e-12) &&
            throw(ArgumentError("windows $i and $j overlap"))
    end
    return exts
end

# Band `i` of `k` of the fillet of radius `rc` joining the wall (radius R, ends ±hw) to the
# cap at z = sgn (hw + rc): the frustum between the arc points at θ = (i-1)π/2k and iπ/2k
# from the wall, on the wall's azimuthal grid, facing the cavity, with its exact cone.
function _corner_band(name, R, hw, rc, sgn, i, k, nphi, phibreaks, material)
    θa, θb = (i - 1) * π / 2k, i * π / 2k
    ρa, za = R - rc + rc * cos(θa), sgn * (hw + rc * sin(θa))
    ρb, zb = R - rc + rc * cos(θb), sgn * (hw + rc * sin(θb))
    φs = isempty(phibreaks) ? range(0.0, 2π; length = nphi + 1) :
         conforming_grid(0.0, 2π, nphi, phibreaks; periodic = true)
    m = Mesh()
    θc = (θa + θb) / 2
    for j in 1:length(φs)-1
        a, b = φs[j], φs[j+1]
        p = [Vec3(ρa * cos(a), ρa * sin(a), za), Vec3(ρa * cos(b), ρa * sin(b), za),
             Vec3(ρb * cos(b), ρb * sin(b), zb), Vec3(ρb * cos(a), ρb * sin(a), zb)]
        φc = (a + b) / 2
        # toward the fillet's center of curvature: into the cavity
        add_polygon!(m, p; facing = -Vec3(cos(θc) * cos(φc), cos(θc) * sin(φc), sgn * sin(θc)))
    end
    # the cone through both rim circles: apex on the axis beyond the cap, opening toward
    # the waist
    zapex = za + (zb - za) * ρa / (ρa - ρb)
    half = atand((ρa - ρb) / abs(zb - za))
    return Component(name, m, material; surface = ConeSurface(half; apex = (0.0, 0.0, zapex), axis = (0.0, 0.0, -sgn)))
end

"""
    equal_area_radius(r, n)

Circumradius of the regular `n`-gon with the area of a circle of radius `r`:
`r √(2π / (n sin(2π/n)))`. Its edge runs from `r π²/3n²` outside the circle at the
vertices to `r π²/6n²` inside at mid-side.
"""
equal_area_radius(r, n) = r * sqrt(2π / (n * sin(2π / n)))

# Default number of sides of a circular window's polygon: about one per wall cell of
# size `hcell` along the outer circumference, a multiple of 8, at least 32.
_window_segments(w::WallWindow, hcell) =
    w.segments > 0 ? w.segments : max(32, 8 * cld(ceil(Int, 2π * (w.width / 2 + w.gap) / hcell), 8))

"""
    _circular_window(box, R, φc, zc, rw, rg, n, hcell) -> (wall, window, gap)

Re-mesh the square `box` of whole wall patches (cut from a `cylinder_mesh` of radius
`R` whose grid contains the box edges) around a circular window centered at azimuth
`φc` (radians) and height `zc`, of radius `rw` with a gap ring out to `rg ≥ rw`. In the
unrolled coordinates (u, v) = (R(φ - φc), z - zc), with every vertex on the cylinder:

- window: polar rings of about `hcell` radial spacing, `n` sides on the edge, halving
  the count inward while the arcs stay below 1.25 `hcell`, and a fan at the center;
- gap: two triangles per side between the window's and the gap's polygons;
- wall: triangles joining the outer polygon to the box boundary, whose vertices are
  exactly those of the box's own grid (no T-junctions with the wall around it), by the
  shorter-diagonal rule.

The window and gap polygons are `equal_area_radius` polygons, so with the exact curved
areas of triangles on the cylinder (R times their (φ, z) area) the components have the
areas π rw² and π(rg² - rw²) to round-off. Every patch is a triangle: a polygon with
more than three vertices on the cylinder is not planar, and the centroid of a folded
patch can sit behind its own facets, which hides the element from the visibility rays.
"""
function _circular_window(box::Mesh, R, φc, zc, rw, rg, n::Integer, hcell)
    uv(p) = (R * (mod(atan(p[2], p[1]) - φc + π, 2π) - π), p[3] - zc)
    xyz((u, v)) = (φ = φc + u / R; Vec3(R * cos(φ), R * sin(φ), zc + v))
    ang(p) = ((u, v) = uv(p); atan(v, u))
    tri!(m, a, b, c) = (g = a + b + c; add_polygon!(m, [a, b, c]; facing = -Vec3(g[1], g[2], 0)))
    # a quad (a, b on the inner ring, d, c on the outer) split along its shorter diagonal
    quad!(m, a, b, c, d) = norm(c - a) <= norm(d - b) ? (tri!(m, a, d, c); tri!(m, a, c, b)) :
                                                         (tri!(m, a, d, b); tri!(m, b, d, c))
    # box boundary, counterclockwise from its point nearest angle -π
    B = maximum(p -> maximum(abs, uv(p)), box.vertices)
    bnd = sort!(unique(p for p in box.vertices if maximum(abs, uv(p)) > B * (1 - 1e-9)); by = ang)
    θ0 = ang(bnd[1])
    ring(r, m) = [xyz((r * cos(θ0 + 2π * k / m), r * sin(θ0 + 2π * k / m))) for k in 0:m-1]
    wall, window, gap = Mesh(), Mesh(), Mesh()
    # window disk
    K = max(1, round(Int, rw / hcell))
    counts = fill(n, K)
    for k in K-1:-1:1
        m = counts[k+1]
        counts[k] = iseven(m) && m ÷ 2 >= 6 && 2π * (k * rw / K) / (m ÷ 2) <= 1.25hcell ? m ÷ 2 : m
    end
    loops = [ring(k == K ? equal_area_radius(rw, n) : k * rw / K, counts[k]) for k in 1:K]
    center = xyz((0.0, 0.0))
    for i in 1:counts[1]
        tri!(window, center, loops[1][i], loops[1][mod1(i + 1, counts[1])])
    end
    for k in 2:K
        inn, out = loops[k-1], loops[k]
        mi, mo = length(inn), length(out)
        for i in 1:mi
            a, b = inn[i], inn[mod1(i + 1, mi)]
            if mo == mi
                quad!(window, a, b, out[mod1(i + 1, mo)], out[i])
            else                       # two outer sides per inner one: three triangles
                p, q, r = out[2i-1], out[2i], out[mod1(2i + 1, mo)]
                tri!(window, a, p, q); tri!(window, a, q, b); tri!(window, b, q, r)
            end
        end
    end
    # gap ring
    outer = loops[K]
    if rg > rw
        gl = ring(equal_area_radius(rg, n), n)
        for i in 1:n
            quad!(gap, outer[i], outer[mod1(i + 1, n)], gl[mod1(i + 1, n)], gl[i])
        end
        outer = gl
    end
    # wall between the outer polygon and the box: advance along whichever loop gives
    # the shorter new edge (both loops are convex about the center and start at the
    # same angle)
    M = length(bnd)
    i = j = 1
    while i <= M || j <= n
        if j > n || (i <= M && norm(bnd[mod1(i + 1, M)] - outer[mod1(j, n)]) <=
                               norm(bnd[mod1(i, M)] - outer[mod1(j + 1, n)]))
            tri!(wall, bnd[i], bnd[mod1(i + 1, M)], outer[mod1(j, n)])
            i += 1
        else
            tri!(wall, bnd[mod1(i, M)], outer[j], outer[mod1(j + 1, n)])
            j += 1
        end
    end
    return wall, window, gap
end

"""
    HohlraumProfile(z, r)

Wall profile of a hohlraum of revolution about `z` (`hohlraum(; profile)`): break points
at heights `z` (mm, nondecreasing, from the bottom end `z[1] = −L/2` to the top end
`z[end] = L/2`) with wall radii `r` (mm, positive). Between consecutive breaks the wall is
a section: a cylinder where the two radii are equal, a conical frustum where they differ,
and a flat ring (a step in radius) where two breaks share a height. A frustraum is two
frusta meeting at the waist ([`frustraum_profile`](@ref)); a multi-section hohlraum is a
radius per break.
"""
struct HohlraumProfile
    z::Vector{Float64}
    r::Vector{Float64}
    function HohlraumProfile(z::AbstractVector{<:Real}, r::AbstractVector{<:Real})
        length(z) == length(r) || throw(ArgumentError("a profile needs one radius per break"))
        length(z) >= 2 || throw(ArgumentError("a profile needs at least two breaks"))
        all(>(0), r) || throw(ArgumentError("a profile's radii must be positive"))
        issorted(z) || throw(ArgumentError("a profile's breaks must go up in z"))
        z[end] > z[1] || throw(ArgumentError("a profile must have a length"))
        abs(z[1] + z[end]) <= 1e-12 * (z[end] - z[1]) ||
            throw(ArgumentError("a profile is centered on its waist: its ends must be at ±length/2"))
        for i in 1:length(z)-1
            z[i] == z[i+1] && r[i] == r[i+1] && throw(ArgumentError("breaks $i and $(i + 1) coincide"))
        end
        (z[1] == z[2] || z[end-1] == z[end]) &&
            throw(ArgumentError("a profile's end is its cap, not a step: start and end with a section"))
        return new(Float64.(z), Float64.(r))
    end
end

"""
    frustraum_profile(; length, waist_radius, end_radius, bottom_end_radius = end_radius)

[`HohlraumProfile`](@ref) of a frustraum: two conical frusta joined at the waist `z = 0`
(radius `waist_radius`, mm), narrowing to `end_radius` at the top end (`z = length/2`)
and `bottom_end_radius` at the bottom. The wall angle from the axis is
`atand((waist_radius − end_radius) / (length/2))`.
"""
function frustraum_profile(; length::Real, waist_radius::Real, end_radius::Real,
                           bottom_end_radius::Real = end_radius)
    length > 0 || throw(ArgumentError("a frustraum needs a positive length"))
    return HohlraumProfile([-length / 2, 0.0, length / 2], [bottom_end_radius, waist_radius, end_radius])
end

# The sections of a profile: (r0, z0, r1, z1) per section, bottom to top.
_sections(p::HohlraumProfile) = [(p.r[i], p.z[i], p.r[i+1], p.z[i+1]) for i in 1:length(p.z)-1]

# A section's kind: :cylinder, :frustum or :ring (a step).
_section_kind(r0, z0, r1, z1) = z0 == z1 ? :ring : r0 == r1 ? :cylinder : :frustum

"""
    hohlraum_counts(; length, radius, leh_radius, h) -> (; nz, nphi, nr_cap, nr_leh)

Mesh counts of `hohlraum` for wall and end-cap elements no larger than `h` (mm) along
either side: `nz` axial and `nphi` azimuthal (a multiple of 8, so the NIF and OMEGA
port symmetries map the mesh onto itself) wall cells, `nr_cap` radial cap rings and
`nr_leh` LEH rings (their sides `≤ h`).
"""
function hohlraum_counts(; length, radius, leh_radius, h)
    h > 0 || throw(ArgumentError("the element size h must be positive"))
    nz = max(2, ceil(Int, length / h))
    nphi = 8 * max(1, ceil(Int, 2π * radius / h / 8))
    nr_cap = max(1, ceil(Int, (radius - leh_radius) / h))
    nr_leh = max(1, ceil(Int, leh_radius / h))
    return (; nz, nphi, nr_cap, nr_leh)
end

# The wall, cap, LEH and fillet counts `hohlraum` meshes with: from the element size `h`
# when given (the straight wall and the flat caps are shorter by the fillets), and the
# fillet bands about one wall cell of arc each unless `nr_corner` is given.
function _hohlraum_grid(; length, radius, leh_radius, h, rc, nz, nphi, nr_cap, nr_leh, nr_corner)
    if h !== nothing
        (; nz, nphi, nr_cap, nr_leh) = hohlraum_counts(; length, radius, leh_radius, h)
        nz = max(2, ceil(Int, (length - 2rc) / h))
        nr_cap = max(1, ceil(Int, (radius - rc - leh_radius) / h))
    end
    hcell = min(2π * radius / nphi, (length - 2rc) / nz)
    nr_corner = rc > 0 ? something(nr_corner, max(2, ceil(Int, (π / 2) * rc / hcell))) : nr_corner
    return (; nz, nphi, nr_cap, nr_leh, nr_corner)
end

"""
    Sombrero(; inner_radius = 4.0, outer_radius = 9.2, offset = 0.623, material = BLACK,
             nr = 4, nphi = 96)

Entrance-hole shield ("sombrero") for `hohlraum(; sombrero)`: a flat annular sheet
around each end of the hohlraum, from `inner_radius` to `outer_radius` (mm), its
mid-plane `offset` mm from the entrance-hole plane (`length/2`) toward the waist
(negative: beyond the hole). `hohlraum` meshes it as an `annulus_mesh` of `nr × nphi`
patches (equal-area rings; the polygon edges lie within `r (1 − cos(π/nphi))` of the
circles, 5 µm at 9.2 mm with 96 sides) carrying its exact `DiskSurface`.

The defaults are representative of a NIF symcap shield, relative to the entrance-hole
washer's inner face: bore 4.0 mm and rim 9.2 mm radius, a 39 µm sheet whose mid-plane
lies 0.623 mm inside that face. Both shields are identical (mirror images). A real sheet
is dimpled, within about ±0.12 mm of its mean plane, and runs past the rim where it is
mounted; the model keeps the flat mean plane.

`material` is an engineering default: `BLACK`, absorbing laser
light and emitting nothing. The shield lies outside the cavity and does not see into it
(beside the wall, the wall hides the hole from it), so its albedo and conversion do not
change the hohlraum's radiation.
"""
struct Sombrero
    inner_radius::Float64
    outer_radius::Float64
    offset::Float64
    material::Material
    nr::Int
    nphi::Int
    function Sombrero(inner_radius, outer_radius, offset, material, nr, nphi)
        0 < inner_radius < outer_radius ||
            throw(ArgumentError("need 0 < inner_radius < outer_radius for a sombrero"))
        nr >= 1 && nphi >= 3 || throw(ArgumentError("a sombrero needs nr ≥ 1 and nphi ≥ 3"))
        new(inner_radius, outer_radius, offset, material, nr, nphi)
    end
end
Sombrero(; inner_radius::Real = 4.0, outer_radius::Real = 9.2, offset::Real = 0.623,
         material::Material = BLACK, nr::Integer = 4, nphi::Integer = 96) =
    Sombrero(inner_radius, outer_radius, offset, material, nr, nphi)

# The two shields of `s` for a hohlraum of half-length `h`, wall radius `R` and entrance
# holes of radius `rl`: one-sided, facing away from the waist.
function _sombreros(s::Sombrero, h, R, rl)
    z = h - s.offset
    z > 0 || throw(ArgumentError("the sombrero offset puts the shields past the waist"))
    # beside the wall the bore must clear it; beyond the hole, the hole
    s.inner_radius > (s.offset >= 0 ? R : rl) ||
        throw(ArgumentError("the sombrero bore (inner_radius) cuts into the hohlraum"))
    return [Component(name, annulus_mesh(s.inner_radius, s.outer_radius, sgn * z; nr = s.nr,
                                         nphi = s.nphi, up = sgn > 0), s.material;
                      surface = DiskSurface(center = (0, 0, sgn * z)))
            for (sgn, name) in ((1, :sombrero_top), (-1, :sombrero_bottom))]
end

"""
    capsule_counts(capsule_radius, capsule_h) -> (; capsule_ntheta, capsule_nphi)

Capsule mesh counts of `hohlraum` for elements no larger than `capsule_h` (mm) along a
meridian or the equator (`capsule_nphi` a multiple of 4).
"""
function capsule_counts(r, ch)
    ch > 0 || throw(ArgumentError("the capsule element size must be positive"))
    nt = max(2, ceil(Int, π * r / ch))
    return (; capsule_ntheta = nt, capsule_nphi = 4 * max(1, ceil(Int, 2π * r / ch / 4)))
end

"""
    hohlraum(; length = 9.4, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9,
             nz = 36, nphi = 48, nr_cap = 6, nr_leh = 4, capsule_ntheta = 24, capsule_nphi = 32,
             h = nothing, capsule_h = nothing, capsule_mesh = :compact, wall = gold(),
             capsule = ablator(), capsule_center = (0,0,0), windows = (), corner_radius = 0,
             nr_corner = nothing, leh_closure = nothing, leh_open_radius = nothing,
             nr_lip = nothing, sombrero = nothing, top_end = nothing, bottom_end = nothing,
             profile = nothing)

Cylindrical hohlraum along `z` (or, with `profile`, a hohlraum of revolution of that
profile), centered on TCC (dimensions in mm), with components

- `:wall` — inner cylinder wall (radiates inward)
- `:cap_top`, `:cap_bottom` — end-cap washers between the LEH and the wall
- `:leh_top`, `:leh_bottom` — laser entrance holes (apertures; their incident flux is
  the x-ray loss through the LEH)
- `:leh_top_lip`, `:leh_bottom_lip` — with `leh_closure`, each hole's closing lip (after
  its hole)
- `:capsule` — spherical capsule (omit with `capsule_radius = 0`)
- `:window_1`, … — diagnostic windows in the wall, one component per `WallWindow` in
  `windows` (or the window's own name): closed windows are opaque wall patches with
  their own material, open ones apertures. A window with a gap is followed by its gap
  component `:window_1_gap`, ….
- `:sombrero_top`, `:sombrero_bottom` — with `sombrero`, the entrance-hole shields
  (last)

Laser light that strikes the outside of the wall or caps (clipping on the LEH edge)
is reported as the `_outside` back-side loss of those components; light that leaves
through an open window is `:missed`. The wall, caps, LEHs, windows and capsule carry
their exact `CylinderSurface`, `DiskSurface` and `SphereSurface`, so the radiation
model uses the true geometry.
The target is declared `closed`, so `RadiationModel` balances its view factors
(`enforce_closure!`) unless given `closed = false`: windows, their gaps and closing lips
are apertures, the sombreros see nothing of the cavity, and with open windows and gaps
the balanced loss and capsule drive agree with a continuum Monte Carlo within 0.1 %
(research notes, "Closure on bored, windowed and CAD cavities").

Each end (`top_end` at `+length/2`, `bottom_end` at `−length/2`) is, by default, an
entrance hole of radius `leh_radius`; given a number, an entrance hole of that radius
(its cap's and hole's ring counts from `h` when given, else `nr_cap` and `nr_leh`; with
`leh_closure`, its open radius `leh_open_radius` scaled by the ratio of the radii);
`:closed`, a cap with no hole (the hole's mesh, and its lip's, merged into the cap, of
the wall's material); `:open`, neither cap nor hole: the cavity is open there, to a
physics package placed against it or to nothing (a halfraum). A target with an open end
is not declared `closed` (its view factors are not balanced: what leaves through the
open end and strikes nothing is lost), and has no fillet (`corner_radius = 0`) or
sombrero. The defaults leave the hohlraum exactly as without these keywords.

With windows, the wall and end-cap meshes contain the windows' edges as grid lines
(`conforming_grid`; the wall's axial grid gets their `z` edges, wall and caps share
the azimuthal grid with their `φ` edges), so each rectangular window (and gap frame) is
tiled exactly by whole patches, and the caps still meet the wall vertex for vertex. A
circular window is cut the same way as a square box around it and its gap, half a wall
cell wider, which is re-meshed to conform to the circles (`_circular_window`) and
stitched to the surrounding grid at the box's own vertices. Without windows the mesh
is the uniform `nz × nphi` grid.

Mesh size: give `h` (mm) to set the wall, cap and LEH counts from an element size
instead (`hohlraum_counts`), and `capsule_h` for the capsule (`capsule_counts`); they
replace the corresponding counts.

`corner_radius` (mm) rounds the wall–cap corners into quarter-circle fillets of that
radius (a hohlraum machined with curved corners): the wall stops `corner_radius` short
of each end and the caps `corner_radius` short of the wall radius. Each fillet is
`nr_corner` bands (default: about one wall cell of arc each, at least 2), components
`:corner_top_1`, … (from the wall to the cap) and `:corner_bottom_1`, …, of the wall
material. Each band is the exact frustum through its two rim circles (`ConeSurface`),
so the fillet is the chord polygon of the arc (sagitta `corner_radius (1 − cos(π/4k))`,
20 µm for 1 mm in 4 bands) and the bands, wall and caps meet vertex for vertex.

`leh_closure`, a `GapClosure`, closes each entrance hole to x-rays from `leh_radius` down
to `leh_open_radius` during the pulse: the annulus between them is a closing lip
(`aperture(...; closure)`), transparent to laser light and lines of sight, which lets
its incident flux escape until `closure.closure` (ns) and across the next
`closure.width` ns turns into plasma that re-emits with the albedo of `closure.after`
(the wall's, e.g. `WallLossAlbedo()`). Engineering default for low-fill NIF hohlraums:
open through the foot, closing over ~0.5 ns at the rise to peak, to
`leh_open_radius = 0.84 leh_radius` (Schneider et al. 2015: ~0.84 r₀ through the peak;
Young et al. 2020 use 85 % static radii; research notes, "Closing entrance holes"). Beams
pass the lip unclipped: the laser keeps its own path open (Chen et al. 2020).

`leh_closure = :computed` (or a `LipClosure`) computes the closure instead: each lip's
plasma edge moves in from `leh_radius` at the sound speed of the gold blown off the lip
and stands a fixed depth inside the rim while the beams pass at high power
(`LipClosure`, calibrated on measured LEH radii; `drive_history` integrates it per hole
from the lip's own temperature and the beams entering that hole). The lip then reaches
down to `leh_open_radius` (default 0.6 `leh_radius`, the floor of the closure) in
`nr_lip` rings (default 6; 1 for a `GapClosure`).

`sombrero` (`true` for the measured `Sombrero()`, or a `Sombrero`) adds an annular
shield around each end, one-sided and facing away from the waist: opaque to beams
(clipping, deposition: power on its waist side is booked as `sombrero_top_outside`) and
lines of sight, hardware for `clearance`. It lies outside the cavity: its elements see
none of the hohlraum's (their view-factor rows and columns are zero), so the enclosure
stays closed, `closed = true` balancing leaves them alone, and the hohlraum's solution
does not change. Off by default.

`profile`, a [`HohlraumProfile`](@ref), gives the wall a radius per break instead of
`length` and `radius` (give neither): a frustraum ([`frustraum_profile`](@ref)) or a
multi-section hohlraum. Each section between two breaks is its own wall component,
`:wall_1`, `:wall_2`, … from the bottom, each with `part = :wall` (so `:wall` names the
whole wall in `elements`, `drive_history(; wall)` and every name-taking analysis) and its
exact surface: a `CylinderSurface`, a `ConeSurface` (a frustum, meshed as planar
trapezoids between parallel chords) or a `DiskSurface` (a flat ring where the radius
steps). All sections share the azimuthal grid, so they meet vertex for vertex, and each
end's cap runs from its hole to that end's radius. With `h` the azimuthal count comes
from the widest radius and each section's axial count from its meridian length;
otherwise `nz` is shared out among the sections by length. Windows are cut into
cylindrical sections only (each window within one); a window on a frustum or across a
break is refused, as are corner fillets. The ends, closing lips, sombreros and the
capsule work as for the cylinder; a sombrero's bore must clear the widest radius. A
profile of one cylindrical section is the cylindrical hohlraum exactly (`:wall` and all).

`capsule_mesh` picks the capsule mesh: `:compact` (the default), a geodesic
`icosphere_mesh` of nearly equilateral triangles with sides of about `capsule_h`
(`icosphere_frequency`; without `capsule_h`, the lat-long meridional size
`π capsule_radius / capsule_ntheta`), or `:latlong` (`sphere_mesh`, `capsule_ntheta ×
capsule_nphi`). The capsule mesh sets the flux's P4: lat-long cells are slivers near
the poles, and the compact mesh reaches the same P4 error with several times fewer
elements (research notes, "Compact capsule mesh"). `:latlong` is the comparison mesh of
that convergence study (`examples/capsule_mesh.jl`, `test/capsule_mesh.jl`).
"""
function hohlraum(; length::Union{Nothing,Real} = nothing, radius::Union{Nothing,Real} = nothing,
                  leh_radius = 1.55, capsule_radius = 0.9,
                  nz::Integer = 36, nphi::Integer = 48, nr_cap::Integer = 6, nr_leh::Integer = 4,
                  capsule_ntheta::Integer = 24, capsule_nphi::Integer = 32,
                  h = nothing, capsule_h = nothing, capsule_mesh::Symbol = :compact,
                  wall::Material = gold(), capsule::Material = ablator(),
                  capsule_center = (0.0, 0.0, 0.0), windows = (), corner_radius::Real = 0.0,
                  nr_corner::Union{Nothing,Integer} = nothing,
                  leh_closure::Union{Nothing,Symbol,GapClosure,LipClosure} = nothing,
                  leh_open_radius::Union{Nothing,Real} = nothing, nr_lip::Union{Nothing,Integer} = nothing,
                  sombrero::Union{Nothing,Bool,Sombrero} = nothing,
                  top_end::Union{Nothing,Real,Symbol} = nothing,
                  bottom_end::Union{Nothing,Real,Symbol} = nothing,
                  profile::Union{Nothing,HohlraumProfile} = nothing)
    if profile !== nothing
        length === nothing && radius === nothing ||
            throw(ArgumentError("a profile sets the length and the radii: give no length or radius with it"))
        secs = _sections(profile)
        if Base.length(secs) == 1 && _section_kind(secs[1]...) === :cylinder
            # one cylindrical section: the cylindrical hohlraum itself
            length, radius = profile.z[end] - profile.z[1], profile.r[1]
            profile = nothing
        end
    end
    length = something(length, 9.4)
    radius = something(radius, 2.875)
    if leh_closure isa Symbol
        leh_closure === :computed || throw(ArgumentError("leh_closure must be :computed, a GapClosure or a LipClosure"))
        leh_closure = LipClosure()
    end
    computed = leh_closure isa LipClosure
    leh_open_radius = something(leh_open_radius, (computed ? 0.6 : 0.84) * leh_radius)
    nr_lip = something(nr_lip, computed ? 6 : 1)
    leh_closure === nothing || 0 < leh_open_radius < leh_radius ||
        throw(ArgumentError("need 0 < leh_open_radius < leh_radius"))
    nr_lip >= 1 || throw(ArgumentError("nr_lip must be at least 1"))
    capsule_mesh in (:latlong, :compact) ||
        throw(ArgumentError("capsule_mesh must be :latlong or :compact, got :$capsule_mesh"))
    windows = windows isa WallWindow ? [windows] : collect(WallWindow, windows)
    profile === nothing ||
        return _profile_hohlraum(profile; leh_radius, capsule_radius, nz, nphi, nr_cap, nr_leh,
                                 capsule_ntheta, capsule_nphi, h, capsule_h, capsule_mesh, wall, capsule,
                                 capsule_center, windows, corner_radius, leh_closure, leh_open_radius,
                                 nr_lip, sombrero, top_end, bottom_end)
    0 < leh_radius < radius || throw(ArgumentError("need 0 < leh_radius < radius"))
    rc = Float64(corner_radius)
    0 <= rc < min(radius - leh_radius, length / 2) ||
        throw(ArgumentError("need 0 ≤ corner_radius < min(radius - leh_radius, length / 2)"))
    capsule_radius < min(radius, length / 2) || throw(ArgumentError("capsule does not fit"))
    (; nz, nphi, nr_cap, nr_leh, nr_corner) =
        _hohlraum_grid(; length, radius, leh_radius, h, rc, nz, nphi, nr_cap, nr_leh, nr_corner)
    ends = (_hohlraum_end(top_end, :top_end; leh_radius, radius, rc, h, nr_cap, nr_leh),
            _hohlraum_end(bottom_end, :bottom_end; leh_radius, radius, rc, h, nr_cap, nr_leh))
    opened = _check_open_ends(ends, rc, sombrero)
    h = length / 2
    hw = h - rc                                     # half-length of the straight wall
    hcell = min(2π * radius / nphi, 2hw / nz)
    exts = _check_windows(windows, radius, hw, hcell)
    zbreaks = [z for e in exts for (_, _, za, zb) in e for z in (za, zb)]
    phibreaks = [φ for e in exts for (φa, Δφ, _, _) in e for φ in (φa, φa + Δφ)]
    wall_mesh = cylinder_mesh(radius, -hw, hw; nz, nphi, inward = true, zbreaks, phibreaks)
    cyl = CylinderSurface(radius)
    wall_mesh, window_comps = _cut_windows(wall_mesh, windows, exts, eachindex(windows), radius, hcell, cyl, wall)
    comps = Component[Component(:wall, wall_mesh, wall; surface = cyl)]
    _end_components!(comps, ends, h, (radius - rc, radius - rc), nphi, phibreaks, wall;
                     leh_closure, leh_radius, leh_open_radius, nr_lip)
    if rc > 0
        k = nr_corner
        k >= 1 || throw(ArgumentError("nr_corner must be at least 1"))
        for (sgn, name) in ((1, :corner_top), (-1, :corner_bottom)), i in 1:k
            push!(comps, _corner_band(Symbol(name, :_, i), radius, hw, rc, sgn, i, k, nphi, phibreaks, wall))
        end
    end
    _capsule_component!(comps, capsule_radius, capsule_mesh, capsule_ntheta, capsule_nphi, capsule_h,
                        capsule, capsule_center)
    append!(comps, window_comps)
    _sombrero_components!(comps, sombrero, h, radius, ends)
    return Target(comps; closed = !opened)
end

# Whether either end is open, refusing what an open end cannot have (a fillet, a sombrero).
function _check_open_ends(ends, rc, sombrero)
    opened = any(e -> e.kind === :open, ends)
    opened && rc > 0 &&
        throw(ArgumentError("an open end has no corner to round: give corner_radius = 0"))
    opened && sombrero !== nothing && sombrero !== false &&
        throw(ArgumentError("a sombrero shields entrance holes: an open end has none"))
    return opened
end

# Cut `windows` (with extents `exts`, numbered `ks`) out of the cylindrical wall mesh of
# radius `R` whose grid contains their edges: the rest of the wall and the windows'
# components (each followed by its gap), on the cylinder `cyl`.
function _cut_windows(wall_mesh, windows, exts, ks, R, hcell, cyl, wall)
    part(name, mesh, open, mat) = open ? aperture(name, mesh; surface = cyl) :
                                  Component(name, mesh, mat === nothing ? wall : mat; surface = cyl)
    window_comps = Component[]
    for (k, w, e) in zip(ks, windows, exts)
        # the mesh conforms to the extents' edges, so patch centroids are well inside or
        # well outside (by half a patch)
        if w.circular
            box, wall_mesh = split_mesh(wall_mesh, c -> _in_window(e[1], c))
            rest, sel, gsel = _circular_window(box, R, deg2rad(w.phi), w.z, w.width / 2,
                                               w.width / 2 + w.gap, _window_segments(w, hcell), hcell)
            wall_mesh = merge_meshes(wall_mesh, rest)
        else
            sel, wall_mesh = split_mesh(wall_mesh, c -> _in_window(e[1], c))
            gsel = nothing
            if w.gap > 0
                gsel, wall_mesh = split_mesh(wall_mesh, c -> _in_window(e[2], c))
            end
        end
        name = w.name === nothing ? Symbol(:window_, k) : w.name
        push!(window_comps, part(name, sel, w.open, w.material))
        w.gap > 0 && push!(window_comps, part(Symbol(name, :_gap), gsel, w.gap_material === nothing,
                                              w.gap_material))
    end
    return wall_mesh, window_comps
end

# The end caps (washers from each end's hole to `rcaps`, top then bottom), the entrance
# holes and their closing lips, appended to `comps`; a closed end's hole (and lip) merged
# into its cap, of the wall's material; an open end with neither cap nor hole.
function _end_components!(comps, ends, h, rcaps, nphi, phibreaks, wall; leh_closure, leh_radius,
                          leh_open_radius, nr_lip)
    push!(comps,
          Component(:cap_top, annulus_mesh(ends[1].radius, rcaps[1], h; nr = ends[1].nr_cap, nphi,
                                           up = false, phibreaks), wall;
                    surface = DiskSurface(center = (0, 0, h))),
          Component(:cap_bottom, annulus_mesh(ends[2].radius, rcaps[2], -h; nr = ends[2].nr_cap, nphi,
                                              up = true, phibreaks), wall;
                    surface = DiskSurface(center = (0, 0, -h))))
    for (sgn, name, e) in ((1, :leh_top, ends[1]), (-1, :leh_bottom, ends[2]))
        e.kind === :open && continue
        disk = DiskSurface(center = (0, 0, sgn * h))
        rl = e.radius
        # the hole takes the cap's azimuths, so its rim meets the cap's inner edge vertex
        # for vertex (half the count left a crack at the rim that rays leak through)
        if leh_closure === nothing
            push!(comps, aperture(name, annulus_mesh(0.0, rl, sgn * h; nr = e.nr_leh, nphi,
                                                     up = sgn < 0, phibreaks); surface = disk))
        else
            # the hole open at the peak, and the lip that closes around it
            ro = rl == leh_radius ? leh_open_radius : leh_open_radius * rl / leh_radius
            push!(comps, aperture(name, annulus_mesh(0.0, ro, sgn * h; nr = max(1, round(Int, e.nr_leh * ro / rl)),
                                                     nphi, up = sgn < 0, phibreaks); surface = disk))
            push!(comps, aperture(Symbol(name, :_lip), annulus_mesh(ro, rl, sgn * h; nr = nr_lip, nphi,
                                                                    up = sgn < 0, phibreaks);
                                  surface = disk, closure = leh_closure))
        end
    end
    for (cap, leh, e) in ((:cap_top, :leh_top, ends[1]), (:cap_bottom, :leh_bottom, ends[2]))
        e.kind === :leh && continue
        k = findfirst(c -> c.name === cap, comps)
        if e.kind === :closed
            hole = filter(c -> c.name in (leh, Symbol(leh, :_lip)), comps)
            c = comps[k]
            comps[k] = Component(cap, merge_meshes(c.mesh, (d.mesh for d in hole)...), c.material;
                                 surface = c.surface)
            filter!(c -> !(c in hole), comps)
        else
            deleteat!(comps, k)
        end
    end
    return comps
end

# The capsule (none for `capsule_radius = 0`), appended to `comps`.
function _capsule_component!(comps, capsule_radius, capsule_mesh, capsule_ntheta, capsule_nphi, capsule_h,
                             capsule, capsule_center)
    capsule_radius > 0 || return comps
    if capsule_mesh === :compact
        capsule_n = icosphere_frequency(capsule_radius, something(capsule_h, π * capsule_radius / capsule_ntheta))
        cm = icosphere_mesh(capsule_radius; n = capsule_n)
    else
        if capsule_h !== nothing
            (; capsule_ntheta, capsule_nphi) = capsule_counts(capsule_radius, capsule_h)
        end
        cm = sphere_mesh(capsule_radius; ntheta = capsule_ntheta, nphi = capsule_nphi)
    end
    cm = translate(cm, Vec3(capsule_center...))
    return push!(comps, Component(:capsule, cm, capsule;
                                  surface = SphereSurface(capsule_radius; center = capsule_center)))
end

# The entrance-hole shields (last), for a hohlraum of half-length `h` whose wall's widest
# radius near the ends is `R`.
function _sombrero_components!(comps, sombrero, h, R, ends)
    sombrero === true && (sombrero = Sombrero())
    sombrero isa Sombrero && append!(comps, _sombreros(sombrero, h, R, max(ends[1].radius, ends[2].radius)))
    return comps
end

# `hohlraum(; profile)` with more than one section, or a section that is not a cylinder.
function _profile_hohlraum(prof::HohlraumProfile; leh_radius, capsule_radius, nz, nphi, nr_cap, nr_leh,
                           capsule_ntheta, capsule_nphi, h, capsule_h, capsule_mesh, wall, capsule,
                           capsule_center, windows, corner_radius, leh_closure, leh_open_radius, nr_lip,
                           sombrero, top_end, bottom_end)
    corner_radius == 0 ||
        throw(ArgumentError("a profile's corners are its breaks: give corner_radius = 0"))
    secs = _sections(prof)
    L = prof.z[end] - prof.z[1]
    hz = prof.z[end]                                # the ends are at ±hz
    rends = (prof.r[end], prof.r[1])                # top, bottom
    rmax = maximum(prof.r)
    0 < leh_radius || throw(ArgumentError("need 0 < leh_radius"))
    # the capsule inside every section: its distance in its meridian half-plane from each
    # section's line, and from the ends
    if capsule_radius > 0
        c = Vec3(capsule_center...)
        q = (hypot(c[1], c[2]), c[3])
        abs(q[2]) + capsule_radius < hz && all(s -> _segment_distance(q, s) > capsule_radius, secs) &&
            _profile_radius(prof, q[2]) > q[1] || throw(ArgumentError("capsule does not fit"))
    end
    # counts: azimuthal from the widest radius, axial per section from its meridian length
    lens = [hypot(r1 - r0, z1 - z0) for (r0, z0, r1, z1) in secs]
    if h !== nothing
        h > 0 || throw(ArgumentError("the element size h must be positive"))
        nphi = 8 * max(1, ceil(Int, 2π * rmax / h / 8))
        ns = [max(1, ceil(Int, l / h)) for l in lens]
    else
        ns = [max(1, round(Int, nz * l / sum(lens))) for l in lens]
    end
    endcounts = map(rends) do re
        h === nothing ? (; nr_cap, nr_leh) :
        (; nr_cap = max(1, ceil(Int, (re - leh_radius) / h)), nr_leh = max(1, ceil(Int, leh_radius / h)))
    end
    ends = (_hohlraum_end(top_end, :top_end; leh_radius, radius = rends[1], rc = 0.0, h, endcounts[1]...),
            _hohlraum_end(bottom_end, :bottom_end; leh_radius, radius = rends[2], rc = 0.0, h, endcounts[2]...))
    for (e, re, key) in zip(ends, rends, (:top_end, :bottom_end))
        e.kind === :open || e.radius < re ||
            throw(ArgumentError("need $key's entrance hole (radius $(e.radius)) narrower than the end ($re)"))
    end
    opened = _check_open_ends(ends, 0.0, sombrero)
    # windows: each within one cylindrical section, cut as in a cylindrical wall
    home = zeros(Int, Base.length(windows))
    exts = Vector{Any}(undef, Base.length(windows))
    for (k, w) in enumerate(windows)
        for (i, (r0, z0, r1, z1)) in enumerate(secs)
            _section_kind(r0, z0, r1, z1) === :cylinder || continue
            e = _window_extents(w, r0, min(2π * r0 / nphi, (z1 - z0) / ns[i]))
            tol = 1e-9 * max(r0, hz)
            z0 - tol <= last(e)[3] && last(e)[4] <= z1 + tol && (home[k] = i; break)
        end
        home[k] == 0 && throw(ArgumentError("window $k does not lie within one cylindrical section: " *
                                            "windows are cut into cylindrical walls only"))
    end
    hcells = [min(2π * r0 / nphi, (z1 - z0) / max(ns[i], 1)) for (i, (r0, z0, r1, z1)) in enumerate(secs)]
    for i in unique(home)
        ks = findall(==(i), home)
        r0, z0, _, z1 = secs[i]
        # `_check_windows` takes a wall centered on z = 0: shift the windows to it
        zc = (z0 + z1) / 2
        shifted = [WallWindow(w.phi, w.z - zc, w.width, w.height, w.circular, w.segments, w.material, w.open,
                              w.name, w.gap, w.gap_material) for w in windows[ks]]
        ex = _check_windows(shifted, r0, (z1 - z0) / 2, hcells[i])
        for (k, e) in zip(ks, ex)
            exts[k] = [(φa, Δφ, za + zc, zb + zc) for (φa, Δφ, za, zb) in e]
        end
    end
    phibreaks = [φ for e in exts for (φa, Δφ, _, _) in e for φ in (φa, φa + Δφ)]
    φs = isempty(phibreaks) ? range(0.0, 2π; length = nphi + 1) :
         conforming_grid(0.0, 2π, nphi, phibreaks; periodic = true)
    comps = Component[]
    window_comps = Component[]
    for (i, (r0, z0, r1, z1)) in enumerate(secs)
        name = Symbol(:wall_, i)
        kind = _section_kind(r0, z0, r1, z1)
        if kind === :cylinder
            ks = findall(==(i), home)
            zbreaks = [z for k in ks for (_, _, za, zb) in exts[k] for z in (za, zb)]
            m = cylinder_mesh(r0, z0, z1; nz = ns[i], nphi, inward = true, zbreaks, phibreaks)
            cyl = CylinderSurface(r0)
            m, wc = _cut_windows(m, windows[ks], exts[ks], ks, r0, hcells[i], cyl, wall)
            append!(window_comps, wc)
            push!(comps, Component(name, m, wall; surface = cyl, part = :wall))
        elseif kind === :frustum
            k = (r1 - r0) / (z1 - z0)
            surf = ConeSurface(atand(abs(k)); apex = (0.0, 0.0, z0 - r0 / k), axis = (0.0, 0.0, sign(k)))
            push!(comps, Component(name, _section_mesh(r0, z0, r1, z1, ns[i], φs), wall; surface = surf,
                                   part = :wall))
        else
            push!(comps, Component(name, _section_mesh(r0, z0, r1, z1, ns[i], φs), wall;
                                   surface = DiskSurface(center = (0, 0, z0)), part = :wall))
        end
    end
    _end_components!(comps, ends, hz, rends, nphi, phibreaks, wall;
                     leh_closure, leh_radius, leh_open_radius, nr_lip)
    _capsule_component!(comps, capsule_radius, capsule_mesh, capsule_ntheta, capsule_nphi, capsule_h,
                        capsule, capsule_center)
    append!(comps, window_comps)
    _sombrero_components!(comps, sombrero, hz, rmax, ends)
    return Target(comps; closed = !opened)
end

# A section of revolution from radius r0 at z0 to r1 at z1 (a frustum, a cylinder or a
# flat ring) in `n` bands on the azimuthal grid `φs`, facing the axis side of the profile
# (the cavity). Each patch is the trapezoid between two parallel chords at the same
# azimuths, so it is planar; the last ring is (r1, z1) exactly, so sections meet vertex
# for vertex.
function _section_mesh(r0, z0, r1, z1, n, φs)
    m = Mesh()
    ring(i) = i == n ? (Float64(r1), Float64(z1)) : (r0 + (r1 - r0) * i / n, z0 + (z1 - z0) * i / n)
    for i in 1:n, j in 1:length(φs)-1
        (ra, za), (rb, zb) = ring(i - 1), ring(i)
        a, b = φs[j], φs[j+1]
        p = [Vec3(ra * cos(a), ra * sin(a), za), Vec3(ra * cos(b), ra * sin(b), za),
             Vec3(rb * cos(b), rb * sin(b), zb), Vec3(rb * cos(a), rb * sin(a), zb)]
        φc = (a + b) / 2
        # the profile's normal toward the cavity: (−dz, dr) in the meridian half-plane
        add_polygon!(m, p; facing = Vec3(-(z1 - z0) * cos(φc), -(z1 - z0) * sin(φc), r1 - r0))
    end
    return m
end

# Distance in the meridian half-plane from the point q = (ρ, z) to the section's line.
function _segment_distance(q, (r0, z0, r1, z1))
    d = (r1 - r0, z1 - z0)
    t = clamp(((q[1] - r0) * d[1] + (q[2] - z0) * d[2]) / (d[1]^2 + d[2]^2), 0.0, 1.0)
    return hypot(q[1] - r0 - t * d[1], q[2] - z0 - t * d[2])
end

# The profile's wall radius at height z (the smallest where a ring steps).
function _profile_radius(p::HohlraumProfile, z)
    best = Inf
    for (r0, z0, r1, z1) in _sections(p)
        z0 <= z <= z1 || continue
        best = min(best, z0 == z1 ? min(r0, r1) : r0 + (r1 - r0) * (z - z0) / (z1 - z0))
    end
    return best
end

# One end of `hohlraum` (keyword `key`): `(; kind, radius, nr_cap, nr_leh)`, `kind` one
# of `:leh` (an entrance hole of `radius`), `:closed` (the cap covers the hole of
# `radius`) and `:open` (neither cap nor hole; `radius` and the counts unused). An end
# of the default radius keeps the hohlraum's counts; another radius takes counts from
# the element size `h` when given.
function _hohlraum_end(e, key::Symbol; leh_radius, radius, rc, h, nr_cap, nr_leh)
    e === nothing && return (; kind = :leh, radius = Float64(leh_radius), nr_cap, nr_leh)
    if e isa Symbol
        e in (:closed, :open) ||
            throw(ArgumentError("$key must be an LEH radius, :closed or :open, got :$e"))
        return (; kind = e, radius = Float64(leh_radius), nr_cap, nr_leh)
    end
    e isa Real || throw(ArgumentError("$key must be an LEH radius, :closed or :open"))
    0 < e < radius - rc || throw(ArgumentError("need 0 < $key < radius - corner_radius"))
    r = Float64(e)
    if r != leh_radius && h !== nothing
        nr_cap = max(1, ceil(Int, (radius - rc - r) / h))
        nr_leh = max(1, ceil(Int, r / h))
    end
    return (; kind = :leh, radius = r, nr_cap, nr_leh)
end

"""
    foil(name, center, normal; width, height = width, n = 20, material = gold(),
         two_sided = true, rotation = 0.0)

Rectangular foil (e.g. a backlighter) with front face along `normal`. `rotation`
(radians) spins it about its normal.
"""
function foil(name, center, normal; width, height = width, n::Integer = 20,
              material::Material = gold(conversion = 0.05), two_sided::Bool = true,
              rotation = 0.0)
    m = rectangle_mesh(width, height; nu = n, nv = max(1, round(Int, n * height / width)))
    Rz = SMatrix{3,3,Float64}(cos(rotation), sin(rotation), 0, -sin(rotation), cos(rotation), 0, 0, 0, 1)
    m = transform(m, rotation_to(normal) * Rz, Vec3(center...))
    return Component(name, m, material; two_sided)
end
