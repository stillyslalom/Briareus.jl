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
             nr_corner = nothing)

Cylindrical hohlraum along `z`, centered on TCC (dimensions in mm), with components

- `:wall` — inner cylinder wall (radiates inward)
- `:cap_top`, `:cap_bottom` — end-cap washers between the LEH and the wall
- `:leh_top`, `:leh_bottom` — laser entrance holes (apertures; their incident flux is
  the x-ray loss through the LEH)
- `:capsule` — spherical capsule (omit with `capsule_radius = 0`)
- `:window_1`, … — diagnostic windows in the wall, one component per `WallWindow` in
  `windows` (or the window's own name): closed windows are opaque wall patches with
  their own material, open ones apertures. A window with a gap is followed by its gap
  component `:window_1_gap`, ….

Laser light that strikes the outside of the wall or caps (clipping on the LEH edge)
is reported as the `_outside` back-side loss of those components; light that leaves
through an open window is `:missed`. The wall, caps, LEHs, windows and capsule carry
their exact `CylinderSurface`, `DiskSurface` and `SphereSurface`, so the radiation
model uses the true geometry.
The target is declared `closed`, so `RadiationModel` balances its view factors
(`enforce_closure!`) unless given `closed = false`.

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

`capsule_mesh` picks the capsule mesh: `:compact` (the default), a geodesic
`icosphere_mesh` of nearly equilateral triangles with sides of about `capsule_h`
(`icosphere_frequency`; without `capsule_h`, the lat-long meridional size
`π capsule_radius / capsule_ntheta`), or `:latlong` (`sphere_mesh`, `capsule_ntheta ×
capsule_nphi`). The capsule mesh sets the flux's P4: lat-long cells are slivers near
the poles, and the compact mesh reaches the same P4 error with several times fewer
elements (research notes, "Compact capsule mesh"). The default was `:latlong` until
2026-10-22; pass it to reproduce earlier baselines.
"""
function hohlraum(; length = 9.4, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9,
                  nz::Integer = 36, nphi::Integer = 48, nr_cap::Integer = 6, nr_leh::Integer = 4,
                  capsule_ntheta::Integer = 24, capsule_nphi::Integer = 32,
                  h = nothing, capsule_h = nothing, capsule_mesh::Symbol = :compact,
                  wall::Material = gold(), capsule::Material = ablator(),
                  capsule_center = (0.0, 0.0, 0.0), windows = (), corner_radius::Real = 0.0,
                  nr_corner::Union{Nothing,Integer} = nothing)
    0 < leh_radius < radius || throw(ArgumentError("need 0 < leh_radius < radius"))
    rc = Float64(corner_radius)
    0 <= rc < min(radius - leh_radius, length / 2) ||
        throw(ArgumentError("need 0 ≤ corner_radius < min(radius - leh_radius, length / 2)"))
    capsule_radius < min(radius, length / 2) || throw(ArgumentError("capsule does not fit"))
    capsule_mesh in (:latlong, :compact) ||
        throw(ArgumentError("capsule_mesh must be :latlong or :compact, got :$capsule_mesh"))
    if h !== nothing
        (; nz, nphi, nr_cap, nr_leh) = hohlraum_counts(; length, radius, leh_radius, h)
        # the straight wall and the flat caps are shorter by the fillets
        nz = max(2, ceil(Int, (length - 2rc) / h))
        nr_cap = max(1, ceil(Int, (radius - rc - leh_radius) / h))
    end
    if capsule_mesh === :compact && capsule_radius > 0
        capsule_n = icosphere_frequency(capsule_radius, something(capsule_h, π * capsule_radius / capsule_ntheta))
    elseif capsule_h !== nothing && capsule_radius > 0
        (; capsule_ntheta, capsule_nphi) = capsule_counts(capsule_radius, capsule_h)
    end
    h = length / 2
    hw = h - rc                                     # half-length of the straight wall
    windows = windows isa WallWindow ? [windows] : collect(WallWindow, windows)
    hcell = min(2π * radius / nphi, 2hw / nz)
    exts = _check_windows(windows, radius, hw, hcell)
    zbreaks = [z for e in exts for (_, _, za, zb) in e for z in (za, zb)]
    phibreaks = [φ for e in exts for (φa, Δφ, _, _) in e for φ in (φa, φa + Δφ)]
    wall_mesh = cylinder_mesh(radius, -hw, hw; nz, nphi, inward = true, zbreaks, phibreaks)
    cyl = CylinderSurface(radius)
    part(name, mesh, open, mat) = open ? aperture(name, mesh; surface = cyl) :
                                  Component(name, mesh, mat === nothing ? wall : mat; surface = cyl)
    window_comps = Component[]
    for (k, (w, e)) in enumerate(zip(windows, exts))
        # the mesh conforms to the extents' edges, so patch centroids are well inside or
        # well outside (by half a patch)
        if w.circular
            box, wall_mesh = split_mesh(wall_mesh, c -> _in_window(e[1], c))
            rest, sel, gsel = _circular_window(box, radius, deg2rad(w.phi), w.z, w.width / 2,
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
    comps = Component[
        Component(:wall, wall_mesh, wall; surface = cyl),
        Component(:cap_top, annulus_mesh(leh_radius, radius - rc, h; nr = nr_cap, nphi, up = false,
                                         phibreaks), wall;
                  surface = DiskSurface(center = (0, 0, h))),
        Component(:cap_bottom, annulus_mesh(leh_radius, radius - rc, -h; nr = nr_cap, nphi, up = true,
                                            phibreaks), wall;
                  surface = DiskSurface(center = (0, 0, -h))),
        aperture(:leh_top, disk_mesh(leh_radius, h; nr = nr_leh, nphi = nphi ÷ 2, up = false);
                 surface = DiskSurface(center = (0, 0, h))),
        aperture(:leh_bottom, disk_mesh(leh_radius, -h; nr = nr_leh, nphi = nphi ÷ 2, up = true);
                 surface = DiskSurface(center = (0, 0, -h))),
    ]
    if rc > 0
        k = something(nr_corner, max(2, ceil(Int, (π / 2) * rc / hcell)))
        k >= 1 || throw(ArgumentError("nr_corner must be at least 1"))
        for (sgn, name) in ((1, :corner_top), (-1, :corner_bottom)), i in 1:k
            push!(comps, _corner_band(Symbol(name, :_, i), radius, hw, rc, sgn, i, k, nphi, phibreaks, wall))
        end
    end
    if capsule_radius > 0
        cm = capsule_mesh === :compact ? icosphere_mesh(capsule_radius; n = capsule_n) :
             sphere_mesh(capsule_radius; ntheta = capsule_ntheta, nphi = capsule_nphi)
        cm = translate(cm, Vec3(capsule_center...))
        push!(comps, Component(:capsule, cm, capsule;
                               surface = SphereSurface(capsule_radius; center = capsule_center)))
    end
    append!(comps, window_comps)
    return Target(comps; closed = true)
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

"""
    disk_foil(name, center, normal; radius, nr = 8, nphi = 32, material, two_sided = true)

Circular foil with front face along `normal`.
"""
function disk_foil(name, center, normal; radius, nr::Integer = 8, nphi::Integer = 32,
                   material::Material = gold(conversion = 0.05), two_sided::Bool = true)
    m = transform(disk_mesh(radius, 0.0; nr, nphi), rotation_to(normal), Vec3(center...))
    return Component(name, m, material; two_sided)
end
