# Ready-made target geometries.

"""
    WallWindow(; phi, z = 0.0, width, height, material = nothing, open = false, name = nothing)

Rectangular diagnostic window cut into a hohlraum wall (see `hohlraum`), centered at
azimuth `phi` (degrees from `+x`) and axial position `z` (mm), `width` mm of arc
length along the wall by `height` mm along the axis.

- A closed window (`open = false`) is an opaque patch of wall that re-emits with its
  own `material` (albedo and laser conversion, e.g. a thin gold-lined or low-Z
  window); `material = nothing` keeps the wall's.
- An open window (`open = true`) is a hole, e.g. a window that has burned through: an
  aperture that laser light passes through and whose incident x-ray flux escapes.

`name` defaults to `:window_<k>` for the `k`-th window.
"""
struct WallWindow
    phi::Float64
    z::Float64
    width::Float64
    height::Float64
    material::Union{Nothing,Material}
    open::Bool
    name::Union{Nothing,Symbol}
    function WallWindow(phi, z, width, height, material, open, name)
        width > 0 && height > 0 || throw(ArgumentError("window width and height must be positive"))
        open && material !== nothing && throw(ArgumentError("an open window has no material"))
        new(phi, z, width, height, material, open, name === nothing ? nothing : Symbol(name))
    end
end
WallWindow(; phi, z = 0.0, width, height, material = nothing, open::Bool = false, name = nothing) =
    WallWindow(phi, z, width, height, material, open, name)

# Azimuthal (radians: start, angular width) and axial (mm) extent of a window on a
# wall of radius R.
_window_extent(w::WallWindow, R) =
    (deg2rad(w.phi) - w.width / 2R, w.width / R, w.z - w.height / 2, w.z + w.height / 2)

# Whether the point `c` (a patch centroid) lies inside the window extent.
function _in_window((φa, Δφ, za, zb), c)
    za < c[3] < zb || return false
    return mod(atan(c[2], c[1]) - φa, 2π) < Δφ
end

function _check_windows(windows, R, h)
    ext = [_window_extent(w, R) for w in windows]
    tol = 1e-9 * max(R, h)
    for (w, (_, Δφ, za, zb)) in zip(windows, ext)
        Δφ <= 2π + 1e-12 || throw(ArgumentError("window wider than the wall circumference"))
        -h - tol <= za && zb <= h + tol || throw(ArgumentError("window extends past the end of the wall"))
    end
    for i in eachindex(ext), j in i+1:lastindex(ext)
        (ai, wi, zai, zbi), (aj, wj, zaj, zbj) = ext[i], ext[j]
        zai < zbj - tol && zaj < zbi - tol || continue
        (mod(aj - ai, 2π) < wi - 1e-12 || mod(ai - aj, 2π) < wj - 1e-12) &&
            throw(ArgumentError("windows $i and $j overlap"))
    end
    return ext
end

"""
    hohlraum(; length = 9.4, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9,
             nz = 36, nphi = 48, nr_cap = 6, capsule_ntheta = 24, capsule_nphi = 32,
             wall = gold(), capsule = ablator(), capsule_center = (0,0,0), windows = ())

Cylindrical hohlraum along `z`, centered on TCC (dimensions in mm), with components

- `:wall` — inner cylinder wall (radiates inward)
- `:cap_top`, `:cap_bottom` — end-cap washers between the LEH and the wall
- `:leh_top`, `:leh_bottom` — laser entrance holes (apertures; their incident flux is
  the x-ray loss through the LEH)
- `:capsule` — spherical capsule (omit with `capsule_radius = 0`)
- `:window_1`, … — diagnostic windows in the wall, one component per `WallWindow` in
  `windows` (or the window's own name): closed windows are opaque wall patches with
  their own material, open ones apertures.

Laser light that strikes the outside of the wall or caps (clipping on the LEH edge)
is reported as the `_outside` back-side loss of those components; light that leaves
through an open window is `:missed`. The wall, caps, LEHs, windows and capsule carry
their exact `CylinderSurface`, `DiskSurface` and `SphereSurface`, so the radiation
model uses the true geometry.

With windows, the wall and end-cap meshes contain the windows' edges as grid lines
(`conforming_grid`; the wall's axial grid gets their `z` edges, wall and caps share
the azimuthal grid with their `φ` edges), so each window is tiled exactly by whole
patches, and the caps still meet the wall vertex for vertex. Without windows the mesh
is the uniform `nz × nphi` grid.
"""
function hohlraum(; length = 9.4, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9,
                  nz::Integer = 36, nphi::Integer = 48, nr_cap::Integer = 6, nr_leh::Integer = 4,
                  capsule_ntheta::Integer = 24, capsule_nphi::Integer = 32,
                  wall::Material = gold(), capsule::Material = ablator(),
                  capsule_center = (0.0, 0.0, 0.0), windows = ())
    0 < leh_radius < radius || throw(ArgumentError("need 0 < leh_radius < radius"))
    capsule_radius < min(radius, length / 2) || throw(ArgumentError("capsule does not fit"))
    h = length / 2
    windows = windows isa WallWindow ? [windows] : collect(WallWindow, windows)
    ext = _check_windows(windows, radius, h)
    zbreaks = [z for (_, _, za, zb) in ext for z in (za, zb)]
    phibreaks = [φ for (φa, Δφ, _, _) in ext for φ in (φa, φa + Δφ)]
    wall_mesh = cylinder_mesh(radius, -h, h; nz, nphi, inward = true, zbreaks, phibreaks)
    window_comps = Component[]
    for (k, (w, e)) in enumerate(zip(windows, ext))
        # the mesh conforms to the window edges, so patch centroids are well inside or
        # well outside (by half a patch)
        sel, wall_mesh = split_mesh(wall_mesh, c -> _in_window(e, c))
        name = w.name === nothing ? Symbol(:window_, k) : w.name
        push!(window_comps, w.open ? aperture(name, sel; surface = CylinderSurface(radius)) :
                            Component(name, sel, w.material === nothing ? wall : w.material;
                                      surface = CylinderSurface(radius)))
    end
    comps = Component[
        Component(:wall, wall_mesh, wall; surface = CylinderSurface(radius)),
        Component(:cap_top, annulus_mesh(leh_radius, radius, h; nr = nr_cap, nphi, up = false,
                                         phibreaks), wall;
                  surface = DiskSurface(center = (0, 0, h))),
        Component(:cap_bottom, annulus_mesh(leh_radius, radius, -h; nr = nr_cap, nphi, up = true,
                                            phibreaks), wall;
                  surface = DiskSurface(center = (0, 0, -h))),
        aperture(:leh_top, disk_mesh(leh_radius, h; nr = nr_leh, nphi = nphi ÷ 2, up = false);
                 surface = DiskSurface(center = (0, 0, h))),
        aperture(:leh_bottom, disk_mesh(leh_radius, -h; nr = nr_leh, nphi = nphi ÷ 2, up = true);
                 surface = DiskSurface(center = (0, 0, -h))),
    ]
    if capsule_radius > 0
        cm = translate(sphere_mesh(capsule_radius; ntheta = capsule_ntheta, nphi = capsule_nphi),
                       Vec3(capsule_center...))
        push!(comps, Component(:capsule, cm, capsule;
                               surface = SphereSurface(capsule_radius; center = capsule_center)))
    end
    append!(comps, window_comps)
    return Target(comps)
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
