# Ready-made target geometries.

"""
    hohlraum(; length = 9.4, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9,
             nz = 36, nphi = 48, nr_cap = 6, capsule_ntheta = 24, capsule_nphi = 32,
             wall = gold(), capsule = ablator(), capsule_center = (0,0,0))

Cylindrical hohlraum along `z`, centered on TCC (dimensions in mm), with components

- `:wall` — inner cylinder wall (radiates inward)
- `:cap_top`, `:cap_bottom` — end-cap washers between the LEH and the wall
- `:leh_top`, `:leh_bottom` — laser entrance holes (apertures; their incident flux is
  the x-ray loss through the LEH)
- `:capsule` — spherical capsule (omit with `capsule_radius = 0`)

Laser light that strikes the outside of the wall or caps (clipping on the LEH edge)
is reported as the `_outside` back-side loss of those components. The wall, caps,
LEHs and capsule carry their exact `CylinderSurface`, `DiskSurface` and `SphereSurface`,
so the radiation model uses the true geometry.
"""
function hohlraum(; length = 9.4, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9,
                  nz::Integer = 36, nphi::Integer = 48, nr_cap::Integer = 6, nr_leh::Integer = 4,
                  capsule_ntheta::Integer = 24, capsule_nphi::Integer = 32,
                  wall::Material = gold(), capsule::Material = ablator(),
                  capsule_center = (0.0, 0.0, 0.0))
    0 < leh_radius < radius || throw(ArgumentError("need 0 < leh_radius < radius"))
    capsule_radius < min(radius, length / 2) || throw(ArgumentError("capsule does not fit"))
    h = length / 2
    comps = Component[
        Component(:wall, cylinder_mesh(radius, -h, h; nz, nphi, inward = true), wall;
                  surface = CylinderSurface(radius)),
        Component(:cap_top, annulus_mesh(leh_radius, radius, h; nr = nr_cap, nphi, up = false), wall;
                  surface = DiskSurface(center = (0, 0, h))),
        Component(:cap_bottom, annulus_mesh(leh_radius, radius, -h; nr = nr_cap, nphi, up = true), wall;
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
