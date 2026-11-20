# This file places backlighter foils on a diagnostic line of sight and checks that the
# line of sight from the foil through the object to the detector is unobstructed.

"""
    backlighter_foil(los; standoff, width = 2.0, height = width, tilt = 0.0,
                     tilt_toward = 0.0, rotation = 0.0, center = (0, 0, 0), n = 8,
                     name = :backlighter, material = foil_material(),
                     two_sided = true) -> Component

Builds a rectangular backlighter foil for the diagnostic line of sight `los` = (θ, φ) in
degrees, the direction from the object toward the detector. The object is at `center`,
TCC by default. The foil is placed on the far side, centered `standoff` mm from `center`
at `center − standoff · port_vector(θ, φ)`. Its front face is heated by the backlighter
beams, and its x rays cross the object to the detector, which the front face is turned
toward. The placement is [`foil_frame`](@ref) with the same keywords. The mesh lies in
the `x`–`y` plane of that frame, with `width` along its `x` axis.

By default the foil normal lies along the line of sight. `tilt` (degrees) rotates the whole
foil about its center, off the line of sight toward `cos(ψ) θ̂ + sin(ψ) φ̂`. The angle
`ψ = tilt_toward` is in degrees from θ̂ of the line of sight toward φ̂, like the `angle` of
a spot. Before the tilt, the `width` edge of the foil lies along θ̂ and `height` along φ̂,
and `rotation` (degrees) turns both about the normal from θ̂ toward φ̂. A foil tilted by β
about its `height` axis (`tilt_toward = 0`) has a projected width `width · cos β` across the
line of sight.

The defaults (a 2 mm square gold foil) are generic engineering values that are not taken
from a facility backlighter design, and real backlighter dimensions and standoffs are
shot-specific. `n` sets the mesh (`n × n` patches for a square foil). The default
`material` is [`foil_material`](@ref)`()`, directly irradiated gold (albedo 0.8,
conversion 0.8); `foil_material(:Ti)` and the other mid-Z elements give a suggested point
inside a suggested range, listed in the table of `foil_material`, and none of these is a
measurement. The x-ray
emission of the foil set by this conversion is broadband radiosity; the He-α line signal
comes from [`kshell_ce`](@ref), fitted to separate data.
"""
function backlighter_foil(los; standoff::Real, width::Real = 2.0, height::Real = width,
                          tilt::Real = 0.0, tilt_toward::Real = 0.0, rotation::Real = 0.0,
                          center = (0.0, 0.0, 0.0), n::Integer = 8, name::Symbol = :backlighter,
                          material::Material = foil_material(), two_sided::Bool = true)
    standoff > 0 || throw(ArgumentError("the backlighter standoff must be positive"))
    width > 0 && height > 0 || throw(ArgumentError("the foil width and height must be positive"))
    f = foil_frame(los; standoff, tilt, tilt_toward, rotation, center)
    m = rectangle_mesh(width, height; nu = n, nv = max(1, round(Int, n * height / width)))
    m = transform(m, f.R, f.origin)
    return Component(name, m, material; two_sided)
end

"""
    LineOfSightReport

The result of `line_of_sight`, with fields:

- `clear`, true when no ray is blocked and the foil covers the object disk;
- `rays` and `blocked`, the rays in the bundle and the rays obstructed;
- `blocked_by`, pairs of hardware part => rays it obstructs, in target order, where a ray
  may count for several parts and rays that would leave the foil through its back count
  under the foil name;
- `covered`, the fraction of the disk sample points whose image ray meets the foil front;
- `samples`, the foil and disk sample counts;
- `los`, the unit vector toward the detector;
- `center` and `radius`, the object disk;
- `detector`, `nothing` or the detector distance used.
"""
struct LineOfSightReport
    clear::Bool
    rays::Int
    blocked::Int
    blocked_by::Vector{Pair{Symbol,Int}}
    covered::Float64
    samples::Tuple{Int,Int}
    los::Vec3
    center::Vec3
    radius::Float64
    detector::Union{Nothing,Float64}
end

function Base.show(io::IO, ::MIME"text/plain", r::LineOfSightReport)
    @printf(io, "LineOfSightReport: %s, %d of %d rays blocked, foil covers %.1f%% of the disk (radius %.3g mm)",
            r.clear ? "clear" : "obstructed", r.blocked, r.rays, 100 * r.covered, r.radius)
    for (p, k) in r.blocked_by
        k > 0 && print(io, "\n  ", p, ": ", k, " rays")
    end
end

# Samples the disk of radius `a` at `c` normal to `d` with its center and `nrings` rings
# of 8k points. Each ring starts on θ̂ of `d`, so the samples include the extremes of
# the disk along θ̂, φ̂ and the diagonals.
function _disk_samples(c::Vec3, d::Vec3, a, nrings)
    eθ, eφ = port_basis(d)
    pts = [c]
    for k in 1:nrings, j in 0:8k-1
        ψ = 2π * j / (8k)
        push!(pts, c + (a * k / nrings) * (cos(ψ) * eθ + sin(ψ) * eφ))
    end
    return pts
end

# Returns the nearest intersection (s, k) of `o + s d`, s > 0, with the triangles `tris`,
# or (Inf, 0) when the ray misses them all.
function _first_hit(t::Target, tris, o::Vec3, d::Vec3)
    best, hit = Inf, 0
    for k in tris
        s = _intersect(o, d, t.tri_v0[k], t.tri_e1[k], t.tri_e2[k], best)
        s < best && ((best, hit) = (s, k))
    end
    return best, hit
end

"""
    line_of_sight(target, los; foil = :backlighter, object = :capsule,
                  center = (0, 0, 0), radius = the object's projected radius,
                  hardware = nothing, level = 0, nrings = 4,
                  detector = nothing) -> LineOfSightReport

Checks whether the diagnostic line of sight `los` is unobstructed from the backlighter
`foil` through the `object` to the detector. `los` is θ, φ in degrees from `center`
toward the detector, or a direction vector. The object is represented by its disk normal
to the line of sight through `center`, of radius `radius`. The default radius is the
projected radius of the object, the largest distance of its vertices from the line of
sight; give `radius` when there is no object (`object = nothing`). The foil and object
parts never obstruct.

`hardware` names the parts or components that may obstruct, and `nothing` takes every
opaque part of the target except the foil and object. Apertures such as entrance holes
and open windows are transparent, as for rays everywhere; a closed window is opaque.

**The bundle.** With `detector = nothing` (the default), the bundle consists of every ray
from a foil sample point through a disk sample point, continued without end past the disk.
These rays are the x-ray paths from the whole emitting area of the foil through the object
toward any detector in their cone, a conservative superset of what one imager accepts. The
samples are taken as follows:

- Foil samples are the vertices of the foil mesh plus `4^level` sub-triangle points per
  triangle. The vertices include the corners and edges, where the extreme rays of the
  bundle start.
- Disk samples are the center and `nrings` rings of `8k` points, starting along θ̂ of the
  line of sight.

A ray counts as obstructed if it meets an opaque triangle of the hardware anywhere past
the foil. A ray that would leave the foil through its back face also counts as
obstructed, under the foil name.

`detector = L` takes the imaging rays instead, with `L` in mm from `center` along the
line of sight, or `Inf` for a distant imager. Each imaging ray is traced from the
detector point through a disk point back to the foil, or along the line of sight for
`Inf`. The segment from the foil to the detector is then checked, unbounded for `Inf`.

**Coverage.** `covered` is the fraction of disk points whose image ray meets the foil front.
The image ray is parallel to the line of sight, or directed toward the detector point when
`detector` is finite. Below 1 the foil does not backlight the whole object. `clear` requires
that no ray is obstructed and that `covered == 1`.

The result is exact for obstructions whose edges meet the bundle at its sampled
extremes, such as straight edges parallel to a foil edge and rims crossing the θ̂ or φ̂
extreme of the bundle. Elsewhere it is a sampled answer, and raising `level` and
`nrings` resolves a near-grazing obstruction.
"""
function line_of_sight(t::Target, los; foil::Symbol = :backlighter,
                       object::Union{Nothing,Symbol} = :capsule, center = (0.0, 0.0, 0.0),
                       radius::Union{Nothing,Real} = nothing, hardware = nothing,
                       level::Integer = 0, nrings::Integer = 4,
                       detector::Union{Nothing,Real} = nothing)
    d = los isa Tuple{Real,Real} ? port_vector(los...) : normalize(Vec3(los...))
    c = Vec3(center...)
    nrings >= 1 || throw(ArgumentError("need nrings ≥ 1"))
    detector === nothing || detector > 0 || throw(ArgumentError("the detector distance must be positive"))
    fi = part_pieces(t, foil)
    oi = object === nothing ? Int[] : part_pieces(t, object)
    if radius === nothing
        object === nothing && throw(ArgumentError("give `radius` when there is no `object`"))
        a = maximum(v -> norm((v - c) - dot(v - c, d) * d),
                    (v for ci in oi for v in t.components[ci].mesh.vertices))
    else
        a = Float64(radius)
    end
    a > 0 || throw(ArgumentError("the object's disk radius must be positive"))
    # List the hardware parts, one entry each, in target order.
    parts = hardware === nothing ?
            unique!([cp.part for (ci, cp) in enumerate(t.components)
                     if cp.kind === :surface && !(ci in fi) && !(ci in oi)]) :
            Symbol[hardware...]
    hw = [setdiff(part_pieces(t, p), fi, oi) for p in parts]
    hwall = Set(ci for cis in hw for ci in cis)
    foil_tris = [k for k in eachindex(t.tri_comp) if t.tri_comp[k] in fi]
    isempty(foil_tris) && throw(ArgumentError("the foil :$foil has no triangles"))

    disk = _disk_samples(c, d, a, nrings)
    # Each ray holds its origin on the foil, direction, tmax and whether it leaves the
    # front.
    rays = Tuple{Vec3,Vec3,Float64,Bool}[]
    far = 4 * (t.scale + norm(c)) + 1
    covered = 0
    for q in disk
        # Trace the image ray through q back toward the foil.
        P = detector === nothing || isinf(detector) ? nothing : c + Float64(detector) * d
        o, back = P === nothing ? (q + far * d, -d) : (P, q - P)
        s, k = _first_hit(t, foil_tris, o, back)
        hit = k > 0 && dot(back, t.tri_normal[k]) < 0
        covered += hit
        if detector !== nothing && k > 0
            p = o + s * back
            push!(rays, P === nothing ? (p, d, Inf, hit) : (p, P - p, 1 - 1e-9, hit))
        end
    end
    nfoil = 0
    if detector === nothing
        pts = Vec3[]
        for k in foil_tris
            v0 = t.tri_v0[k]
            _subtri_vertices!(pts, v0, v0 + t.tri_e1[k], v0 + t.tri_e2[k], level)
        end
        pts = unique(pts)
        nfoil = length(pts)
        # Take the front normal at each point from the nearest foil triangle; on a flat
        # foil every triangle has the same normal.
        normal_at(p) = t.tri_normal[argmin(k -> norm(t.tri_v0[k] + (t.tri_e1[k] + t.tri_e2[k]) / 3 - p),
                                           foil_tris)]
        for p in pts
            nf = normal_at(p)
            for q in disk
                push!(rays, (p, q - p, Inf, dot(q - p, nf) > 0))
            end
        end
    else
        nfoil = length(rays)                      # foil points of the image rays
    end

    skip_all(id) = !(t.tri_comp[id] in hwall)
    counts = zeros(Int, length(parts))
    back = 0
    blocked = 0
    for (o, dir, tmax, front) in rays
        isblocked = !front
        back += !front
        if any_hit(t.bvh, o, dir, tmax, skip_all)
            isblocked = true
            for (j, cis) in enumerate(hw)
                any_hit(t.bvh, o, dir, tmax, id -> !(t.tri_comp[id] in cis)) && (counts[j] += 1)
            end
        end
        blocked += isblocked
    end
    by = [p => k for (p, k) in zip(parts, counts)]
    back > 0 && push!(by, foil => back)
    cov = covered / length(disk)
    return LineOfSightReport(blocked == 0 && covered == length(disk), length(rays), blocked, by,
                             cov, (nfoil, length(disk)), d, c, a,
                             detector === nothing ? nothing : Float64(detector))
end
