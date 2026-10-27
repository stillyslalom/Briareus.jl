# Backlighter foils placed on a diagnostic line of sight, and the check that the line
# of sight from the foil through the object to the detector is unobstructed.

# Unit vector `los` rotated by `tilt` (rad) toward `e` (a unit vector ⟂ los), and the
# rotation that does it, applied to `v`.
function _tilted(v::Vec3, los::Vec3, e::Vec3, tilt)
    k = cross(los, e)                       # rotation axis (unit: los ⟂ e)
    return cos(tilt) * v + sin(tilt) * cross(k, v) + (1 - cos(tilt)) * dot(k, v) * k
end

"""
    backlighter_foil(los; standoff, width = 2.0, height = width, tilt = 0.0,
                     tilt_toward = 0.0, rotation = 0.0, center = (0, 0, 0), n = 8,
                     name = :backlighter, material = gold(conversion = 0.05),
                     two_sided = true) -> Component

Rectangular backlighter foil for the diagnostic line of sight `los` = (θ, φ) in
degrees, the direction from the object (`center`, default TCC) toward the detector. The
foil sits on the far side, centered `standoff` mm from `center` at
`center − standoff · port_vector(θ, φ)`, with its front face (the face the
backlighter beams heat, whose x rays cross the object to the detector) toward the
detector.

By default the foil's normal lies along the line of sight. `tilt` (degrees) turns it
off the line of sight, toward `cos(ψ) θ̂ + sin(ψ) φ̂` with `ψ = tilt_toward` (degrees,
from the line of sight's θ̂ toward φ̂, as a spot's `angle`), rotating the whole foil
about its center. Before the tilt, the foil's `width` edge runs along θ̂ and `height`
along φ̂, turned by `rotation` (degrees) about the normal from θ̂ toward φ̂. A foil
tilted by β about its `height` axis (`tilt_toward = 0`) presents `width · cos β` across
the line of sight.

The defaults (a 2 mm square gold foil) are generic engineering values, not a facility's
backlighter design: real backlighter dimensions and standoffs are shot-specific. `n` sets the mesh (`n × n` patches for a square foil).
"""
function backlighter_foil(los; standoff::Real, width::Real = 2.0, height::Real = width,
                          tilt::Real = 0.0, tilt_toward::Real = 0.0, rotation::Real = 0.0,
                          center = (0.0, 0.0, 0.0), n::Integer = 8, name::Symbol = :backlighter,
                          material::Material = gold(conversion = 0.05), two_sided::Bool = true)
    standoff > 0 || throw(ArgumentError("the backlighter standoff must be positive"))
    width > 0 && height > 0 || throw(ArgumentError("the foil width and height must be positive"))
    d = port_vector(los...)
    eθ, eφ = port_basis(d)
    u = cosd(rotation) * eθ + sind(rotation) * eφ
    v = -sind(rotation) * eθ + cosd(rotation) * eφ
    e = cosd(tilt_toward) * eθ + sind(tilt_toward) * eφ
    β = deg2rad(tilt)
    u, v, nrm = _tilted(u, d, e, β), _tilted(v, d, e, β), _tilted(d, d, e, β)
    m = rectangle_mesh(width, height; nu = n, nv = max(1, round(Int, n * height / width)))
    R = SMatrix{3,3,Float64}(u[1], u[2], u[3], v[1], v[2], v[3], nrm[1], nrm[2], nrm[3])
    m = transform(m, R, Vec3(center...) - standoff * d)
    return Component(name, m, material; two_sided)
end

"""
    LineOfSightReport

Result of `line_of_sight`: `clear` (no ray blocked and the foil covers the object's
disk), `rays` and `blocked` (rays in the bundle and rays obstructed), `blocked_by`
(hardware part => rays it obstructs, in target order; a ray may count for several
parts; rays that would leave the foil through its back count under the foil's name),
`covered` (fraction of the disk's sample points whose image ray meets the foil's
front), `samples` (foil and disk sample counts), `los` (unit vector toward the
detector), `center` and `radius` (the object's disk) and `detector` (`nothing`, or the
detector distance used).
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

# Disk of radius `a` at `c` normal to `d`: center plus `nrings` rings of 8k points, each
# ring starting on θ̂ of `d` (so the samples include the disk's extremes along θ̂, φ̂
# and the diagonals).
function _disk_samples(c::Vec3, d::Vec3, a, nrings)
    eθ, eφ = port_basis(d)
    pts = [c]
    for k in 1:nrings, j in 0:8k-1
        ψ = 2π * j / (8k)
        push!(pts, c + (a * k / nrings) * (cos(ψ) * eθ + sin(ψ) * eφ))
    end
    return pts
end

# Nearest intersection of `o + s d`, s > 0, with the triangles `tris`: (s, k), or (Inf, 0).
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

Whether the diagnostic line of sight `los` (θ, φ in degrees from `center` toward the
detector, or a direction vector) is unobstructed from the backlighter `foil` through
the `object` to the detector. The object is represented by its disk normal to the line
of sight through `center`, of `radius` (default: the largest distance of the object's
vertices from the line of sight, its projected radius; give `radius` when there is no
`object`, `object = nothing`). The foil and object parts never obstruct.

`hardware` names the parts (or components) that may obstruct; `nothing` takes every
opaque part of the target but the foil and object. Apertures (entrance holes, open
windows) are transparent, as for rays everywhere; a closed window is opaque.

**The bundle.** With `detector = nothing` (the default), every ray from a foil sample
point through a disk sample point, continued without end past the disk: the x-ray
paths from the foil's whole emitting area through the object toward any detector in
their cone, a conservative superset of what one imager accepts. Foil samples are the
foil mesh's vertices plus `4^level` sub-triangle points per triangle (vertices
include the corners and edges, where the bundle's extreme rays start); disk samples
are the center and `nrings` rings of `8k` points, starting along the line of sight's
θ̂. A ray counts as obstructed if it meets an opaque triangle of the hardware
anywhere past the foil, or if it would leave the foil through its back face (counted
under the foil's name).

`detector = L` (mm from `center` along the line of sight; `Inf` for a distant
imager) instead takes the imaging rays: from the detector point (or along the line of
sight) through each disk point back to the foil, and checks the segment from the foil
to the detector (unbounded for `Inf`).

**Coverage.** `covered` is the fraction of disk points whose image ray (parallel to
the line of sight, or toward the detector point when `detector` is finite) meets the
foil's front: below 1 the foil does not backlight the whole object. `clear` requires
no obstructed ray and `covered == 1`.

Sampling makes the result exact for obstructions whose edges meet the bundle at its
sampled extremes (a straight edge parallel to a foil edge, a rim crossing the bundle's
θ̂ or φ̂ extreme) and otherwise a sampled answer: raise `level` and `nrings` to resolve
a near-grazing obstruction.
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
    # hardware: parts, one entry each (in target order)
    parts = hardware === nothing ?
            unique!([cp.part for (ci, cp) in enumerate(t.components)
                     if cp.kind === :surface && !(ci in fi) && !(ci in oi)]) :
            Symbol[hardware...]
    hw = [setdiff(part_pieces(t, p), fi, oi) for p in parts]
    hwall = Set(ci for cis in hw for ci in cis)
    foil_tris = [k for k in eachindex(t.tri_comp) if t.tri_comp[k] in fi]
    isempty(foil_tris) && throw(ArgumentError("the foil :$foil has no triangles"))

    disk = _disk_samples(c, d, a, nrings)
    # rays: (origin on the foil, direction, tmax, front?)
    rays = Tuple{Vec3,Vec3,Float64,Bool}[]
    far = 4 * (t.scale + norm(c)) + 1
    covered = 0
    for q in disk
        # the image ray through q, back toward the foil
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
        # front normal at each point: of the nearest foil triangle (a flat foil: any)
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
