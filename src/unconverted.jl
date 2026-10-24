# Unconverted (1ω, 2ω) light: footprint exposure of target hardware.

"""
    UnconvertedColor(name; footprint, edge, focal_shift)

Unconverted-light footprint of one color at the 3ω best-focus plane: a square of side
`footprint` (mm) whose nearest edge lies `edge` mm from the 3ω aim point, displaced in
a port-dependent direction; its focus lies `focal_shift` mm beyond the 3ω focus.
"""
struct UnconvertedColor
    name::Symbol
    footprint::Float64
    edge::Float64
    focal_shift::Float64
end
UnconvertedColor(name; footprint, edge, focal_shift) =
    UnconvertedColor(Symbol(name), footprint, edge, focal_shift)

"""
    nif_unconverted()

NIF 1ω and 2ω footprints (NIF User Guide §4.5, §5.6.1): 1ω ≈ 22 mm square with nearest
edge ≈ 4.8 mm from the 3ω aim point, 2ω ≈ 12 mm square with nearest edge ≈ 2.9 mm. Focal
shifts: 1ω +460 mm (MacFarlane 2003); 2ω ≈ +270 mm (derived from fused-silica
dispersion — see reference/research-notes.md).
"""
nif_unconverted() = [UnconvertedColor(:omega1; footprint = 22.0, edge = 4.8, focal_shift = 460.0),
                     UnconvertedColor(:omega2; footprint = 12.0, edge = 2.9, focal_shift = 270.0)]

"""
    UnconvertedRow

`beam`, `color`, `component`, `exposure` (`:exposed` — hardware inside the footprint
for a known offset direction; `:possible` — inside the footprint for some offset
direction when the direction is unknown; `:clear`), and `distance`: the smallest
transverse distance (mm) from the footprint region to the component (negative inside).
"""
struct UnconvertedRow
    beam::String
    color::Symbol
    component::Symbol
    exposure::Symbol
    distance::Float64
end

"""
    UnconvertedReport

Result of `unconverted_light_check`: per-beam/color/component `rows`, plus per-component
`zone` flags (nearest distance to TCC within `zone_radius`) and `flat` flags (planar
components, which violate the NIF curvature rule inside a 1ω footprint unless dimpled:
1D curvature radius ≤ 2 mm or 2D ≤ 5 mm).
"""
struct UnconvertedReport
    rows::Vector{UnconvertedRow}
    zone::Dict{Symbol,Float64}
    flat::Dict{Symbol,Bool}
    zone_radius::Float64
end

"""
    unconverted_light_check(target, beams; colors = nif_unconverted(), offsets = nothing,
                            hardware = all surfaces, zone_radius = 32.0, level = 1)

Check target hardware against unconverted-light footprints. Footprints are prisms
along each beam through its square footprint at the 3ω focal plane (near TCC the
converging 1ω/2ω beams change size by ≲ 2% per cm, so the prism is a close
approximation; shadowing is ignored, so exposure is conservative).

`offsets` maps beam names to the footprint offset direction: an angle (degrees, from
the port's θ̂ toward φ̂, in the focal plane) or a 3-vector. The per-port directions are
not public; without them every azimuth is considered and hardware is reported as
`:possible` when any offset direction would expose it.

Also reports components within `zone_radius` mm of TCC (where unconverted light
reaches) and planar components (flagged against the NIF curvature rule).
"""
function unconverted_light_check(t::Target, beams::AbstractVector{<:Beam};
                                 colors = nif_unconverted(), offsets = nothing,
                                 hardware = nothing, zone_radius = 32.0, level::Integer = 1)
    # hardware: parts (all their surface pieces) or components
    hw = hardware === nothing ? unique!([c.part for c in t.components if c.kind === :surface]) :
         Symbol[hardware...]
    hw_idx = [hardware === nothing ? [ci for ci in part_pieces(t, n) if t.components[ci].kind === :surface] :
              part_pieces(t, n) for n in hw]
    samples = Dict{Int,Vector{Vec3}}()
    for (ci, cis) in enumerate(hw_idx)
        pts = Vec3[]
        for k in eachindex(t.tri_area)
            t.tri_comp[k] in cis || continue
            v0, e1, e2 = t.tri_v0[k], t.tri_e1[k], t.tri_e2[k]
            _subtri_vertices!(pts, v0, v0 + e1, v0 + e2, level)
        end
        samples[ci] = unique(pts)
    end
    rows = UnconvertedRow[]
    for b in beams, col in colors, ci in eachindex(hw_idx)
        e1, e2 = port_basis(b.port)
        k = direction(b)
        aim = Vec3(ForwardDiff.value.(b.aim))
        dir = offsets === nothing ? nothing : get(offsets, b.name, nothing)
        h = col.footprint / 2
        best = Inf
        for x in samples[ci]
            d = x - aim
            u, v = dot(d, e1), dot(d, e2)
            if dir === nothing
                # union over all offset directions: an annulus from the nearest edge to
                # the far corner of the square
                ρ = hypot(u, v)
                rin, rout = col.edge, hypot(col.edge + col.footprint, h)
                dist = ρ < rin ? rin - ρ : ρ > rout ? ρ - rout : -min(ρ - rin, rout - ρ)
            else
                ô = dir isa Real ? (cosd(dir) * e1 + sind(dir) * e2) : normalize(Vec3(dir...) - dot(Vec3(dir...), k) * k)
                ô⊥ = cross(k, ô)
                c = col.edge + h                      # square center along the offset
                p = dot(d, ô) - c
                q = dot(d, ô⊥)
                ox, oy = abs(p) - h, abs(q) - h
                dist = (ox > 0 || oy > 0) ? hypot(max(ox, 0.0), max(oy, 0.0)) : max(ox, oy)
            end
            best = min(best, dist)
        end
        exposure = best >= 0 ? :clear : dir === nothing ? :possible : :exposed
        push!(rows, UnconvertedRow(b.name, col.name, hw[ci], exposure, best))
    end
    zone = Dict{Symbol,Float64}()
    flat = Dict{Symbol,Bool}()
    for (ci, cis) in enumerate(hw_idx)
        name = hw[ci]
        zone[name] = minimum(norm, samples[ci])
        ns = [t.tri_normal[k] for k in eachindex(t.tri_area) if t.tri_comp[k] in cis && t.tri_area[k] > 0]
        flat[name] = all(n -> abs(dot(n, ns[1])) > 1 - 1e-9, ns)
    end
    return UnconvertedReport(rows, zone, flat, Float64(zone_radius))
end

function Base.show(io::IO, ::MIME"text/plain", r::UnconvertedReport)
    nb = length(unique(x.beam for x in r.rows))
    println(io, "UnconvertedReport: ", nb, " beams")
    println(io, "  component        nearest to TCC   in zone   flat   1ω/2ω exposure (beams)")
    for (c, dmin) in sort(collect(r.zone); by = last)
        rs = [x for x in r.rows if x.component === c]
        ex = join(("$(col): " * string(count(x -> x.color === col && x.exposure !== :clear, rs)) *
                   (any(x -> x.color === col && x.exposure === :exposed, rs) ? " exposed" :
                    any(x -> x.color === col && x.exposure === :possible, rs) ? " possible" : "")
                   for col in unique(x.color for x in rs)), ", ")
        @printf(io, "  %-16s %10.2f mm     %-5s     %-5s  %s\n", c, dmin, dmin <= r.zone_radius,
                r.flat[c], ex)
    end
    any(c -> r.flat[c] && r.zone[c] <= r.zone_radius, keys(r.flat)) &&
        println(io, "  flat surfaces in the unconverted-light zone: check the NIF curvature rule ",
                "(1D radius ≤ 2 mm, 2D ≤ 5 mm, or dimpled)")
end
