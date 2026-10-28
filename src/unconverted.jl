# Unconverted (1ω, 2ω) light: footprint exposure of target hardware.

"""
    UnconvertedColor(name; footprint, edge, focal_shift, absorption = 1.0, conversion = 1.0)

Unconverted-light footprint of one color at the 3ω best-focus plane: a square of side
`footprint` (mm) whose nearest edge lies `edge` mm from the 3ω aim point, displaced in
a port-dependent direction; its focus lies `focal_shift` mm beyond the 3ω focus.

`absorption` is the fraction of this color's incident light a surface absorbs, relative
to the 3ω light (which the deposition takes as fully absorbed), and `conversion` the
x-ray conversion efficiency of the absorbed light relative to a material's 3ω
conversion. [`unconverted_deposition`](@ref) deposits `absorption` times the incident
power; a drive multiplies a material's conversion by `conversion` for this color.
"""
struct UnconvertedColor
    name::Symbol
    footprint::Float64
    edge::Float64
    focal_shift::Float64
    absorption::Float64
    conversion::Float64
end
function UnconvertedColor(name; footprint, edge, focal_shift, absorption = 1.0, conversion = 1.0)
    0 <= absorption <= 1 || throw(ArgumentError("absorption must be in [0, 1]"))
    conversion >= 0 || throw(ArgumentError("conversion must not be negative"))
    return UnconvertedColor(Symbol(name), footprint, edge, focal_shift, absorption, conversion)
end

"""
    nif_unconverted()

NIF 1ω and 2ω footprints (NIF User Guide §4.5, §5.6.1): 1ω ≈ 22 mm square with nearest
edge ≈ 4.8 mm from the 3ω aim point, 2ω ≈ 12 mm square with nearest edge ≈ 2.9 mm. Focal
shifts: 1ω +460 mm (MacFarlane 2003); 2ω ≈ +270 mm (derived from fused-silica
dispersion — see reference/research-notes.md).

Absorption and x-ray conversion relative to 3ω (an *engineering default*, confidence
low): the ratios on gold disks at the lowest intensity measured, ~3 × 10¹³ W/cm² with
~1 ns pulses, held below it.
- 1ω: absorption 0.82 against 0.93 at 3ω (Ze et al., UCRL-86414 (1981), Fig. 4) gives
  0.88; conversion of the absorbed energy 0.50 against 0.57 (Mead et al., LA-UR-86-300
  (1985), Fig. 16; same laboratory) gives 0.88.
- 2ω: absorption equal to 3ω's (Ze, Fig. 4) gives 1.0; conversion 0.72 (small spots)
  to 1.16 (large spots) of 3ω's (Mead, Fig. 16) gives their geometric mean, 0.9.

Unconverted light reaches hardware at ~10¹⁰–10¹² W/cm², 1.5–3.5 decades below these
data. The ratios rise toward 1 as the intensity falls, and inverse-bremsstrahlung
absorption is near complete in that regime. Before a plasma forms, cold metals reflect
1ω almost entirely, and cold silicon transmits it (research notes, "Unconverted
light: absorption and conversion vs. wavelength").
"""
nif_unconverted() = [UnconvertedColor(:omega1; footprint = 22.0, edge = 4.8, focal_shift = 460.0,
                                      absorption = 0.88, conversion = 0.88),
                     UnconvertedColor(:omega2; footprint = 12.0, edge = 2.9, focal_shift = 270.0,
                                      absorption = 1.0, conversion = 0.9)]

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

# ---------------------------------------------------------------- deposition

"""
    unconverted_deposition(target, beams; ratio, colors = nif_unconverted(), offsets = nothing,
                           lens_distance = 7700.0, coarse = 0.5, resolution = 0.02)
        -> NamedTuple

Laser power of the unconverted (1ω, 2ω) light absorbed by the target: one
[`Deposition`](@ref) per color of `colors`, keyed by its name (`:omega1`, `:omega2`),
and `skipped`, the beams whose footprints cannot be placed (no NIF quad name and no
`offsets` entry).

Each beamlet's light of one color is a uniform square beam converging to its own focus.
At the 3ω best-focus plane it is a square of side `footprint` whose nearest edge lies
`edge` mm from the 3ω axis ([`UnconvertedColor`](@ref); [`nif_unconverted`](@ref) gives
the NIF User Guide's 22 mm and 4.8 mm for 1ω, 12 mm and 2.9 mm for 2ω, §4.5 and
Fig. 4-12). It converges to a point `focal_shift` beyond the 3ω focus on the chief ray
from the lens center, `lens_distance` upstream. This is the geometry of
[`talis_check`](@ref)'s 1ω footprint, at the footprint's nominal size where TaLIS uses
its 10 % contour. The square is displaced along φ̂ to the side of the quad axis that
the beamlet's port lies on (Fig. 4-12), a NIF quad-level beam is expanded into its four
beamlets, and `offsets` sets the displacement direction of other beams, as in
[`talis_check`](@ref). The light at a target is the image of the near-field aperture,
so its intensity is taken to be uniform across the square (a flat top; the phase
plate's blur of the edge, ~0.1 mm, is not modeled). The geometry is *derived*
(confidence medium; research notes, "Unconverted light (NIF)").

`ratio` gives, per color, the unconverted power as a fraction of each beam's 3ω
`power`, e.g. `(omega1 = 0.05, omega2 = 0.01)`: a number, or a vector with one entry
per beam. It has no default: by the NIF User Guide (§4.5), most of the unconverted
light is 1ω and shaped pulses leave more of it, with no fraction given. VisRad also
takes it as an input (its "Fraction of Total Beam Power" and "Relative to Blue Cone
Power" models).

Without shadowing, the power a beam puts on a flat triangle is exact: every ray of the
beam passes through its focus, so the triangle's share of the power is the area of its
projection from the focus onto the 3ω focal plane, clipped to the footprint square,
over the square's area, at any tilt and distance from focus. Shadowing is resolved by
rays from the triangle toward the lens. The clipped projection is divided into pieces
no larger than `coarse` mm in the focal plane and, where a shadow edge crosses a piece,
down to `resolution` mm. The side of a surface facing the lens absorbs the color's
`absorption` times the power reaching it: its element or, for an unmodeled back side,
its component's `backside`. The rest is dropped, not reflected (`delivered -
sum(per_beam)` includes it). Light beyond its focus (hundreds of mm past the target)
is not followed. Values are Float64: this is a fixed source, not part of the
differentiable deposition path ([`deposit`](@ref)).
"""
function unconverted_deposition(t::Target, beams::AbstractVector{<:Beam}; ratio,
                                colors = nif_unconverted(), offsets = nothing,
                                lens_distance::Real = 7700.0, coarse::Real = 0.5,
                                resolution::Real = 0.02)
    0 < resolution <= coarse || throw(ArgumentError("need 0 < resolution ≤ coarse"))
    lens_distance > 0 || throw(ArgumentError("lens_distance must be positive"))
    vb = [_value_beam(b) for b in beams]
    table = any(_is_nif, vb) ? _nif_beam_table() : nothing
    bls = [_omega1_beamlets(b, offsets, table) for b in vb]
    tris = findall(t.tri_opaque)
    ε = _ray_eps(t)
    out = map(colors) do col
        haskey(ratio, col.name) || throw(ArgumentError("ratio has no value for $(col.name)"))
        r = ratio[col.name]
        rs = r isa Real ? fill(Float64(r), length(vb)) : Float64.(collect(r))
        length(rs) == length(vb) || throw(DimensionMismatch("ratio for $(col.name): one per beam"))
        all(>=(0), rs) || throw(ArgumentError("ratio for $(col.name) must not be negative"))
        h0 = col.footprint / 2
        c0 = col.edge + h0
        # one job per (beam, beamlet); each sums its triangles in order
        jobs = [(bi, bl) for (bi, bs) in enumerate(bls) if bs !== nothing && rs[bi] > 0 for bl in bs]
        parts = Vector{Vector{Tuple{Int32,Float64}}}(undef, length(jobs))
        Threads.@threads for j in eachindex(jobs)
            _, (k, fb, ô) = jobs[j]
            parts[j] = _uc_beamlet(t, tris, k, fb, ô, c0, h0, Float64(col.focal_shift),
                                   Float64(lens_distance), Float64(coarse), Float64(resolution), ε)
        end
        E = zeros(nelements(t))
        L = zeros(length(t.components))
        per_beam = zeros(length(vb))
        for (j, (bi, _)) in enumerate(jobs)
            P = col.absorption * rs[bi] * vb[bi].power / length(bls[bi])
            for (d, f) in parts[j]
                d > 0 ? (E[d] += P * f) : (L[-d] += P * f)
                per_beam[bi] += P * f
            end
        end
        delivered = sum((rs[bi] * vb[bi].power for bi in eachindex(vb) if bls[bi] !== nothing); init = 0.0)
        col.name => Deposition(E, L, per_beam, delivered)
    end
    skipped = String[vb[bi].name for bi in eachindex(vb) if bls[bi] === nothing]
    return (; out..., skipped)
end

# One beamlet's share (a fraction of its power) per destination: an element (> 0) or a
# component's back side (< 0). The beam is the pyramid from the footprint square (center
# c0 along ô, half-side h0, in the 3ω focal plane through `fb`) to its apex, the focus.
function _uc_beamlet(t::Target, tris, k, fb, ô, c0, h0, f, L, coarse, res, ε)
    ŵ = cross(k, ô)
    A = fb + f * k + c0 * (1 + f / L) * ô
    out = Tuple{Int32,Float64}[]
    sq = (c0 - h0, c0 + h0, -h0, h0)
    for tri in tris
        n = t.tri_normal[tri]
        v0 = t.tri_v0[tri]
        V = (v0, v0 + t.tri_e1[tri], v0 + t.tri_e2[tri])
        # only light still converging (upstream of the focus) is followed
        all(x -> dot(x - fb, k) < f * (1 - 1e-9), V) || continue
        σ = dot(n, A - v0)                           # the same at every point of the plane
        σ == 0 && continue                           # edge-on to every ray
        side_n = σ < 0 ? n : -n                      # the side facing the lens
        e = σ < 0 ? t.tri_front[tri] : t.tri_back[tri]
        dest = e > 0 ? Int32(e) : -Int32(t.tri_comp[tri])
        uv = map(x -> _uc_project(x, A, fb, k, ô, ŵ, f), V)
        (minimum(first, uv) >= sq[2] || maximum(first, uv) <= sq[1] ||
         minimum(last, uv) >= sq[4] || maximum(last, uv) <= sq[3]) && continue
        poly = _clip_box(collect(uv), sq)
        length(poly) < 3 && continue
        lit = 0.0
        for i in 2:length(poly)-1
            lit += _uc_lit(t, (poly[1], poly[i], poly[i+1]), A, fb, ô, ŵ, n, v0, side_n, ε, coarse, res)
        end
        lit > 0 && push!(out, (dest, lit / 4h0^2))
    end
    return out
end

# Focal-plane coordinates (along ô, ŵ from the 3ω focus) of the ray from x to the apex.
function _uc_project(x, A, fb, k, ô, ŵ, f)
    τ = f / (f - dot(x - fb, k))
    p = A + τ * (x - A) - fb
    return (dot(p, ô), dot(p, ŵ))
end

# Sutherland–Hodgman clip of a convex polygon of (u, v) points to the box (u0, u1, v0, v1).
function _clip_box(poly, (u0, u1, v0, v1))
    for (axis, lim, below) in ((1, u0, false), (1, u1, true), (2, v0, false), (2, v1, true))
        isempty(poly) && break
        inside(p) = below ? p[axis] <= lim : p[axis] >= lim
        out = similar(poly, 0)
        for i in eachindex(poly)
            a, b = poly[i], poly[mod1(i + 1, length(poly))]
            ia, ib = inside(a), inside(b)
            ia && push!(out, a)
            if ia != ib
                s = (lim - a[axis]) / (b[axis] - a[axis])
                push!(out, (a[1] + s * (b[1] - a[1]), a[2] + s * (b[2] - a[2])))
            end
        end
        poly = out
    end
    return poly
end

_uc_area(a, b, c) = abs((b[1] - a[1]) * (c[2] - a[2]) - (c[1] - a[1]) * (b[2] - a[2])) / 2

# Lit focal-plane area of a triangle (u, v) of a clipped projection: its corners and
# centroid are mapped back to the target triangle and tested for shadow toward the lens.
function _uc_lit(t, tri, A, fb, ô, ŵ, n, v0, side_n, ε, coarse, res)
    a, b, c = tri
    area = _uc_area(a, b, c)
    area == 0 && return 0.0
    size = max(hypot(b[1] - a[1], b[2] - a[2]), hypot(c[1] - b[1], c[2] - b[2]),
               hypot(a[1] - c[1], a[2] - c[2]))
    size > coarse && return _uc_split(t, tri, A, fb, ô, ŵ, n, v0, side_n, ε, coarse, res)
    m = ((a[1] + b[1] + c[1]) / 3, (a[2] + b[2] + c[2]) / 3)
    seen = count(p -> _uc_visible(t, p, A, fb, ô, ŵ, n, v0, side_n, ε), (a, b, c, m))
    seen == 4 && return area                      # uniformly lit or dark at this size
    seen == 0 && return 0.0
    size <= res && return area * seen / 4
    return _uc_split(t, tri, A, fb, ô, ŵ, n, v0, side_n, ε, coarse, res)
end

function _uc_split(t, (a, b, c), args...)
    mid(p, q) = ((p[1] + q[1]) / 2, (p[2] + q[2]) / 2)
    ab, bc, ca = mid(a, b), mid(b, c), mid(c, a)
    return _uc_lit(t, (a, ab, ca), args...) + _uc_lit(t, (ab, b, bc), args...) +
           _uc_lit(t, (ca, bc, c), args...) + _uc_lit(t, (ab, bc, ca), args...)
end

# Whether the target point on the ray through focal-plane point p sees the lens.
function _uc_visible(t, p, A, fb, ô, ŵ, n, v0, side_n, ε)
    q = fb + p[1] * ô + p[2] * ŵ
    x = A + dot(n, v0 - A) / dot(n, q - A) * (q - A)
    return !any_hit(t.bvh, x + ε * side_n, normalize(x - A))
end
