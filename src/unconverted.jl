# This file models unconverted (1ω, 2ω) light and the exposure of target hardware to its
# footprints.

"""
    UnconvertedColor(name; footprint, edge, focal_shift, absorption = 1.0, conversion = 1.0)

Unconverted-light footprint of one color at the 3ω best-focus plane. The footprint is a
square of side `footprint` (mm) whose nearest edge lies `edge` mm from the 3ω aim point,
displaced in a port-dependent direction. Its focus lies `focal_shift` mm beyond the 3ω
focus.

`absorption` is the fraction of the incident light of this color that a surface absorbs,
relative to the 3ω light, which the deposition takes as fully absorbed. `conversion` is
the x-ray conversion efficiency of the absorbed light relative to the 3ω conversion of a
material. [`unconverted_deposition`](@ref) deposits `absorption` times the incident
power, and [`solve`](@ref) with `unconverted` multiplies the conversion of a material, at
the absorbed flux of this color, by `conversion`.
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

NIF 1ω and 2ω footprints (NIF User Guide §4.5, §5.6.1). The 1ω footprint is a square of
about 22 mm with its nearest edge about 4.8 mm from the 3ω aim point, and the 2ω
footprint a square of about 12 mm with its nearest edge about 2.9 mm away. The focal
shift is +460 mm for 1ω (MacFarlane 2003) and about +270 mm for 2ω, derived from
fused-silica dispersion (see reference/research-notes.md).

The absorption and x-ray conversion relative to 3ω are an *engineering default* with
low confidence. They are the ratios on gold disks at the lowest intensity measured,
about 3 × 10¹³ W/cm² with pulses of about 1 ns, held constant below that intensity.
- 1ω: absorption 0.82 against 0.93 at 3ω (Ze et al., UCRL-86414 (1981), Fig. 4) gives
  0.88; conversion of the absorbed energy 0.50 against 0.57 (Mead et al., LA-UR-86-300
  (1985), Fig. 16; same laboratory) gives 0.88.
- 2ω: absorption equal to that at 3ω (Ze, Fig. 4) gives 1.0. The conversion is 0.72
  (small spots) to 1.16 (large spots) of that at 3ω (Mead, Fig. 16), and the geometric
  mean, 0.9, is used.

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

One row of an unconverted-light check, with fields `beam`, `color`, `component`,
`exposure` and `distance`. `exposure` is `:exposed` for hardware inside the footprint
for a known offset direction, `:possible` for hardware inside the footprint for some
offset direction when the direction is unknown, and `:clear` otherwise. `distance` is
the smallest transverse distance (mm) from the footprint region to the component,
negative inside.
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

Result of `unconverted_light_check`. It contains `rows` per beam, color and component,
and per-component `zone` and `flat` flags. `zone` marks a component whose nearest
distance to TCC is within `zone_radius`. `flat` marks a planar component, which violates
the NIF curvature rule inside a 1ω footprint unless dimpled; the rule requires a 1D
curvature radius ≤ 2 mm or a 2D radius ≤ 5 mm.
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

Check target hardware against unconverted-light footprints. Footprints are prisms along
each beam through its square footprint at the 3ω focal plane. Near TCC the converging
1ω and 2ω beams change size by ≲ 2% per cm, so the prism is a close approximation.
Shadowing is ignored, so the exposure is conservative.

`offsets` maps beam names to the footprint offset direction, given as an angle (degrees,
from θ̂ of the port toward φ̂, in the focal plane) or as a 3-vector. The per-port
directions are not public. Without them every azimuth is considered, and hardware is
reported as `:possible` when any offset direction would expose it.

Also reports components within `zone_radius` mm of TCC (where unconverted light
reaches) and planar components (flagged against the NIF curvature rule).
"""
function unconverted_light_check(t::Target, beams::AbstractVector{<:Beam};
                                 colors = nif_unconverted(), offsets = nothing,
                                 hardware = nothing, zone_radius = 32.0, level::Integer = 1)
    # hardware as parts (all their surface pieces) or components
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
                # The union over all offset directions is an annulus from the nearest edge
                # to the far corner of the square.
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
                           lens_distance = NIF_LENS_DISTANCE, coarse = 0.5, resolution = 0.02)
        -> NamedTuple

Laser power of the unconverted (1ω, 2ω) light absorbed by the target. The result
contains one [`Deposition`](@ref) per color of `colors`, keyed by its name (`:omega1`,
`:omega2`). It also contains `skipped`, the beams whose footprints cannot be placed (no
NIF quad name and no `offsets` entry), and `colors`, the [`UnconvertedColor`](@ref)s
themselves. [`solve`](@ref) and [`flux_response`](@ref) then use the x-ray conversions
of the result with its light.

The light of one color from each beamlet is a uniform square beam converging to its
focus. At the 3ω best-focus plane it is a square of side `footprint` whose nearest edge
lies `edge` mm from the 3ω axis ([`UnconvertedColor`](@ref)). [`nif_unconverted`](@ref)
gives the values of the NIF User Guide, 22 mm and 4.8 mm for 1ω and 12 mm and 2.9 mm
for 2ω (§4.5 and Fig. 4-12). The beam converges to a point `focal_shift` beyond the 3ω
focus on the chief ray from the lens center, `lens_distance` upstream. This is the
geometry of the 1ω footprint of [`talis_check`](@ref), at the nominal size of the
footprint where TaLIS uses its 10% contour. The square is displaced along φ̂ to the side
of the quad axis on which the port of the beamlet lies (Fig. 4-12). A NIF quad-level
beam is expanded into its four beamlets, and `offsets` sets the displacement direction
of other beams, as in [`talis_check`](@ref). The light at a target is the image of the
near-field aperture, so its intensity is taken to be uniform across the square, as a
flat top. The blur of the edge by the phase plate, about 0.1 mm, is not modeled. The
geometry is *derived* (confidence medium; research notes, "Unconverted light (NIF)").

`ratio` gives, per color, the unconverted power as a fraction of the 3ω `power` of each
beam, e.g. `(omega1 = 0.05, omega2 = 0.01)`, as a number or as a vector with one entry
per beam. It has no default. According to the NIF User Guide (§4.5), most of the
unconverted light is 1ω and shaped pulses leave more of it, but the guide does not give
a fraction. VisRad also takes the fraction as an input, in its "Fraction of Total Beam
Power" and "Relative to Blue Cone Power" models.

Without shadowing, the power that a beam puts on a flat triangle is exact. Every ray of
the beam passes through its focus, so the share of the power on the triangle is the area
of its projection from the focus onto the 3ω focal plane, clipped to the footprint
square, over the area of the square. This holds at any tilt and distance from focus.
Shadowing is resolved by rays from the triangle toward the lens. The clipped projection
is divided into pieces of at most `coarse` mm in the focal plane and, where a shadow
edge crosses a piece, down to `resolution` mm. The side of a surface facing the lens
absorbs `absorption` times the power of the color reaching it. The power goes to its
element or, for an unmodeled back side, to the `backside` of its component. The rest is
dropped without reflection, and `delivered - sum(per_beam)` includes it. Light beyond
its focus (hundreds of mm past the target) is not followed. Values are Float64, since
this is a fixed source and lies outside the differentiable deposition path
([`deposit`](@ref)).
"""
function unconverted_deposition(t::Target, beams::AbstractVector{<:Beam}; ratio,
                                colors = nif_unconverted(), offsets = nothing,
                                lens_distance::Real = NIF_LENS_DISTANCE, coarse::Real = 0.5,
                                resolution::Real = 0.02)
    0 < resolution <= coarse || throw(ArgumentError("need 0 < resolution ≤ coarse"))
    lens_distance > 0 || throw(ArgumentError("lens_distance must be positive"))
    vb = [_value_beam(b) for b in beams]
    table = any(_is_nif, vb) ? _nif_beam_table() : nothing
    bls = [_omega1_beamlets(b, offsets, table) for b in vb]
    tris = findall(t.tri_opaque)
    ε = _ray_eps(t)
    _check_color_keys(ratio, Symbol[c.name for c in colors], "ratio")
    out = map(colors) do col
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
    return (; out..., skipped, colors = collect(UnconvertedColor, colors))
end

# Returns the share of one beamlet (a fraction of its power) per destination, which is an
# element (> 0) or the back side of a component (< 0). The beam is the pyramid from the
# footprint square to its apex, the focus. The square has center c0 along ô and half-side
# h0, in the 3ω focal plane through `fb`.
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

# Returns the lit focal-plane area of a triangle (u, v) of a clipped projection. Its
# corners and centroid are mapped back to the target triangle and tested for shadow
# toward the lens.
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

# Whether the target point on the ray through focal-plane point p has a clear line of
# sight to the lens.
function _uc_visible(t, p, A, fb, ô, ŵ, n, v0, side_n, ε)
    q = fb + p[1] * ô + p[2] * ŵ
    x = A + dot(n, v0 - A) / dot(n, q - A) * (q - A)
    return !any_hit(t.bvh, x + ε * side_n, normalize(x - A))
end

# ---------------------------------------------------------------- irradiance histories

"""
    IrradianceHistory

Incident laser irradiance on a target over time, per color. It is the result of
[`unconverted_history`](@ref) (colors `:omega1`, `:omega2`) and of
[`irradiance_history`](@ref) (3ω, `:omega3`).

- `target`, `time` (ns), `colors` (the color names, in order)
- `element` — per color, a `nelements × length(time)` matrix of the irradiance incident
  on each element (W/cm²), which is the power reaching it over its area, at its incidence
- `backside` — per color, a `ncomponents × length(time)` matrix of the irradiance on the
  unmodeled back side of each one-sided component (W/cm²). It is the power over the area
  of the component, a mean that lies below the local peak of the back side where it is
  lit unevenly.
- `absorption` — per color, the fraction of the incident light that a surface absorbs
  (the `absorption` of the `UnconvertedColor`, 1 for 3ω). The absorbed irradiance is
  `absorption[c] .* element[c]`.
- `skipped` — the beams whose unconverted light cannot be placed (see
  [`unconverted_deposition`](@ref))

Between its times the history is taken to be linear, and crossings of a threshold
([`threshold_crossing`](@ref)) are interpolated within a step.
"""
struct IrradianceHistory
    target::Target
    time::Vector{Float64}
    colors::Vector{Symbol}
    element::Dict{Symbol,Matrix{Float64}}
    backside::Dict{Symbol,Matrix{Float64}}
    absorption::Dict{Symbol,Float64}
    skipped::Vector{String}
end

function Base.show(io::IO, ::MIME"text/plain", h::IrradianceHistory)
    nt = length(h.time)
    print(io, "IrradianceHistory: ", nt, " times")
    nt > 0 && @printf(io, " (%.4g–%.4g ns)", h.time[1], h.time[end])
    for c in h.colors
        pk = max(maximum(h.element[c]; init = 0.0), maximum(h.backside[c]; init = 0.0))
        @printf(io, "\n  %-7s peak %.3g W/cm²", c, pk)
    end
    isempty(h.skipped) || print(io, "\n  skipped (cannot be placed): ", join(h.skipped, ", "))
end

# Returns the area of each one-sided component, which is also the area of its back side,
# and 0 for a two-sided component.
_backside_area(t::Target) =
    [c.two_sided ? 0.0 :
     sum(t.elem_area[e] for e in t.comp_elems[ci] if t.elem_side[e] == 1; init = 0.0)
     for (ci, c) in enumerate(t.components)]

# Groups beams by their schedule `sig(b)`, compared by identity. Returns the distinct
# schedules and the beam indices of each, in order of first appearance.
function _schedule_groups(beams, sig)
    sigs = Any[]
    groups = Vector{Int}[]
    for (i, b) in enumerate(beams)
        s = sig(b)
        j = findfirst(k -> k === s, sigs)
        j === nothing ? (push!(sigs, s); push!(groups, [i])) : push!(groups[j], i)
    end
    return sigs, groups
end

# Per-unit-power irradiance (W/cm²) of a deposition, on elements and back sides.
_unit_irradiance(t::Target, d::Deposition, Aback) =
    (1e14 .* d.element ./ t.elem_area,
     [A > 0 ? 1e14 * L / A : 0.0 for (L, A) in zip(d.backside, Aback)])

# Adds `I` (per unit power) times the power history `s` to the matrix `M`.
function _add_outer!(M, I, s)
    for k in eachindex(s)
        s[k] == 0 && continue
        @views M[:, k] .+= s[k] .* I
    end
    return M
end

# A per-color input (`ratio` or `prepulse`, as a Dict or NamedTuple) is keyed exactly by
# the color names as Symbols, as `pulses` is keyed by the beam keys. Any other key, such
# as a String or a misspelled or unknown name, is an error, so that every entry is used.
# With `complete`, every color needs an entry.
function _check_color_keys(d, names, what; complete::Bool = true)
    for n in keys(d)
        n isa Symbol && n in names ||
            throw(ArgumentError("$what has an entry for $(repr(n)), which is not one of the colors " *
                                "$(join(repr.(names), ", ")) (keys are the colors' names, Symbols)"))
    end
    if complete
        for n in names
            haskey(d, n) || throw(ArgumentError("$what has no value for $n"))
        end
    end
    return nothing
end

_fraction_at(r::Real, τ) = Float64(r)
_fraction_at(r, τ) = Float64(r(τ))

# Returns the fraction schedule of beam `b`, at position `i` of `n` beams, from the entry
# `r` of one color in a per-color fraction input. The schedule is `r` itself (a number or
# a function of time, for every beam), `r[i]` of a vector with one entry per beam, or
# `r[key(b)]` of a keyed collection. A beam without an entry is a `KeyError`, as for
# `pulses`. Numbers are converted to Float64, so that beams on equal fractions share a
# schedule, compared by identity.
function _fraction_of(r, key, b, i, n)
    f = if r isa AbstractVector
        length(r) == n || throw(DimensionMismatch("a fraction per beam needs one entry for each of the $n beams, " *
                                                  "not $(length(r))"))
        r[i]
    else
        _schedule_of(r, key, b)
    end
    return f isa Real ? Float64(f) : f
end

# Returns the fraction schedules of each beam, a tuple over the colors `names` of `ratio`.
_fractions_of(ratio, names, key, beams) =
    [Tuple(_fraction_of(ratio[nm], key, b, i, length(beams)) for nm in names) for (i, b) in enumerate(beams)]

# Returns the schedule of beam `b`, `p[key(b)]` of a keyed collection, or else `p` itself
# as one function of time for all beams.
_schedule_of(p, key, b) = p isa Union{AbstractDict,NamedTuple} ? _beam_pulse(p, key, b) : p

# Returns the power of a color per beam (TW) at the times `ts`. It is `r` (a number or a
# function of time) times the 3ω pulse, plus a separate pulse `pre` for the color, if any.
function _color_power(r, pulse, pre, ts, name)
    s = [_fraction_at(r, τ) * pulse(τ) + (pre === nothing ? 0.0 : Float64(pre(τ))) for τ in ts]
    all(>=(0), s) || throw(ArgumentError("the $name power must not be negative (ratio and prepulse)"))
    return s
end

"""
    unconverted_history(target, beams, pulses, times; ratio, prepulse = (;), by = :group,
                        colors = nif_unconverted(), offsets = nothing, lens_distance = NIF_LENS_DISTANCE,
                        coarse = 0.5, resolution = 0.02) -> IrradianceHistory

Irradiance of the unconverted (1ω, 2ω) light on the target at each of `times` (ns), from
the 3ω pulse of each beam and an unconverted fraction that may vary in time.

`pulses` gives the 3ω power history of each beam (TW per beam), keyed as in
[`with_pulses`](@ref). It is one function of time (a [`PulseShape`](@ref)) for all
beams, or a `Dict` or NamedTuple keyed by the `group` of each beam. With `by = :quad`
the key is the quad, with `by = :beam` the name, and otherwise `by(beam)`. For each
color `c` of `colors`, the light of one beam is

    P_c(t) = ratio[c](t) · P₃ω(t) + prepulse[c](t)    (TW per beam)

- `ratio` — per color, the unconverted power as a fraction of the 3ω power of the beam,
  as a number or a function of time (ns), e.g. a `PulseShape(times, fractions)` table
  (zero outside its times). It has no default, since the fraction and its dependence on
  the pulse shape are not public. According to the NIF User Guide (§4.5), most of the
  unconverted light is 1ω and shaped pulses leave more of it. Every color needs an
  entry (0 for none), keyed by its name as a `Symbol`; any other key, such as a `String`
  or an unknown name, is an error. The entry of a color applies to every beam, or gives
  each beam a separate fraction (a number or a function of time). A separate fraction is
  a vector with one entry per beam of `beams`, in order, or is keyed as `pulses` by `by`.
  A beam without an entry is a `KeyError`, as for `pulses`. An example for 1ω on the
  outer cones only is `(omega1 = Dict(:outer44 => 0.05, :outer50 => 0.05,
  :inner23 => 0.0, :inner30 => 0.0), omega2 = 0.01)`.
- `prepulse` — per color (optional), light not proportional to the 3ω pulse, in TW per
  beam. It is one function of time (a `PulseShape`) for all beams or is keyed as
  `pulses`, and it is keyed by color as `ratio` is. An example is the 1ω prepulse, which
  the User Guide (§4.4) puts above 10⁸ W/cm² 3–5 ns before the nominal start of the
  pulse, with no level given. As a *derived* scale, the 22 mm 1ω footprint (4.84 cm² in
  the 3ω focal plane) is at 10⁸ W/cm² when a beamlet carries 0.48 GW (4.8e-4 TW) of 1ω.

The geometry is fixed, with the aims, spots and ports of `beams`, and only the power
varies. The light is linear in the power of each beam, so the whole history costs one
[`unconverted_deposition`](@ref). The beams are grouped by their (pulse, prepulse,
fraction) schedule (equal numbers or the same function), and each group is deposited
once at unit power. The history is that pattern times the power of the group at each
time. A pointing or a target that changes in time needs a deposition per configuration,
which this function does not cover.

The irradiance is the light incident on each surface, with the `absorption` of each
color set to 1, since the threshold checks ([`specular_check`](@ref),
[`preplasma_check`](@ref)) compare incident intensities. Multiplied by the area of each
element (and divided by 10¹⁴ for W/cm² to TW/mm²), `absorption[c] .* element[c]` is the
power that [`unconverted_deposition`](@ref) absorbs at the beam powers of each time. For
a constant fraction the two agree within rounding. `offsets`, `lens_distance`, `coarse`
and `resolution` are those of [`unconverted_deposition`](@ref). The history is sampled
at `times` and is linear between them, so the knots of the pulses and fractions should be
included in `times` to follow them exactly.
"""
Base.@nospecializeinfer function unconverted_history(t::Target, beams::AbstractVector{<:Beam},
                                                     @nospecialize(pulses), @nospecialize(times);
                                                     @nospecialize(ratio), @nospecialize(prepulse = (;)),
                                                     @nospecialize(by = :group),
                                                     colors = nif_unconverted(), offsets = nothing,
                                                     lens_distance::Real = NIF_LENS_DISTANCE, coarse::Real = 0.5,
                                                     resolution::Real = 0.02)
    ts = collect(Float64, times)
    issorted(ts) || throw(ArgumentError("times must be sorted"))
    names = Symbol[c.name for c in colors]
    _check_color_keys(ratio, names, "ratio")
    _check_color_keys(prepulse, names, "prepulse"; complete = false)
    key = _beam_key(by)
    vb = Beam[_value_beam(b) for b in beams]
    pre(b) = Tuple(haskey(prepulse, n) ? _schedule_of(prepulse[n], key, b) : nothing for n in names)
    fr = _fractions_of(ratio, names, key, vb)
    sigs, groups = _schedule_groups(eachindex(vb), i -> (_schedule_of(pulses, key, vb[i]), pre(vb[i]), fr[i]))
    # geometry of the colors, for incident light only
    unit = [UnconvertedColor(c.name; footprint = c.footprint, edge = c.edge, focal_shift = c.focal_shift)
            for c in colors]
    ne, nt = nelements(t), length(ts)
    Aback = _backside_area(t)
    E = Dict(n => zeros(ne, nt) for n in names)
    B = Dict(n => zeros(length(t.components), nt) for n in names)
    skipped = String[]
    for (sig, idx) in zip(sigs, groups)
        pulse, pres, fracs = sig
        d = unconverted_deposition(t, [remake(vb[i]; power = 1.0) for i in idx];
                                   ratio = Dict(n => 1.0 for n in names), colors = unit, offsets,
                                   lens_distance, coarse, resolution)
        append!(skipped, d.skipped)
        for (j, n) in enumerate(names)
            s = _color_power(fracs[j], pulse, pres[j], ts, n)
            Ie, Ib = _unit_irradiance(t, d[n], Aback)
            _add_outer!(E[n], Ie, s)
            _add_outer!(B[n], Ib, s)
        end
    end
    absorption = Dict(c.name => c.absorption for c in colors)
    return IrradianceHistory(t, ts, names, E, B, absorption, skipped)
end

# The unconverted light of a `drive_history`, its `unconverted` argument. Per color it is
# a fraction `ratio` of the 3ω power of each beam, plus an optional `prepulse` power per
# beam. The fraction is a number or a function of time for every beam, one per beam, or
# keyed by `key(beam)`. The prepulse is, per color, one function of time or keyed by
# `key(beam)`. The light is linear in the power of each beam. Beams whose powers have been
# equal at every step so far, on the same fraction and prepulse schedules, therefore
# share one `unconverted_deposition` at unit power. Such a set is an atom of a partition
# of the beams, and it is split when the powers of its beams diverge. The absorbed power
# of a color at a step is the sum over atoms of that pattern times the power of the atom
# in the color, `ratio · P₃ω + prepulse`. For beams on shared pulses the atoms are the
# (pulse, fraction, prepulse) schedules, as in `unconverted_history`.
mutable struct _UnconvertedLight
    ratio::Any
    prepulse::Any                  # `nothing` when there is no prepulse
    key::Any
    colors::Vector{UnconvertedColor}
    kw::NamedTuple
    target::Union{Nothing,Target}
    geometry::Vector{Tuple{String,Vec3,SVector{3,Float64},Float64,String}}
    atoms::Vector{Vector{Int}}
    fracs::Vector{Any}             # per atom, its fraction schedule per color
    pres::Vector{Any}              # per atom, its prepulse schedule per color (`nothing` for none)
    patterns::Vector{Any}          # per atom, its unit-power deposition (`nothing` until computed)
    skipped::Vector{String}
end

function _UnconvertedLight(ratio, prepulse, by, colors, kw)
    cs = collect(UnconvertedColor, colors)
    names = Symbol[c.name for c in cs]
    _check_color_keys(ratio, names, "unconverted")
    prepulse === nothing || _check_color_keys(prepulse, names, "prepulse"; complete = false)
    bad = filter(k -> !(k in (:offsets, :lens_distance, :coarse, :resolution)), collect(keys(kw)))
    isempty(bad) || throw(ArgumentError("unconverted_kw takes unconverted_deposition's offsets, " *
                                        "lens_distance, coarse and resolution, not " * join(bad, ", ")))
    return _UnconvertedLight(ratio, prepulse, _beam_key(by), cs, NamedTuple(kw), nothing,
                             Tuple{String,Vec3,SVector{3,Float64},Float64,String}[], Vector{Int}[], Any[],
                             Any[], Any[], String[])
end

# Returns per-color depositions with the colors with which they were deposited, in the
# form of the `unconverted` argument of `solve`.
_with_colors(d::AbstractDict, colors) = (u = Dict{Symbol,Any}(d); u[:colors] = colors; u)

# Returns the prepulse schedule of a beam per color, as a tuple with `nothing` for a color
# without a prepulse.
_prepulse_of(u::_UnconvertedLight, b) =
    Tuple(haskey(u.prepulse, c.name) ? _schedule_of(u.prepulse[c.name], u.key, b) : nothing for c in u.colors)

# Returns the absorbed unconverted light of each color (a `Deposition`) on target `t` from
# `beams` at time τ. The beams have their 3ω setup powers, before any laser–plasma model.
function _unconverted_at!(u::_UnconvertedLight, t::Target, beams, τ)
    vb = Beam{Float64}[_value_beam(b) for b in beams]
    geo = [(b.name, b.port, SVector{3,Float64}(b.aim), b.defocus, b.quad) for b in vb]
    if isempty(u.atoms)
        u.geometry = geo
        fr = _fractions_of(u.ratio, Symbol[c.name for c in u.colors], u.key, vb)
        if u.prepulse === nothing
            sigs, groups = _schedule_groups(eachindex(vb), i -> fr[i])
            u.atoms, u.fracs, u.pres = groups, Any[sigs...], Any[nothing for _ in groups]
        else
            sigs, groups = _schedule_groups(eachindex(vb), i -> (fr[i], _prepulse_of(u, vb[i])))
            u.atoms, u.fracs, u.pres = groups, Any[first(x) for x in sigs], Any[last(x) for x in sigs]
        end
        u.patterns = Any[nothing for _ in groups]
    elseif !isequal(geo, u.geometry)
        throw(ArgumentError("unconverted light in a drive history needs the beams' pointing fixed in " *
                            "time (their names, ports, aims and defocus), which changed at $τ ns"))
    end
    t === u.target || (u.target = t; fill!(u.patterns, nothing))
    p = Float64[b.power for b in vb]
    atoms, fracs, pres, patterns = Vector{Int}[], Any[], Any[], Any[]
    for (a, f, q, d) in zip(u.atoms, u.fracs, u.pres, u.patterns)
        _, parts = _schedule_groups(a, i -> p[i])
        if length(parts) == 1
            push!(atoms, a)
            push!(fracs, f)
            push!(pres, q)
            push!(patterns, d)
        else
            append!(atoms, [a[g] for g in parts])
            append!(fracs, fill(f, length(parts)))
            append!(pres, fill(q, length(parts)))
            append!(patterns, fill(nothing, length(parts)))
        end
    end
    u.atoms, u.fracs, u.pres, u.patterns = atoms, fracs, pres, patterns
    n, nc = nelements(t), length(t.components)
    E = Dict(c.name => zeros(n) for c in u.colors)
    L = Dict(c.name => zeros(nc) for c in u.colors)
    per_beam = Dict(c.name => zeros(length(vb)) for c in u.colors)
    delivered = Dict(c.name => 0.0 for c in u.colors)
    for (k, a) in enumerate(atoms)
        s = p[a[1]]
        # the fraction of the atom per color, and the power of each color per beam of the atom (TW)
        r = Float64[_fraction_at(f, τ) for f in fracs[k]]
        all(>=(0), r) || throw(ArgumentError("the unconverted fractions must not be negative (at $τ ns)"))
        if pres[k] === nothing
            (s > 0 && any(>(0), r)) || continue
            fs = Float64[r[j] * s for j in eachindex(u.colors)]
        else
            fs = Float64[r[j] * s + (q === nothing ? 0.0 : Float64(q(τ))) for (j, q) in enumerate(pres[k])]
            all(>=(0), fs) || throw(ArgumentError("the unconverted power must not be negative " *
                                                  "(ratio and prepulse, at $τ ns)"))
            any(>(0), fs) || continue
        end
        if patterns[k] === nothing
            d = unconverted_deposition(t, [remake(vb[i]; power = 1.0) for i in a];
                                       ratio = Dict(c.name => 1.0 for c in u.colors), colors = u.colors, u.kw...)
            union!(u.skipped, d.skipped)
            patterns[k] = d
        end
        d = patterns[k]
        for (j, c) in enumerate(u.colors)
            f = fs[j]
            f == 0 && continue
            dc = d[c.name]::Deposition{Float64}
            E[c.name] .+= f .* dc.element
            L[c.name] .+= f .* dc.backside
            per_beam[c.name][a] .+= f .* dc.per_beam
            delivered[c.name] += f * dc.delivered
        end
    end
    return Dict{Symbol,Deposition{Float64}}(c.name => Deposition(E[c.name], L[c.name], per_beam[c.name],
                                                                 delivered[c.name]) for c in u.colors)
end

"""
    irradiance_history(op, beams, pulses, times; by = :group) -> IrradianceHistory

Irradiance of the 3ω light (color `:omega3`) on the target of the
[`DepositionOperator`](@ref) `op` at each of `times` (ns), with the beams following
`pulses`. `pulses` is one function of time for all beams, or keyed by `by` as in
[`with_pulses`](@ref). With the geometry fixed, deposition is linear in the power of
each beam. The history therefore costs one [`deposit`](@ref) per distinct pulse, with
the beams of that pulse at unit power and the others at zero, scaled by the pulse at
each time. The irradiance is the incident light (every ray deposits all its power),
before any laser–plasma model. [`preplasma_check`](@ref) uses it for the time at which
the backlighter beams arrive on their foil.
"""
Base.@nospecializeinfer function irradiance_history(op::DepositionOperator, beams::AbstractVector{<:Beam},
                                                    @nospecialize(pulses), @nospecialize(times);
                                                    @nospecialize(by = :group))
    t = op.target
    ts = collect(Float64, times)
    issorted(ts) || throw(ArgumentError("times must be sorted"))
    key = _beam_key(by)
    vb = Beam[_value_beam(b) for b in beams]
    sigs, groups = _schedule_groups(vb, b -> _schedule_of(pulses, key, b))
    Aback = _backside_area(t)
    E = zeros(nelements(t), length(ts))
    B = zeros(length(t.components), length(ts))
    for (pulse, idx) in zip(sigs, groups)
        d = deposit(op, [remake(b; power = i in idx ? 1.0 : 0.0) for (i, b) in enumerate(vb)])
        s = [Float64(pulse(τ)) for τ in ts]
        all(>=(0), s) || throw(ArgumentError("the pulse power must not be negative"))
        Ie, Ib = _unit_irradiance(t, d, Aback)
        _add_outer!(E, Ie, s)
        _add_outer!(B, Ib, s)
    end
    return IrradianceHistory(t, ts, [:omega3], Dict(:omega3 => E), Dict(:omega3 => B),
                             Dict(:omega3 => 1.0), String[])
end

# Returns the first time at which a linearly interpolated series `I` (at times `ts`)
# reaches `thr`, or Inf if it never does.
function _first_reach(ts, I, thr)
    k = findfirst(>=(thr), I)
    k === nothing && return Inf
    k == 1 && return ts[1]
    return ts[k-1] + (thr - I[k-1]) / (I[k] - I[k-1]) * (ts[k] - ts[k-1])
end

# Returns the irradiance of `colors` summed on the elements and back sides of the surface
# components of `parts`, as rows, and per row the pair (component index, element), with
# element 0 for the back side.
function _part_series(h::IrradianceHistory, parts, colors)
    t = h.target
    cols = colors === nothing ? h.colors : Symbol[colors...]
    all(in(h.colors), cols) || throw(ArgumentError("the history has the colors $(h.colors)"))
    cis = parts === nothing ? collect(eachindex(t.components)) :
          unique!(reduce(vcat, (part_pieces(t, Symbol(p)) for p in parts); init = Int[]))
    labels = Tuple{Int,Int}[]
    for ci in cis
        t.components[ci].kind === :surface || continue
        append!(labels, (ci, e) for e in t.comp_elems[ci])
        t.components[ci].two_sided || push!(labels, (ci, 0))
    end
    S = zeros(length(labels), length(h.time))
    for c in cols, (i, (ci, e)) in enumerate(labels)
        @views S[i, :] .+= e > 0 ? h.element[c][e, :] : h.backside[c][ci, :]
    end
    return S, labels
end

"""
    threshold_crossing(history, threshold; parts = nothing, colors = nothing)
        -> NamedTuple (time, part, element)

The first time (ns) at which the irradiance of `history` (an
[`IrradianceHistory`](@ref)), summed over `colors`, reaches `threshold` (W/cm²) on any
element or one-sided back side of the surface components of `parts`. By default
`colors` are all the colors of the history and `parts` are all parts, and the history is
linear between its times. `part` is the part where the threshold is reached, and
`element` the element (0 for a back side, whose irradiance is the component mean). When
the threshold is never reached, `time = Inf`, `part = :none` and `element = 0`.
"""
function threshold_crossing(h::IrradianceHistory, threshold::Real; parts = nothing, colors = nothing)
    threshold > 0 || throw(ArgumentError("threshold must be positive"))
    S, labels = _part_series(h, parts, colors)
    best, at = Inf, 0
    for i in eachindex(labels)
        τ = _first_reach(h.time, view(S, i, :), threshold)
        τ < best && (best = τ; at = i)
    end
    at == 0 && return (; time = Inf, part = :none, element = 0)
    ci, e = labels[at]
    return (; time = best, part = h.target.components[ci].part, element = e)
end

"""
    specular_check(history; threshold = 1e8, parts = nothing, colors = nothing)
        -> Vector{NamedTuple}

Whether each surface part stays below the irradiance `threshold` (W/cm²) at which it is
taken to stop being a specular reflector, under the light of `history` summed over
`colors`. `history` is an [`IrradianceHistory`](@ref), e.g. of
[`unconverted_history`](@ref). The result has one row per part of `parts` (by default
every surface part), in target order, with the fields `part`, `peak` (the highest
irradiance on any of its elements or back sides, W/cm²), `peak_time` (ns), `crossing`
and `specular`. `crossing` is the first time at which the part reaches `threshold`
([`threshold_crossing`](@ref), `Inf` if never), and `specular` is `crossing == Inf`.

A surface below the threshold stays a cold mirror, and the NIF curvature rule disperses
its specular 1ω reflection. The rule requires a curvature radius ≤ 2 mm in one
dimension, ≤ 5 mm in two, or dimpling (User Guide §5.6.1). Above the threshold a plasma
forms, and the cold-mirror reflection no longer applies. The default 10⁸ W/cm² is an
*engineering default* with low confidence. It is the general guideline of the User
Guide for the onset of weakly ionized plasma on a metal surface (§4.4, citing Spaeth et
al., Fusion Sci. Technol. 69, 25 (2016)), taken as the end of specular reflection. No
facility criterion for specularity sets it. It applies only to metals, and silicon and
plastics need separate values. A fluence-driven change of the surface (melting,
roughening) is not modeled. See the research notes, "Unconverted light over time".
"""
function specular_check(h::IrradianceHistory; threshold::Real = 1e8, parts = nothing, colors = nothing)
    t = h.target
    ps = parts === nothing ? unique!([c.part for c in t.components if c.kind === :surface]) :
         Symbol[parts...]
    return map(ps) do p
        S, _ = _part_series(h, (p,), colors)
        peak, k = isempty(S) ? (0.0, 0) : findmax(vec(maximum(S; dims = 1)))
        crossing = threshold_crossing(h, threshold; parts = (p,), colors).time
        (; part = p, peak, peak_time = k == 0 ? NaN : h.time[k], crossing, specular = crossing == Inf)
    end
end

"""
    preplasma_check(history, arrival; foil = :backlighter, threshold = 1e8, colors = nothing)
        -> NamedTuple (onset, arrival, lead, preplasma)

Whether the foil `foil` forms a pre-plasma before its backlighter beams arrive. `onset`
is the first time (ns) at which the irradiance of `history`, the unconverted light of
[`unconverted_history`](@ref) summed over `colors`, reaches `threshold` (W/cm²) on the
foil ([`threshold_crossing`](@ref)). `arrival` is a time (ns), or the 3ω
[`irradiance_history`](@ref) of the backlighter beams, whose first crossing of the same
threshold on the foil is then the arrival. `lead = arrival - onset` (ns; `-Inf` when the
unconverted light never reaches the threshold), and `preplasma = onset < arrival`.

The default 10⁸ W/cm² is the general guideline of the NIF User Guide for the onset of
weakly ionized plasma on a metal surface (§4.4, citing Spaeth et al., Fusion Sci.
Technol. 69, 25 (2016)). The guide puts the 1ω prepulse above it 3–5 ns before the
nominal start of the pulse. It is a guideline for metals, and no facility criterion for
foils sets it. A foil of another material needs a separate value.
"""
preplasma_check(h::IrradianceHistory, arrival::IrradianceHistory; foil::Symbol = :backlighter,
                threshold::Real = 1e8, colors = nothing) =
    preplasma_check(h, threshold_crossing(arrival, threshold; parts = (foil,)).time; foil, threshold, colors)

function preplasma_check(h::IrradianceHistory, arrival::Real; foil::Symbol = :backlighter,
                         threshold::Real = 1e8, colors = nothing)
    onset = threshold_crossing(h, threshold; parts = (foil,), colors).time
    a = Float64(arrival)
    return (; onset, arrival = a, lead = onset == Inf ? -Inf : a - onset, preplasma = onset < a)
end
