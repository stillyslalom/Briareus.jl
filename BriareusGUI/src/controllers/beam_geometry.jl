# Beam geometry for the viewer: each beam's axis traced from a start point outside the
# target to the first opaque surface it meets (where its light lands), and the beam's
# volume (the best-focus contour carried along the f-cone, as VisRad draws beams),
# each generator line ending where it meets a surface. Pure geometry, no GL.

"""
Ways the viewer draws beams: `:off`, `:rays` (each beam's axis from its start to where
it lands, the aim point marked) and `:volumes` (the beam envelope, cut where its edge
rays meet the target).
"""
const BEAM_STYLES = (:off, :rays, :volumes)

"""
    BeamPath

A beam's axis in the viewer: `start` (`length` mm before the aim point, against the
beam direction), `aim`, `focus` (best focus), `stop` (where the axis meets the first
opaque surface, or `length` mm past the aim when it meets none), `component` (the
landing component's name, `:missed` without one).
"""
struct BeamPath
    start::Vec3
    aim::Vec3
    focus::Vec3
    stop::Vec3
    component::Symbol
end

_hit(t::Target, o, d, tmax = Inf) = Briareus.closest_hit(t.bvh, Vec3(o), Vec3(d), tmax)

"""
    beam_path(t::Target, b::Beam; length = 2 t.scale) -> BeamPath
"""
function beam_path(t::Target, b::Beam; length::Real = 2 * t.scale)
    k = direction(b)
    aim = Vec3(b.aim)
    o = aim - length * k
    τ, tri = _hit(t, o, k)
    focus = aim + b.defocus * k
    tri == 0 && return BeamPath(o, aim, focus, aim + length * k, :missed)
    return BeamPath(o, aim, focus, o + τ * k, t.components[t.tri_comp[tri]].name)
end

beam_paths(t::Target, beams; kw...) = [beam_path(t, b; kw...) for b in beams]

"""
    beam_volume(t::Target, b::Beam; length = 2 t.scale, peak_fraction = 0.05,
                nring = 24, nseg = 6) -> (vertices, faces)

Triangle mesh of a beam's envelope: the best-focus contour at `peak_fraction` of the
peak, grown along the beam by the f-cone (narrowing to best focus, widening after),
from `length` mm before the aim point. Each of `nring` generator lines ends where it
meets an opaque surface (or `length` mm past the aim); `nseg` rings sample each side
of best focus. `faces` are index triples into `vertices`.
"""
function beam_volume(t::Target, b::Beam; length::Real = 2 * t.scale,
                     peak_fraction::Real = 0.05, nring::Integer = 24, nseg::Integer = 6)
    k = direction(b)
    e1, e2 = Briareus.spot_axes(b)
    ρ = envelope_radius(b.spot, b.power; peak_fraction)
    A, B = ρ * b.spot.a, ρ * b.spot.b
    cone = 1 / (2b.fnumber)
    focus = Vec3(b.aim) + b.defocus * k
    s0 = -(length + b.defocus)                   # start: `length` before the aim
    s1 = length - b.defocus                      # end without a hit
    ψs = range(0, 2π; length = nring + 1)[1:end-1]
    # each generator: before focus x = focus + w + s (k - cone u), after: (k + cone u)
    stops = map(ψs) do ψ
        u = cos(ψ) * e1 + sin(ψ) * e2
        w = A * cos(ψ) * e1 + B * sin(ψ) * e2
        d1 = k - cone * u
        o1 = focus + w + s0 * d1
        τ, tri = _hit(t, o1, d1 / norm(d1), -s0 * norm(d1))
        tri != 0 && return s0 + τ / norm(d1)
        d2 = k + cone * u
        τ, tri = _hit(t, focus + w, d2 / norm(d2), s1 * norm(d2))
        return tri != 0 ? τ / norm(d2) : s1
    end
    smax = maximum(stops)
    ss = vcat(range(s0, 0; length = nseg + 1), range(0, max(smax, 0); length = nseg + 1)[2:end])
    smax < 0 && (ss = collect(range(s0, smax; length = nseg + 1)))
    verts = Vec3[]
    for (ψ, st) in zip(ψs, stops), s in ss
        u = cos(ψ) * e1 + sin(ψ) * e2
        w = A * cos(ψ) * e1 + B * sin(ψ) * e2
        sc = min(s, st)
        push!(verts, focus + w + sc * k + cone * abs(sc) * u)
    end
    nsv = Base.length(ss)
    idx(i, j) = (mod1(i, nring) - 1) * nsv + j
    faces = NTuple{3,Int}[]
    for i in 1:nring, j in 1:nsv-1
        a, b_, c, d = idx(i, j), idx(i, j + 1), idx(i + 1, j), idx(i + 1, j + 1)
        push!(faces, (a, b_, d), (a, d, c))
    end
    return verts, faces
end

"""
    beam_hit_distance(path::BeamPath, origin, direction) -> Float64

Shortest distance (mm) between a pick ray and the beam's drawn axis (start → stop).
"""
function beam_hit_distance(p::BeamPath, o, d)
    a, b = Vec3(p.start), Vec3(p.stop)
    u = b - a
    d = Vec3(d) / norm(d)
    w0 = a - Vec3(o)
    uu, ud, dd = dot(u, u), dot(u, d), 1.0
    uw, dw = dot(u, w0), dot(d, w0)
    den = uu * dd - ud^2
    s = den > 1e-12 ? clamp((ud * dw - dd * uw) / den, 0.0, 1.0) : 0.0
    τ = max(dot(a + s * u - Vec3(o), d), 0.0)
    return norm(a + s * u - (Vec3(o) + τ * d))
end

# ---------------------------------------------------------------- unconverted light

# The lens distance (mm upstream of the 3ω focus) of `unconverted_deposition` and
# `talis_check` (their default).
const _LENS_DISTANCE = 7700.0

"""
    unconverted_footprint(harmonic) -> NamedTuple

The unconverted-light footprint drawn for `harmonic` 1 (1ω) or 2 (2ω) at the 3ω
best-focus plane: `edge` (mm from the 3ω axis to its near edge), `size` (side, mm),
`focus` (its focus, mm beyond the 3ω focus) and `lens` (the lens distance, mm before the
3ω focus). These are `nif_unconverted`'s nominal footprints (1ω: 22 mm at 4.8 mm; 2ω:
12 mm at 2.9 mm), the geometry `unconverted_deposition` deposits with; `talis_check`'s
1ω rule uses the same square grown to its 10 % contour (near edge 4.45 mm).
"""
function unconverted_footprint(h::Integer)
    1 <= h <= 2 || throw(ArgumentError("unconverted light is harmonic 1 or 2"))
    uc = nif_unconverted()[h]
    return (edge = uc.edge, size = uc.footprint, focus = uc.focal_shift, lens = _LENS_DISTANCE)
end

const _NIF_TABLE = Ref{Any}(nothing)
function _nif_table()
    _NIF_TABLE[] === nothing && (_NIF_TABLE[] = Briareus._nif_beam_table())
    return _NIF_TABLE[]
end

"""
    has_unconverted(b::Beam) -> Bool

Whether beam `b` has an unconverted-light footprint to draw: a NIF beam or quad (by its
quad name), as for `unconverted_deposition` and `talis_check`.s 1ω rules.
"""
has_unconverted(b::Beam) = Briareus._is_nif(b)

"""
    unconverted_beamlets(b::Beam) -> Vector

The beamlets of beam `b` whose unconverted light `unconverted_deposition` deposits and
`talis_check` checks, as (k, 3ω
focus, ô): the beamlet's direction, its 3ω best focus and the direction its footprint
is displaced in (a NIF quad-level beam gives its four beamlets; empty for a beam
without a footprint, [`has_unconverted`](@ref)).
"""
function unconverted_beamlets(b::Beam)
    has_unconverted(b) || return Tuple{Vec3,Vec3,Vec3}[]
    bl = Briareus._omega1_beamlets(b, nothing, _nif_table())
    return bl === nothing ? Tuple{Vec3,Vec3,Vec3}[] : Tuple{Vec3,Vec3,Vec3}[Tuple(Vec3.(x)) for x in bl]
end

"""
    unconverted_square(beamlet, harmonic, s) -> NTuple{4,Vec3}

Corners of a beamlet's ([`unconverted_beamlets`](@ref)) unconverted-light footprint of
`harmonic` at axial position `s` (mm past its 3ω best focus): the square of
[`unconverted_footprint`](@ref) at `s = 0`, its center on the chief ray from the lens
center and its side shrinking linearly to the unconverted focus, as
`unconverted_deposition` models it.
"""
function unconverted_square((k, f, ô), h::Integer, s::Real)
    fp = unconverted_footprint(h)
    c = (fp.edge + fp.size / 2) * (1 + s / fp.lens)
    a = fp.size / 2 * abs(1 - s / fp.focus)
    q̂ = cross(k, ô)
    m = f + s * k + c * ô
    return (m - a * ô - a * q̂, m + a * ô - a * q̂, m + a * ô + a * q̂, m - a * ô + a * q̂)
end

"""
    unconverted_light(b::Beam, harmonic; length) -> (; vertices, faces, outlines, axes)

What the 3D view draws of beam `b`'s unconverted light of `harmonic` (1 or 2), per
beamlet: the converging square beam (`vertices`, `faces`: index triples) from `length`
mm before the aim to `length` mm past it, not stopped by the target (unconverted light
mostly passes it; the deposition resolves its shadows, the drawing does not), the
footprint square at
the 3ω best-focus plane (`outlines`, segments) and the footprint's center line
(`axes`, segments). Empty for a beam without a footprint.
"""
function unconverted_light(b::Beam, h::Integer; length::Real)
    verts, faces = Vec3[], NTuple{3,Int}[]
    outlines, axes = NTuple{2,Vec3}[], NTuple{2,Vec3}[]
    for bl in unconverted_beamlets(b)
        s0, s1 = -(length + b.defocus), length - b.defocus
        A, B = unconverted_square(bl, h, s0), unconverted_square(bl, h, s1)
        o = Base.length(verts)
        append!(verts, A)
        append!(verts, B)
        for i in 1:4
            j = mod1(i + 1, 4)
            push!(faces, (o + i, o + j, o + 4 + j), (o + i, o + 4 + j, o + 4 + i))
        end
        S = unconverted_square(bl, h, 0.0)
        for i in 1:4
            push!(outlines, (S[i], S[mod1(i + 1, 4)]))
        end
        push!(axes, (sum(A) / 4, sum(B) / 4))
    end
    return (; vertices = verts, faces, outlines, axes)
end
