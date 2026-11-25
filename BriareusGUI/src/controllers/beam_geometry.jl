# Beam geometry for the viewer. The axis of each beam is traced from a start point outside
# the target to the first opaque surface it meets, where its light is deposited. The volume
# of the beam is the best-focus contour carried along the f-cone, as VisRad draws beams,
# with each generator line ending where it meets a surface. The code is pure geometry and
# uses no GL.

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
    # Each generator is x = focus + w + s (k - cone u) before focus and uses (k + cone u)
    # after it.
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

Shortest distance (mm) between a pick ray and the drawn axis of the beam, from start to
stop.
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

"""
    unconverted_footprint(harmonic; mode = :nominal) -> NamedTuple

The unconverted-light footprint drawn for `harmonic` 1 (1ω) or 2 (2ω) at the 3ω
best-focus plane: `edge` (mm from the 3ω axis to its near edge), `size` (side, mm),
`focus` (its focus, mm beyond the 3ω focus) and `lens` (the lens distance, mm before the
3ω focus). With `mode = :nominal` these are the footprints of `nif_unconverted` (1ω: 22 mm
at 4.8 mm; 2ω: 12 mm at 2.9 mm), the geometry `unconverted_deposition` deposits with. For
1ω, `mode = :talis_10pct` gives the 10% contour that `talis_check` measures to (near edge
4.45 mm) and `mode = :talis_worst` that contour grown by the default error set and
required clearance, both from `Briareus.omega1_contour`. The 2ω footprint has only the
nominal mode.
"""
function unconverted_footprint(h::Integer; mode::Symbol = :nominal)
    1 <= h <= 2 || throw(ArgumentError("unconverted light is harmonic 1 or 2"))
    if h == 1
        c = Briareus.omega1_contour(; mode)
        return (edge = c.edge, size = c.size, focus = c.focal_shift, lens = Briareus.NIF_LENS_DISTANCE)
    end
    mode === :nominal || throw(ArgumentError("the 2ω footprint has only the :nominal mode"))
    uc = nif_unconverted()[h]
    return (edge = uc.edge, size = uc.footprint, focus = uc.focal_shift, lens = Briareus.NIF_LENS_DISTANCE)
end

const _NIF_TABLE = Ref{Any}(nothing)
function _nif_table()
    _NIF_TABLE[] === nothing && (_NIF_TABLE[] = Briareus.GUIAPI.nif_beam_table())
    return _NIF_TABLE[]
end

"""
    has_unconverted(b::Beam) -> Bool

Whether beam `b` has an unconverted-light footprint to draw: a NIF beam or quad (by its
quad name), as for `unconverted_deposition` and `talis_check`.s 1ω rules.
"""
has_unconverted(b::Beam) = Briareus.GUIAPI.is_nif(b)

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
    bl = Briareus.GUIAPI.omega1_beamlets(b, nothing, _nif_table())
    return bl === nothing ? Tuple{Vec3,Vec3,Vec3}[] : Tuple{Vec3,Vec3,Vec3}[Tuple(Vec3.(x)) for x in bl]
end

"""
    unconverted_square(beamlet, harmonic, s; mode = :nominal) -> NTuple{4,Vec3}

Corners of a beamlet's ([`unconverted_beamlets`](@ref)) unconverted-light footprint of
`harmonic` at axial position `s` (mm past its 3ω best focus): the square of
[`unconverted_footprint`](@ref) at `s = 0`, its center on the chief ray from the lens
center and its side shrinking linearly to the unconverted focus, as
`unconverted_deposition` models it.
"""
function unconverted_square((k, f, ô), h::Integer, s::Real; mode::Symbol = :nominal)
    fp = unconverted_footprint(h; mode)
    c = (fp.edge + fp.size / 2) * (1 + s / fp.lens)
    a = fp.size / 2 * abs(1 - s / fp.focus)
    q̂ = cross(k, ô)
    m = f + s * k + c * ô
    return (m - a * ô - a * q̂, m + a * ô - a * q̂, m + a * ô + a * q̂, m - a * ô + a * q̂)
end

"""
    unconverted_light(b::Beam, harmonic; length, mode = :nominal) -> (; vertices, faces, outlines, axes)

What the 3D view draws of beam `b`'s unconverted light of `harmonic` (1 or 2), per
beamlet: the converging square beam (`vertices`, `faces`: index triples) from `length`
mm before the aim to `length` mm past it, drawn through the target, since unconverted
light mostly passes it and only the deposition resolves its shadows, the
footprint square at
the 3ω best-focus plane (`outlines`, segments) and the footprint's center line
(`axes`, segments). Empty for a beam without a footprint. `mode` is the contour of the 1ω
footprint ([`unconverted_footprint`](@ref)); the 2ω footprint is always the nominal one.
"""
function unconverted_light(b::Beam, h::Integer; length::Real, mode::Symbol = :nominal)
    mode = h == 1 ? mode : :nominal
    verts, faces = Vec3[], NTuple{3,Int}[]
    outlines, axes = NTuple{2,Vec3}[], NTuple{2,Vec3}[]
    for bl in unconverted_beamlets(b)
        s0, s1 = -(length + b.defocus), length - b.defocus
        A, B = unconverted_square(bl, h, s0; mode), unconverted_square(bl, h, s1; mode)
        o = Base.length(verts)
        append!(verts, A)
        append!(verts, B)
        for i in 1:4
            j = mod1(i + 1, 4)
            push!(faces, (o + i, o + j, o + 4 + j), (o + i, o + 4 + j, o + 4 + i))
        end
        S = unconverted_square(bl, h, 0.0; mode)
        for i in 1:4
            push!(outlines, (S[i], S[mod1(i + 1, 4)]))
        end
        push!(axes, (sum(A) / 4, sum(B) / 4))
    end
    return (; vertices = verts, faces, outlines, axes)
end

# ---------------------------------------------------------------- 1ω closest approach

"""
    closest_approach_surface(approach::Omega1Approach; span = 6.0)
        -> (; outline, vertices, faces, caliper, budget)

What the 3D view draws of a 1ω closest approach (`Briareus.omega1_closest_approach`): the
side or sides of the 1ω beamlet nearest the hardware point, and nothing else of the beam.

- `outline` — segments around those sides of the 10% contour at the focal plane of the
  approach, cut to `span` mm either side of the station of the hardware point
- `vertices`, `faces` — the same sides moved out by the transverse errors, the rule
  tolerance at the hardware point and the required clearance: the plane that hardware
  must stay outside in the worst case (index triples)
- `caliper` — the segment from the hardware point to the closest point of the contour,
  in the plane across the beamlet in which the rule measures
- `budget` — lines of text from the nominal distance to the margin, one per term
"""
function closest_approach_surface(a::Omega1Approach; span::Real = 6.0)
    span > 0 || throw(ArgumentError("the span must be positive"))
    k, f, ô = a.beamlet
    q̂ = cross(k, ô)
    f = f + a.focus_shift * k
    c = a.contour
    b = a.budget
    grow = b.pointing + b.target_offset + b.hardware_tolerance + b.rotation_tol + b.required
    # each side as its outward direction and the direction along it
    sides = ((0x01, -ô, q̂), (0x02, ô, q̂), (0x04, -q̂, ô), (0x08, q̂, ô))
    corner(s, out, along, sgn, g) = begin
        h = c.h0 * abs(1 - s / c.focal_shift)
        f + s * k + c.c0 * (1 + s / c.lens_distance) * ô + (h + g) * out + sgn * (h + g) * along
    end
    outline, verts, faces = NTuple{2,Vec3}[], Vec3[], NTuple{3,Int}[]
    s0, s1 = a.s_star - span, a.s_star + span
    for (bit, out, along) in sides
        a.face_mask & bit == 0 && continue
        p = (corner(s0, out, along, -1, 0.0), corner(s0, out, along, 1, 0.0),
             corner(s1, out, along, 1, 0.0), corner(s1, out, along, -1, 0.0))
        for i in 1:4
            push!(outline, (p[i], p[mod1(i + 1, 4)]))
        end
        o = length(verts)
        push!(verts, corner(s0, out, along, -1, grow), corner(s0, out, along, 1, grow),
              corner(s1, out, along, 1, grow), corner(s1, out, along, -1, grow))
        push!(faces, (o + 1, o + 2, o + 3), (o + 1, o + 3, o + 4))
    end
    line(label, v) = @sprintf("%-22s %+8.3f mm", label, v)
    budget = [line("nominal", a.nominal), line("focus ($(a.focus_shift) mm)", -b.focus_delta),
              line("pointing", -b.pointing), line("target position", -b.target_offset),
              line("hardware tolerance", -b.hardware_tolerance), line("rule tolerance", -b.rotation_tol),
              line("required clearance", -b.required), line("margin", a.margin)]
    return (; outline, vertices = verts, faces, caliper = (a.x_hw, a.x_beam), budget)
end

"""
    reflection_cone(x::Omega1Reflection; length = 10.0, required = 20.0, tolerance = 0.5, n = 32)
        -> (; incident, reflected, port, vertices, faces, half_angle)

What the 3D view draws of the specular 1ω reflection behind a management row
(`Briareus.omega1_reflection`): the `incident` ray into the reflecting point, the
`reflected` axis and the direction to the threatened open `port` (segments of `length`
mm), and the cone about the reflected axis that an open port must stay outside
(`vertices`, `faces`), of half-angle `required` plus twice the rule `tolerance` plus the
normal spread of the facet.
"""
function reflection_cone(x::Omega1Reflection; length::Real = 10.0, required::Real = 20.0,
                         tolerance::Real = 0.5, n::Integer = 32)
    (length > 0 && n >= 3) || throw(ArgumentError("need a positive length and n ≥ 3"))
    half = required + 2 * tolerance + x.spread
    0 < half < 90 || throw(ArgumentError("the cone half-angle must lie between 0 and 90 degrees"))
    p = x.point
    r = x.reflected
    u = normalize(cross(r, abs(r[3]) < 0.9 ? Vec3(0, 0, 1) : Vec3(1, 0, 0)))
    v = cross(r, u)
    verts = Vec3[p]
    for i in 0:n-1
        φ = 2π * i / n
        push!(verts, p + length * (cosd(half) * r + sind(half) * (cos(φ) * u + sin(φ) * v)))
    end
    faces = [(1, 1 + i, 1 + mod1(i + 1, n)) for i in 1:n]
    return (; incident = (p - length * x.incident, p), reflected = (p, p + length * r),
            port = (p, p + length * x.port), vertices = verts, faces, half_angle = half)
end

"""
    daylight_view(approach::Omega1Approach; zoom = 4.0) -> (; center, direction, up, width)

The orthographic view down the beamlet of a 1ω closest approach: the camera looks along
the beamlet, so the 1ω square is seen edge-on in the plane in which the rule measures, and
the gap seen between hardware and the contour is the transverse distance. `center` is
midway along the caliper, `up` the direction across the nearest side, and `width` the
extent of the view (mm), `zoom` times the length of the caliper and at least 1 mm.
"""
function daylight_view(a::Omega1Approach; zoom::Real = 4.0)
    k, _, ô = a.beamlet
    up = a.face_mask & 0x03 != 0 ? cross(k, ô) : ô
    return (; center = (a.x_hw + a.x_beam) / 2, direction = k, up,
            width = max(zoom * norm(a.x_hw - a.x_beam), 1.0))
end

"""
    grazing_view(approach::Omega1Approach; zoom = 4.0) -> (; center, direction, up, width)

The orthographic view along the nearest side of the 1ω beamlet of a closest approach and
across the beam: the side is seen edge-on with the beam running up the view, so the
profile of the hardware is seen against the contour along the beam. `direction` lies in
the side, at right angles to the caliper and to the beamlet.
"""
function grazing_view(a::Omega1Approach; zoom::Real = 4.0)
    k, _, ô = a.beamlet
    along = a.face_mask & 0x03 != 0 ? cross(k, ô) : ô
    return (; center = (a.x_hw + a.x_beam) / 2, direction = along, up = k,
            width = max(zoom * norm(a.x_hw - a.x_beam), 1.0))
end
