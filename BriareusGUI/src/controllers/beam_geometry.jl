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
