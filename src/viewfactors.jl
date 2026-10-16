# Element-to-element view factors with occlusion.

"""
    view_factors(target; near = 3.0, level = 1) -> Matrix{Float64}

Dense view-factor matrix `F[i, j]`: the fraction of diffuse radiation leaving element
`i` that arrives at element `j`. Reciprocity `Aᵢ Fᵢⱼ = Aⱼ Fⱼᵢ` holds exactly by
construction, and rows sum to one for a closed enclosure (counting apertures); see
`closure_error`.

Flat elements use the exact differential-area-to-polygon view factor onto the other
element, integrated over the emitter: over `4^level` sub-triangle points for pairs
closer than `near` element sizes (√2 × the largest centroid-to-vertex distance, the side
of a square element), which stays accurate for elements sharing an edge such as the
wall–cap corner of a hohlraum, and from the centroid beyond that, averaged over the two
directions.

Elements of a component with a `CurvedSurface` are integrated over the curved patch:
quadrature points on the true surface with its normals and areas, point to point
between two curved elements (`4^level` sub-triangle points per triangle when near, one
per triangle out to three times that distance, the projected centroid beyond), and onto
the exact polygon of a planar partner. Elements of a `DiskSurface` are exact polar
sectors, planar receivers whose arcs are divided finely. A curved concave element also
sees itself, which gives the diagonal `Fᵢᵢ`. Flat facets standing in for a curved wall misrepresent the
exchange between nearby elements by an amount of order one relative to its size, which
matters next to a concentrated source such as a laser spot.

Visibility is a single centroid-to-centroid ray between the flat facets.
"""
function view_factors(t::Target; near = 3.0, mid = 12.0, level::Integer = 1, exact_far::Bool = true)
    N = nelements(t)
    F = zeros(N, N)
    C, Nrm, A = t.elem_centroid, t.elem_normal, t.elem_area
    curved = t.elem_curved .!= 0
    # planar elements (flat facets and polar sectors) receive through the exact
    # point-to-polygon kernel; sectors are polygons with finely divided arcs
    sector = [curved[e] && t.components[t.elem_comp[e]].surface isa DiskSurface for e in 1:N]
    planar = .!curved .| sector
    recv = map(1:N) do e
        if sector[e]
            s, sec = t.components[t.elem_comp[e]].surface, element_sector(t, e)
            return _Receiver(_sector_triangles(s, sec), _sector_loop(s, sec))
        end
        return _Receiver([(t.tri_v0[k], t.tri_v0[k] + t.tri_e1[k], t.tri_v0[k] + t.tri_e2[k])
                          for k in t.elem_tris[e]], Vec3[])
    end
    concave = curved .& .!sector
    # point-to-point quadrature on curved surfaces carries an O((h/r)²) bias that does not
    # cancel between pairs, so curved elements keep several points per triangle at every
    # distance: level + 1 near, level out to `mid` element sizes, one per triangle beyond
    qnear = [element_quadrature(t, e, curved[e] && !sector[e] ? level + 1 : level) for e in 1:N]
    qmid = [element_quadrature(t, e, curved[e] ? level : 0) for e in 1:N]
    qfar = [curved[e] ? element_quadrature(t, e, 0) : ([C[e]], [A[e]], [Nrm[e]]) for e in 1:N]
    qself = [concave[e] ? qnear[e] : nothing for e in 1:N]      # same level + 1 rule
    ε = _ray_eps(t)
    # element size: largest distance from the centroid to a vertex, so a long sliver
    # counts as large even when its area is small (√2 × half-diagonal = side for a square)
    size = [sqrt(2) * maximum(max(norm(t.tri_v0[k] - C[e]), norm(t.tri_v0[k] + t.tri_e1[k] - C[e]),
                                  norm(t.tri_v0[k] + t.tri_e2[k] - C[e])) for k in t.elem_tris[e])
            for e in 1:N]
    ctx = (; t, bvh = t.bvh, C, Nrm, A, planar, concave, recv, qnear, qmid, qfar, qself,
           size, ε, near2 = near^2, mid2 = mid^2, level, exact_far)
    # rows have i:N pairs, so they are handed out one at a time, longest first, rather
    # than in equal contiguous chunks (the first chunk would carry most of the work)
    next = Threads.Atomic{Int}(1)
    @sync for _ in 1:min(Threads.nthreads(), N)
        Threads.@spawn while (i = Threads.atomic_add!(next, 1)) <= N
            _vf_row!(F, ctx, i)
        end
    end
    return F
end

# Row i of the view-factor matrix (pairs i < j), behind a function barrier so the
# threaded loop is specialized on the concrete types in `ctx`.
function _vf_row!(F, ctx, i)
    (; t, bvh, C, Nrm, A, planar, concave, recv, qnear, qmid, qfar, qself, size, ε, near2, mid2,
       level, exact_far) = ctx
    N = length(A)
    oi = C[i] + ε * Nrm[i]
    if concave[i]
        F[i, i] = _self_kernel(qself[i]) / A[i]
    end
    for j in i+1:N
        d = C[j] - C[i]
        r2 = dot(d, d)
        r2 == 0 && continue
        # elements must face each other (centroid test)
        dot(Nrm[i], d) > 0 && dot(Nrm[j], d) < 0 || continue
        oj = C[j] + ε * Nrm[j]
        any_hit(bvh, oi, oj - oi, 1 - 1e-9) && continue
        s2 = max(size[i], size[j])^2
        tier = r2 < near2 * s2 ? :near : r2 < mid2 * s2 ? :mid : :far
        # G = Aᵢ Fᵢⱼ = Aⱼ Fⱼᵢ
        if planar[i] && planar[j]
            if tier === :near
                G = (_adaptive_kernel(t, i, recv[j]; mindepth = level + 1) +
                     _adaptive_kernel(t, j, recv[i]; mindepth = level + 1)) / 2
            elseif exact_far
                G = (_pair_kernel(qfar[i], recv[j]) + _pair_kernel(qfar[j], recv[i])) / 2
            else
                G = dot(Nrm[i], d) * -dot(Nrm[j], d) * A[i] * A[j] / (π * r2^2)
            end
        else
            q = tier === :near ? qnear : tier === :mid ? qmid : qfar
            G = !planar[i] && !planar[j] ? _point_kernel(q[i], q[j]) :
                tier === :near ? (planar[j] ? _adaptive_kernel(t, i, recv[j]; mindepth = level + 1) :
                                  _adaptive_kernel(t, j, recv[i]; mindepth = level + 1)) :
                planar[j] ? _pair_kernel(q[i], recv[j]) : _pair_kernel(q[j], recv[i])
        end
        F[i, j] = G / A[i]
        F[j, i] = G / A[j]
    end
    return nothing
end

# One-point rule of a curved element: the projected centroid, its normal and the area.
function _one_point(t::Target, e)
    p, n = surface_point(t.components[t.elem_comp[e]].surface, t.elem_centroid[e])
    return [p], [t.elem_area[e]], [t.elem_curved[e] * n]
end

# Receiving polygon of a planar element: its triangles, and for a polar sector also the
# boundary loop of those triangles (same vertices, same orientation).
struct _Receiver
    tris::Vector{NTuple{3,Vec3}}
    loop::Vector{Vec3}
end

# Boundary of `_sector_triangles(s, sector; m)`, with the same vertices and orientation:
# the outer arc, then the inner arc back (or the center). Interior edges of the triangles
# cancel pairwise, so the contour over this loop is the sum over the triangles whenever
# none of them is clipped.
function _sector_loop(s::DiskSurface, (r0, r1, a, b); m = 16)
    φs = range(a, b; length = m + 1)
    r0 == 0 && return pushfirst!([_sector_point(s, r1, φ) for φ in φs], s.center)
    return [[_sector_point(s, r1, φ) for φ in φs]; [_sector_point(s, r0, φs[k]) for k in m+1:-1:1]]
end

# F(dA at p → receiver). A sector lying entirely in front of p takes the contour over
# its boundary loop (about a third of the edges of its triangles, which would all be
# unclipped and of one orientation, so |Σ| = Σ|·|); otherwise the triangles are summed
# one by one. The two agree to round-off.
function _receiver_vf(p, n, r::_Receiver)
    loop = r.loop
    if !isempty(loop)
        all(v -> dot(n, v - p) > 0, loop) && return _contour_vf(p, n, loop, length(loop))
        all(v -> dot(n, v - p) <= 0, loop) && return 0.0
    end
    f = 0.0
    for tri in r.tris
        f += _point_polygon_vf(p, n, tri)
    end
    return f
end

# Σₚ wₚ F(dAₚ → element j) with the exact point-to-polygon view factor, so only the
# (bounded, continuous) outer integral is done numerically.
function _pair_kernel((pi, wi, ni), r::_Receiver)
    s = 0.0
    for p in eachindex(pi)
        s += wi[p] * _receiver_vf(pi[p], ni[p], r)
    end
    return s
end

# Point-to-point kernel ∫∫ cosθᵢ cosθⱼ / (π r²) over two quadratures.
function _point_kernel((pi, wi, ni), (pj, wj, nj))
    s = 0.0
    for p in eachindex(pi), q in eachindex(pj)
        d = pj[q] - pi[p]
        r2 = dot(d, d)
        ci = dot(ni[p], d)
        cj = -dot(nj[q], d)
        (ci > 0 && cj > 0) || continue
        s += wi[p] * wj[q] * ci * cj / (π * r2^2)
    end
    return s
end

# Self-exchange of a curved concave element, ∫∫ over the patch with itself. Coincident
# points are skipped and the sum rescaled by A²/(A² − Σw²), which is exact for a
# constant kernel (the sphere) and accurate for the bounded kernel of a cylinder.
function _self_kernel((p, w, n))
    s = 0.0
    for a in eachindex(p), b in eachindex(p)
        a == b && continue
        d = p[b] - p[a]
        r2 = dot(d, d)
        ca = dot(n[a], d)
        cb = -dot(n[b], d)
        (ca > 0 && cb > 0) || continue
        s += w[a] * w[b] * ca * cb / (π * r2^2)
    end
    A = sum(w)
    return s * A^2 / (A^2 - sum(abs2, w))
end

"""
    _point_polygon_vf(p, n, tri)

Exact view factor from a differential area at `p` with normal `n` to the triangle
`tri` (clipped to the half-space in front of `p`), by the contour formula
`F = |Σₖ γₖ n·ĝₖ| / 2π`, where `γₖ` is the angle subtended by edge `k` and `ĝₖ` the unit
normal of the plane through `p` and that edge.
"""
function _point_polygon_vf(p::Vec3, n::Vec3, tri)
    h = map(v -> dot(n, v - p), tri)
    all(<=(0), h) && return 0.0
    all(>(0), h) && return _contour_vf(p, n, tri, 3)
    # clipped polygon as an immutable tuple (an MVector here escaped to the heap, one
    # allocation per call)
    z = zero(Vec3)
    poly = (z, z, z, z)
    m = 0
    for k in 1:3
        a, b = tri[k], tri[mod1(k + 1, 3)]
        ha, hb = h[k], h[mod1(k + 1, 3)]
        ha > 0 && (m += 1; poly = Base.setindex(poly, a, m))
        if (ha > 0) != (hb > 0) && ha != hb
            λ = ha / (ha - hb)
            m += 1; poly = Base.setindex(poly, a + λ * (b - a), m)
        end
    end
    m < 3 && return 0.0
    return _contour_vf(p, n, poly, m)
end

# Contour sum over the first `m` vertices of `poly`, each relative vertex formed once.
@inline function _contour_vf(p, n, poly, m)
    s = 0.0
    a1 = poly[1] - p
    a = a1
    for k in 1:m
        b = k == m ? a1 : poly[k+1] - p
        c = cross(a, b)
        nc = norm(c)
        if nc != 0
            s += atan(nc, dot(a, b)) * dot(n, c) / nc
        end
        a = b
    end
    return abs(s) / (2π)
end

"""
    closure_error(F) -> Vector

`1 - Σⱼ Fᵢⱼ` per element: the fraction of each element's emission that leaves the
geometry entirely (zero for every element of a closed enclosure, up to discretization
error).
"""
closure_error(F::AbstractMatrix) = 1 .- vec(sum(F; dims = 2))

# Adaptive outer quadrature for near pairs: emitter cells are subdivided while they lie
# within `ratio` cell sizes of the receiver, so the points concentrate along a shared
# edge, where the point-to-polygon factor has a logarithmic slope (it tends to ½ at the
# edge of two perpendicular faces). Uniform refinement converges slowly there from the
# side of the element that is long across the edge.
function _adaptive_kernel(t::Target, e, rc::_Receiver; mindepth::Integer = 2, maxdepth::Integer = 8, ratio = 2.0)
    s = t.components[t.elem_comp[e]].surface
    sgn = t.elem_curved[e]
    if s isa DiskSurface && sgn != 0
        r0, r1, a, b = element_sector(t, e)
        return _adapt_sector(s, sgn * s.normal, r0, r1, a, b, rc, 0, mindepth, maxdepth, ratio)
    end
    total = 0.0
    for k in t.elem_tris[e]
        v0 = t.tri_v0[k]
        total += _adapt_tri(s, sgn, t.elem_normal[e], v0, v0 + t.tri_e1[k], v0 + t.tri_e2[k],
                            rc, 0, mindepth, maxdepth, ratio)
    end
    return total
end

function _adapt_tri(s, sgn, nflat, a, b, c, rc, depth, mindepth, maxdepth, ratio)
    cen = (a + b + c) / 3
    p, n = (s === nothing || sgn == 0) ? (cen, nflat) : (q = surface_point(s, cen); (q[1], sgn * q[2]))
    h = max(norm(b - a), norm(c - b), norm(a - c))
    if depth < mindepth || (depth < maxdepth && _distance_to(p, rc) < ratio * h)
        ab, bc, ca = (a + b) / 2, (b + c) / 2, (c + a) / 2
        return _adapt_tri(s, sgn, nflat, a, ab, ca, rc, depth + 1, mindepth, maxdepth, ratio) +
               _adapt_tri(s, sgn, nflat, ab, b, bc, rc, depth + 1, mindepth, maxdepth, ratio) +
               _adapt_tri(s, sgn, nflat, ca, bc, c, rc, depth + 1, mindepth, maxdepth, ratio) +
               _adapt_tri(s, sgn, nflat, ab, bc, ca, rc, depth + 1, mindepth, maxdepth, ratio)
    end
    w = (s === nothing || sgn == 0) ? norm(cross(b - a, c - a)) / 2 : _curved_area(s, a, b, c)
    return w * _receiver_vf(p, n, rc)
end

function _adapt_sector(s, n, r0, r1, a, b, rc, depth, mindepth, maxdepth, ratio)
    rm = sqrt((r0^2 + r1^2) / 2)
    p = _sector_point(s, rm, (a + b) / 2)
    h = max(r1 - r0, r1 * (b - a))
    if depth < mindepth || (depth < maxdepth && _distance_to(p, rc) < ratio * h)
        φm = (a + b) / 2
        return _adapt_sector(s, n, r0, rm, a, φm, rc, depth + 1, mindepth, maxdepth, ratio) +
               _adapt_sector(s, n, rm, r1, a, φm, rc, depth + 1, mindepth, maxdepth, ratio) +
               _adapt_sector(s, n, r0, rm, φm, b, rc, depth + 1, mindepth, maxdepth, ratio) +
               _adapt_sector(s, n, rm, r1, φm, b, rc, depth + 1, mindepth, maxdepth, ratio)
    end
    return (r1^2 - r0^2) / 2 * (b - a) * _receiver_vf(p, n, rc)
end

_distance_to(p, r::_Receiver) = minimum(_point_triangle_distance(p, tri...) for tri in r.tris)

# Distance from p to the triangle (a, b, c) (Ericson, Real-Time Collision Detection 5.1.5).
function _point_triangle_distance(p, a, b, c)
    ab, ac, ap = b - a, c - a, p - a
    d1, d2 = dot(ab, ap), dot(ac, ap)
    (d1 <= 0 && d2 <= 0) && return norm(ap)
    bp = p - b
    d3, d4 = dot(ab, bp), dot(ac, bp)
    (d3 >= 0 && d4 <= d3) && return norm(bp)
    vc = d1 * d4 - d3 * d2
    (vc <= 0 && d1 >= 0 && d3 <= 0) && return norm(p - (a + d1 / (d1 - d3) * ab))
    cp = p - c
    d5, d6 = dot(ab, cp), dot(ac, cp)
    (d6 >= 0 && d5 <= d6) && return norm(cp)
    vb = d5 * d2 - d1 * d6
    (vb <= 0 && d2 >= 0 && d6 <= 0) && return norm(p - (a + d2 / (d2 - d6) * ac))
    va = d3 * d6 - d5 * d4
    if va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0
        return norm(p - (b + (d4 - d3) / ((d4 - d3) + (d5 - d6)) * (c - b)))
    end
    den = 1 / (va + vb + vc)
    return norm(p - (a + ab * (vb * den) + ac * (vc * den)))
end
