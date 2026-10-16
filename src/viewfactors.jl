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
    recv = [sector[e] ? _sector_triangles(t.components[t.elem_comp[e]].surface, element_sector(t, e)) :
            [(t.tri_v0[k], t.tri_v0[k] + t.tri_e1[k], t.tri_v0[k] + t.tri_e2[k]) for k in t.elem_tris[e]]
            for e in 1:N]
    concave = curved .& .!sector
    # point-to-point quadrature on curved surfaces carries an O((h/r)²) bias that does not
    # cancel between pairs, so curved elements keep several points per triangle at every
    # distance: level + 1 near, level out to `mid` element sizes, one per triangle beyond
    qnear = [element_quadrature(t, e, curved[e] && !sector[e] ? level + 1 : level) for e in 1:N]
    qmid = [element_quadrature(t, e, curved[e] ? level : 0) for e in 1:N]
    qfar = [curved[e] ? element_quadrature(t, e, 0) : ([C[e]], [A[e]], [Nrm[e]]) for e in 1:N]
    qself = [concave[e] ? element_quadrature(t, e, level + 1) : nothing for e in 1:N]
    ε = _ray_eps(t)
    # element size: largest distance from the centroid to a vertex, so a long sliver
    # counts as large even when its area is small (√2 × half-diagonal = side for a square)
    size = [sqrt(2) * maximum(max(norm(t.tri_v0[k] - C[e]), norm(t.tri_v0[k] + t.tri_e1[k] - C[e]),
                                  norm(t.tri_v0[k] + t.tri_e2[k] - C[e])) for k in t.elem_tris[e])
            for e in 1:N]
    near2 = near^2
    Threads.@threads :dynamic for i in 1:N
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
            any_hit(t.bvh, oi, oj - oi, 1 - 1e-9) && continue
            s2 = max(size[i], size[j])^2
            tier = r2 < near2 * s2 ? :near : r2 < mid^2 * s2 ? :mid : :far
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
    end
    return F
end

# One-point rule of a curved element: the projected centroid, its normal and the area.
function _one_point(t::Target, e)
    p, n = surface_point(t.components[t.elem_comp[e]].surface, t.elem_centroid[e])
    return [p], [t.elem_area[e]], [t.elem_curved[e] * n]
end

# Σₚ wₚ F(dAₚ → element j) with the exact point-to-polygon view factor, so only the
# (bounded, continuous) outer integral is done numerically.
function _pair_kernel((pi, wi, ni), tris)
    s = 0.0
    for p in eachindex(pi)
        f = 0.0
        for tri in tris
            f += _point_polygon_vf(pi[p], ni[p], tri)
        end
        s += wi[p] * f
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
    poly = MVector{4,Vec3}(undef)
    m = 0
    for k in 1:3
        a, b = tri[k], tri[mod1(k + 1, 3)]
        ha, hb = h[k], h[mod1(k + 1, 3)]
        ha > 0 && (m += 1; poly[m] = a)
        if (ha > 0) != (hb > 0) && ha != hb
            λ = ha / (ha - hb)
            m += 1; poly[m] = a + λ * (b - a)
        end
    end
    m < 3 && return 0.0
    s = 0.0
    for k in 1:m
        a = poly[k] - p
        b = poly[mod1(k + 1, m)] - p
        c = cross(a, b)
        nc = norm(c)
        nc == 0 && continue
        s += atan(nc, dot(a, b)) * dot(n, c) / nc
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
function _adaptive_kernel(t::Target, e, tris; mindepth::Integer = 2, maxdepth::Integer = 8, ratio = 2.0)
    s = t.components[t.elem_comp[e]].surface
    sgn = t.elem_curved[e]
    if s isa DiskSurface && sgn != 0
        r0, r1, a, b = element_sector(t, e)
        return _adapt_sector(s, sgn * s.normal, r0, r1, a, b, tris, 0, mindepth, maxdepth, ratio)
    end
    total = 0.0
    for k in t.elem_tris[e]
        v0 = t.tri_v0[k]
        total += _adapt_tri(s, sgn, t.elem_normal[e], v0, v0 + t.tri_e1[k], v0 + t.tri_e2[k],
                            tris, 0, mindepth, maxdepth, ratio)
    end
    return total
end

function _adapt_tri(s, sgn, nflat, a, b, c, tris, depth, mindepth, maxdepth, ratio)
    cen = (a + b + c) / 3
    p, n = (s === nothing || sgn == 0) ? (cen, nflat) : (q = surface_point(s, cen); (q[1], sgn * q[2]))
    h = max(norm(b - a), norm(c - b), norm(a - c))
    if depth < mindepth || (depth < maxdepth && _distance_to(p, tris) < ratio * h)
        ab, bc, ca = (a + b) / 2, (b + c) / 2, (c + a) / 2
        return _adapt_tri(s, sgn, nflat, a, ab, ca, tris, depth + 1, mindepth, maxdepth, ratio) +
               _adapt_tri(s, sgn, nflat, ab, b, bc, tris, depth + 1, mindepth, maxdepth, ratio) +
               _adapt_tri(s, sgn, nflat, ca, bc, c, tris, depth + 1, mindepth, maxdepth, ratio) +
               _adapt_tri(s, sgn, nflat, ab, bc, ca, tris, depth + 1, mindepth, maxdepth, ratio)
    end
    w = (s === nothing || sgn == 0) ? norm(cross(b - a, c - a)) / 2 : _curved_area(s, a, b, c)
    return w * sum(_point_polygon_vf(p, n, tri) for tri in tris)
end

function _adapt_sector(s, n, r0, r1, a, b, tris, depth, mindepth, maxdepth, ratio)
    rm = sqrt((r0^2 + r1^2) / 2)
    p = _sector_point(s, rm, (a + b) / 2)
    h = max(r1 - r0, r1 * (b - a))
    if depth < mindepth || (depth < maxdepth && _distance_to(p, tris) < ratio * h)
        φm = (a + b) / 2
        return _adapt_sector(s, n, r0, rm, a, φm, tris, depth + 1, mindepth, maxdepth, ratio) +
               _adapt_sector(s, n, rm, r1, a, φm, tris, depth + 1, mindepth, maxdepth, ratio) +
               _adapt_sector(s, n, r0, rm, φm, b, tris, depth + 1, mindepth, maxdepth, ratio) +
               _adapt_sector(s, n, rm, r1, φm, b, tris, depth + 1, mindepth, maxdepth, ratio)
    end
    return (r1^2 - r0^2) / 2 * (b - a) * sum(_point_polygon_vf(p, n, tri) for tri in tris)
end

_distance_to(p, tris) = minimum(_point_triangle_distance(p, tri...) for tri in tris)

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
