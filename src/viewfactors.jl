# Element-to-element view factors with occlusion.

"""
    view_factors(target; near = 3.0, level = 1) -> Matrix{Float64}

Dense view-factor matrix `F[i, j]`: the fraction of diffuse radiation leaving element
`i` that arrives at element `j`. Distant pairs use the point kernel
`cosθᵢ cosθⱼ Aⱼ / (π r²)` between centroids. Pairs closer than `near·√max(Aᵢ, Aⱼ)`
integrate the exact differential-area-to-polygon view factor onto `j` over `4^level`
sub-triangle points of `i`, which stays accurate for elements sharing an edge (e.g. the
wall–cap corner of a hohlraum). Visibility is a single
centroid-to-centroid ray against all opaque surfaces. Reciprocity
`Aᵢ Fᵢⱼ = Aⱼ Fⱼᵢ` holds exactly by construction.

Rows sum to one for a closed enclosure (counting apertures); see `closure_error`.
"""
function view_factors(t::Target; near = 3.0, level::Integer = 1, exact_far::Bool = true)
    N = nelements(t)
    F = zeros(N, N)
    C, Nrm, A = t.elem_centroid, t.elem_normal, t.elem_area
    quad = [element_quadrature(t, e, level) for e in 1:N]
    ε = _ray_eps(t)
    near2 = near^2
    Threads.@threads :dynamic for i in 1:N
        oi = C[i] + ε * Nrm[i]
        for j in i+1:N
            d = C[j] - C[i]
            r2 = dot(d, d)
            r2 == 0 && continue
            # elements must face each other (centroid test)
            dot(Nrm[i], d) > 0 && dot(Nrm[j], d) < 0 || continue
            oj = C[j] + ε * Nrm[j]
            any_hit(t.bvh, oi, oj - oi, 1 - 1e-9) && continue
            # G = Aᵢ Fᵢⱼ = Aⱼ Fⱼᵢ, averaged over evaluations from both ends
            if r2 < near2 * max(A[i], A[j])
                G = (_pair_kernel(t, quad[i], Nrm[i], j) + _pair_kernel(t, quad[j], Nrm[j], i)) / 2
            elseif exact_far
                G = (_pair_kernel(t, ((C[i],), (A[i],)), Nrm[i], j) +
                     _pair_kernel(t, ((C[j],), (A[j],)), Nrm[j], i)) / 2
            else
                ci = dot(Nrm[i], d)
                cj = -dot(Nrm[j], d)
                G = ci * cj * A[i] * A[j] / (π * r2^2)
            end
            F[i, j] = G / A[i]
            F[j, i] = G / A[j]
        end
    end
    return F
end

# Near pairs: Σₚ wₚ F(dAₚ → element j) with the exact point-to-polygon view factor, so
# only the (bounded, continuous) outer integral is done numerically.
function _pair_kernel(t::Target, (pi, wi), ni, j)
    s = 0.0
    for p in eachindex(pi)
        f = 0.0
        for k in t.elem_tris[j]
            v0 = t.tri_v0[k]
            f += _point_polygon_vf(pi[p], ni, (v0, v0 + t.tri_e1[k], v0 + t.tri_e2[k]))
        end
        s += wi[p] * f
    end
    return s
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
