# Element-to-element view factors with occlusion.

"""
    view_factors(target; near = 3.0, level = 1, visibility = target.visibility, visibility_grid = 4,
                 storage = Float64, tol, eta, leaf, maxrank) -> Matrix{Float64}, Float32Matrix or HMatrix

Dense view-factor matrix `F[i, j]`: the fraction of diffuse radiation leaving element
`i` that arrives at element `j`. Reciprocity `Aᵢ Fᵢⱼ = Aⱼ Fⱼᵢ` holds exactly by
construction, and rows sum to one for a closed enclosure (counting apertures); see
`closure_error`. `storage = Float32` assembles into a [`Float32Matrix`](@ref) instead
(half the memory; each factor rounded to single precision, products in Float64).
`storage = :hmatrix` assembles a hierarchical [`HMatrix`](@ref): near cluster pairs
dense, far ones of uniform visibility low rank (keywords `tol`, `eta`, `leaf`,
`maxrank`; `Briareus.HMATRIX_DEFAULTS`); on a hohlraum its assembly and storage grow as
N^1.3–1.6 instead of N² (2–3× faster, 5–8× smaller at 17–29k elements).

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

Elements of more than two facets (remeshed CAD, `remesh_mesh`) are not integrated facet
by facet: a flat one receives through its outline (exactly the sum over its facets, to
round-off; piece by piece, each piece its facets in one plane, where the point's tangent
plane cuts it), and near pairs integrate over cells that lump the emitter's facets (side
size / (2^(level+2) √2) before refining toward the receiver, 2^(level+1) on a curved
surface), so a pair costs a few times a two-facet pair whatever the facet count — as
long as the receiver is: each point's factor sums its outline, or its pieces when the
point's tangent plane cuts it, so an element on a finely faceted curved wall (hundreds
of sliver facets, tens of pieces and a jagged outline of hundreds of edges) costs tens
of times more per near pair. Give such parts their fitted surface
(`curved_components`): they are then integrated on the curved surface instead.

Visibility is a single centroid-to-centroid ray between the flat facets; the facets of
a curved surface (cylinder, sphere or cone) do not block a ray between two points inside it,
since the true surface cannot. On curved elements the centroid facing test is relaxed
by the spread of the true normals across the flat facet, so pairs the curved patches
exchange across are kept (a general triangle on a cylinder has its plane tilted axially
by about its size over 2R), and the point kernels apply the exact facing per point.

A pair that sees each other only partly, past the rim of an entrance-hole bore or a gap
or around the capsule, is all or nothing by its centroid ray. `visibility = :sampled`
(the default for targets with rims: `cad_target`, `visrad_target`) estimates the visible fraction instead: up to four probe rays per element (from points
near its extremes to the other's centroid) test whether the pair straddles an occluder;
if they all agree with the centroid ray, the factor is the `:centroid` one bit for bit,
and otherwise rays between the `visibility_grid²` sample points of each element (its
area lumped into a g × g grid of its tangent plane) give the kernel-weighted visible
fraction, which scales the unobstructed factor (reciprocity is kept). On an element
folded across more than two facets (remeshed CAD) the samples are put on its facets, and
its own facets do not occlude them. On a hohlraum with
0.1 mm LEH bores and a 20 µm waist gap at 0.5 mm, it brings the largest row-sum error
before balancing from 1.0% to 0.19% on the walls and from 5.9% to 0.7% on the LEH
apertures, for 2–4× the time; features narrower than a sample cell (the gap's slot) are
not resolved (research notes, "Partial visibility"). Against a Monte Carlo radiosity on
the same facets, it cuts the balanced capsule a₂/a₀ error there (0.51 → 0.31 pp at
0.5 mm, 0.30 → 0.14 pp at 0.25 mm; research notes, "Monte Carlo radiosity reference").
"""
function view_factors(t::Target; near = 3.0, mid = 12.0, level::Integer = 1, exact_far::Bool = true,
                      visibility::Symbol = t.visibility, visibility_grid::Integer = 4,
                      storage = Float64, hmatrix...)
    visibility in (:centroid, :sampled) ||
        throw(ArgumentError("visibility must be :centroid or :sampled, got :$visibility"))
    visibility_grid >= 1 || throw(ArgumentError("visibility_grid must be at least 1"))
    storage in (Float64, Float32, :hmatrix) ||
        throw(ArgumentError("storage must be Float64, Float32 or :hmatrix, got $storage"))
    storage === :hmatrix || isempty(hmatrix) ||
        throw(ArgumentError("keywords $(join(keys(hmatrix), ", ")) apply to storage = :hmatrix only"))
    ctx = _vf_context(t; near, mid, level, exact_far, visibility, visibility_grid)
    storage === :hmatrix && return _hmatrix_view_factors(ctx; hmatrix...)
    N = nelements(t)
    F = storage === Float32 ? Float32Matrix(zeros(Float32, N, N)) : zeros(N, N)
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

# Everything the pair kernel `_vf_pair` reads.
function _vf_context(t::Target; near, mid, level, exact_far, visibility, visibility_grid)
    N = nelements(t)
    C, Nrm, A = t.elem_centroid, t.elem_normal, t.elem_area
    curved = t.elem_curved .!= 0
    # planar elements (flat facets and polar sectors) receive through the exact
    # point-to-polygon kernel; sectors are polygons with finely divided arcs
    sector = [curved[e] && t.components[t.elem_comp[e]].surface isa DiskSurface for e in 1:N]
    planar = .!curved .| sector
    # elements of more than two facets (imported CAD) receive through their outline and
    # emit from cells that lump their facets (see `_PolyReceiver`, `_CellEmitter`)
    multi = [length(t.elem_tris[e]) > 2 && !sector[e] for e in 1:N]
    offset = cumsum([0; [length(c.mesh.triangles) for c in t.components]])
    recv = Union{_Receiver,_PolyReceiver}[_receiver(t, e, sector[e], multi[e] && planar[e], offset) for e in 1:N]
    # (curved elements at the resolution of their point rule `qnear`, flat ones finer: see
    # `_CellEmitter`)
    cells = Any[multi[e] ? _cell_emitter(t, e, curved[e] ? level + 1 : level + 2) : nothing for e in 1:N]
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
    # a flat facet standing in for a curved patch has a normal that departs from the true
    # surface normals across the patch: by Δφ/2 at the edges of a chord rectangle, and
    # also axially for a general triangle on a cylinder (its plane is tilted by ~h/2R).
    # The centroid facing test is relaxed by that spread (`slack`, a sine), so a pair the
    # curved patches exchange across is not culled; the point kernels then apply the
    # exact per-point facing. A visibility ray admitted only by the slack (heading behind
    # its own facet's plane) starts that much farther off the facet, so it clears the
    # facet; every other ray starts ε off it, as for flat elements (lifting all origins
    # would bias them toward the middle of a cavity, and into the capsule's shadow).
    slack = zeros(N)
    for e in 1:N
        curved[e] && !sector[e] || continue
        s = t.components[t.elem_comp[e]].surface
        for k in t.elem_tris[e], v in (t.tri_v0[k], t.tri_v0[k] + t.tri_e1[k], t.tri_v0[k] + t.tri_e2[k])
            slack[e] = max(slack[e], norm(cross(t.elem_curved[e] * last(surface_point(s, v)), Nrm[e])))
        end
    end
    # A segment between two points inside a convex curved surface (cylinder, sphere, cone)
    # cannot cross the surface, so its facets do not occlude such a pair: they would only
    # by faceting (a ray grazing the wall from a small facet is blocked by a larger
    # neighbor, which sits deeper inside the true surface).
    surfs = unique(c.surface for c in t.components if c.surface isa Union{CylinderSurface,SphereSurface,ConeSurface})
    comp_surf = [something(findfirst(==(c.surface), surfs), 0) for c in t.components]
    tri_surf = Int32[comp_surf[c] for c in t.tri_comp]
    inside = BitMatrix(undef, length(surfs), N)
    for e in 1:N, (k, s) in enumerate(surfs)
        inside[k, e] = _inside(s, C[e])
    end
    # partial visibility (`visibility = :sampled`): each element's ray samples and probes
    vis = visibility === :sampled ? [_visibility_samples(t, e, visibility_grid) for e in 1:N] : nothing
    return (; t, bvh = t.bvh, C, Nrm, A, planar, concave, recv, cells, qnear, qmid, qfar, qself,
           size, slack, ε, tri_surf, inside, near2 = near^2, mid2 = mid^2, level, exact_far, vis, multi)
end

# Row i of the view-factor matrix (pairs i < j), behind a function barrier so the
# threaded loop is specialized on the concrete types in `ctx`.
function _vf_row!(F, ctx, i)
    A = ctx.A
    N = length(A)
    if ctx.concave[i]
        F[i, i] = _self_kernel(ctx.qself[i]) / A[i]
    end
    for j in i+1:N
        G = _vf_pair(ctx, i, j)
        G == 0 && continue
        F[i, j] = G / A[i]
        F[j, i] = G / A[j]
    end
    return nothing
end

# Exchange factor G = Aᵢ Fᵢⱼ = Aⱼ Fⱼᵢ of the pair i < j (zero if they do not face or
# see each other); `_vf_exchange` takes any order and the diagonal.
_vf_exchange(ctx, i, j) = i == j ? (ctx.concave[i] ? _self_kernel(ctx.qself[i]) : 0.0) :
                          i < j ? _vf_pair(ctx, i, j) : _vf_pair(ctx, j, i)

# Visibility of the pair i < j: -1 if they do not face each other (or coincide), else
# the visible fraction (0 or 1 by the centroid ray, between with `:sampled`); and the
# centroid separation `d` with its square.
function _vf_visible(ctx, i, j)
    (; t, bvh, C, Nrm, size, slack, ε, tri_surf, inside, vis, multi) = ctx
    d = C[j] - C[i]
    r2 = dot(d, d)
    r2 == 0 && return -1.0, d, r2
    # elements must face each other (centroid test, relaxed on curved surfaces)
    r = sqrt(r2)
    ci, cj = dot(Nrm[i], d) / r, -dot(Nrm[j], d) / r
    ci > -slack[i] && cj > -slack[j] || return -1.0, d, r2
    oi = C[i] + (ε + size[i] * max(0.0, -ci)) * Nrm[i]
    oj = C[j] + (ε + size[j] * max(0.0, -cj)) * Nrm[j]
    skip = k -> (s = tri_surf[k]; s != 0 && inside[s, i] && inside[s, j])
    blocked = any_hit(bvh, oi, oj - oi, 1 - 1e-9, skip)
    vis === nothing && return blocked ? 0.0 : 1.0, d, r2
    # a folded element's own facets do not occlude its samples: the element model
    # radiates from its mean plane, and on a convex part each sample's rays toward
    # the horizon of its own facet would be blocked by that facet or its fold
    # neighbors (all of the reduced capsule rows on a remeshed CAD sphere)
    own(k, e) = multi[e] && (t.tri_front[k] == e || t.tri_back[k] == e)
    skipv = multi[i] || multi[j] ? (k -> skip(k) || own(k, i) || own(k, j)) : skip
    frac = _visible_fraction(bvh, vis[i], vis[j], oi, oj, oi - C[i], oj - C[j], Nrm[i], Nrm[j], blocked, skipv)
    return frac, d, r2
end

function _vf_pair(ctx, i, j)
    frac, d, r2 = _vf_visible(ctx, i, j)
    frac <= 0 && return 0.0
    return _vf_kernel(ctx, i, j, frac, d, r2)
end

# G of the pair i < j that faces and sees itself by the fraction `frac` > 0.
function _vf_kernel(ctx, i, j, frac, d, r2)
    (; t, Nrm, A, planar, recv, cells, qnear, qmid, qfar, size, near2, mid2, level, exact_far) = ctx
    # (`cells` is untyped: assert the kernel's type, or every pair's factor is boxed)
    near_kernel(e, rc) = (cells[e] === nothing ? _adaptive_kernel(t, e, rc; mindepth = level + 1) :
                          _cell_kernel(cells[e], rc))::Float64
    s2 = max(size[i], size[j])^2
    tier = r2 < near2 * s2 ? :near : r2 < mid2 * s2 ? :mid : :far
    # G = Aᵢ Fᵢⱼ = Aⱼ Fⱼᵢ
    if planar[i] && planar[j]
        if tier === :near
            G = (near_kernel(i, recv[j]) + near_kernel(j, recv[i])) / 2
        elseif exact_far
            G = (_pair_kernel(qfar[i], recv[j]) + _pair_kernel(qfar[j], recv[i])) / 2
        else
            G = dot(Nrm[i], d) * -dot(Nrm[j], d) * A[i] * A[j] / (π * r2^2)
        end
    else
        q = tier === :near ? qnear : tier === :mid ? qmid : qfar
        G = !planar[i] && !planar[j] ? _point_kernel(q[i], q[j]) :
            tier === :near ? (planar[j] ? near_kernel(i, recv[j]) : near_kernel(j, recv[i])) :
            planar[j] ? _pair_kernel(q[i], recv[j]) : _pair_kernel(q[j], recv[i])
    end
    frac == 1 || (G *= frac)
    return G::Float64
end

# ------------------------------------------------------------ partial visibility
#
# `visibility = :sampled`. The centroid ray decides a pair that is wholly visible or
# wholly hidden, but a pair that sees past a rim (an LEH bore's, a gap's) is partly
# both. Each element carries ray samples, its facets lumped into a g × g grid of its
# tangent plane (area-weighted mean point and area per cell), and up to four probes, the
# samples farthest along ±e₁ and ±e₂. Probe rays go from each probe of one element to
# the other's centroid. A pair whose probes all agree with its centroid ray keeps the
# centroid's verdict, so the factor is the `:centroid` one bit for bit; otherwise rays
# between every pair of samples give the visible fraction, weighted by the point kernel
# cosθᵢ cosθⱼ / r² (by area / r² if no sample pair faces), and the factor of the
# unobstructed pair is scaled by it. The fraction is computed once per unordered pair
# from segments that do not depend on the direction, so reciprocity holds exactly.

struct _VisSamples
    points::Vector{Vec3}
    weights::Vector{Float64}
    probes::Vector{Int}               # indices into `points`
end

function _visibility_samples(t::Target, e, g)
    C, n = t.elem_centroid[e], t.elem_normal[e]
    e1 = normalize(cross(n, abs(n[1]) < 0.9 ? Vec3(1, 0, 0) : Vec3(0, 1, 0)))
    e2 = cross(n, e1)
    pts, ws = Vec3[], Float64[]
    # sub-triangles of under about half a grid cell
    L = ceil(Int, log2(g)) + 1
    for k in t.elem_tris[e]
        _subtri!(pts, ws, t.tri_v0[k], t.tri_v0[k] + t.tri_e1[k], t.tri_v0[k] + t.tri_e2[k], L)
    end
    u = [dot(p - C, e1) for p in pts]
    v = [dot(p - C, e2) for p in pts]
    u0, v0 = minimum(u), minimum(v)
    du, dv = max(maximum(u) - u0, eps()) / g, max(maximum(v) - v0, eps()) / g
    W = zeros(g, g)
    P = fill(zero(Vec3), g, g)
    for (k, p) in enumerate(pts)
        a = clamp(floor(Int, (u[k] - u0) / du) + 1, 1, g)
        b = clamp(floor(Int, (v[k] - v0) / dv) + 1, 1, g)
        W[a, b] += ws[k]
        P[a, b] += ws[k] * p
    end
    keep = findall(>(0), vec(W))
    points = [vec(P)[c] / vec(W)[c] for c in keep]
    weights = vec(W)[keep]
    # A folded element (remeshed CAD: facets within a few degrees of one plane) has its
    # cell means off its facets, behind them on a convex part (a capsule), where their
    # rays start inside the solid and are blocked by the element's own facets. Each is
    # moved along the normal onto the facet that line meets, as `Target` moves the
    # centroid (`_onto_patch`), or, if it meets none (a cell mean outside a non-convex
    # outline), to the nearest sub-triangle point of its cell. Planar elements, and those
    # of one or two facets (every primitive), keep their samples bit for bit.
    tris = t.elem_tris[e]
    if length(tris) > 2
        h = maximum(norm(p - C) for p in pts)
        if any(p -> abs(dot(p - C, n)) > 1e-9 * h, pts)
            cell = [clamp(floor(Int, (u[k] - u0) / du) + 1, 1, g) + g * (clamp(floor(Int, (v[k] - v0) / dv) + 1, 1, g) - 1)
                    for k in eachindex(pts)]
            for (i, c) in enumerate(keep)
                p, hit = _onto_facets(t, tris, points[i], n, h)
                if !hit
                    ks = findall(==(c), cell)
                    p = pts[ks[argmin([norm(pts[k] - points[i]) for k in ks])]]
                end
                points[i] = p
            end
        end
    end
    probes = unique!([argmax(p -> dot(points[p] - C, d), eachindex(points)) for d in (e1, -e1, e2, -e2)])
    return _VisSamples(points, weights, probes)
end

# `p` moved along `n` onto the nearest of the facets `tris` the line through it meets
# (the rule of `_onto_patch`), and whether it meets one; a move under 1e-9 h is dropped.
function _onto_facets(t::Target, tris, p, n, h)
    best = Inf
    for k in tris
        a, nt = t.tri_v0[k], cross(t.tri_e1[k], t.tri_e2[k])
        den = dot(n, nt)
        den == 0 && continue
        s = dot(a - p, nt) / den
        q = p + s * n
        b, d = a + t.tri_e1[k], a + t.tri_e2[k]
        tol = -1e-9 * dot(nt, nt)
        dot(_tri_cross(a, b, q), nt) >= tol && dot(_tri_cross(b, d, q), nt) >= tol &&
            dot(_tri_cross(d, a, q), nt) >= tol && abs(s) < abs(best) && (best = s)
    end
    isfinite(best) || return p, false
    return abs(best) > 1e-9 * h ? p + best * n : p, true
end

# Visible fraction of a pair whose centroid ray (from `oi` to `oj`) is `blocked` or not;
# `li`, `lj` lift the samples off their facets as the centroid ray's ends are lifted.
function _visible_fraction(bvh, si::_VisSamples, sj::_VisSamples, oi, oj, li, lj, ni, nj, blocked, skip)
    seen(p, q) = !any_hit(bvh, p, q - p, 1 - 1e-9, skip)
    agree = all(k -> seen(si.points[k] + li, oj) != blocked, si.probes) &&
            all(k -> seen(oi, sj.points[k] + lj) != blocked, sj.probes)
    agree && return blocked ? 0.0 : 1.0
    vis, tot, visa, tota = 0.0, 0.0, 0.0, 0.0
    for a in eachindex(si.points), b in eachindex(sj.points)
        p, q = si.points[a] + li, sj.points[b] + lj
        d = q - p
        r2 = dot(d, d)
        r2 > 0 || continue
        wa = si.weights[a] * sj.weights[b] / r2
        w = wa * max(dot(ni, d), 0.0) * max(-dot(nj, d), 0.0) / r2
        s = seen(p, q)
        tot += w; tota += wa
        s && (vis += w; visa += wa)
    end
    return tot > 0 ? vis / tot : tota > 0 ? visa / tota : 0.0
end

# Whether `p` lies inside the (convex) curved surface, to round-off.
_inside(s::CylinderSurface, p) =
    (d = p - s.origin; norm(d - dot(d, s.axis) * s.axis) <= s.radius * (1 + 1e-12))
_inside(s::SphereSurface, p) = norm(p - s.center) <= s.radius * (1 + 1e-12)
# (a cone's nappe bounds a convex solid: points within the half-angle of the axis)
function _inside(s::ConeSurface, p)
    sα, cα = sincosd(s.half_angle)
    z, ρ, _ = _cone_coordinates(s, p)
    return ρ * cα - z * sα <= 1e-12 * norm(p - s.apex)
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

function _receiver(t::Target, e, sector::Bool, poly::Bool, offset)
    if sector
        s, sec = t.components[t.elem_comp[e]].surface, element_sector(t, e)
        return _Receiver(_sector_triangles(s, sec), _sector_loop(s, sec))
    end
    poly && return _poly_receiver(t, e, offset)
    return _Receiver([(t.tri_v0[k], t.tri_v0[k] + t.tri_e1[k], t.tri_v0[k] + t.tri_e2[k])
                      for k in t.elem_tris[e]], Vec3[])
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

# ------------------------------------------------- elements of many facets (CAD)
#
# A remeshed CAD element carries 10–25 facets (slivers of the CAD faces, cut by the
# element's outline), and summing the point-to-polygon factor over each of them made a
# pair cost tens of times a primitive pair. The receiving side uses the element's outline
# instead, and the emitting side lumps its facets into cells; both reduce exactly (the
# receiver to round-off) or in point count to what a two-facet element costs.
#
# Receiver. The contour formula integrates a divergence-free field, so for a point that
# sees the whole element from in front (every vertex in front of its tangent plane, every
# facet facing the same way) the sum over the facets' contours equals the contour over
# the element's boundary edges: the interior edges cancel pairwise. That holds for a
# folded element as well as a planar one. A point whose tangent plane cuts the element
# takes the element piece by piece, a piece being its facets in one plane (to round-off):
# a planar piece cut by the plane is bounded by its edges in front plus one segment on
# the cutting line, whose multiplicities follow from the crossings (with two crossings, the
# segment joins them). That is again the sum over the clipped facets, to round-off; more
# than two crossings fall back to the facets. Straight runs of the outline (vertices
# that neighbors' facets need) are merged into single edges.

struct _Piece
    edges::Vector{NTuple{2,Vec3}}   # boundary edges of the piece's facets
    verts::Vector{Vec3}             # their end points
    normal::Vec3                    # the facets' normal
    point::Vec3                     # a vertex
    radius::Float64                 # largest distance of a vertex from `point`
    tris::Vector{NTuple{3,Vec3}}
end

struct _PolyReceiver
    edges::Vector{NTuple{2,Vec3}}   # boundary edges of the element
    verts::Vector{Vec3}             # every vertex of its facets
    pieces::Vector{_Piece}          # its facets grouped by plane
    tris::Vector{NTuple{3,Vec3}}    # every facet (distances)
    center::Vec3                    # bounding sphere
    radius::Float64
end

# Receiver of flat element `e` of more than two facets; `offset[c]` is the index of the
# component's first triangle in the target, less one.
function _poly_receiver(t::Target, e, offset)
    c = t.elem_comp[e]
    m = t.components[c].mesh
    ks = t.elem_tris[e]
    ids = Dict{Vec3,Int}()
    P = Vec3[]
    vid(v) = get!(() -> (push!(P, v); length(P)), ids, v)
    T = [map(i -> vid(m.vertices[i]), m.triangles[k - offset[c]]) for k in ks]
    center = t.elem_centroid[e]
    radius = maximum(v -> norm(v - center), P)
    # pieces: facets with every vertex within round-off of one plane, facing one way
    tol = 1e-12 * (radius + norm(center))
    groups, gnormal, gpoint = Vector{Int}[], Vec3[], Vec3[]
    for (f, k) in enumerate(ks)
        nk = t.tri_normal[k]
        g = findfirst(eachindex(groups)) do g
            dot(nk, gnormal[g]) > 0 && all(i -> abs(dot(gnormal[g], P[i] - gpoint[g])) <= tol, T[f])
        end
        if g === nothing
            push!(groups, [f]); push!(gnormal, nk); push!(gpoint, P[T[f][1]])
        else
            push!(groups[g], f)
        end
    end
    vec3(es) = [(P[a], P[b]) for (a, b) in es]
    pieces = map(eachindex(groups)) do g
        es = _merge_straight(P, _boundary_edges(T, groups[g]))
        vs = unique!([P[i] for e in es for i in e])
        _Piece(vec3(es), vs, gnormal[g], gpoint[g], maximum(v -> norm(v - gpoint[g]), vs),
               [(P[T[f][1]], P[T[f][2]], P[T[f][3]]) for f in groups[g]])
    end
    return _PolyReceiver(vec3(_merge_straight(P, _boundary_edges(T, eachindex(T)))), P, pieces,
                         [(P[a], P[b], P[d]) for (a, b, d) in T], center, radius)
end

# Directed boundary edges (vertex ids, sorted) of facets `fs` of `T`: edges not cancelled
# by an opposite edge of another facet.
function _boundary_edges(T, fs)
    count = Dict{NTuple{2,Int},Int}()
    for f in fs
        a, b, c = T[f]
        for (u, v) in ((a, b), (b, c), (c, a))
            if get(count, (v, u), 0) > 0
                count[(v, u)] -= 1
            else
                count[(u, v)] = get(count, (u, v), 0) + 1
            end
        end
    end
    return sort!([k for (k, n) in count for _ in 1:n])
end

# Boundary edges with straight runs merged: u→v→w becomes u→w where v is the only edge
# out of and into v and lies on the segment uw to round-off. The contour terms of
# collinear edges add, so the sum is unchanged.
function _merge_straight(P, edges)
    outs, ins = Dict{Int,Vector{Int}}(), Dict{Int,Vector{Int}}()
    for (a, b) in edges
        push!(get!(outs, a, Int[]), b)
        push!(get!(ins, b, Int[]), a)
    end
    merged = true
    while merged
        merged = false
        for v in sort!(collect(keys(outs)))
            haskey(outs, v) && haskey(ins, v) || continue
            o, i = outs[v], ins[v]
            length(o) == 1 && length(i) == 1 || continue
            u, w = i[1], o[1]
            u != w && u != v && w != v && _straight(P[u], P[v], P[w]) || continue
            replace!(outs[u], v => w; count = 1)
            replace!(ins[w], v => u; count = 1)
            delete!(outs, v); delete!(ins, v)
            merged = true
        end
    end
    return sort!([(a, b) for (a, bs) in outs for b in bs])
end

# Whether v lies on the segment ab, to round-off.
function _straight(a, v, b)
    d = b - a
    s = dot(v - a, d)
    return 0 < s < dot(d, d) && norm(cross(v - a, d)) <= 1e-12 * dot(d, d)
end

function _receiver_vf(p, n, r::_PolyReceiver)
    # the bounding sphere wholly in front or behind decides it without the vertices
    hc = dot(n, r.center - p)
    hc < -1.000001 * r.radius && return 0.0
    front, behind = hc > 1.000001 * r.radius, false
    if !front
        front, behind = true, true
        for v in r.verts
            h = dot(n, v - p)
            front &= h > 0
            behind &= h <= 0
        end
    end
    behind && return 0.0
    if front
        length(r.pieces) == 1 && return _edges_vf(p, n, r.edges)
        pos, neg = false, false
        for pc in r.pieces
            σ = dot(pc.normal, p - pc.point)
            pos |= σ > 0
            neg |= σ < 0
        end
        pos && neg || return _edges_vf(p, n, r.edges)
    end
    f = 0.0
    for pc in r.pieces
        f += _piece_vf(p, n, pc)
    end
    return f
end

function _piece_vf(p, n, pc::_Piece)
    front, behind = true, true
    hmax = 0.0
    for v in pc.verts
        h = dot(n, v - p)
        front &= h > 0
        behind &= h <= 0
        hmax = max(hmax, abs(h))
    end
    behind && return 0.0
    front && return _edges_vf(p, n, pc.edges)
    # p's tangent plane parallel to the piece (a coplanar neighbor), or p in the piece's
    # plane (on a neighbor cut from the same CAD facet): the piece is seen edge-on, its
    # factor is zero, and the cutting line is ill-defined or passes through p
    (hmax <= 1e-6 * pc.radius || abs(dot(pc.normal, p - pc.point)) <= 1e-10 * pc.radius) && return 0.0
    return _clipped_vf(p, n, pc)
end

# The piece clipped to the half-space in front of p: its edges in front, cut where they
# cross, and the segment joining the two crossings (exit to entry). Other crossing
# counts (a piece whose outline crosses the line more than twice) sum the facets.
function _clipped_vf(p, n, pc::_Piece)
    s = 0.0
    nx, sx = 0, 0
    X1, X2 = zero(Vec3), zero(Vec3)
    for (a, b) in pc.edges
        ha, hb = dot(n, a - p), dot(n, b - p)
        if ha > 0 && hb > 0
            s += _edge_term(n, a - p, b - p)
        elseif (ha > 0) != (hb > 0)
            X = a + ha / (ha - hb) * (b - a)
            if ha > 0
                s += _edge_term(n, a - p, X - p)          # leaves the half-space at X
                nx == 0 && (sx = 1)
            else
                s += _edge_term(n, X - p, b - p)          # enters it at X
                nx == 0 && (sx = -1)
            end
            nx += 1
            nx == 1 ? (X1 = X) : nx == 2 ? (X2 = X) : nothing
        end
    end
    if nx != 2
        f = 0.0
        for tri in pc.tris
            f += _point_polygon_vf(p, n, tri)
        end
        return f
    end
    # the closing segment runs from the exit to the entry
    s += sx * _edge_term(n, X1 - p, X2 - p)
    return abs(s) / (2π)
end

# Σ over directed edges of the contour terms, |·|/2π.
function _edges_vf(p, n, edges)
    s = 0.0
    for (a, b) in edges
        s += _edge_term(n, a - p, b - p)
    end
    return abs(s) / (2π)
end

# Contour term of the edge a→b (relative to the point): its subtended angle times the
# normal component of the unit normal of the plane through the point and the edge.
# @noinline for the reason given at `_contour_vf`.
@noinline function _edge_term(n, a, b)
    c = cross(a, b)
    nc = norm(c)
    nc == 0 && return 0.0
    return atan(nc, dot(a, b)) * dot(n, c) / nc
end

_distance_to(p, r::_PolyReceiver) = minimum(_point_triangle_distance(p, tri...) for tri in r.tris)

# Whether p lies within d of the receiver (the bounding sphere first).
_closer_than(p, r::_Receiver, d) = _distance_to(p, r) < d
_closer_than(p, r::_PolyReceiver, d) = norm(p - r.center) - r.radius < d && _distance_to(p, r) < d

# Emitter. The adaptive near-pair rule refines toward the receiver from cells of the
# element's tangent plane rather than from its facets: the facets are bisected along
# their longest edges into pieces under half a cell, the pieces lumped by cell (one
# point per cell: the area-weighted mean of the pieces' centroids, with their total area;
# on a curved surface the pieces' exact curved areas, the mean projected onto it), and a
# cell within `ratio` cell sizes of the receiver (a cell's size taken as 2× its side, the
# longest edge of a right triangle of its area) is split into four, as `_adapt_tri` does
# with a triangle, down to cells of side size / (2^maxdepth √2); when splitting, the
# pieces are bisected only to under a child's side (under half, as at the root, cost 70%
# more and changed nothing at 1e-6 on the reference cases). The root cells have side
# size / (2^depth √2): 2·4^depth of them on a square element, whatever the facets'
# number and shape. `view_factors` roots flat elements one level below the two-facet rule
# (depth = level + 2): at the two-facet resolution the one-point cells carry the same
# +2e-4 bias as primitive elements on perpendicular plates (all pairs near); one level
# finer brings it to the 1e-5 of the per-facet rule on CAD elements, for 10-20% more
# time (research notes, "Flat CAD elements"). Curved elements stay at level + 1, the
# resolution of their lumped point rule (`element_quadrature`) for curved partners.
const _Cell = Tuple{Float64,Vec3,Vec3,Vector{NTuple{3,Vec3}},Int,Int}   # weight, point, normal, pieces, key

struct _CellEmitter{S}
    s::S                            # the curved surface, or nothing (flat)
    sgn::Int8
    C::Vec3
    n::Vec3                         # the element normal (flat)
    e1::Vec3
    e2::Vec3
    depth::Int                      # depth of the root cells
    hq::Float64                     # their side
    cells::Vector{_Cell}
end

function _cell_emitter(t::Target, e, depth)
    s = t.components[t.elem_comp[e]].surface
    sgn = t.elem_curved[e]
    S = s === nothing || sgn == 0 ? nothing : s
    C, n = t.elem_centroid[e], t.elem_normal[e]
    r = 0.0
    tris = [(t.tri_v0[k], t.tri_v0[k] + t.tri_e1[k], t.tri_v0[k] + t.tri_e2[k]) for k in t.elem_tris[e]]
    for tri in tris, v in tri
        r = max(r, norm(v - C))
    end
    e1 = normalize(cross(n, abs(n[1]) < 0.9 ? Vec3(1, 0, 0) : Vec3(0, 1, 0)))
    e2 = cross(n, e1)
    hq = r / 2^depth                                     # √2 r / (2^depth √2)
    em = _CellEmitter(S, sgn, C, n, e1, e2, Int(depth), hq, _Cell[])
    acc = Dict{NTuple{2,Int},Tuple{Float64,Vec3,Vector{NTuple{3,Vec3}}}}()
    for (a, b, c) in tris
        _bin_piece!(em, hq, hq / 2, a, b, c) do key, w, q, piece
            w0, q0, ps = get!(() -> (0.0, zero(Vec3), NTuple{3,Vec3}[]), acc, key)
            push!(ps, piece)
            acc[key] = (w0 + w, q0 + w * q, ps)
        end
    end
    for key in sort!(collect(keys(acc)))
        w, wq, ps = acc[key]
        w > 0 && push!(em.cells, _cell(em, w, wq, ps, key...))
    end
    return em
end

# A cell from its pieces' total weight and weighted point sum.
function _cell(em::_CellEmitter, w, wq, ps, i, j)
    em.s === nothing && return (w, wq / w, em.n, ps, i, j)
    q, nq = surface_point(em.s, wq / w)
    return (w, q, em.sgn * nq, ps, i, j)
end

function _cell_kernel(em::_CellEmitter, rc; maxdepth::Integer = 8, ratio = 2.0)
    total = 0.0
    for cell in em.cells
        total += _adapt_cell(em, cell, em.hq, rc, em.depth, maxdepth, ratio)
    end
    return total
end

# h = 2 hq: the longest edge of a right isosceles triangle of the cell's area
function _adapt_cell(em, (w, p, n, ps, i, j), hq, rc, depth, maxdepth, ratio)
    depth < maxdepth && _closer_than(p, rc, 2ratio * hq) &&
        return _refine_cell(em, ps, i, j, hq / 2, rc, depth + 1, maxdepth, ratio)
    return w * _receiver_vf(p, n, rc)
end

# Σ over the four children (side hq, depth `depth`) of cell (i, j), refined further
# where they are close to the receiver. The cell's pieces are bisected to under a child's
# side and go to the child their centroid lies in (or the nearest one, for a part beyond
# the cell). A first pass gives each child's weight and point; only the children that
# are refined again collect their pieces, in a second.
function _refine_cell(em::_CellEmitter, ps, i, j, hq, rc, depth, maxdepth, ratio)
    child((ki, kj)) = 1 + clamp(ki - 2i, 0, 1) + 2clamp(kj - 2j, 0, 1)
    ws = zeros(4)
    wqs = zeros(Vec3, 4)
    for (a, b, c) in ps
        _bin_piece!(em, hq, hq, a, b, c) do key, w, q, piece
            k = child(key)
            ws[k] += w
            wqs[k] += w * q
        end
    end
    total = 0.0
    refine = 0
    for k in 1:4
        ws[k] > 0 || continue
        w, p, n, _ = _cell(em, ws[k], wqs[k], ps, 0, 0)
        if depth < maxdepth && _closer_than(p, rc, 2ratio * hq)
            refine |= 1 << (k - 1)
        else
            total += w * _receiver_vf(p, n, rc)
        end
    end
    refine == 0 && return total
    pcs = _collect_pieces(em, ps, i, j, hq, refine)
    for k in 1:4
        refine & (1 << (k - 1)) != 0 || continue
        total += _refine_cell(em, pcs[k], 2i + (k - 1) % 2, 2j + (k - 1) ÷ 2, hq / 2, rc, depth + 1, maxdepth, ratio)
    end
    return total
end

# The pieces of the children of cell (i, j) (side hq) flagged in `mask`.
function _collect_pieces(em, ps, i, j, hq, mask)
    pcs = ntuple(k -> NTuple{3,Vec3}[], 4)
    for (a, b, c) in ps
        _bin_piece!(em, hq, hq, a, b, c, false) do (ki, kj), _, _, piece
            k = 1 + clamp(ki - 2i, 0, 1) + 2clamp(kj - 2j, 0, 1)
            mask & (1 << (k - 1)) != 0 && push!(pcs[k], piece)
        end
    end
    return pcs
end

# Bisect the piece (a, b, c) along its longest edge until that is under lmax, and hand
# each part to `f(key, weight, point, part)`: its cell (side hq, relative to the
# element's centroid), and, with `weights`, its area (curved on a surface) and its
# centroid (projected).
function _bin_piece!(f, em::_CellEmitter, hq, lmax, a, b, c, weights::Bool = true)
    lab, lbc, lca = sum(abs2, b - a), sum(abs2, c - b), sum(abs2, a - c)
    if max(lab, lbc, lca) > lmax^2
        if lab >= lbc && lab >= lca
            m = (a + b) / 2
            _bin_piece!(f, em, hq, lmax, a, m, c, weights)
            _bin_piece!(f, em, hq, lmax, m, b, c, weights)
        elseif lbc >= lca
            m = (b + c) / 2
            _bin_piece!(f, em, hq, lmax, a, b, m, weights)
            _bin_piece!(f, em, hq, lmax, a, m, c, weights)
        else
            m = (c + a) / 2
            _bin_piece!(f, em, hq, lmax, a, b, m, weights)
            _bin_piece!(f, em, hq, lmax, m, b, c, weights)
        end
        return
    end
    g = (a + b + c) / 3
    key = (floor(Int, dot(g - em.C, em.e1) / hq), floor(Int, dot(g - em.C, em.e2) / hq))
    s = em.s
    if !weights
        f(key, 0.0, g, (a, b, c))
    elseif s === nothing
        f(key, norm(cross(b - a, c - a)) / 2, g, (a, b, c))
    else
        f(key, _curved_area(s, a, b, c), first(surface_point(s, g)), (a, b, c))
    end
    return
end

# Σₚ wₚ F(dAₚ → element j) with the exact point-to-polygon view factor, so only the
# (bounded, continuous) outer integral is done numerically.
function _pair_kernel((pi, wi, ni), r)
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
# @noinline: when inlined into `_point_polygon_vf`, LLVM's x86 register coalescer can drop
# the upper lane of a vectorized vertex around the `atan` call (wrong view factors seen
# under --code-coverage on Julia 1.13 and --track-allocation on 1.12-1.14). No measurable cost.
@noinline function _contour_vf(p, n, poly, m)
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
closure_error(F::AbstractMatrix) = 1 .- _rowsums(F)

# Adaptive outer quadrature for near pairs: emitter cells are subdivided while they lie
# within `ratio` cell sizes of the receiver, so the points concentrate along a shared
# edge, where the point-to-polygon factor has a logarithmic slope (it tends to ½ at the
# edge of two perpendicular faces). Uniform refinement converges slowly there from the
# side of the element that is long across the edge.
function _adaptive_kernel(t::Target, e, rc; mindepth::Integer = 2, maxdepth::Integer = 8, ratio = 2.0)
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
