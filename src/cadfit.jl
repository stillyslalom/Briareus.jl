# Analytic surfaces in imported (CAD) meshes: patches that lie on one plane, cylinder
# or sphere are found by region growing with least-squares fits, and the curved ones are
# attached to their components as `CurvedSurface`s, so the radiation model integrates
# on the true surface as it does for primitive targets.

"""
    SurfaceRegion

Patches of a mesh lying on one analytic surface, as found by [`fit_surfaces`](@ref):

- `kind`: `:plane`, `:cylinder` or `:sphere`
- `surface`: the fitted `CylinderSurface` or `SphereSurface` (`nothing` for a plane)
- `centroid`: area-weighted centroid of the patches (mm)
- `normal`: a plane's unit normal, on the side most of its area faces (zero for
  curved surfaces)
- `patches`: indices of the patches in the mesh
- `area`: total area of the patches (mm², flat facets)
- `rms`, `maxdev`: rms and largest distance of the patches' vertices from the surface
  (mm)
- `facing`: for curved surfaces, `-1` when every patch faces the axis or center (the
  inside of a hohlraum wall), `+1` when every patch faces away from it, `0` when mixed;
  for planes `+1`, or `0` when some patches face the other way
"""
struct SurfaceRegion
    kind::Symbol
    surface::Union{Nothing,CylinderSurface,SphereSurface}
    centroid::Vec3
    normal::Vec3
    patches::Vector{Int}
    area::Float64
    rms::Float64
    maxdev::Float64
    facing::Int
end

function Base.show(io::IO, r::SurfaceRegion)
    p3(x) = @sprintf("(%.4g, %.4g, %.4g)", x[1], x[2], x[3])
    s = r.surface
    geo = r.kind === :plane ? "normal $(p3(r.normal))" :
          r.kind === :cylinder ? @sprintf("radius %.6g mm, axis %s through %s", s.radius, p3(s.axis), p3(s.origin)) :
          @sprintf("radius %.6g mm, center %s", s.radius, p3(s.center))
    side = r.facing == 0 ? ", mixed facing" :
           r.kind === :plane ? "" : r.facing < 0 ? ", facing in" : ", facing out"
    @printf(io, "SurfaceRegion(:%s, %s, %d patches, %.4g mm², rms %.2g mm, max %.2g mm%s)",
            r.kind, geo, length(r.patches), r.area, r.rms, r.maxdev, side)
end

# A fitted surface: plane (point c, unit normal u), cylinder (axis point c, unit axis
# u, radius R) or sphere (center c, radius R).
struct _SurfaceModel
    kind::Symbol
    c::Vec3
    u::Vec3
    R::Float64
end

function _model_distance(m::_SurfaceModel, p)
    d = p - m.c
    m.kind === :plane && return dot(d, m.u)
    m.kind === :cylinder && return norm(d - dot(d, m.u) * m.u) - m.R
    return norm(d) - m.R
end

function _model_normal(m::_SurfaceModel, p)
    m.kind === :plane && return m.u
    d = p - m.c
    m.kind === :cylinder && (d -= dot(d, m.u) * m.u)
    return normalize(d)
end

_model_surface(m::_SurfaceModel) =
    m.kind === :cylinder ? CylinderSurface(m.c, m.u, m.R) :
    m.kind === :sphere ? SphereSurface(m.c, m.R) : nothing

# Least-squares fits to points X (with the patch normals N and areas W for the
# cylinder's axis). `nothing` when the points do not determine the surface, or give a
# radius beyond 1000 times their extent (a plane in disguise).
function _fit_model(kind::Symbol, X::Vector{Vec3}, N::Vector{Vec3}, W::Vector{Float64})
    m = sum(X) / length(X)
    ext = maximum(x -> norm(x - m), X)
    ext > 0 || return nothing
    if kind === :plane
        length(X) >= 3 || return nothing
        C = sum((x - m) * (x - m)' for x in X)
        return _SurfaceModel(:plane, m, Vec3(eigen(Symmetric(Matrix(C))).vectors[:, 1]), Inf)
    elseif kind === :sphere
        length(X) >= 5 || return nothing
        return _fit_sphere(X, m, ext)
    else
        length(X) >= 6 || return nothing
        return _fit_cylinder(X, N, W, m, ext)
    end
end

# Algebraic fit (|x|² = 2c·x + k, linear in c and k), refined by Gauss–Newton on the
# geometric distances |x − c| − R.
function _fit_sphere(X, m, ext)
    M = zeros(SMatrix{4,4,Float64})
    g = zeros(SVector{4,Float64})
    for x in X
        y = x - m
        a = SVector(2y[1], 2y[2], 2y[3], 1.0)
        M += a * a'
        g += a * dot(y, y)
    end
    s = _solve_small(M, g)
    s === nothing && return nothing
    c = Vec3(s[1], s[2], s[3])
    R = sqrt(max(s[4] + dot(c, c), 0.0))
    for _ in 1:20
        (isfinite(R) && 0 < R < 1e3ext) || return nothing
        M = zeros(SMatrix{4,4,Float64})
        g = zeros(SVector{4,Float64})
        for x in X
            d = x - m - c
            l = norm(d)
            l > 0 || continue
            J = SVector(-d[1] / l, -d[2] / l, -d[3] / l, -1.0)
            M += J * J'
            g += J * (l - R)
        end
        δ = _solve_small(M, -g)
        δ === nothing && return nothing
        c += Vec3(δ[1], δ[2], δ[3])
        R += δ[4]
        norm(δ) <= 1e-13 * R && break
    end
    (isfinite(R) && 0 < R < 1e3ext) || return nothing
    return _SurfaceModel(:sphere, m + c, zero(Vec3), R)
end

# Axis from the patch normals (perpendicular to it on a cylinder: the eigenvector of
# Σ w n nᵀ with the least eigenvalue), a circle fit across the axis, then Gauss–Newton
# on the geometric distances over the axis direction, its position and the radius.
function _fit_cylinder(X, N, W, m, ext)
    S = zeros(SMatrix{3,3,Float64})
    for (n, w) in zip(N, W)
        S += w * (n * n')
    end
    a = Vec3(eigen(Symmetric(Matrix(S))).vectors[:, 1])
    e1, e2 = _perp_basis(a)
    M = zeros(SMatrix{3,3,Float64})
    g = zeros(SVector{3,Float64})
    for x in X
        y = x - m
        u, v = dot(y, e1), dot(y, e2)
        q = SVector(2u, 2v, 1.0)
        M += q * q'
        g += q * (u^2 + v^2)
    end
    s = _solve_small(M, g)
    s === nothing && return nothing
    o = s[1] * e1 + s[2] * e2                       # relative to m
    R = sqrt(max(s[3] + s[1]^2 + s[2]^2, 0.0))
    for _ in 1:20
        (isfinite(R) && 0 < R < 1e3ext) || return nothing
        e1, e2 = _perp_basis(a)
        M5 = zeros(SMatrix{5,5,Float64})
        g5 = zeros(SVector{5,Float64})
        for x in X
            d = x - m - o
            h = dot(d, a)
            q = d - h * a
            ρ = norm(q)
            ρ > 0 || continue
            û = q / ρ
            ue1, ue2 = dot(û, e1), dot(û, e2)
            J = SVector(-h * ue1, -h * ue2, -ue1, -ue2, -1.0)
            M5 += J * J'
            g5 += J * (ρ - R)
        end
        δ = _solve_small(M5, -g5)
        δ === nothing && return nothing
        a = normalize(a + δ[1] * e1 + δ[2] * e2)
        o += δ[3] * e1 + δ[4] * e2
        o -= dot(o, a) * a                          # the axis point nearest m
        R += δ[5]
        norm(δ) <= 1e-13 * max(R, 1.0) && break
    end
    (isfinite(R) && 0 < R < 1e3ext) || return nothing
    return _SurfaceModel(:cylinder, m + o, a, R)
end

function _perp_basis(a)
    e1 = normalize(cross(a, abs(a[1]) < 0.9 ? Vec3(1, 0, 0) : Vec3(0, 1, 0)))
    return e1, cross(a, e1)
end

function _solve_small(M, g)
    F = lu(M; check = false)
    issuccess(F) || return nothing
    x = F \ g
    return all(isfinite, x) ? x : nothing
end

# Per-patch data on welded vertices: positions, each patch's vertices, area, unit
# normal, centroid, and the patches sharing a vertex with it.
function _fit_data(m::Mesh, weld)
    id = _welded_ids(m; tol = weld)
    P = Vector{Vec3}(undef, isempty(id) ? 0 : maximum(id))
    for i in reverse(eachindex(id))
        P[id[i]] = m.vertices[i]
    end
    pv = [unique(id[v] for k in r for v in m.triangles[k]) for r in m.patches]
    av = [_patch_area_vector(m, r) for r in m.patches]
    # area as the sum of the facet areas: a patch folded across a curved surface has a
    # shorter area vector
    area = [sum(norm(_tri_cross(triangle_vertices(m, k)...)) for k in r; init = 0.0) / 2 for r in m.patches]
    nrm = [norm(a) > 0 ? a / norm(a) : zero(Vec3) for a in av]
    cen = [_patch_centroid(m, r) for r in m.patches]
    vp = [Int[] for _ in eachindex(P)]
    for (p, vs) in enumerate(pv), v in vs
        push!(vp[v], p)
    end
    nbr = [unique(q for v in pv[p] for q in vp[v] if q != p) for p in eachindex(pv)]
    return (; P, pv, area, nrm, cen, nbr)
end

function _fit_patches(kind, D, ps, scratch)
    X = Vec3[]
    seen = scratch
    for p in ps, v in D.pv[p]
        seen[v] && continue
        seen[v] = true
        push!(X, D.P[v])
    end
    for p in ps, v in D.pv[p]
        seen[v] = false
    end
    isempty(X) && return nothing
    return _fit_model(kind, X, D.nrm[ps], D.area[ps])
end

# Tolerance on the distance of vertices from a model: `rtol` of the radius or of the
# region's extent `ext` if smaller (a gently curved face is not a sphere of any radius
# that bends less than `rtol R` across it), at least `atol`; `atol` for a plane (CAD
# planes are flat to round-off).
_fit_tol(m::_SurfaceModel, ext, rtol, atol) = m.kind === :plane ? atol : max(atol, rtol * min(m.R, ext))

# Diagonal of the bounding box of the patches' vertices.
function _extent(D, ps)
    lo, hi = Vec3(Inf, Inf, Inf), Vec3(-Inf, -Inf, -Inf)
    for p in ps, v in D.pv[p]
        lo, hi = min.(lo, D.P[v]), max.(hi, D.P[v])
    end
    return norm(hi - lo)
end

function _patch_on(m::_SurfaceModel, D, p, tol, cosmax)
    for v in D.pv[p]
        abs(_model_distance(m, D.P[v])) <= tol || return false
    end
    D.area[p] > 0 || return true
    return abs(dot(D.nrm[p], _model_normal(m, D.cen[p]))) >= cosmax
end

# Grow a region of one kind from `seed` over available patches, nearest (to the seed)
# first, refitting as it grows; `nothing` if no start around the seed fits.
function _grow_region(kind, seed, D, avail, W, opts)
    (; rtol, atol, cosmax) = opts
    (; mark, queued, scratch, heap) = W
    W.gid[] += 1
    gid = W.gid[]
    # start from the seed and its neighbors (two rings, then one) facing its way
    near(q) = avail[q] && (D.area[q] == 0 || abs(dot(D.nrm[q], D.nrm[seed])) >= cosmax)
    ring1 = [seed; [q for q in D.nbr[seed] if near(q)]]
    ring2 = unique!([ring1; [q for p in ring1 for q in D.nbr[p] if near(q)]])
    model = nothing
    members = Int[]
    for init in (ring2, ring1, [seed])
        mdl = _fit_patches(kind, D, init, scratch)
        mdl === nothing && continue
        tol = _fit_tol(mdl, _extent(D, init), rtol, atol)
        all(p -> _patch_on(mdl, D, p, tol, cosmax), init) || continue
        model, members = mdl, copy(init)
        break
    end
    model === nothing && return nothing
    return _grow_from!(kind, seed, model::_SurfaceModel, members, D, avail, W, gid, opts)
end

function _push_neighbors!(H, D, p, avail, mark, queued, gid, c0)
    for q in D.nbr[p]
        (avail[q] && mark[q] != gid && queued[q] != gid) || continue
        queued[q] = gid
        _heap_push!(H, norm(D.cen[q] - c0), q, 0)
    end
end

function _grow_from!(kind, seed, model::_SurfaceModel, members::Vector{Int}, D, avail, W, gid, opts)
    (; rtol, atol, cosmax) = opts
    (; mark, queued, scratch, heap) = W
    H = heap
    for p in members
        mark[p] = gid
    end
    empty!(H.key); empty!(H.tri); empty!(H.seed)
    c0 = D.cen[seed]
    for p in members
        _push_neighbors!(H, D, p, avail, mark, queued, gid, c0)
    end
    lo, hi = Vec3(Inf, Inf, Inf), Vec3(-Inf, -Inf, -Inf)
    for p in members, v in D.pv[p]
        lo, hi = min.(lo, D.P[v]), max.(hi, D.P[v])
    end
    tol = _fit_tol(model, norm(hi - lo), rtol, atol)
    nfit = length(members)
    rejected = Int[]
    rounds = 0
    while true
        while !isempty(H.key)
            _, p, _ = _heap_pop!(H)
            mark[p] == gid && continue
            if _patch_on(model, D, p, tol, cosmax)
                mark[p] = gid
                push!(members, p)
                for v in D.pv[p]
                    lo, hi = min.(lo, D.P[v]), max.(hi, D.P[v])
                end
                _push_neighbors!(H, D, p, avail, mark, queued, gid, c0)
                length(members) < 1.25nfit + 4 && continue
            else
                push!(rejected, p)
                continue
            end
            # refit, and reconsider the rejected patches under the new fit
            mdl = _fit_patches(kind, D, members, scratch)
            if mdl !== nothing
                model = mdl
                tol = _fit_tol(model, norm(hi - lo), rtol, atol)
            end
            nfit = length(members)
            for q in rejected
                _heap_push!(H, norm(D.cen[q] - c0), q, 0)
            end
            empty!(rejected)
        end
        rounds += 1
        (isempty(rejected) || nfit == length(members) || rounds > 3) && break
        mdl = _fit_patches(kind, D, members, scratch)
        mdl === nothing || (model = mdl; tol = _fit_tol(model, norm(hi - lo), rtol, atol))
        nfit = length(members)
        for q in rejected
            _heap_push!(H, norm(D.cen[q] - c0), q, 0)
        end
        empty!(rejected)
    end
    # final fit; drop patches that it leaves out of tolerance
    for _ in 1:2
        mdl = _fit_patches(kind, D, members, scratch)
        mdl === nothing && return nothing
        model = mdl
        tol = _fit_tol(model, _extent(D, members), rtol, atol)
        keep = filter(p -> _patch_on(model, D, p, tol, cosmax), members)
        length(keep) == length(members) && break
        for p in members
            mark[p] = 0
        end
        members = keep
        isempty(members) && return nothing
    end
    return members, model
end

# rms and largest distance of the vertices of patches `members` from the model
function _deviation(model, D, members, scratch)
    s2, mx, n = 0.0, 0.0, 0
    for p in members, v in D.pv[p]
        scratch[v] && continue
        scratch[v] = true
        d = abs(_model_distance(model, D.P[v]))
        s2 += d^2
        mx = max(mx, d)
        n += 1
    end
    for p in members, v in D.pv[p]
        scratch[v] = false
    end
    return sqrt(s2 / max(n, 1)), mx
end

function _region(kind, members, model, D, scratch)
    A = sum(D.area[members])
    cen = A > 0 ? sum(D.area[p] * D.cen[p] for p in members) / A : sum(D.cen[members]) / length(members)
    rms, mx = _deviation(model, D, members, scratch)
    sides = [sign(dot(D.nrm[p], _model_normal(model, D.cen[p]))) for p in members if D.area[p] > 0]
    nrm = zero(Vec3)
    if kind === :plane
        w = sum(D.area[p] * dot(D.nrm[p], model.u) for p in members)
        nrm = w < 0 ? -model.u : model.u
        sides .*= w < 0 ? -1 : 1
        facing = all(>(0), sides) ? 1 : 0
    else
        facing = all(>(0), sides) ? 1 : all(<(0), sides) ? -1 : 0
    end
    return SurfaceRegion(kind, _model_surface(model), cen, nrm, sort(members), A, rms, mx, facing)
end

"""
    fit_surfaces(mesh; reference = nothing, kinds = (:plane, :cylinder, :sphere),
                 rtol = 1e-3, atol = 1e-4, angle = 30.0, minpatches = 6, minarea = 0.0,
                 tol = 1e-9) -> Vector{SurfaceRegion}

Regions of `mesh` whose patches lie on one plane, cylinder (any axis) or sphere, by
decreasing area; patches in no region are on none of them (fillets, free-form faces,
patches straddling two faces). Use it on imported CAD (raw facets or after
`remesh_mesh`), then [`curved_components`](@ref) to give the curved regions their exact
surface for the radiation model.

A patch belongs to a surface when every one of its vertices lies within
`max(atol, rtol min(R, L))` mm of it (R the radius, L the region's extent; `atol` for a
plane: CAD planes are flat to round-off) and its normal is within `angle` degrees of the surface normal at its centroid (either
side: windings do not matter). Regions grow from seed patches, by decreasing area, over
patches that share a vertex (welded by position to `tol` times the mesh's extent),
nearest first, refitting by least squares as they grow: planes by principal components,
spheres by an algebraic fit, cylinders by the axis the patch normals are perpendicular
to and a circle fit across it, both refined by Gauss–Newton on the vertex distances.
From each seed all `kinds` are grown and the largest region kept (a plane, then a
cylinder, when within 1% of the largest: a narrow band with vertices only on its two
rims, such as the bore of an entrance hole, fits a sphere as well as a cylinder) among
those with at least `minpatches` patches and `minarea` mm² whose rms distance is under
a quarter of the tolerance: on a true surface the residual is the faceting, well
inside the tolerance, while a fit to a free-form or conical face grows until the
tolerance stops it (rms about 0.4 of it), so such patches stay unfitted. Regions of one kind whose patches all lie on each other's fits (one
surface interrupted by another face) are then merged.

The tolerance must exceed the faceting several times over. Facet vertices of exported
CAD lie on the true surface, but points `remesh_mesh` adds on a facet lie up to w²/8R
inside a curved one (w the facet's width: 0.1 µm for 50 µm strips on a 2.9 mm wall,
5 µm for 0.3 mm strips), more than the tolerance allows once w exceeds ~0.08 R. So
give the facets the mesh was remeshed from as `reference`: the regions are then fitted
on them, and each patch of `mesh` joins the region whose surface all its vertices lie
within that tolerance plus the largest distance of the region's facets from it (and
whose bounding box holds them), with its normal within `angle`. The regions returned
index the patches of `mesh`, with their rms and largest distances from the
(reference) fit.
"""
function fit_surfaces(m::Mesh; reference::Union{Nothing,Mesh} = nothing,
                      kinds = (:plane, :cylinder, :sphere), rtol::Real = 1e-3, atol::Real = 1e-4,
                      angle::Real = 30.0, minpatches::Integer = 6, minarea::Real = 0.0, tol::Real = 1e-9)
    all(in((:plane, :cylinder, :sphere)), kinds) ||
        throw(ArgumentError("kinds must be among :plane, :cylinder and :sphere"))
    opts = (; rtol = Float64(rtol), atol = Float64(atol), cosmax = cosd(angle))
    if reference !== nothing
        regs = fit_surfaces(reference; kinds, rtol, atol, angle, minpatches, minarea, tol)
        return _transfer_regions(regs, reference, m, opts, tol)
    end
    np = npatches(m)
    np == 0 && return SurfaceRegion[]
    D = _fit_data(m, tol)
    W = (; mark = zeros(Int, np), queued = zeros(Int, np), gid = Ref(0), heap = _GrowHeap(),
         scratch = falses(length(D.P)))
    avail = trues(np)
    order = sortperm(D.area; rev = true)
    ordered = [k for k in (:plane, :cylinder, :sphere) if k in kinds]
    found = Tuple{Symbol,Vector{Int},_SurfaceModel}[]
    for s in order
        avail[s] || continue
        best = nothing
        bestA = 0.0
        for kind in ordered
            g = _grow_region(kind, s, D, avail, W, opts)
            g === nothing && continue
            ps, mdl = g
            A = sum(D.area[ps])
            length(ps) >= minpatches && A >= minarea || continue
            # a surface explains its vertices well within the tolerance; a fit that only
            # grew until the tolerance stopped it (a sphere on a free-form face) does not
            _deviation(mdl, D, ps, W.scratch)[1] <= _fit_tol(mdl, _extent(D, ps), opts.rtol, opts.atol) / 4 || continue
            if best === nothing || A > 1.01bestA
                best, bestA = (kind, ps, mdl), A
            end
        end
        best === nothing && continue
        avail[best[2]] .= false
        push!(found, best)
    end
    # merge regions of one kind lying on one surface
    merged = true
    while merged
        merged = false
        for i in eachindex(found), j in i+1:length(found)
            (ki, psi, mi), (kj, psj, mj) = found[i], found[j]
            ki === kj || continue
            ti = _fit_tol(mi, _extent(D, psi), opts.rtol, opts.atol)
            tj = _fit_tol(mj, _extent(D, psj), opts.rtol, opts.atol)
            all(p -> _patch_on(mi, D, p, ti, opts.cosmax), psj) &&
                all(p -> _patch_on(mj, D, p, tj, opts.cosmax), psi) || continue
            ps = [psi; psj]
            mdl = _fit_patches(ki, D, ps, W.scratch)
            mdl === nothing && continue
            t = _fit_tol(mdl, _extent(D, ps), opts.rtol, opts.atol)
            all(p -> _patch_on(mdl, D, p, t, opts.cosmax), ps) || continue
            found[i] = (ki, ps, mdl)
            deleteat!(found, j)
            merged = true
            break
        end
    end
    regions = [_region(k, ps, mdl, D, W.scratch) for (k, ps, mdl) in found]
    return sort!(regions; by = r -> -r.area)
end

_region_model(r::SurfaceRegion) =
    r.kind === :plane ? _SurfaceModel(:plane, r.centroid, r.normal, Inf) :
    r.kind === :cylinder ? _SurfaceModel(:cylinder, r.surface.origin, r.surface.axis, r.surface.radius) :
    _SurfaceModel(:sphere, r.surface.center, zero(Vec3), r.surface.radius)

# Largest distance from the model of the points of the triangle (a, b, c), sampled on a
# barycentric grid of step 1/4 (on a cylinder or sphere the distance is concave across a
# facet whose vertices lie on it, largest inside or on an edge).
function _facet_deviation(m::_SurfaceModel, a, b, c)
    d = 0.0
    for i in 0:4, j in 0:4-i
        d = max(d, abs(_model_distance(m, a + (i / 4) * (b - a) + (j / 4) * (c - a))))
    end
    return d
end

# Regions fitted on `ref` carried over to the patches of `m`, a remeshing of the same
# surface (see `fit_surfaces`).
function _transfer_regions(regs, ref::Mesh, m::Mesh, opts, weld)
    np = npatches(m)
    (np == 0 || isempty(regs)) && return SurfaceRegion[]
    D = _fit_data(m, weld)
    Dr = _fit_data(ref, weld)
    models = [_region_model(r) for r in regs]
    lim = Float64[]
    boxes = Tuple{Vec3,Vec3}[]
    for (r, mdl) in zip(regs, models)
        sag = 0.0
        for p in r.patches, k in ref.patches[p]
            sag = max(sag, _facet_deviation(mdl, triangle_vertices(ref, k)...))
        end
        l = _fit_tol(mdl, _extent(Dr, r.patches), opts.rtol, opts.atol) + 1.25sag
        lo, hi = Vec3(Inf, Inf, Inf), Vec3(-Inf, -Inf, -Inf)
        for p in r.patches, v in Dr.pv[p]
            lo, hi = min.(lo, Dr.P[v]), max.(hi, Dr.P[v])
        end
        push!(lim, l)
        push!(boxes, (lo .- l, hi .+ l))
    end
    members = [Int[] for _ in regs]
    for p in 1:np
        best, score = 0, Inf
        for (k, mdl) in enumerate(models)
            lo, hi = boxes[k]
            all(v -> all(lo .<= D.P[v] .<= hi), D.pv[p]) || continue
            d = maximum(v -> abs(_model_distance(mdl, D.P[v])), D.pv[p])
            d <= lim[k] && d / lim[k] < score || continue
            D.area[p] > 0 && abs(dot(D.nrm[p], _model_normal(mdl, D.cen[p]))) < opts.cosmax && continue
            best, score = k, d / lim[k]
        end
        best > 0 && push!(members[best], p)
    end
    scratch = falses(length(D.P))
    out = [_region(r.kind, ps, mdl, D, scratch) for (r, mdl, ps) in zip(regs, models, members) if !isempty(ps)]
    return sort!(out; by = r -> -r.area)
end

# The patches `ps` of a mesh as a mesh of their own (only the vertices they use).
function _patch_submesh(m::Mesh, ps)
    vmap = zeros(Int, length(m.vertices))
    verts = Vec3[]
    tris = NTuple{3,Int}[]
    patches = UnitRange{Int}[]
    idx(v) = (vmap[v] == 0 && (push!(verts, m.vertices[v]); vmap[v] = length(verts)); vmap[v])
    for p in ps
        t0 = length(tris) + 1
        for k in m.patches[p]
            a, b, c = m.triangles[k]
            push!(tris, (idx(a), idx(b), idx(c)))
        end
        push!(patches, t0:length(tris))
    end
    return Mesh(verts, tris, patches)
end

"""
    curved_components(component; kinds = (:cylinder, :sphere), regions = nothing, kwargs...)
        -> Vector{Component}
    curved_components(components; kwargs...) -> Vector{Component}
    curved_components(target; kwargs...) -> Target

The component split into one component per cylinder or sphere that
[`fit_surfaces`](@ref) finds in its mesh (keywords are passed on to it; or give its
`regions`), each carrying its fitted `CylinderSurface` or `SphereSurface`, and the rest:
`<name>/cyl1`, `<name>/cyl2`, …, `<name>/sph1`, …, by decreasing area, then
`<name>/rest` with the other patches (planes and unrecognized faces). Material, kind
and sidedness are the component's. For a remeshed CAD part pass the facets it was
remeshed from as `reference` (see `fit_surfaces`). A component that already has a surface, or in which
no curved region is found, is returned unchanged.

The radiation model then integrates over the true curved surface, as it does for the
primitive `hohlraum`: each element's sub-triangle points are projected onto the surface
with exact curved-patch areas, whatever the element's shape (the compact,
many-triangle elements of `remesh_mesh` included), and the facets of a cylinder or
sphere no longer occlude pairs of points inside it. Ray casting keeps the facets.
Planar regions are not given a surface: flat elements are already exact there
(`DiskSurface` describes polar-sector elements of primitive end caps, not general
planar CAD faces).

A region of a CAD part is named by its part and kind, so `elements(target, :part)` no
longer exists after the split; collect `elements(target, n)` over the names starting
with the part's.
"""
function curved_components(c::Component; kinds = (:cylinder, :sphere), regions = nothing, kwargs...)
    all(in((:cylinder, :sphere)), kinds) || throw(ArgumentError("kinds must be among :cylinder and :sphere"))
    c.surface === nothing || return [c]
    regs = regions === nothing ? fit_surfaces(c.mesh; kwargs...) : regions
    sel = sort!([r for r in regs if r.kind in kinds && r.surface !== nothing]; by = r -> -r.area)
    isempty(sel) && return [c]
    out = Component[]
    used = falses(npatches(c.mesh))
    nk = Dict{Symbol,Int}()
    for r in sel
        any(used[r.patches]) && throw(ArgumentError("regions overlap"))
        k = nk[r.kind] = get(nk, r.kind, 0) + 1
        name = Symbol(c.name, r.kind === :cylinder ? "/cyl" : "/sph", k)
        push!(out, Component(name, _patch_submesh(c.mesh, r.patches), c.material; kind = c.kind,
                             two_sided = c.two_sided, surface = r.surface))
        used[r.patches] .= true
    end
    rest = findall(!, used)
    isempty(rest) || push!(out, Component(Symbol(c.name, "/rest"), _patch_submesh(c.mesh, rest), c.material;
                                          kind = c.kind, two_sided = c.two_sided))
    return out
end

curved_components(cs::AbstractVector{Component}; kwargs...) =
    reduce(vcat, (curved_components(c; kwargs...) for c in cs); init = Component[])
curved_components(t::Target; kwargs...) = Target(curved_components(t.components; kwargs...))
