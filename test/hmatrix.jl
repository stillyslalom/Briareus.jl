# Hierarchical (H-matrix) view factors, forced on at small N against the dense path.

@testset "hierarchical view factors" begin
    wall = Material(:Au; albedo = WallLossAlbedo(), conversion = 0.8)
    t = hohlraum(; h = 0.45, capsule_h = 0.2, wall)
    N = nelements(t)
    A = t.elem_area
    rel(x, y) = maximum(abs.(x .- y)) / maximum(abs.(y))
    F = view_factors(t)
    # small leaves and a loose tolerance, so that low-rank, zero and dense blocks all occur
    kw = (; tol = 1e-4, eta = 1.0, leaf = 16)
    H = view_factors(t; storage = :hmatrix, kw...)
    @test H isa HMatrix && size(H) == (N, N)
    @test !isempty(H.lowrank) && !isempty(H.dense) && !isempty(H.scattered) && H.info.facing > 0
    @test H.info.masked > 0                       # shafts the capsule's balls decide
    # no cluster mixes surfaces (capsule, wall, a cap with its LEH disc)
    g = Briareus._hsurface_groups(t)
    @test length(unique(g)) == 4
    nodes, perm = Briareus._hcluster_tree(t.elem_centroid, ones(N), 16, g)
    @test all(nd -> nd.kids != (0, 0) || allequal(g[perm[nd.lo:nd.hi]]), nodes)
    @test H.info.bytes < 4 * N^2
    @test occursin("HMatrix", sprint(show, MIME"text/plain"(), H))
    M = Matrix(H)
    @test rel(M, F) < 2e-4                                   # tol 1e-4 per block
    @test isapprox(H[3, 7], M[3, 7]; rtol = 1e-12) && isapprox(H[N, 1], M[N, 1]; rtol = 1e-12)
    # reciprocity: G = A F is stored once per pair
    G = A .* M
    @test maximum(abs.(G .- G')) <= 1e-12 * maximum(abs, G)

    # products against the expanded matrix
    rng = Xoshiro(11)
    x, X, L = rand(rng, N), rand(rng, N, 5), rand(rng, 3, N)
    @test rel(H * x, M * x) < 1e-12
    @test rel(H' * x, M' * x) < 1e-12 && rel(transpose(H) * x, M' * x) < 1e-12
    @test rel(H * X, M * X) < 1e-12 && rel(H' * X, M' * X) < 1e-12
    @test rel(L * H, L * M) < 1e-12
    y = rand(rng, N); y0 = copy(y)
    mul!(y, H, x, 2.0, 0.5)
    @test rel(y, 2 .* (M * x) .+ 0.5 .* y0) < 1e-12
    xd = [ForwardDiff.Dual(x[i], 1.0, x[i]^2) for i in eachindex(x)]
    yd = H * xd
    @test rel(ForwardDiff.value.(yd), M * x) < 1e-12
    @test rel([ForwardDiff.partials(d)[2] for d in yd], M * x .^ 2) < 1e-12

    # the model: balanced through products (the scaling only), iterative
    m = RadiationModel(t; F, closed = true)                   # dense Float64, LU (≡ RadiationModel(t))
    mh = RadiationModel(t; storage = :hmatrix, kw...)
    @test mh.F isa HMatrix && mh.iterative
    @test sprint(show, mh) == "RadiationModel($N elements, hierarchical, iterative)"
    @test maximum(abs, closure_error(mh.F)) < 1e-10
    @test mh.F.lowrank === H.lowrank || mh.F.lowrank == H.lowrank
    Gh = A .* Matrix(mh.F)
    @test maximum(abs.(Gh .- Gh')) <= 1e-12 * maximum(abs, Gh)
    @test_throws ArgumentError RadiationModel(t; storage = :hmatrix, solver = :direct)
    @test_throws ArgumentError RadiationModel(t; F = m.F, storage = :hmatrix)
    @test_throws ArgumentError view_factors(t; tol = 1e-6)          # H-matrix keywords only
    @test RadiationModel(t; F = mh.F, storage = Float64, closed = false).F ≈ Matrix(mh.F)

    q = [remake(b; aim = (0.0, 0.0, sign(b.port[3]) * 4.0)) for b in nif_quads()]
    dep = deposit(DepositionOperator(t, q), q)
    cap = elements(t, :capsule)
    s, sh = solve(m, dep; time = 2.0), solve(mh, dep; time = 2.0)
    @test rel(sh.H, s.H) < 1e-5
    a, ah = legendre_moments(s.H, t, cap; lmax = 4), legendre_moments(sh.H, t, cap; lmax = 4)
    @test abs(ah[1] / a[1] - 1) < 1e-5
    @test abs(ah[3] / ah[1] - a[3] / a[1]) < 1e-5 && abs(ah[5] / ah[1] - a[5] / a[1]) < 1e-5
    eb = energy_balance(sh)
    @test abs(eb.lost) < 1e-8 * eb.source
    W = legendre_weights(t, cap; lmax = 4)
    α = s.albedo
    @test rel(flux_response(mh, cap; albedo = α, weights = W), flux_response(m, cap; albedo = α, weights = W)) < 1e-4
    @test rel(flux_response(mh, cap[1:5]; albedo = α), flux_response(m, cap[1:5]; albedo = α)) < 1e-5
    # derivatives through the iterative solve
    z = [c[3] for c in t.elem_centroid] ./ 5
    f(mod) = p -> sum(solve(mod, (1 .+ p[1] .* z .+ p[2] .* z .^ 2) .* dep.element; albedo = α).H[cap])
    @test rel(ForwardDiff.gradient(f(mh), [0.1, 0.2]), ForwardDiff.gradient(f(m), [0.1, 0.2])) < 1e-5

    # the automatic choice above Briareus.HMATRIX_MODEL[] (on a coarse target: only the
    # storage is checked)
    ts = hohlraum(; h = 1.2, capsule_h = 0.6, wall)
    old = Briareus.HMATRIX_MODEL[]
    try
        Briareus.HMATRIX_MODEL[] = nelements(ts) - 1
        @test RadiationModel(ts).F isa HMatrix
        @test RadiationModel(ts; storage = Float32).F isa Float32Matrix
    finally
        Briareus.HMATRIX_MODEL[] = old
    end
    @test RadiationModel(ts).F isa Matrix{Float64}
end

# Rows that a facing test cut, refined over the column tree into indexed low-rank blocks:
# the entries (read block by block) and the products against the dense view factors.
@testset "hierarchical view factors: cut rows" begin
    wall = Material(:Au; albedo = WallLossAlbedo(), conversion = 0.8)
    t = hohlraum(; h = 0.35, capsule_h = 0.15, wall)
    N, A = nelements(t), t.elem_area
    rel(x, y) = maximum(abs.(x .- y)) / maximum(abs.(y))
    H = view_factors(t; storage = :hmatrix, tol = 1e-4, eta = 1.0, leaf = 16)
    @test !isempty(H.lowrank_index)
    @test sum(size(U, 2) for (_, _, U, _) in H.lowrank_index) > 0
    M = Matrix(H)                                    # through the products
    @test rel(M, view_factors(t)) < 2e-4
    # entries through `getindex` (the stored blocks, not the products), in each indexed block
    ij = [(H.perm[a], H.perm[b]) for (rs, cs, _, _) in H.lowrank_index for a in rs for b in cs]
    @test isapprox([H[i, j] for (i, j) in ij], [M[i, j] for (i, j) in ij]; rtol = 1e-12)
    G = A .* M
    @test maximum(abs.(G .- G')) <= 1e-12 * maximum(abs, G)
    x = rand(Xoshiro(3), N)
    @test rel(H' * x, M' * x) < 1e-12
end

# The H-matrix's vectorized sector contour: its branch-free atan against `Base.atan`, and
# sector factors against the scalar contour (`view_factors`' kernel) from points all
# around the end caps, in front of, behind and cutting the sectors.
@testset "hierarchical view factors: vectorized sector contour" begin
    rng = Xoshiro(5)
    y = abs.(randn(rng, 10^5)) .* exp10.(rand(rng, -4:4, 10^5))
    x = randn(rng, 10^5) .* exp10.(rand(rng, -4:4, 10^5))
    @test all(abs.(Briareus._atan_nb.(y, x) .- atan.(y, x)) .<= 2 .* eps.(atan.(y, x)))
    @test Briareus._atan_nb.([1.0, 0.0, 0.0, 1.0], [0.0, -1.0, 1.0, 1.0]) == [π / 2, Float64(π), 0, π / 4]
    t = hohlraum(; h = 0.6, capsule_h = 0.3)
    rs = [Briareus._receiver(t, e, true, false, nothing) for e in 1:nelements(t)
          if t.components[t.elem_comp[e]].surface isa DiskSurface && t.elem_curved[e] != 0]
    @test length(rs) > 10
    ps = [Vec3(4randn(rng), 4randn(rng), 6randn(rng)) for _ in 1:200]
    ns = [normalize(Vec3(randn(rng), randn(rng), randn(rng))) for _ in 1:200]
    f1 = [Briareus._receiver_vf(p, n, r) for r in rs for (p, n) in zip(ps, ns)]
    f2 = [Briareus._receiver_vf(p, n, Briareus._soa_receiver(r)) for r in rs for (p, n) in zip(ps, ns)]
    @test count(>(0), f1) > length(f1) ÷ 10
    # (absolute below ~1e-6: the contour's terms cancel there)
    @test all(abs.(f2 .- f1) .<= 1e-12 .* abs.(f1) .+ 1e-15)
    # the adaptive rule's distance test by the loop: the triangles' verdict, also at
    # distances just either side of the exact one and in the sector's plane
    qs = [(r, v + Vec3(randn(rng), randn(rng), randn(rng)) * exp10(rand(rng, -4:0)))
          for r in rs for v in r.loop[1:4:end]]
    qs = [qs; [(r, p - dot(p - r.o, r.nr) * r.nr) for (r, p) in qs[1:2:end]]]
    D = [Briareus._distance_to(p, r) for (r, p) in qs]
    @test all(Briareus._closer_than(p, r, D[k] * c) == (D[k] < D[k] * c)
              for (k, (r, p)) in enumerate(qs) for c in (0.5, 1 - 1e-7, 1.0, 1 + 1e-7, 2.0))
end

# The visibility rules (shaft culling without rays, the convex and planar exclusions, the
# exact masks of occupied leaf pairs) on targets with rims, gaps, bores, corner fillets,
# folded CAD elements and partial visibility: every entry as `view_factors` has it, to
# the Float32 rounding and the cross-approximation tolerance.
@testset "hierarchical view factors: visibility rules" begin
    R, Rout, Rleh, h, th, rc = 2.875, 2.975, 1.55, 4.7, 0.05, 0.9
    cavity(m) = split_mesh(m, c -> norm(c[1:2]) < R + 0.01 && abs(c[3]) < h + 1e-6)[1]
    part(name, F, size) = Component(name, cavity(remesh_mesh(orient_mesh(facet_soup(F)); size)),
                                    name === :capsule ? ablator() : gold())
    # (coarser than partial_visibility.jl's: the sampled factors' cost is in the facets)
    cad = Target([part(:top, cad_half(; R, Rout, Rleh, h, t = th, nphi = 180), 1.0),
                  part(:bottom, cad_mirror(cad_half(; R, Rout, Rleh, h, t = th, nphi = 180)), 1.0),
                  part(:capsule, cad_sphere(rc; ntheta = 28, nphi = 56), 0.5),
                  aperture(:leh_top, disk_mesh(Rleh, h; nr = 3, nphi = 24, up = false)),
                  aperture(:leh_bottom, disk_mesh(Rleh, -h; nr = 3, nphi = 24, up = true))]; visibility = :sampled)
    bored = bored_cavity(; hw = 0.8, hc = 0.25)
    for t in (bored, Target(bored.components; visibility = :sampled), cad,
              hohlraum(; h = 0.7, capsule_h = 0.2, capsule_mesh = :latlong),
              hohlraum(; h = 0.7, capsule_h = 0.2, corner_radius = 0.5))
        F = view_factors(t)
        H = view_factors(t; storage = :hmatrix, eta = 1.0, leaf = 16)
        A = t.elem_area
        G, Gh = A .* F, A .* Matrix(H)
        @test maximum(abs.(Gh .- G)) < 1e-5 * maximum(G)
        @test maximum(abs.(closure_error(H) .- closure_error(F))) < 1e-5
    end
end

# Visibility certificates of occupied shafts (a planar part's facets, a cylinder's
# axisymmetric bounds, a convex part's balls and its elements' planes) against the
# visibility ray itself: every segment a certificate decides is decided as the ray
# against that part's facets alone finds it (random segments, and segments crossing a
# washer's plane within 1e-3 to 1e-12 of its rims); the certified masks of every leaf
# pair equal the masks by a ray per pair, and the pairs of empty shafts (those the facing
# slack admits included) get `_vf_pair`'s factor exactly.
@testset "hierarchical view factors: visibility certificates" begin
    B = Briareus
    rng = Xoshiro(11)
    for (t, nseg) in ((bored_cavity(; hw = 0.8, hc = 0.25), 3000), (hohlraum(; h = 0.7, capsule_h = 0.2), 3000),
                      (hohlraum(; h = 0.7, capsule_h = 0.2, capsule_mesh = :latlong), 0))
        ctx0 = B._vf_context(t; near = 3.0, mid = 12.0, level = 1, exact_far = true, visibility = :centroid,
                             visibility_grid = 4)
        ctx, nodes, perm, far, near, opts = B._hsetup(ctx0, 1e-6, 1.0, 16, 32)
        agree = Bool[]
        tried = 0
        for c in eachindex(t.components)
            (opts.cert[c] != 0 || opts.ball[c] != 0) && nseg > 0 || continue
            e0 = findfirst(e -> t.elem_comp[e] != c, 1:nelements(t))     # ends on no element of c
            ks = findall(==(c), t.tri_comp)
            vs = [v for k in ks for v in (t.tri_v0[k], t.tri_v0[k] + t.tri_e1[k], t.tri_v0[k] + t.tri_e2[k])]
            lo, hi = reduce((a, b) -> min.(a, b), vs), reduce((a, b) -> max.(a, b), vs)
            mid, ext = (lo + hi) / 2, (hi - lo) / 2 .+ 0.5
            segs = [(mid + ext .* (2rand(rng, B.Vec3) .- 1), mid + ext .* (2rand(rng, B.Vec3) .- 1)) for _ in 1:nseg]
            k = opts.cert[c]
            if k > 0
                # through the plane near the inner and outer rims
                P = opts.planar[k]
                for _ in 1:nseg
                    ρ = rand(rng, (P.rlo, P.rhi)) * (1 + rand(rng, (-1, 1)) * 10.0^(-rand(rng, 3:12)))
                    φ = 2π * rand(rng)
                    x = P.c + ρ * (cos(φ) * P.u + sin(φ) * P.v)
                    w = normalize(randn(rng, B.Vec3) + 2 * P.n) * (0.1 + rand(rng))
                    push!(segs, (x - rand(rng) * w, x + rand(rng) * w))
                end
            end
            for (oi, oj) in segs
                v = B._hcert_verdict(ctx, opts, c, e0, e0, oi, oj)
                tried += 1
                v < 0 && continue
                hit = any_hit(t.bvh, oi, oj - oi, 1 - 1e-9, id -> t.tri_comp[id] != c)
                push!(agree, (v == 0) == hit)
            end
        end
        nseg > 0 && @test length(agree) > 0.5 * tried && all(agree)
        # leaf pairs and the larger blocks that take masks: certified masks against ray
        # masks, empty shafts against `_vf_pair`
        masks, free = Bool[], Bool[]
        for (s, u) in [near; far]
            a, b = nodes[s], nodes[u]
            rows, cols = a.lo:a.hi, b.lo:b.hi
            B._hleaf(a) && B._hleaf(b) || length(rows) * length(cols) <= opts.mask_size || continue
            occ = B._shaft_occluder(ctx, opts, s, u, true)
            if occ.n > 0
                push!(masks, B._hmask(ctx, opts, perm, rows, cols, occ).v == B._hmask_rays(ctx, perm, rows, cols, rows, cols).v)
            elseif occ.n == 0
                for p in rows, q in cols
                    i, j = minmax(perm[p], perm[q])
                    i == j || push!(free, B._vf_pair_free(ctx, i, j, occ) === B._vf_pair(ctx, i, j))
                end
            end
        end
        @test all(masks) && all(free) && (nseg == 0 || !isempty(masks))
    end
end

# `:sampled` visibility: bundle certificates (a part's verdict on every segment between
# two balls) against the ray of each of a bundle's random segments, its ends in the balls
# (some on their surfaces), against that part's facets alone; and the certified fractions
# of every block that takes them, and the factors of empty shafts, against `_vf_visible`
# and `_vf_pair` (`===`), on the bored cavity (LEH bores, caps, a waist gap with washers).
@testset "hierarchical view factors: sampled visibility bundles" begin
    B = Briareus
    rng = Xoshiro(13)
    t = Target(bored_cavity(; hw = 0.8, hc = 0.25).components; visibility = :sampled)
    ctx0 = B._vf_context(t; near = 3.0, mid = 12.0, level = 1, exact_far = true, visibility = :sampled,
                         visibility_grid = 4)
    ctx, nodes, perm, far, near, opts = B._hsetup(ctx0, 1e-6, 1.0, 16, 32)
    agree, nv = Bool[], [0, 0, 0]
    for c in eachindex(t.components)
        opts.cert[c] != 0 || opts.ball[c] != 0 || continue
        e0 = findfirst(e -> t.elem_comp[e] != c, 1:nelements(t))
        ks = findall(==(c), t.tri_comp)
        vs = [v for k in ks for v in (t.tri_v0[k], t.tri_v0[k] + t.tri_e1[k], t.tri_v0[k] + t.tri_e2[k])]
        lo, hi = reduce((a, b) -> min.(a, b), vs), reduce((a, b) -> max.(a, b), vs)
        mid, ext = (lo + hi) / 2, (hi - lo) / 2 .+ 0.5
        ball(o, r) = o + r * (rand(rng) < 0.3 ? 1.0 : cbrt(rand(rng))) * normalize(randn(rng, B.Vec3))
        for _ in 1:1500
            oi, oj = mid + ext .* (2rand(rng, B.Vec3) .- 1), mid + ext .* (2rand(rng, B.Vec3) .- 1)
            ri, rj = 0.2 * rand(rng), 0.2 * rand(rng)
            v = B._hcert_verdict(ctx, opts, c, e0, e0, oi, oj, ri, rj)
            nv[v+2] += 1
            v < 0 && continue
            for _ in 1:8
                p, q = ball(oi, ri), ball(oj, rj)
                push!(agree, (v == 0) == any_hit(t.bvh, p, q - p, 1 - 1e-9, id -> t.tri_comp[id] != c))
            end
        end
    end
    @test nv[2] > 0 && nv[3] > 0 && all(agree)              # (hidden and clear bundles)
    frac, free = Bool[], Bool[]
    npart = 0
    for (s, u) in [near; far]
        a, b = nodes[s], nodes[u]
        rows, cols = a.lo:a.hi, b.lo:b.hi
        B._hleaf(a) && B._hleaf(b) || length(rows) * length(cols) <= opts.mask_size || continue
        occ = B._shaft_occluder(ctx, opts, s, u, true)
        if occ.n > 0
            F = B._hfrac(ctx, opts, perm, rows, cols, occ).v
            R = [first(B._vf_visible(ctx, minmax(perm[p], perm[q])...)) for p in rows, q in cols]
            push!(frac, all(F .=== R))
            npart += count(f -> 0 < f < 1, R)
        elseif occ.n == 0
            for p in rows, q in cols
                i, j = minmax(perm[p], perm[q])
                i == j || push!(free, B._vf_pair_free(ctx, i, j, occ) === B._vf_pair(ctx, i, j))
            end
        end
    end
    @test npart > 1000 && !isempty(free) && all(frac) && all(free)
end

# The entries a low-rank block's check adds where G may jump (`_htier_entries`): pairs
# off the most common quadrature tier of their row or column, the fewest-off lines
# first. On a latitude–longitude capsule (polar slivers, several sizes larger than their
# neighbours) against the wall: every entry returned is off its row's or column's
# majority tier; every pair that is alone off its line's majority is returned when the
# cap (m + n) is not reached; and a block without entries has a single tier.
@testset "hierarchical view factors: quadrature-tier entries" begin
    B = Briareus
    t = hohlraum(; h = 0.5, capsule_h = 0.12, capsule_mesh = :latlong)
    ctx0 = B._vf_context(t; near = 3.0, mid = 12.0, level = 1, exact_far = true, visibility = :centroid,
                         visibility_grid = 4)
    ctx, nodes, perm, far, near, opts = B._hsetup(ctx0, 1e-6, 1.0, 16, 32)
    tier(p, q) = (i = perm[p]; j = perm[q]; d = ctx.C[j] - ctx.C[i]; B._vf_tier(ctx, i, j, dot(d, d)))
    # (a tie goes to the nearer tier)
    majority(v) = (ts = (:near, :mid, :far); ts[argmax([count(==(x), v) for x in ts])])
    offline, alone, screened = Bool[], Bool[], Bool[]
    nent = 0
    for (s, u) in far
        a, b = nodes[s], nodes[u]
        rows, cols = a.lo:a.hi, b.lo:b.hi
        length(rows) * length(cols) <= 4096 || continue
        T = [tier(p, q) for p in rows, q in cols]
        ent = B._htier_entries(ctx, perm, rows, cols)
        nent += length(ent)
        if isempty(ent)
            push!(screened, all(==(T[1]), T))
            continue
        end
        for (x, y) in ent
            push!(offline, T[x, y] != majority(T[x, :]) || T[x, y] != majority(T[:, y]))
        end
        length(ent) < length(rows) + length(cols) || continue
        for x in axes(T, 1), y in axes(T, 2)
            lone = (count(!=(majority(T[x, :])), T[x, :]) == 1 && T[x, y] != majority(T[x, :])) ||
                   (count(!=(majority(T[:, y])), T[:, y]) == 1 && T[x, y] != majority(T[:, y]))
            lone && push!(alone, (x, y) in ent)
        end
    end
    @test nent > 0 && !isempty(alone) && all(offline) && all(alone)
    # (a block with no entries, screened by its boxes or with no line mixing tiers, has one)
    @test all(screened)
end

# Keyframes of hierarchical models: the interpolated exchange as the two operators.
@testset "hierarchical keyframes" begin
    kw = (length = 10.1, leh_radius = 1.55, capsule_radius = 0.9, h = 0.8, capsule_h = 0.4,
          wall = Material(:Au; albedo = 0.9, conversion = 0.8))
    build = moving_hohlraum(WallMotion(speed = 0.1); radius = 2.875, kw...)
    quads = nif_quads()[1:8]
    kd = Keyframes(build, [0.0, 4.0], quads)
    kh = Keyframes(build, [0.0, 4.0], quads; viewfactor_kw = (; storage = :hmatrix, eta = 1.0, leaf = 16))
    md, mh = model_at(kd, 2.5), model_at(kh, 2.5)
    @test mh.iterative && sprint(show, mh) == "RadiationModel($(nelements(mh.target)) elements, hierarchical, iterative)"
    M = Matrix(mh.F)
    @test maximum(abs.(M .- md.F)) < 1e-5 * maximum(md.F)
    @test maximum(abs, closure_error(mh.F)) < 1e-10
    A = mh.target.elem_area
    @test maximum(abs.(A .* M .- (A .* M)')) < 1e-12 * maximum(A .* M)
    dep = deposit(kd, quads, 2.5)
    rel(x, y) = maximum(abs.(x .- y)) / maximum(abs.(y))
    @test rel(solve(mh, dep).H, solve(md, dep).H) < 1e-5
end
