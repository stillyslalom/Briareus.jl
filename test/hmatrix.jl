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
    m = RadiationModel(t)                                     # dense Float64, LU
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

    # the automatic choice above Briareus.HMATRIX_MODEL[]
    old = Briareus.HMATRIX_MODEL[]
    try
        Briareus.HMATRIX_MODEL[] = N - 1
        @test RadiationModel(t).F isa HMatrix
        @test RadiationModel(t; storage = Float32).F isa Float32Matrix
    finally
        Briareus.HMATRIX_MODEL[] = old
    end
    @test RadiationModel(t).F isa Matrix{Float64}
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
    cad = Target([part(:top, cad_half(; R, Rout, Rleh, h, t = th, nphi = 360), 0.6),
                  part(:bottom, cad_mirror(cad_half(; R, Rout, Rleh, h, t = th, nphi = 360)), 0.6),
                  part(:capsule, cad_sphere(rc; ntheta = 56, nphi = 112), 0.3),
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
