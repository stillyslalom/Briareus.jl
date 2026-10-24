# Large radiation models: Float32 view-factor storage and iterative (factorization-free)
# solves, forced on at small N against the dense Float64 path.

@testset "Float32Matrix products" begin
    rng = Xoshiro(7)
    D = rand(rng, Float32, 300, 200)
    A = Float32Matrix(D)
    Ad = Float64.(D)
    rel(x, y) = maximum(abs.(x .- y)) / maximum(abs.(y))
    @test size(A) == (300, 200) && A[3, 5] === Float64(D[3, 5]) && A == Float32Matrix(copy(D))
    x, z = rand(rng, 200), rand(rng, 300)
    X, Z, L = rand(rng, 200, 7), rand(rng, 300, 5), rand(rng, 4, 300)
    @test rel(A * x, Ad * x) < 1e-14
    @test rel(A' * z, Ad' * z) < 1e-14
    @test rel(transpose(A) * z, Ad' * z) < 1e-14
    @test rel(A * X, Ad * X) < 1e-14
    @test rel(A' * Z, Ad' * Z) < 1e-14
    @test rel(L * A, L * Ad) < 1e-14
    @test rel(A * X[:, 1:1], Ad * X[:, 1:1]) < 1e-14
    y = rand(rng, 300)
    y2 = copy(y)
    mul!(y, A, x, 2.0, 0.5)
    @test rel(y, 2.0 .* (Ad * x) .+ 0.5 .* y2) < 1e-14
    # first-order duals: values and partials through one product
    xd = [ForwardDiff.Dual(x[i], 1.0, x[i]^2) for i in eachindex(x)]
    yd = A * xd
    @test rel(ForwardDiff.value.(yd), Ad * x) < 1e-14
    @test rel([ForwardDiff.partials(d)[2] for d in yd], Ad * x .^ 2) < 1e-14
    # balancing in place keeps the storage
    s = 1 .+ rand(rng, 300)
    S = Float32Matrix(rand(rng, Float32, 300, 300))
    S0 = Float64.(S.data)
    Briareus._balance!(S, s)
    @test S.data == Float32.(s .* S0 .* s')
    @test Briareus._rowsums(S) ≈ vec(sum(Float64.(S.data); dims = 2)) rtol = 1e-14
end

@testset "block GMRES" begin
    rng = Xoshiro(3)
    n = 120
    M = I + 0.3 .* rand(rng, n, n) ./ n
    B = rand(rng, n, 5)
    B[:, 3] .= 0                                   # a zero column stays zero
    op!(Y, X) = mul!(Y, M, X)
    X = Briareus._gmres_block(op!, B, zero(B))
    @test norm(M * X - B) <= 1e-12 * norm(B) && iszero(X[:, 3])
    @test Briareus._gmres_block(op!, B, zero(B); maxiter = 1) === nothing
    # warm start at the solution: no iterations
    Briareus._reset_solve_stats!()
    @test Briareus._gmres_block(op!, B, X) == X && Briareus._SOLVE_STATS.iterations[] == 0
end

@testset "large radiation models (forced on)" begin
    wall = Material(:Au; albedo = WallLossAlbedo(), conversion = 0.8)
    t = hohlraum(; nz = 10, nphi = 16, capsule_ntheta = 8, capsule_nphi = 12, wall)
    q = [remake(b; aim = (0.0, 0.0, sign(b.port[3]) * 4.0)) for b in nif_quads()]
    op = DepositionOperator(t, q)
    dep = deposit(op, q)
    N = nelements(t)
    cap = elements(t, :capsule)
    rel(x, y) = maximum(abs.(x .- y)) / maximum(abs.(y))

    F64 = view_factors(t)
    F32 = view_factors(t; storage = Float32)
    @test F32 isa Float32Matrix && F32.data == Float32.(F64)
    m = RadiationModel(t)                                   # small: dense Float64, LU
    @test m.F isa Matrix{Float64} && !m.iterative
    mi = RadiationModel(t; solver = :iterative)             # Float64, iterative
    @test mi.F isa Matrix{Float64} && mi.iterative && mi.F == m.F
    m32 = RadiationModel(t; storage = Float32)              # Float32 ⇒ iterative
    @test m32.F isa Float32Matrix && m32.iterative
    @test sprint(show, m32) == "RadiationModel($N elements, Float32, iterative)"
    # the automatic threshold
    old = Briareus.LARGE_MODEL[]
    try
        Briareus.LARGE_MODEL[] = N - 1
        ma = RadiationModel(t)
        @test ma.F isa Float32Matrix && ma.iterative && ma.F == m32.F
        @test RadiationModel(t; storage = Float64).iterative          # Float64 above it
        @test !RadiationModel(t; storage = Float64, solver = :direct).iterative
    finally
        Briareus.LARGE_MODEL[] = old
    end
    # given matrices keep their storage unless told otherwise
    @test RadiationModel(t; F = m32.F, closed = false).F === m32.F
    @test RadiationModel(t; F = m32.F, storage = Float64, closed = false).F == Float64.(m32.F.data)
    @test RadiationModel(t; F = m.F, storage = Float32, closed = false).F == Float32Matrix(m.F)
    @test_throws ArgumentError RadiationModel(t; storage = Float16)
    @test_throws ArgumentError RadiationModel(t, m.F; solver = :lu)

    # balancing on Float32 storage: rows to the rounding of the stored factors
    @test maximum(abs, closure_error(m32.F)) < 1e-7
    G = t.elem_area .* m32.F
    @test maximum(abs.(G .- G')) < 1e-6 * maximum(G)

    # solves: temperature-dependent albedo iterated on each model
    s = solve(m, dep; time = 2.0)
    si = solve(mi, dep; time = 2.0)
    s32 = solve(m32, dep; time = 2.0)
    @test maximum(abs.(si.albedo .- s.albedo)) <= 1e-12
    @test rel(si.B, s.B) <= 1e-10 && rel(si.H, s.H) <= 1e-10
    @test si.B ≈ si.source .+ si.albedo .* si.H rtol = 1e-11    # the GMRES residual
    @test rel(s32.H, s.H) <= 1e-7
    a, a32 = legendre_moments(s.H, t, cap; lmax = 4), legendre_moments(s32.H, t, cap; lmax = 4)
    @test abs(a32[1] / a[1] - 1) < 1e-7
    @test abs(a32[3] / a32[1] - a[3] / a[1]) < 1e-7 && abs(a32[5] / a32[1] - a[5] / a[1]) < 1e-7
    @test isempty(mi.cache.keys) && isempty(m32.cache.keys)        # no factorizations
    eb, eb32 = energy_balance(s), energy_balance(s32)
    @test abs(eb32.lost) < 1e-6 * eb32.source && eb32.source == eb.source
    # warm start from a guess
    g = solve(mi, dep; time = 1.0)
    @test rel(solve(mi, 1.3 .* dep.element; time = 2.0, guess = g).H,
              solve(m, 1.3 .* dep.element; time = 2.0, guess = solve(m, dep; time = 1.0)).H) <= 1e-9

    # flux response, contracted and full
    α = albedos(t, 2.0)
    W = legendre_weights(t, cap; lmax = 4)
    R = flux_response(m, cap; time = 2.0, weights = W)
    @test rel(flux_response(mi, cap; time = 2.0, weights = W), R) <= 1e-10
    @test rel(flux_response(m32, cap; time = 2.0, weights = W), R) <= 1e-7
    @test rel(flux_response(mi, cap; albedo = α), flux_response(m, cap; albedo = α)) <= 1e-10

    # derivatives through the iterative solve (dual sources: value and partials as columns)
    f(mod) = p -> sum(solve(mod, p .* dep.element; albedo = α).H[cap])
    p0 = ones(N)
    @test rel(ForwardDiff.gradient(f(mi), p0), ForwardDiff.gradient(f(m), p0)) <= 1e-9
    @test_throws ArgumentError Briareus._iterative_radiosity(mi, α, big.(dep.element), nothing)

    # drive history and a scenario on the large path
    times = range(1.0, 3.0; length = 4)
    h = drive_history(m, op, τ -> q, times)
    h32 = drive_history(m32, op, τ -> q, times)
    @test rel(h32.moments, h.moments) <= 1e-7 && rel(h32.Tr, h.Tr) <= 1e-8
    sc = Scenario(t, q; observe = :capsule, model = m)
    sc32 = Scenario(t, q; observe = :capsule, model = m32)
    @test sc32.model === m32
    @test_throws ArgumentError Scenario(hohlraum(nz = 6, nphi = 8), q; observe = :capsule, model = m)
    x0 = initial_parameters(sc.param)
    obj = symmetry_objective(; modes = (2 => 1.0, 4 => 1.0))
    gr = ForwardDiff.gradient(x -> obj(evaluate(sc, x)), x0)
    @test rel(ForwardDiff.gradient(x -> obj(evaluate(sc32, x)), x0), gr) <= 1e-6
    e, e32 = evaluate(sc, x0), evaluate(sc32, x0)
    @test rel(e32.flux, e.flux) <= 1e-7                           # one solve, not the response
    @test getfield(sc32, :_response).R === nothing
    @test rel(flux_moments(e32), flux_moments(e)) <= 1e-7

    # multigroup radiosity
    sp = spectral_solve(s; edges = [0, 1000, 2000, Inf])
    sp32 = spectral_solve(s32; edges = [0, 1000, 2000, Inf])
    @test rel(sp32.H, sp.H) <= 1e-7

    # keyframe interpolation keeps the storage
    Fi = Briareus._interpolate_exchange(m32.F, m32.F, t.elem_area, t.elem_area, t.elem_area, 0.25)
    @test Fi isa Float32Matrix
    @test Fi ≈ m32.F rtol = 1e-7
end
