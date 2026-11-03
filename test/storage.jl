# Float32 storage against the Float64 matrix it rounds: each entry within 2⁻²⁴ relative,
# so every product within 2⁻²⁴ of the Float64 product's absolute sum; and every product
# path (vector, fused block, BLAS panels, transposed, from the left, α and β, aliased
# operands) equal to the Float64 product with the widened entries.

@testset "Float32 storage against Float64" begin
    rng = Xoshiro(5)
    m, n = 130, 90
    A64 = randn(rng, m, n)
    A = Float32Matrix(A64)
    W = Float64.(A.data)                       # what the storage holds, widened
    @test all(abs.(W .- A64) .<= 2.0^-24 .* abs.(A64))
    @test Float32Matrix(A) === A && A == Float32Matrix(copy(A64))
    rel(x, y) = maximum(abs.(x .- y)) / maximum(abs.(y))

    # against the Float64 original: |F32 x − A x| ≤ 2⁻²⁴ |A| |x| (+ Float64 rounding)
    x = randn(rng, n)
    bound = 2.0^-24 .* (abs.(A64) * abs.(x)) .+ 1e-13 * maximum(abs.(A64) * abs.(x))
    @test all(abs.(A * x .- A64 * x) .<= bound)

    # every path against the widened Float64 product
    z = randn(rng, m)
    X1, Xf, Xp = randn(rng, n, 1), randn(rng, n, 9), randn(rng, n, 60)       # vector, fused, panels
    Zf, Zp = randn(rng, m, 9), randn(rng, m, 60)
    Lp = randn(rng, 7, m)
    @test rel(A' * z, W' * z) < 1e-14 && rel(transpose(A) * z, W' * z) < 1e-14
    @test all(rel(A * X, W * X) < 1e-14 for X in (X1, Xf, Xp))
    @test all(rel(A' * Z, W' * Z) < 1e-14 for Z in (z[:, :], Zf, Zp))
    @test rel(Lp * A, Lp * W) < 1e-14
    # five-argument products: α op(A) X + β Y
    Y1, Y2, Y3 = randn(rng, m, 60), randn(rng, n, 9), randn(rng, 7, n)
    R1, R2, R3 = 2 .* (W * Xp) .- 0.5 .* Y1, 2 .* (W' * Zf) .- 0.5 .* Y2, 2 .* (Lp * W) .- 0.5 .* Y3
    mul!(Y1, A, Xp, 2.0, -0.5); mul!(Y2, A', Zf, 2.0, -0.5); mul!(Y3, Lp, A, 2.0, -0.5)
    @test rel(Y1, R1) < 1e-14 && rel(Y2, R2) < 1e-14 && rel(Y3, R3) < 1e-14
    y = randn(rng, n); y0 = copy(y)
    mul!(y, transpose(A), z, -1.0, 3.0)
    @test rel(y, 3 .* y0 .- W' * z) < 1e-14
    # aliased operands: a square matrix times the matrix the result overwrites
    S = Float32Matrix(randn(rng, Float32, 60, 60))
    Ws = Float64.(S.data)
    X = randn(rng, 60, 5); X0 = copy(X)
    mul!(X, S, X)
    @test rel(X, Ws * X0) < 1e-14
    v = randn(rng, 60); v0 = copy(v)
    mul!(v, S', v)
    @test rel(v, Ws' * v0) < 1e-14
    @test_throws DimensionMismatch mul!(zeros(m), A, zeros(n + 1))
    @test_throws DimensionMismatch mul!(zeros(m, 3), A', zeros(m, 3))
    # duals through the transpose: values and partials as two products
    zd = [ForwardDiff.Dual(z[i], 1.0, z[i]^2) for i in eachindex(z)]
    yd = A' * zd
    @test rel(ForwardDiff.value.(yd), W' * z) < 1e-14 &&
          rel([ForwardDiff.partials(d)[2] for d in yd], W' * z .^ 2) < 1e-14
    # row sums accumulate in Float64; dense balancing as the Float32 one
    @test Briareus._rowsums(A) ≈ vec(sum(W; dims = 2)) rtol = 1e-14
    s = 1 .+ rand(rng, 60)
    D = copy(Ws)
    @test Briareus._balance!(D, s) ≈ s .* Ws .* s' rtol = 1e-15
end
