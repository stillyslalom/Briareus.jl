# Concurrent gradient chunks (`_pointing_gradient!`) reproduce ForwardDiff's serial
# chunk-mode gradient bit for bit, with CBET and a keep-out on the path.

using Briareus, ForwardDiff, Test
using Briareus: _PointingTag, _pointing_gradient!

struct _ChunkStop <: Exception end

@testset "parallel gradient chunks" begin
    t = hohlraum(length = 10.1, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9,
                 h = 1.2, capsule_h = 0.45)
    inner = (:inner23, :inner30)
    q = [remake(b; dlambda = b.group in inner ? 1.0 : 0.0) for b in nif_quads(power = 8.0)]
    p = GroupPointing(q; by = :quad, fields = (:z, :power))         # 96 parameters, 8 chunks
    s = Scenario(t, q; observe = :capsule, parameterization = p, lpi = CBET(leh_z = 5.05),
                 viewfactor_kw = (; closed = true))
    obj = asymmetry_objective(; keepout = ApertureKeepOut(t))
    J(x) = obj(evaluate(s, x))
    x = initial_parameters(p) .+ 0.01 .* sin.(1:nparams(p))
    cfg = ForwardDiff.GradientConfig(J, x, ForwardDiff.Chunk(12),
                                     ForwardDiff.Tag(_PointingTag(), Float64))
    ref = ForwardDiff.DiffResults.GradientResult(similar(x))
    ForwardDiff.gradient!(ref, J, x, cfg, Val(false))
    out = ForwardDiff.DiffResults.GradientResult(similar(x))
    _pointing_gradient!(out, J, x, cfg)
    @test ForwardDiff.DiffResults.gradient(out) == ForwardDiff.DiffResults.gradient(ref)
    @test ForwardDiff.DiffResults.value(out) == ForwardDiff.DiffResults.value(ref) == J(x)
    @test any(!iszero, ForwardDiff.DiffResults.gradient(out))
    # a ragged last chunk (n not a multiple of the chunk size)
    cfg7 = ForwardDiff.GradientConfig(J, x, ForwardDiff.Chunk(7),
                                      ForwardDiff.Tag(_PointingTag(), Float64))
    ForwardDiff.gradient!(ref, J, x, cfg7, Val(false))
    _pointing_gradient!(out, J, x, cfg7)
    @test ForwardDiff.DiffResults.gradient(out) == ForwardDiff.DiffResults.gradient(ref)
    # a chunk's exception reaches the caller unwrapped
    @test_throws _ChunkStop _pointing_gradient!(out, y -> throw(_ChunkStop()), x, cfg)
    # optimize_pointing takes the same path whether or not chunks run concurrently
    lo, hi = bounds(p, x; z = 0.3, power = 0.5)
    runs = map((true, false)) do par
        Briareus._PARALLEL_CHUNKS[] = par
        try
            optimize_pointing(s, obj; x0 = x, lower = lo, upper = hi, iterations = 2)
        finally
            Briareus._PARALLEL_CHUNKS[] = true
        end
    end
    @test runs[1].x == runs[2].x && runs[1].value == runs[2].value
    @test runs[1].evaluations == runs[2].evaluations
end
