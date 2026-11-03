# An objective that differentiates internally: its derivative nests inside the
# optimizer's gradient (`_PointingTag`), whatever the two tags' ForwardDiff counts.

using Briareus, ForwardDiff, Test
using Briareus: _PointingTag, _pointing_gradient!

# a tag with the smallest count, as a tag made at run time gets when the optimizer's
# tag was compiled (and counted) in the precompile session
struct _EarlyTag end
ForwardDiff.tagcount(::Type{ForwardDiff.Tag{_EarlyTag,V}}) where {V} = UInt(0)

@testset "nested derivatives inside the optimizer's gradient" begin
    t = hohlraum(nz = 8, nphi = 12, capsule_ntheta = 6, capsule_nphi = 10)
    p = GroupPointing(nif_quads(power = 8.0); fields = (:z, :power))
    s = Scenario(t, p; observe = :capsule)
    Q(ev) = sum(ev.deposition.element) / 100 + 1e-3 * sum(abs2, ev.flux)
    # d/du sin(Q u) at u = 1 is Q cos Q, by ForwardDiff inside the objective, with its
    # own tag (`derivative`) or one that counts before the optimizer's
    inner(ev) = ForwardDiff.derivative(u -> sin(Q(ev) * u), 1.0)
    function early(ev)
        g(v) = sin(v[1] * Q(ev))              # (the inner dual first: either order nests)
        T = ForwardDiff.Tag{_EarlyTag,Float64}
        cfg = ForwardDiff.GradientConfig(g, [1.0], ForwardDiff.Chunk(1), T())
        return ForwardDiff.gradient(g, [1.0], cfg, Val(false))[1]
    end
    x = initial_parameters(p) .+ 0.01 .* sin.(1:nparams(p))
    function grad(obj)
        J(y) = obj(evaluate(s, y))
        cfg = ForwardDiff.GradientConfig(J, x, ForwardDiff.Chunk(4), ForwardDiff.Tag(_PointingTag(), Float64))
        out = ForwardDiff.DiffResults.GradientResult(similar(x))
        _pointing_gradient!(out, J, x, cfg)
        return copy(ForwardDiff.DiffResults.gradient(out)), ForwardDiff.DiffResults.value(out)
    end
    gq, q = grad(Q)
    @test any(!iszero, gq)
    # the reference: the chain rule on the gradient of Q itself, (cos Q − Q sin Q) ∇Q
    ref = (cos(q) - q * sin(q)) .* gq
    rel(g) = maximum(abs.(g .- ref)) / maximum(abs.(ref))
    gi, vi = grad(inner)
    @test isapprox(vi, q * cos(q); rtol = 1e-12) && rel(gi) <= 1e-12
    ge, ve = grad(early)
    @test isapprox(ve, q * cos(q); rtol = 1e-12) && rel(ge) <= 1e-12
    # the order itself
    PT = ForwardDiff.Tag{_PointingTag,Float64}
    @test ForwardDiff.:≺(PT, ForwardDiff.Tag{_EarlyTag,Float64}) &&
          !ForwardDiff.:≺(ForwardDiff.Tag{_EarlyTag,Float64}, PT) && !ForwardDiff.:≺(PT, PT)
end
