# Inference and compilation guards for the hot paths (CLAUDE.md, "Optimization path").
#
# 1. Functions that take user callables are compiled once for every callable: a second call
#    with new closure types adds no method instance to the body (`Base.@nospecializeinfer`
#    with `@nospecialize`; `@nospecialize` alone still infers, and so instantiates, the
#    body per closure type). Counting instances is deterministic, unlike compile times.
# 2. The optimization path infers concrete types on Float64 and on the optimizer's duals.
# 3. The threaded deposition kernel allocates per task, not per point (function barrier).
# Add a case to 1 for each new function that takes closures, objectives or callbacks.
# The fixture's calls have the precompile workload's argument types (fast with it).

# Method instances of the body of keyword function `f` (its only method), inference-only
# ones included (`Base.bodyfunction` and `Base.specializations` are internal, present
# from Julia 1.10).
function _inf_body_instances(f)
    body = Base.bodyfunction(only(methods(f)))
    body === nothing && error("$f has no keyword body")
    return sum(m -> count(_ -> true, Base.specializations(m)), methods(body))
end

# Bytes allocated by the threaded deposition kernel (warmed), behind a function barrier.
function _inf_kernel_allocs(op, beams)
    T = typeof(first(beams).power)
    frames = [Briareus.BeamFrame(b) for b in beams]
    screens = [Briareus._screen_frame(f, b) for (f, b) in zip(frames, beams)]
    slopes = [Briareus._cone_slopes(op.dirs[bi], frames[bi]) for bi in eachindex(beams)]
    nblock, nback = size(op.plan.boff, 1) - 1, length(op.plan.back)
    E = zeros(T, nelements(op.target))
    PB, LB, SB = zeros(T, length(beams), nblock), zeros(T, length(op.target.components), nback), zeros(T, nback)
    run() = Briareus._deposit_units!(E, PB, LB, SB, op, frames, screens, slopes)
    run()
    return minimum(@allocated(run()) for _ in 1:3)
end

# The optimizer's duals at `x` (its tag and chunk size; see `optimize_pointing`).
function _inf_duals(x)
    N = min(length(x), 12)
    D = ForwardDiff.Dual{ForwardDiff.Tag{Briareus._PointingTag,Float64}}
    return D.(x, Ref(ForwardDiff.Partials(ntuple(i -> Float64(i == 1), N))))
end

@testset "inference guards" begin
    geometry = (; length = 10.1, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9,
                h = 2.0, capsule_h = 0.9)
    t = hohlraum(; geometry...)
    quads = nif_quads(power = 8.0)
    p = GroupPointing(quads; fields = (:z, :r, :power))
    x0 = initial_parameters(p)
    s = Scenario(t, quads; observe = :capsule, parameterization = p,
                 viewfactor_kw = (; closed = true))
    keepout = (ApertureKeepOut(t), SphereKeepOut(t))

    @testset "closures compile once" begin
        tw = hohlraum(; geometry..., wall = Material(:Au; albedo = WallLossAlbedo(), conversion = 0.8))
        mw = RadiationModel(tw; F = s.model.F, closed = true)
        pulse = PulseShape([0.0, 1.0, 2.0], [0.0, 4.0, 4.0])
        drive_history(mw, s.deposition, τ -> with_pulses(quads, pulse, τ), [0.5, 1.0];
                      record = sol -> sol)
        n = _inf_body_instances(Briareus._drive_history)
        # new `beams_at`, `record` and `times` types
        drive_history(mw, s.deposition, τ -> with_pulses(quads, pulse, 0.9τ), 0.5:0.5:1.0;
                      record = sol -> sol.H[1])
        @test _inf_body_instances(Briareus._drive_history) == n

        optimize_pointing(s, symmetry_objective(; keepout); x0, iterations = 1, evaluations = 2)
        n = _inf_body_instances(optimize_pointing)
        # a new objective type and a callback
        optimize_pointing(s, asymmetry_objective(; drive = 1.0); x0, iterations = 1,
                          evaluations = 2, callback = st -> false)
        @test _inf_body_instances(optimize_pointing) == n
    end

    @testset "optimization path infers" begin
        objectives = (symmetry_objective(; drive = 1.0, keepout), asymmetry_objective(; keepout))
        for x in (x0, _inf_duals(x0))   # Float64 and the optimizer's duals
            T = eltype(x)
            beams = @inferred apply(p, x)
            @test beams isa Vector{Beam{T}}
            @test (@inferred deposit(s.deposition, beams)).element isa Vector{T}
            ev = evaluate(s, x)   # `Scenario.lpi` is dynamic by design: one dispatch per call
            @test (@inferred flux_moments(ev)) isa Vector{T}
            @test (@inferred objectives[1](ev)) isa T
            @test (@inferred objectives[2](ev)) isa T
        end
        @test (@inferred flux_response(s.model, s.observed)) isa Matrix{Float64}
    end

    @testset "threaded kernel allocates per task" begin
        # a boxed capture allocates per point (1.8 million here); the tasks'
        # overhead is ~1 kB each
        bound = 4096 * (Threads.nthreads() + 1)
        @test _inf_kernel_allocs(s.deposition, apply(p, x0)) < bound
        @test _inf_kernel_allocs(s.deposition, apply(p, _inf_duals(x0))) < bound
    end
end
