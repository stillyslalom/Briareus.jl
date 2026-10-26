# Inference and compilation guards for the hot paths (CLAUDE.md, "Optimization path").
#
# 1. Functions that take user callables are compiled once for every callable: a second call
#    with new closure types adds no method instance to the body (`Base.@nospecializeinfer`
#    with `@nospecialize`; `@nospecialize` alone still infers, and so instantiates, the
#    body per closure type). Counting instances is deterministic, unlike compile times;
#    JET checks dispatch and inference, not recompilation per closure.
# 2. The optimization path infers concrete types on Float64 and on the optimizer's duals
#    (`@inferred`), and JET finds no runtime dispatch in it (`report_opt`) and no possible
#    error (`report_call`) in it, the solves or the drive history.
# 3. The threaded deposition kernel allocates per task, not per point (function barrier).
# Add a case to 1 for each new function that takes closures, objectives or callbacks,
# and to 2 for each new function on the optimization path.
# The fixture's calls have the precompile workload's argument types (fast with it).
# JET is a test-only dependency (Project.toml [extras]); run standalone without it, its
# checks are skipped with a warning.

const _INF_JET = Base.find_package("JET") !== nothing
_INF_JET && @eval using JET

# Method instances of the body of keyword function `f` (its only method), inference-only
# ones included (`Base.bodyfunction` and `Base.specializations` are internal, present
# from Julia 1.10).
function _inf_body_instances(f)
    body = Base.bodyfunction(only(methods(f)))
    body === nothing && error("$f has no keyword body")
    return sum(m -> count(_ -> true, Base.specializations(m)), methods(body))
end

# Whether JET's analysis `report` (`JET.report_opt` or `JET.report_call`) of `f(args...)`
# reports nothing in Briareus' own frames; prints the reports otherwise. Base, Optim and
# ForwardDiff internals are not ours to fix.
function _inf_jet_clean(report, f, args...)
    r = report(f, Base.typesof(args...); target_modules = (Briareus,))
    ok = isempty(JET.get_reports(r))
    ok || show(stdout, r)
    return ok
end

# Bytes allocated by the threaded deposition kernel (warmed), behind a function barrier.
function _inf_kernel_allocs(op, beams, split = nothing)
    T = typeof(first(beams).power)
    frames = [Briareus.BeamFrame(b) for b in beams]
    screens = [Briareus._screen_frame(f, b) for (f, b) in zip(frames, beams)]
    slopes = [Briareus._cone_slopes(op.dirs[bi], frames[bi]) for bi in eachindex(beams)]
    nblock, nback = size(op.plan.boff, 1) - 1, length(op.plan.back)
    E = zeros(T, nelements(op.target))
    PB, LB, SB = zeros(T, length(beams), nblock), zeros(T, length(op.target.components), nback), zeros(T, nback)
    run() = Briareus._deposit_units!(E, PB, LB, SB, op, frames, screens, slopes, nothing, split)
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
    objectives = (symmetry_objective(; drive = 1.0, keepout), asymmetry_objective(; keepout))
    tw = hohlraum(; geometry..., wall = Material(:Au; albedo = WallLossAlbedo(), conversion = 0.8))
    mw = RadiationModel(tw; F = s.model.F, closed = true)
    pulse = PulseShape([0.0, 1.0, 2.0], [0.0, 4.0, 4.0])
    history(beams_at, times; kw...) = drive_history(mw, s.deposition, beams_at, times; kw...)
    # spherical pointing fields on the capsule (absolute, and relative with split quads)
    ps = GroupPointing(quads; fields = (:theta_s, :dphi_s, :defocus, :power), radius = 0.9)
    xs0 = initial_parameters(ps)
    psr = GroupPointing(nif_beams(power = 2.0)[1:16]; by = :quad, relative = true, radius = 0.9,
                        fields = (:theta_s, :dphi_s, :split_theta, :split_phi))
    # a corona about the capsule (the quads' spots reach it: live rays on both paths)
    corona = CoronaOperator(CoronaProfile(; critical_radius = 0.9, scale_length = 0.2, absorption = 0.9),
                            quads; ncone = 1, lmax = 4, spacing = 0.2, panels = 20)

    @testset "closures compile once" begin
        history(τ -> with_pulses(quads, pulse, τ), [0.5, 1.0]; record = sol -> sol)
        n = _inf_body_instances(Briareus._drive_history)
        # new `beams_at`, `record` and `times` types
        history(τ -> with_pulses(quads, pulse, 0.9τ), 0.5:0.5:1.0; record = sol -> sol.H[1])
        @test _inf_body_instances(Briareus._drive_history) == n

        optimize_pointing(s, objectives[1]; x0, iterations = 1, evaluations = 2)
        n = _inf_body_instances(optimize_pointing)
        # a new objective type and a callback
        optimize_pointing(s, asymmetry_objective(; drive = 1.0); x0, iterations = 1,
                          evaluations = 2, callback = st -> false)
        @test _inf_body_instances(optimize_pointing) == n
    end

    @testset "optimization path infers" begin
        for x in (x0, _inf_duals(x0))   # Float64 and the optimizer's duals
            T = eltype(x)
            beams = @inferred apply(p, x)
            @test beams isa Vector{Beam{T}}
            @test (@inferred deposit(s.deposition, beams)).element isa Vector{T}
            @test (@inferred deposit(s.deposition, beams; absorption = T(0.9))).element isa Vector{T}
            @test (@inferred deposit(corona, beams)).ylm isa Vector{T}
            ev = evaluate(s, x)   # `Scenario.lpi` is dynamic by design: one dispatch per call
            @test (@inferred flux_moments(ev)) isa Vector{T}
            @test (@inferred objectives[1](ev)) isa T
            @test (@inferred objectives[2](ev)) isa T
            c = flux_ylm(ev)
            @test (@inferred mode_rms(c, 6)) isa T
            @test (@inferred mode_rms(c, 6; m = :nonzero)) isa T
        end
        # pointing on the capsule sphere (polar direct drive), with split quads
        for x in (xs0, _inf_duals(xs0))
            T = eltype(x)
            beams = @inferred apply(ps, x)
            @test beams isa Vector{Beam{T}}
            @test (@inferred deposit(s.deposition, beams)).element isa Vector{T}
            @test (@inferred apply(psr, zeros(T, nparams(psr)))) isa Vector{Beam{T}}
        end
        @test (@inferred flux_response(s.model, s.observed)) isa Matrix{Float64}
        # a dual absorption on Float64 beams promotes the deposition
        ηd = only(_inf_duals([0.9]))
        @test (@inferred deposit(s.deposition, quads; absorption = ηd)).element isa Vector{typeof(ηd)}
    end

    if _INF_JET
        @testset "JET" begin
            dep = deposit(s.deposition, quads)
            for x in (x0, _inf_duals(x0))
                beams, ev = apply(p, x), evaluate(s, x)
                # no runtime dispatch in the per-evaluation kernels (the objectives
                # dispatch by design on `Evaluation.scenario` in the uncached moments,
                # as the solves do on albedo models and the model's storage). Julia 1.10's
                # inference leaves per-beam setup dynamic (a comprehension capturing `T`,
                # `Union{Nothing,…}` profile fields) that later versions resolve, so the
                # deposition and keep-out checks run from 1.11.
                @test _inf_jet_clean(JET.report_opt, apply, p, x)
                if VERSION >= v"1.11"
                    @test _inf_jet_clean(JET.report_opt, deposit, s.deposition, beams)
                    @test _inf_jet_clean(JET.report_opt, (o, b, η) -> deposit(o, b; absorption = η),
                                         s.deposition, beams, eltype(x)(0.9))
                    @test _inf_jet_clean(JET.report_opt, Briareus._keepout_penalty, keepout, beams)
                    @test _inf_jet_clean(JET.report_opt, deposit, corona, beams)
                end
                # no possible errors
                @test _inf_jet_clean(JET.report_call, (o, b, η) -> deposit(o, b; absorption = η),
                                     s.deposition, beams, eltype(x)(0.9))
                @test _inf_jet_clean(JET.report_call, deposit, corona, beams)
                @test _inf_jet_clean(JET.report_call, evaluate, s, x)
                @test _inf_jet_clean(JET.report_call, flux_moments, ev)
                @test _inf_jet_clean(JET.report_call, objectives[1], ev)
                @test _inf_jet_clean(JET.report_call, objectives[2], ev)
                @test _inf_jet_clean(JET.report_call, mode_rms, flux_ylm(ev), 6)
            end
            # spherical pointing fields: no dispatch (from 1.11, as the deposition checks
            # above), no possible errors
            for x in (xs0, _inf_duals(xs0))
                if VERSION >= v"1.11"
                    @test _inf_jet_clean(JET.report_opt, apply, ps, x)
                    @test _inf_jet_clean(JET.report_opt, deposit, s.deposition, apply(ps, x))
                    @test _inf_jet_clean(JET.report_opt, apply, psr, zeros(eltype(x), nparams(psr)))
                end
                @test _inf_jet_clean(JET.report_call, apply, ps, x)
                @test _inf_jet_clean(JET.report_call, apply, psr, zeros(eltype(x), nparams(psr)))
            end
            @test _inf_jet_clean(JET.report_call, flux_response, s.model, s.observed)
            @test _inf_jet_clean(JET.report_call, solve, s.model, dep)
            @test _inf_jet_clean(JET.report_call, (m, d) -> solve(m, d; time = 1.0), mw, dep)
            @test _inf_jet_clean(JET.report_call, history,
                                 τ -> with_pulses(quads, pulse, τ), [0.5, 1.0])
            # the history's deposition with the M-band's per-group split (from 1.11, as
            # the deposition checks above)
            if VERSION >= v"1.11"
                @test _inf_jet_clean(JET.report_opt, Briareus._deposit, s.deposition, quads,
                                     Briareus._mband_split(nif_mband(), t, quads))
            end
            @test _inf_jet_clean(JET.report_call, (s, o, x) -> optimize_pointing(s, o; x0 = x),
                                 s, objectives[1], x0)
        end
    else
        @warn "JET is not in the environment: test/inference.jl skips its JET checks"
    end

    @testset "threaded kernel allocates per task" begin
        # a boxed capture allocates per point (1.8 million here); the tasks' overhead
        # is ~1 kB each
        bound = 4096 * (Threads.nthreads() + 1)
        @test _inf_kernel_allocs(s.deposition, apply(p, x0)) < bound
        @test _inf_kernel_allocs(s.deposition, apply(p, _inf_duals(x0))) < bound
        # with the per-group (M-band) split of `drive_history`
        sp = Briareus._mband_split(nif_mband(), t, quads)
        @test _inf_kernel_allocs(s.deposition, quads, sp) < bound
    end
end
