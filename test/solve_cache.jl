# Cached radiosity factorization and the clearance progress/cancel callback.

@testset "radiosity factorization cache" begin
    t = hohlraum(nz = 10, nphi = 16, capsule_ntheta = 8, capsule_nphi = 12,
                 wall = gold(albedo = 0.85, conversion = 0.8))
    q = [remake(b; aim = (0.0, 0.0, sign(b.port[3]) * 4.0)) for b in nif_quads()]
    dep = deposit(DepositionOperator(t, q), q)
    m = RadiationModel(t)
    cache = m.cache
    S = Briareus.conversions(t) .* dep.element ./ t.elem_area
    rel(x, y) = maximum(abs.(x .- y)) / maximum(abs.(y))
    # references: the uncached formulas
    α0 = albedos(t)
    α1 = albedos(t) .* 0.9
    α2 = albedos(t) .* 0.8
    Bref(α) = (I - α .* m.F) \ S
    refs = Dict(α0 => Bref(α0), α1 => Bref(α1), α2 => Bref(α2))

    @test isempty(cache.keys)
    s1 = solve(m, dep)
    @test length(cache.keys) == 1
    s2 = solve(m, dep)
    @test s1.B == s2.B && s1.H == s2.H
    @test rel(s1.B, refs[α0]) <= 1e-12
    @test rel(s1.H, m.F * refs[α0]) <= 1e-12
    @test length(cache.keys) == 1
    # explicit albedo: new key, most recently used first
    sx = solve(m, dep; albedo = α1)
    @test rel(sx.B, refs[α1]) <= 1e-12
    @test cache.keys == [α1, α0]
    @test solve(m, dep; albedo = α0).B == s1.B          # same as the material albedo
    @test cache.keys == [α0, α1]
    # eviction keeps FACTOR_CACHE_SIZE entries and results stay exact
    @test rel(solve(m, dep; albedo = α2).B, refs[α2]) <= 1e-12
    @test length(cache.keys) == Briareus.FACTOR_CACHE_SIZE
    @test !(α1 in cache.keys)
    @test rel(solve(m, dep; albedo = α1).B, refs[α1]) <= 1e-12

    # flux_response shares the cache and matches the direct solve
    obs = elements(t, :capsule)
    n = length(cache.keys)
    R = flux_response(m, obs; albedo = α1)
    @test length(cache.keys) == n && cache.keys[1] == α1
    @test rel(R * dep.element, m.F[obs, :] * refs[α1]) <= 1e-12
    η = Briareus.conversions(t)
    Rref = Matrix((((I - α1 .* m.F)' \ Matrix(m.F[obs, :]')) .* (η ./ t.elem_area))')
    @test rel(R, Rref) <= 1e-12

    Briareus.empty_cache!(m)
    @test isempty(cache.keys) && isempty(cache.facts)

    # concurrent solves on worker threads with more albedos than the cache holds
    αs = [α0, α1, α2]
    tasks = [Threads.@spawn solve(m, dep; albedo = αs[mod1(i, 3)]) for i in 1:24]
    sols = fetch.(tasks)
    @test all(rel(s.B, refs[αs[mod1(i, 3)]]) <= 1e-12 for (i, s) in enumerate(sols))
    @test length(cache.keys) == length(cache.facts) <= Briareus.FACTOR_CACHE_SIZE

    # temperature-dependent albedo: same iteration as without the cache
    tw = hohlraum(nz = 10, nphi = 16, capsule_ntheta = 8, capsule_nphi = 12,
                  wall = Material(:Au; albedo = WallLossAlbedo(), conversion = 0.8))
    mw = RadiationModel(tw)
    depw = deposit(DepositionOperator(tw, q), q)
    sw = solve(mw, depw; time = 2.0)
    sw2 = solve(mw, depw; time = 2.0)
    @test sw.B == sw2.B && sw.albedo == sw2.albedo
    Sw = Briareus.conversions(tw) .* depw.element ./ tw.elem_area
    α = albedos(tw, 2.0)
    B = H = nothing
    for _ in 1:50
        B = (I - α .* mw.F) \ Sw
        H = mw.F * B
        αn = albedos(tw, 2.0, H)
        δ = maximum(abs.(αn .- α))
        δ < 1e-4 && break
        α = αn
    end
    # (iterates after the first are solved by preconditioned GMRES to a 1e-12 residual)
    @test rel(sw.B, B) <= 1e-10 && rel(sw.H, H) <= 1e-10 && maximum(abs.(sw.albedo .- α)) <= 1e-12
    @test length(mw.cache.keys) == 1 && mw.cache.keys[1] == albedos(tw, 2.0)
end

@testset "preconditioned iterative radiosity solves" begin
    q = [remake(b; aim = (0.0, 0.0, sign(b.port[3]) * 4.0)) for b in nif_quads()]
    tw = hohlraum(nz = 12, nphi = 16, capsule_ntheta = 8, capsule_nphi = 12,
                  wall = Material(:Au; albedo = WallLossAlbedo(), conversion = 0.8))
    m = RadiationModel(tw)
    dep = deposit(DepositionOperator(tw, q), q)
    S = Briareus.conversions(tw) .* dep.element ./ tw.elem_area
    rel(x, y) = maximum(abs.(x .- y)) / maximum(abs.(y))
    direct(α, S) = (I - α .* m.F) \ S

    # GMRES with the LU of a nearby albedo against the direct solve
    α0 = albedos(tw, 2.0)
    α1 = albedos(tw, 2.0, m.F * direct(α0, S))
    P = lu(I - α0 .* m.F)
    x = Briareus._gmres_radiosity(m.F, α1, P, S, P \ S)
    @test rel(x, direct(α1, S)) <= 1e-10
    @test norm(S - (I - α1 .* m.F) * x) <= 1e-12 * norm(S)
    # a preconditioner far from the system does not converge in a few iterations
    @test Briareus._gmres_radiosity(m.F, α1, lu(Matrix(1.0I, size(m.F)...)), S, zero(S);
                                    maxiter = 2) === nothing
    # ... and _iterate_solve then factors afresh and keeps that as the preconditioner
    Briareus.empty_cache!(m)
    m.cache.precond[] = lu(rand(Xoshiro(1), size(m.F)...))
    Briareus._reset_solve_stats!()
    @test rel(Briareus._iterate_solve(m, α1, S, nothing), direct(α1, S)) <= 1e-12
    @test Briareus._SOLVE_STATS.lu[] == 1
    @test m.cache.precond[].factors ≈ lu(I - α1 .* m.F).factors

    # temperature-dependent histories: the iterative path against a fresh LU per iterate
    function both(f)
        Briareus.empty_cache!(m)
        Briareus._ITERATIVE_SOLVE[] = false
        a = f()
        Briareus.empty_cache!(m)
        Briareus._ITERATIVE_SOLVE[] = true
        b = f()
        return a, b
    end
    s_old, s_new = both(() -> solve(m, dep; time = 2.0))
    @test maximum(abs.(s_old.albedo .- s_new.albedo)) <= 1e-12
    @test rel(s_new.B, s_old.B) <= 1e-10 && rel(s_new.H, s_old.H) <= 1e-10
    # warm starts (the previous solution as guess and source of the preconditioner)
    g = solve(m, dep; time = 1.0)
    w_old, w_new = both(() -> solve(m, 1.3 .* dep.element; time = 2.0, guess = g))
    @test maximum(abs.(w_old.albedo .- w_new.albedo)) <= 1e-12
    @test rel(w_new.H, w_old.H) <= 1e-10
    # and from a model with the same elements (a moving target's previous time)
    m2 = RadiationModel(tw; F = 0.999 .* m.F)
    Briareus.empty_cache!(m)
    g = solve(m, dep; time = 2.0)
    Briareus._reset_solve_stats!()
    s2 = solve(m2, dep; time = 2.0, guess = g)
    @test Briareus._SOLVE_STATS.lu[] == 0 && m2.cache.precond[] === m.cache.precond[]
    Briareus._ITERATIVE_SOLVE[] = false
    s2_old = solve(m2, dep; time = 2.0, guess = g)
    Briareus._ITERATIVE_SOLVE[] = true
    @test maximum(abs.(s2.albedo .- s2_old.albedo)) <= 1e-12 && rel(s2.H, s2_old.H) <= 1e-10

    # dual-number sources: value and partials solved as columns
    dP = ForwardDiff.Dual{Nothing}.(dep.element, 0.5 .* dep.element, dep.element .^ 2 ./ maximum(dep.element))
    d_old, d_new = both(() -> solve(m, dP; time = 2.0, guess = g))
    @test maximum(abs.(d_old.albedo .- d_new.albedo)) <= 1e-12
    @test rel(ForwardDiff.value.(d_new.H), ForwardDiff.value.(d_old.H)) <= 1e-10
    for k in 1:2
        @test rel(ForwardDiff.partials.(d_new.H, k), ForwardDiff.partials.(d_old.H, k)) <= 1e-10
    end

    # a drive history (warm-started steps, step halving on flux_change)
    op = DepositionOperator(tw, q)
    beams_at = τ -> [remake(b; power = 1.0 + 4τ^2) for b in q]
    h_old, h_new = both(() -> drive_history(m, op, beams_at, 0.25:0.25:2.0; flux_change = 0.02))
    @test h_old.time == h_new.time
    @test maximum(abs.(h_old.Tr .- h_new.Tr)) <= 1e-8
    @test h_new.moments ≈ h_old.moments rtol = 1e-9
end

@testset "clearance progress and cancel" begin
    th = hohlraum(length = 10.1, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9)
    qs = nif_quads()
    qs = [qs[findfirst(b -> b.group === g, qs)] for g in (:inner23, :inner30, :outer44, :outer50)]
    rz = Dict(:inner23 => (1.09, 3.09), :inner30 => (0.6, 3.63), :outer44 => (0.0, 5.2),
              :outer50 => (0.0, 5.01))
    qs = map(qs) do b
        r0, z0 = rz[b.group]
        φ = port_angles(b)[2]
        remake(b; aim = (-r0 * cosd(φ), -r0 * sind(φ), sign(b.port[3]) * z0))
    end
    kw = (; targets = (:wall, :capsule), hardware = (:cap_top, :cap_bottom),
          peak_fraction = 0.01, nangles = 4, nspot = 31)
    ref = clearance(th, qs; kw...)
    @test !ref.canceled && length(ref.rows) == 8
    for batch in (1, 3, 8)
        seen = NamedTuple[]
        tasks = Task[]
        r = clearance(th, qs; kw..., batch,
                      callback = p -> (push!(seen, p); push!(tasks, current_task()); false))
        @test r.rows == ref.rows && !r.canceled
        @test [p.done for p in seen] == [0; collect(batch:batch:3); 4]
        @test seen[end] == (; done = 4, total = 4)
        @test all(==(current_task()), tasks)          # funneled to the calling task
    end
    # cancel after the second beam: rows of the first two beams only
    r = clearance(th, qs; kw..., batch = 1, callback = p -> p.done >= 2)
    @test r.canceled && r.rows == ref.rows[1:4]
    @test occursin("canceled", sprint(show, MIME"text/plain"(), r))
    # cancel before any work, and a cancel after the last batch is ignored
    r = clearance(th, qs; kw..., callback = _ -> true)
    @test r.canceled && isempty(r.rows)
    r = clearance(th, qs; kw..., batch = 2, callback = p -> p.done == 4)
    @test !r.canceled && r.rows == ref.rows
    # the positional constructor (no cancel flag) still works
    r7 = ClearanceReport(ref.rows, ref.contour, ref.transverse, ref.errors, ref.margin,
                         ref.targets, ref.thresholds)
    @test !r7.canceled
end
