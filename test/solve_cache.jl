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
    @test rel(sw.B, B) <= 1e-12 && rel(sw.H, H) <= 1e-12 && sw.albedo == α
    @test length(mw.cache.keys) == 1 && mw.cache.keys[1] == albedos(tw, 2.0)
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
