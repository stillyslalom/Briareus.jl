@testset "moving targets on keyframes" begin
    w = WallMotion(speed = 0.1, start = 1.0, saturation = Inf)
    @test wall_radius(w, 2.875, 0.5) == 2.875
    @test wall_radius(w, 2.875, 5.0) ≈ 2.475
    # the default wall stops: ballistic to the knee, C¹ there, then approaching saturation
    d = WallMotion()
    @test (d.speed, d.start, d.saturation, d.knee) == (0.185, 0.8, 0.45, 0.05)
    tk = d.start + (d.saturation - d.knee) / d.speed
    @test wall_displacement(d, tk) ≈ d.saturation - d.knee
    @test wall_displacement(d, tk + 1e-6) - wall_displacement(d, tk) ≈ d.speed * 1e-6 rtol = 1e-4
    @test wall_displacement(d, 9.0) < d.saturation && wall_displacement(d, 9.0) > d.saturation - 1e-6
    @test wall_displacement(d, 0.5) == 0
    @test_throws ArgumentError WallMotion(speed = 0.1, saturation = 0.04)
    # the fitted defaults reproduce Young 2020's dynamic view-factor radius (N170827)
    y = young2020_radii()
    i = findall(t -> 0.5 <= t <= 6.75, y.t)
    r = [y.wall[1] - y.wall[k] - wall_displacement(WallMotion(start = 0.15 + 0.8), y.t[k]) for k in i]
    @test maximum(abs, r) < 0.035
    # the HYBRID-E waist edge of PRE 109 Fig. 1a: 0.36–0.41 mm in at 6 ns
    f1a = pre109_fig1a_wall()
    waist = findall(z -> z <= 0.7, f1a.z)
    @test all(0.33 .< 3.2 .- f1a.n221204[waist] .< 0.42)
    # N221204's HYDRA capsule radius (Kritcher 2024 Fig. 1): flat to the first shock's
    # breakout, the 6 ns ice–gas interface (0.92 of 0.985 mm), monotone, minimum at bang
    k = kritcher2024_capsule()
    @test all(k.n221204[k.t .<= 3.5] .== 1) && all(diff(k.n221204) .<= 0)
    @test abs(k.n221204[findfirst(==(6.0), k.t)] - 0.92 / 0.985) < 0.03
    @test k.t[argmin(k.n221204)] == 9.5 && minimum(k.n221204) < 0.07
    @test all(k.n210808[k.t .<= 9.0] .<= k.n221204[k.t .<= 9.0] .+ 1e-12)   # thinner shell implodes first
    # N210808's measured inner-cone SBS per quad (Kemp 2026 Fig. 2): the paper's 2.8 and
    # 0.9 kJ per quad, the 23° peak at 5.3 ns before the 30°'s at 5.75, both gone by 7 ns
    s = kemp2026_sbs()
    dt = s.t[2] - s.t[1]
    @test abs(sum(s.inner23) * dt - 2.8) < 0.15 && abs(sum(s.inner30) * dt - 0.9) < 0.05
    @test abs(s.t[argmax(s.inner23)] - 5.3) < 0.03 && abs(s.t[argmax(s.inner30)] - 5.75) < 0.03
    @test all(s.inner23[s.t .>= 6.95] .== 0) && all(s.inner30[s.t .>= 7.0] .== 0) && all(>=(0), s.inner23)

    kw = (length = 10.1, leh_radius = 1.55, capsule_radius = 0.9, h = 0.8, capsule_h = 0.4,
          wall = Material(:Au; albedo = 0.9, conversion = 0.8), corner_radius = 0.3)
    build = moving_hohlraum(w; radius = 2.875, kw...)
    t0, t1 = build(0.0), build(6.0)
    @test nelements(t0) == nelements(t1) && t0.comp_elems == t1.comp_elems
    @test nelements(t0) == nelements(hohlraum(; radius = 2.875, kw...))
    wall0, wall1 = elements(t0, :wall), elements(t1, :wall)
    r(t, e) = hypot(t.elem_centroid[e][1], t.elem_centroid[e][2])
    @test r(t1, wall1[1]) / r(t0, wall0[1]) ≈ 2.375 / 2.875 rtol = 1e-12   # chord centroids scale with R
    leh0, leh1 = elements(t0, :leh_top), elements(t1, :leh_top)
    @test sum(t1.elem_area[leh1]) ≈ sum(t0.elem_area[leh0])                 # holes keep their size
    @test maximum(abs(t1.elem_centroid[e][3]) for e in elements(t1, :cap_top)) ≈
          maximum(abs(t0.elem_centroid[e][3]) for e in elements(t0, :cap_top))

    quads = nif_quads()[1:8]
    kf = Keyframes(build, [0.0, 4.0, 8.0], quads)
    @test length(kf.time) == 3
    # between keyframes: interpolated exchange factors keep reciprocity and closure
    m = model_at(kf, 3.0)
    @test m.target !== model_at(kf, 4.0).target
    @test model_at(kf, 3.0) === m                                           # kept for its time
    A = m.target.elem_area
    @test maximum(abs.(A .* m.F .- (A .* m.F)')) < 1e-12 * maximum(A .* m.F)
    @test maximum(abs.(vec(sum(m.F; dims = 2)) .- 1)) < 1e-9
    @test model_at(kf, -1.0) === model_at(kf, 0.0)                          # nearest keyframe outside
    # only `keep` keyframe models are held
    model_at(kf, 8.0)
    @test length(kf.models) == 2
    # deposition: the keyframes' operators at their times, one built at any other
    @test deposition_at(kf, 4.0) === deposition_at(kf, 4.0)
    @test deposit(kf, quads, 2.0).element ≈ deposit(DepositionOperator(build(2.0), quads), quads).element

    # a still wall reproduces the fixed target's history
    still = moving_hohlraum(WallMotion(speed = 0.0); radius = 2.875, kw...)
    tf = still(0.0)
    times = 0.5:0.5:2.0
    beams_at = τ -> [remake(b; power = 1.0 + τ) for b in quads]
    ref = drive_history(RadiationModel(tf; closed = true), DepositionOperator(tf, quads), beams_at, times)
    dyn = drive_history(Keyframes(still, [0.0, 1.2, 3.0], quads), beams_at, times)
    @test dyn.moments ≈ ref.moments rtol = 1e-10
    @test dyn.Tr ≈ ref.Tr rtol = 1e-10

    # a moving wall's history on keyframes at its own times is the per-time rebuild
    # (two times: each costs a target, its view factors and its deposition operator)
    tk2 = [first(times), last(times)]
    hist = drive_history(Keyframes(build, tk2, quads), beams_at, tk2)
    for (k, τ) in enumerate(tk2)
        tk = build(τ)
        sol = solve(RadiationModel(tk; closed = true), deposit(DepositionOperator(tk, quads), beams_at(τ)); time = τ)
        @test hist.moments[:, k] ≈ legendre_moments(sol.H, tk, elements(tk, :capsule); lmax = 6) rtol = 1e-10
    end

    # a shrinking capsule mesh: the initial radius's counts, the area as R², and none of
    # it without `capsule_shrink`
    for mesh in (:compact, :latlong)
        shrinking = moving_hohlraum(w; radius = 2.875, kw..., capsule_mesh = mesh,
                                    capsule_shrink = ([0.0, 5.0, 9.0], [1.0, 0.8, 0.5]))
        s0, s7 = shrinking(0.0), shrinking(7.0)
        c0, c7 = elements(s0, :capsule), elements(s7, :capsule)
        @test s0.comp_elems == s7.comp_elems
        @test nelements(s0) == nelements(hohlraum(; radius = 2.875, kw..., capsule_mesh = mesh))
        @test sum(s7.elem_area[c7]) / sum(s0.elem_area[c0]) ≈ 0.65^2 rtol = 1e-12
        @test maximum(norm(s7.elem_centroid[e]) for e in c7) < 0.65 * 0.9
        @test sum(s0.elem_area[c0]) == sum(t0.elem_area[elements(t0, :capsule)]) || mesh === :latlong
    end

    # every keyframe must have the same elements
    grow =τ -> hohlraum(; radius = 2.875, kw..., h = τ < 1 ? 0.8 : 0.6)
    @test_throws ArgumentError model_at(Keyframes(grow, [0.0, 2.0], quads), 1.0)
end
