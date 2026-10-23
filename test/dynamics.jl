@testset "moving targets on keyframes" begin
    w = WallMotion(speed = 0.1, start = 1.0)
    @test wall_radius(w, 2.875, 0.5) == 2.875
    @test wall_radius(w, 2.875, 5.0) ≈ 2.475

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
    hist = drive_history(Keyframes(build, collect(times), quads), beams_at, times)
    for (k, τ) in enumerate(times)
        tk = build(τ)
        sol = solve(RadiationModel(tk; closed = true), deposit(DepositionOperator(tk, quads), beams_at(τ)); time = τ)
        @test hist.moments[:, k] ≈ legendre_moments(sol.H, tk, elements(tk, :capsule); lmax = 6) rtol = 1e-10
    end

    # every keyframe must have the same elements
    grow = τ -> hohlraum(; radius = 2.875, kw..., h = τ < 1 ? 0.8 : 0.6)
    @test_throws ArgumentError model_at(Keyframes(grow, [0.0, 2.0], quads), 1.0)
end
