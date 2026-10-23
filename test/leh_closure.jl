@testset "closing entrance holes" begin
    gold9 = Material(:Au; albedo = 0.9, conversion = 0.8)
    kw = (length = 10.1, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9, h = 0.5,
          capsule_h = 0.3, wall = gold9)
    open_t = hohlraum(; kw...)
    lip = GapClosure(closure = 5.0, width = 0.5, after = 0.9)
    t = hohlraum(; kw..., leh_closure = lip)
    names = [c.name for c in t.components]
    @test :leh_top_lip in names && :leh_bottom_lip in names
    @test findfirst(==(:leh_top_lip), names) == findfirst(==(:leh_top), names) + 1
    hole(t, n) = sum(t.elem_area[elements(t, n)])
    @test hole(t, :leh_top) + hole(t, :leh_top_lip) ≈ hole(open_t, :leh_top) rtol = 0.02
    @test hole(t, :leh_top) ≈ π * (0.84 * 1.55)^2 rtol = 0.02
    @test_throws ArgumentError hohlraum(; kw..., leh_closure = lip, leh_open_radius = 1.6)

    # laser light passes the lip: the same deposition as through the open hole
    quads = nif_quads()[1:12]
    d_open = deposit(DepositionOperator(open_t, quads), quads)
    d = deposit(DepositionOperator(t, quads), quads)
    @test d.per_beam ≈ d_open.per_beam rtol = 1e-12
    # and the keep-out rim is the hole's full radius
    @test ApertureKeepOut(t).rims ≈ ApertureKeepOut(open_t).rims
    @test ApertureKeepOut(t).names == [:leh_top, :leh_bottom]

    m, m_open = RadiationModel(t), RadiationModel(open_t)
    cap(t) = elements(t, :capsule)
    # open before closure: the capsule sees what it sees through the open hole
    s0 = solve(m, d; time = 0.0)
    s0_open = solve(m_open, d_open; time = 0.0)
    flux(s, t) = sum(s.H[cap(t)] .* t.elem_area[cap(t)]) / sum(t.elem_area[cap(t)])
    @test flux(s0, t) ≈ flux(s0_open, open_t) rtol = 2e-3
    b0 = energy_balance(s0)
    @test b0.escaped[:leh_top_lip] > 0
    @test abs(b0.absorbed[:leh_top_lip]) < 1e-12

    # closed: the lip re-emits like the cap it extends
    s1 = solve(m, d; time = 6.0)
    @test all(s1.albedo[elements(t, :leh_top_lip)] .≈ 0.9)
    b1 = energy_balance(s1)
    @test b1.escaped[:leh_top_lip] ≈ 0 atol = 1e-12
    lipP = sum(s1.H[e] * t.elem_area[e] for e in elements(t, :leh_top_lip))
    @test b1.absorbed[:leh_top_lip] ≈ 0.1 * lipP rtol = 1e-10
    total = sum(values(b1.absorbed)) + sum(values(b1.escaped)) + sum(values(b1.transmitted)) + b1.lost
    @test total ≈ b1.source rtol = 1e-9
    # (the same laser power on the wall: the smaller hole would clip the beams)
    small = hohlraum(; kw..., leh_radius = 0.84 * 1.55)   # caps reaching the open radius
    w, ws = elements(t, :wall), elements(small, :wall)
    @test t.elem_area[w] == small.elem_area[ws]
    Pw, Ps = zeros(nelements(t)), zeros(nelements(small))
    Pw[w] .= d.element[w]
    Ps[ws] .= d.element[w]
    @test flux(solve(m, Pw; time = 6.0), t) ≈ flux(solve(RadiationModel(small), Ps), small) rtol = 5e-3

    # a drive history resolves the closure with `transition` steps
    hist = drive_history(m, DepositionOperator(t, quads), τ -> quads, [4.0, 6.0]; transition = 0.25)
    @test length(hist.time) > 2
    @test hist.escaped[:leh_top_lip][end] < 1e-12 < hist.escaped[:leh_top_lip][1]
end

@testset "closing entrance holes on a moving wall" begin
    lip = GapClosure(closure = 5.0, width = 0.5, after = 0.9)
    kw = (length = 10.1, leh_radius = 1.55, capsule_radius = 0.9, h = 0.8, capsule_h = 0.4,
          wall = Material(:Au; albedo = 0.9, conversion = 0.8), leh_closure = lip)
    build = moving_hohlraum(WallMotion(speed = 0.1, start = 1.0); radius = 2.875, kw...)
    t0, t1 = build(0.0), build(6.0)
    @test t0.comp_elems == t1.comp_elems && :leh_top_lip in [c.name for c in t1.components]
    lipA(t) = sum(t.elem_area[elements(t, :leh_top_lip)])
    @test lipA(t1) ≈ lipA(t0)
    quads = nif_quads()[1:8]
    hist = drive_history(Keyframes(build, [4.0, 6.0], quads), τ -> quads, [4.0, 6.0]; transition = 0.25)
    @test length(hist.time) > 2
    @test hist.escaped[:leh_top_lip][end] < 1e-12 < hist.escaped[:leh_top_lip][1]
end
