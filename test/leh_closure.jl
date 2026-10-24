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

@testset "computed lip closure (LipClosure)" begin
    m = LipClosure()
    @test m.speed == 0.022 && m.depth == 0.24 && m.hold_time == 0.4 && m.hold == (6.0, 17.0)
    @test_throws ArgumentError LipClosure(hold = (5.0, 5.0))
    @test_throws ArgumentError LipClosure(hold_time = 0.0)
    # free expansion at the sound-speed law: v = speed (T/100 eV)^½, no hold below hold[1]
    ts = collect(0.0:0.25:4.0)
    d = lip_depth(m, ts, fill(400.0, length(ts)), zeros(length(ts)))
    @test d ≈ 2 * 0.022 .* ts rtol = 1e-12
    # held: relaxes to `depth` with time constant hold_time, whatever the temperature
    d = lip_depth(m, ts, fill(300.0, length(ts)), fill(30.0, length(ts)); depth = 0.4)
    @test d ≈ 0.24 .+ 0.16 .* exp.(-ts ./ 0.4) rtol = 1e-12
    # the floor caps it
    @test maximum(lip_depth(m, ts, fill(400.0, length(ts)), zeros(length(ts)); floor = 0.1)) == 0.1

    # hohlraum(; leh_closure = :computed): lips in rings down to 0.6 r₀
    gold9 = Material(:Au; albedo = 0.9, conversion = 0.8)
    kw = (length = 10.1, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9, h = 0.5,
          capsule_h = 0.3, wall = gold9)
    t = hohlraum(; kw..., leh_closure = :computed)
    @test Briareus.element_material(t, first(elements(t, :leh_top_lip))).albedo isa LipClosure
    open_t = hohlraum(; kw...)
    hole(t, n) = sum(t.elem_area[elements(t, n)])
    @test hole(t, :leh_top) ≈ π * (0.6 * 1.55)^2 rtol = 0.02
    @test hole(t, :leh_top) + hole(t, :leh_top_lip) ≈ hole(open_t, :leh_top) rtol = 0.02
    @test_throws ArgumentError hohlraum(; kw..., leh_closure = :prescribed)
    # the prescribed lip keeps its single ring (and 0.84 r₀)
    tg = hohlraum(; kw..., leh_closure = GapClosure(closure = 5.0, after = 0.9))
    @test length(elements(tg, :leh_top_lip)) == length(elements(t, :leh_top_lip)) ÷ 6
    # rings: the open area at an edge radius r is the chord polygon's at r, exactly
    L = only(filter(l -> l.name === :leh_top_lip, Briareus._lip_geometry(t)))
    @test L.r0 ≈ 1.55 && L.rin ≈ 0.6 * 1.55
    for r in (1.55, 1.42, 1.31, 1.0, 0.6 * 1.55)
        s = [Briareus._ring_closed(L.ra[j], L.rb[j], r) for j in eachindex(L.elems)]
        open_area = hole(t, :leh_top) + sum((1 .- s) .* L.area)
        @test open_area / (hole(t, :leh_top) + sum(L.area)) ≈ (r / 1.55)^2 rtol = 1e-10
    end
    # a single solve treats the lip as open; given a closed fraction it re-emits s α
    rm = RadiationModel(t)
    sol = solve(rm, zeros(nelements(t)) .+ 1e-3; time = 6.0)
    @test all(iszero, sol.albedo[elements(t, :leh_top_lip)])
    st = fill(NaN, nelements(t)); st[elements(t, :leh_top_lip)] .= 0.25
    sol = solve(rm, zeros(nelements(t)) .+ 1e-3; time = 6.0, burnthrough = st)
    @test sol.albedo[elements(t, :leh_top_lip)] ≈ 0.25 .* sol.intact[elements(t, :leh_top_lip)]
    @test all(>(0.5), sol.intact[elements(t, :leh_top_lip)])   # the wall's albedo
    @test all(sol.transmission[elements(t, :leh_top_lip)] .≈ 0.75)

    # drive_history integrates it per hole: open at the first time, the closure the
    # standalone integration gives on the history's own lip temperatures and intensities
    quads = nif_quads()
    op = DepositionOperator(t, quads)
    P = PulseShape([0.0, 0.5, 2.0, 3.0, 6.0, 6.5, 9.0], [1.0, 1.0, 1.0, 9.0, 9.0, 1.0, 1.0])  # TW per quad
    times = collect(0.5:0.25:8.0)
    lipe = elements(t, :leh_top_lip)
    lipT(sol) = radiation_temperature(sum(sol.H[lipe] .* t.elem_area[lipe]) / sum(t.elem_area[lipe]))
    h = drive_history(rm, op, τ -> [remake(q; power = P(τ)) for q in quads], times; record = lipT)
    r = h.lip_radius[:leh_top_lip]
    @test r[1] ≈ 1.55
    @test h.lip_radius[:leh_bottom_lip] ≈ r rtol = 1e-6         # symmetric drive
    I = [24 * P(τ) / (π * 1.55^2) for τ in times]               # 24 quads per hole
    @test r ≈ 1.55 .- lip_depth(m, times, h.record, I) atol = 1e-4
    k = findfirst(≈(6.0), times)
    @test r[k] ≈ 1.55 - 0.24 atol = 0.005                       # held at `depth` at peak
    @test r[end] < r[k] - 0.015                                 # free again after the drop
    # the lip's escaped power falls as it closes, and the budget still closes
    esc = h.escaped[:leh_top_lip]
    @test esc[1] > 0
    b = energy_budget(h)
    total = sum(values(b.heating)) + sum(values(b.escaped)) + sum(values(b.transmitted); init = 0.0) +
            b.lost + b.unabsorbed + b.backscattered
    @test total ≈ b.setup rtol = 2e-3

    # calibration: Chen 2020 Fig. 27's measured edge on the five basecamp shots (3.37 mm
    # LEHs), with the lip temperature of the model's coarse basecamp history of N171010
    # (`gold(albedo = WallLossAlbedo(), conversion = nif_conversion())`, h 0.8 / 0.4)
    tT = collect(0.25:0.25:8.0)
    Tl = [88.0, 119.8, 123.0, 123.0, 124.5, 124.1, 124.5, 123.2, 122.5, 143.6, 163.1, 170.8,
          170.9, 179.0, 205.5, 226.1, 249.6, 268.7, 275.2, 279.4, 282.6, 286.9, 288.4, 275.5,
          225.6, 200.6, 183.5, 170.0, 158.4, 147.9, 138.1, 128.7]
    c = chen2020_leh()
    r0 = 1.685
    res, post = Float64[], Float64[]
    for s in unique(c.shot)
        kl = (c.shot .== s) .& (c.kind .== "laser")
        Pc = PulseShape(c.t[kl], c.value[kl])
        tt = [0.0; tT]
        dd = lip_depth(m, tt, [Tl[1]; Tl], [Pc(τ) / 2 / (π * r0^2) for τ in tt])
        for i in findall((c.shot .== s) .& (c.kind .== "gleh"))
            j = searchsortedlast(tt, c.t[i]); w = (c.t[i] - tt[j]) / (tt[j+1] - tt[j])
            rm_ = r0 - ((1 - w) * dd[j] + w * dd[j+1])
            push!(res, c.value[i] - rm_)
            c.t[i] > 6.1 && push!(post, c.value[i] - rm_)
        end
    end
    @test length(res) == 31
    @test sqrt(sum(abs2, res) / length(res)) < 0.035        # 0.033 mm (shot scatter ~0.03)
    @test abs(sum(res) / length(res)) < 0.01
    @test sqrt(sum(abs2, post) / length(post)) < 0.03       # after laser-off (LASNEX: 0.11)
end
