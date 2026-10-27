@testset "inner-beam bubble absorption" begin
    pointing = Dict(:outer50 => (0.0, 5.01), :outer44 => (0.0, 5.20),
                    :inner30 => (0.60, 3.63), :inner23 => (1.09, 3.09))
    function aimed(b)
        r, z = pointing[b.group]
        φ = port_angles(b)[2]
        return remake(b; aim = (-r * cosd(φ), -r * sind(φ), sign(b.port[3]) * z))
    end
    quads = aimed.(nif_quads())
    inner(bs) = [b.group in (:inner23, :inner30) for b in bs]

    # calibration: the tip against Izumi 2021 Fig. 8a (seven shots, ±40 µm), the absorbed
    # fraction against its HYDRA run (Fig. 10a) with the bubble at Izumi's z = ±3 mm
    m3 = BubbleAbsorption(radius = 2.88, turn_on = 0.3)
    m15 = BubbleAbsorption(radius = 2.88, turn_on = 0.3, fill = 0.15)
    tr = izumi2021_bubble()
    meas = findall(==("data"), tr.source)
    @test length(meas) == 28
    res = [bubble_depth(BubbleAbsorption(radius = 2.88, turn_on = 0.3, picket_ratio = tr.energy[i],
                                         fill = tr.fill[i], stagnation = Inf), tr.t[i]) - tr.x_bub[i]
           for i in meas]
    @test sqrt(sum(abs2, res) / length(res)) < 0.045 && maximum(abs, res) < 0.1
    ab = izumi2021_absorption()
    mz = BubbleAbsorption(radius = 2.875, turn_on = 0.3, footprint_z = 3.0)
    ii = findall(inner(quads))
    for i in findall(==("hydra"), ab.source)
        @test sum(bubble_fractions(mz, quads, ab.t[i]).absorbed[ii]) / length(ii) ≈ ab.absorbed[i] atol = 0.03
    end
    @test 0.3 < bubble_transmission(m3, 5.5) < 0.4             # Izumi's own model: 0.29
    @test 0.5 < bubble_transmission(m15, 5.5) < 0.62           # and 0.54
    @test bubble_depth(m3, 5.5) ≈ 0.36 * (5.5 - 0.55)
    @test bubble_speed(m15) / bubble_speed(m3) ≈ 2^0.2
    @test bubble_speed(BubbleAbsorption(radius = 2.88, picket_ratio = 1.5)) ≈ 0.36 * 1.5^0.2
    @test m3.density_exponent == 1.8                            # IB ∝ ρ², S ∝ ρ₀^-0.1
    @test bubble_transmission(m3, 0.4) == 1 && bubble_depth(m3, 0.4) == 0
    @test bubble_transmission(BubbleAbsorption(radius = 2.88, stagnation = Inf), 20.0) == 0

    # stagnation: ballistic (bit for bit) up to the knee, C¹ there, bounded by x_max
    mb = BubbleAbsorption(radius = 2.88, turn_on = 0.3, stagnation = Inf)
    @test bubble_depth(m3, 5.5) == bubble_depth(mb, 5.5)
    @test bubble_depth(m15, 5.5) == bubble_depth(BubbleAbsorption(radius = 2.88, turn_on = 0.3,
                                                                   fill = 0.15, stagnation = Inf), 5.5)
    xmax, w = m3.stagnation, m3.stagnation_width
    @test xmax ≈ 0.80 * 2.88 && w ≈ 0.08 * 2.88
    tk = m3.launch + (xmax - w) / bubble_speed(m3)               # time at the knee
    @test bubble_depth(m3, tk) ≈ xmax - w
    d(m, t) = ForwardDiff.derivative(τ -> bubble_depth(m, τ), t)
    @test d(m3, tk + 1e-9) ≈ bubble_speed(m3) rtol = 1e-6
    @test all(bubble_depth(m3, t) < xmax for t in 0.0:0.5:20.0)
    @test bubble_depth(m3, 50.0) <= xmax
    @test bubble_depth(m3, 50.0) ≈ xmax rtol = 1e-12
    @test bubble_transmission(m3, 50.0) > 0
    @test all(diff([bubble_depth(m3, t) for t in 0.0:0.1:20.0]) .>= 0)

    # the ridge stop's timing: Izumi 2018's 5.8 ns for a 0.3 ns turn-on, ~0.4 ns earlier
    # for the 1.18× capsule (outer radius 1.064 mm), 5.45 ns with a 1.5× picket
    mc = BubbleAbsorption(radius = 2.875, turn_on = 0.3)
    @test mc.collision ≈ 5.8
    @test BubbleAbsorption(radius = 2.875, turn_on = 0.3, capsule_radius = 1.064).collision ≈ 5.4 atol = 0.05
    @test BubbleAbsorption(radius = 2.875, turn_on = 0.3, picket_ratio = 1.5).collision ≈ 5.45 atol = 0.01
    # a wider capsule–wall gap collides later (HYBRID-E: 3.2 mm wall, 1.135 mm capsule)
    @test BubbleAbsorption(radius = 3.2, turn_on = 0.4, capsule_radius = 1.135).collision ≈
          0.4 + 5.5 * (3.2 - 1.135) / (2.875 - 0.908)
    @test m3.collision ≈ 0.3 + 5.5 * (2.88 - 0.908) / (2.875 - 0.908)
    @test ridge_fraction(mc, 5.8 - mc.ridge_rise) == 0 && ridge_fraction(mc, 5.5) == 0
    @test ridge_fraction(mc, 5.8) ≈ 0.5
    @test ridge_fraction(mc, 5.8 + mc.ridge_rise) == 1
    @test ridge_fraction(BubbleAbsorption(radius = 2.88, collision = Inf), 1e3) == 0

    # the clipped fraction against a closed form: a Gaussian spot (1/e radius a) on a
    # vertical, nearly collimated beam keeps exp(-ρ²/a²) of its power outside radius ρ
    b = Beam("v"; theta = 0.0, phi = 0.0, spot = SuperGaussianSpot(0.5), fnumber = 1e6)
    for ρ in (0.2, 0.5, 0.8)
        @test iris_fraction(b, -1.0, ρ; edge = 0.01, nr = 400, nangle = 8, ncone = 1) ≈
              exp(-(ρ / 0.5)^2) atol = 3e-3
    end

    m = BubbleAbsorption(radius = 2.875)
    # before the launch nothing moves
    @test [x.power for x in apply_lpi(m, quads, 0.1)] == [x.power for x in quads]
    # power is conserved, in each hemisphere, and moves from inner to 50° beams only
    for t in (3.0, 5.0, 7.0)
        eff = apply_lpi(m, quads, t)
        for h in (1, -1)
            side(bs) = sum(x.power for x in bs if (x.port[3] > 0) == (h > 0))
            @test side(eff) ≈ side(quads) rtol = 1e-13
        end
        @test all(x.power < y.power for (x, y) in zip(eff, quads) if x.group in (:inner23, :inner30))
        @test all(x.power > y.power for (x, y) in zip(eff, quads) if x.group === :outer50)
        @test all(x.power == y.power for (x, y) in zip(eff, quads) if x.group === :outer44)
    end
    # the absorbed fraction grows with time, from zero to all of the inner power when
    # ballistic, to a plateau below it with the stagnation
    ts = 0.0:0.25:12.0
    absorbed(mm) = [sum(bubble_fractions(mm, quads, t).absorbed) / count(inner(quads)) for t in ts]
    A, Ab = absorbed(m), absorbed(BubbleAbsorption(radius = 2.875, stagnation = Inf))
    @test A[1] == 0
    @test all(diff(A) .>= 0) && all(diff(Ab) .>= 0)
    @test Ab[end] > 0.99
    @test A[end] < 0.98 && A[end] - A[end-4] < 1e-3                 # stagnated
    @test all(A .<= Ab .+ 1e-15)
    fr = bubble_fractions(m, quads, 5.0)
    @test fr.footprint_z[1] ≈ 5.01 - 2.875 / tand(50) rtol = 1e-3    # 50° axes reach the wall
    @test fr.footprint_z[2] ≈ -fr.footprint_z[1]
    @test all(fr.absorbed .<= fr.clipped)
    # a moving wall moves the footprint plane toward the entrance hole
    mw = BubbleAbsorption(radius = 2.875, wall = WallMotion(speed = 0.1))
    @test bubble_fractions(mw, quads, 5.0).footprint_z[1] > fr.footprint_z[1]

    # differentiable in the calibration (speed, opacity) and in the pointing
    function moved(v, τ, dz)
        bs = [x.group === :inner23 ? remake(x; aim = x.aim .+ (0, 0, sign(x.port[3]) * dz)) : x for x in quads]
        eff = apply_lpi(remake(m; speed = v, opacity = τ), bs, 5.0)
        return sum(x.power for x in eff if x.group === :outer50)
    end
    p = [0.35, 0.144, 0.0]
    g = ForwardDiff.gradient(q -> moved(q...), p)
    for i in 1:3
        δ = zeros(3)
        δ[i] = 1e-5
        @test g[i] ≈ (moved((p .+ δ)...) - moved((p .- δ)...)) / 2e-5 rtol = 1e-5
    end
    @test g[1] > 0 && g[2] > 0 && g[3] != 0

    # in a drive history the moved power raises the pole drive late in the pulse
    t = hohlraum(length = 10.1, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9,
                 h = 0.8, capsule_h = 0.4, wall = Material(:Au; albedo = 0.9, conversion = 0.8))
    model, op = RadiationModel(t), DepositionOperator(t, quads)
    times = [1.0, 6.0]
    h0 = drive_history(model, op, τ -> quads, times)
    h1 = drive_history(model, op, τ -> quads, times; lpi = m)
    p2(h) = h.moments[3, :] ./ h.moments[1, :]
    @test p2(h1)[1] ≈ p2(h0)[1] atol = 1e-6
    @test p2(h1)[2] > p2(h0)[2] + 0.02

    # the ridge stop: after t_c the inner power left by the bubble lands on a wall band
    # at the height where each inner axis crosses r = 0.73 R; power is conserved
    tr = 7.0
    eff = apply_lpi(m, quads, tr)
    noridge = BubbleAbsorption(radius = 2.875, collision = Inf)
    d0 = lpi_deposit(noridge, op, apply_lpi(noridge, quads, tr), tr)
    @test d0.element == deposit(op, eff).element                   # no ridge, no change
    @test lpi_deposit(m, op, eff, 1.0).element == deposit(op, eff).element   # before t_c
    d1 = lpi_deposit(m, op, eff, tr)
    @test d1.delivered ≈ sum(x.power for x in quads) rtol = 1e-14
    @test sum(d1.element) + sum(d1.backside) ≈ sum(d1.per_beam) rtol = 1e-12
    Pin = sum(x.power for x in eff if x.group in (:inner23, :inner30))
    ri = [i for (i, x) in enumerate(eff) if x.group in (:inner23, :inner30)]
    @test sum(d1.per_beam[ri]) ≈ Pin + sum(deposit(op, [x.group in (:inner23, :inner30) ? remake(x; power = 0.0) : x
                                                       for x in eff]).per_beam[ri]) rtol = 1e-12
    # the inner beams' wall power moves to the band: up the wall (z_r between the
    # inner spots and the entrance hole's side of the waist), none to the caps or capsule
    we = t.comp_elems[findfirst(c -> c.name === :wall, t.components)]
    other = setdiff(1:nelements(t), we)
    zr = [Briareus._axis_wall_z(eff[i], 0.73 * 2.875) for i in ri]
    @test all(0.3 .< abs.(zr) .< 1.5)
    Δ = d1.element .- deposit(op, [x.group in (:inner23, :inner30) ? remake(x; power = 0.0) : x for x in eff]).element
    @test sum(Δ[we]) ≈ Pin rtol = 1e-3          # all the inner power lands on the wall
    zc = [t.elem_centroid[e][3] for e in we]
    @test sum(Δ[we] .* abs.(zc)) / sum(Δ[we]) < 1.5

    # composes with other models, and at most one may deposit
    @test lpi_deposit((Backscatter(:outer44 => 0.0), m), op, eff, tr).element ≈ d1.element
    @test_throws ArgumentError lpi_deposit((m, m), op, eff, tr)

    # differentiable in the stagnation depth and the collision time (and the pointing)
    zbar = [t.elem_centroid[e][3]^2 for e in 1:nelements(t)]
    function ridge_moment(xs, tc, dz)
        bs = [x.group === :inner30 ? remake(x; aim = x.aim + (sign(x.port[3]) * dz) * Vec3(0, 0, 1)) : x for x in quads]
        mm = remake(m; stagnation = xs, collision = tc)
        dd = lpi_deposit(mm, op, apply_lpi(mm, bs, 6.9), 6.9)
        return sum(dd.element .* zbar)
    end
    q = [2.2, 6.85, 0.0]
    gq = ForwardDiff.gradient(v -> ridge_moment(v...), q)
    for i in 1:3
        δ = zeros(3)
        δ[i] = 1e-5
        @test gq[i] ≈ (ridge_moment((q .+ δ)...) - ridge_moment((q .- δ)...)) / 2e-5 rtol = 1e-4
    end
    @test all(gq .!= 0)
end
