# Thin-foil burnthrough: the Marshak heat front against closed forms and published
# self-similar solutions, the albedo transition, its energy accounting, the factor
# cache, and the per-element tracking in `drive_history`.

@testset "heat front" begin
    f = HeatFront()
    c, a, b = f.coefficient, f.temperature_exponent, f.time_exponent
    # constant temperature: the closed form on any grid
    ts = [0.0, 0.03, 0.1, 0.4, 0.45, 1.0, 2.5, 3.0]
    @test heat_front(f, ts, fill(250.0, length(ts)))[2:end] ≈ c .* 2.5^a .* ts[2:end] .^ b rtol = 1e-13
    @test heat_front(f, ts, zeros(length(ts))) == zeros(length(ts))
    # ramp T = T₀ tᵏ: the integral rule's closed form m = c (1 + k a/b)^(-b) T(t)^a t^b,
    # reached by the trapezoidal integral on a fine grid (the integrand ~ t^(ka/b) has an
    # integrable kink at 0, so the error falls as Δt^(1 + ka/b))
    ramp(k, T0, t) = T0 * t^k
    rule(k) = c * (1 + k * a / b)^(-b)
    for k in (0.123, 0.182)
        ts = range(0, 4; length = 40001)
        m = heat_front(f, ts, ramp.(k, 100.0, ts))
        @test m[end] ≈ rule(k) * (4^k)^a * 4^b rtol = 1e-4
    end
    # against the exact self-similar solutions of the subsonic gold wave: Hammer & Rosen
    # (Phys. Plasmas 10, 1829 (2003), eq. 112) give m₀ = 9.51e-4 at q₀ = 0 (our c) and
    # 7.13e-4 at q₀ = 1 (k = 0.182); Shussman & Heizler (Phys. Plasmas 22, 082109 (2015),
    # eqs. 35-36) give 10.17e-4 at constant T and 8.32e-4 at constant flux (k = 0.123).
    # The rule's ratios m(k)/m(0) are 2.2% and 0.7% high: a single drive integral cannot
    # carry the dependence of the profile shape on the drive history.
    @test rule(0.182) / c ≈ 7.13 / 9.51 rtol = 0.025
    @test rule(0.123) / c ≈ 8.32 / 10.17 rtol = 0.01
    # and the ablated mass at 4 ns for T_s = t^0.182 heV, 2.2-2.4e-3 g/cm² in Hammer &
    # Rosen's HYDRA runs and analytic curve (their Fig. 4b)
    @test 2.2e-3 < rule(0.182) * (4^0.182)^a * 4^b < 2.45e-3
end

@testset "burnthrough albedo" begin
    @test_throws ArgumentError BurnthroughAlbedo(thickness = 0.0)
    @test_throws ArgumentError BurnthroughAlbedo(thickness = 1.0, width = -0.1)
    bt = BurnthroughAlbedo(thickness = 2.0)
    @test areal_mass(bt) ≈ 19.3 * 2e-4
    f = bt.front
    # constant temperature: t_b = (m / c T^a)^(1/b)
    tb300 = (areal_mass(bt) / (f.coefficient * 3^f.temperature_exponent))^(1 / f.time_exponent)
    @test burnthrough_time(bt, [0.0, 0.1, 5.0], [300.0, 300.0, 300.0]) ≈ tb300 rtol = 1e-13
    @test burnthrough_time(bt, [0.0, 0.1], [300.0, 300.0]) == Inf
    # the crossing inside a step solves the trapezoid's quadratic exactly
    r(T) = (T / 100)^(f.temperature_exponent / f.time_exponent)
    times, Tr = [0.0, 2.0], [150.0, 320.0]
    tb = burnthrough_time(bt, times, Tr)
    s = (r(Tr[2]) - r(Tr[1])) / 2
    @test r(Tr[1]) * tb + s * tb^2 / 2 ≈ (areal_mass(bt) / f.coefficient)^(1 / f.time_exponent) rtol = 1e-12
    # a thicker foil burns through later; the estimate ignores an explicit time
    @test burnthrough_time(BurnthroughAlbedo(thickness = 3.6), [0.0, 5.0], [300.0, 300.0]) > tb300
    @test burnthrough_time(BurnthroughAlbedo(thickness = 2.0, burnthrough = 0.1), [0.0, 5.0],
                           [300.0, 300.0]) ≈ tb300
    # the albedo: wall-like before t_b, `after` from t_b + width, a smooth step between
    w = WallLossAlbedo()
    @test albedo(bt, 3.0, 250.0) == albedo(w, 3.0, 250.0)          # intact without a time
    @test albedo(bt, 3.0, 250.0; burnthrough = 3.0) == albedo(w, 3.0, 250.0)
    @test albedo(bt, 3.0, 250.0; burnthrough = 2.5) == 0.0
    @test albedo(bt, 3.0, 250.0; burnthrough = 2.8) ≈ albedo(w, 3.0, 250.0) / 2
    @test albedo(bt, 3.0, 250.0; burnthrough = NaN) == albedo(w, 3.0, 250.0)
    be = BurnthroughAlbedo(thickness = 2.0, burnthrough = 1.0, after = 0.2, width = 0.0)
    @test albedo(be, 0.99, 250.0) == albedo(w, 0.99, 250.0) && albedo(be, 1.0, 250.0) ≈ 0.2
    @test albedo(be, 0.5, 250.0; burnthrough = 0.4) ≈ 0.2          # per-element time wins
    bc = BurnthroughAlbedo(thickness = 1.0, before = 0.8, burnthrough = 2.0)
    @test bc.before isa ConstantAlbedo && albedo(bc, 1.0) == 0.8 && albedo(bc, 3.0) == 0.0

    # --- in a hohlraum: solve, the factor cache, and the energy balance
    kw = (; nz = 12, nphi = 16, capsule_ntheta = 8, capsule_nphi = 12)
    wall = Material(:Au; albedo = WallLossAlbedo(), conversion = 0.8)
    foil(; kwargs...) = Material(:Au; albedo = BurnthroughAlbedo(; thickness = 2.0, kwargs...), conversion = 0.8)
    t = hohlraum(; kw..., wall, windows = [WallWindow(phi = 0, diameter = 2.0, material = foil(burnthrough = 1.0))])
    m = RadiationModel(t; closed = true)
    src = elements(t, :wall)
    laser = zeros(nelements(t))
    laser[src] .= t.elem_area[src] .* 0.02
    ew = elements(t, :window_1)
    early, late = solve(m, laser; time = 0.9), solve(m, laser; time = 2.0)
    # intact: the wall's Marshak albedo at the window's own temperature (to the iteration's
    # tolerance), and nothing transmitted
    wl(s) = albedo.(Ref(WallLossAlbedo()), s.time, radiation_temperature.(s.H[ew]))
    @test early.albedo[ew] == early.intact[ew]
    @test maximum(abs, early.albedo[ew] .- wl(early)) < 1e-4
    @test all(==(0.0), late.albedo[ew]) && maximum(abs, late.intact[ew] .- wl(late)) < 1e-4
    # radiosity from the cache matches the direct solve with the converged albedo
    S = Briareus.conversions(t) .* laser ./ t.elem_area
    for s in (early, late)
        B = (I - s.albedo .* m.F) \ S
        @test maximum(abs.(B .- s.B)) <= 1e-10 * maximum(B)
    end
    # per-element burnthrough times override the model's
    pe = solve(m, laser; time = 2.0, burnthrough = fill(5.0, nelements(t)))
    @test pe.albedo[ew] == pe.intact[ew] && maximum(abs, pe.albedo[ew] .- wl(pe)) < 1e-4
    for s in (early, late)
        bal = energy_balance(s)
        @test bal.source ≈ sum(values(bal.absorbed)) + sum(values(bal.escaped)) +
                           sum(values(bal.transmitted)) + bal.lost rtol = 1e-12
        αb = s.intact[ew]
        Pin = s.H[ew] .* t.elem_area[ew]
        @test bal.absorbed[:window_1] ≈ sum((1 .- αb) .* Pin) rtol = 1e-12
        @test bal.transmitted[:window_1] ≈ sum((αb .- s.albedo[ew]) .* Pin) atol = 1e-14
    end
    @test energy_balance(early).transmitted[:window_1] == 0
    @test energy_balance(late).transmitted[:window_1] > 0.5 * energy_balance(late).absorbed[:window_1]
    # a burned-through window costs drive
    c = elements(t, :capsule)
    @test sum(late.H[c]) < sum(solve(m, laser; time = 2.0, burnthrough = fill(Inf, nelements(t))).H[c])

    # --- drive_history tracks each element's front from its own flux history
    q = [remake(b; aim = (0.0, 0.0, sign(b.port[3]) * 4.0)) for b in nif_quads()]
    pulse = PulseShape([0.0, 0.5, 2.0, 3.0, 5.0], [0.0, 1.0, 1.0, 8.0, 8.0])
    tw = hohlraum(; kw..., wall, windows = [WallWindow(phi = 0, diameter = 2.0, gap = 0.05, material = foil()),
                                             WallWindow(phi = 180, diameter = 2.0, material = foil(burnthrough = 3.5))])
    mw = RadiationModel(tw; closed = true)
    op = DepositionOperator(tw, q)
    times = 0.25:0.125:5.0
    h = drive_history(mw, op, τ -> with_pulses(q, pulse, τ), times)
    tb = h.burnthrough[:window_1]
    @test all(isfinite, tb) && all(t -> 2.0 < t < 5.0, tb)
    @test h.burnthrough[:window_2] == fill(3.5, length(elements(tw, :window_2)))
    # re-solving each time with the final burnthrough times reproduces the march (times
    # past τ do not touch its albedos), and the element's own temperature history gives
    # back its burnthrough time through the standalone estimate
    e1 = elements(tw, :window_1)
    tbv = fill(NaN, nelements(tw))
    tbv[e1] .= tb
    Te = zeros(length(times), length(e1))
    for (k, τ) in enumerate(times)
        s = solve(mw, deposit(op, with_pulses(q, pulse, τ)); time = τ, burnthrough = tbv)
        Te[k, :] .= radiation_temperature.(s.H[e1])
        if k == length(times)
            @test energy_balance(s).transmitted[:window_1] ≈ h.transmitted[:window_1][k] rtol = 1e-10
        end
    end
    model = tw.components[component_index(tw, :window_1)].material.albedo
    @test [burnthrough_time(model, collect(times), Te[:, i]) for i in eachindex(e1)] ≈ tb rtol = 1e-8
    @test h.front[:window_1][end] > 1 && h.front[:window_1][1] < 0.2
    # window loss: what the burned foils transmit plus what escapes through the open gap
    @test h.window_loss ≈ h.transmitted[:window_1] .+ h.transmitted[:window_2] .+ h.escaped[:window_1_gap]
    @test all(==(0), h.transmitted[:window_1][times .< minimum(tb)])
    @test h.transmitted[:window_2][end] > 0 && h.transmitted[:window_2][times .<= 3.5] == zeros(count(times .<= 3.5))
    @test all(>(0), h.escaped[:leh_top]) && issorted(h.Tr[1:5])
end
