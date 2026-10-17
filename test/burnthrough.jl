# Thin-foil burnthrough: the Marshak heat front against closed forms and published
# self-similar solutions, the post-burnthrough transmission, the albedo transition, its
# energy accounting, the factor cache, and the per-element tracking in `drive_history`.

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

@testset "foil transmission" begin
    # Lester et al., Phys. Rev. E 114, L023201 (2026), eq. 2: F_trans = (4/3) F_inc / (τ + 2/3)
    for τ in (2 / 3, 1.0, 3.4, 5.0, 7.2, 20.0)
        @test eddington_transmission(τ) ≈ (4 / 3) / (τ + 2 / 3) rtol = 1e-15
    end
    @test eddington_transmission(5.0) == 4 / 17                # τ ∼ 5 at 2.5 µm: 23.5%
    @test eddington_transmission(0.0) == eddington_transmission(0.5) == 1    # capped below τ = 2/3
    @test_throws DomainError eddington_transmission(-1.0)
    @test ForwardDiff.derivative(eddington_transmission, 5.0) ≈ -12 / 17^2 rtol = 1e-14
    # the default: τ linear in areal mass from τ = 5 at 2.5 µm of solid gold
    bt = BurnthroughAlbedo(thickness = 2.5)
    @test optical_depth(bt) == 5 && bt.transmission == 4 / 17 && bt.after === nothing
    for x in (0.8, 1.7, 3.6)
        m = BurnthroughAlbedo(thickness = x)
        @test optical_depth(m) ≈ 5x / 2.5 rtol = 1e-14
        @test m.transmission ≈ 4 / (3 * 5x / 2.5 + 2) rtol = 1e-14
    end
    @test BurnthroughAlbedo(thickness = 1.7).transmission ≈ 0.328 atol = 1e-3
    @test BurnthroughAlbedo(thickness = 3.6).transmission ≈ 0.169 atol = 1e-3
    @test optical_depth(BurnthroughAlbedo(thickness = 5.0, density = 9.65)) ≈ 5   # same areal mass
    # an explicit transmission, or an explicit `after` instead
    @test BurnthroughAlbedo(thickness = 2.5, transmission = 0.4).transmission == 0.4
    bo = BurnthroughAlbedo(thickness = 2.5, after = 0)
    @test bo.transmission === nothing && bo.after isa ConstantAlbedo
    @test_throws ArgumentError BurnthroughAlbedo(thickness = 2.5, after = 0, transmission = 0.2)
    @test_throws ArgumentError BurnthroughAlbedo(thickness = 2.5, transmission = 1.2)
    @test_throws ArgumentError BurnthroughAlbedo(thickness = 2.5, transmission = -0.1)
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
    # the albedo: wall-like before t_b, (1 - T) α_before from t_b + width (the transmitted
    # part T leaves, the rest re-emits as before), a smooth step s between, along which
    # the transmitted fraction is s T
    w = WallLossAlbedo()
    αw = albedo(w, 3.0, 250.0)
    T = bt.transmission                                             # 2 µm: τ = 4, T = 2/7
    @test T ≈ 2 / 7 rtol = 1e-15
    @test albedo(bt, 3.0, 250.0) == αw                              # intact without a time
    @test albedo(bt, 3.0, 250.0; burnthrough = 3.0) == αw
    @test albedo(bt, 3.0, 250.0; burnthrough = 2.5) == (1 - T) * αw
    @test albedo(bt, 3.0, 250.0; burnthrough = 2.8) ≈ αw * (1 - T / 2) rtol = 1e-15
    @test albedo(bt, 3.0, 250.0; burnthrough = NaN) == αw
    parts(m, tb) = Briareus._burnthrough_split(m, 3.0, 250.0, tb)
    @test parts(bt, 3.0) == (αw, αw, 0.0)                           # nothing transmitted intact
    @test parts(bt, 2.5) == ((1 - T) * αw, αw, T)
    @test all(parts(bt, 2.8) .≈ (αw * (1 - T / 2), αw, T / 2))
    # an explicit `after`: the albedo steps to it, and the loss beyond the intact foil's
    # counts as transmitted. `after = 0` is the open hole, the default before
    # `transmission` (these are the old default's checks).
    bo = BurnthroughAlbedo(thickness = 2.0, after = 0.0)
    @test albedo(bo, 3.0, 250.0; burnthrough = 2.5) == 0.0
    @test albedo(bo, 3.0, 250.0; burnthrough = 2.8) ≈ αw / 2
    @test parts(bo, 2.5) == (0.0, αw, αw) && parts(bo, 3.0) == (αw, αw, 0.0)
    @test albedo(BurnthroughAlbedo(thickness = 2.0, transmission = 1.0), 3.0, 250.0; burnthrough = 2.5) == 0.0
    be = BurnthroughAlbedo(thickness = 2.0, burnthrough = 1.0, after = 0.2, width = 0.0)
    @test albedo(be, 0.99, 250.0) == albedo(w, 0.99, 250.0) && albedo(be, 1.0, 250.0) ≈ 0.2
    @test albedo(be, 0.5, 250.0; burnthrough = 0.4) ≈ 0.2          # per-element time wins
    bc = BurnthroughAlbedo(thickness = 1.0, before = 0.8, burnthrough = 2.0)   # τ = 2, T = 1/2
    @test bc.before isa ConstantAlbedo && albedo(bc, 1.0) == 0.8 && albedo(bc, 3.0) == 0.4
    @test albedo(BurnthroughAlbedo(thickness = 1.0, before = 0.8, burnthrough = 2.0, after = 0), 3.0) == 0.0

    # --- in a hohlraum: solve, the factor cache, and the energy balance
    kw = (; nz = 12, nphi = 16, capsule_ntheta = 8, capsule_nphi = 12)
    wall = Material(:Au; albedo = WallLossAlbedo(), conversion = 0.8)
    foil(; kwargs...) = Material(:Au; albedo = BurnthroughAlbedo(; thickness = 2.0, kwargs...), conversion = 0.8)
    windowed(material) = hohlraum(; kw..., wall, windows = [WallWindow(; phi = 0, diameter = 2.0, material)])
    t = windowed(foil(burnthrough = 1.0))
    m = RadiationModel(t; closed = true)
    src = elements(t, :wall)
    laser = zeros(nelements(t))
    laser[src] .= t.elem_area[src] .* 0.02
    ew = elements(t, :window_1)
    T = 2 / 7                                       # the default for 2 µm (τ = 4)
    # early: intact; mid: halfway through the transition (s = 1/2); late: burned through
    early, mid, late = (solve(m, laser; time) for time in (0.9, 1.2, 2.0))
    # intact: the wall's Marshak albedo at the window's own temperature (to the iteration's
    # tolerance), and nothing transmitted
    wl(s) = albedo.(Ref(WallLossAlbedo()), s.time, radiation_temperature.(s.H[ew]))
    @test early.albedo[ew] == early.intact[ew] && all(==(0.0), early.transmission[ew])
    @test maximum(abs, early.albedo[ew] .- wl(early)) < 1e-4
    # burned through: a fraction T transmitted, the rest re-emitted with the wall's albedo
    @test maximum(abs, late.intact[ew] .- wl(late)) < 1e-4
    @test late.albedo[ew] == (1 - T) .* late.intact[ew] && all(≈(T; rtol = 1e-15), late.transmission[ew])
    @test mid.albedo[ew] ≈ (1 - T / 2) .* mid.intact[ew] rtol = 1e-15
    @test all(≈(T / 2; rtol = 1e-15), mid.transmission[ew])
    @test all(==(0.0), late.transmission[setdiff(1:nelements(t), ew)])
    # radiosity from the cache matches the direct solve with the converged albedo
    S = Briareus.conversions(t) .* laser ./ t.elem_area
    for s in (early, late)
        B = (I - s.albedo .* m.F) \ S
        @test maximum(abs.(B .- s.B)) <= 1e-10 * maximum(B)
    end
    # per-element burnthrough times override the model's
    pe = solve(m, laser; time = 2.0, burnthrough = fill(5.0, nelements(t)))
    @test pe.albedo[ew] == pe.intact[ew] && maximum(abs, pe.albedo[ew] .- wl(pe)) < 1e-4
    # the energy balance closes; the foil transmits s T of its incident power and absorbs
    # (1 - s T)(1 - α_before) of it
    for (s, sT) in ((early, 0.0), (mid, T / 2), (late, T))
        bal = energy_balance(s)
        @test bal.source ≈ sum(values(bal.absorbed)) + sum(values(bal.escaped)) +
                           sum(values(bal.transmitted)) + bal.lost rtol = 1e-12
        αb = s.intact[ew]
        Pin = s.H[ew] .* t.elem_area[ew]
        @test bal.absorbed[:window_1] ≈ sum((1 - sT) .* (1 .- αb) .* Pin) rtol = 1e-12
        @test bal.transmitted[:window_1] ≈ sT * sum(Pin) rtol = 1e-14 atol = 1e-15
        @test bal.absorbed[:window_1] + bal.transmitted[:window_1] ≈ sum((1 .- s.albedo[ew]) .* Pin) rtol = 1e-13
    end
    @test energy_balance(early).transmitted[:window_1] == 0
    @test energy_balance(late).transmitted[:window_1] > 0.5 * energy_balance(late).absorbed[:window_1]
    # a burned-through window costs drive
    c = elements(t, :capsule)
    @test sum(late.H[c]) < sum(solve(m, laser; time = 2.0, burnthrough = fill(Inf, nelements(t))).H[c])

    # the open hole, `after = 0`: the old default, albedo and accounting alike (the
    # intact foil's loss stays absorbed, the excess is transmitted)
    to = windowed(foil(burnthrough = 1.0, after = 0.0))
    mo = RadiationModel(to; F = m.F)
    earlyo, lateo = solve(mo, laser; time = 0.9), solve(mo, laser; time = 2.0)
    @test earlyo.albedo == early.albedo && earlyo.B == early.B      # identical while intact
    @test all(==(0.0), lateo.albedo[ew]) && lateo.transmission[ew] == lateo.intact[ew]
    for s in (earlyo, lateo)
        bal = energy_balance(s)
        @test bal.source ≈ sum(values(bal.absorbed)) + sum(values(bal.escaped)) +
                           sum(values(bal.transmitted)) + bal.lost rtol = 1e-12
        αb = s.intact[ew]
        Pin = s.H[ew] .* t.elem_area[ew]
        @test bal.absorbed[:window_1] ≈ sum((1 .- αb) .* Pin) rtol = 1e-14
        @test bal.transmitted[:window_1] ≈ sum((αb .- s.albedo[ew]) .* Pin) rtol = 1e-14 atol = 1e-15
    end
    # it costs more drive than the default; `transmission = 1` has the same albedo but
    # books the whole loss as transmitted
    @test sum(lateo.H[c]) < sum(late.H[c])
    th = RadiationModel(windowed(foil(burnthrough = 1.0, transmission = 1.0)); F = m.F)
    lateh = solve(th, laser; time = 2.0)
    @test lateh.albedo == lateo.albedo
    @test lateh.B ≈ lateo.B rtol = 1e-14
    bh, bo = energy_balance(lateh), energy_balance(lateo)
    @test bh.absorbed[:window_1] == 0
    @test bh.transmitted[:window_1] ≈ bo.transmitted[:window_1] + bo.absorbed[:window_1] rtol = 1e-13

    # --- drive_history tracks each element's front from its own flux history
    q = [remake(b; aim = (0.0, 0.0, sign(b.port[3]) * 4.0)) for b in nif_quads()]
    pulse = PulseShape([0.0, 0.5, 2.0, 3.0, 5.0], [0.0, 1.0, 1.0, 8.0, 8.0])
    tw = hohlraum(; kw..., wall, windows = [WallWindow(phi = 0, diameter = 2.0, gap = 0.05, material = foil()),
                                             WallWindow(phi = 180, diameter = 2.0, material = foil(burnthrough = 3.5))])
    mw = RadiationModel(tw; closed = true)
    op = DepositionOperator(tw, q)
    times = 0.25:0.125:5.0
    # each step's albedo iteration starts from the previous step's temperatures; iterated
    # tightly, the march and the cold re-solves below reach the same fixed point
    h = drive_history(mw, op, τ -> with_pulses(q, pulse, τ), times; tol = 1e-9)
    @test h.time == collect(times)
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
        s = solve(mw, deposit(op, with_pulses(q, pulse, τ)); time = τ, burnthrough = tbv, tol = 1e-9)
        Te[k, :] .= radiation_temperature.(s.H[e1])
        if k == length(times)
            @test energy_balance(s).transmitted[:window_1] ≈ h.transmitted[:window_1][k] rtol = 1e-7
        end
        # once a window has burned through (its transition over), the history sends
        # exactly the fraction T of the flux on it to `transmitted`
        for (n, tbn) in ((:window_1, maximum(tb)), (:window_2, 3.5))
            τ >= tbn + 0.4 || continue
            en = elements(tw, n)
            Pin = sum(s.H[en] .* tw.elem_area[en])
            @test energy_balance(s).transmitted[n] ≈ T * Pin rtol = 1e-13
            @test h.transmitted[n][k] ≈ T * Pin rtol = 1e-7
        end
    end
    @test h.time[end] >= maximum(tb) + 0.4          # the check above ran for both windows
    model = tw.components[component_index(tw, :window_1)].material.albedo
    @test [burnthrough_time(model, collect(times), Te[:, i]) for i in eachindex(e1)] ≈ tb rtol = 1e-8
    @test h.front[:window_1][end] > 1 && h.front[:window_1][1] < 0.2
    # window loss: what the burned foils transmit plus what escapes through the open gap
    @test h.window_loss ≈ h.transmitted[:window_1] .+ h.transmitted[:window_2] .+ h.escaped[:window_1_gap]
    @test all(==(0), h.transmitted[:window_1][times .< minimum(tb)])
    @test h.transmitted[:window_2][end] > 0 && h.transmitted[:window_2][times .<= 3.5] == zeros(count(times .<= 3.5))
    @test all(>(0), h.escaped[:leh_top]) && issorted(h.Tr[1:5])

    # --- adaptive steps through the burnthrough transitions (width 0.4 ns): window 1
    # burns through at about 2.46 ns, window 2 at its explicit 3.5 ns
    f = τ -> with_pulses(q, pulse, τ)
    coarse = 0.25:0.25:5.0
    hc = drive_history(mw, op, f, coarse)
    @test hc.time == collect(coarse)                    # by default, exactly the given times
    hr = drive_history(mw, op, f, coarse; transition = 0.1)
    inserted = setdiff(hr.time, coarse)
    @test issorted(hr.time) && issubset(coarse, hr.time) && !isempty(inserted)
    # refined only while a transition is under way, and no step advances one by more
    # than the tolerance (1.5 times it after a step's own correction); without
    # refinement a single step takes most of a transition
    tb1 = hr.burnthrough[:window_1]
    @test all(τ -> minimum(tb1) < τ < maximum(tb1) + 0.4 || 3.5 < τ < 3.9, inserted)
    fm = tw.components[component_index(tw, :window_1)].material.albedo
    advance(h) = maximum(Briareus._burnthrough_step(fm, h.time[k], tb) -
                         Briareus._burnthrough_step(fm, h.time[k-1], tb)
                         for tb in [h.burnthrough[:window_1]; 3.5] for k in 2:length(h.time))
    @test advance(hr) <= 0.15 && advance(hc) > 0.5
    @test length(hr.Tr) == length(hr.window_loss) == size(hr.moments, 2) == length(hr.time)
    # against a uniformly fine history: the refined one, interpolated linearly between
    # its steps, follows the window loss through the transitions 5 times more closely
    fine = 0.25:1/32:5.0
    hf = drive_history(mw, op, f, fine)
    interp(h, τ) = (k = clamp(searchsortedlast(h.time, τ), 1, length(h.time) - 1);
                    w = (τ - h.time[k]) / (h.time[k+1] - h.time[k]);
                    (1 - w) * h.window_loss[k] + w * h.window_loss[k+1])
    err(h) = maximum(abs(interp(h, τ) - x) for (τ, x) in zip(hf.time, hf.window_loss)) /
             maximum(hf.window_loss)
    @test err(hr) < 0.015 && err(hr) < err(hc) / 5
    @test maximum(abs.(hr.burnthrough[:window_1] .- hf.burnthrough[:window_1])) < 2e-3
    @test hr.window_loss[end] ≈ hf.window_loss[end] rtol = 1e-3

    # flux_change: the observed flux shape H/⟨H⟩ moves by at most the tolerance per
    # step (checked by re-solving at the steps with the final burnthrough times). The
    # coarse steps move it by up to 0.045 with the default transmission (over 0.05 for
    # the open hole that was the default before), so the tolerance is 0.01.
    hx = drive_history(mw, op, f, coarse; flux_change = 0.01)
    @test issubset(coarse, hx.time) && length(hx.time) > length(coarse)
    obs = elements(tw, :capsule)
    function shape_change(h)
        tbv = fill(NaN, nelements(tw))
        tbv[elements(tw, :window_1)] .= h.burnthrough[:window_1]
        S = map(h.time) do τ
            s = solve(mw, deposit(op, f(τ)); time = τ, burnthrough = tbv)
            s.H[obs] ./ legendre_moments(s.H, tw, obs; lmax = 0)[1]
        end
        return maximum(maximum(abs.(S[k] .- S[k-1])) for k in 2:length(S))
    end
    @test shape_change(hx) <= 0.01 && shape_change(hc) > 0.03
    @test_throws ArgumentError drive_history(mw, op, f, coarse; transition = 0.0)
    @test_throws ArgumentError drive_history(mw, op, f, [1.0, 0.5])
end
