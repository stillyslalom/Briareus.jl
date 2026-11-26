# This file tests the beam rules-of-engagement scan (`roe_scan`, roe.jl). The linear path
# is checked against direct recomputation with the dropped beams at zero power: a fresh
# deposition, a radiosity solve and moments of the solved flux, none of which use the
# tables of the scan or the response matrix of the scenario. Symmetric cases with known
# answers, the tier rules on synthetic criteria, and the exports follow.

using Briareus, LinearAlgebra, Test

# The beams with the indices `drop` at zero power.
_roe_dropped(beams, drop) = [i in drop ? remake(b; power = 0.0) : b for (i, b) in enumerate(beams)]
_roe_unit(r, name) = findfirst(u -> u.name == name, r.units)
_roe_pair(r, a, b) = only(p for p in r.pairs if Set(p.units) == Set((_roe_unit(r, a), _roe_unit(r, b))))

@testset "ROE: capsule in a hohlraum" begin
    t = hohlraum(; length = 10.1, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9,
                 wall = gold(albedo = 0.85, conversion = 0.8), h = 0.7, capsule_h = 0.25)
    quads = [remake(b; aim = (0.0, 0.0, sign(b.port[3]) * (b.group in (:inner23, :inner30) ? 3.3 : 5.1)))
             for b in nif_quads(power = 8.0)]
    s = Scenario(t, quads; observe = :capsule, viewfactor_kw = (; closed = true))
    wall_power(ev) = sum(ev.deposition.element[e] for e in elements(ev.scenario.target, :wall))
    crit = [roe_drive(), roe_mode(2), roe_mode(4), roe_hotspot_velocity(), roe_yield(),
            ROECriterion(:wall_power, wall_power; green = 0.03, amber = 0.06, measure = :loss, label = "fraction")]
    r = roe_scan(s, quads, crit)
    @test length(r.units) == 48 && all(u -> u.kind === :beam, r.units)   # single-beam quads: no quad units
    @test length(r.pairs) == 48 * 47 ÷ 2 && r.linear == [true]
    @test all(u -> u.role === :wall, r.units)

    # The reference of an outage: deposit, solve and project the solved capsule flux.
    obs = s.observed
    Wy = ylm_weights(t, obs; lmax = 1)
    function direct(drop)
        dep = deposit(s.deposition, _roe_dropped(quads, drop))
        H = solve(s.model, dep).H
        a = legendre_moments(H, t, obs; lmax = 4)
        v = hotspot_velocity(Wy * H[obs])
        return (; a, v, wall = sum(dep.element[elements(t, :wall)]))
    end
    d0 = direct(Int[])
    function reference(drop)
        d = direct(drop)
        return [(d0.a[1] - d.a[1]) / d0.a[1], abs(d.a[3] / d.a[1] - d0.a[3] / d0.a[1]),
                abs(d.a[5] / d.a[1] - d0.a[5] / d0.a[1]), norm(d.v - d0.v),
                1 - implosion_metrics(d.v).total / implosion_metrics(d0.v).total,
                (d0.wall - d.wall) / d0.wall]
    end
    close(o) = isapprox(o.scores, reference(o.beams); rtol = 1e-8, atol = 1e-10)
    @test all(close, r.singles[[1, 7, 20, 33, 48]])             # measured 2e-13 relative
    @test all(close, r.pairs[[1, 100, 555, 1128]])
    # the names of the units of a pair are those of its beams
    p = _roe_pair(r, quads[3].name, quads[40].name)
    @test p.beams == [3, 40] && close(p)

    # Every outage evaluated in full agrees with the linear scan.
    rf = roe_scan(s, quads, crit; method = :full, pairs = :none)
    @test rf.linear == [false]
    @test all(isapprox(a.scores, b.scores; rtol = 1e-9, atol = 1e-10) for (a, b) in zip(r.singles, rf.singles))
    @test [o.tier for o in rf.singles] == [o.tier for o in r.singles]

    # A dropped quad of 48 equal ones removes 1/48 of the wall power when every beam
    # lands on the wall, and the drive loss of the quads sums to the whole drive.
    k = findfirst(c -> c.name === :wall_power, r.criteria)
    @test sum(o.scores[k] for o in r.singles) ≈ 1 rtol = 1e-12
    @test sum(o.scores[1] for o in r.singles) ≈ 1 rtol = 1e-10
    # the tier is the worst criterion tier and the binding criterion is one that sets it
    @test all(o.tier === Briareus._ROE_TIERS[maximum(Briareus._roe_rank, o.tiers)] for o in r.singles)
    @test all(o.tiers[o.binding] === o.tier for o in r.singles)

    # Mirror companions pair each quad with one in the other hemisphere, and the two
    # together cancel most of the axial part of mode 1 of each alone. Opposite companions
    # lie near the antipode and cancel most of the whole vector (measured: 7% to 31% of
    # the single-quad velocity remains).
    vel = [roe_hotspot_velocity(; measure = :change)]
    rm = roe_scan(s, quads, vel; pairs = :none, companions = :mirror)
    partner(r, i) = (c = r.companions[i]; c.units[1] == i ? c.units[2] : c.units[1])
    @test all(i -> sign(quads[i].port[3]) == -sign(quads[partner(rm, i)].port[3]), 1:48)
    vz(drop) = hotspot_velocity(evaluate(s, _roe_dropped(quads, drop)))[3]
    @test all(i -> abs(vz([i, partner(rm, i)])) < 0.35 * abs(vz([i])), (1, 9, 30))
    ro = roe_scan(s, quads, vel; pairs = :none, companions = :opposite)
    @test all(i -> dot(quads[i].port, quads[partner(ro, i)].port) < -0.98, 1:48)
    @test all(i -> ro.companions[i].scores[1] < 0.35 * ro.singles[i].scores[1], 1:48)
end

@testset "ROE: face of a package in a halfraum" begin
    hr = hohlraum(; length = 3.0, radius = 1.5, leh_radius = 1.0, capsule_radius = 0.0, h = 0.3,
                  top_end = 1.0, bottom_end = :open)
    smp = foil(:sample, (0, 0, -1.5), (0, 0, 1); width = 1.0, n = 6, two_sided = false,
               material = Material(:sample; albedo = 0.3, conversion = 0.0))
    t = Target([hr.components; smp]; closed = false)
    up = [remake(b; aim = (0.0, 0.0, 1.5)) for b in nif_quads(power = 4.0) if b.port[3] > 0]
    s = Scenario(t, up; observe = :sample)
    sw = remake(s; observe = (:wall, :cap_top))
    crit = [roe_face_flux(), roe_face_nonuniformity(), roe_temperature(),
            roe_temperature(; name = :wall_Tr, case = :wall)]
    r = roe_scan((sample = (s, up), wall = (sw, up)), crit; pairs = :none)
    @test r.cases == [:sample, :wall] && length(r.units) == 24
    A = t.elem_area[s.observed]
    function direct(drop)
        sol = solve(s.model, deposit(s.deposition, _roe_dropped(up, drop)))
        H = sol.H[s.observed]
        return (; mean = sum(A .* H) / sum(A), rms = rms_nonuniformity(H, A),
                Tr = radiation_temperature(sum(A .* H) / sum(A)), wall = wall_temperature(sol, :wall, :cap_top))
    end
    d0 = direct(Int[])
    reference(drop) = (d = direct(drop); [(d0.mean - d.mean) / d0.mean, d.rms - d0.rms, abs(d.Tr - d0.Tr),
                                          abs(d.wall - d0.wall)])
    @test all(isapprox(o.scores, reference(o.beams); rtol = 1e-8, atol = 1e-11) for o in r.singles[[1, 5, 12, 24]])
    @test r.baseline[3] ≈ d0.Tr rtol = 1e-10
    @test r.baseline[4] ≈ d0.wall rtol = 1e-10

    # The drive balance between the sample and the end cap, against the area means of the
    # flux of a direct solve.
    sb = remake(s; observe = (:sample, :cap_top))
    rb = roe_scan(sb, up, [roe_flux_balance(:sample, (:cap_top,))]; pairs = :none)
    ea, eb = elements(t, :sample), elements(t, :cap_top)
    function balance(drop)
        H = solve(sb.model, deposit(sb.deposition, _roe_dropped(up, drop))).H
        Ha = sum(t.elem_area[ea] .* H[ea]) / sum(t.elem_area[ea])
        Hb = sum(t.elem_area[eb] .* H[eb]) / sum(t.elem_area[eb])
        return (Ha - Hb) / (Ha + Hb)
    end
    b0 = balance(Int[])
    @test rb.baseline[1] ≈ b0 rtol = 1e-10
    @test abs(b0) > 0.01                                   # the two faces see different drives
    @test all(isapprox(o.scores[1], abs(balance(o.beams) - b0); rtol = 1e-7, atol = 1e-11) for o in rb.singles[[1, 5, 12, 24]])
    @test any(o -> o.scores[1] > 1e-5, rb.singles)
    # exchanging the two sets changes the sign of the balance and leaves every score
    rx = roe_scan(sb, up, [roe_flux_balance(:cap_top, :sample)]; pairs = :none)
    @test rx.baseline[1] ≈ -b0 rtol = 1e-10
    @test [o.scores[1] for o in rx.singles] ≈ [o.scores[1] for o in rb.singles] rtol = 1e-9
    @test_throws ArgumentError roe_scan(s, up, [roe_flux_balance(:sample, :cap_top)]; pairs = :none)
end

@testset "ROE: foils, symmetric cases and tiers" begin
    # A foil facing +z under four beams at one polar angle, 90° apart in azimuth, each
    # aimed 0.5 mm off center toward its azimuth. A second foil takes two more beams.
    main = foil(:backlighter, (0, 0, 0), (0, 0, 1); width = 6.0, n = 12, two_sided = false)
    side = foil(:second, (20.0, 0, 0), (0, 0, 1); width = 6.0, n = 6, two_sided = false)
    t = Target([main, side])
    az = (0.0, 90.0, 180.0, 270.0)
    four = [Beam("B$i"; theta = 30.0, phi = φ, power = 2.0, spot = SuperGaussianSpot(0.4),
                 aim = (0.5 * cosd(φ), 0.5 * sind(φ), 0.0), quad = i <= 2 ? "QA" : "QB", group = :main)
            for (i, φ) in enumerate(az)]
    two = [Beam("S$i"; theta = 20.0, phi = φ, power = 1.0, spot = SuperGaussianSpot(0.4),
                aim = (20.0, 0.0, 0.0), group = :side) for (i, φ) in enumerate((45.0, 225.0))]
    beams = [four; two]
    s = Scenario(t, beams)
    fe = elements(t, :backlighter)
    Af = sum(t.elem_area[fe])
    los = (0.0, 0.0)
    sig = (; line = :Ge, pulse = 2.0, distance = 500.0, object = nothing)
    crit = [roe_foil_irradiance(:backlighter), roe_foil_centroid(:backlighter),
            roe_foil_nonuniformity(:backlighter), roe_backlighter_signal(los; sig...)]
    r = roe_scan(s, beams, crit; units = :beams, pairs = :beams, companions = :opposite)
    @test [u.role for u in r.units] == [:backlighter, :backlighter, :backlighter, :backlighter, :second, :second]

    # Each of the four beams lands whole on the foil, so one drop removes a quarter of
    # the irradiance, and a beam of the other foil removes nothing.
    @test all(i -> isapprox(r.singles[i].scores[1], 0.25; rtol = 1e-9), 1:4)   # measured 1e-15
    @test all(i -> abs(r.singles[i].scores[1]) < 1e-14, 5:6)
    # The centroid of one beam alone, from a deposition of that beam through the
    # restricted operator: the nominal centroid is the mean of the four, and one drop
    # moves it to the mean of the other three.
    cen(i) = (d = deposit(s.deposition[[i]], [beams[i]]);
              sum(d.element[e] * t.elem_centroid[e] for e in fe) / sum(d.element[fe]))
    c = [cen(i) for i in 1:4]
    @test all(i -> isapprox(r.singles[i].scores[2], norm(sum(c[j] for j in 1:4 if j != i) / 3 - sum(c) / 4);
                            rtol = 1e-9), 1:4)
    @test r.singles[1].scores[2] ≈ norm(c[1]) / 3 rtol = 1e-3       # by symmetry, c₀ = 0
    # Opposed beams are companions, and the pair leaves the centroid where it was.
    partner(i) = (o = r.companions[i]; r.units[o.units[1] == i ? o.units[2] : o.units[1]].name)
    @test [partner(i) for i in 1:6] == ["B3", "B4", "B1", "B2", "S2", "S1"]
    @test all(i -> r.companions[i].scores[2] < 1e-3 * r.singles[i].scores[2], 1:4)
    @test _roe_pair(r, "B1", "B3").scores[2] == r.companions[1].scores[2]
    # nonuniformity and signal against direct evaluation of a fresh deposition
    function reference(drop)
        d0, d = deposit(s.deposition, beams), deposit(s.deposition, _roe_dropped(beams, drop))
        rms(x) = rms_nonuniformity(x.element[fe] ./ t.elem_area[fe], t.elem_area[fe])
        ph(x) = backlighter_signal(t, x, los; sig...).photons
        return [rms(d) - rms(d0), (ph(d0) - ph(d)) / ph(d0)]
    end
    @test all(isapprox(o.scores[3:4], reference(o.beams); rtol = 1e-9, atol = 1e-12)
              for o in [r.singles[[1, 4, 6]]; _roe_pair(r, "B1", "B2")])
    # the conversion efficiency falls with the irradiance, so the signal loses more than a quarter
    @test r.singles[1].scores[4] > 0.25

    # Quads, epochs and beams without power. In the second epoch only B3 and B4 fire.
    late = [b.name in ("B3", "B4") ? remake(b; power = 3.0) : remake(b; power = 0.0) for b in beams]
    idle = [beams; Beam("idle"; theta = 60.0, phi = 10.0, power = 0.0, spot = SuperGaussianSpot(0.4))]
    si = Scenario(t, idle)
    ce = [roe_foil_irradiance(:backlighter), roe_foil_irradiance(:backlighter; name = :late, case = :late)]
    late_i = [late; idle[end]]
    re = roe_scan([:early => (si, idle), :late => (si, late_i)], ce)
    @test [u.name for u in re.units] == ["B1", "B2", "B3", "B4", "S1", "S2", "QA", "QB"]
    @test [u.kind for u in re.units[7:8]] == [:quad, :quad] && re.units[7].beams == [1, 2]
    @test length(re.pairs) == 1 && re.pairs[1].beams == [1, 2, 3, 4]       # the one quad pair
    score(name, k) = re.singles[_roe_unit(re, name)].scores[k]
    @test [score(n, 2) for n in ("B1", "B3", "QA", "QB", "S1")] ≈ [0.0, 0.5, 0.0, 1.0, 0.0] atol = 1e-12
    @test [score(n, 1) for n in ("B1", "QA", "QB")] ≈ [0.25, 0.5, 0.5] rtol = 1e-9
    @test re.pairs[1].scores ≈ [1.0, 1.0] rtol = 1e-12

    # A laser–plasma model switches the scan to full evaluations. Backscatter acts on
    # each beam separately, so the frozen tables give the same scores, and the loss of
    # a backscattered beam is its absorbed share.
    sl = Scenario(t, four; lpi = Backscatter(:main => 0.2))
    half = [i == 1 ? remake(b; group = :other) : b for (i, b) in enumerate(four)]
    sh = Scenario(t, half; lpi = Backscatter(:main => 0.2))
    rl = roe_scan(sh, half, [roe_foil_irradiance(:backlighter)]; units = :beams, pairs = :none)
    rt = roe_scan(sh, half, [roe_foil_irradiance(:backlighter)]; units = :beams, pairs = :none, method = :linear)
    @test rl.linear == [false] && rt.linear == [true]
    @test [o.scores[1] for o in rl.singles] ≈ [1.0, 0.8, 0.8, 0.8] ./ 3.4 rtol = 1e-9
    @test [o.scores[1] for o in rt.singles] ≈ [o.scores[1] for o in rl.singles] rtol = 1e-12
    @test roe_scan(sl, four, [roe_foil_irradiance(:backlighter)]; pairs = :none).linear == [false]

    # Tier rules on synthetic criteria of the set of dropped beams.
    s4 = Scenario(t, four)
    dropped(ev) = findall(b -> b.power == 0, ev.beams)
    count_c(green, amber) = [ROECriterion(:count, ev -> length(dropped(ev)); green, amber, measure = :value)]
    scan(c; kw...) = roe_scan(s4, four, c; units = :beams, pairs = :beams, kw...)
    # every pair red: all four raised to amber, each with three caveats
    ra = scan(count_c(1, 1.5))
    @test all(o -> o.tier === :green, ra.singles) && all(==(:amber), ra.tiers)
    @test all(c -> length(c) == 3, ra.caveats) && only(ra.budgets).limit == 0
    # every pair amber: the units stay green and drop one at a time
    rb = scan(count_c(1, 2))
    @test all(==(:green), rb.tiers) && all(c -> length(c) == 3, rb.caveats)
    @test (b = only(rb.budgets); (b.green, b.pairs, b.safe, b.limit) == (4, 6, 0, 1))
    # every pair green: two together
    rg = scan(count_c(2, 3))
    @test all(isempty, rg.caveats) && (b = only(rg.budgets); (b.safe, b.limit) == (6, 2))
    # one red unit: its pairs are no worse than the unit and add no caveat
    w = [1.0, 1.0, 1.0, 5.0]
    rw = scan([ROECriterion(:weight, ev -> sum(w[dropped(ev)]; init = 0.0); green = 2, amber = 4.5, measure = :value)])
    @test rw.tiers == [:green, :green, :green, :red] && all(isempty, rw.caveats)
    @test (b = only(rw.budgets); (b.green, b.pairs, b.safe, b.limit) == (3, 3, 3, 2))
    # one pair leaves green: an exception of the budget and a caveat of its two units
    rx = scan([ROECriterion(:special, ev -> dropped(ev) == [1, 2] ? 3.0 : 1.0; green = 2, amber = 4, measure = :value)])
    @test only(rx.budgets).exceptions == [("B1", "B2")] && only(rx.budgets).limit == 2
    @test length.(rx.caveats) == [1, 1, 0, 0] && all(==(:green), rx.tiers)
    # pair scopes and roles: with B1 in another role, three pairs lie within a role
    roles = Dict("B1" => :cbi)
    @test length(scan(count_c(2, 3); roles, pair_scope = :within).pairs) == 3
    @test length(scan(count_c(2, 3); roles, pair_scope = :across).pairs) == 3
    @test [u.role for u in scan(count_c(2, 3); roles = :group, pairs = :none).units] == fill(:main, 4)
    @test [u.role for u in scan(count_c(2, 3); roles = b -> Symbol(b.quad), pairs = :none).units] == [:QA, :QA, :QB, :QB]
    rp = roe_scan(s4, four, count_c(2, 3); pairs = :none, extra_pairs = [("B1", "QB")], companions = ["QA" => "QB"])
    @test only(rp.pairs).beams == [1, 3, 4] && rp.companions[_roe_unit(rp, "QA")].beams == [1, 2, 3, 4]
    # measures: a vector value scores by its norm, a NaN score is red
    vec_c = ROECriterion(:vec, ev -> [3.0, 4.0] .* length(dropped(ev)); green = 5, amber = 10, measure = :change)
    nan_c = ROECriterion(:nan, ev -> isempty(dropped(ev)) ? 1.0 : NaN; green = 1, amber = 2, measure = :value)
    rv = scan([vec_c, nan_c]; pairs = :none)
    @test all(o -> o.scores[1] == 5.0 && o.tiers == [:green, :red] && o.binding == 2, rv.singles)
    rel = ROECriterion(:rel, ev -> 4.0 - length(dropped(ev)); green = 0.3, amber = 0.6, measure = :relative)
    @test scan([rel]).pairs[1].scores == [0.5]
    # a hypothetical outage of any units: a pair of a full scan is reproduced exactly, and
    # three beams, or a quad with a beam, drop together
    rf = scan([rel, count_c(1, 2)[1]]; method = :full)
    o2 = roe_outage(rf, s4, four, ["B1", "B3"])
    @test o2.scores == _roe_pair(rf, "B1", "B3").scores && o2.units == [1, 3] && o2.beams == [1, 3]
    o3 = roe_outage(rf, s4, four, ("B1", :B2, "B4"))
    @test o3.scores == [0.75, 3.0] && o3.tiers == [:red, :red] && o3.tier === :red
    rq = roe_scan(s4, four, [rel]; pairs = :none)
    oq = roe_outage(rq, [:default => (s4, four)], ["QA", "B3"])
    @test oq.beams == [1, 2, 3] && oq.units == sort([_roe_unit(rq, "QA"), _roe_unit(rq, "B3")]) && oq.scores == [0.75]
    @test_throws ArgumentError roe_outage(rq, s4, four, ["B9"])
    @test_throws ArgumentError roe_outage(rq, s4, four, String[])
    @test_throws ArgumentError roe_outage(rq, [:other => (s4, four)], ["B1"])

    # errors
    @test_throws ArgumentError ROECriterion(:a, identity; green = 2, amber = 1)
    @test_throws ArgumentError ROECriterion(:a, identity; green = 1, amber = 2, measure = :worst)
    @test_throws ArgumentError scan(ROECriterion[])
    @test_throws ArgumentError scan([rel, rel])
    @test_throws ArgumentError scan([ROECriterion(:a, ev -> 1.0; green = 1, amber = 2, case = :missing)])
    @test_throws ArgumentError scan([rel]; units = :ports)
    @test_throws ArgumentError scan([rel]; companions = :nearest)
    @test_throws ArgumentError scan([rel]; extra_pairs = [("B1", "B9")])
    @test_throws ArgumentError roe_scan([:a => (s4, four), :b => (s, beams)], [rel])
    @test_throws ArgumentError roe_scan(s4, [remake(b; power = 0.0) for b in four], [rel])
    @test_throws ArgumentError roe_mode(0)
    @test_throws Exception scan([ROECriterion(:v, ev -> [1.0, 2.0]; green = 1, amber = 2, measure = :loss)])

    # Tables. The rows run by role, then from critical to auto-drop.
    rows = roe_table(r)
    @test length(rows) == 6 && [x.role for x in rows] == [fill(:backlighter, 4); fill(:second, 2)]
    @test issorted([Briareus._roe_rank(x.tier) for x in rows[1:4]]; rev = true)
    x1 = only(x for x in rows if x.unit == "B1")
    @test x1.beams == ["B1"] && x1.companion == "B3" && x1.scores.irradiance ≈ 0.25
    @test x1.score == r.singles[1].scores[r.singles[1].binding] && x1.decision in ("Critical", "Ask RI", "Auto-Drop")
    @test length(roe_table(re; kind = :quad)) == 2
    mktempdir() do dir
        path = write_roe_csv(joinpath(dir, "roe.csv"), r)
        head, lines = read_csv(path)
        @test head[1:8] == ["role", "unit", "kind", "beams", "tier", "decision", "single", "binding"]
        @test head[12:15] == ["irradiance", "centroid", "irradiance_rms", "backlighter"] && length(lines) == 6
        @test [l[2] for l in lines] == [x.unit for x in rows]
        @test parse.(Float64, [l[12] for l in lines]) ≈ [x.scores.irradiance for x in rows] rtol = 1e-9
        md = read(write_roe_markdown(joinpath(dir, "roe.md"), rx; title = "Shot sheet"), String)
        @test startswith(md, "# Shot sheet") && occursin("## Role: backlighter", md)
        @test occursin("up to 2 green beams together, except 1 of 6 pairs", md)
        @test occursin("Pairs that leave green: B1 + B2.", md) && occursin("| B1 | beam | Auto-Drop |", md)
        @test occursin("B2: amber (special)", md)
        @test occursin("Ask RI (alone: green)", roe_markdown(ra))
        @test count(==('\n'), md) == count(==('\n'), roe_markdown(rx; title = "Shot sheet"))
    end
    @test occursin("4 beams: 0 red, 4 amber, 0 green", sprint(show, ra))
end

@testset "ROE: direct drive" begin
    R = 0.5
    t = Target(Component(:capsule, icosphere_mesh(R; n = 6), Material(:CH); surface = SphereSurface(R)))
    beams = omega_beams(; cones = :all, spot = omega_spot(:SG4), power = 1.0)
    s = Scenario(t, beams; deposition_kw = (; ncone = 1))
    L = 8
    r = roe_scan(s, beams, [roe_illumination_rms(t; lmax = L)]; pairs = :none, companions = :opposite)
    @test length(r.units) == 60
    W = ylm_weights(t, 1:nelements(t); lmax = L)
    σ(drop) = mode_rms(W * (deposit(s.deposition, _roe_dropped(beams, drop)).element ./ t.elem_area), L)
    # the default score is the rise over the full shot
    @test all(isapprox(o.scores[1], σ(o.beams) - σ(Int[]); rtol = 1e-8) for o in r.singles[[1, 17, 60]])
    rv = roe_scan(s, beams, [roe_illumination_rms(t; lmax = L, measure = :value)]; pairs = :none)
    @test all(isapprox(o.scores[1], σ(o.beams); rtol = 1e-9) for o in rv.singles[[1, 17, 60]])
    @test r.baseline[1] ≈ σ(Int[]) rtol = 1e-10
    # OMEGA beams come in opposed pairs, and a companion is the antipode
    @test all(i -> dot(beams[r.companions[i].beams[1]].port, beams[r.companions[i].beams[2]].port) < -0.999, 1:60)
    # a lost beam raises the nonuniformity
    @test all(o -> o.scores[1] > 0, r.singles)

    # The corona operator of three beams against its deposition with one beam off.
    prof = CoronaProfile(; critical_radius = 0.40, scale_length = 0.08, absorption = 1 - exp(-2.4))
    three = beams[[1, 20, 45]]
    op = CoronaOperator(prof, three; ncone = 3)
    s3 = Scenario(t, three; deposition_kw = (; ncone = 1))
    rc = roe_scan(s3, three, [roe_corona_rms(op; lmax = 4)]; pairs = :none)
    σc(drop) = mode_rms(deposit(op, _roe_dropped(three, drop)).ylm, 4)
    @test all(isapprox(o.scores[1], σc(o.beams) - σc(Int[]); rtol = 1e-9, atol = 1e-12) for o in rc.singles)
end
