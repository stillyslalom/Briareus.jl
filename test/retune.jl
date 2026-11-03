# Retuning a drive to a reference: pulses from a cone fraction, pulse tables, and the
# per-time bounded least-squares retune against a full drive history.

@testset "retune to a reference drive" begin
    t = hohlraum(; length = 10.1, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9,
                 wall = Material(:Au; albedo = WallLossAlbedo(), conversion = 0.8),
                 nz = 16, nphi = 24, capsule_ntheta = 8, capsule_nphi = 16)
    rz = Dict(:outer50 => (0.0, 5.01), :outer44 => (0.0, 5.20), :inner30 => (0.60, 3.63),
              :inner23 => (1.09, 3.09))
    quads = map(nif_quads()) do b
        r, z = rz[b.group]
        φ = port_angles(b)[2]
        remake(b; aim = (-r * cosd(φ), -r * sind(φ), sign(b.port[3]) * z))
    end
    # pulses from the total power and the cone fraction
    tt = [0.0, 1.0, 3.0, 4.0]
    total = [0.0, 100.0, 300.0, 0.0]
    cf = [0.3, 0.3, 0.25, 0.25]
    pul = cone_fraction_pulses(quads, tt, total, cf)
    inner = count(b -> b.group in (:inner23, :inner30), quads)
    @test pul[:inner23](3.0) ≈ 0.25 * 300 / inner                  # at the given times
    @test sum(pul[b.group](2.0) for b in quads) ≈ 200
    @test_throws ArgumentError cone_fraction_pulses(quads, tt, total, [0.3, 0.3, 1.2, 0.2])
    # pulse tables round trip
    f = joinpath(mktempdir(), "pulses.csv")
    write_pulse_table(f, pul)
    back = read_pulse_table(f)
    @test Set(keys(back)) == Set(keys(pul))
    @test all(back[k](τ) ≈ pul[k](τ) for k in keys(pul) for τ in (0.5, 2.0, 3.5))
    # keyed pulses
    qp = Dict(Symbol(b.quad) => PulseShape([0.0, 1.0], [1.0, 2.0]) for b in quads)
    @test with_pulses(quads, qp, 0.5; by = :quad)[1].power ≈ 1.5

    # the reference shot, then a dropped bundle, then the retune
    m = RadiationModel(t; closed = true)
    cm = CapsuleModes(m)
    op = DepositionOperator(t, quads)
    times = 0.5:0.5:3.5
    ref = drive_history(m, op, τ -> with_pulses(quads, pul, τ), times; lmax = 4, ylm = true,
                        keep_albedo = true)
    @test length(ref.albedo) == length(ref.time) && length(ref.albedo[1]) == nelements(t)
    keep = findall(b -> nif_bundle(b) != "31", quads)
    kept = quads[keep]
    opk = op[keep]
    amp1(h) = [mode_amplitude(h.ylm[:, k], 1) for k in eachindex(h.time)]
    dropped = drive_history(m, opk, τ -> with_pulses(kept, pul, τ), ref.time; lmax = 4, ylm = true)
    cap(k, τ) = (0.0, 1.5 * pul[first(b.group for b in kept if Symbol(b.quad) == k)](τ) + 1e-9)
    rt = retune_drive(cm, opk, kept, ref; by = :quad, limits = cap)
    @test size(rt.power) == (length(kept), length(ref.time))
    @test all(rt.power[g, k] <= cap(rt.keys[g], rt.time[k])[2] + 1e-9 for g in eachindex(rt.keys)
              for k in eachindex(rt.time))
    @test all(>=(-1e-12), rt.power)
    @test maximum(abs.(rt.drive .- 1)) < 1e-3
    tuned = drive_history(m, opk, τ -> with_pulses(kept, rt.pulses, τ; by = :quad), ref.time;
                          lmax = 4, ylm = true)
    # the per-time linear model is the history's own: same albedos to the solve's tolerance
    @test maximum(abs.(rt.shape .- tuned.ylm[1:25, :] ./ tuned.ylm[1:1, :])) < 1e-4
    # and the dropped bundle's mode 1 is gone, the drive restored
    @test maximum(amp1(dropped)) > 0.01
    @test maximum(amp1(tuned)) < maximum(amp1(dropped)) / 10
    @test maximum(abs.(tuned.Tr .- ref.Tr)) < 0.5
    # staying near nominal powers: a smaller departure, a larger residual
    nominal(k, τ) = pul[first(b.group for b in kept if Symbol(b.quad) == k)](τ)
    near = retune_drive(cm, opk, kept, ref; by = :quad, limits = cap, nominal, stay = 1e-2)
    dep(r) = maximum(abs(r.power[g, k] / nominal(r.keys[g], r.time[k]) - 1)
                     for g in eachindex(r.keys) for k in eachindex(r.time))
    @test dep(near) < dep(rt)
    @test sum(near.residual) > sum(rt.residual)
    # fixed powers (a quad held at zero) are honored
    q0 = rt.keys[1]
    fixed = retune_drive(cm, opk, kept, ref; by = :quad,
                         limits = (k, τ) -> k == q0 ? (0.0, 0.0) : cap(k, τ))
    @test all(iszero, fixed.power[1, :])
    @test_throws ArgumentError retune_drive(cm, opk, kept, dropped)            # no albedos
    @test_throws ArgumentError retune_drive(cm, opk, kept, ref; stay = 1.0)
end

@testset "retune with unconverted light" begin
    # constant albedos and conversion: each step of a history is one linear solve, so the
    # retune's linear prediction is the history's own to rounding
    t = hohlraum(; length = 10.1, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9,
                 wall = Material(:Au; albedo = 0.8, conversion = 0.8),
                 nz = 8, nphi = 12, capsule_ntheta = 6, capsule_nphi = 8)
    rz = Dict(:outer50 => (0.0, 5.01), :outer44 => (0.0, 5.20), :inner30 => (0.60, 3.63),
              :inner23 => (1.09, 3.09))
    quads = map(nif_quads()) do b
        r, z = rz[b.group]
        φ = port_angles(b)[2]
        remake(b; aim = (-r * cosd(φ), -r * sind(φ), sign(b.port[3]) * z))
    end
    pul = cone_fraction_pulses(quads, [0.0, 1.0, 3.0, 4.0], [0.0, 100.0, 300.0, 0.0], [0.3, 0.3, 0.25, 0.25])
    m = RadiationModel(t; closed = true)
    cm = CapsuleModes(m)
    op = DepositionOperator(t, quads)
    times = [0.5, 1.5, 2.5]
    # mock 1ω and 2ω squares near the 3ω axis: NIF's miss the inside of this hohlraum (the
    # wedged lenses keep them ≥ 2.9 mm off the 3ω aim)
    mock = [UnconvertedColor(:omega1; footprint = 2.0, edge = 0.0, focal_shift = 100.0, absorption = 0.88,
                             conversion = 0.88),
            UnconvertedColor(:omega2; footprint = 1.5, edge = 0.3, focal_shift = 60.0, absorption = 1.0, conversion = 0.9)]
    uc = (omega1 = PulseShape([0.0, 4.0], [0.08, 0.03]), omega2 = 0.01)
    pre =(omega1 = Dict(:outer44 => PulseShape([0.0, 4.0], [0.05, 0.02]), :outer50 => τ -> 0.03,
                         :inner23 => τ -> 0.0, :inner30 => τ -> 0.01),)
    ref = drive_history(m, op, τ -> with_pulses(quads, pul, τ), times; lmax = 4, ylm = true, keep_albedo = true,
                        mband = nothing, unconverted = uc, prepulse = pre, colors = mock)
    # the reference retuned to itself, by group: at its own powers (pinned by the limits)
    # the prediction is its drive and shape (measured 9e-16, 2e-16), and the free retune
    # returns those powers to the bounded least-squares solver's tolerance (measured 1e-9;
    # 1e-8 for the same retune of the 3ω-only reference)
    pinned = retune_drive(cm, op, quads, ref; unconverted = uc, prepulse = pre, colors = mock,
                          limits = (k, τ) -> (pul[k](τ), pul[k](τ)))
    @test maximum(abs.(pinned.drive .- 1)) < 1e-12 && maximum(abs.(pinned.shape .- pinned.target_shape)) < 1e-12
    self = retune_drive(cm, op, quads, ref; unconverted = uc, prepulse = pre, colors = mock)
    Pref = [pul[k](τ) for k in self.keys, τ in times]
    @test isapprox(self.power, Pref; rtol = 1e-7) && maximum(abs.(self.drive .- 1)) < 1e-12
    # a reference with unconverted light needs it in the retune; here it raises the drive
    # by 5 % over the 3ω-only history
    @test_throws ArgumentError retune_drive(cm, op, quads, ref)
    ref3 = drive_history(m, op, τ -> with_pulses(quads, pul, τ), times; lmax = 4, ylm = true, keep_albedo = true,
                         mband = nothing)
    @test maximum(abs.(ref.ylm[1, :] ./ ref3.ylm[1, :] .- 1)) > 0.01
    # a dropped bundle, retuned by quad with a prepulse on every beam: the history at the
    # retuned powers is the retune's prediction
    pre2 = (omega1 = PulseShape([0.0, 4.0], [0.04, 0.02]), omega2 = τ -> 0.002)
    ref2 = drive_history(m, op, τ -> with_pulses(quads, pul, τ), times; lmax = 4, ylm = true, keep_albedo = true,
                         mband = nothing, unconverted = uc, prepulse = pre2, colors = mock)
    keep = findall(b -> nif_bundle(b) != "31", quads)
    kept, opk = quads[keep], op[keep]
    cap(k, τ) = (0.0, 1.5 * pul[first(b.group for b in kept if Symbol(b.quad) == k)](τ) + 1e-9)
    rt = retune_drive(cm, opk, kept, ref2; by = :quad, limits = cap, unconverted = uc, prepulse = pre2, colors = mock)
    tuned = drive_history(m, opk, τ -> with_pulses(kept, rt.pulses, τ; by = :quad), times; lmax = 4, ylm = true,
                          mband = nothing, unconverted = uc, prepulse = pre2, colors = mock)
    @test isapprox(rt.shape, tuned.ylm[1:25, :] ./ tuned.ylm[1:1, :]; rtol = 1e-12, atol = 1e-14)
    @test isapprox(rt.drive .* ref2.ylm[1, :], tuned.ylm[1, :]; rtol = 1e-12)
    @test maximum(abs.(rt.drive .- 1)) < 1e-3
    # 1ω on every other quad (a fraction per beam), retuned by cone: each cone's quads
    # carry two fractions, and the history at the retuned powers is the prediction
    half = (omega1 = [isodd(i) ? uc.omega1 : 0.0 for i in eachindex(quads)], omega2 = uc.omega2)
    ref4 = drive_history(m, op, τ -> with_pulses(quads, pul, τ), times; lmax = 4, ylm = true, keep_albedo = true,
                         mband = nothing, unconverted = half, colors = mock)
    uv = (omega1 = half.omega1[keep], omega2 = uc.omega2)
    rv = retune_drive(cm, opk, kept, ref4; by = :group, unconverted = uv, colors = mock)
    tv = drive_history(m, opk, τ -> with_pulses(kept, rv.pulses, τ), times; lmax = 4, ylm = true, mband = nothing,
                       unconverted = uv, colors = mock)
    @test isapprox(rv.shape, tv.ylm[1:25, :] ./ tv.ylm[1:1, :]; rtol = 1e-12, atol = 1e-14)
    @test isapprox(rv.drive .* ref4.ylm[1, :], tv.ylm[1, :]; rtol = 1e-12)
    # (the light is not vacuous: it raises the drive over the 3ω-only history)
    @test all(ref4.ylm[1, :] .> ref3.ylm[1, :])
    # zero fractions: the retune without unconverted light, bit for bit
    r0 = retune_drive(cm, opk, kept, ref3; by = :quad, limits = cap)
    rz0 = retune_drive(cm, opk, kept, ref3; by = :quad, limits = cap, unconverted = (omega1 = 0.0, omega2 = 0.0),
                 colors = mock)
    @test all(isequal(r0[f], rz0[f]) for f in keys(r0) if f !== :pulses) &&
          all(isequal(r0.pulses[k].p, rz0.pulses[k].p) for k in keys(r0.pulses))
    @test_throws ArgumentError retune_drive(cm, opk, kept, ref3; prepulse = pre2)
end

# The face retunes are checked against quantities computed here from their definitions:
# the area-weighted mean and RMS nonuniformity over the face's elements, from forward
# depositions and solves. Two beams symmetric by a half turn about ẑ make a known answer:
# with the second one's pulse at 0.6 of the first's in the reference, the reference
# target recovers 0.6, and the uniform target gives the two equal powers (the RMS at a
# fixed level is convex in the powers and invariant under the turn), at the level asked.
_face_mean(H, A) = sum(A .* H) / sum(A)
_face_rms(H, A) = sqrt(sum(A .* (H .- _face_mean(H, A)) .^ 2) / sum(A)) / _face_mean(H, A)
# the RMS nonuniformity of `p (v1 + ρ v2)` (`v1`, `v2` each beam's distribution at 1 TW)
# at the level `L`, for the power ratio `ρ`
_rms_at_level(v1, v2, A, L, ρ) =
    _face_rms(L / (_face_mean(v1, A) + ρ * _face_mean(v2, A)) .* (v1 .+ ρ .* v2), A)

@testset "retune a foil's laser irradiance" begin
    t = Target(foil(:bl, (0, 0, 0), (0, 0, 1); width = 2.0, n = 24, two_sided = false))
    beam(i, φ, x) = Beam("b$i"; theta = 40, phi = φ, spot = SuperGaussianSpot(0.8), fnumber = 6.7,
                         group = Symbol(:g, i), aim = (x, 0.0, 0.0))
    beams = [beam(1, 0, 0.25), beam(2, 180, -0.25)]          # each spot off center on its side
    op = DepositionOperator(t, beams)
    times = [0.5, 1.0, 1.5]
    pul = Dict(:g1 => PulseShape([0.0, 2.0], [0.2, 0.6]), :g2 => PulseShape([0.0, 2.0], [0.12, 0.36]))
    ref = (; time = times, deposition = [deposit(op, with_pulses(beams, pul, τ)) for τ in times])
    q = FoilIrradiance(:bl; center = (0, 0, 0), radius = 0.6)
    spot = [e for e in elements(t, :bl) if norm(t.elem_centroid[e]) <= 0.6]
    A = t.elem_area[spot]
    irr(P) = deposit(op, [remake(b; power = P[i]) for (i, b) in enumerate(beams)]).element[spot] ./ A
    # the reference's own distribution: the second pulse's scale recovered
    rt = retune_drive(q, op, beams, ref)
    @test rt.keys == [:g1, :g2]
    @test isapprox(rt.power, [pul[k](τ) for k in rt.keys, τ in times]; rtol = 1e-7)
    # the reported mean, RMS and peak-to-mean against the depositions at the retuned powers
    Is = [irr(rt.power[:, k]) for k in eachindex(times)]
    @test isapprox(rt.level, [_face_mean(I, A) for I in Is]; rtol = 1e-12)
    @test isapprox(rt.rms, [_face_rms(I, A) for I in Is]; rtol = 1e-12)
    @test isapprox(rt.peak, [maximum(I) / _face_mean(I, A) for I in Is]; rtol = 1e-12)
    # the uniform target at a prescribed mean: equal powers at that mean, the RMS there
    # below its value at other power ratios on the same mean, and below the reference's
    L(τ) = 0.5τ
    un = retune_drive(q, op, beams, (; time = times); target = :uniform, level = L)
    v1, v2 = irr([1.0, 0.0]), irr([0.0, 1.0])
    @test isapprox(un.power[1, :], un.power[2, :]; rtol = 1e-7)
    @test isapprox([_face_mean(irr(un.power[:, k]), A) for k in eachindex(times)], L.(times); rtol = 1e-10)
    @test isapprox(un.rms[2], _rms_at_level(v1, v2, A, L(1.0), 1.0); rtol = 1e-7)
    @test all(_rms_at_level(v1, v2, A, L(1.0), ρ) > un.rms[2] for ρ in (0.9, 0.97, 1.03, 1.1))
    @test all(un.rms .< rt.rms)
    # a gap in the reference (nothing fires at 1 ns), held off by the limits: not solved,
    # the other times as before
    gap = (; time = times, deposition = [τ == 1.0 ? deposit(op, [remake(b; power = 0.0) for b in beams]) : d
                                         for (τ, d) in zip(times, ref.deposition)])
    rg = retune_drive(q, op, beams, gap; limits = (k, τ) -> τ == 1.0 ? (0.0, 0.0) : (0.0, Inf))
    @test all(iszero, rg.power[:, 2]) && isnan(rg.residual[2]) && rg.level[2] == 0
    @test isapprox(rg.power[:, [1, 3]], rt.power[:, [1, 3]]; rtol = 1e-10)
    @test_throws ArgumentError retune_drive(q, op, beams, gap)                    # an empty target, free powers
    @test_throws ArgumentError retune_drive(q, op, beams, (; time = times))       # no reference drive
    @test_throws ArgumentError retune_drive(q, op, beams, ref; target = :flat)
    @test_throws ArgumentError retune_drive(FoilIrradiance(:bl; center = (5, 5, 0), radius = 0.1), op, beams, ref)
end

@testset "retune a flat face's x-ray drive" begin
    # a halfraum closed at the bottom, whose bottom cap is the observed face, lit through
    # the top hole; constant albedo and conversion, so each step of a history is one linear
    # solve and the retune's linear prediction is the history's own to rounding
    t = hohlraum(; length = 4.0, radius = 1.5, leh_radius = 1.2, capsule_radius = 0.0, nz = 8, nphi = 16,
                 bottom_end = :closed, wall = Material(:Au; albedo = 0.8, conversion = 0.8))
    beams = [Beam("b1"; theta = 35, phi = 0, group = :g1, aim = (-0.3, 0.0, 0.5)),
             Beam("b2"; theta = 35, phi = 180, group = :g2, aim = (0.3, 0.0, 0.5))]
    m = RadiationModel(t; closed = true)
    op = DepositionOperator(t, beams)
    times = [0.5, 1.0, 1.5]
    pul = Dict(:g1 => PulseShape([0.0, 2.0], [1.0, 3.0]), :g2 => PulseShape([0.0, 2.0], [0.6, 1.8]))
    obs = elements(t, :cap_bottom)
    A = t.elem_area[obs]
    hist(P, ts; kw...) = drive_history(m, op, τ -> with_pulses(beams, P, τ), ts; observe = :cap_bottom,
                                       lmax = 2, mband = nothing, kw...)
    ref = hist(pul, times; observed_flux = true, keep_albedo = true)
    @test size(ref.observed_flux) == (length(obs), length(times))
    q = FaceFlux(m, :cap_bottom)
    rt = retune_drive(q, op, beams, ref)
    @test isapprox(rt.power, [pul[k](τ) for k in rt.keys, τ in times]; rtol = 1e-7)
    # the reported mean and RMS against the history at the retuned powers, from each
    # step's incident flux on the face's elements
    tuned = hist(rt.pulses, times; record = sol -> sol.H[obs])
    @test isapprox(rt.level, [_face_mean(H, A) for H in tuned.record]; rtol = 1e-12)
    @test isapprox(rt.rms, [_face_rms(H, A) for H in tuned.record]; rtol = 1e-12)
    @test rt.rms[1] > 1e-3                                   # (the face is not uniform)
    # the uniform target at 1.2 times the reference's mean: equal powers at that mean (from
    # forward solves at 1 TW on each beam), the RMS a minimum along it
    Lr = [_face_mean(ref.observed_flux[:, k], A) for k in eachindex(times)]
    un = retune_drive(q, op, beams, ref; target = :uniform, level = τ -> 1.2 * Lr[findfirst(==(τ), times)])
    flux(P) = solve(m, deposit(op, [remake(b; power = P[i]) for (i, b) in enumerate(beams)]).element).H[obs]
    v1, v2 = flux([1.0, 0.0]), flux([0.0, 1.0])
    @test isapprox(un.power[1, :], un.power[2, :]; rtol = 1e-7)
    @test isapprox(un.power[1, :] .* (_face_mean(v1, A) + _face_mean(v2, A)), 1.2 .* Lr; rtol = 1e-10)
    @test isapprox(un.rms[2], _rms_at_level(v1, v2, A, 1.2 * Lr[2], 1.0); rtol = 1e-7)
    @test all(_rms_at_level(v1, v2, A, 1.2 * Lr[2], ρ) > un.rms[2] for ρ in (0.9, 0.97, 1.03, 1.1))
    @test all(un.rms .< rt.rms)
    @test_throws ArgumentError retune_drive(q, op, beams, hist(pul, times; keep_albedo = true))   # no face flux
    @test_throws ArgumentError retune_drive(q, op, beams, hist(pul, times; observed_flux = true))  # no albedos
end
