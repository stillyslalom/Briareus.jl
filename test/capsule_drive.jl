# CapsuleDrive: the drive handed to a capsule-only hydrodynamics run (ROADMAP "Briareus → CompactLES").

@testset "capsule drive" begin
    t = hohlraum(nz = 10, nphi = 16, capsule_ntheta = 8, capsule_nphi = 12,
                 wall = Material(:Au; albedo = WallLossAlbedo(), conversion = 0.8))
    m = RadiationModel(t)
    q = nif_quads()
    pulse = PulseShape([0.0, 1.0, 5.0, 6.0], [0.0, 2.0, 8.0, 0.0])
    times = [0.5, 2.0, 4.0, 5.5]
    # every beam crosses the axis at z = ±4 and lands toward the waist (equator-hot), or
    # at z = ±1 and lands on the far cap and entrance-hole end (pole-hot)
    run(z) = (qz = [remake(b; aim = (0.0, 0.0, sign(b.port[3]) * z)) for b in q];
              op = DepositionOperator(t, qz);
              (drive_history(m, op, τ -> with_pulses(qz, pulse, τ), times; lmax = 8, wall_history = false), op, qz))
    h, op, qz = run(4.0)
    d = capsule_drive(h, m; provenance = (; pulse = "test ramp", hohlraum = "coarse NIF"))
    # shapes: a uniform grid of ≤ 50 ps over the solved times, l = 0…8
    @test d.time[1] == 0.5 && d.time[end] == 5.5
    @test maximum(diff(d.time)) <= 0.05 + 1e-12 && length(d.time) == 101
    @test size(d.modes) == (9, 101) && all(==(1), d.modes[1, :])
    @test length(d.Tr) == length(d.flux) == length(d.absorbed) == 101
    @test d.radius == 0.9 && d.area ≈ sum(t.elem_area[elements(t, :capsule)])
    @test d.provenance[:pulse] == "test ramp" && d.provenance[:code] == "Briareus"
    @test d.provenance[:observe] == [:capsule] && haskey(d.provenance, :git)
    # round trip: on the history's own times the drive is the history
    dh = capsule_drive(h, m; times = h.time)
    @test dh.flux == h.flux && dh.Tr ≈ h.Tr
    @test dh.modes ≈ h.moments ./ h.moments[1:1, :]
    @test dh.absorbed == h.absorbed[:capsule]
    @test d(2.0).flux ≈ h.flux[2] && d(2.0).modes ≈ h.moments[:, 2] ./ h.moments[1, 2]
    @test d(3.0).flux ≈ (h.flux[2] + h.flux[3]) / 2           # linear between steps
    # resampling the resampled drive at the history times recovers it (piecewise linear)
    @test capsule_drive(h, m; times = [2.0, 4.0]).flux ≈ h.flux[2:3]
    @test capsule_drive(h, m; dt = 1.0).time == collect(0.5:1.0:5.5)
    @test size(capsule_drive(h, m; lmax = 4).modes, 1) == 5
    @test_throws ArgumentError capsule_drive(h, m; lmax = 10)
    @test_throws ArgumentError capsule_drive(h, m; times = [0.0, 1.0])
    @test_throws ArgumentError capsule_drive(h, m; observe = :wall)   # not one sphere
    # a₀ is the capsule flux energy_balance books: (1 − α) a₀ A = absorbed (constant albedo)
    α = albedos(t)[elements(t, :capsule)[1]]
    @test all(k -> isapprox(dh.absorbed[k], (1 - α) * dh.flux[k] * d.area; rtol = 2e-3), 1:4)
    sol = solve(m, deposit(op, with_pulses(qz, pulse, 4.0)); time = 4.0)
    @test energy_balance(sol).absorbed[:capsule] ≈ (1 - α) * d(4.0).flux * d.area rtol = 2e-3
    # sign: P2 > 0 is pole-hot. A field that is larger toward the poles…
    es = elements(t, :capsule)
    H = zeros(nelements(t))
    for e in es
        c = t.elem_centroid[e]
        H[e] = 1 + 0.5 * (c[3] / 0.9)^2
    end
    a = legendre_moments(H, t, es; lmax = 8)
    hs = (; time = [0.0, 1.0], moments = hcat(a, a), absorbed = Dict{Symbol,Vector{Float64}}())
    ds = capsule_drive(hs, t)
    @test ds.modes[3, 1] > 0.1 && all(isnan, ds.absorbed)
    # …and the beams: landing on the caps is pole-hot, toward the waist equator-hot
    @test all(d.modes[3, :] .< 0)
    hp, _, _ = run(1.0)
    @test all(capsule_drive(hp, m).modes[3, :] .> 0)
    @test occursin("CapsuleDrive(101 times", sprint(show, d))
end

@testset "capsule drive: Y_lm (all m)" begin
    t = hohlraum(nz = 10, nphi = 16, nr_cap = 3, nr_leh = 2, capsule_mesh = :latlong,
                 capsule_ntheta = 12, capsule_nphi = 16,
                 wall = Material(:Au; albedo = WallLossAlbedo(), conversion = 0.8))
    m = RadiationModel(t)
    es = elements(t, :capsule)
    # A field of known Schmidt ratios: 1 + 0.05 P2 + 0.1 Ŷ22 + 0.03 Ŷ3,-1 (φ from +x)
    Y22(u) = sqrt(3) / 2 * (u[1]^2 - u[2]^2)
    Y3m1(u) = sqrt(1 / 6) * 1.5 * (5u[3]^2 - 1) * u[2]
    H = zeros(nelements(t))
    for e in es
        u = normalize(t.elem_centroid[e])
        H[e] = 1 + 0.05 * (3u[3]^2 - 1) / 2 + 0.1 * Y22(u) + 0.03 * Y3m1(u)
    end
    a = legendre_moments(H, t, es; lmax = 4)
    c = ylm_weights(t, es; lmax = 4) * H[es]
    hs = (; time = [0.0, 1.0], moments = hcat(a, 2a), ylm = hcat(c, 2c),
          absorbed = Dict{Symbol,Vector{Float64}}())
    d = capsule_drive(hs, t; times = [0.0, 0.5, 1.0])
    @test size(d.ylm) == (25, 3) && all(==(1), d.ylm[1, :])
    # m = 0 ↔ aₗ/a₀, exactly, and the replaced rows agreed to quadrature level
    @test all(l -> d.ylm[mode_index(l, 0), :] == d.modes[l+1, :], 0:4)
    @test all(l -> abs(sqrt((2l + 1) / 4π) * c[mode_index(l, 0)] / a[1] - a[l+1] / a[1]) < 1e-3, 1:4)
    # the ratios recover the field's coefficients (element-constant field: ~1e-3)
    q = d.ylm[:, 2]
    @test q[mode_index(2, 2)] ≈ 0.1 atol = 5e-3
    @test q[mode_index(3, -1)] ≈ 0.03 atol = 2e-3
    @test q[mode_index(2, 0)] ≈ 0.05 atol = 2e-3
    others = setdiff(1:25, mode_index.([0, 2, 2, 3], [0, 0, 2, -1]))
    @test maximum(abs, q[others]) < 4e-3
    # √(Σₘ bₗₘ²) is mode_amplitude (before the m = 0 rows are replaced by modes)
    @test mode_amplitude(c, 2) ≈ sqrt(sum(abs2, sqrt(5 / 4π) * c[mode_index(2, mm)] / a[1] for mm in -2:2))
    # round trip: a₀ Σ bₗₘ Ŷₗₘ rebuilds the field (Ŷₗₘ = √(4π/(2l+1)) Yₗₘ, orthonormal Y)
    Y = zeros(5, 9)
    err = 0.0
    for e in es
        u = normalize(t.elem_centroid[e])
        Briareus._real_ylm!(Y, 4, u[3], atan(u[2], u[1]))
        F = d.flux[2] * sum(q[mode_index(l, mm)] * sqrt(4π / (2l + 1)) * Y[l+1, 5+mm]
                            for l in 0:4 for mm in -l:l)
        err = max(err, abs(F / (1.5 * a[1]) - H[e]))
    end
    @test err < 0.01
    @test d(0.5).ylm ≈ d.ylm[:, 2] && d(0.25).flux ≈ 1.25 * a[1]
    @test capsule_drive(hs, t; ylm = false).ylm === nothing
    @test size(capsule_drive(hs, t; lmax = 2).ylm, 1) == 9
    @test occursin("(all m)", sprint(show, d))
    @test_throws ArgumentError CapsuleDrive(d.time, d.Tr, d.flux, d.modes, d.absorbed, d.radius,
                                            d.area, d.provenance, d.ylm .+ 1)
    @test_throws ArgumentError capsule_drive((; hs..., ylm = nothing), t; ylm = true)
    # A beam pattern rotated by 45° about the axis (the mesh is symmetric under it)
    # rotates the moments: (cₗₘ, cₗ,₋ₘ) turn by mα, and m = 0 stays.
    qb = nif_quads()[[1, 2, 7]]                         # three quads: no azimuthal symmetry
    α = π / 4
    Rz = [cos(α) -sin(α) 0; sin(α) cos(α) 0; 0 0 1]
    rot(b) = Briareus.Beam{Float64}(b.name, Briareus.Vec3(Rz * b.port...), typeof(b.aim)(Rz * b.aim...), b.power, b.defocus, b.spot,
                                    b.fnumber, b.group, b.quad, b.wavelength, b.dlambda)
    pulse = PulseShape([0.0, 1.0, 5.0], [0.0, 2.0, 2.0])
    hist(bs) = drive_history(m, DepositionOperator(t, bs), τ -> with_pulses(bs, pulse, τ),
                             [1.0, 2.0]; lmax = 4, ylm = true, wall_history = false)
    d0 = capsule_drive(hist(qb), m)
    d1 = capsule_drive(hist(rot.(qb)), m)
    @test d1.flux ≈ d0.flux rtol = 1e-8
    @test d1.modes ≈ d0.modes rtol = 1e-6 atol = 1e-9
    b0, b1 = d0.ylm[:, end], d1.ylm[:, end]
    @test abs(b0[mode_index(2, 2)]) + abs(b0[mode_index(2, -2)]) > 2e-3      # a real m = 2 drive (0.5 %)
    for l in 1:4, mm in 1:l
        cp, sp = b0[mode_index(l, mm)], b0[mode_index(l, -mm)]
        @test b1[mode_index(l, mm)] ≈ cp * cos(mm * α) - sp * sin(mm * α) atol = 1e-7
        @test b1[mode_index(l, -mm)] ≈ cp * sin(mm * α) + sp * cos(mm * α) atol = 1e-7
    end
end
