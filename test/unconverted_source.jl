# Per-color x-ray sources in `solve`: a two-plate enclosure solved by hand, the 3ω path
# bit for bit, and the equivalent-3ω weighting for conversions independent of intensity.

@testset "per-color x-ray sources" begin
    uc = nif_unconverted()
    r1, r2 = uc[1].conversion, uc[2].conversion

    # Two directly opposed unit squares 1 mm apart, one element each. View factor between
    # them: the closed form for directly opposed parallel rectangles at X = Y =
    # side/spacing = 1 (Howell's catalog of configuration factors).
    X = 1.0
    fx = 2 / (π * X^2) * (log(sqrt((1 + X^2)^2 / (1 + 2X^2))) + 2X * sqrt(1 + X^2) * atan(X / sqrt(1 + X^2)) -
                          2X * atan(X))
    @test isapprox(fx, 0.19982; atol = 1e-5)
    F = [0.0 fx; fx 0.0]
    # plate 1: Dewald 2008's intensity table; plate 2: a time-dependent conversion
    tab = Conversion(intensities = (1e14, 1e15), efficiencies = (0.95, 0.80))
    m1 = Material(:Au; albedo = 0.6, conversion = tab)
    m2 = Material(:Au; albedo = 0.3, conversion = Conversion(efficiency = 0.5, time_exponent = 0.2))
    sq = rectangle_mesh(1.0, 1.0; nu = 1, nv = 1)
    t2 = Target(Component(:a, sq, m1),
                Component(:b, transform(sq, [1 0 0; 0 -1 0; 0 0 -1], Vec3(0, 0, 1)), m2))
    @test nelements(t2) == 2 && all(≈(1.0), t2.elem_area)
    m = RadiationModel(t2, F)
    P = [3.0, 0.4]                                 # TW on 1 mm²: 3e14 and 4e13 W/cm²
    u = (omega1 = [0.5, 0.1], omega2 = [2.0, 0.0])
    τ = 2.0
    sol = solve(m, P; time = τ, unconverted = u, colors = uc)
    # by hand: η(I) = 0.95 - 0.15 log₁₀(I/1e14), held at 0.95 below 1e14 W/cm²
    η(I) = 0.95 - 0.15 * log10(clamp(I, 1e14, 1e15) / 1e14)
    S1 = η(3e14) * 3.0 + r1 * 0.95 * 0.5 + r2 * η(2e14) * 2.0
    S2 = 0.5 * τ^0.2 * (0.4 + r1 * 0.1)
    α1, α2 = 0.6, 0.3
    B1 = (S1 + α1 * fx * S2) / (1 - α1 * α2 * fx^2)
    B2 = S2 + α2 * fx * B1
    @test isapprox(sol.source, [S1, S2]; rtol = 1e-12)
    @test isapprox(sol.B, [B1, B2]; rtol = 1e-12) && isapprox(sol.H, [fx * B2, fx * B1]; rtol = 1e-12)
    @test sol.laser == P .+ (u.omega1 .+ u.omega2)
    # the equivalent 3ω power, at its own intensity, is another source on plate 1
    Seq = solve(m, P .+ r1 .* u.omega1 .+ r2 .* u.omega2; time = τ).source
    @test abs(Seq[1] - S1) > 1e-3 * S1 && isapprox(Seq[2], S2; rtol = 1e-12)
    # a `Deposition` per color and a `Dict` give the same solution
    dz = zeros(2)
    ud = Dict(:omega1 => Deposition(u.omega1, [0.0, 0.0], [0.6], 0.6), :omega2 => Deposition(u.omega2, dz, [2.0], 2.0))
    @test isequal(solve(m, P; time = τ, unconverted = ud).B, sol.B)
    # derivative in the 3ω power on plate 1: d(ηP)/dP = η + dη/dlnI = η(I) - 0.15/ln 10
    dB = ForwardDiff.derivative(x -> solve(m, [x, 0.4]; time = τ, unconverted = u).B, 3.0)
    dS1 = η(3e14) - 0.15 / log(10)
    @test isapprox(dB, [dS1, α2 * fx * dS1] ./ (1 - α1 * α2 * fx^2); rtol = 1e-12)
    @test_throws ArgumentError solve(m, P; unconverted = (omega1 = u.omega1,))
    @test_throws DimensionMismatch solve(m, P; unconverted = (omega1 = [1.0], omega2 = [1.0]))

    # A hohlraum: the 3ω path bit for bit, and the equivalent 3ω power for conversions
    # independent of intensity (constant gold; nif_gold's time-dependent conversion with a
    # temperature-dependent albedo, iterated)
    th = hohlraum(nz = 8, nphi = 12, capsule_ntheta = 6, capsule_nphi = 8)
    mh = RadiationModel(th)
    rng = Xoshiro(11)
    Ph = 0.1 .* rand(rng, nelements(th))
    uh = (omega1 = 0.01 .* rand(rng, nelements(th)), omega2 = 0.002 .* rand(rng, nelements(th)))
    zero_u = (omega1 = zeros(nelements(th)), omega2 = zeros(nelements(th)))
    s0, sz = solve(mh, Ph), solve(mh, Ph; unconverted = zero_u)
    @test all(isequal(getfield(s0, f), getfield(sz, f)) for f in (:laser, :source, :B, :H, :albedo))
    Peq = Ph .+ r1 .* uh.omega1 .+ r2 .* uh.omega2
    se, sc = solve(mh, Peq), solve(mh, Ph; unconverted = uh)
    @test isapprox(sc.B, se.B; rtol = 1e-12) && isapprox(sc.H, se.H; rtol = 1e-12)
    @test sc.laser == Ph .+ (uh.omega1 .+ uh.omega2)
    mw = RadiationModel(hohlraum(nz = 8, nphi = 12, capsule_ntheta = 6, capsule_nphi = 8, wall = nif_gold()); F = mh.F)
    sew, scw = solve(mw, Peq; time = 3.0), solve(mw, Ph; time = 3.0, unconverted = uh)
    @test isapprox(scw.H, sew.H; rtol = 1e-12) && isapprox(scw.albedo, sew.albedo; rtol = 1e-12)
    # linear at fixed albedo and constant conversion: d H / d(scale of P) = H(P), with
    # the unconverted source fixed (the dual path)
    obs = th.comp_elems[findfirst(c -> c.name === :capsule, th.components)]
    dH = ForwardDiff.derivative(x -> sum(solve(mh, x .* Ph; unconverted = uh).H[obs]), 1.0)
    @test isapprox(dH, sum(s0.H[obs]); rtol = 1e-12)
end

# Two directly opposed unit squares one side length apart: Howell's closed form at X = 1.
const _UC_FX = let X = 1.0
    2 / (π * X^2) * (log(sqrt((1 + X^2)^2 / (1 + 2X^2))) + 2X * sqrt(1 + X^2) * atan(X / sqrt(1 + X^2)) - 2X * atan(X))
end

@testset "per-color flux response and converted share" begin
    uc = nif_unconverted()
    r1, r2 = uc[1].conversion, uc[2].conversion
    fx = _UC_FX
    F = [0.0 fx; fx 0.0]
    tab = Conversion(intensities = (1e14, 1e15), efficiencies = (0.95, 0.80))
    m1 = Material(:Au; albedo = 0.6, conversion = tab)
    m2 = Material(:Au; albedo = 0.3, conversion = Conversion(efficiency = 0.5, time_exponent = 0.2))
    sq = rectangle_mesh(1.0, 1.0; nu = 1, nv = 1)
    t2 = Target(Component(:a, sq, m1),
                Component(:b, transform(sq, [1 0 0; 0 -1 0; 0 0 -1], Vec3(0, 0, 1)), m2))
    m = RadiationModel(t2, F)
    P = [3.0, 0.4]
    u = (omega1 = [0.5, 0.1], omega2 = [2.0, 0.0])
    τ = 2.0
    sol = solve(m, P; time = τ, unconverted = u)
    R = flux_response(m, [1, 2]; time = τ, laser = P, unconverted = u)
    @test R isa Dict{Symbol,Matrix{Float64}} && sort(collect(keys(R))) == [:omega1, :omega2, :omega3]
    @test isequal(R[:omega3], flux_response(m, [1, 2]; time = τ, laser = P))
    Hr(u) = R[:omega3] * P .+ R[:omega1] * u.omega1 .+ R[:omega2] * u.omega2
    @test isapprox(Hr(u), sol.H; rtol = 1e-12)
    # away from the secant: the 2ω flux on plate 1 rises from 2e14 to 2.2e14 W/cm², on
    # the table; the 1ω flux (5.5e13) stays below it and plate 2's conversion depends on
    # time only. The response's error is the source change r₂ [η(2e14) − η(2.2e14)] 2.2
    # on plate 1 through the 2×2 radiosity solve.
    u2 = (omega1 = 1.1 .* u.omega1, omega2 = 1.1 .* u.omega2)
    η(I) = 0.95 - 0.15 * log10(clamp(I, 1e14, 1e15) / 1e14)
    ΔS = r2 * (η(2e14) - η(2.2e14)) * 2.2
    ΔB1 = ΔS / (1 - 0.6 * 0.3 * fx^2)
    ΔH = [fx * 0.3 * fx * ΔB1, fx * ΔB1]
    @test isapprox(Hr(u2) .- solve(m, P; time = τ, unconverted = u2).H, ΔH; rtol = 1e-9)
    # the converted share is the per-color source (no re-emission here), where the 3ω
    # conversion of the total power differs
    @test isequal(sol.converted, sol.source)
    @test spectral_solve(sol).converted == sol.source
    @test abs(Briareus._xray_source(t2, sol.laser, τ)[1] - sol.source[1]) > 1e-3 * sol.source[1]

    # a conversion independent of intensity: the response is the per-color solve at any
    # powers
    th = hohlraum(nz = 8, nphi = 12, capsule_ntheta = 6, capsule_nphi = 8)
    mh = RadiationModel(th)
    rng = Xoshiro(5)
    Ph = 0.1 .* rand(rng, nelements(th))
    uh = (omega1 = 0.01 .* rand(rng, nelements(th)), omega2 = 0.002 .* rand(rng, nelements(th)))
    obs = elements(th, :capsule)
    Rh = flux_response(mh, obs; unconverted = uh)
    Hh(s) = Rh[:omega3] * Ph .+ s .* (Rh[:omega1] * uh.omega1 .+ Rh[:omega2] * uh.omega2)
    Hs(s) = solve(mh, Ph; unconverted = (omega1 = s .* uh.omega1, omega2 = s .* uh.omega2)).H[obs]
    @test isapprox(Hh(1.0), Hs(1.0); rtol = 1e-12) && isapprox(Hh(1.7), Hs(1.7); rtol = 1e-12)
    @test_throws ArgumentError flux_response(mh, obs; unconverted = (omega1 = uh.omega1,))
end
