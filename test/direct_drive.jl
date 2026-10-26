# Hard-sphere direct-drive illumination (polar direct drive, Stage 0; ROADMAP item 15 and
# reference/litreview/polar_direct_drive.md): the deposition operator on a bare sphere
# with collimated beams (`ncone = 1`) against independent references, and the
# capsule-frame metrics `mode_rms` / `mode_power(C; m)` to l = 30.
#
# References (all independent of the operator and of the mode weights):
# - Skupsky & Lee, J. Appl. Phys. 54, 3662 (1983): for identical axisymmetric beams aimed
#   at the centre, σₗ = |aₗ/a₀| Gₗ / √(2l+1) with aₗ the single beam's Legendre
#   coefficients (a 1-D quadrature here) and Gₗ² = Σᵢⱼ Pₗ(Uᵢ·Uⱼ)/N² from the ports alone.
# - Temporal et al., arXiv:1406.7762 (2014): one ring at 49°/131° on R₀ = 0.5 mm, m = 4
#   super-Gaussians; σ = 11.9 % (circular, a = b = R₀) and 5.2 % (elliptical, a = 2b,
#   b = R₀ meridional) centred, minima 1.1 % at Δ/R₀ = 13.7 % and 1.3 % at 7.4 % for a
#   shift Δ perpendicular to the beam toward the equator.
# - Murakami & Nishi, Matter Radiat. Extremes 2, 55 (2017), Eqs. 28–37: ring-averaged
#   absorbed intensity of parallel-ray super-Gaussian beams whose axes cross the polar
#   axis at z R₀, re-solved here by direct quadrature.
# The ring references below integrate the absorbed intensity η_a(ψ) I_L cos ψ over the
# lit azimuths of each latitude (Gauss–Legendre, split where the lit range changes), so
# they converge spectrally; they share no code with Briareus.

# Gauss–Legendre nodes and weights on [-1, 1] (Newton on the three-term recurrence).
function _dd_gauss(n)
    x, w = zeros(n), zeros(n)
    for i in 1:n
        z = cos(π * (i - 0.25) / (n + 0.5))
        dp = 0.0
        for _ in 1:100
            p0, p1 = 1.0, z
            for k in 2:n
                p0, p1 = p1, ((2k - 1) * z * p1 - (k - 1) * p0) / k
            end
            dp = n * (z * p1 - p0) / (z^2 - 1)
            δ = p1 / dp
            z -= δ
            abs(δ) < 1e-15 && break
        end
        x[i], w[i] = z, 2 / ((1 - z^2) * dp^2)
    end
    return x, w
end

function _dd_P(l, μ)
    l == 0 && return one(μ)
    p0, p1 = one(μ), μ
    for n in 1:l-1
        p0, p1 = p1, ((2n + 1) * μ * p1 - n * p0) / (n + 1)
    end
    return p1
end

# One ring (continuous in azimuth) of beams from polar angle j: toward-source direction
# U = (sin j, 0, cos j) at φ = 0, axis through the point p (units of R₀), super-Gaussian
# exp(-[(u/a)² + (v/b)²]^(m/2)) with u azimuthal and v meridional (units of R₀), normal
# absorption η (MN17 Eq. 29; η = 1: everything incident is absorbed). Returns the
# ring-averaged absorbed intensity at latitude μ = cos θ (MN17 Eqs. 28–34).
function _dd_ring(μ, j, p, a, b, m, η, xg, wg)
    st, ct = sqrt(max(0.0, 1 - μ^2)), μ
    A, B = st * sin(j), ct * cos(j)
    B <= -A && return 0.0
    φ0 = B >= A ? π : acos(clamp(-B / A, -1.0, 1.0))
    s = 0.0
    for (x, w) in zip(xg, wg)
        φ = φ0 * x
        r = (st * cos(φ), st * sin(φ), ct)
        c = r[1] * sin(j) + r[3] * cos(j)
        c <= 0 && continue
        d = r .- p
        u = d[2]
        v = d[1] * cos(j) - d[3] * sin(j)
        ηa = η >= 1 ? 1.0 : 1 - (1 - η)^(c^3)
        s += w * ηa * exp(-((u / a)^2 + (v / b)^2)^(m / 2)) * c
    end
    return s * φ0 / (2π)
end

# ∫ f(μ) dμ over [-1, 1] for a ring field, split at the latitudes ±sin j where the lit
# azimuth range changes.
function _dd_integrate(f, j, xg, wg)
    brk = [-1.0, -sin(j), sin(j), 1.0]
    s = 0.0
    for k in 1:3
        lo, hi = brk[k], brk[k+1]
        for (x, w) in zip(xg, wg)
            s += w * (hi - lo) / 2 * f((lo + hi) / 2 + (hi - lo) / 2 * x)
        end
    end
    return s
end

# Legendre coefficients b_l/b_0 (l = 0:lmax) of one ring (MN17 Eqs. 35–36).
function _dd_ring_legendre(j, p, a, b, m, η; lmax = 30, n = 120)
    xg, wg = _dd_gauss(n)
    I = Dict{Float64,Float64}()
    ring(μ) = get!(() -> _dd_ring(μ, j, p, a, b, m, η, xg, wg), I, μ)
    bl = [(2l + 1) / 2 * _dd_integrate(μ -> ring(μ) * _dd_P(l, μ), j, xg, wg) for l in 0:lmax]
    return bl ./ bl[1]
end

# σ (rms / mean over the sphere, all degrees) of a ring and its mirror image in the
# equator: I(μ) + I(-μ).
function _dd_mirrored_sigma(j, p, a, b, m; n = 100)
    xg, wg = _dd_gauss(n)
    tot(μ) = _dd_ring(μ, j, p, a, b, m, 1.0, xg, wg) + _dd_ring(-μ, j, p, a, b, m, 1.0, xg, wg)
    vals = Dict{Float64,Float64}()
    f(μ) = get!(() -> tot(μ), vals, μ)
    m1 = _dd_integrate(f, j, xg, wg) / 2
    m2 = _dd_integrate(μ -> f(μ)^2, j, xg, wg) / 2
    return sqrt(m2 / m1^2 - 1)
end

# Golden-section minimum of a unimodal f on [lo, hi].
function _dd_golden(f, lo, hi; tol = 1e-5)
    g = (sqrt(5) - 1) / 2
    a, b = lo, hi
    c, d = b - g * (b - a), a + g * (b - a)
    fc, fd = f(c), f(d)
    while b - a > tol
        if fc < fd
            b, d, fd = d, c, fc
            c = b - g * (b - a); fc = f(c)
        else
            a, c, fc = c, d, fd
            d = a + g * (b - a); fd = f(d)
        end
    end
    x = (a + b) / 2
    return x, f(x)
end

# A bare sphere of radius R about the origin (the hard-sphere direct-drive target).
_dd_sphere(R, n) = Target(Component(:capsule, icosphere_mesh(R; n), Material(:CH);
                                    surface = SphereSurface(R)))

# Real-harmonic coefficients (l ≤ lmax) of the absorbed intensity of a deposition.
_dd_coeffs(W, t, dep) = W * (dep.element ./ t.elem_area)

@testset "hard-sphere direct drive" begin
    # Two meshes: the element (facet-constant) representation smooths the illumination at
    # O(h²), so references are matched on the finer one and the error is shown to fall
    # as h², with the Richardson extrapolation of the pair closer still.
    R0, L = 0.5, 30
    meshes = map((20, 30)) do n                # 8000 / 18000 elements, edges 28 / 18 µm
        tn = _dd_sphere(R0, n)
        (; n, t = tn, W = ylm_weights(tn, collect(1:nelements(tn)); lmax = L))
    end
    t, W = meshes[2].t, meshes[2].W
    els = collect(1:nelements(t))
    q = (30 / 20)^2
    richardson(x20, x30) = (q .* x30 .- x20) ./ (q - 1)

    @testset "Skupsky–Lee separability (OMEGA 60, SG4)" begin
        sp = omega_spot(:SG4)
        beams = omega_beams(; cones = :all, spot = sp, power = 1.0)
        U = [b.port for b in beams]
        G = [sqrt(abs(sum(_dd_P(l, dot(u, v)) for u in U, v in U))) / length(U) for l in 0:L]
        # single centred beam: I(θ) = I_L(R₀ sin θ) cos θ on the lit hemisphere
        xg, wg = _dd_gauss(2000)
        a = zeros(L + 1)
        for (x, w) in zip(xg, wg)
            θ = π / 4 * (1 + x)
            I = exp(-((R0 * sin(θ) / sp.a)^2)^(sp.order / 2)) * cos(θ) * sin(θ)
            for l in 0:L
                a[l+1] += (2l + 1) / 2 * w * π / 4 * I * _dd_P(l, cos(θ))
            end
        end
        σref = [abs(a[l+1] / a[1]) * G[l+1] / sqrt(2l + 1) for l in 0:L]
        rmsref = sqrt(sum(abs2, σref[2:end]))

        runs = map(meshes) do m
            d = deposit(DepositionOperator(m.t, beams; ncone = 1), beams)
            cc = _dd_coeffs(m.W, m.t, d)
            (; dep = d, c = cc, σ = mode_power(ylm_matrix(cc, L)))
        end
        dep, c, σ = runs[2]
        err20, err30 = (maximum(abs.(r.σ .- σref)) for r in runs)
        rmsx = richardson(mode_rms(runs[1].c, L), mode_rms(c, L))
        # absorbed power: the beam inside the sphere's disk (enclosed fraction)
        @test isapprox(sum(dep.element) / dep.delivered, enclosed_fraction(sp, R0 / sp.a); rtol = 1e-4)
        @test err30 < 2e-5                     # σₗ, absolute in units of the mean
        @test 1.8 < err20 / err30 < 2.8        # O(h²): (30/20)² = 2.25
        @test isapprox(mode_rms(c, L), rmsref; rtol = 0.025)
        @test isapprox(rmsx, rmsref; rtol = 0.005)
        # identities of the spectrum and its m = 0 / m ≠ 0 split
        @test mode_rms(c, L) ≈ sqrt(sum(abs2, σ[2:end])) rtol = 1e-12
        @test mode_rms(c, L; m = :zero)^2 + mode_rms(c, L; m = :nonzero)^2 ≈ mode_rms(c, L)^2 rtol = 1e-12
        C = ylm_matrix(c, L)
        @test mode_power(C; m = :zero) .^ 2 .+ mode_power(C; m = :nonzero) .^ 2 ≈ σ .^ 2 rtol = 1e-12
        @test mode_rms(c, 4) ≈ sqrt(sum(relative_mode_power(c, l) for l in 1:4)) rtol = 1e-12
        # the zonal spectrum is the Legendre spectrum of `legendre_moments`, an independent
        # facet quadrature (level 3 against level 2): it bounds the quadrature error at l ≤ 30
        al = legendre_moments(dep.element ./ t.elem_area, t, els; lmax = L)
        zl = [abs(al[l+1] / al[1]) / sqrt(2l + 1) for l in 0:L]
        @test all(abs.(mode_power(C; m = :zero) .- zl) .< 2e-5)
    end

    @testset "Temporal 2014 single ring" begin
        j = deg2rad(49)
        perp(Δ) = Δ .* (cos(j), 0.0, -sin(j))   # toward the equator, normal to the beam
        circ = (1.0, 1.0)                     # (azimuthal a, meridional b) / R₀
        ell = (2.0, 1.0)
        σref(Δ, (a, b)) = _dd_mirrored_sigma(j, perp(Δ), a, b, 4.0)
        # the reference reproduces the published values (conventions pinned: incident
        # I cos ψ fully absorbed, shift normal to the beam axis, a azimuthal)
        @test round(100σref(0.0, circ); digits = 1) == 11.9
        @test round(100σref(0.0, ell); digits = 1) == 5.2
        Δc, σc = _dd_golden(Δ -> σref(Δ, circ), 0.10, 0.17)
        Δe, σe = _dd_golden(Δ -> σref(Δ, ell), 0.04, 0.11)
        @test abs(Δc - 0.137) < 0.0015 && round(100σc; digits = 1) == 1.1
        @test abs(Δe - 0.074) < 0.0015 && round(100σe; digits = 1) == 1.3

        # Briareus: 8 beams per ring, mirrored; the zonal (m = 0) part of their pattern is
        # the continuous ring's azimuthal average exactly, whatever the number of beams
        N = 8
        ring = [Beam("N$k"; theta = 49, phi = 45k) for k in 0:N-1]
        ring = vcat(ring, [Beam("S$k"; theta = 131, phi = 45k) for k in 0:N-1])
        op = DepositionOperator(t, ring; ncone = 1)
        function σbri(Δ, (a, b))
            spot = SuperGaussianSpot(b * R0, a * R0; order = 4)   # a axis along θ̂ (meridional)
            beams = [remake(bm; spot, aim = Δ * R0 * sign(bm.port[3]) * Briareus.port_basis(bm.port)[1])
                     for bm in ring]
            return mode_rms(_dd_coeffs(W, t, deposit(op, beams)), L; m = :zero)
        end
        cases = [(0.0, circ), (0.0, ell), (Δc, circ), (Δe, ell)]
        got = [σbri(Δ, s) for (Δ, s) in cases]
        want = [σref(Δ, s) for (Δ, s) in cases]
        @test all(abs.(got .- want) .< 1e-4)

        # the metric is differentiable in the pointing (ForwardDiff through `deposit`)
        g = ForwardDiff.derivative(Δ -> σbri(Δ, circ), 0.05)
        h = 1e-4
        fd = (σbri(0.05 + h, circ) - σbri(0.05 - h, circ)) / 2h
        @test isapprox(g, fd; rtol = 1e-6)
    end

    @testset "Murakami & Nishi 2017 two-cone PDD" begin
        # MN17 §4: two cones per hemisphere (j, z, E), a = 0.6, b = 6, η⊥ = 0.9; printed
        # solution j = 46°, 53°, z = -0.65, 0.80, E₂ = 0.40 for B₂ = … = B₁₀ = 0 (Eq. 37,
        # B = Σ Eᵢ bₗ(jᵢ, zᵢ) with bₗ of each ring normalized to its mean).
        a, b, η = 0.6, 6.0, 0.9
        function resid(x; n = 100)
            j1, j2, z1, z2, E2 = x
            b1 = _dd_ring_legendre(deg2rad(j1), (0.0, 0.0, z1), a, a, b, η; lmax = 10, n)
            b2 = _dd_ring_legendre(deg2rad(j2), (0.0, 0.0, z2), a, a, b, η; lmax = 10, n)
            return [b1[l+1] + E2 * b2[l+1] for l in 2:2:10]
        end
        printed = [46.0, 53.0, -0.65, 0.80, 0.40]
        r0 = resid(printed)
        # Newton from the printed values (central-difference Jacobian)
        x = copy(printed)
        for _ in 1:8
            r = resid(x)
            norm(r) < 1e-11 && break
            h = [1e-4, 1e-4, 1e-6, 1e-6, 1e-6]
            J = reduce(hcat, [(resid(x .+ h[k] .* (1:5 .== k)) .- resid(x .- h[k] .* (1:5 .== k))) ./ 2h[k]
                              for k in 1:5])
            x = x .- J \ r
        end
        @test norm(resid(x)) < 1e-10
        # the printed digits do not zero the modes (|B₄| ≈ 0.19); the nearby root found
        # here is j = 48.93°, 51.69°, z = -0.576, 0.769, E₂ = 0.421
        @test abs(r0[2]) > 0.1
        @test isapprox(x, [48.926, 51.685, -0.5756, 0.7687, 0.4208]; atol = 2e-3)

        # Briareus on each ring of the re-solved design, η⊥ = 1 (no absorption model):
        # 12 parallel-ray beams aimed through (0, 0, z R₀), spot 1/e radius 0.6 R₀, order 6
        rings = [(x[1], x[3]), (x[2], x[4])]
        errs = map(rings) do (jd, z)
            beams = [Beam("R$k"; theta = jd, phi = 30k, spot = SuperGaussianSpot(a * R0; order = b),
                          aim = (0.0, 0.0, z * R0)) for k in 0:11]
            bl = map(meshes) do m
                cc = _dd_coeffs(m.W, m.t, deposit(DepositionOperator(m.t, beams; ncone = 1), beams))
                [sqrt(2l + 1) * cc[mode_index(l, 0)] / cc[1] for l in 0:L]
            end
            ref = _dd_ring_legendre(deg2rad(jd), (0.0, 0.0, z), a, a, b, 1.0; lmax = L)
            (maximum(abs.(bl[1] .- ref)), maximum(abs.(bl[2] .- ref)),
             maximum(abs.(richardson(bl...) .- ref)))
        end
        @test all(e -> e[2] < 3e-3, errs)
        @test all(e -> 1.8 < e[1] / e[2] < 2.8, errs)     # O(h²)
        @test all(e -> e[3] < 2e-4, errs)
    end
end
