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
# - Temporal et al., HPLSE 2, e8 (2014), §4: LMJ option A (20 quads at 49°/131°, the
#   hemispheres rotated 18°), elliptical m = 4 quads Δa = 2Δb = 640 µm; minimum over r0
#   of the intrinsic σ0 = 4.6 % (a 2-D quadrature of all 20 beams here).
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

# Temporal et al. 2014 (HPLSE) LMJ option A: 10 quads per ring at 49° and 131°, the
# rings 36° apart in azimuth and the hemispheres rotated 18° (uniform spacing assumed),
# super-Gaussian quads exp(-[(u/Δa)² + (v/Δb)²]^(m/2)) with v meridional, on a sphere of
# radius r0 (all lengths in units of Δb). Returns the toward-source directions with
# their (θ̂, φ̂) frames, and σ0 (rms/mean over the sphere, all modes) of the projected
# intensity I cos ψ fully absorbed, by Gauss–Legendre in μ and the trapezoid rule in φ.
function _dd_lmj_ports()
    pv(θ, φ) = (sind(θ) * cosd(φ), sind(θ) * sind(φ), cosd(θ))
    [(θ, φ, pv(θ, φ), pv(θ + 90, φ), (-sind(φ), cosd(φ), 0.0))
     for k in 0:9 for (θ, φ) in ((49.0, 36k), (131.0, 36k + 18))]
end
function _dd_lmj_sigma(r0, a, b, m; nμ = 200, nφ = 360)
    ports = _dd_lmj_ports()
    xg, wg = _dd_gauss(nμ)
    s1 = s2 = 0.0
    for (μ, wμ) in zip(xg, wg)
        st = sqrt(1 - μ^2)
        for j in 0:nφ-1
            s = (st * cospi(2j / nφ), st * sinpi(2j / nφ), μ)
            I = 0.0
            for (_, _, U, eθ, eφ) in ports
                c = sum(s .* U)
                c > 0 || continue
                u, v = r0 * sum(s .* eφ), r0 * sum(s .* eθ)
                I += exp(-((u / a)^2 + (v / b)^2)^(m / 2)) * c
            end
            s1 += wμ * I / nφ
            s2 += wμ * I^2 / nφ
        end
    end
    return sqrt(2 * s2 / s1^2 - 1)
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
            o = DepositionOperator(m.t, beams; ncone = 1)
            d = deposit(o, beams)
            cc = _dd_coeffs(m.W, m.t, d)
            (; dep = d, c = cc, σ = mode_power(ylm_matrix(cc, L)), op = o)
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

        @testset "absorption weight η_a (MN17 Eq. 29)" begin
            # off is the default, bit for bit (and unchanged from before the option: see the
            # research notes); η⊥ = 1 (η_a ≡ 1) is the same model, to rounding
            op = runs[2].op
            dn = deposit(op, beams; absorption = nothing)
            @test dn.element == dep.element && dn.per_beam == dep.per_beam
            d1 = deposit(op, beams; absorption = 1.0)
            @test isapprox(d1.element, dep.element; rtol = 1e-12)
            @test_throws ArgumentError deposit(op, beams; absorption = 0.0)
            @test_throws ArgumentError deposit(op, beams; absorption = 1.5)
            # a centred beam's ray at polar angle θ from its axis meets the sphere at
            # ψ = θ, so the weighted single-beam pattern stays axisymmetric and
            # separability holds with aₗ of I_L(R₀ sin θ) η_a(θ) cos θ (1-D quadrature)
            η = 0.6
            ηa(μ) = 1 - (1 - η)^(μ^3)
            aη = zeros(L + 1)
            num = den = 0.0
            for (x, w) in zip(xg, wg)
                θ = π / 4 * (1 + x)
                I = exp(-((R0 * sin(θ) / sp.a)^2)^(sp.order / 2)) * cos(θ) * sin(θ)
                num += w * I * ηa(cos(θ))
                den += w * I
                for l in 0:L
                    aη[l+1] += (2l + 1) / 2 * w * π / 4 * I * ηa(cos(θ)) * _dd_P(l, cos(θ))
                end
            end
            σηref = [abs(aη[l+1] / aη[1]) * G[l+1] / sqrt(2l + 1) for l in 0:L]
            rη = map(runs, meshes) do r, m
                d = deposit(r.op, beams; absorption = η)
                (; dep = d, σ = mode_power(ylm_matrix(_dd_coeffs(m.W, m.t, d), L)))
            end
            # σₗ: O(h²) as without the weight (the facets' flat normals now also enter
            # η_a(ψ), so the error at n = 30 is ~3× the unweighted one)
            e20, e30 = (maximum(abs.(r.σ .- σηref)) for r in rη)
            ex = maximum(abs.(richardson(rη[1].σ, rη[2].σ) .- σηref))
            @test e30 < 1e-4 && ex < 1e-5
            @test 1.8 < e20 / e30 < 2.8
            # the weight changes the pattern by far more than the error
            @test maximum(abs.(σηref .- σref)) > 50 * e30
            # absorbed power: the enclosed fraction times the η_a-weighted share
            fη = enclosed_fraction(sp, R0 / sp.a) * num / den
            f20, f30 = (sum(r.dep.element) / r.dep.delivered for r in rη)
            @test isapprox(f30, fη; rtol = 5e-4) && isapprox(richardson(f20, f30), fη; rtol = 2e-5)
            @test sum(rη[2].dep.per_beam) ≈ sum(rη[2].dep.element) rtol = 1e-12
        end
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
        function σbri(Δ, (a, b); absorption = nothing)
            spot = SuperGaussianSpot(b * R0, a * R0; order = 4)   # a axis along θ̂ (meridional)
            beams = [remake(bm; spot, aim = Δ * R0 * sign(bm.port[3]) * Briareus.port_basis(bm.port)[1])
                     for bm in ring]
            return mode_rms(_dd_coeffs(W, t, deposit(op, beams; absorption)), L; m = :zero)
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
        # and so is the absorption-weighted one, in the pointing and in η⊥ itself
        g = ForwardDiff.derivative(Δ -> σbri(Δ, circ; absorption = 0.9), 0.05)
        fd = (σbri(0.05 + h, circ; absorption = 0.9) - σbri(0.05 - h, circ; absorption = 0.9)) / 2h
        @test isapprox(g, fd; rtol = 1e-6)
        g = ForwardDiff.derivative(η -> σbri(0.05, circ; absorption = η), 0.9)
        fd = (σbri(0.05, circ; absorption = 0.9 + h) - σbri(0.05, circ; absorption = 0.9 - h)) / 2h
        @test isapprox(g, fd; rtol = 1e-6)

        @testset "spherical pointing fields" begin
            # The same ring pointed on the sphere (`GroupPointing`'s `:theta_s`): an axis
            # shifted by Δ R₀ toward the equator, normal to the beam, first meets the
            # sphere asin(Δ) further from the pole, so θ_s = 49° + asin(Δ) is the same
            # axis, and collimated beams (`ncone = 1`) deposit the same pattern.
            spot = SuperGaussianSpot(R0, R0; order = 4)
            ringc = [remake(bm; spot) for bm in ring]
            p = GroupPointing(ringc; fields = (:theta_s, :dphi_s), radius = R0)
            @test initial_parameters(p) ≈ [49.0, 0.0] atol = 1e-12   # centred template
            σs(x) = mode_rms(_dd_coeffs(W, t, deposit(op, apply(p, x))), L; m = :zero)
            θc = 49 + asind(Δc)
            @test isapprox(σs([θc, 0.0]), σbri(Δc, circ); rtol = 1e-10)
            @test abs(σs([θc, 0.0]) - σref(Δc, circ)) < 1e-4
            # dσ/dθ_s and dσ/dΔφ_s (all modes) by ForwardDiff against central differences
            σa(x) = mode_rms(_dd_coeffs(W, t, deposit(op, apply(p, x))), L)
            x = [55.0, 3.0]
            gs = ForwardDiff.gradient(σa, x)
            hs = 1e-3
            fds = [(σa(x .+ hs .* (1:2 .== k)) - σa(x .- hs .* (1:2 .== k))) / 2hs for k in 1:2]
            @test isapprox(gs, fds; rtol = 1e-6)
        end
    end

    @testset "Temporal 2014 LMJ intrinsic σ0" begin
        # HPLSE 2, e8 (2014), §4 and Fig. 5: LMJ option A (the 20 quads at 49°/131°),
        # elliptical quads Δb = 320 µm meridional, Δa = 2Δb, m = 4, no beam errors and no
        # PDD shift; the minimum over the capsule radius r0 of the intrinsic σ0 is 4.6 %
        # (its curve's minimum lies near r0 = 350 µm in Fig. 5; the r0 = 320 µm of the
        # text is the optimum of the error-averaged σ̄). Same convention as Te14.
        σref(ρ) = _dd_lmj_sigma(ρ, 2.0, 1.0, 4.0)          # ρ = r0/Δb
        ρmin, σmin = _dd_golden(σref, 1.0, 1.25; tol = 1e-4)
        @test round(100σmin; digits = 1) == 4.6
        @test 340 < 320ρmin < 360
        # the quadrature is converged far below the published digit
        @test abs(_dd_lmj_sigma(ρmin, 2.0, 1.0, 4.0; nμ = 300, nφ = 540) - σmin) < 1e-6

        # Briareus: the same 20 collimated quads on the n = 30 sphere (σ to l = 30)
        quadsL = [Beam("L$k"; theta = θ, phi = φ) for (k, (θ, φ)) in enumerate(_dd_lmj_ports())]
        opL = DepositionOperator(t, quadsL; ncone = 1)
        σL(ρ) = mode_rms(_dd_coeffs(W, t, deposit(opL, [remake(q; spot = SuperGaussianSpot(R0 / ρ, 2R0 / ρ; order = 4))
                                                         for q in quadsL])), L)
        @test abs(σL(ρmin) - σmin) < 1e-4 && abs(σL(1.0) - σref(1.0)) < 1e-4
        @test round(100 * σL(ρmin); digits = 1) == 4.6
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

        # Briareus on each ring of the re-solved design, without absorption (η⊥ = 1) and
        # with MN17's η_a at η⊥ = 0.9 (`absorption`): 12 parallel-ray beams aimed through
        # (0, 0, z R₀), spot 1/e radius 0.6 R₀, order 6
        rings = [(x[1], x[3]), (x[2], x[4])]
        res = map(rings) do (jd, z)
            beams = [Beam("R$k"; theta = jd, phi = 30k, spot = SuperGaussianSpot(a * R0; order = b),
                          aim = (0.0, 0.0, z * R0)) for k in 0:11]
            ops = [DepositionOperator(m.t, beams; ncone = 1) for m in meshes]
            map((nothing, η)) do absorption
                bl = map(ops, meshes) do o, m
                    cc = _dd_coeffs(m.W, m.t, deposit(o, beams; absorption))
                    [sqrt(2l + 1) * cc[mode_index(l, 0)] / cc[1] for l in 0:L]
                end
                ref = _dd_ring_legendre(deg2rad(jd), (0.0, 0.0, z), a, a, b,
                                        something(absorption, 1.0); lmax = L)
                (; bl, err = (maximum(abs.(bl[1] .- ref)), maximum(abs.(bl[2] .- ref)),
                              maximum(abs.(richardson(bl...) .- ref))))
            end
        end
        errs = [r.err for rr in res for r in rr]
        @test all(e -> e[2] < 4e-3, errs)
        @test all(e -> 1.8 < e[1] / e[2] < 2.8, errs)     # O(h²)
        @test all(e -> e[3] < 2e-4, errs)
        # the absorption-weighted design: Briareus's B₂…B₁₀ (Eq. 37) vanish to the mesh error
        B(k) = [res[1][2].bl[k][l+1] + x[5] * res[2][2].bl[k][l+1] for l in 2:2:10]
        Bx = richardson(B(1), B(2))
        @test maximum(abs, B(2)) < 5e-4 && maximum(abs, Bx) < 1e-5
    end
end

@testset "spherical pointing (GroupPointing)" begin
    # Polar-direct-drive pointings are published as points on the target sphere (polar
    # angle, azimuth offset from the port; Barlow et al. 2024 Table I). The reference is
    # the aim set built by hand from those angles, with the mirrored lower hemisphere at
    # 180° − θ_s.
    R = 1.1
    quads = nif_quads(power = 4.0)
    cones = (:inner23, :inner30, :outer44, :outer50)
    θs = Dict(zip(cones, (6.36, 14.95, 52.48, 85.62)))     # Ba24 Table I, as an input set
    Δφs = Dict(zip(cones, (52.02, 160.75, -12.14, 6.90)))
    p = GroupPointing(quads; fields = (:theta_s, :dphi_s, :defocus), radius = R)
    x = reduce(vcat, [[θs[k], Δφs[k], 5.0] for k in p.keys])
    bs = apply(p, x)
    hand = map(quads) do q
        θp, φp = port_angles(q)
        θ = q.port[3] > 0 ? θs[q.group] : 180 - θs[q.group]
        R * port_vector(θ, φp + Δφs[q.group])
    end
    @test all(norm.([b.aim for b in bs] .- hand) .< 1e-14)
    @test all(b -> b.defocus == 5.0, bs) && all(b -> b.power == 4.0, bs)
    # every aim faces its port, so the axis first meets the sphere there
    @test all(b -> dot(b.aim, b.port) > 0, bs)
    # round trip: the angles read back from the new beams' axes
    @test isapprox(initial_parameters(GroupPointing(bs; fields = p.fields, radius = R)), x; atol = 1e-11)
    # relative shifts from that pointing move every aim by the same angles
    pr = GroupPointing(bs; fields = (:theta_s, :dphi_s), radius = R, relative = true)
    @test initial_parameters(pr) == zeros(8)
    @test all(norm.([b.aim for b in apply(pr, zeros(8))] .- hand) .< 1e-14)
    xr = repeat([2.0, -3.0], 4)
    hr = map(q -> R * port_vector(q.port[3] > 0 ? θs[q.group] + 2 : 178 - θs[q.group],
                                  port_angles(q)[2] + Δφs[q.group] - 3), quads)
    @test all(norm.([b.aim for b in apply(pr, xr)] .- hr) .< 1e-13)
    # a shifted centre
    c = (0.1, -0.2, 0.05)
    pc = GroupPointing(quads; fields = (:theta_s, :dphi_s), radius = R, center = c)
    xc = reduce(vcat, [[θs[k], Δφs[k]] for k in pc.keys])
    @test all(norm.([b.aim for b in apply(pc, xc)] .- (hand .+ Ref(Vec3(c...)))) .< 1e-14)

    # split quads on the 192 beams: each beam moves by ±split in θ and Δφ, the signs of
    # its port's offset from its quad's mean port (here from the beam table directly)
    beams = nif_beams()
    ps = GroupPointing(beams; fields = (:theta_s, :dphi_s, :split_theta, :split_phi), radius = R)
    xs = reduce(vcat, [[θs[k], Δφs[k], 1.5, 4.0] for k in ps.keys])
    wrap(a) = mod(a + 180, 360) - 180
    hs = map(beams) do b
        mates = filter(o -> o.quad == b.quad, beams)
        θ(o) = b.port[3] > 0 ? port_angles(o)[1] : 180 - port_angles(o)[1]
        sθ = sign(θ(b) - sum(θ, mates) / 4)
        sφ = sign(wrap(port_angles(b)[2] - port_angles(mates[1])[2]) -
                  sum(o -> wrap(port_angles(o)[2] - port_angles(mates[1])[2]), mates) / 4)
        θb = θs[b.group] + 1.5sθ
        R * port_vector(b.port[3] > 0 ? θb : 180 - θb, port_angles(b)[2] + Δφs[b.group] + 4sφ)
    end
    @test all(norm.([b.aim for b in apply(ps, xs)] .- hs) .< 1e-13)
    # a 2 × 2 quad: two beams on each side in both directions
    @test all(q -> count(i -> beams[i].quad == q && ps.sph[i][3] > 0, eachindex(beams)) == 2 &&
                   count(i -> beams[i].quad == q && ps.sph[i][4] > 0, eachindex(beams)) == 2,
              unique(b.quad for b in beams))
    # bounds in degrees
    lo, hi = bounds(ps, xs; theta_s = 2.0, split_phi = 1.0)
    @test hi[1:4] - lo[1:4] == [4.0, 20.0, 10.0, 2.0]

    @test_throws ArgumentError GroupPointing(quads; fields = (:theta_s, :z), radius = R)
    @test_throws ArgumentError GroupPointing(quads; fields = (:theta_s,))
    @test_throws ArgumentError GroupPointing(omega_beams(); fields = (:split_phi,), radius = R)
    far = [remake(quads[1]; aim = Vec3(5.0, 5.0, 5.0))]
    @test_throws ArgumentError GroupPointing(far; fields = (:theta_s,), radius = R)
end
