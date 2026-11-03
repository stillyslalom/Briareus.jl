# Soft x-ray optics against closed forms: Henke photoabsorption against NIST's tabulated
# μ/ρ at 1 keV (scattering < 0.1 % there), Beer–Lambert transmission from μ, Fresnel
# reflectivity at normal incidence and far above the critical angle, the Névot–Croce factor
# between its limits, and tabulated curves exact on linear data. (dante_response.jl checks
# Henke against NIST at 3 keV and the mirror cutoffs.)

@testset "x-ray optics against closed forms" begin
    # NIST X-Ray Mass Attenuation Coefficients (Hubbell & Seltzer), μ/ρ at 1 keV:
    # Al 1185, C 2211 cm²/g
    @test henke_mass_attenuation(:Al, 1000.0) ≈ 1185 rtol = 0.02
    @test henke_mass_attenuation(:C, 1000.0) ≈ 2211 rtol = 0.02
    # a compound is the mass-weighted mean of its elements (here Lexan, C16H14O3)
    w = (16 * 12.011, 14 * 1.008, 3 * 15.999)
    μ = sum(w .* henke_mass_attenuation.((:C, :H, :O), 2500.0)) / sum(w)
    @test henke_mass_attenuation(:lexan, 2500.0) ≈ μ rtol = 1e-12
    # Beer–Lambert through layers (mm; ρ in g/cm³): exp(−Σ μ/ρ ρ t)
    T = filter_transmission([:lexan => 0.01, :Al => (0.002, 2.7)], 2.5; table = :henke)
    @test T ≈ exp(-(henke_mass_attenuation(:lexan, 2500.0) * 1.20 * 0.001 +
                    henke_mass_attenuation(:Al, 2500.0) * 2.7 * 0.0002)) rtol = 1e-12
    @test_throws ArgumentError henke_mass_attenuation(:Al, 5.0)
    @test_throws ArgumentError henke_mass_attenuation(:CH4, 1000.0)

    # Fresnel: n² = 1 − χ* with χ = r_e λ² N f/π (f = f₁ + i f₂ summed over the formula)
    E = 8000.0
    λ = 1.239841984e-4 / E                               # cm
    f = complex(Briareus._henke_f(:Si, E)...) + 2 * complex(Briareus._henke_f(:O, E)...)
    N = 2.2 * 6.02214076e23 / (28.085 + 2 * 15.999)      # SiO₂ formula units per cm³
    χ = 2.8179403262e-13 * λ^2 * N * f / π
    # at normal incidence r = (1 − n)/(1 + n) ≈ χ*/4: R = |χ|²/16 to O(|χ|)
    @test mirror_reflectivity(:SiO2, 90.0, E) ≈ abs2(χ) / 16 rtol = 1e-4
    # far above the critical angle θc = √(Re χ): r_s ≈ χ*/(4 sin²θ) and |r_p| = |r_s| cos 2θ,
    # to O(θc²/θ²)
    θc = sqrt(real(χ))
    θ = 30θc
    @test mirror_reflectivity(:SiO2, rad2deg(θ), E) ≈ abs2(χ) / (16 * sin(θ)^4) * (1 + cos(2θ)^2) / 2 rtol = 0.005
    # total external reflection below it, bounded by 1
    @test 0.95 < mirror_reflectivity(:SiO2, rad2deg(θc / 4), E) < 1
    # roughness σ: the Névot–Croce factor lies between 1 and the Debye–Waller exp(−4k²σ² sin²θ),
    # and approaches the latter far above θc
    k = 2π / λ
    σ = 0.1                                              # nm
    dw(θ) = exp(-4 * k^2 * sin(θ)^2 * (1e-7σ)^2)
    ratio(θ) = mirror_reflectivity(:SiO2, rad2deg(θ), E; roughness = σ) / mirror_reflectivity(:SiO2, rad2deg(θ), E)
    @test dw(2θc) < ratio(2θc) < 1
    @test ratio(20θc) ≈ dw(20θc) rtol = 0.01
    # a density given overrides the table's: an Si mirror at half density has the χ of half N
    @test mirror_reflectivity(:Si, 90.0, E; density = 2.329 / 2) ≈ mirror_reflectivity(:Si, 90.0, E) / 4 rtol = 1e-4

    # tabulated curves: exact on linear data, held or zero outside
    c = TabulatedCurve([1.0, 2.0, 4.0], [3.0, 5.0, 1.0])
    xs = [1.0, 1.25, 2.0, 3.0, 3.5, 4.0]
    @test c.(xs) ≈ [3.0, 3.5, 5.0, 3.0, 2.0, 1.0] atol = 1e-14
    @test (c(0.5), c(9.0)) == (3.0, 1.0)
    z = TabulatedCurve([1.0, 2.0, 4.0], [3.0, 5.0, 1.0]; extrapolate = :zero)
    @test (z(0.5), z(9.0), z(1.0)) == (0.0, 0.0, 3.0)
    @test_throws ArgumentError TabulatedCurve([1.0, 1.0], [0.0, 1.0])
    @test_throws ArgumentError TabulatedCurve([1.0], [0.0])
    @test_throws ArgumentError TabulatedCurve([1.0, 2.0], [0.0, NaN])
    @test_throws ArgumentError TabulatedCurve([1.0, 2.0], [0.0, 1.0]; extrapolate = :linear)
    # read from a two-column CSV with a header line
    path = tempname() * ".csv"
    write(path, "energy_eV,transmission\n100,0.1\n200,0.5\n400,0.9\n")
    cf = TabulatedCurve(path)
    @test cf.E == [100, 200, 400] && isapprox(cf(300.0), 0.7; rtol = 1e-14) && cf.file == path
    rm(path)
end
