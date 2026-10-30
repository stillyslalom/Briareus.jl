# Synthetic x-ray pinhole imagers (`src/imager.jl`, kind `xray_imager`): the detected
# radiance against a direct quadrature of the response over a Planckian, the image's
# total against the emitter's area, the point-spread function's variance against the
# disk's and the Gaussian's closed forms, gate integrals against closed forms, and the
# built-in definitions.

using Briareus, Test

@testset "x-ray imagers" begin
    R, T = 2.0, 200.0
    disk = Component(:emitter, disk_mesh(R, 0.0; nr = 4, nphi = 32), gold())
    t = Target(disk)
    m = nelements(t)
    Bd = SIGMA_SB * T^4
    sol = RadiationSolution(RadiationModel(t; F = zeros(m, m)), zeros(m), zeros(m), fill(Bd, m),
                            zeros(m), zeros(m), zeros(m), zeros(m), 0.0)
    filters = [:Be => 0.075, :Ti => 0.050]
    sxi(; kw...) = XrayImager(0, 0; pinhole = 0.025, distance = 1656.7, magnification = 2.0, width = 6.0,
                              pixels = 300, filters, detector = (:Si, 0.015), detector_pixel = 0.024,
                              energy = 4000, kw...)
    im = sxi()

    @testset "detected radiance and the image's total" begin
        # the response over a 200 eV Planckian by a direct trapezoid on its own grid
        E = collect(range(50.0, 29_990.0; length = 100_001))
        Rr = [filter_transmission(filters, e / 1000; table = :henke) *
              (1 - filter_transmission([:Si => 0.015], e / 1000; table = :henke)) for e in E]
        P = [e^3 / expm1(e / T) for e in E]
        mR = sum(Rr .* P) / sum(P)
        L = detected_radiance(im, sol)
        @test all(l -> isapprox(l, Bd / π * mR; rtol = 2e-3), L)
        # the blur conserves the total of an object inside the field; the unblurred
        # image covers the 32-gon's area to pixel resolution
        h = 6.0 / 300
        raw = [p > 0 ? L[p] : 0.0 for p in Briareus._imager_pixels(im, t)]
        img = imager_image(im, sol)
        @test sum(img) ≈ sum(raw) rtol = 1e-12
        @test count(>(0), raw) * h^2 ≈ 16R^2 * sin(2π / 32) rtol = 1e-2
        @test imager_image(sxi(center = :emitter), sol) ≈ img rtol = 1e-12   # the disk is at the origin
        @test pinhole_solid_angle(im) ≈ π * 0.025^2 / 4 / 1656.7^2 rtol = 1e-14
        # one group over all energies with the brightness temperature: the gray result
        s = SpectralSolution(sol, [0.0, Inf], zeros(m, 1), fill(Bd, m, 1), zeros(m, 1), zeros(m, 1),
                             fill(T, m), zeros(m), [0.0], false)
        @test detected_radiance(im, s) ≈ L rtol = 1e-12
    end

    @testset "point-spread function" begin
        # the penumbra 37.5 µm, diffraction 1.03 λ distance / d at 4 keV (Schneider 2010's
        # SXI: "diffraction and geometry about equal")
        @test im.penumbra ≈ 0.025 * 1.5 rtol = 1e-14
        @test im.diffraction ≈ 1.03 * 1.239841984e-6 / 4 * 1656.7 / 0.025 rtol = 1e-12
        # the kernel's variance per axis: the disk's D²/16 plus the Gaussian's σ², on a
        # fine sampling
        fine = sxi(width = 0.4, pixels = 400)                  # 1 µm pixels
        k = Briareus._psf(fine)
        r = size(k, 1) ÷ 2
        h = 1e-3
        var = sum(k[i+r+1, j+r+1] * (i * h)^2 for i in -r:r, j in -r:r)
        σ² = (fine.diffraction / 2.3548)^2 + (0.024 / 2)^2 / 12
        @test var ≈ fine.penumbra^2 / 16 + σ² rtol = 1e-2
        @test sum(k) ≈ 1 rtol = 1e-12
    end

    @testset "recording: integrated and gated" begin
        # a signal linear in time at every element: integrals in closed form
        time = collect(0.0:0.5:10.0)
        L0 = detected_radiance(im, sol)
        Ls = L0 * (1 .+ 0.3 .* time')
        integ = imager_record(im, t, time, Ls)
        @test isempty(integ.times) && length(integ.images) == 1
        @test sum(integ.images[1]) ≈ sum(imager_image(im, t, L0)) * (10 + 0.3 * 50) rtol = 1e-12
        # a Gaussian gate of FWHM w inside the history: (1 + 0.3 t0) w √(π / 4 ln 2)
        w, t0 = 1.0, 4.2
        gated = imager_record(sxi(mode = :gated, times = [t0, 6.0], gate = w), t, time, Ls)
        @test gated.times == [t0, 6.0] && length(gated.images) == 2
        @test sum(gated.images[1]) ≈ sum(imager_image(im, t, L0)) * (1 + 0.3t0) * w * sqrt(π / (4log(2))) rtol = 1e-6
        sq = imager_record(sxi(mode = :gated, times = [t0], gate = w, gate_shape = :square), t, time, Ls)
        @test sum(sq.images[1]) ≈ sum(imager_image(im, t, L0)) * (1 + 0.3t0) * w rtol = 1e-3
        @test_throws ArgumentError sxi(mode = :gated)
        @test_throws ArgumentError imager_record(im, t, time, Ls[:, 1:3])
    end

    @testset "definitions" begin
        lo, up = diagnostic("nif_sxi_lower"), diagnostic("nif_sxi_upper")
        @test lo isa XrayImager && lo.mode === :integrated && lo.energy == 4000
        @test isapprox(lo.camera.width, 2048 * 0.024 / 2; rtol = 1e-14) && lo.camera.nx == 1024
        @test lo.camera.dir ≈ -port_vector(161, 326.25) rtol = 1e-14
        @test up.distance == 1626.7 && up.theta == 18 && up.R == lo.R
        issues = filter(i -> basename(i.file) == "nif_gleh2.toml", check_diagnostics())
        @test Set(i.key for i in issues if i.level === :missing) ⊇ Set(["filters"])
        # a group's G-LEH-2: what is missing supplied
        mktempdir() do dir
            write(joinpath(dir, "g.toml"), """
                schema = 1
                id = "group_gleh2"
                extends = "nif_gleh2"
                filters = [{ material = "Ti", thickness_um = 10 }, { material = "polyimide", thickness_um = 0.5 },
                           { material = "Al", thickness_um = 1 }, { material = "lexan", thickness_um = 0.5 }]
                [detector]
                absorber = "Si"
                absorber_um = 10
                [recording]
                times_ns = [6.0, 8.0, 10.0, 12.0]
                gate_ns = 1.0
                """)
            g = diagnostic("group_gleh2"; paths = [diagnostic_paths()[1], dir])
            @test g.mode === :gated && g.times == [6.0, 8.0, 10.0, 12.0] && g.magnification == 1
            @test isapprox(g.camera.width, 1024 * 0.025; rtol = 1e-14)
            @test !any(i -> occursin(dir, i.file) && i.level !== :warning,
                       check_diagnostics(; paths = [diagnostic_paths()[1], dir]))
        end
    end
end
