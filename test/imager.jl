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
        # the response over a 200 eV Planckian by a direct trapezoid on its grid
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
        # The penumbra is 37.5 µm and the diffraction 1.03 λ distance / d at 4 keV. For the
        # SXI, Schneider 2010 gives "diffraction and geometry about equal".
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

    @testset "closing lips on the line of sight" begin
        # An emitting disk seen from the pole through a closing lip, an annular aperture
        # of two rings with open fractions 0.25 and 0.6. The power toward the detector is
        # the radiance times the emitter area less the closed area of the lip rings,
        # (1 − τ) of each lip element. The rendered image, the line-of-sight power and the
        # Dante view take the same open fractions, and the imager scales its image by the
        # same ratio.
        lip = aperture(:lip, annulus_mesh(0.8, 1.5, 1.0; nr = 2, nphi = 32);
                       closure = GapClosure(closure = 5.0))
        tl = Target([disk, lip])
        n = nelements(tl)
        le = elements(tl, :lip)
        τ = zeros(n)
        rmid(e) = hypot(tl.elem_centroid[e][1], tl.elem_centroid[e][2])
        split = sum(rmid, le) / length(le)
        τ[le] .= [rmid(e) < split ? 0.25 : 0.6 for e in le]
        B = [e in le ? 0.0 : Bd for e in 1:n]
        sl = RadiationSolution(RadiationModel(tl; F = zeros(n, n)), zeros(n), zeros(n), B,
                               zeros(n), zeros(n), zeros(n), τ, 0.0)
        emit = sum(tl.elem_area[elements(tl, :emitter)])
        closed = sum((1 - τ[e]) * tl.elem_area[e] for e in le)
        exact = Bd / π * (emit - closed)
        @test 0.05 < closed / emit < 0.5
        cam = OrthoCamera((0, 0, -1); up = (0, 1, 0), width = 6.0, nx = 600)
        px = (6.0 / 600)^2
        @test sum(render(tl, radiance(sl), cam; background = 0.0, sl.transmission)) * px ≈ exact rtol = 5e-3
        @test sum(render(tl, radiance(sl), cam; background = 0.0)) * px ≈ Bd / π * emit rtol = 5e-3
        @test line_of_sight_power(tl, sl.B, (0, 0, 1); sl.transmission, nsub = 64, refine = 1) ≈ exact rtol = 5e-3
        @test line_of_sight_power(tl, sl.B, (0, 0, 1)) ≈ Bd / π * emit rtol = 1e-12
        view = DanteView(0, 0; through = nothing, field_of_view = 6.0, pixels = 600)
        @test dante(view, sl).intensity ≈ exact rtol = 5e-3
        Ll = detected_radiance(im, sl)
        @test sum(imager_image(im, sl)) / sum(imager_image(im, tl, Ll)) ≈ (emit - closed) / emit rtol = 5e-3
        @test imager_image(im, tl, Ll; sl.transmission) == imager_image(im, sl)
        # a record of two times with the lip closing: the integral of L τ(t), linear in t
        τ2 = hcat(ones(n), τ)
        rec = imager_record(sxi(mode = :integrated), tl, [0.0, 2.0], hcat(Ll, Ll); transmission = τ2)
        both = imager_image(im, tl, Ll) .+ imager_image(im, sl)
        @test rec.images[1] ≈ both rtol = 1e-12
        @test_throws ArgumentError imager_record(sxi(mode = :integrated), tl, [0.0, 2.0], hcat(Ll, Ll);
                                                 transmission = ones(n, 3))
    end
end
