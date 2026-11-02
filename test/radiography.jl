# Backlit radiography (`src/radiography.jl`, kind `radiography`): a slab's transmission
# against Beer–Lambert in closed form along each pixel's ray, nested surfaces against the
# shell between them, the open field and the photon count against the chain by hand,
# point projection against a hand sum over the source elements with the chords of an
# axis-aligned box clipped in closed form (its magnification and source blur), a pinhole
# edge's blur against the penumbra, the gates' weights against closed forms, and the
# built-in definition.

using Briareus, LinearAlgebra, Test

# a closed box with outward fronts
function box_mesh(x0, x1, y0, y1, z0, z1)
    m = Mesh()
    v(x, y, z) = Vec3(x, y, z)
    c = v((x0 + x1) / 2, (y0 + y1) / 2, (z0 + z1) / 2)
    for f in ([v(x0, y0, z0), v(x1, y0, z0), v(x1, y1, z0), v(x0, y1, z0)],
              [v(x0, y0, z1), v(x1, y0, z1), v(x1, y1, z1), v(x0, y1, z1)],
              [v(x0, y0, z0), v(x1, y0, z0), v(x1, y0, z1), v(x0, y0, z1)],
              [v(x0, y1, z0), v(x1, y1, z0), v(x1, y1, z1), v(x0, y1, z1)],
              [v(x0, y0, z0), v(x0, y1, z0), v(x0, y1, z1), v(x0, y0, z1)],
              [v(x1, y0, z0), v(x1, y1, z0), v(x1, y1, z1), v(x1, y0, z1)])
        add_polygon!(m, f; facing = sum(f) / 4 - c)
    end
    return m
end

# the length of the segment a → b inside an axis-aligned box (slab clipping)
function clip_length(a, b, lo, hi)
    r = b - a
    s0, s1 = 0.0, 1.0
    for k in 1:3
        if r[k] == 0
            lo[k] <= a[k] <= hi[k] || return 0.0
        else
            t1, t2 = (lo[k] - a[k]) / r[k], (hi[k] - a[k]) / r[k]
            s0, s1 = max(s0, min(t1, t2)), min(s1, max(t1, t2))
        end
    end
    return max(s1 - s0, 0.0) * norm(r)
end

# a target of a one-sided foil (width w, standoff D, facing +z) and `parts`, the foil
# uniformly irradiated at `Ptot` TW
function backlit_scene(parts...; w = 4.0, D = 10.0, n = 8, Ptot = 10.0)
    t = Target([backlighter_foil((0.0, 0.0); standoff = D, width = w, n, two_sided = false), parts...])
    P = zeros(nelements(t))
    fe = elements(t, :backlighter)
    P[fe] .= Ptot .* t.elem_area[fe] ./ sum(t.elem_area[fe])
    return t, Briareus.Deposition(P, Float64[], [Ptot], Ptot)
end

# the object-plane point of pixel (i, j)
function pixel_point(rg, i, j)
    k = rg.camera
    x = (j - 0.5) / k.nx * k.width - k.width / 2
    y = k.height / 2 - (i - 0.5) / k.ny * k.height
    return x * k.right + y * k.up
end

@testset "backlit radiography" begin
    hνJ = 10.279 * 1.602176634e-16
    filters = [:Be => 0.1]
    Tf = exp(-mass_attenuation(:Be, 10.279) * 1.848 * 0.01)
    η = 1 - exp(-henke_mass_attenuation(:Si, 10279.0) * 2.329 * 0.0005)   # 5 µm
    κAl = mass_attenuation(:Al, 10.279) * 2.699 * 0.1                 # 1/mm
    p, L, dpin = 100.0, 1000.0, 0.02
    pin(; kw...) = Radiograph(0, 0; line = :Ge, distance = L, pinhole = dpin, pinhole_distance = p, width = 2.0,
                              pixels = 64, filters, detector = (:Si, 0.005), ce = 0.01, pulse = 1.0,
                              angular = :lambertian, up = Vec3(0, 1, 0), kw...)

    @testset "pinhole: Beer–Lambert, nested parts, the photon count" begin
        t, dep = backlit_scene(Component(:slab, box_mesh(-3, 3, -3, 3, -0.05, 0.05), ablator()))
        rg = pin()
        Pp = Vec3(0, 0, p)
        cosα = [dot(normalize(Pp - pixel_point(rg, i, j)), Vec3(0, 0, 1)) for i in 1:64, j in 1:64]
        r = radiograph(rg, t, dep; object = [:slab => :Al])
        @test r.magnification == (L - p) / p && r.weights == [1.0] && isempty(r.times)
        @test all(isapprox.(r.transmission, exp.(-κAl * 0.1 ./ cosα); rtol = 1e-7))
        # parts not in the object obstruct as hardware: nothing reaches the pinhole
        @test all(iszero, radiograph(rg, t, dep).open)
        # the outer box filled, the inner one removed: the 0.1 mm between them
        tn, depn = backlit_scene(Component(:outer, box_mesh(-3, 3, -3, 3, -0.1, 0.1), ablator()),
                         Component(:inner, box_mesh(-2, 2, -2, 2, -0.05, 0.05), ablator()))
        rn = radiograph(rg, tn, depn; object = [:outer => (:Al, 2.699), :inner => (:Al, -2.699)])
        @test all(isapprox.(rn.transmission, exp.(-κAl * 0.1 ./ cosα); rtol = 1e-7))
        # the open field of a Lambertian foil, R = E_x/(πA) through the pinhole: by hand
        # with cos⁴α off the axis (a flat detector normal to the line of sight)
        Ex = 0.01 * 1e3 * 1.0 * 10.0
        hand = Ex / (π * 16.0) * π * (dpin / 2)^2 / (L - p)^2 * Tf * η / hνJ .* cosα .^ 4
        @test all(isapprox.(r.open[10:55, 10:55], hand[10:55, 10:55]; rtol = 1e-12))
        @test r.images[1] ≈ r.open .* r.transmission rtol = 1e-12
        @test r.ce == 0.01 && r.sigma == 0 && r.line == 10.279
    end

    @testset "pinhole: an opaque edge blurred by the penumbra" begin
        # a half-plane x > 0 (along the image's right) that stops everything; pinhole
        # 50 µm at M = 9: the penumbra D = d(1 + 1/M) on the object plane, diffraction
        # 1.03 λ p/d = 0.25 µm
        dp = 0.05
        rg = pin(pinhole = dp, pixels = 200)
        @test rg.camera.right ≈ Vec3(1, 0, 0)
        t, dep = backlit_scene(Component(:edge, box_mesh(0, 3, -3, 3, -0.005, 0.005), ablator()))
        r = radiograph(rg, t, dep; object = [:edge => (:Au, 1e4)])
        D = dp * (1 + 1 / 9)
        h = 2.0 / 200
        x = [pixel_point(rg, 100, j)[1] for j in 1:200]            # chamber x along the row
        τ = r.transmission[100, :]
        part = findall(v -> 1e-9 < v < 1 - 1e-9, τ)
        @test all(abs.(x[part]) .< D / 2 + h / 2) && maximum(abs.(x[part])) > D / 2 - 1.5h
        # the blur moves no shadow: the row's blocked length is the box's 1 mm in the field
        @test sum(1 .- τ) * h ≈ 1.0 rtol = 1e-6
    end

    @testset "point projection: magnification and the source's blur" begin
        # a 0.4 mm foil 10 mm behind TCC, the detector 500 mm in front: M = 51, and an
        # Al box over x > 0.1 (in chamber x) 20 µm thick
        w, s, Ld = 0.4, 10.0, 500.0
        t, dep = backlit_scene(Component(:edge, box_mesh(0.1, 3, -3, 3, -0.01, 0.01), ablator()); w, D = s)
        rg = Radiograph(0, 0; line = :Ge, distance = Ld, width = 1.0, pixels = 41, ce = 0.01, pulse = 1.0,
                        up = Vec3(0, 1, 0))
        r = radiograph(rg, t, dep; object = [:edge => :Al])
        @test r.magnification ≈ (s + Ld) / s rtol = 1e-14
        # by hand: each foil element a point source (isotropic, E_e/4π sr⁻¹) toward the
        # detector point on the line from the foil's center through the pixel's point, at
        # the ray's incidence on the flat detector (normal +z), attenuated by the box's chord
        fe = elements(t, :backlighter)
        Ee = 0.01 * 1e3 .* dep.element[fe]
        fc = Vec3(0, 0, -s)
        expect = map(Iterators.product(1:41, 1:41)) do (i, j)
            q = pixel_point(rg, i, j)
            Pd = fc + (q - fc) * ((Ld + s) / (q - fc)[3])
            sum(Ee[k] / (4π * sum(abs2, Pd - t.elem_centroid[e])) * normalize(Pd - t.elem_centroid[e])[3] *
                exp(-κAl * clip_length(t.elem_centroid[e], Pd, Vec3(0.1, -3, -0.01), Vec3(3, 3, 0.01)))
                for (k, e) in enumerate(fe)) / hνJ
        end
        @test all(isapprox.(r.images[1], expect; rtol = 1e-10))
        # the open field at the center: a point source at s + L, to the foil's extent
        @test r.open[21, 21] ≈ sum(Ee) / (4π * (s + Ld)^2) / hνJ rtol = 1e-4
        # the source's blur on the object plane, w L/(s + L): partial transmission only
        # within it of the edge (element centroids span w(1 − 1/8))
        xs = [pixel_point(rg, 21, j) ⋅ Vec3(1, 0, 0) for j in 1:41]
        τ = r.transmission[21, :]
        T0 = exp(-κAl * 0.02)
        part = findall(v -> T0 * (1 + 1e-9) < v < 1 - 1e-9, τ)
        @test !isempty(part) && all(abs.(xs[part] .- 0.1) .< w / 2 * Ld / (s + Ld) + 1e-3)   # + the box's side
    end

    @testset "off the axis: the flat detector's obliquity against direct integration" begin
        # pixels at α = 0.25 rad from the line of sight (the last column, 5 pixels), the
        # detector flat and normal to it; a wide uniform foil, so every ray sees it
        α, p5 = 0.25, 100.0
        W = p5 * tan(α) / 0.4
        t, dep = backlit_scene(; w = 100.0, n = 4)
        A = 100.0^2
        Ex = 0.01 * 1e3 * 1.0 * 10.0
        # pinhole: the fluence is R ∫ cosθ_d cosθ_p / D² dA over the pinhole disk (normal
        # +z) from the detector point, by Gauss–Legendre in radius and uniform azimuths
        rg = Radiograph(0, 0; line = :Ge, distance = L, pinhole = dpin, pinhole_distance = p5, width = W,
                        pixels = 5, ce = 0.01, pulse = 1.0, angular = :lambertian, up = Vec3(0, 1, 0))
        q = pixel_point(rg, 3, 5)
        @test isapprox(acos(normalize(Vec3(0, 0, p5) - q)[3]), α; rtol = 1e-12)
        Pp = Vec3(0, 0, p5)
        Pd = Pp + (Pp - q) * ((L - p5) / p5)
        x, wx = Briareus._gauss_legendre(16)
        a = dpin / 2
        G = sum(wx[i] * a / 2 * (ρ = a * (x[i] + 1) / 2; ρ) * 2π / 64 *
                (r = Pd - (Pp + Vec3(ρ * cos(φ), ρ * sin(φ), 0.0)); r[3]^2 / sum(abs2, r)^2)
                for i in eachindex(x), φ in 2π .* (0:63) ./ 64)
        ref = Ex / (π * A) * G / hνJ
        rp = radiograph(rg, t, dep)
        @test isapprox(rp.open[3, 5], ref; rtol = 1e-9) && isapprox(rp.open[3, 3] * cos(α)^4, ref; rtol = 1e-9)
        # point projection: each element's fluence per detector area, E_e/4π times the
        # solid angle of a small detector square about the detector point over its area
        # (two triangles, Van Oosterom–Strackee, the triple product from the edges)
        function Ωtri(R1, R2, R3)
            r1, r2, r3 = norm(R1), norm(R2), norm(R3)
            2atan(abs(dot(R1, cross(R2 - R1, R3 - R1))),
                  r1 * r2 * r3 + dot(R1, R2) * r3 + dot(R1, R3) * r2 + dot(R2, R3) * r1)
        end
        s5, Lp = 10.0, 200.0
        rj = Radiograph(0, 0; line = :Ge, distance = Lp, width = s5 * tan(α) / 0.4,
                        pixels = 5, ce = 0.01, pulse = 1.0, up = Vec3(0, 1, 0))
        tj, dj = backlit_scene(; w = 1.0, n = 2, D = s5)
        fe = elements(tj, :backlighter)
        Ee = 0.01 * 1e3 .* dj.element[fe]
        fc = Vec3(0, 0, -s5)
        qj = pixel_point(rj, 3, 5)
        Pdj = fc + (qj - fc) * ((Lp + s5) / (qj - fc)[3])
        h = 1e-3
        c1, c2, c3, c4 = (Pdj + Vec3(u, v, 0.0) for (u, v) in ((-h, -h), (h, -h), (h, h), (-h, h)))
        refj = sum(Ee[k] / 4π * (Ωtri(c1 - tj.elem_centroid[e], c2 - tj.elem_centroid[e], c3 - tj.elem_centroid[e]) +
                                 Ωtri(c1 - tj.elem_centroid[e], c3 - tj.elem_centroid[e], c4 - tj.elem_centroid[e])) /
                   (2h)^2 for (k, e) in enumerate(fe)) / hνJ
        @test isapprox(acos(normalize(Pdj - fc)[3]), α; rtol = 1e-12)
        @test isapprox(radiograph(rj, tj, dj).open[3, 5], refj; rtol = 1e-8)
    end

    @testset "the CE chain: the open field at TCC is backlighter_signal's" begin
        t, dep = backlit_scene(; w = 1.0, n = 4)
        rg = Radiograph(0, 0; line = :Ge, distance = 500.0, width = 0.5, pixels = 5, pulse = 2.0, prepulse = true)
        r = radiograph(rg, t, dep)
        s = backlighter_signal(t, dep, (0.0, 0.0); line = :Ge, pulse = 2.0, prepulse = true, distance = 500.0,
                               object = nothing)
        # backlighter_signal's fluence is normal to each ray; the radiograph's, per area of
        # the flat detector, is times each element's incidence cosine (isotropic, uniform
        # foil: weights 1/r² per element)
        Pd = Vec3(0, 0, 500.0)
        rs = [Pd - t.elem_centroid[e] for e in elements(t, :backlighter)]
        inc = sum(r -> r[3] / norm(r)^3, rs) / sum(r -> 1 / norm(r)^2, rs)
        @test 1 - 1e-6 < inc < 1 && isapprox(r.open[3, 3], s.photons * inc; rtol = 1e-12)
        @test r.ce ≈ s.ce rtol = 1e-14
        @test r.sigma == s.sigma && r.extrapolated == s.extrapolated
    end

    @testset "gated frames" begin
        t, dep = backlit_scene(; w = 1.0, n = 4)
        g = Radiograph(0, 0; line = :Ge, distance = 500.0, width = 0.5, pixels = 5, pulse = 2.0, ce = 0.01,
                       start = 1.0, mode = :gated, times = [2.0, 1.0], gate = 0.09)
        rg = radiograph(g, t, dep)
        # a Gaussian gate of FWHM w inside the emission records w √(π/4 ln 2)/pulse of it;
        # one at the pulse's start, half that
        full = 0.09 * sqrt(π / (4log(2))) / 2.0
        @test rg.weights[1] ≈ full rtol = 1e-6
        @test rg.weights[2] ≈ full / 2 rtol = 1e-6
        @test rg.times == [2.0, 1.0]
        @test sum(rg.images[1]) ≈ rg.weights[1] * sum(rg.open) rtol = 1e-12
        sq = radiograph(Radiograph(0, 0; line = :Ge, distance = 500.0, width = 0.5, pixels = 5, pulse = 2.0,
                                   ce = 0.01, start = 1.0, mode = :gated, times = [2.0], gate = 0.09,
                                   gate_shape = :square), t, dep)
        @test sq.weights[1] ≈ 0.09 / 2.0 rtol = 1e-3
        @test_throws ArgumentError Radiograph(0, 0; line = :Ge, distance = 500.0, width = 0.5, pulse = 2.0,
                                              mode = :gated)
        @test_throws ArgumentError Radiograph(0, 0; line = :Ge, distance = 500.0, width = 0.5, pulse = 2.0,
                                              pinhole = 0.02)
    end

    @testset "definitions" begin
        issues = filter(i -> basename(i.file) == "nif_inflight_radiography.toml", check_diagnostics())
        @test Set(i.key for i in issues if i.level === :missing) ⊇ Set(["line_of_sight.phi_deg", "filters"])
        # a mock-up completing it: every value below is invented for the test, not a
        # facility value
        mktempdir() do dir
            write(joinpath(dir, "mock.toml"), """
                schema = 1
                id = "mock_inflight_radiography"
                extends = "nif_inflight_radiography"
                public = false
                note = "mock-up for a test; no value is a facility value"
                filters = [{ material = "Be", thickness_um = 100 }]
                [line_of_sight]
                phi_deg = 0.0
                [backlighter]
                pulse_ns = 2.0
                [optic]
                diameter_um = 20
                distance_mm = 100
                [detector]
                absorber = "Si"
                absorber_um = 5
                [view]
                width_mm = 2.0
                pixels = 64
                [recording]
                times_ns = [1.0, 1.5]
                """)
            paths = [diagnostic_paths()[1], dir]
            m = diagnostic("mock_inflight_radiography"; paths)
            @test m isa Radiograph && m.mode === :gated && m.gate == 0.09 && m.line == kshell_line_energy(:Ge)
            @test m.distance ≈ 100 * 9 rtol = 1e-14
        @test m.prepulse
            @test m.filter_transmission ≈ exp(-mass_attenuation(:Be, m.line) * 1.848 * 0.01) rtol = 1e-12
            @test !any(i -> occursin(dir, i.file) && i.level !== :warning, check_diagnostics(; paths))
            # the same definition as a point projection
            write(joinpath(dir, "proj.toml"), """
                schema = 1
                id = "mock_projection"
                extends = "mock_inflight_radiography"
                [optic]
                kind = "projection"
                """)
            issues = check_diagnostics(; paths)
            @test any(i -> occursin("proj.toml", i.file) && i.key == "optic.diameter_um", issues)
        end
    end
end
