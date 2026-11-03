# Hohlraum ends (`hohlraum(; top_end, bottom_end)`): the defaults unchanged bit for bit,
# per-end entrance-hole radii, closed ends, and open ends (a halfraum) against the
# coaxial-disk view factor and an energy balance through an aperture closing the end;
# the face flux objective of a flat observed face.

using Briareus, LinearAlgebra, ForwardDiff, Test

@testset "Hohlraum ends" begin
    same(a, b) = length(a.components) == length(b.components) &&
                 all(((c, d),) -> c.name === d.name && c.mesh.vertices == d.mesh.vertices &&
                                  c.mesh.triangles == d.mesh.triangles && c.mesh.patches == d.mesh.patches &&
                                  c.material == d.material && c.kind === d.kind && c.surface == d.surface,
                     zip(a.components, b.components)) && a.closed == b.closed
    names_of(t) = [c.name for c in t.components]
    kw = (; length = 6.0, radius = 2.0, leh_radius = 1.2, capsule_radius = 0.6, h = 0.5, capsule_h = 0.4)
    # an end given the default radius is the default end, exactly (with a closing lip too)
    @test same(hohlraum(; kw..., top_end = 1.2, bottom_end = 1.2), hohlraum(; kw...))
    lip = GapClosure(closure = 1.0, width = 0.5, after = 0.85)
    @test same(hohlraum(; kw..., leh_closure = lip, bottom_end = 1.2), hohlraum(; kw..., leh_closure = lip))

    # per-end radii: the top hole of 0.8 mm, its cap from it to the wall
    t = hohlraum(; kw..., top_end = 0.8)
    area(n) = sum(t.elem_area[elements(t, n)])
    @test area(:leh_top) ≈ π * 0.8^2 rtol = 1e-12
    @test area(:cap_top) + area(:leh_top) ≈ area(:cap_bottom) + area(:leh_bottom) rtol = 1e-12
    @test maximum(v -> hypot(v[1], v[2]), hohlraum(; kw..., top_end = 0.8).components[2].mesh.vertices) ≈ 2.0
    # a closed end: no hole, its cap the whole end of the wall's material
    t = hohlraum(; kw..., top_end = :closed)
    @test !(:leh_top in names_of(t)) && t.closed
    @test sum(t.elem_area[elements(t, :cap_top)]) ≈ 4π rtol = 1e-12
    # an open end: neither cap nor hole, the target not declared closed
    t = hohlraum(; kw..., bottom_end = :open)
    @test names_of(t) == [:wall, :cap_top, :leh_top, :capsule] && !t.closed
    @test_throws ArgumentError hohlraum(; kw..., top_end = :ajar)
    @test_throws ArgumentError hohlraum(; kw..., top_end = 2.5)
    @test_throws ArgumentError hohlraum(; kw..., top_end = :open, corner_radius = 0.3)
    @test_throws ArgumentError hohlraum(; kw..., top_end = :open, sombrero = true)

    # A halfraum closed at the top and open at the bottom (radius R, length L, no
    # capsule) against the coaxial-disk closed form (r₁ = r₂ = R, L apart): the cap sees
    # the open end with F_dd, and by reciprocity and summation the wall sends
    # A_d (1 − F_dd) / A_w of its emission out through it.
    R, L = 1.5, 3.0
    X = 2 + (L / R)^2
    Fdd = (X - sqrt(X^2 - 4)) / 2
    t = hohlraum(; length = L, radius = R, leh_radius = 0.75, capsule_radius = 0.0, h = 0.2,
                 top_end = :closed, bottom_end = :open)
    F = view_factors(t)
    A = t.elem_area
    escape = 1 .- vec(sum(F; dims = 2))
    cap, wall = elements(t, :cap_top), elements(t, :wall)
    Ad, Aw = π * R^2, 2π * R * L
    @test sum(A[cap] .* escape[cap]) / Ad ≈ Fdd rtol = 2e-3          # measured 0.95e-3
    @test sum(A[wall] .* escape[wall]) / (Ad * (1 - Fdd)) ≈ 1 rtol = 2.5e-3   # 1.3e-3

    # Energy balance: the wall driven, the power lost through the open end equals the
    # power escaping through an aperture disk that closes the end (computed from its own
    # view factors), and each balance closes.
    au = Material(:Au; albedo = 0.8, conversion = 0.8)
    t = hohlraum(; length = L, radius = R, leh_radius = 0.75, capsule_radius = 0.0, h = 0.2,
                 top_end = :closed, bottom_end = :open, wall = au)
    shut = Target([t.components; aperture(:end, disk_mesh(R, -L / 2; nr = 8, nphi = 32, up = true);
                                          surface = DiskSurface(center = (0, 0, -L / 2)))]; closed = false)
    P(t) = [t.components[t.elem_comp[e]].name === :wall ? t.elem_area[e] : 0.0 for e in 1:nelements(t)]
    eo = energy_balance(solve(RadiationModel(t), P(t)))
    es = energy_balance(solve(RadiationModel(shut), P(shut)))
    @test eo.lost / es.escaped[:end] ≈ 1 rtol = 3e-3                # measured 1.45e-3
    @test abs(es.lost) < 3e-3 * es.escaped[:end]
    @test sum(values(eo.absorbed)) + eo.lost ≈ eo.source rtol = 1e-12
end

# The face objective of a flat observed face (an opacity sample at a halfraum's open end)
# against a direct forward solve, its gradient against finite differences, and an
# optimization that lowers it.
@testset "face flux objective" begin
    hr = hohlraum(; length = 3.0, radius = 1.5, leh_radius = 1.0, capsule_radius = 0.0, h = 0.25,
                  top_end = 1.0, bottom_end = :open)
    smp = foil(:sample, (0, 0, -1.5), (0, 0, 1); width = 1.0, n = 8, two_sided = false,
               material = Material(:sample; albedo = 0.3, conversion = 0.0))
    t = Target([hr.components; smp]; closed = false)
    up = [remake(b; aim = (0.0, 0.0, 1.5)) for b in nif_quads(power = 4.0) if b.port[3] > 0]
    p = GroupPointing(up; fields = (:z, :r), mirror = false)
    s = Scenario(t, p; observe = :sample)
    x0 = initial_parameters(p)
    ev = evaluate(s, x0)
    # the value: the RMS of the radiosity solve's incident flux on the face, and its mean
    H = solve(s.model, ev.deposition).H[s.observed]
    A = t.elem_area[s.observed]
    H̄ = sum(A .* H) / sum(A)
    σ = sqrt(sum(A .* (H .- H̄) .^ 2) / sum(A)) / H̄
    shape = face_flux_objective(; clip = 0.0, missed = 0.0, direct = 0.0)
    @test shape(ev) ≈ σ^2 rtol = 1e-12
    @test face_mean_flux(ev) ≈ H̄ rtol = 1e-12
    @test face_flux(ev) ≈ H rtol = 1e-12
    held = face_flux_objective(; clip = 0.0, missed = 0.0, direct = 0.0, flux = 0.8H̄, flux_weight = 2.0)
    @test held(ev) ≈ σ^2 + 2 * (1 / 0.8 - 1)^2 rtol = 1e-12
    # the penalties: the light through the open end and on the sample itself
    full = face_flux_objective()
    @test full(ev) ≈ σ^2 + 10 * (clipped_fraction(ev)^2 + missed_fraction(ev)^2 +
                                 power_fraction(ev, :sample)^2) rtol = 1e-12
    # the term scaled by a reference RMS: the start's nonuniformity weighs `uniformity`
    @test face_nonuniformity(ev) ≈ σ rtol = 1e-12
    scaled = face_flux_objective(; uniformity = 3.0, rms_scale = 2σ, clip = 0.0, missed = 0.0, direct = 0.0)
    @test scaled(ev) ≈ 3 / 4 rtol = 1e-12
    @test_throws ArgumentError face_flux_objective(; uniformity = -1.0)
    @test_throws ArgumentError face_flux_objective(; flux = 0.0)
    @test_throws ArgumentError face_flux_objective(; rms_scale = 0.0)
    # the gradient (forward-mode duals through the flux response) against central differences
    f(x) = full(evaluate(s, x))
    x1 = x0 .+ 0.1 .* (1:length(x0)) ./ length(x0)
    g = ForwardDiff.gradient(f, x1)
    δ = 1e-5
    fd = [(f(x1 .+ δ .* (eachindex(x1) .== i)) - f(x1 .- δ .* (eachindex(x1) .== i))) / 2δ for i in eachindex(x1)]
    @test maximum(abs.(g .- fd)) < 1e-7 * maximum(abs.(g))   # measured 2.4e-10
    # optimizing the shape at a held level lowers the face's RMS (measured 0.61 % → 0.046 %)
    lo, hi = bounds(p, x0; z = 1.0, r = 1.0)
    obj = face_flux_objective(; clip = 0.0, missed = 0.0, direct = 0.0, flux = H̄)
    res = optimize_pointing(s, obj; x0, lower = lo, upper = hi, iterations = 30)
    eo = evaluate(s, res.x)
    Ho = solve(s.model, eo.deposition).H[s.observed]
    @test res.value ≈ obj(eo) rtol = 1e-12
    @test rms_nonuniformity(Ho, A) < 0.5σ
end
