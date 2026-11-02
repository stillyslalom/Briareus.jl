# Hohlraum profiles (`hohlraum(; profile)`): a profile of one cylindrical section is the
# cylindrical hohlraum bit for bit; frustraums and multi-section walls mesh as planar
# sections that meet vertex for vertex; a closed frustraum's exchanges against the
# coaxial-disk closed form, and its drive and LEH loss against Monte Carlo transport.

using Briareus, LinearAlgebra, Statistics, Test

@testset "Hohlraum profiles" begin
    same(a, b) = length(a.components) == length(b.components) &&
                 all(((c, d),) -> c.name === d.name && c.mesh.vertices == d.mesh.vertices &&
                                  c.mesh.triangles == d.mesh.triangles && c.mesh.patches == d.mesh.patches &&
                                  c.material == d.material && c.kind === d.kind && c.surface == d.surface &&
                                  c.part === d.part,
                     zip(a.components, b.components)) && a.closed == b.closed
    names_of(t) = [c.name for c in t.components]

    # one cylindrical section is the cylindrical hohlraum, windows, lips and ends included
    cyl = HohlraumProfile([-3.0, 3.0], [2.0, 2.0])
    kw = (; leh_radius = 1.2, capsule_radius = 0.6, h = 0.5, capsule_h = 0.4)
    win = [WallWindow(; phi = 30, z = 1.0, width = 1.2, height = 0.8, open = true),
           WallWindow(; phi = 200, z = -1.5, diameter = 1.0, gap = 0.1)]
    lip = GapClosure(closure = 1.0, width = 0.5, after = 0.85)
    @test same(hohlraum(; profile = cyl, kw...), hohlraum(; length = 6.0, radius = 2.0, kw...))
    @test same(hohlraum(; profile = cyl, kw..., windows = win, leh_closure = lip, bottom_end = :closed),
               hohlraum(; length = 6.0, radius = 2.0, kw..., windows = win, leh_closure = lip, bottom_end = :closed))
    @test same(hohlraum(; profile = cyl), hohlraum(; length = 6.0, radius = 2.0))

    # what a profile refuses
    @test_throws ArgumentError HohlraumProfile([-1.0, 2.0], [1.0, 1.0])            # not centered
    @test_throws ArgumentError HohlraumProfile([-1.0, 0.0, 1.0], [1.0, 1.0])       # a radius short
    @test_throws ArgumentError HohlraumProfile([-1.0, 1.0, 0.5], [1.0, 1.0, 1.0])  # not sorted
    @test_throws ArgumentError HohlraumProfile([-1.0, -1.0, 1.0], [1.0, 1.5, 1.5]) # ends in a step
    @test_throws ArgumentError HohlraumProfile([-1.0, 1.0], [1.0, 0.0])
    fr = frustraum_profile(; length = 6.0, waist_radius = 2.5, end_radius = 1.5)
    @test_throws ArgumentError hohlraum(; profile = fr, length = 6.0)
    @test_throws ArgumentError hohlraum(; profile = fr, corner_radius = 0.2)
    @test_throws ArgumentError hohlraum(; profile = fr, leh_radius = 1.6)
    @test_throws ArgumentError hohlraum(; profile = fr, leh_radius = 1.0, capsule_radius = 2.4)
    @test_throws ArgumentError hohlraum(; profile = fr, leh_radius = 1.0, windows = [WallWindow(; phi = 0, diameter = 0.5)])

    # a frustraum: two cones, each with its exact area, meeting vertex for vertex at the waist
    t = hohlraum(; profile = fr, leh_radius = 1.0, capsule_radius = 0.8, h = 0.4, capsule_h = 0.3)
    @test names_of(t) == [:wall_1, :wall_2, :cap_top, :cap_bottom, :leh_top, :leh_bottom, :capsule] && t.closed
    @test all(c -> c.part === :wall, t.components[1:2]) && t.components[1].surface isa ConeSurface
    @test elements(t, :wall) == vcat(elements(t, :wall_1), elements(t, :wall_2))
    area(t, n) = sum(t.elem_area[elements(t, n)])
    @test area(t, :wall_1) ≈ π * (1.5 + 2.5) * hypot(1.0, 3.0) rtol = 1e-12
    @test area(t, :wall_1) ≈ area(t, :wall_2) rtol = 1e-12
    @test area(t, :cap_top) + area(t, :leh_top) ≈ π * 1.5^2 rtol = 1e-12
    waist(c) = Set(v for v in c.mesh.vertices if v[3] == 0)
    @test waist(t.components[1]) == waist(t.components[2]) && length(waist(t.components[1])) >= 40
    # the default entrance holes for the keep-out, one per end
    @test ApertureKeepOut(t).names == [:leh_top, :leh_bottom]

    # four sections with a step in radius at ±1 mm: rings face the cavity, a window cut
    # into the wide middle section, one refused on a cone
    st = HohlraumProfile([-3.0, -1.0, -1.0, 1.0, 1.0, 3.0], [1.6, 2.0, 2.2, 2.2, 2.0, 1.6])
    t = hohlraum(; profile = st, leh_radius = 0.8, capsule_radius = 0.6, h = 0.4,
                 windows = [WallWindow(; phi = 90, z = 0.2, diameter = 0.8)])
    @test names_of(t)[1:5] == [:wall_1, :wall_2, :wall_3, :wall_4, :wall_5] && :window_1 in names_of(t)
    @test t.components[2].surface isa DiskSurface && t.components[3].surface isa CylinderSurface
    @test area(t, :window_1) ≈ π * 0.4^2 rtol = 1e-12
    @test area(t, :wall_2) ≈ π * (2.2^2 - 2.0^2) rtol = 1e-12
    @test area(t, :wall_3) + area(t, :window_1) ≈ 2π * 2.2 * 2.0 rtol = 1e-12
    # the lower ring faces up (+z), the upper one down
    nz(t, n) = [t.elem_side[e] * t.tri_normal[first(t.elem_tris[e])][3] for e in elements(t, n)]
    @test all(>(0.999), nz(t, :wall_2)) && all(<(-0.999), nz(t, :wall_4))
    @test_throws ArgumentError hohlraum(; profile = st, leh_radius = 0.8, windows = [WallWindow(; phi = 0, z = -2.0, diameter = 0.5)])
    @test_throws ArgumentError hohlraum(; profile = st, leh_radius = 0.8, windows = [WallWindow(; phi = 0, z = 0.9, diameter = 0.5)])

    # A closed frustraum (waist radius rw, end radius re, length L; ends closed, no
    # capsule) against closed forms built from the coaxial-disk factor F(a → b, d) alone:
    # the ends see each other unobstructed (F_ee = F(re → re, L)); the waist plane splits
    # the cavity, so an end sees the far cone through the waist disk,
    # F(re → rw, L/2) − F_ee, and the near cone with the rest; the near cone's view of the
    # far one is what passes the waist disk less what lands on the far end (reciprocity),
    # and its view of itself the remainder.
    Fdisk(a, b, d) = (A = a / d; B = b / d; X = 1 + (1 + B^2) / A^2; (X - sqrt(X^2 - 4 * (B / A)^2)) / 2)
    L, rw, re = 4.0, 2.0, 1.2
    Ae, Aw, Ac = π * re^2, π * rw^2, π * (re + rw) * hypot(rw - re, L / 2)
    Fee, Few = Fdisk(re, re, L), Fdisk(re, rw, L / 2)
    c1c2 = (Aw - Ae * Few - Ae * (Few - Fee)) / Ac
    exact = (; ee = Fee, near = 1 - Few, far = Few - Fee, c1c2,
             c1c1 = 1 - Ae * (1 - Few) / Ac - Ae * (Few - Fee) / Ac - c1c2)
    p = frustraum_profile(; length = L, waist_radius = rw, end_radius = re)
    t = hohlraum(; profile = p, h = 0.25, leh_radius = 0.6, capsule_radius = 0.0, top_end = :closed,
                 bottom_end = :closed)
    F = view_factors(t)
    A = t.elem_area
    top, bot, c1, c2 = elements(t, :cap_top), elements(t, :cap_bottom), elements(t, :wall_1), elements(t, :wall_2)
    flux(I, J) = sum(A[i] * F[i, j] for i in I, j in J) / sum(A[I])
    @test sum(A[c1]) ≈ Ac rtol = 1e-12
    @test maximum(abs, A .* F .- (A .* F)') < 1e-14 * maximum(A .* F)
    @test flux(top, bot) ≈ exact.ee rtol = 5e-4                   # measured 7.3e-5
    @test flux(top, c2) ≈ exact.near rtol = 2e-3                  # 9.1e-4 (the top's own cone)
    @test flux(top, c1) ≈ exact.far rtol = 2e-4                   # 1.3e-5
    @test flux(c1, c2) ≈ exact.c1c2 rtol = 3e-3                   # 1.5e-3
    @test flux(c1, c1) ≈ exact.c1c1 rtol = 1.5e-3                 # 5.4e-4

    # Drive and loss of a frustraum with entrance holes and a capsule, the wall driven:
    # the energy balance closes, and the LEH loss and capsule drive agree with Monte Carlo
    # transport on the same facets (`reemit = :element`, the element model with exact
    # view factors), as the cylinder of the same length and holes does (measured on the
    # balanced model at h = 0.3: frustraum 1.2e-3 and 1.9e-3, cylinder 1.8e-3 and 2.6e-3,
    # statistical error 6e-4).
    au = Material(:Au; albedo = 0.8, conversion = 0.8)
    t0 = hohlraum(; profile = p, h = 0.3, leh_radius = 0.7, capsule_radius = 0.6, capsule_h = 0.3, wall = au)
    t = Target([Component(c.name, c.mesh, c.material; kind = c.kind, part = c.part) for c in t0.components];
               closed = true)
    P = [t.components[t.elem_comp[e]].part === :wall ? t.elem_area[e] : 0.0 for e in 1:nelements(t)]
    sol = solve(RadiationModel(t), P)
    eb = energy_balance(sol)
    @test sum(values(eb.absorbed)) + sum(values(eb.escaped)) + eb.lost ≈ eb.source rtol = 1e-12
    leh, cap = elements(t, :leh_top, :leh_bottom), elements(t, :capsule)
    A = t.elem_area
    _, batch = montecarlo_incident(t, P, vcat(leh, cap); paths = 4000, batches = 16, reemit = :element)
    nl = length(leh)
    @test sum(A[leh] .* sol.H[leh]) ≈ mean(sum(A[leh] .* batch[b, 1:nl]) for b in axes(batch, 1)) rtol = 5e-3
    @test sum(A[cap] .* sol.H[cap]) ≈ mean(sum(A[cap] .* batch[b, nl+1:end]) for b in axes(batch, 1)) rtol = 5e-3
end
