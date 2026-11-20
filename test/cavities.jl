# Targets with several cavities (assemblies phase C2). The wall temperature of each cavity
# (`wall_temperature`, `drive_history(; cavities)`) is checked against the area-weighted
# incident flux of its wall elements, reduced in this file directly from the solution. The
# history without `cavities` is unchanged, and a Dante field of view that leaves one cavity
# out receives nothing from it.

using Briareus, LinearAlgebra, Test
using Briareus: RadiationSolution

@testset "Several cavities: wall temperatures and the Dante field" begin
    # two closed hohlraums side by side, 12 mm apart along x, as the parts :a and :b with
    # their pieces named `<part>_wall`, … (the GUI's hohlraum part scheme)
    base = hohlraum(; length = 6.0, radius = 2.0, leh_radius = 1.2, capsule_radius = 0.0, h = 0.7)
    placed(p, x) = [translate(Component(Symbol(p, :_, c.name), c.mesh, c.material; kind = c.kind,
                                        two_sided = c.two_sided, surface = c.surface, part = p), (x, 0.0, 0.0))
                    for c in base.components]
    t = Target([placed(:a, -6.0); placed(:b, 6.0)])
    walls(p) = Tuple(Symbol(p, :_, s) for s in (:wall, :cap_top, :cap_bottom))
    cavities = [:a => walls(:a), :b => walls(:b)]
    every = (walls(:a)..., walls(:b)...)
    # the upper beams in a's top entrance hole, the lower in b's bottom one: two cavities
    # driven differently
    q = [remake(b; aim = b.port[3] > 0 ? (-6.0, 0.0, 3.0) : (6.0, 0.0, -3.0)) for b in nif_quads()]
    # (b's beams at a quarter of a's power)
    pulses = Dict(:up => PulseShape([0.0, 1.0, 3.0], [0.0, 2.0, 4.0]),
                  :down => PulseShape([0.0, 1.0, 3.0], [0.0, 0.5, 1.0]))
    side(b) = b.port[3] > 0 ? :up : :down
    m, op = RadiationModel(t), DepositionOperator(t, q)
    run(; kw...) = drive_history(m, op, τ -> with_pulses(q, pulses, τ; by = side), [0.5, 1.5, 2.5];
                                 observe = :a_wall, wall = every, record = sol -> copy(sol.H), kw...)
    h0 = run()
    h = run(; cavities = [cavities; :all => every])

    # The independent reduction finds the wall elements of each cavity from the component
    # names and takes Σ H A / Σ A of the recorded incident flux at each step.
    own(p) = [e for e in 1:nelements(t) if t.components[t.elem_comp[e]].name in walls(p)]
    Tr(H, es) = (sum(H[e] * t.elem_area[e] for e in es) / sum(t.elem_area[e] for e in es) / SIGMA_SB)^0.25
    @test all(isapprox(h.cavity_Tr[p], [Tr(H, own(p)) for H in h.record]; rtol = 1e-12) for p in (:a, :b))
    # two cavities, two temperatures; together they are the history's wall temperature
    @test all(abs.(h.cavity_Tr[:a] .- h.cavity_Tr[:b]) .> 1.0)
    @test h.cavity_Tr[:all] == h.wall_Tr
    # Without `cavities` the history is unchanged.
    @test isempty(h0.cavity_Tr) && all(k -> k === :cavity_Tr || isequal(h0[k], h[k]), keys(h))
    # `wall_temperature` on one solution
    sol = solve(m, deposit(op, q))
    @test wall_temperature(sol, walls(:b)...) ≈ Tr(sol.H, own(:b)) rtol = 1e-12
    # `per_component` gives the energy per piece; the pieces are not summed into the parts.
    # The sum over the pieces of each part equals the by-part history, and every other field
    # is unchanged.
    hc = run(; per_component = true)
    @test haskey(hc.absorbed, :a_cap_top) && !haskey(hc.absorbed, :a) && haskey(h0.absorbed, :a)
    bypart(d, p) = sum(v for (k, v) in d if startswith(String(k), String(p) * "_"))
    @test all(isapprox(bypart(getfield(hc, f), p), getfield(h0, f)[p]; rtol = 1e-12)
              for f in (:laser, :source, :absorbed), p in (:a, :b))
    @test all(k -> k in (:laser, :source, :absorbed, :escaped, :transmitted, :released, :mband) ||
                   isequal(h0[k], hc[k]),
              keys(h0))
    @test_throws ArgumentError run(; cavities = [:a => walls(:a), :a => walls(:b)])
    @test_throws ArgumentError run(; cavities = [walls(:a)])

    # Dante 1's line of sight, each cavity seen through its bottom entrance hole. A
    # blackbody in one cavity at a time (the outside of every one-sided surface is dark).
    d1 = port_angles(:nif, "P143-274")
    n = nelements(t)
    T0 = 250.0
    lit(p) = [c.kind !== :aperture && c.part === p ? SIGMA_SB * T0^4 : 0.0 for c in t.components[t.elem_comp]]
    z = zeros(n)
    bb(B) = RadiationSolution(RadiationModel(t; F = zeros(n, n)), z, z, B, z, z, z, z, 0.0)
    # a field of 8 mm about a's center leaves b out: b gives it nothing
    narrow = DanteView(d1...; through = :a_leh_bottom, field_of_view = 8.0, center = :a, pixels = 400)
    @test dante(narrow, bb(lit(:b))).intensity == 0 && !haskey(visible_area(narrow, t), :b)
    # the views through each hole (no field of view), and a 60 mm field at TCC (a NIF
    # Dante's, research notes) that takes in both cavities
    hole(p) = dante(DanteView(d1...; through = Symbol(p, :_leh_bottom)), bb(lit(p))).intensity
    wide = DanteView(d1...; through = :a_leh_bottom, field_of_view = 60.0, pixels = 1200)
    @test dante(narrow, bb(lit(:a))).intensity ≈ hole(:a) rtol = 1e-3       # measured 1.2e-4
    @test dante(wide, bb(lit(:a) .+ lit(:b))).intensity ≈ hole(:a) + hole(:b) rtol = 1e-3   # 1.4e-4
    @test issubset((:a, :b), keys(visible_area(wide, t)))
    # The built-in Dante 1, with its field at TCC and the Dante 1 channels, sees the whole
    # target. This target does not have a `leh_bottom` to quote Tr over, so only the power
    # is checked.
    fv = dante1()
    @test fv.field_of_view == 60.0 && fv.center == zero(fv.center) && fv.through === :leh_bottom &&
          (fv.theta, fv.phi) == (143.25, 274.0) && length(fv.channels) == 18
    r = dante(fv, bb(lit(:a) .+ lit(:b)))
    @test isnan(r.Tr) && isapprox(r.intensity, hole(:a) + hole(:b); rtol = 5e-3)   # 1.8e-3 (400 pixels)

    # Each channel can have a field of view. Channel 1 has a separate field, and the others
    # use the view field. The blackbodies are given as spectra, with the Planck fractions in
    # the groups of every element.
    E = spectral_edges()
    spec(B) = (G = length(E) - 1; Z = zeros(n, G);
               SpectralSolution(bb(B), E, Z, B .* planck_fractions(T0, E)', Z, Z, fill(T0, n), z, zeros(G), false))
    sa, sb, sab = spec(lit(:a)), spec(lit(:b)), spec(lit(:a) .+ lit(:b))
    chans(f) = [c.number == 1 ? Briareus._with_field(c, f) : Briareus._with_field(c, nothing) for c in dante1_channels()]
    chview(f; kw...) = DanteView(d1...; through = :a_leh_bottom, channels = chans(f), kw...)
    # A channel given the view field records the same values as without one.
    same = dante(chview(60.0; field_of_view = 60.0), sab)
    @test all(k -> isequal(same[k], dante(chview(nothing; field_of_view = 60.0), sab)[k]),
              (:intensity, :Tr, :spectrum, :channels, :color_temperature, :Tr_unfold))
    # With 8 mm about a in the 60 mm field, channel 1 sees nothing of b and records from a
    # the channel 1 signal of the 8 mm view. The rays are the same, and the in-group
    # spectral shape is the Tr_open of each view; the two differ by 1.5e-6 in this case. The
    # other channels see both cavities.
    eight(s) = dante(chview(8.0; field_of_view = 60.0, center = :a), s).channels
    @test eight(sb)[1] == 0 && all(>(0), eight(sb)[2:end])
    @test eight(sa)[1] ≈ dante(chview(nothing; field_of_view = 8.0, center = :a), sa).channels[1] rtol = 1e-4
    @test eight(sa)[1] > 0
    @test eight(sab)[2:end] == dante(chview(nothing; field_of_view = 60.0, center = :a), sab).channels[2:end]
    # The next views see only the hole of a. A channel field that contains the hole changes
    # nothing, and a field of 8 mm about TCC sees none of the hole, which lies 6 mm off TCC.
    hv(f) = dante(chview(f; center = zero(Briareus.Vec3)), sa).channels
    @test hv(60.0) == hv(nothing) && hv(8.0)[1] == 0 && hv(8.0)[2:end] == hv(nothing)[2:end]
end

@testset "Several cavities: per-cavity balancing" begin
    named(c, p) = Component(Symbol(p, :_, c.name), c.mesh, c.material; kind = c.kind,
                            two_sided = c.two_sided, surface = c.surface, part = p)
    rows(F, c) = vec(sum(F[c, c]; dims = 2))
    # A sphere of radius R with a polar hole (its first of `nθ` bands, equal in cos θ),
    # closed by a flat disk aperture: every point of a sphere sees each area A of it with
    # the view factor A / 4πR², so a wall emitting S per unit area with albedo α carries
    # the incident flux H = S w / (1 − α w) everywhere, w = (1 + cos θ₀)/2 its share of
    # the sphere (closed form, every bounce).
    function sphere_cavity(p, R, nθ, nφ, center; albedo)
        m = Briareus.Mesh()
        μs, φs = range(1.0, -1.0; length = nθ + 1), range(0.0, 2π; length = nφ + 1)
        pt(μ, φ) = (s = sqrt(max(0.0, 1 - μ^2)); R * Briareus.Vec3(s * cos(φ), s * sin(φ), μ))
        for i in 2:nθ, j in 1:nφ
            q = i == nθ ? [pt(μs[i], φs[j]), Briareus.Vec3(0, 0, -R), pt(μs[i], φs[j+1])] :
                [pt(μs[i], φs[j]), pt(μs[i+1], φs[j]), pt(μs[i+1], φs[j+1]), pt(μs[i], φs[j+1])]
            Briareus.add_polygon!(m, q; facing = -sum(q) / length(q))
        end
        c = μs[2]
        wall = Component(:wall, m, Material(:w; albedo, conversion = 1.0); surface = SphereSurface(R))
        hole = aperture(:leh, disk_mesh(R * sqrt(1 - c^2), R * c; nr = 2, nphi = nφ, up = false);
                        surface = DiskSurface(center = (0, 0, R * c)))
        return [translate(named(x, p), center) for x in (wall, hole)], (1 + c) / 2
    end
    # Two of them stand side by side, with a black disk in the open above the hole of a,
    # which the wall of a sees through the hole. The disk is an element in no cavity.
    R, α = 1.5, 0.8
    ca, w = sphere_cavity(:a, R, 10, 24, (-4.0, 0.0, 0.0); albedo = α)
    cb, _ = sphere_cavity(:b, R, 10, 24, (4.0, 0.0, 0.0); albedo = α)
    free = translate(Component(:free, disk_mesh(1.0, 0.0; nr = 3, nphi = 16, up = false),
                               Material(:k; albedo = 0.0, conversion = 0.0); surface = DiskSurface()), (-4.0, 0.0, R + 0.5))
    t = Target([ca; cb; free])
    a, b = elements(t, :a), elements(t, :b)
    F0 = view_factors(t)
    @test all(c -> maximum(abs.(rows(F0, c) .- 1)) > 1e-3, (a, b))            # 1.6e-2 on the polar band
    @test open_fraction(t; cavity = :a).fraction == 0 && open_fraction(t; cavity = b).fraction == 0
    @test open_fraction(t; cavity = :a_wall).fraction > 0.05                   # without its hole
    # Each cavity is balanced separately, and its rows sum to one over its elements.
    model(closed) = RadiationModel(t; F = copy(F0), closed)
    m = model([:a, :b])
    @test all(c -> maximum(abs.(rows(m.F, c) .- 1)) < 1e-10, (a, b))
    @test m.F == model([a, b]).F && isapprox(m.F' .* t.elem_area', m.F .* t.elem_area; rtol = 1e-12)
    @test model([1:nelements(t)]).F == model(true).F
    @test_throws ArgumentError model([:a, [first(a)]])
    @test_throws ArgumentError model(:a)
    # The mean incident flux on the wall is checked against the closed form. Per-cavity
    # balancing brings it closer, from 3.3e-3 to 3.7e-4. Balancing the whole target, which
    # includes the rows of the free disk and the view of the disk from the wall, is wrong by
    # 36%.
    wa, wb = elements(t, :a_wall), elements(t, :b_wall)
    P = zeros(nelements(t))
    P[wa] .= t.elem_area[wa]
    Hw = w / (1 - α * w)
    err(m) = (s = solve(m, P); sum(s.H[wa] .* t.elem_area[wa]) / sum(t.elem_area[wa]) / Hw - 1)
    raw, bal = err(model(false)), err(m)
    @test abs(bal) < 1e-3 && abs(bal) < abs(raw) / 5
    @test abs(err(model(true))) > 0.1
    # energy, exactly per cavity: a's source is what a's elements absorb and let out
    # through its hole, Σ (1 − α) H A; b receives nothing
    s = solve(m, P)
    out(es) = sum((1 - s.albedo[e]) * s.H[e] * t.elem_area[e] for e in es)
    @test out(a) ≈ sum(s.source .* t.elem_area) rtol = 1e-12
    @test out(b) == 0
    @test sum(s.H[elements(t, :free)]) > 0                    # (through the hole, outside the cavity)

    # Two halfraums joined by a package, black: each halfraum and the package face that
    # closes its open end form a cavity; the package's outer face is in none. From the
    # face, the coaxial-disk view factors (r₁ = 1.5 to r₂ = 1.0 and 1.5 mm, 2.5 mm apart)
    # give what leaves through the hole, lands on the cap and on the wall.
    h = 0.4
    kw = (; length = 2.5, radius = 1.5, leh_radius = 1.0, capsule_radius = 0.0, h,
          wall = Material(:au; albedo = 0.0, conversion = 1.0))
    nphi = hohlraum_counts(; length = 2.5, radius = 1.5, leh_radius = 1.0, h).nphi
    face = Material(:pk; albedo = 0.0, conversion = 1.0)
    pkg = [Component(:package, cylinder_mesh(1.5, -1.0, 1.0; nz = 5, nphi, inward = false), face;
                     surface = CylinderSurface(1.5)),
           Component(:package_top, disk_mesh(1.5, 1.0; nr = 4, nphi, up = true), face;
                     surface = DiskSurface(center = (0, 0, 1.0))),
           Component(:package_bottom, disk_mesh(1.5, -1.0; nr = 4, nphi, up = false), face;
                     surface = DiskSurface(center = (0, 0, -1.0)))]
    t = Target([[translate(named(x, :top), (0, 0, 2.25)) for x in hohlraum(; kw..., bottom_end = :open).components];
                [translate(named(x, :bot), (0, 0, -2.25)) for x in hohlraum(; kw..., top_end = :open).components]; pkg])
    cavs = [[elements(t, :top); elements(t, :package_top)], [elements(t, :bot); elements(t, :package_bottom)]]
    @test all(c -> open_fraction(t; cavity = c).fraction < 1e-5, cavs)          # 0: the holes meet the caps at shared vertices
    @test open_fraction(hohlraum(h = 1.0)).fraction < 1e-5                       # 2.7e-3 when the holes took half the count
    @test open_fraction(t; cavity = elements(t, :top)).fraction > 0.1             # without the face
    coax(r1, r2, L) = (R1 = r1 / L; R2 = r2 / L; X = 1 + (1 + R2^2) / R1^2; (X - sqrt(X^2 - 4 * (R2 / R1)^2)) / 2)
    Fl, Fe = coax(1.5, 1.0, 2.5), coax(1.5, 1.5, 2.5)
    F0 = view_factors(t)
    m = RadiationModel(t; F = copy(F0), closed = cavs)
    @test all(c -> maximum(abs.(rows(m.F, c) .- 1)) < 1e-10, cavs)
    pf = elements(t, :package_top)
    P = zeros(nelements(t))
    P[pf] .= t.elem_area[pf]
    s = solve(m, P)
    eb = energy_balance(s; by = :component)
    src = sum(P)
    # The results agree with the closed forms within the error of the view factors, and
    # nothing is lost. The raw errors are −6.0e-4, −3.0e-4 and −3.2e-5; the balanced errors
    # are −7.2e-5, +3.1e-4 and −3.5e-5.
    @test eb.escaped[:top_leh_top] / src ≈ Fl rtol = 1e-3
    @test eb.absorbed[:top_cap_top] / src ≈ Fe - Fl rtol = 1e-3
    @test eb.absorbed[:top_wall] / src ≈ 1 - Fe rtol = 1e-3
    @test abs(eb.lost) < 1e-12 * src
    @test eb.absorbed[:top_wall] + eb.absorbed[:top_cap_top] + eb.escaped[:top_leh_top] ≈ src rtol = 1e-12
    @test eb.absorbed[:bot_wall] == eb.absorbed[:package_bottom] == 0
end
