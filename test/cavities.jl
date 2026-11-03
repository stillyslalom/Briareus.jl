# Targets with several cavities (assemblies phase C2): each cavity's wall temperature
# (`wall_temperature`, `drive_history(; cavities)`) against the area-weighted incident
# flux of its wall elements reduced here directly from the solution, the history without
# `cavities` unchanged bit for bit, and a Dante field of view that leaves one cavity out
# receiving nothing from it.

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

    # independent reduction: each cavity's wall elements found from the components'
    # names, Σ H A / Σ A of each step's recorded incident flux
    own(p) = [e for e in 1:nelements(t) if t.components[t.elem_comp[e]].name in walls(p)]
    Tr(H, es) = (sum(H[e] * t.elem_area[e] for e in es) / sum(t.elem_area[e] for e in es) / SIGMA_SB)^0.25
    @test all(isapprox(h.cavity_Tr[p], [Tr(H, own(p)) for H in h.record]; rtol = 1e-12) for p in (:a, :b))
    # two cavities, two temperatures; together they are the history's wall temperature
    @test all(abs.(h.cavity_Tr[:a] .- h.cavity_Tr[:b]) .> 1.0)
    @test h.cavity_Tr[:all] == h.wall_Tr
    # without `cavities` the history is the same, bit for bit
    @test isempty(h0.cavity_Tr) && all(k -> k === :cavity_Tr || isequal(h0[k], h[k]), keys(h))
    # `wall_temperature` on one solution
    sol = solve(m, deposit(op, q))
    @test wall_temperature(sol, walls(:b)...) ≈ Tr(sol.H, own(:b)) rtol = 1e-12
    # `per_component`: the energy per piece, not summed into the parts; each part's sum of
    # its pieces is the by-part history's, and every other field is unchanged
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
    n = nelements(t)
    T0 = 250.0
    lit(p) = [c.kind !== :aperture && c.part === p ? SIGMA_SB * T0^4 : 0.0 for c in t.components[t.elem_comp]]
    z = zeros(n)
    bb(B) = RadiationSolution(RadiationModel(t; F = zeros(n, n)), z, z, B, z, z, z, z, 0.0)
    # a field of 8 mm about a's center leaves b out: b gives it nothing
    narrow = DanteView(143, 274; through = :a_leh_bottom, field_of_view = 8.0, center = :a, pixels = 400)
    @test dante(narrow, bb(lit(:b))).intensity == 0 && !haskey(visible_area(narrow, t), :b)
    # the views through each hole (no field of view), and a 60 mm field at TCC (a NIF
    # Dante's, research notes) that takes in both cavities
    hole(p) = dante(DanteView(143, 274; through = Symbol(p, :_leh_bottom)), bb(lit(p))).intensity
    wide = DanteView(143, 274; through = :a_leh_bottom, field_of_view = 60.0, pixels = 1200)
    @test dante(narrow, bb(lit(:a))).intensity ≈ hole(:a) rtol = 1e-3       # measured 1.2e-4
    @test dante(wide, bb(lit(:a) .+ lit(:b))).intensity ≈ hole(:a) + hole(:b) rtol = 1e-3   # 1.4e-4
    @test issubset((:a, :b), keys(visible_area(wide, t)))
    # the built-in Dante 1 sees the whole target: that field at TCC, Dante 1's channels;
    # this target has no `leh_bottom` to quote Tr over, so only the power
    fv = dante1()
    @test fv.field_of_view == 60.0 && fv.center == zero(fv.center) && fv.through === :leh_bottom &&
          (fv.theta, fv.phi) == (143.0, 274.0) && length(fv.channels) == 18
    r = dante(fv, bb(lit(:a) .+ lit(:b)))
    @test isnan(r.Tr) && isapprox(r.intensity, hole(:a) + hole(:b); rtol = 5e-3)   # 1.8e-3 (400 pixels)
end
