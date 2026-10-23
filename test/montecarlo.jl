# Monte Carlo radiosity reference (`Briareus.montecarlo_incident`) against closed forms.

using Briareus, LinearAlgebra, Test

@testset "Monte Carlo radiosity" begin
    # A unit cube, faces turned inward, each an n × n grid of elements; the top may be open.
    face(name, nrm, c, n; open = false) =
        (m = transform(rectangle_mesh(1.0, 1.0; nu = n, nv = n), rotation_to(nrm), Vec3(c...));
         open ? aperture(name, m) : Component(name, m, gold()))
    cube(top_open, n = 1) =
        Target(face(:bottom, (0, 0, 1), (0, 0, 0), n), face(:top, (0, 0, -1), (0, 0, 1), n; open = top_open),
               face(:xm, (1, 0, 0), (-0.5, 0, 0.5), n), face(:xp, (-1, 0, 0), (0.5, 0, 0.5), n),
               face(:ym, (0, 1, 0), (0, -0.5, 0.5), n), face(:yp, (0, -1, 0), (0, 0.5, 0.5), n))
    stderr_(batch) = (nb = size(batch, 1); [sqrt(sum(abs2, c .- sum(c) / nb) / (nb - 1) / nb) for c in eachcol(batch)])

    # Closed, uniform source: B = H = ηP/A/(1 − α) on every face, in either mode.
    t = cube(false)
    for reemit in (:hit, :element)
        H, _ = Briareus.montecarlo_incident(t, ones(6), 1:6; paths = 20_000, reemit)
        @test maximum(abs, H ./ (0.75 / 0.2) .- 1) < 0.01
    end

    # Open top, source on the bottom: the element model with the closed-form view factors
    # of a cube (opposite squares (Hottel) 0.199825, adjacent ones (1 − 0.199825)/4).
    t = cube(true)
    Fo = 2 / π * (log(sqrt(4 / 3)) + 2 * sqrt(2) * atan(1 / sqrt(2)) - 2 * atan(1.0))
    @test Fo ≈ 0.199825 atol = 1e-6
    opposite = Dict(:bottom => :top, :top => :bottom, :xm => :xp, :xp => :xm, :ym => :yp, :yp => :ym)
    names = [t.components[t.elem_comp[e]].name for e in 1:6]
    F = [i == j ? 0.0 : names[j] === opposite[names[i]] ? Fo : (1 - Fo) / 4 for i in 1:6, j in 1:6]
    P = [n === :bottom ? 1.0 : 0.0 for n in names]
    H = F * ((I - Diagonal(albedos(t)) * F) \ (Briareus.conversions(t) .* P ./ t.elem_area))
    Hmc, batch = Briareus.montecarlo_incident(t, P, 1:6; paths = 64_000, reemit = :element)
    @test all(abs.(Hmc .- H) .< 5stderr_(batch))
    @test maximum(abs, Hmc ./ H .- 1) < 0.01

    # The continuum (re-emission from the hit point) is the limit of finer elements: with
    # one element per face its face averages differ from the element model's by up to 26%,
    # with 8 × 8 by under 1% (measured: 0.6%).
    t = cube(true, 8)
    P = [t.components[t.elem_comp[e]].name === :bottom ? t.elem_area[e] : 0.0 for e in 1:nelements(t)]
    sol = solve(RadiationModel(t), P)
    Hmc, _ = Briareus.montecarlo_incident(t, P, 1:nelements(t); paths = 4_000)
    face_mean(H, n) = (es = elements(t, n); sum(H[es] .* t.elem_area[es]) / sum(t.elem_area[es]))
    @test all(abs(face_mean(Hmc, n) / face_mean(sol.H, n) - 1) < 0.015 for n in (:bottom, :top, :xm, :yp))

    # Closing entrance holes: a closed lip re-emits its albedo, as in `solve`. A small flat
    # hohlraum (the meshes without their curved surfaces) driven on the wall, lips closed
    # at 3 ns, against the same hohlraum with open entrance holes: the capsule's mean
    # incident flux from the element model (sampled visibility) and from the Monte Carlo.
    # Each agrees to its view-factor error (measured 0.3–0.4%); the closure's gain (8.5%)
    # to the Monte Carlo's sampling error (5σ ≈ 0.4%).
    flat(t) = Target([Component(c.name, c.mesh, c.material; kind = c.kind, two_sided = c.two_sided)
                      for c in t.components]; closed = t.closed)
    au = Material(:Au; albedo = 0.8, conversion = 0.8)
    kw = (length = 6.0, radius = 2.0, leh_radius = 1.3, capsule_radius = 0.6, h = 0.6, capsule_h = 0.4, wall = au)
    capflux = map((nothing, GapClosure(closure = 1.0, width = 0.5, after = 0.8))) do lip
        t = flat(hohlraum(; kw..., leh_closure = lip))
        P = [t.components[t.elem_comp[e]].name === :wall ? t.elem_area[e] : 0.0 for e in 1:nelements(t)]
        cap = elements(t, :capsule)
        A = t.elem_area[cap] ./ sum(t.elem_area[cap])
        sol = solve(RadiationModel(t; F = view_factors(t; visibility = :sampled)), P; time = 3.0)
        lip === nothing || @test all(sol.albedo[elements(t, :leh_bottom_lip)] .≈ 0.8)
        _, batch = Briareus.montecarlo_incident(t, P, cap; paths = 8_000, reemit = :element,
                                                albedo = albedos(t, 3.0))
        means = batch * A
        nb = length(means)
        (; sol = sum(sol.H[cap] .* A), mc = sum(means) / nb,
           se = sqrt(sum(abs2, means .- sum(means) / nb) / (nb - 1) / nb))
    end
    (o, c) = capflux
    @test abs(o.mc / o.sol - 1) < 0.006 && abs(c.mc / c.sol - 1) < 0.006
    gain_mc, gain = c.mc / o.mc, c.sol / o.sol
    @test gain > 1.05
    @test abs(gain_mc - gain) < 5gain_mc * hypot(o.se / o.mc, c.se / c.mc)
end
