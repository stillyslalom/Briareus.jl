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
end
