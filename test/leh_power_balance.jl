# X-ray power balance of the N210808 (HYBRID-E) hohlraum against LEH size: fractions of
# the x-ray source escaping through the LEHs and absorbed by the capsule, against a
# Monte Carlo photon walk in the exact geometry (docs/literate/leh_power_balance.jl,
# `monte_carlo` at 10⁸ histories, seed 7, σ 3.5e-5 / 2.8e-5). The source is piecewise
# constant on the h = 0.4 wall mesh, so the references hold for that mesh.

using Briareus, Test

@testset "LEH power balance (N210808)" begin
    αw, αc = albedo(WallLossAlbedo(), 6.0, 300.0), 0.3
    R, L, rc = 3.2, 11.24, 1.13
    reference = Dict(3.0 => (leh = 0.32826, capsule = 0.24170),
                     3.1 => (leh = 0.34421, capsule = 0.23675),
                     3.64 => (leh = 0.42821, capsule = 0.21086))
    # inner cones (33%) at 0.6 ≤ |z| ≤ 1.8 mm, outer cones at 3.4 ≤ |z| ≤ 4.6 mm
    bands = ((0.33, 0.6, 1.8), (0.67, 3.4, 4.6))
    for (d, ref) in sort(collect(reference))
        t = hohlraum(; length = L, radius = R, leh_radius = d / 2, capsule_radius = rc,
                     h = 0.4, capsule_h = 0.4 / 3,
                     wall = Material(:Au; albedo = αw, conversion = 1.0),
                     capsule = Material(:ablator; albedo = αc, conversion = 0.0))
        P = zeros(nelements(t))
        for e in elements(t, :wall)
            pts = [p for k in t.elem_tris[e] for p in (t.tri_v0[k], t.tri_v0[k] + t.tri_e1[k],
                                                        t.tri_v0[k] + t.tri_e2[k])]
            z0, z1 = extrema(p[3] for p in pts)
            for (f, a, b) in bands, (lo, hi) in ((a, b), (-b, -a))
                P[e] += f / 2 * max(0.0, min(z1, hi) - max(z0, lo)) / (b - a) *
                        t.elem_area[e] / (2π * R * (z1 - z0))
            end
        end
        eb = energy_balance(solve(RadiationModel(t; closed = true), P))
        leh = (eb.escaped[:leh_top] + eb.escaped[:leh_bottom]) / eb.source
        @test leh ≈ ref.leh atol = 3e-4
        @test eb.absorbed[:capsule] / eb.source ≈ ref.capsule atol = 3e-4
        # the uniform-flux power balance reads the LEH loss 6.5–7% low
        Aw = 2π * R * L + 2π * (R^2 - d^2 / 4)
        AL = π * d^2 / 2
        @test 0.06 < 1 - AL / ((1 - αw) * Aw + (1 - αc) * 4π * rc^2 + AL) / leh < 0.08
    end
end
