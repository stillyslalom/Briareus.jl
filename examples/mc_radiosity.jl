# Partial visibility against a Monte Carlo radiosity: the capsule's incident flux, a₂/a₀
# and a₄/a₀ from the radiosity solve with centroid and sampled visibility
# (`view_factors(...; visibility = :sampled)`), raw and balanced (`closed`), against
# `Briareus.montecarlo_incident`, which traces diffuse transport on the same facets with
# no view factors. Targets: the bored-LEH cavity of test/partial_visibility.jl (0.1 mm
# LEH bores, 20 µm waist gap) at 0.5 and 0.25 mm, and the flush-LEH `hohlraum()` as a
# control. The source is the CAD tutorial's two Gaussian wall rings, integrated over
# each element. The comparison uses flat elements (the meshes without their curved
# surfaces), the geometry the Monte Carlo traces; the curved models are listed for
# reference. Results are cached in the temp directory (a rerun only prints them).
#
#     julia --project=. -t auto examples/mc_radiosity.jl [paths per capsule element]

using Briareus, Printf, Serialization

const PATHS = isempty(ARGS) ? 100_000 : parse(Int, ARGS[1])
const CACHE = joinpath(tempdir(), "briareus_mc_radiosity")
mkpath(CACHE)

# The bored cavity of test/partial_visibility.jl; `curved = false` drops the surfaces.
function bored_cavity(; hw = 0.5, hc = 0.15, curved = true)
    R, Rout, Rleh, h, g2, rc = 2.875, 2.925, 1.55, 4.7, 0.01, 0.9
    zL = h + 0.1
    nphi = 8 * ceil(Int, 2π * R / hw / 8)
    nz, nr = ceil(Int, (h - g2) / hw), ceil(Int, (R - Rleh) / hw)
    (; capsule_ntheta, capsule_nphi) = capsule_counts(rc, hc)
    au = gold()
    s(x) = curved ? x : nothing
    wall(z0, z1) = cylinder_mesh(R, z0, z1; nz, nphi, inward = true)
    disk(z) = s(DiskSurface(center = (0, 0, z)))
    return Target([Component(:wall_top, wall(g2, h), au; surface = s(CylinderSurface(R))),
                   Component(:wall_bottom, wall(-h, -g2), au; surface = s(CylinderSurface(R))),
                   Component(:cap_top, annulus_mesh(Rleh, R, h; nr, nphi, up = false), au; surface = disk(h)),
                   Component(:cap_bottom, annulus_mesh(Rleh, R, -h; nr, nphi, up = true), au; surface = disk(-h)),
                   Component(:bore_top, cylinder_mesh(Rleh, h, zL; nz = 1, nphi = nphi ÷ 2), au;
                             surface = s(CylinderSurface(Rleh))),
                   Component(:bore_bottom, cylinder_mesh(Rleh, -zL, -h; nz = 1, nphi = nphi ÷ 2), au;
                             surface = s(CylinderSurface(Rleh))),
                   Component(:gap_top, annulus_mesh(R, Rout, g2; nr = 1, nphi, up = false), au; surface = disk(g2)),
                   Component(:gap_bottom, annulus_mesh(R, Rout, -g2; nr = 1, nphi, up = true), au; surface = disk(-g2)),
                   Component(:capsule, sphere_mesh(rc; ntheta = capsule_ntheta, nphi = capsule_nphi), ablator();
                             surface = s(SphereSurface(rc))),
                   aperture(:leh_top, disk_mesh(Rleh, zL; nr = 3, nphi = 24, up = false); surface = disk(zL)),
                   aperture(:leh_bottom, disk_mesh(Rleh, -zL; nr = 3, nphi = 24, up = true); surface = disk(-zL)),
                   aperture(:gap, cylinder_mesh(Rout, -g2, g2; nz = 1, nphi = 180, inward = true);
                            surface = s(CylinderSurface(Rout)))])
end

flat(t::Target) = Target([Component(c.name, c.mesh, c.material; kind = c.kind, two_sided = c.two_sided)
                          for c in t.components]; closed = t.closed)

source(z) = exp(-((abs(z) - 3.9) / 0.8)^2) + 0.6 * exp(-((abs(z) - 2.2) / 0.8)^2)   # TW/mm²
function drive(t)
    P = zeros(nelements(t))
    for e in 1:nelements(t)
        startswith(String(t.components[t.elem_comp[e]].name), "wall") || continue
        points, weights, _ = element_quadrature(t, e, 2)
        P[e] = sum(w * source(p[3]) for (p, w) in zip(points, weights))
    end
    return P
end

ratios(a) = (a[1], 100a[3] / a[1], 100a[5] / a[1])     # a₀, a₂/a₀ (%), a₄/a₀ (%)

function radiosity(t)
    cap = elements(t, :capsule)
    W = legendre_weights(t, cap; lmax = 4)
    P = drive(t)
    rows = Pair{String,Any}[]
    for vis in (:centroid, :sampled)
        sec = @elapsed F = view_factors(t; visibility = vis)
        for closed in (false, true)
            sol = solve(RadiationModel(t; F, closed), P)
            push!(rows, "$vis, closed = $closed" => (ratios(W * sol.H[cap])..., sec))
        end
    end
    return rows
end

function montecarlo(t, reemit)
    cap = elements(t, :capsule)
    W = legendre_weights(t, cap; lmax = 4)
    sec = @elapsed H, batch = Briareus.montecarlo_incident(t, drive(t), cap; paths = PATHS, reemit, batches = 64)
    r = [ratios(W * b) for b in eachrow(batch)]
    nb = length(r)
    err(k) = (x = getindex.(r, k); sqrt(sum(abs2, x .- sum(x) / nb) / (nb - 1) / nb))
    return (ratios(W * H)..., err(1), err(2), err(3), sec)
end

cached(f, key) = (path = joinpath(CACHE, "$key.jls");
                  isfile(path) ? deserialize(path) : (x = f(); serialize(path, x); x))

cases = [("bored, 0.5 mm", () -> bored_cavity(hw = 0.5, hc = 0.15, curved = false),
          () -> bored_cavity(hw = 0.5, hc = 0.15)),
         ("bored, 0.25 mm", () -> bored_cavity(hw = 0.25, hc = 0.1, curved = false),
          () -> bored_cavity(hw = 0.25, hc = 0.1)),
         ("hohlraum()", () -> flat(hohlraum()), () -> hohlraum())]

for (name, flat_target, curved_target) in cases
    key = replace(name, r"[^a-z0-9]+" => "_")
    t = flat_target()
    println("\n", name, ": ", nelements(t), " elements, ", length(elements(t, :capsule)), " on the capsule")
    for reemit in (:hit, :element)
        a0, p2, p4, e0, e2, e4, sec = cached(() -> montecarlo(t, reemit), "$(key)_mc_$(reemit)_$PATHS")
        @printf("  %-34s a₀ %.5f ± %.5f   a₂/a₀ %+.3f ± %.3f   a₄/a₀ %+.3f ± %.3f   (%.0f s)\n",
                "Monte Carlo, reemit = :$reemit", a0, e0, p2, e2, p4, e4, sec)
    end
    for (label, target) in (("flat", flat_target), ("curved", curved_target))
        for (row, (a0, p2, p4, sec)) in cached(() -> radiosity(target()), "$(key)_$label")
            @printf("  %-34s a₀ %.5f             a₂/a₀ %+.3f           a₄/a₀ %+.3f           (F %.1f s)\n",
                    "$label, $row", a0, p2, p4, sec)
        end
    end
end
