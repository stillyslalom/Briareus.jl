# Per-cavity balancing (`RadiationModel(...; closed = cavities)`) on two halfraums joined
# by a package: each halfraum (2.5 mm long, radius 1.5 mm, an entrance hole of 1 mm
# radius at its outer end) and the package face closing its open end form a cavity; the
# package's outer face (absorbing) is in neither. Walls albedo 0.8, package faces 0.3;
# the upper wall emits 1 per unit area, the lower 0.5. Area-mean incident flux on the
# package faces and the walls:
#
# - the default curved path, raw and balanced, at h = 0.4, 0.28, 0.2 and 0.14 mm, against
#   the balanced model at 0.14 mm (convergence in h);
# - flat facets, watertight (the holes' disks of the caps' full azimuthal count), raw and
#   balanced, against the Monte Carlo of the element model with exact view factors
#   (`montecarlo_incident(...; reemit = :element)`, 200,000 paths per receiver), at
#   h = 0.4 and 0.28 mm.
#
# Research notes, "Per-cavity balancing". About 5 minutes on 4 threads.

using Briareus, LinearAlgebra, Printf, Statistics

named(c, p; flat = false) = Component(Symbol(p, :_, c.name), c.mesh, c.material; kind = c.kind,
                                      two_sided = c.two_sided, surface = flat ? nothing : c.surface, part = p)

function pair(h; flat = false)
    au = Material(:au; albedo = 0.8, conversion = 1.0)
    kw = (; length = 2.5, radius = 1.5, leh_radius = 1.0, capsule_radius = 0.0, h, wall = au)
    top = hohlraum(; kw..., bottom_end = :open)
    bot = hohlraum(; kw..., top_end = :open)
    nphi = hohlraum_counts(; length = 2.5, radius = 1.5, leh_radius = 1.0, h).nphi
    nz, nr = max(2, ceil(Int, 2.0 / h)), max(1, ceil(Int, 1.5 / h))
    face = Material(:pk; albedo = 0.3, conversion = 0.0)
    sf(s) = flat ? nothing : s
    pkg = [Component(:package, cylinder_mesh(1.5, -1.0, 1.0; nz, nphi, inward = false),
                     Material(:k; albedo = 0.0, conversion = 0.0); surface = sf(CylinderSurface(1.5))),
           Component(:package_top, disk_mesh(1.5, 1.0; nr, nphi, up = true), face;
                     surface = sf(DiskSurface(center = (0, 0, 1.0)))),
           Component(:package_bottom, disk_mesh(1.5, -1.0; nr, nphi, up = false), face;
                     surface = sf(DiskSurface(center = (0, 0, -1.0))))]
    # flat facets watertight: each hole's disk takes the cap's full azimuthal count (the
    # hohlraum's has half, which leaves cracks that only the curved path closes)
    tight(c) = !flat || c.kind !== :aperture ? c :
               aperture(c.name, disk_mesh(1.0, c.mesh.vertices[1][3]; nr = 2, nphi, up = c.mesh.vertices[1][3] < 0))
    t = Target([[translate(named(tight(x), :top; flat), (0, 0, 2.25)) for x in top.components];
                [translate(named(tight(x), :bot; flat), (0, 0, -2.25)) for x in bot.components]; pkg])
    cavs = [[elements(t, :top); elements(t, :package_top)], [elements(t, :bot); elements(t, :package_bottom)]]
    return t, cavs
end

function drive(t)
    P = zeros(nelements(t))
    P[elements(t, :top_wall)] .= t.elem_area[elements(t, :top_wall)]
    P[elements(t, :bot_wall)] .= 0.5 .* t.elem_area[elements(t, :bot_wall)]
    return P
end
areamean(t, H, n) = (es = elements(t, n); sum(H[es] .* t.elem_area[es]) / sum(t.elem_area[es]))
qs = (:package_top, :package_bottom, :top_wall, :bot_wall)

res = Dict()
for h in (0.4, 0.28, 0.2, 0.14)
    t, cavs = pair(h)
    F0 = view_factors(t)
    P = drive(t)
    for (lab, closed) in (("raw", false), ("balanced", cavs))
        s = solve(RadiationModel(t; F = copy(F0), closed), P)
        res[(h, lab)] = [areamean(t, s.H, q) for q in qs]
    end
    @printf("curved h = %.2f (%d elements): raw %s, balanced %s\n", h, nelements(t),
            string(round.(res[(h, "raw")]; sigdigits = 7)), string(round.(res[(h, "balanced")]; sigdigits = 7)))
end
ref = res[(0.14, "balanced")]
for h in (0.4, 0.28, 0.2, 0.14), lab in ("raw", "balanced")
    @printf("curved h = %.2f %-8s − balanced at 0.14 (relative; %s): %s\n", h, lab, join(qs, ", "),
            string(round.(res[(h, lab)] ./ ref .- 1; sigdigits = 3)))
end

for h in (0.4, 0.28)
    t, cavs = pair(h; flat = true)
    @printf("flat h = %.2f: open fractions of the cavities %s\n", h, string([open_fraction(t; cavity = c).fraction for c in cavs]))
    F0 = view_factors(t)
    P = drive(t)
    rec = [elements(t, :package_top); elements(t, :package_bottom)]
    nt = length(elements(t, :package_top))
    A = t.elem_area[rec]
    faces(H) = [sum(H[1:nt] .* A[1:nt]) / sum(A[1:nt]), sum(H[nt+1:end] .* A[nt+1:end]) / sum(A[nt+1:end])]
    Hm, batch = montecarlo_incident(t, P, rec; paths = 200_000, batches = 40, reemit = :element)
    mc = faces(Hm)
    σ = [std([faces(batch[b, :])[k] for b in axes(batch, 1)]) / sqrt(size(batch, 1)) for k in 1:2]
    @printf("  Monte Carlo, package faces: %s ± %s\n", string(round.(mc; sigdigits = 7)), string(round.(σ; sigdigits = 2)))
    for (lab, closed) in (("raw", false), ("balanced", cavs))
        s = solve(RadiationModel(t; F = copy(F0), closed), P)
        v = [areamean(t, s.H, :package_top), areamean(t, s.H, :package_bottom)]
        @printf("  %-8s − Monte Carlo (relative): %s\n", lab, string(round.(v ./ mc .- 1; sigdigits = 3)))
    end
end
