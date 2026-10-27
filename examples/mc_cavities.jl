# Is view-factor balancing (`closed = true`) right on cavities with rims, windows and CAD
# walls? Continuum Monte Carlo references (`montecarlo_incident`, `reemit =
# :hit`) on watertight flat copies faceted at 0.1 mm, against the radiosity model on the
# default curved path, raw and balanced, centroid and sampled visibility:
#
# - the bored cavity of test/partial_visibility.jl (0.1 mm LEH bores, a 20 µm waist gap
#   closed by apertures) at 0.5 and 0.25 mm, and as flat facets at 0.5 mm, watertight
#   or cracked as the tests build it (24-sided LEH disks in 20-sided bores);
# - the same cavity from CAD-like facets (a 360-sided tessellation, test/cadfit.jl's
#   `cad_half`) through `cad_target(...; remesh, fit = true)`, at 0.5 and 0.25 mm;
# - `hohlraum()` of the same size with an open 1.5 mm window and a closed 1 × 2 mm one,
#   each in an open 50 µm gap, at h = 0.45, 0.28 and 0.17 mm.
#
# Every model shares the capsule mesh of its family (geodesic, 4,500 elements); the
# references use 16,820. Source: the CAD tutorial's two Gaussian wall rings on the
# inner wall, wall α = 0.9, η = 0.75, capsule α 0.3. The loss through the apertures is
# referenced by the same Monte Carlo with the apertures as receivers. Results are cached
# (BRIAREUS_MC_CACHE, default the temp directory); research notes, "Closure on bored,
# windowed and CAD cavities". About 30 min on 4 threads.
#
#     julia --project=. -t auto examples/mc_cavities.jl

using Briareus, Printf, LinearAlgebra, Serialization

const CACHE = get(ENV, "BRIAREUS_MC_CACHE", joinpath(tempdir(), "briareus_mc_cavities"))
mkpath(CACHE)
cached(f, key) = (p = joinpath(CACHE, "$key.jls"); isfile(p) ? deserialize(p) : (x = f(); serialize(p, x); x))

const R, Rout, Rleh, HH, g2, rc = 2.875, 2.925, 1.55, 4.7, 0.01, 0.9
const zL = HH + 0.1
const WALL = gold(albedo = 0.9, conversion = 0.75)
const CAPM = ablator()
capsule_comp(n; curved = true) = Component(:capsule, icosphere_mesh(rc; n), CAPM; surface = curved ? SphereSurface(rc) : nothing)
function openings(; curved = true, ngap = 180, nleh = 24)
    s(x) = curved ? x : nothing
    [aperture(:leh_top, disk_mesh(Rleh, zL; nr = 3, nphi = nleh, up = false); surface = s(DiskSurface(center = (0, 0, zL)))),
     aperture(:leh_bottom, disk_mesh(Rleh, -zL; nr = 3, nphi = nleh, up = true); surface = s(DiskSurface(center = (0, 0, -zL)))),
     aperture(:gap, cylinder_mesh(Rout, -g2, g2; nz = 1, nphi = ngap, inward = true); surface = s(CylinderSurface(Rout)))]
end
# `matched`: bores, LEH disks and the gap aperture share the walls' polygon (watertight facets)
function bored(hw, n; curved = true, matched = !curved)
    nphi = 8 * ceil(Int, 2π * R / hw / 8)
    nz, nr = ceil(Int, (HH - g2) / hw), ceil(Int, (R - Rleh) / hw)
    s(x) = curved ? x : nothing
    wall(z0, z1) = cylinder_mesh(R, z0, z1; nz, nphi, inward = true)
    disk(z) = s(DiskSurface(center = (0, 0, z)))
    bn = matched ? nphi : nphi ÷ 2
    Target([Component(:wall_top, wall(g2, HH), WALL; surface = s(CylinderSurface(R))),
            Component(:wall_bottom, wall(-HH, -g2), WALL; surface = s(CylinderSurface(R))),
            Component(:cap_top, annulus_mesh(Rleh, R, HH; nr, nphi, up = false), WALL; surface = disk(HH)),
            Component(:cap_bottom, annulus_mesh(Rleh, R, -HH; nr, nphi, up = true), WALL; surface = disk(-HH)),
            Component(:bore_top, cylinder_mesh(Rleh, HH, zL; nz = 1, nphi = bn), WALL; surface = s(CylinderSurface(Rleh))),
            Component(:bore_bottom, cylinder_mesh(Rleh, -zL, -HH; nz = 1, nphi = bn), WALL; surface = s(CylinderSurface(Rleh))),
            Component(:gap_top, annulus_mesh(R, Rout, g2; nr = 1, nphi, up = false), WALL; surface = disk(g2)),
            Component(:gap_bottom, annulus_mesh(R, Rout, -g2; nr = 1, nphi, up = true), WALL; surface = disk(-g2)),
            capsule_comp(n; curved), (matched ? openings(; curved, ngap = nphi, nleh = bn) : openings(; curved))...])
end

# CAD-like facets of the same cavity: each half a closed solid (test/cadfit.jl's
# `cad_half`, its waist ring at the half gap), every facet its own vertices
function cad_quad!(F, a, b, c, d, facing)
    dot(cross(b - a, c - a), facing) < 0 && ((a, b, c, d) = (d, c, b, a))
    push!(F, (a, b, c), (a, c, d))
end
cad_pt(r, φ, z) = Vec3(r * cos(φ), r * sin(φ), z)
function cad_half(; nphi, t = 0.1)
    F = NTuple{3,Vec3}[]
    φs = range(0, 2π; length = nphi + 1)
    for j in 1:nphi
        a, b = φs[j], φs[j+1]
        out = cad_pt(1, (a + b) / 2, 0)
        cad_quad!(F, cad_pt(R, a, g2), cad_pt(R, b, g2), cad_pt(R, b, HH), cad_pt(R, a, HH), -out)
        cad_quad!(F, cad_pt(Rout, a, g2), cad_pt(Rout, b, g2), cad_pt(Rout, b, HH + t), cad_pt(Rout, a, HH + t), out)
        cad_quad!(F, cad_pt(Rleh, a, HH), cad_pt(Rleh, b, HH), cad_pt(Rleh, b, HH + t), cad_pt(Rleh, a, HH + t), -out)
        cad_quad!(F, cad_pt(Rleh, a, HH), cad_pt(R, a, HH), cad_pt(R, b, HH), cad_pt(Rleh, b, HH), Vec3(0, 0, -1))
        cad_quad!(F, cad_pt(Rleh, a, HH + t), cad_pt(Rout, a, HH + t), cad_pt(Rout, b, HH + t), cad_pt(Rleh, b, HH + t),
                  Vec3(0, 0, 1))
        cad_quad!(F, cad_pt(R, a, g2), cad_pt(Rout, a, g2), cad_pt(Rout, b, g2), cad_pt(R, b, g2), Vec3(0, 0, -1))
    end
    return F
end
cad_mirror(F) = [(a .* Vec3(1, 1, -1), c .* Vec3(1, 1, -1), b .* Vec3(1, 1, -1)) for (a, b, c) in F]
facet_soup(tris) = (m = Mesh(); foreach(t -> add_polygon!(m, collect(t)), tris); m)
function cad(rm, n; nphi = 360)
    top = cad_half(; nphi)
    cad_target(Any["hohlraum_top" => facet_soup(top), "hohlraum_bottom" => facet_soup(cad_mirror(top)),
                   capsule_comp(n), openings()...]; inside = (0, 0, 0), remesh = rm, material = WALL, closed = true)
end

windows() = [WallWindow(phi = 90, z = 0.0, diameter = 1.5, open = true, gap = 0.05),
             WallWindow(phi = 270, z = 1.0, width = 1.0, height = 2.0, material = gold(albedo = 0.5, conversion = 0.75),
                        gap = 0.05)]
flat(t) = Target([Component(c.name, c.mesh, c.material; kind = c.kind, two_sided = c.two_sided) for c in t.components];
                 closed = t.closed)
windowed(h, ch; curved = true) =
    (t = hohlraum(; length = 2HH, radius = R, leh_radius = Rleh, capsule_radius = rc, h, capsule_h = ch, wall = WALL,
                  windows = windows()); curved ? t : flat(t))

# the two-ring source on the inner wall (r = R, |z| < HH), integrated over each element
source(z) = exp(-((abs(z) - 3.9) / 0.8)^2) + 0.6 * exp(-((abs(z) - 2.2) / 0.8)^2)
function is_wall(t, e)
    t.components[t.elem_comp[e]].kind === :aperture && return false
    p, n = t.elem_centroid[e], t.elem_normal[e]
    return abs(hypot(p[1], p[2]) - R) < 0.02 && abs(n[3]) < 0.5 && abs(p[3]) < HH
end
function drive(t)
    P = zeros(nelements(t))
    for e in 1:nelements(t)
        is_wall(t, e) || continue
        points, weights, _ = element_quadrature(t, e, 2)
        P[e] = sum(w * source(p[3]) for (p, w) in zip(points, weights))
    end
    return P
end
capsule(t) = findall(e -> startswith(String(t.components[t.elem_comp[e]].name), "capsule"), 1:nelements(t))
apertures(t) = findall(e -> t.components[t.elem_comp[e]].kind === :aperture, 1:nelements(t))
ratios(a) = (a[1], 100a[3] / a[1], 100a[5] / a[1])            # a₀, a₂/a₀ (%), a₄/a₀ (%)
se(x, w) = (m = sum(w .* x) / sum(w); sqrt(sum(w .^ 2 .* (x .- m) .^ 2)) / sum(w))

# references: capsule moments (per unit wall power) and the aperture loss (per unit x-ray
# source), batch by batch, over several runs (paths per receiver, seed)
const REFS = Dict("bored" => () -> bored(0.1, 29; curved = false), "windowed" => () -> windowed(0.1, 0.035; curved = false))
const RUNS = Dict("bored" => [(320, 1), (2000, 2), (2000, 3)], "windowed" => [(1500, 1), (2500, 2), (2500, 3)])
const ESC_RUNS = Dict("bored" => [(1000, 1), (25000, 2)], "windowed" => [(6000, 2)])
moments_mc(key, paths, seed) = cached("mc_$(key)_$(paths)_$(seed)") do
    t = REFS[key]()
    P, cap = drive(t), capsule(t)
    _, batch = montecarlo_incident(t, P, cap; paths, batches = 32, seed)
    (; M = (legendre_weights(t, cap; lmax = 4) * batch') ./ sum(P))
end
escape_mc(key, paths, seed) = cached("mcesc_$(key)_$(paths)_$(seed)") do
    t = REFS[key]()
    P, ap = drive(t), apertures(t)
    _, batch = montecarlo_incident(t, P, ap; paths, batches = 32, seed)
    (; esc = (batch * t.elem_area[ap]) ./ sum(conversions(t) .* P))
end
function reference(key)
    M = reduce(hcat, [moments_mc(key, p, s).M for (p, s) in RUNS[key]])
    w = reduce(vcat, [fill(Float64(p), size(moments_mc(key, p, s).M, 2)) for (p, s) in RUNS[key]])
    r = [ratios(c) for c in eachcol(M)]
    E = reduce(vcat, [escape_mc(key, p, s).esc for (p, s) in ESC_RUNS[key]])
    we = reduce(vcat, [fill(Float64(p), length(escape_mc(key, p, s).esc)) for (p, s) in ESC_RUNS[key]])
    return (; a = ratios((M * w) ./ sum(w)), err = ntuple(k -> se(getindex.(r, k), w), 3),
            esc = sum(we .* E) / sum(we), esc_err = se(E, we), open = open_fraction(REFS[key]()).fraction)
end

const MODELS = [("bored", "bored 0.5 mm", () -> bored(0.5, 15), (:centroid, :sampled)),
                ("bored", "bored 0.25 mm", () -> bored(0.25, 15), (:centroid, :sampled)),
                ("bored", "bored 0.5 mm, flat, watertight", () -> bored(0.5, 15; curved = false), (:centroid, :sampled)),
                ("bored", "bored 0.5 mm, flat, cracked", () -> bored(0.5, 15; curved = false, matched = false),
                 (:centroid, :sampled)),
                ("bored", "CAD 0.5 mm", () -> cad(0.5, 15), (:centroid, :sampled)),
                ("bored", "CAD 0.25 mm", () -> cad(0.25, 15), (:sampled,)),
                ("windowed", "windows h 0.45", () -> windowed(0.45, 0.07), (:centroid, :sampled)),
                ("windowed", "windows h 0.28", () -> windowed(0.28, 0.07), (:centroid, :sampled)),
                ("windowed", "windows h 0.17", () -> windowed(0.17, 0.07), (:centroid, :sampled))]
model(name, build, vis) = cached("model_" * replace(name, r"[^A-Za-z0-9.]+" => "_")) do
    t = build()
    P, cap = drive(t), capsule(t)
    W = legendre_weights(t, cap; lmax = 4)
    out = Dict{Any,Any}(:N => nelements(t), :open => open_fraction(t).fraction)
    for v in vis
        F = view_factors(t; visibility = v)
        ce = closure_error(F)
        ap = apertures(t)
        out[(v, :rows)] = (; mean = sum(abs, ce) / length(ce), apertures = maximum(abs, ce[ap]),
                           others = maximum(abs, ce[setdiff(1:nelements(t), ap)]))
        for closed in (false, true)
            sol = solve(RadiationModel(t; F, closed), P)
            b = energy_balance(sol)
            out[(v, closed)] = (ratios((W * sol.H[cap]) ./ sum(P)), sum(values(b.escaped)) / b.source)
        end
    end
    out
end

for key in ("bored", "windowed")
    ref = reference(key)
    @printf("\n%s reference: a₀ %.6f ± %.6f, a₂/a₀ %+.3f ± %.3f, a₄/a₀ %+.3f ± %.3f %%, loss %.5f ± %.5f (open fraction %.1e)\n",
            key, ref.a[1], ref.err[1], ref.a[2], ref.err[2], ref.a[3], ref.err[3], ref.esc, ref.esc_err, ref.open)
    println("  model − reference: a₀ and loss in %, moments in pp")
    for (k, name, build, vis) in MODELS
        k == key || continue
        m = model(name, build, vis)
        @printf("  %s (%d elements, open fraction %.1e)\n", name, m[:N], m[:open])
        for v in vis
            r = m[(v, :rows)]
            @printf("    %-8s |row − 1| mean %.1e, largest %.1e (apertures) %.1e (others)\n", v, r.mean, r.apertures, r.others)
            for closed in (false, true)
                (a0, p2, p4), loss = m[(v, closed)]
                @printf("      %-9s a₀ %+.2f  a₂/a₀ %+.3f  a₄/a₀ %+.3f  loss %+.2f\n", closed ? "balanced" : "raw",
                        100(a0 / ref.a[1] - 1), p2 - ref.a[2], p4 - ref.a[3], 100(loss / ref.esc - 1))
            end
        end
    end
end
