# Capsule mesh study: Legendre moments of the capsule flux versus capsule element count
# for the latitude-longitude capsule (`hohlraum(capsule_mesh = :latlong)`) and the
# geodesic one (`:compact`, `icosphere_mesh`, the default since 2026-10-22), on the hohlraum and the two
# pointings of examples/mesh_convergence.jl with the wall fixed at h = 0.2 mm. P4 is
# extrapolated to zero element size from the two finest compact meshes (error ∝ h²).
# Writes a CSV (path from ARGS[1], default capsule_mesh.csv); a rerun skips rows the
# CSV already has.
#
#     julia --project=. -t auto examples/capsule_mesh.jl [out.csv]

using Briareus, Printf

const GEOM = (; length = 10.1, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9)
const WALL = gold(albedo = 0.85, conversion = 0.8)
const WALL_H = 0.2

function pointed(rz)
    bs = nif_quads(power = 8.0)
    p = GroupPointing(bs; by = :group, fields = (:r, :z))
    x = initial_parameters(p)
    for (g, k) in enumerate(p.keys)
        x[2g-1], x[2g] = rz[k]
    end
    return apply(p, x)
end

const DESIGNS = [
    :nominal => pointed(Dict(:inner23 => (-1.09, 3.09), :inner30 => (-0.6, 3.63),
                             :outer44 => (0.0, 5.2), :outer50 => (0.0, 5.01))),
    :tight => pointed(Dict(:inner23 => (-0.79, 2.09), :inner30 => (0.4, 4.40),
                           :outer44 => (-1.0, 4.36), :outer50 => (-0.93, 4.01))),
]

function measure(kind, ch)
    t = hohlraum(; GEOM..., h = WALL_H, capsule_h = ch, capsule_mesh = kind, wall = WALL)
    ncap = length(elements(t, :capsule))
    t0 = time()
    s = Scenario(t, DESIGNS[1][2]; observe = :capsule, radiation = true,
                 viewfactor_kw = (; closed = true))
    tb = time() - t0
    rows = []
    for (name, beams) in DESIGNS
        ev = evaluate(s, beams)
        a = flux_moments(ev; lmax = 6)
        push!(rows, (; mesh = kind, design = name, capsule_h = ch, capsule = ncap,
                     elements = nelements(t), build_s = tb,
                     P2 = 100a[3] / a[1], P4 = 100a[5] / a[1], P6 = 100a[7] / a[1], drive = a[1]))
    end
    return rows
end

const COLS = (:mesh, :design, :capsule_h, :capsule, :elements, :build_s, :P2, :P4, :P6, :drive)
const SWEEP = [(:latlong, ch) for ch in (0.3, 0.22, 0.16, 0.12, 0.09, 0.07, 0.05, 0.04)] ∪
              [(:compact, ch) for ch in (0.3, 0.22, 0.16, 0.12, 0.09, 0.07, 0.05, 0.04)]

function main(out)
    done = Set{Tuple{String,String}}()
    if isfile(out)
        for l in readlines(out)[2:end]
            f = split(l, ",")
            push!(done, (f[1], f[3]))
        end
    else
        open(io -> println(io, join(COLS, ",")), out, "w")
    end
    for (kind, ch) in SWEEP
        (string(kind), string(ch)) in done && continue
        GC.gc()
        rows = measure(kind, ch)
        open(out, "a") do io
            for r in rows
                println(io, join([getfield(r, c) for c in COLS], ","))
                @printf("%-8s %-8s ch=%.3f Ncap=%6d N=%6d build %6.1fs  P2=%+.4f P4=%+.4f P6=%+.4f drive=%.7g\n",
                        r.mesh, r.design, ch, r.capsule, r.elements, r.build_s, r.P2, r.P4, r.P6, r.drive)
            end
        end
    end
end

main(isempty(ARGS) ? "capsule_mesh.csv" : ARGS[1])
