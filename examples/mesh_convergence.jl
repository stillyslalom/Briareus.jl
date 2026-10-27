# Work-precision study of the hohlraum mesh: capsule flux modes, mean drive and laser
# power accounting versus element size, for the NIF hohlraum of the GUI with NIF quads
# at the example cone aims ("nominal") and at a pointing with the inner cones pulled
# toward the waist and the outer cones past the axis ("tight", as an optimizer leaves
# them). Two sweeps: the wall/cap element size `h` with a fine capsule, then the
# capsule element size with a fine wall. Errors are against the finest mesh of each
# sweep. Writes a CSV (path from ARGS[1], default mesh_convergence.csv).
#
#     julia --project=. -t auto examples/mesh_convergence.jl [out.csv]

using Briareus, Printf

const GEOM = (; length = 10.1, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9)
const WALL = gold(albedo = 0.85, conversion = 0.8)

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

function measure(h, ch)
    t = hohlraum(; GEOM..., h, capsule_h = ch, wall = WALL)
    t0 = time()
    s = Scenario(t, DESIGNS[1][2]; observe = :capsule, radiation = true,
                 viewfactor_kw = (; closed = true))
    tb = time() - t0
    rows = []
    for (name, beams) in DESIGNS
        t1 = time()
        ev = evaluate(s, beams)
        te = time() - t1
        a = flux_moments(ev; lmax = 6)
        c = flux_ylm(ev; lmax = 4)
        push!(rows, (; design = name, h, capsule_h = ch, elements = nelements(t),
                     build_s = tb, eval_s = te,
                     P2 = 100a[3] / a[1], P4 = 100a[5] / a[1], P6 = 100a[7] / a[1],
                     eps4 = 100mode_amplitude(c, 4), drive = a[1],
                     wall = power_fraction(ev, :wall), caps = power_fraction(ev, :cap_top, :cap_bottom),
                     capsule = power_fraction(ev, :capsule), clipped = clipped_fraction(ev)))
    end
    return rows
end

const COLS = (:sweep, :design, :h, :capsule_h, :elements, :build_s, :eval_s, :P2, :P4, :P6, :eps4,
              :drive, :wall, :caps, :capsule, :clipped)

function main(out)
    open(out, "w") do io
        println(io, join(COLS, ","))
        sweeps = [(:wall, [(h, 0.07) for h in (0.8, 0.6, 0.45, 0.35, 0.28, 0.22, 0.17, 0.13)]),
                  (:capsule, [(0.2, ch) for ch in (0.3, 0.22, 0.16, 0.12, 0.09, 0.07, 0.05, 0.04)])]
        for (sweep, pts) in sweeps, (h, ch) in pts
            GC.gc()
            for r in measure(h, ch)
                println(io, join([sweep; [getfield(r, c) for c in COLS[2:end]]], ","))
                flush(io)
                @printf("%-8s %-8s h=%.3f ch=%.3f N=%6d build %6.1fs  P2=%+.4f P4=%+.4f drive=%.6g clip=%.3e\n",
                        sweep, r.design, h, ch, r.elements, r.build_s, r.P2, r.P4, r.drive, r.clipped)
            end
        end
    end
end

main(isempty(ARGS) ? "mesh_convergence.csv" : ARGS[1])
