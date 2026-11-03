# Frustraum cone pointing, derived from the example cylinder pointing (BriareusGUI's
# FRUSTRAUM_CONE_POINTING): each cone's spot centroid on the frustraum wall at the
# fraction of the half-length where the cylinder pointing puts it on the 10.1 mm
# cylinder, where the entrance hole allows it. An engineering default, derived; no
# frustraum pointing beyond the 50° cone's is public (research notes, "Frustraum cone
# pointing: derived (2026-11-03)").
#
#   julia --project=. -t auto examples/frustraum_pointing.jl
#
# Method. A cone's beam axis is fixed by its port; its line through the entrance-hole
# plane is set by the crossing r0 (toward the port azimuth), and its aim (best focus)
# keeps the cylinder aim's axial distance from the hole. The spot centroid is the
# power-weighted mean height of the wall deposition of the cone's upper-hemisphere quads,
# over the half-length. r0 is solved by bisection for the cylinder's fraction; where the
# beams' keep-out clearance from the hole's rim (`ApertureKeepOut` defaults) would fall
# below the cylinder pointing's smallest (its 23.5° cone's), r0 stops at that clearance
# instead, on the side toward the cylinder's fraction. Aims are rounded to 0.01 mm and
# the result checked by ray trace (nspot = 121).

using Briareus, Printf

const CONES = (:inner23, :inner30, :outer44, :outer50)
# the example cylinder pointing (z, r; mm), research notes "Example NIF cone pointing"
const CYL = Dict(:inner23 => (3.09, -1.09), :inner30 => (3.63, -0.60), :outer44 => (5.20, 0.0),
                 :outer50 => (5.01, 0.0))
const Lc, Lf = 10.1, 12.24

cyl = hohlraum(length = Lc, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9, h = 0.2)
fr = hohlraum(profile = frustraum_profile(length = Lf, waist_radius = 4.88, end_radius = 2.282),
              leh_radius = 1.55, capsule_radius = 1.281, h = 0.2)

upper = [b for b in nif_quads(power = 8.0) if b.port[3] > 0]
tanθ = Dict(c => (b = first(b for b in upper if b.group === c); hypot(b.port[1], b.port[2]) / b.port[3])
            for c in CONES)

function pointed(tab)
    p = GroupPointing(upper; fields = (:z, :r))
    return apply(p, parameters(p; (k => (z = tab[k][1], r = tab[k][2]) for k in CONES)...))
end

# a cone's wall-deposition centroid height over the half-length, and its wall share
function fraction(t, deposition, beams, cone, halfL)
    d = deposition([remake(b; power = b.group === cone ? b.power : 0.0) for b in beams])
    w = elements(t, :wall)
    P = d.element[w]
    return sum(P .* [t.elem_centroid[e][3] for e in w]) / sum(P) / halfL,
           sum(P) / sum(b.power for b in beams if b.group === cone)
end

kc, kf = ApertureKeepOut(cyl), ApertureKeepOut(fr)
clearance_of(k, bs, c) = minimum(aperture_clearance(k, [b for b in bs if b.group === c]))
opc, opf = DepositionOperator(cyl, upper), DepositionOperator(fr, upper)

bc = pointed(CYL)
want = Dict(c => fraction(cyl, b -> deposit(opc, b), bc, c, Lc / 2)[1] for c in CONES)
margin = minimum(clearance_of(kc, bc, c) for c in CONES)
@printf "cylinder pointing's smallest clearance: %.4f mm\n" margin

# the aim on the frustraum from the hole-plane crossing r0
aim_of(c, r0) = (za = Lf / 2 + (CYL[c][1] - Lc / 2); (za, r0 - (Lf / 2 - za) * tanθ[c]))
function bisect(g, a, b)
    ga = g(a)
    while b - a > 1e-6
        m = (a + b) / 2
        gm = g(m)
        sign(gm) == sign(ga) ? (a = m; ga = gm) : (b = m)
    end
    return (a + b) / 2
end
tab = Dict(c => aim_of(c, CYL[c][2] + (Lc / 2 - CYL[c][1]) * tanθ[c]) for c in CONES)
for c in CONES
    with(r0) = pointed(merge(tab, Dict(c => aim_of(c, r0))))
    r0 = bisect(r -> fraction(fr, b -> deposit(opf, b), with(r), c, Lf / 2)[1] - want[c], -1.5, 1.0)
    if clearance_of(kf, with(r0), c) < margin
        r0 = bisect(r -> clearance_of(kf, with(r), c) - margin, -1.0, 0.0)
    end
    za, ra = aim_of(c, r0)
    tab[c] = (round(za; digits = 2), round(ra; digits = 2))
end

println("cone     aim z    aim r   fraction: cylinder  frustraum (ray trace, nspot 121)  clearance (mm)")
bf = pointed(tab)
for c in CONES
    fc = fraction(cyl, b -> raytrace_deposition(cyl, b; nspot = 121), [b for b in bc if b.group === c], c, Lc / 2)
    ff = fraction(fr, b -> raytrace_deposition(fr, b; nspot = 121), [b for b in bf if b.group === c], c, Lf / 2)
    @printf "%-8s %6.2f  %6.2f   %8.4f  %8.4f (wall share %.4f)   %.4f\n" c tab[c]... fc[1] ff[1] ff[2] clearance_of(kf, bf, c)
end
