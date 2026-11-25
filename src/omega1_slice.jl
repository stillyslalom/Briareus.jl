# The transverse slice of a 1ω closest approach, as data and as a vector drawing. The slice
# is the plane across the winning beamlet through the hardware point, in which the 1ω
# rules measure.

"""
    omega1_slice(target, approach::Omega1Approach; window = nothing) -> NamedTuple

The plane across the 1ω beamlet of a closest approach ([`omega1_closest_approach`](@ref))
through its hardware point, as 2-D geometry in mm. The coordinates are `(p, q)` along `ô`
and `k × ô` of the beamlet, from the point of its 3ω axis at the station `s_star` of the
worst focal plane.

- `nominal` — the 10% contour of the 1ω footprint at that station, as `(; center, half)`
- `worst` — the same square grown by the transverse errors, the rule tolerance at the
  hardware point and the required clearance: hardware outside it passes
- `hardware` — the cross-section of the opaque triangles of `target` with the plane, as
  segments that reach into the window; a triangle lying in the plane gives its edges
- `caliper` — the hardware point and the closest point of the contour
- `distance`, `margin` — the signed transverse distance (the length of the caliper) and
  the margin of the approach
- `window` — `(pmin, pmax, qmin, qmax)`, the square about the middle of the caliper of
  half-width `window` (default four times the distance, at least 2 mm)
- `beam`, `component`, `s_star`, `focus_shift` — as in the approach
"""
function omega1_slice(t::Target, a::Omega1Approach; window = nothing)
    k, f0, ô = a.beamlet
    q̂ = cross(k, ô)
    f = f0 + a.focus_shift * k
    c = a.contour
    s = a.s_star
    b = a.budget
    grow = b.pointing + b.target_offset + b.hardware_tolerance + b.rotation_tol + b.required
    center = (c.c0 * (1 + s / c.lens_distance), 0.0)
    h = c.h0 * abs(1 - s / c.focal_shift)
    uv(x) = (dot(x - f, ô), dot(x - f, q̂))
    hw, xb = uv(a.x_hw), uv(a.x_beam)
    w = window === nothing ? max(4 * abs(a.distance), 2.0) : Float64(window)
    w > 0 || throw(ArgumentError("window must be positive"))
    mid = ((hw[1] + xb[1]) / 2, (hw[2] + xb[2]) / 2)
    box = (mid[1] - w, mid[1] + w, mid[2] - w, mid[2] + w)
    ε = 1e-9 * max(t.scale, 1.0)
    segs = NTuple{2,NTuple{2,Float64}}[]
    inbox(u, v) = max(u[1], v[1]) >= box[1] && min(u[1], v[1]) <= box[2] &&
                  max(u[2], v[2]) >= box[3] && min(u[2], v[2]) <= box[4]
    keep(x, y) = (u = uv(x); v = uv(y); inbox(u, v) && push!(segs, (u, v)))
    for j in eachindex(t.tri_area)
        (t.tri_opaque[j] && t.tri_area[j] > 0) || continue
        vs = (t.tri_v0[j], t.tri_v0[j] + t.tri_e1[j], t.tri_v0[j] + t.tri_e2[j])
        ds = map(x -> dot(x - f, k) - s, vs)
        if all(d -> abs(d) <= ε, ds)                 # in the plane: its edges
            keep(vs[1], vs[2]); keep(vs[2], vs[3]); keep(vs[3], vs[1])
            continue
        end
        (minimum(ds) > ε || maximum(ds) < -ε) && continue
        pts = Vec3[]
        for (i, l) in ((1, 2), (2, 3), (3, 1))
            di, dl = ds[i], ds[l]
            if abs(di) <= ε
                push!(pts, vs[i])
            elseif abs(dl) > ε && (di > 0) != (dl > 0)
                push!(pts, vs[i] + (di / (di - dl)) * (vs[l] - vs[i]))
            end
        end
        length(pts) == 2 && keep(pts[1], pts[2])
    end
    return (; nominal = (; center, half = h), worst = (; center, half = h + grow), hardware = segs,
            caliper = (hw, xb), distance = a.distance, margin = a.margin, window = box,
            beam = a.beam, component = a.component, s_star = s, focus_shift = a.focus_shift)
end

"""
    omega1_slice_svg(slice; size = 480) -> String

A vector drawing (SVG, `size` pixels square) of an [`omega1_slice`](@ref): the 10% contour
of the 1ω footprint as a solid outline, the worst-case square dashed over a tint, the
cross-section of the hardware in black, the caliper with its length, and a scale bar. `p`
runs to the right and `q` up. The text of the drawing is the same for the same slice:
coordinates are written with three decimals.
"""
function omega1_slice_svg(sl; size::Integer = 480)
    size >= 100 || throw(ArgumentError("size must be at least 100 pixels"))
    p0, p1, q0, q1 = sl.window
    sc = size / (p1 - p0)
    X(p) = (p - p0) * sc
    Y(q) = (q1 - q) * sc
    n(v) = @sprintf("%.3f", v)
    io = IOBuffer()
    println(io, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"", size, "\" height=\"", size,
            "\" viewBox=\"0 0 ", size, " ", size, "\" font-family=\"sans-serif\" font-size=\"12\">")
    println(io, "<rect width=\"", size, "\" height=\"", size, "\" fill=\"#ffffff\"/>")
    function square(sq, attrs)
        (cp, cq), h = sq.center, sq.half
        println(io, "<rect x=\"", n(X(cp - h)), "\" y=\"", n(Y(cq + h)), "\" width=\"", n(2h * sc),
                "\" height=\"", n(2h * sc), "\" ", attrs, "/>")
    end
    square(sl.worst, "fill=\"#D55E00\" fill-opacity=\"0.12\" stroke=\"#D55E00\" stroke-width=\"1\" stroke-dasharray=\"6 4\"")
    square(sl.nominal, "fill=\"none\" stroke=\"#D55E00\" stroke-width=\"2\"")
    line(u, v, attrs) = println(io, "<line x1=\"", n(X(u[1])), "\" y1=\"", n(Y(u[2])), "\" x2=\"", n(X(v[1])),
                                "\" y2=\"", n(Y(v[2])), "\" ", attrs, "/>")
    for (u, v) in sl.hardware
        line(u, v, "stroke=\"#000000\" stroke-width=\"1.5\"")
    end
    hw, xb = sl.caliper
    line(hw, xb, "stroke=\"#0072B2\" stroke-width=\"1.5\"")
    for u in (hw, xb)
        println(io, "<circle cx=\"", n(X(u[1])), "\" cy=\"", n(Y(u[2])), "\" r=\"3\" fill=\"#0072B2\"/>")
    end
    println(io, "<text x=\"", n(X((hw[1] + xb[1]) / 2) + 6), "\" y=\"", n(Y((hw[2] + xb[2]) / 2) - 6),
            "\" fill=\"#0072B2\">", n(sl.distance), " mm</text>")
    # scale bar: the largest of 1, 2 and 5 times a power of ten within a quarter of the view
    span = (p1 - p0) / 4
    bar = 10.0^floor(log10(span))
    bar = 5bar <= span ? 5bar : 2bar <= span ? 2bar : bar
    println(io, "<line x1=\"12.000\" y1=\"", n(size - 14.0), "\" x2=\"", n(12 + bar * sc), "\" y2=\"", n(size - 14.0),
            "\" stroke=\"#000000\" stroke-width=\"2\"/>")
    println(io, "<text x=\"12.000\" y=\"", n(size - 20.0), "\">", @sprintf("%g", bar), " mm</text>")
    println(io, "<text x=\"12.000\" y=\"18.000\">", _xml_escape(sl.beam), " / ", _xml_escape(sl.component), ": 1ω margin ",
            @sprintf("%+.3f", sl.margin), " mm at s = ", n(sl.s_star), " mm (focus ", @sprintf("%+g", sl.focus_shift),
            " mm)</text>")
    print(io, "</svg>\n")
    return String(take!(io))
end

"""
    write_omega1_slice_svg(path, target, approach; window = nothing, size = 480) -> path

Write the drawing of the slice of a 1ω closest approach ([`omega1_slice`](@ref),
[`omega1_slice_svg`](@ref)) to `path`. The drawing is plain text and is made without a graphics context.
"""
function write_omega1_slice_svg(path::AbstractString, t::Target, a::Omega1Approach; window = nothing,
                                size::Integer = 480)
    write(path, omega1_slice_svg(omega1_slice(t, a; window); size))
    return path
end
