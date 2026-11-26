# The chamber port map: the ports of the beams of a shot on two polar disks, one per
# hemisphere, as data and as a vector drawing colored by a value or a tier.

const _PORT_TIER_COLORS = Dict(:green => "#009E73", :amber => "#E69F00", :red => "#D55E00")
const _PORT_NEUTRAL = "#BBBBBB"

"""
    port_map(beams; by = :auto) -> Vector{NamedTuple}

The marks of a chamber port map, one per beam (`by = :beam`) or per quad (`by = :quad`,
the beams that share a `quad` name, at the mean of their port directions). `by = :auto`
takes quads when two beams share a quad name. The ports come from the beams, so any
facility maps. Each mark has:

- `name`, and `beams`, the indices of its beams;
- `theta` and `phi`, the polar angle and azimuth of the port (degrees);
- `hemisphere`, `:upper` for `theta ≤ 90` and `:lower` otherwise;
- `x` and `y`, the position on the unit disk of its hemisphere: the pole is at the
  center and the equator on the rim, at radius `theta / 90` (upper) or
  `(180 − theta) / 90` (lower), with the azimuth counterclockwise from `+x` on both
  disks. The lower disk is the lower hemisphere seen from above, through the chamber.
"""
function port_map(beams::AbstractVector{<:Beam}; by::Symbol = :auto)
    by in (:auto, :beam, :quad) || throw(ArgumentError("by must be :auto, :beam or :quad"))
    quads = String[]
    members = Vector{Int}[]
    for (i, b) in enumerate(beams)
        isempty(b.quad) && continue
        j = findfirst(==(b.quad), quads)
        j === nothing ? (push!(quads, b.quad); push!(members, [i])) : push!(members[j], i)
    end
    if by === :auto
        by = any(m -> length(m) > 1, members) ? :quad : :beam
    end
    if by === :quad
        sum(length, members; init = 0) == length(beams) ||
            throw(ArgumentError("by = :quad needs a quad name on every beam"))
        return [_port_mark(q, m, beams) for (q, m) in zip(quads, members)]
    end
    return [_port_mark(b.name, [i], beams) for (i, b) in enumerate(beams)]
end

function _port_mark(name::AbstractString, members::Vector{Int}, beams)
    p = sum(Vec3(beams[i].port) for i in members)
    n = norm(p)
    n > 1e-9 || throw(ArgumentError("the ports of $name cancel: no mean direction"))
    p /= n
    theta = acosd(clamp(p[3], -1.0, 1.0))
    phi = mod(atand(p[2], p[1]), 360.0)
    hemisphere = theta <= 90 + 1e-9 ? :upper : :lower
    r = hemisphere === :upper ? theta / 90 : (180 - theta) / 90
    return (; name = String(name), beams = members, theta, phi, hemisphere, x = r * cosd(phi), y = r * sind(phi))
end

# A value per mark from a vector in the order of the marks or a dictionary by name.
function _port_lookup(v, marks, what)
    v === nothing && return nothing
    if v isa AbstractDict
        return [get(v, m.name, nothing) for m in marks]
    end
    length(v) == length(marks) || throw(ArgumentError("$what: $(length(v)) entries for $(length(marks)) marks"))
    return collect(Any, v)
end

# A color between light and dark blue for `x` in [0, 1].
function _port_ramp(x::Real)
    a, b = (222, 235, 247), (8, 48, 107)
    c = ntuple(i -> round(Int, a[i] + clamp(x, 0, 1) * (b[i] - a[i])), 3)
    return @sprintf("#%02X%02X%02X", c...)
end

"""
    port_map_svg(beams; by = :auto, tiers = nothing, values = nothing, notes = nothing,
                 badges = nothing, roles = nothing, title = "", unit = "", size = 360,
                 labels = :auto) -> String
    port_map_svg(report::ROEReport, beams; kind = :auto, badges = :auto, kwargs...) -> String

A vector drawing (SVG) of the chamber port map ([`port_map`](@ref)): the upper hemisphere
on the left disk and the lower hemisphere, seen from above, on the right disk, with rings
at polar angles of 30° and 60° from each pole and the equator on the rim.

`tiers` colors each mark green, amber or red (`:green`, `:amber`, `:red`; any other value
is gray), and `values` colors it on a ramp between the smallest and largest value, with
`unit` in the legend. Both are a vector in the order of the marks or a dictionary by mark
name, and a mark without an entry is gray. `notes`, in the same form, is text shown when
the pointer rests on a mark. `badges`, in the same form, is a short text written at the
upper right of each mark, such as the role of its beams, and `roles` a dictionary from
badge to the role that it stands for, listed under the disks. `labels` writes the name
under each mark (`:auto`: for 80 marks or fewer). `size` is the diameter of a disk in
pixels.

With a rules-of-engagement report ([`ROEReport`](@ref)) and the beams of its scan, the
marks are the units of the report (`kind = :quad` or `:beam`; `:auto` takes quads when the
report has quad units), colored by tier, and each note gives the role, the decision and
the binding criterion with its score. With `badges = :auto`, each mark carries the badge
of its role ([`role_badges`](@ref)) when the marks have more than one role, as for a drive
and a backlighter, two halfraums or two foils; `true` writes them always and `false`
never. The text of the drawing is the same for the same input: coordinates are written
with two decimals.
"""
function port_map_svg(beams::AbstractVector{<:Beam}; by::Symbol = :auto, kwargs...)
    return _port_map_svg(port_map(beams; by); kwargs...)
end

function port_map_svg(r::ROEReport, beams::AbstractVector{<:Beam}; kind::Symbol = :auto, badges = :auto, kwargs...)
    kind in (:auto, :beam, :quad) || throw(ArgumentError("kind must be :auto, :beam or :quad"))
    length(beams) == length(r.beams) ||
        throw(ArgumentError("the report has $(length(r.beams)) beams and $(length(beams)) were given"))
    kind === :auto && (kind = any(u -> u.kind === :quad, r.units) ? :quad : :beam)
    is = [i for (i, u) in enumerate(r.units) if u.kind === kind]
    isempty(is) && throw(ArgumentError("the report has no $kind units"))
    marks = [_port_mark(r.units[i].name, r.units[i].beams, beams) for i in is]
    notes = map(is) do i
        o = r.singles[i]
        c = r.criteria[o.binding]
        string("role ", r.units[i].role, ": ", _ROE_DECISION[r.tiers[i]], ", ", c.name, " ",
               _roe_number(o.scores[o.binding]), " (green ≤ ", _roe_number(c.green), ", amber ≤ ",
               _roe_number(c.amber), ")")
    end
    roles = Symbol[r.units[i].role for i in is]
    show = badges === :auto ? length(unique(roles)) > 1 : Bool(badges)
    if show
        bs = role_badges(roles)
        return _port_map_svg(marks; tiers = r.tiers[is], notes, badges = [bs[x] for x in roles],
                             roles = Dict(v => k for (k, v) in bs), kwargs...)
    end
    return _port_map_svg(marks; tiers = r.tiers[is], notes, kwargs...)
end

"""
    role_badges(roles) -> Dict{Symbol,String}

A short badge for each distinct role of `roles` (symbols, such as the parts that absorb
the beams of a scan, [`roe_scan`](@ref)): the initials of the words of its name, split at
`_`, in capitals and with the digits of a word kept, as `HW` for `:hohlraum_wall` and `F2`
for `:foil_2`, at most four characters. A badge that two roles would share takes the
order of the role among those roles, from 1, as a suffix.
"""
function role_badges(roles)
    rs = unique(Symbol[Symbol(x) for x in roles])
    initials(r) = first(join(uppercase(w[1:1]) * replace(w[2:end], r"[^0-9]" => "") for w in split(String(r), '_')
                             if !isempty(w)), 4)
    base = [initials(r) for r in rs]
    out = Dict{Symbol,String}()
    for (k, r) in enumerate(rs)
        b = base[k]
        if count(==(b), base) > 1
            b *= string(count(==(b), base[1:k]))
        end
        out[r] = isempty(b) ? string(k) : b
    end
    return out
end

function _port_map_svg(marks; tiers = nothing, values = nothing, notes = nothing, badges = nothing, roles = nothing,
                       title::AbstractString = "", unit::AbstractString = "", size::Integer = 360, labels = :auto)
    size >= 100 || throw(ArgumentError("size must be at least 100 pixels"))
    tiers === nothing || values === nothing || throw(ArgumentError("give tiers or values, not both"))
    ts = _port_lookup(tiers, marks, "tiers")
    vs = _port_lookup(values, marks, "values")
    ns = _port_lookup(notes, marks, "notes")
    bs = _port_lookup(badges, marks, "badges")
    show_labels = labels === :auto ? length(marks) <= 80 : Bool(labels)
    R = size / 2
    pad, top, bottom = 48, 62, 44
    W, H = 2 * size + 4 * pad, size + top + bottom
    cx = (pad + R, 3 * pad + 3R)
    cy = top + R
    n(v) = @sprintf("%.2f", v)
    # the mark radius: a third of the closest spacing within a hemisphere, within limits
    dmin = Inf
    for a in eachindex(marks), b in (a + 1):length(marks)
        marks[a].hemisphere === marks[b].hemisphere || continue
        dmin = min(dmin, hypot(marks[a].x - marks[b].x, marks[a].y - marks[b].y))
    end
    rm = clamp(0.36 * dmin, 0.012, 0.06) * R
    finite = vs === nothing ? Float64[] : Float64[v for v in vs if v isa Real && isfinite(v)]
    lo, hi = isempty(finite) ? (0.0, 1.0) : extrema(finite)
    function fill(i)
        if ts !== nothing
            return get(_PORT_TIER_COLORS, ts[i], _PORT_NEUTRAL)
        elseif vs !== nothing
            v = vs[i]
            (v isa Real && isfinite(v)) || return _PORT_NEUTRAL
            return _port_ramp(hi > lo ? (v - lo) / (hi - lo) : 0.5)
        end
        return "#0072B2"
    end
    io = IOBuffer()
    println(io, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"", W, "\" height=\"", H, "\" viewBox=\"0 0 ", W,
            " ", H, "\" font-family=\"sans-serif\" font-size=\"11\">")
    println(io, "<rect width=\"", W, "\" height=\"", H, "\" fill=\"#ffffff\"/>")
    isempty(title) || println(io, "<text x=\"", pad, "\" y=\"16\" font-size=\"13\" font-weight=\"bold\">",
                              _xml_escape(title), "</text>")
    for (h, name) in enumerate(("Upper hemisphere", "Lower hemisphere, seen from above"))
        x0 = cx[h]
        println(io, "<text x=\"", n(x0), "\" y=\"", n(top - 22.0), "\" text-anchor=\"middle\">", name, "</text>")
        for (f, θ) in ((1 / 3, 30), (2 / 3, 60), (1.0, 90))
            println(io, "<circle cx=\"", n(x0), "\" cy=\"", n(cy), "\" r=\"", n(f * R), "\" fill=\"none\" stroke=\"",
                    f == 1.0 ? "#444444" : "#CCCCCC", "\" stroke-width=\"1\"/>")
            f < 1 && println(io, "<text x=\"", n(x0 + 3), "\" y=\"", n(cy - f * R - 2), "\" fill=\"#888888\">θ ",
                             h == 1 ? θ : 180 - θ, "°</text>")
        end
        println(io, "<line x1=\"", n(x0 - R), "\" y1=\"", n(cy), "\" x2=\"", n(x0 + R), "\" y2=\"", n(cy),
                "\" stroke=\"#DDDDDD\" stroke-width=\"1\"/>")
        println(io, "<line x1=\"", n(x0), "\" y1=\"", n(cy - R), "\" x2=\"", n(x0), "\" y2=\"", n(cy + R),
                "\" stroke=\"#DDDDDD\" stroke-width=\"1\"/>")
        for (φ, dx, dy, anchor) in ((0, R + 4, 4, "start"), (90, 0, -R - 3, "middle"), (180, -R - 4, 4, "end"),
                                    (270, 0, R + 12, "middle"))
            println(io, "<text x=\"", n(x0 + dx), "\" y=\"", n(cy + dy), "\" text-anchor=\"", anchor,
                    "\" fill=\"#888888\">φ ", φ, "°</text>")
        end
    end
    for (i, m) in enumerate(marks)
        x = cx[m.hemisphere === :upper ? 1 : 2] + m.x * R
        y = cy - m.y * R
        tip = string(m.name, @sprintf(" (θ %.1f°, φ %.1f°)", m.theta, m.phi))
        vs === nothing || vs[i] === nothing || (tip *= string(": ", vs[i] isa Real ? @sprintf("%.4g", vs[i]) : vs[i],
                                                             isempty(unit) ? "" : " " * unit))
        ns === nothing || ns[i] === nothing || (tip *= string(", ", ns[i]))
        println(io, "<g class=\"port\" data-name=\"", _xml_escape(m.name), "\"><title>", _xml_escape(tip),
                "</title><circle cx=\"", n(x), "\" cy=\"", n(y), "\" r=\"", n(rm), "\" fill=\"", fill(i),
                "\" stroke=\"#222222\" stroke-width=\"0.8\"/>",
                show_labels ? string("<text x=\"", n(x), "\" y=\"", n(y + rm + 9), "\" text-anchor=\"middle\" font-size=\"8\">",
                                     _xml_escape(m.name), "</text>") : "",
                bs === nothing || bs[i] === nothing ? "" :
                string("<text x=\"", n(x + 0.7rm), "\" y=\"", n(y - 0.7rm), "\" font-size=\"8\" font-weight=\"bold\" ",
                       "stroke=\"#FFFFFF\" stroke-width=\"2.5\" paint-order=\"stroke\">", _xml_escape(string(bs[i])), "</text>"),
                "</g>")
    end
    yl = H - 14
    if ts !== nothing
        x = pad
        for (t, lab) in ((:red, "Critical"), (:amber, "Ask RI"), (:green, "Auto-Drop"))
            println(io, "<rect x=\"", x, "\" y=\"", yl - 10, "\" width=\"12\" height=\"12\" fill=\"", _PORT_TIER_COLORS[t],
                    "\" stroke=\"#222222\" stroke-width=\"0.8\"/><text x=\"", x + 17, "\" y=\"", yl, "\">", lab, " (",
                    count(==(t), ts), ")</text>")
            x += 110
        end
    elseif vs !== nothing && !isempty(finite)
        for k in 0:9
            println(io, "<rect x=\"", pad + 16k, "\" y=\"", yl - 10, "\" width=\"16\" height=\"12\" fill=\"",
                    _port_ramp(k / 9), "\"/>")
        end
        println(io, "<text x=\"", pad + 166, "\" y=\"", yl, "\">", @sprintf("%.4g", lo), " to ", @sprintf("%.4g", hi),
                isempty(unit) ? "" : " " * _xml_escape(unit), "</text>")
    end
    if roles !== nothing && !isempty(roles)
        key = join((string(b, " ", r) for (b, r) in sort!(collect(roles); by = first)), " · ")
        println(io, "<text x=\"", W - pad, "\" y=\"", yl, "\" text-anchor=\"end\">roles: ", _xml_escape(key), "</text>")
    end
    print(io, "</svg>\n")
    return String(take!(io))
end

"""
    write_port_map_svg(path, beams; kwargs...) -> path
    write_port_map_svg(path, report::ROEReport, beams; kwargs...) -> path

Write the drawing of the chamber port map ([`port_map_svg`](@ref)) to `path`.
"""
function write_port_map_svg(path::AbstractString, args...; kwargs...)
    write(path, port_map_svg(args...; kwargs...))
    return path
end
