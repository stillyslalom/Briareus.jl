# This file tests the headless review outputs: the chamber port map (port_map.jl), the
# HTML packets of the 1ω review and of the rules of engagement, and the command wrappers
# on a mock-up design. Port positions are checked against the published NIF port angles
# and closed-form positions, and packet content against the reports it was built from.

using Briareus, LinearAlgebra, Test

# Whether the tags of an XML or HTML text balance, attribute values are quoted and the
# text between tags has no bare `<`. Script and style bodies are left out.
function _well_formed(text)
    s = replace(text, r"<script>.*?</script>"s => "<script></script>", r"<style>.*?</style>"s => "<style></style>",
                r"<!DOCTYPE[^>]*>" => "", r"<!--.*?-->"s => "")
    stack = String[]
    for m in eachmatch(r"<(/?)([A-Za-z][A-Za-z0-9]*)((?:\s+[A-Za-z_:][-A-Za-z0-9_:.]*=\"[^\"<]*\")*)\s*(/?)>", s)
        name = String(m.captures[2])
        if m.captures[1] == "/"
            (isempty(stack) || pop!(stack) != name) && return false
        elseif m.captures[4] != "/"
            push!(stack, name)
        end
    end
    isempty(stack) || return false
    # every `<` opens one of the tags matched above
    return count(==('<'), s) == length(collect(eachmatch(r"</?[A-Za-z][A-Za-z0-9]*(?:\s+[A-Za-z_:][-A-Za-z0-9_:.]*=\"[^\"<]*\")*\s*/?>", s)))
end

@testset "chamber port map" begin
    # closed-form positions: the pole at the center, the equator on the rim
    spot = SuperGaussianSpot(0.3)
    bs = [Beam("a"; theta = 0.0, phi = 0.0, spot), Beam("b"; theta = 90.0, phi = 0.0, spot),
          Beam("c"; theta = 135.0, phi = 90.0, spot), Beam("d"; theta = 180.0, phi = 0.0, spot),
          Beam("e"; theta = 45.0, phi = 180.0, spot)]
    m = port_map(bs)
    @test [x.name for x in m] == ["a", "b", "c", "d", "e"] && [x.beams for x in m] == [[i] for i in 1:5]
    @test [x.hemisphere for x in m] == [:upper, :upper, :lower, :lower, :upper]
    xy = [(0.0, 0.0), (1.0, 0.0), (0.0, 0.5), (0.0, 0.0), (-0.5, 0.0)]
    @test all(abs(x.x - p[1]) < 1e-12 && abs(x.y - p[2]) < 1e-12 for (x, p) in zip(m, xy))
    @test_throws ArgumentError port_map(bs; by = :quad)
    @test_throws ArgumentError port_map(bs; by = :port)

    # NIF: the quads sit on the four cones of each hemisphere (NIF User Guide angles),
    # and the 192 beams group into the same 48 quads at the same places
    q = port_map(nif_quads())
    cones = [23.5, 30.0, 30.58, 44.5, 50.0]
    @test length(q) == 48 && count(x -> x.hemisphere === :upper, q) == 24
    @test [count(x -> abs(x.theta - c) < 1e-9, q) for c in cones] == [4, 2, 2, 8, 8]
    @test [count(x -> abs(180 - x.theta - c) < 1e-9, q) for c in cones] == [4, 2, 2, 8, 8]
    @test all(x -> abs(hypot(x.x, x.y) - min(x.theta, 180 - x.theta) / 90) < 1e-12, q)
    # chamber ports 1 and 72 of the User Guide table, at (23.5°, 78.75°) and (156.5°, 281.25°)
    byname = Dict(x.name => x for x in q)
    p1, p72 = byname[Briareus.NIF_QUAD_NAMES[1]], byname[Briareus.NIF_QUAD_NAMES[72]]
    @test abs(p1.x - 23.5 / 90 * cosd(78.75)) < 1e-12 && abs(p1.y - 23.5 / 90 * sind(78.75)) < 1e-12
    @test p72.hemisphere === :lower && abs(p72.x - 23.5 / 90 * cosd(281.25)) < 1e-12 &&
          abs(p72.y - 23.5 / 90 * sind(281.25)) < 1e-12
    nb = nif_beams()
    qb = port_map(nb)
    @test length(qb) == 48 && all(x -> length(x.beams) == 4, qb) && length(port_map(nb; by = :beam)) == 192
    @test all(x -> abs(x.theta - byname[x.name].theta) < 0.05 && abs(x.phi - byname[x.name].phi) < 0.05, qb)
    # OMEGA: 60 ports, 30 per hemisphere, no quads
    om = port_map(omega_beams(cones = :all))
    @test length(om) == 60 && count(x -> x.hemisphere === :lower, om) == 30 && all(x -> length(x.beams) == 1, om)

    # the drawing: one group per mark, tiers by name, values on the ramp, text escaped
    tiers = Dict(x.name => (:green, :amber, :red)[mod1(i, 3)] for (i, x) in enumerate(q))
    delete!(tiers, q[1].name)
    svg = port_map_svg(nif_quads(); tiers, title = "A & B <test>")
    @test _well_formed(svg) && startswith(svg, "<svg xmlns=") && endswith(svg, "</svg>\n")
    @test count("<g class=\"port\"", svg) == 48 && occursin("A &amp; B &lt;test&gt;", svg)
    @test count("fill=\"#009E73\"", svg) == 15 + 1 && count("fill=\"#D55E00\"", svg) == 16 + 1 &&
          count("fill=\"#BBBBBB\"", svg) == 1                      # one legend square per tier
    @test occursin("Auto-Drop (15)", svg) && svg == port_map_svg(nif_quads(); tiers, title = "A & B <test>")
    v = port_map_svg(bs; values = [1.0, 2.0, 3.0, NaN, 5.0], unit = "TW", notes = Dict("a" => "note <1>"))
    @test _well_formed(v) && occursin("1 to 5 TW", v) && occursin("a (θ 0.0°, φ 0.0°): 1 TW, note &lt;1&gt;", v)
    @test count("fill=\"#BBBBBB\"", v) == 1 && occursin("fill=\"#DEEBF7\"", v) && occursin("fill=\"#08306B\"", v)
    @test !occursin("<text x=\"", port_map_svg(nb; by = :beam, labels = false)[end-200:end])
    @test_throws ArgumentError port_map_svg(bs; values = [1.0], unit = "TW")
    @test_throws ArgumentError port_map_svg(bs; values = ones(5), tiers = fill(:green, 5))
    path = write_port_map_svg(joinpath(mktempdir(), "ports.svg"), bs)
    @test read(path, String) == port_map_svg(bs)

    # colored by the tiers of a rules-of-engagement report
    main = foil(:backlighter, (0, 0, 0), (0, 0, 1); width = 6.0, n = 8, two_sided = false)
    four = [Beam("B$i"; theta = 30.0, phi = φ, power = 2.0, spot = SuperGaussianSpot(0.4),
                 quad = i <= 2 ? "QA" : "QB") for (i, φ) in enumerate((0.0, 90.0, 180.0, 270.0))]
    r = roe_scan(Scenario(Target([main]), four), four, [roe_foil_irradiance(:backlighter; green = 0.3, amber = 0.6)])
    @test any(u -> u.kind === :quad, r.units)
    sq = port_map_svg(r, four)
    @test _well_formed(sq) && count("<g class=\"port\"", sq) == 2 && occursin("data-name=\"QA\"", sq)
    @test occursin("role backlighter: Ask RI, irradiance 0.5 (green ≤ 0.3, amber ≤ 0.6)", sq)
    sb = port_map_svg(r, four; kind = :beam)
    @test count("<g class=\"port\"", sb) == 4 && occursin("Auto-Drop (4)", sb)
    @test_throws ArgumentError port_map_svg(r, four[1:3])
end

# The ids that the links of a page point to, and the ids the page has.
_page_links(html) = [String(m.captures[1]) for m in eachmatch(r"href=\"#([^\"]+)\"", html)]
_page_ids(html) = [String(m.captures[1]) for m in eachmatch(r"\sid=\"([^\"]+)\"", html)]
# Whether a page stands alone: one script, no linked file, nothing fetched.
_stands_alone(html) = count("<script>", html) == 1 && !occursin("<link", html) && !occursin("src=\"http", html) &&
                      !occursin("href=\"http", html) && !occursin("@import", html) && !occursin("url(", html)

@testset "1ω review packet" begin
    # One NIF beam and two 0.1 mm squares 0.5 mm inside the near and far edges of its 1ω
    # footprint, as in the closest-approach tests. The beam name needs escaping.
    bl = remake(first(nif_beams()); name = "B<1>&")
    k = direction(bl)
    pq = Briareus._nif_quad_port(bl.quad)
    ô = Briareus._own_side(bl.port, pq)
    Rm = hcat(ô, cross(k, ô), k)
    sq(name, c) = Component(name, transform(rectangle_mesh(0.1, 0.1; nu = 2, nv = 2), Rm, c), gold(); two_sided = true)
    t = Target(sq(:own, 3.9 * ô), sq(:other, -3.9 * ô), sq(:far, -12.0 * ô))
    kw = (; capture = (), open_ports = [pq], rtol = 0)
    rv = omega1_review(t, [bl]; threshold = 1.0, kw...)
    ref = talis_check(t, [bl]; kw...)
    same(a, b) = a.beam == b.beam && a.rule === b.rule && a.component === b.component && a.status === b.status &&
                 isequal(a.margin, b.margin)
    @test length(rv.report.rows) == length(ref.rows) && all(same(a, b) for (a, b) in zip(rv.report.rows, ref.rows))
    o1 = [x for x in ref.rows if x.rule in (:clearance_1w, :management_1w, :diagnostic_1w)]
    want = [x for x in o1 if x.status !== :pass || (x.unit === :mm && x.margin < 1.0)]
    @test !isempty(want) && length(want) < length(o1) && length(rv.items) == length(want)
    @test Set((x.row.rule, x.row.component) for x in rv.items) == Set((x.rule, x.component) for x in want)
    @test issorted([x.row.status === :fail ? 0 : x.row.status === :review ? 1 : 2 for x in rv.items])
    @test all(x -> x.approach.margin == omega1_closest_approach(t, [bl], ref, x.row).margin, rv.items)
    @test all(x -> x.slice.margin == x.approach.margin && x.slice.component === x.row.component, rv.items)
    @test occursin("Omega1Review(", sprint(show, rv))
    # the disposition per component: every checked triangle counted once, by area
    d = rv.disposition
    @test d isa Omega1Disposition && d.settings.tight_threshold == 1.0
    @test [c.component for c in rv.components] == [:own, :other, :far]
    area(c) = sum(c[st] for st in (:fail, :review, :specular_clear, :curved_safe, :tight, :clear))
    @test all(c -> isapprox(area(c), sum(t.tri_area[j] for j in eachindex(d.status) if t.components[t.tri_comp[j]].name === c.component);
                            rtol = 1e-12), rv.components)
    @test all(c -> c.margin == minimum(d.margin[j] for j in eachindex(d.status) if t.components[t.tri_comp[j]].name === c.component),
              rv.components)
    @test rv.components[3].worst === :clear && rv.components[1].worst !== :clear
    @test omega1_review(t, [bl]; disposition = false, kw...).disposition === nothing
    @test_throws ArgumentError omega1_review(t, [bl]; threshold = -1.0, kw...)

    html = omega1_review_html(rv; title = "Review <A&B>")
    @test _well_formed(html) && _stands_alone(html) && startswith(html, "<!DOCTYPE html>")
    @test occursin("<title>Review &lt;A&amp;B&gt;</title>", html) && occursin("B&lt;1&gt;&amp;", html) && !occursin("B<1>", html)
    @test issubset(_page_links(html), _page_ids(html)) && allunique(_page_ids(html))
    @test length(_page_links(html)) == length(rv.items) == count("<figure id=\"approach-", html)
    @test count("<svg ", html) == length(rv.items)
    checks = split(split(html, "<table id=\"checks\"")[2], "</table>")[1]
    @test count("<tr>", checks) == length(ref.rows) + 1
    @test all(x -> occursin(">" * Briareus._packet_mm(x.margin) * "</td>", checks), ref.rows)
    @test all(st -> count("class=\"chip $st\"", checks) == count(x -> x.status === st, ref.rows), (:fail, :review, :pass))
    @test all(n -> occursin(Briareus._xml_escape(n), html), ref.notes) && occursin("out of scope: ", html)
    @test html == omega1_review_html(rv; title = "Review <A&B>") && !occursin("<img", html)
    h2 = omega1_review_html(rv; generated = "made on a <date>", images = Dict(1 => "data:image/png;base64,AAAA"))
    @test _well_formed(h2) && count("<img ", h2) == 1 && occursin("made on a &lt;date&gt;", h2)
    # a review with nothing to show still gives a page
    clear = omega1_review(Target(sq(:far, -12.0 * ô)), [bl]; kw...)
    hc = omega1_review_html(clear)
    @test isempty(clear.items) && _well_formed(hc) && occursin("No 1ω check fails", hc) && isempty(_page_links(hc))

    mktempdir() do dir
        paths = write_omega1_review(joinpath(dir, "packet"), rv; json = true)
        @test basename.(paths) == vcat(["review-1w.html", "talis.csv"], ["slice-00$i.svg" for i in 1:length(rv.items)],
                                       ["review-1w.json"])
        @test all(isfile, paths) && read(paths[1], String) == omega1_review_html(rv)
        @test all(p -> _well_formed(read(p, String)), paths[3:end-1])
        js = Briareus.JSON.parsefile(paths[end])
        @test length(js["rows"]) == length(ref.rows) && length(js["approaches"]) == length(rv.items)
        @test js["rows"][1]["beam"] == "B<1>&" && js["threshold_mm"] == 1.0
        @test js["approaches"][1]["margin_mm"] ≈ rv.items[1].approach.margin &&
              js["approaches"][1]["slice"] == "slice-001.svg"
        @test [c["component"] for c in js["disposition"]] == ["own", "other", "far"]
        @test js["counts"]["fail"] == count(x -> x.status === :fail, ref.rows)
        tj = Briareus.JSON.parsefile(write_talis_json(joinpath(dir, "talis.json"), ref))
        @test length(tj["rows"]) == length(ref.rows) && tj["settings"]["tolerance"] == 0.5 &&
              [r["margin"] for r in tj["rows"] if r["margin"] isa Real] ≈ [x.margin for x in ref.rows if isfinite(x.margin)]
        @test length(write_omega1_review(joinpath(dir, "plain"), rv)) == length(paths) - 1
    end
end

@testset "rules-of-engagement packet" begin
    main = foil(:backlighter, (0, 0, 0), (0, 0, 1); width = 6.0, n = 8, two_sided = false)
    side = foil(:second, (20.0, 0, 0), (0, 0, 1); width = 6.0, n = 4, two_sided = false)
    t = Target([main, side])
    four = [Beam("B$i"; theta = 30.0, phi = φ, power = 2.0, spot = SuperGaussianSpot(0.4),
                 quad = i <= 2 ? "QA" : "QB") for (i, φ) in enumerate((0.0, 90.0, 180.0, 270.0))]
    s = Scenario(t, four)
    r = roe_scan(s, four, [roe_foil_irradiance(:backlighter; green = 0.3, amber = 0.6), roe_foil_centroid(:backlighter)];
                 companions = :opposite)
    rows = roe_table(r)
    html = roe_html(r, four; title = "Sheet", notes = ["chosen <by hand>"])
    @test _well_formed(html) && _stands_alone(html) && html == roe_html(r, four; title = "Sheet", notes = ["chosen <by hand>"])
    @test occursin("chosen &lt;by hand&gt;", html) && occursin("<h2>Role: backlighter</h2>", html)
    @test occursin("id=\"ports-quad\"", html) && occursin("id=\"ports-beam\"", html) && count("<svg ", html) == 2
    @test allunique(_page_ids(html)) && issubset(_page_links(html), _page_ids(html))
    @test all(tier -> count("<span class=\"chip $tier\">", html) == count(x -> x.tier === tier, rows), (:green, :amber, :red))
    @test all(x -> occursin("<td>" * x.unit * "</td><td>" * string(x.kind) * "</td>", html), rows)
    @test all(x -> occursin(">" * Briareus._roe_number(x.score) * "</td>", html), rows)
    @test all(b -> occursin(Briareus._xml_escape(b.text), html), r.budgets)
    plain = roe_html(r)
    @test _well_formed(plain) && !occursin("Chamber port map", plain) && !occursin("<svg", plain)

    mktempdir() do dir
        paths = write_roe_packet(joinpath(dir, "roe"), r, four; json = true, notes = ["a note"])
        @test basename.(paths) == ["roe.html", "roe.csv", "roe.md", "ports.svg", "roe.json"] && all(isfile, paths)
        @test read(paths[4], String) == port_map_svg(r, four)
        js = Briareus.JSON.parsefile(paths[5])
        @test [u["unit"] for u in js["units"]] == [x.unit for x in rows] && js["notes"] == ["a note"]
        @test [u["score"] for u in js["units"]] ≈ [x.score for x in rows]
        @test js["counts"]["green"] == count(==(:green), r.tiers) && length(js["criteria"]) == 2
        @test js["criteria"][2]["nominal"] isa Vector && length(js["budgets"]) == length(r.budgets)
        @test basename.(write_roe_packet(joinpath(dir, "bare"), r)) == ["roe.html", "roe.csv", "roe.md"]
    end

    # the criteria chosen from the design
    names(d) = [c.name for c in d.criteria]
    lit = roe_default_criteria(s, four)                    # nothing observed: the lit parts
    @test names(lit) == [:backlighter_irradiance] && occursin("Part backlighter absorbs 100%", lit.notes[1])
    @test !occursin(" %", join(lit.notes)) && occursin("engineering defaults", lit.notes[end])
    two = [four; Beam("S1"; theta = 20.0, phi = 45.0, power = 1.0, spot = SuperGaussianSpot(0.4), aim = (20.0, 0.0, 0.0))]
    @test names(roe_default_criteria(Scenario(t, two), two)) == [:backlighter_irradiance, :second_irradiance]
    @test names(roe_default_criteria(Scenario(t, two), two; share = 0.2)) == [:backlighter_irradiance]
    fo = roe_default_criteria(s, four; foils = :backlighter)
    @test names(fo) == [:backlighter_irradiance, :backlighter_irradiance_rms, :backlighter_centroid]
    @test roe_scan(s, four, fo.criteria; pairs = :none).singles[1].scores[1] ≈ 0.25 rtol = 1e-9
    miss = [remake(b; aim = (0.0, 50.0, 0.0)) for b in four]
    none = roe_default_criteria(Scenario(t, miss), miss)
    @test isempty(none.criteria) && occursin("no criteria", none.notes[1])
    # an observed open face, and an observed closed sphere
    hr = hohlraum(; length = 3.0, radius = 1.5, leh_radius = 1.0, capsule_radius = 0.0, h = 0.5, top_end = 1.0,
                  bottom_end = :open)
    smp = foil(:sample, (0, 0, -1.5), (0, 0, 1); width = 1.0, n = 3, two_sided = false,
               material = Material(:sample; albedo = 0.3, conversion = 0.0))
    th = Target([hr.components; smp]; closed = false)
    up = [remake(b; aim = (0.0, 0.0, 1.5)) for b in nif_quads(power = 4.0) if b.port[3] > 0]
    face = roe_default_criteria(Scenario(th, up; observe = :sample), up)
    @test names(face) == [:face_flux, :face_rms, :Tr] && occursin("open face", face.notes[1])
    hc = hohlraum(; length = 9.0, radius = 2.8, leh_radius = 1.5, capsule_radius = 0.9, h = 1.0, capsule_h = 0.45)
    qs = [remake(b; aim = (0.0, 0.0, sign(b.port[3]) * 3.0)) for b in nif_quads()]
    cap = roe_default_criteria(Scenario(hc, qs; observe = :capsule, viewfactor_kw = (; closed = true)), qs; foils = (:wall,))
    @test names(cap) == [:drive, :P2, :P4, :hotspot_velocity, :wall_irradiance, :wall_irradiance_rms, :wall_centroid]
    @test occursin("close on themselves", cap.notes[1])
end
