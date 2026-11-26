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
    # role badges: none for one role unless asked, one per mark for two roles
    @test !occursin("roles: ", sq) && occursin("roles: B backlighter", port_map_svg(r, four; badges = true))
    rb = role_badges([:hohlraum_wall, :foil_1, :foil_2, :hohlraum_wall, :halfraum_a, :half_raum_a])
    @test rb == Dict(:hohlraum_wall => "HW", :foil_1 => "F1", :foil_2 => "F2", :halfraum_a => "HA",
                     :half_raum_a => "HRA")
    @test role_badges([:capsule, :cbi_crystal_2]) == Dict(:capsule => "C", :cbi_crystal_2 => "CC2")
    second = foil(:second, (0, 0, -40), (0, 0, -1); width = 6.0, n = 8, two_sided = false)
    mixed = [four[1:2]; [Beam("B$i"; theta = 150.0, phi = φ, power = 2.0, spot = SuperGaussianSpot(0.4), quad = "QB",
                              aim = (0, 0, -40)) for (i, φ) in ((3, 180.0), (4, 270.0))]]
    r2 = roe_scan(Scenario(Target([main, second]), mixed), mixed, [roe_foil_irradiance(:backlighter; green = 0.3, amber = 0.6)])
    @test [u.role for u in r2.units if u.kind === :quad] == [:backlighter, :second]
    s2 = port_map_svg(r2, mixed)
    @test _well_formed(s2) && count("font-weight=\"bold\" stroke=\"#FFFFFF\"", s2) == 2 &&
          occursin("roles: B backlighter · S second", s2)
    @test !occursin("roles: ", port_map_svg(r2, mixed; badges = false))
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
    # the 1ω checks by part: one row per part, with the counts and the least margin
    o1 = [x for x in ref.rows if x.rule in Briareus._OMEGA1_RULES]
    bypart = split(split(html, "<table id=\"parts\"")[2], "</table>")[1]
    @test count("<tr>", bypart) == length(unique(x.component for x in o1)) + 1
    @test all(c -> occursin(">" * Briareus._packet_mm(minimum(x.margin for x in o1 if x.component === c && x.unit === :mm)) *
                            "</td>", bypart), unique(x.component for x in o1 if x.unit === :mm))
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
    # the summary by role: one row per role and kind, with the counts of each tier
    roles = split(split(html, "<table id=\"roles\"")[2], "</table>")[1]
    @test count("<tr>", roles) == length(unique((u.role, u.kind) for u in r.units)) + 1
    @test occursin("<td>backlighter</td><td>beam</td>", roles)
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

@testset "design diff" begin
    # Hand-built pairs with known changes.
    plate(name, c; w = 1.0, mat = gold()) = Component(name, translate(rectangle_mesh(w, w; nu = 2, nv = 2), c), mat)
    ball(R) = Component(:ball, icosphere_mesh(R; n = 2), Material(:CH); surface = SphereSurface(R))
    shift = Briareus.Vec3(0.3, 0.0, -0.2)
    ta = Target([plate(:plate, (0.0, 0.0, 2.0)), ball(1.0), plate(:same, (5.0, 0.0, 0.0)), plate(:old, (0.0, 5.0, 0.0)),
                 plate(:coated, (0.0, -5.0, 0.0))])
    tb = Target([plate(:plate, Briareus.Vec3(0.0, 0.0, 2.0) + shift), ball(1.25), plate(:same, (5.0, 0.0, 0.0)),
                 plate(:new, (0.0, 6.0, 0.0)), plate(:coated, (0.0, -5.0, 0.0); mat = Material(:CH))])
    spot = SuperGaussianSpot(0.3)
    ba = [Beam("a"; theta = 30.0, phi = 0.0, spot, power = 1.0), Beam("b"; theta = 150.0, phi = 90.0, spot, power = 2.0),
          Beam("c"; theta = 90.0, phi = 10.0, spot), Beam("gone"; theta = 60.0, phi = 0.0, spot)]
    bb = [remake(ba[1]; aim = (0.1, -0.2, 0.2)), remake(ba[2]; power = 2.5, defocus = 0.4), ba[3],
          Beam("extra"; theta = 70.0, phi = 0.0, spot)]
    d = design_diff(ta, ba, tb, bb; names = ("old", "new"))
    part(n) = only(x for x in d.parts if x.part === n)
    @test [x.part for x in d.parts] == [:plate, :ball, :same, :old, :coated, :new]
    @test part(:plate).status === :moved && part(:plate).shift ≈ shift && part(:plate).distance ≈ norm(shift)
    @test part(:plate).displacement ≈ norm(shift)
    # a sphere that grows about its center: every vertex moves by 0.25 and the centroid stays
    @test part(:ball).status === :changed && part(:ball).displacement ≈ 0.25 && part(:ball).distance < 1e-12
    @test part(:same).status === :same && part(:same).distance == 0 && !part(:same).material
    @test part(:old).status === :removed && part(:new).status === :added
    @test part(:coated).status === :same && part(:coated).material
    beam(n) = only(x for x in d.beams if x.beam == n)
    @test beam("a").changes == [:aim] && beam("a").aim_shift ≈ [0.1, -0.2, 0.2] && beam("a").distance ≈ 0.3
    @test beam("b").changes == [:power, :defocus] && beam("b").power == (2.0, 2.5) && beam("b").defocus == (0.0, 0.4)
    @test beam("c").status === :same && beam("gone").status === :removed && beam("extra").status === :added
    @test isempty(d.talis) && isempty(d.physics)

    # Reports built by hand: one row of each kind of change.
    row(beam, comp, status, margin; rule = :clearance_1w) =
        Briareus.TaLISRow(beam, rule, comp, status, margin, :mm, margin, margin, 0.2, zero(Briareus.Vec3), "")
    rep(rows) = Briareus.TaLISReport(rows, (;), String[], Briareus.Omega1Reflection[])
    ra = rep([row("a", :p, :pass, 1.0), row("b", :p, :fail, -0.5), row("c", :p, :pass, 2.0), row("d", :p, :pass, 2.0),
              row("e", :p, :pass, 3.0), row("f", :p, :pass, 1.0), row("g", :p, :fail, -1.0), row("h", :p, :pass, Inf)])
    rb = rep([row("a", :p, :fail, -0.25), row("b", :p, :pass, 0.5), row("c", :p, :pass, 1.5), row("d", :p, :pass, 2.75),
              row("e", :p, :pass, 3.0), row("g", :p, :fail, -2.0), row("h", :p, :pass, Inf), row("i", :p, :fail, -0.1),
              row("j", :p, :review, 0.1)])
    dt = design_diff(ta, ba, ta, ba; reports = (ra, rb))
    ch = Dict(x.beam => x for x in dt.talis)
    @test [ch[k].change for k in ("a", "b", "c", "d", "e", "f", "g", "h", "i", "j")] ==
          [:newly_failing, :fixed, :worse, :improved, :same, :removed, :worse, :same, :newly_failing, :added]
    @test ch["a"].delta ≈ -1.25 && ch["d"].delta ≈ 0.75 && ch["h"].delta == 0 && isnan(ch["f"].delta)
    @test ch["a"].status == (:pass, :fail) && ch["i"].status == (:absent, :fail)
    @test all(x -> x.status === :same, dt.parts) && all(x -> x.status === :same, dt.beams)

    # The check run on both designs. A plate that moves 0.15 mm toward a 1ω square loses
    # margin, and at most 0.15 mm of it since a distance is 1-Lipschitz in the position.
    # Moved onto the axis of the quad, the plate newly fails the 3ω clearance.
    nb = filter(b -> b.quad == nif_beams()[1].quad, nif_beams())
    bl = nb[1]
    k = direction(bl)
    ô = Briareus._own_side(bl.port, Briareus._nif_quad_port(bl.quad))
    Rm = hcat(ô, cross(k, ô), k)
    c0 = Briareus.Vec3(bl.aim) - 12k
    shield(off) = Component(:shield, transform(rectangle_mesh(1.0, 1.0; nu = 3, nv = 3), Rm, c0 + off * ô), gold(); two_sided = true)
    wall = Component(:wall, transform(rectangle_mesh(6.0, 6.0; nu = 4, nv = 4), Rm, Briareus.Vec3(bl.aim)), gold(); two_sided = true)
    kw = (; capture = (:wall,), hardware = (:shield,))
    t0, t1, t2 = Target([wall, shield(6.0)]), Target([wall, shield(5.85)]), Target([wall, shield(0.0)])
    d1 = design_diff(t0, nb, t1, nb; talis = kw)
    r1 = only(x for x in d1.talis if x.rule === :clearance_1w && x.beam == nb[3].name && x.component === :shield)
    @test r1.status == (:pass, :pass) && r1.margin[1] > 5 && -0.15 - 1e-9 <= r1.delta < -0.05 && r1.change === :worse
    @test only(x for x in d1.parts if x.part === :shield).status === :moved
    d2 = design_diff(t0, nb, t2, nb; talis = kw)
    r2 = only(x for x in d2.talis if x.rule === :clearance_3w && x.beam == bl.name && x.component === :shield)
    @test r2.status == (:pass, :fail) && r2.change === :newly_failing

    # Physics: every power raised by 10% in a hohlraum with a capsule. The radiation is
    # linear in the power at a constant albedo and conversion, so the powers and the flux
    # rise by 10%, the temperature by 1.1^(1/4), and the mode ratios stay.
    cav = hohlraum(; length = 10.1, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9,
                   wall = gold(albedo = 0.85, conversion = 0.8), h = 0.9, capsule_h = 0.3)
    hb = [remake(b; aim = (0.2, 0.0, sign(b.port[3]) * (b.group in (:inner23, :inner30) ? 3.3 : 5.1)))
          for b in nif_quads(power = 8.0)]
    s = Scenario(cav, hb; observe = :capsule, viewfactor_kw = (; closed = true))
    hb2 = [remake(b; power = 1.1 * b.power) for b in hb]
    dp = design_diff(cav, hb, cav, hb2; evaluations = (evaluate(s, hb), evaluate(s, hb2)), foils = (:wall,))
    q(role, name) = only(x for x in dp.physics if startswith(x.role, role) && x.quantity == name)
    @test q("laser", "absorbed power").b ≈ 1.1 * q("laser", "absorbed power").a && q("laser", "absorbed power").a > 0
    @test q("part wall", "absorbed power").delta ≈ 0.1 * q("part wall", "absorbed power").a
    @test q("capsule", "mean flux").b ≈ 1.1 * q("capsule", "mean flux").a rtol = 1e-10
    @test q("capsule", "Tr").b ≈ 1.1^0.25 * q("capsule", "Tr").a rtol = 1e-10
    @test abs(q("capsule", "P2/P0").delta) < 1e-9 && abs(q("capsule", "mode 1").delta) < 1e-9
    @test q("capsule", "hot-spot velocity x").a ≈ hotspot_velocity(evaluate(s, hb))[1] && abs(q("capsule", "hot-spot velocity x").a) > 1
    @test q("foil wall", "absorbed intensity").b ≈ 1.1 * q("foil wall", "absorbed intensity").a
    @test abs(q("foil wall", "intensity RMS over mean").delta) < 1e-9
    @test all(x -> x.changes == [:power], dp.beams)

    # The page and the JSON.
    full = Briareus.DesignDiff(("old", "new <b>"), d.parts, d.beams, dt.talis, dp.physics)
    page = design_diff_html(full; title = "A & B")
    @test _well_formed(page) && occursin("A &amp; B", page) && occursin("new &lt;b&gt;", page)
    @test occursin("newly failing", page) && occursin("id=\"physics\"", page) && !occursin("<td>e</td>", page)
    @test occursin("<td>e</td>", design_diff_html(full; unchanged = true))
    @test page == design_diff_html(full; title = "A & B")
    mktempdir() do dir
        ps = write_design_diff(dir, full; json = true)
        @test basename.(ps) == ["diff.html", "diff.json"] && all(isfile, ps)
        js = Briareus.JSON.parsefile(ps[2])
        @test length(js["talis"]) == length(dt.talis) && js["parts"][1]["status"] == "moved"
        @test js["names"] == ["old", "new <b>"]
    end
    @test occursin("2 newly failing", sprint(show, dt))
end
