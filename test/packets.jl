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
