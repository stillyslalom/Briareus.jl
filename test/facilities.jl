# Facility beam sets against their published counts and symmetries: NIF's 48 quads and
# 192 beams (NIF User Guide Table 4-1, the LMF beam table and chamber-geometry figures),
# OMEGA's 60 beams on a truncated icosahedron, the CPP spots' half widths, and the rows
# that build them.

@testset "facilities" begin
    q = nif_quads()
    @test length(q) == 48
    @test length(nif_beams()) == 192
    angs = port_angles.(q)
    tops = filter(a -> a[1] < 90, angs)
    bots = filter(a -> a[1] > 90, angs)
    # every bottom port is (180-θ, 360-φ) of a top port
    @test all(any(isapprox(180 - b[1], a[1]; atol = 1e-6) && isapprox(mod(360 - b[2], 360), a[2]; atol = 1e-6)
                  for a in tops) for b in bots)
    # four quads at 23.5° and 30°, eight at 44.5° and 50°, in each hemisphere
    @test [count(b -> b.group === g, q) for g in (:inner23, :inner30, :outer44, :outer50)] == [8, 8, 16, 16]
    @test length(tops) == length(bots) == 24
    ob = omega_beams(; cones = :all)
    @test length(ob) == 60
    θs = sort(unique(round.(first.(port_angles.(ob)); digits = 2)))
    @test θs ≈ [21.42, 42.02, 58.85, 81.25, 98.75, 121.15, 137.98, 158.58]
    @test length(omega_beams()) == 40
    # quad names from the LMF chamber-geometry figures: a bijection onto the ports,
    # every bundle one top and one bottom quad, FABS on Q31B (30°) and Q36B (130°, 186°)
    @test length(unique(b.name for b in q)) == 48
    @test all(nif_port(Briareus.NIF_QUAD_NAMES[p]) == p for (p, _, _) in Briareus.NIF_PORTS)
    @test_throws KeyError nif_port("Q99T")
    bundles = nif_bundles(q)
    @test length(bundles) == 24
    @test all(sort([q[i].name[end] for i in v]) == ['B', 'T'] for v in values(bundles))
    @test port_angles(q[findfirst(b -> b.name == "Q36B", q)])[1] ≈ 130
    @test port_angles(q[findfirst(b -> b.name == "Q36B", q)])[2] ≈ 185.62
    @test q[findfirst(b -> b.name == "Q31B", q)].group === :inner30
    @test nif_port("Q11T") == 17 && nif_port("Q11B") == 50
    kinds = Dict{Tuple{Bool,Bool},Int}()
    for v in values(bundles)
        inner = sort([q[i].group in (:inner23, :inner30) for i in v])
        kinds[(inner[1], inner[2])] = get(kinds, (inner[1], inner[2]), 0) + 1
    end
    @test kinds == Dict((false, false) => 12, (false, true) => 8, (true, true) => 4)
    @test all(nif_bundle(b) == nif_bundle(b.quad) for b in nif_beams())
    # the LMF per-beam table: each quad's mean beam direction is its port, which also
    # confirms the quad names read from the chamber-geometry figures
    nb = nif_beams()
    @test length(nb) == 192 && nb[1].name == "B111" && nb[1].quad == "Q11T"
    @test all(acosd(clamp(dot(normalize(sum(x.port for x in nb if x.quad == b.quad)), b.port), -1, 1)) < 0.01
              for b in q)
    @test all(x.group === q[findfirst(y -> y.quad == x.quad, q)].group for x in nb)
    @test all(count(x -> x.quad == b.quad, nb) == 4 for b in q)
end

@testset "facility spots" begin
    # the CPP spots exp(-ρⁿ) fall to half their peak at the tabulated half widths (major
    # along φ̂, minor along θ̂), and `scale` scales both axes
    hwhm = Dict(:inner23 => (0.88585, 0.65579), :inner30 => (0.85838, 0.60591),
                :outer44 => (0.63727, 0.37234), :outer50 => (0.61910, 0.36187))
    half(s, r, a) = exp(-(r / a)^s.order)
    @test all(begin
                  s = nif_spot(c)
                  half(s, hwhm[c][1], s.a) ≈ 0.5 && half(s, hwhm[c][2], s.b) ≈ 0.5 && s.angle == π / 2
              end for c in keys(hwhm))
    s2 = nif_spot(:outer44; scale = 1.5)
    @test s2.a ≈ 1.5 * nif_spot(:outer44).a && s2.b ≈ 1.5 * nif_spot(:outer44).b && s2.order == nif_spot(:outer44).order
    # OMEGA's EIDI-300 plate: 1/e diameters 290 × 212 µm
    e = omega_spot()
    @test (2e.a, 2e.b, e.order) == (0.29, 0.212, 4.1)
    @test omega_spot(:SG4).a == 0.36 && omega_spot(:SG8).order == 4.7
    @test_throws ArgumentError omega_spot(:SG5)
end

@testset "OMEGA ports: a truncated icosahedron about a pentagon axis" begin
    ps = [Vec3(sind(r.theta) * cosd(r.phi), sind(r.theta) * sind(r.phi), cosd(r.theta)) for r in omega_ports()]
    near(p) = minimum(norm(p - x) for x in ps)
    # centrosymmetric, and fivefold about z (the pentagon axis): rotations by 72° and the
    # inversion map the ports onto themselves
    Rz = [cosd(72) -sind(72) 0; sind(72) cosd(72) 0; 0 0 1]
    @test all(near(-p) < 1e-9 && near(Vec3(Rz * p)) < 1e-9 for p in ps)
    # every vertex has three neighbors (one hexagon–hexagon edge, two pentagon edges),
    # the pentagon edges 1.2 times the hexagon–hexagon edge, up to the projection on the sphere
    d = [sort([norm(p - x) for x in ps])[2:5] for p in ps]
    @test all(x -> x[2] ≈ x[3] && x[4] > 1.1x[3], d)
    @test all(x -> 1.15 < x[2] / x[1] < 1.25, d)
    # five beams in each of the two polar cones, ten in each of the other two, per hemisphere
    g = [r.group for r in omega_ports()]
    @test [count(==(Symbol(:cone, k)), g) for k in 1:4] == [10, 10, 20, 20]
    @test count(p -> p[3] > 0, ps) == 30
end

@testset "beams_from_table" begin
    rows = [(name = "a", theta = 30.0, phi = 45.0),
            Dict(:name => :b, :theta => 120.0, :phi => 10.0, :group => :g2, :quad => "Q1",
                 :power => 3.0, :aim => (0.1, 0.2, 0.3), :defocus => 0.5, :fnumber => 7.0,
                 :spot => SuperGaussianSpot(0.2))]
    # defaults: a value, or a function of the row
    bs = beams_from_table(rows; spot = r -> SuperGaussianSpot(r[:theta] / 100), fnumber = 9.0, power = 2.0)
    a, b = bs
    @test a.name == "a" && a.group === :none && a.quad == "" && a.power == 2.0 && a.fnumber == 9.0 &&
          a.spot.a ≈ 0.3 && a.aim == Vec3(0, 0, 0) && a.defocus == 0
    @test all(port_angles(a) .≈ (30.0, 45.0)) && all(port_angles(b) .≈ (120.0, 10.0))
    # the beam direction is the port's opposite: k = -port
    @test direction(a) ≈ -Vec3(sind(30) * cosd(45), sind(30) * sind(45), cosd(30))
    @test b.name == "b" && b.group === :g2 && b.quad == "Q1" && b.power == 3.0 && b.fnumber == 7.0 &&
          b.spot.a == 0.2 && b.aim ≈ Vec3(0.1, 0.2, 0.3) && b.defocus == 0.5
end
