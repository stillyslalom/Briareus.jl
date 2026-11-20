# VisRad interchange (src/visrad.jl): round trips and rows copied from the VisRad 21.1.1
# manual, https://prism-cs.com/Manuals/VisRad/.

using Briareus: _vr_key, _vr_rotation, _vr_reference_angle

@testset "VisRad interchange" begin

vrdir = mktempdir()
vrfile(name, text) = (p = joinpath(vrdir, name); write(p, text); p)
same_spot(s, t; rtol = 1e-9) = isapprox(s.a, t.a; rtol) && isapprox(s.b, t.b; rtol) &&
    isapprox(s.order, t.order; rtol) && (s.a == s.b || isapprox(mod(s.angle - t.angle + π / 2, π) - π / 2, 0; atol = 1e-9))

@testset "names" begin
    @test _vr_key("BL115") == "B115"          # NIF "BL" prefix (laser_beam_data_files.html)
    @test _vr_key("b115") == "B115"
    @test _vr_key("Beam_10") == "Beam_10"     # OMEGA
    @test _vr_key("10") == "Beam_10"          # older "Beam #" column
    @test _vr_key("BL-25") == "Beam_25"
    @test _vr_key("Beam_25_P9") == "Beam_25_P9"
    @test _vr_key("Q31B") == "Q31B"
    @test _vr_key("B01") == "B01"
end

@testset "NIF SST pointing" begin
    # header exactly as in the SST import screenshot (revision_history.html, v21.0.0)
    @test VISRAD_SST_COLUMNS[2] == "Pointing: Beam Aim Point X (mm)"
    @test VISRAD_SST_COLUMNS[5] == "Pointing: Beam Aim Point Beam Defocus (mm)"
    # per-beam round trip, some beams off
    nb = nif_beams()
    bs = [i % 5 == 0 ? remake(b; power = 0.0) :
          remake(b; aim = (0.013i, -0.007i, 0.5 - 0.002i), defocus = 0.003i - 0.2) for (i, b) in enumerate(nb)]
    p = write_visrad_sst(joinpath(vrdir, "beams_sst.csv"), bs)
    @test count(!isempty, readlines(p)) == 1 + count(b -> b.power > 0, bs)
    rb = read_visrad_sst(p; template = nif_beams())
    @test [b.name for b in rb] == [b.name for b in bs]
    @test all(isapprox(a.aim, b.aim; atol = 1e-9) for (a, b) in zip(rb, bs))
    @test all(isapprox(a.defocus, b.defocus; atol = 1e-9) for (a, b) in zip(rb, bs) if b.power > 0)
    @test [a.power for a in rb] == [b.power for b in bs]
    @test all(a.port == b.port && a.spot == b.spot for (a, b) in zip(rb, bs))
    # quads: expanded to four beams on write, aggregated on read
    q = nif_quads()
    qs = [i % 7 == 0 ? remake(b; power = 0.0) : remake(b; aim = (0.1i, -0.05i, 0.3), defocus = 0.01i)
          for (i, b) in enumerate(q)]
    pq = write_visrad_sst(joinpath(vrdir, "quads_sst.csv"), qs)
    rq = read_visrad_sst(pq; template = nif_quads())
    @test all(isapprox(a.aim, b.aim; atol = 1e-9) for (a, b) in zip(rq, qs))
    @test [a.power for a in rq] == [b.power for b in qs]
    # the per-beam file of the quads reads into per-beam beams pointed like their quads
    rbq = read_visrad_sst(pq; template = nif_beams())
    i31 = findfirst(b -> b.name == "Q31B", qs)
    @test all(b.aim ≈ qs[i31].aim for b in rbq if b.quad == "Q31B")
    # a foreign file: preamble, BL names, cm units, extra column, a beam without pointing
    pf = vrfile("foreign_sst.csv", """
        Shot N000000-000-999 pointing export
        Beam,Pointing: Beam Aim Point X (cm),Pointing: Beam Aim Point Y (cm),Pointing: Beam Aim Point Z (cm),Pointing: Beam Aim Point Beam Defocus (cm),Comment
        BL115,0.1,0.2,-0.3,0.05,first
        B116,,,,,no pointing
        """)
    rf = @test_logs min_level = Base.CoreLogging.Warn read_visrad_sst(pf; template = nif_beams())
    i115 = findfirst(b -> b.name == "B115", rf)
    @test rf[i115].aim ≈ [1.0, 2.0, -3.0] && rf[i115].defocus ≈ 0.5 && rf[i115].power == 2.0
    @test rf[i115+1].power == 0                      # B116: no pointing data → off
    @test count(b -> b.power > 0, rf) == 1
    @test_throws ArgumentError read_visrad_sst(vrfile("bad.csv", "Beam,X\nB115,1\n"))
end

@testset "Laser Beam Summary: manual example" begin
    # Rows copied from the v21 summary dialog in laser_summary_table.html
    # (images/laser_summary.jpg, OMEGA). Power and later columns are cut off there.
    pm = vrfile("manual_summary.csv", """
        Cluster,Beam ID,On/Off,Polar Ang,Azm Ang,Pt Ref-CS,Pt Geom,Pt Param-1,Pt Param-2,Pt Param-3,Focus offset,TC Pt Pos-x,TC Pt Pos-y,TC Pt Pos-z
        Backlighter_beams,Beam_10,On,58.85,5.94,Pb_Bi_foil,Cartesian,0,0,0,-0.215,-0.03643,0.09028,-0.004198
        Cone_3,Beam_25,On,158.59,90.00,LEH_top,Cartesian,0,0,0,0,-0.05104,0.01658,-0.02684
        Backlighter_beams,Beam_28,On,58.85,318.06,Pb_Bi_foil,Cartesian,0,0,0,-0.215,-0.08254,-0.05163,-0.004198
        Backlighter_beams,Beam_31,On,42.02,342.00,Pb_Bi_foil,Cartesian,0,0,0,-0.215,-0.08254,0.02682,0.04428
        Cone_2,Beam_40,On,81.25,185.46,LEH_top,Cartesian,0,0,0.01,0,-0.05954,0.01935,-0.03131
        Cone_3,Beam_50,On,137.98,90.00,LEH_top,Cartesian,0,0,0,0,-0.05104,0.01658,-0.02684
        """)
    mb = read_visrad_summary(pm)
    byname = Dict(b.name => b for b in mb)
    @test length(mb) == 6 && all(b -> b.power == 1.0, mb)
    @test all(isapprox.(port_angles(byname["Beam_31"]), (42.02, 342.0); atol = 1e-9))
    @test byname["Beam_10"].defocus ≈ -2.15                   # cm → mm
    # TC Pt Pos is the focus: three beams aimed at the foil origin with focus offset
    # -0.215 cm decode to one aim point (to the table's 4-digit precision)
    foil_aims = [byname[n].aim for n in ("Beam_10", "Beam_28", "Beam_31")]
    @test maximum(norm(a - foil_aims[1]) for a in foil_aims) < 2e-3
    # zero focus offset: TC Pt Pos is the reference point itself
    @test byname["Beam_25"].aim ≈ 10 .* [-0.05104, 0.01658, -0.02684]
    @test byname["Beam_50"].aim ≈ byname["Beam_25"].aim
    # Beam_40 is pointed 0.01 cm along the LEH axis from the LEH origin
    @test norm(byname["Beam_40"].aim - byname["Beam_25"].aim) ≈ 0.1 atol = 2e-3
    # the older dialog (intro/images/laser_summary.jpg): "Beam #", "Pt CS", "Ro"
    po = vrfile("old_summary.csv", """
        Cluster,Beam #,On/Off,Polar Ang,Azm Ang,Pt Ref CS,Pt CS,Pt X-1,Pt X-2,Pt X-3,Focus offset,Power,Sp Prof,Ro
        Cone_2b,16,On,98.75,318.54,end_cap_2,Cartesian,0,0,0.025,0,0.5,Uniform,50
        """)
    @test_throws ArgumentError read_visrad_summary(po)   # pointed in a component frame, no TC columns
    pt = vrfile("old_summary_tc.csv", """
        Cluster,Beam #,On/Off,Polar Ang,Azm Ang,Pt Ref CS,Pt CS,Pt X-1,Pt X-2,Pt X-3,Focus offset,Power,Sp Prof,Ro
        Cone_2b,16,On,98.75,318.54,Target Chamber,Cartesian,0,0,0.025,0.01,0.5,Uniform,50
        Cone_1b,10,Off,58.58,5.94,TCC,Cartesian,0,0,0,0.06,0.5,Uniform,50
        Cone_1b,21,On,121.15,329.94,Target Chamber,Spherical,0.1,90,90,0,T-0.45,Uniform,50
        """)
    ot = @test_logs (:warn, r"time-dependent") match_mode = :any read_visrad_summary(pt)
    @test [b.name for b in ot] == ["Beam_16", "Beam_10", "Beam_21"]
    @test ot[1].aim ≈ [0, 0, 0.25] && ot[1].defocus ≈ 0.1 && ot[1].power == 0.5
    @test ot[2].power == 0 && ot[2].aim == zeros(3)
    @test isapprox(ot[3].aim, [0, 1.0, 0]; atol = 1e-12) && ot[3].power ≈ 0.45
    @test ot[1].spot.a ≈ 0.05 && ot[1].spot.a == ot[1].spot.b   # Ro = 50 µm
    @test ot[1].spot.order == 50                                # "Uniform" → flat top
end

@testset "Laser Beam Summary round trip" begin
    rng = Random.MersenneTwister(7)
    nb = nif_beams()
    bs = map(enumerate(nb)) do (i, b)
        sp = SuperGaussianSpot(0.5 + 0.4rand(rng), 0.3 + 0.2rand(rng); order = 2 + 6rand(rng),
                               angle = π * rand(rng))
        remake(b; aim = 0.5 .* randn(rng, 3), defocus = randn(rng), power = i % 9 == 0 ? 0.0 : 1 + rand(rng),
               spot = sp)
    end
    p = write_visrad_summary(joinpath(vrdir, "summary.csv"), bs)
    rb = read_visrad_summary(p; template = nif_beams())
    @test all(isapprox(a.aim, b.aim; atol = 1e-8) for (a, b) in zip(rb, bs))
    @test all(isapprox(a.defocus, b.defocus; atol = 1e-8) for (a, b) in zip(rb, bs))
    @test all(isapprox(a.power, b.power; atol = 1e-9) for (a, b) in zip(rb, bs))
    @test all(same_spot(a.spot, b.spot) for (a, b) in zip(rb, bs))
    # no template: beams rebuilt from the file's port angles
    rn = read_visrad_summary(p)
    @test [b.name for b in rn] == [b.name for b in bs]
    @test all(isapprox(a.port, b.port; atol = 1e-9) for (a, b) in zip(rn, bs))
    @test all(same_spot(a.spot, b.spot) for (a, b) in zip(rn, bs))
    # the other rotation sense is self-consistent but mirrors the spot
    pt = write_visrad_summary(joinpath(vrdir, "summary_t.csv"), bs; rotation_sense = :target)
    rt = read_visrad_summary(pt; template = nif_beams(), rotation_sense = :target)
    @test all(same_spot(a.spot, b.spot) for (a, b) in zip(rt, bs))
    # a spot with a < b is written as major/minor and read back as the same profile
    b0 = remake(nb[1]; spot = SuperGaussianSpot(0.3, 0.6; order = 4.0, angle = 0.4))
    r0 = only(read_visrad_summary(write_visrad_summary(joinpath(vrdir, "ab.csv"), [b0])))
    x = b0.aim + 0.3 * normalize(cross(b0.port, [0.2, 0.5, 0.7]))
    @test intensity(r0, x) ≈ intensity(b0, x) rtol = 1e-8
    # NIF spots (major axis along φ̂) have rotation 0 in either sense
    @test _vr_rotation(nif_spot(:outer44), :beam) == 0 == _vr_rotation(nif_spot(:outer44), :target)
    # VisRad's rotation is relative to the quad's horizontal: at rotation 0 a NIF beam's
    # major axis lies along φ̂ of its quad's port, which for B111 is turned -1.60° from
    # the beam φ̂ (VisRad 21 deposition on a plate normal to B111 fitted -1.58°); quads and
    # beams of no quad keep their φ̂
    b111 = nb[findfirst(b -> b.name == "B111", nb)]
    ref = _vr_reference_angle(b111)
    @test rad2deg(ref) ≈ -1.60 atol = 0.01
    @test _vr_reference_angle(nif_quads()[1]) == 0 && _vr_reference_angle(omega_beams()[1]) == 0
    @test _vr_rotation(SuperGaussianSpot(0.6, 0.4; angle = π / 2 + ref), :beam, ref) ≈ 0 atol = 1e-12
    # quads expand and aggregate, OMEGA beams keep their names
    q = nif_quads()
    q2 = [remake(b; aim = (0.02i, 0.0, -0.01i), defocus = 0.1) for (i, b) in enumerate(q)]
    rq = read_visrad_summary(write_visrad_summary(joinpath(vrdir, "q.csv"), q2); template = nif_quads())
    @test all(isapprox(a.aim, b.aim; atol = 1e-8) && isapprox(a.power, b.power; rtol = 1e-12)
              for (a, b) in zip(rq, q2))
    ob = omega_beams()
    ro = read_visrad_summary(write_visrad_summary(joinpath(vrdir, "o.csv"), ob); template = omega_beams())
    @test all(same_spot(a.spot, b.spot) && a.power ≈ b.power for (a, b) in zip(ro, ob))
end

@testset "phase plates" begin
    # VisRad's NIF CPP table (HWHM radii) against Briareus's measured-profile fits
    for (plate, cone) in (("Rev3 Inner Cone 2W CPP", :inner23), ("1.0 Rev1a Inner Cone 2W CPP", :inner30),
                          ("Rev3 Outer Cone 2W CPP", :outer44), ("1.0 Rev1 Outer Cone 2W CPP", :outer50))
        s, t = visrad_spot(plate), nif_spot(cone)
        @test isapprox(s.a, t.a; rtol = 0.02) && isapprox(s.b, t.b; rtol = 0.02)
        @test s.angle ≈ t.angle
    end
    s = visrad_spot("OMEGA E-IDI-300 DPP")
    @test s.a ≈ 0.144 && s.b ≈ 0.106 && s.order == 4.3      # 1/e radii in µm in the manual
    @test visrad_spot("Type F").a ≈ 0.315 / log(1 / 0.03)^(1 / 2)
    @test_throws ArgumentError visrad_spot("no such plate")
end

@testset "power files" begin
    q = nif_quads()
    pulses = Dict(:inner23 => PulseShape([0, 1, 5, 6], [0, 4, 8, 0]),
                  :inner30 => PulseShape([0, 2, 6], [0, 6, 0]),
                  :outer44 => PulseShape([0.5, 1, 5.5, 6.5], [0, 8, 8, 0]),
                  :outer50 => PulseShape([0, 3, 6], [0, 10, 0]))
    p = write_visrad_powers(joinpath(vrdir, "powers.csv"), q, pulses)
    @test startswith(readline(p), "time (ns),B")
    raw = read_visrad_powers(p)
    @test length(raw) == 192
    @test raw[:B115](3.0) ≈ pulses[q[findfirst(x -> x.name == "Q11B", q)].group](3.0) / 4
    rq = read_visrad_powers(p; template = nif_quads())
    ts = range(-1, 8; length = 91)
    @test all(isapprox(rq[Symbol(x.name)].(ts), pulses[x.group].(ts); atol = 1e-8) for x in q)
    # the readers key by Symbol, so the pulses compose with with_pulses by beam name
    @test [b.power for b in with_pulses(q, rq, 3.0; by = :beam)] ≈ [b.power for b in with_pulses(q, pulses, 3.0)]
    @test_throws KeyError with_pulses(nif_beams(), rq, 3.0; by = :beam)
    rt = write_visrad_powers(joinpath(vrdir, "rt.csv"), q, rq; by = :beam)
    @test all(isapprox(read_visrad_powers(rt; template = q)[Symbol(x.name)].(ts), rq[Symbol(x.name)].(ts); atol = 1e-8) for x in q)
    # the manual's examples: NIF "BL" names (import_tdep_laser_powers.jpg), OMEGA "Beam_10",
    # negative powers set to zero (import_tdep_beam_powers.html)
    pf = vrfile("tdep.csv", "time (ns),BL115,BL116,Beam_10\n0,0,0,0\n1,2.5,-0.1,0.5\n2,0,0,0\n")
    f = read_visrad_powers(pf)
    @test f[:B115](1.0) == 2.5 && f[:B116](1.0) == 0 && f[:Beam_10](0.5) == 0.25
    fb = read_visrad_powers(pf; template = nif_beams())
    @test keys(fb) == Set([:B115, :B116])
    @test pulse_energy(fb[:B115]) ≈ 2.5
    # seconds in the header; direct Dict form
    ps = vrfile("tdep_s.csv", "time (s),B111\n0,0\n1e-9,1\n2e-9,0\n")
    @test read_visrad_powers(ps)[:B111](1.0) ≈ 1
    pd = write_visrad_powers(joinpath(vrdir, "pd.csv"), Dict(:B111 => PulseShape([0, 1, 2], [0, 3, 0])))
    @test read_visrad_powers(pd)[:B111](0.5) ≈ 1.5
    # standard pulse-shape library file: W versus s
    pl = vrfile("pulse.dat", "# square\n0.0 0.0\n1.0e-10 5.0e11\n1.1e-9 5.0e11\n1.2e-9 0.0\n")
    sh = read_visrad_pulse_shape(pl)
    @test sh.t ≈ [0, 0.1, 1.1, 1.2] && sh(0.5) ≈ 0.5
end

@testset "result tables and fits" begin
    # element table as in surf_elem_table.jpg (whitespace, units in headers)
    pe = vrfile("elems.txt", """
        Z (cm) Azim. ang (deg) Radius (cm) Surf Index Emis T (eV) Rad T (eV) Flux In (TW/cm2) Flux Out (TW/cm2)
        1.0350e-01 1.5000e+01 7.8164e-02 0 6.2594e+01 9.5558e+01 8.5741e+00 1.5785e+00
        1.0350e-01 4.5000e+01 7.8164e-02 1 6.1132e+01 9.3326e+01 7.8007e+00 1.4361e+00
        1.0350e-01 7.5000e+01 7.8164e-02 2 6.1567e+01 9.3991e+01 8.0255e+00 1.4775e+00
        """)
    tb = read_visrad_table(pe)
    @test tb.header == ["Z (cm)", "Azim. ang (deg)", "Radius (cm)", "Surf Index", "Emis T (eV)",
                        "Rad T (eV)", "Flux In (TW/cm2)", "Flux Out (TW/cm2)"]
    @test size(tb.data) == (3, 8) && tb.data[2, 7] == 7.8007
    # power summary as in power_summary_table.jpg, object names as labels
    pp = vrfile("psd.txt", """
        Object\tSurface Area (cm2)\tIncident Rad Pwr\tEmitted Rad Pwr\tPwr Src Term\tLaser Pwr Dep
        hohlraum_barrel\t1.156e-01\t1.131e+00\t1.583e+00\t1.375e+00\t5.000e+00
        capsule\t5.027e-03\t4.894e-02\t3.426e-02\t0.000e+00\t0.000e+00
        Total (all objects)\t1.478e-01\t1.358e+00\t1.650e+00\t1.375e+00\t5.000e+00
        """)
    ts = read_visrad_table(pp)
    @test ts.labels == ["hohlraum_barrel", "capsule", "Total (all objects)"]
    @test ts.header[end] == "Laser Pwr Dep" && ts.data[1, end] == 5.0
    # Legendre and Y_LM fits: VisRad normalization (Y_00 = 1, A_L0 = A_L/√(2L+1))
    rng = Random.MersenneTwister(3)
    θ = acosd.(2 .* rand(rng, 400) .- 1)
    φ = 360 .* rand(rng, 400)
    f = 9.7 .+ 2.5 .* legendre.(2, cosd.(θ)) .+ 0.13 .* legendre.(4, cosd.(θ))
    A = visrad_legendre_fit(θ, f; terms = (0, 2, 4))
    @test A ≈ [9.7, 0, 2.5, 0, 0.13]
    Y = visrad_ylm_fit(θ, φ, f; lmax = 4)
    @test Y[mode_index(0, 0)] ≈ 9.7 && Y[mode_index(2, 0)] ≈ 2.5 / sqrt(5) && Y[mode_index(4, 0)] ≈ 0.13 / 3
    @test abs(1.143 * sqrt(5) - 2.563) / 2.563 < 0.005      # the example pair in the manual
    g = 1 .+ 0.2 .* sind.(θ) .* cosd.(φ) .- 0.1 .* sind.(θ) .^ 2 .* sind.(2 .* φ)
    Yg = visrad_ylm_fit(θ, φ, g; lmax = 2)
    @test Yg[mode_index(1, 1)] ≈ 0.2 / sqrt(3)               # Y_11 = √3 sinθ cosφ
    @test Yg[mode_index(2, -2)] ≈ -0.1 / (sqrt(2) * sqrt(5 / 24) * 3)
end

@testset "benchmark kit and comparison" begin
    t = hohlraum(nz = 8, nphi = 12, capsule_ntheta = 8, capsule_nphi = 12)
    q = filter(b -> b.group in (:inner30, :outer44), nif_quads())
    q = [remake(b; aim = (0, 0, sign(b.port[3]) * 4.0)) for b in q]
    pulses = Dict(:inner30 => PulseShape([0, 1, 3], [0, 6, 6]), :outer44 => PulseShape([0, 3], [0, 9]))
    kit = write_visrad_kit(joinpath(vrdir, "kit"), t, q; pulses, time = 2.0, notes = "test case")
    @test all(isfile, values(kit.files))
    @test kit.results["Laser delivered"] ≈ 6.0 * length(q)          # both pulses at 6 TW at t = 2 ns
    @test read_visrad_summary(kit.files[:summary]; template = q)[1].power ≈ 6.0
    readme = read(kit.files[:readme], String)
    @test occursin("cylinder: radius 2.875 mm", readme) && occursin("pointing_sst.csv", readme)
    # SST and summary in the kit reproduce the beams
    rq = read_visrad_sst(kit.files[:sst]; template = q)
    @test all(isapprox(a.aim, b.aim; atol = 1e-9) for (a, b) in zip(rq, q))
    # a VisRad element export made from the Briareus elements compares to zero
    _, er = read_csv(kit.files[:elements])
    io = IOBuffer()
    println(io, "Radius (cm)  Polar ang (deg)  Azim. ang (deg)  Flux In (TW/cm2)")
    for r in er
        println(io, "0.09  ", r[2], "  ", r[3], "  ", r[5])
    end
    pe = vrfile("kit_elems.txt", String(take!(io)))
    lp = kit.results["Laser Pwr Dep [wall]"]
    rows = compare_visrad(kit.dir; elements = pe,
                          power_summary = Dict("hohlraum_barrel" => Dict("Laser Pwr Dep" => lp)),
                          names = Dict("hohlraum_barrel" => "wall"))
    @test length(rows) == 1 + 1 + 2 + 4 + 1
    @test all(r -> abs(r.difference) < 1e-8 * max(1, abs(r.briareus)), rows)
    @test any(r -> r.quantity == "Laser Pwr Dep [wall]", rows)
end

@testset "workspace (.vrw)" begin
    # the synthetic workspace (vrw_fixture.jl): a primitive hohlraum + foil, a CAD mesh
    # placed through a reference cone, and NIF beam rows; lengths in cm
    tn, tf = TETRA_NODES, TETRA_FACES
    p = vrfile("work.vrw", synthetic_workspace())

    ws = read_visrad_workspace(p)
    @test ws.version == "21.0.0"
    @test [o.kind for o in ws.objects] == [:rectangle, :cylindrical_hohlraum, :cone, :mesh, :sphere]
    @test [o.folder for o in ws.objects] == ["Foil", "", "CAD", "CAD", ""]
    @test ws.clusters == ["Ring_050", "Ring_045"]
    # spherical position (r, θ, φ) in cm
    @test ws.objects[5].position ≈ 3.2 * port_vector(90, 150)
    # the CAD mesh: nodes (cm, centered) placed by its position in the reference cone's
    # frame, which is turned 90° about z (Azimuthal) and raised 2 mm
    cad = ws.objects[4]
    Rz = [0 -1 0; 1 0 0; 0 0 1]
    @test cad.position ≈ Vec3(0, 1, 2)
    @test cad.mesh.vertices ≈ [Vec3(Rz * (10 .* collect(v))) + Vec3(0, 1, 2) for v in tn]
    @test cad.mesh.triangles == [(1, 3, 2), (1, 2, 4), (1, 4, 3), (2, 3, 4)]
    @test cad.albedo == 0.5 && cad.radiosity && !ws.objects[3].radiosity
    # The VisRad element counts for these primitives (power-source table of the
    # workspace they come from): 81 for the rectangle, 216 for the hohlraum
    t = visrad_target(ws)
    @test [c.name for c in t.components] == [:Rectangle_0, :Cyl_Hohlraum_0_wall, :Cyl_Hohlraum_0_cap_top,
        :Cyl_Hohlraum_0_cap_bottom, :Cyl_Hohlraum_0_leh_top, :Cyl_Hohlraum_0_leh_bottom, :CAD_0]
    np = [npatches(c.mesh) for c in t.components]
    @test np[1] == 81 && sum(np[2:6]) == 216 && np[7] == 4
    # materials by object name, a String or a Symbol; a name that is no object is an error
    tm = visrad_target(ws; materials = Dict(:Rectangle_0 => gold()))
    @test tm.components[1].material == gold() && tm.components[7].material == t.components[7].material
    @test_throws ArgumentError visrad_target(ws; materials = Dict(:Rectangle => gold()))
    # foil: 4 × 4 mm at x = 1 mm, tilted 90° so its +z normal turns to +x
    e = t.comp_elems[1]
    @test all(n -> n ≈ Vec3(1, 0, 0), t.elem_normal[e])
    @test sum(t.elem_area[e]) ≈ 16
    @test sum(t.elem_area[i] * t.elem_centroid[i] for i in e) / 16 ≈ Vec3(1, 0, 0)
    # hohlraum: radius 2.55, length 9.35, LEH radius 1.52 mm, centered at x = -4 mm
    w = t.comp_elems[2]
    @test sum(t.elem_area[w]) ≈ 2π * 2.55 * 9.35 rtol = 1e-12        # exact (curved) areas
    @test t.components[5].kind === :aperture
    @test all(norm((t.elem_centroid[i] - Vec3(-4, 0, 0))[1:2]) ≈ 2.55 * cos(π / 12) for i in w)  # 12 chords
    @test all(dot(t.elem_normal[i], t.elem_centroid[i] - Vec3(-4, 0, 0)) < 0 for i in w)  # normals in
    # primitives at a given element size instead of VisRad's grid: same surfaces and places
    rect = only(visrad_component(ws.objects[1]; primitive_h = 0.5))
    @test npatches(rect.mesh) == 8 * 8                                   # 4 mm / 0.5 mm
    tr = Target([rect])
    @test sum(tr.elem_area) ≈ 16 && all(n -> n ≈ Vec3(1, 0, 0), tr.elem_normal)
    hc = visrad_component(ws.objects[2]; primitive_h = 0.5)
    @test [c.name for c in hc] == [c.name for c in t.components[2:6]]
    cnt = hohlraum_counts(; length = 9.35, radius = 2.55, leh_radius = 1.52, h = 0.5)
    @test npatches(hc[1].mesh) == cnt.nz * cnt.nphi && npatches(hc[2].mesh) == cnt.nr_cap * cnt.nphi
    @test npatches(hc[4].mesh) == cnt.nr_leh * cnt.nphi && hc[4].kind === :aperture
    th = Target(hc)
    @test sum(th.elem_area[th.comp_elems[1]]) ≈ 2π * 2.55 * 9.35 rtol = 1e-12
    @test all(dot(th.elem_normal[i], th.elem_centroid[i] - Vec3(-4, 0, 0)) < 0 for i in th.comp_elems[1])
    @test sum(th.elem_area[th.comp_elems[2]]) ≈ sum(t.elem_area[t.comp_elems[3]]) rtol = 1e-12
    sph = only(visrad_component(ws.objects[5]; primitive_h = 0.2))
    @test npatches(sph.mesh) > npatches(only(visrad_component(ws.objects[5])).mesh)
    @test_throws ArgumentError visrad_component(ws.objects[1]; primitive_h = 0.0)
    # CAD mesh: outward front faces; two-sided on request
    @test all(dot(t.elem_normal[i], t.elem_centroid[i] - cad.position) > 0 for i in t.comp_elems[7])
    t2 = visrad_target(ws; select = o -> o.kind === :mesh, two_sided = o -> true)
    @test length(t2.components) == 1 && length(t2.comp_elems[1]) == 8
    @test !t2.closed && visrad_target(ws; select = o -> o.kind === :mesh, closed = true).closed
    # remeshed into elements of about 0.3 mm (its faces are 1 mm right triangles): the
    # same surface, still facing out; `coarsen` and `remesh` exclude each other
    t3 = visrad_target(ws; select = o -> o.kind === :mesh, orient = true, toward = cad.position, remesh = 0.3)
    @test sum(t3.elem_area) ≈ sum(t.elem_area[t.comp_elems[7]]) rtol = 1e-12
    @test nelements(t3) > 4 * 4
    @test all(dot(t3.elem_normal[i], t3.elem_centroid[i] - cad.position) > 0 for i in 1:nelements(t3))
    @test_throws ArgumentError visrad_target(ws; select = o -> o.kind === :mesh, coarsen = 0.3, remesh = 0.3)
    # a tetrahedron has no cylinder or sphere: `fit` (on by default with `remesh`) leaves the
    # part whole, as `fit = false` does
    t4 = visrad_target(ws; select = o -> o.kind === :mesh, orient = true, toward = cad.position, remesh = 0.3,
                       fit = false)
    @test [c.name for c in t4.components] == [c.name for c in t3.components] && t4.elem_area == t3.elem_area
    # beams
    b1, b2, b3 = ws.beams
    @test [b.name for b in ws.beams] == ["B113", "B111", "B121"]
    @test b1.group === :Ring_050 && b3.group === :Ring_045 && b1.fnumber == 20.755
    @test port_angles(b1)[1] ≈ 52.382 && port_angles(b1)[2] ≈ 37.056
    @test b1.aim ≈ Vec3(-0.245, 0.3775, 5.905) && b1.defocus ≈ -0.1
    # power column 1 TW under the "Fraction of Total Beam Power" red/green model with
    # f_G = f_R = 0.05 (green off): the blue cone carries 0.9 TW; under "None" all of it
    @test b1.power ≈ 0.9 && b2.power == 1.0 && b3.power == 0.0
    # Rev1-Outer-2w radii are HWHM: 1/e semi-axes r/(ln 2)^(1/n); rotation 0 puts the
    # major axis along the φ̂ of the beam's quad (from its name without a template)
    s = log(2)^(1 / 4.5)
    @test b1.spot.a ≈ 0.62 / s && b1.spot.b ≈ 0.36 / s && b1.spot.order == 4.5
    @test b1.quad == "Q11T" && b3.quad == "Q12T"
    @test mod(b1.spot.angle - _vr_reference_angle(b1), π) ≈ π / 2 && abs(_vr_reference_angle(b1)) > 0.01
    # a custom super-Gaussian is read at 1/e; pointing in a component's frame (53-field
    # row with a rotation ID) lands on that component's origin
    @test b2.spot.a ≈ 0.3 && b2.spot.order == 2.0
    @test b2.aim ≈ cad.position
    @test ws.rows[2].rot_id == "R1" && ws.rows[1].rot_id == ""
    # template: ports match the NIF table, facility metadata kept, absent beams off
    nifb = nif_beams()
    wt = read_visrad_workspace(p; template = nifb)
    nb = Dict(b.name => b for b in wt.beams)
    @test length(wt.beams) == length(nifb)
    @test nb["B113"].aim ≈ b1.aim && nb["B113"].power ≈ 0.9
    @test nb["B113"].quad == nifb[findfirst(b -> b.name == "B113", nifb)].quad
    @test nb["B113"].spot.angle ≈ b1.spot.angle atol = 1e-4     # the file's ports carry 3 decimals
    @test nb["B121"].power == 0 && nb["B468"].power == 0
    # VisRad's beam model with the file's envelope ratios (per beam)
    @test visrad_profile(ws).envelope["B113"] == 1.53645
    @test_throws ArgumentError read_visrad_workspace(vrfile("not.vrw", "{\"json\": 1}"))

    # written back into the workspace: new pointing, focus, power and on/off, in the
    # frame of each row; everything else unchanged
    edited = [remake(b1; aim = Vec3(0.3, -0.2, 5.5), defocus = 0.25, power = 2.5),
              remake(b2; aim = cad.position + Vec3(0.1, 0.0, 0.0)),        # in CAD_0's frame
              remake(b3; power = 1.5)]                                       # was off
    q = vrfile("back.vrw", "")
    write_visrad_workspace(q, p, edited)
    wb = read_visrad_workspace(q)
    for (x, y) in zip(edited, wb.beams)
        @test y.name == x.name
        @test y.aim ≈ x.aim atol = 1e-6
        @test y.defocus ≈ x.defocus atol = 1e-6
        @test y.power ≈ x.power rtol = 1e-4
    end
    @test wb.rows[3].on_blue == 1 && wb.rows[2].rot_id == "R1"
    @test wb.rows[1].power_blue ≈ 2.5 / 0.9 rtol = 1e-4                  # the total under model 2
    @test wb.rows[2].power_blue == 1.0
    old, new = readlines(p), readlines(q)
    changed = findall(k -> old[k] != new[k], eachindex(old))
    @test length(old) == length(new) && length(changed) == 3              # the three rows only
    @test all(k -> length(split(new[k])) == length(split(old[k])), changed)
    @test_throws ArgumentError write_visrad_workspace(q, p, [remake(b1; port = port_vector(60, 37))])
    @test_logs (:warn, r"no row") write_visrad_workspace(q, p, [edited[1]; Beam("B999"; theta = 10, phi = 0)])

    # VisRad's results read back: the [Power Source Data] section of a computed
    # workspace (vrw_fixture.jl), followed by a table of the same layout in another
    # section, which is not read
    qcad = [10.0, 0.0, 2.5, 7.25]
    qrect = Float64.(1:81)                         # VisRad's row number, to find each row
    tables(; rect = qrect, hohl = zeros(216)) =
        ["Rectangle_0" => rect, "Cyl_Hohlraum_0" => hohl, "CAD_0" => qcad]
    decoy = replace(power_source_section(["Decoy" => [1.0]]), "[Power Source Data]" => "[Beam Capture-Clearance Data]",
                    "[End Power Source Data]" => "[End Beam Capture-Clearance Data]")
    computed(; kw...) = vrfile("computed.vrw", synthetic_workspace() * "\n" * power_source_section(tables(; kw...)) * decoy)
    pc = computed()
    ps = visrad_power_sources(pc)
    @test ps.time == 0 && sort(collect(keys(ps.tables))) == ["CAD_0", "Cyl_Hohlraum_0", "Rectangle_0"]
    tc = ps.tables["CAD_0"]
    @test keys(tc) == (:surface_index, :direct_laser_pwr_source, :uniform_source, :albedo, :x_ray_conv_efficiency,
                       :laser_reflectivity, :pre_reflect_laser_pwr, :n_red_cones_hitting_surf,
                       :n_blue_cones_hitting_surf, :min_angle_overlap_1w, :min_angle_overlap_3w)
    @test tc.surface_index == 0:3 && tc.surface_index isa Vector{Int} && tc.n_blue_cones_hitting_surf == [2, 0, 2, 2]
    @test tc.direct_laser_pwr_source == qcad && tc.albedo == fill(0.7, 4)
    @test ps.units[:direct_laser_pwr_source] == "TW/cm2" && ps.units[:pre_reflect_laser_pwr] == "TW/cm2"
    @test ps.units[:min_angle_overlap_3w] == "deg" && ps.units[:albedo] == ""
    @test length(ps.tables["Cyl_Hohlraum_0"].albedo) == 216
    # not computed (or changed since): no tables
    ps0 = vrfile("stale.vrw", synthetic_workspace() * "\n" * power_source_section([]; up_to_date = false))
    @test_throws ArgumentError visrad_power_sources(ps0)
    @test visrad_power_sources(ps0; stale = :nothing) === nothing
    @test_throws ArgumentError visrad_power_sources(p)                     # no section
    # beams written into a computed workspace: its results belong to the old beams
    qc = joinpath(dirname(pc), "computed_rewritten.vrw")
    write_visrad_workspace(qc, pc, edited)
    @test visrad_power_sources(qc; stale = :nothing) === nothing
    @test occursin("Power sources are not up to date for time", read(qc, String))
    @test all(isapprox(y.aim, x.aim; atol = 1e-6) for (x, y) in zip(edited, read_visrad_workspace(qc).beams))

    # mapped onto the target: power = source × area (TW/cm² × mm² / 100) per element
    wc = read_visrad_workspace(pc)
    tc = visrad_target(wc)
    P = visrad_deposition(wc, tc)
    @test length(P) == nelements(tc) && P == visrad_deposition(pc)
    ec = elements(tc, :CAD_0)
    @test P[ec] ≈ qcad .* tc.elem_area[ec] ./ 100
    @test all(iszero, P[elements(tc, :Cyl_Hohlraum_0_wall, :Cyl_Hohlraum_0_cap_top)])
    # the rectangle's rows run its local x fastest: tilted 90°, local x is chamber -z
    # and local y is chamber y; 9 × 9 elements of 4/9 mm
    er = elements(tc, :Rectangle_0)
    vrow(c) = 9 * round(Int, (c[2] + 2) * 9 / 4 - 0.5) + round(Int, (-c[3] + 2) * 9 / 4 - 0.5) + 1
    @test P[er] ≈ [vrow(tc.elem_centroid[e]) * tc.elem_area[e] / 100 for e in er]
    @test sort(vrow.(tc.elem_centroid[er])) == 1:81
    # two-sided: on the front elements
    t2 = visrad_target(wc; select = o -> o.kind === :mesh, two_sided = true)
    P2 = @test_logs (:warn, r"Rectangle_0") match_mode = :any visrad_deposition(wc, t2)
    @test P2[t2.elem_side .== 1] ≈ qcad .* t2.elem_area[1:4] ./ 100 && all(iszero, P2[t2.elem_side .== -1])
    # VisRad's hohlraum element order is unknown: refused when it has laser power,
    # unless left out of the target; elements not from VisRad are refused
    wh = read_visrad_workspace(computed(; hohl = ones(216)))
    @test_throws ArgumentError visrad_deposition(wh)
    th = visrad_target(wh; select = o -> o.kind in (:rectangle, :mesh))
    @test (@test_logs (:warn, r"Cyl_Hohlraum_0") visrad_deposition(wh, th)) == P[[er; ec]]
    @test_throws ArgumentError visrad_deposition(wc, visrad_target(wc; select = o -> o.kind === :mesh, remesh = 0.3))
    @test_throws ArgumentError visrad_deposition(wc, Target(visrad_component(wc.objects[1]; primitive_h = 0.5)))

    # the writer never overwrites VisRad's results: a computed text workspace (also as the
    # template) or a JSON one is refused unless `overwrite = true`; a JSON template is refused
    @test_throws ArgumentError write_visrad_workspace(pc, p, edited)
    @test_throws ArgumentError write_visrad_workspace(pc, pc, edited)
    qo = vrfile("overwritten.vrw", read(pc, String))
    write_visrad_workspace(qo, pc, edited; overwrite = true)
    @test visrad_power_sources(qo; stale = :nothing) === nothing
    write_visrad_workspace(q, p, edited)                                   # an uncomputed file: rewritten as before
    @test read_visrad_workspace(q).beams[1].aim ≈ edited[1].aim atol = 1e-6
end

@testset "workspace (.vrw): VisRad 21's JSON format" begin
    # The same workspace in the JSON form VisRad 21 saves a computed workspace in
    # (vrw_fixture.jl), its ports at full precision: every object, mesh, folder, beam row
    # and beam is the text reader's, apart from the ports' digits.
    p = vrfile("work.vrw", synthetic_workspace())
    wt = read_visrad_workspace(p)
    full = Dict(r.name => (r.theta + 2.3e-4, r.phi - 4.1e-4, r.fnumber + 3.3e-4) for r in wt.rows)
    qrect = Float64.(1:81)
    qcad = [10.0, 0.0, 2.5, 7.25]
    src = Dict("Rectangle_0" => qrect, "Cyl_Hohlraum_0" => zeros(216), "CAD_0" => qcad)
    pj = vrfile("work_json.vrw", synthetic_json_workspace(wt; sources = src, ports = full))
    @test Briareus._vrw_is_json(pj) && !Briareus._vrw_is_json(p)
    wj = read_visrad_workspace(pj)
    @test wj.version == wt.version && wj.clusters == wt.clusters
    @test [o.name for o in wj.objects] == [o.name for o in wt.objects]
    @test [(o.kind, o.type, o.folder) for o in wj.objects] == [(o.kind, o.type, o.folder) for o in wt.objects]
    @test all((o.albedo, o.conversion, o.reflectivity, o.radiosity, o.laser, o.visible, o.color) ==
              (u.albedo, u.conversion, u.reflectivity, u.radiosity, u.laser, u.visible, u.color) for (o, u) in zip(wj.objects, wt.objects))
    @test all(isapprox(o.position, u.position; atol = 1e-12) && isapprox(o.rotation, u.rotation; atol = 1e-12)
              for (o, u) in zip(wj.objects, wt.objects))
    cj, ct = wj.objects[4].mesh, wt.objects[4].mesh
    @test cj.triangles == ct.triangles && cj.patches == ct.patches && isapprox(cj.vertices, ct.vertices; atol = 1e-12)
    # the text table rounds the ports (and cone_max) to 3 decimals, and lists a beam's power
    # table even when empty (the JSON reader keeps only non-empty ones)
    tabs(r) = filter(!isempty ∘ last, r.tables)
    same(r, s) = tabs(r) == tabs(s) &&
                 all(k -> k in (:theta, :phi, :fnumber, :cone_max, :tables) || isequal(r[k], s[k]), keys(s))
    @test length(wj.rows) == length(wt.rows) && all(same(r, s) for (r, s) in zip(wj.rows, wt.rows))
    @test all((r.theta, r.phi, r.fnumber) == full[r.name] for r in wj.rows)        # the JSON file's digits
    @test all(r.cone_max ≈ atand(1 / (2r.fnumber)) && r.cone_min == 0 for r in wj.rows)
    @test all(x.name == y.name && isapprox(x.aim, y.aim; atol = 1e-12) && x.power == y.power &&
              x.defocus == y.defocus for (x, y) in zip(wj.beams, wt.beams))
    # with a template too, as for text
    @test [b.power for b in read_visrad_workspace(pj; template = nif_beams()).beams] ==
          [b.power for b in read_visrad_workspace(p; template = nif_beams()).beams]
    # the same targets, so VisRad's results map onto either
    tj, tt = visrad_target(wj), visrad_target(wt)
    @test [c.name for c in tj.components] == [c.name for c in tt.components] && tj.elem_area ≈ tt.elem_area
    # power sources: per object in the components' order, the text reader's columns and units
    ps = visrad_power_sources(pj)
    @test ps.time == 0 && sort(collect(keys(ps.tables))) == sort([o.name for o in wt.objects])
    tc = ps.tables["CAD_0"]
    @test keys(tc) == (:surface_index, :direct_laser_pwr_source, :uniform_source, :albedo, :x_ray_conv_efficiency,
                       :laser_reflectivity, :pre_reflect_laser_pwr, :n_red_cones_hitting_surf,
                       :n_blue_cones_hitting_surf, :min_angle_overlap_1w, :min_angle_overlap_3w)
    @test tc.direct_laser_pwr_source == qcad && tc.surface_index == 0:3 && tc.n_blue_cones_hitting_surf == [2, 0, 2, 2]
    @test tc.n_blue_cones_hitting_surf isa Vector{Int}
    @test ps.units[:direct_laser_pwr_source] == "TW/cm2" && ps.units[:min_angle_overlap_3w] == "deg" && ps.units[:albedo] == ""
    # mapped onto the target exactly as the text format's tables
    pt = vrfile("work_text_computed.vrw", synthetic_workspace() * "\n" *
                power_source_section(["Rectangle_0" => qrect, "Cyl_Hohlraum_0" => zeros(216), "CAD_0" => qcad]))
    @test visrad_deposition(wj, tj) == visrad_deposition(read_visrad_workspace(pt), tt)
    # not computed: no tables
    ps0 = vrfile("stale_json.vrw", synthetic_json_workspace(wt))
    @test_throws ArgumentError visrad_power_sources(ps0)
    @test visrad_power_sources(ps0; stale = :nothing) === nothing
    # the writer patches text only, and never writes over a JSON workspace
    @test_throws ArgumentError write_visrad_workspace(vrfile("from_json.vrw", ""), pj, wj.beams)
    @test_throws ArgumentError write_visrad_workspace(pj, p, wt.beams)
    @test Briareus._vrw_is_json(pj)                                            # left as it was
    # JSON that is no VisRad workspace, and a mismatched power-source count
    @test_throws ArgumentError read_visrad_workspace(vrfile("other.json", "{\"json\": 1}"))
    bad = replace(read(pj, String), r"\"Number of target components\": \d+" => "\"Number of target components\": 9")
    @test_throws ArgumentError read_visrad_workspace(vrfile("bad.vrw", bad))
end

@testset "workspace (.vrw): inconsistent CAD windings" begin
    # the tetrahedron with one face wound against the others: read one-sided as stored,
    # visrad_target warns (that face's area share); orient = true and two_sided do not
    txt = synthetic_workspace()
    a = findfirst("Surface element connectivity:", txt)
    b = findnext("End Parameters", txt, last(a))
    conn = txt[last(a)+1:first(b)-1]
    flipped = replace(conn, r"\n(\s+)3(\s+)1(\s+)2(\s+)3(\s+)-1" => s"\n\g<1>3\g<2>1\g<3>3\g<4>2\g<5>-1")
    @test flipped != conn
    w = read_visrad_workspace(vrfile("miswound.vrw", txt[1:last(a)] * flipped * txt[first(b):end]))
    sel = o -> o.kind === :mesh
    f = Briareus._vrw_miswound(only(o for o in w.objects if o.kind === :mesh))
    @test f ≈ (sqrt(3) / 2 * 0.01) / (3 * 0.005 + sqrt(3) / 2 * 0.01) rtol = 1e-6 # the slanted face's share
    @test_logs (:warn, r"1 VisRad CAD meshes are wound inconsistently.*CAD_0 36\.6%") visrad_target(w; select = sel)
    @test_logs visrad_target(w; select = sel, orient = true, toward = w.objects[4].position)
    @test_logs visrad_target(w; select = sel, two_sided = true)
    @test Briareus._vrw_miswound(only(o for o in read_visrad_workspace(vrfile("w.vrw", txt)).objects if o.kind === :mesh)) == 0
end

@testset "workspace (.vrw): required and optional keys" begin
    txt = synthetic_workspace()
    edit(f) = vrfile("keys.vrw", f(txt))
    # in object `name`'s block, the line of `key` replaced (`nothing`: removed)
    function rekey(name, key, value)
        return t -> begin
            a = last(findfirst("Begin Parameters for Object = $name", t))
            b = first(findnext("End Parameters for Object: $name", t, a))
            blk = replace(t[a+1:b-1], Regex("\\n\\s*" * replace(key, "." => "\\.") * "\\s*=[^\\n]*") =>
                          value === nothing ? "" : "\n   $key = $value")
            t[1:a] * blk * t[b:end]
        end
    end
    msg(f) = try f(); "" catch e; e isa ArgumentError ? e.msg : rethrow() end
    # a missing radius: refused when the primitive is meshed, naming key and object
    ws = read_visrad_workspace(edit(rekey("Cyl_Hohlraum_0", "Length 1 value", nothing)))
    m = msg(() -> visrad_target(ws))
    @test occursin("Cyl_Hohlraum_0", m) && occursin("Length 1 value", m)
    # a malformed or missing position, a missing albedo: refused on reading
    m = msg(() -> read_visrad_workspace(edit(rekey("Marker", "Position 2", "abc"))))
    @test occursin("Marker", m) && occursin("Position 2", m) && occursin("abc", m)
    m = msg(() -> read_visrad_workspace(edit(rekey("CAD_0", "Position 1", nothing))))
    @test occursin("CAD_0", m) && occursin("Position 1", m)
    @test occursin("Albedo value", msg(() -> read_visrad_workspace(edit(rekey("Rectangle_0", "Albedo value", nothing)))))
    # a mesh count that is not an integer
    ws = read_visrad_workspace(edit(rekey("Rectangle_0", "Number of mesh pts 1", "9.5")))
    @test occursin("not an integer", msg(() -> visrad_target(ws)))
    # optional keys take their defaults: color, visibility, normals
    ws0 = read_visrad_workspace(vrfile("keys0.vrw", txt))
    ws = read_visrad_workspace(edit(t -> rekey("Rectangle_0", "Normals in/out", nothing)(
                                          rekey("Rectangle_0", "Red value", nothing)(
                                          rekey("Rectangle_0", "Object visibility", nothing)(t)))))
    r, r0 = ws.objects[1], ws0.objects[1]
    @test r.color == (0.7, 0.5, 0.0) && r.visible && r.position == r0.position
    @test visrad_target(ws).elem_normal == visrad_target(ws0).elem_normal
end

@testset "workspace results: element tables, radiosity and modes" begin
    # a CAD capsule (outward sphere, r = 0.5 mm, 1 mm below the chamber center) beside
    # the tetrahedron, which takes VisRad's laser power
    sm = sphere_mesh(0.5; ntheta = 16, nphi = 24)
    nodes = [Tuple(v ./ 10) for v in sm.vertices]
    faces = [(a - 1, b - 1, c - 1) for (a, b, c) in sm.triangles]
    nc = length(faces)
    qcad = [10.0, 0.0, 2.5, 7.25]
    txt = synthetic_workspace(; meshes = [("Capsule_0", nodes, faces, (0.0, 0.0, -0.1))]) * "\n" *
          power_source_section(["Rectangle_0" => zeros(81), "Cyl_Hohlraum_0" => zeros(216), "CAD_0" => qcad,
                                "Capsule_0" => zeros(nc)])
    ws = read_visrad_workspace(vrfile("results.vrw", txt))
    t = visrad_target(ws; select = o -> o.kind === :mesh)
    sol = solve(RadiationModel(t), visrad_deposition(ws, t))
    ec, es = elements(t, :CAD_0), elements(t, :Capsule_0)
    # VisRad's capsule export, rows in reverse order (0-based surface index): an
    # axisymmetric field f₀(1 + a₁P₁ + a₂P₂) about the capsule center at each centroid
    f0, a1, a2 = 0.5, 0.04, 0.1
    μ = [normalize(t.elem_centroid[e] - Vec3(0, 0, -1))[3] for e in es]
    f = f0 .* (1 .+ a1 .* legendre.(1, μ) .+ a2 .* legendre.(2, μ))
    io = IOBuffer()
    println(io, "Radius (cm)  Surf Index  Rad T (eV)  Flux In (TW/cm2)")
    for k in reverse(1:nc)
        println(io, "0.05  ", k - 1, "  ", radiation_temperature(f[k] / 100), "  ", f[k])
    end
    pcap = vrfile("capsule.sed", String(take!(io)))
    # the tetrahedron's export: 1-based indices, reversed, Flux Out 2% above Briareus's
    B = sol.B[ec]
    io = IOBuffer()
    println(io, "Surf Index  Flux Out (TW/cm2)")
    for k in 4:-1:1
        println(io, k, "  ", 102 * B[k])
    end
    ptet = vrfile("tetra.sed", String(take!(io)))

    V = visrad_element_values(ws, t, Dict("Capsule_0" => pcap), r"flux\s*in"i)
    @test V[es] ≈ f ./ 100 && all(isnan, V[ec])                     # TW/cm² → TW/mm², rows by index
    @test visrad_element_values(ws, t, Dict(:CAD_0 => ptet), r"flux\s*out"i)[ec] ≈ 1.02 .* B
    @test visrad_element_values(ws, t, Dict("Capsule_0" => pcap), r"rad\s*t"i)[es] ≈ radiation_temperature.(f ./ 100)

    out = IOBuffer()
    cmp = compare_visrad(ws, sol; elements = Dict("Capsule_0" => pcap, "CAD_0" => ptet), observe = :Capsule_0, io = out)
    row(q) = only(r for r in cmp.rows if r.quantity == q)
    stat(o, q) = only(s for s in cmp.elements if s.object == o && s.quantity == q)
    # modes of VisRad's field by exact solid angles: the Legendre ratios of the field
    @test abs(row("P1/P0 [Capsule_0]").visrad - 100a1) < 0.02
    @test abs(row("P2/P0 [Capsule_0]").visrad - 100a2) < 0.05
    @test abs(row("mode amplitude l=2 [Capsule_0]").visrad - 100a2) < 0.05
    # (centroid samples held constant over 16 rings leak 0.2 pp of the P2 into P4)
    @test all(abs(row("mode amplitude l=$l [Capsule_0]").visrad) < 0.3 for l in (3, 4))
    @test row("Rad T [Capsule_0]").visrad ≈ radiation_temperature(f0 / 100) rtol = 1e-3
    @test row("Incident Rad Pwr [Capsule_0]").briareus ≈ sum(sol.H[es] .* t.elem_area[es]) rtol = 1e-12
    # the tetrahedron: 2% above everywhere, so the rms is 0.02 B over the mean of 1.02 B
    s = stat("CAD_0", "Flux Out")
    A = t.elem_area[ec]
    @test s.n == 4
    @test s.rms ≈ sqrt(sum(A .* (0.02B) .^ 2) / sum(A)) / (sum(A .* 1.02B) / sum(A)) rtol = 1e-12
    @test s.max ≈ maximum(0.02B) / (sum(A .* 1.02B) / sum(A)) rtol = 1e-12
    @test s.correlation ≈ 1
    @test row("Emitted Rad Pwr [CAD_0]").relative ≈ 0.02 rtol = 1e-9
    # laser deposition from the power sources: from VisRad here, so no difference
    @test row("Laser Pwr Dep [CAD_0]").visrad ≈ sum(qcad .* A ./ 100) rtol = 1e-12
    @test row("Laser Pwr Dep [CAD_0]").difference == 0 && stat("CAD_0", "Laser").rms == 0
    @test !any(r -> occursin("Rectangle_0", r.quantity), cmp.rows)               # not in the target
    @test occursin("P2/P0 [Capsule_0]", String(take!(out)))
    # the index base: 1-based found from the largest index; ambiguous tables are refused
    part = (; header = ["Surf Index", "Flux Out (TW/cm2)"], labels = fill("", 3), data = [1.0 1; 2 2; 3 3])
    @test_throws ArgumentError visrad_element_values(ws, t, Dict("CAD_0" => part), r"flux\s*out"i)
    vp = visrad_element_values(ws, t, Dict("CAD_0" => part), r"flux\s*out"i; index_base = 1)[ec]
    @test vp[1:3] ≈ [0.01, 0.02, 0.03] && isnan(vp[4])
    @test_throws ArgumentError visrad_element_values(ws, t, Dict("CAD_0" => (; header = ["Flux Out (TW/cm2)"], labels = ["", ""],
                                                                             data = [1.0; 2;;])), r"flux\s*out"i)
    @test_throws ArgumentError visrad_element_values(ws, t, Dict("Nothing" => ptet), r"flux\s*out"i)
    @test_throws ArgumentError compare_visrad(ws, sol; elements = Dict("Capsule_0" => (; header = ["Surf Index", "Flux In (TW/cm2)"],
                                                                                       labels = ["", ""], data = [0.0 1; 1 1])),
                                              observe = :Capsule_0)
end

end
