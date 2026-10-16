# NIF TaLIS configuration checks (src/talis.jl). Included inside the top-level testset.

@testset "TaLIS" begin
    trow(r, rule, c; beam = nothing) =
        only(x for x in r.rows if x.rule === rule && x.component === c && (beam === nothing || x.beam == beam))
    tan05 = tand(0.5)
    R = 0.05 + 0.025                                   # default pointing + target offset

    # 3ω clearance, analytic: beam straight down onto a foil at z = 0; a plate upstream at
    # z = 2 whose edge runs along x = x0. Worst case = x0 − (1% contour + f-cone at s = −2)
    # − transverse errors − focus/(2F) − r·tan 0.5° at the edge point nearest the axis.
    spot = SuperGaussianSpot(0.3; order = 4.0)
    b = Beam("b"; theta = 0, phi = 0, power = 1.0, spot, fnumber = 8.0)
    x0, w = 1.6, 2.0
    plate = Component(:plate, translate(rectangle_mesh(w, w; nu = 4, nv = 4), (x0 + w / 2, 0, 2.0)), gold();
                      two_sided = true)
    t = Target(foil(:f, (0, 0, 0), (0, 0, 1); width = 4.0, n = 8), plate)
    r = talis_check(t, [b]; capture = (:f,))
    ρ1 = log(100)^(1 / 4) * 0.3                         # 1%-of-peak contour of exp(-ρ⁴)
    nominal = x0 - (ρ1 + 2 / 16)
    worst = nominal - R - 1 / 16 - tan05 * hypot(x0, 2.0)
    c3 = trow(r, :clearance_3w, :plate)
    @test c3.nominal ≈ nominal atol = 1e-9
    @test c3.worst ≈ worst atol = 1e-9
    @test c3.margin ≈ worst - 0.5 atol = 1e-9
    @test c3.point ≈ Vec3(x0, 0, 2.0)
    @test c3.status === :pass && c3.unit === :mm
    # the tolerance term: dropping it, or moving the pivot, changes only r·tan(0.5°)
    c0 = trow(talis_check(t, [b]; capture = (:f,), tolerance = 0), :clearance_3w, :plate)
    @test c0.worst - c3.worst ≈ tan05 * hypot(x0, 2.0) atol = 1e-9
    bA = remake(b; aim = (0.0, 0.0, 0.5))               # aim 0.5 mm up the axis
    cA = trow(talis_check(t, [bA]; capture = (:f,), pivot = :aim), :clearance_3w, :plate)
    @test cA.worst ≈ x0 - (ρ1 + 1.5 / 16) - R - 1 / 16 - tan05 * hypot(x0, 1.5) atol = 1e-9
    cP = trow(talis_check(t, [b]; capture = (:f,), pivot = (x0, 0.0, 2.0)), :clearance_3w, :plate)
    @test cP.worst ≈ nominal - R - 1 / 16 atol = 1e-9    # no allowance at the pivot itself
    cM = trow(talis_check(t, [bA]; capture = (:f,), pivot = :max), :clearance_3w, :plate)
    @test cM.worst <= min(cA.worst, trow(talis_check(t, [bA]; capture = (:f,)), :clearance_3w, :plate).worst) + 1e-12
    # moving the plate in fails the rule, with the blow-by consequence
    tin = Target(foil(:f, (0, 0, 0), (0, 0, 1); width = 4.0, n = 8),
                 translate(plate, (-0.5, 0.0, 0.0)))
    cf = trow(talis_check(tin, [b]; capture = (:f,)), :clearance_3w, :plate)
    @test cf.status === :fail && cf.margin < 0 && occursin("blow-by", cf.note)

    # 3ω capture, analytic: the 10% contour sits 2 − xa − ρ10 from the foil edge at x = 2.
    ρ10 = log(10)^(1 / 4) * 0.3
    for xa in (0.0, 0.6, 1.0)
        rc = talis_check(t, [remake(b; aim = (xa, 0.0, 0.0))]; capture = (:f,))
        cp = trow(rc, :capture_3w, :f)
        @test cp.nominal ≈ 2 - xa - ρ10 atol = 1e-9
        @test cp.worst ≈ 2 - xa - ρ10 - R - 1 / 16 - tan05 * 2 atol = 1e-9
        @test cp.status === (cp.worst >= 0.5 ? :pass : :fail)
        xa > 0 && @test cp.point ≈ Vec3(2, 0, 0) atol = 1e-9
    end
    @test trow(talis_check(t, [remake(b; aim = (1.0, 0.0, 0.0))]; capture = (:f,)), :capture_3w, :f).status === :fail
    @test trow(talis_check(t, [remake(b; aim = (0.9, 0.0, 0.0))]; capture = (:f,)), :capture_3w, :f).status === :pass
    # an axis that misses the substrate fails capture outright
    rmiss = talis_check(t, [remake(b; aim = (3.0, 0.0, 0.0))]; capture = (:f,))
    @test trow(rmiss, :capture_3w, :none).status === :fail
    # incidence: flat foil, beam at θ from its normal; worst = θ + 0.5° tolerance
    big = Target(foil(:f, (0, 0, 0), (0, 0, 1); width = 12.0, n = 12))
    for (θ, st) in ((60.0, :pass), (69.4, :pass), (69.6, :fail), (75.0, :fail))
        ri = talis_check(big, [Beam("i"; theta = θ, phi = 30, spot, fnumber = 8.0)]; capture = (:f,))
        inc = trow(ri, :incidence_3w, :f)
        @test inc.nominal ≈ θ atol = 1e-6
        @test inc.worst ≈ θ + 0.5 atol = 1e-6
        @test inc.margin ≈ 70 - θ - 0.5 atol = 1e-6
        @test inc.status === st && inc.unit === :deg
    end

    # 1ω footprint, analytic: one NIF beamlet aimed at TCC. A 0.1 mm square in the 3ω
    # focal plane with its far edge δ inside the 4.45 mm clear radius on the beamlet's own
    # side; the focus error moves the near edge by 0.027 mm per mm (cone to the 1ω focus).
    nb = nif_beams()
    bl = first(nb)
    k = direction(bl)
    pq = Briareus._nif_quad_port(bl.quad)
    ô = Briareus._own_side(bl.port, pq)
    @test dot(ô, k) ≈ 0 atol = 1e-12
    @test dot(bl.port - pq, ô) > 0                       # own side of the quad axis
    Rm = hcat(ô, cross(k, ô), k)
    δ = 0.5
    sq(name, c) = Component(name, transform(rectangle_mesh(0.1, 0.1; nu = 2, nv = 2), Rm, c), gold();
                            two_sided = true)
    tc = 4.45 - δ - 0.05
    t1 = Target(sq(:own, tc * ô), sq(:other, -tc * ô))
    r1 = talis_check(t1, [bl]; capture = (), open_ports = [pq])
    c0w, h0w = 4.45 + 22.7 / 2, 22.7 / 2
    taper = c0w / 7700 + h0w / 460
    o1 = trow(r1, :clearance_1w, :own)
    @test r1.settings.omega1_size ≈ 22.7
    @test o1.nominal ≈ δ atol = 1e-9
    @test o1.worst ≈ δ - taper - R - tan05 * hypot(4.45 - δ, 0.05) atol = 1e-9
    @test o1.status === :pass                              # 0.5 mm in: clears by > 0.2 mm
    @test trow(r1, :clearance_1w, :other).nominal > 8      # the other side is dark
    # prism model: no taper
    rp = talis_check(t1, [bl]; capture = (), omega1_focus = Inf, lens_distance = Inf)
    @test trow(rp, :clearance_1w, :own).worst ≈ δ - R - tan05 * hypot(4.45 - δ, 0.05) atol = 1e-9
    # the whole quad lights both sides
    q = only(b for b in nif_quads() if b.quad == bl.quad)
    rq = talis_check(t1, [q]; capture = (), open_ports = [pq])
    @test trow(rq, :clearance_1w, :other).nominal < 1
    # an explicit offset direction overrides the NIF convention
    ro = talis_check(t1, [bl]; capture = (), offsets = Dict(bl.name => -ô))
    @test trow(ro, :clearance_1w, :other).nominal ≈ δ atol = 1e-9

    # Uncleared flat plate facing the beamlet: the specular 1ω goes straight back along
    # the port, at the beamlet-to-quad-axis angle from its own quad's FOA.
    pl = Component(:mirror, transform(rectangle_mesh(4.0, 4.0; nu = 2, nv = 2), Rm, 8.0 * ô), gold())
    rm = talis_check(Target(pl), [bl]; capture = (), open_ports = [pq])
    cm = trow(rm, :clearance_1w, :mirror)
    @test cm.status === :review && cm.margin < 0
    mg = trow(rm, :management_1w, :mirror)
    α = acosd(dot(bl.port, pq))
    @test mg.unit === :deg
    @test mg.nominal ≈ α atol = 1e-6
    @test mg.worst ≈ α - 1.0 atol = 1e-6               # less twice the 0.5° tolerance
    @test mg.status === :fail
    # tilted 45° about ô⊥ the reflection leaves at 90° from the beam: passes
    Rt = Rm * [cosd(45) 0 sind(45); 0 1 0; -sind(45) 0 cosd(45)]
    pt = Component(:mirror, transform(rectangle_mesh(4.0, 4.0; nu = 2, nv = 2), Rt, 8.0 * ô), gold())
    mt = trow(talis_check(Target(pt), [bl]; capture = (), open_ports = [pq]), :management_1w, :mirror)
    @test mt.status === :pass && mt.worst > 60

    # Curved surface reflecting toward a port: the upstream pole of a sphere in the 1ω
    # footprint reflects 1ω straight back toward the FOA. R = 6 mm (≥ 5 mm) fails the
    # specular test; R = 4 mm passes on curvature; so does a < 2 mm cylinder.
    function ball(rad, c)
        Component(:ball, translate(sphere_mesh(rad; ntheta = 24, nphi = 48), c), gold();
                  surface = SphereSurface(rad; center = c))
    end
    rb6 = talis_check(Target(ball(6.0, 12.0 * ô)), [bl]; capture = (), open_ports = [-k])
    m6 = trow(rb6, :management_1w, :ball)
    @test m6.status === :fail && m6.worst < 5 - 1.0
    rb4 = talis_check(Target(ball(4.0, 10.0 * ô)), [bl]; capture = (), open_ports = [-k])
    m4 = trow(rb4, :management_1w, :ball)
    @test m4.status === :pass && m4.unit === :mm
    @test m4.margin ≈ 1.0
    cyl(rad) = Component(:stalk, transform(cylinder_mesh(rad, -3.0, 3.0; nz = 6, nphi = 24, inward = false),
                                           hcat(ô, k, cross(ô, k)), 9.0 * ô), gold();
                         surface = CylinderSurface(rad; origin = 9.0 * ô, axis = cross(ô, k)))
    @test trow(talis_check(Target(cyl(1.5)), [bl]; capture = (), open_ports = [-k]), :management_1w, :stalk).status === :pass
    ms = trow(talis_check(Target(cyl(2.5)), [bl]; capture = (), open_ports = [-k]), :management_1w, :stalk)
    @test ms.unit === :deg && ms.status === :fail          # its rim line reflects straight back
    # curvature overrides: dimpled, too-fine dimples, and undetermined meshes
    @test trow(talis_check(Target(pl), [bl]; capture = (), open_ports = [pq],
                           curvature = Dict(:mirror => (:dimpled, 1.0))), :management_1w, :mirror).status === :pass
    @test trow(talis_check(Target(pl), [bl]; capture = (), open_ports = [pq],
                           curvature = Dict(:mirror => (:dimpled, 1.0, 0.1))), :management_1w, :mirror).status === :review
    ub = Component(:ball, translate(sphere_mesh(6.0; ntheta = 24, nphi = 48), 12.0 * ô), gold())   # no surface tag
    @test trow(talis_check(Target(ub), [bl]; capture = (), open_ports = [-k]), :management_1w, :ball).status === :review

    # NIF ignition-style "sombrero": an annular 1ω shield in the LEH plane around a
    # hohlraum. Flat, it mirrors each quad into the conjugate FOA of its own hemisphere and
    # fails; flagged as spherically dimpled (radius < 5 mm) it passes.
    th = hohlraum(length = 10.1, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9,
                  nz = 12, nphi = 24, nr_cap = 2)
    som = Component(:sombrero, annulus_mesh(3.2, 14.0, 5.05; nr = 4, nphi = 48, up = true), gold())
    ts = Target(vcat(th.components, [som]))
    qs = nif_quads()
    rz = Dict(:inner23 => (1.09, 3.09), :inner30 => (0.6, 3.63), :outer44 => (0.0, 5.2), :outer50 => (0.0, 5.01))
    qs = map(filter(b -> b.port[3] > 0, qs)) do b
        r0, z0 = rz[b.group]
        φ = port_angles(b)[2]
        remake(b; aim = (-r0 * cosd(φ), -r0 * sind(φ), z0))
    end
    qs = [qs[findfirst(b -> b.group === g, qs)] for g in (:inner23, :inner30, :outer44, :outer50)]
    rs = talis_check(ts, qs; capture = (:wall, :capsule), level = 1)
    flat = [x for x in rs.rows if x.rule === :management_1w && x.component === :sombrero]
    @test !isempty(flat) && all(x -> x.status === :fail, flat)
    @test all(x -> occursin("flat", x.note), flat)
    rd = talis_check(ts, qs; capture = (:wall, :capsule), level = 1,
                     curvature = Dict(:sombrero => (:dimpled, 2.5)))
    dimp = [x for x in rd.rows if x.rule === :management_1w && x.component === :sombrero]
    @test length(dimp) == length(flat) && all(x -> x.status === :pass, dimp)
    # beams entering through the LEH: no capture rows fail on the wall, and the caps clear
    @test all(x -> x.status === :pass, (x for x in rs.rows if x.rule === :capture_3w))
    @test all(x -> x.margin > 0, (x for x in rs.rows if x.rule === :clearance_3w && x.component === :sombrero))
    @test any(x -> x.rule === :tas_stay_in, rs.rows)

    # Backlighter foil with OMEGA beams: capture, incidence, and a stalk behind the foil.
    ob = filter(b -> b.port[3] > 0, omega_beams(; cones = :all, spot = omega_spot(:SG4), fnumber = 6.7))
    tb = Target(foil(:bl, (0, 0, 0), (0, 0, 1); width = 4.0, n = 16),
                Component(:stalk, transform(cylinder_mesh(0.1, 0.0, 5.0; nz = 4, nphi = 12, inward = false),
                                            [1 0 0; 0 1 0; 0 0 1], Vec3(0, 0, -6.0)), gold()))
    rb = talis_check(tb, ob; capture = (:bl,))
    @test any(x -> occursin("NIF quad", x), rb.notes)                    # no 1ω rules off NIF
    @test !any(x -> x.rule in (:clearance_1w, :tas_stay_in), rb.rows)
    cap = [x for x in rb.rows if x.rule === :capture_3w]
    @test length(cap) == length(ob)
    θof(x) = port_angles(only(b for b in ob if b.name == x.beam))[1]
    @test all(x -> x.status === :pass, (x for x in cap if θof(x) < 45))    # cones 1, 2
    @test all(x -> x.status === :fail, (x for x in cap if θof(x) > 80))    # cone 4 grazes
    # the substrate edge seen orthogonal to the beam, along its tilt: at most the distance
    # to the square's edge in the beam's azimuth × cos θ, less the contour radius
    φof(x) = port_angles(only(b for b in ob if b.name == x.beam))[2]
    @test all(x -> x.nominal <= 2 / max(abs(cosd(φof(x))), abs(sind(φof(x)))) * cosd(θof(x)) -
                                0.36 * log(10)^(1 / 4.2) + 1e-6, cap)
    for x in rb.rows
        x.rule === :incidence_3w || continue
        θb = port_angles(only(b for b in ob if b.name == x.beam))[1]
        @test x.worst ≈ θb + 0.5 atol = 1e-6
        @test x.status === (θb + 0.5 < 70 ? :pass : :fail)
    end
    @test count(x -> x.rule === :incidence_3w && x.status === :fail, rb.rows) ==
          count(b -> port_angles(b)[1] > 69.5, ob) > 0
    st = [x for x in rb.rows if x.rule === :clearance_3w && x.component === :stalk]
    @test all(x -> x.status === :pass && x.worst == Inf, st)           # behind the foil
    # repointing one beam toward the edge fails its capture
    ob2 = copy(ob)
    ob2[1] = remake(ob2[1]; aim = (1.0, 0.0, 0.0))
    rb2 = talis_check(tb, ob2; capture = (:bl,))
    @test trow(rb2, :capture_3w, :bl; beam = ob2[1].name).status === :fail
    @test first(violations(rb2)).status === :fail

    # TAS stay-in zone: a part 20 mm above the lowest upper-beam aim point fails.
    up = translate(Component(:post, disk_mesh(0.5, 0.0; nr = 1, nphi = 8), gold()), (8.0, 8.0, 20.0))
    rt = talis_check(Target(foil(:f, (0, 0, 0), (0, 0, 1); width = 4.0, n = 4), up),
                     [remake(qs[1]; aim = (0.0, 0.0, 0.0))]; capture = (:f,), omega1 = false)
    tas = trow(rt, :tas_stay_in, :post)
    @test tas.status === :fail
    @test tas.margin ≈ -1.0 atol = 1e-9
    @test trow(rt, :tas_stay_in, :f).status === :pass

    # report output
    io = IOBuffer()
    show(io, MIME"text/plain"(), rs)
    txt = String(take!(io))
    @test occursin("ASSUMED", txt) && occursin("r·tan(0.5°)", txt) && occursin("management_1w", txt)
    path = write_talis_csv(tempname() * ".csv", rs)
    @test countlines(path) == 1 + length(rs.rows)
    @test length(TALIS_OUT_OF_SCOPE) >= 5
end
