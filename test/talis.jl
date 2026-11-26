# NIF TaLIS configuration checks (src/talis.jl). Included inside the top-level testset.

@testset "TaLIS" begin
    trow(r, rule, c; beam = nothing) =
        only(x for x in r.rows if x.rule === rule && x.component === c && (beam === nothing || x.beam == beam))
    tan05 = tand(0.5)
    R = 0.05 + 0.025                                   # default pointing + target offset
    # distances are lower bounds of the exact minimum, within the resolution (1e-4 mm)
    lowerbound(v, exact; ε = 1e-4) = exact - ε - 1e-9 <= v <= exact + 1e-9

    # 3ω clearance, analytic: beam straight down onto a foil at z = 0; a plate upstream at
    # z = 2 whose edge runs along x = x0. Worst case = x0 − (1% contour + f-cone at s = −2)
    # − transverse errors − focus/(2F) − r·tan 0.5° at the edge point nearest the axis.
    spot = SuperGaussianSpot(0.3; order = 4.0)
    b = Beam("b"; theta = 0, phi = 0, power = 1.0, spot, fnumber = 8.0)
    x0, w = 1.6, 2.0
    plate = Component(:plate, translate(rectangle_mesh(w, w; nu = 4, nv = 4), (x0 + w / 2, 0, 2.0)), gold();
                      two_sided = true)
    t = Target(foil(:f, (0, 0, 0), (0, 0, 1); width = 4.0, n = 8), plate)
    r = talis_check(t, [b]; capture = (:f,), rtol = 0)
    ρ1 = log(100)^(1 / 4) * 0.3                         # 1%-of-peak contour of exp(-ρ⁴)
    nominal = x0 - (ρ1 + 2 / 16)
    worst = nominal - R - 1 / 16 - tan05 * hypot(x0, 2.0)
    c3 = trow(r, :clearance_3w, :plate)
    @test lowerbound(c3.nominal, nominal)
    @test lowerbound(c3.worst, worst)
    @test c3.margin == c3.worst - 0.5
    @test c3.point ≈ Vec3(x0, 0, 2.0)
    @test c3.status === :pass && c3.unit === :mm
    # Dropping the tolerance term, or moving the pivot, changes only r·tan(0.5°).
    c0 = trow(talis_check(t, [b]; capture = (:f,), tolerance = 0, rtol = 0), :clearance_3w, :plate)
    @test c0.worst - c3.worst ≈ tan05 * hypot(x0, 2.0) atol = 1e-4
    bA = remake(b; aim = (0.0, 0.0, 0.5))               # aim 0.5 mm up the axis
    cA = trow(talis_check(t, [bA]; capture = (:f,), pivot = :aim, rtol = 0), :clearance_3w, :plate)
    @test lowerbound(cA.worst, x0 - (ρ1 + 1.5 / 16) - R - 1 / 16 - tan05 * hypot(x0, 1.5))
    # With the pivot on the edge point, the allowance at that point is zero. Along the edge
    # r·tan(0.5°) grows linearly and the distance only quadratically, so the minimum lies
    # off the pivot (|y| ≈ 0.5° · x0 = 14 µm) and below the vertex value. The reference is
    # the exact edge distance on a 0.1 µm grid.
    cP = trow(talis_check(t, [b]; capture = (:f,), pivot = (x0, 0.0, 2.0), rtol = 0), :clearance_3w, :plate)
    refP = minimum(hypot(x0, y) - (ρ1 + 2 / 16) - R - 1 / 16 - tan05 * abs(y) for y in 0:1e-7:0.1)
    @test refP < nominal - R - 1 / 16 - 5e-5 && lowerbound(cP.worst, refP)
    cM = trow(talis_check(t, [bA]; capture = (:f,), pivot = :max, rtol = 0), :clearance_3w, :plate)
    @test cM.worst <= min(cA.worst, trow(talis_check(t, [bA]; capture = (:f,), rtol = 0), :clearance_3w, :plate).worst) + 1e-4
    # The beam of one quad crosses a plate between its vertices, with the corners ≥ 50 mm
    # off the axis. The minimum over the whole triangle lies on the axis; the vertices do
    # not set it. The exact value is −ρ1 − 2/16 (inside the envelope), less the error stack
    # and 2·tan(0.5°). At the default rtol the bound lies within 1% of its distance from the
    # 0.5 mm rule.
    wide = Component(:plate, translate(rectangle_mesh(100.0, 100.0; nu = 1, nv = 1), (8.0, 8.0, 2.0)), gold();
                     two_sided = true)
    cw = trow(talis_check(Target(foil(:f, (0, 0, 0), (0, 0, 1); width = 4.0, n = 8), wide), [b];
                          capture = (:f,), level = 0), :clearance_3w, :plate)
    exw = -ρ1 - 2 / 16 - R - 1 / 16 - 2tan05
    @test cw.status === :fail && lowerbound(cw.worst, exw; ε = 0.01 * (0.5 - exw))
    # moving the plate in fails the rule, with the blow-by consequence
    tin = Target(foil(:f, (0, 0, 0), (0, 0, 1); width = 4.0, n = 8),
                 translate(plate, (-0.5, 0.0, 0.0)))
    cf = trow(talis_check(tin, [b]; capture = (:f,)), :clearance_3w, :plate)
    @test cf.status === :fail && cf.margin < 0 && occursin("blow-by", cf.note)

    # 3ω capture, analytic: the 10% contour sits 2 − xa − ρ10 from the foil edge at x = 2.
    ρ10 = log(10)^(1 / 4) * 0.3
    for xa in (0.0, 0.6, 1.0)
        rc = talis_check(t, [remake(b; aim = (xa, 0.0, 0.0))]; capture = (:f,), rtol = 0)
        cp = trow(rc, :capture_3w, :f)
        @test lowerbound(cp.nominal, 2 - xa - ρ10)
        @test lowerbound(cp.worst, 2 - xa - ρ10 - R - 1 / 16 - tan05 * 2)
        @test cp.status === (cp.worst >= 0.5 ? :pass : :fail)
        xa > 0 && @test cp.point ≈ Vec3(2, 0, 0) atol = 1e-9
    end
    @test trow(talis_check(t, [remake(b; aim = (1.0, 0.0, 0.0))]; capture = (:f,)), :capture_3w, :f).status === :fail
    @test trow(talis_check(t, [remake(b; aim = (0.9, 0.0, 0.0))]; capture = (:f,)), :capture_3w, :f).status === :pass
    # an axis that misses the substrate fails capture outright
    rmiss = talis_check(t, [remake(b; aim = (3.0, 0.0, 0.0))]; capture = (:f,))
    @test trow(rmiss, :capture_3w, :none).status === :fail
    # For incidence on a flat foil with the beam axis at θ from its normal, the worst angle
    # is θ plus the half-angle of the f-cone, atan(1/16) = 3.576°, with no rule tolerance.
    # The limit therefore falls at θ = 66.42°.
    big = Target(foil(:f, (0, 0, 0), (0, 0, 1); width = 12.0, n = 12))
    hc = atand(1 / 16)
    incs = [(θ, st, trow(talis_check(big, [Beam("i"; theta = θ, phi = 30, spot, fnumber = 8.0)]; capture = (:f,)),
                         :incidence_3w, :f)) for (θ, st) in ((60.0, :pass), (66.3, :pass), (66.5, :fail), (75.0, :fail))]
    @test all(isapprox(i.nominal, θ + hc; atol = 1e-6) && isapprox(i.worst, θ + hc; atol = 1e-6) &&
              isapprox(i.margin, 70 - θ - hc; atol = 1e-6) && i.status === st && i.unit === :deg
              for (θ, st, i) in incs)

    # 1ω footprint, analytic: one NIF beamlet aimed at TCC. A 0.1 mm square in the 3ω
    # focal plane with its far edge δ inside the 4.45 mm clear radius on the side
    # of the beamlet; the focus error moves the near edge by 0.027 mm per mm (cone to the 1ω focus).
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
    r1 = talis_check(t1, [bl]; capture = (), open_ports = [pq], rtol = 0)
    c0w, h0w = 4.45 + 22.7 / 2, 22.7 / 2
    taper = c0w / 7700 + h0w / 460
    o1 = trow(r1, :clearance_1w, :own)
    @test r1.settings.omega1_size ≈ 22.7
    @test lowerbound(o1.nominal, δ)
    @test lowerbound(o1.worst, δ - taper - R - tan05 * hypot(4.45 - δ, 0.05))
    @test o1.status === :pass                              # 0.5 mm in: clears by > 0.2 mm
    @test trow(r1, :clearance_1w, :other).nominal > 8      # the other side is dark
    # The prism model has no taper.
    rp = talis_check(t1, [bl]; capture = (), omega1_focus = Inf, lens_distance = Inf, rtol = 0)
    @test lowerbound(trow(rp, :clearance_1w, :own).worst, δ - R - tan05 * hypot(4.45 - δ, 0.05))
    # a 100 mm plate in the focal plane over the whole footprint, corners ≥ 20 mm out of
    # it: the footprint is found between the vertices; deepest point −h0 at its center
    # (with the focus error, the tapered square of the −1 mm plane)
    wp = Component(:wide, transform(rectangle_mesh(100.0, 100.0; nu = 1, nv = 1), Rm, 50.0 * ô), gold();
                   two_sided = true)
    rw = talis_check(Target(wp), [bl]; capture = (), open_ports = [pq], level = 0, rtol = 0)
    ow = trow(rw, :clearance_1w, :wide)
    @test ow.status === :review && lowerbound(ow.nominal, -h0w)
    # no `diagnostics`: their 1ω rule is reported as not evaluated; given, it is checked
    @test any(n -> occursin("diagnostics to 1ω", n), rw.notes)
    rdg = talis_check(Target(wp), [bl]; capture = (), open_ports = [pq], level = 0, diagnostics = (:wide,))
    @test trow(rdg, :diagnostic_1w, :wide).status === :fail && !any(n -> occursin("diagnostics to 1ω", n), rdg.notes)
    @test rdg.settings.diagnostics == [:wide]
    # the whole quad lights both sides
    q = only(b for b in nif_quads() if b.quad == bl.quad)
    rq = talis_check(t1, [q]; capture = (), open_ports = [pq])
    @test trow(rq, :clearance_1w, :other).nominal < 1
    # an explicit offset direction overrides the NIF convention
    ro = talis_check(t1, [bl]; capture = (), offsets = Dict(bl.name => -ô), rtol = 0)
    @test lowerbound(trow(ro, :clearance_1w, :other).nominal, δ)

    # Uncleared flat plate facing the beamlet: the specular 1ω goes straight back along
    # the port, at the beamlet-to-quad-axis angle from the FOA of its quad.
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
    rt = talis_check(Target(pt), [bl]; capture = (), open_ports = [pq])
    mt = trow(rt, :management_1w, :mirror)
    @test mt.status === :pass && mt.worst > 60
    # the reflection behind each row: the ray, the port it threatens and the projected
    # areas of the plate, which lies whole in the footprint
    xm, xt = omega1_reflection(rm, mg), omega1_reflection(rt, mt)
    @test xm.beam == bl.name && xm.component === :mirror && xm.beamlet_index == 1 && xm.port_index == 1
    @test xm.incident ≈ k && xm.reflected ≈ -k && xm.port == pq && abs(dot(xm.normal, k)) ≈ 1
    @test xm.angle ≈ α atol = 1e-6
    @test xm.spread == 0 && xm.worst == mg.worst && xm.margin == mg.margin && xm.point == mg.point
    @test xm.exposed_area ≈ 16 rtol = 1e-12
    @test xm.offending_area ≈ 16 rtol = 1e-12
    @test xt.exposed_area ≈ 16 * cosd(45) rtol = 1e-12
    @test xt.offending_area == 0 && abs(dot(xt.reflected, k)) < 1e-12
    @test occursin("from 16 mm² of 16 mm²", mg.note)
    # The reflected 1ω power, against the uniform 22 mm square of the beamlet: the plate
    # takes 16/22² of it in the focal plane; 100 mm upstream the square is wider by the
    # cone to the focus 460 mm beyond. A plate over the whole square, with mesh lines on
    # its edges, reflects all of it; the tilted plate reflects none toward the port.
    @test xm.reflected_fraction ≈ 16 / 22^2 rtol = 1e-12
    @test xt.reflected_fraction == 0
    @test omega1_reflected_power(xm, 0.5; reflectivity = 0.9) ≈ 0.45 * 16 / 22^2 rtol = 1e-12
    @test_throws ArgumentError omega1_reflected_power(xm, 0.5; reflectivity = 1.5)
    mirror(w, h, nu, nv, c) = Target(Component(:mirror, transform(rectangle_mesh(w, h; nu, nv), Rm, c), gold()))
    xu = only(talis_check(mirror(4.0, 4.0, 2, 2, 8.0 * ô - 100.0 * k), [bl]; capture = (), open_ports = [pq]).reflections)
    @test xu.reflected_fraction ≈ 16 / (22 * (1 + 100 / 460))^2 rtol = 1e-12
    xw = only(talis_check(mirror(30.0, 26.0, 15, 13, 15.8 * ô), [bl]; capture = (), open_ports = [pq]).reflections)
    @test xw.reflected_fraction ≈ 1 rtol = 1e-12
    @test omega1_reflection(rm, cm) === nothing && length(rm.reflections) == 1

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
    @test omega1_reflection(rb4, m4) === nothing && omega1_reflection(rb6, m6).spread > 0
    @test m4.status === :pass && m4.unit === :mm
    @test m4.margin ≈ 1.0
    # A coarse mesh of that sphere and a port on the reflection from one of its vertices:
    # the reflections of the sampled normals all miss the port by more than 21°, and the
    # spread of the normal over each facet brings the port within reach.
    cb = 12.0 * ô
    coarse = Component(:ball, translate(sphere_mesh(6.0; ntheta = 4, nphi = 8), cb), gold();
                       surface = SphereSurface(6.0; center = cb))
    tcb = Target(coarse)
    vs = unique(v for j in eachindex(tcb.tri_area) for v in (tcb.tri_v0[j], tcb.tri_v0[j] + tcb.tri_e1[j]))
    v = vs[argmax([dot(normalize(x - cb), -k) for x in vs])]
    nvx = normalize(v - cb)
    port = k - 2 * dot(k, nvx) * nvx
    sampled = [(n = normalize(tcb.tri_v0[j] + (tcb.tri_e1[j] + tcb.tri_e2[j]) / 3 - cb); k - 2 * dot(k, n) * n)
               for j in eachindex(tcb.tri_area)]
    @test minimum(acosd(clamp(dot(r, port), -1, 1)) for r in sampled) > 21
    mcb = trow(talis_check(tcb, [bl]; capture = (), open_ports = [port], level = 0), :management_1w, :ball)
    @test mcb.status === :fail && mcb.worst <= -1 + 1e-9
    cyl(rad) = Component(:stalk, transform(cylinder_mesh(rad, -3.0, 3.0; nz = 6, nphi = 24, inward = false),
                                           hcat(ô, k, cross(ô, k)), 9.0 * ô), gold();
                         surface = CylinderSurface(rad; origin = 9.0 * ô, axis = cross(ô, k)))
    @test trow(talis_check(Target(cyl(1.5)), [bl]; capture = (), open_ports = [-k]), :management_1w, :stalk).status === :pass
    ms = trow(talis_check(Target(cyl(2.5)), [bl]; capture = (), open_ports = [-k]), :management_1w, :stalk)
    @test ms.unit === :deg && ms.status === :fail          # its rim line reflects straight back
    # A tapered stalk across the footprint, its axis orthogonal to the beam: the radius of
    # curvature of a cone is its distance from the axis over cos α. The ruling that faces
    # the beam reflects 1ω at 2α from the port. A stalk that stays below 2 mm passes on
    # curvature, on its largest radius, at the base. A wider one reflects into a port
    # placed on that ray from the part of it at or above 2 mm.
    ac = cross(ô, k)
    function taper(r0, r1)
        half = atand((r1 - r0) / 6)
        Component(:stalk, transform(frustum_mesh(r0, -3.0, r1, 3.0; nz = 6, nphi = 24, inward = false),
                                    hcat(ô, k, ac), 9.0 * ô), gold();
                  surface = ConeSurface(half; apex = 9.0 * ô - (3 + r0 / tand(half)) * ac, axis = ac))
    end
    spec(half) = -cosd(2half) * k - sind(2half) * ac
    hs, hw = atand(1 / 6), atand(2 / 6)
    mc = trow(talis_check(Target(taper(0.5, 1.5)), [bl]; capture = (), open_ports = [spec(hs)]), :management_1w, :stalk)
    @test mc.status === :pass && mc.unit === :mm && mc.required == 2
    @test mc.worst ≈ 1.5 / cosd(hs) rtol = 1e-12
    @test mc.margin ≈ 2 - 1.5 / cosd(hs) rtol = 1e-12
    mw = trow(talis_check(Target(taper(1.0, 3.0)), [bl]; capture = (), open_ports = [spec(hw)]), :management_1w, :stalk)
    @test mw.status === :fail && mw.unit === :deg && mw.nominal < 3
    @test dot(mw.point - 9.0 * ô, ac) > -0.5              # found on the part with radius ≥ 2 mm
    # A torus: the radii are the minor radius and |R/cos ψ + r|, the largest over a
    # triangle, unbounded where it crosses a pole of the tube.
    tor = TorusSurface(4.0, 3.0)
    tp(ψ, φ) = Vec3((4 + 3cosd(ψ)) * cosd(φ), (4 + 3cosd(ψ)) * sind(φ), 3sind(ψ))
    @test collect(Briareus._surface_radii(tor, (tp(0, 0), tp(20, 0), tp(20, 15)))) ≈ [3, 4 / cosd(20) + 3]
    @test collect(Briareus._surface_radii(tor, (tp(180, 0), tp(150, 0), tp(150, 15)))) ≈ [4 / cosd(30) - 3, 3]
    @test Briareus._surface_radii(tor, (tp(80, 0), tp(100, 0), tp(100, 15))) == (3.0, Inf)
    @test collect(Briareus._surface_radii(ConeSurface(30.0), (Vec3(1, 0, 2), Vec3(0, 2, 4), Vec3(1, 1, 3)))) ≈
          [4 * tand(30) / cosd(30), Inf]
    function ring(Rj, rn)
        m = Mesh()
        P(ψ, φ) = Vec3((Rj + rn * cosd(ψ)) * cosd(φ), (Rj + rn * cosd(ψ)) * sind(φ), rn * sind(ψ))
        for ψ in 0:30:330, φ in 0:15:345
            a, b2, c2, d = P(ψ, φ), P(ψ, φ + 15), P(ψ + 30, φ + 15), P(ψ + 30, φ)
            out = P(ψ + 15, φ + 7.5) - Vec3(Rj * cosd(φ + 7.5), Rj * sind(φ + 7.5), 0)
            add_polygon!(m, [a, b2, c2]; facing = out)
            add_polygon!(m, [a, c2, d]; facing = out)
        end
        Component(:ring, transform(m, Rm, 14.0 * ô), gold(); surface = TorusSurface(Rj, rn; center = 14.0 * ô, axis = k))
    end
    mr = trow(talis_check(Target(ring(6.0, 1.5)), [bl]; capture = (), open_ports = [-k]), :management_1w, :ring)
    @test mr.status === :pass && mr.unit === :mm && mr.margin ≈ 0.5
    # tube radius 3 mm: the outer half is gentle and its crown reflects straight back
    mg3 = trow(talis_check(Target(ring(4.0, 3.0)), [bl]; capture = (), open_ports = [-k]), :management_1w, :ring)
    @test mg3.status === :fail && mg3.unit === :deg

    # The curvature overrides cover dimpled meshes, too-fine dimples, dimples of unknown
    # spacing (User Guide §5.6.2 needs half-side/radius ≥ 0.25), and undetermined meshes.
    dimpled(c) = trow(talis_check(Target(pl), [bl]; capture = (), open_ports = [pq],
                                  curvature = Dict(:mirror => c)), :management_1w, :mirror)
    @test dimpled((:dimpled, 1.0, 0.25)).status === :pass
    @test dimpled((:dimpled, 1.0, 0.1)).status === :review
    du = dimpled((:dimpled, 1.0))
    @test du.status === :review && occursin("spacing unknown", du.note)
    ub = Component(:ball, translate(sphere_mesh(6.0; ntheta = 24, nphi = 48), 12.0 * ô), gold())   # no surface tag
    @test trow(talis_check(Target(ub), [bl]; capture = (), open_ports = [-k]), :management_1w, :ball).status === :review
    # Meshes without a surface have their radii estimated from the mesh, exact on a
    # regular cylinder. A bare 1.5 mm stalk scatters on its estimate with a quarter added;
    # a bare 2.5 mm one is given the reflection check and reviewed.
    bare(m) = Component(:stalk, transform(m, hcat(ô, k, cross(ô, k)), 9.0 * ô), gold())
    rod(rad) = cylinder_mesh(rad, -3.0, 3.0; nz = 6, nphi = 24, inward = false)
    mb = trow(talis_check(Target(bare(rod(1.5))), [bl]; capture = (), open_ports = [-k]), :management_1w, :stalk)
    @test mb.status === :pass && mb.unit === :mm
    @test mb.worst ≈ 1.25 * 1.5 rtol = 1e-9
    @test mb.margin ≈ 2 - 1.25 * 1.5 rtol = 1e-9
    @test trow(talis_check(Target(bare(rod(2.5))), [bl]; capture = (), open_ports = [-k]), :management_1w, :stalk).status === :review
    # A stalk and a flat flange in one mesh: the stalk scatters, and the flange, edge-on
    # to the beam, reflects it forward, far from the port behind. Judged as one unknown
    # surface, the ruling of the stalk that faces the beam sent the item to review.
    solid = bare(merge_meshes(rod(1.5), annulus_mesh(1.5, 4.0, -3.0; nr = 3, nphi = 24)))
    rs = talis_check(Target(solid), [bl]; capture = (), open_ports = [-k])
    msol = trow(rs, :management_1w, :stalk)
    @test msol.status === :pass && msol.unit === :deg && msol.worst > 150
    @test omega1_reflection(rs, msol).offending_area == 0
    # the estimate itself: cylinders to rounding, spheres within a fifth, a plate flat
    function estimate(m)
        tm = Target(Component(:c, m, gold()))
        tris = [(tm.tri_v0[j], tm.tri_v0[j] + tm.tri_e1[j], tm.tri_v0[j] + tm.tri_e2[j]) for j in eachindex(tm.tri_area)]
        Briareus._mesh_curvature(tris, tm.tri_normal)
    end
    ec = estimate(cylinder_mesh(2.5, -3.0, 3.0; nz = 3, nphi = 12, inward = false))
    @test all(abs(x[1] - 2.5) < 1e-9 && x[2] > 1e9 && x[3] ≈ 15 && x[4] < 1e-9 for x in ec)
    es = estimate(icosphere_mesh(6.0; n = 4))
    @test all(0.8 * 6 < x[1] <= x[2] < 1.2 * 6 && x[4] < 0.1 for x in es)
    @test all(x -> x[1] == Inf && x[3] == 0, estimate(rectangle_mesh(4.0, 4.0; nu = 3, nv = 3)))
    # A check fits only the triangles that it reaches, with the fits of the full estimate,
    # and beams on several threads share them.
    let tl = Target(Component(:c, icosphere_mesh(6.0; n = 4), gold()))
        cxl = (; curv = Dict{Int,Briareus._LazyCurvature}(), curv_lock = ReentrantLock(),
               geo = Dict(1 => Briareus._talis_samples(tl, 1, 0)))
        lc = Briareus._component_curvature(tl, 1, cxl)
        full = Briareus._mesh_curvature(cxl.geo[1].tris, tl.tri_normal; scale = tl.scale)
        ks = 1:3:length(full)
        got = Vector{NTuple{4,Float64}}(undef, length(ks))
        Threads.@threads for i in eachindex(ks)
            got[i] = lc[ks[i]]
        end
        @test got == full[ks]
        @test count(k -> (@atomic lc.state[k]) == 0x02, eachindex(full)) == length(ks)
        @test Briareus._component_curvature(tl, 1, cxl) === lc
    end

    # NIF ignition-style "sombrero": an annular 1ω shield in the LEH plane around a
    # hohlraum. Flat, it mirrors each quad into the conjugate FOA of its hemisphere and
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
                     curvature = Dict(:sombrero => (:dimpled, 2.5, 1.0)))
    dimp = [x for x in rd.rows if x.rule === :management_1w && x.component === :sombrero]
    @test length(dimp) == length(flat) && all(x -> x.status === :pass, dimp)
    # For beams entering through the LEH, every capture row on the wall passes, and the
    # caps clear.
    @test all(x -> x.status === :pass, (x for x in rs.rows if x.rule === :capture_3w))
    @test all(x -> x.margin > 0, (x for x in rs.rows if x.rule === :clearance_3w && x.component === :sombrero))
    @test any(x -> x.rule === :tas_stay_in, rs.rows)

    # LEH edges (`leh`): the 23.5° quads aimed as in the GUI's NIF template cannot clear
    # a 3.1 mm LEH by 0.5 mm; on the caps that is a review, on hardware a failure.
    # Independent reference: the best-focus 1% ellipse (HWHM table of NIF_SPOTS, 1/e
    # radius a·log(100)^(1/n)), grown by the f/7.9 cone, against the exact LEH circle,
    # less the linear error stack and r·tan(0.5°) (r = |rim point| from TCC).
    hl = hohlraum(length = 10.1, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9,
                  nz = 8, nphi = 48, nr_cap = 2)
    q23 = map(filter(b -> b.group === :inner23 && b.port[3] > 0, nif_quads())) do b
        φ = port_angles(b)[2]
        remake(b; aim = (-1.09 * cosd(φ), -1.09 * sind(φ), 3.09))
    end
    n23 = 7.37167
    a23, b23 = (0.88585, 0.65579) ./ log(2)^(1 / n23) .* log(100)^(1 / n23)
    @test a23 ≈ 1.145 atol = 5e-4                       # the 1% semi-axis along φ̂
    rim = [Vec3(1.55cos(ψ), 1.55sin(ψ), 5.05) for ψ in range(0, 2π; length = 2881)[1:end-1]]
    ref = minimum(q23) do b
        θ, φ = port_angles(b)
        eθ = Vec3(cosd(θ)cosd(φ), cosd(θ)sind(φ), -sind(θ))
        eφ = Vec3(-sind(φ), cosd(φ), 0.0)
        k = direction(b)
        minimum(rim) do x
            d = x - Vec3(b.aim)
            u, v, s = dot(d, eφ), dot(d, eθ), dot(d, k)
            minimum(hypot(u - a23 * cos(τ), v - b23 * sin(τ)) for τ in range(0, 2π; length = 2001)) -
                abs(s) / 15.8 - R - 1 / 15.8 - tan05 * norm(x)
        end
    end
    rl = talis_check(hl, q23; capture = (:wall,), hardware = (:cap_top,), leh = (:cap_top,),
                     omega1 = false, rtol = 0)
    lips = [x for x in rl.rows if x.rule === :clearance_3w]
    @test length(lips) == 4 && all(x -> x.component === :cap_top && x.status === :review, lips)
    @test all(x -> occursin("LEH clearance", x.note), lips) && rl.settings.leh == [:cap_top]
    # chords of the 48-gon lie up to 1.55(1 − cos 3.75°) = 3.3 µm inside the circle
    @test minimum(x -> x.worst, lips) ≈ ref atol = 5e-3
    @test minimum(x -> x.margin, lips) ≈ ref - 0.5 atol = 5e-3
    @test ref - 0.5 ≈ -0.4774 atol = 1e-3
    # the same cap given as hardware only, and a shield 50 µm outside the lip, fail
    rh = talis_check(hl, q23; capture = (:wall,), hardware = (:cap_top,), omega1 = false)
    @test all(x -> x.status === :fail, (x for x in rh.rows if x.rule === :clearance_3w))
    shield = Component(:shield, annulus_mesh(1.6, 4.0, 5.1; nr = 2, nphi = 48, up = true), gold())
    rsh = talis_check(Target(vcat(hl.components, [shield])), q23; capture = (:wall,),
                      hardware = (:shield,), leh = (:cap_top, :cap_bottom), omega1 = false)
    @test all(x -> x.status === :fail, (x for x in rsh.rows if x.rule === :clearance_3w && x.component === :shield))
    @test all(x -> x.status === :review, (x for x in rsh.rows if x.rule === :clearance_3w && x.component === :cap_top))

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
    θport(x) = port_angles(only(b for b in ob if b.name == x.beam))[1]
    inc = [x for x in rb.rows if x.rule === :incidence_3w]
    hcone(x) = atand(1 / (2only(b for b in ob if b.name == x.beam).fnumber))
    @test all(x -> isapprox(x.worst, θport(x) + hcone(x); atol = 1e-6) &&
                   x.status === (θport(x) + hcone(x) < 70 ? :pass : :fail), inc)
    @test count(x -> x.status === :fail, inc) == count(b -> port_angles(b)[1] + atand(1 / (2b.fnumber)) >= 70, ob)
    st = [x for x in rb.rows if x.rule === :clearance_3w && x.component === :stalk]
    @test all(x -> x.status === :pass && x.worst == Inf, st)           # behind the foil
    # repointing one beam toward the edge fails its capture
    ob2 = copy(ob)
    ob2[1] = remake(ob2[1]; aim = (1.0, 0.0, 0.0))
    rb2 = talis_check(tb, ob2; capture = (:bl,))
    @test trow(rb2, :capture_3w, :bl; beam = ob2[1].name).status === :fail
    @test first(violations(rb2)).status === :fail

    # In the TAS stay-in zone, a part 20 mm above the lowest upper-beam aim point fails.
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
    @test occursin("70° incidence limit excluded", txt) && occursin("r·tan(0.5°)", txt) && occursin("management_1w", txt)
    path = write_talis_csv(tempname() * ".csv", rs)
    lines = readlines(path)
    @test count(!startswith("#"), lines) == 1 + length(rs.rows)
    # The header states the scope: source, settings, skipped checks and out-of-scope items.
    @test occursin(TALIS_SOURCE, lines[1]) && "# setting tolerance = 0.5" in lines &&
          "# setting pivot = tcc" in lines && "# setting resolution = 0.0001" in lines
    @test count(startswith("# out of scope: "), lines) == length(TALIS_OUT_OF_SCOPE)
    @test count(startswith("# not evaluated: "), lines) == length(rs.notes)
    @test length(TALIS_OUT_OF_SCOPE) >= 5
end

@testset "TaLIS requirement text" begin
    # The requirement of each rule takes its thresholds from the report settings.
    s = (; capture_contour = 0.1, capture_margin = 0.2, max_incidence = 70.0, clearance_contour = 0.01,
         clearance_margin = 0.5, omega1_edge = 1.25, omega1_margin = 3.0, cylinder_radius = 4.0,
         sphere_radius = 2.5, reflection_angle = 10.0, stay_in = 6.0)
    req(rule) = talis_requirement(rule, s)
    @test req(:capture_3w) == "10% contour ≥ 0.2 mm inside the substrate edge"
    @test req(:incidence_3w) == "incidence < 70° to the surface normal"
    @test req(:clearance_3w) == "1% contour > 0.5 mm from the component"
    @test req(:clearance_1w) == "1ω 10% footprint (1.25 mm off the 3ω axis) > 3 mm from the component"
    @test req(:management_1w) ==
          "curvature radius < 4 mm (cylinder) or < 2.5 mm (sphere), or specular 1ω > 10° from open ports"
    @test req(:diagnostic_1w) == "diagnostics ≥ 10 mm (at TCC) to 25 mm (at 1 m) from 1ω light"
    @test req(:tas_stay_in) == "≤ 6 mm beyond the extreme aim points (TAS)"
    @test req(:other) == "other"
    # a report prints the requirement of each rule it checked, from its settings
    b = Beam("b"; theta = 0, phi = 0, power = 1.0, spot = SuperGaussianSpot(0.3; order = 4.0), fnumber = 8.0)
    t = Target(foil(:f, (0, 0, 0), (0, 0, 1); width = 4.0, n = 8))
    r = talis_check(t, [b]; capture = (:f,), rtol = 0)
    shown = sprint(show, MIME("text/plain"), r)
    @test all(occursin(talis_requirement(x.rule, r.settings), shown) for x in r.rows)
end

@testset "1ω closest approach" begin
    # the shared contour
    n, c10, cw = omega1_contour(mode = :nominal), omega1_contour(), omega1_contour(mode = :talis_worst)
    @test (n.edge, n.size, n.focal_shift, n.c0, n.h0) == (4.8, 22.0, 460.0, 15.8, 11.0)
    @test c10.edge == 4.45 && c10.size ≈ 22.7 && c10.c0 ≈ 15.8 && c10.h0 ≈ 11.35
    @test cw.edge ≈ 4.175 && cw.size ≈ 23.25 && cw.c0 == c10.c0 && cw.h0 ≈ c10.h0 + 0.275
    @test omega1_contour(edge = 1.25, size = 6.0, focus = Inf) ==
          (edge = 1.25, size = 6.0, focal_shift = Inf, c0 = 4.25, h0 = 3.0)
    @test_throws ArgumentError omega1_contour(mode = :other)

    # One beamlet and the 0.1 mm square of the analytic 1ω case above, its far edge δ
    # inside the near edge of the footprint. The near edge moves toward the 3ω axis by
    # `taper` per mm upstream of the focus, so the worst focal plane lies 1 mm downstream.
    bl = first(nif_beams())
    k = direction(bl)
    pq = Briareus._nif_quad_port(bl.quad)
    ô = Briareus._own_side(bl.port, pq)
    q̂ = cross(k, ô)
    Rm = hcat(ô, q̂, k)
    δ, tan05, R = 0.5, tand(0.5), 0.075
    sq(name, c) = Component(name, transform(rectangle_mesh(0.1, 0.1; nu = 2, nv = 2), Rm, c), gold();
                            two_sided = true)
    t1 = Target(sq(:own, (4.45 - δ - 0.05) * ô), sq(:other, -(4.45 - δ - 0.05) * ô))
    r1 = talis_check(t1, [bl]; capture = (), open_ports = [pq], rtol = 0)
    row = only(x for x in r1.rows if x.rule === :clearance_1w && x.component === :own)
    a = omega1_closest_approach(t1, [bl], r1, row)
    taper = c10.c0 / 7700 + c10.h0 / 460
    rot = tan05 * hypot(4.45 - δ, 0.05)
    @test a.beam == bl.name && a.component === :own && a.beamlet_index == 1 && a.x_hw == row.point
    @test a.beamlet[1] ≈ k && a.beamlet[3] ≈ ô && norm(a.beamlet[2]) < 1e-12
    @test a.focus_shift == 1.0 && a.face_mask == 0x01
    @test a.s_star ≈ -1 atol = 1e-12
    @test a.nominal ≈ δ atol = 1e-12
    @test a.distance ≈ δ - taper atol = 1e-12
    @test a.x_beam ≈ a.x_hw + (δ - taper) * ô atol = 1e-12
    @test abs(dot(a.x_hw - a.x_beam, k)) < 1e-12 && abs(norm(a.x_hw - a.x_beam) - a.distance) < 1e-12
    @test a.worst ≈ δ - taper - R - rot atol = 1e-12
    @test row.worst - 1e-9 <= a.worst <= row.worst + 1e-4 && a.nominal >= row.nominal - 1e-9
    b = a.budget
    @test (b.pointing, b.target_offset, b.hardware_tolerance, b.required) == (0.05, 0.025, 0.0, 0.2)
    @test b.rotation_tol ≈ rot atol = 1e-12
    @test b.focus_delta ≈ taper atol = 1e-12
    @test a.margin ≈ a.nominal - b.pointing - b.target_offset - b.hardware_tolerance - b.rotation_tol -
                     b.focus_delta - b.required atol = 1e-12
    @test a.margin == a.worst - a.required && a.contour.c0 == c10.c0 && a.contour.lens_distance == 7700
    # the part form finds the same point; the row form without the report takes the defaults
    ap = omega1_closest_approach(t1, [bl], bl.name => :own; rtol = 0)
    @test ap.component === :own && abs(ap.worst - a.worst) < 1e-4
    @test omega1_closest_approach(t1, [bl], row).worst == a.worst
    @test omega1_closest_approach(t1, [bl], 1 => Tuple(row.point)).x_beam == a.x_beam
    # another pivot and error set, through the settings of the report
    r2 = talis_check(t1, [bl]; capture = (), open_ports = [pq], rtol = 0, pivot = (0, 0, 12.0),
                     hardware_tolerance = 0.03, focus = 0.0)
    row2 = only(x for x in r2.rows if x.rule === :clearance_1w && x.component === :own)
    a2 = omega1_closest_approach(t1, [bl], r2, row2)
    @test a2.focus_shift == 0 && a2.budget.focus_delta == 0 && a2.budget.hardware_tolerance == 0.03
    @test a2.budget.rotation_tol ≈ tan05 * norm(row2.point - Vec3(0, 0, 12.0)) atol = 1e-12
    @test row2.worst - 1e-9 <= a2.worst <= row2.worst + 1e-4
    # a diagnostic row carries its required clearance
    rdg = talis_check(t1, [bl]; capture = (), open_ports = [pq], diagnostics = (:other,))
    rowd = only(x for x in rdg.rows if x.rule === :diagnostic_1w)
    ad = omega1_closest_approach(t1, [bl], rdg, rowd)
    @test ad.required == rowd.required && ad.component === :other
    @test rowd.margin - 1e-9 <= ad.margin <= rowd.margin + 0.01 * abs(rowd.margin) + 1e-4

    # Reference at points around a quad: the footprint square of each beamlet and focal
    # plane as four segments, and the distance from the point to the nearest of them in
    # the transverse plane.
    quad = only(b for b in nif_quads() if b.quad == bl.quad)
    qb = remake(quad; aim = (0.4, -0.3, 0.2), defocus = 1.5)
    function reference(x)
        best = (Inf, x, 0, 0.0)
        for (j, (kb, fb, ob)) in enumerate(Briareus._omega1_beamlets(qb, nothing, Briareus._nif_beam_table())),
            sh in (0.0, -1.0, 1.0)
            f = fb + sh * kb
            s = dot(x - f, kb)
            c = f + s * kb + c10.c0 * (1 + s / 7700) * ob
            h = c10.h0 * abs(1 - s / 460)
            u, v = ob, cross(kb, ob)
            cs = [c - h * u - h * v, c + h * u - h * v, c + h * u + h * v, c - h * u + h * v]
            inside = abs(dot(x - c, u)) <= h && abs(dot(x - c, v)) <= h
            ys = [(p = cs[i]; w = cs[mod1(i + 1, 4)] - p; p + clamp(dot(x - p, w) / dot(w, w), 0, 1) * w)
                  for i in 1:4]
            y = ys[argmin([norm(x - y) for y in ys])]
            d = (inside ? -1 : 1) * norm(x - y)
            d < best[1] && (best = (d, y, j, sh))
        end
        return best
    end
    rng = Random.MersenneTwister(7)
    c = Vec3(qb.aim) + 12.0 * quad.port
    pts = [c + Vec3(40 .* (rand(rng, 3) .- 0.5)...) for _ in 1:400]
    as = [omega1_closest_approach(t1, [qb], 1 => p) for p in pts]
    refs = reference.(pts)
    @test all(isapprox(a.distance, r[1]; atol = 1e-10) for (a, r) in zip(as, refs))
    @test all(norm(a.x_beam - r[2]) < 1e-9 && a.beamlet_index == r[3] && a.focus_shift == r[4]
              for (a, r) in zip(as, refs) if a.distance > 1e-6)
    @test all(abs(dot(a.x_hw - a.x_beam, a.beamlet[1])) < 1e-10 for a in as)
    @test any(a.distance < 0 for a in as) && any(count_ones(a.face_mask) == 2 for a in as) &&
          length(unique(a.beamlet_index for a in as)) == 4
    # the mask names the sides: the closest point lies on each side it names
    onside(a) = (kb = a.beamlet[1]; ob = a.beamlet[3]; f = a.beamlet[2] + a.focus_shift * kb;
                 cc = a.contour.c0 * (1 + a.s_star / 7700); h = a.contour.h0 * abs(1 - a.s_star / 460);
                 p = dot(a.x_beam - f, ob) - cc; q = dot(a.x_beam - f, cross(kb, ob));
                 all(abs(v) < 1e-9 for (bit, v) in ((0x01, p + h), (0x02, p - h), (0x04, q + h), (0x08, q - h))
                     if a.face_mask & bit != 0) && a.face_mask != 0)
    @test all(onside, as)

    @test_throws ArgumentError omega1_closest_approach(t1, [bl], first(x for x in r1.rows if x.rule === :tas_stay_in))
    @test_throws ArgumentError omega1_closest_approach(t1, [Beam("b"; theta = 30, phi = 0)], "b" => :own)
    @test_throws ArgumentError omega1_closest_approach(t1, [bl], "none" => :own)
end

@testset "1ω disposition of hardware" begin
    bl = first(nif_beams())
    k = direction(bl)
    pq = Briareus._nif_quad_port(bl.quad)
    ô = Briareus._own_side(bl.port, pq)
    q̂ = cross(k, ô)
    Rm = hcat(ô, q̂, k)
    plate(name, w, R, c) = Component(name, transform(rectangle_mesh(w, w; nu = 2, nv = 2), R, c), gold(); two_sided = true)
    Rt = Rm * [cosd(45) 0 sind(45); 0 1 0; -sind(45) 0 cosd(45)]
    cb, cu = 16.0 * ô - 5.0 * q̂ + 30.0 * k, 20.0 * ô + 5.0 * q̂ - 30.0 * k
    comps = [plate(:mirror, 4.0, Rm, 8.0 * ô),                         # faces the beamlet in the footprint
             plate(:tilted, 3.0, Rt, 8.0 * ô + 7.0 * q̂),               # reflects sideways
             plate(:near, 0.1, Rm, (4.45 - 0.5 - 0.05) * ô),           # 0.5 mm outside the near edge
             plate(:far, 0.1, Rm, -3.9 * ô),                           # the dark side
             Component(:ball, translate(sphere_mesh(4.0; ntheta = 12, nphi = 24), cb), gold();
                       surface = SphereSurface(4.0; center = cb)),
             Component(:lump, translate(sphere_mesh(6.0; ntheta = 12, nphi = 24), cu), gold())]
    t = Target(comps...)
    kw = (; open_ports = [pq, -k])
    d = omega1_disposition(t, [bl]; kw...)
    of(name) = [j for j in eachindex(t.tri_area) if t.components[t.tri_comp[j]].name === name && t.tri_area[j] > 0]
    sts(name) = Set(d.status[of(name)])
    @test length(d.status) == length(t.tri_area) && d.ports == [pq, -k]
    @test sts(:mirror) == Set([:fail]) && sts(:tilted) == Set([:specular_clear])
    @test sts(:near) == Set([:tight]) && sts(:far) == Set([:clear]) && sts(:ball) == Set([:curved_safe])
    @test :review in sts(:lump) && !(:fail in sts(:lump))
    @test all(d.beam[of(:mirror)] .== 1) && all(d.beamlet[of(:mirror)] .== 1) && all(d.margin[of(:far)] .== Inf)
    @test all(d.port[of(:mirror)] .== 2) && all(d.port[of(:far)] .== 0) && all(isnan, d.port_angle[of(:near)])
    @test all(d.R1[of(:ball)] .== 4) && all(d.R2[of(:ball)] .== 4) && all(isnan, d.R1[of(:near)])
    @test all(isinf, d.R1[of(:mirror)]) && !any(d.shadowed)
    @test occursin("fail", sprint(show, d))
    # the rows of talis_check give the same verdicts, margins and angles
    r = talis_check(t, [bl]; capture = (), rtol = 0, kw...)
    row(rule, c) = only(x for x in r.rows if x.rule === rule && x.component === c)
    @test row(:clearance_1w, :near).status === :pass && row(:clearance_1w, :far).status === :pass
    @test row(:management_1w, :mirror).status === :fail && row(:management_1w, :tilted).status === :pass
    @test row(:management_1w, :ball).status === :pass && row(:management_1w, :lump).status === :review
    mn = row(:clearance_1w, :near).margin
    @test mn - 1e-9 <= minimum(d.margin[of(:near)]) <= mn + 2e-4 && 0 < mn < 1
    @test minimum(d.port_angle[of(:mirror)]) ≈ row(:management_1w, :mirror).worst atol = 1e-9
    @test minimum(d.port_angle[of(:tilted)]) ≈ row(:management_1w, :tilted).worst atol = 1e-9
    @test minimum(d.port_angle[of(:lump)][d.status[of(:lump)] .=== :review]) ≈ row(:management_1w, :lump).worst atol = 1e-9
    # a smaller tight threshold clears the near plate; the diagnostic rule fails the far one,
    # 8.35 mm from the footprint where 10 mm are required
    @test Set(omega1_disposition(t, [bl]; tight_threshold = 0.1, kw...).status[of(:near)]) == Set([:clear])
    dd = omega1_disposition(t, [bl]; components = (:far, :mirror), diagnostics = (:far,), kw...)
    @test Set(dd.status[of(:far)]) == Set([:fail]) && Set(dd.status[of(:ball)]) == Set([:unchecked])
    @test all(dd.margin[of(:far)] .< 0) && Set(dd.status[of(:mirror)]) == Set([:fail])
    # a shield upstream of the mirror shadows it and nothing else in the footprint
    ts = Target(comps..., plate(:shield, 5.0, Rm, 8.0 * ô - 3.0 * k))
    ds = omega1_disposition(ts, [bl]; shadows = true, kw...)
    ofs(name) = [j for j in eachindex(ts.tri_area) if ts.components[ts.tri_comp[j]].name === name && ts.tri_area[j] > 0]
    @test all(ds.shadowed[ofs(:mirror)]) && !any(ds.shadowed[ofs(:shield)]) && !any(ds.shadowed[ofs(:tilted)])
    @test Set(ds.status[ofs(:mirror)]) == Set([:fail])                  # shadowing does not change a verdict
    # the whole quad: the verdict of a triangle is the most severe over its beamlets
    quad = only(b for b in nif_quads() if b.quad == bl.quad)
    dq = omega1_disposition(t, [quad, bl]; kw...)
    rank = Briareus._DISPOSITION_RANK
    @test all(rank[a] >= rank[b] for (a, b) in zip(dq.status, d.status)) && all(dq.margin .<= d.margin)
    @test_throws ArgumentError omega1_disposition(t, [bl]; tight_threshold = 6.0)
    # beams without a 1ω footprint leave everything clear
    @test Set(omega1_disposition(t, [Beam("b"; theta = 30, phi = 0)]).status) == Set([:clear])
end

@testset "1ω slice" begin
    bl = first(nif_beams())
    k = direction(bl)
    pq = Briareus._nif_quad_port(bl.quad)
    ô = Briareus._own_side(bl.port, pq)
    q̂ = cross(k, ô)
    # the square of the analytic case, in the focal plane, and a plate along the beam that
    # crosses the slice plane in a segment at q = 0.8
    own = Component(:own, transform(rectangle_mesh(0.1, 0.1; nu = 2, nv = 2), hcat(ô, q̂, k), 3.9 * ô), gold();
                    two_sided = true)
    fin = Component(:fin, transform(rectangle_mesh(1.0, 6.0; nu = 2, nv = 3), hcat(ô, k, cross(ô, k)),
                                    3.0 * ô + 0.8 * q̂), gold(); two_sided = true)
    t = Target(own, fin)
    a = omega1_closest_approach(t, [bl], bl.name => :own; rtol = 0)
    sl = omega1_slice(t, a)
    c10 = omega1_contour()
    taper = c10.c0 / 7700 + c10.h0 / 460
    # at s = -1 the near edge has moved toward the 3ω axis by the taper
    @test a.s_star ≈ -1 atol = 1e-12
    @test sl.nominal.center[1] - sl.nominal.half ≈ 4.45 - taper atol = 1e-12
    @test sl.nominal.center[2] == 0 && sl.worst.center == sl.nominal.center
    b = a.budget
    @test sl.worst.half - sl.nominal.half ≈ b.pointing + b.target_offset + b.rotation_tol + b.required atol = 1e-12
    hw, xb = sl.caliper
    @test hypot(hw[1] - xb[1], hw[2] - xb[2]) ≈ abs(a.distance) atol = 1e-12
    @test xb[1] ≈ 4.45 - taper atol = 1e-12
    @test hw[1] ≈ 3.95 atol = 1e-12
    @test sl.distance == a.distance && sl.margin == a.margin && sl.beam == bl.name && sl.component === :own
    @test sl.window[2] - sl.window[1] ≈ 4 && (sl.window[1] + sl.window[2]) / 2 ≈ (hw[1] + xb[1]) / 2
    # the own square lies in the plane (its edges); the fin crosses it from p = 2.5 to 3.5
    cut = [s for s in sl.hardware if all(u -> abs(u[2] - 0.8) < 1e-9, s)]
    flat = [s for s in sl.hardware if all(u -> abs(u[2]) <= 0.05 + 1e-9, s)]
    @test length(cut) + length(flat) == length(sl.hardware) && !isempty(flat)
    @test minimum(minimum(u[1] for u in s) for s in cut) ≈ 2.5 atol = 1e-9
    @test maximum(maximum(u[1] for u in s) for s in cut) ≈ 3.5 atol = 1e-9
    @test sum(abs(s[1][1] - s[2][1]) for s in cut) ≈ 1.0 atol = 1e-9
    @test all(3.85 - 1e-9 <= u[1] <= 3.95 + 1e-9 for s in flat for u in s)
    @test maximum(u[1] for s in flat for u in s) ≈ 3.95 && minimum(u[1] for s in flat for u in s) ≈ 3.85
    # a narrow window leaves the fin out
    @test all(s -> all(u -> abs(u[2]) < 0.1, s), omega1_slice(t, a; window = 0.3).hardware)
    @test_throws ArgumentError omega1_slice(t, a; window = 0.0)
    # the drawing: the same text for the same slice, with every part in it
    svg = omega1_slice_svg(sl)
    @test svg == omega1_slice_svg(omega1_slice(t, a)) && startswith(svg, "<svg ") && endswith(svg, "</svg>\n")
    @test count("<rect", svg) == 3 && count("<circle", svg) == 2 && count("<line", svg) == length(sl.hardware) + 2
    @test occursin(Briareus.Printf.@sprintf("%.3f mm</text>", a.distance), svg) && occursin("1 mm</text>", svg)
    @test occursin("stroke-dasharray", svg) && occursin(bl.name * " / own", svg)
    # the near edge of the contour is drawn at its place: x = (p - pmin) · size / width
    xe = (4.45 - taper - sl.window[1]) * 480 / 4
    @test occursin(Briareus.Printf.@sprintf("<rect x=\"%.3f\"", xe), svg)
    path = write_omega1_slice_svg(joinpath(mktempdir(), "slice.svg"), t, a; size = 300)
    @test read(path, String) == omega1_slice_svg(sl; size = 300)
    @test_throws ArgumentError omega1_slice_svg(sl; size = 50)
end
