# Frames (`Frame`, `FrameDef`), frame pointing (`FramePointing`) and placement sweeps
# (`placement_sweep`), against closed forms and the code they replaced.

using Briareus, LinearAlgebra, ForwardDiff, Test
using Briareus: Vec3, SVector, SMatrix

@testset "frames" begin
    angles = [(θ, φ, ψ) for θ in (0.0, 23.5, 90.0, 137.0, 180.0) for φ in (0.0, 78.75, 211.0, -40.0)
              for ψ in (0.0, 17.0, 90.0)]
    @testset "frame algebra" begin
        # the z column is the port direction; the x and y columns θ̂ and φ̂ turned by ψ
        @test all(euler_rotation(θ, φ, 0) * Vec3(0, 0, 1) == port_vector(θ, φ) for (θ, φ, _) in angles)
        basis = [Briareus.port_basis(port_vector(θ, φ)) for (θ, φ, _) in angles]
        @test all(begin
                      F = port_frame(θ, φ; roll = ψ)
                      eθ, eφ = basis[i]
                      # at the poles port_basis takes φ from the port vector: compare there by z
                      (θ in (0.0, 180.0) || (F.R[:, 1] ≈ cosd(ψ) * eθ + sind(ψ) * eφ &&
                                             F.R[:, 2] ≈ -sind(ψ) * eθ + cosd(ψ) * eφ)) &&
                          F.R[:, 3] ≈ port_vector(θ, φ) && F.origin == zero(Vec3)
                  end for (i, (θ, φ, ψ)) in enumerate(angles))
        @test all(det(euler_rotation(a...)) ≈ 1 && euler_rotation(a...)' * euler_rotation(a...) ≈ I
                  for a in angles)
        # the GUI's part rotation (radians through cos), to round-off
        Rz(a) = [cos(a) -sin(a) 0; sin(a) cos(a) 0; 0 0 1]
        Ry(a) = [cos(a) 0 sin(a); 0 1 0; -sin(a) 0 cos(a)]
        part_rotation(θ, φ, ψ) = Rz(deg2rad(φ)) * Ry(deg2rad(θ)) * Rz(deg2rad(ψ))
        @test all(isapprox(euler_rotation(a...), part_rotation(a...); atol = 1e-15) for a in angles)
        # the VisRad workspace rotation it replaced, bit for bit
        vRz(a) = SMatrix{3,3,Float64}(cosd(a), sind(a), 0, -sind(a), cosd(a), 0, 0, 0, 1)
        vRy(a) = SMatrix{3,3,Float64}(cosd(a), 0, -sind(a), 0, 1, 0, sind(a), 0, cosd(a))
        @test all(euler_rotation(θ, φ, ψ) == vRz(φ) * vRy(θ) * vRz(ψ) for (θ, φ, ψ) in angles)

        F = euler_frame((0.3, -1.2, 2.5), 37.0, 211.0, 20.0)
        G = port_frame(64.0, 12.0; roll = -33.0, origin = (1.0, 2.0, -0.5))
        p = Vec3(0.7, -0.4, 1.9)
        @test isapprox(compose(F, inv(F)), CHAMBER; atol = 1e-15)
        @test isapprox(compose(inv(F), F), CHAMBER; atol = 1e-15)
        @test to_chamber(inv(F), p) ≈ to_frame(F, p) atol = 1e-15
        @test to_frame(F, to_chamber(F, p)) ≈ p atol = 1e-15
        @test to_chamber(compose(F, G), p) ≈ to_chamber(F, to_chamber(G, p)) atol = 1e-14
        @test to_chamber(F, p; direction = true) == F.R * p
        @test to_frame(F, F.R[:, 2]; direction = true) ≈ Vec3(0, 1, 0) atol = 1e-15
        @test to_chamber(F, (0.0, 0.0, 0.0)) == F.origin
        # duals pass through a Float64 frame: the Jacobian of a point is R
        @test ForwardDiff.jacobian(x -> to_chamber(F, x), [0.1, 0.2, 0.3]) == F.R
        # Frame(origin; z, x): x is the part of `x` normal to z; the default is θ̂ of z
        H = Frame((1, 2, 3); z = (0, 0, 2), x = (1, 1, 5))
        @test H.R ≈ [1/√2 -1/√2 0; 1/√2 1/√2 0; 0 0 1] && H.origin == Vec3(1, 2, 3)
        @test Frame((0, 0, 0); z = port_vector(40.0, 30.0)).R ≈ euler_rotation(40.0, 30.0)
        @test_throws ArgumentError Frame((0, 0, 0); z = (0, 0, 1), x = (0, 0, -3))
        @test_throws ArgumentError Frame((0, 0, 0); z = (0, 0, 0))
    end

    @testset "foil frame and backlighter_foil" begin
        # the pre-frame construction of backlighter_foil, kept as the reference
        function old_foil_vertices(los; standoff, width = 2.0, height = width, tilt = 0.0,
                                   tilt_toward = 0.0, rotation = 0.0, center = (0.0, 0.0, 0.0), n = 8)
            d = port_vector(los...)
            eθ, eφ = Briareus.port_basis(d)
            u = cosd(rotation) * eθ + sind(rotation) * eφ
            v = -sind(rotation) * eθ + cosd(rotation) * eφ
            e = cosd(tilt_toward) * eθ + sind(tilt_toward) * eφ
            β = deg2rad(tilt)
            tl(x) = Briareus._tilted(x, d, e, β)
            u, v, nrm = tl(u), tl(v), tl(d)
            m = rectangle_mesh(width, height; nu = n, nv = max(1, round(Int, n * height / width)))
            R = SMatrix{3,3,Float64}(u[1], u[2], u[3], v[1], v[2], v[3], nrm[1], nrm[2], nrm[3])
            return transform(m, R, Vec3(center...) - standoff * d).vertices
        end
        cases = [(l, tl, tt, r) for l in ((90.0, 78.0), (37.0, 211.0), (0.0, 0.0), (180.0, 10.0))
                 for tl in (0.0, 35.0) for tt in (0.0, 90.0) for r in (0.0, 17.0)]
        kw(c) = (; standoff = 6.5, width = 2.0, height = 1.5, tilt = c[2], tilt_toward = c[3],
                 rotation = c[4], center = (0.1, -0.2, 0.3), n = 4)
        @test all(backlighter_foil(c[1]; kw(c)...).mesh.vertices == old_foil_vertices(c[1]; kw(c)...)
                  for c in cases)
        # the frame: origin at the foil's center, z its front normal, x along its width
        F = foil_frame((37.0, 211.0); standoff = 6.0, tilt = 20.0, tilt_toward = 90.0)
        c = backlighter_foil((37.0, 211.0); standoff = 6.0, tilt = 20.0, tilt_toward = 90.0,
                             width = 2.0, height = 1.0, n = 4)
        t = Target([c])
        loc = [to_frame(F, v) for v in c.mesh.vertices]
        @test maximum(l -> abs(l[3]), loc) < 1e-14
        @test maximum(l -> abs(l[1]), loc) ≈ 1.0 && maximum(l -> abs(l[2]), loc) ≈ 0.5
        @test all(e -> t.elem_normal[e] ≈ F.R[:, 3], findall(==(1), t.elem_side))
        @test part_centroid(t, :backlighter) ≈ F.origin atol = 1e-14
        @test isapprox(component_frame(t, :backlighter), Frame(F.origin; z = F.R[:, 3]); atol = 1e-14)
    end

    @testset "component frames and coordinates" begin
        # a cylinder with its exact surface: the surface axis; without it: the principal
        # axis of its area moments (the symmetry axis of a tube longer than it is wide)
        ax = normalize(Vec3(1, 2, 2))
        Rk = Frame((0, 0, 0); z = ax).R
        tube = transform(cylinder_mesh(0.5, -2.0, 2.0; nz = 8, nphi = 24), Rk, Vec3(1, 1, 1))
        t1 = Target(Component(:tube, tube, gold(); surface = CylinderSurface(0.5; origin = (1, 1, 1), axis = ax)))
        t2 = Target(Component(:tube, tube, gold()))
        @test component_frame(t1, :tube).R[:, 3] == ax
        @test component_frame(t2, :tube).R[:, 3] ≈ ax atol = 1e-12
        @test component_frame(t2, :tube).origin ≈ Vec3(1, 1, 1) atol = 1e-12
        # a sphere has no axis: chamber z at its center
        ts = Target(Component(:s, sphere_mesh(0.7; ntheta = 12, nphi = 24), gold()))
        @test component_frame(ts, :s).R == CHAMBER.R
        @test part_centroid(ts)[:s] ≈ zero(Vec3) atol = 1e-14
        # coordinates: closed forms and round trips through a rotated frame
        @test all(coords_of((1.0, 1.0, √2), :spherical) .≈ (2.0, 45.0, 45.0))
        @test all(coords_of((-1.0, 0.0, 3.0), :cylindrical) .≈ (1.0, 180.0, 3.0))
        F = euler_frame((0.3, -1.2, 2.5), 37.0, 211.0, 20.0)
        p = Vec3(0.7, -0.4, 1.9)
        @test all(point_of(coords_of(p, s; frame = F), s; frame = F) ≈ p for s in (:cartesian, :cylindrical, :spherical))
        @test point_of((2.0, 90.0, 90.0), :spherical; frame = F) ≈ to_chamber(F, (0.0, 2.0, 0.0)) atol = 1e-14
        @test_throws ArgumentError point_of((-1.0, 0.0, 0.0), :cylindrical)
        @test_throws ArgumentError coords_of(p, :polar)
    end

    @testset "frame graph" begin
        defs = Dict(:gp => FrameDef(:tcc; offset = (1.0, 0.0, 0.5), rotation = (30.0, 40.0, 10.0)),
                    :parent => FrameDef(:gp; offset = (0.0, 2.0, 0.0), rotation = (90.0, 0.0, 0.0)),
                    :child => FrameDef(:parent; offset = (0.5, 0.0, -1.0), rotation = (0.0, 0.0, 45.0)),
                    :port => FrameDef(:tcc; rotation = (64.0, 12.0, 5.0)))
        local_frame(d) = euler_frame(d.offset, d.rotation...)
        fr = resolve_frames(defs)
        # the grandchild follows its grandparent exactly, through its parent
        @test fr[:child] == compose(compose(local_frame(defs[:gp]), local_frame(defs[:parent])),
                                    local_frame(defs[:child]))
        @test fr[:port] == port_frame(64.0, 12.0; roll = 5.0) && fr[:tcc] == CHAMBER
        # a base frame (an entrance hole) as the root
        leh = Frame(Diagonal([1.0, -1.0, -1.0]), (0.0, 0.0, -5.0))
        fb = resolve_frames([:hole => FrameDef(:leh_bottom; offset = (0.0, 0.0, 1.0))], Dict(:leh_bottom => leh))
        @test fb[:hole].origin == Vec3(0, 0, -6) && fb[:leh_bottom] == leh
        # a loop names its path
        loop = Dict(:a => FrameDef(:b), :b => FrameDef(:c), :c => FrameDef(:a), :tail => FrameDef(:a))
        err = try resolve_frames(loop); nothing catch e; e end
        @test err isa ArgumentError && occursin(r"loop: ([abc]) → [abc] → [abc] → \1$", err.msg) &&
              all(n -> occursin(n, err.msg), ("a", "b", "c")) && !occursin("tail", err.msg)
        @test_throws ArgumentError resolve_frames(Dict(:a => FrameDef(:nowhere)))
        @test_throws ArgumentError resolve_frames(Dict(:tcc => FrameDef(:tcc)))
        @test_throws ArgumentError resolve_frames([:a => FrameDef(), :a => FrameDef()])
    end

    @testset "placement steps" begin
        snout = Frame((10.0, -2.0, 5.0); z = (0.3, 0.1, 1.0))
        base = Dict(:snout => snout)
        start = FrameDef(:tcc; offset = (1.0, 3.0, -0.5), rotation = (40.0, 70.0, 20.0))
        F0 = euler_frame(start.offset, start.rotation...)
        place(steps...) = resolve_frames(Dict(:a => FrameDef(:tcc; offset = start.offset,
                                                             rotation = start.rotation, steps = collect(steps))), base)[:a]
        report(steps...) = placement_report(Dict(:a => FrameDef(:tcc; offset = start.offset,
                                                                rotation = start.rotation, steps = collect(steps))), base)[:a]
        q = to_chamber(snout, (0.5, 0.0, 1.0))
        a = snout.R[:, 3]
        # each step against its closed form; a turn keeps the axis across the old and new z
        # (the smallest rotation) and the origin
        Fc = place(Coincident(:snout; point = (0.5, 0.0, 1.0)))
        @test isapprox(Fc.origin, q; atol = 1e-14) && Fc.R == F0.R
        Fd = place(Distance(:snout, 2.5; point = (0.5, 0.0, 1.0)))
        @test isapprox(Fd.origin, q + 2.5 * normalize(F0.origin - q); atol = 1e-14) && Fd.R == F0.R
        Fo = place(OnAxis(:snout; distance = -1.5))
        @test Fo.origin ≈ snout.origin - 1.5a atol = 1e-14
        Fl = place(OnAxis(:snout))
        @test Fl.origin ≈ snout.origin + dot(F0.origin - snout.origin, a) * a atol = 1e-14
        k(F, t) = normalize(cross(F.R[:, 3], t))
        Ff = place(Facing(:snout; point = (0.5, 0.0, 1.0)))
        t = normalize(q - F0.origin)
        Q = Ff.R * F0.R'
        @test isapprox(Ff.R[:, 3], t; atol = 1e-14) && isapprox(Q * k(F0, t), k(F0, t); atol = 1e-14) &&
              Ff.origin == F0.origin
        @test place(Facing(:snout; side = :back)).R[:, 3] ≈ normalize(F0.origin - snout.origin) atol = 1e-14
        Fp = place(Parallel(:snout; flip = true))
        @test isapprox(Fp.R[:, 3], -a; atol = 1e-14) && isapprox(Fp.R * F0.R' * k(F0, -a), k(F0, -a); atol = 1e-14)
        Fk = place(Concentric(:snout))
        @test isapprox(Fk.R[:, 3], a; atol = 1e-14) && norm(cross(Fk.origin - snout.origin, a)) < 1e-13 &&
              isapprox(dot(Fk.origin - snout.origin, a), dot(F0.origin - snout.origin, a); atol = 1e-13)
        # antiparallel: half a turn about the frame's own x
        @test resolve_frames(Dict(:a => FrameDef(steps = [Parallel(:tcc; flip = true)])))[:a].R ≈
              Diagonal([1.0, -1.0, -1.0]) atol = 1e-15
        # a roll from chamber x about chamber z
        Fr = resolve_frames(Dict(:a => FrameDef(steps = [Roll(30.0; from = :tcc)])))[:a]
        @test Fr.R ≈ euler_rotation(0, 0, 30) atol = 1e-15
        Fz = place(Parallel(:snout), Roll(-25.0; from = (:snout, :y)))
        e = normalize(snout.R[:, 2] - dot(snout.R[:, 2], a) * a)
        @test Fz.R[:, 1] ≈ cosd(-25) * e + sind(-25) * cross(a, e) atol = 1e-14
        # relative moves are compositions in the frame's own coordinates
        @test place(Offset(0.2, -0.4, 1.1)) ≈ compose(F0, Frame(I(3), (0.2, -0.4, 1.1)))
        @test place(Offset(0.2, -0.4, 1.1; frame = :snout)).origin ≈ F0.origin + snout.R * Vec3(0.2, -0.4, 1.1)
        @test place(Rotate(:z, 35.0)) ≈ compose(F0, euler_frame((0, 0, 0), 0, 0, 35))
        @test place(Rotate(:y, 35.0)) ≈ compose(F0, euler_frame((0, 0, 0), 35, 0, 0))
        @test place(Roll(35.0)) ≈ place(Rotate(:z, 35.0))

        # degrees of freedom: what each step fixes, what stays free, redundancy
        r = report(Distance(:tcc, 2.0), Facing(:snout), Roll(30.0; from = :tcc))
        @test r.fixed == [(1, 0), (0, 2), (0, 1)] && length(r.free_translation) == 2 &&
              isempty(r.free_rotation) && isempty(r.redundant)
        rf = place(r.steps...)
        @test all(u -> abs(dot(u, normalize(rf.origin))) < 1e-14, r.free_translation)
        r = report(Concentric(:snout), Parallel(:snout), OnAxis(:snout; distance = 1.0))
        @test r.fixed == [(2, 2), (0, 0), (1, 0)] && r.redundant == [2] &&
              isempty(r.free_translation) && length(r.free_rotation) == 1 &&
              abs(dot(only(r.free_rotation), a)) ≈ 1
        @test report(Coincident(:tcc)).fixed == [(3, 0)] && length(report().free_rotation) == 3
        # a relative move after a relation leaves a residual; an unmoved relation none
        r = report(OnAxis(:snout; distance = 2.0), Facing(:tcc), Offset(0.5, 0.0, 0.0))
        @test r.residuals[1][1] ≈ 0.5 && r.residuals[2][2] > 1 && r.residuals[3] == (0.0, 0.0)
        @test all(x -> all(<(1e-9), x), report(Concentric(:snout), Roll(10.0; from = :tcc)).residuals)

        # a relation that undoes an earlier one names both steps and the change
        err = try place(Coincident(:tcc), Coincident(:tcc; point = (0.0, 2.0, 0.0))); nothing catch e; e end
        @test err isa ArgumentError &&
              occursin("step 2 (Coincident tcc) conflicts with step 1 (Coincident tcc)", err.msg) &&
              occursin("moves it by 2 mm", err.msg)
        err = try place(Facing(:snout), Parallel(:tcc)); nothing catch e; e end
        @test err isa ArgumentError && occursin("step 2 (Parallel tcc) conflicts with step 1 (Facing snout)", err.msg) &&
              occursin("°", err.msg)
        # moving along the line to the faced point keeps the facing: no conflict
        @test place(Facing(:snout), Distance(:snout, 1.0)).R[:, 3] ≈ normalize(snout.origin - F0.origin) atol = 1e-14
        # undefined directions, bad arguments
        @test_throws ArgumentError place(Coincident(:snout), Facing(:snout))
        @test_throws ArgumentError place(Parallel(:tcc), Roll(0.0; from = (:tcc, :z)))
        @test_throws ArgumentError Distance(:tcc, -1.0)
        @test_throws ArgumentError Facing(:tcc; side = :left)
        @test_throws ArgumentError Parallel(:tcc; axis = :w)
        # a frame follows the frames its steps name; they resolve first, and loops through
        # them are named
        defs = Dict(:mount => FrameDef(:tcc; offset = (0.0, 0.0, 8.0)),
                    :foil => FrameDef(:tcc; offset = (2.0, 0.0, 0.0), steps = [Facing(:mount)]))
        @test resolve_frames(defs)[:foil].R[:, 3] ≈ normalize(Vec3(-2.0, 0, 8)) atol = 1e-14
        err = try resolve_frames(Dict(:a => FrameDef(steps = [Facing(:b)]), :b => FrameDef(steps = [Facing(:a)]))); nothing catch e; e end
        @test err isa ArgumentError && occursin(r"loop: ([ab]) → [ab] → \1$", err.msg)
        @test_throws ArgumentError resolve_frames(Dict(:a => FrameDef(steps = [Facing(:nowhere)])))
    end

    @testset "FramePointing" begin
        quads = nif_quads(power = 8.0)[1:8]
        # chamber frame: the initial parameters reproduce the template aims exactly
        aimed = [remake(b; aim = (0.1i, -0.05i, 0.2 * sign(b.port[3]))) for (i, b) in enumerate(quads)]
        p0 = FramePointing(aimed, CHAMBER; by = :beam, fields = (:u, :v, :w))
        @test all(b.aim == a.aim for (b, a) in zip(apply(p0, initial_parameters(p0)), aimed))
        # per-beam radial frames (x along the port's azimuth, y tangential, z chamber z):
        # (u, v, w) are GroupPointing's (r, phi, z) without mirroring
        radial(b) = (φ = atan(b.port[2], b.port[1]); Frame([cos(φ) -sin(φ) 0; sin(φ) cos(φ) 0; 0 0 1], zero(Vec3)))
        pf = FramePointing(aimed, radial; fields = (:u, :v, :w, :power))
        pg = GroupPointing(aimed; fields = (:r, :phi, :z, :power), mirror = false)
        x = [(0.3, -0.2, 2.5, 7.0)[mod1(i, 4)] + 0.01i for i in 1:nparams(pf)]
        @test nparams(pf) == nparams(pg)
        @test all(isapprox(a.aim, b.aim; atol = 1e-14) && a.power == b.power
                  for (a, b) in zip(apply(pf, x), apply(pg, x)))
        # the non-aim fields go through GroupPointing's code: identical values
        tail = (:defocus, :power, :spot, :dlambda)
        ft, gt = FramePointing(aimed, radial; fields = tail), GroupPointing(aimed; fields = tail)
        xt = repeat([0.5, 7.5, 1.2, 0.8], nparams(ft) ÷ 4)
        @test all(a.power == b.power && a.defocus == b.defocus && a.spot == b.spot && a.dlambda == b.dlambda
                  for (a, b) in zip(apply(ft, xt), apply(gt, xt)))
        @test all(a.aim == b.aim for (a, b) in zip(apply(ft, xt), aimed))   # no aim field: kept exactly
        # a rotated frame: aim = R (u, v, w) + o exactly, and its Jacobian is R
        F = euler_frame((0.3, -1.2, 2.5), 37.0, 211.0, 20.0)
        pr = FramePointing(quads[1:1], F; fields = (:u, :v, :w))
        xr = [0.4, -0.3, 0.2]
        @test apply(pr, xr)[1].aim == F.R * SVector(xr...) + F.origin
        @test ForwardDiff.jacobian(x -> apply(pr, x)[1].aim, xr) == F.R
        # mirrored entrance-hole frames: one local aim per cone gives mirrored aims
        L = 10.0
        top = Frame(SMatrix{3,3,Float64}(I), (0.0, 0.0, L / 2))
        bottom = Frame(Diagonal([1.0, -1.0, -1.0]), (0.0, 0.0, -L / 2))
        pm = FramePointing(nif_quads(), b -> b.port[3] > 0 ? top : bottom; by = :group, fields = (:u, :v, :w))
        xm = collect(Float64, 1:nparams(pm)) ./ 10
        bm = apply(pm, xm)
        g = Dict(k => xm[3i-2:3i] for (i, k) in enumerate(pm.keys))
        @test all(b.aim == (b.port[3] > 0 ? top : bottom).R * SVector(g[b.group]...) +
                           (b.port[3] > 0 ? top : bottom).origin for b in bm)
        @test any(b -> b.port[3] < 0, bm) && any(b -> b.port[3] > 0, bm)
        @test all(b.port[3] > 0 || (b.aim[1], -b.aim[2], -b.aim[3]) == Tuple(bm[findfirst(c -> c.group === b.group && c.port[3] > 0, bm)].aim)
                  for b in bm)
        # relative: zero aim shifts reproduce the template (to round-off); names, bounds
        rel = FramePointing(aimed, F; by = :beam, fields = (:u, :v, :power), relative = true)
        x0 = initial_parameters(rel)
        @test all(iszero, x0[1:3:end]) && all(iszero, x0[2:3:end]) && x0[3:3:end] == [b.power for b in aimed]
        @test all(isapprox(a.aim, b.aim; atol = 1e-14) for (a, b) in zip(apply(rel, x0), aimed))
        @test parameter_names(rel)[1:3] == [aimed[1].name * ".u", aimed[1].name * ".v", aimed[1].name * ".power"]
        lo, hi = bounds(rel, x0; u = 0.5, power = 0.25, power_max = 9.0)
        @test lo[1] == -0.5 && hi[1] == 0.5 && lo[3] == 6.0 && hi[3] == 9.0
        @test parameters(p0, aimed[2].name => (; w = 1.5))[6] == 1.5
        @test_throws ArgumentError FramePointing(aimed, CHAMBER; fields = (:z,))
        @test_throws DimensionMismatch FramePointing(aimed, [CHAMBER])
    end

    @testset "FramePointing in a scenario" begin
        # a backlighter foil heated by four beams aimed in its frame
        los = (90.0, 0.0)
        F = foil_frame(los; standoff = 5.0)
        t = Target([backlighter_foil(los; standoff = 5.0, width = 1.5, n = 6)])
        ports = [(60.0, 20.0), (120.0, 20.0), (60.0, -20.0), (120.0, -20.0)]
        bs = [Beam("b$i"; theta = θ, phi = φ, power = 0.5, spot = SuperGaussianSpot(0.2; order = 4.0),
                   fnumber = 6.7, aim = to_chamber(F, (0.0, 0.0, 0.0)))
              for (i, (θ, φ)) in enumerate(ports)]
        p = FramePointing(bs, F; by = :beam, fields = (:u, :v))
        s = Scenario(t, p; deposition_kw = (; resolution = 0.6))
        front = [e for e in elements(t, :backlighter) if t.elem_side[e] == 1]
        function spread(ev)   # uniform absorbed flux over the front
            P = ev.deposition.element[front]
            A = t.elem_area[front]
            return rms_nonuniformity(P ./ A, A)^2 + (1 - sum(P) / ev.deposition.delivered)^2
        end
        x0 = [0.2, 0.1, -0.2, 0.1, 0.2, -0.1, -0.2, -0.1]
        r = optimize_pointing(s, spread; x0, iterations = 20)
        @test r.value < r.initial_value && r.names == parameter_names(p)
        @test r.lower == x0 .- 1 && r.upper == x0 .+ 1           # the default ±1 mm box
        @test evaluate(s, r.x).deposition.element == deposit(s.deposition, apply(p, r.x)).element
        sub = subset(s, [1, 3])
        @test sub.param isa FramePointing && nparams(sub.param) == 4 &&
              [b.aim for b in apply(sub.param, r.x[[1, 2, 5, 6]])] == [b.aim for b in apply(p, r.x)[[1, 3]]]
    end

    @testset "placement sweeps" begin
        # aims that follow move exactly by value · axis; the others stay
        los = (90.0, 0.0)
        F = foil_frame(los; standoff = 5.0)
        t = Target([backlighter_foil(los; standoff = 5.0, width = 1.5, n = 4),
                    Component(:capsule, sphere_mesh(0.4; ntheta = 6, nphi = 12), ablator())])
        bs = [Beam("b$i"; theta = θ, phi = 20.0, power = 0.5, aim = to_chamber(F, (0.1i, 0.0, 0.0)))
              for (i, θ) in enumerate((60.0, 120.0))]
        vals = [-1.0, 0.0, 0.75]
        sw = placement_sweep(t, bs, :backlighter, vals; frame = F, axis = :y, follow = [2])
        a = F.R[:, 2]
        @test [r.value for r in sw.rows] == vals && !sw.canceled
        @test all(r.beams[2].aim == bs[2].aim + r.value * a && r.beams[1].aim == bs[1].aim for r in sw.rows)
        @test all(isnan(r.min_clearance) && r.los_clear === missing && isnan(r.min_captured) for r in sw.rows)
        # a FramePointing follows with its frames: the parameters keep their meaning
        p = FramePointing(bs, F; by = :beam, fields = (:u, :v))
        d = Briareus._MovedParts(t, Set([:backlighter]), p.template, [true, true], p, 0.75 * a)
        step = d(1.0)
        @test all(isapprox(b.aim, c.aim; atol = 1e-14) for (b, c) in
                  zip(apply(step.pointing, initial_parameters(p)), step.beams))
        @test initial_parameters(step.pointing) ≈ initial_parameters(p) atol = 1e-14

        # a plate beside a beam, swept across it: the clearance changes one for one. Beam
        # along -x onto a target foil at x = -3; the plate (normal y) beside the path at
        # y = 1, the sweep moving it away along +y.
        tgt = foil(:tgt, (-3.0, 0.0, 0.0), (1, 0, 0); width = 2.0, n = 2)
        plate = foil(:plate, (2.0, 1.0, 0.0), (0, 1, 0); width = 1.0, n = 2)
        tp = Target([tgt, plate])
        beam = [Beam("B"; theta = 90.0, phi = 0.0, power = 1.0, spot = SuperGaussianSpot(0.1), fnumber = 8.0)]
        ckw = (; targets = (:tgt,), wings = (), resolution = 1e-7, rtol = 0.0, nangles = 4, nspot = 11)
        sp = placement_sweep(tp, beam, :plate, [0.0, 0.1, 0.3]; axis = :y, clearance_kw = ckw)
        c = [r.min_clearance for r in sp.rows]
        @test c[1] - c[2] ≈ -0.1 atol = 1e-6
        @test c[1] - c[3] ≈ -0.3 atol = 1e-6
        @test all(r.worst_component === :plate && r.worst_beam == "B" for r in sp.rows)

        # a foil approaching an annular aperture on its line of sight (+z): the bundle
        # from the foil's corners (half-width w, at distance D) through the rim of the
        # object's disk (radius a) crosses the aperture plane z = -X at radius
        # (1 - s) a + s w √2, s = X / D; it meets a rim of radius ρ when
        # s* = (ρ - a) / (w √2 - a), at D* = X / s*, a shift Δ* = D0 - D* toward it.
        D0, w, a, X, ρ = 8.0, 1.0, 0.5, 3.0, 1.0
        Dstar = X * (w * √2 - a) / (ρ - a)
        Δstar = D0 - Dstar
        fl = backlighter_foil((0.0, 0.0); standoff = D0, width = 2w, n = 4)
        ap = Component(:aperture_plate, annulus_mesh(ρ, 3.0, -X; nr = 1, nphi = 64), BLACK)
        ta = Target([fl, ap])
        Fa = foil_frame((0.0, 0.0); standoff = D0)
        sl = placement_sweep(ta, Beam{Float64}[], :backlighter, [0.0, Δstar - 1e-6, Δstar + 1e-6];
                             frame = Fa, axis = :z, los = (0.0, 0.0),
                             los_kw = (; object = nothing, radius = a))
        @test [r.los_clear for r in sl.rows] == [true, true, false]
        @test sl.rows[3].blocked_fraction > 0 && sl.rows[3].los.blocked_by[1] == (:aperture_plate => sl.rows[3].los.blocked)

        # a callback cancel returns the completed rows
        seen = Int[]
        sc = placement_sweep(ta, Beam{Float64}[], :backlighter, collect(0.0:0.5:3.0); frame = Fa,
                             callback = st -> (push!(seen, st.done); st.done == 2))
        @test sc.canceled && length(sc.rows) == 2 && seen == [0, 1, 2] && length(sc.values) == 7
        @test placement_sweep(d, [0.0]; callback = st -> true).canceled
        # CSV: unchecked cells empty
        path = tempname() * ".csv"
        write_sweep_csv(path, sl)
        lines = readlines(path)
        @test length(lines) == 4 && startswith(lines[1], "value,min_clearance_mm") &&
              split(lines[4], ",")[7] == "false" && split(lines[4], ",")[2] == ""
        rm(path)
    end

    @testset "sweep with capture and re-optimization" begin
        los = (90.0, 0.0)
        F = foil_frame(los; standoff = 5.0)
        t = Target([backlighter_foil(los; standoff = 5.0, width = 1.5, n = 6)])
        bs = [Beam("b$i"; theta = θ, phi = φ, power = 0.5, spot = SuperGaussianSpot(0.2; order = 4.0),
                   fnumber = 6.7, aim = to_chamber(F, (0.15 * (-1)^i, 0.1, 0.0)))
              for (i, (θ, φ)) in enumerate([(60.0, 20.0), (120.0, -20.0)])]
        p = FramePointing(bs, F; by = :beam, fields = (:u, :v))
        function spread(ev)
            tt = ev.scenario.target
            fr = [e for e in elements(tt, :backlighter) if tt.elem_side[e] == 1]
            P, A = ev.deposition.element[fr], tt.elem_area[fr]
            return rms_nonuniformity(P ./ A, A)^2
        end
        rows = Any[]
        sw = placement_sweep(t, p, :backlighter, [0.0, 1.0]; frame = F, axis = :z, follow = true,
                             capture = :backlighter, capture_kw = (; resolution = 0.6),
                             optimize = (; objective = spread, iterations = 10,
                                         scenario_kw = (; deposition_kw = (; resolution = 0.6))),
                             callback = st -> (push!(rows, st.row); false))
        r1, r2 = sw.rows
        # the whole design translates rigidly: the warm start is already the optimum's
        # neighborhood, and the objective and captured fraction carry over
        @test r1.objective_after <= r1.objective_before
        @test r2.objective_before ≈ r1.objective_after rtol = 1e-6
        @test 0.5 < r1.min_captured <= 1 + 1e-3 && isapprox(r2.min_captured, r1.min_captured; rtol = 1e-3)
        @test length(r1.x) == 4 && rows[1] === nothing && rows[3] === r2
        @test_throws ArgumentError placement_sweep(t, bs, :backlighter, [0.0];
                                                   optimize = (; objective = spread))
    end
end
