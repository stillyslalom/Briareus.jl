# The free-space enclosure of a point: a thick-walled hohlraum built as a closed solid,
# with a shield, a recess behind it and a mount outside.

# axis-aligned box faces (one patch each) facing out of the box, or into it when `inward`
function box_faces(lo, hi; inward = false, open_top = false)
    m = Mesh()
    c = (Vec3(lo...) + Vec3(hi...)) / 2
    x, y, z = (lo[1], hi[1]), (lo[2], hi[2]), (lo[3], hi[3])
    quads = [([Vec3(x[s], y[1], z[1]), Vec3(x[s], y[2], z[1]), Vec3(x[s], y[2], z[2]), Vec3(x[s], y[1], z[2])], Vec3(s == 1 ? -1 : 1, 0, 0)) for s in 1:2]
    append!(quads, [([Vec3(x[1], y[s], z[1]), Vec3(x[2], y[s], z[1]), Vec3(x[2], y[s], z[2]), Vec3(x[1], y[s], z[2])], Vec3(0, s == 1 ? -1 : 1, 0)) for s in 1:2])
    for s in (open_top ? (1,) : (1, 2))
        push!(quads, ([Vec3(x[1], y[1], z[s]), Vec3(x[2], y[1], z[s]), Vec3(x[2], y[2], z[s]), Vec3(x[1], y[2], z[s])], Vec3(0, 0, s == 1 ? -1 : 1)))
    end
    for (q, n) in quads
        add_polygon!(m, q; facing = inward ? -n : n)
    end
    return m
end

# Inner radius r, outer R, cavity |z| ≤ zc, end walls zc ≤ |z| ≤ zo with holes of radius a.
function thick_hohlraum(; r = 2.0, R = 2.3, zc = 4.0, zo = 4.3, a = 1.0, nphi = 24, mount = true)
    au = gold()
    cs = [Component(:inner, cylinder_mesh(r, -zc, zc; nz = 12, nphi, inward = true), au),
          Component(:outer, cylinder_mesh(R, -zo, zo; nz = 12, nphi, inward = false), au),
          Component(:cap_in, merge_meshes(annulus_mesh(a, r, zc; nr = 3, nphi, up = false),
                                          annulus_mesh(a, r, -zc; nr = 3, nphi, up = true)), au),
          Component(:cap_out, merge_meshes(annulus_mesh(a, R, zo; nr = 3, nphi, up = true),
                                           annulus_mesh(a, R, -zo; nr = 3, nphi, up = false)), au),
          Component(:rim, merge_meshes(cylinder_mesh(a, zc, zo; nz = 1, nphi, inward = true),
                                       cylinder_mesh(a, -zo, -zc; nz = 1, nphi, inward = true)), au),
          Component(:capsule, sphere_mesh(0.8; ntheta = 8, nphi = 12), ablator()),
          # a two-sided shield below the capsule: the bottom cap is seen only past it
          Component(:shield, disk_mesh(1.4, -2.0; nr = 2, nphi = 12), au; two_sided = true),
          # a 0.2 mm well 0.3 mm deep on the bottom cap, under the shield: its floor is
          # one 0.2 mm patch that sees out only through the well's mouth
          Component(:well, box_faces((1.4, -0.1, -3.99), (1.6, 0.1, -3.69); inward = true, open_top = true), au),
          Component(:well_out, box_faces((1.4, -0.1, -3.99), (1.6, 0.1, -3.69); open_top = true), au)]
    if mount
        # a block on the axis beyond the top hole: its near face sees the cavity through
        # the hole (view factor ≈ (a/8)² ≈ 0.016 < 0.05), the rest nothing of it
        m = box_faces((-1.0, -1.0, 12.0), (1.0, 1.0, 14.0))
        near, far = split_mesh(m, c -> c[3] < 12.5)
        push!(cs, Component(:mount_near, near, au), Component(:mount, far, au))
    end
    return Target(cs)
end

@testset "enclosure" begin
    t = thick_hohlraum()
    enc = enclosure(t)                     # from the capsule center
    names(tt) = Set(c.name for c in tt.components)
    @test names(enc.target) == Set([:inner, :cap_in, :rim, :capsule, :shield, :well, :well_out])
    # every element of the kept components, both sides of the shield, nothing else
    for c in enc.target.components
        i = component_index(t, c.name)
        @test all(enc.kept[t.comp_elems[i]])
        @test npatches(c.mesh) == npatches(t.components[i].mesh)
    end
    @test count(enc.kept) == nelements(enc.target)
    @test enc.target.elem_centroid == t.elem_centroid[enc.index]
    @test enc.target.elem_normal == t.elem_normal[enc.index]
    shield = t.comp_elems[component_index(t, :shield)]
    @test any(e -> t.elem_normal[e][3] < 0, shield) && any(e -> t.elem_normal[e][3] > 0, shield)
    # the point sees out through the top hole only (the shield covers the bottom one):
    # (1 - cos θ)/2 with tan θ = a/zo = 1/4.3, less the mount's 2 × 2 mm face 12 mm away
    Ωmount = 4 * asin(1 / (1 + 12^2))
    @test enc.escape ≈ (1 - cos(atan(1 / 4.3))) / 2 - Ωmount / 4π rtol = 0.03
    @test enc.back == 1.0                  # every ray starts inside the capsule
    # open area: the two holes (radius 1), from the kept elements' escaping rays
    @test enc.open_area ≈ 2π rtol = 0.1
    # the mount face sees the cavity through the hole with F ≈ 0.016: its flux,
    # about 0.9 F of the cavity's, is below the threshold
    near = t.comp_elems[component_index(t, :mount_near)]
    @test 0.005 < only(enc.flux[near]) < 0.05
    @test all(<(0.05), enc.flux[.!enc.kept])
    @test all(>=(0.05), enc.flux[enc.kept])
    # the outer wall gets no flux at all, the outer caps only what the mount face returns
    # (0.9 × its own × its view of them); the cavity about the reference
    @test all(iszero, enc.flux[t.comp_elems[component_index(t, :outer)]])
    @test maximum(enc.flux[t.comp_elems[component_index(t, :cap_out)]]) < 0.1 * only(enc.flux[near])
    @test 0.3 < minimum(enc.flux[t.comp_elems[component_index(t, :inner)]])
    @test maximum(enc.flux[t.comp_elems[component_index(t, :inner)]]) < 1.5
    # deterministic, and the same set from a point in free space
    @test enclosure(t).kept == enc.kept
    e2 = enclosure(t, (0.0, 1.2, 1.5))
    @test e2.kept == enc.kept && e2.back == 0
    # threshold 0 keeps what is connected at all: the mount face seen through the hole
    e0 = enclosure(t; threshold = 0)
    @test all(e0.kept[near]) && :mount_near in names(e0.target)
    # (and what that face sees in turn, such as the mount's sides, may be; not its back)
    @test !e0.kept[last(t.comp_elems[component_index(t, :mount)])]
    # a bounding sphere stops the rays before the mount
    @test !any(enclosure(t; threshold = 0, radius = 8.0).kept[near])
    # components in, components out
    @test enclosure(t.components).kept == enc.kept
    # a point in nothing
    @test_throws ArgumentError enclosure(t, (0.0, 0.0, 50.0); radius = 5.0)

    # Radiosity: elements the cavity cannot see exchange nothing with it, so the kept
    # target reproduces the capsule flux of the full target; without the mount, to
    # round-off. The mount face returns the little the hole sends it (leak area).
    P = zeros(nelements(t))
    P[t.comp_elems[component_index(t, :inner)]] .= 1.0
    cap(sol, tt) = sum(sol.H[e] * tt.elem_area[e] for e in elements(tt, :capsule)) / sum(tt.elem_area[elements(tt, :capsule)])
    full = solve(RadiationModel(t), P)
    sub = solve(RadiationModel(enc.target), P[enc.index])
    tn = thick_hohlraum(mount = false)
    nomount = solve(RadiationModel(tn), P[1:nelements(tn)])
    @test cap(sub, enc.target) ≈ cap(nomount, tn) rtol = 1e-10
    @test sub.H ≈ nomount.H[enc.index] rtol = 1e-10
    δ = cap(full, t) / cap(sub, enc.target) - 1
    @test 0 < δ < 1e-3
    # all of the leak goes to the mount face (the only dropped element the cavity sees),
    # which would return 0.9 × its flux of it: little against the open holes
    @test enc.leak_area > 0
    @test enc.return_area ≈ 0.9 * only(enc.flux[near]) * enc.leak_area rtol = 1e-12
    @test enc.return_area < 1e-3 * enc.open_area

    # VisRad workspaces: the cavity of the synthetic workspace's hohlraum (centered at
    # x = -4 mm) leaves out the foil and the CAD tetrahedron beside it
    path = joinpath(mktempdir(), "work.vrw")
    write(path, synthetic_workspace())
    ws = read_visrad_workspace(path)
    tv = visrad_target(ws; inside = (-4.0, 0.0, 0.0))
    @test all(startswith(String(c.name), "Cyl_Hohlraum_0") for c in tv.components)
    @test length(tv.components) >= 3
    # CAD is selected on its oriented facets and only the kept facets are remeshed: from
    # a point above the tetrahedron, the faces it sees
    cad = ws.objects[findfirst(o -> o.kind === :mesh, ws.objects)]
    q = cad.position + Vec3(0.3, 0.2, 3.0)
    mesh_only = o -> o.kind === :mesh
    ec = enclosure(visrad_target(ws; select = mesh_only, orient = true), q)
    @test 0 < nelements(ec.target) < 4
    tr = visrad_target(ws; select = mesh_only, orient = true, remesh = 0.3, inside = q)
    @test sum(tr.elem_area) ≈ sum(ec.target.elem_area) rtol = 1e-10
    @test nelements(tr) > 4nelements(ec.target)
    @test_throws ArgumentError visrad_target(ws; coarsen = 0.3, remesh = 0.3, inside = q)
    # `fit` applies to the elements kept (a tetrahedron has no curved region: unchanged)
    tf = visrad_target(ws; select = mesh_only, orient = true, remesh = 0.3, inside = q, fit = true)
    @test [c.name for c in tf.components] == [c.name for c in tr.components] && tf.elem_area == tr.elem_area
end
