# The model tree of assemblies (phase B2): membership separate from placement, nested
# assemblies, subtree show/hide and include/exclude with the mixed state, the tree's
# commands (new, move into, ungroup, rename, remove), the analyses seeing the flattened
# components (the drive unchanged, bit for bit), sessions (round trip, version 9).

const A3 = Briareus.Vec3

# The drive's numbers that must not change when only the tree changes.
_drive_bits(sess) = (r = (build_scenario!(sess); sess.drive.result[]);
                     (r.mean_flux, [m.epsilon for m in r.modes], r.solution.B, r.solution.H,
                      r.evaluation.deposition.element))

@testset "assemblies: membership, nesting, the drive unchanged (no GL)" begin
    s = DesignSession(resolution = :coarse)
    ts, bt, v = s.target, s.beams, s.view
    before = _drive_bits(s)
    fr0, t0 = part_frames(ts), ts.target[]
    # grouping the starter: the hohlraum (placed in TCC, its assembly's frame at the top
    # level) is placed in the new assembly's frame; the capsule, placed in the
    # hohlraum's, keeps its placement
    @test add_assembly!(ts; name = :cavity, members = [:hohlraum, :capsule])
    ps = ts.parts[]
    hp, cp = ps[part_index(ts, :hohlraum)], ps[part_index(ts, :capsule)]
    @test hp.assembly === cp.assembly === :cavity && hp.frame === :cavity && cp.frame === :hohlraum
    @test assemblies(ts) == [ps[end]] && is_assembly(ps[end]) && ps[end].assembly === TOP_LEVEL
    @test assembly_parts(ps, :cavity) == [:hohlraum, :capsule]
    # nothing moved: the same frames, the same components, so the same target
    @test all(n -> part_frames(ts)[n] == fr0[n], keys(fr0))
    @test ts.target[] === t0 && target_summary(ts) == target_summary(TargetSetup(resolution = :coarse))
    # the analyses see the flattened components: the drive is the same, bit for bit
    @test _drive_bits(s) == before

    # nested, and both at the top: a disk beside them
    @test add_part!(ts, :disk; name = :shield, x = 4.0, theta = 90.0, r1 = 0.5)
    @test add_assembly!(ts; name = :outer, members = [:cavity])
    @test ts.parts[][part_index(ts, :cavity)].assembly === :outer
    @test [p.name for p in assembly_children(ts.parts[], TOP_LEVEL)] == [:outer, :shield]
    @test assembly_members(ts.parts[], :outer) == [:cavity, :hohlraum, :capsule]
    # the tree: the assemblies (expanded), the hohlraum's pieces under it, the disk at the top
    v.expanded[] = union(v.expanded[], ("assembly:outer", "assembly:cavity"))
    rows = model_tree(s)
    depth(id) = only(r for r in rows if r.id == id).depth
    @test (depth("assembly:outer"), depth("assembly:cavity"), depth("part:hohlraum"),
           depth("comp:hohlraum_wall"), depth("part:capsule"), depth("part:shield")) == (1, 2, 3, 4, 3, 1)
    ids = [r.id for r in rows]
    @test findfirst(==("assembly:outer"), ids) < findfirst(==("part:hohlraum"), ids) < findfirst(==("part:shield"), ids)
    @test only(r for r in rows if r.id == "assembly:outer").detail == "2 parts"
    # collapsed, its contents are not listed
    toggle_expanded!(v, "assembly:outer")
    rows = model_tree(s)
    @test !any(r -> r.id == "part:capsule", rows) && any(r -> r.id == "part:shield", rows)
    toggle_expanded!(v, "assembly:outer")

    # moving the assembly moves its members and the beams aimed in their frames
    idx = findall(==(:inner23), bt.cluster[])
    set_aim_frame!(s, idx, Symbol("hohlraum.leh_top"))
    a0 = [bt.beams[][i].aim for i in idx]
    @test edit_part!(ts, part_index(ts, :outer), :x, "0.5")
    @test part_frames(ts)[:hohlraum].origin ≈ A3(0.5, 0, 0) atol = 1e-15
    @test part_frames(ts)[:capsule].origin ≈ A3(0.5, 0, 0) atol = 1e-15
    @test all(isapprox.([bt.beams[][i].aim for i in idx], [a + A3(0.5, 0, 0) for a in a0]; atol = 1e-14))
    @test part_frames(ts)[:shield].origin == A3(4, 0, 0)
    # (an assembly sweeps like a part: it is a sweep subject)
    @test :outer in sweep_subjects(s)
    @test edit_part!(ts, part_index(ts, :outer), :x, "0")

    # what an assembly has not: geometry, a role, surface fields
    i = part_index(ts, :outer)
    @test !set_part_role!(ts, i, :foil) && occursin("no surface", ts.error[])
    @test !edit_part!(ts, i, :albedo, "0.5") && occursin("no field", ts.error[])
    # the tree is a tree: no assembly inside itself, no user frame in one
    @test !move_into_assembly!(ts, [:outer], :cavity) && occursin("inside", ts.error[])
    @test !move_into_assembly!(ts, [:cavity], :cavity) && occursin("itself", ts.error[])
    @test add_frame!(ts; name = :mount)
    @test !move_into_assembly!(ts, [:mount], :outer) && occursin("frame", ts.error[])
    @test !move_into_assembly!(ts, [:shield], :nowhere) && occursin("no assembly", ts.error[])
    @test_throws ArgumentError check_tree([PartSpec(:sphere; name = :a, assembly = :b)])
end

@testset "assemblies: show/hide, include/exclude, roles (no GL)" begin
    s = DesignSession(resolution = :coarse)
    ts, v = s.target, s.view
    @test add_part!(ts, :disk; name = :shield, x = 4.0, theta = 90.0, r1 = 0.5)
    @test add_assembly!(ts; name = :cavity, members = [:hohlraum, :capsule])
    v.expanded[] = union(v.expanded[], ("assembly:cavity",))
    row(id) = only(r for r in model_tree(s) if r.id == id)
    @test row("assembly:cavity").on == 1 && row("assembly:cavity").display == "opaque"
    @test Set(node_components(s, "assembly:cavity")) == Set([:capsule; [piece_name(:hohlraum, p) for p in
                                                                        hohlraum_pieces(ts.parts[][1]) if !occursin("lip", String(p))]])
    # partly hidden: the mixed state, on the assembly and on the hohlraum part
    set_display!(v, :capsule, :invisible)
    @test row("assembly:cavity").on == 2 && row("assembly:cavity").display == "mixed"
    @test row("part:hohlraum").on == 1
    set_display!(v, :hohlraum_wall, :transparent)
    @test row("part:hohlraum").display == "opaque"                   # shown, in several modes
    set_display!(v, node_components(s, "assembly:cavity"), :invisible)
    @test row("assembly:cavity").on == 0 && row("part:hohlraum").on == 0
    set_display!(v, node_components(s, "assembly:cavity"), :opaque)

    # excluded as a subtree: its parts leave the calculations and give up their roles
    @test observed_name(ts) === :capsule && :hohlraum_wall in wall_names(ts)
    include_components!(s, [:cavity], false)
    @test ts.excluded[] == Set([:hohlraum, :capsule])
    @test component_names(ts.target[]) == [:shield]
    @test observed_name(ts) === nothing && isempty(wall_names(ts))
    @test occursin("not in calculations", row("assembly:cavity").detail)
    # a part put back alone: the assembly is partly in the calculations
    include_components!(s, [:capsule], true)
    @test observed_name(ts) === :capsule && occursin("partly in calculations", row("assembly:cavity").detail)
    # a piece left out on its own comes back with its assembly
    include_components!(s, [:hohlraum], true)
    @test set_excluded!(ts, :hohlraum_wall, true)
    @test occursin("partly", row("part:hohlraum").detail)
    include_components!(s, [:cavity], true)
    @test isempty(ts.excluded[]) && ts.ghost[] === nothing && row("assembly:cavity").detail == "2 parts"
    @test observed_name(ts) === :capsule && :hohlraum_wall in wall_names(ts)
    # set_excluded! takes an assembly too
    @test set_excluded!(ts, :cavity, true) && ts.excluded[] == Set([:hohlraum, :capsule])
    @test set_excluded!(ts, :cavity, false) && isempty(ts.excluded[])
end

@testset "assemblies: ungroup, rename, remove, sessions (no GL)" begin
    s = DesignSession(resolution = :coarse)
    ts, bt, v = s.target, s.beams, s.view
    @test add_part!(ts, :disk; name = :shield, x = 4.0, theta = 90.0, r1 = 0.5)
    # an assembly placed away from TCC, turned: its members stay where they are
    @test add_assembly!(ts; name = :addon)
    @test edit_part!(ts, part_index(ts, :addon), :x, "1.0") && edit_part!(ts, part_index(ts, :addon), :theta, "90")
    F = part_frames(ts)[:shield]
    @test move_into_assembly!(ts, [:shield], :addon)
    p = ts.parts[][part_index(ts, :shield)]
    @test p.assembly === :addon && p.frame === :addon
    @test isapprox(part_frames(ts)[:shield], F; atol = 1e-14)
    # rename: the member's membership and placement, the beams aimed in it, the selection
    idx = findall(==(:outer50), bt.cluster[])[1:2]
    set_aim_frame!(s, idx, :addon)
    select!(v, "assembly:addon")
    toggle_expanded!(v, "assembly:addon")
    @test edit_part!(ts, part_index(ts, :addon), :name, "mount")
    p = ts.parts[][part_index(ts, :shield)]
    @test p.assembly === :mount && p.frame === :mount && all(==(:mount), bt.aim_frame[][idx])
    @test v.selection[] == "assembly:mount" && "assembly:mount" in v.expanded[]
    # the session stores the tree, and reads it back the same
    state = session_state(s)
    @test only(q for q in state["target"]["parts"] if q["name"] == "shield")["assembly"] == "mount"
    s2 = DesignSession(resolution = :coarse)
    restore!(s2, state)
    @test session_state(s2) == state
    @test [q.assembly for q in s2.target.parts[]] == [q.assembly for q in ts.parts[]]
    # a version-9 session (no assemblies, no `assembly` keys) loads as it is, at the top level
    s9 = DesignSession(resolution = :coarse)
    add_part!(s9.target, :disk; name = :shield, x = 4.0, theta = 90.0, r1 = 0.5)
    old = session_state(s9)
    old["version"] = 9
    foreach(q -> delete!(q, "assembly"), old["target"]["parts"])
    s3 = DesignSession(resolution = :coarse)
    restore!(s3, old)
    @test all(q -> q.assembly === TOP_LEVEL, s3.target.parts[]) && session_state(s3) == session_state(s9)
    # a starter inserted in an assembly's frame from the window joins it; here by hand
    @test insert_starter!(ts, :foil; frame = :mount) && move_into_assembly!(ts, [:foil, :main], :mount)
    @test assembly_parts(ts.parts[], :mount) == [:shield, :foil, :main]
    # ungroup: the members stay where they are, now at the top level and placed in TCC
    Fs = Dict(n => part_frames(ts)[n] for n in (:shield, :foil, :main))
    @test ungroup_assembly!(ts, :mount)
    @test isempty(assemblies(ts)) && v.selection[] == "target"
    @test all(n -> (q = ts.parts[][part_index(ts, n)]; q.assembly === TOP_LEVEL && q.frame === :tcc), keys(Fs))
    @test all(n -> isapprox(part_frames(ts)[n], Fs[n]; atol = 1e-14), keys(Fs))
    @test all(==(:tcc), bt.aim_frame[][idx])                          # its beams keep their aims
    # remove a subtree: everything in it goes; a part outside it that was placed in its
    # frame stays where it is, at the top level's frame
    @test add_assembly!(ts; name = :pack, members = [:foil, :main])
    @test add_part!(ts, :sphere; name = :bead, frame = :pack, z = 1.0, r1 = 0.1)
    @test add_part!(ts, :sphere; name = :pearl, frame = :pack, z = -1.0, r1 = 0.1)
    Fp = part_frames(ts)[:pearl]
    @test replaced_parts(ts, [:pack]) == [:bead, :pearl]
    @test move_into_assembly!(ts, [:bead], :pack) && replaced_parts(ts, [:pack]) == [:pearl]
    @test remove_part!(ts, part_index(ts, :pack))
    @test !any(q -> q.name in (:pack, :foil, :main, :bead), ts.parts[])
    @test ts.parts[][part_index(ts, :pearl)].frame === :tcc && isapprox(part_frames(ts)[:pearl], Fp; atol = 1e-12)
    @test remove_part!(ts, part_index(ts, :pearl))
    @test :shield in component_names(ts.target[])
    # undone like any edit
    record!(s)
    @test add_assembly!(ts; name = :again, members = [:shield])
    record!(s)
    @test undo!(s) && isempty(assemblies(s.target)) && s.target.parts[][part_index(s.target, :shield)].assembly === TOP_LEVEL
end

@testset "assemblies: ungroup a nested assembly, and its refusals (no GL)" begin
    s = DesignSession(resolution = :coarse)
    ts = s.target
    @test add_part!(ts, :disk; name = :shield, x = 4.0, theta = 90.0, r1 = 0.5)
    @test add_part!(ts, :sphere; name = :bead, x = -4.0, r1 = 0.2)
    @test add_assembly!(ts; name = :inner, members = [:shield])
    @test add_assembly!(ts; name = :outer, members = [:inner])
    # placed away from TCC and turned, so the frames compose
    @test edit_part!(ts, part_index(ts, :outer), :x, "1.0") && edit_part!(ts, part_index(ts, :outer), :theta, "90")
    @test edit_part!(ts, part_index(ts, :inner), :z, "0.5")
    F = part_frames(ts)[:shield]
    # a placement step on another part refers to the inner assembly: ungrouping is refused
    ps = copy(ts.parts[])
    k = part_index(ts, :bead)
    ps[k] = Controllers._with(ps[k]; steps = Briareus.PlacementStep[Briareus.OnAxis(:inner; axis = :z)])
    @test Controllers._try_build!(ts; parts = ps)
    @test !ungroup_assembly!(ts, :inner) && occursin("bead", ts.error[]) && part_index(ts, :inner) !== nothing
    ps = copy(ts.parts[])
    ps[k] = Controllers._with(ps[k]; steps = Briareus.PlacementStep[])
    @test Controllers._try_build!(ts; parts = ps)
    # the member joins the parent assembly, placed in its frame, where it was in the chamber
    @test ungroup_assembly!(ts, :inner)
    p = ts.parts[][part_index(ts, :shield)]
    @test part_index(ts, :inner) === nothing && p.assembly === :outer && p.frame === :outer
    @test isapprox(part_frames(ts)[:shield].R, F.R; atol = 1e-14) &&
          isapprox(part_frames(ts)[:shield].origin, F.origin; atol = 1e-14)
    @test assembly_parts(ts.parts[], :outer) == [:shield]
    # a part, or a name of nothing, is not an assembly
    @test !ungroup_assembly!(ts, :shield) && occursin("no assembly", ts.error[])
    @test !ungroup_assembly!(ts, :nowhere) && occursin("no assembly", ts.error[])
end
