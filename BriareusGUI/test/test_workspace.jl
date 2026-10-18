# Session files and undo/redo (no GL).

@testset "session files (no GL)" begin
    s = DesignSession(resolution = :coarse)
    # a design touching every part of the state
    edit_target!(s.target, :length, "10.4")
    add_window!(s.target; phi = 40.0, width = 1.5)
    set_window_open!(s.target, 1, true)
    add_part!(s.target, :disk; z = 7.5)
    edit_part!(s.target, 1, :theta, "20")
    split_quad!(s.beams, "Q31T")
    set_cluster!(s.beams, findall(b -> b.group === :inner23, s.beams.beams[]), :backlighter)
    edit_beam_set!(s.beams, findall(==(:backlighter), s.beams.cluster[]), :z, "3.4")
    set_enabled!(s.beams, 3, false; by = :beam)
    s.optimize.iterations[] = 37
    s.optimize.weights[] = Dict(2 => 1.0, 4 => 0.5)
    s.optimize.fields[] = [:z, :r, :power]
    s.review.peak_fraction[] = 3e-4
    s.talis.tolerance[] = 0.7
    s.drive.spot_radius[] = 0.6
    s.view.beams[] = :volumes
    s.view.hidden[] = Set([:capsule])
    path = joinpath(mktempdir(), "design.briareus.toml")
    save_session(path, s)
    @test occursin("BriareusGUI session", read(path, String))

    # a fresh session restored from the file has the same design, bit for bit
    s2 = DesignSession()
    load_session!(s2, path)
    @test session_state(s2) == session_state(s)
    @test s2.beams.beams[] == s.beams.beams[]                # ports and aims exact
    @test s2.beams.enabled[] == s.beams.enabled[] && s2.beams.cluster[] == s.beams.cluster[]
    @test s2.target.windows[] == s.target.windows[] && s2.target.mesh[] == s.target.mesh[]
    @test [p.values for p in s2.target.parts[]] == [p.values for p in s.target.parts[]]
    @test nelements(s2.target.target[]) == nelements(s.target.target[])
    @test s2.optimize.weights[] == Dict(2 => 1.0, 4 => 0.5) && s2.optimize.iterations[] == 37
    @test s2.view.beams[] === :volumes && s2.view.hidden[] == Set([:capsule])
    @test count(b -> b.quad == "Q31T", s2.beams.beams[]) == 4    # the split quad
    @test !can_undo(s2)                                          # history starts at the load
    # results are not saved: the drive is computed again
    @test s2.drive.scenario[] === nothing

    # other files are refused
    bad = joinpath(mktempdir(), "other.toml")
    write(bad, "format = \"something else\"\n")
    @test_throws ArgumentError load_session!(DesignSession(resolution = :coarse), bad)
    write(bad, "format = \"BriareusGUI session\"\nversion = 99\n")
    @test_throws ArgumentError load_session!(DesignSession(resolution = :coarse), bad)
end

@testset "undo and redo (no GL)" begin
    s = DesignSession(resolution = :coarse)
    @test !can_undo(s) && !can_redo(s) && !undo!(s) && !redo!(s)
    L0 = s.target.params[][:length]
    # each recorded edit is one step
    edit_target!(s.target, :length, "11.0")
    @test record!(s)
    @test !record!(s)                                    # nothing new
    add_window!(s.target; phi = 90.0)
    record!(s)
    b0 = s.beams.beams[][1]
    edit_beams!(s.beams, 1, :power, "5.5"; by = :beam)
    record!(s)
    @test s.beams.beams[][1].power == 5.5
    # view changes (coloring, selection) are not edits
    s.view.beams[] = :volumes
    select!(s.view, "beams")
    @test !record!(s)
    # undo walks back: the power, the window, the length
    @test undo!(s) && s.beams.beams[][1] == b0
    @test undo!(s) && isempty(s.target.windows[])
    @test can_redo(s)
    @test undo!(s) && s.target.params[][:length] == L0
    @test !undo!(s) && !can_undo(s)
    @test s.view.beams[] === :volumes                     # the view stays as it is
    # redo walks forward again, and the target is rebuilt each time
    @test redo!(s) && s.target.params[][:length] == 11.0
    @test redo!(s) && length(s.target.windows[]) == 1
    @test :window_1 in component_names(s.target.target[])
    # a new edit drops the redo steps
    edit_target!(s.target, :capsule_radius, "0.8")
    @test !can_redo(s) || (record!(s); !can_redo(s))
    @test undo!(s) && s.target.params[][:capsule_radius] == 0.9
    @test undo!(s) && isempty(s.target.windows[])
    # an edit that is not recorded yet is undone too (undo records first)
    set_resolution!(s.target, :standard)
    @test undo!(s) && s.target.resolution[] === :coarse
    # stale analyses: a drive computed before an undo of the target is stale after it
    reset_history!(s)
    @test !can_undo(s)
end

@testset "STL mesh parts (no GL)" begin
    # a closed cube (2 × 2 × 2 file units), half its facets wound backward, ASCII STL
    corners = [Briareus.Vec3(x, y, z) for x in (-1, 1), y in (-1, 1), z in (-1, 1)]
    quads = ((1, 3, 4, 2), (5, 6, 8, 7), (1, 2, 6, 5), (3, 7, 8, 4), (1, 5, 7, 3), (2, 4, 8, 6))
    path = joinpath(mktempdir(), "shield.stl")
    open(path, "w") do io
        println(io, "solid shield")
        for (i, q) in enumerate(quads)
            a, b, c, d = corners[collect(q)]
            for tri in ((a, b, c), isodd(i) ? (a, c, d) : (a, d, c))
                println(io, "facet normal 0 0 0\nouter loop")
                foreach(v -> println(io, "vertex ", join(v, " ")), tri)
                println(io, "endloop\nendfacet")
            end
        end
        println(io, "endsolid")
    end
    s = DesignSession(resolution = :coarse)
    ts = s.target
    @test !add_part!(ts, :mesh)                                   # no file
    @test occursin("STL", ts.error[])
    @test add_part!(ts, :mesh; file = path, scale = 0.125, z = 7.5)
    name = ts.parts[][end].name
    t = ts.target[]
    es = elements(t, name)
    # a 0.25 mm cube at z = 7.5, every face out of the solid; with elements of the
    # coarse mesh's 0.6 mm each face is one element, never one across its edges
    @test length(es) == 6
    @test all(sum(t.elem_normal[e] .* (t.elem_centroid[e] - Briareus.Vec3(0, 0, 7.5))) > 0 for e in es)
    @test sum(t.elem_area[es]) ≈ 6 * 0.25^2 rtol = 1e-12
    # at 1 mm each face is remeshed into four elements of about the 0.6 mm size
    @test edit_part!(ts, length(ts.parts[]), :scale, "0.5")
    @test sum(ts.target[].elem_area[elements(ts.target[], name)]) ≈ 6 * 1.0^2 rtol = 1e-12
    @test length(elements(ts.target[], name)) == 24
    @test !edit_part!(ts, length(ts.parts[]), :scale, "-1")       # rejected
    # in a session file
    f = joinpath(mktempdir(), "with_stl.briareus.toml")
    save_session(f, s)
    s2 = DesignSession(resolution = :coarse)
    load_session!(s2, f)
    @test s2.target.parts[][end].file == path && s2.target.parts[][end].kind === :mesh
    @test nelements(s2.target.target[]) == nelements(ts.target[])
end

@testset "x-ray image (no GL)" begin
    s = DesignSession(resolution = :coarse)
    @test s.drive.image[] === nothing
    # looking down the axis through the top entrance hole: computes the drive first
    xray_image!(s; direction = (0, 0, -1), up = (0, 1, 0), width = 6.0, nx = 64)
    @test s.drive.view[] === :image && s.drive.result[] !== nothing
    im = s.drive.image[]
    @test im !== nothing && size(im.values) == (64, 64)
    @test extrema(im.x) == (-3.0, 3.0)
    seen = count(isfinite, im.values)
    # through the top hole only the capsule (0.9 mm) is seen: around it the rays leave
    # through the bottom hole, straight below
    @test 0.9 < seen / (π * 0.9^2 / 36 * 64^2) < 1.1
    @test all(>(0), filter(isfinite, im.values))
    # the image follows the drive: more power, brighter
    peak = maximum(filter(isfinite, im.values))
    edit_beams!(s.beams, 1, :power, string(2 * s.beams.beams[][1].power); by = :beam)
    for b in 2:length(s.beams.beams[])
        edit_beams!(s.beams, b, :power, string(2 * s.beams.beams[][b].power); by = :beam)
    end
    @test s.drive.image[] !== im
    @test maximum(filter(isfinite, s.drive.image[].values)) ≈ 2peak rtol = 1e-6
    set_view!(s.drive, :scene)
    @test s.drive.view[] === :scene
    @test_throws ArgumentError xray_image!(s; direction = (0, 0, 0), width = 1.0)
    @test_throws ArgumentError xray_image!(s; direction = (0, 0, 1), width = 0.0)
end

@testset "part frames (no GL)" begin
    s = DesignSession(resolution = :coarse)
    ts = s.target
    L = ts.params[][:length]
    center(name) = (t = ts.target[]; es = elements(t, name);
                    sum(t.elem_centroid[e] * t.elem_area[e] for e in es) / sum(t.elem_area[es]))
    normal(name) = (t = ts.target[]; t.elem_normal[first(elements(t, name))])
    # shields 0.5 mm outside each entrance hole, from the same placement
    @test add_part!(ts, :disk; name = :shield_top, anchor = :leh_top, z = 0.5, r1 = 2.0)
    @test add_part!(ts, :disk; name = :shield_bottom, anchor = :leh_bottom, z = 0.5, r1 = 2.0)
    @test center(:shield_top) ≈ Briareus.Vec3(0, 0, L / 2 + 0.5) atol = 1e-12
    @test center(:shield_bottom) ≈ Briareus.Vec3(0, 0, -L / 2 - 0.5) atol = 1e-12
    @test normal(:shield_top) ≈ Briareus.Vec3(0, 0, 1) && normal(:shield_bottom) ≈ Briareus.Vec3(0, 0, -1)
    # they follow the hohlraum's length
    @test edit_target!(ts, :length, "11.0")
    @test center(:shield_top)[3] ≈ 6.0 && center(:shield_bottom)[3] ≈ -6.0
    # to the chamber frame, staying put
    i = part_index(ts, :shield_bottom)
    @test set_part_anchor!(ts, i, :tcc)
    p = ts.parts[][i]
    @test p.anchor === :tcc && p.values[:z] ≈ -6.0 && p.values[:theta] ≈ 180
    @test center(:shield_bottom) ≈ Briareus.Vec3(0, 0, -6.0) atol = 1e-12
    @test normal(:shield_bottom) ≈ Briareus.Vec3(0, 0, -1)
    # a frame the template lacks is refused, and so is a template the frames need
    @test !set_part_anchor!(ts, i, :foil) && occursin("backlighter", ts.error[])
    @test !set_target_kind!(ts, :custom) && ts.kind[] === :hohlraum
    # in a session file
    f = joinpath(mktempdir(), "frames.briareus.toml")
    save_session(f, s)
    s2 = DesignSession(resolution = :coarse)
    load_session!(s2, f)
    @test [p.anchor for p in s2.target.parts[]] == [:leh_top, :tcc]
    @test s2.target.target[].elem_centroid == ts.target[].elem_centroid
    # a backlighter's foil frame: z along the foil's back normal, away from TCC
    b = DesignSession(kind = :foil, resolution = :coarse)
    @test add_part!(b.target, :disk; anchor = :foil, z = 0.3, r1 = 1.5, name = :substrate)
    d = b.target.params[][:distance]
    t = b.target.target[]
    es = elements(t, :substrate)
    @test sum(t.elem_centroid[e] for e in es) / length(es) ≈ Briareus.Vec3(d + 0.3, 0, 0) atol = 1e-9
    @test t.elem_normal[first(es)] ≈ Briareus.Vec3(1, 0, 0)
end

@testset "cutaway, exclusion and multi-selection (no GL)" begin
    s = DesignSession(resolution = :coarse)
    v, ts = s.view, s.target
    # the cut: quadrant facing −x −y by default, then a half facing +x
    @test v.cutaway[] && v.cut_width[] == 90 && v.cut_azimuth[] == 225
    @test is_cut(v, :wall, (-1.0, -1.0, 0.0)) && !is_cut(v, :wall, (1.0, -1.0, 0.0))
    @test !is_cut(v, :capsule, (-1.0, -1.0, 0.0))                # not a cut component
    set_cut!(v; width = 180, azimuth = 0, components = [:wall, :capsule])
    @test is_cut(v, :capsule, (1.0, 0.9, 0.0)) && !is_cut(v, :wall, (-1.0, 0.2, 0.0))
    @test_throws ArgumentError set_cut!(v; width = 0)
    v.cutaway[] = false
    @test !is_cut(v, :wall, (1.0, 0.0, 0.0))
    # leaving the capsule out of the calculations: gone from the target, drawn as a ghost,
    # still listed in the tree
    n0 = nelements(ts.target[])
    @test set_excluded!(ts, :capsule, true)
    @test !(:capsule in component_names(ts.target[])) && ts.ghost[] !== nothing
    @test :capsule in all_component_names(ts) && nelements(ts.target[]) < n0
    rows = model_tree(s)
    @test any(r -> r.id == "comp:capsule" && occursin("not in calculations", r.detail), rows)
    select!(v, "comp:capsule")
    edit_target!(ts, :length, "10.5")                              # a rebuild keeps both
    @test v.selection[] == "comp:capsule" && ts.ghost[] !== nothing
    f = joinpath(mktempdir(), "excluded.briareus.toml")
    save_session(f, s)
    s2 = DesignSession(resolution = :coarse)
    load_session!(s2, f)
    @test s2.target.excluded[] == Set([:capsule]) && s2.target.ghost[] !== nothing
    @test s2.view.cut_width[] == 180 && s2.view.cut_components[] == Set([:wall, :capsule])
    @test set_excluded!(ts, :capsule, false) && ts.ghost[] === nothing
    @test :capsule in component_names(ts.target[])
    for c in component_names(ts.target[])
        set_excluded!(ts, c, true)
    end
    @test occursin("every component", ts.error[])                   # the last one stays
    @test length(component_names(ts.target[])) == 1
    for c in collect(ts.excluded[])
        set_excluded!(ts, c, false)
    end
    @test ts.ghost[] === nothing
    # Ctrl adds components to components (and beams to beams), Shift a range in tree order
    select!(v, "comp:wall")
    select_add!(v, "comp:cap_top")
    @test selected_components(s) == [:wall, :cap_top]
    select_add!(v, "beams")                                          # another family: alone
    @test selected_ids(v) == ["beams"]
    order = [r.id for r in model_tree(s)]
    select!(v, "comp:wall")
    select_range!(v, order, "comp:capsule")
    sel = selected_components(s)
    @test first(sel) === :wall && :capsule in sel && :cap_bottom in sel && length(sel) >= 4
    # hidden beams are a view setting, separate from switching them off
    set_beams_hidden!(v, [b.name for b in s.beams.beams[][1:3]], true)
    @test length(v.hidden_beams[]) == 3 && all(s.beams.enabled[])
    set_beams_hidden!(v, [s.beams.beams[][1].name], false)
    @test length(v.hidden_beams[]) == 2
end

@testset "exported results and facilities (no GL)" begin
    s = DesignSession(resolution = :coarse)
    @test_throws ArgumentError export_results(s, mktempdir())        # no drive yet
    build_scenario!(s)
    dir = mktempdir()
    files = export_results(s, dir)
    @test Set(basename.(files)) == Set(["elements.csv", "components.csv", "beams.csv", "legendre.csv",
                                        "ylm.csv", "summary.txt", "session.briareus.toml"])
    @test length(readlines(joinpath(dir, "elements.csv"))) == nelements(s.target.target[]) + 1
    s2 = DesignSession(resolution = :coarse)
    load_session!(s2, joinpath(dir, "session.briareus.toml"))
    @test session_state(s2) == session_state(s)
    # facilities: presets set them, the facility gates the TaLIS rules
    @test s.beams.facility[] === :nif && :talis in session_tasks(s)
    set_facility!(s.beams, :omega)
    @test s.beams.preset[] === :omega_hohlraum && s.beams.facility[] === :omega
    @test !(:talis in session_tasks(s))
    @test !any(r -> r.id == "task:talis", model_tree(s))
    set_facility!(s.beams, :omega_ep)                                  # no presets: beams kept
    @test s.beams.preset[] === :omega_hohlraum && s.beams.facility[] === :omega_ep
    load_preset!(s.beams, :nif_quads)
    @test s.beams.facility[] === :nif
    @test_throws ArgumentError set_facility!(s.beams, :lmj)
end

@testset "coordinate systems and aim references (no GL)" begin
    s = DesignSession(resolution = :coarse)
    ts = s.target
    V = Briareus.Vec3
    center(name) = (t = ts.target[]; es = elements(t, name);
                    sum(t.elem_centroid[e] * t.elem_area[e] for e in es) / sum(t.elem_area[es]))
    # a part in cylindrical coordinates: r, φ, z edit the same position
    @test add_part!(ts, :sphere; name = :post, x = 3.5, r1 = 0.2)
    i = part_index(ts, :post)
    @test set_part_coords!(ts, i, :cylindrical)
    @test Dict(part_position(ts.parts[][i]))[:pos_r] ≈ 3.5
    @test edit_part!(ts, i, :pos_phi, "90")
    @test center(:post) ≈ V(0, 3.5, 0) atol = 1e-9
    @test set_part_coords!(ts, i, :spherical) && edit_part!(ts, i, :pos_theta, "0")
    @test center(:post) ≈ V(0, 0, 3.5) atol = 1e-9
    @test !edit_part!(ts, i, :pos_r, "-1")
    # a part in another part's frame moves and turns with it
    @test add_part!(ts, :disk; name = :shield, anchor = :leh_top, z = 0.5, r1 = 2.0)
    @test add_part!(ts, :sphere; name = :marker, anchor = :shield, z = 1.0, r1 = 0.1)
    L = ts.params[][:length]
    @test center(:marker) ≈ V(0, 0, L / 2 + 1.5) atol = 1e-9
    @test edit_part!(ts, part_index(ts, :shield), :theta, "90")       # tilt the shield
    @test center(:marker) ≈ V(1.0, 0, L / 2 + 0.5) atol = 1e-9
    # renaming the frame's part keeps its children attached; loops are refused
    @test edit_part!(ts, part_index(ts, :shield), :name, "lid")
    @test ts.parts[][part_index(ts, :marker)].anchor === :lid
    @test !set_part_anchor!(ts, part_index(ts, :lid), :marker) && occursin("loop", ts.error[])
    @test !remove_part!(ts, part_index(ts, :lid))                     # its child needs it
    # backlighter-style aiming: beams on a component, following it when it moves
    idx = findall(==(:outer50), s.beams.cluster[])[1:2]
    aim_at_component!(s, idx, :lid)
    p0 = component_points(ts.target[])[:lid]
    @test all(i -> s.beams.beams[][i].aim ≈ p0, idx)
    edit_beam_set!(s.beams, idx, :aim_x, "0.3")                         # an offset from it
    @test edit_part!(ts, part_index(ts, :lid), :z, "1.5")              # move the lid up 1 mm
    p1 = component_points(ts.target[])[:lid]
    @test p1 - p0 ≈ V(0, 0, 1.0) atol = 1e-9
    @test all(i -> s.beams.beams[][i].aim ≈ V(0.3, p1[2], p1[3]), idx)
    @test s.beams.aim_ref[][idx] == [:lid, :lid] && count(!=(:none), s.beams.aim_ref[]) == 2
    f = joinpath(mktempdir(), "refs.briareus.toml")
    save_session(f, s)
    s2 = DesignSession(resolution = :coarse)
    load_session!(s2, f)
    @test s2.beams.aim_ref[] == s.beams.aim_ref[] && s2.beams.beams[] == s.beams.beams[]
    @test s2.target.parts[][part_index(s2.target, :post)].coords === :spherical
    aim_at_component!(s, idx, :none)
    @test all(==(:none), s.beams.aim_ref[])
    @test_throws ArgumentError aim_at_component!(s, idx, :nothing_here)
end

include(joinpath(pkgdir(Briareus), "test", "vrw_fixture.jl"))

@testset "VisRad workspaces in the window (no GL)" begin
    dir = mktempdir()
    src = joinpath(dir, "work.vrw")
    write(src, synthetic_workspace())
    s = DesignSession(resolution = :coarse)
    @test_throws ArgumentError export_visrad_workspace(s, joinpath(dir, "out.vrw"))
    import_visrad_workspace!(s, src)
    ts = s.target
    @test ts.kind[] === :custom && all(p -> p.kind === :visrad, ts.parts[])
    names = component_names(ts.target[])
    # the foil, the hohlraum's five components and the CAD tetrahedron take part; the
    # reference cone (not meshed) is skipped and the marker sphere, which VisRad leaves
    # out of its calculations, is a ghost
    @test :Rectangle_0 in names && :Cyl_Hohlraum_0_wall in names && :CAD_0 in names
    @test :Marker in all_component_names(ts) && !(:Marker in names)
    @test length(s.beams.beams[]) == 3 && s.beams.enabled[] == [true, true, false]
    @test s.beams.cluster[] == [:Ring_050, :Ring_050, :Ring_045]
    # the hohlraum keeps its place (x = -4 mm) and the CAD its winding repaired, outward
    t = ts.target[]
    w = elements(t, :Cyl_Hohlraum_0_wall)
    @test sum(t.elem_centroid[e] * t.elem_area[e] for e in w) / sum(t.elem_area[w]) ≈
          Briareus.Vec3(-4, 0, 0) atol = 1e-6
    # an edit to a beam goes back into a copy of the workspace, the rest unchanged
    edit_beams!(s.beams, 1, :power, "3.0"; by = :beam)
    out = joinpath(dir, "out.vrw")
    export_visrad_workspace(s, out)
    back = read_visrad_workspace(out)
    @test back.beams[1].power ≈ 3.0 rtol = 1e-4
    @test back.beams[3].power == 0
    @test count(k -> readlines(src)[k] != readlines(out)[k], eachindex(readlines(src))) == 3
    # sessions remember the template
    f = joinpath(dir, "vrw.briareus.toml")
    save_session(f, s)
    s2 = DesignSession(resolution = :coarse)
    load_session!(s2, f)
    @test s2.visrad_source[] == abspath(src)
    @test component_names(s2.target.target[]) == names
end

@testset "VisRad import in the background, roles, laser-only drive (no GL)" begin
    dir = mktempdir()
    src = joinpath(dir, "work.vrw")
    write(src, synthetic_workspace())
    # as in a window: the import runs on a worker and is applied from the queue
    s = DesignSession(resolution = :coarse)
    queue = Channel{Any}(Inf)
    s.deliver[] = f -> put!(queue, f)
    s.spawn[] = true
    before = s.target.target[]
    done = Ref(false)
    import_visrad_workspace!(s, src; done = () -> (done[] = true))
    @test !isempty(s.meshing[]) && s.target.target[] === before     # not yet replaced
    t0 = time()
    while !done[] && time() - t0 < 120
        f = take!(queue)
        f()
    end
    @test done[] && isempty(s.meshing[]) && startswith(s.status[], "imported work.vrw")
    @test s.target.kind[] === :custom && :CAD_0 in component_names(s.target.target[])
    @test length(s.beams.beams[]) == 3 && s.visrad_source[] == abspath(src)
    # a failed import says why and keeps the design
    bad = joinpath(dir, "empty.vrw")
    write(bad, "Workspace Format ID  = 11\n")
    t = s.target.target[]
    import_visrad_workspace!(s, bad)
    while isempty(s.status[]) || !isempty(s.meshing[]) || startswith(s.status[], "imported")
        take!(queue)()
    end
    @test occursin("import of empty.vrw failed", s.status[]) && s.target.target[] === t
    s.spawn[] = false
    s.deliver[] = f -> f()
    @test_throws ArgumentError import_visrad_workspace!(s, bad)              # inline: throws

    # roles: the imported rectangle becomes the backlighter foil, and back
    ts = s.target
    i = part_index(ts, :Rectangle_0)
    @test part_role(ts.parts[][i]) === :none
    @test set_part_role!(ts, i, :foil)
    @test ts.parts[][i].name === FOIL_NAME && FOIL_NAME in component_names(ts.target[])
    @test !set_part_role!(ts, part_index(ts, :Cyl_Hohlraum_0), :capsule)      # not a sphere
    @test occursin("is a sphere", ts.error[])
    # a CAD capsule: a sphere fitted to it takes its place in the calculations
    np = length(ts.parts[])
    @test set_part_role!(ts, part_index(ts, :CAD_0), :capsule)
    cap = ts.parts[][end]
    @test length(ts.parts[]) == np + 1 && cap.name === :capsule && cap.kind === :sphere
    @test :CAD_0 in ts.excluded[] && :capsule in component_names(ts.target[])
    tri = only(c for c in ts.ghost[].components if c.name === :CAD_0).mesh.vertices
    @test all(v -> abs(sqrt(sum(abs2, v - part_center(cap))) - cap.values[:r1]) < 1e-6, tri)   # its sphere
    # the excluded marker sphere can be the capsule (the fitted one gives the role up),
    # and stays left out
    m = part_index(ts, :Marker)
    @test set_part_role!(ts, m, :capsule)
    @test :capsule in ts.excluded[] && :capsule in all_component_names(ts) &&
          !(:capsule in component_names(ts.target[])) && ts.parts[][end].name === :sphere_1
    @test set_part_role!(ts, m, :none) && ts.parts[][m].name === :Marker && :Marker in ts.excluded[]
    remove_part!(ts, length(ts.parts[]))
    @test set_excluded!(ts, :CAD_0, false)
    # fit_sphere ignores a minority of points off the sphere (a fill tube)
    pts = [Briareus.Vec3(1.1sin(θ)cos(φ) + 0.2, 1.1sin(θ)sin(φ), 1.1cos(θ)) for θ in 0.1:0.2:3.1 for φ in 0:0.3:6.2]
    append!(pts, [Briareus.Vec3(x, 0.0, 0.0) for x in 1.6:0.1:4.0])
    c0, r0 = fit_sphere(pts)
    @test sqrt(sum(abs2, c0 - Briareus.Vec3(0.2, 0, 0))) < 1e-9 && r0 ≈ 1.1

    # a drive too large for the x-ray view factors computes the laser alone
    RADIATION_MAX_ELEMENTS[] = 50
    try
        @test !radiation_size(ts.target[]).feasible
        build_scenario!(s)
        dr = s.drive
        @test dr.scenario[] !== nothing && dr.scenario[].model === nothing
        @test occursin("laser only", dr.status[]) && occursin("most are in", dr.status[])
        res = dr.result[]
        @test res !== nothing && res.solution === nothing && res.spot !== nothing   # the foil's spot
        @test all(isnan, element_values(res, :radiosity))
        xray_image!(s; direction = (1, 0, 0), width = 4.0)
        @test s.drive.image[] === nothing && startswith(s.status[], "no x-ray image")
        set_objective!(s.optimize, :legendre)
        start_optimize!(s; spawn = false)
        @test occursin("laser only", s.optimize.status[]) && s.optimize.result[] === nothing
    finally
        RADIATION_MAX_ELEMENTS[] = 0
    end
    @test radiation_size(ts.target[]).feasible
    # the foil objective runs on the imported foil once it has its role
    build_scenario!(s)
    @test s.drive.scenario[].model !== nothing
    st = s.optimize
    set_objective!(st, :irradiance)
    st.fields[] = [:r]
    st.by[] = :beam
    st.iterations[] = 2
    start_optimize!(s; spawn = false)
    @test st.result[] !== nothing && st.result[].before.spot !== nothing

    # the target as STL: one named solid per component (ASCII), read back as parts
    f = joinpath(dir, "target.stl")
    export_target_stl(s, f; binary = false)
    solids = read_stl_solids(f)
    @test [Symbol(n) for (n, _) in solids] == component_names(ts.target[])
    @test sum(length(m.triangles) for (_, m) in solids) == length(ts.target[].tri_area)
    fb = joinpath(dir, "target_binary.stl")
    export_target_stl(s, fb)                                                  # binary
    @test length(read_stl(fb).triangles) == length(ts.target[].tri_area)
    # and back in as parts: one per solid of the ASCII file, one for the binary file
    ts2 = TargetSetup(kind = :custom)
    n0 = length(ts2.parts[])
    @test add_stl_parts!(ts2, f; h = 0.5)
    @test length(ts2.parts[]) - n0 == length(solids)
    @test all(p -> p.kind === :mesh && !isempty(p.object), ts2.parts[][n0+1:end])
    @test FOIL_NAME in component_names(ts2.target[])                          # named after its solid
    @test add_stl_parts!(ts2, fb) && isempty(ts2.parts[][end].object)
end

# A synthetic CAD workspace with a part of `n` sliver triangles (long thin strips of a
# tube of radius `radius` cm, as tessellated CAD writes thin shells), placed through a
# reference cone like a STEP import: mounting hardware when the tube is wide.
function sliver_workspace(n::Integer; radius = 0.3)
    txt = synthetic_workspace()
    nodes = NTuple{3,Float64}[]
    faces = NTuple{3,Int}[]
    m = n ÷ 2
    for k in 0:m-1          # two triangles per strip of a cylinder, z ∈ ±0.5 cm
        φ = 2π * k / m
        push!(nodes, (radius * cos(φ), radius * sin(φ), -0.5), (radius * cos(φ), radius * sin(φ), 0.5))
    end
    for k in 0:m-1
        a, b, c, d = 2k, 2k + 1, 2mod(k + 1, m), 2mod(k + 1, m) + 1
        push!(faces, (a, c, b), (b, c, d))
    end
    io = IOBuffer()
    println(io, "Begin Parameters for Object = BIG_CAD")
    println(io, "Object type                 = 13\n")
    # copy the CAD tetrahedron's parameter block, then the big mesh
    blk = match(r"Begin Parameters for Object = CAD_0\n.*?Number of Bounding Volumes = 0\n"s, txt).match
    body = split(blk, '\n'; limit = 3)[3]
    print(io, body)
    println(io, "\n   Gen. mesh object origin: x     = -0.1\n   Gen. mesh object origin: y     = 0\n   Gen. mesh object origin: z     = 0")
    println(io, "Number of vertices = ", length(nodes))
    println(io, "Number of surfaces = ", length(faces), "\n\nNode positions:")
    for (i, v) in enumerate(nodes)
        println(io, "  ", i - 1, "  ", join(v, "  "))
    end
    println(io, "Surface element connectivity:")
    for (i, f) in enumerate(faces)
        println(io, "  ", i - 1, "  ", join(f, "  "), "  -1")
    end
    println(io, "\nEnd Parameters for Object: BIG_CAD\n")
    big = String(take!(io))
    txt = replace(txt, "Number of View Factor Objects =  5" => "Number of View Factor Objects =  6")
    i = findfirst("Begin Parameters for Object = Marker", txt)
    return txt[1:first(i)-1] * big * txt[first(i):end]
end

@testset "a CAD workspace of 10⁵ triangles (no GL)" begin
    dir = mktempdir()
    src = joinpath(dir, "cad.vrw")
    # mounting hardware: 10⁵ slivers on a tube 60 mm across, far beyond the aim points
    write(src, sliver_workspace(100_000; radius = 3.0))
    s = DesignSession(resolution = :coarse)
    t1 = @elapsed import_visrad_workspace!(s, src)          # compiles too
    ts = s.target
    # drawn only: left out of the calculations, from its file's facets (no remeshing)
    @test :BIG_CAD in ts.excluded[] && !(:BIG_CAD in component_names(ts.target[]))
    @test length(ts.ghost[].tri_area) >= 100_000
    @test occursin("drawn only", s.status[]) && occursin("BIG_CAD", s.status[])
    @test t1 < 60
    # again, from the caches (a second import in the window): seconds
    s2 = DesignSession(resolution = :coarse)
    @test (@elapsed import_visrad_workspace!(s2, src)) < 10
    @test (@elapsed model_tree(s2)) < 0.5
    # VisRad's flags alone: in the calculations, too many elements for the x-ray drive
    small = joinpath(dir, "cad_small.vrw")                  # remeshing 10⁵ slivers takes a minute
    write(small, sliver_workspace(20_000))
    RADIATION_MAX_ELEMENTS[] = 500
    try
        s3 = DesignSession(resolution = :coarse)
        import_visrad_workspace!(s3, small; display_only = :none)
        @test :BIG_CAD in component_names(s3.target.target[])
        @test occursin("laser only", s3.status[]) && occursin("BIG_CAD", s3.status[])
        n = nelements(s3.target.target[])
        # left out again: cheap (its facets as they are)
        @test (@elapsed include_components!(s3, [:BIG_CAD], false)) < 10
        @test :BIG_CAD in s3.target.excluded[] && nelements(s3.target.target[]) < n
    finally
        RADIATION_MAX_ELEMENTS[] = 0
    end
end

@testset "retune to a reference drive (no GL)" begin
    s = DesignSession(resolution = :coarse)
    dir = mktempdir()
    beams = s.beams.beams[]
    pul = cone_fraction_pulses(beams, [0.0, 1.0, 3.0, 4.0], [0.0, 100.0, 300.0, 0.0],
                               [0.3, 0.3, 0.25, 0.25])
    table = write_pulse_table(joinpath(dir, "reference.csv"), pul)
    out = joinpath(dir, "retuned.csv")
    retune_to_reference!(s, table, out)                                 # no drive yet
    @test occursin("compute the drive first", s.drive.retune_status[])
    @test !isfile(out)
    build_scenario!(s)
    s.beams.enabled[] = [nif_bundle(b) != "31" for b in beams]          # bundle 31 down
    retune_to_reference!(s, table, out; cap = 1.5)
    @test !s.drive.retune_running[]
    @test startswith(s.drive.retune_status[], "46 quad powers")
    back = read_pulse_table(out)
    @test length(back) == 46
    @test !haskey(back, :Q31T)
    @test all(p -> iszero(p(0.0)) && iszero(p(4.0)), values(back))   # off where the reference is
    # the job: mode 1 of the dropped bundle removed, the drive kept, within the limits
    v = retune_pulses(s.drive.scenario[], beams, s.beams.enabled[], pul; cap = 1.5)
    ε1(h) = maximum(mode_amplitude(h.ylm[:, k], 1) for k in eachindex(h.time))
    @test ε1(v.before) > 0.01
    @test ε1(v.after) < ε1(v.before) / 10
    @test maximum(abs.(v.after.Tr .- v.reference.Tr)) < 0.5
    @test all(v.pulses[k](τ) <= 1.5 * pul[first(b.group for b in beams if Symbol(b.quad) == k)](τ) + 1e-9
              for k in keys(v.pulses) for τ in (1.0, 3.0))
    # per group, from the same table
    vg = retune_pulses(s.drive.scenario[], beams, s.beams.enabled[], pul; by = :group, cap = 1.5)
    @test Set(keys(vg.pulses)) == Set(keys(pul))
    @test_throws ArgumentError retune_pulses(s.drive.scenario[], beams, s.beams.enabled[],
                                             Dict(:x => PulseShape([0.0, 1.0], [1.0, 1.0])))
    @test_throws ArgumentError retune_pulses(s.drive.scenario[], beams, fill(false, length(beams)), pul)
end
