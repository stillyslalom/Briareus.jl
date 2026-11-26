# The `briareus` command's shortcuts and their script (written to temporary folders, never
# the real desktop or menu), the splash's quit protocol, and the icons it ships. The
# splash opens a window for a few seconds unless QT_QPA_PLATFORM=offscreen.
const A = BriareusGUI

@testset "briareus command" begin
    icons = A.icons_dir()
    shim = "/home/x/.julia/bin/briareus"
    log = "/home/x/.local/state/briareus/briareus.log"

    @testset "icons" begin
        @test all(isfile, joinpath.(icons, ("briareus.ico", "briareus.png", "briareus.icns")))
        # ICNS uses big-endian lengths, with the file length in its header and each entry
        # length including the entry header.
        b = read(joinpath(icons, "briareus.icns"))
        len(i) = Int(ntoh(reinterpret(UInt32, b[i:i+3])[1]))
        @test String(b[1:4]) == "icns" && len(5) == length(b)
        offsets = Int[]
        i = 9
        while i < length(b)
            push!(offsets, i)
            i += len(i + 4)
        end
        @test i == length(b) + 1
        @test Set(String(b[o:o+3]) for o in offsets) ==
              Set(["icp4", "icp5", "icp6", "ic07", "ic08", "ic09", "ic10"])
        @test all(b[o+8:o+11] == [0x89, UInt8('P'), UInt8('N'), UInt8('G')] for o in offsets)
        @test isfile(A.SPLASH_QML)
    end

    @testset "the session briareus() starts without one" begin
        # the keywords briareus() passes on construct a session (a renamed keyword once
        # broke the launch with no session while every test built a separate session)
        s = A._launch_session(; starter = :hohlraum, resolution = :coarse, preset = nothing, beams = nothing)
        @test s isa DesignSession && :hohlraum in [p.name for p in s.target.parts[]]
        s = A._launch_session(; starter = :none, resolution = :coarse, preset = :nif_quads, beams = nothing)
        @test s isa DesignSession && isempty(s.target.parts[]) && length(s.beams.beams[]) == 48
    end

    @testset "the app's entry point" begin
        @static if isdefined(Base, Symbol("@main"))
            @test isdefined(A, :main)
        end
        @test A._main(["--from-shortcut", "/no/such/design.briar.toml"]) == 1
    end

    @testset "shortcut scripts" begin
        sh = A.posix_script(shim, log; tool = "/qt/bin/qml")
        @test startswith(sh, "#!/bin/sh\n") && occursin("$(Base.shell_escape(shim)) --from-shortcut >\"\$log\" 2>&1\n", sh)
        @test occursin("export BRIAREUS_SPLASH\n", sh) && occursin("exec /qt/bin/qml $(Base.shell_escape(A.SPLASH_QML)) -- ", sh)
        bare = A.posix_script(shim, log; tool = nothing)
        @test !occursin("qml", bare) && occursin("$(Base.shell_escape(shim)) --from-shortcut", bare)
        cmd = A.windows_script("C:\\x\\briareus.bat", "C:\\x\\Briareus\\briareus.log"; tool = "C:\\qt\\qml.exe")
        @test occursin("\r\ncall \"C:\\x\\briareus.bat\" --from-shortcut >\"%log%\" 2>&1\r\n", cmd)
        @test occursin("start \"\" /b \"C:\\qt\\qml.exe\" ", cmd) && !occursin(r"[^\r]\n", cmd)
    end

    @testset "freedesktop entry" begin
        mktempdir() do dir
            path = withenv("XDG_DATA_HOME" => dir) do
                only(A.desktop_entry("/x/briareus-shortcut.sh", joinpath(icons, "briareus.png")))
            end
            @test path == joinpath(dir, "applications", "briareus.desktop")
            text = read(path, String)
            @test occursin("Exec=\"/x/briareus-shortcut.sh\"\n", text) && occursin("Terminal=false", text)
        end
    end

    @testset "macOS app bundle" begin
        mktempdir() do dir
            script = A.posix_script(shim, log; tool = nothing)
            app = only(A.mac_app(script, joinpath(icons, "briareus.icns"); applications = dir))
            @test app == joinpath(dir, "Briareus.app")
            contents = joinpath(app, "Contents")
            plist = read(joinpath(contents, "Info.plist"), String)
            @test occursin("<key>CFBundleExecutable</key><string>briareus</string>", plist)
            @test read(joinpath(contents, "MacOS", "briareus"), String) == script
            @test read(joinpath(contents, "Resources", "briareus.icns")) ==
                  read(joinpath(icons, "briareus.icns"))
        end
    end

    if Sys.iswindows()
        @testset "Windows shortcut" begin
            mktempdir() do dir
                folder = "'" * replace(dir, "'" => "''") * "'"
                made = A.windows_shortcuts("C:\\x\\briareus-shortcut.cmd", joinpath(icons, "briareus.ico");
                                           folders = folder)
                @test made == [joinpath(dir, "Briareus.lnk")] && isfile(only(made))
            end
        end
    end

    @testset "splash quits on \"quit\"" begin
        if A.qml_tool() !== nothing
            status = tempname() * ".txt"
            write(status, "Testing the splash")
            p = run(pipeline(A.splash_cmd(status); stdout = devnull, stderr = devnull); wait = false)
            sleep(2)
            alive = process_running(p)
            write(status, "quit")
            t = time()
            while process_running(p) && time() - t < 10
                sleep(0.1)
            end
            @test alive && !process_running(p) && p.exitcode == 0
            rm(status; force = true)
        end
    end

    # The script itself runs with a stand-in for the shim of the app. The stand-in reads the
    # status file of the splash from BRIAREUS_SPLASH, sends "quit" as the window would, and
    # exits with the code given. The script passes the code on and removes the status file.
    @testset "shortcut script runs the app" begin
        tool = A.qml_tool()
        mktempdir() do dir
            log = joinpath(dir, "logs", "briareus.log")
            runs = Pair{Int,Bool}[]
            for code in (0, 3), splash in (true, false)
                splash || code == 0 || continue      # a failure without a splash opens the log
                if Sys.iswindows()
                    fake = joinpath(dir, "fake.bat")
                    write(fake, "@echo off\r\necho args %*\r\n" *
                                "if defined BRIAREUS_SPLASH (>\"%BRIAREUS_SPLASH%\" echo quit)\r\n" *
                                "ping -n 2 127.0.0.1 >nul\r\nexit /b $code\r\n")
                    script = joinpath(dir, "briareus-shortcut.cmd")
                    write(script, A.windows_script(fake, log; tool = splash ? tool : nothing))
                    cmd = `cmd /c $script`
                else
                    fake = joinpath(dir, "fake.sh")
                    write(fake, "#!/bin/sh\necho args \"\$@\"\n" *
                                "[ -z \"\$BRIAREUS_SPLASH\" ] || echo quit >\"\$BRIAREUS_SPLASH\"\n" *
                                "sleep 1\nexit $code\n")
                    chmod(fake, 0o755)
                    script = joinpath(dir, "briareus-shortcut.sh")
                    write(script, A.posix_script(fake, log; tool = splash ? tool : nothing))
                    chmod(script, 0o755)
                    cmd = `$script`
                end
                tmp = mkpath(joinpath(dir, "tmp"))
                p = withenv("TEMP" => tmp, "TMPDIR" => tmp) do
                    run(pipeline(ignorestatus(cmd); stdout = devnull, stderr = devnull))
                end
                ok = p.exitcode == code && isempty(readdir(tmp)) &&
                     strip(read(log, String)) == "args --from-shortcut"
                push!(runs, code => ok)
            end
            @test length(runs) == 3 && all(last, runs)
        end
    end
end

# The Help topics (qml/common/help.js) against the QML that shows them: every topic an ⓘ,
# a section heading, the inspector's title, a menu or the Help window names exists; every
# topic is shown by one of them, but for those of the first page of the Help window
# (Getting around), which its page list reaches; each has a one-line hint for the hover
# and a text in short paragraphs; and the topics of the survey's pages are there.
@testset "Help topics" begin
    qml = A.QML_DIR
    js = read(joinpath(qml, "common", "help.js"), String)
    keys_ = Set(m[1] for m in eachmatch(r"key: \"([a-z_]+\.[a-z_]+)\"", js))
    refs = Set{String}()
    for (dir, _, files) in walkdir(qml), f in files
        endswith(f, ".qml") || continue
        src = read(joinpath(dir, f), String)
        for m in eachmatch(r"\b(?:topic|info):\s*([^;\n}]*)", src), k in eachmatch(r"\"([a-z_]+\.[a-z_]+)\"", m[1])
            push!(refs, k[1])
        end
        # a menu or a command that opens the Help window at a topic
        union!(refs, (m[1] for m in eachmatch(r"(?:openAt|br_help)\(\"([a-z_]+\.[a-z_]+)\"\)", src)))
        # the inspector's title ⓘ: its `intro` (a task's is "<task>.intro")
        if f == "Inspector.qml"
            intro = match(r"property string intro: \{(.*?)\n    \}"s, src)[1]
            union!(refs, (m[1] for m in eachmatch(r"return \"([a-z_]+\.[a-z_]+)\"", intro)))
            union!(refs, (t * ".intro" for t in ("drive", "optimize", "clearance", "sweep", "talis")))
        end
    end
    # the Time tab's ⓘ: the chart shown ("time." + app.timeChart; the diagnostics chart is
    # the Diagnostics tab's)
    union!(refs, ("time." * String(c) for c in A.Controllers.TIME_CHARTS if c !== :diagnostics))
    @test isempty(setdiff(refs, keys_))                 # every topic shown exists
    @test all(startswith("start."), setdiff(keys_, refs))   # every topic is shown somewhere
    # Each topic has a hint of one short line and a text without a long unbroken paragraph.
    topics = collect(eachmatch(r"key: \"([a-z_]+\.[a-z_]+)\", title: \"[^\"]*\",\s*tip: \"([^\"]*)\",\s*text: md\(`(.*?)`\)"s, js))
    @test length(topics) == length(keys_)
    @test all(m -> 0 < length(m[2]) <= 64, topics)
    paragraphs(t) = filter(!isempty, strip.(split(t, r"\n\s*\n")))
    @test all(m -> all(p -> length(p) < 560, paragraphs(m[3])), topics)
    survey = ["time.pulses", "time.drive", "diagnostics.chart", "time.energy", "time.retune", "diagnostics.images",
              "diagnostics.record", "diagnostics.config", "diagnostics.folders", "part.windows", "beams.intro", "beams.nodes",
              "beams.table", "optimize.results", "clearance.results", "view.cutaway", "view.display",
              "part.leh_closure"]
    @test issubset(survey, keys_)
end

@testset "the first-frame precompile list (no GL)" begin
    list = joinpath(pkgdir(BriareusGUI), "src", "qt", "first_frame.jl")
    lines = filter(l -> startswith(l, "precompile("), readlines(list))
    @test !isempty(lines) && issorted(lines) && allunique(lines)
    # no anonymous closure named by its numbers, no closure of the GUI, no function of the GUI
    # as the callee (those are explicit statements in the workload)
    @test !any(l -> occursin(r"#\d+#\d+\"", l), lines)
    @test !any(l -> occursin(r"BriareusGUI(\.\w+)*\.var\"#", l), lines)
    @test !any(l -> startswith(l, "precompile(Tuple{typeof(BriareusGUI."), lines)
    # the recorder's filter keeps the whole list (in the repository, where dev/ exists)
    filter_file = joinpath(pkgdir(BriareusGUI), "..", "dev", "gui_first_frame_filter.jl")
    if isfile(filter_file)
        m = Module()
        Base.include(m, filter_file)
        @test Base.invokelatest(m.filter_signatures, lines) == lines
        keep(l) = Base.invokelatest(m.keep_signature, l)
        @test !keep("precompile(Tuple{BriareusGUI.var\"#95#96\"{BriareusGUI.BriareusShell}, Bool})")
        @test !keep("precompile(Tuple{GLMakie.var\"#_#81#83\", Base.Dict{String, String}})")
        @test !keep("precompile(Tuple{typeof(Base._all), BriareusGUI.var\"#f##0#f##1\"{Int64}, Base.Colon})")
        @test !keep("precompile(Tuple{typeof(BriareusGUI.briareus)})")
        @test !keep("# a comment")
        @test keep("precompile(Tuple{typeof(QML.rowcount), QML.ItemModelData{Array{BriareusGUI.TreeItem, 1}}})")
        @test keep("precompile(Tuple{GLMakie.var\"#add_scene!##2#add_scene!##3\"{GLMakie.Screen{QMLMakie.QMLWindow}}, Bool})")
    end
    # the replay counts what compiles and what does not
    f = joinpath(mktempdir(), "list.jl")
    write(f, "# header\nprecompile(Tuple{typeof(Base.sin), Float64})\nprecompile(Tuple{typeof(NoSuchModule.f), Int64})\n" *
             "precompile(Tuple{typeof(Base.no_such_function_here), Int64})\n")
    @test BriareusGUI._replay_precompiles(f) == (; compiled = 1, failed = 2)
    @test BriareusGUI._replay_precompiles(joinpath(dirname(f), "missing.jl")) == (; compiled = 0, failed = 0)

    @testset "arguments" begin
        c = A.parse_command(["design.briar", "--tab", "talis"])
        @test c.command == "" && c.files == ["design.briar"] && c.tab == "talis" && !c.json
        c = A.parse_command(["review-1w", "in.vrw", "--out", "dir", "--threshold", "0.5", "--json"])
        @test (c.command, c.files, c.out, c.threshold, c.json) == ("review-1w", ["in.vrw"], "dir", 0.5, true)
        @test A.parse_command(["review-1w", "in.vrw", "--out", "dir", "--images"]).images && !c.images
        @test A.parse_command(["diff", "a.briar", "b.vrw", "--out", "d"]).files == ["a.briar", "b.vrw"]
        c = A.parse_command(["run", "s.jl", "--out", "x", "talis"])
        @test c.command == "run" && c.files == ["s.jl"] && c.rest == ["--out", "x", "talis"] && c.out === nothing
        @test A.parse_command(String[]).files == String[] && A.parse_command(["--from-shortcut"]).files == String[]
        @test_throws ArgumentError A.parse_command(["talis", "a.briar", "--out"])
        @test_throws ArgumentError A.parse_command(["review-1w", "a.briar", "--threshold", "wide"])
        @test A._main(["--help"]) == 0 && all(t -> occursin(t, A.USAGE), A.TABS)
        mktempdir() do dir
            f = joinpath(dir, "d.briar.toml")
            save_session(f, DesignSession(resolution = :coarse))
            @test A._main([f, "--tab", "nonesuch"]) == 1            # checked before a window opens
            @test A._main(["talis", f]) == 1                         # no --out
            @test A._main(["roe", joinpath(dir, "absent.briar"), "--out", dir]) == 1
        end
    end

    @testset "headless commands" begin
        isdefined(@__MODULE__, :synthetic_workspace) || include(joinpath(pkgdir(Briareus), "test", "vrw_fixture.jl"))
        mktempdir() do dir
            # a saved session of the hohlraum starter: the commands against the core calls
            src = joinpath(dir, "design.briar.toml")
            s0 = DesignSession(resolution = :coarse)
            save_session(src, s0)
            sess = A.load_design(src)
            @test sess isa DesignSession && length(sess.beams.beams[]) == length(s0.beams.beams[])
            @test_throws ArgumentError A.load_design(joinpath(dir, "absent.briar"))
            t, bs, kw = sess.target.target[], active_beams(sess.beams), talis_kwargs(sess.talis, sess.review)
            ref = talis_check(t, bs; kw...)
            out = A.talis_command(sess, joinpath(dir, "out", "talis.csv"); json = true)
            @test basename.(out.paths) == ["talis.csv", "talis.json"] && all(isfile, out.paths)
            @test [(x.beam, x.rule, x.status, x.margin) for x in out.report.rows] ==
                  [(x.beam, x.rule, x.status, x.margin) for x in ref.rows]
            @test length(Briareus.JSON.parsefile(out.paths[2])["rows"]) == length(ref.rows)
            rv = A.review_1w_command(sess, joinpath(dir, "review"); threshold = 2.0, json = true)
            @test rv.review.threshold == 2.0 && length(rv.review.report.rows) == length(ref.rows)
            @test basename(rv.paths[1]) == "review-1w.html" && basename(rv.paths[end]) == "review-1w.json" &&
                  count(p -> endswith(p, ".svg"), rv.paths) == length(rv.review.items)
            html = read(rv.paths[1], String)
            @test occursin("<h2>All checks</h2>", html) && count("<figure id=\"approach-", html) == length(rv.review.items)
            ro = A.roe_command(sess, joinpath(dir, "roe"); json = true)
            @test [c.name for c in ro.criteria] == [:drive, :P2, :P4, :hotspot_velocity]
            @test basename.(ro.paths) == ["roe.html", "roe.csv", "roe.md", "ports.svg", "roe.json"]
            sc = sess.drive.scenario[]
            direct = roe_scan(sc, display_beams(sess.beams), ro.criteria)
            @test ro.report.tiers == direct.tiers && [o.scores for o in ro.report.singles] == [o.scores for o in direct.singles]
            @test occursin("close on themselves", ro.notes[1]) && occursin(ro.notes[1], read(ro.paths[1], String))

            # through the command line: messages, files and exit codes
            code = redirect_stdout(devnull) do
                A._main(["talis", src, "--out", joinpath(dir, "cli", "t.csv"), "--json"])
            end
            @test code == (any(x -> x.status === :fail, ref.rows) ? 2 : 0) && isfile(joinpath(dir, "cli", "t.json"))
            @test redirect_stdout(() -> A._main(["review-1w", src, "--out", joinpath(dir, "cli1w")]), devnull) == 0 &&
                  isfile(joinpath(dir, "cli1w", "review-1w.html")) && !isfile(joinpath(dir, "cli1w", "review-1w.json"))
            @test redirect_stdout(() -> A._main(["roe", src, "--out", joinpath(dir, "cliroe")]), devnull) == 0 &&
                  isfile(joinpath(dir, "cliroe", "roe.html"))

            # two designs: the second with every aim moved 0.2 mm along z
            s1 = DesignSession(resolution = :coarse)
            moved = [remake(b; aim = b.aim + Briareus.Vec3(0.0, 0.0, 0.2)) for b in s1.beams.beams[]]
            set_beams!(s1.beams, moved)
            src1 = joinpath(dir, "moved.briar.toml")
            save_session(src1, s1)
            df = A.diff_command(A.load_design(src), A.load_design(src1), joinpath(dir, "diff"); json = true)
            @test basename.(df.paths) == ["diff.html", "diff.json"] && all(isfile, df.paths)
            @test all(x -> x.status === :same, df.diff.parts)
            @test all(x -> x.changes == [:aim] && abs(x.aim_shift[3] - 0.2) < 1e-9 && abs(x.distance - 0.2) < 1e-9, df.diff.beams)
            @test !isempty(df.diff.talis) && !isempty(df.diff.physics) && isempty(df.notes)
            same = A.diff_command(A.load_design(src), A.load_design(src), joinpath(dir, "diff0"))
            @test all(x -> x.change === :same, same.diff.talis) && all(x -> x.delta == 0, same.diff.physics)
            @test redirect_stdout(() -> A._main(["diff", src, src1, "--out", joinpath(dir, "clidiff")]), devnull) == 0 &&
                  isfile(joinpath(dir, "clidiff", "diff.html"))
            @test A._main(["diff", src, "--out", joinpath(dir, "clidiff")]) == 1

            # the design as a script: run, it gives the target, beams, check and drive of the session
            sp = write_session_script(joinpath(dir, "design.jl"), s1)
            text = read(sp, String)
            @test session_script(s1; title = "design.jl") == text && !occursin("\"window\"", text)
            ms = A.run_script(sp)
            m = (; (n => Base.invokelatest(getglobal, ms, n) for n in (:target, :beams, :talis, :evaluation, :session))...)
            await_target!(s1.target)
            t1 = s1.target.target[]
            @test m.target.elem_area == t1.elem_area && [c.name for c in m.target.components] == [c.name for c in t1.components]
            @test [(b.name, b.aim, b.power, b.spot) for b in m.beams] == [(b.name, b.aim, b.power, b.spot) for b in active_beams(s1.beams)]
            ref1 = talis_check(t1, active_beams(s1.beams); talis_kwargs(s1.talis, s1.review)...)
            @test [(x.beam, x.rule, x.status, x.margin) for x in m.talis.rows] == [(x.beam, x.rule, x.status, x.margin) for x in ref1.rows]
            s1.spawn[] = false
            build_scenario!(s1)
            @test mean_flux(m.evaluation) ≈ mean_flux(evaluate(s1.drive.scenario[], display_beams(s1.beams))) rtol = 1e-12
            @test session_state(m.session)["beams"] == session_state(s1)["beams"]

            # the optimization saved in a design, run without a window for a few evaluations:
            # the optimized session, its SST pointing and the summary carry the optimum
            so = DesignSession(resolution = :coarse)
            so.optimize.max_evaluations[] = 3
            osrc = joinpath(dir, "opt.briar.toml")
            save_session(osrc, so)
            os = A.load_design(osrc)
            start = copy(os.beams.beams[])
            oo = A.optimize_command(os, joinpath(dir, "opt"); json = true)
            @test basename.(oo.paths) == ["optimized.briar", "pointing.csv", "optimize.json"] && all(isfile, oo.paths)
            @test oo.result.applied && os.optimize.max_evaluations[] == 3 && oo.result.value <= oo.result.initial_value
            @test os.beams.beams[][oo.result.indices] == oo.result.beams && any(i -> os.beams.beams[][i].aim != start[i].aim, oo.result.indices)
            back = A.load_design(oo.paths[1])
            @test [b.aim for b in back.beams.beams[]] == [b.aim for b in os.beams.beams[]]
            j = Briareus.JSON.parsefile(oo.paths[3])
            @test j["value"] == oo.result.value && j["objective"] == String(os.optimize.objective[]) &&
                  [Briareus.Vec3(x["aim_mm"]...) for x in j["beams"]] == [b.aim for b in oo.result.beams]
            @test redirect_stdout(() -> A._main(["optimize", osrc, "--out", joinpath(dir, "cliopt")]), devnull) == 0 &&
                  isfile(joinpath(dir, "cliopt", "optimized.briar")) && !isfile(joinpath(dir, "cliopt", "optimize.json"))
            @test A.parse_command(["optimize", osrc, "--out", "o"]).command == "optimize"

            # a foil design with no hohlraum: criteria from the foil, on OMEGA beams, where
            # the TaLIS rules do not apply
            fsrc = joinpath(dir, "foil.briar.toml")
            save_session(fsrc, DesignSession(starter = :foil, resolution = :coarse))
            fs = A.load_design(fsrc)
            fo = A.roe_command(fs, joinpath(dir, "foilroe"))
            @test !isempty(fo.criteria) && all(c -> occursin("irradiance", string(c.name)) || occursin("centroid", string(c.name)), fo.criteria)
            @test isfile(joinpath(dir, "foilroe", "roe.html")) && length(fo.report.units) > 0
            @test_throws ArgumentError A.talis_command(fs, joinpath(dir, "foil.csv"))
            @test A._main(["talis", fsrc, "--out", joinpath(dir, "foil.csv")]) == 1

            # a VisRad workspace mock-up
            vrw = joinpath(dir, "work.vrw")
            write(vrw, synthetic_workspace())
            vs = A.load_design(vrw)
            @test vs.visrad_source[] == vrw || !isempty(vs.beams.beams[])
            vo = A.roe_command(vs, joinpath(dir, "vrwroe"))
            @test !isempty(vo.criteria) && isfile(joinpath(dir, "vrwroe", "roe.html"))
            # an optimization that cannot run gives the reason of the Optimize page
            msg = try
                A.optimize_command(vs, joinpath(dir, "vrwopt"))
                ""
            catch err
                err isa ArgumentError || rethrow()
                err.msg
            end
            @test startswith(msg, "the optimization did not run: capsule symmetry needs a hohlraum with a capsule")

            # a script: Briareus and the controllers in scope, its arguments in ARGS
            script = joinpath(dir, "s.jl")
            write(script, "s = DesignSession(resolution = :coarse)\nwrite(ARGS[1], string(length(nif_quads()), \" \", length(s.beams.beams[])))\n")
            before = copy(ARGS)
            @test A._main(["run", script, joinpath(dir, "s.txt")]) == 0 && read(joinpath(dir, "s.txt"), String) == "48 48"
            @test ARGS == before && A._main(["run", joinpath(dir, "absent.jl")]) == 1
        end
    end
end
