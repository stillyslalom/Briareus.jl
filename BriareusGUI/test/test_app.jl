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
        # ICNS: big-endian lengths, the file's in its header, each entry's own included
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

    # The script itself, with a stand-in for the app's shim: it reads the splash's status
    # file from BRIAREUS_SPLASH, sends "quit" as the window would, and exits with the code
    # given; the script passes the code on and removes the status file.
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
