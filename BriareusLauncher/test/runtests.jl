# The launcher's shortcuts (written to temporary folders, never the real desktop or menu),
# the splash's quit protocol, and the icons it ships. The splash opens a window for about
# two seconds unless QT_QPA_PLATFORM=offscreen.
using Test
using BriareusLauncher
const L = BriareusLauncher

icons = L.icons_dir()
shim = "/home/x/.julia/bin/briareus"

@testset "BriareusLauncher" begin
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
    end

    @testset "freedesktop entry" begin
        mktempdir() do dir
            path = withenv("XDG_DATA_HOME" => dir) do
                only(L.desktop_entry(shim, joinpath(icons, "briareus.png")))
            end
            @test path == joinpath(dir, "applications", "briareus.desktop")
            text = read(path, String)
            @test occursin("Exec=\"$shim\" --from-shortcut\n", text) && occursin("Terminal=false", text)
        end
    end

    @testset "macOS app bundle" begin
        mktempdir() do dir
            app = only(L.mac_app(shim, joinpath(icons, "briareus.icns"); applications = dir))
            @test app == joinpath(dir, "Briareus.app")
            contents = joinpath(app, "Contents")
            plist = read(joinpath(contents, "Info.plist"), String)
            @test occursin("<key>CFBundleExecutable</key><string>briareus</string>", plist)
            @test read(joinpath(contents, "MacOS", "briareus"), String) ==
                  "#!/bin/sh\nexec $shim --from-shortcut\n"
            @test read(joinpath(contents, "Resources", "briareus.icns")) ==
                  read(joinpath(icons, "briareus.icns"))
        end
    end

    if Sys.iswindows()
        @testset "Windows shortcut" begin
            mktempdir() do dir
                folder = "'" * replace(dir, "'" => "''") * "'"
                made = L.windows_shortcuts("C:\\x\\briareus.bat", joinpath(icons, "briareus.ico");
                                           folders = folder)
                @test made == [joinpath(dir, "Briareus.lnk")] && isfile(only(made))
            end
        end
    end

    @testset "splash quits on \"quit\"" begin
        s = L.start_splash("Testing the splash")
        if s.process !== nothing
            sleep(2)
            alive = L.alive(s)
            L.close_splash(s)
            t = time()
            while process_running(s.process) && time() - t < 10
                sleep(0.1)
            end
            @test alive && !process_running(s.process) && s.process.exitcode == 0
        end
        rm(s.status; force = true)
    end
end
