"""
    BriareusLauncher

The `briareus` command (the `[apps]` entry of this package's `Project.toml`, installed as
the BriareusGUI README describes). Loading BriareusGUI takes 15-20 s, most of it loading
GLMakie's and QML's package images, so the command first shows a splash with Qt's own
`qml` tool (within about a second) and closes it when the window shows.

    briareus                       # the window on a fresh design
    briareus design.briareus.toml  # on a saved session
    briareus --shortcut            # add Briareus to the desktop and the applications menu

A shortcut runs `briareus --from-shortcut`, which writes the output to [`log_path`](@ref)
and shows a failure to start in the splash.
"""
module BriareusLauncher

import Qt6Base_jll, Qt6Declarative_jll
using Base64: base64encode

const GUI = Base.PkgId(Base.UUID("27459d3e-0914-41cf-82e2-9495bb7270ba"), "BriareusGUI")
const SPLASH_QML = joinpath(@__DIR__, "splash.qml")

# BriareusGUI's icons, found without loading it
icons_dir() = joinpath(dirname(Base.locate_package(GUI)), "qml", "icons")

"""
    log_path() -> String

Where a launch from a shortcut writes its output: `%LOCALAPPDATA%\\Briareus` on Windows,
`~/Library/Logs/Briareus` on macOS, `\$XDG_STATE_HOME/briareus` elsewhere.
"""
function log_path()
    dir = Sys.iswindows() ? joinpath(get(ENV, "LOCALAPPDATA", homedir()), "Briareus") :
          Sys.isapple() ? joinpath(homedir(), "Library", "Logs", "Briareus") :
          joinpath(get(ENV, "XDG_STATE_HOME", joinpath(homedir(), ".local", "state")), "briareus")
    return joinpath(dir, "briareus.log")
end

# ---------------------------------------------------------------- splash

struct Splash
    process::Union{Nothing,Base.Process}       # nothing: no `qml` tool on this platform
    status::String                             # the file the splash polls
end

function qml_tool()
    bin = joinpath(Qt6Declarative_jll.artifact_dir, "bin")
    for p in (joinpath(bin, "qml.exe"), joinpath(bin, "qml"),
              joinpath(bin, "qml.app", "Contents", "MacOS", "qml"))
        isfile(p) && return p
    end
    return nothing
end

file_url(path) = "file:///" * lstrip(replace(path, '\\' => '/'), '/')

function start_splash(stage::AbstractString)
    status = tempname() * ".txt"
    write(status, stage)
    tool = qml_tool()
    tool === nothing && return Splash(nothing, status)
    # the Qt libraries as a JLL executable finds them
    var = Sys.iswindows() ? "PATH" : Sys.isapple() ? "DYLD_FALLBACK_LIBRARY_PATH" : "LD_LIBRARY_PATH"
    dirs = vcat(Qt6Declarative_jll.LIBPATH_list,
                joinpath(Qt6Declarative_jll.artifact_dir, Sys.iswindows() ? "bin" : "lib"))
    if Sys.isapple()
        # Qt's frameworks: the tool's rpaths name each framework folder of its own artifact
        # (`@loader_path/../lib/QtCore.framework`), so dyld's fallback needs those folders
        frameworks = [joinpath(d, f) for d in dirs if isdir(d) for f in readdir(d) if endswith(f, ".framework")]
        append!(dirs, frameworks)
    end
    haskey(ENV, var) && push!(dirs, ENV[var])
    env = (var => join(dirs, Sys.iswindows() ? ';' : ':'),
           "QT_PLUGIN_PATH" => joinpath(Qt6Base_jll.artifact_dir, "plugins"),
           "QML_IMPORT_PATH" => joinpath(Qt6Declarative_jll.artifact_dir, "qml"),
           "QML_XHR_ALLOW_FILE_READ" => "1")
    icon = joinpath(icons_dir(), "briareus.png")
    cmd = `$tool $SPLASH_QML -- $(file_url(icon)) $(file_url(status))`
    process = try
        run(pipeline(addenv(cmd, env...); stdout = devnull, stderr = devnull); wait = false)
    catch err
        @warn "no splash" exception = err
        nothing
    end
    return Splash(process, status)
end

# A partial write reads as empty, which the splash tolerates for a second.
set_stage(s::Splash, text) = try write(s.status, text) catch end
close_splash(s::Splash) = set_stage(s, "quit")
alive(s::Splash) = s.process !== nothing && process_running(s.process)

# ---------------------------------------------------------------- the command

function (@main)(args)
    "--shortcut" in args && (foreach(p -> println("created ", p), install_shortcut()); return 0)
    detached = "--from-shortcut" in args
    files = filter(a -> !startswith(a, "--"), args)
    log = nothing
    if detached
        mkpath(dirname(log_path()))
        log = open(log_path(), "w")
        redirect_stdio(stdout = log, stderr = log)
    end
    if !isempty(files) && !isfile(files[1])
        println(stderr, "briareus: no such file: ", files[1])
        return 1
    end
    splash = start_splash("Loading the graphics libraries")
    code = 0
    try
        # by name, as `using` loads it: `Base.require(GUI)` (the PkgId) took 36 s, not 16
        gui = Base.require(@__MODULE__, :BriareusGUI)
        set_stage(splash, "Opening the window")
        Base.invokelatest(open_window, gui, files, splash)
    catch err
        code = 1
        @error "Briareus could not start" exception = (err, catch_backtrace())
        if detached
            msg = first(split(sprint(showerror, err), '\n'))
            if alive(splash)
                set_stage(splash, "error: Briareus could not start: $msg. The log is $(log_path()).")
                wait(splash.process)
            else
                open_file(log_path())
            end
        end
    finally
        close_splash(splash)
        sleep(0.5)                         # the splash reads "quit" before the file goes
        rm(splash.status; force = true)
        log === nothing || flush(log)
    end
    return code
end

function open_window(gui, files, splash)
    sess = gui.DesignSession()
    isempty(files) || gui.load_session!(sess, files[1])
    gui.briareus(sess; shown = () -> close_splash(splash))
    return
end

function open_file(path)
    cmd = Sys.iswindows() ? `cmd /c start "" $path` : Sys.isapple() ? `open $path` : `xdg-open $path`
    try run(cmd; wait = false) catch end
    return
end

# ---------------------------------------------------------------- shortcuts

"""
    install_shortcut() -> Vector{String}

Add Briareus to the desktop and the applications menu, as shortcuts that run the
installed `briareus` command with `--from-shortcut`, and return what was created:
on Windows, shortcuts on the desktop and in the Start menu (through `conhost --headless`,
so no console window opens); on macOS, `~/Applications/Briareus.app`; elsewhere, a
freedesktop menu entry. Run again after moving the depot or the clone.
"""
function install_shortcut()
    shim = joinpath(first(DEPOT_PATH), "bin", Sys.iswindows() ? "briareus.bat" : "briareus")
    isfile(shim) || error("no briareus command in $(dirname(shim)): install it first (BriareusGUI README)")
    icons = icons_dir()
    return Sys.iswindows() ? windows_shortcuts(shim, joinpath(icons, "briareus.ico")) :
           Sys.isapple() ? mac_app(shim, joinpath(icons, "briareus.icns")) :
           desktop_entry(shim, joinpath(icons, "briareus.png"))
end

# `folders`: a PowerShell expression for the folders that get a shortcut
function windows_shortcuts(shim, icon; folders = "[Environment]::GetFolderPath('Desktop'), " *
                                                 "[Environment]::GetFolderPath('Programs')")
    q(s) = "'" * replace(s, "'" => "''") * "'"            # a PowerShell string literal
    conhost = joinpath(get(ENV, "SystemRoot", "C:\\Windows"), "System32", "conhost.exe")
    script = """
        \$ProgressPreference = 'SilentlyContinue'
        \$shell = New-Object -ComObject WScript.Shell
        foreach (\$dir in $folders) {
            \$path = Join-Path \$dir 'Briareus.lnk'
            \$link = \$shell.CreateShortcut(\$path)
            \$link.TargetPath = $(q(conhost))
            \$link.Arguments = $(q("--headless \"$shim\" --from-shortcut"))
            \$link.IconLocation = $(q(icon * ",0"))
            \$link.WorkingDirectory = [Environment]::GetFolderPath('MyDocuments')
            \$link.Description = 'Briareus laser-target design'
            \$link.Save()
            \$path
        }
        """
    encoded = base64encode(reinterpret(UInt8, transcode(UInt16, script)))
    return readlines(`powershell -NoProfile -NonInteractive -EncodedCommand $encoded`)
end

function mac_app(shim, icns; applications = joinpath(homedir(), "Applications"))
    app = joinpath(applications, "Briareus.app")
    contents = joinpath(app, "Contents")
    mkpath(joinpath(contents, "MacOS"))
    mkpath(joinpath(contents, "Resources"))
    cp(icns, joinpath(contents, "Resources", "briareus.icns"); force = true)
    write(joinpath(contents, "Info.plist"), """
        <?xml version="1.0" encoding="UTF-8"?>
        <!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
        <plist version="1.0">
        <dict>
            <key>CFBundleName</key><string>Briareus</string>
            <key>CFBundleDisplayName</key><string>Briareus</string>
            <key>CFBundleIdentifier</key><string>io.github.stillyslalom.briareus</string>
            <key>CFBundleExecutable</key><string>briareus</string>
            <key>CFBundleIconFile</key><string>briareus</string>
            <key>CFBundlePackageType</key><string>APPL</string>
            <key>CFBundleShortVersionString</key><string>0.1</string>
            <key>NSHighResolutionCapable</key><true/>
        </dict>
        </plist>
        """)
    exe = joinpath(contents, "MacOS", "briareus")
    write(exe, "#!/bin/sh\nexec $(Base.shell_escape(shim)) --from-shortcut\n")
    chmod(exe, 0o755)
    return [app]
end

function desktop_entry(shim, png)
    dir = joinpath(get(ENV, "XDG_DATA_HOME", joinpath(homedir(), ".local", "share")), "applications")
    mkpath(dir)
    path = joinpath(dir, "briareus.desktop")
    write(path, """
        [Desktop Entry]
        Type=Application
        Name=Briareus
        Comment=Laser-target design
        Exec="$shim" --from-shortcut
        Icon=$png
        Terminal=false
        Categories=Science;Physics;
        """)
    return [path]
end

end # module BriareusLauncher
