# The `briareus` command: the `[apps.briareus]` entry of Project.toml runs `julia -m
# BriareusGUI`, whose `@main` is below. `julia -m` loads the whole package first (15-20 s,
# most of it GLMakie's and QML's package images), so a bare `briareus` in a terminal shows
# nothing until the window. A desktop launch goes through a script that `briareus
# --shortcut` writes (`install_shortcut`): it shows a splash with Qt's own `qml` tool
# within about a second, runs the app with `--from-shortcut` and the splash's status file
# in `BRIAREUS_SPLASH`, and sends the app's output to `log_path()`. The app writes its
# stage to that file, "quit" once the window shows, and "error: <message>" on a failure,
# which the splash then shows with the log's location.

const SPLASH_QML = joinpath(@__DIR__, "splash.qml")
const SPLASH_ENV = "BRIAREUS_SPLASH"

icons_dir() = joinpath(dirname(@__DIR__), "qml", "icons")

"""
    log_path() -> String

Where a launch from a shortcut (`briareus --shortcut`) writes its output:
`%LOCALAPPDATA%\\Briareus` on Windows, `~/Library/Logs/Briareus` on macOS,
`\$XDG_STATE_HOME/briareus` elsewhere.
"""
function log_path()
    dir = Sys.iswindows() ? joinpath(get(ENV, "LOCALAPPDATA", homedir()), "Briareus") :
          Sys.isapple() ? joinpath(homedir(), "Library", "Logs", "Briareus") :
          joinpath(get(ENV, "XDG_STATE_HOME", joinpath(homedir(), ".local", "state")), "briareus")
    return joinpath(dir, "briareus.log")
end

# ---------------------------------------------------------------- the command

function _main(args)
    "--shortcut" in args && (foreach(p -> println("created ", p), install_shortcut()); return 0)
    # the status file of the splash a shortcut's script shows
    splash = "--from-shortcut" in args ? get(ENV, SPLASH_ENV, "") : ""
    stage(text) = isempty(splash) || try write(splash, text) catch end
    files = filter(a -> !startswith(a, "--"), args)
    if !isempty(files) && !isfile(files[1])
        println(stderr, "briareus: no such file: ", files[1])
        return 1
    end
    try
        stage("Opening the window")
        sess = DesignSession()
        isempty(files) || load_session!(sess, files[1])
        briareus(sess; shown = () -> stage("quit"))
    catch err
        @error "Briareus could not start" exception = (err, catch_backtrace())
        # the message on one line (a TOML error's first line is only "TOML Parser error:")
        msg = replace(strip(sprint(showerror, err)), r"\s+" => " ")
        length(msg) > 160 && (msg = first(msg, 159) * "…")
        stage("error: Briareus could not start: $msg. The log is $(log_path()).")
        return 1
    end
    return 0
end

@static if isdefined(Base, Symbol("@main"))         # Julia 1.11; Pkg apps need 1.12
    function (@main)(args)
        return _main(args)
    end
end

# ---------------------------------------------------------------- the splash

function qml_tool()
    bin = joinpath(Qt6Declarative_jll.artifact_dir, "bin")
    for p in (joinpath(bin, "qml.exe"), joinpath(bin, "qml"),
              joinpath(bin, "qml.app", "Contents", "MacOS", "qml"))
        isfile(p) && return p
    end
    return nothing
end

file_url(path) = "file:///" * lstrip(replace(path, '\\' => '/'), '/')

# The variable that finds Qt's libraries for the `qml` tool, its folders (to go before
# the variable's value), and the tool's other variables.
function splash_env()
    var = Sys.iswindows() ? "PATH" : Sys.isapple() ? "DYLD_FALLBACK_LIBRARY_PATH" : "LD_LIBRARY_PATH"
    dirs = vcat(Qt6Declarative_jll.LIBPATH_list,
                joinpath(Qt6Declarative_jll.artifact_dir, Sys.iswindows() ? "bin" : "lib"))
    if Sys.isapple()
        # Qt's frameworks: the tool's rpaths name each framework folder of its own artifact
        # (`@loader_path/../lib/QtCore.framework`), so dyld's fallback needs those folders
        frameworks = [joinpath(d, f) for d in dirs if isdir(d) for f in readdir(d) if endswith(f, ".framework")]
        append!(dirs, frameworks)
    end
    others = ["QT_PLUGIN_PATH" => joinpath(Qt6Base_jll.artifact_dir, "plugins"),
              "QML_IMPORT_PATH" => joinpath(Qt6Declarative_jll.artifact_dir, "qml"),
              "QML_XHR_ALLOW_FILE_READ" => "1"]
    return var, unique(dirs), others
end

# The splash on the status file `status`, as the shortcut's script starts it (tests).
function splash_cmd(status; tool = qml_tool())
    var, dirs, others = splash_env()
    sep = Sys.iswindows() ? ';' : ':'
    haskey(ENV, var) && push!(dirs, ENV[var])
    icon = file_url(joinpath(icons_dir(), "briareus.png"))
    cmd = `$tool $SPLASH_QML -- $icon $(file_url(status))`
    return addenv(cmd, var => join(dirs, sep), others...)
end

# ---------------------------------------------------------------- the shortcut's script

const START_STAGE = "Loading the graphics libraries"

"""
    windows_script(shim, log; tool = qml_tool()) -> String

The batch script a Windows shortcut runs: the splash (when `tool`, Qt's `qml` tool, is
not `nothing`), then the app's `shim` with `--from-shortcut` and its output in `log`.
A failure the app did not report (it did not load) shows in the splash with the exit
code; without a splash, the log opens.
"""
function windows_script(shim, log; tool = qml_tool())
    b(s) = replace(s, "%" => "%%")                 # literal text in a batch script
    lines = ["@echo off",
             "rem Briareus from a shortcut, written by `briareus --shortcut`: a splash while",
             "rem the app loads, and the app's output in the log. Rerun `briareus --shortcut`",
             "rem after updating the app.",
             "chcp 65001 >nul",
             "setlocal",
             "set \"log=$(b(log))\"",
             "if not exist \"$(b(dirname(log)))\" mkdir \"$(b(dirname(log)))\"",
             "set \"BRIAREUS_SPLASH=\""]
    if tool !== nothing
        var, dirs, others = splash_env()
        icon = file_url(joinpath(icons_dir(), "briareus.png"))
        append!(lines, [
            "if not exist \"$(b(tool))\" goto run",
            "set \"BRIAREUS_SPLASH=%TEMP%\\briareus-splash-%RANDOM%%RANDOM%.txt\"",
            ">\"%BRIAREUS_SPLASH%\" echo $START_STAGE",
            "set \"status_url=file:///%BRIAREUS_SPLASH:\\=/%\"",
            "setlocal",
            "set \"$var=$(b(join(dirs, ';')));%$var%\"",
            ("set \"$k=$(b(v))\"" for (k, v) in others)...,
            "start \"\" /b \"$(b(tool))\" \"$(b(SPLASH_QML))\" -- \"$(b(icon))\" \"%status_url%\" >nul 2>&1",
            "endlocal"])
    end
    append!(lines, [
        ":run",
        "call \"$(b(shim))\" --from-shortcut >\"%log%\" 2>&1",
        "set \"code=%ERRORLEVEL%\"",
        "if \"%code%\"==\"0\" goto done",
        "if defined BRIAREUS_SPLASH goto failed",
        "start \"\" \"%log%\"",
        "goto done",
        ":failed",
        "findstr /b /c:\"error:\" \"%BRIAREUS_SPLASH%\" >nul || >\"%BRIAREUS_SPLASH%\" echo error: Briareus could not start (exit code %code%). The log is %log%.",
        "rem the splash reads the message before the file goes",
        "ping -n 3 127.0.0.1 >nul",
        ":done",
        "if defined BRIAREUS_SPLASH del \"%BRIAREUS_SPLASH%\" 2>nul",
        "exit /b %code%"])
    return join(lines, "\r\n") * "\r\n"
end

"""
    posix_script(shim, log; tool = qml_tool()) -> String

The shell script a macOS app bundle or a freedesktop menu entry runs, as
`windows_script` on Windows.
"""
function posix_script(shim, log; tool = qml_tool())
    q = Base.shell_escape
    lines = ["#!/bin/sh",
             "# Briareus from a shortcut, written by `briareus --shortcut`: a splash while the app",
             "# loads, and the app's output in the log. Rerun `briareus --shortcut` after updating",
             "# the app.",
             "log=$(q(log))",
             "mkdir -p \"\$(dirname \"\$log\")\"",
             "BRIAREUS_SPLASH="]
    if tool !== nothing
        var, dirs, others = splash_env()
        icon = file_url(joinpath(icons_dir(), "briareus.png"))
        append!(lines, [
            "if [ -x $(q(tool)) ] && BRIAREUS_SPLASH=\$(mktemp \"\${TMPDIR:-/tmp}/briareus-splash.XXXXXX\"); then",
            "    printf '%s' $(q(START_STAGE)) >\"\$BRIAREUS_SPLASH\"",
            "    (",
            "        export $var=$(q(join(dirs, ':')))\"\${$var:+:\$$var}\"",
            ("        export $k=$(q(v))" for (k, v) in others)...,
            "        exec $(q(tool)) $(q(SPLASH_QML)) -- $(q(icon)) \"file://\$BRIAREUS_SPLASH\"",
            "    ) >/dev/null 2>&1 &",
            "    export BRIAREUS_SPLASH",
            "fi"])
    end
    opener = Sys.isapple() ? "open" : "xdg-open"
    append!(lines, [
        "$(q(shim)) --from-shortcut >\"\$log\" 2>&1",
        "code=\$?",
        "if [ \"\$code\" -ne 0 ]; then",
        "    if [ -n \"\$BRIAREUS_SPLASH\" ]; then",
        "        grep -q '^error:' \"\$BRIAREUS_SPLASH\" ||",
        "            printf 'error: Briareus could not start (exit code %s). The log is %s.' \"\$code\" \"\$log\" >\"\$BRIAREUS_SPLASH\"",
        "        sleep 2    # the splash reads the message before the file goes",
        "    else",
        "        $opener \"\$log\" >/dev/null 2>&1 &",
        "    fi",
        "fi",
        "[ -z \"\$BRIAREUS_SPLASH\" ] || rm -f \"\$BRIAREUS_SPLASH\"",
        "exit \"\$code\""])
    return join(lines, '\n') * '\n'
end

# ---------------------------------------------------------------- shortcuts

"""
    install_shortcut() -> Vector{String}

Add Briareus to the desktop and the applications menu, as shortcuts that run a script
(the splash, then the installed `briareus` command with `--from-shortcut`), and return
what was written: on Windows, the script (beside `log_path()`) and shortcuts on the
desktop and in the Start menu (through `conhost --headless`, so no console window
opens); on macOS, `~/Applications/Briareus.app`; elsewhere, the script and a
freedesktop menu entry. The script names the Qt libraries of the installed app: rerun
after updating the app, or after moving the depot or the clone.
"""
function install_shortcut()
    shim = joinpath(first(DEPOT_PATH), "bin", Sys.iswindows() ? "briareus.bat" : "briareus")
    isfile(shim) || error("no briareus command in $(dirname(shim)): install it first (BriareusGUI README)")
    icons = icons_dir()
    if Sys.iswindows()
        script = joinpath(dirname(log_path()), "briareus-shortcut.cmd")
        mkpath(dirname(script))
        write(script, windows_script(shim, log_path()))
        return [script; windows_shortcuts(script, joinpath(icons, "briareus.ico"))]
    elseif Sys.isapple()
        return mac_app(posix_script(shim, log_path()), joinpath(icons, "briareus.icns"))
    else
        dir = joinpath(get(ENV, "XDG_DATA_HOME", joinpath(homedir(), ".local", "share")), "briareus")
        script = joinpath(dir, "briareus-shortcut.sh")
        mkpath(dir)
        write(script, posix_script(shim, log_path()))
        chmod(script, 0o755)
        return [script; desktop_entry(script, joinpath(icons, "briareus.png"))]
    end
end

# `folders`: a PowerShell expression for the folders that get a shortcut
function windows_shortcuts(script, icon; folders = "[Environment]::GetFolderPath('Desktop'), " *
                                                   "[Environment]::GetFolderPath('Programs')")
    q(s) = "'" * replace(s, "'" => "''") * "'"            # a PowerShell string literal
    conhost = joinpath(get(ENV, "SystemRoot", "C:\\Windows"), "System32", "conhost.exe")
    ps = """
        \$ProgressPreference = 'SilentlyContinue'
        \$shell = New-Object -ComObject WScript.Shell
        foreach (\$dir in $folders) {
            \$path = Join-Path \$dir 'Briareus.lnk'
            \$link = \$shell.CreateShortcut(\$path)
            \$link.TargetPath = $(q(conhost))
            \$link.Arguments = $(q("--headless \"$script\""))
            \$link.IconLocation = $(q(icon * ",0"))
            \$link.WorkingDirectory = [Environment]::GetFolderPath('MyDocuments')
            \$link.Description = 'Briareus laser-target design'
            \$link.Save()
            \$path
        }
        """
    encoded = base64encode(reinterpret(UInt8, transcode(UInt16, ps)))
    return readlines(`powershell -NoProfile -NonInteractive -EncodedCommand $encoded`)
end

# `script`: the text of the bundle's executable
function mac_app(script, icns; applications = joinpath(homedir(), "Applications"))
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
    write(exe, script)
    chmod(exe, 0o755)
    return [app]
end

function desktop_entry(script, png)
    dir = joinpath(get(ENV, "XDG_DATA_HOME", joinpath(homedir(), ".local", "share")), "applications")
    mkpath(dir)
    path = joinpath(dir, "briareus.desktop")
    write(path, """
        [Desktop Entry]
        Type=Application
        Name=Briareus
        Comment=Laser-target design
        Exec="$script"
        Icon=$png
        Terminal=false
        Categories=Science;Physics;
        """)
    return [path]
end
