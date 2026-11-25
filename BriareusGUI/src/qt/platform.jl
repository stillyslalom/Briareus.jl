# Qt Quick shell for the Briareus window (QML.jl + QMLMakie).
#
# Threading: Qt's event loop owns the main thread while the window is open. Calls from QML
# to Julia only change state and return. Background jobs hand their observable updates to
# `deliver`, which queues them; `br_tick` (a QML Timer, ~60 Hz) drains the queue on the
# GUI thread and refreshes what QML shows.
#
# Canvas: the QMLMakie canvas (`MakieArea`) is never destroyed while the window is
# open — destroying one leaves a dangling callback in jlqml that crashes at window
# teardown. The 3D viewer and the flux map are two figures shown in turn on it.
# The 3D viewer's mouse and keyboard gestures reach Makie's events through it
# (canvas/cad_camera.jl).

const QML_DIR = joinpath(@__DIR__, "..", "qml")

# ---------------------------------------------------------------- Qt setup

const _QT_READY = Ref(false)

# GLFW contexts (GLMakie windows and offscreen `colorbuffer` screens) and Qt's GL
# contexts must not share a process: once a GLFW context had existed, Qt's render
# thread crashed in the AMD driver (see HammerheadGUI). The window refuses to open next
# to a GLFW screen.
_glfw_screen_count() = count(s -> !(s isa GLMakie.Screen{QMLMakie.QMLWindow}), GLMakie.ALL_SCREENS)

function _check_no_glfw_screens()
    n = _glfw_screen_count()
    n == 0 && return
    throw(ArgumentError(
        "this Julia session has created $n GLMakie screen" * (n == 1 ? "" : "s") *
        " (a GLMakie window or an offscreen render), and GLMakie windows and the Qt " *
        "window cannot share a process (the graphics driver can crash). Restart Julia " *
        "to open the Briareus window."))
end

function _qt_init!()
    _QT_READY[] && return
    # GLMakie accepts a context only while the GLFW library is initialized, the QMLMakie
    # canvas's included. GLFW.jl 3.4 initialized it on load; 4.x does so on its first
    # call, which a QMLMakie canvas never makes. Initializing creates no GLFW context.
    GLMakie.GLFW.Init() ||
        error("GLFW failed to initialize; the Briareus window's canvases need it")
    if Sys.iswindows()
        # The FluentWinUI3 style plugin depends on a library outside the DLL search
        # path; load it first so the native Windows 11 style is used.
        impl = joinpath(Qt6Declarative_jll.artifact_dir, "bin",
                        "Qt6QuickControls2FluentWinUI3StyleImpl.dll")
        isfile(impl) && Libdl.dlopen(impl)
        # Qt reads msvcrt's copy of the environment, which Julia's ENV (Win32
        # SetEnvironmentVariable) does not update, so set both.
        style = get(ENV, "QT_QUICK_CONTROLS_STYLE", "FluentWinUI3")
        ENV["QT_QUICK_CONTROLS_STYLE"] = style
        ccall((:_putenv_s, "msvcrt"), Cint, (Cstring, Cstring), "QT_QUICK_CONTROLS_STYLE", style)
        ccall((:SetCurrentProcessExplicitAppUserModelID, "shell32"), Clong, (Cwstring,),
              APP_USER_MODEL_ID)
        _ENUM_QT_WINDOWS[] = _enum_qt_windows_ptr()
    end
    _QT_READY[] = true
    return
end

# Glyphs a canvas may show (ASCII, the symbols of the labels, Makie's tick minus). A
# glyph first laid out while the window is open reaches the atlas texture through a
# callback that runs outside Qt's render, where no GL context is current, and draws
# as garbage; glyphs in the atlas before the first render upload with it.
const _CANVAS_GLYPHS = vcat(Char(32):Char(126), collect("°±²³·×θφεΔσ–—−…→≤≥≈"))

# Rendering these glyphs takes seconds, so the warmed atlas is cached on disk next to
# the Makie atlas cache and loaded on the next start, provided no screen uses the
# session's atlas yet and the cache holds every glyph it has.
const _ATLAS_KEY = (2048, 64)                        # GLMakie's atlas
_atlas_cache_path() = Makie.get_cache_path(_ATLAS_KEY...) * ".briareusgui"

function _load_cached_atlas!(atlas)
    path = _atlas_cache_path()
    (isempty(atlas.font_render_callback) && isfile(path)) || return atlas
    cached = try
        Makie.load_texture_atlas(path)
    catch
        return atlas
    end
    all(k -> haskey(cached.mapping, k), keys(atlas.mapping)) || return atlas
    Makie.TEXTURE_ATLASES[_ATLAS_KEY] = cached
    return cached
end

function _warm_glyph_atlas!()
    atlas = _load_cached_atlas!(Makie.get_texture_atlas(_ATLAS_KEY...))
    n = length(atlas.mapping)
    fonts = Makie.theme(:fonts)
    for key in keys(fonts)
        font = Makie.to_font(fonts[key][])
        for c in _CANVAS_GLYPHS
            Makie.insert_glyph!(atlas, c, font)
        end
    end
    if length(atlas.mapping) > n
        path = _atlas_cache_path()
        try
            tmp = path * ".$(getpid()).tmp"
            Makie.store_texture_atlas(tmp, atlas)
            mv(tmp, path; force = true)
        catch
        end
    end
    return
end

# Compile what the canvases first do with a drive or a CAD import before the window opens,
# so that the first drive result, the first coloring by a field and the first parts left out
# of the calculations do not stall the window for most of a second each. The stall is JIT
# compilation, measured on a NIF CAD workspace, and the precompile workload of the package
# does not keep these methods. On the CPU side, a filled flux map and a field colorbar are
# drawn on throwaway canvases with a separate session. On the GL side, at the first frames
# of the window, the transparent ghost plot is drawn from an invisible degenerate triangle.
# The colorbar is shown as a drive field shows it, given a second scale as the next field
# coloring would change it (updates of the colormap and gradient textures), and hidden
# again.
function _warm_canvases!(sh)
    ws = DesignSession(resolution = :coarse)
    mc = map_canvas(ws)
    x, y = collect(1.0:4.0), collect(1.0:3.0)
    _update!(mc.heatmap, x, y, Float32[k + j for k in 1:4, j in 1:3]; colorrange = (1.0, 7.0))
    mc.ax.aspect[] = AxisAspect(2)
    reset_limits!(mc.ax)
    sc = scene_canvas(ws)
    _field_colorbar!(sc.colorbar, 1.3e14)
    _show_colorbar!(sc, true)
    # the canvas of the window
    c = sh.scene
    sh.sess.target.ghost[] === nothing && _update!(c.ghost, _NOTRI...; visible = true)
    if !c.colorbar.blockscene.visible[]
        _field_colorbar!(c.colorbar, 1.3e14)
        _show_colorbar!(c, true)
        ticks = Ref(0)
        once = Ref{Any}(nothing)
        once[] = on(events(c.fig.scene).tick) do _
            # a tick precedes each frame's drawing: a second scale at the second, hidden
            # at the third, unless a drive field colors the body by then
            k = (ticks[] += 1)
            k < 2 && return
            field = _body_coloring(sh.sess) !== nothing
            if k == 2
                field || (_field_colorbar!(c.colorbar, 2.6e14); c.dirty[] = true)
                return
            end
            off(once[])
            field || (_show_colorbar!(c, false); c.dirty[] = true)
        end
    end
    return nothing
end

# A colorbar set as a drive field sets it (`_draw_body!`).
function _field_colorbar!(cb, hi)
    cb.colormap[] = FIELD_COLORMAP
    cb.limits[] = (0.0, hi)
    cb.ticks[] = Makie.automatic
    cb.label[] = Controllers.DRIVE_FIELD_LABELS[:incident]
    return cb
end

# After the window closes Qt has destroyed its GL contexts: forget the canvas's
# screens so GLMakie neither renders to them nor frees them at exit.
function _release_qml_screens!()
    for s in GLMakie.ALL_SCREENS
        s isa GLMakie.Screen{QMLMakie.QMLWindow} && (s.glscreen.context.valid = false)
    end
    filter!(s -> !(s isa GLMakie.Screen{QMLMakie.QMLWindow}), GLMakie.ALL_SCREENS)
    for key in collect(keys(GLMakie.atlas_texture_cache))
        key[2] isa QMLMakie.QMLWindow || continue
        _, callback = pop!(GLMakie.atlas_texture_cache, key)
        Makie.remove_font_render_callback!(key[1], callback)
    end
    return
end

# Convert QML FileDialog URLs to local paths.
function _url_to_path(url::AbstractString)
    s = String(url)
    startswith(s, "file://") || return s
    s = s[8:end]                                   # drop "file://"
    Sys.iswindows() && startswith(s, "/") && occursin(r"^/[A-Za-z]:", s) && (s = s[2:end])
    return _percent_decode(s)
end

function _percent_decode(s::String)
    io = IOBuffer()
    i = firstindex(s)
    while i <= lastindex(s)
        if s[i] == '%' && i + 2 <= lastindex(s)
            write(io, parse(UInt8, s[i+1:i+2]; base = 16))
            i += 3
        else
            write(io, s[i])
            i = nextind(s, i)
        end
    end
    return String(take!(io))
end

# The Briareus icon (docs/make_icons.jl writes it).
const ICON_ICO = joinpath(QML_DIR, "icons", "briareus.ico")
const _ICON_TRIES = Ref(0)

# The taskbar identity of the process (Windows): set before the first window, it keeps
# the window out of julia.exe's taskbar group, whose button shows julia.exe's icon.
const APP_USER_MODEL_ID = "Briareus.DesignWindow"

# Windows: show the icon in the title bar and taskbar. Qt Quick windows have no icon
# property and QML.jl does not expose QGuiApplication's, so the icon is sent to this
# process's top-level Qt windows on the first ticks after the window opens, and set on their
# window class, so that a native window Qt creates later carries it too. Windows are found
# by process instead of by title, since the title changes while a design loads. Returns true
# once a visible window has the icon. Elsewhere the window keeps the platform default.
function _set_window_icon()
    Sys.iswindows() || return true
    empty!(_QT_HWNDS)
    ccall((:EnumWindows, "user32"), Cint, (Ptr{Cvoid}, Int), _ENUM_QT_WINDOWS[], 0)
    isempty(_QT_HWNDS) && return false
    shown = false
    for (which, index, metric) in ((0, -34, 49), (1, -14, 11))  # ICON_SMALL/GCLP_HICONSM/SM_CXSMICON, ICON_BIG/GCLP_HICON/SM_CXICON
        px = ccall((:GetSystemMetrics, "user32"), Cint, (Cint,), metric)
        icon = ccall((:LoadImageW, "user32"), Ptr{Cvoid},
                     (Ptr{Cvoid}, Cwstring, Cuint, Cint, Cint, Cuint),
                     C_NULL, ICON_ICO, 1, px, px, 0x10)      # IMAGE_ICON, LR_LOADFROMFILE
        icon == C_NULL && continue
        for hwnd in _QT_HWNDS
            ccall((:SetClassLongPtrW, "user32"), Csize_t, (Ptr{Cvoid}, Cint, Ptr{Cvoid}),
                  hwnd, index, icon)
            ccall((:SendMessageW, "user32"), Ptr{Cvoid}, (Ptr{Cvoid}, Cuint, Csize_t, Ptr{Cvoid}),
                  hwnd, 0x0080, which, icon)                 # WM_SETICON
            shown |= ccall((:IsWindowVisible, "user32"), Cint, (Ptr{Cvoid},), hwnd) != 0
        end
    end
    return shown
end

const _QT_HWNDS = Ptr{Cvoid}[]

# The EnumWindows callback collects the top-level Qt windows of this process (class
# "Qt…QWindow…").
function _collect_qt_window(hwnd::Ptr{Cvoid}, ::Int)::Cint
    pid = Ref{UInt32}(0)
    ccall((:GetWindowThreadProcessId, "user32"), UInt32, (Ptr{Cvoid}, Ref{UInt32}), hwnd, pid)
    pid[] == getpid() || return 1
    buf = Vector{UInt16}(undef, 256)
    n = ccall((:GetClassNameW, "user32"), Cint, (Ptr{Cvoid}, Ptr{UInt16}, Cint), hwnd, buf, 256)
    cls = transcode(String, buf[1:max(n, 0)])
    startswith(cls, "Qt") && occursin("QWindow", cls) && push!(_QT_HWNDS, hwnd)
    return 1
end

const _ENUM_QT_WINDOWS = Ref(C_NULL)
_enum_qt_windows_ptr() = @cfunction(_collect_qt_window, Cint, (Ptr{Cvoid}, Int))
