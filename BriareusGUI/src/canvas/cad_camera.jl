# CAD-style 3D navigation for the viewer: the camera is a pose (eye, look-at point, up)
# that the canvas sets as the scene's view and projection matrices itself, so nothing
# clips the model: near and far follow the eye's distance to the scene's bounding
# sphere. Mouse gestures (CadCamera docstring) act about points picked by ray casting
# into the target (the canvas supplies `pick`), the way solid modelers do: the wheel
# zooms toward the point under the pointer, drags orbit about the selection (or the
# point under the pointer) and pan the grabbed point with the pointer.

const Vec3d_ = Makie.Vec3d
const _ORBIT_RATE = π              # radians per viewport height of drag
const _ZOOM_STEP = 0.85            # distance factor per wheel notch
const _CLICK_SLOP = 4.0            # pixels a press may move and still be a click
const _DOUBLE_CLICK = 0.4          # seconds between the clicks of a double click

"""
    CadCamera

Camera pose and mouse state of the 3D viewer. Pose: `eye`, `lookat`, `up`, vertical
field of view `fov` (degrees), `ortho` (orthographic projection, at the scale the
perspective view has at the look-at distance). `bounds` is the bounding sphere
`(center, radius)` of everything drawn, which sets the near and far planes; `home` is
the sphere `fit_camera!` shows by default (the target). `turntable` keeps world
`z` up while orbiting (otherwise the rotation is free about screen axes).

Gestures (wired by [`attach_cad_camera!`](@ref)):
- wheel: zoom toward the point under the pointer (it stays under the pointer);
- left or middle drag: orbit about the selection's center, else the point under the
  pointer at the press, else the look-at point;
- right drag, Shift + left/middle drag or Ctrl + middle drag: pan (the grabbed point
  follows the pointer);
- click: `on_click(origin, direction, pixel)` (selection; `click_ctrl` says whether Ctrl
  was held, to add to it); double click: `on_double`;
- `f`: fit the scene; arrow keys orbit 15°.

`pick(origin, direction)` returns the distance to the first surface along the ray, or
`Inf`. `pivot()` returns the rotation center to use (the selection's center) or
`nothing`. `on_change()` runs after every pose change.
"""
mutable struct CadCamera
    eye::Vec3d_
    lookat::Vec3d_
    up::Vec3d_
    fov::Float64
    ortho::Bool
    turntable::Bool
    bounds::Tuple{Vec3d_,Float64}
    home::Tuple{Vec3d_,Float64}
    viewport::NTuple{4,Float64}           # x, y, width, height (pixels)
    pick::Any
    pivot::Any
    on_click::Any
    on_double::Any
    on_change::Any
    # drag state
    mode::Symbol                          # :none, :press, :orbit, :pan
    button::Any
    press_pos::NTuple{2,Float64}
    last_pos::NTuple{2,Float64}
    center::Vec3d_                        # orbit pivot of the drag
    grab::Vec3d_                          # pan point of the drag
    last_click::Float64
    click_ctrl::Bool                      # Ctrl held at the press
end

function CadCamera(; bounds = (Vec3d_(0, 0, 0), 1.0), fov = 30.0)
    c = CadCamera(Vec3d_(1, 1, 1), Vec3d_(0, 0, 0), Vec3d_(0, 0, 1), Float64(fov), false, true,
                  (Vec3d_(bounds[1]...), Float64(bounds[2])),
                  (Vec3d_(bounds[1]...), Float64(bounds[2])), (0.0, 0.0, 800.0, 600.0),
                  (o, d) -> Inf, () -> nothing, (o, d, px) -> nothing, (o, d, px) -> nothing,
                  () -> nothing, :none, nothing, (0.0, 0.0), (0.0, 0.0), Vec3d_(0, 0, 0),
                  Vec3d_(0, 0, 0), -Inf, false)
    orient_camera!(c, :iso)
    return c
end

_normalized(v) = v / norm(v)

# Camera frame: forward (toward the look-at point), right, true up.
function camera_frame(c::CadCamera)
    f = _normalized(c.lookat - c.eye)
    r = cross(f, c.up)
    r = norm(r) < 1e-9 ? _normalized(cross(f, abs(f[3]) < 0.9 ? Vec3d_(0, 0, 1) : Vec3d_(1, 0, 0))) : _normalized(r)
    return f, r, cross(r, f)
end

_aspect(c::CadCamera) = c.viewport[3] / max(c.viewport[4], 1.0)
_tanhalf(c::CadCamera) = tand(c.fov / 2)

"""
    clip_planes(c::CadCamera) -> (near, far)

Near and far distances that hold the whole bounding sphere: in front of it, from its
near side to its far side; inside it, from a small fraction of the look-at distance
(so a zoomed-in close-up is never clipped) to its far side.
"""
function clip_planes(c::CadCamera)
    C, R = c.bounds
    f, _, _ = camera_frame(c)
    zc = dot(C - c.eye, f)                         # depth of the sphere center
    far = zc + 1.02R
    if c.ortho
        return zc - 1.02R, far
    end
    far = max(far, 1e-6)
    d = norm(c.lookat - c.eye)
    near = max(zc - 1.02R, 1e-3 * d, 1e-5 * far)
    return min(near, 0.5far), far
end

"""
    view_matrices(c::CadCamera) -> (view, projection)
"""
function view_matrices(c::CadCamera)
    V = Makie.lookat(c.eye, c.lookat, camera_frame(c)[3])
    near, far = clip_planes(c)
    if c.ortho
        hh = norm(c.lookat - c.eye) * _tanhalf(c)
        hw = hh * _aspect(c)
        P = Makie.orthographicprojection(-hw, hw, -hh, hh, near, far)
    else
        P = Makie.perspectiveprojection(c.fov, _aspect(c), near, far)
    end
    return Makie.Mat4d(V), Makie.Mat4d(P)
end

"""
    pixel_ray(c::CadCamera, px) -> (origin, direction)

The ray through pixel `px` (scene-local pixels, origin at the bottom left).
"""
function pixel_ray(c::CadCamera, px)
    f, r, u = camera_frame(c)
    x = 2px[1] / c.viewport[3] - 1
    y = 2px[2] / c.viewport[4] - 1
    t = _tanhalf(c)
    if c.ortho
        hh = norm(c.lookat - c.eye) * t
        return c.eye + x * hh * _aspect(c) * r + y * hh * u - 10c.bounds[2] * f, f
    end
    return c.eye, _normalized(f + x * t * _aspect(c) * r + y * t * u)
end

"""
    project_point(c::CadCamera, p) -> (x, y)

Scene-local pixel position of the point `p`.
"""
function project_point(c::CadCamera, p)
    f, r, u = camera_frame(c)
    d = p - c.eye
    z = dot(d, f)
    t = _tanhalf(c)
    s = c.ortho ? norm(c.lookat - c.eye) * t : max(z, 1e-12) * t
    x = dot(d, r) / (s * _aspect(c))
    y = dot(d, u) / s
    return ((x + 1) * c.viewport[3] / 2, (y + 1) * c.viewport[4] / 2)
end

# The surface point under a pixel, or the point on the ray at the look-at depth.
function _point_under(c::CadCamera, px)
    o, d = pixel_ray(c, px)
    τ = c.pick(o, d)
    isfinite(τ) && return o + τ * d, true
    f = camera_frame(c)[1]
    s = dot(c.lookat - o, f) / dot(d, f)
    return o + s * d, false
end

"""
    zoom_camera!(c::CadCamera, px, factor)

Scale the eye's distance by `factor` about the point under pixel `px`, which stays
under the pointer (a homothety of the eye and look-at point about it).
"""
function zoom_camera!(c::CadCamera, px, factor::Real)
    p, _ = _point_under(c, px)
    c.eye = p + factor * (c.eye - p)
    c.lookat = p + factor * (c.lookat - p)
    c.on_change()
    return c
end

_rotate(v, axis, θ) = v * cos(θ) + cross(axis, v) * sin(θ) + axis * dot(axis, v) * (1 - cos(θ))

"""
    orbit_camera!(c::CadCamera, dx, dy; center = c.lookat)

Rotate the view about `center` by `dx` radians about world `z` (turntable) or the
screen's up axis, and `dy` radians about the screen's right axis.
"""
function orbit_camera!(c::CadCamera, dx::Real, dy::Real; center = c.lookat)
    _, r, u = camera_frame(c)
    ax1 = c.turntable ? Vec3d_(0, 0, 1) : u
    rot(v) = _rotate(_rotate(v, ax1, -dx), _normalized(_rotate(r, ax1, -dx)), -dy)
    c.eye = center + rot(c.eye - center)
    c.lookat = center + rot(c.lookat - center)
    c.up = rot(u)
    c.on_change()
    return c
end

"""
    pan_camera!(c::CadCamera, dpx, dpy; depth = distance to the look-at point)

Translate the view so a point at `depth` moves by `(dpx, dpy)` pixels with the
pointer.
"""
function pan_camera!(c::CadCamera, dpx::Real, dpy::Real; depth = norm(c.lookat - c.eye))
    _, r, u = camera_frame(c)
    s = 2 * _tanhalf(c) * (c.ortho ? norm(c.lookat - c.eye) : depth) / c.viewport[4]
    δ = -s * (dpx * r + dpy * u)
    c.eye += δ
    c.lookat += δ
    c.on_change()
    return c
end

"""
    fit_camera!(c::CadCamera, center = home center, radius = home radius)

Look at `center` from the current direction, at the distance that shows a sphere of
`radius` whole.
"""
function fit_camera!(c::CadCamera, center = c.home[1], radius = c.home[2])
    f, _, _ = camera_frame(c)
    half = min(c.fov, c.fov * _aspect(c)) / 2
    d = 1.08 * radius / sind(half)
    c.lookat = Vec3d_(center...)
    c.eye = c.lookat - d * f
    c.on_change()
    return c
end

"Standard view directions: the direction the camera looks along, and up."
const STANDARD_VIEWS = Dict(
    :iso => (Vec3d_(1, 1, -0.75), Vec3d_(0, 0, 1)),      # from −x −y, the cutaway side
    :front => (Vec3d_(0, 1, 0), Vec3d_(0, 0, 1)),        # looking along +y
    :back => (Vec3d_(0, -1, 0), Vec3d_(0, 0, 1)),
    :right => (Vec3d_(-1, 0, 0), Vec3d_(0, 0, 1)),       # looking along -x
    :left => (Vec3d_(1, 0, 0), Vec3d_(0, 0, 1)),
    :top => (Vec3d_(0, 0, -1), Vec3d_(0, 1, 0)),
    :bottom => (Vec3d_(0, 0, 1), Vec3d_(0, 1, 0)))

"""
    orient_camera!(c::CadCamera, name::Symbol) / orient_camera!(c, direction, up)

Look along a standard direction (`STANDARD_VIEWS`) or a given one, fitted to the scene.
"""
orient_camera!(c::CadCamera, name::Symbol) = orient_camera!(c, STANDARD_VIEWS[name]...)

function orient_camera!(c::CadCamera, dir, up)
    f = _normalized(Vec3d_(dir...))
    u = Vec3d_(up...)
    abs(dot(f, _normalized(u))) > 0.999 && (u = abs(f[3]) < 0.9 ? Vec3d_(0, 0, 1) : Vec3d_(0, 1, 0))
    c.up = u
    d = max(norm(c.lookat - c.eye), 1e-6)
    c.eye = c.lookat - d * f
    return fit_camera!(c)
end

# ---------------------------------------------------------------- mouse

_mods(ev) = (shift = Makie.Keyboard.left_shift in ev.keyboardstate || Makie.Keyboard.right_shift in ev.keyboardstate,
             ctrl = Makie.Keyboard.left_control in ev.keyboardstate || Makie.Keyboard.right_control in ev.keyboardstate)

function _local(c::CadCamera, mp)
    return (Float64(mp[1]) - c.viewport[1], Float64(mp[2]) - c.viewport[2])
end

_inside(c::CadCamera, p) = 0 <= p[1] <= c.viewport[3] && 0 <= p[2] <= c.viewport[4]

"""
    camera_press!(c, px, button; shift = false, ctrl = false)
    camera_move!(c, px)
    camera_release!(c, px)
    camera_scroll!(c, px, notches)

Gesture handlers (scene-local pixels); [`attach_cad_camera!`](@ref) feeds them from
the scene's events, and tests call them directly.
"""
function camera_press!(c::CadCamera, px, button::Symbol; shift = false, ctrl = false)
    c.press_pos = c.last_pos = (Float64(px[1]), Float64(px[2]))
    c.button = button
    c.click_ctrl = ctrl
    pan = button === :right || (shift && button in (:left, :middle)) || (ctrl && button === :middle)
    c.mode = :press
    p, hit = _point_under(c, px)
    if pan
        c.grab = p
        c.mode = :pan_pending
    else
        q = c.pivot()
        c.center = q !== nothing ? Vec3d_(q...) : hit ? p : c.lookat
        c.mode = :orbit_pending
    end
    return c
end

function camera_move!(c::CadCamera, px)
    c.mode === :none && return false
    x, y = Float64(px[1]), Float64(px[2])
    if c.mode in (:orbit_pending, :pan_pending)
        hypot(x - c.press_pos[1], y - c.press_pos[2]) < _CLICK_SLOP && return true
        c.mode = c.mode === :orbit_pending ? :orbit : :pan
    end
    dx, dy = x - c.last_pos[1], y - c.last_pos[2]
    c.last_pos = (x, y)
    if c.mode === :orbit
        s = _ORBIT_RATE / max(c.viewport[4], 1.0)
        orbit_camera!(c, s * dx, -s * dy; center = c.center)
    elseif c.mode === :pan
        f = camera_frame(c)[1]
        pan_camera!(c, dx, dy; depth = max(dot(c.grab - c.eye, f), 1e-6))
    end
    return true
end

function camera_release!(c::CadCamera, px)
    mode, c.mode = c.mode, :none
    mode in (:orbit_pending, :pan_pending) || return c
    c.button === :left || return c
    o, d = pixel_ray(c, px)
    now = time()
    if now - c.last_click < _DOUBLE_CLICK
        c.last_click = -Inf
        c.on_double(o, d, px)
    else
        c.last_click = now
        c.on_click(o, d, px)
    end
    return c
end

camera_scroll!(c::CadCamera, px, notches::Real) = zoom_camera!(c, px, _ZOOM_STEP^notches)

const _BUTTONS = Dict(Makie.Mouse.left => :left, Makie.Mouse.middle => :middle,
                      Makie.Mouse.right => :right)

"""
    attach_cad_camera!(scene, c::CadCamera)

Drive `scene`'s camera from `c`: set its matrices on every pose change and route the
scene's mouse and keyboard events to the gestures.
"""
function attach_cad_camera!(scene, c::CadCamera)
    ev = Makie.events(scene)
    apply = function ()
        V, P = view_matrices(c)
        cam = Makie.camera(scene)
        f, _, u = camera_frame(c)
        cam.eyeposition[] = Makie.Vec3f(c.eye)
        cam.upvector[] = Makie.Vec3f(u)
        cam.view_direction[] = Makie.Vec3f(f)
        Makie.set_proj_view!(cam, P, V)
    end
    user_change = c.on_change
    c.on_change = () -> (apply(); user_change())
    on(scene, Makie.viewport(scene)) do vp
        o, w = Makie.origin(vp), Makie.widths(vp)
        c.viewport = (Float64(o[1]), Float64(o[2]), Float64(w[1]), Float64(w[2]))
        apply()
    end
    on(scene, ev.mousebutton; priority = 100) do e
        b = get(_BUTTONS, e.button, nothing)
        b === nothing && return Consume(false)
        px = _local(c, ev.mouseposition[])
        if e.action == Makie.Mouse.press
            _inside(c, px) || return Consume(false)
            m = _mods(ev)
            camera_press!(c, px, b; m.shift, m.ctrl)
            return Consume(true)
        elseif e.action == Makie.Mouse.release && c.mode !== :none
            camera_release!(c, px)
            return Consume(true)
        end
        return Consume(false)
    end
    on(scene, ev.mouseposition; priority = 100) do mp
        c.mode === :none && return Consume(false)
        return Consume(camera_move!(c, _local(c, mp)))
    end
    on(scene, ev.scroll; priority = 100) do s
        px = _local(c, ev.mouseposition[])
        _inside(c, px) || return Consume(false)
        s[2] == 0 && return Consume(false)
        camera_scroll!(c, px, sign(s[2]) * min(abs(s[2]), 3))
        return Consume(true)
    end
    on(scene, ev.keyboardbutton; priority = 100) do e
        e.action in (Makie.Keyboard.press, Makie.Keyboard.repeat) || return Consume(false)
        step = deg2rad(15)
        if e.key == Makie.Keyboard.f
            fit_camera!(c)
        elseif e.key == Makie.Keyboard.left
            orbit_camera!(c, -step, 0.0)
        elseif e.key == Makie.Keyboard.right
            orbit_camera!(c, step, 0.0)
        elseif e.key == Makie.Keyboard.up
            orbit_camera!(c, 0.0, -step)
        elseif e.key == Makie.Keyboard.down
            orbit_camera!(c, 0.0, step)
        else
            return Consume(false)
        end
        return Consume(true)
    end
    apply()
    return c
end
