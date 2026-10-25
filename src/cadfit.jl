# Analytic surfaces in imported (CAD) meshes: patches that lie on one plane, cylinder,
# cone, sphere or torus are found by region growing with least-squares fits, and the curved ones
# are attached to their components as `CurvedSurface`s, so the radiation model
# integrates on the true surface as it does for primitive targets.

"""
    SurfaceRegion

Patches of a mesh lying on one analytic surface, as found by [`fit_surfaces`](@ref):

- `kind`: `:plane`, `:cylinder`, `:cone`, `:sphere` or `:torus`
- `surface`: the fitted `CylinderSurface`, `ConeSurface`, `SphereSurface` or
  `TorusSurface` (`nothing` for a plane)
- `centroid`: area-weighted centroid of the patches (mm)
- `normal`: a plane's unit normal, on the side most of its area faces (zero for
  curved surfaces)
- `patches`: indices of the patches in the mesh
- `area`: total area of the patches (mm², flat facets)
- `rms`, `maxdev`: rms and largest distance of the patches' vertices from the surface
  (mm)
- `facing`: for curved surfaces, `-1` when every patch faces the axis, center or spine
  (the inside of a hohlraum wall or of a cone, a fillet in a cavity's corner), `+1` when every patch faces away from it, `0`
  when mixed;
  for planes `+1`, or `0` when some patches face the other way
"""
struct SurfaceRegion
    kind::Symbol
    surface::Union{Nothing,CylinderSurface,ConeSurface,SphereSurface,TorusSurface}
    centroid::Vec3
    normal::Vec3
    patches::Vector{Int}
    area::Float64
    rms::Float64
    maxdev::Float64
    facing::Int
end

function Base.show(io::IO, r::SurfaceRegion)
    p3(x) = @sprintf("(%.4g, %.4g, %.4g)", x[1], x[2], x[3])
    s = r.surface
    geo = r.kind === :plane ? "normal $(p3(r.normal))" :
          r.kind === :cylinder ? @sprintf("radius %.6g mm, axis %s through %s", s.radius, p3(s.axis), p3(s.origin)) :
          r.kind === :cone ? @sprintf("half-angle %.6g°, axis %s from apex %s", s.half_angle, p3(s.axis), p3(s.apex)) :
          r.kind === :torus ? @sprintf("radii %.6g and %.6g mm, axis %s through %s", s.major_radius, s.minor_radius,
                                       p3(s.axis), p3(s.center)) :
          @sprintf("radius %.6g mm, center %s", s.radius, p3(s.center))
    side = r.facing == 0 ? ", mixed facing" :
           r.kind === :plane ? "" : r.facing < 0 ? ", facing in" : ", facing out"
    @printf(io, "SurfaceRegion(:%s, %s, %d patches, %.4g mm², rms %.2g mm, max %.2g mm%s)",
            r.kind, geo, length(r.patches), r.area, r.rms, r.maxdev, side)
end

# A fitted surface: plane (point c, unit normal u), cylinder (axis point c, unit axis
# u, radius R), cone (axis point c, unit axis u, radius R at c and slope k = dρ/dz along
# u: ρ = R + k z, either sign, so a cone tends to a cylinder as k → 0), sphere
# (center c, radius R) or torus (center c, unit axis u, major radius R, minor radius k).
struct _SurfaceModel
    kind::Symbol
    c::Vec3
    u::Vec3
    R::Float64
    k::Float64
end
_SurfaceModel(kind, c, u, R) = _SurfaceModel(kind, c, u, R, 0.0)

function _model_distance(m::_SurfaceModel, p)
    d = p - m.c
    m.kind === :plane && return dot(d, m.u)
    m.kind === :cylinder && return norm(d - dot(d, m.u) * m.u) - m.R
    if m.kind === :cone
        z = dot(d, m.u)
        return (norm(d - z * m.u) - m.R - m.k * z) / sqrt(1 + m.k^2)
    end
    if m.kind === :torus
        z = dot(d, m.u)
        return hypot(norm(d - z * m.u) - m.R, z) - m.k
    end
    return norm(d) - m.R
end

function _model_normal(m::_SurfaceModel, p)
    m.kind === :plane && return m.u
    d = p - m.c
    m.kind === :sphere && return normalize(d)
    z = dot(d, m.u)
    q = d - z * m.u
    ρ̂ = normalize(q)
    m.kind === :torus && return normalize((norm(q) - m.R) * ρ̂ + z * m.u)
    return m.kind === :cone ? (ρ̂ - m.k * m.u) / sqrt(1 + m.k^2) : ρ̂
end

_model_surface(m::_SurfaceModel) =
    m.kind === :cylinder ? CylinderSurface(m.c, m.u, m.R) :
    m.kind === :cone ? ConeSurface(m.c - (m.R / m.k) * m.u, sign(m.k) * m.u, atand(abs(m.k))) :
    m.kind === :sphere ? SphereSurface(m.c, m.R) :
    m.kind === :torus ? TorusSurface(m.c, m.u, m.R, m.k) : nothing

# Radius of curvature across the surface (a cone's at its reference circle, R/cos θ; a
# torus's across its tube, the minor radius).
_model_radius(m::_SurfaceModel) = m.kind === :cone ? m.R * sqrt(1 + m.k^2) : m.kind === :torus ? m.k : m.R

# Least-squares fits to points X (with the patch normals N and areas W for the
# cylinder's and the cone's axis). `nothing` when the points do not determine the
# surface, or give a radius beyond 1000 times their extent (a plane in disguise; for a
# cone, also an apex that far: a cylinder in disguise).
function _fit_model(kind::Symbol, X::Vector{Vec3}, N::Vector{Vec3}, W::Vector{Float64}, Cn = nothing)
    m = sum(X) / length(X)
    ext = maximum(x -> norm(x - m), X)
    ext > 0 || return nothing
    if kind === :torus
        length(X) >= 8 && length(N) >= 3 || return nothing
        return _fit_torus(X, N, W, Cn, m, ext)
    elseif kind === :plane
        length(X) >= 3 || return nothing
        C = sum((x - m) * (x - m)' for x in X)
        return _SurfaceModel(:plane, m, Vec3(eigen(Symmetric(Matrix(C))).vectors[:, 1]), Inf)
    elseif kind === :sphere
        length(X) >= 5 || return nothing
        return _fit_sphere(X, m, ext)
    elseif kind === :cone
        length(X) >= 7 || return nothing
        return _fit_cone(X, N, W, m, ext)
    else
        length(X) >= 6 || return nothing
        return _fit_cylinder(X, N, W, m, ext)
    end
end

# Algebraic fit (|x|² = 2c·x + k, linear in c and k), refined by Gauss–Newton on the
# geometric distances |x − c| − R.
function _fit_sphere(X, m, ext)
    M = zeros(SMatrix{4,4,Float64})
    g = zeros(SVector{4,Float64})
    for x in X
        y = x - m
        a = SVector(2y[1], 2y[2], 2y[3], 1.0)
        M += a * a'
        g += a * dot(y, y)
    end
    s = _solve_small(M, g)
    s === nothing && return nothing
    c = Vec3(s[1], s[2], s[3])
    R = sqrt(max(s[4] + dot(c, c), 0.0))
    for _ in 1:20
        (isfinite(R) && 0 < R < 1e3ext) || return nothing
        M = zeros(SMatrix{4,4,Float64})
        g = zeros(SVector{4,Float64})
        for x in X
            d = x - m - c
            l = norm(d)
            l > 0 || continue
            J = SVector(-d[1] / l, -d[2] / l, -d[3] / l, -1.0)
            M += J * J'
            g += J * (l - R)
        end
        δ = _solve_small(M, -g)
        δ === nothing && return nothing
        c += Vec3(δ[1], δ[2], δ[3])
        R += δ[4]
        norm(δ) <= 1e-13 * R && break
    end
    (isfinite(R) && 0 < R < 1e3ext) || return nothing
    return _SurfaceModel(:sphere, m + c, zero(Vec3), R)
end

# Axis from the patch normals (perpendicular to it on a cylinder: the eigenvector of
# Σ w n nᵀ with the least eigenvalue), a circle fit across the axis, then Gauss–Newton
# on the geometric distances over the axis direction, its position and the radius.
function _fit_cylinder(X, N, W, m, ext)
    S = zeros(SMatrix{3,3,Float64})
    for (n, w) in zip(N, W)
        S += w * (n * n')
    end
    a = Vec3(eigen(Symmetric(Matrix(S))).vectors[:, 1])
    e1, e2 = _perp_basis(a)
    M = zeros(SMatrix{3,3,Float64})
    g = zeros(SVector{3,Float64})
    for x in X
        y = x - m
        u, v = dot(y, e1), dot(y, e2)
        q = SVector(2u, 2v, 1.0)
        M += q * q'
        g += q * (u^2 + v^2)
    end
    s = _solve_small(M, g)
    s === nothing && return nothing
    o = s[1] * e1 + s[2] * e2                       # relative to m
    R = sqrt(max(s[3] + s[1]^2 + s[2]^2, 0.0))
    for _ in 1:20
        (isfinite(R) && 0 < R < 1e3ext) || return nothing
        e1, e2 = _perp_basis(a)
        M5 = zeros(SMatrix{5,5,Float64})
        g5 = zeros(SVector{5,Float64})
        for x in X
            d = x - m - o
            h = dot(d, a)
            q = d - h * a
            ρ = norm(q)
            ρ > 0 || continue
            û = q / ρ
            ue1, ue2 = dot(û, e1), dot(û, e2)
            J = SVector(-h * ue1, -h * ue2, -ue1, -ue2, -1.0)
            M5 += J * J'
            g5 += J * (ρ - R)
        end
        δ = _solve_small(M5, -g5)
        δ === nothing && return nothing
        a = normalize(a + δ[1] * e1 + δ[2] * e2)
        o += δ[3] * e1 + δ[4] * e2
        o -= dot(o, a) * a                          # the axis point nearest m
        R += δ[5]
        norm(δ) <= 1e-13 * max(R, 1.0) && break
    end
    (isfinite(R) && 0 < R < 1e3ext) || return nothing
    return _SurfaceModel(:cylinder, m + o, a, R)
end

# Cone: from two axis directions, the axis of the circle the normals trace (on a cone
# they make a fixed angle with the axis: the least principal direction of the centered
# normals, right for a band around the axis) and the cylinder's (the direction the
# normals are perpendicular to, right for a small patch of a gentle cone), an algebraic
# fit of ρ = R + k z across the axis (u² + v² linear in the axis point, R², 2Rk and k²;
# or without the k² z² term, z being centered, when that fits better: on a band of two
# rims z² is a combination of 1 and z), then Gauss–Newton (step halving when the
# residual grows) on the geometric distances (ρ − R − kz)/√(1+k²) over the axis
# direction, its position, R and k. The better of the two is kept.
function _fit_cone(X, N, W, m, ext)
    wt = sum(W)
    wt > 0 || return nothing
    nbar = sum(w * n for (n, w) in zip(N, W)) / wt
    Sc = zeros(SMatrix{3,3,Float64})
    Su = zeros(SMatrix{3,3,Float64})
    for (n, w) in zip(N, W)
        Sc += w * ((n - nbar) * (n - nbar)')
        Su += w * (n * n')
    end
    best, bestc = nothing, Inf
    for S in (Sc, Su)
        a = Vec3(eigen(Symmetric(Matrix(S))).vectors[:, 1])
        g = _cone_guess(X, m, a, true)
        g1 = _cone_guess(X, m, a, false)
        if g === nothing || (g1 !== nothing && _cone_cost(X, m, g1...) < _cone_cost(X, m, g...))
            g = g1
        end
        g === nothing && continue
        r = _cone_refine(X, m, ext, g...)
        r === nothing && continue
        r[2] < bestc && ((best, bestc) = r)
    end
    return best
end

function _cone_guess(X, m, a, quadratic::Bool)
    e1, e2 = _perp_basis(a)
    M = zeros(SMatrix{5,5,Float64})
    g = zeros(SVector{5,Float64})
    for x in X
        y = x - m
        u, v, z = dot(y, e1), dot(y, e2), dot(y, a)
        q = SVector(2u, 2v, 1.0, z, quadratic ? z^2 : 0.0)
        M += q * q'
        g += q * (u^2 + v^2)
    end
    quadratic || (M += SMatrix{5,5,Float64}(Diagonal(SVector(0.0, 0.0, 0.0, 0.0, 1.0))))
    s = _solve_small(M, g)
    s === nothing && return nothing
    R2 = s[3] + s[1]^2 + s[2]^2
    R2 > 0 || return nothing
    R = sqrt(R2)
    return a, s[1] * e1 + s[2] * e2, R, s[4] / 2R
end

function _cone_cost(X, m, a, o, R, k)
    c = 0.0
    for x in X
        d = x - m - o
        z = dot(d, a)
        c += (norm(d - z * a) - R - k * z)^2
    end
    return c / (1 + k^2)
end

function _cone_refine(X, m, ext, a, o, R, k)
    cost = _cone_cost(X, m, a, o, R, k)
    for _ in 1:30
        (isfinite(R) && isfinite(k) && R > 0) || return nothing
        e1, e2 = _perp_basis(a)
        S = sqrt(1 + k^2)
        M6 = zeros(SMatrix{6,6,Float64})
        g6 = zeros(SVector{6,Float64})
        for x in X
            d = x - m - o
            z = dot(d, a)
            q = d - z * a
            ρ = norm(q)
            ρ > 0 || continue
            û = q / ρ
            ue1, ue2 = dot(û, e1), dot(û, e2)
            r = (ρ - R - k * z) / S
            h = (z + k * ρ) / S
            J = SVector(-ue1 * h, -ue2 * h, -ue1 / S, -ue2 / S, -1 / S, -z / S - r * k / S^2)
            M6 += J * J'
            g6 += J * r
        end
        δ = _solve_small(M6, -g6)
        δ === nothing && return nothing
        step = 1.0
        accepted = false
        a2, o2, R2, k2, c2 = a, o, R, k, cost
        for _ in 1:12
            a2 = normalize(a + step * δ[1] * e1 + step * δ[2] * e2)
            o2 = o + step * δ[3] * e1 + step * δ[4] * e2
            R2, k2 = R + step * δ[5], k + step * δ[6]
            h = dot(o2, a2)                             # back to the axis point nearest m
            o2 -= h * a2
            R2 -= k2 * h
            c2 = _cone_cost(X, m, a2, o2, R2, k2)
            c2 <= cost && (accepted = true; break)
            step /= 2
        end
        accepted || break
        a, o, R, k, cost = a2, o2, R2, k2, c2
        step * norm(δ) <= 1e-13 * max(R, 1.0) && break
    end
    (isfinite(R) && isfinite(k) && R > 0) || return nothing
    R * sqrt(1 + k^2) < 1e3ext || return nothing     # a plane in disguise
    abs(k) * 1e3ext > R || return nothing             # apex beyond 1000 extents: a cylinder
    return _SurfaceModel(:cone, m + o, a, R, k), cost
end

# Torus: candidate axes from the patch normals, each with the axis position that the
# normals' lines come nearest to meeting, a circle fit in the meridian half-plane
# ((ρ − R)² + (z − z₀)² = r²: ρ² + z² linear in R, z₀ and r² − R² − z₀²), and
# Gauss–Newton (step halving when the residual grows) on the geometric distances
# √((ρ − R)² + (z − z₀)²) − r over the axis direction, its position, z₀, R and r; the
# best is kept. On a surface of revolution every normal lies in its meridian plane, so
# the normal lines (patch centroid x, normal n) all meet the axis (a, p₀): with
# ā = p₀ × a, (x × n)·a + n·ā = 0 (a linear complex). Minimizing the weighted squares
# over ā for each a leaves a 3 × 3 form in a whose two least eigenvectors are
# candidates (one is the axis; on a small patch the other can be, its normals also
# nearly meeting the patch's own normal line), with the cone's two (the axis of the
# circle the normals trace, the direction they are perpendicular to). A torus whose
# major or minor radius exceeds 1000 extents is left to the cylinder or the cone.
function _fit_torus(X, N, W, Cn, m, ext)
    Cn === nothing && return nothing
    A = zeros(SMatrix{3,3,Float64})
    B = zeros(SMatrix{3,3,Float64})
    Cm = zeros(SMatrix{3,3,Float64})
    wt = 0.0
    nbar = zero(Vec3)
    for (x, n, w) in zip(Cn, N, W)
        l = cross(x - m, n)
        A += w * (l * l')
        B += w * (l * n')
        Cm += w * (n * n')
        wt += w
        nbar += w * n
    end
    wt > 0 || return nothing
    nbar /= wt
    Sc = zeros(SMatrix{3,3,Float64})
    for (n, w) in zip(N, W)
        Sc += w * ((n - nbar) * (n - nbar)')
    end
    Ci = inv(Cm + (1e-10 * max(tr(Cm), eps())) * I)
    S = A - B * Ci * B'
    E = eigen(Symmetric(Matrix(S + S') / 2))
    # starts (cost, axis, axis point relative to m, z₀, R, r)
    starts = Tuple{Float64,Vec3,Vec3,Float64,Float64,Float64}[]
    for a in (Vec3(E.vectors[:, 1]), Vec3(E.vectors[:, 2]), Vec3(eigen(Symmetric(Matrix(Sc))).vectors[:, 1]),
              Vec3(eigen(Symmetric(Matrix(Cm))).vectors[:, 1]))
        ā = -(Ci * (B' * a))
        o = cross(a, Vec3(ā))                         # the axis point nearest m
        g = _torus_guess(X, m, a, o)
        g === nothing && continue
        push!(starts, (_torus_cost(X, m, a, o, g...), a, o, g...))
    end
    _torus_osculating!(starts, X, N, W, Cn, m, ext)
    sort!(starts; by = first)
    # the three best distinct starts (major radius 2% or axis 0.02 rad apart: on a small
    # patch neighbouring starts of the scan cost nearly the same and lead to the same
    # local minimum), until one fits to round-off
    best, bestc = nothing, Inf
    tried = 0
    exact = (1e-12 * ext)^2 * length(X)
    for (k, (c0, a, o, z0, R, r)) in enumerate(starts)
        (tried >= 3 || bestc <= exact) && break
        any(j -> abs(starts[j][5] - R) <= 0.02R && abs(dot(starts[j][2], a)) >= cos(0.02), 1:k-1) && continue
        tried += 1
        res = _torus_refine(X, m, ext, a, o, z0, R, r)
        res === nothing && continue
        res[2] < bestc && ((best, bestc) = res)
    end
    return best
end

# Starts from the patch's curvature: the shape operator S (dn = −S dx in the tangent
# plane at the weighted mean centroid x₀ and normal n₀, a least-squares fit over the
# patches) gives principal curvatures κ₁, κ₂ and directions; with t₁ across the tube, the
# spine passes through x₀ + n₀/κ₁ along t₂, and the axis lies in the plane of n₀ and t₁
# and meets the normal line at x₀ + n₀/κ₂. That leaves the axis's tilt in that plane,
# set by third-order terms: it is scanned by where the axis crosses the line x₀ + s t₁,
# s = ±2ᵏ/4 extents, k = 0 … 12, for either principal direction across the tube.
function _torus_osculating!(starts, X, N, W, Cn, m, ext)
    n1 = N[argmax(W)]
    wt = 0.0
    x0, n0 = zero(Vec3), zero(Vec3)
    for (x, n, w) in zip(Cn, N, W)
        x0 += w * x
        n0 += w * (dot(n, n1) < 0 ? -n : n)
        wt += w
    end
    (wt > 0 && norm(n0) > 0) || return starts
    x0 /= wt
    n0 = normalize(n0)
    t1, t2 = _perp_basis(n0)
    M = zeros(SMatrix{3,3,Float64})
    g = zeros(SVector{3,Float64})
    for (x, n, w) in zip(Cn, N, W)
        nn = dot(n, n0) < 0 ? -n : n
        dx, dn = x - x0, nn - n0
        u, v = dot(dx, t1), dot(dx, t2)
        q1, q2 = SVector(u, v, 0.0), SVector(0.0, u, v)
        M += w * (q1 * q1' + q2 * q2')
        g -= w * (q1 * dot(dn, t1) + q2 * dot(dn, t2))
    end
    s = _solve_small(M, g)
    s === nothing && return starts
    F = eigen(Symmetric([s[1] s[2]; s[2] s[3]]))
    far = 1e6 * ext
    for i in 1:2
        κ1, κ2 = F.values[i], F.values[3-i]
        abs(κ1) * far > 1 || continue
        d1 = F.vectors[1, i] * t1 + F.vectors[2, i] * t2
        r = 1 / abs(κ1)
        F1 = x0 + n0 / κ1                             # on the spine
        F2 = x0 + n0 * (abs(κ2) * far > 1 ? 1 / κ2 : copysign(far, κ2 + eps()))
        for k in 0:12, sg in (-1, 1)
            P = x0 + (sg * ext * 2.0^k / 4) * d1
            a = normalize(F2 - P)
            c = P + dot(F1 - P, a) * a                 # spine center
            R = norm(F1 - c)
            R > 0 || continue
            o = (P - m) - dot(P - m, a) * a
            z0 = dot(c - m - o, a)
            push!(starts, (_torus_cost(X, m, a, o, z0, R, r), a, o, z0, R, r))
        end
    end
    return starts
end

function _torus_guess(X, m, a, o)
    M = zeros(SMatrix{3,3,Float64})
    g = zeros(SVector{3,Float64})
    for x in X
        d = x - m - o
        z = dot(d, a)
        ρ = norm(d - z * a)
        q = SVector(2ρ, 2z, 1.0)
        M += q * q'
        g += q * (ρ^2 + z^2)
    end
    s = _solve_small(M, g)
    s === nothing && return nothing
    R, z0 = s[1], s[2]
    r2 = s[3] + R^2 + z0^2
    (R > 0 && r2 > 0) || return nothing
    return z0, R, sqrt(r2)
end

function _torus_cost(X, m, a, o, z0, R, r)
    c = 0.0
    for x in X
        d = x - m - o
        z = dot(d, a)
        w, h = norm(d - z * a) - R, z - z0
        c += (sqrt(w^2 + h^2) - r)^2
    end
    return c
end

function _torus_refine(X, m, ext, a, o, z0, R, r)
    cost = _torus_cost(X, m, a, o, z0, R, r)
    λ = 1e-6
    for _ in 1:60
        (isfinite(R) && isfinite(r) && isfinite(z0) && R > 0 && r > 0) || return nothing
        e1, e2 = _perp_basis(a)
        M7 = zeros(SMatrix{7,7,Float64})
        g7 = zeros(SVector{7,Float64})
        for x in X
            d = x - m - o
            z = dot(d, a)
            q = d - z * a
            ρ = norm(q)
            ρ > 0 || continue
            û = q / ρ
            w, h = ρ - R, z - z0
            s = sqrt(w^2 + h^2)
            s > 0 || continue
            ŵ, ĥ = w / s, h / s
            ue1, ue2 = dot(û, e1), dot(û, e2)
            t = ĥ * ρ - ŵ * z
            J = SVector(ue1 * t, ue2 * t, -ŵ * ue1, -ŵ * ue2, -ĥ, -ŵ, -1.0)
            M7 += J * J'
            g7 += J * (s - r)
        end
        # Levenberg–Marquardt: on a small patch the torus is barely determined (its
        # axis and major radius show only in third-order terms), so plain Gauss–Newton
        # steps are huge and step halving stalls
        dg = SVector(M7[1, 1], M7[2, 2], M7[3, 3], M7[4, 4], M7[5, 5], M7[6, 6], M7[7, 7])
        accepted = false
        a2, o2, z02, R2, r2, c2 = a, o, z0, R, r, cost
        δ = zeros(SVector{7,Float64})
        while λ <= 1e12
            δ = _solve_small(M7 + SMatrix{7,7,Float64}(Diagonal(λ * dg .+ 1e-300)), -g7)
            if δ !== nothing
                a2 = normalize(a + δ[1] * e1 + δ[2] * e2)
                o2 = o + δ[3] * e1 + δ[4] * e2
                z02, R2, r2 = z0 + δ[5], R + δ[6], r + δ[7]
                hp = dot(o2, a2)                        # back to the axis point nearest m
                o2 -= hp * a2
                z02 += hp
                c2 = R2 > 0 && r2 > 0 ? _torus_cost(X, m, a2, o2, z02, R2, r2) : Inf
                c2 <= cost && (accepted = true; break)
            end
            λ *= 10
        end
        accepted || break
        λ = max(λ / 10, 1e-12)
        done = cost - c2 <= 1e-12 * cost || norm(δ) <= 1e-13 * max(R, 1.0)
        a, o, z0, R, r, cost = a2, o2, z02, R2, r2, c2
        done && break
    end
    (isfinite(R) && isfinite(r) && R > 0 && r > 0) || return nothing
    (R < 1e3ext && r < 1e3ext) || return nothing
    return _SurfaceModel(:torus, m + o + z0 * a, a, R, r), cost
end

function _perp_basis(a)
    e1 = normalize(cross(a, abs(a[1]) < 0.9 ? Vec3(1, 0, 0) : Vec3(0, 1, 0)))
    return e1, cross(a, e1)
end

function _solve_small(M, g)
    F = lu(M; check = false)
    issuccess(F) || return nothing
    x = F \ g
    return all(isfinite, x) ? x : nothing
end

# Per-patch data on welded vertices: positions, each patch's vertices, area, unit
# normal, centroid, and the patches at each vertex. (Neighbors are visited through the
# vertices, each vertex once per pass: a vertex shared by d patches, the apex of a fan
# of CAD slivers, otherwise costs d² per pass; a fill tube with fans of 800 made a capsule
# part's fit two minutes long.)
function _fit_data(m::Mesh, weld)
    id = _welded_ids(m; tol = weld)
    P = Vector{Vec3}(undef, isempty(id) ? 0 : maximum(id))
    for i in reverse(eachindex(id))
        P[id[i]] = m.vertices[i]
    end
    pv = [unique(id[v] for k in r for v in m.triangles[k]) for r in m.patches]
    av = [_patch_area_vector(m, r) for r in m.patches]
    # area as the sum of the facet areas: a patch folded across a curved surface has a
    # shorter area vector
    area = [sum(norm(_tri_cross(triangle_vertices(m, k)...)) for k in r; init = 0.0) / 2 for r in m.patches]
    nrm = [norm(a) > 0 ? a / norm(a) : zero(Vec3) for a in av]
    cen = [_patch_centroid(m, r) for r in m.patches]
    vp = [Int[] for _ in eachindex(P)]
    for (p, vs) in enumerate(pv), v in vs
        push!(vp[v], p)
    end
    return (; P, pv, area, nrm, cen, vp)
end

# Append to `out` the patches sharing a vertex with the patches `ps` that pass `keep`, in
# the order of first occurrence over p ∈ ps, v ∈ pv[p], q ∈ vp[v], each once, skipping
# those already stamped: patches with `pstamp[q] == s` (stamped as they are added) and
# vertices with `vstamp[v] == s` (stamped as they are visited).
function _append_neighbors!(out, D, ps, keep, pstamp, vstamp, s)
    for p in ps, v in D.pv[p]
        vstamp[v] == s && continue
        vstamp[v] = s
        for q in D.vp[v]
            pstamp[q] == s && continue
            keep(q) || continue
            pstamp[q] = s
            push!(out, q)
        end
    end
    return out
end

function _fit_patches(kind, D, ps, scratch)
    X = Vec3[]
    seen = scratch
    for p in ps, v in D.pv[p]
        seen[v] && continue
        seen[v] = true
        push!(X, D.P[v])
    end
    for p in ps, v in D.pv[p]
        seen[v] = false
    end
    isempty(X) && return nothing
    return _fit_model(kind, X, D.nrm[ps], D.area[ps], kind === :torus ? D.cen[ps] : nothing)
end

# Tolerance on the distance of vertices from a model: `rtol` of the radius or of the
# region's extent `ext` if smaller (a gently curved face is not a sphere of any radius
# that bends less than `rtol R` across it), at least `atol`; `atol` for a plane (CAD
# planes are flat to round-off) and a torus (CAD tori, rolling-ball fillets, have their
# vertices on them to round-off, while within `rtol` of a large one lie strips of the
# cone and the cylinder it blends).
_fit_tol(m::_SurfaceModel, ext, rtol, atol) =
    m.kind === :plane || m.kind === :torus ? atol : max(atol, rtol * min(_model_radius(m), ext))

# Diagonal of the bounding box of the patches' vertices.
function _extent(D, ps)
    lo, hi = Vec3(Inf, Inf, Inf), Vec3(-Inf, -Inf, -Inf)
    for p in ps, v in D.pv[p]
        lo, hi = min.(lo, D.P[v]), max.(hi, D.P[v])
    end
    return norm(hi - lo)
end

function _patch_on(m::_SurfaceModel, D, p, tol, cosmax)
    for v in D.pv[p]
        abs(_model_distance(m, D.P[v])) <= tol || return false
    end
    D.area[p] > 0 || return true
    return abs(dot(D.nrm[p], _model_normal(m, D.cen[p]))) >= cosmax
end

# Grow a region of one kind from `seed` over available patches, nearest (to the seed)
# first, refitting as it grows; `nothing` if no start around the seed fits. `refs` are
# the rms trajectories of rejected growths covering the seed: the growth is then given
# up at a checkpoint where it fits no better than they did (see `fit_surfaces`).
function _grow_region(kind, seed, D, avail, W, opts, refs = nothing)
    (; rtol, atol, cosmax) = opts
    (; mark, queued, scratch, heap) = W
    W.gid[] += 1
    gid = W.gid[]
    # start from the seed and its neighbors (two rings, then one) facing its way
    near(q) = avail[q] && (D.area[q] == 0 || abs(dot(D.nrm[q], D.nrm[seed])) >= cosmax)
    (; pring, vring) = W
    s = W.stamp[] += 1
    pring[seed] = s
    ring1 = _append_neighbors!([seed], D, (seed,), near, pring, vring, s)
    s = W.stamp[] += 1
    pring[ring1] .= s
    ring2 = _append_neighbors!(copy(ring1), D, ring1, near, pring, vring, s)
    model = nothing
    members = Int[]
    for init in (ring2, ring1, [seed])
        # (a torus starts from at most 64 patches: rings around the apex of a fan of CAD
        # slivers hold hundreds, and the torus fit is the costly one)
        kind === :torus && length(init) > 64 && continue
        mdl = _fit_patches(kind, D, init, scratch)
        mdl === nothing && continue
        tol = _fit_tol(mdl, _extent(D, init), rtol, atol)
        all(p -> _patch_on(mdl, D, p, tol, cosmax), init) || continue
        model, members = mdl, copy(init)
        break
    end
    model === nothing && return nothing
    W.growths[] += 1
    return _grow_from!(kind, seed, model::_SurfaceModel, members, D, avail, W, gid, opts, refs)
end

# Queue the available patches sharing a vertex with p, not yet members or queued (each
# vertex is visited once per growth: its patches are queued, members or unavailable from
# then on).
function _push_neighbors!(H, D, p, avail, mark, queued, vgrow, gid, c0)
    for v in D.pv[p]
        vgrow[v] == gid && continue
        vgrow[v] = gid
        for q in D.vp[v]
            (avail[q] && mark[q] != gid && queued[q] != gid) || continue
            queued[q] = gid
            _heap_push!(H, norm(D.cen[q] - c0), q, 0)
        end
    end
end

function _grow_from!(kind, seed, model::_SurfaceModel, members::Vector{Int}, D, avail, W, gid, opts,
                     refs = nothing)
    (; rtol, atol, cosmax) = opts
    (; mark, queued, scratch, heap) = W
    H = heap
    for p in members
        mark[p] = gid
    end
    empty!(H.key); empty!(H.tri); empty!(H.seed)
    c0 = D.cen[seed]
    for p in members
        _push_neighbors!(H, D, p, avail, mark, queued, W.vgrow, gid, c0)
    end
    lo, hi = Vec3(Inf, Inf, Inf), Vec3(-Inf, -Inf, -Inf)
    for p in members, v in D.pv[p]
        lo, hi = min.(lo, D.P[v]), max.(hi, D.P[v])
    end
    tol = _fit_tol(model, norm(hi - lo), rtol, atol)
    nfit = length(members)
    rejected = Int[]
    rounds = 0
    # checkpoints at W.check × 2^k members (k = 0, 1, …): the rms over the members of the
    # current fit, in units of its tolerance (read only: the growth itself is unchanged)
    traj = W.traj
    empty!(traj)
    check = W.check
    while true
        while !isempty(H.key)
            _, p, _ = _heap_pop!(H)
            mark[p] == gid && continue
            W.work[] += 1
            if _patch_on(model, D, p, tol, cosmax)
                mark[p] = gid
                push!(members, p)
                for v in D.pv[p]
                    lo, hi = min.(lo, D.P[v]), max.(hi, D.P[v])
                end
                if length(members) >= check
                    q = _deviation(model, D, members, scratch)[1] / tol
                    while length(members) >= check
                        push!(traj, q)
                        check *= 2
                    end
                    if refs !== nothing && q > _cover_rms(refs, length(traj)) / 2
                        W.givenup[] += 1
                        return nothing
                    end
                    # a torus retracing a sphere, cylinder, cone or plane
                    kind === :torus && (_torus_degenerate(model, norm(hi - lo), tol) ||
                                        _simpler_fits(D, members, opts, scratch; exact = true)) &&
                        (W.degen[] = members; return nothing)
                end
                _push_neighbors!(H, D, p, avail, mark, queued, W.vgrow, gid, c0)
                length(members) < 1.25nfit + 4 && continue
            else
                push!(rejected, p)
                continue
            end
            # refit, and reconsider the rejected patches under the new fit
            mdl = _fit_patches(kind, D, members, scratch)
            if mdl !== nothing
                model = mdl
                tol = _fit_tol(model, norm(hi - lo), rtol, atol)
            end
            nfit = length(members)
            for q in rejected
                _heap_push!(H, norm(D.cen[q] - c0), q, 0)
            end
            empty!(rejected)
        end
        rounds += 1
        (isempty(rejected) || nfit == length(members) || rounds > 3) && break
        mdl = _fit_patches(kind, D, members, scratch)
        mdl === nothing || (model = mdl; tol = _fit_tol(model, norm(hi - lo), rtol, atol))
        nfit = length(members)
        for q in rejected
            _heap_push!(H, norm(D.cen[q] - c0), q, 0)
        end
        empty!(rejected)
    end
    # final fit; drop patches that it leaves out of tolerance
    for _ in 1:2
        mdl = _fit_patches(kind, D, members, scratch)
        mdl === nothing && return nothing
        model = mdl
        tol = _fit_tol(model, _extent(D, members), rtol, atol)
        keep = filter(p -> _patch_on(model, D, p, tol, cosmax), members)
        length(keep) == length(members) && break
        for p in members
            mark[p] = 0
        end
        members = keep
        isempty(members) && return nothing
    end
    return members, model
end

# rms and largest distance of the vertices of patches `members` from the model
function _deviation(model, D, members, scratch)
    s2, mx, n = 0.0, 0.0, 0
    for p in members, v in D.pv[p]
        scratch[v] && continue
        scratch[v] = true
        d = abs(_model_distance(model, D.P[v]))
        s2 += d^2
        mx = max(mx, d)
        n += 1
    end
    for p in members, v in D.pv[p]
        scratch[v] = false
    end
    return sqrt(s2 / max(n, 1)), mx
end

function _region(kind, members, model, D, scratch)
    A = sum(D.area[members])
    cen = A > 0 ? sum(D.area[p] * D.cen[p] for p in members) / A : sum(D.cen[members]) / length(members)
    rms, mx = _deviation(model, D, members, scratch)
    sides = [sign(dot(D.nrm[p], _model_normal(model, D.cen[p]))) for p in members if D.area[p] > 0]
    nrm = zero(Vec3)
    if kind === :plane
        w = sum(D.area[p] * dot(D.nrm[p], model.u) for p in members)
        nrm = w < 0 ? -model.u : model.u
        sides .*= w < 0 ? -1 : 1
        facing = all(>(0), sides) ? 1 : 0
    else
        facing = all(>(0), sides) ? 1 : all(<(0), sides) ? -1 : 0
    end
    return SurfaceRegion(kind, _model_surface(model), cen, nrm, sort(members), A, rms, mx, facing)
end

"""
    fit_surfaces(mesh; reference = nothing, kinds = (:plane, :cylinder, :cone, :sphere, :torus),
                 rtol = 1e-3, atol = 1e-4, angle = 30.0, minpatches = 6, minarea = 0.0,
                 tol = 1e-9) -> Vector{SurfaceRegion}

Regions of `mesh` whose patches lie on one plane, cylinder (any axis), cone, sphere or
torus, by decreasing area; patches in no region are on none of them (free-form faces,
patches straddling two faces). Use it on imported CAD (raw facets or after
`remesh_mesh`), then [`curved_components`](@ref) to give the curved regions their exact
surface for the radiation model.

A patch belongs to a surface when every one of its vertices lies within
`max(atol, rtol min(R, L))` mm of it (R the radius of curvature, R/cos θ on a cone of
half-angle θ; L the region's extent; `atol` for a plane: CAD planes are flat to
round-off) and its normal is within `angle` degrees of the surface normal at its
centroid (either side: windings do not matter). Regions grow from seed patches, by
decreasing area, over patches that share a vertex (welded by position to `tol` times the
mesh's extent), nearest first, refitting by least squares as they grow: planes by
principal components, spheres by an algebraic fit, cylinders by the axis the patch
normals are perpendicular to and a circle fit across it, cones by the axis the normals
make a constant angle with (or the cylinder's, for a small patch) and an algebraic fit of
the radius linear along it, all refined by Gauss–Newton on the vertex distances; tori
by the axis the normals all meet (or, on a small patch, axes through the patch's
centres of curvature) and a circle fit in the meridian plane, refined by
Levenberg–Marquardt. A cone whose apex lies over 1000 extents away is left to the
cylinder, as a surface whose radius of curvature does is left to the plane. From each
seed all `kinds` are grown and the largest region kept, the earlier of plane, cylinder,
cone, sphere, torus when within 1% of the largest: a
narrow band with vertices only on its two rims fits a sphere as well, and a cone too
when its rims differ in radius (the bore of an entrance hole is a cylinder, a conical
shield tessellated as one band a cone: a tessellator puts rings of vertices inside a
spherical zone that bulges beyond its chord tolerance). Only regions with at least
`minpatches` patches and `minarea` mm² whose rms distance is under a quarter of the
tolerance count: on a true surface the residual is the faceting, well inside the
tolerance, while a fit to a free-form face (a variable-radius fillet; a torus without
`:torus`) grows until the tolerance stops it (rms about 0.4 of it), so such patches
stay unfitted. Regions of one kind whose
patches all lie on each other's fits (one surface interrupted by another face) are then
merged.

Tori (CAD fillets: a quarter of the tube where a wall meets a cap) are held to CAD
round-off, since with three more parameters than a sphere a torus also follows free-form
faces and spans neighbouring faces within the relative tolerance: their tolerance is
`atol`, their rms must be under `atol/16`, their vertices must fall on at least four
evenly spaced points of the tube circle (bands of strips with vertices only on their
rims lie on a circle through any three of them), and patches that a plane, cylinder,
cone or sphere fits are left to it. A fillet then becomes one region in place of a
patchwork of small cones and spheres. Pass `kinds` without `:torus` for the fits
before tori (2026-10-24); on the NIF CAD parts that changes only fillet regions
(research notes, "Tori").

On a finely faceted face of none of the `kinds` (a tapered tube without `:cone`, a fine
fillet without `:torus`) every patch would seed such a growth, at a cost quadratic in
the facet count. A rejected growth of at least 512 patches therefore covers the patches near its seed
(within half its radius, while all its patches stay available): a growth of the same
kind from a covered seed is given up once, at 128, 256, … members, its rms reaches
half that of the covering growth at the same size, as it does when it is retracing the
same face. A true surface fits its patches to the faceting and keeps growing.

The tolerance must exceed the faceting several times over. Facet vertices of exported
CAD lie on the true surface, but points `remesh_mesh` adds on a facet lie up to w²/8R
inside a curved one (w the facet's width: 0.1 µm for 50 µm strips on a 2.9 mm wall,
5 µm for 0.3 mm strips), more than the tolerance allows once w exceeds ~0.08 R. So
give the facets the mesh was remeshed from as `reference`: the regions are then fitted
on them, and each patch of `mesh` joins the region whose surface all its vertices lie
within that tolerance plus the largest distance of the region's facets from it (and
whose bounding box holds them), with its normal within `angle`. The regions returned
index the patches of `mesh`, with their rms and largest distances from the
(reference) fit.
"""
function fit_surfaces(m::Mesh; reference::Union{Nothing,Mesh} = nothing,
                      kinds = (:plane, :cylinder, :cone, :sphere, :torus), rtol::Real = 1e-3, atol::Real = 1e-4,
                      angle::Real = 30.0, minpatches::Integer = 6, minarea::Real = 0.0, tol::Real = 1e-9)
    all(in((:plane, :cylinder, :cone, :sphere, :torus)), kinds) ||
        throw(ArgumentError("kinds must be among :plane, :cylinder, :cone, :sphere and :torus"))
    opts = (; rtol = Float64(rtol), atol = Float64(atol), cosmax = cosd(angle))
    if reference !== nothing
        regs = fit_surfaces(reference; kinds, rtol, atol, angle, minpatches, minarea, tol)
        return _transfer_regions(regs, reference, m, opts, tol)
    end
    return _fit_surfaces(m, kinds, opts, minpatches, minarea, tol)[1]
end

# Covers (see `fit_surfaces`): a growth that the rms test rejects with at least
# `_COVER_MIN` patches covers, for its kind, its patches nearer its seed than half the
# region's radius, as long as all of its patches stay available. A growth of that kind
# from a covered seed is given up at the first checkpoint (`_COVER_CHECK` × 2^k
# members) where its rms exceeds half that of each covering growth at the same
# checkpoint. Measured on the NIF CAD parts (research notes, "Cost of fitting faces of
# none of the kinds"): smaller covers change regions that only just pass the rms test.
const _COVER_MIN = 512
const _COVER_CHECK = 128

# Tori: a torus has three more parameters than a sphere or a cylinder and bends both
# ways, so it follows a free-form face (a paraboloid's patches fit tori with rms 0.23 of
# the tolerance) as closely as the other kinds follow fillets, and within the relative
# tolerance of a large torus (r ~ 10 mm) a cone and part of a cylinder beside it; a CAD
# torus (a rolling-ball fillet) has its vertices on it to round-off. Tori are therefore
# accepted at rms under `_TORUS_RMS` times `atol` (6 nm by default; Float32 vertices of
# a few-mm part are rounded by ~0.3 nm), and not on patches that a plane, cylinder,
# cone or sphere fits (a flat ring is the top band of a torus). A torus within twice the
# tolerance of a simpler surface across a region of extent L (a sphere: R ≤ 2 tol; a
# cylinder: spine sagitta L²/8R ≤ 2 tol; a cone or a plane: meridian sagitta
# L²/8r ≤ 2 tol) is left to that kind, and a torus growth that becomes one at a
# checkpoint is given up (it would otherwise retrace every sphere and cylinder of a
# part).
const _TORUS_RMS = 1 / 16
_torus_degenerate(m::_SurfaceModel, L, tol) = m.R <= 2tol || L^2 / 8m.R <= 2tol || L^2 / 8m.k <= 2tol

# Whether the vertices of patches `ps` show the torus: their positions along the tube
# circle (in its meridian half-plane) fall into at least four distinct points (more than
# 10 tol apart) spaced evenly, the largest step under three times the median, past the
# arc they leave empty. Through three points passes a circle, so bands of strips with
# vertices only on their rims (a cylinder and a cone beside it, tessellated as CAD does)
# fit a torus exactly when they have three rims between them, and four rims placed
# symmetrically (two bands at either end of a cylinder) lie on one circle too; a CAD
# torus is tessellated in bands of equal angle.
function _meridian_ok(m::_SurfaceModel, D, ps, tol)
    s = Float64[]
    for p in ps, v in D.pv[p]
        d = D.P[v] - m.c
        z = dot(d, m.u)
        push!(s, m.k * atan(z, norm(d - z * m.u) - m.R))
    end
    sort!(s)
    g = [s[i+1] - s[i] for i in 1:length(s)-1 if s[i+1] - s[i] > 10tol]
    w = 2π * m.k - (s[end] - s[1])                   # the step across ψ = ±π
    w > 10tol && push!(g, w)
    length(g) >= 4 || return false                   # four points on the circle: four steps
    deleteat!(g, argmax(g))                          # the arc the points leave empty
    sort!(g)
    return g[end] <= 3g[(length(g) + 1) ÷ 2]
end

# Whether the bounding boxes of the vertices of patches `ps` and `qs` overlap (to `tol`).
function _boxes_touch(D, ps, qs, tol)
    function box(ps)
        lo, hi = Vec3(Inf, Inf, Inf), Vec3(-Inf, -Inf, -Inf)
        for p in ps, v in D.pv[p]
            lo, hi = min.(lo, D.P[v]), max.(hi, D.P[v])
        end
        return lo, hi
    end
    (l1, h1), (l2, h2) = box(ps), box(qs)
    return all(l1 .<= h2 .+ tol) && all(l2 .<= h1 .+ tol)
end

# Whether a plane, cylinder, cone or sphere fits all the patches `ps` as a region (every
# patch on it, rms under a quarter of the tolerance; with `exact`, under the rms tori
# are held to).
function _simpler_fits(D, ps, opts, scratch; exact::Bool = false)
    L = _extent(D, ps)
    for kind in (:plane, :cylinder, :cone, :sphere)
        mdl = _fit_patches(kind, D, ps, scratch)
        mdl === nothing && continue
        tol = _fit_tol(mdl, L, opts.rtol, opts.atol)
        all(p -> _patch_on(mdl, D, p, tol, opts.cosmax), ps) || continue
        _deviation(mdl, D, ps, scratch)[1] <= (exact ? _TORUS_RMS * opts.atol : tol / 4) && return true
    end
    return false
end

# A seed whose torus region is rejected marks the patches of the region within half its
# radius (all of a growth given up as a simpler surface at a checkpoint, and the seed
# alone when no torus starts there): no torus is grown from them. On a sphere, or on a
# part's planes with no `:plane` among the kinds, every patch would otherwise grow one
# to its first checkpoint.
function _no_torus!(notorus, D, s, ps)
    c0 = D.cen[s]
    r = maximum(p -> norm(D.cen[p] - c0), ps)
    for p in ps
        norm(D.cen[p] - c0) <= r / 2 && (notorus[p] = true)
    end
    notorus[s] = true
    return notorus
end

# the largest rms (in units of the tolerance) of the covering growths at checkpoint k
# (their final rms past their last checkpoint)
_cover_rms(refs, k) = maximum(r -> r[min(k, end)], refs)

# `fit_surfaces` on the mesh itself; also returns the work done (`growths`: growths
# started, `work`: patches tested while growing, `givenup`: growths given up under a
# cover, `covers`) for the tests of its cost.
function _fit_surfaces(m::Mesh, kinds, opts, minpatches, minarea, tol; covermin = _COVER_MIN, check = _COVER_CHECK)
    np = npatches(m)
    np == 0 && return SurfaceRegion[], (; growths = 0, work = 0, givenup = 0, covers = 0)
    D = _fit_data(m, tol)
    nv = length(D.P)
    W = (; mark = zeros(Int, np), queued = zeros(Int, np), gid = Ref(0), heap = _GrowHeap(),
         scratch = falses(nv), vgrow = zeros(Int, nv), pring = zeros(Int, np), vring = zeros(Int, nv),
         stamp = Ref(0), growths = Ref(0), work = Ref(0), givenup = Ref(0), traj = Float64[], check,
         degen = Ref{Union{Nothing,Vector{Int}}}(nothing))
    avail = trues(np)
    order = sortperm(D.area; rev = true)
    ordered = [k for k in (:plane, :cylinder, :cone, :sphere, :torus) if k in kinds]
    found = Tuple{Symbol,Vector{Int},_SurfaceModel}[]
    # covering growths: kind, rms trajectory and covered patches, whether each still
    # stands, the covering growths of each patch and the growths each patch is in
    # seeds near a torus growth that turned out a simpler surface grow no torus
    notorus = falses(np)
    covers = Tuple{Symbol,Vector{Float64},Vector{Int}}[]
    standing = Bool[]
    coveredby = Dict{Int,Vector{Int}}()
    touched = Dict{Int,Vector{Int}}()
    refs = Vector{Float64}[]
    for s in order
        avail[s] || continue
        best = nothing
        bestA = 0.0
        cs = get(coveredby, s, nothing)
        for kind in ordered
            empty!(refs)
            if cs !== nothing
                for g in cs
                    standing[g] && covers[g][1] === kind && push!(refs, covers[g][2])
                end
            end
            kind === :torus && notorus[s] && continue
            W.degen[] = nothing
            g = _grow_region(kind, s, D, avail, W, opts, isempty(refs) ? nothing : refs)
            if g === nothing
                # (a seed whose torus growth fails grows none again, nor do those near it
                # when the growth found a simpler surface)
                kind === :torus && (W.degen[] === nothing ? (notorus[s] = true) : (notorus[W.degen[]] .= true; notorus[s] = true))
                continue
            end
            ps, mdl = g
            A = sum(D.area[ps])
            length(ps) >= minpatches && A >= minarea || (kind === :torus && _no_torus!(notorus, D, s, ps); continue)
            # a surface explains its vertices well within the tolerance; a fit that only
            # grew until the tolerance stopped it (a sphere on a free-form face) does not
            ext = _extent(D, ps)
            tol = _fit_tol(mdl, ext, opts.rtol, opts.atol)
            q = _deviation(mdl, D, ps, W.scratch)[1] / tol
            if kind === :torus
                (_torus_degenerate(mdl, ext, tol) || !_meridian_ok(mdl, D, ps, tol) ||
                 _simpler_fits(D, ps, opts, W.scratch)) &&
                    (_no_torus!(notorus, D, s, ps); continue)
            end
            if kind === :torus ? q * tol > _TORUS_RMS * opts.atol : q > 1 / 4
                kind === :torus && _no_torus!(notorus, D, s, ps)
                if length(ps) >= covermin
                    c0 = D.cen[s]
                    r = maximum(p -> norm(D.cen[p] - c0), ps)
                    push!(covers, (kind, [W.traj; q], [p for p in ps if norm(D.cen[p] - c0) <= r / 2]))
                    push!(standing, true)
                    for p in covers[end][3]
                        push!(get!(Vector{Int}, coveredby, p), length(covers))
                    end
                    for p in ps
                        push!(get!(Vector{Int}, touched, p), length(covers))
                    end
                end
                continue
            end
            if best === nothing || A > 1.01bestA
                best, bestA = (kind, ps, mdl), A
            end
        end
        best === nothing && continue
        avail[best[2]] .= false
        push!(found, best)
        # a rejected growth over patches now taken no longer stands for its seed's
        for p in best[2]
            gs = pop!(touched, p, nothing)
            gs === nothing || (standing[gs] .= false)
        end
    end
    # merge regions of one kind lying on one surface
    merged = true
    apart = Set{Tuple{UInt,UInt}}()                  # torus pairs whose union was tried
    while merged
        merged = false
        for i in eachindex(found), j in i+1:length(found)
            (ki, psi, mi), (kj, psj, mj) = found[i], found[j]
            ki === kj || continue
            if ki === :torus
                # tori fitted on small pieces of one fillet are off by more than `atol`
                # across it: pieces whose boxes touch are merged when their union fits
                (objectid(psi), objectid(psj)) in apart && continue
                push!(apart, (objectid(psi), objectid(psj)))
                _boxes_touch(D, psi, psj, opts.atol) || continue
            else
                ti = _fit_tol(mi, _extent(D, psi), opts.rtol, opts.atol)
                tj = _fit_tol(mj, _extent(D, psj), opts.rtol, opts.atol)
                all(p -> _patch_on(mi, D, p, ti, opts.cosmax), psj) &&
                    all(p -> _patch_on(mj, D, p, tj, opts.cosmax), psi) || continue
            end
            ps = [psi; psj]
            mdl = _fit_patches(ki, D, ps, W.scratch)
            mdl === nothing && continue
            t = _fit_tol(mdl, _extent(D, ps), opts.rtol, opts.atol)
            all(p -> _patch_on(mdl, D, p, t, opts.cosmax), ps) || continue
            ki === :torus && _deviation(mdl, D, ps, W.scratch)[1] > _TORUS_RMS * opts.atol && continue
            found[i] = (ki, ps, mdl)
            deleteat!(found, j)
            merged = true
            break
        end
    end
    regions = [_region(k, ps, mdl, D, W.scratch) for (k, ps, mdl) in found]
    return sort!(regions; by = r -> -r.area), (; growths = W.growths[], work = W.work[], givenup = W.givenup[], covers = length(covers))
end

function _region_model(r::SurfaceRegion)
    r.kind === :plane && return _SurfaceModel(:plane, r.centroid, r.normal, Inf)
    r.kind === :cylinder && return _SurfaceModel(:cylinder, r.surface.origin, r.surface.axis, r.surface.radius)
    r.kind === :sphere && return _SurfaceModel(:sphere, r.surface.center, zero(Vec3), r.surface.radius)
    r.kind === :torus &&
        return _SurfaceModel(:torus, r.surface.center, r.surface.axis, r.surface.major_radius, r.surface.minor_radius)
    # a cone about the axis point level with the region's centroid, as fitted
    s = r.surface
    z = dot(r.centroid - s.apex, s.axis)
    k = tand(s.half_angle)
    return _SurfaceModel(:cone, s.apex + z * s.axis, s.axis, z * k, k)
end

# Largest distance from the model of the points of the triangle (a, b, c), sampled on a
# barycentric grid of step 1/4 (on a cylinder, cone or sphere the distance is concave
# across a facet whose vertices lie on it, largest inside or on an edge).
function _facet_deviation(m::_SurfaceModel, a, b, c)
    d = 0.0
    for i in 0:4, j in 0:4-i
        d = max(d, abs(_model_distance(m, a + (i / 4) * (b - a) + (j / 4) * (c - a))))
    end
    return d
end

# Regions fitted on `ref` carried over to the patches of `m`, a remeshing of the same
# surface (see `fit_surfaces`).
function _transfer_regions(regs, ref::Mesh, m::Mesh, opts, weld)
    np = npatches(m)
    (np == 0 || isempty(regs)) && return SurfaceRegion[]
    D = _fit_data(m, weld)
    Dr = _fit_data(ref, weld)
    models = [_region_model(r) for r in regs]
    lim = Float64[]
    boxes = Tuple{Vec3,Vec3}[]
    for (r, mdl) in zip(regs, models)
        sag = 0.0
        for p in r.patches, k in ref.patches[p]
            sag = max(sag, _facet_deviation(mdl, triangle_vertices(ref, k)...))
        end
        l = _fit_tol(mdl, _extent(Dr, r.patches), opts.rtol, opts.atol) + 1.25sag
        lo, hi = Vec3(Inf, Inf, Inf), Vec3(-Inf, -Inf, -Inf)
        for p in r.patches, v in Dr.pv[p]
            lo, hi = min.(lo, Dr.P[v]), max.(hi, Dr.P[v])
        end
        push!(lim, l)
        push!(boxes, (lo .- l, hi .+ l))
    end
    members = [Int[] for _ in regs]
    for p in 1:np
        best, score = 0, Inf
        for (k, mdl) in enumerate(models)
            lo, hi = boxes[k]
            all(v -> all(lo .<= D.P[v] .<= hi), D.pv[p]) || continue
            d = maximum(v -> abs(_model_distance(mdl, D.P[v])), D.pv[p])
            d <= lim[k] && d / lim[k] < score || continue
            D.area[p] > 0 && abs(dot(D.nrm[p], _model_normal(mdl, D.cen[p]))) < opts.cosmax && continue
            best, score = k, d / lim[k]
        end
        best > 0 && push!(members[best], p)
    end
    scratch = falses(length(D.P))
    out = [_region(r.kind, ps, mdl, D, scratch) for (r, mdl, ps) in zip(regs, models, members) if !isempty(ps)]
    return sort!(out; by = r -> -r.area)
end

# The patches `ps` of a mesh as a mesh of their own (only the vertices they use).
function _patch_submesh(m::Mesh, ps)
    vmap = zeros(Int, length(m.vertices))
    verts = Vec3[]
    tris = NTuple{3,Int}[]
    patches = UnitRange{Int}[]
    idx(v) = (vmap[v] == 0 && (push!(verts, m.vertices[v]); vmap[v] = length(verts)); vmap[v])
    for p in ps
        t0 = length(tris) + 1
        for k in m.patches[p]
            a, b, c = m.triangles[k]
            push!(tris, (idx(a), idx(b), idx(c)))
        end
        push!(patches, t0:length(tris))
    end
    return Mesh(verts, tris, patches)
end

"""
    curved_components(component; kinds = (:cylinder, :cone, :sphere, :torus), regions = nothing,
                      maxradius = Inf, minfraction = 1e-3, kwargs...) -> Vector{Component}
    curved_components(components; kwargs...) -> Vector{Component}
    curved_components(target; kwargs...) -> Target

The component split into one component per cylinder, cone, sphere or torus that
[`fit_surfaces`](@ref) finds in its mesh (other keywords are passed on to it; or give
its `regions`), each carrying its fitted `CylinderSurface`, `ConeSurface`,
`SphereSurface` or `TorusSurface`, and the rest: `<name>/cyl1`, `<name>/cyl2`, …,
`<name>/cone1`, …, `<name>/sph1`, …, `<name>/tor1`, …, by decreasing area, then `<name>/rest` with the other patches (planes
and unrecognized faces). Material, kind, sidedness and `part` are the component's. A
component that one fitted surface covers whole keeps its name (and its mesh) and only
gains the surface. For a remeshed CAD part pass the facets it was remeshed from as
`reference` (see `fit_surfaces`). A component that already has a surface, or in which
no curved region is found, is returned unchanged.

`maxradius` and `minfraction` leave regions in the rest: those whose radius of
curvature (R/cos θ on a cone of half-angle θ, the minor radius on a torus) exceeds
`maxradius` times the part's extent (its bounding-box diagonal), and those of less than
`minfraction` of the part's area. Gently curved CAD faces that are neither planar nor free-form within the fit's
tolerance (radii of tens to hundreds of mm on a hohlraum half) otherwise each become a
component of their own; their flat elements are close to exact (the sagitta of an
element of size h is h²/8R). `minfraction` defaults to 1e-3: small holes and edge faces
of CAD fit as tiny cylinders, cones and spheres (and, before tori, the fillet pieces),
and on a NIF CAD cavity (remeshed at 0.6 and 0.3 mm) the default keeps 68 of 158–159 components, with
the capsule's incident flux within 3e-4 and its Legendre moments within 0.04 percentage
points of fitting everything (research notes, "Fitting CAD by default"). Pass
`minfraction = 0` to fit every region. `maxradius` is off by default.

The radiation model then integrates over the true curved surface, as it does for the
primitive `hohlraum`: each element's sub-triangle points are projected onto the surface
with exact curved-patch areas, whatever the element's shape (the compact,
many-triangle elements of `remesh_mesh` included), and the facets of a cylinder, cone or
sphere no longer occlude pairs of points inside it (a torus is not convex: its facets
still do). Ray casting keeps the facets.
Planar regions are not given a surface: flat elements are already exact there
(`DiskSurface` describes polar-sector elements of primitive end caps, not general
planar CAD faces).

The pieces remember their part (`part`), so the split is computational only:
`elements(target, :part)`, the other name-taking functions and the reports
(`energy_balance`, `component_power`, `write_results`, `report`) see the part, and
each piece stays reachable by its own name (see [`parts`](@ref)).
"""
function curved_components(c::Component; kinds = (:cylinder, :cone, :sphere, :torus), regions = nothing,
                           maxradius::Real = Inf, minfraction::Real = 1e-3, kwargs...)
    all(in((:cylinder, :cone, :sphere, :torus)), kinds) ||
        throw(ArgumentError("kinds must be among :cylinder, :cone, :sphere and :torus"))
    c.surface === nothing || return [c]
    regs = regions === nothing ? fit_surfaces(c.mesh; kwargs...) : regions
    sel = [r for r in regs if r.kind in kinds && r.surface !== nothing]
    if isfinite(maxradius) || minfraction > 0
        m = c.mesh
        lo, hi = Vec3(Inf, Inf, Inf), Vec3(-Inf, -Inf, -Inf)
        for v in m.vertices
            lo, hi = min.(lo, v), max.(hi, v)
        end
        L = norm(hi - lo)
        A = sum(norm(_tri_cross(triangle_vertices(m, k)...)) for k in eachindex(m.triangles); init = 0.0) / 2
        filter!(r -> _model_radius(_region_model(r)) <= maxradius * L && r.area >= minfraction * A, sel)
    end
    sort!(sel; by = r -> -r.area)
    isempty(sel) && return [c]
    if length(sel) == 1 && length(unique(sel[1].patches)) == npatches(c.mesh)
        # one surface covers the whole part: it keeps its name and its mesh
        return [Component(c.name, c.mesh, c.material; kind = c.kind, two_sided = c.two_sided,
                          surface = sel[1].surface, c.part)]
    end
    out = Component[]
    used = falses(npatches(c.mesh))
    nk = Dict{Symbol,Int}()
    for r in sel
        any(used[r.patches]) && throw(ArgumentError("regions overlap"))
        k = nk[r.kind] = get(nk, r.kind, 0) + 1
        name = Symbol(c.name, r.kind === :cylinder ? "/cyl" : r.kind === :cone ? "/cone" : r.kind === :torus ? "/tor" : "/sph", k)
        push!(out, Component(name, _patch_submesh(c.mesh, r.patches), c.material; kind = c.kind,
                             two_sided = c.two_sided, surface = r.surface, c.part))
        used[r.patches] .= true
    end
    rest = findall(!, used)
    isempty(rest) || push!(out, Component(Symbol(c.name, "/rest"), _patch_submesh(c.mesh, rest), c.material;
                                          kind = c.kind, two_sided = c.two_sided, c.part))
    return out
end

curved_components(cs::AbstractVector{Component}; kwargs...) =
    reduce(vcat, (curved_components(c; kwargs...) for c in cs); init = Component[])
curved_components(t::Target; kwargs...) =
    Target(curved_components(t.components; kwargs...); t.closed, t.visibility)
