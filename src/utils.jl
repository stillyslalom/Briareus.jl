# Small numerical utilities shared across the package: orthonormal bases, cubic Hermite
# interpolation, a binary min-heap, Gauss–Legendre rules, the smooth stagnation of a
# ballistic position and nonnegative least squares.

# Returns two unit vectors that complete the unit vector `n` to a right-handed
# orthonormal basis (e1, e2, n). The first is perpendicular to `n` and to the x axis, or
# to the y axis when `n` is within about 26° of the x axis.
@inline function _ortho_basis(n)
    e1 = normalize(cross(n, abs(n[1]) < 0.9 ? Vec3(1, 0, 0) : Vec3(0, 1, 0)))
    return e1, cross(n, e1)
end

# Cubic Hermite interpolation at `t` on the segment from (t0, r0) with slope d0 to
# (t1, r1) with slope d1.
function _hermite(t, t0, r0, d0, t1, r1, d1)
    L = t1 - t0
    s = (t - t0) / L
    return (2s^3 - 3s^2 + 1) * r0 + (s^3 - 2s^2 + s) * L * d0 + (3s^2 - 2s^3) * r1 + (s^3 - s^2) * L * d1
end

# Cubic Hermite interpolation at `x` in the table `y` on the increasing grid `r`, with the
# slopes `dy` at the grid points. Outside the grid the end segment is extended.
@inline function _hermite(r, y, dy, x)
    i = clamp(searchsortedlast(r, x), 1, length(r) - 1)
    h = r[i+1] - r[i]
    t = (x - r[i]) / h
    t2 = t * t
    t3 = t2 * t
    return (2t3 - 3t2 + 1) * y[i] + (t3 - 2t2 + t) * h * dy[i] + (3t2 - 2t3) * y[i+1] +
           (t3 - t2) * h * dy[i+1]
end

# Three-point slopes on a non-uniform grid (exact for quadratics), one-sided at the ends.
function _hermite_slopes(r, y)
    n = length(r)
    d = similar(y)
    d[1] = (y[2] - y[1]) / (r[2] - r[1])
    d[n] = (y[n] - y[n-1]) / (r[n] - r[n-1])
    for k in 2:n-1
        h0, h1 = r[k] - r[k-1], r[k+1] - r[k]
        d[k] = (h0^2 * (y[k+1] - y[k]) + h1^2 * (y[k] - y[k-1])) / (h0 * h1 * (h0 + h1))
    end
    return d
end

# Binary min-heap on a vector of tuples, ordered by their first entry: `_heap_push!` adds a
# tuple and `_heap_pop!` removes and returns the tuple with the least first entry.
function _heap_push!(h, x)
    push!(h, x)
    i = length(h)
    while i > 1 && h[i>>1][1] > h[i][1]
        h[i>>1], h[i] = h[i], h[i>>1]
        i >>= 1
    end
    return h
end

function _heap_pop!(h)
    top, last = h[1], pop!(h)
    isempty(h) && return top
    h[1] = last
    i, n = 1, length(h)
    while 2i <= n
        c = 2i < n && h[2i+1][1] < h[2i][1] ? 2i + 1 : 2i
        h[i][1] <= h[c][1] && break
        h[i], h[c] = h[c], h[i]
        i = c
    end
    return top
end

# The heap of (key, triangle, seed) that grows elements nearest first.
_grow_heap() = Tuple{Float64,Int,Int}[]

# Gauss–Legendre nodes and weights on [-1, 1] (Newton on the three-term recurrence).
function _gauss_legendre(n)
    x, w = zeros(n), zeros(n)
    for i in 1:n
        z = cos(π * (i - 0.25) / (n + 0.5))
        dp = 0.0
        for _ in 1:100
            p0, p1 = 1.0, z
            for k in 2:n
                p0, p1 = p1, ((2k - 1) * z * p1 - (k - 1) * p0) / k
            end
            dp = n * (z * p1 - p0) / (z^2 - 1)
            δ = p1 / dp
            z -= δ
            abs(δ) < 1e-15 && break
        end
        x[i], w[i] = z, 2 / ((1 - z^2) * dp^2)
    end
    return x, w
end

# Returns the position `x` of a ballistic front that stagnates smoothly at `xmax`. The
# position is `x` up to `xmax - w` and approaches `xmax` exponentially over the width `w`
# beyond it. With `xmax = Inf` the front stays ballistic.
function _stagnate(x, xmax, w)
    isinf(xmax) && return x + zero(xmax)
    k = xmax - w
    return x <= k ? x + zero(k) : xmax - w * exp(-(x - k) / w)
end

# Solve the nonnegative least-squares problem argmin ‖A z - b‖ over z ≥ 0 by the
# active-set method of Lawson and Hanson, "Solving Least Squares Problems",
# Prentice-Hall 1974, ch. 23. Rank-deficient passive sets are solved by pivoted QR, which
# gives a basic solution.
function _nnls(A::Matrix{Float64}, b::Vector{Float64})
    m, n = size(A)
    z = zeros(n)
    n == 0 && return z
    P = falses(n)
    tol = 10 * eps() * max(m, n) * max(opnorm(A, 1), 1e-300) * max(norm(b, Inf), 1e-300)
    w = A' * b
    # Columns whose entry left the solution unchanged are skipped until the solution
    # moves. For such a column the least-squares entry came out nonpositive at once, a
    # gain at the rounding level.
    blocked = falses(n)
    for _ in 1:3n
        j = 0
        for i in 1:n
            !P[i] && !blocked[i] && w[i] > tol && (j == 0 || w[i] > w[j]) && (j = i)
        end
        j == 0 && break
        P[j] = true
        entered = true
        while true
            idx = findall(P)
            t = zeros(n)
            t[idx] .= A[:, idx] \ b
            if all(>(0), view(t, idx))
                z = t
                break
            end
            if entered && t[j] <= 0
                P[j] = false
                blocked[j] = true
                break
            end
            entered = false
            α = minimum(z[i] / (z[i] - t[i]) for i in idx if t[i] <= 0)
            z .+= α .* (t .- z)
            for i in idx
                z[i] <= tol && (P[i] = false; z[i] = 0.0)
            end
        end
        blocked[j] || (blocked .= false)
        w = A' * (b - A * z)
    end
    return z
end

# splitmix64 gives a counter-based stream per element, so samples do not depend on threads.
@inline function _splitmix(s::UInt64)
    s += 0x9e3779b97f4a7c15
    z = s
    z = (z ⊻ (z >> 30)) * 0xbf58476d1ce4e5b9
    z = (z ⊻ (z >> 27)) * 0x94d049bb133111eb
    return s, z ⊻ (z >> 31)
end
@inline function _uniform(s::UInt64)
    s, z = _splitmix(s)
    return s, (z >> 11) * 0x1.0p-53
end
