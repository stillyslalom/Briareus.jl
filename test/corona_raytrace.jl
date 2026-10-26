# Independent reference for `CoronaOperator` (test helper, not package API): 3-D ray
# tracing through a corona n_e(x) by RK4, the prototype's model (c) (research notes,
# "Corona model prototype"). It shares only the beam model with Briareus (f-cone
# directions and the best-focus profile); paths, absorption, angular binning and the
# harmonic projection are its own. With the ray parameter τ (ds = n dτ):
#     dx/dτ = v,   dv/dτ = -∇(n_e/n_c)/2,   d(optical depth)/dτ = (ν_ei/c)(n_e/n_c),
# so for an isothermal profile the rate is (ν_c/c)(n_e/n_c)² (Kruer, ch. 5). Rays start
# on a square lattice (spacing `hlat`) in the plane normal to the beam through the centre,
# one lattice per f-cone direction, weighted by the best-focus profile where they cross
# the focal plane × hlat²/ncone. Each RK4 step's absorbed power is put at the step's
# midpoint, cloud-in-cell on an equal-angle (θ, φ) grid, then projected on real Y_lm
# (normalized associated Legendre functions by their own recursion; DFT in φ).

# Isothermal exponential corona n_e/n_c = exp(-(r - rc)/L), zero beyond R.
_rt_ne(rc, L, R, x) = (r = norm(x); r > R ? 0.0 : exp(-(r - rc) / L))
_rt_grad(rc, L, R, x) = (r = norm(x); r > R ? zero(x) : -exp(-(r - rc) / L) / (L * r) * x)

# One ray from x (on the sphere of radius R) along the unit vector v: RK4 steps hτ until
# it leaves the sphere outward; deposits w × (absorbed fraction) on the grid E.
function _rt_ray!(E, nθ, x, v, hτ, w, (rc, L, νc, R))
    grad(y) = -0.5 * _rt_grad(rc, L, R, y)
    rate(y) = νc * _rt_ne(rc, L, R, y)^2
    P, ℓ = 1.0, 0.0
    for _ in 1:10^6
        k1x, k1v, k1l = v, grad(x), rate(x)
        x2 = x + hτ / 2 * k1x
        k2x, k2v, k2l = v + hτ / 2 * k1v, grad(x2), rate(x2)
        x3 = x + hτ / 2 * k2x
        k3x, k3v, k3l = v + hτ / 2 * k2v, grad(x3), rate(x3)
        x4 = x + hτ * k3x
        k4x, k4v, k4l = v + hτ * k3v, grad(x4), rate(x4)
        xn = x + hτ / 6 * (k1x + 2k2x + 2k3x + k4x)
        v += hτ / 6 * (k1v + 2k2v + 2k3v + k4v)
        ℓ += hτ / 6 * (k1l + 2k2l + 2k3l + k4l)
        Pn = exp(-ℓ)
        P > Pn && _rt_cic!(E, nθ, (x + xn) / 2, w * (P - Pn))
        x, P = xn, Pn
        (dot(x, x) > R^2 && dot(x, v) > 0) && return w * (1 - P)
    end
    error("ray did not leave the corona")
end

# Cloud-in-cell deposit of p at the direction of x on the nθ × 2nθ grid (rows θ).
function _rt_cic!(E, nθ, x, p)
    nφ = 2nθ
    sθ = acos(clamp(x[3] / norm(x), -1.0, 1.0)) / π * nθ - 0.5
    sφ = mod(atan(x[2], x[1]), 2π) / 2π * nφ - 0.5
    i0, j0 = floor(Int, sθ), floor(Int, sφ)
    tθ, tφ = sθ - i0, sφ - j0
    i1 = min(i0 + 1, nθ - 1)
    i0 = max(i0, 0)
    j1, j0 = mod(j0 + 1, nφ), mod(j0, nφ)
    E[i0+1, j0+1] += p * (1 - tθ) * (1 - tφ)
    E[i0+1, j1+1] += p * (1 - tθ) * tφ
    E[i1+1, j0+1] += p * tθ * (1 - tφ)
    E[i1+1, j1+1] += p * tθ * tφ
    return E
end

# Fully normalized P̄_lm(μ) (Y_l0 = P̄_l0, Y_lm = √2 P̄_lm cos mφ / sin |m|φ), by the
# standard recursion in l at fixed m.
function _rt_plm(L, μ)
    P = zeros(L + 1, L + 1)
    s = sqrt(1 - μ^2)
    P[1, 1] = 1 / sqrt(4π)
    for m in 1:L
        P[m+1, m+1] = sqrt((2m + 1) / 2m) * s * P[m, m]
    end
    for m in 0:L, l in m+1:L
        a = sqrt((4l^2 - 1) / (l^2 - m^2))
        b = l >= m + 2 ? sqrt(((l - 1)^2 - m^2) / (4(l - 1)^2 - 1)) : 0.0
        P[l+1, m+1] = a * (μ * P[l, m+1] - (l >= m + 2 ? b * P[l-1, m+1] : 0.0))
    end
    return P
end

# Flat real-harmonic coefficients (mode_index order) of the binned power E.
function _rt_project(E, L)
    nθ, nφ = size(E)
    c = zeros((L + 1)^2)
    cs = [cos(m * (j - 0.5) * 2π / nφ) for j in 1:nφ, m in 0:L]
    sn = [sin(m * (j - 0.5) * 2π / nφ) for j in 1:nφ, m in 0:L]
    for i in 1:nθ
        P = _rt_plm(L, cos((i - 0.5) * π / nθ))
        for m in 0:L
            fc = fs = 0.0
            for j in 1:nφ
                fc += E[i, j] * cs[j, m+1]
                fs += E[i, j] * sn[j, m+1]
            end
            for l in m:L
                m == 0 ? (c[l^2+l+1] += P[l+1, 1] * fc) :
                (c[l^2+l+m+1] += sqrt(2) * P[l+1, m+1] * fc; c[l^2+l-m+1] += sqrt(2) * P[l+1, m+1] * fs)
            end
        end
    end
    return c
end

# Absorbed power per steradian of `beams` in the corona (rc, L, νc, R) about the origin:
# flat Y_lm coefficients to degree L and the absorbed power.
function _rt_corona(beams, prof; ncone, hlat, hτ, nθ = 360, lmax = 30)
    R = prof[4]
    work = Tuple{Vec3,Vec3,Float64}[]
    for b in beams
        f = Briareus.BeamFrame(b)
        ds, _ = Briareus.cone_directions(b, ncone)
        n = ceil(Int, R / hlat) + 1
        for d in ds, iu in -n:n, iv in -n:n
            x0 = (iu * hlat) * f.e1 + (iv * hlat) * f.e2
            # entry on the sphere |x| = R of the line x0 + t d
            q = dot(x0, d)
            disc = q^2 - dot(x0, x0) + R^2
            disc > 0 || continue
            xe = x0 + (-q - sqrt(disc)) * d
            # focal-plane crossing of the ray (a sub-beam carries the best-focus profile)
            r = x0 - f.focus
            s = dot(r, f.k) / dot(d, f.k)
            p = r - s * d
            I = Briareus._sg_focal(f.core, dot(p, f.e1), dot(p, f.e2))
            f.has_halo && (I += Briareus._sg_focal(f.halo, dot(p, f.e1), dot(p, f.e2)))
            I > 0 && push!(work, (xe, d, I * hlat^2 / ncone))
        end
    end
    nt = Threads.nthreads()
    grids = [zeros(nθ, 2nθ) for _ in 1:nt]
    absd = zeros(nt)
    Threads.@threads :static for t in 1:nt
        for k in t:nt:length(work)
            xe, d, w = work[k]
            absd[t] += _rt_ray!(grids[t], nθ, xe, d, hτ, w, prof)
        end
    end
    return _rt_project(sum(grids), lmax), sum(absd), length(work)
end
