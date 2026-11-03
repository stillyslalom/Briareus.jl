# Ray casting through the bounding-volume hierarchies (median split and surface-area
# heuristic) against brute force over every triangle, with an intersection written out
# here (plane, then barycentric coordinates) rather than the package's Möller–Trumbore.

# The ray parameter where o + t d meets the triangle (v0, e1, e2), or Inf.
function _brute_t(o, d, v0, e1, e2)
    n = cross(e1, e2)
    dn = dot(n, d)
    dn == 0 && return Inf
    t = dot(n, v0 - o) / dn
    t > 0 || return Inf
    p = o + t * d - v0
    # barycentric coordinates of p in the basis (e1, e2)
    g = [dot(e1, e1) dot(e1, e2); dot(e1, e2) dot(e2, e2)]
    u, v = g \ [dot(p, e1), dot(p, e2)]
    return u >= 0 && v >= 0 && u + v <= 1 ? t : Inf
end

function _brute_closest(o, d, v0, e1, e2, ids; skip = nothing)
    best, id = Inf, Int32(0)
    for k in eachindex(v0)
        skip !== nothing && skip(ids[k]) && continue
        t = _brute_t(o, d, v0[k], e1[k], e2[k])
        t < best && (best = t; id = ids[k])
    end
    return best, id
end

@testset "BVH ray casting against brute force" begin
    # a soup of 300 random triangles in a 4 mm box, and 400 random rays through it
    rng = Xoshiro(11)
    n = 300
    v0 = [Vec3(4 .* rand(rng, 3) .- 2...) for _ in 1:n]
    e1 = [Vec3(0.6 .* randn(rng, 3)...) for _ in 1:n]
    e2 = [Vec3(0.6 .* randn(rng, 3)...) for _ in 1:n]
    ids = Int32.(1000 .+ (1:n))
    median = Briareus.BVH(v0, e1, e2, ids)
    sah = Briareus._sah_bvh(median)
    # from outside the box toward points inside it
    origins = [Vec3(3 .* randn(rng, 3)...) for _ in 1:400]
    rays = [(o, normalize(Vec3(3 .* rand(rng, 3) .- 1.5...) - o)) for o in origins]
    ref = [_brute_closest(o, d, v0, e1, e2, ids) for (o, d) in rays]
    @test count(r -> r[2] != 0, ref) > 200                 # most rays hit something
    same((t, id), (tr, idr)) = id == idr && (id == 0 ? t == Inf : isapprox(t, tr; rtol = 1e-9))
    for b in (median, sah)
        hits = [closest_hit(b, o, d) for (o, d) in rays]
        @test all(same.(hits, ref))
        @test [any_hit(b, o, d) for (o, d) in rays] == [r[2] != 0 for r in ref]
        # tmax just short of and just past the nearest hit
        hit = findall(r -> r[2] != 0, ref)
        @test !any(any_hit(b, rays[i]..., ref[i][1] * (1 - 1e-6)) for i in hit)
        @test all(any_hit(b, rays[i]..., ref[i][1] * (1 + 1e-6)) for i in hit)
        @test all(closest_hit(b, rays[i]..., ref[i][1] * (1 - 1e-6))[2] == 0 for i in hit)
        # skipping the nearest triangle: any_hit then answers for the others
        @test all(any_hit(b, rays[i]..., Inf, id -> id == ref[i][2]) ==
                  (_brute_closest(rays[i]..., v0, e1, e2, ids; skip = id -> id == ref[i][2])[2] != 0) for i in hit)
        # t is in units of |dir|
        i = hit[1]
        @test closest_hit(b, rays[i][1], 2 * rays[i][2])[1] ≈ ref[i][1] / 2 rtol = 1e-12
    end
    # an empty hierarchy hits nothing
    e = Briareus.BVH(Vec3[], Vec3[], Vec3[], Int32[])
    @test closest_hit(e, Vec3(0, 0, 0), Vec3(0, 0, 1)) == (Inf, 0) && !any_hit(e, Vec3(0, 0, 0), Vec3(0, 0, 1))
    @test isempty(Briareus._sah_bvh(e).nodes)
    # an axis-aligned ray (zero direction components) through a known triangle
    tri = Briareus.BVH([Vec3(-1, -1, 0)], [Vec3(2, 0, 0)], [Vec3(0, 2, 0)], Int32[7])
    @test closest_hit(tri, Vec3(-0.5, -0.5, 3), Vec3(0, 0, -1)) == (3.0, 7)
    @test closest_hit(tri, Vec3(0.5, 0.5, 3), Vec3(0, 0, -1)) == (Inf, 0)   # outside the hypotenuse
end

@testset "ray casting on a target" begin
    t = hohlraum(nz = 8, nphi = 12, capsule_ntheta = 6, capsule_nphi = 8)
    rng = MersenneTwister(1)
    ops = findall(t.tri_opaque)
    rays = [(Vec3(randn(rng, 3)...) .* 2, normalize(Vec3(randn(rng, 3)...))) for _ in 1:300]
    # the target's hierarchy holds the opaque triangles, identified by triangle index
    ref = [_brute_closest(o, d, t.tri_v0[ops], t.tri_e1[ops], t.tri_e2[ops], Int32.(ops)) for (o, d) in rays]
    hits = [closest_hit(t.bvh, o, d) for (o, d) in rays]
    @test all(h[1] == r[1] == Inf || isapprox(h[1], r[1]; rtol = 1e-9) for (h, r) in zip(hits, ref))
    @test [any_hit(t.bvh, o, d) for (o, d) in rays] == [r[2] != 0 for r in ref]
    # packets (the deposition operator's visibility rays): the same answers as any_hit,
    # for coherent bundles off a facet and for scattered origins and directions
    agree = true
    cache = zeros(Int32, 16)
    for trial in 1:200
        k = rand(rng, ops)
        side = rand(rng, (-1, 1)) * t.tri_normal[k]
        k0 = normalize(Vec3(randn(rng, 3)...))
        dot(k0, side) > 0 || (k0 = -k0)
        ds = trial % 4 == 0 ? [normalize(Vec3(randn(rng, 3)...)) for _ in 1:16] :
             [normalize(k0 + 0.06 * Vec3(randn(rng, 3)...)) for _ in 1:13]
        os = trial % 4 == 0 ? [Vec3(randn(rng, 3)...) .* 2 for _ in 1:7] :
             [t.tri_v0[k] + u * t.tri_e1[k] + v * t.tri_e2[k] + 1e-6 * side
              for (u, v) in ((rand(rng), rand(rng)) for _ in 1:16) if u + v <= 1]
        m0 = [rand(rng, UInt16) & UInt16((1 << length(ds)) - 1) for _ in os]
        m = Briareus._unblocked!(copy(m0), t.bvh, os, ds, Briareus._safe_inv.(ds), cache)
        for i in eachindex(os), j in eachindex(ds)
            (m0[i] >> (j - 1)) & 1 == 1 || continue
            agree &= ((m[i] >> (j - 1)) & 1 == 1) == !any_hit(t.bvh, os[i], ds[j])
        end
        agree &= all(m .& .~m0 .== 0)
    end
    @test agree
end
