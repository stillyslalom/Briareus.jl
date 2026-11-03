# Per-target caches (a Dante view's footprint, an imager's pixels) hold their targets
# weakly: an entry is released once its target is unreachable.

using Briareus, Test

# A keyframe of a moving target: a target built, both caches filled, and only a weak
# reference to it returned (the targets are made in functions of their own, so that no
# frame of the test keeps one alive).
@noinline function _cache_keyframe(v, im)
    t = hohlraum(nz = 4, nphi = 8, capsule_ntheta = 4, capsule_nphi = 6)
    Briareus._footprint(v, t)
    Briareus._imager_pixels(im, t)
    return WeakRef(t)
end
@noinline _cache_keyframes(v, im, n) = [_cache_keyframe(v, im) for _ in 1:n]

@testset "per-target caches release unreachable targets" begin
    v = dante1()
    im = XrayImager(0, 0; pinhole = 0.025, distance = 1656.7, magnification = 2.0, width = 6.0,
                    pixels = 32, energy = 4000)
    ws = _cache_keyframes(v, im, 3)
    GC.gc(true); GC.gc(true)
    @test all(w -> w.value === nothing, ws)
    @test isempty(v.cache) && isempty(im.cache)
    # a live target keeps its entry, computed once
    t = hohlraum(nz = 4, nphi = 8, capsule_ntheta = 4, capsule_nphi = 6)
    fp = Briareus._footprint(v, t)
    GC.gc(true)
    @test Briareus._footprint(v, t) === fp && length(v.cache) == 1
    # two equal targets built apart are two entries
    t2 = hohlraum(nz = 4, nphi = 8, capsule_ntheta = 4, capsule_nphi = 6)
    Briareus._footprint(v, t2)
    @test length(v.cache) == 2 && t !== t2
end
