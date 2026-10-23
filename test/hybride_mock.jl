# Mock HYBRID-E pointing: builds on the 192 beams and the 48 quads, and has the stated
# symmetry (one rule per cone in each quad's frame, mirrored top to bottom, outer quads
# split across their own ports).

@testset "hybride_mock_pointing" begin
    nb = nif_beams()
    bs = hybride_mock_pointing(nb)
    @test length(bs) == 192
    @test all(b -> b.defocus == 0, bs)
    @test [b.power for b in bs] == [b.power for b in nb]
    @test [b.port for b in bs] == [b.port for b in nb]
    tab = Briareus.HYBRIDE_MOCK_POINTING
    quads = unique(b.quad for b in bs)
    @test length(quads) == 48
    φq = Dict(q => (s = sum(Vec3(b.port...) for b in bs if b.quad == q); atan(s[2], s[1])) for q in quads)
    for q in quads
        mem = [b for b in bs if b.quad == q]
        @test length(mem) == 4
        g = mem[1].group
        p = tab[g]
        er, eφ = Vec3(cos(φq[q]), sin(φq[q]), 0), Vec3(-sin(φq[q]), cos(φq[q]), 0)
        # up–down mirror: the same |z| in both hemispheres
        @test all(b -> b.aim[3] ≈ sign(b.port[3]) * p.z, mem)
        # quad frame: radius along the quad azimuth, the split along its φ̂
        @test all(b -> dot(b.aim, er) ≈ p.r, mem)
        t = [dot(b.aim, eφ) for b in mem]
        if p.s == 0
            @test all(x -> abs(x) < 1e-12, t)                   # inner quads share one aim
        else
            @test sort(t) ≈ [-p.s, -p.s, p.s, p.s]               # two beams each side
            for (b, ti) in zip(mem, t)                           # each crosses to the far side
                δ = mod(atan(b.port[2], b.port[1]) - φq[q] + π, 2π) - π
                @test sign(ti) == -sign(δ)
            end
        end
    end
    # the cones' azimuthal pattern: each cone's aims are its quads' rotated copies
    for g in (:inner23, :inner30, :outer44, :outer50)
        mean_aims = [sum(b.aim for b in bs if b.quad == q) / 4 for q in quads
                     if first(b for b in bs if b.quad == q).group == g]
        @test all(a -> hypot(a[1], a[2]) ≈ abs(tab[g].r) || hypot(a[1], a[2]) < 1e-12, mean_aims)
        @test abs(sum(a[3] for a in mean_aims)) < 1e-12           # top and bottom balance
    end
    # quads as single beams: no split
    qs = hybride_mock_pointing(nif_quads())
    @test length(qs) == 48
    @test all(b -> hypot(b.aim[1], b.aim[2]) ≈ abs(tab[b.group].r) ||
                   hypot(b.aim[1], b.aim[2]) < 1e-12, qs)
end
