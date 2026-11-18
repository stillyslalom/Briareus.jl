# NIF quads and their beams: a quad model split into its facility beams and back, the
# mean aim of a split quad against the closed form of a pointing rule, and a per-beam
# pointing table applied to a quad model through the SST reader.

using Briareus, LinearAlgebra, Test

@testset "NIF quads and their beams" begin
    byname(bs) = sort(bs; by = b -> b.name)
    # the beams of the quads are the facility beams, bit for bit
    @test byname(nif_beams(nif_quads(power = 8.0))) == byname(nif_beams(power = 2.0))

    # a quad model's aim, defocus, power, detuning and group reach its four beams
    quads = [remake(q; aim = (0.1i, -0.2, sign(q.port[3]) * (3.0 + 0.01i)), defocus = 0.05i,
                    power = 4.0 + 0.1i, dlambda = 0.3 * isodd(i), group = Symbol(:g, mod(i, 3)))
             for (i, q) in enumerate(nif_quads())]
    beams = nif_beams(quads)
    @test length(beams) == 192 && [b.quad for b in beams[1:4:end]] == [q.quad for q in quads]
    fields(b) = (b.aim, b.defocus, b.dlambda, b.group, b.wavelength)
    @test all(i -> fields(beams[i]) == fields(quads[cld(i, 4)]), eachindex(beams))
    @test all(i -> beams[i].power == quads[cld(i, 4)].power / 4, eachindex(beams))
    @test all(b -> b.fnumber == 20.7, beams)
    # and back: the quads at their ports, with their values again
    back = nif_quads(beams)
    @test [q.port for q in back] == [q.port for q in quads]
    @test all(i -> isapprox(back[i].aim, quads[i].aim; rtol = 1e-14), eachindex(quads))
    @test all(i -> isapprox(back[i].power, quads[i].power; rtol = 1e-14), eachindex(quads))
    @test all(i -> (back[i].group, back[i].fnumber) == (quads[i].group, 7.9), eachindex(quads))
    @test all(i -> isapprox([back[i].defocus, back[i].dlambda], [quads[i].defocus, quads[i].dlambda];
                            rtol = 1e-14, atol = 1e-15), eachindex(quads))
    # spots: each beam's quad spot, in the quad frame and back
    vq = nif_quads(spot = visrad_nif_spot)
    turned = nif_beams(vq; spot_frame = :quad)
    @test all(i -> turned[i].spot.a == vq[cld(i, 4)].spot.a, eachindex(turned))
    @test all(i -> isapprox(nif_quads(turned; spot_frame = :quad)[i].spot.angle, vq[i].spot.angle;
                            atol = 1e-14), eachindex(vq))
    @test byname(nif_beams(vq; spot_frame = :quad)) == byname(nif_beams(; spot = visrad_nif_spot, spot_frame = :quad, power = 2.0))
    @test nif_quads(nif_beams(); spot = nif_spot)[1].spot == nif_spot(nif_beams()[1].group)

    # a split quad's mean aim: the mock HYBRID-E rule's outer beams split ±0.5 mm along the
    # quad's φ̂ cancel in pairs, leaving r along the quad's mean-port azimuth and ±z
    # (the rule's values, from its docstring)
    rule = Dict(:inner23 => (3.3, -1.3), :inner30 => (4.2, -0.8), :outer44 => (5.7, 0.0),
                :outer50 => (5.8, 0.0))
    nb = nif_beams()
    split = hybride_mock_pointing(nb)
    meanq = nif_quads(split)
    function expected(q)
        port = sum(b.port for b in nb if b.quad == q.quad)
        z, r = rule[q.group]
        φ = atan(port[2], port[1])
        return [r * cos(φ), r * sin(φ), sign(port[3]) * z]
    end
    @test all(q -> isapprox(q.aim, expected(q); atol = 1e-14), meanq)
    outer = [b for b in split if b.group === :outer50]
    @test maximum(norm(b.aim - only(q for q in meanq if q.quad == b.quad).aim) for b in outer) ≈ 0.5 rtol = 1e-12

    # a per-beam pointing table on a quad model: the SST reader with the quads' beams as
    # template keeps each beam's quarter power and detuning and takes its own aim
    path = tempname() * ".csv"
    write_visrad_sst(path, split)
    shot = read_visrad_sst(path; template = nif_beams(quads))
    sq = Dict(b.name => b for b in split)
    @test all(b -> isapprox(b.aim, sq[b.name].aim; atol = 1e-6), shot)
    @test all(i -> (shot[i].power, shot[i].dlambda) == (beams[i].power, beams[i].dlambda), eachindex(shot))
    rm(path)

    @test_throws ArgumentError nif_beams(omega_beams())
    @test_throws ArgumentError nif_quads(omega_beams())
    @test_throws ArgumentError nif_quads([remake(nb[1]; group = :a), remake(nb[2]; group = :b)])
    @test_throws ArgumentError nif_beams(nif_quads(); spot_frame = :cone)
end
