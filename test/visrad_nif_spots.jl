# NIF beams on VisRad's phase plates, their spots in the quad frame: the major axis along
# the quad port's horizontal (φ̂) projected across the beam, by vector geometry; the turn
# against the range the VisRad deposition gave (research notes, "VisRad beam-profile
# probe"); and the default beams unchanged.

using Briareus, LinearAlgebra, Test

@testset "NIF beams on VisRad's phase plates" begin
    # the default is bitwise the beam frame
    @test nif_beams() == nif_beams(spot_frame = :beam)
    @test_throws ArgumentError nif_beams(spot_frame = :cone)
    # the plates by cone, and their spots
    @test sort(collect(keys(VISRAD_NIF_PLATES))) == [:inner23, :inner30, :outer44, :outer50]
    @test all(c -> visrad_nif_spot(c) == visrad_spot(VISRAD_NIF_PLATES[c]), keys(VISRAD_NIF_PLATES))
    # VisRad's Rev3-Inner radii (HWHM 890 × 660 µm, n = 7.6) as 1/e semi-axes
    s = visrad_nif_spot(:inner23)
    @test (s.a, s.b, s.order) == (0.890 / log(2)^(1 / 7.6), 0.660 / log(2)^(1 / 7.6), 7.6)

    beams = nif_beams(; spot = visrad_nif_spot, spot_frame = :quad)
    plain = nif_beams(; spot = visrad_nif_spot)
    @test all(i -> beams[i].spot.a == plain[i].spot.a && beams[i].spot.b == plain[i].spot.b, eachindex(beams))
    # the major (a) axis is the quad port's φ̂ projected normal to the beam
    qport = Dict(q.quad => q.port for q in nif_quads())
    function along_quad(b)
        qφ = normalize(Vec3(-qport[b.quad][2], qport[b.quad][1], 0.0))
        k = -b.port
        proj = normalize(qφ - dot(qφ, k) * k)
        return abs(dot(spot_axes(b)[1], proj))
    end
    @test all(b -> along_quad(b) > 1 - 1e-12, beams)
    @test !all(b -> along_quad(b) > 1 - 1e-6, plain)
    # the turn: 1.6–5.3° (largest on the 23.5° cone), as VisRad's spots are turned
    turn = [rad2deg(abs(b.spot.angle - p.spot.angle)) for (b, p) in zip(beams, plain)]
    @test 1.5 < minimum(turn) < 1.7 && 5.2 < maximum(turn) < 5.4
    @test maximum(turn[i] for i in eachindex(beams) if beams[i].group === :inner23) == maximum(turn)
end
