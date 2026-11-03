# The optimizer's box (`bounds`) for GroupPointing and FramePointing: a group at zero
# power can come back on; unknown fields are refused.

@testset "optimizer bounds" begin
    q = nif_quads(power = 8.0)
    # inner23 dropped (zero power), the others at 8 TW per beam
    q0 = [b.group === :inner23 ? remake(b; power = 0.0) : b for b in q]
    for p in (GroupPointing(q0; fields = (:z, :power)),
              FramePointing(q0, CHAMBER; fields = (:u, :power)))
        x0 = initial_parameters(p)
        lo, hi = bounds(p, x0; power = 0.5)
        ip = 2:2:length(x0)                                   # the power parameters
        k = findfirst(==(:inner23), p.keys)
        others = setdiff(eachindex(p.keys), k)
        @test x0[ip[k]] == 0 && lo[ip[k]] == 0 && hi[ip[k]] == 0.5 * 8.0
        @test lo[ip[others]] == fill(4.0, length(others)) && hi[ip[others]] == fill(12.0, length(others))
        @test bounds(p, x0; power = 0.9, power_max = 3.0)[2][ip[k]] == 3.0
        @test lo[1:2:end] == x0[1:2:end] .- 1 && hi[1:2:end] == x0[1:2:end] .+ 1
    end
    # every group at zero power: the box needs a finite power_max
    pz = GroupPointing([remake(b; power = 0.0) for b in q]; fields = (:power,))
    @test_throws ArgumentError bounds(pz)
    @test bounds(pz; power_max = 9.0) == (zeros(4), fill(9.0, 4))
    # a field without a width is refused, not given the defocus width
    pd = GroupPointing(q; fields = (:defocus,))
    @test_throws ArgumentError Briareus._box_bounds(pd, initial_parameters(pd), (; z = 1.0), 0.5, 0.5, 0.0, Inf)
end
