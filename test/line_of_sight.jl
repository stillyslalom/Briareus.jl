# Line-of-sight power from partly hidden elements, against the visible area in closed
# form: a plate half hidden behind a straight edge (an occluding sheet above it, its edge
# oblique to the mesh), seen straight on and at 30°.

@testset "line of sight: a plate half hidden behind an edge" begin
    plate = Component(:plate, rectangle_mesh(2.0, 2.0; nu = 3, nv = 3), Material(:p; albedo = 0.0, conversion = 0.0))
    β, edge, h = 20.0, 0.2, 0.5          # the edge's normal (degrees from x), its offset, the sheet's height
    R = [cosd(β) -sind(β) 0; sind(β) cosd(β) 0; 0 0 1]
    sheet = Component(:sheet, transform(rectangle_mesh(4.0, 4.0; nu = 1, nv = 1), R, R * Vec3(edge - 2, 0, h)), BLACK)
    t = Target(plate, sheet)
    B = ones(nelements(t))
    on = elements(t, :plate)
    n̂ = Vec3(cosd(β), sind(β), 0)
    # the area of the square |x|, |y| ≤ 1 where n̂ ⋅ p ≥ c (the square clipped by the line)
    function clipped(c)
        sq = [Vec3(-1, -1, 0), Vec3(1, -1, 0), Vec3(1, 1, 0), Vec3(-1, 1, 0)]
        out = Vec3[]
        for i in 1:4
            p, q = sq[i], sq[mod1(i + 1, 4)]
            fp, fq = dot(n̂, p) - c, dot(n̂, q) - c
            fp >= 0 && push!(out, p)
            fp * fq < 0 && push!(out, p + (q - p) * (fp / (fp - fq)))
        end
        return sum(out[i][1] * out[mod1(i + 1, end)][2] - out[mod1(i + 1, end)][1] * out[i][2] for i in eachindex(out)) / 2
    end
    for θ in (0.0, 30.0)
        d = Vec3(sind(θ), 0, cosd(θ))
        # a point p is hidden when p + (h / cos θ) d lies under the sheet
        exact = clipped(edge - h * tand(θ) * n̂[1]) * cosd(θ) / π
        P = line_of_sight_power(t, B, d; elems = on)
        P1 = line_of_sight_power(t, B, d; elems = on, nsub = 1, refine = 1)   # one point per triangle
        Pfine = line_of_sight_power(t, B, d; elems = on, nsub = 128, refine = 1)
        @test abs(P / exact - 1) < 3e-4
        @test abs(Pfine / exact - 1) < 3e-5
        @test abs(P1 / exact - 1) > 10 * abs(P / exact - 1)
    end
    # nothing hidden: every element whole, as before
    @test line_of_sight_power(t, B, (0, 0, 1); elems = elements(t, :sheet)) ≈ 16 / π rtol = 1e-12
end
