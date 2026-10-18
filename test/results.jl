# Exported results: the CSV files of write_results.

@testset "write_results" begin
    t = hohlraum(; length = 6.0, radius = 2.0, leh_radius = 1.0, capsule_radius = 0.6,
                 nz = 8, nphi = 16, capsule_ntheta = 6, capsule_nphi = 12)
    q = [remake(b; aim = (0.0, 0.0, sign(b.port[3]) * 2.5)) for b in nif_quads(power = 4.0)]
    dep = deposit(DepositionOperator(t, q), q)
    sol = solve(RadiationModel(t; closed = true), dep)
    dir = mktempdir()
    files = write_results(dir, t, dep; sol, beams = q, observe = :capsule, lmax = 4)
    @test basename.(files) == ["elements.csv", "components.csv", "beams.csv", "legendre.csv",
                               "ylm.csv", "summary.txt"]
    el = readlines(joinpath(dir, "elements.csv"))
    @test length(el) == nelements(t) + 1 && occursin("incident_Wcm2", el[1])
    hdr = split(el[1], ",")
    col(n) = findfirst(==(n), hdr)
    rows = [split(l, ",") for l in el[2:end]]
    @test sum(parse(Float64, r[col("laser_TW")]) for r in rows) ≈ sum(dep.element) rtol = 1e-9
    @test parse(Float64, rows[5][col("incident_Wcm2")]) ≈ to_Wcm2(sol.H[5]) rtol = 1e-9
    comps = readlines(joinpath(dir, "components.csv"))
    @test length(comps) == length(t.components) + 1
    bl = readlines(joinpath(dir, "beams.csv"))
    @test length(bl) == length(q) + 1
    yl = readlines(joinpath(dir, "ylm.csv"))
    @test length(yl) == 1 + 25 && startswith(yl[2], "0,0,")
    @test occursin("laser_delivered_TW = ", read(joinpath(dir, "summary.txt"), String))
    # without a radiosity solution: the laser files only
    @test basename.(write_results(mktempdir(), t, dep)) == ["elements.csv", "components.csv", "summary.txt"]
end
