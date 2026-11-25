# The boundary with Briareus: the source of the application calls the exported API and
# the names of `Briareus.GUIAPI`, and no other private `Briareus._name`. The replayed
# precompile list is generated and is excluded.

using Test, Briareus

@testset "no private Briareus names outside GUIAPI" begin
    src = normpath(joinpath(@__DIR__, "..", "src"))
    files = [joinpath(root, f) for (root, _, fs) in walkdir(src) for f in fs
             if endswith(f, ".jl") && f != "first_frame.jl"]
    hits = [string(relpath(f, src), ": ", m.match) for f in files
            for m in eachmatch(r"Briareus\._[A-Za-z_0-9!]+", read(f, String))]
    @test isempty(hits)
    # every GUIAPI name that the source uses is bound
    used = unique(m.captures[1] for f in files for m in eachmatch(r"Briareus\.GUIAPI\.([A-Za-z_0-9!]+)", read(f, String)))
    @test !isempty(used) && all(n -> isdefined(Briareus.GUIAPI, Symbol(n)), used)
end
