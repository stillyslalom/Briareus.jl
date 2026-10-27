# Package hygiene (Aqua.jl): method ambiguities and type piracy in Briareus, unbound type
# parameters, undefined exports, stale dependencies, compat bounds for every dependency
# and test extra (Julia's General registry requires them), and tasks left running after
# loading. Aqua is a test-only dependency (Project.toml [extras]).

using Aqua

@testset "Aqua" begin
    Aqua.test_all(Briareus)
end
