using Briareus
using Briareus: closest_hit, any_hit, _point_polygon_vf
using LinearAlgebra
using ForwardDiff
using Random
using Test

const R_SPHERE = 2.0

# BRIAREUS_TEST_SHARD = 1, 2 or 3 runs a third of the suite (CI splits the one-thread job
# three ways); unset, everything runs. Shards 2 and 3 are the files listed here, shard 1
# the rest (basics.jl among them). Files whose helpers others use share a shard:
# meshops.jl, cadfit.jl, partial_visibility.jl, hmatrix.jl; vrw_fixture.jl, visrad.jl.
const SHARD = get(ENV, "BRIAREUS_TEST_SHARD", "")
SHARD in ("", "1", "2", "3") || error("BRIAREUS_TEST_SHARD must be 1, 2 or 3, got \"$SHARD\"")
# BRIAREUS_LONG_TESTS = true also runs the statistical and convergence checks (Monte Carlo
# against linear propagation, a uniformly fine history; ~3 min at four threads). CI's
# one-thread shards run them; local full runs skip them unless asked.
const LONG_TESTS = get(ENV, "BRIAREUS_LONG_TESTS", "false") == "true"
const SHARDS = Dict(f => 2 for f in ("meshops.jl", "stl.jl", "cadfit.jl", "partial_visibility.jl", "hmatrix.jl",
                                     "unconverted_history.jl",
                                     "inference.jl"))
merge!(SHARDS, Dict(f => 3 for f in ("windows.jl", "halfraum.jl", "burnthrough.jl", "uncertainty.jl",
                                     "implosion.jl", "energetics.jl", "spectral.jl", "large_models.jl",
                                     "direct_drive.jl", "corona.jl", "frames.jl", "aqua.jl")))
in_shard(file) = SHARD == "" || get(SHARDS, file, 1) == parse(Int, SHARD)
shard_include(file) = in_shard(file) && include(file)

# Testsets go in included files, not in this block: its body is one top-level expression,
# which keeps the locals of every testset written here alive until the suite ends.
@testset "Briareus" begin

# Docs cross-references, docstrings and @autodocs coverage (cheap; no Documenter).
include("docrefs.jl")

shard_include("basics.jl")

shard_include("talis.jl")
shard_include("vrw_fixture.jl")
shard_include("visrad.jl")
shard_include("solve_cache.jl")
shard_include("windows.jl")
shard_include("halfraum.jl")
shard_include("burnthrough.jl")
shard_include("meshops.jl")
shard_include("stl.jl")
shard_include("cadfit.jl")
shard_include("uncertainty.jl")
shard_include("results.jl")
shard_include("retune.jl")
shard_include("enclosure.jl")
shard_include("leh_power_balance.jl")
shard_include("capsule_mesh.jl")
shard_include("partial_visibility.jl")
shard_include("montecarlo.jl")
shard_include("corners.jl")
shard_include("dynamics.jl")
shard_include("leh_closure.jl")
shard_include("bubble.jl")
shard_include("implosion.jl")
shard_include("hybride_mock.jl")
shard_include("parallel_gradient.jl")
shard_include("cbet_zones.jl")
shard_include("unconverted.jl")
shard_include("unconverted_history.jl")
shard_include("unconverted_source.jl")
shard_include("sombrero.jl")
shard_include("backlighter.jl")
shard_include("backlighter_signal.jl")
shard_include("frames.jl")
shard_include("energetics.jl")
shard_include("capsule_drive.jl")
shard_include("spectral.jl")
shard_include("dante_response.jl")
shard_include("diagnostics.jl")
shard_include("imager.jl")
shard_include("radiography.jl")
shard_include("large_models.jl")
shard_include("direct_drive.jl")
shard_include("corona.jl")
shard_include("hmatrix.jl")
shard_include("parameters.jl")
shard_include("inference.jl")
shard_include("aqua.jl")
end
