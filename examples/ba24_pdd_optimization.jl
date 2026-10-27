# Optimizing Barlow et al. 2024's NIF polar-direct-drive pointing for the vacuum
# illumination nonuniformity σ₁ alone (rms over the sphere relative to the mean, l ≤ 30),
# by gradient (L-BFGS-B with ForwardDiff through `deposit`, from their Optimized Config
# and from random starts) and by two derivative-free baselines already in Optim
# (Nelder–Mead and particle swarm), counting evaluations to converge. Decision 6 of
# reference/litreview/polar_direct_drive.md; the comparison of the design itself with
# Ba24's numbers is the tutorial docs/literate/polar_direct_drive.jl.
#
# Ba24 (PRL 133, 175101, 2024) is the only published genetic-algorithm PDD optimizer. Its
# fitness also holds σ₂ and the ablation pressure of a 4-ns plasma with CBET, which need
# Stages 1–3 of the review's plan; here only its σ₁ term is optimized, CBET-free, over
# its 16 parameters (per cone pair: surface polar angle θ_s, azimuth offset Δφ_s from the
# port, defocus 0–10 mm, power 2–4 TW/quad, Ba24's bounds) with Briareus's own beams:
# LMF-measured CPP spots, f-cone defocus, all four beams of a quad aimed at the quad's
# surface point (Ba24's overlap assumption). Research notes, "Hard-sphere direct-drive
# illumination", "Ba24 σ₁ comparison and optimization".
#
#     julia --project=. -t auto examples/ba24_pdd_optimization.jl [part]
#
# part: lbfgs, multistart, nm (Nelder–Mead), pso (particle swarm), or all (the default;
# about 20 min on 8 threads).
using Briareus, LinearAlgebra, Printf, Random
using Optim

const PART = isempty(ARGS) ? "all" : ARGS[1]
const R, L = 1.1, 30

# Ba24 Table I, the Optimized Config (OC), per cone: polar angle from the pole, azimuth
# offset from the port (the reading the tutorial settles), defocus (mm), power (% of
# 4 TW/quad)
const OC = Dict(:inner23 => (6.36, 52.02, 8.22, 66.26), :inner30 => (14.95, 160.75, 9.07, 66.26),
                :outer44 => (52.48, -12.14, 4.67, 95.05), :outer50 => (85.62, 6.90, 6.83, 82.88))

# The 16 parameters act on the 48 quads (`GroupPointing`'s spherical fields; the quad's
# port is the mean of its four beams'); every beam of a quad takes its quad's aim,
# defocus and a quarter of its power.
struct QuadAims
    quads::GroupPointing{Symbol}
    beams::Vector{Beam{Float64}}
    quad_of::Vector{Int}
end
function QuadAims()
    q = nif_quads()
    beams = nif_beams()
    QuadAims(GroupPointing(q; fields = (:theta_s, :dphi_s, :defocus, :power), radius = R), beams,
             [findfirst(x -> x.quad == b.quad, q) for b in beams])
end
function Briareus.apply(p::QuadAims, x::AbstractVector)
    qb = apply(p.quads, x)
    return [remake(b; aim = qb[j].aim, defocus = qb[j].defocus, power = qb[j].power / 4)
            for (b, j) in zip(p.beams, p.quad_of)]
end

pq = QuadAims()
x_oc = reduce(vcat, [[OC[k][1], OC[k][2], OC[k][3], 4 * OC[k][4] / 100] for k in pq.quads.keys])
# the scenario's deposition operator is built for the beams at the initial parameters
Briareus.initial_parameters(::QuadAims) = x_oc
lower = repeat([0.0, -180.0, 0.0, 2.0], 4)
upper = repeat([90.0, 180.0, 10.0, 4.0], 4)
sphere(n) = Target(Component(:capsule, icosphere_mesh(R; n), Material(:CH); surface = SphereSurface(R)))

# σ₁ on an n = 20 icosphere (8000 elements; σ 0.6 % low against n = 30 at the OC and
# 1.2 % against the extrapolation, a smooth bias the optimum barely feels), rechecked at
# n = 30. The bounds are Ba24's but do not enforce their "visible half of the target"
# rule: an aim behind the limb is still a well-defined axis, which first meets the
# sphere elsewhere.
function setup(n)
    t = sphere(n)
    W = ylm_weights(t, collect(1:nelements(t)); lmax = L)
    s = Scenario(t, pq)
    σ1(ev) = mode_rms(W * (ev.deposition.element ./ t.elem_area), L)
    return (; t, W, s, σ1)
end
t0 = time()
S20 = setup(20)
S30 = setup(30)
@printf("setup (operators and Y_lm weights, n = 20 and 30): %.1f s\n", time() - t0)
σ30(x) = S30.σ1(evaluate(S30.s, x))
σ20(x) = S20.σ1(evaluate(S20.s, x))
@printf("OC: σ₁ = %.3f %% (n = 20), %.3f %% (n = 30); Ba24 6.03 %%\n", 100σ20(x_oc), 100σ30(x_oc))

show_x(x) = join([@sprintf("%s θ %.2f Δφ %.2f d %.2f P %.2f", k, x[4i-3:4i]...)
                  for (i, k) in enumerate(pq.quads.keys)], " | ")

results = []
function run_lbfgs!(results)
    σ20(x_oc); t0 = time()
    marks = (10, 25, 50, 100, 150, 200)
    progress(st) = (st.evaluations in marks &&
                    @printf("  L-BFGS-B: %d evaluations, best σ₁ %.3f %% (%.0f s)\n", st.evaluations,
                            100st.best_value, time() - t0); false)
    r = optimize_pointing(S20.s, S20.σ1; x0 = x_oc, lower, upper, iterations = 200, callback = progress)
    wall = time() - t0
    @printf("L-BFGS-B from the OC: σ₁ %.3f %% -> %.3f %% (n = 20), %.3f %% at n = 30; %d evaluations, %d iterations, %s, %.0f s\n",
            100r.initial_value, 100r.value, 100σ30(r.x), r.evaluations, r.iterations, r.reason, wall)
    println("  ", show_x(r.x))
    push!(results, ("L-BFGS-B from OC", r.value, σ30(r.x), r.evaluations, wall))
end

function run_multistart!(results; nstart = 4, budget = 75)
    rng = Xoshiro(2024)
    best = (value = Inf, x = x_oc)
    total, t0 = 0, time()
    for k in 1:nstart
        x0 = lower .+ rand(rng, length(lower)) .* (upper .- lower)
        r = optimize_pointing(S20.s, S20.σ1; x0, lower, upper, iterations = 150, evaluations = budget)
        total += r.evaluations
        @printf("  start %d: σ₁ %.2f %% -> %.3f %% in %d evaluations (%s)\n", k, 100r.initial_value,
                100r.value, r.evaluations, r.reason)
        r.value < best.value && (best = (value = r.value, x = r.x))
    end
    wall = time() - t0
    @printf("multistart (%d random starts): best σ₁ %.3f %% (n = 20), %.3f %% at n = 30; %d evaluations in all, %.0f s\n",
            nstart, 100best.value, 100σ30(best.x), total, wall)
    println("  ", show_x(best.x))
    push!(results, ("L-BFGS-B, $nstart random starts", best.value, σ30(best.x), total, wall))
end

# Derivative-free baselines on the same objective (Float64 evaluations), in parameters
# scaled to the unit box, capped at `budget` evaluations; the best value so far is
# printed at checkpoints for the evaluations-to-converge comparison.
function run_dfo!(results, name; budget = 1200)
    scale(x) = (x .- lower) ./ (upper .- lower)
    unscale(u) = lower .+ clamp.(u, 0, 1) .* (upper .- lower)
    method = name == "Nelder-Mead" ? NelderMead() :
             ParticleSwarm(; lower = zeros(16), upper = ones(16), n_particles = 24)
    count, bestv, bestu = Ref(0), Ref(Inf), Ref(scale(x_oc))
    marks = (50, 100, 200, 400, 700, 1000, 1500, 3000)
    function f(u)
        count[] += 1
        v = σ20(unscale(u))
        v < bestv[] && (bestv[] = v; bestu[] = copy(u))
        count[] in marks && @printf("  %s: %d evaluations, best σ₁ %.3f %%\n", name, count[], 100bestv[])
        return v
    end
    t0 = time()
    optimize(f, scale(x_oc), method, Optim.Options(; f_calls_limit = budget, iterations = budget))
    wall = time() - t0
    xb = unscale(bestu[])
    @printf("%s from the OC: best σ₁ %.3f %% (n = 20), %.3f %% at n = 30 after %d evaluations, %.0f s\n",
            name, 100bestv[], 100σ30(xb), count[], wall)
    println("  ", show_x(xb))
    push!(results, ("$name from OC", bestv[], σ30(xb), count[], wall))
end

PART in ("all", "lbfgs") && run_lbfgs!(results)
PART in ("all", "multistart") && run_multistart!(results)
PART in ("all", "nm") && run_dfo!(results, "Nelder-Mead")
PART in ("all", "pso") && run_dfo!(results, "particle swarm")

println("\n| run | σ₁ (n = 20) | σ₁ (n = 30) | evaluations | wall (s) |")
println("|---|---|---|---|---|")
for (name, v20, v30, ne, w) in results
    @printf("| %s | %.3f %% | %.3f %% | %d | %.0f |\n", name, 100v20, 100v30, ne, w)
end
