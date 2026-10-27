# Does view-factor balancing (`closed = true`, the `hohlraum()` default) bias the
# capsule drive on the default curved path? A continuum reference by Monte Carlo
# (`Briareus.montecarlo_incident`, `reemit = :hit`) on a finely faceted flat copy of the
# hohlraum (wall h = 0.1 mm, capsule_h = 0.05 mm), against the radiosity model on the
# standard meshes (curved and flat elements, raw and balanced view factors). Source:
# 48 NIF quads at 8 TW at the example cone aims (examples/mesh_convergence.jl's
# nominal design), deposited on each mesh by its own `DepositionOperator` and scaled to
# the reference mesh's total; wall η = 0.8, α = 0.925 and 0.85. Departure: one outer
# quad off (by linearity, the Monte Carlo traces that quad's source alone). Results are
# cached in the temp directory.
#
#     julia --project=. -t auto examples/mc_closure.jl [paths per capsule element at α = 0.925]

using Briareus, Printf, Serialization

const PATHS = isempty(ARGS) ? 45_000 : parse(Int, ARGS[1])
const CACHE = joinpath(tempdir(), "briareus_mc_closure")
mkpath(CACHE)
const GEOM = (; length = 10.1, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9)
const REF = (0.1, 0.05)

function pointed(rz)
    bs = nif_quads(power = 8.0)
    p = GroupPointing(bs; by = :group, fields = (:r, :z))
    x = initial_parameters(p)
    for (g, k) in enumerate(p.keys)
        x[2g-1], x[2g] = rz[k]
    end
    return apply(p, x)
end

const BEAMS = pointed(Dict(:inner23 => (-1.09, 3.09), :inner30 => (-0.6, 3.63),
                           :outer44 => (0.0, 5.2), :outer50 => (0.0, 5.01)))
const DROP = findfirst(b -> b.group === :outer44, BEAMS)
withpower(b, p) = typeof(b)(b.name, b.port, b.aim, p, b.defocus, b.spot, b.fnumber, b.group, b.quad,
                            b.wavelength, b.dlambda)

flat(t::Target) = Target([Component(c.name, c.mesh, c.material; kind = c.kind, two_sided = c.two_sided)
                          for c in t.components]; closed = t.closed)
function target(h, ch, α; curved = true)
    t = hohlraum(; GEOM..., h, capsule_h = ch, wall = gold(albedo = α, conversion = 0.8))
    return curved ? t : flat(t)
end

cached(f, key) = (path = joinpath(CACHE, "$key.jls");
                  isfile(path) ? deserialize(path) : (x = f(); serialize(path, x); x))

# element laser power of the design and of the dropped quad alone
sources(h, ch; curved = true) = cached("src_$(h)_$(ch)_$(curved)") do
    op = DepositionOperator(target(h, ch, 0.9; curved), BEAMS)
    (deposit(op, BEAMS).element,
     deposit(op, [i == DROP ? b : withpower(b, 0.0) for (i, b) in enumerate(BEAMS)]).element)
end

ratios(a) = (a[1], 100a[3] / a[1], 100a[5] / a[1])                 # a₀, a₂/a₀, a₄/a₀ (%)
dep(A, B) = (100 * (A[1] / B[1] - 1), 100A[3] / A[1] - 100B[3] / B[1], 100A[5] / A[1] - 100B[5] / B[1])
se(x) = sqrt(sum(abs2, x .- sum(x) / length(x)) / (length(x) - 1) / length(x))

# Legendre moments of the capsule flux, one column per batch
function montecarlo(h, ch, α, which, paths; batches = 64, seed = 1)
    cached("mc_$(h)_$(ch)_$(α)_$(which)_$(paths)_$(seed)") do
        P, D = sources(h, ch; curved = false)
        t = target(h, ch, α; curved = false)
        cap = elements(t, :capsule)
        _, batch = Briareus.montecarlo_incident(t, which === :P ? P : D, cap; paths, batches, seed)
        legendre_weights(t, cap; lmax = 4) * batch'
    end
end

function model(h, ch, curved)
    cached("model_$(h)_$(ch)_$(curved)") do
        P, D = sources(h, ch; curved)
        k = sum(sources(REF...; curved = false)[1]) / sum(P)
        t0 = target(h, ch, 0.925; curved)
        cap = elements(t0, :capsule)
        W = legendre_weights(t0, cap; lmax = 4)
        F = view_factors(t0)
        rs = vec(sum(F; dims = 2))
        out = Dict{Any,Any}(:N => nelements(t0), :rowerr => (sum(abs, rs .- 1) / length(rs), maximum(abs, rs .- 1)))
        for α in (0.925, 0.85), closed in (false, true)
            m = RadiationModel(target(h, ch, α; curved); F, closed)
            out[(α, closed)] = (W * solve(m, k .* P).H[cap], W * solve(m, k .* (P .- D)).H[cap])
        end
        out
    end
end

println("Reference convergence (α = 0.925, ~32M paths): a₀, a₂/a₀, a₄/a₀")
for (h, ch, n) in ((0.2, 0.05, 4800), (0.1, 0.05, 4800), (0.06, 0.05, 4800), (0.1, 0.035, 2400))
    M = montecarlo(h, ch, 0.925, :P, max(1, PATHS * n ÷ 45_000))
    r = [ratios(c) for c in eachcol(M)]
    @printf("  h %.2f, capsule_h %.3f: a₀ %.5f ± %.5f  a₂/a₀ %+.3f ± %.3f  a₄/a₀ %+.3f ± %.3f\n", h, ch,
            ratios(sum(M; dims = 2) ./ size(M, 2))[1], se(getindex.(r, 1)),
            ratios(sum(M; dims = 2) ./ size(M, 2))[2], se(getindex.(r, 2)),
            ratios(sum(M; dims = 2) ./ size(M, 2))[3], se(getindex.(r, 3)))
end

for (α, paths, seeds) in ((0.925, PATHS, (2, 3)), (0.85, PATHS * 2 ÷ 3, (4, 5)))
    MA = montecarlo(REF..., α, :P, paths; seed = seeds[1])
    MD = montecarlo(REF..., α, :D, paths; seed = seeds[2])
    nb = size(MA, 2)
    A = vec(sum(MA; dims = 2)) / nb
    B = A .- vec(sum(MD; dims = 2)) / nb
    rA = [ratios(MA[:, b]) for b in 1:nb]
    rd = [dep(MA[:, b], MA[:, b] .- MD[:, b]) for b in 1:nb]
    a0, p2, p4 = ratios(A)
    d = dep(A, B)
    println("\nα = $α: a₀ (TW/mm²), a₂/a₀, a₄/a₀ (%); departure (design − one quad off): Δa₀ (%), Δa₂/a₀, Δa₄/a₀ (pp)")
    @printf("  %-28s %.4f ± %.4f  %+.3f ± %.3f  %+.3f ± %.3f  | %+.3f ± %.3f  %+.3f ± %.3f  %+.3f ± %.3f\n",
            "Monte Carlo :hit", a0, se(getindex.(rA, 1)), p2, se(getindex.(rA, 2)), p4, se(getindex.(rA, 3)),
            d[1], se(getindex.(rd, 1)), d[2], se(getindex.(rd, 2)), d[3], se(getindex.(rd, 3)))
    println("  model rows: differences from the Monte Carlo (a₀ in %, moments in pp)")
    for (h, ch) in ((0.45, 0.15), (0.28, 0.09), (0.17, 0.05)), curved in (true, false)
        m = model(h, ch, curved)
        for closed in (false, true)
            mA, mB = m[(α, closed)]
            x, y = ratios(mA), dep(mA, mB)
            @printf("  %-28s %+.2f%%            %+.3f          %+.3f          | %+.3f          %+.3f          %+.3f\n",
                    @sprintf("%s %.2f/%.2f, %s", curved ? "curved" : "flat", h, ch, closed ? "balanced" : "raw"),
                    100 * (x[1] / a0 - 1), x[2] - p2, x[3] - p4, y[1] - d[1], y[2] - d[2], y[3] - d[3])
        end
    end
end
