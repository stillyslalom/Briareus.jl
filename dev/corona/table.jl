# Tables of the (b) vs (c) comparison from the results file written by compare.jl.
#   CORONA_OUT=<file> julia --project=. dev/corona/table.jl
using Serialization, Printf
d = deserialize(ENV["CORONA_OUT"])
ref = d["c_0.01_0.004_720"]["σl"]
ls = (2, 4, 6, 8, 10, 12, 16, 20, 24, 30)
rows = Any[]
for (n, Nb) in ((20, 2048), (30, 4096), (40, 2048))
    haskey(d, "b_n$(n)_Nb$(Nb)") || continue
    r = d["b_n$(n)_Nb$(Nb)"]
    n == 20 && push!(rows, ("hard sphere at r_c (all absorbed)", r["hard sphere r_c, all absorbed"]))
    n == 20 && push!(rows, ("hard sphere at r_c, η_a (Stage 0)", r["hard sphere r_c, eta_a (MN17)"]))
    push!(rows, ("(b) kernels, icosphere n = $n, Nb = $Nb", r["(b) kernel table"]))
end
for h in (0.08, 0.04, 0.02, 0.01)
    k = "c_$(h)_0.004_720"
    haskey(d, k) && push!(rows, ("(c) ray trace, hlat = $h mm", d[k]))
end
for (h, τ) in ((0.04, 0.002), (0.04, 0.008))
    k = "c_$(h)_$(τ)_720"
    haskey(d, k) && push!(rows, ("(c) hlat = $h, hτ = $τ", d[k]))
end
for e in (0.02, 0.04, -0.04)
    k = "p2_$(e)_0.04_0.004"
    haskey(d, k) && push!(rows, ("(c) P2 corona ε = $e", d[k]))
end
g(r, k) = r isa Dict ? r[k] : getproperty(r, Symbol(k))
@printf("%-36s %8s %8s", "model", "absorbed", "σ(l≤30)")
foreach(l -> @printf(" %7s", "σ_$l"), ls); @printf(" %10s\n", "max|Δσ_l|")
for (name, r) in rows
    σl = g(r, "σl")
    @printf("%-36s %8.5f %8.4f", name, g(r, "absorbed"), 100g(r, "σ"))
    foreach(l -> @printf(" %7.4f", 100σl[l+1]), ls)
    @printf(" %10.2e\n", 100maximum(abs.(σl[2:end] .- ref[2:end])))
end
println("(σ in % of the mean; max|Δσ_l| over 1 ≤ l ≤ 30 against (c) at hlat = 0.01, in % of the mean)")
for e in (0.02, 0.04, -0.04)
    k = "p2_$(e)_0.04_0.004"
    haskey(d, k) || continue
    r = d[k]; r0 = d["c_0.04_0.004_720"]
    @printf("P2 ε=%+.2f: Δc20/c00 %+.5f  Δc40/c00 %+.5f  Δσ_2 %+.4f %%  Δabsorbed %+.5f\n", e,
            r["a20"] - r0["a20"], r["a40"] - r0["a40"], 100(r["σl"][3] - r0["σl"][3]), r["absorbed"] - r0["absorbed"])
end
println("\nwork-precision")
for (n, Nb) in ((20, 2048), (30, 4096), (40, 2048))
    haskey(d, "b_n$(n)_Nb$(Nb)") || continue
    r = d["b_n$(n)_Nb$(Nb)"]; b = r["(b) kernel table"]
    @printf("(b) n=%d: %d points, operator %.2f s %.0f MB; kernels %.2f s; map %.2f s %.0f MB (%d pairs); evaluate %.3f s; grad 3 %.2f s; grad 80 %.2f s\n",
            n, r["points"], r["t_op"], r["opbytes"] / 2^20, r["t_K"], b["t_H"], b["bytes"] / 2^20, b["pairs"], b["t_ev"], b["t_g3"], b["t_g80"])
end
for h in (0.08, 0.04, 0.02, 0.01)
    k = "c_$(h)_0.004_720"
    haskey(d, k) && @printf("(c) hlat=%g: %d rays, %.1f s\n", h, d[k]["nr"], d[k]["t"])
end
haskey(d, "cgrad_0.04_0.004") && @printf("(c) gradient (3 ring shifts, hlat 0.04): %.1f s (value %.1f s)\n",
                                         d["cgrad_0.04_0.004"]["tg"], d["cgrad_0.04_0.004"]["t1"])
