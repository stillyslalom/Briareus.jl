# Young et al., PoP 27, 082702 (2020), Sec. II A: the static view-factor model (radii at
# 85 %, albedos 0.8/0.3) against the dynamic one (Fig. 3a's radii and albedos against
# time, on `Keyframes`) for N170827, with a proxy mode-1 source (N170827's delivered
# per-beam powers are not published): a fixed 5 % rms per-quad imbalance, then the
# pieces of the dynamic model one at a time. Research notes, "Young 2020 validation".
#
#     julia --project=. -t auto examples/young2020_dynamic.jl     # ~10 min on 4 threads
using Briareus, LinearAlgebra, Printf, Random
using Briareus: _interp, _hohlraum_grid
const H = parse(Float64, get(ENV, "H", "0.5"))
pub = Dict(:outer50 => (5.01, 0.0), :outer44 => (5.20, 0.0), :inner30 => (3.63, -0.60), :inner23 => (3.09, -1.09))
p = GroupPointing(nif_quads(power = 8.0); fields = (:z, :r))
quads = apply(p, reduce(vcat, [[pub[k][1] * 1.1188, pub[k][2]] for k in p.keys]))
P = young2020_pulses(); side(b) = startswith(String(b.group), "inner") ? :inner : :outer
pul = (; inner = PulseShape(P.inner.t, 4 .* P.inner.p), outer = PulseShape(P.outer.t, 4 .* P.outer.p))
R = young2020_radii(); A = young2020_albedos()
L, R0, r0, c0, rc0 = 11.3, 3.1, 1.82, 0.98, 1.05
wall(τ) = R0 - (R.wall[1] - _interp(R.t, R.wall, τ))
leh(τ) = r0 * _interp(R.t, R.leh, τ) / R.leh[1]
cap(τ) = c0 * _interp(R.t, R.capsule, τ) / R.capsule[1]
g = _hohlraum_grid(; length = L, radius = R0, leh_radius = r0, h = H, rc = rc0, nz = 36, nphi = 48, nr_cap = 6, nr_leh = 4, nr_corner = nothing)
dyn(τ) = hohlraum(; length = L, radius = wall(τ), leh_radius = leh(τ), capsule_radius = cap(τ), corner_radius = rc0 * wall(τ) / R0,
    g..., capsule_ntheta = 20, wall = gold(albedo = TabulatedAlbedo(A.t, A.wall), conversion = 0.9),
    capsule = ablator(albedo = TabulatedAlbedo(A.t, A.capsule)))
stat = hohlraum(; length = L, radius = 0.85R0, leh_radius = 0.85r0, capsule_radius = 0.85c0, corner_radius = 0.85rc0,
    g..., capsule_ntheta = 20, wall = gold(albedo = 0.8, conversion = 0.9))
times = collect(0.25:0.25:7.75)
mult = 1 .+ 0.05 .* randn(MersenneTwister(170827), length(quads))
ba = t -> [remake(b; power = mult[j] * b.power) for (j, b) in enumerate(with_pulses(quads, pul, t; by = side))]
K = young2020_kernel()
m1(h) = reduce(hcat, [100 .* Briareus.mode1_vector(view(h.ylm, :, k)) for k in axes(h.ylm, 2)])
pk = (times .>= 4.5) .& (times .<= 7.0)
report(name, h) = (m = m1(h); n = vec(sqrt.(sum(abs2, m; dims = 1)));
    @printf("%-40s |m| peak %.3f %%, v %.1f km/s, Tr peak %.0f eV\n", name, sum(n[pk]) / count(pk), norm(hotspot_velocity(h; kernel = K)), maximum(h.Tr)))
statm(t) = (m = RadiationModel(t; closed = true); (m, DepositionOperator(t, quads)))
report("static (0.85 radii, 0.8/0.3)", drive_history(statm(stat)..., ba, times; ylm = true, lmax = 1))
sdyn = hohlraum(; length = L, radius = 0.85R0, leh_radius = 0.85r0, capsule_radius = 0.85c0, corner_radius = 0.85rc0,
    g..., capsule_ntheta = 20, wall = gold(albedo = TabulatedAlbedo(A.t, A.wall), conversion = 0.9),
    capsule = ablator(albedo = TabulatedAlbedo(A.t, A.capsule)))
report("static radii, dynamic albedos", drive_history(statm(sdyn)..., ba, times; ylm = true, lmax = 1))
dynr(τ) = hohlraum(; length = L, radius = wall(τ), leh_radius = leh(τ), capsule_radius = cap(τ), corner_radius = rc0 * wall(τ) / R0,
    g..., capsule_ntheta = 20, wall = gold(albedo = 0.8, conversion = 0.9))
kf = Keyframes(dynr, times, quads; closed = true, keep = 2)
report("dynamic radii, static albedos", drive_history(kf, ba, times; ylm = true, lmax = 1))
dync(τ) = hohlraum(; length = L, radius = 0.85R0, leh_radius = 0.85r0, capsule_radius = cap(τ), corner_radius = 0.85rc0,
    g..., capsule_ntheta = 20, wall = gold(albedo = 0.8, conversion = 0.9))
kf2 = Keyframes(dync, times, quads; closed = true, keep = 2)
report("static radii but capsule shrinking, 0.8/0.3", drive_history(kf2, ba, times; ylm = true, lmax = 1))
kfd = Keyframes(dyn, times, quads; closed = true, keep = 2)
report("dynamic (Young Fig. 3a radii and albedos)", drive_history(kfd, ba, times; ylm = true, lmax = 1))
