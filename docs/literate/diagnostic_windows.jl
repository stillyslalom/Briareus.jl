# # Diagnostic windows in the hohlraum wall
#
# To see the capsule from the side during the implosion, a diagnostic needs a line of
# sight through the hohlraum wall. A hole would lose too much drive, so the wall
# carries a thin-foil window instead: a disk of gold about 2 mm across and a couple of
# micrometres thick, backed by a low-Z material that x rays above a few keV pass
# through. Early in the pulse the thin gold re-emits like the rest of the wall. The
# Marshak heat front that carries the x rays into the wall eats through it during the
# drive, and from then on the window lets radiation out: the waist of the capsule
# loses drive, and the poles, relatively, gain it. Two windows facing each other on
# the waist give a positive ``P_2`` and an ``m = 2`` asymmetry, and they cost drive.
# Assembly leaves a thin annular gap around each window that leaks from the start.
#
# This tutorial meshes circular windows with their gaps into the wall, checks the
# heat-front model that decides when a window burns through against published
# solutions and measurements, and then follows a pair of windows through a NIF-like
# drive for three window thicknesses, with and without a gap. It introduces
# `WallWindow` disks and gaps, the `BurnthroughAlbedo` thin-foil model with its
# `HeatFront`, and the window accounting of [`drive_history`](@ref) and
# [`energy_balance`](@ref).

using Briareus
using CairoMakie
using Printf
CairoMakie.activate!(type = "png", px_per_unit = 2)
include(joinpath(pkgdir(Briareus), "docs", "literate", "_plotting.jl"))   # figure helpers
set_theme!(fontsize = 15, linewidth = 2.5)

# ## A circular window on the wall
#
# A [`WallWindow`](@ref) with a `diameter` is a disk drawn on the unrolled wall: arc
# length ``s = R\phi`` along the wall by ``z`` along the axis. Unrolling a cylinder is
# an isometry, so the window's area is ``\pi d^2/4``. A `gap` adds a ring of that
# width around it, its own component (here the default, an open aperture).
#
# [`hohlraum`](@ref) cuts a square box around the window and its gap out of the
# conforming wall grid, half a wall cell wider than the gap, and re-meshes it with
# triangles: polar rings inside the window, one band of triangles across the gap, and a
# band joining the gap's outer edge to the vertices the box shares with the grid
# around it, so the mesh has no T-junctions. Every vertex lies on the cylinder. A
# triangle with its vertices on the cylinder has a curved area of exactly ``R`` times
# its area in ``(\phi, z)``; the radiation model already uses that exact area (see the
# [Integrating sphere](@ref) tutorial for why curved areas matter). The window's edge
# is a regular polygon whose radius is enlarged to [`equal_area_radius`](@ref), so the
# polygon's area is exactly that of the disk. The window and the gap therefore have
# exact areas, while their edges stay within a few micrometres of the circles.

d, gap = 2.0, 0.05
pair = [WallWindow(phi = φ, diameter = d, gap = gap) for φ in (78, 258)]
target = hohlraum(length = 10.1, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9,
                  windows = pair)
area(name) = sum(target.elem_area[elements(target, name)])
exact = (window = π * d^2 / 4, gap = π * ((d / 2 + gap)^2 - (d / 2)^2))
for (name, ref) in ((:window_1, exact.window), (:window_1_gap, exact.gap),
                    (:window_2, exact.window), (:window_2_gap, exact.gap))
    @printf("%-14s %4d elements  area %.12f mm²  exact %.12f  relative error %.1e\n", name,
            length(elements(target, name)), area(name), ref, area(name) / ref - 1)
end
R, L = 2.875, 10.1
@printf("wall + windows + gaps: %.12f mm², 2πRL = %.12f\n",
        area(:wall) + area(:window_1) + area(:window_1_gap) + area(:window_2) + area(:window_2_gap),
        2π * R * L)

# The unrolled mesh around the first window. The box's outer vertices are the grid's
# own, the gap is one band of triangles 50 µm wide, and the window's elements are about
# the size of the wall's:

function unrolled(target, φc)
    fig = Figure(size = (800, 640))
    ax = Axis(fig[1, 1], aspect = DataAspect(), xlabel = "arc length from the window center (mm)",
              ylabel = "z (mm)", title = "Window, gap and wall near φ = $(φc)°")
    colors = Dict(:wall => (:gray85, "wall"), :window_1 => (:goldenrod1, "window"),
                  :window_1_gap => (:dodgerblue3, "gap"))
    uv(p) = Point2f(R * (mod(atan(p[2], p[1]) - deg2rad(φc) + π, 2π) - π), p[3])
    for k in eachindex(target.tri_area)
        name = target.components[target.tri_comp[k]].name
        haskey(colors, name) || continue
        v = target.tri_v0[k]
        pts = uv.((v, v + target.tri_e1[k], v + target.tri_e2[k]))
        all(p -> abs(p[1]) < 1.6 && abs(p[2]) < 1.6, pts) || continue
        poly!(ax, collect(pts); color = colors[name][1], strokecolor = :gray20, strokewidth = 0.6)
    end
    limits!(ax, -1.6, 1.6, -1.6, 1.6)
    Legend(fig[1, 2], [PolyElement(color = c, strokecolor = :gray20, strokewidth = 0.6)
                       for (c, _) in values(colors)], [l for (_, l) in values(colors)])
    return fig
end
unrolled(target, 78)

# ## When a thin gold foil burns through
#
# Radiation of temperature ``T`` drives a subsonic, ablative heat front into a gold
# wall. Hammer and Rosen solved this Marshak wave self-similarly with power-law fits to
# gold's opacity and equation of state (Phys. Plasmas 10, 1829, 2003). At constant
# ``T`` the areal mass behind the front is
#
# ```math
# m_F = 9.51\times10^{-4}\,\left(\frac{T}{100\ \mathrm{eV}}\right)^{1.914}
#       \left(\frac{t}{1\ \mathrm{ns}}\right)^{0.5156}\ \mathrm{g/cm^2}
# ```
#
# (their eq. 112 for a constant drive), the default [`HeatFront`](@ref). A window
# burns through when ``m_F`` reaches its own areal mass, ``19.3\ \mathrm{g/cm^3}``
# times its thickness, 3.86 mg/cm² for 2 µm. A hohlraum's temperature is far from
# constant, so the front is integrated over the drive history:
#
# ```math
# m_F(t) = c\left[\int_0^t \left(\frac{T(t')}{100\ \mathrm{eV}}\right)^{a/b} dt'\right]^b ,
# \qquad c = 9.51\times10^{-4},\ a = 1.914,\ b = 0.5156 .
# ```
#
# This depends only on the drive history and is exact at constant ``T``. For a ramp
# ``T = T_0 t^k`` it gives ``m_F = c\,(1 + ka/b)^{-b}\,T(t)^a t^b``. The exact
# self-similar solutions for ramps are independent references for it.
#
# [`heat_front`](@ref) integrates the drive by the trapezoidal rule. At constant
# temperature on an uneven grid it reproduces the closed form:

front = HeatFront()
(; coefficient, temperature_exponent, time_exponent) = front
ts = [0.0, 0.05, 0.3, 0.31, 1.0, 2.2, 4.0]
m = heat_front(front, ts, fill(300.0, length(ts)))
@printf("constant 300 eV: largest relative error %.1e\n",
        maximum(abs.(m[2:end] ./ (coefficient * 3.0^temperature_exponent .* ts[2:end] .^ time_exponent) .- 1)))

# For ramps, compare with two exact solutions:
# - Hammer and Rosen give ``9.51\times10^{-4}`` at constant temperature and
#   ``7.13\times10^{-4}`` for ``T \propto t^{0.182}`` (their eq. 112, ``q_0 = 1``).
# - Shussman and Heizler's exact solution of the same problem (Phys. Plasmas 22,
#   082109, 2015, eqs. 35–36) gives ``10.17\times10^{-4}`` at constant temperature and
#   ``8.32\times10^{-4}`` at constant absorbed flux, ``T \propto t^{0.123}``.
#
# The two solutions differ by 7% at constant temperature, so the test is the ratio of
# each ramp to its own constant-temperature coefficient. Hammer and Rosen also plot
# the ablated mass at 4 ns for ``T_s = (t/1\,\mathrm{ns})^{0.182}`` hectovolts: about
# 2.2–2.4 mg/cm² from both HYDRA and the analytic curve (their Fig. 4b).

ramp_mass(k, T0, t) = heat_front(front, range(0, t; length = 40001),
                                 T0 .* range(0, t; length = 40001) .^ k)[end]
for (k, exact_ratio, source) in ((0.182, 7.13 / 9.51, "Hammer & Rosen"),
                                 (0.123, 8.32 / 10.17, "Shussman & Heizler"))
    ratio = ramp_mass(k, 100.0, 1.0) / ramp_mass(0.0, 100.0, 1.0)
    @printf("T ∝ t^%.3f: m(k)/m(0) = %.4f, exact %.4f (%s): %+.1f%%\n", k, ratio, exact_ratio,
            source, 100 * (ratio / exact_ratio - 1))
end
@printf("ablated mass at 4 ns, T = t^0.182 heV: %.2f mg/cm² (HYDRA and analytic: 2.2–2.4)\n",
        1e3 * ramp_mass(0.182, 100.0, 4.0))

# The integrated front is 2.2% too heavy for Hammer and Rosen's ramp and 0.7% for
# Shussman and Heizler's. The residual has a known cause: the shape of the
# self-similar profile depends on the drive history, which a single integral of the
# temperature cannot carry. Two percent in mass is four percent in burnthrough time.
#
# ### Against measured burnthrough
#
# The Nova experiments measured burnthrough of gold foils of 0.8–3 µm on the wall of
# a hohlraum near 260 eV. In those experiments the radiation temperature rose as about
# ``t^{0.15}`` under a flat laser pulse (Lindl et al., Phys. Plasmas 11, 339, 2004,
# Fig. 4-6, which compares the data with a heat-front curve for that drive). The times
# below are the half-signal times read from that figure, to about ±0.03 ns.
# [`burnthrough_time`](@ref) gives the heat front's arrival for the same drive, and
# [`BurnthroughAlbedo`](@ref) starts its transition there, reaching half transmission
# `width/2` (0.2 ns by default) later:

nova = [(0.77, 0.33), (1.07, 0.565), (1.55, 0.78), (1.93, 0.86), (2.05, 1.03), (2.27, 1.14),
        (3.0, 1.52)]
tn = range(0, 2.5; length = 2501)
Tnova = 260 .* tn .^ 0.15
println("  Au (µm)   measured (ns)   front arrival   + width/2   constant 260 eV")
offsets = Float64[]
for (x, measured) in nova
    model = BurnthroughAlbedo(thickness = x)
    tb = burnthrough_time(model, tn, Tnova)
    push!(offsets, measured - tb)
    @printf("  %5.2f     %5.2f           %5.2f           %5.2f       %5.2f\n", x, measured, tb,
            tb + model.width / 2, burnthrough_time(model, [0.0, 2.5], [260.0, 260.0]))
end
@printf("measured − front arrival: %.2f ns on average (%.2f–%.2f)\n",
        sum(offsets) / length(offsets), extrema(offsets)...)

# The front tracks the data's thickness dependence: about 0.8 ns from 0.8 to 2.3 µm,
# and 1.2 ns out to 3 µm. The measured half signal comes 0.1–0.25 ns after the front
# arrives, and the default transition width puts half transmission there. A
# constant-temperature estimate would burn the foils through three times too early.
#
# Before the front arrives, the foil absorbs exactly what a thick wall does: in the
# diffusion picture the cold gold ahead of a nonlinear front does not change the flux
# into it, and the measurements show no early transmission. So the thin-foil albedo
# follows the wall's own Marshak scaling, [`WallLossAlbedo`](@ref), until burnthrough.
# Then it falls smoothly to its `after` value, 0 by default: everything that reaches
# the burned-through foil leaves the hohlraum. That is the open-hole limit, an upper
# bound on the loss.

# ## A pair of windows in a NIF-like hohlraum
#
# The hohlraum of the [Freeing a bundle for radiography](@ref) tutorial (10.1 × 5.75 mm,
# 3.1 mm entrance holes, a 0.9 mm capsule) has a gold wall with the Marshak albedo of
# [`WallLossAlbedo`](@ref) and a laser-to-x-ray conversion of 0.8. Its 48 quads are
# pointed as there. Each quad follows a pulse with a picket, a long foot, a second step
# and a 2 ns peak at 9.5 TW per quad (456 TW, 1.6 MJ in all). Two windows 2 mm across
# sit at the waist at φ = 78° and 258°, a pair related by the 180° rotation about the
# axis, which the NIF ports share; their foils are 1.7, 2.5 or 3.6 µm of gold.
# The mesh is coarser than the default to keep the tutorial quick.

wall = Material(:Au; albedo = WallLossAlbedo(), conversion = 0.8)
resolution = (; nz = 24, nphi = 32, capsule_ntheta = 12, capsule_nphi = 24)
build(windows) = hohlraum(; length = 10.1, radius = 2.875, leh_radius = 1.55,
                          capsule_radius = 0.9, wall, windows, resolution...)
foil(thickness) = Material(:Au; albedo = BurnthroughAlbedo(; thickness), conversion = 0.8)
windows(thickness; gap = 0.0) =
    [WallWindow(; phi, diameter = 2.0, gap, material = foil(thickness)) for phi in (78, 258)]

pointing = Dict(:outer50 => (0.0, 5.01), :outer44 => (0.0, 5.20),
                :inner30 => (0.60, 3.63), :inner23 => (1.09, 3.09))
function aimed(b)
    r, z = pointing[b.group]
    φ = port_angles(b)[2]
    return remake(b; aim = (-r * cosd(φ), -r * sind(φ), sign(b.port[3]) * z))
end
quads = aimed.(nif_quads())
pulse = PulseShape([0.0, 0.3, 1.0, 1.4, 4.5, 5.0, 6.0, 7.5, 9.5, 10.0],
                   [0.0, 1.2, 1.2, 0.3, 0.3, 1.0, 1.0, 9.5, 9.5, 0.0])
@printf("laser energy %.2f MJ, peak power %.0f TW\n", 48 * pulse_energy(pulse) / 1e3, 48 * 9.5)

# [`drive_history`](@ref) marches the radiosity through the pulse quasi-statically.
# Each window element tracks its own heat front from its own incident flux, so the
# windows burn through when the flux on them says they should. With `ylm = true` it
# also records the real spherical harmonics of the capsule flux. The three thicknesses
# share one geometry, so they share its view factors and its deposition operator.
#
# The steps are 0.2 ns, but a window's albedo falls to zero over the 0.4 ns `width` of
# its transition, which two steps would jump across. With `transition = 0.05` the march
# inserts steps wherever an element's transition would advance by more than a
# twentieth of the way in one step, planned from each element's known or predicted
# burnthrough time. It returns the times it took in `time`.

times = 0.2:0.2:9.8
beams_at(τ) = with_pulses(quads, pulse, τ)
history(model, op, times) = drive_history(model, op, beams_at, times; lmax = 6, ylm = true,
                                          transition = 0.05)
runs = Dict{Tuple{Float64,Float64},Any}()
shared = Dict{Float64,Any}()                       # view factors and operator per gap
for g in (0.0, 0.05), thickness in (1.7, 2.5, 3.6)
    t = build(windows(thickness; gap = g))
    if !haskey(shared, g)
        shared[g] = (RadiationModel(t; closed = true).F, DepositionOperator(t, quads))
    end
    F, op = shared[g]
    runs[(thickness, g)] = history(RadiationModel(t; F), op, times)
end
@printf("steps: %d given, %d–%d taken\n", length(times),
        extrema(length(h.time) for h in values(runs))...)

# The hohlraum without windows is the reference. It is solved at every time any of the
# windowed runs took, so the differences below compare like with like.

base_target = build(())
base = history(RadiationModel(base_target; closed = true), DepositionOperator(base_target, quads),
               sort(unique(reduce(vcat, [h.time for h in values(runs)]))))
@printf("reference steps: %d\n", length(base.time))

# The change in the capsule's radiation temperature and in ``P_2`` from the hohlraum
# without windows, the ``m = 2`` amplitude ``\sqrt{c_{22}^2 + c_{2,-2}^2}/c_{00}`` of
# the capsule flux, and each window's heat front as a fraction of its foil. Without
# windows the capsule's radiation temperature peaks at the value printed below:

on(h, x) = x[indexin(h.time, base.time)]              # base values at the times of run h
P2(h) = h.moments[3, :] ./ h.moments[1, :]
P4(h) = h.moments[5, :] ./ h.moments[1, :]
m2(h) = hypot.(h.ylm[mode_index(2, 2), :], h.ylm[mode_index(2, -2), :]) ./ h.ylm[1, :]
palette = Makie.wong_colors()
colors = Dict(1.7 => palette[6], 2.5 => palette[2], 3.6 => palette[1])
fig = Figure(size = (820, 700))
axT = Axis(fig[1, 1], ylabel = "ΔTᵣ (eV)", title = "Change in capsule Tᵣ")
axP = Axis(fig[1, 2], ylabel = "ΔP₂/P₀ (%)", title = "Change in P₂ from the windows")
axM = Axis(fig[2, 1], xlabel = "time (ns)", ylabel = "|c₂,±₂|/c₀₀ (%)", title = "m = 2 amplitude")
axF = Axis(fig[2, 2], xlabel = "time (ns)", ylabel = "heat front / foil areal mass",
           title = "Heat front in window 1")
hlines!(axT, [0.0], color = :black)
hlines!(axP, [0.0], color = :black)
for thickness in (1.7, 2.5, 3.6), g in (0.0, 0.05)
    h = runs[(thickness, g)]
    kw = (; color = colors[thickness], linestyle = g == 0 ? :solid : :dash)
    lines!(axT, h.time, h.Tr .- on(h, base.Tr); kw...)
    lines!(axP, h.time, 100 .* (P2(h) .- on(h, P2(base))); kw...)
    lines!(axM, h.time, 100 .* m2(h); kw...)
    g == 0 && lines!(axF, h.time, h.front[:window_1]; kw...)
end
hlines!(axF, [1.0], color = :gray50, linestyle = :dot)
ylims!(axF, 0, 2)
Legend(fig[3, 1:2],
       [[LineElement(color = :black)];
        [LineElement(color = colors[x]) for x in (1.7, 2.5, 3.6)];
        [LineElement(color = :gray30), LineElement(color = :gray30, linestyle = :dash)]],
       ["no windows (reference)", "1.7 µm", "2.5 µm", "3.6 µm", "no gap", "50 µm open gap"],
       orientation = :horizontal, framevisible = false, nbanks = 1)
@printf("capsule Tᵣ without windows: %.0f eV at the peak\n", maximum(base.Tr))
fig

# Without a gap the windows cost nothing until they burn through: an intact window
# re-emits like the wall, and its ``\Delta P_2`` and ``m = 2`` stay below 0.05%, the
# size of the remeshing. Then ``P_2`` jumps by about 5% (the waist goes cold), the
# ``m = 2`` amplitude reaches about 4%, and the capsule's radiation temperature falls
# by 14 eV. Thinner foils burn through earlier: 1.7 µm during the second step of the
# pulse, 2.5 µm on the rise, and 3.6 µm just before the peak. An open 50 µm gap leaks
# from the start. Before burnthrough the gaps alone give about 0.5% of ``P_2``, 0.4% of
# ``m = 2`` and 1.3 eV; at the peak they add 0.6% and 1.5 eV to the windows' cost.
#
# A capsule's shape integrates the drive asymmetry over the pulse. The table therefore
# gives the change in ``P_2``, ``P_4`` and ``m = 2`` weighted by the capsule flux from
# 5 ns on, as well as at the peak. It also gives the drive cost: the radiation
# temperature at the peak, the capsule's absorbed energy over the pulse, and the
# energy that leaves through the windows and gaps (`window_loss`). Each history is
# summed over its own steps, each step weighted by half the time to its neighbors (the
# first and last by a whole step), which on the uniform `times` is the plain sum times
# the step.

function weights(τ)
    w = [(τ[min(k + 1, end)] - τ[max(k - 1, 1)]) / 2 for k in eachindex(τ)]
    w[1] *= 2
    w[end] *= 2
    return w
end
integral(h, x) = sum(weights(h.time) .* x)
function flux_weighted(h, x; from = 5.0)
    w = weights(h.time) .* (h.time .>= from) .* h.moments[1, :]
    return sum(w .* x) / sum(w)
end
function peak_mean(h, x)
    w = weights(h.time) .* (7.6 .<= h.time .<= 9.4)
    return sum(w .* x) / sum(w)
end
println("  Au    gap   burnthrough  ΔTᵣ peak  capsule energy  window loss   ΔP₂     ΔP₂     ΔP₄      m=2")
println(" (µm)  (µm)     (ns)        (eV)       (%)           (kJ)       peak   weighted weighted weighted")
for g in (0.0, 0.05), thickness in (1.7, 2.5, 3.6)
    h = runs[(thickness, g)]
    tb = h.burnthrough[:window_1]
    @printf(" %.1f   %3.0f     %5.2f      %+6.1f      %+6.1f        %6.0f       %+5.2f%%  %+5.2f%%  %+5.2f%%   %4.2f%%\n",
            thickness, 1e3g, sum(tb) / length(tb), peak_mean(h, h.Tr) - peak_mean(base, base.Tr),
            100 * (integral(h, h.flux) / integral(base, base.flux) - 1), integral(h, h.window_loss),
            100 * (peak_mean(h, P2(h)) - peak_mean(base, P2(base))),
            100 * (flux_weighted(h, P2(h)) - flux_weighted(base, P2(base))),
            100 * (flux_weighted(h, P4(h)) - flux_weighted(base, P4(base))),
            100 * flux_weighted(h, m2(h)))
end

# At the peak all three thicknesses have burned through and look alike. What the
# thickness changes is how long the windows have been open. The flux-weighted
# ``\Delta P_2`` rises from 4.1% for 3.6 µm to 5.3% for 1.7 µm, and the capsule loses
# 12–15% of its energy. The pair of windows adds no odd-``m`` harmonics, because it
# shares the ports' 180° symmetry:

h = runs[(2.5, 0.05)]
odd = maximum(abs(h.ylm[mode_index(l, m), k]) / h.ylm[1, k]
              for l in 1:4 for m in -l:l if isodd(m) for k in eachindex(h.time))
@printf("largest odd-m |c_lm|/c₀₀ over the pulse: %.1e\n", odd)

# The capsule flux at the end of the peak, as equal-area (Hammer) maps seen from
# outside, centered between the windows, which face the capsule's waist at φ = 78° and
# 258° (circles). The maps are drawn from the spherical harmonics through ``l = 6``,
# which carry the low modes a capsule responds to, rather than from the coarse capsule
# mesh. The top row shows the flux relative to its own mean, without windows and with
# the 2.5 µm pair and its gaps. Without windows this untuned pointing already drives
# the waist hot, so the bottom row subtracts that map and shows what the windows
# change: on the left without gaps, on the right with them.

function ylm_field(c; lmax = 6)
    C = ylm_matrix(c, lmax)
    Y = zeros(lmax + 1, 2lmax + 1)
    mean = c[1] / sqrt(4π)                          # c₀₀ Y₀₀
    return (θ, φ) -> (Briareus._real_ylm!(Y, lmax, cosd(θ), deg2rad(φ)); 100 * (sum(C .* Y) / mean - 1))
end
φ0 = 168.0                                         # map center, between the windows
deviation_map(h) = sphere_map(ylm_field(h.ylm[:, findfirst(≈(9.4), h.time)]); φ0)
base_map = deviation_map(base)
change_map(h) = (m = deviation_map(h); (m[1], m[2], m[3] .- base_map[3]))
rows = (["no windows" => base_map, "2.5 µm windows, 50 µm gaps" => deviation_map(h)],
        ["2.5 µm windows, no gap" => change_map(runs[(2.5, 0.0)]),
         "2.5 µm windows, 50 µm gaps" => change_map(h)])
limits = [ceil(maximum(maximum(abs, filter(isfinite, m[3])) for (_, m) in row)) for row in rows]
fig = Figure(size = (800, 470))
for (i, row) in enumerate(rows), (j, (title, m)) in enumerate(row)
    ax = sphere_axis(fig[i, j]; title, φ0)
    sphere_heatmap!(ax, m; colormap = :balance, colorrange = (-limits[i], limits[i]))
    mark_directions!(ax, [(90, 78), (90, 258)]; φ0, marker = :circle, markersize = 14,
                     color = :transparent, strokecolor = :black, strokewidth = 2)
end
Colorbar(fig[1, 3]; colormap = :balance, colorrange = (-limits[1], limits[1]),
         height = Relative(0.8), label = "H/⟨H⟩ − 1 (%)")
Colorbar(fig[2, 3]; colormap = :balance, colorrange = (-limits[2], limits[2]),
         height = Relative(0.8), label = "change from no windows (%)")
rowgap!(fig.layout, 6)
fig

# Without windows this untuned pointing drives the waist hot (``P_2/P_0 \approx -12\%``,
# which cone tuning would remove). With the windows open, the waist facing them turns
# cold while the waist between them gains, relative to the mean, and the poles gain
# too. The waist's mean excess drops, which is the positive ``\Delta P_2``, and the
# pattern ``\cos 2(\phi - 78°)`` is the ``m = 2`` asymmetry. The gaps add to it a little.

# ### The window energy balance
#
# At any time, [`energy_balance`](@ref) splits a thin-foil window's loss into the part
# it absorbs as if it were intact, ``(1 - \alpha_\mathrm{before}) H A``, and the excess
# it transmits once burned through, ``(\alpha_\mathrm{before} - \alpha) H A``. With the
# open gap's escape this is `window_loss`. At the end of the peak, for 2.5 µm windows
# with a 50 µm gap:

t = build(windows(2.5; gap = 0.05))
F, op = shared[0.05]
model = RadiationModel(t; F)
tb = fill(NaN, nelements(t))
for n in (:window_1, :window_2)
    tb[elements(t, n)] .= h.burnthrough[n]
end
sol = solve(model, deposit(op, beams_at(9.4)); time = 9.4, burnthrough = tb)
bal = energy_balance(sol)
@printf("source %.1f TW: absorbed %.1f, entrance holes %.1f, windows transmit %.1f, gaps %.1f, lost %.1e\n",
        bal.source, sum(values(bal.absorbed)), bal.escaped[:leh_top] + bal.escaped[:leh_bottom],
        bal.transmitted[:window_1] + bal.transmitted[:window_2],
        bal.escaped[:window_1_gap] + bal.escaped[:window_2_gap],
        bal.source - sum(values(bal.absorbed)) - sum(values(bal.escaped)) - sum(values(bal.transmitted)))
@printf("window loss from the history at 9.4 ns: %.1f TW\n", h.window_loss[findfirst(≈(9.4), h.time)])

# ## What a gray view-factor model leaves out
#
# The model gives the geometry of the loss exactly but not its physics after
# burnthrough. Setting `after = 0` treats a burned-through window as a hole. A real
# foil tamped by its low-Z backing stays hot and goes on re-emitting, and the heated
# backing radiates back. So the losses and asymmetries above are upper bounds once the
# windows burn through: an `after` albedo between 0 and the wall's brackets them. The
# gap is open from the start here, but an assembly gap tens of micrometres wide closes
# as its walls blow off, so the gap rows overstate its late-time leak. The view-factor
# model also has no capsule hydrodynamics: it gives flux asymmetries, not the shape
# they imprint.

# ## What this checks
#
# - A circular window and its gap have the areas ``\pi d^2/4`` and
#   ``\pi[(d/2 + g)^2 - (d/2)^2]`` to round-off. The equal-area edge polygon and the
#   exact curved areas of triangles on the cylinder make them exact, and wall,
#   windows and gaps tile the cylinder exactly.
# - The integrated heat front reproduces the constant-temperature closed form to
#   round-off. For ramps it agrees with the exact self-similar solutions of Hammer
#   and Rosen to 2.2% and of Shussman and Heizler to 0.7%, and with Hammer and
#   Rosen's simulated ablated mass. The cause of the residual is identified.
# - Over 0.8–3 µm, the heat front and the transition width reproduce the measured
#   burnthrough times of gold foils on a 260 eV hohlraum wall to within the scatter
#   of the data. A constant-temperature estimate does not.
# - The thin-foil energy accounting closes to round-off at a single time, and it
#   agrees with the window loss the drive history tracks.
# - A pair of windows related by the 180° rotation about the axis adds no odd-``m``
#   harmonics.
#
# The test suite (`test/windows.jl`, `test/burnthrough.jl`) adds four more checks.
# The view-factor rows of a windowed hohlraum close as well as those without windows.
# Reciprocity holds. Each window element's burnthrough time in the drive history
# matches the standalone estimate from that element's own temperature history. A
# history refined through the transitions follows a uniformly fine one five times more
# closely than the unrefined steps do.
