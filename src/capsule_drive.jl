# The drive handed to a capsule-only hydrodynamics run (ROADMAP 17, LLNL's two-step
# method; reference/litreview/hydro_coupling.md §(c)).

"""
    CapsuleDrive

The radiation drive on the capsule over time, for a capsule-only hydrodynamics run
(LLNL's two-step method: an integrated hohlraum step, then capsule-only runs driven by
its capsule flux; Clark et al., PoP 26, 050601 (2019); Masse et al., PoP 26 (2019)).
Build one with [`capsule_drive`](@ref). Fields:

- `time` — uniform time grid (ns)
- `Tr` — radiation temperature of the capsule-averaged incident flux (eV),
  `radiation_temperature(flux)`
- `flux` — capsule-averaged incident flux `a₀` (TW/mm²; ×1e14 → W/cm²)
- `modes` — `(lmax+1) × nt` Legendre ratios `aₗ/a₀` (row `l+1`, so row 1 is 1), of
  the incident flux `F(θ,t) = a₀(t) Σₗ (aₗ/a₀)(t) Pₗ(cos θ)`, `θ` from the hohlraum axis
  (+z) about the capsule centre. **`P2 > 0` is pole-hot** (Masse 2019's convention).
- `absorbed` — x-ray power absorbed by the capsule (TW), as in `energy_balance`
  (`NaN` where the history does not record it)
- `radius` — the radius (mm) of the sphere the flux was evaluated on (the capsule mesh)
- `area` — the observed elements' area (mm²)
- `provenance` — `Dict{Symbol,Any}`: code, version, git commit, the history's source
  times, observed components, capsule material, mesh size, and whatever the caller adds

The angular structure is that of the flux on the initial capsule radius; applying it
at a moving ablation front far inside that radius changes the geometric smoothing
(hydro_coupling.md §2, pitfall 2). Gray: no spectral (M-band) split yet (ROADMAP 13).
"""
struct CapsuleDrive
    time::Vector{Float64}
    Tr::Vector{Float64}
    flux::Vector{Float64}
    modes::Matrix{Float64}
    absorbed::Vector{Float64}
    radius::Float64
    area::Float64
    provenance::Dict{Symbol,Any}
    function CapsuleDrive(time, Tr, flux, modes, absorbed, radius, area, provenance)
        nt = length(time)
        (length(Tr) == nt && length(flux) == nt && length(absorbed) == nt && size(modes, 2) == nt) ||
            throw(ArgumentError("time, Tr, flux, absorbed and the columns of modes differ in length"))
        size(modes, 1) >= 1 || throw(ArgumentError("modes needs at least the l = 0 row"))
        issorted(time) || throw(ArgumentError("times must be sorted"))
        radius > 0 || throw(ArgumentError("radius must be positive"))
        new(collect(Float64, time), collect(Float64, Tr), collect(Float64, flux),
            Matrix{Float64}(modes), collect(Float64, absorbed), Float64(radius),
            Float64(area), Dict{Symbol,Any}(provenance))
    end
end

"Highest Legendre order a `CapsuleDrive` carries."
_lmax(d::CapsuleDrive) = size(d.modes, 1) - 1

function Base.show(io::IO, d::CapsuleDrive)
    print(io, "CapsuleDrive(", length(d.time), " times on [", first(d.time), ", ", last(d.time),
          "] ns, lmax = ", _lmax(d), ", R = ", d.radius, " mm, peak Tr = ",
          round(maximum(d.Tr); digits = 1), " eV)")
end

"""
    (d::CapsuleDrive)(t) -> (; Tr, flux, modes)

The drive at time `t` (ns), linear in the moments `aₗ` between grid times and held at
the end values outside them.
"""
function (d::CapsuleDrive)(t::Real)
    a = [_interp(d.time, view(d.modes, l + 1, :) .* d.flux, t) for l in 0:_lmax(d)]
    flux = a[1]
    return (; Tr = radiation_temperature(flux), flux, modes = flux > 0 ? a ./ flux : [1; zeros(_lmax(d))])
end

"""
    capsule_drive(history, target; observe = :capsule, dt = 0.05, times = nothing,
                  lmax = size(history.moments, 1) - 1, radius = nothing,
                  provenance = (;)) -> CapsuleDrive

The capsule drive of a [`drive_history`](@ref) result for a capsule-only hydrodynamics
run, on a uniform grid of spacing `≤ dt` (ns; 50 ps by default, hydro_coupling.md §(c))
spanning the history's solved times, or on the given `times` (within them). The
Legendre moments `aₗ` are interpolated linearly between the history's steps and
divided by the interpolated `a₀`; `flux` is `a₀` and `Tr` its radiation temperature,
exactly the history's `flux` and `Tr` at its own times. Run `drive_history(...; lmax = 8)`
for the moments through `l = 8` the hand-over calls for; `lmax` here may only truncate.

`target` (a `Target` or `RadiationModel`; the keyframes' reference target for a moving
hohlraum) supplies the radius of the observed sphere (its `SphereSurface`; pass
`radius` for any other surface), the observed area, and the provenance of the capsule
material and mesh. `observe` must name the components the history observed.
`provenance` (NamedTuple or Dict: hohlraum, pulse, albedo model, …) is merged into the
recorded provenance.
"""
function capsule_drive(h::NamedTuple, target; observe = :capsule, dt::Real = 0.05,
                       times = nothing, lmax::Integer = size(h.moments, 1) - 1,
                       radius = nothing, provenance = (;))
    t = target isa RadiationModel ? target.target : target isa Keyframes ? target.reference : target
    t isa Target || throw(ArgumentError("target must be a Target, RadiationModel or Keyframes"))
    names = observe isa Symbol ? (observe,) : Tuple(observe)
    elems = elements(t, names...)
    isempty(elems) && throw(ArgumentError("no elements in $(names)"))
    0 <= lmax <= size(h.moments, 1) - 1 ||
        throw(ArgumentError("the history carries moments through l = $(size(h.moments, 1) - 1); " *
                            "run drive_history(...; lmax = $lmax)"))
    th = h.time
    length(th) >= 2 || throw(ArgumentError("the history needs at least two times"))
    if times === nothing
        dt > 0 || throw(ArgumentError("dt must be positive"))
        n = max(ceil(Int, (th[end] - th[1]) / dt - 1e-9), 1)
        grid = collect(range(th[1], th[end]; length = n + 1))
    else
        grid = collect(Float64, times)
        issorted(grid) || throw(ArgumentError("times must be sorted"))
        (grid[1] >= th[1] - 1e-12 && grid[end] <= th[end] + 1e-12) ||
            throw(ArgumentError("times must lie within the history's [$(th[1]), $(th[end])] ns"))
    end
    nt = length(grid)
    A = zeros(lmax + 1, nt)
    for l in 0:lmax, k in 1:nt
        A[l+1, k] = _interp(th, view(h.moments, l + 1, :), grid[k])
    end
    flux = A[1, :]
    modes = A ./ reshape(map(x -> x > 0 ? x : Inf, flux), 1, :)
    modes[1, :] .= 1
    absorbed = length(names) == 1 && haskey(h.absorbed, names[1]) ?
               [_interp(th, h.absorbed[names[1]], τ) for τ in grid] : fill(NaN, nt)
    R = radius === nothing ? _observed_radius(t, elems) : Float64(radius)
    area = sum(t.elem_area[elems])
    prov = Dict{Symbol,Any}(
        :code => "Briareus", :version => string(pkgversion(@__MODULE__)), :git => _git_commit(),
        :quantity => "incident x-ray flux on the capsule (gray), Legendre about +z, P2 > 0 pole-hot",
        :units => "time ns, Tr eV, flux TW/mm^2, absorbed TW, radius mm, area mm^2",
        :observe => collect(Symbol, names), :history_times => (th[1], th[end], length(th)),
        :interpolation => "linear in a_l between history steps",
        :capsule_material => string(element_material(t, elems[1]).albedo),
        :elements => length(elems), :target_elements => nelements(t),
        :components => [c.name for c in t.components])
    for (k, v) in pairs(provenance)
        prov[Symbol(k)] = v
    end
    return CapsuleDrive(grid, radiation_temperature.(flux), flux, modes, absorbed, R, area, prov)
end

# The radius of the sphere the observed elements lie on.
function _observed_radius(t::Target, elems)
    s = t.components[t.elem_comp[elems[1]]].surface
    for e in elems
        se = t.components[t.elem_comp[e]].surface
        (se isa SphereSurface && se.radius == s.radius && se.center == s.center) ||
            throw(ArgumentError("the observed elements are not on one sphere; pass `radius`"))
    end
    return s.radius
end

# The package checkout's commit (with "-dirty" for uncommitted changes), or "unknown".
function _git_commit()
    dir = pkgdir(@__MODULE__)
    try
        sha = readchomp(pipeline(`git -C $dir rev-parse --short HEAD`; stderr = devnull))
        dirty = !isempty(readchomp(pipeline(`git -C $dir status --porcelain --untracked-files=no`;
                                            stderr = devnull)))
        return dirty ? sha * "-dirty" : sha
    catch
        return "unknown"
    end
end
