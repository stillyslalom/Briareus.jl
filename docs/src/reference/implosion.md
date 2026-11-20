# Implosion metrics

Reduced models connect drive asymmetry to implosion observables: mode 1 as a hot-spot
velocity, the yield factors of the piston model, P2 swing metrics, and beam-error budgets
in velocity. Sources and confidence are in `reference/research-notes.md`, "Implosion
metrics".

```@autodocs
Modules = [Briareus]
Pages = ["src/implosion.jl"]
Private = false
```

## Capsule drive for a hydrodynamics code

`capsule_drive` turns a `drive_history` into the drive a capsule-only hydrodynamics
run takes (LLNL's two-step method): Tr(t), the capsule-averaged incident flux,
Legendre ratios aₗ/a₀ (P2 > 0 pole-hot) on a uniform ≤ 50 ps grid, the radius the flux
was evaluated on, and its sources. CompactLES takes it in memory through its
`CompactLESBriareusExt` extension (`radiation_drive(drive)`, loaded when both packages
are). The sources are `reference/litreview/hydro_coupling.md` and ROADMAP, "Briareus →
CompactLES".

```@autodocs
Modules = [Briareus]
Pages = ["src/capsule_drive.jl"]
Private = false
```
