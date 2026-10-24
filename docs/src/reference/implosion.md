# Implosion metrics

Reduced models from drive asymmetry to implosion observables: mode 1 as a hot-spot
velocity, the piston model's yield factors, P2 swing metrics, and beam-error budgets in
velocity. Sources and confidence: `reference/research-notes.md`, "Implosion metrics".

```@autodocs
Modules = [Briareus]
Pages = ["src/implosion.jl"]
Private = false
```

## Capsule drive for a hydrodynamics code

`capsule_drive` turns a `drive_history` into the drive a capsule-only hydrodynamics
run takes (LLNL's two-step method): Tr(t), the capsule-averaged incident flux,
Legendre ratios aₗ/a₀ (P2 > 0 pole-hot) on a uniform ≤ 50 ps grid, the radius the flux
was evaluated on, and provenance. CompactLES reads it in memory through its
`CompactLESBriareusExt` extension (`radiation_drive(drive)`, loaded when both packages
are). Sources: `reference/litreview/hydro_coupling.md`; ROADMAP item 17.

```@autodocs
Modules = [Briareus]
Pages = ["src/capsule_drive.jl"]
Private = false
```
