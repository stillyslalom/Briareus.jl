# CLAUDE.md

Briareus.jl — view-factor radiation and laser-pointing design (a VisRad replacement/
complement) for indirect-drive hohlraums and backlighter foils. [README.md](README.md):
scope and facility-data provenance. [ROADMAP.md](ROADMAP.md): backlog and what was done.
`reference/research-notes.md`: every physics default and measured constant, with its source
and confidence; `reference/litreview/`: literature reviews by topic. Read the relevant
notes section before changing a default or a calibrated number.

## Commands

```bash
julia --project=. -t auto -e 'using Pkg; Pkg.test()'    # full suite (~15 min on 4 threads)
BRIAREUS_TEST_SHARD=2 julia --project=. -t auto -e 'using Pkg; Pkg.test()'  # a third (1|2|3)
julia --project=. test/docrefs.jl                        # docs guard, ~3 s
julia --project=docs -t auto docs/make.jl                # docs; runs every tutorial
julia --project=docs -t auto docs/literate/<name>.jl     # one tutorial as a script
julia --project=BriareusGUI -e 'using Pkg; Pkg.test()'   # GUI subpackage
```

`BRIAREUS_NO_WORKLOAD=true` skips the precompile workload (end of `src/Briareus.jl`) for
quick rebuilds; caches don't key on it, so touch a source file before timing first results.
The shard map is at the top of `test/runtests.jl`: a new test file goes in it, and files
whose helpers others use share a shard. Expensive fixtures (deposition operators, view
factors) dominate the suite: build one per geometry and share it across testsets.

## Working rules

- Validate against independent references (closed forms, exact convolutions, converged
  quadrature, Monte Carlo, ray traces converged at nspot ≥ 121), never against the method
  under test. Published shot data used as a check (e.g. PRE 109 Fig. 10, Dante) is never
  also a fit target; say which data a calibration was fitted to.
- New paths must leave existing ones bit for bit unless the change is the point: say so in
  the commit and measure the difference. Screens that skip work (zero weights, culled
  rays) must be exact on Float64.
- Tests: tolerances no tighter than ~1e-12 relative on summed quantities (CI CPUs round
  differently); no `@test` inside loops over elements or points — check the collection
  with one `@test all(...)` (a per-element loop once tripled the suite's test count).
- New public names need a docstring rendered by some `@docs`/`@autodocs` block (new src
  files go in a reference page's `Pages`); `test/docrefs.jl` enforces this and that every
  `@ref` resolves. Docs follow CompactLES.jl: Literate tutorials in `docs/literate/`,
  registered in `docs/make.jl`, each ending in "What this checks" with exact references.
  The docs environment has no ForwardDiff/Optim, and `Vec3` clashes with Makie there.
- Provenance: docstrings say when a value is derived or an engineering default rather than
  a facility criterion; numbers go in the research notes with source and confidence.
  Papers are cite-only (numbers with citations, never copied text or figures).
  Non-public data (anything under `~/Dev/thor`: shot pointing exports, THOR CAD, VisRad
  workspaces) never enters the repository; public examples use mock-ups.
- Windows: Python `open()` defaults to cp1252 and writes CRLF in text mode — pass
  `encoding='utf-8', newline='\n'`, or use the Edit tool. Julia triple-quoted strings
  dedent, so multi-line `replace` patterns built from them miss indented code.

## Conventions

- Units: mm, ns, TW, TW/mm² (×1e14 → W/cm²), eV; angles in degrees at the API.
- Beams: direction fixed by the port (`k = -port`); pointing translates the axis through
  `aim`; best focus at `aim + defocus·k`; the spot `a` axis at `angle` from port θ̂ toward
  φ̂. A `DepositionOperator` depends on geometry and ports only (not materials, aims,
  spots), is valid only while ports are unchanged, and serves every aim — so don't drop
  its points on a profile cut-off at construction.
- Elements are radiating patch sides; one-sided surfaces radiate from the front and their
  back side absorbs laser light as a per-component loss (`<name>_outside`). Apertures are
  transparent with α = η = 0. Patches must be planar (on a cylinder, mesh anything but
  chord rectangles as triangles); remeshed CAD elements are the one folded exception.
- Optimization path (`apply`, `deposit`, the flux response, moments, objectives): generic
  in the number type for ForwardDiff. Threaded loops go through a function barrier (a
  runtime eltype captured in a `Threads.@threads` closure boxes and allocates per point).
  Objectives must be thread-safe (gradient chunks run concurrently). Keep heterogeneous
  tuples out of generators and abstract-typed loops on hot paths (inference blow-ups).
- A long function that takes user callables (closures, objectives, callbacks, schedules
  like `beams_at`) is `Base.@nospecializeinfer` with `@nospecialize` on those arguments;
  `@nospecialize` alone still infers per closure type through the keyword wrapper, and
  each closure then recompiles the whole body. `test/inference.jl` guards this (compiled
  instances of the keyword body) and checks the optimization path on Float64 and the
  optimizer's duals with `@inferred` and JET (`report_opt`, `report_call`, `target_modules =
  (Briareus,)`): add a case there for each such function and each new hot-path function. A
  Float64 constant (`Inf`, `init = 0.0`) returned beside duals is the usual culprit.
- Mode analysis integrates Pₗ/Yₗₘ over each facet's exact solid angle; centroid × A cosθ/R²
  leaks P0 into P2 at 0.5 %.
- Deposition sums the best-focus profile over the f-cone directions; the moment-matched
  `intensity` is for envelopes only (wrong shape far from focus).
- Large radiation models (> 10k elements: `Float32Matrix`; > 15k: `HMatrix`) are never
  densified or factored: go through `mul!`/`*`. `RadiationModel{M}` is parametric, so
  loop over `m.F` only behind a function barrier.
- `Component.part` names the user's part (CAD fits split parts into pieces); name-taking
  APIs go through `parts`/`part_pieces`, and rebuilding a `Component` must pass `c.part`.
- Calibrated tables and constants (VisRad beam profiles, wall loss and re-emission, CBET,
  dynamics) are fitted to specific data recorded in the research notes; don't refit them
  to whatever run is at hand — remeasure against the data they came from.

## Layout

`src/Briareus.jl` includes the files in dependency order (the list there is the map):
roughly geometry, meshes and BVH → targets and CAD fitting → beams, facilities,
deposition → view factors (`viewfactors.jl`, `storage.jl`, `hmatrix.jl`) → radiation →
analysis, imaging, optimization → time dependence (`history.jl`, `spectral.jl`, Dante,
`dynamics.jl`) → clearance, backlighters, uncertainty, implosion → LPI, VisRad and CSV
I/O, CBET. `examples/` holds the facility studies (`n221204.jl` is the shot benchmark);
`BriareusGUI/` is the GUI subpackage.
