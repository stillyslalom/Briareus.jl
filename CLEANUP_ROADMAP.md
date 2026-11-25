# Briareus.jl Cleanup & Refactoring Roadmap

This file is the working copy of the cleanup plan of 2026-11-25, kept on the branch so
that a session can resume from it. The status section comes first; the plan as handed
over follows unchanged except for the status marks. The lasting record of what stayed
open or was declined is item 33 of [reference/ROADMAP.md](reference/ROADMAP.md).

## Status (2026-11-25)

All six phases have been worked through on `claude/optimistic-thompson-pd2w55`. On the
head of the branch the core suite passes (7,701 tests), the GUI suite passes (2,624
tests) and the docs build completes with every tutorial.

Marks used below: **done**, **partial** (with what remains), **declined** (with the
reason in `reference/ROADMAP.md` item 33).

| Item | Status | Note |
|---|---|---|
| 1.1 `_nnls` | done | one copy, in `src/utils.jl` |
| 1.2 races and caches | done | reference models checked without mutable state |
| 1.3 moved-target areas | done | front fraction averaged over the moved areas, with a regression test |
| 1.4 `invokelatest` | done | `_ReferenceModel{L}` |
| 1.5 dead code | done | |
| 2.1 UQ through LPI | done | `sample_modes(...; through_lpi = false)`, opt-in as decided |
| 2.2 lip attenuation in imagers | done | `transmission` keyword through `render`, `imager_image`, `line_of_sight_power` |
| 2.2 BVH over lip triangles | done | `_lip_bvh`; results identical |
| 2.3 bubble frame | done | initial-wall frame kept and documented; awaits the owner's acceptance |
| 2.3 double `linked_bubble` | done | `_resolve_lpi` |
| 2.4 stale claims | done | |
| 3.1 `_StepState`, march split | done | sparse deficits included |
| 3.2 shared utilities | done | VisRad line parser kept separate (declined) |
| 3.3 protocols and types | partial | abstract types, `Keyframes`, workspace beams, BVH guard done; `Scenario`/`Evaluation` parameters declined; CBET schedules evaluated once per apply in place of type parameters |
| 3.4 envelope core | done | exact ellipse distance for keep-outs declined after measurement; keep-outs carry no weight |
| 4.1 splits and include order | partial | all five files split; 31 of 72 forward-reference pairs remain, mutual by design |
| 4.2 exports | partial | `Briareus.Benchmarks` done, `compare_visrad_kit`, spot helpers merged; 514 exports remain, a triage that is the owner's decision |
| 4.3 GUI boundary | done | `Briareus.GUIAPI`, with a test in the GUI suite |
| 5.1 `target_setup.jl` | done | `PartGeometry`; geometry reads on frames answer false or empty |
| 5.2 `shell.jl` | done | dirty domains |
| 5.3 legacy sessions | done | dropped as decided; `SESSION_VERSION` 14 |
| 5.4 first-frame replay | done | recorder not run (it needs a real window) |
| 6.1 design docs, ROADMAP | done | |
| 6.2 research notes | done | `reference/notes/` |
| 6.3 assets and `dev/` | done | |

Open decisions for the owner: whether geometry reads on a frame or assembly `PartSpec`
should throw; whether to triage the export surface; whether the 0.4.0 version bump
stands; whether the initial-wall bubble frame is accepted.

---

## Ground Rules for Execution

- **One commit (or PR) per numbered sub-task** where practical; never mix behavior/physics changes with pure mechanical refactors.
- **Guardrails to run after every phase:**
  - `julia --project=. test/docrefs.jl` (required whenever files are split/renamed or exports change; update `Pages` in `docs/src/reference/*.md`).
  - `julia --project=. -t auto -e 'using Pkg; Pkg.test()'` (or `BRIAREUS_TEST_SHARD=1|2|3` during iteration; `test/inference.jl` and `test/aqua.jl` must stay green).
  - `julia --project=BriareusGUI -t auto -e 'using Pkg; Pkg.test(; julia_args = ["--heap-size-hint=2G"])'` for any core API or GUI change.
  - `julia dev/prose/lint.jl --staged` before every commit.
  - Use `BRIAREUS_NO_WORKLOAD=true` while iterating, then run a clean build without it at the end of Phases 3–5 to confirm `@compile_workload` in `src/Briareus.jl` and `BriareusGUI/src/BriareusGUI.jl`.
- **Per `CLAUDE.md`:** Briareus is `0.x`—do not leave parallel old options or compatibility shims behind when replacing a worse design.

---

## Phase 1: Correctness Bugs, Thread-Safety Races & Dead Code (Zero API Breakage)

High-value fixes that do not change public signatures.

### 1.1 Remove silently shadowed duplicate `_nnls`
- **Problem:** `src/constraints.jl:466` defines `_nnls(A::Matrix{Float64}, b::Vector{Float64})`, which silently shadows `src/dante_response.jl:419-447` (`_nnls(A::AbstractMatrix, b::AbstractVector; tol, maxiter)`) at its only call site `src/dante_response.jl:405` (which passes concrete `Matrix{Float64}, Vector{Float64}` with no kwargs).
- **Action:**
  - Delete `_nnls` from `src/dante_response.jl:419-447` (or keep one canonical `_nnls` in a shared utility file, see 3.2).
  - Verify `test/dante_response.jl` and `test/constraints.jl` remain identical (since `dante_response.jl:405` was already dispatching to `constraints.jl:466`).

### 1.2 Fix data races and unsafe global caches
- **`_CeilingValue.checked` & `_ReferenceModel.checked` (`src/constraints.jl:858, 878-895, 902`):** Mutates an unlocked `Ref{Any}` inside `_check_reference_model` during concurrent ForwardDiff chunk threads (`_pointing_gradient!`, `_pointing_jacobian`).
  - **Fix:** Guard with a lightweight lock, use `@atomic` / `Threads.Atomic` semantics, or better yet, validate reference model compatibility once at `optimize_pointing` / `check_constraints` entry instead of mutating state on every evaluation.
- **`Keyframes.between` (`src/dynamics.jl:327-341`):** `model_at` reads and writes `kf.between[]` outside `kf.lock` (unlike `deposition_at` at `src/dynamics.jl:356`).
  - **Fix:** Wrap the `kf.between[]` lookup and store in `lock(kf.lock) do ... end` (avoiding holding the lock across expensive work if desired, or matching `deposition_at`).
- **`_NORA_TABLE1` (`src/implosion.jl:238-248`):** Add a `ReentrantLock` matching `_HENKE_LOCK` / `_RUBERY_LOCK` / `_KSHELL_LOCK`.
- **`_DIAGNOSTIC_KINDS` (`src/diagnostics.jl:64, 77`):** Guard reads/writes with a `ReentrantLock`.
- **`_VR_NTABS` (`src/beams.jl:556-572`):** Replace `Ref(Dict...)` unlocked read at `L561` with an atomic field wrapper (`mutable struct _AtomicBox{T}; @atomic x::T; end`) or `lock(_VR_NTABS_LOCK)` so the copy-on-write swap has proper acquire/release semantics.
- **`_CHANNEL_CACHE` (`src/diagnostics.jl:856-868`):** Keying solely on `hash(cs)::UInt` risks silent hash collisions and grows unboundedly. Key on the actual normalized channel spec (or `(hash, spec)` with equality check) and cap cache size.

### 1.3 Fix moving-target stale area bug in `_drive_history`
- **Problem:** `src/history.jl:673-678` computes `front[c.name]` area averages using initial `t.elem_area` instead of `sol.model.target.elem_area` (which is bound as `tg` four lines later at `L682`). Also audit `L571` (`msplit` / `_mband_from_split`) and `L669` (`transmitters`) for any other uses of initial `t` when `model_at(τ)` returns a moved `Target`.
- **Action:** Use `tg = sol.model.target` consistently for all per-step area reductions and add a regression test with `moving_hohlraum` + a thin-foil window.

### 1.4 Eliminate the `Base.invokelatest` JET workaround
- **Problem:** `src/picket_fluence.jl:666-672` wraps `_check_model` in `Base.invokelatest` solely so `JET.report_opt` ignores dynamic dispatch on `_ReferenceModel.lpi::Any`.
- **Action:** Remove `Base.invokelatest` by making `_ReferenceModel` concrete or performing reference-model validation once when binding constraints to a scenario (eliminating both the `Ref{Any}` race from 1.2 and the per-evaluation world-age call). Confirm `test/inference.jl` stays JET-clean.

### 1.5 Remove verified dead code
- **Core `src/`:** Remove uncalled private helpers `_gate_lmax` (`src/constraints.jl:1210`), `_one_point` (`src/viewfactors.jl:458`), and `_vr_col` (`src/visrad.jl:1056`).
- **GUI `BriareusGUI/src/`:** Remove `_g` (`BriareusGUI/src/qt/shell.jl:674`), `_is_beam_node` (`BriareusGUI/src/controllers/view.jl:377`), and dead `@qmlfunction` `br_remove_window` (`BriareusGUI/src/qt/shell.jl:1884, 3072`). Remove `br_node_flags` and `br_node_quads` from `@qmlfunction` (`BriareusGUI/src/qt/shell.jl:3066`) since QML never invokes them directly (keep them as plain internal helpers for `br_node_menu`).

---

## Phase 2: Physics Consistency Fixes (Confirm Intent Before Changing Results)

Each item here alters physical outputs or closes a modeling inconsistency. Per `CLAUDE.md`, measure and record any numerical delta in the commit message and `reference/research-notes.md`.

### 2.1 UQ error propagation vs LPI (`src/uncertainty.jl`)
- **Problem:** Both `mode_sensitivity` (`src/uncertainty.jl:32`) and its Monte Carlo check `sample_modes` (`src/uncertainty.jl:289`) call `apply_lpi` *before* perturbing beam pointing/power, dropping nonlinear error propagation through CBET and bubble iris clipping.
- **Action:**
  - Decide explicitly: either (a) add `through_lpi::Bool = false` (or make `true` the default) so `sample_modes` perturbs setup beams *before* `apply_lpi` and `mode_sensitivity` differentiates through `apply_lpi` when present, or (b) keep post-LPI linear sensitivity for speed and have `sample_modes` expose pre-LPI vs post-LPI sampling so users can quantify the LPI-amplified variance.
  - Update docstrings in `src/uncertainty.jl` to state clearly whether errors are pre-LPI (laser delivery) or post-LPI (effective beams).

### 2.2 Closing LEH lip attenuation in `XrayImager` and `line_of_sight_power`
- **Problem:** `DanteView` (`src/dante.jl:254-276`) ray-tests transmitting `LipClosure` apertures and weights radiance by `sol.transmission`, whereas `XrayImager` (`src/imager.jl:123` via `render`) and `line_of_sight_power` (`src/imaging.jl:80-114`) treat all apertures as fully transparent.
- **Action:** Propagate `sol.transmission` for closing-lip apertures through `render` / `imager_image` / `line_of_sight_power` (and optionally thread `_dante_footprint` at `src/dante.jl:259-277` using a small BVH over lip triangles instead of a linear scan). Add a test comparing an SXI `XrayImager` view of a closing LEH against `DanteView`.

### 2.3 `BubbleAbsorption` + `WallMotion` consistency & `PicketLinkedBubble` double deposition
- **Problem:**
  1. `BubbleAbsorption` (`src/bubble.jl:343`) shifts `z_b` with `WallMotion`, but `bubble_depth`, `ρ = m.radius - x` (`L413`), and ridge radius `m.ridge * m.radius` (`L497`) remain anchored to initial `m.radius`. Confirm whether this is intended (initial-wall reference frame) and either make it consistent or document the exact frame explicitly.
  2. `PicketLinkedBubble` (`src/picket_fluence.jl:835-846`) runs `linked_bubble(m, beams)` twice per evaluation (once in `apply_lpi`, once in `_lpi_deposit`).
- **Action:** Pass the linked bubble state from `apply_lpi` into `_lpi_deposit` (or bundle `apply_lpi` + deposition in one helper on `Scenario` / `drive_history`) to cut `PicketLinkedBubble` picket deposition cost in half.

### 2.4 Fix stale physics claims in docs
- Update `README.md:145` (*"Briareus does not trace refraction"*) to reflect `CoronaOperator` (`src/corona.jl`, spherical refraction) vs non-refracting hohlraum `DepositionOperator`.
- Update `ResponseKernel` docstring (`src/implosion.jl:94-95`) which contradicts `young2020_kernel()` in the same file (`L580`).

---

## Phase 3: Internal Deduplication, Protocols & Type Tightening

Refactor internals without fragmenting the public API yet.

### 3.1 Replace the untyped `burnthrough::Vector{Any}` channel
- **Problem:** `_drive_history` (`src/history.jl:466-476`) stuffs `Float64` burnthrough times, `_WallHistory`, `_FoilHistory`, and `LipClosure` closed fractions into one `Vector{Any}` passed via `solve(...; burnthrough = ...)`.
- **Action:**
  - Introduce a concrete internal `_StepState` struct (e.g., `burnthrough::Vector{Float64}`, `wall_history::Union{Nothing,Vector{_WallHistory}}`, `lip_closed::Union{Nothing,Vector{Float64}}`) passed to `solve` (keeping the public `burnthrough` keyword for standalone `solve` calls cleanly typed).
  - While touching `src/history.jl`, check whether `deficits` (`L610, L835-844`) can store sparse/only-`lossy` element subsets (avoiding storing `N`-length vectors for non-wall elements like capsules/apertures) and extract sub-steps out of the 408-line `_drive_history` into small named helpers (`_init_history_state`, `_solve_history_step!`, `_record_history_step!`, `_assemble_history_output`).

### 3.2 Consolidate duplicated utilities (`src/utils.jl` or early base file)
Create a small early-included `src/utils.jl` (included right after `src/geometry.jl`) and consolidate:
- **Cubic Hermite interpolation:** Unify the 3 copies in `src/target.jl:884`, `src/corona.jl:119`, and `dev/corona/Corona.jl:79` (keep uniform/non-uniform methods cleanly typed).
- **Binary min-heap:** Unify `_heap_push!`/`_heap_pop!` duplicated between `src/meshops.jl:871-890` and `src/clearance.jl:832-855`.
- **Lawson–Hanson NNLS:** Move the single canonical `_nnls` (from `src/constraints.jl:466`) to `src/utils.jl` so both `constraints.jl` and `dante_response.jl` call it cleanly.
- **Perpendicular basis helper:** Replace ad-hoc perpendicular vector constructions (`src/geometry.jl:229, 246`, `src/target.jl:1693`, `src/cadfit.jl:179`, `src/uncertainty.jl:213`, `src/viewfactors.jl:331`) with one `_ortho_basis(n::Vec3)`.
- **Material density tables:** Consolidate `_COMPOUNDS` (`src/xray_optics.jl:27-40`), `_XRD_BULK_DENSITY` (`src/dante_response.jl:48`), and `_FILTER_DENSITY` (`src/backlighter_signal.jl:346-351`) into one single source of truth in `src/xray_optics.jl` (this also fixes the forward-reference `xray_optics.jl:115` → `backlighter_signal.jl:346`).
- **CSV parsing:** Reuse `src/io.jl`'s CSV/table parsing inside `src/visrad.jl:30-65` (where compatible) and replace `_nan_table` (`src/cbet.jl:641`) + the 7 bespoke `read_csv` loops in `src/cbet.jl`, `src/backlighter_signal.jl`, and `src/dante_response.jl` with options on `_data_table`.

### 3.3 Formalize abstract protocols & tighten hot struct types
- **Add `abstract type AbstractParameterization end`** (`GroupPointing`, `FramePointing`, `HeldParameters`, `PowerScales`, `PhasedPointing`) and **`abstract type AbstractKeepOut end`** (`ApertureKeepOut`, `SphereKeepOut`, `HardwareKeepOut`, `CaptureKeepIn`, `IncidenceRule`, `Omega1KeepOut`). Replace `hasmethod(parameter_names, ...)` duck-typing (`src/optimize.jl:2053`, `src/parameterizations.jl:51,192`, `src/constraints.jl:986`) with clean dispatch.
- **Tighten abstract fields on hot types** (benchmarking compile time against `test/inference.jl` to avoid specialization bloat):
  - `Scenario{P,M,L}` (`src/optimize.jl:648`) and `Evaluation{T,B,E,S<:Scenario}` (`src/optimize.jl:927`).
  - Parameterize `CBET` (`src/cbet.jl:102-116`) or wrap time-varying scalar schedules in a small concrete union/callable type instead of 8 `::Any` fields.
  - Parameterize `Keyframes` (`src/dynamics.jl:232`, concrete `Vector{Beam{Float64}}`) and `VisRadWorkspace.beams` (`src/visrad_workspace.jl:64`).
  - Add a bounds check or stack-growth guard on `MVector{64,Int32}` traversal stacks in `src/bvh.jl:119, 158`.

### 3.4 Unify the beam-envelope & clearance geometry core
- **Problem:** `clearance.jl`, `keepout.jl`, `talis.jl`, and `design_rules.jl` contain three near-identical envelope structs (`BeamEnvelope`, `_RimFrame`, `_EnvelopeFrame`), mixed ellipse math (first-order Sampson in `keepout.jl` vs exact Eberly in the other three), duplicated 1ω setup (`talis.jl:293-310` vs `design_rules.jl:884-889`), and duplicated p-norm ceiling machinery (`constraints.jl:827` vs `picket_fluence.jl:601`).
- **Action:**
  - Extract a single shared `_BeamConeFrame{T}` (generic in `T` for ForwardDiff) and shared 1ω geometry helper (`_omega1_params`) used by all four files.
  - Evaluate upgrading `ApertureKeepOut` / `SphereKeepOut` from Sampson to exact Eberly distance (using `clearance.jl:_ellipse_closest` + implicit-function/Danskin partials as `HardwareKeepOut` already does), removing the ellipse-distance discrepancy between `keepout.jl` and `clearance.jl`/`talis.jl`.
  - Share the smooth p-norm + `_exact_value` constraint builder between `irradiance_ceiling` and `picket_ceiling`.
  - Clarify/enforce `clearance_floor` vs `keepout.weight` (e.g., construct keep-outs without penalty `weight` when used as constraints, or separate the geometric keep-out object from the quadratic penalty weight).

---

## Phase 4: Core Module Decomposition, Include-Order Cleanup & Export Hygiene (0.x Breaking)

### 4.1 Untangle the 72 forward-reference file pairs & split god files
Reorder `include`s in `src/Briareus.jl` into strict DAG layers (zero backward/cyclic file dependencies) and split the largest multi-responsibility files:

| Current File | Split Into |
|---|---|
| `src/optimize.jl` (2,441 LOC) | `src/pointing.jl` (`GroupPointing`), `src/scenario.jl` (`Scenario`, `Evaluation`, `_mul`), `src/objectives.jl` (all objective factories), `src/optimizer.jl` (`optimize_pointing`, `_lbfgsb`, parallel ForwardDiff, `PointingResult`, `report`) |
| `src/target.jl` (1,773 LOC) | `src/albedo.jl` (all `AlbedoModel`s, `HeatFront`, `RadiusSchedule`, `Conversion`), `src/materials.jl` (`Material`, presets), `src/target.jl` (`Component`, `Target`, quadrature) |
| `src/constraints.jl` (1,381 LOC) | `src/al_solver.jl` (ALGENCAN loop, Phase 1, relaxations, KKT NNLS), `src/constraints.jl` (`Constraint` + built-in constraint factories) |
| `src/visrad.jl` (1,505 LOC) | `src/visrad_io.jl` (phase plates, SST/summary/power/pulse I/O), `src/visrad_kit.jl` (kit writer, Legendre/Ylm fits, comparisons) |
| `src/hmatrix.jl` (2,095 LOC) | `src/hmatrix_tree.jl` (cluster tree + shaft culling + visibility certificates), `src/hmatrix.jl` (`HMatrix`, ACA/SVD assembly, matvec/matmat) |

Break the specific bidirectional include cycles:
- Move `_gauss_legendre` out of `src/corona.jl` into `src/utils.jl` (removes `beams.jl` → `corona.jl`).
- Move `_stagnate` out of `src/bubble.jl` into `src/utils.jl` (removes `dynamics.jl` → `bubble.jl`).
- Move `io.jl` earlier (right after `geometry.jl` / `beams.jl`) so `dynamics.jl`, `implosion.jl`, `xray_optics.jl`, `cbet.jl` never forward-reference `_data_table` / `read_csv`.
- Include `unconverted.jl` (or split `unconverted_types.jl`) before `history.jl`, `retune.jl`, and `talis.jl` so `_UnconvertedLight` and `nif_unconverted` are defined before use.
- Move `_vrw_is_json` dispatch cleanly so `visrad_workspace.jl` and `visrad_json.jl` have a one-way dependency.
- Move `dante1()`/`dante2()` built-in wrappers to `src/diagnostics.jl` (or include the registry before `dante.jl`) to break `dante.jl` ↔ `diagnostics.jl`.

Update `docs/src/reference/*.md` `@autodocs` `Pages` lists and verify with `test/docrefs.jl`.

### 4.2 Trim the 537-symbol flat export surface
- **Create `module Benchmarks` (or `module PublishedData`) inside `Briareus`** (`Briareus.Benchmarks`, not re-exported at top level) and move all **27 digitized paper-figure / shot-fixture functions** out of core exports and physics modules (`src/cbet.jl:623-904`, `src/dynamics.jl:66-155`, `src/bubble.jl:259-271`, `src/implosion.jl:561-618`, `src/lpi.jl:206`, `src/picket_fluence.jl:191-219`, `src/facilities.jl:376`):
  - `n221204_fig9`, `n220919_fig9`, `n210808_fig1`, `pre109_fig1`, `hybride_shots`, `n221204_fig10`, `n221204_pulse`, `zylstra2022_fig2`, `n221204_fig9_beams`, `n221204_cbet`, `hybride_cbet`, `hybride_cbet_iraum`, `hybride_mock_pointing`
  - `young2020_radii`, `young2020_albedos`, `young2020_green`, `young2020_kernel`, `young2020_pulses`, `young2020_dropped_quads`
  - `izumi2021_bubble`, `izumi2021_absorption`, `basecamp_picket_fluence`
  - `chen2020_leh`, `pre109_fig1a_wall`, `kritcher2024_capsule`, `kritcher2024_radius`, `kemp2026_sbs`
  - Update `examples/`, `docs/literate/`, and `test/` to use `Briareus.Benchmarks.<name>` (or `using Briareus.Benchmarks`).
- **Disambiguate overlapping public APIs:**
  - Rename one of the two incompatible `compare_visrad` methods (`src/visrad.jl:1457` kit-directory comparison vs `src/visrad_results.jl:154` workspace-vs-solution comparison), e.g. `compare_visrad_kit` vs `compare_visrad`.
  - Consolidate spot helpers (`nif_spot`, `visrad_spot`, `visrad_nif_spot`) and `VISRAD_PHASE_PLATES` / `VISRAD_NIF_PLATES`.
  - Consider submodules for `VisRad` (`Briareus.VisRad`) and `Diagnostics` (`Briareus.Diagnostics`) if aiming for a tighter `<250`-export core API before `1.0`.

### 4.3 Formalize the `Briareus` ↔ `BriareusGUI` boundary
- **Problem:** `BriareusGUI/src/` directly calls **23 private `Briareus._*` internals** across 34 call sites:
  - Frames/placement: `_relative`, `_step_axis`, `_step_refs`, `_axis_rotation`
  - Targets/beams/pulses: `_profile_radius`, `_beam_key`, `_beam_pulse`, `_is_nif`, `_nif_beam_table`, `_omega1_beamlets`, `_part_series`
  - Drive/unconverted: `_radiative`, `_observed_flux`, `_UnconvertedLight`, `_unconverted_at`
  - Diagnostics/IO: `_BUILTIN_DIAGNOSTICS`, `_read_all`, `_checked`, `_resolve`, `_definitions`, `_imager_pixels`, `_saved_paths`, `_map_csv`
- **Action:** Either promote genuinely useful helpers to documented public API (e.g., `step_references`, `step_axis`, `profile_radius`, `beam_key`, `has_radiation`) or expose a narrow, tested `Briareus.Internal` / `Briareus.GUIAPI` surface and eliminate ad-hoc `Briareus._*` reaches in `BriareusGUI/src/`. Add a test in `BriareusGUI/test/` asserting zero unexpected `Briareus._` references.

---

## Phase 5: `BriareusGUI` De-Crufting

### 5.1 Split `BriareusGUI/src/controllers/target_setup.jl` (3,434 LOC) & bound global caches
- Split into cohesive files under `BriareusGUI/src/controllers/`:
  - `mesh_spec.jl` (`MeshSpec`, `MESH_PRESETS`, `LEH_CLOSURE_FIELDS`)
  - `windows.jl` (`WindowSpec`, window CRUD)
  - `parts.jl` (`PartSpec`, part CRUD, hohlraum sections/profiles)
  - `part_frames.jl` (user frames, placement steps, nudge/turn, `_euler_angles`—consider moving `euler_angles(R)` to core `src/frames.jl`)
  - `roles.jl` (part/piece roles, `observed_components`, `wall_components`, `foil_components`)
  - `cad_cache.jl` (STL/VRW/target caches + background meshing)
  - `starters.jl` (`STARTERS`, `insert_starter!`)
  - `target_setup.jl` (`TargetSetup`, `build_target`, `target_summary`)
- **Bound global caches:** Convert `_VRW_PARTS`, `_VRW_CACHE`, `_STL_CACHE`, and `_STL_SOLIDS` (`BriareusGUI/src/controllers/target_setup.jl:777, 877, 897, 997`) from unbounded `Dict{Any,...}` to small typed LRU caches (matching `_TARGETS`) and expose a `clear_caches!()` helper for tests/memory release.
- **Deduplicate `fit_sphere` (`target_setup.jl:3091-3127`):** Reuse core `src/cadfit.jl:134` `_fit_sphere`'s 4×4 normal-equation solver inside the IRLS loop instead of allocating dense `N×4` matrices via `reduce(vcat, ...)`.
- **Clean up `PartSpec`:** Separate geometry fields from `:frame` and `:assembly` (either distinct spec structs under an abstract `SceneNodeSpec`, or validated sub-structs so frames/assemblies don't carry 14 dead geometry/window/hohlraum fields).

### 5.2 Split `BriareusGUI/src/qt/shell.jl` (3,203 LOC) & tighten QML bridge
- Split into focused files under `BriareusGUI/src/qt/`:
  - `platform.jl` (Qt/GLFW guards, Win32 `ccall`s, `_warm_glyph_atlas!`, `_warm_canvases!`)
  - `models.jl` (the 18 mutable row structs + `JuliaItemModel` diff helpers)
  - `refresh.jl` (`BriareusShell`, `_refresh!`, and panel-specific `_refresh_*!` functions; break the 185-line `_refresh_selection!` into per-kind helpers: `_refresh_part_props!`, `_refresh_beam_props!`, `_refresh_window_props!`, `_refresh_analysis_props!`)
  - `callbacks.jl` (`_with_shell`, `br_tick`, all `br_*` callbacks, `_register_qml_functions`)
  - `window.jl` (`briareus(...)` entry point)
- **Granular dirty flags:** Replace the single `sh.dirty[]` bit (which unconditionally runs all 11 `_refresh_*!` functions and rebuilds the 171-line `model_tree` on every edit) with a small bitmask / set of dirty domains (`:tree`, `:selection`, `:beams`, `:drive`, `:time`, `:optimize`, `:review`, `:talis`, `:sweep`, `:view`) so, e.g., scrubbing a time cursor or toggling a view option doesn't rebuild unrelated panels.

### 5.3 Drop pre-1.0 legacy session migration baggage
- **Problem:** `BriareusGUI/src/controllers/workspace.jl:8-16` maintains conversion paths for **12 legacy versions (`v1`–`v12`)** (`SESSION_VERSION = 13`) plus 75 KB of `v2`/`v7` test fixtures (`BriareusGUI/test/data/session_v2.briar.toml`, `session_v7_*.toml`), contradicting `CLAUDE.md`'s 0.x no-compat-shims policy.
- **Action (confirm with user whether any personal pre-v13 `.briar` files still need loading):**
  - If not needed (or after batch-resaving personal files): reset to a single clean schema version, remove `_TEMPLATE_DEFAULTS`, `_template_parts`, pre-v8 rename maps, and pre-v3/v5/v6/v7 branches from `workspace.jl`, and replace `BriareusGUI/test/data/session_v{2,7}_*.toml` with compact current-schema fixtures.
  - Also split `BriareusGUI/src/controllers/session.jl` (1,889 LOC) and `timeline.jl` (1,606 LOC) along their existing comment sections (`session_drive.jl`, `session_retune.jl`, `timeline_diagnostics.jl`).

### 5.4 Replace the fragile `first_frame.jl` gensym precompile replay
- **Problem:** `BriareusGUI/src/qt/first_frame.jl` (259 KB, 487 lines) replays 234 compiler-generated closure signatures (19 anonymous numbered closures like `BriareusGUI.var"#95#96"`, `Makie.var"#1916#1917"`) via `Core.eval` with a silent `try/catch` (`BriareusGUI/src/BriareusGUI.jl:189-206`), and `dev/gui_first_frame.jl` appends without pruning dead signatures.
- **Action:**
  - Make `dev/gui_first_frame.jl` **overwrite** instead of append-merge, and filter out anonymous numbered closures (`var"#\d+#\d+"`) and unstable `BriareusGUI` internal closures.
  - Move any `BriareusGUI`-owned methods that still appear into explicit calls inside `@compile_workload` in `BriareusGUI/src/BriareusGUI.jl`.
  - Log a debug/dev-mode count of failed `_replay_precompiles` lines so drift is visible instead of silently swallowed.

---

## Phase 6: Reference Docs, `dev/` & Repo Hygiene

### 6.1 Prune stale design docs & refresh `ROADMAP.md`
- Per the rule at the top of `reference/assemblies-design.md`, `reference/frames-design.md`, and `reference/diagnostics-design.md`: delete all already-landed phase write-ups (keep only open Phase E in `assemblies-design.md`, Phase 4 in `frames-design.md`, and open Phase 5 items in `diagnostics-design.md`, or archive completed design specs cleanly).
- Update `reference/ROADMAP.md`:
  - Fix stale claims (e.g., `L762-763` claims GUI screenshots predate the Diagnostics tab and Help window, though `docs/src/assets/gui/diagnostics.png` and `help.png` already exist).
  - Renumber or group items cleanly (removing gaps `1, 3, 5` and out-of-order `37`) and replace the 4-line "Code health" section with any deferred items from this roadmap.

### 6.2 Split `reference/research-notes.md` (15,071 lines / 1.1 MB)
- Split the 121 `##` sections into thematic or chronological files under `reference/notes/` (e.g., `viewfactors_hmatrix.md`, `deposition_visrad.md`, `hohlraum_physics.md`, `optimization_repointing.md`, `diagnostics_spectral.md`, `gui_cad.md`) with a short `reference/research-notes.md` index table of contents so editors, git diffs, and LLM context lookups stay fast.

### 6.3 Minor asset & `dev/` deduplication
- Remove duplicate `docs/src/assets/logo.svg` (identical to `BriareusGUI/src/qml/icons/briareus.svg`) or reference one canonical copy in `docs/make.jl`.
- Audit `dev/corona/Corona.jl` vs `test/corona_raytrace.jl`: either have `test/corona_raytrace.jl` include/call the minimal spherical stepper from `dev/corona/` or clearly isolate the minimal test reference without duplicating `_hermite` and ray-stepping math across three places.
- Clean up the unsourced `_LENS_DISTANCE = 7700.0` constant in `BriareusGUI/src/controllers/beam_geometry.jl:124` (source it from `src/facilities.jl` or document its origin).

---

## Suggested Execution Order & Checkpoints

1. **Checkpoint A (Phase 1):** ~1–2 PRs. Zero API breaks; fixes real races (`constraints.jl`, `dynamics.jl`, `implosion.jl`, `diagnostics.jl`, `beams.jl`), shadowed `_nnls`, moving-target area bug, `invokelatest` hack, and dead code.
2. **Checkpoint B (Phase 2):** User sign-off on 2.1 (`uncertainty.jl` LPI error model) and 2.3 (`BubbleAbsorption` wall frame), then land 2.1–2.4 with measured deltas.
3. **Checkpoint C (Phase 3):** Pure internal cleanup (`_StepState`, `src/utils.jl`, `AbstractParameterization`/`AbstractKeepOut`, type tightening, unified beam-envelope core). Verify `test/inference.jl` and compile/eval benchmarks stay equal or faster.
4. **Checkpoint D (Phase 4):** Coordinated `0.4.0` core cleanup (god-file splits, DAG include order, `Briareus.Benchmarks` submodule, GUI API surface).
5. **Checkpoint E (Phases 5 & 6):** `BriareusGUI` file splits, bounded caches, schema & `first_frame.jl` cleanup, and `reference/` hygiene.
