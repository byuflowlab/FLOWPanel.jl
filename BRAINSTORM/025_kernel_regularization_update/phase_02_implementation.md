# 025 Phase 2 — Implementation and Testing Record

**Date:** 2026-08-20. ~~Default changed to compact-support per Phase 0~~
**Superseded same day: default = GAUSSIAN** (Ryan matched-core-size ruling,
RESET BRIEF (c)); compact and Vatistas selectable. Gaussian radius_inflation
made gradient-aware (review finding 1): Δr/rc = 4.99 / 5.47 / 5.90 at tol
1e-4 / 1e-5 / 1e-6, solving e^(-z)(1+2z)=tol. Nothing committed; no HPC runs
(Phase 3).

## Implementation

- `src/FLOWPanel_elements_fmm.jl`:
  - `@enum FilamentRegularization` (`CompactRegularization` default,
    `GaussianRegularization`, `VatistasRegularization`), module Ref
    `FILAMENT_REGULARIZATION`, setter `set_filament_regularization!`
    (enum or `:compact`/`:gaussian`/`:vatistas`).
  - `_bound_vortex_velocity`: per-family denominator `D` (Phase 1);
    Vatistas branch is the original expression unchanged. Gaussian branch
    evaluated as `g/x²` with `expm1` so the `h → 0` limit is exact; compact
    branch needs no small-`h` guard (`D = B·rc² > 0` on the line).
  - `_bound_vortex_gradient`: per-family `(D, κ)` with `∇D = κ·∇A`
    (Phase 1); Vatistas `κ = A/D` reproduces the shipped `dD_coeff`
    exactly; compact `κ = 2 − rc/h` (clamped at `h ≲ eps·rc`, where
    `∇A → 0` anyway); Gaussian `κ = (1 − h·g′/2g)/g` with series limits
    `D → 2B·rc²`, `κ → 1/2` below `x² < 1e-12`.
  - `radius_inflation(::Type{VortexRing}, rc, tol)` now family-aware:
    compact `rc` (tol-independent), Gaussian gradient-aware `rc·√(2z*)`, `e^(-z*)(1+2z*)=tol` (review finding 1), Vatistas
    `rc·(2/tol)^{1/4}` (legacy). Source/doublet rules unchanged.
- `src/FLOWPanel.jl` `__init__`: env pin `FLOWPANEL_FILAMENT_REG`
  (`compact|gaussian|vatistas`) so frozen drivers (RHPC) can select a family
  with zero code changes. FMM/direct alignment is by construction: both
  backends and `radius_inflation` read the same global.

Callers audited: every filament evaluation routes through the two
`_bound_vortex_*` functions (VortexRing panels in `_induced`, panel-wake
filaments in `src/FLOWPanel_wake.jl:2854-2860`), so the family switch covers
direct kernels, FMM near-field, and wake filaments consistently. The
VortexRing scalar-potential branch (equivalent constant-doublet panel)
already used the compact family — the new default unifies velocity and
potential regularization.

## Tests (all local, julia 1.12.5, ≤4 threads)

New `test/runtests_unit_regularization.jl` (added to `runtests.jl`):
**221/221 pass** — per-family velocity vs independent infinite-filament
closed forms (sign vs singular kernel asserted), gradients vs central finite
differences at h below/at/above rc (convention `g[i,j] = ∂u_i/∂x_j`
verified against the legacy implementation), far-field agreement ≤ 1.01·tol
at each family's `radius_inflation` distance, compact exactness beyond rc
(rtol 1e-14, velocity and gradient), compact C¹ continuity at the support
boundary, on-filament/endpoint guards, `radius_inflation` values, selection
plumbing incl. error path.

Existing suites:

| suite | result |
|---|---|
| runtests_unit_fmm (compact default) | pass (after family-aware `radius_inflation formulas` fixture update) |
| runtests_unit_fmm (`FLOWPANEL_FILAMENT_REG=vatistas`) | pass — legacy behavior exactly reproduced through the env pin |
| runtests_unit_kernel_gradient | pass (Vatistas-only algebraic-convergence fixture pinned to Vatistas; added compact-exactness assertion — stronger property) |
| runtests_unit_liftingbody, runtests_unit_wake, runtests_unit_simulate | pass, unmodified |
| runtests_unit_kutta | pass (one discriminator fixture pinned to Vatistas: it *requires* a family whose kernel depends on the offset at CP distances; the family-independent correctness assertion passed under compact before the pin) |
| runtests_unit_kutta_routeb, runtests_analytical | pass, unmodified |

Fixture updates were intent-preserving: each encoded Vatistas-specific
behavior that the selectable-family change makes explicit.

**Known/pre-existing failure (not this change):** `examples/sweptwing.jl`
dies in GeometricTools loft *plotting* (`PythonCall IndexError`) during mesh
generation, before any solve — identical under the Vatistas pin. Headless
integration coverage substituted via the liftingbody/wake/simulate/kutta
suites.

## Performance regression and fix (cluster A/B, 2026-08-20)

The first selectable-family implementation read `FILAMENT_REGULARIZATION[]`
INSIDE `_bound_vortex_velocity`/`_bound_vortex_gradient` — a non-const global
Ref load plus 3-way branch per EDGE evaluation in the hottest direct loops.
Cluster evidence: the vatistas-pinned control on the new src measured the
body influence pass at **65.0 s vs 43.5–48.5 s across four old-src jobs**
(wake/solve moved <5%, so not node variance; the gaussian arm's 19.2 s
carried the same tax).

Fix — function barrier putting the family in the type domain:
`Val(FILAMENT_REGULARIZATION[])` is read ONCE per `direct!`-level call (the
three entry points whose loops reach the filament kernels:
`FastMultipole.direct!` for `AbstractBody` sources in
`src/FLOWPanel_abstractbody.jl` → `_direct_body!`, for `PanelWake` sources
→ `_direct_panelwake!`, and for `FilamentWrapper` → `_direct_filaments!` in
`src/FLOWPanel_wake.jl`), then crossed as `fam::Val{F}` through
`induced` → `_induced` → `@inline _bound_vortex_*(..., ::Val{F})`, which
compile with zero runtime family branches (`F === ...` folds). Val-less
fallback methods with the old signatures remain for cold call sites (tests,
probes): one Ref read + one dynamic dispatch per CALL, never per edge. The
`@inline` on the Val kernels is load-bearing — without it the family-expanded
bodies exceed the inlining threshold and the call overhead alone costs ~23%.
`radius_inflation`'s Ref read is per-buffer-fill-element (n_sources per call,
not n_pairs) — audited, acceptable, unchanged.

Microbenchmark (local, 1 thread, vatistas pin, 1e6 VortexRing `_induced`
velocity+gradient triangle evaluations = 3e6 edge evals, min of 5):

| variant | time | ns/panel-eval |
|---|---:|---:|
| (a) per-edge Ref read (old behavior, inlined emulation) | 0.077 s | 77.4 |
| (b) Val barrier (shipped fix) | 0.074 s | **73.7** |
| (c) hardcoded Vatistas (ideal floor) | 0.074 s | 74.2 |

(b)/(c) = 0.993; old-vs-new values agree to 2.2e-16. The emulation in (a)
is fully inlined and thus a LOWER bound on the true old cost — the cluster
regression (+34–49%) came from the same Ref-in-kernel pattern in the
non-inlined, branch-bloated real path.

## Known limitation (review finding 3b, 2026-08-20)

Under the Gaussian default the VortexRing velocity is Gaussian-regularized
while its scalar-potential branch (equivalent constant-doublet panel) remains
compact-regularized — ∇φ ≠ u inside the core region for that element (they
agree outside both cores to their contracts). Open Ryan decision, recorded in
the item's Deferred section; unification needs a Gaussian-regularized
doublet-potential derivation.

## Deliberately left open (Phase 3 / Ryan)

- Mature-state cost/error A/B via the 023 harness under the new default
  (the 023 numbers were measured under Vatistas and no longer describe the
  default; radius inflation drops 37.6·rc → rc).
- Solve-conditioning iteration-count comparison on a production-scale case.
- 018 driver compatibility + ledger entry; any commit.

## Clean mature-state benchmarks (2026-08-21, barrier-fixed src, jobs 13247381/82/83)

Mature p018_cs_f1_l3p4 state (181k particles), 64 threads, steady-step means;
thrust trace identical to the Vatistas baseline to 4 digits in every arm
(CF_x −0.07012 → −0.07044 over steps 1035–1039).

| arm | knobs | step | wake infl. | body infl. | solve |
|---|---|---:|---:|---:|---:|
| Vatistas control (barrier src) | production (4/0.4/50, 8/0.4/20) | 183.4 s | 118.7 | **44.6** | 18.1 |
| Gaussian | production knobs | 141.8 s | 109.9 | **17.4** | 12.4 |
| Gaussian | tuned (16/0.6/38, 17/0.7/109) | **72.2 s** | 42.7 | 14.4 | 13.0 |

- Control body pass 44.6 s = back inside the old-src band (43.5–48.5) ⇒ the
  Val-barrier fully removed the dispatch regression (was 65.0 s).
- Gaussian cuts the body pass 2.6x at identical knobs (radius inflation
  37.6·rc → ~5·rc) and, with its own tuned knobs, the FULL STEP lands at
  72.2 s vs the 172–177 s original production baseline: **~2.4x**. A 30-rev
  production run projects to ~21.7 h — inside a 24 h wall (previously timed
  out at rev 28.5 of 30 in 48 h).
- Achieved-error certification (direct reference, job 13247200): wake 1.8e-6
  (target 1e-4), body 1.2e-6 (target 1e-5).
- The solve (13 s) is now 18% of the step — the next lever is 021's
  ILU-GMRES; the wake phase (42.7 s) remains SFS-estimator-dominated (024's
  census + the SFS levers own that).
