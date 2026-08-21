# 024: Far-Field Approximations Inside the Direct Kernels

**Opened:** 2026-08-20 (Ryan directive, staged by agent; spun out of 023)
**Current phase:** Phase 0 (eligibility census) — GATE before any implementation
**Item-level approvals:**
- [ ] Technical Completion
- [ ] Clear-Context Approval
- [ ] User Approval After AI Discussion

## Current status and next actions

> **RESET BRIEF (a) — 2026-08-20 (staging).** Ryan's directive: in the direct
> (near-field) kernel functions, when an evaluation point is far enough from
> the source, evaluate the **singular 1/r kernel for vortex particles** (skip
> the erf/exp regularization math) and the **point source/doublet
> approximation for panels** (skip the full panel-integral). Rationale from
> 023: the near-field direct lists carry most of the influence cost, and the
> profiler shows the transcendental budget (exp/erf/log/atan in
> `zeta_gauserf`/`g_dgdr_gauserf`, `compute_source_dipole`, `_induced`) is a
> large share of the direct kernels. **HARD GATE (Ryan): before implementing,
> measure how many direct-list interactions are far enough for the
> approximation to satisfy the error tolerance.** Census staged as Phase 0.

In-flight jobs: none yet.

## Objective and scope

Reduce near-field direct-evaluation cost in both influence passes by
switching, per pair, to the asymptotic kernel when separation guarantees the
approximation error is below the pass's tolerance (wake 1e-4 relative, body
1e-5 relative — the 023 certification targets). Out of scope until Phase 0
reports and Ryan rules: any kernel code change; changes to MAC/tree logic
(023 owns knob tuning); SFS estimator redesign (noted as a sibling lever
below, Ryan's call).

## Why this should pay (023 evidence)

- 023 nf-split: body pass is 99% near-field (45.6 of 45.7 s) **because
  kerneloffset radius inflation (Δr = `radius_inflation(E, kerneloffset,
  fmm_radius_tolerance)`, ~0.32R at kerneloffset 1e-3) forces geometrically
  FAR pairs into the direct lists** — precisely the pairs that are far enough
  for a point source/doublet swap. Eligibility on the body pass should
  therefore be high by construction.
- Wake pass: the direct particle-particle kernel spends heavily on
  exp/erf (`g_dgdr_gauserf`, `custom_erf`); beyond r/σ ≈ 4–6 the
  regularized kernel equals the singular one to better than the wake
  tolerance, reducing to the bare Biot–Savart 1/r² form.
- The threshold machinery already exists: `radius_inflation` in
  `src/FLOWPanel_abstractbody.jl` (~line 1157) computes exactly "the distance
  at which the singular kernel matches the offset-regularized direct kernel
  within tolerance" — the same logic, reused per-pair instead of per-radius.

## Phase 0 — eligibility census (GATE)

On the mature 023 state (`p018_cs_f1_l3p4` @ step 1034, 181k particles,
warm-started via `benchmark/p018_mature_wake_timing.jl` machinery — NOTE the
023 gotcha: run the full step head incl. `update_TE!` before ANY influence
evaluation, else the stale first wake row NaNs every target):

1. **Analytic thresholds.** For each kernel pair type, derive r_min(tol) as
   the MAX of two independent criteria (Ryan 2026-08-20: the regularization
   must be considered when choosing the singular approximation, not just the
   finite-size/point criterion):
   - *regularization criterion* — distance at which the regularized kernel
     matches the singular one within tol. Particles (gauserf): |g(r/σ)−1| ≤
     tol plus the matching gradient/hessian conditions (wake pass requests
     velocity + hessian at particle targets). Panels: the existing
     `radius_inflation` rules are exactly this — source/doublet regularize
     with compact support (singular beyond `kerneloffset` exactly);
     VortexRing (Vatistas n=2) needs r ≥ rc·(2/tol)^(1/4)
     (`src/FLOWPanel_elements_fmm.jl:982-991`);
   - *finite-size criterion* — point source+dipole vs the full panel
     integral, relative error ~O((a/r)²) with panel radius a.
2. **Direct-list walk.** Instrument (read-only) the FastMultipole direct
   lists for both production passes at BOTH knob sets — production
   (wake 4/0.4/50, body 8/0.4/20) and the 023 perturbed minima
   (wake 16/0.6/24, body 15/0.8/236, pending the p-ladder revision) — and
   histogram pair separations r/σ (particles) and r/a (panels).
3. **Report:** fraction of direct-list pair-interactions eligible at the
   pass tolerance × estimated per-pair speedup of the cheap kernel (micro-
   benchmark both kernel forms in isolation) ⇒ projected s/step saved per
   pass. Include the SFS sibling estimate: `zeta_gauserf` is exponentially
   ZERO beyond ~6σ, so eligible Estr pairs could be skipped outright —
   reported as a separate line, NOT implemented here.

Deliverable: `024_farfield_approximation_direct_kernels/phase_00_census.md`
with the eligibility table; STOP for Ryan's go/no-go.

## Phase 1 — implementation (gated on Phase 0 + Ryan)

Opt-in per-pair branch in the direct kernels:
- `src/FLOWPanel_elements_fmm.jl` panel kernels (`_induced`,
  `compute_source_dipole`, `_bound_vortex_velocity` callers): r² threshold →
  centroid point-source + point-dipole evaluation.
- FLOWVPM particle kernel (`g_dgdr_gauserf` call sites in
  `FLOWVPM_fmm.jl direct!`): r² threshold → g=1, dgdr=0 (singular
  Biot–Savart), consistent hessian branch.
- Threshold plumbed as an explicit tolerance knob, default OFF (exact
  behavior unchanged); direct/FMM backend alignment tests must pass with the
  feature off, and an A/B (CT + achieved-error metric from
  `benchmark/p023_fmm_tune.jl`) certifies it on.

## Contingency chain

- Census shows low eligibility on the wake pass (particles clustered within
  ~5σ — plausible given overlap ratios σ/h ≈ 4) → wake swap is a no-go;
  body pass may still proceed independently.
- Census shows high eligibility but micro-benchmark shows the branch
  mispredicts/vectorization breaks (SIMD direct loops) → restructure as
  per-leaf-pair (box-level) threshold instead of per-pair, using existing
  branch centers/radii — coarser but branch-free inside the loop.
- If the body-pass eligibility is dominated by radius-inflation-forced pairs,
  compare against simply revisiting the kerneloffset/inflation design (023
  lever 3) before adding kernel branches — whichever is simpler wins.

## Deferred

- SFS/Estr pair skipping beyond the zeta support radius (bigger lever than
  the velocity-kernel swap; touches 020-adjacent physics — Ryan's call).
- Changing kerneloffset itself or the `radius_inflation` policy (023 lever).
- GPU direct kernels.
