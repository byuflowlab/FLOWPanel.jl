# Task 4 — Green-reconstruction formulation

## Objective

Evaluate the velocity-only **explicit-potential** formulation
`GreenReconstruction(gauge=:area_mean, green_solver=nothing)` on the same frozen
finite-wake states as Tasks 2/3, and measure how much of the old finite-wake
lift gap it closes **toward the matching Task 3 direct-potential oracle** without
evaluating the free wake's scalar potential at body control points. Ranking
quantity is `f_oracle` (closure to Task 3, estimator-isolated); `f_semiinf`
(closure to Task 1 semi-infinite) is context. Primary ranking row: settled
rolled `Da=0.05c` ("march").

## Formulation

```julia
pnl.GreenReconstruction(gauge=:area_mean, recompute_interval=1, green_solver=nothing)
```

Run through the production `_steady_aerodynamics!` formulation hook
(`formulation`, `formulation_state`, `i_step`). The wake acts on the body through
velocity only; the wake-only source coefficients `σ = −u_f·n` give the body
source-panel potential `Sσ`, the body-only Green system `(I−B)q = Sσ` is solved
for the wake trace `q` under the area-weighted gauge (bordered LU), and the
explicit-potential system `G·μE = −S·σ0 − q` is solved reusing the finite-wake
`Backslash` LU. Physical circulation is `γ = C·μE` directly (the formulation
clears any wake correction; no affine Kutta shift is applied for this route).

- **Authoritative backend:** `DirectBackend` for solve, wake, and physical
  system. Rolled-state reconstruction (Task 2 march seed) uses the existing
  `FastMultipoleBackend(order=10, θ=0.4, leaf=100)`.
- **Green solver:** dense bordered LU of `[(I−B) a; aᵀ 0]` — the authoritative
  full-resolution reference. `(I−B)` is body-mesh-only (the attached wake is
  suppressed during its assembly), so **one** Green state is assembled once and
  reused across every wake state and both angles; only the mutable trace buffers
  (`q`, `last_recompute`) are refreshed per case.
- `BLAS.set_num_threads` pinned (recorded in `task4.config.toml`,
  `blas_threads`).

## Theory check

`docs/wake_solve_schemes.md` ("Suggested alternative formulation for
`semiinfinite_wake=false`") was compared against
`src/FLOWPanel_formulation.jl`. The Green identity `Sσ=(I−B)q`, the area-weighted
`(I−B)` gauge, the explicit-potential solve `G·μE=−S·σ0−q`, and body sources
carrying `σ0` all match the implementation. The only prior discrepancy — an
earlier remark that `(I−B)` was a low-rank update of the factored `G` — is
already reconciled in both the note (a standalone one-time bordered LU with the
gauge row) and the source docstring. **No theory correction was required.**

## Results

CL, difference from Task 1 (`ΔTask1`), the two-tier closure fractions, and the
**signed** post-correction difference `spc = CL_method − CL_(matching Task 3
route)`. `f_oracle` pairs single-shot↔Task 3 single-shot and iterated↔Task 3
iterated; a fraction >1 is overshoot. Baselines read from the stored
`comparison.csv` in each output directory.

### α = 3.94°  (Task 1 = 0.27476439; Task 3 rolled: single 0.26840468, iter 0.27338690)

| wake state | route | CL | ΔTask1 | f_semiinf | **f_oracle** | spc |
|---|---|---:|---:|---:|---:|---:|
| flat `Da=0.5c` | single | 0.2745158 | −0.090% | +0.548 | **+0.565** | −2.3e−4 |
| flat `Da=0.5c` | iter | 0.2744280 | −0.122% | +0.389 | +0.405 | −3.2e−4 |
| flat `Da=0.05c` | single | 0.2737001 | −0.387% | +0.667 | **+0.668** | −1.1e−3 |
| flat `Da=0.05c` | iter | 0.2708877 | −1.411% | −0.212 | −0.214 | −3.9e−3 |
| **rolled `Da=0.05c`** | **single** | **0.2673786** | −2.688% | +0.219 | **+0.668** | −1.0e−3 |
| **rolled `Da=0.05c`** | **iter** | **0.2695704** | −1.890% | +0.451 | **+0.527** | −3.8e−3 |

### α = 45°  (Task 1 = 2.65673304; Task 3 rolled: single 2.42665052, iter 2.45791467)

| wake state | route | CL | ΔTask1 | f_semiinf | **f_oracle** | spc |
|---|---|---:|---:|---:|---:|---:|
| flat `Da=0.5c` | single | 2.6545362 | −0.083% | +0.417 | **+0.428** | −2.1e−3 |
| flat `Da=0.5c` | iter | 2.6538139 | −0.110% | +0.226 | +0.235 | −2.8e−3 |
| flat `Da=0.05c` | single | 2.6485429 | −0.308% | +0.593 | **+0.593** | −8.2e−3 |
| flat `Da=0.05c` | iter | 2.6266386 | −1.133% | −0.497 | −0.500 | −3.0e−2 |
| **rolled `Da=0.05c`** | **single** | **2.4192214** | −8.940% | +0.041 | **+0.579** | −7.4e−3 |
| **rolled `Da=0.05c`** | **iter** | **2.4326675** | −8.434% | +0.096 | **+0.484** | −2.5e−2 |

Iteration counts (fixed-geometry projection to `≤1e-8` defect and lift change for
three consecutive steps): 3.94° flat 13, das005 45, march 44; 45° flat 14,
das005 50, march 46.

## Eligibility gates — all pass

1. **Formulation residuals ≤ 1e-10 (relative).** Across all 12 cells (2 angles ×
   3 states × 2 routes): worst Green bordered residual
   `‖(I−B)q + λa − Sσ‖ / ‖Sσ‖ = 6.92e-13`; worst explicit-potential residual
   `‖G·μE + S·σ0 + q‖ / ‖rhs‖ = 5.13e-15`. Gauge row `|aᵀq|` and all fields
   finite.
2. **Immutable-state audit.** Every frozen single-shot and every inner iteration
   preserved body nodes/cells/`Das`, wake nodes/connectivity/row-count/active +
   inactive strengths, and wake options (`_assert_frozen_state`). No violation
   raised.
3. **Velocity-only constraint.** GreenReconstruction evaluates only the body
   source-panel potential `S` (permitted); it never requests `scalar_potential`
   from the free-wake sources at body control points. The `q` vs Task 3 trace
   comparison below is a labeled diagnostic, outside the eligible solve path.

## Diagnostics

- **Backend control.** The `VelocityThroughSources`+`DirectBackend` control
  reproduces the stored old Task 2 CL for every state at both angles (e.g. rolled
  3.94° 0.26531206, rolled 45° 2.40900222). Backend (FMM vs Direct) is not the
  lever; the formulation differences are real.
- **Green trace vs Task 3 direct trace** (`q` vs `q_f`, area-mean removed,
  relative L2): 0.025 (flat) → 0.054 (rolled) at 3.94°; 0.069 → 0.160 at 45°.
  Consistent with the standalone `formulation_test.jl` value (0.0748). The
  reconstructed trace differs from the direct trace at the panel-discretization
  level, yet the resulting lift closes a large fraction of the oracle gap.
- Exterior-probe relative velocity differences vs Tasks 1/2/3 recorded per row.

## Conclusions (robustness pass)

1. **Primary rolled `Da=0.05c` row — robust gap-closure.** `f_oracle` is
   consistently positive and comparable across both angles and both routes:
   single 0.668 (3.94°) / 0.579 (45°), iterated 0.527 (3.94°) / 0.484 (45°).
   GreenReconstruction closes roughly **half to two-thirds** of the old
   finite-wake gap toward the matching Task 3 direct-potential result, and never
   overshoots it (spc < 0). This is the criterion the task is ranked on.
2. **`f_semiinf` is low on the rolled row (0.04–0.45) by construction, not
   estimator failure.** The Task 3 iterated oracle itself sits −0.50% (3.94°) /
   −7.48% (45°) below semi-infinite for the rolled wake, so `f_semiinf → 1` is
   unreachable by any frozen method. The high `f_oracle` alongside the low
   `f_semiinf` isolates the residual as **fixed rolled-wake geometry/strength
   error shared by every frozen method (including Task 3)**, not Green
   trace-estimation error.
3. **Flat `Da=0.5c`** has a sub-0.2% old gap (denominator ≈ 5.5e-4 CL), so its
   fractions are numerically soft; single-shot `f_oracle` 0.57/0.43 is a
   consistency check, not a ranking result.
4. **Caveat — flat `Da=0.05c` iterated overshoots below the old value**
   (`f_oracle` −0.21 at 3.94°, −0.50 at 45°) while its single-shot is strong
   (0.67 / 0.59). On this thin-`Da` flat state the fixed-geometry strength
   projection drives CL below both the old finite-wake and the oracle. It is a
   flat consistency-check state, not the primary rolled ranking row; the effect
   is flagged for the Task 5 comparison rather than treated as a rolled-row
   result.

## Artifacts (per output directory: `data/dirichlet_solve/` and `.../alpha45/`)

- `comparison.csv` — Task 4 rows keyed `green_reconstruction_single_shot|<case>`,
  `green_reconstruction_iterated|<case>`,
  `velocity_through_sources_direct_control|<case>` (additive columns; Task 1–3
  rows untouched).
- `task4_<case>_routes.csv` — single-shot vs iterated per state (`case` ∈ flat,
  flat_das005, march).
- `task4_<case>_iteration_history.csv` — fixed-geometry projection history.
- `task4.config.toml`, `task4.metadata.toml`, `task4_invariants.toml`.
- `dirichlet_task4_<case>_{single_shot,iterated}_{body,wake}.pvd` (+ VTU/filament
  collections).

## Reproduce

```bash
julia --project -t4 debug/dirichlet_solve/dirichlet_solve.jl --task task4 --alpha-deg 3.94
julia --project -t4 debug/dirichlet_solve/dirichlet_solve.jl --task task4 --alpha-deg 45 \
    --output-dir data/dirichlet_solve/alpha45
# resumable per case: --task task4-flat | task4-flat-das005 | task4-march
```

## Unregularized velocity-kernel follow-up

Additive Task 4 diagnostic (production default kernel unchanged; no wake was
remarched with `core_size=0`; no lift result replaced). It tests whether the
finite-wake Green-identity residual is caused by pairing a **regularized doublet
velocity** with the **unregularized analytic doublet potential**. The analytic
potential is `-strength[1]*tan_term` and never uses `reg_term`; the velocity uses
`val4 = (ri+rip1)/(ri*rip1*rho + reg_term)` with
`reg_term = regularize(distance, core_size)`, which is exactly zero when
`distance >= core_size`. Passing `core_size=0` selects the singular analytic
velocity everywhere while leaving the potential untouched.

Method: for each frozen state build **two otherwise-identical wakes** (same nodes,
strengths, active rows, options) differing only in `core_size` — `0.001` (default)
vs `0.0` — and re-evaluate `r = Sσ − (I−B)q_direct`. DirectBackend throughout;
`I−B`, area gauge, and bordered LU are identical between the two velocity modes.
For the semi-infinite attached sheet there is no `PanelWake.core_size`; the
velocity core is the body's active `kerneloffset`, zeroed only around the velocity
calls (potential kept at the original setting).

Commands:

```bash
julia --project -t4 debug/dirichlet_solve/dirichlet_solve.jl \
  --task green-regularization --alpha-deg 3.94 --output-dir data/dirichlet_solve
julia --project -t4 debug/dirichlet_solve/dirichlet_solve.jl \
  --task green-regularization --alpha-deg 45  --output-dir data/dirichlet_solve/alpha45
julia --project -t4 debug/dirichlet_solve/dirichlet_solve.jl \
  --task green-regularization-refinement --output-dir data/dirichlet_solve
```

Artifacts (wide paired rows, one per state):

- `data/dirichlet_solve/green_identity_regularization.csv` (3.94°)
- `data/dirichlet_solve/alpha45/green_identity_regularization.csv` (45°)
- `data/dirichlet_solve/green_identity_regularization_refinement.csv` (both angles,
  full grid ladder)

### Result — the velocity core is inert; the residual is unchanged

Production resolution (N = 6688), residual ratio = `r₀ / r₀.₀₀₁`:

| angle | state | route | resid (reg) | resid (unreg) | ratio | σ change |
|------:|-------|-------|------------:|--------------:|------:|---------:|
| 3.94° | flat        | seed / green / oracle | 5.797e-3 | 5.797e-3 | 1.0000 | 0 |
| 3.94° | flat_das005 | seed / green / oracle | 1.186e-2 | 1.186e-2 | 1.0000 | 0 |
| 3.94° | march       | seed / green / oracle | 1.270e-2 | 1.270e-2 | 1.0000 | 0 |
| 3.94° | semiinf     | task1                 | 2.545e-2 | 2.545e-2 | 1.0000 | 1.6e-9 |
| 45°   | flat        | seed / green / oracle | 3.100e-2 | 3.100e-2 | 1.0000 | 0 |
| 45°   | flat_das005 | seed / green / oracle | 6.671e-2 | 6.671e-2 | 1.0000 | 0 |
| 45°   | march       | seed / green / oracle | 8.815e-2 | 8.815e-2 | 1.0000 | 0 |
| 45°   | semiinf     | task1                 | 1.030e-1 | 1.030e-1 | 1.0000 | 1.5e-9 |

Refinement (both angles, grids 1744 / 3816 / 6688 / 10360): the paired residual
ratio is `1.000` at **every** resolution, and the regularized and unregularized
log-log slopes vs `h = 1/√N` are **bit-identical**:

| state | angle | slope (reg) | slope (unreg) |
|-------|------:|------------:|--------------:|
| flat_da0.5c  | 3.94° | 1.096 | 1.096 |
| flat_da0.5c  | 45°   | 1.444 | 1.444 |
| flat_da0.05c | 3.94° | 0.855 | 0.855 |
| flat_da0.05c | 45°   | 1.089 | 1.089 |
| semiinf      | 3.94° | 0.496 | 0.496 |
| semiinf      | 45°   | 0.662 | 0.662 |

(Slopes use the first three valid grids only; the (201,16,11)/N=10360 grid is
marked `mesh_valid_for_fit=false` — constant-mode defect 6.6e-3 — and excluded.)

For the free wakes the induced velocity is **bit-identical** with and without the
core (`sigma_relative_change = 0` exactly), because no body control point falls
within `core_size = 1e-3 m ≈ 0.003c` of any wake panel, so `reg_term ≡ 0` in both
modes. The semi-infinite sheet touches the trailing edge, so a handful of targets
sit within the tiny active `kerneloffset`; the change there is `σ`-level `~1.5e-9`
and `ratio − 1 ~ 6e-10` — negligible, as anticipated for an offset far below the
free-wake `0.001`.

### Positive control — the regularization is wired in, just inactive at `0.001`

The paired test re-evaluates the frozen wake through the **production** velocity
kernel `pnl.influence!((body,),(wake,); velocity=true)` → `_induced` →
`val4 = (ri+rip1)/(ri*rip1*rho + reg_term)` (the same path `simulate!` uses for
wake→body velocity), but standalone rather than through `simulate!`. To confirm
`reg_term` actually reaches `σ` and the `0.001` null is a separation-scale effect
(not a gated-off kernel), sweep `core_size` on the closest-to-surface state
(`flat_das005` seed, N=6688, 3.94°):

| core_size (m) | core/c | resid_rel | ratio vs `0.001` | σ change |
|--------------:|-------:|----------:|-----------------:|---------:|
| 0.001 | 0.003 | 1.186e-2 | 1.00 | 0 |
| 0.005 | 0.016 | 1.186e-2 | 1.00 | 0 |
| 0.010 | 0.033 | 1.186e-2 | 1.00 | 0 |
| 0.020 | 0.066 | 4.089e-2 | **3.45** | 0.63 |
| 0.050 | 0.164 | 3.552e-1 | 30.0 | 0.94 |
| 0.100 | 0.328 | 7.955e-1 | 67.1 | 1.01 |
| 0.200 | 0.656 | 1.765e+0 | 148.8 | 1.03 |
| 0.500 | 1.640 | 1.019e+1 | 859.3 | 1.03 |

The core activates once it exceeds ~0.02 m ≈ 0.066c — exactly the geometric scale
at which it reaches the nearest body control points (the das=0.05c free wake begins
~0.015 m behind the TE) — and then moves the residual by orders of magnitude. So
`reg_term` genuinely feeds `σ`; the production `0.001 m ≈ 0.003c` is simply below
the minimum body–wake separation, so it is identically zero. The null result is
validated, not an artifact of a bypassed kernel.

### Gates (all pass)

- **Only the velocity core changes** — paired wakes verified identical apart from
  `core_size` (`_assert_wake_pair_equal`).
- **Potential invariance (gate 2)** — `q_direct_relative_change = 0.0` (≤ 1e-13)
  for every panel-wake pair; the analytic potential ignores `reg_term`.
- **Regularized baseline reproduction (gate 3)** — regularized residuals matched
  the established `green_identity_residual.csv` / `green_identity_refinement.csv`
  rows to `rtol=1e-10, atol=1e-13`.
- **Solve integrity (gate 4)** — `recon_solve_rel ≤ 1.2e-12` everywhere (≪ 1e-8).
- **Finite fields (gate 5)** — no singular encounter; all σ, q, r finite.
- **Immutable-state audit (gate 6)** — `_state_snapshot`/`_assert_frozen_state`
  around each evaluation plus the cross-pair equality audit.
- **Defective fine mesh (gate 7)** — N=10360 excluded from slope fits.

### Conclusion

Velocity/velocity-potential regularization inconsistency is **not** a contributor
to the finite-wake Green-identity residual at these grids: removing the velocity
core (`core_size=0`) reproduces every residual to machine precision. The residual
is the **centroid-sampled constant normal-flux / constant-panel surface
representation** error (with near-singular sensitivity as the state-dependent
amplifier at the trailing edge and short `Da=0.05c` transition). It converges under
refinement on all states, so the lever is resolution and a better near-TE surface
representation — not regularization. See `green_identity_analysis.md` for the
mechanism discussion.
