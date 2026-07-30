# Green-reconstruction identity — refinement + conditioning analysis

## Objective

Localize **where and why** the finite-wake `GreenReconstruction` trace
reconstruction degrades, so the path forward can be chosen (finer mesh /
near-field quadrature vs. regularization vs. accept the trace error since lift
already works per Task 4). Two read-only probes on frozen wakes, DirectBackend
throughout, no time-marching.

- **Probe A (`green-refinement`)** — the discrete Green third-identity residual
  `r = Sσ − (I−B)q_direct` and the reconstructed-vs-direct trace
  `‖q_recon − q_direct‖`, measured across the Task-1 grid ladder
  (1744 → 3816 → 6688 → 10360 panels) on three deterministically rebuildable
  states: flat `Da=0.5c` (far-field control), flat `Da=0.05c` (near-field), and
  the semi-infinite attached sheet (attached / worst near-field). Both angles
  (3.94°, 45°) at each resolution reuse the same body-mesh-only Green LU.
- **Probe B (`green-conditioning`)** — the singular spectrum of `ImB = I − B`
  per resolution, and the smallest-10 right-singular-vector character at
  production (6688): projection onto the constant and trailing-edge
  localization.

## Verification gates (all pass)

- Refinement reproduces the loader-based `green_identity_residual.csv` at 6688,
  flat `Da=0.5c`: **0.03103 at 45°** (expected 0.0310), **0.005804 at 3.94°**
  (expected 0.0058). The deterministic ladder and the loader diagnostic agree.
- The constant is the smallest singular mode: at 6688 the smallest right-singular
  vector has `const_projection = 1.0000`, and the forward defect
  `‖(I−B)·1‖₂/√N = 4.3e-11` (≈ the prior `‖(I−B)·1‖∞ ≈ 7e-10`, up to the
  RMS-vs-∞ norm factor). `svd` reports its value as ~7e-18, a LAPACK
  rank-deficiency floor; the reliable forward matvec (4.3e-11) is the true σ_min.
- Every case: bordered `recon_solve_rel ≤ 8e-13`; finite fields.

## Probe A — the residual converges; it does not plateau

Log–log slope of `r` (green_identity_residual_rel) vs `h = 1/√N`, fit over the
three clean grids **1744 / 3816 / 6688** (the 10360 grid is excluded — see
below):

| state | 3.94° | 45° |
|---|---:|---:|
| flat `Da=0.5c` (far-field) | 1.10 | 1.44 |
| flat `Da=0.05c` (near-field) | 0.86 | 1.09 |
| semi-infinite (attached) | 0.50 | 0.66 |

All slopes are **positive** — the identity residual decays algebraically under
refinement on every state, so the error is **discretization-limited, not a
near-singular plateau**. The far-field flat state converges fastest (~`h^1.1–1.4`);
the two **near-field** states (`Da=0.05c` and the attached sheet) converge
**slower** (~`h^0.5–1.1`), so the wake→body near-field quadrature is the
slower-converging component, but it still converges. The reconstructed-vs-direct
trace `q_recon_vs_direct` tracks `r` closely (its own slopes ~0.5–1.0), i.e. the
trace error scales with the discretization residual, not with the operator
conditioning.

**10360-panel grid excluded — mesh defect.** At `(201,16,11)` the constant-mode
defect `‖(I−B)·1‖₂/√N` jumps to **6.6e-3**, eight orders above the 4.3e-11 seen
at 6688 (and 1.6e-10, 6.2e-14 at the coarser grids). `B` rows no longer sum to
≈1, so that discretization has bad panel influence coefficients (likely
sliver/degenerate panels at the finer round endcap). Its residuals are
contaminated (e.g. 45° flat `Da=0.5c` rises to 0.067, breaking monotonicity) and
are not used in the rates above. **That resolution's panelization should be
checked before it is trusted.**

## Probe B — one isolated near-null mode (the constant); no tiny-σ cluster

`ImB = I − B` singular spectrum (gauge mode = the constant, removed by the area
gauge):

| N | σ_max | σ_min (constant) | σ_2 | cond = σ_max/σ_2 | #(σ<1e-3·σ_max) |
|---:|---:|---:|---:|---:|---:|
| 1744 | 0.958 | ~1e-13* | 1.03e-2 | 92.6 | 1 |
| 3816 | 0.972 | ~2e-10* | 9.76e-3 | 99.6 | 1 |
| 6688 | 0.986 | ~4.3e-11* | 9.45e-3 | 104.3 | 1 |
| 10360 | 0.998 | 6.6e-3 (defective) | 9.25e-3 | 108.0 | 1 |

\*constant-mode defect from the forward matvec (`const_mode_check`), the reliable
σ_min; the smallest mode is the constant (`const_projection = 1.0`).

- **The constant is the only dangerous mode.** `#(σ < 1e-3·σ_max) = 1` at every
  resolution — that one is the constant, which the area-weighted gauge removes.
  There is **no cluster of tiny singular values**: after the constant, σ_2 ≈ 9.5e-3
  and the gauge-fixed conditioning is a **stable ~100** across resolution.
- **Small-mode character (6688).** Modes 2–10 are orthogonal to the constant
  (`const_projection ~ 1e-14`) with σ from 9.5e-3 to 5.8e-2. Their trailing-edge
  localization rises smoothly from **0.3% (mode 2) to 19% (mode 10)** — the
  lowest non-gauge modes are fairly global/smooth; TE concentration grows only in
  the higher modes. So the near-null space is not a TE-localized cluster; it is
  the single global constant.

**⇒ The gauge suffices. The trace is well-determined; regularization is not
needed.**

## Synthesis — predicted vs. observed trace gap

The worst-case bound `‖q_recon − q_direct‖ ≤ ‖r‖ / σ_2` (relative form
`resid_rel · ‖Sσ‖ / (σ_2 · ‖q_direct‖)`, σ_2 at 6688) vs. the observed
`q_recon_vs_direct`:

| angle | state (rolled/near-field) | predicted (‖r‖/σ_2) | observed | ratio |
|---|---|---:|---:|---:|
| 3.94° | flat | 0.492 | 0.0254 | 19.4 |
| 3.94° | flat `Da=0.05c` | 1.08 | 0.0480 | 22.4 |
| 3.94° | rolled `Da=0.05c` (march) | 1.17 | 0.0541 | 21.6 |
| 45° | flat | 0.689 | 0.0690 | 10.0 |
| 45° | flat `Da=0.05c` | 1.94 | 0.160 | 12.1 |
| 45° | rolled `Da=0.05c` (march) | 2.93 | 0.225 | 13.0 |

The observed trace gap is **below the σ_2 worst-case bound everywhere** (ratio
10–40× ⇒ bound respected but loose). The residual `r` is **not** aligned with the
σ_2 near-null direction, so the actual conditioning amplification is mild — the
trace gap is set by `r` projected onto the well-conditioned bulk of the spectrum
and therefore **scales with `r`**, not with `1/σ_2`. This refines the original
"trace gap ≈ ‖r‖/σ_2" hypothesis: the mechanism is discretization residual ×
conditioning, but the amplification factor is modest and refinement-decaying, not
the near-singular worst case.

## Clarifications (two conceptual points)

### What "near-field wake→body quadrature" means

The identity is evaluated at body **control points**. Two terms need the wake's
influence at each control point: `q_direct` (wake scalar **potential**, kernel
`1/r`) and `Sσ` with `σ = −u_f·n` (wake induced **velocity**, kernel `1/r²`).
The panel influences are **analytic constant-panel formulas** (Hess–Smith /
Katz–Plotkin), not numerical quadrature, so the near-field error is *not* an
integration-rule error. Three distinct mechanisms could produce the residual:

1. **Centroid-sampled constant normal-flux representation.** Each panel carries a
   single constant strength and `σ` is sampled from the velocity at one control
   point. The exact identity holds for the continuous field; the discrete
   constant-per-panel / one-sample-per-panel representation of a field that varies
   across the panel is the leading `O(h)` truncation. This is largest where the
   solution varies fastest relative to the panel size — at the trailing edge where
   the shed sheet attaches, and for the short `Da=0.05c` transition panel sitting
   essentially on the surface — which is why the residual splits by state: flat
   `Da=0.5c` (wake half a chord downstream) converges at ~`h^1.1–1.4`, while
   `Da=0.05c` and the attached sheet (wake material on/near the surface) converge
   at only ~`h^0.5`.
2. **Near-singular numerical sensitivity** of the analytic velocity formula when a
   body target sits very close to a wake panel edge (`1/r²` grows steeply).
3. **Velocity/potential inconsistency from the velocity-only regularization core.**
   The analytic doublet **potential** is singular/unregularized, but the doublet
   **velocity** adds a `reg_term = regularize(distance, core_size)` denominator
   term. Using a regularized velocity with an unregularized potential is formally
   inconsistent and could, in principle, break the identity near the wake.

The unregularized velocity-kernel follow-up (`task4.md`, "Unregularized velocity-
kernel follow-up") isolates mechanism 3 and rules it out at these grids: with
`core_size=0` the residual is **bit-identical** to the
regularized `core_size=0.001` result (ratio `1.000`, `σ` change exactly `0` for
the free wakes) because no body control point falls within `1e-3 m ≈ 0.003c` of a
wake panel, so `reg_term` is identically zero either way. The residual is therefore
mechanism 1 (with mechanism 2 as the state-dependent amplifier), **not** the
velocity/potential inconsistency. The lever is resolution and a better constant-
panel surface representation near the trailing edge (higher-order strength, local
refinement), which shrinks `r` at fixed `h`.

### Why the near-null mode does not require (additional) regularization

The worry is legitimate in general: inverting through a σ_min ≈ 4e-11 would
divide some residual component by 4e-11 and blow it up ~10 orders. It does not
happen here because we never invert through that singular value:

1. **Bordered solve, not `(I−B)q = Sσ`.** The gauge appends the area vector `a`:
   `K = [[(I−B), a]; [aᵀ, 0]]`, solving `K·[q; λ] = [Sσ; 0]`. The near-null
   direction of `(I−B)` is the constant (to 13 digits), and `aᵀq = 0` pins
   exactly that direction; `λ` is a Lagrange multiplier, not a free unknown
   floating in the near-null space. `K` is well-conditioned even though `(I−B)`
   is singular — hence the bordered residual `recon_solve_rel ≈ 1e-13` in every
   case. The extra DOF **removes** the rank deficiency; the 4e-11 never enters a
   denominator. This *is* regularization, but **exact deflation of a known null
   space** (the constant), so it adds **zero bias** — unlike Tikhonov / truncated
   SVD / damping, which is the "regularization" that is *not* needed here.
2. **The dangerous subspace is exactly 1-D and the gauge catches all of it.**
   Probe B: `#(σ < 1e-3·σ_max) = 1` at every resolution, that mode has
   `const_projection = 1.0`, and after it `σ_2 ≈ 9.5e-3` (gauge-fixed cond ≈ 100,
   stable in `h`). The surviving residual is amplified by at most ~`1/σ_2 ≈ 100`
   (empirically far less — the 10–40× loose bound). Had there been a *cluster* of
   tiny σ's, the single-mode gauge would leave an ill-conditioned subspace and the
   trace would be under-determined ⇒ real regularization needed; that case is
   ruled out.
3. **The constant is physically inert.** `q` is a potential; the observable is
   its gradient (velocity → pressure → lift). The constant/gauge mode has zero
   gradient, so any residual along it is invisible to forces — the one
   numerically singular direction is also the one that carries no aerodynamics.

Net: error amplification from a near-null mode requires inverting through its tiny
singular value; the gauge/bordering guarantees we never do that for the constant,
there is no second small singular value (cond ≈ 100 after deflation), and the
constant is gradient-free anyway. So the leftover trace error is `O(r)` — a
discretization residual that decays under refinement — not a conditioning blow-up.

## Conclusion and path forward

1. **Method is sound.** This is a forward discretization residual, not a solve or
   gauge failure. `(I−B)` has a single isolated near-null mode (the constant),
   the gauge removes it, and the trace is well-determined (gauge-fixed cond ≈ 100,
   no tiny-σ cluster).
2. **The residual converges under refinement** on every state (no plateau);
   near-field states (`Da=0.05c`, attached sheet) are the slow part at ~`h^0.5–1.1`.
   **The lever is resolution + a better constant-panel surface representation near
   the trailing edge** (centroid-sampled constant normal-flux error), not
   regularization. The follow-up (`task4.md`) confirms the velocity-only core is
   *inert* at these grids (`core_size=0` reproduces the residual to machine
   precision), so velocity/potential inconsistency is not a contributor — consistent
   with Task 4, where the trace error is already benign for lift.
3. **Fix the (201,16,11) mesh** (constant-mode defect 6.6e-3) before using that
   resolution; its Green operator is polluted and its residuals were excluded here.

## Artifacts

- `data/dirichlet_solve/green_identity_refinement.csv` — per (grid, angle, state):
  residual, trace gap, compat defect, σ_min proxy, recon-solve residual.
- `data/dirichlet_solve/green_identity_conditioning.csv` — per grid: σ_max, σ_min,
  σ_2…σ_10, gauge-fixed cond, tiny-σ counts, constant-mode check.
- `data/dirichlet_solve/green_identity_smallmodes.csv` — 6688: smallest-10 modes'
  σ, constant projection, TE localization.

Reproduce:

```bash
julia --project debug/dirichlet_solve/dirichlet_solve.jl --task green-refinement
julia --project debug/dirichlet_solve/dirichlet_solve.jl --task green-conditioning
```
