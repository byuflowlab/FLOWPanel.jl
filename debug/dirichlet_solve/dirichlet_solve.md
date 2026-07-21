# Dirichlet Wake-Coupling Verification

## Agent Instructions

- Read `CLAUDE.md`; then every task-relevant routed policy.
- Read this index; then only the active task file linked below.
- Before starting any task, identify its main goal—the criterion the user will
  ultimately use to judge whether the task is approved—and explicitly verify
  that criterion with the user. Do not begin task work until it is confirmed.
- First task with unchecked `User approved`; no later-task work.
- Substantive solver formulation: complete cited theory section; compare against
  current FLOWPanel code; correct theory before implementation.
- Item files: final working product only—formulation, successful commands,
  artifacts, measurements, conclusions, theory corrections.
- No failure logs, attempted-command histories, superseded approaches, or
  chronological narratives.
- This index: status, item links, terse takeaways only; no duplicated detail.
- Never check `User approved` without explicit user instruction.
- Task 2 onward: frozen finite-wake state copied and verified before each solve;
  only identified transition-panel strengths variable.
- Mechanical plotting, CSV/TOML, VTK, state-loading helpers: no separate theory
  amendment.

## Common Setup

- Watertight AR=4 rectangular NACA 0015
- `n_airfoil=161`; `n_span=13`; `n_endcap=9`
- Chord 1 ft; span 4 ft; $U_\infty=330.2\ \mathrm{ft/s}$;
  $\alpha=3.94^\circ$
- Finite-wake march: `c_per_dt=0.5`; $|D_a|=0.05c$; constant angle; no pitching
- Outputs: `data/dirichlet_solve/`
- Compact comparison: `data/dirichlet_solve/comparison.csv`

## Task 1 — Semi-infinite attached-wake baseline

- [x] Theory checked and documented
- [x] Implementation complete
- Item: [task1.md](task1.md)
- Takeaways:
  - Exact-centroid direct solve; paired Kutta map; finite residual/probes
  - Baseline: 6,688 panels; 13 shedding edges; $C_L=0.2747643938$
  - Four-case proportional grid convergence: 1,744 to 10,360 panels
  - Finite free-wake equivalence: not evaluated
- [x] User approved

## Task 2 — Current FLOWPanel velocity-through-sources route

- [x] Theory checked and documented
- [x] Implementation complete
- Item: [task2.md](task2.md)
- Takeaways:
  - Flat sequence converged by $64c$: $C_L=0.2742139824$; Task 1 difference
    $-5.5041\times10^{-4}$
  - With $|D_a|=0.05c$ and unchanged $0.5c$ free rows, the sequence converged
    by $63.55c$: $C_L=0.2715662450$; Task 1 difference
    $-3.1981\times10^{-3}$
  - With the rolled march corrected to $|D_a|=0.05c$, it settled at
    $40c/U_\infty$: $C_L=0.2653120609$; Task 1 difference
    $-9.4523\times10^{-3}$ ($-3.440\%$)
  - At $45^\circ$, both flat sequences still converged by their terminal
    lengths; their relative baseline differences were $-0.142\%$ and
    $-0.757\%$
  - The corrected $45^\circ$ rolled wake settled at $40c/U_\infty$ but differed
    from its semi-infinite baseline by $-9.325\%$ and had a 0.543 relative
    exterior-probe difference
  - The former $|D_a|=0.5c$ marches understated the discrepancies: $-0.325\%$
    at $3.94^\circ$ and $-1.961\%$ at $45^\circ$
  - All reporting solves preserved frozen body/wake geometry and prescribed
    strengths exactly; no production behavior changed
- [x] User approved

## Task 3 — Direct fixed-wake potential

- [x] Theory checked and documented
- [x] Implementation complete
- Item: [task3.md](task3.md)
- Takeaways:
  - Direct prescribed-panel potential and manual finite-body LU implemented as
    a diagnostic only; all frozen-state and $10^{-10}$ residual gates passed
  - Both terminal flat wakes recover Task 1 lift within $8.0\times10^{-5}$
    relative at $3.94^\circ$ and $5.4\times10^{-5}$ at $45^\circ$ after the
    fixed-geometry strength projection
  - Rolled-wake projection converged and matched the augmented oracle, but
    retained Task 1 differences of $-0.5013\%$ and $-7.4836\%$
  - All six terminal iterations converged without reducing $\omega$; no active
    final filament or vector-potential-only source entered Task 3
- [x] User approved

## Task 4 — Green reconstruction

- [x] Theory checked and documented
- [x] Implementation complete
- Item: [task4.md](task4.md)
- Takeaways:
  - `GreenReconstruction(gauge=:area_mean, green_solver=nothing)` via the
    production formulation hook; dense bordered-LU of body-mesh-only `(I−B)`
    built once and reused across all states/angles; DirectBackend authoritative
  - Theory matched code; no correction needed (the `(I−B)` low-rank remark was
    already reconciled to a standalone bordered LU)
  - All 3 gates pass on all 12 cells (2 angles × 3 states × 2 routes): worst
    Green bordered residual `6.92e-13`, worst explicit-potential residual
    `5.13e-15`; frozen-state and velocity-only preserved
  - **Primary rolled `Da=0.05c`, robust gap-closure** toward the matching Task 3
    oracle: single `f_oracle` 0.668 (3.94°) / 0.579 (45°), iterated 0.527 / 0.484
    — consistent sign & magnitude across both angles and routes; never overshoots
  - `f_semiinf` low on the rolled row (0.04–0.45) is structural (the Task 3
    oracle is itself −0.50%/−7.48% short of semi-infinite there), not estimator
    failure — residual is shared fixed-wake geometry error, not trace estimation
  - Flat `Da=0.5c` gap is sub-0.2% (soft denominator); flat `Da=0.05c` iterated
    overshoots below the old value (`f_oracle` −0.21/−0.50) though its single-shot
    is strong (0.67/0.59) — flagged for the Task 5 comparison
- [x] User approved

## Task 5 — Velocity-only trace correction

- [x] Theory checked and documented
- [x] Implementation complete
- Item: [task5.md](task5.md)
- Takeaways:
  - `TraceCorrected(estimator=:line_integral)` validated on all 36 primary cells
    (2 angles × 3 states × 3 variants × 2 routes); every residual, circulation,
    frozen-state, finite-field, and velocity-only gate passed
  - Straight Simpson ranks first by worst rolled-wake trace error (`4.2872e-11`),
    deformed Simpson is tied (`4.2881e-11`), and trapezoid is third (`1.5284e-7`)
  - Deformed-path depth `0.5–2.0` is immaterial; lift is also quadrature/path
    invariant at shown precision, so the remaining Task 3 lift gap is not trace
    estimator error
  - No theory correction was needed and no self-consistent wake march was run
- [ ] User approved

## Green-identity refinement + conditioning probes

- [x] Theory checked and documented
- [x] Implementation complete
- Item: [green_identity_analysis.md](green_identity_analysis.md)
  (`--task green-refinement`, `--task green-conditioning`)
- Takeaways:
  - The trace degradation is a **forward discretization residual**, not a
    solve/gauge failure. `r = Sσ − (I−B)q_direct` **converges algebraically**
    under refinement on every state (no plateau): far-field flat ~`h^1.1–1.4`,
    near-field (`Da=0.05c`/attached) slower at ~`h^0.5–1.1`
  - `(I−B)` has a **single isolated near-null mode = the constant**
    (`const_projection=1.0`, σ_min≈4e-11 at 6688); the area gauge removes it,
    σ_2≈9.5e-3, gauge-fixed cond ≈ 100 stable across resolution, **no tiny-σ
    cluster** ⇒ gauge suffices, trace well-determined, **no regularization needed**
  - Worst-case `‖r‖/σ_2` **over**predicts the trace gap 10–40× (bound respected
    but loose): `r` avoids the σ_2 direction, so the gap scales with `r`, not `1/σ_2`
  - **Path forward:** resolution + better near-field wake→body quadrature (the
    slow ~`h^0.5` part); trace error benign for lift per Task 4
  - **⚠ (201,16,11) / 10360-panel mesh is defective** — constant-mode defect
    6.6e-3 (8 orders above coarser grids); excluded from rates, needs a mesh check
  - **Unregularized velocity-kernel follow-up** (`task4.md`,
    `--task green-regularization[-refinement]`): removing the doublet-velocity core
    (`core_size=0`) reproduces every residual **bit-for-bit** (ratio 1.000, identical
    slopes) — the core is inert because no body CP is within `1e-3 m ≈ 0.003c` of a
    wake panel. So velocity/potential regularization inconsistency is **not** a
    cause; the residual is centroid-sampled constant-panel representation error
- [x] User approved
