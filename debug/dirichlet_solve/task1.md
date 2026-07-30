# Task 1 — Semi-Infinite Attached-Wake Baseline

## Final Formulation

Theory gate:

- `docs/wake_solve_schemes.md`, “Discrete diagnostic construction: Task 1
  semi-infinite baseline”
- `docs/dirichlet_potential_theory.md`, especially “Symbols and signs” and
  “FLOWPanel's exact-centroid operator”

Implemented theory amendment:

- Prescribed fields and exact system
  $\mathcal G_\infty\mu_\infty=-\mathcal S\sigma_0$
- Implicit $\mathcal B/\mathcal W_\infty\mathcal C$ assembly
- FLOWPanel signs and exact-centroid self limit
- Gauge treatment and $\mathcal Wc=0$
- Required invariants, outputs, and limitations

FLOWPanel conventions:

- Source strength: $\sigma=-\boldsymbol u\boldsymbol\cdot\boldsymbol n$
- Doublet strength: $\mu=\phi_{\rm inside}-\phi_{\rm outside}$
- Dirichlet control points: exact panel centroids
- Doublet self kernel: interior $+\mu/2$ limit; no solver diagonal addition
- Body operator: $\mathcal B$
- Paired Kutta map: $\mathcal C$, upper minus lower
- Lifting operator: $\mathcal G=\mathcal B+\mathcal W\mathcal C$

## Final Implementation

- Pitching-wing watertight geometry and prescribed flow inputs reused.
- Semi-infinite attached-wake system assembled and solved with direct diagnostic
  influence evaluation.
- Watertightness, paired shedding edges, panel indices, and
  $\mathcal C\mathbf1=0$ verified.
- Direct linear residual required finite and below tolerance.
- Transition circulation and steady pressure-derived $C_L$ recorded.
- Common 16-point exterior velocity probe field saved.
- Compact comparison CSV, labeled VTK, configuration, metadata, and one-step
  convergence history written.
- No free wake or frozen-step manifest in Task 1.

## Successful Commands

```bash
JULIA_NUM_THREADS=2 julia --project=. debug/dirichlet_solve/dirichlet_solve.jl task1
JULIA_NUM_THREADS=2 julia --project=. debug/dirichlet_solve/dirichlet_solve.jl task1-grid
JULIA_NUM_THREADS=2 julia --project=. -e 'include("test/runtests_unit_solver.jl")'
JULIA_NUM_THREADS=2 julia --project=. -e 'include("test/runtests_unit_wake.jl")'
JULIA_NUM_THREADS=2 julia --project=. -e 'using Test; import FLOWPanel as pnl; using WriteVTK; include("test/test_helpers.jl"); include("test/runtests_unit_liftingbody.jl")'
git diff --check
```

Results: solver 84/84; wake 82/82; lifting body 40/40; diff check clean.

## Artifacts

- Root: `data/dirichlet_solve/`
- Compact result: `comparison.csv`
- Labeled VTK collection: `dirichlet_task1_semiinfinite_body1.pvd`
- Configuration: `task1.config.toml`
- Metadata: `task1.metadata.toml`
- Direct-solve history: `task1_convergence.csv`
- Exterior probes: `task1_exterior_probes.csv`
- Grid convergence: `task1_grid_convergence.csv`

## Baseline Results

- Body panels: 6,688
- Paired shedding edges: 13
- Watertight mesh
- $\lVert\mathcal C\mathbf1\rVert_\infty=0$
- Absolute residual: $2.1537824326493315\times10^{-13}$
- Relative residual: $2.4807313750444347\times10^{-15}$
- Transition circulation:
  - minimum: $-5.0555235231$
  - maximum: $-2.5732949078$
  - mean: $-4.2824577295$
  - L2 norm: $15.7176054647$
- Steady pressure-derived $C_L=0.2747643938$
- Sixteen finite common exterior probes
- Labeled VTK fields: `normals`, `potential`, `velocity`, `sigma`, `mu`,
  `gauge pressure`, and `F`
- Baseline relative exterior difference: `NaN` (not applicable)
- Finite free-wake equivalence: no conclusion in Task 1

## Grid Results

Proportional sequence:

| Panels | $C_L$ | Relative residual |
|---:|---:|---:|
| 1,744 | 0.2784515 | $1.24\times10^{-15}$ |
| 3,816 | 0.2760592 | $1.95\times10^{-15}$ |
| 6,688 | 0.2747644 | $2.48\times10^{-15}$ |
| 10,360 | 0.2741234 | $3.11\times10^{-15}$ |
