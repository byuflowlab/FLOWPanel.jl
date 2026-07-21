# Task 3 — Direct Fixed-Wake Potential

## Final Formulation

Theory gate: `docs/wake_solve_schemes.md`, “Discrete diagnostic construction:
Task 3 direct fixed-wake potential.”

The diagnostic solves

$$
\mathcal G_\Delta\mu_E=-\mathcal S\sigma_0-q_f,
\qquad
\sigma_0=-U_\infty\cdot n,
$$

where $q_f$ is evaluated directly from the prescribed constant-doublet
`PanelWake` panels at exact body centroids. Wake velocity is not converted into
body sources and no mean is removed from $q_f$. The existing finite-body LU is
applied manually because single-body `solve!` clears preassembled potential.

Task 2 geometry and state construction are reused for all seven lengths in each
flat sequence and for the selected settled rolled state. Selected step,
geometry, options, nodes, active row count, and initial strength hashes match.
The active source tuple contains only the `PanelWake`; every case asserts zero
active final filaments and no vector-potential-only source. Every single shot
and inner solve preserves its frozen state exactly.

On the three terminal geometries, all active row strengths are projected toward
$\mathcal C\mu_E$ at fixed geometry. The maximum rowwise defect and relative
lift change must both be at most $10^{-8}$ for three consecutive iterations.
Relaxation starts at one and follows the documented 5%-decrease rule with a
$1/16$ floor. The 200-iteration nonconvergence status is explicit.

The converged result is checked against

$$
(\mathcal G_\Delta+P\mathcal C)\mu_*=-\mathcal S\sigma_0,
$$

where each column of $P$ is the body-centroid potential of one unit-strength
wake strip. A Woodbury reduction uses the existing LU while the displayed
augmented residual is evaluated directly.

## Successful Commands

```bash
JULIA_NUM_THREADS=2 julia --project=. debug/dirichlet_solve/dirichlet_solve.jl task3-flat
JULIA_NUM_THREADS=2 julia --project=. debug/dirichlet_solve/dirichlet_solve.jl task3-flat-das005
JULIA_NUM_THREADS=2 julia --project=. debug/dirichlet_solve/dirichlet_solve.jl task3-march
JULIA_NUM_THREADS=2 julia --project=. debug/dirichlet_solve/dirichlet_solve.jl task3-flat --alpha-deg=45 --output-dir=data/dirichlet_solve/alpha45
JULIA_NUM_THREADS=2 julia --project=. debug/dirichlet_solve/dirichlet_solve.jl task3-flat-das005 --alpha-deg=45 --output-dir=data/dirichlet_solve/alpha45
JULIA_NUM_THREADS=2 julia --project=. debug/dirichlet_solve/dirichlet_solve.jl task3-march --alpha-deg=45 --output-dir=data/dirichlet_solve/alpha45
JULIA_NUM_THREADS=2 julia --project=. -e 'include("test/runtests_unit_solver.jl")'
JULIA_NUM_THREADS=2 julia --project=. -e 'include("test/runtests_unit_wake.jl")'
JULIA_NUM_THREADS=2 julia --project=. -e 'using Test; import FLOWPanel as pnl; using WriteVTK; include("test/test_helpers.jl"); include("test/runtests_unit_liftingbody.jl")'
JULIA_NUM_THREADS=2 julia --project=. -e 'include("test/runtests_unit_postprocess.jl")'
git diff --check
```

Results: solver 84/84; wake 82/82; lifting body 40/40; post-processing
405/405. Both generated artifact roots passed an additional schema/residual/
terminal-convergence audit. The Task 3 files have no whitespace errors. The
repository-wide `git diff --check` reports only the pre-existing unrelated
`CLAUDE.md:51: new blank line at EOF` in the dirty worktree.

## Artifacts

The normal root is `data/dirichlet_solve/`; the complete 45-degree matrix has
the same layout under `data/dirichlet_solve/alpha45/`.

- `comparison.csv`
- `task3_flat_single_shot.csv`
- `task3_flat_das005_single_shot.csv`
- `task3_march_single_shot.csv`
- `task3_flat_iteration_history.csv`
- `task3_flat_das005_iteration_history.csv`
- `task3_march_iteration_history.csv`
- `task3_invariants.toml`
- `task3.config.toml`
- `task3.metadata.toml`
- terminal single-shot and iterated body/wake PVD collections named
  `dirichlet_task3_{flat,flat_das005,march}_{single_shot,iterated}_{body,wake}.pvd`

## Measurements

### Frozen single shots

| Angle | Geometry | Rows | $C_L$ | $C_L-C_{L,1}$ | $C_L-C_{L,2}$ | Task 1 probe diff. | Task 2 probe diff. | Relative residual |
|---:|---|---:|---:|---:|---:|---:|---:|---:|
| $3.94^\circ$ | $64c$ flat | 127 | 0.2747484319 | $-1.5962\times10^{-5}$ | $+5.3445\times10^{-4}$ | $5.0436\times10^{-6}$ | $3.7932\times10^{-5}$ | $2.5465\times10^{-15}$ |
| $3.94^\circ$ | $63.55c$ short-transition flat | 127 | 0.2747584857 | $-5.9082\times10^{-6}$ | $+3.1922\times10^{-3}$ | $5.2629\times10^{-6}$ | $4.9901\times10^{-5}$ | $2.3762\times10^{-15}$ |
| $3.94^\circ$ | settled rolled wake | 80 | 0.2684046804 | $-6.3597\times10^{-3}$ | $+3.0926\times10^{-3}$ | $1.2441\times10^{-3}$ | $4.8414\times10^{-5}$ | $2.3692\times10^{-15}$ |
| $45^\circ$ | $64c$ flat | 127 | 2.6566372650 | $-9.5780\times10^{-5}$ | $+3.6748\times10^{-3}$ | $5.9160\times10^{-5}$ | $4.2327\times10^{-4}$ | $3.2340\times10^{-15}$ |
| $45^\circ$ | $63.55c$ short-transition flat | 127 | 2.6567063257 | $-2.6719\times10^{-5}$ | $+2.0082\times10^{-2}$ | $5.4814\times10^{-5}$ | $3.8441\times10^{-4}$ | $2.2656\times10^{-15}$ |
| $45^\circ$ | settled rolled wake | 80 | 2.4266505234 | $-2.3008\times10^{-1}$ | $+1.7648\times10^{-2}$ | 0.5427372 | $3.7555\times10^{-4}$ | $2.5347\times10^{-15}$ |

All shorter flat cases are retained in the two per-sequence CSVs. Every
single-shot relative residual is finite and below $10^{-10}$; all $q_f$,
$\mathcal Cq_f$, strengths, lifts, and probe fields are finite.

### Terminal strength iterations and oracles

| Angle | Geometry | Iterations | Final $C_L$ | Final strength defect | Final relative lift change | Oracle residual | Relative $\mu$ difference | Relative $C_L$ difference |
|---:|---|---:|---:|---:|---:|---:|---:|---:|
| $3.94^\circ$ | $64c$ flat | 11 | 0.2747427871 | $3.4250\times10^{-10}$ | $6.4345\times10^{-11}$ | $4.5285\times10^{-15}$ | $1.4341\times10^{-11}$ | $2.2038\times10^{-11}$ |
| $3.94^\circ$ | $63.55c$ short-transition flat | 31 | 0.2747424814 | $4.0837\times10^{-9}$ | $7.4370\times10^{-10}$ | $4.3161\times10^{-15}$ | $1.0377\times10^{-9}$ | $1.8028\times10^{-9}$ |
| $3.94^\circ$ | settled rolled wake | 47 | 0.2733869039 | $4.8396\times10^{-9}$ | $8.8209\times10^{-10}$ | $4.3165\times10^{-15}$ | $1.2223\times10^{-9}$ | $2.1276\times10^{-9}$ |
| $45^\circ$ | $64c$ flat | 12 | 2.6565915984 | $3.5904\times10^{-10}$ | $7.3001\times10^{-12}$ | $8.5665\times10^{-15}$ | $1.9055\times10^{-12}$ | $2.1840\times10^{-12}$ |
| $45^\circ$ | $63.55c$ short-transition flat | 36 | 2.6565895669 | $4.2446\times10^{-9}$ | $7.2084\times10^{-11}$ | $1.0740\times10^{-14}$ | $1.3888\times10^{-10}$ | $1.6861\times10^{-10}$ |
| $45^\circ$ | settled rolled wake | 49 | 2.4579146750 | $3.7610\times10^{-9}$ | $6.4313\times10^{-11}$ | $9.4515\times10^{-15}$ | $1.1486\times10^{-10}$ | $1.3830\times10^{-10}$ |

All six iterations converged at $\omega=1$; the adaptive fallback was not
needed. Inner residuals remained approximately $10^{-15}$. Geometry, row count,
options, nodes, and inactive strength storage stayed exact, and active strengths
changed only between inner solves.

Relative to Task 1, the iterated flat cases differ by about $8.0\times10^{-5}$
or less at $3.94^\circ$ and $5.4\times10^{-5}$ or less at $45^\circ$. The
iterated rolled wake differs by $-0.5013\%$ at $3.94^\circ$ and $-7.4836\%$ at
$45^\circ$.

## Reviewed Conclusions

The methods, generated measurements, invariant hashes, and stated equations
were reviewed together after both angle matrices completed. The direct
fixed-wake potential route removes almost all flat-wake discrepancy from the
semi-infinite reference, including the sensitivity to attached-transition
length seen in Task 2. This is consistent across single shots, fixed-strength
iterations, exterior probes, and the augmented oracle.

For the rolled geometry, direct potential also moves the result toward Task 1,
but it does not eliminate the discrepancy, especially at $45^\circ$. Because
the wake march had already settled, the fixed-geometry strength projection
converged, and the independent oracle agrees, the remaining difference is a
geometry/discretization effect rather than an unconverged strength handoff or
linear solve. Task 3 remains diagnostic-only; no production API or solver
behavior changed.
