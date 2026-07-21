# Task 2 — Finite-Wake Velocity-Through-Sources Verification

## Final Formulation

Theory gate: `docs/wake_solve_schemes.md`, “Discrete diagnostic construction:
Task 2 finite velocity-through-sources wake.”

The implemented finite-wake system is

$$
\mathcal G_\Delta\mu_V=-\mathcal S(\sigma_0+\sigma_f),
\qquad
\sigma_f=-u_f\cdot n,
\qquad
\mathcal G_\Delta=\mathcal B+\mathcal W_\Delta\mathcal C.
$$

- The watertight 6,688-panel Task 1 body, flow, Kutta map, and 16 exterior
  probes are reused exactly.
- The primary flat sequence uses an attached transition panel and fabricated
  free rows of length $0.5c$.
- The second flat sequence and the rolled production march use
  $|D_a|=0.05c$ while retaining $0.5c$ free-wake row spacing per time step.
- The attached transition strength is $\mathcal C\mu_V$; prescribed free rows
  use FLOWPanel's contiguous wake orientation and the same numerical sign.
- Flat free rows carry a fresh Task 1 circulation and remain prescribed.
- The production wake uses default induced-velocity convection at constant
  incidence with no pitching.
- No free-wake scalar potential is added to the body right-hand side.
- No production solver behavior is modified.
- Flat sequences use the terminal $10^{-3}$ relative-change criterion; a miss
  is retained as `nonconverged_at_64c` and does not prevent later cases from
  running. The rolled-wake 0.5% settling criteria and
  `nonconverged_at_80c` classification are unchanged.

Every reporting solve operates on a deep copy. Body nodes/connectivity and
attached geometry, wake nodes/strengths, active row count, and wake options are
asserted bit-exact afterward. Only the copied body's source/doublet values and
the transition circulation implied by its new doublets may change.

## Successful Commands

```bash
JULIA_NUM_THREADS=2 julia --project=. debug/dirichlet_solve/dirichlet_solve.jl task2-flat
JULIA_NUM_THREADS=2 julia --project=. debug/dirichlet_solve/dirichlet_solve.jl task2-flat-das005
JULIA_NUM_THREADS=2 julia --project=. debug/dirichlet_solve/dirichlet_solve.jl task2-march
JULIA_NUM_THREADS=2 julia --project=. debug/dirichlet_solve/dirichlet_solve.jl task1
JULIA_NUM_THREADS=2 julia --project=. debug/dirichlet_solve/dirichlet_solve.jl task2-march --alpha-deg=45 --output-dir=data/dirichlet_solve/alpha45
JULIA_NUM_THREADS=2 julia --project=. -e 'include("test/runtests_unit_solver.jl")'
JULIA_NUM_THREADS=2 julia --project=. -e 'include("test/runtests_unit_wake.jl")'
JULIA_NUM_THREADS=2 julia --project=. -e 'using Test; import FLOWPanel as pnl; using WriteVTK; include("test/test_helpers.jl"); include("test/runtests_unit_liftingbody.jl")'
JULIA_NUM_THREADS=2 julia --project=. -e 'include("test/runtests_unit_postprocess.jl")'
git diff --check
```

Results: solver 84/84; wake 82/82; lifting body 40/40; post-processing
405/405; diff check clean.

## Artifacts

Root: `data/dirichlet_solve/`

- `comparison.csv`
- `task2_flat_length_convergence.csv`
- `task2_flat_das005_length_convergence.csv`
- `task2_march_convergence.csv`
- `task2_flat_frozen_state.toml`
- `task2_flat_das005_frozen_state.toml`
- `task2_march_frozen_state.toml`
- `task2.config.toml`
- `task2.metadata.toml`
- `dirichlet_task2_flat_terminal_body.pvd`
- `dirichlet_task2_flat_terminal_wake.pvd`
- `dirichlet_task2_flat_das005_terminal_body.pvd`
- `dirichlet_task2_flat_das005_terminal_wake.pvd`
- `dirichlet_task2_march_terminal_body.pvd`
- `dirichlet_task2_march_terminal_wake.pvd`

The 45-degree matrix uses the same artifact layout under the isolated root
`data/dirichlet_solve/alpha45/`; the corrected run refreshed its baseline and
rolled-march artifacts while retaining the existing flat-sequence artifacts.

## Measurements

Task 1 was recomputed after Task 2 and remained
$C_{L,\infty}=0.2747643938323718$ with relative direct residual
$2.4807313750\times10^{-15}$. Its comparison row was preserved.

### Fabricated flat wake

| $L/c$ | Rows | $C_L$ | $C_L-C_{L,\infty}$ | Relative difference |
|---:|---:|---:|---:|---:|
| 1 | 1 | 0.2432101847 | -0.0315542091 | -0.114840969 |
| 2 | 3 | 0.2628588527 | -0.0119055412 | -0.043330000 |
| 4 | 7 | 0.2707134454 | -0.0040509484 | -0.014743353 |
| 8 | 15 | 0.2732656785 | -0.0014987154 | -0.005454547 |
| 16 | 31 | 0.2739802786 | -0.0007841153 | -0.002853773 |
| 32 | 63 | 0.2741666448 | -0.0005977490 | -0.002175497 |
| 64 | 127 | 0.2742139824 | -0.0005504114 | -0.002003212 |

The final relative length-to-length change is
$1.7263004809\times10^{-4}\le10^{-3}$, so the sequence is length-converged.
At $64c$, the direct relative residual is $2.5242451747\times10^{-15}$,
transition circulation L2 norm is 15.6930365837, transition/free-row relative
mismatch is 0.0015672509, and the exterior-probe relative L2 difference is
$3.8412533891\times10^{-5}$. All seven frozen states were preserved exactly.

### Fabricated flat wake with $|D_a|=0.05c$

The prescribed free rows remain $0.5c$ long.  The nominal labels retain the
same row counts as the primary sequence; the table reports the actual total
trailing-edge-to-closure extent.

| Nominal $L/c$ | Total extent $/c$ | Rows | $C_L$ | $C_L-C_{L,\infty}$ | Relative difference |
|---:|---:|---:|---:|---:|---:|
| 1 | 0.55 | 1 | 0.2460108956 | -0.0287534983 | -0.104647833 |
| 2 | 1.55 | 3 | 0.2650115934 | -0.0097528004 | -0.035495139 |
| 4 | 3.55 | 7 | 0.2699093080 | -0.0048550858 | -0.017669996 |
| 8 | 7.55 | 15 | 0.2711682416 | -0.0035961522 | -0.013088130 |
| 16 | 15.55 | 31 | 0.2714745256 | -0.0032898683 | -0.011973416 |
| 32 | 31.55 | 63 | 0.2715482839 | -0.0032161099 | -0.011704973 |
| 64 | 63.55 | 127 | 0.2715662450 | -0.0031981489 | -0.011639604 |

The final relative length-to-length change is
$6.6138636288\times10^{-5}\le10^{-3}$, so this sequence is also
length-converged. At the $63.55c$ terminal extent, the direct relative residual
is $2.4264865024\times10^{-15}$, transition circulation L2 norm is
15.5981876606, transition/free-row relative mismatch is 0.0076072780, and the
exterior-probe relative L2 difference is $5.0010898823\times10^{-5}$. All
seven frozen states were preserved exactly.

Relative to the primary $|D_a|=0.5c$ terminal result, shortening only the
attached transition reduces $C_L$ by 0.0026477374, or 0.9656% of the primary
flat value.

### Rolled-up production wake

The production history satisfied both settling criteria at the first allowed
checkpoint, $40c/U_\infty$ (step 80, 80 free rows), with
$|D_a|=0.05c$:

- relative peak-to-peak $C_L$: $2.1545232371\times10^{-4}$;
- normalized linear drift: $2.1138990687\times10^{-4}$;
- frozen $C_L=0.2653120608686893$;
- signed Task 1 difference: $-0.0094523329636825$;
- relative Task 1 difference: $-0.0344015934227970$;
- direct relative residual: $2.3877097022\times10^{-15}$;
- transition circulation L2 norm: 15.2569256732;
- transition/free-row relative mismatch: $6.5082443350\times10^{-6}$;
- exterior-probe relative L2 difference: $1.2523524566\times10^{-3}$.

The former $|D_a|=0.5c$ march gave $C_L=0.2738709687$, or $-0.3252\%$
relative to Task 1.  Correcting the march to $0.05c$ lowers $C_L$ by
$3.1252\%$ relative to that result and increases the semi-infinite discrepancy
to $-3.4402\%$. The settling metrics rule out an unconverged history as the
source of that increase.

Every marching solve preserved its pre-solve body geometry, wake nodes,
strengths, row count, and options. The terminal reporting solve preserved the
frozen state exactly; hashes and options are recorded in the invariant TOML.

## 45-Degree Angle-of-Attack Sensitivity

The complete matrix was repeated at $\alpha=45^\circ$ on the same 6,688-panel
geometry. The new semi-infinite reference is
$C_{L,\infty}=2.6567330445$ with direct relative residual
$1.5562174039\times10^{-14}$. The fixed 3.94-degree lift regression value is
not used at this angle; finite-value, residual, Kutta-map, and frozen-state
assertions remain active.

| Case | Status / terminal metric | Terminal $C_L$ | $C_L-C_{L,\infty}$ | Relative difference | Direct relative residual | Transition/free mismatch | Exterior-probe difference |
|---|---|---:|---:|---:|---:|---:|---:|
| Flat, $|D_a|=0.5c$ | converged at $64c$; final change $1.12284\times10^{-4}$ | 2.6529624339 | -0.0037706106 | -0.001419266 | $4.47419\times10^{-15}$ | 0.001218832 | 0.000452204 |
| Flat, $|D_a|=0.05c$ | converged at $63.55c$; final change $3.60211\times10^{-5}$ | 2.6366242601 | -0.0201087844 | -0.007568989 | $3.48682\times10^{-15}$ | 0.006078884 | 0.000385852 |
| Rolled production wake, $|D_a|=0.05c$ | settled at $40c$; peak-to-peak $2.35957\times10^{-4}$, drift $2.26088\times10^{-4}$ | 2.4090022208 | -0.2477308238 | -0.093246412 | $3.29574\times10^{-15}$ | $6.48021\times10^{-6}$ | 0.542810665 |

Both flat sequences therefore converge under the unchanged $10^{-3}$ terminal
relative-change criterion, and the march settles at its first eligible
checkpoint under both unchanged 0.5% criteria. All seven states in each flat
sequence, every marching solve, and the terminal march reporting state preserve
the recorded frozen-state invariants exactly. Shortening only the transition
panel reduces the terminal flat $C_L$ by 0.6158% relative to the primary flat
result at 45 degrees.

Compared with 3.94 degrees, the magnitudes of the terminal semi-infinite lift
discrepancies are 29.2% smaller for the primary flat wake and 35.0% smaller for
the short-transition flat wake. In contrast, the corrected rolled-wake
discrepancy is 2.71 times as large: 9.325% rather than 3.440%. The corrected
$45^\circ$ result is 7.51% lower in $C_L$ than the former $|D_a|=0.5c$ march,
whose semi-infinite discrepancy was only 1.961%. The transition/free-row
mismatch remains essentially zero at both angles, while the exterior-probe
difference grows from 0.00125 to 0.543. Thus settling and a continuous
circulation handoff do not imply agreement with the semi-infinite exterior
field, especially at high incidence.

## Conclusion

The primary fabricated finite composite approaches a distinct converged value
very near, but not equal to, Task 1. Its residual transition/free-row mismatch
decreases strongly with length, while its remaining signed lift difference is
about $-0.20\%$. With the same half-chord free rows but $|D_a|=0.05c$, the
converged lift difference grows to about $-1.16\%$ and the terminal circulation
handoff mismatch grows to about $0.76\%$. With the rolled march corrected to
$|D_a|=0.05c$, it settles by $40c/U_\infty$ but differs from Task 1 by
$-3.44\%$ at $3.94^\circ$ and $-9.32\%$ at $45^\circ$, despite essentially
continuous transition/free-row circulation handoff. The earlier near-agreement
of $-0.33\%$ used $|D_a|=0.5c$ and was not representative of the production
pitching-wing configuration. The rolled-wake geometry, attached-panel length,
and current velocity-through-sources representation remain measured mechanisms
for later tasks; no production solver behavior changed.
