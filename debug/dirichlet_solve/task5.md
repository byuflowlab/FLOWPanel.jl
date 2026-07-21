# Task 5 — Line-integral trace correction

## Objective and formulation

Evaluate only

```julia
pnl.TraceCorrected(estimator=:line_integral, ...)
```

on the three frozen Task 2 wake states at `α=3.94°` and `45°`, using both a
single-shot route and the Task 3/4 fixed-geometry active-strength projection.
The primary ranking metric is

```text
norm(c - c_oracle) / norm(c_oracle),  c_oracle = C*q_f,
```

where `q_f` is evaluated directly from the same frozen free-wake geometry and
strengths only after the eligible velocity-only solve and physical diagnostics.

The three variants are:

- `trap`: straight lower-to-upper path, trapezoid quadrature;
- `simpson`: straight lower-to-upper path, Simpson quadrature;
- `deformed`: two-segment interior path, Simpson quadrature, `path_depth=1`.

Every independent solve initializes a fresh `TraceCorrectedState` and runs
through `_steady_aerodynamics!` with `DirectBackend` for wake, body solve, and
physical-system evaluation. The solved system is

```text
G*μ̃ + S*(σ0+σ) - W*c = 0,
```

and every downstream consumer uses physical circulation
`γ = C*μ̃ - c` through `_get_wakestrength_Gamma`, including pressure, exterior
probes, transition/free-row comparisons, and VTK output.

The fixed-geometry route projects all active wake rows toward this physical
`γ`, with strength defect and relative lift change both at most `1e-8` for
three consecutive iterations, a 200-iteration cap, and adaptive relaxation no
lower than `1/16`.

## Theory audit

The derivation in `docs/wake_solve_schemes.md` matches
`src/FLOWPanel_formulation.jl` and the affine wake-strength implementation:

- lower-to-upper integration gives `c = q_upper - q_lower`;
- the solve adds `+W*c` to its right-hand side;
- physical circulation is `C*μ̃-c`;
- all estimator probes request wake velocity only.

No theory correction was required. The direct wake scalar-potential oracle is
diagnostic and remains outside the eligible solve path.

## Successful commands

```bash
julia --project test/formulation_test.jl

JULIA_NUM_THREADS=2 julia --project \
  debug/dirichlet_solve/dirichlet_solve.jl \
  --task task5-<trap|simpson|deformed>-<flat|flat-das005|march>-<single|iterated> \
  --alpha-deg 3.94 --output-dir data/dirichlet_solve

JULIA_NUM_THREADS=2 julia --project \
  debug/dirichlet_solve/dirichlet_solve.jl \
  --task task5-<trap|simpson|deformed>-<flat|flat-das005|march>-<single|iterated> \
  --alpha-deg 45 --output-dir data/dirichlet_solve/alpha45

JULIA_NUM_THREADS=2 julia --project \
  debug/dirichlet_solve/dirichlet_solve.jl \
  --task task5-depth-march --alpha-deg 3.94

julia --project debug/dirichlet_solve/dirichlet_solve.jl --task task5-audit
```

The formulation test passed all seven stages: `7/7`, `5/5`, `864/864`,
`4/4`, `3/3`, `8/8`, and `14/14`.

## Primary rolled-wake results and ranking

Relative trace error against the direct oracle:

| Rank | Variant | 3.94° single | 3.94° iterated | 45° single | 45° iterated | Worst |
|---:|---|---:|---:|---:|---:|---:|
| 1 | straight Simpson | 4.2872e-11 | 4.2468e-11 | 2.3062e-11 | 2.3645e-11 | **4.2872e-11** |
| 2 | deformed Simpson | 4.2881e-11 | 4.2568e-11 | 2.3038e-11 | 2.3710e-11 | **4.2881e-11** |
| 3 | straight trapezoid | 1.4811e-7 | 1.4811e-7 | 1.5284e-7 | 1.5283e-7 | **1.5284e-7** |

Straight Simpson wins by the prescribed worst-case criterion, although its
difference from deformed Simpson is negligible. Trapezoid is approximately
`3.6e3` times worse on the worst primary cell.

Rolled-wake lift context for straight Simpson:

| Angle | Route | `CL` | `f_oracle` | Iterations |
|---:|---|---:|---:|---:|
| 3.94° | single | 0.2661161139 | 0.259991 | 0 |
| 3.94° | iterated | 0.2683602493 | 0.377492 | 44 |
| 45° | single | 2.4130610148 | 0.229982 | 0 |
| 45° | iterated | 2.4268700116 | 0.365301 | 46 |

Quadrature changes the trace metric strongly but changes these lift values only
at approximately `1e-9` or less. Therefore the remaining Task 3 lift gap is not
line-integral estimator error. In particular, this almost exact direct trace
closes less of the rolled lift gap than Task 4 Green reconstruction, so Task 4's
larger lift movement was not evidence of a more accurate trailing-edge trace.

The flat robustness cells lead to the same ranking. The short-transition flat
fixed-point route overshoots below the old lift for every variant:
`f_oracle=-0.5973` at 3.94° and `-0.8451` at 45°. Because Simpson and deformed
Simpson have `c` errors of only about `2e-12–5e-12` there, this behavior is a
fixed-geometry coupling result rather than a trace-quadrature failure.

## Deformed-path depth sensitivity

The rolled 3.94° deformed/Simpson results are insensitive to the tested depth:

| `path_depth` | Single trace error | Iterated trace error | Single `CL` | Iterated `CL` |
|---:|---:|---:|---:|---:|
| 0.5 | 4.2867e-11 | 4.2394e-11 | 0.2661161139 | 0.2683602493 |
| 1.0 | 4.2881e-11 | 4.2568e-11 | 0.2661161139 | 0.2683602493 |
| 2.0 | 4.3184e-11 | 4.2393e-11 | 0.2661161139 | 0.2683602493 |

The actual waypoint displacement is `path_depth` times the upper/lower
trailing-edge control-point separation, along the upstream attached-wake
direction. Path spreads between straight and deformed Simpson are at most
`2.4e-13` relative at 3.94° and `6.4e-14` at 45°.

## Acceptance gates and artifact audit

All 36 primary cells and all six depth cells passed.

- Worst primary formulation residual: `4.9635e-15`, below `1e-10`.
- Worst corrected-circulation consistency residual: `9.59e-17`, below `1e-10`.
- All fields and recorded timing values are finite.
- Every independent solve used fresh formulation state.
- Body nodes, cells, shedding maps, attached-wake geometry, wake nodes,
  storage shapes, row count, active/inactive strengths, and wake options were
  exact within each frozen solve. Only solved body strengths and the expected
  affine `c` fields changed; active free-wake strengths changed only between
  outer iterations.
- The eligible path requested no free-wake scalar potential. The oracle was
  evaluated afterward on the same frozen state.
- All outer iterations converged at `ω=1`: flat `13/14`, short-transition flat
  `45/50`, and rolled `44/46` iterations at 3.94°/45°.
- Stored Task 1–4 comparison rows, Task 2 frozen-state records, Task 3/4
  histories and invariants, Task 5 raw route/trace/history CSVs, snapshot TOMLs,
  backend settings, and terminal body/wake VTK artifacts were included in the
  conclusion audit.

Principal artifacts in each angle directory:

- additive `comparison.csv`;
- `task5_<variant>_<case>_<route>_{route,iteration,trace}.csv`;
- `task5_<variant>_<case>_<route>_frozen_state.toml`;
- `task5.config.toml`, `task5.metadata.toml`, `task5_invariants.toml`;
- terminal `dirichlet_task5_*_{body,wake}.pvd` outputs;
- `task5_line_integral_spreads.csv`.

Cross-angle and sensitivity summaries are
`task5_line_integral_primary_cells.csv`,
`task5_line_integral_ranking.csv`, and
`task5_line_integral_path_depth.csv`.

## Audited conclusion

Use the straight Simpson line integral as the Task 5 estimator: it has the
lowest worst-case primary trace error, no meaningful path sensitivity, and the
simplest path. Deformed Simpson is numerically tied and serves as a strong
robustness check. Trapezoid remains accurate in absolute terms but is clearly
third by the requested trace metric. The estimator has therefore been isolated
and validated; no self-consistent wake march was performed.
