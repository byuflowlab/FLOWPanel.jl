# wing_aileron

Solver benchmark study for a wing + aileron (control-surface) configuration.
Compares direct (`Backslash*`) and iterative/Krylov solvers on timing and
accuracy across a mesh-resolution ladder.

## Layout

```
wing_aileron/
├── meshes/     INPUT   — mesh ladder, paired wing_NN_NN.msh + surface_NN_NN.msh
├── results/    OUTPUT  — CSV written by the drivers
├── figures/    OUTPUT  — PNG written by plot/
└── plot/               — plotting scripts
```

The drivers live one level up, in `examples/`, not in this directory.

## Regenerating

All scripts use paths relative to the **repository root**, so run them from
there, not from inside this directory.

### 1. Produce data

```bash
julia --project examples/benchmarks.jl     # -> results/krylov_solvers.csv
julia --project examples/wing_aileron.jl   # -> results/coupled_timing_results.csv
```

`benchmarks.jl` selects the mesh resolution near line 120 — uncomment one
`run_names` pair at a time:

| Pair | Nodes/body | Panels/body |
|---|---|---|
| `wing_13_13` / `surface_13_13` | 330 | 656 |
| `wing_20_20` / `surface_20_20` | 834 | 1664 |
| `wing_25_25` / `surface_25_25` | 1242 | 2480 |
| `wing_40_40` / `surface_40_40` | 3274 | 6544 |
| `wing_60_60` / `surface_60_60` | 7314 | 14624 |

Note: the `25_25` pair is present but not referenced by any driver.

### 2. Produce figures

```bash
julia --project examples/wing_aileron/plot/plot_accuracy.jl
julia --project examples/wing_aileron/plot/plot_krylov_accuracy.jl
python  examples/wing_aileron/plot/plot_nps.py
```

| Script | Reads | Writes |
|---|---|---|
| `plot_accuracy.jl` | `results/krylov_solvers.csv` | `figures/accuracy_{coupled,iterative,mse}.png` |
| `plot_krylov_accuracy.jl` | `results/krylov_solvers.csv` | `figures/krylov_accuracy_{coupled,fmm}_high_nps.png` |
| `plot_nps.py` | `results/krylov_solvers.csv` | `figures/wing_aileron_benchmarks_nps_krylov.png` |
| `plot_2.jl` | `results/second.csv` | (displays only) |
| `plot_5.jl` | `results/times_bench.csv` | (displays only) |
| `plot_6.py` | `results/third.csv` | (displays only) |

## CSV schema

All result CSVs share one schema:

```
solver,nps,AOA,CL,CD,t_build,t_solve,total_time
```

`nps` is the total panel count across bodies. `solver` is one of
`BackslashCoupled`, `BackslashIterative`, `KrylovSolver`, `KrylovCoupled`.

## Known issues

- `examples/wing_aileron.jl` reads `wing.msh` / `surface.msh`, which do not
  exist in `meshes/`. The driver cannot run until those names are pointed at
  real meshes.
- `results/third.csv` is missing the newline after its header row; the header
  and first data row are concatenated on line 1, which `pandas` mis-parses.
- `plot_2.jl`, `plot_5.jl`, `plot_6.py` and their CSVs (`second.csv`,
  `third.csv`, `times_bench.csv`) carry non-descriptive sequential names.
- `results/{benchmarks,benchmarks_updated,runs_plot,times_bench(1)}.csv` have
  no reader script; they appear to be kept for ad-hoc analysis.
