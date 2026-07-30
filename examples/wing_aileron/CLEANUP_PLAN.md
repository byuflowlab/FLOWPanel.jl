# wing_aileron — Cleanup Plan

Executed 2026-07-30. **Status: complete.**

15 files were initially staged into `_to_delete/`. You then rescued 7 figures
back into `figures/`, and approved deletion of the remaining 8. Those 8 are now
removed (staged deletions, not yet committed).

All 8 remain recoverable from git history:

```bash
git checkout dc07f52 -- examples/wing_aileron/<name>   # surface_first, surface_second, wing_first, wing_double
git checkout 7e99049 -- examples/wing_aileron/<name>   # surf501000.msh, wing3.msh
git checkout 4422ea7 -- examples/wing_aileron/<name>   # plotting.jl, benchmarks_11196_panels.png
```

The sections below record why each file was flagged. Entries marked
**[KEPT]** were rescued and are still present in `figures/`.

## How files were assessed

Filesystem mtimes are all identical (checkout artifact), so staleness came from
git history in the `FLOWPanel.jl` submodule plus a repo-wide reference scan.

The whole directory arrived in three commits:

| Commit | Date | Description |
|---|---|---|
| `dc07f52` | 2026-07-08 | plotting and stuff |
| `7e99049` | 2026-07-09 | Clean up |
| `4422ea7` | 2026-07-29 | Massive timing and accuracy plotting save |

`4422ea7` was a bulk save that touched nearly everything, so "last modified
Jul 29" alone is not evidence of use. Each candidate was checked for (a) a real
content diff in `4422ea7`, (b) an inbound reference from any `.jl`/`.py`/`.md`/
`.tex` file in the repo, and (c) byte-identity with another file.

## Producers live outside this directory

This directory is pure data — the generating scripts are siblings:

- `examples/benchmarks.jl` → writes `krylov_solvers.csv`, reads the
  `wing_NN_NN.msh` / `surface_NN_NN.msh` ladder
- `examples/wing_aileron.jl` → writes `coupled_timing_results.csv`

---

## Files flagged (15 originally staged, 2.9 MB)

### Stale meshes — 6 files, ~2.5 MB (all deleted)

Last touched Jul 8–9, never re-edited, referenced by **no** script in the repo.

| File | Size | Nodes | Note |
|---|---|---|---|
| `surf501000.msh` | 961 KB | 10498 | one-off high-res pair with `wing3.msh` |
| `wing3.msh` | 933 KB | 10194 | one-off high-res |
| `surface_second.msh` | 333 KB | 3778 | orphan half-pair — partner `wing_third.msh` deleted in `4422ea7` |
| `surface_first.msh` | 98 KB | 1194 | superseded by the `NN_NN` ladder |
| `wing_first.msh` | 97 KB | 1194 | superseded by the `NN_NN` ladder |
| `wing_double.msh` | 97 KB | 1194 | **byte-identical** to `wing_first.msh` |

### Dead scripts — 1 file (deleted)

| File | Reason |
|---|---|
| `plotting.jl` | Reads `wing_aileron/timing_summary.csv`, **deleted** in `4422ea7`. Cannot run as written. |

### Orphan figures — 8 files (1 deleted, 7 kept)

No surviving script produces these, and none is referenced by
`code/scripts/reports/*.md` or `code/scripts/validation/validation.tex`. All
regenerable.

| File | Reason |
|---|---|
| `cl_aoa_bc_b.png` **[KEPT]** | Jul 9, no producer |
| `nps_aoa_cl.png` **[KEPT]** | Jul 9, no producer |
| `benchmarks_11196_panels.png` | no producer |
| `krylov_accuracy_coupled.png` **[KEPT]** | superseded — `plot_krylov_accuracy.jl` now writes only `_high_nps` variants |
| `krylov_accuracy_fmm.png` **[KEPT]** | superseded, same reason |
| `wing_aileron_benchmarks_nps.png` **[KEPT]** | no producer — `plot_nps.py` now writes only the `_krylov` variant |
| `benchmarks_plots_2500-38600.png` **[KEPT]** | **byte-identical** to `wing_aileron_benchmarks_nps.png` |
| `krylov_timing_2026-07-29_1239.png` **[KEPT]** | **byte-identical** to `wing_aileron_benchmarks_nps_krylov.png`, which is kept |

Verified: no surviving `.jl`/`.py`/`.md`/`.tex` file references any moved file.

---

## Left in place, but worth your judgment (not moved)

Unreferenced but too recent to call stale — possibly staged deliberately.

| File(s) | Situation |
|---|---|
| `wing_25_25.msh`, `surface_25_25.msh` | Added Jul 29. The string `25_25` appears in **no** script — `benchmarks.jl` sweeps 13/20/40/60 only. Staged for an upcoming sweep, or leftover. |
| `benchmarks.csv`, `benchmarks_updated.csv`, `times_bench(1).csv` | No committed reader. `times_bench(1).csv` has a download-collision filename. |
| `runs_plot.csv` | No committed reader, but gained 16 rows in `4422ea7` — actively appended to. |
| `plot_2.jl`+`second.csv`, `plot_5.jl`+`times_bench.csv`, `plot_6.py`+`third.csv` | Each pair self-consistent and runnable, and **all six got real content diffs on Jul 29**. Active despite throwaway names. Rename, don't delete. |

### Two live path bugs found along the way

Independent of cleanup — these scripts point at meshes that do not exist:

1. `deps/FLOWPanel.jl/examples/wing_aileron.jl:110,130` reads `wing.msh` and
   `surface.msh` — neither has ever existed in this directory.
2. `code/scripts/validation/solvers/fmm_solver/wing_aileron.jl:153` reads
   `nasa_wing.msh` and `nasa_surface_spaced_repaired.msh` — both deleted on
   2026-05-08 in commit `1180b86`. This file is currently modified in your
   working tree, so this one probably matters now.

### Data hygiene

`third.csv` is missing the newline after its header row — header and first data
row are concatenated on line 1. `pandas` will silently mis-parse it.

---

## Proposed folder architecture

The directory currently mixes four kinds of thing at one level: input meshes,
output CSVs, output PNGs, and plotting scripts. Splitting by role makes the
pipeline legible (`meshes → results → figures`, driven by `plot/`):

```
examples/wing_aileron/
├── README.md                 # what this study is, how to regenerate it
├── meshes/                   # INPUT — mesh ladder
│   ├── wing_13_13.msh        surface_13_13.msh
│   ├── wing_20_20.msh        surface_20_20.msh
│   ├── wing_40_40.msh        surface_40_40.msh
│   └── wing_60_60.msh        surface_60_60.msh
├── results/                  # OUTPUT of examples/benchmarks.jl
│   └── krylov_solvers.csv
├── figures/                  # OUTPUT of plot/
│   ├── accuracy_coupled.png
│   ├── accuracy_iterative.png
│   ├── accuracy_mse.png
│   ├── krylov_accuracy_coupled_high_nps.png
│   ├── krylov_accuracy_fmm_high_nps.png
│   └── wing_aileron_benchmarks_nps_krylov.png
└── plot/                     # plotting scripts
    ├── plot_accuracy.jl
    ├── plot_krylov_accuracy.jl
    └── plot_nps.py
```

### Cost of adopting it

Every path in the plot scripts is hardcoded relative to the repo root
(`"examples/wing_aileron/krylov_solvers.csv"`), so the moves require edits:

- `benchmarks.jl:120-123, 167` — mesh ladder + output path
- `plot_accuracy.jl:8, 76-77, 101`
- `plot_krylov_accuracy.jl:8, 95-96`
- `plot_nps.py:8, 130`

~12 string edits. Better to define the directory once at the top of each script
than to re-hardcode longer paths.

### One open decision: where meshes live

The repo convention is a shared `examples/data/` store — `cessna.msh`,
`wing_ar4_naca0016_5.msh`, the `dji9443_*` family all sit there. `wing_aileron`
is the only example carrying its own meshes.

- **Option A (shown above): local `meshes/`.** Self-contained and portable;
  smallest diff.
- **Option B: `examples/data/wing_aileron/`.** Matches repo convention, all
  meshes discoverable in one place; larger diff, splits the study across trees.

I lean **A** — this is a self-contained benchmark study rather than a shared
asset, and no other example reads these meshes. B is more consistent if you'd
rather not maintain two conventions.

### Also worth doing

- **Rename the `plot_N` scripts.** `plot_2.jl`, `plot_5.jl`, `plot_6.py` are
  active but the names say nothing. Name them for what they plot; same for
  `second.csv`, `third.csv`, `times_bench.csv`.
- **Consider `.gitignore` for `figures/`.** All six kept PNGs are regenerable
  from `results/` + `plot/`. Committing them is what produced the duplicate and
  orphan figures cleaned up above.

---

## Outcome

Deleted (8): `surf501000.msh`, `wing3.msh`, `surface_second.msh`,
`surface_first.msh`, `wing_first.msh`, `wing_double.msh`, `plotting.jl`,
`benchmarks_11196_panels.png` — 2.5 MB recovered.

Kept (7): the figures marked **[KEPT]** above, now in `figures/`.

Note: two of the kept figures are byte-identical duplicates of others in
`figures/` — `benchmarks_plots_2500-38600.png` = `wing_aileron_benchmarks_nps.png`,
and `krylov_timing_2026-07-29_1239.png` = `wing_aileron_benchmarks_nps_krylov.png`.
