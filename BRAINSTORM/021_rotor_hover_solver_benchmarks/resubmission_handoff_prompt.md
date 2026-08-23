> **SUPERSEDED 2026-08-22 — EXECUTED. Do not run this document again.**
> Both cluster checkouts were synced (whole-tree rsync of `src/` + `benchmark/`,
> md5-verified) and all 8 jobs were resubmitted: 13306549-52 (p2b R1-R4),
> 13306553 (p1-tune-R7), 13306554 (p1-table-R6-multi), 13306556/57 (R6 singles).
> See the `2026-08-22 (later)` entry in `log.md` for provenance and the
> correction that the "never submitted" jobs were actually dependency-cancelled.

# Handoff — sync the 021 fixes to the cluster and resubmit (written 2026-08-22)

## Your task

Get the BRAINSTORM 021 campaign running again. Nothing is submitted; the
previous session diagnosed and fixed the causes but stopped before syncing,
because doing so blind risked mixing source generations.

## State you are inheriting

**Local repo** (`~/Dropbox/research/projects/FLOWPanel.jl`, branch `fastmultipole`,
HEAD `d6bf8b6`): contains UNCOMMITTED fixes for two of the three failure classes,
verified locally (347 unit tests pass; R1 Phase-2b driver smoke passes all three
sections). The working tree is ALSO dirty with unrelated 018/022 work — do not
assume every modified file belongs to 021.

Files carrying the 021 fixes:
- `src/FLOWPanel_fmm.jl`, `src/FLOWPanel_instrumentation.jl`,
  `src/FLOWPanel_solver.jl` — thread `nearfield_cache_max_bytes` /
  `nearfield_cache_max_build_time` from `KrylovSolver` down to
  `FastMultipole.build_nearfield_cache!`. These caps were previously
  UNREACHABLE from FLOWPanel, which hard-blocked R4 (needed 4.47 GiB vs a
  4 GiB default). `FLOWPanel_solver.jl` also restructures the Barba ILU
  `max_pattern_entries` guard to report the TOTAL required, not the running
  subtotal where it tripped.
- `benchmark/rotor_hover_solver_phase2b_nearfield_cache.jl` — caps hoisted to
  the top, applied to all three sections, defaults 32 GiB / 1800 s, overridable
  via `TUNE_MAX_BYTES` / `TUNE_MAX_BUILD_TIME`.
- `benchmark/slurm/p2b_nearfield.sh` — header documents the new caps (no
  behavioural change).
- Seven `benchmark/rotor_hover_solver_*.jl` — ILU `max_pattern_entries`
  2048*N -> 8192*N.

**Cluster** (`ssh orc`; Slurm bins at `/apps/slurm/latest/bin`; use `bash -lc`
for a login shell so `module load julia` works — do NOT wrap heredocs in a
quoted `bash -lc`):
- `~/projects/p2b/FLOWPanel.jl` — plain rsync copy, NOT a git repo. Ryan fixed
  its `Project.toml`/`Manifest.toml` on 08-22 16:56 (GeoIO and PythonPlot removed
  outright; PythonPlot now lives in the root `@v1.12` env). `load_check2.log`
  says `LOAD_OK_WITH_CACHE`. Its `src/` is from Aug 19 and does NOT have the
  fixes above.
- `~/projects/FLOWPanel.jl` — git repo at commit `5615ada`, FOUR commits behind
  local, dirty (`Project.toml`, `src/FLOWPanel.jl`, several examples). STILL has
  GeoIO + PythonPlot as hard `[deps]`, and `import FLOWPanel` STILL FAILS in
  PythonPlot/PythonCall init. This is where the Phase-1 jobs run, so it must be
  fixed the same way p2b was before anything is submitted.
- 8 x BRAINSTORM 018 jobs running (~2 days elapsed), MUST NOT be disturbed.
  Disk healthy: 131 G / 2.0 T.

## Steps

1. **Ask Ryan how to sync**, do not guess. The question is whether to commit the
   021 fixes and `git pull` on the main checkout (clean, but the checkout is
   dirty and 4 commits behind), or rsync only the named files (surgical, but
   risks mixing Aug-19 and Aug-22 source generations — the exact failure that
   already cost this campaign a run). Confirm you may touch
   `~/projects/FLOWPanel.jl` given the 018 jobs share the machine but not that
   checkout.
2. **Fix the main checkout's dependencies** the way p2b was fixed: remove GeoIO
   and PythonPlot from `Project.toml`, resolve, and verify
   `julia --project -e 'import FLOWPanel; println("OK")'` prints OK. Without
   this, every Phase-1 job dies at load in ~7 minutes.
3. **Verify the fixes actually landed** before submitting — on BOTH checkouts:
   `grep -c nearfield_cache_max_bytes src/FLOWPanel_solver.jl` (expect > 0) and
   `grep -c 'max_pattern_entries=8192' benchmark/rotor_hover_solver_phase1_table.jl`.
   Record the src md5 / commit in the run notes (021 has been burned before by
   assuming the cluster source matched local).
4. **Resubmit**, all four p2b rungs together so the ledger stays internally
   consistent (R1/R2 rows measured under the OLD 60 s cap would not be
   comparable to R3/R4 under the new one):
   - p2b R1-R4: see the `sbatch` line in `benchmark/slurm/p2b_nearfield.sh`
     (zen3, exclusive, `--export=ALL,RUNG=Rn,MODE=multi,TUNE_MACS=...`).
     **Investigate why 13242660/62/64/66 never entered Slurm history** — the
     afterok chain looks malformed and silently halved the campaign.
     Consider raising `--mem` above 64G if a rung's cache approaches the cap.
   - `p1-tune-R7` (previously 13206080; failed on GeoIO at 7 min, NOT OOM —
     MaxRSS 1.0 GiB vs ReqMem 256 GiB, so resources were fine).
   - `p1-table-R6-multi` (previously 13206077; needs the 8192*N ILU fix).
   - Check whether any `p1-table-R7-*` was ever submitted — none was as of 08-22.

## Context you should not re-derive

- The FGS unsteady staleness bug (stale FMM tree under rigid rotation, 21% error
  at 80 deg) is CONFIRMED and FIXED as of 08-21, both the branch-center staleness
  and the `simulate!` lifecycle-ordering bug behind the residual 2e-3 plateau.
  018/022 are NOT contaminated — they use dense `pnl.Backslash`.
  `fgs_wake_plateau_handoff_prompt.md` is SUPERSEDED; do not execute it.
- The ILU preconditioner scales ~N^1.5 (487 -> 3691 entries/row R1->R7; apply is
  linear in nnz at ~0.6 ns/nonzero, so ~1 s/apply and ~23 GiB of pattern at R7).
  Measured and written up in `ledger.md` + `ilu_scaling/`. Do not re-measure.
- `log.md` (newest first) is the reliable narrative. The item file's
  `## Current status` now carries a STALE banner pointing there.
- The near-field cache build is SERIAL (only `nearfield_matvec!` is threaded), so
  `max_build_time` is wall-clock, not a per-thread budget. `tune_fmm`'s clamp
  path is dead code when the FIRST trial is infeasible
  (`FastMultipole/src/autotune.jl:109-121`), and its MAC sweep starts at its most
  expensive point (0.3).

## Open decisions for Ryan

- Whether to harden `plot_discretization::Bool = true`
  (`src/liftingline/liftingline.jl:1511`): with PythonPlot demoted to a weakdep
  this throws `UndefVarError: plt` if PythonPlot was never loaded. No in-repo
  caller is affected (all three lifting-line examples import PythonPlot and set
  the flag), so it is latent, not live.
- Whether the root `CondaPkg.toml` should stay now that nothing loads CondaPkg.
- `CLAUDE.md` still says meshes are built with GeometricTools "aliased gt
  everywhere"; GeometricTools is no longer in the Manifest at all.
