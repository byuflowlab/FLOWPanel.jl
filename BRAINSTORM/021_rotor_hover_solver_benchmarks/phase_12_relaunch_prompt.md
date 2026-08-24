# Handoff — BRAINSTORM 021: diagnose the failed HPC jobs, repair, relaunch Phases 1–2

Written 2026-08-24 at a context reset. **Ryan's report: all HPC jobs failed.**

## Your task, in order

1. **Ascertain why they failed.** Nobody has looked at the cluster yet — the
   failure is reported, not diagnosed. Do not assume a cause from the catalogue
   below; use it as a checklist to test against evidence.
2. **Repair the cause**, in the repo where it belongs (launcher, driver, or
   cluster state), and say plainly which of the ten jobs each fix explains.
3. **Relaunch every job Phases 1 and 2 still need.** Get Ryan's go before
   `sbatch` — job submission is inline work and his call, never a subagent's.

Report a compact table: job → exit state → cause → fix → resubmitted id.

## Access

`ssh orc` needs a **live ControlMaster socket**; without one it prompts for 2FA
and will hang. If the socket is dead, ask Ryan to run `! ssh orc` in the session
so the output lands in the conversation and the socket comes up. Do not try to
authenticate yourself.

Read `agent_policies/HPC.md` before touching anything. On the login node,
`/apps/instructions_for_ai_agents/BYU_ORC_AGENTS.md` is the site policy. Slurm
binaries live in `/apps/slurm/latest/bin` (not always on `PATH`).

**Delegation:** route job status, log tailing, and disk reporting to the
`hpc-monitor` subagent, and any VTK sweeping to `hpc-storage`. Both are
read-only/delete-only by design. Keep diagnosis reasoning, repairs, and
submission inline.

## The ten jobs that failed

Submitted 2026-08-22, all were `PENDING` as of 2026-08-22 night. Cluster was at
128/136 nodes.

| Job | id | Launcher | Driver | Walltime |
| --- | --- | --- | --- | --- |
| p2b-nearfield-R1 | 13306549 | `p2b_nearfield.sh` | `rotor_hover_solver_phase2b_nearfield_cache.jl` | 24 h |
| p2b-nearfield-R2 | 13306550 | " | " | 24 h |
| p2b-nearfield-R3 | 13306599 | " | " | 72 h |
| p2b-nearfield-R4 | 13306600 | " | " | 72 h |
| p1-table-R6-multi | 13306554 | `p1_table.sh` | `rotor_hover_solver_phase1_table.jl` | 24 h |
| p1-table-R6-single-a | 13306556 | " | " | 72 h |
| p1-table-R6-single-jacobi | 13306557 (`afterany` on 556) | " | " | 72 h |
| p1-tune-R7-s1 | 13306603 | `p1_tune.sh` `STAGES=1` | `rotor_hover_solver_phase1_tune_hpc.jl` | 72 h |
| p1-tune-R7-s2 | 13306604 (`afterok` on 603) | `STAGES=2` | `rotor_hover_solver_phase1_fgstune.jl` | 72 h |
| p1-tune-R7-s3 | 13306605 (`afterok` on 604) | `STAGES=3` | `rotor_hover_solver_phase1_fgsprecond.jl` | 72 h |

Two cluster checkouts are in play: `~/projects/FLOWPanel.jl` (main; also hosts 8
BRAINSTORM-018 jobs — **do not disturb those**) and `~/projects/p2b/FLOWPanel.jl`
(a plain rsync copy, not a git repo, used by the p2b jobs).

## Diagnosis protocol

`sacct -j <id> --format=JobID,JobName%30,State,ExitCode,Elapsed,Timelimit,ReqMem,MaxRSS,Reason`
plus `seff <id>`, then the per-job `slurm-<name>-<id>.out` / `.err` written into
the **top level of the checkout the job was submitted from**.

Distinguish these before concluding anything:

- `CANCELLED` with `Reason=Dependency` — a chain casualty, not an independent
  failure. **This exact confusion already cost this campaign a wrong diagnosis
  on 08-22** ("never entered Slurm history" — they were dependency-cancelled).
  An `afterok` cascade killed 7 jobs at once that day.
- `TIMEOUT` — walltime, see the ceiling note below.
- `OUT_OF_MEMORY` / non-zero exit with a Julia stacktrace — real job failure.
- Never started / no record — a submission-time rejection; check the submitting
  shell's output, not the cluster.

## Known failure-class catalogue

Test against evidence; do not assume.

1. **`.CondaPkg/lock` EACCES.** Killed a prior generation. Closed 2026-08-22 by
   removing GeoIO + PythonPlot from `Project.toml` on *both* checkouts. Verify
   that is still true on both — an rsync of a local `Project.toml` could have
   put them back.
2. **`afterok` cascades.** Resubmissions must drop cross-job chaining unless the
   jobs share an output directory (the NFS single-writer rule is the only reason
   to chain). The R7 tune chain genuinely needs it (stage 2 consumes stage 1's
   `tune.csv`, stage 3 consumes stage 2's `fgstune_selected.csv`) — but consider
   `afterany` + an input-existence assert instead of `afterok`.
3. **Silent rsync.** The first `rsync` in a `for` loop once transferred nothing
   and exited 0. **Always verify cluster source by per-file md5, never by
   "rsync exited 0".**
4. **Slurm snapshots the batch script at submit time.** Editing a launcher never
   reaches a queued job; that needs cancel + resubmit.
5. **Walltime ceiling.** `m12` MaxTime is `3-00:00:00`. p1-tune scales ~N^1.6
   (R6 = 39 h), so R7 does not fit one job — hence the `STAGES` gate. The R7
   *table* jobs have the same problem (`p1-table-R7-single-a` extrapolates to
   ~100 h, and the ledger already has `R7-single-jacobi` at ~69 h) and **still
   have no splitting strategy**. Per-config splitting is the obvious move.
6. **ILU `max_pattern_entries`.** Was `2048*N`, sized at R3; R6 needs 2,059.5
   entries/row — 0.6% over — which killed `p1-table-R6-multi` (job 13206077)
   after 4 h 09 m. Raised to `8192*N` (clears R7 with ~2.2× headroom). Confirm
   the cluster copies carry the raise.
7. **Near-field cache caps.** Defaults were 4 GiB / 60 s and killed the 08-20
   campaign (R4 needed 4.47 GiB; R3 tripped the tuner's first-trial guard).
   Raised to 32 GiB / 1800 s. The build is **serial**, so the time cap is
   wall-clock. `p2b_nearfield.sh` still requests `--mem=64G` — raise it at
   submit time for any rung whose cache approaches the cap.
8. **Disk.** 117 G / 200 G on 08-22 with 8 018 jobs writing VTK. `HPC.md` says
   sweep at ~100 G *because sweeping late becomes a race*. A full disk can fail
   a job at write time. Launch `hpc-storage`; never hand-roll `rm` (a bad sweep
   killed job 13036477 on 2026-08-04). The VTK protect list is Ryan's file —
   read it, never write it.
9. **Forked launchers dropping env pins.** `p1_tune_s3.sh` was cluster-only and
   had silently lost the `FLOWPANEL_FILAMENT_REG` pin. It is now in the repo and
   deprecated in favour of `STAGES=3`. Check for other cluster-only forks.
10. **Sync scope excludes a real dependency.** THIS is what killed all ten jobs
   on 2026-08-22. The standing rsync ruling was `src/` + `benchmark/`, but
   `benchmark/phase1_case.jl:27` includes
   `examples/dji9443_trailing_edge.jl`, so **`examples/` is a hard dependency of
   every Phase-1 driver**. Removing `import GeoIO` locally without shipping
   `examples/` left both checkouts a full generation stale (45+ files) and every
   driver died at load. Sync `src/` + `benchmark/` + `examples/*.jl`.
11. **Never rsync `benchmark/results/`.** The cluster tree is the authoritative
   accepted data (561 files, per-rung R1–R6); the local one is far thinner. A
   whole-tree `benchmark/` rsync silently clobbers months of results. Exclude it
   explicitly, along with `Project.toml` (local still lists PythonPlot, which
   the cluster deliberately dropped), `Manifest.toml` and `examples/data/`.
12. **Hardware pinned in the header, not the submit line (2026-08-24).** All
   seven 021 launchers now carry
   `--partition=m12 --constraint=zen3 --exclusive --mem=0` plus a runtime assert
   on `SLURM_CPUS_ON_NODE == 128` and `HARDWARE_TAG=orc-m12-zen3-2x64c-512g`.
   Before this, `p1_table.sh`'s documented submit line omitted
   `--constraint`/`--exclusive` entirely and `ReqMem` drifted 64/128/192/256 G
   across jobs meant to be mutually comparable. **Do not override these at
   `sbatch` time** — override `--job-name` and `--time` only.

## Before resubmitting

**Sync source and verify by md5.** Ryan's standing ruling is a **whole-tree
rsync of `src/` + `benchmark/`** to both checkouts — the local tree is dirty with
mixed 018/021/022 work and `benchmark/` is untracked in full, so a partial sync
mixes source generations. Then per-file md5 diff + `import FLOWPanel` on both.

**Commit first if Ryan agrees.** The tree is still ~103 modified / 34 untracked,
and both cluster checkouts are rsync copies with no git anchor — a multi-day
campaign running on source that exists in exactly one place. Logical chunks.

**Stage scripts APPEND, they do not checkpoint-skip.** Never re-run a stage
whose CSVs already exist under `results/phase1/multi/<RUNG>/` — it duplicates
rows. Check what the failed jobs managed to write before they died, and clear or
quarantine it deliberately.

## What landed locally since the jobs were submitted — and what it does NOT change

A Phase 3 solver-lifecycle fix landed on 2026-08-23 (local, uncommitted; see
`log.md` 2026-08-23 and the amended `decision_rules.md`). Warm-start history now
advances once per completed top-level step instead of once per raw `_solve!`, and
**a one-body tuple `solve!` now breaks after the first block solve** instead of
re-solving the same system every outer iteration.

**Verified 2026-08-24: this does not change any measured quantity in the ten
failed jobs — but the earlier justification for that was wrong and is corrected
here.** Their five drivers — `phase1_table`, `phase1_tune_hpc`, `phase1_fgstune`,
`phase1_fgsprecond`, `phase2b_nearfield_cache` — call either raw
`pnl._solve!(rotor, solver)` or the public single-body `pnl.solve!(rotor, solver)`,
and **no 1-tuple call appears in any of them.**

Both of those paths *were* in fact touched by the commit, contrary to the
original claim that "the raw kernel is untouched":

- the public single-body `solve!` gained the step brackets — Dirichlet
  `src/FLOWPanel_solver.jl:414-438`, Neumann `:445-462`: `finalize_step::Bool=true`
  plus `begin_step_solution!` / `note_step_solve!` / `finalize_step_solution!`;
- the Dirichlet method swapped `influence!(body, body, backend; scalar_potential=true)`
  for the new `_source_influence!` seam;
- the FGS `_solve!` and `_krylov_launch!` kernels lost their per-kernel
  `save_solution!` / `x_prev` persistence — that is exactly where the
  once-per-step relocation lives.

The conclusion survives on a different argument. `_source_influence!`
(`src/FLOWPanel_solver.jl:384-395`) dispatches to the *original* `influence!`
call whenever `solver.S === nothing`, and `Backslash`'s
`assemble_source_potential` kwarg defaults to `false`; grep confirms no driver in
`benchmark/` opts in (the sole `assemble_source_potential!` call is
`fm051_solve_pricing.jl:106`, local-only and outside this campaign). So the
Dirichlet assembly is the same work as before. And for a cold
(`warmstart=false`, `solution_history_length=0`) single-body solve the step
commit is the same length-N copy that `_krylov_launch!` used to do, just
relocated inside the same timed region. So relaunching against the new source is
safe, and the pre-existing accepted rows stay comparable — **conditional on
`assemble_source_potential` staying `false`. The day a driver opts in, Phase-1
timings stop being comparable across that boundary.**

**Rows that ARE affected — handled 2026-08-24, not part of the relaunch.**
`benchmark/rotor_hover_solver_phase1.jl:268` (the local `phase1_costcheck`
driver) times `pnl.solve!((rotor,), (solver,))`; `..._smoke.jl:290-292`
(`backslash_iterative`) does the same. Those are now one solve rather than two,
so their `t_solve_min` and `alloc_solve` roughly halve. Done:

- the stale note strings were corrected (`phase1.jl:291`,
  `smoke.jl:278-283,292`);
- Ryan's ruling was to **quarantine** the pre-fix rows: they moved to
  `runs_pre_phase3.csv` alongside each live `runs.csv` (17 rows in
  `results/phase1/multi/`, 1 each in `results/smoke/{multi,single}/`), recorded in
  `ledger.md` under "Data quarantine". Nothing was deleted.

`..._phase1_availability.jl:128,133` also uses the 1-tuple path but records no
timings and no `notes` column, so nothing there needed correcting.

**One new operational guard will fire.** The unsteady driver now refuses to
append to an `unsteady.csv` whose first line is not the current (Phase 3) header.
`benchmark/results/phase2/multi/unsteady.csv` carries the **old** header, so any
Phase-2 unsteady rerun will error out by design and must be pointed at a fresh or
versioned path. This affects `p2_unsteady.sh` only — none of the ten failed jobs.

## What "all required jobs for Phases 1–2" means

Establish the outstanding matrix from evidence, not from this list — check which
CSVs already exist under `results/phase1/multi/<rung>/` and `results/phase2/` on
the cluster, and read `ledger.md` for what has been accepted. As of the last
known state:

- **Phase 1** — ladder frozen R1–R7; R1–R6 tuned. Outstanding: the R6 table rows
  (multi + the two single-mode configs) and the whole R7 row (tune chain, then
  the table, which still needs a splitting strategy). Cost-ceiling rule from
  `ledger.md`: `backslash_ldiv` drops at R6–R7 both modes; `krylov_gmres` drops
  from single mode at R6–R7. Drops are per config, not per rung.
- **Phase 2** — 2a table jobs (`p2_table.sh`) were **prepared 2026-08-18 and
  never submitted**; 2b near-field cache R1–R4 are the four failed jobs. The
  unsteady arm (`p2_unsteady.sh`) was also prepared and never submitted, and now
  carries the CSV-schema note above.

Submit commands are documented in each launcher's header block — read them
rather than reconstructing the `--export` syntax (note: `CONFIGS` is
colon-separated because commas collide with `--export`).

## Context you must NOT re-derive

- **The filament-family question is CLOSED.** Channel A = the filament kernel
  (`src/FLOWPanel_elements_fmm.jl:938-1067`), reached only on velocity/gradient
  passes where filaments exist → moves **numbers**, wake-carrying runs only.
  Channel B = `radius_inflation` (`:1130-1142`) → moves the near/far split, i.e.
  **cost only, never values**. Phase 1 and 2b are **family-invariant** (measured:
  R1 bit-identical under both, same iteration count), so vatistas-era and
  Gaussian-era Phase-1 rows are mutually comparable. All six launchers now
  default to gaussian.
- **The family-relabel TODO carried at the top of `log.md` is moot if these jobs
  produced no accepted results** — it existed only because the ten jobs
  snapshot-ran as vatistas. Confirm and close it out.
- FGS unsteady staleness + `simulate!` lifecycle ordering: confirmed and fixed
  2026-08-21. ILU scales ~N^1.5 (not ~N) — written up in `ledger.md`; do not
  re-measure.
- `filament_reg` is recorded by `assert_and_banner()` (`benchmark/common.jl`) and
  required by `validate_runs_csv`.
- Entry point for narrative is `log.md` (**newest first**). Do NOT read the item
  file's `## Current status` (stale, carries a banner), and do not read
  `phase_01_consistency.md` (73 KB) or `phase_02_single_step_benchmarks.md`
  (41 KB) end-to-end — grep them.

## Also open, not blocking the relaunch

- `phase3_analysis.jl`'s trajectory-identity guard flags a warm arm DIVERGED
  above 1e-8 relative checksum difference; the post-fix Phase-3 smoke arms match
  on particle count but differ 7.9e-8 / 2.4e-7. Needs Ryan's call before Phase 3
  arms run. Threshold left untouched.
- Phase 3 itself is technically unblocked but `NOT STARTED` with an unticked
  approval box; `benchmark/slurm/p3_warmstart.sh` is written and never submitted.
  Launching it needs Ryan's explicit go, separately from this relaunch.
- `_U_boundvortex` (`src/FLOWPanel_elements_fmm.jl:1559-1615`) is hardcoded
  Vatistas while `_U_boundvortex_gradient` (`:1627`) honours the global. Latent —
  it fires only on the `semiinfinite_wake=true` path, which 021 does not use.
