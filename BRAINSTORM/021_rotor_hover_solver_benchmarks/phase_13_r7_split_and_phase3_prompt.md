# Handoff — BRAINSTORM 021: R7 table splitting strategy, then (maybe) Phase 3

Written 2026-08-24 at a context reset, immediately after relaunching 13 jobs.
**The 08-22 failure is diagnosed, repaired and closed — do not re-open it.**

## Your task, in order

1. **Design the R7 table splitting strategy.** This is the real work. It is a
   design question with no accepted answer yet; see the section below for the
   measured cost data you need and the constraint that binds.
2. **Check on the 13 relaunched jobs** and report a compact table. They were
   `PENDING` with a projected start of 2026-08-24T17:49 (cluster busy) when the
   session ended.
3. **Phase 3, only if Ryan says so.** It is technically unblocked but has an
   unticked approval box, and it carries one unresolved numeric question (below).
   `benchmark/slurm/p3_warmstart.sh` is written, hardware-pinned and never
   submitted. Launching it needs Ryan's explicit, separate go.

## Access

`ssh orc` needs a **live ControlMaster socket**; without one it prompts for 2FA
and hangs. If the socket is dead, ask Ryan to run `! ssh orc` in the session so
the output lands in the conversation and the socket comes up. Do not try to
authenticate yourself. The socket was live for the whole 08-24 session.

Read `agent_policies/HPC.md` before touching anything. On the login node,
`/apps/instructions_for_ai_agents/BYU_ORC_AGENTS.md` is the site policy. Slurm
binaries live in `/apps/slurm/latest/bin` (not on `PATH`). Julia needs a login
shell: `ssh orc 'bash -lc "module load julia && ..."'`.

**Delegation:** route job status, log tailing and disk reporting to
`hpc-monitor`, and any VTK sweeping to `hpc-storage`. Keep diagnosis reasoning,
repairs and submission inline. **Job submission is Ryan's call — never a
subagent's, and never without his go.**

---

## 1. The R7 table splitting problem (your main task)

### The constraint

`m12` MaxTime is `3-00:00:00`. There is no partition with a longer ceiling.
`p1-table-R7-single-a` extrapolates to **~100 h** and the ledger already carries
`R7-single-jacobi` at **~69 h** — the latter fits only barely, with no margin for
a slow node or a restart.

### The measured cost data you need

From `slurm-p1-table-R6-multi-13206077.out` (the 08-20 job that died on the ILU
limit after 4 h 09 m), R6 = 212,108 panels, **multi** mode, 64 threads:

| Config | setup [s] | solve [s] | niter |
| --- | --- | --- | --- |
| krylov_gmres | 1.38 | 1364.29 | 82 |
| krylov_jacobi | 8830.7 | 443.73 | 22 |

Note the shape: `krylov_jacobi`'s cost is almost entirely **setup**, and
`krylov_gmres`'s is almost entirely **solve**. A split that assumes uniform
per-config cost will be wrong. `ledger.md` also records that the Barba ILU
preconditioner scales **~N^1.5, not ~N** (measured, do not re-measure) — at R6 it
needs 2,059.5 pattern entries/row, and `max_pattern_entries` is now `8192*N`
(clears R7 with ~2.2x headroom).

The frozen ladder doubles N per rung: R6 = 212,108, R7 = 419,276.

Cost-ceiling drops from `ledger.md` (per config, not per rung):
`backslash_ldiv` drops at R6–R7 in **both** modes; `krylov_gmres` drops from
**single** mode at R6–R7. So the R7 config sets are:

- **multi**: `krylov_gmres`, `krylov_jacobi`, `krylov_ilu`, `fgs`, `fgmres_fgs`
- **single**: `krylov_jacobi`, `krylov_ilu`, `fgs`, `fgmres_fgs`

### What is already known about how to split

`p1_table.sh` takes `CONFIGS` as a **colon-separated** list at `--export`
(commas collide with `--export` syntax), so **per-config splitting works today
with no code change** — it is how the 08-24 relaunch scoped R6-multi down to
`krylov_ilu:fgs:fgmres_fgs` and put `krylov_jacobi` in its own job. That is the
obvious move and the one Ryan has implicitly already accepted at R6.

The open questions are the ones per-config splitting does **not** answer:

- Is one config per job enough at R7, or does a single config still exceed 72 h?
  `krylov_gmres` multi at R6 was 1364 s/solve with `k_reps=5`; work the R7
  projection from the ladder's N-doubling and the measured scaling, and say
  plainly which configs still do not fit.
- If a single config does not fit, the lever is `K_REPS` (`p1_table.sh` defaults
  to 5) or splitting the setup from the solve. **Check `decision_rules.md` for
  the campaign's timing protocol before touching `k_reps`** — min-of-k is part of
  the accepted measurement contract, not a free knob.
- `krylov_jacobi`'s 8830 s **setup** at R6 is the thing that will bite: it is
  serial-ish preconditioner construction, so it does not shrink with threads.
- **Stage scripts APPEND, they do not checkpoint-skip.** Any split must be
  designed so a partial job leaves a re-runnable remainder, and so a re-run never
  duplicates rows. Before submitting anything, check what exists under
  `benchmark/results/phase1/<mode>/R7/` and quarantine deliberately.

### Do not forget

R7's table is **downstream of the R7 tune chain** (`13396783/784/785`), which
produces `tune.csv` → `fgstune_selected.csv` → `fgsprecond.csv` under
`benchmark/results/phase1/multi/R7/`. The table cannot run until those land.
Design the split now; submit it after the chain completes.

---

## 2. The 13 jobs relaunched 2026-08-24

All hardware-pinned (`m12`/`zen3`/`--exclusive`/`--mem=0`) and validated with
`sbatch --test-only` before submission.

| Job | id | Checkout | Notes |
| --- | --- | --- | --- |
| p1-table-R6-multi | 13396780 | main | `CONFIGS=krylov_ilu:fgs:fgmres_fgs` only |
| p1-table-R6-single-a | 13396781 | main | same 3 configs, 1 thread |
| p1-table-R6-single-jacobi | 13396782 | main | `krylov_jacobi` alone |
| p1-tune-R7-s1 | 13396783 | main | `STAGES=1` |
| p1-tune-R7-s2 | 13396784 | main | `STAGES=2`, **afterany** on 783 |
| p1-tune-R7-s3 | 13396785 | main | `STAGES=3`, **afterany** on 784 |
| p2-table-R2-multi | 13396805 | main | 6 configs |
| p2-table-R3-multi | 13396807 | main | 6 configs |
| p2-table-R4-multi | 13396809 | main | 6 configs |
| p2-table-R5-multi | 13396811 | main | 6 configs |
| p2b-nearfield-R2 | 13396813 | p2b | raised caps |
| p2b-nearfield-R3 | 13396815 | p2b | raised caps |
| p2b-nearfield-R4 | 13396817 | p2b | raised caps |

**R6-multi runs only 3 configs on purpose:** `krylov_gmres` and `krylov_jacobi`
were already accepted on 2026-08-20 (job 13206077, before it died on the ILU
limit). The stage scripts append, so re-running them would duplicate rows.

**Phase 2a runs 6 configs, not 9:** the `*_nfcache` variants consume 2b's
`tune_cached.csv`, which R2–R5 do not have yet. Once 2b R2–R4 land, the nfcache
variants become runnable — that is a follow-on round, not a re-run.

**Watch for:** `p1-table-R6-single-a` (13396781) is a genuine timeout risk.
Multi-mode `krylov_gmres` alone took 1364 s/solve at R6; single-threaded
`krylov_ilu`/`fgs`/`fgmres_fgs` on 212k panels against a hard 72 h ceiling is not
obviously safe. If it TIMEOUTs, that is data for the R7 split design, not a bug.

Eight BRAINSTORM-018 jobs (`13396725`–`13396732`) were RUNNING and untouched.
**Do not disturb them.**

---

## 3. Phase 3 — needs Ryan's go, and has one open number

- `phase3_analysis.jl`'s trajectory-identity guard flags a warm arm DIVERGED
  above a **1e-8** relative checksum difference. The post-fix Phase-3 smoke arms
  match on particle count but differ **7.9e-8 / 2.4e-7**. **Needs Ryan's call
  before Phase 3 arms run.** The threshold was deliberately left untouched.
- `benchmark/slurm/p3_warmstart.sh` is written, now hardware-pinned, and has
  never been submitted.
- Phase 3 itself is `NOT STARTED` with an unticked approval box in the item file.

---

## Context you must NOT re-derive

- **The 08-22 failure is closed.** Single cause: the cluster's
  `examples/dji9443_trailing_edge.jl` still had `import GeoIO` at line 2 after
  GeoIO was dropped from `Project.toml`, because the standing sync ruling was
  `src/` + `benchmark/` only. `benchmark/phase1_case.jl:27` includes that file,
  so all five drivers died at load in ≤3 min. Two jobs were `afterok` casualties,
  not failures. Repaired, md5-verified (0 mismatches over 236 files on both
  checkouts), and both checkouts now load clean. Full write-up: `log.md`
  2026-08-24 (newest first); catalogue items 10–12 in
  `phase_12_relaunch_prompt.md`.
- **The sync recipe is now fixed**, and getting it wrong is destructive:
  `rsync -a` of `src/ benchmark/ examples/`, **no `--delete`**, excluding
  `Project.toml`, `Manifest.toml`, `examples/data/` and **`benchmark/results/`**.
  That last exclusion is the dangerous one — the cluster results tree is the
  authoritative accepted data (561 files, per-rung R1–R6) and the local tree is
  far thinner. Verify by **per-file md5, never by "rsync exited 0"**.
- **Hardware is pinned in the SBATCH headers of all seven 021 launchers** and
  must not be overridden at `sbatch` time. Override `--job-name` and `--time`
  only. Audit finding: every job ever run for this campaign was already
  zen3/exclusive/128-CPU, so the ladder *is* mutually comparable; only `ReqMem`
  had drifted. Node: 2 sockets x 64 cores, ThreadsPerCore=1, 512 GiB.
- **Chain with `afterany`, never `afterok`** (Ryan, 2026-08-24). `p1_tune.sh`
  carries a stage-input assert that names the missing CSV and exits 1.
- **The filament-family question is CLOSED.** Channel A = the filament kernel
  (`src/FLOWPanel_elements_fmm.jl:938-1067`), velocity/gradient passes where
  filaments exist → moves **numbers**, wake-carrying runs only. Channel B =
  `radius_inflation` (`:1130-1142`) → **cost only, never values**. Phase 1 and 2b
  are family-invariant (measured: R1 bit-identical, same iteration count). All
  launchers default to gaussian.
- **The family-relabel TODO is CLOSED** (2026-08-24): the ten jobs produced zero
  rows, so there was nothing to relabel.
- **The one-body `solve!` early break** (`src/FLOWPanel_solver.jl:2000-2011`)
  landed 08-23: a single-body tuple now does one block solve, not two. Note
  strings fixed at `benchmark/rotor_hover_solver_phase1.jl:291` and
  `..._smoke.jl:278-283,292`; pre-fix rows quarantined to `runs_pre_phase3.csv`
  (17 rows phase1, 1 each smoke multi/single) and recorded in `ledger.md` under
  "Data quarantine". This affects **no** HPC driver — none uses the 1-tuple path.
- FGS unsteady staleness + `simulate!` lifecycle ordering: confirmed and fixed
  2026-08-21. ILU scales ~N^1.5 (not ~N) — in `ledger.md`; do not re-measure.
- `filament_reg` is recorded by `assert_and_banner()` (`benchmark/common.jl`) and
  required by `validate_runs_csv`.
- Entry point for narrative is `log.md` (**newest first**). Do NOT read the item
  file's `## Current status` (stale, carries a banner), and do not read
  `phase_01_consistency.md` (73 KB) or `phase_02_single_step_benchmarks.md`
  (41 KB) end-to-end — grep them.

## Also open, not blocking

- **The tree is still uncommitted.** ~103 modified / 34 untracked, and both
  cluster checkouts are rsync copies with no git anchor — a multi-day campaign
  running on source that exists in exactly one place. Ryan declined to commit on
  2026-08-24 ("relaunch now"). Worth re-raising once the jobs land.
- Disk: `~/projects` was **104 G** against a 200 G cap on 2026-08-24 (filesystem
  itself 139 G / 2.0 T). Above `HPC.md`'s ~100 G sweep trigger. Launch
  `hpc-storage`; never hand-roll `rm` (a bad sweep killed job 13036477 on
  2026-08-04). The VTK protect list is Ryan's file — read it, never write it.
- Phase 2 unsteady (`p2_unsteady.sh`) is blocked by design: the driver refuses to
  append to an `unsteady.csv` whose first line is not the Phase-3 header, and
  `benchmark/results/phase2/multi/unsteady.csv` carries the old one. Point it at
  a fresh or versioned path.
- `_U_boundvortex` (`src/FLOWPanel_elements_fmm.jl:1559-1615`) is hardcoded
  Vatistas while `_U_boundvortex_gradient` (`:1627`) honours the global. Latent —
  fires only on `semiinfinite_wake=true`, which 021 does not use.
