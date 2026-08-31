# BRAINSTORM 021 — Phase 18: context-reset handoff

**Written 2026-08-25 ~22:45 UTC while 7 jobs run. Read this file and, if you
need the campaign's rules, `decision_rules.md`. You should NOT need to read
`ledger.md`, `log.md`, or the item file to start. Do not read
`phase_01_consistency.md` (73 KB) or `phase_02_single_step_benchmarks.md`
(41 KB) end-to-end — grep them.**

Predecessor context: `phase_17_harvest_prompt.md` (the harvest brief this
session executed). Its §4 trap list and §7 quarantine list are still current
and still authoritative — re-read them before you publish any table.

---

## 1. What happened this session, in one paragraph

The 2026-08-24 launch of 13 jobs all stalled ~12 h inside Julia package
loading and four of them burned their entire walltime producing nothing. Root
cause found and fixed (§3). All wrappers were repaired, verified in
production, and the campaign relaunched. Separately, a new quality gate in the
merge script revealed that **every Phase 2 tuning descent that had completed
was hitting its wall-clock timer rather than converging** — none was a tuned
optimum. The per-rung timers were raised 5x (Ryan's call) and all five tuning
jobs were cancelled and relaunched, resuming from their checkpoints.

---

## 2. What is in flight RIGHT NOW

All `physics2` / `standby` / `--constraint=zen3 --exclusive --mem=500G
--requeue`, 64 threads, launched from `~/flowpanel-021/FLOWPanel.jl` on `orc`.

| job id | name | walltime | started | lands ~ |
| --- | --- | --- | --- | --- |
| 13477933 | p2-tune-R1 | 24 h | 08-25 22:32 | 08-26 22:30 |
| 13477934 | p2-tune-R2 | 2 d | 08-25 22:32 | 08-27 |
| 13477935 | p2-tune-R3 | 4 d | 08-25 22:32 | 08-29 |
| 13477936 | p2-tune-R4 | 6 d | 08-25 22:32 | 08-31 |
| 13477937 | p2-tune-R5 | 7 d | 08-25 22:32 | 09-01, **needs a 2nd cycle** |
| 13469157 | p1-agree-R6 | 2 d | 08-25 10:25 | 08-27 10:25 |
| 13469158 | p1-agree-R7 | 3 d | 08-25 10:25 | 08-28 10:25 |

**Nothing is queued or pending.** Every `sbatch` is Ryan's call (this session's
submissions were explicitly authorized by him, one time, for the relaunches).

---

## 3. The root cause that was fixed — do not reintroduce it

**Symptom.** All 13 jobs of the 08-24 launch printed only their bash banner
(290 bytes) and sat for ~12 h. `p2-tune-R1/R2` and `p1-agree-R1/R2` (12 h
limits) died having produced nothing. Survivors all began working within
~2 minutes of those four being killed.

**Cause.** Julia auto-precompilation racing on one NFS depot. Job 13469151
(phys-1-7) won the per-package pidfile lock and held it; the rest blocked
*before executing any user code*. Recorded verbatim in
`slurm-p1-agree-R3-13469154.err`:

```
Being precompiled by another machine (hostname: phys-1-7, pid: 2905490)
```

with a final wait of `4.32e7 ms` = 12.0 h, plus an `EACCES` writing the
`BangBangDataFramesExt` pidfile and a stale pidfile for `Meshes`.
`p2b_nearfield.sh` already carried a `JULIA_PKG_PRECOMPILE_AUTO=0` guard; the
two wrappers this campaign actually used never got it.

**Fix, now in all nine 021 wrappers** (`p1_agreement`, `p2_tune`, `p1_table`,
`p2_table`, `p1_tune`, `p1_tune_s3`, `p2_unsteady`, `p3_warmstart`,
`p2b_nearfield`):

1. `export JULIA_PKG_PRECOMPILE_AUTO=0`.
2. A **serialized precompile** using `flock -E 99 -w 3600` on
   `$HOME/.julia/flowpanel-021-precompile.lock`. `/home` is NFSv4.1 with
   `local_lock=none` (verified), so flock serializes across nodes. If `flock`
   is missing it **skips the shared build and says so loudly** rather than
   letting every job precompile at once.
3. `--compiled-modules=existing` on all 11 driver invocations. A job that
   cannot write a cache file can never take the write lock. Verified on Julia
   1.12: a never-precompiled package still loads and runs, and a **stale**
   cache is ignored and recompiled in memory — it does not silently run old
   code. Warm-cache cost is nil (3.7 s vs 5.0 s on the real stack).
4. Observability: 5-minute timestamped heartbeat, `JULIA_DEBUG=loading`,
   start/exit stamps. Measured in production: the lock serialized three jobs
   at 26–30 s each.

**Two lessons that cost this session hours — do not repeat them:**

- **`sacct TotalCPU` is not gathered on this cluster.** A COMPLETED job
  reports `00:00:00`. It is NOT evidence a job burned no CPU. Use `ps` TIME on
  the node (in **seconds**, not centiseconds).
- **Julia block-buffers stdout to a file.** An empty `.out` is not evidence a
  job is stalled. A healthy 14 h job and a wedged one looked identical.

---

## 4. The finding that matters most — Phase 2 was never converging

The new merge gate (§5) fired immediately on real data:

```
! R1 budget   0 GiB: descent TIMED OUT — n_candidates=6
! R1 budget  16 GiB: descent TIMED OUT — n_candidates=16
! R2 budget   0 GiB: descent TIMED OUT — n_candidates=6
! R3 budget   0 GiB: descent TIMED OUT — n_candidates=5
```

The descent admits up to `max_iters=20` **accepted** moves plus every rejected
candidate it explores; the timers allowed ~6. **Cause:** `TUNE_MAX_SECONDS`
was sized in Phase 16 to sit "below the typical preemption window" — a salvage
device from when the tuner had no memo and a preemption cost a whole budget.
Phase 16's candidate-level trace removed that premise, but the numbers were
never re-sized, and were then inherited by the **real-solve** driver where a
candidate is a full solve x `TUNE_REPS`, not one cheap apply.

**Ryan's call 2026-08-25: raise 5x.** Now in
`benchmark/rotor_hover_solver_phase2_tune.jl`:

| rung | R1 | R2 | R3 | R4 | R5 | R6 | R7 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| s | 18000 | 36000 | 72000 | 108000 | 216000 | 432000 | 432000 |
| h | 5 | 10 | 20 | 30 | 60 | 120 | 120 |

Walltimes in §2 were set to `4 x timer` where the partition allows. **R5 is
the exception**: 4 x 60 h = 240 h exceeds the 7-day partition cap, so R5 is a
deliberate **two-cycle** job — when it hits walltime, resubmit and it resumes.
Ryan accepted this rather than capping R5's timer.

**Do NOT change `TUNE_REPS`** (5 at R1–R4, 2 at R5+). Unlike `max_seconds`, it
**is** part of the trace provenance guard, so changing it makes every existing
trace refuse to load and discards all banked candidates. Ryan asked about
dropping to 2 and this is why the answer was no. (`t_warm` is min-of-reps —
the fastest — and 2 is the floor that still discards first-touch cost.)

**70 candidates are banked** across 10 trace files and were confirmed to
replay cleanly on relaunch (no provenance/schema errors; trace line count
unchanged at 70). Superseded rows are at
`results/phase2/multi/R{1,2,3}/tune_phase2_pre20260825_oldtimer.csv.bak`.

---

## 5. Code changed this session — all synced and md5-verified

Eleven files, all local edits synced to `~/flowpanel-021/FLOWPanel.jl` and
verified md5-identical on both sides.

- **9 wrappers** in `benchmark/slurm/` — §3 fixes.
- **`benchmark/rotor_hover_solver_phase2_tune.jl`** — 5x `_MAXSEC`.
- **`benchmark/p021_merge.jl`** — two changes:
  1. **`residual_backend` added to the agreement key.** The O(N^2)
     cross-check writes the same configs at the same rung into the same file
     (`RESIDUAL_BACKEND` is a CSV column, not a path component). Keying on
     `(rung, config)` alone reported R1's 9 legitimate rows as duplicates —
     the NFS-clobber signature — on correct data. Now keyed on
     `(rung, config, residual_backend)`; a real same-backend duplicate is
     still caught. Validated on the live R1 file: 9 rows, zero warnings.
  2. **Quality gate.** `tune_timed_out=true`, `cache_capped=true`, and
     `bc_certified=false` are now reported problems, so `STRICT=1` refuses a
     clean merge. **Below-floor rows are exempt** — they carry
     `cache_capped=true` and `bc_certified=false` by construction and would
     otherwise fire on every legitimate below-floor answer at R6/R7.

Gated rows are still written to the merged CSV — the gate blocks the *claim*
of a clean merge, it does not delete data.

**Everything is uncommitted.** The working tree was already ~100 modified /
34 untracked before this session (`benchmark/` is untracked in full) and both
cluster checkouts are rsync copies with no git anchor. This was flagged on
2026-08-22 and is unchanged. Committing in logical chunks is worth raising.

---

## 6. Results so far

**Phase 1 is complete R1–R5.** Seven solvers agree, reference-free:

| rung | n_panels | configs | max pairwise relL2 |
| --- | --- | --- | --- |
| R1 | 8,016 | 7 (+2 direct) | 1.73e-5 |
| R3 | 28,752 | 7 | 1.90e-5 |
| R4 | 58,192 | 7 | — |
| R5 | 108,240 | 5 | — |

R3 detail (all `bc_certified=true`): `backslash_ref` 45.75 s, `krylov_ilu`
35.70 s / 11 it, `fgs` 34.47 s, `fgmres_fgs` 39.48 s / 10 it, `krylov_jacobi`
44.51 s / 16 it, `krylov_gmres` 118.67 s / 69 it, `backslash_coupled`
140.59 s. `CT_laplace` spans 5.9e-7 absolute across all seven.
`backslash_coupled` is a consistent `bc_rel_l2` outlier (4.6e-5 at R3 and R4,
four decades worse than `backslash_ref`) while landing on the same solution.
`CT_laplace` climbs monotonically R1→R5 (0.06041 / 0.06102 / 0.06124 /
0.06143) — **not mesh-converged**. `krylov_gmres` scales ~N^1.4 R3→R5 with
iteration count nearly flat (69→76).

**Phase 2 has no publishable rows.** All prior rows were timed-out and are
set aside; the relaunched descents are the first that can converge.

---

## 7. What to do next

1. **Wait.** Nothing needs submitting until a job lands.
2. **When `p1-agree-R6`/`R7` land** (08-27 / 08-28): Phase 1 is complete.
   Expect `krylov_ilu` at R7 to possibly not fit (~23 GiB of pattern entries)
   — that is a RESULT, not a failure.
3. **When a `p2-tune` job lands**: check `tune_timed_out` FIRST. The merge
   gate will tell you. If it timed out again even at 5x, that is a finding,
   not something to paper over.
4. **`p2-tune-R5` will need a second submission** after it hits its 7-day
   walltime. Resubmit the same line; it resumes from the trace.
5. **Then**: the campaign's actual open question — does the real-solve winner
   AGREE with the apply proxy at R1–R5? That decides whether R6/R7 can be
   tuned by proxy at all. Deferred, and Ryan's call.
6. **Step C (`p2_table.sh`) is still ungated and unrun.** It needs that rung's
   Step A landed and the pre-W3 `phase2.csv` moved aside first.

**Harvest command** (from `~/flowpanel-021/FLOWPanel.jl`):

```bash
julia --project benchmark/p021_merge.jl
RUNGS=R1:R2:R3:R4:R5 BUDGETS=0.0:16.0:128.0:500.0 STRICT=1 \
  julia --project benchmark/p021_merge.jl
```

**Monitoring gotchas:** `.err` files now need **grepping, not reading** —
`JULIA_DEBUG=loading` emits a bounded burst of "Rejecting cache file" lines at
startup (normal: Julia rejecting the shipped stdlib caches built with
non-default flags). Real errors are at the tail. `squeue`/`sacct` need
`ssh orc 'bash -lc "..."'` — a plain non-interactive ssh has no Slurm on PATH.

---

## 8. Is Phase 3 ready for work?

**Not ready to launch. Ready to plan. And there is one piece of real work that
is unblocked right now.**

- **Approval gate.** Phase 3 is `NOT STARTED` with its approval box unticked.
  "Completing a phase does not authorize the next" — starting it needs Ryan's
  explicit go. `phase_03_launch_prompt.md` (139 lines) is the planning brief
  and `phase_03_warmstart.md` (~3.5 KB) is the spec; read both in full.
- **Its numeric baselines are downstream of what is running.** Phase 3's cold
  single-step baselines come from `phase2.csv` — the Step C output, which has
  never been run and is itself gated on Step A tuning landing. That is now
  1–7 days out. Frozen knobs from `results/phase1/multi/*` are fine and
  family-invariant.
- **The "FGS iteration capture" blocker is CLOSED — do not re-do it.**
  `phase_03_launch_prompt.md:103` calls it "pending callback plumbing (Phase
  3)"; that document was written 2026-08-22 and **superseded the next day**.
  The 08-23 work implemented it end-to-end: `FGSSolver.niter`
  (`src/FLOWPanel_solver.jl:1383`), a counting callback with an explicit
  off-by-one correction (`:1674-1709`), and the `SolveStepStats` /
  `note_step_solve!` / `publish_step_stats!` lifecycle under the section header
  `PER-STEP SOLVE STATISTICS ... (BRAINSTORM 021 Phase 3)` at `:53`. It is
  wired into `benchmark/rotor_hover_solver_unsteady.jl:241`, which reads
  `s.niter` for BOTH Krylov and FGS. See `phase_19_fgs_niter_followups.md`.
- **Ready:** the driver exists (`benchmark/rotor_hover_solver_unsteady.jl`,
  reuses `examples/rotor_hover_pressure_comparison.jl` with
  `RHPC_SETUP_ONLY=1`); its launcher `benchmark/slurm/p2_unsteady.sh` exists
  but has **never been submitted**, so no unsteady HPC rows exist at all.
  Phase 0's warmstart plumbing landed and is unit-tested. Phase 3 is
  Gaussian-from-birth, so it carries no filament rerun debt.
- **Open questions for Ryan when planning:** whether Phase 3 runs at R1 only
  or across rungs; how many steps reach the wake-developed regime; and whether
  per-step cost under the two filament families matters to the baseline (the
  existing `unsteady_family_ab_20260822.csv` timing columns are
  noise-confounded and **must not be cited**).

**Recommendation:** while the campaign runs, plan Phase 3 and build the FGS
iteration-capture plumbing. Do not submit any Phase 3 job.

---

## 9. Still Ryan's call — do not decide

- Every `sbatch`, including Step C and R5's second cycle.
- Starting Phase 3.
- The R6/R7 tuning objective (real-solve vs apply proxy).
- Any notebook entry — **offer, don't write**. One is owed for this session:
  the stall, the root cause, the fix, the timer finding, and the R1–R5
  agreement tables. Ask for the day header and how much detail per section.
- Committing the working tree.
