# Reset prompt: 021 champion adoption — resume reruns live + harvest + FGS-scaling question (2026-09-21, supersedes 20260920)

Copy everything below the line into a fresh agent session started in
`~/Dropbox/research/projects/FLOWPanel.jl`.

---

## Where this stands

BRAINSTORM 021 is CLOSED AT PROMOTION (champion = dagteam + f32full +
`numactl --interleave=0-3 --cpunodebind=0-3` @ j16/BLAS 1, R4 cold solve
4.475 s = 2.26×; story `fgs_acceleration_summary_20260919.md`). Read first:
`CLAUDE.md`, `agent_policies/WORKFLOW.md` + `TESTING.md` + `HPC.md` before
corresponding work. Prior context: `champion_adoption_reset_prompt_20260920.md`
and the two provenance files (each now has a "Resume rerun 2026-09-21"
addendum).

The 2026-09-21 session harvested the first pass of both campaigns and
launched warm-start resume reruns.

### Harvested (first pass, 13777133 + 13778533)

R4 thread scaling, median accepted solve time (s):

| j | FGS-dagteam | FGS-colored | krylov_ilu (CAVEAT: tuner budget-500 warm times, not certified measurement) |
|---|---|---|---|
| 1 | 33.07 | 35.58 | 249.9 |
| 8 | 5.97 | 9.07 | 15.66 |
| 16 | 4.50 | 7.26 | 12.03 |
| 32 | 4.42 | 6.91 | 9.16 |
| 64 | 6.42 | 12.30 | 7.70 |

- Both FGS families plateau at j*=16 by the ruling (<~10%/doubling above j*,
  in `thread_scaling_provenance_20260919.md`); both REGRESS at 64 (0.69× /
  0.56×). krylov_ilu still scales at 64 (+19%/doubling) but never beats
  dagteam anywhere. **No R4 krylov_ilu_nfcache data yet** (see live jobs).
- R1–R2 (13778533): 8/14 arms completed (R1 j4–64, R2 j16–64), full tables
  in the harvest (R1 winner krylov_ilu_nfcache 0.165–0.31 s; FGS plateaus
  ≤j16 there too). Six arms timed out at the 24 h wall mid-p2tune (R1 j1/j2,
  R2 j1/j2/j4/j8). **R2 j64 is anomalous** (fgs 1.20→10.7 s vs j32; treat as
  contention pathology, discard from scaling reads pending repeat).
- Duplicate-row scan of all 19 cluster `tune_phase2.csv`: CLEAN (the "16" vs
  "16.0" resume-bug artifact did not land in this dataset). Minor: j16
  thread-scaling `tune_phase2.csv` has a duplicated header row — make harvest
  scripts skip repeated headers.
- Expected finding dagteam@j1 STATUS_fgs_calibrate=FAILED did NOT occur —
  calibrate passed at j1 (datum: dagteam calibrates single-threaded).

### Root cause found + fixed: ilu_measure FAILED (all 13777133 arms)

`run_r4_thread_scaling.slurm.sh` hard-coded `TUNE_SEED_B0=10:0.6:6`; P=10 is
error-tolerance-violating (>1e-6) at R4 under FM `f4d6b671` (budget-500
winner needed P=12). The tuner warns-and-continues on a failed budget →
STATUS_ilu_tune=ok but no budget-0 row → `rotor_hover_solver_phase2.jl:103`
hard-errors → STATUS_ilu_measure=FAILED at every j. Fix (commit `946cec2` on
`fastmultipole`, tag `campaign/p021-resume-source-20260921`): default B0
seed = the budget>0 seed `15:0.55:32` (env-overridable) + warm-start resume
in BOTH launchers (`RESUME_FROM_JOB_ID=<old id>` reuses the old run dirs,
skips `STATUS_*=ok` stages, tuner stages always re-run and rely on row-level
(rung,budget,julia_threads) resume; p2/measure keep the ok-skip because they
append to phase2.csv; old logs/pins go to `logs.before.<new id>/`;
COMPLETED removed on resume). r12 launcher wall 24→48 h. Provenance commit
`9b1f411`; both pushed to orc branch `p021-thread-scaling-20260919`.

## LIVE JOBS (submitted 2026-09-21, all 11 tasks started immediately on m12; judge by outputs, never sacct)

Data root unchanged: `/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/`
(SAME run dirs as the first pass — resume reuses them).

- **13829231** — thread-scaling resume, array 0–4 → j=1/8/16/32/64, run dirs
  `thread-scaling-j<J>-13777133/` (note: OLD job id in dir name). Per arm:
  FGS arm skips; budget-0 knob descent (≤10 h cap) + krylov_ilu AND
  krylov_ilu_nfcache measurement at R4 (`CONFIGS=krylov_ilu,krylov_ilu_nfcache
  MEM_BUDGETS=500`; nfcache reads budget-500 knobs, uncached reads the new
  budget-0 knobs). Julia threads = j, **BLAS = 1** (champion convention,
  HARDWARE_TAG `...-blas1`). Worktree
  `/home/rander39/campaigns/p021-thread-scaling-20260919/` now at exec
  `0c002ae` (tag `campaign/p021-thread-scaling-exec-20260921`, partial pick:
  only its own launcher; pins.toml updated).
- **13829232** — R1–R2 resume, array 0,1,7,8,9,10 → R1 j1/j2 + R2
  j1/j2/j4/j8, run dirs `r12-champion-<rung>-j<j>-13778533/`. fgstune +
  fgsprecond skip (certified staircases reused); p2tune resumes row-level;
  then p2. **BLAS = j** (this campaign's convention). Worktree
  `/home/rander39/campaigns/p021-r12-champion-20260919/` at exec `e2b1410`
  (tag `campaign/p021-r12-champion-exec-20260921`; pins.toml updated).

## NEXT ACTIONS

1. **Monitor via `hpc-monitor`; harvest via `harvester`** when arrays land.
   Priority deliverable (Ryan asked explicitly): the R4 table above extended
   with certified krylov_ilu and **krylov_ilu_nfcache** columns vs j from
   13829231 (`<rundir>/ilu/phase2.csv`; knobs in `<rundir>/ilu/tune_phase2.csv`
   — budget-0 row must now exist). Verify budget-0 rows landed at all five j.
   Then merge the six resumed R1–R2 arms into the R1–R2 tables (per-run
   `phase2/phase2.csv`; check for duplicate rows anyway — p2 ran fresh so
   there should be none).
2. **Plateau verdict discussion with Ryan** (deliberately DEFERRED until 1
   lands): FGS j*=16 is solid; open question Ryan posed back: should
   krylov_ilu('s family) carry weight in pruning at all given it never wins
   on absolute time? If not, everything above 32-thread caps is prunable.
   Ryan leans practical; present the completed table and let him rule.
3. **Ryan-gated, discussion first: why doesn't FGS scale past ~16 threads?**
   Ryan flagged this as the likely next investigation AFTER 1+2 and his
   approval. Existing evidence (suggestive, no direct measurement yet):
   - krylov_ilu keeps scaling to j=64 in the SAME jobs/placement/machine, so
     it is FGS-specific, not a node ceiling.
   - Both FGS variants REGRESS beyond 32 (0.69×/0.56×) — active contention,
     not just saturation; R2 j64 blow-up (~9×) on the smaller rung points the
     same way (worse for small rungs = less work per sync).
   - BLAS is measured-inert for FGS (per-leaf gemvs below BLAS threading
     thresholds, 2026-09-19 ruling) — so it is not BLAS oversubscription.
   - Placement is load-bearing (socket-membind collapsed dagteam to 1.09×;
     champion = interleave 0-3): consistent with memory-bandwidth/NUMA
     sensitivity; j=64 uses ALL of socket 0's cores (no headroom for GC/OS).
   - Structural suspects: Gauss-Seidel sweep dependency chains (dagteam DAG
     width / coloring class sizes bound usable parallelism), fine-grained
     tasks at leaf≈100, per-sweep synchronization cost growing with j.
   - NOT yet measured: sweep-phase timing breakdown vs j, DAG width/level
     occupancy, bandwidth counters (likwid/perf), GC/spin time. A cheap first
     probe: rerun the j-ladder A/B fixture with per-phase timers, or a
     j16-vs-j64 `perf stat` pair on the existing harness. Discuss scope with
     Ryan before building anything.
4. **Ryan-gated follow-ons** (unchanged): R3+ re-runs (parked on plateau
   verdict); optional R1–R2 f32full arms (resubmit 13829232's launcher env
   with `FGS_DAGTEAM_PRECISION=f32full` AFTER f64 staircases certify);
   optional zen3 BLAS A/B rider at the R4 champion (expected null).
5. **Standing Ryan-gated ledger**: origin pushes (both repos + v21/v22/v23 +
   campaign tags incl. the three new `*-20260921` tags), notebook entries
   owed (v21, diagnostics ladder, v22 chunked, NUMA, v23 promotion, + these
   campaigns once harvested — offer, don't write), WeakKeyDict/warmstart fix
   (`_publish_block_gs_status!`, from `7fbd68a`), `:dagteam` unit test
   missing in `test/runtests_unit_solver.jl`, 3 RECENT p018 runs awaiting
   archive approval, hpc-storage archive-pass report collection.

## Traps (2026-09-21 update; older ones in the 20260920 prompt still bind)

- Resume run dirs carry the OLD job ids (13777133/13778533) — do not look
  for dirs named 13829231/13829232.
- `logs.before.<new job id>/` inside each resumed run dir holds the FIRST
  pass's logs/pins — harvest scripts must not double-count them.
- STATUS_ilu_tune=ok does NOT prove all budgets landed (warn-and-continue) —
  always check tune_phase2.csv rows directly. Same for any tuner stage.
- BLAS conventions DIFFER by campaign: R4 thread-scaling = BLAS 1; R1–R2 =
  BLAS j. The FGS family is measured BLAS-insensitive but the ILU
  factorization was never BLAS-swept — flag before any cross-campaign ILU
  comparison.
- Tolerances are per order AND per environment — never carry a tolerance;
  every point re-staircases. f32full certified at R4/zen3 only.
- The krylov_ilu column in the table above is superseded the moment 13829231
  lands — replace it, don't average it.
- Thread-scaling worktree branch predates the r12 launcher — its exec is a
  partial pick; don't "fix" the missing file there.
- Local runs ≤4 threads; `ssh orc` needs a live ControlMaster socket; Slurm
  CLI needs a login shell (`bash -lc` + `source /etc/profile`); judge runs by
  outputs, never sacct.

## House rules (binding, unchanged)

HPC submission Ryan-gated; monitoring via `hpc-monitor`; storage via
`hpc-storage` (400 G cap); notebook writes Ryan-gated (offer, don't write);
dated status/provenance files in BRAINSTORM/021; each campaign in its own
worktree, never a shared live checkout while jobs run.
