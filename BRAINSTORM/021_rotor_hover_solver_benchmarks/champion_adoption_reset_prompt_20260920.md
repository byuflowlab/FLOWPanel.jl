# Reset prompt: 021 champion adoption — live campaigns + harvest (2026-09-20, supersedes 20260919g)

Copy everything below the line into a fresh agent session started in
`~/Dropbox/research/projects/FLOWPanel.jl`.

---

## Where this stands

BRAINSTORM 021 FGS acceleration is CLOSED AT PROMOTION (champion = dagteam +
f32full + `numactl --interleave=0-3 --cpunodebind=0-3` @ j16/BLAS 1, R4 cold
solve 4.475 s = 2.26×; story `fgs_acceleration_summary_20260919.md`, config
`benchmark/retained_r4_champion.toml`). Read first: `CLAUDE.md`,
`agent_policies/WORKFLOW.md` + `TESTING.md` + `HPC.md` (before corresponding
work).

The 2026-09-19 session executed Ryan's adoption directive. Key finding first:
**no production chain outside `benchmark/` uses FGS** (repo-wide inventory;
018/022/026/032 use other solvers) — so "re-run the FGS runs" = re-running
the 021 benchmark suite itself, and no physics-null was needed.

**Done and committed** (branch `fastmultipole`: `db4f053`, `d3ae35e`,
`4a8a852`, + provenance commits; pushed to the orc unified repo as branch
`p021-thread-scaling-20260919`; origin pushes remain Ryan-pending):

1. **General policy (Ryan-ruled 2026-09-19)** in
   `benchmark/rotor_hover_solver_phase2_tune.jl`: memory ladder
   0/16/32/64/128/500 GiB, `machine_max_threads(b) = min(64, b ÷ 2)`
   (~2 GiB/thread; budget 0 uncapped); tuning identity = (rung, budget,
   julia_threads) for rows, resume, and trace files; budgets whose cap < the
   running thread count are skipped. `PHASE2_OUTDIR` override in tuner +
   `rotor_hover_solver_phase2.jl` (concurrent-job isolation; the 2026-08-18
   NFS append hazard). **Fixed a pre-existing silent bug**: row-level resume
   compared "16" vs "16.0" and never matched → standby requeues re-tuned and
   DUPLICATED landed budgets — check historical cluster `tune_phase2.csv`
   files for duplicate (rung, budget) rows at harvest time.
2. **FGS order plumbing**: `FGS_SWEEP_ORDER` / `FGS_DAGTEAM_PRECISION` env in
   `benchmark/phase1_case.jl` (`fgs_order_kw()`), splatted into every
   `FGSSolver`/`FGSPreconditioner` ctor in fgstune, fgsprecond, phase2;
   defaults unchanged (lexicographic/f64); order recorded in phase2 config
   strings. `*margin_verify.jl` NOT plumbed (not in the pipeline).
3. **BLAS ruling (measured)**: BLAS thread count is inert for the FGS family
   (local R1@j4 dagteam: BLAS 1/2/4 within 0.3% — per-leaf gemvs sit below
   BLAS threading thresholds). The champion was never BLAS-swept (BLAS 1 by
   harness convention); expected null at R4 too (leaf=100). R1–R2 re-run uses
   uniform BLAS=j per task (backslash honesty; no cross-BLAS tolerance carry
   since the whole task shares one environment).

## LIVE JOBS (submitted 2026-09-19, judge by outputs, never sacct)

Data root: `/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/`.
Per-arm verdicts = `STATUS_*` files; task-level = `COMPLETED` marker.

- **13777133** — R4 thread-scaling array `_0–_4` → j=1/8/16/32/64
  (`thread-scaling-j<J>-13777133/`). Per j: FGS dagteam-f32full + colored
  calibrate (evaluator ≤1e-6) + A/B trials; krylov_ilu apply-knob re-descent
  at budgets 0+500 then krylov_ilu(+nfcache) measurement (retune REQUIRED:
  FastMultipole apply path changed 4c0f1b8f→f4d6b671 since its 08-25 knobs).
  Champion placement, BLAS 1, `HARDWARE_TAG=orc-m12-zen3-socket0-ilv0-3-blas1`.
  Worktree `/home/rander39/campaigns/p021-thread-scaling-20260919/` (exec
  `28eef244`). Provenance `thread_scaling_provenance_20260919.md`.
  dagteam@j1 is uncovered by any gate — STATUS_fgs_calibrate=FAILED there is
  a FINDING, not a job failure.
- **13778533** — R1–R2 champion re-run array `_0–_13` → (R1,R2) ×
  j=1/2/4/8/16/32/64 (`r12-champion-<rung>-j<j>-13778533/`). Per task:
  fgstune (dagteam **f64** knob descent + tolerance staircase = the recorded
  new accuracy) → fgsprecond (SWEEP_LADDER_1E6=1) → phase2 tuner (full
  ladder, caps prune) → phase2 full CONFIGS. BLAS=j,
  `HARDWARE_TAG=orc-m12-zen3-socket0-ilv0-3`. Own worktree
  `/home/rander39/campaigns/p021-r12-champion-20260919/` (exec `f85ac255`;
  separate because 13777133 runs from the other worktree — concurrency rule).
  Provenance `r12_champion_rerun_provenance_20260919.md`.

## NEXT ACTIONS

1. **Monitor via `hpc-monitor`** (read-only). When arrays finish, **harvest
   via `harvester`**: scaling curves (median accepted solve time vs j) for
   FGS-dagteam / FGS-colored / krylov_ilu / krylov_ilu_nfcache from 13777133
   (`fgs-trials/results/ab_summary.toml`, `ilu/tune_phase2.csv` +
   `phase2.csv` per run dir); R1–R2 tables from 13778533 run dirs. Also scan
   harvested tune CSVs for the historical duplicate-row artifact (item 1
   above).
2. **Plateau verdict to Ryan** (his decision input for pruning thread tiers
   from the general ladder): decision rule in
   `thread_scaling_provenance_20260919.md` — a family plateaus at j* if
   points above improve <~10%/doubling; tiers whose caps exceed every
   family's j* are pruning candidates.
3. **Ryan-gated follow-ons**: R3+ re-runs (parked on the plateau verdict);
   optional f32full arms for R1–R2 (resubmit 13778533's launcher with
   `FGS_DAGTEAM_PRECISION=f32full` AFTER the f64 staircases certify);
   optional zen3 BLAS A/B rider at the R4 champion (qos=test, AB trials j16
   BLAS 4 vs 1) — expected null, only if Ryan wants certainty.
4. **Standing Ryan-gated ledger** (unchanged): origin pushes (both repos +
   v21/v22/v23 + the two new campaign tag sets), notebook entries owed (v21,
   diagnostics ladder, v22 chunked, NUMA, v23 promotion; + these two
   campaigns once harvested — offer, don't write), WeakKeyDict/warmstart
   one-line fix (`_publish_block_gs_status!`, from `7fbd68a`), 3 RECENT p018
   runs awaiting archive approval, hpc-storage archive-pass report
   collection. Also missing: a `:dagteam` unit test in
   `test/runtests_unit_solver.jl` (lex/colored/chunked covered, dagteam not).

## Traps (unchanged from 20260919g, plus new ones)

- Tolerances are per order AND per environment (mesh, rung, threads, BLAS,
  machine) — never carry 3.43e-7 anywhere; every point re-staircases.
- f32full is certified at R4/zen3 ONLY; off-R4 first rung is dagteam+f64;
  the independent evaluator (not the internal residual) is the accuracy
  authority.
- Placement is load-bearing (socket-membind collapsed dagteam to 1.09×);
  "node 0" ≠ "socket 0" on NPS4 EPYC (socket 0 = nodes 0-3).
- NEW: `HARDWARE_TAG` is part of the tuning-trace hard guard — the
  socket0-ilv tags exist so these rows/traces never mix with historical
  both-socket phase-1/2 data; don't "clean them up".
- NEW: phase-1 knob CSVs have NO sweep_order column; `stage3_winner`/
  `staircase_for` select on rung+knobs only → per-run `BENCH_CASE_ROOT`
  isolation is the correctness boundary. Never point two sweep orders at one
  case root.
- NEW: julia block-buffers stdout to files — a killed job's empty log means
  lost buffer, not "did nothing". zsh does NOT word-split unquoted `$var`
  env lists; `echo ===` breaks zsh too.
- Local runs ≤4 threads; `ssh orc` needs a live ControlMaster socket; Slurm
  CLI needs a login shell (`bash -lc` + `source /etc/profile`).
- RigidWakeBody shedding from CONSTRUCTED cells if touching rotor drivers.

## House rules (binding, unchanged)

HPC submission Ryan-gated; monitoring via `hpc-monitor`; storage via
`hpc-storage` (400 G cap); notebook writes Ryan-gated (offer, don't write);
dated status/provenance files in BRAINSTORM/021; each campaign in its own
worktree, never a shared live checkout while jobs run.
