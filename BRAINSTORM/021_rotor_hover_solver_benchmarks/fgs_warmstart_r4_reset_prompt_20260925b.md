# Reset prompt: 021 warm-start R4 campaign — finish babysit, harvest + results (2026-09-25b; supersedes fgs_warmstart_r4_reset_prompt_20260925.md)

Copy everything below the line into a fresh agent session started in
`~/Dropbox/research/projects/FLOWPanel.jl`.

---

## Mission

The 021 warm-start R4 campaign is DEPLOYED and NEARLY DONE (Ryan's launch
go was given 2026-09-24). Remaining work, in order:

1. **Finish babysitting Job 1 = 13890195** on orc. As of 2026-09-25
   ~20:15 UTC it was RUNNING and healthy with **10/14 legs complete**
   (see State). Arm a background squeue poll (~30 min interval) that
   fires when the job leaves the queue; on wake, verify the terminal
   state independently (a dead ssh socket also ends a poll loop).
   Job 2 = 13889502 is COMPLETE and healthy — nothing to do until harvest.
2. When Job 1 is terminal: judge each leg by its STATUS sentinel AND the
   CSV `solved` column (never sacct exit codes — and note the driver now
   hard-fails a leg on any `solved=false` row).
3. **Harvest BOTH jobs**: `benchmark/fgs_r4_warmstart_harvest.jl` run from
   the deploy tree `/home/rander39/campaigns/p021-fgs-warmstart-20260924/FLOWPanel.jl`
   (harvest is leg-aware: Window A = non-restarted steps 1..36; Window B =
   restarted rows, global step = restart_step + local step).
4. **Write `BRAINSTORM/021_rotor_hover_solver_benchmarks/fgs_warmstart_r4_results_20260925.md`**:
   per-arm Window A/B tables (niter + time-to-target mean±spread/median),
   setup costs in their own columns (never amortized), per-step cost traces
   vs step index, solution-agreement report (cross-arm strength deltas from
   the `*_strength_snapshots.bin`; the known ~2e-3 FGS-vs-Krylov wake-on
   fixed-point discrepancy is REPORTED, not chased), CT traces, and the
   compete/no-compete recommendation. **Ryan rules** — present medians,
   spreads, and cumulative-cost crossover curves, don't pre-judge.
   REQUIRED reporting note: solver warm-start histories are not serialized
   in checkpoints, so each restarted warm leg's first (order+1) steps are
   effectively cold — inside Window B's transient-included stats by design.
5. Also fold Job 2 into the deliverable: re-issue the 2026-09-22 R4 cold
   table with an "FGS-dagteam+backoff (new default)" column
   (run dirs `fgs-cold-newdefault-j{1,8,16,32,64}-13889502`; j64 cold min
   3.24 s, exactly the expected 3.2–3.3 s).
6. Offer (don't write) the notebook entry: FGS Stages 1+2, gate-0, dagedge
   campaign + verdict, default adoption, warm-start campaign + winB restart
   bug. Ask Ryan for verbosity per topic.

Required reads first: `CLAUDE.md`, `agent_policies/HPC.md`, cluster
`BYU_ORC_AGENTS.md`, and
`fgs_warmstart_r4_provenance_20260924.md` (provenance of record — fully
updated through the 13890195 submission). `ssh orc` needs a live
ControlMaster socket (`ssh orc -fN` if cold; never retry into 2FA).
Local runs ≤4 threads. Slurm needs a LOGIN shell: `ssh orc 'bash -lc "…"'`.

## State (2026-09-25, ~20:15 UTC)

- **Job 1 = 13890195 RUNNING on m12-3-12**, elapsed 15:00 of 36 h,
  healthy (~2–2.5 h/leg, steady). **10 legs "ok"**: fgs_cold_ckpt,
  ilu_nfcache_cold_ckpt, fgs_prev_winA, fgs_proj1_winA, fgs_proj2_winA,
  ilu_nfcache_prev_winA, ilu_nfcache_proj1_winA, fgs_cold_winB,
  fgs_prev_winB, fgs_proj1_winB. **Running**: fgs_proj2_winB. Remaining
  after it: the ilu_nfcache winB family (cold/prev/proj1) → projected
  finish ~2026-09-26 02:00–06:00 UTC, well inside walltime.
- **The winB restart fix (`0aaceaf`) is holding in production**: three
  winB legs (which exercise the restarted-solve path that used to die)
  have completed "ok".
- **Open oddity (non-blocking)**: monitor found 0 `.vtp` files in the
  expected ckpt VTK dirs `/home/rander39/projects/FLOWPanel.jl/data/fgs_wsr4_R4_ckpt_{fgs,ilu}/`,
  yet both ckpt legs completed "ok" AND winB restarts from those
  checkpoints succeed — the checkpoints evidently live under a different
  path (likely inside the run dir). Confirm actual ckpt/VTK location at
  harvest time before applying archive-first storage policy; do not
  treat the empty dirs as a failure.
- Run dir (sentinels, leg logs, CSVs):
  `/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/fgs-wsr4-13890195/`.
- **Job 2 = 13889502 COMPLETE + healthy** (all 5 rungs COMPLETED +
  certified; j64 cold min 3.24 s).
- Everything else (fix verification, tag `campaign/p021-fgs-warmstart-20260924`
  recut at `0aaceaf`, three deployment defects fixed, dead run dirs
  13889222/13889501) is unchanged from the 20260925 reset prompt and
  fully recorded in `fgs_warmstart_r4_provenance_20260924.md`.

## Traps

- Task/leg logs are OUTPUT-BUFFERED — judge liveness by sstat CPU, STATUS
  sentinels, and CSV outputs, never by log tails. "fatal: not a git
  repository" noise in leg logs is expected (archive-export mode; CSV
  `commit` columns read "unknown"; provenance = pins + manifests).
- sacct FAILED ≠ real failure and COMPLETED ≠ success — judge by outputs.
- Cold = zero-initial-guess in the SAME process (Ryan 2026-09-23).
  Pre-sim setup costs excluded from per-step comparisons, reported in the
  setup columns.
- A `RESUME_FROM_JOB_ID` path exists in the launcher (landed legs skip by
  STATUS) if Job 1 dies mid-sequence — resubmit with that env var set to
  13890195 rather than re-marching finished legs.
- If a winB leg fails, its family ckpt is on disk and restartable —
  debug locally at R1/NT=4 first (smoke harness:
  `benchmark/run_r4_fgs_warmstart_smoke.sh`, STAGES env selects legs).
- 018 NT-ladder jobs run concurrently — disk alarms are their VTK; launch
  hpc-storage, don't touch their queue. Ckpt output is new disk load;
  archive-first policy applies AFTER harvest.
- Local laptop `runtests_benchmark_cold.jl` failure is pre-existing (BLAS
  pin) — don't chase.

## Standing gates (surface, don't act)

- Notebook entry (see mission step 6): offer once, verbosity per Ryan.
- Origin pushes: branches + tags `campaign/p021-fgs-stage2-20260923`,
  `campaign/p021-fgs-dagedge-20260924`, `campaign/p021-fgs-warmstart-20260924`
  (three repos) after `gh auth login -h github.com`.
- dagedge harvest (13879622/13879625, both out of queue) is a SEPARATE
  task (`fgs_dagedge_benchmark_provenance_20260924.md`).
- Ryan may still veto the four submission defaults used (champion cold
  tolerance 3.4309419310610173e-7; 36 h walltime; `ilu_nfcache_proj1`
  kept; sequential legs on one exclusive node).
