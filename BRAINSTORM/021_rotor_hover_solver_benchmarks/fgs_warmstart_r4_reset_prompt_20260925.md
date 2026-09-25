# Reset prompt: 021 warm-start R4 campaign — babysit Job 1, then harvest + results (2026-09-25; supersedes fgs_warmstart_r4_reset_prompt_20260924b.md)

Copy everything below the line into a fresh agent session started in
`~/Dropbox/research/projects/FLOWPanel.jl`.

---

## Mission

The 021 warm-start R4 campaign is FIXED, DEPLOYED, and RUNNING (Ryan's go
was given 2026-09-24: "permission to launch on hpc once bugs are fixed").
Remaining work, in order:

1. **Babysit Job 1 = 13890195** on orc (hourly; hpc-monitor or quick
   squeue). Job 2 = 13889502 is COMPLETE and healthy.
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

## State (2026-09-25, ~07:00 UTC)

- **The winB restart blocker is ROOT-CAUSED, FIXED, VERIFIED** (commit
  `0aaceaf` on `fastmultipole`): `simulate_warmstart!` replayed rigid
  kinematics but never mirrored them into persistent solver FMM state via
  `transform_body_solvers!` — FGS trees and the primed persistent Krylov
  plan+nfcache stayed at construction pose (wrong operator at the first
  restarted solve; both families died; NOT NT-dependent). Fix accumulates
  the net rigid delta across the whole replay and applies it once after
  normals/control-points refresh. Forward NT=4 control was clean; fixed
  winB reproduces the forward checksum to ~1e-13. Attribution: fix + old
  f32 ckpt also passes ⇒ f32 was replay-exactness only. Leg smoke 7/7
  PASS; suites solver 513/513 + warmstart all PASS. Also fixed: driver
  forces `FLOWPANEL_PARTICLE_PRECISION=f64`; driver fails legs on
  `solved=false`; `_publish_block_gs_status!` guarded for immutable
  solvers (pre-existing WeakKeyDict crash). NO unit pin for the defect-4
  mechanism (warmstart suite only uses noop/Backslash) — the leg smoke is
  the regression guard.
- **Tag** `campaign/p021-fgs-warmstart-20260924` recut at `0aaceaf`;
  orc deploy rsync-refreshed, manifest VERIFIED, `pins.toml` updated.
  Local-only (origin push still owed after `gh auth login`).
- **Three deployment defects found at submission and fixed** (full detail
  in the provenance "Staged submissions" section):
  1. deploy `data` symlink was NESTED inside the extracted real `data/`
     dir (archive extraction raced the ln -s) → Das preflight refused;
     now a real symlink → `/home/rander39/projects/FLOWPanel.jl/data`;
     the content manifest now EXCLUDES `data/` (live-mutable by design —
     runs append to `rotor_hover_pressure_comparison.metadata.toml`).
  2. campaign env Manifest had been resolved under juliaup default 1.12.7
     vs the launchers' `module load julia/1.11.7` → rebuilt from the
     dagedge 1.11.7 env, dev-repointed at the warmstart deploy trees.
  3. env lacked `VSPGeom` (fixture example's first import; dagedge chain
     never loads it) → added v0.6.6; full driver import chain audited.
- **Job 2 = 13889502 COMPLETE + healthy** (all 5 rungs COMPLETED +
  certified; j64 cold min 3.24 s).
- **Job 1 = 13890195 RUNNING on m12-3-12** since ~05:15 UTC, healthy at
  last check (~06:50): first leg `fgs_cold_ckpt` ~58/108 steps at
  ~55 s/step, ~19-core avg CPU, 29.6 GB RSS. Projection ~12–17 h total
  for all 14 sequential legs (walltime 36 h — fits). Leg order: 2 ckpt
  (108 steps each, SAVE_VTK) → 5 winA (36) → 7 winB (36, gated on family
  ckpt STATUS). Run dir (sentinels, leg logs, CSV):
  `/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/fgs-wsr4-13890195/`.
  Ckpt VTK (step counter = .vtp count):
  `/home/rander39/projects/FLOWPanel.jl/data/fgs_wsr4_R4_ckpt_{fgs,ilu}/`.
  Failed earlier submissions 13889222 (preflight) and 13889501 (VSPGeom)
  left dead run dirs — ignore them; their forensics are in the provenance.

## Traps

- Task/leg logs are OUTPUT-BUFFERED — judge liveness by sstat CPU, STATUS
  sentinels, and VTK/CSV outputs, never by log tails. "fatal: not a git
  repository" noise in leg logs is expected (archive-export mode; CSV
  `commit` columns read "unknown"; provenance = pins + manifests).
- sacct FAILED ≠ real failure and COMPLETED ≠ success — judge by outputs.
- Cold = zero-initial-guess in the SAME process (Ryan 2026-09-23).
  Pre-sim setup costs excluded from per-step comparisons, reported in the
  setup columns.
- A `RESUME_FROM_JOB_ID` path exists in the launcher (landed legs skip by
  STATUS) if Job 1 dies mid-sequence — resubmit with that env var set to
  13890195 rather than re-marching finished legs.
- If a winB leg fails, its family ckpt VTK is on disk and restartable —
  debug locally at R1/NT=4 first (smoke harness:
  `benchmark/run_r4_fgs_warmstart_smoke.sh`, STAGES env selects legs).
- 018 NT-ladder jobs run concurrently — disk alarms are their VTK; launch
  hpc-storage, don't touch their queue. Two ckpt VTK trees (~R4, 108 steps
  each) are new disk load; archive-first policy applies AFTER harvest.
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
