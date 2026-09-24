# Reset prompt: 021 :dagedge benchmark — babysit, harvest, analyze (2026-09-24, supersedes fgs_dagedge_benchmark_reset_prompt_20260924.md)

Copy everything below the line into a fresh agent session started in
`~/Dropbox/research/projects/FLOWPanel.jl`.

---

## Mission

Babysit, harvest, and analyze the 021 :dagedge benchmark campaign — **both
jobs are already SUBMITTED** (2026-09-24, the only two Ryan-pre-approved
submissions; resubmission of a failed/timed-out job via the recorded resume
path is in-scope, anything else needs a fresh Ryan go):

| job | run | state at handoff | run dir (under `/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/`) |
|---|---|---|---|
| **13879622** | performance (`RUN_MODE=perf`, uninstrumented) | PENDING on m12; test-only ETA 2026-09-24T11:59, backfill likely | `fgs-dagedge-perf-13879622` |
| **13879625** | profile (`RUN_MODE=profile`, instrumented) | PENDING, `--dependency=afterany:13879622` | `fgs-dagedge-profile-13879625` |

Read `CLAUDE.md`, `agent_policies/HPC.md`, top-level `BYU_ORC_AGENTS.md`
before any HPC work. Delegate: job status/log tails → `hpc-monitor`;
CSV/TOML scraping → `harvester`; archiving → `hpc-storage` (judge by
outputs, never sacct). `ssh orc` needs a live ControlMaster socket
(`ssh orc -fN` first if cold; never retry into 2FA).

## What's being tested (one paragraph)

`FGSSolver(sweep_order=:dagedge)` — edge-level partial pulls on a static
per-worker schedule — vs the Stage-2 champion dagteam+backoff (3.24 s/solve
@ j=64, R4 27×3 fixed cold work). Gate-0 predicts a 4.25× sweep bound
(θ=4KB), roughly flat in j; expected end-to-end ~1.7–1.8 s @ j64.
**Decision threshold: < ~3× sweep gain ⇒ scheduler overhead ate the prize**
— the profile run must then attribute it. Full context:
`fgs_dagedge_benchmark_provenance_20260924.md` (campaign record, binding
measurement caveats), `fgs_dagedge_design_20260924.md`,
`fgs_lshortening_gate0_20260924.md`,
`fgs_scalability_stage2_results_20260924.md`.

## Campaign facts (all committed; provenance file is authoritative)

- Tag `campaign/p021-fgs-dagedge-20260924` in all three repos: FLOWPanel
  `a418f8d` (harness), FastMultipole `90a60cc3` (:dagedge impl), FLOWVPM
  `8d4a3b4`. Deployed rsync-mode (git archive, manifests VERIFIED) to
  `/home/rander39/campaigns/p021-fgs-dagedge-20260924/` (ARCHIVER_SKIP;
  never edit these trees while either job is queued/running). Env at
  `.../env`, pins at `.../pins.toml`.
- Launcher `benchmark/run_r4_fgs_dagedge.slurm.sh`; driver
  `benchmark/fgs_r4_dagedge.jl` (arms batched per (j, placement) process;
  cold = zero-initial-guess solves).
- Perf arms: `primary-p{1,2,3}` (paired dagteam vs dagedge:4096 @ j64
  backoff, in-process order alternates by block), `ladder-j{16,32}` (both
  executors), `theta-j64` (dagedge θ ∈ {0,4096,16384}, one process).
- Profile arms: `prof-j64-b{1,2}`, `prof-j{16,32}` — same pairs,
  `DAGEDGE_DIAG=1`.
- Per stage dir: `results/dagedge_solves.csv` (one row per solve; arm,
  sweep_order, theta, plan_* and sim_* schedule metadata, diag_* columns
  with −1 sentinels), `results/dagedge_summary.toml` (per-arm medians +
  schedule metadata), `residual_history_arm*.csv`, `STATUS_<name>` files,
  `COMPLETED` sentinel with failed_count.
- Resume a failed/timed-out job: resubmit the same sbatch line from the
  deployed FLOWPanel tree top level with `RESUME_FROM_JOB_ID=<old id>`
  added to `--export` (STATUS_*=ok stages skip):
  `sbatch -p m12 --export=ALL,RUN_MODE=<perf|profile>,RESUME_FROM_JOB_ID=<id>,COLD_PROJECT=/home/rander39/campaigns/p021-fgs-dagedge-20260924/env,CAMPAIGN_PINS=/home/rander39/campaigns/p021-fgs-dagedge-20260924/pins.toml,COLD_DATA_ROOT=/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910 benchmark/run_r4_fgs_dagedge.slurm.sh`

## Harvest & analysis (write `fgs_dagedge_benchmark_results_20260924.md` — or actual date — in BRAINSTORM/021)

Primary questions, in order:

1. **Ranking (perf run ONLY):** paired dagedge-vs-dagteam medians per
   primary block @ j64; did 3.24 → ~1.7–1.8 s materialize? Report paired
   deltas per block (same process ⇒ same placement/compile state), then
   pooled median.
2. **Flat-in-j:** dagedge solve time across j16/32/64 vs dagteam's plateau
   (backoff ladder baseline 4.155/3.478/3.262 s).
3. **θ probe:** dagedge @ θ=0/4KB/16KB @ j64 (in-process comparison; the
   primary pairs are the ranking-grade numbers).
4. **Profile attribution (never pooled with perf rows):** busy_lower vs
   idle vs wait/reduce shares of team-time (durations overlap — never sum
   as elapsed); does idle collapse vs dagteam @ j64? Per-task busy
   (busy_lower / n_lower) vs the R2-extrapolated ~0.3 µs/task at ~31.6k
   tasks/sweep; busy_max/min imbalance (static-schedule skew is the known
   risk; recorded fallback = work-stealing deques; next optimization = NUMA
   first-touch repack of split Lmat blocks). lockmgmt must be identically 0
   for dagedge.
5. Sanity: schedule metadata (plan_ntasks, sim_edge_L ≈ 68–69 MB/sweep at
   θ=4KB @ R4 expected from gate-0) — build-time simulation outputs, not
   measurements. Cross-executor deltas (tripwire 1e-5; R2 measured 3.4e-7).

Binding caveats: fixed-work gates (27 iterations, certified evaluator, 1e-8
repeat); diag −1 sentinels; instrumentation shifts timing (Stage 2: 8.5%);
failures are findings (STATUS_*=FAILED, job continues).

## Traps

- Judge runs by outputs (STATUS_*/COMPLETED), never sacct.
- Never edit the campaign trees or live checkouts while jobs are queued or
  running; laptop FLOWPanel HEAD `069b173` already matches the deployed
  content plus the filled provenance.
- The θ-probe arms share one process — construction-order effects possible;
  don't rank from them.
- These runs write no VTK; storage pressure is unlikely, but 018 NT-ladder
  GPU jobs (13878882, 13879081–94) are running/queued concurrently — if
  disk alarms fire, that's their VTK, launch `hpc-storage`.
- Local runs ≤4 threads; laptop `runtests_benchmark_cold.jl` failure is
  pre-existing/environmental (BLAS pin) — don't chase.

## Standing Ryan gates (surface, don't act)

- backoff as :dagteam production default (recommended); Stage-3
  cancellation (recommended). If dagedge wins, a NEW recommendation
  (dagedge as production default + any follow-up ladder) goes to Ryan — no
  further submissions are pre-approved.
- Notebook entry for Stages 1+2 + gate-0 + dagedge prototype + this
  campaign: OFFER once harvested, don't write.
- Origin push of branches + tags (`campaign/p021-fgs-stage2-20260923`,
  `campaign/p021-fgs-dagedge-20260924`, all three repos) after Ryan runs
  `gh auth login -h github.com`.

## House rules (binding)

Dated status/results files in BRAINSTORM/021; notebook Ryan-gated;
delegation per CLAUDE.md; monitoring cadence ≥60 s between scheduler
queries; no new HPC submissions beyond the recorded resume path.
