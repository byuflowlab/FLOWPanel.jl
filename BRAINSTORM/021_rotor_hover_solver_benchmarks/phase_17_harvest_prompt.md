# BRAINSTORM 021 — Phase 17: harvest the Phase 1 + Phase 2 HPC campaign

**Written 2026-08-24 at launch, for a fresh agent picking this up after the
jobs land. Read this file and `ledger.md`'s "Phase 16" section; you should not
need to read anything else to start.**

Your job is to HARVEST and REPORT. Nothing here authorizes an `sbatch`.

---

## 1. What is in flight

Launched 2026-08-24 ~19:15 from `~/flowpanel-021/FLOWPanel.jl` on `orc`.
All 13 jobs are `physics2` / `standby` / `--constraint=zen3 --exclusive
--mem=500G --requeue`, 64 threads.

| job id | name | what it does | walltime |
| --- | --- | --- | --- |
| 13469147 | p2-tune-R1 | Phase 2 real-solve tuning, budgets 0:16:128:500 | 12 h |
| 13469148 | p2-tune-R2 | " | 12 h |
| 13469149 | p2-tune-R3 | " | 24 h |
| 13469150 | p2-tune-R4 | " | 2 d |
| 13469151 | p2-tune-R5 | " | 3 d |
| 13469152 | p1-agree-R1 | Phase 1 agreement, 7 configs incl. dense backslash | 12 h |
| 13469153 | p1-agree-R2 | " | 12 h |
| 13469154 | p1-agree-R3 | " | 24 h |
| 13469155 | p1-agree-R4 | " | 24 h |
| 13469156 | p1-agree-R5 | 5 configs, dense dropped (8N^2 = 87 GiB) | 2 d |
| 13469157 | p1-agree-R6 | " (335 GiB) | 2 d |
| 13469158 | p1-agree-R7 | " (1.31 TiB) | 3 d |
| 13469159 | p1-agree-R1-direct | `RESIDUAL_BACKEND=direct` O(N^2) BC cross-check, configs `krylov_gmres,fgs`, **chained `afterany:13469152`** | 12 h |

At launch: all 5 tuning jobs RUNNING, agreement R1-R3 RUNNING, R4 pending on
Resources, R5-R7 on Priority, direct job held on Dependency.

**Preemption is verified to requeue automatically** — `PreemptType=preempt/qos`,
`PreemptMode=REQUEUE` at both cluster and partition level, `Requeue=1` on the
jobs. A requeued job keeps its ID and reruns the batch script from the top;
`Restarts=` in `scontrol show job` tells you how often it happened.
`GraceTime=0` and `PreemptExemptTime=0`, so kills are immediate with no flush
window — which is why the tuner checkpoints per candidate rather than on a
signal.

## 2. Where the output goes

```
~/flowpanel-021/FLOWPanel.jl/benchmark/results/
  phase2/multi/<RUNG>/tune_phase2.csv            <- Step A, one row per budget
  phase2/multi/<RUNG>/tune_trace_<RUNG>_b<B>.csv <- tuning checkpoint, not a result
  phase1/multi/<RUNG>/agreement.csv              <- Step B, one row per config
  phase1/multi/<RUNG>/agreement_spread.csv       <- ensemble spread, per rung
  phase1/multi/<RUNG>/agreement_x_<RUNG>_<cfg>.bin  <- persisted solutions
```

One job writes one file — that is Ryan's ruling of 2026-08-24 and the whole
reason for the per-rung directories. Do not "tidy" them into a shared file.

## 3. How to harvest

```bash
cd ~/flowpanel-021/FLOWPanel.jl
julia --project benchmark/p021_merge.jl
RUNGS=R1:R2:R3:R4:R5 BUDGETS=0.0:16.0:128.0:500.0 STRICT=1 \
  julia --project benchmark/p021_merge.jl
```

Writes `tune_phase2_merged.csv`, `agreement_merged.csv`,
`agreement_spread_merged.csv` at the `<mode>` level. It ASSERTS one row per
`(rung, budget)` and per `(rung, config)` and names anything missing or
duplicated. A missing row is a detectable absence — report it, do not paper
over it. The merged files are DERIVED: regenerate, never append.

Per repo policy, route the log tailing / status polling to the `hpc-monitor`
subagent and the CSV scraping to `harvester`. Keep the physics reasoning and
the conclusions inline.

## 4. What to check in the harvested tables — the traps

1. **`tune_timed_out`**. A `true` means the descent was cut short by
   `TUNE_MAX_SECONDS` and the knobs are BEST-SO-FAR, not a minimum. Do not
   publish a timed-out row as a tuned optimum. Note the subtlety recorded in
   the ledger: a preempted-and-resumed descent gets a FRESH `max_seconds` on
   each requeue, so this flag is the only thing separating "converged" from
   "ran out of clock", and a heavily-preempted rung can be MORE converged than
   an uninterrupted one. The descent is still bounded by `max_iters=20`
   accepted moves, so it cannot run away.
2. **`resumed_from_trace n_replayed=N`** in the `notes` column marks a row whose
   descent was resumed. The replayed timings came from the pre-preemption node.
   All nodes are spec-identical pinned zen3 and the winner is independently
   re-timed, so this affects which path the descent walked, never a published
   number — but say so if you quote such a row.
3. **`cache_capped=true`** means the winner was stopped by the memory budget,
   not by cost. Never read a capped winner as an unconstrained optimum.
4. **`bc_certified` / `bc_rel_l2`**. This is the error metric, not a reference
   comparison. An uncertified row is not a result.
5. **`DUPLICATE of previous budget`** in `notes` means that budget produced the
   identical configuration as the one below it and adds no new operating point.
   Do not present it as independent evidence on the memory curve.
6. **`krylov_ilu` at R7** may not fit (~23 GiB of pattern entries). That is a
   RESULT, not a failure — record it as such.

## 5. Failure modes and what they mean

- **Trace provenance error** (`trace provenance disagrees with this run`) — a
  job was relaunched with different `TUNE_REPS`, `TUNE_ABANDON_FACTOR` or
  hardware than the trace was measured under. This is a deliberate guard, not a
  bug: a memoized `t` is a timing. Either restore the original setting or move
  the trace aside to start that descent fresh. Do NOT delete it silently.
- **Schema guard** (`has a pre-W2/W3 schema; move it aside`) — an old CSV is in
  the way. Move aside with a dated `.bak`, never delete.
- **Below-floor row** — a budget cheaper than the rung's minimum uncached
  configuration writes an explanatory row rather than failing. That is the
  intended answer for a 16 GiB laptop at R6/R7.
- **A rung produced nothing** — check `Restarts=` and the `slurm-*.err` in the
  repo root. Under standby, repeated preemption is expected; the tuner resumes
  per candidate, but **Phase 1 agreement has NO candidate-level resume** and
  loses the config in flight. At R6/R7 one config is many hours, so that is the
  real exposure if those rungs thrash.

## 6. Hazards specific to re-running anything

**Per-rung directories isolate RUNGS, not REPEATS of a rung.** Launching a
second `p1-agree-R3` while the first is live puts two writers back on one NFS
file — the exact 2026-08-18 clobber the campaign was restructured to remove.
The same applies to any `RESIDUAL_BACKEND=direct` re-run, which is why 13469159
is chained rather than concurrent: `RESIDUAL_BACKEND` is a CSV COLUMN, not a
path component. If you need to re-run a rung, let the first finish or cancel it
first.

## 7. Quarantined — never mix into a result table

Local Mac smoke artifacts and superseded rows, all deliberately kept:

- cluster: `phase1/multi/R1/tune_superseded_13447582.csv.bak` (the mis-tuned
  leaf=43 apply-proxy row), `phase2/multi/tune_trace_R1_b0.0_login01depcheck.csv.bak`
- local: `phase2/multi/tune_phase2_localsmoke{2,3,4}.csv`,
  `tune_phase2_localsmoke5_traceresume.csv`,
  `tune_trace_R1_b0.0_localsmoke5.csv`, `phase2_w3_localsmoke.csv`,
  `agreement{,_spread}_localsmoke.csv`, `agreement_pre20260824.csv.bak`,
  `tune_cached.csv`

**No local timing is evidence.** Ryan's Mac scatters 22-39% at fixed knobs
against a ~15% effect. Local runs verify control flow only.

## 8. Still Ryan's call — do not decide

- **Every `sbatch`.** Including Step C (`p2_table.sh`), which is gated on that
  rung's Step A landing and needs the pre-W3 `phase2.csv` moved aside first.
- **R6/R7 tuning objective.** Real-solve tuning is affordable to R5; R6 is
  8-16 h per budget and R7 does not fit. Deferred until R1-R5 land — the first
  question to answer with the harvest is whether the real-solve winner AGREES
  with the apply proxy at R1-R5, because that decides whether R6/R7 can be
  tuned by proxy at all.
- **`TUNE_REPS` at R5+.** Currently the driver default of 2 (the wrapper no
  longer forces 5). Changing it must happen BEFORE a resume — it is in the
  trace provenance guard.
- Any notebook entry. Offer, don't write.

## 9. The question the campaign exists to answer

Phase 1: do the solvers agree on the SOLUTION, reference-free, up to R7?
Phase 2: given a machine with X GiB, which solver is fastest, and does the
near-field cache pay for itself? Report both as .md tables with absolute dates.
