# Reset prompt: 021 FGS scalability Stage 2 — babysit + harvest (2026-09-23, supersedes fgs_scalability_stage1_reset_prompt_20260923.md)

Copy everything below the line into a fresh agent session started in
`~/Dropbox/research/projects/FLOWPanel.jl`.

---

## Where this stands

BRAINSTORM 021 is CLOSED AT PROMOTION; active follow-on = FGS scalability
diagnostic. Read `CLAUDE.md` first and `agent_policies/HPC.md` before HPC
work. Stage 1 (job 13858983) is COMPLETE and analyzed
(`fgs_scalability_stage1_results_20260923.md`): both effects reproduced and
fully localized to the dagteam near-field sweep (`nonself_product`
anti-scales 2.0→2.7→4.5 s over j=16/32/64); placement settled (socket-0
champion by +4.4 s); worker cap16 @ j=64 is the best known operating point
(3.41 s). Commits `a835265` (Stage-1 results) and `90f7452` (Stage-2 prep)
on `fastmultipole`, local only.

**Stage 2 (mechanism) SUBMITTED 2026-09-23: job 13875511** (`p021-fgs-stage2`,
m12 zen3 exclusive 128c/500G, **12 h wall**, PENDING at submission; test-only
estimated start 2026-09-24T15:41, likely earlier via backfill). Ryan
pre-approved this launch (2026-09-23), resume path included. Provenance +
pins: `fgs_scalability_stage2_provenance_20260923.md` (tag
`campaign/p021-fgs-stage2-20260923` in all three repos; rsync deployment to
`/home/rander39/campaigns/p021-fgs-stage2-20260923/`, ARCHIVER_SKIP — never
touch). Run dir:
`/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/fgs-stage2-13875511`.

Stage 2 measures (a) per-worker drain-loop aggregates — discriminator:
per-task busy inflating with j → bandwidth/locality; flat busy with growing
idle+lockmgmt shares → scheduling — and (b) a bounded-backoff idle-policy
A/B at j=64 (recovery isolates idle lock hammering). 17 stages + 4
conditional j=32 idle pairs (only if backoff beats spin at j=64;
`backoff_verdict.txt`).

## TASK — babysit job 13875511, then harvest + analyze

1. Monitor via `hpc-monitor` only. Judge by outputs (STATUS_* files, expect
   17 or 21, all ok; `COMPLETED` with `failed_count=N`), never sacct; task
   logs output-buffered. Individual FAILEDs are findings, not aborts.
2. On COMPLETED: rsync the run dir's CSVs/TOMLs/STATUS/verdict files locally
   (Stage-1 pattern: ~5 MB) and run
   `julia --project=benchmark -t 1 benchmark/fgs_stage2_analysis.jl <local_copy>`
   — FIRST real execution; expect to debug against the actual CSVs (Stage 1's
   script needed 2 fixes on first contact). Sync `analysis/` back to the run
   dir. Read `backoff_verdict.txt`.
3. **Gates (binding):** attribute aggregates only if the aggdiag-vs-anchor
   overhead gate passes (≤5% each rung). Worker durations are shares of
   team-time, never summed as elapsed. Plateau vs regression stay separate
   conclusions. An empty ready queue does not by itself distinguish
   dependency starvation from all-work-already-running (Stage-3 territory).
4. Report the discriminator verdict + backoff A/B to Ryan; write a dated
   results note in BRAINSTORM/021; offer (don't write) a notebook entry.
   **Stage 3 / any production adoption of backoff or cap16 is Ryan-gated.**
5. Resume path if the job dies: resubmit with `RESUME_FROM_JOB_ID=13875511`
   in the `--export` list (same COLD_PROJECT/CAMPAIGN_PINS/COLD_DATA_ROOT,
   from the deploy tree top level); STATUS_*=ok stages skip.

## Owed (carried)

- **032 ledger mirror** of the p032-rootomit archive lines — offered, awaiting
  Ryan (`archiver_campaigns_support_20260922.md`).
- **Origin push** of branches+tags in all three repos once Ryan re-auths
  GitHub (`gh auth login -h github.com`). Now also includes Stage-2 commits
  `90f7452` (FLOWPanel), `053c8de7` (FastMultipole) and the
  `campaign/p021-fgs-stage2-20260923` tags.
- Archiver **T5 pre-existing failure** (exit 9 vs 8) awaits Ryan's ruling.
- Stage-1 run dir `fgs-stage1-13858983` and the timed-out
  `p021-r4-thread-scaling` dir are archive-eligible once quiet (storage flow).
- Standing Ryan-gated ledger unchanged
  (`fgs_scalability_reset_prompt_20260922.md`).

## Traps (prior traps bind; key ones)

- Fixed-work rows: `solved=false`/`eligible=false` BY CONSTRUCTION — gate on
  certified-accepted + iterations==27 + 1e-8 repeat.
- `diag_*` = −1 marks uninstrumented rows; never pool instrumented with
  uninstrumented. New `diag_dagteam_*` aggregate columns are zero on
  pre-Stage-2 pins and for non-dagteam sweeps.
- `lockmgmt_ns` excludes idle-streak lock churn (folded into `idle_ns` by
  design); `busy_max/min_ns` are sums of per-outer-iteration extrema.
- The spin baseline binary now carries qhint maintenance — compare against
  THIS job's anchors, not Stage-1 absolute numbers; paired A/Bs are
  same-binary.
- `ssh orc` needs a live ControlMaster socket (`ssh orc -fN` if cold) AND
  `bash -lc "..."` for slurm/module commands. Local runs NEVER >4 threads.
  Never edit source while a job uses its deployment (manifest-verified).
- Run-dir names case-sensitive; judge runs by outputs, never sacct.
- `test/runtests_benchmark_cold.jl` fails on the laptop (BLAS pin,
  pre-existing/environmental) — not a regression signal.

## House rules (binding)

Monitoring `hpc-monitor`; harvesting `harvester`; storage `hpc-storage`;
notebook writes Ryan-gated (offer, don't write); dated status/provenance
files in BRAINSTORM/021. HPC submission approved for STAGE 2 job 13875511
only (its resume path included); Stage 3 and anything else Ryan-gated.
