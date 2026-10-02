# Reset prompt: FGS acceleration — harvest gate-2 rev b, decide, then implement (2026-09-19b)

Copy everything below the line into a fresh agent session started in
`~/Dropbox/research/projects/FLOWPanel.jl`.

---

## Where this stands

You are mid-way through implementing
`BRAINSTORM/021_rotor_hover_solver_benchmarks/fgs_acceleration_recommendation_20260918.md`
(authoritative spec — read it in full, plus companion
`thread_efficiency_top5_20260918.md`; where they disagree the recommendation
wins). The gate-2 HPC benchmark has now run TWICE: rev a (job 13773451) was
**invalidated by a numactl binding bug**; the corrected rev b (job
**13773494**) has **just completed cleanly, unharvested** (`== done ==` in
its slurm out). Also read `CLAUDE.md`, `agent_policies/WORKFLOW.md`,
`agent_policies/TESTING.md`, and `agent_policies/HPC.md` before the
corresponding work.

Goal and acceptance (unchanged): accelerate the prepared cold FGS body solve
at R4. Baseline = colored @ j16 = **10.116 s** (26 iters). Minimum
**≤ 6.744 s** (1.5×), design target **≤ 5.058 s** (2×), planning range
4.5–5.9 s. Accuracy gate: authoritative BC relative L2 ≤ 1e-6 via the
independent evaluator (`benchmark/fgs_cold_README.md` in FLOWPanel).

## What happened to rev a (read `fgs_acceleration_status_20260919.md`)

Every arm ran under `numactl --cpunodebind=0 --membind=0`. On the m12
EPYC 7763 (NPS4) that is ONE NUMA domain — 16 CPUs, 2 of 8 socket-0 memory
channels — not socket 0 (= nodes 0-3, 64 CPUs). Consequences: all t32 arms
oversubscribed (33 spin threads on 16 cores → ~25 s/sweep garbage),
bandwidth capped at ~27-33 GB/s actual, first-touch comparison nullified
(membind forced all pages to node 0), interleave controls void. qos=test was
cleared as a suspect (AllocCPUS=128, full node). Rev a DID establish:
handoff kill switch passes (1.81/2.18/2.92 µs at t4/8/16, budget <5.8 µs),
F32 convert-on-load kernel works (~1.7× stream reduction), structural
checks exact on cluster inputs.

## Rev b (job 13773494) — harvest this FIRST

Fix (commit **`b1c6c7af`** on dev branch, `fgs_acceleration_provenance_20260919b.md`
is the provenance): team arms bound `--cpunodebind=0-3 --membind=0-3` (true
socket 0); serial 1-core baselines keep node-0 binding (correct for 1
thread); interleave controls `--interleave=0-3 --cpunodebind=0-3` at t16 AND
t32; per-arm `numastat -p` snapshot 60 s in (`numastat_<arm>.txt`); dag arm
priced from corrected serial B + t4 handoff h — its sim has NO
shared-bandwidth cap (W workers at B each), so sanity-check aggregate W×B
against the measured socket ceiling before believing makespans.

Dead first attempt 13773493: FAILED in 1 s — submitted WITHOUT the
CENSUS/EDGES env overrides (driver default path doesn't exist on orc). No
results, no code difference. 13773494 was submitted with
`CENSUS=/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/diag-v15-13694724/j64-b1/results/gemv_census.csv`
and matching `EDGES=.../dependency_edges.csv`, md5-verified at submission
(`927a9b9988305febf52d448f0127a318` / `40f7f2d22bc096ae85769777e2e2e020`).

Unharvested outputs on orc:

- `/home/rander39/wt-p021-fgs-gate2/slurm-p021-fgs-replay-gate2-13773494.{out,err}`
- `/home/rander39/wt-p021-fgs-gate2/benchmark/replay_gate2_13773494/`
  (topology.txt, per-arm logs, numastat snapshots, dag.log)

Harvest to
`BRAINSTORM/021_rotor_hover_solver_benchmarks/fgs_r4_followup_evidence_20260914/replay-gate2-13773494/`
(rev a's harvest sits beside it as `replay-gate2-13773451/`). Log format:
each arm log has `useful BW (min)` (F64-equivalent = 231.891 GB / 81-sweep
stream) and `projected 81-sweep coefficient-stream time`.

## IMMEDIATE TASK 1: evaluate decision rules and report go/no-go

- **Handoff kill switch**: ladder avg < ~5.8 µs (86,508 handoffs/solve;
  10 µs costs 0.865 s). Check t32 is now sane (64 cores available).
- **Projected solve time**: T ≈ 3.0 s remainder + stream + H, stream =
  81 × best sweep min. Reject candidates that cannot reach 6.744 s with
  margin; design for ≤ 5.058 s. (231.891 GB F64 / 115.945 GB F32 per solve;
  minimum needs 61.9 GB/s F64-equiv or 62 F64-equiv for F32; halving needs
  112.7 / 56.3 actual. NUMA-diagnostic ceilings: 29.4 one core, 74
  serial-touch + parallel consume, 164.4 one socket affine.)
- **First-touch**: serial vs owner is now a REAL comparison (membind spans
  nodes 0-3); verify with the numastat snapshots; sub-page tiles can defeat
  ownership. Compare against interleave controls at t16/t32.
- **Source-major vs pull-DAG**: compare best rowpar arms against dag.log's
  priced schedule, remembering its missing bandwidth cap and the 2.849
  work/span bound; spec rule = promote split only if a complete inclusive
  schedule comparison at equal precision/placement wins with margin.
- Anomalies: if rowpar scales worse than the synthetic NUMA bench even with
  correct binding, diagnose before recommending.

Write `fgs_acceleration_status_20260919b.md` (or c if the name is taken) in
BRAINSTORM/021: handoff µs by team size, best useful BW per arm, numastat
placement verdicts, projected solve times, schedule pick, go/no-go.
**Report to Ryan and STOP before touching the production solver path.**

## What exists (all committed on the dev branch)

FastMultipole dev branch `p021-fgs-accel-20260918`, HEAD **`b1c6c7af`**
(chain: `c18e4b46` Manifest pin → `b0946c36` gates → `a6492f49` rev-a driver
→ `b1c6c7af` rev-b binding fix). Checkouts:

- local worktree `/private/tmp/fastmultipole-p021-fgs-accel-20260918`
- orc worktree `/home/rander39/wt-p021-fgs-gate2` at `b1c6c7af` (synced by
  pushing ref `p021-fgs-gate2b` to the orc clone via ssh — the https origin
  needs Ryan's token, origin push still Ryan-pending — then `git merge
  --ff-only p021-fgs-gate2b` in the worktree; repeat that dance for future
  syncs, a checked-out branch can't be pushed to directly)

Artifacts on the branch: gate-1 harness `test/fgs_rowpar_gate1_test.jl`
(shadow row-tiled executor vs production lex `gs_sweep!`; PASSING locally,
bitwise identical incl. transformed fixtures); gate-2 benchmark
`benchmark/fgs_sequence_replay.jl` (modes serial/rowpar/handoff/dag,
structural checks exact: 1,068 leaves, 2,862,850,032 B/sweep, 48,167 lower
edges, unit path 279, work/span 2.849); HPC driver
`benchmark/fgs_sequence_replay_orc.slurm.sh` (rev b).

**Gate-1 finding:** two-system `FastGaussSeidel((sysA,sysB),(sysA,sysB))`
fails to CONSTRUCT at `c18e4b46` (BoundsError in nonself-matrix fill) —
pre-existing; R4 is single-system so gates proceed, but production must
revisit or declare the limitation.

## TASK 2 (after Ryan's go): production implementation

In FastMultipole on the dev branch (NOT the dirty live checkouts — both
FLOWPanel `fastmultipole` and FastMultipole `flowpanel-20260817` carry
unrelated uncommitted work):

1. Persistent adaptive row-parallel execution of the existing source-major
   nonself cache (new `sweep_order` value or equivalent), consumer-aligned
   placement. Keep lex sequence, far-field refresh, RHS semantics
   (`+= old` then `-= new`, never a delta). No allocations/task creation in
   the leaf loop. Gate-1 shadow executor defines the schedule; gate-2 rowpar
   team is the reference coordination pattern (epoch/done atomics +
   GC.safepoint in spin loops).
2. Float32 nonself storage / Float64 arithmetic: explicit container change —
   `Matrices{TF}` (containers.jl:1092) couples coefficient and product/RHS
   types; `FastGaussSeidel` couples self/nonself precision. Convert on load,
   accumulate F64, never round strengths to F32. Fallback = selective
   per-block F64 retention, never a relaxed gate.
3. Gate order binding: correctness (gate-1 harness vs the REAL
   implementation, incl. multi-system) → numerical gate (independent F32
   certification, external evaluator) → end-to-end interleaved A/B vs
   unchanged champion, ≥1.5× accepted throughput, full campaign ceremony
   (tagged worktrees `campaign/<item>-<slug>-YYYYMMDD`, Manifest pins,
   provenance before submission).

## Code anchors (verified at c18e4b46)

FastMultipole `src/solve.jl`: `gs_sweep!` :1214 (lex branch :1277–1302),
`compute_nonself_products!` :923, `scatter_nonself_influence!` :948,
`nonself_influence_matrices` :143 (fill :256–321), `residual!` :1803
(**shared scratch — threading races without private scratch**),
`color_leaves` :1126 (symmetrized, NOT the pull graph), outer loop :1380+
(FLOWPanel sets `final_update=false`, `src/FLOWPanel_solver.jl:1814–1824`);
`src/containers.jl:1092` `Matrices{TF}`; `nearfield_cache.jl:438–478`
private-buffer parallel-assembly pattern.

## Known traps (each cost a prior agent time)

- **"node 0" ≠ "socket 0" on NPS4 EPYC** — the rev-a lesson; verify
  placement with numastat, never assume.
- sbatch of the replay driver NEEDS the CENSUS/EDGES env overrides (the
  rev-b-first-attempt lesson).
- Row-tiled kernels: mathematically equivalent to BLAS, not guaranteed
  bit-identical. Apple BLAS WAS bit-identical in gate 1; re-check on the
  HPC OpenBLAS build, else certify accuracy.
- Source-affine first touch is WRONG for row-worker consumption.
- Chunked v22 LOST end-to-end (14.264 s, 44 iters) despite faster sweeps;
  its 38.8 active threads measured a Jacobi-lagged schedule.
- Colored and lex have separately calibrated tolerances/orderings — never
  mix their iteration counts.
- The reverse flag repeats forward order — don't silently change it.
- Split design (if promoted): backward products must not overwrite the
  frozen upper accumulator while targets read it; initialize Ux^0 for warm
  starts.
- Sub-range BLAS gemv on strided views: fine; one strided dot per row: never.
- M2 local runs: E-cores drag balanced tiles at t4 — machinery-validation
  only.
- Judge runs by outputs, never sacct exit status (13773494 shows COMPLETED,
  but rev a also "completed" while invalid — read the logs).

## House rules (binding)

- Local runs ≤ 4 threads; full-scale timing only on HPC.
- HPC submission Ryan-gated (rev b was approved 2026-09-19 "go ahead"; new
  submissions need fresh approval). Monitoring via `hpc-monitor` subagent;
  `ssh orc` needs a live ControlMaster socket (else 2FA — stop and ask Ryan
  to run `ssh orc -fN`). Slurm CLI on orc needs
  `source /etc/profile; module load slurm`.
- Notebook writes and anything outward-facing are Ryan-gated. Ryan-pending
  ledger: origin pushes (merged branches, v21/v22 tags, AND
  `p021-fgs-accel-20260918`), 4 notebook entries, WeakKeyDict/warmstart fix.
- Dated status/provenance files in BRAINSTORM/021 per existing convention.
- The orc worktree `/home/rander39/wt-p021-fgs-gate2` is this thread's own;
  don't touch other campaigns' worktrees or the live clones.

## Context: scaling expectations (discussed with Ryan 2026-09-19)

Prediction on record, to be revisited when a larger rung exists: per-sweep
cost linear in N (~49 KB/DOF F64 nearfield cache, fixed leaf/MAC), handoff
and remainder constant fractions, so T(N) ≈ (N/58k) × T(R4) IF outer
iterations stay ~mesh-independent (second-kind operator argument; the
uncounted risk term — needs cold solves on an R5-type rung to measure).
Multi-system constructor bug becomes blocking for multi-rotor growth; leaf
size likely wants retuning upward at larger N (spec rank 4, report
separately).

## Suggested first moves

1. Read the recommendation + top-5 note end to end, then
   `fgs_acceleration_status_20260918.md`, `..._status_20260919.md`, and
   `..._provenance_20260919{,b}.md`.
2. Harvest 13773494 (scp slurm logs + results dir as above).
3. Tabulate against the decision rules; write the status file; report
   go/no-go + schedule pick to Ryan. Stop there until he answers.
