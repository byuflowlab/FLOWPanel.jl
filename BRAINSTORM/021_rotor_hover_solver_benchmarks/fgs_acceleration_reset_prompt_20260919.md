# Reset prompt: FGS acceleration — harvest gate 2, decide, then implement (2026-09-19)

Copy everything below the line into a fresh agent session started in
`~/Dropbox/research/projects/FLOWPanel.jl`.

---

## Where this stands

You are mid-way through implementing
`BRAINSTORM/021_rotor_hover_solver_benchmarks/fgs_acceleration_recommendation_20260918.md`
(authoritative spec — read it in full, plus companion
`thread_efficiency_top5_20260918.md`; where they disagree the recommendation
wins). The previous session (status:
`fgs_acceleration_status_20260918.md`) built and validated both decision
gates; the gate-2 HPC benchmark has **just completed, unharvested**. Also
read `CLAUDE.md`, `agent_policies/WORKFLOW.md`, `agent_policies/TESTING.md`,
and `agent_policies/HPC.md` before the corresponding work.

Goal and acceptance (unchanged): accelerate the prepared cold FGS body solve
at R4. Baseline = colored @ j16 = **10.116 s** (26 iters). Minimum
**≤ 6.744 s** (1.5×), design target **≤ 5.058 s** (2×), planning range
4.5–5.9 s. Accuracy gate: authoritative BC relative L2 ≤ 1e-6 via the
independent evaluator (`benchmark/fgs_cold_README.md` in FLOWPanel).

## What exists (all committed on the dev branch)

FastMultipole dev branch `p021-fgs-accel-20260918`, HEAD **`a6492f49`**
(off Manifest pin `c18e4b46`). Checkouts:

- local worktree `/private/tmp/fastmultipole-p021-fgs-accel-20260918`
- orc worktree `/home/rander39/wt-p021-fgs-gate2` (same SHA; branch was
  pushed by direct `git push ssh://orc/home/rander39/projects/FastMultipole`
  because the https origin needs Ryan's token — origin push still pending
  with Ryan's other pushes)

Three artifacts on that branch:

1. **Gate-1 harness** `test/fgs_rowpar_gate1_test.jl` — shadow row-tiled
   executor of the source-major schedule vs production `:lexicographic`
   `gs_sweep!`. PASSING locally: **bitwise identical** at 1–4 tiles /
   4 threads, 3 sweeps, nonzero starts, plain + rigidly transformed
   (`transform_solver!`) fixtures. Runs standalone
   (`julia -t4 --project=. -e 'include("test/fgs_rowpar_gate1_test.jl")'`
   from the worktree).
   **Finding:** two-system `FastGaussSeidel((sysA,sysB),(sysA,sysB))` fails
   to CONSTRUCT at `c18e4b46` (BoundsError in the nonself-matrix fill,
   `reshape(view(...),524,76) at [523:526, 40:76]`) — pre-existing; R4 is
   single-system so gates proceed, but the production implementation must
   revisit or declare the limitation.
2. **Gate-2 benchmark** `benchmark/fgs_sequence_replay.jl` — replays the
   exact R4 sequence from the saved census/edges with production RHS
   semantics (`+= old` then `-= new`, never a delta). Structural checks all
   exact vs spec: 1,068 leaves, 2,862,850,032 B/sweep, m_j = Σ n_i per
   source, 48,167 lower edges, unit path 279, byte work/span 2.849. Modes:
   `serial`, `rowpar` (persistent spin team + coordinator, adaptive
   `--small-bytes` serial policy, `--first-touch serial|owner`, `--f32` =
   Float32 storage/Float64 accumulate), `handoff` (zero-work sync
   microbench), `dag` (pull-DAG structural + priced critical-path
   list-schedule sim with backward-filler model; `--dag-bandwidth GB/s`,
   `--dag-handoff-us`). "Useful BW" is reported in F64-equivalent bytes so
   F32 arms compare directly; the projected 81-sweep stream seconds are the
   unambiguous numbers.
3. **HPC driver** `benchmark/fgs_sequence_replay_orc.slurm.sh` — the run
   described next.

## IMMEDIATE TASK 1: harvest job 13773451 and issue the go/no-go

Gate-2 ran to completion 2026-09-19 on m12 (zen3, exclusive, qos=test):
job **13773451**, clean finish (`== done ==`, empty stderr). Provenance:
`fgs_acceleration_provenance_20260919.md` (SHAs, input md5s, matrix).
Unharvested outputs on orc:

- `/home/rander39/wt-p021-fgs-gate2/slurm-p021-fgs-replay-gate2-13773451.{out,err}`
- `/home/rander39/wt-p021-fgs-gate2/benchmark/replay_gate2_13773451/`
  (26 files: `topology.txt`, per-arm logs `serial_f64.log`,
  `serial_f32.log`, `handoff_t{4,8,16,32}.log`,
  `rowpar_{f64,f32}_{serial,owner}_t{4,8,16,32}.log`,
  `rowpar_*_interleave_t32.log`, `dag.log`)

Matrix: serial {F64,F32} 1-core; handoff ladder; rowpar precision ×
first-touch × ladder (all socket-0 numactl-pinned); interleave control;
dag priced with the node's measured serial B and t4 handoff h. 12 timed
sweeps/arm, min reported, OPENBLAS_NUM_THREADS=1.

Harvest (scp the results dir + log tails; small text files — inline or
`harvester` subagent), then evaluate the spec's decision rules:

- **Handoff kill switch**: ladder average must be < ~5.8 µs (86,508
  handoffs/solve; 10 µs costs 0.865 s). Local M2 gave 2.52 µs at t4.
- **Projected solve time**: T ≈ 3.0 s remainder + stream + H, stream =
  81 × best sweep min. Reject candidates that cannot reach 6.744 s with
  margin; design for ≤ 5.058 s. (Reference points: 231.891 GB F64 /
  115.945 GB F32 per solve; halving needs 112.7 / 56.3 GB/s; measured
  ceilings 29.4 one core, 164.4 one socket affine.)
- **Source-major vs pull-DAG**: compare rowpar arms against `dag.log`'s
  priced schedule (its last run: forward makespan saturates ~W=8 at
  ~1.74 s/81-sweep lower-stream + absorbed backward — but re-read it with
  the node's actual B/h, and remember the spec: promote split only if a
  complete inclusive schedule comparison at equal precision/placement wins
  with margin; the 2.849 work/span bound caps one-worker-per-pull
  parallelism).
- Sanity-check owner vs serial first-touch vs interleave against the NUMA
  findings (74 vs 164 GB/s; verify with the topology log; sub-page tiles
  can defeat ownership).

**Report the verdict to Ryan before touching the production solver path**
(status file `<topic>_status_20260919.md` in BRAINSTORM/021, naming
convention as before). Include: handoff µs by team size, best useful
bandwidths per arm, projected solve times, the schedule pick, and any
anomalies (e.g. if rowpar scales worse than the synthetic NUMA bench,
diagnose before recommending).

## TASK 2 (after Ryan's go): production implementation

In FastMultipole on the dev branch (NOT the dirty live checkouts — both
FLOWPanel `fastmultipole` and FastMultipole `flowpanel-20260817` carry
unrelated uncommitted work; never commit/revert/build on top of it):

1. Persistent adaptive row-parallel execution of the existing source-major
   nonself cache (new `sweep_order` value or equivalent), consumer-aligned
   placement. Keep lex sequence, far-field refresh, RHS semantics. No
   allocations/task creation inside the leaf loop. The gate-1 shadow
   executor defines the schedule; the gate-2 rowpar team is the reference
   coordination pattern (epoch/done atomics + GC.safepoint in spin loops).
2. Float32 nonself storage / Float64 arithmetic: explicit container change —
   `Matrices{TF}` (containers.jl:1092) couples coefficient and product/RHS
   types; `FastGaussSeidel` couples self/nonself precision. Convert on load,
   accumulate F64, never round strengths to F32. Fallback = selective
   per-block F64 retention, never a relaxed gate.
3. Gate order is binding: correctness (gate-1 harness against the REAL
   implementation, incl. the multi-system question) → numerical gate
   (independent F32 certification, external evaluator rules) → end-to-end
   interleaved A/B vs unchanged champion, ≥1.5× accepted throughput. The
   A/B is a full campaign: tagged worktrees
   (`campaign/<item>-<slug>-YYYYMMDD`), Manifest pins, provenance before
   submission.

## Code anchors (verified at c18e4b46)

FastMultipole `src/solve.jl`: `gs_sweep!` :1214 (lex branch :1277–1302),
`compute_nonself_products!` :923, `scatter_nonself_influence!` :948,
`nonself_influence_matrices` :143 (fill :256–321),
`residual!` :1803 (**shared scratch — threading races without private
scratch**), `color_leaves` :1126 (symmetrized, NOT the pull graph), outer
loop :1380+ (FLOWPanel sets `final_update=false`,
`src/FLOWPanel_solver.jl:1814–1824`); `src/containers.jl:1092` `Matrices{TF}`;
`nearfield_cache.jl:438–478` private-buffer parallel-assembly pattern.

## Known traps (each cost a prior agent time)

- Row-tiled kernels: mathematically equivalent to BLAS, not guaranteed
  bit-identical. Apple BLAS WAS bit-identical in gate 1; re-check on the
  HPC OpenBLAS build, else certify accuracy.
- Source-affine first touch is WRONG for row-worker consumption; verify
  placement with numastat, don't assume.
- Chunked v22 LOST end-to-end (14.264 s, 44 iters) despite faster sweeps;
  its 38.8 active threads measured a Jacobi-lagged schedule.
- Colored and lex have separately calibrated tolerances/orderings — never
  mix their iteration counts.
- The reverse flag repeats forward order — don't silently change it.
- Split design (if promoted): backward products must not overwrite the
  frozen upper accumulator while targets read it; initialize Ux^0 for warm
  starts.
- Sub-range BLAS gemv on strided views: fine; one strided dot per row: never.
- M2 local runs: E-cores drag balanced tiles at t4 (t2 P-cores beat serial)
  — local perf numbers are machinery-validation only.

## House rules (binding)

- Local runs ≤ 4 threads; full-scale timing only on HPC.
- HPC submission Ryan-gated (gate-2 was explicitly approved 2026-09-19;
  new submissions need fresh approval). Monitoring via `hpc-monitor`
  subagent; `ssh orc` needs a live ControlMaster socket (else 2FA — stop
  and ask Ryan to run `ssh orc -fN`).
- Judge runs by outputs, never by sacct exit status.
- Notebook writes and anything outward-facing are Ryan-gated. Ryan-pending
  ledger items now include: origin pushes (merged branches, v21/v22 tags,
  AND `p021-fgs-accel-20260918`), 4 notebook entries, WeakKeyDict/warmstart
  one-line fix.
- Dated status/provenance files in BRAINSTORM/021 per existing convention.
- The orc worktree `/home/rander39/wt-p021-fgs-gate2` is this thread's own;
  don't touch other campaigns' worktrees or the live clones.

## Suggested first moves

1. Read the recommendation + top-5 note end to end, then the 20260918
   status and 20260919 provenance files.
2. Harvest 13773451 (scp `replay_gate2_13773451/` + slurm logs to
   `fgs_r4_followup_evidence_20260914/replay-gate2-13773451/` or similar).
3. Tabulate against the decision rules; write the status file; report
   go/no-go + schedule pick to Ryan. Stop there until he answers.
