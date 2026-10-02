# Reset prompt: FGS acceleration — harvest gates 2c+2d, pick precision & schedule, then production (2026-09-19c)

Copy everything below the line into a fresh agent session started in
`~/Dropbox/research/projects/FLOWPanel.jl`.

---

## Where this stands

You are continuing
`BRAINSTORM/021_rotor_hover_solver_benchmarks/fgs_acceleration_recommendation_20260918.md`
(authoritative spec — read it in full; where the companion
`thread_efficiency_top5_20260918.md` disagrees, the recommendation wins).
Also read `CLAUDE.md`, `agent_policies/WORKFLOW.md`, `agent_policies/TESTING.md`,
and `agent_policies/HPC.md` before the corresponding work, plus (in order):
`fgs_acceleration_status_20260919.md` (rev a invalidation),
`fgs_acceleration_status_20260919b.md` (rev b verdict),
`fgs_acceleration_provenance_20260919c.md` (the two jobs you are harvesting).

Goal and acceptance (unchanged): accelerate the prepared cold FGS body solve
at R4. Baseline = colored @ j16 = **10.116 s** (26 iters). Minimum
**≤ 6.744 s** (1.5×), design target **≤ 5.058 s** (2×). Accuracy gate:
authoritative BC relative L2 ≤ 1e-6 via the independent evaluator
(`benchmark/fgs_cold_README.md` in FLOWPanel). Model: T ≈ 3.0 s remainder +
81-sweep stream + H; 231.891 GB F64 / 115.945 GB F32 per solve; useful BW
reported F64-equivalent; minimum needs 61.9 GB/s F64-equiv, halving 112.7.

## Rev b results (job 13773494, harvested — the standing verdict)

True socket-0 binding fixed everything rev a broke. Handoff kill switch PASS
(2.0–4.6 µs across t4–t32). **F64 rowpar fails the minimum everywhere**
(best 54.8 GB/s, owner t4 → T ≈ 7.24 s = 1.40×). **F32 convert-on-load
passes**: best = interleave 0-3 @ t16, 81.6 GB/s F64-equiv → stream 2.843 s
→ T ≈ 5.84 s = 1.73×. Interleave beats owner touch at F32 (81.6 vs 74.8;
sub-page tiles defeat ownership at halved bytes). Priced dag sim's 1.76 s
lower stream implies ~132 GB/s aggregate — beyond demonstrated rates, so the
split was NOT promoted on rev b evidence; gate-2d now measures it for real.
Anomaly on record: real kernels plateau ~55 GB/s actual vs 164.4 synthetic
affine (leaf-serial short bursts, ~2.68 MB mean work/leaf); and rev-b F32
moved only ~41 GB/s ACTUAL bytes vs F64's 55 → the convert kernel was the
limiter — which motivated gate 2c. All rev-b numastat snapshots were empty
(single 60 s shot raced process exit); fixed in the 2c/2d drivers.

## The two jobs to harvest FIRST (Ryan approved 2026-09-19)

Both submitted from `/home/rander39/wt-p021-fgs-gate2` at FastMultipole
**`e904e763`** (branch `p021-fgs-accel-20260918`; chain `c18e4b46` →
`b0946c36` → `a6492f49` → `b1c6c7af` → `e904e763`), CENSUS/EDGES env
overrides md5-verified (`927a…318` / `40f7…020`). qos=test runs them
SEQUENTIALLY (MaxJobsPerUser=1): 13773580 was RUNNING (m12-2-13),
13773581 auto-starts after it.

- **Gate-2c, job 13773580** (`benchmark/fgs_replay_gate2c_f32full_orc.slurm.sh`):
  Ryan's full-F32 question — `--f32-full` = F32 storage AND
  state/accumulate/LU, no convert cost. Serial×3 precisions, handoff t16,
  f32full rowpar owner t4–t32 + interleave t16/t32, f32conv champions rerun
  in-job. Early live numbers: serial F64 27.8 / F32conv 34.0 /
  **F32full 52.2** GB/s — the convert cost is real and large.
- **Gate-2d, job 13773581** (`benchmark/fgs_replay_gate2d_dagteam_orc.slurm.sh`):
  REAL split dual-layout executor (`--mode dagteam`): target-major lower
  pulls driven by readiness counters over the 48,167 directed lower edges,
  source-major backward upper products as filler tasks, target-owned
  boundary reduction (frozen u^s never overwritten mid-sweep). F64 ladder
  t4–t32 + socket control, F32conv/F32full at t16/t32, in-job rowpar
  references at matching precision/placement — the spec's inclusive
  equal-precision comparison. Local M2 smoke passed with exact structural
  invariants (1,068 tasks, 3 roots, 48,167 edges, lower/upper bytes
  1,513,294,744 / 1,349,555,288).

Outputs: `~/wt-p021-fgs-gate2/slurm-p021-fgs-replay-gate2{c,d}-*.{out,err}`
and `~/wt-p021-fgs-gate2/benchmark/replay_gate2c_13773580/`,
`.../replay_gate2d_13773581/`. Harvest to
`BRAINSTORM/021_rotor_hover_solver_benchmarks/fgs_r4_followup_evidence_20260914/`
as `replay-gate2c-13773580/` and `replay-gate2d-13773581/` (rev a/b harvests
sit beside them). Judge by outputs, never sacct status.

## IMMEDIATE TASK: harvest both, decide, report, STOP

1. Harvest both jobs (scp slurm logs + results dirs as above). Check the
   numastat snapshots actually captured this time (poll-and-keep-last fix).
2. **Precision pick (2c)**: does f32full beat f32conv at equal placement on
   the socket? Tabulate the ladder; if the best f32full stream approaches
   2.06 s (112.7 GB/s F64-equiv), the 2× design target is in measured
   reach. Adoption is CONDITIONAL on the independent evaluator at 1e-6
   (F32 eps is 1.2e-7 — only one decade of headroom; GS accumulates over
   81 sweeps and the internal residual can flatter the rounded operator).
   Fallback ladder: f32full → f32conv (certified rev-b GO at T≈5.84 s) →
   selective per-block F64. Never relax the gate.
3. **Schedule pick (2d)**: dagteam vs rowpar best-vs-best AND at matched
   team size/precision/placement. Sanity-check dagteam against the 2.849
   work/span bound. Promote the split ONLY on a win with margin enough to
   justify the rebuild (spec rule); otherwise source-major stands.
4. Write `fgs_acceleration_status_20260920.md` (or next free date-name) in
   BRAINSTORM/021: tables, precision pick, schedule pick, updated projected
   solve time, go/no-go for production. **Report to Ryan and STOP before
   touching the production solver path.**

If a job died: diagnose from its logs; resubmission needs fresh Ryan
approval (sbatch from the worktree with the CENSUS/EDGES overrides —
omitting them kills the job in 1 s; that was rev-b's dead first attempt
13773493).

## After Ryan's go: production implementation (TASK 2 of the 20260919b prompt)

In FastMultipole dev branch `p021-fgs-accel-20260918` (NOT the dirty live
checkouts): persistent adaptive row-parallel execution of the source-major
nonself cache (or split, if 2d overturns), winning precision mode,
consumer-aligned/interleave placement per 2c. Keep lex sequence, far-field
refresh, RHS semantics (`+= old` then `-= new`, never a delta). No
allocations/task creation in the leaf loop. Gate order: correctness
(gate-1 harness `test/fgs_rowpar_gate1_test.jl` vs the REAL implementation,
incl. multi-system) → numerical (independent evaluator) → end-to-end
interleaved A/B vs unchanged champion, ≥1.5× accepted throughput, full
campaign ceremony (tagged worktrees `campaign/<item>-<slug>-YYYYMMDD`,
Manifest pins, provenance before submission).

## Code anchors (verified at c18e4b46; replay additions at e904e763)

FastMultipole `src/solve.jl`: `gs_sweep!` :1214 (lex branch :1277–1302),
`compute_nonself_products!` :923, `scatter_nonself_influence!` :948,
`nonself_influence_matrices` :143 (fill :256–321), `residual!` :1803
(**shared scratch — threading races without private scratch**),
`color_leaves` :1126 (symmetrized, NOT the pull graph), outer loop :1380+
(FLOWPanel sets `final_update=false`, `src/FLOWPanel_solver.jl:1814–1824`);
`src/containers.jl:1092` `Matrices{TF}` couples coefficient and product/RHS
types; `nearfield_cache.jl:438–478` private-buffer parallel-assembly
pattern. Replay harness: `benchmark/fgs_sequence_replay.jl` (modes
serial/rowpar/handoff/dag/dagteam; `--f32` convert, `--f32-full`).

## Known traps (each cost a prior agent time)

- "node 0" ≠ "socket 0" on NPS4 EPYC — socket 0 = nodes 0-3; verify
  placement, never assume.
- sbatch of any replay driver NEEDS the CENSUS/EDGES env overrides.
- qos=test: MaxJobsPerUser=1 — approved jobs queue sequentially.
- Row-tiled kernels: not guaranteed bit-identical to BLAS; Apple BLAS was
  bit-identical in gate 1, re-check on cluster OpenBLAS else certify
  accuracy.
- Source-affine first touch is WRONG for row-worker consumption; sub-page
  tiles can defeat owner touch (rev b measured it at F32).
- Chunked v22 LOST end-to-end despite faster sweeps (Jacobi-lagged
  schedule); colored and lex have separately calibrated tolerances — never
  mix iteration counts.
- The reverse flag repeats forward order — don't silently change it.
- Split design: backward products must not overwrite the frozen upper
  accumulator mid-sweep; initialize Ux^0 for warm starts.
- Two-system `FastGaussSeidel((sysA,sysB),(sysA,sysB))` fails to CONSTRUCT
  at `c18e4b46` (BoundsError, pre-existing) — production must fix or
  declare before multi-rotor use.
- M2 local runs: machinery validation only (≤4 threads house rule).
- Judge runs by outputs, never sacct exit status.
- zsh does not word-split unquoted `$var` — spell out CLI flag lists.

## House rules (binding)

- Local runs ≤ 4 threads; full-scale timing only on HPC.
- HPC submission Ryan-gated (13773580/13773581 were approved 2026-09-19;
  new submissions need fresh approval). Monitoring via `hpc-monitor`
  subagent; `ssh orc` needs a live ControlMaster socket (else 2FA — stop
  and ask Ryan to run `ssh orc -fN`). Slurm CLI on orc: `source
  /etc/profile` (the `module load slurm` step errors, sbatch/squeue work
  regardless).
- Notebook writes and anything outward-facing are Ryan-gated. Ryan-pending
  ledger: origin pushes (merged branches, v21/v22 tags, AND
  `p021-fgs-accel-20260918`), 4 notebook entries, WeakKeyDict/warmstart fix.
- Dated status/provenance files in BRAINSTORM/021 per existing convention.
- The orc worktree `/home/rander39/wt-p021-fgs-gate2` is this thread's own;
  don't touch other campaigns' worktrees or the live clones. Local worktree:
  `/private/tmp/fastmultipole-p021-fgs-accel-20260918`. Sync dance: push a
  fresh ref to the orc clone via ssh (`git push orc HEAD:refs/heads/<ref>`),
  then in the worktree `git fetch /home/rander39/projects/FastMultipole
  <ref> && git merge --ff-only FETCH_HEAD` (a checked-out branch can't be
  pushed to directly; https origin push needs Ryan's token).

## Suggested first moves

1. Read the recommendation end to end, then the three status/provenance
   files listed at the top.
2. Check both jobs finished (hpc-monitor); harvest both evidence dirs.
3. Tabulate 2c (precision) and 2d (schedule) against the decision rules;
   write the status file; report picks + go/no-go to Ryan; stop.
