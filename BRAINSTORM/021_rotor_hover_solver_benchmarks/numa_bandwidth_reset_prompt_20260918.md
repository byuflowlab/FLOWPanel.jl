# NUMA-placement investigation: context-reset handoff (2026-09-18)

Follows the v22 chunked A/B wrap (commit `3a94f1a`; results in
`fgs_r4_followup_evidence_20260914/chunked-v22-13749231/analysis/ab_summary.md`,
supersedes `fgs_r4_context_reset_20260918b.md` — that task is DONE).

**Ryan's directive (2026-09-18):** Float32 nearfield storage stays staged
and approved, but FIRST investigate bigger speedups via NUMA-aware
placement, taking the FASTEST path. **Submission is PRE-APPROVED by Ryan
for this investigation** (microbenchmark job and, if confirmed, the
in-situ rerun) — no further ask needed to submit those two jobs. Anything
beyond them (code changes to FastMultipole/FLOWPanel, new A/B campaigns)
remains Ryan-gated.

## The hypothesis (why we're here)

v22 verdict: chunked loses everywhere (best j64 14.264 s vs colored@j16
10.116 s; iterations 27→44). BUT the parallel mechanism worked: nearfield
span/iteration shrank 0.346 → 0.277 s at 38.8 avg active threads. The
products themselves sped up at most ~1.3× on 64 threads. Decomposition
(v20 counters, lex, j64, per iteration): nonself products 0.292 s,
scatter 0.032 s, leaf solves 0.020 s → chain 0.346 s/it.

Bandwidth math: each inner sweep streams the full 2.86 GB influence cache
(v20 census: 1,068 leaves), ×3 inner sweeps/iter ≈ 8.6 GB/iter.
- lex serial: 8.6/0.292 ≈ **29 GB/s from ONE core** — near a zen3 core's
  DRAM limit; the kernel is bandwidth-bound (~2 flops/byte dgemv stream).
- chunked 64-thread: ~37 GB/s aggregate — barely above one core, FAR
  below the socket's ~160–200 GB/s.

Node topology (captured in `chunked-v22-13749231/numactl_hardware.txt` +
`j64-trials/numactl_show.txt`): m12 = dual-socket zen3, **NPS4 → 8 NUMA
nodes × 16 CPUs × ~64 GB**; job pinned `cpubind 0–3` (one socket, 4 NUMA
nodes), `membind 0–7`, **policy: default = first-touch**. If the
influence cache is first-touched by one thread, ALL 2.86 GB sits on one
NUMA node (~40–50 GB/s local ceiling) and 64 threads starve on one
memory controller — which fits every measured number. The huge busy-time
(10.7 CPU-s per 0.28 s span) is barrier spin while memory-starved
stragglers finish.

If placement is the cap: interleaving across nodes 0–3 → ~4× delivered
bandwidth → products ~0.08–0.10 s/it → chunked j64 total ≈ 7–8 s even at
44 iterations (beats 10.116), multiplicative with Float32 (→ ~4–5 s).

## The fastest path (staged; execute in order)

**Step 0 — local code check (free, ~minutes).** Read the nearfield
influence-cache assembly in FastMultipole
(`/Users/ryan/Dropbox/research/projects/FastMultipole`, branch
`flowpanel-20260817` @ `c18e4b46`; start from `NearfieldInfluenceCache` /
the cache build used when `cache_leaf_lu=true`, near `src/solve.jl` and
`src/containers.jl`). Determine: is the cache assembled single-threaded
(→ all pages one node; hypothesis maximally likely) or `@threads`
(→ pages scattered by builder — interleave gains less, chunk-affine
first-touch still wins). Record the answer + file:line in the experiment
doc BEFORE submitting. Local runs: ≤4 threads (Ryan's global rule).

**Step 1 — standalone microbenchmark (~3 min compute; THE go/no-go).**
No FLOWPanel/FastMultipole involvement — pure Julia + BLAS dgemv:
- Mimic the cache: ~1,000 Float64 matrices of ~100×550 (match v20 census
  scale: 2.86 GB total; exact shape uncritical, total bytes + block size
  are what matter), plus x/y vectors.
- Arms (each: warm up, then median of ≥10 timed full passes over all
  matrices; report GB/s = bytes/second):
  a. single thread, single-thread first-touch (lex analogue; expect
     ~29 GB/s);
  b. 64 threads static-partitioned, single-thread first-touch (chunked
     v22 analogue; hypothesis predicts ~40–50 GB/s);
  c. same as (b) under `numactl --interleave=0-3` (predicts ~4× (b));
  d. 64 threads, **parallel chunk-affine first-touch** (each thread
     allocates+fills its own partition; predicts ≥ (c)).
- BLAS threads = 1; `Threads.@threads :static`; pin like the campaign
  (the launcher already does `physcpubind 0–63`, `cpubind/nodebind 0–3` —
  copy the v22 sbatch header from
  `benchmark/run_r4_chunked_ab.slurm.sh`); one exclusive m12 node,
  ~30 min wall request. Capture `numastat -p $$` (or
  `/proc/self/numa_maps` summary) inside arms (b) and (d) to SHOW page
  placement, not just infer it.
- Script + sbatch: NEW small files, fine to keep in a scratch dir or
  `benchmark/` (this is a diagnostic, not an official campaign — no
  worktree/pin ceremony needed for a standalone microbenchmark that
  touches no repo code; keep outputs small, CSV + logs only, to
  `~/projects/FLOWPanel.jl/data/p021-cold-20260910/numa-bench-<job>/`).
- **Gate:** (c) or (d) ≥ ~3× over (b) → placement CONFIRMED → step 2.
  <1.5× → placement REFUTED → STOP, report; next levers are Float32
  (already approved) and iteration reduction; per-chunk timers/uncore
  counters are follow-ups for Ryan.

**Step 2 — in-situ confirmation (ONLY if step 1 confirms; ~1 h job).**
Rerun ONLY the j64 activity pair under `numactl --interleave=0-3`,
reusing the pinned v22 setup unchanged: worktrees + env at
`/home/rander39/campaigns/p021-r4-chunked-20260918-v22/`, certified
config `<v22 run>/calibrate/results/chunked_selected.toml` (tolerance
5.348662427942506e-7, chunks=64) — no recalibration, activity stage only
(v22's took ~14 min; drop calibrate/trials stages or write a small
launcher variant that invokes the driver's activity mode directly —
see `AB_MODE` handling in `benchmark/fgs_r4_chunked_ab.jl` /
`test/runtests_r4_chunked_ab_driver.jl`). Include a lex arm as control
(expect ~no change — single reader). Compare nearfield span/it vs
baselines 0.277 (chunked) / 0.346 (lex).
- **Gate:** chunked span/it ≤ ~0.15 s → hypothesis confirmed in the real
  kernel → write up + hand Ryan the decision on step 3 (chunk-affine
  first-touch assembly in FastMultipole = code change, Ryan-gated;
  natural affinity exists via the static chunk→thread map).
- This IS campaign-adjacent: judge by outputs, harvest evidence to
  `fgs_r4_followup_evidence_20260914/numa-insitu-<job>/` SHA256-verified,
  cite the v22 tags (unchanged) in a short provenance note.

**Wrap:** results → a new `numa_placement_findings_<date>.md` in this
directory + 021 item Current status/decision log + log.md + memory
`project_021_solver_benchmarks.md`. OFFER notebook entry (now 4 owed:
v21 A/B, diagnostics ladder, v22 chunked, NUMA). Never write the notebook
without Ryan's approval.

## Cluster gotchas (all bit recent sessions)

- `ssh orc` needs the live ControlMaster socket — auth failure = STOP,
  never trigger MFA.
- Non-login ssh shells have NO slurm on PATH: use
  `/apps/slurm/latest/bin/squeue` (absolute); `bash -lc` injects an ANSI
  banner into stdout.
- Judge jobs by OUTPUTS, never sacct exit status; ≥300 s sacct spacing.
- module `julia/1.11.7-6bmogfl` on login nodes; `JULIA_PKG_PRECOMPILE_AUTO=0`.
- /home was 359.1 G of 400 G cap after the 09-17 storage cycle — this
  work writes only KBs; no storage cycle needed.
- m12 queue has been fast (v22 started 2 min after submit);
  `slurm-availability` probe only if it isn't.

## Known issues / do-not-touch

- FLOWPanel full suite has ONE known NEW failure, unrelated to v22 and to
  this work: `_publish_block_gs_status!` (`src/FLOWPanel_solver.jl:2331`,
  from commit `7fbd68a`) uses the first solver as a WeakKeyDict key; the
  warmstart test's immutable `WarmstartNoopSolver` errors
  (`test/runtests_unit_warmstart.jl:36`). One-line fix (`ismutable` guard
  or `objectid` key) flagged to Ryan — NOT in this task's scope. The
  historical Kutta `:jump` failure did not fire on 2026-09-18.
  FastMultipole full suite: PASS (1,394,073 tests, 2026-09-18).
- Pre-existing dirty files from other campaigns (018 mechanism tests,
  026, data/, examples/) — leave them alone.

## Ryan-pending (do not act)

- Origin pushes: merged branches + v21/v22 campaign tags.
- Notebook entries (4 owed after this — offer at wrap).
- RECENT-VTK archiving approval (208 GiB / 3 runs, listed in
  `fgs_r4_context_reset_20260918b.md`).
- The WeakKeyDict/warmstart one-line fix.
- Step 3 (chunk-affine assembly code change) and Float32 sequencing after
  this investigation reports.
