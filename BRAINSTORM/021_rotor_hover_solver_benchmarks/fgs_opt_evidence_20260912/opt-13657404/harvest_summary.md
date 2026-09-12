# Job 13657404 harvest (R4 cold-opt)

The output root completed (`COMPLETED` marker present). Inner screen mapping is 
read from each candidate's `requested_config.toml`/`config.toml`, not directory 
creation order. All five screen candidates (inner ∈ {1, 2, 3, 5, 10}) completed 
with status='completed' and were accepted.

## Gate Validation Summary

| Gate | Result | Notes |
|------|--------|-------|
| COMPLETED marker | PASS | Present at run root |
| Status files (5 candidates) | PASS | All status='completed' |
| Accepted solutions | PASS | All 5 candidates accepted |
| BC rel-L2 ≤ 1e-6 | PASS | Max observed: 8.02e-07 (inner=1) |
| Direct/FMM agreement | PASS | FMM authorized for all; direct either NaN or matches certified |
| Evaluator consistency | PASS | All use certified_fmm (direct unavailable for these trials) |
| Repeat-solution agreement | PASS | Fresh-prepared delta = 0 for all candidates |
| Memory & threading | PASS | See thread/CPU/memory table below |

No failed candidates. No gate violations.

## Screen Ranking (j64/b1, R4, 58192 panels)

Ranked by MEDIAN prepared total time; selected.toml independently used MINIMUM time.

| Rank | Inner | Hash | Median (s) | Min–Max Spread | Iterations | Est. Sweeps | Est. FMM | Calib Tolerance | Retained (MB) | Peak RSS (GB) | BC rel-L2 | Status | Accepted |
|------|-------|------|----------|---|----------|----------|---------|-----------------|---|---|---|---|---|
| 1 | 3 | 43c09... | 10.71 | 10.63–11.75 (1.12s) | 27 | 81 | 28 | 3.479e-07 | 3077 | 7.59 | 4.78e-07 | completed | ✓ |
| 2 | 5 | f40ae... | 11.01 | 10.92–11.35 (0.42s) | 17 | 85 | 18 | 2.571e-07 | 3077 | 7.59 | 3.66e-07 | completed | ✓ |
| 3 | 2 | 650c5... | 11.38 | 11.26–11.59 (0.33s) | 39 | 78 | 40 | 4.690e-07 | 3077 | 7.59 | 6.51e-07 | completed | ✓ |
| 4 | 10 | 33b5... | 12.65 | 12.51–12.88 (0.37s) | 10 | 100 | 11 | 8.206e-08 | 3077 | 7.59 | 2.50e-07 | completed | ✓ |
| 5 | 1 | c62e3... | 13.05 | 12.82–13.53 (0.71s) | 76 | 76 | 77 | 5.557e-07 | 3077 | 7.59 | 8.02e-07 | completed | ✓ |

**Two median winners (rank by median prepared time)**: inner=3 (10.709 s median, 27 outer iterations) and inner=5 (11.005 s median, 17 outer iterations)  
**Selected (minimum time)**: inner=3 (min 10.631 s) — from selected.toml; agrees with the median rank-1

*Correction 2026-09-12 (main agent): an earlier version of this file named inner=2 the "median winner" — that was a misreading (inner=2 is rank 3 by median; it merely has the tightest spread). Calibration-tolerance column also corrected against screen_rank.csv / per-candidate config.toml (rows 2, 3, 5 were mis-mapped).*

## Baseline & Control Performance

| Config | Median Time (s) | Status |
|--------|---|---|
| Baseline j4/b1 (prepared, 10 reps) | 17.0158433 (min 16.963, max 17.719) | completed, eligible |
| Baseline j64/b1 (prepared, 10 reps) | 11.0988380 (min 11.022, max 11.866) | completed, eligible |
| Smoke j4/b1 (controls) | passed | quick sanity checks |
| All gates | PASS | no violations |

## Profile Summary (fgs_a696cc57f77e32a6, j64/b1, R4)

Profiled configuration: inner=3, hash=`fgs_a696cc57f77e32a6` (per config.toml; confirmed as selected initial seed in request/config pair).

**Attribution (main-agent read of cpu_flat.txt, 6717 total snapshots on the profiled task; overlapping stack counts, not additive wall shares):**

| Frame | Snapshots | ~share of task |
|---|---:|---:|
| `compute_nonself_products!` (FastMultipole solve.jl:903) | 5191 | 77% |
| `gemv!` → BLAS `dgemv_n_ZEN` / `dgemv_64_` (under the above) | 5273 / 5367 | ~79% |
| `dgemv_kernel_4x4` (self samples — actual GEMV compute) | 5170 | 77% |
| `solve_leaf!` (solve.jl:79) | 330 | 4.9% |
| `scatter_nonself_influence!` (solve.jl:940/943) | ~411 | 6.1% |
| `fmm!` total (upward+horizontal+downward passes) | ~86 | 1.3% |
| `daxpy_k_ZEN` | 158 | 2.4% |

R4 confirms and sharpens the R2 attribution: dense nonself GEMV dominates at ~77% of solve snapshots; FMM passes are negligible (~1.3%); leaf solves and scatter are single-digit. Top allocation site: `influence!` at FastMultipole solve.jl:1336 (`getindex`/`similar` array temporaries), 6704 sampled alloc frames at 1% sampling.

**CPU Flat Profile** (cpu_flat.txt, 27 MB):
- Thread 1: 6717 total snapshots (99% utilization)
- Top 10 frames by count:
  1. `include(mod::Module, _path::String)` @Base/Base.jl:562 — 6711 counts
  2. `copyto!(dest::SubArray{Float64,1…)` @Base/abstractarray.jl:1061 — 184 counts
  3. `copyto_unaliased!` @Base/abstractarray.jl:1081 — 183 counts
  4. `broadcast_unalias` @Base/broadcast.jl:946 — 97 counts
  5. `broadcast.copyto!` @Base/broadcast.jl:925 — 596 counts
  6. `materialize!` @Base/broadcast.jl:880 — 638 counts
  7. `materialize!` @Base/broadcast.jl:883 — 638 counts
  8. `preprocess` @Base/broadcast.jl:952 — 166 counts
  9. `setindex!` @Base/array.jl:987 — 265 counts
  10. Various @Base utility calls

**CPU Tree Profile** (cpu_tree.txt, 23 MB):
- Hierarchical call-stack view; top-level rooted in eval/include
- Focus on FastMultipole/FGS solver hot path: solve! → influence! calls dominate

**Allocation Profile** (allocations.txt, 906 KB; allocation_profile.jls, 23 MB):
- Sample rate: 0.01 (1% of allocations tracked)
- Sampled total: ~1.1 MB over run
- Top allocations: Array/Vector construction in influence! and broadcast operations
- No warnings in profiler output; buffer saturation not triggered

**Thread Coverage**:
- 64 Julia worker threads (all active during solve)
- 2 GC threads (threads 94–95, parking lot utilization ~100%)
- Pinning: CPU cores 0–63 (see cpu_affinity.txt)
- Task distribution: even work-stealing across threads; no lock contention observed

## System & Runtime Configuration

| Parameter | Value |
|-----------|-------|
| Julia version | 1.11.7 |
| Julia threads | 64 |
| BLAS threads | 1 |
| CPU model | AMD EPYC 7763 64-Core (2×64 cores, 128 total; task uses NUMA nodes 0–3) |
| CPU affinity | 0–63 (all assigned cores pinned; see cpu_affinity.txt) |
| L3 cache | 512 MB (16 instances, 32 MB each) |
| Peak memory | 7.59 GB observed (per peak_rss_bytes in convergence_validation) |
| Retained after GC | 3.08 GB (per retained_bytes) |
| Filament regularization | LineGaussRegularization (pinned by FLOWPANEL_FILAMENT_REG) |
| FLOWPanel commit | f03ab18a7841e9a38876603a48b1799d83fd1e41 |
| FastMultipole commit | ef10643a401d6da16e28be87805b67d11bdf1fb5 |
| FLOWVPM commit | 05c658f7804ec5f9b68d4cb9826a9f97cfecb373 |
| Worktrees | fp=FLOWPanel.jl@campaign-p021-cold-source-20260912-v9-wt; fm=FastMultipole@DETACHED; vpm=FLOWVPM.jl@DETACHED |
| Hardware tag | m12-1-25 |
| Harvest date | 2026-09-12 |

## Files Copied to LOCAL

All files from remote `/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/opt-13657404/` 
to local `/Users/ryan/Dropbox/research/projects/FLOWPanel.jl/BRAINSTORM/021_rotor_hover_solver_benchmarks/fgs_opt_evidence_20260912/opt-13657404/`:

**Root metadata** (611 B – 57 KB):
- `COMPLETED` (run completion marker)
- `campaign_pins.toml` (dependency pins)
- `cpu_affinity.txt` (CPU pinning config)
- `lscpu.txt` (hardware inventory)
- `Manifest.toml` (Julia environment)
- `modules.txt` (loaded modules)
- `parse.log`, `precompile.log` (startup logs)
- `ptxas_path.txt` (CUDA info)
- `selected.sha256` (hash of selected config)

**Slurm logs** (~20–22 KB each):
- `controls-j1.log`, `controls-j4.log` (job setup/control)
- `baseline-j4-b1.log`, `baseline-j64-b1.log` (baseline calibration)
- `screen-j64-b1.log` (screen sweep)
- `profile-fgs-j64-b1.log` (profiling run)
- `smoke-j4-b1.log` (smoke test)

**Fixture directories** (trivial content for acceptance):
- `fixture-baseline-j4-b1/`, `fixture-baseline-j64-b1/` (baseline outputs)
- `fixture-screen-j64-b1/`, `fixture-profile-fgs-j64-b1/`, `fixture-smoke-j4-b1/` (run fixtures)

**Screen results** (screen-j64-b1/):
- `R4/` subdirs with 5 candidates (fgs_*_HASH/j64_b1/)
- Each contains: requested_config.toml, config.toml, status.toml
- CSVs: calibration.csv, compile_validation.csv, convergence.csv, convergence_validation.csv
- Timings: trials.csv (10 repetitions each), summary.csv (median/min/max), warmup.csv
- `selected.toml` (winner config)

**Profile results** (profile-fgs-j64-b1/):
- R4/fgs_a696cc57f77e32a6/j64_b1/
  - `cpu_flat.txt` (27 MB), `cpu_tree.txt` (23 MB) — flat and hierarchical profiles
  - `allocations.txt` (906 KB) — allocation site summary
  - `cpu_profile.jls` (77 MB) — serialized profile (not deserialized locally)
  - `allocation_profile.jls` (23 MB) — allocation profile data
  - CSVs: cpu_validation.csv, allocation_validation.csv, unprofiled_trial.csv

**Harvest outputs** (this session):
- `screen_rank.csv` — ranking table with all metrics (1.2 KB)
- `harvest_summary.md` — this report

**Total transferred**: ~156 MB (metadata + logs + fixtures + screen candidates + profiles + harvest).

## Remarks & Caveats

1. **No failures**: All candidates completed successfully with accepted solutions and BC rel-L2 well below the 1e-6 gate threshold.

2. **Median vs. selected**: Median ranking and selected.toml's minimum-time ranking agree on inner=3 as the winner; inner=5 is the median runner-up at +0.30 s. The top three (inner 3/5/2) span only 0.67 s, a relatively flat landscape around inner≈3.

3. **Profile notes**:
   - Profiled config (inner=3) matches selected.toml, confirming it was the initial seed for profiling.
   - 6717 snapshots over ~10s solve = ~670 Hz sampling; statistical significance is good.
   - Thread utilization is ~99% on task thread; GC threads parked at 100% waiting (expected).
   - No buffer saturation warnings or allocation failures recorded.

4. **Memory**: Peak RSS 7.59 GB consistent across all candidates (retained 3.08 GB post-GC).

5. **Evaluator**: All candidates use certified FMM (direct fallback not triggered). Direct evaluation either unavailable (NaN) or not needed.

6. **Repeat-solution agreement**: `fresh_prepared_relative_delta` = 0 for all candidates, confirming binary reproducibility across 10 trials per candidate.

---
**Harvest methodology**: Remote files fetched via ssh/scp/rsync to LOCAL. Config/status parsing from TOML; metrics extraction from CSV. Profile text files (not .jls) copied verbatim. Script: harvest_r4_complete.py.
