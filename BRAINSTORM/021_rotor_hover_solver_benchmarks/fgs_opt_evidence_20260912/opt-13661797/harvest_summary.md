# Job 13661797 harvest (R4 P/MAC screen, stage-B, base inner=3)

## Job Facts

| Field | Value |
|---|---|
| Job ID | 13661797 |
| Node | m12-1-17 |
| Elapsed | 2h40m45s |
| Stage | screen_profile (baseline-j4-b1, baseline-j64-b1, screen-j64-b1, profile-fgs-j64-b1) |
| Panel count | 58,192 |
| Rung | R4 |
| Base config | inner=3, leaf=100, P=8, MAC=0.4 (stage-A median winner, opt-13660643) |
| Screen set | `P:6,10;MAC:0.3,0.5` (`screen-j64-b1/provenance.toml` `screen_set`) |
| COMPLETED marker | `completed` (verbatim contents of `COMPLETED`) |

## Baselines

Seed configuration (both baselines, from each `config.toml`): inner=3, leaf=100, P=8, MAC=0.4, tolerance=3.479128881193055e-7.

| Baseline | Median (s) | Min (s) | Max (s) | Spread (s) | Repetitions | Eligible |
|---|---|---|---|---|---|---|
| baseline-j4-b1 (R4/fgs_a696cc57f77e32a6/j4_b1) | 17.2611866 | 17.2117292 | 17.8971683 | 0.68543908 | 10 | true |
| baseline-j64-b1 (R4/fgs_a696cc57f77e32a6/j64_b1) | 11.2115955 | 11.1075987 | 11.5877532 | 0.48015455 | 10 | true |

Both baseline `status.toml` = `completed`.

## Screen Ranking (j64/b1, R4, 58192 panels)

5-candidate screen: base inner=3/leaf=100/P=8/MAC=0.4 plus one-factor neighbors P∈{6,10}, MAC∈{0.3,0.5}. Ranked by MEDIAN prepared seconds ascending; 2 candidates (P=6, MAC=0.5) failed calibration.

| Rank | Dir | Inner | Leaf | P | MAC | Tolerance | Median (s) | Min (s) | Max (s) | Spread (s) | Reps | Eligible | Outer Iter | Authoritative rel-L2 | FMM certified | Rel. solution delta | Status |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | fgs_43c0957fd1108d49 | 3 | 100 | 8 | 0.4 | 3.479128881193055e-07 | 11.0144373 | 10.9351946 | 11.243838 | 0.308643431 | 10 | true | 28 | 4.77952513e-07 | true | 0 | completed |
| 2 | fgs_1986c75bfc7c003c | 3 | 100 | 10 | 0.4 | 3.479130058413085e-07 | 11.9685691 | 11.7167291 | 12.1239953 | 0.407266215 | 10 | true | 28 | 4.49123195e-07 | true | 0 | completed |
| 3 | fgs_06bf28cf18db65be | 3 | 100 | 8 | 0.3 | 4.924606105994385e-07 | 14.9706686 | 14.4195788 | 16.9086529 | 2.4890741 | 10 | true | 27 | 6.09081208e-07 | true | 0 | completed |
| — | fgs_b139402530c377cb | 3 | 100 | 6 | 0.4 | 0.0 (requested; not calibrated) | — | — | — | — | — | — | — | — | — | — | failed |
| — | fgs_26c6c360f62b7979 | 3 | 100 | 8 | 0.5 | 0.0 (requested; not calibrated) | — | — | — | — | — | — | — | — | — | — | failed |

Failed candidates' first error line (`status.toml`, `error` field, first line — both candidates share the identical message):
`FGS staircase has no certified crossing with a decreasing successor; capped candidate`

`screen-j64-b1/selected.toml` (harness-selected candidate) = inner=3, leaf=100, P=8, MAC=0.4, tolerance=3.479128881193055e-7 — matches rank 1 above.

## Gates

| Gate | Result |
|---|---|
| Worst `authoritative_rel_l2` across the 3 completed candidates + 2 baselines | 6.09081208e-07 (candidate fgs_06bf28cf18db65be, P=8/MAC=0.3) |
| `fmm_certified` values observed | `true` for all 3 completed candidates and both baselines (no `false` observed; not applicable to the 2 failed candidates, which have no `convergence_validation.csv`) |
| `relative_solution_delta` values observed | `0` for all 3 completed candidates and both baselines |
| `status.toml` not "completed" | 2 (`fgs_b139402530c377cb` P=6, `fgs_26c6c360f62b7979` MAC=0.5) — both `status = "failed"` |
| `eligible = false` | none found among the 3 completed candidates' + 2 baselines' `summary.csv` prepared rows (all `true`) |

## Tolerance Map

Calibrated tolerance per completed candidate, read from each `config.toml` (`tolerance` field). Failed candidates have no calibrated tolerance (`requested_config.toml` carries the placeholder `tolerance = 0.0`, never resolved).

| Inner | Leaf | P | MAC | Tolerance | Status |
|---|---|---|---|---|---|
| 3 | 100 | 8 | 0.3 | 4.924606105994385e-07 | completed |
| 3 | 100 | 8 | 0.4 | 3.479128881193055e-07 | completed |
| 3 | 100 | 8 | 0.5 | — (calibration failed) | failed |
| 3 | 100 | 6 | 0.4 | — (calibration failed) | failed |
| 3 | 100 | 10 | 0.4 | 3.479130058413085e-07 | completed |

## Profile

Profiled configuration (`profile-fgs-j64-b1/R4/fgs_a696cc57f77e32a6/j64_b1/config.toml`): inner=3, leaf=100, P=8, MAC=0.4, tolerance=3.479128881193055e-7, cache_leaf_lu=true, max_iterations=300, rlx=1.0, sweep_order="lexicographic", rung="R4" — same hash/config as the seed baselines and screen rank-1 candidate.

Total sample count (`cpu_flat.txt` header for the profiled Task): **Thread 1, Task 0x00007f51c7e5c010, "Total snapshots: 6740. Utilization: 98%"**.

Top 5 flat frames by the `Count` column within Thread 1's block (`cpu_flat.txt`, lines up to the start of Thread 2 at line 4443), ties broken by file order (stable sort):

| Count | Overhead | File | Line | Function |
|---|---|---|---|---|
| 6734 | 0 | — | — | `[any unknown stackframes]` |
| 6727 | 0 | @Base/Base.jl | 562 | `include(mod::Module, _path::String…)` |
| 6727 | 0 | @Base/boot.jl | 430 | `eval` |
| 6727 | 0 | @Base/client.jl | 531 | `_start()` |
| 6727 | 0 | @Base/client.jl | 323 | `exec_options(opts::Base.JLOptions)` |

(Caveat: the `Count` column is inclusive/whole-stack count, so the top rows (after the `[any unknown stackframes]` placeholder) are tied at 6727 — the same linear call chain from program entry to `solve!`, not differentiating hot compute. A 6th row, `@Base/essentials.jl:1057 #invokelatest#2`, is also tied at 6727 but excluded from the top-5 by file order.)

Top self-time frame (`Overhead` column, i.e. exclusive time), `cpu_flat.txt`, Thread 1 block: `dgemv_kernel_4x4` @ `…lia/libopenblas64_.so` — Count 5158, Overhead 5158 (largest self-time entry by a wide margin; next is `setindex!`@array.jl:987 at Count 306/Overhead 291).

Top 3 allocation sites from `allocations.txt` (sample_rate=0.01, sampled_bytes=1132640), by leading count column, innermost application frame (all three share the identical call site):

| Sampled count | Innermost frames (from stack, in order) |
|---|---|
| 6072 | `similar`@abstractarray.jl:822 [inlined] → `getindex`@array.jl:938 [inlined] → `influence!(...)`@FastMultipole solve.jl:1336 |
| 3280 | (same site) `influence!(...)`@FastMultipole solve.jl:1336 |
| 2680 | (same site) `influence!(...)`@FastMultipole solve.jl:1336 |

`unprofiled_trial.csv` row (single trial, profiling excluded):

```
setup_seconds,solve_seconds,total_seconds,allocated_bytes,gc_seconds,retained_bytes,process_peak_rss_bytes,iterations,estimated_inner_sweeps,estimated_fmm_passes,work_count_kind,formulation_subsolves,solved,fmm_rel_l2,fmm_rel_max,fmm_certified,fmm_seconds,epsilon_requested,direct_rel_l2,evaluator_delta,direct_seconds,authoritative_evaluator,authoritative_rel_l2,accepted,eligible
0,17.38661,17.38661,642601024,0.097015564,3221617552,4542525440,27,81,28,estimates; -1 means unavailable,0,true,4.77952513e-07,2.32095749e-05,true,5.1273107,1.22339544e-09,NaN,NaN,0,certified_fmm,4.77952513e-07,true,true
```

## Provenance

`selected.sha256` (root of opt-13661797):
```
092cd29526572e89658406597b378a91ef08d242bebf971e978f501b1342f2f4  /home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/opt-13657404/smoke-j4-b1/selected.toml
```
(Note: this hash/path points at an earlier job's — 13657404 — `smoke-j4-b1/selected.toml`, not a file inside this job's own tree; reported verbatim, not reconciled.)

`screen-j64-b1/screen_bases.toml` (path: `.../opt-13661797/screen-j64-b1/screen_bases.toml`) — one saved base, comment header:

> "BRAINSTORM 021 R4 tuning stage B (handoff step 7): P/MAC neighbors around retained stage-A base inner=3/leaf=100 (opt-13660643 median winner, 11.1603485 s prepared j64/b1). Tolerance is its stage-A calibrated value (harness resets to 0 and recalibrates every roster point)."

| Base | inner | leaf | P | MAC | tolerance |
|---|---|---|---|---|---|
| 1 | 3 | 100 | 8 | 0.4 | 3.479128881193055e-07 |

Package pins (`screen-j64-b1/provenance.toml` and root `campaign_pins.toml`, consistent):
- FLOWPanel: sha `f03ab18a7841e9a38876603a48b1799d83fd1e41`, tag `campaign/p021-cold-exec-20260912-v9`
- FastMultipole: sha `ef10643a401d6da16e28be87805b67d11bdf1fb5`, tag `campaign/p021-cold-exec-20260910-v1`
- FLOWVPM: sha `05c658f7804ec5f9b68d4cb9826a9f97cfecb373`, tag `campaign/p021-cold-exec-20260910-v1`

Other `provenance.toml` fields (screen-j64-b1): hostname=`m12-1-17`, julia_version=`1.11.7`, julia_threads=64, requested_blas_threads=1, blas=`LBTConfig([ILP64] libopenblas64_.so)`, filament_reg=`LineGaussRegularization`, timing_scope=`frozen _solve!; reset and BC diagnostics excluded; no formulation subsolves`, minimum_reps=10, screen_base_sha256=`e09ecb6f8437060d2034526fcf033a89e6e884c3c882af47f533adf89a6d83b5`, manifest_sha256=`03cdb054fc315cb9a48cbb0ddf0331235a1b41484b4ac7889f4787b2cb64ce67`.
