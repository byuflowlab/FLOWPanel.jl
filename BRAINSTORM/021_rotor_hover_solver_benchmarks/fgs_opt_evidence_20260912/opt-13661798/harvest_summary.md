# Job 13661798 harvest (R4 P/MAC screen, stage-B, base inner=5)

## Job Facts

| Field | Value |
|---|---|
| Job ID | 13661798 |
| Node | m12-1-25 |
| Elapsed | 2h37m38s |
| Stage | screen_profile (baseline-j4-b1, baseline-j64-b1, screen-j64-b1, profile-fgs-j64-b1) |
| Panel count | 58,192 |
| Rung | R4 |
| Base config | inner=5, leaf=100, P=8, MAC=0.4 (stage-A median runner-up, opt-13660643) |
| Screen set | `P:6,10;MAC:0.3,0.5` (`screen-j64-b1/provenance.toml` `screen_set`) |
| COMPLETED marker | `completed` (verbatim contents of `COMPLETED`) |

## Baselines

Seed configuration (both baselines, from each `config.toml`): inner=3, leaf=100, P=8, MAC=0.4, tolerance=3.479128881193055e-7 (note: the baselines in this job use the *original* seed base, inner=3, not the stage-B inner=5 base being screened here).

| Baseline | Median (s) | Min (s) | Max (s) | Spread (s) | Repetitions | Eligible |
|---|---|---|---|---|---|---|
| baseline-j4-b1 (R4/fgs_a696cc57f77e32a6/j4_b1) | 16.8442106 | 16.7297441 | 17.5003216 | 0.77057752 | 10 | true |
| baseline-j64-b1 (R4/fgs_a696cc57f77e32a6/j64_b1) | 11.0814568 | 11.0223721 | 11.8720526 | 0.849680472 | 10 | true |

Both baseline `status.toml` = `completed`.

## Screen Ranking (j64/b1, R4, 58192 panels)

5-candidate screen: base inner=5/leaf=100/P=8/MAC=0.4 plus one-factor neighbors P∈{6,10}, MAC∈{0.3,0.5}. Ranked by MEDIAN prepared seconds ascending; 2 candidates (P=6, MAC=0.5) failed calibration.

| Rank | Dir | Inner | Leaf | P | MAC | Tolerance | Median (s) | Min (s) | Max (s) | Spread (s) | Reps | Eligible | Outer Iter | Authoritative rel-L2 | FMM certified | Rel. solution delta | Status |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | fgs_f40aebe0fdcb61ea | 5 | 100 | 8 | 0.4 | 2.5712316195808637e-07 | 11.125801 | 11.0531718 | 11.609372 | 0.556200234 | 10 | true | 18 | 3.65840408e-07 | true | 0 | completed |
| 2 | fgs_360414f33341980a | 5 | 100 | 10 | 0.4 | 2.5712327255007005e-07 | 11.7918432 | 11.6013373 | 12.0987384 | 0.497401122 | 10 | true | 18 | 3.30384916e-07 | true | 0 | completed |
| 3 | fgs_28b71b3b4615f119 | 5 | 100 | 8 | 0.3 | 2.2715585509969438e-07 | 15.3764733 | 15.3612248 | 15.8215412 | 0.460316416 | 10 | true | 18 | 2.56096621e-07 | true | 0 | completed |
| — | fgs_799d55a10c28378b | 5 | 100 | 6 | 0.4 | 0.0 (requested; not calibrated) | — | — | — | — | — | — | — | — | — | — | failed |
| — | fgs_64a1ce9afb77f9ce | 5 | 100 | 8 | 0.5 | 0.0 (requested; not calibrated) | — | — | — | — | — | — | — | — | — | — | failed |

Failed candidates' first error line (`status.toml`, `error` field, first line — both candidates share the identical message):
`FGS staircase has no certified crossing with a decreasing successor; capped candidate`

`screen-j64-b1/selected.toml` (harness-selected candidate) = inner=5, leaf=100, P=8, MAC=0.4, tolerance=2.5712316195808637e-7 — matches rank 1 above.

## Gates

| Gate | Result |
|---|---|
| Worst `authoritative_rel_l2` across the 3 completed candidates + 2 baselines | 3.65840408e-07 (candidate fgs_f40aebe0fdcb61ea, P=8/MAC=0.4, inner=5) — note this is the *largest* value observed in this job's set, since the two baselines (inner=3 seed) both report 4.77952513e-07, which is larger still; worst overall across all 5 dirs with a `convergence_validation.csv` is 4.77952513e-07 (both baselines) |
| `fmm_certified` values observed | `true` for all 3 completed candidates and both baselines (no `false` observed; not applicable to the 2 failed candidates, which have no `convergence_validation.csv`) |
| `relative_solution_delta` values observed | `0` for all 3 completed candidates and both baselines |
| `status.toml` not "completed" | 2 (`fgs_799d55a10c28378b` P=6, `fgs_64a1ce9afb77f9ce` MAC=0.5) — both `status = "failed"` |
| `eligible = false` | none found among the 3 completed candidates' + 2 baselines' `summary.csv` prepared rows (all `true`) |

## Tolerance Map

Calibrated tolerance per completed candidate, read from each `config.toml` (`tolerance` field). Failed candidates have no calibrated tolerance (`requested_config.toml` carries the placeholder `tolerance = 0.0`, never resolved).

| Inner | Leaf | P | MAC | Tolerance | Status |
|---|---|---|---|---|---|
| 5 | 100 | 8 | 0.3 | 2.2715585509969438e-07 | completed |
| 5 | 100 | 8 | 0.4 | 2.5712316195808637e-07 | completed |
| 5 | 100 | 8 | 0.5 | — (calibration failed) | failed |
| 5 | 100 | 6 | 0.4 | — (calibration failed) | failed |
| 5 | 100 | 10 | 0.4 | 2.5712327255007005e-07 | completed |

## Profile

Profiled configuration (`profile-fgs-j64-b1/R4/fgs_a696cc57f77e32a6/j64_b1/config.toml`): inner=3, leaf=100, P=8, MAC=0.4, tolerance=3.479128881193055e-7, cache_leaf_lu=true, max_iterations=300, rlx=1.0, sweep_order="lexicographic", rung="R4" — note this profiled config is the seed base (inner=3), matching the baselines in this job, not the inner=5 base being screened.

Total sample count (`cpu_flat.txt` header for the profiled Task): **Thread 1, Task 0x00007f6848824010, "Total snapshots: 6727. Utilization: 98%"**.

Top 5 flat frames by the `Count` column within Thread 1's block (`cpu_flat.txt`, lines up to the start of Thread 2 at line 4586), ties broken by file order (stable sort):

| Count | Overhead | File | Line | Function |
|---|---|---|---|---|
| 6725 | 0 | — | — | `[any unknown stackframes]` |
| 6722 | 0 | …-1-dot-11/src/julia.h | 2157 | `jl_apply` |
| 6721 | 0 | @Base/Base.jl | 562 | `include(mod::Module, _path::String…)` |
| 6721 | 0 | @Base/boot.jl | 430 | `eval` |
| 6721 | 0 | @Base/client.jl | 531 | `_start()` |

(Caveat: the `Count` column is inclusive/whole-stack count, so the top rows (after the `[any unknown stackframes]` placeholder) are tied at 6721/6722 — the same linear call chain from program entry to `solve!`, not differentiating hot compute.)

Top self-time frame (`Overhead` column, i.e. exclusive time), `cpu_flat.txt`, Thread 1 block: `dgemv_kernel_4x4` @ `…lia/libopenblas64_.so` — Count 5106, Overhead 5106 (largest self-time entry by a wide margin; next is `setindex!`@array.jl:987 at Count 298/Overhead 273).

Top 3 allocation sites from `allocations.txt` (sample_rate=0.01, sampled_bytes=1098952), by leading count column, innermost application frame (all three share the identical call site):

| Sampled count | Innermost frames (from stack, in order) |
|---|---|
| 11600 | `similar`@abstractarray.jl:822 [inlined] → `getindex`@array.jl:938 [inlined] → `influence!(...)`@FastMultipole solve.jl:1336 |
| 6072 | (same site) `influence!(...)`@FastMultipole solve.jl:1336 |
| 6072 | (same site) `influence!(...)`@FastMultipole solve.jl:1336 |

`unprofiled_trial.csv` row (single trial, profiling excluded):

```
setup_seconds,solve_seconds,total_seconds,allocated_bytes,gc_seconds,retained_bytes,process_peak_rss_bytes,iterations,estimated_inner_sweeps,estimated_fmm_passes,work_count_kind,formulation_subsolves,solved,fmm_rel_l2,fmm_rel_max,fmm_certified,fmm_seconds,epsilon_requested,direct_rel_l2,evaluator_delta,direct_seconds,authoritative_evaluator,authoritative_rel_l2,accepted,eligible
0,16.835992,16.835992,642646688,0.094418815,3221617552,4540772352,27,81,28,estimates; -1 means unavailable,0,true,4.77952513e-07,2.32095749e-05,true,5.17447229,1.22339544e-09,NaN,NaN,0,certified_fmm,4.77952513e-07,true,true
```

## Provenance

`selected.sha256` (root of opt-13661798):
```
092cd29526572e89658406597b378a91ef08d242bebf971e978f501b1342f2f4  /home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/opt-13657404/smoke-j4-b1/selected.toml
```
(Note: this hash/path points at an earlier job's — 13657404 — `smoke-j4-b1/selected.toml`, not a file inside this job's own tree; reported verbatim, not reconciled.)

`screen-j64-b1/screen_bases.toml` (path: `.../opt-13661798/screen-j64-b1/screen_bases.toml`) — one saved base, comment header:

> "BRAINSTORM 021 R4 tuning stage B (handoff step 7): P/MAC neighbors around retained stage-A base inner=5/leaf=100 (opt-13660643 median runner-up, 11.2826077 s prepared j64/b1). Tolerance is its stage-A calibrated value (harness resets to 0 and recalibrates every roster point)."

| Base | inner | leaf | P | MAC | tolerance |
|---|---|---|---|---|---|
| 1 | 5 | 100 | 8 | 0.4 | 2.5712316195808637e-07 |

Package pins (`screen-j64-b1/provenance.toml` and root `campaign_pins.toml`, consistent):
- FLOWPanel: sha `f03ab18a7841e9a38876603a48b1799d83fd1e41`, tag `campaign/p021-cold-exec-20260912-v9`
- FastMultipole: sha `ef10643a401d6da16e28be87805b67d11bdf1fb5`, tag `campaign/p021-cold-exec-20260910-v1`
- FLOWVPM: sha `05c658f7804ec5f9b68d4cb9826a9f97cfecb373`, tag `campaign/p021-cold-exec-20260910-v1`

Other `provenance.toml` fields (screen-j64-b1): hostname=`m12-1-25`, julia_version=`1.11.7`, julia_threads=64, requested_blas_threads=1, blas=`LBTConfig([ILP64] libopenblas64_.so)`, filament_reg=`LineGaussRegularization`, timing_scope=`frozen _solve!; reset and BC diagnostics excluded; no formulation subsolves`, minimum_reps=10, screen_base_sha256=`1d0f6290ecae06d6679d364230ace7fd2d9e8f65cc2cd1af5c9c330f68cf4ba3`, manifest_sha256=`03cdb054fc315cb9a48cbb0ddf0331235a1b41484b4ac7889f4787b2cb64ce67`.
