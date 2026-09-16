# Job 13660643 harvest (R4 leaf screen, stage-A)

## Job Facts

| Field | Value |
|---|---|
| Job ID | 13660643 |
| Node | m12-1-25 |
| Elapsed | 3h52m52s |
| Stage | screen_profile (baseline-j4-b1, baseline-j64-b1, screen-j64-b1, profile-fgs-j64-b1) |
| Panel count | 58,192 |
| Rung | R4 |
| Screen set | `leaf:25,50,200` (per screen-j64-b1 provenance.toml `screen_set`) |
| COMPLETED marker | `completed` (verbatim contents of `COMPLETED`) |

## Baselines

Seed configuration (both baselines): inner=3, leaf=100, P=8, MAC=0.4, tolerance=3.479128881193055e-7 (from each `config.toml`; identical to screen candidate rank 1).

| Baseline | Median (s) | Min (s) | Max (s) | Spread (s) | Repetitions | Eligible |
|---|---|---|---|---|---|---|
| baseline-j4-b1 (R4/fgs_a696cc57f77e32a6/j4_b1) | 17.129249 | 16.9065676 | 17.6554496 | 0.748882056 | 10 | true |
| baseline-j64-b1 (R4/fgs_a696cc57f77e32a6/j64_b1) | 10.9684716 | 10.9033167 | 11.3608665 | 0.457549848 | 10 | true |

Both baseline `status.toml` = `completed`.

## Screen Ranking (j64/b1, R4, 58192 panels)

8-candidate saved-base screen: bases inner=3 and inner=5 at leaf=100, plus leaf ∈ {25,50,200} one-factor neighbors of each. Ranked by MEDIAN prepared seconds ascending.

| Rank | Dir | Inner | Leaf | P | MAC | Tolerance | Median (s) | Min (s) | Max (s) | Spread (s) | Reps | Eligible | Outer Iter | Authoritative rel-L2 | FMM certified | Rel. solution delta |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | fgs_43c0957fd1108d49 | 3 | 100 | 8 | 0.4 | 3.479128881193055e-07 | 11.1603485 | 11.0767292 | 11.7159075 | 0.63917833 | 10 | true | 28 | 4.77952513e-07 | true | 0 |
| 2 | fgs_f40aebe0fdcb61ea | 5 | 100 | 8 | 0.4 | 2.5712316195808637e-07 | 11.2826077 | 11.1205288 | 11.3245494 | 0.204020612 | 10 | true | 18 | 3.65840408e-07 | true | 0 |
| 3 | fgs_c86061a8bfbba8db | 3 | 50 | 8 | 0.4 | 2.3864048149973856e-07 | 11.5759369 | 11.4392816 | 11.7244017 | 0.285120095 | 10 | true | 30 | 4.91488859e-07 | true | 0 |
| 4 | fgs_e69901af1bbc5516 | 3 | 25 | 8 | 0.4 | 3.80292896308757e-07 | 11.7758337 | 11.5304376 | 13.3462421 | 1.81580454 | 10 | true | 29 | 6.68094917e-07 | true | 0 |
| 5 | fgs_8a04b989bea2b26c | 5 | 50 | 8 | 0.4 | 1.2891762136325558e-07 | 11.9747802 | 11.8696902 | 12.0783038 | 0.208613622 | 10 | true | 20 | 3.7764383e-07 | true | 0 |
| 6 | fgs_03b8211c1747f926 | 5 | 25 | 8 | 0.4 | 1.5209757150978087e-07 | 12.3897912 | 12.3238828 | 12.8369401 | 0.513057279 | 10 | true | 20 | 4.49557946e-07 | true | 0 |
| 7 | fgs_46208a8714edb343 | 5 | 200 | 8 | 0.4 | 3.9670158724694127e-07 | 12.7605998 | 12.6818926 | 12.8982438 | 0.216351215 | 10 | true | 17 | 4.41290605e-07 | true | 0 |
| 8 | fgs_5708e3b50dbdfd42 | 3 | 200 | 8 | 0.4 | 4.252870058799867e-07 | 12.8517574 | 12.6593747 | 12.8891094 | 0.229734718 | 10 | true | 27 | 4.99975993e-07 | true | 0 |

All 8 candidate `status.toml` = `completed`; all 8 `summary.csv` prepared rows have `eligible = true`.

`screen-j64-b1/selected.toml` (harness-selected candidate) = inner=3, leaf=100, P=8, MAC=0.4, tolerance=3.479128881193055e-7 — matches rank 1 above.

## Gates

| Gate | Result |
|---|---|
| Worst `authoritative_rel_l2` across the 8 candidates + 2 baselines | 6.68094917e-07 (candidate fgs_e69901af1bbc5516, inner=3/leaf=25) |
| `fmm_certified` values observed | `true` for all 8 candidates and both baselines (no `false` observed) |
| `relative_solution_delta` values observed | `0` for all 8 candidates and both baselines |
| `status.toml` not "completed" | none found (all 10 result dirs: 8 candidates + 2 baselines = `completed`) |
| `eligible = false` | none found (all 10 `summary.csv` prepared rows = `true`) |

## Tolerance Map

Calibrated tolerance per candidate, read from each `config.toml` (`tolerance` field):

| Inner | Leaf | Tolerance |
|---|---|---|
| 3 | 25 | 3.80292896308757e-07 |
| 3 | 50 | 2.3864048149973856e-07 |
| 3 | 100 | 3.479128881193055e-07 |
| 3 | 200 | 4.252870058799867e-07 |
| 5 | 25 | 1.5209757150978087e-07 |
| 5 | 50 | 1.2891762136325558e-07 |
| 5 | 100 | 2.5712316195808637e-07 |
| 5 | 200 | 3.9670158724694127e-07 |

## Profile

Profiled configuration (`profile-fgs-j64-b1/R4/fgs_a696cc57f77e32a6/j64_b1/config.toml`): inner=3, leaf=100, P=8, MAC=0.4, tolerance=3.479128881193055e-7, cache_leaf_lu=true, max_iterations=300, rlx=1.0, sweep_order="lexicographic", rung="R4" — same hash/config as the seed baselines.

Total sample count (`cpu_flat.txt` header for the profiled Task): **Thread 1, Task 0x00007f16798a0010, "Total snapshots: 6741. Utilization: 98%"**. (Other threads in the file are GC/idle worker threads with 0–1 snapshots each and 100% "utilization" while parked; not the profiled compute task.)

Top 5 flat frames by the `Count` column within Thread 1's block (`cpu_flat.txt`, lines 1–4500), ties broken by file order (stable sort):

| Count | Overhead | File | Line | Function |
|---|---|---|---|---|
| 6738 | 1 | — | — | `[any unknown stackframes]` |
| 6730 | 0 | @Base/Base.jl | 562 | `include(mod::Module, _path::String…)` |
| 6730 | 0 | @Base/boot.jl | 430 | `eval` |
| 6730 | 0 | @Base/client.jl | 531 | `_start()` |
| 6730 | 0 | @Base/client.jl | 323 | `exec_options(opts::Base.JLOptions)` |

(Caveat: the `Count` column is inclusive/whole-stack count, so the top rows are all the same linear call chain from program entry to `solve!`, tied at 6730 — not differentiating hot compute. For reference, the file's `Overhead` column — self-time — separately shows `dgemv_kernel_4x4` at 5198/5198 as the largest self-time frame, with `setindex!`@array.jl:987 at 239/221 next; reported here as raw supplementary data, not as a substitute top-5.)

Top 3 allocation sites from `allocations.txt` (sample_rate=0.01, sampled_bytes=1157344), by leading count column, innermost application frame:

| Sampled count | Innermost frames (from stack, in order) |
|---|---|
| 6704 | `similar` @abstractarray.jl:822 [inlined] → `getindex` @array.jl:938 [inlined] → `influence!(...)` @FastMultipole solve.jl:1336 |
| 6704 | (same site) `influence!(...)` @FastMultipole solve.jl:1336 |
| 6072 | (same site) `influence!(...)` @FastMultipole solve.jl:1336 |

`unprofiled_trial.csv` row (single trial, profiling excluded):

```
setup_seconds,solve_seconds,total_seconds,allocated_bytes,gc_seconds,retained_bytes,process_peak_rss_bytes,iterations,estimated_inner_sweeps,estimated_fmm_passes,work_count_kind,formulation_subsolves,solved,fmm_rel_l2,fmm_rel_max,fmm_certified,fmm_seconds,epsilon_requested,direct_rel_l2,evaluator_delta,direct_seconds,authoritative_evaluator,authoritative_rel_l2,accepted,eligible
0,17.2803827,17.2803827,642597760,0.085711203,3221617552,4551806976,27,81,28,estimates; -1 means unavailable,0,true,4.77952513e-07,2.32095749e-05,true,5.17661351,1.22339544e-09,NaN,NaN,0,certified_fmm,4.77952513e-07,true,true
```

## Provenance

`selected.sha256` (root of opt-13660643):
```
092cd29526572e89658406597b378a91ef08d242bebf971e978f501b1342f2f4  /home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/opt-13657404/smoke-j4-b1/selected.toml
```
(Note: this hash/path points at the *prior* job's — 13657404 — `smoke-j4-b1/selected.toml`, not a file inside this job's own tree; reported verbatim, not reconciled.)

`screen-j64-b1/screen_bases.toml` (path: `.../opt-13660643/screen-j64-b1/screen_bases.toml`) — two saved bases, tolerances as calibrated at the time the bases file was written:

| Base | inner | leaf | tolerance |
|---|---|---|---|
| 1 | 3 | 100 | 3.479128881193055e-07 |
| 2 | 5 | 100 | 2.5712316195808637e-07 |

Comment header in the file states these are "Two median winners from opt-13657404 inner screen; tolerances are their calibrated values (harness resets to 0 and recalibrates every roster point)."

Package pins (`screen-j64-b1/provenance.toml` and root `campaign_pins.toml`, consistent):
- FLOWPanel: sha `f03ab18a7841e9a38876603a48b1799d83fd1e41`, tag `campaign/p021-cold-exec-20260912-v9`
- FastMultipole: sha `ef10643a401d6da16e28be87805b67d11bdf1fb5`, tag `campaign/p021-cold-exec-20260910-v1`
- FLOWVPM: sha `05c658f7804ec5f9b68d4cb9826a9f97cfecb373`, tag `campaign/p021-cold-exec-20260910-v1`

Other `provenance.toml` fields (screen-j64-b1): hostname=`m12-1-25`, julia_version=`1.11.7`, julia_threads=64, requested_blas_threads=1, blas=`LBTConfig([ILP64] libopenblas64_.so)`, filament_reg=`LineGaussRegularization`, timing_scope=`frozen _solve!; reset and BC diagnostics excluded; no formulation subsolves`, minimum_reps=10, screen_base_sha256=`7b343e2a0715a22305a08227265ea25aab7c3f5e740d3221772aace300244943`, manifest_sha256=`03cdb054fc315cb9a48cbb0ddf0331235a1b41484b4ac7889f4787b2cb64ce67`.
