# Job 13663309 harvest (R4 finalist confirmation)

## Job Facts

| Field | Value |
|---|---|
| Job ID | 13663309 |
| Node | m12-1-17 |
| Elapsed | 2h56m27s |
| Stage | screen_profile (baseline-j4-b1, baseline-j64-b1, screen-j64-b1, profile-fgs-j64-b1) |
| Panel count | 58,192 |
| Rung | R4 |
| CONFIG_FILE | `confirm_r4_20260912.toml` (two configs: finalist A inner=3, finalist B inner=5; leaf=100/P=8/MAC=0.4 both) |
| COMPLETED marker | `completed` (verbatim contents of `COMPLETED`) |

## Confirmation table

Per (stage, config), read from each config-dir's `config.toml` (inner/tolerance), `summary.csv` (prepared row), data-row count of `convergence.csv` (outer_iterations), and `convergence_validation.csv` row 1 (authoritative_rel_l2 / fmm_certified / relative_solution_delta). All `status.toml` = `completed`.

| Stage | Inner | Tolerance | Median (s) | Min (s) | Max (s) | Spread (s) | Reps | Eligible | Outer Iter | Authoritative rel-L2 | FMM certified | Rel. solution delta | Status |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| baseline-j4-b1 | 3 | 3.479128881193055e-07 | 16.897097 | 16.7756771 | 17.5015255 | 0.725848343 | 10 | true | 28 | 4.77952513e-07 | true | 0 | completed |
| baseline-j4-b1 | 5 | 2.5712316195808637e-07 | 14.6531677 | 14.5911795 | 15.2713502 | 0.680170702 | 10 | true | 18 | 3.65840408e-07 | true | 0 | completed |
| baseline-j64-b1 | 3 | 3.479128881193055e-07 | 11.4016262 | 11.2975517 | 11.9151186 | 0.617566841 | 10 | true | 28 | 4.77952513e-07 | true | 0 | completed |
| baseline-j64-b1 | 5 | 2.5712316195808637e-07 | 11.3566549 | 11.1579443 | 11.9980002 | 0.840055903 | 10 | true | 18 | 3.65840408e-07 | true | 0 | completed |
| screen-j64-b1 | 3 | 3.479128881193055e-07 | 10.8483714 | 10.7304166 | 11.3270652 | 0.596648574 | 10 | true | 28 | 4.77952513e-07 | true | 0 | completed |
| screen-j64-b1 | 5 | 2.5712316195808637e-07 | 11.1291186 | 10.977071 | 11.6023793 | 0.625308245 | 10 | true | 18 | 3.65840408e-07 | true | 0 | completed |

(Identical to `confirm_table.csv` in this directory.)

## Gates

| Gate | Result |
|---|---|
| Worst `authoritative_rel_l2` across all 6 rows | 4.77952513e-07 (inner=3, all three stages) |
| `fmm_certified` values observed | `true` for all 6 rows (no `false` observed) |
| `relative_solution_delta` values observed | `0` for all 6 rows |
| `status.toml` not "completed" | none found (all 6 result dirs = `completed`) |
| `eligible = false` | none found (all 6 `summary.csv` prepared rows = `true`) |

## Screen recalibrated tolerances

Read from `screen-j64-b1/R4/fgs_*/j64_b1/config.toml` (`tolerance` field, the harness-recalibrated value at screen time):

| Finalist | Inner | Recalibrated tolerance (screen-j64-b1) | Input tolerance (confirm_r4_20260912.toml) | Equal? |
|---|---|---|---|---|
| A | 3 | 3.479128881193055e-07 | 3.479128881193055e-07 | yes |
| B | 5 | 2.5712316195808637e-07 | 2.5712316195808637e-07 | yes |

Both finalists' recalibrated tolerances are exactly identical (bit-for-bit as printed) to the input tolerances.

## Profiles (profile-fgs-j64-b1)

Profiled config dirs: `R4/fgs_a696cc57f77e32a6/j64_b1` (inner=3) and `R4/fgs_6db6c2eac3662cd9/j64_b1` (inner=5); both config.toml otherwise leaf=100, P=8, MAC=0.4, cache_leaf_lu=true, max_iterations=300, rlx=1.0, sweep_order="lexicographic", rung="R4".

### Config: inner=3 (fgs_a696cc57f77e32a6)

- Profiled task total snapshot count (`cpu_flat.txt`, task with largest header total): **6726** (task `0x00007fce18e28010`, utilization 99%).
- Top self-time frame by the `Overhead` column: `dgemv_kernel_4x4` (…lia/libopenblas64_.so), Count=5065, Overhead=5065.
- Top 3 allocation sites (`allocations.txt`, sample_rate=0.01, sampled_bytes=1160168), by leading count, innermost application frames:
  1. count=6704: `similar`@array.jl:372[inlined] → `similar`@abstractarray.jl:822[inlined] → `getindex`@array.jl:938[inlined] → `influence!`@solve.jl:1336
  2. count=3880: same site (`influence!`@solve.jl:1336)
  3. count=2680: same site (`influence!`@solve.jl:1336)
- `unprofiled_trial.csv` row (single trial, profiling excluded):
```
setup_seconds,solve_seconds,total_seconds,allocated_bytes,gc_seconds,retained_bytes,process_peak_rss_bytes,iterations,estimated_inner_sweeps,estimated_fmm_passes,work_count_kind,formulation_subsolves,solved,fmm_rel_l2,fmm_rel_max,fmm_certified,fmm_seconds,epsilon_requested,direct_rel_l2,evaluator_delta,direct_seconds,authoritative_evaluator,authoritative_rel_l2,accepted,eligible
0,16.7666963,16.7666963,642687104,0.091191494,3221617552,4533350400,27,81,28,estimates; -1 means unavailable,0,true,4.77952513e-07,2.32095749e-05,true,5.33970831,1.22339544e-09,NaN,NaN,0,certified_fmm,4.77952513e-07,true,true
```

### Config: inner=5 (fgs_6db6c2eac3662cd9)

- Profiled task total snapshot count: **6842** (task `0x00007fce18e28010`, utilization 99%).
- Top self-time frame: `dgemv_kernel_4x4` (…lia/libopenblas64_.so), Count=5318, Overhead=5318.
- Top 3 allocation sites (`allocations.txt`, sample_rate=0.01, sampled_bytes=843768), innermost application frames:
  1. count=3880: `similar`@array.jl:372[inlined] → `similar`@abstractarray.jl:822[inlined] → `getindex`@array.jl:938[inlined] → `influence!`@solve.jl:1336
  2. count=3880: same site (`influence!`@solve.jl:1336)
  3. count=2280: `zeros`@array.jl:589[inlined] → `zeros`@array.jl:585 → `(::FastMultipole.var"#451#457"{Int64, DataType})(::Int64)`@none:0
- `unprofiled_trial.csv` row:
```
setup_seconds,solve_seconds,total_seconds,allocated_bytes,gc_seconds,retained_bytes,process_peak_rss_bytes,iterations,estimated_inner_sweeps,estimated_fmm_passes,work_count_kind,formulation_subsolves,solved,fmm_rel_l2,fmm_rel_max,fmm_certified,fmm_seconds,epsilon_requested,direct_rel_l2,evaluator_delta,direct_seconds,authoritative_evaluator,authoritative_rel_l2,accepted,eligible
0,11.3042064,11.3042064,97124440,0.016285182,3221617552,5110079488,17,85,18,estimates; -1 means unavailable,0,true,3.65840408e-07,1.49999495e-05,true,3.26854223,1.22339544e-09,NaN,NaN,0,certified_fmm,3.65840408e-07,true,true
```

(Caveat: `cpu_flat.txt`'s `Count` column is inclusive/whole-stack, so the linear entry-point chain ties at high values there; `Overhead` (self-time) is the differentiating column and is reported here, consistent with the format template's caveat.)

## Provenance

`selected.sha256` (root of opt-13663309), verbatim:
```
677c13773dcf3c74a84356753426dc6ad2d2f00845b86659d6bde418eebe151f  /home/rander39/campaigns/p021-cold-opt-20260912-v9/confirm_r4_20260912.toml
```

`screen-j64-b1/screen_bases.toml` — two saved bases (comment: "Two median winners from opt-13657404 inner screen; tolerances are their calibrated values (harness resets to 0 and recalibrates every roster point)."):

| Base | inner | leaf | tolerance |
|---|---|---|---|
| 1 | 3 | 100 | 3.479128881193055e-07 |
| 2 | 5 | 100 | 2.5712316195808637e-07 |

Package pins (`screen-j64-b1/provenance.toml` and root `campaign_pins.toml`, consistent):
- FLOWPanel: sha `f03ab18a7841e9a38876603a48b1799d83fd1e41`, tag `campaign/p021-cold-exec-20260912-v9`
- FastMultipole: sha `ef10643a401d6da16e28be87805b67d11bdf1fb5`, tag `campaign/p021-cold-exec-20260910-v1`
- FLOWVPM: sha `05c658f7804ec5f9b68d4cb9826a9f97cfecb373`, tag `campaign/p021-cold-exec-20260910-v1`

Other `provenance.toml` fields (screen-j64-b1): hostname=`m12-1-17`, julia_version=`1.11.7`, julia_threads=64, requested_blas_threads=1, blas=`LBTConfig([ILP64] libopenblas64_.so)`, filament_reg=`LineGaussRegularization`, timing_scope=`frozen _solve!; reset and BC diagnostics excluded; no formulation subsolves`, minimum_reps=10, screen_base_sha256=`7b343e2a0715a22305a08227265ea25aab7c3f5e740d3221772aace300244943`, manifest_sha256=`03cdb054fc315cb9a48cbb0ddf0331235a1b41484b4ac7889f4787b2cb64ce67`.
