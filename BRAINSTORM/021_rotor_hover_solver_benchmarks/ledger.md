# 021 ledger — single running results table

One row per accepted result set; superseded rows are struck through, never deleted. Raw CSVs
live in `benchmark/results/` / `data/`; this table holds the headline numbers.

## Frozen mesh ladder (set in Phase 1, then immutable)

**RE-FROZEN 2026-08-13 evening on CAPPED (captess4) meshes per Ryan's formulation
ruling: the campaign case matches `examples/rotor_hover_pressure_comparison.jl`
(Dirichlet — watertight capped meshes, Union{ConstantSource,VortexRing}, DBC=true,
kerneloffset_panel=R·1e-10 / targets=1e-3, rho 1.179, RPM 6000).** The same-day
uncapped freeze (struck through below) is superseded; the uncapped meshes remain on
disk for Neumann diagnostics. n_panels = blade 4(U−1)(W−1) + cap cells (measured;
verified watertight via the TE-detector topology assert, full TE traced with the
RHPC end_node-anchored recipe + root clip at r/R 0.1 + cap-wrap guard).
Constant-aspect-ratio family anchored on the workhorse 45_145: scale (U−1, W−1)
together by √2 per rung so N = 4(U−1)(W−1) doubles per rung and panel aspect ratio
stays fixed (the old 30_xx/45_xxx candidates drift AR 2.76→4.18 — rejected for that
reason). Tess_W odd (script-enforced ⇒ even chordwise panel count ⇒ upper/lower
control points mirror across the section plane). Uncapped, stock hub.
**OpenVSP rounding constraint (measured 2026-08-13)**: the chordwise panel count is
rounded UP to a multiple of 8, so honest tessellation requires W−1 ≡ 0 (mod 8)
(every prior DJI mesh satisfies this). The original R2 31_103 / R4 63_205 picks
violated it (produced 12,480 / 51,584 cells under mislabeled names — regenerated as
33_105 / 65_209 with U re-balanced to hold the anchor AR (W−1)/(U−1) ≈ 3.27; AR now
3.25–3.29 across all rungs). All rungs generated 2026-08-13 and verified: exact
N = 4(U−1)(W−1), watertight-topology check via TE detection, full-span shedding
(U−1 edges per blade) traced on the constructed body.
**Cost ceiling (set 2026-08-13 from the measured R1–R3 local cost check, phase_01
Log): matrix-ful (Backslash) configs run R1–R5 and drop out at R6–R7.** Measured
(uncapped Neumann) 4-thread assembly/LU scale as ~N^2.0 / ~N^2.9; at the capped
counts the projections are R5 = 94 GiB dense G (fits a large HPC node), R6 =
360 GiB (does not) — the R1–R5 / R6–R7 split is unchanged; Dirichlet R1–R3 costs
re-measured in the capped floor runs. Drop rows marked per config, not per rung.

| Rung | Tess (U_W) | Mesh file | n_panels (blade + caps) | AR (W−1)/(U−1) |
| --- | --- | --- | --- | --- |
| R1 | 23_73 | `dji9443_20260813_23_73_capped_captess4.msh` | 8,016 (6,336 + 1,680) | 3.27 |
| R2 | 33_105 | `dji9443_20260813_33_105_capped_captess4.msh` | 15,760 (13,312 + 2,448) | 3.25 |
| R3 | 45_145 | `dji9443_20260813_45_145_capped_captess4.msh` | 28,752 (25,344 + 3,408) | 3.27 |
| R4 | 65_209 | `dji9443_20260813_65_209_capped_captess4.msh` | 58,192 (53,248 + 4,944) | 3.25 |
| R5 | 89_289 | `dji9443_20260813_89_289_capped_captess4.msh` | 108,240 (101,376 + 6,864) | 3.27 |
| R6 | 125_409 | `dji9443_20260813_125_409_capped_captess4.msh` | 212,108 (202,368 + 9,740) | 3.29 |
| R7 | 177_577 | `dji9443_20260813_177_577_capped_captess4.msh` | 419,276 (405,504 + 13,772) | 3.27 |

<details><summary>Superseded same-day uncapped (Neumann) freeze — struck through, meshes kept on disk</summary>

| Rung | Tess (U_W) | Mesh file | n_panels | AR |
| --- | --- | --- | --- | --- |
| ~~R1~~ | ~~23_73~~ | ~~`dji9443_20260813_23_73_uncapped.msh`~~ | ~~6,336~~ | ~~3.27~~ |
| ~~R2~~ | ~~33_105~~ | ~~`dji9443_20260813_33_105_uncapped.msh`~~ | ~~13,312~~ | ~~3.25~~ |
| ~~R3~~ | ~~45_145~~ | ~~`dji9443_20260813_45_145_uncapped.msh`~~ | ~~25,344~~ | ~~3.27~~ |
| ~~R4~~ | ~~65_209~~ | ~~`dji9443_20260813_65_209_uncapped.msh`~~ | ~~53,248~~ | ~~3.25~~ |
| ~~R5~~ | ~~89_289~~ | ~~`dji9443_20260813_89_289_uncapped.msh`~~ | ~~101,376~~ | ~~3.27~~ |
| ~~R6~~ | ~~125_409~~ | ~~`dji9443_20260813_125_409_uncapped.msh`~~ | ~~202,368~~ | ~~3.29~~ |
| ~~R7~~ | ~~177_577~~ | ~~`dji9443_20260813_177_577_uncapped.msh`~~ | ~~405,504~~ | ~~3.27~~ |

</details>

## Frozen per-solver settings (Phase 1 calibration)

| Solver config | Settings | Matched RMS residual achieved |
| --- | --- | --- |
| — | — | — |

## Results

| Date | Phase | Solver config | Rung | Threads (Julia/BLAS) | Setup [s] (asm/fact/tree/precond) | Per-step [s] | Iters | RMS / max residual | Mem | Notes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 2026-08-17 | 2b local cache A/B | FGS uncached | R1 | 4/1 | total 11.588; leaf factor 0 | 1.19617 | 11 | certified BC rel-L2 2.83536e-7 | 249,430,448 B | Phase-1 frozen τ=1e-6 knobs; cold min-of-5; leaf pass 0.0129974 s / 11,554,848 B alloc |
| 2026-08-17 | 2b local cache A/B | FGS cached | R1 | 4/1 | total 9.656; leaf factor 0.00977 | 0.455358 | 11 | certified BC rel-L2 2.83536e-7 | 260,097,472 B | 2.63x solve; leaf pass 0.00046275 s / 0 B alloc; cache 10,656,608 B |
| 2026-08-17 | 2b local cache A/B | FGMRES+FGS uncached | R1 | 4/1 | total 11.633; leaf factor 0 | 3.30494 | 1 | certified BC rel-L2 8.64628e-7 | 257,079,680 B | Phase-1 frozen τ=1e-6 knobs; cold min-of-5; leaf pass 0.0103698 s / 11,554,848 B alloc |
| 2026-08-17 | 2b local cache A/B | FGMRES+FGS cached | R1 | 4/1 | total 9.657; leaf factor 0.01110 | 2.22390 | 1 | certified BC rel-L2 8.64628e-7 | 267,746,704 B | 1.49x solve; leaf pass 0.000673542 s / 0 B alloc; cache 10,656,608 B |
| 2026-08-19 | 2b local nearfield-cache A/B (local-smoke) | krylov_gmres cache_tree, nearfield uncached | R1 | 4/4 | — | 130.969 (cold min-of-3) | 83 | certified BC rel-L2 8.47471e-9 | — | Phase-1 frozen apply knobs p17/mac0.5/leaf21; atol 1e-10 / rtol 1e-8; isolated near-field pass 1.43309 s (min-of-7) |
| 2026-08-19 | 2b local nearfield-cache A/B (local-smoke) | krylov_gmres cache_tree + cache_nearfield | R1 | 4/4 | nearfield-cache build 11.892 | 23.7209 (cold min-of-3) | 83 | certified BC rel-L2 8.47471e-9 | cache 272,389,024 B | **5.52x solve, 265x isolated near-field** (0.0054061 s); break-even 8.3 applies; solution shift rel-L2 3.5e-15; apply exactness 4.3e-16; dense-G baseline 8N²=514 MB |
| 2026-08-19 | 2b nfcache configs (local-smoke) | krylov_gmres_nfcache | R1 | 4/4 | ctor 0.14 | 33.794 (cold min-of-5) | 59 | certified BC rel-L2 9.673e-7 | cache 810,094,272 B/solve | cached shared knobs p14/mac0.5/leaf326 (tune_cached.csv); cache build 31.7 s INSIDE every solve (per-solve state) — build-dominated |
| 2026-08-19 | 2b nfcache configs (local-smoke) | krylov_ilu_nfcache | R1 | 4/4 | ILU 3.81 | 33.050 (cold min-of-5) | 7 | certified BC rel-L2 8.295e-7 | cache 810,094,272 B/solve | same shared knobs; build-dominated — ILU's 59→7 iter saving worth <1 s at cached-apply cost |
| 2026-08-19 | 2b nfcache configs (local-smoke) | fgmres_fgs_nfcache | R1 | 4/4 | precond 12.27 | 34.115 (cold min-of-5) | 1 | certified BC rel-L2 5.884e-7 | cache 810,094,272 B/solve | Stage-3 winner precond knobs unchanged; build-dominated ⇒ lever = staged rigid-motion cache persistence |

## Phase-1 finding: the Barba ILU preconditioner scales ~N^1.5, not ~N (2026-08-22)

Measured on this machine (4 threads, geometry-only probe over the frozen
ladder; leaf_size=10, MAC=1.0 — the ILU knobs every 021 benchmark script uses).
`nnz` equals the direct-list pattern entry count exactly at every rung.

| Rung | N | Pattern entries (= nnz) | Entries/row | Construct [s] | Apply `ldiv!` [ms] | Apply [ns/nnz] |
| --- | --- | --- | --- | --- | --- | --- |
| R1 | 8,016 | 3,906,740 | 487.4 | 3.80 | 2.31 | 0.59 |
| R2 | 15,760 | 8,002,490 | 507.8 | 7.12 | 4.85 | 0.61 |
| R3 | 28,752 | 19,150,596 | 666.1 | 16.03 | 11.30 | 0.59 |
| R4 | 58,192 | 48,433,784 | 832.3 | 53.53 | 29.52 | 0.61 |
| R5 | 108,240 | 130,925,412 | 1,209.6 | 215.39 | 92.21 | 0.70 |
| R6 | 212,108 | 436,830,804 | 2,059.5 | — | ~284 (projected) | — |
| R7 | 419,276 | 1,547,431,218 | 3,690.7 | — | ~1,000 (projected) | — |

R6/R7 pattern sizes are measured; their apply times are projected at the
measured ~0.65 ns/nnz (apply is linear in nnz to within 2% across R1–R4, 18% at
R5), because building the factors needs more RAM than this machine has.

**Scaling.** Entries per row grow as ~sqrt(N), so nnz ~ N^1.5. Both costs follow
it: construction R1→R5 is N x13.5 -> time x56.7 (exponent 1.55), apply is
N x13.5 -> time x39.9 (exponent 1.44). Apply is linear in nnz at a near-constant
~0.6 ns/nonzero, so the exponent is inherited from the pattern, not from the
triangular solve.

**Consequence for the campaign.** An FMM operator apply grows ~N log N while
this preconditioner's apply grows ~N^1.5, so the ILU share of every GMRES
iteration rises without bound along the ladder. At R7 one apply is ~1 s and the
pattern alone is ~23 GiB before factors. The ILU rows in the Phase-1 table
should therefore be read as rung-local results, not as a solver that carries to
production mesh sizes; the crossover against the unpreconditioned/Jacobi arms is
itself a Phase-1 deliverable worth stating explicitly.

**Operational note.** `max_pattern_entries` was 2048*N in all seven benchmark
scripts, set at R3. R6 needs 2,059.5 entries/row — 0.6% over — which is what
killed `p1-table-R6-multi` (job 13206077) after 4 h 09 m. Raised to 8192*N
(clears R7 with ~2.2x headroom); the limit is a guard only, nothing is allocated
from it. The guard now reports the TOTAL entries required and entries/row
instead of the running subtotal at which it tripped, so it can be sized in one
shot. NOTE: the unit test "ILU(0) construction ladder remains linear"
(test/runtests_unit_solver.jl:581) does not actually enforce linearity — it
allows 2.5x nnz growth across ~2.1x N steps (i.e. up to N^1.28) on 128–512-panel
spheres, a regime that does not reproduce the real meshes' N^1.5.

## FINDING: stock `tune_fmm` picks a ~2x-too-expensive leaf_size (2026-08-24)

**The new FastMultipole is NOT slower. The tuner is choosing badly.** This was
nearly mis-attributed as a FastMultipole regression, so the evidence is recorded
in full.

**Symptom.** R1 multi under FM `d8258a7d` came out 1.5-2.5x slower than under
`a9b734a` on every FMM path, at bit-comparable accuracy:

| config | solve old | solve new | delta | niter o/n |
| --- | --- | --- | --- | --- |
| krylov_gmres | 14.59 s | 36.26 s | +148% | 59/59 |
| krylov_ilu | 1.718 s | 4.333 s | +152% | 7/7 |
| krylov_jacobi | 4.774 s | 11.92 s | +150% | 19/19 |
| fgs (target_1e-6) | 0.587 s | 1.115 s | +90% | 8/7 |
| fgmres_fgs | 0.746 s | 1.081 s | +45% | 1/1 |
| **backslash_ldiv (no FMM in solve)** | **0.0378 s** | **0.0362 s** | **-4%** | — |

Iteration counts identical, certified BC agreeing to 4 significant figures
(`krylov_gmres` 9.674e-7 -> 9.675e-7), and the one config whose solve contains no
FMM apply unchanged. Phase 2a R1 reproduced it independently. Only knob that
moved: apply `leaf_size` 9 -> 60 (`p` and `mac` identical at 17 / 0.5).

**Diagnostic** (`benchmark/fm_leaf_ab.jl`, job 13407693, R1, current FM, fixed
p=17 mac=0.5):

| leaf | t_solve_min | niter |
| --- | --- | --- |
| **9** | **15.51 s** | 57 |
| 21 | 20.65 s | 57 |
| 40 | 26.10 s | 55 |
| 60 | 32.64 s | 57 |
| 71 | 34.95 s | 56 |

Monotone in leaf, 2.25x across the range at constant iteration count. **At
matched knobs (leaf 9) the two FM generations agree to +6% (14.59 -> 15.51 s).**
That +6% is the true generation delta; the remaining ~140% was the knob.
(NB: the `bc_rel_l2` column of that diagnostic is invalid — wrong arguments were
passed to `bc_error!` — so accuracy parity rests on the campaign rows above, not
on the sweep.)

**Root cause.** Stock `FastMultipole.tune_fmm` optimizes a SINGLE evaluation
against an error tolerance. It does not measure the cost of an iterative solve,
so it will happily accept a large leaf that meets the accuracy target while being
~2x more expensive per apply — which is precisely the quantity Phase 1 exists to
measure. Under `a9b734a` it happened to land on leaf 9; under `d8258a7d` it lands
on 60-71.

**Fix available in-tree.** `FastMultipole.tune_fmm_perturb`
(`src/autotune.jl:337`) is a coordinate-descent tuner over
`(expansion_order, multipole_acceptance, leaf_size_source)` that selects on
MEASURED wall-clock (`reps=2`, `improve_tol=0.02`, `leaf_factor=1.5`) subject to
an `error_tolerance` constraint. It requires a starting point and errors if the
seed misses the tolerance, so the robust composition is: stock `tune_fmm` to find
an accuracy-valid seed, then `tune_fmm_perturb` from that seed to descend on
cost. That needs no prior knobs, so it also works for R7, which has never been
tuned.

**Campaign impact.** All Phase-1 tune jobs were cancelled ~40 min in rather than
spend up to 9 node-days producing knobs that are ~2x off; 21 dependent
table/2a/2b jobs were HELD (not cancelled) so nothing consumes them. The R1 rows
above stand as the demonstration of the effect and must not be treated as
accepted results.

## GENERATION BREAK: new FMM + isolated checkout (2026-08-24)

**Every row above this line is superseded.** The campaign restarts the full
ladder under a new FastMultipole and a new, isolated checkout. Old rows are kept
for the record, never deleted, but must not be plotted on the same axes as
anything produced after this date.

**What changed and why.** The 2026-08-24 FLOWPanel source (021 Phase-3 solver
lifecycle + the Backslash `S` seam) references nine FastMultipole symbols that do
not exist in `a9b734a` — `NEARFIELD_CACHE_DEFAULT_MAX_BYTES` (used at
`src/FLOWPanel_solver.jl:504,644`, i.e. the solver every 021 driver runs),
`FmmPlan`, `NearfieldInfluenceCache`, `transform_plan!`, `transform_solver!`,
`direct_rectangular!`, `RectangularGaussianErfVortex`,
`RectangularPanelInfluence`, `load_cuda_radix_lifecycle!`. Continuing on the old
FMM was therefore not an option: it was a runtime failure, not a preference.

**Provenance, which was the real hazard.** The old pin `fm_commit =
a9b734ad…-dirty` was an *uncommitted* worktree whose `src/autotune.jl` blob
(`47beaa59`) existed in no git history anywhere — verified with `git cat-file`
and `--find-object`. Worse, when the new FastMultipole content was first copied
onto the cluster it was laid on top of the same commit, so two entirely different
FMM generations would both have recorded the string `a9b734ad…-dirty`. All three
repos were then committed and pushed, and the campaign moved to clean clones, so
`fm_commit` now resolves to a real object.

| | before | after |
| --- | --- | --- |
| FastMultipole | `a9b734a`-dirty (unversioned blob) | **`d8258a7d`, clean** |
| FLOWPanel | working tree via rsync | **`321473f`, clean** |
| FLOWVPM | `a950790` (Jun 30, dirty) | **`4494c25`, clean** |
| checkout | `~/projects/FLOWPanel.jl` (shared with 018/023) | **`~/flowpanel-021/`, isolated** |
| hardware | m12/zen3, flags typed by hand | m12 / physics2 / m12pws, pinned in-header |

The old pinned state is archived twice: `fm_a9b734a_dirty_021pinned_20260824.patch`
on the cluster, and a verified 5-file copy taken before the overwrite.

**Knobs are re-tuned, not reused.** `tune_fmm` tunes FMM parameters, so its
output is FM-dependent by construction. Reusing knobs selected under `a9b734a`
would leave the ladder inconsistent in exactly the dimension this break is meant
to fix, so the tune chain is re-run for every rung R1–R7.

**Hardware pin.** Three partitions are spec-identical (2 sockets x 64 zen3 cores,
`ThreadsPerCore=1`, `CPUTot=128`, `RealMemory=524288 MB`): `m12` (public, 3 d),
`physics2` (7 d) and `m12pws` (28 d), the latter two via `--qos=standby`. All
seven launchers pin `--constraint=zen3 --cpus-per-task=64 --mem=250G` in the
SBATCH header and assert, at runtime, both `SLURM_CPUS_ON_NODE == 64` (this job
holds the 64 CPUs its 64 Julia threads expect) and
`getconf _NPROCESSORS_ONLN == 128` (the node is the pinned zen3 part). Hardware
can no longer depend on what someone typed at the submit line. `HARDWARE_TAG` is
the stable `orc-m12-zen3-2x64c-512g`; pre-break rows carry `m12-2-30` or
`zen3-exclusive`, which denote the same verified hardware.

**HALF-NODE ALLOCATION (Ryan ruling 2026-08-24), and its cost.** Jobs take 64
CPUs, not a whole 128-core node. Julia runs 64 threads either way, so
`--exclusive` was wasting half of every node and forcing each job to wait for a
fully free one; at 64 CPUs two jobs pack per node and 24 of 32 started
immediately instead of queuing ~15 h. **The trade-off is recorded deliberately:
these partitions are `OverSubscribe=YES:4`, so a co-tenant can now share the node
and contend for memory bandwidth.** The FMM/ILU workload is bandwidth-bound, so
co-tenancy is a real (if usually small) contaminant on every timing in this
generation — uniform across rungs, so the ladder stays internally comparable, but
NOT zero, and the 1-thread `MODE=single` arms are the most exposed. Any
cross-generation comparison against the pre-break `--exclusive` rows must account
for it.

**Gotcha, cost one full fan-out.** Do not test CPU allocation with `nproc` in
these scripts: GNU `nproc` honours `OMP_NUM_THREADS`, and Slurm presets that to 1
whenever `--cpus-per-task` is used, so `nproc` reports 1 while the job's real
affinity is the full `0-63` (verified via `/proc/self/status Cpus_allowed_list`).
Use `SLURM_CPUS_ON_NODE`. Equally, `--ntasks=64 --cpus-per-task=64` requests
64x64 = 4096 CPUs; the correct spec for one 64-thread process is
`--nodes=1 --ntasks=1 --cpus-per-task=64`.

**Walltimes are right-sized per job, not blanket.** On a priority-0 `standby`
QOS, Slurm backfills work that fits an existing gap, so an oversized request
schedules far worse than an accurate one. Limits are set from the measured
2026-08-18/19 timings with ~2x margin (R1 table 3.6 min, R2 12 min, R3 20 min,
R4 48 min, R5 2 h 08) rather than a uniform 7 days.

**Resume support (required by preemption).** `standby` is preemptible and the
launchers set `--requeue`, so a preempted job restarts from scratch. The three
appending drivers now skip work that already landed —
`rotor_hover_solver_phase1_table.jl` and `rotor_hover_solver_phase2.jl` on
`(rung, config, row_kind)`, and `..._phase2b_nearfield_cache.jl` on
`(rung, row_kind)` driving its existing `SKIP_AB` / `TUNE_CACHED` gates. Each row
is flushed as its measurement completes, so a row in the CSV means that
measurement is finished and is never repeated. Verified against real data: the
parser recovers R6's 2 partial rows and R1's complete 7 including the `fgs` dual
row.

## Data quarantine: pre-Phase-3 `phase1_costcheck` timings (2026-08-24)

The 021 Phase-3 solver-lifecycle fix (2026-08-23) added a one-body early break to
the tuple `solve!` (`src/FLOWPanel_solver.jl:2000-2011`): a single-body tuple has
no coupling to iterate, so the block-GS outer loop now exits after the first
block solve instead of re-solving the identical system every outer iteration.
Every `phase1_costcheck` / `phase0_smoke backslash_iterative` row measured before
that date therefore carries **two** inner solves in `t_solve_min` and
`alloc_solve_bytes`, and is not comparable to anything measured after it.

Those rows were moved out of the live files (rows preserved verbatim with header,
nothing deleted):

| Was | Now | Rows |
| --- | --- | --- |
| `benchmark/results/phase1/multi/runs.csv` | `benchmark/results/phase1/multi/runs_pre_phase3.csv` | 17 (all `phase1_costcheck`, `backslash_ref`, R1–R3) |
| `benchmark/results/smoke/multi/runs.csv` | `benchmark/results/smoke/multi/runs_pre_phase3.csv` | 1 (`backslash_iterative_multi`) |
| `benchmark/results/smoke/single/runs.csv` | `benchmark/results/smoke/single/runs_pre_phase3.csv` | 1 (`backslash_iterative_single`) |

The stale note strings that described the two-solve behaviour were corrected at
the same time: `benchmark/rotor_hover_solver_phase1.jl:291` and
`benchmark/rotor_hover_solver_smoke.jl:278-283,292`. No analysis script globs
`runs*.csv`, so the quarantine filename is inert; `smoke` rewrites its `runs.csv`
with `open(path, "w")` on every run, which is why the smoke rows had to be moved
rather than left in place.

Only the tuple path is affected. The other smoke configs and every
`phase1_table` / `phase1_tune` row go through the single-body `solve!` or raw
`_solve!`, which the same commit re-bracketed with the step lifecycle but did not
change numerically for a cold (`warmstart=false`, `solution_history_length=0`)
solve. The Dirichlet source-potential seam `_source_influence!`
(`src/FLOWPanel_solver.jl:384-395`) falls back to the original
`influence!(body, body, backend; scalar_potential=true)` whenever
`solver.S === nothing`, and `Backslash`'s `assemble_source_potential` kwarg
defaults to `false` — no benchmark driver opts in. **If a driver ever does opt
in, Phase-1 timings stop being comparable across that boundary.**

---

## Near-field cache feasibility map (W0, 2026-08-24)

`benchmark/nfcache_feasibility.jl` → `benchmark/results/phase2/multi/nfcache_feasibility.csv`.
MEASURED via `FastMultipole.estimate_nearfield_cache` (exact size-pass
arithmetic, nothing built, nothing solved) at p=17/MAC=0.5 — `bytes` is
independent of expansion order and hardware, so these numbers are exact and
locally valid. Cache size in GiB, with the ratio to a dense `8N²` operator in
parentheses.

| rung | N | dense 8N² [GiB] | leaf 9 | leaf 21 | leaf 45 | leaf 100 | leaf 225 | leaf 500 | leaf 1100 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| R1 | 8,016 | 0.5 | 0.20 (0.423×) | 0.25 (0.530×) | 0.31 (0.642×) | 0.41 (0.857×) | 0.64 (1.337×) | 0.92 (1.932×) | 0.96 (2.001×) |
| R2 | 15,760 | 1.9 | 0.52 (0.280×) | 0.64 (0.344×) | 0.80 (0.433×) | 1.13 (0.608×) | 1.68 (0.906×) | 2.49 (1.347×) | 3.58 (1.933×) |
| R3 | 28,752 | 6.2 | 1.29 (0.210×) | 1.62 (0.263×) | 2.00 (0.325×) | 2.81 (0.457×) | 3.87 (0.629×) | 6.16 (1.000×) | 9.26 (1.504×) |
| R4 | 58,192 | 25.2 | 3.16 (0.125×) | 3.96 (0.157×) | 5.03 (0.199×) | 7.19 (0.285×) | 10.76 (0.426×) | 15.42 (0.611×) | 27.19 (1.078×) |
| R5 | 108,240 | 87.3 | 8.18 (0.094×) | 9.97 (0.114×) | 12.83 (0.147×) | 17.11 (0.196×) | 25.37 (0.291×) | 37.54 (0.430×) | 64.15 (0.735×) |
| R6 | 212,108 | 335.2 | 23.83 (0.071×) | 27.50 (0.082×) | 33.47 (0.100×) | 45.51 (0.136×) | 67.00 (0.200×) | 103.31 (0.308×) | 157.28 (0.469×) |
| R7 | 419,276 | 1309.8 | 78.23 (0.060×) | 89.32 (0.068×) | 102.56 (0.078×) | 128.68 (0.098×) | 180.82 (0.138×) | 284.89 (0.218×) | 426.88 (0.326×) |

**Reading it.**

1. **The cache/dense ratio FALLS with N** — 0.53 → 0.068 at the kernel leaf
   across R1→R7. The dense degeneracy that motivated the memory axis is an
   R1–R2 artifact; matrix-free character genuinely returns at scale.
2. **A cache has a FLOOR.** Shrinking the leaf raises the block count without
   shrinking total near-field area, so bytes never fall below the small-leaf
   value: 0.20 GiB at R1 but **78 GiB at R7**, 23.8 GiB at R6. A budget under a
   rung's floor admits NO cached configuration at all. The prior estimate
   ("naive R7 ≈ 32 GiB", `phase_15_caching_and_objective_prompt.md`)
   underestimated by **2.4-2.8×** (78.2 GiB at leaf 9, 89.3 GiB at the kernel
   leaf 21).

   Cross-check on the one previously recorded datum: the handoff's "R4 ≈ 4.5
   GiB" sits between this map's leaf 30 (4.40 GiB) and leaf 45 (5.03 GiB) —
   the map reproduces an independently measured point.
3. **The cache asymptotes at ~2× dense, not 1×** (R1 leaf ≥750 → 2.00×). The
   direct list holds both (i,j) and (j,i), so symmetric blocks are stored
   twice. A dense solver is therefore strictly more memory-efficient at the
   dense end — a real point on the memory axis, not a rounding artifact.
4. **Budget-ladder pruning.** Budgets at or above ~2× a rung's dense size all
   admit the same saturated optimum; only budgets between the floor and that
   ceiling are distinct. R1's window is 0.2–1.0 GiB (narrow); R7's is
   78–427 GiB.

**Anchor check:** leaf 340 → 0.803 GiB reproduces the recorded 862104272 B at
leaf 342 (`nearfield_cache_ab.csv`) exactly.

## Superseded by the W5 corrections (2026-08-24)

| Row | Status |
| --- | --- |
| `nearfield_cache_ab.csv` `solve_krylov_gmres_cache_{on,off}` — the 5.5× cold-solve speedup | **UPPER BOUND.** The A/B had no warm-up rep, so the first arm (uncached) absorbed all the JIT, biasing it toward the cached arm. Warm-up added; re-measure before quoting. |
| `fm_leaf_ab.jl` `bc_rel_l2` column (~1.02–1.09) | **INVALID, cause corrected.** Not `bc_error!` arguments: the local `reset!` omitted the frozen-potential restore, and the Dirichlet rhs is `-rotor.potential`, which a solve clobbers — so every solve after the first ran against a stale b. Now uses `reset_cold!`. **Timings unaffected**; the R2 2.18× leaf result stands. |

## D4 additivity — cold = build + warm (R1, 2026-08-24)

Ryan asked whether the warm test alone suffices, since "a cold solve would
merely add the cache build time". MEASURED at R1, budget 2 GiB, p16/MAC0.5/
leaf21, min-of-3 each (`phase2.csv` `row_kind="additivity"`):

| quantity | [s] |
| --- | --- |
| `t_cold` (plan+cache rebuilt inside every timed solve) | 16.059 |
| `t_build` (cache build, serial) | 11.800 |
| `t_warm` (warm cache, cold solution) | 4.223 |
| `t_build + t_warm` | 16.024 |
| **relative discrepancy** | **0.22%** |

**Confirmed.** The warm measurement plus the build column recovers the cold
number, so the cached rows need only ONE timed pass, and ruling 12's cold
convention is preserved in the sense that matters (the SOLUTION is cold:
`reset_cold!` zeroes strengths, warmstart off, `niter=59` unchanged).

Consequence for the pre-W3 rows: folding the build into `t_solve` inflated the
reported per-step cost of every `*_nfcache` row by ~2.9x at R1
(4.22 s warm vs ~16 s cold). Those rows must not be compared against the new
ones. LOCAL timings — cluster re-measurement required before publication.

## Local artifacts quarantined (2026-08-24) — do NOT read as campaign data

| File | Why |
| --- | --- |
| `benchmark/results/phase2/multi/tune_phase2_localsmoke.csv` | Renamed out of the way. These are R1 knobs from a Mac smoke with `TUNE_REPS=1` and a 240 s cap that TIMED OUT — local timings scatter 22-39% against the effect being measured, so they are NOT campaign knobs. `rotor_hover_solver_phase2.jl` reads `tune_phase2.csv`; leaving these in place would have silently promoted them. |
| `benchmark/results/phase2/multi/phase2_w3_localsmoke.csv` | The W3 smoke rows (warm + additivity). Evidence that the driver works, not benchmark data. |
| `benchmark/results/phase2/multi/phase2.csv` | PRE-W3 schema (no `mem_budget_gib` / `nfcache_build_time` columns) and its `*_nfcache` rows fold the cache build into `t_solve`. The new driver's header guard REFUSES to append to it — move it aside before the retuned campaign runs. |
| `benchmark/results/phase2/multi/tune_cached.csv` | Accuracy-only objective (stock `tune_fmm`). Superseded by `tune_phase2.csv`; kept as evidence of record, not read by anything now. |

## Reference-free verification + total-memory axis — R1 local smoke (2026-08-24)

Mechanism checks for the two corrections. LOCAL timings are indicative only;
the `bc_rel_l2` and memory figures are not timings and are exact.

**Reference-free BC metric** (`rotor_hover_solver_phase1_agreement.jl`, FMM
evaluation certified at `0.1 x target_rel x rms_b`):

| config | bc_rel_l2 | certified | CT (laplace, bernoulli, kj) |
| --- | --- | --- | --- |
| krylov_gmres | 9.673e-7 | true | 0.060406, 0.053168, 0.057221 |
| fgs | 4.607e-8 | true | 0.060406, 0.053168, 0.057221 |

Cross-validation: `krylov_gmres` reads **9.673e-7** here and **9.673e-7** in
the independent Phase 2 tuner at the same rung and knobs — two different
drivers, same number.

Ensemble spread (reference-free, and it correctly picked up a solution
persisted by an EARLIER process):

| | relative L2 |
| --- | --- |
| fgs, deviation from ensemble mean | 6.206e-6 |
| krylov_gmres, deviation from ensemble mean | 6.206e-6 |
| **max pairwise (fgs vs krylov_gmres)** | **1.241e-5** |

Both solvers are individually correct to below the 1e-6 BC criterion and agree
with each other to 1.2e-5 in solution — no reference solve, no dense G, no
ordering constraint.

**BUG FOUND AND FIXED while smoking this**: the first run reported
`bc_rel_l2 = 0.985` for a converged `krylov_gmres` — "BC completely
unsatisfied". Cause: `steady!` runs post-processing monitors that OVERWRITE
`rotor.velocity` (and leave `core_size` at the TARGETS offset), and
`bc_error!` recomputes the BC source strengths from `body.velocity`, so sigma
no longer matched the frozen kinematic state. The driver now restores the
frozen velocity and the panel core_size before evaluating. Any BC number this
driver produced before this fix is meaningless. (Same family as the
`fm_leaf_ab.jl` stale-`b` defect above — post-processing clobbering solve-time
state is a recurring hazard in these drivers.)

**Total-memory budget axis** (`rotor_hover_solver_phase2_tune.jl`). A budget is
body + solver state + FMM plan + cache. Cached-or-not is an OUTCOME of what
fits:

| budget [GiB] | cached | t_warm [s] | cache [MB] | mem predicted [MB] | mem measured [MB] | note |
| --- | --- | --- | --- | --- | --- | --- |
| 0.05 | **no** | 127.3 | 0 | 14.7 | 15.3 | no cache fits -> tuned uncached automatically |
| 0.4 | yes | 8.3 | 272 | 287.1 | 289.2 | build 14.1 s |
| 4 | yes | 8.3 | 272 | 287.1 | 289.2 | flagged DUPLICATE - same operating point |

Predicted vs measured total agree to **0.7%**, so the a-priori accounting used
to reject over-budget candidates is sound. The 0.05 GiB row is the behaviour
Ryan asked for: a rung/budget that prohibits a cache is tuned without one, via
the same code path. The two lower rows show the tradeoff the axis exists to
expose — **19x the memory for ~15x the speed** at R1.

## Campaign memory ladder = machine classes (Ryan 2026-08-24)

Ryan: *"it might be more insightful to choose memory in terms of what is
typically available: 16G for a laptop, 128 for a workstation, and 500 for a
supercomputer."* `MEM_BUDGETS = 0:16:128:500` (GiB TOTAL memory: body + solver
state + FMM plan + cache).

| budget | meaning |
| --- | --- |
| **0** | cache FORBIDDEN — the matrix-free method endpoint, not a machine. Memory unconstrained. Required: at R1-R5 every real machine can afford a cache, so without it no budget would ever return an uncached winner and the uncached configs would have no tuned knobs. |
| **16** | laptop |
| **128** | workstation |
| **500** | supercomputer node (the pinned `--mem=500G` on 512 GB zen3) |

What a machine class MEANS varies by rung, which is the finding, not a defect
(derived from the W0 map + measured non-cache overhead):

| rung | N | cache floor | saturation | laptop 16G | workstation 128G | node 500G |
| --- | --- | --- | --- | --- | --- | --- |
| R1 | 8,016 | 0.22 G | 1.0 G | saturated | saturated | saturated |
| R2 | 15,760 | 0.55 G | 3.6 G | saturated | saturated | saturated |
| R3 | 28,752 | 1.35 G | 9.3 G | saturated | saturated | saturated |
| R4 | 58,192 | 3.27 G | 27.3 G | **constrained** | saturated | saturated |
| R5 | 108,240 | 8.38 G | 64.4 G | **constrained** | saturated | saturated |
| R6 | 212,108 | 24.2 G | 157.7 G | **UNCACHED** | constrained | saturated |
| R7 | 419,276 | 79.0 G | 427.7 G | **UNCACHED** | constrained | saturated |

Reading: below R4 the whole near field fits on anything, so caching is
unconstrained everywhere and the three machine rows coincide. From R6 a laptop
**cannot cache at all** — the floor exceeds 16 GiB — and is forced matrix-free,
which the tuner now handles automatically by tuning that point uncached.
R6/R7 saturate below 500 GiB, so the full memory axis is reachable on the
pinned hardware without touching the 1 TB (`cs`) or 2 TB (`m11-2`) partitions.

Why not a per-rung geometric ladder: a flat `4:8:...:500` sweep was checked
against the map first and gives only **13 distinct operating points out of 56
descents** (R1-R2 saturate at every budget; R7 is uncached at five of eight).
The machine ladder is 4 points x 7 rungs = 28 descents and its duplicates ARE
the result.

Watch for: `krylov_ilu`'s pattern is ~23 GiB of entries at R7 (see the
`max_pattern_entries` note in `rotor_hover_solver_phase1_agreement.jl`), so
some configs will not fit some machine classes at all. The driver reports that
as the infeasible end of the axis rather than failing — which is itself a
publishable column ("which solvers are even AVAILABLE on this machine").

R1 smoke of the two-point ladder (local, indicative timings only):

| budget | cached | t_warm [s] | cache | mem total |
| --- | --- | --- | --- | --- |
| 0 (cache forbidden) | no | 93.9 | — | 14.7 MB |
| 16 (laptop) | yes | 4.83 | 272 MB | 287 MB |

## CORRECTION: the cache "floor" is softer than first stated (2026-08-24)

The W0 map sampled leaves 9..1100, and I read its smallest-leaf value as a
floor, writing that "shrinking the leaf raises block count without shrinking
total near-field area, so bytes never fall below the small-leaf value". That
was an extrapolation past the sampled range, not a measurement. Checked
directly:

| rung | leaf 9 | leaf 7 | leaf 5 | leaf 4 | leaf 3 | leaf 2 | change 9->2 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| R6 | 23.83 G | 22.91 G | 21.86 G | 21.04 G | 20.19 G | 19.31 G | **-19%** |
| R7 | 78.23 G | 75.40 G | 72.92 G | 71.22 G | 69.46 G | 67.00 G | **-14%** |

Cache bytes DO keep falling below leaf 9, though the curve is clearly
asymptoting toward a genuine MAC-set geometric minimum near ~19 GiB (R6) and
~67 GiB (R7).

**What survives:** a 16 GiB laptop still cannot cache at R6 or R7 — the
cheapest cache is ~19.7 GiB and ~67.8 GiB including overhead. **What changes:**
the R6 margin is 23%, not the ~50% the quoted 24 GiB implied. R7 is unaffected
in substance.

**The real lesson, and a code fix.** Total memory is body + solver state +
plan + cache, and the FMM **plan grows as the leaf shrinks** while the cache
shrinks. The total therefore has an INTERIOR minimum, so pricing a rung's floor
at any single leaf is wrong in both directions — it can overstate the floor by
missing a cheaper small-leaf point and understate it by ignoring plan growth.
`rotor_hover_solver_phase2_tune.jl`'s floor check now SCANS `FLOOR_LEAVES`
(default 2:4:9:21) and takes the minimum of the actual total, reporting both
the cheapest uncached and the cheapest cached configuration.

Note the W0 map reports CACHE bytes only, not total; the per-rung "floor"
column earlier in this ledger should be read as cache-at-leaf-9 plus a scaled
overhead estimate, not as a minimised total.

## CORRECTION: the cache "floor" is softer than first stated (2026-08-24)

The W0 map sampled leaves 9..1100, and I read its smallest-leaf value as a
floor, writing that "shrinking the leaf raises block count without shrinking
total near-field area, so bytes never fall below the small-leaf value". That
was an extrapolation past the sampled range, not a measurement. Checked
directly:

| rung | leaf 9 | leaf 7 | leaf 5 | leaf 4 | leaf 3 | leaf 2 | change 9->2 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| R6 | 23.83 G | 22.91 G | 21.86 G | 21.04 G | 20.19 G | 19.31 G | **-19%** |
| R7 | 78.23 G | 75.40 G | 72.92 G | 71.22 G | 69.46 G | 67.00 G | **-14%** |

Cache bytes DO keep falling below leaf 9, though the curve is clearly
asymptoting toward a genuine MAC-set geometric minimum near ~19 GiB (R6) and
~67 GiB (R7).

**What survives:** a 16 GiB laptop still cannot cache at R6 or R7 — the
cheapest cache is ~19.7 GiB and ~67.8 GiB including overhead. **What changes:**
the R6 margin is 23%, not the ~50% the quoted 24 GiB implied. R7 is unaffected
in substance.

**The real lesson, and a code fix.** Total memory is body + solver state +
plan + cache, and the FMM **plan grows as the leaf shrinks** while the cache
shrinks. The total therefore has an INTERIOR minimum, so pricing a rung's floor
at any single leaf is wrong in both directions — it can overstate the floor by
missing a cheaper small-leaf point and understate it by ignoring plan growth.
`rotor_hover_solver_phase2_tune.jl`'s floor check now SCANS `FLOOR_LEAVES`
(default 2:4:9:21) and takes the minimum of the actual total, reporting both
the cheapest uncached and the cheapest cached configuration.

Note the W0 map reports CACHE bytes only, not total; the per-rung "floor"
column earlier in this ledger should be read as cache-at-leaf-9 plus a scaled
overhead estimate, not as a minimised total.

---

## Phase 16 — candidate-level tuning resume, one CSV per job (2026-08-24)

Three things blocked the Phase 1/2 HPC launch as specified. Ryan ruled on all
three; this records what was implemented and what the results are allowed to
mean.

### 1. Interruption-transparent tuning

**Problem.** Under `physics2/standby` the tuner wrote a budget's row only after
that budget's whole descent finished, and `tune_fmm_perturb` did not memoize
across runs. A preempted R5 job therefore lost 4-8 h of measurement,
repeatedly. `TUNE_MAX_SECONDS` did not fix this — it only converted total loss
into an honestly-labelled *un-converged* answer.

**Ruling (Ryan).** Checkpoint every benchmark during tuning so a re-launched
job picks up where the interrupted one stopped, then just use `--requeue`.

**Implementation.** `tune_fmm_perturb` gained `memo=` (a caller-supplied dict,
seeded from and written into) and `on_measure=` (a hook firing once per
FRESHLY measured candidate). The Phase 2 driver persists every candidate to
`tune_trace_<rung>_b<budget>.csv`, flushed per line, and replays it into the
memo on restart. FastMultipole `4c0f1b8f`.

**One behavioural fix was required, not optional.** `benchmark()` returned
early on a memo hit *before* the `t_best_ok[]` update, so a preloaded memo
would have replayed with the abandonment cutoff stuck at `Inf`, letting
newly-measured candidates escape abandonment. A resumed descent would then have
been both slower than and DIVERGENT from an uninterrupted one. A memo hit now
tightens the cutoff exactly as a fresh measurement does. This is what makes
replay faithful rather than merely fast; the new tests fail 66 assertions
without it.

**Honest limitation.** A replayed `t` was measured on the pre-preemption node.
All jobs are pinned `--constraint=zen3 --exclusive --mem=500G` on
spec-identical hardware, and the winner is independently re-timed by the driver
and again by `rotor_hover_solver_phase2.jl`, so this affects WHICH PATH the
descent walks, never a published number. The trace records its provenance
(rung, budget, reps, abandon factor, hardware tag, fm commit) and REFUSES to
load under a mismatch — a memoized `t` is a timing and must not be replayed
across a change of reps, abandonment rule or hardware.

**`max_seconds` now means something slightly different.** The timeout check
deliberately skips memo hits, so a replay cannot be cut short by it. The
backstop is therefore per-JOB, not per-DESCENT: a repeatedly-preempted R5 can
accumulate more total tuning time than `TUNE_MAX_SECONDS` suggests. That is the
intended consequence of the ruling, but it is a change in what the knob means.

### 2. One CSV per job

**Problem.** Step A's jobs all appended to one `tune_phase2.csv` per rung and
Step B's 7 jobs to one *flat* `agreement.csv`, each holding the file open for
hours over NFS. This is the exact non-atomic append that destroyed R1's rows on
2026-08-18.

**Ruling (Ryan).** Every job writes its own CSV; synthesize afterwards.

**Implementation.** Step A runs one job per rung with all four budgets
sequential (`PER_RUNG_DIR=1` was already honoured), which is exactly one writer
per file and keeps budget 0 first as `phase2.jl` requires. `p2_tune.sh`'s header
had documented one-job-per-BUDGET — four writers on one file, the very hazard —
and was rewritten. `phase1_agreement.jl` received the `PER_RUNG_DIR` branch it
had never been given, copied verbatim from `phase1_case.jl`; its `CONFIGS` split
was widened to `r"[:,]"` so the colon form the wrappers promote no longer
silently selects nothing. New `benchmark/p021_merge.jl` concatenates the
per-rung files and asserts one row per `(rung, budget)` and per `(rung,
config)`, so a missing row is a DETECTABLE absence rather than a silent hole.

### 3. Unrelated fix found while here

`p2_tune.sh` forced `TUNE_REPS=5` unconditionally, overriding the driver's
per-rung default of 2 at R5+ — this would have more than doubled every R5
descent. It now exports `TUNE_REPS` only when explicitly set.

### Verification

Local runs verify CONTROL FLOW ONLY — Ryan's Mac scatters 22-39% at fixed knobs
against a ~15% effect, so no local timing here is evidence.

| test | result |
| --- | --- |
| FastMultipole suite (incl. 240 new memo/on_measure assertions) | PASS, 332,265 assertions, 0 failures |
| Full trace replay, R1 budget 0 | 22 candidates replayed free; identical winner p=17 mac=0.65 leaf=14 |
| Truncated trace (10 of 22) | prefix replayed byte-identical, then continued live |
| Provenance guard (reps / abandon / hardware mismatch) | all three error before any solve |
| Per-job CSV isolation, 2 concurrent rungs | two separate complete CSVs |
| `p021_merge.jl` | flags removed rows, duplicates, header drift; `STRICT=1` errors |
| Cluster dep stack, R1 driver end-to-end on login node | runtime-clean; bc_rel_l2 9.67e-7 certified |

The truncated-trace test also exposes why the *original* acceptance criterion
could not be run locally. Candidate 11, re-measured live, read 96.093 s in the
reference and 93.666 s on resume — pure machine scatter. Comparing a resumed
descent against a freshly re-measured uninterrupted one is therefore NOT a test
of the memo; only replay of STORED timings is deterministic. Path identity at
every cut point is asserted instead on an analytic cost surface in
`test/autotune_cost_test.jl`, which is immune to machine noise.

An unplanned demonstration: the login-node dep check was SIGKILLed by the CPU
cgroup mid-descent, and its trace file survived with the measured candidate
intact — the checkpoint working under exactly the failure mode it exists for.

### State at handoff

Sync complete and audited: FastMultipole fast-forwarded to `4c0f1b8f` on the
cluster (`src/autotune.jl` md5-identical both sides); `benchmark/` 54/54 files
md5-identical, no orphans either side. Superseded `phase1/multi/R1/tune.csv`
(job 13447582, the mis-tuned leaf=43 apply-proxy row) moved aside; all 7
`bcache_R*.bin` intact. Disk 227 G / 2.0 T. Nothing from 021 in the queue.

Submissions are written out in `benchmark/slurm/SUBMIT_021_phase12.txt`.
**Every `sbatch` is Ryan's call and none has been submitted.** R6/R7 tuning
objective remains deferred per §7 of the plan.

---

## 2026-08-28 — Phase 1 HARVESTED, R1–R7 complete

`benchmark/p021_merge.jl` rerun on the cluster (script md5-identical both sides,
`efa1489e…`). The pre-existing merged tables were dated 08-25 16:20 and silently
predated R7 by a day; they are now regenerated. **Phase 1 portion merged CLEAN** —
44 agreement rows, 52 spread rows, every expected row present exactly once, no
duplicates, header consistent across all seven rungs. Full tables:
`phase1_agreement_table_R1_R7.md`.

`CONFIGS` deliberately left UNSET for the merge: the phase-1 expectation is then
"each config that appeared, appeared once", which accepts R7's reduced 4-config
roster without a spurious missing-`krylov_jacobi` error while still catching
duplicates.

| rung | panels | n_configs | max pairwise relL2 |
| --- | --- | --- | --- |
| R1 | 8,016 | 7 (+2 direct) | 1.73e-5 |
| R2 | 15,760 | 7 | 1.87e-5 |
| R3 | 28,752 | 7 | 1.90e-5 |
| R4 | 58,192 | 7 | 2.54e-5 |
| R5 | 108,240 | 5 | 2.93e-5 |
| R6 | 212,108 | 5 | 3.76e-5 |
| R7 | 419,276 | 4 | 4.31e-5 |

**All 44 rows `bc_certified=true`.** Agreement degrades smoothly, 2.5× over a 52×
panel count.

Three results:

- **`krylov_gmres` is the outlier at every rung, and that is the stopping rule, not
  a defect.** Its BC error is pinned at 9.0–9.7e-7 across all seven rungs — it halts
  just under the 1e-6 target — while `fgs` overshoots to 2e-8–2.5e-7 and `krylov_ilu`
  to 3e-7–9e-7. The ensemble spread is therefore largely gmres's own tolerance band,
  so the growth in the table above measures the *target*, not divergence between
  solvers.
- **`krylov_ilu` FIT at R7** — 24 iterations, 3788 s. The anticipated ~23 GiB
  pattern non-fit did not materialise; that contingency (`phase_17` §4 trap 6) is
  closed.
- **`backslash_coupled`'s BC of 4.6e-5 confirmed flat across all four rungs it ran**
  (4.78 / 4.67 / 4.61 / 4.55e-5, slightly *decreasing* with N), consistent with
  `phase_18` §6: same solution, four decades worse BC. Not re-diagnosed here.

The R1 `RESIDUAL_BACKEND=direct` cross-check (job 13469159) validates the FMM BC
metric against O(N²): **0.00 % on `krylov_gmres`, 0.96 % on `fgs`**.

Merge also reported three Phase 2 problems, all expected and none new: R6/R7
`tune_phase2.csv` missing (that tuning was never submitted, per §7 of the plan),
and R3 budget 0 `tune_timed_out=true` (best-so-far, not an optimum).

### OPEN DEFECT — R7 `CT_laplace` collapse

| rung | CT_laplace | CT_bernoulli | CT_kj |
| --- | --- | --- | --- |
| R1 | 0.060406 | 0.053168 | 0.057221 |
| R2 | 0.061307 | 0.053654 | 0.057736 |
| R3 | 0.061021 | 0.053402 | 0.058153 |
| R4 | 0.061241 | 0.053424 | 0.058611 |
| R5 | 0.061430 | 0.053414 | 0.059026 |
| R6 | 0.062278 | 0.053391 | 0.059453 |
| R7 | **−0.001488** | 0.053203 | 0.059882 |

At R7 `CT_laplace` reads −0.001488 against ~0.062 at R1–R6 — a sign flip and a 40×
collapse. **It is not a usable number and must not be published as one.**

What is established:

1. **It is common-mode, not solver disagreement.** All four R7 configs agree on it
   to 3.4e-5 absolute (−0.0014964 / −0.0014956 / −0.0014964 / −0.0014963 for
   `krylov_ilu` / `fgs` / `fgmres_fgs` / `krylov_gmres`). Every solver reproduces the
   same value, so the solve is consistent.
2. **The wake and shedding are NOT implicated.** The `RigidWakeBody`
   shedding-from-raw-mesh-cells invariant was the first hypothesis — its signature is
   exactly a CT collapse with no error raised — and it is **ruled out**: `CT_kj`
   (0.059882) and `CT_bernoulli` (0.053203) both sit exactly on their R1→R6 trends at
   R7. A shedding failure would drag all three force methods down together, since KJ
   is computed from the bound circulation directly.
3. **BC error is unaffected** — all four R7 rows certify at 2.2e-7…9.3e-7. The
   Dirichlet surface condition is satisfied regardless, which is consistent with a
   defect in post-processing rather than in the solved strengths.

This localises the defect to **`PressureLaplace` recovery on the R7 mesh**
(419,276 panels, `dji9443_20260813_..._capped_captess4`). Cause not investigated.
See `agent_policies/MONITORS.md` for the PressureLaplace constraints — the natural
first check is whether one of them is violated at R7 resolution but not at R6.

Note this does **not** taint the Phase 1 agreement deliverable, which is a statement
about solver-to-solver consistency and stands on its own at R7.

### R7 provenance — must be restated wherever the R7 table is published

`krylov_gmres` came from job `13469158` (cancelled mid-run on 2026-08-25 after
`krylov_jacobi` spent 9 h 20 m with no result); `krylov_ilu`, `fgs`, `fgmres_fgs`
came from `13482349`. Same mesh, same commit `c665634…-dirty`, same knobs
`tuned(17, 0.5, 21)`. Legitimate under the reference-free design, which has no
privileged config and no config ordering — but it is split provenance and is stated
as such. `krylov_jacobi` is dropped at **R7 only**; it remains selectable R1–R6.
