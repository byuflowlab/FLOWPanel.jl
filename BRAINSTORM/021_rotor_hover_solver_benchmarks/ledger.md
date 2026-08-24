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
