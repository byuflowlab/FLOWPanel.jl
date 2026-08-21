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
