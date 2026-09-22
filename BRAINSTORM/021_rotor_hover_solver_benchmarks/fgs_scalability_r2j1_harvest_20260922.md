# 021: R2-j1 harvest (13829232_7) — owed-item resolution (2026-09-22)

Closes the owed R2-j1 item from `fgs_scalability_stage1_reset_prompt_20260922b.md`
(carried since `fgs_scalability_reset_prompt_20260922.md` §"R1–R2 resume
merge"). Job 13829232_7 finished 2026-09-22 (queue drained ~18 h into its
48 h wall; `phase2/phase2.csv` mtime Sep 22 06:44).

Source: `/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/r12-champion-R2-j1-13778533/phase2/phase2.csv`
(21 data rows; `tune_phase2.csv` 6 rows). Harvested via `harvester`.

## R2, j=1 — per-config median t_solve_min (s), fastest first

| config | t_solve (s) | rows | bc_rel_l2 |
|---|---|---|---|
| backslash_ldiv | 0.064 | 1 | 5.96e-09 |
| fgs | 2.634 | 1 | 3.29e-07 |
| fgmres_fgs_nfcache | 3.684 | 5 | 8.44e-07 |
| krylov_ilu_nfcache | 15.889 | 5 | 3.32e-07 |
| fgmres_fgs | 25.894 | 1 | 8.44e-07 |
| krylov_gmres_nfcache | 102.567 | 5 | 9.53e-07 |
| krylov_ilu | 241.029 | 1 | 3.32e-07 |
| krylov_jacobi | 333.092 | 1 | 8.78e-07 |
| krylov_gmres | 1545.106 | 1 | 9.53e-07 |

Checks: 0 rows excluded; no duplicate config+rep rows; all 21 rows
bc_rel_l2 ≤ 1e-6 measured (range 5.96e-09 – 9.53e-07), independent of the
`bc_certified` self-certification flag. BLAS convention: R1–R2 runs use
BLAS = j (here 1), unlike R4 thread-scaling (BLAS 1 at all j).

## Read

- Completes the R1–R2 resume-merge matrix: with R1-j32 (resolved in
  `fgs_scalability_stage0_audit_20260922.md` §5), all 13829232 arms are
  now harvested. Both formerly-missing arms were case-sensitivity/glob
  misses or still-running, not failures.
- At R2 j=1 the story matches the small-rung pattern: `backslash_ldiv`
  untouchable (0.064 s), FGS the best iterative (2.63 s), nfcache variants
  pay their cache economics with no thread parallelism to amortize them
  (krylov_ilu_nfcache 15.9 s — it only overtakes FGS from j≥8).
- No new owed items from this run.
