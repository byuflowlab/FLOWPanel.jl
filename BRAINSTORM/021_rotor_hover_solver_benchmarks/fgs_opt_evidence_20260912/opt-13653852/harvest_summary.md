# Job 13653852 harvest

The output root completed (`COMPLETED`; Slurm 13653852 and all recorded
stages exit 0). Inner screen mapping is read from each candidate's
`requested_config.toml`/`config.toml`, not directory creation order. All five
screen candidates were completed and eligible; no screen arm failed.

The ranked prepared-only medians are in `screen_rank.csv`. The two fastest
accepted candidates are inner=3 (1.39697637 s) and inner=5 (1.44921144 s),
with the frozen seed inner=10 at 1.52587803 s. Relative to the seed, inner=3
is 8.45% faster and inner=5 is 5.02% faster by median prepared solve time.
Each candidate has five repetitions and zero recorded fresh-prepared relative
delta; all calibration and final validation gates pass.

Gate values for the screen are bounded by: authoritative certified-FMM
relative L2 <= 7.78557547e-7; evaluator disagreement <= 5.29506755e-9; and
repeat delta = 0, against limits 1e-6, 1e-7, and 1e-8 respectively. The
reported `fmm_rel_max` is retained in each candidate's
`convergence_validation.csv`; it is not the evaluator-disagreement gate.

The seed baselines are under `baseline-j4-b1/.../j4_b1` and
`baseline-j64-b1/.../j64_b1`; their prepared medians are 1.82464498 s and
1.51121875 s respectively, both completed and eligible.

Profile evidence is retained under
`profile-fgs-j64-b1/R2/fgs_4ea88e39ddb33b21/j64_b1/`: `cpu_flat.txt`,
`cpu_tree.txt`, and `allocations.txt` (no `.jls` files). The profile banner
records Julia 1.11.7, 64 Julia threads, BLAS=1, m12-3-5, and the pinned
FLOWPanel/FastMultipole/FLOWVPM commits. The 3,857-snapshot flat profile
contains 2,715 samples at `LinearAlgebra.gemv!` and 2,697 at
`compute_nonself_products!`, versus 296 at `solve_leaf!`; it also shows
`dgemv_64_` with 2,712 samples. This preserves the R2 attribution that
nonself dense GEMV remains dominant. The profile log reports warnings for
groups with no samples, but the main thread profile is full (100% utilization).
These are overlapping stack counts, not additive wall-time shares. This R2
screen is preliminary five-trial evidence, not a final cross-rung speedup claim.
