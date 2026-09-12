# R1 CPU pilot evidence — 2026-09-12 UTC

These CSVs were harvested without running Julia locally. Raw authoritative results remain under `/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/pilot-13653115/` on ORC; the complete text snapshot is `/tmp/p021-cold-13653115-evidence/`.

`13653115_timing_trials.csv` contains all 60 measured samples; `timing_summaries` contains min/median/max/spread. `all_compile_validation`, `all_warmup`, `all_fresh_warmup`, and `all_convergence` retain excluded first-call/warmup and independent convergence evidence. `validation`, `status`, `provenance`, and `config_hashes` summarize original files; original TOMLs remain authoritative.

NaN direct-evaluator fields in timed trials mean the independent direct crosscheck was outside that scope. Use warmup/convergence validation for direct comparisons. `fmm_rel_max` is a residual infinity norm, not FMM/direct disagreement (`evaluator_delta`).

Original job 13653115 failed only at profile startup; profiles are being continued in a separate immutable v7 generation, job 13653450. No original output was altered.
