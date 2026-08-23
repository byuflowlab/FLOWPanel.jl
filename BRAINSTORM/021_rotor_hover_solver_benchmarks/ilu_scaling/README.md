# ILU pattern/cost scaling probe (2026-08-22)

`ilu_ladder_20260822.csv` — Barba direct-list ILU pattern size, construction
time, and apply (`ldiv!`) time across the frozen 021 mesh ladder, at the ILU
knobs every benchmark script uses (leaf_size=10, MAC=1.0).

Method: geometry-only probe — the body is built from `benchmark/phase1_case.jl`
truncated before its O(N^2) frozen-`b` assembly (which is not needed for the
pattern, the factors, or the apply), so R6/R7 pattern sizes are reachable on a
16 GB machine. Apply is min-of-k `ldiv!` after a warmup.

`measured=pattern_only` rows have no factors built: at R6/R7 the pattern alone
is 6.5 / 23.1 GiB before values, beyond this machine. Their apply times are
projected in `ledger.md` at the measured ~0.65 ns/nnz, not measured.

Finding and consequences: see the "Phase-1 finding: the Barba ILU
preconditioner scales ~N^1.5, not ~N" section of `../ledger.md`.
