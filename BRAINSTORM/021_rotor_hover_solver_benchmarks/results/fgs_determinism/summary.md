# R1 FGS determinism attribution and fix

Run date: 2026-08-17. Fixture: R1, expansion order 6, MAC 0.3,
leaf size 150, five inner iterations, ten fixed outer iterations. Each matrix
coordinate used five cold solves with one solver, three fresh constructors,
ten representative-leaf factor/solves, one warmup, and five timed cold solves.
The driver set `GOTO_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, and
`OMP_NUM_THREADS` and asserted the effective BLAS count after LAPACK.

## Attribution

| Component | Pre-fix 4T/BLAS-1T | First difference | Verdict |
|---|---:|---:|---|
| Fresh constructors, tree topology/geometry | exact | none | not causal |
| M2L/direct lists and packed matrices | exact | none | not ctor variation |
| Upward pass | exact | none | not causal |
| Fixed-input M2L | differs in all repeats | packed expansion index 2353 or 3361 | root cause |
| Fixed-input L2L and L2B | exact | none | not causal |
| Representative leaf LU/solve | exact | none | not causal |
| Iterate strengths | differs in all repeats | outer state 2, first strength as early as 1 | downstream effect |
| Residual sequence | differs in all repeats | outer iteration 3--5 | downstream effect |

The raw interaction list was deterministic but not owner-major. Threaded
`assign_m2l!` assumes all interactions for a target are contiguous; without
that premise, separate assignments concurrently accumulated into the same
target expansion. M2M/L2L ownership remains a verified static invariant:
their assignments contain complete parent branches.

## Fix

Immediately after list construction, before matrices or index maps, apply
FastMultipole's stable counting sort by source and then by target. Stability
makes the result lexicographic `(target, source)`, gives every threaded M2L
target one owner, and determinizes within-owner order. The same canonical
ordering is applied to the direct list before dependent matrix/index maps.

On the actual 1,760-entry R1 M2L list, the selected FastMultipole two-pass
counting sort took 8.208 us minimum versus 25.791 us for generic
lexicographic `sort!` (3.14x faster, identical output).

## Post-fix verdict and timing

| Julia threads | BLAS threads | Exact repeated solves/ctors/LU | Fixed-10 minimum (s) |
|---:|---:|---:|---:|
| 1 | 1 | PASS | 1.000965 |
| 1 | 4 | PASS | 1.122893 |
| 4 | 1 | PASS | 0.894550 |
| 4 | 4 | PASS | 0.999269 |

The selected Julia-4T/BLAS-1T pre-fix replay was nondeterministic and took
0.886844 s minimum; the post-fix minimum is 0.894550 s, a 0.87% slowdown.
This passes the 5% performance gate. Equality is required within, not across,
fixed thread configurations.

## Margin verification

The stopping tolerance is the geometric mean of the BC-crossing internal
residual and the first later finite residual that is strictly lower. R3 at
tau 1e-6 initially lacked such a successor; its staircase was extended rather
than using a fallback. All 18 R1--R3 selections pass certified BC error <= tau
at Julia 4T/BLAS 1T. Worst achieved error/tau ratios were 0.482 (R1), 0.597
(R2), and 0.760 (R3).

## Provenance and verification

- FLOWPanel: `b2510715b6e299bddbd900b12b878ec6badcfa11-dirty`.
- FastMultipole baseline: `5adde3bda9b62458c349ddbda0bfa77c7e14351d-dirty`;
  the determinism hunk is limited to canonical list ordering plus its
  premise-guarded regression. Existing campaign fixes remain uncommitted and
  `scripts/20250404_prediction_accuracy.jl` was untouched.
- FastMultipole `solve_test.jl`, Julia 4T/BLAS 1T: 582,335 checks passed.
- FLOWPanel solver history: 15/15 passed. Solver unit suite: 200/200 passed.
- No LU caching, colored Gauss--Seidel, or batched nonself update was attempted.

Raw files: `probe_pre_j4_b1.csv`, `probe_j{1,4}_b{1,4}.csv`,
`timing_*.csv`, and `ordering_*.csv` in this directory. Margin evidence:
`benchmark/results/phase1/multi/fgstune_margin_verify.csv` and its three
per-rung sidecars.
