# 033 B-I3 — ILU setup + prime attribution at R4 j64 (2026-10-01)

Attribution of the warm-start-R4 ILU arm's one-time costs (t_setup ≈ 80 s,
t_prime ≈ 28 s at R4 j64) into (phase, seconds, threaded?), per the B-I3
staging (Ryan 2026-09-29). Out-of-scope ruling honored: `ILUZero.ilu0` is
NOT parallelized — its cost is recorded as ILU's accepted residual.

## Evidence

1. **Campaign harvest**: the 021 warm-start R4 campaign CSV persisted the
   ILU ctor's stats dict per restart window —
   `/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/fgs-wsr4-13890195/unsteady.csv`
   (job 13890195, R4 = 58,192 panels, j64, CONFIG=krylov_ilu_nfcache,
   `ilu leaf10/mac1.0`, nfcache_max_gib=500, apply knobs P12/MAC0.55/leaf48),
   8 restart windows, <4% spread: t_setup 80.63 (79.00–82.54), t_prime
   27.97 (27.63–28.24), tree 1.89, lists 0.28, assembly 12.10,
   factorization 44.22, t_ctor 0.13.
2. **Fresh instrumented run (gap-closers)**: job `fp033-bi3-iluattr-r4j64`
   **13948848** (campaign `p033-bi2-20261001`, pins in
   `~/campaigns/p033-bi2-20261001/pins.toml`; same knobs/env as the campaign
   arm, N_STEPS=2). Two instrumentation-only additions (FLOWPanel `6e15e6a`):
   `pattern_time` bracket in `_ilu_direct_pattern` (sizing +
   diagonal-coverage pass) and `prime_nfcache_build` persisted from the
   nfcache's own `build_time`. R1 smokes + full unit-solver tests pass
   (513/513).

## Setup split — R4 j64 (fresh run 13948848; campaign numbers in parens)

| phase | what it is | t [s] | threaded? |
|---|---|---:|---|
| t_setup | make_config_solver total | 79.11 (80.63) | — |
| └ t_precond = ILUPreconditioner ctor | | 75.05 (76.44) | — |
|   ├ tree | 2× FastMultipole.Tree (Barba) | 1.92 (1.89) | partially |
|   ├ lists | build_interaction_lists | 0.30 (0.28) | no |
|   ├ pattern | sizing + diagonal-coverage pass, O(entries) | **8.40** (n/a) | no |
|   ├ assembly | kernel eval + `sparse()` CSC | 11.80 (12.10) | kernel loop yes; `sparse()` no |
|   ├ factorization | `ILUZero.ilu0` + pivot check | **43.15** (44.22) | **no — accepted residual (ruling)** |
|   └ other ≈ 9.5 | geometry refreshes, row-nnz scan, pivot/stats incl. `Base.summarysize` | ~9.5 | no |
| └ t_ctor (KrylovSolver ctor) | workspace alloc | 0.14 (0.13) | no |
| └ outside the two @elapsed | env/knob parsing, GC attribution | ~3.9 | — |

## Prime split — R4 j64 (t_prime 27.91 fresh; 27.97 campaign)

| phase | t [s] | threaded? |
|---|---:|---|
| near-field cache build (`prime_nfcache_build`) | **14.60** | yes (`@spawn` chunk queue) |
| plan build + priming GMRES iterations | 13.31 | plan partially; GMRES applies via threaded FMM backend |

## Remaining-headroom statement (feeds Z1 fairness note)

Of ILU's ~107 s one-time cost at R4 j64 (79.1 setup + 27.9 prime):

- **43.2 s (40%) is `ilu0` factorization** — serial by the accepted-residual
  ruling; not addressable within 033's scope.
- **~18 s is serial pattern/bookkeeping** (8.4 pattern pass + ~9.5
  stats/refresh residue) — in principle threadable, claimed by no one.
- The threaded shares (assembly kernels, nfcache build) are already
  parallel; prime's other half (13.3 s) is plan build + an unavoidable
  priming solve.

Fairness note for Tier-2: ILU's effective one-time cost floor under the
ruling is ≈ 43 s factorization + ~28 s prime ≈ 71 s even if every serial
bookkeeping second were threaded away. FGS post-B-I2 full ctor at R4 j64 is
22.0 s (certified, job 13948843) — FGS setup is now ≈ 3.6× cheaper than
ILU's measured 79 s setup (and cheaper than ILU's theoretical floor).
