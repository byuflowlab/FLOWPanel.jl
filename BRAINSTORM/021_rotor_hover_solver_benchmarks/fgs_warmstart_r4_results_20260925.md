# 021 warm-start R4 campaign — results (harvested 2026-09-26)

Campaign per `fgs_warmstart_r4_provenance_20260924.md` (pins: tag
`campaign/p021-fgs-warmstart-20260924`, FLOWPanel `0aaceaf1`, FastMultipole
`745af760`, FLOWVPM `8d4a3b4d`; deploy tree
`/home/rander39/campaigns/p021-fgs-warmstart-20260924/FLOWPanel.jl`).

- **Job 1 = 13890195** (warm-start A/B, m12-3-12, zen3 exclusive):
  COMPLETED 2026-09-26 00:34 UTC, 19:17 elapsed, **14/14 legs "ok"**, all
  655 CSV rows `solved=true`, `bcerr_certified=true` on every row,
  `nsolves=1` everywhere. The winB restart fix (`0aaceaf`) held through all
  seven restarted legs. Harvest audit flags exactly ONE anomalous step:
  `fgs_cold` winB has 1 step with bcerr above the certified tolerance
  (still certified-column true; see harvest summary).
- **Job 2 = 13889502** (cold R4 FGS under the new dagteam+backoff default):
  all 5 rungs COMPLETED + evaluator-certified.
- **Follow-up job 13899021** (Ryan-directed `ilu_nfcache_proj2` arm, tag
  `campaign/p021-fgs-warmstart-20260926`, FLOWPanel `b6a09244`; see the
  2026-09-26 addendum in the provenance file): COMPLETED 2026-09-26,
  2:28 elapsed, both legs "ok" (winA + winB restarted from the existing
  ILU step-108 checkpoint), all rows `solved=true`, zero audit flags.
  Tables below include this arm (harvest re-run 2026-09-26).
- Raw harvest outputs (backing data for every table below):
  `fgs_warmstart_r4_results_20260925/harvest_{summary.md,traces.csv,solution_deltas.csv}`
  (copies of the run-dir originals in
  `<data root>/p021-cold-20260910/fgs-wsr4-13890195/`).

Arms: `fgs_*` = FGSSolver (dagteam+backoff default), `ilu_nfcache_*` =
KrylovSolver ILU + near-field cache. Warm-start modes: `cold` =
zero-initial-guess, `prev` = previous-step strengths, `proj1`/`proj2` =
polynomial projection of order 1/2. Windows (transients INCLUDED, Ryan
2026-09-24): **A** = steps 1–36 from scratch; **B** = fourth revolution,
steps 109–144, restarted from the family checkpoint at step 108.

**Required reporting note (Ryan 2026-09-24):** solver warm-start histories
are not serialized in checkpoints, so each restarted warm leg's first
(order+1) steps are effectively cold — inside Window B's
transient-included stats by design. The refill is directly visible in the
niter traces below.

## Window A (startup, steps 1–36)

Time-to-target = `t_solve` (s), the per-step solve time to the certified
tolerance. Setup costs are NEVER amortized into per-step numbers.

| arm | t_solve mean ± spread (median) [s] | niter_first mean ± spread (med) | t_project med [s] | t_setup / t_prime [s] | unconverged | bcerr>tol |
|---|---|---|---|---|---|---|
| fgs_cold | 14.62 ± 6.7 (14.44) | 28.25 ± 2 (28) | 0 | 315.35 / 0 | 0 | 0 |
| fgs_prev | 12.58 ± 9.3 (12.40) | 19.89 ± 14 (20) | 0.00103 | 317.51 / 0 | 0 | 0 |
| fgs_proj1 | 11.70 ± 10 (11.35) | 16.25 ± 19 (15) | 0.00111 | 317.14 / 0 | 0 | 0 |
| fgs_proj2 | 11.15 ± 11 (10.68) | 13.94 ± 23 (12) | 0.00118 | 319.51 / 0 | 0 | 0 |
| ilu_nfcache_cold | 11.10 ± 2.6 (11.06) | 16.25 ± 1 (16) | 0 | 81.67 / 28.22 | 0 | 0 |
| ilu_nfcache_prev | 10.22 ± 2.8 (10.19) | 11.78 ± 6 (12) | 0 | 79.41 / 28.24 | 0 | 0 |
| ilu_nfcache_proj1 | 9.92 ± 3.6 (9.90) | 10.31 ± 7 (10) | 0.00033 | 81.86 / 28.13 | 0 | 0 |
| ilu_nfcache_proj2 | 9.12 ± 3.0 (9.11) | 9.64 ± 8 (9) | 0.00034 | 82.24 / 27.67 | 0 | 0 |

## Window B (developed wake, steps 109–144, restarted)

| arm | t_solve mean ± spread (median) [s] | niter_first mean ± spread (med) | t_project med [s] | t_setup / t_prime [s] | unconverged | bcerr>tol |
|---|---|---|---|---|---|---|
| fgs_cold | 12.97 ± 5.1 (12.88) | 21.89 ± 1 (22) | 0 | 315.35 / 0 | 0 | **1** |
| fgs_prev | 10.08 ± 8.7 (9.86) | 9.50 ± 14 (9) | 0.00101 | 317.51 / 0 | 0 | 0 |
| fgs_proj1 | 8.85 ± 9.2 (8.50) | 4.19 ± 19 (3) | 0.00117 | 317.14 / 0 | 0 | 0 |
| fgs_proj2 | 8.63 ± 9.7 (8.25) | 3.22 ± 20 (2) | 0.00116 | 319.51 / 0 | 0 | 0 |
| ilu_nfcache_cold | 10.64 ± 1.3 (10.65) | 15.00 ± 0 (15) | 0 | 81.67 / 28.22 | 0 | 0 |
| ilu_nfcache_prev | 9.27 ± 3.3 (9.24) | 7.08 ± 9 (7) | 0 | 79.41 / 28.24 | 0 | 0 |
| ilu_nfcache_proj1 | 8.73 ± 3.4 (8.68) | 5.11 ± 11 (5) | 0.00034 | 81.86 / 28.13 | 0 | 0 |
| ilu_nfcache_proj2 | 8.25 ± 3.1 (8.22) | 4.42 ± 12 (4) | 0.00032 | 82.24 / 27.67 | 0 | 0 |

## Per-step traces (full traces in `harvest_traces.csv`)

`niter_first`, first 12 steps of each window — the Window-B history refill
is exactly the designed (order+1)-step cold transient:

| arm | winA steps 1–12 | winB steps 109–120 |
|---|---|---|
| fgs_cold | 27 28 28 29 29 29 29 29 29 29 29 29 | 21 21 21 21 22 22 22 22 22 22 22 22 |
| fgs_prev | 27 24 25 24 22 21 19 16 13 16 17 18 | 21 **11** 10 10 10 10 10 10 10 10 10 10 |
| fgs_proj1 | 27 28 21 23 21 20 19 18 17 15 14 14 | 21 21 **4** 5 4 3 2 2 3 3 3 2 |
| fgs_proj2 | 27 28 21 24 21 19 18 17 15 14 12 10 | 21 21 **4** 4 2 5 1 2 1 2 2 2 |
| ilu_nfcache_cold | 16 16 16 16 17 17 17 17 17 17 17 17 | 15 15 15 15 15 15 15 15 15 15 15 15 |
| ilu_nfcache_prev | 16 14 14 13 13 12 11 10 10 10 10 11 | 15 **7** 7 7 7 7 7 7 7 7 7 7 |
| ilu_nfcache_proj1 | 16 14 13 13 12 12 11 10 10 9 9 9 | 15 **7** 5 5 5 5 5 5 5 5 4 4 |
| ilu_nfcache_proj2 | 16 14 13 14 13 11 10 10 9 8 8 8 | 15 **7** 5 4 4 4 4 4 3 3 4 4 |

Projection engages at step order+2 as designed; in the developed wake
(winB) FGS proj1/proj2 settle at 2–5 iterations vs 22 cold — an ~8×
iteration reduction — while per-step time only drops ~35% (12.88→8.25 s
median): the FMM/near-field product per iteration dominates less than
setup phases at low iteration counts.

## Cumulative cost and crossover (setup + Σ t_solve+t_project)

Cumulative through Window A (36 steps), setup INCLUDED:

| arm | cumulative @ step 36 [s] |
|---|---|
| fgs_cold | 841.6 |
| fgs_prev | 770.5 |
| fgs_proj1 | 738.4 |
| fgs_proj2 | 720.9 |
| ilu_nfcache_cold | 509.4 |
| ilu_nfcache_prev | 475.4 |
| ilu_nfcache_proj1 | 467.0 |
| ilu_nfcache_proj2 | 438.3 |

- Every ilu_nfcache arm is ahead of every FGS arm from **step 1** (FGS
  setup 315–320 s vs ILU setup+prime 108–110 s) and stays ahead through
  step 144 of this campaign.
- **Steady-state crossover, best-vs-best (Window-B medians, t_solve +
  t_project): NONE.** With the order-matched head-to-head now measured,
  `ilu_nfcache_proj2` (8.216 s/step) is *faster* than `fgs_proj2`
  (8.246 s/step) by 0.031 s/step median (0.384 s/step by means, 8.250 vs
  8.634) *and* carries the ~210 s cheaper setup — fgs_proj2 never breaks
  even. The per-step median gap is within noise (~0.4%, spreads ±3.1 vs
  ±9.7), so the honest statement is: per-step it is a statistical tie,
  and the setup difference then decides every horizon in ILU's favor.
- The previous **~480-step crossover** was an artifact of the original
  slate's order asymmetry (fgs_proj2 vs ilu_nfcache_*proj1*, 8.68 s/step);
  it is superseded by the row above. With `prev` warm-start the sign was
  already reversed (fgs_prev 0.62 s/step slower) with no crossover.
- Note the ILU numbers here exclude nothing: its near-field cache build is
  inside the reported `t_setup`/`t_prime` columns.

## Cross-arm solution agreement (rel-L2 strength deltas vs `fgs_cold`)

From `*_strength_snapshots.bin` (per-step full strength vectors; table in
`harvest_summary.md`):

| arm | winA mean / max | winB mean / max |
|---|---|---|
| fgs_prev | 3.6e-05 / 2.0e-04 | 1.0e-05 / 1.5e-05 |
| fgs_proj1 | 8.4e-05 / 2.4e-04 | 4.4e-05 / 6.3e-05 |
| fgs_proj2 | 4.6e-05 / 1.4e-04 | 3.1e-05 / 6.9e-05 |
| ilu_nfcache_cold | 6.9e-05 / 2.4e-04 | 6.0e-04 / 1.1e-03 |
| ilu_nfcache_prev | 7.4e-05 / 4.0e-04 | 6.0e-04 / 1.1e-03 |
| ilu_nfcache_proj1 | 6.4e-05 / 1.9e-04 | 6.0e-04 / 1.1e-03 |
| ilu_nfcache_proj2 | 6.7e-05 / 3.8e-04 | 6.0e-04 / 1.1e-03 |

Warm-starting does not move the answer: within a solver family every
warm mode agrees with its cold reference to ≤4e-04. The larger (~6e-04 to
1.1e-03) FGS-vs-ILU Window-B deltas carry BOTH the two family checkpoints'
divergence and the known ~2e-3 FGS-vs-Krylov wake-on fixed-point
discrepancy (`rigid_motion_tree_reuse_item.md` §5) — REPORTED here, not
chased, per the campaign charter.

## CT traces

- Window A: all eight arms (incl. ilu_nfcache_proj2, CT@36 = 0.0509585)
  agree to 1e-5 absolute — CT@36 spans 0.050952–0.050959 (spread 0.01%). Checkpoint legs continue to step 108
  (CT@108 = 0.014940 fgs / 0.014939 ilu, agreeing to 5e-6 — the low value
  is the usual mid-transient CT dip of this fixture, not a solver effect).
- **Caveat: every restarted (winB) row logs CT = −0** — the restart path
  does not reattach the CT reduction (all other columns, including bcerr
  certification and strength checksums/snapshots, are live). Window-B
  physics agreement is therefore certified via the strength-snapshot
  deltas above rather than CT. Worth a one-line driver fix if winB CT is
  ever needed directly.
- Full CT traces (non-restarted legs) are in `harvest_traces.csv`.

## Job 2 — the 2026-09-22 R4 cold table, re-issued

New column measured under identical contract (cold isolated solves,
`COLD_PREPARED_ONLY=1`, champion placement, BLAS=1, certified evaluator;
run dirs `fgs-cold-newdefault-j{1,8,16,32,64}-13889502`). Values are
median accepted solve time in seconds (Job 2 minima in parentheses).

| j | FGS-dagteam (13777133) | **FGS-dagteam+backoff (new default)** | FGS-colored | krylov_ilu (budget-0) | krylov_ilu_nfcache (budget-500) |
|---|---|---|---|---|---|
| 1 | 33.07 | **32.79** (32.74) | 35.58 | FAILED | FAILED |
| 8 | 5.97 | **5.88** (5.85) | 9.07 | 325.7 | 4.25 |
| 16 | 4.50 | **4.27** (4.25) | 7.26 | 158.5 | 3.56 |
| 32 | 4.42 | **3.53** (3.52) | 6.91 | 83.8 | 2.73 |
| 64 | 6.42 | **3.27** (3.24) | 12.30 | 40.3 | 2.41 |

The j64 min of 3.24 s lands exactly in the 3.2–3.3 s projection from the
Stage-1 cap16 A/B (3.41 s median there). Backoff removes the 64-thread
regression outright (6.42 → 3.27 s, ratio vs j32 now 0.93× instead of
1.45×) and the 16→32 plateau (1.03× → 1.21× gain). krylov_ilu_nfcache
remains the fastest cold R4 solve at every j ≥ 8 (same caveats as the
original table: ~8.5 GB cache + ~9.8 GB solver state vs the far lighter
dagteam).

## Compete/no-compete data (decision is Ryan's)

Presented per the campaign charter — medians, spreads, and crossover;
no pre-judgment:

- **Per-step, developed wake (winB medians):** with the order-matched arm
  landed, ilu_nfcache_proj2 8.22 s and fgs_proj2 8.25 s are a statistical
  tie (means favor ILU, 8.25 vs 8.63); the FGS spread (±9.7) is ~3× the
  ILU spread (±3.1).
- **Cumulative:** ILU's ~210 s setup advantage plus per-step parity means
  fgs_proj2 has NO break-even horizon against ilu_nfcache_proj2 under
  these measurements — the previous ~480-step crossover was vs the
  order-1 ILU arm and is superseded.
- **Iteration budget:** warm-start compresses FGS niter_first far more than
  ILU's (22→2 vs 15→4), so any future per-iteration cost reduction (e.g.
  the Stage-2 sweep-width work) leverages FGS disproportionately.
- **Memory:** ilu_nfcache carries ~18 GB of cache+state at R4; FGS is far
  lighter (unchanged from the 2026-09-22 caveats).
- Solution quality is equivalent across all arms including the follow-up
  (tables above); the only audit flag across all rows is one bcerr>tol
  step in fgs_cold winB.

## Status / gates

- Harvest COMPLETE for both jobs + the 13899021 follow-up (this file is
  the deliverable of record; backing copies refreshed 2026-09-26 after the
  follow-up harvest).
- Storage: run dirs + ckpt trees (`fgs_wsr4_R4_ckpt_{fgs,ilu}`, 847 M each;
  checkpoints confirmed to live INSIDE those run dirs as nested .vtm/.vtu —
  the "0 .vtp" monitor oddity was a wrong-depth glob) handed to the
  archive-first flow 2026-09-26.
- Notebook entry OFFERED to Ryan (FGS Stages 1+2, gate-0, dagedge, default
  adoption, warm-start campaign + winB restart fix) — verbosity per topic
  pending his call.
- Origin pushes of ALL campaign tags (now incl.
  `campaign/p021-fgs-warmstart-20260926`) still owed after `gh auth login`.
- Ryan may still veto the four submission defaults used (champion cold
  tolerance 3.4309419310610173e-7; 36 h walltime; ilu_nfcache_proj1 kept;
  sequential legs on one exclusive node).
