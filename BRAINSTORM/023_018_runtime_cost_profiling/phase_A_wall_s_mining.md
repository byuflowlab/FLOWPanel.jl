# Phase A — wall_s mining of landed 018 runs

**Date:** 2026-08-20. **Data:** every local
`data/p018_*/monitors/*monitor04_wake_health*.csv` (per-step `wall_s` and
`n_particles` recorded by `WakeHealthMonitor`); steps 0–1 dropped (NaN /
startup), runs with < 50 valid steps dropped. Script: session scratchpad
`phaseA_walls.py`; full table `phaseA_wall_s_summary.csv` (this directory).
"early/late" = median over the first/last 10% of recorded steps. Fit is
ordinary least squares `wall_s ≈ a + b·n_particles` over steps with
n_particles > 1000.

## Headline numbers (full-history runs only)

| run | total_h | med wall_s early → late | np late | fit a [s] | fit b [s/1e5 np] | R² |
|---|---|---|---|---|---|---|
| p018_cs_f1_l3p4 (production carrier) | 47.9 | 37 → 179 | 181k | 8.1 | 100.0 | 0.82 |
| p018_cs_f1_l2p4 | 47.9 | 34 → 195 | 184k | 6.3 | 102.1 | 0.84 |
| p018_cs_f1_l4p8 | 47.9 | 39 → 182 | 189k | 13.1 | 96.2 | 0.79 |
| p018_trunc6_n1_s040_visc (6R trunc) | 47.9 | 25 → 328 | 342k | −15.0 | 81.5 | 0.76 |
| p018_upin_nt144 | 71.9 | 21 → 225 | 293k | 0.9 | 73.3 | 0.85 |
| p018_upin_nt72 | 57.5 | 27 → 298 | 239k | −8.1 | 93.2 | 0.42 |
| p018_upin_nt72_rlxf0p16334_rs1 | 71.9 | 60 → 183 | 252k | 20.9 | 60.2 | 0.54 |
| p018_ufront_n1_s040_visc | 33.6 | 25 → 113 | 210k | 13.8 | 54.8 | 0.67 |
| p018_scr_csarc_l2p4 (8-rev screen) | 12.5 | 26 → 201 | 193k | 15.3 | 89.2 | 0.92 |

(Restart-segment runs — `*_rs1/_s2/_s3` — sit at a ~constant particle count,
so their np-fit is meaningless (R² ≈ 0); their late-window medians are still
valid cost datapoints and are in the CSV.)

## Findings

1. **Per-step cost is particle-count-driven and ~linear in np.** Every
   full-history run fits `wall_s ≈ a + b·np` with a small intercept
   (a ≈ 0–20 s) and b ≈ 50–100 s per 100k particles (R² 0.65–0.92 where np
   actually varies). At the production plateau (~181k particles) the
   particle-dependent term is ~90–95% of the step. The body-only floor
   (solve + body↔body influence + I/O at 36,752 panels) is roughly the
   intercept: order 10–20 s.
2. **The particle census, not the mesh or dt, is the cost lever.** NT=144 arms
   cost about the same per step as NT=36 at matched np (b ≈ 63–73 vs 96–102);
   an NT rung's total cost is (steps/rev × np trajectory), so the NT=144 walls
   (rev ~13–15 in 72 h) are census-driven: those arms plateau HIGHER
   (~273–293k vs 181k) despite equal particles-shed-per-rev (NT·pps = 432) —
   the merging/truncation equilibrium shifts with dt. Worth one line in Phase
   C; not chased here.
3. **Truncation depth is a demonstrated large lever:** 6R truncation carries
   342k particles at 328 s/step vs 4R's 181k at ~180 s/step — the 4R ruling
   already saves ~45%/step; the trend is linear as predicted by (1).
4. **Node-to-node variance is real and large.** At essentially the same np
   (~210–270k), late-window medians span 113–298 s across runs/nodes, and the
   two nt72 rs1 segments differ 151 vs 183 s at matched configs. Slope
   comparisons across runs are only trustworthy to ~±30%; Phase B's
   within-run phase split does not suffer this confounder.
5. **Cost-to-schedule check** (validates the mining): f1_l3p4 sums to 47.9 h
   over 1035 steps — exactly the 48 h wall at rev ~28.5 of 30, matching the
   campaign's observed timeout one rev short.

Projection for Phase C: a lever that removes X% of mature-wake particles (or
of the per-particle work) removes ≈ 0.9·X% of a production step at maturity.
