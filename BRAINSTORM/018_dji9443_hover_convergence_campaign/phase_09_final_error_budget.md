# Phase 9 — Final Production Run + Error Budget

**Objective:** the study's headline numbers: CT̄ ± CI and Γ̄(r/R) ± CI at the
converged operating point (Das\*, NT\*, σ\*, N=4, 4R, production formulation
from Phase 7), over ≥30 settled revolutions, with the error budget closed.

## Cases

| tag | knobs | schedule |
| --- | --- | --- |
| `p018_final` | converged settings via env at submission | extend the existing final-settings run (Phase 6/7/8 twin) by restart segments to ≥30 settled revs; only run fresh if no twin exists at the exact final settings |

## Deliverables

1. CT̄ ± bootstrap CI (M1) and Γ̄(r/R) ± CI (M2) over the full settled window;
   settling plot (CT̄ vs window start) as in Phase 1.
2. `error_budget.md` completed — every term filled from its phase, disclosures
   included.
3. Convergence statement per axis (one quantitative sentence each).
4. ONLY THEN: comparison against CT_exp = 0.072 and the CCBlade BEM bracket
   (0.060–0.071), with the explicit statement that no setting was selected
   using the target (ruling 3).
5. If evidence gathered during the campaign suggests the relaxation filter
   would help (e.g. residual ripple limits the CI), RECOMMEND it for follow-up
   — do not run it (ruling 2).
6. Item file: phase-gate table updated to COMPLETE; INDEX.md Outcome cell
   updated (concise); decision log entry.

## Exit criteria

All deliverables written; item ready for clear-context review.

## Log

(append dated entries here)

## 2026-08-13 — rlxf sensitivity ladder (budget term 9): CLOSED by double ignition

The 3-point ladder {0.3 = `_s3`, 0.15 = 13157881, 0.075 = 13157882}
(production-carrier clones, warm from `_s2`@1619 = rev 45, banners +
warmstart verified) produced NO scorable reduced-rlxf window: **both
reduced rungs ignited from an identical healthy warm state.**

| rung | job | onset (max_u>100 / min_sr<0.06) | onset rev | +revs from handoff | peak max_u | end state |
| --- | --- | --- | --- | --- | --- | --- |
| rlxf 0.15 | 13157881 | step 1734 / 1735 | 48.17–48.19 | +3.2 | 47,492 m/s @1737 | min_sigma pinned at 9.407e-5 viscous floor, min_sr 0.020; scancel'd @1772 |
| rlxf 0.075 | 13157882 | step 1684 / 1695 | 46.78–47.08 | +2.0 | 1.12e6 m/s @1706 | wake self-annihilated 217k→11k particles; ran to 2087, exit 0, garbage forces |
| rlxf 0.3 (`_s3`) | 13157752 | no ignition | — | — | 49.8 m/s | healthy (min_sr ≥0.117) through step 1829+ |

Both onsets are abrupt (min_sr 0.125→<0.06 within 1–4 steps) with
min_sigma pinned at the √(2νΔt)=9.41e-5 m viscous floor — the 019
Γ-side-ignition signature. Dose-response: quartering rlxf ignites ~1.2
revs sooner than halving. **Verdict: relaxation at stock 0.3 is
load-bearing for stability at σ=0.04R/N=1 viscous; the term-9 downward
slope is unmeasurable at this operating point.** The ignition boundary
itself is the ladder's result (bounds rlxf ≥ ~0.3 on this carrier).
Carrying mode is a RYAN DECISION (options: reframe as model definition;
probe rlxf>0.3; joint σ–rlxf probe at a more stable operating point). No
further submissions without Ryan. Evidence: local
`data/p018_rlxf0p{15,075}_n1_s040_visc/monitors/` (full wake_health +
force + inventory CSVs); cluster VTK swept to newest-10 restart sets
(6325+5782 MB freed, retention rules; CSVs/TOML/monitors sacrosanct).

**VTK deletion log 2026-08-13 ~19:30 MDT (rlxf upward pair):** live sweep
— rlxf0p45 15,616MB (kept 1966–1975), rlxf0p6 15,595MB (kept 1961–1970);
newest-10 + monitors/CSVs retained. Health: both arms rev ~55, no
ignition ~10 revs in.

## 2026-08-14 — rlxf UPWARD pair harvested (13159912/13): stable, slope UNRESOLVED

Both COMPLETED ~19h, healthy tails (exit 0 trusted): rlxf 0.45 min_sr
0.153 / max_u 29.6 / n_p 221.2k; rlxf 0.6 min_sr 0.201 / max_u 32.8 /
n_p 231.0k. Both reached step 2087 = rev 58, exactly matching `_s3`.
**Stability dose-response confirmed upward: tail min_sr rises with rlxf
(0.3: 0.149 → 0.45: 0.153 → 0.6: 0.201), and neither upward rung shows
any ignition precursor over +13 revs — against the reduced rungs'
ignition at +2.0/+3.2 revs.**

Matched revs 45–57 (identical windows, same warm seed `_s2`@1619):

| rlxf | run | CT̄ | 95% CI | per-rev std | ΔCT vs 0.3 | M2 vs `_s3` |
| --- | --- | --- | --- | --- | --- | --- |
| 0.3 | `_s3` | 0.076806 | [0.075467, 0.078458] | 0.002606 | — | — |
| 0.45 | rlxf0p45 | 0.075475 | [0.074381, 0.076574] | 0.001979 | −1.73% | ε_max 0.980% PASS |
| 0.6 | rlxf0p6 | 0.075884 | [0.074293, 0.076978] | 0.002038 | −1.20% | ε_max 0.708% PASS |

**Verdict: the term-9 local slope is UNRESOLVED at 13-rev windows** —
deltas are NON-MONOTONE in rlxf (0.45 sits below both neighbors), all
three CIs overlap heavily, and quiet-limit regressions over the same
window (0.3: 0.074295, 0.45: 0.072214, 0.6: 0.073603; n_quiet=3 each)
are equally non-monotone. The observed |Δ| ~1.2–1.7% is the same order
as the burst-sampling scatter established in phase_15 (episode means
swing ±3–8%), so CT differences here cannot be attributed to rlxf.
Γ̄(r/R) is genuinely insensitive (M2 PASS ≤1% both rungs). Compare 006
legacy downward slope ~+0.0011 CT per halving (~1.5%): same order as
the noise floor of this window size ⇒ a measured local sensitivity
would need either much longer settled windows or quiet-epoch-only
scoring. **Sum for Ryan: rlxf ∈ [0.3, 0.6] moves CT̄ by less than the
burst noise and Γ̄ by <1%, improves stability monotonically, and rlxf <
0.3 ignites. Options for carrying term 9 unchanged from phase_09
§2026-08-13 (model-definition reframe now has direct support).**
