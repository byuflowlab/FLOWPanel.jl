# 018 Decision Rules and Metrics

## M1 — settled-mean thrust coefficient

$\bar{C_T}$ = cycle-mean of `CT_bernoulli` over the last **≥15 settled
revolutions** (post-withdrawal). Settledness: successive 5-rev block means of
the averaging window drift < ~0.3% of $\bar{C_T}$ (no monotone trend); extend
the run by restart segments until this holds. Report a moving-block bootstrap
95% CI on per-rev means as supporting evidence. Executor discretion on the
exact test is allowed (Ryan 2026-07-30) — **the requirement is that the mean
settles**; record whatever test was used in the phase file.

**Limit-cycle defense (publication text, and why the old gate is retired):**
the self-sustaining hover wake exhibits bounded limit cycles with periods up
to ~9–9.5 revs (BRAINSTORM/006). A 5-rev window with 0.5% rel-std / 2% p-p
tolerances cannot pass by construction on such a signal and no configuration
has ever passed it. Cycle-means over ≥1 full period are the defensible
observable; quote them with their CI, never a 2–5-rev window.

## M2 — circulation distribution (co-equal observable, ruling 9)

$\bar\Gamma(r/R)$ = per-section cycle-mean of `circulation_te` — the
**trailing-edge μ jump** (upper−lower wake strength), the physically
meaningful bound-circulation estimator — with `circulation_slice` as
cross-check, from the BoundCirculationMonitor CSV, averaged over the same
window as M1, blades averaged (axisymmetry confirmed in 006). Comparison
metric between two runs: $\varepsilon_\Gamma$ = max and RMS of
$|\Delta\bar\Gamma|/\max_r|\bar\Gamma|$ over $0.3 \le r/R \le 0.95$; when the
runs' section grids differ (mesh phases 10/11), interpolate to the common
grid first. **Every convergence axis must pass BOTH M1 and M2** — a rung that
converges CT̄ but not $\bar\Gamma(r/R)$ is not converged.

## Parameter-convergence thresholds (per axis, successive rungs)

- $|\Delta\bar{C_T}| \le 0.5\%$ AND $\varepsilon_{\Gamma,\max} \le 1\%$
  ⇒ converged at the coarser rung.
- Das axis: "smallest that still converges" ≡ lower edge of the log plateau —
  successive-doubling deltas ≤ 0.5% and consistent with the ~0.2%/doubling
  wing law (014). Below-plateau (amplified) points are by definition
  non-converged.
- σ axis: fit $\bar{C_T}(\sigma) = C_{T\infty} + A\sigma^p$ over the 3 rungs;
  accept if $p > 0$ and $|\bar{C_T}(L1) - C_{T\infty}| \le 0.75\%$; the gap
  enters the error budget. Non-monotone or unstable ⇒ contingency chain.
- Null demonstrations: truncation $\le 0.3\%$; merging $\le 0.25\%$
  ($\varepsilon_\Gamma \le 0.5\%$); N spot-check $\le 0.5\%$.
- Green: $|\Delta\bar{C_T}| \le 1\%$ ⇒ velocity stays production, Δ is a
  budget term; else Green becomes production and Phase 9 reruns green.

## Reporting rules

- Every quoted number carries: case tag, job id, averaging window (revs),
  and CI. Deltas are computed between matched windows.
- Convergence and agreement with experiment are separate claims;
  $C_{T,\mathrm{exp}}$ enters only in Phase 9 after the budget closes.
- Reused recorded results (legacy-σ η ladder, naive NT=72, 4R/6R null) are
  labeled *corroborating* — they used `overlap_pps` shedding and are not
  ladder points of this campaign.
