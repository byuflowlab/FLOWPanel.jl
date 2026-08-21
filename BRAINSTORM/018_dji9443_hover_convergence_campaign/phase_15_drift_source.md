# Phase 15 — CT drift-source diagnosis (EXECUTING 2026-08-12)

**Status: the required plan discussion with Ryan was HELD 2026-08-12 (second
session; ledger §"Phase-15 PLAN DISCUSSION HELD"). Implementation is
authorized and IN PROGRESS. Session decisions: ΣΓ⃗ per band added to the
monitor; FOUR conditionals authorized (subtasks (c)–(f) below); burst
statistics REPORT-ONLY; σ stays 0.04R (under-damping near σ_stab disfavored:
n1_s050 at 25% margin drifts too, min_sr↔CT null, and detrended
r(s_k, min_sr) has no consistent sign across arms). Merge cancellation
(H4b) is ELEVATED to leading candidate for H1's Σ|Γ| decline — merging is
the only always-on Σ|Γ| sink at conserved ΣΓ⃗ and the decline is
uniform-sign across bands; the ΣΓ⃗ columns are its passive discriminator
and subtask (f) its active test.**

## Objective

The viscous uniform-d_front carrier (`p018_ufront_n1_s040_visc`, 45 revs) is
stable (min_sr arrest confirmed) but NOT converged: CT̄ carries a monotone
drift — 2.1% over revs 15–29, still 1.1% over revs 30–44. Since the main
objective is a converged simulation, and 15 extra revs reduced but did not
arrest the drift, we diagnose the drift's SOURCE rather than buy open-ended
extensions. Ryan's framing (2026-08-12): the aperiodic wake fluctuation is
accepted as (likely) physics; the non-stationary MEAN is the convergence
blocker to attack.

## Established facts feeding the design (ledger §2026-08-12)

- Gross particle count is NOT still filling: n_p peaks 228k @ rev 15, dips
  to 208k @ rev 30, then +0.17%/rev over revs 30–44 (parallels the late CT
  rise). Whatever drifts lives in the circulation/σ DISTRIBUTION, not the
  headcount.
- min_sr fluctuates aperiodically (band 0.09–0.17, ~5–10-rev episodes) in
  ALL stable arms; per-rev min_sr↔CT correlation is NULL (detrended r≈+0.08
  at lag 0, no consistent lag structure across arms). min_sr is a pure
  stability tripwire; CT scoring need not condition on it.
- Both n_p and CT are non-monotone over 30 revs — H4 (slow mode, not
  secular drift) has nonzero prior.

## Hypotheses and kill conditions

| # | hypothesis | diagnostic | prediction if TRUE | killed if |
| --- | --- | --- | --- | --- |
| H1 | Far-wake Σ\|Γ\| distribution still equilibrating (disk inflow tracks it) | banded inventory monitor | band-resolved Σ\|Γ\| trends on the ~15-rev CT-drift timescale; arrests when CT arrests | band Σ\|Γ\| stationary while CT still drifts |
| H2 | Truncation-boundary secular artifact (deletion at 4R removes downwash the disk should feel; balance shifts as far-wake rollup evolves) | truncation-depth A/B, 4R vs 6R (ONE dedicated job, ~1 day, NT=36 carrier) — **gated on H1's monitor implicating the outermost band; do not buy blind** | drift rate changes with depth | drift rate depth-independent |
| H3 | Viscous σ-distribution secular evolution (CoreSpreading grows aged cores; merges inherit them; effective far-wake resolution degrades slowly) | σ statistics per band (free rider on the same monitor) | σ quantiles trend with the drift while Σ\|Γ\| does not | σ statistics stationary |
| H4 | Not drift: an undersampled ultra-low-frequency mode (period ≳ 30 revs; "monotone" = one flank) | (a) free: AR/spectral screen on the existing 45-rev CT series; (b) decisive: `_s3` extension to ~60 revs (~22 h) | trend reverses on longer horizon; spectral power at ≳20-rev periods | trend continues monotone through ~60 revs |

**Merge-rate hypothesis (H4b, deferred by Ryan 2026-08-12):** a secular
shift in the merge/deletion balance (cloud density → merge rate →
resolution) could also drive slow CT evolution. Plumbing merge statistics
out of FLOWVPM is more than a monitor-side change, so it is SKIPPED in v1 —
but if H1–H4 all come back killed, or the σ-band statistics (H3) trend
without CoreSpreading being able to explain the magnitude, CHECK THIS
HYPOTHESIS before inventing new ones.

## Instrumentation design (Ryan-approved 2026-08-12)

ONE monitor-side src change, bundled, appended LAST in the monitor tuple so
existing CSV names/indices are untouched (same discipline as
WakeHealthMonitor; verify bit-identity on physics):

1. **`WakeInventoryMonitor`** — per step, per band: particle count, Σ|Γ|,
   σ statistics. **Bands (decided): 5 in z/R — 0–0.5, 0.5–1, 1–2, 2–3,
   3–4 — plus a radial split (r/R ≶ 1) on the nearest band** to separate
   root/tip contributions.
2. **σ statistics (decided): include quantiles** (recommend p5/p50/p95 per
   band, plus mean; quantiles are robust where max is not).
3. **min_sr attribution columns** in WakeHealthMonitor (or alongside):
   **p1-percentile of σ/σ_shed** (field-level contraction measure — if IT
   oscillates, the fluctuation is physical, not order-statistic churn) and
   **argmin particle position** (does the same region persist through a
   dip?).
4. **NO merge statistics in v1** (see H4b above).

**Sequencing (Ryan-decided): the instrumentation GATES the pinned dt
ladder** — implement + pass the pre-submission gate (formulation_test, full
runtests since src changed, cluster md5 verify) BEFORE the NT=72/144 rungs
are submitted, so every ladder rung carries the diagnostics for free and
the NT=72 rung doubles as the H1/H3 measurement run. Dedicated drift jobs
only if the monitors demand them (H2's A/B, H4's extension are the only
candidates).

## Subtasks

- [x] (a) DONE 2026-08-12 (subagent; full numbers ledger §2026-08-12
      "Phase-15 subtask (a)"): H4-as-posed INDISTINGUISHABLE from drift at
      35 revs (sinusoid wins AIC by ~1, loses BIC; band power p=0.25; no
      coherent slow mode across arms). TWO structural findings that should
      reshape the plan discussion: (1) per-rev CT̄ is burst-rectification —
      r(mean, within-rev std)=0.61, quiet-limit CT≈0.0730 ⇒ candidate
      convergence statistics = quiet-epoch CT level + burst amplitude,
      SEPARATELY; (2) **restart confound** — the rev-30 block step sits
      exactly at the stitch and within-rev std collapses across it ⇒ the
      "+2.3% offset" revision carries a restart caveat; verify with one
      restart placed mid-quiet-epoch (does the boundary reset burst
      state?). Mean-based H4 discrimination would need +30–45 revs (to rev
      75–90).
- [x] (a2) DONE 2026-08-12 (subagent; full table ledger §2026-08-12 "H1
      two-epoch VTK check"): **H1 SUPPORTED** — rev 30 → rev 45, total
      Σ|Γ| −11.4% (63 sd) with count +3%; all bands −8 to −16.5%; biggest
      movers band 2–3R and the near-disk outboard cell (z<0.5R, r>1R,
      −95% = tightening slipstream, plausible CT link). H3 partial
      (σ_p5 +5–12% but never without Σ|Γ| motion). **H2 concretized: hard
      deletion boundary measured at x=3.49R / r=1.50R** — outer monitor
      band must end there (4R+ empty by construction), and the
      banner-vs-measured depth discrepancy (4R vs 3.49R) must be
      reconciled against the driver before designing the H2 A/B.
      Two-epoch caveat: a trend claim still needs the per-step monitor
      (or a third epoch).
- [ ] (b) Monitor implementation per the design above + pre-submission
      gate. **IN PROGRESS 2026-08-12** (ledger §"Phase-15 instrumentation
      IMPLEMENTED"): `WakeInventoryMonitor` (count, Σ|Γ|, ΣΓ⃗, σ
      mean/p5/p50/p95 per cell; outer band at the RECONCILED 3.5R
      boundary — banner "4R" = cylinder length incl. 0.5R upstream; z
      bands + r/R≶1 split; `outside` catch-all so the partition conserves
      totals) + `WakeHealthMonitor(attribution=true)` p1/argmin columns;
      unit tests green; formulation_test + NT=6 bit-identity smoke pair +
      full runtests + cluster md5 pending. Burst statistics live in the
      ANALYSIS layer: `scripts/p018_burst_stats.py` (report-only;
      validated on B0 + n1_s040; s* pooled per carrier FAMILY only).
- [ ] (c) AUTHORIZED (Ryan 2026-08-12): truncation-depth A/B (H2) —
      `TRUNCATION_DEPTH_R=6` ⇒ 5.5R downstream vs production 3.5R, NT=36
      carrier, ~1 day. Submit after (b) deploys.
      **FIRST READ 2026-08-14 (job 13157751 TIMEOUT at 48 h, step
      1062/1080 = rev 29.5 — matched window fully covered; local harvest
      `data/p018_trunc6_n1_s040_visc/monitors/`): H2 NOT KILLED — drift
      is depth-DEPENDENT.** Matched 15–29: trunc6 CT̄ 0.074411, block
      drift **−3.56% monotone DOWN** (0.0760→0.0740→0.0733) vs carrier
      3.5R CT̄ 0.075298, **+2.08% monotone UP**; late window 22–29:
      0.073619 vs 0.075691 = **−2.7% level effect**; trunc6 quiet-rev
      mean 0.072068 / quiet-limit regression 0.070970 vs carrier quiet
      limit ≈0.0730; M2 vs carrier 15–29: ε_max 1.87%, ε_rms 0.71%.
      **CONFOUND (disclosed): the 5.5R domain never finished filling in
      the window** — n_p still rising 6–9k/rev through rev 27 (carrier
      saturates ~rev 12), so the downward drift mixes fill transient
      with depth response; the clean claims are (i) the carrier's upward
      drift is NOT depth-independent, (ii) deeper truncation LOWERS CT.
      **Deconfound attempt: `p018_trunc6_n1_s040_visc_s2` = job
      13179250** (warm @1061, SETTLE 36 ⇒ rev 44) — **IGNITED rev ~34.5
      (2026-08-15): n_p 344k→13k, min_sigma pinned 9.49e-5, exit 0 with
      garbage. NOT a seam artifact — the cold 5.5R run's |Γ|/σ² had
      already departed the carrier envelope pre-handoff (5.3e3 vs 1.6e3
      at rev 29; carrier stays O(1e3) through rev 44) ⇒ the 3.5R
      deletion boundary is LOAD-BEARING for stability** (aged-tail
      removal pre-empts Γ-ignition; ledger §2026-08-15). H2
      settled-window quantification is BLOCKED at 6R with production
      knobs; the cold-run read above is the H2 result of record.
      Options if Ryan wants more depth resolution (his call, staged in
      the item-file STATUS REPORT): a 4.5R arm (shorter aged tail, may
      stay stable), or accept the cold-run read + carry truncation as a
      one-sided level term (−1.2..−2.7%, fill-confounded).
- [ ] (d) AUTHORIZED (Ryan 2026-08-12): `n1_s040_visc_s3` extension to
      ~60 revs (H4 / more quiet-epoch support). Submit after (b) deploys.
- [ ] (e) AUTHORIZED (Ryan 2026-08-12): restart-seam burst check.
      PROTOCOL REVISION (2026-08-12, constrained by retention): quiet
      epochs on the carrier are revs 31–34 and 39–41 (s_k table via
      `p018_burst_stats.py`; family s*=0.0014), but retained restart
      sets exist only at rev ~30 (step 1079) and rev ~45 (step 1620) —
      no mid-quiet set survives. Read the seam from (d) instead, FREE:
      `_s3` restarts at step 1619 (rev ~45; CSV shows 1620 = S+1, VTU retained 1610–1619); if
      within-rev std collapses across THAT seam like it did at rev 30
      (0.007→0.0005), the restart truncates burst state — confound
      confirmed without quiet placement. Only if ambiguous: protect a
      mid-quiet VTK set from `_s3` (revs 46+ quiet epoch) and run one
      dedicated short restart from it.
- [ ] (f) AUTHORIZED (Ryan 2026-08-12): merge-off A/B (H4b active test) —
      warm restart of `n1_s040_visc_s2` from the SAME step 1619 as (d)
      with `MERGE_PARTICLES=false` (driver knob, sets merge every=0),
      ~5 revs (180 steps; n_p headroom 215k of 500k). **(d) is its stock
      control**: same restart step, only the merge knob differs. Score:
      per-band Σ|Γ| slope (both arms carry the new inventory monitor) +
      CT level over the matched 5 revs. Σ|Γ| decline + CT drift stop ⇒
      merging owns the drift; unchanged ⇒ H4b killed, H1-as-convection
      stands. Diagnostic only, never production.

## Related, queued post-ladder (note in item file)

Numerical-vs-physical discriminator for the min_sr/CT fluctuation: compare
band, episode timescale, and CT per-rev std across NT 36/72/144 — a
numerical mode should move with Δt, a physical wake mode should not. Comes
free with the ladder; record the verdict in the ledger.

**VTK deletion log 2026-08-12 ~23:15 MDT (200G budget sweep):** older-step
VTK of live conditionals swept via `p018_vtk_sweeper.sh` — mergeoff 873MB
(kept 1640–1649), trunc6 1313MB (kept 133–142), `_s3` 915MB (kept
1641–1650); newest-10 restartable sets + all monitors/CSVs retained.

## 2026-08-13 — (f) merge-off A/B harvested: H4b KILLED

`p018_mergeoff_n1_s040_visc` (13157753) COMPLETED cleanly (3h56, exit 0
trusted: tail min_sr 0.140, global min 0.117, peak max_u 57.5, no
ignition), reaching step 1727 = rev 47.97 (+3 revs from the rev-45
handoff). Merge-off verified genuine at launch (n_p excess vs `_s3`
growing ~+44/step). Matched window vs `_s3` = steps 1620–1727 (revs
45–48):

- **Banded Σ|Γ| slope (monitor05, %/rev of band mean):** total mergeoff
  **+3.104** vs `_s3` **+3.001** — trajectory UNCHANGED (Δ = +0.10%/rev
  on a +3%/rev signal). Dominant bands match to ≤0.2%/rev
  (z1p0_2p0 −0.53 vs −0.73; z2p0_3p0 +10.40 vs +10.35; z3p0_3p5 −3.98
  vs −4.07). Meanwhile the particle sink is real: n_p slope +2.00 vs
  +1.22 %/rev — merging removes ~0.8%/rev of particles with no
  measurable Σ|Γ| footprint.
- **CT matched revs 45–47:** mergeoff 0.075284 vs `_s3` 0.075118 =
  **+0.22%** (within per-rev std ~2.7%). M2 (45–48): ε_max 1.56%,
  ε_rms 0.75%.
- **Verdict per the pre-registered rule ("unchanged ⇒ H4b killed"):
  H4b KILLED — merge-cancellation does not own the Σ|Γ| decline.**
  Caveats: 3-rev window (resolution ~0.1%/rev slope difference — still
  well below the drift-scale signal H4b was to explain); the window
  itself sits in a Σ|Γ|-rising (burst) phase. Side-benefit for budget
  term 8: merging ΔCT +0.22% over 3 revs (target ≤0.25% — consistent,
  not yet a settled-window claim).

## 2026-08-13 — (d)+(e) `_s3` harvested (13157752, COMPLETED 18h15, healthy)

Exit 0 TRUSTED: tail min_sr 0.1488 (recovering), max_u 22.6, n_p 222,584;
CT CSVs written. Coverage steps 1620–2087 = revs 45–58.

**(d) H4 trend through rev 58 — grand mean drifts ON, quiet baseline
STATIONARY:**
- M1 revs 45–57: CT̄ 0.076806 [0.075467, 0.078458], 5-rev blocks
  0.075435 → 0.077829, block drift **+3.117% monotone** (vs _s2's 30–44:
  0.073632, +1.107%) ⇒ grand-mean drift CONTINUES and STEEPENS — no
  arrest, no reversal.
- Burst/quiet decomposition (stitched `_s2+_s3@1619`, revs 30–57, n=27):
  quiet-limit regression **0.072975** vs `_s2`-alone (30–44) **0.072979**
  — the quiet baseline is UNCHANGED across 30→58; quiet-rev means
  0.073332 vs 0.073171. Burst rectification r=+0.611, incidence 0.74,
  quiet-series drift +1.676%±1.380% (half-split +0.21%). **Reading: the
  apparent CT drift is burst rectification (episodes growing more
  frequent/energetic, e.g. rev-53 mean 0.082433, ptp 0.021), not a
  moving baseline — quiet limit stays ≈0.0730.** H4's aperiodic mode
  persists at 45–58 with ~5–10-rev episode spacing.

**Seam check (rev-45 stitch): CLEAN.** Per-step CFz continuous across
1619→1620 (consecutive-step ramp +87%/+73%/+48% — normal blade-passage
ramp, no discontinuity); rev-45 within-rev ptp 0.00794 within neighbor
range (0.0096–0.0225) — NO std collapse like the rev-30 concern.

**(e) banded inventory (monitor05, revs 45–58, %/rev of band mean):**
total Σ|Γ| **+0.451 %/rev RISING** (4.188→4.908) — the earlier-feared
monotone Σ|Γ| decline is NOT present in this window. Bands: rin +0.53,
z0p5_1p0 +0.38, z1p0_2p0 −0.91, z2p0_3p0 +0.98, z3p0_3p5 +0.78, outside
+3.13. With H4b killed (merge-off A/B above), the Σ|Γ| story reads as
slow far-wake population growth, not merge cancellation.

**VTK deletion log 2026-08-13 ~16:00 MDT:** post-harvest sweeps —
`_s3` 19,539MB freed (kept 2078–2087, the `_s4`/rlxf-comparison chain
source), mergeoff 3488MB freed (kept 1718–1727). CSVs/TOML/monitors
retained cluster-side and harvested locally.

**VTK deletion log 2026-08-13 ~19:30 MDT:** live sweep (disk alert) —
trunc6 22,389MB freed (kept 692–701); newest-10 + monitors/CSVs retained.
