# Phase 2 — IGE first light + headline comparison

**Objective:** CT_IGE at h/R = 1.0 with the best-guess ground (disc 4R,
panel_length 0.15R, truncation radius 3R, particle policy `none`), same
carrier knobs and mesh rungs as Phase 1; report CT_IGE/CT_OGE per rung.

## Cases

| tag | mesh | knobs | job | status |
|---|---|---|---|---|
| p022_ige_fine | 45_185_ct4 | carrier + GROUND_ENABLE=true GROUND_H_R=1.0 GROUND_RADIUS_R=4.0 GROUND_PANEL_LENGTH_R=0.15 TRUNC_RADIUS_R=3.0 GROUND_PARTICLE_POLICY=none, NREVS=30 | — | pending |
| p022_ige_coarse | 56_57 | same | — | pending |

## Decision (pre-registered)

- Health gates per `decision_rules.md`, PLUS ground-specific: GS
  `gs_nonconverged = 0`; tangency RMS bounded; below-ground census reported.
- Headline: CT_IGE/CT_OGE per rung with quadrature error bars. Coarse-rung
  ratio is the leading indicator while fine runs; the fine-rung ratio is the
  reported number.
- Sanity anchor: momentum theory h/R=1 ⇒ ratio ≈ 1.07. Ratio < 1 or > ~1.2
  ⇒ diagnose (truncation amputating the wall jet? ground disc too small?
  below-ground Γ accumulation?) before reporting.
- IGE settling: if per-rev means still trend at end-of-run, extend/chain
  (mind the warm-start-with-ground driver guard: rerun longer cold instead).

## Exit criteria

CT_IGE and the IGE/OGE ratio (both rungs) in `ledger.md`; headline paragraph
in the item's Current status; offer Ryan a notebook entry at this milestone.

## Log

- 2026-08-18 (staging): cases defined; awaiting Phase 0 PASS then submission.
- 2026-08-19: **p022_ige_coarse (13207682) blew up and was cancelled** (Ryan's
  call) at 16 h 28 m, step 636/1007 (rev ≈ 17.64). Full signature table in
  `ledger.md`. Headlines: CF_x −0.074 → −371.9 in one step; ground tangency
  RMS(U·n) 0.22 → 0.842 and max|U·n| 4.9 → 20.17 at step 630; below-ground
  sum|Γ| 0.38 → 9.859; walltime/step 134 s → 633 s. Steps ≤ 635 are clean.

  Leading suspect is the below-ground particle population under the leave-be
  policy (ruling 2): the census grew monotonically 456 → 6362 over steps
  240–630 while sum|Γ| stayed ~0.35, i.e. accumulation was tolerated for ~390
  steps before something tipped. This is exactly the trigger the item's
  contingency chain names for pulling Phase 4's cull ablation earlier —
  **but the second half of that trigger ("AND CT drifts") was not met**: CT was
  stationary at −0.078 ± 0.004 right up to the event. So this reads as a sudden
  instability, not a slow corruption, and the cull ablation is indicated but
  not yet proven to be the fix.

  **RESOLVED 2026-08-20: NOT mesh-independent.** p022_ige_fine (13207681)
  passed the equivalent rev and was healthy at step 885/1007 (rev ≈ 24.6):
  tangency RMS 0.22, below-ground census PLATEAUED at ~4600 (vs the coarse
  rung's monotone climb), sum|Γ| stable at ~0.27. The fine rung tolerates the
  same leave-be policy the coarse rung died under. Also p022_oge_fine
  (13207679) COMPLETED: **CT_OGE fine = 0.07480 ± 0.00152** (final 10 revs,
  in-hover, all finite) — the Phase-2 denominator exists.

  (Superseded original note:) **Not established as mesh-independent.** p022_ige_fine (13207681) was at step
  445 with comparable tangency (RMS ≈ 0.24, 2714 below-ground) and zero |CF_x|
  excursions when checked on 2026-08-19; it is still running and reaches the
  equivalent rev around step 636 of its own schedule. Watch it there before
  concluding anything about the mechanism.

  Phase 2 gate is **NOT passed**: the coarse rung has no harvestable M1 window
  (any window reaching rev 18 includes the blow-up). A preliminary matched read
  over revs 7–17 is recorded in `ledger.md` (ratio 1.045) and is explicitly not
  an M1 harvest. The fine rung is still the live path to the Phase 2 headline.
- 2026-08-19 (ops): OGE coarse (13207680) will also miss its window — measured
  ~87 s/step means the 1007-step schedule needs 25–31 h against a 24 h wall, so
  it lands around rev 25 of 28. Recorded in `ops_reference.md`; Phase 5 sizes
  its coarse runs at 48 h as a result.
