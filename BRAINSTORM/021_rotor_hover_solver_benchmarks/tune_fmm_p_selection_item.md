# tune_fmm p-request semantics (epsilon construction) — staged fix item

**Status: STAGED 2026-08-19 — do NOT execute without Ryan's go.** Scoped
during the cached-nearfield feedback round; deliberately NOT explored further
in-session. **Interim ruling (Ryan, same day): the tuner's full p is KEPT
for now** — the worst-case point genuinely needs it under the requested
max-abs criterion, and the ~2× apply cost (R2: 77.9 ms at p=15 vs
38.8–43.4 ms at p=6–8) is accepted as insurance for that point. The HPC
retune adopts tuner knobs AS-IS (leaf, MAC, p); no descending-p step. This
item is therefore purely a FUTURE improvement to the epsilon-construction
semantics, not a pending campaign dependency.

## Symptom (measured 2026-08-19, local 4T)

`tune_fmm` with the campaign request
`error_tolerance = PowerAbsolutePotential(0.1·1e-6·rms_b)` (ABSOLUTE epsilon
derived from the 1e-6 relative BC target, 0.1 safety included) selects
p=14 (R1, cached-tuned leaf 342) / p=15 (R2, cached-tuned leaf 196), while
the campaign-facing REL-L2 apply error vs a p=20 reference at the SAME
leaf/MAC meets the 1e-6 relative target at **p=6** (R2: 3.1e-7; p=8:
9.8e-9). Cost: R2 apply 77.9 ms at p=15 vs 38.8–43.4 ms at p=6–8
(~1.8–2×). See "Corrected diagnosis" below: the algorithm satisfies its
max-abs request nearly minimally — the mismatch is between that request and
the campaign metric.

## Corrected diagnosis (max-error check, R2, 2026-08-19 — Ryan-requested)

Measured MAX pointwise absolute apply error at leaf=196 vs the requested
epsilon 1.18e-9: p=6 → 4.90e-7 (**414× OVER**), p=8 → 2.86e-8 (24× over),
p=15 → 5.02e-11 (meets, 24× margin ⇒ minimal compliant order ≈ p=12–13 at
~1.4 decades/Δp=2). **So tune_fmm is nearly correct for the max-abs
criterion it was GIVEN — only ~1 decade of bound conservatism.** The
apparent 4-decade over-selection lives in the REQUEST construction:

1. the 0.1 safety factor (1 decade, deliberate);
2. **max-abs-per-point vs the campaign rel-L2 BC intent** (~1.5–2 decades):
   the max is set by a few worst points near close m2l pairs while BC
   certification averages in L2 (p=6: max 414× over epsilon yet rel-L2
   3.1e-7, inside the 1e-6 target);
3. **normalization**: epsilon scales with rms_b (0.0118) but the operator is
   linear, so an "absolute" tolerance is implicitly relative to the strength
   magnitude present at tune time; the O(1) probe potential makes the
   effective pointwise-relative request ≈1.2e-9, ~3 decades tighter than the
   stated 1e-6 relative intent.

## Candidate directions (pick after Ryan's go; smallest first)

- **Fix the campaign epsilon convention** (likely the whole fix): construct
  the PowerAbsolutePotential epsilon from the metric actually certified
  (rel-L2 BC), not a per-point max at rms_b scale — or use a relative error
  model if FastMultipole exposes one suited to the certified metric.
- **Measured-verify step in tune_fmm** (robust regardless of convention):
  after the bound-driven pick, descend p at the tuned leaf with actual
  applies until the requested criterion is violated; return the last
  passing p.
- **Bound calibration** is NOT the lever (only ~1 decade of slack there).

## Evidence pointers

- Sweep tables: `phase_02_single_step_benchmarks.md` Log 2026-08-19
  "(later)" entry (R1 p-sweep @ leaf 342; R2 p-sweep @ leaf 196, m2l=516).
- Phase-1 bound-margin history: `phase_01_consistency.md` (sliver-fix floor
  work, ~50× margin at p=10; tuner degenerates if requested tol < floor).
- Related landed code: cached-path tuning (`tune_nearfield_cache`,
  FastMultipole 204188a + 1ec0af9), `estimate_nearfield_cache` (0ef4e83).

## Non-goals

- No changes to the campaign's frozen Phase-1 knobs or certification metric.
- Not blocking the HPC retune: that uses tuner (leaf, MAC) + descending
  measured p-sweep + per-run BC certification, independent of this fix.
