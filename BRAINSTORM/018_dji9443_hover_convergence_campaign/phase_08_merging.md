# Phase 8 — Merging Null + Speedup

**Objective:** demonstrate particle merging doesn't noticeably affect the
solution at converged settings, and report the speedup it buys.

## Case (needs Phase 4's σ\*)

| tag | knobs | time |
| --- | --- | --- |
| `p018_nomerge` | final settings + MERGE_PARTICLES=false (pass final-settings env at submission) | 48 h+ (unmerged count grows the N-body cost; budget ~1.5× the merged twin) |

**Particle-cap watch:** without merging the count can exceed the driver's
hard-coded `max_particles=500_000` at σ\* — if the count trends toward the cap,
raise the constant (`examples/rotor_hover_pressure_comparison.jl:343`,
one-line; redeploy + md5-verify) and resubmit. Do not silently let the run die
at the cap.

## Deliverables

1. |ΔCT̄| vs the matched merge-ON final-settings run — target ≤ 0.25%;
   ε_Γ ≤ 0.5%. → error-budget term 8.
2. **Speedup report:** wall time ratio, settled particle-count ratio, and the
   per-step timing trend (merge-ON vs OFF). Cross-reference
   `data/rotor_hover_pressure_comparison/particles_merging_savings.csv`
   (legacy recording) if useful.
3. Confirm `MERGE_R_FACTOR` scaling held (r_merge/σ ≈ 0.138 at σ\*: factor
   ≈ 0.0053 for L1) — an over-aggressive merge at refined σ would undo the
   Phase-4 refinement (012 caveat).

If the null FAILS (Δ > 0.25%): halve MERGE_R_FACTOR and retest before
concluding; if it still fails, production runs go merge-OFF and the speedup
claim is dropped (report honestly).

## Exit criteria

Null + speedup numbers in the ledger.

## Log

(append dated entries here)
