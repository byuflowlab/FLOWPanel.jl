# Phase 4 — ground representation convergence ladders

**Objective:** show CT_IGE is converged in the ground-representation knobs.
One axis at a time off the Phase-2 first-light carrier **plus the Phase-3
production particle policy** (order swapped with the old Phase 3 by Ryan
2026-08-20 — the policy must be settled before these ladders mean anything).
Coarse mesh rung first for screening; promote the selected values to a
fine-rung confirm run.

## Ladders (pre-registered)

| axis | rungs (carrier value in bold) | ENV knob |
|---|---|---|
| ground disc radius | 2R, **4R**, 6R | `GROUND_RADIUS_R` |
| ground panel length | 0.3R, **0.15R**, 0.075R | `GROUND_PANEL_LENGTH_R` |
| truncation radius | 1.5R, **3R**, 4.5R | `TRUNC_RADIUS_R` |

Notes:
- Disc radius and truncation radius interact (the wall jet must fit inside
  both); if the disc-radius ladder fails to plateau, widen truncation first.
- CONSTRAINT: TRUNC_RADIUS_R ≤ GROUND_RADIUS_R always (the truncation
  cylinder is the horizontal domain boundary — particles must be deleted
  before leaving the paneled ground support; the driver warns if violated).
  The 4.5R truncation rung therefore pairs with the 6R disc, not the 4R
  carrier disc.
- Ground panel count grows ~ (radius/panel_length)²: 4R/0.075R ≈ O(10k)
  triangles — still cheap for FlatGroundSolver (O(N)) but the FMM influence
  cost grows; watch wall time.
- Truncation depth stays 4R (particles below ground are the leave-be policy's
  business, Phase 4).

## Decision (pre-registered)

Per `decision_rules.md`: a knob is converged when doubling (halving for
panel length) moves M1 cycle-mean < 0.5% and within combined cycle-std.

## Exit criteria

Selected converged ground prescription recorded in the item file + `ledger.md`;
fine-rung confirm run matches the coarse-rung selection within tolerance.

## Log

- 2026-08-18 (staging): ladders defined; blocked on Phase 2 first light.
- 2026-08-20: renumbered from Phase 3 to Phase 4 (swap with particle policy,
  Ryan); now additionally blocked on Phase 3's policy verdict.
