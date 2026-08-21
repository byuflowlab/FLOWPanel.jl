# Phase 3 — below-ground particle policy (refined, Ryan 2026-08-20)

**Objective:** stop particles from passing through the ground and stabilize
the IGE runs, per Ryan's refined two-step policy. Promoted from Phase 4 to
Phase 3 (order swapped with ground convergence, Ryan 2026-08-20) after the
p022_ige_coarse blow-up implicated below-ground particle accumulation
(monotone census 456 → 6362 over steps 240–630, then sum|Γ| 0.38 → 9.86 and
CF_x −0.074 → −371.9 in one step; see `phase_02_ige_first_light.md`).

## The refined policy (supersedes ruling 2's pure leave-be)

1. **Near-ground vertical velocity cutoff.** For particles whose axial
   (vertical) coordinate is within `GROUND_DAMP_BAND_R * R` of the ground
   (start value **R/10**), scale the ground-ward axial velocity component by
   the **linear** factor f = d/band (d = height above the plane; f = 0 at and
   below it). Receding motion is never damped, so strays can recover.
   Implemented as a driver-local `propagate!` overwrite that edits the stored
   particle U just before the Euler position update (`_euler` convects with
   exactly that U). Discrete-step guarantee: with dt = 2.78e-4 s and band =
   R/10, a particle cannot cross the plane below ~43 m/s axial velocity.
   Not foolproof by design ⇒ measure (A/B below).
2. **Truncation floored at the ground.** `GROUND_PARTICLE_POLICY=truncate`
   raises the truncation-cylinder floor to coincide with the ground: cylinder
   length becomes `ground_x + 0.5R` (≈ 1.46R at h/R = 1) instead of 4R, so
   anything that still crosses is deleted at the next maintenance pass.

Diagnostics: `*_ground_diagnostics.csv` gains `n_damped`, `n_inband` columns
and — with `*_gs_convergence.csv` — is now written **incrementally** per step
(2026-08-20 lesson: post-run dumps lose everything on a walled/cancelled run;
13207682's census survived only in the slurm log).

## Cases

| tag | mesh | policy knobs | job | wall | status |
|---|---|---|---|---|---|
| (A/B "without") p022_ige_coarse | 56_57 | damp 0, none | 13207682 (cancelled @636) | — | census in slurm log, steps ≤635 clean |
| p022_ige_coarse_damp | 56_57 | GROUND_DAMP_BAND_R=0.1, GROUND_PARTICLE_POLICY=none | 13246557 | 48 h | submitted 2026-08-20 |
| p022_ige_coarse_trunc | 56_57 | GROUND_DAMP_BAND_R=0.1, GROUND_PARTICLE_POLICY=truncate | 13246558 | 48 h | submitted 2026-08-20 |

All other knobs identical to p022_ige_coarse (carrier + ground first-light).
`_damp` keeps `none` deliberately so the below-ground census still counts
pass-throughs — that is the A/B observable. `_trunc` is the production
candidate (damping + floor).

## Decision (pre-registered)

- **A/B pass-through comparison:** below-ground census (count and sum|Γ|)
  of `_damp` vs the cancelled run's log over the matched window (revs ≤ 17.5,
  steps ≤ 630). Damping is EFFECTIVE if it cuts the below-ground count by
  ≥ 10× over that window; PARTIALLY effective in (2×, 10×); ineffective below
  2× (then the linear band is the wrong tool — escalate to Ryan before
  inventing alternatives).
- **Stability:** does `_damp` (and `_trunc`) survive past rev 17.6, where
  13207682 ignited? Survival to end-of-schedule with stationary per-rev CT =
  the accumulation hypothesis is supported AND the fix works; a blow-up at
  the same rev despite an empty below-ground census = the hypothesis is
  wrong, and the mechanism hunt reopens (do not tune further knobs blindly).
- **CT effect:** ΔCT̄(damp vs pre-blow-up window of 13207682) and
  ΔCT̄(trunc vs damp) reported with cycle-std; a policy that shifts CT̄ by
  more than the combined std is a modeling decision for Ryan, not an
  automatic accept.
- Production policy for Phases 4/5 = whatever passes here (expected:
  damp band R/10 + truncate). Record as a standing-ruling candidate.

## Exit criteria

Pass-through A/B numbers + survival verdict + chosen production policy in
`ledger.md` and the item's Current status; Phase 4/5 arms updated to carry
the chosen policy knobs.

## Log

- 2026-08-18 (staging, as old Phase 4): cull-vs-none ablation defined.
- 2026-08-20: Ryan refined the policy (velocity cutoff + truncate) and swapped
  this phase ahead of ground convergence; driver implements
  `GROUND_DAMP_BAND_R` (linear, ground-ward only) and
  `GROUND_PARTICLE_POLICY=truncate`; diagnostics CSVs made incremental.
