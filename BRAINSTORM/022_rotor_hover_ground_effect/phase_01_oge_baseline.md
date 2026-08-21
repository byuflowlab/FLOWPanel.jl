# Phase 1 — OGE baseline at the 022 operating point

**Objective:** CT_OGE at RPM 6000 / ρ 1.16 / R 0.1195 with the 018 carrier
knobs, from `examples/rotor_hover_ground_effect.jl` with `GROUND_ENABLE=false`
(matched-settings contract, ruling 4). Two mesh rungs; 56_57 lands early as a
cross-check while 45_185_ct4 (headline) runs.

## Cases

| tag | mesh | knobs | job | status |
|---|---|---|---|---|
| p022_oge_fine | 45_185_ct4 | carrier (see item Fixed operating point), GROUND_ENABLE=false, TRUNC_RADIUS_R=1.5, NREVS=30 | — | pending |
| p022_oge_coarse | 56_57 | same | — | pending |

Carrier ENV (both): `NT=36 OVERLAP=2.75 P_PER_STEP=12 DAS_UNIFORM_DSIGMA=3.4
NWAKEROWS=1 MERGE_R_FACTOR=0.0055 RELAX_RLXF=0.3 CORE_SPREADING_ACTIVE=true
WAKE_CORE_BETA=1e9 PARTICLE_SHEDDING=sigma_overlap SETTLE_REVS=20
RPM=6000 RHO=1.16 ROTOR_R=0.1195`.

## Decision (pre-registered)

- Harvest M1 per `decision_rules.md` (health gates first; window in hover,
  per-rev blocks stationary). CT_OGE_fine is the Phase-2 denominator.
- Cross-check: |CT_coarse − CT_fine|/CT_fine reported (expected few %, per
  018 mesh experience); a >10% gap = mesh/shedding pathology, investigate
  before Phase 2 harvest.
- Context anchors (not acceptance): 018 quiet-limit CT ≈ 0.0730 at
  RPM 5400/ρ 1.179; CT is dimensionless so the operating-point change should
  move it only through Re/startup differences (expect same ballpark).

## Exit criteria

CT_OGE (cycle-mean ± std) recorded in `ledger.md` for both rungs, health
gates green.

## Log

- 2026-08-18 (staging): cases defined; awaiting Phase 0 PASS then submission.
