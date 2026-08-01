# Phase 7 — GreenReconstruction Delta

**Objective:** measure what the wake→body transfer scheme does at converged
settings. Motivation: item 015 (wing) showed `GreenReconstruction` removes 93%
of the small-Das Dirichlet differential — the default `VelocityThroughSources`
equivalent-source transfer is itself a major error source there. The rotor
measurement at legacy settings (+0.93% at η=0.2) was below cycle scatter and
unresolved. Reporting ΔCT̄ and ΔΓ̄(r/R) at converged settings is a **required
deliverable regardless of size**.

## Cases (need Phase 4's σ\*)

| tag | knobs | time |
| --- | --- | --- |
| `p018_green` | final settings + RHPC_FORMULATION=green (pass final-settings env at submission) | 48 h+ (~1.5× velocity cost; warm start if validated) |
| `p018_green_coarse` (optional) | B0 carrier + green | 24–48 h — run if the final-settings Δ is >0.5%, to show the Δ trend with σ |

Memory: green build peak ~30 GB at this mesh (on top of the particle wake) —
64G per ruling 5; `GreenReconstruction(gauge=:area_mean,
recompute_interval=1)` as in 2e.

## Decision

|ΔCT̄| ≤ 1% ⇒ velocity remains production; Δ (CT and Γ) → error-budget term 7.
> 1% ⇒ **Green becomes the production formulation**: Phase 9 runs green, and
the velocity-based ladders are re-examined for whether any rung decision would
flip (spot-rerun the decisive rung in green if so — flag to Ryan before
spending it).

## Exit criteria

ΔCT̄ and ΔΓ̄(r/R) on matched windows in the ledger; production formulation
decided.

## Log

(append dated entries here)
