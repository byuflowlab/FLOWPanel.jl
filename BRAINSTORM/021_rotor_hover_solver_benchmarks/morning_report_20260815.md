# 021 morning report — autonomous session 2026-08-14 evening (for Ryan, 2026-08-15)

All three work-queue items advanced; every judgment is logged in
`phase_01_consistency.md` (Log entries "2026-08-14 evening"). Nothing was
written to the lab notebook. Local runs stayed ≤4 threads.

## 1. PENDING CHECK resolved: tune_fmm transfers SAFELY to FGS — but keep production FGS knobs

Empirical check (`benchmark/rotor_hover_solver_phase1_fgs_check.jl` →
`fgs_check.csv`): R1–R3 × {fgs, fgmres_fgs} × three knob candidates (prod
p10/MAC0.4/leaf150; per-rung tuned; hybrid = tuned p+MAC with leaf 150), all
judged by DIRECT-evaluated true residual vs the frozen b.

- **All 18 runs meet rel 1e-6** — the replayed source-topology target tree
  introduces NO new residual floor; internal convergence is honest.
- But tune_fmm's small leaf is bad for FGS *cost*: leaf is also the
  Gauss–Seidel block size (R2: 12 outer iterations at leaf 12 vs 6 at
  leaf 150; 20–65% more wall time). The tuner has no model of that role.
- Production knobs fastest or tied at every rung; hybrid within ±3%.

**PROPOSED FGS-family freeze** (in `decision_rules.md`, pending you):
`fgs` = p10/MAC0.4/leaf150, inner 20, rlx 1.0, shrink=true, tolerance =
1e-6·rms(b) absolute (its max-abs metric), maxit 300 (≤12 used), same knobs
R4–R7, guard = the per-run direct-evaluated residual (an O(N²) apply — no
dense G — feasible at every rung incl. R6–R7). `fgmres_fgs` = same
preconditioner knobs (sweeps 1, inner 2), tuned apply backend, rtol 1e-6 /
atol 1e-14 / memory 50 / itmax 500 (9–11 iterations used). Do NOT re-tune
FGS-internal knobs per rung. FGS cost scaling is excellent: 2.6 / 5.7 /
20.0 s (R1/R2/R3, multi @4T) ≈ N^1.15, vs dense LU 2.8 / 18 / 105 s.

## 2. backslash_coupled: NOT x≈0, NOT a layout bug — root-caused, fixed test-gated, roster 7/7 clean

The rel_rms = 1.0 availability rows were the signature of **x = 2·x_ref
exactly** (G·2x_ref − b = b). Probe proved two independent causes:

1. **Script state leak**: the coupled solvers treat entry `body.potential`
   as an EXTERNAL potential and add it to their internally assembled source
   potential; the availability script's `reset!` left the b-capture's
   assembly in `rotor.potential` → doubled RHS. With the potential zeroed,
   `BackslashCoupled` matched the Backslash reference to machine precision
   (rel 5.4e-15) — the coupled DBC solution-column layout is CORRECT.
2. **Real latent src bug**: `calc_bc_dirichlet` ACCUMULATED
   (`RHS .-= potential`) into the coupled solvers' persistent, never-zeroed
   `rhs` — every solve after the first doubled the Dirichlet rows. Hits
   `BackslashCoupled` AND `KrylovCoupled` in any repeated-solve use
   (time-marching). **Fixed** to assignment semantics
   (`src/FLOWPanel_solver.jl`), regression testset "Coupled Dirichlet rhs
   idempotence" added (premise-verified: with the old code temporarily
   restored, ‖x₂‖/‖x₁‖ = 2.0), full `runtests_unit_solver.jl` green.

Availability rerun of record: **backslash_coupled PASS rel 5.43e-15**; all
seven roster configs clean on the Dirichlet campaign case.

## 3. Agreement tables vs Backslash at matched residual (exit criterion 1) — R1–R3 COMPLETE, all 7 configs

New stage `benchmark/rotor_hover_solver_phase1_agreement.jl` →
`agreement.csv`: all 7 configs through the production `steady!` path
(skeleton = rotor_hover_convergence.jl), frozen/PROPOSED settings, shared
tuned post-processing backend, CT from PressureLaplace(edge_difference,
material_derivative)+Force, steady PressureBernoulli+Force, and
KuttaJoukowski (fresh monitors per config; steady-Bernoulli moving-body
caveat accepted because the metric is the cross-config DELTA).

Headline vs provisional thresholds (relL2 ≤ 1e-5, relmax ≤ 1e-4, |dCT| ≤ 0.1%):

| config | relL2 (R1 / R2 / R3) | relmax (R1 / R2 / R3) | worst dCT % | verdict |
| --- | --- | --- | --- | --- |
| krylov_jacobi | 1.1e-6 / 1.2e-6 / 1.6e-6 | 9.4e-6 / 5.2e-6 / 3.4e-5 | 4.5e-5 | PASS |
| krylov_ilu | 1.6e-6 / 6.2e-7 / 1.1e-6 | 3.3e-6 / 1.2e-6 / 2.0e-6 | 2.8e-5 | PASS |
| fgs | 1.8e-8 / 1.8e-8 / 1.1e-8 | 2.4e-8 / 3.5e-8 / 3.1e-8 | 1.0e-6 | PASS (best) |
| fgmres_fgs | 2.1e-6 / 1.8e-6 / 4.0e-6 | 1.2e-5 / 1.3e-5 / 1.9e-5 | 2.7e-4 | PASS |
| krylov_gmres | 1.2e-5 / 1.5e-5 / 1.6e-5 | 4.3e-5 / 5.2e-5 / 5.4e-5 | 7.0e-4 | relL2 MARGINAL |
| backslash_coupled | 1.3e-5 / 1.3e-5 / 1.3e-5 | 1.2e-5 / 1.2e-5 / 1.2e-5 | 8.3e-4 | FLAG (below) |

- **CT agreement is uniformly excellent**: worst |dCT| = 8.3e-4 % — 120×
  inside the 0.1% threshold, across all three monitors.
- **krylov_gmres relL2 1.2–1.5e-5 vs the provisional 1e-5**: this is not a
  solver defect — it is exactly the residual→solution amplification
  κ_eff·rtol with κ_eff ≈ 12–15 (its residual honestly sits at 1e-6).
  **PROPOSED threshold calibration** (decision_rules contract says calibrate
  from observed sensitivity): relL2 ≤ 3e-5, relmax ≤ 1e-4 (unchanged), dCT ≤
  0.1% (unchanged) at the rel-1e-6 matched-residual target. Alternative if
  you prefer the 1e-5 solution bar: tighten the matched target to 3e-7 (adds
  ~1 gmres restart-cycle of iterations; all other configs already comply).
- **backslash_coupled inside steady! — ROOT-CAUSED (same night): a
  formulation-gauge inconsistency, not a solver bug.** Probe chain proved:
  G/LU factors, entry velocity, and BC sources bit-identical between the
  direct call and the steady! path; the entire 4.8e-5 difference is the rhs,
  and it equals `steady!`'s `apply_freestream!` writing the freestream
  scalar potential φ∞ = U∞·x into `body.potential` pre-solve. The
  single-body Dirichlet wrapper ZEROES that workspace (perturbation-potential
  gauge — the production convention), while the coupled solvers honor entry
  potential as an external field (φ∞ enters the Dirichlet rhs). Both
  self-consistent; difference scales with |U∞| (here magVinf = 1e-4) and
  vanishes at exact hover. **Needs your ruling on which convention the two
  paths should share.** Benign at campaign thresholds (CT deltas 8e-4 %).

R3 operational notes: the environment kills long background jobs on this
machine (three kills; 16 GB RAM is also tight for two live 6.6 GB dense Gs),
so R3 ran in per-config chunks with the reference solution persisted to
disk (`agreement_xref_R3.bin`). Also the **ILU pattern guard legitimately
trips at R3** — the leaf-10/MAC-1.0 Barba direct list exceeds the default
512 entries/panel as this sliver-heavy mesh refines; the agreement run used
max_pattern_entries=2048N. That density growth is a Phase 2 knob
consideration for R4+ (pattern nnz/panel is not N-independent on this mesh).

Caveat (R3 only): the PressureLaplace CG hit its 1000-iteration cap at R3
(relative residual ~3e-6 vs its internal target) — Laplace CT at R3 is
slightly under-converged, identically for every config (shared monitor
recipe), so deltas remain comparable; flagged for the HPC sweep (raise
itmax).

## Housekeeping / open items for you

1. **Rulings wanted**: (a) tolerance freeze incl. the FGS-family block and
   MAC≥0.6 rejection (decision_rules PROPOSED); (b) agreement-threshold
   calibration relL2 → 3e-5 (above); (c) which Dirichlet external-potential
   convention the single-body and coupled paths should share inside
   steady!/simulate! (φ∞ discarded vs honored — see the root-caused item).
2. The FastMultipole working tree still carries the three uncommitted fixes
   (FGS tree replay + two tune_fmm) — plus the unrelated dirty
   `scripts/20250404_prediction_accuracy.jl` to segregate before committing.
3. New FLOWPanel changes this session (uncommitted, test-gated):
   `calc_bc_dirichlet` assignment fix + regression testset; agreement/
   fgs-check benchmark scripts; availability-script `reset!` fix.
4. Environment note: long (>~45 min) background jobs on this machine get
   killed mid-run (hit twice at R3 backslash_coupled; 16 GB RAM is also
   tight for two live 6.6 GB dense Gs) — R3 agreement was run in short
   per-config chunks; HPC remains the right home for R4+.
5. Next actions (my ordering): backslash_coupled steady! probe → freeze
   Phase 1 tables on your rulings → R4–R7 TUNE + agreement on HPC → Phase 2.
