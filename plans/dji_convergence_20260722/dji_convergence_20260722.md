# Adaptive DJI 9443 Circulation and CT Convergence Study

Date: 2026-07-22

## Goal

Determine whether the re-fit DJI 9443 airfoils improve bound circulation, then
seek numerical convergence of hover thrust coefficient `CT` while separating
the effects of mesh, topology/formulation, timestep, wake core, target offset,
mu-gradient reconstruction, and Green reconstruction.

Experimental reference: `CT_exp = 0.072` at 5400 RPM. Numerical convergence
means a relative `CT` change no larger than 1% under refinement or a matched
perturbation. Convergence and agreement with experiment are separate claims.

## Global directives

- Use 5400 RPM for every case in this study.
- Run adaptively: inspect each completed batch before choosing the next cases.
- Stop after every numbered phase and obtain Ryan's explicit approval before
  beginning any work, deployment, or submission from the next phase.
- Allow at most three active study jobs, each using 64 Julia/BLAS threads and
  32 GB.
- Retain at most four full ParaView histories from this new HPC study.
- Preserve unrelated local changes, remote changes, results, and histories.
- Never derive shedding from raw mesh cells. Construct the no-shedding body
  first, then trace on its rewound nodes and cells.

## How to resume

Read only:

1. this dashboard;
2. the current phase plan and the top **Current snapshot** section of its log;
3. the repository policies named by that phase.

Read older log entries or completed phase documents only to audit evidence or
diagnose a problem. At a phase boundary, copy only the decisions needed by the
next phase into that phase's log snapshot.

## Progress

Current phase: **Phase 2c — COMPLETE; decision recorded, awaiting Ryan's Phase 3 approval.** Phase 2b was approved by Ryan on 2026-07-24 (attribution: chordwise n_airfoil section-resolution discretization difference; oracle converged to 0.11% at 120 chordwise panels/section). **Phase 2c decision on the actual DJI mesh (30 spanwise; n_airfoil ∈ {81, 97, 121}): FINER MESH DESIRABLE.** The trustworthy signals support the attribution — the Neumann reference is flat/converged (Δ ≤ 0.30%/rung) and the outboard (|r/R|≥0.7) Dir-Neu gap shrinks monotonically (3.33% → 2.95% → 2.69%) with Dirichlet climbing toward Neumann — but it is **not converged** (outboard still 2.69% at n=121, > 1%). The full-∫Γ gap is non-monotonic (3.76% → 1.01% → 2.98%) solely because the **n_airfoil=97 capped/Dirichlet mesh is a single-rung inboard outlier** (rotor-hover 56_57 style); attribution NOT challenged. **Spanwise addendum (Ryan's 60_97 meshes, 2026-07-24): spanwise is CONVERGED** — the outlier-free indicators (Neumann solve, outboard integral) move ≤0.3% from 30→60 span; the 30_97 capped mesh was a defect, and the 60_97 Dirichlet value rejoins the family, making the corrected Dir–Neu ladder monotonic (3.76% → 3.30% → 2.98%). 30-spanwise is adequate; regenerate/replace 30_97 capped (or use 60_97). **Extended-chordwise addendum (Ryan's 45-span 145/185/201/249, 2026-07-24): the Dir–Neu gap CONVERGES monotonically toward the Neumann reference under chordwise refinement (the oracle mechanism, now on the real DJI mesh), but SLOWLY** — ∫Γ gap 3.76→3.30→2.98→**2.73%** at n=145 (Neumann flat, Dirichlet climbing ~0.1%/rung); still >1%. Only n=145 solved locally (dense Backslash; 17 GB RAM caps it); **185/201/249 submitted to HPC** (Ryan approved 2026-07-24; BYU `orc` Slurm job **12890990**, 64 thr / 32 GB / 12 h). On completion, pull the 6 raws back and re-run `PHASE2C_MODE=extend` locally to fill the pending rows. Phase 3 blocked until Ryan approves.

### [Phase 1 — Circulation Audit](phase_01_circulation_audit.md)

Log: [Phase 1 notes](../../logs/dji_convergence_20260722/phase_01_circulation_audit.md)

- [x] Setup and topology validation
- [x] Matched steady solves
- [x] Analysis and approval report

### [Phase 2 — TE Adequacy](phase_02_te_adequacy.md)

Log: [Phase 2 notes](../../logs/dji_convergence_20260722/phase_02_te_adequacy.md)

- [x] Geometry and operator audit
- [x] Circulation-observable sensitivity
- [x] Off-collocation and kernel-offset checks
- [x] Conditional thickness/local-refinement study
- [x] Adequacy decision and approval report

### [Phase 2b — Formulation/Topology Attribution](phase_02b_formulation_attribution.md)

Log: [Phase 2b notes](../../logs/dji_convergence_20260722/phase_02b_formulation_attribution.md)

- [x] Simple generated mesh oracle (3-rung ladder, matched semi-infinite wake)
- [x] Gate A field-circulation diagnostics (loop-circ = Γ_TE to 3e-5 → not extraction)
- [x] Conditional branch diagnostics (extraction branch closed; formulation channel excluded via Gate A0)
- [x] Conditional thickness/TE closure screen (gap grows with thickness → thinness does NOT explain DJI 5%)
- [x] Refinement-dimension separation (chordwise n_airfoil is the entire convergence lever)
- [ ] ~~Conditional DJI 2x2 bridge~~ (not needed; oracle isolated the mechanism)
- [x] Attribution decision and approval report (approved 2026-07-24)

### [Phase 2c — DJI Mesh Convergence](phase_02c_dji_mesh_convergence.md)

Log: [Phase 2c notes](../../logs/dji_convergence_20260722/phase_02c_dji_mesh_convergence.md)

- [x] Setup and topology validation (all 6 meshes)
- [x] Matched steady solves (3 rungs × 2 formulations)
- [x] Convergence table and decision report (FINER MESH DESIRABLE; n=97 capped outlier)

### [Phase 3 — HPC Driver](phase_03_hpc_driver.md)

Log: [Phase 3 notes](../../logs/dji_convergence_20260722/phase_03_hpc_driver.md)

- [ ] Parameterized driver
- [ ] Slurm launcher and deployment safeguards
- [ ] Local verification and approval report

### [Phase 4 — Green Gate](phase_04_green_gate.md)

Log: [Phase 4 notes](../../logs/dji_convergence_20260722/phase_04_green_gate.md)

- [ ] Existing regression tests
- [ ] Rotor frozen-wake oracle
- [ ] Green solver-route validation
- [ ] Accuracy classification and approval report

### [Phase 5 — CT Convergence](phase_05_ct_convergence.md)

Log: [Phase 5 notes](../../logs/dji_convergence_20260722/phase_05_ct_convergence.md)

- [ ] Steady preflights and cold baselines
- [ ] Plateau continuation
- [ ] Regularization, gradient, and Green perturbations
- [ ] Timestep refinement
- [ ] Final 57-series confirmations and approval report

### [Phase 6 — Final Report](phase_06_final_report.md)

Log: [Phase 6 notes](../../logs/dji_convergence_20260722/phase_06_final_report.md)

- [ ] Factor ranking and experimental comparison
- [ ] Retention cleanup and manifest
- [ ] Reviewed final report
