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

**End use (Ryan, 2026-07-24):** the **Dirichlet formulation is the production
path** for modeling the DJI rotor. This study exists to verify that Dirichlet
converges to the right answer and to determine the mesh resolution and settings
required to do so. A later simulation with **Green reconstruction of the
wake-induced potential** (Phase 4 gate) is expected to close much of the
remaining validation gap vs experiment. The Neumann solves are the reference
oracle only — a trustworthy Dirichlet configuration is the primary deliverable.

## Global directives

- Use 5400 RPM for every case in this study.
- Run adaptively: inspect each completed batch before choosing the next cases.
- Stop after every numbered phase and obtain Ryan's explicit approval before
  beginning any work, deployment, or submission from the next phase.
- Allow at most **six** active study jobs (raised from three by Ryan on
  2026-07-25), each using 64 Julia/BLAS threads. Memory is per-case: 32 GB
  suffices for the steady solves; the Phase 2e unsteady cases on the 36,752-panel
  ct4 mesh need 64 GB (Green peaks ~30 GB).
- Retain at most four full ParaView histories from this new HPC study.
- Preserve unrelated local changes, remote changes, results, and histories.
- Never derive shedding from raw mesh cells. Construct the no-shedding body
  first, then trace on its rewound nodes and cells.

## Tools

- **`scripts/generate_dji9443_mesh.sh`** — agents may generate additional DJI
  9443 meshes from the command line without opening the OpenVSP/Gmsh GUIs
  (usage examples in the script header; `--help` for full options). Positional
  form: `SPANWISE CHORDWISE ROOT_CAP TIP_CAP [CAP_TESS]` with cap modes
  `none|flat|round` selectable independently (legacy combined modes `capped` =
  flat root + round tip, `uncapped`, `flat`, `round`). CHORDWISE must be odd.
  Writes `examples/data/dji9443_YYYYMMDD_<span>_<chord>_<caps>.msh`; existing
  files are protected unless `--force`. Runs on Ryan's Mac only (needs local
  OpenVSP at `/Users/ryan/OpenVSP-3.51.1-MacOS/vsp` and `gmsh` on PATH;
  overridable via `DJI9443_VSP_BIN`/`DJI9443_GMSH_BIN`/`DJI9443_VSP3`) —
  generate locally, then `scp` to the cluster if needed.

## How to resume

Read only:

1. this dashboard;
2. the current phase plan and the top **Current snapshot** section of its log;
3. the repository policies named by that phase.

Read older log entries or completed phase documents only to audit evidence or
diagnose a problem. At a phase boundary, copy only the decisions needed by the
next phase into that phase's log snapshot.

## Progress

Current phase: **Phase 2e — Unsteady Hover CT Convergence (opened 2026-07-25; all six required cases RUNNING on HPC as jobs 12894164/12894165/12894166/12894168/12894169/12894170 since 02:39-02:48).** Goal: a converged unsteady hover CT on the Phase 2d production mesh (`dji9443_20260725_45_185_capped_captess4.msh`, 5400 RPM) with a `PanelParticleWake` at `nwakerows=1`, run both **without** Green reconstruction (`VelocityThroughSources`) and **with** it (`GreenReconstruction`), converging **timestep (NT)** and **wake truncation depth**. Convergence criterion (Ryan): CT settles to a small amplitude with a mean that changes little over 5 revolutions. Phase 2e **front-runs elements of Phase 4 (Green gate) and Phase 5 (CT convergence)** at Ryan's direction; mitigation is `test/formulation_test.jl` preflight plus always reporting Green CT side-by-side with the matched velocity run. HPC: 64 threads, **64 GB**/job (the ct4 mesh is 36,752 panels, so Green peaks ~30 GB — 24 GB would OOM), ≤6 active study jobs. Cluster access is `ssh orc`. — Prior status: **Phase 2d COMPLETE (2026-07-25) — awaiting Ryan's approval before Phase 3.** Root cause of the Phase 2c tip artifact: the **round CapUMinTess=3 tip cap** of the 20260723 `capped` meshes (Dirichlet tip Γ fragile ±10% under that recipe; eliminated: extraction, conditioning, regularization — architecturally inert for Dirichlet, since the potential kernel is unregularized — wake attachment, first-order cap geometry, spanwise sampling). Fix implemented & F-verified: **round CapUMinTess=4 tip caps** (`generate_dji9443_mesh.sh <span> <chord> flat round 4`; ct4≡ct2≡flat ≤0.5%): six-rung ladder (97→249) has a strictly monotone tip Γ approaching Neumann from below (no overshoot), all Dir–Neu gap bands monotone-decreasing (full ∫Γ 2.90→2.11% vs the matched tip-capped Neumann reference — Ryan's root-open/tip-capped idea, which closed the geometry share of the tip gap), and Dir ∫Γ changes ≤0.22%/rung from n=121. **Production prescription: flat root + round-ct4 tip, 45 span, n_airfoil ≥145 (185 preferred); do NOT use the 20260723 round-ct3 `capped` meshes for Dirichlet solves.** Full report: `data/dji_convergence_20260722/phase_02d_tip_cap_discrepancy/phase_02d_report.md`. **Appendix G added 2026-07-25 (Ryan follow-up on flow tangency; local-only, no HPC):** a Morino solve never imposes tangency (σ = −U_ext·n; the solver enforces the interior potential), so it was measured for the first time. **Decisive control: through the identical evaluation path, the Neumann solution satisfies U·n=0 to rms 4.8e-7·ΩR while Dirichlet sits at 3.1e-1·ΩR — five orders worse.** That validates the measurement machinery (near-singular sliver pairs, self-pair offset conditioning, ½∇μ half-jump all reproduce an exactly-tangent solution to ~1e-6) and establishes the Dirichlet residual as a genuine property of the Morino solution. Bulk tangency ~1–3% of local blade speed at n=145, concentrated in tip-cap apex slivers (top 1% of panels = 93% of Σ(U·n)²; worst r/R=0.998, |U|=991 m/s vs 67 m/s tip speed); mid-span best; TE panels NOT a hot spot. ∇μ reconstruction **exonerated** (tri/tri_robust/quad identical to 4 digits; half-jump 6 orders down). Under chordwise-only refinement the residual **grows** (inboard rms ×2.75 over 97→249), driven mainly by **panel aspect ratio** (median inboard AR 7.0→17.8; residual ~linear in AR; within-AR-band growth only +29% vs +175% aggregate). **Dir–Neu velocity gap monotone-decreasing at all 7 stations on all 4 rungs, converging at least as fast as the ∫Γ gap** (−41% at r/R=0.5 vs −19% for ∫Γ over 97→185), though its absolute level is higher at the tip (10% of local |U|) — **velocity is the more demanding metric; a converged Γ does not certify a converged velocity field.** Dirichlet velocity self-converges ~first order. Thrust contamination (exactly ΔP = ½ρ(U·n)², reference-free) 0.84% of this solve's CT at n=145 → 3.24% at n=249, but only ~1e-4 in absolute CT (offending panels hold 0.009% of wetted area). **Production prescription unaffected**; actionable recommendation is **refine spanwise alongside chordwise** to control AR. **A kernel-offset "fix" (koff R·1e-10→R·1e-3, an apparent 80× improvement) was tested and REJECTED**: the same change destroys the Neumann control's tangency (4.97e-7→4.88e-2), so it is an inconsistency with the solve operator, not accuracy — retained in the report as G2b so it is not re-proposed. Method note: residuals were recomputed from archived `_panels.csv` strengths with **no re-solve** (O(n) memory), which is what put n=201/249 (dense G 12.8/19.6 GB) back within local reach; ordering gate 0.0 and faithfulness gate 2.4e-6 vs a real solve. — Prior status: Phase 2c COMPLETE. Phase 2b was approved by Ryan on 2026-07-24 (attribution: chordwise n_airfoil section-resolution discretization difference; oracle converged to 0.11% at 120 chordwise panels/section). **Phase 2c decision on the actual DJI mesh (30 spanwise; n_airfoil ∈ {81, 97, 121}): FINER MESH DESIRABLE.** The trustworthy signals support the attribution — the Neumann reference is flat/converged (Δ ≤ 0.30%/rung) and the outboard (|r/R|≥0.7) Dir-Neu gap shrinks monotonically (3.33% → 2.95% → 2.69%) with Dirichlet climbing toward Neumann — but it is **not converged** (outboard still 2.69% at n=121, > 1%). The full-∫Γ gap is non-monotonic (3.76% → 1.01% → 2.98%) solely because the **n_airfoil=97 capped/Dirichlet mesh is a single-rung inboard outlier** (rotor-hover 56_57 style); attribution NOT challenged. **Spanwise addendum (Ryan's 60_97 meshes, 2026-07-24): spanwise is CONVERGED** — the outlier-free indicators (Neumann solve, outboard integral) move ≤0.3% from 30→60 span; the 30_97 capped mesh was a defect, and the 60_97 Dirichlet value rejoins the family, making the corrected Dir–Neu ladder monotonic (3.76% → 3.30% → 2.98%). 30-spanwise is adequate; regenerate/replace 30_97 capped (or use 60_97). **Extended-chordwise addendum (Ryan's 45-span 145/185/201/249, 2026-07-24): the Dir–Neu gap CONVERGES monotonically toward the Neumann reference under chordwise refinement (the oracle mechanism, now on the real DJI mesh), but SLOWLY** — ∫Γ gap 3.76→3.30→2.98→**2.73%** at n=145 (Neumann flat, Dirichlet climbing ~0.1%/rung); still >1%. Extended ladder now COMPLETE to n=249 (145 local; 185/201/249 on HPC, BYU `orc` Slurm job 12892031 — after fixing a grad_mu crash on the finest mesh and a stale-cluster-src monitor bug). **FINAL (corrected after per-station tip inspection): attribution DIRECTION supported, but the DJI gap is NOT converged (~2.2% and still decreasing).** The trustworthy signal is the **inboard (root) gap**: smooth monotone convergence 4.05→2.16% with Dirichlet climbing toward the flat Neumann reference — the Phase 2b mechanism — still ~2% at the finest rung. The apparent **outboard "convergence to ~0.28% at 185/201" was a TIP-CAP ARTIFACT, not real convergence** (capped Dirichlet overshoots Neumann at r/R≳0.9, non-monotonic in n_airfoil, absent from uncapped Neumann). grad_mu ruled out (dji45_145 bit-identical :quad vs :tri). **Ryan states the tip caps use an identical recipe for all meshes** (so NOT two recipes) → the tip artifact is now the subject of **Phase 2d**, which investigates it with high confidence (open hypotheses: cap↔sharp-TE-corner interaction, Morino tip conditioning, fixed-regularization-vs-shrinking-panels, and which μ side carries it). Phase 3 blocked until Ryan approves. **Phase 2d reviewed & re-scoped (2026-07-24):** Ryan clarified the end goal — Dirichlet is the production formulation (Neumann is only the oracle), so the tip must be made trustworthy, not merely documented-around. Scope: **A + C are the required core** (A localizes the anomaly and decomposes μ; C tests H5, the one hypothesis with a general code-fix payoff), B kept (nearly free), **D and E conditional — run only if A+C leave the cause ambiguous**. Decision rule reweighted: outcomes (i) mesh-recipe fix and (ii) regularization fix are preferred; (iii) "interior stations only" is a last-resort fallback that must still name which meshes give trustworthy Dirichlet tips. Agents may generate control meshes themselves via `scripts/generate_dji9443_mesh.sh` (see Tools). **Success criterion (Ryan): fix-and-verify — Phase 2d is successful only when a converged Dirichlet solution and its settings are demonstrated (step F: fix implemented, targeted code changes in scope — a major code overhaul needs Ryan's permission first — ladder re-run to a physical/monotonic tip and a ≤1% or demonstrably-converging Dir–Neu gap); diagnosis alone does not close the phase.**

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
- [x] Convergence table and decision report (extended to n=249; attribution direction supported, NOT converged ~2.2%; tip 185/201 = artifact → Phase 2d)

### [Phase 2d — Tip-Cap Circulation Discrepancy](phase_02d_tip_cap_discrepancy.md)

Log: [Phase 2d notes](../../logs/dji_convergence_20260722/phase_02d_tip_cap_discrepancy.md)

- [x] A — instrumented solve (145c vs 185c): field-vs-extraction + μ decomposition + conditioning (core)
- [x] B — tip-corner geometry/aspect-ratio diff across {145,185,201,249} (cheap, keep)
- [x] C — regularization (kerneloffset/kernelcutoff) sensitivity sweep (core)
- [x] D — simple-case reproduction — skipped (not needed; E2 isolated the cause on the DJI meshes)
- [x] E — causal confirmation via generated control meshes (E2, Ryan-approved)
- [x] F — fix implemented (round-ct4 mesh recipe) + six-rung ladder re-run verified
- [x] Decision report + Phase 2c pointer update
- [x] G — **Appendix (Ryan, 2026-07-25, after the phase closed):** flow-tangency residual,
  its localization and cause, Dir–Neu velocity-field agreement, Cp/CT consequence

### [Phase 2e — Unsteady Hover CT Convergence](phase_02e_unsteady_ct.md)

Log: [Phase 2e notes](../../logs/dji_convergence_20260722/phase_02e_unsteady_ct.md)

Front-runs Phase 4/5 elements at Ryan's direction (2026-07-25).

- [x] Phase docs + dashboard entry
- [x] Driver edits (`RHPC_MESH` honored, auto TE seeds, `RHPC_FORMULATION`, metadata)
- [x] Slurm launcher `examples/run_dji9443_hover_ct_hpc.slurm.sh`
- [x] Local preflight (formulation tests + both-formulation smoke)
- [x] Deploy, md5 verification, submission (jobs 12894164/12894165/12894166)
- [x] Batch 1 launched (velocity/green at NT=36 d4; velocity NT=72 d4)
- [x] Batch 1b launched (green NT=72 d4; velocity/green NT=36 d6) — cap now 6
- [ ] Batch 2 (adaptive: truncation depth, matched Green confirmation)
- [ ] Decision report and approval

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
