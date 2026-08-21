# Item 020 — Theory Remediation and Phase-2 Reopening

**Execution record (2026-08-12): COMPLETE through Phase 2.** The corrected
screen crossed the Euler control's step-237 death healthy, then failed on an
explicit non-finite-ratio guard after a tail runaway at step 243. Literal
rescoring leaves “threshold moves” failed and “fidelity collapse remains” not
confirmed by its registered Boolean. See `phase_02_evidence_pack.md`; stop at
the revised gate before predictive fixtures or Phase 3.

## Summary

- Do not approve the current Phase-2 gate or begin closure implementation.
- Preserve Leg 1 and Leg 3 as supporting evidence, but invalidate the mechanistic conclusion from Stage B: the landed `euler_exp` in FLOWVPM freezes the stretching vector, not the velocity gradient. Under aligned strain its amplification approaches $5/3$ instead of the claimed physical $e^{2\Delta tZ}$.
- Treat the fixed-filter closure in `phase_01_theory.md` as a candidate model, not a uniquely derived or proven-unconditionally-stable closure.
- Target a predictive model: stability is required, but constants may not be chosen merely to hold $\sigma=\sigma_0$. A clean negative result is acceptable.

## Theory Corrections

- Distinguish molecular viscosity $\nu$, modeled SGS viscosity $\nu_{\mathrm{sgs}}$, and total $\nu_{\mathrm{eff}}=\nu+\nu_{\mathrm{sgs}}$; remove the present double-counting ambiguity.
- Scope the diffusion derivation to a locally uniform coefficient, isolated Gaussian blob, and localized/no-transfer approximation. Global vorticity conservation permits strength redistribution and does not uniquely require every $\dot{\Gamma}_p=0$.
- Replace the unconditional stability claim with a conditional result for the scalar Kelvin-filament model. The relation $Z\sim\Gamma_v/(2\pi\sigma^2)$ and its cap are modeling assumptions, especially because the registered $\Gamma_{\mathrm{implied}}$ criterion failed.
- Correct the conservation claims: fixed-strength core spreading preserves circulation and linear impulse, but angular impulse changes through its $\sigma^2\Gamma$ contribution. Replace T4 accordingly.
- Define the strain invariant explicitly. Compare aligned stretching, $s_+=\max(\hat\Gamma^TD\hat\Gamma,0)$, against $\sqrt{2D:D}$; $\kappa=1/\sqrt5$ is valid only for the former normalization and is a cutoff calibration, not a literature-derived physical constant.
- Relabel the reported 1.9–3.3-decade gap as a scenario estimate until the $r_c/c=0.01$–$0.05$ bracket and vortex-profile coefficient are independently supported.
- Correct the Squire/Bhagwat–Leishman discussion: Squire supports eddy viscosity proportional to circulation, while the empirical coefficient has a broad model-dependent range. Squire's original result does not justify $\kappa$.

## Phase-2R Work

- Replace the prototype with a frozen-gradient geometric update for the live $f=0$ formulation:
  - Compute $v=\exp(L\Delta t)\Gamma_n$ and $r=\lVert v\rVert/\lVert\Gamma_n\rVert$.
  - Update $\Gamma_{n+1}=v\,r^{-3/5}$ and $\sigma_{n+1}=\sigma_n r^{-1/5}$.
  - Apply the existing SFS forcing as a separately identified first-order split.
  - Integrate $y=\sigma^2$ and diffusion together:

    $$
    y_{n+1}=e^{-2Z\Delta t}y_n+
    \frac{\nu_{\mathrm{eff}}}{Z}\left(1-e^{-2Z\Delta t}\right),
    $$

    using the $Z\to0$ limit $y_{n+1}=y_n+2\nu_{\mathrm{eff}}\Delta t$.
- Verify aligned and general constant-gradient cases against matrix-exponential/reference integration, positivity for large $\Delta tZ$, first-order convergence with time-varying gradients, and off-state bit identity.
- Rerun only the decisive $\sigma/R=0.02$ Stage-B screen with the corrected integrator and otherwise identical settings. Retain the existing Euler run as the control.
- Rescore Stage B from survival, $\sigma/\sigma_{\mathrm{shed}}$, $M$, $|\Gamma|/\sigma^2$, and velocity histories. Update `phase_02_evidence_pack.md` before reconsidering the gate.

## Predictive Closure Gate and Tests

- Evaluate candidate strain estimators first on prescribed axial strain, pure rotation, simple shear, Lamb–Oseen/Squire core growth, and a strained Burgers vortex. Reject any estimator that produces material artificial spreading in an isolated coherent vortex or needs case-specific retuning.
- Fix constants only from these independent fixtures. Require analytic cases within 5% and literature-envelope cases within 15%; do not use hover stability or $C_T$ to fit them.
- If no fixed-constant candidate passes, close 020 with a negative predictive verdict and route a robust $\sigma_0$ cutoff to a separately labeled numerical model.
- If a candidate passes, implement it as a default-off `ViscousScheme` with per-particle $\nu_{\mathrm{sgs}}$, explicit strain-estimator metadata, and the corrected integrator as its required companion.
- Phase-3 tests must cover the corrected frozen-gradient map, the exact discrete fixed point, laminar recovery, circulation and linear impulse, the expected angular-impulse change, isolated-blob dissipation, off-state identity, ignition fixtures, and timestep convergence. Do not require monotonic whole-field enstrophy for arbitrary overlapping particles with spatially varying diffusivity without first proving that property.
- Only after those tests pass should Phase 4 repeat the 019 regime subset and validate fixed constants on DJI9443 hover plus a second flow family.

## Assumptions

- Phase 2 is reopened before any closure implementation.
- The closure objective is predictive modeling rather than a stability-only cutoff.
- Existing Phase-2 Leg 1 and Leg 3 data remain supporting evidence but are rescored under the corrected theory.
- The existing dirty worktrees in FLOWPanel and FLOWVPM are preserved; remediation must not overwrite unrelated user changes.
