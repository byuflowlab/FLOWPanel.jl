# Suddenly Started Wing / Wagner Validation

## Summary

Create a maintained `examples/suddenly_started_wing.jl` validation case reproducing the executable VortexLattice setup: \(AR=100\), \(c=U_\infty=1\), \(\alpha=5^\circ\), \(t^*=0:1/8:7\), an open-tip NACA 0012 wing, Neumann vortex-ring panels, and a rolled-up panel wake. Compare \(C_L(t)/C_{L,\mathrm{steady}}\) against

\[
\Phi(2t^*)=1-0.165e^{-0.09t^*}-0.335e^{-0.6t^*}.
\]

The script's actual timestep is \(1/8\); its `1/16` comment is stale.

## Implementation

- Generate a cosine-spaced, closed-trailing-edge NACA 0012 contour, extrude it uniformly across the span, triangulate the surface, and omit both tip endcaps.
- Construct `RigidWakeBody{VortexRing,1,Float64,false}` with `watertight=false`. Follow the required two-stage construction: build a no-shedding body, derive trailing-edge shedding from its rewound cells, then rebuild the final body.
- Use a fixed frame and no-op maneuver. Because FLOWPanel solves at \(t=0\), initialize that sample with zero freestream and apply the full angled freestream for \(t>0\); omit \(t=0\) from comparisons. Store dimensional forces and normalize afterward by the fixed \(q_\infty S\).
- Use `PanelWake` with all wake rows retained, `include_final_filament=false`, `core_size=10^{-3}c`, induced wake convection enabled, and attached-wake displacement \(0.25U_\infty\Delta t\), matching VortexLattice's `eta=0.25`.
- Recover unsteady lift using `PressureBernoulli(unsteady=true)` followed by `ForceMonitor`; project force onto the wind-axis lift direction. Keep Kutta-Joukowski force optional as a circulatory-only diagnostic.
- For each mesh, run a separate identical steady body with a semi-infinite wake and normalize its transient history by that mesh's steady \(C_L\), matching the VortexLattice script.
- Expose keyword/environment controls for mesh size, timestep, solver/backend, VTK output, output path, and single-case versus convergence execution. No package-level public API changes are required.

## Refinement and Outputs

- Spanwise ladder at 21 airfoil points: 12, 24, 48, then 96 strips as needed.
- Chord/surface ladder at the accepted span resolution: 21, 41, then 81 airfoil points as needed.
- Declare successive histories converged when the maximum pointwise relative difference over \(t^*=0.125\ldots7\) is at most 2%. Also report relative \(L_2\) change, raw steady \(C_L\), Wagner RMS error, and maximum Wagner error.
- Confirm joint convergence by doubling the accepted span resolution at the accepted chord resolution; continue alternating the unresolved direction if this changes the curve by more than 2%.
- After spatial convergence, rerun the accepted grid at \(\Delta t^*=1/16\) and report timestep sensitivity without changing the primary \(1/8\) comparison.
- Use a direct/backslash coarse case as a reference cross-check, then fixed-tolerance Krylov/FMM settings for scalable refinement. Benchmark each rung; runs expected to exceed 20-30 local minutes move to a prepared Slurm launcher, which remains for user submission.
- Write reproducible per-case CSV/VTK data and an aggregate convergence CSV under `data/suddenly_started_wing/`. Produce:
  - a Wagner comparison plot containing every refinement, with the accepted FLOWPanel curve emphasized;
  - a convergence/error plot showing successive-grid and Wagner errors versus panel count.

## Verification

- Add lightweight tests for NACA geometry, open tips, continuous trailing-edge shedding, Neumann body type, timestep/time coordinates, sudden-start initialization, Wagner values, and wake capacity.
- Add a tiny two-step simulation smoke test checking finite post-start lift, monitor ordering, and no wake overflow.
- Run the new focused test plus the lifting-body, wake, simulation, and post-processing unit suites.
- Run the coarse case first, inspect force sign/time alignment and wake geometry, then execute the adaptive refinement ladder and the final half-timestep check.
- Preserve all existing unrelated working-tree changes; this work should be confined to the new validation example/test, its optional launcher, and generated validation outputs.

## Progress (2026-07-21)

- [x] Implement example geometry/formulation, corrected startup, wake, pressure/force recovery, steady normalization, environment controls, resume, CSV/VTK, and plots.
- [x] Add focused smoke/contract tests and run focused plus lifting/wake/simulate/postprocess validation.
- [x] Complete FMM interaction-backend span cases 12, 24, 48, and 96 plus the coarse direct cross-check. All body linear solves through 3,840 panels used dense Backslash.
- [x] Diagnose the coarse direct/FMM disagreement with a shortened case through \(t^*=2.25\), including direct-like FMM settings and split wake/system/pressure backends. The dominant trigger is FMM bound-body influence on wake nodes (`backend_system`), whose small numerical perturbations are amplified by the underresolved free wake; this is not an FMM linear-solve discrepancy or a Bernoulli-pressure discrepancy.
- [ ] Continue span refinement at 192 strips on HPC, then further if required, until both maximum absolute and relative-L2 changes are at most 0.02 or a formulation fault is established.
- [ ] Run airfoil refinement 21 to 41 to 81 at the accepted span.
- [ ] Double the accepted span for joint confirmation; alternate unresolved dimensions if needed.
- [ ] Run the accepted spatial grid at \(\Delta t^*=1/16\).
- [ ] Document the finite-thickness NACA versus flat VLM model difference, refresh final CSV/plots, and state the Wagner validation result.

## Current blockers and cautions

The 48 to 96 curve change is max absolute 0.24457 and relative \(L_2\) 0.071762. The coarse direct/FMM relative \(L_2\) difference is 0.27197. `CLsteady` remains highly mesh-dependent and nonphysical (6.98 even at 96 strips versus roughly 0.55 from thin-airfoil theory), so convergence must not be claimed.

The coarse histories agree to about \(10^{-8}\) through \(t^*=1.5\), then diverge after FMM tree subdivision becomes active. In the shortened diagnostic, standard FMM and `SSW_FMM_THETA=0` were bit-for-bit identical and differed from Direct by 0.05810 at \(t^*=2.25\), while `SSW_FMM_LEAF=10000` matched Direct bit-for-bit. Hybrid pressure backends followed the dynamics backend, excluding Bernoulli evaluation as the source. Split dynamics showed that FMM `backend_system` with Direct `backend_wake` reproduced nearly all of the difference, while Direct `backend_system` with FMM `backend_wake` stayed near Direct. This supports amplification of small tree partition/direct-summation perturbations by the underresolved wake; it does not by itself establish an FMM kernel defect. The exposed `SSW_FMM_ORDER`, `SSW_FMM_THETA`, and `SSW_FMM_LEAF` controls support follow-up sensitivity runs.

Long runs belong on HPC and agents must not submit them. Verify the site's Julia module and the 20,000-panel dense cutoff before user submission. See `plans/20260721_suddenly_started_wing_handoff.md` for detailed restart state.
