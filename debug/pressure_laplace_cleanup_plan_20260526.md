# PressureLaplace Cleanup, Mixable Modes, and Stage 2 Rerun

## Summary

Refactor `PressureLaplace` so `gradient_mode` controls whether the raw
induced-velocity Hessian is requested and used. `acceleration_form` remains the
RHS formulation choice. Remove dead scratch state, preserve the existing public
keywords, and rerun Stage 2 with all four combinations.

## Key Changes

- Treat `gradient_mode` as the derivative-source selector:
  - `:raw_hessian`: request/use `body.velocity_gradient` and add the rigid-body
    `[Omega]_x` contribution where needed.
  - `:surface_velocity`: do not request raw Hessian; reconstruct surface velocity
    gradients only when needed for material derivative.
- Make `acceleration_form=:material_derivative` use the selected gradient source:
  - raw mode: build `Du/Dt` from `body.velocity_gradient + [Omega]_x`.
  - surface mode: build `Du/Dt` from `compute_surface_velocity_gradient!`.
  - feed the resulting acceleration to `_pressure_rhs_from_acceleration!`.
- Make `acceleration_form=:lamb_vector` use mode-specific vorticity:
  - raw mode: compute vorticity as curl of the raw induced-velocity gradient.
  - surface mode: use the existing `body.induced_vorticity` field, populated
    through FMM/direct `extra_outputs=3`.
- Remove dead/obsolete code:
  - remove `tangential` from `PressureLaplace`, constructor allocation, docs, and
    checks.
  - keep `acceleration` because material derivative will actively use it.
  - remove the old edge-difference material-derivative RHS helper and replace its
    tests.

## Simulate/Replay Dependency Handling

- Update monitor traits to reflect the mode/form combination:
  - `monitor_requires_body_hessian(m::PressureLaplace) = m.gradient_mode == :raw_hessian`
  - `monitor_requires_induced_vorticity(m::PressureLaplace) = m.acceleration_form == :lamb_vector && m.gradient_mode == :surface_velocity`
- This makes `simulate!` use its existing `velocity_gradient=Tuple(requires_hessian(...))`
  path for raw-hessian variants.
- This makes `simulate!` request `extra_outputs=3` only for
  `lamb_vector + surface_velocity`, so `body.induced_vorticity` is available
  without computing the raw Hessian.
- Update replay auto-recompute with the same trait behavior:
  - raw-hessian variants request `:velocity_gradient`
  - surface Lamb requests `:induced_vorticity`
  - surface material derivative requests neither Hessian nor induced vorticity

## Tests

- Update `test/runtests_unit_postprocess.jl`:
  - all four `gradient_mode` x `acceleration_form` combinations construct and run.
  - raw-hessian variants require body Hessian.
  - surface-material variant requires no Hessian and no induced vorticity.
  - surface-Lamb variant requires induced vorticity but no Hessian.
  - material derivative tests verify gradient-driven acceleration feeds
    `_pressure_rhs_from_acceleration!`.
  - Lamb-vector tests verify raw mode uses curl of `body.velocity_gradient`, while
    surface mode uses `body.induced_vorticity`.
- Update simulate/replay dependency tests:
  - raw-hessian PressureLaplace causes `velocity_gradient=true`.
  - surface-Lamb PressureLaplace causes `extra_outputs=3`.
  - raw-Lamb PressureLaplace does not request `extra_outputs=3`.
- Run:
  - `julia --project=. -e 'include("test/runtests_unit_postprocess.jl")'`
  - `julia --project=. -e 'include("test/runtests_unit_simulate.jl")'`
  - if replay dependency tests exist/are added:
    `julia --project=. -e 'include("test/runtests_unit_replay.jl")'`

## Stage 2 Rerun

Rerun Stage 2 after tests pass:

```bash
env RUN_STAGE2_LAPLACE=true RUN_STAGE1_LOG=true RUN_BERNOULLI=false RUN_LAPLACE_MD=false RUN_LAPLACE_LAMB=false RUN_KJ=true NSTEPS_REPLAY=4 STAGE1_LOG_PATH=debug/logs/rotor_hover_stage2_gradient_active_20260526_final4.csv julia --project=. examples/rotor_hover_pressure_comparison_replay.jl
```

Append to `debug/debug_rotor_hover_20260526.md`:

- command and CSV path
- CT for all four variants
- pressure RMS over `q_rotor`
- panelwise cancellation values
- interpretation of whether the split follows `gradient_mode`, `acceleration_form`,
  or both

## Assumptions

- The repo field corresponding to "body vorticity" is `body.induced_vorticity`.
- Public `PressureLaplace` keywords remain unchanged.
- `:surface_velocity + :lamb_vector` intentionally uses `body.induced_vorticity`
  from `extra_outputs=3`, not curl of the reconstructed surface gradient.
