# Monitors Policy

Read this for any monitor, simulation-output, pressure-recovery, or force work.
Source of truth: `src/FLOWPanel_simulate_monitors.jl` (docstrings),
`src/FLOWPanel_simulate_monitors_fieldprobe.jl`, `src/FLOWPanel_replay.jl`;
tests in `test/runtests_unit_simulate.jl` and `test/runtests_unit_postprocess.jl`.
Theory (pressure Poisson formulation): `docs/pressure_poisson.md`.

## Dependency Contract

Monitors are callables `(systems, wakes, frames, uinf, i_step, dt) -> nothing`
invoked in tuple order each timestep by `simulate!`. They declare data contracts
via traits:

- `monitor_provides(m)` — tuple of symbols (`:P`, `:F`, `:sectional_F`) this monitor registers
- `monitor_requires(m)` — symbols that must be provided by an *earlier* monitor in the tuple
- `monitor_requires_body_hessian(m)` — `true` flips `body.needs_velocity_gradient[]` in `simulate!` so the FastMultipole Hessian is populated into `body.velocity_gradient`
- `monitor_requires_induced_vorticity(m)` — `true` makes the solver accumulate volumetric induced vorticity into `body.induced_vorticity`

`audit_monitors(monitors)` validates ordering at the start of `simulate!` and
throws `ArgumentError` on the first unmet dependency. Data flows through a
`MonitorContext` (`monitor_register!` / `monitor_field`), so provided fields are
per-body registered values, not just flags. If you change what a monitor reads
or writes, update its traits and add tests in both the simulate and postprocess
suites.

Note: `monitor_requires(::ForceMonitor)` is **instance-dependent** —
`source=:pressure` requires `:P`; `source=:context_force` requires `:F`.

## Monitor Inventory (provides / requires)

- `PressureBernoulli(rho; unsteady, correct_kuttacondition, clip, backend)` — provides `:P`. Steady or unsteady Bernoulli at control points; unsteady `∂φ/∂t` uses zero/BE/variable-step-BDF2 startup and subtracts the moving-control-point term. Only scalar-potential sources enter the φ history; vector-potential-only wake sources trigger a one-time warning.
- `PressureLaplace(bodies, rho; ...)` — provides `:P`; `monitor_requires_body_hessian` and `monitor_requires_induced_vorticity` are `true` only for the gradient/acceleration modes that need them. See constraints below.
- `ForceMonitor(nt, i_system; i_frame, normalization, correct_kuttacondition, select, source, verbose)` — provides `:F`; requires `:P` (`source=:pressure`, default) or `:F` (`source=:context_force`, e.g. after `DragPolarMonitor`). Writes `distributed_force`, integrates force/moment in `frames[i_frame]` coordinates, stores 3×nt histories in `.force`/`.moment`. Optional `select` predicate `cp_frame -> Bool` zeroes excluded panels (mask cached while `ncells` unchanged; body assumed rigid in the selected frame).
- `SpanwiseLoadingMonitor(nbins, i_system; components, span_axis, normalization, binning, select, ...)` — requires `:F`, provides `:sectional_F`. Bins panel forces along a span axis; `binning=:control_point` or `:span_overlap`.
- `DragPolarMonitor(nt, i_system, polar, chord; inflow_method, ...)` — requires `(:F, :sectional_F)`, **re-registers `:F`** as inviscid + panel-distributed profile drag. Follow with `ForceMonitor(...; source=:context_force)` to integrate totals. `inflow_method`: `:surface`, `:probe`, or `:both` (default: probe drives drag, surface recorded for the `inflow_angle_diff` comparison).
- `KuttaJoukowskiForce(body, nt, i_system; rho, backend, normalization, select)` — no contract symbols; independent cross-check `F = ρ Σ γ (V × Δs)` at edge midpoints via a `FastMultipole.ProbeSystem`, body-relative velocity, force-on-body sign matching `ForceMonitor`. Requires `ConstantDoublet`/`VortexRing` kernel.
- `SurfaceVorticityForce(body, nt, i_system; ...)` — provides `:F` (diagnostic). Reconstructs `κ = -n × ∇ₛμ` from the strength gradient and integrates `ρ (V_cp × κ) dS`; the minus sign follows FLOWPanel's stored strength convention (opposite much of the literature).
- `BoundCirculationMonitor(body, nt, i_system; i_frame, radial_dimension, R)` — memory-only rotor-section bound circulation histories (TE `upper − lower` strength convention, plus slice estimate).
- `CylindricalFieldProbeMonitor` (`src/FLOWPanel_simulate_monitors_fieldprobe.jl`) — off-body field probes; no contract symbols.

## Normalization Callables

For `ForceMonitor`-family, called as `normalization(CF, CM, systems, frames, uinf)`:
- `WingNormalization(rho, Sref, Lref)` — `0.5 ρ |U∞|² Sref` (moments `… Lref`); default
- `NoNormalization()` — pass-through dimensional
- `RotorNormalization(rho, D, i_frame)` — `ρ n² D⁴` / `ρ n² D⁵`, `n` (rev/s) = `frames[i_frame].ω / 2π`
- `RotorNormalization2(rho, D, i_frame)` — disk-area dynamic-pressure scale from `ω` (rad/s)

For `SpanwiseLoadingMonitor`: `NoSectionalNormalization()`,
`FreestreamSectionalNormalization(rho, Lref)`, `RotorSectionalNormalization(rho, R, i_frame)`.

## PressureLaplace Constraints

- Construct with the **actual body objects**, in the same order later passed to `simulate!` — it preallocates per-body operators, preconditioners, and Krylov CG workspaces. Body count is checked at call time; identity is not.
- **Never pass `dt` at construction**; `simulate!` supplies it at monitor-call time. (A deprecated `(rho, dt)` constructor exists — do not use it.)
- Solves a sparse panel-centered surface pressure Poisson equation (FV Laplacian from shared-edge weights `w_ij = ℓ_ij/d_ij`, gauge-fixed by pinning `reference_panel` to `reference_pressure`; CG via Krylov.jl with `JacobiPressurePreconditioner` by default; `IncompleteCholesky`/`AMG` variants are reserved, not implemented).
- `gradient_mode` must be selected explicitly for production (`:edge_difference`, `:corrected_hessian`, or `:surface_velocity`); omitting it warns and falls back to `:edge_difference`. **No formulation has yet passed the complete unsteady pitching-wing acceptance gate.** `:corrected_hessian` (alias `:raw_hessian`) is diagnostic-only — the exterior hypersingular surface limit is not implemented.
- `acceleration_form=:material_derivative` is the default; `:lamb_vector` is a deprecated diagnostic (its `lamb_vorticity` and `kappa_basis` options are for sensitivity studies; rotor-hover lamb CT moved ~2.5e-3 between `:quad`/`:tri` bases).
- Unsteady mode finite-differences the inertial surface trace with zero/BE/variable-step-BDF2 startup; steady configurations allocate no history arrays.
- Operator, preconditioner, and workspace are reused across steps by default; set `rebuild_every_step=true` for deforming geometry.

## Replay

`replay(path, run_name; monitors, reconstruct, recompute, steps, ...)` in
`src/FLOWPanel_replay.jl` loads saved VTK body/wake/frame state and runs
monitors for selected saved steps **without** solving, wake propagation,
shedding, or kinematics. Triangular VTU cells only. State not recoverable from
VTK requires a replay manifest or a `reconstruct` callback, otherwise it throws.
Tests: `test/runtests_unit_replay.jl`.
