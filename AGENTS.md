# AGENTS.md

Operational guidance for coding agents working in this repository. This file is intentionally complementary to `CLAUDE.md`: use `CLAUDE.md` for architecture and domain background, and use this file for task routing, verification, and repo-specific working habits.

## First Reads

- Read `CLAUDE.md` first for module load order, key types, monitor contracts, and solver terminology.
- Read `README.md` only for package-level positioning and public-facing examples.
- Before editing, inspect `git status --short`; this repository often has a dirty worktree with user debug scripts and generated artifacts.

## Task Routing

- Panel kernels, element conventions, and sign behavior: `src/FLOWPanel_elements.jl`
- FMM/direct backend dispatch and acceleration behavior: `src/FLOWPanel_fmm.jl`, `src/FLOWPanel_elements_fmm.jl`
- Body data structures and boundary-condition setup: `src/FLOWPanel_abstractbody.jl`, `src/FLOWPanel_nonliftingbody.jl`, `src/FLOWPanel_liftingbody.jl`
- Linear solves and coupled solve orchestration: `src/FLOWPanel_solver.jl`, `src/FLOWPanel_system.jl`
- Post-processing fields, forces, and pressure recovery: `src/FLOWPanel_postprocess.jl`, `src/FLOWPanel_simulate_monitors.jl`
- Unsteady stepping, wake evolution, and monitor execution order: `src/FLOWPanel_simulate.jl`, `src/FLOWPanel_wake.jl`, `src/FLOWPanel_simulate_monitors.jl`
- Reference-frame and kinematic issues: `src/FLOWPanel_frames.jl`
- Lifting-line coupling: `src/FLOWPanel_liftingline.jl`, `src/liftingline/*`
- Public docs and theory notes: `docs/`, especially `docs/pressure_poisson.md`, `docs/kernel_gradients.md`, and `docs/dirichlet_potential_theory.md`
- Reproducible usage examples: `examples/`

## Verification Matrix

- Broad regression: `julia --project -e 'include("test/runtests.jl")'`
- Solver changes: `julia --project -e 'include("test/runtests_unit_solver.jl")'`
- FMM or induced-velocity changes: `julia --project -e 'include("test/runtests_unit_fmm.jl")'`
- Kernel gradient or Hessian-sensitive changes: `julia --project -e 'include("test/runtests_unit_kernel_gradient.jl")'`
- Body assembly / geometry bookkeeping changes: `julia --project -e 'include("test/runtests_unit_body.jl")'`
- Lifting-body or wake changes: `julia --project -e 'include("test/runtests_unit_liftingbody.jl")'`, `julia --project -e 'include("test/runtests_unit_wake.jl")'`
- Simulation / monitor ordering / time-marching changes: `julia --project -e 'include("test/runtests_unit_simulate.jl")'`, `julia --project -e 'include("test/runtests_unit_postprocess.jl")'`
- Warm-start changes: `julia --project -e 'include("test/runtests_unit_warmstart.jl")'`
- FGS convergence-history changes: `julia --project -e 'include("test/runtests_unit_fgs_history.jl")'`
- Lifting-line changes: `julia --project -e 'include("test/runtests_liftingline.jl")'`
- Analytical consistency checks: `julia --project -e 'include("test/runtests_analytical.jl")'`

When a change is localized, run the narrowest relevant test first, then expand to the nearest higher-level suite.

## Examples As Integration Checks

- Use examples only after unit coverage is in place or when reproducing behavior that unit tests do not capture.
- Good smoke tests:
  - `julia --project examples/sphere.jl`
  - `julia --project examples/duct.jl`
  - `julia --project examples/sweptwing.jl`
  - `julia --project examples/rotor_hover_pressurelaplace.jl`
- Prefer examples that exercise the subsystem you changed instead of long-running showcase cases.

## Repo Realities

- The top level contains many exploratory scripts, debug notes, images, VTK outputs, and generated meshes. Do not treat every top-level file as maintained source.
- Untracked files are common here. Avoid cleanup unless the user explicitly asks for it.
- `examples/data/` and many geometry assets are large reference inputs, not places for casual edits.
- `docs/src/generate_examples*.jl` mirrors example behavior into the documentation; if an example interface changes, docs generation may also need updates.

## Editing Boundaries

- Primary maintained code lives in `src/`, `test/`, `examples/`, and selected `docs/` files.
- Be conservative around:
  - `.CondaPkg/`, `.julia_depot/`, `.mplconfig/`
  - `*.pvd`, `*.vtu`, `*.vtp`, generated images, and ad hoc debug outputs
  - one-off top-level debug scripts unless the user asked about that specific workflow
- Do not revert unrelated user changes in a dirty worktree.

## Common Failure Modes

- Sign conventions are easy to break. Normals come from `GeometricTools`, while panel formulas follow Hess & Smith / Katz & Plotkin conventions with explicit sign handling.
- Small changes to control-point offsets, panel characteristic lengths, or normal orientation can cause large pressure or boundary-condition regressions.
- FMM and direct backends should stay behaviorally aligned; do not validate one path only when touching shared induced-field logic.
- Monitor ordering matters. `audit_monitors` enforces dependency contracts; changes to what a monitor reads or writes usually require tests in both simulate and postprocess suites.
- `PressureLaplace` depends on runtime body compatibility, cached geometry signatures, and Hessian-backed velocity gradients. Regressions often show up as subtle pressure differences rather than obvious crashes.
- Multi-body and matrix-free solver paths can fail differently from single-body direct solves; if you touch solve orchestration, test both narrow solver units and at least one higher-level simulation path.

## Working Style

- State which subsystem you are changing and which verification command covers it.
- Prefer the smallest reproducer that still exercises the bug.
- If behavior differs between a unit test and an example, trust the smaller reproducer first and then explain the gap.
- Keep new operational guidance in this file; keep deep architecture notes in `CLAUDE.md` or dedicated docs.
