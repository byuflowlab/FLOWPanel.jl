# Workflow Policy

Read this before diagnosing or modifying repository code. For test selection see
`agent_policies/TESTING.md`; for monitors and pressure/force recovery see
`agent_policies/MONITORS.md`; for HPC/Slurm see `agent_policies/HPC.md`.

## Task Routing

- Panel kernels, element conventions, and sign behavior: `src/FLOWPanel_elements.jl`
- FMM/direct backend dispatch and acceleration behavior: `src/FLOWPanel_fmm.jl`, `src/FLOWPanel_elements_fmm.jl`
- Body data structures and boundary-condition setup: `src/FLOWPanel_abstractbody.jl`, `src/FLOWPanel_nonliftingbody.jl`, `src/FLOWPanel_liftingbody.jl`
- Linear solves and coupled solve orchestration: `src/FLOWPanel_solver.jl`, `src/FLOWPanel_system.jl`
- Post-processing fields, forces, and pressure recovery: `src/FLOWPanel_postprocess.jl`, `src/FLOWPanel_simulate_monitors.jl`
- Unsteady stepping, wake evolution, and monitor execution order: `src/FLOWPanel_simulate.jl`, `src/FLOWPanel_wake.jl`, `src/FLOWPanel_simulate_monitors.jl`
- Off-body field probes: `src/FLOWPanel_simulate_monitors_fieldprobe.jl`
- Replay of saved VTK output through monitors: `src/FLOWPanel_replay.jl`
- Warm-starting simulations: `src/FLOWPanel_warmstart.jl`
- Run metadata (TOML manifests): `src/FLOWPanel_metadata.jl`
- Reference-frame and kinematic issues: `src/FLOWPanel_frames.jl`
- Lifting-line coupling: `src/FLOWPanel_liftingline.jl`, `src/liftingline/*`
- Module wiring / load order: the `for header_name in [...]` include loop in `src/FLOWPanel.jl` (~line 78) is the authoritative list
- Public docs and theory notes: `docs/`, especially `docs/pressure_poisson.md`, `docs/kernel_gradients.md`, and `docs/dirichlet_potential_theory.md`
- Reproducible usage examples: `examples/`

## Repo Realities

- The top level contains many exploratory scripts, debug notes, images, VTK outputs, and generated meshes. Do not treat every top-level file as maintained source.
- Untracked files are common. Before editing, inspect `git status --short`; do not clean up or revert unrelated user changes in a dirty worktree unless explicitly asked.
- `examples/data/` and many geometry assets are large reference inputs, not places for casual edits.
- `docs/src/generate_examples*.jl` mirrors example behavior into the documentation; if an example interface changes, docs generation may also need updates.

## Editing Boundaries

- Primary maintained code lives in `src/`, `test/`, `examples/`, and selected `docs/` files.
- Be conservative around:
  - `.CondaPkg/`, `.julia_depot/`, `.mplconfig/`
  - `*.pvd`, `*.vtu`, `*.vtp`, generated images, and ad hoc debug outputs
  - one-off top-level debug scripts unless the user asked about that specific workflow

## Common Failure Modes

- Sign conventions are easy to break. Normals come from `GeometricTools` (right-hand rule), while panel formulas follow Hess & Smith / Katz & Plotkin conventions with explicit sign flips in the element kernels.
- Small changes to control-point offsets, panel characteristic lengths, or normal orientation can cause large pressure or boundary-condition regressions.
- FMM and direct backends should stay behaviorally aligned; do not validate one path only when touching shared induced-field logic.
- Monitor ordering matters. `audit_monitors` enforces dependency contracts; changes to what a monitor reads or writes usually require tests in both simulate and postprocess suites. See `agent_policies/MONITORS.md`.
- `RigidWakeBody` shedding must come from the constructed (re-wound) cells, never the raw mesh — see the invariant in `CLAUDE.md` and the `RigidWakeBody` docstring in `src/FLOWPanel_liftingbody.jl`.
- `PressureLaplace` depends on runtime body compatibility, cached geometry, and gradient-mode selection. Regressions often show up as subtle pressure differences rather than crashes.
- Multi-body and matrix-free solver paths can fail differently from single-body direct solves; if you touch solve orchestration, test both narrow solver units and at least one higher-level simulation path.

## Working Style

- State which subsystem you are changing and which verification command covers it.
- Prefer the smallest reproducer that still exercises the bug.
- If behavior differs between a unit test and an example, trust the smaller reproducer first and then explain the gap.
