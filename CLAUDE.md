# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository. It is the routing hub: read the matching policy file below **before** starting the corresponding kind of work.

## What This Is

FLOWPanel.jl is a 3D panel method solver for low-speed (inviscid, incompressible) aerodynamics, written in Julia. It supports non-lifting bodies (source/doublet panels), lifting bodies (vortex-ring panels with rigid wake), and unsteady simulation with free wakes. Computation is accelerated via the `FastMultipole` package (FMM) or direct N-body evaluation.

## Required Policy Reads (by task)

- Diagnosing or modifying repository code → `agent_policies/WORKFLOW.md` (subsystem/source-file routing, repo realities, editing boundaries, failure modes)
- Selecting or running tests / validation → `agent_policies/TESTING.md` (verification matrix, example smoke tests)
- Monitors, simulation output, pressure or force recovery, replay → `agent_policies/MONITORS.md` (dependency contracts, monitor inventory, PressureLaplace constraints)
- Long simulations, HPC, or Slurm work → `agent_policies/HPC.md` (allocation, threading, assets, checkpoints, output retention, submission boundary)

## Quick Commands

```bash
# Run all tests
julia --project -e 'include("test/runtests.jl")'

# Run an example
julia --project examples/sphere.jl
```

For anything narrower, use the matrix in `agent_policies/TESTING.md`.

## Orientation

Source lives in `src/` as `FLOWPanel_<subsystem>.jl` files; the authoritative load order is the `for header_name in [...]` include loop in `src/FLOWPanel.jl`. Roughly: element kernels → N-body backends (`DirectBackend`, `FastMultipoleBackend`) → body types (`NonLiftingBody`, `RigidWakeBody`, parameterized `{E,N,TF,DBC}`) → solvers (`BackslashNeumann`, `BackslashDirichlet`, `KrylovSolver`, `FGSSolver`) → frames/kinematics → post-processing → wake → monitors → `simulate!` → warm-start and replay. Multi-body solving uses `solve!(bodies::Tuple, solvers::Tuple)`. Meshes are built with `GeometricTools` (aliased `gt` everywhere). Forward/reverse AD is supported through `ImplicitAD`. For which file owns what, see the routing table in `agent_policies/WORKFLOW.md`.

`simulate!(body, frames, maneuver!, Vinf, t_range)` time-marches: propagate kinematics → solve body+wake → wake convection/rollup → call each monitor in tuple order.

## Cross-Cutting Conventions

- Panel influence functions follow Hess & Smith / Katz & Plotkin sign conventions; GeometricTools defines normals with the right-hand rule (opposite K&P) — element kernels account for this with explicit sign flips. Sign handling is the most common source of regressions.
- `strength[:, j]` stores the strength of the j-th element type for all panels (ncells × N matrix).
- `Uinf` passed to `solve` is a 3×ncells matrix (freestream at each control point).
- The `fields` vector in each body tracks which solution fields have been computed (used by post-processing).

## Critical Invariant: RigidWakeBody Shedding

Compute shedding from the *constructed* body's cells, never the raw mesh: with `ensure_winding=true` (default) the constructor re-winds `cells` in place, so shedding computed from raw mesh cells attaches the wake at the wrong edges with **no error** — the body silently sheds almost no circulation (observed: rotor-hover CT collapsed ~3.6×, 0.0505→0.014). Build a `noshedding` body first, run `calc_shedding_from_seed` on *its* `.nodes`/`.cells`, then rebuild with the shedding. See the `RigidWakeBody` docstring in `src/FLOWPanel_liftingbody.jl` and `examples/rotor_hover_convergence.jl`.

## Response Preferences

When you respond to prompts, if you ever need to ask me questions where I decide between possible options (not including asking permissions to perform shell commands), and you suspect any of the option will be token-heavy, also include a brief estimate of the token count required as a fraction of my 5-hour limit.

When in plan mode, default to preparing a plan with enough context that a new agent with clean context doesn't have to waste tokens by reading other files. Save the plan to file, tell me where it is saved, and stop so I can run with a new agent with clear context. If you don't think this would save any tokens in the long run, ask me if I would like you to prepare a plan for you to implement without clearing context.
