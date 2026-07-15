# Testing Policy

Read this when planning or running validation. All commands run from the repo top level.

## Commands

```bash
# All tests from shell
julia --project -e 'include("test/runtests.jl")'

# A specific test file
julia --project -e 'include("test/runtests_unit_solver.jl")'
```

From the Julia REPL (with `--project`): `] test` for the full suite, or
`include("test/<file>.jl")` for one file.

## Verification Matrix

Pick the narrowest test that covers the change first, then expand to the nearest
higher-level suite.

- Broad regression: `test/runtests.jl`
- Solver changes: `test/runtests_unit_solver.jl`
- FMM or induced-velocity changes: `test/runtests_unit_fmm.jl`
- Kernel gradient or Hessian-sensitive changes: `test/runtests_unit_kernel_gradient.jl`
- Body assembly / geometry bookkeeping changes: `test/runtests_unit_body.jl`
- Lifting-body or wake changes: `test/runtests_unit_liftingbody.jl`, `test/runtests_unit_wake.jl`
- Simulation / monitor ordering / time-marching changes: `test/runtests_unit_simulate.jl`, `test/runtests_unit_postprocess.jl`
- Replay changes: `test/runtests_unit_replay.jl`
- Warm-start changes: `test/runtests_unit_warmstart.jl`
- FGS convergence-history changes: `test/runtests_unit_fgs_history.jl`
- Lifting-line changes: `test/runtests_liftingline.jl`
- Pitching-wing experimental data handling: `test/runtests_unit_pitching_wing_exp.jl`
- Analytical consistency checks: `test/runtests_analytical.jl`
- Example-level pitching-wing regressions (long-running; run only when the change targets those examples): `test/runtests_example_pitching_wing.jl`, `test/runtests_example_pitching_wing_convergence.jl`, `test/runtests_example_pitching_wing_pressure_comparison.jl`

`test/dirichlet_potential_test.jl` is an ad hoc diagnostic script, not part of the
matrix. `test/test_helpers.jl` is shared setup included by the suites.

## Examples As Integration Checks

Use examples only after unit coverage is in place or when reproducing behavior
that unit tests do not capture. Prefer examples that exercise the subsystem you
changed instead of long-running showcase cases. Good smoke tests:

```bash
julia --project examples/sphere.jl
julia --project examples/duct.jl
julia --project examples/sweptwing.jl
julia --project examples/rotor_hover_pressurelaplace.jl
```
