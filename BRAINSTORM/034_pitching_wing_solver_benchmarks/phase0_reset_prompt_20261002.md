# 034 Phase 0 reset prompt — 2026-10-02

You are starting **Phase 0 (harness adaptation + availability smoke)** of BRAINSTORM
item 034: replicating 021's solver-benchmark campaign on the pitching-wing case.
The item is STAGED — nothing has been run, no code has been changed.

## Read first (in this order, nothing else up front)

1. `BRAINSTORM/034_pitching_wing_solver_benchmarks.md` — control doc: RESET BRIEF,
   standing rulings 1–11, Phase gates table. Rulings are binding.
2. `BRAINSTORM/034_pitching_wing_solver_benchmarks/decision_rules.md` — binding
   metrics/protocols (cold protocol, BC certification + its Dirichlet-only caveat,
   iteration counting, CSV schema).
3. `agent_policies/WORKFLOW.md` and `agent_policies/TESTING.md` before editing or
   testing repo code (per CLAUDE.md).
Do NOT read 021's full item file inline; `decision_rules.md` here already carries what
Phase 0 needs. Use `brainstorm-scout` if you need more 021 context.

## Phase 0 deliverables (gate = clear-context subagent review; Ryan approval separate)

1. **Solver injection.** Add a `solver` kwarg (or solver factory, since the static
   polar builds its own body) to `prepare_pitching_wing` / the static-polar path in
   `examples/pitching_wing.jl`, threaded to the two hardcoded sites:
   `solver = pnl.Backslash(wing)` at `:1008` (unsteady; consumed as
   `body_solvers=(solver,)` at `:1160`/`:1166`) and `body_solvers=pnl.Backslash(body)`
   at `:867` (static polar). Default must preserve `Backslash` bit-identically —
   verify with the existing smoke test `test/runtests_example_pitching_wing.jl`
   (config at `:28-47`: `n_cycles=0.01, n_span=1, n_airfoil=21, n_endcap=5`,
   DirectBackend).
2. **Mesh ladder.** Define 3–4 rungs spanning ~2k → ~30k panels via
   (`n_span`, `n_airfoil`, `n_endcap`); cells =
   `2*n_sec*n_span + 4*(n_chord-2)*(n_endcap-1)` (`examples/pitching_wing.jl:93-95`,
   defaults 13/161/9 → 6688). Document the chosen knob triples and computed counts in
   the ledger. The ladder freezes at the END of Phase 1, so mark it provisional.
3. **LHS-rebuild measurement (headline question).** Instrument a short `simulate!`
   run to MEASURE whether G is rebuilt/refactored per step under the pitching
   maneuver (`pitching_maneuver!` at `:971` writes only the frame angular rate; the
   newest `PanelWake` row couples to the unknowns, so constancy is a hypothesis, not
   a given). Deliver: (a) what the code actually does today per step
   (rebuild/refactor/reuse, with file:line evidence), (b) measured per-step cost of
   whatever rebuild occurs on the smallest rung. Related prior art:
   `BRAINSTORM/021_rotor_hover_solver_benchmarks/rigid_motion_tree_reuse_item.md`
   (`transform_tree!` — COMPLETE, tree-level reuse under rigid motion).
4. **Dirichlet confirmation.** Confirm the wing body under the pinned default
   formulation (`VelocityThroughSources`) is Dirichlet so `bc_error!`
   (`benchmark/common.jl`) applies as-is. If Neumann, STOP and flag — the metric
   section of `decision_rules.md` must be amended via the control doc's decision log.
5. **Availability smoke (021 W1–W6 analog).** All four arms run to completion on the
   smallest rung with a certified BC check where applicable: `backslash`,
   `krylov_gmres` (`pnl.KrylovSolver`, `src/FLOWPanel_solver.jl:1013`),
   `krylov_ilu_nfcache` (021 production iterative winner; 033's threaded repack),
   `fgs` (`pnl.FGSSolver`, `:1547`; 021/033 champion config dagteam+backoff as the
   seed). Usage patterns: `examples/sweptwing_solverbenchmark.jl:114`,
   `examples/suddenly_started_wing.jl:513`. Harness pieces go in flat `benchmark/`
   next to the 021 drivers; REUSE `benchmark/common.jl` (banner, `bc_error!`,
   `solver_state_bytes`, CSV validation) — do not rewrite it.

## Hard constraints

- Local runs ≤ 4 threads, smoke only — nothing local is publishable.
- Wake pinned `:panel`; formulation pinned `VelocityThroughSources`; no multi-body;
  no particle wake.
- Judge from harness-written CSVs, never stdout; banner beside every CSV.
- No HPC submissions in Phase 0 (everything here is smoke-scale). Any run > 20 min
  goes to HPC under campaign worktree rules — that starts in Phase 1 at the earliest.
- Write results to `BRAINSTORM/034_pitching_wing_solver_benchmarks/ledger.md` (dated,
  append-only), narrative to `log.md` (newest first), and create
  `phase_00_harness.md` in the item dir for the phase record. Update the control
  doc's `## Current status` and RESET BRIEF in place when Phase 0 state changes.
- Notebook entries, INDEX checkbox ticks, commits, and pushes are Ryan-gated — offer,
  don't do.

## Exit

Phase 0 ends with: solver kwarg merged locally + default-path smoke green, ladder
table in the ledger, the LHS-rebuild question answered with a measurement, Dirichlet
confirmed, 4/4 arms smoking on rung 1, and a clear-context subagent review of the
phase. Then stop and report to Ryan before Phase 1.
