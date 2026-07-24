# Phase 5 Log — CT Convergence

Plan: [Phase 5 — CT Convergence](../../plans/dji_convergence_20260722/phase_05_ct_convergence.md)

Dashboard: [DJI convergence progress](../../plans/dji_convergence_20260722/dji_convergence_20260722.md)

## Current snapshot

Status: **Blocked on Phase 4 approval**

- No study HPC cases exist and no new job has been submitted.
- No ParaView history from this study has been retained or deleted.
- Autonomous HPC submission is authorized within the explicitly approved
  Phase 5 scope, subject to the queue and retention limits.
- The first batch cannot be selected until Phase 1's old-control decision,
  Phase 2's geometry decision, Phase 3's verified interface, and Phase 4's
  Green classification are summarized below.

Required handoff:

- Accepted geometry/TE remedy: _pending_
- Old 40-series HPC control required: _pending_
- Verified driver/launcher and deployment state: _pending_
- Green classification and accepted route: _pending_

## Established numerical context

Prior steady mu-gradient histories:

| Method | Old 40-series CT | Old 56/57-series CT |
|---|---:|---:|
| default quad-growing | `0.05048051147416609` | `0.051489958562886455` |
| triangular | `0.05038462967921074` | `0.05128457315128889` |
| robust triangular | `0.050291047437732786` | `0.051158092178601434` |
| quad, growth disabled | `0.0970424573` | `0.0476993392` |

Therefore compare normal quad-growing with triangular first. Do not include
quad-without-growth in the routine sweep.

The existing particle-wake history at
`data/rotor_hover_pressure_comparison/` uses an old capped Dirichlet mesh,
`NT=36`, 6000 RPM, and wake core/target offset `1e-3`; it settles near
`CT ~= 0.062`. It provides scale and settings context only, not a matched
5400 RPM re-fit comparison.

Prior finite panel-wake in-run Bernoulli/Kutta-Joukowski values had known
force-recovery inconsistencies. Treat circulation and formulation residuals as
essential cross-checks.

## Active job ledger

No study jobs yet.

| Job ID | Case | Mesh/formulation | Key settings | Restart | Submitted | State | Output/error logs |
|---|---|---|---|---|---|---|---|

## Completed case ledger

No study cases yet.

| Case | Result | CT mean/variability | Circulation | Wall/RSS | Full VTK? | Decision/next action |
|---|---|---|---|---|---|---|

## ParaView retention ledger

Full histories retained: **0 / 4**

| Action/date | Exact path | Summary retained | Warm-start check | Reason |
|---|---|---|---|---|

## Dated entries

Immediately record exact submission commands, settings, job IDs/states,
elapsed time, peak memory, restart provenance, logs, numerical results,
failures/root causes, and every retained or removed history. Keep the snapshot
and ledgers current after each batch.
