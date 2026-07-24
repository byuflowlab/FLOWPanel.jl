# Phase 3 Log — HPC Driver

Plan: [Phase 3 — HPC Driver](../../plans/dji_convergence_20260722/phase_03_hpc_driver.md)

Dashboard: [DJI convergence progress](../../plans/dji_convergence_20260722/dji_convergence_20260722.md)

## Current snapshot

Status: **Blocked on Phase 2b approval**

- No study driver or launcher has been implemented.
- No new HPC files have been deployed and no study jobs have been submitted.
- At the 2026-07-22 read-only check, the `orc` queue was empty, the remote
  checkout was dirty with prior formulation/panel-wake work, and the four new
  meshes were absent.
- Phase 3 implements and verifies deployment interfaces but does not submit
  production cases.
- Phase 2b has been inserted to determine why Dirichlet solves produce lower
  integrated circulation than Neumann solves. Do not begin Phase 3 until that
  handoff is complete or Ryan explicitly cancels Phase 2b.

Phase 2b handoff: _pending_

## Established environment

- Cluster alias: `orc`
- Remote checkout: `/home/rander39/projects/FLOWPanel.jl`

Non-interactive SSH may lack Julia and Slurm on `PATH`. Prior successful jobs
used a login shell and the Julia 1.11.7 fallback documented in
`plans/20260721_panel_wake_progress.md`; recheck rather than assume.

Relevant existing untracked work to inspect without overwriting:

- `examples/rotor_panel_wake_study.jl`
- `examples/run_dji9443_panel_wake_convergence.jl`
- associated Slurm launchers

Reuse ideas cautiously; prior force histories are not automatically valid for
this study.

## Working records

### Implementation and validation

| Date | Command or file | Purpose/result |
|---|---|---|

### Deployment manifest

| Local file | Remote file | Local checksum | Remote checksum | Result |
|---|---|---|---|---|

### Cluster checks

| Date | Queue state | Julia/Slurm environment | Retained study histories |
|---|---|---|---|

### Decisions and next-phase handoff

No decisions yet.

## Dated entries

Append exact local validation and deployment commands, changed files and why,
checksums, failures/root causes, and artifact paths here. Keep the snapshot
above current after every material action.
