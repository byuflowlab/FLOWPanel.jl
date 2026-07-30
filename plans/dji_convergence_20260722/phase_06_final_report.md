# Phase 6 — Final Report

## Purpose and required context

Produce the reviewed factor ranking, experimental comparison, and retention
manifest without launching an optional follow-on study.

Read the repository instructions, the dashboard, the top snapshot in the
[Phase 6 log](../../logs/dji_convergence_20260722/phase_06_final_report.md),
and:

- `agent_policies/MONITORS.md`
- `agent_policies/HPC.md`

Do not read earlier phase documents. The Phase 5 log must identify the accepted
cases, unresolved dimensions, artifact locations, and histories still required
for warm starts.

## Analysis

Rank matched relative changes in final two-revolution `CT` and circulation for:

- re-fit geometry;
- combined cap/formulation;
- mesh;
- timestep;
- wake core;
- target offset;
- mu-gradient method;
- Green reconstruction.

Treat a change below the larger of 1% or twice revolution-block variability as
unresolved noise.

For every final candidate, report:

```text
signed relative error = (CT - 0.072) / 0.072
absolute relative error = abs(CT - 0.072) / 0.072
```

State numerical convergence and experimental agreement separately. Do not
tune a numerical parameter solely to force agreement with the experiment.

## Retention

Preserve compact CSV/TOML summaries and logs for every case. Before stripping a
history:

- verify its completion marker and compact summary;
- confirm it is not needed for a warm start;
- delete only its `.pvd`, `.vtu`, `.vtp`, and associated ParaView directory;
- record the exact path, preserved artifacts, and reason in the log.

Finish with the two final 57-series histories plus no more than two
diagnostically important checkpoints, for a maximum of four full histories.
Do not remove pre-existing local results or unrelated cluster histories.

## Deliverables and exit gate

- Reviewed final report with factor ranking, convergence claims, experimental
  errors, and limitations.
- Compact results index and complete retention manifest.
- Updated Phase 6 log and completed dashboard checkboxes.

Present the report and retention manifest for explicit approval before any
optional cleanup or follow-on study outside this plan.

