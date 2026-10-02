# Reset prompt: 021 champion adoption + FGS re-runs (2026-09-19g, supersedes f)

Copy everything below the line into a fresh agent session started in
`~/Dropbox/research/projects/FLOWPanel.jl`.

---

## Where this stands

BRAINSTORM 021 FGS acceleration is CLOSED AT PROMOTION (Ryan, 2026-09-19).
Read first: `CLAUDE.md`, `agent_policies/WORKFLOW.md` + `TESTING.md` +
`HPC.md` (before corresponding work), then in BRAINSTORM/021:
`fgs_acceleration_summary_20260919.md` (one-page story: what was tried, what
won, why, accuracy certification, evidence map) and
`fgs_acceleration_status_20260919d.md` (full arc incl. the PROMOTED
section). Provenance of the deciding jobs =
`fgs_acceleration_provenance_20260919d.md`.

**Champion (promoted):** `sweep_order=:dagteam` + `dagteam_precision=:f32full`
+ `numactl --interleave=0-3 --cpunodebind=0-3` @ 16 Julia threads / BLAS 1 →
R4 prepared cold solve 4.475 s = 2.26× vs the v21 accepted colored@j16
baseline (10.116 s). Config + calibrated tolerances (dagteam 3.43e-7,
colored fallback 5.22e-7) = `benchmark/retained_r4_champion.toml`; the
placement is PART of the champion (~1.6× by itself) and lives in launchers,
not the library. Accuracy certified per solve by the independent evaluator
(BC rel-L2 ≤ 1e-6; jobs 13773687/13773689, 240 accepted trials).

**Code state (all committed, NOT pushed to origin):** live FastMultipole
`../FastMultipole` branch `flowpanel-20260817` @ `f4d6b671` (dagteam
executor + multi-system fill fix; the live FLOWPanel Manifest loads it;
promotion smoke passed). FLOWPanel `fastmultipole`: `8aa5511` plumbing,
`894ba2d` v23 A/B assets, `7dc7a17` champion config, `363e29a` wrap docs +
evidence. **Defaults are deliberately UNCHANGED** (lex + f64): f32full's
accuracy is certified at R4 specifically, placement can't be a code default,
and tolerances are calibrated per order/environment — the champion is
adopted per launcher, never silently.

## IMMEDIATE TASK: adopt the champion and re-run the FGS runs

Ryan's directive (2026-09-19): FGS runs must now be RE-RUN on the promoted
configuration. Work order:

1. **Inventory FGS consumers** (delegate to `code-scout`): every launcher,
   driver, and example that constructs `FGSSolver` (or selects an fgs config
   in the cold harness / simulate! pathways) — production rotor-hover
   chains, ground-effect (022) drivers, benchmark launchers
   (`run_r4_diagnostics.slurm.sh` etc.), anything on the standing-runs
   list. Table: script → mesh/rung → current solver settings → whether its
   results are cited anywhere (ledger/notebook).
2. **Adopt the champion per consumer**: pass
   `sweep_order=:dagteam, dagteam_precision=:f32full` and add the numactl
   placement to the launcher (single-socket zen3 form
   `--interleave=0-3 --cpunodebind=0-3`, 16 threads, BLAS 1 for the body
   solve; on other node types re-derive placement, don't copy blindly).
   **Recalibrate the FGS tolerance per mesh/rung by the standard staircase**
   — NEVER carry 3.43e-7 to a different mesh, rung, thread count, or BLAS.
   For non-R4 or non-zen3 cases where certification is untested, the safe
   first rung is dagteam + f64 (mathematically equivalent to lex, still
   ~1.6× with placement), with f32full promoted only after that case's own
   evaluator/staircase certifies it.
3. **Physics-null validation BEFORE production re-runs** (023+025
   precedent: CT nulls +0.023%/+0.0038% when solver internals changed):
   run a short backend-matched discriminator pair — previous solver config
   vs champion — on a standing rotor-hover case; require CT (and Cd) null
   within the established budget. Local smoke ≤4 threads first; HPC
   submission Ryan-gated. If the null fails at f32full, drop that consumer
   to dagteam+f64 and report.
4. **Re-run the affected FGS runs**: ASK RYAN which standing runs to
   re-execute (candidates: current 018/022/026/032 chains that spend
   body-solve time, any benchmark whose solve-time numbers are cited).
   Official re-runs = full campaign ceremony (annotated tags
   `campaign/<item>-<slug>-YYYYMMDD`, worktrees from tags — the v23
   worktrees at `/home/rander39/campaigns/p021-r4-dagteam-20260919-v23/`
   are reusable pins for FastMultipole/FLOWPanel if unchanged — Manifest
   pins, provenance BEFORE submission, outputs to the data root).
5. **Ledger (Ryan-gated, standing):** origin pushes (both repos + v21/v22/
   v23 campaign tags), 5 notebook entries owed for 021 (v21 colored,
   diagnostics ladder, v22 chunked, NUMA, v23 gates+promotion),
   WeakKeyDict/warmstart one-line fix (`_publish_block_gs_status!`, from
   `7fbd68a`), 3 RECENT p018 runs awaiting `--include-recent --only`
   archive approval. Also collect the hpc-storage archive-pass report
   (~213 GB of finished p018 runs, running at session end) and append its
   ledger line.

## Follow-on tracks (separate, need Ryan's direction before starting)

The 3.0 s non-sweep remainder now dominates the 4.475 s champion solve —
further large wins come from iterations/remainder, not bandwidth: (a) joint
leaf/MAC/P/inner retune around the new cost balance (staircase per point;
report separately from same-iterate engineering — P8/MAC0.4/leaf100/inner3
was tuned for serial 29 GB/s sweeps); (b) safeguarded Anderson acceleration
of the outer FGS map (spec §algorithm-changing, rank 1 there).

## Traps

- f32full perturbs the operator: the internal residual can flatter it — the
  independent evaluator (or a physics null) is the only accuracy authority.
- Tolerances are per order AND per environment (cluster dagteam 3.43e-7 vs
  local 8.9e-8 for the same config); never mix iteration counts or carry
  staircase tolerances across sweep orders/machines.
- Placement is load-bearing: socket-membind collapsed dagteam to 1.09×.
  "node 0" ≠ "socket 0" on NPS4 EPYC (socket 0 = nodes 0-3).
- Never index per-thread scratch by threadid() (Julia 1.12 interactive
  pool); spinning dagteam workers must never be alive across the FMM call.
- dagteam vs lex is tolerance-equivalent NOT bitwise; per-run determinism IS
  bitwise at any thread count.
- Gate-1's multi-system oracle must stay EXACTLY one sweep (block-GS on the
  1/r fixture amplifies ~5e5/sweep; divergence there is fixture-class).
- RigidWakeBody shedding: compute from the CONSTRUCTED body's cells
  (CLAUDE.md critical invariant) if touching rotor drivers.
- Local runs ≤ 4 threads; `ssh orc` needs a live ControlMaster socket; Slurm
  CLI needs a login shell (`bash -lc` + `source /etc/profile`); judge runs
  by outputs (COMPLETED markers/TOMLs), never sacct.

## House rules (binding, unchanged)

HPC submission Ryan-gated; monitoring via `hpc-monitor`; storage via
`hpc-storage` (cap 400 G across /home/rander39); notebook writes Ryan-gated
(offer, don't write); dated status/provenance files in BRAINSTORM/021; each
agent works in its own worktree for campaigns, never a shared live checkout
while jobs run.
