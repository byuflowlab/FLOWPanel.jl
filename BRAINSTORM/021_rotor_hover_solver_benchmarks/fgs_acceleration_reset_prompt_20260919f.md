# Reset prompt: FGS dagteam PROMOTED — wrap-up and follow-ons (2026-09-19f)

Copy everything below the line into a fresh agent session started in
`~/Dropbox/research/projects/FLOWPanel.jl`.

---

## Where this stands

BRAINSTORM 021 FGS acceleration is **DONE through promotion**. Read first:
`CLAUDE.md`, `agent_policies/WORKFLOW.md` + `TESTING.md` + `HPC.md` (before
corresponding work), then in BRAINSTORM/021:
`fgs_acceleration_status_20260919d.md` (the full arc: gate-1 closure,
plumbing, numerical gate, A/B verdict, PROMOTED section) and
`fgs_acceleration_provenance_20260919d.md` (v23 pins/jobs/decision rules).
The spec was `fgs_acceleration_recommendation_20260918.md`. The one-page
campaign story (what was tried, what won, why, accuracy certification,
evidence map) = `fgs_acceleration_summary_20260919.md`.

**Final verdict (2026-09-19, Ryan promoted):** dagteam split executor +
f32full + `numactl --interleave=0-3 --cpunodebind=0-3` @ j16/BLAS1 =
**4.475 s median R4 solve = 2.26× vs the v21 accepted colored@j16 baseline
(10.116 s)**; numerical gate certified (evaluator BC rel-L2 ≤ 1e-6, dagteam
tolerance 3.43e-7, 27 iters, bitwise repeatable). Side finding: interleave
placement alone gives colored@j16 7.243 s (free 1.40×) — adopted for the
fallback rung. Jobs 13773687 (numgate) + 13773689 (A/B), both COMPLETED
clean; evidence harvested to
`fgs_r4_followup_evidence_20260914/dagteam-{numgate-13773687,v23-13773689}/`.

## Promotion state (already executed — do NOT redo)

- Live FastMultipole `../FastMultipole` branch `flowpanel-20260817`
  fast-forwarded `c18e4b46` → `f4d6b671` (dagteam executor, gate tests,
  multi-system fill fix; additive). FLOWPanel live env loads it.
- FLOWPanel branch `fastmultipole` carries: `8aa5511` (FGSSolver
  `dagteam_precision` plumbing, forwarded only when sweep_order===:dagteam),
  `894ba2d` (v23 A/B driver `benchmark/fgs_r4_dagteam_ab.jl` + launchers +
  smoke test), and the champion config `benchmark/retained_r4_champion.toml`
  (calibrated tolerances for dagteam champion AND colored fallback in its
  header; placement documented there — it is part of the champion).
- Campaign infra: `/home/rander39/campaigns/p021-r4-dagteam-20260919-v23/`
  (worktrees, env, pins), tags `campaign/p021-r4-dagteam-{source,exec,fm}-
  20260919-v23`, orc refs `p021-r4-dagteam-v23` (FLOWPanel) and
  `p021-fgs-accel-20260918-v23` (FastMultipole).
- Gate-1 test `test/fgs_dagteam_gate1_test.jl` is 28/28; the old
  multi-system recovery oracle was replaced by a one-sweep fixed-point
  oracle — it MUST stay exactly 1 sweep (block-GS on the 1/r fixture
  amplifies ~5e5 per sweep; divergence there is fixture-class, not a bug).

## Remaining work (in rough priority order)

1. **Ryan-pending ledger** (all Ryan-gated): origin pushes (FLOWPanel
   `fastmultipole` incl. 8aa5511/894ba2d/champion-config commit;
   FastMultipole `flowpanel-20260817`@f4d6b671 + branch
   `p021-fgs-accel-20260918` + v21/v22/v23 campaign tags), **5 notebook
   entries owed for 021** (v21 colored A/B, diagnostics ladder, v22 chunked,
   NUMA, v23 dagteam gates+promotion), WeakKeyDict/warmstart one-line fix
   (pre-existing FLOWPanel failure from `7fbd68a`), 3 RECENT p018 runs
   awaiting `--include-recent --only` archive approval.
2. **Storage report**: an hpc-storage archive pass (~213 GB of finished p018
   runs) was running at session end — collect its final before/after +
   STALE/VERIFY counts and append the ledger line if the session died before
   it reported.
3. **Adopt the champion in consumers**: production launchers/scripts that
   run R4 FGS solves should load `retained_r4_champion.toml` and the
   documented numactl placement; the colored fallback should get interleave
   too. Check `benchmark/run_r4_diagnostics.slurm.sh`-style launchers and
   any simulate! pathways that construct FGSSolver at R4 defaults.
4. **Spec follow-ons (rank 4/5, separate tracks, need Ryan's direction):**
   joint leaf/MAC/P/inner retune around the NEW cost balance (staircase per
   point; report separately from same-iterate engineering — the old
   P8/MAC0.4/leaf100/inner3 was tuned for serial 29 GB/s sweeps), and the
   algorithm-changing track (safeguarded Anderson first). The 3.0 s
   non-sweep remainder now dominates (measured 4.475 total); further big
   wins come from remainder/iteration count, not bandwidth.
5. Housekeeping: local scratch env at `<scratchpad>/dagteam_env` is
   session-specific (gone after reset; rebuild via Pkg.develop live repos if
   needed — live FastMultipole now suffices, no worktree needed). The local
   worktree `/private/tmp/fastmultipole-p021-fgs-accel-20260918` can be
   retired once branches are pushed (it holds no unique state:
   f4d6b671 is merged; keep until origin push lands).

## Traps

- Never index per-thread scratch by threadid() (Julia 1.12 interactive pool).
- Spinning dagteam workers must never be alive across the FMM call
  (per-inner-sweep-block spawn/stop pattern in solve_dagteam.jl).
- dagteam vs lex is tolerance-equivalent NOT bitwise; per-run determinism IS
  bitwise at any thread count.
- Colored and dagteam tolerances are separately calibrated per environment —
  never carry one environment's staircase tolerance to another, and never
  mix iteration counts across sweep orders.
- `numactl` placement is part of the champion: without interleave the
  dagteam socket-membind control collapsed to 1.09×.
- Local runs ≤ 4 threads; `ssh orc` needs a live ControlMaster socket;
  Slurm CLI needs a login shell (`bash -lc` + `source /etc/profile`).
- Judge cluster runs by outputs (COMPLETED markers, TOMLs), never sacct.

## House rules (binding, unchanged)

HPC submission Ryan-gated; monitoring via `hpc-monitor`; storage via
`hpc-storage`; notebook writes Ryan-gated (offer, don't write); dated
status/provenance files in BRAINSTORM/021 per convention; campaign ceremony
per ~/.claude/CLAUDE.md for any new official runs.
