# Reset prompt: 021 FGS scalability diagnostic (staged) + resume-harvest tail (2026-09-22, supersedes champion_adoption_reset_prompt_20260921.md)

Copy everything below the line into a fresh agent session started in
`~/Dropbox/research/projects/FLOWPanel.jl`.

---

## Where this stands

BRAINSTORM 021 is CLOSED AT PROMOTION (champion = dagteam + f32full +
`numactl --interleave=0-3 --cpunodebind=0-3` @ j16/BLAS 1; story
`fgs_acceleration_summary_20260919.md`). Read first: `CLAUDE.md`, then
`agent_policies/WORKFLOW.md`/`TESTING.md`/`HPC.md` before corresponding
work. Prior context: `champion_adoption_reset_prompt_20260921.md` (its
traps still bind) and the two provenance files with 2026-09-21 addenda.

The 2026-09-21/22 session verified and harvested the warm-start resume
reruns (13829231 R4 thread-scaling, 13829232 R1–R2) and delivered Ryan's
priority table. Ryan then directed the next investigation and staged its
plan (see NEXT TASK).

## Harvest results (2026-09-21/22, certified rows)

### R4 thread scaling — FINAL (median accepted solve time, s)

| j | FGS-dagteam | FGS-colored | krylov_ilu (certified, budget-0 knobs) | krylov_ilu_nfcache (budget-500 knobs) |
|---|---|---|---|---|
| 1 | 33.07 | 35.58 | FAILED | FAILED |
| 8 | 5.97 | 9.07 | 325.7 | **4.25** |
| 16 | 4.50 | 7.26 | 158.5 | **3.56** |
| 32 | 4.42 | 6.91 | 83.8 | **2.73** |
| 64 | 6.42 | 12.30 | 40.3 | **2.41** |

- Source: `<data root>/thread-scaling-j<J>-13777133/ilu/{tune_phase2,phase2}.csv`,
  data root `/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/`.
  Budget-0 rows verified present at j=8/16/32/64 (P=15, MAC=0.55,
  leaf=32/32/32/21; j=8 descent overran its 10 h cap, `tune_timed_out=true`,
  which is why that arm ran ~16 h — the cap only stops new candidates).
- **Headline: krylov_ilu_nfcache beats the dagteam champion at every
  measured j ≥ 8 at R4 and is the only config still scaling at 64.**
  Caveats for any comparison: nfcache t_solve EXCLUDES its one-time
  near-field cache build (~94 s @ j8, recorded in `nfcache_build_time`) and
  carries ~8.5 GB cache + ~9.8 GB solver state; dagteam is far lighter.
  This overturns the earlier "krylov_ilu never wins anywhere" premise
  (that was based on superseded budget-500 warm tuner times).
- **j=1 ILU FAILED again despite the `946cec2` seed fix**: the budget-0
  descent at j1 recorded its only candidate as t=Inf/success=false
  (consistent with a per-candidate timeout — one uncached j1 solve is
  enormous), so `rotor_hover_solver_phase2.jl` hard-errored at measure.
  R4 ILU@j1 stays missing unless rerun with a longer candidate cap
  (Ryan-gated; arguably academic).
- Diagnostic note: task logs are OUTPUT-BUFFERED — logs freeze for ~10 h
  mid-descent while the process computes at full CPU. Judge liveness by
  CPU/ps or output files, never by log mtime.

### R1–R2 resume merge (13829232)

- 5 of 6 arms harvested clean (R1 j1/j2, R2 j2/j4/j8 merged with
  first-pass rows; per-run `phase2/phase2.csv`; duplicate-row checks
  clean). Cross-rung story: `backslash_ldiv` untouchable at small rungs;
  FGS plateaus by j8–j16; krylov_ilu_nfcache wins R1 (0.164 s @ j64), R2
  from j8 up, and R4 from j8 up. R2 j64 column remains the known
  contention pathology (fgs ~10.7 s) — discard from scaling reads.
- **OWED (small): R2-j1 (13829232_7, run dir
  `r12-champion-r2-j1-13778533/`) was still running at handoff**
  (p2tune done 09-21 16:15, p2 pending; 48 h wall, plenty left). Check it
  (via `hpc-monitor`), harvest its `phase2/phase2.csv` (via `harvester`),
  and fill the R2 j1 column.
- **OWED (small): R1 j32 came back empty in the merge sweep** — recheck
  `r12-champion-r1-j32-13778533/phase2/phase2.csv` (first-pass arm; was it
  parsed / does it exist?) and fill or explain the gap.

## NEXT TASK (staged by Ryan 2026-09-21): execute the FGS scalability diagnostic

**Plan = `fgs_scalability_diagnostic_plan_20260921c.md` (plan C, the
consolidated synthesis). Execute from it; it is self-contained.** Plans
`..._20260921.md` (A, Ryan's revision) and `..._20260921b.md` (B, prior
agent's standalone rev 3) are provenance only. Review verdict on C
(2026-09-22): strongest of the three; adds anti-p-hacking confirmation
batches, empty-queue disambiguation (starved vs saturated), interaction
honesty, per-solve data schema + shipped analysis script; correctly
rejects whole-solve frequency normalization. Two OPTIONAL small patches
Ryan saw but did not rule on — ask or fold in at Stage 0:
1. Add the recorded R2-j64 blow-up (~9×, worse for small rungs) to the
   Stage 0 evidence table as zero-cost prior weight for H2 (sync cost).
2. C's "verify nfcache records before citing": already satisfied — the
   numbers above are certified rows (bc_certified=true @ 1e-6, cold
   isolated solves on warm cache) in the thread-scaling run dirs; cite
   those paths instead of re-verifying.

Plan-C execution order: Stage 0 desk audit (local, ≤4 threads, small) →
measurement contract → Stage 1 single-allocation reproduce+decompose
(HPC submission is Ryan-gated: prepare launcher + worktrees, get
approval) → Stage 2 targeted interventions → Stage 3 only for named
ambiguities. Campaign rules apply: annotated-tag-pinned worktrees for
FLOWPanel + all dev deps (FastMultipole), Manifest dev-paths at the
worktrees, outputs to the consolidated data root, pins recorded in a
provenance file before submitting.

## Also owed / standing (Ryan-gated ledger, carried forward)

- Plateau/pruning verdict with Ryan: FGS j*=16 is solid; the old question
  "does krylov_ilu carry weight in pruning given it never wins" is now
  MOOT in its premise (nfcache wins R4 j≥8) — re-pose as: what does
  production adopt at R4 (dagteam vs nfcache, speed-vs-memory), and what
  is prunable above 32 threads?
- Origin pushes (both repos + v21/v22/v23 + campaign tags incl. the three
  `*-20260921` tags); notebook entries owed (v21, diagnostics ladder, v22
  chunked, NUMA, v23 promotion, + these harvests — offer, don't write);
  WeakKeyDict/warmstart fix (`_publish_block_gs_status!`, from `7fbd68a`);
  `:dagteam` unit test missing in `test/runtests_unit_solver.jl`; 3 RECENT
  p018 runs awaiting archive approval; hpc-storage archive-pass report.
- Optional Ryan-gated reruns: R4 ILU j=1 with longer candidate cap; R3+
  re-runs; R1–R2 f32full arms; zen3 BLAS A/B rider.

## Traps (new this session; 20260921 prompt's traps still bind)

- Resume run dirs carry OLD job ids (13777133/13778533); first-pass logs
  live in `logs.before.<new id>/` inside each run dir — never double-count.
- Logs are output-buffered (see above) — never diagnose a hang from log
  mtime; check process CPU.
- STATUS_ilu_tune=ok does not prove budgets landed — check tune_phase2.csv
  rows (and skip duplicated header rows when parsing).
- BLAS conventions differ: R4 thread-scaling = BLAS 1; R1–R2 = BLAS j.
- nfcache-vs-FGS comparisons must state the cache-build exclusion and
  memory footprint every time — otherwise the table misleads.
- `ssh orc` needs a live ControlMaster socket (else 2FA); Slurm CLI needs
  a login shell; login banners pollute command output — filter them; judge
  runs by outputs, never sacct. Local runs ≤4 threads.

## House rules (binding, unchanged)

HPC submission Ryan-gated; monitoring via `hpc-monitor`; harvesting via
`harvester`; storage via `hpc-storage` (400 G cap); notebook writes
Ryan-gated (offer, don't write); dated status/provenance files in
BRAINSTORM/021; each campaign in its own annotated-tag-pinned worktree,
never a shared live checkout while jobs run.
