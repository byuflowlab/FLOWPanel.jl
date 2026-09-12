# 026 reset prompt — smokes running detached; commit plan RULED (2026-09-12)

You are picking up BRAINSTORM 026 (resolution-preserving particle splitting)
in `/Users/ryan/Dropbox/research/projects/FLOWPanel.jl` (panel solver) +
`/Users/ryan/Dropbox/research/projects/FLOWVPM.jl` (VPM core). Read
`CLAUDE.md` and the routing policies it names before touching code.

**Predecessor prompt** (full design/state context — read it first):
`adaptive_elongation_reset_prompt_20260911.md` in this directory. Everything
there still holds (documentation map, §19 + 2026-09-09 layer contents, key
physics facts, unit-test verification status, open questions, ground rules).
This prompt only records what changed on 2026-09-12.

## What changed today

### 1. Driver smokes are RUNNING in a detached background script

The three verification smokes (task 1 of the predecessor prompt) are running
sequentially at 4 threads via a `nohup`-detached runner, launched ~12:40 on
2026-09-12:

- Runner: PID 78730 (zsh, PPID 1), script + logs in
  `/private/tmp/claude-502/-Users-ryan-Dropbox-research-projects-FLOWPanel-jl/74eae7d8-c925-45f8-83ab-41a854dafb9e/scratchpad/`
  (`run_smokes.sh`, `smokes_driver.log`, `smoke_a_adaptive.log`,
  `smoke_b_mergeoverlap.log`, `smoke_c_legacy.log`).
- Smoke A = sane-fraction adaptive default (`WAKE_SPLIT_STRETCH=true
  WAKE_SPLIT_FRAC_COMPRESS=0.5 WAKE_SPLIT_FRAC_ELONGATE=0.3 SIGMA_CEIL=1.0
  NREVS=0.2 RUN_MONITORS=false`, 467 steps), B = A + `MERGE_OVERLAP=3.5`,
  C = A + `WAKE_SPLIT_ELONGATE_OVERLAP=NaN` (legacy pair2 arm). B and C were
  trimmed to `NREVS=0.1` (~233 steps) on Ryan's approval 2026-09-12:
  elongation onset is ~step 103 (measured), so ~233 steps suffices for B's
  churn/retention check and C's 2-children check; only A needs full length
  (compression onset ~step 400 + clean teardown are the unproven parts).
  Wall rate ~150 steps/h at 4 threads (the driver's bracketed step timer
  undercounts; don't trust it for ETAs) → A ~3 h, B/C ~1.5 h each.
  Compare B-vs-A particle retention at B's final step using A's log at the
  same step.
- Expected healthy signatures (seen in two prior partial runs, both killed
  externally, never by the physics): adaptive elongation children = 3× events
  (m=3 at f_elong=0.3), compress events firing, zero capacity skips, zero
  errors; benign warnings only (early-step empty-fmm, PressureLaplace
  deprecation). A prior 4-thread run reached step 406/467 healthy with 224
  split events. Smoke B should retain MORE particles than A (overlap gate
  merges strictly less) with no split/merge churn; smoke C should emit 2
  children per elongate event.
- Pass criteria: clean exit (`SMOKE X exit: 0` lines in
  `smokes_driver.log`, ending `ALL SMOKES DONE`), sane counters as above.

**IMPORTANT — harness kill bug**: this session's harness (Claude Code
2.1.231) repeatedly killed its OWN tracked background tasks (3 kills:
02:22, 07:09, 10:46 on 09-12), each coinciding with the session becoming
active again; the OS and other processes were ruled out (no kernel
memory-kill events, untracked jobs survived every kill). Mitigation: run
long local jobs with `nohup ... & disown` (PPID 1) and poll their logs;
do NOT use the harness's `run_in_background` for anything that must
survive.

Machine context: the 4-thread `addendum_052e2a_realsim` job (another
session's) finished ~12:30 on 09-12, freeing the ≤4-threads-total budget;
that's why the smokes run at 4 threads now. Four stale Aug-29 watcher loops
(`watch_13508968.sh`) were killed on Ryan's order.

### 2. Commit sequencing RULED by Ryan (AskUserQuestion, 2026-09-12)

**Two commits, one per repo, FLOWVPM first; metadata.toml EXCLUDED.**
Commit ONLY after all three smokes pass. Both layers (§19 fractional gating
+ 2026-09-09 adaptive elongation) go together in each commit — they are
interleaved within the same files, and the as-is tree is what the green
suites (FLOWVPM 1011/1011, wake 730/730, replay 148/148) validated.

Commit 1 — FLOWVPM (branch `flowpanel`):
- `src/FLOWVPM_resolution_split.jl` (fractional gating, adaptive elongation
  kernel, circ bookkeeping fix)
- `src/FLOWVPM_timeintegration.jl` (attempted pre-clamp Δσ² accumulators)
- `test/runtests_resolution_split.jl` (t10 adaptive elongation, t11 circ)
- `examples/p026_ring_split_test.jl`
- do NOT commit untracked `examples/p026_ring_split_test_out/`

Commit 2 — FLOWPanel (branch `flowpanel`):
- `src/FLOWPanel_wake.jl` (5-field rsplit warm-start persistence w/ legacy
  6-field migration; exposure removal; heal-comment update)
- `src/FLOWPanel_gpu_wake.jl` (device-path comment: fractional gating means
  no trigger fires device-side; host-mirror-only)
- `test/runtests_unit_wake.jl`
- `examples/rotor_hover_pressure_comparison.jl` (fraction knobs,
  `WAKE_SPLIT_ELONGATE_OVERLAP`, `WAKE_SPLIT_ELONGATE_M_MAX`,
  `MERGE_OVERLAP`)
- BRAINSTORM 026 docs: `particle_splitting_design.md` (§20/§20.7 append),
  untracked `splitting_theory.md` + the reset prompts in this directory
- EXCLUDE `data/rotor_hover_pressure_comparison.metadata.toml` (~7.9k lines
  of dev-run provenance; leave uncommitted)

The working tree also carries OTHER items' uncommitted work (018, 021, 030,
031 files) — do not touch or commit any of it. All four FLOWPanel source/test
files above were diff-verified as 026-only.

## Your tasks

1. Check the smokes (poll the scratchpad logs above; the runner is
   detached, no notification will arrive). If a smoke fails or the runner
   died, diagnose; relaunch with the same nohup pattern (≤4 threads total,
   check first that no other big job took the machine).
2. On three clean exits, verify the counters/signatures listed above (A vs
   B particle retention; C = 2 children/elongate) and report to Ryan.
3. Make the two commits exactly as ruled (no further approval needed for
   the breakdown itself; the smoke-pass gate is the only condition).
   Sensible messages referencing 026 + the §19/§20.7 design-doc sections.
4. Then the item returns to the §19 pending queue (predecessor prompt task
   3): commit-7 campaign arms (§8.4 cap030/cap018, §9 s020v matrix),
   Ryan-gated, fractions to re-derive in fraction space, per-arm
   `WAKE_SPLIT_ELONGATE_OVERLAP`/`MERGE_OVERLAP` choices —
   MERGE_OVERLAP=3.5 vs absolute A/B is the natural first discriminator
   (σ-pump finding, theory doc §4). Campaign launches: worktrees + tags per
   global policy, Ryan-gated.

## Ground rules (carry-over)

Local runs ≤ 4 threads TOTAL on the machine (check for other jobs first).
Commit only per the ruling above; anything beyond it needs Ryan. Campaign
launches Ryan-gated. Theory doc rewritten in place; design doc append-only.
Known unrelated pre-existing failure: `runtests_unit_warmstart.jl` FIRST
testset (immutable WarmstartNoopSolver as WeakKeyDict key) — not ours.
