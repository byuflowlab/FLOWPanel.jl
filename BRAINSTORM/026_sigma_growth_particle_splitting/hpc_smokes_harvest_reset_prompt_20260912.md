# 026 reset prompt — smoke A PASSED; B/C running on HPC silo; harvest → commit → DELETE SILO (2026-09-12 evening)

You are picking up BRAINSTORM 026 (resolution-preserving particle splitting)
in `/Users/ryan/Dropbox/research/projects/FLOWPanel.jl` (+ FLOWVPM.jl,
FastMultipole siblings). Read `CLAUDE.md` + the policies it names.

**Predecessor prompts** (read in this order, same directory):
1. `adaptive_elongation_reset_prompt_20260911.md` — full design/state context.
2. `smoke_verify_commit_reset_prompt_20260912.md` — smoke definitions, pass
   criteria, and the RULED commit plan. Everything there holds except as
   updated below.

## State as of ~19:30 2026-09-12

### Smoke A (adaptive default, 467 steps): PASSED locally
Exit 0, clean teardown, ~5.1 h at 4 threads. Counters: 2173 elongate events
→ 6519 children (exactly 3 per event, zero violations — m=3 at f_elong=0.3),
786 compress events (first ~log-line 400), viscous=0, zero capacity skips,
zero mech_disabled, zero errors; only known-benign warnings. Log:
`/private/tmp/claude-502/-Users-ryan-Dropbox-research-projects-FLOWPanel-jl/74eae7d8-c925-45f8-83ab-41a854dafb9e/scratchpad/smoke_a_adaptive.log`.
Per-step particle counts snapshotted (before local B could overwrite VTK) to
`.../74eae7d8-.../scratchpad/smoke_a_particle_counts.csv`.
**Key reference numbers: A @ step 233 = 24275 particles; A @ 466 = 37438.**

### Smokes B and C: running on an HPC silo (Ryan-approved)
Ryan approved a NEW rsync silo on orc **on the explicit condition that the
silo is DELETED immediately after the runs are harvested. This is mandatory
— do not leave the silo behind.**

- Silo: `orc:~/silo_p026_smokes_20260912/{FLOWPanel.jl,FLOWVPM.jl,FastMultipole}`
  (sibling layout preserves the Manifest's relative dev paths; instantiated
  with cluster Julia 1.12.6, local is 1.12.5).
- Jobs (attempt 3, running since ~19:10 on m9, 16 threads / 48G / 4h,
  `--qos=normal`): **B = 13662583** (`MERGE_OVERLAP=3.5`,
  `RUN_NAME=smoke_b_mergeoverlap`), **C = 13662584**
  (`WAKE_SPLIT_ELONGATE_OVERLAP=NaN`, `RUN_NAME=smoke_c_legacy`). Both
  NREVS=0.1 (~233 steps), otherwise smoke-A env. Expected done ~1–1.5 h from
  start. Slurm logs: `silo.../FLOWPanel.jl/slurm-p026-smoke-{b,c}-<jobid>.out`.
- Two earlier submissions FAILED at t=3 s (13662572/73: `module: command not
  found`; 13662580/81: `set -u` vs `/etc/profile` HISTCONTROL). Fixed order
  in the scripts now on the silo: `set -eo pipefail` → `source /etc/profile`
  → `module load julia` → `set -u`. Judge jobs by outputs, not sacct state.
- `ssh orc` needs the live ControlMaster socket (2FA otherwise); slurm needs
  `export PATH=$PATH:/apps/slurm/latest/bin` in non-login shells.

### Local hedge run
The local detached runner (PID 78730, script/logs in
`.../74eae7d8-.../scratchpad/`) is still running smoke B locally at 4
threads (then C). Once the HPC pair passes harvest, **kill the local runner
and its julia child** (redundant; frees the machine). If HPC fails, the
local runs are the fallback evidence.

### Commit status (ruling updated by events)
- **FLOWPanel side: already in history.** Another session committed
  multi-item WIP snapshot `004ce84` (2026-09-12 15:02, at Ryan's request at
  052e.2a closure) containing ALL 026 FLOWPanel files AND the
  metadata.toml the ruling had excluded. Reported to Ryan; he did not
  object. Leave history alone. The metadata.toml is again modified in the
  working tree (local smoke writes it) — leave uncommitted.
- **FLOWVPM commit still pending, gated on B+C passing.** Branch
  `flowpanel`, exactly: `src/FLOWVPM_resolution_split.jl`,
  `src/FLOWVPM_timeintegration.jl`, `test/runtests_resolution_split.jl`,
  `examples/p026_ring_split_test.jl`. Do NOT commit untracked
  `examples/p026_ring_split_test_out/`. Message references 026 + design doc
  §19/§20.7. No further approval needed once smokes pass.

## Your tasks
1. Check jobs 13662583/13662584 (poll slurm + silo logs; monitors from the
   previous session are dead after reset). If a job died, diagnose from its
   logs; the fixed scripts are in the silo root (`smoke_b.slurm.sh`,
   `smoke_c.slurm.sh`) — resubmit from `silo.../FLOWPanel.jl`.
2. Harvest & verify:
   - B: clean exit line `SMOKE B exit: 0`; healthy counters (same grammar as
     A); **particle retention at its final step (~233) MUST EXCEED A's 24275**
     — read `NumberOfPoints` from
     `silo.../FLOWPanel.jl/data/smoke_b_mergeoverlap/smoke_b_mergeoverlap_wake1_particles/*.233.vtp`
     (or the highest-index .vtp); no split/merge churn.
   - C: clean exit; every elongate event emits exactly **2 children**
     (legacy pair2), i.e. children == 2×elongate on compress-free lines.
   - Copy the two slurm .out logs and the B particle-count evidence to the
     session scratchpad BEFORE deleting anything.
3. **Delete the silo** (`rm -rf ~/silo_p026_smokes_20260912` on orc) —
   mandatory, immediately after harvest. Also kill the local runner (task 1
   note above).
4. Make the FLOWVPM commit (above). Report A/B/C verification numbers to
   Ryan.
5. Item then returns to the §19 pending queue (commit-7 campaign arms —
   Ryan-gated; see predecessor prompt task 4).

## Ground rules (carry-over)
Local runs ≤4 threads TOTAL. HPC harness kill bug: never use harness
`run_in_background` for must-survive local jobs — `nohup ... & disown` and
poll. Campaign launches Ryan-gated, worktrees+tags per global policy. Theory
doc rewritten in place; design doc append-only. Known unrelated failure:
`runtests_unit_warmstart.jl` first testset — not ours.

## UPDATE ~19:45 — local smoke B finished; two findings

**Local B: exit 0, clean.** Zero errors, zero capacity skips. BUT it ran the
FULL 467 steps: `required_revs = max(nrevs, schedule_revs)` in the driver
(examples/rotor_hover_pressure_comparison.jl:83) — the default freestream
schedule (~467 steps at NT=36) dominates, so **NREVS=0.1 never shortened B/C**.
The HPC jobs will also run 467 steps (should still fit 4h walltime at 16
threads; watch for TIMEOUT — if hit, judge by outputs at the last written
step). The trimmed-length rationale is moot; comparisons use matching steps.

**Local B counters (467 steps): elongate=549 children=1647 (3/event),
compress=99, viscous=0. Particle retention @233 = 8879 (vs A 24275), @467 =
9978 (vs A 37438).** This FAILS the written pass criterion "B retains MORE
than A". Before calling it a bug: `merge_particles!` with sigma_relative
uses `r_pair = r_merge*sigma_min` (FLOWVPM_merging.jl:562), i.e. merge
radius SCALES WITH σ. The "merges strictly less" prediction assumed σ near
shed value (0.286σ < 0.02R absolute); in a σ-growth run, large-σ particles
get proportionally large merge radii → aggressive late-wake merging may be
the overlap gate working as designed (the σ-pump discriminator, theory doc
§4). OPEN QUESTION for Ryan/next agent: extract σ distribution at ~step 233
from B's VTK (or rerun a short probe) to decide designed-behavior vs bug,
and whether the churn guard (children re-merging: Φ_merge=3.5 vs child
target overlap) held. Do NOT mark smoke B passed/failed without this call —
report both numbers to Ryan and let him rule.

Local runner has moved on to local smoke C (4 threads, ~5h). HPC C
(13662584) will finish first; local C then becomes redundant — kill the
runner once HPC C is harvested (unless HPC C died).

## UPDATE ~20:0x — HPC smoke B finished, exit 0, MATCHES local B

HPC B (13662583): 467/467 steps, zero errors/skips, elongate=571
children=1713 (3/event), compress=85, viscous=0; np@233=8876, np@467=10065.
Local-vs-HPC agreement (8879 vs 8876 @233) confirms the low-retention
behavior is systematic, not threading/machine. Evidence copied to session-A
scratchpad (`hpc_smoke_b_13662583.log`, `hpc_smoke_b_particle_counts.csv`,
469 rows). HPC B is HARVESTED — only C's harvest still gates silo deletion.
The retention pass/fail ruling (designed σ-scaling vs bug) remains with
Ryan (see previous update).

## FINAL — all three smokes PASSED; commits done; silo deleted (2026-09-12 ~20:30)

Ryan's rulings (this session): **B = PASS** (low retention is the overlap
gate's σ-scaled merge radius acting on the aged wake — r_pair = σ_min/3.5
grows past the legacy absolute 0.02R once σ > ~2x shed value; prediction
"B retains more" was derived at shed-σ only and was wrong, not the code).
**C = PASS** (ruled on the partial log; job 13662584 cancelled early at
Ryan's order at ~step 300).

| smoke | arm | steps | elongate (children/event) | compress | np@233 | np@467 | exit |
|---|---|---|---|---|---|---|---|
| A (local, 4t) | adaptive default | 467/467 | 2173 (3.000) | 786 | 24275 | 37438 | 0 |
| B (local, 4t) | +MERGE_OVERLAP=3.5 | 467/467 | 549 (3.000) | 99 | 8879 | 9978 | 0 |
| B (HPC 13662583, 16t) | same | 467/467 | 571 (3.000) | 85 | 8876 | 10065 | 0 |
| C (HPC 13662584, 16t) | +ELONGATE_OVERLAP=NaN | ~300 (cancelled) | 354 (2.000, 0 violations) | 1 | — | — | ruled pass |

All arms: viscous=0, zero capacity skips, zero mech_disabled, zero errors.
Evidence in session-A scratchpad (`74eae7d8-...`): smoke_a_adaptive.log,
smoke_a_particle_counts.csv, smoke_b_mergeoverlap.log,
hpc_smoke_b_13662583.log, hpc_smoke_b_particle_counts.csv,
hpc_smoke_c_13662584_partial.log.

Done this session: FLOWVPM commit **f51f4ee** on `flowpanel` (4 ruled files;
untracked ring-test out dir excluded). FLOWPanel side was already in WIP
snapshot 004ce84. Silo `~/silo_p026_smokes_20260912` DELETED per Ryan's
condition. Local runner killed (local C redundant). metadata.toml left
uncommitted.

Item returns to the §19 pending queue: commit-7 campaign arms (Ryan-gated;
see predecessor prompt task 4). Open analysis thread carried forward:
σ-distribution check distinguishing healthy σ-scaled thinning from
split-merge churn under MERGE_OVERLAP=3.5 (feeds the campaign's
MERGE_OVERLAP discriminator).
