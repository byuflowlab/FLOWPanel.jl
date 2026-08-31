# Handoff — BRAINSTORM 021: fix the cost tuner's noise floor, retune every rung, relaunch

Written 2026-08-24 at a context reset. **Nothing is queued. Nothing is running.**
The cluster is clean and waiting for you.

## Ryan's plan for you, in order

1. **Make the tuner adjustment in PLAN MODE first.** Ryan wants to see the plan
   before you change code. It is a small change — see "The fix" below.
2. **Retune every rung** R1–R7 with the corrected tuner. Do not reuse any
   existing knobs (reasons below — this is not optional caution, the old ones are
   known-suspect).
3. **Then release/relaunch the rest of the campaign**, making sure every job uses
   the new tuning parameters.

**Job submission is inline work and Ryan's call — get his go before `sbatch`,
never route it to a subagent.**

## The state you are inheriting

Clean slate, deliberately:

- **No 021 jobs queued or running.** Everything was cancelled.
- **All `benchmark/results/**/*.csv` purged** on the cluster. The seven
  `bcache_R*.bin` were KEPT — they are frozen-`b` direct N² assemblies,
  independent of both FMM generation and tuning, and R7's is expensive to rebuild.
- Local repo has the two-stage tuner already wired in (see below); it is synced to
  the cluster and md5-verified.

## The problem, fully diagnosed — do NOT re-derive any of this

### Stock `tune_fmm` optimizes the wrong thing

`FastMultipole.tune_fmm` tunes a **single evaluation** against an error
tolerance. It never measures iterative-solve cost, so it accepts any knob set
that is accurate enough — including a `leaf_size` that is ~2× more expensive per
apply, which is exactly what Phase 1 exists to measure.

Measured at R1 on the current FM (job 13407693, fixed p=17 / mac=0.5,
`benchmark/fm_leaf_ab.jl`):

| leaf | 9 | 21 | 40 | 60 | 71 |
| --- | --- | --- | --- | --- | --- |
| `t_solve_min` | **15.51 s** | 20.65 | 26.10 | 32.64 | 34.95 |

Monotone, 2.25× across the range, at constant iteration count (~57).

### There is NO FastMultipole regression

At matched knobs (leaf 9) the old and new FM agree to **+6%** (14.59 → 15.51 s).
The apparent 2.5× slowdown was entirely the knob. The 051 work is fine. This was
nearly mis-filed as an FM performance regression — do not re-open it.

### The old ladder is also suspect

Old-campaign picks were leaf **9 / 59 / 58 / 32 / 25 / 20** across R1–R6. R1's 9
was luck, not a better tuner — the same accuracy-only objective produced 59 and
58 at R2/R3. **So parts of the previously-accepted ladder may carry the same ~2×
cost error.** This is inferred from R1's sweep, NOT measured at R2. If you want it
settled, `benchmark/fm_leaf_ab.jl` with `RUNG=R2` against its old leaf 59 costs
~15 min on one node and tells you whether the published ladder has a systematic
problem. Worth doing — it matters beyond this campaign.

This is also why **every rung must be retuned** rather than reusing old knobs.

## The fix (what plan mode is for)

`benchmark/rotor_hover_solver_phase1_tune_hpc.jl` already runs two stages:
stock `tune_fmm` for an accuracy-valid seed, then
`FastMultipole.tune_fmm_perturb` (`src/autotune.jl:337`) to coordinate-descend on
**benchmarked wall-clock** under the same tolerance. `PERTURB_TUNE=0` reverts to
seed-only (records variant `tuned_seed_only`).

That machinery is correct and **already probes leaf in both directions**. Its
trace at R1:

```
start:  P=17 MAC=0.5 leaf=68  t=0.662 s
iter 1: leaf=102 t=0.943 s          <- upward probe, worse
iter 1: leaf=45  t=0.444 s   -> move
iter 2: leaf=30  t=0.450 s          <- NOT taken
minimum: leaf=45
```

**It stopped on measurement noise.** Leaf 30 measured 1.4% *worse* than 45 when
the truth there is clearly better (the sweep predicts ≈0.40 s). With `reps=2` the
benchmark cannot resolve the real difference. Note this was a measured
*regression*, not a below-threshold gain — so lowering `improve_tol` does NOT
help; only better sampling does.

Result: leaf 68 → 45 is a genuine improvement (~35 → 28 s) but still ~1.8× off
the ~15.5 s available at leaf 9.

**The change Ryan approved:**

- **`reps=2` → 7** in the `tune_fmm_perturb` call. This is the actual fix.
- Add a **`TUNE_MAX_SECONDS` wall-clock guard** on the stage so a pathological
  rung cannot stall a 3-day job. Free at R1 (<1 s per evaluation) but real at R7
  (419k panels, 7 reps × ~12 candidates).
- **Consider scaling reps by rung** — 7 at R1–R4, fewer at R5–R7 where each
  evaluation is long enough that noise is proportionally smaller.

Do NOT write a custom leaf ladder. That was considered and rejected: `perturb`
already does bidirectional coordinate descent over `(p, mac, leaf)` and correctly
rejected `p=16` and `mac=0.55` on accuracy. The only defect is the noise floor.

**Validation target:** rerun R1's tune. A working tuner should walk
68 → 45 → 30 → 20 → 13 → 9 and land near **leaf 9**, and R1's
`krylov_gmres/standard` `t_solve_min` should come back near **15.5 s**. If it
stalls at 45 again, `reps` is still too low or the wall-clock guard fired — check
the descent trace in the job's `.out` (`grep -A30 'Begin FastMultipole.tune_fmm_perturb'`).

Open question nobody has answered: the sweep never sampled **below leaf 9**, so
the true optimum may be lower. FMM cost is U-shaped in leaf (large = too much
direct work, small = tree/translation overhead), and at R1/p=17 we are on the
"too large" branch across the entire sampled range. The corrected descent should
find the turn on its own; just do not assume 9 is the floor.

## Relaunching after validation

Submit per rung: tune → table, chained **`afterok`** (see trap 1). Then Phase 2a
R1–R5 and Phase 2b R1–R4 the same way. A prior fan-out script lives at
`~/flowpanel-021/fanout2.sh` on the cluster — its structure, walltimes and config
sets are correct, but it predates the tuner fix, so re-read it rather than running
it blind.

Cost-ceiling drops from `ledger.md`, per config not per rung: `backslash_ldiv`
drops at R6–R7 in **both** modes; `krylov_gmres` additionally drops from
**single** mode at R6–R7. R2–R5 use all six configs in both modes.

Right-sized walltimes (from measured 2026-08-18/19 timings, ~2× margin). Do NOT
use a blanket 7 days — a priority-0 `standby` QOS backfills accurate requests far
better than padded ones:

| | tune | table multi | table single | p2a | p2b |
| --- | --- | --- | --- | --- | --- |
| R1 | 4 h | 2 h | 4 h | 2 h | 2 h |
| R2 | 6 h | 3 h | 8 h | 3 h | 4 h |
| R3 | 12 h | 4 h | 12 h | 6 h | 8 h |
| R4 | 1 d | 8 h | 1 d | 12 h | 16 h |
| R5 | 2 d | 12 h | 2 d | 1 d | — |
| R6 | 3 d | 1 d | 3 d | — | — |
| R7 | 3 d × 3 stages | 2 d | 7 d | — | — |

Phase 2a runs **six** configs, not nine: the `*_nfcache` variants consume 2b's
`tune_cached.csv`, which those rungs do not have yet — a follow-on round.

R7 has never been tuned and needs its 3-stage chain (`STAGES=1`, `2`, `3`)
because the full tune does not fit one job. Ryan's ruling: **R7 does not gate the
other rungs** — launch the rest without waiting for it.

### Three traps that have already cost real time

1. **`afterany` treats a CANCELLED upstream as SATISFIED.** 21 held jobs had to be
   cancelled because releasing them would have run tables with no knobs at all.
   Use `afterok` for tune→table.
2. **Resume support makes stale rows sticky.** The drivers now skip any
   `(rung, config, row_kind)` already in the CSV — correct for preemption, but it
   means bad rows are silently inherited and never recomputed. **Purge a rung's
   CSVs before re-running it.** Keep `bcache_R*.bin`.
3. **`nproc` lies here.** GNU `nproc` honours `OMP_NUM_THREADS`, which Slurm
   presets to 1 under `--cpus-per-task`, so it reports 1 while real affinity is
   `0-127`. Use `SLURM_CPUS_ON_NODE`. Also `--ntasks=64 --cpus-per-task=64`
   requests 64×64 = **4096** CPUs.

## Environment — verified working, do not re-derive

- **Run everything from `~/flowpanel-021/`** on the cluster. Isolated from
  `~/projects/`, which is shared with 018/023 — two agents collided there on
  08-24. Layout: `FLOWPanel.jl` + sibling `FastMultipole` (symlink to the `.jl`
  clone; the Manifest dev path is `../FastMultipole`) + `FLOWVPM.jl`.
- Clones are clean: FLOWPanel `321473f`, FastMultipole `d8258a7d`, FLOWVPM
  `4494c25`. `fm_commit` resolves with **no `-dirty`** — the old `a9b734ad-dirty`
  pin was an unversioned blob in no git history, and two different FMM
  generations briefly claimed that same string.
- `benchmark/` and `Manifest.toml` are **untracked upstream** and must be rsynced
  in after any clone. Leave `Project.toml` alone — `PythonPlot` is an inert
  `[weakdeps]` entry; the 08-22 killer was `GeoIO`, a hard dep, already gone.
- **Hardware pin, in the SBATCH header of all seven launchers — do not override:**
  `--constraint=zen3 --exclusive --mem=500G`, 128 CPUs allocated, `THREADS=64`.
  Exclusive is deliberate: these partitions are `OverSubscribe=YES:4`, and a
  64-CPU cgroup can bind all threads to one socket, so a half-node allocation
  changes timings at the same thread count. Override `--partition`, `--qos`,
  `--time`, `--job-name` only.
- Three **spec-identical** zen3 pools (2×64c, `ThreadsPerCore=1`, 512 GB):
  `m12` (public, 3 d), `physics2` (7 d, 10 nodes), `m12pws` (28 d, 1 node) — the
  latter two via `--qos=standby`, preemptible; launchers set `--requeue`, which is
  safe because the drivers resume.
- BLAS threads are explicitly pinned AND asserted in both modes in
  `benchmark/common.jl`. Multi mode previously let OpenBLAS default to the machine
  core count — 128 BLAS threads against Julia's 64 on an exclusive node.
- `ssh orc` needs a live ControlMaster socket; if dead, ask Ryan to run
  `! ssh orc`. Slurm binaries live in `/apps/slurm/latest/bin`. Julia needs a
  login shell: `ssh orc 'bash -lc "module load julia && ..."'`.

## Also open, not blocking

- **R7's `p1-table-R7-single` has a 7-day request and unknown runtime.** This is
  the unsolved R7 splitting problem. `m12pws`'s 28-day ceiling is the escape
  hatch; per-config splitting via the colon-separated `CONFIGS` export is the
  mechanism.
- `benchmark/fm_leaf_ab.jl`'s `bc_rel_l2` column is **invalid** (wrong args to
  `bc_error!`); its timings are sound. Accuracy parity rests on the campaign rows
  (`krylov_gmres` certified BC 9.674e-7 → 9.675e-7), not on that sweep.
- All three repos are committed and pushed, but still carry untracked
  working-tree files; `benchmark/` is essentially untracked and exists in one
  place plus the cluster rsync.
- Phase 2 unsteady (`p2_unsteady.sh`) is blocked by design — the driver refuses to
  append to an `unsteady.csv` whose header is not current. Point it at a fresh path.
- Phase 3 (`benchmark/slurm/p3_warmstart.sh`) written, never submitted, needs
  Ryan's separate go. Its trajectory-identity guard flags DIVERGED above 1e-8
  while post-fix smoke arms differ 7.9e-8 / 2.4e-7 — needs Ryan's call first.
- Narrative entry point is `log.md` (**newest first**). Do NOT read the item
  file's `## Current status` (stale); do not read `phase_01_consistency.md`
  (73 KB) or `phase_02_single_step_benchmarks.md` (41 KB) end-to-end — grep them.
