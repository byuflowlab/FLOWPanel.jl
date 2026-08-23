# 021 Session Log

Newest first. Narrative only — results go to the phase files and `ledger.md`.

## Dated entries

### 2026-08-22 (night) — Gaussian rerun executed: scope was 3 artifacts, not the campaign

Executed the pre-Phase-3 rerun. The scope question ("what changes under
Gaussian?") resolved to **two mechanisms**, and knowing which one is active in a
given driver decides everything.

**Channel A — the filament kernel.** Only `_bound_vortex_velocity` /
`_bound_vortex_gradient` (`src/FLOWPanel_elements_fmm.jl:938-1067`) branch on the
family, reached from vortex-ring panel edges (`:862,867`), panel-wake trailing
filaments (`src/FLOWPanel_wake.jl:2891-2897`), and the semi-infinite gradient
helper (`:1627`) — and **only on velocity/gradient passes**. Potential-only
evaluation routes through `ConstantDoublet` and never reaches a filament.
Vortex particles / FLOWVPM are entirely unaffected.

**Channel B — FMM buffer radii (the non-obvious one).**
`radius_inflation(::Type{VortexRing}, core_size, tol)`
(`src/FLOWPanel_elements_fmm.jl:1130-1142`) reads the same global but dispatches
on the body's ELEMENT TYPE, not on whether a filament is evaluated. It feeds the
source-buffer radius (`src/FLOWPanel_abstractbody.jl:1200`,
`src/FLOWPanel_wake.jl:355,2815`), moving the near/far split. So it can change
timings in a run with no wake at all.

**Measured, geometry only (no solve), both ends of the ladder:**

| offset | Δ(buffer) vat−gau | vs smallest panel radius, R1 | vs smallest, R7 |
| --- | --- | --- | --- |
| `core_size_panel = R·1e-10` (all solve + tuning runs) | 3.773e-10 m | 1.7e-5 | 1.5e-4 |
| `core_size_targets = 1e-3` (`steady!` target pass) | 3.171e-2 m | 1413× | 12321× |

(R1 panel radius min 2.24e-5 max 4.43e-3 m; R7 min 2.57e-6 max 5.66e-4 m.)

**This is why the FMM tuning runs do NOT need redoing** (Ryan's question). Every
tuning artifact — `tune.csv`, the FGS τ-ladder, `fgsprecond.csv`, the p2b
cached-economics tune — runs at the panel offset, where the family moves the
buffer by 0.0017% of the smallest panel at R1 and 0.015% at R7. The optimal
(p, MAC, leaf) is family-invariant. Corroborated by an R1 A/B through
`FastMultipoleBackend` at tuned knobs (so the near/far split was genuinely
exercised): rms_b, ‖x‖, Σx bit-identical and the SAME 118 iterations.
CAVEAT, stated honestly: the panel-offset perturbation GROWS ~9× from R1 to R7
as panels shrink, and 1.5e-4 relative is not machine-zero. A marginal MAC
comparison sitting within 0.015% of threshold could in principle flip at R7.
Bit-identity is measured at R1 only; R7 rests on the scaling argument.

#### What was rerun

| # | Artifact | Result |
| --- | --- | --- |
| 1 | `results/phase1/multi/agreement.csv` | **NOT sensitive — prediction refuted by measurement, NO rerun needed.** See the correction below; the file is restored as the record of evidence. |
| 2 | `results/phase2/multi/unsteady.csv` | Rerun R1 backslash+fgs under Gaussian; now 12 rows, all `filament_reg=GaussianRegularization`. |
| 3 | `rigid_motion_tree_reuse/fgs_staleness.csv`, arm `fgs_geometry_order_fixed_8step` | Rerun under Gaussian. **Reproduced BIT-FOR-BIT** — 171/171 comparable cells identical, divergence 1.740e-7 → 1.886e-7 flat, growth 1.1×. |

Result 3 is the useful one: it **resolves the provenance indeterminacy**. The
arm was already Gaussian (as suspected — the in-repo default became Gaussian
mid-day 2026-08-20 and this arm ran 08-21), so the Phase-3 trust gate
(`phase_02...md:313`, "NO Phase-3 unsteady FGS row is trusted until the
discriminator runs") is certified in the campaign's kernel. No further action.

#### CORRECTION — Channel B is a COST channel, not an ACCURACY channel

Predicted, from the 12321× buffer difference at the targets offset, that
`agreement.csv` would move. **It does not.** Rerunning R1 (all 7 configs) under
Gaussian reproduces `CT_laplace` to 8 decimals against rows dated 2026-08-14 —
which are definitively Vatistas, since the kernel was hardcoded and not
selectable until 08-20:

| config | CT_lap gaussian | CT_lap vatistas (08-14) | Δ% |
| --- | --- | --- | --- |
| backslash_ref | 0.06040613 | 0.06040613 | 0.0000 |
| backslash_coupled | 0.06040565 | 0.06040565 | 0.0000 |
| krylov_gmres | 0.06040636 | 0.06040636 | 0.0000 |
| krylov_jacobi | 0.06040614 | 0.06040614 | 0.0000 |
| krylov_ilu | 0.06040613 | 0.06040613 | 0.0000 |
| fgs | 0.06040613 | 0.06040613 | 0.0000 |
| fgmres_fgs | 0.06040607 | 0.06040607 | 0.0000 |

**Why the prediction was wrong.** `radius_inflation` changes WHICH panels are
classified near vs far — it does not change the VALUE the FMM returns, because
both the near and far paths evaluate the same influence to the same tolerance.
So Channel B moves **cost**, never numbers. Only Channel A (the filament kernel
actually evaluated on a velocity/gradient pass) moves numbers, and that requires
filaments to exist — i.e. a wake.

Revised rule, and the one to carry forward:

- **Channel A → numbers**, wake-carrying runs only. Confirmed by the 328 vs 329
  particle split below.
- **Channel B → timings only**, any `VortexRing`-typed body, wake or not.

Consequence: the rerun scope was even smaller than planned — only the
wake-carrying unsteady artifacts ever needed it. R2/R3 of `agreement` were
cancelled once R1 established the null; `agreement.csv` is restored unchanged as
the record, and the Gaussian R1 confirmation is kept as
`agreement_gaussian_verify_R1_20260822.csv`.

Superseded originals are in `quarantine_vatistas_20260822/` beside each artifact,
with a README; they carry no `filament_reg` column and must not be pooled with
post-08-22 rows.

#### Same-code family A/B (new measurement) — `unsteady_family_ab_20260822.csv`

Old-vs-new could NOT isolate the family: the 08-18 rows predate three days of
code changes (and `nsim` is derived from `maneuver!` deltas, so `N_STEPS=5` now
yields 6 recorded steps vs the old 5). So both families were run on TODAY's code:

| metric | gaussian | vatistas | vat−gau |
| --- | --- | --- | --- |
| Σ t_solve, backslash | 26.432 s | 23.723 s | **−10.25%** |
| Σ t_solve, fgs | 30.333 s | 30.851 s | **+1.71%** |
| particles at final step | 328 | 329 | +1 |

**The robust finding is the particle count: 328 vs 329.** That is deterministic
and load-independent, and it proves the family changes the unsteady trajectory —
so 025's CT null (+0.023%) does not license treating wake-carrying 021 work as
family-invariant.

**The timing columns are NOT trustworthy as stated** and are recorded as
indicative only: they are single runs (the campaign's own ruling 5/7 requires
min-of-k, k≥5, after warmup), the two families ran in separate batches so
machine thermal/load state is confounded with the family, and the signs
DISAGREE between configs (−10% backslash vs +2% fgs) — which is what noise
domination looks like, since a pure near/far-cost effect should move both the
same way. A proper measurement belongs on HPC under the min-of-k protocol.

#### Instrumentation fix — `filament_reg` now recorded (SCHEMA AMENDMENT)

Root cause of the whole indeterminacy: nothing recorded the family. Fixed
centrally in `assert_and_banner()` (`benchmark/common.jl`) so every 021 driver
inherits it — banner text, returned NamedTuple, `RUNS_CSV_COLUMNS`, and the
`validate_runs_csv` required-non-empty list. Pattern lifted from
`benchmark/fm051_pass_parity.jl:922` (previously the ONLY recorder in the tree).
Row-construction sites updated in `rotor_hover_solver_phase1.jl`,
`rotor_hover_solver_smoke.jl`, and the three rerun drivers.
`decision_rules.md` amended (the schema is pinned there).
Verified live: the new unsteady rows carry `filament_reg=GaussianRegularization`.

#### Latent bug found, NOT fixed

`_U_boundvortex` (`src/FLOWPanel_elements_fmm.jl:1559-1615`) is **hardcoded
Vatistas** and ignores the global, while its gradient counterpart
`_U_boundvortex_gradient` (`:1627`) honours it. Only fires on the
`semiinfinite_wake=true` path, which 021 does not use — latent, not live.
Separate fix; flagged here so it is not rediscovered.

## OPEN TODO — relabel the 10 in-flight jobs after they finish (Ryan, 2026-08-22)

The ten jobs queued on 2026-08-22 (13306549/50/54/56/57, 13306599/600,
13306603/604/605) were submitted BEFORE the launcher default was flipped to
Gaussian. Slurm snapshots the batch script at submit time, so they will run and
self-report as `FLOWPANEL_FILAMENT_REG=vatistas`.

**Ryan's ruling: let them run.** Their numbers are family-invariant by
construction (Phase 1 / 2b are potential-only at the panel kernel offset —
measured bit-identical, see the 2026-08-22 evening entry), so the label is a
provenance artifact, not a correctness problem. Cancelling would have cost queue
position on a cluster at 128/136 nodes for zero numerical change.

**ACTION WHEN THEY FINISH — do not skip:** update the family attribution on
their outputs so the campaign reads as uniformly Gaussian.

- [ ] 13306549 p2b-nearfield-R1
- [ ] 13306550 p2b-nearfield-R2
- [ ] 13306599 p2b-nearfield-R3
- [ ] 13306600 p2b-nearfield-R4
- [ ] 13306554 p1-table-R6-multi
- [ ] 13306556 p1-table-R6-single-a
- [ ] 13306557 p1-table-R6-single-jacobi
- [ ] 13306603 p1-tune-R7-s1
- [ ] 13306604 p1-tune-R7-s2
- [ ] 13306605 p1-tune-R7-s3

For each: relabel the `filament_reg` attribution on the produced rows
(vatistas -> gaussian) and annotate in `ledger.md` that the row is
family-invariant by measurement rather than by the pin it ran under. Note these
runs predate the `filament_reg` CSV column reaching the cluster, so the
attribution lives in `banner.txt` / the launcher snapshot, not in the rows
themselves — check `scontrol show job <id>` output or the job's `slurm-*.out`
banner to confirm what each actually ran under before relabelling.

### 2026-08-22 (evening) — Campaign moved to Gaussian regularization; MEASURED inert for Phase 1/2b

**Ryan's ruling (2026-08-22):** the campaign uses Gaussian filament
regularization from here on; the earlier vatistas-era rungs are to be RERUN
under Gaussian later. **TODO carried forward — see "Gaussian rerun backlog"
below.**

**Measured first, before acting.** Slurm snapshots the batch script at submit
time, so editing a launcher cannot reach an already-queued job — applying the
switch to the 10 in-flight jobs would have meant cancelling and resubmitting all
of them. That is only worth doing if the flag actually changes the answer, so it
was tested rather than assumed. R1, GMRES to rtol 1e-12, tuned knobs, 2 threads,
identical in every respect but `FLOWPANEL_FILAMENT_REG`:

| reg | rms_b | norm(x) | sum(x) | niter |
| --- | --- | --- | --- | --- |
| vatistas | 0.011505099103103303 | 6.437344805815382 | 11.41656294891941 | 118 |
| gaussian | 0.011505099103103303 | 6.437344805815382 | 11.41656294891941 | 118 |

Bit-identical to the last digit, same iteration count. **Root cause: the
regularization is unreachable in this case.** `phase1_case.jl` builds a frozen
single-step Dirichlet solve — RHS is a pure body->body direct source-potential
assembly, `semiinfinite_wake=false`, and there is no time-marching, so no wake
is ever shed and no filament kernel is evaluated. Every one of the 10 queued
jobs (p1_tune, p1_table, p2b_nearfield) includes `phase1_case.jl`, so the flag
is inert for all of them.

**Therefore the 10 queued jobs were NOT cancelled** — resubmitting would have
cost a queue cycle and changed no number. The launcher DEFAULTS were flipped to
`gaussian` in all six `benchmark/slurm/*.sh` (md5-verified identical on both
checkouts), so the switch applies from the next submission onward and, more
importantly, to the wake-carrying `p2_unsteady.sh` arm where it is NOT inert.

**Consequence for the ledger:** the existing Phase-1 rows are unaffected by this
switch — they are filament-free by construction, so vatistas-era and Gaussian-era
Phase-1 rows ARE mutually comparable. The regularization only partitions the
evidence for wake-carrying work (Phase 2 unsteady, and campaigns 018/022/025).

#### Gaussian rerun backlog (opened 2026-08-22)

To revisit once the current queue drains. Nothing here is urgent, and per the
measurement above most of it may be a no-op:

- [ ] Phase-2 unsteady arms (`p2_unsteady.sh`) — the ONE place in 021 where the
      regularization is known to be live. Any pre-08-22 unsteady rows were taken
      under vatistas and need rerunning under Gaussian before being pooled with
      post-08-22 rows.
- [ ] Confirm the inertness result holds at a wake-carrying Phase-2 state (the
      08-21 FGS work used a "first particle-carrying R1 state"); the R1 probe
      above only establishes it for the frozen wake-free case.
- [ ] Phase-1 R1-R6 rungs + `tune.csv` / `bcache_R*.bin`: rerun NOT required on
      regularization grounds (filament-free), but note in the ledger that they
      were produced pre-switch so the provenance is explicit.

### 2026-08-22 (later still) — Walltime audit: R7 tune restructured into stages; p2b R3/R4 raised

Before any of the 8 resubmitted jobs started, an elapsed-vs-N fit on the
COMPLETED rungs showed two of them could not finish in the walltime requested.
The fit validates: it predicts R5-single at 30-57 h (mid 38 h) against an actual
44.6 h, so mid-range numbers are good to roughly +-40%.

| pipeline | measured rungs (h) | exponent in N | R6/R7 extrapolation |
| --- | --- | --- | --- |
| p1-tune | R1 0.29, R3 1.29, R4 4.32, R5 10.46, R6 39.05 | 1.18-1.96 | **R7 87-148 h (mid ~114)** |
| p1-table multi | R1 0.08, R2 0.22, R3 0.36, R4 0.87, R5 2.32 | 0.77-1.58 | R6 4-7 h |
| p1-table single (all cfg) | R1 0.81, R2 3.37, R3 6.38, R4 15.37, R5 44.6 (a+gmres) | 1.06-2.12 | R6 ~61 h across the two split jobs |

`n_panels` doubles per rung (R1 8,016 -> R7 419,276) and the `m12` partition
MaxTime is **3-00:00:00**. So R7's tune at ~114 h does not fit a single job at
ANY allowed walltime, let alone the 48 h first requested. This is structural,
not a sizing slip, and it extends to the rest of the R7 row:
`p1-table-R7-single-a` extrapolates to ~100 h and the ledger already had
`R7-single-jacobi` at ~69 h ~= its limit. **The R7 row needs a splitting
strategy before it is submitted, not a bigger `--time`.**

**Actions taken.** Cancelled 13306551/52 (p2b R3/R4) and 13306553 (p1-tune-R7)
while all were still PENDING, so nothing was wasted.

- p2b R3/R4 resubmitted as 13306599/13306600 at `3-00:00:00`. Reason: raising
  the per-trial cache-build cap 60 s -> 1800 s raised the tuner's worst-case
  cost ~30x across a 6-MAC sweep, but the script's 24 h walltime had not been
  raised to match. R1/R2 left at 24 h. NOTE these are the weakest estimates in
  the campaign — no p2b job has completed at ANY rung, so there is no measured
  reference, only the R1 local smoke.
- `p1_tune.sh` gained a `STAGES` gate (colon/comma separated subset of `1:2:3`,
  **default `1:2:3`** so every other rung's invocation is unchanged). R7
  resubmitted as three chained jobs 13306603 (s1 tune_fmm) -> 13306604
  (s2 fgstune) -> 13306605 (s3 fgsprecond), `afterok`, 72 h and 256G each.
  `afterok` is correct here — unlike the cascade that killed 7 jobs on 08-20,
  this is a genuine data dependency (s2 consumes s1's `tune.csv`). Verified
  `results/phase1/multi/R7/` is empty first: the stage scripts APPEND rather
  than checkpoint-skip, so re-running a completed stage duplicates rows.

**Latent bug found and fixed: `p1_tune_s3.sh` was missing the filament pin.**
It is a pre-08-20 fork of `p1_tune.sh` that never received the BRAINSTORM/025
`FLOWPANEL_FILAMENT_REG=vatistas` export, so any stage-3 recovery run since
08-20 would have silently used the new Gaussian default and produced rows not
comparable with the rest of Phase 1. It was also **cluster-only** — never in the
local repo, and it survived the 08-22 sync only because the `benchmark/` rsync
carried no `--delete`. Now pulled into the repo, pin added, and marked
DEPRECATED in favour of `STAGES=3`; forking the launcher again is what caused
the drift, so future splits go through the gate.

### 2026-08-22 (later) — Fixes synced to both cluster checkouts; all 8 jobs resubmitted

**Sync method (Ryan's ruling):** rsync the WHOLE local `src/` (with `--delete`)
and `benchmark/` to both `~/projects/FLOWPanel.jl` and
`~/projects/p2b/FLOWPanel.jl`, rather than committing or cherry-picking files.
Rationale: the local tree is dirty with mixed 018/022 work and all of
`benchmark/` is untracked, so any partial sync would have left a mixed
Aug-19/Aug-22 generation inside `src/` — the exact failure that already cost
this campaign a run. Whole-tree rsync makes the cluster source byte-for-byte
verifiable against local.

**Provenance (verified post-sync, both checkouts):**

| item | value |
| --- | --- |
| local branch / HEAD | `fastmultipole` @ `d6bf8b6` + uncommitted 021 fixes |
| `src/FLOWPanel_fmm.jl` | `4128f687d3bf1d161e0cdabea640dac9` |
| `src/FLOWPanel_instrumentation.jl` | `a08ac333b20d17da523d6d298024b493` |
| `src/FLOWPanel_solver.jl` | `114428fa7c6aa70528fb2f77a1004dca` |
| `benchmark/rotor_hover_solver_phase1_table.jl` | `3f9d5b8ff6c104e7f16b6091dbbca9df` |
| `benchmark/rotor_hover_solver_phase2b_nearfield_cache.jl` | `58c2e32fd93ca8617bd259ee74f0d8cb` |
| `benchmark/slurm/p2b_nearfield.sh` | `7b9d338ddcbeb352e881c86c93ada974` |
| full `src/*.jl` per-file md5 list | diffs CLEAN vs local on both checkouts |
| `grep -c nearfield_cache_max_bytes src/FLOWPanel_solver.jl` | 7 (was 0) |
| `grep -c 'max_pattern_entries=8192' .../phase1_table.jl` | 1 (was 0) |
| `import FLOWPanel` | `LOAD_OK_MAIN`, `LOAD_OK_P2B` (FLOWPanel recompiled, 240 deps cached) |

The first rsync pass silently transferred nothing and was caught only by the
md5/grep check — a direct vindication of the "verify on the cluster, never
assume" rule. Re-run individually, it landed. Cluster-generated flat-layer CSVs
under `benchmark/results/phase1/multi/` were NOT clobbered (identical size+mtime,
rsync skipped them); the per-rung `R1`–`R6` dirs were untouched and the p2b
R1–R4 knob CSVs + `bcache_R*.bin` prerequisites are intact.

**Failure class 3 (EACCES) is closed by dependency surgery, not by lock cleanup.**
Ryan removed GeoIO and PythonPlot from `Project.toml` on both checkouts on 08-22,
so `import FLOWPanel` no longer runs PythonCall's `__init__` -> `CondaPkg.resolve()`
and there is no `.CondaPkg/lock` to contend for. Deps lists on the two checkouts
are now identical to local's. The new `FLOWPanel_gpu_influence.jl` /
`FLOWPanel_gpu_wake.jl` (task 051) add no dependencies and are inert unless armed.

**Correction to the 08-22 entry above:** 13242660/62/64/66 were NOT
"never submitted". They are in `sacct` as `CANCELLED`, Reason `Dependency`,
00:00:00 elapsed — their `afterok` parent failed and Slurm cancelled the chain.
Same for `p1-table-R6-single-a/-jacobi` (13206078/79) and ALL THREE
`p1-table-R7-*` (13206081/82/83). So the campaign lost 7 jobs to dependency
cascade on top of the 3 that actually failed. The afterok chain was not
malformed; it was doing what afterok does. Resubmission therefore drops the
cross-job chaining wherever the jobs do not share an output directory.

**Resubmitted 2026-08-22, all zen3 + `--exclusive` (128-core / 512 GB nodes):**

| job | id | notes |
| --- | --- | --- |
| p2b-nearfield-R1..R4 | 13306549/50/51/52 | from the p2b checkout, four INDEPENDENT jobs (no afterok), `--mem=256G`, `MODE=multi`, `TUNE_MACS=0.3:...:0.8` |
| p1-tune-R7 | 13306553 | retry of 13206080, `--mem=256G`, `--time=48:00:00` |
| p1-table-R6-multi | 13306554 | retry of 13206077, needs the 8192*N ILU fix; 5 configs (drops `backslash_ldiv` per the R6/R7 drop-out recipe) |
| p1-table-R6-single-a | 13306556 | `krylov_ilu:fgs:fgmres_fgs`, 192G/3-00:00:00 |
| p1-table-R6-single-jacobi | 13306557 | `krylov_jacobi`, 128G/3-00:00:00, `afterany:13306556` |

`--mem` raised to 256G on the p2b rungs (from the script's 64G default) because
the near-field cache cap is now 32 GiB; under `--exclusive` this is free and does
not affect timing comparability. The two R6 `-single` jobs remain serially
chained — they share `results/phase1/single/R6/` and the NFS single-writer rule
applies — but via `afterany`, so an upstream failure no longer silently deletes
downstream work. `p1-table-R7-*` cannot go back in until 13306553 produces R7's
knob CSVs.

Not yet resubmitted / still open: `p1-table-R7-multi`, `p1-table-R7-single-a`,
`p1-table-R7-single-jacobi` (blocked on p1-tune-R7).

### 2026-08-22 — Cluster campaign was 0/8; three independent causes found and two fixed

Live `sacct` contradicted the docs: of the eight Phase-2b jobs treated as in
flight, **four were never submitted** (13242660/62/64/66 absent from Slurm
history — the afterok chain needs checking before relaunch) and four failed.
Phase-1 was also incomplete: `p1-tune-R5/R6` and `p1-table-R5-multi` COMPLETED,
but `p1-tune-R7` (13206080) and `p1-table-R6-multi` (13206077) FAILED. The four
p2b failures were NOT one cause, contrary to the single error text first seen:

1. **Near-field cache caps (R3, R4).** R4 needed 4.47 GiB for 207,344 blocks
   against a 4 GiB cap — a 12% overshoot. R3 tripped the tuner's first-trial
   guard at leaf 58 / MAC 0.3. Root cause of the abort-rather-than-degrade
   behaviour: `tune_fmm` clamps back to `last_feasible_leaf`, but that is
   initialised to `nothing` and `&&` short-circuits, so when the FIRST trial is
   infeasible the `error()` fires and the clamp path below it is dead code
   (`FastMultipole/src/autotune.jl:109-121`). The leaf ladder only ever grows,
   so the message's own advice (start lower) is the only route the code supports.
   Compounding it, the MAC sweep starts at 0.3, its most expensive point — at R1
   the near field is 85% of the full matrix there vs 53% at MAC 0.5.
   FIXED: caps threaded from `KrylovSolver` down to `build_nearfield_cache!`
   (they were previously **unreachable from FLOWPanel** — `FLOWPanel_fmm.jl`
   called the builder with no `max_bytes`, so R4 was hard-blocked); driver
   defaults raised to 32 GiB / 1800 s, per-rung overridable via
   `TUNE_MAX_BYTES` / `TUNE_MAX_BUILD_TIME`.
   NOTE the cache build is **serial** (`nearfield_cache.jl` `_build_nearfield_cache`
   has no `@threads`; the threaded code is `nearfield_matvec!`, the evaluation),
   so the time cap is wall-clock, not a per-thread budget — measured R1 19.5 s at
   4 threads vs 20.0 s at 1. `estimate_nearfield_cache` overshoots actual build
   time by a consistent ~1.65x.
   The old 60 s cap was distorting the result, not just aborting jobs: with caps
   raised the R1 cached-economics tuner picks leaf=275 (0.69 GiB, 57.6 s build)
   rather than the kernel-tuned leaf=21 — the answer sat just past the old cap.

2. **ILU pattern cap (p1-table-R6).** `max_pattern_entries=2048*N`, set at R3,
   undershot R6 by 0.6% (needs 2,059.5 entries/row). Raised to 8192*N across all
   seven benchmark scripts. See the new ILU scaling finding in `ledger.md` +
   `ilu_scaling/` — the pattern grows ~N^1.5, which is a Phase-1 result in its
   own right, not just a cap constant. The guard now reports the total required
   rather than the subtotal at which it tripped.

3. **Cluster filesystem EACCES (p2b R1, R2; p1-tune-R7).** NOT a repo bug.
   `open(".../p2b/FLOWPanel.jl/.CondaPkg/lock") — permission denied` for the p2b
   jobs, and `open("~/.julia/compiled/v1.12/Pkg/*.ji.pidfile") — permission
   denied` for the R7 tune. FLOWPanel hard-depends on PythonPlot
   (`Project.toml:30`), so `import FLOWPanel` always runs PythonCall's `__init__`
   -> `CondaPkg.resolve()` even though the benchmark path never plots;
   `JULIA_PKG_PRECOMPILE_AUTO=0` does not prevent this, as it is module init and
   not precompilation. Ryan is clearing the stale lock state cluster-side
   (2026-08-22 ruling); submit scripts deliberately left untouched. Same depot as
   the in-flight 018 jobs, so worth watching there too.

Verification: 347 unit tests pass (`runtests_unit_solver.jl` +
`runtests_unit_fmm.jl`), including the ILU guard and 021 rigid-motion testsets;
R1 Phase-2b driver smoke passes all three sections and writes both CSVs.

Still open: resubmission (waiting on the EACCES cleanup), the four never-submitted
p2b jobs, and reconciling the item file's `## Current status` / phase-gate table,
which still read "Phase 2 NOT STARTED" and are stale since 08-13.

### 2026-08-21 — Wake-on FGS plateau diagnosed and fixed (H3 lifecycle ordering)

Dense-LU arbitration on the first particle-carrying R1 state identified FGS
as the incorrect path (FGS solution error 2.073e-3 and true relative residual
1.941e-2; fresh Krylov 5.337e-9 and 9.542e-10). The underlying H3 bug was a
one-step control-point lag: `simulate!` transformed persistent solver target
buffers before recalculating moved control points. Reordering geometry refresh
before `transform_body_solvers!` removes the jump. Fixed 0°–80° divergence is
flat at 1.740e-7–1.893e-7; dense wake-on agreement is 1.879e-7. Added the
triaxial panel rigid-motion gate and simulation ordering regression. Details:
`phase_02_single_step_benchmarks.md` 2026-08-21 entry.

### 2026-08-20 — Rigid-motion tree/cache reuse executed; FGS unsteady staleness confirmed + fixed

`rigid_motion_tree_reuse_item.md` executed (Ryan's go via launch prompt).
The pre-registered FGS staleness discriminator CONFIRMED the bug: per-step μ
divergence vs a fresh-Krylov reference grows ~10⁶× with rotation angle
(1.7e-7 at 0° → 21% at 80°, R1, 8 steps) — pre-fix unsteady FGS results are
untrustworthy. FastMultipole gained `transform_tree!`/`transform_plan!`/
`transform_solver!` (commits eea944d/087bf4a/d714544; note Ryan's concurrent
645cc96 swept in the transform_tree! src half) with exact scalar
nearfield-cache persistence under rigid motion and loud v1 refusals for
direction-carrying outputs; FLOWPanel (uncommitted) mirrors per-step rigid
kinematics deltas into persistent solver state inside `simulate!` and adds
`KrylovSolver(persistent_plan=true)` (deferred nearfield-cache commit 4).
All suites green. Details + tables: `phase_02_single_step_benchmarks.md` Log
2026-08-20 (tree-reuse entry) and the item file's "Execution results".

### 2026-08-20 — Phase-2b HPC chains launched

R1–R4 per-rung afterok chains submitted (p2b-nearfield → p2-table-nf:
13242659–13242666, zen3 exclusive, multi 64T) from a SEPARATE cluster
checkout `~/projects/p2b/` (Phase-1 chains still in flight — main checkout
untouched). Full detail + five ops judgment calls (p2b checkout, FLOWVPM
symlink, FM test/-dir ship requirement, JULIA_PKG_PRECOMPILE_AUTO=0, R5–R7
deferred): `phase_02_single_step_benchmarks.md` Log 2026-08-20.

### 2026-08-19 (evening) — nfcache benchmark configs + rigid-motion staging

Per Ryan: three nfcache variants added to the phase2 driver
(gmres/ilu/fgmres_fgs; Jacobi skipped), shared cached-tuned knobs per rung
via new `tune_cached.csv` plumbing; cache build counts inside per-step
(per-solve state, matching the Krylov tree-rebuild convention). R1 smoke:
**all three ≈33–34 s cold, BUILD-DOMINATED (~94% = the 31.7 s per-solve
cache build; preconditioner worth <1 s)** — points at the newly STAGED
`rigid_motion_tree_reuse_item.md` (transform trees under rigid motion,
shared FmmPlan+FGS machinery, cross-timestep cache persistence) as the real
lever. That item also carries a **CORRECTNESS FLAG: FGS unsteady staleness**
(ctor-frozen trees under a rotating body — far field silently degrades with
rotation angle; steady campaign unaffected; verify before any Phase-3
unsteady FGS row). Detail: `phase_02_single_step_benchmarks.md` Log
2026-08-19 "(evening)" entry.

### 2026-08-19 (later) — Cached-path autotuning + build-time guard (Ryan feedback round)

Ryan rulings: FmmPlan cache ownership confirmed (Tree-field rejected);
tune=true refusal superseded. Landed FM 0ef4e83
(`estimate_nearfield_cache` + `max_build_time` guard — estimate from one
kernel sample, never time the build to discover it) + 204188a (tune on the
cached path: provided cache tunes expansion_order ONLY; `tune_fmm
tune_nearfield_cache=true` builds throwaway caches per trial, build cost
excluded, leaf changes stop at the bytes/build-time caps). R1 smoke:
cached-economics knobs **p=14/MAC0.5/leaf=342** vs kernel-tuned
p17/MAC0.5/leaf21 — leaf optimum ~16× larger, cache 862 MB / est build 34 s,
uncapped. All suites green 4T. Detail + rulings:
`phase_02_single_step_benchmarks.md` Log 2026-08-19 "(later)" entry.

### 2026-08-19 — Near-field influence-matrix cache IMPLEMENTED (scheduled autonomous session)

Executed the staged plan: FastMultipole commits 72d3f3d
(`NearfieldInfluenceCache` + builder + deterministic owner-partitioned
`nearfield_matvec!` + standalone `direct!` form + tests) and 02f071c
(`fmm!`/`FmmPlan` integration); FLOWPanel `cache_nearfield` KrylovSolver
surface implemented but UNCOMMITTED per the launch constraint. V0
(`_induced_wake` linearity) PASSED (8.6e-16 on the shedding diamond, 4.3e-16
at R1 scale) — no refusal path needed. R1 local smoke A/B: isolated
near-field **265×**, cold krylov_gmres solve 131→23.7 s (**5.5×**, same 83
iters, certified BC 8.47e-9 both, solution shift 3.5e-15), build 11.9 s /
272 MB / break-even 8.3 applies. All suites green at 4T. Commit 4
(`persistent_plan` cross-solve reuse) deferred pending Ryan. Full detail +
9 flagged judgment calls: `phase_02_single_step_benchmarks.md` Log
2026-08-19 (IMPLEMENTED entry); rows in `ledger.md`; CSV
`benchmark/results/phase2/multi/nearfield_cache_ab.csv`.

### 2026-08-19 — Near-field influence-matrix cache lever STAGED (Ryan)

Plan-only session: new Stage 2b lever (packed dense near-field blocks keyed
to a Tree, general over `direct!`, timestep-reusable; performance-first)
staged as a self-contained implementation plan at
`nearfield_matrix_cache_plan.md`; pointer + design summary in
`phase_02_single_step_benchmarks.md` Log 2026-08-19. No code touched.

### 2026-08-18 (eve) — Phase 2 prep executed concurrent with Phase 1 HPC chains

Ryan authorized concurrent Phase 2 prep at full scope (harness + ALL staged
src levers — the selection constitutes the sign-off Part C2/C3 of the FGS
performance plan were gated on) + a corrected unsteady driver via RHPC
setup-only include. Everything landed and smoke-validated locally; full
detail in `phase_02_single_step_benchmarks.md` Log 2026-08-18. Headlines:
Phase 2a componentized driver + unsteady driver + p2 Slurm launchers (not
submitted) + profiling harness + analysis/TikZ pipeline all working at R1;
`FmmPlan`/`cache_tree` per-solve Krylov tree cache (bitwise-verified,
235/235 solver units); `sweep_order=:colored` FGS sweeps with the batching
theorem proven bitwise (2216/2216) — **but multicolor GS measurably diverges
on a synthetic case where lexicographic converges, so campaign convergence
A/B is a hard adoption gate**; scheduling A/B per Ryan: default
`Threads.@threads` beats `:static` ~4–10%, default kept. Nothing rsynced to
the cluster; Phase 1 evidence untouched (one pure code-motion extraction
from table.jl, smoke-verified). Legacy FGS `reverse_pass` found to be
never-actually-reversed (counter-vs-branch quirk) — preserved bit-identical,
flagged for Ryan. FGS per-step iteration capture in `simulate!` flagged as a
Phase 3 prerequisite. FastMultipole full suite at session end (4T, incl. the
two new test files): exit 0, 73 testset summaries, zero failures/errors.

### 2026-08-17 — Leaf-LU caching implemented and local R1 cache-only A/B complete

Implemented the deferred D3 improvement across FastMultipole and FLOWPanel.
Original self-influence blocks remain untouched for residual matvecs; checked
LU factors live in one distinct contiguous buffer and are reused through
three-argument `ldiv!`. Default is cached; `cache_leaf_lu=false` is the
low-memory/control path. Wrapper keywords and metadata are covered by tests.

FastMultipole's full package suite passes (including 153 cache
invariants/trajectory assertions); FLOWPanel's full regression suite also
passes (focused solver 217/217 and history suites 43/43). Local R1
4T/BLAS-1 min-of-5 at frozen Phase-1 settings preserved iterations and
certified BC errors exactly while improving leaf pass 28.1x, FGS 2.63x, and
FGMRES+FGS 1.49x. The factor cache is 10.16 MiB and built in ~10--11 ms.
Raw evidence is under `benchmark/results/phase2/leaf_lu_cache/multi/`.

No HPC results or cache-enabled production retuning are claimed: Phase 1 has
not frozen R4--R7 settings, so those prerequisites remain pending. No notebook
entry was written.

### 2026-08-17 — FGS nondeterminism attributed to M2L target overlap; fixed and verified

Executed the staged determinism plan on R1 with a true
`{Julia 1T,4T} x {BLAS 1T,4T}` subprocess matrix. The pre-fix 4T/BLAS-1T
replay had exact constructors, tree geometry/topology, lists, packed matrices,
upward/L2L/L2B passes, and representative leaf LU solves, but fixed-input M2L
differed in every repeat. Root cause: FGS fed the raw non-owner-major M2L list
to `assign_m2l!`, whose partitioning requires contiguous target ownership;
separate tasks could `+=` into one target expansion. Canonical stable counting
sorts by source then target now run before matrices/index maps. A premise-guarded
FastMultipole regression requires multiple leaves, nonempty M2L, and exact
repeated cold residual/strength histories.

Post-fix, all four thread coordinates pass five repeated cold solves and three
fresh constructors bit-for-bit. Selected 4T/BLAS-1T fixed-ten minimum changed
0.886844 -> 0.894550 s (+0.87%, gate PASS). FastMultipole counting sort was
3.14x faster than generic lexicographic `sort!` on the 1,760-pair R1 list.
Corrected post-crossing geometric-mean margins pass all 18 R1--R3 checks; R3
tau 1e-6 correctly forced a staircase extension before verification. Evidence
and full attribution table: `results/fgs_determinism/summary.md`. FastMultipole
solve tests 582,335/582,335; FLOWPanel history 15/15; solver units 200/200.
Stopped before LU caching or any performance redesign. No notebook entry.

### 2026-08-13 (later) — First gate review DECLINED on evidence-record consistency; remediated + rerun

The clear-context gate review (subagent, authorized by Ryan) **verified all W6 substance**
(code correctness incl. permutation/equilibration algebra, tests green, tables ≡ CSVs,
acceptance claims recomputed) but **DECLINED** on six findings — chief (Medium): the CSV
of record labeled the cold benchmark solves `warmstart=true`, because `run_config!`
captured `solver_settings` after the `simulate!` leg's `body_solvers_for` mutation
(ruling 10 makes that column authoritative; the timed solves were traced genuinely cold).
Plus: stale control-doc header ("STAGED — not started"), stale phase_00 header/W6 section,
a wrong setup-breakdown claim ("tree+lists ~0.5 s"; actual 0.15 s + ~0.4 s untracked ctor
overhead), and the 08-12 log entry's misordered approve/rescind paragraphs. One
informational note (untested missing-diagonal top-up path — dead under Barba
`self_induced=true`, guarded loudly; no action).

Remediation, same day: harness captures `solver_settings` before the sim leg (flip
recorded in `notes`); all four doc findings fixed; **both smoke modes rerun** (agent
judgment: reviewer deemed rerun optional, but clean CSVs beat an annotation for a
publishable campaign). Per Ryan (back online): the reruns ran **concurrently** (5
threads total) with **multi at 4 threads** — accepted as CSVs of record with caveats
noted in the phase_00 tables; iteration counts/residuals reproduced exactly, wall times
carry contention (e.g. plain GMRES single 281 s vs 223 s quiet). Ryan rulings: proceed
into Phase 1 on the re-review's approval; no notebook entry yet.

**Second gate review: APPROVED** (same day). All six findings verified fixed/no-action;
both 08-13 tables cross-checked cell-for-cell against the new CSVs; acceptance holds in
both modes; unit solver suite green. Two informational notes, no action: the sim-leg
note in `notes` is appended unconditionally (cosmetic; `solver_settings` authoritative
and correct), and "lists 0.04 s" matches multi (single shows 0.058 s, within "~").
Phase 0 gate: clear-context approved; Ryan's re-approval checkbox pending his review.
**Phase 1 begins** per Ryan's 2026-08-13 ruling.

### 2026-08-13 — W6 complete: ILU implemented (Ryan), reviewed, smoke-run PASS both modes

Ryan implemented W6 directly on the evening of 2026-08-12 — `phase_00_ilu_plan.md`'s
Stages A–C collapsed into one sitting; the plan file is marked superseded. Design
decisions made by the implementation: pattern = FastMultipole **Barba direct interaction
list** (dedicated ctor-time tree, not the staged radius-cutoff fallback); factorization =
**`ILUZero.jl`** (new dep). This session (agent): reviewed the implementation for
correctness and performance — verdict CORRECT (assembly ≡ `_G!` semantics; permutation +
equilibration algebra checked; ILUZero internals audited at source). Fixes applied:
ctor deduplicated onto `_ilu_direct_pattern` (helper now returns timings); assembly
threaded over direct-list pairs (mirrors `_G!`; `induced` is pure). Two crash bugs in the
harness's new code fixed at launch time: `\"` inside `$()` interpolation (1.12 parser
rejects) and empty-generator `sum` in the `setup_total` fallback. Unit + full suites
green.

Smoke runs relaunched from scratch (both modes, quiet machine; the harness truncates
`runs.csv` per launch, so the 2026-08-12 corrected-rerun CSVs are replaced on disk —
carried-over configs reproduce to ~1–2%, and the old numbers persist in
`phase_00_availability.md`'s tables). **W6 acceptance PASS**: ILU-GMRES 290→15 iters,
223.4→11.51 s (single) and 55.2→2.86 s (multi) vs plain GMRES; setup ~1.0–1.2 s,
break-even <0.02 solves; displaces FGMRES+FGS as best iterative config at smoke scale.
Full tables + stats in `phase_00_availability.md` Log 2026-08-13. Side observation: multi
FGS hit 500 iterations this run (prior rerun: 267) — thread-order nondeterminism,
consistent with the pre-registered Phase 2 finding.

Phase 0 status: technically complete incl. W6; pending clear-context review + Ryan's
gate re-approval. Ryan authorized (2026-08-12) autonomous continuation incl. a subagent
clear-context gate review, then Phase 1.

### 2026-08-12 — Rulings 11/12; clear-context review findings; remediation

Morning: Ryan's GMRES directives adopted as ruling 11 (right-preconditioning honesty —
landed in code, Jacobi now routed `N=`; Phase-1-calibrated tolerances; FMM-floor
measurement; preconditioner exploration with sparse near-field ILU as priority — theory doc
`theory_nearfield_ilu.md` written) and ruling 12 (warmstart benefit = headline deliverable;
Phase 2 strictly cold, Phase 3 owns the warmstart matrix incl. break-even step counts).

Afternoon: a clear-context review declined Phase 0 approval with 5 findings. Verified
assessment — all confirmed, one scope correction:

1. (High, CONFIRMED) **Warm-state FGS/block-GS evidence**: the smoke timing loop never
   reset `body.strength`; FGS seeds from current strengths and breaks on the residual check
   *before* sweeping, so every timed standalone-FGS rep after warmup was a no-op residual
   check. The "FGS converges in ONE outer iteration / 0.047 s" headline was an artifact and
   had propagated to the phase log, control doc, INDEX, and agent memory. Block-GS's
   1-outer history was the same warm-entry artifact. **Scope correction**: config 1f
   (FGMRES+FGS) was NOT warm-biased — Krylov never seeds from `body.strength` and
   `FGSPreconditioner.ldiv!` zeroes strengths every apply (the linearity contract) — so
   "1f ≫ plain GMRES, 26 vs 290 iterations" survives; 1c/1d timings were also genuinely
   cold per rep.
2. (High, CONFIRMED) Config 1d CSVs predate the ruling-11a right-routing change; the
   left-routed "Jacobi hurts, 362 iters" conclusion must be re-derived.
3. (Medium, CONFIRMED) CSV provenance gaps (`julia_version`, tolerances, FMM knobs,
   `t_rhs`, banner stdout-only).
4. (Medium, CONFIRMED) `Base.summarysize(solver)` had inconsistent object boundaries
   (Krylov holds the body; Backslash doesn't).
5. (Medium, CONFIRMED) FastMultipole callback was an uncommitted dev-checkout edit.

Remediation landed same day: cold-state reset before every timed/history/alloc solve
(`min_of_k` gained an untimed `setup!`; FGS cold tripwire assert); provenance columns
`fm_commit, julia_version, solver_settings, backend_settings` + populated `t_rhs` +
`banner.txt` (schema updated in `decision_rules.md`); memory metric finalized as
`solver_state_bytes` (summarysize minus referenced bodies — comparable boundary);
FastMultipole callback committed as `5adde3b` (dev checkout, local commit) + portable diff
at `patches/fastmultipole_callback.patch`.

**Rerun complete same day** (single mode clean; the first multi rerun was discarded —
Ryan had a concurrent job on the machine — and repeated on a quiet machine). Corrected
tables + headline findings in `phase_00_availability.md` Log 2026-08-12; all four doc
surfaces replaced. Biggest corrections: cold FGS = 35.2 s / 204 iterations (not 1/0.047 s);
"Jacobi hurts" re-confirmed honestly; NEW measured finding — FGS multi-thread solve is
2.2× slower with more iterations (267 vs 204; threaded sweep-order convergence
degradation), now isolated to the solve itself. Untouched configs reproduced to ~1%.
Phase 0 back to technically complete. **Ryan approved the Phase 0 gate later the same day**
(recorded in the control-doc decision log and gate table).

**Evening: that Phase 0 approval was RESCINDED by Ryan** — ILU pulled forward from
Phase 2b into Phase 0 as **W6** (develop now: theory doc → implementation design →
implementation/testing; roster gains config (g) GMRES+ILU). Self-contained execution
plan for a fresh agent staged at `phase_00_ilu_plan.md`; W1–W5 evidence and remediation
stand unchanged. Next action: execute W6 — resolved in the 2026-08-13 entry above.

### 2026-08-11 — Phase 0 started; W1 dropped; D2/D3 ruled

Ryan gave the Phase 0 go-ahead. Planning exploration (2026-08-07) verified the solver
inventory and produced three corrections now folded into the control doc:

1. The rotor pipeline is **not** blocked by the Kutta gate — `rotor_hover.jl` uses the legacy
   linear Kutta condition; the `FLOWPanel_kutta.jl:495-503` hard error only fires for the
   opt-in 015 pressure-continuity closure (Dirichlet paired bodies).
2. `BackslashCoupled` availability bug: ctor installs a dummy identity `lu!`
   (`FLOWPanel_solver.jl:989`) awaiting `update_G=true`, but `solve_formulation!` never passes
   `update_G` (default false) — inside `simulate!`, roster config 1a silently solves against
   the identity. Fixed in Phase 0 (WP-A).
3. `Backslash` stale-`Glu` corruption: ctor's `lu!(G)` aliases `G` as `Glu.factors`; a later
   `_solve!(update_G=true)` factorizes into a *local* LU, leaving `solver.Glu` with new
   factors + old pivots — garbage for direct consumers (Kutta Route A, GreenReconstruction,
   DirectWakePotential). Fixed in Phase 0 (WP-D).

**Ryan rulings (recorded in control-doc decision log):**

- **W1 dropped.** The campaign runs the legacy linear Kutta condition throughout; upgrade to
  the 015 pressure-continuity closure only if accuracy/convergence later demand it. The
  kutta.jl hard error stays. Phase 0 W1 is rescoped to the two bug fixes above.
- **D2 — FGS convergence history via FastMultipole callback**: add a non-breaking
  `callback=nothing` kwarg to `FastGaussSeidel` `solve!` in the dev checkout
  (`../FastMultipole/src/solve.jl`), invoked once per outer iteration with
  (iteration, residual) so FLOWPanel records exact production-loop timestamps (ruling 4).
  Cross-repo edit — noted here per the repo boundary.
- **D3 — leaf-LU caching deferred to Phase 2b.** Discovery: FastGaussSeidel leaf blocks are
  dense and re-factorized every sweep (`mat \ rhs`, FastMultipole `solve.jl:934`) — the
  control doc's "leaf LUs built once" was wrong (corrected). This inflates FGSSolver and every
  FGS-preconditioner apply; Phase 1's 1e/1f calibration numbers must be read with that caveat.

**Archived W1/Kutta design (revivable as its own item if the closure is ever needed):**
option (a) solver-generic inner solve — `Backslash` mutable + `Glu` write-back;
`KuttaRuntime{TS<:Backslash}` widened to `AbstractSolver`; `_kutta_inner_solve!(rt, body)`
replacing the direct `_kutta_ldiv!` (Backslash = ldiv path; Krylov = production `_solve!` with
forced warmstart from `rt.current.mu`, essential at ~10–50 trials/step; FGS = seed
`strength[:,2]`, finite-only status); validator relaxed to accept
Backslash|KrylovSolver|FGSSolver, still rejecting coupled/tuple≠1; Route B for matrix-free =
`_assemble_W!` + version bump (no G rebuild); FGS × TEAnchoredAttachment needs per-step
FastGaussSeidel rebuild (stale leaf matrices) — the contested piece. Correctness test design:
unsteady wing, Krylov/FGS vs Backslash at Phase 1 provisional thresholds.

Implementation plan for this session: `/Users/ryan/.claude/plans/work-on-brainstorm-item-bright-lobster.md`
(work packages A→B→C→E→F→D). Also noted: single-body tuple block-GS always runs a second
identical solve to observe Δ=0 — config 1b timing must record solve counts (harness notes).

**Same-day completion:** all work packages landed and Phase 0 exit criteria met (results and
per-config smoke tables in `phase_00_availability.md` Log, 2026-08-11 entry). Beyond the
planned scope, the smoke run exposed and fixed a third availability bug: FastMultipole's
`strength_to_value` had no overload for `RigidWakeBody{<:Any,1}`, so FGSSolver had never
actually worked on the uncapped Neumann rotor body type (the commented FGSSolver block in
`rotor_hover.jl` was aspirational). Fixed in `src/FLOWPanel_liftingbody.jl`, test-gated.
Cross-repo edits to `../FastMultipole` this session: the D2 `callback` kwarg only
(`src/solve.jl`); FastMultipole solve tests pass. Removed stdout `@show workspace.stats` from
Krylov solves (aligned with judge-from-CSVs). Full FLOWPanel suite green (45 testsets).
Phase 0 gate approval + notebook entry pend Ryan.

### 2026-08-07 — Item opened and staged

Item created from Ryan's 2026-08-06 brief plus a solver-subsystem inventory (control doc
"Solver matrix"). Key discoveries baked into Phase 0: Kutta closure hard-errors on
Krylov/FGS; FGS-as-preconditioner and Krylov `x0` don't exist; no convergence-history capture
anywhere. Ryan's decisions recorded in the control-doc decision log (benchmark scope, HPC,
dual threading modes, NACA 0015 wing). No code touched. Next: Phase 0 on go-ahead.
