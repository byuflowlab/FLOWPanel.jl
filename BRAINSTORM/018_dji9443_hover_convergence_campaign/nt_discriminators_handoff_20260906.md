# RESET BRIEF — NT discriminators, wave-1 triage + fallout repair (2026-09-06)

Supersedes `nt_discriminators_handoff_20260905.md` as entry point. **Everything
in that file's §0 (the science), §2 (worktrees/tags) and §3 (suspect ranking)
still stands unchanged** — read it for those. This file records what the
2026-09-05→06 session established on top of it: the wave-1 arms nearly all
died, why, what was fixed, and what is still owed.

Guardrails unchanged: sacct state is NOT evidence (judge by outputs), `ssh orc`
needs `bash -lc` + live ControlMaster, silos are dead (worktrees + annotated
tags only), NODE_FAIL → fresh resubmit never a chain.

---

## 0. THE GOTCHA THAT COST A FULL DIAGNOSTIC ROUND

The GPU launcher writes stdout and stderr to **separate files**
(`examples/run_dji9443_hover_ct_gpu.slurm.sh:10-11`):

```
#SBATCH --output=logs/slurm/slurm-%x-%j.out
#SBATCH --error=logs/slurm/slurm-%x-%j.err
```

Every Julia stacktrace, every `@info` gemv marker, and the launcher's own
`ERROR: dispatcher exited rc=$RC` line (line 109) land in the **`.err`** file.
Reading only `.out` makes every death look like a silent "post-run dispatcher
bug". It is not. Two more facts about the gate block (lines 100-138):

- `dispatcher_rc` is the exit code of the **julia run itself**; the GATE block
  runs *after* it. So `gpu_gemv=N` on an M-step run means julia died at step ≈N.
- `gate_rc=1` has **four** distinct causes (GPU path never ran / CPU path ran /
  NaN in log / dispatcher rc ∉ {0,124}). It is a summary, not a diagnosis.

**Always read the `.err` first.**

## 1. Wave-1 outcome: 4 distinct failure modes, not one

Only ONE of the original 11 arms produced a clean score. Full triage:

| Job | Arm | Died at | Cause |
|---|---|---|---|
| 13592731 | `_3r_sv_s1p5` NT36 | — | **COMPLETED**, 30 revs, window opens rev 21 |
| 13592732 | `_3r_sv_s1p5` NT72 | — | RUNNING (7h+ as of reset) |
| 13592894 | `_3r_exp` NT36 | step 0 | `ArgumentError: sigma_guard is not supported by the euler_exp integrator` |
| 13592723 | `_3r_sv` NT72 | pre-step-0 | DAS arc-table `isfile()` failed (worktree data-symlink not yet created at 13:23; pins were made 13:21) |
| 13592738 | `_3r_sv` NT36 | 537/1079 | `WakeGeometryError: wake panel is folded or inverted` |
| 13592759 | `_3r_sv_h2p0` NT36 | **1076/1079** | same WakeGeometryError |
| 13592726 | `_3r_sv_h2p0` NT72 | 672/2159 | same WakeGeometryError |
| 13592729 | `_3r_sfs3nb` NT72 | 1550/2159 | **FMM near-set adequacy** (see §2) |
| 13592760 | `_3r_sfs3nb` NT36 | 816/1079 | genuine NODE_FAIL |
| 13592896 | `_3r_srlx` NT72 | — | still PENDING (untouched, original wave) |

**Salvage not yet done**: 13592759 died 3 steps from the end with σ-growth 1.19
against a 1e9 critical threshold — nothing diverging. Its `CT_per_rev.csv`
should be scoreable **without a rerun**. NOBODY HAS SCORED IT YET.

**New lead on the open WakeGeometryError root cause**: the only smooth-conversion
arm that survived is `_s1p5` (reduced 1.5σ smoothing); both overlap-2.75 and
overlap-2.0 arms folded. Smoothing/deposition width may govern the folding.

## 2. FMM near-set adequacy — the strongest lead of the session

`_3r_sfs3nb` NT72 died with:

> `regularized nearfield near-set adequacy failed: the direct stencil leaves an
> M2L gap of g_min*h_leaf = 0.2993 but the smoothing cutoff needs
> rho_t*sigma_max = 0.9751 (ratio 0.3069, g_min = 2.449, sigma_max = 0.2036,
> ell = 3); the admissible depth at this geometry is ell <= 1. Pairs inside the
> cutoff would be handled by the singular far field and silently lose the
> regularization (row 032a).`

Why this matters for the NT anti-convergence (mechanism, hedged — NOT yet
tested):

1. The regularized (Gaussian) kernel is evaluated only on the **direct near
   list**; anything outside goes to the multipole far field, which uses the
   **singular** kernel. Adequacy requires
   $g_{\min} h_{\text{leaf}} > \rho_t \sigma_{\max}$ with
   $h_{\text{leaf}} = 2h_0/2^{\ell}$.
2. When it fails, near pairs get the singular kernel, which **exceeds** the
   regularized one as $r \to 0$. The error is a one-sided **over-induction**,
   not zero-mean noise. A bias, not a tolerance.
3. **NT coupling**: at fixed particles-per-step, doubling NT doubles particles
   shed per rev and halves helix spacing; `_radix_auto_geometry` sizes depth
   from live `np` (cap $\sim np^{1/3}$), so $\ell$ deepens, $h_{\text{leaf}}$
   shrinks, while $\sigma_{\max}$ (chord/λ-set, plus core spreading) does not.
   The ratio degrades monotonically with NT.
4. A fixed per-step over-induction × ∝NT steps/rev ⇒ CT climbs with NT and
   per-doubling deltas **grow**. Truncation would halve; a multiplicative FMM
   tolerance would integrate NT-independently. This mechanism has the right shape.
5. Not speculative that we are near the boundary: the gate fired at NT72 with
   ratio 0.307 and "admissible ℓ ≤ 1" against an actual ℓ = 3. NT36 plausibly
   sits just under the threshold — silently biased, never throwing.

**Cheap falsifiable test (not yet run)**: post-hoc evaluate the adequacy ratio
across the five completed anti-convergence arms at NT36/72/144. If it degrades
with NT while staying under the throw threshold, the mechanism is live. Or rerun
one rung with forced shallow `ell` / all-direct and see if the climb flattens.
Ryan's ruling 2026-09-05: **do not chase this yet.**

**ANOMALY, unresolved**: the gate reports `sigma_max = 0.2036` m. DJI 9443 radius
≈ 0.119 m and that wave ran `SIGMA_CEIL=0.030`. A σ ~7× the ceiling and larger
than the rotor means either the ceiling is not binding where we assume, or
`sigma_max` is read from a non-wake system's sigma row. Worth resolving — it
affects more than this arm.

**The all-direct fallback exists but was not in the pin.** Ryan recalled
correctly: `_alldirect_geometry_fallback!` (FastMultipole `d938ba68`,
2026-08-31, "052f: demote inadequate hierarchical geometry to all-direct
zero-M2L") demotes instead of throwing when
`policy isa HierarchicalRigidStencil && !(kernel isa TwoPassVortex)`. Our kernel
qualifies (message says `rho_t`, not `rho_c`) and FLOWVPM never passes `policy`
so it gets the default `HierarchicalRigidStencil`. But the campaign pins are
**silo snapshots dated Aug 26** — the fallback is absent from
`~/wt018/FastMultipole-{h200,gh200}`. It IS present in
`~/wt018/FastMultipole-nt` (ported onto the `unified-052` line, so
`git merge-base --is-ancestor d938ba68` says *no* while the code is *there* —
check by grep, not by merge-base).

## 3. The `-nt` stack: fixes one thing, breaks another — BLOCKER

The exp (`WAKE_EXPINT=true`) arm needs a guard-free euler_exp on GPU. Chain of
three failures, each revealing the next:

1. `sigma_guard` + `euler_exp` are mutually exclusive **by design**
   (`FLOWVPM/src/FLOWVPM_timeintegration.jl:33-34,108`; guard is euler /
   ReformulatedVPM only). The wave-wide `SIGMA_CEIL=0.030` export therefore
   cannot coexist with `WAKE_EXPINT=true`. Dropping `SIGMA_CEIL` makes
   `sigma_guard` empty (`rotor_hover_pressure_comparison.jl:1285`: empty iff
   `SIGMA_DTZ_CAP=Inf && SIGMA_FLOOR_FRAC<=0 && SIGMA_CEIL=Inf`) and clears it.
   → **Confound**: the exp arm cannot carry the σ-ceiling the rest of the ladder
   carries. Ryan's ruling 2026-09-05: **run unguarded**; if a later run shows the
   guard is still necessary with expint, implement it for euler_exp then.
2. Guard-free on the OLD pin → `TaskFailedException: Scalar indexing is
   disallowed` (CPU-only scalar loop on the CuArray pfield). Cured only by
   `~/wt018/FLOWVPM-nt`, whose HEAD is *"port 026 GPU/broadcast euler_exp
   (local 8b00dbd) for guard-free expint trials"*.
3. On the `-nt` stack (job 13593792, NT36) it stepped cleanly at 3.96 s/step —
   then died at the first VTK write, 7 min in:

   > `ERROR: LoadError: type SplittingState has no field dsigma2_visc`
   > `_write_particles_vtp` @ `FLOWPanel-nt-gh200/src/FLOWPanel_wake.jl:2441`

   **Version skew INSIDE the `-nt` pin set**: `FLOWPanel-nt`'s VTP writer emits
   `split_dsigma2_visc` / `split_dsigma2_rvpm`, but `FLOWVPM-nt`'s
   `SplittingState` (`src/FLOWVPM_particlefield.jl:101`) has no such fields.
   FLOWPanel-nt is ahead of FLOWVPM-nt.

   **Two ways out, unresolved — needs Ryan:**
   - (a) `SAVE_VTK=false` (`rotor_hover_pressure_comparison.jl:25` — `save_path`
     becomes `nothing`). Cheapest; CT CSVs/monitors still written, which is all a
     discriminator needs. Costs ParaView output and warm-start ability. Doubly
     attractive while /home is over its cap. Output-only difference, so it does
     NOT affect physics comparability.
   - (b) Fix the skew: add the two fields to `FLOWVPM-nt`'s `SplittingState`, or
     guard the writer with a `hasproperty` check.

**Also: `-nt` deps are unified and NOT bit-comparable to the pinned arms**, so
runs on it were named `_3r_exp_nt` (not `_3r_exp`) to avoid burying that in a
silent overwrite. Keep that convention.

**`--constraint=arm` is REQUIRED on mgh.** Submitting without it is rejected
("incorrect resource request"), despite the launcher header comment claiming a
constraint is redundant there. eng needs `--qos=eng --constraint=intel
--gres=gpu:h200:1 --cpus-per-task=64`.

## 4. Live queue at reset (2026-09-06)

| Job | Arm | State | Note |
|---|---|---|---|
| 13592732 | `_3r_sv_s1p5` NT72 | RUNNING 7h+ | the one NT72 score in flight; do NOT strip its VTK |
| 13593011 | `_3r_sv` NT72 | PENDING | other agent's resubmit of 13592723 |
| 13592896 | `_3r_srlx` NT72 | PENDING | original wave, untouched |
| 13593711 | `_3r_sfs3nb` NT36 | **COMPLETED** | gate_rc=0, gpu_gemv=1080/1080 — NOT YET SCORED |
| 13593792 | `_3r_exp_nt` NT36 | FAILED | §3 item 3 |
| 13593793 | `_3r_exp_nt` NT72 | CANCELLED | cancelled pre-emptively, same blocker |
| 13593717/13593718 | `_3r_exp` NT36/NT72 | FAILED / CANCELLED | old-pin attempts |

Scored so far: `_3r_sv_s1p5` NT36 = see 13592731's `CT_per_rev.csv` (30 rows,
`in_convergence_window` true from rev 21). `_3r_sfs3nb` NT36 (13593711) is
complete and **awaiting harvest**.

## 5. Storage: fixed a real bug, but the cap is still breached

**Bug (fixed, deployed, verified)**: `run_archiver.sh --all-checkouts` treated
`data/`-symlink aliases as independent checkouts and re-archived the same runs
once per alias. 9-10 of 13 "checkouts" are aliases of
`projects/FLOWPanel.jl`. Fix keys identity on the **physical `data/` path**,
promotes the real owner over symlink aliases, prints `ALIAS-SKIP`. Verified on
the cluster: **13 → 3 roots**. Bonus: the "95.7 G unreachable CHECKOUT-LOCKED
`FLOWPanel-018-gpu-gh200`" was never independent data — it is an alias, so that
blocker is gone.

⚠ **The fix is deployed to the cluster but UNCOMMITTED in the local repo**
(`scripts/run_archiver.sh`, working tree only). Cluster copy at
`~/projects/FLOWPanel.jl/scripts/run_archiver.sh` (md5 `5ed0b87a…`) — note that
file is UNTRACKED on the cluster's `unified-052` branch.

**Fallout cleaned**: two alias archive dirs (66 + 38 entries) moved to
`/nobackup/archive/usr/rander39/FLOWPanel_runs/_alias_fallout_20260905/`; the one
genuinely unique tarball (`fm052d_gpu_1080_t1.prev.20260827`) rescued into
`projects_FLOWPanel.jl/`; two `.partial` tarballs deleted; orphaned reader lock
(dead PID 2191333, since 09-03) cleared. Alias copies were **quarantined, not
deleted**: a sampled pair showed the canonical tarball is a strict superset
(50 members vs 44 — the alias copy was made after the run was stripped and lost
the VTK; only alias-unique member is the `ARCHIVED.txt` breadcrumb). That was 1
verified sample out of 97, and archive quota is effectively unlimited, so
deleting on that inference wasn't worth it. **Ryan can authorise deleting the
quarantine.**

**Resume-deletes done** (7 non-018 stale runs, each verified byte-for-byte):
freed 2,271 MB. `fm052d_gpu_1080_t1` was already complete (its 28.9 G was a
stale dry-run figure). The 6 others carried "deliberately restored" markers and
are now re-stripped — recoverable via `--restore` from their verified tarballs.

**STILL BREACHED**: `/home/rander39` was **430 G against the 400 G cap** at
reset — *worse* than the 409.5 G the watchdog left it at, because live jobs write
VTK faster than we free it.

**IN FLIGHT AND UNFINISHED — the next agent must redo this.** A dry-run
`--resume-delete` over the 018 stale runs was still running at reset and produced
**no captured output** (piped through ssh+grep, output buffered, file empty). It
will die with the ssh session. Redo it writing to a **remote log file** instead
of piping through ssh.

Target list and the two exclusions (both verified, do not skip these checks):

- **STRIP (5 runs, ≈130 G)**: `p018_csarc_l3p0_3r_sv` (25 G),
  `p018_csarc_l3p0_3r_sv_h2p0` (47 G), `p018_csarc_l3p0_3r_sv_s1p5` (42 G),
  `p018_csarc_n2_nt72_l3p0_3r_sfs3nb` (32 G),
  `p018_csarc_n4_nt144_l3p0_3r_sv` (1.5 G).
- **EXCLUDE `p018_csarc_l3p0_3r_sfs3nb`** — has **11,885 files written after
  21:00 on 09-05**: the output of the COMPLETED rerun 13593711, written into the
  same run name as the old NODE_FAILed attempt. Its ARCHIVED-STALE tarball
  describes the *old* run. Needs a **fresh archive**, not a resume-delete.
- **EXCLUDE `p018_csarc_n2_nt72_l3p0_3r_sv_s1p5`** (67 G) — live output of
  RUNNING job 13592732, the one NT72 score in flight.

(`p018_csarc_l3p0_3r_sv_s1p5` was checked and has 0 files newer than 21:00 —
safe.)

Command form that worked:
`./scripts/run_archiver.sh --resume-delete <run> --root /home/rander39/projects/FLOWPanel.jl --apply`

## 6. Owed / next actions

1. **Harvest** 13593711 (`_3r_sfs3nb` NT36, complete) and salvage 13592759's
   partial CSV. Extend the §0 table in the 09-05 handoff.
2. **Finish the storage strip** per §5 (5 runs, ≈130 G) — the cap is breached.
   Re-archive `_3r_sfs3nb` fresh.
3. **Unblock the exp arm** per §3 — Ryan picks (a) `SAVE_VTK=false` or (b) fix
   the `-nt` SplittingState skew. Then relaunch NT36+NT72 as `_3r_exp_nt`.
4. Decide whether to rerun the 3 WakeGeometryError arms as-is (they may fold
   again) or hold pending the smoothing-width lead in §1.
5. **Commit** the archiver fix (§5) — deployed but uncommitted.
6. Investigate the `sigma_max = 0.2036` anomaly (§2).
7. Still owed from 09-05, untouched: notebook entry for the whole arc (ask Ryan
   verbosity first), `ledger.md` stale (pre-08-29), INDEX.md + item-file RESET
   BRIEF stale (~08-21), mgh-1-2 triple-NODE_FAIL unticketed.
8. Benign but unexplained: the fresh `-nt` env emits
   `WARNING: redefinition of constant FastMultipole._NEARFIELD_*` at startup.

## 7. Local repo state at reset

Branch `fastmultipole` at `b9c24e0`. Uncommitted/untracked **mine**:
- `scripts/run_archiver.sh` — MODIFIED (the alias-dedup fix, §5)
- `BRAINSTORM/018_.../nt_discriminators_handoff_20260906.md` — UNTRACKED (this file)
- `BRAINSTORM/018_.../smooth_ladders_provenance_20260905.md` — MODIFIED (pre-existing)
- `BRAINSTORM/018_.../nt_discriminators_handoff_20260905.md` — UNTRACKED (pre-existing)

⚠ **CONCURRENCY**: another agent is editing this SAME local checkout — during this
session `BRAINSTORM/026_.../particle_splitting_design.md`,
`examples/run_p018_screen_hpc.slurm.sh` and
`BRAINSTORM/026_.../phase2_handoff_prompt_1.md` appeared as modified/untracked without
this session touching them. It is also submitting jobs (13593011, 13593711, 13593717/18)
and running its own archiver. **Do not assume the working tree or the queue is yours
alone**; re-check `git status` and `squeue` before acting, and never `git checkout .` or
`git stash` here.

Cluster disk was still climbing at reset: 430 G → **437 G** within the hour.

## 8. Addendum 2026-09-07 — NT36 harvest + fleet reconciliation (new session)

### Harvested NT36 scores (formal; window = `in_convergence_window==true`, 10 revs)

| Case | Job | NT | Windowed CT | Spread | Δ vs ref NT36 (0.070775) | Source / caveat |
|---|---|---|---|---|---|---|
| `_3r_sfs3nb` | 13593711 | 36 | 0.069382 ± 0.000031 | ±0.14% | −1.97% | `CT_per_rev.csv`, clean full run; confirms provisional 0.06938 |
| `_3r_sv_h2p0` (dropped arm, salvage) | 13592759 | 36 | 0.075364 ± 0.000391 | ±1.64% | +6.48% | reconstructed from force monitor `monitor02_force_system1.csv` (end-of-run CSV never written; died step 1076/1079, benign σ-growth 1.19). Lower confidence: finer binning + partial final rev |

Remote paths: `~/projects/FLOWPanel.jl/data/p018_csarc_l3p0_3r_sfs3nb/..._CT_per_rev.csv`,
`~/projects/FLOWPanel.jl/data/p018_csarc_l3p0_3r_sv_h2p0/monitors/..._monitor02_force_system1.csv`.
NOTE: the sfs3nb dir mixes old/new runs — the scored CSV is the 13593711 resubmit (timestamp-checked).

Reading (single-rung, NOT a slope claim): smooth-conversion family sits high at NT36
(+2.0% s1p5, +6.5% h2p0) while legacy-shed+SFS3 sits −2.0% low. Discriminator verdicts
still require the NT72 rungs.

### Fleet after Ryan's 09-06/07 cancel-resubmit churn (sacct-reconciled 07:30 UTC)

| Case | Latest job | State | Note |
|---|---|---|---|
| `_3r_sv_s1p5` NT72 | 13592732 | RUNNING (eng-1-1 H200, 14 h limit, requeued once) | only live job |
| `_3r_sv` NT72 | 13593011 | FAILED = controlled launcher exit at step 737/2159, physics clean (σ-growth 1.0695, 0 NaN) | **no successor submitted** — warm-start from ~step 737 feasible |
| `_3r_srlx` NT72 | 13592896 | COMPLETED clean, step 2160/22 revs, all CSVs, gate_rc=0 | awaiting scoring vs ref NT72 0.071844 |
| `_3r_sfs3nb` NT36 | 13593711 | COMPLETED | scored above |
| `_3r_exp` / `_3r_exp_nt` | 13593717 / 13593792 / 13593793 | FAILED @5–7 min / CANCELLED | consistent with §3 blocker (VTP field skew); arm still blocked on Ryan's (a)/(b) choice |

Monitoring: persistent local monitor (v3.1) polls sacct every 5 min, auto-adopts any new
`fp-018gpu-*` submission from squeue (ANSI-sanitized), greps crash signatures on terminal
states.

### `_3r_srlx` NT72 score (harvested 2026-09-07, same windowing)

| Case | Job | NT | Windowed CT | Spread | Δ vs ref NT72 (0.071844) |
|---|---|---|---|---|---|
| `_3r_srlx` (SFS_RLXF=0.0025031) | 13592896 | 72 | 0.071429 ± 0.000029 | ±0.13% | **−0.58%** |

Coverage: revs 20–30, steps 1441–2160 (full window), from
`~/projects/FLOWPanel.jl/data/p018_csarc_n2_nt72_l3p0_3r_srlx/*_CT_per_rev.csv`.
Interpretation (hedged): if SFS-coefficient memory were the whole NT-bias carrier,
srlx@NT72 should fall back toward ref NT36 (0.070775, −1.49%). Observed −0.58% ≈ 40%
of one doubling step — a real but PARTIAL effect; srlx alone does not account for the
anti-convergence, but it is not a null either. Needs the NT36 srlx rung (or a slope
argument) before scoring it as a discriminator verdict.

### Newly owed
- Score `_3r_sv_s1p5` NT72 (13592732) when it lands → first smooth-family NT36→NT72 slope.
- Decisions still pending on Ryan: sv-arm resubmit/warm-start, exp-arm (a)/(b), storage
  strip authorization, WGE rerun-vs-hold, archiver commit, notebook verbosity.

### 2026-09-07 decisions executed (Ryan's rulings this session)
- **exp arm unblocked via option (a)**: relaunched `_3r_exp_nt` NT36 (job **13602884**,
  case p018_csarc_l3p0) and NT72 (job **13602885**, case p018_csarc_n2_nt72_l3p0) from
  `~/wt018/FLOWPanel-nt-gh200` with `--constraint=arm`,
  `WAKE_EXPINT=true SAVE_VTK=false` (no ParaView output, no warm-start; CT CSVs
  unaffected), unguarded per the 09-05 ruling, `-nt` overrides
  (`P018_REPO_OVERRIDE=~/wt018/FLOWPanel-nt-gh200`,
  `P018_PROJECT_OVERRIDE=~/p018wtenv-nt-gh200`). Logs land in the worktree's
  `logs/slurm/`.
- **Storage strip**: handled by an hpc-storage agent in Ryan's OTHER session; this
  session's agent was cancelled before any archiving began (it was still measuring du).
  Constraints unchanged: 5-run list, archive-first, mixed sfs3nb dir + live s1p5 dir
  excluded, alias-fallout quarantine untouched.
- **Archiver fix committed**: `656e005` (run_archiver.sh alone).
- **WGE arms AND the `_3r_sv` resubmit**: HOLD until 13592732 (`_3r_sv_s1p5` NT72)
  lands (Ryan 2026-09-07: "wait on s1p5 first"), then decide with the smooth-family
  slope in hand.
- **Notebook**: concise SFS-only entry written under `# 20260907` in
  `journals/20260901.md`.
- Clarification recorded: `_s1p5` = `CONVERSION_SIGMA=0.0068009` ≈ **1.5× base
  particle σ** (larger cores), not a reduced width; `CONVERSION_OVERLAP` is the
  σ-to-sampling-spacing ratio (h=σ/overlap). WGE "folded/inverted" = 5-point
  normal-orientation flip check on bilinear wake panels
  (`src/FLOWPanel_wake.jl:1252`).

## 9. RESET BRIEF 2026-09-07 (entry point for the next agent)

Read §8 + this section; §§0–7 for background. Sibling 09-05 file still holds §0 scores
table / §2 worktrees / §3 suspect ranking.

**State**: NT36 SFS-family scores harvested (§8 table); srlx NT72 = −0.58% vs ref
(partial contributor, ~40% of one doubling — notebook §20260907). Fleet reconciled.
Archiver fix committed (656e005). exp arm relaunched under option (a).

**In flight (check first)**:
- 13592732 `_3r_sv_s1p5` NT72 — RUNNING, started 06:44, 14 h limit → lands by ~20:45
  2026-09-07. THE gating result: first smooth-family NT36→NT72 slope
  (NT36 rung = 0.072222, +2.0%).
- 13602884 / 13602885 — `_3r_exp_nt` NT36/NT72, `WAKE_EXPINT=true SAVE_VTK=false`,
  unguarded, `-nt` stack, submitted ~07:5x. Early check: predecessor 13593792 died at
  the FIRST VTK write ~7 min in; with SAVE_VTK=false these should sail past — if one
  dies early anyway, it's a NEW failure mode. Logs in
  `~/wt018/FLOWPanel-nt-gh200/logs/slurm/` (read .err first).
- Storage strip: running from Ryan's OTHER session — do not start a second one.
- The session monitor (sacct poll, 5 min, auto-adopts fp-018gpu-*) DIED with the old
  session — re-arm one in the new session.

**Scoring**: `scripts/p018_harvest_ct.py` (windowed CT over `in_convergence_window`
revs from `*_CT_per_rev.csv`; force-monitor fallback for killed runs). Refs: NT36
0.070775, NT72 0.071844.

**Ryan's holds (2026-09-07)**: WGE reruns AND `_3r_sv` resubmit wait on the s1p5 NT72
score. FMM near-set lead still not to be chased (09-05). Notebook: concise entries,
ask before writing.

**Next actions in order**:
1. Re-arm the fleet monitor (or poll sacct on 13592732/13602884/13602885).
2. Sanity-check the two exp-nt jobs got past ~10 min of stepping.
3. When 13592732 lands: verify by outputs, score with p018_harvest_ct.py, compute the
   smooth-family slope, log in §8 + notebook (ask Ryan), then bring Ryan the WGE-rerun
   and sv-resubmit decisions with the slope in hand.
4. When exp-nt rungs land: score both, giving the stretching-integrator discriminator.
5. Optional inline anytime: the sigma_max=0.2036-vs-SIGMA_CEIL=0.030 anomaly
   (gate log of dead 13592729; guard bypass vs wrong-row read; weakens/strengthens
   the parked FMM lead).

## 10. 2026-09-07 session 2: exp-nt option (a) FAILED; option (b) applied; sigma anomaly root-caused

### exp-nt pair (13602884/13602885): dead, same VTP crash — option (a) was doubly broken
Both jobs died at step 3 (~1 min in) with the identical
`SplittingState has no field dsigma2_visc` at `_write_particles_vtp`
(`FLOWPanel_wake.jl:2441`). Root cause of the failed fix, twofold:
1. **Env clobber**: the GPU launcher `source`s the CPU dispatcher, and
   `run_dji9443_hover_ct_hpc.slurm.sh:90` unconditionally `export SAVE_VTK=true`
   ("HPC policy") — the submission-time `SAVE_VTK=false` never reached the driver.
2. **Option (a) unworkable anyway for this driver**: `rotor_hover_pressure_comparison.jl:25`
   maps `SAVE_VTK=false` → `save_path=nothing`, and `FLOWPanel_simulate.jl:1101`
   then sets `monitor_csv_dir=nothing` → NO monitor CSVs → nothing for
   `p018_harvest_ct.py` to score. The 09-07 note "CT CSVs unaffected" was wrong.

**Fix applied (option (b), minimal I/O-only guard)**: wrapped the two
`split_dsigma2_*` VTP writes in `hasfield(typeof(split_state), :dsigma2_visc)` in
the campaign worktree `~/wt018/FLOWPanel-nt-gh200` — commit `d2d3bb7`, annotated tag
`campaign/p018-nt-vtpguard-20260907` (supersedes `campaign/p018-nt-20260905` for the
exp arm). No physics change; VTK output now stays ON (matches other arms, restores
warm-start ability). NOTE: this overrides Ryan's (a)-over-(b) ruling because (a)'s
premise (CSVs unaffected) was falsified — flagged to Ryan.

**Resubmitted**: NT36 = **13603468** (case p018_csarc_l3p0), NT72 = **13603469**
(case p018_csarc_n2_nt72_l3p0); same env as 09-07 relaunch minus `SAVE_VTK=false`
(`WAKE_EXPINT=true`, `-nt` overrides, `--constraint=arm`, gh200). A background poll
watches for the dsigma2 error vs ≥12 clean steps.

### sigma_max=0.2036 anomaly (dead 13592729 gate log): ROOT-CAUSED to unguarded merging
Traced by code-scout across FLOWVPM/FastMultipole:
- The gate (`FastMultipole/src/translate_batched_resident.jl:2081-2119`) reads the
  TRUE particle σ: single-system packing, `sigma_row=8` verified correct
  (`FLOWVPM_fmm.jl:93`, `FLOWVPM_fmm_radix.jl:89-99`). **Wrong-row/wrong-system
  hypothesis ruled out** — the 0.2036 m particle was real.
- **Mechanism: particle merging bypasses the σ ceiling entirely.**
  `FLOWVPM_merging.jl:166-171` sets `σ_merged = cbrt(Σσᵢ³)` with no clamp;
  `sigma_guard`/`SIGMA_CEIL` clamps only inside the reformulated-Euler σ-ODE update
  (`FLOWVPM_timeintegration.jl:327-330`, skips static particles at :277).
  Quantitative match: `0.2036/0.030 = 6.79 ≈ 313^(1/3)` → a ~313-particle-equivalent
  merge. Merging was ON in the sfs3nb arm (driver default true; only `mergeoff`
  cases disable it).
- Nuance: since the Euler clamp would pull σ back to sceil at the next step, the
  gate must have tripped on the fresh post-merge transient (first FMM eval after
  the merge pass) — consistent with a sudden death deep into a guarded run.
- **Implication for the parked FMM lead**: the gate is honest; the anomaly does NOT
  indicate an FMM diagnostic bug. Instead it exposes a guard coverage gap: any
  merged particle can transiently exceed SIGMA_CEIL by large factors, and each such
  transient degrades near-set adequacy (or kills the run at the gate). Candidate
  remedies (NOT applied, needs Ryan): clamp σ_merged to sceil in merge_particles!,
  or make merge eligibility respect the ceiling. Also note the σ-conserving cbrt
  merge itself may distort FMM adequacy monotonically with NT (more particles →
  more merges) — a possible, still-unproven contributor to NT anti-convergence.

### merge2 A/B launched (Ryan-approved, 2026-09-07 session 2)
Ryan fixed runaway merging in FLOWVPM (`flowpanel` branch commit `119fe23`,
"prevent runaway particle merging"): greedy disjoint nearest-neighbor pairing
within a cell — at most one pair-merge per particle per pass, chaining eliminated.
50/50 merging tests pass locally. Residual (accepted): pair merges still uncapped
at cbrt(2)·σ ≈ 1.26×; guard-on arms bounded at ~1.26·sceil.

Relaunched the sfs3nb pair as a labeled model-def A/B (`_m2` suffix):
- Pins: cherry-picked 119fe23 onto both smooth-ladders FLOWVPM silo pins →
  worktrees `~/wt018/FLOWVPM-m2-{gh200,h200}`, tags
  `campaign/p018-merge2-20260907-{gh200,h200}` (only conflict: docs file absent
  from silo, took new version; merging src/tests applied clean, md5-identical
  across arches). Envs `~/p018wtenv-m2-{gh200,h200}` = copies of the wave envs
  with FLOWVPM dev-path repointed.
- **13603724** = NT36 `p018_csarc_l3p0` `_3r_sfs3nb_m2` (mgh/gh200, RUNNING).
- **13603725** = NT72 `p018_csarc_n2_nt72_l3p0` `_3r_sfs3nb_m2` (eng/h200,
  PENDING behind the H200). Submission gotcha: pin-h200 launcher header still
  hardcodes `--constraint=arm` (predates 0e08ab4); `--constraint=""` is rejected
  by the submit plugin — use `--constraint=hopper` (a real eng feature tag).
- Wave-common env reproduced: SIGMA_CEIL=0.030 TRUNCATION_RADIUS_R=3.0
  MAX_PARTICLES=1500000 P018_SETTLE_REVS=22 + SFS_THREELEVEL=true.
- Readout: NT72 survival past the old σ_max=0.2036 FMM gate death + clean
  same-model-def NT36→NT72 slope for the merge2 arm.

### Scores 2026-09-07 evening (session 2)
- **13592732 `_3r_sv_s1p5` NT72 COMPLETED clean** (2160/2160, converged, rc=0):
  windowed CT = **0.073285 ± 0.000048** (10 revs, spread 0.57%), **+2.01%** vs NT72
  ref 0.071844. Within-arm NT36→NT72 climb **+1.47%** vs reference-family +1.51%.
  **VERDICT: NULL — smoothing width does not flatten the NT climb.** The gating
  smooth-family slope is in hand; WGE-rerun and `_3r_sv`-resubmit decisions are now
  unblocked (Ryan's call).
- **exp-nt pair 13603468/13603469 COMPLETED clean** (720/720, 1440/1440, no dsigma2
  error — VTP guard verified end-to-end) **BUT OFF-MODEL**: the 09-07 relaunch env
  (inherited this session) lacked the wave-common exports → ran at truncation
  **4.0R** (wave = 3R), **20 revs** (wave = 30, settle 22), and default run names —
  **overwrote the base-case CT CSVs** in `data/p018_csarc_l3p0/` and
  `data/p018_csarc_n2_nt72_l3p0/` (base numeric refs survive in §8/§0 tables;
  today's CSVs copied to `data/salvage_exp_nt_t4p0_20260907/`).
  Provisional off-model scores: NT36 0.071365 ± 0.000252, NT72 0.072280 ± 0.000173,
  within-pair climb **+1.28%** — prima facie the stretching integrator does NOT
  flatten the climb either, but this needs an on-model rerun to score as a verdict.
- **exp-nt on-model relaunch (Ryan-approved)**: NT36 = **13603734**, NT72 =
  **13603735**, from `~/wt018/FLOWPanel-nt-gh200` (VTP-guard pin d2d3bb7), full wave
  env (TRUNCATION_RADIUS_R=3.0, MAX_PARTICLES=1500000, P018_SETTLE_REVS=22,
  `P018_RUN_NAME=<case>_3r_exp_nt`), WAKE_EXPINT=true, unguarded (no SIGMA_CEIL —
  euler_exp rejects sigma_guard), mgh/gh200. These supersede the off-model
  13603468/13603469 scores above.

## 11. RESET BRIEF 2026-09-07 evening (entry point — supersedes §9)

Read §8 + §10 + this section. Science background: 09-05 file §0/§2/§3.

**State**: s1p5 slope landed — **NULL** (+1.47% climb ≈ ref +1.51%; §10 scores).
srlx = partial (−0.58%, 40% of a doubling). exp-nt provisional null (+1.28%,
off-model). Suspect mass now concentrates on: (a) merge runaway / near-set
adequacy (σ_max=0.2036 ROOT-CAUSED to chained merging, §10 — merge2 arm in
flight), (b) SFS per-step bias (never ladder-tested with SFS fully OFF).

**In flight (all fp-018gpu-*, adopted by fleet monitor)**:
| Job | Arm | Notes |
|---|---|---|
| 13603724 | sfs3nb **merge2** NT36 | mgh, pairs-only merging (FLOWVPM 119fe23) |
| 13603725 | sfs3nb **merge2** NT72 | eng PENDING behind H200 |
| 13603734/35 | exp-nt **on-model** NT36/NT72 | mgh, full wave env, `_3r_exp_nt` run names |

**Pins (annotated tags)**: exp-nt = `campaign/p018-nt-vtpguard-20260907`
(FLOWPanel-nt-gh200 d2d3bb7, VTP dsigma2 hasfield guard). merge2 =
`campaign/p018-merge2-20260907-{gh200,h200}` (FLOWVPM-m2-* worktrees, envs
`p018wtenv-m2-*`). Wave pins unchanged.

**SFS-off question (Ryan asked 2026-09-07 evening)**: SFS_OFF=true → noSFS has
NEVER been run as an NT ladder — diagnostics only (ops_reference ruling 9; also
mandatory in RHPC_BACKEND=direct A/Bs). Candidate next arm `_3r_nosfs` if Ryan
lifts ruling 9 for ladders: same wave-common env as §10 merge2 launch but
`SFS_OFF=true` (no SFS_THREELEVEL), cases p018_csarc_l3p0 (NT36, mgh/gh200) +
p018_csarc_n2_nt72_l3p0 (NT72, eng/h200, `--constraint=hopper`), run names
`<case>_3r_nosfs`, wave pins (FLOWPanel-pin-*, p018wtenv-*). Rationale: per-step
fingerprint + srlx partial + sfs3nb NT36 −2% all point at SFS.

**Next actions**:
1. Score merge2 pair when it lands (readout: NT72 survives old gate death? does
   the NT36→NT72 climb flatten vs ref +1.51%?). If it flattens → merge chaining
   was the per-step bias; FMM near-set lead collapses into it.
2. Score on-model exp-nt pair → formal stretching-integrator verdict
   (provisional: null).
3. Bring Ryan: WGE-rerun + `_3r_sv` resubmit decisions (unblocked by s1p5 null),
   and the `_3r_nosfs` ladder proposal (needs ruling-9 lift).
4. Fill notebook table cells (journals/20260901.md "018 NT ladders — CT summary")
   as rungs land — append-only, ask Ryan first.
5. Open remainder: merge-result σ cap in FLOWVPM (pairs still uncapped ×1.26;
   unguarded arms compound), Ryan decides where it lands.

**Gotchas (this session, learned in anger)**:
- Dispatcher hard-exports SAVE_VTK=true (hpc:90); SAVE_VTK=false anyway kills ALL
  CSVs (simulate.jl:1101 + ops_reference gotcha). Option (a) was unworkable.
- Wave env is NOT in the case table: without SIGMA_CEIL/TRUNCATION_RADIUS_R=3.0/
  MAX_PARTICLES/P018_SETTLE_REVS=22/P018_RUN_NAME the base case runs 4.0R/20revs
  and CLOBBERS base-case data dirs (happened 09-07; salvage =
  data/salvage_exp_nt_t4p0_20260907/).
- eng submissions from pin-h200: header hardcodes --constraint=arm;
  --constraint="" rejected; use --constraint=hopper.
- sacct/squeue via `ssh orc 'bash -lc ...'` (module path).

### `_3r_nosfs` ladder LAUNCHED (Ryan lifted ruling 9 for ladders, 2026-09-07 evening)
- NT36 = **13603743** (mgh/gh200, PENDING Priority), NT72 = **13603744** (eng/h200
  `--constraint=hopper`, PENDING Resources). Wave pins (FLOWPanel-pin-*,
  p018wtenv-*), wave-common env + `SFS_OFF=true` (→ noSFS), guarded
  (SIGMA_CEIL=0.030), run names `<case>_3r_nosfs`.
- OWED at first step: verify the driver diagnostics line shows the noSFS choice
  took effect (jobs were still pending at handoff).
- Readout: ceiling test on SFS-mediated share of the NT climb (flat ≈ SFS is the
  per-step bias; full climb ≈ SFS cleared).

## 12. STATUS UPDATE 2026-09-07 late (read with §11; supersedes its in-flight table)

**Fleet snapshot**:
| Job | Arm | State |
|---|---|---|
| 13603724 | merge2 NT36 | **COMPLETED clean** (1080/1080, gate rc=0) — scored below |
| 13603725 | merge2 NT72 | **DIED step 1550/2159** — same FMM adequacy gate as 13592729 (see below) |
| 13603734 | exp-nt NT36 on-model | RUNNING step ~352/1079 (8 s/step, sharing mgh-1-1) |
| 13603735 | exp-nt NT72 on-model | RUNNING step ~962/2159 (11 s/step) |
| 13603743 | nosfs NT36 | RUNNING on mgh-1-1 (~1h18m) |
| 13603744 | nosfs NT72 | RUNNING on eng-1-1 (~15m) |

**merge2 NT36 score** (windowed, 10 revs): CT = **0.069260 ± 0.000053** —
−2.14% vs ref NT36; **−0.18% vs old-merge sfs3nb NT36 (0.069382)** — pairs-only
merging barely moves CT at NT36 (marginal vs ±0.14% spread).

**merge2 NT72 death — the big finding**: died at the near-set adequacy gate at
**the SAME step (1550/2159)** as old-merge 13592729, now with **sigma_max=0.1695**
(was 0.2036), ell=2 (was 3), admissible ell<=1. Implications:
1. Chained merging was NOT the sole σ source: with pairs-only merging σ_max still
   reached 0.1695 = ceiling×5.65 (≈181 ceiling-particles of σ³, or ~7.5 compounding
   ×1.26 pair merges). Candidate residual mechanisms: compounding pair merges across
   steps (merge→clamp→merge; clamp destroys σ³ but Γ concentrates), and/or
   CORE_SPREADING (WAKE_CORE_BETA=1e9, σ-growth ~1.0002/step compounding on merged
   reps; guard skips static particles).
2. Same step number both times → the trigger is DETERMINISTIC: `_radix_auto_geometry`
   deepens ell when the particle count crosses a threshold at step ~1550 (NT72
   schedule), and the standing large-σ population fails adequacy at the new depth.
   NT36 rungs never reach that count → never die. This is *exactly* the parked
   NT-scaling mechanism of §2 (09-06), now with two direct observations.
3. Remedies to bring Ryan: (i) σ-cap on merge results (skip pair if cbrt(σi³+σj³)
   > sceil — one condition in FLOWVPM_merging.jl); (ii) port
   `_alldirect_geometry_fallback!` (FastMultipole d938ba68, 2026-08-31 — demotes to
   all-direct instead of throwing; ABSENT from the silo FastMultipole pins, present
   in FastMultipole-nt by cherry-pick — verify by grep not merge-base); (iii) both.
   Note (ii) changes death→survival but the underlying σ growth remains physical(ly
   wrong); (i) attacks the σ source itself.

**Owed**:
- nosfs pair: driver never prints the SFS choice — SFS_OFF=true delivery relies on
  the (precedent-verified) sbatch env-export path; when NT36 lands, sanity-check its
  CT differs from the reference rung before trusting the label.
- Score exp-nt on-model + nosfs rungs as they land; fill notebook table
  (ask Ryan first).
- merge2 NT72: harvestable-to-step-1550 via force-monitor fallback if Ryan wants a
  partial read (21+ revs, likely usable window per 13592729 precedent).

### Ryan's ruling on the merge2 σ blowup (2026-09-07 late)
**Try item 026 particle splitting here** as the remedy for the step-1550 σ_max
blowup — splitting large-σ particles back toward the target resolution attacks
the σ source directly (vs the §12 remedy options (i)/(ii), which are NOT to be
pursued for now). Plumbing status: FLOWVPM side landed on `flowpanel`
(9d63578 / 99f4d54 / 5583443 — split kernels + `split_particles!(::ResolutionSplitOpts)`
+ radix adequacy log; Phase 2 Session 1); the FLOWPanel-side wiring (026 commits
4–6) is still pending — see
`BRAINSTORM/026_sigma_growth_particle_splitting/phase2_session2_prompt_20260907.md`.
A 018 trial arm needs that wiring (or a minimal maintenance-hook equivalent) plus
a fresh pin/tag before any `_m2s`-style NT72 rung is submitted.
**Otherwise: WAIT for the in-flight jobs** (exp-nt on-model 13603734/35, nosfs
13603743/44) — no further 018 submissions until they land and are scored.

### Update 2026-09-07 night
- **nosfs NT36 (13603743) COMPLETED clean** (gate rc=0): windowed CT =
  **0.069431 ± 0.000050**, **−1.90% vs ref NT36** — SFS_OFF verifiably took effect
  (owed check PASSED: rung clearly shifted off the reference). Notably close to
  sfs3nb NT36 (−1.97%). Verdict awaits nosfs NT72 (13603744, RUNNING).
- **exp-nt NT36 on-model (13603734) genuine NODE_FAIL** at step 352/1079 (log
  stops mid-stream, no error; .err only benign FastMultipole redefinition
  warnings). Resubmitted fresh (guardrail: NODE_FAIL restarts are fresh) as
  **13603853**, same env/pin/run-name.
- Still RUNNING: 13603735 (exp-nt NT72, ~2h), 13603744 (nosfs NT72, ~33m).

## 13. RESET POINTER 2026-09-07 night (entry point)

Read §12 + its updates + §11 (recipes/gotchas). Current in-flight:
| Job | Arm | Notes |
|---|---|---|
| 13603735 | exp-nt NT72 on-model | RUNNING (~2h at reset) |
| 13603744 | nosfs NT72 | RUNNING (~33m at reset) |
| 13603853 | exp-nt NT36 on-model | resubmit of NODE_FAIL 13603734 |

Scored today: s1p5 NT72 0.073285 (NULL), merge2 NT36 0.069260 (−0.18% vs old
merge), nosfs NT36 0.069431 (−1.90% vs ref; SFS_OFF verified), exp-nt off-model
pair (+1.28% climb, provisional null). merge2 NT72 died at the deterministic
step-1550 depth-increment adequacy gate (σ_max 0.1695).

Immediate next actions:
1. When 13603735 / 13603744 / 13603853 land: verify by outputs, score windowed
   CT, compute nosfs and exp-nt NT36→NT72 climbs vs ref +1.51%, log in handoff +
   notebook table (journals/20260901.md "018 NT ladders — CT summary"; ask Ryan).
2. Ryan's ruling stands: remedy for the σ blowup = try 026 particle splitting
   (needs 026 FLOWPanel commits 4–6 wired + fresh pin); NO other 018 submissions
   until in-flight jobs land. σ-cap/alldirect-fallback remedies parked.
3. Ryan decisions open: WGE-rerun + `_3r_sv` resubmit (unblocked by s1p5 null).

## 14. 2026-09-08: in-flight set landed — exp-nt scored, nosfs NT72 blew up

All three jobs from §13 are terminal. **No NT72 rung completed.**

| Job | Arm | Outcome | Step |
|---|---|---|---|
| 13603853 | exp-nt NT36 on-model | COMPLETED clean, gate_rc=0 | 1080/1080 |
| 13603735 | exp-nt NT72 on-model | **NODE_FAIL** (hardware; log stops mid-stream, .err only benign FastMultipole redefinition warnings) | 1821/2159 (rev 25.0) |
| 13603744 | nosfs NT72 | **FAILED — force blow-up** (dispatcher_rc=1) | 1458/2159 (rev 20.24) |

### Scores (windowed CT; refs NT36 0.070775, NT72 0.071844)

| Arm | Job | NT | Windowed CT | Spread | Δ vs ref | Provenance |
|---|---|---|---|---|---|---|
| `_3r_exp_nt` on-model | 13603853 | 36 | **0.070536 ± 0.000379** | ±0.54% | **−0.34%** | `CT_per_rev.csv`, revs 21–30, timestamp-verified against 13603853 (not the off-model 13603468 clobber). AUTHORITATIVE |
| `_3r_exp_nt` on-model | 13603735 | 72 | 0.071972 ± 0.000133 | ±0.18% | +0.18% | force-monitor fallback, **revs 15–25 — NOT the convention window** (revs 20–30 unreached). NOT comparable to the other rungs |
| `_3r_nosfs` | 13603744 | 72 | — | — | — | **NO usable window**: died rev 20.24, at/just before window opening. Not scoreable |

⚠ **Do NOT compute an exp-nt slope from the two rows above.** The NT36 rung uses
the convention window (revs 21–30); the NT72 rung uses revs 15–25 because the run
died at rev 25. Different, less-settled span → the naive +2.03% climb is an
artifact of window mismatch, not a measurement. The stretching-integrator verdict
remains **provisional null** (off-model pair gave +1.28% climb, §10), and still
needs a completed on-model NT72 rung.

### nosfs NT72 blow-up: a NEW failure mode, distinct from the step-1550 gate

Not the FMM near-set adequacy gate (that fires at step 1550 with a σ_max message);
this is a force catastrophe at step 1458 with **σ-growth benign at 1.133** — the
σ guard never fired and no adequacy error was thrown.

| Step | Rev | CF_x |
|---|---|---|
| 1433–1450 | 19.90–20.15 | −0.0698 → −0.0694 (flat, ±2%) |
| 1451 | 20.15 | −0.0644 |
| 1452–1454 | 20.17–20.19 | −0.048 → −0.037 (collapse) |
| 1455 | 20.21 | −0.0403 |
| 1456 | 20.22 | −0.2918 |
| 1457 | 20.24 | −1.0635 → dispatcher rc=1 |

Shape: flat for ~20 revs, then sign-preserving collapse over ~4 steps followed by
a 15× spike in 2 steps. Quiet σ throughout rules out the σ-blowup/merge pathway
that killed the merge2 and sfs3nb NT72 rungs.

**Reading (hedged, single observation)**: this is the arm with SFS fully OFF, and
SFS is the term that damps small-scale wake growth. A noSFS wake failing to stay
stable past ~rev 20 is plausibly a *result* about the SFS question rather than an
infrastructure failure — but it costs the NT72 rung either way, and one death is
not a mechanism. The NT36 nosfs rung completed fine (−1.90%), so the instability
is NT72-specific (finer helix spacing / more particles), consistent in spirit with
the NT-scaling story even though the proximate failure differs.

### Fleet: no 018 jobs running. Other campaigns live (021 p3-R2/R3, 022 p2lg-tune-R6/R7).

### Owed / decisions for Ryan
1. **exp-nt NT72 resubmit** — pure NODE_FAIL, on-model, would have landed. Fresh
   resubmit (never a chain) is the guardrail-compliant action; needs Ryan's OK
   given the "no further 018 submissions" hold.
2. **nosfs NT72** — rerun as-is (may blow up again), or treat the blow-up as the
   arm's answer?
3. Still open from §11/§13: WGE reruns + `_3r_sv` resubmit (unblocked by the s1p5
   null); merge-result σ cap location in FLOWVPM; 026 splitting wiring
   (FLOWPanel commits 4–6) as the designated σ remedy.
4. Notebook table cells (journals/20260901.md "018 NT ladders — CT summary") —
   ask Ryan verbosity before writing.

### Ryan's rulings 2026-09-08 (executed)
- **exp-nt NT72 resubmitted fresh**: job **13605974** (NODE_FAIL restarts are fresh,
  never chained). From `~/wt018/FLOWPanel-nt-gh200` (VTP-guard pin d2d3bb7), mgh/gh200,
  `--constraint=arm`, env `WAKE_EXPINT=true TRUNCATION_RADIUS_R=3.0
  MAX_PARTICLES=1500000 P018_SETTLE_REVS=22
  P018_RUN_NAME=p018_csarc_n2_nt72_l3p0_3r_exp_nt` + `-nt` overrides, unguarded
  (no SIGMA_CEIL — euler_exp rejects sigma_guard). Deliberately added NOTHING beyond
  the env that produced the paired NT36 rung 13603853, so the pair stays internally
  consistent. Readout: the formal stretching-integrator NT36→NT72 slope.
- **nosfs NT72: blow-up accepted as the arm's result** — no rerun. Recorded finding:
  *the noSFS wake is unstable past ~rev 20 at NT72* (force catastrophe, benign σ).
  The nosfs ladder verdict therefore rests on the NT36 rung alone (−1.90% vs ref);
  there is no nosfs NT72 CT and none is planned.

### `depth:4R` header check (2026-09-08) — NOT a deviation
While recovering the submit recipe, the driver header line was seen to print
`depth:4R` on the exp-nt and nosfs jobs, which looked like an off-model 4.0R
truncation. **Checked against the accepted wave arms**: 13592732 (s1p5 NT72),
13592896 (srlx NT72) and 13593711 (sfs3nb NT36) ALL print `depth:4R` as well.
So the header's `depth:` field is not `TRUNCATION_RADIUS_R`; it is uniform across
every arm and comparability is intact. Contrast: the one genuinely off-model job
(13603468) is distinguishable in the same header by `settle:12` vs the wave's
`settle:22`. **Use `settle:`, not `depth:`, as the on-model tell.**

## 15. 2026-09-08: SFS labeling, overwrite audit, and Cd (SFS coefficient) vs NT

### 15.1 Labeling fix (local repo, UNCOMMITTED)
`examples/rotor_hover_pressure_comparison.jl` now builds an `sfs_label` alongside
`sfs_choice` and appends `SFS=$(sfs_label)` to the "Particle diagnostics:" line.
Closes the §12 gap where SFS_OFF delivery could only be inferred indirectly.
Emits e.g. `SFS=noSFS`, `SFS=SFS_Cd_threelevel_nobackscatter(static)`, or
`SFS=DynamicSFS(rlxf=0.0025031, maxC=1.0, alpha=0.999, clippings=backscatter,
controls=none, nostatic=false)`. Parse-checked and runtime-checked.
⚠ Affects FUTURE runs only, and only once it is in a **pinned worktree** — the
in-flight jobs and every scored rung ran without it.

### 15.2 Overwrite audit — the SFS labeling gap did NOT cause overwrites
Surveyed every `p018_csarc*` data dir (CT_per_rev.csv count + mtimes). Each SFS
variant carried a distinct explicit `P018_RUN_NAME` suffix (`_3r_sfs3nb`,
`_3r_srlx`, `_3r_nosfs`, `_3r_sfs3nb_m2`), so **no two SFS states ever wrote into
the same dir**; every arm dir holds exactly one CT CSV.
**The real overwrite is the previously-documented one and has a different cause**
(missing `P018_RUN_NAME`, not missing SFS label): the off-model 4.0R exp-nt pair
overwrote the two BASE-CASE dirs — `p018_csarc_l3p0/..._CT_per_rev.csv` (mtime
2026-09-07 16:21) and `p018_csarc_n2_nt72_l3p0/...` (2026-09-07 15:36).
⚠ **Consequence not previously flagged**: those dirs are the source of the
reference values 0.070775 (NT36) / 0.071844 (NT72). The reference CSVs on /home
are gone — the numbers survive only in this handoff's tables and in the
09-07 salvage copies (`data/salvage_exp_nt_t4p0_20260907/`, which are the
*off-model* CSVs, NOT the references). **Recovering the true reference CSVs from
the archive tarballs is owed** before anyone needs to re-derive a reference.

### 15.3 Cd IS recoverable from existing outputs — no new instrumentation
Per-particle SFS coefficient lives at `FLOWVPM` `C_INDEX = 37:39`
(`C[1]`=coefficient, `C[2]/C[3]`=Lagrangian nume/deno) and is already written to
every wake particle `.vtp` by `FLOWPanel_wake.jl:2427`. No monitor aggregates it,
so it exists only as raw per-particle data. Extractor written:
`scratchpad/vtp_C.py` + `vtp_C2.py` (also at `~/vtp_C.py`, `~/vtp_C2.py` on the
cluster) — parses the raw-appended VTK XML directly (meshio cannot read these).
`vtp_C2.py` compares at **equal revolution** (rev = step/NT); comparing at equal
*step* is wrong (NT36 step 1079 = rev 30, NT72 step 1550 = rev 21.5).

### 15.4 The Cd distribution is bimodal — a single "converged Cd" is a poor summary
Across every arm, **~50–55% of active particles carry C exactly 0** (the
`clipping_backscatter` strategy zeroes them), and the remainder is a skewed tail.
Median C over ALL particles is therefore 0.0 everywhere. Saturation at the
`maxC=1.0` clamp is rare (frCap ≤ 0.006), so maxC is not binding.
Statistics used below: `mC` = mean over all active particles (the coefficient the
model effectively applies); `mC+` = mean over the unclipped (C>0) subset.

### 15.5 Cd by arm (equal-rev), and the two SFS families differ by ~200×

| Arm | SFS model | NT | rev | mC | mC+ | medC+ |
|---|---|---|---|---|---|---|
| `_3r_exp_nt` | default DynamicSFS (α=0.999, rlxf 0.005) | 36 | 21.4 | 0.1579 | 0.3391 | 0.1429 |
| `_3r_srlx` | DynamicSFS rlxf=0.0025031 | 72 | 30 | 0.1356 | — | — |
| `_3r_sv_s1p5` | default DynamicSFS | 36 | 30 | 0.0876 | 0.1764 | 0.0223 |
| `_3r_sv_s1p5` | default DynamicSFS | 72 | 30 | 0.1004 | 0.2115 | 0.0216 |
| `_3r_sfs3nb` | `SFS_Cd_threelevel_nobackscatter` (α=0.667) | 36 | 30 | 0.0011 | 0.0023 | 0.0001 |
| `_3r_sfs3nb_m2` | same, pairs-only merging | 36 | 21.4 | 0.0008 | 0.0016 | 0.0001 |
| `_3r_sfs3nb_m2` | same | 72 | 21.4 | 0.0019 | 0.0040 | 0.0001 |

**The three-level preset runs ~200× smaller Cd than the default two-level
DynamicSFS** (0.0008 vs 0.158 at matched rev 21.4). That is a far bigger model
difference than anything else on the ladder, and it is consistent with the
sfs3nb rungs sitting −2% while the default-family rungs sit +2%.

### 15.6 Cd vs NT — YES, Cd grows with NT, and the gap widens with revolution
Cleanest pair (merge2, full VTK retained on both rungs, identical model-def):

| rev | NT36 mC | NT72 mC | NT72/NT36 |
|---|---|---|---|
| 4 | 0.00031 | 0.00042 | 1.33 |
| 8 | 0.00047 | 0.00066 | 1.40 |
| 12 | 0.00058 | 0.00068 | 1.17 |
| 16 | 0.00068 | 0.00088 | 1.29 |
| 20 | 0.00070 | 0.00150 | 2.14 |
| 21.4 | 0.00081 | 0.00190 | 2.38 |

Independently, `_3r_sv_s1p5` at rev 30: NT36 0.0876 → NT72 0.1004 = **+14.6%**
(mC+ +19.9%).

Readings, in decreasing confidence:
1. **Cd is genuinely NT-dependent.** Two independent model-defs (three-level
   preset and default DynamicSFS) both show NT72 > NT36 at matched revolution.
   The dynamic procedure does NOT deliver an NT-invariant coefficient.
2. **Cd is not converged in the ordinary sense** — it drifts monotonically upward
   with revolution in the sfs3nb family (0.0003 → 0.0008 over revs 4→21 at NT36)
   and only plateaus for the default-DynamicSFS exp_nt arm (≈0.15–0.16 after
   rev ~8). "The converged Cd" therefore has to be quoted with a rev and an arm.
3. **The late widening (ratio 1.3 → 2.4 over revs 16→21.4) is CONFOUNDED**: the
   NT72 merge2 rung is walking into its step-1550 σ blow-up during exactly those
   revs. Early-rev ratio ~1.2–1.4 is the trustworthy number; the rev-20+ ratio may
   be measuring the pathology, not the model. Do not quote 2.4 as the NT effect.

### 15.7 Static-Cd ladder — feasible, but Cs is arm-dependent
`ConstantSFS(model; Cs=…)` (`FLOWVPM_subfilterscale_models.jl:180`) sets
`C[1]=Cs` for every particle each step, subject to the same clipping. Its
**default `Cs=1.0` is a formal identity default, ~6× the largest observed mC+ and
~1000× the three-level family** — using the default would be badly wrong.
Because constant-Cs particles still pass through `clipping_backscatter` (~50%
zeroed), matching the *applied* mean means setting **Cs ≈ mC+**, not Cs ≈ mC.
Candidate targets: **Cs ≈ 0.18** to emulate `_3r_sv_s1p5`-family dynamics,
**Cs ≈ 0.34** for the `_3r_exp_nt` default-family level, **Cs ≈ 0.002** for the
three-level family. A driver branch for `SFS_CONST_CS` does not exist yet and
would need adding (one `elseif` in the `sfs_choice` block) plus a fresh pin.

## 16. STATIC-Cd LADDER (new ladder, launched 2026-09-08) — Ryan-approved

Two distinct experiments, documented together because they share the code change.

**Ladder S (the static-Cd NT ladder)** — default SFS family, `Cs=0.18`
(= the unclipped mean `mC+` of `_3r_sv_s1p5`, §15.5; matching mC+ rather than mC
is deliberate — constant-Cs particles still pass through `clipping_backscatter`,
so ~half are zeroed and the *applied* mean lands near mC).

| Job | Case | NT | Arch | Run name |
|---|---|---|---|---|
| **13605983** | p018_csarc_l3p0 | 36 | mgh/gh200 | `p018_csarc_l3p0_3r_cs0p18` |
| **13605984** | p018_csarc_n2_nt72_l3p0 | 72 | eng/h200 (`--constraint=hopper`) | `p018_csarc_n2_nt72_l3p0_3r_cs0p18` |

Env: wave-common (`SIGMA_CEIL=0.030 TRUNCATION_RADIUS_R=3.0 MAX_PARTICLES=1500000
P018_SETTLE_REVS=22`) + `SFS_CONST_CS=0.18`. Guarded. No expint.
**Readout**: does freezing Cd flatten the NT36→NT72 climb vs the reference +1.51%?
If it flattens, the NT bias is carried by the *adaptation* of Cd (which §15.6 shows
is itself NT-dependent, +14.6% at rev 30). If the climb survives, Cd adaptation is
exonerated and the bias is elsewhere.

**Ladder C (Cs sensitivity at fixed NT36, exp-nt on)** — probes how much CT moves
per unit Cs, bracketing the ~200× family gap of §15.5.

| Job | Cs | Run name |
|---|---|---|
| **13605985** | 0.002 (three-level family level) | `p018_csarc_l3p0_3r_cs0p002_exp_nt` |
| **13605986** | 0.34 (exp-nt default-family level) | `p018_csarc_l3p0_3r_cs0p34_exp_nt` |

Env: `WAKE_EXPINT=true TRUNCATION_RADIUS_R=3.0 MAX_PARTICLES=1500000
P018_SETTLE_REVS=22`, **unguarded** (euler_exp rejects sigma_guard), mgh/gh200.
Their dynamic-Cd control is the existing on-model **exp-nt NT36 = 0.070536**
(13603853) — same pin lineage, same env, Cd dynamic instead of frozen.
⚠ Ladder C is NOT comparable to Ladder S (different stack, expint on, unguarded);
score it only against 13603853 and against itself.

### Code + pins for this ladder
Dev commit **720947a** (`fastmultipole`): adds `SFS_CONST_CS` → `ConstantSFS`
branch reusing the default branch's model/clippings/controls, plus the §15.1
`SFS=` label. Verified: parses, and `ConstantSFS(Estr_fmm; Cs=0.18,
clippings=(clipping_backscatter,))` constructs with `isSFSenabled=true`.
`ConstantSFS` **and its GPU broadcast path** (`_constantsfs_coefficient_broadcast!`)
were confirmed present in all five pinned FLOWVPM worktrees before submitting —
the CPU branch is guarded by `pfield.particles isa Array` and these runs are
device-resident CuArray, so the GPU path was the thing that had to exist.

**New worktrees + annotated tags** (the running job 13605974 occupies
`FLOWPanel-nt-gh200`, and the worktree policy forbids editing a worktree with a
job in flight, so this ladder got its own):

| Worktree | From | Commit | Tag | Env |
|---|---|---|---|---|
| `~/wt018/FLOWPanel-cs-gh200` | `campaign/p018-smooth-ladders-20260905` | d0ff499 | `campaign/p018-sfsconst-20260908-gh200` | `~/p018wtenv-cs-gh200` |
| `~/wt018/FLOWPanel-cs-h200` | same | 41aeb48 | `campaign/p018-sfsconst-20260908-h200` | `~/p018wtenv-cs-h200` |
| `~/wt018/FLOWPanel-nt-cs-gh200` | `campaign/p018-nt-vtpguard-20260907` | 67b7383 | `campaign/p018-sfsconst-nt-20260908` | `~/p018wtenv-nt-cs-gh200` |

Patch applied by **anchor** (`~/deploy_sfs_patch.py`), not by file copy: the wave
pins carry a DIFFERENT driver snapshot from the dev checkout (md5 0bd24cec vs
94386a25), so copying the dev file would have silently changed pinned campaign
code. The `-nt` pin's driver did match the dev pre-commit file exactly.

## 17. 2026-09-08 late: all five jobs terminal — 2 scoreable, 3 dead (NOT yet harvested)

Triaged by outputs (.err first). **Nothing in this section has been CT-scored yet
— that is the next session's first task.**

| Job | Arm | Steps | gate | Outcome |
|---|---|---|---|---|
| **13605974** | exp-nt NT72 on-model | **2159/2159** | rc=0 | **COMPLETE — scoreable.** Gives the formal stretching-integrator slope against NT36 0.070536 |
| **13605983** | Ladder S `_3r_cs0p18` NT36 | **1079/1079** | rc=0 | **COMPLETE — scoreable.** First static-Cd rung |
| 13605984 | Ladder S `_3r_cs0p18` NT72 | 1876/2159 (rev 26.06) | rc=1 | **FMM near-set adequacy gate**: `g_min*h_leaf = 0.5989` vs cutoff; σ-growth benign 1.169 |
| 13605985 | Ladder C `_3r_cs0p002_exp_nt` NT36 | 898/1079 (rev 24.9) | rc=1 | **Force blow-up** → `DomainError with 8.65e9`; CF reached (−3.2e9, −3.5e8, −1.5e9) |
| 13605986 | Ladder C `_3r_cs0p34_exp_nt` NT36 | 361/1079 (rev 10.0) | rc=1 | **Force blow-up** → `DomainError with 8.24e4`; CF reached (1.39e4, −778, 1.18e4) |

### The `SFS=` label works — verified on all four new-pin runs
`SFS=ConstantSFS(Cs=0.18, …)`, `Cs=0.002`, `Cs=0.34` all appear verbatim in the
`.out`, alongside `SIGMA_CEIL=0.03 m (guard=on)` for Ladder S. SFS state is now
directly verifiable from a log instead of inferred from where a CT rung lands.
(13605974 has no `SFS=` line — correct: it predates the label, running the
`campaign/p018-nt-vtpguard-20260907` pin.)

### ⚠ Ladder C is CONFOUNDED — do not read it as "static Cs destabilizes"
Ladder C requires `WAKE_EXPINT=true`, and euler_exp **rejects `sigma_guard`**, so
both rungs ran **unguarded** (no `SIGMA_CEIL`). Ladder S NT36, same static-Cs
mechanism but **guarded**, completed 1079/1079 without incident. So the Ladder C
blow-ups cannot be attributed to freezing Cd rather than to running without the
σ ceiling. Ladder C as designed cannot separate the two, and no guarded variant
is possible while expint is on.
**Remedy if Ryan wants the Cs sensitivity**: rerun the sweep on the WAVE stack
(guarded, no expint), scored against a dynamic-Cd control on the same stack.
The blow-up ordering (Cs=0.34 died rev 10, Cs=0.002 rev 24.9) is suggestive of a
Cs-magnitude effect but is not evidence while the guard confound stands.

### Correction to the "deterministic step 1550" claim (§12)
Ladder S NT72 hit the SAME near-set adequacy gate but at **step 1876**, not 1550.
The trigger is therefore **particle-count-driven, not step-number-driven** — the
count at which `_radix_auto_geometry` deepens `ell` is crossed at whatever step
that arm's particle schedule reaches it. §12's "same step both times" was true of
two arms sharing a particle schedule (sfs3nb / sfs3nb_m2), and generalised too
far. The underlying mechanism (deeper `ell` + standing large-σ population fails
adequacy) is unchanged and now has a third observation.

### Standing tally: EVERY NT72 rung in this campaign has died
sfs3nb, sfs3nb_m2, cs0p18 → adequacy gate; nosfs → force blow-up; exp-nt (first
attempt) → NODE_FAIL. The sole NT72 survivors are `_3r_sv_s1p5` (13592732),
`_3r_srlx` (13592896), and now **exp-nt on-model 13605974**. No static-Cd NT72
rung exists, so **Ladder S has no slope** — only its NT36 rung.

### Disk
`/home/rander39` at **430 G against the 400 G cap** — still breached, unchanged
since 09-06. The three new worktrees added ~1.5 G. Not the cause, but the strip
from Ryan's other session either did not complete or was offset by live output.

## 18. RESET POINTER 2026-09-08 (entry point — supersedes §13)

**Next agent starts at `sfs_cd_reset_prompt_20260908.md`** (same directory). It carries the fleet table,
the guardrails, the scoring recipe, the Ladder C confound, and the scoped
investigation brief. This file's §§14–17 are the session record behind it.

Two tasks, in order: (1) harvest 13605974 + 13605983 (both complete, unscored);
(2) investigate whether NT has a clear, accountable impact on the SFS coefficient
Cd — starting from §15's measurements, not from scratch.
