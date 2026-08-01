# Phase 2e Log — Unsteady Hover CT Convergence

Plan: [Phase 2e — Unsteady Hover CT Convergence](../../plans/dji_convergence_20260722/phase_02e_unsteady_ct.md)

Dashboard: [DJI convergence progress](../../plans/dji_convergence_20260722/dji_convergence_20260722.md)

## Current snapshot

Status (2026-07-27): **BASELINES COMPLETE — first real Phase 2e numbers.
Batch 2 all OOM'd and has been resubmitted with more memory.**

**Headline result — unfiltered baselines, 5400 RPM, production mesh, NT=36,
truncation 4R, both `COMPLETED` (~10 h each, 719 steps, 20 revs):**

| formulation | 10-rev cycle-mean CT | std | 5-rev window |
|---|---:|---:|---:|
| velocity | **0.06148** | ±0.00152 (2.5%) | 0.06052 ± 0.00134 |
| green | **0.06205** | ±0.00188 (3.0%) | 0.06072 ± 0.00155 |

- **Green − velocity = +0.93%** on the 10-rev mean (+0.34% on the last 5). That
  is **smaller than the ±2.5–3% cycle scatter**, so Green reconstruction is
  *not resolved* as an improvement by this pair — and it is far short of the
  large gap-closing effect the study hoped for. The +1.5% seen mid-run at revs
  8–15 shrank as the wake developed; it was transient, not the settled answer.
- **Neither converged.** Per-rev spread 2.7–2.8% (tol 0.5%), p-p 4.4–5.7%
  (tol 2%), and revs 10→19 drift *downward* 0.0628 → 0.0596 (velocity). The
  wake is still settling at rev 19; this is not a stationary limit cycle yet.
- vs references: steady Dirichlet rigid/semi-infinite 0.0505–0.0515; prior
  6000-RPM unsteady 0.0617–0.0622 (our 5400-RPM value lands on it, as expected
  for a near-RPM-independent CT); **experiment 0.072 — still ~15% above us.**

- Batch 1: 12894164 (`p2e_vel_nt36_d4`) and 12894165 (`p2e_green_nt36_d4`) are
  the **filter-off baselines**. The NT=72 pair and the depth-6R pair were
  cancelled 08:15 (see the dated entry).
- Batch 2 (submitted 08:55): 12894470/12894471/12894472/12894473 — a **2x2
  factorial in (particle sigma, Das eta)** at relaxation-filter 2.0R, green.
- **TRUNCATION DEPTH: ANSWERED, null** (|dCT| <= 0.22%, 4R vs 6R, our own data).
- **Judge filtered cases on the 10-rev CYCLE-MEAN**, not the 5-rev criterion --
  the filtered wake settles into a ~9-rev limit cycle, so the 5-rev test cannot
  pass by construction and is retained only as a diagnostic.
- **6 active study jobs = the cap** (raised from 3 by Ryan on 2026-07-25).
  Do not submit more until one finishes.
- Blow-up fallback for the sigma-fine cases: resubmit with `OVERLAP=4`.
- No ParaView history retained or deleted by this phase yet.
- **Cluster access is currently BLOCKED (2026-07-27).** `ssh orc` only works
  while an SSH **ControlMaster socket** exists
  (`~/.ssh/config`: `ControlPersist yes`,
  `ControlPath ~/.ssh/master-%r@%h:%p.socket`). Auth is *not* passwordless — BYU
  requires keyboard-interactive 2FA. The socket has since expired
  (`ls ~/.ssh/master-*.socket` → none), so every connection now returns
  `Permission denied (keyboard-interactive)` and no agent can reach the cluster.
  **Fix: Ryan runs `! ssh orc -fN` once** to re-establish the master; agent
  connections then multiplex over it. Wrap remote `squeue`/`sbatch` in
  `bash -lc "..."`; parse only via `scripts/p2e_status.sh` sentinels.
- **Memory request raised 24 GB → 64 GB** (see the 2026-07-25 entry): the
  production mesh is 36,752 triangular panels, not the ~16.6k the plan assumed,
  and the Green case peaks near 30 GB.

Handoff facts carried in from earlier phases:

- **Mesh (Phase 2d production recipe):**
  `examples/data/dji9443_20260725_45_185_capped_captess4.msh` — flat root cap +
  round `CapUMinTess=4` tip cap, 45 spanwise, n_airfoil=185. Never use the
  20260723 round-ct3 `capped` meshes for Dirichlet solves (tip-circulation
  artifact).
- **Steady 5400-RPM Dirichlet CT baselines:** ≈0.0505–0.0515 (semi-infinite
  rigid wake).
- **Experimental reference:** `CT_exp = 0.072` at 5400 RPM.
- **Only converged unsteady particle-wake point to date:** 6000 RPM, NT=36,
  10 revs, `CT_bernoulli ≈ 0.0617–0.0622` (σ≈1.5e-4), from
  `data/rotor_hover_pressure_comparison/rotor_hover_pressure_comparison_CT_vs_rev.csv`.
- Phase 2e front-runs Phase 4 (Green gate) and Phase 5 (CT convergence)
  elements at Ryan's direction; Green results are always reported next to the
  matched velocity run.

## Working records

| Date | Command/file | Purpose/result |
|---|---|---|
| 2026-07-25 | `plans/dji_convergence_20260722/phase_02e_unsteady_ct.md` | Phase 2e plan created (config, run matrix, convergence criterion, deliverables) |
| 2026-07-25 | `logs/dji_convergence_20260722/phase_02e_unsteady_ct.md` | This log created |
| 2026-07-25 | `examples/rotor_hover_pressure_comparison.jl` | Driver edits: `RHPC_MESH` honored (hard-coded 40_40 override removed) + `45_185_ct4`/`45_145_ct4` keys + `RHPC_MESH_FILE` escape hatch; auto TE seeds for non-legacy meshes; `RHPC_FORMULATION` ∈ {velocity, green}; per-rev block stats + convergence verdict; `<run>_CT_per_rev.csv` and `<run>_case_metadata.toml` |
| 2026-07-25 | `examples/run_dji9443_hover_ct_hpc.slurm.sh` | New launcher, one case per job (case tag as `$1`), 64 threads / 64 GB / 24 h |
| 2026-07-25 | `scripts/deploy_phase02e.sh` | Deployment + md5 verification helper; default remote corrected to the working `orc` alias and the remote call wrapped in `bash -lc` |
| 2026-07-25 | `julia --project -t 6 test/formulation_test.jl` | **PASS** — all 10 stages, incl. Stage 5/6 Green machinery and Green-vs-TraceCorrected consistency (`‖γ_trace−γ_green‖/‖γ‖ = 3.15e-13`) |
| 2026-07-25 | driver smoke, 40_40, `RHPC_FORMULATION=velocity`, NT=6, 6 revs | exit 0; finite CT; no unpaired-shedding warning; per-rev CSV + metadata written; window mean CT 0.129096 |
| 2026-07-25 | driver smoke, 40_40, `RHPC_FORMULATION=green`, NT=6, 6 revs | exit 0; Green validation passed on the particle wake under FMM; finite CT; window mean CT 0.130069 (+0.75% vs velocity, coarse smoke settings only — not a physics claim) |
| 2026-07-25 | ct4 body-construction preflight (scratchpad `ct4_preflight.jl`) | Auto TE seeds `[17743,17559,9191]` / `[8553,8369,1]`; **ncells 36752**; body is `RigidWakeBody{_,2,_,true}`, `semiinfinite_wake=false`; **0 unpaired shedding edges** on both blades (40 edges each); shedding roots at r/R = ∓0.111 (requested ≥0.1) |

## Active job ledger

Batch 1 submitted 2026-07-25 02:39:36 (3 jobs = the study cap; queue was empty).

| Job ID | Case | Key settings | Restart | Submitted | State | Logs |
|---|---|---|---|---|---|---|
| 12894164 | `p2e_vel_nt36_d4` | velocity, NT=36, depth 4R, rlxf 0.3, 64 thr / 64 GB / 24 h, node m12-1-29 | cold | 2026-07-25 02:39 | RUNNING | `logs/slurm/slurm-fp-p2e-vel-nt36-d4-12894164.{out,err}` |
| 12894165 | `p2e_green_nt36_d4` | green (`:area_mean`, interval 1, dense), NT=36, depth 4R, rlxf 0.3, 64 thr / 64 GB / 24 h, node m12-2-1 | cold | 2026-07-25 02:39 | RUNNING | `logs/slurm/slurm-fp-p2e-green-nt36-d4-12894165.{out,err}` |
| 12894166 | `p2e_vel_nt72_d4` | velocity, NT=72, depth 4R, rlxf 0.15, 64 thr / 64 GB / **48 h**, node m12-2-2 | cold | 2026-07-25 02:39 | **CANCELLED 08:15** (NT=72 deferred; reached rev 8.4, dir KEPT for warm start) | `logs/slurm/slurm-fp-p2e-vel-nt72-d4-12894166.{out,err}` |
| 12894168 | `p2e_green_nt72_d4` | green, NT=72, depth 4R, rlxf 0.15, 64 thr / 64 GB / **48 h**, node m12-2-9 | cold | 2026-07-25 02:48 | **CANCELLED 08:15** (NT=72 deferred; reached rev 8.3, dir KEPT for warm start) | `logs/slurm/slurm-fp-p2e-green-nt72-d4-12894168.{out,err}` |
| 12894169 | `p2e_vel_nt36_d6` | velocity, NT=36, **depth 6R**, rlxf 0.3, 64 thr / 64 GB / **36 h**, node m12-2-14 | cold | 2026-07-25 02:48 | **CANCELLED 08:15** (depth question dropped; reached rev 12.9, dir DELETED) | `logs/slurm/slurm-fp-p2e-vel-nt36-d6-12894169.{out,err}` |
| 12894170 | `p2e_green_nt36_d6` | green, NT=36, **depth 6R**, rlxf 0.3, 64 thr / 64 GB / **36 h**, node m12-2-14 | cold | 2026-07-25 02:48 | **CANCELLED 08:15** (depth question dropped; reached rev 12.9) | `logs/slurm/slurm-fp-p2e-green-nt36-d6-12894170.{out,err}` |

**Batch 2 (2×2 factorial, near-wake relaxation filter), submitted 2026-07-25 08:55.**
All: green, filter 2.0R, NT=36, truncation 4R, `SETTLE_REVS=20`, `CONVERGENCE_REVS=10`.

| Job ID | Case | σ config | `Das` η | Wall | State |
|---|---|---|---|---|---|
| 12894470 | `p2e_green_f2p0_sigC_das0p2` | stock (ov 3, pps 2) | 0.2 | 24 h | RUNNING |
| 12894471 | `p2e_green_f2p0_sigC_das4p0` | stock | 4.0 | 24 h | RUNNING |
| 12894472 | `p2e_green_f2p0_sigF_das0p2` | fine (ov 2, pps 8, merge 0.0025R) | 0.2 | 48 h / 64 GB | **OUT_OF_MEMORY** at step 56 (14 min) → resubmitted as 12894481 |
| 12894473 | `p2e_green_f2p0_sigF_das4p0` | fine | 4.0 | 48 h / 64 GB | **CANCELLED** pre-emptively (same config, same fate) → resubmitted as 12894482 |
| **12894481** | `p2e_green_f2p0_sigF_das0p2` | fine, σ unchanged | 0.2 | 48 h / **256 GB** | RUNNING (resubmit) |
| **12894482** | `p2e_green_f2p0_sigF_das4p0` | fine, σ unchanged | 4.0 | 48 h / **256 GB** | RUNNING (resubmit) |

## Batch 3 — staged studies A and B (submitted 2026-07-27)

Carrier: unfiltered production mesh, NT=36, truncation 4R, `SETTLE_REVS=12`
(719 steps, matching the completed baselines), **velocity** formulation,
128 GB / 24 h. Staged in `BRAINSTORM/006` "Staged studies (2026-07-27)".

| Job ID | Case | Study | Filter | `Das` η | State |
|---|---|---|---|---|---|
| 12918580 | `p2e_filt0p5` | B — filter survivability | 0.5R | 0.2 | RUNNING |
| 12918581 | `p2e_filt1p0` | B — filter survivability | 1.0R | 0.2 | RUNNING |
| 12918582 | `p2e_das1p0` | A — `Das` convergence | off | 1.0 | RUNNING |
| 12918583 | `p2e_das4p0` | A — `Das` convergence | off | 4.0 | RUNNING |

**`p2e_das0p2` deliberately NOT submitted** — it is an exact re-run of the
completed velocity baseline 12894164 (unfiltered, NT=36, 4R, η=0.2, 719 steps;
the driver's hard-coded `init_Das_eta_kinematic` was 0.2, identical to the env
default). That ladder point already exists at **CT = 0.06148 ± 0.00152**
(10-rev cycle-mean). Saves a slot with no loss of information.

**Held for the next free slots:** `p2e_das0p5`, `p2e_das2p0`. The ladder was
started at η ∈ {0.2 (baseline), 1.0, 4.0} — a 20× span across three points. If
1.0 and 4.0 agree, `Das` is converged by η=1.0 and the intermediate rungs are
unnecessary; if they disagree, η=2.0 becomes the informative fill and η=0.5
brackets the low end.

## Completed case ledger

No study cases yet.

| Case | Formulation | NT | Depth | CT mean (final 5 revs) | Per-rev spread / p-p | Converged? | Wall/RSS | Artifacts |
|---|---|---:|---:|---:|---|---|---|---|

## ParaView retention ledger

Full histories currently on disk: **8** (2 running baselines + 2 archived NT=72
+ 4 Batch-2) — above the ≤4 cap; trim at analysis time, keeping the CT CSVs.
Cluster disk 192 G / 2.0 T, so there is no pruning pressure.

| Action/date | Exact path | Summary retained | Warm-start check | Reason |
|---|---|---|---|---|
| deleted 2026-07-25 | `data/p2e_vel_nt36_d6`, `data/p2e_green_nt36_d6` (9.2 GB) | `data/p2e_harvest/{vel,green}-nt36-d6_per_rev.csv` | n/a | truncation-depth question permanently dropped; result was a null |
| kept 2026-07-25 | `data/p2e_vel_nt72_d4`, `data/p2e_green_nt72_d4` (12.3 GB) | `data/p2e_harvest/*nt72*_per_rev.csv` | intended restart point ~rev 8.3 | Ryan: "may come back to NT=72 later"; keeping the history avoids re-paying the cold wake fill. Cluster is 190 G / 2.0 T so this costs nothing there; delete if space is ever needed |
| destroyed 2026-07-27 | `data/p2e_green_f2p0_{sigC_das0p2,sigC_das4p0,sigF_das0p2}` | `logs/slurm/*.out` CT series (extracted) | n/a | wiped unintentionally by the retry launch — see the dated entry |
| rescued → local 2026-07-27 | `data/p2e_OOM_sigF_das4p0_12894487` (1.3 GB) | full VTK, 129 steps, 64 k particles | n/a | copied to the local repo `data/` for ParaView; **cluster copy deleted** after transfer. Delete the local copy once viewed — it lives in Dropbox |

## Dated entries

Immediately record exact submission commands, settings, job IDs/states, elapsed
time, peak memory, restart provenance, logs, numerical results, failures and
root causes, md5 verification of deployed assets, and every retained or removed
history. Keep the snapshot and ledgers current after each batch.

### 2026-07-29 — 014 η×dt matrix + `Das` log law RE-RUN on the corrected mesh (LOCAL) — all findings stand

Ryan directed a re-run after the two topological mesh defects were fixed (entry
below). All local, `julia --project -t 6`, `DirectBackend`. No HPC time.

**Trap handled first.** `_ssw_case_tag` encoded no mesh version and
`SSWETA_RESUME` defaults to `true`, so a plain re-run would have silently reused
every defective-mesh history and produced identical output. Two mitigations:
`SSWETA_RESUME=false` on the first pass, and a permanent
`const SSW_MESH_REVISION = 2` now appended to every case tag
(`examples/suddenly_started_wing.jl`). Verified by `CLsteady`:
**0.38018584 (old mesh) → 0.38005447 (corrected)** on all 24 matrix rows.

Also found: the pre-existing `data/ssw_eta_convergence/analysis.csv` had been
clobbered down to a single row by the old η=8/16 probe invocation, so 014's F1
ratios were hand-derived from a console table with **no surviving artifact**.
Fixed by (a) a new `das_analysis.csv` writer in `examples/ssw_eta_convergence.jl`
computing the fixed-`Das` diagonals, fixed-η ladders and log-law fits directly,
and (b) running the probes before the matrix so the matrix's summary survives.
New `SSWETA_LOGLAW_EXTRA` env folds out-of-matrix probes into the fit only.

```bash
# (a) large-Das probes, main root, run FIRST (~10 s)
SSWETA_ETAS=8,16 SSWETA_DTS=0.25 SSWETA_MESHCHECK=false SSWETA_RESUME=false \
  julia --project -t 6 examples/ssw_eta_convergence.jl
# (b) 24-case matrix + mesh check (~35 min)
SSWETA_RESUME=false SSWETA_LOGLAW_EXTRA=8:0.25,16:0.25 \
  julia --project -t 6 examples/ssw_eta_convergence.jl
# (c) long-window addendum, t*=20, ns=48 (~62 min; 1226-1244 s per case)
SSWETA_ETAS=0.25,1,4 SSWETA_DTS=0.0625 SSWETA_NSPAN=48 SSWETA_TEND=20 \
  SSWETA_MESHCHECK=false SSWETA_RESUME=false \
  SSWETA_OUTPUT=data/ssw_eta_convergence/long_fine \
  julia --project -t 6 examples/ssw_eta_convergence.jl
```

**Verdict: F1–F4 and the addendum all stand; numbers moved <1%.**
F1 fixed-`Das` diagonal ratios 0.39–0.47 (was 0.36–0.44); the sharpest statement
is new — the tail η-spread is dt-invariant (2.637% at dt*=1/4 vs 2.692% at
dt*=1/32). F2 dt-converged log slope **+0.205%/doubling, RMS resid 0.0035%** (was
a hand-derived +0.23%). F4 rotor +36.8% vs a predicted +0.48% over the same 2.36
`Das`-doublings — ~77×. Addendum tails 0.992681 / 0.996590 / 0.997892 at
`Das` = c/64 / c/16 / c/4 (was 0.9930 / 0.9966 / 0.9978).

One statement corrected rather than confirmed: the earlier "+0.23%/doubling from
0.004c to 4c" was read off a fixed-dt row, along which η — and hence the dt
error — also varies. The per-doubling increment along that row actually runs
1.11% → 0.21%, and a naive fit to it has RMS residual 0.303% against a 0.396%
slope. The log law is a dt-converged statement. Full detail in 014.

Not re-run: F2's wake-core robustness probes. `wake_core_over_c` has no ENV hook
and is absent from the case tag, so core variants collide in the resume cache.

### 2026-07-29 — Wake-representation consistency study (LOCAL, not HPC) — CONVERGED; hypothesis eliminated

Executed `plans/20260729_wake_representation_consistency.md`. All runs local,
`julia --project -t 6`, `DirectBackend` (mandatory: FMM inflates the wake-panel
radius by `|Das|`, `src/FLOWPanel_liftingbody.jl:996`, which would sweep its own
error alongside the variable of interest). No HPC time used.

**Code:** `src/FLOWPanel_wake.jl` — `PanelWake` gains `freestream_convection::Bool`
(default `false`, behaviour-preserving) plus a new leading branch in
`propagate!` translating all rows by `dt·freestream`; manifest round-trip in
`src/FLOWPanel_replay.jl`; unit test in `test/runtests_unit_wake.jl`
("PanelWake freestream convection mode", 8 assertions) covering exact rigid
translation, untouched rows beyond `nwakes[]`, and both legacy modes unchanged.
`test/runtests_unit_wake.jl` 90/90 pass.

**Driver:** `examples/suddenly_started_wing.jl` — `SSW_MODE=wake_consistency`,
`SSW_FREESTREAM_CONVECTION`, `SSW_WAKE_LENGTHS`, `SSW_WAKE_MODES`. Case tag gains
`_t<t_end>` / `_fc` suffixes only when non-default, so pre-existing histories
still resume. Mode defaults to AR=6 and `SAVE_VTK=false`.

Command:
```
SSW_MODE=wake_consistency SSW_BACKEND=direct SSW_WAKE_LENGTHS=1,2,4,8,16,32,64 \
  SAVE_VTK=false SSW_VERBOSE=false julia --project -t 6 \
  examples/suddenly_started_wing.jl
```
Elapsed: 3.4 / 0.5 / 1.3 / 3.7 / 12.3 / 55.4 / 292 s per case (rolled-up) and
0.21 / 0.48 / 1.13 / 3.28 / 11.5 / 51.3 / 285 s (straight sheet); whole sweep
under 15 min. Data: `data/suddenly_started_wing/wake_consistency.csv` and
`wake_consistency.png`; AR=3 cross-check in `data/suddenly_started_wing_ar3/`.

**Result:** the freestream-convected (straight-sheet) `PanelWake` converges to
the semi-infinite reference `cl=0.42122` monotonically — rel err −19.49% →
−0.0136% over `L/c` 1 → 64, asymptotic order ≈2.2 (successive ratios 0.213 /
0.220 / 0.236). The rolled-up control plateaus at −0.96%, correctly, since a
drooped wake is a different geometry. AR=3 reproduces it (order ≈2.1).
**The two wake representations agree; this candidate for the +37% `Das`
sensitivity is eliminated.** Full verdict and consequences in
`BRAINSTORM/014_first_wake_row_offset_selection.md`.

**Two mesh defects found and fixed en route** (both topological — symmetric node
cloud, asymmetric triangulation): TE panels were not chordwise reflections
(upper/lower TE control points 0.8% chord apart), and **all 480 cells violated
y→−y symmetry**, giving a ~10% spanwise asymmetry in shed circulation that did
NOT refine away (1.01e-1 / 1.13e-1 / 1.10e-1 at n_span 12/24/48). Fixed via an
XOR of the chordwise and spanwise half-tests; asymmetry now 1.3e-14. Guard
`assert_ssw_mesh_symmetry` added. **`examples/ssw_eta_convergence.jl` shares this
mesh**, so the F1–F4 η×dt matrix and `Das` log law in 014 were computed on the
defective mesh — they were re-run the same day, see the next entry.

### 2026-07-25 — Phase 2e opened

Ryan directed Phase 2e: converged unsteady hover CT on the Phase 2d production
mesh with `nwakerows=1`, with and without Green reconstruction, converging
timestep (`NT`) and wake truncation depth. HPC budget as directed: 64 threads,
24 GB per job. Phase docs written; dashboard updated.

### 2026-07-25 — Memory request raised 24 GB → 64 GB (measured)

The plan estimated ~16.6k panels for `dji9443_20260725_45_185_capped_captess4.msh`
(45 span × ~184 chordwise × 2 blades). That counts **quads**; the mesh is
triangulated, so the body actually has **36,752 panels** (18,380 nodes). Dense
operator footprint at that size, from `src/FLOWPanel_solver.jl:221-233` and
`src/FLOWPanel_formulation.jl:494-525`:

| Operator | Size |
|---|---:|
| `Backslash` G, `lu!` in place (one allocation) | 10.06 GB |
| `GreenReconstruction(:area_mean)` build: B **and** bordered K live simultaneously | 20.13 GB |
| … retained after the build (K/LU only; B becomes collectible on return) | 10.06 GB |
| **velocity case, steady** | **~10 GB** |
| **green case, peak during `initialize_formulation`** | **~30 GB** |

24 GB would OOM every Green case and leaves a velocity case with little
headroom for the particle wake and FMM trees. The launcher therefore requests
`--mem=64G`; a velocity-only case may be submitted leaner with
`sbatch --mem=24G`. Both dense factorizations are built **once**
(`initialize_formulation` is called before the time loop,
`src/FLOWPanel_simulate.jl:739`), so the per-step Green cost is only an
O(N²) triangular solve — cheap. `green_solver=nothing` (dense) is therefore the
right starting route; a `KrylovSolver` route would trade that 10 GB for
per-step matrix-free products and is only worth it if memory, not time, binds.

### 2026-07-25 — Local preflight complete; deployment blocked on interactive SSH

`test/formulation_test.jl` passes in full. Both formulations run end-to-end
through the modified driver on the cheap 40_40 mesh. The ct4 production body
constructs cleanly with auto-detected TE seeds and **zero unpaired shedding
edges**, which is the condition `_validate_formulation_common`
(`src/FLOWPanel_formulation.jl:342-346`) warns on for Green.

`ssh -o BatchMode=yes rander39@ssh.rc.byu.edu` returns
`Permission denied (keyboard-interactive)` — the cluster requires interactive
2FA, so no agent can scp, md5-verify, or inspect `squeue`. Deployment is handed
to Ryan via `scripts/deploy_phase02e.sh`.

Local md5 of the assets to deploy (record the remote values beside these after
`deploy_phase02e.sh` runs):

| md5 (local) | file | status |
|---|---|---|
| `9fab48aa70f196be58a9f06d4d4291c3` | `examples/data/dji9443_20260725_45_185_capped_captess4.msh` | untracked → **must be copied** |
| `0b073ffc082f1f251e3ddf5511424bab` | `examples/rotor_hover_pressure_comparison.jl` | modified → **must be copied** |
| `4ef1ce01bdf65a598ce97f01e15afb60` | `examples/run_dji9443_hover_ct_hpc.slurm.sh` | new → **must be copied** |
| `0f9017b72f56e0ea5e5b79edf3cbb4e7` | `examples/dji9443_trailing_edge.jl` | clean vs HEAD → verify only |
| `156c2ae79df53f8a7376691633c5101e` | `src/FLOWPanel_formulation.jl` | clean vs HEAD → verify only |
| `46f5042323ec22565c59b877e6a670b6` | `src/FLOWPanel_simulate.jl` | clean vs HEAD → verify only |
| `704dd2f56720e37330352dba14bdb216` | `src/FLOWPanel_liftingbody.jl` | clean vs HEAD → verify only |

No `src/` file was modified by this phase. The verify-only entries exist
because of the Phase 2c stale-cluster-src monitor bug: if any differs, reconcile
the cluster checkout before submitting.

### 2026-07-25 02:39 — Batch 1 deployed and SUBMITTED

Ryan granted explicit permission to submit and pointed at the working SSH
alias: **`ssh orc`** reaches `login04` non-interactively (the earlier failure
was using the raw `rander39@ssh.rc.byu.edu` host, which demands 2FA). Recorded
here so no future agent repeats the wrong host.

Deployment verification (local vs cluster md5, all seven match):

| md5 | file | action |
|---|---|---|
| `9fab48aa70f196be58a9f06d4d4291c3` | `examples/data/dji9443_20260725_45_185_capped_captess4.msh` | already on cluster; verified |
| `0b073ffc082f1f251e3ddf5511424bab` | `examples/rotor_hover_pressure_comparison.jl` | copied |
| `4ef1ce01bdf65a598ce97f01e15afb60` | `examples/run_dji9443_hover_ct_hpc.slurm.sh` | copied |
| `0f9017b72f56e0ea5e5b79edf3cbb4e7` | `examples/dji9443_trailing_edge.jl` | verified only |
| `156c2ae79df53f8a7376691633c5101e` | `src/FLOWPanel_formulation.jl` | verified only |
| `46f5042323ec22565c59b877e6a670b6` | `src/FLOWPanel_simulate.jl` | verified only |
| `704dd2f56720e37330352dba14bdb216` | `src/FLOWPanel_liftingbody.jl` | verified only |

Cluster checkout on branch `fastmultipole` at `5615ada` — same commit as local,
and no `src/` file is modified locally, so the Phase 2c stale-cluster-src
failure mode is excluded. `logs/slurm/` created before submission (Slurm opens
log paths before the script runs). `squeue -u rander39` was empty, so all three
Batch 1 jobs fit under the 3-job study cap.

Submitted commands (from `/home/rander39/projects/FLOWPanel.jl`):

```
sbatch --job-name=fp-p2e-vel-nt36-d4   examples/run_dji9443_hover_ct_hpc.slurm.sh p2e_vel_nt36_d4
sbatch --job-name=fp-p2e-green-nt36-d4 examples/run_dji9443_hover_ct_hpc.slurm.sh p2e_green_nt36_d4
sbatch --job-name=fp-p2e-vel-nt72-d4 --time=48:00:00 \
       examples/run_dji9443_hover_ct_hpc.slurm.sh p2e_vel_nt72_d4
```

→ jobs **12894164 / 12894165 / 12894166**, all RUNNING on partition `m12`,
account `sn72`. Launcher banners confirm the intended settings per case
(mesh `45_185_ct4`, RPM 5400, correct formulation/NT/depth/rlxf). Empty `.err`
files at first check.

Note: `squeue`/`sbatch` are not on the non-interactive `PATH`; wrap remote
Slurm calls in `bash -lc "..."`.

### 2026-07-25 02:48 — Job cap raised 3 → 6; full required matrix launched

Ryan raised the concurrent-job cap from three to six. Updated in the dashboard
global directives, the Phase 2e plan, and the launcher header.

The three added cases are **not adaptive choices** — they are the rungs the
Phase 2e decision rule requires regardless of Batch 1's outcome, since
parameter convergence is evaluated *separately per formulation*:

| Job | Case | Role |
|---|---|---|
| 12894168 | `p2e_green_nt72_d4` | NT 36→72 refinement, **green** (partner of 12894166) |
| 12894169 | `p2e_vel_nt36_d6` | depth 4R→6R refinement, **velocity** |
| 12894170 | `p2e_green_nt36_d6` | depth 4R→6R refinement, **green** |

With these, all six required matrix cells are in flight; the 2×2×2 design is
complete (formulation × NT × depth, sharing the nt36_d4 baselines). Genuinely
adaptive items — warm-start settle extension, the Krylov-Green cost check, a
depth-8R probe — remain deferred to Batch 2 and now have no free slots, so they
wait for a completion.

Wall-time allocation reflects expected cost rather than a flat default:
NT=72 doubles the step count (48 h), and depth 6R retains ~1.5× the particles
of 4R for ~1.5× the per-step FMM cost (36 h).

Launch-safety note: this went out ~9 minutes into Batch 1, before any Julia
stdout had appeared (package load + the 1.35e9-pair G assembly + LU precede the
first print, and stdout is block-buffered to file). That was judged acceptable
because all three new cases reuse the identical driver, mesh, and launcher —
any systemic failure would surface in Batch 1 first and `scancel` makes the
extra jobs cheaply reversible. Jobs 12894169/12894170 landed on the same node
(m12-2-14); Slurm enforces the 64-core/64 GB allocations, but shared memory
bandwidth may slow both relative to the solo jobs — compare wall times before
reading anything physical into them.

### 2026-07-25 03:07 — Health check (step rate RETRACTED — see the 03:53 entry)

**The step-rate figures in this entry are wrong by ~3× and the "~1 h wall"
conclusion is withdrawn.** `ls data/<case>/<case>_body1/*.vtu | wc -l` counts
about **3 files per step index**, not one, so the counts below are file counts,
not step counts. Corrected progress and cost are in the 03:53 entry. The
health conclusion (all six alive, no errors) and the Green-vs-velocity cost
ratio (a ratio, so unaffected by the scaling error) both stand.

Original entry follows.

#### 2026-07-25 03:07 — Health check: all six healthy

All six jobs RUNNING, all `.err` files 0 bytes. Step counts from the on-disk
body VTU series (the reliable progress signal — Julia stdout is block-buffered
to file, so the `.out` files still showed only the launcher banner at 27 min,
which is **not** evidence of a stall):

| Job | Case | Elapsed | Body steps written |
|---|---|---:|---:|
| 12894164 | `p2e_vel_nt36_d4` | 27 min | 354 |
| 12894165 | `p2e_green_nt36_d4` | 27 min | 342 |
| 12894166 | `p2e_vel_nt72_d4` | 27 min | 402 |
| 12894168 | `p2e_green_nt72_d4` | 19 min | 294 |
| 12894169 | `p2e_vel_nt36_d6` | 19 min | 273 |
| 12894170 | `p2e_green_nt36_d6` | 19 min | 255 |

Green costs no more than velocity per step at NT=36 d4 (342 vs 354 *files* in
the same 27 min, ≈3% apart) — a ratio, so unaffected by the file-vs-step
scaling error. The dense Green route's per-step cost really is just the O(N²)
triangular solve, as predicted. No case for switching to a `KrylovSolver` route
on cost grounds.

A seventh unrelated job (12894176, `fm024bcpu`) is also in Ryan's queue; it is
not a Phase 2e case and does not count against the study cap.

Methodology note: an earlier check piped the remote output through `tail -40`,
which silently clipped the `squeue` header and three of the six job rows and
made it look as though half the batch had vanished. Count jobs explicitly
(`squeue --noheader | wc -l`) rather than reading a truncated listing.

### 2026-07-25 03:53 — Corrected progress and cost; wall-time requests are appropriate

All six still RUNNING, all `.err` files 0 bytes, ~1:13 elapsed for the first
three and ~1:04 for the second three.

**Correct way to read progress.** Use the highest step *index*, not the file
count: the `<case>_body1/` directory holds `.vtm` + multiple `.vtu` pieces per
step (~3 `.vtu` per index, exactly `3*(maxidx+1)`), so a raw `*.vtu | wc -l`
overstates progress ~3×.

```
maxidx=$(ls data/<case>/<case>_body1/*.vtu | sed 's/.*\.\([0-9]*\)\.vtu/\1/' | sort -n | tail -1)
```

The `.out` files lag by only ~20 steps (block buffering), so `step N/total` in
the log is a fine secondary check — it agreed with `maxidx` to within 20.

| Job | Case | Elapsed | Step (maxidx) | of | Fraction |
|---|---|---:|---:|---:|---:|
| 12894164 | `p2e_vel_nt36_d4` | 1:13 | 219 | 719 | 30% |
| 12894165 | `p2e_green_nt36_d4` | 1:13 | 212 | 719 | 29% |
| 12894166 | `p2e_vel_nt72_d4` | 1:13 | 258 | 1439 | 18% |
| 12894168 | `p2e_green_nt72_d4` | 1:04 | 235 | 1439 | 16% |
| 12894169 | `p2e_vel_nt36_d6` | 1:04 | 201 | 719 | 28% |
| 12894170 | `p2e_green_nt36_d6` | 1:04 | 200 | 719 | 28% |

**The step rate is decaying as the wake fills**, which the earlier single-point
estimate could not see. For `p2e_vel_nt36_d4`: 118 steps in the first 27.5 min
(4.3 steps/min), then 102 steps in the next 46 min (2.2 steps/min). Expected —
particle count grows until the truncation cylinder reaches steady state, and
FMM cost tracks it.

Revised projections (remaining steps at a continued-decay rate ~1.5–2/min):

| Case family | Projected wall | Requested | Verdict |
|---|---|---|---|
| NT=36, d4 | ~6 h | 24 h | comfortable |
| NT=36, d6 | ~9 h (≈1.5× particles) | 36 h | comfortable |
| NT=72, d4 | ~24 h (2× steps × ~2× particles/rev) | 48 h | adequate, watch it |

So **the original wall-time requests were appropriate, not 20× conservative**.
NT=72 deserves monitoring: at NT=72 the shedding rate per revolution doubles,
so the steady-state particle count is ~2× the NT=36 case and per-step cost
roughly doubles on top of the doubled step count. If either NT=72 case
approaches its limit, resume with `RESTART_STEP` rather than a cold rerun.

Mid-transient CT values in the logs (step ~199, rev ~5.5) are inside the
freestream withdraw phase, so they are **not** hover CT and must not be
compared across cases yet — the buffered tails also show different cases at
different steps, which would make any such comparison doubly meaningless.

### 2026-07-25 04:55 — Rate decay tracked; all six on course

All six RUNNING, `.err` all 0 bytes, no case finished yet (2:15 / 2:07 elapsed).

| Case | step @03:53 | step @04:55 | of | rate over the hour |
|---|---:|---:|---:|---:|
| `p2e_vel_nt36_d4` | 219 | 306 | 719 | 1.40 /min |
| `p2e_green_nt36_d4` | 212 | 300 | 719 | 1.42 /min |
| `p2e_vel_nt36_d6` | 201 | 295 | 719 | 1.52 /min |
| `p2e_green_nt36_d6` | 200 | 295 | 719 | 1.53 /min |
| `p2e_vel_nt72_d4` | 258 | 380 | 1439 | 1.97 /min |
| `p2e_green_nt72_d4` | 235 | 354 | 1439 | 1.92 /min |

Decay for `p2e_vel_nt36_d4`: 4.3 → 2.2 → 1.4 steps/min. It should **flatten
soon**: step 306 = rev 8.5, and the 4.5R-deep truncation cylinder fills in
roughly 6–12 revs at the hover convection speed, after which particle count —
and therefore FMM cost — is capped. Treat the projections below as valid only
if the next interval confirms the flattening.

Projected completion (assuming the rate settles near 1.2–1.4 /min):
NT=36 cases ~5–6 h more (total ~8 h, limit 24 h); NT=72 ~17 h more if its rate
settles near 1.0 /min (total ~20 h, limit 48 h). NT=72 remains the only case
without comfortable margin.

Green-vs-velocity per-step cost stays within ~2% at matched NT and depth
(300 vs 306, 295 vs 295, 354 vs 380 — the last pair started 9 min apart).

### 2026-07-25 08:15 — Batch 1 truncated; TRUNCATION DEPTH ANSWERED (null)

Ryan cancelled the NT=72 pair (revisit later) and dropped the truncation-depth
question, freeing 4 slots. `scancel 12894166 12894168 12894169 12894170`.

Before deleting anything, per-rev CT was harvested from the `.out` logs of all
six Batch 1 jobs into `data/p2e_harvest/*_per_rev.csv` on the cluster.

**Truncation depth 4R vs 6R, production mesh, 5400 RPM, matched revs:**

| rev | vel d4 | vel d6 | Δ | green d4 | green d6 | Δ |
|---:|---:|---:|---:|---:|---:|---:|
| 7 | 0.06430 | 0.06428 | −0.02% | 0.06514 | 0.06512 | −0.03% |
| 8 | 0.06093 | 0.06107 | +0.22% | 0.06181 | 0.06181 | +0.00% |
| 9 | 0.06078 | 0.06087 | +0.16% | 0.06167 | 0.06171 | +0.05% |
| 10 | 0.06284 | 0.06284 | +0.01% | 0.06383 | 0.06392 | +0.13% |

**|ΔCT| ≤ 0.22%** — the documented 6000-RPM/40_40 null reproduces on the
production mesh. Deeper truncation does not raise thrust. Question closed with
our own data rather than only the prior.

**Unplanned finding — NT may be an upward lever.** The killed NT=72 run reached
**rev 7 at CT 0.0696** vs NT=36's **0.0643** at the same rev (**+8%**). Still
transient (freestream withdraw ends rev 6.5), so not conclusive, but it points
*opposite* to the 005 `nt_hi` result (NT 18→36 gave −6.6%). Worth revisiting;
the NT=72 run directories were **kept intact** (12.3 GB) precisely so a return
can warm-start from ~rev 8 rather than pay the cold fill again.

Retention: the two d6 run directories were deleted (9.2 GB, question dead, CT
harvested). Cluster disk is 192 G / 2.0 T used, so no pruning pressure. Current
full histories = 2 running baselines + 2 archived NT=72 + 4 new Batch 2 = 8,
above the study's ≤4 cap; trim at analysis time, keeping the CT CSVs.

### 2026-07-25 08:55 — Batch 2 submitted: relaxation filter × σ × Das (2×2)

**σ was measured, not assumed.** From job 12894164 step 481 (56,626 particles),
particle `sigma/R` = **0.022 (min) / 0.145 (median) / 0.379 (max)** — the tip
vortex is smeared over a third of the radius.

Mechanism: `_shed_particles!` (`src/FLOWPanel_wake.jl:618-622`) sets
`sigma = dist·overlap/p_per_step`, and the shed segment spans node row 1 → row 2
(`src/FLOWPanel_simulate.jl:1248-1251`). Because `update_TE!` puts row 1 at
TE+`Das` and `shed_wake!` shifts row 1 → row 2 each step, **the `Das` offset
cancels between the rows** and `dist` ≈ one step's TE travel (≈2πr/nt), which is
why σ grows with radius. **σ is therefore independent of `Das`** — an earlier
claim in this session that η=4 would give σ=1.05R was wrong and is retracted.

Ryan's targets: σ ≈ 0.01R, `overlap` ≥ 2. Cost of σ=0.01R: overlap 3→2 buys a
free 1.5×, the remaining ~10× needs `p_per_step ≈ 20` ⇒ ~560k particles (over
`max_particles=500_000`) and ~10× FMM cost ≈ 4–5 days/job. Not reachable.
Selected instead: **overlap 2 + p_per_step 8** ⇒ σ_med ≈ 0.024R (6× finer),
~230k particles, ~4× cost, 48 h wall. `MERGE_R_FACTOR` cut 0.02 → 0.0025 so
merging cannot undo the refinement.

Design is a **2×2 factorial** in (σ, η) at filter 2.0R, so every comparison is
single-variable — corrected from an earlier draft whose "σ effect" contrast
differed in both σ and η:

- 12894470 σ stock, η 0.2 — **anchor**, gives the pure filter effect vs baseline 12894165
- 12894472 σ fine, η 0.2 — σ effect (vs 12894470)
- 12894471 σ stock, η 4.0 — η effect (vs 12894470)
- 12894473 σ fine, η 4.0 — interaction

**Blow-up contingency (Ryan, standing instruction):** if a σ-fine case diverges,
resubmit it with `OVERLAP=4` (σ_med ≈ 0.048R, still 3× finer than stock). The
launcher honors an inherited `OVERLAP`, so
`OVERLAP=4 sbatch ... p2e_green_f2p0_sigF_das0p2` is sufficient.

Driver changes (md5 `437288127018ad6c518d1e7ba55bd5c6`, verified on cluster):
new envs `DAS_ETA_KINEMATIC` / `DAS_MIN_DISPLACEMENT_R` /
`SHED_WITH_INDUCED_VELOCITY` (all defaulting to the previous hard-coded values);
per-rev CSV now covers **all** complete revolutions with an
`in_convergence_window` flag; `CT_cycle_mean` / `CT_cycle_std` /
`CT_cycle_std_rel` added as the headline metric for limit-cycle configurations.

Local validation (40_40, NT=6, 3 revs): filter constructed with
`point = [0.238, 0, 0]` (= 2.0R) and `normal = [1, 0, 0]` ⇒ relaxation **off**
for `0 ≤ x < 2.0R`, on beyond — the requested behavior; all-revs CSV and all new
metadata fields written. Two earlier local smokes died to memory contention from
an unrelated Julia process on Ryan's laptop (a second session), not to a code
fault — the heavy σ-fine config is exercised on the cluster instead.

### 2026-07-25 08:55 — Monitoring gotcha: the login MOTD contaminates scalars

A divergence alarm fired seconds after Batch 2 launched, reporting
"DIVERGENCE DETECTED in 10140316000 sigma-fine log(s)". **False alarm, caused by
this log's own tooling, not by the jobs.** The monitor extracted a count with
`ssh ... | tr -dc '0-9'`, which concatenates every digit in the BYU login
**MOTD banner** (announcement URLs, version numbers) into the "count". All four
`.err` files were 0 bytes and the jobs were ~20 s old.

This is the second time the MOTD has corrupted a remote reading (the first was a
`tail -40` that clipped the `squeue` header and three job rows, making it look
as though half of Batch 1 had vanished). **Rule for this study: never parse a
bare `ssh` result. Emit a sentinel and extract it**, e.g.

```
ssh orc 'bash -lc "... && echo MARK_N=$(squeue -u rander39 --noheader | grep -c fp-p2e)"' \
  | grep -o 'MARK_N=[0-9]*' | cut -d= -f2
```

Also: a heredoc `ssh orc 'bash -s'` does **not** get a login shell, so `squeue`
is missing from `PATH`; either use `bash -lc` or
`export PATH=/apps/slurm/latest/bin:$PATH` inside the script.

Verified with sentinel parsing at 08:58: `MARK_N=6`, `MARK_BAD=0` — all six jobs
healthy, no divergence.

**Then it happened a third way.** The same sentinel probe wrapped in a shell
*function* returned `MARK_N=0` while the very same output block listed six
RUNNING jobs — `$(...)` collapsing somewhere in the nested
single-quote/double-quote/`bash -lc` layering. Three different corruption modes
in one session is enough: all remote parsing now lives in a script **on the
cluster**, `scripts/p2e_status.sh`, invoked bare:

```
ssh orc /home/rander39/projects/FLOWPanel.jl/scripts/p2e_status.sh | grep '^P2E_'
```

It emits only sentinel-prefixed lines (`P2E_NJOBS=`, `P2E_NBAD=`, `P2E_JOB=`,
`P2E_STEP=`, `P2E_CT=`), so nothing the caller reads can be contaminated by the
MOTD or by quoting. Do not go back to inline `ssh '...'` pipelines for status.

### 2026-07-25 09:08 — σ-fine cases OOM at 64 GB; resubmitted at 256 GB

`sacct`: 12894472 `State=OUT_OF_MEMORY`, `ExitCode 0:125`, elapsed 00:14:14,
`ReqMem 64G`; `.err` carries
`slurmstepd: error: Detected 1 oom_kill event`. Died at step 56 of 954.

**Diagnosis — it is not a stability blow-up and not particle count.** Particle
counts read from the VTP headers:

| case | pps | step | particles |
|---|---:|---:|---:|
| `sigF_das0p2` (OOM) | 8 | 56 | **33,519** |
| `sigC_das0p2` (healthy) | 2 | 80 | 11,500 |
| `sigF_das4p0` | 8 | 30 | 17,386 |
| `vel_nt36_d4` baseline (healthy) | 2 | 527 | **53,999** |

The killed job held **fewer** particles (33.5 k) than the healthy baseline
(54.0 k), so per-particle storage cannot explain it — implied ~1.3 MB/particle
is absurd. Something in the σ-fine configuration allocates non-linearly
(suspect: an FMM near-field or merge-stage transient that grows as σ shrinks).
Not chased further; the practical fix is headroom.

**Ryan's `OVERLAP=4` fallback deliberately NOT applied here.** That instruction
targeted a *stability* blow-up. This is memory, and with `OverlapPPS`,
`p_per_step` fixes the particle count while `overlap` only scales σ — so
overlap 2→4 would not reduce memory, it would merely coarsen σ 0.024R → 0.048R
and surrender half the refinement Ryan asked for. Node `m12-1-8` has
**RealMemory = 512 GB** and we had requested only 64 GB, so this was simply an
under-request.

Action: `scancel 12894473` (same config, ~17 k particles and climbing — would
have hit the same wall and wasted a 48 h slot), then both σ-fine cases
resubmitted **unchanged except `--mem=256G`**:

```
sbatch --job-name=fp-p2e-f2-sigF-das0p2 --time=48:00:00 --mem=256G \
       examples/run_dji9443_hover_ct_hpc.slurm.sh p2e_green_f2p0_sigF_das0p2
sbatch --job-name=fp-p2e-f2-sigF-das4p0 --time=48:00:00 --mem=256G \
       examples/run_dji9443_hover_ct_hpc.slurm.sh p2e_green_f2p0_sigF_das4p0
```
→ 12894481 / 12894482, both RUNNING; queue back to 6.

**Escalation if 256 GB also fails.** Memory is evidently not linear in particle
count, and the count still has to grow ~7× before the truncation cylinder
saturates (~rev 8). If these OOM again, drop `P_PER_STEP` 8 → 4 (σ_med 0.048R,
half the particles) rather than raising memory a third time — and only then
consider `OVERLAP=4`.

Also resolved: the η=4.0 cases that had written no step at the 5-minute check
were merely slower to start, not stalled — 12894470/12894471 reached steps
69/28 by 09:08. No action taken.

### 2026-07-25 09:33 — CORRECTION: the σ-fine failure IS a stability blow-up

**The 09:08 diagnosis above was wrong and is retracted.** I concluded "not a
stability blow-up, just an under-request" from the absence of NaN/Inf and from
the OOM'ing job holding fewer particles than a healthy baseline. Both
observations were real; the inference was not.

The 256 GB resubmit (12894481) failed in 16 min with `ExitCode 1:0`, and its
last force print was:

```
CF = (53.46, 0.1716, 0.004559)        # CT ~ 53, versus a normal ~0.06
```

Max |CT| over that log is **90.17**. The run had already diverged numerically;
the exception was `OutOfMemoryError()` raised inside
**`merge_particles!` (`FLOWVPM_merging.jl:454`, a `resize!` of a hash bucket)**
— the exploded particle cloud wrecked the merge spatial hash. So the original
64 GB OOM was the *same* divergence, one stage earlier. Memory was a symptom
throughout; **`OVERLAP=2` is simply unstable for this wake.**

Why the detector missed it: divergence happened in large-but-**finite**
arithmetic, and the monitor only grepped for `NaN|Inf`. `scripts/p2e_status.sh`
now also flags on magnitude (hover CT is O(0.06), so **|CT| > 1 ⇒ diverged**)
and emits `P2E_DIVERGED=<log> reason=...`. Verified: it flags 12894481.

**Ryan's `OVERLAP=4` fallback was the correct instruction all along** — this is
exactly the stability failure it was meant for. Applied now as the launcher
default for σ-fine cases:

| | σ_med/R | vs stock | particles | status |
|---|---:|---|---:|---|
| stock (ov 3, pps 2) | 0.145 | 1× | ~54 k | healthy |
| ov 2, pps 8 | 0.024 | 6× finer | ~230 k | **DIVERGED** |
| **ov 4, pps 8** | **0.048** | **3× finer** | ~230 k | resubmitted |

`overlap` scales σ only; `p_per_step` still fixes the particle count, so the
resolution in particle *spacing* is unchanged — only the core overlap doubles,
which is what buys stability. `MERGE_R_FACTOR` 0.0025 → **0.0066** to hold
`r_merge/σ` at the stock ratio (0.02/0.145 = 0.138).

Action: `scancel 12894482` (same unstable config), resubmit both σ-fine cases
with the new launcher (md5 `a27dad12423aac9fcf4ef6f36af6ba65`, verified) at
`--mem=128G` → **12894486 / 12894487**, both RUNNING. Queue back to 6.

Standing note for the monitor: `P2E_NBAD` has a **baseline of 1** from the
retained 12894481 log; only a rise above that indicates a new divergence.

**If ov 4 also diverges**, the next lever is `P_PER_STEP` 8 → 4 (σ_med 0.048R at
overlap 2, or 0.097R at overlap 4) — i.e. back off the refinement rather than
push overlap higher, since overlap ≥4 starts to defeat the purpose of refining σ.

### 2026-07-25 10:11 — OVERLAP=4 holding; early Das signal is large

All 6 RUNNING, no new divergence (`P2E_NBAD=1` is the retained 12894481 log).

**The stability fix is working.** σ-fine at overlap 4 has reached steps **168**
(12894486) and **118** (12894487); overlap 2 diverged and died at step **56–58**.
Passing 3× that point without a magnitude excursion is good evidence the wake is
stable at σ_med ≈ 0.048R.

Baseline per-rev CT, revs 10–15 (filter OFF, both at step ~610/719):

| rev | 10 | 11 | 12 | 13 | 14 | 15 |
|---|---:|---:|---:|---:|---:|---:|
| velocity | 0.06284 | 0.06359 | 0.06295 | 0.06194 | 0.06088 | 0.06162 |
| green | 0.06383 | 0.06449 | 0.06389 | 0.06302 | 0.06170 | 0.06214 |

Green holds a steady **+1.5%** over velocity across all six revs. Both oscillate
in the 0.061–0.064 band rather than flattening — i.e. **the unfiltered
configuration also looks like a limit cycle, not a plateau**, which supports
using the cycle-mean as the headline metric for every case, not just filtered
ones.

**Early `Das` signal (transient revs 0–4, NOT hover — do not quote as a result):**

| rev | 0 | 1 | 2 | 3 | 4 |
|---|---:|---:|---:|---:|---:|
| η=0.2 (12894470) | 0.06519 | 0.04812 | 0.04565 | 0.05042 | 0.06274 |
| η=4.0 (12894471) | 0.08181 | 0.05561 | 0.05116 | 0.05725 | 0.07184 |
| Δ | +25.5% | +15.6% | +12.1% | +13.5% | **+14.5%** |

η=4 sits **12–25% above** η=0.2 at every matched transient rev — larger than the
relaxation-filter effect (+7%) and larger than the oracle-wing analogue (+13%
over a 4× `das` range). These revs are inside the freestream ramp/hold/withdraw,
so they are not hover CT and the magnitude will change; but the *direction and
scale* are already clear. If this survives into the settle window it makes
first-wake-row placement the **dominant** discretization uncertainty in the
whole unsteady model — bigger than mesh, formulation, or relaxation — and the
correct conclusion would be that `Das` needs convergence study, not that a
particular η validates against experiment.

### 2026-07-27 — Cluster access lost; Batch 2 state as of the last successful probe

`ssh orc` now fails with `Permission denied (keyboard-interactive)`. Root cause
is **not** a credential change: `Host orc` multiplexes over a ControlMaster
socket, that socket has expired, and the underlying auth is 2FA. My earlier
statement in this log that "`ssh orc` works non-interactively" was true only
while the socket was live and has been corrected in the snapshot above.
**Nothing further can be read or submitted until Ryan runs `! ssh orc -fN`.**

Last successful probe was 2026-07-25 11:17. State at that moment:

| Job | Case | State | Note |
|---|---|---|---|
| 12894164 | `p2e_vel_nt36_d4` | RUNNING 8:38 | step 633/719 — should have completed |
| 12894165 | `p2e_green_nt36_d4` | RUNNING 8:38 | step 631/719 — should have completed |
| 12894470 | `f2_sigC_das0p2` | RUNNING 2:25 | step 311/1007 |
| 12894471 | `f2_sigC_das4p0` | RUNNING 2:25 | step 225/1007 |
| 12894486 | `f2_sigF_das0p2` (ov 4) | RUNNING 1:42 | step 188 — **past the step-56 divergence point** |
| **12894487** | `f2_sigF_das4p0` (ov 4) | **OUT_OF_MEMORY** at 1:35:32 | 128 GB |

**New finding: σ-fine + η=4 is the failing combination.** With `OVERLAP=4`, the
η=0.2 σ-fine case survived well past where `OVERLAP=2` died, but the η=4 σ-fine
case still OOM'd at 128 GB after 1.5 h. Since η=4 displaces the whole attached
wake row ~0.70R downstream, combining that with 4× the particles plausibly
spreads the cloud far enough to blow up the merge hash again. Whether this was
another divergence-driven OOM or a genuine memory ceiling **cannot be
determined without log access** — check `max|CT|` in
`logs/slurm/slurm-fp-p2e-f2-sigF-das4p0-12894487.out` first (the
`P2E_DIVERGED` magnitude test in `scripts/p2e_status.sh` will answer it).

**On resuming access, in order:**
1. `ssh orc /home/rander39/projects/FLOWPanel.jl/scripts/p2e_status.sh` — full state.
2. `sacct --starttime 2026-07-25T02:00 ... | grep p2e` — final states of all jobs.
3. Collect `data/p2e_*/*_case_metadata.toml` and `*_CT_per_rev.csv` — the
   completed baselines should now carry `CT_cycle_mean` / `CT_cycle_std`.
4. Diagnose 12894487 by magnitude before deciding whether to retry it at
   `P_PER_STEP=4` (the pre-agreed next lever) rather than more memory.

### 2026-07-27 — Batch 2 post-mortem: the filter breaks particle-count saturation

Access restored (Ryan re-established the ControlMaster socket). Final states:

| Job | Case | State | Elapsed |
|---|---|---|---|
| 12894164 | `p2e_vel_nt36_d4` | **COMPLETED** | 10:06:07 |
| 12894165 | `p2e_green_nt36_d4` | **COMPLETED** | 10:04:17 |
| 12894470 | `f2_sigC_das0p2` (64 GB) | OUT_OF_MEMORY | 04:37:25 |
| 12894471 | `f2_sigC_das4p0` (64 GB) | OUT_OF_MEMORY | 06:36:08 |
| 12894486 | `f2_sigF_das0p2` (ov 4, 128 GB) | OUT_OF_MEMORY | 02:52:54 |
| 12894487 | `f2_sigF_das4p0` (ov 4, 128 GB) | OUT_OF_MEMORY | 01:35:32 |

**RETRACTED 2026-07-27 — these OOMs ARE divergence.** This entry originally
read "these OOMs are NOT divergence", inferred from max|CT| of **0.1295**
(12894471) and **0.09747** (12894486) with final forces −0.073 / −0.068. That
inference was invalid: Julia's stdout is **block-buffered to file and a SIGKILL
discards the buffer**, so a killed run's log stops tens of steps *before* the
failure — a sane maximum over that truncated prefix says nothing about the end
state. Ryan's ParaView inspection of `p2e_OOM_sigF_das4p0` shows the wake blew
up, and the OOM historically throws from `merge_particles!` (confirmed by
12894481's stack trace at `FLOWVPM_merging.jl:454`, with CT already at 90).
Correct chain: **instability → particles scatter → merge hash allocation
explodes → OOM**. Memory is the symptom. Magnitude tests are only valid on
flushed (completed) logs.

**Particle count never saturates under the filter** (consistent with a
scattering, diverging cloud escaping both the truncation cylinder and the
merge — not, as first written, an independent memory ceiling):

| case | step | particles |
|---|---:|---:|
| baseline green (completed) | 719 | **43,906** — saturated |
| `sigC_das0p2` (filter) | 421 | 55,856 — still climbing |
| `sigC_das4p0` (filter) | 420 | 55,707 — still climbing |
| `sigF_das0p2` (filter, pps 8) | 219 | 118,417 — still climbing |

The unfiltered wake levels off at ~44 k because the `GlobalCylinder` sweeps
particles out. With relaxation disabled inside 2R the particles are never
realigned with the local flow, so they scatter and **linger inside the
cylinder** instead of convecting through it. Count keeps growing and the FMM
octree spans a larger, sparser volume. Prior filter work (BRAINSTORM/006) never
hit this because it ran on the 7,288-panel mesh with a far smaller wake.

**Retry (submitted 2026-07-27): 12915034/35/36/37.** Two changes:
1. `SETTLE_REVS=12` (719 steps) instead of 20 — matches the baselines exactly,
   so the filter comparison is like-for-like on run length, and still leaves
   ~12 revs of hover for a 10-rev cycle-mean.
2. Memory well above what failed: **256 GB** for σ-stock, **384 GB** for σ-fine
   (node `RealMemory` is 512 GB).

The launcher now honors an inherited `SETTLE_REVS`/`CONVERGENCE_REVS`
(md5 `c05c2e210d298b7694bafbc7e51cdc29`, verified on cluster), passed through
with `--export=ALL,SETTLE_REVS=12`.

**Note on the baselines' `converged=false` metadata:** those jobs started at
02:39, *before* the `CT_cycle_*` metadata fields were added at ~08:50, so they
carry only the old 5-rev window. The 10-rev cycle-means above were recomputed
from the full 20-rev `ForceMonitor` series in their logs.

### 2026-07-27 — DATA LOSS: three OOM'd histories destroyed by the retry launch

Ryan asked to view the OOM'd cases in ParaView. **Three of the four were already
gone.** The retry jobs 12915034/35/36 reuse the same `RUN_NAME`, and the
launcher's first action was `rm -rf "data/$RUN_NAME"` — so as each retry
started it wiped its own predecessor's output:

| case | OOM'd job | VTK status |
|---|---|---|
| `sigC_das0p2` | 12894470 | **destroyed** by 12915034 |
| `sigC_das4p0` | 12894471 | **destroyed** by 12915035 |
| `sigF_das0p2` | 12894486 | **destroyed** by 12915036 |
| `sigF_das4p0` | 12894487 | **PRESERVED** — its retry 12915037 was still PENDING |

This was avoidable and is on me: I wrote the comment documenting that `rm -rf`
and still resubmitted without moving the OOM'd output aside first, which would
have cost nothing.

Rescue of the survivor, done in the PENDING window:
`scontrol hold 12915037` → `mv data/p2e_green_f2p0_sigF_das4p0
data/p2e_OOM_sigF_das4p0_12894487` → `scontrol release 12915037`.
Contents verified: 129 body steps (0–128), PVD manifest with 129 entries, 129
particle `.vtp`, 383 filament files, **64,000 particles** at the final step.
Copied to the local repo `data/` for ParaView.

**What is NOT lost:** the full per-revolution CT histories of all three wiped
cases survive in `logs/slurm/*.out`, and the particle counts and max|CT| values
quoted in the 2026-07-27 post-mortem were already extracted from them. The
quantitative conclusions stand; only the wake *geometry* is gone.

**Retention policy (Ryan, 2026-07-27): "in general, you can wipe old data. We
don't have space to keep all paraview files."** The launcher therefore stays
**destructive by default** (md5 `05f00b3132f67d0d024a4c04e9fb2f71`), with
`RHPC_KEEP_PREV=true` as an opt-in that moves the old run to `data/<name>.prev`.

**Rule for this study: before resubmitting a case tag, harvest the CT series
from `logs/slurm/*.out` (cheap, always worth it); set `RHPC_KEEP_PREV=true`
only when the previous attempt FAILED and its wake geometry has not yet been
inspected.** The VTK is expendable once the CT history is out.

Open items for Batch 1 to settle, to be judged from the logs rather than
assumed:

- **Wall time.** No 36,752-panel unsteady case has been run before, so
  24 h for 720 steps (NT=36, ~20 revs) is an estimate. If a case times out,
  resume with `RESTART_STEP`/`RESTART_NAME`/`RESTART_PATH` rather than a cold
  rerun.
- **`MAGVINF_PEAK=5.0`.** The staged-startup recipe was validated at 6000 RPM,
  where 5.0 m/s ≈ 1.1 × the hover induced velocity. At 5400 RPM `v_h` falls
  ~10%, so the pulse is slightly stronger relative to self-induction. Kept at
  the validated value rather than silently rescaled; if the startup misbehaves,
  4.5 m/s is the RPM-matched value.
- **`RELAX_RLXF` at NT=72.** Stock corrected-Pedrizzetti `rlxf` is 0.3
  (`FLOWVPM.jl:170`). The launcher halves it to 0.15 at NT=72 to hold the
  physical relaxation rate fixed as `dt` halves.

### 2026-07-27 23:40 — Batch 3 launched (studies A and B); configs verified at the header

Four jobs, all `velocity`, mesh `45_185_ct4`, 5400 RPM, NT=36, depth 4R,
rlxf 0.3, overlap 3, pps 2, merge_r 0.02, `SETTLE_REVS=12` (719 steps — matched
to the completed baselines), 128 GB / 24 h. Verified from each job's own banner
line in `logs/slurm/slurm-fp-p2e-*-129185*.out`, not from the submit command:

| Job | Case | Study | `relax_filter` | `das_eta` |
|---|---|---|---|---|
| 12918580 | `p2e_filt0p5` | B | 0.5R | 0.2 |
| 12918581 | `p2e_filt1p0` | B | 1.0R | 0.2 |
| 12918582 | `p2e_das1p0` | A | off | 1.0 |
| 12918583 | `p2e_das4p0` | A | off | 4.0 |

Study A's η=0.2 rung is the completed velocity baseline 12894164
(**CT = 0.06148 ± 0.00152**), which is why no `p2e_das0p2` job was submitted.

**Verified: the `Das` sweep is σ-clean (measured, was previously an assumption).**
BRAINSTORM/006 asserts that with `nwakerows=1` particle σ is independent of the
`Das` offset. That is not obvious for a *rotating* blade: the trailing filament
whose length sets σ (`sigma = dist*OVERLAP/P_PER_STEP`,
`src/FLOWPanel_wake.jl:621`) runs between consecutive shed rows
(`:1251`), and a node displaced tangentially by `Das` sits at a larger radius
(`sqrt(r²+Das²)`) and so sweeps a longer arc per step. At η=4, `Das ≈ 0.70R`, a
naive estimate of that effect is +35% on σ at the 75% station — which would have
confounded the whole sweep.

Measured instead, from step 9 of the two live jobs (particle `sigma` decoded
from the `.vtp`; helper `vtp_sigma.py` in the session scratchpad):

| case | n | σ/R min | σ/R median | σ/R max | median σ/R, 0.6–0.95R |
|---|---:|---:|---:|---:|---:|
| `das1p0` (η=1.0) | 1050 | 0.0174 | 0.0685 | 0.1111 | 0.0847 |
| `das4p0` (η=4.0) | 1103 | 0.0174 | 0.0693 | 0.1152 | 0.0845 |

A 4× change in η moves the median σ by **+1.2%** and the outboard median by
**−0.2%** (min identical, max +3.7%) — i.e. the offset does cancel between
consecutive steps, the curvature residual is negligible, and **η moves one thing
only**. Study A is a clean single-variable sweep as staged.

Also running: the Batch-2 σ-stock filtered retries at 256 GB, 12915034
(`f2_sigC_das0p2`, step 349) and 12915035 (`f2_sigC_das4p0`, step 268), 3 h 06 in.
Their 64 GB predecessors died at steps 421 and 420, so **step ~421 is the
decision point** — passing it is the first evidence that the 2.0R filter can
survive at production mesh size given headroom, and failing it again confirms
divergence for the σ-stock branch as it already did for σ-fine
(12915036/12915037 both re-failed at 384 GB).

All six study slots are occupied, so staged studies C (unfiltered refined-σ) and
E (`shed_with_induced_velocity`) remain queued behind the first free slot.

### 2026-07-28 01:45 — Study A has already flattened: CT is η-converged by η=1.0, and η=0.2 was under-converged by ~16%

Matched per-revolution CT (revs 0–4, thrust convention) computed directly from
each case's `monitors/*_monitor02_force_system1.csv`, 36 steps/rev. All velocity
formulation, production mesh, NT=36, 4R, stock σ. The η=0.2 column is the
**completed unfiltered baseline** 12894164, so this is a like-for-like ladder.

| rev | η=0.2 (12894164) | η=1.0 (12918582) | η=4.0 (12918583) | 0.2→1.0 | 1.0→4.0 |
|---:|---:|---:|---:|---:|---:|
| 0 | 0.06336 | 0.07886 | 0.08151 | +24.5% | +3.4% |
| 1 | 0.04698 | 0.05534 | 0.05569 | +17.8% | +0.6% |
| 2 | 0.04446 | 0.05119 | 0.05102 | +15.1% | −0.3% |
| 3 | 0.04863 | 0.05663 | 0.05661 | +16.5% | −0.04% |
| 4 | 0.06025 | 0.07034 | 0.07054 | +16.8% | +0.3% |

**Two findings.**

1. **The ladder is flat between η=1 and η=4**: a 4× change in the first-wake-row
   offset moves CT by −0.3% … +0.6% from rev 2 onward, against +15–18% for the
   0.2→1.0 step. That is study A's acceptance criterion ("CT-vs-η flattens")
   met, with **η = 1.0 the converged value** and η=0.2 — the driver's
   long-standing hard-coded default, used in *every* prior 006 number — sitting
   ~16% low. η=0.2 is the outlier, not η=4.
2. **The elevation is not the tangent-vector geometric artifact.** The
   23:40 entry showed the legacy straight-vector `Das` lands the first row at
   radius `r√(1+θ²)`: 0.015R too far out at η=1.0 but **0.22R at η=4.0**. If
   that radial displacement were driving CT up, η=4 would sit well above η=1.
   It does not — they agree to <0.5%. The +16% is therefore a real `Das`-length
   effect, and the arc correction (below), while still the right construction,
   does not explain or threaten this result.

**Why this matters for BRAINSTORM/006.** The settled unfiltered baseline is
CT = 0.06148 ± 0.00152. If the ~16% transient elevation carries through to the
settled hover value, the η-converged CT lands near **0.071** — inside 006's
0.068–0.072 acceptance band, and reached *without* the near-disk relaxation
filter that destabilises this mesh. That would make first-wake-row placement,
not wake resolution, the missing magnitude lever.

**Do not quote 0.071 as a result.** Revs 0–4 are inside the staged freestream
ramp/withdrawal (`SETTLE_REVS=12`), where all cases are still transient. The
deliverable is the 10-rev cycle-mean at revs 10–19, due when 12918582/12918583
complete (~10 h from 23:40). The 16% is a *lead*, measured in the same transient
window that the 2026-07-27 staging note flagged as "12–25% above".

**Held rungs η ∈ {0.5, 2.0} are now the informative ones.** The staging note said
"if 1.0 and 4.0 agree, `Das` is converged by η=1.0 and the intermediate rungs are
unnecessary." They agree — so η=2.0 is indeed unnecessary, but **η=0.5 is now the
valuable fill**: the entire jump happens between 0.2 and 1.0, and η=0.5 locates
its knee. Recommend η=0.5 for the next free slot ahead of η=2.0.

### 2026-07-28 01:45 / 02:34 — Batch 2 retry outcome: ALL FOUR filtered cases diverge, at any memory

| Job | Case | Mem | Outcome |
|---|---|---|---|
| 12915034 | `f2_sigC_das0p2` | 256 GB | **FAILED** — diverged, `max|CT| = 285.7` |
| 12915035 | `f2_sigC_das4p0` | 256 GB | **FAILED** at 5:00:45 — diverged, `max|CT| = 188.8` |
| 12915036 | `f2_sigF_das0p2` | 384 GB | FAILED at 2:41:23 |
| 12915037 | `f2_sigF_das4p0` | 384 GB | FAILED at 1:41:10 — `max|CT| = 374.3` |

12915035 was the first filtered case to show divergence **directly in flushed
per-rev CT** rather than by inference: `r7 = −0.0915`, `r9 = −9.089`. No
buffering caveat applies — unambiguous blow-up.

**The 2026-07-27 "3× memory changed nothing ⇒ divergence, not a memory ceiling"
reading is CONFIRMED and general.** Every filtered case fails by divergence
regardless of headroom, σ, or `Das` η.

**RETRACTED (written 01:45, corrected 02:34):** an earlier version of this entry
recorded that `sigC_das0p2` had "passed its old death step" at 256 GB and
concluded that its original 64 GB OOM was "at least partly a genuine memory
ceiling", and that "filter survivability is η-dependent". **Both claims were
wrong.** The case was still running when I read it (step 441, past the 421 where
its 64 GB attempt died) and it diverged shortly afterwards. Extra memory bought
it ~20 extra steps, not survival.

**Methodological lesson, and it is the same one the 2026-07-27 entry already
recorded in a different guise:** a running job that has passed a previous failure
point is *not* evidence of survival — only a terminal state is. The earlier
lesson was "magnitude tests are only valid on flushed (completed) logs"; the
generalisation is **do not draw conclusions from any in-progress job's prefix**,
whether the metric is max|CT| or elapsed steps.

### 2026-07-28 — `Das` arc construction implemented (code change, unreleased)

Prompted by Ryan's question of whether `Das` should follow the blade's arc
rather than a fixed vector from the TE. It should; the legacy construction is
first-order only. Implemented behind a switch, default on.

- `src/FLOWPanel_frames.jl`: `_rigid_back_displacement` (exact rigid backward
  path `Rot(ω̂,−|ω|τ)·d − d − v·τ`) and `accumulate_Das_arc!` (frame-tree walk
  mirroring `kinematic_velocity!`).
- `src/FLOWPanel_simulate.jl`: `_accumulate_Das_kinematic!(...; arc)`, preserving
  the legacy min-displacement-floor semantics; threaded through
  `initialize_Das!` and `simulate!`; recorded in run metadata.
- `src/FLOWPanel_warmstart.jl`: same path, so restarts rebuild `Das` identically.
- `examples/rotor_hover_pressure_comparison.jl`: env `DAS_KINEMATIC_ARC`
  (default `true`), written to `<run>_case_metadata.toml` as `das_kinematic_arc`.

Tests: new `kinematic Das arc construction` testset (25 assertions) in
`test/runtests_unit_simulate.jl` — radius preserved and azimuth exactly `−θ` up
to θ = 0.7 rad, pure translation bit-identical, `arc=false` reproducing the
legacy path exactly, arc-vs-tangent relative difference measured at exactly
`θ/2`. Regressions green: warm-start 14/14, replay 96/96, lifting body 42/42,
wake 82/82. `runtests_unit_liftingbody.jl`'s "initialize_Das! matches simulate
inline logic" was updated — it pinned the helper to the tangent construction, so
it now tests the legacy path explicitly plus the new default. The pre-existing
FLOWVPM SFS error in `runtests_unit_simulate.jl` reproduces on HEAD.

**Not deployed to the cluster** — every running job predates it, and at the
production η=0.2 the correction is 6e-4 R, so no existing 006 number changes.
Per Ryan, the in-flight straight-vector ladder runs to completion as the A/B
reference against a later arc-based ladder.

**Incident: local baseline data partly overwritten.** A driver smoke test was
launched with `RHPC_RUN_NAME=...`, but the driver reads **`RUN_NAME`**
(`examples/rotor_hover_pressure_comparison.jl:22`). The guessed variable was
ignored, so the run wrote into `data/rotor_hover_pressure_comparison/` and was
killed after 42 steps. Overwritten: both `monitors/*.csv` (360-step originals →
42 smoke steps), `metadata.toml`, all four `.pvd` manifests, and VTK steps 0–41.
**Intact:** `rotor_hover_pressure_comparison_CT_vs_rev.csv` (the settled ≈0.0624
history), VTK steps 42–359, all `particle_*.csv` analysis outputs,
`spanwise_loading_replay/`. Recoverable from Dropbox version history; no other
directory touched. Correct smoke invocation is `RUN_NAME=<name>`.

### 2026-07-28 11:26 — η=0.5 rung submitted; local source now DIVERGES from cluster (intentional)

| Job | Case | Study | Filter | `Das` η | State |
|---|---|---|---|---|---|
| 12920967 | `p2e_das0p5` | A — `Das` convergence | off | 0.5 | RUNNING (m12-1-7, 128 GB) |

Banner-verified: `relax_filter:offR das_eta:0.5 overlap:3 pps:2 settle:12`,
velocity, NT=36, 4R — matching the rest of the ladder. Chosen over the staged
η=2.0 because the ladder is already flat for η≥1: the entire +16% lives between
0.2 and 1.0, so **η=0.5 locates the knee** and makes study A a 4-point
convergence study (0.2 / 0.5 / 1.0 / 4.0) rather than a two-point inference.

**Cluster source is deliberately NOT updated with the arc `Das` change**, so the
whole ladder shares one construction. Record the divergence so nobody
misdiagnoses it (this study has been bitten before by cluster-vs-local drift):

| file | local (arc) | cluster (legacy tangent) |
|---|---|---|
| `examples/rotor_hover_pressure_comparison.jl` | `ff11fe2211ea04fabe857d73f80a7ec3` | `437288127018ad6c518d1e7ba55bd5c6` |
| `src/FLOWPanel_simulate.jl` | `6f35bb10d5029b2895aff1e6330105a0` | `46f5042323ec22565c59b877e6a670b6` |
| `src/FLOWPanel_frames.jl` | `4e5c1e6dd3a86130b36d5f60d918a3fa` | `2b29bc6402751f745b335f3190089962` |

Deploy the arc change only when starting a *new* ladder; re-verify md5s then.
Every Batch-3 job (12918580–83) and 12920967 ran on the cluster hashes above.

### 2026-07-28 11:46 — STUDIES A AND B DELIVER: η=1.0 gives CT = 0.0713, inside the 006 band; shallow filters survive

Three Batch-3 jobs reached terminal state, all **`COMPLETED` exit 0:0** (sacct):
12918580 `filt0p5` (9:05), 12918581 `filt1p0` (10:21), 12918582 `das1p0` (9:43).

**Ignore their `P2E_DIVERGED reason=nan_inf` flags** — the known-good completed
baselines 12894164/12894165 carry the identical flag. `nan_inf` fires on the
all-NaN `CT_kj` column and is a false positive for completed runs; only
`reason=magnitude` has discriminated real blow-ups in this study.

10-rev cycle-means, revs 10–19, computed from each case's
`monitors/*_monitor02_force_system1.csv` (36 steps/rev; mean of per-rev means,
std across the 10 rev-means). Production mesh, 5400 RPM, velocity, NT=36,
truncation 4R:

| case | `Das` η | filter | CT | std | rel | per-rev p-p |
|---|---:|---:|---:|---:|---:|---:|
| `p2e_vel_nt36_d4` (baseline) | 0.2 | off | 0.06148 | 0.00152 | 2.47% | 7.79% |
| **`p2e_das1p0`** | **1.0** | **off** | **0.07133** | **0.00159** | **2.23%** | **7.15%** |
| `p2e_filt0p5` | 0.2 | 0.5R | 0.06616 | 0.00109 | 1.64% | 4.34% |
| `p2e_filt1p0` | 0.2 | 1.0R | 0.06578 | 0.00170 | 2.59% | 7.96% |

Method check: the baseline recomputes to **0.06148 ± 0.00152**, identical to the
value logged on 2026-07-27, so the extraction is validated against a known result.

Per-rev series, revs 10→19:
- `das1p0`: 0.07372 0.07278 0.07141 0.06975 0.07176 0.07255 0.07070 0.06862 0.06983 0.07217
- `filt0p5`: 0.06536 0.06717 0.06682 0.06588 0.06502 0.06492 0.06779 0.06744 0.06615 0.06500
- `filt1p0`: 0.06476 0.06708 0.06689 0.06678 0.06542 0.06355 0.06249 0.06629 0.06773 0.06686

**Study A — ANSWERED.** `Das` η is the missing magnitude lever. η=1.0 gives
**CT = 0.0713 ± 0.0016**, **+16.0%** over the η=0.2 baseline and **inside 006's
0.068–0.072 acceptance band**, at 4R truncation, **unfiltered**, on the
production mesh. The settled result matches the +15–18% predicted from the revs
0–4 transient, so the transient reading held up. η=4.0 (12918583, still running,
rev 16) is tracking 0.0703–0.0733 — consistent, confirming the ladder is flat for
η ≥ 1 in the settled state too, not just in the transient.

**Study B — ANSWERED, and positively.** Its acceptance question was "does any
filter depth both complete 719 steps and raise the 10-rev cycle-mean above
0.0615?" **Yes, both.** 0.5R → 0.06616, 1.0R → 0.06578, both `COMPLETED`. The
2.0R filter's instability is therefore a *depth* effect, not a filter effect:
shallow filters are survivable at production mesh size where 2.0R is not.

**006 acceptance: the magnitude half is met, the stability half is NOT.** η=1.0
lands in band but its per-rev spread is 2.23% and p-p 7.15%, against tolerances
of 0.5% and 2%. This is a limit cycle at the right level, i.e. exactly the
"near-magnitude-but-oscillating" partial result 006's acceptance anticipates.

**The named next lever is the combination.** `filt0p5` has the **lowest ripple of
any case run at this mesh** (1.64% / 4.34%, vs 2.47% / 7.79% unfiltered) *and*
raises CT — the two levers are independent and both point the right way. The
obvious next configuration is **η=1.0 + `RELAX_FILTER_DOWNSTREAM_R=0.5`**, which
is the first candidate to satisfy both halves of 006 simultaneously. Caveat to
carry: the filter added +0.0047 at η=0.2, so if the effects are additive this
lands near 0.076 — above the band. That would itself be informative (it would
make the filter a magnitude knob to back off rather than a free stabiliser).

Caveat on construction: all of these ran the **legacy tangent `Das`** (cluster
hashes in the 11:26 entry). At η=1.0 the arc correction is 0.015R, and η=1.0 vs
η=4.0 agreement across a 0.015R→0.22R placement change argues the result is
insensitive to it — but the winning configuration should eventually be
re-confirmed under `DAS_KINEMATIC_ARC=true`.

### 2026-07-28 12:00 — Arc `Das` deployed to the cluster; arc A/B and the A+B combination submitted

**First cluster deployment of the arc `Das` construction.** Until now every
Phase 2e number, including the winning η=1.0 result, ran the legacy tangent
construction; the arc code existed only locally.

Deployed and md5-verified (local == cluster, all five):

| file | md5 |
|---|---|
| `src/FLOWPanel_frames.jl` | `4e5c1e6dd3a86130b36d5f60d918a3fa` |
| `src/FLOWPanel_simulate.jl` | `6f35bb10d5029b2895aff1e6330105a0` |
| `src/FLOWPanel_warmstart.jl` | `2cadde7f75ef08d84bd3a40f5516747c` |
| `examples/rotor_hover_pressure_comparison.jl` | `ff11fe2211ea04fabe857d73f80a7ec3` |
| `examples/run_dji9443_hover_ct_hpc.slurm.sh` | `f939cda70229f56a496123765bff7519` |

**Reproducibility guard.** The driver defaults `DAS_KINEMATIC_ARC=true`, so
deploying it would have silently flipped the construction for every subsequent
submission and contaminated the η ladder. The launcher therefore now exports
`DAS_KINEMATIC_ARC=false` globally (line 81, before the case block); only cases
that deliberately test the arc set it back to true. The banner echoes
`das_arc:<value>` so every future log states which construction it ran.

Two new case tags (**new tags on purpose** — the launcher does
`rm -rf data/$RUN_NAME`, so re-running `p2e_das1p0` would have destroyed the
0.0713 result, the same mechanism that wiped three histories on 2026-07-27):

| Job | Case | η | filter | `das_arc` | Purpose |
|---|---|---:|---:|---|---|
| 12921071 | `p2e_das1p0_arc` | 1.0 | off | **true** | A/B confirmation of the winning config under the correct construction |
| 12921072 | `p2e_das1p0_filt0p5` | 1.0 | 0.5R | false | A+B combination — first candidate to meet **both** halves of 006 |

Banner-verified on both. The combination deliberately runs **tangent** so it is
directly comparable to both parents (`p2e_das1p0` 0.07133 and `p2e_filt0p5`
0.06616), keeping one variable at a time.

Running jobs were unaffected by the deployment (Julia loads source at job start):
12918583 `das4p0` (12:23 elapsed) and 12920967 `das0p5` (37 min) both launched
against the pre-deploy source and are legacy-tangent ladder points, as intended.

Four jobs active, two slots free.

**What each outcome would mean.**
- *Arc ≈ 0.0713*: the construction correction is immaterial at η=1.0 and the
  headline result stands as reported. Expected, given η=1.0 and η=4.0 agree to
  <0.5% despite a 0.015R vs 0.22R placement difference.
- *Arc materially different*: the tangent ladder measured placement error rather
  than offset length, and the whole ladder needs re-running on the arc
  construction. This is the case that would retract the +16% attribution.
- *Combination in band with reduced ripple*: 006 is satisfied on both halves.
- *Combination near 0.076*: the levers are additive, the filter becomes a
  magnitude knob to back off rather than a free stabiliser, and the target
  configuration is η=1.0 with a filter shallower than 0.5R.

### 2026-07-28 12:36 — `das4p0` completes: settled ladder confirmed flat for η ≥ 1

12918583 `p2e_das4p0` **COMPLETED** (12:48:06, 20 revs). Settled study-A ladder,
10-rev cycle-means (revs 10–19), production mesh, 5400 RPM, velocity, NT=36, 4R,
unfiltered, legacy tangent `Das`:

| η | CT | std | rel | p-p | in 0.068–0.072? |
|---:|---:|---:|---:|---:|---|
| 0.2 | 0.06148 | 0.00152 | 2.47% | 7.79% | no (−11%) |
| 0.5 | *running* (12920967) | | | | |
| 1.0 | 0.07133 | 0.00159 | 2.23% | 7.15% | **yes** |
| 4.0 | 0.07190 | 0.00101 | 1.40% | 4.19% | **yes** (top edge) |

**η=1.0 → 4.0 is +0.80%, below either run's own cycle scatter** (2.23% / 1.40%).
The ladder is flat for η ≥ 1 in the *settled* state, not merely in the transient
where it was first spotted. Study A's acceptance — "CT-vs-η flattens, naming a
converged η and an uncertainty band" — is met: **converged η = 1.0**, converged
**CT ≈ 0.0713–0.0719**, uncertainty band ±0.0016 (the cycle scatter, which
dominates the η sensitivity above η=1).

`das0p5` interim (rev 5, still inside the freestream schedule), matched-rev
comparison at rev 4: η=0.2 → 0.06025, **η=0.5 → 0.06805**, η=1.0 → 0.07034,
η=4.0 → 0.07054. η=0.5 captures **~77%** of the 0.2→1.0 jump, so the knee sits
**below η=0.5** and the curve is already near-saturated there. Transient values —
the settled number governs.

**Observation, not a result:** η=4.0 has the lowest ripple of any unfiltered case
at this mesh (1.40% / 4.19% vs 2.47% / 7.79% at η=0.2), hinting that a longer
first-row offset also smooths the limit cycle. One run, and η=4 is the
geometrically distorted end under the tangent construction, so this is a thing to
watch in the arc A/B (12921071), not a claim.

### 2026-07-28 13:00 — What `Das` actually is, and the principled argument for η=1 (corrects a mechanism claim)

Ryan challenged an explanation I had given verbally — that raising η weakens the
self-induced downwash by moving the near wake away from the blade — on the
grounds that the body is solved with a doublet wake on `Das` regardless. **He is
right; that explanation was wrong.** Verified in the source:

- `src/FLOWPanel_elements_fmm.jl:1044-1066`: the body's shedding element is a
  doublet panel spanning `TE1 → TE2 → TE2+Db → TE1+Da` (two triangles,
  `vw = v + Da`). It is **attached at the trailing edge**. There is no void
  between the TE and TE+`Das` at any η.
- `src/FLOWPanel_simulate.jl:1021` (`update_TE!`):
  `nodes[:, 1, i_shed] .= v1 .+ view(Das, :, i_shed)` — the `PanelWake`'s **row 1
  sits at TE + `Das`**, i.e. exactly where the attached panel ends.
- `src/FLOWPanel_wake.jl` contains **no** reference to `Das`; row 2 is one
  convection step beyond row 1, which is why particle σ measured
  `Das`-independent (23:40 entry).

Resulting chain: TE →[attached doublet panel, length `Das`]→ TE+`Das`
→[`PanelWake` row, one step of travel]→ particles.

**Corrected mechanism.** Raising η does not remove near-wake influence; it moves
the **boundary between two representations**, lengthening the flat attached
doublet sheet and pushing the start of the particle wake downstream by `Das`.

**Principled argument for η=1 (replaces "it's where the curve flattens").** The
particle wake advances one step per timestep, so the stretch between the TE and
where particles legitimately begin is exactly one step of TE travel — i.e.
**η=1 by construction**. The attached panel should cover that gap and no more.
- η<1 seeds particles *inside* the stretch the attached panel already spans, so
  they sit closer to the blade than the shedding process places them and
  over-induce downwash. **η=0.2 was therefore an inconsistency, not a
  conservative choice**, and the +16% recovers from it rather than tuning toward
  a target.
- η>1 hands wake to a flat panel that should have rolled up as particles.

This matters because it is dt-relative by construction (`Das = η·dt·|V_te|`): at
η=1 the flat-panel region shrinks automatically with dt. **The consistency
argument predicts η=1 is correct at any dt**, which is exactly what 12921559
(NT=72, η=1.0) tests. Scale reference: at NT=36 one step of TE travel at 0.75R
is **1.02 local chords** (2.0 chords at the tip, 10° azimuth), so the flat-panel
region is currently chord-scale — ample motivation for the refinement.

**Residual caveat.** Flatness over η ∈ [1,4] shows CT is insensitive to *where*
the panel/particle handoff sits across a 4× range. That is reassuring but is
insensitivity to the partition, not proof that either representation is correct.
The dt refinement is the test that can distinguish them.

`Das` in local chords (NT=36, RPM 5400, R=0.119 m, chord from
`examples/rotor_hover_scan/processed/dji9443_brainstorm_item003.csv`):

| η | azimuth | r/R=0.25 | 0.50 | 0.75 | 0.90 | tip |
|---:|---:|---:|---:|---:|---:|---:|
| 0.2 | 2° | 0.03 | 0.09 | 0.20 | 0.31 | 0.40 |
| 0.5 | 5° | 0.08 | 0.24 | 0.51 | 0.77 | 1.00 |
| 1.0 | 10° | 0.17 | 0.47 | 1.02 | 1.55 | 2.00 |
| 4.0 | 40° | 0.66 | 1.90 | 4.08 | 6.18 | 7.98 |

### 2026-07-28 12:55 — NT=72 at the converged η submitted (study D, reframed)

| Job | Case | NT | η | rlxf | filter | Wall |
|---|---|---:|---:|---:|---|---|
| 12921559 | `p2e_nt72_das1p0` | 72 | 1.0 | 0.15 | off | 48 h |

Banner-verified (`NT:72 rlxf:0.15 das_eta:1.0 das_arc:false`). Launcher md5
`c6e5ba83dbe37269bb4f1b09c8dc5a0d`, deployed and verified on the cluster.

**Reframes study D.** The staged plan proposed scaling η ∝ NT to hold `Das`
fixed in physical units. That isolates dt from `Das`, but it is *not* the test
of the scheme we actually run: the scheme's specification is "shed one timestep
back", i.e. **η fixed at 1.0** with `Das` shrinking as dt shrinks. That is what
this job does.
- CT stays ≈0.0713 ⇒ η=1 is genuinely dt-converged and the headline result holds.
- CT falls toward 0.0615 ⇒ the plateau was an artifact of chord-scale dt and the
  +16% attribution must be retracted.

### 2026-07-28 13:30 — MEASURED: `Das` is set at 40% RPM; σ law verified; the core-overlap hypothesis is REFUTED

Measurement on `data/p2e_das1p0` settled step 700 (41,778 particles), σ decoded
from the `.vtp`. Script: `das_vs_sigma.py` (session scratchpad).

**Correction — `Das` is 2.5× smaller than the 13:00 chord table states.**
`initialize_Das!` runs once at `t=0`, where
`spinup_fraction(0) = SPINUP_START_FRACTION = 0.4` (launcher line 61), and
`rotate_Das!` only *rotates* the stored vector — it never rescales. So `Das`
keeps its 40%-RPM magnitude for the whole run while the rotor spins to full RPM.
Expressed in steps of TE travel **at operating RPM**, `η_eff = 0.4·η`:

| nominal η | η_eff | `Das` at 0.75R | in local chords |
|---:|---:|---:|---:|
| 0.2 | 0.08 | 1.25 mm | 0.08 |
| 0.5 | 0.20 | 3.12 mm | 0.20 |
| 1.0 | 0.40 | 6.23 mm | 0.41 |
| 4.0 | 1.60 | 24.9 mm | 1.63 |

**Supersedes the 13:00 table** (which used full RPM). Note the CT plateau
η ∈ [1,4] spans η_eff ∈ [0.4, 1.6], bracketing one step of TE travel — an
observation, not a selection criterion, since both ends give the same CT.

**σ law verified against data.** Predicted trailing-particle
`σ/r = Δψ·OVERLAP/PPS` vs measured medians:

| r/R band | measured σ/R (p50) | predicted σ/R |
|---|---:|---:|
| 0.30–0.50 | 0.1179 | 0.1047 |
| 0.50–0.70 | 0.1758 | 0.1571 |
| 0.70–0.85 | 0.2175 | 0.2029 |
| 0.85–1.00 | 0.2264 | 0.2422 |

Within ~10%; overall median consistent with the 0.145 logged for this mesh.

**HYPOTHESIS REFUTED.** `Das/σ = 0.4·η/(OVERLAP/PPS) = 0.267·η`, exactly
radius-independent (both scale as `r`): 0.053 (η=0.2) → 0.133 (0.5) → 0.267
(1.0) → 1.067 (4.0). If "the newly shed particle's core straddles the trailing
edge" were the controlling pathology, CT should improve until `Das/σ ≳ 1`
(η≈4). CT is flat from η=1, where the core still exceeds the offset **3.7×**.
The core-overlap criterion therefore does **not** explain the η plateau and is
withdrawn as an explanation. Together with the two retractions above (the
"near-wake pushed out of influence" mechanism and the "particles seeded inside
the panel span" inconsistency), **no proposed mechanism for the +16% survives
scrutiny; the effect remains empirical.**

**Real finding, independent of the η question.** The shed particle core is
**~1.5 local chords** at 0.75R (σ = 23.4 mm vs chord 15.3 mm) and exceeds the TE
offset at every η run so far. Newly shed particles engulf the blade section that
produced them. This is `σ = 1.5 ×` one step of TE travel by construction
(`OVERLAP/PPS = 3/2`), so it is inherited directly from the timestep. This is a
sharper statement of near-wake under-resolution than the σ/R = 0.145 figure and
belongs to **item 012**.

**Caveat now attached to 12921559 (NT=72).** At fixed η, halving dt halves
`Das` *and* σ along with dt, so the run refines three things at once (as the
staged note warned). `Das/σ = 0.267` is invariant under it. If CT holds at
0.0713 that is a strong simultaneous-convergence result; if it moves, the run
alone cannot attribute the change to dt vs `Das` vs σ.

### 2026-07-28 (later) — BRAINSTORM 014 attribution cases STAGED (not yet submitted): floor test, nwakerows handoff trades, un-frozen `Das`

Four new case tags added to `examples/run_dji9443_hover_ct_hpc.slurm.sh`, all on
the unfiltered carrier (velocity, NT=36, 4R, `SETTLE_REVS=12`, 719 steps,
tangent `Das`) for like-for-like comparison against the completed ladder
(η=0.2 → 0.06148, η=1.0 → 0.07133, η=4.0 → 0.07190):

| tag | config delta | question |
|---|---|---|
| `p2e_das0p2_nofloor` | η=0.2, `DAS_MIN_DISPLACEMENT_R=1e-6` | **Highest priority.** Was the "+16% η effect" partly the 0.01R floor (which sets `Das` over ~68% of the span at η=0.2)? Pre-registered reads in BRAINSTORM/014. |
| `p2e_nrows2_dassmall` | η=0.1, floor off, `NWAKEROWS=2` | Tiny rigid sheet + one free row ⇒ handoff ≈ 1.0 step (between the flat η=1.0/η=4.0 points). CT ≈ 0.0713 ⇒ only total handoff matters; CT → 0.0615 ⇒ the rigid `Das` row itself carries the effect. |
| `p2e_nrows4_das1p0` | η=1.0, `NWAKEROWS=4` | Extent insensitivity at fixed rigid sheet (handoff ≈ 3.4 steps). |
| `p2e_das1p0_refresh` | η=1.0, `DAS_REFRESH=true` | `Das` re-derived from *current* kinematics each step ⇒ the frozen 0.4× spin-up factor disappears (η_eff = 1.0, 2.5× the frozen length). Does 0.0713 survive making η physically meaningful? |

Driver support added the same day (`examples/rotor_hover_pressure_comparison.jl`):
`NWAKEROWS` env (was hard-coded 1; the log line and metadata literal that would
have silently lied are now interpolated) and `DAS_REFRESH` env (new
`set_Das_refresh` kwarg in `simulate!`; when true the driver skips its own
`initialize_Das!` so simulate!'s per-step re-derivation is the only writer —
pre-initializing would double-accumulate). Unit-tested
(`test/runtests_unit_simulate.jl`, "set_Das_refresh" testset). NOT yet deployed
to the cluster: cluster src must be synced (md5 check per the 2026-07-25
lesson) before submitting `p2e_nrows*` or `p2e_das1p0_refresh`;
`p2e_das0p2_nofloor` needs only the launcher (env hook already existed in the
deployed driver).

### 2026-07-28 (later still) — 014 discriminators SUBMITTED: `p2e_das0p2_nofloor` (12925945), `p2e_nrows2_dassmall` (12925947)

Cluster synced first (md5 check): only the three files edited today differed
(`src/FLOWPanel_simulate.jl` — verified by diff to be exactly the
`set_Das_refresh` hunks — plus the driver and this launcher); pushed all three.
Both jobs at `--mem=24G` (velocity-only). Six active jobs = cap;
`p2e_nrows4_das1p0` and `p2e_das1p0_refresh` queue for free slots.

Pre-registered reads (from BRAINSTORM/014 + the 2026-07-28 oracle results):
- `p2e_das0p2_nofloor`: CT→0.0713 ⇒ floor artifact; ~0.0615 ⇒ Das itself; between ⇒ both.
- `p2e_nrows2_dassmall`: CT≈0.0713 ⇒ only total handoff distance matters (rigid
  sheet has no special status); CT→0.0615 ⇒ the rigid `Das` row carries the effect.
- In-flight 12921559 (NT=72, η=1.0): oracle log-law-only prediction 0.0708–0.0712;
  a fall toward 0.0615 locates the rotor's amplified regime near `Das`≈0.2 chords.

### 2026-07-28 21:20 — `das0p5` completes; floor test, arc A/B, combination and `nwakerows=2` all report

Ryan submitted two BRAINSTORM/014 experiments directly (12925945, 12925947) and
added banner fields `das_min_R` / `das_refresh` / `nwakerows`.

**`p2e_das0p5` COMPLETED** (09:42:18). 10-rev cycle-mean revs 10–19:
**CT = 0.06942 ± 0.00185** (rel 2.66%, p-p 7.30%). Study A ladder, settled:

| η | CT | fraction of the 0.2→1.0 jump |
|---:|---:|---:|
| 0.2 | 0.06148 | — |
| 0.5 | **0.06942** | 81% |
| 1.0 | 0.07133 | 100% |
| 4.0 | 0.07190 | flat (+0.8%) |

Knee below η=0.5, confirming the 77% transient estimate.

Everything below is **matched revs 5–10 and still transient** (`SETTLE_REVS=12`);
settled cycle-means follow on completion.

**FLOOR TEST ANSWERED — and it strengthens the η attribution.** 12925945
(`p2e_das0p2_nofloor`: η=0.2, `das_min_R=1e-6`, nwakerows 1):

| rev | floor ON (0.01R) | floor OFF |
|---:|---:|---:|
| 5 | 0.06952 | 0.05786 |
| 8 | 0.06093 | 0.05149 |
| 10 | 0.06284 | 0.05411 |

Removing the floor **lowers** CT ~14% — the *opposite* of the risk flagged when
the confound was found. The floor was propping the η=0.2 end **up**, so it
masked part of the η sensitivity rather than manufacturing it. Floor-free, the
ladder runs ~0.053 → 0.0713, i.e. **~+30%**, not +16%. **The confound is removed
from 014's critical path and η is unambiguously the lever** — while the
model-form uncertainty concern gets *worse*, not better.

**ARC A/B: immaterial, as predicted.** 12921071 vs 12921072's tangent parent:
rev 5 0.08037 vs 0.08073, rev 8 0.07020 vs 0.07005, rev 10 0.07342 vs 0.07372 —
within 0.5% throughout. **CT = 0.0713 stands under the corrected arc
construction**, consistent with the η=1 vs η=4 insensitivity to a 0.015R→0.22R
placement change.

**COMBINATION OVERSHOOTS, as forecast in the 12:00 entry.**
`p2e_das1p0_filt0p5` runs ~5% above η=1.0 alone (rev 10: 0.07769 vs 0.07372),
revs 10–15 at 0.0745–0.0784 — **above the 0.068–0.072 band**. The levers are
additive ⇒ the 0.5R filter is a **magnitude knob to back off**, not a free
stabiliser. 006's stability half will need a *weaker* filter at η=1.0, not this
one.

**`nwakerows=2` — a clean NEGATIVE, and the most informative result here.**
12925947 (`p2e_nrows2_dassmall`: η=0.1, floor off, **nwakerows 2**) gives
**0.0366 / 0.0350 / 0.0332** at revs 5 / 7 / 10 — roughly **half** the baseline.

Its *total* handoff distance is **larger** than η=1.0's: `Das`(η=0.1) ≈ 0.04
steps of TE travel plus one full free row ≈ 1.0 step, versus 0.4 steps for
(η=1.0, nrows=1). Yet CT is 0.033 vs 0.0713.

⇒ **CT does not depend on the total distance to the particle handoff.** A free
convecting row is **not** interchangeable with a longer rigid `Das`. This
refutes the "`Das` is merely a way to buy panel-wake extent" branch of 014's
Proposal 1, and points at **Proposal 2**: as `Das → 0` the attached doublet
panel degenerates and something essential to the circulation — most plausibly
the **Kutta enforcement** — degrades with it.

The floor test and the `nwakerows` test are mutually consistent: both say CT
collapses as `Das` shrinks, **regardless of what is downstream of it**. That is
the first coherent, non-refuted signal about what η is actually doing, and it
relocates 014's centre of gravity from "how much near wake is a flat sheet" to
"what does the attached panel do for the Kutta condition".

### 2026-07-28 21:15 — `p2e_das0p5` COMPLETE: 10-rev cycle-mean CT = 0.06942 ± 0.00185

Job 12920967, 719 steps, exit clean. Matched-window (revs 10–19) cycle-mean from
`p2e_das0p5_CT_per_rev.csv`: **0.06942 ± 0.00185**. The η ladder is now
0.2 → 0.06148, **0.5 → 0.06942**, 1.0 → 0.07133, 4.0 → 0.07190: 81% of the
0.2→1.0 jump is realized by η=0.5. In the physical-`Das` frame (frozen,
`Das`@0.75R ≈ 0.21 → 0.41 local chords), the 0.5→1.0 doubling is **+2.75%** —
an order of magnitude steeper than the oracle's quasi-steady log slope
(0.1–0.4%/doubling) — while 1.0→4.0 (+0.8% over two doublings) matches it. So
the rotor's amplified regime extends to `Das` ≈ 0.2–0.4 local chords and the
log plateau begins near `Das` ≈ 0.4 chords (η ≥ 1 frozen). Whether the
amplifier is the floor (19.7% of span clamped at η=0.5 vs 3.6% at η=1) or a
genuine near-field effect is exactly what 12925945 (nofloor) and 12925947
(nrows2) are running to decide.

### 2026-07-28 21:42 — combination completes: levers are cleanly ADDITIVE; it overshoots; no config meets the stability gate

`p2e_das1p0_filt0p5` **COMPLETED** (09:29:01). 10-rev cycle-mean revs 10–19:
**CT = 0.07625 ± 0.00140** (rel 1.84%, p-p 5.38%).

**Superposition confirmed to 0.3%.** η=1.0 alone 0.07133; the 0.5R filter
contributes +0.00468 at η=0.2 (0.06148→0.06616). Predicted additive result
0.07601 vs measured **0.07625**. The 12:00 entry's advance prediction ("near
0.076 — above the band") is confirmed, so **the 0.5R filter is a magnitude knob,
not a free stabiliser**, exactly as flagged.

Complete settled table, against 006's gates (rel ≤ 0.5%, p-p ≤ 2%):

| case | η | filter | CT | rel | p-p | in band | stable |
|---|---:|---:|---:|---:|---:|:--:|:--:|
| `p2e_vel_nt36_d4` | 0.2 | off | 0.06148 | 2.47% | 7.79% | no | no |
| `p2e_das0p5` | 0.5 | off | 0.06942 | 2.66% | 7.30% | no | no |
| `p2e_das1p0` | 1.0 | off | 0.07133 | 2.23% | 7.15% | **yes** | no |
| `p2e_das4p0` | 4.0 | off | 0.07190 | **1.40%** | **4.19%** | **yes** | no |
| `p2e_filt0p5` | 0.2 | 0.5R | 0.06616 | 1.64% | 4.34% | no | no |
| `p2e_das1p0_filt0p5` | 1.0 | 0.5R | 0.07625 | 1.84% | 5.38% | no (over) | no |

**NO configuration tested meets the stability gate.** Best is `das4p0` at
1.40% / 4.19% — still 2.8× and 2.1× the tolerances. The combination improved
ripple over η=1.0 alone but did not beat `filt0p5`, while pushing magnitude out
the top of the band.

**Discipline note.** Because the two levers superpose cleanly, a configuration
near η≈0.35 + 0.5R filter would land ~0.070. That is **tuning two knobs to a
known target** and must not be reported as a result — it is exactly what
BRAINSTORM/014 exists to guard against. Record it only as evidence of
additivity.

**Notable:** `das4p0` is in band **and** has the lowest ripple of any unfiltered
case, with no filter at all. If the stability half is ever met it may come from
the η side rather than from damping — worth remembering, though η=4 is the far
end of the ladder and 014's verdict governs whether any η is defensible.

### 2026-07-28 22:01 — arc A/B CLOSED (immaterial); Γ_TE diagnostic locates the η effect in the bound circulation

**`p2e_das1p0_arc` COMPLETED** (09:49:38). 10-rev cycle-mean revs 10–19:
**CT = 0.07108 ± 0.00184** (rel 2.59%, p-p 8.27%) vs the tangent parent
`p2e_das1p0` **0.07133 ± 0.00159**. Difference **−0.35%**, far inside either
run's ±2.2–2.6% cycle scatter.

**Arc-vs-tangent question CLOSED: immaterial.** The corrected arc construction
does not change the result at η=1.0, so the headline CT ≈ 0.0713 stands under
the correct construction, and only the offset's *length* matters, not its
*direction*. (The arc code remains the right construction on correctness
grounds; it is simply not a lever here.)

**Γ_TE diagnostic** (read-only, from `monitor03_bound_circulation`, no new
compute; settled steps ≥ 360). Span-averaged |Γ_TE| against CT:

| η | mean Γ_TE | Γ ratio | CT | CT ratio |
|---:|---:|---:|---:|---:|
| 0.2 | 0.18540 | 1.000 | 0.06148 | 1.000 |
| 0.5 | 0.21138 | 1.140 | 0.06942 | 1.129 |
| 1.0 | 0.22193 | 1.197 | 0.07133 | 1.160 |
| 4.0 | 0.22550 | 1.216 | 0.07190 | 1.169 |

Γ tracks CT in near-lockstep ⇒ **the η sensitivity acts through the blade's own
bound circulation**, not through the wake acting downstream on unchanged
loading. Spanwise, the effect is **concentrated inboard** (Γ(0.25R) +48% from
η=0.2→4, while Γ(0.75R) saturates by η=0.5 and slightly declines at η=4) — i.e.
where `Das/chord` is smallest and the attached panel is nearest degenerate.

**Explicitly NOT established:** this does not separate Kutta enforcement from
wake-induced downwash, since altered induced velocity would also change
effective incidence → Γ → CT.

**Next and cheapest discriminator: a STEADY rigid-wake `Das` sweep** (item 001's
machinery, seconds per case). No wake evolution ⇒ `Das` changes only the
attached panel geometry. Γ moving ~20% there ⇒ Kutta/geometry, and 014 becomes a
question answerable entirely in cheap steady solves; Γ insensitive ⇒ the route is
induced velocity. **Run before any further unsteady η runs.** Detail in
`BRAINSTORM/014`.

### 2026-07-28 23:40 — arc A/B and A+B combination COMPLETE; nrows4 + refresh submitted

- **`p2e_das1p0_arc` (12921071): 10-rev cycle-mean CT = 0.07108 ± 0.00184** vs
  tangent 0.07133 ± 0.00159 — Δ −0.35%, well inside cycle scatter. **The η=1.0
  magnitude result is robust to the arc-vs-tangent `Das` construction.**
- **`p2e_das1p0_filt0p5` (12921072): CT = 0.07625 ± 0.00140** — the additive
  caveat came true almost exactly (0.0713 + 0.0047 ≈ 0.0760): the 0.5R filter
  is a **magnitude knob (+0.005), not a free stabilizer**, and the combination
  OVERSHOOTS the 0.068–0.072 band. Ripple improvement is marginal (per-rev std
  1.8% vs 2.2%; within-rev p-p ~6.3% vs 7.2%). 006's stability half is NOT
  solved by this combination; treat relaxation filtering as part of the CT
  error budget, not a stacking lever.
- Submitted with the freed slots: `p2e_nrows4_das1p0` (12927923),
  `p2e_das1p0_refresh` (12927924).
- `p2e_nrows2_dassmall` (12925947) is slow but healthy (julia ~12 cores active
  on m12-2-21; log block-buffering makes step 402 stale) — expect it hours
  behind `p2e_das0p2_nofloor` (12925945, step 609/719 and on schedule).

### 2026-07-28 (Ryan ruling) — relaxation policy

Production configs: **stock relaxation, no near-disk filter** (rule). Filter and
`rlxf` adjustments are diagnostics/convergence probes only, revisited after
near-disk resolution (012). Detail in BRAINSTORM/006 "Ruling 2026-07-28".

### 2026-07-29 04:00 — floor-free and `nwakerows=2` settle; true η sensitivity is +37%

`p2e_das0p2_nofloor` **COMPLETED** (08:34:52) and `p2e_nrows2_dassmall`
**COMPLETED** (11:45:53). 10-rev cycle-means, revs 10–19:

| case | CT | std | rel | p-p |
|---|---:|---:|---:|---:|
| baseline η=0.2 (floor 0.01R) | 0.06148 | 0.00152 | 2.47% | 7.79% |
| **η=0.2, floor OFF (1e-6)** | **0.05215** | 0.00078 | 1.49% | 5.10% |
| **η=0.1, nrows=2, floor OFF** | **0.03362** | 0.00046 | 1.35% | 3.71% |

- **Floor-free η sensitivity = +36.8%** (0.05215 → 0.07133), not the +16%
  measured with the floor active. The floor inflated the η=0.2 end by 15.2% and
  masked over half the dependence. **All ±16% statements are superseded by ±37%.**
- **`nwakerows` negative confirmed settled:** larger total handoff distance
  (≈1.04 vs 0.40 steps) yet less than half the CT. `Das` length, not wake extent,
  is the controlling variable.
- Ripple: low-CT cases have lower *relative* ripple, but absolute std tracks CT
  magnitude, so most of that is proportional. **η=4.0 remains the only case with
  both high CT and low relative ripple.** Nothing meets 006's 0.5% / 2% gates.

**`p2e_nt72_das1p0` is healthy**, at step ~1001 of 1438 after 15 h (~6 h left).
Its `P2E_CT` rev count lags actual progress — judge it by particle-file count,
not the digest.

**Two new jobs launched by Ryan** (banners to verify when they start):
12927923 `fp-p2e-nrows4-das1p0` (nwakerows=4 at η=1.0, extends the decoupling
trend) and 12927924 `fp-p2e-das1p0-refresh` (un-freeze `Das`, addressing the
"frozen at t=0 40%-RPM" open question).

### 2026-07-29 — banners verified for Ryan's two new jobs, with predictions fixed in advance

| Job | Case | η | rows | `das_refresh` | floor | Started |
|---|---|---:|---:|---|---|---|
| 12927923 | `p2e_nrows4_das1p0` | 1.0 | **4** | false | 0.01R | 00:21 |
| 12927924 | `p2e_das1p0_refresh` | 1.0 | 1 | **true** | 0.01R | 01:29 |

(At η=1.0 the 0.01R floor is active over only ~3.6% of the shedding span, so
these are effectively floor-free.)

**`nrows4_das1p0` — the clean complement to `nrows2_dassmall`.** It holds `Das`
**fixed** at η=1.0 (0.40 steps) and adds free rows, giving a total handoff
distance of ≈3.40 steps versus 0.40 for `p2e_das1p0`.
- **CT ≈ 0.0713 (unchanged) ⇒ free rows are inert and `Das` length alone
  controls CT.** Combined with `nrows2` (small `Das` + extra row → CT collapsed),
  this would close Proposal 1 completely: the panel/particle handoff distance is
  irrelevant, only the rigid attached panel matters.
- **CT rises toward/past 0.0719 ⇒ wake extent does contribute**, and the
  `nrows2` collapse was driven by its tiny `Das` rather than by rows being inert.

**`das1p0_refresh` — tests the frozen-`Das` question.** Recomputing `Das` from the
*current* `|V_te|` each step removes the 0.4 spin-up factor, so at operating RPM
its `Das` ≈ **1.0 steps** of TE travel rather than the frozen 0.40 — i.e. it
should behave like frozen η≈2.5.
- **CT ≈ 0.0719 (the saturated plateau) ⇒ consistent**, and it confirms that the
  frozen-vs-refreshed distinction acts purely through the effective `Das` length,
  with `η_eff = 0.4·η` the right bookkeeping.
- **CT materially different from the frozen ladder at matched effective `Das`
  ⇒ refreshing does something beyond rescaling**, e.g. through the per-step
  variation of `Das` during spin-up, and the ladder cannot be reduced to a single
  effective length.

Both sit in the saturated region of the floor-free ladder, so neither is expected
to move CT much; their value is in **discriminating mechanisms**, not magnitude.

### 2026-07-29 — `das1p0_refresh` CANCELLED as INVALID: `Backslash` caches an LU built with `Das = 0`

12927924 cancelled at ~7.5 h / step 444. **The run was invalid from step 0** and
its numbers must not be used.

**Symptom.** CT ≈ 0.0076 at revs 10–11 (vs 0.0713 for the matched frozen case),
with bound circulation collapsed by **21.7×** (mean |Γ_TE| 0.010484 vs 0.227854
over steps 360–430). Γ collapsing with CT localises this to the **solve**, not
force integration.

**Cause.** `Backslash` factorizes once at construction
(`src/FLOWPanel_solver.jl:230`, `Glu = lu!(G)`), and `G` includes the attached
wake panel's influence, so **`G` depends on `Das`**. This assumption is stated
in `src/FLOWPanel_formulation.jl:357`: one-time assembled operators assume
"rigidly invariant body+`Das` geometry (the same assumption `Backslash` reuse
already makes)."

The driver **skips** `initialize_Das!` when `DAS_REFRESH=true`
(`examples/rotor_hover_pressure_comparison.jl:409`), and `Backslash(rotor)` is
constructed *after* that block (line 416). **So `G` was factorized with
`Das = 0`** — a zero-length attached panel, i.e. no Kutta panel at all — and
every subsequent solve reused that stale factorization regardless of the
per-step refresh.

The skip itself was well-motivated: the accumulate helpers are `+=`-based, so
pre-initialising *and* refreshing would double-accumulate. The per-step path
correctly calls `_das_zero!` first (`FLOWPanel_simulate.jl:865`). The defect is
only that skipping the pre-init leaves `Das = 0` for the one operation that
bakes it in permanently.

**What a valid version needs** (any of):
1. Retain the pre-init so `G` is factorized with a representative `Das`, and
   keep the per-step `_das_zero!` + refresh (accepting that `G` is then slightly
   stale whenever `Das` changes).
2. Rebuild/refactorize `G` when `Das` changes — correct but expensive.
3. Use a solver that does not cache a `Das`-dependent factorization.

Even with (1), `Backslash` reuse is only defensible if `Das` is near-constant
after spin-up; during spin-up `Das` varies with ω by construction.

**The 014 question "does un-freezing `Das` change the η sensitivity?" remains
OPEN and untested.** Nothing was learned about it.

**Generalisable lesson:** any run that varies `Das` (or body geometry) *during*
a simulation is incompatible with cached-operator solvers unless the operator is
rebuilt. Check solver caching before adding a per-step geometry knob.

### 2026-07-29 — `nrows4_das1p0` completes; free wake rows contribute +4.2% (vs `Das`'s +37%)

12927923 `p2e_nrows4_das1p0` **COMPLETED** (09:36:18). η=1.0, `nwakerows=4`,
unfiltered, floor 0.01R (inactive at η=1.0). 10-rev cycle-mean revs 10–19:

| case | η | rows | CT | std | rel | p-p |
|---|---:|---:|---:|---:|---:|---:|
| `p2e_das1p0` | 1.0 | 1 | 0.07133 | 0.00159 | 2.23% | 7.15% |
| **`p2e_nrows4_das1p0`** | 1.0 | **4** | **0.07431** | 0.00096 | **1.30%** | **4.52%** |

**+4.2% from three extra free panel rows at fixed `Das`.** This is a wake-*extent*
effect (a longer panel sheet), not a subdivision effect, and it should saturate as
rows → many. Contrast with `Das` itself: **+37%** across the floor-free ladder. So
`Das` dominates by roughly an order of magnitude, but free rows are **not inert** —
my pre-registered "CT ≈ 0.0713 unchanged" prediction (12:00 entry) is **wrong**;
the alternative branch ("wake extent does contribute") is what occurred.

Note 0.07431 sits **above** 006's 0.068–0.072 band, and this case has the lowest
relative ripple of any high-CT run at this mesh (1.30% / 4.52%).

Also completed/closed since the last entry: 12921071 `das1p0_arc` **0.07108**
(arc immaterial), 12921072 `das1p0_filt0p5` **0.07625** (overshoot, levers
additive), 12927924 `das1p0_refresh` **CANCELLED as invalid** (`Backslash`
factorized `G` with `Das = 0`).

### 2026-07-29 — 014 discriminators HARVESTED: the floor was propping CT UP; total-handoff hypothesis DEAD; panel-sheet extent has its own +4% lever

10-rev cycle-means (revs 10–19), production mesh, unfiltered, stock relaxation:

| case | job | CT | vs pre-registered reads |
|---|---|---|---|
| `p2e_das0p2_nofloor` (η=0.2, floor off) | 12925945 | **0.05215 ± 0.00078** | Neither pre-registered branch: CT went DOWN, not up. The 0.01R floor was *increasing* `Das` over 68% of the span (kinematic η=0.2 `Das` < 0.01R there), so the floored 0.06148 was floor-PROPPED. Consistent with monotone CT-vs-`Das`: the clean η=0.2 point is 0.0522, making the true η 0.2→1.0 sensitivity **+37%**, not +16%. |
| `p2e_nrows2_dassmall` (η=0.1, floor off, 2 rows) | 12925947 | **0.03362 ± 0.00046** | **"Only total handoff distance matters" is DEAD**: handoff ≈ 1.0 step (inside the flat η=1–4 zone) yet CT collapsed far below even the η=0.2-nofloor value. The rigid `Das` sheet itself (and/or a near-TE free-row trajectory pathology — the free row convects with strong noisy induced velocity just behind the TE) governs. Caveat: two changes vs nofloor (η 0.2→0.1 AND +1 free row); contributions not separable from this run alone. |
| `p2e_nrows4_das1p0` (η=1.0, 4 rows) | 12927923 | **0.07431 ± 0.00096** | NOT extent-insensitive: +4.2% over the 1-row 0.07133, above cycle scatter AND above the band. More panel-sheet before particle conversion ⇒ more CT. |
| `p2e_das1p0_refresh` | 12927924 | CANCELLED at step 398 | — |

**Reading (sobering):** CT rises monotonically with how much near wake is carried
as panel sheet rather than particles — rigid `Das` length (steep at small `Das`),
free-row extent (+4.2% for 3 rows), limiting at 001's semi-infinite rigid sheet
(0.110). Experiment (0.072) sits inside that range, so the η=1 in-band value is
at risk of being an interpolation artifact of the sheet/particle split rather
than a converged prediction. The needed demonstration is a **plateau in BOTH
`Das` and `nwakerows`** (does the nrows ladder 1→4→8→… asymptote, and where?),
or a principled near-wake treatment. NT=72 (12921559, step 1213/1439) still
pending as the dt-refinement check with the F5 pre-registration (0.0708–0.0712).

### 2026-07-29 (later) — final-filament review CLEAN; wing sheet/particle ladder RUNNING; rotor nwakerows-up ladder SUBMITTED

- **Final-filament check (Ryan):** verified in source that the sheet-end
  filament is ACTIVE in every particle-shedding configuration. Rotor driver:
  `include_final_filament` defaulted `true` (`get_sources` includes the
  `FilamentWrapper`), paired deliberately with `unsteady_filament=false` +
  `NoShed` (filament cancels the current last row). Wing precedents
  (pitching_wing particle branch, task1 consistency): default `true` with
  `unsteady_filament=true` + OverlapPPS. The `include_final_filament=false`
  instances are all pure-panel potential oracles, where the last row's own ring
  edge IS the physical starting vortex. No fix needed.
- **New wing experiment** (`examples/ssw_sheet_particle_split.jl` + particle
  mode in `suddenly_started_wing.jl`): panel_rows ladder {2,4,8,16,32,64} +
  pure-panel control at AR=6, η=1, dt*=1/8, t*=20, plus a rotor-like-σ
  contrast wave (OverlapPPS overlap 3.0 at rows 2/8/32). Running locally.
- **Bug found & fixed while smoking it:** `_shed_particles!` admitted
  zero-strength particles; their |Γ|=0 NaNs corrected-Pedrizzetti relaxation
  and poisons the whole field. Exact zeros arise generically on symmetric wings
  (mid-span trailing difference) and impulsive starts (zero-strength first
  row) — the rotor never hits exact zeros, which is why this never surfaced.
  Guard added in `src/FLOWPanel_wake.jl` (`Γ == 0 && return nothing`);
  wake unit tests pass; rotor results unaffected (no exact zeros there).
- **Submitted** rotor sheet-extent asymptote ladder at η=1.0 (unfiltered, stock
  relaxation): `p2e_nrows8_das1p0` (12943025, 24G), `p2e_nrows16_das1p0`
  (12943026, 24G), `p2e_nrows36_das1p0` (12943027, 32G; one revolution of
  sheet). Convergence claim requires CT(N) to plateau; N=1 → 0.0713,
  N=4 → 0.0743 so far. NT=72 (12921559) still running (~step 1400/1439).

### 2026-07-29 (evening) — NT=72 harvested: 0.07337, OUTSIDE the F5 pre-registration; matched-physical-units NT=72 submitted

`p2e_nt72_das1p0` (12921559, 1439 steps, COMPLETED): 10-rev cycle-mean
**CT = 0.07337 ± 0.00063** (max within-rev p-p 0.0045) vs NT=36's
0.07133 ± 0.00159. The F5 log-law-only prediction (0.0708–0.0712) FAILED, in
the **+** direction — but attribution is confounded exactly as 006 study D
warned: at fixed η, NT 36→72 halves frozen `Das` (0.41→0.20 local chords, edge
of the amplified region), halves σ/chord (1.5→0.75 — the wing ladder's most
sensitive range), and halved `rlxf`. Notable: ripple improved sharply (per-rev
std 0.9% vs 2.2%).

Submitted the clean isolation `p2e_nt72_das2p0_ov6` (**12950996**, 48 h): NT=72
with `Das` pinned (η=2.0 ⇒ 0.41 chords again) and σ pinned (OVERLAP=6 ⇒ same
σ, halved particle spacing). If IT reproduces ~0.0713, dt is converged and the
naive NT=72 shift was the `Das`/σ physical changes; if it lands near 0.0734,
the timestep itself was the mover and NT=36 is under-resolved in time.

### 2026-07-30 — sheet-extent ladder COMPLETE: non-monotone, peak at N=8, and N=36 is the most stable config measured; N=72 submitted

10-rev cycle-means at η=1.0 (velocity, NT=36, 4R, unfiltered, stock relaxation):

| `nwakerows` | sheet extent | CT | rel std |
|---:|---:|---:|---:|
| 1 | 1/36 rev | 0.07133 ± 0.00159 | 2.23% |
| 4 | 1/9 rev | 0.07431 ± 0.00096 | 1.29% |
| 8 | 2/9 rev | 0.07506 ± 0.00160 | 2.13% |
| 16 | 4/9 rev | 0.07304 ± 0.00203 | 2.78% |
| 36 | 1 rev | **0.07049 ± 0.00047** | **0.67%** |

**The "CT rises monotonically with sheet fraction" reading (based on N=1→4) is
FALSIFIED by the fuller ladder**: CT peaks at N=8 (0.0751) then declines through
N=16 to N=36 (0.0705). It does NOT head toward the semi-infinite 0.110.
Consistent with the wing σ-ladder's sign structure: chord-scale-σ particles
OVER-induce very near the body and under-induce farther away — as the
conversion point moves aft, the over-inducing near-body particles are replaced
by sheet (CT drops back), and the remaining particle error shrinks.

**N=36 (one revolution of sheet) is the best configuration measured to date on
BOTH 006 axes**: CT = 0.0705, inside the 0.068–0.072 band, with rel std 0.67%
(tol 0.5%) and within-rev p-p 3.28% (tol 2%) — 3× quieter than N=1. Submitted
`p2e_nrows72_das1p0` (**12955430**, two revolutions of sheet) to test whether
0.0705 is the asymptote: N=36→72 agreement within scatter would be the
convergence demonstration this study needs.

Still running: `p2e_sigF_nofilt` (12943696, ~12 h) and the matched-physical-
units NT=72 (12950996, ~10 h).

### 2026-07-30 (Ryan ruling) — N=36 not adopted; targeted wing study directed

The nwakerows ladder is non-monotone (peak N=8, decline to N=36) — not a
convergent metric — so N=36's in-band low-ripple result is NOT adopted.
Directed next step: a targeted straight-wing steady-conditions study to isolate
the cause of the sheet/particle variability (measure induced-velocity fields,
not just CL, at rotor-like σ/chord). Detail + in-flight-job harvest table
(12955430 nrows72, 12943696 sigF_nofilt, 12950996 nt72match) in
BRAINSTORM/014, "Ruling 2026-07-30".

## 2026-07-31 — Final Phase 2e harvest (handoff to BRAINSTORM/018)

- `p2e_nrows72_das1p0` (12955430, COMPLETED): CT̄(revs 10–19) = 0.06931 ± 0.00154 (2.23%). The nwakerows tail continues falling (36 → 0.07049, 72 → 0.06931) — non-monotone-ladder rejection confirmed.
- `p2e_sigF_nofilt` (12943696, OUT_OF_MEMORY at step 403/719): ReqMem was 32G, MaxRSS 33.4 GB; CF tail smooth ≈ 0.066 with no |CT|>1 ⇒ genuine memory ceiling, NOT the divergence pattern. σ-refined unfiltered wake was STABLE while it ran.
- `p2e_nt72_das2p0_ov6` (12950996, COMPLETED 2026-07-31 after 41 h): CT̄(revs 10–19) = 0.06852 ± 0.00157 (5-rev block drift −1.14%, window not fully settled). This is the clean dt isolation (Das AND σ pinned at NT=36 physical values; metadata-verified OverlapPPS 6.0/2, nwakerows=1): **−3.9% vs the NT=36 partner `p2e_das1p0` (0.07133)** ⇒ NT=36 is under-resolved in time and true dt refinement moves CT DOWN — the naive NT=72's 0.07337 was dominated by its halved-Das confound, as F5 suspected. Phase 3 of 018 goes to the full NT ladder incl. NT=144.

The publishable convergence campaign continues as **BRAINSTORM/018** (`p018_*` case tags in the same launcher); this ledger closes for new cases.
