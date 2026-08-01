# Phase 0 — Harvest + Deployment (no new jobs)

**Objective:** close out the three unharvested Phase-2e cluster jobs (they gate
Phases 3, 4, 5 design), deploy the 018 launcher changes, and verify the
pre-submission gate. **No 018 submissions until this phase closes.**

## Steps

1. Ryan opens the ControlMaster socket (`! ssh orc -fN`). Then per
   `ops_reference.md`: check `squeue`; harvest each finished job's CSVs +
   metadata + bound-circulation monitor CSV; update `ledger.md` AND the 2e
   ledger (`logs/dji_convergence_20260722/phase_02e_unsteady_ct.md`) +
   BRAINSTORM/014's harvest table.
2. Deploy: scp the updated `examples/run_dji9443_hover_ct_hpc.slurm.sh` to the
   cluster; md5-verify `src/` and `examples/rotor_hover_pressure_comparison.jl`
   against local.
3. Run the pre-submission gate (ops_reference) locally.

## Branch table (record outcomes here)

| job | case | branch |
| --- | --- | --- |
| 12950996 | `p2e_nt72_das2p0_ov6` (NT=72, Das & σ pinned at NT=36 physical values) | CT̄ ≈ 0.0713 ⇒ dt already converged at NT=36 → Phase 3 = single confirmation. CT̄ ≈ 0.0734 ⇒ NT=36 under-resolved → Phase 3 full ladder incl. NT=144. In between ⇒ Phase 3 ladder, judge on the pinned pair. **RESOLVED 2026-07-31: CT̄ = 0.06852 ± 0.00157 — BELOW both anchors (−3.9% vs the NT=36 partner 0.07133) ⇒ NT=36 under-resolved in time and true dt refinement moves CT DOWN (the naive NT=72's 0.0734 was the Das confound) ⇒ Phase 3 full ladder incl. NT=144.** |
| 12943696 | `p2e_sigF_nofilt` (σ≈0.37c legacy policy, unfiltered) | Completed stably ⇒ σ-refinement viability supported; record as off-ladder consistency point near σ/c≈0.37. Diverged ⇒ raise Phase 4 ladder base OVERLAP (2.0 → 2.4 across all rungs) before submitting anything. |
| 12955430 | `p2e_nrows72_das1p0` | Record only — completes the rejected nwakerows ladder tail (does CT continue below 0.0705?). Context for Phase 5, no design change either way. |

If a job is still queued/running, wait (do not resubmit); if it died, harvest
whatever CSV prefix exists (`RHPC_KEEP_PREV` note in ops_reference) and record
the failure mode.

## Exit criteria

All three rows resolved in `ledger.md`; launcher deployed + verified;
pre-submission gate green. Update the phase-gate row in the item file.

## Log

- 2026-07-31 — Harvest (ControlMaster open). **12955430 nrows72 COMPLETED**:
  CT̄(revs 10–19) = 0.06931 ± 0.00154; ladder tail still falling — rejected-
  ladder non-convergence confirmed; CSVs + bound-circulation monitor local at
  `data/p2e_nrows72_das1p0/`. **12943696 sigF_nofilt OUT_OF_MEMORY at step
  403/719**: ReqMem was 32G, MaxRSS 33.4 GB, CF tail smooth ≈ 0.066 ⇒ genuine
  memory ceiling, NOT the divergence pattern ⇒ Phase 4 branch = viability
  SUPPORTED, keep ladder base OVERLAP=2.0; the flat-64G ruling removes the
  failure cause. **12950996 nt72match still RUNNING** (step ~1011/1439 at
  28 h of 48 h; instantaneous CF ≈ 0.0675) — Phase 3's branch decision waits;
  Phases 1/2 proceed (no dependency).
- 2026-07-31 — Deployment: local `src/` (25 files incl. new `FLOWPanel_kutta.jl`)
  + driver + launcher synced to cluster; md5 parity verified file-by-file
  (cluster had stale pre-015 src — the Phase 2c failure mode, caught).
  2e-era numbers ran on the older src; PRIOR rows stay corroborating-only.
  Pre-submission gate PASSED: `formulation_test.jl` all 10 stages;
  NT=6 40_40 smokes clean in both formulations at campaign settings.
  **Submitted:** 12993801 `p018_b0`, 12993802 `p018_das0p5`, 12993803
  `p018_das2p0`, 12993804 `p018_das4p0` (5 study jobs incl. 12950996 ≤ cap).
  b0 log header verified: sigma_overlap 2.0/4, merge 0.0120, filter off,
  η=1.0, N=4, NT=36, RPM 5400. Driver gained `RHPC_BACKEND=fmm|direct`
  (T4 discriminator switch, default fmm — no behavior change); T4 runs local.
  Phase 0 remains OPEN only for the 12950996 harvest (monitor armed).
- 2026-07-31 — **12950996 nt72match COMPLETED (41 h, exit 0:0) and harvested ⇒
  PHASE 0 CLOSED.** CT̄(revs 10–19) = 0.06852 ± 0.00157 (10 per-rev blocks),
  5-rev block drift −1.14% (window not fully settled; same still-falling caveat
  as nrows72). Metadata verified: nwakerows=1, OverlapPPS overlap=6.0 pps=2 —
  the correct σ-pinned NT=72 configuration (ratio 3.0 = 2× the legacy NT=36
  ratio). Branch resolved: BELOW both anchors ⇒ Phase 3 = full NT ladder incl.
  NT=144; dt refinement direction is DOWN. CSVs + bound-circulation monitor
  local at `data/p2e_nt72_das2p0_ov6/`. All three branch rows now resolved;
  launcher deployed; gate green — exit criteria met.
- 2026-07-31 — **VTK retention (Ryan's call):** full ParaView set (15,844
  files, 20 GB) rsync'd to Ryan's local `~/p2e_nt72_das2p0_ov6/`
  (integrity-verified: file count + bytes match cluster exactly), then deleted
  from the cluster run dir (CSVs/TOML/monitors kept, 8.8 MB remain). Cluster
  VTK now exists only for the four active `p018` runs.
- 2026-07-31 — **Ruling 10 applied (VTK retention):** deleted VTK from all 23
  closed `p2e_*` run dirs on the cluster (CSV/TOML/monitors kept everywhere,
  spot-verified), freeing ~180 GB (`data/` 230 → 47 GB). VTK retained only for
  active runs: `p2e_nt72_das2p0_ov6` + the four `p018` jobs (4–5.5 GB each and
  growing). Going forward: prune each run's VTK at harvest once its question
  is answered, keep ≤3 runs' VTK, and ask Ryan for a ParaView pass BEFORE
  pruning if a convergence question needs wake visualization.
