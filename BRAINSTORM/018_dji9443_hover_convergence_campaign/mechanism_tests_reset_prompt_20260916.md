# RESET PROMPT — 018 mechanism tests: monitor wave-1 tail + geometry analyses (2026-09-16)

You are a clean-context agent on BRAINSTORM item 018 (DJI-9443 hover
convergence, NT-ladder non-convergence: CT̄ climbs +1.27%/NT-doubling in the
g25 guarded pair). Wave 1 is mostly scored; your jobs are (A) monitor/score
the three in-flight runs, (B) run the approved free analyses 1–3 below, and
(C) prepare — but do NOT submit — GPU tests 4–6, which run only on Ryan's
explicit approval. Authoritative context, read in order:
`mechanism_tests_status_20260916.md` (this directory — wave-1 results,
verdicts, incident log), `mechanism_tests_provenance_20260915.md` (pins,
arms, env), and skim `gamma_distribution_status_20260915.md` (the original
Γ signature) and `mechanism_tests_reset_prompt_20260915.md` (launch charter:
scoring recipes, ORC gotchas). Standing rules: never tune to CT_exp;
notebook/ledger writes need Ryan's approval (stale since 2026-09-08; backlog
now includes wave-1 launch + warmstart incident + gate PASS + node-fail
resubmits + wave-1 partial verdicts + the late-onset finding below).

## In-flight jobs (poll `sacct`; prior session's watchers died)

`sacct -j 13733312,13733313,13733314 -X --format=JobID,JobName%24,State,Elapsed`

| job | test | run name | status at reset |
|---|---|---|---|
| 13733312 | T1a NT36 rev 30→60 | `p018_csarc_l3p0_3r_g25_s2` | RUNNING mgh-1-1 since ~19:20 MDT 09-16; banner VERIFIED (NT:36, rlxf:0.3, settle:52, resume step 1079, DynamicSFS rlxf=0.005, guard=on, mechtests repo). Wall 08:00. |
| 13733313 | T1b NT72 rev 30→60 | `p018_csarc_n2_nt72_l3p0_3r_srlx_g25_s2` | RUNNING mgh-1-2; banner VERIFIED (NT:72, rlxf:0.16334, settle:52, resume step 2159, DynamicSFS rlxf=0.0025031, guard=on). Wall 16:00. |
| 13733314 | T5b inviscid NT72 (resubmit after NODE_FAIL 13712097) | `p018_csarc_n2_nt72_l3p0_3r_srlx_nv_g25` | PENDING (Resources; starts when T1a's node frees). **Banner-verify on start**: `visc:false`, NO CoreSpreading line, NT:72, rlxf 0.16334, DynamicSFS rlxf=0.0025031, guard=on. Wall 14:00. |

History (do not re-diagnose): 13712154 T1-chain failed only because job-env
python3 lacks numpy (gate itself PASSED when run manually: max rel ΔCF
2.1e-5 vs tol 2e-4, wake state identical — restart machinery is validated);
13712097 T5b was a genuine NODE_FAIL (partial dir preserved at
`data/…_nv_g25_nodefail13712097`, delete after 13733314 lands cleanly).

## Wave-1 verdict summary (details in mechanism_tests_status_20260916.md)

- **T2 (merge cadence)**: far-wake Σ|Γ| divergence is a merge-cadence
  artifact (−98% of the rev-29 excess) but **NOT the CT carrier** (climb
  +1.23% vs baseline +1.27%). Γ signature and CT climb are DECOUPLED.
- **T4 (ConstantSFS 0.14)**: far-wake divergence 3× WORSE, but CT climb
  shrinks to +0.85%/doubling → SFS dynamics carry ~1/3 (dissipation-mismatch
  caveat: uniform Cs vs 54%-clipped dynamic model, counts −40%).
- **T5a (inviscid NT36)**: far wake +29%, CT +0.15%; pair verdict needs T5b.
- **NEW (2026-09-16): the CT climb is LATE-ONSET.** NT72−NT36 CT̄ by window:
  −3.8% (revs 8–12), +0.5% (11–15), −0.3% (16–20), +1.2% (21–25), +1.4%
  (26–30). Not a fixed evaluation offset; it develops revs ~15–30, two-sided
  (NT36 drifts down, NT72 up). Surviving suspect class: wake-borne, slowly
  developing, but NOT total far-wake strength ⇒ prime candidate is **wake
  geometry/structure** (induced inflow), plus SFS dynamics, plus whatever T1
  says about equilibration.

## Task A — monitor and score the in-flight runs

1. Banner-verify 13733314 as it starts (log
   `orc:~/wt018/FLOWPanel-mechtests/logs/slurm/slurm-fp-018mech-t5b-nv72-13733314.out`).
   Wrong banner ⇒ scancel immediately.
2. When runs land: judge health from monitors CSVs + the last `GATE:` line
   (gpu_gemv>0, nan_lines=0), exit 0 ≠ health.
3. Score PROMPTLY (sweeper culls landed runs to newest-36 VTPs; T1b has 4320
   steps):
   - Γ time series: `python3 ~/p018_mech_tests_20260915/p018_gamma_ts2.py <NT> <run>`
     (parameterized copy; NOTE: its `REVS` list and the `step >= 30*nt` cap
     are hardcoded for 30-rev runs — for the `_s2` runs edit a copy to cover
     revs 30–60, cap 60*nt).
   - Phase-avg per-bin ratios: `p018_gamma_dist_avg2.py RUN36 RUN72`
     (argv-parameterized; `WINDOWS=[15,23,29]` hardcoded — for T1 edit to add
     late windows, e.g. 44 and 59). Pair T1a×T1b; pair T5a×T5b
     (`p018_csarc_l3p0_3r_nv_g25` × the T5b run).
   - CT̄ windows: `cd ~/projects/FLOWPanel.jl && python3 scripts/p018_analyze.py
     m1 --revs A B <run>` — T5b: 21 25 / 26 30 / 21 30; T1: also 31 40 /
     41 50 / 51 60 (both rungs).
4. T1 verdict question: does NT72's far wake saturate by rev ~45–60 and does
   the CT climb SHRINK in late windows (transient-length mechanism) or
   persist (persistent carrier)? Given the late-onset finding, this is now
   the central arbiter.
5. Update `mechanism_tests_status_20260916.md` (append a dated section) or
   write `mechanism_tests_status_<date>.md`.

## Task B — free analyses (approved by Ryan 2026-09-16 chat: "do 1-3")

All on existing data; write scripts to `orc:~/p018_mech_tests_20260915/`,
never cat VTP bytes, print summaries only. Reader helper:
`sys.path.insert(0, "/home/rander39"); from vtp_C import read_vtp` (see
`p018_gamma_ts2.py` for the pattern; fields: `Points`, `gamma` (Nx3);
axial axis is **x**, R=0.11995 m). Baseline pair dirs under
`~/projects/FLOWPanel.jl/data/`: `p018_csarc_l3p0_3r_g25` (NT36) /
`p018_csarc_n2_nt72_l3p0_3r_srlx_g25` (NT72), full VTP sets on disk (also
tarred in `/nobackup/archive/usr/rander39/FLOWPanel_runs/`). Compare late
revs (e.g. phase-avg over rev 29→30, 12 phases, matched blade phase:
step = r0*NT + (NT//12)*j).

1. **Wake geometry profiling** (the big one): per rung, phase-averaged
   Γ-weighted structure of the near/mid wake — (a) slipstream contraction:
   Γ-weighted radial centroid r̄(x) in x-slices (e.g. Δx = 0.25R, x ∈
   [0, 3R]); (b) tip-vortex trajectory: in each x-slice, r and x of the
   |Γ|-densest cluster (or Γ-weighted r percentiles p25/p50/p75 as a robust
   proxy); (c) axial spacing of the tip spiral near first passages if
   extractable. Readout: does NT72's wake sit geometrically different from
   NT36's within the first 1–2 tip passages (x ≲ 1R)? That region sets
   induced inflow at the disk.
2. **Radial rebinning of Γ**: same avg machinery, but bin Σ|Γ| in cylindrical
   r/R (bins e.g. 0–0.2,…,1.2–1.4R) restricted to x < 1.5R, NT72/NT36
   ratios. Total near Σ|Γ| is NT-invariant — look for radial REdistribution.
3. **Spanwise loading split**: check what the monitors carry
   (`data/<run>/monitors/` — force monitor is
   `<run>_monitor02_force_system1.csv`; look for sectional/spanwise files or
   per-panel outputs). If sectional data exists: where along the blade does
   the late +1.4% live (tip-concentrated vs uniform)? If none exists, report
   that and move on (adding a monitor is a wave-2 code change, not yours).

Deliverable for B: dated results section (tables) + interpretation w.r.t.
the geometry hypothesis, in the status file.

## Task C — GPU tests, prep only, submit ONLY on Ryan's explicit approval

4. **Wake-transplant cross-restart** (strongest discriminator): restart NT72
   stepping from the NT36 rev-30 wake state (and mirror), run 2–3 revs.
   CT snaps to stepping-rung value within ~a rev ⇒ per-step near-field
   evaluation carrier; starts at donor value and drifts slowly ⇒ wake-evolution
   carrier. FEASIBILITY FIRST: arms differ in pps (12 vs 6) and nwakerows
   (1 vs 2) — audit `src/FLOWPanel_warmstart.jl` reconstruction
   (`RESTART_STEP/NAME/PATH` env; quadruplet `<name>_body1.<S>.vtu` +
   `<name>_wake1.{1,2}.<S>.vts` + `<name>_wake1_particles.<S>.vtp`) for
   cross-arm validity (particle state is arm-agnostic; panel-wake rows may
   need regeneration; NT36 arm has only wake1.1). Present findings + proposed
   sbatch lines to Ryan.
5. **Single-knob restart perturbations**: from existing rev-30 states, ~5-rev
   restarts flipping ONE of: MERGE_EVERY (done, =T2), pps, nwakerows,
   RELAX_RLXF (Pedrizzetti particle-Γ realignment factor, banner `rlxf`;
   srlx scaling 0.16334 = 1−(1−0.3)^½ matches compounded per-rev relaxation
   but the operation is nonlinear so residual NT-dependence is possible),
   SFS_RLXF (DynamicSFS coefficient under-relaxation, same ½-power scaling).
   Caution: pps/nwakerows are case-arm constants — per the standing ops rule,
   new case arms (examples commit), never env-overrides of unconditional
   exports. Watch the late-window CT response.
6. **Per-rung matched ConstantSFS**: Cs=0.145 (NT36) / 0.131 (NT72) — if the
   climb returns to ~1.27% the T4 reduction was SFS *level*; if it stays
   ~0.85% it's the *dynamics/fluctuations*. Reuse the T4 submit pattern
   (SubmitLine of 13712094 in sacct) with per-rung SFS_CONST_CS.

Wave-2 framing for Ryan: T3 (merging off) is pre-approved and triggered by
T2-positive, but mainly confirms the (already-explained) Γ artifact; flag
particle-count feasibility at NT72 (MAX_PARTICLES=1.5e6, baseline peaked
~340k WITH per-step merging). The CT-carrier axis (4–6 + T1 readout) is
likely higher value — Ryan decides.

## Ops (compressed; full list in the 0915 charter + ops_reference.md)

- `ssh orc 'bash -lc "…"'` for slurm/python; heredoc `ssh orc 'bash -ls' <<'EOF'`
  for quotes/$; MOTD glues to line 1 — lead with `echo`.
- Job-env python3 has NO numpy — any in-job scripted check must load a
  numpy-capable python first (login node python3 is fine).
- mgh = 2 nodes × 1 GH200. Submit from `~/wt018/FLOWPanel-mechtests`
  (`logs/slurm/` is relative), `--constraint=arm`. Campaign tag
  `campaign/p018-mech-tests-20260915` (FLOWPanel `5cdf058`); env
  `~/p018wtenv-expguard-gh200`. Run data in `~/projects/FLOWPanel.jl/data/<run>/`.
- Never cat VTP/CSV bytes; `vtk_protect_list.txt` is Ryan's — never write.
  Local checkouts are DIRTY with other sessions' work — read-only; ≤4 local
  threads; cluster campaign work only in `~/wt018/*` worktrees.
- Scoring scripts/logs live in `orc:~/p018_mech_tests_20260915/`
  (gate_check.py, t1a.slurm.sh, t1chain.slurm.sh, p018_gamma_ts2.py,
  p018_gamma_dist_avg2.py, gamma_ts_wave1.log, gamma_avg_wave1.log,
  ct_windows_wave1.log, ct_windows_baselines.log, sfs C extraction in
  `~/p018_gamma_dist_20260915/`).
