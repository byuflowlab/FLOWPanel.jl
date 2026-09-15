# 026 reset prompt — de-risk harvest → wave-2 gating (2026-09-15)

You are picking up BRAINSTORM 026 (resolution-preserving particle splitting)
in `/Users/ryan/Dropbox/research/projects/FLOWPanel.jl` (+ siblings
FLOWVPM.jl, FastMultipole). Read `CLAUDE.md` and the policies it names
(HPC.md before any cluster work; `ssh orc` needs a live ControlMaster
socket — ask Ryan to run `! ssh orc echo ok` if 2FA blocks you).

## What happened (2026-09-12 → 09-15, all Ryan-authorized)

1. **Task 1 verdict (smoke B σ-check)**: the 2.7× retention drop under
   MERGE_OVERLAP=3.5 is the overlap gate consuming freshly split
   CROSS-FAMILY children (67% of merges are split×split pairs at
   σ ≈ 0.5–0.9 σ_shed, 98% fresh σ≈σ₀, ~100 merges/step), NOT aged-wake
   thinning; sibling re-merge 0–3%, at-release <1%, σ-pump ~5–9%. Regime
   caveat: that smoke sheds σ = 0.26R where the overlap gate is uniformly
   ~3.7× wider than absolute; the campaign s9 arms (σ ≈ 0.0381R) have the
   opposite young-particle polarity, so B does not transfer — the A/B
   stands. Evidence: `sigma_check_evidence_20260912/` (this dir; merge-event
   log gz + analysis scripts) and `data/smokeB_vtk_322-467_20260912.tar.gz`
   (B's only surviving VTK — smoke C overwrote steps 0–321 in the shared
   run dir). Full ruling record: `campaign_prep_20260912.md` §5.

2. **GPU splitting DONE + verified** (was "device accumulators only" gap):
   FLOWVPM `edc9d95`+`bf88806` (ResolutionSplitState storage matches the
   pfield array type; broadcast accumulation in euler/euler_exp/rk3 +
   3 CoreSpreading twins, pre-clamp semantics, NaN-safe; device-safe
   lifecycle hooks — shedding calls add_particle ON the device field);
   FLOWPanel `a804a95`+`035f50b` (`_gpu_sync_rsplit!` seam sync with
   WIDENED H2D vs stale device tails; device-canonical state; verification
   case defs scr_p026{gpuv_split,gpuv_splitmerge,cpuv_split} + campaign
   case scr_p026s9_exp_split_mo35; NaN-gate-safe splitting banner). Tests:
   rsplit 1055/1055 (t12 parity ×44), wake/replay clean. CUDA smokes:
   13689273 gate PASS (all 3 mechanisms fired on GPU, 7.6 vs 176 s/step
   CPU ≈ 23×); CPU twin 13689275 split-rate parity (34 lines both, per-step
   events within ~10%); 13689274 (MERGE_OVERLAP=3.5) died on the euler_exp
   substep-budget guard during an elongate storm — physics of the ignition
   case, not a code failure. All commits PUSHED.

3. **De-risk wave LAUNCHED and COMPLETED** (Ryan rulings 2026-09-14:
   3-arm de-risk; cold-start GPU; f_visc=0.587 enabled; SIGMA_FLOOR_FRAC=0.1
   with 0.25 fallback). Tag `campaign/p026-derisk-20260914` on all three
   repos (FLOWVPM bf88806 / FLOWPanel 035f50b / FastMultipole ac7230a6),
   worktrees + env at `orc:~/campaigns/p026-derisk-20260914/`, provenance =
   `p026_derisk_20260914_provenance.md` (this dir, mirrored on orc). Jobs
   (m13h H200, all COMPLETED overnight, **gate_rc=0**, banners verified —
   f_visc=0.587 f_comp=0.73 f_elong=0.3, Φ_t 2.75 cap / 2.4 exp):
   - 13691080 `scr_p026sp_nt144_cap030` NREVS=20 → ran 1295 steps
   - 13691081 `scr_p026s9_exp_split` NREVS=12, vatistas → 323 steps
   - 13691082 `scr_p026s9_exp_split_mo35` (MERGE_OVERLAP=3.5 twin) → 323 steps

## YOUR TASKS

### Task A — harvest the de-risk trio (start here)

Outputs live in the CAMPAIGN WORKTREE
`orc:~/campaigns/p026-derisk-20260914/FLOWPanel.jl/data/<case>/` and logs in
`.../logs/slurm/`. Per the provenance data policy: MOVE each run dir to the
shared root `~/projects/FLOWPanel.jl/data/` and leave a symlink behind,
then harvest (delegate scraping to `harvester`/`hpc-monitor`):

1. **cap030 cliff question (CHECK FIRST)**: the arm ran 1295 steps at
   NREVS=20 (≈64.75 steps/rev — NOT the 144 steps/rev I assumed; the
   original cliff bracket "steps 2200–2248" may sit at rev ≈ 34.7 in this
   driver's stepping). Determine whether the FMM adequacy stress / cliff
   window was actually reached (check `sigma_max / current_adequacy_limit`
   log lines, s/step trend, split telemetry). If not covered, propose a
   chained extension (RESTART_STEP from the retained state, NREVS≈36) —
   Ryan-gated.
2. exp_split vs exp_split_mo35 **merge A/B**: particle counts, σ census,
   split/merge telemetry, CT trend, s/step; did the overlap gate reproduce
   the churn/σ-pump signatures (cf. Task 1 verdict) at campaign σ?
   WAKE_HEALTH is on — health CSVs are in the run dirs.
3. Floor 0.1: any sign the lower guard floor destabilized the healthy
   phase (compare against the ef bracket history / expguard verdict
   +0.39% CT memory).
4. Score against the acceptance block in the provenance file; report
   PASS/FAIL to Ryan. **If PASS → he ruled the remaining 11 arms launch
   with floor 0.1 kept; if FAIL → floor back to 0.25.** Wave-2 submission
   itself remains Ryan-gated (present the plan + cost first; reuse the
   same tag/worktrees if code is unchanged, else new tag).

### Task B — hygiene (quick)

- Verification worktree `orc:~/wt026gpu/` (FLOWPanel@035f50b, FLOWVPM@
  bf88806 + env): smoke run dirs wrote THROUGH its `data` symlink into the
  shared root (`scr_p026gpuv_split`, `scr_p026gpuv_splitmerge`,
  `scr_p026cpuv_split` + logs in the worktree). After confirming smoke data
  is harvested/kept, ask Ryan whether to remove the wt026gpu worktrees
  (`git worktree remove`) — they were verification-only.
- `data/p026_restart_gpu40_s950` on orc was re-extracted from
  `/nobackup/archive/.../p026_restart_gpu40_s950.tar.zst` (archiver had
  stripped it); leave for the archiver.
- Local: `data/rotor_hover_pressure_comparison.metadata.toml` stays
  uncommitted (perpetually rewritten). BRAINSTORM files from this arc are
  uncommitted — bundle them into the next approved commit
  (campaign_prep_20260912.md, provenance, evidence dir, this file).

### Owed / parked

- **Notebook**: Ryan said "not yet" (2026-09-14) — the whole 026 arc
  (σ-check verdict, GPU splitting, de-risk wave) is unlogged; offer an
  entry at the next milestone. Also owed from earlier arcs: expguard
  three-arm result, Cd transient, rlxf derivation, Ladder C forensics.
- Remaining 11 campaign arms (§8.4 cap018 + 12−2 s020v arms minus the two
  run) per `campaign_prep_20260912.md` §2 — after Task A gating.
- f_visc=0.587 rationale (4^(1/3)−1 count-matched tetra4 analog) was
  Ryan-picked from options; it fired rarely (viscous=1–2 events/smoke) —
  keep an eye on it in harvest.

## Ground rules (carry-over)

Local ≤4 threads; long local jobs via `nohup ... & disown` + log polling
(harness kill bug), never harness run_in_background for the sim itself.
Slurm in non-login shells needs `bash -lc` (PATH); MOTD banners contaminate
ssh output — filter (`grep -E '^13...'` / `-oE 'COMPLETED|FAILED|...'`).
sacct state is not evidence; judge runs by outputs (.err before .out).
Commits, campaign launches, and notebook writes are Ryan-gated. Design doc
append-only; theory doc rewritten in place. Known unrelated failure:
`runtests_unit_warmstart.jl` first testset. The dispatcher wipes
`data/$RUN_NAME` on cold runs — never pre-place symlinks there.
