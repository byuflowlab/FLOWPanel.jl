# Handoff prompt — 026 §16 expint-fails reruns, scoring leg (written 2026-09-05 ~11:45 MDT)

Continue BRAINSTORM 026: babysit and score the four `scr_p026ef_*` s020v rerun arms on orc, record the verdict as design doc **§17**, apply Ryan's linegauss ruling if it fires, and report. Read `BRAINSTORM/026_sigma_growth_particle_splitting/particle_splitting_design.md` §15–§16 (RK3 verdict + shortlist + launch record), memory `project_026_particle_splitting.md`, and `agent_policies/HPC.md` first. Delegate CSV scraping to `harvester`, status to `hpc-monitor`.

**Context (one paragraph):** Phase 1b closed with §15: RK3 falsified (no arrest, no delay; both arms died PARTICLE OVERFLOW; 5–10× cost) — the lever is euler_exp, GPU-capable since Task 1. §16 shortlisted "expint-fails" candidates; the star is the 020 Phase-2R σ/R=0.02 screen where the corrected map still blew up (resolution-loss regime, pre-onset M≈28). Ryan approved cold reruns of that pair on the current stack, plus LineGauss twins, all retaining full VTP series so we finally have warm-startable states for this event.

**In flight (launched 2026-09-05 ~08:47, m12 --qos=normal, 12 h wall, from `~/wt026/FLOWPanel.jl` with `P018_REPO=$PWD`, worktree project julia 1.12.6):**

| arm | case | job | knobs delta |
|---|---|---|---|
| ctrl | `scr_p026ef_ctrl_s020v` | 13591760 | — (euler, Vatistas) |
| exp | `scr_p026ef_exp_s020v` | 13591761 | `WAKE_EXPINT=true` |
| ctrl-lg | `scr_p026ef_ctrl_s020v_lg` | 13591762 | `FLOWPANEL_FILAMENT_REG=linegauss` |
| exp-lg | `scr_p026ef_exp_s020v_lg` | 13591763 | linegauss + `WAKE_EXPINT=true` |

All four banners VERIFIED (expint/regularization correct). Cases are clones of `scr_p019_s020v` (OVERLAP 2.4, PPS 21, MERGE_R 0.00275, NWAKEROWS 1, DSIGMA 3.4, CoreSpreading β=1e9, WAKE_HEALTH_DTZ+ATTRIBUTION) — added to the dispatcher `examples/run_p018_screen_hpc.slurm.sh` both locally (UNCOMMITTED) and in the cluster worktree. 323 steps each. At 11:17 MDT: steps 135–169, ~78–108 s/step and rising with particle count. Expect ctrl ignitions ~step 200–250 (original refs: euler ctrl u>100 @213, corrected-map @242; dep stack moved 08-24, timing may shift — that re-measurement is the point). Ignited arms likely die by PARTICLE OVERFLOW (500k) like the RK3 arms. ETA: deaths early-mid afternoon; a surviving exp arm finishes ~16:30–18:30 MDT. Logs: `~/wt026/FLOWPanel.jl/logs/slurm/slurm-fp-p026ef-{ctrl,exp,ctrl-lg,exp-lg}-135917{60,61,62,63}.{out,err}`; run dirs `~/projects/FLOWPanel.jl/data/scr_p026ef_*`.

**Do, in order:**

1. **Babysit** (prior session's monitors die with it — re-arm: poll each run's `monitors/*_monitor04_wake_health_system1.csv` last step + grep .err for `ERROR|DomainError|OutOfMemory` every ~20 min; squeue needs the full path `/apps/slurm/latest/bin/squeue` in non-interactive ssh). Judge by outputs, never sacct.
2. **Score as §13/§15**: per arm — max_u, max γ/σ², min σ, max dtZ, first step max_u>100 (+ first dtZ>2/3); ctrl-vs-exp ignition timing per family. CT_bernoulli (`data/<run>/<run>_CT_vs_rev.csv`) exists ONLY for arms that finish; crashed arms have no CT file (written post-`simulate!`, line ~1468 of the driver) — use monitor02 CFz same-step windows for load agreement instead.
3. **Rule per family**: does the current-stack expint arm still blow up (expint-FAILS confirmed on current code → validates §16 candidate 1 as the Phase-2 splitting motivation case) or does it arrest (the 08-12→now stack changes rescued it → §16 needs a new candidate)? Record as design doc **§17** (§13 format + pacing note + which VTP steps bracket each ignition for future warm starts). Report to Ryan BEFORE further 026 work.
4. **Ryan's ruling on record (2026-09-05, act on it):** *if the LineGauss pair showcases the failure (exp-lg blows up), the campaign default filament regularization changes to linegauss from then on.* If it fires, propose the concrete change (dispatcher line 43 pin + wherever else Vatistas is pinned) — implementation itself may still warrant his sign-off on scope.
5. After scoring, surface (don't act): archive the four `scr_p026ef_*` dirs once scored (hpc-storage); the still-pending items below.

**Also live right now:**

- **hpc-storage agent pass PENDING**: approved to archive the 11 `scr_p026ph1*` dirs (8 Phase-1 + 2 rk3 + smoke) but blocked on the checkout lock held by another running archiver (p018_csarc pass, PID 1023011 on login04). It waits and auto-resumes; /home was 525 GiB vs 400 G cap before either pass. If the session reset killed the agent, relaunch hpc-storage with the same 11-run `--only` list; it must NOT touch `p026_restart_*`, `scr_p019_*`, `scr_p020*`, `scr_p026ef_*`, or the protect list.

**Gotchas (all live):**

- `ssh orc` needs the ControlMaster socket; on 2FA ask Ryan for `! ssh orc echo ok`.
- Non-interactive ssh has no slurm on PATH (`squeue: command not found` caused one false "both jobs done" alarm) — use `/apps/slurm/latest/bin/squeue`; never scancel on a monitor's word alone.
- Step lines in .out begin with a TAB (`grep -a 'step .* at time'`).
- Wake-health header: `step,time,n_particles,max_u,min_sigma,min_sigma_ratio,max_gamma_over_sigma2,wall_s,max_dtZ,p1_sigma_ratio,argmin_*`.
- The archiver ate `data/p026_restart_*` once (restored); suggest protect-listing to Ryan — protect list is HIS file, never edit it.
- Do NOT score these reruns against the original 019/020 verdicts (dep stack moved 08-24); score ctrl-vs-exp within each family.
- Never tick notebook checkboxes; notebook writes need Ryan's approval.

**Pending Ryan decisions (surface when convenient, don't act):** notebook entry for Phase 0/1/1b (+§16/§17) — deferred twice, ask verbosity when offering; protect-listing `p026_restart_*`/`scr_p026ph1*`; keep or delete `~/wt026` worktrees + env-gpu after 026 closes; merge RK3 wiring + the `scr_p026ef_*` dispatcher cases to mainline (commit currently missing locally); GPU-preference memory saved 2026-09-05 (future long runs default GPU where sensible — see `feedback_gpu_launch_preference.md`).
