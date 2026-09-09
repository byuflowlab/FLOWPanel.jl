# Handoff prompt — 026 Phase 2: does particle splitting fix the s020v expint-fails event? (written 2026-09-05 ~22:20 MDT)

Continue BRAINSTORM 026 into Phase 2. Read `particle_splitting_design.md` **§9–§10 and §14–§17** (phase structure, split-geometry design state, RK3 verdict, expint-fails shortlist + launch, and the fresh §17 harvest/ruling), memory `project_026_particle_splitting.md`, and `agent_policies/HPC.md` first. Delegate status to `hpc-monitor`, CSV scraping to `harvester`.

**Context (one paragraph):** §17 (2026-09-05) validated §16 candidate 1 on the current stack: in the s020v resolution-loss regime the exponential integrator does NOT arrest — Vatistas expint ignited *earlier* than its euler ctrl (210 vs 225, died +1 step on NaN), LineGauss expint delayed 55 steps but still ignited (285) and died (Inf @300). All four `scr_p026ef_*` arms retained full VTP series with ignition bracketed at ±1 step (ctrl 224/226, exp 209/211, ctrl-lg 229/231, exp-lg 284/286) — warm-startable states for cheap splitting iteration. Ryan's linegauss ruling FIRED (exp-lg blew up): campaign default filament regularization is now linegauss — already implemented in BOTH dispatcher copies (`examples/run_p018_screen_hpc.slurm.sh:41`, local + `~/wt026`; cluster commit `54ad90e` on branch `p026-phase1b-rk3`, tag `campaign/p026ef-rk3-20260905`).

**In flight (babysit first):** `scr_p026ef_rk3_s020v_lg` — job **13593720**, m12/normal, **36 h wall** (RK3 is 5–10× euler cost, §15), submitted 2026-09-05 ~22:15 MDT from `~/wt026/FLOWPanel.jl` at the tag above. Same knobs as ctrl-lg + `WAKE_INTEGRATOR=rk3`. The prior session's monitor dies with it — re-arm (poll `monitors/scr_p026ef_rk3_s020v_lg_monitor04_wake_health_system1.csv` last step + grep .err every ~20–30 min). **Verify the banner** (`WAKE_INTEGRATOR=rk3`, `FLOWPANEL_FILAMENT_REG=linegauss`) — mandatory, not yet done (job was PD at handoff). Logs: `~/wt026/FLOWPanel.jl/logs/slurm/slurm-fp-p026ef-rk3-lg-13593720.{out,err}`.

**Do, in order:**

1. **Verify banner** of 13593720, then babysit. Reference for comparison: euler ctrl-lg ignited @230 (u>100), died @260 (overflow); expint-lg @285/@300.
2. **Score the RK3 arm as §13/§17** (max_u, max γ/σ², min σ, max dtZ, first u>100, first dtZ>2/3, death step+cause; loads = monitor02 **CFx** — thrust axis is x; CFz is near-zero lateral, never average it). Record as **§18**. This is a *confound check*: **if RK3 arrests the event, the failure is time-integration accuracy, not resolution loss — the Phase-2 splitting motivation weakens; STOP and report to Ryan before any splitting work.** If RK3 also ignites (like §15 on the gpu40/LG events), the resolution-loss reading stands and splitting proceeds.
3. **Prep the splitting test** (Ryan's framing: "let's see if particle splitting fixes it"): re-read §9 (shrink-side splitting validation case) and §10 (Phase 2 definition); check whether any split-geometry implementation exists yet (grep FLOWVPM worktree for split; Phase 2 was ON HOLD at §13, so likely design-only — if so, the deliverable is an implementation plan + validation-case wiring, not a launch). Use the s020v warm-start VTPs (brackets above) to iterate cheaply near ignition instead of 5–6 h cold runs; propose the arm matrix (ctrl-lg vs split-lg minimum; keep pairs backend-matched) to Ryan BEFORE launching.
4. After the RK3 verdict, surface (don't act): archive all five `scr_p026ef_*` dirs via hpc-storage (Ryan said archive after his review).

**Not skipped, still pending (surface when convenient, don't act):** notebook entry for Phase 0/1/1b + §16/§17 (deferred 3×; ask verbosity when offering); **local repo has uncommitted state** — dispatcher (ef cases + linegauss default + rk3 case) and design-doc §17 — plus the standing decision to merge the RK3 wiring + dispatcher cases to mainline (cluster branch `p026-phase1b-rk3` has them committed; local does not); protect-listing `p026_restart_*`/`scr_p026ph1*`; the blocked hpc-storage pass over the 11 `scr_p026ph1*` dirs; `~/wt026` worktree retention after 026 closes.

**Gotchas (all live, inherited from prompt_3):**

- `ssh orc` needs the ControlMaster socket; on 2FA ask Ryan for `! ssh orc echo ok`.
- Non-interactive ssh: slurm needs full paths (`/apps/slurm/latest/bin/squeue`); **sacct state is unreliable both ways** — this session sacct showed a RUNNING job as COMPLETED and produced a false "ctrl stopped" alarm plus a bogus subagent kill-theory. Judge by outputs only; never scancel on a monitor's word.
- Step lines in .out begin with a TAB (`grep -a 'step .* at time'`).
- Wake-health header: `step,time,n_particles,max_u,min_sigma,min_sigma_ratio,max_gamma_over_sigma2,wall_s,max_dtZ,p1_sigma_ratio,argmin_*`. Euler overflow deaths are intra-step: last CSV row can show far fewer than 500k particles.
- No `*_CT_vs_rev.csv` unless an arm finishes all steps (written post-`simulate!`).
- Do NOT score against original 019/020 verdicts (dep stack moved 08-24); compare within the current `scr_p026ef_*` family only.
- Never tick notebook checkboxes; notebook writes need Ryan's approval. Harvester subagents can botch column choice — spot-check any load table against raw monitor02 headers.
