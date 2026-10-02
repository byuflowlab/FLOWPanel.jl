# 026 reset prompt — §22.3 rerun slate IN FLIGHT; watch arms, verdicts, harvest (2026-09-18b)

You are picking up BRAINSTORM 026 in
`/Users/ryan/Dropbox/research/projects/FLOWPanel.jl` (branch
**fastmultipole**) + sibling `/Users/ryan/Dropbox/research/projects/FLOWVPM.jl`
(branch **flowpanel**). Read `CLAUDE.md` and the policies it names
(`agent_policies/HPC.md` before touching ORC). Predecessor:
`merge_redesign_landed_reset_prompt_20260918.md` — its tasks are DONE.
This file supersedes it. Durable record for everything below:
`rerunslate_provenance_20260918.md` (pins, slate, submissions, banner
verdicts, outcomes so far) — trust it over this summary if they diverge.

## State (all as of 2026-09-18 evening)

1. **Committed + tagged**: §22 merge redesign FLOWVPM `8d4a3b4`; FLOWPanel
   doc bundle `8eff54e` + RUN_NAME_OVERRIDE launcher hook `6e52628`
   (both repos tagged `campaign/p026-rerunslate-20260918`, pushed to orc
   BY NAME — branches NOT pushed, divergence ruling still parked with
   Ryan). FMM pin unchanged `ac7230a6`. NOT pushed to github origin.
2. **Smoke dir** `data/smoke_mergesigma2m_20260917/`: VTK deleted (Ryan-
   approved), CSVs/monitors kept, 4.3 MB.
3. **Rerun slate LAUNCHED** (Ryan-approved 5 arms; NT144 cap-ladder rung
   PARKED). Campaign worktrees `orc:~/campaigns/p026-derisk-20260914/`
   at the tags, dirty=0, env Manifest verified. Run names are
   `scr_p026s9r2_*` (RUN_NAME_OVERRIDE; avoids wave-2 harvest collision).
   Slurm logs: campaign worktree `logs/slurm/slurm-fp-052-scr-gpu-<jobid>.out`.

| arm | job | case (run=…r2_<case-suffix>) | mech | state at handoff |
|---|---|---|---|---|
| A1 | 13758586 (m13h) | scr_p026s9_explg_fs | fs, uncapped | **DIED step 328/467**: DomainError dt·\|L\|=3461.6 (euler_exp budget) — same class as wave-2 explg_floor/split deaths; redesign delayed (+33–54 steps) but did NOT prevent exp ignition. Precursors: adequacy ratio warnings 310–328 with count-driven limit shrunk to 0.00446; 87k merge events. Banner was PASS 9/9. |
| A2 | 13763819 (eng) | scr_p026s9_explg_floor | floor only | **DIED step 274/467**: DomainError dt·|L|=5251.8 — EXACT same step as wave-2 twin (2133 @274), larger magnitude. Zero redesign effect on floor-only ⇒ Γ-side ignition channel independent of merge σ-pump. Banner was PASS 9/9. GOTCHA: traceback is in logs/slurm/…<jobid>.err, NOT .out (.out ends at GATE/dispatcher_rc=1). |
| A3 | 13763820 (eng) | scr_p026s9_ctrllg_fs | fs | RUNNING since ~18:58, banner PASS 9/9. Ctrl-channel probe: death ~260–290 again ⇒ ctrl blow-up is pump-independent. 12 h wall (wave-2 hit 8 h wall). |
| A4 | 13763821 (eng) | scr_p026s9_explg_fs (run …_fs_cap3) | fs + WAKE_SPLIT_SIGMA_MAX=0.0071 + SIGMA_CEIL=0.0071 (k=3 × shed σ 0.00237) | PENDING. Tripwire arm vs A1. Given A1's death, watch whether the cap ARRESTS the ignition (binding events) or the arm dies anyway ⇒ σ-cap can't stop Γ-side. |
| A5 | 13763822 (eng) | scr_p026s9_explg_split | split only | PENDING. Replicate of wave-2 DomainError @280 class. |

  IMPORTANT death-character correction (Ryan ParaView review + monitor04
  analysis; full section in the provenance file): the DomainError deaths
  are NOT field-wide blow-ups. The guard is max-over-particles
  (dt·|L|>2048 budget, `FLOWVPM_timeintegration.jl:647-713`); A1's
  global metrics stayed healthy while the extreme tail ignited (min_σ
  pinned at viscous floor from ~315, max Γ/σ² 826→9.1e4, max_u 49→919,
  steps 315→327) — localized Γ-side ignition near the root/fountain
  region. A1 also rode the 500k particle cap from ~step 318 (wave-2
  explg_fs died OF that cap @295 — the new merge law kept A1 alive at
  the ceiling). Use "localized Γ-ignition tripping the substep guard",
  not "blow-up".

  History note: A2–A5 were first submitted to m13h as 13758587–90,
  cancelled while PENDING and resubmitted on eng
  (`-p eng --qos=eng --gres=gpu:h200:1`, Ryan: "we get priority on eng").
4. **A1 last-50-step VTK** (steps 278–327 + monitors/merge_events/
   metadata, 3.6 G) downloaded to local `~/scr_p026s9r2_explg_fs_last50steps/`
   for Ryan's ParaView pass (.pvd files reference all steps — open the
   file series in the subdirs, not the .pvd). Ryan is examining the
   death; expect rulings from that. A2's identical-step death sharpens
   the picture: the exp/floor DomainError ignition is NOT merge-driven;
   the interesting remaining discriminators are A4 (does the k=3 σ-cap
   arrest it?) and A3 (ctrl channel).
5. The predecessor session's Slurm monitor died with it — re-establish
   your own watch (hpc-monitor subagent or a Monitor probing
   `sacct -n -X -j 13763819,13763820,13763821,13763822 -o JobID,State`
   via `ssh orc 'bash -lc ...'`; anchor awk on the job-id field — MOTD
   contaminates ssh output and non-login shells lack sbatch/sacct).

## YOUR TASKS (in order)

1. **Watch the four eng arms** (judge by outputs, never sacct state).
   On each start: banner verification per wave-2 checklist (12 revs/468
   steps, linegauss, mechanism knobs per arm — A3/A4 fs = floor 0.1 +
   split trio 0.587/0.73/0.3; A5 split-only = floor OFF; A4 additionally
   SIGMA_CEIL=0.0071 and clamp=[…, 0.0071] instead of Inf — merge log +
   sigma telemetry lines, r2 run name). Record verdicts in the
   provenance Banner table.
2. **On each terminal arm**: diagnose by outputs (log tail, first error
   signature, CT/monitor CSV state), append to the provenance Outcomes
   table. Key discriminations are in the arm table above.
3. **When all four are terminal**: harvest (delegate `harvester`) — MOVE
   run dirs from campaign worktree `data/` to shared root
   `orc:~/projects/FLOWPanel.jl/data/` + symlink back (NEVER pre-place
   symlinks). Assemble the slate verdict: per-event merge σ law vs
   survival; whether caps are needed and at what k; ctrl-channel
   status. Ledger line + offer (don't write) the notebook entry.
4. **Storage**: five VTK-writing runs — launch `hpc-storage` for a cycle
   once arms finish (A1's run dir is 11 G in the worktree; wave-2 dirs
   already in shared root).

## Owed / parked (carried)

- Notebook entry for the whole 026 arc (Ryan "not yet" ×2) — offer only.
- orc branch-divergence ruling (tags carry pins; do NOT force-push).
- Push of fastmultipole/flowpanel branches + tags to github origin —
  Ryan-gated.
- scr_p026gpuv_split archive retry via hpc-storage once quiet ≥24 h.
- NT144 cap-ladder rung: parked until this slate reads out.
- 021 silo cleanup owed (021 package §7); 018/022 queue jobs are other
  sessions'.

## Ground rules (carry-over)

Local ≤4 threads. `ssh orc` needs a live ControlMaster socket (ask Ryan
to run `! ssh orc echo ok` if 2FA blocks). Use `bash -lc` on orc for any
Slurm command. MOTD contaminates ssh output — filter. For any ORC
execution read `BYU_ORC_AGENTS.md` first. Commits, submissions beyond
this slate, and notebook writes are Ryan-gated.
