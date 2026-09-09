# Handoff prompt — 026 Phase 2: finalize the splitting implementation plan with Ryan (written 2026-09-06 ~00:15 MDT)

Continue BRAINSTORM 026 Phase 2. This session's job is to **talk through and finalize the implementation plan with Ryan**, then (on his approval) implement it. Read first: the draft plan at `~/.claude/plans/shimmering-wobbling-sunrise.md` (context-complete — the implementing agent needs no further exploration), design doc `particle_splitting_design.md` **§18** (fresh RK3 ruling) and skim §3a/§3b/§4/§6/§7 for the spec, memory `project_026_particle_splitting.md`.

**Where things stand:**

- **§18 (2026-09-06): RK3 confound check ruled NO ARREST** — `scr_p026ef_rk3_s020v_lg` (13593720) ignited @74 (vs ctrl-lg 230, exp-lg 285), died @92 PARTICLE OVERFLOW. Order-of-accuracy falsified as the lever; the §17 resolution-loss reading stands; **Phase-2 splitting proceeds as motivated**. Recorded in the design doc; babysitters all shut down.
- **Implementation plan DRAFTED, NOT approved.** Ryan rejected plan-approval twice to discuss specifics (that's the current conversation thread, mid-stream). Rulings already made and folded into the plan:
  1. Implement BOTH mechanisms (§3a tetra-4 viscous split AND §3b tri-3 stretch split), OFF by default, env knobs.
  2. Mechanism B direction = **time-averaged stretching axis, averaged since the last split** (reset-on-split sign-invariant running director sum + weight in `SplittingState`), `STRETCH_AVG` enum member, STRENGTH fallback below `coherence_min`.
  3. Existing `SplitParticles`/`_do_split!` are experimental — do NOT build on them; new `ResolutionSplit` policy + `split_particles_026!` + `ResolutionSplitOptions`.
  4. **NO lineage_id/split_generation state** — fully random orientations (serial pass-2 loop, default RNG), not reproducible across warm starts, accepted for now, REVISIT note in plan.
  5. **RBF reset site untouched** — axis/weight NOT cleared there (flow history survives re-projection; site unreachable in production with β=1e9).
- **Last exchange before reset:** Ryan asked for a summary/definition of every `ResolutionSplitOptions` field and said "let's discuss first." The summary was delivered (see plan Stage 3.3 for the struct); two design points were flagged for his eyes, still UNANSWERED: (a) `coherence_min`=0.5 + STRENGTH-fallback is the only soft physics judgment in the struct; (b) severity ranking under `max_fraction` mixes grow and shrink candidates in one budget — offered a per-side quota if he wants shrink (ignition-critical) events guaranteed a share. Open decision points D1–D5 are listed at the plan's end (all env-settable defaults). Resume the discussion there.

**In flight on orc (check early):** detached worker resume-deleting the 8 Ryan-approved `p018_csarc_*` ARCHIVED-STALE runs (~253 GiB; each re-verified vs tarball before delete). Log: `/home/rander39/resume_delete_p018csarc_20260906.log` (ends with `DONE <timestamp>` + final `du -sm /home/rander39`; was mid-run-1 at handoff). Disk was 421G/400G cap before it. The 6 restore-marker STALE runs (`fp052d5*_xverify*`, `p026_restart_*`) are intentional residue — leave alone. NOTE: the hpc-storage agent refuses relayed authorization for irreversible deletes (correctly) — that's why the worker was launched from the main session on Ryan's direct in-session approval.

**Still pending (surface, don't act without Ryan):**

- The five `scr_p026ef_*` run dirs: Ryan authorized archiving 2026-09-05, but the archiver skips them as RECENT until quiet >24 h. Re-run hpc-storage once they age out (~1 day). Warm-start bracket VTPs (ctrl 224/226, exp 209/211, ctrl-lg 229/231, exp-lg 284/286) sit OUTSIDE the newest-5 retention — ask hpc-storage to keep those steps on /home, or extract from tarball later.
- Notebook entry for Phase 0/1/1b + §16–§18 (deferred 4×; ask verbosity when offering).
- Local FLOWPanel repo uncommitted state: dispatcher (ef cases + linegauss default + rk3 case), design-doc §17+§18, this handoff + prompt_1. Standing decision: merge RK3 wiring + dispatcher cases to mainline (cluster branch `p026-phase1b-rk3` has them; local does not). A campaign worktree + pinned tags are REQUIRED before any implementation commits turn into launches (no-silos policy).
- Protect-listing `p026_restart_*`/`scr_p026ph1*`; `~/wt026` retention.

**Gotchas (inherited, all live):** `ssh orc` needs the ControlMaster socket (2FA otherwise; ask Ryan for `! ssh orc echo ok`). Non-interactive ssh: full slurm paths. sacct state unreliable both ways — judge by outputs. Step lines in .out begin with a TAB. Wake-health CSV early rows contain NaN in max_dtZ — guard awk with `$9+0`. Loads = monitor02 **CFx** (thrust axis x; never average CFz). Do not score vs pre-08-24 verdicts. Never tick notebook checkboxes; notebook writes need Ryan's approval.
