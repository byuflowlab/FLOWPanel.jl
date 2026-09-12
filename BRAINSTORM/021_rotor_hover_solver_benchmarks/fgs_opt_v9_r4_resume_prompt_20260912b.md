# BRAINSTORM 021: v9 R4 diagnostics resume prompt B (2026-09-12, second reset)

Entry point for the next agent. Supersedes
`fgs_opt_v9_r4_resume_prompt_20260912.md` (its steps 1–3 are DONE; its
"Session gotchas" section remains valid — read it). The governing handoff is
still `fgs_opt_r4_diagnostics_handoff_20260912.md` (binding: stopping
boundary, gates, diagnostics sequence). Pins and both job launches are in
`fgs_opt_v9_provenance_20260912.md`. Read `~/.claude/CLAUDE.md`, repo
`CLAUDE.md`, `agent_policies/HPC.md`; on ORC read
`/apps/instructions_for_ai_agents/BYU_ORC_AGENTS.md` before cluster ops.

## Mission (unchanged)

Finish the bounded near/far tuning (stage B), run finalist confirmation and
reprofile, and hand Ryan validated tables, text profiles, and a ranked list of
>=5% implementation opportunities (measured vs hypothesis distinguished).
**Do not implement solver optimizations. No notebook entry without Ryan's
approval.** Config-only reruns reuse clean v9; any harness code change needs a
new pinned generation (v10) — never move v9 tags or edit a running worktree.

## State at reset

**Job 13657404 (R4 all-stage) COMPLETE and harvested.** All gates PASS.
Durable evidence:
`BRAINSTORM/021_rotor_hover_solver_benchmarks/fgs_opt_evidence_20260912/opt-13657404/`
(start with `harvest_summary.md`, `screen_rank.csv`). Key results, j64/b1
prepared medians (58,192 panels, 10 reps):

- Inner screen: inner=3 10.709 s (27 outer iters) < inner=5 11.005 s (17) <
  inner=2 11.383 s < inner=10 12.650 s < inner=1 13.049 s. Median and
  min-time (selected.toml) rankings agree: winner inner=3, runner-up inner=5.
- Baselines (seed P8/MAC0.4/leaf100/inner3): j4/b1 17.016 s, j64/b1 11.099 s.
- R4 CPU profile (inner=3 seed, 6717 snapshots, read independently of R2):
  `compute_nonself_products!` → BLAS `dgemv` **~77%** of snapshots
  (`dgemv_kernel_4x4` 5170 self-samples); `scatter_nonself_influence!` ~6%;
  `solve_leaf!` ~5%; all `fmm!` passes ~1.3%; `daxpy` ~2.4%. Top allocation
  site: `influence!` array temporaries (FastMultipole solve.jl:1336), 1%
  sampling. R2's dense-nonself-GEMV attribution is CONFIRMED and sharper.
- Provenance note: `harvest_summary.md` was harvested by a subagent and then
  CORRECTED by the main agent (it had named inner=2 "median winner" — wrong;
  mis-mapped calibration tolerances; missed `baseline-*/`+`smoke-*/` dirs,
  since pulled). Trust the corrected file and `screen_rank.csv`.

**Job 13660643 (step-7 stage A: leaf axis) RUNNING at reset.** Submitted
2026-09-12 ~13:25, 8 h walltime, node m12-1-25 (same as job 1 — matched).
Output: `/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/opt-13660643/`.
Launch env recorded in provenance file: `COLD_OPT_STAGE=screen_profile`,
`CONFIG_FILE=.../opt-13657404/smoke-j4-b1/selected.toml`,
`COLD_OPT_SCREEN_BASE_FILE=/home/rander39/campaigns/p021-cold-opt-20260912-v9/bases_r4_leaf_20260912.toml`
(md5 `be538fc013863860e9a181f538076649`; two bases = inner=3 and inner=5
winners with calibrated tolerances), `COLD_OPT_SCREEN_SET=leaf:25,50,200`,
`COLD_MIN_REPS=10`. Roster = 8 candidates (2 bases + 6 leaf neighbors;
leaf=100 covered by the bases; harness resets tolerance to 0 and recalibrates
each point). Progress at reset (~2 h in): baselines done, screen 2/8 status
files, third candidate in flight (~20 min/candidate pace); then a stock
~35 min profile of the CONFIG_FILE seed. The old background poller died with
the context reset — start your own (plain Bash ssh loop on the `COMPLETED`
marker + sacct terminal state, >=60 s between scheduler queries; judge by
outputs, not sacct).

## Ryan's decisions this session (binding)

1. **Stage B runs as a 2-way parallel split** (approved): one job per
   retained stage-A base, each node self-anchored (its own base rerun + stock
   baselines) so all ranking comparisons stay within-node. Node-to-node
   variance ~ effect size (~2–3%), so never rank candidates across nodes
   without a shared anchor.
2. More parallel HPC submissions are welcome where they don't compromise
   comparisons: confirmation as 2 jobs on separate nodes (independent
   replications), reprofile concurrent with confirmation.

## Next actions

1. Poll 13660643 to completion. Validate gates per handoff (BC rel-L2 <=1e-6,
   certified evaluator or explicit direct fallback, repeat agreement <=1e-8,
   convergence; map configs via config.toml, not dir order; failed candidates
   are findings). Harvest to
   `fgs_opt_evidence_20260912/opt-13660643/` mirroring the opt-13657404
   layout (delegate scraping to `harvester` but AUDIT its tables against the
   raw CSVs — it made real errors last time — and make sure it pulls
   `baseline-*/` result dirs, not just logs). Rank all 8 by median prepared
   time; retain top two by median.
2. **Stage B (2 parallel jobs):** for each retained base, write a one-config
   bases TOML (full config table incl. its newly calibrated tolerance from
   stage A output `config.toml`) to
   `/home/rander39/campaigns/p021-cold-opt-20260912-v9/bases_r4_<slug>_20260912.toml`,
   append launch settings + bases md5 to the provenance file BEFORE sbatch,
   verify worktree still clean at `f03ab18`, then submit each with the job-2
   env template (provenance file) changing only:
   `COLD_OPT_SCREEN_BASE_FILE=<that bases file>` and
   `COLD_OPT_SCREEN_SET="P:6,10;MAC:0.3,0.5"` (5 candidates/job: base + 4
   one-factor neighbors; ~2.5–3 h; 6 h walltime is enough). Same
   CONFIG_FILE (job-1 smoke selected.toml). sbatch needs a login shell:
   `ssh orc 'bash -lc "... sbatch ..."'` (`/apps/slurm/latest/bin`). Use
   `--test-only` first. No full Cartesian sweep; retain two by median
   overall (within-node comparisons + each job's baseline anchor).
3. **Confirmation (2 parallel jobs, separate nodes):** confirm finalists vs
   the seed baseline with >=10 unprofiled trials in alternating batches.
   Stock `cold_run` times a multi-config CONFIG_FILE sequentially in one
   process — so a `confirm.toml` with `configs=[seed, finalist1, finalist2]`
   (each with its calibrated positive tolerance; distinct configs only, dup
   cold_ids collide) makes each baseline process time B→C1→C2 as consecutive
   10-trial batches in one matched environment; two simultaneous jobs on
   different nodes = two independent replications. Use
   `COLD_OPT_STAGE=screen_profile` with `CONFIG_FILE=confirm.toml`; the
   screen stage must run — give it a trivial roster (bases file = finalists,
   `SCREEN_SET` naming one axis at its base value dedupes neighbors into the
   base) which doubles as an extra recalibrated candidate batch. CAVEAT:
   check how `rotor_hover_solver_phase2_profile.jl` handles a multi-config
   CONFIG_FILE before submitting (it profiles CONFIG_FILE at stage end);
   if it errors or profiles only configs[1], that's fine — note it.
   Repeat ambiguous comparisons; never loosen gates to force a winner.
4. **Reprofile the retained configuration** (if it differs from the already
   profiled inner=3/leaf=100 seed) in its own clean allocation, concurrent
   with step 3: `screen_profile` job with single-config
   `CONFIG_FILE=<retained calibrated config>`, trivial screen roster as
   above, `COLD_OPT_PROFILE_REPS=20` to match job-1 profile density. The
   profile stage profiles CONFIG_FILE with `COLD_INVESTIGATION=1`.
   (`attribution` stage can NOT reprofile an arbitrary config — the launcher
   overwrites CONFIG_FILE with the smoke selected.toml in that path.)
5. Harvest stage B + confirmation + reprofile into
   `fgs_opt_evidence_20260912/opt-<jobid>/` dirs; then deliver the final
   diagnostics package: validated config/timing tables, text profiles,
   ranked >=5% implementation opportunities with supporting R4
   stacks/allocations (measured vs hypothesis), extra measurements needed.
   Current measured picture the ranking must reflect: dense nonself GEMV
   ~77%, influence! allocation temporaries, scatter ~6%, leaf solves ~5%,
   FMM ~1.3%. Then STOP. Offer a notebook entry (ask detail level); do not
   write it without approval.

## Cost/scale notes

- Screen pace at R4: ~20 min/candidate (calibration solves dominate, timed
  10 reps are ~2 min); baselines ~14 min each; profile ~35 min at 20 reps.
- Storage: campaign outputs are small (no VTK); local profile text is the
  bulk (~150 MB/harvest). Home disk was ~192 GiB used this morning.
- All ssh via the `orc` alias (needs live ControlMaster; if cold, ask Ryan
  to run `ssh orc -fN`; never retry into MFA).
