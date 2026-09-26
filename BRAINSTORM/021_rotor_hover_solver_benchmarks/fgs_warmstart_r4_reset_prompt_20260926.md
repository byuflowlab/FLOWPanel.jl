# Reset prompt: 021 warm-start R4 — babysit + harvest ilu_nfcache_proj2 (2026-09-26; supersedes fgs_warmstart_r4_reset_prompt_20260925b.md)

Copy everything below the line into a fresh agent session started in
`~/Dropbox/research/projects/FLOWPanel.jl`.

---

## Mission

The 021 warm-start R4 campaign is HARVESTED and its deliverable written
(`fgs_warmstart_r4_results_20260925.md` + backing
`fgs_warmstart_r4_results_20260925/harvest_*`, committed `b6a0924`). Ryan
then directed ONE follow-up arm: **`ilu_nfcache_proj2`** (quadratic
extrapolation warm-start for the ILU family), closing the best-vs-best
asymmetry — the original slate compared fgs_proj2 (order 2) against
ilu_nfcache_proj1 (order 1). Remaining work, in order:

1. **Babysit job 13899021** on orc (submitted 2026-09-26 09:14 MDT,
   started RUNNING immediately on m12-2-11, `--time=08:00:00`; expected
   ~2.5–3 h total). It runs exactly two legs, sequential:
   `ilu_nfcache_proj2:winA` (steps 1–36 from scratch) then
   `ilu_nfcache_proj2:winB` (steps 109–144 restarted from the existing
   `fgs_wsr4_R4_ckpt_ilu` step-108 checkpoint). It APPENDS to the existing
   run dir `/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/fgs-wsr4-13890195/`
   (RESUME_FROM_JOB_ID=13890195; the 14 landed legs skip by STATUS).
   Judge by `STATUS_ilu_nfcache_proj2_win{A,B}` sentinels + CSV
   `solved` column, never sacct.
2. **Re-run the harvest** (same command as before):
   `cd /home/rander39/campaigns/p021-fgs-warmstart-20260924/FLOWPanel.jl &&
   julia --project benchmark/fgs_r4_warmstart_harvest.jl
   /home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/fgs-wsr4-13890195`
   — it is arm-agnostic (keys on config/warmstart/order) and will pick up
   the new rows and snapshots automatically.
3. **Update the results doc** `fgs_warmstart_r4_results_20260925.md`:
   add the ilu_nfcache_proj2 row to the Window A and Window B tables, the
   cumulative table, the agreement table, and re-run the crossover math
   (fgs_proj2 vs ilu_nfcache_proj2 is now THE head-to-head; previous
   estimate vs proj1 was ~480 steps). Refresh the backing copies in
   `fgs_warmstart_r4_results_20260925/`. Update the notebook 20260926
   entry's tables the same way (day header exists; append a correction
   note, never rewrite — see notebook policy). Commit.
4. Surface (don't act): compete/no-compete decision is Ryan's; origin
   pushes of ALL campaign tags (now incl.
   `campaign/p021-fgs-warmstart-20260926`) owed after `gh auth login`;
   Ryan may still veto the submission defaults (champion tol
   3.4309419310610173e-7, proj1 kept, sequential legs, exclusive node).

Required reads first: `CLAUDE.md`, `agent_policies/HPC.md`, cluster
`BYU_ORC_AGENTS.md`, `fgs_warmstart_r4_provenance_20260924.md`
(**including its 2026-09-26 addendum** — the follow-up's full provenance).
`ssh orc` needs a live ControlMaster socket (`ssh orc -fN` if cold; never
retry into 2FA). Local runs ≤4 threads. Slurm needs a LOGIN shell:
`ssh orc 'bash -lc "…"'`.

## State (2026-09-26, ~09:20 MDT)

- **Job 13899021 RUNNING** on m12-2-11 (2 legs, 8 h walltime). Slurm log:
  `$DEPLOY/FLOWPanel.jl/logs/slurm/r4-fgs-wsr4-13899021.out`,
  `DEPLOY=/home/rander39/campaigns/p021-fgs-warmstart-20260924`.
- **Pins**: FLOWPanel `b6a09244` = tag `campaign/p021-fgs-warmstart-20260926`
  (adds the one ARM_TABLE line); FastMultipole `745af760` + FLOWVPM
  `8d4a3b4d` unchanged (tag `…-20260924`). Deploy tree refreshed
  minimally (one file + manifest line), whole tree `sha256sum -c`
  re-verified clean, `pins.toml` updated. R1 smoke of the new arm PASS
  (3/3 legs; winB niter 7→4→3→3, extrap order 2).
- **Prior results (for comparison when updating the doc)** — winB medians:
  fgs_proj2 8.245 s / 2 iters; ilu_nfcache_proj1 8.680 s / 5 iters;
  ilu_nfcache_prev 9.244 s / 7. WinA medians: fgs_proj2 10.68;
  ilu_nfcache_proj1 9.90. Setup: FGS 315–320 s, ILU ~110 s (incl. prime).
  Expected: ilu proj2 gains ~0.2–0.5 s/step over its proj1.
- **Storage (concurrent, both may still be running)**: two detached
  archiver `--apply` workers launched 2026-09-26 on orc login node —
  logs `/home/rander39/archiver_worker_{main,ntladder}_20260926.log`.
  Scope: `fgs_wsr4_R4_ckpt_fgs` (main checkout) + three finished
  p018-ntladder arms (~325 G would-free; /home was at 624 G of the 400 G
  cap). MUST collect both workers' `STALE_COUNT=`, `VERIFY_FAIL_COUNT=`
  and exit status before declaring the cycle done. RECENT items
  deliberately NOT archived: `fgs_wsr4_R4_ckpt_ilu` (13899021 restarts
  from it — leave until that job is done) and one p018 ntladder arm
  (Ryan's call). `p021-cold-20260910` reads NO-VTK to the archiver —
  correct, benchmark run dirs hold only CSVs/snapshots (policy: stay on
  /home); their VTK is the separate ckpt trees.

## Traps

- Leg logs are OUTPUT-BUFFERED — judge liveness by sstat CPU, STATUS
  sentinels, CSV rows; never by log tails. "fatal: not a git repository"
  noise is expected (archive-export deploy; CSV `commit` = "unknown").
- sacct FAILED ≠ real failure, COMPLETED ≠ success — judge by outputs.
- The harvest overwrites `harvest_{summary.md,traces.csv,solution_deltas.csv}`
  in the run dir — fine, they're committed at the pre-proj2 state in
  `fgs_warmstart_r4_results_20260925/`.
- Restarted (winB) rows log CT = −0 (monitor not reattached on restart) —
  known, reported in the results doc; don't chase.
- Window-B stats INCLUDE the (order+1)-step history-refill transient by
  design (Ryan 2026-09-24) — for proj2 that's 3 effectively-cold steps.
- If 13899021 dies mid-sequence: resubmit the same submit line (it's in
  the provenance addendum) — STATUS-landed legs skip automatically.
- 018 NT-ladder jobs run concurrently — disk alarms are theirs; launch
  hpc-storage, never touch their queue.
- Do NOT edit the live deploy tree while 13899021 runs.
