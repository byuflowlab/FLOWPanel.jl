# 021 FGS scalability — Stage 1 provenance (2026-09-22)

Plan: `fgs_scalability_diagnostic_plan_20260921c.md` (plan C). Stage 0:
`fgs_scalability_stage0_audit_20260922.md`. Handoff executed:
`fgs_scalability_stage1_reset_prompt_20260922b.md`. HPC submission
pre-approved by Ryan for Stage 1 ONLY.

## Campaign pins (annotated tag `campaign/p021-fgs-stage1-20260922` in all three repos)

| repo | branch | pinned commit | tag verified annotated |
|---|---|---|---|
| FLOWPanel.jl | `fastmultipole` | `c8621af789eccbe4b8234da9cd4339ade50fcb77` | yes |
| FastMultipole | `flowpanel-20260817` | `649405cba0db01b82fbf8336e8795be1d01ef1fa` | yes |
| FLOWVPM.jl | `flowpanel` | `8d4a3b4d3012c42fc7d078629234c105b1e570f7` | yes |

- FLOWPanel `c8621af` = Stage-1 commit on top of `4de873c` (driver
  `benchmark/fgs_r4_dagteam_stage1.jl`, launcher
  `benchmark/run_r4_fgs_stage1.slurm.sh`, analysis
  `benchmark/fgs_stage1_analysis.jl`, `dagteam_workers` plumbing in
  `FGSSolver`/metadata/`fgs_cold_common.jl`, plan C + reset prompts).
- FastMultipole `649405cb` = Stage-1 commit on top of `f4d6b671`
  (`build_dagteam_plan` `nworkers` kwarg, `FastGaussSeidel`
  `dagteam_workers`, coordinator coarse timers
  `:dagteam_{spawn,join,wait,reduce}_ns` — all subsets of
  `:nonself_product_ns`).
- FLOWVPM `8d4a3b4` = loaded commit of the live checkout (new merge law,
  production default per Ryan 2026-09-19); FLOWVPM is untouched by Stage 1
  and pinned only to complete the triple.

## Local verification (2026-09-22, laptop, ≤4 threads)

- `test/runtests_unit_solver.jl`: all testsets pass (exit 0).
- 4-thread dagteam worker-cap smoke (scratchpad script): 37/37 pass —
  caps 0/1/2/4/99 build team size `clamp(cap,1,nthreads)` (`plan.xg/yb`
  length checked), solutions match the lexicographic reference to 1e-6 and
  are **bit-identical across team sizes**; diagnostics dict carries the 4
  new keys, their sum ≤ `:nonself_product_ns` when :dagteam, all zero for
  :lexicographic; solver metadata exposes `dagteam_workers`.
- `cold_check_config` whitelist smoke: `dagteam_workers` accepted on
  dagteam (0/16), rejected for negative, Bool, and non-dagteam configs.
- `benchmark/fgs_stage1_analysis.jl`: parse-checked (`Meta.parseall`);
  first execution will be against the harvested Stage-1 CSVs.

## Calibrated inputs (verified on orc 2026-09-22)

`dagteam_selected.toml` present for j ∈ {1,16,32,64} under
`/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/thread-scaling-j<J>-13777133/fgs-calibrate/results/`
(202 B each, mtimes 2026-09-19). CALIB_JOB=13777133.

## Deployment (fill at submission)

- **RULING (Ryan 2026-09-22)**: GitHub creds dead on this machine, so the
  origin push is deferred (still owed once Ryan re-auths). Deployment
  switches to the harness's `deployment = "rsync"` mode (content-manifest
  verified at runtime by `cold_packages`), shipping the TAGGED content into
  fresh dirs under `/home/rander39/campaigns/p021-fgs-stage1-20260922/`.
  See `fgs_scalability_stage1_reset_prompt_20260922c.md` TASK 1 for the
  full mechanics.
- [ ] Branches + tags pushed to origin (deferred; owed after re-auth)
- [ ] Tagged triple rsync-deployed to
      `/home/rander39/campaigns/p021-fgs-stage1-20260922/` (fresh dirs;
      manifest sha256s recorded here when created)
- [ ] Campaign env (COLD_PROJECT) with Manifest dev-paths at the deployed
      trees
- [ ] CAMPAIGN_PINS toml with `deployment = "rsync"` +
      `content_manifest(_sha256)` per package (schema consumed by
      `cold_packages`, `benchmark/fgs_cold_common.jl:272`)
- [ ] zen3 availability confirmed via slurm-availability (`--cpus 64
      --mem-gb 500 --eta`) immediately before submission
- [ ] `sbatch --export` submission of `benchmark/run_r4_fgs_stage1.slurm.sh`
      from the FLOWPanel campaign worktree top level with
      COLD_PROJECT / CAMPAIGN_PINS /
      COLD_DATA_ROOT=`/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910`
      (job id recorded here)

## Measurement caveats (carried from the reset prompt, binding on analysis)

- Fixed-work rows have `solved=false`/`eligible=false` BY CONSTRUCTION
  (tolerance=0); gates are certified-accepted + iterations==27 + 1e-8
  repeat. Never filter Stage-1 fixed rows on eligible/solved.
- Accepted solves run ONE extra fmm+influence+residual vs fixed work.
- The worker cap does NOT pin active workers to specific cores (tasks float
  inside the 64-core cpuset); the plan's "matching nested CPU set" is
  approximated by the champion cpuset.
- `diag_*` = −1 marks uninstrumented rows; never pool instrumented and
  uninstrumented rows.
- Plateau (16→32) and regression (32→64) are SEPARATE conclusions, each
  gated on reproduction; no interventions on a non-reproduced effect.
