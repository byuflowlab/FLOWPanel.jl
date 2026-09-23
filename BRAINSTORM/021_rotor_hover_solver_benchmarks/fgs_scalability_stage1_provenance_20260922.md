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
- [x] Tagged triple rsync-deployed 2026-09-22 to
      `/home/rander39/campaigns/p021-fgs-stage1-20260922/{FLOWPanel.jl,FastMultipole,FLOWVPM.jl}`
      (fresh dirs, did not exist before; content shipped via
      `git archive <tag> | ssh orc tar -x`, so exactly the tagged tracked
      content — no dirty working-tree state). Per-repo sha256sum manifests
      generated from the same local export, shipped to
      `<deploy>/MANIFEST.<name>.sha256`, and **verified on orc**:
      manifest self-hashes match both sides and `sha256sum --quiet -c`
      passes in all three trees. R4 mesh family
      `examples/data/dji9443_20260813_*_capped_captess4.msh` confirmed
      present in the deployed FLOWPanel tree.
      | manifest | files | sha256 |
      |---|---|---|
      | MANIFEST.FLOWPanel.jl.sha256 | 3190 | `dfcf84e3c0033d5f3cf4504d3956951922c4af70c307ab1005027ec530b17aa5` |
      | MANIFEST.FastMultipole.sha256 | 6478 | `b9aa6449a40fbe2d64ebf467e25efa8ba981e340859d184c002801c0bfc234e4` |
      | MANIFEST.FLOWVPM.jl.sha256 | 188 | `b4d0fa8d4a66cc76483f8c29c7261593bad7daad67d331d3af37f222215ebd77` |
- [x] Campaign env (COLD_PROJECT =
      `/home/rander39/campaigns/p021-fgs-stage1-20260922/env`) built with
      julia/1.11.7-6bmogfl: Project.toml copied from the
      p021-r12-champion-20260919 env, `Pkg.develop` on the three deployed
      trees + `Pkg.instantiate()`; Manifest dev-paths confirmed to resolve
      to the deploy dirs (satisfies the `cold_packages` realpath gate).
- [x] CAMPAIGN_PINS at
      `/home/rander39/campaigns/p021-fgs-stage1-20260922/pins.toml` with
      `deployment = "rsync"` + `content_manifest(_sha256)` per package
      (schema consumed by `cold_packages`,
      `benchmark/fgs_cold_common.jl:276`); tag/sha per the pins table above.
- [x] zen3 availability confirmed via slurm-availability (`--cpus 128
      --mem-gb 500 --time 48:00:00 --eta`, 2026-09-23T01:44Z): m12
      `access=normal`, 5 nodes fit and idle, estimated start
      2026-09-23T03:27Z under `--qos=normal`.
- [x] **Submitted 2026-09-22 (job 13858983, m12, PENDING at submission)**:
      `sbatch --export=ALL,COLD_PROJECT=.../env,CAMPAIGN_PINS=.../pins.toml,COLD_DATA_ROOT=/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910 benchmark/run_r4_fgs_stage1.slurm.sh`
      from the deployed FLOWPanel tree top level
      (`/home/rander39/campaigns/p021-fgs-stage1-20260922/FLOWPanel.jl`,
      `logs/slurm/` pre-created). CALIB_JOB default 13777133; the four
      calibrated `dagteam_selected.toml` files verified present 2026-09-22.
      Run dir will be
      `data/p021-cold-20260910/fgs-stage1-13858983`.

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
