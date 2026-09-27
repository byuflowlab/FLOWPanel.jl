# 033 HPC wave 1 — provenance (submitted 2026-09-26, orc time 2026-09-27 UTC)

Campaign: **A-R3** (cooperative-executor paired j-ladder) + **B-I1** (threaded
setup j64/R5 re-measure). Authorized by Ryan 2026-09-26 ("launch the hpc runs").

## Pins (annotated tag `campaign/p033-hpcwave1-20260926` in all three repos)

| repo | tag commit | worktree HEAD (tag + data-symlink commit) |
|---|---|---|
| FLOWPanel.jl | `6ecd0b3` (033 session 2) | `99aa273` |
| FastMultipole | `b4c35f67` (033 A-R2+B-R2, on `flowpanel-20260817`) | `e46e91d7` |
| FLOWVPM.jl | `8d4a3b4` (unchanged dependency pin) | `5983c34` |

Tags pushed to the `orc` remote (cluster clones). GitHub `origin` push FAILED
locally (https credentials unavailable in the session shell) — **tags are NOT
on origin yet**; push from an interactive shell when convenient.

Pre-commit gate: FLOWPanel solver tests 513/513 PASS; FastMultipole full suite
PASS (+ dagteam/dagedge gates) at default toggles, 2026-09-26 local.

## Cluster layout

- Worktrees: `~/campaigns/p033-hpcwave1-20260926/{FLOWPanel.jl,FastMultipole,FLOWVPM.jl}`
  via `scripts/prep_campaign_worktree.sh` (extracted from the tag; not present
  on the orc clones' `unified-052` branch). `data/` symlinked to
  `~/projects/FLOWPanel.jl/data`. No uncommitted state.
- Env: `~/campaigns/p033-hpcwave1-20260926/env/` — copied from
  `p021-fgs-warmstart-20260924/env`, Manifest dev-paths rewritten to the p033
  worktrees (3 substitutions verified); sha256 in `MANIFEST.env.sha256`;
  `pins.toml` alongside. `ARCHIVER_SKIP` set.
- Launcher: `~/campaigns/p033-hpcwave1-20260926/launch_p033.sh` (generic:
  `HARNESS=coop|setup`, `J`, harness knobs via `--export=ALL`). BLAS pinned
  single-threaded at load (`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1` +
  `BENCH_BLAS_THREADS=1`); `THREADING_MODE=multi EXPECT_JULIA_THREADS=$J`.
- Outputs: `~/projects/FLOWPanel.jl/data/p033_hpcwave1_20260926/<jobname>-<jobid>/`
  (outside the worktree). Slurm logs: worktree `logs/slurm/`.

## Jobs (all m12 zen3, `--qos=normal`, 64 tasks, 96 G; submitted 2026-09-27 UTC)

| job | id | harness / knobs | walltime |
|---|---|---|---|
| fp033-ar3-j1 | 13901206 | coop MODE=time RUNG=R4 ROUNDS=5 WMAX=1 CACHE_B=1, J=1 | 6 h |
| fp033-ar3-j8 | 13901207 | coop … WMAX=4, J=8 | 4 h |
| fp033-ar3-j16 | 13901208 | coop … WMAX=4, J=16 | 4 h |
| fp033-ar3-j32 | 13901209 | coop … WMAX=4, J=32 | 4 h |
| fp033-ar3-j64 | 13901210 | coop … WMAX=4, J=64 | 4 h |
| fp033-bi1-r4j64 | 13901211 | setup RUNG=R4 ARMS=old,new AB_K=3 CERT_SOLVE=1, J=64 | 6 h |
| fp033-bi1-r5j64 | 13901212 | setup RUNG=R5 ARMS=old,new AB_K=3 CERT_SOLVE=1, J=64 | 8 h |
| fp033-bi1-r4j1 | 13901213 | setup RUNG=R4 ARMS=old,new AB_K=2 CERT_SOLVE=0, J=1 | 6 h |

At submission +2 min: 7 RUNNING, fp033-bi1-r4j1 PENDING.

## Analysis contracts

- A-R3 judge: paired same-process w∈{1,2,4} per j from each job's CSVs;
  A-G stop gate = team-of-4 gain <~15% paired wall clock at j64
  (uninstrumented paired data). Also harvest in-situ h_w vs A-R1/A-R2 locals.
- B-I1 judge: old vs new setup at R4/R5 j64; feed the measured new setup into
  the Tier-2 cumulative-crossover arithmetic vs krylov_ilu_nfcache (~110 s
  setup, warm-start R4 doc). R5 = the "larger mesh" owed from B-T1's spec
  deviation.
- Slurm exit status is unreliable (standing ruling) — judge from output CSVs.
