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

## Wave 1 FAILED — diagnosis and wave 1b resubmission (2026-09-28)

All 8 jobs (13901206–13) terminated at walltime (13901213 FAILED at ~4 h)
with **zero CSVs**: no user code ever ran. Root cause: all 8 jobs launched
simultaneously and raced to precompile the campaign env (fresh cache slug
`_2g3HB` in `~/.julia/compiled/v1.12/`). Julia's pidfile locks are created
mode 0444; on the VAST NFS home an existing 0444 pidfile surfaces as
`EACCES` (not `EEXIST`), so one early crash left 21 stale read-only
pidfiles (StructUtils → … → FLOWPanel) that wedged every job in
precompilation until walltime. First error (13901213):
`IOError: open(".../StructUtils/vyJDy_2g3HB.ji.pidfile", 194, 292): permission denied (EACCES)`.
Quota/permissions ruled out (823 G/2 T, login-node write test OK).

Remediation (2026-09-28): deleted the 21 stale `*_2g3HB.ji.pidfile`;
added `launch_p033_precompile.sh` (campaign dir, alongside `launch_p033.sh`)
— single m12 job that clears stale pidfiles (>60 min old), runs
`Pkg.precompile()` on the campaign env, and verifies
`using FLOWPanel; using FastMultipole`. Wave 1b resubmitted identically
(same knobs/walltimes/partition/qos as the table above) with
`--dependency=afterok:<precompile>` so no job precompiles concurrently.

| job | wave-1b id | note |
|---|---|---|
| fp033-precompile | 13905250 | 2 h cap, gates all arms (afterok) |
| fp033-ar3-j1 | 13905251 | = 13901206 knobs |
| fp033-ar3-j8 | 13905252 | = 13901207 |
| fp033-ar3-j16 | 13905253 | = 13901208 |
| fp033-ar3-j32 | 13905254 | = 13901209 |
| fp033-ar3-j64 | 13905255 | = 13901210 |
| fp033-bi1-r4j64 | 13905256 | = 13901211 |
| fp033-bi1-r5j64 | 13905257 | = 13901212 |
| fp033-bi1-r4j1 | 13905258 | = 13901213 |

Output dirs get new `<jobname>-<jobid>` names; the eight empty wave-1 dirs
(`*-139012??`) can be deleted at harvest time. Lesson for future waves:
never let N simultaneous jobs cold-precompile a fresh env on this
filesystem — always gate the wave on one precompile job.

## Wave 1b precompile FAILED — julia version mismatch; wave 1c (2026-09-28)

Precompile 13905250 FAILED at 3.5 min with the real underlying failure of
wave 1 exposed: compute nodes' default `julia` is now **1.12.7** (juliaup;
`command -v julia` finds it, so the launcher's 1.11.7 spack fallback never
fires), but the campaign Manifest was resolved with **1.11.7**. Under
1.12.7, HDF5_jll's `libhdf5.so` needs symbol version `CURL_4`, which
1.12.7's bundled `libcurl.so.4` lacks → `InitError` → HDF5 → FLOWPanel
precompile fails. So wave 1 was doubly doomed: the pidfile race wedged it,
and even without the race every arm would have died here.

Remediation: both launchers now export the spack julia-1.11.7 bin dir at
the FRONT of PATH unconditionally (no `command -v` guard). Cancelled the
held 1b arms (13905251–58; 13905250 already failed) and resubmitted as
wave 1c, same knobs/dependency structure:

| job | wave-1c id |
|---|---|
| fp033-precompile | 13905268 |
| fp033-ar3-j1 | 13905269 |
| fp033-ar3-j8 | 13905270 |
| fp033-ar3-j16 | 13905271 |
| fp033-ar3-j32 | 13905272 |
| fp033-ar3-j64 | 13905273 |
| fp033-bi1-r4j64 | 13905274 |
| fp033-bi1-r5j64 | 13905275 |
| fp033-bi1-r4j1 | 13905276 |

Lesson: pin the julia binary explicitly in every campaign launcher —
juliaup's default channel drifts (1.11.7 → 1.12.7 sometime before
2026-09-27) and `command -v julia` fallbacks silently pick up the drift.

## Wave 1c results + B-I1 partial rerun (2026-09-28)

All 8 wave-1c arms COMPLETED with outputs (8–67 min each; the original
4–8 h walltimes were sized for the precompile wedge, not real work).

- **A-R3 landed clean**: full j-ladder harvested to
  `BRAINSTORM/033_ar3_20260928/`; results in `033_ar3_results_20260928.md`;
  A-R3 + A-G ticked (worker + clear-context review). A-G verdict: j64 w=4
  paired uninstrumented gain −18.85% vs ≥+15% required ⇒ **Track A
  implementation stops**.
- **B-I1 partially degraded by a submission bug**: the wave-1c resubmission
  passed knobs via `--export=ALL,ARMS=old,new,...`; sbatch splits
  `--export` on commas, so every setup job got `ARMS=old` (the `new` token
  was dropped). The pass-loop ran old-arm only; at j64 the new-setup
  numbers still landed via the `CERT_SOLVE=1` full-ctor path (1 pass:
  R4 30.0 s, R5 83.5 s, matrices + solves certified identical), but the
  j1 sanity arm (CERT_SOLVE=0) produced no new-arm data at all.
  **Gotcha for all future submissions: never pass comma-valued env via
  `--export` lists — use `VAR=x sbatch` env inheritance.**
- B-I1 rerun with correct env passing (no gate needed, cache warm):
  fp033-bi1-r4j64b **13905507**, fp033-bi1-r5j64b **13905508**,
  fp033-bi1-r4j1b **13905509** (same knobs as the table above).
- Banner-vs-pin SHA note: job banners print worktree HEADs
  (fp 99aa273, fm e46e91d7, vpm 5983c34), each exactly one commit ahead of
  its `campaign/p033-hpcwave1-20260926` tag commit (fp 6ecd0b3,
  fm b4c35f67, vpm 8d4a3b4) — the extra commit is only the standard
  `data/ -> consolidated data root` symlink commit (verified via
  `git log --stat`; source identical to pins). Provenance intact.
