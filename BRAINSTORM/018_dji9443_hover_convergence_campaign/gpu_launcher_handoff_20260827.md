# 018 GPU launcher — handoff (2026-08-27)

Status: **GPU launcher built and smoke-validated on GH200 and H200.** Written by
the agent that built it, for the next agent. Ryan has seen the smoke results.

## What exists

**Launcher** (in this repo, uncommitted):
- `examples/run_dji9443_hover_ct_gpu.slurm.sh` — GPU sbatch payload. Args
  `<arch> <case_tag>` (arch ∈ gh200|h200|h100; gh200 default matches the
  in-file `#SBATCH` header). Wraps the CPU dispatcher so the ~103 `p018_*`
  case blocks stay single-sourced. Sets the 052 GPU env (`VPM_ARRAYTYPE=cuarray`,
  `FLOWPANEL_GPU_INFLUENCE=cuda`, `RHPC_SOLVER_S=true`, `RHPC_SOLVER_S_GPU=true`,
  `FASTMULTIPOLE_FORCE_CUDA_LOAD=1`; sources `FM052_GPU_ENV` from the FLOWVPM
  silo when present). Runs julia under an internal timeout (wall − 10 min) so
  the post-run gate always executes; timeout exit 124 = acceptable partial run
  (probe semantics). Gate: `source_influence_s_gpu_gemv` > 0 across
  .out+.err (the @info markers land on **stderr**), CPU `source_influence_s_gemv`
  == 0, no `NaN` outside the CT table (CT-KJ columns are NaN by design when
  `RUN_KJ=false`).
- `examples/run_dji9443_hover_ct_hpc.slurm.sh` — 3 backward-compatible
  parameterizations: `THREADS="${P018_THREADS:-64}"`,
  `EXPECTED_REPO="${P018_REPO:-...}"`, final invocation
  `${P018_JULIA:-julia} --project="${P018_PROJECT:-.}"`; the x86 spack-julia
  PATH fallback is skipped when `P018_JULIA` is set. Unset ⇒ CPU behavior
  identical (running Phase-17 CPU jobs unaffected).

**Cluster state (orc)** — all deployed and verified by md5:
- Silo trees `~/{FLOWPanel,FLOWVPM,FastMultipole}-018-gpu-{gh200,h200}` —
  copies of the 052 silos with `src/` **rsynced byte-identical to the local Mac
  working repos** (FLOWPanel main incl. commit 3e5aa08 device-wake fixes +
  uncommitted 052c work; FLOWVPM d871d47 incl. `sigma_guard`; FastMultipole
  branch flowpanel-20260817). The 052 trees/envs are untouched (052c in flight).
- Envs `~/p018env-{gh200,h200}` — copies of `fm052env-<arch>` with the three
  dev paths repointed to the 018 silos (Manifest lines ~456-482).
- GH200 arch facts: partition `mgh` **requires `-C arm`**, no module loads
  (ARM node has no x86 module tree), julia = `~/julia/julia-1.11.7/bin/julia`
  (aarch64), `JULIA_DEPOT_PATH=~/fm052depot-gh200`. H200 via
  `--partition=eng --qos=eng --constraint=hopper` (must override the in-file
  arm constraint; empty `--constraint=""` is rejected). `test` QOS caps at 1 h.
- Both launchers also copied to `orc:~/projects/FLOWPanel.jl/examples/`.

**Submission recipe** (validated):
```bash
# GH200 (started immediately both times; non-preemptible)
cd ~/FLOWPanel-018-gpu-gh200 && sbatch --job-name=fp-018gpu-<tag> \
  --partition=mgh --gres=gpu:gh200:1 --constraint=arm --qos=normal --no-requeue \
  --cpus-per-task=72 --mem=192G --time=02:00:00 \
  --output=logs/slurm/slurm-%x-%j.out --error=logs/slurm/slurm-%x-%j.err \
  --export=ALL,P018_RUN_NAME=<run_name> \
  examples/run_dji9443_hover_ct_gpu.slurm.sh gh200 <p018_case>
# H200 (eng; preempts standby jobs — courtesy consideration)
cd ~/FLOWPanel-018-gpu-h200 && sbatch ... --partition=eng --qos=eng \
  --gres=gpu:h200:1 --constraint=hopper --cpus-per-task=64 ... h200 <p018_case>
```

## Smoke results (jobs 13501825 gh200 / 13501826 h200, case p018_csarc_l3p4, cold, 719 steps ≈ 20 revs)

| metric | GH200 | H200 | CPU ref (64 cores, tuned) |
|---|---|---|---|
| time-marching wall | 3547 s | 3883 s | — |
| mean s/step | 4.90 | 5.36 | 94.9 |
| s/step last 100 (mature, ~last revs) | 5.96 | 6.45 | 94.9 |
| gpu gemv / cpu gemv | 720 / 0 | 720 / 0 | — |
| NaN (outside CT table) | 0 | 0 | — |

CT̄ per rev vs the CPU `p018_csarc_l3p4` reference, same revs (sanity band —
band of ±0.36 %, comparable to the campaign's arm-to-arm scatter):

| rev | GH200 | H200 | CPU | gh−cpu | h−cpu |
|---|---|---|---|---|---|
| 15 | 0.070858 | 0.070503 | 0.070736 | +0.17 % | −0.33 % |
| 18 | 0.070906 | 0.071050 | 0.071163 | −0.36 % | −0.16 % |
| 20 | 0.070913 | 0.071138 | 0.070895 | +0.02 % | +0.34 % |

**≈ 16–19× speedup at the campaign operating point.** A 30-rev NT36 arm
(1080 steps) extrapolates to ~2 h wall vs the current 48 h.

Harvested artifacts: `data/p018_gpusmoke_20260827/` (CT CSVs, case metadata,
gz'd slurm logs, CPU-ref CT CSV).

Slurm marks both smoke jobs FAILED — cosmetic: their spooled (pre-fix) gate
grepped stdout only. The deployed launcher has the fixed gate (.out+.err, CT
table excluded); it has **not yet been proven in-job**.

## Debug history (why 3 attempts)

1. 13501799/13501804: `module load cuda` fails on ARM mgh → removed for gh200.
2. 13501811 + h200 twin: driver refuses `RHPC_SOLVER_S_GPU` without
   `RHPC_SOLVER_S=true` (driver:1050); then MethodError `_euler(...;
   sigma_guard)` — FLOWPanel main is ahead of the 052 FLOWVPM tree → built the
   FLOWVPM/FastMultipole silos + repointed envs.
3. 13501825/13501826: clean run.

Notes: gemv only starts firing mid-run (0 at step 261, 720 by end — begins when
the Dirichlet/source-mode path activates); `NREVS` cannot shorten a run (the
rev schedule wins) — bound smokes by wall/timeout instead; sacct is
authoritative; strip ANSI from orc output.

## Open items for the next agent

1. **Campaign migration** (needs Ryan's direction): which arms move to GPU;
   walltimes (NT36 30-rev fits ~3 h on GH200); `_s2` restart chaining through
   the GPU launcher (P018_RUN_NAME mechanics unchanged, warm-start VTK loading
   has a GPU staging path in FLOWPanel_warmstart.jl); VTK sweeper paths for the
   silo trees; only 2 GH200 nodes exist — H200/eng is the overflow (preempts
   others' standby jobs) and m13h/qos=gpu the deep queue.
2. **Silo refresh discipline**: silos are snapshots of the local repos as of
   2026-08-27 ~00:00; resync `src/` + driver + launchers (rsync from local)
   whenever local moves, and re-run a 1-case smoke. Depot precompile reruns
   automatically in-job after src changes (first step slower).
3. **Gate**: prove the fixed gate on the next job; consider also gating on the
   FLOWVPM wake path (no per-step marker exists for the device wake — the
   `device-resident CuArray` banner + step time are the current evidence).
4. **052c interaction**: smoke ran unguarded (no `SIGMA_DTZ_CAP`); long
   (1080-step) runs should adopt the 052c sigma guard once it lands, and silos
   should be resynced afterward.
5. Ledger/notebook entries for the campaign are NOT written (needs Ryan's
   approval; draft from this doc + `data/p018_gpusmoke_20260827/`).
6. **NT144 cliff root cause + fix track (2026-08-27):** the cliff is rVPM
   sigma growth of a single runaway particle crossing the radix-FMM ell=3
   adequacy limit (0.0389) → rebuild at ell=2 → quasi-dense near field. Real
   fix: particle splitting — designed as standalone item
   `../026_sigma_growth_particle_splitting/particle_splitting_design.md`
   (extends the existing `FLOWVPM_splitting.jl`, growth-side triggers + two
   split geometries). Interim band-aid: absolute σ ceiling in the rVPM update
   (env `P018`-style knob), deployed to unblock the NT144 arms.
7. **Metadata timing (Ryan, 2026-08-27):** the run `.metadata.toml` is only
   written at end-of-run, so cancelled/killed/running jobs leave no settings
   record (hampered the NT144-cliff diagnosis — radix settings were
   unrecoverable from the cancelled runs). Fix after the sigma-growth bug:
   write the settings/config portion of the metadata at run START, appending
   or updating the completion fields (wall time, convergence, counters) at
   end-of-run.
