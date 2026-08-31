# 018 GPU — NT144 performance cliff: findings + handoff (2026-08-27)

> **RESOLVED (same day, later session):** root cause was NOT the counting
> sort — it was rVPM sigma growth of one runaway particle crossing the
> radix-FMM ell=3 adequacy limit, forcing an ell=2 rebuild. Confirmed
> quantitatively and fixed (σ-cap band-aid + 3R/λ2.0 ladder relaunch). See
> `gpu_3r_ladder_status_20260827.md` and BRAINSTORM item 026. The document
> below is the historical investigation record.

Status: **root cause NOT confirmed.** A sharp, reproducible performance cliff
blocks all NT144 GPU arms. Mechanism narrowed to one leading hypothesis;
the test that would confirm it has not yet been run successfully.
Companion doc: `gpu_launcher_handoff_20260827.md` (launcher/silo/env state).

---

## 1. The cliff (verified)

The NT144 arm runs at 6–8 s/step, then **abruptly and permanently** drops to
~52–54 s/step. Measured on GH200 (job 13502380), from the wake-health monitor
of `p018_csarc_n5_nt144_l2p4_s2gpu`:

| step | n_particles | wall (s) |
|---|---|---|
| 2251 | 266,240 | 8.2 |
| 2252 | 266,250 | 6.0 |
| 2253 | 266,288 | 24.3 |
| 2254 | 266,347 | 54.4 |
| 2255 | 266,365 | 54.0 |
| 2260 | 266,482 | 51.9 |

**Threshold ≈ 266,300 particles**, resolved to ~100 particles. Not gradual
degradation — a step change, then flat.

### Verified properties

- **Persistent, not process state.** Killed and restarted from a post-cliff
  snapshot twice (13502552 from step 2278, 13502600 from 2288). Both ran at
  52.2 s/step **from their first step** — no fast phase. A fresh process does
  not clear it.
- **Not hardware / not GH200-specific-so-far.** Same-case A/B (job 13502563,
  `p018_csarc_n2_nt72_l3p0` on GH200) vs the live H200 job 13502383, compared
  at **matched step numbers**: GH200 2.49–2.61 s/step, H200 2.81–2.92 s/step.
  GH200 is marginally *faster*. The node had also completed a full 30-rev NT72
  arm at 6.2 s/step hours earlier.
- **Not memory.** Device: 95 GB free, pool reserved 78 GB / used 15.9 GB,
  steady across the transition. No OOM, no pool exhaustion, no ECC errors,
  no throttling (clocks at max 1980 MHz, 28 °C).
- **Not the physics.** Across the transition `n_particles` grows smoothly
  (~30–60/step), `max_u` steady ~12 m/s, `min_sigma` smooth. No blow-up.
- **Not GPU contention.** Sole compute process on the GPU.

### Profiling evidence (perf, read-only, no ptrace stop)

`perf record -F 199 -p <pid>` on the live slow process:

- ~72 % of samples in one tight loop in `libcuda.so` (addresses clustered
  within ~100 bytes) + `vdso` `clock_gettime` — a **driver busy-wait spin**.
- Call graph names the callers: **`cuStreamSynchronize`** and
  **`cuMemcpyHtoDAsync_v2`**. Host→device transfers plus synchronisation,
  *not* kernel compute.
- Exactly **one** Julia thread at 99.9 % CPU, other 182 asleep; node load 2.6
  of 72 cores. (Argues against GC: Julia 1.11 would light ~half the pool.)
- `nvidia-smi` reports "100 % utilization" but **power draw is only
  113–165 W of a 900 W limit** — the GPU is near-idle and being poked by
  small transfers. Do not trust the utilization number here.
- The GPU source-influence solve is still healthy (`source_influence_s_gpu_gemv`
  ~3 ms, firing once per step). The regression is elsewhere.

## 2. Why only NT144 (verified)

Particle counts at matched revolutions (N differs by design; the NT72 N=2 vs
N=3 arms agree within 0.6 %, so **NT drives this, not N**):

| rev | NT36 (N=1) | NT72 (N=3) | NT144 (N=5) |
|---|---|---|---|
| 2 | 31,558 | 34,134 | 40,302 |
| 6 | 112,680 | 120,218 | 138,219 |
| 10 | 184,385 | 196,222 | 222,550 |
| 12 | 197,478 | 213,996 | 246,462 |
| 14 | 199,340 | 219,514 | 261,491 |
| 20 | 188,028 | 216,966 | — |
| 28 | 183,865 | 226,495 | — |

- NT36 peaks ~199k (rev 14) then settles 180–190k. **Never reaches the cliff.**
- NT72 plateaus 217–227k (max 226,846 at rev 30). **Never reaches the cliff.**
- NT144 **never plateaus** — still climbing steeply; crosses ~266.3k at
  rev ~15.6 and only goes further above it.

So NT36/NT72 are structurally safe; NT144 operates above the threshold for the
back half of its run. `NT × PPS = 432` does **not** hold particle count
constant in practice (NT144 is ~21 % above NT36 at rev 10) — worth knowing for
the campaign independently of this bug.

Incidental cross-validation: at rev 10 the CPU NT144 run reports 222,550
particles vs the GPU run's 222,549 — 1 particle apart after 1440 steps.

## 3. Leading hypothesis (NOT confirmed)

`FastMultipole/src/translate_batched_cuda.jl:171–181`:

```julia
const RADIX_CUDA_COUNTING_SORT = Ref(true)
const RADIX_CUDA_COUNTING_SORT_MAX_ELL = Ref(6)
_cuda_counting_sort_enabled(ell) =
    radix_setting(:RADIX_CUDA_COUNTING_SORT) && ell <= radix_setting(:RADIX_CUDA_COUNTING_SORT_MAX_ELL)
_cuda_counting_sort_ready(ctx, ell) =
    _cuda_counting_sort_enabled(ell) && length(ctx.counting_histogram) == 1 << (3ell)
```

At `:6668` a failed `_cuda_counting_sort_ready` falls back to a generic
`sortperm!` (`:6302`). The comment states "Falling back is always safe" — i.e.
**silent**. Two ways to trip it: `ell > MAX_ELL` (=6), or the
construction-time histogram no longer matching the current depth.

Supporting: `8^6 = 262,144` sits right against the measured ~266,300 cliff;
sharpness and permanence both fit a discrete depth increment.

**Unresolved tension:** on restart the cache is rebuilt at the current particle
count, so a pure histogram-size mismatch should self-heal — yet restarts are
slow immediately. That points to `ell > MAX_ELL` rather than the size check,
but this is unverified.

**Critically: `ell` is never logged.** No `ell`/depth marker exists anywhere in
the run output (`depth=4` in the logs is `TRUNCATION_DEPTH_R`, unrelated).
This is the single biggest gap.

## 4. Failed bisection attempt — read before retrying

Job 13502600 tried `RADIX_CUDA_COUNTING_SORT_MAX_ELL=8` from the post-cliff
snapshot. **The result (52.5 s/step) is NOT a disproof — the harness did not
take effect.**

Mechanism built (deployed, uncommitted):
- `examples/run_dji9443_hover_ct_gpu.slurm.sh` — added
  `P018_JULIA_OVERRIDE` (inert when unset; both arch branches).
- `orc:~/p018_julia_radix.sh` — shim: `exec julia -L ~/p018_radix_preload.jl "$@"`.
- `orc:~/p018_radix_preload.jl` — maps `P018_RADIX_*` env → `set_radix_setting!`.

Verified the shim WAS used: the live process argv was
`julia -L /home/rander39/p018_radix_preload.jl --project=... -t 72 examples/rotor_hover_pressure_comparison.jl`.
Yet the preload emitted **zero output** — including its unconditional echo
loop. Job stderr is 32 lines, all "redefinition of constant". Unexplained.

**Related finding, independent of this campaign:** those warnings are the
FastMultipole CUDA settings `Ref`s being **re-created** because
`translate_batched_cuda.jl` is included more than once per process. Any value
set before that re-include is **reset to its default**. This makes the
documented `set_radix_setting!` surface unreliable from a preload. Pre-existing
(present in every prior GPU run incl. the validated smokes). Relevant to
tasks 047/052 settings hardening.

## 5. Recommended next step

Instrument rather than guess. In the silo FastMultipole
(`~/FastMultipole-018-gpu-gh200`), add a rate-limited `@info` near
`translate_batched_cuda.jl:6668` printing `ell`, `MAX_ELL`,
`length(ctx.counting_histogram)`, `1 << (3ell)`, and which branch was taken.
Run ~20 steps from the post-cliff snapshot (step 2288). That answers the `ell`
question directly and tells us whether the fallback is even involved.
**This is science-adjacent source; Ryan's approval needed** (why it stopped here).

If confirmed and `ell = 7`: raising `MAX_ELL` to 7 costs a `8^7 = 2,097,152`
entry histogram (~25 MB) — trivial on 98 GB. Would unblock all three NT144 arms.
If not the sort: next candidates are CUDA graph-capture invalidation (several
settings are "captured"; losing capture means per-step launches + more sync,
which fits the `cuStreamSynchronize` profile), `CUDA_NEARFIELD_SHAPE`, and
`CUDA_NEARFIELD_SUBSORT` (documented as running "during the (uncaptured)
**host refresh**" — would directly explain the HtoD traffic).

## 6. Job ledger (2026-08-27)

| job | arm | outcome |
|---|---|---|
| 13502379 | NT72 λ2.4 restart → `_s2gpu` @1162 | **DONE, rev 30**, 1:57:40, 6.2 s/step |
| 13502381 | `csarc_l3p0` NT36 N=1 λ3.0 | **DONE, 30 revs**, CT̄ rev30 = 0.070884 |
| 13502382 | `csarc_l4p0` NT36 N=1 λ4.0 | **DONE, 30 revs**, CT̄ rev30 = 0.071048 |
| 13502383 | `n2_nt72_l3p0` | RUNNING, healthy 5.6–6.9 s/step |
| 13502384 | `n2_nt72_l4p0` | PENDING (eng) |
| 13502385/6 | `n4_nt144_l{3p0,4p0}` | PENDING — **will hit the cliff ~rev 15.6** |
| 13502567 | NT144 restart `_s4h200` @2288, H200/eng | PENDING — **architecture test** |
| 13502380 | NT144 `_s2gpu` @1162 | cancelled at the cliff (snapshots kept) |
| 13502552 | NT144 `_s3gpu` @2278 | cancelled — slow from step 1 |
| 13502563 | GH200 x-check of `n2_nt72_l3p0` | cancelled — GH200 healthy |
| 13502600 | bisect `MAX_ELL=8` | cancelled — **harness inert, invalid** |

NT144 progress banked: rev ~15.9 (step 2288, `_s3gpu`, gh200 silo + staged in
h200 silo). Latest complete restartable step **2288**.

Slurm "FAILED" on the three completed arms is **cosmetic** — the gate tripped
on the driver's by-design `Bernoulli vs KJ: NaN` summary line (KJ disabled),
which is not `|`-separated and slipped the CT-table filter. Gate fixed and
mirrored (md5 `3e660960b6f8fe5c9e22396a73357a0a`); substantive checks all
passed (gpu_gemv 997–1080, cpu_gemv 0, rc 0).

## 7. Uncommitted local changes

- `examples/run_dji9443_hover_ct_hpc.slurm.sh` — 6 new Phase-17b case blocks
  (`csarc_{l3p0,l4p0,n2_nt72_l*,n4_nt144_l*}`, λ=3.0/4.0, N=1/2/4 by NT).
  Mirrored md5-identical to CPU tree + both silos (`b909b52ee53f53cbccaecd73db0fc1cf`).
- `examples/run_dji9443_hover_ct_gpu.slurm.sh` — gate NaN fix + inert
  `P018_JULIA_OVERRIDE`. Mirrored to both silos.
- Ledger/notebook entries NOT written (needs Ryan's approval).

## 8. Local artifacts

`~/p018_nt144_cliff_2250_2260/` (409 MB) — ParaView files for steps 2250–2260
of `_s2gpu`, bracketing the cliff (2 clean steps before, 6 after): 11 particle
`.vtp`, 22 wake `.vts`, 55 body/tw `.vtu`. Open the file series directly; the
included `.pvd` files index the full run and will complain about absent steps.
Ryan was inspecting these for anything visually obvious — result not yet known.
