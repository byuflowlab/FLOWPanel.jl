# Handoff — BRAINSTORM 021: near-field caching in Phase 2, and the tuner's objective

Written 2026-08-24 at a context reset, superseding
`phase_14_cost_tuner_relaunch_prompt.md` (still accurate on environment, traps
and walltimes — read it for those; its *diagnosis* is revised below).

**Nothing is queued. Nothing is running.** Cluster disk 183 G / 2.0 T.

---

## 1. Ryan's ruling — do this first

> Make sure caching is performed in **both tuning and benchmarks for Phase 2**.
> **Ignore the caching cost while tuning.** **Measure it for benchmarks, but keep
> it separate from the per-step solve cost.**

Concretely that means:

- Phase 2's tuning must run with the near-field cache active, and the cache
  build time must NOT enter the tuning objective. This is legitimate because
  the panels-on-self operator has frozen geometry — the cache is built once a
  priori and reused across every solve iteration and every timestep.
- Phase 2's benchmark rows must run cached too, and must report the cache build
  cost as its **own column**, never folded into `t_solve` / per-step cost.
  A reader must be able to see both "what one step costs" and "what the one-off
  cache build cost" without them being mixed.
- Ryan also confirmed (and the code agrees): **only panels-on-self is worth
  caching.** Wake→body geometry moves every step, so those influence matrices
  change every step and are never reused inside a solve. The campaign already
  freezes `b` (`bcache_R*.bin`), which is exactly this fact.

### What exists today (verified 2026-08-24, do not re-derive)

- `benchmark/rotor_hover_solver_phase2.jl` — the `*_nfcache` config variants
  (`krylov_gmres_nfcache`, `krylov_ilu_nfcache`, `fgmres_fgs_nfcache`) already
  run cached: `nfcache_setup()` at line 183 builds an `FmmPlan` and calls
  `FastMultipole.build_nearfield_cache!` (line 193), then passes
  `cache_nearfield=true` (lines 218, 301). The non-`nfcache` configs do not.
- `benchmark/rotor_hover_solver_phase2b_nearfield_cache.jl:195` — the cached
  tune (`TUNE_CACHED=1`) calling stock `FastMultipole.tune_fmm` with
  `tune_nearfield_cache=true`, writing `tune_cached.csv`, which phase2 reads at
  line 164. Its header comment already says "build cost excluded", so the
  ignore-cost-while-tuning half of the ruling is ALREADY the policy there.
- **The gap:** that cached tune uses **stock `tune_fmm`**, i.e. the
  accuracy-only cost-model objective — precisely the defect fixed for the
  uncached path on 2026-08-24 (see §2). And `tune_fmm_perturb` has **no cache
  support at all**, so the cost-descent tuner cannot currently be pointed at
  cached economics.
- **Phase 1 caches nothing**, by design. Grep of `phase1_case.jl`,
  `phase1_knobs.jl`, `rotor_hover_solver_phase1.jl` for
  `cache_nearfield`/`nfcache` returns nothing.

### Open question Ryan has NOT ruled on

Whether the **Phase 1** ladder should also be tuned/run cached. His ruling names
Phase 2 only. Phase 1 currently spans all seven rungs uncached; Phase 2's cached
variants are a subset. **Ask before changing Phase 1.**

### Feasibility you should measure before committing

Cache bytes scale roughly as `panels x leaf`. The one recorded datum is
**R4 (58,192 panels) ~ 4.5 GiB**; naive scaling puts R7 (419,276 panels) near
**~32 GiB** at the same leaf, and higher if the tuned leaf grows — which it
will, since caching makes larger leaves cheaper. That is fine on a 500 GB node
but at/above the campaign's current 32 GiB cap, and **the cache build is
serial**, so `max_build_time` may bind before memory does.
`FastMultipole.estimate_nearfield_cache(plan.target_tree, plan.source_tree,
plan.direct_list, plan.derivatives_switches, (rotor,))` returns
`(bytes, est_build_time)` without building anything — use it to find where the
cache caps out across R1–R7 before deciding the cached arm's rung coverage.
`tune_fmm` returns a third value `tune_info.cache_capped` flagging when leaf
growth was clamped by the caps rather than by the cost optimum; propagate it.

---

## 2. State of the cost-tuner work (landed, synced, committed)

### The diagnosis, revised

`phase_14`'s claim was that the descent stalled on **noise** and `reps=2 -> 7`
was the fix. That was half right. The dominant defect was **bias**:
`tune_fmm_perturb` timed each candidate through FastMultipole's `Cache` path,
which rebuilds both octrees AND the interaction lists inside every timed call.
A steady solve does not pay that per apply — `src/FLOWPanel_fmm.jl` builds one
`FmmPlan` and reuses its trees/lists across every iteration. Build cost RISES as
leaf shrinks, so the objective carried a leaf-dependent penalty pushing the
answer toward large leaves.

### What shipped

`FastMultipole/src/autotune.jl`, `tune_fmm_perturb` (commits `fd52a2ad`,
`150b3fc3`, `40b52902` on branch `flowpanel-20260817`):

- **`tree_amortization::Real=1`** — applies sharing one tree/list build;
  candidate cost is `t_build/tree_amortization + t_apply`.
  `1` = legacy `Cache` path verbatim, correct when geometry moves every step
  (**this is the value for tuning a PARTICLE field**). `Inf` = charge nothing,
  correct for the frozen panels-on-self operator — **Ryan's ruling, and the
  driver default**. Finite `n` = a build reused over `n` applies. For any value
  `!= 1` the plan is built once per candidate, timed separately, released.
- **`abandon_factor=1.3`** — stops a trial once its running min passes 1.3x the
  fastest complete, error-satisfying candidate so far. Safe because acceptance
  needs `t < t0*(1-improve_tol)` and `t0` is never below that best; the
  threshold only tightens. Flagged `abandoned=true` in `history`; a candidate
  abandoned at the build stage has NOT had its error tolerance verified.
- **`max_seconds=Inf`** — absolute backstop, best-so-far + `@warn` + `timed_out`.
- Third return value `info = (timed_out, t_elapsed, n_candidates, n_abandoned,
  t_best)`. Two-value destructuring still works.
- `FMMPLAN_STRUCTURAL_KWARGS` in `src/fmm.jl` splits caller kwargs between
  `FmmPlan` (no catch-all `optargs`) and the per-apply `fmm!`.

`benchmark/rotor_hover_solver_phase1_tune_hpc.jl`:

- `TUNE_REPS` (5 at R1–R4, 2 at R5–R7, **asserted <= 5**, Ryan's cap),
  `TUNE_TREE_AMORTIZATION` (**default `Inf`**), `TUNE_ABANDON_FACTOR` (1.3),
  `TUNE_MAX_SECONDS` (900/1800/3600/7200/14400/21600/43200 s for R1–R7).
- `tune.csv` gained `tune_reps,tree_amortization,tune_timed_out` plus a header
  guard refusing to append to a pre-fix schema.
- `benchmark/phase1_case.jl`: hardcoded R1–R3 `TUNED` fallbacks now labelled as
  pre-fix accuracy-only picks, not valid campaign knobs.
- `benchmark/p023_fmm_tune.jl`: left at `tree_amortization=1` deliberately (023
  prices a mature-wake unsteady step, geometry moves every step) with a comment.
- `benchmark/slurm/fm_leaf_ab.sh`: NEW wrapper, reuses `p1_tune.sh`'s pinned
  hardware header verbatim, 2 h walltime. Untracked like the rest of
  `benchmark/`.

---

## 3. THE BLOCKER — the tuner's objective still disagrees with a real solve

**This is the open technical problem. Do not relaunch the ladder until it is
settled.**

At R1, p=17 / MAC=0.5:

| leaf | 9 | 21 | 40 | 60 | 71 |
| --- | --- | --- | --- | --- | --- |
| real solve, `t_solve_min` [s] | **15.51** | 20.65 | 26.10 | 32.64 | 34.95 |
| per iteration (niter ~57) | 0.272 | 0.362 | 0.458 | 0.573 | 0.613 |

Monotone down to leaf 9. But the tuner (job 13447582, `tree_amortization=Inf`,
`reps=5`, same p/MAC) measured **leaf 43 = 0.434 s and leaf 29 = 0.435 s** —
flat — and stopped after one move. Full trace in
`~/flowpanel-021/FLOWPanel.jl/slurm-p1-tune-R1-13447582.out`.

The tuner minimizes **one `fmm!` apply**; Phase 1 publishes **solve wall clock**.
That proxy is only worth anything if it is monotone in the same direction, and
below leaf ~40 it is not.

### Ruled out — do NOT re-investigate

- **`tune=true` overhead.** The trials pass it, real solves do not, and it adds
  a counting pass over `direct_list` whose length grows as leaf shrinks.
  Measured: **<=3%, and only at leaf 9.** Not the cause.
- **Per-candidate plan allocation / NUMA first-touch.** Hypothesis was that a
  fresh `FmmPlan` per candidate gives each one a persistent placement penalty
  that `min`-of-reps cannot average out. Tested: 4 independent plans at
  identical knobs vs 1 plan re-measured 4 times — **both scatter the same
  amount**, so it is not plan identity.

### Still live

Comparing the tuner against the solve-implied curve at matched leaf gives
+3% / -9% / +9% at leaf 64 / 43 / 29 — roughly +-9% scatter, exactly enough to
erase the true ~15% gap between leaf 43 and 29 and stall a descent with
`improve_tol=0.02`. Whether that is resolution or a real divergence between
apply cost and solve cost is **unknown**.

### CRITICAL instrument warning

**Ryan's laptop cannot measure this.** Measured 2026-08-24: min-of-5 apply
timings on the local Mac scatter **22–39%** at fixed knobs, against a 15% effect.
Any local timing "verification" is meaningless. **All timing evidence must come
from an exclusive cluster node.** (An earlier local `tree_amortization`
1/50/Inf comparison was reported in this session as verification; it was inside
that noise band and should be disregarded — the `Inf` code path is verified to
RUN correctly, its effect on the answer is not verified.)

### Two ways forward — Ryan has NOT chosen

He was asked and the question was interrupted, so **ask again before building**:

**(a) Time a real solve per candidate.** Replace the apply proxy with the actual
objective. `tune_fmm_perturb` would take a caller-supplied cost function; the
driver passes a closure building a `KrylovSolver` at the candidate knobs and
timing `_solve!`, exactly as `benchmark/fm_leaf_ab.jl:28-35` already does. Keeps
the descent, abandonment and accuracy-rejection machinery. Removes the proxy
question outright and makes the objective identical to what the campaign
publishes — and makes cached-vs-uncached fall out for free, which matters for §1.
Cost with `reps=2` and ~20 candidates, roughly halved by early abandonment:
R1 ~10–25 min, R2 ~25–55 min, R3 1–2 h, R4 2–4 h, R5 4–8 h, R6 8–16 h — all
inside the budgeted tune walltimes. **R7 is the only rung this does not fit.**

**(b) Debug the proxy.** Add per-rep printing (not just the min) and a fixed
diagnostic leaf ladder at the winning p/MAC to the **3.5-minute** R1 tune job.
If leaf 29's five reps are tight around 0.435 the measurement is real and apply
cost genuinely is not what the solve pays; if they scatter it is resolution.
Cheaper, and a working proxy is the only thing that scales to R7.

A hybrid is attractive: (a) for R1–R4 where it is affordable, which also
*validates* whether the proxy agrees — and if it does, use the proxy for R5–R7
with evidence behind it.

---

## 4. Result worth keeping: the PUBLISHED ladder is mis-tuned

Job 13443627, `fm_leaf_ab.jl` at R2 (15,760 panels), p=17 / MAC=0.5:

| leaf | 9 | 21 | 40 | 60 | 71 |
| --- | --- | --- | --- | --- | --- |
| `t_solve_min` [s] | **36.71** | 46.72 | 54.58 | 73.17 | 80.20 |
| niter | 62 | 63 | 62 | 62 | 60 |

Monotone, **2.18x at constant iteration count**. R2's old campaign pick was leaf
**59** — ~73 s where ~37 s was available. This was previously only inferred from
R1; it is now measured. R3 (old pick leaf 58) is affected by the same argument.

**Consequence:** apply knobs are shared across configs within a rung, so the
inflation hits every FMM-based config alike, but `backslash_ldiv` runs no FMM in
its solve and is untouched. The published rungs therefore biased the comparison
**against** the iterative FMM solvers by up to ~2x. Any backslash-vs-Krylov
crossover conclusion from those rows must be re-read against the retuned ladder.

(`fm_leaf_ab.jl`'s `bc_rel_l2` column is INVALID — wrong args to `bc_error!` —
and read ~1.02–1.09 here. Only its timings are sound.)

---

## 5. Cluster state

- `~/flowpanel-021/` — FLOWPanel `7af4271` (branch `fastmultipole`),
  FastMultipole `40b5290` (branch `flowpanel-20260817`), both clean, both
  md5/hash-verified against local. `fm_commit` resolves `v2.3.0-222-g40b5290`,
  **no `-dirty`**.
  - The FLOWPanel banner's `commit` field always shows `-dirty` on this clone.
    That is structural: `_git_describe` (`benchmark/common.jl:95-104`) flags on
    `git status --porcelain`, which counts untracked files, and `benchmark/` is
    untracked there by design. Not a problem; do not chase it.
  - Local FastMultipole has one newer commit, `b332bb05` (051 docs), added by a
    concurrent session. Docs only, not synced, harmless.
- **Syncing tracked code**: use `git bundle`, not a GitHub push and not an
  rsync of tracked files. `git bundle create X.bundle <base>..<branch>`, `scp`
  it, then `git fetch /tmp/X.bundle <branch> && git merge --ff-only FETCH_HEAD`.
  Identical hashes on both sides, clone stays clean, nothing goes public.
  Untracked `benchmark/` files go by `rsync` + `md5sum` check as before.
- Results tree: only `benchmark/results/phase1/multi/R1/tune.csv` exists (from
  job 13447582, `tree_amortization=Inf` -> p=17/MAC=0.45... see the row; it is
  the SUPERSEDED-pending-§3 result). All seven `bcache_R*.bin` intact. Purge a
  rung's CSVs before re-running it — the drivers skip rows already present.
- **Use `--partition=physics2 --qos=standby`.** The default `m12`/normal left
  both jobs PENDING on Priority; moving them to physics2/standby started them
  within seconds (8 idle nodes). `scontrol update jobid=N partition=physics2
  qos=standby` works on a pending job and preserves the job ID — no resubmit.
  Everything else in the SBATCH header is pinned hardware: do not override
  `--constraint=zen3 --exclusive --mem=500G`.
- Job submission is inline work and **Ryan's call — get his go before `sbatch`,
  never route it to a subagent.**

## 6. Reading order for a fresh agent

1. This file.
2. `log.md` (**newest first**) — entries for 2026-08-24 carry the bias finding,
   the abandonment design, the `Inf` ruling, and the R2 measurement.
3. `phase_14_cost_tuner_relaunch_prompt.md` for environment, the three traps
   (`afterany` vs `afterok`, sticky resume rows, `nproc` lying) and the
   right-sized walltime table. Its diagnosis section is superseded by §2 here.
4. Do NOT read the item file's `## Current status` (stale), and do not read
   `phase_01_consistency.md` (73 KB) or `phase_02_single_step_benchmarks.md`
   (41 KB) end-to-end — grep them.
