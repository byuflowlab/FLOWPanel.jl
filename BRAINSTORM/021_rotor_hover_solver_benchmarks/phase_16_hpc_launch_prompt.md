# Handoff — BRAINSTORM 021: launch the Phase 1–2 HPC campaign

Written 2026-08-24 at a context reset. **Supersedes
`phase_15_caching_and_objective_prompt.md`** (read its §5 for cluster state and
`phase_14_cost_tuner_relaunch_prompt.md` for environment/traps/walltimes; their
*diagnoses* are settled below).

**Your job: sync the code to the cluster and launch the Phase 1 and Phase 2
jobs described in §4.** Everything is implemented and locally smoke-tested.
Nothing has ever run on the cluster. **Nothing is queued. Nothing is running.**

> **`sbatch` IS RYAN'S CALL.** Prepare the submissions, show him the exact
> command lines and the expected cost, and get his explicit go before
> submitting. Never route submission to a subagent.

---

## 1. What changed, in one page

Ryan's rulings, all from 2026-08-24:

1. **Caching happens in BOTH Phase 2 tuning and Phase 2 benchmarks.** The cache
   build cost is **excluded from the tuning objective** and **measured
   separately in the benchmarks**, never folded into per-step cost.
2. **The tuning objective is a REAL SOLVE**, not one `fmm!` apply, for R1–R5.
   (R6–R7 deferred — see §6.)
3. **Memory is a swept AXIS, expressed as machines**: 0 (cache forbidden) /
   16 (laptop) / 128 (workstation) / 500 (node) GiB of **TOTAL** memory.
4. **Error is measured against the BOUNDARY CONDITION, not a reference**, and
   evaluated through the FMM (not `DirectBackend`, which is O(N²)).
5. **Phase 1 is descoped to verification**: solve the ladder up to R7 and
   confirm the solvers agree on the SOLUTION. Phase 2 owns tuning.

### Why the memory axis exists

An uncapped cached tune **provably walks to a dense operator** — measured at R1,
the cached optimum sat at 1.3–1.7× dense with `cache_capped=false`, so the caps
never bound; the cost optimum genuinely wanted it. That is structural: if
assembly is free (and ruling 1 makes it free for tuning), matrix-ful always
beats matrix-free wherever memory allows. So a memory cap does not *guard* the
answer, it *determines* it. Hence: sweep it, and sweep TOTAL memory (body +
solver state + FMM plan + cache), which puts `backslash_ldiv`'s dense 8N² on the
same axis. Whether a candidate is cached is an **outcome** of what fits.

---

## 2. Files (all implemented, all locally smoke-tested)

| File | Role |
| --- | --- |
| `benchmark/nfcache_feasibility.jl` | **W0 map.** `estimate_nearfield_cache` over a leaf grid; nothing built, nothing solved. Deterministic → valid locally. Already run for R1–R7; results in `ledger.md`. You should not need to re-run it. |
| `benchmark/rotor_hover_solver_phase2_tune.jl` | **Phase 2 tuner.** Real-solve objective over the machine-class memory ladder. Writes `tune_phase2.csv`. |
| `benchmark/rotor_hover_solver_phase2.jl` | **Phase 2 benchmarks.** Reads `tune_phase2.csv`; `*_nfcache` rows are WARM with the build in its own column. |
| `benchmark/rotor_hover_solver_phase1_agreement.jl` | **Phase 1 verification.** BC metric + reference-free ensemble spread. |
| `benchmark/slurm/p2_tune.sh` | NEW wrapper for the tuner. One job per (rung, budget). |
| `benchmark/slurm/p1_agreement.sh` | NEW wrapper for Phase 1 verification. |
| `FastMultipole/src/autotune.jl` | `tune_fmm_perturb(...; cost=)` — the caller-supplied objective. |
| `FastMultipole/test/autotune_cost_test.jl` | 13 tests for it, wired into `runtests.jl`. |

Both new wrappers reuse `p1_tune.sh`'s pinned-hardware preamble **verbatim**
(zen3 / `--exclusive` / `--mem=500G` / 128-CPU asserts / 64 threads / gaussian
filament reg). `benchmark/` is untracked on the cluster by design.

---

## 3. The settings, exactly

**Mesh ladder** (frozen, unchanged): R1–R7, 8,016 → 419,276 panels.

**Memory ladder** — `MEM_BUDGETS = 0:16:128:500`, GiB of TOTAL memory:

| budget | meaning |
| --- | --- |
| **0** | cache FORBIDDEN — the matrix-free endpoint, not a machine. Memory unconstrained. **REQUIRED at every rung and must run FIRST**: it is what tunes the uncached configs, and at R1–R5 every real machine can afford a cache so no positive budget ever returns an uncached winner. `phase2.jl` errors without it. |
| **16** | laptop |
| **128** | workstation |
| **500** | supercomputer node (the pinned `--mem=500G` on 512 GB zen3) |

What a class MEANS varies by rung — that is the finding, not a defect:

| rung | laptop 16G | workstation 128G | node 500G |
| --- | --- | --- | --- |
| R1–R3 | saturated | saturated | saturated |
| R4–R5 | **constrained** | saturated | saturated |
| R6–R7 | **cannot cache at all** | constrained | saturated |

R6/R7 saturate at ~158 and ~428 GiB, both under 500, so the whole axis is
reachable on the pinned hardware. Do **not** chase the 1 TB (`cs`) or 2 TB
(`m11-2`) partitions: that changes hardware and destroys timing comparability
(`m11-2` is zen2 outright). Saturated/duplicate rows are flagged by the driver.

**Other knobs**: `TUNE_REPS` 5 at R1–R4 / 2 above (asserted ≤5, Ryan's cap);
`TUNE_ABANDON_FACTOR=1.3`; `FLOOR_LEAVES=2:4:9:21`;
`NFCACHE_MAX_BUILD_TIME=Inf` (the build is SERIAL — consider capping it at
R6/R7); `target_rel=1e-6`, BC certified at `0.1×` that.

---

## 4. What to launch

Order matters only in that **budget 0 must land before `phase2.jl` runs**, and
Phase 2 tuning must land before Phase 2 tables. Budgets are independent → one
job per (rung, budget), which also suits preemptible standby.

**Step A — Phase 2 tuning, R1–R5 × 4 budgets = 20 jobs.**

```
sbatch --job-name=p2-tune-R1-b0 --time=4:00:00 \
  --export=ALL,RUNG=R1,MEM_BUDGETS=0 benchmark/slurm/p2_tune.sh
```
…and so on for `MEM_BUDGETS` ∈ {0,16,128,500} × `RUNG` ∈ {R1..R5}.

Per-budget walltime guide (from `phase_15` §3, one real-solve descent):
R1 10–25 min · R2 25–55 min · R3 1–2 h · R4 2–4 h · R5 4–8 h. Request ~3× the
guide. Whole step ≈ 38–77 h of node time, fully parallel across jobs.

**Step B — Phase 1 verification, R1–R7 = 7 jobs** (independent of A; can run
concurrently).

```
sbatch --job-name=p1-agree-R7 --time=2-00:00:00 \
  --export=ALL,RUNG=R7 benchmark/slurm/p1_agreement.sh
```

- **Drop `backslash_ref` and `backslash_coupled` above R4** via `CONFIGS` —
  they assemble a dense 8N² G (25 GiB at R4, 87 at R5, 335 at R6, **1.31 TiB**
  at R7). Nothing requires them any more; run them where they fit as an
  independent non-iterative cross-check.
- `krylov_ilu`'s pattern is ~23 GiB of entries at R7 — it may not fit some
  machine classes. That is a result, not a failure.
- Run `RESIDUAL_BACKEND=direct` **once at R1 or R2** as a cross-check that the
  FMM BC metric matches the exact O(N²) one. Do not attempt it above R3.

**Step C — Phase 2 tables**, after A lands for a rung. `p2_table.sh` /
`rotor_hover_solver_phase2.jl`, `MEM_BUDGETS=16:128:500`. **The pre-W3
`phase2.csv` must be moved aside first** — the new schema is header-guarded and
the driver will refuse to append (see §5).

---

## 5. Traps — read before submitting

1. **Quarantined local artifacts.** Everything below is LOCAL smoke output and
   must NOT reach the cluster results tree or be read as campaign data:
   `tune_phase2_localsmoke{,2,3,4}.csv`, `phase2_w3_localsmoke.csv`,
   `agreement_localsmoke.csv`, `agreement_spread_localsmoke.csv`,
   `phase2_prefix_w3.csv.bak`, `agreement_pre20260824.csv.bak`,
   `tune_cached.csv` (accuracy-only objective, superseded).
   **Local timings scatter 22–39% at fixed knobs on Ryan's Mac** against a ~15%
   effect. No local timing is evidence. Ever.
2. **Schema guards will refuse to append** to pre-2026-08-24 `phase2.csv` and
   `agreement.csv`. Move them aside; do not delete.
3. **Preemptible by default — already in the headers.** `p2_tune.sh`,
   `p1_agreement.sh` and (patched 2026-08-24) `p2_table.sh` all carry
   `--partition=physics2 --qos=standby --requeue` in the SBATCH header, so you
   do not have to type the override. `m12`/normal left the 08-24 jobs PENDING
   on Priority indefinitely; physics2/standby started them within seconds. The
   three partitions are spec-identical (2x64 zen3, 512 GB), so this does not
   affect comparability. If you submit any OTHER wrapper (`p1_tune.sh`,
   `p1_table.sh`, ...), it still defaults to m12/normal — pass the override, or
   `scontrol update jobid=N partition=physics2 qos=standby` on a pending job
   (preserves the job ID).

   **But know what standby costs here.** The tuner writes a budget's CSV row
   only AFTER that budget's descent finishes; nothing is written mid-descent
   and `tune_fmm_perturb` does not memoize across runs. So **resume
   granularity is one whole budget** — a preemption at R5 can throw away 4-8 h
   of descent, repeatedly. Mitigations, in order of preference:
   - R1-R3 (<=2 h per budget): standby is straightforwardly right.
   - R4-R5: set **`TUNE_MAX_SECONDS` below the typical preemption window**.
     It is a wall-clock backstop that returns best-so-far, WRITES the row, and
     flags `tune_timed_out=true` — converting "lose everything" into "keep a
     partial answer that is honestly labelled". A timed-out descent must never
     be read as converged (the driver warns).
   - If preemption still bites: `m12`/normal and accept the queue, or `m12pws`
     (28 d MaxTime, 1 node, standby).
4. **Never override** `--constraint=zen3 --exclusive --mem=500G`. Override
   `--time` and `--job-name` freely.
5. **Chain with `afterany`, never `afterok`** (Ryan's ruling): an `afterok`
   cascade turns one failure into N silent CANCELLED jobs, misdiagnosed twice
   on this campaign already.
6. **`nproc` lies under Slurm** — use `SLURM_CPUS_ON_NODE`. Already handled in
   the wrappers.
7. **Drivers APPEND and resume per row.** Purge a rung's CSVs before re-running
   it; keep `bcache_R*.bin`.
8. **Sync tracked code by `git bundle`**, not a push and not an rsync of
   tracked files: `git bundle create X.bundle <base>..<branch>`, `scp`, then
   `git fetch /tmp/X.bundle <branch> && git merge --ff-only FETCH_HEAD`.
   Untracked `benchmark/` goes by rsync + `md5sum`. **Two repos must sync:
   FLOWPanel.jl AND FastMultipole** (`autotune.jl` + the new test file). After
   a partial sync, verify the WHOLE dep stack — import-clean ≠ runtime-clean.
9. `ssh orc` needs a live ControlMaster socket (2FA otherwise).
10. The FLOWPanel banner always shows `-dirty` on the cluster clone because
    `benchmark/` is untracked there. Structural; do not chase it.

---

## 6. Still needs Ryan — do not decide these yourself

- **R6–R7 tuning objective.** Real-solve tuning is affordable to R5. R6 is
  8–16 h *per budget*; R7 does not fit. Ryan deferred this until R1–R5 land.
  If the real-solve winner agrees with the apply proxy at R1–R5, the proxy can
  be used above with evidence behind it — check that first.
- **Anything past 500 GiB** (would require different hardware).
- **Every `sbatch`.**

---

## 7. Verification already done (do not redo)

- FastMultipole: 13/13 new `autotune_cost_test.jl`, 54/54 near-field cache
  tests, at `-t 4`.
- **D4 additivity, R1**: `t_cold` 16.059 s vs `t_build + t_warm` 16.024 s —
  **0.22%**. A cold solve really is just the build plus a warm solve, so the
  cached rows need only ONE timed pass. Ryan asked for this explicitly.
- **Build genuinely excluded from the objective**, structurally (the warm-up
  solve that builds the cache is not wrapped in `@elapsed`; only the rep loop
  is) and empirically (R1: build 11.7–14.2 s, `t_warm` 4.8 s — if the build
  were included `t_warm` could not be below 11.7 s).
- BC metric cross-validates: `krylov_gmres` at R1 reads 9.673e-7 from both the
  Phase 1 agreement driver and the independent Phase 2 tuner.
- Reference-free ensemble spread works across processes: R1 `fgs` vs
  `krylov_gmres` max pairwise 1.24e-5, both BC-certified below 1e-6.

## 8. Reading order

1. This file.
2. `ledger.md`, last ~5 sections (W0 map, machine ladder, additivity, the
   floor correction, quarantined artifacts).
3. `log.md`, newest first — the two 2026-08-24 entries.
4. `phase_14_cost_tuner_relaunch_prompt.md` for environment/traps/walltimes
   ONLY; its diagnosis is superseded.
5. Do NOT read the item file's `## Current status` (stale), and grep rather
   than read `phase_01_consistency.md` (73 KB) or
   `phase_02_single_step_benchmarks.md` (41 KB).
