# BRAINSTORM 021 — Phase 21: context-reset handoff

**Written 2026-08-27 14:10 MDT. Read this file first. `phase_20_context_reset_prompt.md`
§2–§5 are still current and still authoritative for the Phase 3 BC guard —
read them before touching anything in Phase 3. `phase_18_context_reset_prompt.md`
§2/§3 remain the reference for cluster env and the precompilation root cause.
Campaign rules: `decision_rules.md`. Harvest traps: `phase_17_harvest_prompt.md`
§4 and §7. You should NOT need `ledger.md`, `log.md`, or the item file to start.
Do not read `phase_01_consistency.md` (73 KB) or `phase_02_single_step_benchmarks.md`
(41 KB) end-to-end — grep them.**

---

## 0. START HERE — two things need Ryan's decision

1. **Disk is at 360 G of the 400 G cap** (`df -h /home/rander39`, 08-27 14:05).
   That is past the 300 G escalation line, and the `fp-018gpu-*` chain is
   actively writing. Policy: escalate retention to 72 timesteps, sweep with the
   `hpc-storage` subagent (restart-set-aware, deletes ParaView output only).
   **Offered to Ryan 08-27, not yet authorized. Ask before sweeping.**
2. **A notebook entry on the Phase 3 BC-guard rebuild has been offered twice and
   never written** (carried from `phase_20` §6.2). Nothing has gone into
   `~/Dropbox/research/notebooks/journals/` for this work. The finding is not
   "a threshold moved" — it is that a validity guard was asking a question the
   solver never answers. Offer, ask how much to write, never write unprompted.

---

## 1. Cluster access — read this before delegating any status check

**`squeue`/`sacct` are NOT on the non-interactive `PATH`.** They are shell
functions defined in the interactive profile, wrapping `/apps/slurm/latest/bin/`.
A plain `ssh orc 'squeue'` returns "command not found", which looks exactly like
a broken cluster. Use either form:

```bash
ssh orc '/apps/slurm/latest/bin/squeue -u rander39 -o "%.12i %.16j %.2t %.11M %.11l %R"'
ssh orc 'bash -lc "squeue -u rander39 ..."'
```

**The `hpc-monitor` subagent got this wrong twice on 08-25**, then — once given
the path — reported UTC/MDT-mixed timestamps and declared a stalled R7 "complete".
Its raw numbers were right; its arithmetic and verdicts were not. If you delegate,
**hand it the binary path up front and check its table against itself** (start
time vs mtime must be consistent). Short targeted inline `ssh` commands proved
cheaper and more reliable than the subagent for this particular job.

**The cluster runs MDT, and `squeue` `StartTime` is MDT.** File mtimes are MDT
too. `date -u` differs by +6 h. The `phase_20` job table's start times were
simply wrong (it said 22:32; actual 16:33).

**Connection gotcha (unchanged):** `ssh orc` needs a live ControlMaster socket or
it hits 2FA and hangs. If it stalls, ask Ryan to open a session.

---

## 2. Job state as of 2026-08-27 14:05 MDT — verify, do not trust

### Phase 1 agreement: COMPLETE for every rung

R1–R6 all completed 08-25. **R7 completed 08-26 08:54** (`13482349`, 10 h 18 m).

R7 landed as a **4-config ensemble**:

| config | t_steady (s) | niter | bc_rel_l2 | certified | dev_from_mean relL2 |
| --- | --- | --- | --- | --- | --- |
| krylov_gmres | 7315.3 | 86 | 9.05e-7 | true | 3.20e-5 |
| krylov_ilu | 3788.0 | 24 | 9.27e-7 | true | 1.06e-5 |
| fgs | 3166.4 | (n/a) | 2.19e-7 | true | 1.06e-5 |
| fgmres_fgs | 3123.6 | 13 | 6.40e-7 | true | 1.11e-5 |

**Ensemble max pairwise 4.31e-5 relL2, `n_configs=4`, all certified.** gmres is
the outlier; the other three agree to ~1.1e-5.

Those wall times are **not benchmarks** — `p1_agreement.sh` says so in its header
(Phase 1 was descoped 08-24; timing lives in Phase 2).

### Phase 2 tuning

| job | rung | state | budgets in `tune_phase2.csv` |
| --- | --- | --- | --- |
| 13477933 | R1 | COMPLETED 08-26 03:41 (11 h 07 m) | 0 / 16 / 128 / 500 |
| 13477934 | R2 | COMPLETED 08-26 06:27 (13 h 53 m) | 0 / 16 / 128 / 500 |
| 13477935 | R3 | COMPLETED 08-27 07:26 (1 d 14 h 52 m) | 0 / 16 / 128 / 500 |
| 13477936 | R4 | **RUNNING — REQUEUED**, restarted 08-27 06:32 on phys-1-2 | 0 only |
| 13477937 | R5 | RUNNING, 1 d 21 h 31 m, phys-1-7 | 0 / 16 / 128 |

**R4 was preempted and requeued** (standby QOS + `--requeue`), losing ~1.5 days
of wall clock; a campus power outage hit CTB/CB nodes around then and may be the
cause. Its 6-day timer restarted with it, so it is not at risk of timing out —
it is just far behind, and the tuner resumes candidates from its trace files
rather than starting cold.

Sample winners (`mem_budget_gib`, `p`, `mac`, `leaf`, `t_solve_warm`):

- R1 @500: `p=10 mac=0.5 leaf=72`, 1.25 s, niter 59
- R3 @128: `p=13 mac=0.5 leaf=27`, 6.34 s, niter 69
- R5 @128: `p=17 mac=0.6 leaf=32`, 34.8 s, niter 76

### Not 021, same account

`fp-022m-4r-oge` (1 d 09 h of 3 d) and ~8 `fp-018gpu-*` jobs (2 running, rest
pending on Priority/Resources).

---

## 3. What this session changed

**Ryan's ruling, 2026-08-25:** at R7, drop `krylov_jacobi`; keep vanilla
`krylov_gmres`, `krylov_ilu`, `fgs`, `fgmres_fgs`. **Scope is R7 only — jacobi
stays selectable at R1–R6**, where it converges in reasonable time.

Why: job `13469158` spent **9 h 20 m on `krylov_jacobi` with no result** after
`krylov_gmres` had solved the same rung in 7315 s / 86 iters. Ryan cancelled it.

Sequence executed (all authorized in the moment):

1. `scancel 13469158`.
2. Archived the existing gmres artifacts so they could not silently rejoin the
   ensemble, then **restored them** when Ryan ruled to keep vanilla gmres.
3. Resubmitted as `13482349` with
   `--export=ALL,RUNG=R7,CONFIGS=krylov_ilu:fgs:fgmres_fgs`.
4. On completion the job appended its 3 rows to the restored `agreement.csv` and
   the spread step picked up all 4 `.bin` files — the 4-config table in §2.

**Local repo edits (NOT on the cluster, `benchmark/` is untracked so `git diff`
shows nothing):**

- `benchmark/slurm/p1_agreement.sh:19` — R7 submit example now carries the
  explicit roster plus a "krylov_jacobi IS NOT SUPPORTED AT R7" note with the
  evidence and the colon-separator warning. Passes `bash -n`.
- `benchmark/slurm/SUBMIT_021_phase12.txt:60` — R7 split out of the shared
  R5–R7 `CONFIGS` line, same rationale. R5/R6 keep jacobi.

---

## 4. Traps

- **The spread step ignores `CONFIGS`.** `rotor_hover_solver_phase1_agreement.jl:492-497`
  globs every `agreement_x_<rung>_*.bin` in the results dir. A leftover solution
  from a cancelled run rejoins the ensemble silently. `agreement.csv` **appends**
  (`:61`), so stale rows persist too. Check the dir before relaunching any rung.
- **`CONFIGS` needs COLONS, not commas** (`:310-314`). A comma inside
  `--export=ALL,...` breaks sbatch's own parsing and the driver's `[:,]` split
  then selects **nothing** — silently.
- **R7's provenance is split across two job IDs.** `krylov_gmres` came from
  `13469158`, the other three from `13482349`. Legitimate — the design has no
  privileged reference and states config order does not matter, and it is the
  same mesh, commit (`c665634-dirty`), and knobs `tuned(17,0.5,21)` — but
  **this must be stated when R7's table is published** (`phase_17` §4 trap list).
- **`bcache_R7.bin` is config-independent** and was deliberately kept through the
  relaunch. It is a speed-only cache and moves no reported metric.
- **The rotor is DIRICHLET, not Neumann** — see `phase_20` §4, still current.
- **`tol_abs = 1.77e-8` is CORRECT** for the Phase 3 winner knob set — see
  `phase_20` §4, still current. Do not "correct" it to 9.93e-8.
- **`niter` is not a homogeneous column** — Krylov iterations, FGS sweeps, and
  outer FGMRES counts share it. `fgs` writes `-1`. Do not "fix" it.

---

## 5. Recommended next actions, in order

> **NEW RULING 2026-09-05 (Ryan), read before any measurement work:** the item is
> **Bernoulli-only** for force recovery from here on — steady `PressureBernoulli` +
> Force. `PressureLaplace` and `KuttaJoukowski` are out of the forward path (pull back
> in for diagnosis only). Retires the `PressureLaplace` CG `itmax` per-rung knob;
> demotes the R7 `CT_laplace` defect to non-blocking (still logged in `ledger.md`).
> Full text: `decision_rules.md`, "Force recovery — Bernoulli-only".

1. **Resolve the two items in §0** — disk sweep authorization, notebook entry.
2. **Harvest Phase 1.** Every rung R1–R7 now has `agreement.csv` +
   `agreement_spread.csv`. This is a complete, publishable deliverable and
   nothing is blocking it. Run it through `phase_17_harvest_prompt.md` §4/§7
   before publishing any table, and carry the R7 split-provenance note.
3. **Phase 3 launch is unblocked for R1–R3** (Step A complete there). Per
   `phase_20` §6.3: needs `PHASE=phase3` (which auto-enables the monitors `CT`
   depends on) and the frozen knobs. On the first real rows, confirm
   `bcerr_certified` is true before trusting any BC verdict, and re-read
   `t_bcerr / t_solve` — the 62 % figure is from a cold-start smoke and should
   fall with a developed wake (`phase_20` §5). **Launching is Ryan's call.**
4. **Watch R4.** It is the long pole now and has already been preempted once.
5. **Still open, deliberately not done:** `benchmark/slurm/p2_unsteady.sh` has
   never been submitted; the whole tree is uncommitted (~100 modified,
   `benchmark/` untracked in full) and both cluster checkouts are rsync copies
   with no git anchor — flagged 08-22, still true, and this session added two
   more uncommitted edits (§3). Ryan's call, separate conversation.

---

## 6. Standing limits (unchanged)

- **No `sbatch` / `scancel` / `scontrol`** without Ryan asking in the moment.
  (The 08-25 cancel/relaunch was explicitly requested.)
- **No sync to `~/flowpanel-021/FLOWPanel.jl`** while jobs run — this is why §3's
  edits are local-only. R4/R5 are still running from that checkout.
- **No git commit.**
- **Notebook entries and BRAINSTORM item edits: offer, don't write.**
