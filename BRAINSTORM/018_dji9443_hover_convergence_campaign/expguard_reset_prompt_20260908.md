# RESET PROMPT — 018: euler_exp sigma guard built, cluster port + relaunch pending (2026-09-08)

You are picking up BRAINSTORM item 018 (DJI 9443 hover CT convergence campaign).
**Entry point: this file.** It supersedes `sfs_cd_reset_prompt_20260908.md`, whose
three tasks are all COMPLETE (results below — do NOT re-derive them). Read the
predecessor only for campaign background; use `brainstorm-scout` for anything
older, never read whole BRAINSTORM files inline.

## Where this stopped

Ryan's instruction was: **"add the sigma floor/cap to the exp path and relaunch."**
The code is written, tested, committed and tagged **locally in both repos**.
The **cluster port has not been applied and nothing has been submitted.**
Ryan then called a break for a context reset. Pick up at "What to do next".

## Standing guardrails (violating these has cost this campaign whole days)

- **sacct state is NOT evidence.** Judge every run by its outputs.
- **Read the `.err` before the `.out`** — Julia stacktraces go only to `.err`.
  (This is how the Ladder C mislabel below was caught.)
- **`ssh orc` needs `bash -lc` + a live ControlMaster** (2FA otherwise). The
  socket drops occasionally; a dropped socket is NOT a hard blocker, just retry.
  numpy IS available on the cluster (`/apps/python/3.12.2`, numpy 1.26.4) —
  a subagent wrongly reported otherwise on 09-08.
- **NODE_FAIL → fresh resubmit, never a chain.**
- **Worktrees + annotated tags only; no silos.** Never edit a worktree that has a
  queued or running job. Campaign worktrees carry no uncommitted tracked state.
- **Another agent shares the local checkout and the queue.** Re-check
  `git status` and `squeue` before acting; never `git checkout .` or `git stash`,
  and stage only your own files when committing.
- Local jobs: never more than 4 threads (HPC exempt).
- **No 018 submissions beyond the three arms Ryan authorised below.**

## Ryan's decisions on 2026-09-08 (binding)

1. **Floor value: `SIGMA_FLOOR_FRAC=0.25`** → floor = 0.25 x tip_sigma_default
   = 0.25 x 0.00476 = **0.00119 m**. Chosen from the sigma low-tail data: it
   clamps ~60 of 180k particles (3e-4) in healthy runs but engages on 2.2% of
   the field in cs0p002 at rev 24.9.
2. **Relaunch all three arms** (he declined the "hold for disk" option):
   - `_3r_cs0p002_exp_nt` NT36, guarded
   - `_3r_cs0p34_exp_nt` NT36, guarded
   - `_3r_exp_nt` NT36, guarded — the control (dynamic Cd, same stack; it
     already ran healthy UNguarded, so this shows whether the guard perturbs a
     healthy run).

## COMPLETE — results from the predecessor's three tasks (do NOT re-derive)

### Task 1: the two rungs are scored
Windowed CT = mean over revs 21-30 (`in_convergence_window == true`), matched
windows, CSV timestamps and job->run-name mapping both verified.

| Job | Arm | NT | Windowed CT | SE | Coverage | vs ref |
|---|---|---|---|---|---|---|
| 13605974 | `_3r_exp_nt` on-model | 72 | **0.071914** | 0.000047 | 10 revs, steps 1441-2160 | +0.10% (ref 0.071844) |
| 13605983 | Ladder S `_3r_cs0p18` | 36 | **0.071099** | 0.000022 | 10 revs, steps 721-1080 | +0.46% (ref 0.070775) |

**exp-nt formal slope**: paired with the on-model NT36 rung 0.070536 (13603853),
climb = **+1.95%** vs the reference family's +1.51%. `WAKE_EXPINT` is NOT a null
on the climb — levels match reference within +-0.5% at both rungs but the climb
is ~30% steeper. **Retire the "provisional null" label** (it came from the
off-model +1.28% pair).

### Task 2: NT has NO accountable impact on Cd in the scoring window
- The Cd formula contains **no dt at all** except through `rlxf`. Source
  docstring (`FLOWVPM_subfilterscale.jl:714`): `rlxf = dt/T`, T = averaging
  window. Fixed rlxf => T halves per NT doubling.
- **rlxf was never NT-scaled**: `sfs_rlxf` defaults to 0.005 in both drivers and
  no launcher sets `SFS_RLXF` at any NT — while the *vortex* relaxation IS
  exact-rate scaled (`RELAX_RLXF` 0.3 / 0.16334 / 0.08539).
- **The fix has already been run**: exact rate `1-(1-0.005)^(36/72) = 1-sqrt(0.995)
  = 0.0025031` is bit-for-bit the `_3r_srlx` value (as `0.16334 = 1-sqrt(0.7)`).
  srlx is the identity at NT36, so its NT72 number is a valid corrected slope:
  climb **+1.51% -> +0.92%**, i.e. Δt-invariant rlxf removes **39%** of the climb.
- **Cd itself converges between rungs.** Clean pair `_3r_exp_nt` (both rungs
  COMPLETE, full VTK, healthy to rev 30), mC ratio NT72/NT36 by rev:
  8 -> 1.243, 12 -> 1.128, 16 -> 1.094, 20 -> 1.035, 24 -> 1.016, 28 -> 0.995.
  **The NT dependence of Cd is a startup transient that decays to 1.00 by rev 28**
  — i.e. it is ~0 across the CT scoring window (revs 21-30). The previously
  quoted "1.2-1.4x" was measured at revs 4-16, inside the transient.
  So **Cd is not the carrier of the CT climb.**
- Extractor validated on the constant-Cs control (13605983): C takes exactly two
  values {0, 0.180000007} (float32), fr0 = 0.525/0.544 — the ~53% zeros are
  `clipping_backscatter` and are **Cs-independent**.
- Matched-rev wake state confirms the ladder design (`NT x P_PER_STEP = 432`):
  median sigma ratio 1.02-1.04, particle counts within ~7%.
- **Caveat**: mean `nume` is negative while mean C is positive, so
  `mean(nume)/mean(deno)` is NOT a valid decomposition of mean Cd. Per-particle
  ratio statistics were not run.
- **Recommended (not yet ruled on by Ryan)**: adopt exact-rate `SFS_RLXF` for NT
  rungs; do NOT carry "Cd NT-drift" in the CT error budget; quote Cd with its rev.

### Task 3: Ladder C dies by sigma COLLAPSE, and the label was wrong
- **The deaths are not force blow-ups.** Both `.err` files show
  `DomainError ... dt*|L| bound ... exceeds euler_exp broadcast substep budget`
  in `_euler_exp_broadcast!`. The numbers quoted through the handoffs as forces
  — **8.6e9 and 8.2e4 — are dt*||grad u||**, not forces. Correct this wherever it
  appears (it also affects the `_3r_nosfs` NT72 death label).
- **sigma never runs away.** p99 flat at ~0.0103-0.0116 in both dead runs right
  to the death step, indistinguishable from the healthy control; the 0.030
  ceiling is never crossed. **|Gamma| is what diverges** (cs0p002 |G|max:
  3.5e-3 -> 6.9e-2 -> 2.2e-1 -> 3.9e7 over steps 673/880/897/898), and in the
  ignition particles **sigma collapses** (0.0090 -> 0.0004 -> ~0).
- **Mechanism, from the code.** The campaign runs the default
  `ReformulatedVPM(f=0, g=1/5)`. On the standard `euler` path the SFS term drops
  out of Z identically and `dsigma/dt = -(1/5) sigma S` exactly with
  `S = ghat . grad u . ghat` — so sigma growth requires S<0 (Ryan's
  anti-stretching prediction, **confirmed**: top-1% sigma is enriched in S<0,
  0.60-0.68 vs a 0.39-0.44 baseline, **including in the healthy control** — so
  anti-stretching is ordinary and bounded, and is NOT what killed Ladder C).
  Ladder C runs `euler_exp`, whose geometric step is different: with
  `q = exp(dt L) G` and `r = |q|/|G|`, it applies `|G| ~ r^(2g)`,
  `sigma ~ r^(-g)`, so **`|G| sigma^2` is an exact invariant**. Gamma runaway and
  sigma collapse are the same event. sigma->0 makes induced gradients diverge,
  which raises r for neighbours: a closed feedback loop. `euler_exp` rejected
  `sigma_guard`, so there was **no floor** to break it.
- **The FMM-adequacy deaths are the opposite end of the same invariant** (sigma
  ceiling). Two distinct failure modes.
- **Consequence for 026**: particle splitting attacks *oversize* sigma, so it
  does NOT address the Ladder C / `_3r_nosfs` mode. This cuts against the
  standing "splitting is the remedy" plan for that half. Ryan has not ruled.
- **Ladder C is de-confounded — but not by `cs0p18`** (which differs in BOTH
  guard and integrator). The right control is **`_3r_exp_nt`**: same `euler_exp`,
  same unguarded status, healthy to rev 30, differing only in dynamic vs frozen
  Cd. On that comparison static Cs is the distinguishing factor (n=2).
- **A sigma floor will probably NOT save cs0p34.** One step before death its
  sigma low tail is completely healthy (min 0.00164, **zero** particles below
  0.25*sigma0) while |G|max goes 0.0045 -> 0.805 -> 2.71 -> 831 in six steps.
  Only cs0p002 has a slow sigma-collapse precursor (2.2% of the field below
  0.25*sigma0 by rev 24.9 vs 0.03% healthy). The remaining Gamma source in
  `euler_exp` is the additive SFS Lie split
  `G -= dt*C*SFS*sigma^3/zeta0`, which is not norm-preserving and which no sigma
  guard touches; Cs=0.34 is ~2.2x the dynamic mean (0.158). Ryan was told this
  and chose to run cs0p34 anyway — if it dies identically, that confirms the
  SFS-additive route.

Sigma low-tail reference numbers (sigma0 shed = 0.00476 m):

| run | step | rev | sigma min | p0.1 | p1 | n < 0.25*sigma0 |
|---|---|---|---|---|---|---|
| cs0p18 (healthy) | 1079 | 29.97 | 0.00070 | 0.00140 | 0.00191 | 61 |
| exp_nt36 (healthy) | 1079 | 29.97 | 0.00065 | 0.00144 | 0.00190 | 67 |
| cs0p002 (dies 898) | 897 | 24.92 | 0.00003 | 0.00060 | 0.00100 | 4401 |
| cs0p34 (dies 361) | 360 | 10.00 | 0.00164 | 0.00235 | 0.00272 | 0 |

## DONE — code, tests, commits, tags (LOCAL ONLY)

**Design**: the guard clamps the geometric **gain ratio r**, not sigma post-hoc,
so Gamma, sigma and M[9] all derive from the clamped r and `|G| sigma^2` stays
exact. Clamping sigma alone would floor the core size while leaving the
circulation amplification unguarded — which is the actual failure. Bounds:
`ceil => r >= (sigma/ceil)^(1/g)`; `floor => r <= (sigma/floor)^(1/g)`;
`dtz_cap => r <= exp(dtz_cap/g)` (since `dt*Z = g*log(r)` on this path).
An empty `sigma_guard` still reproduces the unguarded step **bit-exactly** (the
unguarded expression is kept on its own branch, not reassociated).

| repo | branch | commit | tag |
|---|---|---|---|
| FLOWVPM.jl | `flowpanel` | **21eeaaa** | `campaign/p018-expguard-20260908` |
| FLOWPanel.jl | `fastmultipole` | **7dee1ab** | `campaign/p018-expguard-20260908` |

Files changed: `FLOWVPM/src/FLOWVPM_timeintegration.jl` (guard on both the scalar
and broadcast/GPU euler_exp paths + new `_exp_ratio_bounds` helper),
`FLOWVPM/test/runtests_expint.jl` (+16 tests), `FLOWVPM/test/runtests_gpu_fmm.jl`
(guarded CPU-vs-broadcast parity), `FLOWPanel/src/FLOWPanel_wake.jl` (forward the
guard; rungekutta3 still rejects it), `FLOWPanel/examples/rotor_hover_ground_effect.jl`
(same removal in that file's local propagate copy, so 022 does not diverge).

**Tests, all green**: expint 28 pre-existing + 16 new; gpu_fmm euler_exp
loop-vs-broadcast equivalence 411/411 with the guard active (bounds respected,
invariant conserved); dsigma2 accumulators 44/44; FLOWPanel `runtests_unit_wake.jl`
(714 free-wake tests, expint kwarg testset) and `runtests_unit_simulate.jl`.
No test anywhere asserted the old rejection.

Driver defaults are unchanged: `SIGMA_FLOOR_FRAC` / `SIGMA_CEIL` / `SIGMA_DTZ_CAP`
still default to off, so an unset environment reproduces prior behaviour exactly.

## NOT DONE — this is where you start

### The cluster is on DIVERGENT branches — do not just push
- Cluster `~/projects/FLOWVPM.jl` is on branch **`unified-052` @ 3315b22**, and
  **lacks the local 026 resolution-split / dsigma2 accumulator lines**. Cluster
  `~/projects/FLOWPanel.jl` is at **4e6b5b7** and has no rungekutta3 branch in
  `propagate!`. Both remotes are public GitHub (`byuflowlab`) — pushing there is
  outward-facing; **do not push without asking Ryan.**
- **Verified**: the euler_exp *arithmetic* is identical on cluster and local
  (the broadcast block is byte-identical; the CPU block differs only by the 026
  lines). So the guard ports cleanly.
- A ready port script is already uploaded at **`~/port_guard.py`** on the cluster
  (source of truth also in this session's scratchpad as `port_guard.py`). It
  takes two args — the FLOWVPM `src/FLOWVPM_timeintegration.jl` path and the
  FLOWPanel `src/FLOWPanel_wake.jl` path — applies the cluster-variant anchors,
  and prints `PORT OK`. Every replacement is asserted, so it fails loudly rather
  than half-applying. **It has NOT been run.**

### Steps remaining
1. Create fresh branches + worktrees on the cluster so **no existing worktree is
   touched** (wt052/campaign-052 and wt026 share these repos; 021 and 022 jobs
   are live in the queue). Suggested:
   `git -C ~/projects/FLOWVPM.jl branch p018-expguard-20260908 3315b22` then
   `git -C ~/projects/FLOWVPM.jl worktree add ~/wt018/FLOWVPM-expguard p018-expguard-20260908`
   (and the same for FLOWPanel from `4e6b5b7` -> `~/wt018/FLOWPanel-expguard`).
2. Run `~/port_guard.py` against the files **inside the new worktrees**, commit
   there, and create annotated tags (`campaign/p018-expguard-20260908` on the
   cluster side too) so the pins are citable.
3. Point the campaign Julia environment's Manifest dev-paths at the new
   worktrees (not at `~/projects/*`).
4. Smoke-test on the cluster before submitting — at minimum load the package and
   run a few steps with the guard armed, since the ported CPU block was NOT
   compiled locally in its cluster form.
5. Submit the three arms with `SIGMA_FLOOR_FRAC=0.25` and `SIGMA_CEIL=0.030`
   (Ryan said "floor/cap"; 0.030 is the wave-stack ceiling). **Confirm the
   launcher case names** in `examples/run_dji9443_hover_ct_hpc.slurm.sh` — the
   run names are `p018_csarc_l3p0_3r_cs0p002_exp_nt`,
   `p018_csarc_l3p0_3r_cs0p34_exp_nt`, `p018_csarc_l3p0_3r_exp_nt`, but the case
   labels were not verified this session. Give the guarded reruns **new run
   names** so they do not clobber the existing data dirs (the base-case dirs were
   already clobbered once, on 09-07).
6. Record the pins (tag + SHA) in a provenance file before submitting.

## Open items for Ryan (do not act unprompted)
- **DISK: `/home/rander39` is at 548 G against the 400 G cap** — up from the
  430 G quoted in the previous handoff. Filesystem quota is 2.0 T (27% used), so
  submissions will not fail, but this is 148 G over the project's own cap. A
  strip belongs to Ryan's OTHER session — **do not start a second archiver
  without checking with him.**
- Whether to adopt exact-rate `SFS_RLXF` as campaign policy (Task 2).
- Whether 026 splitting still counts as the remedy now that it only addresses the
  ceiling half of the failure (Task 3).
- WGE-arm reruns and the `_3r_sv` resubmit (still unblocked, still unruled).
- Where a merge-result sigma cap lands in FLOWVPM (pair merges uncapped at x1.26).
- Recovering the clobbered base-case reference CSVs (0.070775 / 0.071844) from
  the archive tarballs — still owed.

## Notebook
`~/Dropbox/research/notebooks/journals/20260901.md` has a `# 20260908` header.
**Added this session** (approved by Ryan, table only, no prose):
`## 018 NT ladders — full attempt inventory` — 18 arms with knob, NT36/NT72/NT144,
climb and status.
**Still owed** for the 09-07/09-08 arc: the Cd transient finding, the rlxf
derivation, and the Ladder C sigma-collapse forensics. **ASK Ryan for verbosity
before writing**; append-only, under the existing `# 20260908` header.

## Reusable tooling (cluster)
- `~/vtp_C.py` — `read_vtp(path, want)` parses the raw-appended VTK XML directly
  (**meshio cannot read these files**). Note the wake VTPs contain **no `static`
  array**, so all particles in them are active wake particles.
- `~/cd_nt.py` — Cd-vs-NT tables + constant-Cs extractor validation.
- `~/sig.py` — sigma/|Gamma|/stretching-S percentiles and tail enrichment.
- `~/loc.py` — spatial + age localisation of the high-|Gamma| population.
- `~/port_guard.py` — the cluster port described above.
