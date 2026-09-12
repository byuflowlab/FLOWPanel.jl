# RESET PROMPT — 026 splitting launch prep + 018 expguard results (2026-09-11)

**Entry point: this file.** It carries (a) the finished 018 `expguard` campaign
and its verdict, and (b) a half-finished 026 launch prep whose numeric half is
**STALE by design** — the splitting syntax and options changed after it was
written. Read `expguard_reset_prompt_20260908.md` (in the 018 item dir) only for
campaign background. Use `brainstorm-scout` for anything older; never read a
whole BRAINSTORM file inline.

Ryan's last instruction: *"the splitting syntax and options just changed. Prep
for a context reset."* He was mid-way through clarifying a set of launch
questions (below) and has **not ruled on any of them**.

---

## PART 0 — Standing guardrails (violating these has cost this campaign whole days)

- **sacct state is NOT evidence.** Judge every run by its outputs.
- **Read the `.err` before the `.out`** — Julia stacktraces go only to `.err`.
- **`ssh orc` needs `bash -lc` + a live ControlMaster** (2FA otherwise). A
  dropped socket is not a hard blocker; retry. numpy IS available on the
  cluster (`python3`, numpy 1.26.4).
- **NODE_FAIL → fresh resubmit, never a chain.**
- **Worktrees + annotated tags only; no silos.** Never edit a worktree with a
  queued or running job. Campaign worktrees carry no uncommitted tracked state.
- **Another agent shares the local checkout and the queue** — and is actively
  editing the splitting code right now (see Part 2). Re-check `git status` and
  `squeue` before acting; never `git checkout .` or `git stash`; stage only
  your own files when committing.
- Local jobs: never more than 4 threads (HPC exempt).
- **No 026/018 submissions Ryan has not authorised.** The three 018 expguard
  arms were authorised and are DONE. **Nothing for 026 is authorised yet.**

---

## PART 1 — DONE: 018 expguard campaign (do NOT re-derive)

### Pins (all tagged `campaign/p018-expguard-20260908`)

| repo | cluster worktree | commit | base |
|---|---|---|---|
| FLOWVPM.jl | `~/wt018/FLOWVPM-expguard` | `7468712` | `3315b22` (`unified-052`) |
| FLOWPanel.jl | `~/wt018/FLOWPanel-expguard` | `d5dd772` | `67b7383` (`p018cs-fp-nt-20260908`) |
| FastMultipole | `~/wt018/FastMultipole-expguard` | `3da58a1a` | unmodified |

Local counterparts: FLOWVPM `21eeaaa` (branch `flowpanel`), FLOWPanel `7dee1ab`
(branch `fastmultipole`), same tag name. Julia env:
`~/p018wtenv-expguard-gh200` (copy of `~/p018wtenv-nt-cs-gh200`, all three
dev-paths repointed). Full provenance:
`BRAINSTORM/018_dji9443_hover_convergence_campaign/expguard_provenance_20260908.md`.

**Base-choice note:** the 09-08 handoff said to base FLOWPanel on `4e6b5b7`
(`unified-052`). That was wrong — `sacct` WorkDir shows the arms ran from
`67b7383`, which differs from `4e6b5b7` in the driver *and* 150 lines of
`src/FLOWPanel_wake.jl`. The worktree was torn down and rebuilt on `67b7383`
(a strict superset of the control's `d2d3bb7`). Smoke job 13610775 passed
16/16 guard unit tests on the ported code, CUDA functional.

### Results — all three arms terminal

Guard armed at `SIGMA_FLOOR_FRAC=0.25` (floor 0.00119 m) + `SIGMA_CEIL=0.030`.
Windowed CT = mean over revs 21–30 (`in_convergence_window == true`).

| job | arm | outcome | windowed CT |
|---|---|---|---|
| 13610777 | `_3r_exp_nt_g25` control, dynamic Cd | COMPLETED 1079/1079 | **0.070808** |
| 13610778 | `_3r_cs0p002_exp_nt_g25` | COMPLETED 1079/1079 | **0.069317** |
| 13610779 | `_3r_cs0p34_exp_nt_g25` | **DIED step ~354**, `dt*|L| = 92786.7` | — |

Unguarded predecessors: control 13603853 windowed CT **0.070536** (COMPLETED);
cs0p002 13605985 died step 898 (rev 24.9); cs0p34 13605986 died step 361.

### Verdict (robust, three-arm, n=1 per arm but the contrast is clean)

1. **The σ floor arrests the geometric-collapse channel.** cs0p002 went from
   dying at step 898 to completing all 1079 steps with a scoreable CT.
2. **The σ floor does nothing to the SFS-additive channel.** cs0p34 died at
   step ~354 vs 361 unguarded — same signature, ~2% earlier (noise). Confirms
   the prediction made from the code before the run.
3. **The guard is NOT a null on a healthy run.** Control moved
   0.070536 → 0.070808 = **+0.39%**. The NT climb under study is +1.95%, so
   the guard is worth ~20% of the effect size. **This must either enter the CT
   error budget or drive the floor down.** Ryan has not ruled.

### The mechanism (established from the code, confirmed by the runs)

On `euler_exp` with `ReformulatedVPM(f=0, g=1/5)` the step is a Lie split:

- **Geometric piece.** With frozen gradient $L$, $q = e^{\Delta t L}\Gamma_0$,
  gain ratio $r = |q|/|\Gamma_0|$: $\Gamma \leftarrow q\,r^{-3g}$ and
  $\sigma \leftarrow \sigma r^{-g}$. At $g=1/5$, $1-5g=0$, so
  $|\Gamma|\sigma^2$ is an **exact invariant**. Γ runaway and σ collapse are
  the same event.
- **SFS piece**, applied *after*: $\Gamma \mathrel{-}= \Delta t\,C\,\mathrm{SFS}\,\sigma^3/\zeta_0$.
  Additive, not norm-preserving — the **only** operation in the integrator that
  can change $P = |\Gamma|\sigma^2$. No σ guard touches it.

**Ignition loop:** a kick raising $|\Gamma|$ at fixed σ raises $P$ permanently
(the geometric step conserves it thereafter) → larger induced gradients
($\nabla u \sim \Gamma/\sigma^4$) → larger $r$ → but $r>1$ *forces*
$\sigma \leftarrow \sigma r^{-g}$ → smaller σ amplifies $\nabla u$ again.
Closed feedback, finite-time blow-up, terminating at the `dt*|L|` substep-budget
`DomainError` in `_euler_exp_broadcast!`.

**Why SFS is the igniter, not the sustainer:** the kick scales as $\sigma^3$, so
once σ collapses it shuts itself off and the geometric loop runs alone. Hence
cs0p34 ($C_s=0.34 \approx 2.2\times$ the dynamic mean 0.158) dies with **no σ
precursor**, while cs0p002 (kick ≈ 0) takes the slow geometric route with a
2.2%-of-field σ-collapse precursor. The guard clamps $r$ (not σ post-hoc), so it
breaks the geometric link only — exactly matching the observed outcomes.

**Label correction carried forward:** the numbers 8.6e9 / 8.2e4 / 92786.7 quoted
through these handoffs as "forces" are **dt·‖∇u‖**, not forces.

---

## PART 2 — STALE: everything below depends on the splitting API

**Ryan: "the splitting syntax and options just changed."** As of this writing
the following carry **uncommitted** modifications from the other agent:

```
FLOWVPM.jl:   src/FLOWVPM_resolution_split.jl
              src/FLOWVPM_timeintegration.jl
              test/runtests_resolution_split.jl
              examples/p026_ring_split_test.jl
FLOWPanel.jl: examples/rotor_hover_pressure_comparison.jl
              examples/run_p018_screen_hpc.slurm.sh
              BRAINSTORM/026_sigma_growth_particle_splitting/particle_splitting_design.md
```

FLOWVPM HEAD `8b0b70d`, FLOWPanel HEAD `8dce66c` — but **the working tree is
ahead of both**. Treat every claim in this Part as a hypothesis to re-verify.

### YOUR FIRST TASK: re-read the API, then re-derive

Read, in this order:
1. `FLOWVPM.jl/src/FLOWVPM_resolution_split.jl` — the options struct, the
   trigger predicate(s), and the three emission kernels (child count, child σ,
   child Γ for each mechanism).
2. The `# --- BRAINSTORM 026` block in
   `FLOWPanel.jl/examples/rotor_hover_pressure_comparison.jl` (was ~line 710–785)
   — the env-knob surface and what defaults to what.
3. The `scr_p026sp_*` / `scr_p026s9_*` arm table in
   `examples/run_p018_screen_hpc.slurm.sh` (was ~line 245–280).
4. The newest `# §` section of `particle_splitting_design.md` (was §19).

Then re-derive the two cost models below against whatever the API now is.

### What the API looked like BEFORE the change (for diffing only — do not trust)

`ResolutionSplitOpts` fields: `f_visc`, `f_comp`, `f_elong` (per-mechanism
growth fractions vs the particle's own `sigma_0`, `NaN` disables each);
`sigma_min`/`sigma_max` (emission **clamps**, not triggers, defaulting to the
052c `sigma_guard` floor and `SIGMA_CEIL`); `enable_viscous_split`,
`enable_stretch_split`; `viscous_offset_ratio` 1.3503, `compress_offset_ratio`
0.6, `elongate_offset_ratio` 0.5; `use_stretch_axis`, `axis_coherence_min`.
Driver env: `WAKE_SPLIT_VISCOUS` / `_STRETCH` / `_FRAC_VISCOUS` /
`_FRAC_COMPRESS` / `_FRAC_ELONGATE` / `_SIGMA_MIN` / `_SIGMA_MAX` /
`_STRETCH_AXIS` / `_*_OFFSET_RATIO` / `_EVERY` / `_VERBOSE`.

Geometry as read on 09-09:

| mechanism | fires when | children | σ_child | Γ_child |
|---|---|---|---|---|
| viscous `tetra4` | $\sqrt{\sigma_0^2+\Delta_{visc}}/\sigma_0 > 1+f_{visc}$ | 4 | $0.630\,\sigma_p$ | $\Gamma_p/4$ |
| compress `tri3` | $\sqrt{\sigma_0^2+\Delta_{rvpm}}/\sigma_0 > 1+f_{comp}$ | 3 | $\sigma_p/\sqrt3$ | $\Gamma_p/3$ |
| elongate `pair2` | $\sqrt{\sigma_0^2+\Delta_{rvpm}}/\sigma_0 < 1-f_{elong}$ | 2 | $\sigma_p$ **unchanged** | $\Gamma_p/2$ |

### Two findings from that read worth re-testing — they change the item's plan

**(a) Splitting now attacks BOTH halves of the invariant.** The `pair2`
elongation regime (Ryan's 09-07 ruling) halves $|\Gamma|$ per particle at
*fixed* σ, which is a direct attack on the $|\Gamma|\sigma^2$ runaway of Part 1.
**This retires the standing note that "026 splitting only addresses oversize σ
and therefore cannot help the Ladder C / `_3r_nosfs` mode."** If the new API
keeps an elongation regime with σ_child ≈ σ_parent, the retraction stands.
Ryan has not been told this yet — it is worth raising, because it reopens
"is splitting the remedy?" for the collapse half.

**(b) Splitting may be CPU-only, silently.** The Δσ² accumulators were updated
only in the *scalar CPU* integrator paths (`_euler_cpu_reformulated!`,
`_euler_exp_cpu!`, and the rk3 CPU path). Neither `_euler_broadcast_reformulated!`
nor `_euler_exp_broadcast!` touched them. On a device-resident wake
(`VPM_ARRAYTYPE=cuarray`, which the GPU launcher sets by default) the
accumulators would stay zero forever and **no trigger would ever fire, with no
error** — `_heal_unseeded_rsplit_slots!` repairs `sigma_0` on the host mirror
but not the accumulators. The old design doc said the same ("device-resident
integrator twins still do not maintain `dvisc`/`drvpm`; splitting remains
host-mirror-only").
**Re-verify against the edited `FLOWVPM_timeintegration.jl`.** If it still
holds: these arms must run CPU (the `run_p018_screen_hpc.slurm.sh` dispatcher,
not the GPU launcher), this cuts against Ryan's standing GPU-default
preference, and a fail-fast guard is worth proposing.

### The cost models (method is API-independent; re-run the arithmetic)

Both rest on one property: **`sigma_0` resets to the child σ on every split**,
so repeated firing is a *multiplicative ladder*, not a linear one.

**Elongation (`pair2`, σ unchanged, ×2/fire).** A particle descending to the
0.25σ₀ floor fires $n = \ln(0.25)/\ln(1-f_{elong})$ times, costing $2^n$:

| $f_{elong}$ | fires | multiplier | healthy field (60 of 180k) | ignited tail (4401 of 180k) |
|---|---|---|---|---|
| 0.2 | 6.2 | 74× | +4.4k | +325k (≈2× field) |
| 0.3 | 3.9 | 15× | +0.9k (+0.5%) | +66k (+37%) |
| 0.4 | 2.7 | 6.5× | +0.4k | +29k (+16%) |
| 0.5 | 2.0 | 4× | +0.2k | +18k (+10%) |

Population counts are measured, from the σ low-tail census: healthy runs put
~60–67 of 180k particles below 0.25σ₀ at rev 30; cs0p002 put 4401 there by
rev 24.9. σ₀(shed) = 0.0047597 m.

**Compression (`tri3`, σ_c = σ_p/√3, ×3/fire).** After a fire σ becomes
$(1+f_{comp})\sigma_p/\sqrt3$, a *net reduction* iff
$$f_{comp} < \sqrt3 - 1 \approx 0.732$$
Below that bound each fire more than undoes the growth that triggered it and the
σ_max emission clamp becomes effectively unreachable — which was §19's stated
goal ("splits fire before particles sit long at the clamp"). **Re-derive this
bound from the new σ_child rule**; it is the single most useful number for
picking $f_{comp}$ and it falls straight out of the geometry.

Note the bulk σ growth in these runs is **core spreading** (`CORE_SPREADING_ACTIVE=true`,
`WAKE_CORE_BETA=1e9`), which landed in the `dvisc` accumulator and stayed
unsplit because every dispatcher arm set only `WAKE_SPLIT_STRETCH=true`. That
is what keeps the cost numbers small. If `f_visc` is enabled the arithmetic
changes a lot.

### Values I was about to recommend (NOT ruled on, and API-dependent)

`f_comp = 0.5` (net 0.87× per fire, inside the self-limiting bound, fires on the
~1% anti-stretching tail) and `f_elong = 0.3` (acts at 0.7σ, well before the
0.25 floor engages, ~free on a healthy run). Treat as a starting point only.

---

## PART 3 — The launch that did NOT happen

"The corresponding jobs" was read as the **commit-7 arms**, still Ryan-gated:

- **§9 matrix, 12 arms**: `scr_p026s9_{ctrl,exp,ctrllg,explg}_{floor,split,fs}`
  — shrink-side validation. `floor` arms need `SIGMA_FLOOR_FRAC`; `split` arms
  need the elongation fraction; `fs` arms need both. Γ̂-fallback comparison arms
  add `WAKE_SPLIT_STRETCH_AXIS=false`.
- **§8.4 cap arms, 2**: `scr_p026sp_nt144_cap030` / `_cap018` at NT144 —
  ceiling side, need the compression fraction; their `WAKE_SPLIT_SIGMA_MAX`
  (0.030 / 0.018) is now an emission clamp, not a trigger.

Both families live in `examples/run_p018_screen_hpc.slurm.sh` (**modified —
re-read**). All still carry `<re-derive>` placeholders in their comments.

### Open questions for Ryan (he asked to clarify these before answering)

1. **Which arms?** All 14, the 12 §9 arms first, or a 3-arm de-risk
   (`exp_{floor,split,fs}`) to confirm the triggers actually fire at the chosen
   fractions before committing the matrix.
2. **Cold vs warm-start.** The §9 arms have s020v warm-start VTPs bracketing
   ignition; everything above was costed as cold 30-rev runs. Warm-starting
   near ignition would make the fraction choice far cheaper to iterate and
   changes how much needs settling up front.
3. **Does `f_visc` stay off?** Every dispatcher arm enables only the stretch
   mechanism. Enabling viscous splitting changes the cost model substantially.
4. **How hard should the +0.39% guard perturbation drive the floor?** Either it
   is a reason to drop the floor below 0.25, or it simply enters the CT error
   budget and 0.25 stands. This also decides whether a floor-sensitivity pair
   (0.25 vs 0.15) is worth 4 extra arms.
5. **The two fractions themselves**, once re-derived against the new API.

Ryan had rejected the question set to add clarification and had not yet said
what he wanted clarified. **Ask him that first** — do not re-ask the old
questions verbatim.

---

## PART 4 — Other open items (unchanged, do not act unprompted)

- **DISK: clear as of 2026-09-11 — `/home/rander39` is at 191 G** against the
  400 G cap (it was 548 G on 09-08; Ryan's other session evidently archived).
  There is headroom for a 12–14 arm CPU launch, but re-measure before
  submitting and do not start a second archiver without checking with him.
- **Queue as of 2026-09-11**: no 018 or 026 jobs. Live: 13653115
  `p021-cold-pilot`, 13593021 `p2lg-tune-R7` (both item 021, not yours).
- Whether to adopt exact-rate `SFS_RLXF` for NT rungs (removes 39% of the CT
  climb: +1.51% → +0.92%).
- WGE-arm reruns and the `_3r_sv` resubmit — unblocked, unruled.
- Where a merge-result σ cap lands in FLOWVPM (pair merges uncapped at ×1.26).
- Recovering the clobbered base-case reference CSVs (0.070775 / 0.071844) from
  the archive tarballs — still owed.
- `WAKE_EXPINT` is **not** a null on the NT climb (+1.95% vs the reference
  family's +1.51%); the "provisional null" label is retired.
- Cd is **not** the carrier of the CT climb — its NT dependence is a startup
  transient that decays to 1.00 by rev 28, i.e. ~0 across the scoring window.

## Notebook (append-only, ASK before writing, ASK how verbose)

`~/Dropbox/research/notebooks/journals/20260901.md` has a `# 20260908` header
with `## 018 NT ladders — full attempt inventory` already written. **Still
owed**: the Cd transient finding, the rlxf derivation, the Ladder C
σ-collapse forensics, and now the expguard three-arm result (Part 1) — which is
the cleanest single result of the arc and should probably lead.

## Reusable tooling (cluster)

`~/vtp_C.py` (`read_vtp(path, want)` — meshio CANNOT read these files; wake VTPs
carry no `static` array so all particles in them are active), `~/cd_nt.py`,
`~/sig.py` (σ/|Γ|/stretching-S percentiles + tail enrichment), `~/loc.py`,
`~/port_guard.py`, `~/expguard_smoke.jl` + `~/expguard_smoke.slurm.sh`.
