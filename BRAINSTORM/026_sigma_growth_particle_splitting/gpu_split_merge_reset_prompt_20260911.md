# RESET PROMPT — 026 GPU splitting/merging implementation (2026-09-11)

**Entry point: this file.** Ryan directed (2026-09-11): *implement the
splitting and merging approaches on GPU*, as the prerequisite for the gated
026 campaign arms (his standing preference is GPU-default for long
particle-heavy runs). You are picking up in
`/Users/ryan/Dropbox/research/projects/FLOWPanel.jl` (panel solver) +
`/Users/ryan/Dropbox/research/projects/FLOWVPM.jl` (VPM core). Read both
repos' `CLAUDE.md` + the routing policies before touching code. Design/history
record: `particle_splitting_design.md` §19–§21 (append-only); theory:
`splitting_theory.md` (single living draft, rewrite in place). Launch
context if needed: `split_launch_reset_prompt_20260911.md` (same dir).
Use `brainstorm-scout` for anything older; never read a whole BRAINSTORM
file inline.

## Standing guardrails

- **Local machine is macOS — no CUDA.** All GPU execution/verification goes
  through the cluster (`ssh orc`, needs `bash -lc` + live ControlMaster
  socket; 2FA otherwise; retry a dropped socket). Local: CPU-only dev +
  unit tests, ≤ 4 threads.
- **Another agent may share the local checkouts** (the adaptive-elongation
  session had uncommitted layers on BOTH repos and pending commit
  sequencing — see `adaptive_elongation_reset_prompt_20260911.md`). Check
  `git status` before acting; never `git checkout .` / `git stash`; stage
  only your own files; coordinate commit boundaries with Ryan.
- **No 026 campaign submissions** — launches stay Ryan-gated. HPC smokes of
  the GPU code path are fine (they are verification, not campaign arms),
  but say what you're submitting.
- sacct state is not evidence; read the `.err` before the `.out`.
- Worktrees + annotated tags for anything campaign-grade; no silos.

## Ruled context (do not re-litigate)

- §21 (design doc): arm settings ADOPTED — `f_elong=0.3`, `f_comp=0.73`,
  adaptive elongation at default Φ_t (= shedding `OVERLAP`),
  `MERGE_OVERLAP=3.5` vs production-absolute as the first merge A/B.
- §19/§20.7: per-mechanism fractional triggers on ATTEMPTED (pre-clamp)
  Δσ² accumulators; σ bounds are clamps; adaptive in-line elongation kernel
  (m = clamp(round(λ_att·σ₀/σ_c), 2, m_max), σ_c = σ_p, Γ/m); sigma_guard
  is a permanent partner of splitting.
- Splitting attacks both halves of the |Γ|σ² invariant (elongation divides
  Γ at fixed σ) — it is a live candidate remedy for the geometric-collapse
  channel (018 expguard verdict, `split_launch_reset_prompt_20260911.md`
  Part 1).

## THE GAP (verified 2026-09-11 against the live trees)

Storage-type dispatch in `FLOWVPM.jl/src/FLOWVPM_timeintegration.jl`
(`pfield.particles isa Array` → scalar CPU; else broadcast/GPU — `_euler`
~l.246, `_euler_exp` ~l.457). **Every `_rsplit_accumulate!` /
`_rsplit_accumulate_dsigma2!` call sits in a CPU scalar branch; the
broadcast twins have none and there is no fail-fast**, so on
`VPM_ARRAYTYPE=cuarray` the accumulators stay zero and no trigger ever
fires, silently.

CPU accumulation sites to twin (line numbers from the 09-11 tree — re-check):

| site | file:line | what it accumulates |
|---|---|---|
| `_euler_cpu_reformulated!` | timeintegration.jl:319, 347 | axis sample (dt·S); rVPM Δσ² → `drvpm` |
| `_euler_exp_cpu!` | timeintegration.jl:549, 578 | axis sample; rVPM Δσ² (pre-clamp — note `_exp_ratio_bounds` clamping happens AFTER; capture the unclamped ratio) |
| `update_particle_states_cpu_reformulated!` (rk3) | timeintegration.jl:1087, 1119 | axis sample (rsplit_sample gate); rVPM Δσ² |
| CoreSpreading euler | FLOWVPM_viscous.jl:167 | `dvisc += 2ν·dt` |
| CoreSpreading euler_exp blended | FLOWVPM_viscous.jl:195 | `dvisc += diffusion` |
| CoreSpreading rk3 | FLOWVPM_viscous.jl:216 | `dvisc += aux2·M[7]` |

Broadcast twins needing the additions: `_euler_broadcast_reformulated!`
(:354), `_euler_exp_broadcast!` (:631),
`update_particle_states_broadcast_reformulated!` (:1126);
`_corespreading_{euler,eulerexp,rk3}_broadcast!` (viscous.jl:281/302/320).
GPU ext: `FLOWVPM.jl/ext/FLOWVPMCUDAExt.jl`.

## What ALREADY works on GPU (don't rebuild it)

`FLOWPanel.jl/src/FLOWPanel_gpu_wake.jl` (052 stage A seam): on a
device-backed wake, `_apply_particle_maintenance_device!` (~l.130) does
D2H → host maintenance → H2D on a cached host-mirror ParticleField
(`_gpu_sync_mirror_from_device!` / `_gpu_sync_device_from_mirror!`,
contiguous live-prefix `copyto!`; `_gpu_copy_side_buffers!` currently syncs
only the filament edge graph). **`MergeParticles` AND the `ResolutionSplit`
split pass already run there** (policies:
`FLOWPanel_wake.jl:1630`/`1659`, applied :1732/:1747;
`_heal_unseeded_rsplit_slots!` :1769 seeds `sigma_0` for device-shed
particles). The `MERGE_OVERLAP` knob is a MergeParticles config
(`sigma_relative=true, r=1/Φ_merge`) — it should already work on GPU via
the mirror; **verify, don't reimplement**.

The deliberate design note at `FLOWPanel_gpu_wake.jl:62–75` says
`resolution_split` state is canonical on the MIRROR and the device field
keeps `resolution_split === nothing`. That design is what your work
replaces/amends: with fractional gating, host-only state means zero
accumulation and dead triggers.

## Recommended shape (starting point — you own the design)

1. **Accumulate on device, inline in the broadcast twins.** The state is 7
   scalars/particle (`sigma_0`, `axis`(3), `weight`, `dvisc`, `drvpm`).
   Every update is elementwise-expressible: the sign-aligned axis add is
   `sgn = ifelse.(dot-row .< 0, -1, 1)` per column; Δσ² attributions are
   row broadcasts. Make `ResolutionSplitState` backing arrays match the
   pfield's array type (parameterize the struct or allocate via
   `similar(pfield.particles, ...)`), so the same broadcast code serves
   both backends. Accumulate PRE-CLAMP values (in `_euler_exp_broadcast!`
   compute the unclamped r-derived Δσ² before applying the guard bounds).
2. **Extend the maintenance seam sync.** Add the rs arrays (live prefix
   only) to `_gpu_sync_mirror_from_device!`/`_gpu_sync_device_from_mirror!`
   so the mirror's split pass sees real accumulators and children's fresh
   state flows back. State becomes canonical on the DEVICE field; drop the
   mirror-canonical special case and update the l.62–75 comment.
3. **Keep the split pass itself on the host mirror** (serial loop, RNG
   draws, add_particle) — that seam is proven and split events are rare.
   Device-side `add_particle` shedding still misses the init hook →
   `_heal_unseeded_rsplit_slots!` must now also ensure healed slots have
   zeroed accumulators (it only seeds `sigma_0` today; verify what device
   shedding leaves in the rs arrays).
4. A simpler fallback if (1) fights the type system: a seam-owned bundle of
   device arrays updated by the ext, synced into the host
   ResolutionSplitState at each maintenance pass. Fine too — pick what's
   maintainable.

Watch: `remove_particle` swap-with-last semantics — where do removals
happen on device-backed runs (mirror only, or device too)? The rs arrays
must mirror whatever swaps the particle matrix undergoes on the side where
they live. Also `_GPU_PFIELD_MIRRORS` is keyed by `objectid` — wake
rebuilds call `clear_gpu_pfield_mirrors!`.

## Verification

- **Local (CPU)**: `FLOWVPM.jl/test/runtests_resolution_split.jl` must stay
  green (1011/1011 as of 09-11, includes t10 adaptive + t11 circulation).
  FLOWPanel `runtests_unit_wake.jl` (730/730), `runtests_unit_replay.jl`
  (148/148). Known pre-existing unrelated failure: first testset of
  `runtests_unit_warmstart.jl` (WeakKeyDict/immutable) — not yours.
- **New parity tests**: CPU-vs-broadcast accumulator parity on an
  Array-backed field (the broadcast path runs on plain Arrays too — you can
  test the new code paths WITHOUT a GPU by calling the broadcast twins
  directly), N steps, all three integrators + CoreSpreading variants;
  trigger-decision parity; split-event count parity (RNG makes geometry
  differ — compare counts/σ/Γ-totals, not positions).
- **Cluster (CUDA)**: the dispatcher already has
  `scr_p026ph1b_expgpu_smoke` (exercises `_euler_exp_broadcast!` +
  `_corespreading_eulerexp_broadcast!` end-to-end,
  `examples/run_p018_screen_hpc.slurm.sh` ~l.230) — extend/clone it with
  `WAKE_SPLIT_STRETCH=true WAKE_SPLIT_FRAC_COMPRESS=0.73
  WAKE_SPLIT_FRAC_ELONGATE=0.3` and confirm: nonzero split counters in the
  log, sane σ census, and a companion `MERGE_OVERLAP=3.5` smoke. Compare
  split-event rates against a CPU run of the same arm (backend-matched
  discriminator hygiene).
- Acceptance: GPU smoke fires both mechanisms with counters comparable to
  CPU, no per-step cost blowup from the added sync (budget: the seam already
  moves ~2×67 MB per maintenance pass; rs adds ~7/45 more rows).

## After it works

1. Propose commit sequencing to Ryan (this stacks on the possibly-still-
   uncommitted §19 + §20.7 layers — coordinate; do not commit the other
   agent's files without his say).
2. Report back for launch: the 14 gated arms
   (`scr_p026s9_*`, `scr_p026sp_nt144_cap{030,018}`) can then run GPU with
   the §21 settings. Open launch questions Ryan has still not ruled:
   all-14 vs 3-arm de-risk, cold vs warm-start (s020v VTPs bracket
   ignition), whether `f_visc` stays off (core spreading feeds `dvisc`;
   enabling it changes the cost model a lot), and whether the +0.39%
   guard perturbation drives the floor below 0.25. He also still owes what
   he wanted clarified about the old question set — ask, don't re-ask the
   old questions verbatim.

## Notebook

Append-only, ASK Ryan before writing anything, ASK how verbose. Owed from
earlier arcs (see `split_launch_reset_prompt_20260911.md`): expguard
three-arm result should lead; plus Cd transient, rlxf derivation, Ladder C
forensics.
