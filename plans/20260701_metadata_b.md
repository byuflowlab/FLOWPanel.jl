# Improvement 1 is superseded by `plans/20260701_metadata_c.md`; Improvements 2-4 are absorbed there.

# Review & Improve: `plans/20260701_metadata.md` (replay/metadata serialization coverage)

## Context

The completed plan extended `<run>.metadata.toml` (de)serialization so more particle-wake
settings round-trip instead of silently reconstructing to a wrong default. I reviewed the
implementation (working-tree diff on `src/FLOWPanel_replay.jl`, `src/FLOWPanel_metadata.jl`,
`test/runtests_unit_replay.jl`) against the actual FLOWVPM/FLOWPanel source.

**Review verdict — the implemented work is correct.** Every tag string, constructor signature,
and struct field was verified against source:
- SFS singletons exist by those names; `SFS_default === SFS_none` (alias), all are immutable
  `<: SubFilterScale`, so the new `_vpm_optarg_manifest(::FLOWVPM.SubFilterScale)` dispatch and the
  `===` identity checks are sound (`FLOWVPM.jl/src/FLOWVPM.jl:221-278`, `FLOWVPM_subfilterscale.jl:89,136,168,245`).
- `ParticleStrengthExchange{R}(nu; recalculate_vols=true)` matches the serialize/deserialize branches (`FLOWVPM_viscous.jl:231`).
- `RelativeMinGamma(fraction)`, `SplitParticles(every, opts, verbose)`, and
  `CylindricalFieldProbeMonitor` fields all match (`FLOWPanel_wake.jl:689,772`, `FLOWPanel_simulate_monitors_fieldprobe.jl:23`).
- `julia --project -e 'include("test/runtests_unit_replay.jl")'` → **64/64 pass**.

The rest of this plan is the **substantial improvements** found during review. Decisions confirmed
with the user: (Q1) unrecognized/`unsupported` tags on read should **warn and keep loading** (non-breaking);
(Q2) prefer **reflection-based tagging** — record the actual binding name and resolve it blindly at
read time via `getfield(FLOWVPM, …)`, so a previously-written type that FLOWVPM later drops errors on
the *FLOWVPM side* (a clean `UndefVarError`) rather than needing a hand-maintained registry.

---

## Improvement 1 — Reflection-based singleton tagging (replaces the hardcoded SFS chains)

The recurring failure mode (a new FLOWVPM SFS gets added but someone forgets to update *both*
enum chains → silent wrong default) is structural. The just-added SFS if-chains extend the
enumeration but don't remove the fragility. Replace them with two generic helpers, then delete the
per-singleton branches.

Add to `src/FLOWPanel_replay.jl`:
- `_singleton_tag(mod, obj) -> Union{String,Nothing}`: scan `names(mod; all=true)`, guard with
  `isdefined(mod, n)` and skip gensym names (contain `#`), return `"FLOWVPM.<n>"` for the first
  `n` whose `getfield(mod, n) === obj` (deterministic: `names` is sorted; aliases like
  `SFS_default`/`SFS_none`/`noSFS` all resolve to the same object on read, so any match is safe).
  Return `nothing` for user-constructed (non-`const`) objects.
- `_resolve_singleton(tag) -> Any`: split `"FLOWVPM.<name>"`; only the `FLOWVPM` module is in scope,
  so `getfield(FLOWVPM, Symbol(name))`. A tag whose binding no longer exists throws `UndefVarError`
  here — this is the intended "error on the FLOWVPM side" (Q2).

Rewire:
- `_sfs_manifest(sfs)`: `t = _singleton_tag(FLOWVPM, sfs); t === nothing ? _metadata_unsupported_dict(typeof(sfs)) : Dict("type" => t)`.
  Keep the explicit `sfs === nothing` branch. **Delete** the six hardcoded `elseif sfs === FLOWVPM.SFS_x` branches.
- `_deserialize_sfs(meta)`: read `stype`; if `"nothing"`/`"unsupported"`/missing → warn + fallback
  (Improvement 2); else `return _resolve_singleton(stype)`. **Delete** the six hardcoded tag branches.
  (Legacy `"FLOWVPM.SFS_default"` resolves for free — the const still exists.)
- `_vpm_optarg_manifest(::FLOWVPM.SubFilterScale) = _sfs_manifest(val)` stays (already added, correct).

Net effect: SFS coverage is complete for all current *and future* named singletons with zero
per-type maintenance; only genuinely un-nameable (inline-constructed) SFS fall to `unsupported`.
Leave the relaxation (`_relaxation_manifest`) and viscous singleton-vs-data handling as-is for now,
but note in a comment that `_singleton_tag`/`_resolve_singleton` are the reusable mechanism.

## Improvement 2 — Loud read-side fallback (Q1: warn, keep loading)

The deserializers still silently substitute a *physically different* default when a tag is
`"unsupported"` or unknown — the same corruption class the plan set out to kill, on the read side.
For a restart/reproduction tool this means replay silently runs a different simulation.

In each of `_deserialize_sfs`, `_deserialize_viscous`, and `_deserialize_wake_shedding`, replace the
bare `else -> <default>` with a warn-then-default: e.g.
`@warn "replay: unrecognized SFS tag \"$stype\"; falling back to SFS_none (result may differ)"`.
Non-breaking: existing manifests with `unsupported` markers still load, but now emit a signal.
(`_deserialize_particle_policy`'s `else -> nothing` is fine as-is — a dropped trim policy is inert,
not a silent physical substitution — but a `@debug` there is optional.)

## Improvement 3 — Fix the dead-ternary kernel bug in `_deserialize_viscous` (`src/FLOWPanel_replay.jl:550`)

Pre-existing, sits three lines above the new PSE branch, same silent-wrong-default class:
```julia
kernel_fn = kernel == "FLOWVPM.zeta_fmm" ? FLOWVPM.zeta_fmm : FLOWVPM.zeta_fmm  # both branches identical
```
Any non-default CoreSpreading kernel is silently coerced to `zeta_fmm`. Resolve it generically with
the Improvement-1 helper: `kernel_fn = _resolve_singleton(String(get(meta, "kernel", "FLOWVPM.zeta_fmm")))`
(kernels are `const` functions in FLOWVPM, so they tag/resolve the same way). This also makes the
CoreSpreading `kernel` field a real round-trip.

## Improvement 4 — Restore dropped test coverage; add reflection round-trip checks

The diff *replaced* `CoreSpreading`+`SFS_Cd_twolevel_nobackscatter` in the main testset with
PSE+`SFS_Cs_nobackscatter`, so the CoreSpreading round-trip and the (previously sole) two-level SFS
round-trip are no longer asserted anywhere. In `test/runtests_unit_replay.jl`:
- Add a small loop asserting reflection round-trips each named singleton `===` its original:
  `SFS_none, SFS_Cs_nobackscatter, SFS_Cd_twolevel_nobackscatter, SFS_Cd_twolevel_nobackscatter_projection,
  SFS_Cd_twolevel_backscatter_signed, SFS_Cd_threelevel_nobackscatter` (reuse the existing
  `_write_metadata_toml`→`replay` helper already in the `SFS_Cd_threelevel` sub-block).
- Re-add one `CoreSpreading` round-trip asserting the reconstructed `kernel` (guards Improvement 3).
- Assert an inline/unknown SFS tag triggers the warn path and returns the documented default
  (`@test_logs (:warn,) …` or `@test_warn`), documenting the Q1 behavior.

## Minor notes (low priority, optional)

- `_deserialize_particle_policy`'s explicit `elseif ptype == "SplitParticles"; return nothing`
  duplicates the `else`; keep only if it carries a comment explaining the intentional partial
  coverage, else drop it.
- `CylindricalFieldProbeMonitor` metadata also cannot recover `azimuth_range` (same reason as
  axial/radial — baked into `positions_local`) and drops `backend`; record-only, acceptable — just
  extend the plan's noted caveat.

---

## Verification

- `julia --project -e 'include("test/runtests_unit_replay.jl")'` — expect all existing + new
  assertions to pass (baseline is 64/64 after the reviewed work).
- Smoke: in a REPL, build a wake with each SFS singleton and with a `CoreSpreading(...; kernel)`,
  call `FLOWPanel._wake_manifest_dict(wake, 1)`, confirm `pfield_optargs["SFS"]["type"]` is the
  reflected `"FLOWVPM.<name>"` (not `"unsupported"`) and `viscous["kernel"]` round-trips.
- Confirm the pre-existing, unrelated `test/runtests_unit_warmstart.jl` `basis=:quad` failure is
  unchanged (do not fix here).
