# Replay metadata: reflection-based FLOWVPM tagging (no lists) + migration function

Standalone plan — executable without prior conversation context. Supersedes Improvement 1 of
`plans/20260701_metadata_b.md`; carries its Improvements 2–4 forward inline.

## Context

FLOWPanel writes each run's settings to `<run>.metadata.toml` and can reconstruct a simulation
from it (`src/FLOWPanel_replay.jl`, "replay"). Particle-wake settings that reference FLOWVPM
objects (SFS models, viscous-scheme kernels) are serialized as string tags like
`"FLOWVPM.SFS_Cs_nobackscatter"`. Today both directions are **hand-maintained if-chains**:
`_sfs_manifest` (~`src/FLOWPanel_replay.jl:177`) enumerates six SFS singletons on write, and
`_deserialize_sfs` (~`src/FLOWPanel_replay.jl:569`) enumerates the same six on read. The failure
mode is structural: when FLOWVPM adds a singleton and the chains aren't updated, the manifest
records `unsupported` and replay **silently substitutes a physically different default** —
corrupted reproductions with no signal.

**Decision (user-confirmed):** replace the chains with pure reflection — record the actual FLOWVPM
binding name at write time and resolve it back via `getfield` at read time. **No stored
registries/lists anywhere in the runtime path.** Backward compatibility is handled by a one-shot
`migrate_metadata_toml` function, not runtime alias tables. Unrecognized/`unsupported` tags on
read must **warn and keep loading** (never silently substitute); a tag whose FLOWVPM binding no
longer exists must **error loudly**.

**Empirical findings (verified against the installed FLOWVPM at
`/Users/ryan/Dropbox/research/projects/FLOWVPM.jl`) — the reflection scheme is sound:**
- Scanning `names(FLOWVPM; all=true)` for bindings `isa FLOWVPM.SubFilterScale` finds exactly the
  6 SFS singletons: `SFS_Cs_nobackscatter` (`ConstantSFS`), `SFS_Cd_twolevel_nobackscatter`,
  `SFS_Cd_twolevel_nobackscatter_projection`, `SFS_Cd_twolevel_backscatter_signed`,
  `SFS_Cd_threelevel_nobackscatter` (`DynamicSFS`), and the alias group
  `SFS_default === SFS_none === noSFS` (`NoSFS`).
- `names` returns sorted symbols, so first-match tagging is deterministic (the no-SFS group tags
  as the alphabetically-first alias, `SFS_default`). Any alias resolves to the same object on
  read, so which alias gets written is immaterial.
- Immutable egal: a freshly-constructed `ConstantSFS` with fields identical to
  `SFS_Cs_nobackscatter` is `===` the const, so equivalent user-constructed SFS objects tag
  correctly rather than falling to `unsupported`. Only genuinely novel (differently-parameterized)
  objects are un-nameable.
- All SFS binding names are **unexported** FLOWVPM internals. Accepted trade-off of the no-list
  design: if FLOWVPM renames a binding, old files error loudly at read; the fix is a reactive
  one-line entry in the migration function's remap dict (empty today), never a maintained runtime
  list.
- `getfield` resolution must be type-validated, else a wrong/corrupted tag resolves to an
  arbitrary FLOWVPM function/type and fails confusingly downstream.
- `FLOWVPM.CoreSpreading` stores its kernel in the `zeta` field. The write side
  (`_viscous_manifest`, ~`src/FLOWPanel_replay.jl:161`) currently hardcodes
  `"kernel" => "FLOWVPM.zeta_fmm"` regardless of the actual kernel, and the read side
  (~`src/FLOWPanel_replay.jl:550`) has a dead ternary
  (`kernel == "FLOWVPM.zeta_fmm" ? FLOWVPM.zeta_fmm : FLOWVPM.zeta_fmm`) — the kernel never
  round-trips.

Files touched: `src/FLOWPanel_replay.jl`, `test/runtests_unit_replay.jl`,
`plans/20260701_metadata_b.md` (note only). Baseline: `test/runtests_unit_replay.jl` passes 64/64.

---

## Step 1 — Mark the old plan superseded

Add one line at the top of `plans/20260701_metadata_b.md`: its Improvement 1 is superseded by
`plans/20260701_metadata_c.md` (this plan); Improvements 2–4 are absorbed here.

## Step 2 — Reflection helpers in `src/FLOWPanel_replay.jl`

Add near the other manifest helpers:

- `_singleton_tag(mod::Module, obj)::Union{String,Nothing}` — iterate `names(mod; all=true)`;
  skip names containing `'#'` (gensyms); guard `isdefined(mod, n)`; return `"$(nameof(mod)).$n"`
  for the first `n` with `getfield(mod, n) === obj`; return `nothing` if no binding matches
  (caller falls back to `_metadata_unsupported_dict`). Deterministic because `names` is sorted.
- `_resolve_singleton(tag::AbstractString, expected::Type)` — require prefix `"FLOWVPM."`
  (else `ArgumentError` quoting the tag); `sym = Symbol(chopprefix(tag, "FLOWVPM."))`;
  if `!isdefined(FLOWVPM, sym)` throw `ArgumentError` quoting the tag, suggesting a FLOWVPM
  version mismatch, and pointing at `migrate_metadata_toml`; `x = getfield(FLOWVPM, sym)`;
  if `!(x isa expected)` throw `ArgumentError` quoting the tag, `typeof(x)`, and `expected`;
  return `x`.

## Step 3 — Rewire the SFS manifest/deserializer

- `_sfs_manifest(sfs)`: keep the `sfs === nothing` branch (writes `"type" => "nothing"`); replace
  the six hardcoded `elseif sfs === FLOWVPM.SFS_x` branches with:
  `t = _singleton_tag(FLOWVPM, sfs); t === nothing ? _metadata_unsupported_dict(typeof(sfs)) : Dict{String,Any}("type" => t)`.
- `_deserialize_sfs(meta)`: if the tag is `"nothing"`, `"unsupported"`, or missing → `@warn`
  ("replay: … falling back to SFS_none (result may differ)") and return `FLOWVPM.SFS_none`;
  otherwise `return _resolve_singleton(stype, FLOWVPM.SubFilterScale)`. Delete the six hardcoded
  branches. (Existing files tagged `"FLOWVPM.SFS_default"` etc. resolve for free — those bindings
  still exist.)
- The `_vpm_optarg_manifest(val::FLOWVPM.SubFilterScale) = _sfs_manifest(val)` dispatch
  (~`src/FLOWPanel_replay.jl:239`) stays. Leave `_relaxation_manifest` and the FLOWPanel-native
  tag branches (`OverlapPPS`, `RelativeMinGamma`, …) as-is — those are parameterized constructors,
  not singletons; add a comment noting `_singleton_tag`/`_resolve_singleton` are the reusable
  mechanism for any future singleton-tag need.

## Step 4 — Loud read-side fallbacks (kill the remaining silent substitutions)

In `_deserialize_viscous` (~line 544) and `_deserialize_wake_shedding` (~line 484), replace each
bare `else → <default>` with a warn-then-default, e.g.
`@warn "replay: unrecognized viscous tag \"$vtype\"; falling back to Inviscid (result may differ)"`.
(`_deserialize_sfs` is covered by Step 3.) Leave `_deserialize_particle_policy`'s `else → nothing`
silent — a dropped trim policy is inert, not a physical substitution; drop or comment its
redundant `elseif ptype == "SplitParticles"; return nothing` branch while there.

## Step 5 — Fix the CoreSpreading kernel round-trip via the new helpers

- Read (~`src/FLOWPanel_replay.jl:550`): replace the dead ternary with
  `kernel_fn = _resolve_singleton(String(get(meta, "kernel", "FLOWVPM.zeta_fmm")), Function)`.
- Write (`_viscous_manifest` CoreSpreading branch): replace the hardcoded
  `"kernel" => "FLOWVPM.zeta_fmm"` with the reflected `_singleton_tag(FLOWVPM, viscous.zeta)`
  (fall back to the literal `"FLOWVPM.zeta_fmm"` only if the tag is `nothing`).

## Step 6 — `migrate_metadata_toml` (backward compatibility)

New function in `src/FLOWPanel_replay.jl`:

```julia
migrate_metadata_toml(path::AbstractString; backup=true) -> path
```

- Parse the TOML. For each wake entry:
  - Move legacy top-level `"viscous"`/`"SFS"`/`"relaxation"` keys under `"pfield_optargs"`
    (the read path currently special-cases legacy placement at ~`src/FLOWPanel_replay.jl:457-459`;
    migration makes files self-consistent so that fallback can eventually be dropped — don't drop
    it in this change).
  - Canonicalize every `"FLOWVPM.*"` tag **list-free** by round-tripping through reflection:
    `_singleton_tag(FLOWVPM, _resolve_singleton(tag, Any))` (e.g. `"FLOWVPM.noSFS"` →
    `"FLOWVPM.SFS_default"`). Tags that fail to resolve: first consult `_LEGACY_TAG_REMAP`, a
    `Dict{String,String}` defined next to this function and **empty today** — the single
    designated place a mapping ever lives if FLOWVPM renames a binding; it is *not* consulted by
    the runtime read path. Still-unresolvable tags are collected and reported in one error listing
    them all.
- `backup=true` writes `<path>.bak` before overwriting. Idempotent on an already-current file.

## Step 7 — Tests in `test/runtests_unit_replay.jl`

The prior diff replaced the only `CoreSpreading` and two-level-SFS round-trip tests; restore
coverage reflection-natively (reuse the existing `_write_metadata_toml` → `replay` helper already
used in the `SFS_Cd_threelevel` sub-block):

- **Self-updating SFS round-trip loop**: for every binding in `names(FLOWVPM; all=true)` whose
  value `isa FLOWVPM.SubFilterScale`, assert
  `_deserialize_sfs(_sfs_manifest(getfield(FLOWVPM, n))) === getfield(FLOWVPM, n)`. Covers all
  current and future singletons with no test edits.
- `CoreSpreading` round-trip asserting the reconstructed viscous scheme's `zeta === FLOWVPM.zeta_fmm`
  (guards Step 5, both sides).
- `@test_logs (:warn,) …` for an `"unsupported"`/unknown SFS entry returning `SFS_none`
  (documents the warn-and-continue behavior).
- `_resolve_singleton` error paths: `"FLOWVPM.SFS_does_not_exist"` → `ArgumentError` mentioning
  the tag and migration; `"FLOWVPM.run_vpm!"` with `expected=FLOWVPM.SubFilterScale` →
  `ArgumentError` (type mismatch); `"NotVPM.SFS_none"` → `ArgumentError` (bad prefix).
- Migration test: write a legacy-style TOML fixture (top-level `"SFS"`/`"viscous"` keys, alias tag
  `"FLOWVPM.noSFS"`); run `migrate_metadata_toml`; assert keys moved under `"pfield_optargs"`,
  tags canonicalized, `.bak` created, a second run is a no-op, and the migrated file replays.

## Minor note (optional)

`CylindricalFieldProbeMonitor` metadata cannot recover `azimuth_range` (baked into
`positions_local`, like axial/radial) and drops `backend`; record-only is acceptable — extend the
existing caveat comment.

## Verification

- `julia --project -e 'include("test/runtests_unit_replay.jl")'` — baseline 64/64 plus all new
  assertions pass.
- REPL smoke: build a wake with each SFS singleton and a `CoreSpreading(nu, sgm0, kernel)`, call
  `FLOWPanel._wake_manifest_dict(wake, 1)`, confirm `pfield_optargs["SFS"]["type"]` is the
  reflected `"FLOWVPM.<name>"` (not `"unsupported"`) and `viscous["kernel"]` round-trips.
- Run `migrate_metadata_toml` on a scratchpad copy of a real run's metadata (e.g.
  `data/rotor_hover_pressure_comparison/rotor_hover_pressure_comparison.metadata.toml`) and
  confirm it still replays.
- The pre-existing, unrelated `test/runtests_unit_warmstart.jl` `basis=:quad` failure must remain
  unchanged (do not fix here).
