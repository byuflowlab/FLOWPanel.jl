# 018 expguard relaunch — campaign provenance (2026-09-08)

Guarded re-run of the three `euler_exp` NT36 arms with the `sigma_guard` gain-ratio
clamp active on the `euler_exp` path. Authorised by Ryan 2026-09-08
(`SIGMA_FLOOR_FRAC=0.25`, "floor/cap", relaunch all three arms).

## Pins (annotated tags, all named `campaign/p018-expguard-20260908`)

| repo | worktree (cluster) | commit | base | note |
|---|---|---|---|---|
| FLOWVPM.jl | `~/wt018/FLOWVPM-expguard` | `7468712` | `3315b22` (`unified-052`) | guard on scalar + broadcast `euler_exp` |
| FLOWPanel.jl | `~/wt018/FLOWPanel-expguard` | `d5dd772` | `67b7383` (`p018cs-fp-nt-20260908`) | forwards `sigma_guard` to `euler_exp` |
| FastMultipole | `~/wt018/FastMultipole-expguard` | `3da58a1a` | `3da58a1a` (`unified-052`) | unmodified |

Local (development) counterparts, from which the cluster port was derived:

| repo | branch | commit | tag |
|---|---|---|---|
| FLOWVPM.jl | `flowpanel` | `21eeaaa` | `campaign/p018-expguard-20260908` |
| FLOWPanel.jl | `fastmultipole` | `7dee1ab` | `campaign/p018-expguard-20260908` |

The cluster port was applied with `~/port_guard.py` (FLOWVPM) and an equivalent
anchored replacement (FLOWPanel, whose `src/FLOWPanel_wake.jl` at `67b7383`
differs from the `4e6b5b7` variant the script targets but carries the identical
`euler_exp` branch). The ported FLOWVPM file was diffed against local `21eeaaa`:
the only differences are the BRAINSTORM-026 resolution-split / `dsigma2`
accumulator lines that the cluster branch does not carry. In the local file the
guard bounds are evaluated from the `sig_before` local; in the cluster file from
`get_sigma(p)[]` read at the same point (before `get_sigma(p)[] *= rc^(-g)`),
i.e. the same value.

**Base choice:** `67b7383` is the exact FLOWPanel tree that ran the three
unguarded arms (jobs 13603853, 13605985, 13605986 — `sacct` WorkDir
`~/wt018/FLOWPanel-nt-gh200` and `~/wt018/FLOWPanel-nt-cs-gh200`; `67b7383` is a
strict superset of the control's `d2d3bb7`, adding only the `SFS_CONST_CS`
driver branch). Basing on `unified-052` (`4e6b5b7`) would have changed the
driver and `src/FLOWPanel_wake.jl` and confounded the comparison.

## Julia environment

`~/p018wtenv-expguard-gh200` — copy of `~/p018wtenv-nt-cs-gh200` with all three
dev-paths repointed at the worktrees above. Depot `~/fm052depot-gh200`,
julia `~/julia/julia-1.11.7/bin/julia` (aarch64).

## Arms

Case tag `p018_csarc_l3p0` for all three; knobs by environment.

Common: `WAKE_EXPINT=true`, `P018_SETTLE_REVS=22`, `SIGMA_FLOOR_FRAC=0.25`
(floor = 0.25 x 0.0047597 = 0.00119 m), `SIGMA_CEIL=0.030`,
`P018_REPO_OVERRIDE=~/wt018/FLOWPanel-expguard`,
`P018_PROJECT_OVERRIDE=~/p018wtenv-expguard-gh200`.

| run name | extra env | unguarded predecessor | predecessor fate |
|---|---|---|---|
| `p018_csarc_l3p0_3r_exp_nt_g25` | (dynamic Cd) | 13603853 | COMPLETED, healthy to rev 30 — control for guard perturbation |
| `p018_csarc_l3p0_3r_cs0p002_exp_nt_g25` | `SFS_CONST_CS=0.002` | 13605985 | died step 898 (rev 24.9), slow sigma-collapse precursor |
| `p018_csarc_l3p0_3r_cs0p34_exp_nt_g25` | `SFS_CONST_CS=0.34` | 13605986 | died step 361 (rev 10.0), no sigma precursor |

New run names (`_g25` suffix) so the existing data dirs are not clobbered.

## Smoke test

Job 13610775 (gh200, 40 min): guard unit checks on the ported CPU scalar path
(empty guard bit-exact, floor/ceil/dtz_cap bounds binding exactly,
`|Gamma| sigma^2` conserved, unknown-key throw) plus a full precompile of the
FLOWPanel/CUDA stack in the new environment.

## Submitted

| job | run name | state at submit |
|---|---|---|
| 13610777 | `p018_csarc_l3p0_3r_exp_nt_g25` | RUNNING mgh-1-1; banner `guard=on`, floor 0.00119 m, ceil 0.03 m, `DynamicSFS(rlxf=0.005)` |
| 13610778 | `p018_csarc_l3p0_3r_cs0p002_exp_nt_g25` | RUNNING mgh-1-2; banner `guard=on`, `ConstantSFS(Cs=0.002)` |
| 13610779 | `p018_csarc_l3p0_3r_cs0p34_exp_nt_g25` | PENDING (Resources) |

All three: 1079 steps, ~3-5 s/step, 24 h wall.
