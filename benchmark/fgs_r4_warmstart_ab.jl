#=##############################################################################
BRAINSTORM 021 — warm-start R4 head-to-head, per-arm entry
(fgs_warmstart_r4_reset_prompt_20260924.md, Job 1).

FGS (dagteam+backoff production default) vs krylov_ilu_nfcache
(persistent_plan), warm-started, wake-on, R4, NT=36, scored over BOTH windows
INCLUDING their transients (Ryan 2026-09-24): Window A = steps 1-36 (startup
revolution, from the very first step), Window B = steps 109-144 (fourth
revolution). No settling exclusion anywhere (SKIP_STEPS=0).

Checkpoint + restart layout (Ryan 2026-09-24, so revs 2-3 are simulated once
per solver FAMILY instead of once per arm — see the WSR4_LEG block below):
the two cold arms each march revs 1-3 with VTK on (leg=ckpt: their own
Window A + the family restart source), warm arms march rev 1 from scratch
(leg=winA), and every arm's Window B is a rev-4 leg restarted from its
family's checkpoint (leg=winB; symmetric treatment cold and warm).

One process per (arm, leg). ARM selects the configuration; everything else
is shared. Cold = zero-initial-guess solves in the SAME process (Ryan
2026-09-23), which is exactly WARMSTART=cold here.

All simulation/measurement mechanics live in rotor_hover_solver_unsteady.jl
(the ruling-12 Phase-3 vehicle): per-step t_solve/niter_first/nsolves/solved,
per-step certified BC residual (the fixed-accuracy convergence contract's
per-step spot check), per-step CT identity, per-step t_project (warm-start
guess-construction cost), per-step solved-strength snapshots (cross-solver
solution-agreement deliverable), and the setup-cost split (ILU factorization
lives for the whole run: built once, exactly rigid-invariant for this
Dirichlet body — design check 2026-09-24 — so there are no per-step
factorization events; the plan+nfcache build is pulled out of step 1 by the
priming solve and reported as t_setup_prime).

Pre-sim setup costs (nfcache/plan build, ILU factorization, solver
construction) are EXCLUDED from all per-step cost comparisons and reported
once per arm in the setup columns (Ryan 2026-09-24). No amortization curves.

Required env:
  ARM           one of the ARM_TABLE keys below
  FGS_TOL_ABS   explicit FGS stopping tolerance (FGS arms; re-staircased for
                this fixture/environment — never carried from the cold table)
Optional env (defaults set here; the unsteady driver documents the rest):
  RUNG (R4), N_STEPS (144), NT (36), KNOBS_* / KNOBS_BUDGET (apply knobs),
  FGS_KNOBS_FILE (benchmark/retained_r4_champion.toml), FGS_PRECISION (f64 —
  f32full is certified COLD R4 only, do not enable without re-certification),
  OUTDIR_OVERRIDE (campaign data root), SNAPSHOT_STRENGTHS (1)

Writes STATUS_<ARM> (running/ok/FAILED) and, on success, COMPLETED_<ARM> in
the output directory. Judge arms by these outputs, never sacct.

Local smoke (R1, all arms, few steps, <=4 threads):
  bash benchmark/run_r4_fgs_warmstart_smoke.sh
=###############################################################################

const ARM_TABLE = Dict(
    # arm name          => (CONFIG,               WARMSTART, WARMSTART_ORDER)
    "fgs_cold"          => ("fgs",                "cold",    1),
    "fgs_prev"          => ("fgs",                "prev",    1),
    "fgs_proj1"         => ("fgs",                "extrap",  1),
    "fgs_proj2"         => ("fgs",                "extrap",  2),
    "ilu_nfcache_cold"  => ("krylov_ilu_nfcache", "cold",    1),
    "ilu_nfcache_prev"  => ("krylov_ilu_nfcache", "prev",    1),
    # optional arm (Ryan may strike): extrapolated x0 reusing the SAME shared
    # extrapolation coefficients as FGS's project_solution!, for comparability
    "ilu_nfcache_proj1" => ("krylov_ilu_nfcache", "extrap",  1),
    # added 2026-09-26 (Ryan): quadratic ILU arm for a proj2-vs-proj2
    # head-to-head with fgs_proj2 (closes the best-vs-best asymmetry)
    "ilu_nfcache_proj2" => ("krylov_ilu_nfcache", "extrap",  2),
)

arm = get(ENV, "ARM", "")
haskey(ARM_TABLE, arm) ||
    error("ARM must be one of $(sort(collect(keys(ARM_TABLE)))); got $(repr(arm))")
config_arm, warmstart_arm, ws_order_arm = ARM_TABLE[arm]

# --- legs (Ryan 2026-09-24: checkpoint + restart so revs 2-3 are not
# re-simulated per arm; TWO checkpoints, one per solver family) -------------
#   ckpt — cold arms only: march revs 1-3 (3*NT steps) from scratch with
#          SAVE_VTK=true. Doubles as the cold arm's Window A data (its first
#          NT steps) and as the family's shared restart source.
#   winA — warm arms: march rev 1 (NT steps) from scratch, VTK off.
#   winB — ALL arms (cold included, for symmetric treatment): restart from
#          the FAMILY checkpoint at step 3*NT and march rev 4. REPORTING
#          NOTE (Ryan 2026-09-24): solver warm-start histories are not
#          serialized in the checkpoint, so the first (order+1) steps of a
#          restarted warm leg are effectively cold — the history-fill
#          transient lands inside Window B's transient-included stats and
#          must be flagged in the harvest, never silently excluded.
#   full — the original single-march 4-revolution behavior (kept for
#          flexibility/smokes; not used by the staged campaign).
leg = get(ENV, "WSR4_LEG", "full")
leg in ("full", "ckpt", "winA", "winB") ||
    error("WSR4_LEG must be full, ckpt, winA, or winB; got $(repr(leg))")
family = config_arm == "fgs" ? "fgs" : "ilu"
leg == "ckpt" && warmstart_arm != "cold" &&
    error("WSR4_LEG=ckpt is generated by the family's COLD arm only (got $arm)")

_setdefault!(k, v) = haskey(ENV, k) || (ENV[k] = v)
# p033 thread-scaling campaign: 018-ported wake-physics defaults (inert unless
# WAKE_ENV_018=1; see the file header for provenance)
include(joinpath(@__DIR__, "p033_wake_env_018.jl"))
ENV["CONFIG"] = config_arm
ENV["WARMSTART"] = warmstart_arm
ENV["WARMSTART_ORDER"] = string(ws_order_arm)
_setdefault!("RUNG", "R4")
_setdefault!("NT", "36")
nt_arm = parse(Int, ENV["NT"])
ckpt_name = "fgs_wsr4_$(ENV["RUNG"])_ckpt_$family"
if leg == "ckpt"
    _setdefault!("N_STEPS", string(3nt_arm))   # revs 1-3
    ENV["SAVE_VTK"] = "true"                   # the restart source
    _setdefault!("RUN_NAME", ckpt_name)
    ENV["RESTART_STEP"] = "-1"
elseif leg == "winA"
    _setdefault!("N_STEPS", string(nt_arm))    # rev 1, from the very first step
    _setdefault!("SAVE_VTK", "false")
    _setdefault!("RUN_NAME", "fgs_wsr4_$(ENV["RUNG"])_$(arm)_winA")
    ENV["RESTART_STEP"] = "-1"
elseif leg == "winB"
    _setdefault!("N_STEPS", string(4nt_arm))   # global schedule; marches rev 4
    _setdefault!("SAVE_VTK", "false")
    _setdefault!("RUN_NAME", "fgs_wsr4_$(ENV["RUNG"])_$(arm)_winB")
    _setdefault!("RESTART_STEP", string(3nt_arm))
    _setdefault!("RESTART_NAME", ckpt_name)
    _setdefault!("RESTART_PATH", joinpath("data", ckpt_name))
else
    # 4 full revolutions in one march; windows are cut at harvest
    _setdefault!("N_STEPS", string(4nt_arm))
    _setdefault!("SAVE_VTK", "false")
    _setdefault!("RUN_NAME", "fgs_wsr4_$(ENV["RUNG"])_$arm")
    ENV["RESTART_STEP"] = get(ENV, "RESTART_STEP", "-1")
end
# no settling exclusion anywhere (Ryan 2026-09-24); the summary printed by the
# unsteady driver then covers all steps, and windows are cut at harvest
ENV["SKIP_STEPS"] = "0"
# The fixture example defaults the particle .vtp series to Float32
# (visualization-focused). Campaign checkpoints are RESTART SOURCES: force the
# package-level f64 default so winB continuations are replay-faithful (the
# warm-start loader warns on a Float32 series). Explicit env still wins.
_setdefault!("FLOWPANEL_PARTICLE_PRECISION", "f64")
_setdefault!("PHASE", "phase3wsr4")          # phase3* => CT monitor on
_setdefault!("SNAPSHOT_STRENGTHS", "1")
if config_arm == "fgs" && !haskey(ENV, "FGS_P")
    _setdefault!("FGS_KNOBS_FILE",
                 joinpath(@__DIR__, "retained_r4_champion.toml"))
    # default f64 BOTH sides: f32full is certified cold-R4/zen3 only
    _setdefault!("FGS_PRECISION", "f64")
    haskey(ENV, "FGS_TOL_ABS") || error("FGS arms need an explicit FGS_TOL_ABS")
end

armleg = "$(arm)_$(leg)"
outdir_status = get(ENV, "OUTDIR_OVERRIDE", joinpath(@__DIR__, "results",
    ENV["PHASE"]))
mkpath(outdir_status)
_status(s) = write(joinpath(outdir_status, "STATUS_$armleg"), s * "\n")
_status("running")

try
    include(joinpath(@__DIR__, "rotor_hover_solver_unsteady.jl"))
    # A leg whose steps hit itmax is NOT ok even though the march completed:
    # judge convergence, not just survival (2026-09-24 smoke: the ilu winB leg
    # read ok while every restarted step was unconverged at CF ~1e7).
    n_unconverged = count(!, timed_formulation.solved)
    n_unconverged == 0 || error("$(n_unconverged) of $(length(timed_formulation.solved)) " *
        "steps unconverged (solved=false) — failing leg $armleg")
    _status("ok")
    write(joinpath(outdir_status, "COMPLETED_$armleg"),
          "completed $(time_string())\n")
catch err
    _status("FAILED")
    rethrow()
end
