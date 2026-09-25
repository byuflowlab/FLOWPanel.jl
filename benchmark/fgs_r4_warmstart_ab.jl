#=##############################################################################
BRAINSTORM 021 — warm-start R4 head-to-head, per-arm entry
(fgs_warmstart_r4_reset_prompt_20260924.md, Job 1).

FGS (dagteam+backoff production default) vs krylov_ilu_nfcache
(persistent_plan), warm-started, wake-on, R4, 4 revolutions (144 steps at
NT=36), scored over BOTH windows INCLUDING their transients (Ryan
2026-09-24): Window A = steps 1-36 (startup revolution, from the very first
step), Window B = steps 109-144 (fourth revolution). No settling exclusion
anywhere (SKIP_STEPS=0). Windows are harvest-side cuts of the per-step CSV;
this driver just runs the full schedule.

One process per arm (arms are separate Slurm tasks). ARM selects the
configuration; everything else is shared. Cold = zero-initial-guess solves in
the SAME process (Ryan 2026-09-23), which is exactly WARMSTART=cold here.

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
)

arm = get(ENV, "ARM", "")
haskey(ARM_TABLE, arm) ||
    error("ARM must be one of $(sort(collect(keys(ARM_TABLE)))); got $(repr(arm))")
config_arm, warmstart_arm, ws_order_arm = ARM_TABLE[arm]

_setdefault!(k, v) = haskey(ENV, k) || (ENV[k] = v)
ENV["CONFIG"] = config_arm
ENV["WARMSTART"] = warmstart_arm
ENV["WARMSTART_ORDER"] = string(ws_order_arm)
_setdefault!("RUNG", "R4")
_setdefault!("NT", "36")
# 4 full revolutions; Window A = rev 1 (steps 1-36), Window B = rev 4
# (steps 109-144), both harvested WITH their transients
_setdefault!("N_STEPS", string(4 * parse(Int, ENV["NT"])))
# no settling exclusion anywhere (Ryan 2026-09-24); the summary printed by the
# unsteady driver then covers all steps, and windows are cut at harvest
ENV["SKIP_STEPS"] = "0"
_setdefault!("PHASE", "phase3wsr4")          # phase3* => CT monitor on
_setdefault!("SNAPSHOT_STRENGTHS", "1")
_setdefault!("RUN_NAME", "fgs_wsr4_$(ENV["RUNG"])_$arm")
if config_arm == "fgs" && !haskey(ENV, "FGS_P")
    _setdefault!("FGS_KNOBS_FILE",
                 joinpath(@__DIR__, "retained_r4_champion.toml"))
    # default f64 BOTH sides: f32full is certified cold-R4/zen3 only
    _setdefault!("FGS_PRECISION", "f64")
    haskey(ENV, "FGS_TOL_ABS") || error("FGS arms need an explicit FGS_TOL_ABS")
end
ENV["RESTART_STEP"] = get(ENV, "RESTART_STEP", "-1")  # full schedule, no restart

outdir_status = get(ENV, "OUTDIR_OVERRIDE", joinpath(@__DIR__, "results",
    ENV["PHASE"]))
mkpath(outdir_status)
_status(s) = write(joinpath(outdir_status, "STATUS_$arm"), s * "\n")
_status("running")

try
    include(joinpath(@__DIR__, "rotor_hover_solver_unsteady.jl"))
    _status("ok")
    write(joinpath(outdir_status, "COMPLETED_$arm"),
          "completed $(time_string())\n")
catch err
    _status("FAILED")
    rethrow()
end
