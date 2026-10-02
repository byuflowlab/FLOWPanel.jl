#=##############################################################################
BRAINSTORM 021 — unsteady campaign driver (Phase 2 full-timestep share,
ruling 7; Phase 3 warmstart runs reuse this file, ruling 12).

Construction is NOT duplicated: this driver includes
examples/rotor_hover_pressure_comparison.jl (RHPC — the campaign's reference
unsteady driver, Ryan 2026-08-18) with RHPC_SETUP_ONLY=1, inheriting the
corrected DJI conventions (scale by radial dim 2, rotation about x, staged-
startup freestream pulse, RigidWakeBody shedding-from-constructed-cells,
stock relaxation) and the production wake. examples/rotor_hover.jl and its
stale Phantom-era lines are not touched or used.

What this driver adds:
  - CONFIG env selects the roster solver (backslash | krylov_gmres |
    krylov_jacobi | krylov_ilu | fgs | fgmres_fgs) built at the Phase-1
    frozen per-rung knobs (phase1_knobs.jl; FGS family needs the rung's
    fgstune/fgsprecond CSVs);
  - per-step instrumentation with NO src changes:
      t_step_total  — wall time between consecutive maneuver! calls
                      (maneuver! runs once at the start of every step);
      t_solve       — a TimedFormulation wrapper around the production
                      VelocityThroughSources formulation times
                      solve_formulation! (the body solve inside
                      _steady_aerodynamics!);
      niter         — iterations of the LAST inner solve of the step
                      (KrylovSolver.niter and, since Phase 3, FGSSolver.niter;
                      -1 for Backslash, a direct solve with no iteration
                      count). Kept for continuity with Phase 1/2; it is a
                      diagnostic, NOT the warmstart comparison metric;
      niter_first   — iterations of the FIRST inner solve of the step. This is
                      the solve that consumes the step-to-step initial guess,
                      so it is the Phase 3 warmstart headline metric;
      nsolves       — per-body solves the step took (1 for a one-body run;
                      -1 = unavailable for direct/generic solvers). Anything
                      other than 1 here means the timestep re-solved and the
                      warmstart comparison is invalid;
      solved        — AND of the step's inner-solver convergence flags;
  - per-step CSV rows (unsteady.csv) + banner.

Phase 3 (warmstart, ruling 12) adds, all defaulting to the Phase-2 behavior
so existing cold rows stay reproducible:
  - WARMSTART=cold|prev|extrap selects the initial guess for the iterative
    solvers (backslash is the null control and accepts cold only);
  - WARMSTART_ORDER sets the extrapolation order, using the SAME shared
    coefficients (pnl._extrapolation_coefficients) on both the Krylov and FGS
    sides so the two are comparable by construction;
  - RESTART_STEP/RESTART_PATH/RESTART_NAME resume from a wake-developed
    checkpoint via simulate_warmstart! instead of re-marching the 13-rev
    RHPC startup schedule in every arm;
  - PHASE routes output (phase2 default; the Phase-3 launcher sets phase3).

Accuracy guard (spec: "guard against silent accuracy drift"): the per-step
`solved` flag catches a warm guess that lets a step stop short, and the
end-of-run particle count + strength checksum let cold and warm arms of the
same config be compared for trajectory identity — both converge to the same
tolerance, so any divergence means the guess changed the answer. A certified
BC error is NOT recorded here: the campaign's frozen BC metric
(decision_rules.md) is defined against a static single-step RHS, whereas the
unsteady RHS moves every step, so reusing that column would misreport it.

Monitors default OFF (RUN_MONITORS=false: the metric is solver cost);
opt back in with RUN_MONITORS=true [BERNOULLI_ONLY=true] for CT
cross-checks. VTK off by default.

Run (local smoke; published rows are HPC-only per ruling 5):
  RUNG=R1 CONFIG=fgs N_STEPS=8 EXPECT_JULIA_THREADS=4 THREADING_MODE=multi \
    julia --project -t 4 benchmark/rotor_hover_solver_unsteady.jl
N_STEPS truncates the RHPC schedule (smoke); leave unset for the full
staged-startup + settle schedule (Phase 3 wake-developed regime).
=###############################################################################

include(joinpath(@__DIR__, "common.jl"))
using TOML

banner = assert_and_banner()

# --- ladder + env plumbing (keep LADDER in sync with phase1_case.jl / ledger.md)
const LADDER = Dict(
    "R1" => ("dji9443_20260813_23_73_capped_captess4.msh",    8016),
    "R2" => ("dji9443_20260813_33_105_capped_captess4.msh",  15760),
    "R3" => ("dji9443_20260813_45_145_capped_captess4.msh",  28752),
    "R4" => ("dji9443_20260813_65_209_capped_captess4.msh",  58192),
    "R5" => ("dji9443_20260813_89_289_capped_captess4.msh", 108240),
    "R6" => ("dji9443_20260813_125_409_capped_captess4.msh", 212108),
    "R7" => ("dji9443_20260813_177_577_capped_captess4.msh", 419276),
)
rung = get(ENV, "RUNG", "")
haskey(LADDER, rung) || error("RUNG must be one of $(sort(collect(keys(LADDER)))); got \"$rung\"")
msh_ladder, n_expected = LADDER[rung]
config = lowercase(get(ENV, "CONFIG", "backslash"))

# --- Phase 3 warmstart knobs (defaults reproduce the Phase-2 cold behavior) --
warmstart = lowercase(get(ENV, "WARMSTART", "cold"))
warmstart in ("cold", "prev", "extrap") ||
    error("WARMSTART must be cold, prev, or extrap; got $(repr(warmstart))")
warmstart_order = parse(Int, get(ENV, "WARMSTART_ORDER", "1"))
warmstart == "extrap" && warmstart_order < 1 &&
    error("WARMSTART=extrap needs WARMSTART_ORDER >= 1; got $warmstart_order")
# effective polynomial order handed to the solvers: 0 == plain previous-solution
ws_order = warmstart == "extrap" ? warmstart_order : 0
config == "backslash" && warmstart != "cold" &&
    error("CONFIG=backslash is the warmstart null control (a direct solve is " *
          "guess-independent); it accepts WARMSTART=cold only, got $(repr(warmstart))")

# --- Per-step BC check (Ryan 2026-08-25: measure the BOUNDARY CONDITION,
# not the solution) -----------------------------------------------------------
# The validity claim for a warm start is "this arm still satisfies the BCs to
# the tolerance it promised", and that is directly measurable — no need to
# compare solutions between arms and then argue about how a residual tolerance
# maps to a solution difference.
#
# `bc_error!` (benchmark/common.jl, 021 ruling 3) is exactly this instrument and
# needs only ONE influence pass: it reloads σ from the BC and the solved
# doublets together, so the control-point potential IS the BC residual
# φ_σ + G_μ x = G x − b. No reference solution, no dense G, at any rung.
#
# ACCURACY. Ryan asked for the measurement to be ~10x sharper than the thing
# measured. Rather than guess an expansion order, `bc_error!` REQUESTS an
# absolute FMM error target (`PowerAbsolutePotential`) and reports whether every
# M2L pair certified it. Requesting `BCERR_SAFETY` (default 0.1) x the arm's own
# stopping tolerance makes the instrument 10x sharper BY CONSTRUCTION, and
# `bcerr_certified=false` says loudly when it failed to be.
#
# The pass sits OUTSIDE the `t_solve` timer — it is measurement, not solve cost.
bcerr_safety = parse(Float64, get(ENV, "BCERR_SAFETY", "0.1"))
(isfinite(bcerr_safety) && bcerr_safety > 0) ||
    error("BCERR_SAFETY must be positive and finite; got $bcerr_safety")
bcerr_max_p = parse(Int, get(ENV, "BCERR_MAX_P", "20"))
# COST. Certifying an error target 10x below tol_abs makes the FMM climb its
# expansion order; the pass measured ~62% of t_solve at R1 (t_bcerr column). It
# is outside the t_solve timer, so the headline cost metric is unaffected, and
# Ryan ruled 2026-08-25 that the extra WALL time is acceptable (run twice, or
# recover it through replay if it ever matters). BCERR_EVERY therefore defaults
# to 1 — measure EVERY step — and exists only as an escape hatch; skipped steps
# write empty cells rather than a fabricated value.
bcerr_every = parse(Int, get(ENV, "BCERR_EVERY", "1"))
bcerr_every >= 1 || error("BCERR_EVERY must be >= 1; got $bcerr_every")

restart_step = parse(Int, get(ENV, "RESTART_STEP", "-1"))
restart_name = get(ENV, "RESTART_NAME", "")
restart_path = get(ENV, "RESTART_PATH", "")
# steps to drop from wake-developed averages: the warmstart history needs
# filling after a restart (Krylov order p and FGS order p both need p+1 solves)
skip_steps = parse(Int, get(ENV, "SKIP_STEPS", "3"))

_rung_sub = get(ENV, "PER_RUNG_DIR", "0") == "1" ? rung : ""
knobs_dir = joinpath(@__DIR__, "results", "phase1",
                     get(ENV, "KNOBS_MODE", banner.threading_mode), _rung_sub)
outdir = if !isempty(get(ENV, "OUTDIR_OVERRIDE", ""))
    # campaign runs write outside the (worktree) checkout — consolidated data
    # root discipline, see agent_policies/HPC.md
    isabspath(ENV["OUTDIR_OVERRIDE"]) ||
        error("OUTDIR_OVERRIDE must be an absolute path")
    ENV["OUTDIR_OVERRIDE"]
else
    joinpath(@__DIR__, "results", get(ENV, "PHASE", "phase2"),
             banner.threading_mode, _rung_sub)
end
mkpath(outdir)

target_rel = 1e-6
include(joinpath(@__DIR__, "phase1_knobs.jl"))

# tuned shared Krylov apply knobs: csv-first from knobs_dir/tune.csv, hardcoded
# R1–R3 local fallback (keep in sync with phase1_case.jl TUNED)
function tuned_knobs(rung)
    tuned = Dict("R1" => (17, 0.5, 21), "R2" => (17, 0.5, 12),
                 "R3" => (18, 0.5, 18))
    parsed = read_rows(joinpath(knobs_dir, "tune.csv"))
    if parsed !== nothing
        cols, rows = parsed
        for c in rows
            length(c) >= length(cols) || continue
            c[cols["rung"]] == rung || continue
            c[cols["variant"]] == "tuned" || continue
            tuned[rung] = (parse(Int, c[cols["expansion_order"]]),
                           parse(Float64, c[cols["multipole_acceptance"]]),
                           parse(Int, c[cols["leaf_size"]]))   # latest wins
        end
    end
    haskey(tuned, rung) || error("No tuned apply knobs for $rung — run the tune stage")
    return tuned[rung]
end

# Phase 2 owns tuning; Phase 3 owns amortization (decision_rules ruling 2), so
# under PHASE=phase3* the apply knobs come from Phase 2's cost-tuned LineGauss
# winners in results/phase2/<mode>/<rung>/tune_phase2.csv — NEVER from the
# Phase-1 tune.csv path above, whose Gaussian-era hardcoded fallback would
# otherwise engage silently (verified 2026-09-07: the campaign worktree has no
# Phase-1 knob CSVs at all, only bcache_*.bin). Phase 2 winners are keyed by
# memory budget, and the budgets are machine classes whose cached=true rows
# were tuned WITH the nearfield cache this driver does not build — so the
# budget is an explicit submit-line choice (KNOBS_BUDGET, in GiB, matching
# mem_budget_gib), not a default. The uncached-consistent row is budget 0
# (matrix-free endpoint), which as of 2026-09-07 exists for R3/R5 only.
function phase2_knobs(rung)
    budget = get(ENV, "KNOBS_BUDGET", "")
    isempty(budget) && error("PHASE=phase3 requires KNOBS_BUDGET (GiB, a " *
        "mem_budget_gib value from Phase 2's tune_phase2.csv): the winner " *
        "knobs differ per budget and cached-row winners assume a nearfield " *
        "cache this driver does not build, so the choice must be explicit")
    path = joinpath(@__DIR__, "results", "phase2",
                    get(ENV, "KNOBS_MODE", banner.threading_mode), rung,
                    "tune_phase2.csv")
    # a 0-byte file is a still-running tune job that has created but not yet
    # written its CSV — same meaning as missing, and read_rows would crash on it
    parsed = isfile(path) && filesize(path) > 0 ? read_rows(path) : nothing
    parsed === nothing && error("phase3 knobs: $path is missing or empty — " *
        "Phase 2 Step A has not landed for $rung")
    cols, rows = parsed
    sel = nothing
    avail = String[]
    for c in rows
        length(c) >= length(cols) || continue
        c[cols["rung"]] == rung || continue
        # below-floor placeholder rows carry bc_certified=false by construction
        c[cols["bc_certified"]] == "true" || continue
        push!(avail, c[cols["mem_budget_gib"]])
        c[cols["mem_budget_gib"]] == budget || continue
        sel = c   # latest wins (resume-replay rows carry identical knobs)
    end
    sel === nothing && error("phase3 knobs: no bc_certified row for $rung " *
        "at KNOBS_BUDGET=$budget GiB in $path (budgets present: " *
        "$(join(unique(avail), ", ")))")
    return (parse(Int, sel[cols["expansion_order"]]),
            parse(Float64, sel[cols["multipole_acceptance"]]),
            parse(Int, sel[cols["leaf_size"]]))
end

# --- RHPC setup-only include (env defaults only where the user set nothing) --
setdefault!(k, v) = haskey(ENV, k) || (ENV[k] = v)
ENV["RHPC_SETUP_ONLY"] = "1"
ENV["RHPC_MESH_FILE"] = msh_ladder
# The warmstart suffix keeps the three guess types of one Phase-3 job from
# sharing an output name (and, with SAVE_VTK=true, from appending into each
# other's VTK/TOML). It also keeps an arm from ever writing into the shared
# checkpoint directory it restarts FROM.
setdefault!("RUN_NAME", "rotor_hover_solver_unsteady_$(rung)_$(config)" *
            (warmstart == "cold" ? "" : "_$(warmstart)$(ws_order)"))
# Monitors stay OFF by default (the metric is solver cost) EXCEPT under
# PHASE=phase3, where the Bernoulli force monitor supplies the per-step CT that
# is now the arms' identity metric. Bernoulli-only keeps the added cost to the
# cheap pressure path; it lands outside `t_solve` either way. Phase 2 rows are
# unaffected, so their reproducibility story is untouched.
if startswith(get(ENV, "PHASE", "phase2"), "phase3")
    setdefault!("RUN_MONITORS", "true")
    setdefault!("BERNOULLI_ONLY", "true")
end
setdefault!("RUN_MONITORS", "false")
setdefault!("SAVE_VTK", "false")
setdefault!("NT", "36")

include(joinpath(pnl.examples_path, "rotor_hover_pressure_comparison.jl"))

rotor.ncells == n_expected ||
    error("Rung $rung: expected $n_expected panels, got $(rotor.ncells)")

# truncate the schedule for smokes; unset = full RHPC staged-startup schedule
n_steps_env = parse(Int, get(ENV, "N_STEPS", "0"))
if n_steps_env > 0
    t_range = range(0.0, step=dt, length=n_steps_env + 1)
end

# --- roster solver apply knobs: Phase-2 winners under phase3, Phase-1 frozen
# knobs otherwise (see phase2_knobs above) ------------------------------------
p_t, mac_t, leaf_t = if haskey(ENV, "KNOBS_P")
    # explicit apply knobs on the submit line (021 warm-start campaign / local
    # smokes without tune CSVs); all three required, recorded in the CSV
    (parse(Int, ENV["KNOBS_P"]), parse(Float64, ENV["KNOBS_MAC"]),
     parse(Int, ENV["KNOBS_LEAF"]))
elseif startswith(get(ENV, "PHASE", "phase2"), "phase3")
    phase2_knobs(rung)
else
    tuned_knobs(rung)
end
backend_apply = pnl.FastMultipoleBackend(; expansion_order=p_t,
    multipole_acceptance=mac_t, leaf_size=leaf_t)
krylov_rtol = target_rel
# Fixed-accuracy stopping for the warm-start campaign (PHASE=phase3*):
# Krylov.jl's rtol is relative to the INITIAL residual, so a warm x0 shrinks
# the target with it — same iterations, deeper solve, no time-to-target gain
# by construction. rtol_rhs anchors the target to ||b|| instead (identical
# for cold solves, where r0 = b); Phase-2 rows keep the historical semantics.
krylov_rtol_rhs = startswith(get(ENV, "PHASE", "phase2"), "phase3") ?
    krylov_rtol : 0.0
krylov_kw = (; itmax=500, atol=1e-14, rtol=krylov_rtol,
             rtol_rhs=krylov_rtol_rhs, memory=50,
             backend=backend_apply,
             # Phase 3: warmstart=false is the cold baseline; order 0 reuses the
             # previous solution, order>=1 extrapolates from the rolling history.
             warmstart=(warmstart != "cold"), warmstart_order=ws_order)

function make_config_solver(config)
    config == "backslash" && return pnl.Backslash(rotor), "dense LU"
    config == "krylov_gmres" && return pnl.KrylovSolver(rotor; method=:gmres,
        krylov_kw...), "rtol=$krylov_rtol;atol=1e-14;memory=50"
    config == "krylov_jacobi" && return pnl.KrylovSolver(rotor; method=:gmres,
        krylov_kw...,
        preconditioner=pnl.FastMultipole.JacobiPreconditioner((rotor,);
            cell_size=R/4)), "jacobi cell_size=R/4"
    config == "krylov_ilu" && return pnl.KrylovSolver(rotor; method=:gmres,
        krylov_kw...,
        preconditioner=pnl.ILUPreconditioner(rotor; leaf_size=10,
            multipole_acceptance=1.0,
            max_pattern_entries=8192 * rotor.ncells)), "ilu leaf10/mac1.0"
    if config == "krylov_ilu_nfcache"
        # 021 warm-start R4 campaign (fgs_warmstart_r4_reset_prompt_20260924):
        # same ILU as krylov_ilu, plus the persistent plan + dense near-field
        # cache. Setup costs are campaign deliverables reported SEPARATELY
        # from per-step cost (Ryan 2026-09-24), so each piece is timed here
        # and the plan/cache build is pulled out of step 1 by a priming solve
        # (prime_persistent_plan! below). The ILU factorization is built ONCE
        # here and never refactorized: transform_solver_geometry! only
        # transforms the plan (+cache), and for this Dirichlet body the
        # operator is exactly invariant under rigid motion, so the factors
        # stay valid by the same argument as the near-field cache
        # (design check 2026-09-24; see the KrylovSolver docstring).
        nfcache_gib = parse(Float64, get(ENV, "NFCACHE_MAX_GIB", "500"))
        (isfinite(nfcache_gib) && nfcache_gib > 0) ||
            error("NFCACHE_MAX_GIB must be positive and finite")
        t_precond = @elapsed P = pnl.ILUPreconditioner(rotor; leaf_size=10,
            multipole_acceptance=1.0, max_pattern_entries=8192 * rotor.ncells)
        st = P.stats
        t_ctor = @elapsed solver = pnl.KrylovSolver(rotor; method=:gmres,
            krylov_kw..., preconditioner=P, cache_tree=true,
            cache_nearfield=true, persistent_plan=true,
            nearfield_cache_max_bytes=round(Int, nfcache_gib * 1024^3))
        global setup_split = (; t_precond, t_ctor,
            ilu_tree=st["tree_time"], ilu_lists=st["interaction_list_time"],
            ilu_assembly=st["assembly_time"],
            ilu_factorization=st["factorization_time"])
        return solver, "ilu leaf10/mac1.0;cache_nearfield;persistent_plan;" *
            "nfcache_max_gib=$nfcache_gib"
    end
    if config in ("fgs", "fgmres_fgs")
        # 021 warm-start R4 campaign: FGS_KNOBS_FILE points at a champion
        # config TOML (e.g. benchmark/retained_r4_champion.toml) instead of
        # the Phase-1 fgsprecond winner CSVs. Tolerances are per rung AND per
        # environment (never carried silently), so this path REQUIRES an
        # explicit FGS_TOL_ABS on the submit line, recorded in the CSV knobs.
        # Sweep executor/precision come from the FGSSolver defaults
        # (dagteam+backoff, f64 since 2026-09-24) unless FGS_PRECISION is set;
        # the champion file's own dagteam_precision (f32full, certified COLD
        # R4 only) is deliberately NOT read.
        knobs_file = get(ENV, "FGS_KNOBS_FILE", "")
        if haskey(ENV, "FGS_P")
            # fully explicit FGS knobs (local smokes / ad-hoc arms); all four
            # required, plus the explicit tolerance
            winner = (; p=parse(Int, ENV["FGS_P"]),
                        mac=parse(Float64, ENV["FGS_MAC"]),
                        leaf=parse(Int, ENV["FGS_LEAF"]),
                        inner=parse(Int, ENV["FGS_INNER"]), sweeps=1)
            haskey(ENV, "FGS_TOL_ABS") || error("FGS_P requires FGS_TOL_ABS")
            tol_abs = parse(Float64, ENV["FGS_TOL_ABS"])
            (isfinite(tol_abs) && tol_abs > 0) || error("Invalid FGS_TOL_ABS")
        elseif !isempty(knobs_file)
            isfile(knobs_file) || error("FGS_KNOBS_FILE not found: $knobs_file")
            cf = TOML.parsefile(knobs_file)
            cf["rung"] == rung || error("FGS_KNOBS_FILE rung $(cf["rung"]) != $rung")
            winner = (; p=Int(cf["P"]), mac=Float64(cf["MAC"]),
                        leaf=Int(cf["leaf"]), inner=Int(cf["inner"]),
                        sweeps=1)
            haskey(ENV, "FGS_TOL_ABS") || error("FGS_KNOBS_FILE requires an " *
                "explicit FGS_TOL_ABS (tolerances are per rung AND per " *
                "environment — re-staircase, never carry)")
            tol_abs = parse(Float64, ENV["FGS_TOL_ABS"])
            (isfinite(tol_abs) && tol_abs > 0) || error("Invalid FGS_TOL_ABS")
        else
            winner = stage3_winner()
            winner === nothing && error("no fgsprecond winner for $rung")
            sc = staircase_for(winner.p, winner.mac, winner.leaf, winner.inner)
            i_cross = findfirst(t -> t[4] <= target_rel, sc)
            i_cross === nothing && error("winner staircase never crosses 1e-6")
            tol_abs = margin_tol(sc, i_cross)
        end
        fgs_precision = Symbol(get(ENV, "FGS_PRECISION", "f64"))
        fgs_precision in (:f64, :f32conv, :f32full) || error("Invalid FGS_PRECISION")
        knobs = "p=$(winner.p);mac=$(winner.mac);leaf=$(winner.leaf);" *
                "inner=$(winner.inner);precision=$(fgs_precision)"
        # Phase 3: FGS owns the history/extrapolation machinery already; cold
        # keeps history off entirely so the buffer costs nothing (ruling 8).
        fgs_history = warmstart == "cold" ? 0 : ws_order + 1
        config == "fgs" && return pnl.FGSSolver(rotor;
            expansion_order=winner.p, multipole_acceptance=winner.mac,
            leaf_size=winner.leaf, inner_iterations=winner.inner,
            max_iterations=300, tolerance=tol_abs, rlx=1.0, shrink=true,
            recenter=false, reverse_pass=false,
            dagteam_precision=fgs_precision,
            # SOLVER_VERBOSE=1 turns on FGS's per-iteration residual log and its
            # projected-vs-actual strength dump (diagnostic only; off by default)
            verbose=(get(ENV, "SOLVER_VERBOSE", "0") == "1"),
            solution_history_length=fgs_history,
            project_solution=(warmstart != "cold"),
            project_solution_order=ws_order),
            knobs * ";tol_abs=$tol_abs"
        P = pnl.FGSPreconditioner(rotor; sweeps=winner.sweeps,
            inner_iterations=winner.inner, rlx=1.0, expansion_order=winner.p,
            multipole_acceptance=winner.mac, leaf_size=winner.leaf,
            shrink=true, recenter=false)
        return pnl.KrylovSolver(rotor; method=:fgmres, krylov_kw...,
            preconditioner=P), knobs * ";sweeps=$(winner.sweeps);rtol=$krylov_rtol"
    end
    error("Unknown CONFIG=$(repr(config)); use backslash, krylov_gmres, " *
          "krylov_jacobi, krylov_ilu, fgs, or fgmres_fgs")
end

setup_split = nothing   # filled by configs with a decomposed setup cost
t_setup = @elapsed solver_config, solver_knobs = make_config_solver(config)

# --- persistent-plan priming (021 warm-start campaign, Ryan 2026-09-24) ------
# For persistent_plan configs the FMM plan + near-field cache build lazily on
# the FIRST solve, which would land the build inside step 1's t_solve — and
# step 1 is inside scoring Window A (transients are scored, setup is not). One
# priming solve here moves the build into the reported setup cost. The body
# and solver state are snapshotted and fully restored so the arm's dynamics
# and warm-start semantics are untouched (the priming solution never seeds
# anything).
t_prime = 0.0
if solver_config isa pnl.KrylovSolver && solver_config.persistent_plan
    _strength0 = copy(rotor.strength)
    _velocity0 = copy(rotor.velocity)
    _potential0 = copy(rotor.potential)
    _core0 = rotor.core_size
    # a zero RHS would let Krylov return x=0 without ever applying the
    # operator, so the plan/cache would NOT build. The Dirichlet rhs is
    # -body.potential (zero pre-march), so prime under a synthetic unit
    # potential; Neumann paths read the control-point velocity instead.
    if pnl.has_dirichlet_bc(rotor)
        rotor.potential .= 1.0
    else
        rotor.velocity .= 0
        rotor.velocity[1, :] .= 1.0
    end
    t_prime = @elapsed pnl._solve!(rotor, solver_config)
    solver_config.niter > 0 || error("priming solve did no iterations — " *
        "the plan/cache build cannot be assumed complete")
    rotor.strength .= _strength0
    rotor.velocity .= _velocity0
    rotor.potential .= _potential0
    rotor.core_size = _core0
    # wipe every trace of the priming solution from the warm-start state
    solver_config.have_x_prev = false
    solver_config.x_history_nsaved = 0
    fill!(solver_config.x_prev, 0)
    fill!(solver_config.x_history, 0)
    fill!(solver_config.x0_scratch, 0)
    fill!(solver_config.workspace.x, 0)
    solver_config.niter = 0
    solver_config.solved = false
    pnl.begin_step_solution!(solver_config)
    solver_config.kop.plan_slot[] !== nothing ||
        error("persistent_plan priming solve left no live plan")
    if solver_config.cache_nearfield
        solver_config.kop.plan_slot[][1].nearfield_cache[] !== nothing ||
            error("cache_nearfield set but the priming solve built no cache")
        # B-I3 (033, 2026-10-01): the cache records its own threaded build
        # wall time — persist it so t_prime decomposes into (nfcache build,
        # plan build + priming Krylov iterations)
        _nfc_bt = solver_config.kop.plan_slot[][1].nearfield_cache[].build_time
        setup_split === nothing ||
            (global setup_split = (; setup_split..., prime_nfcache_build=_nfc_bt))
        println("  near-field cache build (threaded): " *
                "$(round(_nfc_bt; digits=2)) s of the priming solve")
    end
    println("persistent-plan priming solve: $(round(t_prime; digits=2)) s " *
            "(plan + near-field cache build; EXCLUDED from per-step costs)")
end
body_solvers = (solver_config,)
println("CONFIG=$config ($solver_knobs), setup $(round(t_setup; digits=2)) s; " *
        "WARMSTART=$warmstart" * (warmstart == "extrap" ? " order=$ws_order" : "") *
        (restart_step >= 0 ? "; RESTART_STEP=$restart_step" : ""))

# Backslash is a direct solve and has no iteration count. Krylov and (since
# Phase 3) FGS both expose niter/solved on the solver object.
#
# `niter` is the LAST inner solve of the step, kept for continuity with the
# Phase 1/2 rows. The Phase 3 warm-start headline metric is `niter_first` — the
# FIRST per-body solve of the step, which is the one that actually consumes the
# step-to-step initial guess. `nsolves` is how many per-body solves the step
# took; for a one-body run it must be 1 on every row, and anything else means
# the timestep re-solved and the warm-start comparison is invalid.
# Generic/direct solvers report -1 (unavailable), never a fabricated success.
_niter_of(s) = s isa pnl.KrylovSolver || s isa pnl.FGSSolver ? s.niter : -1
_niter_first_of(s) = pnl.step_niter_first(s)
_nsolves_of(s) = pnl.step_nsolves(s)
_solved_of(s) = pnl.step_solved(s)

# --- src-free per-step instrumentation --------------------------------------
"Times solve_formulation! of the wrapped production formulation (the body
solve inside _steady_aerodynamics!) and records the solver's iterations."
struct TimedFormulation{F} <: pnl.AbstractSolveFormulation
    inner::F
    t_solve::Vector{Float64}
    niter::Vector{Int}         # last inner solve of the step (Phase 1/2 legacy)
    niter_first::Vector{Int}   # FIRST inner solve of the step (Phase 3 headline)
    nsolves::Vector{Int}       # per-body solves the step took
    solved::Vector{Bool}       # AND of the step's inner convergence flags
    t_project::Vector{Float64} # warm-start guess-construction cost of the
                               # step's last raw solve (== its only solve when
                               # nsolves==1); 0 for cold/prev, NaN unavailable
    # Per-step BC residual statistics (2026-08-25). `bcerr_tol` is the arm's own
    # effective ABSOLUTE acceptance level, so `bcerr_max <= bcerr_tol` is the
    # arm's own promise, checked independently and ~10x more sharply.
    bcerr_max::Vector{Float64}
    bcerr_min::Vector{Float64}
    bcerr_q1::Vector{Float64}
    bcerr_med::Vector{Float64}
    bcerr_q3::Vector{Float64}
    bcerr_rms::Vector{Float64}
    bcerr_tol::Vector{Float64}
    bcerr_eps::Vector{Float64}    # absolute FMM error target requested
    bcerr_cert::Vector{Bool}      # did every M2L pair certify that target?
    t_bcerr::Vector{Float64}      # wall time of the measurement pass (NOT in t_solve)
    n_particles::Vector{Int}      # per step, REPORTED not asserted (chaotic wake)
    tol_of::Base.RefValue{Any}    # (solver, body) -> effective absolute tolerance
    bc_kw::Base.RefValue{Any}     # frozen bc_error! knobs (incl. the `every` stride)
    istep::Base.RefValue{Int}     # steps seen, for the BCERR_EVERY stride
    phi::Vector{Float64}          # per-panel BC residual scratch
    vel_pre::Matrix{Float64}      # body.velocity snapshot taken BEFORE the solve
    # per-step solved-strength snapshots (021 warm-start campaign): the
    # cross-solver solution-agreement deliverable needs the actual per-step
    # solutions, not checksums — FGS and Krylov converge to slightly different
    # fixed points once free wake exists (~2e-3, rigid_motion_tree_reuse_item
    # §5), and the harvest quantifies that delta per step across arms.
    snapshots::Vector{Vector{Float64}}
    snapshot_on::Bool
end
pnl.initialize_formulation(f::TimedFormulation, args...) =
    pnl.initialize_formulation(f.inner, args...)
pnl.formulation_prewake!(f::TimedFormulation, state, systems_tuple) =
    pnl.formulation_prewake!(f.inner, state, systems_tuple)
function pnl.solve_formulation!(f::TimedFormulation, state, systems,
        systems_tuple, wakes_tuple, body_solvers; kwargs...)
    # `bc_error!`'s contract is that `body.velocity` holds the apparent velocity
    # at the control points on entry — it re-derives σ from it. Snapshot it here,
    # before the `:solve` stage runs, and restore it for the measurement. Phase 1
    # does the same thing with its frozen velocity
    # (`phase1_agreement.jl`: `rotor.velocity .= frozen_velocity` before the
    # call); unsteady has no frozen velocity, so the snapshot is the analogue.
    f.vel_pre .= systems_tuple[1].velocity
    t0 = time_ns()
    out = pnl.solve_formulation!(f.inner, state, systems, systems_tuple,
                                 wakes_tuple, body_solvers; kwargs...)
    push!(f.t_solve, (time_ns() - t0) / 1e9)
    # Read AFTER the formulation returns: the step statistics are published at
    # the step boundary inside it (021 Phase 3), not by the raw kernel.
    push!(f.niter, _niter_of(body_solvers[1]))
    push!(f.niter_first, _niter_first_of(body_solvers[1]))
    push!(f.nsolves, _nsolves_of(body_solvers[1]))
    push!(f.solved, _solved_of(body_solvers[1]))
    push!(f.t_project, pnl.step_t_project(body_solvers[1]))

    # --- per-step BC residual -------------------------------------------------
    body = systems_tuple[1]
    f.snapshot_on && push!(f.snapshots,
        copy(vec(body.strength[:, pnl._fgs_solved_strength_index(body)])))
    f.istep[] += 1
    kw0 = f.bc_kw[]
    if (f.istep[] - 1) % kw0.every != 0
        # subsampled: keep the vectors index-aligned with the step, but never
        # fabricate a value for a step that was not measured
        for v in (f.bcerr_max, f.bcerr_min, f.bcerr_q1, f.bcerr_med, f.bcerr_q3,
                  f.bcerr_rms, f.bcerr_tol, f.bcerr_eps, f.t_bcerr)
            push!(v, NaN)
        end
        push!(f.bcerr_cert, false)
        push!(f.n_particles, length(wakes_tuple) > 0 && wakes_tuple[1] !== nothing ?
              pnl.FLOWVPM.get_np(wakes_tuple[1].pfield) : -1)
        return out
    end
    x = vec(body.strength[:, pnl._fgs_solved_strength_index(body)])
    tol = f.tol_of[](body_solvers[1], body)
    body.velocity .= f.vel_pre          # bc_error!'s entry contract (restored by it)
    kw = f.bc_kw[]
    # rms_b=tol with target_rel=1 makes the requested FMM error `safety * tol`:
    # the instrument is sharper than the acceptance level by construction.
    e = bc_error!(body, x; rms_b=tol, target_rel=1.0, safety=kw.safety,
                  max_expansion_order=kw.max_p, multipole_acceptance=kw.mac,
                  leaf_size=kw.leaf, backend=:fmm, phi_out=f.phi)
    ap = sort!(abs.(f.phi))
    np_ = length(ap)
    qi(q) = ap[clamp(round(Int, q * (np_ - 1)) + 1, 1, np_)]
    push!(f.bcerr_max, ap[end]);       push!(f.bcerr_min, ap[1])
    push!(f.bcerr_q1, qi(0.25));       push!(f.bcerr_med, qi(0.50))
    push!(f.bcerr_q3, qi(0.75))
    push!(f.bcerr_rms, e.rel_l2 * tol)   # rel_l2 was normalized by rms_b=tol
    push!(f.bcerr_tol, tol)
    push!(f.bcerr_eps, e.epsilon_requested)
    push!(f.bcerr_cert, e.error_success)
    push!(f.t_bcerr, e.t_eval)
    push!(f.n_particles, length(wakes_tuple) > 0 && wakes_tuple[1] !== nothing ?
          pnl.FLOWVPM.get_np(wakes_tuple[1].pfield) : -1)
    return out
end

# Each solver family promises a different thing; `bcerr_tol` records what THIS
# arm actually promised, in absolute units, so `bcerr_max <= bcerr_tol` is a
# like-for-like check across configs.
#   FGSSolver     — absolute max-abs residual bound (`solver.tolerance`)
#   KrylovSolver  — Krylov.jl stops at atol + rtol*||b||, and for this Dirichlet
#                   body ||b|| is the source-potential norm, which moves every
#                   step; `body.potential` holds it at the solve boundary
#   Backslash     — direct; no iterative tolerance, reported as NaN
function _effective_tol(solver, body)
    solver isa pnl.FGSSolver && return solver.tolerance
    # ||b|| must come from the solver's OWN rhs vector, written at the solve
    # boundary and untouched afterwards: `body.potential` is overwritten by
    # the solve itself, so reading it here collapsed the recorded tolerance
    # to atol and flagged every Krylov row (caught by the 2026-09-24 smoke).
    solver isa pnl.KrylovSolver && return solver.rtol_rhs > 0 ?
        solver.atol + solver.rtol_rhs * LinearAlgebra.norm(solver.rhs) :
        krylov_kw.atol + krylov_kw.rtol * LinearAlgebra.norm(solver.rhs)
    return NaN
end

# bc_error! evaluates at the arm's tuned mac/leaf under a REQUESTED absolute
# error target; max_expansion_order is the dynamic-P cap it may climb to.
_bc_kw = (; safety=bcerr_safety, max_p=bcerr_max_p, mac=mac_t, leaf=leaf_t,
           every=bcerr_every)

pnl.has_dirichlet_bc(rotor) ||
    error("per-step BC check uses bc_error!, which implements the DIRICHLET BC " *
          "metric (control-point potential = φ_σ + G_μ x). This rotor reports " *
          "Neumann; use true_residual! instead and re-derive the tolerance.")

snapshot_on = get(ENV, "SNAPSHOT_STRENGTHS", "0") == "1"
timed_formulation = TimedFormulation(formulation, Float64[], Int[], Int[], Int[],
    Bool[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[],
    Float64[], Float64[], Float64[], Bool[], Float64[], Int[],
    Ref{Any}(_effective_tol), Ref{Any}(_bc_kw), Ref(0),
    zeros(rotor.ncells), zeros(3, rotor.ncells),
    Vector{Float64}[], snapshot_on)

# maneuver! runs once at the start of every step: consecutive-call deltas are
# the full-step wall times (the tail of the last step is closed at sim end)
step_stamps = UInt64[]
function maneuver_timed!(frames, systems, wakes, t)
    push!(step_stamps, time_ns())
    return maneuver!(frames, systems, wakes, t)
end

# --- time march --------------------------------------------------------------
# Phase 3 restarts every arm from ONE shared wake-developed checkpoint, so the
# only variable across arms is the initial guess (and re-marching the 13-rev
# RHPC startup schedule per arm is avoided). RESTART_STEP=-1 keeps the Phase-2
# cold-start path exactly. The checkpoint needs no protect-list entry: the VTK
# sweeper's standing newest-36-timestep retention always leaves the last step
# on disk, and simulate_warmstart! resolves a negative restart_step to the last
# PVD entry — so an arm simply restarts from the newest surviving step.
march! = restart_step >= 0 ? pnl.simulate_warmstart! : pnl.simulate!
_rpath = isempty(restart_path) ? save_path : restart_path
_rname = isempty(restart_name) ? run_name : restart_name
# simulate_warmstart! forwards start_step=restart_step+1 and APPENDS to the
# output it is pointed at, so an arm writing under the checkpoint's own
# path/name would grow (and eventually invalidate) the very state the other
# arms restart from. Refuse rather than silently corrupt the shared source.
restart_step >= 0 && _rpath == save_path && _rname == run_name &&
    error("restart source ($_rpath/$_rname) is this run's own output target; " *
          "set RUN_NAME (or RESTART_NAME/RESTART_PATH) so the arm writes " *
          "somewhere other than the shared checkpoint")
restart_kw = restart_step >= 0 ?
    (; restart_step, restart_path=_rpath, restart_name=_rname) : (;)

sim_wall = @elapsed march!(systems, wakes, frames, maneuver_timed!,
    Uinf, t_range;
    restart_kw...,
    set_Das_eta_kinematic=set_Das_refresh ? init_Das_eta_kinematic : NaN,
    set_Das_min_kinematic_displacement,
    set_Das_kinematic_arc,
    set_Das_refresh,
    monitors,
    formulation=timed_formulation,
    body_solvers, backend, backend_wake,
    wakerow_no_hessian_to_particles,
    body_hessian_to_particles,
    body_gradient_core_size,
    body_on_wake,
    panel_wake_on_particles,
    particle_hessian_self,
    particle_relax,
    bound_strength_rlx,
    verbose=true,
    path=save_path, name=run_name,
)
push!(step_stamps, time_ns())   # close the last step

# --- per-step CT (identity metric, Ryan 2026-08-25) --------------------------
# Chaos makes exact trajectory identity the wrong thing to demand: tiny
# differences compound over many steps, so `n_particles` is REPORTED per step
# (to show the arms stay mostly consistent), never asserted. CT is the physical
# quantity the campaign cares about and the best identity signal available.
# `ForceMonitor` carries a RotorNormalization, so `force` is already CT.
_fmons = filter(m -> m isa pnl.ForceMonitor, collect(monitors))
ct_series = if isempty(_fmons)
    @warn "no ForceMonitor active (RUN_MONITORS=$(ENV["RUN_MONITORS"])) — CT " *
          "identity unavailable; set RUN_MONITORS=true BERNOULLI_ONLY=true"
    Float64[]
else
    f_ct = _fmons[1].force
    [-f_ct[axial_dimension, i] for i in 1:size(f_ct, 2)]
end
_at(v, i) = (i >= 1 && i <= length(v)) ?
    (v[i] isa AbstractFloat && isnan(v[i]) ? "" : v[i]) : ""
# a step skipped by BCERR_EVERY has no verdict at all — blank, not `false`
_measured(i) = i >= 1 && i <= length(timed_formulation.bcerr_max) &&
               !isnan(timed_formulation.bcerr_max[i])
_at_cert(i) = _measured(i) ? timed_formulation.bcerr_cert[i] : ""
# The BC pass sits INSIDE the maneuver!-to-maneuver! window, so `t_step_total`
# and any solve-share derived from it are inflated by the measurement even
# though `t_solve` is not (its timer closes before the pass). Full-timestep
# share is a campaign deliverable (Ryan 2026-08-06), so emit the corrected
# figure alongside the raw one rather than leaving the caller to know this.
_t_step_net(i) = (tt = (step_stamps[i+1] - step_stamps[i]) / 1e9;
                  _measured(i) ? tt - timed_formulation.t_bcerr[i] : tt)

# --- per-step CSV ------------------------------------------------------------
write(joinpath(outdir, "banner.txt"), banner.text * "\n")
csv_path = joinpath(outdir, "unsteady.csv")
header = "rung,mesh_file,n_panels,config,step,t,t_step_total,t_solve,niter," *
    "niter_first,nsolves,t_project," *
    "solved,warmstart,warmstart_order,restart_step,skip_steps," *
    "n_particles_end,strength_checksum,t_setup,t_setup_prime,setup_detail," *
    "solver_knobs,apply_p,apply_mac,apply_leaf," *
    "nt,n_steps,run_monitors,threading_mode,julia_threads,blas_threads," *
    "bcerr_max,bcerr_min,bcerr_q1,bcerr_med,bcerr_q3,bcerr_rms,bcerr_tol," *
    "bcerr_eps,bcerr_certified,t_bcerr,t_step_net,n_particles,CT," *
    "commit,fm_commit,filament_reg,date,notes"
fresh = !isfile(csv_path) || filesize(csv_path) == 0
# Phase 3 widened the schema. Appending new-width rows under an old-width header
# would silently misparse every column past `niter`, so refuse instead: the
# caller must point at a fresh/versioned output path.
if !fresh
    existing_header = open(io -> readline(io), csv_path)
    existing_header == header || error(
        "CSV schema mismatch at $csv_path: the existing header is not the " *
        "current (Phase 3) schema, so appending would mix row widths under one " *
        "header. Write to a fresh/versioned path (set PHASE or move the file " *
        "aside) rather than appending.\n  existing: $existing_header\n  expected: $header")
end
io = open(csv_path, "a")
fresh && println(io, header)

nsim = length(step_stamps) - 1
nsolve = length(timed_formulation.t_solve)
nsim == nsolve || @warn "step count mismatch: $nsim maneuver! deltas vs " *
    "$nsolve solves — check simulate!'s per-step call pattern" nsim nsolve
n_particles_end = pnl.FLOWVPM.get_np(wakes[1].pfield)
# Trajectory-identity guard: cold and warm arms of the same config solve to the
# same tolerance, so they must agree on both of these. Divergence means the warm
# guess changed the answer rather than just reaching it faster.
strength_checksum = sum(abs, rotor.strength)
note = "RHPC setup-only include; t_solve = solve_formulation! wall (incl. " *
    "the per-step Dirichlet source-potential assembly inside solve!); " *
    "t_step from maneuver! deltas; niter now recorded for FGS as well as " *
    "Krylov (Phase 3); niter=-1 means a direct solve. niter is the LAST " *
    "inner solve of the step (legacy diagnostic); niter_first is the FIRST " *
    "and is the warmstart comparison metric; nsolves is the per-body solve " *
    "count of the step (must be 1 for a one-body run; -1 = unavailable); " *
    "solved is the AND of the step's inner convergence flags. Wake-developed " *
    "averages should drop the first skip_steps rows. bcerr_* are order " *
    "statistics of the per-panel Dirichlet BC residual |phi_sigma + G_mu x| " *
    "from ONE bc_error! pass per step; bcerr_tol is THIS arm's own effective " *
    "absolute acceptance level, so bcerr_max <= bcerr_tol is the arm's promise " *
    "checked independently. The pass requests an absolute FMM error of " *
    "bcerr_eps = BCERR_SAFETY * bcerr_tol, so the instrument is ~10x sharper " *
    "than the acceptance level by construction; bcerr_certified=false means " *
    "some M2L pair could not certify that bound and the row is UNCERTIFIED. " *
    "The pass is EXCLUDED from t_solve, but NOT from t_step_total (it lies " *
    "inside the maneuver!-to-maneuver! window); use t_step_net = t_step_total " *
    "- t_bcerr for full-timestep share. Its own cost is t_bcerr (~62% of t_solve " *
    "at R1 — it certifies to 10x below tol_abs, so the FMM climbs its order). " *
    "BCERR_EVERY>1 subsamples it; skipped steps write EMPTY bcerr_* cells, " *
    "never a fabricated value. n_particles and " *
    "CT are per-step identity REPORTING, not assertions: the wake is chaotic, " *
    "so small differences compound and exact trajectory identity is not " *
    "expected. CT requires RUN_MONITORS=true (default under PHASE=phase3)."
for i in 1:min(nsim, nsolve)
    println(io, join(_csv_cell.([rung, msh_ladder, rotor.ncells, config, i,
        t_range[i], (step_stamps[i+1] - step_stamps[i]) / 1e9,
        timed_formulation.t_solve[i], timed_formulation.niter[i],
        timed_formulation.niter_first[i], timed_formulation.nsolves[i],
        _at(timed_formulation.t_project, i),
        timed_formulation.solved[i], warmstart, ws_order, restart_step, skip_steps,
        i == nsim ? n_particles_end : "",
        i == nsim ? strength_checksum : "", i == 1 ? t_setup : "",
        i == 1 ? t_prime : "",
        i == 1 && setup_split !== nothing ?
            join(("$k=$v" for (k, v) in pairs(setup_split)), ";") : "",
        solver_knobs, p_t, mac_t, leaf_t, parse(Int, ENV["NT"]), nsim,
        ENV["RUN_MONITORS"], banner.threading_mode, banner.julia_threads,
        banner.blas_threads,
        _at(timed_formulation.bcerr_max, i), _at(timed_formulation.bcerr_min, i),
        _at(timed_formulation.bcerr_q1, i), _at(timed_formulation.bcerr_med, i),
        _at(timed_formulation.bcerr_q3, i), _at(timed_formulation.bcerr_rms, i),
        _at(timed_formulation.bcerr_tol, i), _at(timed_formulation.bcerr_eps, i),
        _at_cert(i), _at(timed_formulation.t_bcerr, i), _t_step_net(i),
        _at(timed_formulation.n_particles, i), _at(ct_series, i),
        banner.commit, banner.fm_commit, banner.filament_reg, time_string(),
        i == 1 ? note : ""]), ","))
end
close(io)

# --- per-step solution snapshots (021 warm-start campaign) --------------------
# Raw Float64, ncells x nsteps, column-major (one column per step, in step
# order), with a sidecar TOML naming the shape — the harvest computes per-step
# cross-arm solution deltas from these.
if snapshot_on && !isempty(timed_formulation.snapshots)
    snap_base = joinpath(outdir, "$(run_name)_strength_snapshots")
    open(snap_base * ".bin", "w") do sio
        foreach(s -> write(sio, s), timed_formulation.snapshots)
    end
    open(snap_base * ".toml", "w") do sio
        TOML.print(sio, Dict("ncells" => rotor.ncells,
            "nsteps" => length(timed_formulation.snapshots),
            "eltype" => "Float64", "layout" => "column-major, one column per step",
            "strength_column" => pnl._fgs_solved_strength_index(rotor),
            "config" => config, "warmstart" => warmstart,
            "warmstart_order" => ws_order, "run_name" => run_name,
            # leg discrimination for the warm-start campaign harvest: a
            # restarted leg's column i is GLOBAL step restart_step + i
            "restart_step" => restart_step,
            "n_steps_env" => n_steps_env))
    end
    println("strength snapshots: $(length(timed_formulation.snapshots)) steps " *
            "-> $(snap_base).bin")
end

# Net of the measurement, for the same reason as `t_step_net` above.
t_bcerr_total = sum(x -> isnan(x) ? 0.0 : x, timed_formulation.t_bcerr; init=0.0)
sim_wall_net = sim_wall - t_bcerr_total
t_solve_share = sum(timed_formulation.t_solve) / max(sim_wall_net, eps())
n_unconverged = count(!, timed_formulation.solved)
println("\n$rung/$config unsteady [warmstart=$warmstart" *
        (warmstart == "extrap" ? "/order$ws_order" : "") * "]: " *
        "$nsim steps, sim wall $(round(sim_wall_net; digits=1)) s net " *
        "($(round(sim_wall; digits=1)) s raw, incl. $(round(t_bcerr_total; digits=1)) s " *
        "of BC measurement), solve share $(round(100 * t_solve_share; digits=1))% " *
        "(of net), " *
        "$(n_particles_end) particles at end, checksum $strength_checksum. " *
        "Rows appended to $csv_path")

# Loud, because this is the failure mode the accuracy guard exists to catch: a
# warm guess that lets steps stop before they have actually converged.
n_unconverged == 0 ||
    @warn "$n_unconverged of $nsim steps did NOT meet the solver tolerance — " *
          "this row set is not a valid warmstart comparison" config warmstart

# One-body solve-count guard (021 Phase 3). This run is a 1-tuple, so a valid
# warmstart comparison needs exactly ONE per-body solve per step: more means the
# timestep re-solved the same system and `niter_first` is no longer measuring a
# step-to-step warm guess. `nsolves == -1` is "unavailable" (Backslash and other
# generic solvers), which is NOT a passed check — it just isn't applicable.
# Deliberately not generalized to multibody: genuine block-GS outer counts vary
# per timestep there.
if nsolve > 0
    win_all = (min(skip_steps, nsolve) + 1):nsolve   # wake-developed rows
    ns = timed_formulation.nsolves[win_all]
    if isempty(ns) || first(ns) < 0
        println("nsolves: unavailable for CONFIG=$config (direct/generic solver) " *
                "— the one-body solve-count check does not apply")
    else
        bad = count(!=(1), ns)
        bad == 0 ?
            println("nsolves == 1 on all $(length(ns)) wake-developed rows " *
                    "(one-body step contract holds)") :
            @warn "$bad of $(length(ns)) wake-developed rows have nsolves != 1 " *
                  "— the timestep solved the body more than once, so " *
                  "niter_first is NOT a step-to-step warmstart measurement and " *
                  "this row set is an INVALID warmstart comparison" config warmstart extrema(ns)
    end
end

# Wake-developed window, excluding the post-restart history-fill transient.
if nsim > skip_steps
    win = (skip_steps + 1):min(nsim, nsolve)
    iters = timed_formulation.niter[win]
    firsts = timed_formulation.niter_first[win]
    solves = timed_formulation.t_solve[win]
    if !isempty(firsts) && first(firsts) >= 0
        # niter_first is the headline: it is the solve that consumes the
        # step-to-step initial guess. niter (last inner solve) is legacy.
        println("wake-developed window (steps $(first(win))-$(last(win))): " *
                "niter_first mean $(round(sum(firsts)/length(firsts); digits=2)) " *
                "[min $(minimum(firsts)), max $(maximum(firsts))], " *
                "t_solve mean $(round(sum(solves)/length(solves); digits=4)) s" *
                (isempty(iters) || first(iters) < 0 ? "" :
                 "; legacy niter (last inner solve) mean " *
                 "$(round(sum(iters)/length(iters); digits=2))"))
    end
end
