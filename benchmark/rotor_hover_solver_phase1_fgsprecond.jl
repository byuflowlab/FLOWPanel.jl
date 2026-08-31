#=##############################################################################
BRAINSTORM 021 Phase 1 Stage 3 — FGS-preconditioner selection via the
τ-ladder (Ryan 2026-08-15 ruling 3; phase_01_tuning_plan.md Stage 3).

For each τ rung of fgstune_selected.csv (Stage 2 output, same rung): wrap the
τ-tuned FGS config as an FGSPreconditioner with a FIXED sweep count
(= the ladder's precond_sweeps; fixed counts + zero seed keep each apply a
linear map — the clean FGMRES regime; tolerance-stopped applies are banned),
then measure end-to-end preconditioned-FGMRES wall time to BC 1e-6, with the
outer Krylov apply on the rung's tuned apply knobs (tune_fmm output). The
winning τ fixes the SINGLE FGS knob set shared between solver and
preconditioner roles (ruling 3).

Judgment (logged): the preconditioner is built with shrink=true to match the
solver-role FGS trees (ruling 3 shared-knob intent); Phase 0's fgs_check used
the ctor default shrink=false for the preconditioner.

Run AFTER fgstune on the same rung:
  RUNG=R1 EXPECT_JULIA_THREADS=4 THREADING_MODE=multi \
    julia --project -t 4 benchmark/rotor_hover_solver_phase1_fgsprecond.jl
Appends fgsprecond.csv under benchmark/results/phase1/<mode>/.
=###############################################################################

include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "phase1_case.jl"))

target_rel = 1e-6
safety = 0.1
p_t, mac_t, leaf_t = TUNED[rung]

bc_fmm(x) = bc_error!(rotor, x; rms_b, target_rel, safety,
                      max_expansion_order=20, multipole_acceptance=mac_t,
                      leaf_size=leaf_t, backend=:fmm)

# ---- read the Stage 2 ladder (latest row per τ for this rung) ----
sel_path = joinpath(results_dir, "fgstune_selected.csv")
isfile(sel_path) || error("Run fgstune first: $sel_path not found")
ladder = Dict{Float64,NamedTuple}()
let lines = readlines(sel_path)
    cols = split(lines[1], ",")
    ix = Dict(String(c) => i for (i, c) in enumerate(cols))
    for line in lines[2:end]
        isempty(strip(line)) && continue
        cells = split(line, ",")
        cells[ix["rung"]] == rung || continue
        tau = parse(Float64, cells[ix["tau"]])
        ladder[tau] = (;                       # later rows overwrite → latest
            p=parse(Int, cells[ix["p"]]),
            mac=parse(Float64, cells[ix["mac"]]),
            leaf=parse(Int, cells[ix["leaf"]]),
            inner=parse(Int, cells[ix["inner"]]),
            sweeps=parse(Int, cells[ix["precond_sweeps"]]))
    end
end
isempty(ladder) && error("No fgstune_selected.csv rows for rung $rung")

# Ryan ruling 2026-08-17: the shared FGS knob set = the τ=1e-6-tuned config;
# only the preconditioner SWEEP COUNT is selected end-to-end. With
# SWEEP_LADDER_1E6=1 the τ-ladder is rebuilt as {τ => 1e-6 config with
# sweeps = (that config's own staircase crossing at τ) − 1}.
if get(ENV, "SWEEP_LADDER_1E6", "0") == "1"
    w = ladder[1e-6]
    st_cols, st_rows = let lines = readlines(joinpath(results_dir, "fgstune_staircase.csv"))
        (Dict(String(c) => i for (i, c) in enumerate(split(lines[1], ","))),
         [split(l, ",") for l in lines[2:end] if !isempty(strip(l))])
    end
    pts = Dict{Int,Float64}()          # iter => bc (latest rows win)
    for c in st_rows
        c[st_cols["rung"]] == rung || continue
        parse(Int, c[st_cols["p"]]) == w.p || continue
        parse(Float64, c[st_cols["mac"]]) == w.mac || continue
        parse(Int, c[st_cols["leaf"]]) == w.leaf || continue
        parse(Int, c[st_cols["inner"]]) == w.inner || continue
        pts[parse(Int, c[st_cols["iter"]])] = parse(Float64, c[st_cols["bc_rel_l2"]])
    end
    its = sort(collect(keys(pts)))
    empty!(ladder)
    for tau in (1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6)
        i = findfirst(k -> pts[k] <= tau, its)
        i === nothing && continue
        ladder[tau] = (; w.p, w.mac, w.leaf, w.inner, sweeps=its[i] - 1)
    end
    println("sweep ladder of the 1e-6 config p=$(w.p)/$(w.mac)/$(w.leaf)/$(w.inner): " *
            join(["τ=$t→$(ladder[t].sweeps)" for t in sort(collect(keys(ladder)); rev=true)], ", "))
end

header = "rung,mesh_file,n_panels,tau,p,mac,leaf,inner,sweeps,t_setup," *
    "t_solve_min,k_reps,niter,converged,bc_rel_l2,bc_certified,meets_target," *
    "apply_p,apply_mac,apply_leaf,rms_b,radius_tol,threading_mode," *
    "julia_threads,commit,fm_commit,date,notes"
csv_path = joinpath(results_dir, "fgsprecond.csv")
fresh = !isfile(csv_path) || filesize(csv_path) == 0
io = open(csv_path, "a")
fresh && println(io, header)

backend_apply = pnl.FastMultipoleBackend(; expansion_order=p_t,
    multipole_acceptance=mac_t, leaf_size=leaf_t)

best_tau = NaN; best_t = Inf
for tau in sort(collect(keys(ladder)); rev=true)     # coarse → fine
    cfg = ladder[tau]
    cfg.sweeps >= 1 || (println("τ=$tau: sweeps=$(cfg.sweeps) < 1, skipped"); continue)
    println("--- τ=$tau: precond p=$(cfg.p) MAC=$(cfg.mac) leaf=$(cfg.leaf) " *
            "inner=$(cfg.inner) sweeps=$(cfg.sweeps) ---")
    t_setup = @elapsed begin
        P = pnl.FGSPreconditioner(rotor; sweeps=cfg.sweeps,
            inner_iterations=cfg.inner, rlx=1.0, expansion_order=cfg.p,
            multipole_acceptance=cfg.mac, leaf_size=cfg.leaf,
            shrink=true, recenter=false)
        solver = pnl.KrylovSolver(rotor; method=:fgmres, itmax=500,
            atol=1e-14, rtol=target_rel, memory=50, backend=backend_apply,
            preconditioner=P)
    end
    t_solve, _ = min_of_k(() -> pnl._solve!(rotor, solver);
                          k=k_reps, warmup=1, setup! = reset_cold!)
    x = copy(vec(rotor.strength[:, solution_column]))
    e = bc_fmm(x)
    meets = e.error_success && e.rel_l2 <= target_rel
    meets && t_solve < best_t && (global best_t = t_solve; global best_tau = tau)
    println("  -> niter=$(solver.niter) t=$(round(t_solve; digits=2))s " *
            "bc=$(e.rel_l2) $(meets ? "MEETS 1e-6" : "MISSES 1e-6")")
    println(io, join(_csv_cell.([rung, msh_name, rotor.ncells, tau, cfg.p,
        cfg.mac, cfg.leaf, cfg.inner, cfg.sweeps, t_setup, t_solve, k_reps,
        solver.niter, solver.solved, e.rel_l2, e.error_success, meets,
        p_t, mac_t, leaf_t, rms_b, pnl.FMM_RADIUS_TOL[],
        banner.threading_mode, banner.julia_threads, banner.commit,
        banner.fm_commit, time_string(),
        "fgmres rtol=1e-6 atol=1e-14 memory=50; apply=tuned; fixed-sweep " *
        "linear applies (zero seed); shrink=true (shared-knob judgment); " *
        "t excludes N^2 source assembly (restored)"]), ","))
    flush(io)
    P = nothing; solver = nothing; GC.gc()
end
close(io)

println("\n$rung preconditioner τ-ladder done. Best: τ=$best_tau " *
        "($(round(best_t; digits=2)) s end-to-end). Rows in $csv_path")
