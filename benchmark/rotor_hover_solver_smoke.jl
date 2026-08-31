#=##############################################################################
BRAINSTORM 021 Phase 0 W5 — solver-roster smoke harness.

Runs all seven roster configurations (control doc ruling 1) on one small DJI9443
mesh in a `rotor_hover.jl`-style setup, and emits schema-valid CSVs
(`runs.csv` + per-run `history_<run_id>.csv` sidecars) under
`benchmark/results/smoke/<threading_mode>/`. Availability proof + CSV
plumbing only — published numbers come from Phases 2+ on HPC.

Roster:
  a  backslash_coupled   BackslashCoupled((rotor,)), monolithic
  b  backslash_iterative (Backslash(rotor),) through the tuple block-GS path
  c  krylov_gmres        KrylovSolver, unpreconditioned
  d  krylov_jacobi       KrylovSolver + block-Jacobi preconditioner
  e  fgs                 FGSSolver
  f  fgmres_fgs          KrylovSolver(:fgmres) + FGSPreconditioner
  g  krylov_ilu          KrylovSolver(:gmres) + FMM-direct-list ILU(0)

Run (never mix threading modes in one comparison — ruling 6):
  EXPECT_JULIA_THREADS=1 THREADING_MODE=single julia --project -t 1 benchmark/rotor_hover_solver_smoke.jl
  EXPECT_JULIA_THREADS=6 THREADING_MODE=multi  julia --project -t 6 benchmark/rotor_hover_solver_smoke.jl
(64-thread multi mode is the same script on HPC; local runs stay ≤ 6 threads.)

ENV knobs: SMOKE_MESH (mesh path), N_STEPS (default 3), K_REPS (default 3),
KRYLOV_ITMAX, KRYLOV_RTOL, HARDWARE_TAG.
=###############################################################################

import FLOWPanel as pnl
include(joinpath(pnl.examples_path, "helper_functions.jl"))
include(joinpath(pnl.examples_path, "dji9443_trailing_edge.jl"))
include(joinpath(@__DIR__, "common.jl"))

using FLOWPanel.FastMultipole.StaticArrays
using LinearAlgebra: norm, cross, lu!

################################################################################
# Banner + output layout
################################################################################

banner = assert_and_banner()

results_dir = joinpath(@__DIR__, "results", "smoke", banner.threading_mode)
mkpath(results_dir)
write(joinpath(results_dir, "banner.txt"), banner.text * "\n")

################################################################################
# Case parameters (mirrors examples/rotor_hover.jl)
################################################################################

magVinf = 0.0001
AOA     = 0.0
rho     = 1.225
RPM     = 6000
Vinf    = magVinf * [0.0, -cosd(AOA), sind(AOA)]
R       = 0.119
nt      = 36
dt      = 60 / RPM / nt
n_steps = parse(Int, get(ENV, "N_STEPS", "3"))
k_reps  = parse(Int, get(ENV, "K_REPS", "3"))   # smoke only; published runs use k ≥ 5 (ruling 5)

# Smoke-scale Krylov settings: the unpreconditioned operator is ill-conditioned
# on this case (residual 17293 → 56 over 200 GMRES iterations observed), so
# smoke uses a loose rtol to finish locally. Phase 1 calibrates the real,
# frozen tolerances (matched true-residual levels per decision_rules.md).
krylov_itmax = parse(Int, get(ENV, "KRYLOV_ITMAX", "400"))
krylov_rtol  = parse(Float64, get(ENV, "KRYLOV_RTOL", "1e-4"))
krylov_atol  = 1e-6
t_range = range(0.0, step=dt, length=n_steps)

core_size = R * 1e-3
kernelcutoff = R * 1e-13
init_Das_eta_kinematic = 0.2
wake_core_size = 1e-3

msh_file = get(ENV, "SMOKE_MESH",
    joinpath(pnl.examples_path, "data", "dji9443_20260722_40_41_uncapped.msh"))
watertight = occursin("capped", basename(msh_file)) && !occursin("uncapped", basename(msh_file))
watertight && error("Smoke harness assumes an uncapped (Neumann) mesh; got $msh_file")

println("Mesh: $msh_file")
te_indices_1, te_indices_2 = find_dji9443_trailing_edge_indices(msh_file; watertight)

msh = pnl.read_gmsh(msh_file)
nodes0, cells0 = pnl.meshes2nodes_cells(msh)
nodes0 .*= R / maximum(nodes0[1, :])

# Shedding from the *constructed* body's cells (CLAUDE.md critical invariant):
# build a noshedding body first, trace shedding on its nodes/cells, rebuild.
kernel = pnl.VortexRing
DBC = false
rotor0 = pnl.RigidWakeBody{kernel}(nodes0, cells0, pnl.noshedding;
            core_size, kernelcutoff, semiinfinite_wake=false, watertight, DBC)
shedding1 = pnl.calc_shedding_from_seed(
    rotor0.nodes, rotor0.cells, te_indices_1[1], te_indices_1[2];
    end_node=te_indices_1[3], normal_jump_tol=0.2, max_turn_angle=pi/3)
shedding2 = pnl.calc_shedding_from_seed(
    rotor0.nodes, rotor0.cells, te_indices_2[1], te_indices_2[2];
    end_node=te_indices_2[3], normal_jump_tol=0.2, max_turn_angle=pi/3)
const ROTOR_NODES = copy(rotor0.nodes)
const ROTOR_CELLS = copy(rotor0.cells)
const ROTOR_SHEDDING = [shedding1, shedding2]

"Fresh rotor body + reference frame, Das initialized (pre-solver, per rotor_hover.jl)."
function build_case()
    rotor = pnl.RigidWakeBody{kernel}(ROTOR_NODES, ROTOR_CELLS, ROTOR_SHEDDING;
                core_size, kernelcutoff, semiinfinite_wake=false, watertight,
                ensure_winding=true, DBC)
    frames = pnl.ReferenceFrame(rotor;
        origin = SVector{3}(0.0, 0.0, 0.0),
        v = SVector{3}(0.0, 0.0, 0.0),
        ω_axis = SVector{3}(1.0, 0.0, 0.0),
        ω = -2*pi * RPM/60,
        R = SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
        name = "vehicle",
        child_index = Int[],
        dependent_index = [1])
    pnl.initialize_Das!((rotor,), frames, t -> Vinf, t_range[1], t_range[2] - t_range[1];
        set_Das_eta_kinematic=init_Das_eta_kinematic)
    return rotor, frames
end

"Freeze the rotor-like RHS: apparent velocity Vinf - ω×r at each control point."
function set_frozen_velocity!(rotor)
    ω_vec = SVector{3}(-2*pi * RPM/60, 0.0, 0.0)
    pnl.calc_normals!(rotor)
    pnl.calc_controlpoints!(rotor)
    for i in 1:rotor.ncells
        r = SVector{3}(rotor.controlpoints[1, i], rotor.controlpoints[2, i], rotor.controlpoints[3, i])
        rotor.velocity[:, i] .= Vinf .- cross(ω_vec, r)
    end
    return nothing
end

backend_sim = pnl.FastMultipoleBackend()
const SOLVE_SUMMARY = Dict{String,NamedTuple}()

################################################################################
# Per-config runner
################################################################################

"""
Run one roster config end-to-end: timed setup, min-of-k frozen-RHS isolated
solve, history capture, true residual (shared evaluator, DirectBackend per
ruling 3), then an `n_steps`-step `simulate!` as the availability proof.
Writes one runs.csv row + a history sidecar; returns the final Bernoulli CT.
"""
function run_config!(io_runs, config::String;
        make_solver,          # rotor -> (solver, setup_times::NamedTuple, notes::String)
        solve_once!,          # (rotor, solver) -> nothing  (one isolated solve)
        capture_history! =nothing,  # (rotor, solver, hist) -> nothing, or nothing
        iterations_of =(solver, hist) -> length(hist),
        body_solvers_for =(solver) -> solver,
        cold_tripwire =nothing)     # (solver, hist) -> nothing; error on warm-state artifacts

    println("\n=== config $config ===")
    rotor, frames = build_case()
    set_frozen_velocity!(rotor)

    solver, setup, notes = make_solver(rotor)

    # COLD-STATE RESET, run untimed before EVERY solve below: FGS seeds from
    # the current body.strength (and block-GS converges on Δstrength), so
    # without this every rep after warmup is a warm no-op, not a solve
    # (2026-08-12 review finding 1).
    reset_cold! = () -> (rotor.strength .= 0; nothing)

    # First ever application includes cold compilation/cache effects and is
    # reported separately from the warmed min-of-k steady solve.
    reset_cold!()
    t_cold_first = @elapsed solve_once!(rotor, solver)

    # isolated frozen-RHS solve: min over k_reps after 1 warmup (rulings 5/7)
    t_solve_min, _ = min_of_k(() -> solve_once!(rotor, solver);
                              k=k_reps, warmup=1, setup! = reset_cold!)
    reset_cold!()
    alloc_solve = @allocated solve_once!(rotor, solver)

    # RHS-assembly cost, timed separately (ruling 7 / review finding 3)
    b = zeros(rotor.ncells)
    t_rhs, _ = min_of_k(() -> pnl.assemble_rhs!(b, rotor); k=k_reps, warmup=1)

    # convergence history (021 ruling 4) — from a cold solve
    hist = pnl.ConvergenceHistory()
    if capture_history! !== nothing
        reset_cold!()
        capture_history!(rotor, solver, hist)
    end
    cold_tripwire === nothing || cold_tripwire(solver, hist)

    # true residual through the shared evaluator (021 ruling 3)
    x = copy(vec(rotor.strength[:, 1]))
    r = zeros(rotor.ncells)
    rms, rmax = pnl.true_residual!(r, rotor, x, b; backend=pnl.DirectBackend())
    println("  true residual: rms = $rms, max = $rmax")

    # capture the benchmark solve's iteration count NOW — the simulate! below
    # runs warmstarted solves that overwrite solver.niter
    iterations = iterations_of(solver, hist)
    setup_total = haskey(setup, :t_precond) && setup.t_precond !== nothing ?
        setup.t_precond :
        sum((v for v in values(setup) if v !== nothing); init=0.0)
    SOLVE_SUMMARY[config] = (; iterations, t_steady=t_solve_min,
        t_cold_first, setup_total, rms_residual=rms, max_residual=rmax)

    run_id = "$(config)_$(banner.threading_mode)"
    length(hist) > 0 && write_history_csv(joinpath(results_dir, "history_$(run_id).csv"), hist)

    # Capture the benchmark-leg solver settings NOW: body_solvers_for below may
    # mutate the solver for the simulate! availability leg (e.g. warmstart=true),
    # and the CSV row must describe the timed cold solves (ruling 10).
    solver_settings_str = settings_string(pnl._solver_metadata_dict(solver))

    # availability proof: n_steps of the full unsteady pipeline
    wake = pnl.PanelWake(rotor; nwakerows=n_steps + 1, core_size=wake_core_size)
    pressure = pnl.PressureBernoulli(rho; unsteady=true, allow_partial=true,
                    correct_kuttacondition=false, backend=backend_sim)
    force = pnl.ForceMonitor(length(t_range), 1; i_frame=1,
                    normalization=pnl.RotorNormalization(rho, 2*R, 1),
                    correct_kuttacondition=false, verbose=false)
    maneuver!(frames, systems, wakes, t) = nothing
    t_sim = @elapsed pnl.simulate!((rotor,), (wake,), frames, maneuver!, t -> Vinf, t_range;
        set_Das_eta_kinematic=NaN,
        monitors=(pressure, force),
        body_solvers=body_solvers_for(solver),
        backend=backend_sim, verbose=false,
        path=joinpath(results_dir, "sim_$(config)"), name="smoke_$(config)")
    CT = force.force[2, end]
    println("  simulate!($(n_steps) steps): $(round(t_sim; digits=2)) s, CT[end] = $CT")

    write_run_row!(io_runs;
        run_id, phase="phase0_smoke", solver_config=config,
        mesh_file=basename(msh_file), n_panels=rotor.ncells,
        threading_mode=banner.threading_mode,
        julia_threads=banner.julia_threads, blas_threads=banner.blas_threads,
        t_assembly=get(setup, :t_assembly, nothing),
        t_factorize=get(setup, :t_factorize, nothing),
        t_tree=get(setup, :t_tree, nothing),
        t_precond=get(setup, :t_precond, nothing),
        t_rhs, t_solve_min, k_reps,
        iterations,
        rms_residual=rms, max_residual=rmax,
        mem_state_bytes=solver_state_bytes(solver), alloc_solve_bytes=alloc_solve,
        commit=banner.commit, fm_commit=banner.fm_commit,
        julia_version=banner.julia_version, hardware_tag=banner.hardware_tag,
        filament_reg=banner.filament_reg,
        solver_settings=solver_settings_str,
        backend_settings=settings_string(pnl._backend_metadata_dict(backend_sim)),
        notes=notes * " | setup_total $(setup_total)s | cold_first $(t_cold_first)s" *
              " | sim $(n_steps) steps $(round(t_sim; digits=2))s CT $(CT)" *
              " | sim leg ran warmstarted (body_solvers_for)")

    return CT
end

################################################################################
# The seven roster configs
################################################################################

runs_path = joinpath(results_dir, "runs.csv")
CTs = Dict{String,Float64}()

open(runs_path, "w") do io
    println(io, runs_csv_header())

    # (a) BackslashCoupled — monolithic; ctor is trivial, the assembly+lu!
    # happens on the first solve (t_build; combined assembly+factorization)
    CTs["backslash_coupled"] = run_config!(io, "backslash_coupled";
        make_solver = rotor -> begin
            solver = pnl.BackslashCoupled((rotor,))
            pnl.solve!((rotor,), solver)   # triggers the one-time build
            (solver, (; t_assembly=solver.t_build),
             "t_assembly = coupled t_build (assembly+factorize combined)")
        end,
        solve_once! = (rotor, solver) -> pnl.solve!((rotor,), solver),
        iterations_of = (solver, hist) -> nothing,
        body_solvers_for = solver -> solver)

    # (b) Backslash through the tuple block-GS path (what simulate! uses for
    # per-body solvers). NOTE: with a single body there is no coupling to
    # iterate, so the outer loop breaks after the first block solve
    # (src/FLOWPanel_solver.jl:2000, 021 Phase 3) — 1 inner solve per call and
    # exactly one ConvergenceHistory entry. Before that fix this path ran a
    # second identical solve to observe Δ=0, so pre-2026-08-23 rows carry
    # roughly double t_solve_min/alloc_solve.
    CTs["backslash_iterative"] = run_config!(io, "backslash_iterative";
        make_solver = rotor -> begin
            G = zeros(rotor.ncells, rotor.ncells)
            t_assembly = @elapsed pnl._G!(G, rotor, rotor; core_size=rotor.core_size_panel)
            t_factorize = @elapsed lu!(G)
            solver = pnl.Backslash(rotor)
            (solver, (; t_assembly, t_factorize),
             "tuple block-GS: 1 inner solve per call on a single body (one-body " *
             "early break; 021 Phase 3)")
        end,
        solve_once! = (rotor, solver) -> pnl.solve!((rotor,), (solver,); backend=(backend_sim,)),
        capture_history! = (rotor, solver, hist) ->
            pnl.solve!((rotor,), (solver,); backend=(backend_sim,), history=hist),
        iterations_of = (solver, hist) -> length(hist),
        body_solvers_for = solver -> (solver,))

    # (c) Krylov GMRES, unpreconditioned
    CTs["krylov_gmres"] = run_config!(io, "krylov_gmres";
        make_solver = rotor -> begin
            solver = pnl.KrylovSolver(rotor; method=:gmres, itmax=krylov_itmax,
                atol=krylov_atol, rtol=krylov_rtol, backend=backend_sim)
            (solver, (;), "")
        end,
        solve_once! = (rotor, solver) -> pnl.solve!(rotor, solver),
        capture_history! = (rotor, solver, hist) -> begin
            solver.record_history = true
            pnl.solve!(rotor, solver)
            hist.metric = solver.history.metric
            append!(hist.iter, solver.history.iter)
            append!(hist.t_ns, solver.history.t_ns)
            append!(hist.residual_internal, solver.history.residual_internal)
            hist.t0_ns = solver.history.t0_ns
            solver.record_history = false
        end,
        iterations_of = (solver, hist) -> solver.niter,
        # warmstart across timesteps for the availability run (off during the
        # cold isolated benchmark above; Phase 3 studies this properly)
        body_solvers_for = solver -> (solver.warmstart = true; (solver,)))

    # (d) Krylov GMRES + block-Jacobi preconditioner
    CTs["krylov_jacobi"] = run_config!(io, "krylov_jacobi";
        make_solver = rotor -> begin
            t_precond = @elapsed P = pnl.FastMultipole.JacobiPreconditioner((rotor,); cell_size=R/4)
            solver = pnl.KrylovSolver(rotor; method=:gmres, itmax=krylov_itmax,
                atol=krylov_atol, rtol=krylov_rtol, backend=backend_sim, preconditioner=P)
            (solver, (; t_precond), "jacobi cell_size = R/4")
        end,
        solve_once! = (rotor, solver) -> pnl.solve!(rotor, solver),
        capture_history! = (rotor, solver, hist) -> begin
            solver.record_history = true
            pnl.solve!(rotor, solver)
            hist.metric = solver.history.metric
            append!(hist.iter, solver.history.iter)
            append!(hist.t_ns, solver.history.t_ns)
            append!(hist.residual_internal, solver.history.residual_internal)
            hist.t0_ns = solver.history.t0_ns
            solver.record_history = false
        end,
        iterations_of = (solver, hist) -> solver.niter,
        # warmstart across timesteps for the availability run (off during the
        # cold isolated benchmark above; Phase 3 studies this properly)
        body_solvers_for = solver -> (solver.warmstart = true; (solver,)))

    # (e) FGSSolver (knobs from the stubbed block in examples/rotor_hover.jl)
    CTs["fgs"] = run_config!(io, "fgs";
        make_solver = rotor -> begin
            t_tree = @elapsed solver = pnl.FGSSolver(rotor;
                max_iterations=500, tolerance=1e-6, rlx=1.0,
                expansion_order=10, multipole_acceptance=0.4, leaf_size=150,
                shrink=true, recenter=false, inner_iterations=20,
                reverse_pass=false, verbose=false)
            (solver, (; t_tree), "")
        end,
        solve_once! = (rotor, solver) -> pnl.solve!(rotor, solver),
        capture_history! = (rotor, solver, hist) -> pnl._solve_history!(rotor, solver, hist),
        iterations_of = (solver, hist) -> length(hist),
        body_solvers_for = solver -> (solver,),
        # cold tripwire (2026-08-12 review finding 1): a genuinely cold FGS
        # solve must START far from converged — a first-iteration residual at
        # tolerance means the timed reps were warm no-ops, not solves
        cold_tripwire = (solver, hist) -> begin
            length(hist) >= 1 || error("fgs: empty convergence history")
            hist.residual_internal[1] > 100 * solver.tolerance || error(
                "fgs: warm-state artifact — first cold-solve residual " *
                "$(hist.residual_internal[1]) is already at tolerance $(solver.tolerance)")
        end)

    # (f) FGMRES + FGS preconditioner (config 1f, implemented Phase 0 W2)
    CTs["fgmres_fgs"] = run_config!(io, "fgmres_fgs";
        make_solver = rotor -> begin
            t_precond = @elapsed P = pnl.FGSPreconditioner(rotor;
                sweeps=1, inner_iterations=2, expansion_order=10,
                multipole_acceptance=0.4, leaf_size=150)
            solver = pnl.KrylovSolver(rotor; method=:fgmres, itmax=krylov_itmax,
                atol=krylov_atol, rtol=krylov_rtol, backend=backend_sim, preconditioner=P)
            (solver, (; t_precond), "fgs precond sweeps=1 inner=2")
        end,
        solve_once! = (rotor, solver) -> pnl.solve!(rotor, solver),
        capture_history! = (rotor, solver, hist) -> begin
            solver.record_history = true
            pnl.solve!(rotor, solver)
            hist.metric = solver.history.metric
            append!(hist.iter, solver.history.iter)
            append!(hist.t_ns, solver.history.t_ns)
            append!(hist.residual_internal, solver.history.residual_internal)
            hist.t0_ns = solver.history.t0_ns
            solver.record_history = false
        end,
        iterations_of = (solver, hist) -> solver.niter,
        # warmstart across timesteps for the availability run (off during the
        # cold isolated benchmark above; Phase 3 studies this properly)
        body_solvers_for = solver -> (solver.warmstart = true; (solver,)))

    # (g) GMRES + dedicated Barba direct-list ILU(0)
    CTs["krylov_ilu"] = run_config!(io, "krylov_ilu";
        make_solver = rotor -> begin
            P = pnl.ILUPreconditioner(rotor; leaf_size=10,
                multipole_acceptance=1.0)
            st = P.stats
            solver = pnl.KrylovSolver(rotor; method=:gmres,
                itmax=krylov_itmax, atol=krylov_atol, rtol=krylov_rtol,
                backend=backend_sim, preconditioner=P)
            setup = (; t_tree=st["tree_time"],
                t_assembly=st["assembly_time"],
                t_factorize=st["factorization_time"],
                t_precond=st["total_time"])
            notes = "Barba leaf=10 MAC=1.0 nnz=$(st["nnz"]) " *
                    "nnz/N=$(st["nnz_per_panel"]) factor_nnz=$(st["factor_nnz"]) " *
                    "list_time=$(st["interaction_list_time"])s " *
                    "retained=$(st["retained_bytes"])B"
            (solver, setup, notes)
        end,
        solve_once! = (rotor, solver) -> pnl.solve!(rotor, solver),
        capture_history! = (rotor, solver, hist) -> begin
            solver.record_history = true
            pnl.solve!(rotor, solver)
            hist.metric = solver.history.metric
            append!(hist.iter, solver.history.iter)
            append!(hist.t_ns, solver.history.t_ns)
            append!(hist.residual_internal, solver.history.residual_internal)
            hist.t0_ns = solver.history.t0_ns
            solver.record_history = false
        end,
        iterations_of = (solver, hist) -> solver.niter,
        body_solvers_for = solver -> (solver.warmstart = true; (solver,)))
end

################################################################################
# Schema self-validation + CT sanity
################################################################################

nrows = validate_runs_csv(runs_path)
println("\nruns.csv: $nrows schema-valid rows at $runs_path")

println("\nCT[end] sanity (availability proof, NOT converged physics):")
for (config, CT) in sort(collect(CTs))
    println("  $(rpad(config, 22)) CT = $CT")
end

plain = SOLVE_SUMMARY["krylov_gmres"]
ilu = SOLVE_SUMMARY["krylov_ilu"]
saved_per_solve = plain.t_steady - ilu.t_steady
break_even = saved_per_solve > 0 ? ilu.setup_total / saved_per_solve : Inf
println("\nILU versus plain GMRES at matched requested tolerance:")
println("  iterations: $(plain.iterations) -> $(ilu.iterations)")
println("  steady solve: $(plain.t_steady)s -> $(ilu.t_steady)s")
println("  ILU total setup: $(ilu.setup_total)s")
println("  ILU cold first solve: $(ilu.t_cold_first)s")
println("  amortization break-even: $break_even solves")
println("  true residual RMS: $(plain.rms_residual) vs $(ilu.rms_residual)")
