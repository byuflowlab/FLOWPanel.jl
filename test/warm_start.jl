using Pkg
Pkg.activate("/Users/ryan/Dropbox/research/projects/FLOWPanel.jl/")
using FLOWPanel
using FLOWPanel.StaticArrays
using LinearAlgebra
using BSON
using PythonPlot
FLOWPanel_dir = normpath(joinpath(splitdir(pathof(FLOWPanel))[1], ".."))
include("auxiliary_functions.jl")

# set number of threads to 1
BLAS.set_num_threads(1)

# create functions
function prepare_wing(;
        kernel=ConstantSource(),
        AR=10, span=1.0,
        nc=10, ns=30,
        freestream=SVector{3}(1.0,0,0),
    )
    # create wing
    wing = create_wing_structured(kernel; span, AR, nc, ns)

    # apply freestream
    apply_freestream!(wing, freestream)

    return wing
end

function warm_start_wing_gmres_mf(wing, tolerance, max_iterations, expansion_order, leaf_size, multipole_threshold, reuse_tree::Bool, freestream=SVector{3}(1.0,0,0); St, Δt=0.001, f=1.0)

    # get St parameters
    U = norm(freestream)

    # reset wing
    reset_potential_velocity!(wing)
    apply_freestream!(wing, freestream) # apply freestream
    FLOWPanel.set_unit_strength!(wing)

    # create solver
    scheme = FLOWPanel.Scheme{FLOWPanel.DirectNeumann, FlowTangency}
    solver = FLOWPanel.MatrixFreeSolver(wing, scheme; expansion_order, leaf_size, multipole_threshold, reuse_tree)
    print(" Created solver.")

    # solve
    this_tolerance = tolerance
    solve!(wing, solver; verbose=false, tolerance=this_tolerance, max_iterations, rtol=0.0)
    cold_start_strengths = deepcopy(wing.strengths)
    print(" Finished cold start in $(solver.solver.stats.niter) iterations.")

    # warm start
    ts_solve = zeros(length(St))
    resids_max = zeros(length(St))
    resids_mean = zeros(length(St))
    dummy_velocity = zeros(SVector{3,Float64}, length(wing.panels))
    for (i,st) in enumerate(St)
        L = st * U / f
        Δv = SVector{3}(0.0,0.0,2*pi^2 * f * st * U * Δt)
        reset_potential_velocity!(wing)
        apply_freestream!(wing, freestream + Δv)

        # set initial strength
        for i in eachindex(cold_start_strengths)
            wing.strengths[i] = cold_start_strengths[i]
        end
        FLOWPanel.grid_2_panels_strength!(wing)

        # solve
        this_solver = deepcopy(solver)
        ts_solve[i] = @elapsed solve!(wing, this_solver; verbose=false, tolerance=this_tolerance, max_iterations, rtol=0.0, history=false)
        print(" Solved St=$st in $(solver.solver.stats.niter) iterations.")

        # get residual
        resid_max, resid_mean = get_residual(dummy_velocity, wing, freestream + Δv)
        resids_max[i] = resid_max
        resids_mean[i] = resid_mean
        print(" Got resid.")
    end

    # cold start comparison
    ts_solve_cold = zeros(length(St))
    resids_max_cold = zeros(length(St))
    resids_mean_cold = zeros(length(St))
    for (i,st) in enumerate(St)
        L = st * U / f
        Δv = SVector{3}(0.0,0.0,2*pi^2 * f * st * U * Δt)
        reset_potential_velocity!(wing)
        apply_freestream!(wing, freestream + Δv)

        # set unit strengths
        FLOWPanel.set_unit_strength!(wing)

        # solve
        ts_solve_cold[i] = @elapsed solve!(wing, solver; verbose=false, tolerance=this_tolerance, max_iterations, rtol=0.0, history=false)
        print(" Cold solved St=$st in $(solver.solver.stats.niter) iterations.")

        # get residual
        resid_max, resid_mean = get_residual(dummy_velocity, wing, freestream + Δv)
        resids_max_cold[i] = resid_max
        resids_mean_cold[i] = resid_mean
        print(" Got resid.")
    end

    return ts_solve, resids_max, resids_mean, ts_solve_cold, resids_max_cold, resids_mean_cold
end

function guess_next_strength!(panels::FLOWPanel.AbstractPanels, history, dt, i_step)
    h = dt
    h2 = h^2
    for i_panel in eachindex(panels.panels)

        # estimate derivatives
        σ0 = panels.panels[i_panel].strength[1]
        σm1 = history[i_panel,1]
        σm2 = history[i_panel,2]
        σm3 = history[i_panel,3]
        if i_step > 3
            dσdt = (11/6*σ0 - 3*σm1 + 3/2*σm2 - 1/3*σm3) / h
            d2σdt2 = (2*σ0 - 5*σm1 + 4*σm2 - σm3) / h2
        elseif i_step == 3
            dσdt = (3/2*σ0 - 2*σm1 + 1/2*σm2) / h
            d2σdt2 = (σ0 - 2*σm1 + σm2) / h2
        elseif i_step == 2
            dσdt = (σ0 - σm1) / h
            d2σdt2 = 0.0
        elseif i_step == 1
            dσdt = 0.0
            d2σdt2 = 0.0
        end

        # extrapolate
        σ1 = σ0 + dσdt * h + d2σdt2 * h2 / 2

        # store guess in panels.strengths (not panels.panels)
        panels.strengths[i_panel] = SVector{1}(σ1)

    end

    # update strength history
    update_strength_history!(history, panels)

    # update panels with new guess
    grid_2_panels_strength!(panels)
end

function update_strength_history!(history, panels::FLOWPanel.AbstractPanels)
    for i in size(history,2):-1:2
        history[:,i] .= history[:,i-1]
    end
    for i in eachindex(panels.panels)
        history[i,1] = panels.panels[i].strength[1]
    end
end

function warm_start_wing_gmres_mf_cycle(wing, tolerance, max_iterations, expansion_order, leaf_size, multipole_threshold, reuse_tree::Bool, freestream=SVector{3}(1.0,0,0); St, Δt=0.001, f=1.0, nt=10, t0=0.0, smart_guess=false)

    # get St parameters
    U = norm(freestream)

    # other parameters
    scheme = FLOWPanel.Scheme{FLOWPanel.DirectNeumann, FlowTangency}
    this_tolerance = tolerance
    strength_history = fill(NaN, length(wing.panels), 3)

    # warm start
    ts_solve = zeros(nt, length(St))
    resids_max = zeros(nt, length(St))
    resids_mean = zeros(nt, length(St))
    niter = zeros(Int,nt,length(St))
    histories = Matrix{Vector{Float64}}(undef,nt,length(St))
    dummy_velocity = zeros(SVector{3,Float64}, length(wing.panels))
    for (i,st) in enumerate(St)
        L = st * U / f

        # t=0
        v = pi*f*L*cos(2*pi*f*(t0/f+0))
        Δv = SVector{3}(0,0,v)

        # reset wing
        reset_potential_velocity!(wing)
        apply_freestream!(wing, freestream + Δv) # apply freestream
        FLOWPanel.set_unit_strength!(wing)

        # reset strength history
        smart_guess && strength_history .= NaN

        # create solver
        solver = FLOWPanel.MatrixFreeSolver(wing, scheme; expansion_order, leaf_size, multipole_threshold, reuse_tree)
        print(" Created solver.")

        # solve
        solve!(wing, solver; verbose=false, tolerance=this_tolerance, max_iterations, rtol=0.0)
        if smart_guess
            update_strength_history!(strength_history, wing)
        end
        print(" Finished cold start in $(solver.solver.stats.niter) iterations.")

        for i_t in 1:nt
            # velocity
            v = pi*f*L*cos(2*pi*f*(i_t*Δt + t0/f))
            Δv = SVector{3}(0.0,0.0,v)
            reset_potential_velocity!(wing)
            apply_freestream!(wing, freestream + Δv)

            # Taylor series guess
            smart_guess && guess_next_strength!(wing, strength_history, Δt, i_t)

            # solve
            solver = FLOWPanel.MatrixFreeSolver(wing, scheme; expansion_order, leaf_size, multipole_threshold, reuse_tree)
            ts_solve[i_t,i] = @elapsed solve!(wing, solver; verbose=false, tolerance=this_tolerance, max_iterations, rtol=0.0, history=true)
            print(" Solved St=$st in $(solver.solver.stats.niter) iterations.")
            niter[i_t,i] = solver.solver.stats.niter
            histories[i_t,i] = deepcopy(solver.solver.stats.residuals)

            # residual
            resid_max, resid_mean = get_residual(dummy_velocity, wing, freestream + Δv)
            resids_max[i_t,i] = resid_max
            resids_mean[i_t,i] = resid_mean
            print(" Got resid.")
        end
        println()
    end

    print("\n\tMF-GMRES-a Cold Start:")
    # cold start comparison
    ts_solve_cold = zeros(nt, length(St))
    resids_max_cold = zeros(nt, length(St))
    resids_mean_cold = zeros(nt, length(St))
    niter_cold = zeros(Int,nt,length(St))
    histories_cold = Matrix{Vector{Float64}}(undef,nt, length(St))
    for (i,st) in enumerate(St)
        L = st * U / f
        for i_t in 1:nt
            # velocity
            v = pi*f*L*cos(2*pi*f*(i_t*Δt + t0/f))
            Δv = SVector{3}(0.0,0.0,v)
            reset_potential_velocity!(wing)
            apply_freestream!(wing, freestream + Δv)
            FLOWPanel.set_unit_strength!(wing)

            # solve
            this_solver = FLOWPanel.MatrixFreeSolver(wing, scheme; expansion_order, leaf_size, multipole_threshold, reuse_tree)
            print(" Created solver.")
            ts_solve_cold[i_t,i] = @elapsed solve!(wing, this_solver; verbose=false, tolerance=this_tolerance, max_iterations, rtol=0.0, history=true)
            print(" Solved St=$st in $(this_solver.solver.stats.niter) iterations.")
            niter_cold[i_t,i] = this_solver.solver.stats.niter
            histories_cold[i_t,i] = deepcopy(this_solver.solver.stats.residuals)

            # residual
            resid_max, resid_mean = get_residual(dummy_velocity, wing, freestream + Δv)
            resids_max_cold[i_t,i] = resid_max
            resids_mean_cold[i_t,i] = resid_mean
            print(" Got resid.")
        end
    end
    println()

    return ts_solve, resids_max, resids_mean, niter, ts_solve_cold, resids_max_cold, resids_mean_cold, niter_cold, histories, histories_cold
end

function warm_start_wing_fgs(wing, tolerance, max_iterations, expansion_order, leaf_size, multipole_threshold, relaxation, reuse_tree::Bool, freestream=SVector{3}(1.0,0,0); St, Δt=0.001, f=1.0)

    # get St parameters
    U = norm(freestream)

    # reset wing
    reset_potential_velocity!(wing)
    apply_freestream!(wing, freestream) # apply freestream
    FLOWPanel.set_unit_strength!(wing)

    # create solver
    scheme = FLOWPanel.Scheme{FLOWPanel.DirectNeumann, FlowTangency}
    solver = FLOWPanel.FastGaussSeidel(wing, scheme; expansion_order, leaf_size, multipole_threshold, reuse_tree)
    print(" Created solver.")

    # solve
    this_tolerance = tolerance
    solve!(wing, solver; verbose=false, tolerance=this_tolerance, max_iterations, relaxation, save_history=true)
    cold_start_strengths = deepcopy(wing.strengths)
    print(" Finished cold start in $(length(solver.convergence_history)) iterations.")

    # warm start
    resids_max = zeros(length(St))
    resids_mean = zeros(length(St))
    dummy_velocity = zeros(SVector{3,Float64}, length(wing.panels))
    ts_solve = zeros(length(St))
    for (i,st) in enumerate(St)
        L = st * U / f
        Δv = SVector{3}(0.0,0.0,2*pi^2 * f * st * U * Δt)
        reset_potential_velocity!(wing)
        apply_freestream!(wing, freestream)
        apply_freestream!(wing, Δv)

        # set initial strength
        for i in eachindex(cold_start_strengths)
            wing.strengths[i] = cold_start_strengths[i]
        end
        grid_2_panels_strength!(wing)

        # solve
        ts_solve[i] = @elapsed solve!(wing, solver; verbose=false, tolerance=this_tolerance, max_iterations, relaxation, save_history=true)
        print(" Solved St=$st in $(length(solver.convergence_history)) iterations.")

        # get residual
        resid_max, resid_mean = get_residual(dummy_velocity, wing, freestream + Δv)
        resids_max[i] = resid_max
        resids_mean[i] = resid_mean
        print(" Got resid.")
    end

    return ts_solve, resids_max, resids_mean
end

function warm_start_wing_fgs_cycle(wing, tolerance, max_iterations, expansion_order, leaf_size, multipole_threshold, relaxation, reuse_tree::Bool, freestream=SVector{3}(1.0,0,0); St, Δt=0.001, f=1.0, nt=10, t0=0.0, smart_guess=false)

    # get St parameters
    U = norm(freestream)

    # get other parameters
    this_tolerance = tolerance
    scheme = FLOWPanel.Scheme{FLOWPanel.DirectNeumann, FlowTangency}

    # warm start
    strength_history = fill(NaN, length(wing.panels), 3)
    ts_solve = zeros(nt, length(St))
    resids_max = zeros(nt, length(St))
    resids_mean = zeros(nt, length(St))
    niter = zeros(Int, nt, length(St))
    histories = Matrix{Vector{Float64}}(undef,nt, length(St))
    dummy_velocity = zeros(SVector{3,Float64}, length(wing.panels))
    for (i,st) in enumerate(St)
        L = st * U / f

        # t=0
        v = pi*f*L*cos(2*pi*f*(t0/f))
        Δv = SVector{3}(0,0,v)

        # reset wing
        reset_potential_velocity!(wing)
        apply_freestream!(wing, freestream + Δv) # apply freestream
        FLOWPanel.set_unit_strength!(wing)

        # create solver
        solver = FLOWPanel.FastGaussSeidel(wing, scheme; expansion_order, leaf_size, multipole_threshold, reuse_tree)
        print(" Created solver.")

        # solve
        solve!(wing, solver; verbose=false, tolerance=this_tolerance, max_iterations, relaxation, save_history=false)
        print(" Finished cold start.")

        for i_t in 1:nt
            # velocity
            v = pi*f*L*cos(2*pi*f*(i_t*Δt + t0/f))
            Δv = SVector{3}(0.0,0.0,v)
            reset_potential_velocity!(wing)
            apply_freestream!(wing, freestream + Δv)

            # Taylor series guess
            println("before:")
            @show wing.panels[1].strength[1]
            smart_guess && guess_next_strength!(wing, strength_history, Δt, i_t)
            println("smart guess:")
            @show wing.panels[1].strength[1]

            # solve
            ts_solve[i_t,i] = @elapsed solve!(wing, solver; verbose=false, tolerance=this_tolerance, max_iterations, relaxation, save_history=false)
            println("Solved.")
            @show wing.panels[1].strength[1]
            print(" Solved St=$st.")

            # get residual
            resid_max, resid_mean = get_residual(dummy_velocity, wing, freestream + Δv)
            resids_max[i_t,i] = resid_max
            resids_mean[i_t,i] = resid_mean
            print(" Got resid.")

        end
    end

    # get history with L2 norm
    for (i,st) in enumerate(St)
        L = st * U / f

        # t=0
        v = pi*f*L*cos(2*pi*f*(t0/f))
        Δv = SVector{3}(0,0,v)

        # reset wing
        reset_potential_velocity!(wing)
        apply_freestream!(wing, freestream + Δv) # apply freestream
        FLOWPanel.set_unit_strength!(wing)

        # reset strength history
        strength_history .= NaN

        # create solver
        solver = FLOWPanel.FastGaussSeidel(wing, scheme; expansion_order, leaf_size, multipole_threshold, reuse_tree)
        print(" Created history solver.")

        # solve
        solve!(wing, solver; verbose=false, tolerance=this_tolerance, max_iterations, relaxation, save_history=true)
        print(" Finished history cold start in $(length(solver.convergence_history)) iterations.")

        for i_t in 1:nt
            # velocity
            v = pi*f*L*cos(2*pi*f*(i_t*Δt + t0/f))
            Δv = SVector{3}(0.0,0.0,v)
            reset_potential_velocity!(wing)
            apply_freestream!(wing, freestream + Δv)

            # Taylor series guess
            smart_guess && guess_next_strength!(wing, strength_history, Δt, i_t)

            # solve
            solve!(wing, solver; verbose=false, tolerance, max_iterations=100, relaxation, save_history=true)
            niter[i_t,i] = length(solver.convergence_history)
            histories[i_t,i] = deepcopy(solver.convergence_history)
            print(" Finished history warm start in $(length(solver.convergence_history)) iterations.")

        end
    end


    print("\n\tFGS-a Cold Start:")
    # cold start comparison
    ts_solve_cold = zeros(nt, length(St))
    resids_max_cold = zeros(nt, length(St))
    resids_mean_cold = zeros(nt, length(St))
    niter_cold = zeros(Int, nt, length(St))
    histories_cold = Matrix{Vector{Float64}}(undef, nt, length(St))
    for (i,st) in enumerate(St)
        L = st * U / f
        for i_t in 1:nt
            # velocity
            v = pi*f*L*cos(2*pi*f*(i_t*Δt + t0/f))
            Δv = SVector{3}(0.0,0.0,v)
            reset_potential_velocity!(wing)
            apply_freestream!(wing, freestream + Δv)
            FLOWPanel.set_unit_strength!(wing)

            # solve
            solver = FLOWPanel.FastGaussSeidel(wing, scheme; expansion_order, leaf_size, multipole_threshold, reuse_tree)
            print(" Created solver.")
            ts_solve_cold[i_t,i] = @elapsed solve!(wing, solver; verbose=false, tolerance=this_tolerance, max_iterations, relaxation, save_history=false)
            print(" Solved St=$st.")

            # residual
            resid_max, resid_mean = get_residual(dummy_velocity, wing, freestream + Δv)
            resids_max_cold[i_t,i] = resid_max
            resids_mean_cold[i_t,i] = resid_mean
            print(" Got resid.")
        end
    end

    for (i,st) in enumerate(St)
        L = st * U / f
        for i_t in 1:nt
            # velocity
            v = pi*f*L*cos(2*pi*f*(i_t*Δt + t0/f))
            Δv = SVector{3}(0.0,0.0,v)
            reset_potential_velocity!(wing)
            apply_freestream!(wing, freestream + Δv)
            FLOWPanel.set_unit_strength!(wing)

            # solve
            solver = FLOWPanel.FastGaussSeidel(wing, scheme; expansion_order, leaf_size, multipole_threshold, reuse_tree)
            print(" Created history solver.")
            solve!(wing, solver; verbose=false, tolerance=this_tolerance, max_iterations=100, relaxation, save_history=true)
            print(" Solved history St=$st in $(length(solver.convergence_history)) iterations.")
            niter_cold[i_t,i] = length(solver.convergence_history)
            histories_cold[i_t,i] = deepcopy(solver.convergence_history)

        end
    end
    println()

    return ts_solve, resids_max, resids_mean, niter, ts_solve_cold, resids_max_cold, resids_mean_cold, niter_cold, histories, histories_cold
end

function warm_start_wing(wing::FLOWPanel.AbstractPanels, tolerance, expansion_order_a, leaf_size_a, multipole_threshold_a, expansion_order_b, leaf_size_b, multipole_threshold_b, relaxation, nc, ns;
        freestream=SVector{3}(1.0,0,0),
        benchmark_dir,
        St, Δt, f
    )

    # benchmark MF-GMRES-a
    print("\n\tMF-GMRES-a...")
    reuse_tree = false
    max_iterations = 30
    ts_gmres_mf_a, resids_max_gmres_mf_a, resids_mean_gmres_mf_a = warm_start_wing_gmres_mf(wing, tolerance, max_iterations, expansion_order_a, leaf_size_a, multipole_threshold_a, reuse_tree, freestream; St, Δt, f)

    # benchmark FGS-a
    print("\n\tFast Gauss-Seidel-a...")
    max_iterations=500
    reuse_tree = false
    ts_fgs_a, resids_max_fgs_a, resids_mean_fgs_a = warm_start_wing_fgs(wing, tolerance, max_iterations, expansion_order_a, leaf_size_a, multipole_threshold_a, relaxation, reuse_tree, freestream; St, Δt, f)

    return (ts_gmres_mf_a, resids_max_gmres_mf_a, resids_mean_gmres_mf_a), (ts_fgs_a, resids_max_fgs_a, resids_mean_fgs_a)
end

function warm_start_wing_cycle(wing::FLOWPanel.AbstractPanels, tolerance, expansion_order_a, leaf_size_a, multipole_threshold_a, expansion_order_b, leaf_size_b, multipole_threshold_b, relaxation, nc, ns;
        freestream=SVector{3}(1.0,0,0),
        benchmark_dir,
        St, Δt, f, nt, t0, smart_guess
    )

    # benchmark MF-GMRES-a
    print("\n\tMF-GMRES-a...")
    reuse_tree = false
    max_iterations = 30
    ts_gmres_mf_a, resids_max_gmres_mf_a, resids_mean_gmres_mf_a, niter_gmres_mf_a, ts_gmres_mf_a_cold, resids_max_gmres_mf_a_cold, resids_mean_gmres_mf_a_cold, niter_gmres_mf_a_cold, histories_gmres_mf_a, histories_gmres_mf_a_cold = warm_start_wing_gmres_mf_cycle(wing, tolerance, max_iterations, expansion_order_a, leaf_size_a, multipole_threshold_a, reuse_tree, freestream; St, Δt, f, nt, t0, smart_guess)

    # benchmark FGS-a
    print("\n\tFast Gauss-Seidel-a...")
    max_iterations=100
    reuse_tree = false
    ts_fgs_a, resids_max_fgs_a, resids_mean_fgs_a, niter_fgs_a, ts_fgs_a_cold, resids_max_fgs_a_cold, resids_mean_fgs_a_cold, niter_fgs_a_cold, histories_fgs_a, histories_fgs_a_cold = warm_start_wing_fgs_cycle(wing, tolerance, max_iterations, expansion_order_a, leaf_size_a, multipole_threshold_a, relaxation, reuse_tree, freestream; St, Δt, f, nt, t0, smart_guess)

    return (ts_gmres_mf_a, resids_max_gmres_mf_a, resids_mean_gmres_mf_a, niter_gmres_mf_a), (ts_gmres_mf_a_cold, resids_max_gmres_mf_a_cold, resids_mean_gmres_mf_a_cold, niter_gmres_mf_a_cold), (ts_fgs_a, resids_max_fgs_a, resids_mean_fgs_a, niter_fgs_a), (ts_fgs_a_cold, resids_max_fgs_a_cold, resids_mean_fgs_a_cold, niter_fgs_a_cold), (histories_gmres_mf_a, histories_gmres_mf_a_cold), (histories_fgs_a, histories_fgs_a_cold)
end

function warm_start_wing(ms, benchmark_dir;
        St, Δt=0.001, f=1.0,
        replace_file=false,
        freestream=SVector{3}(1.0,0,0),
        kernel=ConstantSource(),
        nc0=5, ns0=50, span=1.0, AR=10,
        tolerances = [1e-3, 1e-6, 1e-9],
        ps_a = fill((5,7,12), length(ms)),
        ls_a = fill((20,30,100), length(ms)),
        mts_a = fill((0.65,0.4,0.4), length(ms)),
        ps_b = [(1,1,1), (3,1,1), (4,1,1), (4,10,6), (4,8,14), (4,8,14), (4,8,14), (4,8,14)],
        ls_b = [(200,200,200), (100,800,800), (100,800,800), (100,400,800), (100,300,400), (100,300,400), (100,300,400), (100,300,400)],
        mts_b = [(0.4,0.4,0.4), (0.55,0.6,0.6), (0.5,0.6,0.6), (0.5,0.6,0.4), (0.5,0.4,0.45), (0.5,0.4,0.45), (0.5,0.4,0.45), (0.5,0.4,0.45)],
        relaxation = 1.4,
    )
    # check if folders exist
    !isdir(joinpath(FLOWPanel_dir, benchmark_dir)) && mkdir(joinpath(FLOWPanel_dir, benchmark_dir))
    for dir in ("warm_start_gmres_mf_a", "benchmark_fgs_a")
        if !isdir(joinpath(FLOWPanel_dir, benchmark_dir, dir))
            mkdir(joinpath(FLOWPanel_dir, benchmark_dir, dir))
        end
    end

    # check if file exists
    file = joinpath(FLOWPanel_dir, benchmark_dir, "assembled_benchmarks.bson")
    if isfile(file) && !replace_file
        BSON.@load file ws_gmres_mf_a ws_fgs_a n_panels tolerances ps_a ls_a mts_a ps_b ls_b mts_b nc0 ns0 span AR relaxation St Δt f
    else
        # preallocate benchmarks
        ws_gmres_mf_a = zeros(3, length(St), length(tolerances), length(ms))
        ws_fgs_a = zeros(3, length(St), length(tolerances), length(ms))
        n_panels = fill(-1, length(ms))

        # populate with NaNs
        ws_gmres_mf_a .= NaN
        ws_fgs_a .= NaN
    end

    # run benchmarks
    for (i_m, (m, p_a, l_a, mt_a, p_b, l_b, mt_b)) in enumerate(zip(ms, ps_a, ls_a, mts_a, ps_b, ls_b, mts_b))
        # create wing
        nc, ns = nc0*m, ns0*m
        println("\n#------- m = $m; n_panels = $(ns0*nc0*2*m^2) -------#")
        print("\n\tpreparing warm start wing...")
        wing = prepare_wing(; kernel, span, AR, nc, ns, freestream)
        print(" Done.\n")
        n_panels[i_m] = length(wing.panels)

        # loop over tolerances
        for (i_tol,tol) in enumerate(tolerances)

            println("\n#--- m = $m; tolerance: $tol ---#\n")

            # run benchmarks
            gmres_mf_a, fgs_a = warm_start_wing(wing, tol, p_a[i_tol], l_a[i_tol], mt_a[i_tol], p_b[i_tol], l_b[i_tol], mt_b[i_tol], relaxation, nc, ns; freestream, benchmark_dir, St, f, Δt)

            # store results
            for i in 1:3
                ws_gmres_mf_a[i,:,i_tol,i_m] .= gmres_mf_a[i]
                ws_fgs_a[i,:,i_tol,i_m] .= fgs_a[i]
            end

            # store bson
            BSON.@save file ws_gmres_mf_a ws_fgs_a n_panels tolerances ps_a ls_a mts_a ps_b ls_b mts_b nc0 ns0 span AR relaxation St Δt f
        end
    end

    return ws_gmres_mf_a, ws_fgs_a
end

function warm_start_wing_cycle(ms, benchmark_dir;
        St, Δt=0.001, f=1.0, nt, t0, smart_guess,
        replace_file=false,
        freestream=SVector{3}(1.0,0,0),
        kernel=ConstantSource(),
        nc0=5, ns0=50, span=1.0, AR=10,
        tolerances = [1e-3, 1e-6, 1e-9],
        ps_a = fill((5,7,12), length(ms)),
        ls_a = fill((20,30,100), length(ms)),
        mts_a = fill((0.65,0.4,0.4), length(ms)),
        ps_b = [(1,1,1), (3,1,1), (4,1,1), (4,10,6), (4,8,14), (4,8,14), (4,8,14), (4,8,14)],
        ls_b = [(200,200,200), (100,800,800), (100,800,800), (100,400,800), (100,300,400), (100,300,400), (100,300,400), (100,300,400)],
        mts_b = [(0.4,0.4,0.4), (0.55,0.6,0.6), (0.5,0.6,0.6), (0.5,0.6,0.4), (0.5,0.4,0.45), (0.5,0.4,0.45), (0.5,0.4,0.45), (0.5,0.4,0.45)],
        relaxation = 1.4,
    )
    # check if folders exist
    !isdir(joinpath(FLOWPanel_dir, benchmark_dir)) && mkdir(joinpath(FLOWPanel_dir, benchmark_dir))
    for dir in ("warm_start_gmres_mf_a", "benchmark_fgs_a")
        if !isdir(joinpath(FLOWPanel_dir, benchmark_dir, dir))
            mkdir(joinpath(FLOWPanel_dir, benchmark_dir, dir))
        end
    end

    # check if file exists
    file = joinpath(FLOWPanel_dir, benchmark_dir, "assembled_benchmarks_cycle.bson")
    if isfile(file) && !replace_file
        BSON.@load file ws_gmres_mf_a ws_gmres_mf_a_cold ws_fgs_a ws_fgs_a_cold histories_gmres_mf_a histories_gmres_mf_a_cold histories_fgs_a histories_fgs_a_cold n_panels tolerances ps_a ls_a mts_a ps_b ls_b mts_b nc0 ns0 span AR relaxation St Δt f smart_guess
    else
        # preallocate benchmarks
        ws_gmres_mf_a = zeros(nt, 4, length(St), length(tolerances), length(ms))
        ws_gmres_mf_a_cold = zeros(nt, 4, length(St), length(tolerances), length(ms))
        ws_fgs_a = zeros(nt, 4, length(St), length(tolerances), length(ms))
        ws_fgs_a_cold = zeros(nt, 4, length(St), length(tolerances), length(ms))
        histories_gmres_mf_a = Matrix{Matrix{Vector{Float64}}}(undef, length(tolerances), length(ms))
        for i in eachindex(histories_gmres_mf_a)
            histories_gmres_mf_a[i] = [[NaN] for _ in 1:nt, _ in 1:length(St)]
        end
        histories_gmres_mf_a_cold = deepcopy(histories_gmres_mf_a)
        histories_fgs_a = deepcopy(histories_gmres_mf_a)
        histories_fgs_a_cold = deepcopy(histories_gmres_mf_a)
        n_panels = fill(-1, length(ms))

        # populate with NaNs
        ws_gmres_mf_a .= NaN
        ws_fgs_a .= NaN
    end

    # run benchmarks
    for (i_m, (m, p_a, l_a, mt_a, p_b, l_b, mt_b)) in enumerate(zip(ms, ps_a, ls_a, mts_a, ps_b, ls_b, mts_b))
        # create wing
        nc, ns = nc0*m, ns0*m
        println("\n#------- m = $m; n_panels = $(ns0*nc0*2*m^2) -------#")
        print("\n\tpreparing warm start wing cycle...")
        wing = prepare_wing(; kernel, span, AR, nc, ns, freestream)
        print(" Done.\n")
        n_panels[i_m] = length(wing.panels)

        # loop over tolerances
        for (i_tol,tol) in enumerate(tolerances)

            println("\n#--- m = $m; tolerance: $tol ---#\n")

            # run benchmarks
            gmres_mf_a, gmres_mf_a_cold, fgs_a, fgs_a_cold, histories_gmres, histories_fgs = warm_start_wing_cycle(wing, tol, p_a[i_tol], l_a[i_tol], mt_a[i_tol], p_b[i_tol], l_b[i_tol], mt_b[i_tol], relaxation, nc, ns; freestream, benchmark_dir, St, f, Δt, nt, t0, smart_guess)

            # store results
            for i in 1:4
                ws_gmres_mf_a[:,i,:,i_tol,i_m] .= gmres_mf_a[i]
                ws_fgs_a[:,i,:,i_tol,i_m] .= fgs_a[i]
                ws_gmres_mf_a_cold[:,i,:,i_tol,i_m] .= gmres_mf_a_cold[i]
                ws_fgs_a_cold[:,i,:,i_tol,i_m] .= fgs_a_cold[i]
                histories_gmres_mf_a[i_tol,i_m] = histories_gmres[1]
                histories_gmres_mf_a_cold[i_tol,i_m] = histories_gmres[2]
                histories_fgs_a[i_tol,i_m] = histories_fgs[1]
                histories_fgs_a_cold[i_tol,i_m] = histories_fgs[2]
            end

            # store bson
            BSON.@save file ws_gmres_mf_a ws_gmres_mf_a_cold ws_fgs_a ws_fgs_a_cold histories_gmres_mf_a histories_gmres_mf_a_cold histories_fgs_a histories_fgs_a_cold n_panels tolerances ps_a ls_a mts_a ps_b ls_b mts_b nc0 ns0 span AR relaxation St Δt f ms smart_guess
        end
    end

    return ws_gmres_mf_a, ws_fgs_a, ws_gmres_mf_a_cold, ws_fgs_a_cold, histories_gmres_mf_a, histories_gmres_mf_a_cold, histories_fgs_a, histories_fgs_a_cold
end

function plot_warm_start_cycle_history(file, save_name; i_tol=2, i_m=1, i_t=2)

    # load data
    BSON.@load file ws_gmres_mf_a ws_gmres_mf_a_cold ws_fgs_a ws_fgs_a_cold histories_gmres_mf_a histories_gmres_mf_a_cold histories_fgs_a histories_fgs_a_cold n_panels tolerances ps_a ls_a mts_a ps_b ls_b mts_b nc0 ns0 span AR relaxation St Δt f

    # prepare colors
    n_colors = 6
    cmap = PythonPlot.get_cmap("viridis_r", n_colors)

    # generate figure
    fig = figure("history")
    fig.clf()
    fig.add_subplot(421, ylabel="L-2 norm residual")
    fig.add_subplot(423, ylabel="L-2 norm residual")
    fig.add_subplot(425, ylabel="L-2 norm residual")
    fig.add_subplot(427, xlabel="iteration", ylabel="L-2 norm residual")
    fig.add_subplot(422)
    fig.add_subplot(424)
    fig.add_subplot(426)
    fig.add_subplot(428, xlabel="iteration")
    axs = fig.get_axes()

    # loop over St
    for (i_st, st) in enumerate(St)
        # plot GMRES

        ## get time
        t = ws_gmres_mf_a[i_t,1,i_st,i_tol,i_m]
        t_cold = ws_gmres_mf_a_cold[i_t,1,i_st,i_tol,i_m]

        ## get histories
        history = histories_gmres_mf_a[i_tol,i_m][i_t,i_st]
        history_cold = histories_gmres_mf_a_cold[i_tol,i_m][i_t,i_st]

        i_ax = i_st-1
        axs[i_ax].plot(1:length(history), history, color=cmap(2))
        axs[i_ax].plot(1:length(history_cold), history_cold, color=cmap(2), "--")

        # plot FGS

        ## get time
        t = ws_fgs_a[i_t,1,i_st,i_tol,i_m]
        t_cold = ws_fgs_a_cold[i_t,1,i_st,i_tol,i_m]

        ## get histories
        history = histories_fgs_a[i_tol,i_m][i_t,i_st]
        history_cold = histories_fgs_a_cold[i_tol,i_m][i_t,i_st]

        i_ax += 4
        axs[i_ax].plot(1:length(history), history, color=cmap(2))
        axs[i_ax].plot(1:length(history_cold), history_cold, color=cmap(2), "--")

    end
    axs[7].legend(["warm started", "cold started"])

    for i in 0:7
        # axs[i].set_ylim([1,24])
        # axs[i].set_yticks([5,10,15,20])
        # axs[i].set_xscale("log")
        axs[i].set_xlim([1,25])
        axs[i].set_yscale("log")
        axs[i].set_ylim([1e-7, 1e1])
    end
    tight_layout()

    savefig(save_name*".png")
end

function plot_warm_start_cycle_history_individual(file, save_name; i_tol=2, i_m=1, i_t=2)

    # load data
    BSON.@load file ws_gmres_mf_a ws_gmres_mf_a_cold ws_fgs_a ws_fgs_a_cold histories_gmres_mf_a histories_gmres_mf_a_cold histories_fgs_a histories_fgs_a_cold n_panels tolerances ps_a ls_a mts_a ps_b ls_b mts_b nc0 ns0 span AR relaxation St Δt f

    # prepare colors
    n_colors = 6
    cmap = PythonPlot.get_cmap("viridis_r", n_colors)

    # loop over St
    for (i_st, st) in enumerate(St)
        # plot GMRES

        ## generate figure
        fig = figure("history_individual")
        fig.clf()
        fig.add_subplot(111, ylabel="L-2 norm residual")
        ax = fig.get_axes()[0]
        i_st == 4 && ax.set_xlabel("iteration")

        ## get time
        t = ws_gmres_mf_a[i_t,1,i_st,i_tol,i_m]
        t_cold = ws_gmres_mf_a_cold[i_t,1,i_st,i_tol,i_m]

        ## get histories
        history = histories_gmres_mf_a[i_tol,i_m][i_t,i_st]
        history_cold = histories_gmres_mf_a_cold[i_tol,i_m][i_t,i_st]

        ax.plot(1:length(history), history, color=cmap(2))
        ax.plot(1:length(history_cold), history_cold, color=cmap(2), "--")
        ax.set_xlim([1,25])
        ax.set_yscale("log")
        ax.set_ylim([1e-6, 1e1])
        tight_layout()
        savefig(save_name*"_gmres_mf_a_st$(i_st).png")

        # plot FGS

        ## generate figure
        fig = figure("history_individual")
        fig.clf()
        fig.add_subplot(111)
        ax = fig.get_axes()[0]
        i_st == 4 && ax.set_xlabel("iteration")

        ## get time
        t = ws_fgs_a[i_t,1,i_st,i_tol,i_m]
        t_cold = ws_fgs_a_cold[i_t,1,i_st,i_tol,i_m]

        ## get histories
        history = histories_fgs_a[i_tol,i_m][i_t,i_st]
        history_cold = histories_fgs_a_cold[i_tol,i_m][i_t,i_st]

        ax.plot(1:length(history), history, color=cmap(2))
        ax.plot(1:length(history_cold), history_cold, color=cmap(2), "--")
        ax.set_xlim([1,22])
        ax.set_yscale("log")
        ax.set_ylim([1e-7, 1e0])
        i_st == 4 && ax.legend(["warm started", "cold started"])
        tight_layout()
        savefig(save_name*"_fgs_a_st$(i_st).png")

    end
end

function plot_warm_start_cycle(file, save_name; i_tol=2, i_m=1)

    # load data
    BSON.@load file ws_gmres_mf_a ws_gmres_mf_a_cold ws_fgs_a ws_fgs_a_cold n_panels tolerances ps_a ls_a mts_a ps_b ls_b mts_b nc0 ns0 span AR relaxation St Δt f ms

    # prepare colors
    n_colors = 6
    cmap = PythonPlot.get_cmap("viridis_r", n_colors)

    # generate figure
    fig = figure("iterations per timestep")
    fig.clf()
    fig.add_subplot(421, ylabel="iterations")
    fig.add_subplot(423, ylabel="iterations")
    fig.add_subplot(425, ylabel="iterations")
    fig.add_subplot(427, xlabel="time step", ylabel="iterations")
    fig.add_subplot(422)
    fig.add_subplot(424)
    fig.add_subplot(426)
    fig.add_subplot(428, xlabel="time step")
    axs = fig.get_axes()

    # loop over St
    for (i_st, st) in enumerate(St)
        # plot GMRES
        n_iter = ws_gmres_mf_a[:,4,i_st,i_tol,i_m]
        n_iter_cold = ws_gmres_mf_a_cold[:,4,i_st,i_tol,i_m]
        i_ax = i_st-1
        axs[i_ax].plot(1:length(n_iter), n_iter, color=cmap(2))
        axs[i_ax].plot(1:length(n_iter_cold), n_iter_cold, color=cmap(2), "--")

        # plot FGS
        n_iter = ws_fgs_a[:,4,i_st,i_tol,i_m]
        n_iter_cold = ws_fgs_a_cold[:,4,i_st,i_tol,i_m]
        i_ax += 4
        axs[i_ax].plot(1:length(n_iter), n_iter, color=cmap(2))
        axs[i_ax].plot(1:length(n_iter_cold), n_iter_cold, color=cmap(2), "--")
    end
    axs[7].legend(["warm started", "cold started"])

    for i in 0:7
        axs[i].set_ylim([0,32])
        #axs[i].set_yticks([5,10,15,20,25])
    end
    tight_layout()

    savefig(save_name*".png")
end

function plot_warm_start_cycle_individual(file, save_name; i_tol=2, i_m=1)

    # load data
    BSON.@load file ws_gmres_mf_a ws_gmres_mf_a_cold ws_fgs_a ws_fgs_a_cold n_panels tolerances ps_a ls_a mts_a ps_b ls_b mts_b nc0 ns0 span AR relaxation St Δt f ms

    # prepare colors
    n_colors = 6
    cmap = PythonPlot.get_cmap("viridis_r", n_colors)

    # generate figure
    fig = figure("iterations per timestep individual")

    # loop over St
    for (i_st, st) in enumerate(St)
        # plot GMRES
        fig.clf()
        fig.add_subplot(111,ylabel="iterations")
        ax=fig.get_axes()[0]
        i_st == 4 && ax.set_xlabel("time step")
        n_iter = ws_gmres_mf_a[:,4,i_st,i_tol,i_m]
        n_iter_cold = ws_gmres_mf_a_cold[:,4,i_st,i_tol,i_m]
        ax.plot(1:length(n_iter), n_iter, color=cmap(2))
        ax.plot(1:length(n_iter_cold), n_iter_cold, color=cmap(2), "--")
        ax.set_ylim([0,32])
        ax.set_yticks([5,10,15,20])
        tight_layout()
        savefig(save_name*"_gmres_mf_a_st$(i_st).png")

        # plot FGS
        fig.clf()
        fig.add_subplot(111)
        ax=fig.get_axes()[0]
        i_st == 4 && ax.set_xlabel("time step")
        n_iter = ws_fgs_a[:,4,i_st,i_tol,i_m]
        n_iter_cold = ws_fgs_a_cold[:,4,i_st,i_tol,i_m]
        ax.plot(1:length(n_iter), n_iter, color=cmap(2))
        ax.plot(1:length(n_iter_cold), n_iter_cold, color=cmap(2), "--")
        ax.set_ylim([0,32])
        ax.set_yticks([5,10,15,20])
        i_st == 4 && ax.legend(["warm started", "cold started"])
        tight_layout()
        savefig(save_name*"_fgs_a_st$(i_st).png")
    end

end

#--- final warm start plots ---#

function plot_warm_start_cycle_history_final(file, smart_file, save_name; i_tol, i_m, i_t, i_tol_smart, i_m_smart, i_t_smart)

    # load data
    BSON.@load file ws_gmres_mf_a ws_gmres_mf_a_cold ws_fgs_a ws_fgs_a_cold histories_gmres_mf_a histories_gmres_mf_a_cold histories_fgs_a histories_fgs_a_cold n_panels tolerances ps_a ls_a mts_a ps_b ls_b mts_b nc0 ns0 span AR relaxation St Δt f

    # prepare colors
    n_colors = 6
    cmap = PythonPlot.get_cmap("viridis_r", n_colors)

    # generate figure
    fig = figure("history_final")
    fig.clf()
    fig.add_subplot(141, xlabel="iteration", ylabel="L-2 norm residual")
    fig.add_subplot(142, xlabel="iteration")#, ylabel="L-2 norm residual")
    fig.add_subplot(143, xlabel="iteration")#, ylabel="L-2 norm residual")
    fig.add_subplot(144, xlabel="iteration")#, ylabel="L-2 norm residual")
    axs = fig.get_axes()

    # loop over St
    maxiter = 0
    for i in eachindex(St)
        maxiter = max(maxiter, length(histories_gmres_mf_a_cold[i_tol,i_m][i_t,i]))
    end

    for (i_st, st) in enumerate(St)
        # plot GMRES

        ## get time
        t_cold = ws_gmres_mf_a_cold[i_t,1,i_st,i_tol,i_m]
        t_cold_fgs = ws_fgs_a_cold[i_t,1,i_st,i_tol,i_m]

        ## get histories
        history_cold = histories_gmres_mf_a_cold[i_tol,i_m][i_t,i_st]
        history_cold_fgs = histories_fgs_a_cold[i_tol,i_m][i_t,i_st]

        i_ax = i_st-1
        axs[i_ax].plot(1:length(history_cold), history_cold, color=cmap(2))
        axs[i_ax].plot(1:length(history_cold_fgs), history_cold_fgs, color=cmap(4))

    end
    axs[3].legend(["MF-GMRES", "FGS"])

    for i in 0:3
        # axs[i].set_ylim([1,24])
        # axs[i].set_yticks([5,10,15,20])
        # axs[i].set_xscale("log")
        axs[i].set_xlim([1,maxiter])
        axs[i].set_yscale("log")
        axs[i].set_ylim([1e-7, 1e1])
    end
    tight_layout()

    savefig(save_name*"_cold.png")

    #--- warm start ---#
    fig = figure("history_warm")
    fig.clf()
    fig.add_subplot(141, xlabel="iteration", ylabel="L-2 norm residual")
    fig.add_subplot(142, xlabel="iteration")#, ylabel="L-2 norm residual")
    fig.add_subplot(143, xlabel="iteration")#, ylabel="L-2 norm residual")
    fig.add_subplot(144, xlabel="iteration")#, ylabel="L-2 norm residual")
    axs = fig.get_axes()

    # loop over St
    maxiter = 0
    for i in eachindex(St)
        maxiter = max(maxiter, length(histories_gmres_mf_a[i_tol,i_m][i_t,i]))
    end
    for (i_st, st) in enumerate(St)
        # plot GMRES

        ## get time
        t = ws_gmres_mf_a[i_t,1,i_st,i_tol,i_m]
        t_fgs = ws_fgs_a[i_t,1,i_st,i_tol,i_m]

        ## get histories
        history = histories_gmres_mf_a[i_tol,i_m][i_t,i_st]
        history_fgs = histories_fgs_a[i_tol,i_m][i_t,i_st]

        i_ax = i_st-1
        axs[i_ax].plot(1:length(history), history, color=cmap(2))
        axs[i_ax].plot(1:length(history_fgs), history_fgs, color=cmap(4))

    end
    axs[3].legend(["MF-GMRES", "FGS"])

    for i in 0:3
        # axs[i].set_ylim([1,24])
        # axs[i].set_yticks([5,10,15,20])
        # axs[i].set_xscale("log")
        axs[i].set_xlim([1,maxiter])
        axs[i].set_yscale("log")
        axs[i].set_ylim([1e-7, 1e1])
    end
    tight_layout()

    savefig(save_name*"_warm.png")

    #--- smart warm start ---#

    # load data
    BSON.@load smart_file ws_gmres_mf_a ws_gmres_mf_a_cold ws_fgs_a ws_fgs_a_cold histories_gmres_mf_a histories_gmres_mf_a_cold histories_fgs_a histories_fgs_a_cold n_panels tolerances ps_a ls_a mts_a ps_b ls_b mts_b nc0 ns0 span AR relaxation St Δt f

    @show tolerances[i_tol_smart], size(ws_gmres_mf_a_cold,5)
    fig = figure("history_smart")
    fig.clf()
    fig.add_subplot(141, xlabel="iteration", ylabel="L-2 norm residual")
    fig.add_subplot(142, xlabel="iteration")#, ylabel="L-2 norm residual")
    fig.add_subplot(143, xlabel="iteration")#, ylabel="L-2 norm residual")
    fig.add_subplot(144, xlabel="iteration")#, ylabel="L-2 norm residual")
    axs = fig.get_axes()

    # loop over St
    maxiter = 0
    for i in eachindex(St)
        maxiter = max(maxiter, length(histories_gmres_mf_a[i_tol_smart,i_m_smart][i_t_smart,i]))
    end
    for (i_st, st) in enumerate(St)
        # plot GMRES

        ## get time
        t = ws_gmres_mf_a[i_t,1,i_st,i_tol_smart,i_m_smart]
        t_fgs = ws_fgs_a[i_t,1,i_st,i_tol_smart,i_m_smart]

        ## get histories
        history = histories_gmres_mf_a[i_tol_smart,i_m_smart][i_t_smart,i_st]
        history_fgs = histories_fgs_a[i_tol_smart,i_m_smart][i_t_smart,i_st]

        i_ax = i_st-1
        axs[i_ax].plot(1:length(history), history, color=cmap(2))
        axs[i_ax].plot(1:length(history_fgs), history_fgs, color=cmap(4))

    end
    axs[3].legend(["MF-GMRES", "FGS"])

    for i in 0:3
        # axs[i].set_ylim([1,24])
        # axs[i].set_yticks([5,10,15,20])
        # axs[i].set_xscale("log")
        axs[i].set_xlim([1,maxiter])
        axs[i].set_yscale("log")
        axs[i].set_ylim([1e-7, 1e1])
    end
    tight_layout()

    savefig(save_name*"_smart.png")
end

function plot_warm_start_final_time_justfgs(file, smart_file, save_name; i_tol, i_m, i_t, i_tol_smart, i_m_smart, i_t_smart, i_St)

    # load data
    BSON.@load file ws_gmres_mf_a ws_gmres_mf_a_cold ws_fgs_a ws_fgs_a_cold histories_gmres_mf_a histories_gmres_mf_a_cold histories_fgs_a histories_fgs_a_cold n_panels tolerances ps_a ls_a mts_a ps_b ls_b mts_b nc0 ns0 span AR relaxation St Δt f

    # prepare colors
    n_colors = 6
    cmap = PythonPlot.get_cmap("viridis_r", n_colors)

    # generate figure
    fig = figure("history_final")
    fig.clf()
    fig.add_subplot(111, xlabel="timestep", ylabel="time cost, seconds")
    ax = fig.get_axes()[0]

    # loop over St
    maxiter = 0
    for i in eachindex(St)
        maxiter = max(maxiter, length(histories_gmres_mf_a_cold[i_tol,i_m][i_t,i]))
    end

    TFGS = nothing

    maxt= 0
    for i in eachindex(St)
        maxt = max(maxt, maximum(ws_gmres_mf_a_cold[:,1,i,i_tol,i_m]))
        maxt = max(maxt, maximum(ws_fgs_a_cold[:,1,i,i_tol,i_m]))
    end
    for (i_st, st) in enumerate(St)
        if i_st == i_St
            # plot GMRES

            ## get time
            t_cold = ws_gmres_mf_a_cold[:,1,i_st,i_tol,i_m]
            t_cold_fgs = ws_fgs_a_cold[:,1,i_st,i_tol,i_m]

            ## get histories
            niter_cold = ws_gmres_mf_a_cold[:,4,i_st,i_tol,i_m]
            niter_cold_fgs = ws_fgs_a_cold[:,4,i_st,i_tol,i_m]

            ## other
            i_ax = i_st-1
            ax.plot(0:length(niter_cold_fgs), vcat([t_cold_fgs[1]], t_cold_fgs), color=cmap(4))

            TFGS = t_cold_fgs[1]
        end
    end
    #axs[3].legend(["MF-GMRES", "FGS"])

    for i in 0:3
        # axs[i].set_ylim([1,24])
        # axs[i].set_yticks([5,10,15,20])
        # axs[i].set_xscale("log")
    end
    tight_layout()

    #savefig(save_name*"_cold.png")

    # loop over St
    maxiter = 0
    for i in eachindex(St)
        maxiter = max(maxiter, length(histories_gmres_mf_a[i_tol,i_m][i_t,i]))
    end
    maxt= 0
    for i in eachindex(St)
        maxt = max(maxt, maximum(ws_gmres_mf_a[:,1,i,i_tol,i_m]))
        maxt = max(maxt, maximum(ws_fgs_a[:,1,i,i_tol,i_m]))
    end
    for (i_st, st) in enumerate(St)
        if i_st == i_St
            # plot GMRES

            ## get time
            t = ws_gmres_mf_a[:,1,i_st,i_tol,i_m]
            t_fgs = ws_fgs_a[:,1,i_st,i_tol,i_m]

            ## get histories
            niter = ws_gmres_mf_a[:,4,i_st,i_tol,i_m]
            niter_fgs = ws_fgs_a[:,4,i_st,i_tol,i_m]

            ax.plot(0:length(niter_fgs), vcat(TFGS,t_fgs), "--", color=cmap(4))
        end
    end

    # load data
    BSON.@load smart_file ws_gmres_mf_a ws_gmres_mf_a_cold ws_fgs_a ws_fgs_a_cold histories_gmres_mf_a histories_gmres_mf_a_cold histories_fgs_a histories_fgs_a_cold n_panels tolerances ps_a ls_a mts_a ps_b ls_b mts_b nc0 ns0 span AR relaxation St Δt f

    # loop over St
    maxiter = 0
    for i in eachindex(St)
        maxiter = max(maxiter, length(histories_gmres_mf_a[i_tol_smart,i_m_smart][i_t_smart,i]))
    end
    maxt= 0
    for i in eachindex(St)
        maxt = max(maxt, maximum(ws_gmres_mf_a[:,1,i,i_tol_smart,i_m_smart]))
        maxt = max(maxt, maximum(ws_fgs_a[:,1,i,i_tol_smart,i_m_smart]))
    end
    for (i_st, st) in enumerate(St)
        if i_st == i_St
            # plot GMRES

            ## get time
            t = ws_gmres_mf_a[:,1,i_st,i_tol_smart,i_m_smart]
            t_fgs = ws_fgs_a[:,1,i_st,i_tol_smart,i_m_smart]

            ## get histories
            niter = ws_gmres_mf_a[:,4,i_st,i_tol_smart,i_m_smart]
            niter_fgs = ws_fgs_a[:,4,i_st,i_tol_smart,i_m_smart]

            i_ax = i_st-1
            ax.plot(0:length(niter_fgs), vcat(TFGS,t_fgs), ":", color=cmap(4))
        end
    end

    tight_layout()

    savefig(save_name*"_smart_justfgs_St$(i_St).pdf")
end

function plot_warm_start_final_time(file, smart_file, save_name; i_tol, i_m, i_t, i_tol_smart, i_m_smart, i_t_smart)

    # load data
    BSON.@load file ws_gmres_mf_a ws_gmres_mf_a_cold ws_fgs_a ws_fgs_a_cold histories_gmres_mf_a histories_gmres_mf_a_cold histories_fgs_a histories_fgs_a_cold n_panels tolerances ps_a ls_a mts_a ps_b ls_b mts_b nc0 ns0 span AR relaxation St Δt f

    # prepare colors
    n_colors = 6
    cmap = PythonPlot.get_cmap("viridis_r", n_colors)

    # generate figure
    fig = figure("history_final")
    fig.clf()
    fig.add_subplot(141, xlabel="timestep", ylabel="time cost, seconds")
    fig.add_subplot(142, xlabel="timestep")#, ylabel="L-2 norm residual")
    fig.add_subplot(143, xlabel="timestep")#, ylabel="L-2 norm residual")
    fig.add_subplot(144, xlabel="timestep")#, ylabel="L-2 norm residual")
    axs = fig.get_axes()

    # loop over St
    maxiter = 0
    for i in eachindex(St)
        maxiter = max(maxiter, length(histories_gmres_mf_a_cold[i_tol,i_m][i_t,i]))
    end

    maxt= 0
    for i in eachindex(St)
        maxt = max(maxt, maximum(ws_gmres_mf_a_cold[:,1,i,i_tol,i_m]))
        maxt = max(maxt, maximum(ws_fgs_a_cold[:,1,i,i_tol,i_m]))
    end
    for (i_st, st) in enumerate(St)
        # plot GMRES

        ## get time
        t_cold = ws_gmres_mf_a_cold[:,1,i_st,i_tol,i_m]
        t_cold_fgs = ws_fgs_a_cold[:,1,i_st,i_tol,i_m]

        ## get histories
        niter_cold = ws_gmres_mf_a_cold[:,4,i_st,i_tol,i_m]
        niter_cold_fgs = ws_fgs_a_cold[:,4,i_st,i_tol,i_m]

        ## other
        i_ax = i_st-1
        axs[i_ax].plot(1:length(niter_cold), t_cold, color=cmap(2))
        axs[i_ax].plot(1:length(niter_cold_fgs), t_cold_fgs, color=cmap(4))

    end
    axs[3].legend(["MF-GMRES", "FGS"])

    for i in 0:3
        # axs[i].set_ylim([1,24])
        # axs[i].set_yticks([5,10,15,20])
        # axs[i].set_xscale("log")
        axs[i].set_ylim([0,maxiter])
    end
    tight_layout()

    savefig(save_name*"_cold.png")

    #--- warm start ---#
    fig = figure("history_warm")
    fig.clf()
    fig.add_subplot(141, xlabel="timestep", ylabel="time cost, seconds")
    fig.add_subplot(142, xlabel="timestep")#, ylabel="L-2 norm residual")
    fig.add_subplot(143, xlabel="timestep")#, ylabel="L-2 norm residual")
    fig.add_subplot(144, xlabel="timestep")#, ylabel="L-2 norm residual")
    axs = fig.get_axes()

    # loop over St
    maxiter = 0
    for i in eachindex(St)
        maxiter = max(maxiter, length(histories_gmres_mf_a[i_tol,i_m][i_t,i]))
    end
    maxt= 0
    for i in eachindex(St)
        maxt = max(maxt, maximum(ws_gmres_mf_a[:,1,i,i_tol,i_m]))
        maxt = max(maxt, maximum(ws_fgs_a[:,1,i,i_tol,i_m]))
    end
    for (i_st, st) in enumerate(St)
        # plot GMRES

        ## get time
        t = ws_gmres_mf_a[:,1,i_st,i_tol,i_m]
        t_fgs = ws_fgs_a[:,1,i_st,i_tol,i_m]

        ## get histories
        niter = ws_gmres_mf_a[:,4,i_st,i_tol,i_m]
        niter_fgs = ws_fgs_a[:,4,i_st,i_tol,i_m]

        i_ax = i_st-1
        axs[i_ax].plot(1:length(niter), t, color=cmap(2))
        axs[i_ax].plot(1:length(niter_fgs), t_fgs, color=cmap(4))

    end
    axs[3].legend(["MF-GMRES", "FGS"])

    for i in 0:3
        # axs[i].set_ylim([1,24])
        # axs[i].set_yticks([5,10,15,20])
        # axs[i].set_xscale("log")
        axs[i].set_ylim([0,maxiter])
    end
    tight_layout()

    savefig(save_name*"_warm.png")

    #--- smart warm start ---#

    # load data
    BSON.@load smart_file ws_gmres_mf_a ws_gmres_mf_a_cold ws_fgs_a ws_fgs_a_cold histories_gmres_mf_a histories_gmres_mf_a_cold histories_fgs_a histories_fgs_a_cold n_panels tolerances ps_a ls_a mts_a ps_b ls_b mts_b nc0 ns0 span AR relaxation St Δt f

    @show tolerances[i_tol_smart], size(ws_gmres_mf_a_cold,5)
    fig = figure("history_smart")
    fig.clf()
    fig.add_subplot(141, xlabel="timestep", ylabel="time cost, seconds")
    fig.add_subplot(142, xlabel="timestep")#, ylabel="L-2 norm residual")
    fig.add_subplot(143, xlabel="timestep")#, ylabel="L-2 norm residual")
    fig.add_subplot(144, xlabel="timestep")#, ylabel="L-2 norm residual")
    axs = fig.get_axes()

    # loop over St
    maxiter = 0
    for i in eachindex(St)
        maxiter = max(maxiter, length(histories_gmres_mf_a[i_tol_smart,i_m_smart][i_t_smart,i]))
    end
    maxt= 0
    for i in eachindex(St)
        maxt = max(maxt, maximum(ws_gmres_mf_a[:,1,i,i_tol_smart,i_m_smart]))
        maxt = max(maxt, maximum(ws_fgs_a[:,1,i,i_tol_smart,i_m_smart]))
    end
    for (i_st, st) in enumerate(St)
        # plot GMRES

        ## get time
        t = ws_gmres_mf_a[:,1,i_st,i_tol_smart,i_m_smart]
        t_fgs = ws_fgs_a[:,1,i_st,i_tol_smart,i_m_smart]

        ## get histories
        niter = ws_gmres_mf_a[:,4,i_st,i_tol_smart,i_m_smart]
        niter_fgs = ws_fgs_a[:,4,i_st,i_tol_smart,i_m_smart]

        i_ax = i_st-1
        axs[i_ax].plot(1:length(niter), t, color=cmap(2))
        axs[i_ax].plot(1:length(niter_fgs), t_fgs, color=cmap(4))

    end
    axs[3].legend(["MF-GMRES", "FGS"])

    for i in 0:3
        # axs[i].set_ylim([1,24])
        # axs[i].set_yticks([5,10,15,20])
        # axs[i].set_xscale("log")
        axs[i].set_ylim([0,maxiter])
    end
    tight_layout()

    savefig(save_name*"_smart.png")
end

function plot_warm_start_final(file, smart_file, save_name; i_tol, i_m, i_t, i_tol_smart, i_m_smart, i_t_smart)

    # load data
    BSON.@load file ws_gmres_mf_a ws_gmres_mf_a_cold ws_fgs_a ws_fgs_a_cold histories_gmres_mf_a histories_gmres_mf_a_cold histories_fgs_a histories_fgs_a_cold n_panels tolerances ps_a ls_a mts_a ps_b ls_b mts_b nc0 ns0 span AR relaxation St Δt f

    # prepare colors
    n_colors = 6
    cmap = PythonPlot.get_cmap("viridis_r", n_colors)

    # generate figure
    fig = figure("history_final")
    fig.clf()
    fig.add_subplot(141, xlabel="timestep", ylabel="iterations to convergence")
    fig.add_subplot(142, xlabel="timestep")#, ylabel="L-2 norm residual")
    fig.add_subplot(143, xlabel="timestep")#, ylabel="L-2 norm residual")
    fig.add_subplot(144, xlabel="timestep")#, ylabel="L-2 norm residual")
    axs = fig.get_axes()

    # loop over St
    maxiter = 0
    for i in eachindex(St)
        maxiter = max(maxiter, length(histories_gmres_mf_a_cold[i_tol,i_m][i_t,i]))
    end

    for (i_st, st) in enumerate(St)
        # plot GMRES

        ## get time
        t_cold = ws_gmres_mf_a_cold[i_t,1,i_st,i_tol,i_m]
        t_cold_fgs = ws_fgs_a_cold[i_t,1,i_st,i_tol,i_m]

        ## get histories
        niter_cold = ws_gmres_mf_a_cold[:,4,i_st,i_tol,i_m]
        niter_cold_fgs = ws_fgs_a_cold[:,4,i_st,i_tol,i_m]

        ## other
        i_ax = i_st-1
        axs[i_ax].plot(1:length(niter_cold), niter_cold, color=cmap(2))
        axs[i_ax].plot(1:length(niter_cold_fgs), niter_cold_fgs, color=cmap(4))

    end
    axs[3].legend(["MF-GMRES", "FGS"])

    for i in 0:3
        # axs[i].set_ylim([1,24])
        # axs[i].set_yticks([5,10,15,20])
        # axs[i].set_xscale("log")
        axs[i].set_ylim([0,maxiter])
    end
    tight_layout()

    savefig(save_name*"_cold.png")

    #--- warm start ---#
    fig = figure("history_warm")
    fig.clf()
    fig.add_subplot(141, xlabel="timestep", ylabel="iterations to convergence")
    fig.add_subplot(142, xlabel="timestep")#, ylabel="L-2 norm residual")
    fig.add_subplot(143, xlabel="timestep")#, ylabel="L-2 norm residual")
    fig.add_subplot(144, xlabel="timestep")#, ylabel="L-2 norm residual")
    axs = fig.get_axes()

    # loop over St
    maxiter = 0
    for i in eachindex(St)
        maxiter = max(maxiter, length(histories_gmres_mf_a[i_tol,i_m][i_t,i]))
    end
    for (i_st, st) in enumerate(St)
        # plot GMRES

        ## get time
        t = ws_gmres_mf_a[i_t,1,i_st,i_tol,i_m]
        t_fgs = ws_fgs_a[i_t,1,i_st,i_tol,i_m]

        ## get histories
        niter = ws_gmres_mf_a[:,4,i_st,i_tol,i_m]
        niter_fgs = ws_fgs_a[:,4,i_st,i_tol,i_m]

        i_ax = i_st-1
        axs[i_ax].plot(1:length(niter), niter, color=cmap(2))
        axs[i_ax].plot(1:length(niter_fgs), niter_fgs, color=cmap(4))

    end
    axs[3].legend(["MF-GMRES", "FGS"])

    for i in 0:3
        # axs[i].set_ylim([1,24])
        # axs[i].set_yticks([5,10,15,20])
        # axs[i].set_xscale("log")
        axs[i].set_ylim([0,maxiter])
    end
    tight_layout()

    savefig(save_name*"_warm.png")

    #--- smart warm start ---#

    # load data
    BSON.@load smart_file ws_gmres_mf_a ws_gmres_mf_a_cold ws_fgs_a ws_fgs_a_cold histories_gmres_mf_a histories_gmres_mf_a_cold histories_fgs_a histories_fgs_a_cold n_panels tolerances ps_a ls_a mts_a ps_b ls_b mts_b nc0 ns0 span AR relaxation St Δt f

    @show tolerances[i_tol_smart], size(ws_gmres_mf_a_cold,5)
    fig = figure("history_smart")
    fig.clf()
    fig.add_subplot(141, xlabel="timestep", ylabel="iterations to convergence")
    fig.add_subplot(142, xlabel="timestep")#, ylabel="L-2 norm residual")
    fig.add_subplot(143, xlabel="timestep")#, ylabel="L-2 norm residual")
    fig.add_subplot(144, xlabel="timestep")#, ylabel="L-2 norm residual")
    axs = fig.get_axes()

    # loop over St
    maxiter = 0
    for i in eachindex(St)
        maxiter = max(maxiter, length(histories_gmres_mf_a[i_tol_smart,i_m_smart][i_t_smart,i]))
    end
    for (i_st, st) in enumerate(St)
        # plot GMRES

        ## get time
        t = ws_gmres_mf_a[i_t,1,i_st,i_tol_smart,i_m_smart]
        t_fgs = ws_fgs_a[i_t,1,i_st,i_tol_smart,i_m_smart]

        ## get histories
        niter = ws_gmres_mf_a[:,4,i_st,i_tol_smart,i_m_smart]
        niter_fgs = ws_fgs_a[:,4,i_st,i_tol_smart,i_m_smart]

        i_ax = i_st-1
        axs[i_ax].plot(1:length(niter), niter, color=cmap(2))
        axs[i_ax].plot(1:length(niter_fgs), niter_fgs, color=cmap(4))

    end
    axs[3].legend(["MF-GMRES", "FGS"])

    for i in 0:3
        # axs[i].set_ylim([1,24])
        # axs[i].set_yticks([5,10,15,20])
        # axs[i].set_xscale("log")
        axs[i].set_ylim([0,maxiter])
    end
    tight_layout()

    savefig(save_name*"_smart.png")
end


# set save path
warm_start_dir = "warm_start_20240629_2"
smart_dir = "warm_start_20240629_2_smart_guess"

# get results
St = [1e-3, 1e-2, 1e-1, 1e0]
nt = 10
ms = 2:2
ms_smart = 2:2
t0 = 0.25
smart_guess = true
# ws = warm_start_wing(ms, warm_start_dir; St, replace_file=false)

# ws_cycle = warm_start_wing_cycle(ms, warm_start_dir; St, replace_file=true, nt, t0, smart_guess)

# plot
cycle_dir = joinpath(FLOWPanel_dir, warm_start_dir)
smart_cycle_dir = joinpath(FLOWPanel_dir, smart_dir)
file = joinpath(cycle_dir, "assembled_benchmarks_cycle.bson")
smart_file = joinpath(smart_cycle_dir, "assembled_benchmarks_cycle.bson")
i_m, i_tol, i_t = 2, 2, 2
i_m_smart, i_tol_smart, i_t_smart = 1, 2, 2

save_name = joinpath(cycle_dir, "cycle_m$(i_m)_itol$(i_tol)_time")
# plot_warm_start_cycle(file, save_name; i_tol, i_m)
# plot_warm_start_cycle_individual(file, save_name; i_tol, i_m)

plot_warm_start_final_time(file, smart_file, save_name; i_tol, i_m, i_t, i_tol_smart, i_m_smart, i_t_smart)
plot_warm_start_final_time_justfgs(file, smart_file, save_name; i_tol, i_m, i_t, i_tol_smart, i_m_smart, i_t_smart, i_St=1)
plot_warm_start_final_time_justfgs(file, smart_file, save_name; i_tol, i_m, i_t, i_tol_smart, i_m_smart, i_t_smart, i_St=4)

save_name = joinpath(cycle_dir, "cycle_m$(i_m)_itol$(i_tol)_it$(i_t)_history")
#plot_warm_start_cycle_history(file, save_name; i_tol, i_m, i_t)
#plot_warm_start_cycle_history_individual(file, save_name; i_tol, i_m, i_t)

#plot_warm_start_cycle_history_final(file, smart_file, save_name; i_tol, i_m, i_t, i_tol_smart, i_m_smart, i_t_smart)
