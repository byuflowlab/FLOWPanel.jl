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

function benchmark_wing_lu(wing, freestream=SVector{3}(1.0,0,0))

    # reset wing
    reset_potential_velocity!(wing)
    apply_freestream!(wing, freestream) # apply freestream
    FLOWPanel.set_unit_strength!(wing)

    # create solver
    scheme = FLOWPanel.Scheme{FLOWPanel.DirectNeumann, FlowTangency}
    solver, t_aic, t_lu, t_alloc = FLOWPanel.LUDecomposition_benchmark(wing, scheme)
    print(" Created solver.")

    # solve
    t_solve = @elapsed solve!(wing, solver)
    print(" Finished solve.")

    return t_aic+t_alloc, t_lu, t_solve
end

function benchmark_wing_gmres(wing, tolerance, max_iterations, freestream=SVector{3}(1.0,0,0); cold_start=true)

    # reset wing
    reset_potential_velocity!(wing) # zero velocity
    apply_freestream!(wing, freestream) # apply freestream
    cold_start && FLOWPanel.set_unit_strength!(wing) # cold start initial guess

    # create solver
    scheme = FLOWPanel.Scheme{FLOWPanel.DirectNeumann, FlowTangency}
    solver, t_aic, t_alloc = FLOWPanel.IterativeSolver_benchmark(wing, scheme)
    print(" Created solver.")

    # solve
    t_solve = @elapsed solve!(wing, solver; verbose=false, tolerance, max_iterations, rtol=0.0)
    print(" Finished solve.")

    return t_aic+t_alloc, t_solve
end

function benchmark_wing_gmres_mf(wing, tolerance, max_iterations, expansion_order, leaf_size, multipole_threshold, reuse_tree::Bool, freestream=SVector{3}(1.0,0,0); cold_start=true)

    # reset wing
    reset_potential_velocity!(wing)
    apply_freestream!(wing, freestream) # apply freestream
    cold_start && FLOWPanel.set_unit_strength!(wing)

    # create solver
    scheme = FLOWPanel.Scheme{FLOWPanel.DirectNeumann, FlowTangency}
    solver, t_alloc = FLOWPanel.MatrixFreeSolver_benchmark(wing, scheme; expansion_order, leaf_size, multipole_threshold, reuse_tree)
    print(" Created solver.")

    # solve
    t_solve = @elapsed solve!(wing, solver; verbose=false, tolerance, max_iterations, rtol=0.0)
    print(" Finished solve.")

    return t_alloc, t_solve
end

function benchmark_wing_fgs(wing, tolerance, max_iterations, expansion_order, leaf_size, multipole_threshold, relaxation, reuse_tree::Bool, freestream=SVector{3}(1.0,0,0); cold_start=true)

    # reset wing
    reset_potential_velocity!(wing)
    apply_freestream!(wing, freestream) # apply freestream
    cold_start && FLOWPanel.set_unit_strength!(wing)

    # create solver
    scheme = FLOWPanel.Scheme{FLOWPanel.DirectNeumann, FlowTangency}
    solver, t_solver, t_fmm = FLOWPanel.FastGaussSeidel_benchmark(wing, scheme; expansion_order, leaf_size, multipole_threshold, reuse_tree)
    print(" Created solver.")

    # solve
    this_tolerance = tolerance
    t_solve = @elapsed solve!(wing, solver; verbose=false, tolerance=this_tolerance, max_iterations, relaxation)
    print(" Finished solve.")

    return t_solver, t_fmm, t_solve
end

function benchmark_wing(wing::FLOWPanel.AbstractPanels, dummy_velocity, tolerance, expansion_order_a, leaf_size_a, multipole_threshold_a, expansion_order_b, leaf_size_b, multipole_threshold_b, relaxation, nc, ns, t_aic_lu, t_lu_lu, t_solve_lu, resid_max_lu, resid_mean_lu;
        freestream=SVector{3}(1.0,0,0),
        benchmark_dir="benchmarks_20240625",
    )

    # benchmark GMRES
    print("\tGMRES...")
    max_iterations = 100
    t_aic_gmres, t_solve_gmres = benchmark_wing_gmres(wing, tolerance, max_iterations, freestream)

    ## get residual
    resid_max_gmres, resid_mean_gmres = get_residual(dummy_velocity, wing, freestream)

    ## save vtk
    vtk(joinpath(FLOWPanel_dir, benchmark_dir, "benchmark_gmres", "benchmark_gmres_nc$(nc)_ns$(ns)_atol$(tolerance).vts"), wing)

    # benchmark MF-GMRES-a
    print("\n\tMF-GMRES-a...")
    reuse_tree = false
    t_alloc_gmres_mf_a, t_solve_gmres_mf_a = benchmark_wing_gmres_mf(wing, tolerance, max_iterations, expansion_order_a, leaf_size_a, multipole_threshold_a, reuse_tree, freestream)

    ## get residual
    resid_max_gmres_mf_a, resid_mean_gmres_mf_a = get_residual(dummy_velocity, wing, freestream)

    ## save vtk
    vtk(joinpath(FLOWPanel_dir, benchmark_dir, "benchmark_gmres_mf_a", "benchmark_gmres_mf_a_nc$(nc)_ns$(ns)_atol$(tolerance)_p$(expansion_order_a)_l$(leaf_size_a)_$(multipole_threshold_a).vts"), wing)

    #=
    # benchmark MF-GMRES-b
    print("\n\tMF-GMRES-b...")
    reuse_tree = true
    t_alloc_gmres_mf_b, t_solve_gmres_mf_b = benchmark_wing_gmres_mf(wing, tolerance, max_iterations, expansion_order_b, leaf_size_b, multipole_threshold_b, reuse_tree, freestream)

    ## get residual
    resid_max_gmres_mf_b, resid_mean_gmres_mf_b = get_residual(dummy_velocity, wing, freestream)

    ## save vtk
    vtk(joinpath(FLOWPanel_dir, benchmark_dir, "benchmark_gmres_mf_b", "benchmark_gmres_mf_b_nc$(nc)_ns$(ns)_atol$(tolerance)_p$(expansion_order_b)_l$(leaf_size_b)_$(multipole_threshold_b).vts"), wing)
    =#
    t_alloc_gmres_mf_b, t_solve_gmres_mf_b, resid_max_gmres_mf_b, resid_mean_gmres_mf_b = NaN, NaN, NaN, NaN

    # benchmark FGS-a
    print("\n\tFast Gauss-Seidel-a...")
    max_iterations=500
    reuse_tree = false
    t_solver_fgs_a, t_fmm_fgs_a, t_solve_fgs_a = benchmark_wing_fgs(wing, tolerance, max_iterations, expansion_order_a, leaf_size_a, multipole_threshold_a, relaxation, reuse_tree, freestream)

    ## get residual
    resid_max_fgs_a, resid_mean_fgs_a = get_residual(dummy_velocity, wing, freestream)

    ## save vtk
    vtk(joinpath(FLOWPanel_dir, benchmark_dir, "benchmark_fgs_a", "benchmark_fgs_a_nc$(nc)_ns$(ns)_atol$(tolerance)_p$(expansion_order_a)_l$(leaf_size_a)_$(multipole_threshold_a)_relax$(relaxation).vts"), wing)

    # benchmark FGS-b
    #=
    print("\n\tFast Gauss-Seidel-b...")
    max_iterations=500
    reuse_tree = true
    t_solver_fgs_b, t_fmm_fgs_b, t_solve_fgs_b = benchmark_wing_fgs(wing, tolerance, max_iterations, expansion_order_b, leaf_size_b, multipole_threshold_b, relaxation, reuse_tree, freestream)

    ## get residual
    resid_max_fgs_b, resid_mean_fgs_b = get_residual(dummy_velocity, wing, freestream)

    ## save vtk
    vtk(joinpath(FLOWPanel_dir, benchmark_dir, "benchmark_fgs_b", "benchmark_fgs_b_nc$(nc)_ns$(ns)_atol$(tolerance)_p$(expansion_order_b)_l$(leaf_size_b)_$(multipole_threshold_b)_relax$(relaxation).vts"), wing)
    println()
    =#

    t_solver_fgs_b, t_fmm_fgs_b, t_solve_fgs_b, resid_max_fgs_b, resid_mean_fgs_b = NaN, NaN, NaN, NaN, NaN

    return (t_aic_lu, t_lu_lu, t_solve_lu, resid_max_lu, resid_mean_lu), (t_aic_gmres, t_solve_gmres, resid_max_gmres, resid_mean_gmres), (t_alloc_gmres_mf_a, t_solve_gmres_mf_a, resid_max_gmres_mf_a, resid_mean_gmres_mf_a), (t_alloc_gmres_mf_b, t_solve_gmres_mf_b, resid_max_gmres_mf_b, resid_mean_gmres_mf_b), (t_solver_fgs_a, t_fmm_fgs_a, t_solve_fgs_a, resid_max_fgs_a, resid_mean_fgs_a), (t_solver_fgs_b, t_fmm_fgs_b, t_solve_fgs_b, resid_max_fgs_b, resid_mean_fgs_b)
end

function benchmark_wing(ms, benchmark_dir;
        replace_file=false, omit_lu=false,
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
    for dir in ("benchmark_lu", "benchmark_gmres", "benchmark_gmres_mf_a", "benchmark_gmres_mf_b", "benchmark_fgs_a", "benchmark_fgs_b")
        if !isdir(joinpath(FLOWPanel_dir, benchmark_dir, dir))
            mkdir(joinpath(FLOWPanel_dir, benchmark_dir, dir))
        end
    end

    # check if file exists
    file = joinpath(FLOWPanel_dir, benchmark_dir, "assembled_benchmarks.bson")
    if isfile(file) && !replace_file
        BSON.@load file bms_lu bms_gmres bms_gmres_mf_a bms_gmres_mf_b bms_fgs_a bms_fgs_b n_panels tolerances ps_a ls_a mts_a ps_b ls_b mts_b nc0 ns0 span AR relaxation
    else
        # preallocate benchmarks
        bms_lu = zeros(5, length(tolerances), length(ms))
        bms_gmres = zeros(4, length(tolerances), length(ms))
        bms_gmres_mf_a = zeros(4, length(tolerances), length(ms))
        bms_gmres_mf_b = zeros(4, length(tolerances), length(ms))
        bms_fgs_a = zeros(5, length(tolerances), length(ms))
        bms_fgs_b = zeros(5, length(tolerances), length(ms))
        n_panels = fill(-1, length(ms))

        # populate with NaNs
        bms_lu .= NaN
        bms_gmres .= NaN
        bms_gmres_mf_a .= NaN
        bms_gmres_mf_b .= NaN
        bms_fgs_a .= NaN
        bms_fgs_b .= NaN
    end

    # run benchmarks
    for (i_m, (m, p_a, l_a, mt_a, p_b, l_b, mt_b)) in enumerate(zip(ms, ps_a, ls_a, mts_a, ps_b, ls_b, mts_b))

        # create wing
        nc, ns = nc0*m, ns0*m
        println("\n#------- m = $m; n_panels = $(ns0*nc0*2*m^2) -------#")
        print("\n\tpreparing wing...")
        wing = prepare_wing(; kernel, span, AR, nc, ns, freestream)
        print(" Done.\n")
        n_panels[i_m] = length(wing.panels)
        dummy_velocity = zeros(SVector{3,Float64}, length(wing.panels))

        # benchmark LU
        print("\n\tLU Decomposition...")
        if omit_lu
            # bms_lu[:,i_tol,i_m] .= bm_lu
            t_aic_lu, t_lu_lu, t_solve_lu, resid_max_lu, resid_mean_lu = bms_lu[:,1,i_m]
        else
            t_aic_lu, t_lu_lu, t_solve_lu = benchmark_wing_lu(wing, freestream)

            ## get residual
            resid_max_lu, resid_mean_lu = get_residual(dummy_velocity, wing, freestream)

            ## save vtk
            vtk(joinpath(FLOWPanel_dir, benchmark_dir, "benchmark_lu", "benchmark_lu_nc$(nc)_ns$(ns).vts"), wing)
        end
        println()

        # loop over tolerances
        for (i_tol,tol) in enumerate(tolerances)

            println("\n#--- m = $m; tolerance: $tol ---#\n")

            # run benchmarks
            bm_lu, bm_gmres, bm_gmres_mf_a, bm_gmres_mf_b, bm_fgs_a, bm_fgs_b = benchmark_wing(wing, dummy_velocity, tol, p_a[i_tol], l_a[i_tol], mt_a[i_tol], p_b[i_tol], l_b[i_tol], mt_b[i_tol], relaxation, nc, ns, t_aic_lu, t_lu_lu, t_solve_lu, resid_max_lu, resid_mean_lu; freestream, benchmark_dir)

            # store results
            bms_lu[:,i_tol,i_m] .= bm_lu
            bms_gmres[:,i_tol,i_m] .= bm_gmres
            bms_gmres_mf_a[:,i_tol,i_m] .= bm_gmres_mf_a
            bms_gmres_mf_b[:,i_tol,i_m] .= bm_gmres_mf_b
            bms_fgs_a[:,i_tol,i_m] .= bm_fgs_a
            bms_fgs_b[:,i_tol,i_m] .= bm_fgs_b

            # store bson
            BSON.@save joinpath(FLOWPanel_dir, benchmark_dir, "assembled_benchmarks.bson") bms_lu bms_gmres bms_gmres_mf_a bms_gmres_mf_b bms_fgs_a bms_fgs_b n_panels tolerances ps_a ls_a mts_a ps_b ls_b mts_b nc0 ns0 span AR relaxation
        end
    end

    return bms_lu, bms_gmres, bms_gmres_mf_a, bms_gmres_mf_b, bms_fgs_a, bms_fgs_b, n_panels

end

# function benchmark_convergence(

function plot_benchmark_1b(file, save_name)

    # load data
    BSON.@load file bms_lu bms_gmres bms_gmres_mf_a bms_gmres_mf_b bms_fgs_a bms_fgs_b n_panels tolerances ps_a ls_a mts_a ps_b ls_b mts_b nc0 ns0 span AR

    # prepare colors
    n_colors = 6
    cmap = PythonPlot.get_cmap("viridis_r", n_colors)

    # loop over tolerances
    for (i_tol, tolerance) in enumerate(tolerances)

        fig = figure("benchmark_wing_1b_$i_tol")
        fig.clf()
        fig.add_subplot(211, ylabel="total cost, seconds")
        fig.add_subplot(212, ylabel=L"L^2"*" norm error, meters per second", xlabel="number of panels")
        axs = fig.get_axes()

        # plot total cost for each solver

        ## lu decomposition
        axs[0].plot(n_panels, bms_lu[1,i_tol,:] + bms_lu[2,i_tol,:] + bms_lu[3,i_tol,:], color=cmap(0), "-", marker=".")
        axs[1].plot(n_panels, bms_lu[5,i_tol,:], color=cmap(0), "-", marker=".")

        ## gmres
        axs[0].plot(n_panels, bms_gmres[1,i_tol,:] + bms_gmres[2,i_tol,:], color=cmap(1), "-", marker="*")
        axs[1].plot(n_panels, bms_gmres[4,i_tol,:], color=cmap(1), "-", marker="*")

        ## gmres-mf-a
        axs[0].plot(n_panels, bms_gmres_mf_a[1,i_tol,:] + bms_gmres_mf_a[2,i_tol,:], color=cmap(2), "-", marker="x")
        axs[1].plot(n_panels, bms_gmres_mf_a[4,i_tol,:], color=cmap(2), "-", marker="x")

        ## gmres-mf-b
        # axs[0].plot(n_panels, bms_gmres_mf_b[1,i_tol,:] + bms_gmres_mf_b[2,i_tol,:], color=cmap(3), "--", marker="+")
        # axs[1].plot(n_panels, bms_gmres_mf_b[4,i_tol,:], color=cmap(3), "--", marker="+")

        ## fgs-a
        axs[0].plot(n_panels, bms_fgs_a[1,i_tol,:] + bms_fgs_a[2,i_tol,:] + bms_fgs_a[3,i_tol,:], color=cmap(4), "-", marker="1")
        axs[1].plot(n_panels, bms_fgs_a[5,i_tol,:], color=cmap(4), "-", marker="1")

        ## fgs-b
        # axs[0].plot(n_panels, bms_fgs_b[1,i_tol,:] + bms_fgs_b[2,i_tol,:] + bms_fgs_b[3,i_tol,:], color=cmap(5), "--", marker="2")
        # axs[1].plot(n_panels, bms_fgs_b[5,i_tol,:], color=cmap(5), "--", marker="2")

        # set axes
        axs[0].set_yscale("log")
        axs[1].set_yscale("log")
        axs[1].set_ylim([1e-10, 1e0])
        axs[0].set_xscale("log")
        axs[1].set_xscale("log")

        # create legend
        # axs[1].legend(["LU decomposition", "GMRES", "FMM+GMRES (a)", "FMM+GMRES (b)", "Fast Gauss-Seidel (a)", "Fast Gauss-Seidel (b)"])
        axs[1].legend(["LU decomposition", "GMRES", "MF-GMRES", "Fast Gauss-Seidel"])

        # tight layout
        fig.tight_layout()

        # save figure
        savefig(save_name*"_$i_tol.pdf")
    end
end

# set save path
# benchmark_dir = "benchmarks_20240630_10aoa_new_fmm_params_12"
benchmark_dir = "benchmarks_20240629_10aoa_1_12"

# get results
ms = vcat(collect(1:8),12)
# ms = [12]
# bm = benchmark_wing(vcat(collect(1:8),12), benchmark_dir; omit_lu=false, freestream=SVector{3}(1.0,0,-tan(10*pi/180)))
# benchmark_wing(ms, benchmark_dir; replace_file=false, omit_lu=true,
#         freestream=SVector{3}(1.0,0,-tan(10*pi/180)), kernel=ConstantSource(),
#         nc0=5, ns0=50, span=1.0, AR=10,
#         tolerances = [1e-3, 1e-6, 1e-9],
#         ps_a = fill((5,7,12), length(ms)),
#         ls_a = fill((20,30,100), length(ms)),
#         mts_a = fill((0.6,0.3,0.35), length(ms)),
#         ps_b = [(1,1,1), (3,1,1), (4,1,1), (4,10,6), (4,8,14), (4,8,14), (4,8,14), (4,8,14)],
#         ls_b = [(200,200,200), (100,800,800), (100,800,800), (100,400,800), (100,300,400), (100,300,400), (100,300,400), (100,300,400)],
#         mts_b = [(0.4,0.4,0.4), (0.55,0.6,0.6), (0.5,0.6,0.6), (0.5,0.6,0.4), (0.5,0.4,0.45), (0.5,0.4,0.45), (0.5,0.4,0.45), (0.5,0.4,0.45)],
#         relaxation = 1.4,
#     )
# # plot results
this_path = joinpath(FLOWPanel_dir, benchmark_dir)
file = joinpath(this_path, "assembled_benchmarks.bson")
save_name = joinpath(this_path, "assembled_benchmark")
plot_benchmark_1b(file, save_name)

