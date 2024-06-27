using Pkg
Pkg.activate("/Users/ryan/Dropbox/research/projects/FLOWPanel.jl/")
using FLOWPanel
using FLOWPanel.StaticArrays
include("auxiliary_functions.jl")

function benchmark_wing_solver(;
        kernel=ConstantSource(),
        solver_type=FLOWPanel.LUDecomposition_benchmark,
        AR=10, span=1.0,
        nc=10, ns=30,
        freestream=SVector{3}(1.0,0,0),
    )
    # create wing
    wing = prepare_wing(; kernel, span, AR, nc, ns, freestream)

    # create solver
    solver, t_solver_prep, t_solver_alloc = solver_type(wing, FLOWPanel.Scheme{FLOWPanel.DirectNeumann, FlowTangency})

    return wing, solver, t_solver_prep, t_solver_alloc
end

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

function benchmark_wing(vtk_base_name::String;
            solver_type = LUDecomposition_benchmark,
            multipliers = range(1,10),
            nc0=5, ns0=50,
            solver_kwargs...
       )

    ts_pre = zeros(length(multipliers))
    ts_alloc = zeros(length(multipliers))
    ts_solve = zeros(length(multipliers))
    for (i,m) in enumerate(multipliers)
        println("Benchmark: m=$m")
        wing, solver, t_pre, t_alloc = benchmark_wing_solver(; solver_type, nc=nc0*m, ns=ns0*m)
        t_solve = @elapsed solve!(wing, solver; solver_kwargs...)
        @show t_pre, t_alloc, t_solve
        ts_pre[i] = t_pre
        ts_alloc[i] = t_alloc
        ts_solve[i] = t_solve

        # save vtk file
        vtk("benchmark_$(vtk_base_name)_nc$(nc0*m)_ns$(ns0*m).vts", wing)
    end
    return ts_pre, ts_alloc, ts_solve
end

function benchmark_wing_iterative(vtk_base_name::String;
            solver_type = FLOWPanel.IterativeSolver_benchmark,
            multipliers = range(1,8),
            tols = [1e-3, 1e-6, 1e-9],
            nc0=5, ns0=50,
            solve_kwargs...
       )

    ts_pre = zeros(length(multipliers))
    ts_alloc = zeros(length(multipliers))
    ts_solve = zeros(length(multipliers), length(tols))

    for (i,m) in enumerate(multipliers)
        println("\nBenchmark: m=$m")
        wing, solver, t_pre, t_alloc = benchmark_wing_solver(; solver_type, nc=nc0*m, ns=ns0*m)
        @show t_pre t_alloc
        ts_pre[i] = t_pre
        ts_alloc[i] = t_alloc

        for (j,tol) in enumerate(tols)
            println("\ntol = $tol")
            t_solve = @elapsed solve!(wing, solver; tolerance=tol, max_iterations=200, solve_kwargs...)
            @show t_solve
            ts_solve[i,j] = t_solve
            # save vtk file
            vtk("benchmark_$(vtk_base_name)_nc$(nc0*m)_ns$(ns0*m)_tol$(tol).vts", wing)
        end

    end
    return ts_pre, ts_alloc, ts_solve
end

#=
#--- benchmark LU decomposition ---#
ts_pre_lu, ts_alloc_lu, ts_solve_lu = benchmark_wing("lu"; solver_type=FLOWPanel.LUDecomposition_benchmark, multipliers=range(1,8))

#--- benchmark matrix-powered GMRES ---#
ts_pre_gmres, ts_alloc_gmres, ts_solve_gmres = benchmark_wing_iterative("matrix_GMRES"; solver_type=FLOWPanel.IterativeSolver_benchmark, tols=[1e-3, 1e-6, 1e-9], multipliers=range(1,8))

#--- benchmark matrix-free GMRES ---#
ts_pre_gmres_mf, ts_alloc_gmres_mf, ts_solve_gmres_mf = benchmark_wing_iterative("matrix_free_GMRES"; solver_type=FLOWPanel.MatrixFreeSolver_benchmark, tols=[1e-3, 1e-6, 1e-9], multipliers=range(1,8))

#--- save file ---#
n_panels = range(1,10) .^2 .* 5*50*2
BSON.@save "benchmarks.bson" ts_pre_lu ts_alloc_lu ts_solve_lu ts_pre_gmres ts_alloc_gmres ts_solve_gmres ts_pre_gmres_mf ts_alloc_gmres_mf ts_solve_gmres_mf n_panels
=#
