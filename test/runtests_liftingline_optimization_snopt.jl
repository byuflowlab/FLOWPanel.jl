#=##############################################################################
# DESCRIPTION
    Unit tests of lifting line solver
=###############################################################################

using Test
import Printf: @printf
import FLOWPanel as pnl
import FLOWPanel: mean, norm, dot, cross

import ForwardDiff: Dual, Partials, value, partials, jacobian
import Snopt
import SNOW: minimize, Options, SNOPT, ForwardAD, ComplexStep

try
    verbose
catch
    global verbose = true
end
v_lvl = 0


@testset verbose=verbose "Lifting Line Optimization Tests w/ SNOPT" begin

    @testset "AOA vs twist optimization equivalence" begin

        # --------------- OPTIMIZATION TEST: Optimum AOA -----------------------
        if verbose
            println("\n"*"\t"^(v_lvl)*"Optimization test: Max L/D w.r.t. AOA")
        end

        model_cache = Dict()
        model_cache["fcalls"] = 0

        function f!(g, x; optimize_aoa=true, cache=model_cache)

            # Fetch design variables
            aoa = x[1]

            if optimize_aoa
                alpha = aoa
                twist_root = 0.0
                twist_tip = 0.0
            else
                alpha = 0.0
                twist_root = aoa
                twist_tip = aoa
            end

            magUinf = 49.7

            chord_root = 0.2
            chord_tip = 0.2

            sweep = 0.0

            dihedral = 0.0

            airfoil_path = joinpath(pnl.examples_path, "data")

            nelements = 40

            # Evaluate model
            out = pnl.run_liftingline(; 
                                    alpha,
                                    magUinf,
                                    nelements,
                                    verbose     = false,

                                    rho         = 1.225,
                                    b           = 98*0.0254,                    # (m) wing span
                                    
                                    # Chord distribution (nondim y-position, nondim chord)
                                    chord_distribution = [
                                    #   2*y/b   c/b
                                        0.0     chord_root
                                        1.0     chord_tip
                                    ],
                                    
                                    # Twist distribution (nondim y-position, twist)
                                    twist_distribution = [
                                    #   2*y/b   twist (deg)
                                        0.0     twist_root
                                        1.0     twist_tip
                                    ],
                                    
                                    # Sweep distribution (nondim y-position, sweep)
                                    sweep_distribution = [
                                    #   2*y/b   sweep (deg)
                                        0.0     sweep
                                        1.0     sweep
                                    ],
                                    
                                    # Dihedral distribution (nondim y-position, dihedral)
                                    dihedral_distribution = [
                                    #   2*y/b   dihedral (deg)
                                        0.0     dihedral
                                        1.0     dihedral
                                    ],
                                    
                                    # Airfoil contour distribution (nondim y-position, polar, airfoil type)
                                    airfoil_distribution = [
                                    #    2*y/b  polar file                        airfoil type
                                        (0.00, "rae101-Re1p7e6-smooth180-3.csv",  pnl.SimpleAirfoil),
                                        (1.00, "rae101-Re1p7e6-smooth180-3.csv",  pnl.SimpleAirfoil)
                                    ],
                
                                    airfoil_path,
                                    element_optargs = (;    path = airfoil_path,
                                                            plot_polars = !true,
                                                            extrapolatepolar = false,   # Whether to extrapolate the 2D polars to ±180 deg
                                                            verbose=false,
                                                        ),
                                    
                                    
                                    # ------------------ SOLVER PARAMETERS -----------------------------------------
                                    deltasb         = 1.0,                          # Blending distance, deltasb = 2*dy/b
                                    deltajoint      = 1.0,                          # Joint distance, deltajoint = dx/c
                                    
                                    sigmafactor     = 0.0,                          # Dragging line amplification factor (set to -1.0 for rebust post-stall method)
                                    sigmaexponent   = 4.0,                          # Dragging line amplification exponent (no effects if `sigmafactor==0.0`)
                                    
                                                                                    # Nonlinear solver
                                    solver          = pnl.SimpleNonlinearSolve.SimpleDFSane(),             # Indifferent to initial guess, but somewhat not robust
                                    
                                    solver_optargs  = (; 
                                                        abstol = 1e-12,  
                                                        maxiters = 800,
                                                        ),
                                    
                                    align_joints_with_Uinfs = false,                # Whether to align joint bound vortices with the freestream
                                    
                                    use_Uind_for_force = true,                      # Whether to use Uind as opposed to selfUind for force postprocessing
                                                                                    # (`true` for more accurate spanwise cd distribution, but worse integrated CD)
                                    
                                    stability_derivatives = !true,                   # Whether to calculate stability derivatives
                
                                    cache
                                    )

            cache["fcalls"] += 1

            # Fetch results
            L = out.L
            D = out.D
            
            # Objective (minimization)
            obj = -L/D

            # Constraints (dummies, we defined two or SNOW won't work)
            # g[1] = (-x[1] - 999) / 1e-3
            # g[2] = (-x[2] - 999) / 1e-3
            # g[2] = (-x[1] - 999) / 1e-3
            # g[1] = (-x[1] - 1) / 1e-3
            # g[2] = (-x[2] - 1) / 1e-3
            # g[1] = -L
            # g[2] = -L

            return obj
        end

        x0 = [0.0]             # Initial guess

        lx = [-10.0]           # Lower bounds on x
        ux = [10.0]            # Upper bounds on x

        ng = 0                      # Number of constraints
        lg = -Inf*ones(ng)          # Lower bounds on g
        ug = zeros(ng)              # Upper bounds on g

        snopt_options = Dict(
                "Major iterations limit" => 500,
            )
        solver = SNOPT(options=snopt_options)
        options = Options(; solver, derivatives=ForwardAD())
        # options = Options(; solver, derivatives=ComplexStep())
        # options = Options(; solver)

        model_cache = Dict()
        model_cache["fcalls"] = 0

        t_snopt = @elapsed begin
            xopt, fopt, info = minimize(f!, deepcopy(x0), ng, lx, ux, lg, ug, options)
        end

        aoaopt_snopt = sum(xopt)
        fopt_snopt = value(f!(zeros(ng), xopt))
        fcalls_snopt = model_cache["fcalls"]

        # Brute force answer
        model_cache = Dict()
        model_cache["fcalls"] = 0

        t_brute = @elapsed begin
            aoas = range(lx[1], ux[1], length=200)
            fs = [value(f!(zeros(2), [aoa, 0])) for aoa in aoas]
        end

        aoaopt_brute_i = argmin(fs)
        aoaopt_brute = aoas[aoaopt_brute_i]
        fopt_brute = fs[aoaopt_brute_i]
        fcalls_brute = model_cache["fcalls"]

        if verbose
            println()
            @printf "%s%15.15s   %-10s %-10s %-10s %-10s\n" "\t"^(v_lvl+1) "" "Opt AOA" "Max L/D" "Time" "Function calls"
            @printf "%s%15.15s   %-10.4f %-10.4f %-3.0f secs  %4.3g\n" "\t"^(v_lvl+1) "Brute force" aoaopt_brute -fopt_brute t_brute fcalls_brute
            @printf "%s%15.15s   %-10.4f %-10.4f %-3.0f secs  %4.3g\n" "\t"^(v_lvl+1) "Ipopt" aoaopt_snopt -fopt_snopt t_snopt fcalls_snopt
        end

        @testset "Max L/D w.r.t. alpha" begin
            @test info == :Solved_To_Acceptable_Level || info == :Solve_Succeeded
            @test isapprox(aoaopt_snopt, aoaopt_brute, atol=0.05)
            @test isapprox(fopt_snopt, fopt_brute, atol=0.01)
        end




        # # --------------- OPTIMIZATION TEST: Optimum twist ---------------------
        # if verbose
        #     println("\n"*"\t"^(v_lvl)*"Optimization test: Max L/D w.r.t. twist")
        # end

        # function f2!(args...; optargs...)
        #     return f!(args...; optargs..., optimize_aoa=false)
        # end

        # model_cache = Dict()
        # model_cache["fcalls"] = 0

        # t_snopt = @elapsed begin
        #     xopt, fopt, info = minimize(f2!, deepcopy(x0), ng, lx, ux, lg, ug, options)
        # end

        # twistopt_snopt = sum(xopt)
        # fopt_snopt = value(f!(zeros(ng), xopt))
        # fcalls_snopt = model_cache["fcalls"]

        # # Brute force answer
        # model_cache = Dict()
        # model_cache["fcalls"] = 0

        # t_brute = @elapsed begin
        #     twists = range(lx[1], ux[1], length=200)
        #     fs = [value(f2!(zeros(2), [twist, 0])) for twist in twists]
        # end

        # twistopt_brute_i = argmin(fs)
        # twistopt_brute = twists[twistopt_brute_i]
        # fopt_brute = fs[twistopt_brute_i]
        # fcalls_brute = model_cache["fcalls"]

        # if verbose
        #     println()
        #     @printf "%s%15.15s   %-10s %-10s %-10s %-10s\n" "\t"^(v_lvl+1) "" "Opt twist" "Max L/D" "Time" "Function calls"
        #     @printf "%s%15.15s   %-10.4f %-10.4f %-3.0f secs  %4.3g\n" "\t"^(v_lvl+1) "Brute force" twistopt_brute -fopt_brute t_brute fcalls_brute
        #     @printf "%s%15.15s   %-10.4f %-10.4f %-3.0f secs  %4.3g\n" "\t"^(v_lvl+1) "Ipopt" twistopt_snopt -fopt_snopt t_snopt fcalls_snopt
        # end

        # # @testset "Max L/D w.r.t. twist" begin
        # #     @test info == :Solved_To_Acceptable_Level || info == :Solve_Succeeded
        # #     @test isapprox(twistopt_snopt, twistopt_brute, atol=0.05)
        # #     @test isapprox(fopt_snopt, fopt_brute, atol=0.01)
        # # end


    end

end

nothing
