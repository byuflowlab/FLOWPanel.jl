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


# NOTE: Use CSDA with DFSane, or ForwardDiff with Broyden


@testset verbose=verbose "Lifting Line Optimization Tests w/ SNOPT" begin

    @testset "AOA vs twist optimization equivalence" begin

        # --------------- OPTIMIZATION TEST: Optimum AOA -----------------------
        if verbose
            println("\n"*"\t"^(v_lvl)*"Optimization test: Max L/D w.r.t. AOA")
        end

        model_cache = Dict()

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

            sweep = 20.0

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
                                    solver        = pnl.SimpleNonlinearSolve.SimpleDFSane(),         # Indifferent to initial guess, but somewhat not robust   <---- NOT COMPATIBLE WITH FORWARDDIFF
                                    # solver        = pnl.SimpleNonlinearSolve.SimpleTrustRegion(),    # Trust region needs a good initial guess, but it converges very reliably
                                    # solver          = pnl.NonlinearSolve.SimpleBroyden(),              # This seems to converge well while being compatible with ForwardDiff
                                    
                                    solver_optargs  = (; 
                                                        abstol = 1e-12,             # <-- tight tolerance to converge derivatives
                                                        maxiters = 800,
                                                        ),
                                    
                                    align_joints_with_Uinfs = !false,                # Whether to align joint bound vortices with the freestream
                                    
                                    use_Uind_for_force = true,                      # Whether to use Uind as opposed to selfUind for force postprocessing
                                                                                    # (`true` for more accurate spanwise cd distribution, but worse integrated CD)
                                    
                                    stability_derivatives = !true,                   # Whether to calculate stability derivatives
                
                                    cache
                                    )

            # Fetch results
            L = out.L
            D = out.D
            
            # Objective (minimization)
            obj = -L/D

            # Constraints (none)
            # g[1] = 

            return obj
        end


        x0 = [2.0]             # Initial guess

        lx = [-10.0]           # Lower bounds on x
        ux = [10.0]            # Upper bounds on x

        ng = 0                      # Number of constraints
        lg = -Inf*ones(ng)          # Lower bounds on g
        ug = zeros(ng)              # Upper bounds on g

        snopt_options = Dict(
                "Major iterations limit" => 50,
            )
        solver = SNOPT(options=snopt_options)

        # options = Options(; solver, derivatives=ForwardAD())
        options = Options(; solver, derivatives=ComplexStep())

        model_cache = Dict()

        t_snopt = @elapsed begin
            xopt, fopt, info = minimize(f!, deepcopy(x0), ng, lx, ux, lg, ug, options)
        end

        aoaopt_snopt = sum(xopt)
        fopt_snopt = value(f!(zeros(ng), xopt))
        fcalls_snopt = model_cache["fcalls"]

        # Brute force answer
        model_cache = Dict()

        t_brute = @elapsed begin
            aoas = range(lx[1], ux[1], length=200)
            fs = [value(f!(zeros(ng), [aoa])) for aoa in aoas]
        end

        aoaopt_brute_i = argmin(fs)
        aoaopt_brute = aoas[aoaopt_brute_i]
        fopt_brute = fs[aoaopt_brute_i]
        fcalls_brute = model_cache["fcalls"]

        fopt_aoa = fopt_snopt

        if verbose
            println()
            @printf "%s%15.15s   %-10s %-10s %-10s %-10s\n" "\t"^(v_lvl+1) "" "Opt AOA" "Max L/D" "Time" "Function calls"
            @printf "%s%15.15s   %-10.4f %-10.4f %-3.0f secs  %4.3g\n" "\t"^(v_lvl+1) "Brute force" aoaopt_brute -fopt_brute t_brute fcalls_brute
            @printf "%s%15.15s   %-10.4f %-10.4f %-3.0f secs  %4.3g\n" "\t"^(v_lvl+1) "SNOPT" aoaopt_snopt -fopt_snopt t_snopt fcalls_snopt
        end

        @testset "Max L/D w.r.t. alpha" begin
            @test info in [:Solved_To_Acceptable_Level, :Solve_Succeeded, "Finished successfully: optimality conditions satisfied"]
            @test isapprox(aoaopt_snopt, aoaopt_brute, atol=0.05)
            @test isapprox(fopt_snopt, fopt_brute, atol=0.01)
        end


        # --------------- TEST DERIVATIVES OVER TWIST --------------------------
        if verbose
            println("\n"*"\t"^(v_lvl)*"AD test over twist")
        end

        function f2!(args...; optargs...)
            return f!(args...; optargs..., optimize_aoa=false)
        end

        g = zeros(0)

        # Compute derivative through finite difference
        global daoa = 0.001
        global cache_m1 = Dict()
        global cache_diff = Dict()
        global cache_p1 = Dict()
        global f_m1 = f2!(g, x0 .- daoa; cache=cache_m1)
        global f_diff = f2!(g, x0; cache=cache_diff)
        global f_p1 = f2!(g, x0 .+ daoa; cache=cache_p1)

        global dfdx_diff = (f_p1 - f_m1) / (2*daoa)

        # Compute derivatives through complex-step
        global h = eps()^2
        global cache_csda = Dict()
        global f_csda = f2!(g, [Complex(x0[1], h)]; cache=cache_csda)
        global dfdx_csda = imag(f_csda) / h
        f_csda = real(f_csda)

        # Compute derivatives through Dual
        global cache_dual = Dict()
        global f_dual = f2!(g, [Dual(x0[1], Partials((1.0,)))]; cache=cache_dual)
        global dfdx_dual = partials(f_dual)[1]
        f_dual = value(f_dual)

        f_ref = f_diff
        dfdx_ref = dfdx_diff

        # Tests
        @testset "Derivative verification" begin
            for (lbl, (f, dfdx), (atol1, atol2)) in (
                                    ("Finite Difference", (f_diff, dfdx_diff), (1e-10, 1e-6)),
                                    ("CSDA", (f_csda, dfdx_csda), (1e-10, 5e-4)),
                                    ("Dual", (f_dual, dfdx_dual), (1e-10, 1e-6)),
                                )

                if verbose
                    println("\t"^(v_lvl+1)*"Comparing $(lbl)...")
                end

                @testset "$(lbl) f" begin @test isapprox(f, f_ref; atol=atol1) end
                @testset "$(lbl) dfdx" begin @test isapprox(dfdx, dfdx_ref; atol=atol2) end
                
            end
        end

        # --------------- OPTIMIZATION TEST: Optimum twist ---------------------
        if verbose
            println("\n"*"\t"^(v_lvl)*"Optimization test: Max L/D w.r.t. twist")
        end

        model_cache = Dict()

        t_snopt = @elapsed begin
            xopt, fopt, info = minimize(f2!, deepcopy(x0), ng, lx, ux, lg, ug, options)
        end

        twistopt_snopt = sum(xopt)
        fopt_snopt = value(f!(zeros(ng), xopt))
        fcalls_snopt = model_cache["fcalls"]

        # Brute force answer
        model_cache = Dict()

        t_brute = @elapsed begin
            twists = range(lx[1], ux[1], length=200)
            fs = [value(f2!(zeros(ng), [twist])) for twist in twists]
        end

        twistopt_brute_i = argmin(fs)
        twistopt_brute = twists[twistopt_brute_i]
        fopt_brute = fs[twistopt_brute_i]
        fcalls_brute = model_cache["fcalls"]

        fopt_twist = fopt_snopt

        if verbose
            println()
            @printf "%s%15.15s   %-10s %-10s %-10s %-10s\n" "\t"^(v_lvl+1) "" "Opt twist" "Max L/D" "Time" "Function calls"
            @printf "%s%15.15s   %-10.4f %-10.4f %-3.0f secs  %4.3g\n" "\t"^(v_lvl+1) "Brute force" twistopt_brute -fopt_brute t_brute fcalls_brute
            @printf "%s%15.15s   %-10.4f %-10.4f %-3.0f secs  %4.3g\n" "\t"^(v_lvl+1) "SNOPT" twistopt_snopt -fopt_snopt t_snopt fcalls_snopt
        end

        @testset "Max L/D w.r.t. twist" begin
            @test info in [:Solved_To_Acceptable_Level, :Solve_Succeeded, "Finished successfully: optimality conditions satisfied"]
            @test isapprox(twistopt_snopt, twistopt_brute, atol=0.05)
            @test isapprox(fopt_snopt, fopt_brute, atol=0.05)
        end


        @testset "AOA and twist equivalence" begin
            @test isapprox(aoaopt_snopt, twistopt_snopt, atol=0.05)
            @test isapprox(fopt_aoa, fopt_twist, atol=1e-4)
        end


    end

end

nothing
