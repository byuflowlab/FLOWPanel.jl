#=##############################################################################
# DESCRIPTION
    Unit tests of lifting line solver
=###############################################################################

using Test
import Printf: @printf
import FLOWPanel as pnl
import FLOWPanel: mean, norm, dot, cross

import ForwardDiff: Dual, Partials, value, partials, jacobian
import SNOW: minimize, Options, IPOPT, ForwardAD, ComplexStep

import ForwardDiff

try
    verbose
catch
    global verbose = true
end
v_lvl = 0


@testset verbose=verbose "Lifting Line Optimization Tests" begin

    @testset "run_liftingline compatibility with AD" begin

        # --------------- run_liftingline TESTS --------------------------------
        if verbose
            println("\n"*"\t"^(v_lvl)*"run_liftingline AD compatibility tests")
        end

        alpha = 4.2
        beta = 1.0

        X0 = [-0.2 * 98*0.0254 * 0.25, 0, 0]

        sweep_distribution = [
        #   2*y/b   sweep (deg)
            0.0     10.0
            1.0     10.0
        ]

        dihedral_distribution = [
        #   2*y/b   dihedral (deg)
            0.0     4.0
            1.0     4.0
        ]

        optargs = (; alpha, beta, X0, 
                    sweep_distribution, dihedral_distribution, verbose=false)

        # Wall-clock benchmark for reference
        out = pnl.run_liftingline(; optargs...)
        cache = out.cache
        pnl.run_liftingline(; cache, optargs...)

        t = @elapsed pnl.run_liftingline(; cache, optargs...)

        if verbose
            println("\n"*"\t"^(v_lvl+1)*"run_liftingline run time: $(Int(round(t*1000))) ms\n")
        end
        

        # Compute stability derivatives through finite difference
        global dalpha = 0.001
        global out_m1 = pnl.run_liftingline(; optargs..., alpha=alpha-dalpha, stability_derivatives=false)
        global out_1  = pnl.run_liftingline(; optargs..., alpha=alpha,        stability_derivatives=false)
        global out_p1 = pnl.run_liftingline(; optargs..., alpha=alpha+dalpha, stability_derivatives=false)

        dbeta = 0.001
        out_m2 = pnl.run_liftingline(; optargs..., alpha, beta=beta-dbeta,  stability_derivatives=false)
        out_2  = pnl.run_liftingline(; optargs..., alpha, beta=beta,        stability_derivatives=false)
        out_p2 = pnl.run_liftingline(; optargs..., alpha, beta=beta+dbeta,  stability_derivatives=false)

        out_diff = (;   CD=out_1.CD, CY=out_1.CY, CL=out_1.CL,
                        Cl=out_1.Cl, Cm=out_1.Cm, Cn=out_1.Cn,
                        dCDdα=(out_p1.CD-out_m1.CD)/(2*dalpha), dCYdα=(out_p1.CY-out_m1.CY)/(2*dalpha), dCLdα=(out_p1.CL-out_m1.CL)/(2*dalpha),
                        dCldα=(out_p1.Cl-out_m1.Cl)/(2*dalpha), dCmdα=(out_p1.Cm-out_m1.Cm)/(2*dalpha), dCndα=(out_p1.Cn-out_m1.Cn)/(2*dalpha),
                        dCDdβ=(out_p2.CD-out_m2.CD)/(2*dbeta), dCYdβ=(out_p2.CY-out_m2.CY)/(2*dbeta), dCLdβ=(out_p2.CL-out_m2.CL)/(2*dbeta),
                        dCldβ=(out_p2.Cl-out_m2.Cl)/(2*dbeta), dCmdβ=(out_p2.Cm-out_m2.Cm)/(2*dbeta), dCndβ=(out_p2.Cn-out_m2.Cn)/(2*dbeta),
                    )

        # Compute stability derivatives through complex-step
        global h = eps()^2
        global out_csdaα = pnl.run_liftingline(; optargs..., alpha=Complex(alpha, h), stability_derivatives=false)
        out_csdaβ = pnl.run_liftingline(; optargs...,   beta=Complex(beta, h), stability_derivatives=false)

        out_csda = (;   CD=real(out_csdaα.CD), CY=real(out_csdaα.CY), CL=real(out_csdaα.CL),
                        Cl=real(out_csdaα.Cl), Cm=real(out_csdaα.Cm), Cn=real(out_csdaα.Cn),
                        dCDdα=imag(out_csdaα.CD)/h, dCYdα=imag(out_csdaα.CY)/h, dCLdα=imag(out_csdaα.CL)/h,
                        dCldα=imag(out_csdaα.Cl)/h, dCmdα=imag(out_csdaα.Cm)/h, dCndα=imag(out_csdaα.Cn)/h,
                        dCDdβ=imag(out_csdaβ.CD)/h, dCYdβ=imag(out_csdaβ.CY)/h, dCLdβ=imag(out_csdaβ.CL)/h,
                        dCldβ=imag(out_csdaβ.Cl)/h, dCmdβ=imag(out_csdaβ.Cm)/h, dCndβ=imag(out_csdaβ.Cn)/h,
                    )

        # Built-in Dual computation of stability derivatives
        global out_dual = pnl.run_liftingline(; optargs..., alpha, stability_derivatives=true)

        # ForwardDiff computation of stability derivatives
        function f(x; stability_derivatives=false, cache=Dict())
            out = pnl.run_liftingline(; optargs..., alpha=x[1], beta=x[2], 
                                            stability_derivatives, cache)

            (; CD, CY, CL, Cl, Cm, Cn,
            dCDdα, dCYdα, dCLdα, dCldα, dCmdα, dCndα,
            dCDdβ, dCYdβ, dCLdβ, dCldβ, dCmdβ, dCndβ) = out

            return [CD, CY, CL, Cl, Cm, Cn,
            dCDdα, dCYdα, dCLdα, dCldα, dCmdα, dCndα,
            dCDdβ, dCYdβ, dCLdβ, dCldβ, dCmdβ, dCndβ]
        end

        global cache_f1 = Dict()
        f1(x) = f(x; stability_derivatives=false, cache=cache_f1)

        out_fd1_p = f1([alpha, beta])
        out_fd1_j = jacobian(f1, [alpha, beta])

        out_fd1 = (;    CD=out_fd1_p[1], CY=out_fd1_p[2], CL=out_fd1_p[3], Cl=out_fd1_p[4], Cm=out_fd1_p[5], Cn=out_fd1_p[6],
                        dCDdα=out_fd1_j[1, 1], dCYdα=out_fd1_j[2, 1], dCLdα=out_fd1_j[3, 1], dCldα=out_fd1_j[4, 1], dCmdα=out_fd1_j[5, 1], dCndα=out_fd1_j[6, 1],
                        dCDdβ=out_fd1_j[1, 2], dCYdβ=out_fd1_j[2, 2], dCLdβ=out_fd1_j[3, 2], dCldβ=out_fd1_j[4, 2], dCmdβ=out_fd1_j[5, 2], dCndβ=out_fd1_j[6, 2])

        # ForwardDiff computation of stability derivatives embedded in the built-in Dual computation
        global cache_f2 = Dict()
        f2(x) = f(x; stability_derivatives=true, cache=cache_f2)

        out_fd2_p = f2([alpha, beta])
        out_fd2_j = jacobian(f2, [alpha, beta])

        out_fd2 = (;    CD=out_fd2_p[1], CY=out_fd2_p[2], CL=out_fd2_p[3], Cl=out_fd2_p[4], Cm=out_fd2_p[5], Cn=out_fd2_p[6],
                        dCDdα=out_fd2_j[1, 1], dCYdα=out_fd2_j[2, 1], dCLdα=out_fd2_j[3, 1], dCldα=out_fd2_j[4, 1], dCmdα=out_fd2_j[5, 1], dCndα=out_fd2_j[6, 1],
                        dCDdβ=out_fd2_j[1, 2], dCYdβ=out_fd2_j[2, 2], dCLdβ=out_fd2_j[3, 2], dCldβ=out_fd2_j[4, 2], dCmdβ=out_fd2_j[5, 2], dCndβ=out_fd2_j[6, 2])

        # Reference values
        # out_ref = out_diff
        out_ref = out_csda
    
        # Tests
        for (lbl, out, atol) in (
                                ("Finite Difference", out_diff, 5e-6),
                                ("CSDA", out_csda, 1e-12),
                                ("Built-in Duals", out_dual, 5e-6),
                                ("ForwardDiff w/o built-in Duals", out_fd1, 5e-6),
                                ("ForwardDiff with built-in Duals", out_fd2, 5e-6),
                            )

            
            if verbose
                println("\t"^(v_lvl+1)*"Comparing $(lbl)...")
            end
            
            @testset "$lbl" begin
                for sym in (:CD, :CY, :CL, :Cl, :Cm, :Cn,
                            :dCDdα, :dCYdα, :dCLdα, :dCldα, :dCmdα, :dCndα,
                            :dCDdβ, :dCYdβ, :dCLdβ, :dCldβ, :dCmdβ, :dCndβ)
        
                    @testset "$sym" begin
                        @test isapprox(getproperty(out, sym), getproperty(out_ref, sym); atol)
                    end
        
                end
            end
            
        end
        
    end


    @testset "AOA optimization: Max L/D" begin

        # --------------- OPTIMIZATION TEST: Optimum AOA -----------------------
        if verbose
            println("\n"*"\t"^(v_lvl)*"Optimization test: Max L/D w.r.t. AOA")
        end

        model_cache = Dict()
        model_cache["fcalls"] = 0

        function f!(g, x; optimize_aoa=true, cache=model_cache)

            # Fetch design variables
            aoa = x[1] + x[2]    # We defined two variables, otherwise SNOW won't work

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

            sweep = 45.0

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
                                    solver          = pnl.optimization_solver,
                                    
                                    solver_optargs  = (; 
                                                        abstol = 1e-9,  
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
            g[1] = -x[1] - 999
            g[2] = -x[2] - 999

            return obj
        end

        x0 = [4.2; 0.0]             # Initial guess

        lx = [-10.0, 0.0]         # Lower bounds on x
        ux = [10.0, 0.0]           # Upper bounds on x

        ng = 2                      # Number of constraints
        lg = -Inf*ones(ng)          # Lower bounds on g
        ug = zeros(ng)              # Upper bounds on g

        ipopt_options = Dict(
                "max_iter" => 50,
                "nlp_scaling_method" => "none",
            )
        solver = IPOPT(ipopt_options)
        options = Options(; solver, derivatives=ForwardAD())
        # options = Options(; solver, derivatives=ComplexStep())

        model_cache = Dict()
        model_cache["fcalls"] = 0

        t_ipopt = @elapsed begin
            xopt, fopt, info = minimize(f!, x0, ng, lx, ux, lg, ug, options)
        end

        aoaopt_ipopt = sum(xopt)
        fopt_ipopt = value(f!(zeros(ng), xopt))
        fcalls_ipopt = model_cache["fcalls"]

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
            @printf "%s%15.15s   %-10.4f %-10.4f %-3.0f secs  %4.3g\n" "\t"^(v_lvl+1) "Ipopt" aoaopt_ipopt -fopt_ipopt t_ipopt fcalls_ipopt
        end

        @testset "Max L/D w.r.t. alpha" begin
            @test info == :Solved_To_Acceptable_Level || info == :Solve_Succeeded
            @test isapprox(aoaopt_ipopt, aoaopt_brute, atol=0.05)
            @test isapprox(fopt_ipopt, fopt_brute, atol=0.02)
        end


    end

end

nothing
