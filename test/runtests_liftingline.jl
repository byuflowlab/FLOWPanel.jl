#=##############################################################################
# DESCRIPTION
    Unit tests of lifting line solver
=###############################################################################

using Test
import Printf: @printf
import FLOWPanel as pnl
import FLOWPanel: mean, norm, dot, cross

import ForwardDiff: Dual, Partials, value, partials, jacobian

try
    verbose
catch
    global verbose = true
end
v_lvl = 0


@testset verbose=verbose "Lifting Line Tests" begin

    @testset "Swept wing validation" begin

        # --------------- SWEPT WING TESTS -----------------------------------------
        if verbose
            println("\n"*"\t"^(v_lvl)*"Swept-wing test of lifting line")
        end

        CLexp = 0.238
        CDexp = 0.005

        CLvsp = 2.572866962740000263e-01
        CDvsp = 3.355815780000000192e-03
        Cmvsp = -3.758078951200000128e-01

        airfoil_path  = joinpath(pnl.examples_path, "data")


        # ----------------- SIMULATION PARAMETERS --------------------------------------
        AOA             = 4.2                           # (deg) angle of attack
        magUinf         = 49.7                          # (m/s) freestream velocity
        Uinf            = magUinf*[cosd(AOA), 0, sind(AOA)] # Freestream

        rho             = 1.225                         # (kg/m^3) air density

        # ------------------ INPUT PARAMETERS ----------------------------------

        # High-level parameters
        b               = 98*0.0254                     # (m) wing span

        # Discretization parameters
        nelements       = 40                            # Number of stripwise elements per semi-span
        discretization  = [                             # Multi-discretization of wing (seg length, ndivs, expansion, central)
                            (1.00,  nelements, 1/10, false),
                        ]
        symmetric       = true                          # Whether the wing is symmetric

        # Chord distribution (nondim y-position, nondim chord)
        chord_distribution = [
        #   2*y/b   c/b
            0.0     0.2
            1.0     0.2
        ]

        # Twist distribution (nondim y-position, twist)
        twist_distribution = [
        #   2*y/b   twist (deg)
            0.0     0.0
            1.0     0.0
        ]

        # Sweep distribution (nondim y-position, sweep)
        sweep_distribution = [
        #   2*y/b   sweep (deg)
            0.0     45.0
            1.0     45.0
        ]

        # Dihedral distribution (nondim y-position, dihedral)
        dihedral_distribution = [
        #   2*y/b   sweep (deg)
            0.0     0.0
            1.0     0.0
        ]

        # Span-axis distribution: chordwise point about which the wing is twisted, swept, and dihedralized (nondim y-position, nondim chord-position)
        spanaxis_distribution = [
        #   2*y/b   x/c
            0.0     0.25
            1.0     0.25
        ]

        # Airfoil contour distribution (nondim y-position, polar, airfoil type)
        airfoil_distribution = [
        #    2*y/b  polar file            airfoil type
            (0.00, "rae101-Re1p7e6.csv",     pnl.SimpleAirfoil),
            (1.00, "rae101-Re1p7e6.csv",     pnl.SimpleAirfoil)
        ]

        # ------------------ SOLVER PARAMETERS -------------------------------------
        deltasb         = 1.0                           # Blending distance, deltasb = 2*dy/b

        # deltajoint    = 0.15                          # Joint distance, deltajoint = dx/c
        deltajoint      = 1.0
                                                        # Nonlinear solver
        solver          = pnl.SimpleNonlinearSolve.SimpleDFSane()
        solver_optargs  = (; abstol = 1e-9)

        Dhat            = Uinf/norm(Uinf)               # Drag direction
        Shat            = [0, 1, 0]                     # Span direction
        Lhat            = cross(Dhat, Shat)             # Lift direction

        X0              = [0.0 * chord_distribution[1, 2]*b, 0, 0] # Center about which to calculate moments
        lhat            = Dhat                          # Rolling direction
        mhat            = Shat                          # Pitching direction
        nhat            = Lhat                          # Yawing direction

        cref            = mean(chord_distribution[:, 2]*b) # Reference chord
        nondim          = 0.5*rho*magUinf^2*b*cref      # Normalization factor


        # ------------------ GENERATE LIFTING LINE ---------------------------------
        ll = pnl.LiftingLine(
                                        airfoil_distribution; 
                                        b, chord_distribution, twist_distribution,
                                        sweep_distribution, dihedral_distribution,
                                        spanaxis_distribution,
                                        discretization,
                                        symmetric,
                                        deltasb, deltajoint,
                                        plot_discretization = false,
                                        element_optargs = (; 
                                            path = airfoil_path, 
                                            plot_polars = false,
                                            verbose = false
                                            )
                                        )

        # Test unit vectors
        @testset "Unit vectors" begin
            for unit in (ll.tangents, ll.spans, ll.normals,
                            ll.swepttangents, ll.lines, ll.sweptnormals)

                @test maximum(abs.(norm.(eachcol(unit)) .- 1)) <= eps(10.0)

            end
        end


        # Test sweep calculation
        @testset "Sweep calculation" begin
            @test abs(-pnl.calc_sweep(ll, 1) - sweep_distribution[1, 2]) < 1e-6
            @test abs(pnl.calc_sweep(ll, ll.nelements) - sweep_distribution[end, 2]) < 1e-6
        end


        # Test precomputation of velocity geometric matrix
        ll.Gammas .= 1.0

        Us_lazy = repeat(Uinf, 1, ll.nelements)
        pnl.selfUind!(ll, Us_lazy; precomputed=false)

        Us_precomputed = repeat(Uinf, 1, ll.nelements)
        pnl.selfUind!(ll, Us_precomputed; precomputed=true)

        @testset "G matrix precomp" begin
            @test maximum(abs.(Us_lazy - Us_precomputed)) <= eps(5e2)
        end


        # ------------------ CALL LINEAR SOLVER ------------------------------------

        # First run for JIT compilation (so the later run will be faster after 
        # compilation)
        pnl.solve_linear(ll, Uinf)

        # Timed run
        t_linear = @elapsed pnl.solve_linear(ll, Uinf)


        # ------------------ POSTPROCESS LINEAR SOLUTION ---------------------------
        # Force per stripwise element using Kutta-Joukowski theorem
        Fkj = pnl.calcfield_Fkj(ll, rho; fieldname="Fkj")

        # Integrated force
        Ftot = pnl.calcfield_Ftot(ll; F_fieldname="Fkj")

        # Decompose integrated force into lift and drag
        LDS = pnl.calcfield_LDS(ll, Lhat, Dhat, Shat; F_fieldname="Fkj")

        L = LDS[:, 1]
        D = LDS[:, 2]

        # # Loading distribution (force per unit span) using Kutta-Joukowski theorem
        # fs = pnl.calcfield_fkj(ll, rho)

        # lds = pnl.decompose(fs, Lhat, Dhat)

        # l = lds[1, :]
        # d = lds[2, :]
        
        # Integrated moment
        Mtot = pnl.calcfield_Mtot(ll, X0, rho; F_fieldname="Fkj", cm_fieldname=nothing)
        
        # Moment decomposed into axes
        lmn = pnl.calcfield_lmn(ll, lhat, mhat, nhat)
        roll, pitch, yaw = collect(eachcol(lmn))

        # Force coefficients
        CL_linear = sign(dot(L, Lhat)) * norm(L) / nondim
        CD_linear = sign(dot(D, Dhat)) * norm(D) / nondim

        # cl = l / (nondim/b)
        # cd = d / (nondim/b)

        # Cl = sign(dot(roll, lhat)) * norm(roll) / (nondim*cref)
        Cm_linear = sign(dot(pitch, mhat)) * norm(pitch) / (nondim*cref)
        # Cn = sign(dot(yaw, nhat)) * norm(yaw) / (nondim*cref)

        # Error
        err_CL_linear = abs(CL_linear-CLexp)/CLexp
        err_CD_linear = abs(CD_linear-CDexp)/CDexp
        err_Cm_linear = abs(Cm_linear-Cmvsp)/abs(Cmvsp)

        # ------------------ CALL NONLINEAR SOLVER ---------------------------------
        Uinfs = repeat(Uinf, 1, ll.nelements)

        # First run for JIT compilation (so the later run will be faster after 
        # compilation)
        pnl.solve(ll, Uinfs; solver, solver_optargs)

        # Timed run
        t_nonlinear = @elapsed begin 
            result, solver_cache = pnl.solve(ll, Uinfs; solver, solver_optargs)
        end

        # Check solver success
        # pnl.SimpleNonlinearSolve.SciMLBase.successful_retcode(result)

        # ------------------ POSTPROCESSING NONLINEAR SOLUTION ---------------------

        # Calculate force and moment coefficients
        calcs = pnl.calc_forcemoment_coefficients(ll, Uinfs, Uinf, 
                                                    rho, cref, b;
                                                    X0)
                                                    
        # Unpack calculations
        (; CD, CY, CL) = calcs
        (; Cl, Cm, Cn) = calcs
        (; Dhat, Shat, Lhat) = calcs
        (; lhat, mhat, nhat) = calcs
        (; q, Aref, bref, cref) = calcs


        # Coefficients
        CL_nonlinear = CL
        CD_nonlinear = CD
        Cm_nonlinear = Cm

        # Error
        err_CL_nonlinear = abs(CL_nonlinear-CLexp)/CLexp
        err_CD_nonlinear = abs(CD_nonlinear-CDexp)/CDexp
        err_Cm_nonlinear = abs(Cm_nonlinear-Cmvsp)/abs(Cmvsp)


        if verbose
            println()
            @printf "%s%15.15s %-7s %-7s %-7s %-10s %-10s %-10s %-10s\n" "\t"^(v_lvl+1) "" "CL" "CD" "Cm" "Time" "CL error" "CD error" "Cm error"
            @printf "%s%15.15s %-7s %-7s %-7s %-10s %-10s %-10s %-10s\n" "\t"^(v_lvl+1) "Experimental" CLexp CDexp "-" "-" "-" "-" "-"
            @printf "%s%15.15s %-7.4f %-7.4f %-7.4f %-10s %-10s %-10s %-10s\n" "\t"^(v_lvl+1) "OpenVSP" CLvsp CDvsp Cmvsp "-" "-" "-" "-"
            @printf "%s%15.15s %-7.4f %-7.4f %-7.4f %-3.0f msecs  %4.3g﹪     %4.3g﹪     %4.3g﹪\n" "\t"^(v_lvl+1) "Linear LL" CL_linear CD_linear Cm_linear t_linear*1000 err_CL_linear*100 err_CD_linear*100 err_Cm_linear*100
            @printf "%s%15.15s %-7.4f %-7.4f %-7.4f %-3.0f msecs  %4.3g﹪     %4.3g﹪     %4.3g﹪\n" "\t"^(v_lvl+1) "Nonlinear LL" CL_nonlinear CD_nonlinear Cm_nonlinear t_nonlinear*1000 err_CL_nonlinear*100 err_CD_nonlinear*100 err_Cm_nonlinear*100
        end

        # Test result
        @testset "Weber validation" begin
            @test err_CL_linear <= 0.03
            @test err_CD_linear <= 0.50
            @test err_Cm_linear <= 0.12
            @test err_CL_nonlinear <= 0.04
            @test err_CD_nonlinear <= 0.11
            @test err_Cm_nonlinear <= 0.06
        end


        # ------------------ TEST STABILITY DERIVATIVES ----------------------------

        # NOTE: Here we define a sideslip angle in addition to the angle of attack

        alpha           = 4.2                           # (deg) angle of attack
        beta            = 2.0                           # (deg) sideslip angle

        # NOTE: In the dual numbers we will implicitely defined the first partial to be 
        #       the derivative w.r.t. angle of attack and the second partial to be
        #       the derivative w.r.t. sideslip angle

        alpha           = Dual(alpha, Partials((1.0, 0.0))) # Convert angle of attack into dual number for automatic differentiation
        beta            = Dual(beta,  Partials((0.0, 1.0))) # Convert sideslip angle into dual number for automatic differentiation
        NumType         = typeof(alpha)                 # Number type for LiftingLine

        Uinf            = magUinf*pnl.direction(; alpha, beta) # Freestream vector

        # NOTE: We need to manually define the direction of lift, drag, and span, as 
        #       well as roll, pitching, and yawing moment since we now have a sideslip
        #       angle

        Dhat            = pnl.direction(; alpha)        # Drag direction
        Shat            = [0, 1, 0]                     # Span direction
        Lhat            = cross(Dhat, Shat)             # Lift direction

        lhat            = Dhat                          # Rolling direction
        mhat            = Shat                          # Pitching direction
        nhat            = Lhat                          # Yawing direction

        # Add dihedral so the sideslip force is non-zero
        dihedral_distribution = [
        #   2*y/b   dihedral (deg)
            0.0     5.0
            1.0     5.0
        ]

        # Redefine lifting line with Dual numbers
        ll = pnl.LiftingLine{NumType}(
                                        airfoil_distribution; 
                                        b, chord_distribution, twist_distribution,
                                        sweep_distribution, dihedral_distribution,
                                        spanaxis_distribution,
                                        discretization,
                                        symmetric,
                                        deltasb, deltajoint,
                                        plot_discretization = false,
                                        element_optargs = (; 
                                            path = airfoil_path, 
                                            plot_polars = false,
                                            verbose = false
                                            )
                                        )

        # Freestream velocity at each stripwise element
        Uinfs = repeat(Uinf, 1, ll.nelements)

        # Run solver
        result, solver_cache = pnl.solve(ll, Uinfs; 
                                            debug=true,             # `true` returns the residual rms
                                            aoas_initial_guess=alpha, 
                                            solver, solver_optargs)


        # Calculate force and moment coefficients
        calcs = pnl.calc_forcemoment_coefficients(ll, Uinfs, Uinf, 
                                                    rho, cref, b;
                                                    X0, 
                                                    Dhat, Shat, Lhat,
                                                    lhat, mhat, nhat)

        # Unpack calculations
        (; CD, CY, CL) = calcs                      # Drag, side, and lift forces
        (; Cl, Cm, Cn) = calcs                      # Roll, pitch, and yawing moment
        (; q, Aref, bref, cref) = calcs             # Reference dynamic pressure, area, span, and chord

        dCLdα = partials(CL)[1]
        dCDdα = partials(CD)[1]
        dCmdα = partials(Cm)[1]

        dCLdβ = partials(CL)[2]
        dCDdβ = partials(CD)[2]
        dCmdβ = partials(Cm)[2]

        CL = value(CL)
        CD = value(CD)
        Cm = value(Cm)

        #Tests
        @testset "Stability derivatives" begin
            @test abs(CL - 0.24709439885542692) <= 1e-12
            @test abs(CD - 0.005540765716982027) <= 1e-12
            @test abs(Cm - -0.35783803349227894) <= 1e-12
            @test abs(dCLdα - 0.06773567828194929) <= 1e-12
            @test abs(dCDdα - 0.001944238708320056) <= 1e-12
            @test abs(dCmdα - -0.09941038196468542) <= 1e-12
            @test abs(dCLdβ - 0.0002089112131300865) <= 1e-12
            @test abs(dCDdβ - -6.872400820022926e-6) <= 1e-12
            @test abs(dCmdβ - -0.00033435257221903946) <= 1e-12
        end

    end



    @testset "run_liftingline compatibility with AD" begin

        # --------------- run_liftingline TESTS --------------------------------
        if verbose
            println("\n"*"\t"^(v_lvl)*"run_liftingline AD compatibility tests")
        end

        alpha = 4.2
        beta = 1.0

        X0 = zeros(3)

        dihedral_distribution = [
        #   2*y/b   dihedral (deg)
            0.0     4.0
            1.0     4.0
        ]

        optargs = (; alpha, beta, X0, dihedral_distribution, verbose=false)

        # Wall-clock benchmark for reference
        out = pnl.run_liftingline(; optargs...)
        cache = out.cache
        pnl.run_liftingline(; cache, optargs...)

        t = @elapsed pnl.run_liftingline(; cache, optargs...)

        if verbose
            println("\n"*"\t"^(v_lvl+1)*"run_liftingline run time: $(Int(round(t*1000))) ms")
        end
        

        # Reference values (these where obtained through built-in Duals)
        out_ref = (;    CD = 0.0055338499113409614, CY = 0.0002269891577923591, CL = 0.24681427783627455, 
                        Cl = -0.009055479511985452, Cm = -0.3571370514170437, Cn = 0.0003528888094935135, 
                        dCDdα = 0.0019359109713538113, dCYdα = -4.73088519324675e-6, dCLdα = 0.06765598314629265, 
                        dCldα = -0.00136491717526742, dCmdα = -0.09912479914111964, dCndα = 6.05371898560579e-5, 
                        dCDdβ = -2.8762128067772757e-6, dCYdβ = 0.00022690799615432357, dCLdβ = 6.364513730187299e-5, 
                        dCldβ = -0.009053779138291685, dCmdβ = -0.00010373311929720024, dCndβ = 0.00035263241876742245)

        # Compute stability derivatives through finite difference
        dalpha = 0.001
        out_m1 = pnl.run_liftingline(; optargs..., alpha=alpha-dalpha, stability_derivatives=false)
        out_1  = pnl.run_liftingline(; optargs..., alpha=alpha,        stability_derivatives=false)
        out_p1 = pnl.run_liftingline(; optargs..., alpha=alpha+dalpha, stability_derivatives=false)

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

        # Built-in Dual computation of stability derivatives
        out_dual = pnl.run_liftingline(; optargs..., alpha, stability_derivatives=true)

        # ForwardDiff computation of stability derivatives
        function f(x; stability_derivatives=false)
            out = pnl.run_liftingline(; optargs..., alpha=x[1], beta=x[2], stability_derivatives)

            (; CD, CY, CL, Cl, Cm, Cn,
            dCDdα, dCYdα, dCLdα, dCldα, dCmdα, dCndα,
            dCDdβ, dCYdβ, dCLdβ, dCldβ, dCmdβ, dCndβ) = out

            return [CD, CY, CL, Cl, Cm, Cn,
            dCDdα, dCYdα, dCLdα, dCldα, dCmdα, dCndα,
            dCDdβ, dCYdβ, dCLdβ, dCldβ, dCmdβ, dCndβ]
        end

        f1(x) = f(x; stability_derivatives=false)

        out_fd1_p = f1([alpha, beta])
        out_fd1_j = jacobian(f1, [alpha, beta])

        out_fd1 = (;    CD=out_fd1_p[1], CY=out_fd1_p[2], CL=out_fd1_p[3], Cl=out_fd1_p[4], Cm=out_fd1_p[5], Cn=out_fd1_p[6],
                        dCDdα=out_fd1_j[1, 1], dCYdα=out_fd1_j[2, 1], dCLdα=out_fd1_j[3, 1], dCldα=out_fd1_j[4, 1], dCmdα=out_fd1_j[5, 1], dCndα=out_fd1_j[6, 1],
                        dCDdβ=out_fd1_j[1, 2], dCYdβ=out_fd1_j[2, 2], dCLdβ=out_fd1_j[3, 2], dCldβ=out_fd1_j[4, 2], dCmdβ=out_fd1_j[5, 2], dCndβ=out_fd1_j[6, 2])

        # ForwardDiff computation of stability derivatives embedded in the built-in Dual computation
        f2(x) = f(x; stability_derivatives=true)

        out_fd2_p = f2([alpha, beta])
        out_fd2_j = jacobian(f2, [alpha, beta])

        out_fd2 = (;    CD=out_fd2_p[1], CY=out_fd2_p[2], CL=out_fd2_p[3], Cl=out_fd2_p[4], Cm=out_fd2_p[5], Cn=out_fd2_p[6],
                        dCDdα=out_fd2_j[1, 1], dCYdα=out_fd2_j[2, 1], dCLdα=out_fd2_j[3, 1], dCldα=out_fd2_j[4, 1], dCmdα=out_fd2_j[5, 1], dCndα=out_fd2_j[6, 1],
                        dCDdβ=out_fd2_j[1, 2], dCYdβ=out_fd2_j[2, 2], dCLdβ=out_fd2_j[3, 2], dCldβ=out_fd2_j[4, 2], dCmdβ=out_fd2_j[5, 2], dCndβ=out_fd2_j[6, 2])

        # # Print dual result for reference
        # (; CD, CY, CL, Cl, Cm, Cn,
        #     dCDdα, dCYdα, dCLdα, dCldα, dCmdα, dCndα,
        #     dCDdβ, dCYdβ, dCLdβ, dCldβ, dCmdβ, dCndβ)  = out_dual

        # @show (; CD, CY, CL, Cl, Cm, Cn,
        #     dCDdα, dCYdα, dCLdα, dCldα, dCmdα, dCndα,
        #     dCDdβ, dCYdβ, dCLdβ, dCldβ, dCmdβ, dCndβ) 
    
        for (lbl, out, atol) in (
                                ("Finite Difference", out_diff, 1e-6),
                                ("Built-in Duals", out_dual, eps(10.0)),
                                ("ForwardDiff w/o built-in Duals", out_fd1, 1e-13),
                                ("ForwardDiff with built-in Duals", out_fd2, eps(10.0)),
                            )
            
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


end

nothing
