#=##############################################################################
# DESCRIPTION
    Unit tests of solvers
=###############################################################################

using Test
import Printf: @printf
import FLOWPanel as pnl

try
    verbose
catch
    global verbose = true
end
v_lvl = 0

solvers_to_test = (
                        ( "Backslash", (solver=pnl.solve_backslash!, solver_optargs=()) ),
                        ( "LUdiv", (solver=pnl.solve_ludiv!, solver_optargs=()) ),
                        # ( "GMRES tol=1e-8", (solver=pnl.solve_gmres!, solver_optargs=(atol=1e-8, rtol=1e-8,)) ),
                        # ( "GMRES tol=1e-2", (solver=pnl.solve_gmres!, solver_optargs=(atol=1e-2, rtol=1e-2,)) ),
                  )


@testset verbose=verbose "Solver Tests" begin

    # --------------- SWEPT WING TESTS -----------------------------------------
    if verbose
        println("\n"*"\t"^(v_lvl)*"Swept-wing solver tests")
    end

    CLexp = 0.238
    CDexp = 0.005

    airfoil_path  = joinpath(pnl.examples_path, "data")

    # ----------------- SIMULATION PARAMETERS
    AOA           = 4.2                            # (deg) angle of attack
    magVinf       = 30.0                           # (m/s) freestream velocity
    Vinf          = magVinf*[cos(AOA*pi/180), 0, sin(AOA*pi/180)] # Freestream

    rho           = 1.225                          # (kg/m^3) air density

    # ----------------- GEOMETRY DESCRIPTION
    b             = 98*0.0254                      # (m) span length
    ar            = 5.0                            # Aspect ratio b/c_tip
    tr            = 1.0                            # Taper ratio c_tip/c_root
    twist_root    = 0                              # (deg) twist at root
    twist_tip     = 0                              # (deg) twist at tip
    lambda        = 45                             # (deg) sweep
    gamma         = 0                              # (deg) dihedral
    airfoil       = "airfoil-rae101.csv"           # Airfoil contour file

    # span_n = 60
    span_n = 8
    span_NDIVS_l    = [(1.0, span_n, 10.0, false)]
    span_NDIVS_r    = [(1.0, span_n, 10.0, false)]

    # rfl_n = 10
    rfl_n = 8
    rfl_NDIVS = [ (0.25, rfl_n,   10.0, false),
                  (0.50, rfl_n,    1.0, true),
                  (0.25, rfl_n, 1/10.0, false)]


    # ----------------- GENERATE BODY
    # Generate body
    if verbose; println("\t"^(v_lvl+2)*"Generating body..."); end;

    bodytype = pnl.RigidWakeBody{pnl.VortexRing}

    bodyoptargs_l = (
                    CPoffset=0.000000005,
                    kerneloffset=1e-8,
                    kernelcutoff=1e-14,
                    characteristiclength=(args...)->b/ar
                  )

    # Same arguments but negative CPoffset since the normals are flipped
    bodyoptargs_r = (
                    CPoffset=-bodyoptargs_l.CPoffset,
                    kerneloffset=bodyoptargs_l.kerneloffset,
                    kernelcutoff=bodyoptargs_l.kernelcutoff,
                    characteristiclength=bodyoptargs_l.characteristiclength
                  )


    # Loft left side of the wing from left to right
    wing_left = pnl.simplewing(b, ar, tr, twist_root, twist_tip, lambda, gamma;
                                bodytype=bodytype, bodyoptargs=bodyoptargs_l,
                                airfoil_root=airfoil, airfoil_tip=airfoil,
                                airfoil_path=airfoil_path,
                                rfl_NDIVS=rfl_NDIVS,
                                delim=",",
                                span_NDIVS=span_NDIVS_l,
                                b_low=-1.0, b_up=0.0
                               )

    # Loft right side of the wing from right to left
    wing_right = pnl.simplewing(b, ar, tr, twist_root, twist_tip, lambda, gamma;
                                bodytype=bodytype, bodyoptargs=bodyoptargs_r,
                                airfoil_root=airfoil, airfoil_tip=airfoil,
                                airfoil_path=airfoil_path,
                                rfl_NDIVS=rfl_NDIVS,
                                delim=",",
                                span_NDIVS=span_NDIVS_r,
                                b_low=1.0, b_up=0.0,
                               )

    # Put both sides together to make a wing with symmetric discretization
    bodies = [wing_left, wing_right]
    names = ["L", "R"]

    body = pnl.MultiBody(bodies, names)

    if verbose; println("\t"^(v_lvl+2)*"Number of panels:\t$(body.ncells)"); end;

    # Freestream at every control point
    Uinfs = repeat(Vinf, 1, body.ncells)
    Das = repeat(Vinf/magVinf, 1, body.nsheddings)
    Dbs = repeat(Vinf/magVinf, 1, body.nsheddings)


    if verbose
        println()
        @printf "%s%15.15s %-7s %-7s %-12s\t%-10s\n" "\t"^(v_lvl+1) "Solver" "CL" "CD" "Time" "CL error"
        @printf "%s%15.15s %-7s %-7s %-12s\t%-10s\n" "\t"^(v_lvl+1) "Experimental" CLexp CDexp "-" "-"
    end

    for (lbl, (optargs)) in solvers_to_test
        @test begin

            # ----------------- CALL SOLVER
            # Solve body
            t = @elapsed pnl.solve(body, Uinfs, Das, Dbs; optargs...)

            # ----------------- POST PROCESSING ------------------------------------------------
            # Calculate velocity away from the body
            Us = pnl.calcfield_U(body, body; fieldname="Uoff",
                                    offset=0.02, characteristiclength=(args...)->b/ar)

            # Calculate pressure coeffiecient
            Cps = pnl.calcfield_Cp(body, magVinf; U_fieldname="Uoff")

            # Calculate the force of each panel
            Fs = pnl.calcfield_F(body, magVinf, rho)
            # Calculate total force of the vehicle decomposed as lfit, drag, and sideslip
            Dhat = Vinf/pnl.norm(Vinf)        # Drag direction
            Shat = [0, 1, 0]              # Span direction
            Lhat = pnl.cross(Dhat, Shat)      # Lift direction

            LDS = pnl.calcfield_LDS(body, Lhat, Dhat)

            L = LDS[:, 1]
            D = LDS[:, 2]

            # Force coefficients
            nondim = 0.5*rho*magVinf^2*b^2/ar   # Normalization factor
            CL = sign(pnl.dot(L, Lhat)) * pnl.norm(L) / nondim
            CD = sign(pnl.dot(D, Dhat)) * pnl.norm(D) / nondim
            err = abs(CL-CLexp)/CLexp

            if verbose
                @printf "%s%15.15s %-7.4f %-7.4f %4.2f seconds\t%4.3g﹪\n" "\t"^(v_lvl+1) lbl CL CD t err*100
            end

            res = err <= 0.6

            # Test result
            res
        end
    end

end
