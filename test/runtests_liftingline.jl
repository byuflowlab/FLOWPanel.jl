#=##############################################################################
# DESCRIPTION
    Unit tests of lifting line solver
=###############################################################################

using Test
import Printf: @printf
import FLOWPanel as pnl
import FLOWPanel: mean, norm, dot, cross

try
    verbose
catch
    global verbose = true
end
v_lvl = 0


@testset verbose=verbose "Lifting Line Tests" begin

    # --------------- SWEPT WING TESTS -----------------------------------------
    if verbose
        println("\n"*"\t"^(v_lvl)*"Swept-wing test of lifting line")
    end

    CLexp = 0.238
    CDexp = 0.005

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
        0.0     1/5.0
        1.0     1/5.0
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


    # ------------------ GENERATE LIFTING LINE ---------------------------
    liftingline = pnl.LiftingLine{Float64}(
                                    airfoil_distribution; 
                                    b, chord_distribution, twist_distribution,
                                    sweep_distribution, dihedral_distribution,
                                    spanaxis_distribution,
                                    discretization,
                                    symmetric,
                                    plot_discretization = false,
                                    element_optargs = (; 
                                        path = airfoil_path, 
                                        plot_polars = false
                                        )
                                    )



    # ------------------ CALL SOLVER ----------------------------------
    t = @elapsed pnl.solve_linear(liftingline, Uinf)


    # ------------------ POSTPROCESSING -------------------------------
    # Force per stripwise element using Kutta-Joukowski theorem
    Fkj = pnl.calcfield_Fkj(liftingline, rho; fieldname="Fkj")

    # Integrated force
    Ftot = pnl.calcfield_Ftot(liftingline; F_fieldname="Fkj")

    # Integrated force decomposed into lift and drag
    Dhat = Uinf/norm(Uinf)    # Drag direction
    Shat = [0, 1, 0]              # Span direction
    Lhat = cross(Dhat, Shat)  # Lift direction

    LDS = pnl.calcfield_LDS(liftingline, Lhat, Dhat, Shat; F_fieldname="Fkj")

    L = LDS[:, 1]
    D = LDS[:, 2]

    # # Loading distribution (force per unit span) using Kutta-Joukowski theorem
    # fs = pnl.calcfield_fkj(liftingline, rho)

    # lds = pnl.decompose(fs, Lhat, Dhat)

    # l = lds[1, :]
    # d = lds[2, :]

    # Force coefficients
    cref = mean(chord_distribution[:, 2]*b)
    nondim = 0.5*rho*magUinf^2*b*cref   # Normalization factor

    CL = sign(dot(L, Lhat)) * norm(L) / nondim
    CD = sign(dot(D, Dhat)) * norm(D) / nondim
    err = abs(CL-CLexp)/CLexp

    # cl = l / (nondim/b)
    # cd = d / (nondim/b)


    if verbose
        println()
        @printf "%s%15.15s %-7s %-7s %-12s\t%-10s\n" "\t"^(v_lvl+1) "" "CL" "CD" "Time" "CL error"
        @printf "%s%15.15s %-7s %-7s %-12s\t%-10s\n" "\t"^(v_lvl+1) "Experimental" CLexp CDexp "-" "-"
        @printf "%s%15.15s %-7.4f %-7.4f %-3.0f msecs   \t%4.3g﹪\n" "\t"^(v_lvl+1) "Linear LL" CL CD t*1000 err*100
    end

    res = err <= 0.03

    # Test result
    @test res

end
