#=##############################################################################
# DESCRIPTION
    Unit tests of other solvers
=###############################################################################

using Test
import Printf: @printf
import CSV
import DataFrames: DataFrame

import FLOWPanel as pnl

try
    verbose
catch
    global verbose = true
end
v_lvl = 0

solvers_to_test = Any[
                        ( "Backslash", (solver=pnl.solve_backslash!, solver_optargs=()) ),
                        ( "LUdiv", (solver=pnl.solve_ludiv!, solver_optargs=()) ),
                  ]

# Add GPU test cases if GPU hardware is available
# NOTE: Here we test with Float64 instead of Float32 since the later leads to
#       an error larger than the tolerance
try

    import Pkg; Pkg.add("CUDA")
    import CUDA

    if CUDA.functional()

        using CUDA
        push!(solvers_to_test, ( "LUdiv + GPU (CUDA 64)",
                                (solver=pnl.solve_ludiv!, GPUArray=CuArray{Float64}, solver_optargs=()) ) )

    end

catch e
    @warn "CUDA device not found: $e"
end

try

    import Pkg; Pkg.add("Metal")
    import Metal

    if Metal.functional()

        using Metal
        push!(solvers_to_test, ( "LUdiv + GPU (Metal 64)",
                                (solver=pnl.solve_ludiv!, GPUArray=MtlArray{Float32}, solver_optargs=()) ) )

    end
catch e
    @warn "Metal device not found: $e"
end

try

    import Pkg; Pkg.add("AMDGPU")
    import AMDGPU

    if AMDGPU.functional()

        using AMDGPU
        push!(solvers_to_test, ( "LUdiv + GPU (AMD 64)",
                                (solver=pnl.solve_ludiv!, GPUArray=ROCArray{Float64}, solver_optargs=()) ) )

    end
catch e
    @warn "AMD device not found: $e"
end



@testset verbose=verbose "Solver Tests 2" begin

    # --------------- SWEPT WING TESTS -----------------------------------------
    if verbose
        println("\n"*"\t"^(v_lvl)*"Duct test (least-squares solver + ∇μ scheme)")
    end

    Lref = pnl.norm([-30.49874883049934, 0.0, 348.60229430138884])
    Dref = pnl.norm([17.96575492217547, 0.0, 1.5717998873754226])

    airfoil_path  = joinpath(pnl.examples_path, "data")

    # ----------------- SIMULATION PARAMETERS
    AOA             = 5                         # (deg) angles of attack to evaluate
    magVinf         = 30.0                      # (m/s) freestream velocity
    rho             = 1.225                     # (kg/m^3) air density

    # ----------------- GEOMETRY DESCRIPTION
    # Read duct contour (table in figure 7.4 of Lewis 1991)
    filename        = joinpath(pnl.examples_path, "data", "naca662015.csv")
    contour         = CSV.read(filename, DataFrame)

    aspectratio     = 0.6                       # Duct trailing edge aspect ratio l/d
    d               = 2*0.835                   # (m) duct diameter

    # ----------------- SOLVER PARAMETERS
    # Discretization
    NDIVS_theta     = 60                        # Number of azimuthal panels
    n_rfl           = 10                        # This controls the number of chordwise panels

    NDIVS_rfl_up = [                            # Discretization of airfoil upper surface
                    (0.25, n_rfl,   10.0, false),
                    (0.50, n_rfl,    1.0, true),
                    (0.25, n_rfl,    0.1, false)]

    NDIVS_rfl_lo = NDIVS_rfl_up                 # Discretization of airfoil lower surface

    # Solver: Vortex-ring least-squares
    bodytype        = pnl.RigidWakeBody{pnl.VortexRing, 2} # Elements and wake model

    # ----------------- GENERATE BODY
    # Generate body
    if verbose; println("\t"^(v_lvl+2)*"Generating body..."); end;

    # Re-discretize the contour of the body of revolution according to NDIVS
    xs, ys = pnl.gt.rediscretize_airfoil(contour[:, 1], contour[:, 2],
                                            NDIVS_rfl_up, NDIVS_rfl_lo;
                                            verify_spline=false)

    # Make sure that the contour is closed
    ys[end] = ys[1]

    # Scale contour by duct length
    xs *= d*aspectratio
    ys *= d*aspectratio

    # Move contour to the radial position
    ys .+= d/2

    # Collect points that make the contour of the body of revolution
    points = hcat(xs, ys)

    # Generate body of revolution
    body = pnl.generate_revolution_liftbody(bodytype, points, NDIVS_theta;
                                            bodyoptargs = (
                                                            CPoffset=1e-14,
                                                            kerneloffset=1e-8,
                                                            kernelcutoff=1e-14,
                                                            characteristiclength=(args...)->d*aspectratio
                                                )
                                            )

    if verbose; println("\t"^(v_lvl+2)*"Number of panels:\t$(body.ncells)"); end;

    # Freestream vector
    Vinf = magVinf*[cos(AOA*pi/180), 0, sin(AOA*pi/180)]

    # Freestream at every control point
    Uinfs = repeat(Vinf, 1, body.ncells)

    # Unitary direction of semi-infinite vortex at points `a` and `b` of each
    # trailing edge panel
    Das = repeat(Vinf/magVinf, 1, body.nsheddings)
    Dbs = repeat(Vinf/magVinf, 1, body.nsheddings)

    if verbose
        println()
        @printf "%s%19.19s %-7s %-7s %-12s\t%-10s\t%-10s\n" "\t"^(v_lvl+1) "Solver" "Lift" "Drag" "Time" "L error" "D error"
        @printf "%s%19.19s %-7.1f %-7.1f %-12s\t%-10s\t%-10s\n" "\t"^(v_lvl+1) "Reference" Lref Dref "-" "-" "-"
    end

    for (lbl, (optargs)) in solvers_to_test
        @test begin

            # ----------------- CALL SOLVER
            # Solve body
            t = @elapsed pnl.solve(body, Uinfs, Das, Dbs; optargs...)

            # ----------------- POST PROCESSING ------------------------------------------------
            # Calculate surface velocity U on the body
            Us = pnl.calcfield_U(body, body)

            # Calculate surface velocity U_∇μ due to the gradient of the doublet strength
            # UDeltaGamma = pnl.calcfield_Ugradmu(body)
            UDeltaGamma = pnl.calcfield_Ugradmu(body; sharpTE=true, force_cellTE=false)

            # Add both velocities together
            pnl.addfields(body, "Ugradmu", "U")

            # Calculate pressure coefficient (based on U + U_∇μ)
            Cps = pnl.calcfield_Cp(body, magVinf)

            # Calculate the force of each panel (based on Cp)
            Fs = pnl.calcfield_F(body, magVinf, rho)

            # Calculate total force of the vehicle decomposed as lift, drag, and sideslip
            Dhat = Vinf/pnl.norm(Vinf)        # Drag direction
            Shat = [0, 1, 0]                  # Span direction
            Lhat = pnl.cross(Dhat, Shat)      # Lift direction

            LDS = pnl.calcfield_LDS(body, Lhat, Dhat)

            L = pnl.norm(LDS[:, 1])
            D = pnl.norm(LDS[:, 2])

            Lerr = abs(L-Lref)/Lref
            Derr = abs(D-Dref)/Dref

            if verbose
                @printf "%s%19.19s %-7.1f %-7.1f %4.2f seconds\t%4.3g﹪\t%4.3g﹪\n" "\t"^(v_lvl+1) lbl L D t Lerr*100 Derr*100
            end

            Lres = Lerr <= 0.005
            Dres = Derr <= 0.005

            # Test result
            Lres && Dres
        end
    end

end
