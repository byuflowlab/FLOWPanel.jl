#=##############################################################################
# DESCRIPTION
    Unit tests of the integration of Meshes.jl with OpenVSP and SolidWorks+Gmsh
=###############################################################################

using Test
import Printf: @printf
import GeoIO

import FLOWPanel as pnl
import FLOWPanel: norm, dot, cross

try
    verbose
catch
    global verbose = true
end
v_lvl = 0


@testset verbose=verbose "OpenVSP + Meshes.jl test" begin

    if verbose
        println("\n"*"\t"^(v_lvl)*"OpenVSP + Meshes.jl integration test")
    end

    # --------------- SWEPT WING TESTS -----------------------------------------
    if verbose
        println("\n"*"\t"^(v_lvl+1)*"Swept wing test on OpenVSP grid")
    end

    # ----------------- SIMULATION PARAMETERS ----------------------------------
    AOA             = 4.2                           # (deg) angle of attack
    magVinf         = 30.0                          # (m/s) freestream velocity
    Vinf            = magVinf*[cos(AOA*pi/180), 0, sin(AOA*pi/180)] # Freestream

    rho             = 1.225                         # (kg/m^3) air density

    # ----------------- GEOMETRY DESCRIPTION -----------------------------------

    # VTK files to read
    meshfile        = joinpath(pnl.examples_path, "data", "sweptwing-openvsp.msh")

    # Span directions used to orient each trailing edge
    spandir         = [0, 1, 0]

    # Whether to flip control points against the direction of normals
    flip            = true

    # Reference dimensions
    b               = 98*0.0254                     # (m) span length
    ar              = 5.0                           # Aspect ratio b/c_tip

    # ----------------- SOLVER SETTINGS ----------------------------------------

    # Solver: direct linear solver for open bodies
    # bodytype = pnl.RigidWakeBody{pnl.VortexRing}    # Elements and wake model

    # Solver: least-squares solver for watertight bodies
    bodytype        = pnl.RigidWakeBody{pnl.VortexRing, 2} # Elements and wake model


    # ----------------- GENERATE BODY ------------------------------------------
    # Read VTK as a Meshes object
    msh = GeoIO.load(meshfile)
    msh = msh.geometry

    # Wrap Meshes object into a Grid object from GeometricTools
    grid = pnl.gt.GridTriangleSurface(msh)

    # Define trailing edge
    X1 = [1.74244, -1.2446, 0.0]
    X2 = [0.49784, 0.0, 0.0]
    X3 = [1.74244, 1.2446, 0.0]

    trailingedge1 = stack([X1 + val*(X2-X1) for val in range(0, 1, length=10001)])
    trailingedge2 = stack([X2 + val*(X3-X2) for val in range(0, 1, length=10001)])
    trailingedge = hcat(trailingedge1, trailingedge2[:, 2:end])

    # Sort TE points from left to right
    trailingedge = sortslices(trailingedge; dims=2, by = X -> pnl.dot(X, spandir))

    # Estimate span length (used as a reference length)
    spantips = extrema(X -> pnl.dot(X, spandir), eachcol(trailingedge))
    span = spantips[2] - spantips[1]

    # Generate shedding matrix
    shedding = pnl.calc_shedding(grid, trailingedge; tolerance=0.0001*span)

    # Generate paneled body
    body = bodytype(grid, shedding; CPoffset=(-1)^flip * 1e-14)

    if verbose; println("\t"^(v_lvl+2)*"Number of panels:\t$(body.ncells)"); end;

    # ----------------- CALL SOLVER --------------------------------------------

    # Freestream at every control point
    Uinfs = repeat(Vinf, 1, body.ncells)

    # Unitary direction of semi-infinite vortex at points `a` and `b` of each
    # trailing edge panel
    Das = repeat(Vinf/magVinf, 1, body.nsheddings)
    Dbs = repeat(Vinf/magVinf, 1, body.nsheddings)

    pnl.solve(body, Uinfs, Das, Dbs)
    # pnl.solve(body, Uinfs, Das, Dbs; GPUArray=CUDA.CuArray{Float32})

    # ----------------- POST PROCESSING ----------------------------------------
    # Calculate surface velocity U on the body
    Us = pnl.calcfield_U(body, body)

    # Calculate surface velocity U_∇μ due to the gradient of the doublet strength
    UDeltaGamma = pnl.calcfield_Ugradmu(body)

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

    L = LDS[:, 1]
    D = LDS[:, 2]

    # Force coefficients
    nondim = 0.5*rho*magVinf^2*b^2/ar   # Normalization factor
    CL = sign(dot(L, Lhat)) * norm(L) / nondim
    CD = sign(dot(D, Dhat)) * norm(D) / nondim

    CLexp = 0.238
    CDexp = 0.005

    # Test lift and drag
    CLerr = (CLexp - CL) / CLexp
    CDerr = (CDexp - CD) / CDexp

    if verbose
        println()
        @printf "%s%19.19s\t%-7s %-7s %-10s\t%-10s\n" "\t"^(v_lvl+1) "" "CL" "CD" "CL error" "CD error"
        @printf "%s%19.19s\t%-7.3f %-7.4f %-10s\t%-10s\n" "\t"^(v_lvl+1) "Experimental" CLexp CDexp "-" "-"
        @printf "%s%19.19s\t%-7.3f %-7.4f %4.3g﹪\t\t%4.3g﹪\n" "\t"^(v_lvl+1) "OpenVSP Grid" CL CD CLerr*100 CDerr*100
    end

    @test CLerr <= 0.05 && CDerr <= 1.5

end
