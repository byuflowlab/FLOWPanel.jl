#=##############################################################################
# DESCRIPTION
    Unit tests of the integration of Meshes.jl
=###############################################################################

using Test
import Printf: @printf

import Pkg; Pkg.add("GeoIO")    # This is here as opposed than in Project.toml
import GeoIO                    # because GeoIO is only compatible with Julia
                                # v1.9+, and we don't want to make the whole
                                # package follow the same compatibility
                                # constraint only because of this unit test
import FLOWPanel as pnl

try
    verbose
catch
    global verbose = true
end
v_lvl = 0


@testset verbose=verbose "Meshes.jl test" begin

    if verbose
        println("\n"*"\t"^(v_lvl)*"Meshes.jl integration test")
    end

    # --------------- SWEPT WING TESTS -----------------------------------------
        if verbose
            println("\n"*"\t"^(v_lvl+1)*"Swept wing test")
        end

    # --------------------------------------------------------------------------
    # --------------- SOLUTION IN STRUCTURED GRID ------------------------------
    # --------------------------------------------------------------------------
    if verbose; println("\t"^(v_lvl+2)*"Solving on structured grid..."); end;

    airfoil_path    = joinpath(pnl.examples_path, "data") # Where to find airfoil contours

    # ----------------- SIMULATION PARAMETERS --------------------------------------
    AOA             = 4.2                           # (deg) angle of attack
    magVinf         = 30.0                          # (m/s) freestream velocity
    Vinf            = magVinf*[cos(AOA*pi/180), 0, sin(AOA*pi/180)] # Freestream

    rho             = 1.225                         # (kg/m^3) air density

    # ----------------- GEOMETRY DESCRIPTION ---------------------------------------
    b               = 98*0.0254                     # (m) span length
    ar              = 5.0                           # Aspect ratio b/c_tip
    tr              = 1.0                           # Taper ratio c_tip/c_root
    twist_root      = 0                             # (deg) twist at root
    twist_tip       = 0                             # (deg) twist at tip
    lambda          = 45                            # (deg) sweep
    gamma           = 0                             # (deg) dihedral
    airfoil         = "airfoil-rae101.csv"          # Airfoil contour file


    # ----- Chordwise discretization

    n_rfl           = 8                             # Control number of chordwise panels
    NDIVS_rfl = [ (0.25, n_rfl,   10.0, false),
                  (0.50, n_rfl,    1.0, true),
                  (0.25, n_rfl, 1/10.0, false)]

    # ----- Spanwise discretization
    n_span          = 15                            # Number of spanwise panels on each side of the wing
    NDIVS_span_l    = [(1.0, n_span, 10.0, false)]  # Discretization of left side
    NDIVS_span_r    = [(1.0, n_span, 10.0, false)]  # Discretization of right side

    # ----------------- GENERATE BODY
    # Generate body
    bodytype = pnl.RigidWakeBody{pnl.VortexRing}    # Elements and wake model

    # Arguments for lofting the left side of the wing
    bodyoptargs_l = (
                    CPoffset=1e-14,                 # Offset control points slightly in the positive normal direction
                    characteristiclength=(args...)->b/ar,   # Characheristic length for control point offset
                    kerneloffset=1e-8,              # Offset of kernel to avoid singularities
                    kernelcutoff=1e-14              # Cutoff of kernel to avoid singularities
                  )

    # Same arguments but negative CPoffset since the normals are flipped
    bodyoptargs_r = (
                    CPoffset=-bodyoptargs_l.CPoffset,
                    characteristiclength=bodyoptargs_l.characteristiclength,
                    kerneloffset=bodyoptargs_l.kerneloffset,
                    kernelcutoff=bodyoptargs_l.kernelcutoff
                  )

    # Loft left side of the wing from left to right
    wing_left = pnl.simplewing(b, ar, tr, twist_root, twist_tip, lambda, gamma;
                                bodytype=bodytype, bodyoptargs=bodyoptargs_l,
                                airfoil_root=airfoil, airfoil_tip=airfoil,
                                airfoil_path=airfoil_path,
                                rfl_NDIVS=NDIVS_rfl,
                                delim=",",
                                span_NDIVS=NDIVS_span_l,
                                b_low=-1.0, b_up=0.0
                               )

    # Loft right side of the wing from right to left
    wing_right = pnl.simplewing(b, ar, tr, twist_root, twist_tip, lambda, gamma;
                                bodytype=bodytype, bodyoptargs=bodyoptargs_r,
                                airfoil_root=airfoil, airfoil_tip=airfoil,
                                airfoil_path=airfoil_path,
                                rfl_NDIVS=NDIVS_rfl,
                                delim=",",
                                span_NDIVS=NDIVS_span_r,
                                b_low=1.0, b_up=0.0,
                               )

    # Put both sides together to make a wing with symmetric discretization
    bodies = [wing_left, wing_right]
    names = ["L", "R"]

    body = pnl.MultiBody(bodies, names)

    # ----------------- CALL SOLVER --------------------------------------------

    # Freestream at every control point
    Uinfs = repeat(Vinf, 1, body.ncells)

    # Unitary direction of semi-infinite vortex at points `a` and `b` of each
    # trailing edge panel
    Das = repeat(Vinf/magVinf, 1, body.nsheddings)
    Dbs = repeat(Vinf/magVinf, 1, body.nsheddings)

    # Solve body (panel strengths) giving `Uinfs` as boundary conditions and
    # `Das` and `Dbs` as trailing edge rigid wake direction
    pnl.solve(body, Uinfs, Das, Dbs)

    # ----------------- POST PROCESSING ----------------------------------------
    # Calculate surface velocity U on the body
    Us = pnl.calcfield_U(body, body)

    # Calculate surface velocity U_∇μ due to the gradient of the doublet strength
    UDeltaGamma = pnl.calcfield_Ugradmu(body)
    # UDeltaGamma = pnl.calcfield_Ugradmu_cell(body)
    # UDeltaGamma = pnl.calcfield_Ugradmu_node(body)

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


    # ----------------- OUTPUTS TO COMPARE -------------------------------------
    # Panel strengths
    strength_str = [deepcopy(wing.strength) for wing in bodies]

    # Lift and drag
    Lstr = pnl.norm(LDS[:, 1])
    Dstr = pnl.norm(LDS[:, 2])

    # List of every neighbor of every cell (converted to linear index)
    neighbors_str = [[[ pnl.gt.neighbor(b.grid, ni, ci; preserveEdge=true)
                                                                for ni in 1:3]
                                                                for ci in 1:b.ncells]
                                                                for b in body.bodies]

    lincs = [LinearIndices( Tuple(collect( 1:(d != 0 ? d : 1)   for d in b.grid._ndivscells )) )
                                                                for b in body.bodies]

    neighbors_str = [[[ prod(ncoor .!= 0) ? linc[ncoor...] : 0  for ncoor in neighbors]
                                                                for neighbors in cellneigh]
                                                                for (cellneigh, linc) in zip(neighbors_str, lincs)]

    # Surface velocity from ∇𝜇 scheme
    Ugradmu_str = deepcopy(collect.(pnl.get_field(body, "Ugradmu")["field_data"]))


    # --------------------------------------------------------------------------
    # --------------- SOLUTION IN UNSTRUCTURED MESH ----------------------------
    # --------------------------------------------------------------------------
    if verbose; println("\t"^(v_lvl+2)*"Solving on unstructured grid..."); end;


    # ----------------- GEOMETRY DESCRIPTION -----------------------------------
    # VTK files to read
    meshfile1        = joinpath(pnl.examples_path, "data", "sweptwing-left000.vtu")
    meshfile2        = joinpath(pnl.examples_path, "data", "sweptwing-right000.vtu")

    # Span directions used to orient each trailing edge
    spandir1         = [0, 1, 0]
    spandir2         = [0, -1, 0]

    # Whether to flip the normals of each mesh (the user must check that the
    # normals that get outputed with `save(body, ...; debug=true)` are pointing
    # out of the body)
    flip1            = false
    flip2            = true

    # ----------------- GENERATE BODY ------------------------------------------
    # Read VTK as a Meshes object
    msh1 = GeoIO.load(meshfile1)
    msh1 = msh1.geometry

    # Wrap Meshes object into a Grid object from GeometricTools
    grid1 = pnl.gt.GridTriangleSurface(msh1)

    # Define trailing edge
    X1 = [1.74244,-1.2446,1.38094e-07]
    X2 = [0.49784,-1.38178e-16,1.38094e-07]
    trailingedge = stack([X1 + val*(X2-X1) for val in range(0, 1, length=1000)])

    # Sort TE points from left to right
    trailingedge = sortslices(trailingedge; dims=2, by = X -> pnl.dot(X, spandir1))

    # Estimate span length (used as a reference length)
    spantips = extrema(X -> pnl.dot(X, spandir1), eachcol(trailingedge))
    span = spantips[2] - spantips[1]

    # Generate shedding matrix
    shedding = pnl.calc_shedding(grid1, trailingedge; tolerance=0.001*span)

    # Generate paneled body
    wing1 = bodytype(grid1, shedding; CPoffset=(-1)^flip1 * 1e-14)

    # Repeat the same process for the other side of the wing
    msh2 = GeoIO.load(meshfile2)
    msh2 = msh2.geometry
    grid2 = pnl.gt.GridTriangleSurface(msh2)

    X1 = [1.74244,1.2446,1.38094e-07]
    X2 = [0.49784,1.38178e-16,1.38094e-07]
    trailingedge = stack([X1 + val*(X2-X1) for val in range(0, 1, length=1000)])
    trailingedge = sortslices(trailingedge; dims=2, by = X -> pnl.dot(X, spandir2))

    spantips = extrema(X -> pnl.dot(X, spandir2), eachcol(trailingedge))
    span = spantips[2] - spantips[1]
    shedding = pnl.calc_shedding(grid2, trailingedge; tolerance=0.001*span)

    wing2 = bodytype(grid2, shedding; CPoffset=(-1)^flip2 * 1e-14)

    # Put both sides together to make a wing with symmetric discretization
    bodies = [wing1, wing2]
    names = ["1", "2"]

    body = pnl.MultiBody(bodies, names)

    # ----------------- CALL SOLVER --------------------------------------------
    pnl.solve(body, Uinfs, Das, Dbs)


    # ----------------- POST PROCESSING ----------------------------------------
    # Calculate surface velocity U on the body
    Us = pnl.calcfield_U(body, body)

    # Calculate surface velocity U_∇μ due to the gradient of the doublet strength
    UDeltaGamma = pnl.calcfield_Ugradmu(body)
    # UDeltaGamma = pnl.calcfield_Ugradmu_cell(body)
    # UDeltaGamma = pnl.calcfield_Ugradmu_node(body)

    # Add both velocities together
    pnl.addfields(body, "Ugradmu", "U")

    # Calculate pressure coefficient (based on U + U_∇μ)
    Cps = pnl.calcfield_Cp(body, magVinf)

    # Calculate the force of each panel (based on Cp)
    Fs = pnl.calcfield_F(body, magVinf, rho)

    # Calculate total force of the vehicle decomposed as lift, drag, and sideslip
    LDS = pnl.calcfield_LDS(body, Lhat, Dhat)


    # ----------------- OUTPUTS TO COMPARE -------------------------------------
    # Panel strengths
    strength_uns = [deepcopy(wing.strength) for wing in bodies]

    # Lift and drag
    Luns = pnl.norm(LDS[:, 1])
    Duns = pnl.norm(LDS[:, 2])

    # List of every neighbor of every cell
    neighbors_uns = [[[pnl.gt.neighbor(b.grid, ni, ci; preserveEdge=true)[1] for ni in 1:3]
                                                                             for ci in 1:b.ncells]
                                                                             for b in body.bodies]

    # Surface velocity from ∇𝜇 scheme
    Ugradmu_uns = deepcopy(collect.(pnl.get_field(body, "Ugradmu")["field_data"]))

    bodies_uns = [deepcopy(b) for b in bodies]

    # --------------------------------------------------------------------------
    # --------------- COMPARE SOLUTIONS ----------------------------------------
    # --------------------------------------------------------------------------

    # Test solution error
    for wi in 1:2

        str = strength_str[wi]
        uns = strength_uns[wi]

        if verbose
            println("\t"^(v_lvl+1)*"Maximum Gamma discrepancy of wing $(wi):\t$(maximum(abs.(str - uns)))")
        end

        @test prod(isapprox.(uns, str, atol=1e-10))
    end

    # Test neighborhood of every cell
    @test begin
        prod(
            prod(
                prod(
                    ni in neighbors_str[bi][ci] for ni in neighbors_uns[bi][ci]
                    ) for ci in 1:body.ncells
                ) for (bi, body) in enumerate(body.bodies)
            )
    end

    #= Uncomment this to debug
    for (bi, body) in enumerate(body.bodies)
        for ci in 1:body.ncells
            res = prod( ni in neighbors_str[bi][ci] for ni in neighbors_uns[bi][ci])

            if !res
                println("bi = $bi\tci = $ci\t$(neighbors_str[bi][ci]) != $(neighbors_uns[bi][ci])")
            end
        end
    end
    =#

    # Test ∇𝜇 scheme
    if verbose
        println("\t"^(v_lvl+1)*"Maximum ∇𝜇 discrepancy:\t\t\t$(maximum(pnl.norm.(Ugradmu_str - Ugradmu_uns)))")
    end

    @test prod(isapprox.(pnl.norm.(Ugradmu_str - Ugradmu_uns), 0, atol=1e-9))

    # Test lift and drag
    Lerr = (Lstr - Luns) / Lstr
    Derr = (Dstr - Duns) / Dstr

    if verbose
        println()
        @printf "%s%19.19s\t%-7s %-7s %-10s\t%-10s\n" "\t"^(v_lvl+1) "Grid" "Lift" "Drag" "L error" "D error"
        @printf "%s%19.19s\t%-7.1f %-7.1f %-10s\t%-10s\n" "\t"^(v_lvl+1) "GeometricTools" Lstr Dstr "-" "-"
        @printf "%s%19.19s\t%-7.1f %-7.1f %4.3g﹪\t%4.3g﹪\n" "\t"^(v_lvl+1) "Meshes.jl" Luns Duns Lerr*100 Derr*100
    end

    @test Lerr <= 1e-9 && Derr <= 1e-9

end
