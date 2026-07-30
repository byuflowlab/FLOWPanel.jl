# ENV["KMP_DUPLICATE_LIB_OK"] = get(ENV, "KMP_DUPLICATE_LIB_OK", "TRUE")
# ENV["MPLCONFIGDIR"] = get(ENV, "MPLCONFIGDIR", joinpath(@__DIR__, "..", ".mplconfig"))
# ENV["XDG_CACHE_HOME"] = get(ENV, "XDG_CACHE_HOME", joinpath(@__DIR__, "..", ".cache"))

import FLOWPanel as pnl
# include(joinpath(pnl.examples_path, "helper_functions.jl"))
import FLOWPanel: norm, dot, cross

# import Meshes
import GeoIO

import PythonPlot as plt
import PythonPlot: @L_str

# import CUDA                               # Uncomment this to use GPU (if available)

function plot_Cp(wing, spanposs, b, c, rho, magVinf; 
    clearme=true, 
    spanpos_offset=0.0, 
    colormap = plt.get_cmap("RdBu",15),
    colormap_indices=[2, 5, 8, 11],
    x_cp_upper = [nothing],
    x_cp_lower = [nothing],
    slicetol = 0.02 * b
)
    
    colors = [colormap(i) for i in colormap_indices]
    qinf = 0.5 * rho * magVinf^2

    figname = "simplewing"
    fig = plt.figure(figname, figsize=(8,6))
    if clearme
        plt.clf()
        ax = fig.add_subplot(111, xlabel=L"x/c", ylabel=L"C_P")
    else
        axs = fig.get_axes()
        ax = length(axs) > 0 ? axs[0] : fig.add_subplot(111, xlabel=L"x/c", ylabel=L"C_P")
    end

    stls = ["-", "--", "-.", ":"]

    for (i, spanpos) in enumerate(spanposs)
        points, slicePs = pnl.slice_scalarfield(wing, :P, 2, spanpos * b * (1-spanpos_offset), slicetol)
        Cps = slicePs ./ qinf
        ax.plot(points[1, :] ./ c, Cps, label="2y/b=$(round(spanpos, digits=2))", stls[i], color=colors[i])
        if x_cp_upper[i] != nothing
            ax.scatter(x_cp_upper[i][:, 1], x_cp_upper[i][:, 2], marker="^", label="Exp. (upper)", color=colors[i])
        end
        if x_cp_lower[i] != nothing
            ax.scatter(x_cp_lower[i][:, 1], x_cp_lower[i][:, 2], marker="v", label="Exp. (lower)", color=colors[i])
        end
    end

    ax.invert_yaxis()
    ax.legend()
    plt.tight_layout()
end

function scatter_Cp(wing, spanposs, b, c, rho, magVinf; 
    clearme=true, 
    spanpos_offset=0.0, 
    colormap = plt.get_cmap("RdBu",15),
    colormap_indices=[2, 5, 8, 11],
    x_cp_upper = [nothing],
    x_cp_lower = [nothing],
    slicetol = 0.02 * b,
    labels = ["2y/b=$(round(spanpos, digits=2))" for spanpos in spanposs],
    ss = [50 for _ in spanposs],
    mrks = ["o", "s", "^", "v"]
)
    
    colors = [colormap(i) for i in colormap_indices]
    qinf = 0.5 * rho * magVinf^2

    figname = "simplewing"
    fig = plt.figure(figname, figsize=(8,6))
    if clearme
        plt.clf()
        ax = fig.add_subplot(111, xlabel=L"x/c", ylabel=L"C_P")
        ax.invert_yaxis()
    else
        axs = fig.get_axes()
        ax = length(axs) > 0 ? axs[0] : fig.add_subplot(111, xlabel=L"x/c", ylabel=L"C_P")
    end

    # stls = ["-", "--", "-.", ":"]

    # points, slicePs = pnl.slice_scalarfield(wing, :P, 2, spanpos * b * (1-spanpos_offset), slicetol)
    plot_data = [pnl.slice_scalarfield(wing, :P, 2, spanpos * b * (1-spanpos_offset), slicetol) for spanpos in spanposs]
    for (i, spanpos) in enumerate(spanposs)
        points, slicePs = plot_data[i]
        Cps = slicePs ./ qinf
        ax.scatter(points[1, :] ./ c, Cps, label=labels[i], color=colors[i], marker=mrks[i], s=ss[i])
        if x_cp_upper[i] != nothing
            ax.scatter(x_cp_upper[i][:, 1], x_cp_upper[i][:, 2], marker="^", label="Exp. (upper)", color=colors[i])
        end
        if x_cp_lower[i] != nothing
            ax.scatter(x_cp_lower[i][:, 1], x_cp_lower[i][:, 2], marker="v", label="Exp. (lower)", color=colors[i])
        end
    end

    ax.legend()
    plt.tight_layout()

    return plot_data
end

function run(; AOA=5.0, magVinf=56.0, endplates=false, meshfile="",
        kernel = Union{pnl.ConstantSource, pnl.VortexRing}
    )

    println("\n#===== SIMPLE WING SIMULATION =====#\n")
    println("\tAOA: $(AOA) degrees")
    println("\tmagVinf: $(magVinf) m/s\n")

    run_name        = "wing_capped"             # Name of this run

    save_path       = joinpath("data", run_name) # Where to save outputs
    paraview        = true                      # Whether to visualize with Paraview

    # ----------------- SIMULATION PARAMETERS --------------------------------------
    rho             = 1.225                     # (kg/m^3) air density
    AR = 4.0                               # Aspect ratio of the wing (b/c)
    c = 2.0                               # (m) root chord length
    b = AR * c                            # (m) span length


    # ----------------- GEOMETRY DESCRIPTION ---------------------------------------
    
    stretch = 1.0
    b = 2.7 * stretch
    c = 0.76
    AR = b/c
    # trailingedgefile= joinpath(read_path, "zeroebwb-TE.msh") # Gmsh file with trailing edge

    # offset          = [0, 0, 0]                 # Offset to center the mesh
    # rotation        = RotZ(-90*pi/180)*RotX(90*pi/180) # Rotation to align mesh
    scaling         = 1.0                      # Factor to scale original mesh to
                                                # the approximate dimensions of the
                                                # ZEROEe BWB subscale model

    spandir         = [0, 1, 0]                 # Span direction used to orient the trailing edge
    flip            = false                     # Whether to flip control points against the direction of normals
                                                # NOTE: use `flip=true` if the normals
                                                #       point inside the body

    Sref            = c * b                       # (m^2) reference area

    # ----------------- SOLVER SETTINGS -------------------------------------------

    # Body and wake model
    # kernel = pnl.ConstantSource               # Kernel type to use
    # kernel = pnl.ConstantDoublet               # Kernel type to use
    # kernel = Union{pnl.ConstantSource, pnl.ConstantDoublet}               # Kernel type to use
    # kernel = Union{pnl.ConstantSource, pnl.VortexRing}               # Kernel type to use
    # kernel = pnl.VortexRing               # Kernel type to use

    # body type
    # bodytype = pnl.NonLiftingBody{kernel}    # Elements and wake model
    DBC = pnl.kernel_dim(kernel) == 2
    bodytype = pnl.RigidWakeBody{kernel, pnl.kernel_dim(kernel), Float64, DBC}    # Elements and wake model

    # Processing
    clip_Cp         = 1 - 342.0/magVinf         # Clip pressure coefficients that are lower than this threshold


    # ----------------- GENERATE BODY ----------------------------------------------
    # Read Gmsh mesh
    println("\tLoading mesh...")
    msh = GeoIO.load(meshfile).geometry
    println("\tDone.\n")

    # Transform the original mesh: Translate, rotate, and scale
    # msh = msh |> Meshes.Scale(scaling)

    # Uncomment this to do 10 smoothing iterations on the mesh
    # msh = msh |> Meshes.TaubinSmoothing(10)

    # get trailing edge line
    nte = 10000
    trailingedge = zeros(3, nte)
    trailingedge[1, :] .= c 
    trailingedge[2, :] .= range(-b/2, stop=b/2, length=nte)
    trailingedge[2, :] .= range(0.0, stop=b, length=nte)
    trailingedge[3, :] .= 0.0

    # Freestream vector
    Vinf = magVinf*[cos(AOA*pi/180), 0, sin(AOA*pi/180)]

    # Generate paneled body
    println("\tGenerating body...")
    nodes, cells = pnl.meshes2nodes_cells(msh)
    nodes[2, :] .*= stretch # stretch aspect ratio
    # cells = hcat(cells[:, 1:2044], cells[:, 2046:end]) # remove 1 panel for Neumann BC
    
    if bodytype == pnl.NonLiftingBody{pnl.ConstantSource}
        body = bodytype(nodes, cells; cp_outer=iseven(flip))
    elseif bodytype <: pnl.RigidWakeBody
        shedding = pnl.calc_shedding(nodes, cells, trailingedge; tolerance=0.001*b)
        # shedding2 = pnl.calc_shedding_from_seed(nodes, cells, 396, 364; bbox=nothing, end_node=nothing, normal_jump_tol=0.2, max_turn_angle=pi/3, debug=false)
        body = bodytype(nodes, cells, shedding; cp_outer=true, kerneloffset=1e-3, ensure_winding=true, semiinfinite_wake=true)
        shedding = pnl.calc_shedding(body.nodes, body.cells, trailingedge; tolerance=0.001*b)
        body = bodytype(nodes, cells, shedding; cp_outer=true, kerneloffset=1e-3, ensure_winding=true, semiinfinite_wake=true)
        
        body.Das[1] .= repeat(Vinf/magVinf, 1, size(body.Das[1], 2))
    else
        error("Unsupported body type")
    end

    Uinfs = repeat(Vinf, 1, body.ncells)

    println("Number of panels:\t$(body.ncells)")

    #------------------- GENERATE END PLATES ----------------------------------------------
    
    left_center = [c*0.5, 0.0, 0.0]
    left_normal = [0.0, 1.0, 0.0]
    left_radius = c * 20
    left_plate = pnl.FlatGround(left_center, left_normal, left_radius; panel_length=c*0.1)
    right_center = [c*0.5, b, 0.0]
    right_normal = [0.0, -1.0, 0.0]
    right_radius = c * 20
    right_plate = pnl.FlatGround(right_center, right_normal, right_radius; panel_length=c*0.1)

    #------------------- SOLVE BODY ----------------------------------------------
    println("\tSolving body...")

    # Unitary direction of semi-infinite vortex at points `a` and `b` of each
    # trailing edge panel
    # body.Das .= repeat(Vinf/magVinf, 1, body.nsheddings+1)

    # select backend for n-body calculation
    backend = pnl.FastMultipoleBackend(
            expansion_order=7,
            multipole_acceptance=0.4,
            leaf_size=20
        )
    # backend = pnl.DirectBackend()
        
    # Solve body (panel strengths) giving `Uinfs` as boundary conditions and
    @time begin
        # global solver = pnl.Backslash(body; least_squares=true)
        # solver = pnl.KrylovSolver(body;
        #     method=:gmres,
        #     itmax=40,
        #     atol=1e-4,
        #     rtol=1e-4,
        #     # elprescribe=Tuple{Int,Float64}[],   # No prescribed strengths
        #     backend=pnl.FastMultipoleBackend(
        #                 expansion_order=7,
        #                 multipole_acceptance=0.4,
        #                 leaf_size=10
        #             )
        # )

        # single body solve
        bodies = body
        solvers = pnl.Backslash(body)

        # endplate solve
        # bodies = (left_plate, right_plate, body)
        pnl.apply_freestream!(bodies, Vinf)
        if endplates
            solvers = (pnl.FlatGroundSolver(left_plate), pnl.FlatGroundSolver(right_plate), pnl.Backslash(body))
            pnl.solve!(bodies, solvers; backend)
        else
            solver = pnl.Backslash(body)
            pnl.solve!(body, solver; backend)
        end

    end

    # ----------------- POST PROCESSING ----------------------------------------
    println("\tPost processing...")

    # Calculate surface velocity U on the body
    pnl.calcfield_U!((body,), Vinf; backend=backend)

    # Calculate gauge pressure (based on U + U_∇μ)
    pnl.calcfield_P!(body, magVinf, rho; correct_kuttacondition=false)

    # Calculate the force of each panel (based on P)
    pnl.calcfield_F!(body)

    # --------- Integrated forces: lift and induced drag

    # Calculate total force of the vehicle decomposed as lift, drag, and sideslip
    Dhat = Vinf/norm(Vinf)        # Drag direction
    Shat = [0, 1, 0]              # Span direction
    Lhat = cross(Dhat, Shat)      # Lift direction

    LDS = pnl.calcfield_LDS(body, Lhat, Dhat)

    L = LDS[:, 1]
    D = LDS[:, 2]

    # Force coefficients
    nondim = 0.5*rho*magVinf^2*Sref   # Normalization factor
    CL = sign(dot(L, Lhat)) * norm(L) / nondim
    CD = sign(dot(D, Dhat)) * norm(D) / nondim

    @show CL
    @show CD


    # ----------------- VISUALIZATION ----------------------------------------------
    # if paraview
        mkpath(save_path)
        str = pnl.write_vtk(joinpath(save_path, run_name*"_AOA$(AOA)"), body, 0, 0.0; overwrite=true)
        str = pnl.write_vtk(joinpath(save_path, run_name*"_AOA$(AOA)_leftplate"), left_plate, 0, 0.0; overwrite=true)
        str = pnl.write_vtk(joinpath(save_path, run_name*"_AOA$(AOA)_rightplate"), right_plate, 0, 0.0; overwrite=true)

        # Call Paraview
        # run(`paraview --data=$(str)`)
    # end

    return CL, CD, bodies, b, c, rho, magVinf
end

read_path       = joinpath(pnl.examples_path, "data") # Where to read Gmsh files from
# meshfile        = joinpath(read_path, "wing_ar4_naca0016_refined.msh")    # Gmsh file to read
# meshfile        = joinpath(read_path, "wing_ar4_naca0016_5.msh")    # Gmsh file to read
# meshfile        = joinpath(read_path, "naca0012_nc133.msh")    # Gmsh file to read
meshfile        = joinpath(read_path, "naca0012_nc101_nw26_refined.msh")    # Gmsh file to read
# shedding_points = [1972, 1900, 0] .+ 1 # Nodes from which to seed shedding (1-based indexing)
# meshfile        = joinpath(read_path, "naca0012_a.msh")    # Gmsh file to read
# meshfile        = joinpath(read_path, "naca0012_nc70_nocaps.msh")    # Gmsh file to read
# meshfile        = joinpath(read_path, "naca0012_nc70_cfd.msh")    # Gmsh file to read
# meshfile        = joinpath(read_path, "naca0012_nc70_cfd_refined.msh")    # Gmsh file to read
# meshfile        = joinpath(read_path, "naca0012_nc70_refineLTE.msh")    # Gmsh file to read
# meshfile        = joinpath(read_path, "naca0012_nc70_refineLE.msh")    # Gmsh file to read

kernel = Union{pnl.ConstantSource, pnl.VortexRing}

# CL0, _ = run(;AOA=0.0)
CL6, _, bodies, b, c, rho, magVinf = run(;AOA=7.0, meshfile, kernel)
# CL6, _, bodies, b, c, rho, magVinf = run(;AOA=5.88, meshfile)
# CL10, _ = run(;AOA=10.0)

# slope = (CL10 - CL0) / (10.0 - 0.0) * 180.0 / π

x_cp_upper = [
      0.0                    -0.6549295774647887
      0.0021008403361344537  -1.732394366197183
      0.004201680672268907   -2.23943661971831
      0.0063025210084033615  -2.471830985915493
      0.014705882352941176   -2.3661971830985915
      0.018907563025210083   -2.176056338028169
      0.04831932773109243    -1.7112676056338028
      0.07563025210084033    -1.436619718309859
      0.09873949579831932    -1.3309859154929577
      0.12394957983193276    -1.1408450704225352
      0.14705882352941177    -1.0985915492957745
      0.17436974789915966    -0.9507042253521126
      0.28361344537815125    -0.7394366197183099
      0.31932773109243695    -0.6971830985915493
      0.38235294117647056    -0.5704225352112676
      0.4432773109243697     -0.5070422535211268
      0.5084033613445378     -0.4225352112676056
      0.5756302521008403     -0.3380281690140845
      0.6386554621848739     -0.2535211267605634
      0.703781512605042      -0.2112676056338028
      0.76890756302521       -0.1267605633802817
      0.8361344537815125     -0.06338028169014084
      0.8949579831932772      0.02112676056338028
      0.9600840336134453      0.1267605633802817
      1.0                     0.16901408450704225
  ]
  
x_cp_lower = [
      0.0021008403361344537   0.0
      0.004201680672268907    0.8028169014084507
      0.008403361344537815    0.9929577464788732
      0.02310924369747899     0.8028169014084507
      0.03991596638655462     0.6126760563380281
      0.06092436974789916     0.4436619718309859
      0.10504201680672269     0.2535211267605634
      0.16176470588235292     0.1267605633802817
      0.2310924369747899      0.06338028169014084
      0.2899159663865546      0.02112676056338028
      0.35084033613445376     0.0
      0.40756302521008403     0.0
      0.4852941176470588     -0.02112676056338028
      0.5714285714285714     -0.02112676056338028
      0.6491596638655462      0.0
      0.7668067226890756      0.02112676056338028
      0.8760504201680672      0.08450704225352113
      0.9243697478991596      0.1056338028169014
      1.0                     0.14788732394366197
  ]

# plot_Cp(bodies[3], [0.5], b, c, rho, magVinf; 
#     clearme=true, x_cp_lower=[x_cp_lower], x_cp_upper=[x_cp_upper],
#     slicetol=0.02*b)

# data_str = scatter_Cp(bodies[3], [0.5], b, c, rho, magVinf; 
#     clearme=false, 
#     # spanpos_offset=0.0, 
#     colormap = plt.get_cmap("RdBu",15),
#     colormap_indices=[11],#[2, 5, 8, 11],
#     x_cp_upper = [nothing],
#     x_cp_lower = [nothing],
#     slicetol = 0.02 * b,
#     # labels = ["2y/b=$(round(spanpos, digits=2))" for spanpos in spanposs]
#     labels = ["structured"],
#     ss = [20], mrks = ["x"]
# )
