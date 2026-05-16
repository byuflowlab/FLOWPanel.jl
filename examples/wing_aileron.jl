import FLOWPanel: norm, dot, cross
import Meshes
import GeoIO
import FLOWPanel as pnl
import GeometricTools as gt
using LinearAlgebra: diag
using StaticArrays: SVector

## May use kerneloffset for vortexring bodies

function generate_body(
    meshfile::String,
    chord::Float64,
    span::Float64,
    bodytype::Type{<:pnl.RigidWakeBody};
    translate::NTuple{3,Float64} = (0.0, 0.0, 0.0),
    scaling::Float64 = 1.0,
    flip::Int64 = 1,
    Vinf::AbstractVector{<:Real} = zeros(3),
    kerneloffset::Float64 = 1.0e-6,
    firstnode = -1,
    secondnode = -1
)

    magVinf = norm(Vinf)

    # --- Load + transform mesh ---
    msh = GeoIO.load(meshfile).geometry
    msh = msh |> Meshes.Scale(scaling)
    msh = msh |> Meshes.Translate(translate)

    grid = pnl.gt.GridTriangleSurface(msh)

    CPoffset = 1e-6
    empty_shedding = zeros(Int, 6, 0)

    # --- Step 1: build body WITHOUT relying on TE shedding ---
    body = bodytype(
        grid,
        [empty_shedding];
        CPoffset,
        kerneloffset,
        flip_normals = false
    )

    # --- Step 2: find trailing edge nodes (YOU PROVIDE THIS FUNCTION) ---
    if firstnode == -1 || secondnode == -1
        # Replace this with your TE detection routine
        firstnode, secondnode = find_trailing_edge_nodes(body)
    end

    # --- Step 3: build shedding properly from finalized topology ---
    shedding = pnl.calc_shedding_from_seed(
        body.nodes,
        body.cells,
        firstnode,
        secondnode
    )

    body = bodytype(
        body.nodes,
        body.cells,
        [shedding];
        CPoffset,
        kerneloffset,
        flip_normals = false,
        ensure_winding = false
    )

    # --- initialize wake doublets ---
    for i in eachindex(body.Das)
        body.Das[i] .= repeat(Vinf / magVinf, 1, size(body.Das[i], 2))
    end

    return body
end

function postprocess!(bodies, Vinf, rho, chords, span, scaling::Float64=1.0)
    Dhat = Vinf / norm(Vinf)
    Shat = [0, 1, 0]
    Lhat = cross(Dhat, Shat)
    Sref = 0.0
    for chord in chords
        Sref += chord * span * scaling * scaling
    end
    
    pnl.calcfield_U!(bodies, Vinf)
    pnl.calcfield_P!(bodies, magVinf, rho)
    pnl.calcfield_F!(bodies)
    LDS = pnl.calcfield_LDS!(zeros(3,3), bodies, Lhat, Dhat, cross(Lhat, Dhat))

    # Force coefficients
    @show nondim = 0.5 * rho * norm(Vinf)^2 * Sref
    CL = sign(dot(LDS[:,1], Lhat)) * norm(LDS[:,1]) / nondim
    CD = sign(dot(LDS[:,2], Dhat)) * norm(LDS[:,2]) / nondim

    return CL, CD
end

# """
# Assumes freestream has already been applied.
# """
# function flow_tangency_residuals(bodies::Tuple)
#     for body in bodies
#         pnl.calc_normals!(body)
#         pnl.calc_controlpoints!(body; off=1e-10)
#     end

#     pnl.influence!(bodies, bodies; scalar_potential=false, velocity=true)

#     res = zeros(length(bodies))
#     for (i, body) in enumerate(bodies)
#         r = 0.0
#         for (vel, normal) in zip(eachcol(body.velocity), eachcol(body.normals))
#             vx, vy, vz = vel
#             nx, ny, nz = normal
#             r1 = vx * nx + vy * ny + vz * nz
#             r += r1 * r1
#         end
#         res[i] = r / size(body.normals, 2)
#     end

#     return res
# end

# function flow_potential_residuals(bodies::Tuple; cp_off=-1e-10)
#     for body in bodies
#         pnl.calc_normals!(body)
#         pnl.calc_controlpoints!(body; off=cp_off)
#         body.potential .= 0.0
#     end

#     pnl.influence!(bodies, bodies; scalar_potential=true, velocity=false)
#     res = zeros(length(bodies))
#     for (i, body) in enumerate(bodies)
#         r = 0.0
#         for potential in body.potential
#             r += potential * potential
#         end
#         res[i] = r / length(body.potential)
#     end

#     return res
# end

run_names       = ["wing.msh", "surface.msh"]
file_path       = "examples"
paraview        = true                      # Whether to visualize with Paraview
out_file        = joinpath(pnl.examples_path, "wing_aileron", "coupled_timing_results.csv")

function find_trailing_edge_nodes(body; tol=1e-8)
    nodes = body.nodes
    cells = body.cells

    xs = nodes[1, :]

    xmax = maximum(xs)
    # Find the node with the maximum x-coordinate (trailing edge candidate)
    candidates = findall(x -> abs(x - xmax) < tol, xs)
    # If there are multiple candidates, break ties using the y-coordinate
    @show candidates[end], candidates[end-1]
    return candidates[end], candidates[end-1]
end

# ----------------- SIM SETUP ---------------------------------------------------
run_names = ["wing.msh", "surface.msh"]
files = [joinpath(pnl.examples_path, "wing_aileron", name) for name in run_names]
# nodes1 = [788, 788]
# nodes2 = [768, 768]
# nodes1 = [-1, -1]
# nodes2 = [-1, -1]

# ----------------- SIMULATION PARAMETERS --------------------------------------
m              = 0.0254
AOAs = [
-9.58790170132325  
# -7.841209829867704 
-5.988657844990541 
# -4.136105860113403 
# -2.3894139886577754
-0.536862003780719 
#  1.3156899810964084
 3.168241965973529 
#  5.02079395085066  
 6.926275992438619 
#  8.672967863894144 
10.472589792060482 
# 12.325141776937613 
14.230623818525515 
# 15.077504725897917 
# 16.400756143667294 
18.465028355387517 
]

# AOAs = [-8.0]

magVinf         = 117.3 * 12 * m                      # (m/s) freestream velocity
rho             = 1.225                     # (kg/m^3) air density

# ----------------- GEOMETRY DESCRIPTION ---------------------------------------
c_body1 = 10.0
b = 60.0                        # (m) span length
c_body2 = 2.0 
AR_body1 = b / c_body1                             # (m) span length
AR_body2 = b / c_body2                             # (m) span length

chords = [c_body1, c_body2]
ARs = [AR_body1, AR_body2] 
Sref = b * (c_body1 + c_body2)
scaling = m
trs = [
    (0.0, 0.0, 0.0),
    (9.8489*m, 0.0, -0.12*m)
]
# ----------------- SOLVER SETTINGS -------------------------------------------

# Body and wake model
kernel = Union{pnl.ConstantSource, pnl.ConstantDoublet}               # Kernel type to use
# kernel = Union{pnl.ConstantSource, pnl.VortexRing}               # Kernel type to use
# body type
bodytype = pnl.RigidWakeBody{kernel}

# Processing
clip_Cp         = 1 - 342.0/magVinf         # Clip pressure coefficients that are lower than this threshold
paraview        = false

for (i, AOA) in enumerate(AOAs)
    Vinf = magVinf * [cosd(AOA), 0.0, sind(AOA)]

    bodies = tuple([generate_body(file, chord, b, bodytype; translate=tr, scaling=scaling, Vinf=Vinf)
                    for (file, chord, tr) in zip(files, chords, trs)]...)
    
    #------------------- SOLVE BODY ----------------------------------------------
 
    for body in bodies
        body.velocity .= 0.0
        pnl.apply_freestream!(body, Vinf)
    end

    backend = pnl.DirectBackend()
    solver = pnl.BackslashCoupled(bodies)
    println("Solving body...")

    nps = sum(b.ncells for b in bodies)

    t_build, t_solve = pnl.solve!(bodies, solver; update_G=true)

    CL, CD = postprocess!(bodies, Vinf, rho, chords, b, scaling)

    open(out_file, "a") do io
        if i == 1
            write(io, "BackslashCoupled\n")
        end
        write(io,
            "$AOA,$CL,$CD,$(t_build),$(t_solve)\n"
        )
    end

    if i == 1 && paraview
        filestr1 = pnl.write_vtk(joinpath("examples", "wing_val"), bodies[1], 0, 0.0)
        files1 = split(filestr1, ", ")
        pvd1 = first(filter(f -> endswith(f, ".pvd"), files1))

        filestr2 = pnl.write_vtk(joinpath("examples", "surface_val"), bodies[2], 0, 0.0)
        files2 = split(filestr2, ", ")
        pvd2 = first(filter(f -> endswith(f, ".pvd"), files2))

        run(`paraview $pvd1 $pvd2`, wait=false)
    end
end

# for (i, AOA) in enumerate(AOAs)
#     Vinf = magVinf * [cosd(AOA), 0.0, sind(AOA)]

#     bodies = tuple([generate_body(file, chord, b, bodytype; translate=tr, scaling=scaling, Vinf=Vinf, firstnode=firstnode, secondnode=secondnode)
#                     for (file, chord, tr, firstnode, secondnode) in zip(files, chords, trs, nodes1, nodes2)]...)

#     for body in bodies
#         body.velocity .= 0.0
#         pnl.apply_freestream!(body, Vinf)
#     end

#     solver2 = (pnl.Backslash(bodies[1]),pnl.Backslash(bodies[2]))

#     println("Solving bodies part 2...")

#     t_build2, t_solve2 = pnl.solve!(bodies, solver2)
#     CL, CD = postprocess!(bodies, Vinf, rho, chords, b, scaling)

#     open(out_file, "a") do io
#         if i == 1
#             write(io, "BackslashIterative\n")
#         end
#         write(io,
#             "$AOA,$CL,$CD,$(t_build2),$(t_solve2)\n"
#         )
#     end

#     if i == 1 && paraview
#         filestr1 = pnl.write_vtk(joinpath("examples", "wing_val"), bodies[1], 0, 0.0)
#         files1 = split(filestr1, ", ")
#         pvd1 = first(filter(f -> endswith(f, ".pvd"), files1))

#         filestr2 = pnl.write_vtk(joinpath("examples", "surface_val"), bodies[2], 0, 0.0)
#         files2 = split(filestr2, ", ")
#         pvd2 = first(filter(f -> endswith(f, ".pvd"), files2))

#         run(`paraview $pvd1 $pvd2`, wait=false)
#     end

#     println("Resetting bodies...")
# end

# for (i, AOA) in enumerate(AOAs)
#     Vinf = magVinf * [cosd(AOA), sind(AOA), 0.0]
#     bodies = tuple([generate_body(file, chord, b, bodytype; translate=tr, scaling=scaling, Vinf=Vinf, firstnode=firstnode, secondnode=secondnode)
#                     for (file, chord, tr, firstnode, secondnode) in zip(files, chords, trs, nodes1, nodes2)]...)

#     for body in bodies
#         body.velocity .= 0.0
#         pnl.apply_freestream!(body, Vinf)
#     end

#     backend3 = fill(pnl.FastMultipoleBackend(), length(bodies))
#     solver3 = (pnl.KrylovSolver(bodies[1]), pnl.KrylovSolver(bodies[2]))

#     t_build3, t_solve3 = pnl.solve!(bodies, solver3)
#     CL, CD = postprocess!(bodies, Vinf, rho, chords, b)

#     open(out_file, "a") do io
#         if i == 1
#             write(io, "KrylovSolver\n")
#         end
#         write(io,
#             "$AOA,$CL,$CD,$(t_build3),$(t_solve3)\n"
#         )
#     end

#     println("Resetting bodies...")
# end

# for (i, AOA) in enumerate(AOAs)
#     Vinf = magVinf * [cosd(AOA), 0.0, sind(AOA)]
#     bodies = tuple([generate_body(file, chord, b, bodytype; translate=tr, scaling=scaling, Vinf=Vinf, firstnode=firstnode, secondnode=secondnode)
#                     for (file, chord, tr, firstnode, secondnode) in zip(files, chords, trs, nodes1, nodes2)]...)

#     for body in bodies
#         body.velocity .= 0.0
#         pnl.apply_freestream!(body, Vinf)
#     end

#     solver4 = pnl.KrylovCoupled(bodies)
#     t_build3, t_solve3 = pnl.solve!(bodies, solver4)
#     CL, CD = postprocess!(bodies, Vinf, rho, chords, b)

#     open(out_file, "a") do io
#         if i == 1
#             write(io, "KrylovCoupled\n")
#         end
#         write(io,
#             "$AOA,$CL,$CD,$(t_build3),$(t_solve3)\n"
#         )
#     end

#     println("Resetting bodies...")
# end