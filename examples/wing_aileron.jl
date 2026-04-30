import FLOWPanel as pnl
import FLOWPanel: norm, dot, cross

import Meshes
import GeoIO

# import CUDA                               # Uncomment this to use GPU (if available)

# run_names = ["nasa_wing.msh", "nasa_surface_spaced.msh"]
# run_names = ["nasa_surface_spaced.msh"]
run_names = ["0010.msh"]
file_path       = "examples"
paraview        = true                      # Whether to visualize with Paraview
out_file = joinpath(pnl.examples_path, "wing_aileron", "coupled_timing_results.csv")

files = [joinpath(pnl.examples_path, "wing_aileron", name) for name in run_names]

# ----------------- SIMULATION PARAMETERS --------------------------------------
m              = 0.0254
AOA             = 10.0                      # (deg) freestream angle of attack
magVinf         = 100.0                      # (m/s) freestream velocity
rho             = 1.225                     # (kg/m^3) air density

# ----------------- GEOMETRY DESCRIPTION ---------------------------------------
c_body1 = 10.0
b = 60.0                           # (m) span length
c_body2 = 2.0
AR_body1 = b / c_body1                             # (m) span length
AR_body2 = b / c_body2                             # (m) span length

chords = [c_body1]
ARs = [AR_body1] 

scaling = 1.0
# ----------------- SOLVER SETTINGS -------------------------------------------

# Body and wake model
kernel = Union{pnl.ConstantSource, pnl.ConstantDoublet}               # Kernel type to use
# kernel = pnl.ConstantSource
# body type
# bodytype = pnl.NonLiftingBody{kernel}    # Elements and wake model
bodytype = pnl.RigidWakeBody{kernel}

# Processing
clip_Cp         = 1 - 342.0/magVinf         # Clip pressure coefficients that are lower than this threshold

function postprocess!(bodies, Vinf, magVinf, rho, backend, chords, span)
    Dhat = Vinf / norm(Vinf)
    Shat = [0, 1, 0]
    Lhat = cross(Dhat, Shat)
    Sref = 0.0
    for chord in chords
        Sref += chord * span
        @show Sref
    end

    pnl.calcfield_U!(bodies, Vinf; backend)
    pnl.apply_freestream!(bodies, Vinf)
    pnl.calcfield_Cp!(bodies, magVinf; correct_kuttacondition=fill(true, length(bodies)))
    pnl.calcfield_F!(bodies, magVinf, rho)
    LDS = pnl.calcfield_LDS!(zeros(3,3), bodies, Lhat, Dhat, cross(Lhat, Dhat))

    # Force coefficients
    @show nondim = 0.5 * rho * magVinf^2 * Sref
    CL = sign(dot(LDS[:,1], Lhat)) * norm(LDS[:,1]) / nondim
    CD = sign(dot(LDS[:,2], Dhat)) * norm(LDS[:,2]) / nondim

    return CL, CD
end

# ----------------- GENERATE BODY ---------------------------------------------

function generate_body(
    meshfile::String,
    chord::Float64,
    span::Float64,
    bodytype::Type{<:pnl.NonLiftingBody},
    scaling::Float64 = 1.0,
    flip::Int64 = 1,
    Vinf::AbstractVector{<:Real} = zeros(3)
)
    magVinf = norm(Vinf)

    # Read Gmsh mesh
    msh = GeoIO.load(meshfile).geometry

    # Transform the mesh: scale
    msh = msh |> Meshes.Scale(scaling)

    # Wrap into Grid object
    grid = pnl.gt.GridTriangleSurface(msh)

    # Create trailing edge line
    nte = 200
    trailingedge = zeros(3, nte)
    trailingedge[1, :] .= chord
    trailingedge[2, :] .= range(-span/2, stop=span/2, length=nte)
    trailingedge[3, :] .= 0.0

    # 6enerate the paneled body
    body = bodytype(grid; CPoffset = (-1)^flip * 1e-14)

    pnl.apply_freestream!(body, Vinf)

    return body
end

function generate_body(
    meshfile::String,
    chord::Float64,
    span::Float64,
    bodytype::Type{<:pnl.RigidWakeBody},
    scaling::Float64 = 1.0,
    flip::Int64 = 1,
    Vinf::AbstractVector{<:Real} = zeros(3)
)
    magVinf = norm(Vinf)

    # Read Gmsh mesh
    msh = GeoIO.load(meshfile).geometry

    # Transform the mesh: scale
    msh = msh |> Meshes.Scale(scaling)

    # Wrap into Grid object
    grid = pnl.gt.GridTriangleSurface(msh)

    # Create trailing edge line
    offset = minimum(grid._nodes[1, :])
    nte = 200
    trailingedge = zeros(3, nte)
    trailingedge[1, :] .= chord + offset
    trailingedge[2, :] .= range(-span/2, stop=span/2, length=nte)
    trailingedge[3, :] .= 0.0

    # Generate the paneled body
    body = bodytype(grid; CPoffset = (-1)^flip * 1e-14)
    shedding = pnl.calc_shedding(
        body.nodes,
        body.cells,
        trailingedge;
        tolerance = 0.001 * span
    )
    body = bodytype(body.nodes, body.cells, [shedding];
                    CPoffset = (-1)^flip * 1e-14,
                    ensure_winding=false)

    # initialize wake doublets
    for i in eachindex(body.Das)
        body.Das[i] .= repeat(Vinf/magVinf, 1, size(body.Das[i],2))
    end

    pnl.apply_freestream!(body, Vinf)

    return body
end

# function benchmarks(out, bodies, solver; backend=DirectBackend())
#     @time begin
#         t_build, t_solve = pnl.solve!(bodies, solver; backend, update_G=true)
#     end

#     nps = sum(b.ncells for b in bodies)

#     open(out, "a") do io 
#         write(io, 
#         "$(nps),$(t_build),$(t_solve)\n"
#         )
#         flush(io)
#     end
#     return t_build, t_solve
# end

Vinf = magVinf * [cosd(AOA), sind(AOA), 0.0]

body = generate_body(files[1], chords[1], b, bodytype, scaling, 1, Vinf)
# bodies = tuple([generate_body(file, chord, b, bodytype, scaling, 1, Vinf)
                # for (file, chord) in zip(files, chords)]...)

#------------------- SOLVE BODY ----------------------------------------------
backend = pnl.DirectBackend()
# backend = pnl.FastMultipoleBackend()
solver = pnl.BackslashCoupled((body,))
# solver = pnl.Backslash(body)
# solver = pnl.FGSSolver(body)
println("Solving body...")
# body = body

# benchmarks(out_file, bodies, solver; backend)
pnl.solve!((body,), solver)

# pnl.solve!(bodies, solver)
println("Strength column 1:")
println("  max = ", maximum(body.strength[:, 1]))
println("  min = ", minimum(body.strength[:, 1]))

println("Strength column 2:")
println("  max = ", maximum(body.strength[:, 2]))
println("  min = ", minimum(body.strength[:, 2]))

# ----------------- POST PROCESSING ----------------------------------------
println("Post processing...")

CL, CD = postprocess!((body,), Vinf, magVinf, rho, backend, chords, b)
println("CL =$(CL)")
println("CD =$(CD)")

###
# LDS[:, 1] = Lift
# LDS[:, 2] = Drag
# LDS[:, 3] = Side
###

# Force coefficients
# nondim = 0.5 * rho * magVinf^2 * sum(Sref)
# CL = sign(dot(LDS[:,1], Lhat)) * norm(LDS[:,1]) / nondim
# CD = sign(dot(LDS[:,2], Dhat)) * norm(LDS[:,2]) / nondim

# @show CL
# @show CD

println("done")

#=
# ----------------- VISUALIZATION ----------------------------------------------
# if paraview
    str = save_path*"/"

    # Save wing as a VTK
    str *= pnl.save(body, run_name; path=save_path)
    # str *= pnl.save(body, run_name; path=save_path, out_wake=false)

    # Call Paraview
    # run(`paraview --data=$(str)`)
# end
=#
