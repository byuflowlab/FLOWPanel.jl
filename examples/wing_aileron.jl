import FLOWPanel: norm, dot, cross
import Meshes
import GeoIO
import FLOWPanel as pnl
import GeometricTools as gt
using LinearAlgebra: diag
using StaticArrays: SVector

function generate_body(
    meshfile::String,
    chord::Float64,
    span::Float64,
    bodytype::Type{<:pnl.RigidWakeBody},
    scaling::Float64 = 1.0,
    flip::Int64 = 1,
    Vinf::AbstractVector{<:Real} = zeros(3),
    firstnode=-1,
    secondnode=-1
)
    magVinf = norm(Vinf)

    # Read Gmsh mesh
    msh = GeoIO.load(meshfile).geometry

    # Transform the mesh: scale
    msh = msh |> Meshes.Scale(scaling)

    # Wrap into Grid object
    grid = pnl.gt.GridTriangleSurface(msh)

    # Generate TE shedding matrix
    shedding = zeros(Int, 6, 0)

    # Generate the paneled body
    CPoffset = 1e-6
    body = bodytype(grid, [shedding]; CPoffset, flip_normals=false)
    # pnl.write_vtk("spaced_nasa", body)

    # Recompute shedding from the finalized cell winding used by `body`.
    if firstnode == -1 || secondnode == -1
        @warn "firstnode and secondnode not provided; TE shedding will be disabled. This may cause inaccurate results for lifting bodies."
    else
        shedding = pnl.calc_shedding_from_seed(
            body.nodes,
            body.cells,
            firstnode, secondnode
        )
        body = bodytype(
            body.nodes,
            body.cells,
            [shedding];
            CPoffset,
            flip_normals=false,
            ensure_winding=false
        )
    end

    # initialize wake doublets
    for i in eachindex(body.Das)
        body.Das[i] .= repeat(Vinf/magVinf, 1, size(body.Das[i],2))
    end

    pnl.apply_freestream!(body, Vinf)

    return body
end

function postprocess!(bodies, Vinf, rho, backend, chords, span)
    Dhat = Vinf / norm(Vinf)
    Shat = [0, 1, 0]
    Lhat = cross(Dhat, Shat)
    Sref = 0.0
    for chord in chords
        Sref += chord * span
    end
    
    pnl.calcfield_U!(bodies, Vinf; backend)
    # pnl.apply_freestream!(bodies, Vinf)
    pnl.calcfield_Cp!(bodies, magVinf; correct_kuttacondition=fill(true, length(bodies)))
    pnl.calcfield_F!(bodies, magVinf, rho)
    LDS = pnl.calcfield_LDS!(zeros(3,3), bodies, Lhat, Dhat, cross(Lhat, Dhat))

    # Force coefficients
    nondim = 0.5 * rho * norm(Vinf)^2 * Sref
    CL = sign(dot(LDS[:,1], Lhat)) * norm(LDS[:,1]) / nondim
    CD = sign(dot(LDS[:,2], Dhat)) * norm(LDS[:,2]) / nondim

    return CL, CD
end

"""
Assumes freestream has already been applied.
"""
function flow_tangency_residuals(bodies::Tuple)
    for body in bodies
        pnl.calc_normals!(body)
        pnl.calc_controlpoints!(body; off=1e-10)
    end

    pnl.influence!(bodies, bodies; scalar_potential=false, velocity=true)

    res = zeros(length(bodies))
    for (i, body) in enumerate(bodies)
        r = 0.0
        for (vel, normal) in zip(eachcol(body.velocity), eachcol(body.normals))
            vx, vy, vz = vel
            nx, ny, nz = normal
            r1 = vx * nx + vy * ny + vz * nz
            r += r1 * r1
        end
        res[i] = r / size(body.normals, 2)
    end

    return res
end

function flow_potential_residuals(bodies::Tuple; cp_off=-1e-10)
    for body in bodies
        pnl.calc_normals!(body)
        pnl.calc_controlpoints!(body; off=cp_off)
        body.potential .= 0.0
    end

    pnl.influence!(bodies, bodies; scalar_potential=true, velocity=false)
    res = zeros(length(bodies))
    for (i, body) in enumerate(bodies)
        r = 0.0
        for potential in body.potential
            r += potential * potential
        end
        res[i] = r / length(body.potential)
    end

    return res
end

run_names = ["nasa_wing.msh", "nasa_surface_spaced_repaired.msh"]
file_path       = "examples"
paraview        = true                      # Whether to visualize with Paraview
out_file = joinpath(pnl.examples_path, "wing_aileron", "coupled_timing_results.csv")

files = [joinpath(pnl.examples_path, "wing_aileron", name) for name in run_names]
nodes1 = [42, 43]
nodes2 = [34, 35]
nodes1 = [42, 19]
nodes2 = [34, 3]

# ----------------- SIMULATION PARAMETERS --------------------------------------
m              = 0.0254
AOA             = 0.0                      # (deg) freestream angle of attack
magVinf         = 117.3 * m * 12                      # (m/s) freestream velocity
rho             = 1.225                     # (kg/m^3) air density

# ----------------- GEOMETRY DESCRIPTION ---------------------------------------
c_body1 = 10 * m
b = 60 * m                            # (m) span length
c_body2 = 2 * m
AR_body1 = b / c_body1                             # (m) span length
AR_body2 = b / c_body2                             # (m) span length

chords = [c_body1, c_body2]
ARs = [AR_body1, AR_body2] 
Sref = b * (c_body1 + c_body2)

scaling = 1.0
# ----------------- SOLVER SETTINGS -------------------------------------------

# Body and wake model
kernel = Union{pnl.ConstantSource, pnl.ConstantDoublet}               # Kernel type to use
# body type
bodytype = pnl.RigidWakeBody{kernel}

# Processing
clip_Cp         = 1 - 342.0/magVinf         # Clip pressure coefficients that are lower than this threshold

Vinf = magVinf * [cosd(AOA), sind(AOA), 0.0]

bodies = tuple([generate_body(file, chord, b, bodytype, scaling, 1, Vinf, firstnode, secondnode)
                for (file, chord, firstnode, secondnode) in zip(files, chords, nodes1, nodes2)]...)

#------------------- SOLVE BODY ----------------------------------------------
solver = pnl.BackslashCoupled(bodies)
println("Solving bodies...")

@show nps = sum(b.ncells for b in bodies)

t_build, t_solve = pnl.solve!(bodies, solver; update_G=true)

write_header = !isfile(out_file) || filesize(out_file) == 0

open(out_file, "a") do io
    if write_header
        write(io, "solver,nps,t_build,t_solve,res,pot\n")
    end

    write(io,
        "BackslashCoupled,$(t_build),$(t_solve)\n"
    )
end

println("Resetting bodies...")

for body in bodies
    body.velocity .= 0.0
    pnl.apply_freestream!(body, Vinf)
end

backend2 = fill(pnl.DirectBackend(), length(bodies))
solver2 = (pnl.Backslash(bodies[1]), pnl.Backslash(bodies[2]))

println("Solving bodies part 2...")

t_build2, t_solve2 = pnl.solve!(bodies, solver2)

open(out_file, "a") do io
    write(io,
        "BackslashIterate,$(t_build2),$(t_solve2)\n"
    )
end

println("Resetting bodies...")

for body in bodies
    body.velocity .= 0.0
    pnl.apply_freestream!(body, Vinf)
end

backend3 = fill(pnl.FastMultipoleBackend(), length(bodies))
solver3 = (pnl.KrylovSolver(bodies[1]), pnl.KrylovSolver(bodies[2]))

println("Solving bodies part 2...")

t_build3, t_solve3 = pnl.solve!(bodies, solver3)

open(out_file, "a") do io
    write(io,
        "Krylove,$(t_build3),$(t_solve3)\n"
    )
end

println("done")