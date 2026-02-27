import FLOWPanel as pnl
import FLOWPanel: norm, dot, cross

import Meshes
import GeoIO

# import CUDA                               # Uncomment this to use GPU (if available)


run_name        = "wing_capped"             # Name of this run

save_path       = run_name                  # Where to save outputs
paraview        = true                      # Whether to visualize with Paraview
read_path       = joinpath(pnl.examples_path, "data") # Where to read Gmsh files from

# ----------------- SIMULATION PARAMETERS --------------------------------------
AOA             = 15.0                      # (deg) freestream angle of attack
magVinf         = 1.0                      # (m/s) freestream velocity
rho             = 1.225                     # (kg/m^3) air density
AR = 4.0                               # Aspect ratio of the wing (b/c)
c = 2.0                               # (m) root chord length
b = AR * c                            # (m) span length


# ----------------- GEOMETRY DESCRIPTION ---------------------------------------
meshfile        = joinpath(read_path, "wing_ar4_naca0016_refined.msh")    # Gmsh file to read
# meshfile        = joinpath(read_path, "wing_ar4_naca0016_5.msh")    # Gmsh file to read
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
kernel = Union{pnl.ConstantSource, pnl.ConstantDoublet}               # Kernel type to use

# body type
# bodytype = pnl.NonLiftingBody{kernel}    # Elements and wake model
bodytype = pnl.RigidWakeBody{kernel}    # Elements and wake model

# Processing
clip_Cp         = 1 - 342.0/magVinf         # Clip pressure coefficients that are lower than this threshold


# ----------------- GENERATE BODY ----------------------------------------------
# Read Gmsh mesh
msh = GeoIO.load(meshfile).geometry

# Transform the original mesh: Translate, rotate, and scale
msh = msh |> Meshes.Scale(scaling)

# Uncomment this to do 10 smoothing iterations on the mesh
# msh = msh |> Meshes.TaubinSmoothing(10)

# Wrap Meshes object into a Grid object from GeometricTools
grid = pnl.gt.GridTriangleSurface(msh)

# get trailing edge line
nte = 10000
trailingedge = zeros(3, nte)
trailingedge[1, :] .= c 
trailingedge[2, :] .= range(-b/2, stop=b/2, length=nte)
trailingedge[3, :] .= 0.0

# Generate TE shedding matrix
# TE_indices = [161, 129, 97, 65, 3, 1, 268, 300, 332, 364, 396]
# shedding = pnl.calc_shedding(grid, TE_indices, trailingedge; tolerance=0.001*b)
shedding = pnl.calc_shedding(grid, trailingedge; tolerance=0.001*b)

# Freestream vector
Vinf = magVinf*[cos(AOA*pi/180), 0, sin(AOA*pi/180)]

# Freestream at every control point
Uinfs = repeat(Vinf, 1, body.ncells)

# Generate paneled body
if bodytype == pnl.NonLiftingBody{pnl.ConstantSource}
    body = bodytype(grid; CPoffset=(-1)^flip * 1e-14)
elseif bodytype <: pnl.RigidWakeBody
    body = bodytype(grid, shedding; CPoffset=(-1)^flip * 1e-14)
    body.Das .= repeat(Vinf/magVinf, 1, body.nsheddings)
    body.Dbs .= repeat(Vinf/magVinf, 1, body.nsheddings)
else
    error("Unsupported body type")
end

println("Number of panels:\t$(body.ncells)")

#------------------- SOLVE BODY ----------------------------------------------
println("Solving body...")

# Unitary direction of semi-infinite vortex at points `a` and `b` of each
# trailing edge panel
# body.Das .= repeat(Vinf/magVinf, 1, body.nsheddings)
# body.Dbs .= repeat(Vinf/magVinf, 1, body.nsheddings)

# select backend for n-body calculation
backend = pnl.FastMultipoleBackend(
        expansion_order=7,
        multipole_acceptance=0.4,
        leaf_size=100000
    )
# backend = pnl.DirectBackend()
    
# Solve body (panel strengths) giving `Uinfs` as boundary conditions and
@time begin
    # global solver = pnl.Backslash(body; least_squares=true)
    # solver = pnl.KrylovSolver(body;
    #     method=:gmres,
    #     itmax=20,
    #     atol=1e-4,
    #     rtol=1e-4,
    #     # elprescribe=Tuple{Int,Float64}[],   # No prescribed strengths
    #     backend=pnl.FastMultipoleBackend(
    #                 expansion_order=7,
    #                 multipole_acceptance=0.4,
    #                 leaf_size=10
    #             )
    # )
    solver = pnl.BackslashDirichlet(body)
    # solver = pnl.Backslash(body; least_squares=false)
    pnl.solve2!(body, Uinfs, solver; backend)
end

# ----------------- POST PROCESSING ----------------------------------------
println("Post processing...")

# Calculate surface velocity U on the body
@time Us = pnl.calcfield_U(body, body; backend)

# Calculate surface velocity U_∇μ due to the gradient of the doublet strength
Gammai = kernel == pnl.VortexRing ? 1 : kernel == Union{pnl.ConstantSource, pnl.ConstantDoublet} ? 2 : 0
# UDeltaGamma = pnl.calcfield_Ugradmu(body; Gammai)
UDeltaGamma = pnl.calcfield_Ugradmu(body; Gammai=2)

# Add both velocities together
pnl.addfields(body, "Ugradmu", "U")

# Calculate pressure coefficient (based on U + U_∇μ)
@time Cps = pnl.calcfield_Cp(body, magVinf)

# Calculate the force of each panel (based on Cp)
@time Fs = pnl.calcfield_F(body, magVinf, rho)

# --------- Integrated forces: lift and induced drag

# Calculate total force of the vehicle decomposed as lift, drag, and sideslip
Dhat = Vinf/norm(Vinf)        # Drag direction
Shat = [0, 1, 0]              # Span direction
Lhat = cross(Dhat, Shat)      # Lift direction

LDS = pnl.calcfield_LDS(body, Lhat, Dhat)

L = LDS[:, 1]
D = LDS[:, 2]

# Force coefficients
nondim = 0.5*rho*magVinf^2*b^2/AR   # Normalization factor
CL = sign(dot(L, Lhat)) * norm(L) / nondim
CD = sign(dot(D, Dhat)) * norm(D) / nondim

@show CL
@show CD


# ----------------- VISUALIZATION ----------------------------------------------
# if paraview
    str = save_path*"/"

    # Save wing as a VTK
    str *= pnl.save(body, run_name; path=save_path)
    # str *= pnl.save(body, run_name; path=save_path, out_wake=false)

    # Call Paraview
    # run(`paraview --data=$(str)`)
# end