import FLOWPanel as pnl
import FLOWPanel: norm, dot, cross

import FLOWPanel.GeometricTools.Meshes
import GeoIO

# import CUDA                               # Uncomment this to use GPU (if available)


run_name        = "wing_capped"             # Name of this run

save_path       = run_name                  # Where to save outputs
paraview        = true                      # Whether to visualize with Paraview
read_path       = joinpath(pnl.examples_path, "data") # Where to read Gmsh files from

# ----------------- SIMULATION PARAMETERS --------------------------------------
magVinf         = 1.0                      # (m/s) freestream velocity
rho             = 1.225                     # (kg/m^3) air density


# ----------------- GEOMETRY DESCRIPTION ---------------------------------------
this_dir = @__DIR__
meshfile        = joinpath(this_dir, "data", "phantom_3_mod2.msh")    # Gmsh file to read
meshfile_te     = joinpath(this_dir, "data", "phantom_3_mod2_te.msh")    # Gmsh file to read

# offset          = [0, 0, 0]                 # Offset to center the mesh
# rotation        = RotZ(-90*pi/180)*RotX(90*pi/180) # Rotation to align mesh
scaling         = 1.0                      # Factor to scale original mesh to
                                            # the approximate dimensions of the
                                            # ZEROEe BWB subscale model

spandir         = [1, 0, 0]                 # Span direction used to orient the trailing edge
flip            = false                     # Whether to flip control points against the direction of normals
                                            # NOTE: use `flip=true` if the normals
                                            #       point inside the body

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
shedding = pnl.calc_shedding(grid._nodes, pnl.grid2cells(grid), trailingedge; tolerance=0.001*b)

# Freestream vector
Vinf = magVinf*[cos(AOA*pi/180), 0, sin(AOA*pi/180)]

# Generate paneled body
if bodytype == pnl.NonLiftingBody{pnl.ConstantSource}
    body = bodytype(grid; CPoffset=(-1)^flip * 1e-14)
elseif bodytype <: pnl.RigidWakeBody
    body = bodytype(grid, shedding; CPoffset=(-1)^flip * 1e-14)
    body.Das[1] .= repeat(Vinf/magVinf, 1, size(body.Das[1], 2))
else
    error("Unsupported body type")
end
