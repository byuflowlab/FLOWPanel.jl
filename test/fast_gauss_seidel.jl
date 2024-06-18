using Pkg
Pkg.activate("/Users/ryan/Dropbox/research/projects/FLOWPanel.jl/")
using FLOWPanel
using FLOWPanel.StaticArrays
import FLOWPanel.FastMultipole as fmm
import FLOWPanel.WriteVTK
using LinearAlgebra
using BSON
include("auxiliary_functions.jl")
include("benchmark_wing.jl")

function simple_grid(;
        kernel=ConstantSource(),
        invert_normals=false
    )
    grid = zeros(3,3,2,1)
    xs = range(0,1,length=3)
    ys = range(0,1,length=2)
    zs = [0,0,1.0]

    for j in 1:2
        y = ys[j]
        for i in 1:3
            x = xs[i]
            z = zs[i]
            grid[:,i,j,1] .= SVector{3}(x,y,z)
        end
    end

    return PanelArray(grid, kernel; invert_normals)
end

#--- Fast Gauss-Seidel: small system w/o expansion ---#

panels = simple_grid()
apply_freestream!(panels, SVector{3}(0,0,-1.0))
vtk("test_simple_grid.vts", panels)

# solve using LU decomposition
scheme = FLOWPanel.Scheme{FLOWPanel.DirectNeumann, FLOWPanel.FlowTangency}
solver = LUDecomposition(panels, scheme)
solve!(panels, solver)
grid_2_panels_strength!(panels)
vtk("test_simple_grid_lu.vts", panels)

# solve using Fast Gauss Seidel
FLOWPanel.set_unit_strength!(panels.panels)
solver = FastGaussSeidel(panels, scheme; leaf_size=1)
solve!(panels, solver)
vtk("test_simple_grid_fgs.vts", panels)

#--- Fast Gauss-Seidel: larger system w/o expansion ---#

# create sphere
n_phi = 10
n_theta = 7
sphere = create_sphere_unstructured(ConstantSource(); n_phi, n_theta)
apply_freestream!(sphere, SVector{3}(1.0,0,0))

# solve using LU decomposition
solver_lu = LUDecomposition(sphere, scheme)
solve!(sphere, solver_lu)
grid_2_panels_strength!(sphere)

vtk("test_sphere_lu.vtu", sphere)
strengths_lu = deepcopy(sphere.strengths)
strengths_lu_vec = zeros(length(strengths_lu))
for i in eachindex(strengths_lu)
    strengths_lu_vec[i] = strengths_lu[i][1]
end

# solve using FGS (no expansions)
FLOWPanel.set_unit_strength!(sphere)
FLOWPanel.reset_potential_velocity!(sphere)
apply_freestream!(sphere, SVector{3}(1.0,0,0))
solver_fgs = FastGaussSeidel(sphere, scheme; leaf_size=20)
solve!(sphere, solver_fgs; relaxation=0.0, max_iterations=100, multipole_threshold=0.0)

vtk("test_sphere_fgs.vtu", sphere)
strengths_fgs = deepcopy(sphere.strengths)
strengths_fgs_vec = zeros(length(strengths_fgs))
for i in eachindex(strengths_fgs)
    strengths_fgs_vec[i] = strengths_fgs[i][1]
end

# solve using FGS (with expansions)
FLOWPanel.set_unit_strength!(sphere)
FLOWPanel.reset_potential_velocity!(sphere)
apply_freestream!(sphere, SVector{3}(1.0,0,0))
solver_fgs = FastGaussSeidel(sphere, scheme; leaf_size=2)
solve!(sphere, solver_fgs; relaxation=0.0, max_iterations=100, multipole_threshold=0.5)

vtk("test_sphere_fgs_with_expansions.vtu", sphere)
strengths_fgs = deepcopy(sphere.strengths)
strengths_fgs_vec = zeros(length(strengths_fgs))
for i in eachindex(strengths_fgs)
    strengths_fgs_vec[i] = strengths_fgs[i][1]
end

