using Pkg
Pkg.activate("/Users/ryan/Dropbox/research/projects/FLOWPanel.jl/")
using FLOWPanel
using FLOWPanel.StaticArrays
import FLOWPanel.FastMultipole as fmm
import FLOWPanel.WriteVTK
using LinearAlgebra
using BSON
include("auxiliary_functions.jl")
include("benchmark_wing_1a.jl")

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

# #--- Fast Gauss-Seidel: small system w/o expansion ---#
#
# panels = simple_grid()
# apply_freestream!(panels, SVector{3}(0,0,-1.0))
# vtk("test_simple_grid.vts", panels)
#
# # solve using LU decomposition
scheme = FLOWPanel.Scheme{FLOWPanel.DirectNeumann, FLOWPanel.FlowTangency}
# solver = LUDecomposition(panels, scheme)
# solve!(panels, solver)
# grid_2_panels_strength!(panels)
# vtk("test_simple_grid_lu.vts", panels)
#
# # solve using Fast Gauss Seidel
# FLOWPanel.set_unit_strength!(panels.panels)
# solver = FastGaussSeidel(panels, scheme; leaf_size=1)
# solve!(panels, solver)
# vtk("test_simple_grid_fgs.vts", panels)
#
# #--- Fast Gauss-Seidel: larger system w/o expansion ---#
#
# # create sphere
# n_phi = 10
# n_theta = 8
# sphere = create_sphere_unstructured(ConstantSource(); n_phi, n_theta)
# apply_freestream!(sphere, SVector{3}(1.0,0,0))
#
# # solve using LU decomposition
# solver_lu = LUDecomposition(sphere, scheme; save_lu=false)
# solve!(sphere, solver_lu)
# grid_2_panels_strength!(sphere)
#
# vtk("test_sphere_lu.vtu", sphere)
# strengths_lu = deepcopy(sphere.strengths)
# strengths_lu_vec = zeros(length(strengths_lu))
# for i in eachindex(strengths_lu)
#     strengths_lu_vec[i] = strengths_lu[i][1]
# end

# # solve using FGS (no expansions)
# # FLOWPanel.set_unit_strength!(sphere)
# FLOWPanel.reset_potential_velocity!(sphere)
# apply_freestream!(sphere, SVector{3}(1.0,0,0))
# solver_fgs = FastGaussSeidel(sphere, scheme; leaf_size=20)
# solve!(sphere, solver_fgs; tolerance=1e-10, relaxation=0.0, max_iterations=100, multipole_threshold=0.0)
#
# vtk("test_sphere_fgs.vtu", sphere)
# strengths_fgs = deepcopy(sphere.strengths)
# strengths_fgs_vec = zeros(length(strengths_fgs))
# for i in eachindex(strengths_fgs)
#     strengths_fgs_vec[i] = strengths_fgs[i][1]
# end

# # solve using FGS (with expansions)
# # FLOWPanel.set_unit_strength!(sphere)
# FLOWPanel.reset_potential_velocity!(sphere)
# apply_freestream!(sphere, SVector{3}(1.0,0,0))
# println("Building solver...")
# solver_fgs = FastGaussSeidel(sphere, scheme; leaf_size=30, expansion_order=7)
# println("solving...")
# solve!(sphere, solver_fgs; tolerance=1e-5, relaxation=0.0, max_iterations=100, multipole_threshold=0.4)
# println("Save vtk...")
# vtk("test_sphere_fgs_with_expansions.vtu", sphere)
# println("done.")
# strengths_fgs = deepcopy(sphere.strengths)
# strengths_fgs_vec = zeros(length(strengths_fgs))
# for i in eachindex(strengths_fgs)
#     strengths_fgs_vec[i] = strengths_fgs[i][1]
# end

# LU decomposition for comparison
# wing_lu = prepare_wing(;AR=1, span=1/20, nc=2, ns=1, freestream=SVector{3}(1.0,0,0))
# solver_lu = LUDecomposition(wing_lu, scheme; save_lu=false)
# solve!(wing_lu, solver_lu)
# fmm.direct!(wing_lu)
# vtk("test_wing_lu.vts", wing_lu)
# resid_lu_max, resid_lu_mean = get_residual(wing_lu)

# # don't reuse tree each iteration
# wing = prepare_wing(;nc=20, ns=100, freestream=SVector{3}(1.0,0,0))
# println("Building solver...")
# solver_fgs_wing = FastGaussSeidel(wing, scheme; leaf_size=30, expansion_order=7, multipole_threshold=0.4, reuse_tree=false)
# println("solving...")
# @time solve!(wing, solver_fgs_wing; tolerance=1e-3, relaxation=0.0, max_iterations=100)
# println("Save vtk...")
# fmm.direct!(wing)
# vtk("test_wing_fgs_a.vts", wing)
# # check residual
# resid_a = get_residual(wing)

# gmres
# wing = prepare_wing(;nc=15, ns=150, freestream=SVector{3}(1.0,0,0))
# println("Building solver...")
# solver_gmres = IterativeSolver(wing, scheme)
# println("solving...")
# @time solve!(wing, solver_gmres; tolerance=1e-3, max_iterations=20)
# println("Save vtk...")
# # fmm.direct!(wing)
# vtk("test_gmres.vts", wing)
# # check residual
# # resid_gmres = get_residual(wing)
# # resid_gmres_2 = get_residual2(wing)
# resid_gmres_3, velocity = get_residual(wing; apply_induced=true)

# mf-gmres
# println("Building solver...")
# reset_potential_velocity!(wing)
# apply_freestream!(wing, SVector{3}(1.0,0,0))
# solver_gmres_mf = MatrixFreeSolver(wing, scheme; expansion_order=10, leaf_size=50, reuse_tree=false, multipole_threshold=0.4)
# println("solving...")
# @time solve!(wing, solver_gmres_mf; tolerance=1e-3, max_iterations=20)
# println("Save vtk...")
# # fmm.direct!(wing)
# vtk("test_gmres_mf.vts", wing)
# # check residual
# # resid_gmres = get_residual(wing)
# # resid_gmres_2 = get_residual2(wing)
# resid_gmres_mf_3, velocity = get_residual(wing; apply_induced=true)

# mf-gmres
# reset_potential_velocity!(wing)
# apply_freestream!(wing, SVector{3}(1.0,0,0))
# println("Building solver...")
# solver_gmres_mf2 = MatrixFreeSolver(wing, scheme; expansion_order=10, leaf_size=50, reuse_tree=false, multipole_threshold=0.4)
# println("solving...")
# @time solve!(wing, solver_gmres_mf; tolerance=1e-3, max_iterations=20)
# println("Save vtk...")
# # fmm.direct!(wing)
# vtk("test_gmres_mf.vts", wing)
# # check residual
# # resid_gmres = get_residual(wing)
# # resid_gmres_2 = get_residual2(wing)
# # resid_gmres_mf2_3, velocity = get_residual(wing; apply_induced=true)

# reuse tree/use fmm influence matrices
# wing_fgs_b = prepare_wing(; AR=1, span=1/20, nc=20, ns=100, freestream=SVector{3}(1.0,0,0))
m=2
wing_fgs_b = prepare_wing(; nc=5*m, ns=50*m, freestream=SVector{3}(1.0,0,0))
# wing_fgs_b = prepare_wing(; nc=5, ns=1, AR=1, freestream=SVector{3}(1.0,0,0))
println("Building solver...")
solver_fgs_wing = FastGaussSeidel(wing_fgs_b, scheme; leaf_size=5, expansion_order=3, multipole_threshold=0.0, reuse_tree=true)
# solver_fgs_wing = FastGaussSeidel(wing_fgs_b, scheme; leaf_size=20, expansion_order=5, multipole_threshold=0.65, reuse_tree=true)
println("solving...")
@time solve!(wing_fgs_b, solver_fgs_wing; tolerance=1e-5, relaxation=1.4, max_iterations=50, verbose=true, save_history=true)
resid_b_max, resid_b_mean, velocity = get_residual(wing_fgs_b)
@show resid_b_max
FLOWPanel.set_unit_strength!(wing_fgs_b)
@time solve!(wing_fgs_b, solver_fgs_wing; tolerance=1e-5, relaxation=1.4, max_iterations=50, verbose=true, save_history=false)
resid_b_max, resid_b_mean, velocity = get_residual(wing_fgs_b)
@show resid_b_max
FLOWPanel.set_unit_strength!(wing_fgs_b)
@time solve!(wing_fgs_b, solver_fgs_wing; tolerance=1e-5, relaxation=1.4, max_iterations=50, verbose=true, save_history=true, error_style=:L2)
resid_b_max, resid_b_mean, velocity = get_residual(wing_fgs_b)
@show resid_b_max
FLOWPanel.set_unit_strength!(wing_fgs_b)
@time solve!(wing_fgs_b, solver_fgs_wing; tolerance=1e-5, relaxation=1.4, max_iterations=50, verbose=true, error_style=:L2,  save_history=false)
resid_b_max, resid_b_mean, velocity = get_residual(wing_fgs_b)
@show resid_b_max
println()
fmm.direct!(wing_fgs_b)
println("Save vtk...")
vtk("test_wing_fgs_b.vts", wing_fgs_b)
# check residual
@show resid_b_max

