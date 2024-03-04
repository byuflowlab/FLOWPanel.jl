# @testset "solver: lu decomposition" begin
    # using LinearAlgebra

#####
##### source panels
#####

# create panels
sphere_structured = create_sphere_structured(ConstantSource(); n_phi=30, n_theta=22, d_theta=5*pi/180)
sphere_unstructured = create_sphere_unstructured(ConstantSource(); n_phi=30, n_theta=27)

# apply freestream
apply_freestream!(sphere_structured, SVector{3}(1.0,0,0))
apply_freestream!(sphere_unstructured, SVector{3}(1.0,0,0))

# solve panels
solver_structured = LUDecomposition(sphere_structured, FLOWPanel.Scheme{FLOWPanel.DirectNeumann,FlowTangency})
solve!(sphere_structured, solver_structured)
grid_2_panels_strength!(sphere_structured)
tree_structured = FLOWFMM.fmm!(sphere_structured; expansion_order=14, n_per_branch=50, theta=0.4)
# vtk("structured_source_sphere_solved", sphere_structured)

for (velocity,panel) in zip(sphere_structured.velocity, sphere_structured.panels)
    @test isapprox(dot(panel.normal, velocity), 0.0; atol=1e-6)
end

solver_unstructured = LUDecomposition(sphere_unstructured, FLOWPanel.Scheme{FLOWPanel.DirectNeumann,FlowTangency})
solve!(sphere_unstructured, solver_unstructured)
grid_2_panels_strength!(sphere_unstructured)
tree_unstructured = FLOWFMM.fmm!(sphere_unstructured; expansion_order=14, n_per_branch=50, theta=0.4)
# vtk("unstructured_source_sphere_solved", sphere_unstructured)

for (velocity,panel) in zip(sphere_unstructured.velocity, sphere_unstructured.panels)
    @test isapprox(dot(panel.normal, velocity), 0.0; atol=1e-12)
end

# #####
# ##### pure dipole panels
# #####

# # create panels
# dipole_sphere_structured = create_sphere_structured(ConstantNormalDoublet(); n_phi=30, n_theta=22, d_theta=5*pi/180)
# dipole_sphere_unstructured = create_sphere_unstructured(ConstantNormalDoublet(); n_phi=30, n_theta=27)

# # apply freestream
# apply_freestream!(dipole_sphere_structured, SVector{3}(1.0,0,0))
# apply_freestream!(dipole_sphere_unstructured, SVector{3}(1.0,0,0))

# # solve panels
# dipole_solver_structured = LUDecomposition(dipole_sphere_structured, FLOWPanel.Scheme{FLOWPanel.DirectNeumann,FlowTangency})
# solve!(dipole_sphere_structured, dipole_solver_structured, SVector{3,Float64}(0,0,0))
# grid_2_panels_strength!(dipole_sphere_structured)
# tree_structured = FLOWFMM.fmm!(dipole_sphere_structured; expansion_order=14, n_per_branch=50, theta=0.4)
# vtk("structured_dipole_sphere_solved", dipole_sphere_structured)

# dipole_solver_unstructured = LUDecomposition(dipole_sphere_unstructured, FLOWPanel.Scheme{FLOWPanel.DirectNeumann,FlowTangency})
# solve!(dipole_sphere_unstructured, dipole_solver_unstructured, SVector{3,Float64}(0,0,0))
# grid_2_panels_strength!(dipole_sphere_unstructured)
# tree_unstructured = FLOWFMM.fmm!(dipole_sphere_unstructured; expansion_order=14, n_per_branch=50, theta=0.4)
# vtk("unstructured_dipole_sphere_solved", dipole_sphere_unstructured)

# errs = Float64[]
# for (velocity,panel) in zip(dipole_sphere_unstructured.velocity, dipole_sphere_unstructured.panels)
#     push!(errs, dot(panel.normal, velocity))
# end




# displ = rand(SVector{3,Float64})
# points = [
#     SVector{3}(0.0,0.0,-0.2) + displ,
#     SVector{3}(2,1,0.3) + displ,
#     SVector{3}(-1,2,0.1) + displ
# ]
# meshcells = [WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_TRIANGLE, SVector{3}(1,2,3))]
# target = [5.5, -4.2, 5.3] + displ
# centroid = (points[1] + points[2] + points[3])/3
# normal = cross(points[2]-points[1],points[3]-points[1])
# area = norm(normal)/2
# normal /= area*2
# strength = 1.0
# panels_dipole = PanelArray(points, meshcells, ConstantNormalDoublet())
# this_potential_dipole, this_velocity_dipole, _ = induced(panels_dipole.panels[1].control_point, panels_dipole.panels[1], ConstantNormalDoublet(); toggle_potential=true, toggle_velocity=true, toggle_hessian=false)

# # component in the normal/tangential direction
# vn_dipole = dot(this_velocity_dipole, normal) * normal
# vt_dipole = this_velocity_dipole - vn_dipole

# panels_source = PanelArray(points, meshcells, ConstantSource())
# this_potential_source, this_velocity_source, _ = induced(panels_source.panels[1].control_point, panels_source.panels[1], ConstantSource(); toggle_potential=true, toggle_velocity=true, toggle_hessian=false)

# # component in the normal/tangential direction
# vn_source = dot(this_velocity_source, normal) * normal
# vt_source = this_velocity_source - vn_source

# end

# @testset "solver: fast jacobi" begin

# # create panels
# sphere_structured = create_sphere_structured(ConstantSource(); n_phi=30, n_theta=22, d_theta=5*pi/180)
# sphere_unstructured = create_sphere_unstructured(ConstantSource(); n_phi=30, n_theta=27)

# # apply freestream
# apply_freestream!(sphere_structured, SVector{3}(1.0,0,0))
# apply_freestream!(sphere_unstructured, SVector{3}(1.0,0,0))

# # solve panels
# solver_structured = FastJacobi(sphere_structured, FLOWPanel.Scheme{FLOWPanel.DirectNeumann,FlowTangency}; n_per_branch=16, max_iter=200, epsilon=1e-6)
# solve!(sphere_structured, solver_structured; verbose=true, update_influence_matrices=false)
# vtk("structured_source_sphere_solved_fastjacobi", sphere_structured)

# for (velocity,panel) in zip(sphere_structured.velocity, sphere_structured.panels)
#     @test isapprox(dot(panel.normal, velocity), 0.0; atol=1e-6)
# end

# structured_solution_saved = deepcopy(solver_structured.strength_history[:,1])

# # perturb strength solution slightly
# solver_structured.previous_dt .= [0.0,-0.1,-0.2]
# solver_structured.strength_history[:,1] .+= solver_structured.strength_history[:,1] .* 0.05
# solver_structured.strength_history[:,2] .+= solver_structured.strength_history[:,2] .* 0.1
# solver_structured.strength_history[:,3] .+= solver_structured.strength_history[:,2] .* 0.3
# for i in 1:length(sphere_structured.velocity)
#     sphere_structured.velocity[i] = zero(eltype(sphere_structured.velocity))
# end
# sphere_structured.potential .= 0.0
# apply_freestream!(sphere_structured, SVector{3}(1.0,0,0))
# solve!(sphere_structured, solver_structured, 0.1)
# vtk("structured_source_sphere_solved_fastjacobi_taylorguess", sphere_structured)

# for (velocity,panel) in zip(sphere_structured.velocity, sphere_structured.panels)
#     @test isapprox(dot(panel.normal, velocity), 0.0; atol=1e-6)
# end

# solver_unstructured = FastJacobi(sphere_unstructured, FLOWPanel.Scheme{FLOWPanel.DirectNeumann,FlowTangency}; n_per_branch=15, max_iter=200, epsilon=1e-6)
# solve!(sphere_unstructured, solver_unstructured; verbose=false)
# vtk("unstructured_source_sphere_solved_fastjacobi", sphere_unstructured)

# for (velocity,panel) in zip(sphere_unstructured.velocity, sphere_unstructured.panels)
#     @test isapprox(dot(panel.normal, velocity), 0.0; atol=1e-6)
# end

# end

# @testset "taylor series guess" begin

# dummy_strength_history = zeros(10,3)
# previous_dt = MVector{3,Float64}(0.0, -0.1, -0.3)
# dt = 0.1
# fit_order = Val(2)
# for i in axes(dummy_strength_history,1)
#     for j in eachindex(previous_dt)
#         dummy_strength_history[i,j] = sin(i*0.13 + previous_dt[j])
#     end
# end

# test_next = zeros(size(dummy_strength_history,1))
# for i in axes(dummy_strength_history,1)
#     test_next[i] = sin(i*0.13 + dt)
# end

# FLOWPanel.guess_strengths!(dummy_strength_history, previous_dt, fit_order, dt)

# for i in eachindex(test_next)
#     @test isapprox(test_next[i], dummy_strength_history[i,end], atol=5e-1)
# end

# # using PythonPlot
# # fig = figure("test_taylor")
# # clf()
# # fig.add_subplot(111,xlabel="x",ylabel="y")
# # ax = fig.get_axes()[0]
# # xs = range(-0.3, stop=0.5, length=50)
# # ys_true = [sin.(i*0.13 .+ xs) for i in 1:10]
# # ax.plot(xs, ys_true[1], "-k")
# # ax.plot(previous_dt, sin.(0.13 .+ previous_dt), "ko")
# # ax.plot(dt, dummy_strength_history[1,end], "k>")
# # ax.plot(xs, ys_true[5], "-b")
# # ax.plot(previous_dt, sin.(0.13*5 .+ previous_dt), "bo")
# # ax.plot(dt, dummy_strength_history[5,end], "b>")
# # ax.plot(xs, ys_true[10], "-g")
# # ax.plot(previous_dt, sin.(0.13*10 .+ previous_dt), "go")
# # ax.plot(dt, dummy_strength_history[10,end], "g>")

# end

# @testset "fast jacobi: intelligent guess" begin

# # create panels
# sphere_structured = create_sphere_structured(ConstantSource(); n_phi=40, n_theta=30, d_theta=5*pi/180)

# # instantiate solver
# solver_structured = FastJacobi(sphere_structured, FLOWPanel.Scheme{FLOWPanel.DirectNeumann,FlowTangency}; n_per_branch=16, max_iter=200, epsilon=1e-7, fit_order=1, n_previous_steps=2)

# # unsteady simulation
# nsteps = 20
# ts_solve = zeros(nsteps)
# for (i,fs) in enumerate(range(1.0, stop=2, length=nsteps))
#     # println("Begin Vinf=$fs:")

#     # apply freestream
#     apply_freestream!(sphere_structured, SVector{3}(fs,0,0))
    
#     # solve panels
#     tsolve = @elapsed solve!(sphere_structured, solver_structured, 0.001; verbose=true, update_influence_matrices=false)
#     ts_solve[i] = tsolve

#     # check boundary condition
#     for (velocity,panel) in zip(sphere_structured.velocity, sphere_structured.panels)
#         @test isapprox(dot(panel.normal, velocity), 0.0; atol=1e-6)
#     end

#     # save vtk
#     # vtk("structured_sphere_fast_solve.$i.vts", sphere_structured)

#     # reset for next time
#     reset_potential_velocity!(sphere_structured)
# end

# @test ts_solve[1] > sum(ts_solve[2:end]) / (length(ts_solve)-1) * 2

# # using PythonPlot

# # cla()
# # plot(1:100,ts_solve,"-x")
# # xlabel("timestep")
# # ylabel("time cost, seconds")
# # suptitle("linear extrapolation")

# end

# @testset "solver: fast gauss-seidel" begin

# create panels
sphere_structured = create_sphere_structured(ConstantSource(); n_phi=40, n_theta=48, d_theta=5*pi/180)
sphere_structured_lu = create_sphere_structured(ConstantSource(); n_phi=40, n_theta=48, d_theta=5*pi/180)
sphere_unstructured = create_sphere_unstructured(ConstantSource(); n_phi=30, n_theta=27)

# apply freestream
apply_freestream!(sphere_structured, SVector{3}(1.0,0,0))
apply_freestream!(sphere_structured_lu, SVector{3}(1.0,0,0))
apply_freestream!(sphere_unstructured, SVector{3}(1.0,0,0))

# solve panels
solver_structured, tree = FastGaussSeidel(sphere_structured, FLOWPanel.Scheme{FLOWPanel.DirectNeumann,FlowTangency}; n_per_branch=16, theta=0.5, max_iter=50, epsilon=1e-10)
solver_lu = LUDecomposition(sphere_structured_lu, FLOWPanel.Scheme{FLOWPanel.DirectNeumann,FlowTangency})
solve!(sphere_structured_lu, solver_lu)
vtk("structured_source_sphere_solved_ludebug", sphere_structured)

solve!(sphere_structured, solver_structured; verbose=true, update_influence_matrices=false)
vtk("structured_source_sphere_solved_fastgaussseidel", sphere_structured)
#=
for (velocity,panel) in zip(sphere_structured.velocity, sphere_structured.panels)
    @test isapprox(dot(panel.normal, velocity), 0.0; atol=1e-6)
end

structured_solution_saved = deepcopy(solver_structured.strength_history[:,1])

# perturb strength solution slightly
solver_structured.previous_dt .= [0.0,-0.1,-0.2]
solver_structured.strength_history[:,1] .+= solver_structured.strength_history[:,1] .* 0.05
solver_structured.strength_history[:,2] .+= solver_structured.strength_history[:,2] .* 0.1
solver_structured.strength_history[:,3] .+= solver_structured.strength_history[:,2] .* 0.3
for i in 1:length(sphere_structured.velocity)
    sphere_structured.velocity[i] = zero(eltype(sphere_structured.velocity))
end
sphere_structured.potential .= 0.0
apply_freestream!(sphere_structured, SVector{3}(1.0,0,0))
solve!(sphere_structured, solver_structured, 0.1)
vtk("structured_source_sphere_solved_fastgaussseidel_taylorguess", sphere_structured)

for (velocity,panel) in zip(sphere_structured.velocity, sphere_structured.panels)
    @test isapprox(dot(panel.normal, velocity), 0.0; atol=1e-6)
end

solver_unstructured = FastJacobi(sphere_unstructured, FLOWPanel.Scheme{FLOWPanel.DirectNeumann,FlowTangency}; n_per_branch=15, max_iter=200, epsilon=1e-6)
solve!(sphere_unstructured, solver_unstructured; verbose=false)
vtk("unstructured_source_sphere_solved_fastgaussseidel", sphere_unstructured)

for (velocity,panel) in zip(sphere_unstructured.velocity, sphere_unstructured.panels)
    @test isapprox(dot(panel.normal, velocity), 0.0; atol=1e-6)
end

# end
=#