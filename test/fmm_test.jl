# @testset "single panel expansions" begin

# single tri panel
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

# P = 13

# expansion_order = Val(P)
# branch = FastMultipole.SingleBranch(
#     1:1,
#     0,
#     1:0,
#     0,
#     SVector{3}(0.0,0.0,0.0) + displ,
#     0.4,
#     FastMultipole.initialize_expansion(P),
#     FastMultipole.initialize_expansion(P),
#     zeros(2, (P+1)*(P+1)),
#     zeros(2,4),
#     ReentrantLock(),
# )

# FastMultipole.B2M!_tripanel(panels_dipole, branch, 1:1, zeros(2,(P+1)^2), Val{P}(), FastMultipole.UniformNormalDipolePanel())

# target_potential_dipole = zeros(4)
# irregular_harmonics_dipole = FastMultipole.initialize_harmonics(P,Float64)
# multipole_expansion_dipole = branch.multipole_expansion
# FastMultipole.M2B!(target_potential_dipole, target, displ, irregular_harmonics_dipole, multipole_expansion_dipole, P)

# this_potential, _, _ = induced(target, panels_dipole.panels[1], ConstantNormalDoublet())

# @test isapprox(target_potential_dipole[1], this_potential; atol=1e-12)

# # reset and repeat for source panel
# branch.multipole_expansion .= 0.0
# source_panels = PanelArray(points, meshcells, ConstantSource())

# FastMultipole.B2M!_tripanel(source_panels, branch, 1:1, zeros(2,(P+1)^2), Val{P}(), FastMultipole.UniformSourcePanel())

# target_potential = zeros(4)
# irregular_harmonics = FastMultipole.initialize_harmonics(P,Float64)
# multipole_expansion = branch.multipole_expansion
# FastMultipole.M2B!(target_potential, target, displ, irregular_harmonics, multipole_expansion, P)

# probe = ProbeSystem([SVector{3}(target)]; scalar_potential=true, velocity=true, velocity_gradient=true)
# direct!(probe, source_panels)

# @test isapprox(target_potential[1], probe.scalar_potential[1]; atol=1e-12)

#--- tri panel ---#

points = [
    SVector{3}(0,0,0.0),
    rand(SVector{3}),
    rand(SVector{3})
]
translate_vector = rand(SVector{3,Float64})
for i in eachindex(points)
    points[i] += translate_vector
end
meshcells = [WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_TRIANGLE, SVector{3}(1,2,3))]
tri_doublet = PanelArray(points, meshcells, ConstantNormalDoublet())
tri_source = PanelArray(points, meshcells, ConstantSource())
probe = ProbeSystem([SVector{3}(4.0,4.0,4.0)]; scalar_potential=true, velocity=true, velocity_gradient=true)

# # doublet panel: direct
# direct!(probe, tri_doublet)
# potential_doublet_check = deepcopy(probe.scalar_potential[1])
# velocity_doublet_check = deepcopy(probe.velocity[1])
# gradient_doublet_check = deepcopy(probe.velocity_gradient[1])
# reset!(probe)

# # doublet panel: fmm
# fmm!(probe, tri_doublet; expansion_order=10, n_per_branch_source=1, n_per_branch_target=1, multipole_acceptance_criterion=1.0)
# potential_doublet_fmm = probe.scalar_potential[1]
# velocity_doublet_fmm = deepcopy(probe.velocity[1])
# gradient_doublet_fmm = deepcopy(probe.velocity_gradient[1])
# reset!(probe)

# source panel: direct
direct!(probe, tri_source)
potential_source_check = probe.scalar_potential[1]
velocity_source_check = deepcopy(probe.velocity[1])
gradient_source_check = deepcopy(probe.velocity_gradient[1])
reset!(probe)

# source panel: fmm
fmm!(probe, tri_source; expansion_order=20, n_per_branch_source=1, n_per_branch_target=1, multipole_acceptance_criterion=1.0)
potential_source_fmm = probe.scalar_potential[1]
velocity_source_fmm = deepcopy(probe.velocity[1])
gradient_source_fmm = deepcopy(probe.velocity_gradient[1])
reset!(probe)

# @show potential_doublet_check - potential_doublet_fmm
# @show velocity_doublet_check - velocity_doublet_fmm
# @show gradient_doublet_check - gradient_doublet_fmm
# gradient_doublet_check2 = test_gradient(probe.position[1], tri_doublet.panels[1], ConstantNormalDoublet())
# @show gradient_doublet_check2 - gradient_doublet_check
@show gradient_source_check - gradient_source_fmm


# test
# @test isapprox(potential_doublet_check, potential_doublet_fmm; atol=1e-12)
@test isapprox(potential_source_check, potential_source_fmm; atol=1e-12)
# @test isapprox(velocity_doublet_check, velocity_doublet_fmm; atol=1e-12)
@test isapprox(velocity_source_check, velocity_source_fmm; atol=1e-12)
# @test isapprox(gradient_doublet_check, gradient_doublet_fmm; atol=1e-12)
@test isapprox(gradient_source_check, gradient_source_fmm; atol=1e-12)


#--- quad panel ---#

points = [
    SVector{3}(0,0,0.0),
    rand(SVector{3}),
    rand(SVector{3}),
    rand(SVector{3})
]
points[4] = points[2]*0.9 + points[3]*1.2
translate_vector = rand(SVector{3,Float64})
for i in eachindex(points)
    points[i] += translate_vector
end
meshcells = [WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_QUAD, SVector{4}(1,2,3,4))]
quad_doublet = PanelArray(points, meshcells, ConstantNormalDoublet())
quad_source = PanelArray(points, meshcells, ConstantSource())
probe = ProbeSystem([SVector{3}(4.0,4.0,4.0)]; scalar_potential=true, velocity=true, velocity_gradient=true)

# doublet panel: direct
direct!(probe, quad_doublet)
potential_doublet_check = probe.scalar_potential[1]
velocity_doublet_check = deepcopy(probe.velocity[1])
gradient_doublet_check = deepcopy(probe.velocity_gradient[1])
reset!(probe)

# doublet panel: fmm
fmm!(probe, quad_doublet; expansion_order=14, n_per_branch_source=1, n_per_branch_target=1, multipole_acceptance_criterion=1.0)
potential_doublet_fmm = probe.scalar_potential[1]
velocity_doublet_fmm = deepcopy(probe.velocity[1])
gradient_doublet_fmm = deepcopy(probe.velocity_gradient[1])
reset!(probe)

# source panel: direct
direct!(probe, quad_source)
potential_source_check = probe.scalar_potential[1]
velocity_source_check = deepcopy(probe.velocity[1])
gradient_source_check = deepcopy(probe.velocity_gradient[1])
reset!(probe)

# source panel: fmm
direct!(probe, quad_source)
potential_source_fmm = probe.scalar_potential[1]
velocity_source_fmm = deepcopy(probe.velocity[1])
gradient_source_fmm = deepcopy(probe.velocity_gradient[1])
reset!(probe)
@show potential_doublet_check potential_doublet_fmm

# test
@test isapprox(potential_doublet_check - potential_doublet_fmm, 0.0; atol=1e-12)
@test isapprox(potential_source_check - potential_source_fmm, 0.0; atol=1e-12)
@test isapprox(velocity_doublet_check - velocity_doublet_fmm, 0.0; atol=1e-12)
@test isapprox(velocity_source_check - velocity_source_fmm, 0.0; atol=1e-12)
@test isapprox(gradient_doublet_check - gradient_doublet_fmm, 0.0; atol=1e-12)
@test isapprox(gradient_source_check - gradient_source_fmm, 0.0; atol=1e-12)


# end

@testset "fmm sphere" begin

#####
##### uniform source panels
#####
sphere_structured = create_sphere_structured(ConstantSource(); n_phi=30, n_theta=20, d_theta=5*pi/180)
sphere_unstructured = create_sphere_unstructured(ConstantSource(); n_phi=30, n_theta=26)

# compute potential directly
FastMultipole.direct!(sphere_unstructured)
FastMultipole.direct!(sphere_structured)
potential_unstructured_direct = deepcopy(sphere_unstructured.potential)
potential_structured_direct = deepcopy(sphere_structured.potential)

vtk("unstructured_sphere_source_direct", sphere_unstructured)
vtk("structured_sphere_source_direct", sphere_structured)

# reset potential
sphere_unstructured.potential .= 0.0
sphere_structured.potential .= 0.0

# compute potential using the FMM
tree_unstructured = FastMultipole.fmm!(sphere_unstructured; expansion_order=14, n_per_branch=50, multipole_acceptance_criterion=0.4)
tree_structured = FastMultipole.fmm!(sphere_structured; expansion_order=14, n_per_branch=50, multipole_acceptance_criterion=0.4)

vtk("unstructured_sphere_source_fmm", sphere_unstructured)
vtk("structured_sphere_source_fmm", sphere_structured)

potential_unstructured_fmm = deepcopy(sphere_unstructured.potential)
potential_structured_fmm = deepcopy(sphere_structured.potential)

for i in eachindex(potential_unstructured_direct)
    @test isapprox(potential_unstructured_direct[i], potential_unstructured_fmm[i]; atol=1e-12)
end

for i in eachindex(potential_structured_direct)
    @test isapprox(potential_structured_direct[i], potential_structured_fmm[i]; atol=1e-12)
end

#####
##### uniform normal dipole panels
#####
sphere_structured = create_sphere_structured(ConstantNormalDoublet(); n_phi=30, n_theta=20, d_theta=5*pi/180)
sphere_unstructured = create_sphere_unstructured(ConstantNormalDoublet(); n_phi=30, n_theta=26)

# compute potential directly
FastMultipole.direct!(sphere_unstructured)
FastMultipole.direct!(sphere_structured)
potential_unstructured_direct = deepcopy(sphere_unstructured.potential)
potential_structured_direct = deepcopy(sphere_structured.potential)

vtk("unstructured_sphere_dipole_direct", sphere_unstructured)
vtk("structured_sphere_dipole_direct", sphere_structured)

# reset potential
sphere_unstructured.potential .= 0.0
sphere_structured.potential .= 0.0

# compute potential using the FMM
tree_unstructured = FastMultipole.fmm!(sphere_unstructured; expansion_order=15, n_per_branch=50, multipole_acceptance_criterion=0.4)
tree_structured = FastMultipole.fmm!(sphere_structured; expansion_order=15, n_per_branch=50, multipole_acceptance_criterion=0.4)

vtk("unstructured_sphere_dipole_fmm", sphere_unstructured)
vtk("structured_sphere_dipole_fmm", sphere_structured)

potential_unstructured_fmm = deepcopy(sphere_unstructured.potential)
potential_structured_fmm = deepcopy(sphere_structured.potential)

for i in eachindex(potential_unstructured_direct)
    @test isapprox(potential_unstructured_direct[i], potential_unstructured_fmm[i]; atol=1e-12)
end

for i in eachindex(potential_structured_direct)
    @test isapprox(potential_structured_direct[i], potential_structured_fmm[i]; atol=1e-12)
end

end
