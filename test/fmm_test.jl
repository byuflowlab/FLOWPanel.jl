@testset "single panel expansions" begin

#--- tri panel ---#

points = [
    SVector{3}(0,0,0.0),
    rand(SVector{3}),
    rand(SVector{3}),
]
translate_vector = rand(SVector{3,Float64})
for i in eachindex(points)
    points[i] += translate_vector
end
for point in points[1:3]
    push!(points, -point)
end
# points = SVector{3, Float64}[[0.044818005017491114, 0.933353287277165, 0.5805599818745412], [0.5660318005528742, 1.5201600447305135, 1.4714386799673222], [0.23572470404325396, 1.4590156788192123, 0.9711482573058853], [-0.044818005017491114, -0.933353287277165, -0.5805599818745412], [-0.5660318005528742, -1.5201600447305135, -1.4714386799673222], [-0.23572470404325396, -1.4590156788192123, -0.9711482573058853]]
dist = 4.0
meshcells = [WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_TRIANGLE, SVector{3}(1,2,3)), WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_TRIANGLE, SVector{3}(4,5,6))]
tri_doublet = PanelArray(points, [meshcells[1]], ConstantNormalDoublet())
tri_source = PanelArray(points, meshcells, ConstantSource())
probe = ProbeSystem([SVector{3}(1.0,1.0,1.0) * dist]; scalar_potential=true, velocity=true, velocity_gradient=true)

# doublet panel: direct
direct!(probe, tri_doublet)
potential_doublet_check = probe.scalar_potential[1]
velocity_doublet_check = deepcopy(probe.velocity[1])
gradient_doublet_check = deepcopy(probe.velocity_gradient[1])
reset!(probe)

# doublet panel: fmm
st_now, tt_now = fmm!(probe, tri_doublet; expansion_order=20, leaf_size_source=1, leaf_size_target=1, multipole_threshold=1.0)
potential_doublet_fmm = deepcopy(probe.scalar_potential[1])
velocity_doublet_fmm = deepcopy(probe.velocity[1])
gradient_doublet_fmm = deepcopy(probe.velocity_gradient[1])
reset!(probe)

# source panel: direct
direct!(probe, tri_source)
potential_source_check = deepcopy(probe.scalar_potential[1])
velocity_source_check = deepcopy(probe.velocity[1])
gradient_source_check = deepcopy(probe.velocity_gradient[1])
reset!(probe)

# check again
pcheck, vcheck, gcheck = induced(probe.position[1], tri_source.panels[1], ConstantSource())

gcheck2 = ForwardDiff.jacobian((x) -> induced(x, tri_source.panels[1], ConstantSource())[2], probe.position[1])

# source panel: fmm
source_tree, target_tree = fmm!(probe, tri_source; expansion_order=20, leaf_size_source=2, leaf_size_target=2, multipole_threshold=1.0)
potential_source_fmm = deepcopy(probe.scalar_potential[1])
velocity_source_fmm = deepcopy(probe.velocity[1])
gradient_source_fmm = deepcopy(probe.velocity_gradient[1])
reset!(probe)

# test
@test isapprox(potential_doublet_check, potential_doublet_fmm; atol=1e-12)
@test isapprox(potential_source_check, potential_source_fmm; atol=2e-12)
@test isapprox(velocity_doublet_check, velocity_doublet_fmm; atol=1e-12)
@test isapprox(velocity_source_check, velocity_source_fmm; atol=1e-12)
@test isapprox(gradient_doublet_check, gradient_doublet_fmm; atol=1e-12)
@test isapprox(gradient_source_check, gradient_source_fmm; atol=1e-12)

#--- quad panel ---#

points = [
    SVector{3}(0,0,0.0),
    rand(SVector{3}),
    rand(SVector{3}),
    rand(SVector{3})
]
points[3] = points[2]*0.9 + points[4]*1.2
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
fmm!(probe, quad_doublet; expansion_order=20, leaf_size_source=1, leaf_size_target=1, multipole_threshold=1.0)
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

# test
@test isapprox(potential_doublet_check, potential_doublet_fmm; atol=1e-12)
@test isapprox(potential_source_check, potential_source_fmm; atol=1e-12)
@test isapprox(velocity_doublet_check, velocity_doublet_fmm; atol=1e-12)
@test isapprox(velocity_source_check, velocity_source_fmm; atol=1e-12)
@test isapprox(gradient_doublet_check, gradient_doublet_fmm; atol=1e-12)
@test isapprox(gradient_source_check, gradient_source_fmm; atol=1e-12)

end

#####
##### debug VortexLattice
#####

points = [
    SVector{3}(0.25,0,0),
    SVector{3}(1.0,0,0),
    SVector{3}(1.0,1.0,0),
    SVector{3}(0.25,1.0,0)
]
meshcells = [WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_QUAD, SVector{4}(1,2,3,4))]
quad_doublet = PanelArray(points, meshcells, ConstantNormalDoublet())
probe = ProbeSystem([SVector{3}(0.5, 0.5, 2.0)]; scalar_potential=true, velocity=true, velocity_gradient=true)

# doublet panel: direct
direct!(probe, quad_doublet)
potential_doublet_check = probe.scalar_potential[1]
velocity_doublet_check = deepcopy(probe.velocity[1])
gradient_doublet_check = deepcopy(probe.velocity_gradient[1])
reset!(probe)

# doublet panel: fmm
fmm!(probe, quad_doublet; expansion_order=20, leaf_size_source=1, leaf_size_target=1, multipole_threshold=1.0)
potential_doublet_fmm = probe.scalar_potential[1]
velocity_doublet_fmm = deepcopy(probe.velocity[1])
gradient_doublet_fmm = deepcopy(probe.velocity_gradient[1])
reset!(probe)

@show velocity_doublet_fmm velocity_doublet_check


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
tree_unstructured = FastMultipole.fmm!(sphere_unstructured; expansion_order=14, leaf_size=50, multipole_threshold=0.4)
tree_structured = FastMultipole.fmm!(sphere_structured; expansion_order=14, leaf_size=50, multipole_threshold=0.4)

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
tree_unstructured = FastMultipole.fmm!(sphere_unstructured; expansion_order=15, leaf_size=50, multipole_threshold=0.4)
tree_structured = FastMultipole.fmm!(sphere_structured; expansion_order=15, leaf_size=50, multipole_threshold=0.4)

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
