@testset "single panel expansions" begin

# single tri panel
displ = rand(SVector{3,Float64})
points = [
    SVector{3}(0.0,0.0,-0.2) + displ,
    SVector{3}(2,1,0.3) + displ,
    SVector{3}(-1,2,0.1) + displ
]
meshcells = [WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_TRIANGLE, SVector{3}(1,2,3))]
target = [5.5, -4.2, 5.3] + displ
centroid = (points[1] + points[2] + points[3])/3
normal = cross(points[2]-points[1],points[3]-points[1])
area = norm(normal)/2
normal /= area*2
strength = 1.0
panels_dipole = PanelArray(points, meshcells, ConstantNormalDoublet())

P = 13

expansion_order = Val(P)
branch = FLOWFMM.SingleBranch(
    1:1,
    0,
    1:0,
    0,
    SVector{3}(0.0,0.0,0.0) + displ,
    0.4,
    FLOWFMM.initialize_expansion(P),
    FLOWFMM.initialize_expansion(P),
    zeros(2, (P+1)*(P+1)),
    zeros(2,4),
    ReentrantLock(),
)

FLOWFMM.B2M!_tripanel(panels_dipole, branch, 1:1, zeros(2,(P+1)^2), Val{P}(), FLOWFMM.UniformNormalDipolePanel())

target_potential_dipole = zeros(4)
irregular_harmonics_dipole = FLOWFMM.initialize_harmonics(P,Float64)
multipole_expansion_dipole = branch.multipole_expansion
FLOWFMM.M2B!(target_potential_dipole, target, displ, irregular_harmonics_dipole, multipole_expansion_dipole, P)

this_potential, _, _ = induced(target, panels_dipole.panels[1], ConstantNormalDoublet(); toggle_potential=true, toggle_velocity=false, toggle_hessian=false)

@test isapprox(target_potential_dipole[1], this_potential; atol=1e-12)

# reset and repeat for source panel
branch.multipole_expansion .= 0.0
source_panels = PanelArray(points, meshcells, ConstantSource())

FLOWFMM.B2M!_tripanel(source_panels, branch, 1:1, zeros(2,(P+1)^2), Val{P}(), FLOWFMM.UniformSourcePanel())

target_potential = zeros(4)
irregular_harmonics = FLOWFMM.initialize_harmonics(P,Float64)
multipole_expansion = branch.multipole_expansion
FLOWFMM.M2B!(target_potential, target, displ, irregular_harmonics, multipole_expansion, P)

@test isapprox(target_potential[1], 0.022849841377239076; atol=1e-12)

end

@testset "fmm sphere" begin

#####
##### uniform source panels
#####
sphere_structured = create_sphere_structured(ConstantSource(); n_phi=30, n_theta=20, d_theta=5*pi/180)
sphere_unstructured = create_sphere_unstructured(ConstantSource(); n_phi=30, n_theta=26)

# compute potential directly
FLOWFMM.direct!(sphere_unstructured)
FLOWFMM.direct!(sphere_structured)
potential_unstructured_direct = deepcopy(sphere_unstructured.potential)
potential_structured_direct = deepcopy(sphere_structured.potential)

vtk("unstructured_sphere_source_direct", sphere_unstructured)
vtk("structured_sphere_source_direct", sphere_structured)

# reset potential
sphere_unstructured.potential .= 0.0
sphere_structured.potential .= 0.0

# compute potential using the FMM
tree_unstructured = FLOWFMM.fmm!(sphere_unstructured; expansion_order=14, n_per_branch=50, theta=0.4)
tree_structured = FLOWFMM.fmm!(sphere_structured; expansion_order=14, n_per_branch=50, theta=0.4)

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
FLOWFMM.direct!(sphere_unstructured)
FLOWFMM.direct!(sphere_structured)
potential_unstructured_direct = deepcopy(sphere_unstructured.potential)
potential_structured_direct = deepcopy(sphere_structured.potential)

vtk("unstructured_sphere_dipole_direct", sphere_unstructured)
vtk("structured_sphere_dipole_direct", sphere_structured)

# reset potential
sphere_unstructured.potential .= 0.0
sphere_structured.potential .= 0.0

# compute potential using the FMM
tree_unstructured = FLOWFMM.fmm!(sphere_unstructured; expansion_order=15, n_per_branch=50, theta=0.4)
tree_structured = FLOWFMM.fmm!(sphere_structured; expansion_order=15, n_per_branch=50, theta=0.4)

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
