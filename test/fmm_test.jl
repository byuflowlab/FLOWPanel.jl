@testset "fmm" begin

sphere_structured = create_sphere_structured(ConstantSource(); n_phi=30, n_theta=20, d_theta=5*pi/180)
sphere_unstructured = create_sphere_unstructured(ConstantSource(); n_phi=30, n_theta=26)

# compute potential directly
FLOWFMM.direct!(sphere_unstructured)
FLOWFMM.direct!(sphere_structured)
potential_unstructured_direct = deepcopy(sphere_unstructured.potential)
potential_structured_direct = deepcopy(sphere_structured.potential)

vtk("unstructured_sphere_direct", sphere_unstructured)
vtk("structured_sphere_direct", sphere_structured)

# reset potential
sphere_unstructured.potential .= 0.0
sphere_structured.potential .= 0.0

# compute potential using the FMM
tree_unstructured = FLOWFMM.fmm!(sphere_unstructured; expansion_order=14, n_per_branch=50, theta=0.4)
tree_structured = FLOWFMM.fmm!(sphere_structured; expansion_order=14, n_per_branch=50, theta=0.4)

vtk("unstructured_sphere_fmm", sphere_unstructured)
vtk("structured_sphere_fmm", sphere_structured)

potential_unstructured_fmm = deepcopy(sphere_unstructured.potential)
potential_structured_fmm = deepcopy(sphere_structured.potential)

for i in eachindex(potential_unstructured_direct)
    @test isapprox(potential_unstructured_direct[i], potential_unstructured_fmm[i]; atol=1e-12)
end

for i in eachindex(potential_structured_direct)
    @test isapprox(potential_structured_direct[i], potential_structured_fmm[i]; atol=1e-12)
end

end