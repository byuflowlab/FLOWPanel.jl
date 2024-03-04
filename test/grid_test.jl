@testset "create grids" begin

sphere_structured = create_sphere_structured(ConstantSource(); n_phi=20, n_theta=10, d_theta=5*pi/180)
sphere_unstructured = create_sphere_unstructured(ConstantSource(); n_phi=20, n_theta=18)

#####
##### apply freestream
#####
freestream = SVector{3}(1.0,0,0)
apply_freestream!(sphere_structured, freestream)
apply_freestream!(sphere_unstructured, freestream)

vtk("unstructured_sphere", sphere_unstructured)
vtk("structured_sphere", sphere_structured)

end
