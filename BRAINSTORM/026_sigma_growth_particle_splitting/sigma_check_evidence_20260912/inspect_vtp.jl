using ReadVTK
f = "/Users/ryan/Dropbox/research/projects/FLOWPanel.jl/data/rotor_hover_pressure_comparison/rotor_hover_pressure_comparison_wake1_particles/rotor_hover_pressure_comparison_wake1_particles.233.vtp"
vtk = VTKFile(f)
pd = get_point_data(vtk)
println("arrays: ", keys(pd))
println("npoints: ", vtk.n_points)
