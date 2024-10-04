using Pkg
this_dir = @__DIR__
Pkg.activate(normpath(joinpath(this_dir,"..")))

using FLOWPanel
using FLOWPanel.FastMultipole
using StaticArrays
using BenchmarkTools

function create_checkerboard_structured(kernel; lx=1.0, nx=3, invert_normals=false)
    corner_grid = zeros(SVector{3,Float64}, nx+1, nx+1, 1)
    y = -lx/2
    for iy in 1:nx+1
        x = -lx/2
        for ix in 1:nx+1
            corner_grid[ix,iy,1] = SVector{3}(x,y,0.0)
            x += lx/nx
        end
        y += lx/nx
    end

    return PanelArray(corner_grid, kernel; invert_normals)
end

# input arguments
lx, nx = 1.0, 3

# create probes
probe_positions = [
    SVector{3}(-lx/2,-lx/2,0),
    SVector{3}(-lx/2,0,0),
    SVector{3}(-lx/2,lx,0),
    SVector{3}(-lx/2+lx/3, 0, 0),
    SVector{3}(-lx/2+lx/3, -lx/2+lx/3, 0)
]
probes_doublet = ProbeSystem(probe_positions; velocity=true)
probes_vortexring = ProbeSystem(probe_positions; velocity=true)

# doublet panels
checkerboard_doublet = create_checkerboard_structured(ConstantNormalDoublet())
test_doublet!(checkerboard_doublet) = direct!(checkerboard_doublet)

test_doublet!(checkerboard_doublet)
vtk("checkerboard_doublet", checkerboard_doublet)
v_doublet = deepcopy(checkerboard_doublet.velocity)

direct!(probes_doublet, checkerboard_doublet)
FastMultipole.visualize_bodies("doublet_on_probes", probes_doublet)

# @btime test_doublet!($checkerboard_doublet)

# vortex ring panels
checkerboard_vortexring = create_checkerboard_structured(VortexRing())
test_vortexring!(checkerboard_vortexring) = direct!(checkerboard_vortexring)

test_vortexring!(checkerboard_vortexring)
vtk("checkerboard_vortexring", checkerboard_vortexring)
v_vortexring = deepcopy(checkerboard_vortexring.velocity)

direct!(probes_vortexring, checkerboard_vortexring)
FastMultipole.visualize_bodies("vortexring_on_probes", probes_vortexring)

# @btime test_vortexring!($checkerboard_vortexring)

for iv in eachindex(v_doublet)
    for i in 1:3
        @show isapprox(v_vortexring[iv][i], v_doublet[iv][i]; atol=1e-12)
    end
end