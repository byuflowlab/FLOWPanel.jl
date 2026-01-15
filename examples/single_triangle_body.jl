# Example: construct a NonLiftingBody with a single triangular panel
# Vertices: [0,0,0], [1,0,0], [0,1,0]

using FLOWPanel
using FLOWPanel.gt.Meshes
import FLOWPanel.PyPlot as plt

function make_single_triangle_body()
    pnl = FLOWPanel

    # Create vertices for a single triangle
    vertices = [
        Point(0.0, 0.0, 0.0),
        Point(1.0, 0.0, 0.0),
        Point(0.0, 1.0, 0.0)
    ]

    # Create a single triangle using Meshes.jl
    triangle = [connect((1,2,3))]
    mesh = SimpleMesh(vertices, triangle)

    # Wrap it as a GridTriangleSurface
    grid = FLOWPanel.gt.GridTriangleSurface(mesh)

    # Construct NonLiftingBody passing `cells` explicitly to avoid calling gt-based grid2cells
    body = pnl.NonLiftingBody{pnl.ConstantSource}(grid)

    return body
end

# generate body
body = make_single_triangle_body()

# set unit strength
body.strength[1,1] = 1.0

# compute velocity at a line of points above the triangle
npoints = 100
dz = range(0.0, stop=10.0, length=npoints)
points = hcat([ [0.33; 0.33; z] for z in dz ]...)  # points above centroid of triangle
points = hcat([[0.33; 0.33; 0.0] .+ [z; z; z] for z in dz ]...)  # points in y from the centroid

# set "solved" flag
FLOWPanel._solvedflag(body, true)

# direct backend
backend = pnl.DirectBackend()
out_direct = zeros(3, npoints)
pnl.Uind!(body, points, out_direct, backend)

# fast multipole backend
backend_fmm = pnl.FastMultipoleBackend(
                                    expansion_order=20,
                                    multipole_acceptance=1.0,
                                    leaf_size=1,
                                )
out_fmm = zeros(3, npoints)
pnl.Uind!(body, points, out_fmm, backend_fmm)

# plot results
out_direct_z = out_direct[1, :]
out_fmm_z = out_fmm[1, :]
plt.figure("verify panel")
plt.cla()
plt.plot(dz, abs.(out_direct_z - out_fmm_z), label="error")
plt.yscale("log")
plt.xlabel("z (above triangle)")
plt.ylabel("Uind_z (direct - FMM)")
# plt.plot(dz, out_direct_z, label="direct")
# plt.plot(dz, out_fmm_z, label="FMM", "--")
# plt.legend()
plt.tight_layout()
plt.show()
