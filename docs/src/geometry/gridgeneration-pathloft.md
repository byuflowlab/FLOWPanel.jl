# Path-Lofting Method

The following function in GeometricTools facilitates the path lofting of a
surface geometry through a set of cross sections and a path.
This method will use an Akima spline to extrapolate in between the sections and
path points that are given by the user, which resembles the flexibility of NURBs
but without weights to make the process simpler for the user.

```@docs
FLOWPanel.GeometricTools.surface_pathloft
```

## Example
First, we generate a set of conical sections for the loft:

```julia
import GeometricTools as gt
import FLOWMath as math
import PyPlot as plt

import LinearAlgebra: det
import PyPlot: @L_str

# Circular-to-rectangular transition
areadist        = [ # (xpos, area/R^2)      # Area distribution
                    (0.0, pi),
                    (1.0, pi)
                  ]
ar              = 2.0                       # Rectangle's aspect ratio
flatbottom      = false                     # If true, makes the lower inner surface flat, if false, the geometry is symmetric

sections_params = [ # (x-pos (0==LE, 1==TE), conic shape parameter rho, aspect ratio y/z)
                                (0.0,  0.4142, 1.0),
                                (0.125, 0.4142, 1.0),
                                (0.25, 0.4142, 1.0),
                                (0.8,  0.82,   ar),
                                (0.9,  0.82,   ar),
                                (1.0,   0.82,   ar),
                            ]

# npoints         = 20                       # Number of points per conic segment per section
npoints         = 100

# Spline area/R^2 distribution
area_spl = math.Akima([xa[1] for xa in areadist], [xa[2] for xa in areadist])

# Generate deforming sections from conic contours
sections = [ # (x, yzpoints)
            ]

# Plot area/R^2 distribution
fig = plt.figure(figsize=[7, 5]*5/9)
ax = fig.gca()
ax.plot(getindex.(areadist, 1), getindex.(areadist, 2)/pi, "ok")
ax.plot(range(0, 1, length=100), area_spl(range(0, 1, length=100))/pi, "-r")
ax.set_title("Area distribution spline")
ax.set_xlabel(L"x/d")
ax.set_ylabel(L"\mathrm{Area} / \pi d^2")

# Generate and plot conic sections
fig = plt.figure(figsize=[7, 5]*5/9)
ax = fig.gca()
ax.set_aspect("equal")

for (xpos, rho, ar) in sections_params

    # Contour points
    Ps = [ [1, 0], [0, 1], [-1, 0], [0, -1] ]

    # Control points
    CPs = [ [1, 1], [-1, 1], [-1, -1], [1, -1] ]

    # Apply aspect ratio
    Ps = [ [x, y/ar] for (x, y) in Ps ]
    CPs = [ [x, y/ar] for (x, y) in CPs ]

    # Conic shape parameter for each control point
    rhos = [ rho for i in 1:length(CPs) ]

    # Parametric probing points along each segment of the compound conic
    ss = [ collect(range(0, 1, length=npoints)) for i in 1:length(CPs) ]

    points = gt.conic_cross_section(Ps, CPs, rhos, ss)

    # Calculate original area of this compound conic
    area_org = 0
    for (i, X) in enumerate(points)

        Y = points[i%length(points) + 1]

        area_org += det(hcat(X, Y))/2

    end

    # Scale contour by target area distribution
    area_trg = area_spl(xpos)
    points *= sqrt(area_trg / area_org)

    # Verify new area
    area_new = 0
    for (i, X) in enumerate(points)

        Y = points[i%length(points) + 1]

        area_new += det(hcat(X, Y))/2

    end
    println("Section x=$(xpos)")
    println("\tOriginal area: $(area_org)")
    println("\tTarget area: $(area_trg)")
    println("\tActual area: $(area_new)")

    # Move contour down to have a straight bottom surface
    zmin = flatbottom ? minimum(z for (y, z) in points) : -1.0
    for X in points
        X[2] -= (zmin + 1.0)
    end

    # Plot contour
    ys = [y for (y, z) in points]
    zs = [z for (y, z) in points]
    ax.plot(ys, zs, ".-", label=L"x = "*"$(xpos)", alpha=0.8)

    # Convert points from array of arrays to a matrix
    points_matrix = hcat(ys, zs)

    push!(sections, (; x=xpos, points=points_matrix))
end
```

Now we define the path and use the sections to generate the path loft:


```julia
save_path = "pathloft-test000"

xpos_NDIVS = 7*10
sec_NDIVS = 20*10

sec_NDIVS = [[(1.0, sec_NDIVS, 1.0, true)] for i in 1:6]

sections = [(pos, points) for (pos, points) in sections]


path_Xs = [[1, 0, 0], [3.0, 0, 0], [4.0, 0.5, 0], [6.0, 1, 1]]
path_normals = [[1, 0, 0], [1, 0, 0], [1, 1, 0]/sqrt(2), [1, 1, 1]/sqrt(3)]
path_twists = [0, 45, 0, -45]

path = collect(zip(path_Xs, path_normals, path_twists))

grid = gt.surface_pathloft(sections, sec_NDIVS, xpos_NDIVS, path;
                                    save_path=save_path, paraview=true)
```
