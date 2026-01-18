# Space Transformations

A Cartesian grid can be used as the starting point to define a more complex structured grid through any arbitrary non-linear transformation:

```@docs
FLOWPanel.GeometricTools.transform!
FLOWPanel.GeometricTools.lintransform!
```
GeometricTools provides the following non-linear [orthogonal space transformations](https://en.wikipedia.org/wiki/Orthogonal_coordinates): `cylindrical3D(X)`, `cylindrical2D(X)`, `spherical3D(X)`, `parabolic3D(X)`, `paraboloidal3D(X)`, `elliptic3D(X; a=1)`, `prolate3D(X; a=1)`, `oblate3D(X; a=1)`, `bipolar3D(X; a=1)`, `toroidal3D(X; a=1)`, and `conical3D(X; b=2, c=1)`.

For linear transformations, `rotation_matrix(yaw::Real, pitch::Real, roll::Real)` returns the rotation matrix of such angles, and `axis_rotation(r::Array{Float64, 1}, angle_deg::Float64)` returns the transformation matrix of rotation around an arbitrary axis of unit vector `r`.

The user can also define an arbitrary space transformation. The airfoil example in the Looped Grid section shows the example of transforming a quasi-one-dimensional line into an two-dimensional airfoil contour.

## Example — 2D circular grid

To demonstrate the use of space transformations, in this example we generate a two-dimensional circular grid.
First we define the boundaries of the circular section $r_\text{min}, r_\text{max}, \theta_\text{min}, \theta_\text{max}$ as a cartesian grid of minimum and maximum bounds $P_\text{min}=(r_\text{min}, \theta_\text{min})$ and $P_\text{max}=(r_\text{max}, \theta_\text{max})$, respectively.
We then applying a cylindrical transformation on the grid.
The cylindrical transformation takes the current $x,y$ coordinates of the grid as $r,\theta$ for the new grid. In this case, $r_\text{min}=R/4$, $r_\text{max}=R$, $\theta_\text{min}=0^\circ$, $\theta_\text{max}=270^\circ$, as shown below.

```julia
import GeometricTools as gt

file_name = "circulargrid00"                # Output file
paraview = true                             # Whether to visualize the grid in Paraview

R = 1.0                   # Radius of circle

P_min = [R/4, 0*pi/180]   # Lower boundaries r, theta
P_max = [R, 270*pi/180]   # Upper boundaries r, theta

NDIVS = [15, 60]          # 15 radius divisions, 60 angle divisions

# Generate the Cartesian grid
grid = gt.Grid(P_min, P_max, NDIVS)

# Convert grid to cylindrical
gt.transform!(grid, gt.cylindrical2D)

# Visualization
if paraview

    # Output a vtk file
    gt.save(grid, file_name; format="vtk")

    # Call paraview
    run(`paraview --data=$file_name.vtk`)

else
    # Use PythonPlot instead of Paraview
    gt.plot(grid; labelnodes=!true, labelcells=!true, labelndivs=true)
end;
```
```@raw html
<img src="../../assets/images/circ01.png" alt="Pic here" style="width: 600px;"/>
<br><br>
```

Notice that `NDIVS` can still be used to defined sections of refinement:
```julia
import GeometricTools as gt

file_name = "circulargrid01"                # Output file
paraview = true                             # Whether to visualize the grid in Paraview

NDIVS = [   # r sections
            [(1.0, 15, 5.0, true)],
            # theta sections
            [(1/3, 15, 1/3, false),
             (1/3, 15, 3.0, true),
             (1/3, 15, 3.0, false)]
        ]

# Generates the grid as cartesian
grid = gt.Grid(P_min, P_max, NDIVS)

# Converts to cylindrical
gt.transform!(grid, gt.cylindrical2D)

# Visualization
if paraview

    # Outputs a vtk file
    gt.save(grid, file_name; format="vtk")

    # Calls paraview
    run(`paraview --data=$file_name.vtk`)
else
    # Use PythonPlot instead of Paraview
    gt.plot(grid; labelnodes=!true, labelcells=!true, labelndivs=true)
end;
```
```@raw html
<img src="../../assets/images/circ02.png" alt="Pic here" style="width: 600px;"/>
```
