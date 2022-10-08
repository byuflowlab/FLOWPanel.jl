# Looped Grid

For some applications it is convenient to define a grid that loops on itself, making the last cell in a dimension close back with the first cell. The `Grid` object accepts an optional argument indicating the dimension that loops on itself as follows:

`Grid(P_min, P_max, NDIVS, loop_dim)`

where `loop_dim` is an integer indicating the dimension to loop.

Here is an example of a cylindrical grid $(r, \theta)$ where the nodes at $\theta=0^\circ$ and $\theta=360^\circ$ end up overlapping (notice nodes 1, 2, and 3 overlapping with 16, 17, and 18).

```julia
import GeometricTools as gt

R = 1.0                   # Radius of circle

P_min = [R/4, 0*pi/180]   # Lower boundaries r, theta
P_max = [R, 360*pi/180]   # Upper boundaries r, theta

NDIVS = [2, 5]          # 15 radius divisions, 60 angle divisions

# Generates the grid as cartesian
grid = gt.Grid(P_min, P_max, NDIVS)

# Converts to cylindrical
gt.transform!(grid, gt.cylindrical2D)

gt.plot(grid; labelnodes=true, labelcells=!true, labelndivs=!true, fontsize=8)
```
```@raw html
<img src="../../assets/images/geometry-loopedgrid00.png" alt="Pic here" style="width: 600px;"/>
<br><br><br><br>
```

Using `loop_dim=2`, now the $\theta$-coordinate loops on itself, merging the nodes that were overlapping:


```julia
loop_dim = 2              # Loop the theta dimension

# Generates the grid as cartesian
grid = gt.Grid(P_min, P_max, NDIVS, loop_dim)

# Converts to cylindrical
gt.transform!(grid, gt.cylindrical2D)

gt.plot(grid; labelnodes=true, labelcells=!true, labelndivs=!true, fontsize=8)
```
```@raw html
<img src="../../assets/images/geometry-loopedgrid01.png" alt="Pic here" style="width: 600px;"/>
```

## Example — Airfoil Contour

We will now both a space transformation and a looped grid to generate a two-dimensional airfoil section. The airfoil section is generated through the following procedure:

1. Read a collection of points from a CSV file describing an airfoil contour ("Original Airfoil Geometry").
2. Split the contour into upper and lower surfaces.
3. Spline through the points generating two analytic curve parameterized from inputs 0 to 1.
4. Generate a 1D line that represents a structured grid ("Original Grid").
5. Apply a space transformation to the 1D line to curve it into the shape of the airfoil contours ("Transformed Grid").

```julia
import FLOWPanel as pnl
import GeometricTools as gt
                                # Airfoil data path
airfoil_path = joinpath(pnl.def_data_path, "airfoils");

# ----------------- READ AND PARAMETERIZE AIRFOIL -----------------------------
# Read airfoil contour
x,y = gt.readcontour(joinpath(airfoil_path, "S809.txt"); header_len=2)

gt.plot_airfoil(x,y; style="--.k", title_str="Original Airfoil Geometry")
```
```@raw html
<img src="../../assets/images/geometry-loopedgrid-airfoil00.png" alt="Pic here" style="width: 900px;"/>
<br><br><br><br>
```

```julia
# Separate upper and lower surfaces to make the contour injective in x
upper, lower = gt.splitcontour(x,y)

# Parameterize both surfaces independently
fun_upper = gt.parameterize(upper[1], upper[2], zeros(size(upper[1])); inj_var=1, s=1e-6)
fun_lower = gt.parameterize(lower[1], lower[2], zeros(size(lower[1])); inj_var=1, s=1e-6)


# ----------------- CREATE GRID OF AIRFOIL CONTOUR ----------------------------
# Create the grid as a quasi-1D line
#   `X[1]` is the arc-length around the airfoil contour (between 0 and 1),
#   `X[2]` is a dummy value.

P_min = [0, 0]             # Lower boundaries of arclength and dummy
P_max = [1, 0]             # Upper boundaries of arclength and dummy
NDIVS = [10, 0]            # 100 arclength divisions, 0 dummys  divisions
loop_dim = 1               # Loop the arclength dimension

grid = gt.Grid(P_min, P_max, NDIVS, loop_dim)

gt.plot(grid; labelnodes=true, labelcells=true, labelndivs=!true,
                fontsize=8, fig_name="org", title_str="Original Grid")
```
```@raw html
<img src="../../assets/images/geometry-loopedgrid-airfoil01.png" alt="Pic here" style="width: 600px;"/>
<br><br><br><br>
```

```julia
# Create space transformation function
function my_space_transform(X)
    if X[1]<0.5
        return fun_upper(1-2*X[1])[1:2]
    else
        return fun_lower(2*(X[1]-0.5))[1:2]
    end
end

# Transforms the quasi-1D line into the 2D airfoil contour
gt.transform!(grid, my_space_transform)

gt.plot(grid; labelnodes=true, labelcells=true, labelndivs=!true,
                fontsize=8, title_str="Transformed grid");
```
```@raw html
<img src="../../assets/images/geometry-loopedgrid-airfoil02.png" alt="Pic here" style="width: 600px;"/>
```
* Green = Cell index
* Black = Node index
