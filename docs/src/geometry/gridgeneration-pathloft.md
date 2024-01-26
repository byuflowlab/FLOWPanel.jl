# Path-Lofting Method

The following function in GeometricTools lofts a set of contours along a path.
This method uses an Akima spline to extrapolate in between the contours
and path points supplied by the user, which resembles the flexibility of NURBs
but without weights in order to make the process simpler for the user.

```@docs
FLOWPanel.GeometricTools.surface_pathloft
```

## Example
In this example we read two conic sections from a CSV and define a loft path
to generate a lofted surface grid. This animation shows the path of the loft (
points in black, normals in green, path spline in red):

```@raw html
<center>
<img src="http://edoalvar2.groups.et.byu.net/public/FLOWPanel/pathloft-example_4-small.gif" alt="Vid here", width="70%"/>
</center>
```

```julia
import GeometricTools as gt

import FLOWPanel: examples_path
import CSV
import DataFrames: DataFrame


save_path       = "pathloft-example"                # Where to save grid VTK
contour_path    = joinpath(examples_path, "data")   # Where to read the conic contour from


# ------------ Define contours of loft -----------------------------------------

# Read conic contours for loft
points1 = CSV.read(joinpath(contour_path, "conicsection1.csv"), DataFrame)
points1 = Matrix(points1)

points2 = CSV.read(joinpath(contour_path, "conicsection2.csv"), DataFrame)
points2 = Matrix(points2)

# Define loft sections
sections = [ # (non-dimensional arclength position along path, contour)
                (0.000,  points1),
                (0.125, points1),
                (0.250, points1),
                (0.800,  points2),
                (0.900,  points2),
                (1.000,  points2)
            ]


# ------------ Define loft path ------------------------------------------------

# Points along path
path_Xs = [
            [1.0, 0, 0],
            [2.0, 0, 0],
            [3.0, 0, 0],
            [4.0, 0.5, 0],
            [6.0, 1, 1]
          ]

# Section normals along path
path_normals = [
                [1, 0, 0],
                [1, 0, 0],
                [1, 0, 0],
                [1, 1, 0]/sqrt(2),
                [1, 1, 1]/sqrt(3)
               ]

# Twist of section contours along path
path_twists = [0, -30, -45, 0, -45]

# Collect the path
path = collect(zip(path_Xs, path_normals, path_twists))


# ------------ Generate lofted surface grid ----------------------------------

# Define discretization
section_NDIVS   = 200                               # Discretization of section contours
path_NDIVS      = 100                               # Discretization of path

# Generate loft
grid = gt.surface_pathloft(sections, path,
                            section_NDIVS, path_NDIVS; save_path=save_path, paraview=true)
```

`surface_pathloft(...; verify_spline=true, ...)` automatically generates the
following plots for the user to verify that the spline is correctly fitting
the raw contours and path.
Here are the contours:

```@raw html

<center>
  <table>
      <tr>
          <td>
              <img src="../../assets/images/pathloft-example-contour000.png" alt="Pic here" width="100%">
          </td>
          <td>
              <img src="../../assets/images/pathloft-example-contour001.png" alt="Pic here" width="100%">
          </td>
          <td>
              <img src="../../assets/images/pathloft-example-contour002.png" alt="Pic here" width="100%">
          </td>
      </tr>
      <tr>
          <td>
              <img src="../../assets/images/pathloft-example-contour003.png" alt="Pic here" width="100%">
          </td>
          <td>
              <img src="../../assets/images/pathloft-example-contour004.png" alt="Pic here" width="100%">
          </td>
          <td>
              <img src="../../assets/images/pathloft-example-contour005.png" alt="Pic here" width="100%">
          </td>
      </tr>
  </table>
</center>

<br><br>
```

The path is visualized showing the $x-y$ and $x-z$ coordinates of the path
points and normals, as shown below:

```@raw html
<center>
  <table>
      <tr>
          <td>
              <img src="../../assets/images/pathloft-example-path000.png" alt="Pic here" width="100%">
          </td>
      </tr>
  </table>
</center>
```

Also the contour twist is shown along the dimensional length of the path:

```@raw html
<center>
  <table>
      <tr>
          <td>
              <img src="../../assets/images/pathloft-example-path001.png" alt="Pic here" width="100%">
          </td>
      </tr>
  </table>
</center>
```

`surface_pathloft(...; paraview=true, ...)` calls ParaView after generating the
grid to visualize the loft:

```@raw html
<center>
  <img src="http://edoalvar2.groups.et.byu.net/public/FLOWPanel/pathloft-example_5-small.gif" alt="Vid here", width="70%"/>
  <img src="http://edoalvar2.groups.et.byu.net/public/FLOWPanel/pathloft-example_6-small.gif" alt="Vid here", width="70%"/>
</center>
```
