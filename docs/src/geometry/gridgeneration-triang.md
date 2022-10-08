# Grid Triangulation

Since some element types require planar panels, the quadrilateral grids generated through the previous methods need to be transformed into triangular grids that ensure planar panels. This is done by splitting each quadrilateral panel into two triangular panel through the following function:

```@docs
FLOWPanel.GeometricTools.GridTriangleSurface
```

The indexing pattern of a triangulated lofted surface is shown below

```@raw html
<center>
  <table>
      <tr>
          <th>
              <center>Node Index</center>
          </th>
          <th>
              <center>Cell Index</center>
          </th>
      </tr>
      <tr>
          <th>
              <img src="../../assets/images/loft-triang-nodeindex00.png" alt="Pic here" width="450px">
          </th>
          <th>
              <img src="../../assets/images/loft-triang-cellindex00.png" alt="Pic here" width="450px">
          </th>
      </tr>
  </table>
</center>

<br><br>

<center>
  <table>
      <tr>
          <th>
              <center>First Coordinate</center>
          </th>
          <th>
              <center>Second Coordinate</center>
          </th>
      </tr>
      <tr>
          <th>
              <img src="../../assets/images/loft-triang-cellcoordinate03.png" alt="Pic here" width="450px">
          </th>
          <th>
              <img src="../../assets/images/loft-triang-cellcoordinate02.png" alt="Pic here" width="450px">
          </th>
      </tr>
  </table>
</center>

<br><br>

<center>
  <table>
      <tr>
          <th>
              <center>Cell (1, 1)</center>
          </th>
          <th>
              <center>Cell (31, 1)</center>
          </th>
          <th>
              <center>Cell (30, 2)</center>
          </th>
      </tr>
      <tr>
          <th>
              <img src="../../assets/images/loft-triang-cellnodes00.0000.png" alt="Pic here" width="450px">
          </th>
          <th>
              <img src="../../assets/images/loft-triang-cellnodes00.0002.png" alt="Pic here" width="450px">
          </th>
          <th>
              <img src="../../assets/images/loft-triang-cellnodes00.0005.png" alt="Pic here" width="450px">
          </th>
      </tr>
      <tr>
          <th>
              <center>Cell (2, 1)</center>
          </th>
          <th>
              <center>Cell (32, 1)</center>
          </th>
          <th>
              <center>Cell (1, 2)</center>
          </th>
      </tr>
      <tr>
          <th>
              <img src="../../assets/images/loft-triang-cellnodes00.0001.png" alt="Pic here" width="450px">
          </th>
          <th>
              <img src="../../assets/images/loft-triang-cellnodes00.0003.png" alt="Pic here" width="450px">
          </th>
          <th>
              <img src="../../assets/images/loft-triang-cellnodes00.0004.png" alt="Pic here" width="450px">
          </th>
      </tr>
  </table>
</center>

<br><br>
```

It is possible to define a normal vector on a planar panel, which is then used to define the panel coordinate system as explained in [Panel Definition](@ref). The normal vector of each panel is shown here below in the left, while the orthonormal bases $\left( \hat{\mathbf{t}},\, \hat{\mathbf{o}},\, \hat{\mathbf{n}} \right)$ are shown in the right ($\hat{\mathbf{t}} = \mathrm{red}$, $\hat{\mathbf{o}} = \mathrm{yellow}$, $\hat{\mathbf{n}} = \mathrm{green}$).

```@raw html
<center>
  <table>
      <tr>
          <th>
              <center>Normals</center>
          </th>
          <th>
              <center>Panel Coordinate System</center>
          </th>
      </tr>
      <tr>
          <th>
              <img src="../../assets/images/loft-triang-normals00.png" alt="Pic here" width="450px">
          </th>
          <th>
              <img src="../../assets/images/loft-triang-vectors03.png" alt="Pic here" width="450px">
          </th>
      </tr>
  </table>
</center>
```
