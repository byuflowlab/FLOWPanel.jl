```@raw html
<center>
  <img src="http://edoalvar2.groups.et.byu.net/public/FLOWPanel/zeroebwb009.jpg" alt="Pic here" style="width: 100%;"/>
</center>
```

In this example we analyze a blended wing body that is generated from CAD
(SolidWorks).
We show how to import the NURBS geometry into [Gmsh](https://gmsh.info) to generate an
unstructured watertight mesh with arbitrary regions of refinement, which is
then read by FLOWPanel using [Meshes.jl](https://juliageometry.github.io/MeshesDocs).

To run this tutorial you will have to install a few dependencies for
processing unstructured grids:
```julia
import Pkg

Pkg.add(Meshes)
Pkg.add(GeoIO)
Pkg.add(Rotations)
```

# CAD Model
1. Download original SolidWorks CAD file from GrabCAD:
    [grabcad.com/library/airbuszeroe-blended-wing-body-concept-1](https://grabcad.com/library/airbuszeroe-blended-wing-body-concept-1)
2. Remove unnecessary parting lines and ducted fan array
3. Added a parting line along the leading edge around which to customize refinement
4. Export as a STEP file

!!! info "STEP File"
    The resulting `.STEP` file is available here:
    [LINK](https://edoalvar2.groups.et.byu.net/public/FLOWPanel/zeroebwb.STEP)
    (`right click → save as...`).

