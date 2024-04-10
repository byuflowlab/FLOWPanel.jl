# Unstructured Mesh Generation

We now open the NURBS from the STEP file in
[Gmsh](https://gmsh.info/#Download) and generate an
unstructured surface mesh with refinement along the leading edge of
the wing.


> **NOTE:** *the `CurveList`s described below are valid for the file
> [`zeroebwb.STEP`](https://edoalvar2.groups.et.byu.net/public/FLOWPanel/zeroebwb.STEP)
> and Gmsh v4.8.4*

0. Import STEP file
    * `File > Merge`

1. Create physical group with all NURBS surfaces
    * Make sure that geometric surfaces are visible: `Tools > Options > Geometry > Visibility > Surfaces`
    * Create physical group: `Modules > Geometry > Physical Groups > Add > Surface`
    * Press `CTRL + Right Click` to box-select
    ![pic](http://edoalvar2.groups.et.byu.net/public/FLOWPanel/gmsh000.png)

2. Create criterion for refinement at LE
    * `Modules > Mesh > Define > Size fields`
    * Define leading edge lines: `Distance` with `CurvesList` set to 33,25,53,45 and 1000 for `NumPointsPercurve`
    * Define transition lines: `Distance` with `CurvesList` set to 18,19,56,54 and 1000 for `NumPointsPercurve`
    * Define upper-bound lines of leading edge: `Distance` with `CurvesList` set to 30, 22, 55, 47 and 1000 for `NumPointsPercurve`
    * Define lower-bound lines of leading edge: `Distance` with `CurvesList` set to 35, 28, 51, 43 and 1000 for `NumPointsPercurve`
    * Define cell size based on distance to leading edge: `MathEval` with the formula `1.5*(F1^1.2 + 500/F2^0.5)/2 + 30`
    * Define cell size based on thickness of leading edge: `MathEval` with the formula `3*((5*F1+F3+F4)/7)^1.5 + 5`
    * Mix both cell sizes: `Min` with `FieldsList` set to 5,6, select "Set as background field" and click Apply. Defined this way, it refines the mesh in the proximity of the LE lines based on LE thickness, while transitioning the refinement close to the transition lines.

3. Set mesher settings
    * Set size factor: `Tools > Options > Mesh > General > Element size factor` and set to `0.1`
    * Select `MeshAdapt` for the 2D algorithm

4. Mesh the NURBS surface
    * Discretize 1D curves: `Modules > Mesh > 1D`
    * Discretize 2D surfaces: `Modules > Mesh > 2D`
    * Smooth out the discretization by clicking `Modules > Mesh > Smooth 2D` a few times
    ![pic](http://edoalvar2.groups.et.byu.net/public/FLOWPanel/gmsh001.png)

5. Export mesh as `.msh` in ASCII format v4: `File > Export` and select `Mesh - Gmsh MSH (*.msh)` in the dropdown menu



