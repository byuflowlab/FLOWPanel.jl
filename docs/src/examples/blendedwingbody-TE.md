# Export Trailing Edge

Along with the mesh, we also need to export the trailing edge as a line of points
that we can use to identify the trailing edge in the mesh.

First, close the Gmsh session from the previous section and launch Gmsh again to
clean up the workspace. Then follow these steps using the same STEP file that
you were processing before:

0. Import STEP file
    * `File > Merge`

1. Create physical group with only the TE curves
    * Make sure that geometric curves are visible: `Tools > Options > Geometry > Visibility > ✓ Curves`
    * Create physical group: `Modules > Geometry > Physical Groups > Add > Curve`
    * Select all the curves that define the trailing edge
    ![pic](http://edoalvar2.groups.et.byu.net/public/FLOWPanel/gmsh02.png)

2. Discretize the TE
    * Set size factor to something small: `Tools > Options > Mesh > General > Element size factor > 0.01`
    * Discretize 1D curves: `Modules > Mesh > 1D`

3. *(Optional)* Verify discretization
    * Make nodes visible: `Tools > Options > Mesh > Visibility > ✓ Nodes`
    * Zoom into the trailing edge and visually confirm that the line is finelly discretized
    ![pic](http://edoalvar2.groups.et.byu.net/public/FLOWPanel/gmsh03.png)

4. Export discretized curve as `.msh` in ASCII format v4: `File > Export` and select `Mesh - Gmsh MSH (*.msh)` in the dropdown menu

!!! info "`.msh` File"
    The resulting `.msh` file is available here:
    [LINK](https://github.com/byuflowlab/FLOWPanel.jl/raw/master/examples/data/zeroebwb-TE.msh)
    (`right click → save as...`).

