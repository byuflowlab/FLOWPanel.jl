# Export Trailing Edge

For each trailing edge in the mesh we need to generate a collection of
ordered points that we can use to recognize such edge in the `.msh` file.
We will use Gmsh to process the STEP file with NURBS and create such lines
of points along each trailing edge:

0. Import STEP file
    * `File > Merge`

1. Set discretization settings
    * Set size factor to something small: `Tools > Options > Mesh > General > Element size factor > 0.01`
    ![pic](http://edoalvar2.groups.et.byu.net/public/FLOWPanel/gmsh011.png)

Then for each trailing edge:

1. Create physical group with the TE curve
    * Make sure that geometric curves are visible: `Tools > Options > Geometry > Visibility > ✓ Curves`
    * Create physical group: `Modules > Geometry > Physical Groups > Add > Curve`
    * Select the curve that define the trailing edge
    ![pic](http://edoalvar2.groups.et.byu.net/public/FLOWPanel/gmsh012.png)

2. Discretize the TE:
    * `Modules > Mesh > 1D`

3. *(Optional)* Verify discretization
    * Make nodes visible: `Tools > Options > Mesh > Visibility > ✓ Nodes`
    * Zoom into the trailing edge and visually confirm that the line is finelly discretized
    ![pic](http://edoalvar2.groups.et.byu.net/public/FLOWPanel/gmsh013.png)
    * Hide nodes (otherwise you won't be able to select other curves in steps
        1 and 5): `Tools > Options > Mesh > Visibility > □ Nodes`

4. Export discretized curve as `.msh` in ASCII format v4: `File > Export`
    and select `Mesh - Gmsh MSH (*.msh)` in the dropdown menu

5. Delete physical group
    * `Modules > Geometry > Physical Groups > Remove > Curve`
    * Select the curve that we just exported

*Repeat steps 1 through 5 for each trailing edge.*

!!! info "Trailing Edge Files"
    The resulting TE `.msh` files are available here:
    * [`cessna-TE-leftwing.msh`](https://github.com/byuflowlab/FLOWPanel.jl/raw/master/examples/data/cessna-TE-leftwing.msh)
    * [`cessna-TE-rightwing.msh`](https://github.com/byuflowlab/FLOWPanel.jl/raw/master/examples/data/cessna-TE-rightwing.msh)
    * [`cessna-TE-leftelevator.msh`](https://github.com/byuflowlab/FLOWPanel.jl/raw/master/examples/data/cessna-TE-leftelevator.msh)
    * [`cessna-TE-rightelevator.msh`](https://github.com/byuflowlab/FLOWPanel.jl/raw/master/examples/data/cessna-TE-rightelevator.msh)
    * [`cessna-TE-rudder.msh`](https://github.com/byuflowlab/FLOWPanel.jl/raw/master/examples/data/cessna-TE-rudder.msh)
    (`right click → save as...`)

