output_name = "cessna"
data_path = joinpath(module_path, "..", "resources", "data")



# ---------------------- CAD Model ---------------------------------------------
open(joinpath(output_path, output_name*"-openvsp.md"), "w") do fout

    println(fout, """
    ```@raw html
    <center>
      <img src="https://edoalvar2.groups.et.byu.net/public/FLOWPanel/cessna002.jpg" alt="Pic here" style="width: 100%;"/>
    </center>
    ```

    In this example we analyze a Cessna 210 aircraft that is generated using
    [OpenVSP](https://openvsp.org).
    We show how to polish and mesh the OpenVSP model, and the mesh is
    then read by FLOWPanel using [Meshes.jl](https://juliageometry.github.io/MeshesDocs).

    To run this tutorial you will have to install a few dependencies for
    processing unstructured grids:
    ```julia
    import Pkg

    Pkg.add(Meshes)
    Pkg.add(GeoIO)
    Pkg.add(Rotations)
    ```
    """)

    println(fout, "# OpenVSP Model")

    println(fout, """
    We will use the Cessna 210 model that is available in the [OpenVSP Ground
    School](https://vspu.larc.nasa.gov/example-file-download) as the starting
    point (you can also download the `.vsp3` file directly from
    [here](https://github.com/byuflowlab/FLOWPanel.jl/raw/master/examples/data/Cessna-210.vsp3),
    `right click → save as...`).
    Follow the steps below to polish the model, mesh it, and export the mesh.

    ## Polish the Model
    If you try to run the model as-is through the OpenVSP aero solver, you might
    notice that VSPAERO crashes.
    This is an indication that the model needs some clean up and polish before
    we can do any numerical analysis with it.

    We recommend doing the following modifications to the model,

    1. Round up the wing tips to avoid sharp edges
    2. Sharpen the wing trailing edge
    3. Merge the prop spinner with the main body by translating it to close the gap
    4. Adequately adapt the mesh resolution:
        * Refinement towards wing tips
        * Refinement towards both wing LE and TE
        * Avoid abrupt changes in cell sizes

    ```@raw html
    <center>
        <span style="font-size: 1.1em; color:black;"><b>
            Before
        </b></span>
        <br>
        <img src="https://edoalvar2.groups.et.byu.net/public/FLOWPanel/openvsp000.png" alt="Pic here" style="width: 100%;"/>
    </center>
    <br>
    <center>
        <span style="font-size: 1.1em; color:black;"><b>
            After
        </b></span>
        <br>
        <img src="https://edoalvar2.groups.et.byu.net/public/FLOWPanel/openvsp001.png" alt="Pic here" style="width: 100%;"/>
    </center>
    ```

    !!! info "OpenVSP Files"
        The `.vsp3` files are available here:
        [[before]](https://github.com/byuflowlab/FLOWPanel.jl/raw/master/examples/data/Cessna-210.vsp3)
        [[after]](https://github.com/byuflowlab/FLOWPanel.jl/raw/master/examples/data/cessna.vsp3)
        (`right click → save as...`).

    ## Quality Check
    Make sure that the OPENVSP aero solver (VSPAERO) can successfully perform
    analysis with the mesh:
    * Confirm that the analysis completes without returning `nan`s
    * Confirm that the analysis returns reasonable values for ``C_L`` and ``C_D``
    * Visualize the wake and make sure that the body is not improperly
        interfering with the wake (e.g.,, wake penetrating into the body)
    * Visualize ``C_p`` and make sure that there aren't any regions where the
        solution blows up

    If VSPAERO can't successfully do numerical analysis with this mesh, it is
    likely that FLOWPanel will also fail.

    ## Export

    * Export mesh: `File > Export > Gmsh (.msh)`
    * Export NURBS: `File > Export > Untrimmed STEP (.stp)`, set length unit to `MM`, and set tolerance to `1e-12`=

    OpenVSP exports `.msh` files in an older format that is not compatible with
    [Meshes.jl](https://juliageometry.github.io/MeshesDocs), so here we
    recommend using [Gmsh](https://gmsh.info) to convert the `.msh` file to a
    newer version:

    * Open the `.msh` file in Gmsh
    * Re-export `.msh` in ASCII format v4: `File > Export` and select
        `Mesh - Gmsh MSH (*.msh)` in the dropdown menu

    !!! info "Mesh Files"
        The resulting `.msh` and `.stp` files are available here:
        [`cessna.msh`](https://edoalvar2.groups.et.byu.net/public/FLOWPanel/cessna.msh)
        [`cessna.stp`](https://edoalvar2.groups.et.byu.net/public/FLOWPanel/cessna.stp)
        (`right click → save as...`).
    """)
end


# ---------------------- Export Trailing Edge ----------------------------------
open(joinpath(output_path, output_name*"-TE.md"), "w") do fout

    println(fout, "# Export Trailing Edge")

    println(fout, """

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
    """)
end


# ---------------------- Import and Solve --------------------------------------
open(joinpath(output_path, output_name*"-aero.md"), "w") do fout


    println(fout, """
    ```@raw html
    <center>
      <img src="https://edoalvar2.groups.et.byu.net/public/FLOWPanel/cessna000.jpg" alt="Pic here" style="width: 100%;"/>
    </center>
    ```
    """)

    println(fout, "# Import Mesh and Solve")

    println(fout, """
    Here we import the mesh into FLOWPanel using
    [Meshes.jl](https://juliageometry.github.io/MeshesDocs), identify the
    trailing edge, and run the watertight solver.


    !!! info "Default OpenVSP files"
        We have pre-generated and uploaded an OpenVSP mesh to this example so that
        you can run this section without needing to complete the previous
        sections. However, if you would like to use your own mesh, simply change
        `read_path`, `meshfile`, and `trailingedges` to point to
        your files.
    """)

    println(fout, "```julia")

    open(joinpath(pnl.examples_path, "cessna2.jl"), "r") do fin
        for l in eachline(fin)
            # if contains(l, "COMPARISON TO EXPERIMENTAL DATA")
            #     break
            # end

            println(fout, l)
        end
    end

    println(fout, "```")


    println(fout, """
    ```@raw html
    <span style="font-size: 0.9em; color:gray;"><i>
        Number of panels: 24,000. <br>
        Run time: ~100 seconds on a Dell Precision 7760 laptop (no GPU). <br>
    </i></span>
    <br><br>
    ```

    |                         | VSPAERO | FLOWPanel |
    | ----------------------: | :-----: | :-------: |
    | Lift ``C_L``            | 0.499   | 0.554     |
    | Drag ``C_D``            | 0.0397  | 0.0496    |


    ```@raw html
    <center>
      <img src="https://edoalvar2.groups.et.byu.net/public/FLOWPanel/cessna002.jpg" alt="Pic here" style="width: 100%;"/>
      <img src="https://edoalvar2.groups.et.byu.net/public/FLOWPanel/cessna005.jpg" alt="Pic here" style="width: 100%;"/>
    </center>
    <br><br>
    <br><br>
    ```

    """)


    println(fout, """

    !!! details "Tip"
        You can also automatically run this example
        with the following command:
        ```julia
        import FLOWPanel as pnl

        include(joinpath(pnl.examples_path, "cessna.jl"))
        ```
    """
    )

    println(fout, """

    !!! details "Checklist for importing meshes"
        1. **Check whether normals point into the body:** Using the flag
            `debug=true` in `pnl.save(body, run_name; path=save_path, debug=true)`
            will output the control points of the body along with the associated
            normal vector of each panel.
                We recommend opening the body and control points in ParaView and
            visualizing the normals with the Glyph filter.
                Whenever the normals are pointing into the body, the user needs
            to flip the offset of the
            control points with `CPoffset=-1e-14` or any other negligibly small
            negative number. This won't flip the normals outwards, but it will flip
            the zero-potential domain from outwards back to inside the body
            (achieved by shifting the control points slightly into the body).
            If you pull up the solution in ParaView and realize that the surface
            velocity is much smaller than the freestream everywhere along the
            aircraft, that's an indication that the normals are point inwards
            and you need to set `CPoffset` to be negative.
        2. **Check that the trailing edge was correctly identified:**
            `pnl.save(body, run_name; path=save_path)` automatically outputes the
            wake.
                We recommend opening the body and wake in ParaView and visually
            inspecting that the wake runs along the trailing edge line that you
            defined under `trailingedge`.
                If not successful, increase the resolution of `trailingedge` and tighten
            the tolerance to something small like
            `pnl.calc_shedding(grid, trailingedge; tolerance=0.0001*span)`.
        3. **Choose the right solver for the geometry:**
            Use the least-squares solver with watertight bodies
            (`bodytype = pnl.RigidWakeBody{pnl.VortexRing, 2}`), and the direct
            linear solver with open bodies
            (`bodytype = pnl.RigidWakeBody{pnl.VortexRing}`). The least-squares
            solver runs much faster in GPU
            (`pnl.solve(body, Uinfs, Das, Dbs; GPUArray=CUDA.CuArray{Float32})`),
            but it comes at the price of sacrificing accuracy (single precision
            numbers as opposed to double).

    !!! tip "Visualization"
        To help you practice in ParaView, we have uploaded the solution files
        of this simulation along with the ParaView state file (`.pvsm`) that
        we used to generate the visualizations shown above:
        [DOWNLOAD](https://edoalvar2.groups.et.byu.net/public/FLOWPanel/cessna.zip)

        To open in ParaView: `File → Load State → cessna.pvsm` then
        select "Search files under specified directory" and point it to the
        folder with the outputs of FLOWPanel.
    """)

end



# ---------------------- GPU and CPU Acceleration ------------------------------
open(joinpath(output_path, output_name*"-vspgeom.md"), "w") do fout

    println(fout, "# VSPGeom.jl")

    println(fout, """

    [VSPGeom.jl](https://flow.byu.edu/VSPGeom.jl) is a Julia package for reading
    and processing the geometry of OpenVSP, which we can also use with
    FLOWPanel.
    For instance, in OpenVSP we can generate a CFD-quality surface mesh that
    gets exported as an STL file.
    The STL can then be read following [these instructions](https://flow.byu.edu/VSPGeom.jl/dev/howto/#STL-Files) and wrapped as
    a [Meshes.jl](https://juliageometry.github.io/MeshesDocs) object that is
    passed to FLOWPanel.

    ```@raw html
    <center>
      <img src="https://edoalvar2.groups.et.byu.net/public/FLOWPanel/openvsp004.png" alt="Pic here" style="width: 100%;"/>
    </center>
    ```
    """)

end
