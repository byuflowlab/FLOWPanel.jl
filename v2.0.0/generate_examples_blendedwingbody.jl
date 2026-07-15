output_name = "blendedwingbody"
data_path = joinpath(module_path, "..", "resources", "data")



# ---------------------- CAD Model ---------------------------------------------
open(joinpath(output_path, output_name*"-cad.md"), "w") do fout

    println(fout, """
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
    """)

    println(fout, "# CAD Model")

    println(fout, """
    1. Download the original SolidWorks CAD file from GrabCAD:
        [grabcad.com/library/airbuszeroe-blended-wing-body-concept-1](https://grabcad.com/library/airbuszeroe-blended-wing-body-concept-1)
        ![pic](http://edoalvar2.groups.et.byu.net/public/FLOWPanel/zeroebwb-sw000.png)
    2. Remove unnecessary parting lines and ducted fan array
        ![pic](http://edoalvar2.groups.et.byu.net/public/FLOWPanel/zeroebwb-sw002.png)
    3. Add a parting line along the leading edge around which we will later
        customize the mesh refinement
        ![pic](http://edoalvar2.groups.et.byu.net/public/FLOWPanel/zeroebwb-sw001.png)
    4. Export as a STEP file

    !!! info "STEP File"
        The resulting `.STEP` file is available here:
        [LINK](https://github.com/byuflowlab/FLOWPanel.jl/raw/master/examples/data/zeroebwb.STEP)
        (`right click → save as...`).
    """)
end


# ---------------------- Unstructured Meshing ----------------------------------
open(joinpath(output_path, output_name*"-gmsh.md"), "w") do fout

    println(fout, "# Unstructured Mesh Generation")

    println(fout, """

    We now open the NURBS from the STEP file in
    [Gmsh](https://gmsh.info/#Download) and generate an
    unstructured surface mesh with refinement along the leading edge of
    the wing.

    > **NOTE:** *the `CurveList`s described below are valid for the file
    > [`zeroebwb.STEP`](https://github.com/byuflowlab/FLOWPanel.jl/raw/master/examples/data/zeroebwb.STEP)
    > and Gmsh v4.8.4*

    0. Import STEP file
        * `File > Merge`

    1. Create physical group with all NURBS surfaces
        * Make sure that geometric surfaces are visible: `Tools > Options > Geometry > Visibility > Surfaces`
        * Click the `Y` button in the bottom left corner to get a top view
        * Create physical group: `Modules > Geometry > Physical Groups > Add > Surface`
        * Press `CTRL + Right Click` to box-select **only half of the airframe** (we will export half the mesh
            and will later mirror it to obtain a mesh that is exactly symmetric)
        ![pic](http://edoalvar2.groups.et.byu.net/public/FLOWPanel/gmsh020.png)

    2. Create criterion for refinement at LE
        * `Modules > Mesh > Define > Size fields`
        * Define leading edge lines: `Distance` with `CurvesList` set to `33, 25, 53, 45` and `1000` for `NumPointsPercurve`
        * Define transition lines: `Distance` with `CurvesList` set to `18, 19, 56, 54` and `1000` for `NumPointsPercurve`
        * Define upper-bound lines of leading edge: `Distance` with `CurvesList` set to `30, 22, 55, 47` and `1000` for `NumPointsPercurve`
        * Define lower-bound lines of leading edge: `Distance` with `CurvesList` set to `35, 28, 51, 43` and `1000` for `NumPointsPercurve`
        * Define cell size based on distance to leading edge: `MathEval` with the formula `1.5*(F1^1.2 + 500/F2^0.5)/2 + 30`
        * Define cell size based on thickness of leading edge: `MathEval` with the formula `3*((5*F1+F3+F4)/7)^1.5 + 5`
        * Mix both cell sizes: `Min` with `FieldsList` set to `5, 6` select "Set as background field" and click Apply. Defined this way, it refines the mesh in the proximity of the LE lines based on LE thickness, while transitioning the refinement close to the transition lines.

    3. Set mesher settings
        * Set size factor: `Tools > Options > Mesh > General > Element size factor` and set to `0.1`
        * Select `MeshAdapt` for the 2D algorithm

    4. Mesh the NURBS surface
        * Discretize 1D curves: `Modules > Mesh > 1D`
        * Discretize 2D surfaces: `Modules > Mesh > 2D`
        * Smooth out the discretization by clicking `Modules > Mesh > Smooth 2D` a few times
        ![pic](http://edoalvar2.groups.et.byu.net/public/FLOWPanel/gmsh021.png)

    5. Export mesh as `.msh` in ASCII format v4: `File > Export` and select `Mesh - Gmsh MSH (*.msh)` in the dropdown menu

    !!! info "`.msh` File"
        The resulting `.msh` file is available here:
        [LINK](https://github.com/byuflowlab/FLOWPanel.jl/raw/master/examples/data/zeroebwb.msh)
        (`right click → save as...`).
    """)
end


# ---------------------- Export Trailing Edge ----------------------------------
open(joinpath(output_path, output_name*"-TE.md"), "w") do fout

    println(fout, "# Export Trailing Edge")

    println(fout, """

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
    """)
end


# ---------------------- Import and Solve --------------------------------------
open(joinpath(output_path, output_name*"-aero.md"), "w") do fout


    println(fout, """
    ```@raw html
    <center>
      <img src="http://edoalvar2.groups.et.byu.net/public/FLOWPanel/zeroebwb005.jpg" alt="Pic here" style="width: 100%;"/>
    </center>
    ```
    """)

    println(fout, "# Import Mesh and Solve")

    println(fout, """
    Here we import the mesh into FLOWPanel using
    [Meshes.jl](https://juliageometry.github.io/MeshesDocs), identify the
    trailing edge, and run the watertight solver.


    !!! info "Default Gmsh files"
        We have pre-generated and uploaded a Gmsh mesh to this example so that
        you can run this section without needing to complete the previous
        sections. However, if you would like to use your own mesh, simply change
        `read_path`, `meshfile`, and `trailingedgefile` to point to
        your files.
    """)

    println(fout, "```julia")

    open(joinpath(pnl.examples_path, "blendedwingbody.jl"), "r") do fin
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
        Number of panels: 22,000. <br>
        Run time: ~90 seconds on a Dell Precision 7760 laptop, no GPU (~10 seconds with GPU). <br>
    </i></span>
    <br><br>
    ```


    ```@raw html
    <center>
      <img src="http://edoalvar2.groups.et.byu.net/public/FLOWPanel/zeroebwb009.jpg" alt="Pic here" style="width: 100%;"/>
      <img src="http://edoalvar2.groups.et.byu.net/public/FLOWPanel/zeroebwb004.jpg" alt="Pic here" style="width: 100%;"/>
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

        include(joinpath(pnl.examples_path, "blendedwingbody.jl"))
        ```
    """
    )

    println(fout, """

    !!! tip "Checklist for importing meshes"
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


    """)

end


# ---------------------- GPU and CPU Acceleration ------------------------------
open(joinpath(output_path, output_name*"-gpucpu.md"), "w") do fout

    println(fout, "# GPU and CPU Acceleration")

    println(fout, """

    ## CPU Multi-Threading

    The kernels and solvers implemented in FLOWPanel are parallelized (threaded)
    in CPU by default.
    However, in order to activate the CPU parallelization, the user needs to
    [launch Julia with multi-threading activated](https://docs.julialang.org/en/v1/manual/multi-threading/#Starting-Julia-with-multiple-threads).
    For instance, to launch Julia with 4 threads:
    ```bash
    \$ julia --threads 4
    ```

    You can then verify that the 4 threads became available:

    ```julia-repl
    julia> Threads.nthreads()
    4
    ```

    ## Porting to GPU

    The solver can be seamlessly ported to GPU by indicating the type of
    array to be used internally.
    [The Julia GPU interface](https://juliagpu.org/) is the same for any GPU hardware and platform
    (NVIDIA CUDA, AMD ROCm, and Mac Metal), however, we have only tested NVIDIA
    GPUs.

    For an NVIDIA GPU, first import the CUDA package before running the code of
    the previous section,
    ```julia-repl
    julia> import CUDA
    ```
    check that the GPU hardware is ready to be used,
    ```julia-repl
    julia> CUDA.functional()
    true
    ```
    and instead of letting the solver default its internal arrays to CPU, change
    the solver call from
    ```julia-repl
    julia> pnl.solve(body, Uinfs, Das, Dbs)
    ```
    to
    ```julia-repl
    julia> pnl.solve(body, Uinfs, Das, Dbs; GPUArray=CUDA.CuArray{Float32})
    ```

    For AMD GPU:
    ```julia-repl
    julia> import AMDGPU
    julia> AMDGPU.functional()
    true
    julia> AMDGPU.functional(:MIOpen)
    true
    julia> pnl.solve(body, Uinfs, Das, Dbs; GPUArray=AMDGPU.ROCArray{Float32})
    ```

    For Metal GPU:
    ```julia-repl
    julia> import Metal
    julia> Metal.functional()
    true
    julia> pnl.solve(body, Uinfs, Das, Dbs; GPUArray=Metal.MtlArray{Float32})
    ```


    !!! info "GPU"
        We have only tested NVIDIA GPUs
    """)

end
