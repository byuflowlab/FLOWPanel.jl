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

# OpenVSP Model
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

