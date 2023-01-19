# Fluid Domain

The previous script uses the flag `fluiddomain` to indicate whether
to generate the fluid domain around the duct or not. When
`fluiddomain=true`, the fluid domain is generated through the function call
`generate_fluiddomain(body, AOA, Vinf, d, aspectratio, save_path)`, which
creates a grid around the duct where the velocity, pressure, and potential
field are probed.
This grid is then saved as an
[XDMF](https://www.xdmf.org/index.php/XDMF_Model_and_Format) file that can
be visualized in [Paraview](https://www.paraview.org/).

`generate_fluiddomain` is defined under
[examples/duct_postprocessing.jl](https://github.com/byuflowlab/FLOWPanel.jl/blob/master/examples/duct_postprocessing.jl)
as follows:

```julia
function generate_fluiddomain(body, AOA, Vinf, d, aspectratio, save_path;
                                halfdomain=false, # Whether to cover only one side of the duct
                                gridname="fluidomain",
                                num=nothing,
                                verbose=true,
                                v_lvl=0
                                )

    if verbose; println("\t"^(v_lvl)*"Generating fluid domain..."); end;

    # ---------------- GENERATE FLUID DOMAIN GRID ------------------------------
    # Bounds of grid
    Pmax = d*[aspectratio*1.5,    aspectratio*0.005,  0.5*1.35] # Upper bound
    Pmin = d*[-aspectratio*0.35, -Pmax[2], -Pmax[3]]            # Lower bound

    halfdomain ? Pmin[3] = 0 : nothing

    # Grid discretization
    dx         = 0.005*d*aspectratio                            # Cell size
    dy, dz     = dx, dx

    NDIVS = ceil.(Int, (Pmax .- Pmin) ./ [dx, dy, dz]) # Divisions in each dimension

    # Generate grid
    @time grid  = pnl.gt.Grid(Pmin, Pmax, NDIVS) # Grid

    if verbose; println("\t"^(v_lvl+1)*"Grid size:\t\t$(NDIVS)"); end;
    if verbose; println("\t"^(v_lvl+1)*"Number of nodes :\t$(grid.nnodes)"); end;

    # Translate and rotate grid to align with freestream
    O = zeros(3)
    Oaxis = pnl.gt.rotation_matrix2(0, AOA, 0)
    pnl.gt.lintransform!(grid, Oaxis, O)

    # Targets where to probe the velocity
    targets = grid.nodes
    ntargets = size(targets, 2)

    # Array where to store potential and velocity
    phis = zeros(ntargets)
    Us = repeat(Vinf, 1, ntargets)

    # Calculate potential and velocity fields
    @time pnl.phi!(body, targets, phis)
    @time pnl.Uind!(body, targets, Us)

    # Save fields
    pnl.gt.add_field(grid, "phi", "scalar", phis, "node")
    pnl.gt.add_field(grid, "U", "vector", collect(eachcol(Us)), "node")

    # Output fluid domain
    @time vtks = pnl.gt.save(grid, gridname; path=save_path, num=num)

    return vtks
end
```
The fluid domain can be processed to visualize streamlines and contours of
velocity/pressure/potential fields as shown below (the Paraview state that
generated these images is available here: [LINK](https://edoalvar2.groups.et.byu.net/public/FLOWPanel/ducthill.pvsm)
<- right click and "save as")

```@raw html
    <table width=100%>
    <tr>
        <td>
            <center>
                <img src="../../assets/images/duct-hill-aoa15-0001-viz-glyph00.png" alt="Pic here" style="width:100%">
            </center>
        </td>
    </tr>
    <tr>
        <td>
            <center>
                <img src="../../assets/images/duct-hill-aoa15-0001-viz-streamline00.png" alt="Pic here" style="width:100%">
            </center>
        </td>
    </tr>
    <tr>
        <td>
            <center>
                <img src="../../assets/images/duct-hill-aoa15-0001-viz-u00.png" alt="Pic here" style="width:100%">
            </center>
        </td>
    </tr>
    <tr>
        <td>
            <center>
                <img src="../../assets/images/duct-hill-aoa15-0001-viz-cp00.png" alt="Pic here" style="width:100%">
            </center>
        </td>
    </tr>
    <tr>
        <td>
            <center>
                <img src="../../assets/images/duct-hill-aoa15-0001-viz-phi00.png" alt="Pic here" style="width:100%">
            </center>
        </td>
    </tr>
</table>
```

