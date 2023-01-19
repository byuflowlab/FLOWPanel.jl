output_name = "duct"
data_path = joinpath(module_path, "..", "resources", "data")



# ---------------------- AOA SWEEP ---------------------------------------------
open(joinpath(output_path, output_name*"-aoasweep.md"), "w") do fout

    println(fout, """
    ```@raw html
    <center>
      <img src="../../assets/images/duct-hill-aoa15-0001-viz-glyph00.png" alt="Pic here" style="width: 100%;"/>
    </center>
    ```

    In this example we solve the flow around a fan duct (*aka* annular airfoil
    or engine cowl) at an angle of attack leading to asymmetric flow.
    The geometry is generated as a body of revolution that is watertight.
    """)
    println(fout, "# AOA Sweep")

    println(fout, """
    Here we generate the duct as a body of revolution and then we run a sweep
    of inflow angle of attack
    """)

    println(fout, "```julia")

    open(joinpath(pnl.examples_path, "duct.jl"), "r") do fin
        for l in eachline(fin)
            # if contains(l, "COMPARISON TO EXPERIMENTAL DATA")
            #     break
            # end

            println(fout, l)
        end
    end

    println(fout, "```")

    println(fout, """
    (see the complete example under
    [examples/duct.jl](https://github.com/byuflowlab/FLOWPanel.jl/blob/master/examples/duct.jl)
    )
    """)


    println(fout, """
    ```@raw html
    <center>
        <br><b>No AOA (symmetric flow)</b><br>
        <img src="../../assets/images/duct-hill00-Cp-AOA0.png" alt="Pic here" style="width:100%;"/>
    </center>

    <center>
        <br><b>5° angle of attack</b><br>
        <img src="../../assets/images/duct-hill00-Cp-AOA5.png" alt="Pic here" style="width:100%;"/>
    </center>

    <center>
        <br><b>15° angle of attack</b><br>
        <img src="../../assets/images/duct-hill00-Cp-AOA15.png" alt="Pic here" style="width:100%;"/>
    </center>
    ```

    In the plots above, the upper and lower sides correspond to the slice
    shown here below:

    ```@raw html
    <center>
      <img src="../../assets/images/duct-hill-aoa15-slice02.png" alt="Pic here" style="width: 100%;"/>
    </center>
    ```
    """)
end



# ---------------------- FLUID DOMAIN ------------------------------------------
open(joinpath(output_path, output_name*"-fluiddomain.md"), "w") do fout

    println(fout, """
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
    """)

    println(fout, "```julia")

    open(joinpath(pnl.examples_path, "duct_postprocessing.jl"), "r") do fin

        start = false

        for l in eachline(fin)
            start = start || contains(l, "function generate_fluiddomain")

            if start
                println(fout, l)
            end

        end
    end

    println(fout, "```")

    println(fout, """
    The fluid domain can be processed to visualize streamlines and contours of
    velocity/pressure/potential fields as shown below (the Paraview state that
    generated these images is available here: [LINK](https://edoalvar2.groups.et.byu.net/public/FLOWPanel/ducthill.pvsm)
    <- right click and "save as")
    """)


    println(fout, """
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
    """)
end





# ---------------------- LEAST-SQUARE SOLVER -----------------------------------
open(joinpath(output_path, output_name*"-leastsquare.md"), "w") do fout

    println(fout, """
    # Least-Square Solver

    It is well known that a purely doublet (or vortex ring) solver encounters
    difficulties when the geometry is closed (*i.e.*, watertight).
    The flow field around a body is given by the vortex filaments that make the
    edges of the panel, and the
    strength of each filament is simply the difference between adjacent panels.
    Thus, in the absence of an open edge (like in a watertight geometry), the
    strengths of the vortex-ring elements become irrelevant, and the problem
    is purely determined by the difference in strength between adjacent panels.
    This leads to an overdetermined problem where one of the original degrees of
    freedom (panels strengths) has become redundant.

    Let \$n\$ the number of panels. The problem is well defined for a open
    geometry formulating the solver as
    ```math
    \\begin{align*}
        G \\Gamma = -b
    ,\\end{align*}
    ```
    where \$b \\in \\mathbb{R}^{n}\$ is the normal of the freestream condition
    at each panel that needs to be canceled ("no-flow-through" boundary
    condition), \$G \\in \\mathbb{R}^{n\\times n}\$ contains the geometric
    information of the panels and \$\\Gamma \\in \\mathbb{R}^{n}\$ is the
    strength of each vortex-ring panel.
    However, for a watertight geometry, \$G\$ is no longer full rank and the
    problem becomes ill-conditioned.
    Due to numerical roundoff, the system of equations can still be inverted
    but the numerical solution ends up giving panel strengths (vortex-ring
    circulations) that are in the order of \$10^{16}\$ and large numerical
    noise.

    In order to circumvent this issue, we can transform the original
    problem into a least-square problem as follows.
    Since one of the panel strengths is redundant in a watertight geometry,
    we can simply pick an arbitrary panel and prescribe an arbitrary strength.
    Then, \$G\$ has become a \$n\\times n-1\$ matrix, \$\\Gamma\$ is a
    vector of length \$n-1\$, while \$b\$ is still a vector of length \$n\$.
    To formulate the least-square problem, we substract the velocity \$b_p\$
    induced by the "prescribed" panel  to the right-hand side,
    ```math
    \\begin{align*}
        G \\Gamma = -b - b_p
    ,\\end{align*}
    ```
    and we solve the problem as
    ```math
    \\begin{align*}
        \\Gamma = - \\left( G^{t} G \\right)^{-1} G^{t} \\left( b + b_p \\right)
    ,\\end{align*}
    ```
    where the superscript \$t\$ denotes the transposed matrix.

    Either solver (*i.e.*, the original vortex-ring solver or the vortex-ring
    least-square one) is automatically called whenever the function
    `FLOWPanel.solve(body, Uinfs, Das, Dbs)` is called.
    The solver to be used is identified based on the body type.
    A body type `bodytype = pnl.RigidWakeBody{pnl.VortexRing}` corresponds to
    the original vortex-ring solver, while the least-square solver is called
    by declaring `bodytype = pnl.RigidWakeBody{pnl.VortexRing, 2}`.

    Even though both solvers lead to roughly the same flow field solution, the
    numerical noise of the ill-condition problem is evident when visualizing the
    potential field:

    ```@raw html
    <center>
        <br><b>Vortex-Ring Solver</b><br>
        <img src="../../assets/images/duct-hill0110-viz-phi00.png" alt="Pic here" style="width:70%;"/>
    </center>
    <br>
    <center>
        <br><b>Vortex-Ring Least-Square Solver</b><br>
        <img src="../../assets/images/duct-hill0111-viz-phi00.png" alt="Pic here" style="width:70%;"/>
    </center>
    ```
    """)

end
