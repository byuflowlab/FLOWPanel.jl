# ------------- CENTERBODY EXAMPLE ---------------------------------------------

output_name = "centerbody.md"
data_path = joinpath(module_path, "..", "resources", "data")

header = """
# Centerbody
"""


open(joinpath(output_path, output_name), "w") do fout

    println(fout, header)

    println(fout, """
    ```@raw html
    <center>
      <img src="../../assets/images/centerbody-viz00.png" alt="Pic here" style="width: 100%;"/>
    </center>
    ```

    In this example we solve the flow around a body of revolution resembling
    the centerbody (hub) of a ducted fan.
    """)

    # -------- Source Solver --------------------------
    println(fout, "## Source Elements")

    println(fout, """
    First we run this example with source elements, which are especially
    accurate for non-lifting bodies.

    """)

    println(fout, "```julia")

    open(joinpath(pnl.examples_path, "centerbody.jl"), "r") do fin
        for l in eachline(fin)
            if contains(l, "COMPARISON TO EXPERIMENTAL DATA")
                break
            end

            println(fout, l)
        end
    end

    println(fout, "```")

    println(fout, """
    (see the complete example under
    [examples/centerbody.jl](https://github.com/byuflowlab/FLOWPanel.jl/blob/master/examples/centerbody.jl)
    )
    """)


    # -------- Comparison to Experimental --------------------------

    println(fout, """
    ## Slice

    FLOWPanel provides the following function to obtain the solution field
    along a slice along a body:

    ```@docs
    FLOWPanel.slicefield
    ```

    Now we process the solution to plot the surface velocity along a slice
    of the body of revolution.

    """)

    println(fout, "```julia")

    open(joinpath(pnl.examples_path, "centerbody_postprocessing.jl"), "r") do fin
        for l in eachline(fin)
            if contains(l, "Beautify")
                break
            end

            println(fout, l)
        end
    end

    println(fout, "```")


    println(fout, """
    ```@raw html
    <center>
        <br><b>Surface velocity</b><br>
        <img src="../../assets/images/centerbody-lewis00-velocity-source.png" alt="Pic here" style="width: 60%;"/>
    </center>
    ```
    """)

    # -------- Vortex Ring Solver --------------------------

    println(fout, """

    ## Vortex Ring Elements

    While source elements are physically adequate to model a non-lifting body,
    in some circumstances it may be benefitial to use all vortex ring elements.
    A thick body with only vortex ring elements leads to a surface velocity
    that is inaccurate at the exact surface of the body, but that
    approximates the physical solution away from the surface. For this
    reason, we probe the velocity used to calculate Cp slightly away from
    the body.

    Here we repeat the example using only vortex ring elements.

    """)

    println(fout, "```julia")

    open(joinpath(pnl.examples_path, "centerbody_vortexring.jl"), "r") do fin
        for l in eachline(fin)
            if contains(l, "COMPARISON TO EXPERIMENTAL DATA")
                break
            end

            println(fout, l)
        end
    end

    println(fout, "```")


    println(fout, """
    ```@raw html
    <center>
        <br><b>Surface velocity</b><br>
        <img src="../../assets/images/centerbody-lewis00-velocity-vortexring.png" alt="Pic here" style="width: 60%;"/>
    </center>
    ```
    """)

end
