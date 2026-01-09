# ------------- LIFTING LINE WEBER EXAMPLE -------------------------------------

output_name = "ll-a50k27"
data_path = joinpath(module_path, "..", "resources", "data")




# -------- AOA SWEEP -----------------------------------------------------------
open(joinpath(output_path, output_name*".md"), "w") do fout

    println(fout, """
    In this example we use the lifting line method to solve the aerodynamics of 
    a \$36^\\circ\$ swept-back wing with a cambered NACA 64(1)-612 airfoil from the 
    NACA RM A50K27 Flying Wing​ report.
    """)

    println(fout, "```julia")

    input_name = "liftingline_a50k27.jl"
    break_flag = "COMPARISON TO EXPERIMENTAL DATA"

    open(joinpath(pnl.examples_path, input_name), "r") do fin
        for (li, l) in enumerate(eachline(fin))
            if contains(l, break_flag)
                break
            end

            if (
                    contains(l, "plotformat.jl") || contains(l, "@L_str") ||
                    contains(l, "import CSV") || contains(l, "import DataFrames") ||
                    contains(l, "save_outputs    =") ||
                    contains(l, "n64_1_A612-Re0p5e6-smooth180-4.csv") ||
                    contains(l, "distributions = []") || 
                    contains(l, "output_distributions=distributions") ||
                    contains(l, "sweepname = run_name") ||
                    contains(l, "plots_path") || contains(l, "extraplots_path")
                )
                nothing
            else
                println(fout, l)
            end
        end
    end

    println(fout, "```")

    println(fout, """
    ```@raw html
    <span style="font-size: 0.9em; color:gray;"><i>
        Run time: ~7 seconds, evaluated 202 AOAs with 86% success rate. <br>
    </i></span>
    <br><br>
    ```
    (see the complete example under
    [examples/liftingline_a50k27.jl](https://github.com/byuflowlab/FLOWPanel.jl/blob/master/examples/liftingline_a50k27.jl)
    to see how to postprocess the spanwise loading that is plotted below)

    ```@raw html
    <center>
        <br><br><b>Wing Polar</b><br>
        <img src="../../assets/images/ll-a50k27-sweep-CLCDCm.png" alt="Pic here" style="width: 100%;"/>

        <br><br><b>Wing Polar Post-Stall</b><br>
        <img src="../../assets/images/ll-a50k27-sweep-CLCDCm-zoomout.png" alt="Pic here" style="width: 100%;"/>
    </center>
    ```

    """)

    println(fout, """

    !!! details "Tip"
        You can also automatically run this example and generate these plots
        with the following command:
        ```julia
        import FLOWPanel as pnl

        include(joinpath(pnl.examples_path, "liftingline_a50k27.jl"))
        ```
    """
    )

end
