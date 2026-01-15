# ------------- LIFTING LINE WEBER EXAMPLE -------------------------------------

output_name = "ll-stabilityderivatives"
data_path = joinpath(module_path, "..", "resources", "data")




# -------- 4.2deg AOA ----------------------------------------------------------
open(joinpath(output_path, output_name*".md"), "w") do fout

    println(fout, "# Stability Derivatives")

    println(fout, """
    In this example we show how to calculate stability derivatives using automatic 
    differentiation. This example uses Weber's wing with 5deg dihedral, and
    the nonlinear lifting line method.
    """)

    println(fout, "```julia")

    input_name = "liftingline_stabilityderivatives.jl"
    break_flag = "COMPARISON TO EXPERIMENTAL DATA"

    open(joinpath(pnl.examples_path, input_name), "r") do fin
        for (li, l) in enumerate(eachline(fin))
            if contains(l, break_flag)
                break
            end

            if (
                    contains(l, "plotformat.jl") || contains(l, "@L_str") ||
                    contains(l, "import CSV") || contains(l, "import DataFrames") ||
                    contains(l, "save_outputs    =")
                )
                nothing
            else
                println(fout, l)
            end
        end
    end

    println(fout, "```")

    println(fout, """
    (see the complete example under
    [examples/liftingline_stabilityderivatives.jl](https://github.com/byuflowlab/FLOWPanel.jl/blob/master/examples/liftingline_stabilityderivatives.jl)
    to see how to postprocess the spanwise loading that is plotted below)

    """)

    println(fout, """

    !!! details "Tip"
        You can also automatically run this example and generate these plots
        with the following command:
        ```julia
        import FLOWPanel as pnl

        include(joinpath(pnl.examples_path, "liftingline_stabilityderivatives.jl"))
        ```
    """
    )

end

