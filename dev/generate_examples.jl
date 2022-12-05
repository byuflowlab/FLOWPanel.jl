import FLOWPanel as pnl

module_path = splitdir(@__FILE__)[1]      # Path to this module
output_path = joinpath(module_path, "examples") # Path where to store  examples as markdown

# ------------- SWEPT WING EXAMPLE ---------------------------------------------

input_name = "sweptwing.jl"
output_name = "sweptwing.md"

header = """
# Swept Wing
"""

break_flag = "COMPARISON TO EXPERIMENTAL DATA"

open(joinpath(output_path, output_name), "w") do fout

    println(fout, header)

    println(fout, """
    ```@raw html
    <center>
      <img src="../../assets/images/sweptwing-viz00.png" alt="Pic here" style="width: 100%;"/>
    </center>
    ```
    """)

    println(fout, "```julia")

    open(joinpath(pnl.examples_path, input_name), "r") do fin
        for l in eachline(fin)
            if contains(l, break_flag)
                break
            end

            println(fout, l)
        end
    end

    println(fout, "```")

    println(fout, """
    (see the complete example under
    [examples/sweptwing.jl](https://github.com/byuflowlab/FLOWPanel.jl/blob/master/examples/sweptwing.jl)
    to see how to postprocess the solution to calculate the slices of pressure
    distribution and spanwise loading that is ploted here below)

    ```@raw html
    <center>
        <br><b>Chordwise pressure distribution</b>
        <img src="../../assets/images/sweptwing000-chordpressure.png" alt="Pic here" style="width: 100%;"/>

        <br><br><b>Pressure difference</b>
        <img src="../../assets/images/sweptwing000-deltapressure.png" alt="Pic here" style="width: 100%;"/>

        <br><br><b>Spanwise loading distribution</b>
        <img src="../../assets/images/sweptwing000-loading.png" alt="Pic here" style="width: 100%;"/>
    </center>
    ```
    """)

    open(joinpath(module_path, "..", "resources", "data", "sweptwing000-CLCD.md"), "r") do fin
        for l in eachline(fin)
            if contains(l, break_flag)
                break
            end

            println(fout, l)
        end
    end
end
