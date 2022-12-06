import FLOWPanel as pnl

module_path = splitdir(@__FILE__)[1]      # Path to this module
output_path = joinpath(module_path, "examples") # Path where to store  examples as markdown

# ------------- SWEPT WING EXAMPLE ---------------------------------------------

output_name = "sweptwing.md"

header = """
# Swept Wing
"""


open(joinpath(output_path, output_name), "w") do fout

    println(fout, header)

    # -------- 4.2deg AOA --------------------------
    println(fout, "## \$4.2^\\circ\$ Angle of Attack")

    println(fout, """
    ```@raw html
    <center>
      <img src="../../assets/images/sweptwing-viz00.png" alt="Pic here" style="width: 100%;"/>
    </center>
    ```
    """)

    println(fout, "```julia")

    input_name = "sweptwing.jl"
    break_flag = "COMPARISON TO EXPERIMENTAL DATA"

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
    distribution and spanwise loading that is plotted here below)

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

    # -------- AOA Sweep --------------------------
    println(fout, "## AOA Sweep")

    println(fout, """
        \nUsing the wing defined in the previous section, we now sweep the angle
        of attack.
    """)

    println(fout, "```julia")

    input_name = "sweptwing_aoasweep.jl"
    break_flag = "COMPARISON TO EXPERIMENTAL DATA"

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
    [examples/sweptwing_aoasweep.jl](https://github.com/byuflowlab/FLOWPanel.jl/blob/master/examples/sweptwing_aoasweep.jl)
    to see how to postprocess the solution as plotted here below)

    ```@raw html
    <center>
        <br><b>Spanwise loading distribution</b>
        <img src="../../assets/images/sweptwing000-sweep-loading.png" alt="Pic here" style="width: 100%;"/>

        <br><br><b>Lift and induced drag</b>
        <img src="../../assets/images/sweptwing000-sweep-CLCD.png" alt="Pic here" style="width: 100%;"/>

        <br><br><b>Pitching moment</b><br>
        <img src="../../assets/images/sweptwing000-sweep-Cm.png" alt="Pic here" style="width: 50%;"/>
    </center>
    ```
    """)
end
