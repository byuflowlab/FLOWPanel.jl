# ------------- LIFTING LINE WEBER EXAMPLE -------------------------------------

output_name = "ll-weber"
data_path = joinpath(module_path, "..", "resources", "data")




# -------- 4.2deg AOA ----------------------------------------------------------
open(joinpath(output_path, output_name*"-4p2aoa.md"), "w") do fout

    println(fout, """
    In this example we solve the aerodynamics of Weber's \$45^\\circ\$ swept-back 
    wing at an angle of attack of \$4.2^\\circ\$ using a the lifting line solver
    with a rigid wake model.
    """)

    println(fout, "# \$4.2^\\circ\$ Angle of Attack")

    println(fout, "```julia")

    input_name = "liftingline_weber.jl"
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
    [examples/liftingline_weber.jl](https://github.com/byuflowlab/FLOWPanel.jl/blob/master/examples/liftingline_weber.jl)
    to see how to postprocess the spanwise loading that is plotted below)

    ```@raw html
    <center>
        <br><br><b>Spanwise loading distribution</b>
        <img src="../../assets/images/ll-weber-loading.png" alt="Pic here" style="width: 100%;"/>
    </center>
    ```

    """)

    open(joinpath(data_path, "ll-weber-CLCD.md"), "r") do fin
        for l in eachline(fin)
            if contains(l, break_flag)
                break
            end

            println(fout, l)
        end
    end

    println(fout, """

    !!! details "Tip"
        You can also automatically run this example and generate these plots
        with the following command:
        ```julia
        import FLOWPanel as pnl

        include(joinpath(pnl.examples_path, "liftingline_weber.jl"))
        ```
    """
    )

end




# -------- AOA Sweep -----------------------------------------------------------
open(joinpath(output_path, output_name*"-aoasweep.md"), "w") do fout

    println(fout, "# AOA Sweep")

    println(fout, """
        \nUsing the wing defined in the previous section, we now sweep the angle
        of attack.
    """)

    println(fout, "```julia")

    input_name = "liftingline_weber.jl"
    start_flag = "AOA SWEEP"
    break_flag = "COMPARISON TO EXPERIMENTAL SWEEP"

    open(joinpath(pnl.examples_path, input_name), "r") do fin

        started = false

        for (li, l) in enumerate(eachline(fin))
            
            if contains(l, start_flag)
                started = true
            end

            if contains(l, break_flag)
                break
            end

            if (
                    !started || 
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
        Run time: ~5 seconds, evaluated 208 AOAs with 47% success rate. <br>
    </i></span>
    <br><br>
    ```
    (see the complete example under
    [examples/liftingline_weber.jl](https://github.com/byuflowlab/FLOWPanel.jl/blob/master/examples/liftingline_weber.jl)
    to see how to postprocess the solution as plotted here below)

    ```@raw html
    <center>
        <br><b>Spanwise loading distribution</b>
        <img src="../../assets/images/ll-weber-sweep-loading.png" alt="Pic here" style="width: 100%;"/>

        <br><br><b>Wing Polar</b><br>
        <img src="../../assets/images/ll-weber-sweep-CLCDCm.png" alt="Pic here" style="width: 100%;"/>
    </center>
    ```

    Notice that the spanwise drag distribution shows really good agreement at 
    low angles of attack but it starts to deviate at the tip for ``\\alpha \\geq 6.3^\\circ``.
    This also seen in the polar curves.
    To get more accurate ``C_D`` predictions (at the trade of a less accurate 
    spanwise drag distribution), we recommend using `use_Uind_for_force = false` as follows:
    """)


    println(fout, "```julia")

    input_name = "liftingline_weber.jl"
    start_flag = "AOA SWEEP 2"
    break_flag = "# --------- Load distribution 2"

    open(joinpath(pnl.examples_path, input_name), "r") do fin

        started = false

        for (li, l) in enumerate(eachline(fin))
            
            if contains(l, start_flag)
                started = true
            end

            if contains(l, break_flag)
                break
            end

            if (
                    !started || 
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
    <center>
        <br><b>Spanwise loading distribution</b>
        <img src="../../assets/images/ll-weber-sweep-loading2.png" alt="Pic here" style="width: 100%;"/>

        <br><br><b>Wing Polar</b><br>
        <img src="../../assets/images/ll-weber-sweep-CLCDCm2.png" alt="Pic here" style="width: 100%;"/>
    </center>
    ```
    """)


    println(fout, """

    We have added the option of superimpossing a dragging line in order to make 
    the lifting line method more robust post-stall, which is activated using 
    `sigmafactor=-1.0`. Re-running the sweep with that option increases the
    success rate from 47% to 63%.
    Furthermore, we can use `align_joints_with_Uinfs = true` to further increase
    the success rate to 78%.

    Here are the polars zoomed out to post-stall using `sigmafactor=-1.0` and
    `align_joints_with_Uinfs=true`:

    ```@raw html
    <center>
        <br><br><b>Wing Polar post-stall</b><br>
        <img src="../../assets/images/ll-weber-sweep-CLCDCm-zoomout2.png" alt="Pic here" style="width: 100%;"/>
    </center>
    ```

    ```@raw html
    <center>
        <br><br><b>sigmafactor=-1.0</b><br>
        <img src="../../assets/images/ll-weber2-sweep-CLCDCm-zoomout2.png" alt="Pic here" style="width: 100%;"/>
    </center>
    ```

    ```@raw html
    <center>
        <br><br><b>sigmafactor=-1.0, align_joints_with_Uinfs=true</b><br>
        <img src="../../assets/images/ll-weber3-sweep-CLCDCm-zoomout2.png" alt="Pic here" style="width: 100%;"/>
    </center>
    ```
    """)



    println(fout, """
    !!! details "Tip"
        You can also automatically run this example and generate these plots
        with the following command:
        ```julia
        import FLOWPanel as pnl

        include(joinpath(pnl.examples_path, "liftingline_weber.jl"))
        
        ```
    """)
end

