# ------------- SWEPT WING EXAMPLE ---------------------------------------------

output_name = "sweptwing"
data_path = joinpath(module_path, "..", "resources", "data")




# -------- 4.2deg AOA ----------------------------------------------------------
open(joinpath(output_path, output_name*"-4p2aoa.md"), "w") do fout

    println(fout, """
    ```@raw html
    <center>
      <img src="../../assets/images/sweptwing-viz00.png" alt="Pic here" style="width: 100%;"/>
    </center>
    ```

    In this example we solve the flow around a \$45^\\circ\$ swept-back wing at an
    angle of attack of \$4.2^\\circ\$ using a rigid wake model.
    """)

    println(fout, "# \$4.2^\\circ\$ Angle of Attack")

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

    open(joinpath(data_path, "sweptwing000-CLCD.md"), "r") do fin
        for l in eachline(fin)
            if contains(l, break_flag)
                break
            end

            println(fout, l)
        end
    end

end




# -------- AOA Sweep -----------------------------------------------------------
open(joinpath(output_path, output_name*"-aoasweep.md"), "w") do fout

    println(fout, "# AOA Sweep")

    println(fout, """
        \nUsing the wing defined in the previous section, we now sweep the angle
        of attack.
    """)

    println(fout, "```julia")

    input_name = "sweptwing_aoasweep.jl"
    break_flag = "COMPARISON TO EXPERIMENTAL DATA"

    open(joinpath(pnl.examples_path, input_name), "r") do fin
        for (li, l) in enumerate(eachline(fin))
            if contains(l, break_flag)
                break
            end

            if li>=13
                println(fout, l)
            end
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




# -------- Linear solver benchmark ---------------------------------------------
open(joinpath(output_path, output_name*"-solver.md"), "w") do fout

    println(fout, "# Solver Benchmark")

    println(fout, """
    \nThe problem of solving the panel strengths that satisfy the
    no-flow-through condition poses a linear system of equation of the form

    ```math
    \\begin{align*}
            A y = b
    ,\\end{align*}
    ```

    where \$A\$ is the matrix containing the geometry of the panels and wake,
    \$y\$ is the vector of panel strengths, and \$b\$ is the vector of boundary
    conditions.
    This is trivially solved as

    ```math
    \\begin{align*}
            y = A^{-1}b
    ,\\end{align*}
    ```

    however, depending on the size of \$A\$ (which depends on the number of
    panels) it can become inefficient or even unfeasible to explicitely
    calculate the inverse of \$A\$.
    Multiple linear solvers are available in FLOWPanel that avoid
    explicitely inverting \$A\$, which are described and benchmarked as
    follows.

    > The complete code is available at """*"""
    [examples/sweptwing_solverbenchmark.jl](https://github.com/byuflowlab/FLOWPanel.jl/blob/master/examples/sweptwing_solverbenchmark.jl) """*"""
    but you should also be able to copy and paste these lines after running the """*"""
    first section of this example.

    """)



    println(fout, """
    ## Backslash operator `\\`

    Most programming languages implement an operator `\\` that directly
    calculates the matrix-vector product \$A^{-1}b\$.
    This is more efficient than directly inverting \$A\$ and then multiplying
    by \$b\$, without loosing any accuracy.
    This linear solver is available under this function:
    ```@docs
    FLOWPanel.solve_backslash!
    ```
    and is used as follows:
    """)


    println(fout, "```julia")
    open(joinpath(pnl.examples_path, "sweptwing_solverbenchmark.jl"), "r") do fin
        for (li, l) in enumerate(eachline(fin))
            if 12 <= li <= 65+6
                println(fout, l)
            end
        end
    end
    println(fout, "```")

    println(fout, "\n")
    open(joinpath(data_path, "sweptwing000-backslash.md"), "r") do fin
        for l in eachline(fin)
            println(fout, l)
        end
    end



    println(fout, """
    ## LU decomposition

    Pre-calculating and re-using the
    [LU decomposition](https://en.wikipedia.org/wiki/LU_decomposition) of \$A\$
    is advantageous when the linear system needs to be solved for multiple
    boundary conditions.

    ```@docs
    FLOWPanel.solve_ludiv!
    ```

    Running the solver:
    """)


    println(fout, "```julia")
    open(joinpath(pnl.examples_path, "sweptwing_solverbenchmark.jl"), "r") do fin
        for (li, l) in enumerate(eachline(fin))
            if 77+6 <= li <= 79+6
                println(fout, l)
            end
        end
    end
    println(fout, "```")

    open(joinpath(data_path, "sweptwing000-ludiv.md"), "r") do fin
        for l in eachline(fin)
            println(fout, l)
        end
    end



    println(fout, """
    ## Iterative Krylov Solver

    Iterative Krylov subspace solvers converge to the right solution
    rather than directly solving the system of equations.
    This allows the user to trade off accuracy for computational speed by
    tuning the tolerance of the solver.

    The [generalized minimal residual](https://en.wikipedia.org/wiki/Generalized_minimal_residual_method)
    (GMRES) method provided by
    [Krylov.jl](https://juliasmoothoptimizers.github.io/Krylov.jl/stable/solvers/unsymmetric/#GMRES)
    is available in FLOWPanel.

    ```@docs
    FLOWPanel.solve_gmres!
    ```
    """)


    println(fout, "Running the solver with tolerance \$10^{-8}\$:")
    println(fout, "```julia")
    open(joinpath(pnl.examples_path, "sweptwing_solverbenchmark.jl"), "r") do fin
        for (li, l) in enumerate(eachline(fin))
            if 91+6 <= li <= 97+6
                println(fout, l)
            end
        end
    end
    println(fout, "```")
    open(joinpath(data_path, "sweptwing000-gmres8.md"), "r") do fin
        for l in eachline(fin)
            println(fout, l)
        end
    end

    println(fout, "Running the solver with tolerance \$10^{-2}\$:")
    println(fout, "```julia")
    open(joinpath(pnl.examples_path, "sweptwing_solverbenchmark.jl"), "r") do fin
        for (li, l) in enumerate(eachline(fin))
            if 110+6 <= li <= 116+6
                println(fout, l)
            end
        end
    end
    println(fout, "```")
    open(joinpath(data_path, "sweptwing000-gmres2.md"), "r") do fin
        for l in eachline(fin)
            println(fout, l)
        end
    end

    println(fout, "Running the solver with tolerance \$10^{-1}\$:")
    println(fout, "```julia")
    open(joinpath(pnl.examples_path, "sweptwing_solverbenchmark.jl"), "r") do fin
        for (li, l) in enumerate(eachline(fin))
            if 130+6 <= li <= 136+6
                println(fout, l)
            end
        end
    end
    println(fout, "```")
    open(joinpath(data_path, "sweptwing000-gmres1.md"), "r") do fin
        for l in eachline(fin)
            println(fout, l)
        end
    end

end
