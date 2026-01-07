# AOA Sweep
    
Using the wing defined in the previous section, we now sweep the angle
    of attack.

```julia
# ----------------- AOA SWEEP --------------------------------------------------

# Sequence of sweeps to run
# NOTE: To help convergence and speed it up the sweep, we recommend starting
#       each sweep at 0deg AOA since the sweep steps over points using the last
#       converged solution as the initial guess
sweep1 = [0, 2.1, 4.2, 6.3, 8.4, 10.5]      # Same AOAs than Weber's experiment
sweep2 = range(0, -50, step=-0.5)           # Sweep from 0 into deep negative stall (-50deg)
sweep3 = range(0, 50, step=0.5)             # Sweep from 0 into deep positive stall (50deg)


@time wingpolar = pnl.run_polarsweep(ll,
                            magUinf, rho, X0, cref, b;
                            aoa_sweeps = (sweep1, sweep2, sweep3),
                            solver,
                            solver_optargs,
                            align_joints_with_Uinfs, 
                            use_Uind_for_force
                        )


```
(see the complete example under
[examples/liftingline_weber.jl](https://github.com/byuflowlab/FLOWPanel.jl/blob/master/examples/liftingline_weber.jl)
to see how to postprocess the solution as plotted here below)

```@raw html
<center>
    <br><b>Spanwise loading distribution</b>
    <img src="../../assets/images/ll-weber-sweep-loading.png" alt="Pic here" style="width: 100%;"/>

    <br><br><b>Wing Polar</b><br>
    <img src="../../assets/images/ll-weber-sweep.png" alt="Pic here" style="width: 100%;"/>
</center>
```

!!! details "Tip"
    You can also automatically run this example and generate these plots
    with the following command:
    ```julia
    import FLOWPanel as pnl

    include(joinpath(pnl.examples_path, "liftingline_weber.jl"))
    
    ```

