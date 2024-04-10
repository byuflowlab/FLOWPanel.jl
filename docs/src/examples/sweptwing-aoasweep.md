# AOA Sweep
    
Using the wing defined in the previous section, we now sweep the angle
    of attack.

```julia
AOAs = [0, 2.1, 4.2, 6.3, 8.4, 10.5, 12, 14, 16] # (deg) angles of attack
Xac = [0.0*b/ar, 0, 0]                 # (m) aerodynamic center for moment calculation

# Results are stored in these arrays
Ls, Ds = [], []                         # Lift and drag at each angle of attack
Lhats, Dhats = [], []                   # Direction of lift and drag at each AOA

rolls, pitchs, yaws = [], [], []        # Rolling, pitching, and yawing moment
lhats, mhats, nhats = [], [], []        # Direction of roll, pitch, and yaw

ls, ds = [], []                         # Load and drag distributions
spanposs = []                           # Spanwise positions for load distributions


# ----------------- AOA SWEEP --------------------------------------------------
for AOA in AOAs

    Vinf = magVinf*[cos(AOA*pi/180), 0, sin(AOA*pi/180)] # Freestream

    # ----------------- CALL SOLVER --------------------------------------------
    # Freestream at every control point
    Uinfs = repeat(Vinf, 1, body.ncells)
    Das = repeat(Vinf/magVinf, 1, body.nsheddings)
    Dbs = repeat(Vinf/magVinf, 1, body.nsheddings)

    # Solve body
    @time pnl.solve(body, Uinfs, Das, Dbs)

    # ----------------- POST PROCESSING ----------------------------------------
    # Calculate velocity away from the body
    Us = pnl.calcfield_U(body, body; characteristiclength=(args...)->b/ar)

    # Calculate surface velocity U_∇μ due to the gradient of the doublet strength
    UDeltaGamma = pnl.calcfield_Ugradmu(body)

    # Add both velocities together
    pnl.addfields(body, "Ugradmu", "U")

    # Calculate pressure coeffiecient
    Cps = pnl.calcfield_Cp(body, magVinf)

    # Calculate the force of each panel
    Fs = pnl.calcfield_F(body, magVinf, rho)

    # Integrated force decomposed into lift and drag
    Dhat = Vinf/pnl.norm(Vinf)    # Drag direction
    Shat = [0, 1, 0]              # Span direction
    Lhat = pnl.cross(Dhat, Shat)  # Lift direction

    LDS = pnl.calcfield_LDS(body, Lhat, Dhat, Shat)

    L = LDS[:, 1]
    D = LDS[:, 2]

    push!(Ls, L)
    push!(Ds, D)
    push!(Lhats, Lhat)
    push!(Dhats, Dhat)

    # Integrated moment decomposed into rolling, pitching, and yawing moments
    lhat = Dhat                   # Rolling direction
    mhat = Shat                   # Pitching direction
    nhat = Lhat                   # Yawing direction

    lmn = pnl.calcfield_lmn(body, Xac, lhat, mhat, nhat)
    roll, pitch, yaw = collect(eachcol(lmn))

    push!(rolls, roll)
    push!(pitchs, pitch)
    push!(yaws, yaw)
    push!(lhats, lhat)
    push!(mhats, mhat)
    push!(nhats, nhat)

    # Calculate loading distribution
    fs, spanpos = pnl.calcfield_sectionalforce(wing_right; spandirection=[0, 1, 0])
    lds = pnl.decompose(fs, Lhat, Dhat)

    l = lds[1, :]
    d = lds[2, :]

    push!(spanposs, spanpos)
    push!(ls, l)
    push!(ds, d)
end




```
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

!!! details "Tip"
    You can also automatically run this example and generate these plots
    with the following command:
    ```julia
    import FLOWPanel as pnl

    include(joinpath(pnl.examples_path, "sweptwing.jl"))
    include(joinpath(pnl.examples_path, "sweptwing_aoasweep.jl"))
    
    ```

