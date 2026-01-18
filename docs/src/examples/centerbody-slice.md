# Slice

FLOWPanel provides the following function to obtain the solution field
along a slice along a body:

```@docs
FLOWPanel.slicefield
```

Now we process the solution to plot the surface velocity along a slice
of the body of revolution.


```julia
import PythonPlot as plt
import LaTeXStrings: @L_str
include(joinpath(pnl.examples_path, "plotformat.jl"))

# ----------------- COMPARISON TO EXPERIMENTAL DATA ----------------------------
#=
    NOTE: Here we take a slice of the body and plot the velocity distribution
    along the slice.
=#

# Get a slice of the body
position        = 0.0        # Position of slice (slice along origin)
direction       = [0, 1, 0]  # Direction of slice (slice along the xz-plane)
row             = false      # If true, it slices along azimuth; centerline if false

slicepoints, sliceCps = pnl.slicefield(body, "Cp", position, direction, row)
slicepoints, sliceUs = pnl.slicefield(body, "U", position, direction, row)

# Plot experimental surface velocity distribution (figure 4.6 in Lewis 1991)
fig = plt.figure(figsize=[7, 5*0.8]*2/3)
ax = fig.gca()

filename = joinpath(pnl.examples_path, "data",
                                "centerbody-lewis-fig4p6.csv")
VoVinf_lewis = CSV.read(filename, DataFrame)

ax.plot(VoVinf_lewis[:, 1], VoVinf_lewis[:, 2], "ok",
                            markersize=5, label="Experimental")

# Plot surface velocity distribution of FLOWPanel
ax.plot(slicepoints[1, :], pnl.norm.(sliceUs)/magVinf, "-", color="cyan",
                            linewidth=2.0, alpha=0.9, label="FLOWPanel")

# Plot contour of centerbody
ax2 = ax.twinx()
xs = vcat(slicepoints[1, :], reverse(slicepoints[1, :]), slicepoints[1, 1])
ys = vcat(slicepoints[3, :], -reverse(slicepoints[3, :]), slicepoints[3, 1])
ax2.plot(xs, ys, "-k", alpha=0.25)

```
```@raw html
<center>
    <br><b>Surface velocity</b><br>
    <img src="../../assets/images/centerbody-lewis00-velocity-source.png" alt="Pic here" style="width: 60%;"/>
</center>
```

