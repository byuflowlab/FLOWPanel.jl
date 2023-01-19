import PyPlot as plt
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

# Beautify the plot
xlims = [-0.1, 1.5]
ylims = [-0.5, 1.5]
ax.set_xlim(xlims)
ax.set_xticks(0:0.5:1.5)
ax.set_ylim(ylims)
ax.set_yticks(0:0.5:1.5)

ax2.set_aspect(1.0)
ax2.set_ylim([-0.62, 1])
ax2.set_yticks([])
ax2.plot(xlims, zeros(2), ":k", alpha=0.25, linewidth=1)

ax.set_xlabel(L"$x$-position (m)")
ax.set_ylabel(L"Velocity $u/u_\infty$")
ax.legend(loc="best", frameon=false, fontsize=8)

for a in [ax, ax2]
    a.spines["right"].set_visible(false)
    a.spines["top"].set_visible(false)
end

fig.tight_layout()
