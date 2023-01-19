import PyPlot as plt
import LaTeXStrings: @L_str
include(joinpath(pnl.examples_path, "plotformat.jl"))

"""
Compare solution to experimental surface pressure (figures 4.8 and 6.4 in Lewis
1991)
"""
function plot_Cp(body, AOA)

    fig = plt.figure("AOA $(AOA)", figsize=[7, 5*0.8]*2/3 .* [2, 1])
    axs = fig.subplots(1, 2)

    for (ax, upperside) in zip(axs, [true, false])  # Iterate over duct side

        #=
            NOTE: Here we take a slice of the body and plot the velocity
            distribution along the slice
        =#

        # Get a slice of the body
        position        = 0.0        # Position of slice (slice along origin)
        direction       = [0, 1, 0]  # Direction of slice (slice along the xz-plane)
        row             = false      # Slices along azimuth if true; centerline if false
        # upperside     = true       # +z slice if true, -z if false

        # Data filter to get only +z or only -z data points
        slicefilter(x, i) = (x[3] >= 0) == upperside

        slicepoints, sliceCps = pnl.slicefield(body, "Cp",
                                                position, direction, row;
                                                filter=slicefilter)

        side = upperside ? "upper" : "lower"

        # Plot experimental surface pressure (figures 4.8 and 6.4 in Lewis 1991)
        if AOA==0
            fname = "lewis1991-fig4p8a.csv"
        elseif AOA==5 || AOA==15
            fname = "lewis1991-fig6p4$(side)-aoa$(ceil(Int, AOA)).csv"
        else
            error("Experimental data for AOA $(AOA) not available")
        end

        fname = joinpath(pnl.examples_path, "data", fname)
        Cp_exp = CSV.read(fname, DataFrame)

        ax.plot(Cp_exp[:, 1], Cp_exp[:, 2], "ok",
                                    markersize=5, label="Experimental")

        # Plot surface pressure of FLOWPanel
        ax.plot(slicepoints[1, :]/(d*aspectratio), sliceCps, "-", color="cyan",
                                    linewidth=2.0, alpha=0.9, label="FLOWPanel")

        # Beautify the plot
        xlims = [-0.1, 1.1]
        ylims = [1.0, -2.5]
        ax.set_xlim(xlims)
        ax.set_xticks(0:0.25:1.0)
        ax.set_ylim(ylims)
        ax.set_yticks(ylims[1]:-1.0:ylims[2])

        ax.set_xlabel(L"x/c")
        ax.set_ylabel(L"Pressure $C_p$")
        ax.legend(loc="best", frameon=false, fontsize=10)

        ax.spines["right"].set_visible(false)
        ax.spines["top"].set_visible(false)

        ax.set_title(side*" side", x=0.5, y=0.8)

    end

    fig.tight_layout()

    return fig, axs
end


function generate_fluiddomain(body, AOA, Vinf, d, aspectratio, save_path;
                                halfdomain=false, # Whether to cover only one side of the duct
                                gridname="fluidomain",
                                num=nothing,
                                verbose=true,
                                v_lvl=0
                                )

    if verbose; println("\t"^(v_lvl)*"Generating fluid domain..."); end;

    # ---------------- GENERATE FLUID DOMAIN GRID ------------------------------
    # Bounds of grid
    Pmax = d*[aspectratio*1.5,    aspectratio*0.005,  0.5*1.35] # Upper bound
    Pmin = d*[-aspectratio*0.35, -Pmax[2], -Pmax[3]]            # Lower bound

    halfdomain ? Pmin[3] = 0 : nothing

    # Grid discretization
    dx         = 0.005*d*aspectratio                            # Cell size
    dy, dz     = dx, dx

    NDIVS = ceil.(Int, (Pmax .- Pmin) ./ [dx, dy, dz]) # Divisions in each dimension

    # Generate grid
    @time grid  = pnl.gt.Grid(Pmin, Pmax, NDIVS) # Grid

    if verbose; println("\t"^(v_lvl+1)*"Grid size:\t\t$(NDIVS)"); end;
    if verbose; println("\t"^(v_lvl+1)*"Number of nodes :\t$(grid.nnodes)"); end;

    # Translate and rotate grid to align with freestream
    O = zeros(3)
    Oaxis = pnl.gt.rotation_matrix2(0, AOA, 0)
    pnl.gt.lintransform!(grid, Oaxis, O)

    # Targets where to probe the velocity
    targets = grid.nodes
    ntargets = size(targets, 2)

    # Array where to store potential and velocity
    phis = zeros(ntargets)
    Us = repeat(Vinf, 1, ntargets)

    # Calculate potential and velocity fields
    @time pnl.phi!(body, targets, phis)
    @time pnl.Uind!(body, targets, Us)

    # Save fields
    pnl.gt.add_field(grid, "phi", "scalar", phis, "node")
    pnl.gt.add_field(grid, "U", "vector", collect(eachcol(Us)), "node")

    # Output fluid domain
    @time vtks = pnl.gt.save(grid, gridname; path=save_path, num=num)

    return vtks
end
