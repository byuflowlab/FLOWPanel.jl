using CSV
using DataFrames
using Plots

# -------------------------------
# Load + clean data
# -------------------------------
function load_data(filepath::String)
    df = CSV.read(filepath, DataFrame)

    # Ensure correct types (helps avoid subtle bugs)
    df.AOA = Float64.(df.AOA)
    df.CL  = Float64.(df.CL)

    return df
end

# -------------------------------
# Plot AOA vs CL
# -------------------------------
function plot_aoa_cl(df::DataFrame)
    plt = plot()

    for solver in unique(df.solver)
        sub = df[df.solver .== solver, :]
        println("solver: $solver")

        # Sort by AOA for smooth curves
        sort!(sub, :AOA)

        scatter!(
            sub.AOA,
            sub.CL,
            label = solver,
            marker = :circle,
            color = :auto,
        )
    end

    xlabel!("Angle of Attack (deg)")
    ylabel!("Lift Coefficient (CL)")
    # title!("AOA vs CL Comparison")
    return plt
end

function plot_nps_cl(df::DataFrame)
    plt = plot()

    for nps in unique(df.nps)
        sub = df[df.nps .== nps, :]
        println("nps: $nps")

        # Sort by AOA for smooth curves
        sort!(sub, :AOA)

        scatter!(
            sub.AOA,
            sub.CL,
            label = "nps = $nps",
            marker = :circle,
            color = :auto
        )
    end

    xlabel!("Angle of Attack (deg)")
    ylabel!("Lift Coefficient (CL)")
    # title!("AOA vs CL Comparison by Panel Count")

    return plt
end

# -------------------------------
# Optional: Lift slope calculation
# -------------------------------
function compute_lift_slope(df::DataFrame; aoa_range=(-5,5))
    println("Lift slope (dCL/dα) near 0 AoA:")

    for solver in unique(df.solver)
        sub = df[df.solver .== solver, :]

        # Filter near linear region
        sub = sub[(sub.AOA .>= aoa_range[1]) .& (sub.AOA .<= aoa_range[2]), :]
        sort!(sub, :AOA)

        # Simple linear fit
        A = [sub.AOA ones(length(sub.AOA))]
        coeffs = A \ sub.CL
        slope = coeffs[1]

        println("  $solver → $(round(slope, digits=4)) per deg")
    end
end


function mse_by_panel_count(df::DataFrame, CL_exp::Matrix{Float64})
    # Convert experimental data → DataFrame
    df_exp = DataFrame(AOA = CL_exp[:,1], CL = CL_exp[:,2])
    sort!(df_exp, :AOA)

    # Interpolator for experimental curve
    itp = LinearInterpolation(df_exp.AOA, df_exp.CL, extrapolation_bc=Line())

    results = DataFrame(nps=Int[], mse=Float64[])

    for nps in unique(df.nps)
        sim = df[df.nps .== nps, :]
        sort!(sim, :AOA)

        # Interpolate experimental CL at sim AOA points
        exp_cl = itp.(sim.AOA)

        mse = mean((sim.CL .- exp_cl).^2)

        push!(results, (nps=nps, mse=mse))
    end

    return sort(results, :nps)
end

# --- Experimental data ---
CL_exp = [
 -9.58790170132325  -0.34195250659630627
 -7.841209829867704 -0.20211081794194818
 -5.988657844990541 -0.05435356200527597
 -4.136105860113403  0.09340369393140424
 -2.3894139886577754 0.2358839050132029
 -0.536862003780719  0.38100263852242966
  1.3156899810964084 0.5261213720316649
  3.168241965973529  0.6686015831134591
  5.02079395085066   0.8084432717678132
  6.926275992438619  0.9535620052770526
  8.672967863894144  1.0907651715039615
 10.472589792060482  1.2279683377308745
 12.325141776937613  1.3546174142480258
 14.230623818525515  1.4759894459102951
 15.077504725897917  1.528759894459108
 16.400756143667294  1.3308707124010595
 18.465028355387517  1.2807387862796877
]

aoa_exp = CL_exp[:,1]
cl_exp  = CL_exp[:,2]

# -------------------------------
# MAIN
# -------------------------------
filepath = "examples/wing_aileron/second.csv" 

df = load_data(filepath)

# Plot CL
plt1 = plot_nps_cl(df)
scatter!(plt1, aoa_exp, cl_exp, label="Experimental", marker=:square, color=:black, grid=false)
display(plt1)

# Plot CD (optional)
# plt2 = plot_aoa_cd(df)
# display(plt2)

# Compute lift slope (optional)
compute_lift_slope(df)