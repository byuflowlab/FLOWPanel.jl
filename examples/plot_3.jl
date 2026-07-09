using CSV
using DataFrames
using Statistics
using Interpolations
using Plots

# -----------------------------
# Load simulation data
# -----------------------------
df = CSV.read("examples/wing_aileron/second.csv", DataFrame)

# Clean data
df = dropmissing(df, [:AOA, :CL, :nps])

# -----------------------------
# Experimental data
# -----------------------------
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

df_exp = DataFrame(AOA = CL_exp[:,1], CL = CL_exp[:,2])
sort!(df_exp, :AOA)

# Interpolator
itp = LinearInterpolation(df_exp.AOA, df_exp.CL, extrapolation_bc=Line())

# -----------------------------
# Compute error at each AOA
# -----------------------------
function error_by_aoa(df, itp)
    results = DataFrame(AOA=Float64[], nps=Int[], sq_error=Float64[])

    for nps in unique(df.nps)
        sim = df[df.nps .== nps, :]
        sort!(sim, :AOA)

        exp_cl = itp.(sim.AOA)
        sq_err = (sim.CL .- exp_cl).^2

        for i in eachindex(sim.AOA)
            push!(results, (
                AOA = sim.AOA[i],
                nps = nps,
                sq_error = sq_err[i]
            ))
        end
    end

    return results
end

err_df = error_by_aoa(df, itp)

using CSV
# println("Error by AOA:")
# println(first(err_df, 10))
sort!(err_df, [:nps, :AOA])
CSV.write("error_by_aoa.csv", err_df)

# -----------------------------
# Plot: error vs AOA (grouped by nps)
# -----------------------------
# plt = plot()

# for nps in unique(err_df.nps)
#     sub = err_df[err_df.nps .== nps, :]
#     sort!(sub, :AOA)

#     plot!(
#         sub.AOA,
#         sub.sq_error,
#         marker = :circle,
#         label = "nps = $nps"
#     )
# end

# xlabel!("Angle of Attack (deg)")
# ylabel!("Squared Error")
# title!("Error vs AOA by Panel Count")

# display(plt)