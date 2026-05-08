using DelimitedFiles
using Plots

dir = @__DIR__

# --- Load your CSV (skip header lines manually) ---
data = readlines(joinpath(dir, "wing_aileron/coupled_timing_results.csv"))

# Split sections
idx_coupled = findfirst(contains("BackslashCoupled"), data)
idx_iter = findfirst(contains("BackslashIterative"), data)

coupled_lines = data[idx_coupled+1 : idx_iter-1]
iter_lines    = data[idx_iter+1 : end]

# Parse into arrays
function parse_block(lines)
    aoa = Float64[]
    cl  = Float64[]
    cd  = Float64[]
    
    for line in lines
        vals = split(strip(line), ",")
        if length(vals) >= 3
            push!(aoa, parse(Float64, vals[1]))
            push!(cl,  parse(Float64, vals[2]))
            push!(cd,  parse(Float64, vals[3]))
        end
    end
    return aoa, cl, cd
end

aoa_c, cl_c, _ = parse_block(coupled_lines)
aoa_i, cl_i, _ = parse_block(iter_lines)

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

# --- Plot ---
plot(aoa_c, cl_c, label="Coupled", lw=2, grid=false)
plot!(aoa_i, cl_i, label="Iterative", lw=2, ls=:dash)
scatter!(aoa_exp, cl_exp, label="Experimental", ms=4)

xlabel!("AOA (deg)")
ylabel!("CL")
title!("AOA vs CL Comparison")

savefig("second_check.png")