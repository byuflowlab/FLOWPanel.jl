## Postprocessing plots for the rotor_hover and rotor_Ground examples
# Author : Timothy Harlow
# Created : April 15, 2026

using DelimitedFiles
using Plots
using LaTeXStrings

CT_exp = 0.072 # Experimental CT
save_path = "/Users/Timothy/Library/CloudStorage/OneDrive-BrighamYoungUniversity/Capstone/Plots/FLOWPanel"

## Convergence Study
p_per_step = [1, 2, 5, 10, 20] # Particles per step
p_per_rev  = p_per_step * 36  # nt = 36 for these runs

data_p1 = readdlm("CT_v_t_hover_RPM5400_pps1nt36_overlap2.0_kerneloff0.00012.csv", ',')
data_p2 = readdlm("CT_v_t_hover_RPM5400_pps2_nt36_overlap4.0_kerneloff0.00012.csv", ',')
data_p3 = readdlm("CT_v_t_hover_RPM5400_pps5_nt36_overlap10.0_kerneloff0.00012.csv", ',')
data_p4 = readdlm("CT_v_t_hover_RPM5400_pps10_nt36_overlap20.0_kerneloff0.00012.csv", ',')
data_p5 = readdlm("CT_v_t_hover_RPM5400_pps20_nt36_overlap40.0_kerneloff0.00012.csv", ',')

t_p = data_p1[2, :]
CT_data_p = [data_p1[2,:] data_p2[2,:] data_p3[2,:] data_p4[2,:] data_p5[2,:]]

di = 36 * 2 # `nt` for these runs
CT_p = sum(CT_data_p[end-di+1 : end, :]; dims=1) / length(t_p[end-di+1 : end])
# percent_error_p = abs.((CT_p .- CT_exp) ./ CT_exp) * 100

plt_CT_ppr = plot(p_per_rev, CT_p';
    title   = "Coefficient of Thrust "*L"(n_t = 36)",
    xlabel  = "No. of particles per revolution "*L"n_p",
    ylabel  = L"C_T",
    # ylims   = [0.067, 0.07],
    markers = :circle,
    label   = L"n_p")
# savefig(plt_CT_ppr, savepath*"/ppr_ct.png")

plt_conv_pps = plot(1:length(t_p), CT_data_p[:,3];
    title  = L"n_p = 10",
    xlabel = "time step",
    ylabel = L"C_T",
    label  = "VPM")
plot!(plt_conv_pps, 1:length(t_p), fill(CT_exp, length(t_p)),
    label  = "Experiment", linestyle=:dash)
# savefig(plt_conv_pps, savepath*"/pps_convergence.png")

## Kernel offset
kerneloffset = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2] # Kernel offset

data_k1 = readdlm("CT_v_t_hover_RPM5400_pps2_nt36_overlap2.0_kerneloff0.0012.csv", ',')
data_k2 = readdlm("CT_v_t_hover_RPM5400_pps2_nt36_overlap2.0_kerneloff0.00012.csv", ',')
data_k3 = readdlm("CT_v_t_hover_RPM5400_pps2_nt36_overlap2.0_kerneloff1.2e-5.csv", ',')
data_k4 = readdlm("CT_v_t_hover_RPM5400_pps2_nt36_overlap2.0_kerneloff1.2e-6.csv", ',')
data_k5 = readdlm("CT_v_t_hover_RPM5400_pps2_nt36_overlap2.0_kerneloff1.2e-7.csv", ',')

t_k = data_k1[2, :]
CT_data_k = [data_k1[2,:] data_k2[2,:] data_k3[2,:] data_k4[2,:] data_k5[2,:]]

di = 36 * 2 # `nt` for these runs
CT_k = sum(CT_data_k[end-di+1 : end, :]; dims=1) / length(t_k[end-di+1 : end])
# percent_error_p = abs.((CT_p .- CT_exp) ./ CT_exp) * 100

plt_CT_ppr = plot(kerneloffset, CT_k';
    title   = "Coefficient of Thrust "*L"(n_t = 36, np = 2)",
    xlabel  = "Kernel offset • R ",
    xscale  = :log10,
    ylabel  = L"C_T",
    ylims   = [0.0095, 0.0115],
    markers = :circle,
    label   = L"kernel offset")
# savefig(plt_CT_ppr, savepath*"/ppr_ct.png")

plt_conv_pps = plot(1:length(t_k), CT_data_k[:,3];
    title  = L"k_o = 10",
    xlabel = "time step",
    ylabel = L"C_T",
    label  = "VPM")
plot!(plt_conv_pps, 1:length(t_p), fill(CT_exp, length(t_p)),
    label  = "Experiment", linestyle=:dash)
# savefig(plt_conv_pps, savepath*"/pps_convergence.png")