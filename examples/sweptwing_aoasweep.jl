#=##############################################################################
# DESCRIPTION
    AOA sweep on 45deg swept-back wing.

# AUTHORSHIP
  * Author    : Eduardo J. Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Dec 2022
  * License   : MIT License
=###############################################################################


AOAs = [0, 2.1, 4.2, 6.3, 8.4, 10.5, 12, 14, 16] # (deg) angles of attack
Xac = [0.25*b/ar, 0, 0]                 # (m) aerodynamic center for moment calculation

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




# ----------------- COMPARISON TO EXPERIMENTAL DATA ----------------------------

# --------- Load distribution
nondim = 0.5*rho*magVinf^2*b/ar   # Normalization factor

cls = ls / nondim
cds = ds / nondim

fig4 = plt.figure(figsize=[7*2, 5*1*0.8]*2/3)
axs = fig4.subplots(1, 2)


for (axi, (ax, vals_exp)) in enumerate(zip(axs, [cls_web, cds_web]))

    first = true

    for (AOA, pos, cl, cd) in zip(AOAs, spanposs, cls, cds)
        rowi = findfirst(a -> a==AOA, alphas_web)
        if rowi != nothing && AOA in (axi==1 ? [2.1, 4.2, 6.3, 8.4] : [4.2, 6.3, 8.4])

            # Filter out NaNs
            ys = vals_exp[rowi, :]
            xs = [val for (vali, val) in enumerate(y2b_web) if !isnan(ys[vali])]
            ys = [val for (vali, val) in enumerate(ys) if !isnan(ys[vali])]

            # Plot experimental
            for f in [-1, 1]
                ax.plot(f*xs, ys, "o--k",
                            label=("Experimental"^(f==1))^first,
                            linewidth=0.5, markersize=5, alpha=1.0)
            end

            # Plot FLOWPanel
            ax.plot(pos*2/b, axi==1 ? cl : cd, "-", label="FLOWPanel"^first,
                            color="steelblue", markersize=8, linewidth=1)

            first = false
        end

    end

    xlims = [0, 1]
    ax.set_xlim(xlims)
    ax.set_xticks(xlims[1]:0.2:xlims[end])
    ax.set_xlabel(L"Span position $2y/b$")

    if axi==1
        ylims = [0, 0.6]
        ax.set_ylim(ylims)
        ax.set_yticks(ylims[1]:0.2:ylims[end])
        ax.set_ylabel(L"Sectional lift $c_\ell$")

        ax.legend(loc="best", frameon=false, fontsize=6)

        ax.annotate(L"\alpha=8.4^\circ", [0.38, 0.53], xycoords="data", fontsize=ANNOT_SIZE, color="black", alpha=0.6)
        ax.annotate(L"\alpha=6.3^\circ", [0.38, 0.402], xycoords="data", fontsize=ANNOT_SIZE, color="black", alpha=0.6)
        ax.annotate(L"\alpha=4.2^\circ", [0.38, 0.275], xycoords="data", fontsize=ANNOT_SIZE, color="black", alpha=0.6)
        ax.annotate(L"\alpha=2.1^\circ", [0.38, 0.145], xycoords="data", fontsize=ANNOT_SIZE, color="black", alpha=0.6)
    else
        ylims = [-0.04, 0.12]
        ax.set_ylim(ylims)
        ax.set_yticks(ylims[1]:0.04:ylims[end])
        ax.set_ylabel(L"Sectional drag $c_d$")


        ax.annotate(L"\alpha=8.4^\circ", [0.25, 0.030], xycoords="data", fontsize=ANNOT_SIZE, color="black", alpha=0.6, rotation=-10)
        ax.annotate(L"\alpha=4.2^\circ", [0.25, -0.005], xycoords="data", fontsize=ANNOT_SIZE, color="black", alpha=0.6, rotation=-5)

        ax.annotate(L"\alpha=6.3^\circ", [0.5, 0.035], xycoords="data", fontsize=ANNOT_SIZE, color="black", alpha=0.6)

        ax.annotate("", [0.4, 0.0145], xycoords="data",
                    xytext=[0.5, 0.035], textcoords="data",
                    arrowprops=Dict(:facecolor=>"black", :linewidth=>0, :alpha=>0.4,
                                    :shrink=>0, :width=>1.0, :headwidth=>5.0, :headlength=>7))
    end

    ax.spines["right"].set_visible(false)
    ax.spines["top"].set_visible(false)
end

fig4.tight_layout()

# --------- Integrated forces: lift and drag
nondim = 0.5*rho*magVinf^2*b^2/ar   # Normalization factor

CLs = sign.(dot.(Ls, Lhats)) .* norm.(Ls) / nondim
CDs = sign.(dot.(Ds, Dhats)) .* norm.(Ds) / nondim

# VSPAERO CL and CD
data_vsp = CSV.read(vsp_file, DataFrame; skipto=397, limit=419-397+1)
alphas_vsp = [val for val in data_vsp[1, 2:end]]
CDi_vsp = [val for val in data_vsp[3, 2:end]]
CDtot_vsp = [val for val in data_vsp[6, 2:end]]
CL_vsp = [val for val in data_vsp[11, 2:end]]
CMy_vsp = [val for val in data_vsp[16, 2:end]]

fig5 = plt.figure(figsize=[7*2, 5*1*0.75]*2/3)
axs = fig5.subplots(1, 2)

ax = axs[1]
ax.plot(alphas_web, CLs_web, "-ok", label="Experimental")
ax.plot(alphas_vsp, CL_vsp, ":v", color=color_vsp, alpha=0.9, label="VSPAERO")
ax.plot(AOAs, CLs, ":^", label="FLOWPanel", color="steelblue", markersize=8, alpha=0.7)

ylims = [0, 0.8]
ax.set_ylim(ylims)
ax.set_yticks(ylims[1]:0.2:ylims[end])
ax.set_ylabel(L"Lift coefficient $C_L$")

ax.legend(loc="lower right", fontsize=10, frameon=false, reverse=true)

ax = axs[2]
ax.plot(alphas_web, CDs_web, "-ok", label="Experimental")
ax.plot(alphas_vsp, CDtot_vsp, ":v", color=color_vsp, alpha=0.9, label="VSPAERO")
ax.plot(AOAs, CDs, ":^", label="FLOWPanel", color="steelblue", markersize=8, alpha=0.7)

ylims = [0, 0.04]
ax.set_ylim(ylims)
ax.set_yticks(ylims[1]:0.01:ylims[end])
ax.set_ylabel(L"Drag coefficient $C_D$")

for ax in axs
    xlims = [0, 12]
    xticks = xlims[1]:2:xlims[2]
    ax.set_xlim(xlims)
    ax.set_xticks(xticks)
    ax.set_xticklabels(["$val"*L"^\circ" for val in xticks])
    ax.set_xlabel(L"Angle of attack $\alpha$")

    ax.spines["right"].set_visible(false)
    ax.spines["top"].set_visible(false)
end

fig5.tight_layout()

# --------- Integrated moment: Pitching moment
nondim = 0.5*rho*magVinf^2*b*(b/ar)^2 # Normalization factor

Cls = sign.(dot.(rolls, lhats)) .* norm.(rolls) / nondim
Cms = sign.(dot.(pitchs, mhats)) .* norm.(pitchs) / nondim
Cns = sign.(dot.(yaws, nhats)) .* norm.(yaws) / nondim

fig6 = plt.figure(figsize=[7*1, 5*1*0.75]*2/3)
ax = fig6.gca()

ax.plot(alphas_vsp, CMy_vsp, ":v", color=color_vsp, alpha=0.9, label="VSPAERO")
ax.plot(AOAs, Cms, ":o", label="FLOWPanel", color="steelblue", markersize=8, alpha=0.7)

xlims = [0, 16]
xticks = xlims[1]:2:xlims[2]
ax.set_xlim(xlims)
ax.set_xticks(xticks)
ax.set_xticklabels(["$val"*L"^\circ" for val in xticks])
ax.set_xlabel(L"Angle of attack $\alpha$")

ylims = [-1.2, 0.2]
ax.set_ylim(ylims)
ax.set_yticks(ylims[1]:0.2:ylims[end])
ax.set_ylabel(L"Pitching moment $C_m$")

ax.spines["right"].set_visible(false)
ax.spines["top"].set_visible(false)

ax.legend(loc="best", frameon=false, fontsize=10, reverse=true)

fig6.tight_layout()


# --------- Save figures
if save_outputs
    fig4.savefig(joinpath(fig_path, "$(run_name)-sweep-loading.png"),
                                                dpi=300, transparent=true)
    fig5.savefig(joinpath(fig_path, "$(run_name)-sweep-CLCD.png"),
                                                dpi=300, transparent=true)
    fig6.savefig(joinpath(fig_path, "$(run_name)-sweep-Cm.png"),
                                                dpi=300, transparent=true)
end
