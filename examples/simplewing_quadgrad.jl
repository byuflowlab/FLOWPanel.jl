#=##############################################################################
# DESCRIPTION
#   Straight wing at an angle of attack of 8.4 deg. This example solves the
#   triangular-panel simple wing once, then compares Bernoulli pressure and
#   loading from two doublet-gradient surface-velocity reconstructions:
#   the raw triangle stencil and the standalone quad agglomerate μ-difference
#   stencil.
#
# AUTHORSHIP
#   * Created   : Jun 2026
#   * License   : MIT License
=###############################################################################

import FLOWPanel as pnl
include(joinpath(pnl.examples_path, "helper_functions.jl"))
using GeometricTools: PyPlot as plt
import LinearAlgebra
import Printf
import Statistics
import WriteVTK

const _norm = LinearAlgebra.norm
const _dot = LinearAlgebra.dot
const _cross = LinearAlgebra.cross

envflag(name, default=false) =
    lowercase(get(ENV, name, string(default))) in ("1", "true", "yes", "on")

airfoil_path = joinpath(pnl.examples_path, "data")

save_vtk = envflag("FLOWPANEL_SIMPLEWING_QUADGRAD_VTK", true)
save_plots = envflag("FLOWPANEL_SIMPLEWING_QUADGRAD_PLOTS", true)
snap_cp_slices = envflag("FLOWPANEL_SIMPLEWING_QUADGRAD_SNAP_CP", true)

# ----------------- SIMULATION PARAMETERS --------------------------------------
AOA = 8.4
aoa_tag = replace(string(AOA), "." => "p")
magVinf = 30.0
Vinf = magVinf * [cos(AOA * pi / 180), 0, sin(AOA * pi / 180)]
rho = 1.225

# ----------------- GEOMETRY DESCRIPTION ---------------------------------------
b = 98 * 0.0254
ar = 5.0
c = b / ar
tr = 1.0
twist_root = 0
twist_tip = 0
lambda = 0
gamma = 0
airfoil = "airfoil-rae101.csv"

n_rfl = 16
NDIVS_rfl = [(0.25, n_rfl, 10.0, false),
             (0.50, n_rfl, 1.0, true),
             (0.25, n_rfl, 1 / 10.0, false)]

n_span = 60
NDIVS_span = [(1.0, n_span, 1.0, false)]
n_loading_bins = 2 * n_span

# ----------------- GENERATE BODY ----------------------------------------------
println("Generating wing...")

bodyoptargs = (;
    kerneloffset=1e-6,
    kernelcutoff=1e-12,
    semiinfinite_wake=true,
    DBC=false,
)

kernel = pnl.VortexRing
bodytype = pnl.RigidWakeBody{kernel}
run_name = "simplewing_quadgrad"
save_path = joinpath("data", run_name)

@time body = simplewing_mirrored(b, ar, tr, twist_root, twist_tip, lambda, gamma;
                                 bodytype=bodytype, bodyoptargs=bodyoptargs,
                                 airfoil_root=airfoil, airfoil_tip=airfoil,
                                 airfoil_path=airfoil_path,
                                 rfl_NDIVS=NDIVS_rfl,
                                 delim=",",
                                 span_NDIVS=NDIVS_span,
                                 verify_spline=false,
                                 verify_rflspline=false)
@show typeof(body)

wake_direction = reshape(Vinf ./ magVinf, :, 1)
for i in eachindex(body.Das)
    body.Das[i] .= repeat(wake_direction, 1, size(body.Das[i], 2))
end

println("Number of panels:\t$(body.ncells)")

# ----------------- SOLVE ONCE --------------------------------------------------
println("Solving body once...")

backend = pnl.FastMultipoleBackend()
Dhat = Vinf / _norm(Vinf)
Shat = [0, 1, 0]
Lhat = _cross(Dhat, Shat)

Sref = b^2 / ar
c_ref = b / ar
force_q = 0.5 * rho * magVinf^2 * Sref
sectional_qc = 0.5 * rho * magVinf^2 * c_ref

frames = pnl.ReferenceFrame(body)
@time pnl.steady!(body, frames, Vinf;
    body_solvers=pnl.Backslash(body),
    backend=backend,
    monitors=(),
    path=nothing,
    verbose=false)

# ----------------- RECONSTRUCT AND POSTPROCESS --------------------------------
function set_wake_direction!(body_case, wake_direction)
    for i in eachindex(body_case.Das)
        body_case.Das[i] .= repeat(wake_direction, 1, size(body_case.Das[i], 2))
    end
    return body_case
end

function postprocess_case(solved_body, label, tag; grad_mu_options)
    body_case = set_wake_direction!(deepcopy(solved_body), wake_direction)
    pnl.calcfield_U!(body_case, Vinf; backend=backend,
                     grad_mu_options=grad_mu_options)

    pressure = zeros(body_case.ncells)
    pnl.calcfield_P!(pressure, body_case, body_case.velocity, magVinf, rho, nothing;
                     correct_kuttacondition=false)

    force = zeros(3, body_case.ncells)
    pnl.calcfield_F!(force, body_case, pnl.calc_areas(body_case), body_case.normals, pressure;
                     correct_kuttacondition=false)
    total_force = pnl.calcfield_Ftot(body_case, force)
    Cp = pnl.calcfield_Cp(pressure, magVinf, rho)
    CL = _dot(total_force, Lhat) / force_q
    CD = _dot(total_force, Dhat) / force_q
    return (; label, tag, body=body_case, pressure, Cp, force, total_force, CL, CD)
end

# Two surface-velocity reconstructions on the SAME γ solve: raw triangle stencil
# and standalone quad-pitch μ-difference. The solve is not repeated — only the
# post-solve ∇μ recovery differs.
#
# Every knob that affects the ∇μ reconstruction is spelled out explicitly below; none
# is left to a silent default. Notes:
#   * `scale` (the ±½∇μ half-jump weight) is fixed at 0.5 inside calcfield_U! and is
#     not exposed at this call site.
#   * `basis=:tri` uses the raw per-triangle 1-ring stencil; the quad pairing/growth
#     options below are then inert (shown for completeness/visibility).
#   * triangle robustness is OFF (tri_robust=false), so its donor target is inert.
tri_grad_mu_options = (;
    basis = :tri,
    tri_robust = false,
    tri_robust_ar_threshold = 10.0,
    tri_robust_max_depth = 4,
    tri_robust_target_healthy = 6,
)
quad_grad_mu_options = (;
    basis = :quad,
    quad_grow = false,                 # agglomerate-stencil growth OFF
    quad_grow_stop = :cond,            # if grown: stop on conditioning (:cond) vs fixed rings (:depth)
    quad_grow_cond_max = 1e3,          # LS condition-number threshold for :cond
    quad_grow_max_depth = 4,           # ring cap for quad agglomerate growth
    quad_normal_dot_min = cos(pi / 4),  # fold gate: never pair/difference across a sharper crease
)
cases = [
    postprocess_case(body, "Tri gradient", "tri";
        grad_mu_options=tri_grad_mu_options),
    postprocess_case(body, "Quad mu_diff", "quad_mudiff";
        grad_mu_options=quad_grad_mu_options),
]
case_tri = cases[1]

struct QuadGradVTKFields
    pressure::Vector{Float64}
    Cp::Vector{Float64}
    force::Matrix{Float64}
end

function pnl.write_vtk_fields!(vtk, m::QuadGradVTKFields, body, i_system::Int,
                               i_step::Int, field_names::pnl.VTKFieldNameAllocator,
                               i_monitor::Int)
    vtk["pressure", WriteVTK.VTKCellData()] = m.pressure
    vtk["Cp", WriteVTK.VTKCellData()] = m.Cp
    vtk["F", WriteVTK.VTKCellData()] = m.force
    return nothing
end

function finite_check(case)
    all(isfinite, case.body.velocity) || error("Non-finite velocity in $(case.label).")
    all(isfinite, case.pressure) || error("Non-finite pressure in $(case.label).")
    all(isfinite, case.force) || error("Non-finite force in $(case.label).")
    return nothing
end

foreach(finite_check, cases)

function column_norms(A)
    return [_norm(view(A, :, i)) for i in axes(A, 2)]
end

function relative_norm_delta(a, b)
    denom = max(_norm(b), eps(eltype(float.(b))))
    return _norm(a .- b) / denom
end

function te_jump_summary(case)
    Ujumps = Float64[]
    Cpjumps = Float64[]
    for shedding in case.body.shedding
        for (pi, _nia, _nib, pj, _nja, _njb) in eachcol(shedding)
            pj == -1 && continue
            push!(Ujumps, _norm(case.body.velocity[:, pi] .- case.body.velocity[:, pj]))
            push!(Cpjumps, abs(case.Cp[pi] - case.Cp[pj]))
        end
    end
    isempty(Ujumps) && return (; n=0, mean_U=NaN, max_U=NaN, mean_Cp=NaN, max_Cp=NaN)
    return (; n=length(Ujumps),
            mean_U=Statistics.mean(Ujumps), max_U=maximum(Ujumps),
            mean_Cp=Statistics.mean(Cpjumps), max_Cp=maximum(Cpjumps))
end

println("\n#===== INTEGRATED FORCE COEFFICIENTS =====#")
for case in cases
    Printf.@printf("%-16s CL = %.8f   CD = %.8f\n", case.label, case.CL, case.CD)
end

# Pairwise comparisons: each non-tri reconstruction vs the raw triangle baseline,
# plus mu_diff vs grad_avg (the two quad recoveries against each other).
comparisons = [(cases[k], case_tri) for k in 2:length(cases)]
length(cases) >= 3 && push!(comparisons, (cases[3], cases[2]))

println("\n#===== FIELD DIFFERENCES (a vs b) =====#")
for (a, b) in comparisons
    dU_colnorm = column_norms(a.body.velocity .- b.body.velocity)
    Printf.@printf("%-16s vs %-16s  rel||dU||=%.6e  max|dU|=%.6e m/s  rel||dCp||=%.6e  max|dCp|=%.6e\n",
                   a.label, b.label,
                   relative_norm_delta(a.body.velocity, b.body.velocity),
                   maximum(dU_colnorm),
                   relative_norm_delta(a.Cp, b.Cp),
                   maximum(abs.(a.Cp .- b.Cp)))
end

println("\n#===== TRAILING-EDGE ADJACENT JUMPS =====#")
for case in cases
    te = te_jump_summary(case)
    Printf.@printf("%-16s n = %d   mean |dU_TE| = %.8e   max |dU_TE| = %.8e   mean |dCp_TE| = %.8e   max |dCp_TE| = %.8e\n",
                   case.label, te.n, te.mean_U, te.max_U, te.mean_Cp, te.max_Cp)
end

function spanwise_loading(case, nbins, b, sectional_qc)
    cps = case.body.controlpoints
    span_coord = vec(Shat' * cps)
    smin, smax = minimum(span_coord), maximum(span_coord)
    width = (smax - smin) / nbins
    bin_force = zeros(3, nbins)
    counts = zeros(Int, nbins)
    for p in 1:case.body.ncells
        ibin = clamp(floor(Int, (span_coord[p] - smin) / width) + 1, 1, nbins)
        bin_force[:, ibin] .+= case.force[:, p]
        counts[ibin] += 1
    end

    populated = findall(!=(0), counts)
    isempty(populated) && error("No populated spanwise loading bins for $(case.label).")
    for ibin in 1:nbins
        counts[ibin] == 0 || continue
        nearest = populated[argmin(abs.(populated .- ibin))]
        bin_force[:, ibin] .= bin_force[:, nearest]
    end

    centers = [smin + (ibin - 0.5) * width for ibin in 1:nbins]
    cl = [_dot(bin_force[:, ibin] ./ width, Lhat) / sectional_qc for ibin in 1:nbins]
    cd = [_dot(bin_force[:, ibin] ./ width, Dhat) / sectional_qc for ibin in 1:nbins]
    return (; y2b=2 .* centers ./ b, cl, cd)
end

loadings = [spanwise_loading(case, n_loading_bins, b, sectional_qc) for case in cases]

# ----------------- VISUALIZATION ----------------------------------------------
if save_vtk
    mkpath(save_path)
    for case in cases
        pnl.write_vtk(joinpath(save_path, run_name * "_$(case.tag)_AOA$(aoa_tag)"),
                      case.body, 0, 0.0;
                      monitors=(QuadGradVTKFields(case.pressure, case.Cp, case.force),),
                      overwrite=true)
    end
    println("Wrote VTK files to $(save_path)/")
end

function place_figure_legend_right!(fig, axs)
    handles = Any[]
    labels = String[]
    for ax in axs
        axis_handles, axis_labels = ax.get_legend_handles_labels()
        for (handle, label) in zip(axis_handles, axis_labels)
            label == "_nolegend_" && continue
            label in labels && continue
            push!(handles, handle)
            push!(labels, label)
        end

        legend = ax.get_legend()
        legend !== nothing && legend.remove()
    end

    isempty(handles) && return nothing
    fig.legend(handles, labels; loc="center left", bbox_to_anchor=(0.80, 0.5),
               fontsize=10, frameon=false)
    return nothing
end

function plot_Cp_overlays(cases, spanposs, b, c;
                          slicetol=0.03 * b,
                          snap_to_span_row=false,
                          dim_span=2,
                          stls=["-", "--"],
                          xlims=[-0.1, 1.1],
                          ylims=[1.0, -1.5])
    npos = length(spanposs)
    fig = plt.figure("simplewing_quadgrad_chordwise_loading",
                     figsize=(9.0, 3.2 * ceil(npos / 2)))
    plt.clf()
    axs_py = fig.subplots(ceil(Int, npos / 2), 2)
    axs = [axs_py[i, j] for j in 1:size(axs_py, 2), i in 1:size(axs_py, 1)]

    for (axi, (ax, spanpos)) in enumerate(zip(axs, spanposs))
        target_y = spanpos * b / 2
        ax.set_title("2y/b=$(round(spanpos, digits=2))")

        for (ci, case) in enumerate(cases)
            body_case = case.body
            cps = body_case.controlpoints
            nrm = body_case.normals
            target_sign = sign(spanpos)
            on_requested_side(y) = target_sign == 0 || y == 0 || sign(y) == target_sign

            if snap_to_span_row
                side_idx = [p for p in 1:body_case.ncells if on_requested_side(cps[dim_span, p])]
                isempty(side_idx) && error("No panels found on requested side for spanwise $target_y for $(case.label).")

                yrows = unique(cps[dim_span, side_idx])
                row_y = yrows[argmin(abs.(yrows .- target_y))]
                rowtol = max(100eps(Float64), 1e-10 * b)
                idx = [p for p in side_idx if abs(cps[dim_span, p] - row_y) <= rowtol]
                isempty(idx) && error("No panels found within row tolerance $rowtol of snapped spanwise $row_y for $(case.label).")

                outside_slicetol = abs(row_y - target_y) > slicetol
                println("Cp slice snap: requested 2y/b=$(spanpos), y=$(target_y), ",
                        "snapped y=$(row_y), snapped 2y/b=$(2 * row_y / b), ",
                        "panels=$(length(idx))",
                        outside_slicetol ? " (outside slicetol=$(slicetol))" : "")
            else
                idx = [p for p in 1:body_case.ncells
                       if abs(cps[dim_span, p] - target_y) <= slicetol &&
                          on_requested_side(cps[dim_span, p])]
                isempty(idx) && error("No panels found within tolerance $slicetol of spanwise $target_y for $(case.label).")
            end

            points = cps[:, idx]
            xoc = points[1, :] ./ c
            upper = nrm[3, idx] .>= 0
            lower = .!upper
            ord_u = sortperm(xoc[upper])
            ord_l = sortperm(xoc[lower])

            lines = ax.plot(xoc[upper][ord_u], case.Cp[idx][upper][ord_u],
                            stls[min(ci, length(stls))];
                            label=case.label, linewidth=1.5, clip_on=false)
            color = lines[1].get_color()
            ax.plot(xoc[lower][ord_l], case.Cp[idx][lower][ord_l],
                    stls[min(ci, length(stls))];
                    label="_nolegend_", color=color, linewidth=1.5, clip_on=false)
        end

        ax.set_xlim(xlims)
        ax.set_ylim(ylims)
        ax.set_xticks(0:0.2:1)
        ax.set_yticks(ylims[1]:-0.5:ylims[2])
        if axi >= length(axs) - 1
            ax.set_xlabel("x/c")
        end
        if isodd(axi)
            ax.set_ylabel("Cp")
        end
        ax.spines["right"].set_visible(false)
        ax.spines["top"].set_visible(false)
    end

    place_figure_legend_right!(fig, axs)
    fig.subplots_adjust(left=0.11, right=0.74, bottom=0.11, top=0.95,
                        wspace=0.18, hspace=0.30)
    return fig, axs
end

function plot_spanwise_loading(loadings, labels;
                               stls=["-", "--"],
                               ylims=[0.0, 0.8])
    fig = plt.figure("simplewing_quadgrad_spanwise_loading", figsize=(7, 4))
    plt.clf()
    ax = fig.add_subplot(111)
    for (li, loading) in enumerate(loadings)
        ax.plot(loading.y2b, loading.cl, stls[min(li, length(stls))];
                label=labels[li], linewidth=1.5, clip_on=false)
    end
    ax.set_xlabel("2y/b")
    ax.set_ylabel("c_l")
    ax.set_xlim([-1, 1])
    ax.set_ylim(ylims)
    ax.set_xticks(-1:0.25:1)
    ax.legend(loc="best", fontsize=10, frameon=false)
    ax.spines["right"].set_visible(false)
    ax.spines["top"].set_visible(false)
    return fig, ax
end

if save_plots
    spanposs_cps = [0.5, 0.8, 0.9, 1.0]
    fig1, _ = plot_Cp_overlays(Tuple(cases), spanposs_cps, b, c;
                               stls=["-", "--", ":"],
                               snap_to_span_row=snap_cp_slices)
    fig1.savefig(joinpath(@__DIR__, "..", "simplewing_quadgrad_chordwise_loading.png"),
                 dpi=150, bbox_inches="tight")
    println("Saved chordwise Cp comparison plot to simplewing_quadgrad_chordwise_loading.png")

    fig2, _ = plot_spanwise_loading(
        Tuple(loadings),
        ["$(case.label) (CL=$(round(case.CL, digits=3)))" for case in cases];
        stls=["-", "--", ":"])
    fig2.tight_layout()
    fig2.savefig(joinpath(@__DIR__, "..", "simplewing_quadgrad_spanwise_loading.png"),
                 dpi=150)
    println("Saved spanwise loading comparison plot to simplewing_quadgrad_spanwise_loading.png")
end
