#=##############################################################################
# DESCRIPTION
#   Straight wing at an angle of attack of 8.4 deg. This example builds a
#   mirrored full-span lifting body, solves it with `steady!`, and compares the
#   standard +y-half-mirrored geometry against a -y-half-mirrored discretization.
#
# AUTHORSHIP
#   * Author    : Eduardo J. Alvarez
#   * Email     : Edo.AlvarezR@gmail.com
#   * Created   : Dec 2022
#   * License   : MIT License
=###############################################################################

import FLOWPanel as pnl
include(joinpath(pnl.examples_path, "helper_functions.jl"))
using GeometricTools: PyPlot as plt
import LinearAlgebra

const _norm = LinearAlgebra.norm
const _dot = LinearAlgebra.dot
const _cross = LinearAlgebra.cross

envflag(name, default=false) = lowercase(get(ENV, name, string(default))) in ("1", "true", "yes", "on")

airfoil_path    = joinpath(pnl.examples_path, "data") # Where to find airfoil contours

paraview        = true                          # Whether to write VTK output
snap_cp_slices  = envflag("FLOWPANEL_SIMPLEWING_SNAP_CP", true) # Whether to snap Cp slices to control-point rows

# ----------------- SIMULATION PARAMETERS --------------------------------------
AOA             = 8.4                           # (deg) angle of attack
aoa_tag         = replace(string(AOA), "." => "p")
magVinf         = 30.0                          # (m/s) freestream velocity
Vinf            = magVinf*[cos(AOA*pi/180), 0, sin(AOA*pi/180)] # Freestream

rho             = 1.225                         # (kg/m^3) air density

# ----------------- GEOMETRY DESCRIPTION ---------------------------------------
b               = 98*0.0254                     # (m) span length
ar              = 5.0                           # Aspect ratio b/c_tip
c               = b / ar                        # (m) root chord length
tr              = 1.0                           # Taper ratio c_tip/c_root
twist_root      = 0                             # (deg) twist at root
twist_tip       = 0                             # (deg) twist at tip
lambda          = 0                             # (deg) sweep
gamma           = 0                             # (deg) dihedral
airfoil         = "airfoil-rae101.csv"          # Airfoil contour file

# ----- Chordwise discretization
# n_rfl           = 8                             # Control number of chordwise panels
n_rfl         = 16                            # <-- uncomment this for finer discretization

NDIVS_rfl = [ (0.25, n_rfl,   10.0, false),
              (0.50, n_rfl,    1.0, true),
              (0.25, n_rfl, 1/10.0, false)]

# ----- Spanwise discretization
# n_span          = 15                            # Number of spanwise panels on each side of the wing
n_span        = 60                            # <-- uncomment this for finer discretization
NDIVS_span      = [(1.0, n_span, 1.0, false)]
n_loading_bins  = 2*n_span


# ----------------- GENERATE BODY ----------------------------------------------
println("Generating wing...")

bodyoptargs = (;
    kerneloffset=1e-6,
    kernelcutoff=1e-12,
    semiinfinite_wake=true,
    DBC=false,
)

kernel = pnl.VortexRing

run_name = if kernel <: pnl.ConstantDoublet
    "simplewing_doublet"
elseif kernel <: pnl.VortexRing
    "simplewing_vortexring"
else
    error("Neumann RigidWakeBody supports only ConstantDoublet or VortexRing kernels; got $(kernel)")
end
save_path = joinpath("data", run_name)

bodytype = pnl.RigidWakeBody{kernel}

function simplewing_mirrored_from_negative(b, ar, tr, twist_root, twist_tip, lambda, gamma;
                                           bodytype, bodyoptargs=(;),
                                           airfoil_root, airfoil_tip, airfoil_path,
                                           rfl_NDIVS, span_NDIVS, delim=",",
                                           mirror_tol=100eps(Float64),
                                           reference_nodes=nothing)
    half = simplewing(b, ar, tr, twist_root, twist_tip, lambda, gamma;
                      bodytype=bodytype, bodyoptargs=bodyoptargs,
                      airfoil_root=airfoil_root, airfoil_tip=airfoil_tip,
                      airfoil_path=airfoil_path,
                      rfl_NDIVS=rfl_NDIVS,
                      delim=delim,
                      span_NDIVS=span_NDIVS,
                      b_low=-1.0, b_up=0.0,
                      verify_spline=false,
                      verify_rflspline=false)

    half_nodes = half.nodes
    half_cells = half.cells
    mirror_index = Vector{Int}(undef, size(half_nodes, 2))
    nodes = copy(half_nodes)
    for ni in axes(half_nodes, 2)
        if abs(half_nodes[2, ni]) <= mirror_tol
            mirror_index[ni] = ni
        else
            nodes = hcat(nodes, [half_nodes[1, ni], -half_nodes[2, ni], half_nodes[3, ni]])
            mirror_index[ni] = size(nodes, 2)
        end
    end

    half_centers_y = [sum(half_nodes[2, half_cells[:, ci]]) / 3 for ci in axes(half_cells, 2)]
    neg_order = sort(collect(axes(half_cells, 2)); by=ci -> half_centers_y[ci])
    pos_order = sort(collect(axes(half_cells, 2)); by=ci -> -half_centers_y[ci])

    cells = Matrix{Int}(undef, 3, 2 * size(half_cells, 2))
    out_ci = 0
    for ci in neg_order
        out_ci += 1
        cells[:, out_ci] .= half_cells[:, ci]
    end
    for ci in pos_order
        out_ci += 1
        cells[:, out_ci] .= reverse(mirror_index[half_cells[:, ci]])
    end

    te_nodes = Int[]
    for col in eachcol(half.shedding[1])
        pi, nia, nib = col[1], col[2], col[3]
        push!(te_nodes, half_cells[nia, pi])
        push!(te_nodes, half_cells[nib, pi])
    end
    full_te_nodes = unique(vcat(te_nodes, mirror_index[te_nodes]))

    if !isnothing(reference_nodes)
        size(reference_nodes) == size(nodes) ||
            error("Cannot reindex negative-half mirror: reference node size $(size(reference_nodes)) differs from generated node size $(size(nodes)).")

        old_to_reference = zeros(Int, size(nodes, 2))
        reference_used = falses(size(reference_nodes, 2))
        for old_i in axes(nodes, 2)
            match_i = 0
            for ref_i in axes(reference_nodes, 2)
                reference_used[ref_i] && continue
                if maximum(abs.(view(nodes, :, old_i) .- view(reference_nodes, :, ref_i))) <= mirror_tol
                    match_i = ref_i
                    break
                end
            end
            match_i != 0 ||
                error("Cannot reindex negative-half mirror: generated node $old_i has no matching reference node within $(mirror_tol).")
            old_to_reference[old_i] = match_i
            reference_used[match_i] = true
        end

        all(reference_used) ||
            error("Cannot reindex negative-half mirror: not all reference nodes were matched.")

        nodes = copy(reference_nodes)
        cells = old_to_reference[cells]
        full_te_nodes = unique(old_to_reference[full_te_nodes])
    end

    sort!(full_te_nodes; by=ni -> nodes[2, ni])
    shedding = pnl.calc_shedding(nodes, cells, full_te_nodes, zeros(eltype(nodes), 3, 0))

    watertight, _ = pnl.iswatertight(nodes, cells)
    final_bodyoptargs = merge((ensure_winding=false,), bodyoptargs)
    return bodytype(nodes, cells, [shedding]; watertight, final_bodyoptargs...)
end

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


# ----------------- CALL SOLVER AND MONITORS -----------------------------------
println("Solving body (combined Bernoulli + Laplace monitor stack)...")

backend = pnl.FastMultipoleBackend()
Dhat = Vinf/_norm(Vinf)
Shat = [0, 1, 0]
Lhat = _cross(Dhat, Shat)

Sref = b^2 / ar
c_ref = b / ar
normalization = pnl.WingNormalization(rho, Sref, c_ref)
sectional_qc = 0.5*rho*magVinf^2*c_ref

frames = pnl.ReferenceFrame(body)

pressure_bernoulli = pnl.PressureBernoulli(rho)
force_bernoulli = pnl.ForceMonitor(1, 1; i_frame=-1, normalization=normalization,
    correct_kuttacondition=false, verbose=false)
spanwise_bernoulli = pnl.SpanwiseLoadingMonitor(n_loading_bins, 1;
    components=(lift=Lhat, drag=Dhat),
    span_axis=Shat,
    per_length=true,
    normalization=pnl.NoSectionalNormalization())

pressure_laplace = pnl.PressureLaplace((body,), rho;
    reference_panel=1, reference_pressure=0.0, verbose=false,
    unsteady=false,
    gradient_mode=:surface_velocity,
    acceleration_form=:lamb_vector)
force_laplace = pnl.ForceMonitor(1, 1; i_frame=-1, normalization=normalization,
    correct_kuttacondition=false, verbose=false)
spanwise_laplace = pnl.SpanwiseLoadingMonitor(n_loading_bins, 1;
    components=(lift=Lhat, drag=Dhat),
    span_axis=Shat,
    per_length=true,
    normalization=pnl.NoSectionalNormalization())

monitors = (pressure_bernoulli, force_bernoulli, spanwise_bernoulli,
            pressure_laplace,   force_laplace,   spanwise_laplace)

@time pnl.steady!(body, frames, Vinf;
    body_solvers=pnl.Backslash(body),
    backend=backend,
    monitors=monitors,
    path=paraview ? save_path : nothing,
    name=run_name*"_bernoulli_AOA$(aoa_tag)",
    verbose=false)

F_bernoulli = force_bernoulli.force[:, 1]
CL_bernoulli = _dot(F_bernoulli, Lhat)
CD_bernoulli = _dot(F_bernoulli, Dhat)

F_laplace = force_laplace.force[:, 1]
CL_laplace = _dot(F_laplace, Lhat)
CD_laplace = _dot(F_laplace, Dhat)

function solve_bernoulli_loading!(body_case, label)
    for i in eachindex(body_case.Das)
        body_case.Das[i] .= repeat(wake_direction, 1, size(body_case.Das[i], 2))
    end

    frames_case = pnl.ReferenceFrame(body_case)
    pressure_case = pnl.PressureBernoulli(rho)
    force_case = pnl.ForceMonitor(1, 1; i_frame=-1, normalization=normalization,
        correct_kuttacondition=false, verbose=false)
    spanwise_case = pnl.SpanwiseLoadingMonitor(n_loading_bins, 1;
        components=(lift=Lhat, drag=Dhat),
        span_axis=Shat,
        per_length=true,
        normalization=pnl.NoSectionalNormalization())
    monitors_case = (pressure_case, force_case, spanwise_case)

    println("Solving body ($(label) Bernoulli loading stack)...")
    @time pnl.steady!(body_case, frames_case, Vinf;
        body_solvers=pnl.Backslash(body_case),
        backend=backend,
        monitors=monitors_case,
        path=nothing,
        verbose=false)

    F_case = force_case.force[:, 1]
    return (; pressure=pressure_case, force=force_case, spanwise=spanwise_case,
            CL=_dot(F_case, Lhat), CD=_dot(F_case, Dhat))
end

body_negative_mirror = simplewing_mirrored_from_negative(
    b, ar, tr, twist_root, twist_tip, lambda, gamma;
    bodytype=bodytype, bodyoptargs=bodyoptargs,
    airfoil_root=airfoil, airfoil_tip=airfoil,
    airfoil_path=airfoil_path,
    rfl_NDIVS=NDIVS_rfl,
    delim=",",
    span_NDIVS=NDIVS_span,
    reference_nodes=body.nodes)
negative_mirror_bernoulli = solve_bernoulli_loading!(
    body_negative_mirror, "negative-half mirror")


# ----------------- DIAGNOSTICS ------------------------------------------------
println("\n#===== INTEGRATED FORCE COEFFICIENTS =====#")
@show CL_bernoulli CD_bernoulli
@show CL_laplace CD_laplace
@show negative_mirror_bernoulli.CL negative_mirror_bernoulli.CD

println("\n#===== MIRRORED DISCRETIZATION COMPARISON =====#")
println("+y half mirrored CL: $(CL_bernoulli)")
println("-y half mirrored CL: $(negative_mirror_bernoulli.CL)")
println("Absolute CL difference: $(abs(CL_bernoulli - negative_mirror_bernoulli.CL))")

normals = body.normals
Udotn = sum(body.velocity .* normals, dims=1)
resid = maximum(abs.(Udotn))
println("Max flow tangency residual (+y half mirrored): $resid")


# ----------------- VISUALIZATION ----------------------------------------------
if paraview
    println("Wrote primary steady! VTK files to $(save_path)/")
    pnl.write_vtk(joinpath(save_path, run_name * "_negative_half_mirror_bernoulli_AOA$(aoa_tag)"),
                  body_negative_mirror, 0, 0.0;
                  monitors=(negative_mirror_bernoulli.pressure,), i_system=1,
                  overwrite=true)
    println("Saved VTK for negative-half-mirror case to $(save_path)/")
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

function plot_Cp_overlays(cases, spanposs, b, c, rho, magVinf;
                          slicetol=0.03*b,
                          snap_to_span_row=false,
                          dim_span=2,
                          stls=["-", "--"],
                          xlims=[-0.1, 1.1],
                          ylims=[1.0, -1.5])
    qinf = 0.5 * rho * magVinf^2
    npos = length(spanposs)
    fig = plt.figure("simplewing_Cps_mirrored", figsize=(9.0, 3.2*ceil(npos/2)))
    plt.clf()
    axs_py = fig.subplots(ceil(Int, npos/2), 2)
    axs = [axs_py[i, j] for j in 1:size(axs_py, 2), i in 1:size(axs_py, 1)]

    for (axi, (ax, spanpos)) in enumerate(zip(axs, spanposs))
        target_y = spanpos * b / 2
        ax.set_title("2y/b=$(round(spanpos, digits=2))")

        for (ci, case) in enumerate(cases)
            body_case = case.body
            pressure = case.pressure
            cps = body_case.controlpoints
            nrm = body_case.normals
            target_sign = sign(spanpos)
            on_requested_side(y) = target_sign == 0 || y == 0 || sign(y) == target_sign

            if snap_to_span_row
                side_idx = [p for p in 1:body_case.ncells if on_requested_side(cps[dim_span, p])]
                isempty(side_idx) && error("No panels found on requested side for spanwise $target_y for $(case.label).")

                yrows = unique(cps[dim_span, side_idx])
                row_y = yrows[argmin(abs.(yrows .- target_y))]
                rowtol = max(100eps(Float64), 1e-10*b)
                idx = [p for p in side_idx if abs(cps[dim_span, p] - row_y) <= rowtol]
                isempty(idx) && error("No panels found within row tolerance $rowtol of snapped spanwise $row_y for $(case.label).")

                outside_slicetol = abs(row_y - target_y) > slicetol
                println("Cp slice snap: requested 2y/b=$(spanpos), y=$(target_y), ",
                        "snapped y=$(row_y), snapped 2y/b=$(2*row_y/b), ",
                        "panels=$(length(idx))",
                        outside_slicetol ? " (outside slicetol=$(slicetol))" : "")
            else
                idx = [p for p in 1:body_case.ncells
                       if abs(cps[dim_span, p] - target_y) <= slicetol &&
                          on_requested_side(cps[dim_span, p])]
                isempty(idx) && error("No panels found within tolerance $slicetol of spanwise $target_y for $(case.label).")
            end

            points = cps[:, idx]
            Cps = pressure[idx] ./ qinf
            xoc = points[1, :] ./ c
            upper = nrm[3, idx] .>= 0
            lower = .!upper
            ord_u = sortperm(xoc[upper])
            ord_l = sortperm(xoc[lower])

            lines = ax.plot(xoc[upper][ord_u], Cps[upper][ord_u], stls[min(ci, length(stls))];
                            label=case.label, linewidth=1.5, clip_on=false)
            color = lines[1].get_color()
            ax.plot(xoc[lower][ord_l], Cps[lower][ord_l], stls[min(ci, length(stls))];
                    label="_nolegend_", color=color, linewidth=1.5, clip_on=false)
        end

        ax.set_xlim(xlims)
        ax.set_ylim(ylims)
        ax.set_xticks(0:0.2:1)
        ax.set_yticks(ylims[1]:-0.5:ylims[2])
        if axi >= length(axs)-1
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

function monitor_sectional_coefficients(spanwise_monitor, b, sectional_qc)
    return spanwise_monitor.bin_center .* 2 ./ b,
           spanwise_monitor.load_components ./ sectional_qc
end

function plot_monitor_loading(spanwise_monitors, labels, b, sectional_qc;
                              stls=["-", "--"],
                              ylims=[0.0, 0.8])
    lift_i = findfirst(==(:lift), spanwise_monitors[1].component_names)
    lift_i === nothing && error("SpanwiseLoadingMonitor has no :lift component.")

    fig = plt.figure("simplewing_loading_mirrored_discretizations", figsize=(7, 4))
    plt.clf()
    ax = fig.add_subplot(111)

    for (mi, monitor) in enumerate(spanwise_monitors)
        y2b, coeffs = monitor_sectional_coefficients(monitor, b, sectional_qc)
        ax.plot(y2b, coeffs[lift_i, :], stls[min(mi, length(stls))];
                label=labels[mi], linewidth=1.5, clip_on=false)
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

spanposs_cps = [0.5, 0.8, 0.9, 1.0]
fig1, _ = plot_Cp_overlays((
        (; body=body, pressure=pressure_bernoulli.pressure[1], label="+y half mirrored"),
        (; body=body_negative_mirror, pressure=negative_mirror_bernoulli.pressure.pressure[1], label="-y half mirrored"),
    ),
    spanposs_cps, b, c, rho, magVinf;
    snap_to_span_row=snap_cp_slices)
fig1.savefig(joinpath(@__DIR__, "..", "simplewing_Cps_mirrored.png"),
             dpi=150, bbox_inches="tight")
println("Saved Cp overlay plot to simplewing_Cps_mirrored.png")

fig2, _ = plot_monitor_loading(
    (spanwise_bernoulli, negative_mirror_bernoulli.spanwise),
    ("+y half mirrored (CL=$(round(CL_bernoulli, digits=3)))",
     "-y half mirrored (CL=$(round(negative_mirror_bernoulli.CL, digits=3)))"),
    b, sectional_qc)
fig2.tight_layout()
fig2.savefig(joinpath(@__DIR__, "..", "simplewing_loading_mirrored_discretizations.png"),
             dpi=150)
println("Saved mirrored spanwise lift plot to simplewing_loading_mirrored_discretizations.png")
