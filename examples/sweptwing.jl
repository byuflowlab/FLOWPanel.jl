#=##############################################################################
# DESCRIPTION
    45deg swept-back wing at an angle of attack of 4.2deg. This wing has an
    aspect ratio of 5.0, a RAE 101 airfoil section with 12% thickness, and no
    dihedral, twist, nor taper. This test case matches the experimental setup
    of Weber, J., and Brebner, G., "Low-Speed Tests on 45-deg Swept-Back Wings,
    Part I," Tech. rep., 1951.

# AUTHORSHIP
  * Author    : Eduardo J. Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Dec 2022
  * License   : MIT License
=###############################################################################

import FLOWPanel as pnl
include(joinpath(pnl.examples_path, "helper_functions.jl"))
using GeometricTools: PyPlot as plt
import LinearAlgebra
const _norm  = LinearAlgebra.norm
const _dot   = LinearAlgebra.dot
const _cross = LinearAlgebra.cross

envflag(name, default=false) = lowercase(get(ENV, name, string(default))) in ("1", "true", "yes", "on")

run_name        = "sweptwing000"                # Name of this run; grid tag appended below

save_path       = joinpath("data", run_name)    # Where to save outputs; grid tag appended below
airfoil_path    = joinpath(pnl.examples_path, "data") # Where to find airfoil contours

paraview        = envflag("FLOWPANEL_SWEPTWING_VTK", true) # Whether to write VTK
load_vtk        = envflag("FLOWPANEL_SWEPTWING_LOAD_VTK", false) # Whether to load saved VTK for Cp-only plotting
snap_cp_slices  = envflag("FLOWPANEL_SWEPTWING_SNAP_CP", false) # Whether to snap Cp slices to control-point rows

# ----------------- SIMULATION PARAMETERS --------------------------------------
AOA             = 4.2                           # (deg) angle of attack
aoa_tag         = replace(string(AOA), "." => "p")
magVinf         = 30.0                          # (m/s) freestream velocity
Vinf            = magVinf*[cos(AOA*pi/180), 0, sin(AOA*pi/180)] # Freestream

rho             = 1.225                         # (kg/m^3) air density

# ----------------- GEOMETRY DESCRIPTION ---------------------------------------
b               = 98*0.0254                     # (m) span length
ar              = 5.0                           # Aspect ratio b/c_tip
tr              = 1.0                           # Taper ratio c_tip/c_root
twist_root      = 0                             # (deg) twist at root
twist_tip       = 0                             # (deg) twist at tip
lambda          = 45                            # (deg) sweep
gamma           = 0                             # (deg) dihedral
airfoil         = "airfoil-rae101.csv"          # Airfoil contour file

# ----- Chordwise discretization
n_rfl           = 24                            # Control number of chordwise panels
NDIVS_rfl = [ (0.25, n_rfl,   10.0, false),
              (0.50, n_rfl,    1.0, true),
              (0.25, n_rfl, 1/10.0, false)]

# ----- Spanwise discretization (full span, single loft)
# Uniform distribution: with `central=true, expansion=20` the *root* panels
# end up coarsest and the tips finest, which makes the inner Cp slice
# (2y/b ≈ 0.04) too under-resolved to plot cleanly. Uniform spacing keeps
# panel size constant across the span.
n_span_full     = 24                            # Number of spanwise panels across full span
NDIVS_span      = [(1.0, n_span_full, 1.0, true)]

grid_tag        = "nrf$(n_rfl)_nspan$(n_span_full)"
run_name        = "sweptwing000_" * grid_tag
save_path       = joinpath("data", run_name)
load_run_name   = get(ENV, "FLOWPANEL_SWEPTWING_LOAD_RUN_NAME", run_name)
load_path       = get(ENV, "FLOWPANEL_SWEPTWING_LOAD_PATH", joinpath("data", load_run_name))

# ----------------- GENERATE BODY ----------------------------------------------
println("Generating body...")

bodytype = pnl.RigidWakeBody{pnl.VortexRing, 1, Float64, false}

bodyoptargs = (;)

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
                      b_low=-1.0, b_up=0.0)

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
                                )
wake_direction = reshape(Vinf ./ magVinf, :, 1)
for i in eachindex(body.Das)
    body.Das[i] .= repeat(wake_direction, 1, size(body.Das[i], 2))
end

# Freestream at every control point
Uinfs = repeat(Vinf, 1, body.ncells)

println("Number of panels:\t$(body.ncells)")

function saved_body_vtu_path(path, body_name, idx=0)
    return joinpath(path, body_name, "$(body_name).$(idx).vtu")
end

function load_vtu_cell_scalar(path, body_name, field_name, expected_ncells; idx=0)
    vtu_path = saved_body_vtu_path(path, body_name, idx)
    isfile(vtu_path) ||
        error("Required body VTU file not found: $(vtu_path)")

    vtk = pnl.ReadVTK.VTKFile(vtu_path)
    cell_data = pnl.ReadVTK.get_cell_data(vtk)
    field_name in keys(cell_data) ||
        error("Required cell field '$(field_name)' missing from $(vtu_path). Available cell fields: $(collect(keys(cell_data))).")

    values = Vector{Float64}(pnl.ReadVTK.get_data(cell_data[field_name]))
    length(values) == expected_ncells ||
        error("Loaded field '$(field_name)' from $(vtu_path) has $(length(values)) values, expected $(expected_ncells).")

    return values
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

if load_vtk
    println("Loading saved VTK pressure for Cp-only plotting from $(load_path)/")
    pnl.calc_normals!(body)
    pnl.calc_controlpoints!(body)

    body_negative_mirror = simplewing_mirrored_from_negative(
        b, ar, tr, twist_root, twist_tip, lambda, gamma;
        bodytype=bodytype, bodyoptargs=bodyoptargs,
        airfoil_root=airfoil, airfoil_tip=airfoil,
        airfoil_path=airfoil_path,
        rfl_NDIVS=NDIVS_rfl,
        delim=",",
        span_NDIVS=NDIVS_span,
        reference_nodes=body.nodes)
    for i in eachindex(body_negative_mirror.Das)
        body_negative_mirror.Das[i] .= repeat(wake_direction, 1, size(body_negative_mirror.Das[i], 2))
    end
    pnl.calc_normals!(body_negative_mirror)
    pnl.calc_controlpoints!(body_negative_mirror)

    positive_body_name = load_run_name * "_bernoulli_AOA$(aoa_tag)_body1"
    negative_body_name = load_run_name * "_negative_half_mirror_bernoulli_AOA$(aoa_tag)"
    pressure_bernoulli_loaded = load_vtu_cell_scalar(load_path, positive_body_name,
                                                     "gauge pressure", body.ncells)
    pressure_negative_loaded = load_vtu_cell_scalar(load_path, negative_body_name,
                                                   "gauge pressure", body_negative_mirror.ncells)

    include(joinpath(pnl.examples_path, "sweptwing_postprocessing.jl"))

    if envflag("FLOWPANEL_SWEPTWING_PLOTS", true)
        side = 1
        spanposs_cps = side*parse.(Float64, keys(weber_Cps["$AOA"]))[[2, 4, 5, 7]]
        xLE_fn = y -> abs(y) * tan(lambda * pi / 180)
        fig1, axs = plot_Cps(body, pressure_bernoulli_loaded, spanposs_cps, b, rho, magVinf;
                                    xscaling=ar/b, AOA=AOA,
                                    xlims=[-0.1, 1.1], ylims=[1.0, -1.5], stl="-",
                                    slicetol=0.013*b, xLE_fn=xLE_fn,
                                    snap_to_span_row=snap_cp_slices,
                                    show_axis_legend=false,
                                    plot_vsp=false,
                                    plot_optargs=(label="+y half mirrored",))
        plot_Cps(body_negative_mirror, pressure_negative_loaded,
                 spanposs_cps, b, rho, magVinf;
                 _fig=fig1, _axs=axs,
                 xscaling=ar/b, AOA=AOA,
                 xlims=[-0.1, 1.1], ylims=[1.0, -1.5], stl="--",
                 slicetol=0.013*b, xLE_fn=xLE_fn,
                 snap_to_span_row=snap_cp_slices,
                 show_axis_legend=false,
                 plot_exp=false, plot_vsp=false,
                 plot_optargs=(label="-y half mirrored",))
        for ax in axs
            ax.set_xlim([-0.1, 1.1])
            ax.set_ylim([1.0, -1.5])
        end
        place_figure_legend_right!(fig1, axs)
        fig1.subplots_adjust(left=0.11, right=0.74, bottom=0.11, top=0.95,
                             wspace=0.18, hspace=0.30)
        fig1.savefig(joinpath(@__DIR__, "..", "sweptwing_Cps_mirrored.png"),
                     dpi=150, bbox_inches="tight")
        println("Saved Cp overlay plot to sweptwing_Cps_mirrored.png")
    end

    exit()
end


# ----------------- CALL SOLVER AND MONITORS -----------------------------------
println("Solving body (combined Bernoulli + Laplace monitor stack)...")

backend = pnl.FastMultipoleBackend()
Dhat = Vinf/_norm(Vinf)        # Drag direction
Shat = [0, 1, 0]               # Span direction
Lhat = _cross(Dhat, Shat)      # Lift direction

Sref = b^2 / ar
c_ref = b / ar
normalization = pnl.WingNormalization(rho, Sref, c_ref)
sectional_qc = 0.5*rho*magVinf^2*c_ref

frames = pnl.ReferenceFrame(body)

# Bernoulli post-processing stack
pressure_bernoulli = pnl.PressureBernoulli(rho)
force_bernoulli = pnl.ForceMonitor(1, 1; i_frame=-1, normalization=normalization,
    correct_kuttacondition=false, verbose=false)
spanwise_bernoulli = pnl.SpanwiseLoadingMonitor(n_span_full, 1;
    components=(lift=Lhat, drag=Dhat),
    span_axis=Shat,
    per_length=true,
    normalization=pnl.NoSectionalNormalization())

# Laplace post-processing stack on the *same* body. Both pressure monitors keep
# their own storage (`.pressure` vs `.p`) and VTK auto-disambiguates the field
# names, so the two stacks coexist in a single solve. The Bernoulli trio runs
# first so each force/spanwise monitor reads the pressure/force from its own
# immediately-preceding monitor in the tuple ordering below.
pressure_laplace = pnl.PressureLaplace((body,), rho;
    reference_panel=1, reference_pressure=0.0, verbose=false,
    unsteady=false,
    gradient_mode=:corrected_hessian,
    acceleration_form=:material_derivative)
force_laplace = pnl.ForceMonitor(1, 1; i_frame=-1, normalization=normalization,
    correct_kuttacondition=false, verbose=false)
spanwise_laplace = pnl.SpanwiseLoadingMonitor(n_span_full, 1;
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

function solve_bernoulli_loading!(body_case, label)
    frames_case = pnl.ReferenceFrame(body_case)
    pressure_case = pnl.PressureBernoulli(rho)
    force_case = pnl.ForceMonitor(1, 1; i_frame=-1, normalization=normalization,
        correct_kuttacondition=false, verbose=false)
    spanwise_case = pnl.SpanwiseLoadingMonitor(n_span_full, 1;
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

# `WingNormalization` already divides by 0.5 ρ |Vinf|² Sref, so force_laplace.force
# is in coefficient form. Project onto Lhat/Dhat to recover CL/CD.
F_lap = force_laplace.force[:, 1]
CL_laplace = _dot(F_lap, Lhat)
CD_laplace = _dot(F_lap, Dhat)


# ----------------- TRAILING-EDGE VELOCITY/PRESSURE DIAGNOSTICS ---------------
function reconstructed_half_jump(body)
    jump = zeros(3, body.ncells)
    pnl.compute_mu_gradient!(jump, body.controlpoints, body.normals,
        body.cells, body.neighbor,
        view(body.strength, :, pnl.get_Gammai(body)),
        view(body.shedding_full, 1:2, :);
        scale=0.5,
        nodes=body.nodes)
    return jump
end

function lower_te_panels(body)
    panels = Int[]
    for shedding in body.shedding
        for col in eachcol(shedding)
            pi, pj = col[1], col[4]
            if pj == -1
                push!(panels, pi)
            else
                lower = body.controlpoints[3, pi] <= body.controlpoints[3, pj] ? pi : pj
                push!(panels, lower)
            end
        end
    end
    panels = unique(panels)
    sort!(panels; by=p -> body.controlpoints[2, p])
    return panels
end

function adjacent_jump_summary(values, panels)
    length(panels) >= 2 || return (; n=0, mean=NaN, rms=NaN, max=NaN)
    jumps = Float64[]
    size(values, 1) == 1 && (values = reshape(vec(values), 1, :))
    for k in 1:(length(panels)-1)
        i, j = panels[k], panels[k+1]
        push!(jumps, _norm(values[:, j] .- values[:, i]))
    end
    return (; n=length(jumps),
             mean=sum(jumps) / length(jumps),
             rms=sqrt(sum(abs2, jumps) / length(jumps)),
             max=maximum(jumps))
end

function print_jump_summary(label, values, panels)
    s = adjacent_jump_summary(values, panels)
    println(rpad(label, 24), " n=$(s.n)",
            " mean=$(round(s.mean, sigdigits=5))",
            " rms=$(round(s.rms, sigdigits=5))",
            " max=$(round(s.max, sigdigits=5))")
end

function report_te_diagnostics(body, Vinf, p_bernoulli, p_laplace)
    jump = reconstructed_half_jump(body)
    pv_velocity = body.velocity .- jump
    pv_induced_velocity = pv_velocity .- repeat(Vinf, 1, body.ncells)
    panels = lower_te_panels(body)
    gamma = reshape(collect(view(body.strength, :, pnl.get_Gammai(body))), 1, :)
    bernoulli = reshape(collect(p_bernoulli), 1, :)
    laplace = reshape(collect(p_laplace), 1, :)

    println("\n#===== LOWER-TE ADJACENT JUMP DIAGNOSTICS =====#")
    @show length(panels)
    print_jump_summary("gamma", gamma, panels)
    print_jump_summary("PV induced velocity", pv_induced_velocity, panels)
    print_jump_summary("half-jump velocity", jump, panels)
    print_jump_summary("final velocity", body.velocity, panels)
    print_jump_summary("Bernoulli pressure", bernoulli, panels)
    print_jump_summary("Laplace pressure", laplace, panels)

    return (; panels, jump, pv_velocity, pv_induced_velocity)
end

te_diagnostics = report_te_diagnostics(body, Vinf, pressure_bernoulli.pressure[1],
                                       pressure_laplace.p[1])


# ----------------- VISUALIZATION ----------------------------------------------
if paraview
    println("Wrote VTK files to $(save_path)/")
end


# ----------------- COMPARISON TO EXPERIMENTAL DATA ----------------------------
include(joinpath(pnl.examples_path, "sweptwing_postprocessing.jl"))

save_outputs = false

fig_path = joinpath(pnl.examples_path, "..", "docs", "resources", "images")
outdata_path = joinpath(pnl.examples_path, "..", "docs", "resources", "data")

# `plot_Cps` uses the unstructured-mesh-compatible `slice_scalarfield` primitive
# and works on this body. `plot_deltaCps` / `plot_loading` still depend on
# structured-grid helpers (`slicefield` / `calcfield_sectionalforce`) and remain
# disabled until those are ported.
make_plots_cps = envflag("FLOWPANEL_SWEPTWING_PLOTS", true)
make_plots_loading = envflag("FLOWPANEL_SWEPTWING_PLOTS", true)

function monitor_sectional_coefficients(spanwise_monitor, b, sectional_qc)
    return spanwise_monitor.bin_center .* 2 ./ b,
           spanwise_monitor.load_components ./ sectional_qc
end

function report_spanwise_symmetry(label, spanwise_monitor, b, sectional_qc)
    lift_i = findfirst(==(:lift), spanwise_monitor.component_names)
    lift_i === nothing && error("SpanwiseLoadingMonitor has no :lift component.")

    _, coeffs = monitor_sectional_coefficients(spanwise_monitor, b, sectional_qc)
    counts_match = spanwise_monitor.counts == reverse(spanwise_monitor.counts)
    center_symmetry = maximum(abs.(spanwise_monitor.bin_center .+ reverse(spanwise_monitor.bin_center)))
    lift_symmetry = maximum(abs.(coeffs[lift_i, :] .- reverse(coeffs[lift_i, :])))

    println("$label mirrored bin counts: $counts_match")
    println("$label max mirrored bin-center sum: $(center_symmetry)")
    println("$label max mirrored sectional lift diff: $(lift_symmetry)")

    return (; counts_match, center_symmetry, lift_symmetry)
end

function cell_scalar_to_nodes(body, cell_values)
    length(cell_values) == body.ncells ||
        error("Expected one cell value per panel; got $(length(cell_values)) values for $(body.ncells) panels.")

    areas = if hasproperty(body, :areas) && length(getproperty(body, :areas)) == body.ncells
        collect(getproperty(body, :areas))
    else
        pnl.calc_areas(body)
    end

    values = zeros(promote_type(eltype(cell_values), eltype(areas)), size(body.nodes, 2))
    weights = zeros(eltype(areas), size(body.nodes, 2))
    for ci in 1:body.ncells
        area = areas[ci]
        value = cell_values[ci]
        for ni in view(body.cells, :, ci)
            values[ni] += area * value
            weights[ni] += area
        end
    end

    unused = findall(iszero, weights)
    isempty(unused) || error("Cannot interpolate cell values to nodes; $(length(unused)) nodes have no incident panels.")

    values ./= weights
    return values
end

function plot_monitor_loading(spanwise_monitors, labels, b, sectional_qc;
                              to_plot=collect(1:length(spanwise_monitors[1].component_names)),
                              stls=["-", "--"],
                              AOA=nothing,
                              xlims=[-1, 1],
                              ylims=([0.0, 0.8, 0.2], [-0.02, 0.08, 0.02], [-0.1, 0.1, 0.05]))
    fig = plt.figure(figsize=[7, 5*0.75]*2/3 .* [length(to_plot), 1])
    axs = fig.subplots(1, length(to_plot))
    axs = length(to_plot)==1 ? [axs] : [axs[i] for i in 1:length(to_plot)]

    for (axi, (ax, pi)) in enumerate(zip(axs, to_plot))
        if AOA != nothing && pi != 3
            vals_exp = (cls_web, cds_web[2:end, :])[pi]
            rowi = findfirst(a -> a==AOA, alphas_web)
            if rowi != nothing
                for f in [-1, 1]
                    ax.plot(f*y2b_web, vals_exp[rowi, :], "o--k",
                            label="Experimental"^(f==1), linewidth=0.5, clip_on=true)
                end
            else
                println("Experimental data at AOA=$(AOA) not found; valid AOAs are $(alphas_web).")
            end
        end

        for (mi, monitor) in enumerate(spanwise_monitors)
            y2b, coeffs = monitor_sectional_coefficients(monitor, b, sectional_qc)
            ax.plot(y2b, coeffs[pi, :], stls[min(mi, length(stls))];
                    label=labels[mi], linewidth=1.5, clip_on=false)
        end

        if xlims!=nothing; ax.set_xlim(xlims); end
        if ylims!=nothing; ax.set_ylim(ylims[axi][1:2]); end

        if xlims!=nothing; ax.set_xticks(xlims[1]:0.25:xlims[2]); end
        if ylims!=nothing; ax.set_yticks(ylims[axi][1]:ylims[axi][3]:ylims[axi][2]); end

        ax.set_xlabel(L"2y/b")
        coeff_symbols = (lift=L"c_l", drag=L"c_d", sideslip=L"c_s")
        ax.set_ylabel(get(coeff_symbols, spanwise_monitors[1].component_names[pi], L"c"))

        if axi==1
            ax.legend(loc="best", fontsize=10, frameon=false, reverse=true)
        end

        ax.spines["right"].set_visible(false)
        ax.spines["top"].set_visible(false)
    end

    return fig, axs
end

# --------- Summary (printed before plotting so we get numbers even if plots fail)
CLexp = CLs_web[2]
CDexp = CDs_web[2]

println("\n#===== INTEGRATED CL/CD =====#")
@show CL_bernoulli CL_laplace CLexp
@show CD_bernoulli CD_laplace CDexp
println("Bernoulli vs Laplace CL diff: $(round(abs(CL_laplace-CL_bernoulli), sigdigits=4))")
println("Bernoulli vs Laplace CD diff: $(round(abs(CD_laplace-CD_bernoulli), sigdigits=4))")
println("Bernoulli CL error: $(round(abs(CL_bernoulli-CLexp)/CLexp*100, digits=2))%")
println("Laplace   CL error: $(round(abs(CL_laplace-CLexp)/CLexp*100, digits=2))%")

println("\n#===== SPANWISE LOADING MONITORS =====#")
@show length(spanwise_bernoulli.bin_center) length(spanwise_laplace.bin_center) n_span_full
@show all(isfinite, spanwise_bernoulli.load_components)
@show all(isfinite, spanwise_laplace.load_components)

symmetry_bernoulli = report_spanwise_symmetry("Bernoulli", spanwise_bernoulli, b, sectional_qc)
symmetry_laplace = report_spanwise_symmetry("Laplace", spanwise_laplace, b, sectional_qc)
@assert symmetry_bernoulli.counts_match
@assert symmetry_laplace.counts_match
@assert symmetry_bernoulli.center_symmetry <= 1e-12
@assert symmetry_laplace.center_symmetry <= 1e-12
@assert symmetry_bernoulli.lift_symmetry <= 1e-8

# Build & solve the negative-half mirrored body up-front so it's available for
# both the Cp overlay (step 0 diagnostic) and the loading comparison.
body_negative_mirror = simplewing_mirrored_from_negative(
    b, ar, tr, twist_root, twist_tip, lambda, gamma;
    bodytype=bodytype, bodyoptargs=bodyoptargs,
    airfoil_root=airfoil, airfoil_tip=airfoil,
    airfoil_path=airfoil_path,
    rfl_NDIVS=NDIVS_rfl,
    delim=",",
    span_NDIVS=NDIVS_span,
    reference_nodes=body.nodes)
for i in eachindex(body_negative_mirror.Das)
    body_negative_mirror.Das[i] .= repeat(wake_direction, 1, size(body_negative_mirror.Das[i], 2))
end
negative_mirror_bernoulli = solve_bernoulli_loading!(
    body_negative_mirror, "negative-half mirror")

if paraview
    Gi = pnl.get_Gammai(body)
    gamma_pos = collect(view(body.strength, :, Gi))
    gamma_neg = collect(view(body_negative_mirror.strength, :, Gi))

    size(body.nodes) == size(body_negative_mirror.nodes) ||
        error("Cannot compute nodewise gamma difference: node array sizes differ ($(size(body.nodes)) vs $(size(body_negative_mirror.nodes))).")
    node_mismatch = maximum(abs.(body.nodes .- body_negative_mirror.nodes))
    println("Max node coordinate mismatch for gamma diagnostic: $(node_mismatch)")
    node_mismatch <= 1e-10 ||
        error("Cannot compute nodewise gamma difference: max node coordinate mismatch $(node_mismatch) exceeds 1e-10.")

    gamma_pos_nodes = cell_scalar_to_nodes(body, gamma_pos)
    gamma_neg_nodes = cell_scalar_to_nodes(body_negative_mirror, gamma_neg)
    gamma_diff_nodes = gamma_pos_nodes .- gamma_neg_nodes

    abs_gamma_diff_nodes = abs.(gamma_diff_nodes)
    println("Max abs(gamma_difference_node): $(maximum(abs_gamma_diff_nodes))")
    println("Mean abs(gamma_difference_node): $(sum(abs_gamma_diff_nodes) / length(abs_gamma_diff_nodes))")
    println("RMS abs(gamma_difference_node): $(sqrt(sum(abs2, gamma_diff_nodes) / length(gamma_diff_nodes)))")

    triangle_cells = [collect(c) .- 1 for c in eachcol(body.cells)]
    pnl._write_vtk_points_or_lines(joinpath(save_path, "gamma_difference_nodes"), body.nodes;
        cells=triangle_cells,
        point_data=(
            Dict("field_name" => "gamma_pos_node", "field_data" => gamma_pos_nodes),
            Dict("field_name" => "gamma_neg_node", "field_data" => gamma_neg_nodes),
            Dict("field_name" => "gamma_difference_node", "field_data" => gamma_diff_nodes),
        ),
        cell_data=(
            Dict("field_name" => "gamma_pos_cell", "field_data" => gamma_pos),
        ),
        override_cell_type=pnl.WriteVTK.VTKCellTypes.VTK_TRIANGLE)
    println("Saved nodewise gamma-difference VTU to $(save_path)/gamma_difference_nodes.vtu")
end

if make_plots_cps
    side = 1
    spanposs_cps = side*parse.(Float64, keys(weber_Cps["$AOA"]))[[2, 4, 5, 7]]
    # 45° sweep: LE x at spanwise position y is |y|*tan(λ).
    xLE_fn = y -> abs(y) * tan(lambda * pi / 180)
    fig1, axs = plot_Cps(body, pressure_bernoulli.pressure[1], spanposs_cps, b, rho, magVinf;
                                xscaling=ar/b, AOA=AOA,
                                xlims=[-0.1, 1.1], ylims=[1.0, -1.5], stl="-",
                                slicetol=0.013*b, xLE_fn=xLE_fn,
                                snap_to_span_row=snap_cp_slices,
                                show_axis_legend=false,
                                plot_vsp=false,
                                plot_optargs=(label="+y half mirrored",))
    plot_Cps(body_negative_mirror, negative_mirror_bernoulli.pressure.pressure[1],
             spanposs_cps, b, rho, magVinf;
             _fig=fig1, _axs=axs,
             xscaling=ar/b, AOA=AOA,
             xlims=[-0.1, 1.1], ylims=[1.0, -1.5], stl="--",
             slicetol=0.013*b, xLE_fn=xLE_fn,
             snap_to_span_row=snap_cp_slices,
             show_axis_legend=false,
             plot_exp=false, plot_vsp=false,
             plot_optargs=(label="-y half mirrored",))
    for ax in axs
        ax.set_xlim([-0.1, 1.1])
        ax.set_ylim([1.0, -1.5])
    end
    place_figure_legend_right!(fig1, axs)
    fig1.subplots_adjust(left=0.11, right=0.74, bottom=0.11, top=0.95,
                         wspace=0.18, hspace=0.30)
    fig1.savefig(joinpath(@__DIR__, "..", "sweptwing_Cps_mirrored.png"),
                 dpi=150, bbox_inches="tight")
    println("Saved Cp overlay plot to sweptwing_Cps_mirrored.png")
end

if make_plots_loading
    fig2, axs2 = plot_monitor_loading((spanwise_bernoulli, spanwise_laplace),
                                      ("Bernoulli", "Laplace"),
                                      b, sectional_qc; AOA=AOA)
    fig2.tight_layout()
    fig2.savefig(joinpath(@__DIR__, "..", "sweptwing_loading.png"), dpi=150)
    println("Saved spanwise loading plot to sweptwing_loading.png")
    if paraview
        pnl.write_vtk(joinpath(save_path,
                      run_name * "_negative_half_mirror_bernoulli_AOA$(aoa_tag)"),
                      body_negative_mirror, 0, 0.0;
                      monitors=(negative_mirror_bernoulli.pressure,), i_system=1,
                      overwrite=true)
        println("Saved VTK for negative-half-mirror case to $(save_path)/")
    end

    println("\n#===== MIRRORED DISCRETIZATION COMPARISON =====#")
    @show CL_bernoulli negative_mirror_bernoulli.CL CLexp

    fig3, axs3 = plot_monitor_loading(
        (spanwise_bernoulli, negative_mirror_bernoulli.spanwise),
        ("+y half mirrored (CL=$(round(CL_bernoulli, digits=3)))",
         "-y half mirrored (CL=$(round(negative_mirror_bernoulli.CL, digits=3)))"),
        b, sectional_qc;
        to_plot=[1],
        stls=["-", "--"],
        AOA=AOA,
        ylims=([0.0, 0.4, 0.1],))
    fig3.tight_layout()
    fig3.savefig(joinpath(@__DIR__, "..", "sweptwing_loading_mirrored_discretizations.png"), dpi=150)
    println("Saved mirrored spanwise lift plot to sweptwing_loading_mirrored_discretizations.png")
end
