## Bagai/Leishman-style finite free-wake relaxation probe for DJI9443 hover.
##
## This is a diagnostic harness, not a supported wake model. It reuses the
## DJI9443 setup from rotor_hover_force_method_audit.jl, fills a PanelWake with
## a finite helical sheet as the initial condition, then pseudo-implicitly
## relaxes the wake geometry and TE-jump strengths together.
##
## Quick smoke:
##   RUN_NAME=/private/tmp/bl_freewake_smoke SAVE_VTK=false MAX_ITER=2 WAKE_REVS=0.25 WAKE_ROWS_PER_REV=12 julia --project examples/rotor_hover_prescribed_helical_wake.jl
##   RUN_NAME=/private/tmp/bl_nested_smoke SAVE_VTK=false ITERATION_MODE=nested_pitch MAX_OUTER_ITER=2 MAX_INNER_ITER=3 WAKE_REVS=0.25 WAKE_ROWS_PER_REV=12 julia --project examples/rotor_hover_prescribed_helical_wake.jl
##   RUN_NAME=/private/tmp/bl_debug_smoke SAVE_VTK=true DEBUG_VTK_EVERY_ITER=true ITERATION_MODE=nested_pitch MAX_OUTER_ITER=1 MAX_INNER_ITER=2 WAKE_REVS=0.027777777777777776 WAKE_ROWS_PER_REV=36 julia --project examples/rotor_hover_prescribed_helical_wake.jl
##
## Intended diagnostic:
##   RUN_NAME=/private/tmp/bl_freewake_default SAVE_VTK=false MAX_ITER=20 julia --project examples/rotor_hover_prescribed_helical_wake.jl

import FLOWPanel as pnl
using FLOWPanel.FastMultipole.StaticArrays
import GeoIO
using LinearAlgebra: cross, dot, norm
using Printf: @sprintf

run_name  = get(ENV, "RUN_NAME", "rotor_hover_prescribed_helical_wake")
save_path = isabspath(run_name) ? run_name : joinpath("data", run_name)
isdir(save_path) || mkpath(save_path)

save_vtk          = get(ENV, "SAVE_VTK", "true") == "true"
debug_vtk_every_iter = get(ENV, "DEBUG_VTK_EVERY_ITER", "false") == "true"
wake_revs         = parse(Float64, get(ENV, "WAKE_REVS", "1.25"))
rows_per_rev      = parse(Int, get(ENV, "WAKE_ROWS_PER_REV", "36"))
max_iter          = parse(Int, get(ENV, "MAX_ITER", "20"))
iteration_mode    = Symbol(get(ENV, "ITERATION_MODE", "nested_pitch"))
inner_geometry_mode = Symbol(get(ENV, "INNER_GEOMETRY_MODE", "streamline"))
streamline_substeps = parse(Int, get(ENV, "STREAMLINE_SUBSTEPS", "1"))
max_outer_iter    = parse(Int, get(ENV, "MAX_OUTER_ITER", string(max_iter)))
max_inner_iter    = parse(Int, get(ENV, "MAX_INNER_ITER", "8"))
initial_only      = get(ENV, "INITIAL_ONLY", "false") == "true"
initial_strength_mode = Symbol(get(ENV, "INITIAL_STRENGTH_MODE", "candidate_constant"))
wake_relax        = parse(Float64, get(ENV, "WAKE_RELAX", "0.15"))
strength_relax    = parse(Float64, get(ENV, "STRENGTH_RELAX", "0.25"))
pseudo_dt_factor  = parse(Float64, get(ENV, "WAKE_PSEUDO_DT_FACTOR", "1.0"))
max_node_step_R   = parse(Float64, get(ENV, "MAX_NODE_STEP_R", "0.03"))
gamma_tol         = parse(Float64, get(ENV, "GAMMA_TOL", "5e-3"))
node_tol_R        = parse(Float64, get(ENV, "NODE_TOL_R", "2e-3"))
wake_core_size    = parse(Float64, get(ENV, "WAKE_CORE_SIZE", "1e-3"))
initial_inflow    = parse(Float64, get(ENV, "INITIAL_INFLOW", "0.08"))
axial_advance_ratio = parse(Float64, get(ENV, "AXIAL_ADVANCE_RATIO", "0.0"))
ct_tol            = parse(Float64, get(ENV, "CT_TOL", "2e-4"))
sample_r_min      = parse(Float64, get(ENV, "SAMPLE_R_MIN", "0.25"))
sample_r_max      = parse(Float64, get(ENV, "SAMPLE_R_MAX", "0.95"))
near_wake_mode    = Symbol(get(ENV, "NEAR_WAKE_MODE", "das_offset"))
fountain_min_upstream_rows = parse(Int, get(ENV, "FOUNTAIN_MIN_UPSTREAM_ROWS", "2"))
fountain_tol_R    = parse(Float64, get(ENV, "FOUNTAIN_TOL_R", "0.05"))

# ----- fixed operating point / geometry (DJI9443 40_40, matches audit) -------
msh_file = joinpath(pnl.examples_path, "data", "dji9443_new_40_40.msh")
isfile(msh_file) || error("Mesh file not found: $(msh_file)")
mesh_tag = splitext(basename(msh_file))[1]

# 0-based ParaView TE seed point IDs (see rotor_hover_pressure_comparison.jl)
te_indices_1 = [1614, 1574, 45]   .+ 1
te_indices_2 = [3324, 3284, 1755] .+ 1

rho     = 1.179
RPM     = parse(Float64, get(ENV, "RPM", "5400"))
R       = 0.119
shedding_r_over_R = 0.1

kerneloffset_panel   = R * 1e-10
kerneloffset_targets = parse(Float64, get(ENV, "KERNELOFFSET", "1e-3"))
kernelcutoff         = R * 1e-13
init_Das_eta_kinematic = 0.2
set_Das_min_kinematic_displacement = 0.01 * R

axial_dimension  = 1
radial_dimension = 2
omega_axis       = SVector{3}(-1.0, 0.0, 0.0)
default_helix_winding_sign = sign(omega_axis[axial_dimension])
default_axial_wake_sign = -sign(omega_axis[axial_dimension])
helix_winding_sign = parse(Float64, get(ENV, "HELIX_WINDING_SIGN",
    string(default_helix_winding_sign)))
axial_wake_sign = parse(Float64, get(ENV, "WAKE_AXIAL_SIGN",
    string(default_axial_wake_sign)))
omega            = 2 * pi * RPM / 60
dt               = 60 / RPM / 36
wake_row_dt      = 60 / RPM / rows_per_rev
pseudo_dt        = pseudo_dt_factor * wake_row_dt
tip_speed        = abs(omega) * R
# See docs/src/prescribed_helical_wake.md for the pitch-ratio derivation.
initial_helix_pitch_ratio = axial_advance_ratio + initial_inflow
Vinf             = [i == axial_dimension ? axial_wake_sign * axial_advance_ratio * tip_speed : 0.0
                    for i in 1:3]
Uinf(t)          = Vinf

msh = GeoIO.load(msh_file).geometry
base_nodes, base_cells = pnl.meshes2nodes_cells(msh)
base_nodes .*= R / maximum(base_nodes[radial_dimension, :])

kernel = Union{pnl.ConstantSource, pnl.VortexRing}
DBC    = kernel == pnl.VortexRing ? false : true

function make_shedding_bbox(nodes, seed_nodes)
    radial_midpoint = sum(nodes[radial_dimension, seed_nodes]) / length(seed_nodes)
    radial_sign = sign(radial_midpoint)
    lower = [minimum(nodes[i, :]) for i in 1:size(nodes, 1)]
    upper = [maximum(nodes[i, :]) for i in 1:size(nodes, 1)]
    padding = max(sqrt(eps(eltype(nodes))) * R, R * 1e-6)
    lower .-= padding
    upper .+= padding
    radial_cutoff = shedding_r_over_R * R
    radial_sign > 0 ? (lower[radial_dimension] = radial_cutoff - padding) :
                      (upper[radial_dimension] = -radial_cutoff + padding)
    return (SVector{3}(lower...), SVector{3}(upper...))
end

const BASE_ROTOR = pnl.RigidWakeBody{kernel}(base_nodes, base_cells, pnl.noshedding;
    kerneloffset=kerneloffset_panel, kerneloffset_panel, kerneloffset_targets,
    kernelcutoff, semiinfinite_wake=false, watertight=true, DBC)
const SHEDDING1 = pnl.calc_shedding_from_seed(BASE_ROTOR.nodes, BASE_ROTOR.cells,
    te_indices_1[1], te_indices_1[2];
    bbox=make_shedding_bbox(BASE_ROTOR.nodes, te_indices_1[1:2]),
    normal_jump_tol=0.2, max_turn_angle=pi/3, debug=false)
const SHEDDING2 = pnl.calc_shedding_from_seed(BASE_ROTOR.nodes, BASE_ROTOR.cells,
    te_indices_2[1], te_indices_2[2];
    bbox=make_shedding_bbox(BASE_ROTOR.nodes, te_indices_2[1:2]),
    normal_jump_tol=0.2, max_turn_angle=pi/3, debug=false)

function build_rotor_and_frames()
    rotor = pnl.RigidWakeBody{kernel}(copy(BASE_ROTOR.nodes), copy(BASE_ROTOR.cells),
        [copy(SHEDDING1), copy(SHEDDING2)];
        kerneloffset=kerneloffset_panel, kerneloffset_panel, kerneloffset_targets,
        kernelcutoff, semiinfinite_wake=false, watertight=true,
        ensure_winding=true, DBC)

    frames = pnl.ReferenceFrame(rotor;
        origin=SVector{3}(0.0, 0.0, 0.0), v=SVector{3}(0.0, 0.0, 0.0),
        ω_axis=omega_axis, ω=omega,
        R=SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
        name="vehicle", child_index=Int[], dependent_index=[1])

    pnl.initialize_Das!((rotor,), frames, Uinf, 0.0, dt;
        set_Das_eta_kinematic=init_Das_eta_kinematic,
        set_Das_min_kinematic_displacement)
    return rotor, frames
end

rotor, frames = build_rotor_and_frames()
n_wake_rows = max(1, ceil(Int, wake_revs * rows_per_rev))
wake = pnl.PanelWake(rotor; nwakerows=n_wake_rows, core_size=wake_core_size,
    unsteady_filament=false)
backend = pnl.FastMultipoleBackend(;
    expansion_order=parse(Int, get(ENV, "FMM_EXPANSION_ORDER", "8")),
    multipole_acceptance=parse(Float64, get(ENV, "FMM_ACCEPTANCE", "0.4")),
    leaf_size=parse(Int, get(ENV, "FMM_LEAF_SIZE", "20")),
)

rotation_x(theta) = SMatrix{3,3}(
    1.0, 0.0,         0.0,
    0.0, cos(theta), -sin(theta),
    0.0, sin(theta),  cos(theta),
)

function te_vertices(system, shedding)
    out = zeros(eltype(system.nodes), 3, size(shedding, 2) + 1)
    for i_shed in axes(shedding, 2)
        i_panel = shedding[1, i_shed]
        idx = shedding[3, i_shed] # nib
        out[:, i_shed] .= view(system.nodes, :, system.cells[idx, i_panel])
    end
    i_panel = shedding[1, end]
    idx = shedding[2, end] # final nia
    out[:, end] .= view(system.nodes, :, system.cells[idx, i_panel])
    return out
end

function validate_near_wake_mode(mode::Symbol)
    mode in (:te, :das_offset) ||
        error("Invalid NEAR_WAKE_MODE=$(mode). Expected te or das_offset.")
    return mode
end

function validate_iteration_mode(mode::Symbol)
    mode in (:coupled_relax, :nested_pitch) ||
        error("Invalid ITERATION_MODE=$(mode). Expected coupled_relax or nested_pitch.")
    return mode
end

function validate_inner_geometry_mode(mode::Symbol)
    mode in (:streamline, :length_preserving) ||
        error("Invalid INNER_GEOMETRY_MODE=$(mode). Expected streamline or length_preserving.")
    return mode
end

function validate_initial_strength_mode(mode::Symbol)
    mode in (:zero, :candidate_latest, :candidate_constant) ||
        error("Invalid INITIAL_STRENGTH_MODE=$(mode). Expected zero, candidate_latest, or candidate_constant.")
    return mode
end

function effective_near_wake_mode()
    validate_iteration_mode(iteration_mode)
    return iteration_mode == :nested_pitch ? :das_offset : validate_near_wake_mode(near_wake_mode)
end

function near_wake_vertices(system, i_surf; mode)
    validate_near_wake_mode(mode)
    out = te_vertices(system, system.shedding[i_surf])
    if mode == :das_offset
        out .+= system.Das[i_surf]
    end
    return out
end

function pin_near_wake_row!(wake, system; mode)
    for i_surf in eachindex(wake.nodes)
        wake.nodes[i_surf][:, 1, :] .= near_wake_vertices(system, i_surf; mode)
    end
    return nothing
end

function populate_helical_nodes!(wake, system; pitch_inflow_ratio, mode)
    vi = max(pitch_inflow_ratio, 0.0) * tip_speed
    for i_surf in eachindex(wake.nodes)
        seed = near_wake_vertices(system, i_surf; mode)
        nodes = wake.nodes[i_surf]
        fill!(nodes, zero(eltype(nodes)))
        for jrow in 1:(wake.nwakes[] + 1)
            age_revs = (jrow - 1) / rows_per_rev
            theta = 2 * pi * age_revs
            q = rotation_x(theta * helix_winding_sign)
            dx = axial_wake_sign * vi * theta / max(abs(omega), eps())
            for jcol in axes(nodes, 3)
                r = q * SVector{3}(seed[1, jcol], seed[2, jcol], seed[3, jcol])
                nodes[:, jrow, jcol] .= r .+ SVector{3}(dx, 0.0, 0.0)
            end
        end
    end
    return nothing
end

function candidate_strengths(system)
    out = Vector{Vector{Float64}}(undef, length(system.shedding))
    for i_surf in eachindex(system.shedding)
        out[i_surf] = [pnl._get_wakestrength_Gamma(system, ished, i_surf)
                       for ished in axes(system.shedding[i_surf], 2)]
    end
    return out
end

function strength_stats(values)
    total = 0.0
    maxval = 0.0
    n = 0
    for surf_values in values, gamma in surf_values
        agamma = abs(gamma)
        total += agamma
        maxval = max(maxval, agamma)
        n += 1
    end
    return (mean = n == 0 ? 0.0 : total / n, max = maxval)
end

function relax_strengths!(wake, candidates; relax)
    gamma_num = 0.0
    gamma_den = 0.0
    relaxed = Vector{Vector{Float64}}(undef, length(wake.strength))
    for i_surf in eachindex(wake.strength)
        str = wake.strength[i_surf]
        relaxed[i_surf] = similar(candidates[i_surf])
        for ished in axes(str, 3)
            current = str[1, 1, ished]
            candidate = candidates[i_surf][ished]
            next_gamma = (1 - relax) * current + relax * candidate
            relaxed[i_surf][ished] = next_gamma
            gamma_num += abs(next_gamma - current)
            gamma_den += abs(current)
            str[1, :, ished] .= next_gamma
            if size(str, 1) > 1
                str[2:end, :, ished] .= zero(eltype(str))
            end
        end
    end
    return relaxed, gamma_num / max(gamma_den, eps())
end

function advance_latest_strengths!(wake, candidates; relax)
    gamma_num = 0.0
    gamma_den = 0.0
    active_rows = wake.nwakes[]
    for i_surf in eachindex(wake.strength)
        str = wake.strength[i_surf]
        if active_rows >= 2
            for jrow in active_rows:-1:2
                str[:, jrow, :] .= view(str, :, jrow - 1, :)
            end
        end
        for ished in axes(str, 3)
            current = str[1, 1, ished]
            candidate = candidates[i_surf][ished]
            next_gamma = (1 - relax) * current + relax * candidate
            gamma_num += abs(next_gamma - current)
            gamma_den += abs(current)
            str[1, 1, ished] = next_gamma
        end
        if size(str, 1) > 1
            str[2:end, 1:active_rows, :] .= zero(eltype(str))
        end
        if active_rows + 1 <= size(str, 2)
            str[:, active_rows + 1:end, :] .= zero(eltype(str))
        end
    end
    return gamma_num / max(gamma_den, eps())
end

function pin_latest_strengths!(wake, candidates)
    gamma_num = 0.0
    gamma_den = 0.0
    active_rows = wake.nwakes[]
    for i_surf in eachindex(wake.strength)
        str = wake.strength[i_surf]
        if active_rows >= 2
            for jrow in active_rows:-1:2
                str[:, jrow, :] .= view(str, :, jrow - 1, :)
            end
        end
        for ished in axes(str, 3)
            current = str[1, 1, ished]
            candidate = candidates[i_surf][ished]
            gamma_num += abs(candidate - current)
            gamma_den += abs(current)
            str[1, 1, ished] = candidate
        end
        if size(str, 1) > 1
            str[2:end, 1:active_rows, :] .= zero(eltype(str))
        end
        if active_rows + 1 <= size(str, 2)
            str[:, active_rows + 1:end, :] .= zero(eltype(str))
        end
    end
    return gamma_den == 0 ? (gamma_num == 0 ? 0.0 : Inf) : gamma_num / gamma_den
end

snapshot_strengths(wake) = [copy(str) for str in wake.strength]

function max_strength_delta(wake, snapshot)
    delta = 0.0
    for i in eachindex(wake.strength)
        delta = max(delta, maximum(abs.(wake.strength[i] .- snapshot[i])))
    end
    return delta
end

function assert_strength_unchanged(wake, snapshot; context="")
    delta = max_strength_delta(wake, snapshot)
    if delta > 64 * eps(Float64)
        error("Wake strength mutated during pinned-strength inner loop$(isempty(context) ? "" : " ($(context))"): max abs delta $(delta)")
    end
    return delta
end

snapshot_body_strength(system) = copy(system.strength)

function max_body_strength_delta(system, snapshot)
    return maximum(abs.(system.strength .- snapshot))
end

function assert_body_strength_unchanged(system, snapshot; context="")
    delta = max_body_strength_delta(system, snapshot)
    if delta > 64 * eps(Float64)
        error("Body strength mutated during pinned-strength inner loop$(isempty(context) ? "" : " ($(context))"): max abs delta $(delta)")
    end
    return delta
end

rowget(r, name::Symbol, default) = name in propertynames(r) ? getproperty(r, name) : default

function max_row_value(rows, name::Symbol)
    vals = [rowget(r, name, NaN) for r in rows]
    vals = filter(isfinite, vals)
    return isempty(vals) ? NaN : maximum(vals)
end

function classify_inner_failure(rows; tol_R)
    isempty(rows) && return "no_inner_iterations"
    last = rows[end]
    if rowget(last, :capped_node_fraction, 0.0) > 0.5
        return "capped-step limited"
    elseif :fountain in propertynames(last) && last.fountain.rows_with_upstream_motion > 0
        return "fountain/upstream directed"
    end
    residuals = [rowget(r, :max_target_residual_R, NaN) for r in rows]
    residuals = filter(isfinite, residuals)
    length(residuals) < 4 && return "residual floor"
    tail = residuals[max(1, end - min(9, length(residuals) - 1)):end]
    diffs = diff(tail)
    sign_changes = sum(sign(diffs[i]) != sign(diffs[i - 1]) for i in 2:length(diffs)
        if abs(diffs[i]) > 0.05 * tol_R && abs(diffs[i - 1]) > 0.05 * tol_R)
    if sign_changes >= 3
        return "oscillatory"
    elseif maximum(tail) - minimum(tail) < 0.1 * max(sum(tail) / length(tail), tol_R)
        return "residual floor"
    else
        return "not settled"
    end
end

function wake_strength_stats(wake)
    total = 0.0
    maxval = 0.0
    n = 0
    active_rows = wake.nwakes[]
    for str in wake.strength, jrow in 1:active_rows, ished in axes(str, 3)
        agamma = abs(str[1, jrow, ished])
        total += agamma
        maxval = max(maxval, agamma)
        n += 1
    end
    return (mean = n == 0 ? 0.0 : total / n, max = maxval)
end

function set_wake_strengths_from_candidates!(wake, candidates; mode)
    validate_initial_strength_mode(mode)
    active_rows = wake.nwakes[]
    fill!.(wake.strength, 0.0)
    mode == :zero && return wake_strength_stats(wake)

    for i_surf in eachindex(wake.strength)
        str = wake.strength[i_surf]
        for ished in axes(str, 3)
            gamma = candidates[i_surf][ished]
            if mode == :candidate_latest
                str[1, 1, ished] = gamma
            elseif mode == :candidate_constant
                str[1, 1:active_rows, ished] .= gamma
            end
        end
        if size(str, 1) > 1
            str[2:end, :, :] .= zero(eltype(str))
        end
        if active_rows + 1 <= size(str, 2)
            str[:, active_rows + 1:end, :] .= zero(eltype(str))
        end
    end
    return wake_strength_stats(wake)
end

function wake_row_inflow_ratio(wake)
    total = 0.0
    n = 0
    for (isurf, nodes) in enumerate(wake.nodes)
        vel = wake.velocity[isurf]
        for jrow in 2:(wake.nwakes[] + 1), jcol in axes(nodes, 3)
            y = nodes[2, jrow, jcol]
            z = nodes[3, jrow, jcol]
            rr = sqrt(y^2 + z^2) / R
            if sample_r_min <= rr <= sample_r_max
                induced_axial = abs(vel[axial_dimension, jrow, jcol] - Vinf[axial_dimension])
                total += induced_axial
                n += 1
            end
        end
    end
    return n == 0 ? NaN : max(total / n / tip_speed, 0.0)
end

downstream_component(v) = axial_wake_sign * v[axial_dimension]
signed_anchor_distance(nodes, jrow, jcol, anchor_col) =
    axial_wake_sign * (nodes[axial_dimension, jrow, jcol] -
                       nodes[axial_dimension, 1, anchor_col])

function initial_wake_direction_stats(wake)
    total = 0.0
    minval = Inf
    maxval = -Inf
    n = 0
    upstream = 0
    for nodes in wake.nodes
        size(nodes, 2) >= 2 || continue
        for jcol in axes(nodes, 3)
            ds = signed_anchor_distance(nodes, 2, jcol, jcol)
            total += ds
            minval = min(minval, ds)
            maxval = max(maxval, ds)
            upstream += ds <= 0 ? 1 : 0
            n += 1
        end
    end
    return (
        mean_R = n == 0 ? NaN : total / n / R,
        min_R = isfinite(minval) ? minval / R : NaN,
        max_R = isfinite(maxval) ? maxval / R : NaN,
        upstream_nodes = upstream,
        n_nodes = n,
    )
end

function write_initial_vtk!(wake, system)
    tag = basename(run_name)
    pnl.write_vtk(joinpath(save_path, tag * "_initial_body1"), system, 0, 0.0;
        overwrite=true, compress=true)
    pnl.write_vtk(joinpath(save_path, tag * "_initial_wake"), wake, 0, 0.0;
        overwrite=true, compress=true)
    return nothing
end

function write_current_paraview_state!(system, wake; vtk_iteration::Int,
        vtk_series_suffix::AbstractString, overwrite::Bool=false, final::Bool=false)
    wake_snapshot = snapshot_strengths(wake)
    body_snapshot = snapshot_body_strength(system)
    if final
        body_name = joinpath(save_path, "$(basename(run_name))_body1")
        wake_name = joinpath(save_path, "$(basename(run_name))_wake")
        filament_name = joinpath(save_path, "$(basename(run_name))_wake_filaments")
        idx = 0
        t = 0.0
        vtk_overwrite = true
    else
        body_name = joinpath(save_path,
            "$(basename(run_name))_body1_$(vtk_series_suffix)")
        wake_name = joinpath(save_path,
            "$(basename(run_name))_wake_$(vtk_series_suffix)")
        filament_name = joinpath(save_path,
            "$(basename(run_name))_wake_filaments_$(vtk_series_suffix)")
        idx = vtk_iteration
        t = vtk_iteration
        vtk_overwrite = overwrite
    end

    pnl.write_vtk(body_name, system, idx, t; overwrite=vtk_overwrite, compress=true)
    pnl.write_vtk(wake_name, wake, idx, t;
        overwrite=vtk_overwrite, compress=true, filament_name=filament_name)
    assert_strength_unchanged(wake, wake_snapshot;
        context="current-state VTK write $(vtk_series_suffix) step $(vtk_iteration)")
    assert_body_strength_unchanged(system, body_snapshot;
        context="current-state VTK write $(vtk_series_suffix) step $(vtk_iteration)")
    return nothing
end

function fountain_flow_diagnostic(wake)
    node_ux_negative = 0
    panel_ux_negative = 0
    nodes_upstream_of_anchor = 0
    panels_crossing_anchor_plane = 0
    rows_with_upstream_motion = 0
    min_anchor_distance = Inf
    min_signed_velocity = Inf
    min_panel_signed_velocity = Inf
    min_raw_x = Inf
    min_raw_ux = Inf
    x_tol = max(fountain_tol_R * R, 1e-10 * R, sqrt(eps(Float64)) * R)

    for i_surf in eachindex(wake.nodes)
        nodes = wake.nodes[i_surf]
        vel = wake.velocity[i_surf]
        for jrow in 2:(wake.nwakes[] + 1)
            row_signed_velocity_total = 0.0
            row_n = 0
            for jcol in axes(nodes, 3)
                ds = signed_anchor_distance(nodes, jrow, jcol, jcol)
                ux = vel[axial_dimension, jrow, jcol]
                signed_velocity = downstream_component(view(vel, :, jrow, jcol))
                min_anchor_distance = min(min_anchor_distance, ds)
                min_signed_velocity = min(min_signed_velocity, signed_velocity)
                min_raw_x = min(min_raw_x, nodes[axial_dimension, jrow, jcol])
                min_raw_ux = min(min_raw_ux, ux)
                node_ux_negative += ux < 0 ? 1 : 0
                nodes_upstream_of_anchor += ds < -x_tol ? 1 : 0
                row_signed_velocity_total += signed_velocity
                row_n += 1
            end
            rows_with_upstream_motion += row_n > 0 &&
                row_signed_velocity_total / row_n < -x_tol ? 1 : 0
        end

        for jrow in 1:wake.nwakes[], jcol in axes(wake.strength[i_surf], 3)
            ds = (
                signed_anchor_distance(nodes, jrow, jcol, jcol),
                signed_anchor_distance(nodes, jrow + 1, jcol, jcol),
                signed_anchor_distance(nodes, jrow + 1, jcol + 1, jcol + 1),
                signed_anchor_distance(nodes, jrow, jcol + 1, jcol + 1),
            )
            panels_crossing_anchor_plane += minimum(ds) < -x_tol && maximum(ds) > x_tol ? 1 : 0
            ux = 0.25 * (
                vel[axial_dimension, jrow, jcol] +
                vel[axial_dimension, jrow + 1, jcol] +
                vel[axial_dimension, jrow + 1, jcol + 1] +
                vel[axial_dimension, jrow, jcol + 1]
            )
            signed_velocity = axial_wake_sign * ux
            min_panel_signed_velocity = min(min_panel_signed_velocity, signed_velocity)
            panel_ux_negative += ux < 0 ? 1 : 0
        end
    end

    fountain = nodes_upstream_of_anchor > 0 || panels_crossing_anchor_plane > 0
    return (
        fountain=fountain,
        node_ux_negative=node_ux_negative,
        panel_ux_negative=panel_ux_negative,
        nodes_upstream_of_anchor=nodes_upstream_of_anchor,
        panels_crossing_anchor_plane=panels_crossing_anchor_plane,
        rows_with_upstream_motion=rows_with_upstream_motion,
        min_anchor_distance_R=isfinite(min_anchor_distance) ? min_anchor_distance / R : NaN,
        min_signed_velocity_tip=isfinite(min_signed_velocity) ? min_signed_velocity / tip_speed : NaN,
        min_panel_signed_velocity_tip=isfinite(min_panel_signed_velocity) ? min_panel_signed_velocity / tip_speed : NaN,
        min_raw_x_R=isfinite(min_raw_x) ? min_raw_x / R : NaN,
        min_raw_ux_tip=isfinite(min_raw_ux) ? min_raw_ux / tip_speed : NaN,
    )
end

function empty_node_relaxation_stats()
    return (
        residual_mean_R=NaN,
        residual_max_R=NaN,
        capped_node_fraction=NaN,
        applied_mean_R=NaN,
        applied_max_R=NaN,
        segment_length_error_mean_R=NaN,
        segment_length_error_max_R=NaN,
        row_stats=NamedTuple[],
        mean_R=NaN,
        max_R=NaN,
    )
end

function snapshot_streamwise_segment_lengths(wake)
    lengths = Vector{Matrix{Float64}}(undef, length(wake.nodes))
    for i_surf in eachindex(wake.nodes)
        nodes = wake.nodes[i_surf]
        lengths[i_surf] = zeros(Float64, wake.nwakes[], size(nodes, 3))
        for jseg in 1:wake.nwakes[], jcol in axes(nodes, 3)
            dx = SVector{3}(
                nodes[1, jseg + 1, jcol] - nodes[1, jseg, jcol],
                nodes[2, jseg + 1, jcol] - nodes[2, jseg, jcol],
                nodes[3, jseg + 1, jcol] - nodes[3, jseg, jcol],
            )
            lengths[i_surf][jseg, jcol] = norm(dx)
        end
    end
    return lengths
end

function segment_length_error_stats(wake, target_lengths; start_row::Int=2)
    total = 0.0
    maxerr = 0.0
    n = 0
    n_node_rows = wake.nwakes[] + 1
    row_total = zeros(Float64, n_node_rows)
    row_max = zeros(Float64, n_node_rows)
    row_count = zeros(Int, n_node_rows)
    for i_surf in eachindex(wake.nodes)
        nodes = wake.nodes[i_surf]
        for jrow in start_row:(wake.nwakes[] + 1), jcol in axes(nodes, 3)
            dx = SVector{3}(
                nodes[1, jrow, jcol] - nodes[1, jrow - 1, jcol],
                nodes[2, jrow, jcol] - nodes[2, jrow - 1, jcol],
                nodes[3, jrow, jcol] - nodes[3, jrow - 1, jcol],
            )
            err = abs(norm(dx) - target_lengths[i_surf][jrow - 1, jcol])
            total += err
            maxerr = max(maxerr, err)
            row_total[jrow] += err
            row_max[jrow] = max(row_max[jrow], err)
            row_count[jrow] += 1
            n += 1
        end
    end
    rows = NamedTuple[]
    for jrow in start_row:n_node_rows
        row_count[jrow] == 0 && continue
        push!(rows, (
            wake_node_row=jrow,
            segment_length_error_mean_R=row_total[jrow] / row_count[jrow] / R,
            segment_length_error_max_R=row_max[jrow] / R,
        ))
    end
    return (
        mean_R = n == 0 ? 0.0 : total / n / R,
        max_R = maxerr / R,
        row_stats = rows,
    )
end

function unit_or(v, fallback)
    vnorm = norm(v)
    if vnorm > sqrt(eps(Float64))
        return v / vnorm
    end
    fnorm = norm(fallback)
    return fnorm > sqrt(eps(Float64)) ? fallback / fnorm : SVector{3}(1.0, 0.0, 0.0)
end

function slerp_unit(a, b, fraction)
    c = clamp(dot(a, b), -1.0, 1.0)
    theta = acos(c)
    if theta < sqrt(eps(Float64))
        return a
    end
    denom = sin(theta)
    out = (sin((1 - fraction) * theta) / denom) * a +
          (sin(fraction * theta) / denom) * b
    return unit_or(out, b)
end

function relax_wake_nodes_length_preserving!(wake, system, target_lengths;
        relax, max_step, start_row::Int=2, mode=effective_near_wake_mode())
    pin_near_wake_row!(wake, system; mode=mode)
    residual_total = 0.0
    residual_max = 0.0
    applied_total = 0.0
    applied_max = 0.0
    capped_nodes = 0
    n = 0
    n_node_rows = wake.nwakes[] + 1
    row_residual_total = zeros(Float64, n_node_rows)
    row_residual_max = zeros(Float64, n_node_rows)
    row_applied_total = zeros(Float64, n_node_rows)
    row_applied_max = zeros(Float64, n_node_rows)
    row_capped_nodes = zeros(Int, n_node_rows)
    row_node_count = zeros(Int, n_node_rows)

    for i_surf in eachindex(wake.nodes)
        nodes = wake.nodes[i_surf]
        vel = wake.velocity[i_surf]
        for jrow in start_row:(wake.nwakes[] + 1), jcol in axes(nodes, 3)
            upstream = SVector{3}(nodes[1, jrow - 1, jcol],
                                  nodes[2, jrow - 1, jcol],
                                  nodes[3, jrow - 1, jcol])
            current = SVector{3}(nodes[1, jrow, jcol],
                                 nodes[2, jrow, jcol],
                                 nodes[3, jrow, jcol])
            current_vec = current - upstream
            L = target_lengths[i_surf][jrow - 1, jcol]
            L <= sqrt(eps(Float64)) && continue
            u_here = SVector{3}(vel[1, jrow, jcol],
                                vel[2, jrow, jcol],
                                vel[3, jrow, jcol])
            current_dir = unit_or(current_vec, u_here)
            target_dir = unit_or(u_here, current_dir)
            target = upstream + L * target_dir
            residual = target - current
            residual_norm = norm(residual)
            is_capped = residual_norm > max_step
            fraction = relax * (is_capped ? max_step / residual_norm : 1.0)
            new_dir = slerp_unit(current_dir, target_dir, clamp(fraction, 0.0, 1.0))
            new_node = upstream + L * new_dir
            disp = new_node - current
            disp_norm = norm(disp)
            nodes[:, jrow, jcol] .= new_node
            residual_total += residual_norm
            residual_max = max(residual_max, residual_norm)
            applied_total += disp_norm
            applied_max = max(applied_max, disp_norm)
            capped_nodes += is_capped ? 1 : 0
            row_residual_total[jrow] += residual_norm
            row_residual_max[jrow] = max(row_residual_max[jrow], residual_norm)
            row_applied_total[jrow] += disp_norm
            row_applied_max[jrow] = max(row_applied_max[jrow], disp_norm)
            row_capped_nodes[jrow] += is_capped ? 1 : 0
            row_node_count[jrow] += 1
            n += 1
        end
    end
    pin_near_wake_row!(wake, system; mode=mode)
    wake.nwakes[] = n_wake_rows
    wake.overflowed[] = true
    length_stats = segment_length_error_stats(wake, target_lengths; start_row=start_row)
    length_rows = Dict(row.wake_node_row => row for row in length_stats.row_stats)
    rows = NamedTuple[]
    for jrow in start_row:n_node_rows
        row_node_count[jrow] == 0 && continue
        length_row = length_rows[jrow]
        push!(rows, (
            wake_node_row=jrow,
            residual_mean_R=row_residual_total[jrow] / row_node_count[jrow] / R,
            residual_max_R=row_residual_max[jrow] / R,
            capped_node_fraction=row_capped_nodes[jrow] / row_node_count[jrow],
            applied_mean_R=row_applied_total[jrow] / row_node_count[jrow] / R,
            applied_max_R=row_applied_max[jrow] / R,
            segment_length_error_mean_R=length_row.segment_length_error_mean_R,
            segment_length_error_max_R=length_row.segment_length_error_max_R,
        ))
    end
    applied_mean_R = n == 0 ? 0.0 : applied_total / n / R
    applied_max_R = applied_max / R
    return (
        residual_mean_R = n == 0 ? 0.0 : residual_total / n / R,
        residual_max_R = residual_max / R,
        capped_node_fraction = n == 0 ? 0.0 : capped_nodes / n,
        applied_mean_R = applied_mean_R,
        applied_max_R = applied_max_R,
        segment_length_error_mean_R = length_stats.mean_R,
        segment_length_error_max_R = length_stats.max_R,
        row_stats = rows,
        mean_R = applied_mean_R,
        max_R = applied_max_R,
    )
end

function snapshot_wake_nodes(wake)
    return [copy(nodes) for nodes in wake.nodes]
end

function restore_wake_nodes!(wake, nodes_snapshot)
    for i_surf in eachindex(wake.nodes)
        wake.nodes[i_surf] .= nodes_snapshot[i_surf]
    end
    return nothing
end

function _frame_material_velocity_at_point(frames, i_system::Int, x_global)
    velocity = zero(x_global)
    identity = SMatrix{3, 3, Float64, 9}(
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
    )
    return _frame_material_velocity_at_point(
        frames, i_system, x_global, 1, zero(SVector{3, Float64}), identity, velocity)
end

function _frame_material_velocity_at_point(frames, i_system::Int, x_global,
        i_frame::Int, dx_parent_to_global, R_parent_to_global, velocity)
    frame = frames[i_frame]
    origin_global = R_parent_to_global * frame.x + dx_parent_to_global
    v_global = R_parent_to_global * frame.v
    omega_global = R_parent_to_global * frame.ω_axis * frame.ω
    if i_system in frame.dependent_index
        velocity += v_global + cross(omega_global, x_global - origin_global)
    end
    dx_parent_to_global = origin_global
    R_parent_to_global = R_parent_to_global * frame.R
    for child in frame.child_index
        velocity = _frame_material_velocity_at_point(
            frames, i_system, x_global, child,
            dx_parent_to_global, R_parent_to_global, velocity)
    end
    return velocity
end

function probe_effective_velocities!(out, points, system, wake, frames, source_nodes)
    nprobe = length(points)
    probes = pnl.FastMultipole.ProbeSystem(nprobe, Float64)
    zero_v = zero(SVector{3, Float64})
    zero_h = zero(SMatrix{3, 3, Float64, 9})
    @inbounds for k in 1:nprobe
        probes.position[k] = points[k]
        probes.gradient[k] = zero_v
        probes.scalar_potential[k] = 0.0
        probes.hessian[k] = zero_h
    end

    current_nodes = snapshot_wake_nodes(wake)
    old_offset = system.kerneloffset
    try
        restore_wake_nodes!(wake, source_nodes)
        wake_sources = pnl.get_sources(wake)
        if length(wake_sources) > 0
            pnl.influence!((probes,), wake_sources, backend;
                precalc=true, scalar_potential=false, gradient=true,
                hessian=(false,))
        end
        system.kerneloffset = system.kerneloffset_targets
        pnl.influence!((probes,), (system,), backend;
            precalc=false, scalar_potential=false, gradient=true,
            hessian=(false,),
            direct_conditioning=pnl._self_panel_kerneloffset_conditioning())
    finally
        system.kerneloffset = old_offset
        restore_wake_nodes!(wake, current_nodes)
    end

    uinf = SVector{3, Float64}(Vinf[1], Vinf[2], Vinf[3])
    @inbounds for k in 1:nprobe
        material_velocity = _frame_material_velocity_at_point(frames, 1, points[k])
        out[k] = probes.gradient[k] + uinf - material_velocity
    end
    return out
end

function relax_wake_nodes_streamline!(wake, system, frames;
        dt_row, substeps::Int=1, start_row::Int=2, mode=effective_near_wake_mode())
    substeps >= 1 || error("STREAMLINE_SUBSTEPS must be >= 1; got $(substeps).")
    pin_near_wake_row!(wake, system; mode=mode)
    source_nodes = snapshot_wake_nodes(wake)
    candidate_nodes = snapshot_wake_nodes(wake)

    displacement_total = 0.0
    displacement_max = 0.0
    n = 0
    n_node_rows = wake.nwakes[] + 1
    row_displacement_total = zeros(Float64, n_node_rows)
    row_displacement_max = zeros(Float64, n_node_rows)
    row_node_count = zeros(Int, n_node_rows)

    for jrow in start_row:n_node_rows
        points = SVector{3, Float64}[]
        refs = Tuple{Int, Int}[]
        for i_surf in eachindex(candidate_nodes)
            nodes = candidate_nodes[i_surf]
            for jcol in axes(nodes, 3)
                push!(points, SVector{3, Float64}(
                    nodes[1, jrow - 1, jcol],
                    nodes[2, jrow - 1, jcol],
                    nodes[3, jrow - 1, jcol],
                ))
                push!(refs, (i_surf, jcol))
            end
        end

        velocities = Vector{SVector{3, Float64}}(undef, length(points))
        midpoint_velocities = similar(velocities)
        dt_sub = dt_row / substeps
        current_points = copy(points)
        for _ in 1:substeps
            probe_effective_velocities!(velocities, current_points, system, wake, frames, source_nodes)
            midpoints = [current_points[k] + 0.5 * dt_sub * velocities[k]
                         for k in eachindex(current_points)]
            probe_effective_velocities!(
                midpoint_velocities, midpoints, system, wake, frames, source_nodes)
            for k in eachindex(current_points)
                current_points[k] = current_points[k] + dt_sub * midpoint_velocities[k]
            end
        end

        for (k, (i_surf, jcol)) in enumerate(refs)
            old_node = SVector{3, Float64}(
                wake.nodes[i_surf][1, jrow, jcol],
                wake.nodes[i_surf][2, jrow, jcol],
                wake.nodes[i_surf][3, jrow, jcol],
            )
            new_node = current_points[k]
            disp = norm(new_node - old_node)
            candidate_nodes[i_surf][:, jrow, jcol] .= new_node
            displacement_total += disp
            displacement_max = max(displacement_max, disp)
            row_displacement_total[jrow] += disp
            row_displacement_max[jrow] = max(row_displacement_max[jrow], disp)
            row_node_count[jrow] += 1
            n += 1
        end
    end

    for i_surf in eachindex(wake.nodes)
        wake.nodes[i_surf] .= candidate_nodes[i_surf]
    end
    pin_near_wake_row!(wake, system; mode=mode)
    wake.nwakes[] = n_wake_rows
    wake.overflowed[] = true

    rows = NamedTuple[]
    for jrow in start_row:n_node_rows
        row_node_count[jrow] == 0 && continue
        row_mean_R = row_displacement_total[jrow] / row_node_count[jrow] / R
        row_max_R = row_displacement_max[jrow] / R
        push!(rows, (
            wake_node_row=jrow,
            residual_mean_R=row_mean_R,
            residual_max_R=row_max_R,
            capped_node_fraction=0.0,
            applied_mean_R=row_mean_R,
            applied_max_R=row_max_R,
        ))
    end
    mean_R = n == 0 ? 0.0 : displacement_total / n / R
    max_R = displacement_max / R
    return (
        residual_mean_R=mean_R,
        residual_max_R=max_R,
        capped_node_fraction=0.0,
        applied_mean_R=mean_R,
        applied_max_R=max_R,
        segment_length_error_mean_R=NaN,
        segment_length_error_max_R=NaN,
        row_stats=rows,
        mean_R=mean_R,
        max_R=max_R,
    )
end

function relax_wake_nodes!(wake, system; relax, pseudo_dt, max_step,
        start_row::Int=2, mode=effective_near_wake_mode())
    pin_near_wake_row!(wake, system; mode=mode)
    residual_total = 0.0
    residual_max = 0.0
    applied_total = 0.0
    applied_max = 0.0
    capped_nodes = 0
    n = 0
    n_node_rows = wake.nwakes[] + 1
    row_residual_total = zeros(Float64, n_node_rows)
    row_residual_max = zeros(Float64, n_node_rows)
    row_applied_total = zeros(Float64, n_node_rows)
    row_applied_max = zeros(Float64, n_node_rows)
    row_capped_nodes = zeros(Int, n_node_rows)
    row_node_count = zeros(Int, n_node_rows)
    for i_surf in eachindex(wake.nodes)
        nodes = wake.nodes[i_surf]
        vel = wake.velocity[i_surf]
        for jrow in start_row:(wake.nwakes[] + 1), jcol in axes(nodes, 3)
            # A steady wake has nonzero convection velocity. Relax toward the
            # age-marched position from the upstream row so node motion has a
            # fixed point, instead of adding U*dt to every node indefinitely.
            upstream = SVector{3}(nodes[1, jrow - 1, jcol],
                                  nodes[2, jrow - 1, jcol],
                                  nodes[3, jrow - 1, jcol])
            current = SVector{3}(nodes[1, jrow, jcol],
                                 nodes[2, jrow, jcol],
                                 nodes[3, jrow, jcol])
            u_up = SVector{3}(vel[1, jrow - 1, jcol],
                              vel[2, jrow - 1, jcol],
                              vel[3, jrow - 1, jcol])
            u_here = SVector{3}(vel[1, jrow, jcol],
                                vel[2, jrow, jcol],
                                vel[3, jrow, jcol])
            target = upstream + 0.5 * (u_up + u_here) * pseudo_dt
            residual = target - current
            residual_norm = norm(residual)
            is_capped = residual_norm > max_step
            capped = is_capped ? residual * (max_step / residual_norm) : residual
            disp = relax * capped
            nodes[:, jrow, jcol] .+= disp
            disp_norm = norm(disp)
            residual_total += residual_norm
            residual_max = max(residual_max, residual_norm)
            applied_total += disp_norm
            applied_max = max(applied_max, disp_norm)
            capped_nodes += is_capped ? 1 : 0
            row_residual_total[jrow] += residual_norm
            row_residual_max[jrow] = max(row_residual_max[jrow], residual_norm)
            row_applied_total[jrow] += disp_norm
            row_applied_max[jrow] = max(row_applied_max[jrow], disp_norm)
            row_capped_nodes[jrow] += is_capped ? 1 : 0
            row_node_count[jrow] += 1
            n += 1
        end
    end
    pin_near_wake_row!(wake, system; mode=mode)
    wake.nwakes[] = n_wake_rows
    wake.overflowed[] = true
    rows = NamedTuple[]
    for jrow in start_row:n_node_rows
        row_node_count[jrow] == 0 && continue
        push!(rows, (
            wake_node_row=jrow,
            residual_mean_R=row_residual_total[jrow] / row_node_count[jrow] / R,
            residual_max_R=row_residual_max[jrow] / R,
            capped_node_fraction=row_capped_nodes[jrow] / row_node_count[jrow],
            applied_mean_R=row_applied_total[jrow] / row_node_count[jrow] / R,
            applied_max_R=row_applied_max[jrow] / R,
        ))
    end
    applied_mean_R = n == 0 ? 0.0 : applied_total / n / R
    applied_max_R = applied_max / R
    return (
        residual_mean_R = n == 0 ? 0.0 : residual_total / n / R,
        residual_max_R = residual_max / R,
        capped_node_fraction = n == 0 ? 0.0 : capped_nodes / n,
        applied_mean_R = applied_mean_R,
        applied_max_R = applied_max_R,
        row_stats = rows,
        # Backward-compatible names; these are applied, capped displacements.
        mean_R = applied_mean_R,
        max_R = applied_max_R,
    )
end

function evaluate_wake_velocity_with_pinned_strengths!(system, wake)
    pnl.reset!(wake)
    pnl.apply_freestream!(wake, Vinf)
    wake_probes = pnl.get_probes(wake)
    wake_sources = pnl.get_sources(wake)
    if length(wake_sources) > 0
        pnl.influence!(wake_probes, wake_sources, backend;
            precalc=true, scalar_potential=false, velocity=true,
            velocity_gradient=false)
    end
    pnl._set_kerneloffsets!((system,), :kerneloffset_targets)
    pnl.influence!(wake_probes, (system,), backend;
        precalc=false, scalar_potential=false, velocity=true,
        velocity_gradient=false,
        direct_conditioning=pnl._self_panel_kerneloffset_conditioning())
    return nothing
end

function run_steady_with_explicit_wake!(system, wake, frames, iteration; write_state::Bool,
        update_trailing_edges::Bool=false, tag_suffix::AbstractString="",
        vtk_series::Bool=false, vtk_iteration::Int=iteration,
        vtk_series_suffix::AbstractString="iter", vtk_series_overwrite::Bool=false,
        vtk_final::Bool=false)
    pressure = pnl.PressureBernoulli(rho;
        unsteady=false, correct_kuttacondition=false, backend=backend, file=false)
    force = pnl.ForceMonitor(1, 1; i_frame=1,
        normalization=pnl.RotorNormalization(rho, 2 * R, 1),
        correct_kuttacondition=false, verbose=false, file=false)
    pressure_laplace_lamb = pnl.PressureLaplace((system,), rho;
        unsteady=false, acceleration_form=:lamb_vector, file=false)
    force_laplace_lamb = pnl.ForceMonitor(1, 1; i_frame=1,
        normalization=pnl.RotorNormalization(rho, 2 * R, 1),
        correct_kuttacondition=false, verbose=false, file=false)
    vorticity_force = pnl.SurfaceVorticityForce(system, 1, 1; rho=rho, i_frame=1,
        normalization=pnl.RotorNormalization(rho, 2 * R, 1),
        correct_kuttacondition=false, verbose=false, file=false)
    kj = pnl.KuttaJoukowskiForce(system, 1, 1; rho=rho, backend=backend,
        i_frame=1, normalization=pnl.RotorNormalization(rho, 2 * R, 1),
        verbose=false, file=false)
    monitors = (pressure, force, pressure_laplace_lamb, force_laplace_lamb,
        vorticity_force, kj)

    systems_tuple = (system,)
    wakes_tuple = (wake,)
    pnl.audit_monitors(monitors)
    for sys in systems_tuple
        sys.needs_velocity_gradient[] = any(pnl.monitor_requires_body_hessian, monitors)
        pnl.calc_normals!(sys)
        pnl.calc_controlpoints!(sys)
    end

    pnl._steady_aerodynamics!(systems_tuple, systems_tuple, wakes_tuple, frames, Vinf,
        (pnl.Backslash(system),);
        backend_wake=backend, backend_solve=backend, backend_system=backend,
        needs_induced_vorticity=any(pnl.monitor_requires_induced_vorticity, monitors),
        update_trailing_edges)

    ctx = pnl.MonitorContext()
    pnl.monitor_set_time!(ctx, Float64(iteration - 1))
    for monitor in monitors
        pnl._run_monitor!(monitor, ctx, systems_tuple, wakes_tuple, frames, Vinf, 0, 1.0)
    end

    if write_state
        if vtk_series
            tag = "$(basename(run_name))_body1_$(vtk_series_suffix)"
            pnl.write_vtk(joinpath(save_path, tag), system, vtk_iteration, vtk_iteration;
                monitors=monitors, i_system=1, overwrite=vtk_series_overwrite,
                compress=true)
            pnl.write_vtk(joinpath(save_path, "$(basename(run_name))_wake_$(vtk_series_suffix)"),
                wake, vtk_iteration, vtk_iteration;
                overwrite=vtk_series_overwrite, compress=true,
                filament_name=joinpath(save_path,
                    "$(basename(run_name))_wake_filaments_$(vtk_series_suffix)"))
        elseif vtk_final
            pnl.write_vtk(joinpath(save_path, "$(basename(run_name))_body1"),
                system, 0, 0.0;
                monitors=monitors, i_system=1, overwrite=true, compress=true)
            pnl.write_vtk(joinpath(save_path, "$(basename(run_name))_wake"),
                wake, 0, 0.0; overwrite=true, compress=true)
        else
            suffix = isempty(tag_suffix) ? "" : "_$(tag_suffix)"
            tag = "$(basename(run_name))_iter$(lpad(iteration, 2, '0'))$(suffix)"
            pnl.write_vtk(joinpath(save_path, tag * "_body1"), system, iteration - 1, iteration - 1;
                monitors=monitors, i_system=1, overwrite=iteration == 1, compress=true)
            pnl.write_vtk(joinpath(save_path, tag * "_wake"), wake, iteration - 1, iteration - 1;
                overwrite=iteration == 1, compress=true)
        end
    end

    return (
        CT_bernoulli = -force.force[axial_dimension, 1],
        CT_laplace_lamb = -force_laplace_lamb.force[axial_dimension, 1],
        CT_vorticity = -vorticity_force.force[axial_dimension, 1],
        CT_kj = -kj.force[axial_dimension, 1],
        wake_row_inflow = wake_row_inflow_ratio(wake),
    )
end

function write_iteration_table(rows)
    ntget(r, name::Symbol, default) = name in propertynames(r) ? getproperty(r, name) : default
    csv_path = joinpath(save_path, "iteration_table.csv")
    open(csv_path, "w") do io
        println(io, "iteration,outer_iteration,inner_iteration,wake_revs,rows_per_rev,nwake_panels,CT_bernoulli,CT_laplace_lamb,CT_vorticity,CT_kj,delta_CT,mean_abs_gamma,max_abs_gamma,mean_abs_gamma_candidate,max_abs_gamma_candidate,delta_gamma_rel,strength_checksum_delta,mean_node_step_R,max_node_step_R,mean_target_residual_R,max_target_residual_R,capped_node_fraction,mean_applied_step_R,max_applied_step_R,mean_segment_length_error_R,max_segment_length_error_R,wake_row_inflow_ratio,fountain_node_ux_negative,fountain_panel_ux_negative,fountain_nodes_upstream_of_anchor,fountain_panels_crossing_anchor_plane,fountain_rows_with_upstream_motion,fountain_min_anchor_distance_R,fountain_min_signed_velocity_tip,fountain_min_panel_signed_velocity_tip,fountain_min_raw_x_R,fountain_min_raw_ux_tip,converged,stopping_reason")
        for r in rows
            diag = ntget(r, :fountain, (
                node_ux_negative=0, panel_ux_negative=0,
                nodes_upstream_of_anchor=0, panels_crossing_anchor_plane=0,
                rows_with_upstream_motion=0, min_anchor_distance_R=NaN,
                min_signed_velocity_tip=NaN, min_panel_signed_velocity_tip=NaN,
                min_raw_x_R=NaN, min_raw_ux_tip=NaN,
            ))
            println(io, @sprintf("%d,%d,%d,%.8f,%d,%d,%.8f,%.8f,%.8f,%.8f,%.8f,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8f,%.8e,%.8e,%.8e,%.8e,%.8f,%d,%d,%d,%d,%d,%.8e,%.8e,%.8e,%.8e,%.8e,%s,%s",
                r.iteration, ntget(r, :outer_iteration, r.iteration),
                ntget(r, :inner_iteration, 0),
                r.wake_revs, r.rows_per_rev, r.nwake_panels,
                r.CT_bernoulli, r.CT_laplace_lamb, r.CT_vorticity, r.CT_kj,
                r.delta_CT,
                r.mean_abs_gamma, r.max_abs_gamma,
                r.mean_abs_gamma_candidate, r.max_abs_gamma_candidate,
                r.delta_gamma_rel, ntget(r, :strength_checksum_delta, NaN),
                r.mean_node_step_R, r.max_node_step_R,
                ntget(r, :mean_target_residual_R, r.mean_node_step_R),
                ntget(r, :max_target_residual_R, r.max_node_step_R),
                ntget(r, :capped_node_fraction, NaN),
                ntget(r, :mean_applied_step_R, r.mean_node_step_R),
                ntget(r, :max_applied_step_R, r.max_node_step_R),
                ntget(r, :mean_segment_length_error_R, NaN),
                ntget(r, :max_segment_length_error_R, NaN),
                r.wake_row_inflow_ratio,
                diag.node_ux_negative, diag.panel_ux_negative,
                diag.nodes_upstream_of_anchor, diag.panels_crossing_anchor_plane,
                diag.rows_with_upstream_motion, diag.min_anchor_distance_R,
                diag.min_signed_velocity_tip, diag.min_panel_signed_velocity_tip,
                diag.min_raw_x_R, diag.min_raw_ux_tip,
                string(r.converged), r.stopping_reason))
        end
    end
    return csv_path
end

function write_row_residual_table(rows)
    ntget(r, name::Symbol, default) = name in propertynames(r) ? getproperty(r, name) : default
    csv_path = joinpath(save_path, "row_residual_table.csv")
    open(csv_path, "w") do io
        println(io, "iteration,outer_iteration,inner_iteration,wake_node_row,residual_mean_R,residual_max_R,capped_node_fraction,applied_mean_R,applied_max_R,segment_length_error_mean_R,segment_length_error_max_R")
        for r in rows
            for row in ntget(r, :row_residuals, NamedTuple[])
                println(io, @sprintf("%d,%d,%d,%d,%.8e,%.8e,%.8f,%.8e,%.8e,%.8e,%.8e",
                    r.iteration, ntget(r, :outer_iteration, r.iteration),
                    ntget(r, :inner_iteration, 0), row.wake_node_row,
                    row.residual_mean_R, row.residual_max_R,
                    row.capped_node_fraction, row.applied_mean_R,
                    row.applied_max_R,
                    ntget(row, :segment_length_error_mean_R, NaN),
                    ntget(row, :segment_length_error_max_R, NaN)))
            end
        end
    end
    return csv_path
end

function main()
    validate_iteration_mode(iteration_mode)
    validate_inner_geometry_mode(inner_geometry_mode)
    validate_initial_strength_mode(initial_strength_mode)
    near_mode = effective_near_wake_mode()
    wake.nwakes[] = n_wake_rows
    wake.overflowed[] = true
    populate_helical_nodes!(wake, rotor; pitch_inflow_ratio=initial_helix_pitch_ratio,
        mode=near_mode)
    fill!.(wake.strength, 0.0)
    pin_near_wake_row!(wake, rotor; mode=near_mode)
    initial_stats = initial_wake_direction_stats(wake)
    write_initial_vtk!(wake, rotor)

    println("\nBagai/Leishman-style PanelWake relaxation probe - $(mesh_tag), $(RPM) RPM")
    println("  iteration mode: $(iteration_mode)")
    if iteration_mode == :nested_pitch
        println("  inner geometry mode: $(inner_geometry_mode), streamline substeps: $(streamline_substeps)")
    end
    println(@sprintf("  wake %.2f revs, %d rows/rev, %d panel rows, core %.3e m",
        wake_revs, rows_per_rev, wake.nwakes[], wake_core_size))
    println(@sprintf("  near wake: %s, helix winding sign %.1f, axial wake sign %.1f",
        string(near_mode), helix_winding_sign, axial_wake_sign))
    println(@sprintf("  axial advance ratio %.4f, freestream axial speed %.4f m/s (%.4f tip)",
        axial_advance_ratio, Vinf[axial_dimension],
        abs(Vinf[axial_dimension]) / tip_speed))
    println(@sprintf("  induced seed ratio %.4f, total initial helix pitch ratio %.4f",
        initial_inflow, initial_helix_pitch_ratio))
    println("  initial strength mode: $(initial_strength_mode)")
    println(@sprintf("  initial row-2 anchor distance: mean %.4e R, min %.4e R, max %.4e R, upstream nodes %d/%d",
        initial_stats.mean_R, initial_stats.min_R, initial_stats.max_R,
        initial_stats.upstream_nodes, initial_stats.n_nodes))
    if initial_stats.upstream_nodes > 0
        println("  WARNING: initialized wake has row-2 nodes upstream of their local TE anchors.")
    end
    println(@sprintf("  relax: wake %.3f, strength %.3f, pseudo dt %.3e s, max node step %.4f R",
        wake_relax, strength_relax, pseudo_dt, max_node_step_R))
    println(@sprintf("  anchor-relative fountain tolerance: %.4e R", fountain_tol_R))
    if debug_vtk_every_iter
        println("  WARNING: DEBUG_VTK_EVERY_ITER=true writes current body/wake VTK for every inner and outer iteration without recomputing pressure or force monitor fields; short runs can create many files.")
    end
    println("  Wrote initial body/wake VTK snapshots under $(save_path)/")
    println("  References: steady rigid wake CT ~= 0.0505; VPM ~= 0.062; BEM ~= 0.068; experiment ~= 0.072")
    if initial_only
        println("  INITIAL_ONLY=true: stopping after initial VTK write.")
        return nothing
    end

    if iteration_mode != :nested_pitch && initial_strength_mode != :zero
        seed_result = run_steady_with_explicit_wake!(rotor, wake, frames, 0;
            write_state=false, update_trailing_edges=true)
        seed_candidates = candidate_strengths(rotor)
        seed_candidate_stats = strength_stats(seed_candidates)
        seed_stats = set_wake_strengths_from_candidates!(wake, seed_candidates;
            mode=initial_strength_mode)
        println(@sprintf("  Seeded wake strengths from first solve: CT(Bern) %.5f, candidate mean %.3e, wake mean %.3e, wake max %.3e",
            seed_result.CT_bernoulli, seed_candidate_stats.mean, seed_stats.mean,
            seed_stats.max))
    end

    rows = NamedTuple[]
    previous_ct = NaN
    previous_growth = (ct = NaN, gamma = NaN, node = NaN)
    unstable_streak = 0
    converged = false
    stopping_reason = "max_iter"
    global_iteration = 0

    if iteration_mode == :nested_pitch
        println("\n  outer inner  CT(Bern)  CT(Lamb)  CT(vort)  CT(KJ)    dCT      dGamma   resMax/R  capFrac  stepMax/R  gammaMean  candMean  rowInflow  fountain status")
        println("  " * "-"^160)

        for outer in 1:max_outer_iter
            vtk_outer_suffix = "iter$(lpad(outer, 4, '0'))"
            last_inner_iteration = 0
            inner_converged = false
            global_iteration += 1
            wake.nwakes[] = n_wake_rows
            pin_near_wake_row!(wake, rotor; mode=near_mode)
            result = run_steady_with_explicit_wake!(rotor, wake, frames, global_iteration;
                write_state=false, update_trailing_edges=false)
            candidates = candidate_strengths(rotor)
            candidate_stats = strength_stats(candidates)
            delta_gamma_rel = pin_latest_strengths!(wake, candidates)
            relaxed_stats = wake_strength_stats(wake)
            pinned_strengths = snapshot_strengths(wake)
            pinned_body_strength = snapshot_body_strength(rotor)
            fixed_lengths = snapshot_streamwise_segment_lengths(wake)
            node_stats = empty_node_relaxation_stats()
            strength_delta = 0.0
            diag = fountain_flow_diagnostic(wake)

            if debug_vtk_every_iter
                wake.nwakes[] = n_wake_rows
                write_current_paraview_state!(rotor, wake;
                    vtk_iteration=0, vtk_series_suffix=vtk_outer_suffix,
                    overwrite=true)
                assert_strength_unchanged(wake, pinned_strengths;
                    context="outer $(outer) debug write")
                assert_body_strength_unchanged(rotor, pinned_body_strength;
                    context="outer $(outer) debug write")
            end

            for inner in 1:max_inner_iter
                last_inner_iteration = inner
                global_iteration += 1
                wake.nwakes[] = n_wake_rows
                if inner_geometry_mode == :streamline
                    node_stats = relax_wake_nodes_streamline!(wake, rotor, frames;
                        dt_row=wake_row_dt, substeps=streamline_substeps,
                        start_row=2, mode=near_mode)
                else
                    evaluate_wake_velocity_with_pinned_strengths!(rotor, wake)
                    node_stats = relax_wake_nodes_length_preserving!(wake, rotor, fixed_lengths;
                        relax=wake_relax, max_step=max_node_step_R * R,
                        start_row=2, mode=near_mode)
                end
                strength_delta = assert_strength_unchanged(wake, pinned_strengths;
                    context="outer $(outer), inner $(inner)")
                evaluate_wake_velocity_with_pinned_strengths!(rotor, wake)
                diag = fountain_flow_diagnostic(wake)
                row_inflow = wake_row_inflow_ratio(wake)
                delta_ct = isfinite(previous_ct) ? abs(result.CT_bernoulli - previous_ct) : Inf

                if diag.fountain
                    stopping_reason = "fountain_flow"
                    push!(rows, (
                        iteration=global_iteration, outer_iteration=outer,
                        inner_iteration=inner, wake_revs=wake_revs,
                        rows_per_rev=rows_per_rev,
                        nwake_panels=sum(size(s, 3) for s in wake.strength) * wake.nwakes[],
                        CT_bernoulli=result.CT_bernoulli,
                        CT_laplace_lamb=result.CT_laplace_lamb,
                        CT_vorticity=result.CT_vorticity, CT_kj=result.CT_kj,
                        delta_CT=delta_ct,
                        mean_abs_gamma=relaxed_stats.mean,
                        max_abs_gamma=relaxed_stats.max,
                        mean_abs_gamma_candidate=candidate_stats.mean,
                        max_abs_gamma_candidate=candidate_stats.max,
                        delta_gamma_rel=delta_gamma_rel,
                        strength_checksum_delta=strength_delta,
                        mean_node_step_R=node_stats.mean_R,
                        max_node_step_R=node_stats.max_R,
                        mean_target_residual_R=node_stats.residual_mean_R,
                        max_target_residual_R=node_stats.residual_max_R,
                        capped_node_fraction=node_stats.capped_node_fraction,
                        mean_applied_step_R=node_stats.applied_mean_R,
                        max_applied_step_R=node_stats.applied_max_R,
                        mean_segment_length_error_R=node_stats.segment_length_error_mean_R,
                        max_segment_length_error_R=node_stats.segment_length_error_max_R,
                        row_residuals=node_stats.row_stats,
                        wake_row_inflow_ratio=row_inflow,
                        fountain=diag,
                        converged=false, stopping_reason=stopping_reason,
                    ))
                    if debug_vtk_every_iter
                        write_current_paraview_state!(rotor, wake;
                            vtk_iteration=inner, vtk_series_suffix=vtk_outer_suffix)
                    else
                        run_steady_with_explicit_wake!(rotor, wake, frames, global_iteration;
                            write_state=true, update_trailing_edges=false,
                            tag_suffix="fountain_outer$(outer)_inner$(inner)")
                    end
                    assert_strength_unchanged(wake, pinned_strengths;
                        context="outer $(outer), fountain write")
                    if debug_vtk_every_iter
                        assert_body_strength_unchanged(rotor, pinned_body_strength;
                            context="outer $(outer), fountain write")
                    end
                    println(@sprintf("  %5d %5d  %8.5f  %8.5f  %8.5f  %8.5f  %7.4f  %8.3e  %9.4f  %7.3f  %9.4f  %9.3e  %8.3e  %9.4f  %8d %s",
                        outer, inner, result.CT_bernoulli, result.CT_laplace_lamb,
                        result.CT_vorticity, result.CT_kj,
                        isfinite(delta_ct) ? delta_ct : NaN, delta_gamma_rel,
                        node_stats.residual_max_R, node_stats.capped_node_fraction,
                        node_stats.applied_max_R, relaxed_stats.mean, candidate_stats.mean,
                        isfinite(row_inflow) ? row_inflow : NaN,
                        diag.rows_with_upstream_motion, "fountain"))
                    break
                end

                inner_converged = node_stats.residual_max_R < node_tol_R
                push!(rows, (
                    iteration=global_iteration, outer_iteration=outer,
                    inner_iteration=inner, wake_revs=wake_revs,
                    rows_per_rev=rows_per_rev,
                    nwake_panels=sum(size(s, 3) for s in wake.strength) * wake.nwakes[],
                    CT_bernoulli=result.CT_bernoulli,
                    CT_laplace_lamb=result.CT_laplace_lamb,
                    CT_vorticity=result.CT_vorticity, CT_kj=result.CT_kj,
                    delta_CT=delta_ct,
                    mean_abs_gamma=relaxed_stats.mean,
                    max_abs_gamma=relaxed_stats.max,
                    mean_abs_gamma_candidate=candidate_stats.mean,
                    max_abs_gamma_candidate=candidate_stats.max,
                    delta_gamma_rel=delta_gamma_rel,
                    strength_checksum_delta=strength_delta,
                    mean_node_step_R=node_stats.mean_R,
                    max_node_step_R=node_stats.max_R,
                    mean_target_residual_R=node_stats.residual_mean_R,
                    max_target_residual_R=node_stats.residual_max_R,
                    capped_node_fraction=node_stats.capped_node_fraction,
                    mean_applied_step_R=node_stats.applied_mean_R,
                    max_applied_step_R=node_stats.applied_max_R,
                    mean_segment_length_error_R=node_stats.segment_length_error_mean_R,
                    max_segment_length_error_R=node_stats.segment_length_error_max_R,
                    row_residuals=node_stats.row_stats,
                    wake_row_inflow_ratio=row_inflow,
                    fountain=diag,
                    converged=false,
                    stopping_reason=inner_converged ? "inner_converged" : "pitch_relax",
                ))
                println(@sprintf("  %5d %5d  %8.5f  %8.5f  %8.5f  %8.5f  %7.4f  %8.3e  %9.4f  %7.3f  %9.4f  %9.3e  %8.3e  %9.4f  %8d %s",
                    outer, inner, result.CT_bernoulli, result.CT_laplace_lamb,
                    result.CT_vorticity, result.CT_kj,
                    isfinite(delta_ct) ? delta_ct : NaN,
                    delta_gamma_rel, node_stats.residual_max_R,
                    node_stats.capped_node_fraction, node_stats.applied_max_R,
                    relaxed_stats.mean, candidate_stats.mean,
                    isfinite(row_inflow) ? row_inflow : NaN,
                    diag.rows_with_upstream_motion,
                    inner_converged ? "inner-converged" : "pitch-relax"))
                if debug_vtk_every_iter
                    write_current_paraview_state!(rotor, wake;
                        vtk_iteration=inner, vtk_series_suffix=vtk_outer_suffix)
                    assert_strength_unchanged(wake, pinned_strengths;
                        context="outer $(outer), inner $(inner) debug write")
                    assert_body_strength_unchanged(rotor, pinned_body_strength;
                        context="outer $(outer), inner $(inner) debug write")
                end
                inner_converged && break
            end

            stopping_reason == "fountain_flow" && break

            delta_ct = isfinite(previous_ct) ? abs(result.CT_bernoulli - previous_ct) : Inf
            converged = outer > 1 &&
                delta_ct < ct_tol &&
                delta_gamma_rel < gamma_tol &&
                node_stats.residual_max_R < node_tol_R
            stopping_reason = converged ? "converged" : "max_iter"

            status = converged ? "converged" : (inner_converged ? "outer-updated" : "inner-max")
            println(@sprintf("  %5d %5s  %8.5f  %8.5f  %8.5f  %8.5f  %7.4f  %8.3e  %9.4f  %7.3f  %9.4f  %9.3e  %8.3e  %9.4f  %8d %s",
                outer, "post", result.CT_bernoulli, result.CT_laplace_lamb,
                result.CT_vorticity, result.CT_kj,
                isfinite(delta_ct) ? delta_ct : NaN, delta_gamma_rel,
                node_stats.residual_max_R, node_stats.capped_node_fraction,
                node_stats.applied_max_R, relaxed_stats.mean, candidate_stats.mean,
                isfinite(wake_row_inflow_ratio(wake)) ? wake_row_inflow_ratio(wake) : NaN,
                diag.rows_with_upstream_motion, status))

            if debug_vtk_every_iter
                write_current_paraview_state!(rotor, wake;
                    vtk_iteration=0, vtk_series_suffix=vtk_outer_suffix, final=true)
                assert_strength_unchanged(wake, pinned_strengths;
                    context="outer $(outer) final write")
                assert_body_strength_unchanged(rotor, pinned_body_strength;
                    context="outer $(outer) final write")
            elseif save_vtk
                run_steady_with_explicit_wake!(rotor, wake, frames, global_iteration;
                    write_state=true, update_trailing_edges=false,
                    tag_suffix="outer$(outer)_post")
                assert_strength_unchanged(wake, pinned_strengths;
                    context="outer $(outer) final write")
            end
            previous_ct = result.CT_bernoulli
            (converged || outer >= max_outer_iter) && break
        end
    else
        println("\n  iter  CT(Bern)  CT(Lamb)  CT(vort)  CT(KJ)    dCT      dGamma   resMax/R  capFrac  stepMax/R  gammaMean  candMean  rowInflow  fountain status")
        println("  " * "-"^151)

        for iter in 1:max_iter
            wake.nwakes[] = n_wake_rows
            result = run_steady_with_explicit_wake!(rotor, wake, frames, iter;
                write_state=save_vtk || debug_vtk_every_iter,
                tag_suffix="iter$(iter)",
                vtk_series=debug_vtk_every_iter,
                vtk_iteration=iter,
                vtk_series_overwrite=iter == 1)
            candidates = candidate_strengths(rotor)
            candidate_stats = strength_stats(candidates)
            relaxed_strengths, delta_gamma_rel = relax_strengths!(wake, candidates; relax=strength_relax)
            relaxed_stats = strength_stats(relaxed_strengths)
            node_stats = relax_wake_nodes!(wake, rotor;
                relax=wake_relax, pseudo_dt=pseudo_dt, max_step=max_node_step_R * R,
                mode=near_mode)
            diag = fountain_flow_diagnostic(wake)

            delta_ct = isfinite(previous_ct) ? abs(result.CT_bernoulli - previous_ct) : Inf
            row_inflow = result.wake_row_inflow

            converged = iter > 1 &&
                delta_ct < ct_tol &&
                delta_gamma_rel < gamma_tol &&
                node_stats.residual_max_R < node_tol_R
            if converged
                stopping_reason = "converged"
            elseif diag.fountain
                stopping_reason = "fountain_flow"
            end

            ct_abs = abs(result.CT_bernoulli)
            grows_ct = isfinite(previous_growth.ct) && ct_abs > 1.35 * max(previous_growth.ct, eps())
            grows_gamma = isfinite(previous_growth.gamma) && relaxed_stats.mean > 1.35 * max(previous_growth.gamma, eps())
            grows_node = isfinite(previous_growth.node) && node_stats.residual_max_R > 1.35 * max(previous_growth.node, eps())
            if grows_ct || grows_gamma || grows_node || !isfinite(result.CT_bernoulli) || !isfinite(delta_gamma_rel)
                unstable_streak += 1
            else
                unstable_streak = 0
            end
            unstable = unstable_streak >= 3
            if unstable
                stopping_reason = "unstable_or_divergent"
            end

            push!(rows, (
                iteration=iter, wake_revs=wake_revs, rows_per_rev=rows_per_rev,
                nwake_panels=sum(size(s, 3) for s in wake.strength) * wake.nwakes[],
                CT_bernoulli=result.CT_bernoulli,
                CT_laplace_lamb=result.CT_laplace_lamb,
                CT_vorticity=result.CT_vorticity, CT_kj=result.CT_kj,
                delta_CT=delta_ct,
                mean_abs_gamma=relaxed_stats.mean, max_abs_gamma=relaxed_stats.max,
                mean_abs_gamma_candidate=candidate_stats.mean,
                max_abs_gamma_candidate=candidate_stats.max,
                delta_gamma_rel=delta_gamma_rel,
                mean_node_step_R=node_stats.mean_R, max_node_step_R=node_stats.max_R,
                mean_target_residual_R=node_stats.residual_mean_R,
                max_target_residual_R=node_stats.residual_max_R,
                capped_node_fraction=node_stats.capped_node_fraction,
                mean_applied_step_R=node_stats.applied_mean_R,
                max_applied_step_R=node_stats.applied_max_R,
                row_residuals=node_stats.row_stats,
                wake_row_inflow_ratio=row_inflow,
                fountain=diag,
                converged=converged, stopping_reason=stopping_reason,
            ))

            status = diag.fountain ? "fountain" : (converged ? "converged" : (unstable ? "unstable" : "relaxing"))
            println(@sprintf("  %4d  %8.5f  %8.5f  %8.5f  %8.5f  %7.4f  %8.3e  %9.4f  %7.3f  %9.4f  %9.3e  %8.3e  %9.4f  %8d %s",
                iter, result.CT_bernoulli, result.CT_laplace_lamb,
                result.CT_vorticity, result.CT_kj,
                isfinite(delta_ct) ? delta_ct : NaN, delta_gamma_rel,
                node_stats.residual_max_R, node_stats.capped_node_fraction,
                node_stats.applied_max_R, relaxed_stats.mean, candidate_stats.mean,
                isfinite(row_inflow) ? row_inflow : NaN,
                diag.rows_with_upstream_motion, status))

            if diag.fountain
                run_steady_with_explicit_wake!(rotor, wake, frames, iter;
                    write_state=true, tag_suffix="fountain")
            end

            previous_ct = result.CT_bernoulli
            previous_growth = (ct = ct_abs, gamma = relaxed_stats.mean, node = node_stats.residual_max_R)
            (converged || unstable || diag.fountain) && break
        end
    end

    iter_limit = iteration_mode == :nested_pitch ? max_outer_iter : max_iter
    if !isempty(rows) && rows[end].stopping_reason == "max_iter" && length(rows) >= iter_limit
        rows[end] = merge(rows[end], (stopping_reason="max_iter",))
    end

    csv_path = write_iteration_table(rows)
    row_csv_path = write_row_residual_table(rows)
    last = rows[end]
    println("\nWrote $(csv_path)")
    println("Wrote $(row_csv_path)")
    if save_vtk || debug_vtk_every_iter || last.stopping_reason == "fountain_flow"
        println("Wrote body/wake VTK snapshots under $(save_path)/")
    end

    println("\nVerdict:")
    if iteration_mode == :nested_pitch
        max_strength_delta = max_row_value(rows, :strength_checksum_delta)
        max_length_error_R = max_row_value(rows, :max_segment_length_error_R)
        println(@sprintf("  Inner-loop invariants: max strength checksum delta %.3e, max streamwise segment length error %.3e R.",
            max_strength_delta, max_length_error_R))
    end
    if last.stopping_reason == "fountain_flow"
        diag = last.fountain
        println(@sprintf("  Fountain flow detected: %d downstream nodes upstream of their local anchors, %d anchor-plane crossing panels, %d downstream rows with mean signed velocity < 0.",
            diag.nodes_upstream_of_anchor, diag.panels_crossing_anchor_plane,
            diag.rows_with_upstream_motion))
        println(@sprintf("  Stopped before accepting convergence or CT trends. Minimum anchor distance/R %.4e, minimum signed velocity/tip %.4e.",
            diag.min_anchor_distance_R, diag.min_signed_velocity_tip))
    elseif last.stopping_reason == "unstable_or_divergent"
        println(@sprintf("  Unstable/divergent: CT %.5f, delta_gamma %.3e, max target residual %.4f R (applied step %.4f R).",
            last.CT_bernoulli, last.delta_gamma_rel, last.max_target_residual_R,
            last.max_applied_step_R))
    elseif last.converged && 0.060 <= last.CT_bernoulli <= 0.080
        println(@sprintf("  Converged at CT %.5f (Laplace Lamb %.5f), in the target VPM/BEM/experiment neighborhood.",
            last.CT_bernoulli, last.CT_laplace_lamb))
    elseif last.converged && last.CT_bernoulli < 0.060
        println(@sprintf("  Converged at CT %.5f (Laplace Lamb %.5f), below the target VPM/BEM/experiment neighborhood.",
            last.CT_bernoulli, last.CT_laplace_lamb))
    elseif last.converged
        println(@sprintf("  Converged at CT %.5f (Laplace Lamb %.5f), above the target VPM/BEM/experiment neighborhood.",
            last.CT_bernoulli, last.CT_laplace_lamb))
    else
        iter_label = iteration_mode == :nested_pitch ? "MAX_OUTER_ITER" : "MAX_ITER"
        println(@sprintf("  Not converged within %s=%d: CT %.5f (Laplace Lamb %.5f), delta_gamma %.3e, max target residual %.4f R (applied step %.4f R, capped %.1f%%).",
            iter_label, iter_limit, last.CT_bernoulli, last.CT_laplace_lamb,
            last.delta_gamma_rel, last.max_target_residual_R,
            last.max_applied_step_R, 100 * last.capped_node_fraction))
        if iteration_mode == :nested_pitch
            println("  Inner residual classification: $(classify_inner_failure(rows; tol_R=node_tol_R)).")
        end
        if :fountain in propertynames(last) && last.fountain.rows_with_upstream_motion > 0
            println(@sprintf("  Warning metric: %d downstream rows had mean signed velocity < 0, but no anchor-relative wake crossing was detected.",
                last.fountain.rows_with_upstream_motion))
        end
    end
end

main()
