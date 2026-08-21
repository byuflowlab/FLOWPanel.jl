## TE downwash probe (BRAINSTORM/018 Phase 16 F1b, Route B).
## Replays the last saved steps of a SETTLED run and measures the induced
## drift velocity at the trailing-edge stations (or at the sheet midpoint
## TE + 0.5*Das), averaged over the replayed steps and over blades, writing
## the per-station drift table CSV consumed by DAS_ARC_TABLE.
##
## The drift is the INERTIAL fluid velocity u = u_induced + Uinf (NOT minus
## the kinematic blade velocity — the frozen-path construction handles blade
## kinematics via the backward swept arc; see _set_Das_station_arc!).
## Components are stored in the local TE basis so the table co-rotates with
## the geometry by construction:
##   e_t = normalized kinematic TE tangent (-omega_hat x d_perp)
##   e_r = radial/spanwise unit
##   e_n = e_t x e_r
##
## Env knobs:
##   P018_PROBE_SOURCE    run dir under data/ (default p018_cs_l3p4_rs1)
##   P018_PROBE_RUN_NAME  saved run_name inside that dir (default = SOURCE)
##   P018_PROBE_LOCATION  te | mid  (default te)
##   P018_PROBE_STEPS     "lo:hi" saved-step range, or "auto" = last 10
##   P018_PROBE_CLAMP_DIR "x,y,z" GLOBAL direction; any station drift with a
##                        positive component along it gets that component
##                        zeroed (Ryan 2026-08-18: user-specified clamp;
##                        expected inactive on a developed wake). Optional.
##   P018_PROBE_OUT       output CSV path (default derived)

import FLOWPanel as pnl
using LinearAlgebra: norm, dot, cross
using StaticArrays: SVector, SMatrix
import FastMultipole

ENV["RHPC_SETUP_ONLY"] = "true"
include(joinpath(pnl.examples_path, "rotor_hover_pressure_comparison.jl"))

probe_source = get(ENV, "P018_PROBE_SOURCE", "p018_cs_l3p4_rs1")
probe_run_name = get(ENV, "P018_PROBE_RUN_NAME", probe_source)
probe_loc = Symbol(lowercase(get(ENV, "P018_PROBE_LOCATION", "te")))
probe_loc in (:te, :mid) || error("P018_PROBE_LOCATION must be te or mid")
steps_spec = get(ENV, "P018_PROBE_STEPS", "auto")
out_csv = get(ENV, "P018_PROBE_OUT",
    joinpath("data", "$(probe_source)_te_downwash_$(probe_loc).csv"))
clamp_spec = get(ENV, "P018_PROBE_CLAMP_DIR", "")

source_path = joinpath("data", probe_source)
isdir(source_path) || error("source run dir not found: $(source_path)")

# --- step selection ---
probe_steps = if steps_spec == "auto"
    body_dir = joinpath(source_path, "$(probe_run_name)_body1")
    isdir(body_dir) || error("no body1 VTU dir at $(body_dir)")
    idxs = sort(unique(parse(Int, m.captures[1]) for f in readdir(body_dir)
        for m in (match(r"\.(\d+)\.vtu$", f),) if m !== nothing))
    length(idxs) >= 2 || error("found $(length(idxs)) VTU steps in $(body_dir)")
    idxs[max(1, end - 9):end]
else
    lo, hi = parse.(Int, split(steps_spec, ':'))
    collect(lo:hi)
end
println("Probing $(probe_source) run_name=$(probe_run_name) at steps " *
    "$(first(probe_steps))..$(last(probe_steps)) location=$(probe_loc)")

clamp_dir = if isempty(clamp_spec)
    nothing
else
    v = SVector{3}(parse.(Float64, split(clamp_spec, ','))...)
    v / norm(v)
end

# --- the monitor ---
# NOTE: probe/accumulator sizes are allocated LAZILY on the first call from
# the REPLAYED body's dimensions (2026-08-18 lesson: sizing from the
# setup-only rotor breaks when the setup env mesh differs from the saved
# run's mesh — BoundsError in job 13206092).
mutable struct TEDownwashProbe{TB}
    i_system::Int
    backend::TB
    loc::Symbol
    probes::Any
    nsta::Int          # stations per shedding (assumed equal across blades)
    nshedsets::Int
    sum_u_loc::Matrix{Float64}   # 3 x nsta, local basis, WAKE-only sources
    sum_u_glob::Matrix{Float64}  # 3 x nsta, global, WAKE-only sources
    sum_u_all::Matrix{Float64}   # 3 x nsta, local basis, ALL sources (diag:
                                 # includes the body's own bound near field,
                                 # which poisons a point sample — never the
                                 # table value)
    sum_r::Vector{Float64}       # station radius accum (for r/R)
    nsamples::Int                # (steps x blades) accumulated
end

function _tdp_ensure_sized!(m::TEDownwashProbe, body)
    nsta = size(body.Das[1], 2)
    nprobes = sum(size(body.Das[k], 2) for k in eachindex(body.Das))
    if m.nsta != nsta || m.probes === nothing
        m.nsamples == 0 || error("TE probe: replayed body station count " *
            "changed mid-replay ($(m.nsta) -> $(nsta))")
        m.nsta = nsta
        m.nshedsets = length(body.Das)
        m.probes = FastMultipole.ProbeSystem(nprobes, Float64)
        m.sum_u_loc = zeros(3, nsta)
        m.sum_u_glob = zeros(3, nsta)
        m.sum_u_all = zeros(3, nsta)
        m.sum_r = zeros(nsta)
    end
    return nothing
end

pnl.monitor_requires(::TEDownwashProbe) = ()
pnl.monitor_provides(::TEDownwashProbe) = ()

function (m::TEDownwashProbe)(systems, wakes, frames, uinf, i_step, dt)
    systems_tuple = systems isa Tuple ? systems : (systems,)
    body = systems_tuple[m.i_system]
    _tdp_ensure_sized!(m, body)
    zero_v = zero(SVector{3, Float64})
    zero_h = zero(SMatrix{3, 3, Float64, 9})

    # positions
    idx = 0
    for (k, shed) in enumerate(body.shedding)
        nshed = size(shed, 2)
        Das = body.Das[k]
        for j in 1:(nshed + 1)
            node_idx = j <= nshed ? body.cells[shed[3, j], shed[1, j]] :
                                    body.cells[shed[2, nshed], shed[1, nshed]]
            p = SVector{3}(body.nodes[1, node_idx], body.nodes[2, node_idx],
                           body.nodes[3, node_idx])
            if m.loc == :mid
                p += 0.5 * SVector{3}(Das[1, j], Das[2, j], Das[3, j])
            end
            idx += 1
            m.probes.position[idx] = p
            m.probes.gradient[idx] = zero_v
            m.probes.scalar_potential[idx] = 0.0
            m.probes.hessian[idx] = zero_h
        end
    end

    # PRIMARY: induced velocity from the WAKE only. The convection drift that
    # shapes the frozen path is the wake+freestream field; a point probe at
    # or near the TE would otherwise be dominated by the body's own bound-
    # sheet near field (measured at 0.5-0.9 tip speed in the shakeout — an
    # unrepresentative point sample of a field that averages out along the
    # path). The all-sources value is recorded as a diagnostic column.
    wake_sources = pnl._collect_wake_sources(wakes)
    grad_wake = Vector{SVector{3, Float64}}(undef, idx)
    if isempty(wake_sources)
        fill!(grad_wake, zero_v)
    else
        pnl.influence!((m.probes,), wake_sources, m.backend;
            precalc=false, scalar_potential=false, gradient=true, hessian=(false,))
        for i in 1:idx
            grad_wake[i] = m.probes.gradient[i]
            m.probes.gradient[i] = zero_v
        end
    end

    # DIAGNOSTIC: all sources (kernel-offset guard per DragPolarMonitor)
    all_sources = (systems_tuple..., wake_sources...)
    old_offsets = [sys.kerneloffset for sys in systems_tuple]
    try
        for (i, sys) in pairs(systems_tuple)
            sys.kerneloffset = i == m.i_system ? sys.kerneloffset_panel :
                                                 sys.kerneloffset_targets
        end
        pnl.influence!((m.probes,), all_sources, m.backend;
            precalc=false, scalar_potential=false, gradient=true, hessian=(false,))
    finally
        for (sys, offset) in zip(systems_tuple, old_offsets)
            sys.kerneloffset = offset
        end
    end

    # accumulate in the local TE basis
    uinf_sv = SVector{3, Float64}(uinf[1], uinf[2], uinf[3])
    axis_hat = frames[1].ω_axis / norm(frames[1].ω_axis)
    origin = frames[1].x
    idx = 0
    for (k, shed) in enumerate(body.shedding)
        nshed = size(shed, 2)
        for j in 1:(nshed + 1)
            node_idx = j <= nshed ? body.cells[shed[3, j], shed[1, j]] :
                                    body.cells[shed[2, nshed], shed[1, nshed]]
            p = SVector{3}(body.nodes[1, node_idx], body.nodes[2, node_idx],
                           body.nodes[3, node_idx])
            idx += 1
            u = grad_wake[idx] + uinf_sv
            u_all = m.probes.gradient[idx] + uinf_sv
            d = p - origin
            d_perp = d - axis_hat * dot(axis_hat, d)
            rj = norm(d_perp)
            rj > 0 || continue
            e_r = d_perp / rj
            vte = -cross(axis_hat, d_perp)
            e_t = vte / norm(vte)
            e_n = cross(e_t, e_r)
            m.sum_u_loc[1, j] += dot(u, e_t)
            m.sum_u_loc[2, j] += dot(u, e_r)
            m.sum_u_loc[3, j] += dot(u, e_n)
            m.sum_u_glob[:, j] .+= u
            m.sum_u_all[1, j] += dot(u_all, e_t)
            m.sum_u_all[2, j] += dot(u_all, e_r)
            m.sum_u_all[3, j] += dot(u_all, e_n)
            m.sum_r[j] += rj
        end
        m.nsamples += 1
    end
    return nothing
end

# sizes allocated lazily from the REPLAYED body on first call
probe_monitor = TEDownwashProbe(1, backend_wake, probe_loc,
    nothing, 0, 0, zeros(3, 0), zeros(3, 0), zeros(3, 0), zeros(0), 0)

result = pnl.replay(source_path, probe_run_name;
    monitors     = (probe_monitor,),
    Uinf         = Uinf,
    backend      = backend,
    backend_wake = backend_wake,
    steps        = probe_steps,
    verbose      = true,
)

# --- write the table ---
n = probe_monitor.nsamples
n > 0 || error("no samples accumulated")
nsta = probe_monitor.nsta
u_loc = probe_monitor.sum_u_loc ./ n
u_glob = probe_monitor.sum_u_glob ./ n
u_all = probe_monitor.sum_u_all ./ n
r_avg = probe_monitor.sum_r ./ n

# diagnostic: data-derived mean downstream direction (mid-span stations),
# REPORT-ONLY (Ryan 2026-08-18: the clamp is user-specified, never derived)
mid = [j for j in 1:nsta if 0.4 <= r_avg[j] / R <= 0.8]
u_mid = isempty(mid) ? zeros(3) : vec(sum(u_glob[:, mid]; dims=2)) ./ length(mid)
downstream = norm(u_mid) > 0 ? u_mid ./ norm(u_mid) : u_mid

isdir(dirname(out_csv)) || mkpath(dirname(out_csv))
open(out_csv, "w") do io
    println(io, "# TE downwash probe: source=$(probe_source) steps=" *
        "$(first(probe_steps)):$(last(probe_steps)) location=$(probe_loc) " *
        "samples=$(n) (steps x blades)")
    println(io, "# basis: u_t along kinematic TE tangent, u_r radial/spanwise, " *
        "u_n = e_t x e_r; u = WAKE-induced + Uinf inertial drift (body's own " *
        "bound field EXCLUDED — poisons a point sample; uall_* diagnostic " *
        "columns include it)")
    println(io, "# diagnostic mean downstream dir (mid-span, global): " *
        "$(round.(downstream, sigdigits=4))")
    clamp_dir === nothing ||
        println(io, "# clamp dir (global): $(round.(clamp_dir, sigdigits=4))")
    # The clamp operates on the LOCAL components (that is what the driver
    # consumes). The local n-axis is the axial direction (e_n = e_t x e_r,
    # sign set by the rotation sense), so a clamp direction must be
    # (anti)parallel to the rotation axis for the mapping to be exact —
    # error out otherwise rather than clamp approximately.
    n_sign = 0.0
    if clamp_dir !== nothing
        axis_hat0 = frames[1].ω_axis / norm(frames[1].ω_axis)
        a = dot(clamp_dir, axis_hat0)
        abs(abs(a) - 1) < 1e-6 || error(
            "P018_PROBE_CLAMP_DIR must be (anti)parallel to the rotation " *
            "axis (got |cos| = $(abs(a))); the local-basis table only maps " *
            "an axial clamp exactly")
        # sign of e_n along +axis: e_n = e_t x e_r; for this construction
        # e_n = -axis_hat * sign convention resolved empirically from the
        # mid-span diagnostic: use the measured mapping.
        j0 = isempty(mid) ? 1 : mid[1]
        en_along_axis = abs(u_loc[3, j0]) > 0 ?
            sign(dot(SVector{3}(u_glob[1, j0], u_glob[2, j0], u_glob[3, j0]),
                     axis_hat0) / u_loc[3, j0]) : 1.0
        # clamp fires when the GLOBAL drift has positive projection on
        # clamp_dir <=> local u_n * en_along_axis * a > 0
        n_sign = en_along_axis * a
    end
    println(io, "r_over_R,u_t,u_r,u_n,u_mag,u_gx,u_gy,u_gz,clamped," *
        "uall_t,uall_r,uall_n")
    for j in 1:nsta
        ut, ur, un = u_loc[1, j], u_loc[2, j], u_loc[3, j]
        clamped = 0
        if clamp_dir !== nothing && un * n_sign > 0
            clamped = 1
            un = 0.0
        end
        println(io, join((round(r_avg[j] / R, digits=6),
            round(ut, sigdigits=8), round(ur, sigdigits=8),
            round(un, sigdigits=8),
            round(sqrt(ut^2 + ur^2 + un^2), sigdigits=8),
            round(u_glob[1, j], sigdigits=8), round(u_glob[2, j], sigdigits=8),
            round(u_glob[3, j], sigdigits=8), clamped,
            round(u_all[1, j], sigdigits=8), round(u_all[2, j], sigdigits=8),
            round(u_all[3, j], sigdigits=8)), ','))
    end
end
println("Wrote $(out_csv) ($(nsta) stations, $(n) samples)")
for j in 1:nsta
    umag = norm(u_loc[:, j])
    println("  r/R $(round(r_avg[j]/R, digits=3)): u_t $(round(u_loc[1,j], sigdigits=3)) " *
        "u_r $(round(u_loc[2,j], sigdigits=3)) u_n $(round(u_loc[3,j], sigdigits=3)) " *
        "|u| $(round(umag, sigdigits=3)) m/s " *
        "($(round(umag/(ω_full*R), sigdigits=2)) tip) " *
        "[all-src |u| $(round(norm(u_all[:,j]), sigdigits=3))]")
end
