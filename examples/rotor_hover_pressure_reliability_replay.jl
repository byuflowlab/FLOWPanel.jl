#=##############################################################################
# DESCRIPTION (T2+T3 of plans/20260709_pressure_monitor_reliability.md)
#   Lamb-vector cases in this audit are deprecated diagnostic channels only;
#   they are not validation candidates or production pressure methods.
#   Single replay pass over the settled 2.0R rotor-hover cycle that recomputes
#   the pressure monitors with lamb-vector ingredient swaps, a projected-du
#   material-derivative variant, and a steady/unsteady Bernoulli pair, so the
#   0.0665 (lamb) <-> 0.0713 (Du/Dt) <-> 0.0685 (Bernoulli) CT spread can be
#   attributed term by term:
#     laplace_lamb           baseline (omega = induced_vorticity = wake + kappa)
#     laplace_lamb_nobound   omega = induced - kappa   (wake-only vorticity)
#     laplace_lamb_boundonly omega = kappa             (no wake vorticity)
#     laplace_lamb_hesscurl  omega = Hessian curl + kappa (T4 predicts == baseline
#                            up to FMM error on the extra outputs)
#     laplace_matderiv       baseline material-derivative (raw edge du)
#     laplace_matderiv_proj  du projected off the edge-averaged normal (T1
#                            predicts this collapses matderiv onto lamb modulo
#                            the omega injection)
#     bernoulli              steady twin (matches the original runs)
#     bernoulli_unsteady     + rho*dphi/dt (T3: rotating-frame unsteadiness;
#                            first step carries the FD transient - drop it)
#
#   All PressureLaplace variants use unsteady=false (as the originals), rho =
#   1.179, gauge at panel 1. Loaded body velocity/potential from the VTUs are
#   reused; induced_vorticity and the velocity Hessian are recomputed per step
#   (they are not saved).
#
# ENV knobs:
#   RELIABILITY_MODE  "gate" (default: 3 steps, check baselines vs
#                     data/rotor_hover_replay2p0_forces/replay_CT_vs_rev.csv),
#                     "full" (steps 1152:1223 + CSV + cycle-mean table)
#   RELIABILITY_STEPS "a:b" override of the replayed steps
#   RELIABILITY_OUT   output dir (default data/rotor_hover_pressure_reliability)
#
# Usage: julia -t auto --project examples/rotor_hover_pressure_reliability_replay.jl
=###############################################################################

import FLOWPanel as pnl
using Printf

const ROTOR_R = 0.119
const RHO = 1.179
const NT_PER_REV = 36
const AXIAL_DIMENSION = 1
const EXPERIMENT_CT = 0.072

const DATA_PATH = joinpath("data", "rotor_hover_relaxfilter2p0_ws")
const RUN_NAME = "rotor_hover_relaxfilter2p0_ws"
const REF_CSV = joinpath("data", "rotor_hover_replay2p0_forces", "replay_CT_vs_rev.csv")

mode = get(ENV, "RELIABILITY_MODE", "gate")
out_dir = get(ENV, "RELIABILITY_OUT", joinpath("data", "rotor_hover_pressure_reliability"))
steps_str = get(ENV, "RELIABILITY_STEPS", mode == "gate" ? "1152:1154" : "1152:1223")
lo, hi = parse.(Int, split(steps_str, ":"))
selected = collect(lo:hi)

# ---------------- monitor tuple -------------------------------------------------

const VARIANTS = (
    :laplace_lamb,
    :laplace_lamb_nobound,
    :laplace_lamb_boundonly,
    :laplace_lamb_hesscurl,
    :laplace_matderiv,
    :laplace_matderiv_proj,
    :bernoulli,
    :bernoulli_unsteady,
    :laplace_lamb_trikappa,   # must stay last: KappaSwap mutates induced_vorticity
)

# Swap the quad-basis kappa inside body.induced_vorticity (simulate!/replay
# default) for the tri-basis kappa, to measure the lamb form's sensitivity to
# the bound-vorticity discretization. Mutates the body, so this variant runs
# after all others.
struct KappaSwap end
pnl.monitor_provides(::KappaSwap) = ()
function pnl._run_monitor!(::KappaSwap, ctx::pnl.MonitorContext, systems, wakes,
        frames::AbstractVector{<:pnl.ReferenceFrame}, uinf, i_step::Int, dt::Real,
        t=nothing)
    for body in systems
        pnl._subtract_bound_surface_vorticity!(body.induced_vorticity, body;
            grad_mu_options=(; basis=:quad))
        pnl._add_bound_surface_vorticity_into!(body.induced_vorticity, body;
            grad_mu_options=(; basis=:tri))
    end
    return nothing
end

function make_monitors(systems, nt)
    normalization() = pnl.RotorNormalization(RHO, 2 * ROTOR_R, 1)
    fm() = pnl.ForceMonitor(nt, 1; i_frame=1, normalization=normalization(),
        correct_kuttacondition=false, verbose=false, file=false, vtk_fields=())
    pl(; kwargs...) = pnl.PressureLaplace(systems, RHO;
        unsteady=false, reference_panel=1, reference_pressure=0.0,
        verbose=false, file=false, vtk_fields=(), kwargs...)
    bern(unsteady) = pnl.PressureBernoulli(RHO;
        unsteady, allow_partial=unsteady, correct_kuttacondition=false,
        file=false, vtk_fields=())

    pressure_monitors = (
        pl(; acceleration_form=:lamb_vector),
        pl(; acceleration_form=:lamb_vector, lamb_vorticity=:no_bound),
        pl(; acceleration_form=:lamb_vector, lamb_vorticity=:bound_only),
        pl(; acceleration_form=:lamb_vector, lamb_vorticity=:hessian_curl),
        pl(; acceleration_form=:material_derivative),
        pl(; acceleration_form=:material_derivative, project_edge_du=true),
        bern(false),
        bern(true),
    )
    monitors = Any[]
    force_monitors = Any[]
    for pm in pressure_monitors
        f = fm()
        push!(monitors, pm, f)
        push!(force_monitors, f)
    end
    # tri-kappa lamb variant, isolated at the end (KappaSwap mutates state)
    f_tri = fm()
    push!(monitors, KappaSwap(), pl(; acceleration_form=:lamb_vector), f_tri)
    push!(force_monitors, f_tri)
    return Tuple(monitors), force_monitors
end

# ---------------- replay pass ----------------------------------------------------

println("rotor_hover_pressure_reliability_replay: mode=$(mode), steps $(lo):$(hi) ($(length(selected)) steps)")
force_monitors = nothing
factory = (systems, wakes, frames, t_range) -> begin
    monitors, force_monitors_local = make_monitors(systems, length(t_range))
    global force_monitors = force_monitors_local
    return monitors
end
elapsed = @elapsed pnl.replay(DATA_PATH, RUN_NAME;
    monitor_factory=factory, steps=selected, recompute=(:auto,), verbose=true)
@printf("replay done in %.1f s (%.2f s/step)\n", elapsed, elapsed / length(selected))

cts = Dict(v => -force_monitors[i].force[AXIAL_DIMENSION, :]
           for (i, v) in enumerate(VARIANTS))

# ---------------- gate: baselines must reproduce the saved-pressure replay -------

function reference_cts()
    isfile(REF_CSV) || return nothing
    lines = readlines(REF_CSV)
    header = Symbol.(split(lines[1], ","))
    rows = Dict{Int, Dict{Symbol, Float64}}()
    for line in lines[2:end]
        vals = parse.(Float64, split(line, ","))
        rows[Int(vals[1])] = Dict(zip(header, vals))
    end
    return rows
end

function run_gate(cts, selected)
    ref = reference_cts()
    if ref === nothing
        println("(no reference CSV at $(REF_CSV); gate skipped)")
        return
    end
    println("\n=== gate: recomputed baselines vs saved-pressure replay (rtol tolerance 1e-3) ===")
    pairs = ((:laplace_lamb, :CT_laplace_lamb),
             (:laplace_matderiv, :CT_laplace_matderiv),
             (:bernoulli, :CT_bernoulli))
    worst = 0.0
    for (k, step) in enumerate(selected)
        haskey(ref, step) || continue
        for (v, refcol) in pairs
            got = cts[v][k]
            want = ref[step][refcol]
            dev = abs(got - want) / max(abs(want), 1e-12)
            worst = max(worst, dev)
            k <= 3 && @printf("  step %d %-18s recomputed=%.6f saved=%.6f rel.dev=%.2e\n",
                              step, String(v), got, want, dev)
        end
    end
    @printf("  worst relative deviation over %d steps: %.3e  %s\n",
            length(selected), worst, worst < 1e-3 ? "PASS" : "FAIL (diagnose before full run)")
end

run_gate(cts, selected)

# ---------------- outputs ---------------------------------------------------------

mkpath(out_dir)
csv_path = joinpath(out_dir, "T2_T3_variant_CT_vs_step.csv")
open(csv_path, "w") do io
    println(io, "step,revolution," * join("CT_" .* String.(VARIANTS), ","))
    for (k, step) in enumerate(selected)
        println(io, "$step,$(step / NT_PER_REV)," *
                    join((@sprintf("%.8f", cts[v][k]) for v in VARIANTS), ","))
    end
end
println("Wrote $(csv_path)")

println("\n=== cycle statistics over steps $(lo):$(hi) (bernoulli_unsteady drops step 1: FD transient) ===")
@printf("%-26s %10s %10s %12s\n", "variant", "mean CT", "std", "gap to $(EXPERIMENT_CT)")
for v in VARIANTS
    vals = cts[v]
    v == :bernoulli_unsteady && (vals = vals[2:end])
    m = sum(vals) / length(vals)
    s = sqrt(sum((x - m)^2 for x in vals) / max(length(vals) - 1, 1))
    @printf("%-26s %10.5f %10.5f %+12.5f\n", String(v), m, s, m - EXPERIMENT_CT)
end
