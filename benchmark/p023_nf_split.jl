#!/usr/bin/env julia

# BRAINSTORM/023: measure the nearfield-vs-farfield split of the two dominant
# influence phases directly, by timing fmm! with its nearfield/farfield
# switches on the frozen mature state (p018_cs_f1_l3p4 @ ~step 1034).
# Motivation: the flat/tree profiles show nearfield direct kernels dominating,
# but one 244k-sample FastMultipole threaded block resists attribution;
# wall-clock ablation is unambiguous.
#
# For each of the two production calls (wake pfield -> all targets;
# body -> all targets), times min-of-REPS of: full, nearfield-only
# (farfield=false), farfield-only (nearfield=false). Also with/without
# postcalc on the pfield call to price the Dynamic-SFS estimator (Estr_fmm!)
# that rides on it.

import FLOWPanel as pnl
using Printf

include(joinpath(@__DIR__, "p018_mature_wake_timing.jl"))

const REPS = parse(Int, get(ENV, "SPLIT_REPS", "2"))

function split_main()
    println("BRAINSTORM/023 nearfield/farfield split — $(RESTART_RUN), $(Threads.nthreads()) threads, min-of-$(REPS)")
    cfg = canonical_setup!()
    restart_step = parse(Int, get(ENV, "RESTART_STEP", "-1"))
    restart_step < 0 && (restart_step = Base.invokelatest(default_restart_step))
    systems_tuple, wakes_tuple, _ = Base.invokelatest(warmstart_restore!, cfg, restart_step)
    println("  restored $(particle_count(wakes_tuple)) particles at step $(restart_step)")
    flush(stdout)

    needs_grad = any(pnl.monitor_requires_body_hessian, cfg.monitors)
    for sys in systems_tuple
        sys.needs_velocity_gradient[] = needs_grad
    end
    pnl._set_core_sizes!(systems_tuple, :core_size_targets)

    FV_PF = pnl.FLOWVPM.ParticleField
    _, targets, wake_sources = pnl._sa_collect(systems_tuple, wakes_tuple)
    pfield_sources = Tuple(s for s in wake_sources if s isa FV_PF)
    no_vort = pnl._induced_vorticity_extra_outputs(targets, false)
    wake_hessian = Tuple(pnl.requires_hessian(sys) for sys in targets)
    body_hessian = Tuple((sys isa FV_PF && !cfg.body_hessian_to_particles) ?
        false : pnl.requires_hessian(sys) for sys in targets)

    wake_call(; kwargs...) = pnl.influence!(targets, pfield_sources, cfg.backend_wake;
        precalc=true, postcalc=true, scalar_potential=false, velocity=true,
        velocity_gradient=wake_hessian, extra_outputs=no_vort, kwargs...)
    wake_call_nopost(; kwargs...) = pnl.influence!(targets, pfield_sources, cfg.backend_wake;
        precalc=true, postcalc=false, scalar_potential=false, velocity=true,
        velocity_gradient=wake_hessian, extra_outputs=no_vort, kwargs...)
    body_call(; kwargs...) = pnl.influence!(targets, systems_tuple, cfg.backend;
        precalc=false, scalar_potential=false, velocity=true,
        velocity_gradient=body_hessian, extra_outputs=no_vort,
        direct_conditioning=pnl._self_panel_core_size_conditioning(), kwargs...)

    function timeit(f; kwargs...)
        t = Inf
        for _ in 1:REPS
            t = min(t, elapsed_s(() -> f(; kwargs...)))
        end
        return t
    end

    println("\n| call | full [s] | nearfield-only [s] | farfield-only [s] | near+far−full [s] |")
    println("|---|---:|---:|---:|---:|")
    for (label, f) in (("wake pfield (with SFS postcalc)", wake_call),
                       ("wake pfield (no postcalc)", wake_call_nopost),
                       ("body -> all targets", body_call))
        t_full = timeit(f)
        t_near = timeit(f; farfield=false)
        t_far = timeit(f; nearfield=false)
        @printf("| %s | %.2f | %.2f | %.2f | %+.2f |\n",
            label, t_full, t_near, t_far, t_near + t_far - t_full)
        flush(stdout)
    end
    println("\nDone.")
end

if abspath(PROGRAM_FILE) == (@__FILE__) || isempty(PROGRAM_FILE)
    split_main()
end
