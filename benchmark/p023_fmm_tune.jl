#!/usr/bin/env julia

# BRAINSTORM/023 stage: FMM parameter tuning on the MATURE production state
# (p018_cs_f1_l3p4 @ step ~1034, ~181k particles), per Ryan 2026-08-20:
#   1) FastMultipole.tune_fmm for ~1e-4 RELATIVE velocity error due to the
#      wake pass and ~1e-5 relative error due to the panel (body) pass;
#   2) scripted perturbation descent from the tuned point until the measured
#      cost reaches a minimum (FastMultipole.tune_fmm_perturb, added to
#      FastMultipole/src/autotune.jl for reuse);
#   3) verification: cost AND achieved relative error of production knobs,
#      the tune_fmm point, and the perturbed minimum, all measured through
#      the PRODUCTION influence path (_sa_wake_influence!/_sa_body_influence!)
#      against a high-accuracy FMM reference.
#
# Reuses the restore/setup machinery of p018_mature_wake_timing.jl (included
# below; its main() does not fire on include). Run via
# benchmark/slurm/p023_fmm_tune.sh with the verbatim production env.
#
# Env knobs: TUNE_WAKE_REL (1e-4), TUNE_BODY_REL (1e-5), TUNE_REPS (2),
# TUNE_MAX_ITERS (8), TUNE_REF_P (16), TUNE_REF_MAC (0.3),
# TUNE_PASSES ("wake,body").

import FLOWPanel as pnl
import FastMultipole
using Printf

include(joinpath(@__DIR__, "p018_mature_wake_timing.jl"))

const WAKE_REL = parse(Float64, get(ENV, "TUNE_WAKE_REL", "1e-4"))
const BODY_REL = parse(Float64, get(ENV, "TUNE_BODY_REL", "1e-5"))
const TUNE_REPS = parse(Int, get(ENV, "TUNE_REPS", "2"))
const TUNE_MAX_ITERS = parse(Int, get(ENV, "TUNE_MAX_ITERS", "8"))
const REF_P = parse(Int, get(ENV, "TUNE_REF_P", "16"))
const REF_MAC = parse(Float64, get(ENV, "TUNE_REF_MAC", "0.3"))
const PASSES = split(get(ENV, "TUNE_PASSES", "wake,body"), ",")
const TUNE_CSV = joinpath(RESULTS_DIR, "p023_fmm_tune.csv")
const HIST_CSV = joinpath(RESULTS_DIR, "p023_fmm_tune_history.csv")

function tune_main()
    threads = Threads.nthreads()
    println("BRAINSTORM/023 FMM tuning — restart run $(RESTART_RUN), threads $(threads)")
    println("  targets: wake rel $(WAKE_REL), body rel $(BODY_REL); reps $(TUNE_REPS), max_iters $(TUNE_MAX_ITERS)")
    flush(stdout)

    # TUNE_RADIUS_TOL: override pnl.FMM_RADIUS_TOL[] (default 1e-6) — the
    # radius-inflation tolerance driving Δr = core_size*(2/tol)^(1/4) for
    # VortexRing-bearing bodies. 1e-6 is 10-100x tighter than the pass
    # certification targets; this measures what a pass-consistent tolerance
    # buys (Ryan 2026-08-20, "75x seems excessive"). Achieved error is still
    # certified against the direct reference below.
    if haskey(ENV, "TUNE_RADIUS_TOL")
        pnl.FMM_RADIUS_TOL[] = parse(Float64, ENV["TUNE_RADIUS_TOL"])
        println("  FMM_RADIUS_TOL overridden to $(pnl.FMM_RADIUS_TOL[])")
    end

    cfg = canonical_setup!()
    restart_step = parse(Int, get(ENV, "RESTART_STEP", "-1"))
    restart_step < 0 && (restart_step = Base.invokelatest(default_restart_step))
    systems_tuple, wakes_tuple, _ =
        Base.invokelatest(warmstart_restore!, cfg, restart_step)
    println("  restored $(particle_count(wakes_tuple)) particles at step $(restart_step)")
    flush(stdout)

    # production per-step state the influence phases see: monitors demand no
    # hessian/vorticity in this config; core_sizes in target mode
    needs_grad = any(pnl.monitor_requires_body_hessian, cfg.monitors)
    for sys in systems_tuple
        sys.needs_velocity_gradient[] = needs_grad
    end
    pnl._set_core_sizes!(systems_tuple, :core_size_targets)

    FV_PF = pnl.FLOWVPM.ParticleField
    wake_probes, targets, wake_sources = pnl._sa_collect(systems_tuple, wakes_tuple)
    pfield_sources = Tuple(s for s in wake_sources if s isa FV_PF)
    pfield = only(pfield_sources)

    # production derivative requests (mirror _sa_wake_influence!/_sa_body_influence!)
    wake_hessian = Tuple(pnl.requires_hessian(sys) for sys in targets)
    body_hessian = Tuple((sys isa FV_PF && !cfg.body_hessian_to_particles) ?
        false : pnl.requires_hessian(sys) for sys in targets)
    conditioning = pnl._self_panel_core_size_conditioning()

    # --- production-path phase evaluators (cost + target-velocity capture) ---
    # Each evaluation replays the production step head: maneuver!-consistent
    # kinematics were set once below; reset + freestream + kinematic velocity
    # + update_TE! mirror _steady_aerodynamics! exactly. update_TE! is
    # REQUIRED: after the restore's replayed shed_wake!, the first wake row is
    # stale/sentinel until update_TE! rebuilds it from the body TE + Das —
    # skipping it NaNs every target through the panel-wake sources (root cause
    # of jobs 13245865/91/13246014's all-NaN fields; localized by local
    # bisection 2026-08-20).
    t_step = cfg.t_range[restart_step + 2]
    Base.invokelatest(cfg.maneuver, cfg.frames, systems_tuple, wakes_tuple, t_step)
    uinf = Base.invokelatest(cfg.Uinf, t_step)
    function zero_targets!()
        pnl._sa_reset_freestream_kinematic!(systems_tuple, wakes_tuple, cfg.frames, uinf)
        for (sys, w) in zip(systems_tuple, wakes_tuple)
            !isnothing(w) && pnl.update_TE!(w, sys)
        end
    end
    function wake_phase!(backend)
        zero_targets!()
        pnl._sa_wake_influence!(targets, wake_sources, backend;
            needs_induced_vorticity=false,
            wakerow_no_hessian_to_particles=cfg.wakerow_no_hessian_to_particles,
            panel_wake_on_particles=cfg.panel_wake_on_particles,
            particle_hessian_self=cfg.particle_hessian_self)
    end
    function body_phase!(backend)
        zero_targets!()
        pnl._sa_body_influence!(targets, systems_tuple, backend;
            needs_induced_vorticity=false,
            body_on_wake=cfg.body_on_wake,
            body_hessian_to_particles=cfg.body_hessian_to_particles,
            body_gradient_core_size=cfg.body_gradient_core_size)
    end
    function snapshot()
        body_v = reduce(hcat, (copy(sys.velocity) for sys in systems_tuple))
        np = pfield.np
        U = Matrix{Float64}(undef, 3, np)
        for i in 1:np
            U[:, i] .= pnl.FLOWVPM.get_U(pfield, i)
        end
        return hcat(body_v, U)
    end
    # The prelude writes -V_kinematic (rotation apparent velocity, ~tip speed)
    # into body.velocity; subtract that deterministic baseline so error
    # metrics act on the WAKE-INDUCED field only ("error due to the wake",
    # Ryan 2026-08-20), not a denominator inflated by rigid-body motion.
    zero_targets!()
    base0 = snapshot()
    wake_field() = snapshot() .- base0
    relerr(ref, v, mask) = sqrt(sum(abs2, view(v, :, mask) .- view(ref, :, mask)) /
        sum(abs2, view(ref, :, mask)))
    function measure(phase!, backend, ref, mask)
        t = Inf
        for _ in 1:TUNE_REPS
            t = min(t, elapsed_s(() -> phase!(backend)))
        end
        v = wake_field()
        new_bad = count(j -> mask[j] && !all(isfinite, view(v, :, j)), axes(v, 2))
        new_bad > 0 && println("   WARNING: candidate produced $(new_bad) NEW non-finite columns")
        return t, relerr(ref, v, mask), v
    end
    # Diagnose and mask non-finite target columns (job 13245957: the reference
    # snapshot contained NaN, poisoning both FastMultipole's relative error
    # method and the naive RMS; localize them, then compute all metrics over
    # the finite columns only)
    function finite_mask(ref)
        mask = vec(all(isfinite, ref; dims=1))
        n_bad = count(!, mask)
        n_body = size(ref, 2) - pfield.np
        if n_bad > 0
            bad = findall(!, mask)
            bad_body = count(<=(n_body), bad)
            println("   WARNING: $(n_bad) non-finite target columns " *
                "($(bad_body) body CPs, $(n_bad - bad_body) particles); masking them")
            for j in first(bad, 5)
                if j <= n_body
                    println("     body CP column $(j)")
                else
                    i = j - n_body
                    try
                        P = pnl.FLOWVPM.get_particle(pfield, i)
                        println("     particle $(i): X=$(P.X) sigma=$(P.sigma) " *
                            "|Gamma|=$(sqrt(sum(abs2, P.Gamma)))")
                    catch
                        println("     particle $(i)")
                    end
                end
            end
        end
        return mask
    end

    io = open(TUNE_CSV, "w")
    println(io, "pass,variant,expansion_order,multipole_acceptance,leaf_size," *
        "phase_t_s,rel_err,meets_target,threads,n_particles,restart_step")
    io_h = open(HIST_CSV, "w")
    println(io_h, "pass,iter,expansion_order,multipole_acceptance,leaf_size,trial_t_s,error_success,accepted")

    for pass in PASSES
        pass = strip(pass)
        is_wake = pass == "wake"
        tol = is_wake ? WAKE_REL : BODY_REL
        phase! = is_wake ? wake_phase! : body_phase!
        tune_targets = targets
        tune_sources = is_wake ? pfield_sources : systems_tuple
        hessian = is_wake ? wake_hessian : body_hessian
        fmm_kwargs = is_wake ?
            (; scalar_potential=false, gradient=true, hessian) :
            (; scalar_potential=false, gradient=true, hessian,
               direct_conditioning=conditioning)
        prod_backend = is_wake ?
            pnl.FastMultipoleBackend(4, 0.4, 50) :   # RHPC wake defaults
            pnl.FastMultipoleBackend(8, 0.4, 20)     # RHPC body defaults

        println("\n#### pass=$(pass): tolerance $(tol) (relative velocity) ####")
        flush(stdout)

        # Reference for achieved-error verification: DirectBackend = exact.
        # (A p=16/mac=0.3 FMM reference produced a FULLY non-finite target
        # field — all 218,059 columns NaN, job 13245991 — while p=4 production
        # knobs are clean across four jobs; high expansion order at this scale
        # is itself broken. Documented via the fmm_p16_diag row below.)
        println("-- reference evaluation (DirectBackend, exact) --")
        t_ref = elapsed_s(() -> phase!(pnl.DirectBackend()))
        ref = wake_field()
        println("   reference took $(round(t_ref; digits=1)) s")
        flush(stdout)

        # Relative target implemented as absolute = rel x RMS reference
        # velocity (021 convention). The cluster FastMultipole's
        # PowerRelativeGradient path evaluates its tolerance to NaN in
        # translate.jl (observed job 13245865), spuriously rejecting every
        # MAC, so the relative error method is unusable at this version.
        mask = finite_mask(ref)
        rms_ref = sqrt(sum(abs2, view(ref, :, mask)) / count(mask))
        error_tolerance = FastMultipole.PowerAbsoluteGradient(tol * rms_ref)
        println("   RMS reference velocity $(rms_ref) -> absolute tolerance $(tol * rms_ref)")

        # TUNE_KNOWN_WAKE / TUNE_KNOWN_BODY = "p,mac,leaf": reuse an
        # already-tuned point (skip stages 1-2, which cost ~31 min on the
        # body pass) and go straight to verification.
        known = get(ENV, is_wake ? "TUNE_KNOWN_WAKE" : "TUNE_KNOWN_BODY", "")
        local p1, mac1, leaf1, p2, mac2, leaf2
        if !isempty(known)
            ps, ms, ls = split(known, r"[:,]")
            p1 = p2 = parse(Int, ps)
            mac1 = mac2 = parse(Float64, ms)
            leaf1 = leaf2 = parse(Int, ls)
            println("known tuned point from ENV: p=$(p1) mac=$(mac1) leaf=$(leaf1); skipping tune+perturb")
        else

        # 1) tune_fmm at the target tolerance (older cluster tune_fmm returns
        # a 2-tuple; index rather than destructure)
        t_tune = @elapsed tune_result = FastMultipole.tune_fmm(
            tune_targets, tune_sources;
            error_tolerance,
            multipole_acceptances=range(0.3, 0.8, step=0.1),
            max_expansion_order=20,
            verbose=true, fmm_kwargs...)
        tuned = tune_result[1]
        p1, mac1 = tuned.expansion_order, tuned.multipole_acceptance
        leaf1 = tuned.leaf_size_source[1]
        println("tune_fmm ($(round(t_tune; digits=1)) s): p=$(p1) mac=$(mac1) leaf=$(leaf1)")
        flush(stdout)

        # 2) scripted perturbation descent until the cost minimum
        # tree_amortization is left at its default 1 DELIBERATELY: 023 prices a
        # mature-wake unsteady step, whose geometry moves every step, so the
        # tree/interaction-list build IS paid per apply and belongs in the
        # objective. A steady iterative solve reusing one FmmPlan should pass
        # tree_amortization=<iteration count> instead (see BRAINSTORM 021,
        # 2026-08-24, and the tune_fmm_perturb docstring).
        (best, history) = FastMultipole.tune_fmm_perturb(
            tune_targets, tune_sources;
            tuned...,
            error_tolerance,
            max_expansion_order=20,
            reps=TUNE_REPS, max_iters=TUNE_MAX_ITERS,
            verbose=true, fmm_kwargs...)
        p2, mac2 = best.expansion_order, best.multipole_acceptance
        leaf2 = best.leaf_size_source isa Integer ? best.leaf_size_source :
            best.leaf_size_source[1]
        for h in history
            lf = h.leaf_size_source isa Integer ? h.leaf_size_source : h.leaf_size_source[1]
            println(io_h, "$(pass),$(h.iter),$(h.expansion_order)," *
                "$(h.multipole_acceptance),$(lf),$(h.t),$(h.error_success),$(h.accepted)")
        end
        flush(stdout)

        end  # !isempty(known)

        # 3) verification through the production influence path.
        # TUNE_P_LADDER (comma list of fixed expansion orders) appends a
        # fixed-p ladder at the perturbed minimum's (MAC, leaf): the descent's
        # P-1 rejections come from FastMultipole's conservative worst-pair
        # dynamic-P model, while production uses FIXED p judged on measured
        # global error — lower fixed p may still meet the target (Ryan
        # 2026-08-20, "the expansion orders look too high").
        variants = Any[
            ("production", prod_backend),
            ("tuned", pnl.FastMultipoleBackend(p1, mac1, leaf1)),
            ("perturbed_min", pnl.FastMultipoleBackend(p2, mac2, leaf2)),
            ("fmm_p16_diag", pnl.FastMultipoleBackend(REF_P, REF_MAC, 50)),
        ]
        for pstr in split(get(ENV, "TUNE_P_LADDER", ""), r"[:,]"; keepempty=false)
            p = parse(Int, strip(pstr))
            push!(variants, ("pladder_p$(p)", pnl.FastMultipoleBackend(p, mac2, leaf2)))
        end
        for (variant, backend) in variants
            t, err, _ = measure(phase!, backend, ref, mask)
            meets = err <= tol
            @printf("  %-14s p=%2d mac=%.3f leaf=%4d: %8.2f s  rel_err=%.3e %s\n",
                variant, backend.expansion_order, backend.multipole_acceptance,
                backend.leaf_size, t, err, meets ? "" : "[misses target]")
            println(io, "$(pass),$(variant),$(backend.expansion_order)," *
                "$(backend.multipole_acceptance),$(backend.leaf_size),$(t),$(err)," *
                "$(meets),$(threads),$(pfield.np),$(restart_step)")
            flush(io); flush(stdout)
        end
    end
    close(io); close(io_h)
    println("\nWrote $(TUNE_CSV) and $(HIST_CSV)\nDone.")
end

if abspath(PROGRAM_FILE) == (@__FILE__) || isempty(PROGRAM_FILE)
    tune_main()
end
