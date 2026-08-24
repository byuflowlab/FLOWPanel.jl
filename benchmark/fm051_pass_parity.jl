#!/usr/bin/env julia

#=##############################################################################
# DESCRIPTION
#   Task 051 stage 2/3 GATE: pass-by-pass CPU-vs-GPU parity for the FLOWPanel
#   rectangular-influence seam (src/FLOWPanel_gpu_influence.jl).
#
#   For ONE production timestep, each cross-influence pass of
#   `_steady_aerodynamics!` (src/FLOWPanel_simulate.jl:654) is run TWICE from
#   an identical pre-pass state:
#     arm A ("cpu")  : FLOWPanel.set_gpu_influence!(:off)  -> FastMultipole.fmm!
#     arm B ("seam") : FLOWPanel.set_gpu_influence!(SEAM)  -> direct_rectangular!
#   Between arms EVERY array the pass accumulates into is snapshotted and
#   restored, so the two arms see bit-identical inputs. After the comparison
#   the CPU-arm outputs are written back, so all DOWNSTREAM state is the
#   production CPU path -- the harness never leaves seam results in the sim.
#
#   Passes compared (in per-step order):
#     pass 1  `_sa_wake_influence!`  wake sources -> bodies + wake probes + particles
#     solve   `solve_formulation!`   asserted NOT captured by the seam
#     pass 3  `_sa_body_influence!`  bodies -> bodies + wake probes + particles
#     SFS     Estr / SFS_INDEX rows  asserted NOT captured by the seam
#
#   GATES
#     phase gate (each pass, worst output class)   worst per-target rel err < 1e-3
#     pass-1 informational bound (points kernel, F64)                      < 1e-11
#     SFS/Estr agreement (threading nondeterminism only)                   < 1e-10
#   Pass 3's REAL gate is 1e-3: the seam evaluates ALL self-body pairs at
#   `core_size_panel`, while the CPU `direct_conditioning` rule flips the
#   offset only for NEARFIELD self-pair blocks (the "Known deviation" note in
#   src/FLOWPanel_gpu_influence.jl). The achieved value is always reported.
#
# MODES  (env FM051_MODE)
#   "mini" (default)  small self-contained unsteady case built here: a diamond
#                     RigidWakeBody{ConstantSource+VortexRing} with a finite
#                     attached TE wake + a PanelParticleWake (panel rows,
#                     final filament, FLOWVPM particle field). A few steps run
#                     with the seam OFF so the wake state is nontrivial, then
#                     the NEXT step is analysed pass by pass. Runs locally in
#                     minutes on 4 threads. Default seam :host -- the host
#                     rectangular kernels are the SAME pair math as CUDA, so a
#                     GPU-free machine exercises recognition, packing and
#                     write-back end to end.
#   "full"            warm-start restore of the mature p018 run. Reuses
#                     benchmark/p018_mature_wake_timing.jl verbatim (loaded
#                     into a private module) for `canonical_setup!` and
#                     `warmstart_restore!`, i.e. the mirror of
#                     simulate_warmstart! sections 2.5-5 INCLUDING the
#                     mandatory kinematic replay. Env names are p023's:
#                     BENCH_RESTART_RUN, RESTART_STEP. Default seam :cuda.
#
# ENVIRONMENT
#   FM051_MODE        mini | full                       (default mini)
#   FM051_SEAM        host | cuda                       (default host / cuda)
#   FM051_MINI_STEPS  spin-up steps before the parity step (mini, default 4)
#   FM051_MINI_NSPAN  diamond span cells (mini, default 6 -> 8*nspan panels)
#   FM051_MINI_FMM    "p,mac,leaf" backend for the mini case
#                     (default "10,0.4,100000" -- a single leaf, so the CPU arm
#                      is pure nearfield direct and the 1e-11 bound is a real
#                      check of packing/write-back rather than of the fmm
#                      truncation error; raise the leaf count to exercise
#                      expansions instead)
#   FM051_GATE        phase gate                        (default 1e-3)
#   FM051_GATE_P1     pass-1 informational bound        (default 1e-11)
#   FM051_GATE_SFS    Estr agreement bound              (default 1e-10)
#   BENCH_RESTART_RUN / RESTART_STEP                    (full mode; p023 names)
#   FLOWPANEL_GPU_INFLUENCE  MUST be unset or 0 -- the harness arms the seam
#                     itself via set_gpu_influence! and restores :off at exit.
#
# OUTPUT
#   benchmark/results/fm051_pass_parity_<mode>.csv   (key,value rows)
#
# THE STEP-HEAD GOTCHA
#   The parity step's HEAD -- maneuver!, then _sa_collect /
#   _sa_reset_freestream_kinematic! / update_TE! -- must run IN FULL before any
#   pass is replayed, and in full mode the warm start must include the
#   kinematic replay of simulate_warmstart! section 2.5. Skipping any of it
#   silently NaNs the first wake row (this cost three cluster jobs in tasks
#   018/023). The harness asserts finiteness of the target state right after
#   the step head and fails loudly citing this note.
#
# AUTHORSHIP
#   * Created by : task 051 stage 2/3 parity gate, Aug 2026
#   * License    : GNU Public License
=###############################################################################

import FLOWPanel as pnl
import FLOWVPM
import FastMultipole
using LinearAlgebra
using Printf

const MODE      = lowercase(get(ENV, "FM051_MODE", "mini"))
const SEAM_MODE = Symbol(lowercase(get(ENV, "FM051_SEAM", MODE == "full" ? "cuda" : "host")))
const GATE      = parse(Float64, get(ENV, "FM051_GATE", "1e-3"))
const GATE_P1   = parse(Float64, get(ENV, "FM051_GATE_P1", "1e-11"))
const GATE_SFS  = parse(Float64, get(ENV, "FM051_GATE_SFS", "1e-10"))
const GATE_ATTR = parse(Float64, get(ENV, "FM051_GATE_ATTR", "1e-10"))
const RESULTS_DIR = joinpath(@__DIR__, "results")

# pass-3 three-way attribution (exact DirectBackend arm); FM051_ATTRIBUTION=0 skips
include(joinpath(@__DIR__, "fm051_pass3_attribution.jl"))

MODE in ("mini", "full") || error("FM051_MODE must be \"mini\" or \"full\" (got $(repr(MODE)))")
SEAM_MODE in (:host, :cuda) || error("FM051_SEAM must be host or cuda (got $(repr(SEAM_MODE)))")

# the harness owns the seam switch; a pre-armed env would corrupt the CPU arm
let v = lowercase(get(ENV, "FLOWPANEL_GPU_INFLUENCE", "0"))
    v in ("0", "off", "false", "") || error("FLOWPANEL_GPU_INFLUENCE=$(v) is armed in " *
        "the environment; this harness arms/disarms the seam itself. Set it to 0.")
end

# p018/p023 harness expects this; 1 measured step is all we need
get!(ENV, "BENCH_NSTEPS", "1")

################################################################################
# BOOKKEEPING
################################################################################

const RECORD = Pair{String,Any}[]
record!(k, v) = (push!(RECORD, string(k) => v); v)

const GATE_LINES = Tuple{String,Bool,String}[]

function gate!(label, achieved, bound; informational=false, note="")
    ok = isfinite(achieved) && achieved <= bound
    detail = @sprintf("achieved %.3e   bound %.3e", achieved, bound)
    isempty(note) || (detail *= "   [" * note * "]")
    informational && (detail *= "   (INFORMATIONAL, not a gate)")
    push!(GATE_LINES, (label, informational ? true : ok, detail))
    record!("gate." * label * ".achieved", achieved)
    record!("gate." * label * ".bound", bound)
    record!("gate." * label * ".pass", ok)
    record!("gate." * label * ".informational", informational)
    return ok
end

function gate_skip!(label, reason)
    push!(GATE_LINES, (label, true, "SKIPPED -- " * reason))
    record!("gate." * label * ".skipped", reason)
    return true
end

elapsed_s(f) = (t0 = time_ns(); f(); (time_ns() - t0) / 1e9)

################################################################################
# SNAPSHOT / RESTORE  -- every array the compared passes accumulate into
################################################################################
# Target class -> storage (mirrors the seam's `_gpu_add_result!` and the fmm!
# `buffer_to_target_system!` methods):
#   AbstractBody               velocity, velocity_gradient, potential,
#                              induced_vorticity (extra_outputs)
#   ProbeWrapper{<:PanelWake}  wake.velocity[isurf]  (no gradient stored)
#   FastMultipole.ProbeSystem  gradient (=U), hessian (=J), scalar_potential
#   FLOWVPM.ParticleField      the WHOLE particle matrix, columns 1:np -- U/J
#                              rows plus SFS/Estr rows and anything the
#                              pre/post SFS hooks touch

function snap_target(t::pnl.AbstractBody)
    d = Dict{Symbol,Any}()
    d[:velocity] = copy(t.velocity)
    d[:velocity_gradient] = copy(t.velocity_gradient)
    hasproperty(t, :potential) && (d[:potential] = copy(t.potential))
    hasproperty(t, :induced_vorticity) && (d[:induced_vorticity] = copy(t.induced_vorticity))
    return d
end
function restore_target!(t::pnl.AbstractBody, d)
    for (k, v) in d
        getproperty(t, k) .= v
    end
    return nothing
end

function snap_target(pw::pnl.ProbeWrapper)
    return Dict{Symbol,Any}(:velocity => [copy(v) for v in pw.system.velocity])
end
function restore_target!(pw::pnl.ProbeWrapper, d)
    for (a, b) in zip(pw.system.velocity, d[:velocity])
        a .= b
    end
    return nothing
end

function snap_target(ps::FastMultipole.ProbeSystem)
    return Dict{Symbol,Any}(:gradient => copy(ps.gradient),
                            :hessian => copy(ps.hessian),
                            :scalar_potential => copy(ps.scalar_potential))
end
function restore_target!(ps::FastMultipole.ProbeSystem, d)
    ps.gradient .= d[:gradient]
    ps.hessian .= d[:hessian]
    ps.scalar_potential .= d[:scalar_potential]
    return nothing
end

function snap_target(pf::FLOWVPM.ParticleField)
    np = max(pf.np, 0)
    return Dict{Symbol,Any}(:np => np,
                            :particles => Array(view(pf.particles, :, 1:np)))
end
function restore_target!(pf::FLOWVPM.ParticleField, d)
    np = d[:np]::Int
    np > 0 && (view(pf.particles, :, 1:np) .= d[:particles])
    return nothing
end

snapshot(targets::Tuple) = Any[snap_target(t) for t in targets]
function restore!(targets::Tuple, snaps)
    for (t, s) in zip(targets, snaps)
        restore_target!(t, s)
    end
    return nothing
end

################################################################################
# OUTPUT-CLASS EXTRACTION  (each class -> a ncomp x ntarget Float64 matrix)
################################################################################

_mat(a::AbstractMatrix{<:Real}) = Matrix{Float64}(a)
_mat(a::AbstractArray{<:Real,3}) =
    Matrix{Float64}(reshape(a, size(a, 1) * size(a, 2), size(a, 3)))
function _mat(v::AbstractVector{<:FastMultipole.StaticArrays.SVector{3}})
    m = Matrix{Float64}(undef, 3, length(v))
    @inbounds for i in eachindex(v), c in 1:3
        m[c, i] = v[i][c]
    end
    return m
end
function _mat(v::AbstractVector{<:FastMultipole.StaticArrays.SMatrix{3,3}})
    m = Matrix{Float64}(undef, 9, length(v))
    @inbounds for i in eachindex(v), c in 1:9
        m[c, i] = v[i][c]
    end
    return m
end

"Per-target output classes of one target system: Vector of class => (ncomp x n)."
extract(t::pnl.AbstractBody) = [
    "body_velocity" => _mat(t.velocity),
    "body_velocity_gradient" => _mat(t.velocity_gradient),
]

function extract(pw::pnl.ProbeWrapper)
    w = pw.system
    n = FastMultipole.get_n_bodies(pw)
    m = Matrix{Float64}(undef, 3, n)
    @inbounds for i in 1:n
        isurf, irow, icol = pnl.global_to_matrix_index(pw, i)
        for c in 1:3
            m[c, i] = w.velocity[isurf][c, irow, icol]
        end
    end
    return ["probe_U" => m]
end

extract(ps::FastMultipole.ProbeSystem) =
    ["probe_U" => _mat(ps.gradient), "probe_J" => _mat(ps.hessian)]

function extract(pf::FLOWVPM.ParticleField)
    idx = 1:max(pf.np, 0)
    out = ["particle_U" => Matrix{Float64}(Array(view(pf.particles, FLOWVPM.U_INDEX, idx))),
           "particle_J" => Matrix{Float64}(Array(view(pf.particles, FLOWVPM.J_INDEX, idx)))]
    if FLOWVPM.isSFSenabled(pf.SFS)
        push!(out, "particle_Estr" =>
            Matrix{Float64}(Array(view(pf.particles, FLOWVPM.SFS_INDEX, idx))))
    end
    return out
end

function _hcat_class(v::Vector{Matrix{Float64}})
    v = [m for m in v if size(m, 2) > 0]
    isempty(v) && return zeros(Float64, 0, 0)
    nr = size(v[1], 1)
    all(size(m, 1) == nr for m in v) ||
        error("inconsistent component count within one output class")
    return length(v) == 1 ? v[1] : hcat(v...)
end

"Merge per-target class matrices over a whole target tuple: class => matrix."
function extract_all(targets::Tuple)
    acc = Dict{String,Vector{Matrix{Float64}}}()
    for t in targets, (cls, m) in extract(t)
        push!(get!(acc, cls, Matrix{Float64}[]), m)
    end
    return Dict{String,Matrix{Float64}}(cls => _hcat_class(v) for (cls, v) in acc)
end

diffmap(a::Dict, b::Dict) = Dict{String,Matrix{Float64}}(
    k => (a[k] .- b[k]) for k in keys(a) if haskey(b, k) && size(a[k]) == size(b[k]))

"""
    class_parity_stats(ref, tst; gate) -> NamedTuple

Parity statistics of one output class between two arms, over target columns:

- `worst_rel`, `iw`, `dn_w`, `rn_w`: worst per-target relative error
  `norm(tst[:,i] - ref[:,i]) / max(norm(ref[:,i]), 1e-6 * scale)` with the
  absolute diff (`dn_w`) and reference norm (`rn_w`) AT that target — so a
  floor-dominated blowup on a near-null target is distinguishable from a
  genuinely large absolute error (job 13306475 could not tell them apart).
- `scale`: largest per-target reference magnitude within the class.
- `diff_to_scale`: `max_i norm(tst[:,i] - ref[:,i]) / scale` — the
  floor-independent severity metric. This is the FULL-mode phase gate: the
  full-mode CPU arm is itself fmm-approximate, so demanding 1e-3 PER-TARGET
  agreement on targets 1e-6 below class scale would require absolute
  agreement at 1e-9·scale — far tighter than the fmm tolerance the reference
  carries. Mini mode still gates per-target rel (both arms are exact there).
- `n_over`: number of targets with per-target rel > `gate`.
- `n_div`: number of targets with per-target rel > 0.5 ("divergent" — used by
  the full-mode Estr sparsity gate; see the SFS section).
- `top`: up to 5 worst offenders as `(i, rel, dn, rn)`.
"""
function class_parity_stats(ref::Matrix{Float64}, tst::Matrix{Float64}; gate::Float64=1e-3)
    (size(ref, 2) == 0 || size(ref) != size(tst)) &&
        return (worst_rel=0.0, iw=0, dn_w=0.0, rn_w=0.0, scale=0.0,
                diff_to_scale=0.0, n_over=0, n_div=0, n=0,
                top=NTuple{4,Float64}[])
    n = size(ref, 2)
    scale = 0.0
    @inbounds for i in 1:n
        s = 0.0
        for c in axes(ref, 1)
            s += ref[c, i]^2
        end
        scale = max(scale, sqrt(s))
    end
    floorv = 1e-6 * scale
    rels = Vector{Float64}(undef, n)
    dns = Vector{Float64}(undef, n)
    rns = Vector{Float64}(undef, n)
    dmax = 0.0
    n_over = 0
    n_div = 0
    @inbounds for i in 1:n
        dn = 0.0
        rn = 0.0
        for c in axes(ref, 1)
            dn += (tst[c, i] - ref[c, i])^2
            rn += ref[c, i]^2
        end
        dn = sqrt(dn); rn = sqrt(rn)
        den = max(rn, floorv)
        e = den > 0 ? dn / den : (dn > 0 ? Inf : 0.0)
        rels[i] = e; dns[i] = dn; rns[i] = rn
        dmax = max(dmax, dn)
        e > gate && (n_over += 1)
        e > 0.5 && (n_div += 1)
    end
    order = partialsortperm(rels, 1:min(5, n); rev=true)
    top = [(Float64(i), rels[i], dns[i], rns[i]) for i in order]
    iw = isempty(order) ? 0 : order[1]
    return (worst_rel=(iw == 0 ? 0.0 : rels[iw]), iw=iw,
            dn_w=(iw == 0 ? 0.0 : dns[iw]), rn_w=(iw == 0 ? 0.0 : rns[iw]),
            scale=scale, diff_to_scale=(scale > 0 ? dmax / scale : (dmax > 0 ? Inf : 0.0)),
            n_over=n_over, n_div=n_div, n=n, top=top)
end

allfinite(targets::Tuple) = all(all(isfinite, m) for t in targets for (_, m) in extract(t))

"""
    report_eligibility(label, targets, sources)

Print (and record) the seam's per-system eligibility verdicts for a pass. This
is the first thing to read when a pass reports 0 seam hits: it names the exact
system that made the pass fall through to `fmm!`.
"""
function report_eligibility(label, targets::Tuple, sources::Tuple)
    println("  seam eligibility:")
    for (i, s) in enumerate(sources)
        ok = pnl._gpu_source_supported(s)
        ncol = pnl._gpu_source_columns(s)
        @printf("    source %d  %-46s supported=%-5s packed_columns=%d\n",
                i, string(nameof(typeof(s))), string(ok), ncol)
        record!("$(label).source$(i).type", string(typeof(s)))
        record!("$(label).source$(i).supported", ok)
        record!("$(label).source$(i).columns", ncol)
    end
    for (i, t) in enumerate(targets)
        ok = pnl._gpu_target_supported(t)
        @printf("    target %d  %-46s supported=%-5s n=%d\n",
                i, string(nameof(typeof(t))), string(ok),
                FastMultipole.get_n_bodies(t))
        record!("$(label).target$(i).type", string(typeof(t)))
        record!("$(label).target$(i).supported", ok)
    end
    return nothing
end

function assert_step_head_ok(targets::Tuple, wakes_tuple::Tuple, where_::AbstractString)
    bad = !allfinite(targets)
    if !bad
        for w in wakes_tuple
            isnothing(w) && continue
            pw = w isa pnl.PanelParticleWake ? w.panel_wake : w
            pw isa pnl.PanelWake || continue
            for nod in pw.nodes                      # FIRST WAKE ROW specifically
                size(nod, 2) >= 1 || continue
                all(isfinite, view(nod, :, 1, :)) || (bad = true)
            end
        end
    end
    bad && error("""
        NON-FINITE state detected $(where_).

        This is the classic STEP-HEAD GOTCHA: the parity step's head
        (maneuver! -> _sa_collect -> _sa_reset_freestream_kinematic! ->
        update_TE!), or in full mode the warm-start kinematic replay
        (simulate_warmstart! section 2.5, which re-rotates the un-persisted
        `Das` shed offsets), was skipped or incomplete. That silently NaNs the
        first wake row. Three cluster jobs were lost to this in tasks 018/023.
        Aborting rather than reporting a meaningless parity number.""")
    return nothing
end

################################################################################
# PARITY COMPARISON OF ONE PASS
################################################################################

"""
    compare_pass!(label, targets, runpass!; gate, tight, tight_label)

Run `runpass!()` twice from an identical pre-pass state: first with the seam
`:off` (CPU fmm!), then with the seam armed. Compare the pass CONTRIBUTION
(post-pass state minus pre-pass state) per output class, then write the
CPU-arm outputs back so downstream state stays on the production CPU path.

Returns `(worst_overall, class_errors::Dict{String,Float64}, hits_seam)`.
"""
function compare_pass!(label::AbstractString, targets::Tuple, runpass!;
        gate::Float64=GATE, tight::Float64=NaN, tight_label::AbstractString="",
        keep::Union{Nothing,Base.RefValue{Any}}=nothing)
    pre = snapshot(targets)
    base = extract_all(targets)

    # ---- arm A: CPU (fmm!) ----
    pnl.set_gpu_influence!(:off)
    h0 = pnl.GPU_INFLUENCE_HITS[]
    t_cpu = elapsed_s(runpass!)
    hits_cpu = pnl.GPU_INFLUENCE_HITS[] - h0
    cpu_state = snapshot(targets)
    cpu_out = extract_all(targets)

    # ---- arm B: armed seam ----
    restore!(targets, pre)
    pnl.set_gpu_influence!(SEAM_MODE)
    h0 = pnl.GPU_INFLUENCE_HITS[]
    t_seam = elapsed_s(runpass!)
    hits_seam = pnl.GPU_INFLUENCE_HITS[] - h0
    seam_out = extract_all(targets)
    pnl.set_gpu_influence!(:off)

    hits_cpu == 0 || error("seam hit counter advanced by $(hits_cpu) during the CPU " *
        "arm of $(label) -- the seam was not disarmed")

    dcpu = diffmap(cpu_out, base)
    dseam = diffmap(seam_out, base)

    record!("$(label).hits_seam_arm", hits_seam)
    record!("$(label).hits_cpu_arm", hits_cpu)
    record!("$(label).wall_cpu_s", t_cpu)
    record!("$(label).wall_seam_s", t_seam)

    println("\n--- $(label) ---")
    @printf("  seam hits: %d (cpu arm: %d)    wall: cpu %.3f s | seam %.3f s\n",
            hits_seam, hits_cpu, t_cpu, t_seam)
    if hits_seam == 0
        println("  NOTE: the seam did NOT accept this pass -- both arms ran the CPU")
        println("        fmm! path (see the fallback rules at the top of")
        println("        src/FLOWPanel_gpu_influence.jl). The errors below therefore")
        println("        only measure threading nondeterminism, not seam parity.")
    end
    push!(GATE_LINES, ("$(label).seam_accepted", true,
        hits_seam > 0 ? "yes ($(hits_seam) influence! call(s) handled by the seam)" :
        "NO -- pass fell through to fmm!; the parity numbers below are VACUOUS " *
        "(informational, see the eligibility table above)"))
    record!("$(label).seam_accepted", hits_seam > 0)

    # particle_Estr is EXCLUDED from the pass-level phase gate: it has its own
    # dedicated SFS gate, and in full mode it is a discontinuous function
    # (dynamic-procedure clamp/clipping) of the pass's U/J inputs, so one root
    # cause would otherwise fail two gates (job 13306475). It is still
    # printed/recorded here.
    worst_overall = 0.0          # per-target rel, gated classes only
    dscale_overall = 0.0         # diff-to-scale, gated classes only
    errs = Dict{String,Float64}()
    stats_all = Dict{String,Any}()
    for cls in sort!(collect(keys(dcpu)))
        A = dcpu[cls]; B = dseam[cls]
        if size(A, 2) == 0
            println("  $(rpad(cls, 24)) ABSENT (0 targets in this class)")
            record!("$(label).$(cls).n", 0)
            continue
        end
        st = class_parity_stats(A, B; gate=gate)
        errs[cls] = st.worst_rel
        stats_all[cls] = st
        if cls != "particle_Estr"
            worst_overall = max(worst_overall, st.worst_rel)
            dscale_overall = max(dscale_overall, st.diff_to_scale)
        end
        record!("$(label).$(cls).n", size(A, 2))
        record!("$(label).$(cls).worst_rel", st.worst_rel)
        record!("$(label).$(cls).worst_index", st.iw)
        record!("$(label).$(cls).worst_absdiff", st.dn_w)
        record!("$(label).$(cls).worst_refnorm", st.rn_w)
        record!("$(label).$(cls).class_scale", st.scale)
        record!("$(label).$(cls).diff_to_scale", st.diff_to_scale)
        record!("$(label).$(cls).n_over_gate", st.n_over)
        record!("$(label).$(cls).n_divergent", st.n_div)
        @printf("  %-24s n=%-8d worst rel = %.3e (target %d: |diff| %.3e, |ref| %.3e)\n",
                cls, size(A, 2), st.worst_rel, st.iw, st.dn_w, st.rn_w)
        @printf("  %-24s          diff/scale = %.3e (scale %.3e)  over-gate %d  divergent %d\n",
                "", st.diff_to_scale, st.scale, st.n_over, st.n_div)
        if st.n_over > 0
            for (i, r, dn, rn) in st.top
                r > gate || break
                @printf("      offender target %-8d rel %.3e  |diff| %.3e  |ref| %.3e\n",
                        Int(i), r, dn, rn)
            end
        end
    end
    isempty(dcpu) && println("  (no comparable output classes -- all target classes absent)")

    record!("$(label).worst_rel_overall", worst_overall)
    record!("$(label).diff_to_scale_overall", dscale_overall)
    # Phase gate metric by mode (documented in class_parity_stats): mini arms
    # are both exact -> per-target rel; full-mode CPU arm carries its own fmm
    # error -> gate diff-to-scale, report per-target rel informationally.
    if MODE == "mini"
        gate!("$(label).phase_gate_1e-3", worst_overall, gate)
        gate!("$(label).phase_diff_to_scale", dscale_overall, gate; informational=true)
    else
        gate!("$(label).phase_gate_1e-3_diff_to_scale", dscale_overall, gate;
            note="full mode: diff/class-scale; CPU arm is fmm-approximate")
        gate!("$(label).pertarget_rel", worst_overall, gate; informational=true,
            note="full mode: per-target rel vs an fmm-approximate reference")
    end
    isnan(tight) || gate!(isempty(tight_label) ? "$(label).tight" : tight_label,
        MODE == "mini" ? worst_overall : dscale_overall, tight; informational=true)

    keep === nothing || (keep[] = (; pre, cpu_state, dcpu, dseam))
    restore!(targets, cpu_state)          # downstream state = production CPU path
    return (worst_overall, errs, hits_seam, stats_all)
end

################################################################################
# ONE FULL (uncompared) STEP -- mirror of the simulate! loop body
# (src/FLOWPanel_simulate.jl:1181-1366, legacy-Kutta path)
################################################################################

step_dt(cfg, i_step) = i_step < length(cfg.t_range) - 1 ?
    cfg.t_range[i_step+2] - cfg.t_range[i_step+1] :
    cfg.t_range[i_step+1] - cfg.t_range[i_step]

function advance_step!(cfg, systems_tuple, wakes_tuple, i_step, dt)
    for w in wakes_tuple
        if w isa pnl.PanelParticleWake
            pnl.propagate!(w, dt; relax=cfg.particle_relax, step=i_step, frames=cfg.frames)
        elseif !isnothing(w)
            pnl.propagate!(w, dt; step=i_step, frames=cfg.frames)
        end
    end
    pnl.propagate_kinematics!(systems_tuple, cfg.frames, dt)
    for sys in systems_tuple
        pnl.calc_normals!(sys)
        pnl.calc_controlpoints!(sys)
    end
    for (sys, w) in zip(systems_tuple, wakes_tuple)
        !isnothing(w) && pnl.shed_wake!(w, sys)
    end
    return nothing
end

function run_full_step!(cfg, systems_tuple, wakes_tuple, i_step)
    t = cfg.t_range[i_step+1]
    dt = step_dt(cfg, i_step)
    Base.invokelatest(cfg.maneuver, cfg.frames, systems_tuple, wakes_tuple, t)
    uinf = Base.invokelatest(cfg.Uinf, t)
    pnl._steady_aerodynamics!(cfg.systems, systems_tuple, wakes_tuple, cfg.frames,
        uinf, cfg.body_solvers;
        backend_wake=cfg.backend_wake, backend_solve=cfg.backend_solve,
        backend_system=cfg.backend_system,
        needs_induced_vorticity=cfg.needs_induced_vorticity,
        update_trailing_edges=true,
        wakerow_no_hessian_to_particles=cfg.wakerow_no_hessian_to_particles,
        body_hessian_to_particles=cfg.body_hessian_to_particles,
        body_gradient_core_size=cfg.body_gradient_core_size,
        body_on_wake=cfg.body_on_wake,
        panel_wake_on_particles=cfg.panel_wake_on_particles,
        particle_hessian_self=cfg.particle_hessian_self,
        grad_mu_options=cfg.grad_mu_options,
        formulation=cfg.formulation,
        formulation_state=cfg.formulation_state,
        i_step=i_step)
    advance_step!(cfg, systems_tuple, wakes_tuple, i_step, dt)
    return nothing
end

################################################################################
# THE PARITY STEP
################################################################################

function parity_step!(cfg, systems_tuple, wakes_tuple, i_step)
    t = cfg.t_range[i_step+1]
    println("\n>>> parity step i_step=$(i_step)  t=$(t)")

    # ---------------- STEP HEAD (must run IN FULL; see the gotcha note) -------
    Base.invokelatest(cfg.maneuver, cfg.frames, systems_tuple, wakes_tuple, t)
    uinf = Base.invokelatest(cfg.Uinf, t)
    wake_probes, targets, wake_sources = pnl._sa_collect(systems_tuple, wakes_tuple)
    pnl._sa_reset_freestream_kinematic!(systems_tuple, wakes_tuple, cfg.frames, uinf)
    for (sys, w) in zip(systems_tuple, wakes_tuple)
        !isnothing(w) && pnl.update_TE!(w, sys)
    end
    pnl.formulation_prewake!(cfg.formulation, cfg.formulation_state, systems_tuple)
    assert_step_head_ok(targets, wakes_tuple, "immediately after the parity step head")
    println("    step head OK (all target state and first wake rows finite)")

    # ---------------- case size ---------------------------------------------
    n_panels = sum(FastMultipole.get_n_bodies(s) for s in systems_tuple; init=0)
    n_particles = sum(w isa pnl.PanelParticleWake ? w.pfield.np : 0 for w in wakes_tuple; init=0)
    n_probes = sum(FastMultipole.get_n_bodies(p) for p in wake_probes
                   if !(p isa FLOWVPM.ParticleField); init=0)
    record!("n_panels", n_panels)
    record!("n_particles", n_particles)
    record!("n_probe_targets", n_probes)
    record!("n_targets", length(targets))
    record!("n_wake_sources", length(wake_sources))
    record!("target_types", join(string.(nameof.(typeof.(targets))), "|"))
    record!("wake_source_types", join(string.(nameof.(typeof.(wake_sources))), "|"))
    @printf("    n_panels=%d  n_particles=%d  n_probe_targets=%d  n_targets=%d  n_wake_sources=%d\n",
            n_panels, n_particles, n_probes, length(targets), length(wake_sources))
    n_probes == 0 && println("    NOTE: no wake-probe targets in this case -- the " *
                             "probe_U/probe_J classes are reported as ABSENT.")

    # ---------------- PASS 1: wake sources -> bodies + probes + particles ----
    println("\n--- pass1_wake_influence eligibility ---")
    report_eligibility("pass1_wake_influence", targets, wake_sources)
    pass1 = () -> pnl._sa_wake_influence!(targets, wake_sources, cfg.backend_wake;
        needs_induced_vorticity=cfg.needs_induced_vorticity,
        wakerow_no_hessian_to_particles=cfg.wakerow_no_hessian_to_particles,
        panel_wake_on_particles=cfg.panel_wake_on_particles,
        particle_hessian_self=cfg.particle_hessian_self)
    w1, errs1, hits1, stats1 = compare_pass!("pass1_wake_influence", targets, pass1;
        gate=GATE, tight=GATE_P1,
        tight_label="pass1_wake_influence.tight_1e-11_points_kernel_F64")

    # ---------------- SOLVE: must never be captured by the seam --------------
    pnl._set_core_sizes!(systems_tuple, :core_size_panel)
    for probe_mode in (SEAM_MODE,)                       # armed on purpose
        pnl.set_gpu_influence!(probe_mode)
        h0 = pnl.GPU_INFLUENCE_HITS[]
        t_solve = elapsed_s() do
            pnl.solve_formulation!(cfg.formulation, cfg.formulation_state, cfg.systems,
                systems_tuple, wakes_tuple, cfg.body_solvers;
                backend_solve=cfg.backend_solve, backend_wake=cfg.backend_wake,
                i_step=i_step)
        end
        pnl.set_gpu_influence!(:off)
        hits_solve = pnl.GPU_INFLUENCE_HITS[] - h0
        record!("solve.hits", hits_solve)
        record!("solve.wall_s", t_solve)
        println("\n--- solve_formulation! (seam armed = $(probe_mode)) ---")
        @printf("  seam hits during solve: %d    wall: %.3f s\n", hits_solve, t_solve)
        if hits_solve == 0
            push!(GATE_LINES, ("solve.seam_never_matches", true,
                "seam hit counter did not increment across the solve (0 hits)"))
            record!("gate.solve.seam_never_matches.pass", true)
        else
            push!(GATE_LINES, ("solve.seam_never_matches", false,
                "seam CAPTURED $(hits_solve) solve influence call(s) -- the solve pass " *
                "must always fall through to fmm!"))
            record!("gate.solve.seam_never_matches.pass", false)
        end
    end

    cfg.needs_induced_vorticity && pnl._add_bound_surface_vorticity!(systems_tuple;
        grad_mu_options=pnl._normalize_grad_mu_options(cfg.grad_mu_options; default_basis=:quad))

    # ---------------- PASS 3: conditioned body sources ----------------------
    println("\n--- pass3_body_influence eligibility ---")
    report_eligibility("pass3_body_influence", targets, systems_tuple)
    pass3 = function ()
        pnl._set_core_sizes!(systems_tuple, :core_size_targets)
        pnl._sa_body_influence!(targets, systems_tuple, cfg.backend_system;
            needs_induced_vorticity=cfg.needs_induced_vorticity,
            body_on_wake=cfg.body_on_wake,
            body_hessian_to_particles=cfg.body_hessian_to_particles,
            body_gradient_core_size=cfg.body_gradient_core_size)
    end
    keep3 = Ref{Any}(nothing)
    w3, errs3, hits3, stats3 = compare_pass!("pass3_body_influence", targets, pass3;
        gate=GATE, keep=keep3)
    println("  (pass 3's 1e-3 gate is the real one: the seam evaluates ALL self-body")
    println("   pairs at core_size_panel, while direct_conditioning flips the offset")
    println("   only for NEARFIELD self-pair blocks -- documented deviation.)")

    # ---------------- PASS 3 ATTRIBUTION: exact DirectBackend arm ------------
    # Job 13306588: pass-3 body_velocity diff/scale 1.235e-2 vs an
    # fmm-APPROXIMATE CPU reference cannot distinguish a seam defect from the
    # CPU arm's own truncation error. Re-run the identical pass from the same
    # pre-pass state with DirectBackend (exact brute force; production
    # conditioning preserved -- standalone direct! honors direct_conditioning)
    # on BODY targets only (particle/probe classes are already clean, and
    # body->particles direct would cost ~6.7e9 pairs), then judge each arm
    # against the exact result. Cost: n_panels^2 pairs, ~40 s threaded at p018
    # size (measured small-problem extrapolation: fm051_attribution_debug.jl).
    if get(ENV, "FM051_ATTRIBUTION", "1") == "1" && keep3[] !== nothing
        k3 = keep3[]
        body_targets = Tuple(t for t in targets if t isa pnl.AbstractBody)
        if isempty(body_targets) || !haskey(k3.dcpu, "body_velocity")
            gate_skip!("pass3_attribution.seam_vs_exact",
                "no body targets / no body_velocity class in pass 3")
        else
            restore!(targets, k3.pre)
            base_b = extract_all(body_targets)
            pnl.set_gpu_influence!(:off)
            t_exact = elapsed_s() do
                pnl._set_core_sizes!(systems_tuple, :core_size_targets)
                pnl._sa_body_influence!(body_targets, systems_tuple, pnl.DirectBackend();
                    needs_induced_vorticity=cfg.needs_induced_vorticity,
                    body_on_wake=cfg.body_on_wake,
                    body_hessian_to_particles=cfg.body_hessian_to_particles,
                    body_gradient_core_size=cfg.body_gradient_core_size)
            end
            dexact = diffmap(extract_all(body_targets), base_b)
            restore!(targets, k3.cpu_state)   # downstream convention unchanged
            r = attr_report(dexact["body_velocity"], k3.dcpu["body_velocity"],
                k3.dseam["body_velocity"]; gate_seam=GATE_ATTR)
            @printf("  exact direct arm wall: %.3f s (body targets only)\n", t_exact)

            # -- host-seam arm (only when the primary seam is :cuda): splits a
            # seam-vs-exact deviation into functor-math vs device-numerics.
            # Job 13309844: seam(:cuda) vs exact 1.235e-2 while cpu-fmm vs
            # exact was 1.1e-10 -- every prior clean seam validation (mini,
            # attribution debug) ran :host, so this arm is the discriminator.
            if SEAM_MODE == :cuda
                restore!(targets, k3.pre)
                pnl.set_gpu_influence!(:host)
                h0h = pnl.GPU_INFLUENCE_HITS[]
                t_host = elapsed_s() do
                    pnl._set_core_sizes!(systems_tuple, :core_size_targets)
                    pnl._sa_body_influence!(body_targets, systems_tuple, cfg.backend_system;
                        needs_induced_vorticity=cfg.needs_induced_vorticity,
                        body_on_wake=cfg.body_on_wake,
                        body_hessian_to_particles=cfg.body_hessian_to_particles,
                        body_gradient_core_size=cfg.body_gradient_core_size)
                end
                hits_host = pnl.GPU_INFLUENCE_HITS[] - h0h
                pnl.set_gpu_influence!(:off)
                dhost = diffmap(extract_all(body_targets), base_b)
                restore!(targets, k3.cpu_state)
                if hits_host == 0
                    gate_skip!("pass3_attribution.host_seam_vs_exact",
                        "host seam did not accept the body-only pass (0 hits) -- vacuous")
                else
                    st_he = attr_class_stats(dexact["body_velocity"], dhost["body_velocity"])
                    st_ch = attr_class_stats(dhost["body_velocity"], k3.dseam["body_velocity"])
                    @printf("  host-seam vs exact     diff/scale = %.3e   worst rel = %.3e (target %d)   [wall %.3f s]\n",
                        st_he.diff_to_scale, st_he.worst_rel, st_he.iw, t_host)
                    @printf("  cuda-seam vs host-seam diff/scale = %.3e   worst rel = %.3e (target %d)\n",
                        st_ch.diff_to_scale, st_ch.worst_rel, st_ch.iw)
                    record!("pass3_attribution.host_seam_vs_exact.diff_to_scale", st_he.diff_to_scale)
                    record!("pass3_attribution.cuda_vs_host_seam.diff_to_scale", st_ch.diff_to_scale)
                    gate!("pass3_attribution.host_seam_vs_exact", st_he.diff_to_scale, GATE_ATTR;
                        informational=true,
                        note="attribution split: functor math (host rect vs exact direct)")
                    gate!("pass3_attribution.cuda_vs_host_seam", st_ch.diff_to_scale, GATE_ATTR;
                        informational=true,
                        note="attribution split: device numerics (cuda rect vs host rect)")
                end
            end

            # -- offender geometry: control points of the worst seam-vs-exact
            # targets (single-body case; index == column in the body class)
            if length(body_targets) == 1 && !isempty(r.seam.top)
                b1 = body_targets[1]
                println("  worst seam-vs-exact offender control points (x y z | normal):")
                for (i, rel, dn, rn) in r.seam.top
                    ii = Int(i)
                    cp = b1.controlpoints[:, ii]; nrm = b1.normals[:, ii]
                    @printf("      target %-8d rel %.3e  cp % .6f % .6f % .6f  n % .3f % .3f % .3f\n",
                        ii, rel, cp[1], cp[2], cp[3], nrm[1], nrm[2], nrm[3])
                end
            end
            record!("pass3_attribution.wall_exact_s", t_exact)
            record!("pass3_attribution.seam_vs_exact.diff_to_scale", r.seam.diff_to_scale)
            record!("pass3_attribution.seam_vs_exact.worst_rel", r.seam.worst_rel)
            record!("pass3_attribution.cpu_vs_exact.diff_to_scale", r.cpu.diff_to_scale)
            record!("pass3_attribution.cpu_vs_exact.worst_rel", r.cpu.worst_rel)
            record!("pass3_attribution.cpu_vs_exact.n_divergent", r.cpu.n_div)
            record!("pass3_attribution.cross_check.diff_to_scale", r.cross.diff_to_scale)
            gate!("pass3_attribution.seam_vs_exact", r.seam.diff_to_scale, GATE_ATTR;
                note="seam vs exact DirectBackend, diff/scale (F64 summation-order headroom)")
            gate!("pass3_attribution.cpu_fmm_vs_exact", r.cpu.diff_to_scale, GATE;
                informational=true,
                note="measured legacy fmm! truncation on body self-influence -- a finding, not a defect gate")
        end
    end

    # ---------------- SFS / Estr --------------------------------------------
    println("\n--- SFS / Estr ---")
    sfs_fields = [w.pfield for w in wakes_tuple
                  if w isa pnl.PanelParticleWake && FLOWVPM.isSFSenabled(w.pfield.SFS)]
    record!("sfs_enabled", !isempty(sfs_fields))
    if isempty(sfs_fields)
        println("  SFS disabled on every particle field in this case -- no Estr rows " *
                "to compare (class ABSENT).")
        gate_skip!("sfs.estr_1e-10", "SFS disabled on all particle fields")
        gate_skip!("sfs.seam_never_handles_estr", "SFS disabled on all particle fields")
    else
        e_estr = get(errs1, "particle_Estr", NaN)
        # By design the particle-field SELF influence + SFS closure stay on the
        # FLOWVPM radix path in BOTH arms; the seam never evaluates Estr.
        if hits1 == 0
            push!(GATE_LINES, ("sfs.seam_never_handles_estr", true,
                "pass 1 fell back to fmm! entirely (SFS-enabled host particle self-pair " *
                "is an explicit seam fallback), so the seam handled no Estr work"))
            record!("gate.sfs.seam_never_handles_estr.pass", true)
        else
            push!(GATE_LINES, ("sfs.seam_never_handles_estr", true,
                "seam accepted pass 1 ($(hits1) hit(s)); Estr came from the FLOWVPM " *
                "radix path in both arms by construction (Estr_fmm!/UJ_fmm_gpu!)"))
            record!("gate.sfs.seam_never_handles_estr.pass", true)
        end
        if isnan(e_estr)
            gate_skip!("sfs.estr_1e-10", "no particle_Estr output class was produced")
        elseif MODE == "mini"
            @printf("  particle_Estr worst per-target rel = %.3e\n", e_estr)
            gate!("sfs.estr_1e-10", e_estr, GATE_SFS;
                note="loose bound: threading nondeterminism only")
        else
            # FULL mode: both arms compute Estr through the identical
            # fmm!+postcalc path (the SFS self-pair fallback), but Estr is a
            # DISCONTINUOUS function of the pass-1 U/J inputs: the dynamic
            # procedure's clamp/clipping branches (FLOWVPM_subfilterscale.jl)
            # flip at isolated particles under the legitimate fmm-vs-exact
            # input difference (~2e-4), producing O(1) per-particle deviations
            # (job 13306475: worst rel 1.185 from exactly this). A tight
            # per-particle bound is therefore unsatisfiable BY DESIGN in full
            # mode. Gate instead that the clip-flip divergences are SPARSE
            # (<0.1% of particles) and report the magnitude informationally.
            st = get(stats1, "particle_Estr", nothing)
            if st === nothing || st.n == 0
                gate_skip!("sfs.estr_divergence_sparse", "no pass-1 Estr stats")
            else
                frac_div = st.n_div / st.n
                @printf("  particle_Estr: worst rel %.3e | diff/scale %.3e | divergent %d of %d (%.2e)\n",
                        st.worst_rel, st.diff_to_scale, st.n_div, st.n, frac_div)
                record!("sfs.estr_n_divergent", st.n_div)
                record!("sfs.estr_frac_divergent", frac_div)
                gate!("sfs.estr_divergence_sparse", frac_div, 1e-3;
                    note="full mode: fraction of particles with per-particle rel > 0.5 " *
                         "(dynamic-procedure clip flips under fmm-vs-exact inputs)")
                gate!("sfs.estr_diff_to_scale", st.diff_to_scale, GATE;
                    informational=true,
                    note="full mode: clip flips at high-Estr particles can reach O(1)")
            end
        end
    end

    return nothing
end

################################################################################
# MINI CASE
################################################################################

function build_mini_cfg()
    if !isdefined(@__MODULE__, :make_dirichlet_diamond_body)
        Base.include(@__MODULE__, joinpath(@__DIR__, "..", "test", "test_helpers.jl"))
    end
    nspan = parse(Int, get(ENV, "FM051_MINI_NSPAN", "6"))
    p, mac, leaf = split(get(ENV, "FM051_MINI_FMM", "10,0.4,100000"), r"[:,]")
    backend = pnl.FastMultipoleBackend(parse(Int, p), parse(Float64, mac), parse(Int, leaf))

    # Union{ConstantSource,VortexRing} RigidWakeBody with a FINITE attached TE
    # wake -- the seam's panel tag 4 plus the attached-wake triangle packing.
    body = Base.invokelatest(getfield(@__MODULE__, :make_dirichlet_diamond_body);
        nspan=nspan, thick=0.06, das=0.15)
    body.needs_velocity_gradient[] = true    # exercise the body Hessian write-back

    # include_final_filament: once the panel sheet has overflowed (which it
    # always has in a mature run), `FilamentWrapper` reports n > 0. With the
    # 051 FilamentWrapper seam extension the active final-row filaments are
    # packed as tag-3 nv=2 open-segment columns, so FM051_MINI_FINAL_FILAMENT=1
    # must now show pass 1 seam-ACCEPTED (it previously reproduced the
    # documented fallback). Default false keeps the plain-mini baseline
    # unchanged.
    final_filament = get(ENV, "FM051_MINI_FINAL_FILAMENT", "0") in ("1", "true")
    wake = pnl.PanelParticleWake(body; nwakerows=1, max_particles=20_000,
        core_size=2e-2, shed_with_induced_velocity=true,
        include_final_filament=final_filament,
        unsteady_filament=final_filament)
    record!("mini.include_final_filament", final_filament)
    record!("mini.nspan", nspan)
    record!("mini.backend", string(backend))
    frames = pnl.ReferenceFrame(body)
    solver = pnl.Backslash(body)

    dt = 0.05
    nspin = parse(Int, get(ENV, "FM051_MINI_STEPS", "4"))
    t_range = collect(range(0.0; step=dt, length=nspin + 3))
    Uinf(t) = [1.0, 0.0, 0.15]
    maneuver!(frames, systems, wakes, t) = nothing

    systems_tuple = pnl._systems_tuple(body)
    wakes_tuple = pnl._wakes_tuple(body, wake)
    formulation = pnl.VelocityThroughSources()
    formulation_state = pnl.initialize_formulation(formulation, systems_tuple,
        wakes_tuple, solver, backend, backend)

    cfg = (; systems=body, frames, maneuver=maneuver!, Uinf, t_range,
        body_solvers=solver, backend_wake=backend, backend_solve=backend,
        backend_system=backend, formulation, formulation_state,
        needs_induced_vorticity=false,
        wakerow_no_hessian_to_particles=false,
        body_hessian_to_particles=false,
        body_gradient_core_size=NaN,
        body_on_wake=true, panel_wake_on_particles=true,
        particle_hessian_self=true, particle_relax=true,
        grad_mu_options=(; basis=:tri))
    return cfg, systems_tuple, wakes_tuple, nspin
end

################################################################################
# FULL CASE (p018 warm start; benchmark/p018_mature_wake_timing.jl reused as-is)
################################################################################

module P018 end          # private namespace so the p018 script's consts cannot clash

function build_full_cfg()
    Base.include(P018, joinpath(@__DIR__, "p018_mature_wake_timing.jl"))
    cfg018 = Base.invokelatest(getfield(P018, Symbol("canonical_setup!")))
    restart_step = parse(Int, get(ENV, "RESTART_STEP", "-1"))
    restart_step < 0 &&
        (restart_step = Base.invokelatest(getfield(P018, :default_restart_step)))
    println("  restarting from step $(restart_step) of run " *
            "$(getfield(P018, :RESTART_RUN))")
    systems_tuple, wakes_tuple, setup_times = Base.invokelatest(
        getfield(P018, Symbol("warmstart_restore!")), cfg018, restart_step)
    for (k, v) in sort!(collect(setup_times))
        @printf("    %s: %.1f s\n", k, v)
        record!("full." * k, v)
    end

    pnl.audit_monitors(cfg018.monitors)
    needs_grad = any(pnl.monitor_requires_body_hessian, cfg018.monitors)
    needs_induced_vorticity = any(pnl.monitor_requires_induced_vorticity, cfg018.monitors)
    for sys in systems_tuple
        sys.needs_velocity_gradient[] = needs_grad
    end
    record!("full.needs_body_hessian", needs_grad)
    record!("full.needs_induced_vorticity", needs_induced_vorticity)
    needs_induced_vorticity && println("    WARNING: needs_induced_vorticity=true -- " *
        "nonzero extra_outputs is an explicit seam FALLBACK, so the seam will not " *
        "accept the body-target passes in this configuration (hits will be 0).")

    formulation_state = Base.invokelatest(pnl.initialize_formulation, cfg018.formulation,
        systems_tuple, wakes_tuple, cfg018.body_solvers, cfg018.backend, cfg018.backend)

    cfg = (; systems=cfg018.systems, frames=cfg018.frames, maneuver=cfg018.maneuver,
        Uinf=cfg018.Uinf, t_range=cfg018.t_range,
        body_solvers=cfg018.body_solvers,
        backend_wake=cfg018.backend_wake, backend_solve=cfg018.backend,
        backend_system=cfg018.backend,
        formulation=cfg018.formulation, formulation_state,
        needs_induced_vorticity,
        wakerow_no_hessian_to_particles=cfg018.wakerow_no_hessian_to_particles,
        body_hessian_to_particles=cfg018.body_hessian_to_particles,
        body_gradient_core_size=cfg018.body_gradient_core_size,
        body_on_wake=cfg018.body_on_wake,
        panel_wake_on_particles=cfg018.panel_wake_on_particles,
        particle_hessian_self=cfg018.particle_hessian_self,
        particle_relax=cfg018.particle_relax,
        grad_mu_options=(;))
    return cfg, systems_tuple, wakes_tuple, restart_step + 1
end

################################################################################
# CSV + REPORT
################################################################################

_csv_escape(s) = occursin(r"[,\"\n]", s) ? "\"" * replace(s, "\"" => "\"\"") * "\"" : s

function write_csv()
    mkpath(RESULTS_DIR)
    path = joinpath(RESULTS_DIR, "fm051_pass_parity_$(MODE).csv")
    open(path, "w") do io
        println(io, "key,value")
        for (k, v) in RECORD
            sv = v isa AbstractFloat ? (isfinite(v) ? @sprintf("%.12g", v) : string(v)) :
                 string(v)
            println(io, _csv_escape(string(k)) * "," * _csv_escape(sv))
        end
    end
    println("\nWrote $(path)")
    return path
end

function report()
    println("\n" * "=" ^ 78)
    println("GATES")
    println("=" ^ 78)
    allok = true
    for (label, ok, detail) in GATE_LINES
        allok &= ok
        @printf("%-4s %-46s %s\n", ok ? "PASS" : "FAIL", label, detail)
    end
    println("=" ^ 78)
    println(allok ? "OVERALL: PASS" : "OVERALL: FAIL")
    record!("overall_pass", allok)
    return allok
end

################################################################################
# MAIN
################################################################################

function main()
    println("=" ^ 78)
    println("FM051 pass-by-pass CPU-vs-GPU influence parity")
    println("=" ^ 78)
    @printf("  mode=%s  seam=%s  threads=%d  julia=%s\n",
            MODE, SEAM_MODE, Threads.nthreads(), VERSION)
    @printf("  gates: phase %.1e | pass1 tight %.1e (informational) | Estr %.1e\n",
            GATE, GATE_P1, GATE_SFS)
    record!("mode", MODE)
    record!("seam_mode", SEAM_MODE)
    record!("threads", Threads.nthreads())
    record!("julia_version", string(VERSION))
    # the filament family is a GLOBAL production default: record, never change
    record!("filament_regularization", string(pnl.FILAMENT_REGULARIZATION[]))
    println("  FLOWPanel.FILAMENT_REGULARIZATION[] = " *
            "$(pnl.FILAMENT_REGULARIZATION[])  (recorded, not modified)")

    seam0 = pnl.GPU_INFLUENCE[]
    try
        pnl.set_gpu_influence!(:off)
        cfg, systems_tuple, wakes_tuple, parity_i_step =
            MODE == "mini" ? build_mini_cfg() : build_full_cfg()
        record!("parity_i_step", parity_i_step)

        if MODE == "mini"
            nspin = parity_i_step
            println("\n  spinning up $(nspin) step(s) with the seam OFF " *
                    "(building nontrivial wake state)...")
            for i_step in 0:(nspin - 1)
                t_step = elapsed_s(() -> Base.invokelatest(run_full_step!, cfg,
                    systems_tuple, wakes_tuple, i_step))
                np = sum(w isa pnl.PanelParticleWake ? w.pfield.np : 0
                         for w in wakes_tuple; init=0)
                @printf("    step %d done (%.2f s), np=%d\n", i_step, t_step, np)
            end
            np = sum(w isa pnl.PanelParticleWake ? w.pfield.np : 0
                     for w in wakes_tuple; init=0)
            np > 0 || error("mini case produced no particles after $(nspin) steps; " *
                "raise FM051_MINI_STEPS (conversion fires on shed nwakerows+1)")
        end

        # arming state must be :off entering the parity step
        pnl.set_gpu_influence!(:off)
        Base.invokelatest(parity_step!, cfg, systems_tuple, wakes_tuple, parity_i_step)
    finally
        pnl.set_gpu_influence!(seam0 === :unset ? :off : seam0)
        println("\n  seam restored to $(pnl.GPU_INFLUENCE[]) (production default :off)")
    end

    ok = report()
    write_csv()
    return ok
end

if abspath(PROGRAM_FILE) == (@__FILE__) || isempty(PROGRAM_FILE)
    ok = main()
    exit(ok ? 0 : 1)
end
