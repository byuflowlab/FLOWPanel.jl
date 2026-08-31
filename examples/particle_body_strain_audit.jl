#!/usr/bin/env julia
# Body/panel-wake strain omission audit (BRAINSTORM item 008 follow-up).
#
# The rotor_hover_pressure_comparison run advects particles with the TOTAL
# velocity (body + panel-wake rows + particles) but stretches them with the
# particle-only velocity gradient (BODY_HESSIAN_TO_PARTICLES=false,
# PANEL_WAKE_HESSIAN_TO_PARTICLES=false). This script sizes the omitted
# stretching terms offline on a saved step:
#
#   1. Rebuild the example systems (RHPC_SETUP_ONLY=true), warmstart-load the
#      solved body + wake state at STEP.
#   2. Evaluate body->particle and panel-wake-row->particle velocity gradients
#      with the same influence machinery the live run would have used.
#   3. Gate: body+wakerow+particle-particle induced velocity at the particles
#      must reproduce the saved particle `velocity`.
#   4. Report per-particle stretching-rate ratios ||S(J_omitted, Γ)|| /
#      ||S(J_particle, Γ)|| banded by axial position (wake age proxy).
#
# ENV: STEP (350), NBANDS (10), AXIAL_DIM (1). Run from repo root:
#   STEP=350 julia --project=. -t auto examples/particle_body_strain_audit.jl

haskey(ENV, "RHPC_SETUP_ONLY") || (ENV["RHPC_SETUP_ONLY"] = "true")
@assert ENV["RHPC_SETUP_ONLY"] == "true"

step = parse(Int, get(ENV, "STEP", "350"))
nbands = parse(Int, get(ENV, "NBANDS", "10"))
axial_dim = parse(Int, get(ENV, "AXIAL_DIM", "1"))

import Printf: @printf
import Statistics: mean, median, quantile
import LinearAlgebra: norm

include(joinpath(@__DIR__, "rotor_hover_pressure_comparison.jl"))
import FLOWVPM

const RUN = "rotor_hover_pressure_comparison"
const DATA = joinpath("data", RUN)

# --- load saved state at `step` into the constructed objects ---------------
pnl._load_body_vtk!(rotor, DATA, "$(RUN)_body1", step)
pnl._load_panel_particle_wake_vtk!(wake_rotor, DATA, "$(RUN)_wake1", step)
pf = wake_rotor.pfield
np = pf.np
@printf("Loaded step %d: np=%d particles, body strength norm=%.6g\n",
        step, np, norm(rotor.strength))

Xp = copy(pf.particles[FLOWVPM.X_INDEX, 1:np])          # 3×np positions
Gp = copy(pf.particles[FLOWVPM.GAMMA_INDEX, 1:np])      # 3×np strengths
J_saved = copy(pf.particles[FLOWVPM.J_INDEX, 1:np])     # 9×np particle-only J
U_saved = copy(pf.particles[FLOWVPM.U_INDEX, 1:np])     # 3×np total velocity

probes = pnl._collect_wake_probes((wake_rotor,))
wake_sources = pnl._collect_wake_sources((wake_rotor,))
panel_sources = Tuple(s for s in wake_sources if !(s isa FLOWVPM.ParticleField))
pfield_sources = Tuple(s for s in wake_sources if s isa FLOWVPM.ParticleField)

function zero_probe_state!(pf, np)
    pf.particles[FLOWVPM.U_INDEX, 1:np] .= 0
    pf.particles[FLOWVPM.J_INDEX, 1:np] .= 0
    return nothing
end

grab_J(pf, np) = copy(pf.particles[FLOWVPM.J_INDEX, 1:np])
grab_U(pf, np) = copy(pf.particles[FLOWVPM.U_INDEX, 1:np])

# --- pass 1: body -> particles (velocity + gradient) ------------------------
# Mirrors the body-on-wake influence in _steady_aerodynamics! (core_size
# :core_size_targets, precalc=false) but with the Hessian ENABLED for the
# particle target.
zero_probe_state!(pf, np)
pnl._set_core_sizes!((rotor,), :core_size_targets)
pnl.influence!(probes, (rotor,), backend;
    precalc=false, scalar_potential=false, velocity=true,
    velocity_gradient=Tuple(true for _ in probes),
    direct_conditioning=pnl._self_panel_core_size_conditioning())
J_body = grab_J(pf, np)
U_body = grab_U(pf, np)

# --- pass 2: panel-wake rows -> particles ------------------------------------
zero_probe_state!(pf, np)
if length(panel_sources) > 0
    pnl.influence!(probes, panel_sources, backend_wake;
        precalc=true, scalar_potential=false, velocity=true,
        velocity_gradient=Tuple(true for _ in probes))
end
J_wrow = grab_J(pf, np)
U_wrow = grab_U(pf, np)

# --- pass 3 (gate): particles -> particles velocity -------------------------
# Total induced velocity must reproduce the saved particle `velocity`
# (freestream is withdrawn in settled hover).
zero_probe_state!(pf, np)
pnl.influence!(pfield_sources, pfield_sources, backend_wake;
    precalc=true, postcalc=true, scalar_potential=false, velocity=true,
    velocity_gradient=(true,))
J_pp = grab_J(pf, np)
U_pp = grab_U(pf, np)

U_total = U_body .+ U_wrow .+ U_pp
gate_vel = norm(U_total .- U_saved) / (norm(U_saved) + 1e-30)
gate_J = norm(J_pp .- J_saved) / (norm(J_saved) + 1e-30)
@printf("GATE  |U_body+U_wrow+U_pp - U_saved|/|U_saved| = %.4e\n", gate_vel)
@printf("GATE  |J_pp - J_saved|/|J_saved|               = %.4e\n", gate_J)

# --- stretching-term ratios --------------------------------------------------
# Live FLOWVPM transposed branch: S(J, Γ). The omitted per-particle strength
# rates are S(J_body,Γ) and S(J_wrow,Γ); the retained one is S(J_saved,Γ).
function rates(J, G, np)
    S = Matrix{Float64}(undef, 3, np)
    @inbounds for i in 1:np
        Jv = view(J, :, i); Gv = view(G, :, i)
        S[1,i] = Jv[1]*Gv[1] + Jv[2]*Gv[2] + Jv[3]*Gv[3]
        S[2,i] = Jv[4]*Gv[1] + Jv[5]*Gv[2] + Jv[6]*Gv[3]
        S[3,i] = Jv[7]*Gv[1] + Jv[8]*Gv[2] + Jv[9]*Gv[3]
    end
    return S
end

S_part = rates(J_saved, Gp, np)
S_body = rates(J_body, Gp, np)
S_wrow = rates(J_wrow, Gp, np)

nS_part = vec(sqrt.(sum(abs2, S_part; dims=1)))
nS_body = vec(sqrt.(sum(abs2, S_body; dims=1)))
nS_wrow = vec(sqrt.(sum(abs2, S_wrow; dims=1)))
floorS = 1e-6 * median(nS_part)
ratio_body = nS_body ./ max.(nS_part, floorS)
ratio_wrow = nS_wrow ./ max.(nS_part, floorS)

@printf("\nGlobal Frobenius ratios: |S_body|/|S_part|=%.4f  |S_wrow|/|S_part|=%.4f\n",
        norm(S_body)/norm(S_part), norm(S_wrow)/norm(S_part))

# --- axial banding (same convention as particle_overlap_age_diag.jl) --------
axial = vec(Xp[axial_dim, :])
ord = sortperm(axial)
bands = [ord[(1 + div((b-1)*np, nbands)):div(b*np, nbands)] for b in 1:nbands]

rows = NamedTuple[]
@printf("\n%-4s %7s  [%9s, %9s]  %9s %9s  %9s %9s  %9s\n",
        "band", "count", "ax_lo", "ax_hi", "body_med", "body_q90",
        "wrow_med", "wrow_q90", "bodyfrob")
for (b, S) in enumerate(bands)
    frob_body = norm(@view S_body[:, S]) / (norm(@view S_part[:, S]) + 1e-30)
    r = (; run_name=RUN, step, band=b, count=length(S),
        ax_lo=minimum(axial[S]), ax_hi=maximum(axial[S]),
        ratio_body_med=median(ratio_body[S]), ratio_body_q90=quantile(ratio_body[S], 0.9),
        ratio_wrow_med=median(ratio_wrow[S]), ratio_wrow_q90=quantile(ratio_wrow[S], 0.9),
        ratio_body_frob=frob_body,
        Ubody_med=median(vec(sqrt.(sum(abs2, U_body[:, S]; dims=1)))))
    push!(rows, r)
    @printf("%-4d %7d  [%9.4f, %9.4f]  %9.4f %9.4f  %9.4f %9.4f  %9.4f\n",
            b, r.count, r.ax_lo, r.ax_hi, r.ratio_body_med, r.ratio_body_q90,
            r.ratio_wrow_med, r.ratio_wrow_q90, r.ratio_body_frob)
end

outdir = joinpath("data", "overlap_age_diag")
mkpath(outdir)

# Dump the omitted gradients for reuse (e.g. particle_overlap_age_diag.jl
# EXTRA_J mode adds them to the saved particle-only J).
import Serialization
extraj_path = joinpath(outdir, "$(RUN)_step$(step)_extraJ.bin")
open(extraj_path, "w") do io
    Serialization.serialize(io, (; step, np, J_body, J_wrow, U_body, U_wrow))
end
@printf("Wrote %s\n", extraj_path)

out = joinpath(outdir, "$(RUN)_step$(step)_body_strain.csv")
open(out, "w") do io
    names = propertynames(rows[1])
    println(io, join(names, ","))
    for r in rows
        println(io, join((getproperty(r, n) for n in names), ","))
    end
end
@printf("\nWrote %s\n", out)
