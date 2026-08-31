#=##############################################################################
BRAINSTORM 021 Phase 1 — shared campaign-case construction (RHPC/Dirichlet).

Factored 2026-08-15 from the identical blocks in
rotor_hover_solver_phase1{,_availability,_fgs_check,_solvetime}.jl so the
BC-metric-era scripts (bcerror validation, FGS τ-ladder tuner, preconditioner
selection, tuned-knob solvetime) cannot drift from each other. The four
pre-existing scripts are left untouched (their CSVs are evidence of record).

Include AFTER common.jl. Expects `RUNG` in ENV. Defines (non-exhaustive):
  rung, msh_name, n_expected, rotor, frames, R, rho, RPM, Vinf, dt,
  frozen_velocity, potential_frozen, b, rms_b, solution_column,
  reset_cold!(), LADDER, TUNED, PROD, banner, results_dir, k_reps

Contracts (phase_01_consistency.md):
- Dirichlet solution = strength column 2; col 1 = BC-set sources.
- Frozen b = −potential after DIRECT source-potential assembly
  (assemble_rhs! contract); reset_cold!() restores the cached assembly so
  timed solves exclude the N² source assembly.
- reset_cold!() gives a genuine cold solve (FGS seeds from body.strength).
- Body core_size is left at the PANEL offset after construction; scripts
  that run steady! target passes must restore it before solver ctors.
=###############################################################################

import FLOWPanel as pnl
import FastMultipole
include(joinpath(pnl.examples_path, "dji9443_trailing_edge.jl"))

using FLOWPanel.FastMultipole.StaticArrays
using LinearAlgebra: norm, cross

banner = assert_and_banner()
# PER_RUNG_DIR=1 (HPC): isolate each rung's CSVs in results/phase1/<mode>/<rung>/
# — concurrent per-rung jobs on different nodes append over NFS, and
# non-atomic appends CLOBBER shared files (observed 2026-08-18: R1's
# staircase rows destroyed by R2–R4 writers). Merge per-rung dirs for
# analysis; local sequential runs keep the flat layout.
_rung_sub = get(ENV, "PER_RUNG_DIR", "0") == "1" ? get(ENV, "RUNG", "") : ""
results_dir = joinpath(@__DIR__, "results", "phase1", banner.threading_mode, _rung_sub)
mkpath(results_dir)
# knobs_dir: where tuning outputs (tune.csv, fgstune_*, fgsprecond.csv) and
# the frozen-b cache are read from. Tuning runs in multi mode only; single-
# mode table jobs pass KNOBS_MODE=multi to reuse them (b and knobs are
# mode-independent — pure functions of geometry/BC and the tuner's choices).
knobs_dir = joinpath(@__DIR__, "results", "phase1",
                     get(ENV, "KNOBS_MODE", banner.threading_mode), _rung_sub)
mkpath(knobs_dir)

const LADDER = Dict(
    "R1" => ("dji9443_20260813_23_73_capped_captess4.msh",    8016),
    "R2" => ("dji9443_20260813_33_105_capped_captess4.msh",  15760),
    "R3" => ("dji9443_20260813_45_145_capped_captess4.msh",  28752),
    "R4" => ("dji9443_20260813_65_209_capped_captess4.msh",  58192),
    "R5" => ("dji9443_20260813_89_289_capped_captess4.msh", 108240),
    "R6" => ("dji9443_20260813_125_409_capped_captess4.msh", 212108),
    "R7" => ("dji9443_20260813_177_577_capped_captess4.msh", 419276),
)
# tune_fmm output per rung — the shared Krylov apply knobs (ruling 3).
# Precedence: the knobs_dir tune.csv's latest "tuned"-variant row for the
# rung (produced by the tune stage ON THE BENCHMARK HARDWARE) wins; the
# hardcoded R1–R3 values (local Mac tune.csv 2026-08-14) are a FALLBACK
# only. Apply knobs are hardware-dependent in cost; certification is
# per-run either way. Scripts needing TUNED on an untuned rung fail loudly.
# WARNING (2026-08-24): the hardcoded R1–R3 triples predate the cost-tuner
# fix — they came from the accuracy-only objective, which is blind to solve
# cost and picks leaf_size ~2x too large (R1: leaf 21 here vs ~9 measured
# fastest). They are a last-resort fallback, NOT valid campaign knobs. Every
# rung is retuned from scratch; if you see this fallback fire, the rung's
# tune stage did not run.
const TUNED = Dict{String,NTuple{3,Any}}(
    "R1" => (17, 0.5, 21),
    "R2" => (17, 0.5, 12),
    "R3" => (18, 0.5, 18),
)
let rung_ = get(ENV, "RUNG", "")
    if rung_ != ""
        tune_csv = joinpath(knobs_dir, "tune.csv")
        if isfile(tune_csv)
            lines = readlines(tune_csv)
            cols = Dict(String(c) => i for (i, c) in enumerate(split(lines[1], ",")))
            for l in lines[2:end]
                c = split(l, ",")
                length(c) >= length(cols) || continue
                c[cols["rung"]] == rung_ || continue
                c[cols["variant"]] == "tuned" || continue
                TUNED[rung_] = (parse(Int, c[cols["expansion_order"]]),
                                parse(Float64, c[cols["multipole_acceptance"]]),
                                parse(Int, c[cols["leaf_size"]]))   # latest wins
            end
        end
        haskey(TUNED, rung_) || @warn "No tuned apply knobs for $rung_ — run the " *
            "tune stage first; scripts touching TUNED[rung] will fail" rung_
        haskey(TUNED, rung_) && println("TUNED[$rung_] = $(TUNED[rung_]) " *
            "(csv-first precedence; hardcoded R1-R3 are fallback)")
    end
end
# Phase 0 production FGS knobs (fgs_check.csv fastest-or-tied at every rung)
const PROD = (10, 0.4, 150)

rung = get(ENV, "RUNG", "")
haskey(LADDER, rung) || error("RUNG must be one of $(sort(collect(keys(LADDER)))); got \"$rung\"")
msh_name, n_expected = LADDER[rung]
k_reps = parse(Int, get(ENV, "K_REPS", "1"))

# ---- Dirichlet case (RHPC: rotor_hover_pressure_comparison conventions) ----
magVinf = 0.0001; rho = 1.179; RPM = 6000; R = 0.119
Vinf = magVinf * [1.0, 0.0, 0.0]
dt = 60 / RPM / 36
core_size_panel = R * 1e-10
core_size_targets = 1e-3
kernelcutoff = R * 1e-13
radial_dimension = 2
shedding_r_over_R = 0.1

msh_file = joinpath(pnl.examples_path, "data", msh_name)
te_indices_1, te_indices_2 = find_dji9443_trailing_edge_indices(msh_file; watertight=true)
msh = pnl.read_gmsh(msh_file)
nodes0, cells0 = pnl.meshes2nodes_cells(msh)
nodes0 .*= R / maximum(nodes0[radial_dimension, :])

kernel = Union{pnl.ConstantSource, pnl.VortexRing}
body_kwargs = (; core_size=core_size_panel, core_size_panel,
                 core_size_targets, kernelcutoff,
                 semiinfinite_wake=false, watertight=true, DBC=true)
rotor0 = pnl.RigidWakeBody{kernel}(nodes0, cells0, pnl.noshedding; body_kwargs...)

function make_shedding_bbox(nodes, seed_nodes, radial_dimension, R, cutoff_r_over_R)
    radial_midpoint = sum(nodes[radial_dimension, seed_nodes]) / length(seed_nodes)
    radial_sign = sign(radial_midpoint)
    lower = [minimum(nodes[i, :]) for i in 1:size(nodes, 1)]
    upper = [maximum(nodes[i, :]) for i in 1:size(nodes, 1)]
    padding = max(sqrt(eps(eltype(nodes))) * R, R * 1e-6)
    lower .-= padding; upper .+= padding
    if radial_sign > 0
        lower[radial_dimension] = cutoff_r_over_R * R - padding
    else
        upper[radial_dimension] = -cutoff_r_over_R * R + padding
    end
    return (pnl.SVector{3}(lower...), pnl.SVector{3}(upper...))
end
function clip_shedding_root(nodes, shedding, cells, radial_dimension, R, clip_r_over_R)
    keep = Int[]
    for j in axes(shedding, 2)
        p, nia, nib = shedding[1, j], shedding[2, j], shedding[3, j]
        na, nb = cells[nia, p], cells[nib, p]
        mid = (nodes[radial_dimension, na] + nodes[radial_dimension, nb]) / 2
        abs(mid) / R >= clip_r_over_R && push!(keep, j)
    end
    return shedding[:, keep]
end
sheddings = map((te_indices_1, te_indices_2)) do te
    bbox = make_shedding_bbox(rotor0.nodes, te[1:2], radial_dimension, R, 0.0)
    full = pnl.calc_shedding_from_seed(rotor0.nodes, rotor0.cells, te[1], te[2];
        bbox, end_node=te[3], normal_jump_tol=0.2, max_turn_angle=pi/3)
    clip_shedding_root(rotor0.nodes, full, rotor0.cells, radial_dimension, R, shedding_r_over_R)
end
rotor = pnl.RigidWakeBody{kernel}(rotor0.nodes, rotor0.cells, collect(sheddings);
            ensure_winding=true, body_kwargs...)
rotor.ncells == n_expected ||
    error("Rung $rung: expected $n_expected panels, got $(rotor.ncells)")

frames = pnl.ReferenceFrame(rotor;
    origin=SVector{3}(0.0,0.0,0.0), v=SVector{3}(0.0,0.0,0.0),
    ω_axis=SVector{3}(-1.0,0.0,0.0), ω=2*pi*RPM/60,
    R=SMatrix{3,3}(1.0,0,0, 0,1.0,0, 0,0,1.0),
    name="vehicle", child_index=Int[], dependent_index=[1])
pnl.initialize_Das!((rotor,), frames, t -> Vinf, 0.0, dt; set_Das_eta_kinematic=0.2)
ω_vec = SVector{3}(-2*pi*RPM/60, 0.0, 0.0)
pnl.calc_normals!(rotor); pnl.calc_controlpoints!(rotor)
for i in 1:rotor.ncells
    rr = SVector{3}(rotor.controlpoints[:, i]...)
    rotor.velocity[:, i] .= Vinf .- cross(ω_vec, rr)
end
frozen_velocity = copy(rotor.velocity)

# frozen Dirichlet b (direct source-potential assembly) + cached potential:
# rhs = -potential is a pure function of geometry + BC, so the assembly is
# restored (not recomputed) before every solve. CACHE_B=1 (HPC, large rungs)
# persists the assembly to disk — the O(N²) direct pass at R7 costs tens of
# minutes per script invocation otherwise.
# SKIP_B=1 (021 W0, 2026-08-24): geometry-only consumers — e.g. the
# near-field-cache feasibility map, which needs trees/direct lists and never
# solves — skip the assembly entirely. rms_b is then NaN so any script that
# derives a tolerance from it fails loudly rather than solving against b=0.
b_skipped = get(ENV, "SKIP_B", "0") == "1"
pnl.set_strengths!(rotor)
rotor.potential .= 0
bcache_path = joinpath(knobs_dir, "bcache_$(rung).bin")
if b_skipped
    t_rhs_assembly = 0.0
    @warn "SKIP_B=1: frozen b NOT assembled for $rung — rms_b is NaN and this " *
          "case MUST NOT be solved (geometry-only consumers only)"
elseif get(ENV, "CACHE_B", "0") == "1" && isfile(bcache_path)
    v = Vector{Float64}(undef, rotor.ncells)
    read!(bcache_path, v)
    rotor.potential .= v
    t_rhs_assembly = 0.0
    println("frozen b loaded from $bcache_path")
else
    t_rhs_assembly = @elapsed pnl.influence!(rotor, rotor, pnl.DirectBackend();
        scalar_potential=true, velocity=false)
    if get(ENV, "CACHE_B", "0") == "1"
        write(bcache_path, Vector{Float64}(rotor.potential))
        println("frozen b cached to $bcache_path")
    end
end
potential_frozen = copy(rotor.potential)
b = zeros(rotor.ncells)
b_skipped || pnl.assemble_rhs!(b, rotor)
rms_b = b_skipped ? NaN : norm(b) / sqrt(length(b))
println("$rung Dirichlet case: $(rotor.ncells) panels, rms(b) = $rms_b " *
        "(direct source assembly $(round(t_rhs_assembly; digits=1)) s)\n")

solution_column = 2
"Restore the frozen pre-solve state (cold solve; FGS seeds from body.strength)."
function reset_cold!()
    pnl.set_strengths!(rotor)                       # col 1 = BC sources
    rotor.strength[:, solution_column] .= 0
    rotor.velocity .= frozen_velocity
    rotor.potential .= potential_frozen             # rhs = -potential contract
    return nothing
end

rotor.core_size = rotor.core_size_panel   # ctor-time trees at panel offset
