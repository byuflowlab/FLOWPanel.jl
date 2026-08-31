## Rotor-hover CT force-recovery audit (BRAINSTORM/001).
##
## Produces ONE steady solved rotor state per configuration and evaluates the
## thrust coefficient CT every way the codebase allows, then sweeps the cheap
## knobs. Covers the four low-effort implementation checks in
## BRAINSTORM/001_rotor_hover_low_effort_checks.md:
##
##   1. backend            : DirectBackend vs FastMultipoleBackend (same state)
##   2. force/pressure path : Bernoulli, PressureLaplace, SurfaceVorticityForce,
##                            KuttaJoukowskiForce (KJ is a lower-confidence
##                            cross-check)
##   3. Kutta condition     : correct_kuttacondition = false vs true
##   4. rigid wake          : semiinfinite_wake = false vs true
##
## Acceptance: the implementation-error track survives only if some variant
## moves CT materially or two equivalent evaluations disagree on the same
## solved state. If everything agrees and CT stays low, the gap is wake physics.
##
##   julia --project examples/rotor_hover_force_method_audit.jl

import FLOWPanel as pnl
using FLOWPanel.FastMultipole.StaticArrays
using LinearAlgebra: norm
using Printf: @sprintf

run_name  = get(ENV, "RUN_NAME", "rotor_hover_force_method_audit")
save_path = get(ENV, "SAVE_VTK", "true") == "true" ? joinpath("data", run_name) : nothing
save_path !== nothing && (isdir(save_path) || mkpath(save_path))

# ----- fixed operating point / geometry (DJI9443 40_40, matches convergence) --
msh_file = joinpath(pnl.examples_path, "data", "dji9443_new_40_40.msh")
isfile(msh_file) || error("Mesh file not found: $(msh_file)")
mesh_tag = splitext(basename(msh_file))[1]

# 0-based ParaView TE seed point IDs (see rotor_hover_pressure_comparison.jl)
te_indices_1 = [1614, 1574, 45]   .+ 1
te_indices_2 = [3324, 3284, 1755] .+ 1

magVinf = 0.0001
AOA     = 0.0
rho     = 1.179        # NASA paper
RPM     = parse(Float64, get(ENV, "RPM", "5400"))
R       = 0.119
shedding_r_over_R = 0.1

core_size_panel   = R * 1e-10
core_size_targets = parse(Float64, get(ENV, "CORE_SIZE", get(ENV, "KERNELOFFSET", "1e-3")))
kernelcutoff         = R * 1e-13
init_Das_eta_kinematic = 0.2
set_Das_min_kinematic_displacement = 0.01 * R

# DJI9443 geometry is rotated relative to the typical rotor convention.
axial_dimension  = 1
radial_dimension = 2
omega_axis       = SVector{3}(-1.0, 0.0, 0.0)
Vinf             = magVinf * [cosd(AOA), sind(AOA), 0.0]
Uinf(t)          = Vinf
omega            = 2 * pi * RPM / 60
dt               = 60 / RPM / 36   # nominal step for Das initialization

msh = pnl.read_gmsh(msh_file)
base_nodes, base_cells = pnl.meshes2nodes_cells(msh)
base_nodes .*= R / maximum(base_nodes[radial_dimension, :])

kernel = Union{pnl.ConstantSource, pnl.VortexRing}
DBC    = kernel == pnl.VortexRing ? false : true  # Union{source,ring} -> Dirichlet (matches convergence example)

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

# Mirror the convergence example exactly: derive the shedding from a noshedding
# RigidWakeBody's nodes/cells (the watertight constructor may reorient cells),
# then rebuild with the shedding + ensure_winding. Compute shedding once and
# reuse for both wake flags (geometry-only).
const BASE_ROTOR = pnl.RigidWakeBody{kernel}(base_nodes, base_cells, pnl.noshedding;
    core_size=core_size_panel, core_size_panel, core_size_targets,
    kernelcutoff, semiinfinite_wake=false, watertight=true, DBC)
const SHEDDING1 = pnl.calc_shedding_from_seed(BASE_ROTOR.nodes, BASE_ROTOR.cells,
    te_indices_1[1], te_indices_1[2];
    bbox=make_shedding_bbox(BASE_ROTOR.nodes, te_indices_1[1:2]),
    normal_jump_tol=0.2, max_turn_angle=pi/3, debug=false)
const SHEDDING2 = pnl.calc_shedding_from_seed(BASE_ROTOR.nodes, BASE_ROTOR.cells,
    te_indices_2[1], te_indices_2[2];
    bbox=make_shedding_bbox(BASE_ROTOR.nodes, te_indices_2[1:2]),
    normal_jump_tol=0.2, max_turn_angle=pi/3, debug=false)

function build_rotor_and_frames(semiinfinite_wake::Bool)
    rotor = pnl.RigidWakeBody{kernel}(copy(BASE_ROTOR.nodes), copy(BASE_ROTOR.cells),
        [copy(SHEDDING1), copy(SHEDDING2)];
        core_size=core_size_panel, core_size_panel, core_size_targets,
        kernelcutoff, semiinfinite_wake, watertight=true, ensure_winding=true, DBC)

    frames = pnl.ReferenceFrame(rotor;
        origin=SVector{3}(0.0, 0.0, 0.0), v=SVector{3}(0.0, 0.0, 0.0),
        ω_axis=omega_axis, ω=omega,
        R=SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
        name="vehicle", child_index=Int[], dependent_index=[1])

    pnl.initialize_Das!((rotor,), frames, Uinf, 0.0, dt;
        set_Das_eta_kinematic=init_Das_eta_kinematic,
        set_Das_min_kinematic_displacement)

    # Finite wake uses Das as the (short) displacement to the next wake row.
    # The semi-infinite kernel scales by SEMIINFINITE_LENGTH internally and
    # requires unit-norm direction columns, so renormalize for that case.
    if semiinfinite_wake
        for Das in rotor.Das
            for j in axes(Das, 2)
                n = norm(@view Das[:, j])
                n > 0 && (Das[:, j] ./= n)
            end
        end
    end
    return rotor, frames
end

# ----- one clean solve -> CT from ONE force-recovery method -------------------
# Each method gets its own steady! solve. The linear solve (γ) is identical
# across methods, so this compares force recovery on the same solved state
# WITHOUT cross-monitor interference (e.g. PressureLaplace flipping
# needs_velocity_gradient or methods mutating shared body.velocity / body.P).
function solve_one(; backend, semiinfinite_wake::Bool, method::Symbol, kutta::Bool)
    rotor, frames = build_rotor_and_frames(semiinfinite_wake)
    norm_rotor = pnl.RotorNormalization(rho, 2 * R, 1)

    if method === :bernoulli
        pressure = pnl.PressureBernoulli(rho;
            unsteady=false, correct_kuttacondition=kutta, backend=backend, file=false)
        force = pnl.ForceMonitor(1, 1; i_frame=1, normalization=norm_rotor,
            correct_kuttacondition=kutta, verbose=false, file=false)
        monitors, force_monitor = (pressure, force), force
    elseif method === :laplace
        pressure = pnl.PressureLaplace((rotor,), rho; unsteady=false, file=false)
        force = pnl.ForceMonitor(1, 1; i_frame=1, normalization=norm_rotor,
            correct_kuttacondition=false, verbose=false, file=false)
        monitors, force_monitor = (pressure, force), force
    elseif method === :vorticity
        force = pnl.SurfaceVorticityForce(rotor, 1, 1; rho=rho, i_frame=1,
            normalization=norm_rotor, correct_kuttacondition=kutta,
            verbose=false, file=false)
        monitors, force_monitor = (force,), force
    elseif method === :kj
        force = pnl.KuttaJoukowskiForce(rotor, 1, 1; rho=rho, backend=backend,
            i_frame=1, normalization=norm_rotor, verbose=false, file=false)
        monitors, force_monitor = (force,), force
    else
        error("unknown method $(method)")
    end

    bname = backend isa pnl.DirectBackend ? "direct" : "fmm"
    tag = "$(mesh_tag)_$(bname)_$(semiinfinite_wake ? "semiinf" : "finite")_$(method)_$(kutta ? "kutta" : "nokutta")"
    pnl.steady!((rotor,), frames, Vinf;
        body_solvers=(pnl.Backslash(rotor),),
        backend=backend, monitors=monitors,
        path=save_path, name=tag, verbose=false)

    return -force_monitor.force[axial_dimension, 1]
end

# ----- sweep ------------------------------------------------------------------
fmm = pnl.FastMultipoleBackend(; expansion_order=8, multipole_acceptance=0.4, leaf_size=20)

println("\nRotor-hover CT force-recovery audit — $(mesh_tag), $(RPM) RPM")
println("Reference CT ≈ 0.072 (exp), 0.068 (BEM); particle-wake ≈ 0.062; prior steady rigid-wake Bernoulli ≈ 0.0505")

# QUICK mode: one Bernoulli/FMM/finite solve to validate the harness (≈0.0505).
if get(ENV, "AUDIT_QUICK", "false") == "true"
    ct = solve_one(; backend=fmm, semiinfinite_wake=false, method=:bernoulli, kutta=false)
    println("\n[QUICK] fmm finite bernoulli kutta=off  CT = $(round(ct, digits=5))  (expect ≈ 0.0505)")
    exit()
end

# kutta knob only affects Bernoulli and SurfaceVorticity; Laplace/KJ run once.
function method_kuttas(method)
    method in (:bernoulli, :vorticity) ? (false, true) : (false,)
end

results = NamedTuple[]
# Primary FMM sweep over wake type and method (Check 2, 3, 4).
for semiinf in (false, true), method in (:bernoulli, :laplace, :vorticity, :kj), kutta in method_kuttas(method)
    ct = solve_one(; backend=fmm, semiinfinite_wake=semiinf, method=method, kutta=kutta)
    push!(results, (backend=:fmm, wake=semiinf ? "semiinf" : "finite",
                    method=method, kutta=kutta, ct=ct))
end
# Backend cross-check (Check 1): Direct vs FMM at finite wake, kutta off, all methods.
for method in (:bernoulli, :laplace, :vorticity, :kj)
    ct = solve_one(; backend=pnl.DirectBackend(), semiinfinite_wake=false, method=method, kutta=false)
    push!(results, (backend=:direct, wake="finite", method=method, kutta=false, ct=ct))
end

# ----- report -----------------------------------------------------------------
ksym(k) = k ? "on" : "off"
println("\nCT by (backend, wake, method, kutta):")
println(@sprintf("  %-7s %-8s %-10s %-6s  %s", "backend", "wake", "method", "kutta", "CT"))
println("  " * "-"^45)
for r in results
    println(@sprintf("  %-7s %-8s %-10s %-6s  %.5f",
        string(r.backend), r.wake, string(r.method), ksym(r.kutta), r.ct))
end

reldiff(a, b) = abs(a - b) / max(abs(a), abs(b), eps())
ctval(; backend, wake, method, kutta) = first(r.ct for r in results
    if r.backend==backend && r.wake==wake && r.method==method && r.kutta==kutta)

# Check 2: force/pressure-method spread on the same solved state (FMM, finite, kutta off).
println("\nCheck 2 — force-method spread on the same solved state (fmm, finite, kutta=off):")
b = ctval(backend=:fmm, wake="finite", method=:bernoulli, kutta=false)
l = ctval(backend=:fmm, wake="finite", method=:laplace,   kutta=false)
v = ctval(backend=:fmm, wake="finite", method=:vorticity, kutta=false)
k = ctval(backend=:fmm, wake="finite", method=:kj,        kutta=false)
println(@sprintf("  bernoulli=%.5f  laplace=%.5f  vorticity=%.5f  kj(*)=%.5f", b, l, v, k))
println(@sprintf("  max spread (bern/lapl/vort) = %.1f%%", 100*maximum(reldiff(x,y) for x in (b,l,v) for y in (b,l,v))))

# Check 3: Kutta-condition sensitivity.
println("\nCheck 3 — correct_kuttacondition sensitivity (fmm, finite):")
for method in (:bernoulli, :vorticity)
    off = ctval(backend=:fmm, wake="finite", method=method, kutta=false)
    on  = ctval(backend=:fmm, wake="finite", method=method, kutta=true)
    println(@sprintf("  %-10s off=%.5f  on=%.5f  Δ=%.1f%%", string(method), off, on, 100*reldiff(off,on)))
end

# Check 4: finite vs semi-infinite rigid wake.
println("\nCheck 4 — finite vs semi-infinite rigid wake (fmm, kutta=off):")
for method in (:bernoulli, :laplace, :vorticity, :kj)
    fin = ctval(backend=:fmm, wake="finite",  method=method, kutta=false)
    sem = ctval(backend=:fmm, wake="semiinf", method=method, kutta=false)
    println(@sprintf("  %-10s finite=%.5f  semiinf=%.5f  Δ=%.1f%%", string(method), fin, sem, 100*reldiff(fin,sem)))
end

# Check 1: backend agreement on the same solved state.
println("\nCheck 1 — Direct vs FMM relative diff (finite, kutta=off):")
for method in (:bernoulli, :laplace, :vorticity, :kj)
    f = ctval(backend=:fmm,    wake="finite", method=method, kutta=false)
    d = ctval(backend=:direct, wake="finite", method=method, kutta=false)
    println(@sprintf("  %-10s fmm=%.5f  direct=%.5f  Δ=%.3f%%", string(method), f, d, 100*reldiff(f,d)))
end

# ----- CSV --------------------------------------------------------------------
if save_path !== nothing
    csv_path = joinpath(save_path, "audit_table.csv")
    open(csv_path, "w") do io
        println(io, "backend,wake,method,kutta,CT")
        for r in results
            println(io, @sprintf("%s,%s,%s,%s,%.8f",
                string(r.backend), r.wake, string(r.method), ksym(r.kutta), r.ct))
        end
    end
    println("\nWrote $(csv_path)")
end
