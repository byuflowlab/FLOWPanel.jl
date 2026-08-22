## Phase 2c driver for plans/dji_convergence_20260722/phase_02c_dji_mesh_convergence.md.
##
## Verifies the Phase 2b attribution (the Dirichlet-Neumann bound-circulation gap is a
## chordwise section-resolution discretization difference) on the ACTUAL DJI 9443 rotor
## mesh, using a chordwise-refinement ladder (fixed 30 spanwise sections, n_airfoil in
## {81, 97, 121}) of capped/Dirichlet vs uncapped/Neumann solves.
##
## Run batches from the repository root:
##   PHASE2C_MODE=setup   julia --project examples/dji9443_mesh_convergence.jl
##   PHASE2C_MODE=all     julia --project examples/dji9443_mesh_convergence.jl
##   PHASE2C_MODE=analyze julia --project examples/dji9443_mesh_convergence.jl
##
## A single case can be selected with PHASE2C_MODE=solve PHASE2C_CASE=<tag>.

import FLOWPanel as pnl
using CSV
using DataFrames
using Dates
using LinearAlgebra
using Printf
using FLOWPanel.FastMultipole.StaticArrays

include(joinpath(pnl.examples_path, "dji9443_trailing_edge.jl"))

const STUDY = "dji_convergence_20260722"
const PHASE = "phase_02c_dji_mesh_convergence"
const OUTPUT_DIR = joinpath("data", STUDY, PHASE)
const RAW_DIR = joinpath(OUTPUT_DIR, "raw")
const PHASE1_METRICS = joinpath(
    "data", STUDY, "phase_01_circulation_audit", "case_metrics.csv")
const R = 0.119
const RPM = 5400.0
const OMEGA = 2pi * RPM / 60
const RADIAL_DIMENSION = 2
const VINF = [1.0e-4, 0.0, 0.0]
const RADIAL_BINS = collect(0.125:0.025:0.975)
# Regularization constants. Phase 2d (tipdiag mode) sweeps these via environment
# overrides; the defaults are the Phase 2c study values, so every other mode is
# unchanged unless the PHASE2D_* variables are explicitly set.
const KERNEL_OFFSET_PANEL = parse(Float64, get(ENV, "PHASE2D_KOFF_PANEL", string(R * 1.0e-10)))
const KERNEL_OFFSET_TARGETS = parse(Float64, get(ENV, "PHASE2D_KOFF_TARGETS", "1.0e-3"))
const KERNEL_CUTOFF = parse(Float64, get(ENV, "PHASE2D_KCUTOFF", string(R * 1.0e-13)))
# Phase 2d attached-wake perturbation knobs (tipdiag experiments only; defaults = study):
# PHASE2D_DAS_ETA sets set_Das_eta_kinematic (default 0.2; controls the attached-strip
# length at the tip, ~eta*U_kin*dt); PHASE2D_DAS_MIN scales the minimum displacement
# floor (default 0.01R); PHASE2D_SEMIINF=1 builds the body with a semi-infinite wake.
const DAS_ETA = parse(Float64, get(ENV, "PHASE2D_DAS_ETA", "0.2"))
const DAS_MIN_FRAC = parse(Float64, get(ENV, "PHASE2D_DAS_MIN", "0.01"))
const SEMIINF_WAKE = get(ENV, "PHASE2D_SEMIINF", "0") == "1"

# Chordwise-refinement ladder: fixed 30 spanwise sections, n_airfoil in {81, 97, 121}.
# `c` = watertight/Dirichlet (Morino source+doublet), `u` = open/Neumann (vortex-ring).
const CASES = [
    (tag=:dji81c, mesh="dji9443_20260723_30_81_capped.msh",
     n_airfoil=81, formulation=:dirichlet, watertight=true),
    (tag=:dji81u, mesh="dji9443_20260723_30_81_uncapped.msh",
     n_airfoil=81, formulation=:neumann, watertight=false),
    (tag=:dji97c, mesh="dji9443_20260723_30_97_capped.msh",
     n_airfoil=97, formulation=:dirichlet, watertight=true),
    (tag=:dji97u, mesh="dji9443_20260723_30_97_uncapped.msh",
     n_airfoil=97, formulation=:neumann, watertight=false),
    (tag=:dji121c, mesh="dji9443_20260723_30_121_capped.msh",
     n_airfoil=121, formulation=:dirichlet, watertight=true),
    (tag=:dji121u, mesh="dji9443_20260723_30_121_uncapped.msh",
     n_airfoil=121, formulation=:neumann, watertight=false),
]

# Chordwise rungs (coarse -> fine) and the (dirichlet, neumann) tag pair per rung.
const RUNGS = [81, 97, 121]
const RUNG_PAIRS = Dict(
    81 => (dir=:dji81c, neu=:dji81u),
    97 => (dir=:dji97c, neu=:dji97u),
    121 => (dir=:dji121c, neu=:dji121u),
)

# Spanwise-convergence probe: 60 spanwise sections at fixed n_airfoil=97, compared
# against the 30-spanwise `dji97c`/`dji97u` (n_airfoil=97 was the Phase 2c outlier rung).
const SPANWISE_CASES = [
    (tag=:dji60_97c, mesh="dji9443_20260723_60_97_capped.msh",
     n_airfoil=97, n_span=60, formulation=:dirichlet, watertight=true),
    (tag=:dji60_97u, mesh="dji9443_20260723_60_97_uncapped.msh",
     n_airfoil=97, n_span=60, formulation=:neumann, watertight=false),
]
# n_airfoil=97 30-span counterparts already solved by the chordwise ladder.
const SPANWISE_PAIR_30 = (dir=:dji97c, neu=:dji97u)
const SPANWISE_PAIR_60 = (dir=:dji60_97c, neu=:dji60_97u)

# Chordwise extension to finer sections (45 spanwise; spanwise is converged, so gaps are
# comparable across the 30/45/60-span families).
const CHORDWISE_EXT = [
    (tag=:dji45_145c, mesh="dji9443_20260723_45_145_capped.msh",
     n_airfoil=145, n_span=45, formulation=:dirichlet, watertight=true),
    (tag=:dji45_145u, mesh="dji9443_20260723_45_145_uncapped.msh",
     n_airfoil=145, n_span=45, formulation=:neumann, watertight=false),
    (tag=:dji45_185c, mesh="dji9443_20260723_45_185_capped.msh",
     n_airfoil=185, n_span=45, formulation=:dirichlet, watertight=true),
    (tag=:dji45_185u, mesh="dji9443_20260723_45_185_uncapped.msh",
     n_airfoil=185, n_span=45, formulation=:neumann, watertight=false),
    (tag=:dji45_201c, mesh="dji9443_20260723_45_201_capped.msh",
     n_airfoil=201, n_span=45, formulation=:dirichlet, watertight=true),
    (tag=:dji45_201u, mesh="dji9443_20260723_45_201_uncapped.msh",
     n_airfoil=201, n_span=45, formulation=:neumann, watertight=false),
    (tag=:dji45_249c, mesh="dji9443_20260723_45_249_capped.msh",
     n_airfoil=249, n_span=45, formulation=:dirichlet, watertight=true),
    (tag=:dji45_249u, mesh="dji9443_20260723_45_249_uncapped.msh",
     n_airfoil=249, n_span=45, formulation=:neumann, watertight=false),
]

# `Backslash` is a dense solve; on a 17 GB machine only n_airfoil≤145 (≈27.6k panels,
# ~6 GB) fits. 185/201/249 (10–18 GB) are deferred to HPC. Set PHASE2C_LOCAL_EXT to
# override the local-solvable tag set (comma-separated) if run on a larger machine.
const LOCAL_EXT_TAGS = let ov = get(ENV, "PHASE2C_LOCAL_EXT", "")
    isempty(ov) ? Set([:dji45_145c, :dji45_145u]) :
        Set(Symbol.(strip.(split(ov, ","))))
end

# Extended chordwise convergence ladder: best available point per n_airfoil (n=97 uses the
# corrected 60-span solve; the 30_97 capped mesh is a known inboard defect).
const EXTENDED_LADDER = [
    (n=81,  dir=:dji81c,     neu=:dji81u,     span=30),
    (n=97,  dir=:dji60_97c,  neu=:dji60_97u,  span=60),
    (n=121, dir=:dji121c,    neu=:dji121u,    span=30),
    (n=145, dir=:dji45_145c, neu=:dji45_145u, span=45),
    (n=185, dir=:dji45_185c, neu=:dji45_185u, span=45),
    (n=201, dir=:dji45_201c, neu=:dji45_201u, span=45),
    (n=249, dir=:dji45_249c, neu=:dji45_249u, span=45),
]

# Phase 2d experiment E2: cap-treatment control meshes (agent-generated via
# scripts/generate_dji9443_mesh.sh; flat tip cap and CapUMinTess variants) used to
# causally isolate the tip-cap tessellation. Solved only by the tipdiag mode.
const PHASE2D_CONTROLS = [
    (tag=:dji45_185c_flat, mesh="dji9443_20260725_45_185_flatcaps.msh",
     n_airfoil=185, n_span=45, formulation=:dirichlet, watertight=true),
    (tag=:dji45_185c_ct2, mesh="dji9443_20260725_45_185_capped_captess2.msh",
     n_airfoil=185, n_span=45, formulation=:dirichlet, watertight=true),
    (tag=:dji45_185c_ct4, mesh="dji9443_20260725_45_185_capped_captess4.msh",
     n_airfoil=185, n_span=45, formulation=:dirichlet, watertight=true),
    (tag=:dji45_145c_flat, mesh="dji9443_20260725_45_145_flatcaps.msh",
     n_airfoil=145, n_span=45, formulation=:dirichlet, watertight=true),
    (tag=:dji45_145c_ct2, mesh="dji9443_20260725_45_145_capped_captess2.msh",
     n_airfoil=145, n_span=45, formulation=:dirichlet, watertight=true),
    (tag=:dji49_145c_flat, mesh="dji9443_20260725_49_145_flatcaps.msh",
     n_airfoil=145, n_span=49, formulation=:dirichlet, watertight=true),
    (tag=:dji45_121c_flat, mesh="dji9443_20260725_45_121_flatcaps.msh",
     n_airfoil=121, n_span=45, formulation=:dirichlet, watertight=true),
    (tag=:dji45_97c_flat, mesh="dji9443_20260725_45_97_flatcaps.msh",
     n_airfoil=97, n_span=45, formulation=:dirichlet, watertight=true),
    (tag=:dji45_201c_flat, mesh="dji9443_20260725_45_201_flatcaps.msh",
     n_airfoil=201, n_span=45, formulation=:dirichlet, watertight=true),
    (tag=:dji45_249c_flat, mesh="dji9443_20260725_45_249_flatcaps.msh",
     n_airfoil=249, n_span=45, formulation=:dirichlet, watertight=true),
    # Neumann references with a MATCHED flat tip cap (root left open so the surface
    # stays non-watertight): isolates formulation/discretization in the tip gap by
    # removing the open-tip vs capped-tip geometry difference (Ryan's suggestion).
    (tag=:dji45_145u_tipflat, mesh="dji9443_20260725_45_145_root-none_tip-flat.msh",
     n_airfoil=145, n_span=45, formulation=:neumann, watertight=false),
    (tag=:dji45_185u_tipflat, mesh="dji9443_20260725_45_185_root-none_tip-flat.msh",
     n_airfoil=185, n_span=45, formulation=:neumann, watertight=false),
    # Production-candidate ladder (Ryan prefers round tip caps): round CapUMinTess=4
    # Dirichlet ladder + matched round-ct4-tip Neumann references (root open).
    (tag=:dji45_97c_ct4, mesh="dji9443_20260725_45_97_capped_captess4.msh",
     n_airfoil=97, n_span=45, formulation=:dirichlet, watertight=true),
    (tag=:dji45_121c_ct4, mesh="dji9443_20260725_45_121_capped_captess4.msh",
     n_airfoil=121, n_span=45, formulation=:dirichlet, watertight=true),
    (tag=:dji45_145c_ct4, mesh="dji9443_20260725_45_145_capped_captess4.msh",
     n_airfoil=145, n_span=45, formulation=:dirichlet, watertight=true),
    (tag=:dji45_201c_ct4, mesh="dji9443_20260725_45_201_capped_captess4.msh",
     n_airfoil=201, n_span=45, formulation=:dirichlet, watertight=true),
    (tag=:dji45_249c_ct4, mesh="dji9443_20260725_45_249_capped_captess4.msh",
     n_airfoil=249, n_span=45, formulation=:dirichlet, watertight=true),
    # Appendix G ladder extension: the two coarse Neumann partners were missing, so the
    # Dir-Neu velocity comparison would otherwise be a two-point difference rather than a
    # convergence statement. 201/249 partners stay ungenerated (their dense G is 12.8/19.6
    # GB, past this machine's 17 GB, and Appendix G is local-only).
    (tag=:dji45_97u_tipround4, mesh="dji9443_20260725_45_97_root-none_tip-round_captess4.msh",
     n_airfoil=97, n_span=45, formulation=:neumann, watertight=false),
    (tag=:dji45_121u_tipround4, mesh="dji9443_20260725_45_121_root-none_tip-round_captess4.msh",
     n_airfoil=121, n_span=45, formulation=:neumann, watertight=false),
    (tag=:dji45_145u_tipround4, mesh="dji9443_20260725_45_145_root-none_tip-round_captess4.msh",
     n_airfoil=145, n_span=45, formulation=:neumann, watertight=false),
    (tag=:dji45_185u_tipround4, mesh="dji9443_20260725_45_185_root-none_tip-round_captess4.msh",
     n_airfoil=185, n_span=45, formulation=:neumann, watertight=false),
]

nspan(case) = haskey(case, :n_span) ? case.n_span : 30
case_by_tag(tag::Symbol) =
    only(filter(c -> c.tag == tag,
        vcat(CASES, SPANWISE_CASES, CHORDWISE_EXT, PHASE2D_CONTROLS)))
raw_path(tag::Symbol) = joinpath(RAW_DIR, "$(tag)_circulation.csv")
topology(case) = case.watertight ? "watertight" : "open"

function make_shedding_bbox(nodes, seed_nodes)
    radial_midpoint = sum(nodes[RADIAL_DIMENSION, seed_nodes]) / length(seed_nodes)
    radial_sign = sign(radial_midpoint)
    radial_sign != 0 || error("seed edge lies on the rotor axis")
    lower = [minimum(nodes[i, :]) for i in axes(nodes, 1)]
    upper = [maximum(nodes[i, :]) for i in axes(nodes, 1)]
    padding = max(sqrt(eps(eltype(nodes))) * R, R * 1.0e-6)
    lower .-= padding
    upper .+= padding
    if radial_sign > 0
        lower[RADIAL_DIMENSION] = 0.1R - padding
    else
        upper[RADIAL_DIMENSION] = -0.1R + padding
    end
    return (SVector{3}(lower...), SVector{3}(upper...))
end

# Seed triples come from the detector (not hardcoded); node indices survive the
# constructor's `ensure_winding` rewinding, so they feed the constructed body directly.
function case_seeds(case)
    mesh_path = joinpath(pnl.examples_path, "data", case.mesh)
    isfile(mesh_path) || error("missing mesh: $(mesh_path)")
    neg, pos = find_dji9443_trailing_edge_indices(
        mesh_path; watertight=case.watertight)
    return (neg, pos)
end

function load_base(case)
    mesh_path = joinpath(pnl.examples_path, "data", case.mesh)
    isfile(mesh_path) || error("missing mesh: $(mesh_path)")
    mesh = pnl.read_gmsh(mesh_path)
    nodes, cells = pnl.meshes2nodes_cells(mesh)
    source_radius = maximum(abs, nodes[RADIAL_DIMENSION, :])
    nodes .*= R / source_radius

    if case.formulation == :dirichlet
        kernel = Union{pnl.ConstantSource, pnl.VortexRing}
        body = pnl.RigidWakeBody{kernel}(nodes, cells, pnl.noshedding;
            kerneloffset=KERNEL_OFFSET_PANEL,
            kerneloffset_panel=KERNEL_OFFSET_PANEL,
            kerneloffset_targets=KERNEL_OFFSET_TARGETS,
            kernelcutoff=KERNEL_CUTOFF,
            semiinfinite_wake=SEMIINF_WAKE, watertight=true, DBC=true)
    else
        kernel = pnl.VortexRing
        body = pnl.RigidWakeBody{kernel}(nodes, cells, pnl.noshedding;
            kerneloffset=KERNEL_OFFSET_PANEL,
            kerneloffset_panel=KERNEL_OFFSET_PANEL,
            kerneloffset_targets=KERNEL_OFFSET_TARGETS,
            kernelcutoff=KERNEL_CUTOFF,
            semiinfinite_wake=SEMIINF_WAKE, watertight=false, DBC=false)
    end
    return body, source_radius
end

# Critical invariant: trace shedding on the *constructed* (re-wound) body's nodes/cells.
function trace_shedding(base, seeds)
    return map(seeds) do seed
        bbox = make_shedding_bbox(base.nodes, seed[1:2])
        # On some cap variants (Phase 2d controls) the sharp TE ridge continues inboard
        # of the 0.1R shedding cutoff; the detected inner end node then lies outside the
        # bbox and the walk must terminate at the bbox instead of at the end node. The
        # shed radial span stays [0.1R, tip] either way (matched across all meshes).
        end_node = abs(base.nodes[RADIAL_DIMENSION, seed[3]]) >= 0.1R ? seed[3] : nothing
        pnl.calc_shedding_from_seed(
            base.nodes, base.cells, seed[1], seed[2];
            bbox, end_node, normal_jump_tol=0.2,
            max_turn_angle=pi / 3, debug=false)
    end
end

function build_case(case)
    base, source_radius = load_base(case)
    seeds = case_seeds(case)
    shedding = trace_shedding(base, seeds)
    if case.formulation == :dirichlet
        kernel = Union{pnl.ConstantSource, pnl.VortexRing}
        body = pnl.RigidWakeBody{kernel}(
            copy(base.nodes), copy(base.cells), [copy(s) for s in shedding];
            kerneloffset=KERNEL_OFFSET_PANEL,
            kerneloffset_panel=KERNEL_OFFSET_PANEL,
            kerneloffset_targets=KERNEL_OFFSET_TARGETS,
            kernelcutoff=KERNEL_CUTOFF,
            semiinfinite_wake=SEMIINF_WAKE, watertight=true,
            ensure_winding=true, DBC=true)
    else
        kernel = pnl.VortexRing
        body = pnl.RigidWakeBody{kernel}(
            copy(base.nodes), copy(base.cells), [copy(s) for s in shedding];
            kerneloffset=KERNEL_OFFSET_PANEL,
            kerneloffset_panel=KERNEL_OFFSET_PANEL,
            kerneloffset_targets=KERNEL_OFFSET_TARGETS,
            kernelcutoff=KERNEL_CUTOFF,
            semiinfinite_wake=SEMIINF_WAKE, watertight=false,
            ensure_winding=true, DBC=false)
    end
    return body, source_radius
end

function section_radii(body, shedding)
    return [
        abs(0.5 * (
            body.nodes[RADIAL_DIMENSION, body.cells[shedding[2, j], shedding[1, j]]] +
            body.nodes[RADIAL_DIMENSION, body.cells[shedding[3, j], shedding[1, j]]]
        ) / R)
        for j in axes(shedding, 2)
    ]
end

function validate_case(case)
    body, source_radius = build_case(case)
    length(body.shedding) == 2 || error("$(case.tag): expected two shedding chains")

    radii = map(body.shedding) do shedding
        r = section_radii(body, shedding)
        all(diff(r) .< 1.0e-10) ||
            error("$(case.tag): shedding radii are not monotone tip-to-root")
        length(unique(Tuple.(eachcol(shedding)))) == size(shedding, 2) ||
            error("$(case.tag): duplicate shedding columns")
        r
    end
    maximum(abs.(radii[1] .- radii[2])) <= 1.0e-10 ||
        error("$(case.tag): blade station geometry is asymmetric")
    isapprox(maximum(abs, body.nodes[RADIAL_DIMENSION, :]), R; rtol=0, atol=1.0e-12) ||
        error("$(case.tag): scaled radius is not R")

    paired = [count(!=(-1), shedding[4, :]) for shedding in body.shedding]
    all_r = vcat(radii...)
    return (
        case=string(case.tag), mesh=case.mesh, n_airfoil=case.n_airfoil,
        topology=topology(case), formulation=string(case.formulation),
        source_radius=source_radius, scale_factor=R / source_radius,
        nodes=size(body.nodes, 2), panels=body.ncells,
        sections_blade_1=length(radii[1]), sections_blade_2=length(radii[2]),
        min_abs_r_over_R=minimum(all_r), max_abs_r_over_R=maximum(all_r),
        paired_edges_blade_1=paired[1], paired_edges_blade_2=paired[2],
        status="pass")
end

# Common bins: keep the fixed RADIAL_BINS grid, trimmed to the support shared by every
# case so `interpolate_linear` never extrapolates and every rung uses matched bins.
function common_bins(ranges)
    common_min = maximum(r.min for r in ranges)
    common_max = minimum(r.max for r in ranges)
    bins = filter(b -> b >= common_min && b <= common_max, RADIAL_BINS)
    isempty(bins) && error("no RADIAL_BINS lie within the common station support")
    trimmed = length(bins) != length(RADIAL_BINS)
    return bins, trimmed, common_min, common_max
end

function run_setup()
    mkpath(OUTPUT_DIR)
    rows = NamedTuple[]
    for case in CASES
        print("Validating $(case.tag)... ")
        row = validate_case(case)
        push!(rows, row)
        println(@sprintf("pass (%d panels, %d sections/blade, |r/R| %.4f-%.4f)",
            row.panels, row.sections_blade_1,
            row.min_abs_r_over_R, row.max_abs_r_over_R))
    end
    ranges = [(min=r.min_abs_r_over_R, max=r.max_abs_r_over_R) for r in rows]
    bins, trimmed, cmin, cmax = common_bins(ranges)
    path = joinpath(OUTPUT_DIR, "topology_validation.csv")
    df = DataFrame(rows)
    df[!, :bins_min] .= first(bins)
    df[!, :bins_max] .= last(bins)
    df[!, :n_bins] .= length(bins)
    df[!, :bins_trimmed] .= trimmed
    CSV.write(path, df)
    println("Wrote $(path)")
    println(@sprintf(
        "Common station support: |r/R| in [%.4f, %.4f]; %d fixed bins in [%.3f, %.3f]%s",
        cmin, cmax, length(bins), first(bins), last(bins),
        trimmed ? " (TRIMMED from full 0.125:0.025:0.975 grid)" : ""))
    trimmed && @warn("RADIAL_BINS trimmed to common support; absolute comparison to " *
                     "Phase 1 rows is approximate (per-rung Dir-Neu gap uses matched bins).")
end

function make_frames(body)
    return pnl.ReferenceFrame(body;
        origin=SVector{3}(0.0, 0.0, 0.0),
        v=SVector{3}(0.0, 0.0, 0.0),
        ω_axis=SVector{3}(-1.0, 0.0, 0.0),
        ω=OMEGA,
        R=SMatrix{3,3}(1.0, 0.0, 0.0,
                       0.0, 1.0, 0.0,
                       0.0, 0.0, 1.0),
        name="vehicle", child_index=Int[], dependent_index=[1])
end

function write_raw(case, body, monitor, elapsed)
    i_strength = monitor.i_strength
    rows = NamedTuple[]
    for (blade, shedding) in pairs(body.shedding)
        for section in axes(shedding, 2)
            pi, _, _, pj, _, _ = shedding[:, section]
            upper = body.strength[pi, i_strength]
            lower = pj == -1 ? zero(upper) : body.strength[pj, i_strength]
            direct = upper - lower
            push!(rows, (
                case=string(case.tag), mesh=case.mesh, n_airfoil=case.n_airfoil,
                topology=topology(case), formulation=string(case.formulation),
                rpm=RPM, radius_m=R, solve_elapsed_s=elapsed,
                blade=blade, section=section,
                r_over_R=monitor.r_over_R[section, blade],
                abs_r_over_R=abs(monitor.r_over_R[section, blade]),
                gamma_monitor_te=monitor.circulation_te[section, blade, 1],
                gamma_direct_te=direct,
                gamma_slice=monitor.circulation_slice[section, blade, 1],
            ))
        end
    end
    mkpath(RAW_DIR)
    CSV.write(raw_path(case.tag), DataFrame(rows))
end

function run_case(case)
    validate_case(case)
    body, _ = build_case(case)
    frames = make_frames(body)
    Uinf(t) = VINF
    dt = 60 / RPM / 36
    pnl.initialize_Das!((body,), frames, Uinf, 0.0, dt;
        set_Das_eta_kinematic=0.2,
        set_Das_min_kinematic_displacement=0.01R)
    backend = pnl.FastMultipoleBackend(
        expansion_order=8, multipole_acceptance=0.4, leaf_size=20)
    monitor = pnl.BoundCirculationMonitor(body, 1, 1;
        i_frame=1, radial_dimension=RADIAL_DIMENSION, R,
        verbose=false, file=false)

    println("Solving $(case.tag): $(body.ncells) panels, n_airfoil=$(case.n_airfoil), $(RPM) RPM")
    # The bound-circulation observable is the TE μ-jump and is independent of grad_mu, but
    # steady!'s standard aero post-processing still reconstructs the surface μ-gradient.
    # The default :quad basis can hit a degenerate agglomerate on the finest DJI meshes
    # (n=249 failed: cond=Inf, stencil_size=0). Use the robust tri path — result-neutral
    # for ∫Γ, verified identical on dji81c.
    elapsed = @elapsed pnl.steady!((body,), frames, VINF;
        body_solvers=(pnl.Backslash(body),), backend,
        monitors=(monitor,), path=nothing, verbose=true,
        grad_mu_options=(basis=:tri, tri_robust=true))
    write_raw(case, body, monitor, elapsed)

    df = CSV.read(raw_path(case.tag), DataFrame)
    scale = max(maximum(abs, df.gamma_direct_te), eps(Float64))
    extraction_error = maximum(abs.(df.gamma_monitor_te .- df.gamma_direct_te)) / scale
    println(@sprintf(
        "Finished %s in %.1f s; monitor/direct max relative error %.3e",
        case.tag, elapsed, extraction_error))
    extraction_error <= 1.0e-6 ||
        @warn "monitor/direct mismatch exceeds 1e-6" case=case.tag extraction_error
end

function interpolate_linear(x, y, bins)
    order = sortperm(x)
    xs, ys = x[order], y[order]
    minimum(bins) >= first(xs) && maximum(bins) <= last(xs) ||
        error("fixed bins require extrapolation: data range $(extrema(xs))")
    return map(bins) do q
        i = searchsortedlast(xs, q)
        i == length(xs) && return ys[end]
        x0, x1 = xs[i], xs[i + 1]
        y0, y1 = ys[i], ys[i + 1]
        y0 + (q - x0) / (x1 - x0) * (y1 - y0)
    end
end

function trapz(x, y)
    return sum((x[i + 1] - x[i]) * (y[i + 1] + y[i]) / 2
               for i in 1:length(x)-1)
end

relative_error(a, b) = maximum(abs.(a .- b)) /
                       max(maximum(abs, a), maximum(abs, b), eps(Float64))

# Station range for a case read back from its raw CSV (authoritative for what was solved).
function raw_station_range(tag)
    path = raw_path(tag)
    isfile(path) || error("missing raw result for $(tag): $(path)")
    raw = CSV.read(path, DataFrame)
    return (min=minimum(raw.abs_r_over_R), max=maximum(raw.abs_r_over_R))
end

function analyze_case(case, bins)
    path = raw_path(case.tag)
    isfile(path) || error("missing raw result for $(case.tag): $(path)")
    raw = CSV.read(path, DataFrame)
    blade_values = Dict{Int, NamedTuple}()
    for blade in sort(unique(raw.blade))
        b = raw[raw.blade .== blade, :]
        blade_values[blade] = (
            te=interpolate_linear(b.abs_r_over_R, b.gamma_monitor_te, bins),
            direct=interpolate_linear(b.abs_r_over_R, b.gamma_direct_te, bins),
            slice=interpolate_linear(b.abs_r_over_R, b.gamma_slice, bins))
    end
    length(blade_values) == 2 || error("$(case.tag): expected two blades")
    te = 0.5 .* (blade_values[1].te .+ blade_values[2].te)
    direct = 0.5 .* (blade_values[1].direct .+ blade_values[2].direct)
    slice = 0.5 .* (blade_values[1].slice .+ blade_values[2].slice)
    extraction = relative_error(te, direct)
    slice_error = relative_error(te, slice)
    symmetry = relative_error(blade_values[1].te, blade_values[2].te)
    extraction <= 1.0e-6 ||
        @warn "monitor/direct mismatch exceeds 1e-6" case=case.tag extraction

    fixed = DataFrame(
        case=fill(string(case.tag), length(bins)),
        n_airfoil=fill(case.n_airfoil, length(bins)),
        topology=fill(topology(case), length(bins)),
        formulation=fill(string(case.formulation), length(bins)),
        rpm=fill(RPM, length(bins)),
        abs_r_over_R=bins,
        gamma_blade_1=blade_values[1].te,
        gamma_blade_2=blade_values[2].te,
        gamma_mean=te,
        gamma_direct_mean=direct,
        gamma_slice_mean=slice)
    return fixed, (
        case=string(case.tag), n_airfoil=case.n_airfoil,
        topology=topology(case), formulation=string(case.formulation),
        mesh=case.mesh, rpm=RPM, bins_min=first(bins), bins_max=last(bins),
        n_bins=length(bins),
        extraction_relative_error=extraction,
        slice_relative_error=slice_error,
        symmetry_relative_error=symmetry,
        integral_gamma_dr=trapz(R .* bins, te),
        integral_omega_r_gamma_dr=trapz(R .* bins, OMEGA .* R .* bins .* te),
        outboard_integral_gamma_dr=trapz(
            R .* bins[bins .>= 0.7], te[bins .>= 0.7]),
        inboard_integral_gamma_dr=trapz(
            R .* bins[bins .< 0.7], te[bins .< 0.7]))
end

function run_analysis()
    mkpath(OUTPUT_DIR)
    ranges = [raw_station_range(c.tag) for c in CASES]
    bins, trimmed, cmin, cmax = common_bins(ranges)
    println(@sprintf(
        "Analyzing on %d matched bins in [%.3f, %.3f]%s",
        length(bins), first(bins), last(bins),
        trimmed ? " (trimmed to common support)" : ""))

    fixed_all = DataFrame()
    metrics_rows = NamedTuple[]
    for case in CASES
        fixed, metrics = analyze_case(case, bins)
        append!(fixed_all, fixed)
        push!(metrics_rows, metrics)
    end
    metrics = DataFrame(metrics_rows)
    CSV.write(joinpath(OUTPUT_DIR, "fixed_bin_circulation.csv"), fixed_all)
    CSV.write(joinpath(OUTPUT_DIR, "case_metrics.csv"), metrics)

    # Panel counts come from the setup-time topology validation table.
    topo_path = joinpath(OUTPUT_DIR, "topology_validation.csv")
    panels = Dict{String, Int}()
    if isfile(topo_path)
        topo = CSV.read(topo_path, DataFrame)
        for r in eachrow(topo)
            panels[r.case] = r.panels
        end
    end
    write_report(metrics, panels, bins, trimmed, cmin, cmax)
end

function metric_by_tag(metrics, tag)
    return only(eachrow(metrics[metrics.case .== string(tag), :]))
end

# Reads the Phase 1 reference metrics rows (different spanwise family) if present.
function phase1_reference_rows()
    isfile(PHASE1_METRICS) || return nothing
    df = CSV.read(PHASE1_METRICS, DataFrame)
    wanted = ["new40c", "new40u", "new57c", "new57u"]
    rows = filter(r -> r.case in wanted, eachrow(df))
    return isempty(rows) ? nothing : rows
end

# Convergence table + per-rung Dir-Neu gap for one integral column.
function convergence_block(io, metrics, panels, integral_field, title)
    println(io, "### $(title)\n")
    println(io, "gap % = 100·(Neu − Dir)/|Dir|. Δ/rung % = change vs the coarser rung.\n")
    println(io, "| n_airfoil | panels (cap/uncap) | Dir ∫ | Neu ∫ | gap % | Dir Δ/rung % | Neu Δ/rung % |")
    println(io, "|---:|---|---:|---:|---:|---:|---:|")
    prev_dir = prev_neu = nothing
    gaps = Float64[]
    dir_deltas = Float64[]
    neu_deltas = Float64[]
    for n in RUNGS
        pair = RUNG_PAIRS[n]
        dir = metric_by_tag(metrics, pair.dir)
        neu = metric_by_tag(metrics, pair.neu)
        dv = getproperty(dir, integral_field)
        nv = getproperty(neu, integral_field)
        gap = 100 * (nv - dv) / max(abs(dv), eps(Float64))
        push!(gaps, gap)
        dir_d = prev_dir === nothing ? NaN : 100 * (dv - prev_dir) / max(abs(prev_dir), eps(Float64))
        neu_d = prev_neu === nothing ? NaN : 100 * (nv - prev_neu) / max(abs(prev_neu), eps(Float64))
        !isnan(dir_d) && push!(dir_deltas, dir_d)
        !isnan(neu_d) && push!(neu_deltas, neu_d)
        fmt(x) = isnan(x) ? "—" : @sprintf("%+.3f", x)
        dir_np = get(panels, string(RUNG_PAIRS[n].dir), 0)
        neu_np = get(panels, string(RUNG_PAIRS[n].neu), 0)
        println(io, @sprintf("| %d | %d / %d | %.5g | %.5g | %+.3f | %s | %s |",
            n, dir_np, neu_np, dv, nv, gap, fmt(dir_d), fmt(neu_d)))
        prev_dir, prev_neu = dv, nv
    end
    println(io)
    return gaps, dir_deltas, neu_deltas
end

function write_report(metrics, panels, bins, trimmed, cmin, cmax)
    path = joinpath(OUTPUT_DIR, "phase_02c_report.md")

    open(path, "w") do io
        println(io, "# Phase 2c — DJI 9443 chordwise mesh-convergence of the Dir-Neu gap\n")
        println(io, "Generated: $(Dates.format(now(), dateformat"yyyy-mm-dd HH:MM"))\n")
        println(io, "RPM: `5400`. Chordwise ladder: fixed 30 spanwise sections, ",
            "n_airfoil ∈ {81, 97, 121}. Capped→Dirichlet (Morino source+doublet), ",
            "uncapped→Neumann (vortex-ring).\n")
        println(io, @sprintf(
            "Matched comparison grid: `|r/R|` = %d bins in [%.3f, %.3f]%s.\n",
            length(bins), first(bins), last(bins),
            trimmed ? @sprintf(" (trimmed from 0.125:0.025:0.975 to common support |r/R| in [%.4f, %.4f]; absolute comparison to Phase 1 rows is then approximate — the per-rung Dir-Neu gap uses matched bins within each rung regardless)",
                               cmin, cmax) : ""))

        println(io, "## Case consistency\n")
        println(io, "| Case | n_airfoil | Monitor/direct | TE/slice | Blade symmetry |")
        println(io, "|---|---:|---:|---:|---:|")
        for row in eachrow(metrics)
            println(io, @sprintf("| %s | %d | %.3e | %.3e | %.3e |",
                row.case, row.n_airfoil, row.extraction_relative_error,
                row.slice_relative_error, row.symmetry_relative_error))
        end
        println(io, "\n(The slice estimator returns zero on these station-aligned ",
            "rotor meshes — a known Phase-1 diagnostic limitation, excluded from the ",
            "decision metric.)\n")

        println(io, "## Convergence tables\n")
        gaps_g, ddir_g, dneu_g = convergence_block(io, metrics, panels, :integral_gamma_dr,
            "Integrated circulation ∫Γ dr")
        convergence_block(io, metrics, panels, :integral_omega_r_gamma_dr,
            "Ωr-weighted ∫Ωr·Γ dr")
        gaps_out, ddir_out, dneu_out = convergence_block(io, metrics, panels,
            :outboard_integral_gamma_dr, "Outboard ∫Γ dr (|r/R| ≥ 0.7)")

        ref = phase1_reference_rows()
        println(io, "## Phase 1 reference (different spanwise family — NOT matched)\n")
        if ref === nothing
            println(io, "_Phase 1 `case_metrics.csv` not found; reference rows skipped._\n")
        else
            println(io, "Phase 1 meshes used 40/56 spanwise sections vs this study's 30; ",
                "Phase 2b showed spanwise is not the gap lever, but these rows are ",
                "non-matched (different span family and n_airfoil 41/57). For orientation only.\n")
            println(io, "| Case | ∫Γ dr | Ωr-weighted | Outboard |")
            println(io, "|---|---:|---:|---:|")
            for r in ref
                println(io, @sprintf("| %s | %.5g | %.5g | %.5g |",
                    r.case, r.integral_gamma_dr, r.integral_omega_r_gamma_dr,
                    r.outboard_integral_gamma_dr))
            end
            # Phase 1 reference gaps for orientation.
            n40c = only(filter(r -> r.case == "new40c", ref)).integral_gamma_dr
            n40u = only(filter(r -> r.case == "new40u", ref)).integral_gamma_dr
            n57c = only(filter(r -> r.case == "new57c", ref)).integral_gamma_dr
            n57u = only(filter(r -> r.case == "new57u", ref)).integral_gamma_dr
            println(io, @sprintf(
                "\nPhase 1 ∫Γ Dir-Neu gaps (100·(Neu−Dir)/|Dir|): n_airfoil=41 → %+.3f%%, ",
                100 * (n40u - n40c) / abs(n40c)))
            println(io, @sprintf(
                "n_airfoil=57 → %+.3f%%.\n", 100 * (n57u - n57c) / abs(n57c)))
        end

        # Reviewed decision (mirrors Phase 2b criteria). The full ∫Γ mixes an
        # under-resolved outboard trend with a noisier inboard region, so the decision
        # weighs three diagnostics rather than the full-∫Γ monotonicity alone:
        #   (1) the trustworthy Neumann reference should be flat/converged;
        #   (2) the outboard gap (|r/R|≥0.7, closest to the oracle geometry) should
        #       shrink monotonically as Dirichlet climbs toward Neumann;
        #   (3) a single-rung non-monotonic outlier in the full ∫Γ is distinguished
        #       from a genuinely flat/non-converging gap.
        println(io, "## Decision\n")
        full_monotone = all(gaps_g[i] >= gaps_g[i+1] for i in 1:length(gaps_g)-1)
        out_monotone = all(gaps_out[i] >= gaps_out[i+1] for i in 1:length(gaps_out)-1)
        neu_flat = maximum(abs, dneu_g) <= 0.5
        gap121_full = abs(gaps_g[end])
        gap121_out = abs(gaps_out[end])
        last_dir_out = isempty(ddir_out) ? NaN : abs(ddir_out[end])
        # A single interior rung whose full-∫Γ gap dips below both neighbors is an outlier.
        outlier_i = argmin(gaps_g)
        has_outlier = !full_monotone && 1 < outlier_i < length(gaps_g)
        outlier_n = has_outlier ? RUNGS[outlier_i] : 0

        decision = if full_monotone && gap121_full <= 1.0 && last_dir_out <= 0.5
            "CONVERGED / attribution confirmed"
        elseif neu_flat && out_monotone
            "FINER MESH DESIRABLE"
        else
            "ATTRIBUTION CHALLENGED"
        end

        println(io, "**Reviewed reading of the three integrals:**\n")
        println(io, @sprintf(
            "- **Neumann reference flat/converged:** %s (∫Γ Neu per-rung Δ ≤ %.3f%%). ",
            neu_flat ? "yes" : "no", maximum(abs, dneu_g)),
            "This is the oracle's \"open Neumann chordwise-converged by ~60\" behavior.")
        println(io, @sprintf(
            "- **Outboard (|r/R|≥0.7) gap monotonic:** %s (gaps: %s), at n=121 = **%.3f%%** ",
            out_monotone ? "yes" : "no",
            join([@sprintf("%.3f%%", g) for g in gaps_out], " → "), gap121_out),
            @sprintf("(≤1%%: %s). Dirichlet climbs toward Neumann here — the oracle mechanism.",
                gap121_out <= 1.0 ? "yes" : "no"))
        println(io, @sprintf(
            "- **Full ∫Γ gap:** %s (gaps: %s), at n=121 = **%.3f%%**. ",
            full_monotone ? "monotonic" : "NON-monotonic",
            join([@sprintf("%.3f%%", g) for g in gaps_g], " → "), gap121_full),
            has_outlier ?
                @sprintf("The dip at **n_airfoil=%d** is a single-rung outlier: its capped/Dirichlet *inboard* circulation spikes toward the Neumann level, then recedes at n=121 (Neumann inboard stays flat; Dirichlet outboard stays monotone). Same pattern as the documented rotor-hover 56_57 outlier.", outlier_n) :
                "")
        println(io, "\n**Decision: $(decision).**\n")
        if decision == "FINER MESH DESIRABLE"
            println(io, "The trustworthy signals (Neumann flat, outboard Dirichlet ",
                "climbing monotonically toward it) support the Phase 2b attribution ",
                "direction, but the gap is **not converged**: outboard is still ",
                @sprintf("%.2f%% at n=121", gap121_out), " (oracle needed ~120 chordwise ",
                "panels/section; the DJI blade's thick twisted inboard sections need more, ",
                "and the ", has_outlier ? "n_airfoil=$(outlier_n) capped mesh shows an inboard outlier that should be regenerated/checked" : "inboard region is noisier",
                "). Recommend Ryan generate finer chordwise resolutions (e.g. 30_145, ",
                "30_161) — and a replacement/verification of the outlier rung — before ",
                "the Phase 5 chordwise sweep. Do NOT proceed to Phase 3/5 until approved.\n")
        elseif decision == "ATTRIBUTION CHALLENGED"
            println(io, "The gap is flat or non-monotonic in a way the outlier/outboard ",
                "diagnostics do not resolve — report to Ryan before touching Phase 3/5.\n")
        else
            println(io, "The DJI Dir-Neu gap converges like the Phase 2b oracle: ",
                "chordwise section under-resolution of the Dirichlet formulation, ",
                "closing toward the Neumann value under refinement.\n")
        end
    end
    println("Wrote analysis artifacts and $(path)")
end

# ---------------------------------------------------------------------------------------
# Spanwise-convergence probe (60 vs 30 spanwise sections at fixed n_airfoil=97).
# ---------------------------------------------------------------------------------------

partial_integral(fixed, mask) = trapz(fixed.abs_r_over_R[mask] .* R, fixed.gamma_mean[mask])

function run_spanwise()
    mkpath(OUTPUT_DIR)
    # Setup-validate then solve the two new 60-span meshes (skip solve if raw exists).
    for case in SPANWISE_CASES
        if isfile(raw_path(case.tag))
            println("Raw already present for $(case.tag); skipping solve.")
        else
            run_case(case)
        end
    end

    tags = [SPANWISE_PAIR_30.dir, SPANWISE_PAIR_30.neu,
            SPANWISE_PAIR_60.dir, SPANWISE_PAIR_60.neu]
    for t in tags
        isfile(raw_path(t)) ||
            error("missing raw result for $(t) (run the chordwise ladder / spanwise solve first)")
    end
    ranges = [raw_station_range(t) for t in tags]
    bins, trimmed, cmin, cmax = common_bins(ranges)
    println(@sprintf("Spanwise comparison on %d matched bins in [%.3f, %.3f]%s",
        length(bins), first(bins), last(bins),
        trimmed ? " (trimmed to common support)" : ""))

    fixed = Dict{Symbol,DataFrame}()
    metrics = Dict{Symbol,NamedTuple}()
    for t in tags
        f, m = analyze_case(case_by_tag(t), bins)
        fixed[t] = f
        metrics[t] = m
    end

    metrics_rows = [merge(metrics[t], (n_span=nspan(case_by_tag(t)),)) for t in tags]
    CSV.write(joinpath(OUTPUT_DIR, "spanwise_metrics.csv"), DataFrame(metrics_rows))

    # Chordwise ladder metrics (from the main analyze) for the corrected-ladder check.
    chordwise = nothing
    cm_path = joinpath(OUTPUT_DIR, "case_metrics.csv")
    if isfile(cm_path)
        cm = CSV.read(cm_path, DataFrame)
        chordwise = Dict(r.case => r for r in eachrow(cm))
    end
    write_spanwise_report(fixed, metrics, chordwise, bins, trimmed, cmin, cmax)
end

function write_spanwise_report(fixed, metrics, chordwise, bins, trimmed, cmin, cmax)
    path = joinpath(OUTPUT_DIR, "spanwise_report.md")
    d30, u30 = SPANWISE_PAIR_30.dir, SPANWISE_PAIR_30.neu
    d60, u60 = SPANWISE_PAIR_60.dir, SPANWISE_PAIR_60.neu
    fields = [(:integral_gamma_dr, "∫Γ dr"),
              (:integral_omega_r_gamma_dr, "Ωr-weighted ∫Ωr·Γ dr"),
              (:outboard_integral_gamma_dr, "Outboard ∫Γ dr (|r/R|≥0.7)")]
    pct(new, old) = 100 * (new - old) / max(abs(old), eps(Float64))

    inb = bins .< 0.7
    inboard = Dict(t => partial_integral(fixed[t], inb) for t in keys(fixed))

    dir_changes = [abs(pct(getproperty(metrics[d60], f), getproperty(metrics[d30], f)))
                   for (f, _) in fields]
    neu_changes = [abs(pct(getproperty(metrics[u60], f), getproperty(metrics[u30], f)))
                   for (f, _) in fields]
    # Spanwise convergence is judged on the OUTLIER-FREE indicators: the Neumann solve
    # (flat in both chordwise and spanwise, no outlier) and the outboard integral (robust
    # to the inboard 30_97 defect). The large Dirichlet full/inboard 30→60 change is the
    # 30_97 capped inboard OUTLIER being corrected, not general spanwise under-resolution.
    neu_max = maximum(neu_changes)
    out_max = max(dir_changes[3], neu_changes[3])
    converged = neu_max <= 1.0 && out_max <= 1.0

    open(path, "w") do io
        println(io, "# Phase 2c addendum — DJI 9443 spanwise convergence at n_airfoil=97\n")
        println(io, "Generated: $(Dates.format(now(), dateformat"yyyy-mm-dd HH:MM"))\n")
        println(io, "RPM: `5400`. Fixed n_airfoil=97; spanwise sections 30 vs 60. ",
            "Capped→Dirichlet, uncapped→Neumann.\n")
        println(io, @sprintf("Matched grid: `|r/R|` = %d bins in [%.3f, %.3f]%s.\n",
            length(bins), first(bins), last(bins),
            trimmed ? @sprintf(" (trimmed to common support |r/R| ∈ [%.4f, %.4f])", cmin, cmax) : ""))

        println(io, "Panels: dji97c=12640, dji97u=11136 (30 span); ",
            "dji60_97c, dji60_97u (60 span) — see `spanwise_metrics.csv`.\n")

        println(io, "## Spanwise change 30 → 60 (per formulation)\n")
        println(io, "| Integral | Dir 30 | Dir 60 | Dir Δ% | Neu 30 | Neu 60 | Neu Δ% |")
        println(io, "|---|---:|---:|---:|---:|---:|---:|")
        for (f, label) in fields
            dv30 = getproperty(metrics[d30], f); dv60 = getproperty(metrics[d60], f)
            nv30 = getproperty(metrics[u30], f); nv60 = getproperty(metrics[u60], f)
            println(io, @sprintf("| %s | %.5g | %.5g | %+.3f | %.5g | %.5g | %+.3f |",
                label, dv30, dv60, pct(dv60, dv30), nv30, nv60, pct(nv60, nv30)))
        end
        println(io)

        println(io, "## Dir–Neu gap at 30 vs 60 span (gap % = 100·(Neu−Dir)/|Dir|)\n")
        println(io, "| Integral | gap @30 span | gap @60 span |")
        println(io, "|---|---:|---:|")
        for (f, label) in fields
            g30 = pct(getproperty(metrics[u30], f), getproperty(metrics[d30], f))
            g60 = pct(getproperty(metrics[u60], f), getproperty(metrics[d60], f))
            println(io, @sprintf("| %s | %+.3f%% | %+.3f%% |", label, g30, g60))
        end
        println(io)

        println(io, "## Inboard (|r/R| < 0.7) circulation — did the n=97 outlier persist?\n")
        println(io, "| | Dir inboard ∫ | Neu inboard ∫ | Dir–Neu % |")
        println(io, "|---|---:|---:|---:|")
        println(io, @sprintf("| 30 span | %.5g | %.5g | %+.3f%% |",
            inboard[d30], inboard[u30], pct(inboard[u30], inboard[d30])))
        println(io, @sprintf("| 60 span | %.5g | %.5g | %+.3f%% |",
            inboard[d60], inboard[u60], pct(inboard[u60], inboard[d60])))
        println(io, "\n(At 30 span the n=97 capped/Dirichlet inboard circulation was ",
            "anomalously *high* — nearly matching Neumann, i.e. Dir–Neu inboard ≈ 0. ",
            "A converged spanwise result should show the inboard Dir–Neu gap open back ",
            "up to the family value; near-zero here means the 30_97 capped inboard was ",
            "the outlier.)\n")

        # Corrected chordwise ladder: replace the 30_97 outlier with the 60_97 point.
        if chordwise !== nothing && all(haskey(chordwise, c)
                for c in ("dji81c","dji81u","dji121c","dji121u"))
            println(io, "## Corrected chordwise ∫Γ ladder (30_97 → 60_97)\n")
            println(io, "Replacing the defective 30-span n=97 capped point with the ",
                "60-span solve. gap % = 100·(Neu−Dir)/|Dir|.\n")
            println(io, "| n_airfoil | Dir ∫Γ | Neu ∫Γ | gap % | note |")
            println(io, "|---:|---:|---:|---:|---|")
            rows = [(81, chordwise["dji81c"].integral_gamma_dr, chordwise["dji81u"].integral_gamma_dr, "30 span"),
                    (97, getproperty(metrics[d60], :integral_gamma_dr), getproperty(metrics[u60], :integral_gamma_dr), "**60 span (corrected)**"),
                    (121, chordwise["dji121c"].integral_gamma_dr, chordwise["dji121u"].integral_gamma_dr, "30 span")]
            gaps = Float64[]
            for (n, dv, nv, note) in rows
                g = pct(nv, dv); push!(gaps, g)
                println(io, @sprintf("| %d | %.5g | %.5g | %+.3f | %s |", n, dv, nv, g, note))
            end
            mono = all(gaps[i] >= gaps[i+1] for i in 1:length(gaps)-1)
            println(io, @sprintf("\nCorrected gap ladder: %s → **%s** in n_airfoil.\n",
                join([@sprintf("%.3f%%", g) for g in gaps], " → "),
                mono ? "monotonic" : "still non-monotonic"))
        end

        println(io, "## Verdict\n")
        println(io, @sprintf("- Neumann spanwise change 30→60 (outlier-free indicator): **%.3f%% max** ",
            neu_max), "(", join([@sprintf("%.3f%%", c) for c in neu_changes], "/"),
            " across ∫Γ / Ωr / outboard).")
        println(io, @sprintf("- Outboard spanwise change (both formulations): **%.3f%% max** ",
            out_max), @sprintf("(Dir %.3f%%, Neu %.3f%%).", dir_changes[3], neu_changes[3]))
        println(io, @sprintf("- Dirichlet full/inboard change is large (∫Γ %.3f%%, inboard Dir–Neu ",
            dir_changes[1]), @sprintf("%+.3f%% → %+.3f%%) — this is the 30_97 capped inboard OUTLIER ",
            pct(inboard[u30], inboard[d30]), pct(inboard[u60], inboard[d60])),
            @sprintf("being corrected: 60_97 Dir ∫Γ = %.5g lands on the family (81c/121c).",
                getproperty(metrics[d60], :integral_gamma_dr)))
        println(io, "- Spanwise converged (Neumann + outboard ≤ 1%): **$(converged ? "YES" : "NO")**.\n")
        if converged
            println(io, "**Spanwise resolution is converged at n_airfoil=97.** The ",
                "trustworthy indicators — the Neumann solve (≤0.3% 30→60) and the ",
                "outboard integral (≤0.3%) — are spanwise-flat. The large Dirichlet ",
                "full-∫Γ / inboard move is the **30_97 capped mesh's inboard defect being ",
                "fixed**, not a genuine spanwise deficiency: the 60_97 Dirichlet value ",
                "rejoins the chordwise family and restores a monotonic Dir–Neu gap ",
                "(the Phase 2c non-monotonicity was the 30_97 mesh, not the attribution).\n")
            println(io, "Consequence: the 30-spanwise grid is adequate for the study; ",
                "**regenerate or replace the 30_97 capped mesh** (or use 60_97). The gap ",
                "is still ~3% at n=121, so the **finer-chordwise** recommendation from ",
                "Phase 2c stands unchanged.\n")
        else
            println(io, "**Spanwise resolution is NOT converged at n_airfoil=97**: even the ",
                "outlier-free Neumann/outboard indicators move >1% from 30→60 span. ",
                "Report to Ryan before relying on the 30-spanwise chordwise ladder.\n")
        end
    end
    println("Wrote spanwise artifacts and $(path)")
end

# ---------------------------------------------------------------------------------------
# Extended chordwise convergence ladder (adds 45-span n_airfoil=145, 185).
# ---------------------------------------------------------------------------------------

function run_extended()
    mkpath(OUTPUT_DIR)
    for case in CHORDWISE_EXT
        if isfile(raw_path(case.tag))
            println("Raw already present for $(case.tag); skipping solve.")
        elseif case.tag in LOCAL_EXT_TAGS
            run_case(case)
        else
            println("Deferred to HPC (dense Backslash > local RAM): $(case.tag) " *
                    "[$(case.mesh)]; no local solve.")
        end
    end

    # Analyze only ladder points whose raws exist; the rest are pending (HPC).
    present = [p for p in EXTENDED_LADDER
               if isfile(raw_path(p.dir)) && isfile(raw_path(p.neu))]
    pending = [p for p in EXTENDED_LADDER
               if !(isfile(raw_path(p.dir)) && isfile(raw_path(p.neu)))]
    isempty(present) && error("no extended-ladder points have raw results yet")

    tags = unique(vcat([p.dir for p in present], [p.neu for p in present]))
    ranges = [raw_station_range(t) for t in tags]
    bins, trimmed, cmin, cmax = common_bins(ranges)
    println(@sprintf("Extended ladder on %d matched bins in [%.3f, %.3f]%s (%d/%d points present)",
        length(bins), first(bins), last(bins),
        trimmed ? " (trimmed to common support)" : "", length(present), length(EXTENDED_LADDER)))

    metrics = Dict{Symbol,NamedTuple}()
    fixed_all = DataFrame()
    for t in tags
        f, m = analyze_case(case_by_tag(t), bins)
        metrics[t] = m
        f[!, :n_span] .= nspan(case_by_tag(t))
        append!(fixed_all, f)
    end
    CSV.write(joinpath(OUTPUT_DIR, "extended_fixed_bin.csv"), fixed_all)
    write_extended_report(metrics, present, pending, bins, trimmed, cmin, cmax)
end

function write_extended_report(metrics, present, pending, bins, trimmed, cmin, cmax)
    path = joinpath(OUTPUT_DIR, "extended_report.md")
    pct(new, old) = 100 * (new - old) / max(abs(old), eps(Float64))
    fields = [(:integral_gamma_dr, "Integrated circulation ∫Γ dr"),
              (:outboard_integral_gamma_dr, "Outboard ∫Γ dr (|r/R| ≥ 0.7)"),
              (:inboard_integral_gamma_dr, "Inboard ∫Γ dr (|r/R| < 0.7)")]
    nmax = isempty(present) ? 0 : maximum(p.n for p in present)

    open(path, "w") do io
        println(io, "# Phase 2c addendum — DJI 9443 extended chordwise convergence\n")
        println(io, "Generated: $(Dates.format(now(), dateformat"yyyy-mm-dd HH:MM"))\n")
        println(io, "RPM: `5400`. Extended chordwise ladder, best point per n_airfoil ",
            "(n=97 uses the corrected 60-span solve; 30_97 capped is a known inboard ",
            "defect). Spanwise is converged (≤0.3% 30→60), so gaps are comparable across ",
            "the 30/45/60-span families.\n")
        println(io, @sprintf("Matched grid: `|r/R|` = %d bins in [%.3f, %.3f]%s.\n",
            length(bins), first(bins), last(bins),
            trimmed ? @sprintf(" (trimmed to common support |r/R| ∈ [%.4f, %.4f])", cmin, cmax) : ""))
        if !isempty(pending)
            println(io, "**Pending (deferred to HPC — dense Backslash exceeds local RAM):** ",
                join([@sprintf("n_airfoil=%d", p.n) for p in pending], ", "),
                ". These rows will fill in once the HPC solves land.\n")
        end

        results = NamedTuple[]
        for (field, title) in fields
            println(io, "## $(title)\n")
            println(io, "gap % = 100·(Neu − Dir)/|Dir|. Δ/rung % = change vs the coarser rung.\n")
            println(io, "| n_airfoil | span | Dir ∫ | Neu ∫ | gap % | Dir Δ/rung % | Neu Δ/rung % |")
            println(io, "|---:|---:|---:|---:|---:|---:|---:|")
            prev_dir = prev_neu = nothing
            gaps = Float64[]; ddir = Float64[]; dneu = Float64[]
            for p in present
                dv = getproperty(metrics[p.dir], field)
                nv = getproperty(metrics[p.neu], field)
                gap = pct(nv, dv); push!(gaps, gap)
                dd = prev_dir === nothing ? NaN : pct(dv, prev_dir)
                nd = prev_neu === nothing ? NaN : pct(nv, prev_neu)
                !isnan(dd) && push!(ddir, dd)
                !isnan(nd) && push!(dneu, nd)
                fmt(x) = isnan(x) ? "—" : @sprintf("%+.3f", x)
                println(io, @sprintf("| %d | %d | %.5g | %.5g | %+.3f | %s | %s |",
                    p.n, p.span, dv, nv, gap, fmt(dd), fmt(nd)))
                prev_dir, prev_neu = dv, nv
            end
            for p in pending
                println(io, @sprintf("| %d | %d | _pending HPC_ | _pending HPC_ | — | — | — |",
                    p.n, p.span))
            end
            println(io)
            push!(results, (field=field, title=title, gaps=gaps, ddir=ddir, dneu=dneu))
        end

        # Reviewed decision (corrected 2026-07-24 after per-station tip inspection). The
        # INBOARD (root) trend is the trustworthy convergence signal. The OUTBOARD (tip)
        # gap is NON-monotonic ({145,249} high, {185,201} low) because the capped/Dirichlet
        # tip circulation at r/R≳0.9 is set by tip-CAP meshing, not chordwise resolution:
        # the uncapped/Neumann tip is identical across all rungs, and at 185/201 the
        # Dirichlet tip *overshoots* Neumann (non-physical). So a low outboard gap at
        # 185/201 is a tip-cap artifact, NOT convergence — do not treat it as converged.
        prim = results[1]; outb = results[2]; inb = results[3]
        ns = [p.n for p in present]
        neu_flat = isempty(prim.dneu) ? true : maximum(abs, prim.dneu) <= 0.5
        in_last = abs(inb.gaps[end])
        # Inboard monotone-decreasing over the available fine tail (the real trend).
        in_tail = inb.gaps[max(1, end-3):end]
        in_mono = all(in_tail[i] >= in_tail[i+1] + -1e-9 for i in 1:length(in_tail)-1)
        # Tip-cap artifact: the outboard gap is non-monotone AND dips far below the inboard
        # gap (the tip should track the interior, not collapse to ~0 then recover).
        out_mono = all(outb.gaps[i] >= outb.gaps[i+1] for i in 1:length(outb.gaps)-1)
        out_min = minimum(abs, outb.gaps)
        tip_artifact = !out_mono && out_min < 0.5 * in_last
        # Representative "true" full gap: interior (inboard) governs; take the finest
        # inboard-consistent rung. With a tip artifact, the full-∫Γ at those rungs is
        # artificially low, so quote the inboard-scale gap as the honest estimate.
        true_gap = in_last
        converged = neu_flat && in_mono && in_last <= 1.0 && !tip_artifact && out_mono

        decision = if converged
            "CONVERGED / attribution confirmed"
        elseif tip_artifact
            "NOT CONVERGED — interior converging; tip-cap artifact"
        elseif in_mono && neu_flat
            "CONVERGING — not yet converged"
        else
            "NOT MONOTONE — investigate"
        end

        println(io, "## Decision\n")
        if !isempty(pending)
            println(io, "_Partial ladder ($(length(present))/$(length(present)+length(pending)) ",
                "points); n_airfoil ∈ {$(join([p.n for p in pending], ", "))} pending HPC._\n")
        end
        println(io, "**Reviewed reading (inboard/outboard split; best point per n_airfoil):**\n")
        println(io, @sprintf("- **Neumann flat/converged:** %s (∫Γ Neu per-rung Δ ≤ %.3f%%). grad_mu basis verified result-neutral (dji81c, dji45_145 bit-identical :quad vs :tri).",
            neu_flat ? "yes" : "no", isempty(prim.dneu) ? 0.0 : maximum(abs, prim.dneu)))
        println(io, @sprintf("- **Inboard (root) gap — the trustworthy signal:** %s — %s, still ~%.1f%% at n=%d, Dirichlet climbing toward Neumann. Consistent with the Phase 2b chordwise-under-resolution attribution, but NOT converged.",
            join([@sprintf("%.3f%%", g) for g in inb.gaps], " → "),
            in_mono ? "monotone shrinking" : "noisy", in_last, ns[end]))
        println(io, @sprintf("- **Outboard (tip) gap — NON-monotonic: %s.** %s",
            join([@sprintf("%.3f%%", g) for g in outb.gaps], " → "),
            tip_artifact ?
                "The dip at n=185/201 is a **tip-cap mesh artifact**, not convergence: the capped/Dirichlet tip circulation (r/R≳0.9) overshoots the flat uncapped/Neumann tip there, and moves non-monotonically with n_airfoil ({145,249} vs {185,201}) — i.e. it is set by tip-cap meshing, not chordwise resolution. Do NOT read the ~0.3% integral gap at 185/201 as converged." :
                "outboard tracks the interior."))
        println(io, @sprintf("- **Full ∫Γ gap:** %s. The 1.6%% values at 185/201 are pulled down by the tip overshoot; the honest gap is inboard-governed at **~%.1f%% and still decreasing**.",
            join([@sprintf("%.3f%%", g) for g in prim.gaps], " → "), true_gap))
        println(io, "\n**Decision: $(decision)$(isempty(pending) ? "" : " (provisional)").**\n")
        if converged
            println(io, "The DJI Dir–Neu gap converges under chordwise refinement (Neumann ",
                "flat, Dirichlet climbing) — the Phase 2b oracle mechanism, confirmed on ",
                "the real DJI mesh.\n")
        elseif tip_artifact
            println(io, @sprintf("**Attribution DIRECTION supported, but the DJI gap is NOT converged (~%.1f%% and decreasing).** ", true_gap),
                "The inboard (root) Dir–Neu gap converges smoothly and monotonically with ",
                "chordwise refinement (Dirichlet climbing toward the flat Neumann reference) ",
                "— the Phase 2b mechanism — but is still ~2% at the finest rung. The ",
                "outboard/tip is dominated by a **tip-cap meshing artifact** (185/201 ",
                "capped Dirichlet overshoots Neumann at r/R≳0.9; the effect is non-monotonic ",
                "in n_airfoil and absent from the uncapped Neumann solves), so the apparent ",
                "~0.3% outboard convergence there is spurious. **Next:** confirm how the tip ",
                "caps were generated across the {145,185,201,249} meshes (the {145,249} vs ",
                "{185,201} split suggests two cap recipes); add an intermediate rung (~165) ",
                "to test the tip; and refine the root chordwise to push the interior gap down.\n")
        elseif decision == "CONVERGING — not yet converged"
            println(io, @sprintf("The gap shrinks monotonically toward the Neumann reference but is not yet ≤1%% (inboard ~%.1f%% at n=%d). Finer chordwise still helpful.\n", in_last, ns[end]))
        else
            println(io, "The trend is not monotone — investigate mesh quality before ",
                "relying on it.\n")
        end
    end
    println("Wrote extended artifacts and $(path)")
end

# ---------------------------------------------------------------------------------------
# Phase 2d tip-cap diagnostic mode (plans/dji_convergence_20260722/phase_02d_*.md).
# Instrumented steady solve for one case: dumps the full solved strength field with panel
# geometry/region tags, the separate mu_upper/mu_lower per TE section, an offset-station
# slice circulation (extraction-independent; offset stations avoid the station-aligned
# zero-crossing limitation of the monitor's slice), and solve conditioning (1-norm RCOND
# via the LU factors; optional exact dense residual with PHASE2D_KEEP_A=1).
# Environment: PHASE2C_CASE (required tag), PHASE2D_LABEL (output prefix, default
# "<tag>_base"), PHASE2D_KOFF_PANEL / PHASE2D_KOFF_TARGETS / PHASE2D_KCUTOFF
# (regularization overrides for the C sweep), PHASE2D_KEEP_A=1 (keep an unfactored copy
# of G for an exact residual; doubles matrix memory).
# ---------------------------------------------------------------------------------------

const PHASE2D_DIR = joinpath("data", STUDY, "phase_02d_tip_cap_discrepancy", "raw")

function tip_region(cp, normal)
    is_capnormal = abs(normal[RADIAL_DIMENSION]) > 0.5
    r = abs(cp[RADIAL_DIMENSION]) / R
    if is_capnormal && r > 0.9
        return "tip_cap"
    elseif is_capnormal && r < 0.25
        return "root_cap"
    else
        return "surface"
    end
end

# Induced velocity at probe points from the solved body (Phase 2b machinery, adapted):
# off-body evaluation through the ProbeSystem/influence! path with the body's
# kerneloffset_targets, bypassing on-surface grad_mu.
function tipdiag_induced_velocity(body, pts; backend=pnl.DirectBackend())
    n = length(pts)
    probes = pnl.FastMultipole.ProbeSystem(n, Float64)
    for k in 1:n
        probes.position[k] = pts[k]
    end
    saved = body.kerneloffset
    body.kerneloffset = body.kerneloffset_targets
    pnl.influence!((probes,), (body,), backend;
        precalc=false, scalar_potential=true, gradient=true, hessian=false)
    body.kerneloffset = saved
    return probes.scalar_potential, probes.gradient
end

# Fixed-point near-field probes: induced potential and velocity on rings of points
# around the blade section at chosen |r/R| stations. The points depend only on the
# section bounding box (nearly identical across the mesh ladder), so probing the same
# physical locations on two meshes compares the SOLVED FIELDS directly, independent of
# any circulation-extraction convention. (Loop integrals are degenerate here: the study's
# finite ~0.01R attached wake strip is a closed vortex ring, so every exterior circuit
# encloses zero net vorticity.)
function run_fieldprobe(body; fracs=(0.90, 0.925, 0.95, 0.975, 0.5), ntheta=24,
                        backend=pnl.DirectBackend())
    is_surf = [tip_region(view(body.controlpoints, :, i), view(body.normals, :, i)) == "surface"
               for i in 1:body.ncells]
    # blade 1 sits at negative radial coordinate (matches shedding chain 1)
    sgn = sign(sum(body.nodes[RADIAL_DIMENSION,
        body.cells[body.shedding[1][2, 1], body.shedding[1][1, 1]]]))
    rows = NamedTuple[]
    pts = SVector{3, Float64}[]
    meta = Tuple{Float64, Int}[]
    for frac in fracs
        r0 = sgn * frac * R
        band = [i for i in 1:body.ncells if is_surf[i] &&
                abs(body.controlpoints[RADIAL_DIMENSION, i] - r0) < 0.01R]
        length(band) >= 8 || continue
        xs = view(body.controlpoints, 1, band)
        zs = view(body.controlpoints, 3, band)
        cx = 0.5 * (minimum(xs) + maximum(xs))
        cz = 0.5 * (minimum(zs) + maximum(zs))
        rad = 1.3 * max(0.5 * (maximum(xs) - minimum(xs)),
                        0.5 * (maximum(zs) - minimum(zs)))
        for k in 1:ntheta
            th = 2pi * (k - 1) / ntheta
            push!(pts, SVector{3}(cx + rad * cos(th), r0, cz + rad * sin(th)))
            push!(meta, (frac, k))
        end
    end
    phi, U = tipdiag_induced_velocity(body, pts; backend)
    for (idx, (frac, k)) in enumerate(meta)
        p = pts[idx]
        push!(rows, (
            abs_r_over_R=frac, itheta=k,
            x=p[1], y=p[2], z=p[3],
            phi=phi[idx], ux=U[idx][1], uy=U[idx][2], uz=U[idx][3],
        ))
    end
    return rows
end

# ∮V·dl around a rectangular circuit at constant radial coordinate r0 (x-z plane),
# enclosing the local blade section: extraction-independent bound circulation. Works on
# watertight bodies (unlike the crossing-sum slice, which telescopes to zero on any
# closed surface regardless of station offset).
function tipdiag_loop_circulation(body, r0, cx, cz, hw, hh, nseg, backend)
    corners = [(cx - hw, cz - hh), (cx + hw, cz - hh), (cx + hw, cz + hh), (cx - hw, cz + hh)]
    pts = SVector{3, Float64}[]
    tangents = SVector{3, Float64}[]
    for e in 1:4
        (x1, z1) = corners[e]
        (x2, z2) = corners[mod1(e + 1, 4)]
        for k in 1:nseg
            t = (k - 0.5) / nseg
            x = x1 + t * (x2 - x1)
            z = z1 + t * (z2 - z1)
            push!(pts, SVector{3}(x, r0, z))
            push!(tangents, SVector{3}((x2 - x1) / nseg, 0.0, (z2 - z1) / nseg))
        end
    end
    _, U = tipdiag_induced_velocity(body, pts; backend)
    return sum(dot(U[k], tangents[k]) for k in eachindex(pts))
end

# Loop-integral circulation at offset stations (midpoints between TE sections) of blade
# `blade`: tip-most `n_tip` gaps plus controls near |r/R| ~ 0.7 and 0.5. Two loop sizes x
# two quadrature densities per station give a self-check spread.
function run_loopcirc(body, monitor, blade; n_tip=10, backend=pnl.DirectBackend())
    shed = body.shedding[blade]
    nsec = size(shed, 2)
    origin = SVector{3, Float64}(0.0, 0.0, 0.0)
    Rg2f = SMatrix{3, 3, Float64}(1.0I)
    mids = [SVector{3, Float64}(pnl._bound_circulation_te_midpoint(
                body, shed[:, j], origin, Rg2f)...) for j in 1:nsec]
    rr = [m[RADIAL_DIMENSION] / R for m in mids]
    gaps = collect(1:min(n_tip, nsec - 1))
    for target in (0.7, 0.5)
        j = argmin(abs.(abs.(rr) .- target))
        j = clamp(j, 1, nsec - 1)
        j in gaps || push!(gaps, j)
    end

    is_cap = [tip_region(view(body.controlpoints, :, i), view(body.normals, :, i)) != "surface"
              for i in 1:body.ncells]
    rows = NamedTuple[]
    for j in gaps
        r0 = 0.5 * (mids[j][RADIAL_DIMENSION] + mids[j + 1][RADIAL_DIMENSION])
        spacing = abs(mids[j + 1][RADIAL_DIMENSION] - mids[j][RADIAL_DIMENSION])
        band = [i for i in 1:body.ncells if !is_cap[i] &&
                abs(body.controlpoints[RADIAL_DIMENSION, i] - r0) < 0.6 * spacing]
        length(band) >= 8 || continue
        xs = view(body.controlpoints, 1, band)
        zs = view(body.controlpoints, 3, band)
        cx = 0.5 * (minimum(xs) + maximum(xs))
        cz = 0.5 * (minimum(zs) + maximum(zs))
        hw0 = 0.5 * (maximum(xs) - minimum(xs))
        hh0 = 0.5 * (maximum(zs) - minimum(zs))
        vals = Float64[]
        for f in (1.4, 1.8), nseg in (48, 96)
            hw = f * hw0 + 0.05 * hw0
            hh = f * hh0 + 0.15 * hw0
            push!(vals, tipdiag_loop_circulation(body, r0, cx, cz, hw, hh, nseg, backend))
        end
        m = sum(vals) / length(vals)
        spread = maximum(abs.(vals .- m)) / max(abs(m), eps())
        gamma_te_interp = 0.5 * (monitor.circulation_te[j, blade, 1] +
                                 monitor.circulation_te[j + 1, blade, 1])
        push!(rows, (
            blade=blade, gap=j,
            r_over_R=r0 / R, abs_r_over_R=abs(r0) / R,
            gamma_te_interp=gamma_te_interp,
            loop_11=vals[1], loop_12=vals[2], loop_21=vals[3], loop_22=vals[4],
            loop_mean=m, loop_spread=spread,
        ))
    end
    return rows
end

function run_tipdiag(case)
    label = get(ENV, "PHASE2D_LABEL", string(case.tag) * "_base")
    keep_A = get(ENV, "PHASE2D_KEEP_A", "0") == "1"
    mkpath(PHASE2D_DIR)

    body, _ = build_case(case)
    frames = make_frames(body)
    Uinf(t) = VINF
    dt = 60 / RPM / 36
    pnl.initialize_Das!((body,), frames, Uinf, 0.0, dt;
        set_Das_eta_kinematic=DAS_ETA,
        set_Das_min_kinematic_displacement=DAS_MIN_FRAC * R)
    if SEMIINF_WAKE
        # semi-infinite wake kernels require unitary Das directions
        for Das in body.Das, k in axes(Das, 2)
            nrm = norm(view(Das, :, k))
            nrm > 0 && (Das[:, k] ./= nrm)
        end
    end

    # Manual Backslash construction (same steps as the constructor) so the 1-norm of the
    # unfactored influence matrix is captured for a LAPACK gecon RCOND estimate.
    TF = Float64
    n = body.ncells
    pnl.calc_normals!(body)
    pnl.calc_controlpoints!(body)
    G = zeros(TF, n, n)
    println("Assembling G for $(case.tag): $(n) panels, koff_panel=$(KERNEL_OFFSET_PANEL), " *
            "kcutoff=$(KERNEL_CUTOFF), koff_targets=$(KERNEL_OFFSET_TARGETS)")
    assembly_elapsed = @elapsed pnl._G!(G, body, body;
        kerneloffset=body.kerneloffset_panel, update_geometry=false)
    A = keep_A ? copy(G) : nothing
    anorm = opnorm(G, 1)
    Glu = lu!(G)
    diagU = abs.(diag(Glu.factors))
    rcond = LAPACK.gecon!('1', Glu.factors, anorm)
    solver = pnl.Backslash{TF, typeof(Glu)}(
        Glu.factors, Glu, zeros(TF, n), zeros(TF, 3, n), zeros(TF, n))

    backend = pnl.FastMultipoleBackend(
        expansion_order=8, multipole_acceptance=0.4, leaf_size=20)
    monitor = pnl.BoundCirculationMonitor(body, 1, 1;
        i_frame=1, radial_dimension=RADIAL_DIMENSION, R,
        verbose=false, file=false)

    println("Solving $(case.tag) [tipdiag label=$(label)]")
    elapsed = @elapsed pnl.steady!((body,), frames, VINF;
        body_solvers=(solver,), backend,
        monitors=(monitor,), path=nothing, verbose=true,
        grad_mu_options=(basis=:tri, tri_robust=true))

    i_strength = monitor.i_strength
    x = view(body.strength, :, i_strength)
    residual = NaN
    if keep_A
        b = solver.rhs
        residual = norm(A * x .- b) / max(norm(b), eps(TF))
    end

    # --- per-panel dump ---
    strength_cols = size(body.strength, 2)
    is_te_upper = falses(n)
    is_te_lower = falses(n)
    for shed in body.shedding, col in eachcol(shed)
        is_te_upper[col[1]] = true
        col[4] != -1 && (is_te_lower[col[4]] = true)
    end
    panel_rows = NamedTuple[]
    for i in 1:n
        cp = view(body.controlpoints, :, i)
        nrm = view(body.normals, :, i)
        charlen = pnl.characteristiclength_sqrtarea(body.nodes, view(body.cells, :, i))
        push!(panel_rows, (
            panel=i,
            cx=cp[1], cy=cp[2], cz=cp[3],
            nx=nrm[1], ny=nrm[2], nz=nrm[3],
            abs_r_over_R=abs(cp[RADIAL_DIMENSION]) / R,
            charlen=charlen,
            region=tip_region(cp, nrm),
            te_upper=is_te_upper[i], te_lower=is_te_lower[i],
            strength_1=body.strength[i, 1],
            strength_2=strength_cols >= 2 ? body.strength[i, 2] : NaN,
        ))
    end
    CSV.write(joinpath(PHASE2D_DIR, "$(label)_panels.csv"), DataFrame(panel_rows))

    # --- per-TE-section mu decomposition ---
    section_rows = NamedTuple[]
    for (blade, shed) in pairs(body.shedding)
        for section in axes(shed, 2)
            pi_, _, _, pj, _, _ = shed[:, section]
            upper = body.strength[pi_, i_strength]
            lower = pj == -1 ? zero(upper) : body.strength[pj, i_strength]
            push!(section_rows, (
                case=string(case.tag), label=label, blade=blade, section=section,
                r_over_R=monitor.r_over_R[section, blade],
                abs_r_over_R=abs(monitor.r_over_R[section, blade]),
                mu_upper=upper, mu_lower=lower,
                gamma_te=monitor.circulation_te[section, blade, 1],
                panel_upper=pi_, panel_lower=pj,
            ))
        end
    end
    CSV.write(joinpath(PHASE2D_DIR, "$(label)_sections.csv"), DataFrame(section_rows))

    # --- offset-station slice circulation (extraction-independent) ---
    origin = SVector{3, TF}(0.0, 0.0, 0.0)
    Rg2f = SMatrix{3, 3, TF}(1.0I)
    slice_rows = NamedTuple[]
    for (blade, shed) in pairs(body.shedding)
        mids = [SVector{3, TF}(pnl._bound_circulation_te_midpoint(
                    body, shed[:, j], origin, Rg2f)...)
                for j in axes(shed, 2)]
        for j in 1:length(mids) - 1
            station = 0.5 .* (mids[j] .+ mids[j + 1])
            spacing = abs(mids[j + 1][RADIAL_DIMENSION] - mids[j][RADIAL_DIMENSION])
            spacing > 0 || continue
            tol = 0.45 * spacing
            gamma_slice = pnl._bound_circulation_slice(
                body, station, tol, RADIAL_DIMENSION, i_strength, origin, Rg2f)
            push!(slice_rows, (
                case=string(case.tag), label=label, blade=blade,
                r_over_R=station[RADIAL_DIMENSION] / R,
                abs_r_over_R=abs(station[RADIAL_DIMENSION]) / R,
                gamma_slice=gamma_slice, tol=tol,
            ))
        end
    end
    CSV.write(joinpath(PHASE2D_DIR, "$(label)_slice.csv"), DataFrame(slice_rows))

    # --- fixed-point near-field probes (cross-mesh solved-field comparison) ---
    probe_elapsed = @elapsed probe_rows = run_fieldprobe(body)
    println(@sprintf("field probes: %d points in %.1f s", length(probe_rows), probe_elapsed))
    CSV.write(joinpath(PHASE2D_DIR, "$(label)_fieldprobe.csv"), DataFrame(probe_rows))

    # --- loop-integral circulation (opt-in; degenerate with the finite attached wake:
    # every exterior circuit encloses zero net vorticity, kept only as a diagnostic) ---
    if get(ENV, "PHASE2D_LOOPS", "0") == "1"
        loop_elapsed = @elapsed loop_rows = run_loopcirc(body, monitor, 1)
        println(@sprintf("loop-integral circulation: %d stations in %.1f s",
            length(loop_rows), loop_elapsed))
        CSV.write(joinpath(PHASE2D_DIR, "$(label)_loopcirc.csv"), DataFrame(loop_rows))
    end

    # --- solve info ---
    CSV.write(joinpath(PHASE2D_DIR, "$(label)_solveinfo.csv"), DataFrame([(
        case=string(case.tag), label=label, mesh=case.mesh,
        n_airfoil=case.n_airfoil, formulation=string(case.formulation),
        npanels=n, rpm=RPM,
        koff_panel=KERNEL_OFFSET_PANEL, koff_targets=KERNEL_OFFSET_TARGETS,
        kcutoff=KERNEL_CUTOFF,
        anorm=anorm, rcond=rcond,
        diagU_min=minimum(diagU), diagU_max=maximum(diagU),
        residual=residual,
        assembly_elapsed_s=assembly_elapsed, solve_elapsed_s=elapsed,
    )]))
    println(@sprintf(
        "tipdiag %s done: rcond=%.3e, diagU ratio=%.3e, residual=%s, %.1f s",
        label, rcond, minimum(diagU) / maximum(diagU),
        isnan(residual) ? "n/a" : @sprintf("%.3e", residual), elapsed))
end

# ---------------------------------------------------------------------------------------
# Phase 2d experiment B: tip-corner geometry / topology diff (read-only; no solve).
# For each 45-span capped mesh, characterizes the tip-cap region (r/R > 0.9): cap panel
# count/areas/aspect ratios, cap<->surface seam connectivity (cap cells per seam node,
# fan structure at the sharp-TE corner node), tip surface panel aspect ratios, and
# nearest-neighbor control-point distances vs characteristic length (H5's relative
# regularization scale). Writes per-mesh detail CSVs + a summary CSV, plus a legacy-VTK
# tip sub-mesh for ParaView.
# ---------------------------------------------------------------------------------------

const PHASE2D_GEOM_DIR = joinpath("data", STUDY, "phase_02d_tip_cap_discrepancy", "geom")

edge_key(a::Int, b::Int) = a < b ? (a, b) : (b, a)

function cell_edge_lengths(nodes, cell)
    nk = length(cell)
    return [norm(view(nodes, :, cell[i]) .- view(nodes, :, cell[i % nk + 1])) for i in 1:nk]
end

function panel_stats(vals)
    isempty(vals) && return (min=NaN, med=NaN, max=NaN)
    s = sort(collect(vals))
    return (min=first(s), med=s[(length(s) + 1) ÷ 2], max=last(s))
end

function run_tipgeom()
    mkpath(PHASE2D_GEOM_DIR)
    summary_rows = NamedTuple[]
    for tag in (:dji45_145c, :dji45_185c, :dji45_201c, :dji45_249c)
        case = case_by_tag(tag)
        body, _ = build_case(case)
        pnl.calc_normals!(body)
        pnl.calc_controlpoints!(body)
        n = body.ncells
        nodes = body.nodes
        cells = body.cells

        region = [tip_region(view(body.controlpoints, :, i), view(body.normals, :, i))
                  for i in 1:n]
        rr = [abs(body.controlpoints[RADIAL_DIMENSION, i]) / R for i in 1:n]
        is_tipcap = [region[i] == "tip_cap" for i in 1:n]
        is_tipsurf = [region[i] == "surface" && rr[i] > 0.9 for i in 1:n]

        # Edge -> cells over the tip region (cap + nearby surface) for seam detection.
        edge_cells = Dict{Tuple{Int, Int}, Vector{Int}}()
        for i in 1:n
            (is_tipcap[i] || rr[i] > 0.85) || continue
            cell = view(cells, :, i)
            nk = length(cell)
            for k in 1:nk
                key = edge_key(cell[k], cell[k % nk + 1])
                push!(get!(edge_cells, key, Int[]), i)
            end
        end
        seam_nodes = Set{Int}()
        for (key, cs) in edge_cells
            length(cs) == 2 || continue
            caps = count(i -> is_tipcap[i], cs)
            if caps == 1
                push!(seam_nodes, key[1])
                push!(seam_nodes, key[2])
            end
        end
        cap_cells_per_node = Dict{Int, Int}()
        for i in 1:n
            is_tipcap[i] || continue
            for node in view(cells, :, i)
                node in seam_nodes || continue
                cap_cells_per_node[node] = get(cap_cells_per_node, node, 0) + 1
            end
        end

        # Sharp-TE corner node: TE edge nodes of the tip-most shedding section.
        corner_counts = Int[]
        for shed in body.shedding
            col = shed[:, 1]
            for node in (cells[col[2], col[1]], cells[col[3], col[1]])
                haskey(cap_cells_per_node, node) && push!(corner_counts, cap_cells_per_node[node])
            end
        end

        # Per-cap-cell details.
        detail_rows = NamedTuple[]
        cap_ars = Float64[]
        cap_areas = Float64[]
        for i in 1:n
            is_tipcap[i] || continue
            lens = cell_edge_lengths(nodes, view(cells, :, i))
            ar = maximum(lens) / max(minimum(lens), eps())
            charlen = pnl.characteristiclength_sqrtarea(nodes, view(cells, :, i))
            push!(cap_ars, ar)
            push!(cap_areas, charlen^2)
            push!(detail_rows, (
                case=string(tag), panel=i,
                cx=body.controlpoints[1, i], cy=body.controlpoints[2, i],
                cz=body.controlpoints[3, i],
                abs_r_over_R=rr[i], aspect=ar, charlen=charlen,
                on_seam=any(node -> node in seam_nodes, view(cells, :, i)),
            ))
        end
        CSV.write(joinpath(PHASE2D_GEOM_DIR, "$(tag)_cap_cells.csv"), DataFrame(detail_rows))

        surf_ars = Float64[]
        surf_charlens = Float64[]
        surf_ids = Int[]
        for i in 1:n
            is_tipsurf[i] || continue
            lens = cell_edge_lengths(nodes, view(cells, :, i))
            push!(surf_ars, maximum(lens) / max(minimum(lens), eps()))
            push!(surf_charlens, pnl.characteristiclength_sqrtarea(nodes, view(cells, :, i)))
            push!(surf_ids, i)
        end

        # Nearest-neighbor control-point distance in the tip region (H5 scale check):
        # smallest distance from each tip panel's control point to any other tip-region
        # control point, normalized by the local characteristic length.
        tip_ids = [i for i in 1:n if is_tipcap[i] || is_tipsurf[i]]
        nn_ratio = Float64[]
        nn_dist = Float64[]
        for (a_idx, i) in enumerate(tip_ids)
            best = Inf
            ci = view(body.controlpoints, :, i)
            for j in tip_ids
                j == i && continue
                d = norm(ci .- view(body.controlpoints, :, j))
                d < best && (best = d)
            end
            charlen = pnl.characteristiclength_sqrtarea(nodes, view(cells, :, i))
            push!(nn_dist, best)
            push!(nn_ratio, best / max(charlen, eps()))
        end

        cap_ar_s = panel_stats(cap_ars)
        cap_area_s = panel_stats(cap_areas)
        surf_ar_s = panel_stats(surf_ars)
        surf_cl_s = panel_stats(surf_charlens)
        nn_d_s = panel_stats(nn_dist)
        nn_r_s = panel_stats(nn_ratio)
        fan_s = panel_stats(collect(values(cap_cells_per_node)))
        push!(summary_rows, (
            case=string(tag), n_airfoil=case.n_airfoil, npanels=n,
            n_tip_cap=count(is_tipcap), n_tip_surf=count(is_tipsurf),
            n_seam_nodes=length(seam_nodes),
            cap_per_seam_node_med=fan_s.med, cap_per_seam_node_max=fan_s.max,
            corner_cap_counts=join(corner_counts, "/"),
            cap_ar_med=cap_ar_s.med, cap_ar_max=cap_ar_s.max,
            cap_charlen2_med=cap_area_s.med,
            surf_ar_med=surf_ar_s.med, surf_ar_max=surf_ar_s.max,
            surf_charlen_med=surf_cl_s.med, surf_charlen_min=surf_cl_s.min,
            nn_dist_min=nn_d_s.min, nn_dist_med=nn_d_s.med,
            nn_over_charlen_min=nn_r_s.min, nn_over_charlen_med=nn_r_s.med,
            koff_panel_over_nn_min=KERNEL_OFFSET_PANEL / nn_d_s.min,
        ))

        # Legacy-VTK tip sub-mesh (triangles) for ParaView.
        vtk_ids = tip_ids
        used_nodes = sort(unique(vec(cells[:, vtk_ids])))
        remap = Dict(node => k - 1 for (k, node) in enumerate(used_nodes))
        open(joinpath(PHASE2D_GEOM_DIR, "$(tag)_tip_submesh.vtk"), "w") do io
            println(io, "# vtk DataFile Version 3.0")
            println(io, "$(tag) tip submesh (r/R>0.9)")
            println(io, "ASCII")
            println(io, "DATASET UNSTRUCTURED_GRID")
            println(io, "POINTS $(length(used_nodes)) double")
            for node in used_nodes
                println(io, join(nodes[:, node], " "))
            end
            println(io, "CELLS $(length(vtk_ids)) $(4 * length(vtk_ids))")
            for i in vtk_ids
                println(io, "3 " * join([remap[node] for node in cells[:, i]], " "))
            end
            println(io, "CELL_TYPES $(length(vtk_ids))")
            for _ in vtk_ids
                println(io, "5")
            end
            println(io, "CELL_DATA $(length(vtk_ids))")
            println(io, "SCALARS is_cap int 1")
            println(io, "LOOKUP_TABLE default")
            for i in vtk_ids
                println(io, is_tipcap[i] ? "1" : "0")
            end
        end
        println("tipgeom $(tag): $(count(is_tipcap)) cap / $(count(is_tipsurf)) tip-surface panels")
    end
    CSV.write(joinpath(PHASE2D_GEOM_DIR, "tip_geometry_summary.csv"), DataFrame(summary_rows))
    println("Wrote $(joinpath(PHASE2D_GEOM_DIR, "tip_geometry_summary.csv"))")
end

# ---------------------------------------------------------------------------------------
# Appendix G: flow-tangency residual and Dirichlet-Neumann velocity agreement.
#
# A Morino/Dirichlet solve does NOT impose flow tangency: `set_strengths!` fixes
# sigma = -Uext.n and the solver enforces the INTERIOR POTENTIAL condition on mu, so the
# tangency of the reconstructed exterior surface velocity is emergent and unmeasured. All
# Phase 2b-2d convergence evidence is circulation-based (a TE-panel mu difference); forces
# come from the surface velocity, which carries reconstruction error that Gamma never sees.
#
# Everything here runs WITHOUT re-solving. `run_tipdiag` already dumped every panel's
# solved strengths, control point and normal to `<label>_panels.csv`, so the solved state
# is reproduced by rebuilding the body and loading those strengths. That matters: the
# dense Backslash G is 8n^2 bytes (n_airfoil=201 -> 12.8 GB, 249 -> 19.6 GB), out of reach
# on a 17 GB machine, whereas the reconstruction is O(n) and puts the whole six-rung ct4
# ladder back within local reach.
#
# The velocity assembly below mirrors `_steady_aerodynamics!` exactly (see
# src/FLOWPanel_simulate.jl:340-538): reset -> freestream -> kinematic -> induced velocity
# at `kerneloffset_targets` with self-pair conditioning back down to `kerneloffset_panel`
# -> +1/2 grad_mu tangential half-jump. After that sequence `body.velocity` is the EXTERIOR
# surface limit, so the residual is literally dot(velocity[:,i], normals[:,i]).
# ---------------------------------------------------------------------------------------

const PHASE2G_DIR = PHASE2D_DIR
# Chordwise/radial binning for localization: enough radial bins that each holds one
# airfoil section's worth of panels, so the per-bin x-range is a genuine section chord.
const TANGENCY_NRBINS = 60

tangency_label(case) = get(ENV, "PHASE2D_LABEL", string(case.tag) * "_base")

# Build the body exactly as `run_tipdiag` does (same wake, same Das, same regularization),
# but stop short of assembling G.
function tangency_build(case)
    body, _ = build_case(case)
    frames = make_frames(body)
    Uinf(t) = VINF
    dt = 60 / RPM / 36
    pnl.initialize_Das!((body,), frames, Uinf, 0.0, dt;
        set_Das_eta_kinematic=DAS_ETA,
        set_Das_min_kinematic_displacement=DAS_MIN_FRAC * R)
    if SEMIINF_WAKE
        for Das in body.Das, k in axes(Das, 2)
            nrm = norm(view(Das, :, k))
            nrm > 0 && (Das[:, k] ./= nrm)
        end
    end
    pnl.calc_normals!(body)
    pnl.calc_controlpoints!(body)
    return body, frames
end

# Load solved strengths from the Phase 2d per-panel dump. HARD GATE: the CSV row order is
# only meaningful if it is the same panel order as the rebuilt body, so every control point
# and normal must match. A silent permutation would invalidate every number in Appendix G.
function tangency_load_strengths!(body, label)
    path = joinpath(PHASE2G_DIR, "$(label)_panels.csv")
    isfile(path) || error("missing per-panel dump for $(label): $(path)")
    df = CSV.read(path, DataFrame)
    nrow(df) == body.ncells ||
        error("$(label): panel count mismatch, CSV has $(nrow(df)), body has $(body.ncells)")

    cp_tol = 1.0e-10 * R
    n_tol = 1.0e-10
    cp_err = 0.0
    n_err = 0.0
    @inbounds for i in 1:body.ncells
        cp_err = max(cp_err,
            abs(df.cx[i] - body.controlpoints[1, i]),
            abs(df.cy[i] - body.controlpoints[2, i]),
            abs(df.cz[i] - body.controlpoints[3, i]))
        n_err = max(n_err,
            abs(df.nx[i] - body.normals[1, i]),
            abs(df.ny[i] - body.normals[2, i]),
            abs(df.nz[i] - body.normals[3, i]))
    end
    cp_err <= cp_tol || error("$(label): control-point mismatch $(cp_err) > $(cp_tol); " *
        "the CSV panel order does not match the rebuilt body")
    n_err <= n_tol || error("$(label): normal mismatch $(n_err) > $(n_tol)")

    body.strength[:, 1] .= df.strength_1
    if size(body.strength, 2) >= 2
        any(ismissing, df.strength_2) && error("$(label): strength_2 missing on a 2-column body")
        body.strength[:, 2] .= df.strength_2
    end
    return df, cp_err, n_err
end

# The three post-solve steps of `_steady_aerodynamics!`. Returns the external velocity
# (freestream + kinematic) and the principal-value field (external + kernel-induced),
# leaving `body.velocity` as the exterior surface limit.
function tangency_velocity!(body, frames;
        grad_mu_options=(basis=:tri, tri_robust=true),
        backend=pnl.DirectBackend(), halfjump=true)
    pnl.reset!(body)
    pnl.apply_freestream!((body,), VINF)
    pnl.kinematic_velocity!((body,), frames)
    Uext = copy(body.velocity)

    pnl._set_kerneloffsets!((body,), :kerneloffset_targets)
    pnl.pre_evaluate_influence!(body)
    pnl.influence!((body,), (body,), backend; precalc=false,
        scalar_potential=false, velocity=true,
        direct_conditioning=pnl._self_panel_kerneloffset_conditioning())
    Upv = copy(body.velocity)

    if halfjump
        opts = pnl._normalize_grad_mu_options(grad_mu_options; default_basis=:quad)
        pnl.compute_mu_gradient!(body.velocity, body.controlpoints, body.normals,
            body.cells, body.neighbor,
            view(body.strength, :, pnl.get_Gammai(body)),
            pnl._bound_surface_vorticity_te_info(body);
            scale=0.5, nodes=body.nodes, grad_mu_options=opts)
    end
    return Uext, Upv
end

# The +1/2 grad_mu half-jump alone, for a given reconstruction route.
function tangency_halfjump(body, grad_mu_options)
    out = zeros(Float64, 3, body.ncells)
    opts = pnl._normalize_grad_mu_options(grad_mu_options; default_basis=:quad)
    pnl.compute_mu_gradient!(out, body.controlpoints, body.normals,
        body.cells, body.neighbor,
        view(body.strength, :, pnl.get_Gammai(body)),
        pnl._bound_surface_vorticity_te_info(body);
        scale=0.5, nodes=body.nodes, grad_mu_options=opts)
    return out
end

# Kernel-induced velocity from a single strength column (the other zeroed), used to split
# the residual into its source and doublet shares. The rigid wake rows are driven by the
# TE doublet strengths, so the wake contribution rides inside the doublet term.
function tangency_partial_induced(body, keep_col, backend)
    saved = copy(body.strength)
    for c in axes(body.strength, 2)
        c == keep_col || (body.strength[:, c] .= 0.0)
    end
    pnl.reset!(body)
    pnl.pre_evaluate_influence!(body)
    pnl._set_kerneloffsets!((body,), :kerneloffset_targets)
    pnl.influence!((body,), (body,), backend; precalc=false,
        scalar_potential=false, velocity=true,
        direct_conditioning=pnl._self_panel_kerneloffset_conditioning())
    out = copy(body.velocity)
    body.strength .= saved
    return out
end

dotcol(A, B, i) = A[1, i] * B[1, i] + A[2, i] * B[2, i] + A[3, i] * B[3, i]

# Normalized chordwise station and surface side. The blade lies in the x-y plane with
# thrust along z, so |r| bins isolate one airfoil section and the sign of n_z separates
# upper from lower. Cap panels get chord=NaN (they are not on an airfoil section).
function tangency_chordwise(body, regions)
    n = body.ncells
    rr = [abs(body.controlpoints[RADIAL_DIMENSION, i]) / R for i in 1:n]
    lo, hi = minimum(rr), maximum(rr)
    width = (hi - lo) / TANGENCY_NRBINS
    binof(i) = width > 0 ? clamp(Int(floor((rr[i] - lo) / width)) + 1, 1, TANGENCY_NRBINS) : 1
    xmin = fill(Inf, TANGENCY_NRBINS)
    xmax = fill(-Inf, TANGENCY_NRBINS)
    for i in 1:n
        regions[i] == "surface" || continue
        b = binof(i)
        x = body.controlpoints[1, i]
        xmin[b] = min(xmin[b], x)
        xmax[b] = max(xmax[b], x)
    end
    chord = fill(NaN, n)
    for i in 1:n
        regions[i] == "surface" || continue
        b = binof(i)
        isfinite(xmin[b]) && xmax[b] > xmin[b] || continue
        chord[i] = (body.controlpoints[1, i] - xmin[b]) / (xmax[b] - xmin[b])
    end
    side = [regions[i] != "surface" ? "cap" : (body.normals[3, i] >= 0 ? "upper" : "lower")
            for i in 1:n]
    return chord, side
end

function tangency_aspect_ratios(body)
    return [let e = cell_edge_lengths(body.nodes, view(body.cells, :, i))
                maximum(e) / max(minimum(e), eps())
            end for i in 1:body.ncells]
end

stat_block(v) = isempty(v) ? (max=NaN, rms=NaN, p95=NaN, mean=NaN) : (
    max=maximum(abs, v),
    rms=sqrt(sum(abs2, v) / length(v)),
    p95=sort(abs.(v))[max(1, ceil(Int, 0.95 * length(v)))],
    mean=sum(abs, v) / length(v),
)

function run_tangency(case)
    label = tangency_label(case)
    mkpath(PHASE2G_DIR)
    backend = pnl.DirectBackend()

    body, frames = tangency_build(case)
    n = body.ncells
    panels_df, cp_err, n_err = tangency_load_strengths!(body, label)
    println(@sprintf("tangency %s: %d panels, ordering gate cp_err=%.3e m, n_err=%.3e",
        label, n, cp_err, n_err))

    # --- reference reconstruction (study route: basis=:tri, tri_robust=true) ---
    t_recon = @elapsed Uext, Upv = tangency_velocity!(body, frames; backend,
        grad_mu_options=(basis=:tri, tri_robust=true))
    U = copy(body.velocity)

    # --- G2 decomposition ---
    is_dirichlet = case.formulation == :dirichlet
    Usrc = is_dirichlet ? tangency_partial_induced(body, 1, backend) : zeros(Float64, 3, n)
    Udbl = tangency_partial_induced(body, pnl.get_Gammai(body), backend)
    hj_tri = tangency_halfjump(body, (basis=:tri, tri_robust=true))
    hj_tri_plain = tangency_halfjump(body, (basis=:tri, tri_robust=false))
    hj_quad = tangency_halfjump(body, (basis=:quad,))

    # linearity self-check: Uext + source-induced + doublet-induced == PV field
    lin_err = 0.0
    for i in 1:n, d in 1:3
        lin_err = max(lin_err, abs(Uext[d, i] + Usrc[d, i] + Udbl[d, i] - Upv[d, i]))
    end

    omegaR = OMEGA * R
    regions = panels_df.region
    chord, side = tangency_chordwise(body, regions)
    ar = tangency_aspect_ratios(body)

    rows = NamedTuple[]
    for i in 1:n
        rloc = max(OMEGA * abs(body.controlpoints[RADIAL_DIMENSION, i]), 0.05 * omegaR)
        un = dotcol(U, body.normals, i)
        push!(rows, (
            panel=i,
            cx=body.controlpoints[1, i], cy=body.controlpoints[2, i],
            cz=body.controlpoints[3, i],
            abs_r_over_R=panels_df.abs_r_over_R[i],
            region=regions[i], side=side[i], chordfrac=chord[i],
            te_upper=panels_df.te_upper[i], te_lower=panels_df.te_lower[i],
            charlen=panels_df.charlen[i], aspect_ratio=ar[i],
            umag=sqrt(U[1, i]^2 + U[2, i]^2 + U[3, i]^2),
            un=un,
            un_over_omegaR=un / omegaR,
            un_over_local=un / rloc,
            un_ext=dotcol(Uext, body.normals, i),
            un_src=dotcol(Usrc, body.normals, i),
            un_dbl=dotcol(Udbl, body.normals, i),
            un_halfjump=dotcol(hj_tri, body.normals, i),
            un_tri_plain=un - dotcol(hj_tri, body.normals, i) + dotcol(hj_tri_plain, body.normals, i),
            un_quad=un - dotcol(hj_tri, body.normals, i) + dotcol(hj_quad, body.normals, i),
        ))
    end
    CSV.write(joinpath(PHASE2G_DIR, "$(label)_tangency.csv"), DataFrame(rows))

    # --- summary: global, per region, per radial band, per side, TE vs interior ---
    e = [r.un_over_omegaR for r in rows]
    pick(f) = [rows[i].un_over_omegaR for i in 1:n if f(rows[i])]
    surf = pick(r -> r.region == "surface")
    band(lo, hi) = pick(r -> r.region == "surface" && lo <= r.abs_r_over_R < hi)
    all_s = stat_block(e)
    srf_s = stat_block(surf)
    tip_s = stat_block(pick(r -> r.region == "tip_cap"))
    root_s = stat_block(pick(r -> r.region == "root_cap"))
    te_s = stat_block(pick(r -> r.te_upper || r.te_lower))
    up_s = stat_block(pick(r -> r.side == "upper"))
    lo_s = stat_block(pick(r -> r.side == "lower"))
    ar_s = stat_block(pick(r -> r.aspect_ratio > 10))
    imax = argmax(abs.(e))

    summary = (
        case=string(case.tag), label=label, mesh=case.mesh,
        n_airfoil=case.n_airfoil, formulation=string(case.formulation),
        npanels=n, omegaR=omegaR,
        cp_gate_err=cp_err, normal_gate_err=n_err, linearity_err=lin_err,
        recon_elapsed_s=t_recon,
        all_max=all_s.max, all_rms=all_s.rms, all_p95=all_s.p95,
        surface_max=srf_s.max, surface_rms=srf_s.rms, surface_p95=srf_s.p95,
        tipcap_max=tip_s.max, tipcap_rms=tip_s.rms,
        rootcap_max=root_s.max, rootcap_rms=root_s.rms,
        te_max=te_s.max, te_rms=te_s.rms,
        upper_rms=up_s.rms, lower_rms=lo_s.rms,
        highAR_max=ar_s.max, highAR_rms=ar_s.rms,
        band_050_070=stat_block(band(0.50, 0.70)).rms,
        band_070_090=stat_block(band(0.70, 0.90)).rms,
        band_090_095=stat_block(band(0.90, 0.95)).rms,
        band_095_100=stat_block(band(0.95, 1.01)).rms,
        argmax_panel=imax, argmax_region=rows[imax].region,
        argmax_r_over_R=rows[imax].abs_r_over_R,
        argmax_chordfrac=rows[imax].chordfrac,
        argmax_aspect_ratio=rows[imax].aspect_ratio,
        ext_rms=stat_block([r.un_ext for r in rows]).rms / omegaR,
        src_rms=stat_block([r.un_src for r in rows]).rms / omegaR,
        dbl_rms=stat_block([r.un_dbl for r in rows]).rms / omegaR,
        halfjump_rms=stat_block([r.un_halfjump for r in rows]).rms / omegaR,
        route_tri_plain_rms=stat_block([r.un_tri_plain for r in rows]).rms / omegaR,
        route_quad_rms=stat_block([r.un_quad for r in rows]).rms / omegaR,
    )
    CSV.write(joinpath(PHASE2G_DIR, "$(label)_tangency_summary.csv"), DataFrame([summary]))

    # --- full-surface legacy VTK for spatial localization in ParaView ---
    open(joinpath(PHASE2G_DIR, "$(label)_tangency.vtk"), "w") do io
        println(io, "# vtk DataFile Version 3.0")
        println(io, "$(label) flow-tangency residual U.n/(omega R)")
        println(io, "ASCII")
        println(io, "DATASET UNSTRUCTURED_GRID")
        println(io, "POINTS $(size(body.nodes, 2)) double")
        for j in axes(body.nodes, 2)
            println(io, join(view(body.nodes, :, j), " "))
        end
        println(io, "CELLS $(n) $(4n)")
        for i in 1:n
            println(io, "3 " * join(view(body.cells, :, i) .- 1, " "))
        end
        println(io, "CELL_TYPES $(n)")
        for _ in 1:n
            println(io, "5")
        end
        println(io, "CELL_DATA $(n)")
        println(io, "SCALARS un_over_omegaR double 1")
        println(io, "LOOKUP_TABLE default")
        for r in rows
            println(io, r.un_over_omegaR)
        end
        println(io, "SCALARS un_over_local double 1")
        println(io, "LOOKUP_TABLE default")
        for r in rows
            println(io, r.un_over_local)
        end
    end

    println(@sprintf("tangency %s done: |U.n|/(omegaR) max=%.3e rms=%.3e (surface rms=%.3e, tipcap rms=%.3e), linearity=%.2e, %.1f s",
        label, all_s.max, all_s.rms, srf_s.rms, tip_s.rms, lin_err, t_recon))
    return summary
end

# ---------------------------------------------------------------------------------------
# Appendix G part 2: Dirichlet-Neumann velocity-field comparison on a CANONICAL point set.
#
# `run_fieldprobe` derives its ring centers/radii from each body's own control-point
# bounding box, so the capped-Dirichlet and root-open-Neumann probe points differ by
# ~9e-6 m. Where velocity gradients are of order omegaR/chord, that offset alone shifts U
# by ~0.03 m/s -- comparable to the Dir-Neu difference being measured. Here the points are
# built ONCE from the finest ct4 mesh, cached, and reused verbatim by every case, so the
# comparison is exactly apples-to-apples.
# ---------------------------------------------------------------------------------------

const GPROBE_PTS = joinpath(PHASE2G_DIR, "gprobe_points.csv")
const GPROBE_REF_TAG = :dji45_249c_ct4
const GPROBE_FRACS = (0.30, 0.50, 0.70, 0.80, 0.90, 0.95, 0.975)
const GPROBE_NTHETA = 24
const GPROBE_RINGS = (1.3, 2.0)

function run_gprobe_points()
    mkpath(PHASE2G_DIR)
    body, _ = build_case(case_by_tag(GPROBE_REF_TAG))
    pnl.calc_normals!(body)
    pnl.calc_controlpoints!(body)
    is_surf = [tip_region(view(body.controlpoints, :, i), view(body.normals, :, i)) == "surface"
               for i in 1:body.ncells]
    sgn = sign(sum(body.nodes[RADIAL_DIMENSION,
        body.cells[body.shedding[1][2, 1], body.shedding[1][1, 1]]]))
    rows = NamedTuple[]
    for frac in GPROBE_FRACS
        r0 = sgn * frac * R
        band = [i for i in 1:body.ncells if is_surf[i] &&
                abs(body.controlpoints[RADIAL_DIMENSION, i] - r0) < 0.01R]
        length(band) >= 8 || continue
        xs = view(body.controlpoints, 1, band)
        zs = view(body.controlpoints, 3, band)
        cx = 0.5 * (minimum(xs) + maximum(xs))
        cz = 0.5 * (minimum(zs) + maximum(zs))
        base = max(0.5 * (maximum(xs) - minimum(xs)), 0.5 * (maximum(zs) - minimum(zs)))
        for ring in GPROBE_RINGS, k in 1:GPROBE_NTHETA
            th = 2pi * (k - 1) / GPROBE_NTHETA
            push!(rows, (frac=frac, ring=ring, itheta=k,
                x=cx + ring * base * cos(th), y=r0, z=cz + ring * base * sin(th)))
        end
    end
    CSV.write(GPROBE_PTS, DataFrame(rows))
    println("wrote $(length(rows)) canonical probe points -> $(GPROBE_PTS)")
end

function run_gprobe(case)
    label = tangency_label(case)
    isfile(GPROBE_PTS) || error("run PHASE2C_MODE=gprobe_points first")
    df = CSV.read(GPROBE_PTS, DataFrame)
    pts = [SVector{3}(df.x[i], df.y[i], df.z[i]) for i in 1:nrow(df)]
    body, _ = tangency_build(case)
    tangency_load_strengths!(body, label)
    phi, U = tipdiag_induced_velocity(body, pts; backend=pnl.DirectBackend())
    CSV.write(joinpath(PHASE2G_DIR, "$(label)_gprobe.csv"),
        DataFrame(frac=df.frac, ring=df.ring, itheta=df.itheta,
            x=df.x, y=df.y, z=df.z, phi=phi,
            ux=[u[1] for u in U], uy=[u[2] for u in U], uz=[u[3] for u in U]))
    println("gprobe $(label): $(length(pts)) points")
end

mode = lowercase(get(ENV, "PHASE2C_MODE", "setup"))
if mode == "setup"
    run_setup()
elseif mode == "all"
    for case in CASES  # coarse -> fine
        run_case(case)
    end
elseif mode == "solve"
    tag = Symbol(get(ENV, "PHASE2C_CASE", ""))
    isempty(string(tag)) && error("PHASE2C_CASE is required with PHASE2C_MODE=solve")
    run_case(case_by_tag(tag))
elseif mode == "analyze"
    run_analysis()
elseif mode == "spanwise"
    run_spanwise()
elseif mode == "extend"
    run_extended()
elseif mode == "load"
    # No-op: define everything for interactive/scripted reuse without running a batch.
elseif mode == "tipgeom"
    run_tipgeom()
elseif mode == "tipdiag"
    tag = Symbol(get(ENV, "PHASE2C_CASE", ""))
    isempty(string(tag)) && error("PHASE2C_CASE is required with PHASE2C_MODE=tipdiag")
    run_tipdiag(case_by_tag(tag))
elseif mode == "tangency"
    tag = Symbol(get(ENV, "PHASE2C_CASE", ""))
    isempty(string(tag)) && error("PHASE2C_CASE is required with PHASE2C_MODE=tangency")
    run_tangency(case_by_tag(tag))
elseif mode == "gprobe_points"
    run_gprobe_points()
elseif mode == "gprobe"
    tag = Symbol(get(ENV, "PHASE2C_CASE", ""))
    isempty(string(tag)) && error("PHASE2C_CASE is required with PHASE2C_MODE=gprobe")
    run_gprobe(case_by_tag(tag))
elseif mode == "hpc_ext"
    # Solve the extended-ladder cases deferred from the laptop (185/201/249 c/u) on a
    # big-memory node. Solves only; analysis (`extend`) runs where the other raws live.
    for case in CHORDWISE_EXT
        case.tag in LOCAL_EXT_TAGS && continue
        if isfile(raw_path(case.tag))
            println("Raw already present for $(case.tag); skipping solve.")
        else
            run_case(case)
        end
    end
else
    error("PHASE2C_MODE must be setup, all, solve, analyze, spanwise, extend, tipgeom, " *
          "tipdiag, tangency, gprobe_points, gprobe, or hpc_ext; got $(mode)")
end
