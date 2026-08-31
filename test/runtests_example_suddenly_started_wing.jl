using Test
import FLOWPanel as pnl
using LinearAlgebra: norm

if !isdefined(@__MODULE__, :SSWConfig)
    old_no_plot = get(ENV, "SSW_NO_PLOT", nothing)
    ENV["SSW_NO_PLOT"] = "true"
    try
        include(joinpath(pnl.examples_path, "suddenly_started_wing.jl"))
    finally
        if old_no_plot === nothing
            delete!(ENV, "SSW_NO_PLOT")
        else
            ENV["SSW_NO_PLOT"] = old_no_plot
        end
    end
end
if !isdefined(@__MODULE__, :ssw_convert_rows)
    include(joinpath(pnl.examples_path, "ssw_representation_probe.jl"))
end
if !isdefined(@__MODULE__, :_sps_case_row)
    include(joinpath(pnl.examples_path, "ssw_sheet_particle_split.jl"))
end

@testset "Suddenly-started wing shedding configuration and tags" begin
    base = SSWConfig(wake_model=:particle, backend_kind=:direct, save_vtk=false,
        verbose=false)
    overlap_pps = _ssw_shedding_method(base)
    @test overlap_pps isa pnl.OverlapPPS
    @test overlap_pps.overlap == 1.3
    @test overlap_pps.p_per_step == 2

    sigma_pps_config = ssw_with(base; shed_method=:sigma_pps,
        sigma_over_c=0.08, pps_n=4)
    sigma_pps = _ssw_shedding_method(sigma_pps_config)
    @test sigma_pps isa pnl.SigmaPPS
    @test sigma_pps.sigma ≈ 0.08 * base.c
    @test sigma_pps.p_per_step == 4

    sigma_overlap_config = ssw_with(base; shed_method=:sigma_overlap,
        sigma_over_c=0.3, pps_overlap=1.7)
    sigma_overlap = _ssw_shedding_method(sigma_overlap_config)
    @test sigma_overlap isa pnl.SigmaOverlap
    @test sigma_overlap.sigma ≈ 0.3 * base.c
    @test sigma_overlap.overlap ≈ 1.7

    @test_throws ArgumentError _ssw_shedding_method(
        ssw_with(base; shed_method=:unknown))
    @test_throws ArgumentError _ssw_shedding_method(
        ssw_with(base; shed_method=:sigma_pps))
    @test_throws ArgumentError _ssw_shedding_method(
        ssw_with(base; pps_n=0))

    tags = Set(_ssw_case_tag(c) for c in (
        base,
        ssw_with(base; pps_n=3),
        ssw_with(base; pps_overlap=1.4),
        ssw_with(base; max_particles=100_000),
        ssw_with(base; wake_core_over_c=1e-4),
        sigma_pps_config,
        sigma_overlap_config,
    ))
    @test length(tags) == 7
    @test occursin("_core0p0001", _ssw_case_tag(
        ssw_with(base; wake_core_over_c=1e-4)))

    withenv(
        "SSW_WAKE_MODEL" => "particle",
        "SSW_PANEL_ROWS" => "32",
        "SSW_SHED_METHOD" => "sigma_pps",
        "SSW_SIGMA_OVER_C" => "0.15",
        "SSW_PPS" => "8",
        "SSW_OVERLAP" => "2.1",
        "SSW_WAKE_CORE" => "0.01",
        "SSW_MAX_PARTICLES" => "12345",
    ) do
        parsed = _ssw_config_from_env()
        @test parsed.wake_model == :particle
        @test parsed.panel_rows == 32
        @test parsed.shed_method == :sigma_pps
        @test parsed.sigma_over_c == 0.15
        @test parsed.pps_n == 8
        @test parsed.pps_overlap == 2.1
        @test parsed.wake_core_over_c == 0.01
        @test parsed.max_particles == 12_345
    end
end

@testset "Suddenly-started wing settledness classification" begin
    t = collect(0.0:0.125:2.0)
    settled = ssw_settled_stats(fill(1.0, length(t)), t)
    @test settled.settled
    @test settled.drift == 0
    @test settled.ripple == 0

    drifting = ssw_settled_stats(1 .+ 0.002 .* t, t)
    @test !drifting.settled
    rippling = ssw_settled_stats(1 .+ 0.002 .* (-1.0) .^ eachindex(t), t)
    @test !rippling.settled
end

@testset "Suddenly-started wing row conversion ledgers and field limit" begin
    config = SSWConfig(AR=2.0, n_span=2, n_airfoil=21,
        backend_kind=:direct, save_vtk=false, verbose=false)
    body = build_suddenly_started_wing(config)
    wake = pnl.PanelWake(body; nwakerows=1, include_final_filament=false,
        core_size=1e-4)
    y = range(-1.0, 1.0; length=3)
    for column in 1:3
        wake.nodes[1][:, 1, column] .= (1.0, y[column], 0.0)
        wake.nodes[1][:, 2, column] .= (1.5, y[column], 0.0)
    end
    wake.strength[1][1, 1, :] .= (0.7, 0.4)
    wake.nwakes[] = 1

    converted = ssw_convert_rows(wake, 1, 1, pnl.SigmaPPS(1e-4, 128);
        max_particles=4096)
    @test converted.family_counts[:trailing] == 256
    @test converted.family_counts[:unsteady] == 256
    @test converted.family_counts[:tip] == 256
    @test converted.actual_gamma ≈ converted.expected_gamma atol=1e-12
    @test converted.actual_impulse ≈ converted.expected_impulse atol=1e-12

    sheet = ssw_subset_wake(body, wake, 1, 1)
    probes = [SVector(0.5, 0.0, 0.4), SVector(1.25, 0.0, 0.5),
        SVector(2.0, 0.25, 0.4)]
    error = ssw_compare_fields(probes, sheet, converted.pfield)
    @test error.velocity_max < 1e-3
    @test error.gradient_max < 1e-3
end

@testset "Suddenly-started wing relaxation recorder preserves stock update" begin
    relaxation, recorder = ssw_recording_relaxation()
    pfield = pnl.FLOWVPM.ParticleField(1, Float64;
        fmm=pnl.FLOWVPM.FMM(autotune_reg_error=false),
        relaxation=relaxation)
    pnl.FLOWVPM.add_particle(pfield, zeros(3), [1.0, 0.0, 0.0], 0.1)
    pnl.FLOWVPM.get_J(pfield, 1) .= [0.0, 0.0, 0.0, 0.0, 0.0, -1.0,
        0.0, 1.0, 0.0]
    reference = deepcopy(pfield)
    pnl.FLOWVPM.relax_correctedpedrizzetti(
        pnl.FLOWVPM.relaxation_correctedpedrizzetti.rlxf, reference, 1)
    relaxation(pnl.FLOWVPM.get_particle(pfield, 1))
    @test pnl.FLOWVPM.get_Gamma(pfield, 1) ≈
        pnl.FLOWVPM.get_Gamma(reference, 1)
    @test ssw_relaxation_stats(recorder).samples == 1
end

@testset "Suddenly-started wing admissibility classification" begin
    config = SSWConfig(AR=6.0, n_span=24, dt_star=0.125, panel_rows=8)
    artifact(values) = (; mean=values)
    control = (; tail_CL=(; mean=1.0, drift=0.0, ripple=0.0, settled=true),
        tail_gamma=(; drift=0.0, ripple=0.0, settled=true),
        gamma_artifact=artifact([1.0, 0.8]),
        loading_artifact=artifact([1.0, 0.8]), config, wake_rows=160,
        tag="control")
    result = (; tail_CL=(; mean=1.004, drift=5e-4, ripple=2e-3, settled=true),
        tail_gamma=(; drift=5e-4, ripple=2e-3, settled=true),
        gamma_artifact=artifact([1.005, 0.795]),
        loading_artifact=artifact([1.003, 0.797]), config, wake_rows=8,
        n_particles=100, tag="case")
    eligible = _sps_case_row(result, control, :sigma_overlap, 0.3, 2, 1.3, 1.0)
    @test eligible.eligible
    @test eligible.admissible
    @test eligible.h_trailing_over_sigma ≈ 0.125 / 0.3
    @test eligible.h_unsteady_over_sigma ≈ 0.25 / (2 * 0.3)
    mechanism = _sps_case_row(result, control, :sigma_pps, 0.3, 2, 1.3, 1.0)
    @test !mechanism.eligible
    @test !mechanism.admissible
    @test mechanism.h_trailing_over_sigma ≈ 0.125 / (2 * 0.3)
    @test mechanism.h_unsteady_over_sigma ≈ 0.25 / (2 * 0.3)
    bad = merge(result, (; tail_CL=(; mean=1.006, drift=5e-4,
        ripple=2e-3, settled=true)))
    @test !_sps_case_row(bad, control, :sigma_overlap, 0.3, 2, 1.3, 1.0).admissible
end

@testset "Suddenly-started wing matched-geometry controller" begin
    pfield = pnl.FLOWVPM.ParticleField(2, Float64;
        fmm=pnl.FLOWVPM.FMM(autotune_reg_error=false))
    pnl.FLOWVPM.add_particle(pfield, [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0], 0.1)
    pfield.particles[pnl.FLOWVPM.U_INDEX, 1] .= [2.0, 3.0, 0.0]
    # Mimic the position immediately after the normal Euler update.
    pnl.FLOWVPM.get_X(pfield, 1) .+= 0.1 .* [2.0, 3.0, 0.0]
    controller = SSWMatchedGeometry([1.0, 0.0, 0.0])
    pnl.apply_particle_policy!(controller, pfield,
        pnl.ParticleMaintenanceContext(nothing, 1, 0.1))
    @test pnl.FLOWVPM.get_X(pfield, 1) ≈ [1.1, 0.0, 0.0]
    @test pnl.FLOWVPM.get_Gamma(pfield, 1) ≈ [0.0, 1.0, 0.0]

    # On subsequent updates the stored trajectory advances by U∞ regardless
    # of the induced velocity retained in the normal Euler step.
    pnl.FLOWVPM.get_X(pfield, 1) .+= 0.1 .* [-4.0, 2.0, 1.0]
    pnl.apply_particle_policy!(controller, pfield,
        pnl.ParticleMaintenanceContext(nothing, 2, 0.1))
    @test pnl.FLOWVPM.get_X(pfield, 1) ≈ [1.2, 0.0, 0.0]
    @test ssw_matched_geometry_stats(controller).samples == 2
end

function _ssw_boundary_edges(cells)
    counts = Dict{Tuple{Int,Int},Int}()
    for cell in eachcol(cells)
        for (a, b) in ((cell[1], cell[2]), (cell[2], cell[3]), (cell[3], cell[1]))
            edge = a < b ? (a, b) : (b, a)
            counts[edge] = get(counts, edge, 0) + 1
        end
    end
    return [edge for (edge, count) in counts if count == 1]
end

@testset "Suddenly-started wing geometry and formulation" begin
    contour = naca0012_contour(21)
    @test size(contour, 2) == 2
    @test minimum(contour[:, 1]) ≈ 0.0 atol=100eps()
    @test maximum(contour[:, 1]) ≈ 1.0 atol=100eps()
    @test maximum(contour[:, 2]) ≈ -minimum(contour[:, 2]) rtol=1e-13
    @test contour[1, 1] ≈ 1.0 atol=100eps()
    @test contour[1, 2] ≈ 0.0 atol=1e-12
    @test count(i -> isapprox(contour[i, 1], 0.0; atol=100eps()), axes(contour, 1)) == 1

    config = SSWConfig(n_span=2, n_airfoil=21, save_vtk=false,
        backend_kind=:direct, verbose=false)
    nodes, cells = suddenly_started_wing_mesh(config.c, config.AR * config.c;
        n_span=config.n_span, n_airfoil=config.n_airfoil)
    n_section = size(contour, 1)
    @test size(nodes) == (3, n_section * (config.n_span + 1))
    @test size(cells) == (3, 2 * n_section * config.n_span)
    @test collect(extrema(nodes[2, :])) ≈
        [-config.AR * config.c / 2, config.AR * config.c / 2]

    watertight, open_cells = pnl.iswatertight(nodes, cells; return_open_cells=true)
    @test !watertight
    @test !isempty(open_cells)
    boundary_edges = _ssw_boundary_edges(cells)
    @test length(boundary_edges) == 2 * n_section
    @test all(edge -> begin
        y = nodes[2, collect(edge)]
        all(isapprox.(abs.(y), config.AR * config.c / 2; atol=100eps()))
    end, boundary_edges)

    body = build_suddenly_started_wing(config)
    @test body isa pnl.RigidWakeBody{pnl.VortexRing,1,Float64,false}
    @test !pnl.has_dirichlet_bc(body)
    @test !body.watertight
    @test !body.semiinfinite_wake
    @test length(body.shedding) == 1
    @test size(only(body.shedding), 2) == config.n_span
    @test body.nsheddings == config.n_span

    fmm_config = ssw_with(config; backend_kind=:fmm,
        fmm_expansion_order=8, fmm_multipole_acceptance=0.25,
        fmm_leaf_size=64)
    backend = _ssw_backend(fmm_config)
    @test backend.expansion_order == 8
    @test backend.multipole_acceptance == 0.25
    @test backend.leaf_size == 64

    shedding = only(body.shedding)
    te_midpoints = map(eachcol(shedding)) do edge
        panel, nia, nib = edge[1:3]
        node_ids = body.cells[[nia, nib], panel]
        (body.nodes[:, node_ids[1]] + body.nodes[:, node_ids[2]]) / 2
    end
    @test all(p -> isapprox(p[1], config.c; atol=1e-8 * config.c), te_midpoints)
    @test issorted(getindex.(te_midpoints, 2)) || issorted(getindex.(te_midpoints, 2); rev=true)
end

@testset "Suddenly-started wing time, inflow, Wagner function, and wake capacity" begin
    config = SSWConfig(n_span=2, n_airfoil=21, t_end_star=0.25,
        dt_star=0.125, save_vtk=false, backend_kind=:direct, verbose=false)
    time = collect(ssw_time_range(config))
    @test time ≈ [0.0, 0.125, 0.25] .* config.c ./ config.U
    @test diff(time) ≈ fill(config.dt_star * config.c / config.U, 2)

    sim = prepare_suddenly_started_wing(config)
    @test sim.Uinf(0.0) == zero(sim.Uinf(0.0))
    post_start = sim.Uinf(time[2])
    @test norm(post_start) ≈ config.U
    @test post_start[1] ≈ config.U * cosd(config.alpha_deg)
    @test post_start[2] == 0
    @test abs(post_start[3]) ≈ config.U * sind(config.alpha_deg)

    @test ssw_wagner(0.0) ≈ 0.5
    @test ssw_wagner(1.0) ≈ 1 - 0.165exp(-0.09) - 0.335exp(-0.6)
    @test ssw_wagner(Inf) == 1.0
    @test issorted(ssw_wagner.(0.0:0.25:7.0))

    @test length(sim.t_range) == length(time)
    @test sim.t_range ≈ time
    @test sim.wake isa pnl.PanelWake
    @test size(only(sim.wake.nodes), 2) - 1 == length(time) - 1
    @test !sim.wake.include_final_filament
    @test sim.wake.core_size ≈ config.c * config.wake_core_over_c
    @test !sim.wake.overflowed[]
    @test pnl.monitor_provides(sim.pressure) == (:P,)
    @test pnl.monitor_requires(sim.force) == (:P,)
    @test pnl.audit_monitors((sim.pressure, sim.force)) !== nothing
end

@testset "Suddenly-started wing tiny unsteady smoke" begin
    config = SSWConfig(n_span=2, n_airfoil=21, t_end_star=0.125,
        dt_star=0.125, save_vtk=false, backend_kind=:direct, verbose=false,
        backslash_max_panels=10_000)
    sim = prepare_suddenly_started_wing(config)

    pnl.simulate!((sim.wing,), (sim.wake,), sim.frames, sim.maneuver!, sim.Uinf, sim.t_range;
        body_solvers=(sim.solver,),
        backend=sim.backend,
        monitors=sim.monitors,
        path=nothing,
        set_Das_eta_freestream=NaN,
        grad_mu_options=(; basis=:tri),
        verbose=false,
    )

    @test all(isfinite, sim.force.force[:, 2])
    @test norm(sim.force.force[:, 2]) > 0
    @test sim.wake.nwakes[] == 1
    @test !sim.wake.overflowed[]
end
