using LinearAlgebra: cross, dot, norm

if !isdefined(@__MODULE__, :pitching_wing_mesh)
    include(joinpath(pnl.examples_path, "pitching_wing.jl"))
end

@testset "Pitching-wing section interpolation bins" begin
    bin_eta = [0.75, 0.25, 0.5]
    section_eta = [0.4, 0.5, 0.1, 0.9]

    @test _section_interpolation_bins(bin_eta, section_eta) == [
        (0.25, 0.5),
        (0.5, 0.75),
        (0.25, 0.25),
        (0.75, 0.75),
    ]

    io = IOBuffer()
    _print_section_interpolation_bins(io, bin_eta, section_eta)
    output = String(take!(io))
    @test occursin("requested η=0.4000, interpolation bins η=[0.2500, 0.5000]", output)
    @test occursin("requested η=0.5000, interpolation bins η=[0.5000, 0.7500]", output)
    @test occursin("requested η=0.1000, interpolation bins η=[0.2500, 0.2500]", output)
    @test occursin("requested η=0.9000, interpolation bins η=[0.7500, 0.7500]", output)
end

@testset "Pitching-wing restart configuration and histories" begin
    sim = prepare_pitching_wing(;
        n_cycles=0.01,
        n_span=1,
        n_airfoil=21,
        n_endcap=5,
        include_static_polar=false,
        save_vtk=false,
        backend=pnl.DirectBackend(),
    )
    @test sim.setup.c_per_dt == 0.05
    @test sim.setup.das_chord_fraction == 0.05
    @test sim.setup.dt * sim.setup.U / sim.setup.c ≈ 0.05
    @test all(Das -> all(col -> norm(col) ≈ 0.05 * sim.setup.c, eachcol(Das)),
        sim.wing.Das)

    default_policy = only(sim.wake.particle_maintenance.trim_policies)
    @test default_policy isa pnl.FrameBox
    @test default_policy.i_frame == 1
    @test default_policy.xmin == SVector(-Inf, -Inf, -Inf)
    @test default_policy.xmax == SVector(
        sim.setup.c * (1 - sim.setup.pivot_chord_fraction) + 2 * sim.setup.b,
        Inf, Inf)
    @test sim.setup.wake_length_spans == 2.0
    @test sim.setup.wake_length == 2 * sim.setup.b

    custom = prepare_pitching_wing(;
        n_cycles=0.01, n_span=1, n_airfoil=21, n_endcap=5,
        include_static_polar=false, save_vtk=false, backend=pnl.DirectBackend(),
        wake_length_spans=3.5)
    custom_policy = only(custom.wake.particle_maintenance.trim_policies)
    @test custom_policy.xmax[1] ≈
        custom.setup.c * (1 - custom.setup.pivot_chord_fraction) + 3.5 * custom.setup.b
    @test custom.setup.wake_length ≈ 3.5 * custom.setup.b

    untrimmed = prepare_pitching_wing(;
        n_cycles=0.01, n_span=1, n_airfoil=21, n_endcap=5,
        include_static_polar=false, save_vtk=false, backend=pnl.DirectBackend(),
        wake_length_spans=nothing)
    @test isempty(untrimmed.wake.particle_maintenance.trim_policies)
    @test isnothing(untrimmed.setup.wake_length)
    @test isnothing(untrimmed.setup.wake_downstream_boundary)
    untrimmed_config = _pitching_wing_config(untrimmed.setup)
    @test untrimmed_config["wake"]["wake_length_spans"] == "disabled"
    @test untrimmed_config["wake"]["wake_length_m"] == "disabled"

    for invalid in (0.0, -1.0, NaN, Inf, -Inf)
        @test_throws ArgumentError prepare_pitching_wing(;
            n_cycles=0.01, n_span=1, n_airfoil=21, n_endcap=5,
            include_static_polar=false, save_vtk=false,
            backend=pnl.DirectBackend(), wake_length_spans=invalid)
    end

    static_body = build_pitching_wing_body(1.0, 2.0;
        n_span=1, n_airfoil=21, n_endcap=5, semiinfinite_wake=true)
    set_wake_Das!(static_body, [2.0, 0.0, 0.0])
    @test all(Das -> all(col -> norm(col) ≈ 1.0, eachcol(Das)), static_body.Das)

    @test _prompt_pitching_wing_restart(IOBuffer("y\n"), IOBuffer())
    @test !_prompt_pitching_wing_restart(IOBuffer("n\n"), IOBuffer())
    prompt_output = IOBuffer()
    @test _prompt_pitching_wing_restart(IOBuffer("maybe\nYES\n"), prompt_output)
    @test occursin("Please enter y or n", String(take!(prompt_output)))
    eof_error = try
        _prompt_pitching_wing_restart(IOBuffer(), IOBuffer())
        nothing
    catch err
        err
    end
    @test eof_error isa ArgumentError
    @test occursin("check_existing=false", sprint(showerror, eof_error))
    no_prompt_state = (; any=true)
    @test !_choose_pitching_wing_restart(false, no_prompt_state, IOBuffer(), IOBuffer())

    path = mktempdir()
    stale = joinpath(path, "stale.txt")
    write(stale, "old")
    _clear_pitching_wing_output!(path)
    @test isdir(path)
    @test !isfile(stale)

    config = _pitching_wing_config(sim.setup)
    @test config["wake"]["wake_length_spans"] == 2.0
    @test config["wake"]["wake_length_m"] == 2 * sim.setup.b
    @test config["wake"]["downstream_boundary_from_pivot_m"] ==
        sim.setup.wake_downstream_boundary
    config_path = _write_pitching_wing_config(path, config)
    @test config_path == joinpath(path, PITCHING_WING_CONFIG_FILE)
    @test _validate_pitching_wing_config(path, config) == TOML.parsefile(config_path)
    mismatched = deepcopy(config)
    mismatched["time"]["c_per_dt"] = 0.02
    @test_throws ArgumentError _validate_pitching_wing_config(path, mismatched)
    mismatched_wake = deepcopy(config)
    mismatched_wake["wake"]["wake_length_spans"] = 3.0
    @test_throws ArgumentError _validate_pitching_wing_config(path, mismatched_wake)
    rm(config_path)
    legacy_error = try
        _validate_pitching_wing_config(path, config)
        nothing
    catch err
        err
    end
    @test legacy_error isa ArgumentError
    @test occursin("legacy", sprint(showerror, legacy_error))

    section_eta = [0.25, 0.75]
    period = 0.5
    times = collect(0.0:0.1:0.4)
    alpha = collect(1.0:5.0)
    CL = collect(11.0:15.0)
    sections = [21.0 22 23 24 25; 31.0 32 33 34 35]
    _write_unsteady_csv(path, times, period, alpha, CL, section_eta, sections)
    history = _load_pitching_wing_history(path, times, period, section_eta, 2)
    @test history.time == times[1:3]
    @test history.alpha == alpha[1:3]
    @test history.CL == CL[1:3]
    @test history.section == sections[:, 1:3]

    merged_CL = fill(NaN, length(times))
    merged_sections = fill(NaN, length(section_eta), length(times))
    merged_CL[1:3] .= history.CL
    merged_sections[:, 1:3] .= history.section
    merged_CL[4:5] .= [16.0, 17.0]
    merged_sections[:, 4:5] .= [26.0 27.0; 36.0 37.0]
    plot_data = _pitching_wing_convergence_data(times, period, merged_CL,
        section_eta, merged_sections)
    @test plot_data.cycles == times ./ period
    @test plot_data.CL == [11.0, 12.0, 13.0, 16.0, 17.0]
    @test [series.eta for series in plot_data.sections] == section_eta
    @test plot_data.sections[1].cl == [21.0, 22.0, 23.0, 26.0, 27.0]
    @test plot_data.sections[2].cl == [31.0, 32.0, 33.0, 36.0, 37.0]

    bad_header = read(_pitching_wing_csv_path(path), String)
    write(_pitching_wing_csv_path(path), replace(bad_header, "cl_eta_0p7500" => "cl_eta_0p5000"))
    @test_throws ArgumentError _load_pitching_wing_history(path, times, period, section_eta, 2)
    _write_unsteady_csv(path, times .+ 0.01, period, alpha, CL, section_eta, sections)
    @test_throws ArgumentError _load_pitching_wing_history(path, times, period, section_eta, 2)
end

function _reflection_maps(nodes, cells; atol=1e-13)
    reflected_node = Vector{Int}(undef, size(nodes, 2))
    for i in axes(nodes, 2)
        target = (nodes[1, i], nodes[2, i], -nodes[3, i])
        reflected_node[i] = something(findfirst(j ->
            isapprox(nodes[1, j], target[1]; atol, rtol=0) &&
            isapprox(nodes[2, j], target[2]; atol, rtol=0) &&
            isapprox(nodes[3, j], target[3]; atol, rtol=0), axes(nodes, 2)), 0)
    end

    cell_lookup = Dict(
        Tuple(sort(collect(cells[:, i]))) => i for i in axes(cells, 2))
    reflected_cell = map(axes(cells, 2)) do i
        mapped_nodes = reflected_node[cells[:, i]]
        any(iszero, mapped_nodes) && return 0
        return get(cell_lookup, Tuple(sort(mapped_nodes)), 0)
    end
    return reflected_node, reflected_cell
end

function _pitching_wing_static_cl(alpha_deg, backend)
    chord = 1.0
    span = 4.0
    rho = 1.0
    body = build_pitching_wing_body(chord, span;
        n_span=2, n_airfoil=21, n_endcap=5, semiinfinite_wake=true)
    freestream = _uinf_from_alpha(1.0, alpha_deg)
    lift, _, _ = _lift_drag_span_directions(freestream)
    set_wake_Das!(body, freestream)

    pressure = pnl.PressureBernoulli(rho; correct_kuttacondition=false)
    force = pnl.ForceMonitor(1, 1;
        i_frame=-1,
        normalization=pnl.WingNormalization(rho, span * chord, chord),
        correct_kuttacondition=false,
        file=false,
    )
    pnl.steady!(body, pnl.ReferenceFrame(body), freestream;
        body_solvers=pnl.Backslash(body),
        backend,
        monitors=(pressure, force),
        path=nothing,
        name="pitching_wing_symmetry",
        verbose=false,
    )
    return dot(SVector{3}(force.force[:, 1]), lift)
end

@testset "Pitching-wing rounded end caps and Dirichlet body" begin
    chord = 1.0
    span = 4.0
    nodes, cells = pitching_wing_mesh(chord, span;
        n_span=2, n_airfoil=21, n_endcap=5)

    @test size(nodes) == (3, 114)
    @test size(cells) == (3, 224)

    reflected_node, reflected_cell = _reflection_maps(nodes, cells)
    @test all(!iszero, reflected_node)
    @test all(!iszero, reflected_cell)
    @test all(i -> reflected_node[reflected_node[i]] == i, axes(nodes, 2))
    @test all(i -> reflected_cell[reflected_cell[i]] == i, axes(cells, 2))

    @test pnl.iswatertight(nodes, cells) == (true, Int[])
    area2 = [norm(cross(
        nodes[:, cells[2, i]] - nodes[:, cells[1, i]],
        nodes[:, cells[3, i]] - nodes[:, cells[1, i]],
    )) for i in axes(cells, 2)]
    @test minimum(area2) > 100 * eps(Float64) * chord^2

    @test maximum(nodes[2, :]) > span / 2
    @test minimum(nodes[2, :]) < -span / 2
    for xcorner in (0.0, chord)
        corner_nodes = findall(
            i -> isapprox(nodes[1, i], xcorner; atol=100eps()), axes(nodes, 2))
        @test maximum(abs.(nodes[2, corner_nodes])) <= span / 2 + 100eps()
    end

    upper_cap = findall(i -> nodes[2, i] > span / 2, axes(nodes, 2))
    xmid = nodes[1, upper_cap[argmax(nodes[2, upper_cap])]]
    section = findall(i -> isapprox(nodes[1, i], xmid; atol=100eps()), axes(nodes, 2))
    root = filter(i -> isapprox(nodes[2, i], span / 2; atol=100eps()), section)
    outer = section[argmax(nodes[2, section])]
    @test length(root) == 2
    @test isapprox(nodes[3, outer], sum(nodes[3, root]) / 2; atol=100eps())
    @test isapprox(nodes[2, outer] - span / 2,
        (maximum(nodes[3, root]) - minimum(nodes[3, root])) / 2; atol=100eps())

    body = build_pitching_wing_body(chord, span;
        n_span=2, n_airfoil=21, n_endcap=5, semiinfinite_wake=true)
    @test body.nnodes == 114
    @test body.ncells == 224
    @test body.watertight
    @test pnl.has_dirichlet_bc(body)
    @test size(body.strength) == (body.ncells, 2)
    @test pnl.strength_names(body) == ("sigma", "mu")
    @test body.nsheddings == 2
    @test body.semiinfinite_wake

    pnl.calc_normals!(body)
    pnl.calc_controlpoints!(body)

    body_reflected_node, body_reflected_cell =
        _reflection_maps(body.nodes, body.cells)
    @test all(!iszero, body_reflected_node)
    @test all(!iszero, body_reflected_cell)
    for i in axes(body.cells, 2)
        j = body_reflected_cell[i]
        @test body.controlpoints[1, i] ≈ body.controlpoints[1, j] atol=1e-12 rtol=0
        @test body.controlpoints[2, i] ≈ body.controlpoints[2, j] atol=1e-12 rtol=0
        @test body.controlpoints[3, i] ≈ -body.controlpoints[3, j] atol=1e-12 rtol=0
        @test body.normals[1, i] ≈ body.normals[1, j] atol=1e-12 rtol=0
        @test body.normals[2, i] ≈ body.normals[2, j] atol=1e-12 rtol=0
        @test body.normals[3, i] ≈ -body.normals[3, j] atol=1e-12 rtol=0
    end

    # Each upper/lower shedding pair must use mirrored off-edge vertices at the
    # same span station.  Check the winding-normalized cells from which shedding
    # was constructed, rather than relying on raw mesh-local edge indices.
    point_line_distance(point, edge_a, edge_b) =
        norm(cross(point - edge_a, edge_b - edge_a)) / norm(edge_b - edge_a)
    for (pi, nia, nib, pj, nja, njb) in eachcol(body.shedding[1])
        third_i = only(setdiff(axes(body.cells, 1), (nia, nib)))
        third_j = only(setdiff(axes(body.cells, 1), (nja, njb)))
        vertex_i = body.nodes[:, body.cells[third_i, pi]]
        vertex_j = body.nodes[:, body.cells[third_j, pj]]
        upper_panel, lower_panel, upper_vertex, lower_vertex =
            vertex_i[3] > vertex_j[3] ?
            (pi, pj, vertex_i, vertex_j) : (pj, pi, vertex_j, vertex_i)

        @test upper_vertex[1:2] ≈ lower_vertex[1:2]
        @test upper_vertex[3] ≈ -lower_vertex[3]
        @test body.controlpoints[1:2, upper_panel] ≈
              body.controlpoints[1:2, lower_panel]

        edge_a = body.nodes[:, body.cells[nia, pi]]
        edge_b = body.nodes[:, body.cells[nib, pi]]
        @test point_line_distance(body.controlpoints[:, upper_panel], edge_a, edge_b) ≈
              point_line_distance(body.controlpoints[:, lower_panel], edge_a, edge_b)
    end

    right_cap = findall(
        i -> body.controlpoints[2, i] > span / 2, axes(body.controlpoints, 2))
    left_cap = findall(
        i -> body.controlpoints[2, i] < -span / 2, axes(body.controlpoints, 2))
    @test all(>(0), body.normals[2, right_cap])
    @test all(<(0), body.normals[2, left_cap])

    freestream = [cosd(5.0), 0.0, sind(5.0)]
    set_wake_Das!(body, freestream)
    freestream_direction = freestream / norm(freestream)
    @test all(Das -> all(col -> isapprox(col, freestream_direction), eachcol(Das)), body.Das)
    pnl.apply_freestream!(body, freestream)
    solver = pnl.Backslash(body)
    pnl.solve!(body, solver; backend=pnl.DirectBackend())
    @test all(isfinite, body.strength)
    @test any(!iszero, view(body.strength, :, 1))
    @test any(!iszero, view(body.strength, :, 2))
end

@testset "Pitching-wing zero-angle aerodynamic symmetry" begin
    cl0_direct = _pitching_wing_static_cl(0.0, pnl.DirectBackend())
    @test abs(cl0_direct) <= 1e-10

    fmm = pnl.FastMultipoleBackend(
        expansion_order=8, multipole_acceptance=0.4, leaf_size=40)
    cl0_fmm = _pitching_wing_static_cl(0.0, fmm)
    @test abs(cl0_fmm) <= 1e-6

    cl_positive = _pitching_wing_static_cl(3.0, pnl.DirectBackend())
    cl_negative = _pitching_wing_static_cl(-3.0, pnl.DirectBackend())
    @test cl_positive ≈ -cl_negative atol=1e-10 rtol=1e-10
end
