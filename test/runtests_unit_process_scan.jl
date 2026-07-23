using Test

include(joinpath(@__DIR__, "..", "examples", "rotor_hover_scan_new", "process_scan.jl"))

@testset "section review output layout" begin
    mktempdir() do output_dir
        cst_path = cst_output_path(output_dir)
        paths = section_review_paths(3; output_dir)
        @test dirname(paths.directory) == dirname(cst_path)
        @test basename(paths.directory) == "section_03_airfoil_review"
        @test all(dirname(path) == paths.directory for path in
                  (paths.raw, paths.airfoil, paths.leading_edge, paths.raw_selection,
                   paths.normalization, paths.blunt, paths.taper, paths.final))
        @test basename(paths.leading_edge) == "section_03_leading_edge_row_guide.png"
        @test basename(paths.raw_selection) == "section_03_checkpoint_1_raw_selection.png"
        @test basename(paths.normalization) ==
              "section_03_checkpoint_1_normalization_and_canonicalization.png"
        @test basename(paths.blunt) == "section_03_checkpoint_2_blunt_cst.png"
        @test basename(paths.taper) == "section_03_checkpoint_3_taper.png"
        @test basename(paths.final) == "section_03_checkpoint_4_final_cst.png"
    end
end

@testset "scan-processing input parsers" begin
    @test parse_coordinate_pair("1.25, -0.03") == (1.25, -0.03)
    @test parse_coordinate_pair("1.25 -0.03") == (1.25, -0.03)
    @test parse_coordinate_pair(""; default=(1.0, 0.1)) == (1.0, 0.1)
    @test_throws ErrorException parse_coordinate_pair("1.0")
    @test_throws ErrorException parse_coordinate_pair("a,b")

    @test parse_exclusion_rows("", 40) == Int[]
    @test parse_exclusion_rows("12", 40) == [12]
    @test parse_exclusion_rows("24,15-18,12,16", 40) == [12, 15, 16, 17, 18, 24]
    @test_throws ErrorException parse_exclusion_rows("1,,2", 40)
    @test_throws ErrorException parse_exclusion_rows("8-3", 40)
    @test_throws ErrorException parse_exclusion_rows("-1", 40)
    @test_throws ErrorException parse_exclusion_rows("0", 40)
    @test_throws ErrorException parse_exclusion_rows("41", 40)
    @test_throws ErrorException parse_exclusion_rows("5", 40; le_idx=5)
    @test_throws ErrorException parse_exclusion_rows("1-9", 40)
end

@testset "initial section skip prompt" begin
    for answer in ("y\n", "YES\n")
        @test prompt_skip_section(3; input=IOBuffer(answer), output=IOBuffer())
    end
    for answer in ("n\n", "No\n", "\n")
        @test !prompt_skip_section(3; input=IOBuffer(answer), output=IOBuffer())
    end
    output = IOBuffer()
    @test prompt_skip_section(3; input=IOBuffer("maybe\ny\n"), output)
    prompt_text = String(take!(output))
    @test occursin("Skip section 3? [y/N]:", prompt_text)
    @test occursin("Please answer y or n.", prompt_text)

    airfoils = Any[]
    all_params = Any[]
    @test !record_processed_section!(airfoils, all_params, nothing)
    @test isempty(airfoils)
    @test isempty(all_params)
    first_result = (:first_airfoil, [1.0, 2.0])
    second_result = (:second_airfoil, [3.0, 4.0])
    @test record_processed_section!(airfoils, all_params, first_result)
    @test !record_processed_section!(airfoils, all_params, nothing)
    @test record_processed_section!(airfoils, all_params, second_result)
    @test airfoils == [:first_airfoil, :second_airfoil]
    @test all_params == [[1.0, 2.0], [3.0, 4.0]]
end

function suggestion_slice()
    x_surface = [0.85, 0.9, 0.95, 0.97]
    upper = 0.08 .- 0.02 .* x_surface
    lower = -0.06 .+ 0.01 .* x_surface
    nodes = hcat([0.0, 0.0], [1.0, 0.0],
                 [[x, y] for (x, y) in zip(x_surface, upper)]...,
                 [[x, y] for (x, y) in zip(x_surface, lower)]...,
                 [0.98, 0.03], [1.02, 0.0], [0.98, -0.03])
    return MeshSlice(0.5, nodes)
end

@testset "square trailing-edge suggestion" begin
    ms = suggestion_slice()
    suggestion = suggest_te_anchors(ms, 1, [11, 12, 13])
    @test suggestion.xstar ≈ 1.02
    @test suggestion.upper[1] ≈ 1.02
    @test suggestion.upper[2] ≈ 0.08 - 0.02 * 1.02 atol=1e-12
    @test suggestion.lower[1] ≈ 1.02
    @test suggestion.lower[2] ≈ -0.06 + 0.01 * 1.02 atol=1e-12

    fallback = suggest_te_anchors(ms, 1, Int[])
    @test fallback.xstar ≈ 1.0

    sparse = MeshSlice(0.5, hcat([0.0, 0.0], [1.0, 0.0], [0.9, 0.04],
                                 [0.9, -0.04], [0.8, 0.03], [0.8, -0.03]))
    no_suggestion = suggest_te_anchors(sparse, 1, Int[])
    @test no_suggestion.upper === nothing
    @test no_suggestion.lower === nothing
end

@testset "normalization and canonical anchors" begin
    nodes = hcat([0.0, 0.0], [0.5, 0.12], [0.8, 0.15], [0.8, 0.05])
    upper = (1.05, 0.30)
    lower = (0.95, 0.10)
    ms = AnchoredMeshSlice(0.5, nodes, upper, 1, lower, Int[], 0.8)
    normalized = normalized_section(ms)
    midpoint = reshape(0.5 .* (collect(upper) .+ collect(lower)), 2, 1)
    mx, my = transform_points(midpoint, normalized.le_point,
                              normalized.twist, normalized.chord)
    @test mx[1] ≈ 1.0 atol=1e-14
    @test my[1] ≈ 0.0 atol=1e-14
    @test normalized.canonical_x == [1.0, 1.0]
    @test normalized.canonical_y == [normalized.h_te, -normalized.h_te]
    @test 0.5 .* (normalized.canonical_y .+ reverse(normalized.canonical_y)) == [0.0, 0.0]
    projected_half_distance = abs(normalized.anchor_y[1] - normalized.anchor_y[2]) / 2
    full_half_distance = norm(collect(upper) .- collect(lower)) / (2 * normalized.chord)
    @test normalized.h_te ≈ projected_half_distance
    @test normalized.h_te < full_half_distance
    @test_throws ErrorException normalized_section(
        AnchoredMeshSlice(0.5, nodes, lower, 1, upper, Int[], 0.8))
    @test_throws ErrorException normalized_section(
        AnchoredMeshSlice(0.5, nodes, upper, 1, upper, Int[], 0.8))
end

function synthetic_anchored_section(n=31; noise=0.0)
    x = collect(range(0.0, 0.97; length=n))
    camber = @. 0.02 * x * (1 - x)
    thickness = @. 0.12 * sqrt(x) * (1 - x) + 0.02 * x
    perturbation = @. noise * sin(37x)
    upper = camber .+ thickness .+ perturbation
    lower = camber .- thickness .- perturbation
    nodes = hcat(vcat(x', upper'), vcat(x', lower'))
    return AnchoredMeshSlice(0.6, nodes, (1.0, 0.02), 1, (1.0, -0.02), Int[], 0.8)
end

@testset "blunt CST reconstruction" begin
    ms = synthetic_anchored_section(35; noise=2e-4)
    normalized = normalized_section(ms)
    blunt = fit_blunt_cst(ms, normalized)
    @test all(isfinite, blunt.params)
    @test blunt.params[end] ≈ 2 * normalized.h_te
    @test blunt.yu[end] == normalized.h_te
    @test blunt.yl[end] == -normalized.h_te
    @test blunt.camber[end] == 0.0
    @test_throws ErrorException fit_blunt_cst(synthetic_anchored_section(8))
end

@testset "cosine grid and camber-preserving taper" begin
    blunt = fit_blunt_cst(synthetic_anchored_section(35))
    x = blunt.x_common
    @test x[1] == 0.0
    @test x[end] == 1.0
    @test all(diff(x) .> 0.0)
    for blend_start in (0.8, 0.7)
        sharp = taper_blunt_geometry(blunt, blend_start)
        forward = x .<= blend_start
        @test sharp.yu[forward] == blunt.yu[forward]
        @test sharp.yl[forward] == blunt.yl[forward]
        @test sharp.camber ≈ blunt.camber atol=10eps(Float64) rtol=0.0
        @test minimum(sharp.half_thickness) >= 0.0
        @test sharp.half_thickness[end] == 0.0
    end
    @test quadratic_thickness_weight(0.0) == 1.0
    @test quadratic_thickness_weight(1.0) == 0.0
    @test quadratic_thickness_weight_derivative(0.0) == 0.0
    @test quadratic_thickness_weight_derivative(1.0) == -2.0
end

@testset "final sharp CST and diagnostics" begin
    blunt = fit_blunt_cst(synthetic_anchored_section(35; noise=1e-4))
    sharp = taper_blunt_geometry(blunt, 0.8)
    final = fit_final_cst(sharp)
    @test all(isfinite, final.params)
    @test final.params[end] == 0.0
    @test final.yu[end] == 0.0
    @test final.yl[end] == 0.0
    diagnostics = reconstruction_diagnostics(blunt, sharp, final)
    @test all(isfinite, diagnostics.camber_error)
    @test isfinite(diagnostics.final_upper_stats.rms)
    @test isfinite(diagnostics.final_lower_stats.maxabs)
    af = Airfoil(synthetic_anchored_section(35; noise=1e-4))
    @test af.y_u[end] == 0.0
    @test af.y_l[end] == 0.0
    mktempdir() do directory
        path = joinpath(directory, "section.csv")
        write_cst_table(path, [af], [final.params])
        lines = readlines(path)
        expected_header = vcat(["r", "twist_rad", "chord", "axial", "tangential"],
                               ["wu$i" for i in 1:8], ["wl$i" for i in 1:8],
                               ["leading_edge_weight", "dz"])
        @test split(lines[1], ',') == expected_header
        @test parse(Float64, split(lines[2], ',')[end]) == 0.0
    end
end

@testset "headless diagnostic figures" begin
    ms = synthetic_anchored_section(35)
    normalized = normalized_section(ms)
    blunt = fit_blunt_cst(ms, normalized)
    sharp = taper_blunt_geometry(blunt, 0.8)
    final = fit_final_cst(sharp)
    diagnostics = reconstruction_diagnostics(blunt, sharp, final)
    raw_selection, normalization = plot_selection_checkpoint(ms, normalized, nothing, 1)
    figures = (
        plot_le_row_selection(MeshSlice(ms.r, ms.nodes), 1),
        raw_selection,
        normalization,
        plot_blunt_diagnostics(ms, normalized, blunt, diagnostics, 1),
        plot_taper_diagnostics(blunt, sharp, 0.8, 1),
        plot_final_diagnostics(sharp, final, diagnostics, 1),
    )
    show_interactive_figures(figures...; pause_seconds=0.001)
    mktempdir() do directory
        for (i, figure) in enumerate(figures)
            path = joinpath(directory, "diagnostic_$i.png")
            figure.savefig(path)
            @test isfile(path)
            plt.close(figure)
        end
    end
end
