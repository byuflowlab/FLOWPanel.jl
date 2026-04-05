using Test
import FLOWPanel as pnl
import GeometricTools as gt

@testset verbose=true "Abstract Body Geometry And NonLiftingBody" begin
    @testset "calc_normals!" begin
        normals = zeros(3, 1)
        pnl.calc_normals!(NODES_1TRI, CELLS_1TRI, normals)
        @test isapprox(norm(normals[:, 1]), 1.0; atol=1e-12)
        @test isapprox(abs(normals[3, 1]), 1.0; atol=1e-12)

        flipped = reshape(Int[1, 3, 2], 3, 1)
        flipped_normals = zeros(3, 1)
        pnl.calc_normals!(NODES_1TRI, flipped, flipped_normals)
        @test isapprox(flipped_normals[:, 1], -normals[:, 1]; atol=1e-12)

        body = make_nonlifting(pnl.ConstantSource)
        pnl.calc_normals!(body)
        @test size(body.normals) == (3, 2)
        @test all(isapprox.(vec(sqrt.(sum(body.normals .^ 2; dims=1))), 1.0; atol=1e-12))
    end

    @testset "calc_areas!" begin
        areas = zeros(1)
        pnl.calc_areas!(NODES_1TRI, CELLS_1TRI, areas)
        @test areas[1] == 0.5

        scaled_nodes = Float64[
            0 3 0;
            0 0 4;
            0 0 0;
        ]
        scaled_areas = zeros(1)
        pnl.calc_areas!(scaled_nodes, CELLS_1TRI, scaled_areas)
        @test scaled_areas[1] == 6.0

        equi_nodes = Float64[
            0 2 1;
            0 0 sqrt(3);
            0 0 0;
        ]
        equi_areas = pnl.calc_areas(equi_nodes, CELLS_1TRI)
        @test isapprox(equi_areas[1], sqrt(3); atol=1e-12)

        body = make_nonlifting(pnl.ConstantSource)
        areas_body = pnl.calc_areas(body)
        @test length(areas_body) == body.ncells
    end

    @testset "calc_controlpoints!" begin
        normals = zeros(3, 1)
        pnl.calc_normals!(NODES_1TRI, CELLS_1TRI, normals)
        cps = zeros(3, 1)

        pnl.calc_controlpoints!(NODES_1TRI, CELLS_1TRI, cps, normals; off=0.0)
        centroid = vec(sum(NODES_1TRI; dims=2) ./ 3)
        @test isapprox(cps[:, 1], centroid; atol=1e-12)

        pnl.calc_controlpoints!(NODES_1TRI, CELLS_1TRI, cps, normals; off=0.25, characteristiclength=pnl.characteristiclength_unitary)
        @test isapprox(cps[:, 1] - centroid, 0.25 .* normals[:, 1]; atol=1e-12)

        cps_bbox = zeros(3, 1)
        pnl.calc_controlpoints!(NODES_1TRI, CELLS_1TRI, cps_bbox, normals; off=0.5, characteristiclength=pnl.characteristiclength_bbox)
        @test isapprox(norm(cps_bbox[:, 1] - centroid), 0.5 * sqrt(2); atol=1e-12)
    end

    @testset "characteristic lengths" begin
        panel = view(CELLS_1TRI, :, 1)
        @test pnl.characteristiclength_unitary(NODES_1TRI, panel) == 1
        @test isapprox(pnl.characteristiclength_bbox(NODES_1TRI, panel), sqrt(2); atol=1e-12)
        @test isapprox(pnl.characteristiclength_maxdist(NODES_1TRI, panel), 1.0; atol=1e-12)
        @test isapprox(pnl.characteristiclength_sqrtarea(NODES_1TRI, panel), sqrt(0.5); atol=1e-12)
    end

    @testset "calc_tangents! and calc_obliques!" begin
        body = make_nonlifting(pnl.ConstantSource)
        pnl.calc_normals!(body)
        tangents = pnl.calc_tangents(body)
        obliques = pnl.calc_obliques(body)
        tangents2 = pnl.calc_tangents(body)
        obliques2 = pnl.calc_obliques(body)

        for i in 1:body.ncells
            @test abs(dot(tangents[:, i], body.normals[:, i])) < 1e-12
            @test abs(dot(obliques[:, i], body.normals[:, i])) < 1e-12
            @test abs(dot(tangents[:, i], obliques[:, i])) < 1e-12
            @test isapprox(norm(tangents[:, i]), 1.0; atol=1e-12)
            @test isapprox(norm(obliques[:, i]), 1.0; atol=1e-12)
        end

        @test tangents == tangents2
        @test obliques == obliques2
    end

    @testset "rotate! reset! get_cell grid2cells" begin
        body = make_nonlifting(pnl.ConstantSource)
        original_nodes = copy(body.nodes)
        pnl.rotate!(body, 0, 0, 90)
        expected = gt.rotation_matrix2(0, 0, -90) * original_nodes
        @test isapprox(body.nodes, expected; atol=1e-10)
        @test isapprox(body.Oaxis, gt.rotation_matrix2(0, 0, -90); atol=1e-10)

        body2 = make_nonlifting(pnl.ConstantSource)
        pnl.rotate!(body2, 0, 0, 0; translation=[1.0, 0.0, 0.0])
        @test isapprox(body2.nodes, NODES_2TRI .+ [1.0, 0.0, 0.0]; atol=1e-12)

        body3 = make_nonlifting(pnl.ConstantSource)
        pnl.rotate!(body3, 0, 0, 90)
        pnl.rotate!(body3, 0, 0, 90)
        body4 = make_nonlifting(pnl.ConstantSource)
        pnl.rotate!(body4, 0, 0, 180)
        @test isapprox(body3.nodes, body4.nodes; atol=1e-10)

        body4.velocity .= 1.0
        body4.potential .= 2.0
        body4.Cp .= 3.0
        body4.F .= 4.0
        pnl.reset!(body4)
        @test all(body4.velocity .== 0)
        @test all(body4.potential .== 0)
        @test all(body4.Cp .== 0)
        @test all(body4.F .== 0)

        @test pnl.get_cell(body4, 1) == (CELLS_2TRI[1, 1], CELLS_2TRI[2, 1], CELLS_2TRI[3, 1])

        tri_grid = make_basic_triangle_surface()
        grid_cells = pnl.grid2cells(tri_grid)
        @test size(grid_cells, 1) == 3
        @test size(grid_cells, 2) == tri_grid.ncells
        @test minimum(grid_cells) >= 1
        @test maximum(grid_cells) <= size(tri_grid._nodes, 2)
    end

    @testset "NonLiftingBody construction" begin
        body = make_nonlifting(pnl.ConstantSource)
        @test body.nnodes == 4
        @test body.ncells == 2
        @test size(body.strength) == (2, 1)
        @test body.CPoffset == 1e-14
        @test body.watertight == false
        @test length(body.vtk_cells) == 2

        body_doublet = make_nonlifting(pnl.ConstantDoublet)
        @test size(body_doublet.strength) == (2, 1)

        body_union = make_nonlifting(Union{pnl.ConstantSource, pnl.ConstantDoublet})
        @test size(body_union.strength) == (2, 2)

        @test pnl.has_dirichlet_bc(make_nonlifting(pnl.ConstantSource; DBC=true)) == true
        @test pnl.has_dirichlet_bc(make_nonlifting(pnl.ConstantSource)) == false

        tri_grid = make_basic_triangle_surface()
        grid_body = pnl.NonLiftingBody{pnl.ConstantSource}(tri_grid)
        @test size(grid_body.neighbor) == (3, grid_body.ncells)
        @test grid_body.nnodes == size(tri_grid._nodes, 2)
        @test grid_body.ncells == tri_grid.ncells
    end
end
