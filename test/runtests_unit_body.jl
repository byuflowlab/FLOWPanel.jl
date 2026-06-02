using Test
import FLOWPanel as pnl

function shared_edge_direction(cells::AbstractMatrix{Int}, ci::Int, a::Int, b::Int)
    n1, n2, n3 = cells[1, ci], cells[2, ci], cells[3, ci]
    if (n1 == a && n2 == b) || (n2 == a && n3 == b) || (n3 == a && n1 == b)
        return 1
    elseif (n1 == b && n2 == a) || (n2 == b && n3 == a) || (n3 == b && n1 == a)
        return -1
    end
    error("Cell $ci does not contain edge ($a, $b).")
end

const NODES_TET = Float64[
    1  -1  -1   1;
    1  -1   1  -1;
    1   1  -1  -1;
]
const CELLS_TET = Int[
    1  1  1  2;
    3  4  2  4;
    2  3  4  3;
]

const NODES_NONMANIFOLD = Float64[
    0 1 0 0 0;
    0 0 1 -1 0;
    0 0 0 0 1;
]
const CELLS_NONMANIFOLD = Int[
    1 2 1;
    2 1 2;
    3 4 5;
]

@testset verbose=true "Abstract Body Geometry And NonLiftingBody" begin
    @testset "ensure_consistent_winding" begin
        flipped_cells = copy(CELLS_TET)
        flipped_cells[:, 2] = flipped_cells[[1, 3, 2], 2]
        flipped_cells[:, 4] = flipped_cells[[1, 3, 2], 4]

        consistent_cells = pnl.ensure_consistent_winding(NODES_TET, flipped_cells; watertight=true)
        edge_to_cells = pnl._calc_edge_to_cells(consistent_cells)
        for ((a, b), refs) in edge_to_cells
            if length(refs) == 2
                c1 = refs[1][1]
                c2 = refs[2][1]
                @test shared_edge_direction(consistent_cells, c1, a, b) == -shared_edge_direction(consistent_cells, c2, a, b)
            end
        end

        normals = pnl.calc_normals(NODES_TET, consistent_cells)
        body_centroid = vec(sum(NODES_TET; dims=2) ./ size(NODES_TET, 2))
        for ci in axes(consistent_cells, 2)
            centroid = vec(sum(view(NODES_TET, :, consistent_cells[:, ci]); dims=2) ./ 3)
            @test dot(normals[:, ci], centroid - body_centroid) > 0
        end

        flipped_normals_cells = pnl.ensure_consistent_winding(NODES_TET, flipped_cells; watertight=true, flip_normals=true)
        flipped_normals = pnl.calc_normals(NODES_TET, flipped_normals_cells)
        @test isapprox(flipped_normals, -normals; atol=1e-12)
    end

    @testset "iswatertight" begin
        @test pnl.iswatertight(NODES_TET, CELLS_TET) == (true, Int[])
        @test pnl.iswatertight(NODES_OCT, CELLS_OCT) == (true, Int[])

        @test pnl.iswatertight(NODES_2TRI, CELLS_2TRI) == (false, Int[])
        @test pnl.iswatertight(NODES_2TRI, CELLS_2TRI; return_open_cells=true) == (false, [1, 2])

        @test pnl.iswatertight(NODES_NONMANIFOLD, CELLS_NONMANIFOLD) == (false, Int[])
        @test pnl.iswatertight(NODES_NONMANIFOLD, CELLS_NONMANIFOLD; return_open_cells=true) == (false, [1, 2, 3])

        open_nodes, open_cells = make_basic_triangle_surface()
        @test pnl.iswatertight(open_nodes, open_cells; return_open_cells=true) == (false, [1, 2])

        closed_nodes, closed_cells = make_octa_triangle_surface()
        closed_grid_raw = pnl.iswatertight(closed_nodes, closed_cells; return_open_cells=true)
        @test closed_grid_raw == (true, Int[])

        open_body = make_nonlifting(pnl.ConstantSource)
        open_body_raw = pnl.iswatertight(open_body.nodes, open_body.cells; return_open_cells=true)
        @test pnl.iswatertight(open_body; return_open_cells=true) == open_body_raw

        closed_body = pnl.NonLiftingBody{pnl.ConstantSource}(copy(NODES_OCT), copy(CELLS_OCT); watertight=true)
        closed_body_raw = pnl.iswatertight(closed_body.nodes, closed_body.cells; return_open_cells=true)
        @test closed_body_raw == (true, Int[])
        @test pnl.iswatertight(closed_body; return_open_cells=true) == closed_body_raw
    end

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

        pnl.calc_controlpoints!(NODES_1TRI, CELLS_1TRI, cps, normals)
        centroid = vec(sum(NODES_1TRI; dims=2) ./ 3)
        @test isapprox(cps[:, 1], centroid; atol=1e-12)

        body_default = pnl.NonLiftingBody{pnl.ConstantSource}(copy(NODES_1TRI), copy(CELLS_1TRI))
        pnl.calc_normals!(body_default)
        pnl.calc_controlpoints!(body_default)
        @test body_default.characteristiclength === pnl.characteristiclength_sqrtarea
        @test isapprox(body_default.controlpoints[:, 1], centroid; atol=1e-12)
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

    @testset "rotate! reset! get_cell no grid2cells" begin
        body = make_nonlifting(pnl.ConstantSource)
        original_nodes = copy(body.nodes)
        pnl.rotate!(body, 0, 0, 90)
        expected = [0 1 0; -1 0 0; 0 0 1] * original_nodes
        @test isapprox(body.nodes, expected; atol=1e-10)
        @test isapprox(body.Oaxis, [0 1 0; -1 0 0; 0 0 1]; atol=1e-10)

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
        body4.velocity_gradient .= 3.0
        body4.induced_vorticity .= 4.0
        pnl.reset!(body4)
        @test all(body4.velocity .== 0)
        @test all(body4.potential .== 0)
        @test all(body4.velocity_gradient .== 0)
        @test all(body4.induced_vorticity .== 0)

        @test pnl.get_cell(body4, 1) == (CELLS_2TRI[1, 1], CELLS_2TRI[2, 1], CELLS_2TRI[3, 1])

        @test !isdefined(pnl, :grid2cells)
    end

    @testset "NonLiftingBody construction" begin
        body = make_nonlifting(pnl.ConstantSource)
        @test body.nnodes == 4
        @test body.ncells == 2
        @test size(body.strength) == (2, 1)
        @test body.watertight == false
        @test length(body.vtk_cells) == 2

        body_doublet = make_nonlifting(pnl.ConstantDoublet)
        @test size(body_doublet.strength) == (2, 1)

        body_union = make_nonlifting(Union{pnl.ConstantSource, pnl.ConstantDoublet})
        @test size(body_union.strength) == (2, 2)

        @test pnl.has_dirichlet_bc(make_nonlifting(pnl.ConstantSource; DBC=true)) == true
        @test pnl.has_dirichlet_bc(make_nonlifting(pnl.ConstantSource)) == false

        tri_nodes, tri_cells = make_basic_triangle_surface()
        grid_constructor_body = pnl.NonLiftingBody{pnl.ConstantSource}(tri_nodes, tri_cells)
        @test size(grid_constructor_body.neighbor) == (3, grid_constructor_body.ncells)
        @test grid_constructor_body.nnodes == size(tri_nodes, 2)
        @test grid_constructor_body.ncells == size(tri_cells, 2)
        @test grid_constructor_body.watertight == false
        @test !(:P in fieldnames(typeof(grid_constructor_body)))
        @test !(:F in fieldnames(typeof(grid_constructor_body)))

        closed_nodes, closed_cells = make_octa_triangle_surface()
        closed_grid_body = pnl.NonLiftingBody{pnl.ConstantSource}(closed_nodes, closed_cells; watertight=true)
        @test closed_grid_body.watertight == true

        flipped_cells = copy(CELLS_TET)
        flipped_cells[:, 3] = flipped_cells[[1, 3, 2], 3]
        flipped_body = pnl.NonLiftingBody{pnl.ConstantSource}(copy(NODES_TET), flipped_cells; watertight=true)
        pnl.calc_normals!(flipped_body)
        body_centroid = vec(sum(flipped_body.nodes; dims=2) ./ flipped_body.nnodes)
        for ci in 1:flipped_body.ncells
            centroid = vec(sum(view(flipped_body.nodes, :, flipped_body.cells[:, ci]); dims=2) ./ 3)
            @test dot(flipped_body.normals[:, ci], centroid - body_centroid) > 0
        end

        raw_body = pnl.NonLiftingBody{pnl.ConstantSource}(copy(NODES_TET), flipped_cells; watertight=true, ensure_winding=false)
        @test raw_body.cells == flipped_cells
    end

    @testset "calcfield_Cp" begin
        rho = 1.0
        uinf = 10.0
        qinf = 0.5 * rho * uinf^2

        @test pnl.calcfield_Cp([qinf], uinf, rho) == [1.0]
        @test pnl.calcfield_Cp([0.0], uinf, rho) == [0.0]
        @test pnl.calcfield_Cp([-3.0 * qinf], uinf, rho) == [-3.0]
    end
end
