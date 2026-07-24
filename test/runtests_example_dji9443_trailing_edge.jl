using Test
import FLOWPanel as pnl
import GeoIO

include(joinpath(pnl.examples_path, "dji9443_trailing_edge.jl"))

@testset verbose=true "DJI 9443 trailing-edge detection" begin
    cases = [
        ("dji9443_20260722_40_41_capped.msh", true,
         ([3333, 3293, 1715], [1619, 1579, 1]), 39),
        ("dji9443_20260722_40_41_uncapped.msh", false,
         ([3162, 3122, 1601], [1562, 1522, 1]), 39),
        ("dji9443_20260722_57_57_capped.msh", true,
         ([6573, 6517, 3355], [3219, 3163, 1]), 56),
        ("dji9443_20260722_57_57_uncapped.msh", false,
         ([6330, 6274, 3193], [3138, 3082, 1]), 56),
    ]

    @testset "tracked meshes and shedding construction" begin
        for (filename, watertight, expected, nsections) in cases
            mesh_file = joinpath(pnl.examples_path, "data", filename)
            @test find_dji9443_trailing_edge_indices(mesh_file; watertight) == expected

            mesh = GeoIO.load(mesh_file).geometry
            nodes, cells = pnl.meshes2nodes_cells(mesh)
            kernel = watertight ?
                     Union{pnl.ConstantSource, pnl.VortexRing} :
                     pnl.VortexRing
            base = pnl.RigidWakeBody{kernel}(
                nodes, cells, pnl.noshedding;
                watertight, DBC=watertight, semiinfinite_wake=false)
            shedding = map(expected) do indices
                pnl.calc_shedding_from_seed(
                    base.nodes, base.cells, indices[1], indices[2];
                    end_node=indices[3], normal_jump_tol=0.2,
                    max_turn_angle=pi / 3)
            end
            @test all(size(chain) == (6, nsections) for chain in shedding)

            body = pnl.RigidWakeBody{kernel}(
                copy(base.nodes), copy(base.cells), collect(shedding);
                watertight, DBC=watertight, semiinfinite_wake=false,
                ensure_winding=true)
            @test body.nsheddings == 2nsections
        end
    end

    @testset "uniform scaling" begin
        filename, watertight, expected, _ = cases[1]
        mesh = GeoIO.load(joinpath(pnl.examples_path, "data", filename)).geometry
        nodes, cells = pnl.meshes2nodes_cells(mesh)
        for scale in (1.0e-3, 7.5, 1.0e3)
            @test _find_dji9443_trailing_edge_indices(
                scale .* nodes, cells; watertight) == expected
        end
    end

    @testset "topology mismatch" begin
        for (filename, watertight, _, _) in cases
            mesh_file = joinpath(pnl.examples_path, "data", filename)
            error = try
                find_dji9443_trailing_edge_indices(
                    mesh_file; watertight=!watertight)
                nothing
            catch err
                err
            end
            @test error isa ErrorException
            @test occursin("Topology mismatch", sprint(showerror, error))
        end
    end
end
