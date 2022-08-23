#=##############################################################################
# DESCRIPTION
    Unit tests of elements
=###############################################################################

using Test
import Printf: @printf
import FLOWPanel as pnl

verbose = true
v_lvl = 0



@testset verbose=verbose "Grid Tests" begin

    # --------------- PANEL DEFINITION TEST ------------------------------------
    # This test randomly compares the outputs of get_cell(...) and
    # get_cell_t!(...) to test consistency between cell definitions

    @test begin

        R = 1                            # (m) radius of sphere
        P_min = [0.15, 0, 0]             # Lower bounds of (theta, phi, dummy)
        P_max = [pi-P_min[1], 2*pi, 0]   # Upper bounds of (theta, phi, dummy)
        NDIVS = 1*[15, 30, 0]            # Number of divisions (cells) of (theta, phi, dummy)
        loop_dim = 2                     # Coordinate to loop (1==theta)

        # Generates parametric (theta, phi) grid
        grid = pnl.gt.Grid(P_min, P_max, NDIVS, loop_dim)

        # Transforms the grid into a spherical cartesian space
        my_transform(X) = pnl.gt.spherical3D(vcat(R, X[1:2]))
        pnl.gt.transform!(grid, my_transform)

        # Splits the quad cells into triangular cells
        dimsplit = 1
        triang_grid = pnl.gt.GridTriangleSurface(grid, dimsplit)

        # Pre-allocate memory for get_cell_t!(...)
        lin = LinearIndices(triang_grid._ndivsnodes)
        ndivscells = vcat(triang_grid._ndivscells...)

        tri_out = zeros(Int, 3)
        quadcoor = zeros(Int, 3)
        quad_out = zeros(Int, 4)

        # Randomized test
        function this_test()
            coor = [ndiv!=0 ? rand(1:ndiv) : 1 for ndiv in triang_grid._ndivscells]
            A = pnl.gt.get_cell(triang_grid, coor)
            B = pnl.gt.get_cell_t!(tri_out, quadcoor, quad_out, triang_grid, coor, lin, ndivscells)

            return prod(A.==B)
        end

        # Test result
        prod(this_test() for i in 1:10*triang_grid.ncells)

    end

end
