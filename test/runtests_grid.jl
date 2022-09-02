#=##############################################################################
# DESCRIPTION
    Unit tests of elements
=###############################################################################

using Test
import Printf: @printf
import FLOWPanel as pnl

try
    verbose
catch
    global verbose = true
end
v_lvl = 0



@testset verbose=verbose "Grid Tests" begin

    R = 1.5                          # (m) radius of sphere
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

    # Generate body object
    body = pnl.NonLiftingBody{pnl.ConstantSource}(triang_grid)

    # Unitary vectors
    tangents = pnl._calc_tangents(body)
    obliques = pnl._calc_obliques(body)
    normals = pnl._calc_normals(body)

    # Pre-allocate memory for get_cell_t!(...)
    lin = LinearIndices(triang_grid._ndivsnodes)
    ndivscells = vcat(triang_grid._ndivscells...)
    cin = CartesianIndices(Tuple(collect( 1:(d != 0 ? d : 1) for d in triang_grid._ndivscells)))

    tri_out = zeros(Int, 3)
    tricoor = zeros(Int, 3)
    quadcoor = zeros(Int, 3)
    quad_out = zeros(Int, 4)

    # --------------- PANEL DEFINITION TEST ------------------------------------
    # This test randomly compares the outputs of get_cell(...) and
    # get_cell_t!(...) to test consistency between cell definitions
    @test begin

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



    # --------------- ORTHOGONAL UNITARY VECTORS (TEST GeometricTools) ---------
    @test begin

        nodes = triang_grid.orggrid.nodes
        crss = zeros(3)

        result_t = true
        result_o = true
        result_n = true
        result_ortt = true
        result_orto = true
        result_ortn = true
        result_crss = true

        for pi in 1:triang_grid.ncells

            panel = pnl.gt.get_cell_t!(tri_out, tricoor, quadcoor, quad_out,
                                                 triang_grid, pi, lin, ndivscells, cin)

            # Tangent, oblique, and normal vectors
            t1, t2, t3 = pnl.gt._calc_t1(nodes, panel), pnl.gt._calc_t2(nodes, panel), pnl.gt._calc_t3(nodes, panel)
            o1, o2, o3 = pnl.gt._calc_o1(nodes, panel), pnl.gt._calc_o2(nodes, panel), pnl.gt._calc_o3(nodes, panel)
            n1, n2, n3 = pnl.gt._calc_n1(nodes, panel), pnl.gt._calc_n2(nodes, panel), pnl.gt._calc_n3(nodes, panel)

            # Test that vectors are unitary
            result_t *= abs(sqrt(t1^2 + t2^2 + t3^2) - 1) < 1e1*eps()
            result_o *= abs(sqrt(o1^2 + o2^2 + o3^2) - 1) < 1e1*eps()
            result_n *= abs(sqrt(n1^2 + n2^2 + n3^2) - 1) < 1e1*eps()

            # Test that vectors are orthogonal
            result_ortt *= pnl.dot((t1,t2,t3), (o1,o2,o3)) < 1e1*eps()
            result_orto *= pnl.dot((o1,o2,o3), (n1,n2,n3)) < 1e1*eps()
            result_ortn *= pnl.dot((n1,n2,n3), (t1,t2,t3)) < 1e1*eps()

            # Test that vectors make a basis
            crss .= 0; pnl.cross!(crss, (t1,t2,t3), (o1,o2,o3))
            result_crss *= pnl.norm(crss .- [n1,n2,n3]) < 1e1*eps()
        end

        # Test result
        result_t && result_o && result_n && result_ortt && result_orto && result_ortn && result_crss
    end


    # --------------- ORTHOGONAL UNITARY VECTORS (TEST Body)--------------------
    # Test that vectors are unitary
    @test prod(abs.( pnl.norm.(eachcol(tangents)) .- 1 ) .< 1e1*eps())
    @test prod(abs.( pnl.norm.(eachcol(obliques)) .- 1 ) .< 1e1*eps())
    @test prod(abs.( pnl.norm.(eachcol(normals))  .- 1 ) .< 1e1*eps())

    # Test that vectors are orthogonal
    @test prod(pnl.dot.(eachcol(tangents), eachcol(obliques)) .< 1e1*eps())
    @test prod(pnl.dot.(eachcol(obliques), eachcol(normals)) .< 1e1*eps())
    @test prod(pnl.dot.(eachcol(normals), eachcol(tangents)) .< 1e1*eps())

    # Test that vectors make a orthonormal basis
    @test prod(pnl.norm.( pnl.cross.(eachcol(tangents), eachcol(obliques)) .- eachcol(normals) ) .< 1e1*eps())

end
