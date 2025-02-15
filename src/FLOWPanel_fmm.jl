# using Revise

################################################################################
# FMM WRAPPER
################################################################################
module FMM

    using StaticArrays
    import LinearAlgebra
    using LinearAlgebra: cross, norm, dot, mul!, lu!, LU
    using WriteVTK
    import Krylov
    # import LinearOperators

    # https://github.com/byuflowlab/FastMultipole
    using FastMultipole

    const ONE_OVER_4PI = 1/4/pi

    for header_name in ["types", "panel", "kernel", "geometry", "fmm",
                        "solve", "freestream", "vtk"]

      include(joinpath("fmm", header_name*".jl"))

    end

end # END OF MODULE


################################################################################
# FMM-BASED SOLVER
################################################################################
function solvefmm(body::AbstractBody{VortexRing},
                    Uinfs::AbstractMatrix{T1},
                    Das::AbstractMatrix{T2},
                    Dbs::AbstractMatrix{T3},
                    args...;
                    optargs...
                    ) where {T1, T2, T3}

    T = promote_type(T1, T2, T3)

    # Convert body to FMM body
    body_fmm = body2fmmbody(body)

    # Call FMM-based solver
    solver = _solvefmm!(body, body_fmm, Uinfs, Das, Dbs, args...; optargs...)

    # Fetch strengths
    Gamma = [P.strength[1] for P in body_fmm.panels]

    # Store strengths as field Gamma
    add_field(body, "Gamma", "scalar", Gamma, "cell")

    return solver
end

function _solvefmm!(body::AbstractBody, body_fmm::FMM.UnstructuredGrid,
                    Uinfs::AbstractMatrix,
                    Das::AbstractMatrix,
                    Dbs::AbstractMatrix,
                    args...; optargs...)

    # Error cases
    if size(Uinfs) != (3, body.ncells)
        error("Invalid Uinfs;"*
              " expected size (3, $(body.ncells)), got $(size(Uinfs))")
    elseif size(Das) != (3, body.nsheddings)
        error("Invalid Das;"*
              " expected size (3, $(body.nsheddings)), got $(size(Das))")
    elseif size(Dbs) != (3, body.nsheddings)
        error("Invalid Dbs;"*
              " expected size (3, $(body.nsheddings)), got $(size(Dbs))")
    end

    # Solve linear system of equations
    solver = _solvefmm!(body_fmm, Uinfs, Das, Dbs, args...; optargs...)

    # Set strengths calculated through FMM-based solver back in the body
    for (Pi, P) in enumerate(body_fmm.panels)
        set_strength(body, Pi, P.strength)
    end

    _solvedflag(body, true)
    add_field(body, "Uinf", "vector", collect(eachcol(Uinfs)), "cell")
    add_field(body, "Da", "vector", collect(eachcol(Das)), "system")
    add_field(body, "Db", "vector", collect(eachcol(Dbs)), "system")

    return solver

end

function _solvefmm!(body_fmm::FMM.UnstructuredGrid,
                    Uinfs::AbstractMatrix,
                    Das::AbstractMatrix,
                    Dbs::AbstractMatrix;
                    solver=FMM.MatrixFreeSolver,
                    scheme=FMM.Scheme{FMM.DirectNeumann, FMM.FlowTangency},
                    set_initial_guess=FMM.set_unit_strength!,
                    settings_fmm=(; expansion_order = 10,
                                    leaf_size = 50,
                                    multipole_threshold = 0.4),
                    settings_solver=(;
                                    tolerance=1e-6,
                                    max_iterations=1000,
                                    history=true
                                    )
                    )

    # TODO: Do something with Das, Dbs, and TE once FMMVortexFilament kernel is
    #       implemented!

    # Set freestream at each FMM panel
    for (i, Uinf) in zip(eachindex(body_fmm.velocity), eachcol(Uinfs))
        body_fmm.velocity[i] -= body_fmm.velocity[i] # Cancel out current velocity
        body_fmm.velocity[i] += Uinf                 # Set freestream instead
    end

    # Set initial guess of strengths
    set_initial_guess(body_fmm)

    # Initiate solver
    _solver = solver(body_fmm, scheme; reuse_tree=false, settings_fmm...)

    # Call solver
    FMM.solve!(body_fmm, _solver; settings_solver...)

    return _solver

end

################################################################################
# Body TO FMM Body CONVERTER
################################################################################
const SVector = FMM.SVector
const WriteVTK = FMM.WriteVTK


bodytype2fmmkernel(::AbstractBody{<:VortexRing}) = FMM.ConstantNormalDoublet

function body2fmmbody(body)

    args = allocate_body2fmmbody(body)
    return body2fmmbody!(body, args...)

end

function body2fmmbody!(body,
                        body_fmm::FMM.UnstructuredGrid{<:Any, TF, NK, <:Any},
                        normals, controlpoints
                        ) where {TF, NK}

    # Fetch data from FMM body
    (points_vtk, meshcells_vtk,
        normals_fmm, controlpoints_fmm,
        panels_fmm,
        strengths_fmm) = (getproperty(body_fmm, p) for p in (:points, :meshcells,
                                                            :normals, :control_points,
                                                            :panels, :strengths))

    # Update VTK points
    _get_vtkpoints!(body, points_vtk)

    # Update normals and control points
    calc_normals!(body, normals)
    calc_controlpoints!(body, controlpoints, normals)

    # Format normals and control points into FMM inputs
    for (i, (N, CP)) in enumerate(zip(eachcol(normals), eachcol(controlpoints)))
        normals_fmm[i] = SVector(N...)
        controlpoints_fmm[i] = SVector(CP...)
    end

    # Update FMM panels
    for (i, meshcell) in enumerate(meshcells_vtk)

        vertices = SVector(Tuple(points_vtk[j] for j in meshcell.connectivity))
        controlpoints = controlpoints_fmm[i]
        normal = normals_fmm[i]
        strength = SVector{NK, TF}(get_strength(body, i)...)
        radius = FMM.get_radius(controlpoints, vertices...)

        panels_fmm[i] = FMM.Panel(vertices, controlpoints, normal, strength, radius)
    end

    # Update extra data of UnstructuredGrid object
    for (i, panel) in enumerate(panels_fmm)
        strengths_fmm[i] = panel.strength
    end

    return body_fmm

end

function allocate_body2fmmbody(body)

    # Identify FMM kernel from body type
    kernel_fmm = bodytype2fmmkernel(body)

    # Convert body geometry to VTK format
    points_vtk = _get_vtkpoints(body)
    meshcells_vtk = _get_vtkcells(body)

    # Allocate more memory
    normals = calc_normals(body)
    controlpoints = calc_controlpoints(body, normals)

    normals_fmm = zeros(SVector{size(normals, 1)}, size(normals, 2))
    controlpoints_fmm = zeros(SVector{size(normals, 1)}, size(controlpoints, 2))

    n_cells = length(meshcells_vtk)
    TF = eltype(eltype(points_vtk))
    NK = FMM.kernel_multiplicity(kernel_fmm)
    NS = meshcells_vtk[1].ctype.nodes
    panels_fmm = Vector{FMM.Panel{TF, NK, NS}}(undef, n_cells)

    strengths_fmm = zeros(SVector{NK, TF}, n_cells)
    potential_fmm = zeros(TF, n_cells)
    velocity_fmm = zeros(SVector{3, TF}, n_cells)

    # TODO: Update this with Kutta condition! Leaving it empty for now
    wake_points_vtk = SVector{3, TF}[]

    # Create FMM body
    body_fmm = FMM.UnstructuredGrid{kernel_fmm, TF, NK, NS}(points_vtk, meshcells_vtk,
                                                            controlpoints_fmm, normals_fmm,
                                                            strengths_fmm,
                                                            potential_fmm, velocity_fmm,
                                                            panels_fmm, wake_points_vtk)

    return (body_fmm, normals, controlpoints)
end


"""
Convert nodes into VTK points
"""
function _get_vtkpoints(body::AbstractBody; TF=Float64)

    points = zeros(SVector{3, TF}, body.nnodes)
    _get_vtkpoints!(body, points)

    return points
end

function _get_vtkpoints!(body::AbstractBody, points::AbstractVector{<:SVector})

    @assert length(points)==body.nnodes ""*
        "In-place target output `points` has a different length than needed"*
        " ($(length(points)) != $(body.nnodes))"

    for (i, X) in enumerate(eachcol(body.grid._nodes))
        points[i] = SVector(X...)
    end

end

function _get_vtkpoints!(multibody::MultiBody, points::AbstractVector{<:SVector})

    @assert length(points)==multibody.nnodes ""*
        "In-place target output `points` has a different length than needed"*
        " ($(length(points)) != $(multibody.nnodes))"

    offset = 0

    for body in multibody.bodies

        _get_vtkpoints!(body, view(points, (1:body.nnodes) .+ offset))

        offset += body.nnodes
    end

end


"""
Return cells in VTK format (WriteVTK.MeshCell)
"""
function _get_vtkcells(body::AbstractBody; offset=0)

    # Pre-allocate memory for panel calculation
    (tri_out, tricoor, quadcoor, quad_out,
            lin, ndivscells, cin) = gt.generate_getcellt_args!(body.grid)

    # Format mesh connectivity in VTK format
    meshcells = Vector{WriteVTK.MeshCell{WriteVTK.VTKCellType, SVector{3,Int64}}}(undef, body.ncells)

    for ni in 1:body.ncells

        # Fetch panel
        panel = gt.get_cell_t!(tri_out, tricoor, quadcoor, quad_out,
                                            body.grid, ni, lin, ndivscells, cin)
        # Offset point indices
        panel .+= offset

        # Convert panel to VTK format
        meshcells[ni] = WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_TRIANGLE, SVector{3}(panel...))

    end

    return meshcells
end

function _get_vtkcells(multibody::MultiBody)

    meshcells = WriteVTK.MeshCell{WriteVTK.VTKCellType, SVector{3,Int64}}[]

    offset = 0
    for body in multibody.bodies

        meshcells = vcat(meshcells, _get_vtkcells(body; offset=offset))

        offset += body.nnodes
    end

    return meshcells
end
