#=##############################################################################
# DESCRIPTION
    Definition of 2D body types

# AUTHORSHIP
  * Created by  : Eduardo J. Alvarez
  * Email       : Edo.AlvarezR@gmail.com
  * Date        : Jan 2025
  * License     : MIT License
=###############################################################################

################################################################################
# 2D BODY TYPE
################################################################################

struct Body2D{E, N} <: AbstractBody{E, N}

    # User inputs
    grid::gt.GridSurface2D                       # Paneled geometry

    # Properties
    nnodes::Int                               # Number of nodes
    ncells::Int                               # Number of cells
    fields::Array{String, 1}                  # Available fields (solutions)
    Oaxis::Array{<:Number,2}                  # Coordinate system of original grid
    O::Array{<:Number,1}                      # Position of CS of original grid

    # Internal variables
    strength::Matrix                          # strength[i,j] is the stength of the i-th panel with the j-th element type
    CPoffset::Float64                         # Control point offset in normal direction
    kerneloffset::Float64                     # Kernel offset to avoid singularities
    kernelcutoff::Float64                     # Kernel cutoff to avoid singularities
    characteristiclength::Function            # Characteristic length of each panel


    function Body2D{E, N}(
                            grid;
                            nnodes=grid.nnodes, ncells=grid.ncells,
                            fields=Array{String,1}(),
                            Oaxis=Array(1.0I, 3, 3), O=zeros(3),
                            strength=zeros(grid.ncells, N),
                            CPoffset=1e-14,
                            kerneloffset=1e-8,
                            kernelcutoff=1e-14,
                            characteristiclength=characteristiclength_unitary
                            ) where {E, N}

        return new(grid,
                    nnodes, ncells,
                    fields,
                    Oaxis, O,
                    strength,
                    CPoffset, kerneloffset, kernelcutoff, characteristiclength)
    end
end


save(body::Body2D, args...; optargs...) = save_base(body, args...; optargs...)
