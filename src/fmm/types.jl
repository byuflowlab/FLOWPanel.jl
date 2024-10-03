abstract type AbstractPanels{TK,TF,NK,NS} end

abstract type AbstractKernel{M} end

abstract type AbstractRotatedKernel{M} <: AbstractKernel{M} end

abstract type AbstractUnrotatedKernel{M} <: AbstractKernel{M} end


@inline kernel_multiplicity(::AbstractKernel{M}) where M = M
