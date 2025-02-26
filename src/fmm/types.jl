abstract type AbstractKernel{M} end

abstract type AbstractRotatedKernel{M} <: AbstractKernel{M} end

abstract type AbstractUnrotatedKernel{M} <: AbstractKernel{M} end

abstract type AbstractPanels{TK<:AbstractKernel, TF, NK, NS} end

@inline kernel_multiplicity(::Type{<:AbstractKernel{M}}) where M = M

