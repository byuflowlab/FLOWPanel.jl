abstract type AbstractKernel{M} end

abstract type AbstractRotatedKernel{M} <: AbstractKernel{M} end

abstract type AbstractUnrotatedKernel{M} <: AbstractKernel{M} end

abstract type AbstractPanels{TK<:AbstractKernel, TF, NK, NS} end

@inline kernel_multiplicity(::AbstractKernel{M}) where M = M

function kernel_multiplicity(type::Type)
    if type <: AbstractKernel
        return type.parameters[1]
    else
        error("Requested kernel multiplicity over invalid type $(type)")
    end
end
