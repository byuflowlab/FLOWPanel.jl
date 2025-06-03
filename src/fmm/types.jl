abstract type AbstractKernel{M} end

abstract type AbstractRotatedKernel{M} <: AbstractKernel{M} end

abstract type AbstractUnrotatedKernel{M} <: AbstractKernel{M} end

abstract type AbstractPanels{TK<:AbstractKernel, TF, NK, NS} end

@inline kernel_multiplicity(::Type{<:AbstractKernel{M}}) where M = M

struct Source{M} <: AbstractRotatedKernel{M} end
const UniformSource = Source{1}

struct NormalDoublet{M} <: AbstractRotatedKernel{M} end
const ConstantNormalDoublet = NormalDoublet{1}

struct SourceNormalDoublet{M} <: AbstractRotatedKernel{M} end
const UniformSourceNormalDoublet = SourceNormalDoublet{2}

struct Vortex{M} <: AbstractUnrotatedKernel{M} end
const VortexRing = Vortex{1}