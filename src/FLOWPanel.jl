module FLOWPanel

using StaticArrays
using LinearAlgebra: cross, norm, dot
using WriteVTK
using FLOWFMM

const ONE_OVER_4PI = 1/4/pi

include("types.jl")
include("panel.jl")
include("kernel.jl")
include("geometry.jl")
include("fmm.jl")
include("freestream.jl")
include("vtk.jl")

export Panel, PanelArray
export ConstantSource, ConstantNormalDoublet, ConstantSourceNormalDoublet
export induced, apply_freestream!
export vtk

end # END OF MODULE
