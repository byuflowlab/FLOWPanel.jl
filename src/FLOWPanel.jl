module FLOWPanel

using StaticArrays
using LinearAlgebra: cross, norm, dot, mul!
using WriteVTK
using FLOWFMM

const ONE_OVER_4PI = 1/4/pi

include("types.jl")
include("panel.jl")
include("kernel.jl")
include("geometry.jl")
include("fmm.jl")
include("solve.jl")
include("freestream.jl")
include("vtk.jl")

export Panel, PanelArray
export ConstantSource, ConstantNormalDoublet, ConstantSourceNormalDoublet
export induced, apply_freestream!, reset_potential_velocity!
export FlowTangency, LUDecomposition, FastJacobi, FastGaussSeidel, solve!, grid_2_panels_strength!
export vtk

end # END OF MODULE
