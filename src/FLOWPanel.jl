module FLOWPanel

using StaticArrays
using LinearAlgebra: cross, norm, dot, mul!
using WriteVTK
using FastMultipole
import Krylov
import LinearOperators

const ONE_OVER_4PI = 1/4/pi
const TEST_FLAG = true


if TEST_FLAG
    FastMultipole.direct!(stuff::String) = "hello"
else
    FastMultipole.direct!(stuff::String) = "world"
end

include("types.jl")
include("panel.jl")
include("kernel.jl")
include("geometry.jl")
include("fmm.jl")
include("solve.jl")
include("freestream.jl")
include("vtk.jl")

export Panel, PanelArray
export ConstantSource, ConstantNormalDoublet, ConstantSourceNormalDoublet, VortexRing
export induced, apply_freestream!, reset_potential_velocity!
export FlowTangency, LUDecomposition, FastGaussSeidel, solve!, grid_2_panels_strength!
export vtk

end # END OF MODULE
