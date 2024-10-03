module FMM

using StaticArrays
using LinearAlgebra: cross, norm, dot, mul!, lu!, LU
using WriteVTK
import Krylov
import LinearOperators

# https://github.com/byuflowlab/FastMultipole
using FastMultipole

const ONE_OVER_4PI = 1/4/pi

for header_name in ["types", "panel", "kernel", "geometry", "fmm",
                    "solve", "freestream", "vtk"]

  include(joinpath("fmm", header_name*".jl"))
  
end

end # END OF MODULE
