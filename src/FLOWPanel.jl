"""
  Three-dimensional panel method for high-Reynolds aerodynamics.

  # AUTHORSHIP
    * Author    : Eduardo J. Alvarez
    * Email     : Edo.AlvarezR@gmail.com
    * Created   : Jun 2018 originally as MyPanel.jl
    * License   : MIT License
"""
module FLOWPanel

# ------------ GENERIC MODULES -------------------------------------------------
import Dierckx
import ForwardDiff
import PyPlot
const plt = PyPlot

# ------------ FLOW LAB MODULES ------------------------------------------------
# GeometricTools from https://github.com/byuflowlab/GeometricTools.jl
import GeometricTools
const gt = GeometricTools


# ------------ GLOBAL VARIABLES AND DATA STRUCTURES ----------------------------
const module_path = splitdir(@__FILE__)[1]      # Path to this module
                                                # Default path to data files
const def_data_path = joinpath(module_path, "..", "docs", "data")
                                            # Default path to airfoil data files
const def_rfl_path = joinpath(def_data_path, "airfoils")

const RType = Union{Float64,                    # Concrete real types
                    Int64,
                    ForwardDiff.Dual{Nothing,Float64,3},
                    ForwardDiff.Dual{Nothing,Int64,3}
                    }

# Discretization parameter type
const ndivstype = Union{Float64, gt.multidisctype, Nothing}

# ------------ HEADERS ---------------------------------------------------------
for header_name in ["elements", "G", "abstractbody", "utils"]
  include("FLOWPanel_"*header_name*".jl")
end



# ------------ USEFUL FUNCTIONS ------------------------------------------------
dot(A, B) = sum(a*b for (a,b) in zip(A, B))
norm(A) = sqrt(mapreduce(x->x^2, +, A))
function cross!(out, A, B)
    out[1] = A[2]*B[3] - A[3]*B[2]
    out[2] = A[3]*B[1] - A[1]*B[3]
    out[3] = A[1]*B[2] - A[2]*B[1]
end
cross(A,B) = (out = zeros(3); cross!(out, A, B); return out)
mean(xs) = sum(xs)/length(xs)

end # END OF MODULE
