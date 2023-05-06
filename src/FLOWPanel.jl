"""
  Three-dimensional panel method for high-Reynolds aerodynamics.

  # AUTHORSHIP
    * Created by    : Eduardo J. Alvarez
    * Email         : Edo.AlvarezR@gmail.com
    * Date          : Jun 2018 originally as MyPanel.jl
    * License       : MIT License
"""
module FLOWPanel

export  solve, save, Uind!, phi!,
        get_ndivscells, get_ndivsnodes,
        get_cart2lin_cells, get_cart2lin_nodes,
        get_field, get_fieldval, add_field,
        calc_normals!, calc_normals,
        calc_tangents!, calc_tangents,
        calc_obliques!, calc_obliques,
        calc_controlpoints!, calc_controlpoints,
        calc_areas!, calc_areas

# ------------ GENERIC MODULES -------------------------------------------------
import Dierckx
import PyPlot
import PyPlot: @L_str
const plt = PyPlot
import LinearAlgebra as LA
import LinearAlgebra: I
import Krylov

# ------------ FLOW LAB MODULES ------------------------------------------------
# GeometricTools from https://github.com/byuflowlab/GeometricTools.jl
import GeometricTools
const gt = GeometricTools
import ImplicitAD as IAD
import ImplicitAD: ForwardDiff as FD, ReverseDiff as RD


# ------------ GLOBAL VARIABLES AND DATA STRUCTURES ----------------------------
const module_path = splitdir(@__FILE__)[1]      # Path to this module
                                                # Default path to data files
const def_data_path = joinpath(module_path, "..", "docs", "resources", "data")
                                            # Default path to airfoil data files
const def_rfl_path = joinpath(def_data_path, "airfoils")
                                                # Path to examples
const examples_path = joinpath(module_path, "..", "examples")

# const RType = Union{Float64,                    # Concrete real types
#                     Int64,
#                     # ForwardDiff.Dual{Nothing,Float64,3},
#                     # ForwardDiff.Dual{Nothing,Int64,3}
#                     }

# Discretization parameter type
const ndivstype = Union{Float64, gt.multidisctype, Nothing}

# Identity matrix
const Im = Array(1.0I, 3, 3)

# ------------ HEADERS ---------------------------------------------------------
for header_name in ["elements", "linearsolver",
                    "abstractbody", "nonliftingbody",
                    "abstractliftingbody", "liftingbody",
                    "multibody",
                    "utils", "postprocess", "monitor"
                    ]
  include("FLOWPanel_"*header_name*".jl")
end

end # END OF MODULE
