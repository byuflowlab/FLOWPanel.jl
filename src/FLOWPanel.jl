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
import LinearAlgebra as LA
import LinearAlgebra: I
import Krylov
import NLsolve
import NonlinearSolve
import SimpleNonlinearSolve
import Requires: @require

# ------------ FLOW LAB MODULES ------------------------------------------------
# GeometricTools from https://github.com/byuflowlab/GeometricTools.jl
import GeometricTools
const gt = GeometricTools
import GeometricTools: Meshes
import ImplicitAD as IAD
import ImplicitAD: ForwardDiff as FD, ReverseDiff as RD
import FLOWMath as math

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

# Shedding matrix for a RigidWakeBody without shedding
const noshedding = zeros(Int, 6, 0)

# ------------ HEADERS ---------------------------------------------------------
for header_name in ["elements", "linearsolver",
                    "abstractbody", "nonliftingbody",
                    "abstractliftingbody", "liftingbody",
                    "multibody",
                    "liftingline",
                    "utils", "postprocess",
                    # "fmm"
                    ]
  include("FLOWPanel_"*header_name*".jl")
end



function __init__()

    # Conditionally load monitors if PythonPlot is available
    try
        @require PythonPlot="274fc56d-3b97-40fa-a1cd-1b4a50311bf9" begin

            import .PythonPlot as plt
            import .PythonPlot: pyconvert

            for header_name in ["monitor"]
              include("FLOWPanel_"*header_name*".jl")
            end

        end

    catch e
        @warn "PythonPlot is not available; monitors will not be loaded"
    end


    # Conditionally load FMM solver if FastMultipole is available
    try
        @require FastMultipole="ce07d0d3-2b9f-49ba-89eb-12c800257c85" begin

            for header_name in ["fmm"]
              include("FLOWPanel_"*header_name*".jl")
            end

        end

    catch e
        @warn "FastMultipole is not available; FMM solvers will not be loaded"
    end

end



end # END OF MODULE
