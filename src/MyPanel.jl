"""
  Three-dimensional panel method for high-Reynolds aerodynamics.

  # AUTHORSHIP
    * Author    : Eduardo J. Alvarez
    * Email     : Edo.AlvarezR@gmail.com
    * Created   : Jun 2018
    * License   : AGPL-3.0
"""
module MyPanel

# ------------ GENERIC MODULES -------------------------------------------------
import Dierckx
import PyPlot; const plt = PyPlot
import ForwardDiff
using Statistics
using LinearAlgebra

# ------------ FLOW LAB MODULES ------------------------------------------------

# GeometricTools from https://github.com/byuflowlab/GeometricTools.jl
import GeometricTools
const gt = GeometricTools


# ------------ GLOBAL VARIABLES AND DATA STRUCTURES ----------------------------
const module_path = splitdir(@__FILE__)[1]      # Path to this module
                                                # Default path to data files
const def_data_path = joinpath(module_path, "../data")
                                            # Default path to airfoil data files
const def_rfl_path = joinpath(def_data_path, "airfoils")

const RType = Union{Float64,                    # Concrete real types
                    Int64,
                    ForwardDiff.Dual{Nothing,Float64,3},
                    ForwardDiff.Dual{Nothing,Int64,3}
                    }

# Structure of implemented fields
const FIELDS = Dict(
  # Declare fields as follow:
  # "MyField" => Dict( "field_type"  => "scalar" or "vector",
  #                     "entry_type" => "node", "cell", or "system")

  # Flag indicating whether the paneled body has been solved
  "solved"    => Dict(  "field_type"  => "scalar",
                        "entry_type"  => "system"),

  # Freestream velocity at every control point
  "Vinf"      => Dict(  "field_type"  => "vector",
                        "entry_type"  => "cell"),

  # Constant source strength at every panel
  "sigma"     => Dict(  "field_type"  => "scalar",
                        "entry_type"  => "cell"),

  # Constant doublet strength at every panel
  "mu"        => Dict(  "field_type"  => "scalar",
                        "entry_type"  => "cell"),

  # Vortex ring strength at every panel
  "Gamma"     => Dict(  "field_type"  => "scalar",
                        "entry_type"  => "cell"),

  # Freestream direction at every trailing edge point
  "D"         => Dict(  "field_type"  => "vector",
                        "entry_type"  => "system"),

  # Strength of wake elements
  "Gammawake" => Dict(  "field_type"  => "scalar",
                        "entry_type"  => "system"),
)

# Discretization parameter type
const ndivstype = Union{Float64, gt.multidisctype, Nothing}

# ------------ HEADERS ---------------------------------------------------------
for header_name in ["solver", "abstractbody", "utils"]
  include("MyPanel_"*header_name*".jl")
end


end # END OF MODULE
