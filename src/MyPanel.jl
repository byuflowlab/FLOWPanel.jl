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
import PyPlot
plt = PyPlot

# ------------ FLOW LAB MODULES ------------------------------------------------

# GeometricTools from https://github.com/byuflowlab/GeometricTools.jl
import GeometricTools
const gt = GeometricTools


# ------------ GLOBAL VARIABLES AND DATA STRUCTURES ----------------------------
const module_path = splitdir(@__FILE__)[1]      # Path to this module
const data_path = joinpath(module_path, "../data")  # Path to data files

# Structure of implemented fields
const FIELDS = Dict(
  # Declare fields as follow:
  # "MyField" => Dict( "field_type"  => "scalar" or "vector",
  #                     "entry_type" => "node" or "cell",
  #                     "field_data" => data)
)


# ------------ HEADERS ---------------------------------------------------------
for header_name in ["solver", "abstractbody", "utils"]
  include("MyPanel_"*header_name*".jl")
end


end # END OF MODULE
