#=##############################################################################
# DESCRIPTION
    Non-lifting paneled body type definition.
# AUTHORSHIP
  * Author    : Eduardo J. Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Jun 2018
  * License   : AGPL-3.0
=###############################################################################


################################################################################
# NON-LIFTING BODY TYPE
################################################################################
struct NonLiftingBody <: AbstractBody

  # User inputs
  grid::gt.GridTriangleSurface              # Paneled geometry

  # Properties
  nnodes::Int64                             # Number of nodes
  ncells::Int64                             # Number of cells
  fields::Array{String, 1}                  # Available fields (solutions)

  # Internal variables

  NonLiftingBody( grid,
                  nnodes=grid.nnodes, ncells=grid.ncells,
                  fields=Array{String,1}()
         ) = new( grid,
                  nnodes, ncells,
                  fields
         )
end


"""
  `get_field(self::NonLiftingBody, field_name::String)`

Returns the requested field.
"""
function get_field(self::NonLiftingBody, field_name::String)
  if !(field_name in self.fields)
    error("Field $field_name not found! Available fields: $(self.fields)")
  end

  return self.grid.field[field_name]
end

"""
  `get_field(self::NonLiftingBody, field_name::String, i::Int64)`

Returns the requested field value. Give it `_check=false` to skip checking logic
for faster operation.
"""
function get_field_val(self::NonLiftingBody, field_name::String, i::Int64;
                      _check::Bool=true)
  if _check
    if i<1
      error("Invalid index $i.")
    elseif FIELDS[field_name]["entry_type"]=="node" && i>self.nnodes
      error("Invalid index $i. Maximum is $(self.nnodes).")
    elseif FIELDS[field_name]["entry_type"]=="cell" && i>self.ncells
      error("Invalid index $i. Maximum is $(self.ncells).")
    end
  end

  return self.grid.field[field_name]["field_data"][i]
end


"""
  `add_field(self::NonLiftingBody, field_name::String, field_data)`

Adds a new field to the body.
"""
function add_field(self::NonLiftingBody, field_name::String, field_data)

  # ERROR CASES
  if !(field_name in keys(FIELDS))
    error("Invalid field $field_name. Implemeted fields are: $(keys(FIELDS))")
  end

  gt.add_field(self.grid, field_name, FIELDS[field_name]["field_type"],
                field_data, FIELDS[field_name]["field_type"])

  if !(field_name in self.fields)
    push!(self.fields, field_name)
  end

  nothing
end
##### INTERNAL FUNCTIONS  ######################################################


##### END OF NON-LIFTING BODY ##################################################
