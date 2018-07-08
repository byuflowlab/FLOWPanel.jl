#=##############################################################################
# DESCRIPTION
    Abstract body type definition.
# AUTHORSHIP
  * Author    : Eduardo J. Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Jun 2018
  * License   : AGPL-3.0
=###############################################################################


################################################################################
# ABSTRACT BODY TYPE
################################################################################
"""
  Implementations of AbstractBody are expected to have the following fields.
  * ` `             : .

  and the following functions
  ```julia
  ```
"""
abstract type AbstractBody end


for header_name in ["nonliftingbody"]
  include("MyPanel_"*header_name*".jl")
end

# Implementations of AbstractGrid
const BodyTypes = Union{NonLiftingBody}



##### COMMON FUNCTIONS  ########################################################


##### END OF ABSTRACT BODY #####################################################
