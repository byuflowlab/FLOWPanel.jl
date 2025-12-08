#=##############################################################################
# DESCRIPTION

    Nonlinear lifting line method. This method is formulated based on the 
    following references:

    * Martinez-Tossas, L. A., Allaerts, D., Branlard, E., and Churchfield, M. J. 
        (2025), ASME. J. Fluids Eng., "A Solution Method for the Filtered 
        Lifting Line Theory."

    * Cory D. Goates and Douglas F. Hunsaker (2023), Journal of Aircraft, 
        "Modern Implementation and Evaluation of Lifting-Line Theory for Complex 
        Geometries" 

    * Jackson T. Reid (2020), PhD Dissertation, "A General Approach to 
        Lifting-Line Theory, Applied to Wings With Sweep"


# AUTHORSHIP
  * Created by  : Eduardo J. Alvarez
  * Email       : Edo.AlvarezR@gmail.com
  * Date        : Nov 2025
  * License     : MIT License
=###############################################################################

#=##############################################################################
# DESCRIPTION
    Nonlinear solver for lifting line method with sweep and dihedral.


# AUTHORSHIP
  * Created by  : Eduardo J. Alvarez
  * Email       : Edo.AlvarezR@gmail.com
  * Date        : Nov 2025
  * License     : MIT License
=###############################################################################

import CSV
import DataFrames: DataFrame
import AirfoilPrep

import LaTeXStrings: @L_str

# ------------ HEADERS ---------------------------------------------------------
for header_name in ["stripwise", "liftingline", 
                    "linearsolver", "nonlinearsolver", 
                    "postprocess", "utils"
                    ]

  include(joinpath("liftingline", header_name*".jl"))

end
# ------------------------------------------------------------------------------