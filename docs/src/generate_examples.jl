import FLOWPanel as pnl

module_path = splitdir(@__FILE__)[1]      # Path to this module
output_path = joinpath(module_path, "examples") # Path where to store  examples as markdown

include("generate_examples_sweptwing.jl")
include("generate_examples_centerbody.jl")
include("generate_examples_duct.jl")
include("generate_examples_blendedwingbody.jl")
include("generate_examples_cessna.jl")
include("generate_examples_ll_weber.jl")
include("generate_examples_ll_a50k27.jl")
