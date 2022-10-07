using Documenter
using FLOWPanel
const pnl = FLOWPanel

makedocs(
    sitename = "FLOWPanel.jl",
    format = Documenter.HTML(;
                                sidebar_sitename = false,
                                assets = ["assets/favicon.ico"],
                            ),
    pages = [
                "Home"              => "index.md",
                "Potential Flow"    => "potentialflow.md",
                "Elements"          => [
                                        "elements/paneldefinition.md",
                                        "Constant Source"  => "elements/constantsource.md",
                                        "Constant Doublet" => "elements/constantdoublet.md",
                                        "Semi-Infinite Doublet" => "elements/semiinfdoublet.md",
                                        "Non-Planar Semi-Infinite Doublet" => "elements/semiinfnonplanardoublet.md",
                                       ]
                # "API Reference"     => ["api.md",
                #                         "api-elements.md",
                #                         "api-abstractbody.md"
                #                        ]
            ]
)



# # Documenter can also automatically deploy documentation to gh-pages.
# # See "Hosting Documentation" and deploydocs() in the Documenter manual
# # for more information.
deploydocs(
    repo = "github.com/byuflowlab/FLOWPanel.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
    # devbranch = "main"
)
