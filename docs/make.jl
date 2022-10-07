using Documenter
using FLOWPanel
const pnl = FLOWPanel

makedocs(
    sitename="FLOWPanel",
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
