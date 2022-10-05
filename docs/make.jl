using Documenter
import FLOWPanel
pnl = FLOWPanel

makedocs(
    sitename="FLOWPanel",
    pages = [
                "Home" => "index.md",
                "API Reference" => ["api.md", "api-elements.md"]
            ]
)
