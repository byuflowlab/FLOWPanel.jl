#=##############################################################################
# DESCRIPTION
    Panel System data structure.

# AUTHORSHIP
  * Created by  : Ryan Anderson
  * Email       : Ry.M.Anderson@gmail.com
  * Date        : Feb 2026
  * License     : MIT License
=###############################################################################

"""
    PanelSystem

Compact container for surface mesh geometry, wake indexing, and associated panel
data used by higher-level coupled workflows.
"""
struct PanelSystem{TF}
    nodes::Matrix{TF}
    panels::Vector{MeshCell{VTKCellType, Vector{Int64}}}
    wake_index::Matrix{Int64}
    data::Matrix{TF}
end
