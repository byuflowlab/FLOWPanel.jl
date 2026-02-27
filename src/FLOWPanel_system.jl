#=##############################################################################
# DESCRIPTION
    Panel System data structure.

# AUTHORSHIP
  * Created by  : Ryan Anderson
  * Email       : Ry.M.Anderson@gmail.com
  * Date        : Feb 2026
  * License     : MIT License
=###############################################################################

struct PanelSystem{TF}
    nodes::Matrix{TF}
    panels::Vector{MeshCell{VTKCellType, Vector{Int64}}}
    wake_index::Matrix{Int64}
    data::Matrix{TF}
end

