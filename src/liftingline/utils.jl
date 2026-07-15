#=##############################################################################
# DESCRIPTION
    Lifting line auxiliary functions

# AUTHORSHIP
  * Created by  : Eduardo J. Alvarez
  * Email       : Edo.AlvarezR@gmail.com
  * Date        : Nov 2025
  * License     : MIT License
=###############################################################################

function calc_area(ll::LiftingLine{R}) where {R}

    nodes = ll.grid.nodes
    linearindices = ll.linearindices

    area = zero(R)

    for ei in 1:ll.nelements

        TEa = view(nodes, :, linearindices[ei, 1])    # TE at a-side
        TEb = view(nodes, :, linearindices[ei+1, 1])  # TE at b-side
        LEa = view(nodes, :, linearindices[ei, 2])    # LE at a-side
        LEb = view(nodes, :, linearindices[ei+1, 2])  # LE at b-side

        area += 0.5*norm(cross(TEa-LEa, LEb-LEa))
        area += 0.5*norm(cross(LEb-TEb, TEa-TEb))

    end

    return area
end

function calc_noseposition(ll::LiftingLine)

    nodes = ll.grid.nodes
    lin = ll.linearindices

    # LE on b-side of root element
    return nodes[:, lin[ll.root_i, 2]] 

end


function save(self::LiftingLine{R}, filename::AbstractString; 
                    format="vtk", 
                    horseshoe_suffix="_horseshoe",
                    midpoint_suffix="_midp", 
                    controlpoint_suffix="_cp", 
                    planar_suffix="_planar", 
                    wake_length_factor=0.25,
                    debug=false,
                    optargs...) where R

    str = ""

    # Function for fetching values out of Dual numbers
    value(x::Number) = x
    value(x::FD.Dual) = FD.value(x)
    value(x::AbstractArray) = x
    value(x::AbstractArray{<:FD.Dual}) = FD.value.(x)

    # ------------- OUTPUT HORSESHOES ------------------------------------------

    npoints = size(self.horseshoes, 2)                  # Number of points in a horseshoe
    nhorseshoes = size(self.horseshoes, 3)

    # Flatten the tensor of horseshoe point into a matrix of points
    horseshoes = self.horseshoes
    horseshoes = reshape(horseshoes, ( size(horseshoes, 1), npoints*nhorseshoes ))

    # Connect the horseshoe points (0-indexed for VTK format)
    lines = [ collect( (0:npoints-1) .+ (ei-1)*npoints ) for ei in 1:nhorseshoes]

    # Determine a characteristic length for the wake
    wake_length = norm( maximum(horseshoes; dims=2) - minimum(horseshoes; dims=2) )
    wake_length *= wake_length_factor

    # Add fictitious wake end points to the horseshoes
    horseshoes = hcat(horseshoes, zeros(3, 2*nhorseshoes))

    for ei in 1:nhorseshoes

        # TE points
        Ap = horseshoes[:, lines[ei][1]+1]
        Bp = horseshoes[:, lines[ei][end]+1]

        # Indices of wake end points
        ai = npoints*nhorseshoes + 2*(ei-1) + 1
        bi = npoints*nhorseshoes + 2*(ei-1) + 2

        # Add wake end points
        horseshoes[:, ai] .= Ap + wake_length*self.Dinfs[:, 1, ei]
        horseshoes[:, bi] .= Bp + wake_length*self.Dinfs[:, 2, ei]

        # Add points to the VTK line definition
        lines[ei] = vcat(ai-1, lines[ei], bi-1)

    end

    # Format the points as a vector of vectors 
    horseshoes = eachcol(value(horseshoes))

    horseshoes_data = [
                            Dict(   "field_name"  => "Gamma",
                                    "field_type"  => "scalar",
                                    "field_data"  => value(self.Gammas))

                            Dict(   "field_name"  => "sigma",
                                    "field_type"  => "scalar",
                                    "field_data"  => value(self.sigmas))

                            Dict(   "field_name"  => "angleofattack",
                                    "field_type"  => "scalar",
                                    "field_data"  => value(self.aoas))
                                    
                            Dict(   "field_name"  => "claero",
                                    "field_type"  => "scalar",
                                    "field_data"  => value(self.claeros))
                        ]

    str *= gt.generateVTK(filename*horseshoe_suffix, horseshoes; cells=lines,
                                cell_data=horseshoes_data, num=debug ? 0 : nothing,
                                override_cell_type=4, optargs...)


    # ------------- OUTPUT EFFECTIVE LIFTING LINE HORSESHOES -------------------
    if debug
        for ei in 1:self.nelements

            # Flatten the tensor of horseshoe point into a matrix of points
            horseshoes = view(self.effective_horseshoes, :, :, :, ei)
            horseshoes = reshape(horseshoes, ( size(horseshoes, 1), npoints*nhorseshoes ))

            # Connect the horseshoe points (0-indexed for VTK format)
            lines = [ collect( (0:npoints-1) .+ (ni-1)*npoints ) for ni in 1:nhorseshoes]

            # Determine a characteristic length for the wake
            wake_length = norm( maximum(horseshoes; dims=2) - minimum(horseshoes; dims=2) )
            wake_length *= wake_length_factor

            # Add fictitious wake end points to the horseshoes
            horseshoes = hcat(horseshoes, zeros(3, 2*nhorseshoes))

            for ni in 1:nhorseshoes

                # TE points
                Ap = horseshoes[:, lines[ni][1]+1]
                Bp = horseshoes[:, lines[ni][end]+1]

                # Indices of wake end points
                ai = npoints*nhorseshoes + 2*(ni-1) + 1
                bi = npoints*nhorseshoes + 2*(ni-1) + 2

                # Add wake end points
                horseshoes[:, ai] .= Ap + wake_length*self.Dinfs[:, 1, ni]
                horseshoes[:, bi] .= Bp + wake_length*self.Dinfs[:, 2, ni]

                # Add points to the VTK line definition
                lines[ni] = vcat(ai-1, lines[ni], bi-1)

            end

            # Format the points as a vector of vectors 
            horseshoes = eachcol(value(horseshoes))

            this_str = gt.generateVTK(filename*horseshoe_suffix, horseshoes; 
                                        cells=lines, num=ei,
                                        override_cell_type=4, optargs...)

            if ei==1
                str *= replace(this_str, ".$ei." => ".*.")
            end

            # Output associated midpoint
            midpoints = [value(self.midpoints[:, ei])]

            midpoints_data = [
                                    Dict(   "field_name"  => "tangent",
                                            "field_type"  => "vector",
                                            "field_data"  => [value(self.tangents[:, ei])]),

                                    Dict(   "field_name"  => "span",
                                            "field_type"  => "vector",
                                            "field_data"  => [value(self.spans[:, ei])]),

                                    Dict(   "field_name"  => "normal",
                                            "field_type"  => "vector",
                                            "field_data"  => [value(self.normals[:, ei])]),

                                    Dict(   "field_name"  => "swepttangent",
                                            "field_type"  => "vector",
                                            "field_data"  => [value(self.swepttangents[:, ei])]),

                                    Dict(   "field_name"  => "line",
                                            "field_type"  => "vector",
                                            "field_data"  => [value(self.lines[:, ei])]),

                                    Dict(   "field_name"  => "sweptnormal",
                                            "field_type"  => "vector",
                                            "field_data"  => [value(self.sweptnormals[:, ei])]),

                                    Dict(   "field_name"  => "angleofattack",
                                            "field_type"  => "scalar",
                                            "field_data"  => [value(self.aoas[ei])]),

                                    Dict(   "field_name"  => "U",
                                            "field_type"  => "vector",
                                            "field_data"  => [value(self.Us[:, ei])])
                                ]

            this_str = gt.generateVTK(filename*midpoint_suffix, midpoints;
                                        num=ei, 
                                        point_data=midpoints_data, optargs...)

            if ei==1
                str *= replace(this_str, ".$ei." => ".*.")
            end

        end


    end

    # ------------- OUTPUT MIDPOINTS -------------------------------------------
    midpoints = eachcol(value(self.midpoints))

    midpoints_data = [
                            Dict(   "field_name"  => "tangent",
                                    "field_type"  => "vector",
                                    "field_data"  => eachcol(value(self.tangents))),

                            Dict(   "field_name"  => "span",
                                    "field_type"  => "vector",
                                    "field_data"  => eachcol(value(self.spans))),

                            Dict(   "field_name"  => "normal",
                                    "field_type"  => "vector",
                                    "field_data"  => eachcol(value(self.normals))),

                            Dict(   "field_name"  => "swepttangent",
                                    "field_type"  => "vector",
                                    "field_data"  => eachcol(value(self.swepttangents))),

                            Dict(   "field_name"  => "line",
                                    "field_type"  => "vector",
                                    "field_data"  => eachcol(value(self.lines))),

                            Dict(   "field_name"  => "sweptnormal",
                                    "field_type"  => "vector",
                                    "field_data"  => eachcol(value(self.sweptnormals))),
                                    
                            Dict(   "field_name"  => "angleofattack",
                                    "field_type"  => "scalar",
                                    "field_data"  => value(self.aoas)),

                            Dict(   "field_name"  => "U",
                                    "field_type"  => "vector",
                                    "field_data"  => eachcol(value(self.Us)))
                        ]

    str *= gt.generateVTK(filename*midpoint_suffix, midpoints; 
                                num = debug ? 0 : nothing,
                                point_data=midpoints_data, optargs...)

    # ------------- OUTPUT CONTROL POINTS --------------------------------------
    controlpoints = eachcol(value(self.controlpoints))

    controlpoints_data = [
                            Dict(   "field_name"  => "normal",
                                    "field_type"  => "vector",
                                    "field_data"  => eachcol(value(self.normals))),

                            Dict(   "field_name"  => "angleofattack",
                                    "field_type"  => "scalar",
                                    "field_data"  => value(self.aoas))
                        ]

    str *= gt.generateVTK(filename*controlpoint_suffix, controlpoints; 
                                point_data=controlpoints_data, optargs...)

    #  ------------- OUTPUT PLANAR GEOMETRY ------------------------------------
    if R <: FD.Dual

        # Fetch original Dual data
        nodes_org = self.grid.nodes
        fields_org = Dict(field_name => field["field_data"] for (field_name, field) in self.grid.field)

        # Replace Dual data with Real data
        self.grid.nodes = value(nodes_org)

        for (field_name, field_data) in fields_org

            field_type = self.grid.field[field_name]["field_type"]
            self.grid.field[field_name]["field_data"] = field_type == "vector" ? value.(field_data) : value(field_data)

        end

        # Save grid with Real data
        str *= gt.save(self.grid, filename*planar_suffix; format, optargs...)

        # Restore Dual data
        self.grid.nodes = nodes_org

        for (field_name, field_data) in fields_org
            self.grid.field[field_name]["field_data"] = field_data
        end

    else
        str *= gt.save(self.grid, filename*planar_suffix; format, optargs...)
    end

    return str
end



#=
THE FOLLOWING FUNCTIONS WHERE COPY/PASTED FROM FLOWMath.jl AND MODIFIED TO IMPROVE 
PERFORMANCE FOR BIG DATASETS 
=#

(spline::math.Akima)(y::AbstractVector, x::AbstractVector) = broadcast!(spline, y, x)

akima!(x, y, xpt, ypt; delta=0.0, eps=1e-30) = math.Akima(x, y, delta, eps)(ypt, xpt)

function interp2d!(interp1d!, xdata, ydata, fdata, 
                        xpt::AbstractArray{Tx}, ypt::AbstractArray{Ty}; 
                        yinterp = Array{promote_type(Tx, Ty)}(undef, length(ydata), length(xpt)),
                        output = Array{promote_type(Tx, Ty)}(undef, length(xpt), length(ypt)),
                        optargs...
                        ) where {Tx, Ty}
    
    for i in 1:ny
         interp1d!(xdata, view(fdata, :, i), xpt, view(yinterp, i, :); optargs...)
    end
    for i in 1:nxpt
         interp1d!(ydata, view(yinterp, :, i), ypt, view(output, i, :); optargs...)
    end

    return output
end


function interp3d!(interp1d!, xdata, ydata, zdata, fdata, 
                        xpt::AbstractArray{Tx}, ypt::AbstractArray{Ty}, 
                        zpt::AbstractArray{Tz};
                        zinterp = Array{promote_type(Tx, Ty, Tz)}(undef, length(zdata), length(xpt), length(ypt)),
                        output = Array{promote_type(Tx, Ty, Tz)}(undef, length(xpt), length(ypt), length(zpt)),
                        optargs...
                        ) where {Tx, Ty, Tz}

    nz = length(zdata)
    nxpt = length(xpt)
    nypt = length(ypt)
    nzpt = length(zpt)
    
    for i in 1:nz
         interp2d!(interp1d!, xdata, ydata, view(fdata, :, :, i), xpt, ypt; output=view(zinterp, i, :, :), optargs...)
    end
    for j in 1:nypt
        for i in 1:nxpt
            interp1d!(zdata, view(zinterp, :, i, j), zpt, view(output, i, j, :))
        end
    end

    return output
end


function interp4d!(interp1d!, xdata, ydata, zdata, tdata, fdata, 
                        xpt::AbstractArray{Tx}, ypt::AbstractArray{Ty}, 
                        zpt::AbstractArray{Tz}, tpt::AbstractArray{Tt};
                        tinterp = Array{promote_type(Tx, Ty, Tz, Tt)}(undef, length(tdata), length(xpt), length(ypt), length(zpt)),
                        output = Array{promote_type(Tx, Ty, Tz, Tt)}(undef, length(xpt), length(ypt), length(zpt), length(tpt)),
                        optargs...
                        ) where {Tx, Ty, Tz, Tt}
    nt = length(tdata)
    nxpt = length(xpt)
    nypt = length(ypt)
    nzpt = length(zpt)
    ntpt = length(tpt)

    for i in 1:nt
        interp3d!(
            interp1d!, xdata, ydata, zdata, view(fdata, :, :, :, i), xpt, ypt, zpt;
            output = view(tinterp, i, :, :, :), optargs...
        )
    end
    for k in 1:nzpt
        for j in 1:nypt
            for i in 1:nxpt
                interp1d!(tdata, view(tinterp, :, i, j, k), tpt, view(output, i, j, k, :))
            end
        end
    end

    return output
end

# COPIED FROM DuctAPE.jl AND MODIFIED
function interp5d!(interp1d!, x1data, x2data, x3data, x4data, x5data, fdata, 
                        x1pt::AbstractArray{Tx1}, x2pt::AbstractArray{Tx2}, 
                        x3pt::AbstractArray{Tx3}, x4pt::AbstractArray{Tx4}, 
                        x5pt::AbstractArray{Tx5};
                        x5interp = Array{promote_type(Tx1, Tx2, Tx3, Tx4, Tx5)}(undef, length(x5data), length(x1pt), length(x2pt), length(x3pt), length(x4pt)),
                        output = Array{promote_type(Tx1, Tx2, Tx3, Tx4, Tx5)}(undef, length(x1pt), length(x2pt), length(x3pt), length(x4pt), length(x5pt)),
                        optargs...
                        ) where {Tx1, Tx2, Tx3, Tx4, Tx5}
    nd = length(x5data)
    nx1pt = length(x1pt)
    nx2pt = length(x2pt)
    nx3pt = length(x3pt)
    nx4pt = length(x4pt)
    nx5pt = length(x5pt)

    for i in 1:nd
        interp4d!(
            interp1d,
            x1data,
            x2data,
            x3data,
            x4data,
            view(fdata, :, :, :, :, i),
            x1pt,
            x2pt,
            x3pt,
            x4pt;
            output = view(x5interp, i, :, :, :),
            optargs...
        )
    end
    for ell in 1:nx4pt
        for k in 1:nx3pt
            for j in 1:nx2pt
                for i in 1:nx1pt
                    interp1d!(
                        x5data, x5interp[:, i, j, k, ell], x5pt, view(output, i, j, k, ell, :)
                    )
                end
            end
        end
    end

    return output
end



# COPIED FROM DuctAPE.jl WITHOUT MODIFICATIONS
function interp5d(
    interp1d, x1data, x2data, x3data, x4data, x5data, fdata, x1pt, x2pt, x3pt, x4pt, x5pt
)
    nd = length(x5data)
    nx1pt = length(x1pt)
    nx2pt = length(x2pt)
    nx3pt = length(x3pt)
    nx4pt = length(x4pt)
    nx5pt = length(x5pt)

    R = promote_type(eltype(x1pt), eltype(x2pt), eltype(x3pt), eltype(x4pt), eltype(x5pt))
    x5interp = Array{R}(undef, nd, nx1pt, nx2pt, nx3pt, nx4pt)
    output = Array{R}(undef, nx1pt, nx2pt, nx3pt, nx4pt, nx5pt)

    for i in 1:nd
        x5interp[i, :, :, :] .= math.interp4d(
            interp1d,
            x1data,
            x2data,
            x3data,
            x4data,
            fdata[:, :, :, :, i],
            x1pt,
            x2pt,
            x3pt,
            x4pt,
        )
    end
    for ell in 1:nx4pt
        for k in 1:nx3pt
            for j in 1:nx2pt
                for i in 1:nx1pt
                    output[i, j, k, ell, :] .= interp1d(
                        x5data, x5interp[:, i, j, k, ell], x5pt
                    )
                end
            end
        end
    end

    return output
end