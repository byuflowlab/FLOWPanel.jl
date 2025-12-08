#=##############################################################################
# DESCRIPTION
    Lifting line auxiliary functions

# AUTHORSHIP
  * Created by  : Eduardo J. Alvarez
  * Email       : Edo.AlvarezR@gmail.com
  * Date        : Nov 2025
  * License     : MIT License
=###############################################################################



function save(self::LiftingLine, filename::AbstractString; 
                    format="vtk", 
                    horseshoe_suffix="_horseshoe",
                    midpoint_suffix="_midp", 
                    controlpoint_suffix="_cp", 
                    planar_suffix="_planar", 
                    wake_length_factor=0.25,
                    debug=false,
                    optargs...)

    str = ""

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
    horseshoes = eachcol(horseshoes)

    horseshoes_data = [
                            Dict(   "field_name"  => "Gamma",
                                    "field_type"  => "scalar",
                                    "field_data"  => self.Gammas)

                            Dict(   "field_name"  => "sigma",
                                    "field_type"  => "scalar",
                                    "field_data"  => self.sigmas)

                            Dict(   "field_name"  => "angleofattack",
                                    "field_type"  => "scalar",
                                    "field_data"  => self.aoas)
                                    
                            Dict(   "field_name"  => "claero",
                                    "field_type"  => "scalar",
                                    "field_data"  => self.claeros)
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
            horseshoes = eachcol(horseshoes)

            str *= gt.generateVTK(filename*horseshoe_suffix, horseshoes; 
                                        cells=lines, num=ei,
                                        override_cell_type=4, optargs...)

            # Output associated midpoint
            midpoints = [self.midpoints[:, ei]]

            midpoints_data = [
                                    Dict(   "field_name"  => "tangent",
                                            "field_type"  => "vector",
                                            "field_data"  => [self.tangents[:, ei]]),

                                    Dict(   "field_name"  => "span",
                                            "field_type"  => "vector",
                                            "field_data"  => [self.spans[:, ei]]),

                                    Dict(   "field_name"  => "normal",
                                            "field_type"  => "vector",
                                            "field_data"  => [self.normals[:, ei]]),

                                    Dict(   "field_name"  => "swepttangent",
                                            "field_type"  => "vector",
                                            "field_data"  => [self.swepttangents[:, ei]]),

                                    Dict(   "field_name"  => "line",
                                            "field_type"  => "vector",
                                            "field_data"  => [self.lines[:, ei]]),

                                    Dict(   "field_name"  => "sweptnormal",
                                            "field_type"  => "vector",
                                            "field_data"  => [self.sweptnormals[:, ei]]),

                                    Dict(   "field_name"  => "angleofattack",
                                            "field_type"  => "scalar",
                                            "field_data"  => [self.aoas[ei]]),

                                    Dict(   "field_name"  => "U",
                                            "field_type"  => "vector",
                                            "field_data"  => [self.Us[:, ei]])
                                ]

            str *= gt.generateVTK(filename*midpoint_suffix, midpoints;
                                        num=ei, 
                                        point_data=midpoints_data, optargs...)

        end
    end

    # ------------- OUTPUT MIDPOINTS -------------------------------------------
    midpoints = eachcol(self.midpoints)

    midpoints_data = [
                            Dict(   "field_name"  => "tangent",
                                    "field_type"  => "vector",
                                    "field_data"  => eachcol(self.tangents)),

                            Dict(   "field_name"  => "span",
                                    "field_type"  => "vector",
                                    "field_data"  => eachcol(self.spans)),

                            Dict(   "field_name"  => "normal",
                                    "field_type"  => "vector",
                                    "field_data"  => eachcol(self.normals)),

                            Dict(   "field_name"  => "swepttangent",
                                    "field_type"  => "vector",
                                    "field_data"  => eachcol(self.swepttangents)),

                            Dict(   "field_name"  => "line",
                                    "field_type"  => "vector",
                                    "field_data"  => eachcol(self.lines)),

                            Dict(   "field_name"  => "sweptnormal",
                                    "field_type"  => "vector",
                                    "field_data"  => eachcol(self.sweptnormals)),
                                    
                            Dict(   "field_name"  => "angleofattack",
                                    "field_type"  => "scalar",
                                    "field_data"  => self.aoas),

                            Dict(   "field_name"  => "U",
                                    "field_type"  => "vector",
                                    "field_data"  => eachcol(self.Us))
                        ]

    str *= gt.generateVTK(filename*midpoint_suffix, midpoints; 
                                num = debug ? 0 : nothing,
                                point_data=midpoints_data, optargs...)

    # ------------- OUTPUT CONTROL POINTS --------------------------------------
    controlpoints = eachcol(self.controlpoints)

    controlpoints_data = [
                            Dict(   "field_name"  => "normal",
                                    "field_type"  => "vector",
                                    "field_data"  => eachcol(self.normals)),

                            Dict(   "field_name"  => "angleofattack",
                                    "field_type"  => "scalar",
                                    "field_data"  => self.aoas)
                        ]

    str *= gt.generateVTK(filename*controlpoint_suffix, controlpoints; 
                                point_data=controlpoints_data, optargs...)

    #  ------------- OUTPUT PLANAR GEOMETRY ------------------------------------

    str *= gt.save(self.grid, filename*planar_suffix; format, optargs...)

    return str
end
