#=##############################################################################
# DESCRIPTION
    Linear solver for lifting line method (aka Weissinger VLM).

# AUTHORSHIP
  * Created by  : Eduardo J. Alvarez
  * Email       : Edo.AlvarezR@gmail.com
  * Date        : Nov 2025
  * License     : MIT License
=###############################################################################

function solve_linear(self::LiftingLine, Uinf::AbstractVector, 
                                            args...; optargs...) 
    solve_linear(self, repeat(Uinf, 1, self.nelements), args...; optargs...)
end

function solve_linear(self::LiftingLine, Uinfs::AbstractMatrix;
                        solver=solve_ludiv!, solver_optargs=(),
                        addfields=true, raise_warn=false,
                        optargs...)

    # Update semi-infinite wake to align with freestream
    calc_Dinfs!(self, Uinfs)

    # Compute geometric matrix (left-hand-side influence matrix)
    self.G .= 0
    _G_U!(self; optargs...)

    # Calculate boundary conditions (right-hand side of system of equations)
    self.RHS .= 0
    calc_bc_noflowthrough!(self.RHS, Uinfs, self.normals)

    # Solve system of equations
    solver(self.Gammas, self.G, self.RHS; solver_optargs...)

    # Calculate velocity at lifting-line midpoints
    self.Us .= Uinfs
    Uind!(self, self.midpoints, self.Us)

    if addfields
        gt.add_field(self.grid, "Uinf", "vector", collect(eachcol(Uinfs)), "cell"; raise_warn)
        gt.add_field(self.grid, "Gamma", "scalar", self.Gammas, "cell"; raise_warn)
        gt.add_field(self.grid, "angleofattack", "scalar", self.aoas, "cell"; raise_warn)
    end

    return nothing

end

function _G_U!(self::LiftingLine; 
                    ground_normal=self.ground_normal, 
                    ground_position=self.ground_position,
                    optargs...)

    # Ground plane arguments
    if isfinite(norm(ground_position))
        ground = (ground_position..., ground_normal...)
    else
        ground = ()
    end

    # Build geometric matrix from panel contributions
    elements = 1:self.nelements
    chunks = collect(Iterators.partition(elements, max(length(elements) ÷ Threads.nthreads(), 3*Threads.nthreads())))

    Threads.@threads for chunk in chunks      # Distribute panel iteration among all CPU threads

        for ei in chunk                       # Iterate over elements

            U_vortexring(
                              view(self.horseshoes, :, :, ei),   # All nodes
                              1:4,                               # Indices of nodes that make this panel (closed ring)
                              1.0,                               # Unitary strength
                              ground...,
                              self.controlpoints,                # Targets
                              view(self.G, :, ei);               # Velocity of ei-th panel on every CP
                              dot_with=self.normals,             # Normal of every CP
                              offset=self.kerneloffset,          # Offset of kernel to avoid singularities
                              cutoff=self.kernelcutoff,          # Kernel cutoff to avoid singularities
                              optargs...
                             )

         end
    end


    # Add wake contributions
    TE = [1, size(self.horseshoes, 2)]          # Indices of TE nodes in each horseshoe

    Threads.@threads for chunk in chunks        # Distribute wake iteration among all CPU threads

        for ei in chunk                          # Iterate horseshoes

            da1 = self.Dinfs[1, 1, ei]
            da2 = self.Dinfs[2, 1, ei]
            da3 = self.Dinfs[3, 1, ei]
            db1 = self.Dinfs[1, 2, ei]
            db2 = self.Dinfs[2, 2, ei]
            db3 = self.Dinfs[3, 2, ei]

            U_semiinfinite_horseshoe(
                              view(self.horseshoes, :, :, ei),   # All nodes
                              TE,                                # Indices of nodes that make the shedding edge
                              da1, da2, da3,                     # Semi-infinite direction da
                              db1, db2, db3,                     # Semi-infinite direction db
                              1.0,                               # Unitary strength
                              ground...,
                              self.controlpoints,                # Targets
                              view(self.G, :, ei);               # Velocity of wake panel on every CP
                              dot_with=self.normals,             # Normal of every CP
                              offset=self.kerneloffset,          # Offset of kernel to avoid singularities
                              cutoff=self.kernelcutoff,          # Kernel cutoff to avoid singularities
                              optargs...
                             )

        end

    end
end