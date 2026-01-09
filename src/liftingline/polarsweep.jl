#=##############################################################################
# DESCRIPTION
    Method for running sweeps generating CL, CD, and CM vs AOA lookup tables

# AUTHORSHIP
  * Created by  : Eduardo J. Alvarez
  * Email       : Edo.AlvarezR@gmail.com
  * Date        : Jan 2026
  * License     : MIT License
=###############################################################################

function run_polarsweep(liftingline::LiftingLine, 
                        magUinf::Number,                    # (m/s) freestream velocity
                        rho::Number,                        # (kg/m^3) air density
                        X0::AbstractVector,                 # (m) center about which to calculate moments
                        cref::Number,                       # (m) reference chord
                        bref::Number;                       # (m) reference span
                        Aref = cref*bref,                   # (m^2) reference area

                    # ------- SWEEP PARAMETERS ---------------------------------
                    aoa_sweeps = [range(0, -50, step=-0.5), range(0.5, 50, step=0.5)],  # Set of AOA sweeps to run
                    fixed_parameters = Dict(),                                          # Other parameters that are fixed during this sweep

                    # ------- SOLVER SETTINGS ----------------------------------
                    solver = SimpleNonlinearSolve.SimpleDFSane(),                       # Indifferent to initial guess, but somewhat not robust
                    # solver = SimpleNonlinearSolve.SimpleTrustRegion()                 # Trust region needs a good initial guess, but it converges very reliably
                    solver_optargs  = (; abstol = 1e-9, maxiters = 200),
                    align_joints_with_Uinfs = false,                                    # Whether to align joint bound vortices with the freestream
                    reattempt_align_method = false,                                     # Re-attempt joint alingment if first option fails (`true` for more robust polars, but might take longer)
                    use_Uind_for_force = false,                                         # Whether to use Uind as opposed to selfUind for force postprocessing
        
                    # ------- OUTPUTS ------------------------------------------
                    output_distributions = nothing,                                     # Give it an array and it will output spanwise distributions here
                    # File outputs
                    sweepname = "polarsweep",                                           # Name of this sweep
                    masterdatabase_file = "alldata.csv",                                # File where to append this sweep

                    save_path = nothing,                                                # If provided, it saves all outputs under this path
                    sweep_path = nothing,                                               # If provided, it stores outputs under subfolders

                    polars_path=!isnothing(sweep_path) ? joinpath(sweep_path, "polars") : save_path,        # Where to save polars
                    plots_path=!isnothing(sweep_path) ? joinpath(sweep_path, "plots") : save_path,     # Where to save plots
                    extraplots_path=!isnothing(sweep_path) ? joinpath(sweep_path, "plots-extra") : save_path,    # Where to save extra plots

                    # Plotting
                    plot_sweep=true,        # Plots to generate
                    lbl = "FLOWPanel LL",
                    stl = ".-",
                    fmt = (; color="steelblue", markersize=3, alpha=0.6),
        
                    # Interactions with user
                    prompt = true,                                                      # Whether to prompt the user
                    verbose=true, tab_lvl=0,                                            # Settings for verbosing
                    verbose_sweep0=true, verbose_sweep=true

    )

    nondim = 0.5*rho*magUinf^2*Aref      # Normalization factor

    ll = liftingline

    ############################################################################
    #                          INITIALIZATION
    ############################################################################

    header = vcat("aoa (deg)", keys(fixed_parameters), "CD", "CY", "CL", "Cl", "Cm", "Cn")
    header = join(header, ",")

    # Initialize sweep path
    if !isnothing(sweep_path) && !ispath(sweep_path)
        
        # Create sweep path
        mkdir(sweep_path)
        mkdir(plots_path)
        mkdir(plots_path)
        mkdir(polars_path)

        # Initialize master database file
        open(joinpath(sweep_path, masterdatabase_file), "w") do f
            println(f, header)
        end
        
    end


    
    ############################################################################
    #                          RUN AOA SWEEP
    ############################################################################
    solver_cache = Dict()

    _alphas, _CDs, _CYs, _CLs, _Cls, _Cms, _Cns = ([] for i in 1:7)
    _cds, _cys, _cls = ([] for i in 1:3)
    succeeded = []
    
    if verbose_sweep0
        println("\t"^(tab_lvl) * "Running AOA sweep...")
    end
        
    for AOAs in aoa_sweeps
    
        success_previous = false
        last_successfull_aoas = nothing
        last_successfull = 0
        
        for AOA in AOAs

            if verbose_sweep
                println("\t"^(tab_lvl+1)*"\n##################################")
                @show AOA
            end
            
            Uinf  = magUinf*[cosd(AOA), 0, sind(AOA)] # Freestream
                
            Uinfs = repeat(Uinf, 1, liftingline.nelements)
        
            
            # ------------------ CALL NONLINEAR SOLVER ------------------------
            success = false
            for (tryi, attempt) in enumerate([align_joints_with_Uinfs, !align_joints_with_Uinfs][1:2^reattempt_align_method])
        
                if !attempt
                    jointerize!(ll; deltajoint=ll.deltajoint)
                end
        
                if tryi==2
                    println("\nalign_joints_with_Uinfs=$(!attempt) case did not converge. Trying $(attempt) now.\n" )
                end
                
                aoas_initial_guess = !isnothing(last_successfull_aoas) && last_successfull <= 3 ? last_successfull_aoas : AOA
                
                result, solver_cache = solve(liftingline, Uinfs; debug=!true, 
                                                    aoas_initial_guess, align_joints_with_Uinfs=attempt, 
                                                    solver, solver_optargs, solver_cache)
        
                # Check solver success
                success = SimpleNonlinearSolve.SciMLBase.successful_retcode(result)
                
                if verbose_sweep
                    @show solver_cache[:fcalls]
                    # display(result)
                    @show success
                end
                
                if success
                    break
                end
            end
            
            
            # ------------------ POSTPROCESSING -------------------------------
            
            Dhat            = Uinf/norm(Uinf)               # Drag direction
            Shat            = [0, 1, 0]                     # Span direction
            Lhat            = cross(Dhat, Shat)             # Lift direction
            
            lhat            = Dhat                          # Rolling direction
            mhat            = Shat                          # Pitching direction
            nhat            = Lhat                          # Yawing direction
            
            # NOTE: Coefficients must be evaluated on using the velocity from 
            #       the effective horseshoes as shown below, which is automatically
            #       computed by the solver already
            # ll.Us .= Uinfs
            # selfUind!(ll)
            
            # Calculate stripwise coefficients
            calcfield_cl(liftingline)
            calcfield_cd(liftingline)
            calcfield_cm(liftingline)
            
            # Convert velocity to effective swept velocity
            # NOTE: Forces must use the velocity from the original horseshoes for
            #       best accuracy, as done here
            if use_Uind_for_force
                ll.Us .= Uinfs
                Uind!(ll, ll.midpoints, ll.Us)
            end
            calc_UΛs!(ll, ll.Us)
            
            # Force per stripwise element integrating lift and drag coefficient
            calcfield_F(liftingline, rho)
            
            # Integrated force
            Ftot = calcfield_Ftot(liftingline)
            
            # Integrated force decomposed into lift and drag
            LDS = calcfield_LDS(liftingline, Lhat, Dhat, Shat)
            
            L = LDS[:, 1]
            D = LDS[:, 2]
            S = LDS[:, 3]
            
            if !isnothing(output_distributions)
                # Loading distribution (force per unit span)
                fs = calcfield_f(liftingline)
                
                lds = decompose(fs, Lhat, Dhat)
                
                l = lds[1, :]
                d = lds[2, :]
                s = lds[3, :]
            
                cl = l / (nondim/bref)
                cd = d / (nondim/bref)
                cy = s / (nondim/bref)
            end
    
            # Integrated moment
            Mtot = calcfield_Mtot(ll, X0, rho)
            
            # Moment decomposed into axes
            lmn = calcfield_lmn(ll, lhat, mhat, nhat)
            roll, pitch, yaw = collect(eachcol(lmn))
            
            # Coefficients
            CL = sign(dot(L, Lhat)) * norm(L) / nondim
            CD = sign(dot(D, Dhat)) * norm(D) / nondim
            CY = sign(dot(S, Shat)) * norm(S) / nondim
            
            Cl = sign(dot(roll, lhat)) * norm(roll) / (nondim*cref)
            Cm = sign(dot(pitch, mhat)) * norm(pitch) / (nondim*cref)
            Cn = sign(dot(yaw, nhat)) * norm(yaw) / (nondim*cref)

            if verbose_sweep
                @show CL
                @show CD
                @show Cm
            end
    
            # Detect blown up case
            if success && abs(CL) >= 500
                success = false
            end
    
            push!(succeeded, success)
            push!(_alphas, AOA)
            push!(_CDs, success ? CD : NaN)
            push!(_CYs, success ? CY : NaN)
            push!(_CLs, success ? CL : NaN)
            push!(_Cls, success ? Cl : NaN)
            push!(_Cms, success ? Cm : NaN)
            push!(_Cns, success ? Cn : NaN)

            if !isnothing(output_distributions)
                push!(_cds, success ? cd : NaN)
                push!(_cys, success ? cy : NaN)
                push!(_cls, success ? cl : NaN)
            end
        
            if success
                last_successfull_aoas = deepcopy(ll.aoas)
                last_successfull = 0
            else
                last_successfull += 1
            end
            
            success_previous = success
        
        end
    
    end

    # Sort results by AOA while filtering unconverged points
    sorted_indices = sort(collect(enumerate(_alphas)), by= x -> x[2])
    sorted_indices = [i for (i, a) in sorted_indices if succeeded[i]]
    
    AOAs = _alphas[sorted_indices]
    CDs = _CDs[sorted_indices]
    CYs = _CYs[sorted_indices]
    CLs = _CLs[sorted_indices]
    Cls = _Cls[sorted_indices]
    Cms = _Cms[sorted_indices]
    Cns = _Cns[sorted_indices]

    if !isnothing(output_distributions)
        cds = _cds[sorted_indices]
        cys = _cys[sorted_indices]
        cls = _cls[sorted_indices]
    end

    success_rate = sum(succeeded) / length(succeeded)
    
    if verbose_sweep0
        println("\t"^(tab_lvl) * "Success rate: $(Int(round(100*success_rate)))%")
    end

    # Output spanwise distributions
    if !isnothing(output_distributions)

        ypos = (ll.ypositions[2:end] .+ ll.ypositions[1:end-1]) / 2
        spanwise_distributions = (; AOA=AOAs, yposition=ypos, cd=cds, cy=cys, cl=cls)

        push!(output_distributions, spanwise_distributions)
    end


    files_to_write = []

    # Write sweep to individual polar file
    if !isnothing(polars_path)
        push!(files_to_write, (joinpath(polars_path, sweepname*".csv"), "w"))
    end

    # Write sweep to master database file
    if !isnothing(sweep_path)
        push!(files_to_write, (joinpath(sweep_path, masterdatabase_file), "a"))
    end
    
    # Write sweep to files
    for (file, flag) in files_to_write
        open(file, flag) do f

            if flag=="w"
                println(f, header)
            end
            
            for (aoa, CD, CY, CL, Cl, Cm, Cn) in zip(AOAs, CDs, CYs, CLs, Cls, Cms, Cns)

                entry = vcat(aoa, values(fixed_parameters), CD, CY, CL, Cl, Cm, Cn)
                entry = join(entry, ",")
                println(f, entry)

            end

        end
    end
    
    
    ##############################################################################
    #                          PLOT SWEEP
    ##############################################################################
    if plot_sweep
        
        # Compare database and getters
        fig = plt.figure(figsize = [7*2, 0.75*5*2]*2/3)
        axs = fig.subplots(2, 2)

        axs = [axs[j, i] for i in 1:size(axs, 1), j in 1:size(axs, 2)]
        
        ax = axs[1]
        ax.plot(AOAs, CLs, stl; label=lbl, fmt...)
        
        ax.set_ylabel(L"Lift coeff. $C_L$")
        # ax.set_ylim([-0.8, 1.0])
        # ax.set_yticks(range(-0.5, 1.1, step=0.5))
        
        ax = axs[2]
        ax.plot(AOAs, CDs, stl; label=lbl, fmt...)
        
        ax.set_ylabel(L"Drag coeff. $C_D$")
        # ax.set_ylim([-0.005, 0.075])
        # ax.set_yticks(range(0, 0.08, step=0.025))
        
        ax = axs[3]
        ax.plot(CDs, CLs, stl; label=lbl, fmt...)
        
        ax.set_xlabel(L"Drag coeff. $C_D$")
        # ax.set_xlim([-0.005, 0.075])
        # ax.set_xticks(range(0, 0.08, step=0.025))
        
        ax.set_ylabel(L"Lift coeff. $C_L$")
        # ax.set_ylim([-0.8, 1.0])
        # ax.set_yticks(range(-0.5, 1.1, step=0.5))
        
        ax = axs[4]
        ax.plot(AOAs, Cms, stl; label=lbl, fmt...)
        
        ax.set_ylabel(L"Pitching moment $C_m$")
        # ax.set_ylim([-1.2, 1.0])
        # ax.set_yticks(range(-1.0, 1.1, step=0.5))
        
        for (axi, ax) in enumerate(axs)
            
            if axi != 3
                ax.set_xlabel(L"Angle of attack ($^\circ$)")
                # ax.set_xlim([-18, 18])
                # ax.set_xlim([-30, 40])
                # ax.set_xticks(range(-30, 40, step=10))
            end
                
            ax.spines["top"].set_visible(false)
            ax.spines["right"].set_visible(false)
            
            ax.legend(loc="best", frameon=false, fontsize=8)
            
        end
            
        fig.tight_layout()

        if !isnothing(plots_path)
            for ext in [".png", ".pdf"]
                fig.savefig(joinpath(plots_path, sweepname*"-sweep"*ext), transparent=true, dpi=300)
            end
        end
            
        ax = axs[1]
        ax.set_yticks(range(-1.0, 1.0, step=0.5))
        ax.set_ylim([-1.2, 1.2])
        
        ax = axs[2]
        ax.set_yticks(range(0, 0.1, step=0.02))
        ax.set_ylim([0, 0.1])
        
        ax = axs[3]
        ax.set_xticks(range(0, 0.1, step=0.02))
        ax.set_xlim([0, 0.1])
        
        ax.set_yticks(range(-1.0, 1.0, step=0.5))
        ax.set_ylim([-1.2, 1.2])
        
        ax = axs[4]
        ax.set_yticks(range(-1.6, 1.6, step=0.5))
        ax.set_ylim([-1.6, 1.6])
        
        for (axi, ax) in enumerate(axs)
            
            if axi != 3
                ax.set_xticks(range(-15, 15, step=5))
                ax.set_xlim([-16, 16])
            end
            
        end
            
        if !isnothing(extraplots_path)
            for ext in [".png", ".pdf"]
                fig.savefig(joinpath(extraplots_path, sweepname*"-sweep2"*ext), transparent=true, dpi=300)
            end
        end

    end

    return (; AOA=AOAs, CD=CDs, CY=CYs, CL=CLs, Cl=Cls, Cm=Cms, Cn=Cns)

end