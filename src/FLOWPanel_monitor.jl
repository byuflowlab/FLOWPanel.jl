#=##############################################################################
# DESCRIPTION
    Functions to plot and monitor specific metrics.

# AUTHORSHIP
  * Created by  : Eduardo J. Alvarez
  * Email       : Edo.AlvarezR@gmail.com
  * Date        : May 2023
  * License     : MIT License
=###############################################################################

"""
    monitor_slice(body::Union{NonLiftingBody, AbstractLiftingBody},
                        controlpoints, fieldname, sliceposs;
                        slicenormal=[0, 1, 0],
                        row=false,
                        filter=(x, i)->true,
                        proc_xfun=(points, vals)->getindex.(eachcol(points), 1),
                        proc_yfun=(points, vals)->vals,
                        # ------------- PLOT OUTPUTS ---------------------------
                        disp_plot=true,
                        _fig=nothing, _axs=nothing,
                        ttl=nothing, xscaling=1,
                        stl="-", xlims=nothing, ylims=nothing,
                        plot_optargs=[],
                        figext=".png",
                        figname="",
                        fignum=nothing,
                        # ------------- OTHER OUTPUTS --------------------------
                        save_path=nothing,
                        filepref="slice",
                        num=nothing,
                        out=[],
                        output_csv=true,
                        output_vtk=true
                        )

Plots and outputs a slice of field `fieldname`.

# ARGUMENTS
* `sliceposs::Tuple`            :   Position of each slice.
* `slicenormal::Vector`         :   Unit vector normal to the slicing plane.
* `row::Bool`                   :   Slice along first grid dimension if true;
                                    second if false.
* `filter::Function`            :   Filter points along the slice according to
                                    this logic.
* `proc_xfun::Function`         :   Process slice and output x-values for plots
* `proc_yfun::Function`         :   Process slice and output y-values for plots

NOTE: Current implementation of `find_i` does not work on MultiBody.
"""
function monitor_slice(body::Union{NonLiftingBody, AbstractLiftingBody},
                        controlpoints, fieldname, sliceposs;
                        slicenormal=[0, 1, 0],
                        row=false,
                        filter=(x, i)->true,
                        proc_xfun=(points, vals)->getindex.(eachcol(points), 1),
                        proc_yfun=(points, vals)->vals,
                        # ------------- PLOT OUTPUTS ---------------------------
                        disp_plot=true,
                        _fig=nothing, _axs=nothing,
                        ttl=nothing, xscaling=1,
                        stl="-", xlims=nothing, ylims=nothing,
                        plot_optargs=[],
                        figext=".png",
                        figname="",
                        fignum=nothing,
                        savefig=true,
                        # ------------- OTHER OUTPUTS --------------------------
                        save_path=nothing,
                        filepref="slice",
                        num=nothing,
                        out=[],
                        output_csv=true,
                        output_vtk=true
                        )

    # Create path
    if !isnothing(save_path) && !isfile(save_path)
      mkpath(save_path)
    end

    _num = num==nothing ? "" : ".$(num)"
    _fignum = fignum==nothing ? "" : ".$(fignum)"

    # Create figures
    if disp_plot
        npos = length(sliceposs)
        fig = _fig==nothing ? plt.figure(figname, figsize=[7, 5*0.75]*2/3 .* [2, ceil(npos/2)]) : _fig
        axs = _axs==nothing ? fig.subplots(ceil(Int, npos/2), 2) : _axs
        axs = _axs==nothing ? [axs[i, j] for j in 1:size(axs, 2), i in 1:size(axs, 1)] : axs
    else
        fig = _fig
        axs = _axs
    end

    for (axi, (ax, pos)) in enumerate(zip(axs, sliceposs))

        # Position of grid row/columns that are the closest to target
        points, vals = slicefield(body, controlpoints, fieldname, pos,
                                        slicenormal, row; filter=filter)
        # Process slice
        xs = proc_xfun(points, vals)
        ys = proc_yfun(points, vals)

        # Plot processed slice
        if disp_plot
            ax.plot(xs, ys, stl; plot_optargs...)

            if xlims!=nothing; ax.set_xlim(xlims); end;
            if ylims!=nothing; ax.set_ylim(ylims); end;

            ax.spines["right"].set_visible(false)
            ax.spines["top"].set_visible(false)

            # ax.set_title(ttl==nothing ? "Slice "*"$(pos)" : ttl)
        end

        # Output processed slice to array
        push!(out, [pos, xs, ys])

        if !isnothing(save_path)

            # Output processed slice to CSV
            if output_csv
                open(joinpath(save_path, filepref*"-$(axi)"*_num*".csv"), "w") do f
                    println(f, "pos,", pos)
                    for (x, y) in zip(xs, ys)
                        println(f, x, ",", y)
                    end
                end
            end

            # Output raw slice as a vtk
            if output_vtk

                nd = ndims(vals)
                field_type = nd==1 ? "scalar" : nd==2 ? "vector" : error("Invalid field with dims $(nd)")
                field_data = nd==1 ? vals : collect(eachcol(vals))

                gt.generateVTK(filepref*"-$(axi)", collect(eachcol(points));
                                    lines=[collect(0:size(points, 2))],
                                    point_data=[Dict(   "field_name"=>fieldname,
                                                        "field_type"=>field_type,
                                                        "field_data"=>vals)],
                                    num=num, path=save_path)
            end

        end

    end

    # Save plots
    if disp_plot

        fig.tight_layout()

        if !isnothing(save_path) && savefig
            fig.savefig(joinpath(save_path, filepref*_fignum*figext);
                                                    dpi=300, transparent=true)
        end

    end

    return fig, axs

end

function monitor_slice(body::AbstractBody, fieldname, sliceposs; optargs...)

    normals = calc_normals(body)
    controlpoints = calc_controlpoints(body, normals)

    return monitor_slice(body::AbstractBody, controlpoints, fieldname,
                                                        sliceposs; optargs...)
end




function monitor_loading(body::Union{NonLiftingBody, AbstractLiftingBody}, Lhat, Dhat, b;
                        spandirection=[0, 1, 0], dimspan=2, dimchord=1,
                        to_plot=[1, 2, 3], yscalings=[1.0, 1.0, 1.0],
                        F_fieldname="F",
                        # ------------- PLOT OUTPUTS ---------------------------
                        disp_plot=true,
                        _fig=nothing, _axs=nothing,
                        stl="-", xlims=[-1, 1], ylims=nothing,
                        xlbl=L"Span position $2y/b$",
                        ylbls="Sectional ".*[L"lift $\ell$", L"drag $d$", L"sideslip $s$"].*" (N/m)",
                        plot_optargs=[],
                        figext=".png",
                        figname="",
                        fignum=nothing,
                        savefig=true,
                        # ------------- OTHER OUTPUTS --------------------------
                        save_path=nothing,
                        filepref="loading",
                        num=nothing,
                        out=[],
                        output_csv=true
                        )

    # Calculate sectional force along the span
    fs, spanpos = calcfield_sectionalforce(body; spandirection=spandirection,
                                                    dimspan=dimspan,
                                                    dimchord=dimchord,
                                                    F_fieldname=F_fieldname)

    # Decompose the force into lift, drag, and sideslip
    lds = decompose(fs, Lhat, Dhat)

    # Create path
    if !isnothing(save_path) && !isfile(save_path)
      mkpath(save_path)
    end

    _num = num==nothing ? "" : ".$(num)"
    _fignum = fignum==nothing ? "" : ".$(fignum)"

    # Plot loading
    if disp_plot
        fig = _fig==nothing ? plt.figure(figname, figsize=[7, 5*0.75]*2/3 .* [length(to_plot), 1]) : _fig
        axs = _axs==nothing ? fig.subplots(1, length(to_plot)) : _axs
        axs = _axs==nothing ? length(to_plot)==1 ? [axs] : [axs[i, j] for j in 1:size(axs, 2), i in 1:size(axs, 1)] : axs
    else
        fig = _fig
        axs = _axs
    end

    xs = spanpos*2/b
    yss = []

    for (axi, (ax, pi)) in enumerate(zip(axs, to_plot))

        ys = lds[pi, :]*yscalings[pi]
        push!(yss, ys)

        # Plot processed data
        if disp_plot
            ax.plot(xs, ys, stl; plot_optargs...)

            if xlims!=nothing; ax.set_xlim(xlims); end;
            if ylims!=nothing; ax.set_ylim(ylims[axi]); end;

            ax.set_xlabel(xlbl)
            ax.set_ylabel(ylbls[axi])

            ax.spines["right"].set_visible(false)
            ax.spines["top"].set_visible(false)
        end

    end

    # Push processed data to array
    push!(out, [xs, yss])

    # Output processed data to CSV
    if !isnothing(save_path) && output_csv
        open(joinpath(save_path, filepref*_num*".csv"), "w") do f
            println(f, "posscaled"*prod([",$(var)scaled" for var in ["l", "d", "s"][to_plot]]))
            for xy in zip(xs, yss...)
                print(f, xy[1])
                for y in xy[2:end]
                    print(f, ",", y)
                end
                println(f)
            end
        end
    end

    # Save plots
    if disp_plot

        # fig.tight_layout()

        if !isnothing(save_path) && savefig
            fig.savefig(joinpath(save_path, filepref*_fignum*figext);
                                                    dpi=300, transparent=true)
        end

    end

    return fig, axs

end


function monitor_loading(multibody::MultiBody, args...;
                            _fig=nothing, _axs=nothing,
                            filepref="loading", optargs...)

    fig, axs = _fig, _axs

    for (bi, body) in enumerate(multibody.bodies)

        fig, axs = monitor_loading(body, args...; _fig=fig, _axs=axs,
                                    filepref=filepref*"-b$(bi)", optargs...)
    end

    return fig, axs

end
