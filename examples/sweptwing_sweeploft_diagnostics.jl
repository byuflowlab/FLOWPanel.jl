#=##############################################################################
# DESCRIPTION
#   Post-convergence diagnostics for the sweep-loft swept-wing study.
#   Loads the final-level VTK state written by examples/sweptwing_sweeploft.jl
#   (NO re-solving) for both mirror-symmetric diagonal variants (diagAC/diagBD)
#   and investigates why they converge to different CL:
#
#     1. Quad-level comparison of gamma / velocity / Cp between the two
#        triangulations of the IDENTICAL node set (the two variants share the
#        same quads; cells 2q-1, 2q are the two triangles of quad q in both).
#     2. Spanwise distribution of the lift difference (which spanwise stations
#        drive the CL gap).
#     3. Lower-TE adjacent-jump diagnostics (same as examples/sweptwing.jl).
#     4. Diagonal-edge orientation relative to the freestream per variant.
#
#   Writes findings to data/sweptwing_sweeploft/findings_data.md
=###############################################################################

ENV["FLOWPANEL_SWEEPLOFT_RUN"] = "false"
include(joinpath(@__DIR__, "sweptwing_sweeploft.jl"))

import Statistics: mean, median

# ----------------- LOAD SAVED STATE -------------------------------------------
function load_vtu_cell_field(path, body_name, field_name; idx=0)
    vtu_path = joinpath(path, body_name, "$(body_name).$(idx).vtu")
    isfile(vtu_path) || error("Body VTU not found: $(vtu_path)")
    vtk = pnl.ReadVTK.VTKFile(vtu_path)
    cell_data = pnl.ReadVTK.get_cell_data(vtk)
    field_name in keys(cell_data) ||
        error("Field '$(field_name)' missing from $(vtu_path); available: $(collect(keys(cell_data)))")
    return Array{Float64}(pnl.ReadVTK.get_data(cell_data[field_name]))
end

"""
Rebuild the final mesh of a variant and populate it with the saved solution
fields (gamma, velocity). Returns `(body, pressure)` where `pressure` is the
saved Bernoulli gauge pressure per cell.
"""
function load_variant(n_ch, n_span, variant_name, swap_diagonals)
    body = sweeploft_wing(n_ch, n_span; mirror_diagonals=true, swap_diagonals)
    body_name = run_name * "_" * variant_name * "_body1"

    gamma = vec(load_vtu_cell_field(save_path, body_name, "gamma"))
    velocity = load_vtu_cell_field(save_path, body_name, "velocity")
    pressure = vec(load_vtu_cell_field(save_path, body_name, "gauge pressure"))
    F = load_vtu_cell_field(save_path, body_name, "F")

    length(gamma) == body.ncells ||
        error("$(variant_name): loaded gamma has $(length(gamma)) cells, mesh has $(body.ncells)")
    size(velocity, 2) == body.ncells || (velocity = collect(velocity'))

    body.strength[:, pnl.get_Gammai(body)] .= gamma
    body.velocity .= velocity
    pnl.calc_normals!(body)
    pnl.calc_controlpoints!(body)
    return body, pressure, F
end

# ----------------- TE DIAGNOSTICS (from examples/sweptwing.jl) ----------------
function reconstructed_half_jump(body)
    jump = zeros(3, body.ncells)
    pnl.compute_mu_gradient!(jump, body.controlpoints, body.normals,
        body.cells, body.neighbor,
        view(body.strength, :, pnl.get_Gammai(body)),
        view(body.shedding_full, 1:2, :);
        scale=0.5,
        nodes=body.nodes)
    return jump
end

function lower_te_panels(body)
    panels = Int[]
    for shedding in body.shedding
        for col in eachcol(shedding)
            pi, pj = col[1], col[4]
            if pj == -1
                push!(panels, pi)
            else
                lower = body.controlpoints[3, pi] <= body.controlpoints[3, pj] ? pi : pj
                push!(panels, lower)
            end
        end
    end
    panels = unique(panels)
    sort!(panels; by=p -> body.controlpoints[2, p])
    return panels
end

function adjacent_jump_summary(values, panels)
    length(panels) >= 2 || return (; n=0, mean=NaN, rms=NaN, max=NaN)
    jumps = Float64[]
    size(values, 1) == 1 && (values = reshape(vec(values), 1, :))
    for k in 1:(length(panels)-1)
        i, j = panels[k], panels[k+1]
        push!(jumps, LA.norm(values[:, j] .- values[:, i]))
    end
    return (; n=length(jumps),
             mean=sum(jumps)/length(jumps),
             rms=sqrt(sum(abs2, jumps)/length(jumps)),
             max=maximum(jumps))
end

function te_diagnostics(io, label, body, pressure)
    jump = reconstructed_half_jump(body)
    panels = lower_te_panels(body)
    gamma = reshape(collect(view(body.strength, :, pnl.get_Gammai(body))), 1, :)
    pv_induced = (body.velocity .- jump) .- repeat(Vinf, 1, body.ncells)
    println(io, "\n### Lower-TE adjacent jumps: $(label)")
    for (name, vals) in (("gamma", gamma),
                         ("PV induced velocity", pv_induced),
                         ("half-jump velocity", jump),
                         ("final velocity", body.velocity),
                         ("Bernoulli pressure", reshape(collect(pressure), 1, :)))
        s = adjacent_jump_summary(vals, panels)
        println(io, "- ", rpad(name, 22), " n=$(s.n)",
            " mean=$(round(s.mean, sigdigits=5))",
            " rms=$(round(s.rms, sigdigits=5))",
            " max=$(round(s.max, sigdigits=5))")
    end
end

# ----------------- QUAD-LEVEL VARIANT COMPARISON ------------------------------
"Average a per-cell scalar onto quads (cells 2q-1, 2q are quad q's triangles)."
quad_average(vals::AbstractVector) = (vals[1:2:end] .+ vals[2:2:end])/2
quad_average(vals::AbstractMatrix) = (vals[:, 1:2:end] .+ vals[:, 2:2:end])/2

"Quad centroid y from the body's cells (works for either variant; same quads)."
function quad_y(body)
    yq = zeros(body.ncells ÷ 2)
    for q in eachindex(yq)
        nis = unique(vcat(body.cells[:, 2q-1], body.cells[:, 2q]))
        yq[q] = mean(body.nodes[2, ni] for ni in nis)
    end
    return yq
end

function quad_comparison(io, bodyAC, pAC, bodyBD, pBD)
    gAC = quad_average(vec(collect(view(bodyAC.strength, :, pnl.get_Gammai(bodyAC)))))
    gBD = quad_average(vec(collect(view(bodyBD.strength, :, pnl.get_Gammai(bodyBD)))))
    uAC = quad_average(bodyAC.velocity)
    uBD = quad_average(bodyBD.velocity)
    cpAC = quad_average(pAC) ./ (0.5*rho*magVinf^2)
    cpBD = quad_average(pBD) ./ (0.5*rho*magVinf^2)

    dg = abs.(gAC .- gBD)
    du = vec(sqrt.(sum(abs2, uAC .- uBD; dims=1)))
    dcp = abs.(cpAC .- cpBD)

    println(io, "\n### Quad-level diagAC vs diagBD (same node set, same quads)")
    println(io, "- |Δgamma|/max|gamma|: mean=", round(mean(dg)/maximum(abs.(gAC)), sigdigits=4),
        " max=", round(maximum(dg)/maximum(abs.(gAC)), sigdigits=4))
    println(io, "- |Δvelocity|/Vinf:    mean=", round(mean(du)/magVinf, sigdigits=4),
        " max=", round(maximum(du)/magVinf, sigdigits=4))
    println(io, "- |ΔCp|:               mean=", round(mean(dcp), sigdigits=4),
        " max=", round(maximum(dcp), sigdigits=4))
    return (; dg, du, dcp)
end

# ----------------- SPANWISE LIFT-DIFFERENCE DISTRIBUTION ----------------------
function spanwise_lift(body, F, Lhat, nbins)
    yq = view(body.controlpoints, 2, :)
    edges = range(-semispan, semispan; length=nbins + 1)
    lift = zeros(nbins)
    for ci in 1:body.ncells
        bin = clamp(searchsortedlast(edges, yq[ci]), 1, nbins)
        lift[bin] += LA.dot(view(F, :, ci), Lhat)
    end
    centers = (edges[1:end-1] .+ edges[2:end])/2
    return collect(centers), lift
end

# ----------------- DIAGONAL ORIENTATION ----------------------------------------
"""
Mean |cos(angle)| between each quad's diagonal edge and the freestream
direction, mid-span quads only (|y| < semispan - L_t with L_t = c).
"""
function diagonal_alignment(body, n_ch)
    Dhat = Vinf/LA.norm(Vinf)
    npts = 2*n_ch
    aligns = Float64[]
    for q in 1:(body.ncells ÷ 2)
        tri = view(body.cells, :, 2q-1)
        # diagonal = edge shared by the two triangles of the quad
        tri2 = view(body.cells, :, 2q)
        shared = intersect(tri, tri2)
        length(shared) == 2 || continue
        e = body.nodes[:, shared[2]] .- body.nodes[:, shared[1]]
        ymid = (body.nodes[2, shared[1]] + body.nodes[2, shared[2]])/2
        abs(ymid) < semispan - c || continue
        push!(aligns, abs(LA.dot(e, Dhat))/LA.norm(e))
    end
    return (; mean=mean(aligns), median=median(aligns))
end

# ----------------- MAIN --------------------------------------------------------
import DataFrames
conv = CSV.read(joinpath(save_path, "convergence.csv"), DataFrames.DataFrame)
finals = Dict(v => last(conv[conv.variant .== v, :]) for v in ("diagAC", "diagBD"))
n_ch_f, n_span_f = finals["diagAC"].n_ch, finals["diagAC"].n_span
println("Final level: $(n_ch_f)x$(n_span_f) ($(finals["diagAC"].panels) panels)")

bodyAC, pAC, FAC = load_variant(n_ch_f, n_span_f, "diagAC", false)
bodyBD, pBD, FBD = load_variant(n_ch_f, n_span_f, "diagBD", true)

Dhat = Vinf/LA.norm(Vinf)
Lhat = LA.cross(Dhat, [0.0, 1.0, 0.0])

open(joinpath(save_path, "findings_data.md"), "w") do io
    println(io, "# Sweep-loft diagnostics (final level $(n_ch_f)x$(n_span_f), ",
        finals["diagAC"].panels, " panels)")
    println(io, "\nCL diagAC = ", finals["diagAC"].CL,
        "   CL diagBD = ", finals["diagBD"].CL,
        "   ΔCL = ", finals["diagAC"].CL - finals["diagBD"].CL)

    aAC = diagonal_alignment(bodyAC, n_ch_f)
    aBD = diagonal_alignment(bodyBD, n_ch_f)
    println(io, "\n### Mid-span diagonal alignment with freestream |cos|")
    println(io, "- diagAC: mean=", round(aAC.mean, sigdigits=4),
        " median=", round(aAC.median, sigdigits=4))
    println(io, "- diagBD: mean=", round(aBD.mean, sigdigits=4),
        " median=", round(aBD.median, sigdigits=4))

    quad_comparison(io, bodyAC, pAC, bodyBD, pBD)

    nbins = min(n_span_f, 48)
    yc, liftAC = spanwise_lift(bodyAC, FAC, Lhat, nbins)
    _, liftBD = spanwise_lift(bodyBD, FBD, Lhat, nbins)
    dlift = liftAC .- liftBD
    println(io, "\n### Spanwise lift difference diagAC - diagBD (", nbins, " bins)")
    println(io, "- total ΔL = ", round(sum(dlift), sigdigits=5),
        " (", round(100*sum(dlift)/sum(liftBD), sigdigits=3), "% of diagBD lift)")
    imax = argmax(abs.(dlift))
    println(io, "- largest bin Δ at 2y/b = ", round(2yc[imax]/b, sigdigits=3),
        ": ", round(dlift[imax], sigdigits=4))
    println(io, "- 2y/b vs Δlift:")
    for i in 1:nbins
        println(io, "  ", round(2yc[i]/b, digits=3), ", ", round(dlift[i], sigdigits=5))
    end
    CSV.write(joinpath(save_path, "spanwise_lift_difference.csv"),
        DataFrames.DataFrame(y2b=2 .* yc ./ b, liftAC=liftAC, liftBD=liftBD,
                             dlift=dlift))

    te_diagnostics(io, "diagAC", bodyAC, pAC)
    te_diagnostics(io, "diagBD", bodyBD, pBD)
end

println(read(joinpath(save_path, "findings_data.md"), String))

# ----------------- QUAD-BASED ∇μ EXPERIMENT -----------------------------------
# Rebuild the surface velocity manually from the loaded γ (no re-solve):
#   U = U∞ + U_PV (influence! FMM pass) + ½∇μ
# with ∇μ from (a) the existing triangle scheme (validation against the VTU),
# (b) the structured-quad LS gradient, (c) the node-based P1 gradient.

"""
Recompute `U∞ + U_PV` into body.velocity via the same influence! pass solve!
uses post-solve (src/FLOWPanel_simulate.jl:409-426), and return it along with
the triangle-based half-jump field Jtri (compute_mu_gradient!, scale=0.5).
body.strength must already hold the loaded γ.
"""
function manual_velocity!(body)
    Gi = pnl.get_Gammai(body)
    gamma = copy(body.strength[:, Gi])
    pnl.reset!(body)                      # zeroes velocity/potential, not strength
    body.strength[:, Gi] .= gamma
    pnl.apply_freestream!(body, Vinf)
    pnl._set_kerneloffsets!((body,), :kerneloffset_targets)
    pnl.influence!((body,), (body,), pnl.FastMultipoleBackend();
        precalc=false, scalar_potential=false, velocity=true,
        velocity_gradient=(false,),
        direct_conditioning=pnl._self_panel_kerneloffset_conditioning())
    U0 = copy(body.velocity)              # U∞ + PV

    Jtri = zeros(3, body.ncells)
    pnl.compute_mu_gradient!(Jtri, body.controlpoints, body.normals,
        body.cells, body.neighbor, view(body.strength, :, Gi),
        view(body.shedding_full, 1:2, :); scale=0.5, nodes=body.nodes)
    return U0, Jtri
end

"Quad-average a 3xncells triangle field with triangle-area weights."
function quad_field(vals, areas)
    Nq = size(vals, 2) ÷ 2
    out = zeros(3, Nq)
    for q in 1:Nq
        a1, a2 = areas[2q-1], areas[2q]
        out[:, q] .= (a1 .* view(vals, :, 2q-1) .+ a2 .* view(vals, :, 2q)) ./ (a1 + a2)
    end
    return out
end

"""
Calibrate the half-jump sign: pick s ∈ ±1 such that s*0.5*gradmu best matches
the triangle half-jump on smooth mid-span quads (excluding TE rows and tips).
"""
function calibrate_jump_sign(geom, gradmu, Jtri_quad, n_ch, n_span)
    npts = 2*n_ch
    acc = 0.0
    for k in 3:(n_span - 2), j in 3:(npts - 2)
        q = (k - 1)*npts + j
        abs(geom.centroid[2, q]) < semispan - c || continue
        acc += LA.dot(0.5 .* view(gradmu, :, q), view(Jtri_quad, :, q))
    end
    return acc >= 0 ? 1.0 : -1.0
end

function reconstruction_comparison()
    n_ch, n_span = n_ch_f, n_span_f
    println("\n#===== QUAD-BASED ∇μ RECONSTRUCTION ($(n_ch)x$(n_span)) =====#")

    Dhat = Vinf/LA.norm(Vinf)
    Lhat = LA.cross(Dhat, [0.0, 1.0, 0.0])

    out = Dict{String, Any}()
    for (name, body) in (("diagAC", bodyAC), ("diagBD", bodyBD))
        U0, Jtri = manual_velocity!(body)
        geom = quad_geometry(body)

        # Validation: manual tri velocity must reproduce the VTU velocity
        Utri = U0 .+ Jtri
        vtu_vel = name == "diagAC" ?
            load_vtu_cell_field(save_path, run_name*"_diagAC_body1", "velocity") :
            load_vtu_cell_field(save_path, run_name*"_diagBD_body1", "velocity")
        size(vtu_vel, 2) == body.ncells || (vtu_vel = collect(vtu_vel'))
        valid = maximum(abs.(Utri .- vtu_vel))/magVinf
        println("$(name): max|manual tri velocity - VTU velocity|/Vinf = ", valid)

        Jtri_quad = quad_field(Jtri, geom.areas)
        grad_quad = quad_mu_gradient(geom, n_ch, n_span)
        grad_p1 = p1_mu_gradient(body, n_ch, n_span)
        s_quad = calibrate_jump_sign(geom, grad_quad, Jtri_quad, n_ch, n_span)
        s_p1 = calibrate_jump_sign(geom, grad_p1, Jtri_quad, n_ch, n_span)
        println("$(name): jump signs quad=$(s_quad) p1=$(s_p1)")

        res = Dict(
            "tri"  => quad_velocity_CL(geom, U0, 2 .* Jtri_quad, 1.0, Lhat, Dhat),
            "quad" => quad_velocity_CL(geom, U0, grad_quad, s_quad, Lhat, Dhat),
            "p1"   => quad_velocity_CL(geom, U0, grad_p1, s_p1, Lhat, Dhat))
        out[name] = (; body, geom, U0, Jtri, grad_quad, grad_p1, s_quad, res, valid)
    end

    println("\nCL by reconstruction method (quad-integrated):")
    @printf("%8s %12s %12s %12s\n", "variant", "tri", "quad", "p1")
    for name in ("diagAC", "diagBD")
        r = out[name].res
        @printf("%8s %12.6f %12.6f %12.6f\n", name,
            r["tri"].CL, r["quad"].CL, r["p1"].CL)
    end
    @printf("%8s %12.6f %12.6f %12.6f\n", "|AC-BD|",
        abs(out["diagAC"].res["tri"].CL - out["diagBD"].res["tri"].CL),
        abs(out["diagAC"].res["quad"].CL - out["diagBD"].res["quad"].CL),
        abs(out["diagAC"].res["p1"].CL - out["diagBD"].res["p1"].CL))
    println("(solver CLs: diagAC=", finals["diagAC"].CL,
        " diagBD=", finals["diagBD"].CL, "; experiment 0.238)")

    println("\nCross-variant velocity agreement, mean |U_AC - U_BD| / Vinf per method:")
    for m in ("tri", "quad", "p1")
        dU = out["diagAC"].res[m].U .- out["diagBD"].res[m].U
        du = vec(sqrt.(sum(abs2, dU; dims=1)))
        @printf("  %-4s mean=%.5g  max=%.5g\n", m, mean(du)/magVinf, maximum(du)/magVinf)
    end

    # ---- Decomposition: is the residual ΔCL carried by U_PV or by ∇μ? ----
    # Cross the inputs: CL(U0 of variant v, quad gradient of variant w).
    # The quads are shared, so gradient fields are directly transplantable.
    println("\nΔCL decomposition with the quad gradient, CL(U0_v, grad_w):")
    AC, BD = out["diagAC"], out["diagBD"]
    s = AC.s_quad
    @printf("%14s %12s %12s\n", "", "grad_AC", "grad_BD")
    for (vn, v) in (("U0_AC", AC), ("U0_BD", BD))
        cl1 = quad_velocity_CL(v.geom, v.U0, AC.grad_quad, s, Lhat, Dhat).CL
        cl2 = quad_velocity_CL(v.geom, v.U0, BD.grad_quad, s, Lhat, Dhat).CL
        @printf("%14s %12.6f %12.6f\n", vn, cl1, cl2)
    end

    # ---- PV-free (Dirichlet/Katz-Plotkin-style) velocity ----
    # Exterior tangential velocity from the doublet sheet alone:
    # U_t = (I - n̂n̂ᵀ)U∞ + full ∇μ (twice the half-jump). No panel-induced
    # PV evaluation anywhere; tests whether the near-singular on-surface PV
    # evaluation is what carries the remaining diagonal dependence.
    println("\nPV-free velocity (tangential U∞ + full ∇μ):")
    pvfree = Dict{String, Any}()
    for (name, v) in (("diagAC", AC), ("diagBD", BD))
        Nq = length(v.geom.A)
        Uinf_t = zeros(3, Nq)
        for q in 1:Nq
            nhat = view(v.geom.normal, :, q)
            Uinf_t[:, q] .= Vinf .- LA.dot(Vinf, nhat) .* nhat
        end
        pvfree[name] = quad_velocity_CL(v.geom, Uinf_t, v.grad_quad, 2*v.s_quad,
                                        Lhat, Dhat)
        @printf("  %s: CL = %.6f  CD = %.6f\n", name,
            pvfree[name].CL, pvfree[name].CD)
    end
    @printf("  |ΔCL| PV-free = %.6f\n",
        abs(pvfree["diagAC"].CL - pvfree["diagBD"].CL))
    # PV-free with the OTHER variant's gradient: any remaining split here is
    # pure γ/quad-geometry sensitivity with no velocity evaluation involved.
    let
        Nq = length(AC.geom.A)
        Uinf_t = zeros(3, Nq)
        for q in 1:Nq
            nhat = view(AC.geom.normal, :, q)
            Uinf_t[:, q] .= Vinf .- LA.dot(Vinf, nhat) .* nhat
        end
        cl_cross = quad_velocity_CL(AC.geom, Uinf_t, BD.grad_quad, 2*AC.s_quad,
                                    Lhat, Dhat).CL
        @printf("  PV-free CL(geom_AC, grad_BD) = %.6f (vs grad_AC %.6f)\n",
            cl_cross, pvfree["diagAC"].CL)
    end

    # ---- Where does the quad-method ΔCL live chordwise? ----
    npts = 2*n_ch
    regions = ("TE rows (j=1,npts)", "LE rows (j=n_ch,n_ch+1)", "rest")
    contrib = Dict(r => 0.0 for r in regions)
    for q in 1:length(AC.geom.A)
        j = mod1(q, npts)
        r = (j == 1 || j == npts) ? regions[1] :
            (j == n_ch || j == n_ch + 1) ? regions[2] : regions[3]
        dF = -(AC.res["quad"].p[q] - BD.res["quad"].p[q]) * AC.geom.A[q] .*
             view(AC.geom.normal, :, q)
        contrib[r] += LA.dot(dF, Lhat)/(0.5*rho*magVinf^2*Sref)
    end
    println("\nQuad-method ΔCL (AC-BD) by chordwise region:")
    for r in regions
        @printf("  %-26s %+10.6f\n", r, contrib[r])
    end

    open(joinpath(save_path, "findings_recon_table.md"), "w") do io
        println(io, "\n## Quad-based ∇μ experiment ($(n_ch)x$(n_span), no re-solve)")
        println(io, "\nManual `U∞+PV+½∇μ_tri` reproduces the VTU velocity to ",
            "max ", max(out["diagAC"].valid, out["diagBD"].valid), "·V∞.")
        println(io, "\n| method | CL diagAC | CL diagBD | &#124;ΔCL&#124; | mean &#124;ΔU&#124;/V∞ |")
        println(io, "|---|---|---|---|---|")
        for m in ("tri", "quad", "p1")
            dU = out["diagAC"].res[m].U .- out["diagBD"].res[m].U
            du = mean(vec(sqrt.(sum(abs2, dU; dims=1))))/magVinf
            println(io, "| ", m, " | ", round(out["diagAC"].res[m].CL, digits=6),
                " | ", round(out["diagBD"].res[m].CL, digits=6),
                " | ", round(abs(out["diagAC"].res[m].CL - out["diagBD"].res[m].CL), sigdigits=3),
                " | ", round(du, sigdigits=3), " |")
        end
        println(io, "\nSolver CLs: diagAC=", finals["diagAC"].CL, ", diagBD=",
            finals["diagBD"].CL, "; experiment 0.238.")
    end
    return out
end

recon = reconstruction_comparison()
