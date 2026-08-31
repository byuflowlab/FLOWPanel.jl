#=##############################################################################
# DESCRIPTION
#   Weber 45deg swept-back wing (AR 5, RAE 101 12%, no taper/twist/dihedral)
#   meshed by lofting ALONG THE SWEEP AXIS instead of along y.
#
#   Since the wing is a pure sheared extrusion z = c*ztilde((x - |y| tanλ)/c),
#   any vertical cut plane is an exact analytic section: a cut along direction
#   (cosφ, sinφ, 0) through a leading-edge point gives the unit airfoil with
#   in-plane chord c_cut(φ) = c / (cosφ - sinφ tanλ) and z UNSCALED (= c*ztilde).
#   φ = 0 is the streamwise RAE 101; φ = -λ is the sweep-normal section
#   (chord c cosλ, i.e. a ~17%-thick relative section).
#
#   The cut angle blends smoothly (smoothstep over transition length L_t) from
#   streamwise at the centerline and the tip to sweep-normal in the mid-span,
#   so the planform is identical to the y-lofted mesh and every node lies
#   exactly on the true surface; only the grid lines differ.
#
#   The script runs a CL convergence ladder (refine until |ΔCL| < tol or the
#   panel cap), logs to data/sweptwing_sweeploft/convergence.csv, overwrites
#   VTK in the same directory each level (so the last/converged level's
#   ParaView state survives), and prints a comparison against the recorded
#   converged values of the old y-lofted +y/-y mirrored meshes.
#
#   Environment knobs:
#     FLOWPANEL_SWEEPLOFT_CASES="24:24,36:36,54:54,80:80,120:120"
#       Case format is n_ch:n_span where n_ch is the number of chordwise
#       panels per side (loop has 2*n_ch panels) and n_span the number of
#       spanwise panels across the full span (must be even).
#     FLOWPANEL_SWEEPLOFT_MAX_PANELS=60000
#     FLOWPANEL_SWEEPLOFT_CL_TOL=5e-4
#     FLOWPANEL_SWEEPLOFT_BACKSLASH_MAX_PANELS=30000  (Krylov+FMM above this)
=###############################################################################

import FLOWPanel as pnl
import FLOWMath
import LinearAlgebra
import CSV
import DataFrames: DataFrame
using DelimitedFiles: readdlm
using Printf

const LA = LinearAlgebra

run_name      = "sweptwing_sweeploft"
save_path     = joinpath("data", run_name)
airfoil_file  = joinpath(pnl.examples_path, "data", "airfoil-rae101.csv")

# ----------------- SIMULATION PARAMETERS (match examples/sweptwing.jl) -------
AOA           = 4.2                         # (deg) angle of attack
magVinf       = 30.0                        # (m/s)
Vinf          = magVinf*[cos(AOA*pi/180), 0, sin(AOA*pi/180)]
rho           = 1.225                       # (kg/m^3)

b             = 98*0.0254                   # (m) span
ar            = 5.0                         # aspect ratio b/c
lambda        = 45.0                        # (deg) sweep
c             = b/ar                        # (m) chord (untapered)
semispan      = b/2
lambda_rad    = lambda*pi/180

Sref          = b^2/ar
c_ref         = b/ar
core_size  = 1e-10                       # matches the old convergence study

CLexp         = 0.238                       # Weber experiment
CL_old_posy   = 0.2394                      # old +y-mirrored y-loft, 40x56
CL_old_negy   = 0.2671                      # old -y-mirrored y-loft, 40x56

# ----------------- AIRFOIL (no GeometricTools) -------------------------------
"""
Read the airfoil contour CSV (TE -> upper -> LE -> lower -> TE) and return
`(xu, zu, xl, zl)` with x ascending on each side.
"""
function load_airfoil(path)
    data = readdlm(path, ',', Float64; skipstart=1)
    x, z = data[:, 1], data[:, 2]
    ile = argmin(x)
    xu, zu = reverse(x[1:ile]), reverse(z[1:ile])
    xl, zl = x[ile:end], z[ile:end]
    issorted(xu) && issorted(xl) ||
        error("Airfoil contour in $(path) is not TE->upper->LE->lower->TE.")
    return (; xu, zu, xl, zl)
end

"""
Sample the unit airfoil at `n_ch` cosine-clustered panels per side (Akima
interpolation of the contour) and return the closed section loop
`(xloop, zloop)` of length `2*n_ch`, ordered TE -> upper -> LE -> lower with a
single shared TE node (j=1) and a single shared LE node (j=n_ch+1).
"""
function section_loop(airfoil, n_ch)
    xs = [(1 - cos(pi*i/n_ch))/2 for i in 0:n_ch]
    zu = FLOWMath.akima(airfoil.xu, airfoil.zu, xs)
    zl = FLOWMath.akima(airfoil.xl, airfoil.zl, xs)
    zu[1] = zl[1] = 0.0                     # sharp LE point
    zu[end] = zl[end] = 0.0                 # sharp TE point
    xloop = vcat(reverse(xs), xs[2:end-1])
    zloop = vcat(reverse(zu), zl[2:end-1])
    return xloop, zloop
end

# ----------------- SWEEP-LOFT MESH GENERATOR ---------------------------------
smoothstep(t) = (t = clamp(t, 0.0, 1.0); 3t^2 - 2t^3)

"""
Cut-plane angle at spanwise position |y|: 0 (streamwise) at the centerline and
the tip, -λ (sweep-normal) in the mid-span, blended by smoothsteps over
transition length `L_t` at each end.
"""
function cut_angle(yabs, semispan, L_t, lambda_rad)
    s = smoothstep(yabs/L_t) * smoothstep((semispan - yabs)/L_t)
    return -lambda_rad*s
end

"""
    sweeploft_wing(n_ch, n_span; L_t=c, core_size=1e-10)

Build the full-span sweep-lofted RigidWakeBody. `n_ch` chordwise panels per
side, `n_span` (even) spanwise panels across the full span;
ncells = 4*n_ch*n_span.
"""
function sweeploft_wing(n_ch::Int, n_span::Int; L_t::Real=c,
                        core_size::Real=1e-10, airfoil=AIRFOIL,
                        mirror_diagonals::Bool=true, swap_diagonals::Bool=false)
    iseven(n_span) || error("n_span must be even for a symmetric mesh; got $n_span")

    xloop, zloop = section_loop(airfoil, n_ch)
    npts = 2*n_ch
    nst = n_span + 1
    ys = range(-semispan, semispan; length=nst)

    nodes = zeros(3, npts*nst)
    nid(k, j) = (k - 1)*npts + j
    for (k, y) in enumerate(ys)
        yabs, sgn = abs(y), (y < 0 ? -1.0 : 1.0)
        # Section built in the +y half then mirrored by flipping y
        phiabs = cut_angle(yabs, semispan, L_t, lambda_rad)
        c_cut = c/(cos(phiabs) - sin(phiabs)*tan(lambda_rad))
        xLE = yabs*tan(lambda_rad)
        for j in 1:npts
            xi = c_cut*xloop[j]
            nodes[1, nid(k, j)] = xLE + xi*cos(phiabs)
            nodes[2, nid(k, j)] = sgn*(yabs + xi*sin(phiabs))
            nodes[3, nid(k, j)] = c*zloop[j]
        end
    end

    # TE points must be strictly monotonic in y or the sections fan over
    te_y = [nodes[2, nid(k, 1)] for k in 1:nst]
    all(diff(te_y) .> 0) ||
        error("TE points not monotonic in y; increase L_t (= $(L_t)).")

    # Quads between consecutive sections, split along a diagonal.
    # Loop (A, D, C, B) = ((k,j), (k+1,j), (k+1,j+1), (k,j+1)) winds outward
    # (on the upper surface: AD ~ +y, DC ~ -x, cross = +z).
    # With `mirror_diagonals` the negative-y strips use the opposite diagonal
    # so the triangulation (not just the nodes) is mirror-symmetric in y,
    # matching the symmetric construction of the old mirrored meshes.
    cells = Matrix{Int}(undef, 3, 4*n_ch*n_span)
    ci = 0
    for k in 1:n_span, j in 1:npts
        jp = mod1(j + 1, npts)
        A, B, C, D = nid(k, j), nid(k, jp), nid(k + 1, jp), nid(k + 1, j)
        if (mirror_diagonals && k <= n_span ÷ 2) ⊻ swap_diagonals
            cells[:, ci += 1] .= (B, A, D)      # diagonal B-D (mirror image)
            cells[:, ci += 1] .= (C, B, D)
        else
            cells[:, ci += 1] .= (A, D, C)      # diagonal A-C
            cells[:, ci += 1] .= (A, C, B)
        end
    end

    # Orientation guard: signed volume must be positive for outward normals
    vol = 0.0
    for ci in axes(cells, 2)
        a = view(nodes, :, cells[1, ci])
        bb = view(nodes, :, cells[2, ci])
        cc = view(nodes, :, cells[3, ci])
        vol += LA.dot(a, LA.cross(bb, cc))/6
    end
    vol > 0 || error("Mesh signed volume $(vol) <= 0; normals point inward.")

    te_nodes = [nid(k, 1) for k in 1:nst]   # already sorted by y (monotonic)
    shedding = pnl.calc_shedding(nodes, cells, te_nodes, zeros(Float64, 3, 0))
    watertight, _ = pnl.iswatertight(nodes, cells)

    body = pnl.RigidWakeBody{pnl.VortexRing, 1, Float64, false}(
        nodes, cells, [shedding];
        watertight, ensure_winding=false, core_size)

    wake_direction = reshape(Vinf ./ magVinf, :, 1)
    for i in eachindex(body.Das)
        body.Das[i] .= repeat(wake_direction, 1, size(body.Das[i], 2))
    end
    return body
end

"""
Verify every node lies exactly on the analytic sheared-extrusion surface:
z must equal c*ztilde((x - |y| tanλ)/c) on the matching side.
"""
function max_surface_residual(body, airfoil, n_ch)
    npts = 2*n_ch
    resid = 0.0
    for ni in axes(body.nodes, 2)
        x, y, z = body.nodes[1, ni], body.nodes[2, ni], body.nodes[3, ni]
        xt = (x - abs(y)*tan(lambda_rad))/c
        xt = clamp(xt, 0.0, 1.0)
        j = mod1(ni, npts)
        side = j <= n_ch + 1 ? (airfoil.xu, airfoil.zu) : (airfoil.xl, airfoil.zl)
        zref = c*FLOWMath.akima(side[1], side[2], [xt])[1]
        if j == 1 || j == n_ch + 1                  # sharp TE/LE nodes
            zref = 0.0
        end
        resid = max(resid, abs(z - zref))
    end
    return resid
end

# ----------------- QUAD-BASED ∇μ SURFACE-VELOCITY RECONSTRUCTION --------------
# The CL split between the two triangulations of the same node set enters
# through the triangle-based ½∇μ half-jump term (see findings.md). The methods
# below recompute ∇μ from diagonal-independent data: the structured quad grid
# underlying both triangulations, or P1 nodal interpolation.

"""
    quad_mesh(n_ch, n_span) -> Matrix{Int} (4 x Nq)

Node indices of each quad, columns ordered (A, D, C, B) = ((k,j), (k+1,j),
(k+1,j+1), (k,j+1)), matching the outward winding of `sweeploft_wing`. Quad q
corresponds to triangle cells 2q-1 and 2q for BOTH diagonal variants, and
q = (k-1)*npts + j with npts = 2*n_ch chordwise quads per spanwise strip.
"""
function quad_mesh(n_ch::Int, n_span::Int)
    npts = 2*n_ch
    nid(k, j) = (k - 1)*npts + j
    quads = Matrix{Int}(undef, 4, npts*n_span)
    q = 0
    for k in 1:n_span, j in 1:npts
        jp = mod1(j + 1, npts)
        quads[:, q += 1] .= (nid(k, j), nid(k + 1, j), nid(k + 1, jp), nid(k, jp))
    end
    return quads
end

"""
Per-quad geometry and strength from the triangulated body: area (sum of the two
triangles), area-weighted centroid and unit normal, and area-weighted average
γ of the two triangles composing each quad.
"""
function quad_geometry(body)
    Nq = body.ncells ÷ 2
    areas = pnl.calc_areas(body)
    gamma = view(body.strength, :, pnl.get_Gammai(body))

    A = zeros(Nq)
    centroid = zeros(3, Nq)
    normal = zeros(3, Nq)
    mu = zeros(Nq)
    tri_centroid = zeros(3)
    for q in 1:Nq
        for (t, w) in ((2q - 1, areas[2q-1]), (2q, areas[2q]))
            tri_centroid .= 0
            for ni in view(body.cells, :, t)
                tri_centroid .+= view(body.nodes, :, ni)
            end
            tri_centroid ./= 3
            A[q] += w
            centroid[:, q] .+= w .* tri_centroid
            normal[:, q] .+= w .* view(body.normals, :, t)
            mu[q] += w * gamma[t]
        end
        centroid[:, q] ./= A[q]
        normal[:, q] ./= LA.norm(view(normal, :, q))
        mu[q] /= A[q]
    end
    return (; A, centroid, normal, mu, areas)
end

"""
    quad_mu_gradient(geom, n_ch, n_span) -> 3 x Nq

Tangential least-squares gradient of μ on the structured quad grid. Stencil of
quad (k,j): chordwise neighbors (k, j∓1) and spanwise neighbors (k∓1, j),
one-sided at the TE seam (no differencing between j=1 and j=npts, where μ
jumps by the shed circulation) and at the tips. Diagonal-independent by
construction.
"""
function quad_mu_gradient(geom, n_ch::Int, n_span::Int)
    npts = 2*n_ch
    qid(k, j) = (k - 1)*npts + j
    grad = zeros(3, npts*n_span)
    for k in 1:n_span, j in 1:npts
        q = qid(k, j)
        neighbors = Int[]
        j > 1      && push!(neighbors, qid(k, j - 1))   # never wrap the TE seam
        j < npts   && push!(neighbors, qid(k, j + 1))
        k > 1      && push!(neighbors, qid(k - 1, j))
        k < n_span && push!(neighbors, qid(k + 1, j))
        m = length(neighbors)
        @assert m >= 2 "Quad ($k,$j) has fewer than 2 stencil neighbors"

        nhat = view(geom.normal, :, q)
        # local tangent basis from the first neighbor offset
        dx1 = geom.centroid[:, neighbors[1]] .- geom.centroid[:, q]
        t1 = dx1 .- LA.dot(dx1, nhat) .* nhat
        t1 ./= LA.norm(t1)
        t2 = LA.cross(nhat, t1)

        M = zeros(m, 2)
        rhs = zeros(m)
        for (row, nq) in enumerate(neighbors)
            dx = geom.centroid[:, nq] .- geom.centroid[:, q]
            M[row, 1] = LA.dot(dx, t1)
            M[row, 2] = LA.dot(dx, t2)
            rhs[row] = geom.mu[nq] - geom.mu[q]
        end
        g = M \ rhs
        grad[:, q] .= g[1] .* t1 .+ g[2] .* t2
    end
    return grad
end

"""
Per-triangle gradient of a linear (P1) field given nodal values, projected on
the triangle plane. Verified exact for linear fields.
"""
function tri_p1_gradient(p1, p2, p3, f1, f2, f3)
    n = LA.cross(p2 .- p1, p3 .- p1)
    A2 = LA.norm(n)
    g = f1 .* (p2 .- p3) .+ f2 .* (p3 .- p1) .+ f3 .* (p1 .- p2)
    return LA.cross(g, n ./ A2) ./ A2
end

"""
    p1_mu_gradient(body, n_ch, n_span) -> 3 x Nq

Node-based P1 gradient: γ is area-weight-averaged onto nodes (with separate
upper/lower registers at the TE-seam nodes so the μ jump across the TE is
never smeared), each triangle gets the exact linear-shape-function gradient,
and quads get the area-weighted average of their two triangles.
Diagonal-independent because nodal values are shared by construction.
"""
function p1_mu_gradient(body, n_ch::Int, n_span::Int)
    npts = 2*n_ch
    nst = n_span + 1
    Nq = body.ncells ÷ 2
    areas = pnl.calc_areas(body)
    gamma = view(body.strength, :, pnl.get_Gammai(body))

    # cell -> chordwise strip j (upper surface: j <= n_ch)
    strip(ci) = mod1((ci + 1) ÷ 2, npts)
    is_upper(ci) = strip(ci) <= n_ch
    te_node = falses(size(body.nodes, 2))
    for k in 1:nst
        te_node[(k - 1)*npts + 1] = true
    end

    nnodes = size(body.nodes, 2)
    val = zeros(nnodes); wt = zeros(nnodes)            # regular register
    val_lo = zeros(nnodes); wt_lo = zeros(nnodes)      # lower-side register (TE nodes)
    for ci in 1:body.ncells
        a, g = areas[ci], gamma[ci]
        for ni in view(body.cells, :, ci)
            if te_node[ni] && !is_upper(ci)
                val_lo[ni] += a*g; wt_lo[ni] += a
            else
                val[ni] += a*g; wt[ni] += a
            end
        end
    end
    mu_node = val ./ max.(wt, eps())
    mu_node_lo = ifelse.(wt_lo .> 0, val_lo ./ max.(wt_lo, eps()), mu_node)

    grad = zeros(3, Nq)
    for q in 1:Nq
        Aq = 0.0
        for t in (2q - 1, 2q)
            n1, n2, n3 = body.cells[1, t], body.cells[2, t], body.cells[3, t]
            f(ni) = (te_node[ni] && !is_upper(t)) ? mu_node_lo[ni] : mu_node[ni]
            gt = tri_p1_gradient(view(body.nodes, :, n1), view(body.nodes, :, n2),
                                 view(body.nodes, :, n3), f(n1), f(n2), f(n3))
            grad[:, q] .+= areas[t] .* gt
            Aq += areas[t]
        end
        grad[:, q] ./= Aq
    end
    return grad
end

"""
Per-quad exterior surface velocity `U∞+PV (area-weighted from the 2 triangles)
+ jump_sign*0.5*gradmu`, Bernoulli gauge pressure, and CL/CD (ForceMonitor
formula F = Σ -p A n̂, WingNormalization).
"""
function quad_velocity_CL(geom, U0, gradmu, jump_sign, Lhat, Dhat)
    Nq = length(geom.A)
    U0_is_quad = size(U0, 2) == Nq          # accept tri-level or quad-level U0
    U = zeros(3, Nq)
    p = zeros(Nq)
    F = zeros(3)
    for q in 1:Nq
        if U0_is_quad
            U[:, q] .= view(U0, :, q)
        else
            a1, a2 = geom.areas[2q-1], geom.areas[2q]
            U[:, q] .= (a1 .* view(U0, :, 2q-1) .+ a2 .* view(U0, :, 2q)) ./ (a1 + a2)
        end
        U[:, q] .+= (jump_sign*0.5) .* view(gradmu, :, q)
        p[q] = 0.5*rho*(magVinf^2 - LA.dot(view(U, :, q), view(U, :, q)))
        F .-= p[q] .* geom.A[q] .* view(geom.normal, :, q)
    end
    F ./= 0.5*rho*magVinf^2*Sref
    return (; U, p, CL=LA.dot(F, Lhat), CD=LA.dot(F, Dhat))
end

# ----------------- SOLVE ONE CASE ---------------------------------------------
function solve_case(n_ch::Int, n_span::Int;
                    backslash_max_panels::Int=30000,
                    mirror_diagonals::Bool=true, swap_diagonals::Bool=false,
                    path=nothing, name=run_name)
    body = sweeploft_wing(n_ch, n_span; mirror_diagonals, swap_diagonals)

    solver = body.ncells <= backslash_max_panels ?
        pnl.Backslash(body) :
        pnl.KrylovSolver(body; backend=pnl.FastMultipoleBackend(),
                         method=:gmres, atol=1e-9, rtol=1e-9, itmax=200)
    solver_name = body.ncells <= backslash_max_panels ? "backslash" : "krylov"

    Dhat = Vinf/LA.norm(Vinf)
    Shat = [0.0, 1.0, 0.0]
    Lhat = LA.cross(Dhat, Shat)

    pressure = pnl.PressureBernoulli(rho; correct_kuttacondition=false)
    force = pnl.ForceMonitor(1, 1; i_frame=-1,
        normalization=pnl.WingNormalization(rho, Sref, c_ref),
        correct_kuttacondition=false, verbose=false)
    spanwise = pnl.SpanwiseLoadingMonitor(n_span, 1;
        components=(lift=Lhat, drag=Dhat),
        span_axis=Shat, per_length=true,
        normalization=pnl.NoSectionalNormalization())

    elapsed = @elapsed pnl.steady!(body, pnl.ReferenceFrame(body), Vinf;
        body_solvers=solver,
        backend=pnl.FastMultipoleBackend(),
        monitors=(pressure, force, spanwise),
        path=path,
        name=name,
        verbose=false)

    CL = LA.dot(force.force[:, 1], Lhat)
    CD = LA.dot(force.force[:, 1], Dhat)
    y = view(body.controlpoints, 2, :)
    lift = [LA.dot(view(force.distributed_force, :, i), Lhat) for i in 1:body.ncells]

    variant = mirror_diagonals ? (swap_diagonals ? "diagBD" : "diagAC") : "uniform"
    result = (; variant, n_ch, n_span, panels=body.ncells,
              shedding=size(body.shedding[1], 2), solver=solver_name,
              CL, CD, CL_error_pct=100*(CL - CLexp)/CLexp,
              lift_positive=sum(lift[y .> 0]), lift_negative=sum(lift[y .< 0]),
              elapsed)
    return result, body, (; pressure, force, spanwise, Lhat, Dhat)
end

# ----------------- CONVERGENCE DRIVER -----------------------------------------
function parse_case_list(s::AbstractString)
    out = Tuple{Int, Int}[]
    for item in split(s, ",")
        isempty(strip(item)) && continue
        parts = split(strip(item), ":")
        length(parts) == 2 || error("Expected n_ch:n_span case, got $(item)")
        push!(out, (parse(Int, parts[1]), parse(Int, parts[2])))
    end
    isempty(out) && error("No cases parsed from $(s)")
    return out
end

function print_result(r)
    @printf("%8s %5d %6d %8d %6d %10s %+10.6f %+9.3f %+11.6f %8.2f\n",
        r.variant, r.n_ch, r.n_span, r.panels, r.shedding, r.solver, r.CL,
        r.CL_error_pct, r.CD, r.elapsed)
    flush(stdout)
end

function main()
    cases = parse_case_list(get(ENV, "FLOWPANEL_SWEEPLOFT_CASES",
        "24:24,36:36,54:54,80:80,120:120"))
    max_panels = parse(Int, get(ENV, "FLOWPANEL_SWEEPLOFT_MAX_PANELS", "60000"))
    cl_tol = parse(Float64, get(ENV, "FLOWPANEL_SWEEPLOFT_CL_TOL", "5e-4"))
    backslash_max_panels = parse(Int,
        get(ENV, "FLOWPANEL_SWEEPLOFT_BACKSLASH_MAX_PANELS", "30000"))

    rm(save_path; force=true, recursive=true)
    mkpath(save_path)

    # The two mirror-symmetric triangulations of the SAME node set: the +y-half
    # quad diagonal is A-C ("diagAC") or B-D ("diagBD"), with the -y half always
    # the mirror image. At 36x36 they give CL ~ 0.2655 / 0.2423 respectively,
    # reproducing the old +y/-y mirrored-mesh split on identical nodes.
    variants = [(name="diagAC", swap=false), (name="diagBD", swap=true)]

    println("Sweep-loft convergence ladder: ", cases)
    @printf("%8s %5s %6s %8s %6s %10s %10s %9s %11s %8s\n",
        "variant", "n_ch", "n_span", "panels", "shed", "solver", "CL", "CLerr%",
        "CD", "time_s")

    results = []
    final = Dict{String, Any}()
    CL_prev = Dict{String, Float64}()
    for (n_ch, n_span) in cases
        npanels = 4*n_ch*n_span
        if npanels > max_panels
            println("Skipping $(n_ch)x$(n_span) ($(npanels) panels > cap $(max_panels)).")
            continue
        end
        deltas = Float64[]
        for v in variants
            r, body, mon = solve_case(n_ch, n_span; backslash_max_panels,
                                      mirror_diagonals=true, swap_diagonals=v.swap,
                                      path=save_path, name=run_name*"_"*v.name)
            print_result(r)
            push!(results, r)
            final[v.name] = (; r, body, mon)
            push!(deltas, abs(r.CL - get(CL_prev, v.name, Inf)))
            CL_prev[v.name] = r.CL
            CSV.write(joinpath(save_path, "convergence.csv"), DataFrame(results))
        end
        if all(d -> d < cl_tol, deltas)
            println("Converged: max |ΔCL| over variants = ", maximum(deltas),
                " < ", cl_tol)
            break
        end
    end

    isempty(results) && error("No cases were run.")

    println("\n#===== CONVERGED CL COMPARISON =====#")
    for v in variants
        haskey(final, v.name) || continue
        CLconv = final[v.name].r.CL
        @printf("sweep-loft %s (%d panels): CL = %.6f  (vs +y %.4f: %+.4f | vs -y %.4f: %+.4f | vs exp %.3f: %+.2f%%)\n",
            v.name, final[v.name].r.panels, CLconv,
            CL_old_posy, CLconv - CL_old_posy,
            CL_old_negy, CLconv - CL_old_negy,
            CLexp, 100*(CLconv - CLexp)/CLexp)
    end
    println("\nVTK of the final (converged) meshes is in $(save_path)/")

    return results, final
end

const AIRFOIL = load_airfoil(airfoil_file)

if !isinteractive() && get(ENV, "FLOWPANEL_SWEEPLOFT_RUN", "true") == "true"
    results, final = main()
end
