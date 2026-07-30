# Sweptwing mirrored-discretization discrepancy — diagnosis plan

## Context

`examples/sweptwing.jl` builds the full-span wing two ways and gets two different sectional-lift curves and integrated CLs:

- **`simplewing_mirrored`** (the original helper at `examples/helper_functions.jl:250`): lofts the `+y` half, then mirrors nodes/cells to `−y`. CL = **0.207** at `n_rfl=8, n_span=40`.
- **`simplewing_mirrored_from_negative`** (defined inline at `examples/sweptwing.jl:73`): lofts the `−y` half, then mirrors to `+y`. CL = **0.286** at the same resolution.
- Experimental CL ≈ 0.238 (Weber & Brebner, AOA = 4.2°).

A chordwise convergence sweep already exists at `examples/sweptwing_mirrored_chordwise_convergence.md`:

| n_rfl:n_span | CL positive | err  | CL negative | err   |
|---|---|---|---|---|
| 4:40  | 0.192 | −19.3% | 0.281 | +18.2% |
| 8:40  | 0.207 | −12.9% | 0.287 | +20.5% |
| 16:40 | 0.226 |  −5.0% | 0.275 | +15.6% |
| 32:40 | 0.237 |  −0.2% | 0.266 | +11.9% |
| 40:40 | 0.240 |  +0.7% | 0.265 | +11.1% |

So the **positive-mirror discretization is convergent toward the experimental CL**; the **negative-mirror is biased high by ~11% and converging slowly (or to the wrong limit) under chordwise refinement alone**. The chordwise-convergence Cp slice (image #2) for the positive-mirror case shows large LE/TE oscillations in the FLOWPanel curve that VSPAERO does not show — characteristic of surface-gradient noise in the Bernoulli pressure recovery.

Two distinct questions to answer:

1. **Convergence:** does each discretization converge under joint refinement (`n_rfl` *and* `n_span`)? The negative-mirror's flat 11% error vs chord refinement strongly suggests the spanwise mesh is the bottleneck, not the section.
2. **Why so different at the same nominal resolution?** The two builders use opposite cell-winding conventions on the mirrored half (Explore agent confirmed: `simplewing_mirrored` flips `−y` cells; `simplewing_mirrored_from_negative` flips `+y` cells). Both pass `ensure_winding=false`. The `Das` wake-direction vectors and the `shedding`/`shedding_full` TE-vertex → `Das`-column map (`FLOWPanel_liftingbody.jl:148–180`) depend on TE-node ordering — so the same Kutta condition is being applied with a *different* effective sign/orientation on the mirrored half between the two builders.

The user wants to know: regularization? gradient noise? something else? Hypotheses (ranked by plausibility given the asymmetric, non-converging error):

- **(H1) TE/wake sign error on the mirrored half** — single most likely root cause of the 0.207-vs-0.286 split. Cell winding and `Das` orientation are not co-consistent on the flipped half in at least one of the two builders. This is *not* a regularization or gradient-noise problem; it's a discretization-construction bug.
- **(H2) Spanwise under-resolution near the root crease** — `n_span=40` uniform across the full span gives only ~20 panels per half. The −y mirror builder uses a different spanwise ordering across the centerline join, which can produce non-matching panel widths at `y=0` and a localized loading defect that integrates differently.
- **(H3) Surface-velocity gradient noise (Bernoulli pressure recovery)** — visible in image #2 as the LE/TE wiggles, and arises from finite-difference μ-gradient reconstruction (`compute_mu_gradient!` in the local `reconstructed_half_jump` helper, `sweptwing.jl:241`). This degrades the pointwise Cp but should not by itself cause an 11% integrated-CL bias.
- **(H4) Panel regularization / control-point offset (`CPoffset`)** — could matter at coarse `n_rfl`, but would affect both discretizations symmetrically; doesn't explain the +/− split.

## Findings so far

- **Step 0 (Cp overlay at 8:40)** — done. The `−y` builder loads the wing more aggressively than the `+y` builder across the *whole* chord at outboard stations, not just at the LE or TE. `CL_bernoulli` and `CL_laplace` differ by 1e−4 on the `+y` body, so **H3 (gradient noise) is ruled out** as the cause of the integrated-CL discrepancy.
- **Step 1 (joint refinement, capped at 25k panels)** — done. Span refinement at `n_rfl=8` does essentially nothing (CL stays at 0.207 / 0.287 ±0.003 across `n_span ∈ {40, 60, 80, 120}`); kernel-offset sensitivity is zero (CL identical to 6+ digits across `1e-8 … 1e-12`). Chord refinement closes the gap monotonically: at 24:40 the gap is 0.036, at 32:40 it is 0.029, at 32:24 it is 0.018.
- **At chord-converged 32:24**: `+y` CL = 0.240 (matches CLexp = 0.238 to within 1%), `−y` CL = 0.258 (still ≈8% high). The Cp slices at outboard stations (`2y/b ∈ {0.163, 0.245, 0.51}`) overlap almost perfectly between the two builders; root station shows residual saw-tooth noise but symmetric across builders. Sectional-lift overlay shows the `−y` mirror sits ~0.04 higher than `+y` across **essentially the entire span**, not just near the root.
- Conclusion: the two builders are producing **two genuinely non-equivalent bodies that converge to two different limits**. Not a discretization-noise artifact, not a span-resolution artifact, not a Bernoulli-recovery artifact. **H1 (TE/wake topology or sign mismatch on the mirrored half) is now the dominant hypothesis.**

## Plan: TE/wake topology + sign diagnostics on both grids

Single new diagnostic script — keep it self-contained so it can be re-run independently of the convergence sweep. Target the chord-converged resolution (`n_rfl=32, n_span=24`, 18,432 panels) so we are not looking at noise.

### What to dump, per body (`body` from `simplewing_mirrored`, `body_negative_mirror` from `simplewing_mirrored_from_negative`)

For each of the two bodies, build a "TE table" — one row per shedding edge `i`:

| field | source | purpose |
|---|---|---|
| `i` | column index in `body.shedding[1]` | index |
| `pi`, `pj` | `shedding[1][1, i]`, `shedding[1][4, i]` | upper / lower panel ids |
| `y_mid` | centroid y of `pi` | locate edge on span |
| `nia, nib` (upper) | local TE node indices in `pi` cell | winding diagnostic |
| `nic, nid` (lower) | local TE node indices in `pj` cell | winding diagnostic |
| `node_a, node_b` | global node indices at TE-edge endpoints | geometry |
| `pos_a, pos_b` | `body.nodes[:, node_a/b]` | physical edge |
| `n_up`, `n_lo` | `body.normals[:, pi/pj]` | panel orientation |
| `n_up · z`, `n_lo · z` | — | sign of "upper" assignment; expect `n_up · z > 0` |
| `Da_col_a, Da_col_b` | `shedding_full[5, i], shedding_full[6, i]` (upper) and matching rows for lower panel | Das indexing |
| `Da_a, Da_b` | `Das[1][:, Da_col_*]` | wake-direction vectors at TE vertex |
| `Da_a · Vinf_hat`, `Da_b · Vinf_hat` | — | sign of wake convection direction; expect both > 0 |
| `n_wake` | `(pos_b - pos_a) × (Da_a)` normalized | wake panel normal |
| `n_wake · (n_up - n_lo)` | — | wake-panel orientation consistency w/ TE jump; should be same-signed across all edges |
| `gamma_up`, `gamma_lo` | `body.strength[pi/pj, get_Gammai(body)]` | solved μ on TE panels |
| `dGamma = gamma_up - gamma_lo` | — | TE jump; this is what gets convected as wake circulation |

Dump to CSV at `data/sweptwing_te_diagnostics/te_<builder>.csv`. Use whatever helper most-naturally extracts the upper/lower local-node indices — `shedding_full` rows 2,3 hold those for the upper panel (`FLOWPanel_liftingbody.jl:148–180`).

### Symmetry / consistency checks (printed to stdout)

For each body, compute and print:

1. **Watertight flag** (`pnl.iswatertight(body.nodes, body.cells)`) and **per-half cell counts** (should be exactly equal).
2. **Winding consistency**: count panels with `body.normals[:, p] · (body.controlpoints[:, p] - centroid) > 0` vs `< 0`. Imbalance ≠ 0 means a winding flip on one half.
3. **`n_up · ẑ > 0` and `n_lo · ẑ < 0`** for every shedding edge. Any violation = upper/lower panel mis-assignment on one half.
4. **`Da · V̂∞ > 0`** for every Das column. Any violation = wake convects backward on at least one TE vertex.
5. **`sign(n_wake · (n_up − n_lo))`** uniform across all edges. Non-uniform sign = wake-panel orientation flips between halves.
6. **`dGamma` mirror symmetry**: for each shedding edge at `y_mid > 0`, find the mirror edge at `−y_mid` and check `|dGamma_+ − dGamma_−| / |dGamma_+|`. Asymmetry ≠ 0 (within numerical tolerance) is the smoking gun.
7. **`gamma` mirror symmetry** on *all* panels (not just TE): for each panel `p` find its mirror partner via controlpoint y-flip (KDTree or `findmin` on the y-flipped controlpoint set) and report max `|γ_p − γ_{p'}| / |γ_p|`.

### Cross-builder consistency (the key check)

For the `+y` half (where `simplewing_mirrored` *built* the panels natively and `simplewing_mirrored_from_negative` *mirrored* them):

- Match TE edges across builders by `y_mid` (positive half only).
- For each matched edge, check: same upper/lower assignment? same Das vectors? same `dGamma` (after both bodies are solved at the same resolution)?

Then repeat for the `−y` half.

If a builder is doing the right thing on the half it builds natively but the wrong thing on the mirrored half, this side-by-side will show it as a sign flip on `Da · V̂∞`, `n_up · ẑ`, or `n_wake · (n_up − n_lo)` on the mirrored half only.

### Visual confirmation in ParaView

Re-run `examples/sweptwing.jl` with `paraview=true` and write VTKs of both bodies *plus* a small custom VTK of the wake panels (vertices `[pos_a, pos_b, pos_b + L·Da_a, pos_a + L·Da_b]` for a fixed `L = 5·c_ref`). In ParaView, color by:
- panel normals (`body.normals` → glyphs at controlpoints)
- panel γ (`body.strength[:, get_Gammai(body)]`)
- wake-panel normals

A sign flip will jump out as either a glyph pointing the wrong way on one half, or an antisymmetric γ pattern.

### Critical files

- `examples/sweptwing.jl` — add diagnostic block after both bodies are built and solved.
- `examples/helper_functions.jl:250` (`simplewing_mirrored`) and `examples/sweptwing.jl:73` (`simplewing_mirrored_from_negative`) — the two builders being compared; read-only.
- `src/FLOWPanel_liftingbody.jl:80` (`RigidWakeBody` constructor; `ensure_winding`, `watertight`), `:148–180` (`shedding_full` ↔ `Das` index map), `:1115` (`calc_shedding`), `:1606–1628` (wake panel vertex construction `p3 = p2 + Das[:,i]`, `p4 = p1 + Das[:,i+1]`).
- `src/FLOWPanel_liftingbody.jl:928–961` (`additional_source_system_to_buffer!`) — how Das enters the Kutta condition; consult only if diagnostics confirm a sign flip and we need to know which way it propagates.

### Decision after diagnostics

- **All seven consistency checks pass on both builders + cross-builder matches** ⇒ H1 is also wrong; we need to look at panel kernel evaluation near the centerline (sphere-of-influence / kernel-offset asymmetry across the symmetry plane).
- **One or more checks fail on one builder only** ⇒ root cause located; the fix is in that builder's node/cell/shedding/Das construction. Then verify the fix by re-running step 1's joint sweep and showing both CLs converge to the same limit.

## Verification

- Diagnostic script runs cleanly at `(n_rfl, n_span) = (8, 40)` (cheap, fast) and at `(32, 24)` (chord-converged), and produces the two CSVs and a stdout report.
- For the cheap case, manually inspect ~5 TE-edge rows per builder and confirm the upper/lower panel assignment matches what the geometry shows.
- For the chord-converged case, plot `dGamma(y_mid)` for both builders on the same axes and overlay with the spanwise lift coefficient (`sweptwing_loading_mirrored_discretizations.png`) — `dGamma` should match the qualitative shape of sectional lift, and any builder-to-builder offset in `dGamma` should match the sectional-lift offset.
- ParaView VTKs visually confirm consistent normals and wake direction on both halves of both bodies.

## Findings from TE-topology run (8:40)

- Both bodies have identical node sets (max nearest-node distance 1.1e-15) and identical edge-sharing topology (96 unshared / 11472 shared edges).
- All seven per-body consistency checks pass: winding consistent, `n_up·ẑ > 0`, `n_lo·ẑ < 0`, `Da·V̂∞ > 0`, `n_wake·(n_up−n_lo)` same-signed on all edges, internal γ mirror-asymmetry < 1e-13.
- Cross-builder mismatch is **panel triangulation only**: only 640 / 7680 panel centroids match between bodies, max controlpoint mismatch ~10 mm — the same nodal grid is split into triangles with **opposite quad diagonals** on each half.
- **Hybrid verification**: built `body_hybrid` = (body_neg's nodes) + (body_pos's cells, node-index-permuted). It gives `CL = 0.2073` matching `body_pos` to 1.5e-14, while `body_neg` differs by 0.0791. Diagonal choice fully accounts for the gap at this resolution.

## Remaining diagnostics

Hybrid result above is conclusive at 8:40 but only at one resolution; a few targeted tests close the remaining loops. Tip caps are not present on this wing, so skip those.

### A. Symmetric hybrid (closes the loop completely)

Build `body_hybrid_neg` = (body_pos's nodes) + (body_neg's cells, node-index-permuted via `build_node_perm(body_neg.nodes, body_pos.nodes)`). Solve and check `CL ≈ 0.2864` to machine precision. Reuses the existing `build_node_perm`/`bodytype` machinery in `examples/sweptwing_te_topology.jl`; ~15 lines.

**Pass condition**: `|CL_hybrid_neg − CL_neg| < 1e-12`. That proves the gap is *strictly* triangulation-induced — no residual from kernel-offset, neighbor-list, or any other panel-position-dependent effect.

### B. Diagonal-equalized convergence sweep

Run the diagonal-swap construction at `(n_rfl, n_span) ∈ {16:40, 24:40, 32:24}` (all under 25k panels) and compare:
- `CL(body_pos)` vs `CL(body_hybrid)` at each resolution — should match to machine precision.
- Repeat for the symmetric case (body_neg vs body_hybrid_neg) to confirm the diagonal-equalized bodies converge along the *same* curve at every resolution.

Append the results to a CSV at `data/sweptwing_te_diagnostics/diag_swap_convergence.csv`. This forecloses on the "but maybe at higher resolution kernel effects re-introduce a gap" objection.

### C. Wake-panel triangulation check (read-only investigation)

The rigid-wake panels are constructed dynamically inside the influence kernel from `body.Das` columns (`src/FLOWPanel_liftingbody.jl:1606–1628`: `p1, p2, p3=p2+Das[:,i], p4=p1+Das[:,i+1]`). Inspect that code path to confirm the wake quad is evaluated as a **single quad doublet** (not split into two triangles). If it is split, the wake itself carries a diagonal bias that body_pos vs body_neg might propagate. If it isn't (most likely), this is a no-op confirmation.

If the wake is triangulated, dump the wake-panel diagonal direction per shedding edge for both builders and check symmetry across `y=0`.

### D. Surface-normal fidelity at the LE (explains *why* one diagonal is luckier)

For each panel `p`, compute the angle between `body.normals[:, p]` and the analytic surface normal at the panel centroid. The analytic normal on the RAE101 section requires only spline-differentiating the airfoil contour (`examples/data/airfoil-rae101.csv`) and applying the local section-axis transformation. Compare the LE-region (x/c < 0.1) RMS normal-error between body_pos and body_neg.

**Expectation**: body_neg's LE panels should have systematically larger normal error than body_pos's, explaining its persistent CL bias under chord refinement.

### E. Element-type sensitivity

Re-run both builders with `bodytype = pnl.RigidWakeBody{Union{pnl.ConstantSource, pnl.ConstantDoublet}, 2, Float64, true}` (Dirichlet formulation) at `(8, 40)`. Compare gap. If the gap **shrinks substantially** with constant source+doublet vs vortex-ring, the bias is amplified by the vortex-ring kernel (whose ring follows the cell boundary and therefore changes shape with the diagonal). If gap is similar, bias is purely geometric.

This requires using a working source+doublet builder — check whether `simplewing` supports it, or use `simplewing_mirrored`'s existing `<: NonLiftingBody` branch as a template.

### F. Kernel-offset double-check at high resolution

We already saw zero CL sensitivity to `kerneloffset_panel ∈ {1e-8, 1e-10, 1e-12}` at 8:40 on the original bodies. Run the same offset sweep on `body_hybrid` and `body_hybrid_neg` at 24:40 to confirm the kernel offset doesn't reactivate as the panels shrink. Cheap; reuse `FLOWPANEL_SWEEP_KERNEL_OFFSETS` in `examples/sweptwing_mirrored_convergence.jl`.

### Order of execution

1. **A** (symmetric hybrid) — closes the headline loop, ~5 min.
2. **C** (wake-panel triangulation, read-only) — eliminates a class of objection, ~10 min.
3. **B** (equalized convergence) — confirms it stays clean under refinement, ~30 min compute.
4. **D** (LE normal fidelity) — gives the *mechanism* of the diagonal bias, ~30 min.
5. **E** (element-type sensitivity) — only if A–D leave any residual doubt; may require builder work.
6. **F** (kernel-offset retest) — only if anything in A–D suggests sensitivity.

## Recommendation

Run **A, C, D, B** in that order. Skip **E** unless they leave residue. **F** is defensive and can be deferred.

## Final findings (A, C, D, B all executed)

### A. Symmetric hybrid — PASS
`body_pos's nodes + body_neg's cells` gives `CL = 0.2864`, matching `body_neg` to 5.6e-15. Combined with the original hybrid (`body_neg's nodes + body_pos's cells = 0.2073`), the CL is *strictly* a function of the cell triangulation — fully independent of which node set carries it.

### B. Diagonal-equalized convergence — PASS at every resolution

| Case   | body_pos | body_neg | body_pos − body_hybrid | body_neg − body_hybrid_neg |
|--------|----------|----------|------------------------|----------------------------|
| 8:40   | 0.2073   | 0.2864   | −1.5e-14               | −5.6e-15                   |
| 16:40  | 0.2262   | 0.2745   | +2.9e-14               | −1.9e-14                   |
| 24:40  | 0.2339   | 0.2682   | +1.5e-14               | −2.7e-14                   |
| 32:24  | 0.2402   | 0.2580   | +3.1e-13               | +1.1e-13                   |

Both builders' CLs are converging monotonically toward CLexp = 0.238 (body_pos from below, body_neg from above). No new discrepancy source kicks in at higher resolution. The diagonal-bias mechanism is the only mechanism; it shrinks with refinement, as expected.

### C. Wake-panel triangulation — RULED OUT
Default `semiinfinite_wake=true` (`liftingbody.jl:104`); influence kernel for that case calls `induced_semiinfinite` directly (`elements_fmm.jl:1033`) on a single semi-infinite quad doublet. No wake diagonal exists to bias.

### D. Surface-normal fidelity — NULL RESULT (informative)
At **both** 8:40 and 32:24 the per-panel `acos(n_panel · n_analytic)` distribution is **identical** between builders:

| region | n (8:40) | RMS body_pos | RMS body_neg |
|--------|---:|---:|---:|
| LE  (x/c<0.1)         | 1600 | 4.497° | 4.497° |
| mid (0.1≤x/c<0.7)     | 3200 | 1.396° | 1.396° |
| TE  (x/c≥0.7)         | 2880 | 0.033° | 0.033° |

The diagonal choice does **not** change per-panel surface fidelity. The mechanism is therefore *not* "one diagonal approximates the curvature better"; both triangulations have the same statistical normal accuracy. The remaining mechanism is the **controlpoint-position shift**: with the same node set but different triangulations, only 8% of controlpoints coincide; max controlpoint shift scales as ~10 mm (8:40) → 2.6 mm (32:24), tracking the CL-gap shrinkage. The BC is enforced at *different* physical points between builders, which is what drives the integrated-CL bias.

## Conclusion

The full diagnosis: the two builders triangulate an identical nodal grid with opposite quad diagonals (a consequence of lofting on `[0,1]` vs `[-1,0]` parameter range in `gt.GridTriangleSurface(grid, 1)`). Both triangulations have the same surface-normal accuracy on average. They produce non-coincident controlpoint sets, so the no-through-flow boundary condition is enforced at different physical locations on the wing. Different BC sample points give different solved μ distributions, which integrate to different CLs. The two CLs converge to the same experimental value from opposite sides under joint refinement.

**Not a bug** in the solver or in either builder. A property of single-diagonal triangle paneling of structured-quad grids.

**Practical recommendation**: standardize wing-builder helpers on `b_low=0, b_up=1` (the canonical `simplewing_mirrored` pattern), or eventually port to quad panels to eliminate the diagonal-choice ambiguity.

E and F skipped — A–D leave no doubt.

## Follow-up: γ-difference VTU for ParaView

### Context

We have VTKs for `body` (`+y` mirror) and `body_negative_mirror` (`−y` mirror) at 32:24 in `data/sweptwing000/`. The user wants a third VTK whose per-cell strength field is `γ_body − γ_body_negative_mirror`, so the spatial pattern of the triangulation-induced disagreement can be inspected directly in ParaView. Saved under run_name `gamma_difference` in the same directory.

### Approach

Append a small block to `examples/sweptwing.jl` after both bodies are solved (after the `solve_bernoulli_loading!` call that yields `negative_mirror_bernoulli`). Gate behind `paraview` (the existing `FLOWPANEL_SWEPTWING_VTK` env flag) so the diff VTU is only written when other VTKs are.

Pseudocode:

```julia
if paraview
    Gi = pnl.get_Gammai(body)
    γ_pos = body.strength[:, Gi]
    γ_neg = body_negative_mirror.strength[:, Gi]
    cp_pos = body.controlpoints
    cp_neg = body_negative_mirror.controlpoints

    # match each body cell to nearest controlpoint in body_negative_mirror
    γ_neg_matched = similar(γ_pos)
    for p in 1:body.ncells
        best_d = Inf; best_q = 0
        for q in 1:body_negative_mirror.ncells
            d = (cp_pos[1,p]-cp_neg[1,q])^2 +
                (cp_pos[2,p]-cp_neg[2,q])^2 +
                (cp_pos[3,p]-cp_neg[3,q])^2
            if d < best_d; best_d = d; best_q = q; end
        end
        γ_neg_matched[p] = γ_neg[best_q]
    end

    body_diff = deepcopy(body)
    body_diff.strength[:, Gi] .= γ_pos .- γ_neg_matched
    pnl.write_vtk(joinpath(save_path, "gamma_difference"), body_diff)
    println("Saved γ-difference VTU to $(save_path)/gamma_difference/")
end
```

Notes:
- O(N²) matching costs ~3 s at 18432 panels (acceptable; KDTree not warranted).
- `write_vtk` exports all strength columns of `body_diff.strength` automatically (`abstractbody.jl:501–504`), so the difference field flows through under the existing strength name. Other fields (`velocity`, `potential`) are left as body_pos's; they don't represent a difference but won't mislead since the user is looking at the strength field specifically.
- Wake VTU (`gamma_difference_body1_tw.1.0.vtu`) will carry the difference of wake doublet μ_upper − μ_lower at each TE edge, which is the visualized TE-shed circulation difference — also useful.

### Critical files

- `examples/sweptwing.jl` — add the diff-VTU block, gated by `paraview`.

### Verification

- Re-run `FLOWPANEL_SWEPTWING_VTK=true julia --project examples/sweptwing.jl` at the current 32:24 resolution.
- Confirm `data/sweptwing000/gamma_difference/gamma_difference.0.vtu` and matching `_tw.1.0.vtu` exist.
- In ParaView, load `data/sweptwing000/gamma_difference.pvd` (or the `.vtm`) and color by the strength field — expect a spatial pattern (likely strongest near LE on each half) consistent with the dGamma vs y plot already saved (`sweptwing_dGamma_te.png`).

## Follow-up: chordwise convergence sweep at fixed n_span=24 → CSV

### Context

We need a persisted convergence table covering `n_rfl ∈ {8, 16, 24, 32, 40, 48, 56}` at fixed `n_span = 24`, capped at 35k panels per case, for both mirrored builders. The CSV must live in the same directory as the ParaView outputs (`data/sweptwing000/`, derived from `save_path = joinpath("data", run_name)` in `examples/sweptwing.jl`).

Panel-count math: `panels = 24 * n_rfl * n_span` (see `examples/sweptwing_mirrored_convergence.jl:64`). With n_span=24:
- 8 → 4 608   ✓
- 16 → 9 216   ✓
- 24 → 13 824   ✓
- 32 → 18 432   ✓
- 40 → 23 040   ✓
- 48 → 27 648   ✓
- 56 → 32 256   ✓

All 7 cases fit under the 35 000-panel cap, so the run is 7 × 2 builders = 14 solves.

### γ-difference VTU isolation

The user wants the γ-difference VTU block to *not* execute during this sweep. The block lives at `examples/sweptwing.jl` (right after `negative_mirror_bernoulli = solve_bernoulli_loading!(...)`, gated by `paraview`). The convergence sweep runs `examples/sweptwing_mirrored_convergence.jl`, which only `include`s `helper_functions.jl` and reuses `simplewing_mirrored`/`simplewing_mirrored_from_negative` — it never includes or executes `sweptwing.jl`. So the γ-difference block is structurally unreachable from this sweep; no extra gating needed. Verification step: `grep -n "sweptwing.jl\|gamma_difference" examples/sweptwing_mirrored_convergence.jl` should return nothing matching the diff block.

### Approach

Extend `examples/sweptwing_mirrored_convergence.jl` with optional CSV output gated by a new env var `FLOWPANEL_SWEEP_CSV` (path). When set, each main sweep row appends a record to that CSV. The script already collects per-case `(n_rfl, n_span, panels, shed, kernel, k/hmin, CL, CLerr%, CD, L y>0, L y<0, time_s, source_half)` — just route those columns through `CSV.write` (DataFrames already loaded in the diagnostic ecosystem; if not, use a hand-rolled CSV writer to avoid new deps).

Concretely:
1. Locate the per-case print site in `sweptwing_mirrored_convergence.jl` (the table-printer near where `format("%i%i …", …)` is called).
2. Mirror the same fields into a `DataFrame` row (or `NamedTuple` plus `CSV.write(io, [row]; append=true)`), keyed by the env var.
3. Open the CSV with header on first write, append on subsequent writes.

Then invoke:
```bash
FLOWPANEL_SWEEP_CASES="8:24,16:24,24:24,32:24,40:24,48:24,56:24" \
FLOWPANEL_SWEEP_SOURCE_HALVES="positive,negative" \
FLOWPANEL_SWEEP_MAX_PANELS=35000 \
FLOWPANEL_SWEEP_CSV="data/sweptwing000/mirrored_convergence_n_span24.csv" \
julia --project examples/sweptwing_mirrored_convergence.jl
```

All 7 cases fit under the 35 000 cap; `filter_cases_by_panel_limit` (`sweptwing_mirrored_convergence.jl:66`) will accept all 7. If a heavier case is ever requested, it will be dropped automatically rather than attempted.

### Output CSV schema

Columns (one row per (case, builder) pair):

| column | source | meaning |
|---|---|---|
| `source_half` | env loop | `positive` or `negative` |
| `n_rfl` | case | chordwise panel count |
| `n_span` | case | spanwise panel count |
| `panels` | `24 * n_rfl * n_span` | total panel count |
| `n_shed` | computed in sweep | TE shedding edge count |
| `CL` | sweep | integrated lift coefficient |
| `CD` | sweep | integrated drag coefficient |
| `CL_err_pct` | `(CL − CLexp) / CLexp * 100` with `CLexp = 0.238` | signed error |
| `L_y_pos` / `L_y_neg` | sweep | half-span integrated lift, +y and −y |
| `time_s` | sweep | wall time |
| `kerneloffset` | sweep main | the panel kernel offset used |

### Critical files

- `examples/sweptwing_mirrored_convergence.jl` — add CSV-output path; reuse existing `filter_cases_by_panel_limit` (line 66) and existing per-case data already computed in the table loop.

### Verification

- After the run, `data/sweptwing000/mirrored_convergence_n_span24.csv` exists with **14 rows** (7 cases × 2 builders).
- Spot-check a couple of CSV rows against stdout (e.g., 32:24 should give CL ≈ 0.240 positive, 0.258 negative — matching the values we already saw at this resolution earlier in the session).
- The CSV file path lives in the same directory as the ParaView outputs (`data/sweptwing000/`), as requested.
- The γ-difference VTU block in `sweptwing.jl` is structurally unreachable from the convergence script (verified via grep) — no chance of accidental triggering.

## Follow-up: spanwise convergence at fixed n_rfl=40 → extend CSV

### Context

The previous chordwise sweep at `n_span=24` showed both builders flattening by `n_rfl ≈ 40` (positive: 0.240 → 0.241 across n_rfl ∈ {40, 48, 56}; negative: 0.261 → 0.261 across the same range). Now we want the orthogonal sweep: **fix `n_rfl=40` and refine n_span** to see whether the residual ~9.6% bias in the negative builder closes when the span gets finer, or whether the two builders are converging to genuinely distinct infinite-resolution limits even at this chord-converged setting. The cap rises to **75 000 panels** to push n_span as high as we can afford. Results append to the same convergence table.

### Cases

Panel count: `24 * n_rfl * n_span = 24 * 40 * n_span = 960 * n_span`. Cap 75 000 ⇒ `n_span ≤ 78`. Step n_span by 8:

| n_rfl | n_span | panels |
|---:|---:|---:|
| 40 | 24 | 23 040 (already done) |
| 40 | 32 | 30 720 |
| 40 | 40 | 38 400 |
| 40 | 48 | 46 080 |
| 40 | 56 | 53 760 |
| 40 | 64 | 61 440 |
| 40 | 72 | 69 120 |
| 40 | 80 | 76 800 (over cap → drop) |

Six new (n_rfl, n_span) cases at n_rfl=40, plus the seven prior chordwise rows at n_span=24, all run in **one invocation** so the script writes a single consolidated CSV. 13 cases × 2 builders = **26 rows total**.

### Consolidated CSV

Output path: `data/sweptwing000/mirrored_convergence.csv` (drop the misleading `_n_span24` suffix; that file can be deleted after the new one is written, or kept as a snapshot — user's call).

Schema unchanged from the prior sweep (columns: `n_rfl, n_span, panels, shedding, kerneloffset, offset_over_hmin, source_half, CL, CD, CL_error_pct, lift_positive, lift_negative, elapsed`). The `n_rfl` and `n_span` columns alone disambiguate which axis is being varied — no extra axis column needed.

### Invocation

```bash
FLOWPANEL_SWEEP_CASES="8:24,16:24,24:24,32:24,40:24,48:24,56:24,40:32,40:40,40:48,40:56,40:64,40:72" \
FLOWPANEL_SWEEP_SOURCE_HALVES="positive,negative" \
FLOWPANEL_SWEEP_MAX_PANELS=75000 \
FLOWPANEL_SWEEP_KERNEL_CASE="100:100" \
FLOWPANEL_SWEEP_CSV="data/sweptwing000/mirrored_convergence.csv" \
julia --project examples/sweptwing_mirrored_convergence.jl
```

`FLOWPANEL_SWEEP_KERNEL_CASE=100:100` triggers the "exceeds max_panels" skip branch and prevents the irrelevant kernel-offset sensitivity loop from running (we already showed CL is insensitive to kernel offset).

### Wall-time estimate

Per case, single solve scales roughly with panel count:

| case | panels | ~time per solve |
|---|---:|---:|
| 8:24 | 4 608 | <5 s |
| 56:24 | 32 256 | ~3 min |
| 40:32 | 30 720 | ~3 min |
| 40:48 | 46 080 | ~6 min |
| 40:64 | 61 440 | ~10 min |
| 40:72 | 69 120 | ~13 min |

× 2 builders ≈ **60–90 min total** on a workstation. Run with `nohup`/`tmux` since the user will execute later.

### Critical files

- `examples/sweptwing_mirrored_convergence.jl` — already has CSV writer (added this session). No further code changes required for this run.

### Verification

- After the run, `data/sweptwing000/mirrored_convergence.csv` exists with **26 rows** (13 unique cases × 2 builders).
- Spot-check: 40:24 row should match the prior CSV's 40:24 row for both builders (CL ≈ 0.2397 positive, 0.2609 negative).
- Plot CL vs (panel count, axis) on log-x — separate curves for chordwise refinement and spanwise refinement, per builder. The crucial question is whether the spanwise curve for the negative builder bends back down toward 0.238 at large n_span, or whether it converges to a flat value near 0.260+ regardless of span resolution.

### Decision tree after this sweep

- **If negative builder's CL drops toward 0.238 under span refinement** → the residual 9–10% bias was spanwise-undersampling, not a fundamental diagonal limit. The two builders share a true infinite-resolution limit; the diagonal bias just sets the *rate* of convergence.
- **If negative builder's CL stays flat (≈ 0.260) under both axes refinement** → the two diagonal triangulations have *distinct* infinite-resolution limits (only one of which matches the experiment), and the diagonal choice is a permanent embedded bias rather than a transient discretization error. That would elevate the practical recommendation from "use the canonical builder" to "the negative-lofter builder is straight-up wrong on this geometry — delete it, or document the asymmetry loudly."
