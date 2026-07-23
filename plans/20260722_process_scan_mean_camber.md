# Preserve Mean Camber While Sharpening Scanned Airfoils

## Summary

Modify only the new scan-processing workflow in `examples/rotor_hover_scan_new/process_scan.jl` so each section is reconstructed as a noise-smoothed blunt airfoil anchored at user-confirmed square trailing-edge points. Evaluate that fit on a common cosine grid, preserve its mean camber exactly while smoothly tapering only thickness to zero, and fit the resulting sharp geometry with CST for the existing CSV output.

The existing CST formulation supports a blunt trailing edge through the `x * dz_half` term in `half_cst`. At `x = 1`, the class-function contribution vanishes and the upper/lower ordinates become `±dz/2`. The final sharp fit remains compatible with the existing format by using `dz = 0`.

Normalize using the midpoint of the confirmed trailing-edge points, then canonicalize their normalized positions to exactly `(1,+h_TE)` and `(1,-h_TE)`. Although midpoint-based normalization already makes their mean ordinate zero, explicit canonicalization removes any streamwise tilt or floating-point drift and guarantees that both the blunt and sharp mean camber lines end at `(1,0)`.

Interactive flow per section: the user selects the LE row, lists the rounded-TE rows to exclude, and confirms (or overrides) an automatically suggested square-TE coordinate pair per surface. There are exactly **two approval stops**: one after selection/canonicalization, one after the full reconstruction. Any rejection re-prompts with the previous answers as defaults and re-runs the entire pipeline, so no invariant can be bypassed.

## Existing Implementation Context

**No repository exploration is needed.** The only files touched are `examples/rotor_hover_scan_new/process_scan.jl` (654 lines, self-contained script, modified) and `test/runtests_unit_process_scan.jl` (new). Everything an implementer needs about the script is inventoried below; open only that script when editing it. The required policy reads (`CLAUDE.md`, `agent_policies/WORKFLOW.md`, `agent_policies/TESTING.md`) are repository rules, not exploration; nothing else in `src/` or `test/` is relevant — the script only calls FLOWPanel for mesh I/O (`pnl.meshes2nodes_cells`, `pnl.NonLiftingBody`, `pnl.write_vtk`), all in the untouched alignment stage.

### Function inventory of `examples/rotor_hover_scan_new/process_scan.jl` (line numbers pre-change)

| Lines | Item | Role | Fate |
|---|---|---|---|
| 11–28 | constants | paths, `plane_idx`/`radial_idx`, `NSLICES=1`, `CST_NCOEFFS=8`, `CST_RECONSTRUCTION_POINTS=201`, `SHARP_TE_BLEND_START=0.8`, `N_RFL_SECS=100`, `USE_COSINE_RFL_BINS=true` | keep |
| 39–50 | `MeshSlice`, `AnchoredMeshSlice` | slice structs; `nodes` is 2×N (row 1 axial x, row 2 thickness z) | `AnchoredMeshSlice` fields change (§1) |
| 60–70 | `Airfoil` struct | `r,twist,chord,axial,tangential,x_u,y_u,x_l,y_l` (LE→TE) | keep |
| 74–186 | `require_indices`, `load_body`, `write_body_vtk`, `fit_plane`, `alignment_from_manual_points`, `estimate_slice_width`, `generate_slices` | mesh load/align/slice stage | untouched |
| 188–211 | `rotated_section(ms)` | LE→(0,0), TE-midpoint→(1,0) normalization via translate/rotate/scale; snaps TE row x's to `x_te`; chordwise mask `x_le ≤ x ≤ x_te`; returns `(twist, chord, x, y, original_indices)` | rewrite per §2 (anchors are coordinates; also return transformed anchors) |
| 213–230 | `add_unique!`, `chordwise_bin_edges` | cosine bin edges over `[xmin,xmax]` | keep |
| 232–275 | `sectioned_selig_indices(x, y; …, upper_te_pos, lower_te_pos)` | bins TE→LE, assigns each point upper/lower by continuity from seeds `y[upper_te_pos]`/`y[lower_te_pos]`, picks per-bin extreme, returns Selig-ordered indices (upper TE→LE, then lower LE→TE) | keep logic; seeds/required entries come from appended anchor positions (§2) |
| 277–285 | `required_selig_indices` | maps the three anchor rows into filtered positions | refactor (§2) |
| 287–297 | `selig_ordered_slice` | no callers | delete |
| 299–326 | `Airfoil(ms::AnchoredMeshSlice)` | normalize→bin→split at `argmin(x)`→`sharpen_trailing_edge!`; computes `axial=-chord_center[2]`, `tangential=-chord_center[1]` from `chord_center = 0.5*(le_point+te_point)` with `te_point` from TE **row indices** | rewrite: `te_point = m_TE` from coordinate anchors; sharp geometry comes from §5 pipeline |
| 328–348 | `sharpen_trailing_edge!` | current linear collapse to TE midpoint (destroys camber aft of 0.8) | delete after removing callers |
| 350–389 | `selig_coordinates`, `write_dat`, `write_xy_dat`, `write_indexed_dat`, `read_dat` | I/O helpers; `write_indexed_dat` writes `# row axial thickness` rows 1-based | keep; reuse for all new file output |
| 391–402 | `prompt_slice_row(label, nrows)` | loop-until-valid stdin prompt pattern | extend pattern: add default-value display, Enter-accepts (§1) |
| 404–433 | `select_airfoil_anchors` | writes `section_XX_raw.dat`, shows numbered scatter (labels only if ≤400 pts, `ax.set_aspect("equal")`), prompts 3 row indices | rewrite per §1 prompt order |
| 435–467 | `count_upstream_of_le`, `select_airfoil_anchors_validated` | warn on points with x<0 after rotation; `[r]emove / [c]hoose` prompt | update to coordinate anchors; fold warning into checkpoint 1 (§2) |
| 469–506 | `approve_airfoil` | approve loop with edit-file-and-reload rejection path | replaced by two-stop flow (§3, §7) — remove edit-and-reload |
| 508–553 | `bernstein`, `half_cst`, `cst_coordinates`, `fit_cst_half`, `fit_cst` | CST core: `y = C*S + x*dz_half + LEW-term`; `fit_cst_half` masks `1e-10 < x < 1-1e-10`; `fit_cst` sets `dz = yu[end]-yl[end]`, returns `[w_upper; w_lower; 0.0; dz]` | keep unchanged |
| 555–567 | `coords_to_cst`, `airfoil_to_cst` | Selig→halves→`fit_cst` | keep |
| 569–603 | `approve_cst` | second approve loop reading the editable `.dat` | replaced by §7 |
| 605–618 | `write_cst_table` | CSV writer, schema `r,twist_rad,chord,axial,tangential,wu1..8,wl1..8,leading_edge_weight,dz` | keep unchanged |
| 620–654 | `main` + `PROGRAM_FILE` guard | align → slice → per-section `approve_airfoil`+`approve_cst` → CSV | per-section loop becomes `approve_processed_airfoil` (§8) |

Interaction conventions to imitate: figures via `fig, ax = PythonPlot.subplots()` … `PythonPlot.show()` then `plt.close(fig)`; prompts via `print(...)` + `lowercase(strip(readline(stdin)))` loops; files named with `@sprintf("section_%02d…", index)` under `joinpath(OUTPUT_DIR, "airfoils")`.

- The target directory is currently untracked in a dirty worktree. Preserve all unrelated changes and do not clean or reset anything.
- Do not modify the older duplicate `examples/rotor_hover_scan/process_scan.jl`; this task targets only `rotor_hover_scan_new`.
- The target script already imports `FLOWPanel`, `GeoIO`, `DelimitedFiles`, `LinearAlgebra`, `Printf`, `PythonPlot`, and `Statistics`; no new package is needed. `SHARP_TE_BLEND_START = 0.8` becomes a per-section default, overridable at the final approval stop.
- `leading_edge_weight` is exported as zero in the CST parameter vector; do not change that behavior.
- The production CSV schema must remain:
  `r,twist_rad,chord,axial,tangential,wu1...wu8,wl1...wl8,leading_edge_weight,dz`.
  Downstream code expects this exact layout and the same CST formula. The per-section blend start is baked into the geometry and is NOT added to the CSV.

## Implementation Changes

### 1. Selection data and interactive input

Replace the row-index trailing-edge fields in `AnchoredMeshSlice` with:

```julia
struct AnchoredMeshSlice
    r::Float64
    nodes::Matrix{Float64}
    upper_te::NTuple{2,Float64}
    le_idx::Int
    lower_te::NTuple{2,Float64}
    excluded_indices::Vector{Int}
    blend_start::Float64            # per-section taper start, default SHARP_TE_BLEND_START
end
```

Keep `MeshSlice` and `Airfoil` conceptually unchanged. `Airfoil.x_u/y_u/x_l/y_l` remain leading-edge-to-trailing-edge vectors and will hold the approved sharp target geometry.

Continue writing `section_XX_raw.dat` with 1-based row numbers before prompting. **Prompt order** (chosen so exclusions are known before the TE suggestion is computed):

1. Leading-edge row.
2. Rounded-TE rows to exclude (blank allowed). These are the scan rows forming the rounded trailing edge; they are removed from all fitting and replaced by the square upper/lower TE anchors below.
3. Upper square-TE raw `(x,z)` coordinate — default is the suggested pair (§1b); Enter accepts.
4. Lower square-TE raw `(x,z)` coordinate — same behavior.

All prompts support **defaults**: a small prompt helper takes a label and an optional previous/suggested value, displays it in the prompt, and returns it on blank input. On any rejection (either checkpoint), re-run the same prompt sequence pre-filled with the previously accepted answers, so fixing one input does not require re-entering the rest. The full normalization/fit pipeline always re-runs after re-prompting.

Input helpers (small, pure, unit-testable):

- Parse one coordinate pair from one line, accepting whitespace- or comma-separated values.
- Parse exclusion rows from blank input or expressions such as `12,15-21,24`; expand inclusive ranges, sort, deduplicate.
- Reject malformed expressions, descending ranges, zero/negative rows, or rows beyond the current slice.
- Prohibit exclusion of `le_idx`.
- Error if exclusions leave fewer than `4*CST_NCOEFFS` retained rows in the slice (guards against a wildcard range wiping the section).
- Require distinct upper/lower TE coordinates and a nonzero distance from the LE to their midpoint.

### 1b. Suggested square-TE anchors

After the exclusion list is accepted, compute one suggested `(x,z)` pair per surface from the retained (non-excluded) points, entirely in raw slice coordinates:

1. Provisional chord frame: direction from the selected LE point to the retained point farthest from it; rotate/translate all retained points into this frame (reuse the rotation pattern of `rotated_section`).
2. Suggested TE plane `x*`: the minimum provisional-`x` among the *excluded* rows (the rounded TE begins where the excluded points begin). If no rows were excluded, use the maximum provisional-`x` among retained rows.
3. Split retained points near the TE (provisional `x/c ≥ 0.85`, tunable local constant) into upper/lower by sign of provisional `y`; require at least 2 points per surface, otherwise fall back to no suggestion (prompt then has no default and requires manual entry — print why).
4. Fit a least-squares line through each surface's near-TE points, evaluate at `x*`, and transform the two resulting points back to raw coordinates. These are the suggested anchors.

Draw the suggested anchors and the `x*` plane on the selection figure so the user judges them visually before accepting. The suggestion is a *default*, never silently applied: the prompt always shows the numeric pair and waits for Enter or an override.

### 2. Normalization and canonical trailing-edge anchors

For each accepted selection:

1. Define the raw trailing-edge midpoint:
   `m_TE = 0.5 * (p_upper + p_lower)`.
2. Define the chord vector from the selected LE point to `m_TE`.
3. Use its norm as the physical chord and its angle as the twist.
4. Translate, rotate, and divide by chord exactly as the current `rotated_section` does, mapping:

   - LE to `(0,0)`.
   - `m_TE` to `(1,0)`.

5. Transform the two TE anchor coordinates with the same operation.
6. Confirm the transformed upper ordinate is greater than the transformed lower ordinate. Treat reversed ordering as invalid selection rather than silently swapping surface labels.
7. Compute the chord-normal half-thickness:
   `h_TE = abs(y_upper_transformed - y_lower_transformed)/2`.
8. Require `h_TE > eps(Float64)`.
9. Replace the transformed anchors with the canonical points:

   - Upper: `(1,+h_TE)`
   - Lower: `(1,-h_TE)`

Using the difference in transformed y ordinates intentionally retains only the endpoint separation projected normal to the chord. Do not use half of the full tilted endpoint distance, since that would convert manual streamwise placement error into excess thickness.

Print a non-fatal warning (and annotate checkpoint 1) if either transformed anchor's `|x − 1|` exceeds `0.02`: this indicates a mis-entered or mis-suggested coordinate, but the user decides at the checkpoint.

Record diagnostics for plotting:

- Originally transformed upper/lower anchors.
- Canonical anchors.
- Each original anchor's x-offset from `1`.
- Each anchor's total displacement during canonicalization.
- `2h_TE`.

Remove `excluded_indices` before applying the existing chordwise `[0,1]` filter. Keep a mapping from retained normalized points back to raw row indices for diagnostics. Append the two canonical anchors separately; do not give them fake raw row numbers.

Refactor `required_selig_indices` and `sectioned_selig_indices` so explicit positions for the appended upper and lower anchors replace the old lookup of TE row indices (append the anchor coordinates to the working `x`/`y` arrays and pass their positions as the required TE entries and continuity seeds). Continue requiring the selected LE point on both ordered surfaces.

Update `Airfoil(ms::AnchoredMeshSlice)` to compute `axial`/`tangential` from the raw `m_TE` (midpoint of the coordinate anchors) instead of TE row indices. Update `count_upstream_of_le` similarly; keep the existing remove-vs-reselect prompt for upstream points.

### 3. First approval stop: selection and canonicalization

After input validation and normalization, show and save a two-panel figure:

- Raw-coordinate panel:

  - All scan points in a neutral color.
  - Retained points and excluded rounded-TE rows distinguished.
  - Raw row numbers.
  - Selected LE.
  - Suggested TE plane `x*` and the accepted upper/lower square-TE coordinates.
  - TE segment, midpoint, and LE-to-midpoint chord line.

- Normalized-coordinate panel:

  - Retained normalized scan points.
  - Points dropped by the chordwise `[0,1]` filter, in a subdued style (these are typically the rounded-TE bulge aft of `x=1`).
  - Retained points with normalized `x` greater than both transformed anchor x-positions **highlighted** as likely missed rounded-TE rows the user may still want to exclude.
  - Originally transformed TE anchor points.
  - Canonical points at `(1,±h_TE)`.
  - Lines connecting original and canonical positions.
  - LE and chord line.

Annotate normalized TE thickness, original x-offsets, and correction distances. Use equal aspect ratio in both panels.

Prompt for approval. Rejection re-prompts the selection sequence with previous answers as defaults. Save the approved figure as `section_XX_checkpoint_1_selection.png` beside the section data.

### 4. Intermediate blunt CST reconstruction

Remove the call to the current `sharpen_trailing_edge!`; delete that function after all callers are removed.

Use the existing chordwise binning and surface-continuity logic to select ordered upper and lower points from retained scan data, with the canonical anchors forced into the appropriate surfaces.

Build LE-to-TE arrays satisfying:

- Both surfaces begin at the selected LE.
- Both surfaces end at `x=1`.
- Upper ends at `+h_TE`.
- Lower ends at `-h_TE`.

Fit each surface with the existing eight-coefficient CST least-squares model:

- Set `dz = 2h_TE`.
- Pass `+dz/2` to the upper half fit and `-dz/2` to the lower half fit.
- Continue excluding exact endpoints from the coefficient solve; CST enforces those endpoints analytically.
- Error clearly if either surface has fewer than `CST_NCOEFFS + 2` interior points after binning (rank-deficient fit).
- Retain `leading_edge_weight = 0`.

Evaluate both fitted surfaces on one shared cosine grid:

```julia
θ = range(0.0, π; length=CST_RECONSTRUCTION_POINTS)
x_common = @. 0.5 * (1.0 - cos(θ))
```

This grid must run monotonically from exactly `0` to exactly `1`.

### 5. Mean camber and thickness taper

On `x_common`, calculate:

```julia
camber = 0.5 .* (y_upper_blunt .+ y_lower_blunt)
half_thickness = 0.5 .* (y_upper_blunt .- y_lower_blunt)
```

Before tapering, assert:

- `camber[end] ≈ 0`.
- `half_thickness[end] ≈ h_TE`.
- No value of `half_thickness` is materially negative; allow only a small floating-point tolerance and error if the fitted surfaces genuinely cross.

Define a dedicated thickness multiplier using the section's `blend_start` (default `SHARP_TE_BLEND_START`):

- `w = 1` for `x <= blend_start`.
- For the remaining chord:
  `s = (x-blend_start)/(1-blend_start)`.
- Use:
  `w = 1 - (10s^3 - 15s^4 + 6s^5)`.

Construct:

```julia
y_upper_sharp = camber .+ w .* half_thickness
y_lower_sharp = camber .- w .* half_thickness
```

Explicitly set both final ordinates to `0.0` at `x_common[end]`. This only removes floating-point residue because the blunt camber endpoint is analytically zero.

The quintic multiplier preserves thickness exactly through `blend_start`, reaches zero at `1.0c`, and has zero first and second derivatives at both ends of the taper interval.

### 6. Final sharp CST fit

Fit the sharp cosine-grid upper and lower arrays with the existing CST machinery:

- Use the same `CST_NCOEFFS = 8`.
- Obtain `dz = y_upper_sharp[end] - y_lower_sharp[end]`, which must be exactly zero.
- Continue exporting `leading_edge_weight = 0`.
- Reconstruct the final fit with `CST_RECONSTRUCTION_POINTS = 201`.

The final CST representation analytically places both surfaces at `(1,0)`. The final least-squares fit may approximate the interior sharp target and therefore may not reproduce its mean camber exactly away from the endpoints; this discrepancy must be exposed in the final diagnostics.

### 7. Second (final) approval stop: reconstruction diagnostics

The blunt fit, taper, and final CST are deterministic once selection is approved, so they share **one** approval stop. Show all three diagnostic figures (as three windows or sequential figures per the script's existing `PythonPlot.show()` pattern), and save each PNG regardless of the approval outcome so a rejected attempt can still be inspected:

**Figure A — blunt fit (`section_XX_checkpoint_2_blunt_cst.png`), two panels:**

- Geometry panel: selected upper/lower scan samples in separate colors; canonical square-TE anchors; blunt upper/lower CST curves; excluded raw points in a subdued style; a view close enough to make the blunt TE visible.
- Residual panel: upper residuals `y_scan - y_fit(x_scan)`, lower residuals, zero reference line, annotated RMS and maximum absolute residual per surface.

**Figure B — camber preservation and taper (`section_XX_checkpoint_3_taper.png`), two panels:**

- Geometry panel: blunt CST surfaces; tapered sharp target surfaces; blunt mean camber; sharp mean camber visibly overlaid; vertical marker at `x/c = blend_start`.
- Distribution panel: original full thickness `2*half_thickness`; tapered full thickness `2*w*half_thickness`; camber difference `0.5*(y_upper_sharp+y_lower_sharp) - camber`; zero reference line.
- Annotate maximum absolute camber change, final TE thickness, and maximum geometry change forward of `blend_start`.

**Figure C — final CST (`section_XX_checkpoint_4_final_cst.png`), multi-panel:**

- Whole-airfoil geometry: sharp target points and final upper/lower CST reconstruction.
- Mean-camber comparison: sharp target camber, final reconstructed CST camber, and their difference.
- Surface residuals: upper/lower target-minus-fit residuals with RMS and maximum absolute values.
- Trailing-edge close-up: both target surfaces, both reconstructed surfaces, marker at `(1,0)` demonstrating coincident endpoints.

Residual statistics are approval diagnostics, not automatic pass/fail thresholds. The prompt accepts:

- `y` — approve; proceed to output.
- `n` — reject; return to the selection sequence with previous answers as defaults (checkpoint 1 runs again).
- `b` — enter a new `blend_start` in `[0,1)` for this section (default keeps the current value), then re-run only the taper and final fit (§5–§6) and re-show this stop.

### 8. Output and main-loop refactor

Integrate the two approval stops into one per-section routine that returns both:

- The final sharp `Airfoil`.
- Its approved final CST parameter vector.

Update `main` so it no longer calls the old separate `approve_airfoil` followed by `approve_cst`. Instead:

```julia
af, params = approve_processed_airfoil(slice, i)
push!(airfoils, af)
push!(all_params, params)
```

Preserve these outputs:

- `processed/airfoils/section_XX_raw.dat`: all indexed raw scan points.
- `processed/airfoils/section_XX.dat`: approved sharp cosine-grid target in Selig order.
- `processed/dji9443_medium_cst.csv`: approved final CST parameters and unchanged metadata/schema.

Write `section_XX.dat` only after the final approval. Use the existing `selig_coordinates` ordering: trailing edge down the upper surface to the LE, followed by the lower surface back to the trailing edge, without duplicating the LE.

Blunt CST coefficients are diagnostic only and must not be added to the CSV.

## Interfaces and Compatibility

- No FLOWPanel public API changes.
- No new package dependencies.
- No changes to the production CST CSV columns, coefficient count, reconstruction convention, or metadata.
- Final CSV rows must have `dz = 0`.
- The interactive workflow changes from three indexed anchors plus two approvals to:

  - One indexed LE.
  - One rounded-edge exclusion list.
  - Two TE anchors (suggested defaults, overridable).
  - Two visual approval stops (four diagnostic PNGs).
  - Optional per-section blend-start override at the final stop.

- Remove the current edit-and-reload behavior from rejected approvals. Every rejection re-enters the prompt sequence (with defaults) and re-runs the full pipeline, so edits cannot bypass canonicalization or camber-preservation invariants.

## Automated Test Plan

Add a standalone focused test file `test/runtests_unit_process_scan.jl`. Include the target script; its existing `PROGRAM_FILE == @__FILE__` guard prevents `main()` from running. Do not register this test in the already-modified broad `test/runtests.jl`. Run headless:

```bash
MPLBACKEND=Agg julia --project -e 'include("test/runtests_unit_process_scan.jl")'
```

(The include imports PythonPlot; the Agg backend avoids any display requirement.)

Cover:

- Exclusion parser:

  - Blank input.
  - Individual rows.
  - Mixed comma/range expressions.
  - Sorting and deduplication.
  - Malformed, descending, negative, zero, and out-of-range input.
  - Conflict with the selected LE.
  - Too-few-retained-rows error.

- TE suggestion:

  - Synthetic section with straight near-TE surfaces and a fake rounded cap marked excluded: suggested anchors must land on the extrapolated lines at the excluded rows' minimum provisional x.
  - No-exclusion fallback uses the maximum retained provisional x.
  - Fewer than 2 near-TE points on a surface yields no suggestion (sentinel), not an error.

- Normalization and canonicalization:

  - Use deliberately tilted synthetic TE input.
  - Verify the supplied midpoint maps to `(1,0)`.
  - Verify canonical endpoints are exactly `(1,+h_TE)` and `(1,-h_TE)`.
  - Verify their midpoint and equal-distance center are `(1,0)`.
  - Verify `h_TE` comes from normal projection rather than full tilted distance.
  - Verify reversed upper/lower labels and degenerate anchors fail clearly.

- Blunt CST:

  - Fit noisy synthetic upper/lower samples.
  - Verify finite coefficients.
  - Verify reconstructed TE ordinates are exactly `±h_TE`.
  - Verify the reconstructed mean-camber endpoint is zero.
  - Verify the too-few-interior-points error triggers.

- Cosine grid and taper:

  - Shared, monotonic grid containing exact endpoints.
  - Geometry unchanged for all `x <= blend_start` (test at 0.8 and one non-default value, e.g. 0.7).
  - Mean camber unchanged to floating-point tolerance.
  - Nonnegative tapered thickness.
  - Zero thickness at `x=1`.
  - Quintic multiplier and its first derivative match at both blend boundaries; also check the analytically expected zero second derivatives.

- Final CST:

  - Finite coefficients.
  - `dz == 0`.
  - Both reconstructed surfaces end at `(1,0)`.
  - Residual and camber-error diagnostic calculations return finite values.

## Interactive Verification

Run the script for its configured one-section case:

```bash
julia --project examples/rotor_hover_scan_new/process_scan.jl
```

During the smoke test:

1. Enter the LE row and a rounded-TE exclusion list; accept the suggested TE anchors with Enter.
2. Reject the first stop once and verify the re-prompt shows the previous answers as defaults; change only the exclusion list and confirm the rest carries over.
3. At the final stop, use `b` once to change `blend_start` and confirm only Figures B/C change; then approve.
4. Confirm the raw `.dat`, sharp `.dat`, four checkpoint PNGs, and CST CSV are written.
5. Inspect the final CSV header and verify its schema is unchanged and its row has `dz = 0`.

## Acceptance Criteria

- The confirmed TE midpoint defines normalized `(1,0)`.
- Canonical blunt anchors are exactly `(1,±h_TE)`.
- Rounded-edge scan rows do not participate in binning or fitting, and checkpoint 1 highlights likely missed rounded-TE rows.
- The suggested TE anchors are visibly drawn before confirmation and never applied without user acceptance.
- The blunt CST smooths scan noise while enforcing the canonical TE.
- The sharp target preserves blunt-fit camber pointwise and modifies only thickness.
- Geometry is unchanged through `blend_start` and smoothly closes thereafter.
- The sharp target and final CST both end at `(1,0)`.
- Both approval stops expose sufficient geometry and error information for informed approval, and rejections re-prompt with defaults.
- Existing CSV consumers require no changes.

## Assumptions

- TE coordinates (typed or suggested) are in the raw two-dimensional `(axial, thickness)` coordinate system displayed in the selection plot.
- Exclusion indices refer to 1-based rows in `section_XX_raw.dat`; comma-separated indices and inclusive ranges are supported.
- "Square trailing edge" means perpendicular to the normalized chord, so TE thickness is derived from endpoint separation projected onto the chord-normal direction.
- The user-labelled upper point must transform above the lower point.
- Residual statistics are informational and user-approved; no new hard fit-error threshold is introduced.
- A rejected stop re-runs the complete selection and reconstruction of that section, with prior answers offered as prompt defaults.
