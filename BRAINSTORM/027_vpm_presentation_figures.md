# 027 — VPM presentation figures

- [ ] Technical Completion
- [ ] Clear-Context Approval

## RESET BRIEF

**Status (2026-08-29, session 2):** five TikZ figures built through 4 feedback rounds,
all current and linked in the notebook (`journals/20260821.md`, `# 20260829` entry): F1
composite + F3a–d. F1/F3a are at r3, F3b/F3c/F3d at r4 (σ annotations + |ω|(r) profiles
moved from 3c into 3d; rVPM σ̇ equation added to 3b) — NONE formally approved right now;
the earlier F3c approval was superseded by the r4 edit. Notebook links still point at the same PNG filenames (contents
regenerated), so no notebook edit is pending. F2 (transport equation) deferred by Ryan.
**What this is:** concept figures for Ryan's presentation explaining how the VPM works.
Ryan supplies figure requests/feedback; agent maintains the ledger, delegates big
revisions to Opus subagents (small ones inline), reviews renders, sends PNGs to Ryan.
Tangential to the CT mission (presentation support).
**NEXT ACTION (handed to fresh agent):** the F3s experimental internet images —
Ryan will direct picks from the shortlist in section "F3s" below; then download, file
into `~/Dropbox/research/notebooks/img/20260829_vpm_presentation/`, record attribution.
**Conventions:** PowerPoint/Keynote target → standalone TikZ `.tex` + 300-dpi white PNG
(`magick -density 300 <fig>.pdf <fig>.png`); ALL figures live in
`~/Dropbox/research/notebooks/img/20260829_vpm_presentation/` (this item = ledger only).
**Style rulings (accumulated, binding — details in Style section + Log):** soft radial
fades = continuous vorticity; compact ball-shaded spheres = discrete particles; no prose
captions in-image but subcaptions (F1) and vector labels (ω, u, Γ_p, σ) allowed;
transition-arrow annotations are derivative equations (ω̇=(ω·∇)u, Γ̇_p=(Γ_p·∇)u,
ω̇=ν∇²ω, σ̇=ν/σ); uniform left↔right separation, F3c is the spacing reference; σ is a
PARTICLE quantity — σ/σ′ annotations and the |ω|(r) profile plots belong to F3d only,
F3c is a bare field picture (r4 ruling); F3c/F3d canvases pinned to identical widths via
the SAME invisible `\path (-1.4135,0) (10.0965,0);` markers now present in BOTH .tex files
(centers register at equal rendered width — keep when editing; verify with `sips -g
pixelWidth`, currently 336.230 vs 336.246 bp).
**Kernel fact:** F1c plots FLOWVPM's `zeta_gauserf` ζ(ρ)=(2π)^{-3/2}e^{-ρ²/2}
(`FLOWVPM.jl/src/FLOWVPM_kernel.jl:51`), NOT `zeta_gaus` or the panel filament reg.

## Figure ledger

| ID | File | Title | Medium | Status |
|----|------|-------|--------|--------|
| F1 | `fig_what_is_a_particle.tex` | What is a vortex particle? (composite a/b/c, subcaptions) | TikZ | r3 revised, awaiting Ryan |
| F2 | — | Vorticity transport equation | TBD | deferred — Ryan undecided |
| F3a | `fig_operator_splitting_a.tex` | Stretching: tube, ω̇=(ω·∇)u | TikZ | r3 revised, awaiting Ryan |
| F3b | `fig_operator_splitting_b.tex` | Stretching: particle, Γ̇_p=(Γ_p·∇)u + σ̇ | TikZ | r4: rVPM σ̇ equation added below arrow, awaiting Ryan |
| F3c | `fig_operator_splitting_c.tex` | Diffusion: continuous, ω̇=ν∇²ω | TikZ | r4: σ arrows + \|ω\|(r) profiles removed, awaiting Ryan |
| F3d | `fig_operator_splitting_d.tex` | Diffusion: particle, σ̇=ν/σ | TikZ | r4: \|ω\|(r) profiles moved in from 3c, awaiting Ryan |
| F3s | — | Experimental supplements (web images) | curated list | shortlist below, awaiting Ryan picks |

## Specs

### F1 — What is a vortex particle?

Three-panel composite (left→right tells the discretization story):
- (a) continuous vorticity field ω(x,t): smooth shaded swirl/patch;
- (b) same field overlaid with particles (positions x_p, arrows Γ_p), caption
  ω(x) ≈ Σ_p Γ_p ζ_σ(x−x_p) — "discretizing the vorticity field";
- (c) single-particle zoom: radially shaded blob (shading ∝ |ω|), σ arrow from center to
  characteristic radius, Γ_p vector through it; inset pgfplots panel of the Gaussian
  ζ_σ(r) vs r/σ with σ marked. Use the FLOWVPM Gaussian normalization
  (verify against FLOWVPM source; state normalization in a comment in the .tex).

### F3 — Operator splitting

- (a) **Vortex stretching** (inviscid substep), fluid-element picture: material vortex
  tube before/after stretching — length ↑, cross-section ↓, ω ↑, **circulation Γ
  conserved** (annotate Γ = const). Callout: in VPM, dΓ_p/dt = (Γ_p·∇)u.
- (b) **Viscous diffusion** (viscous substep), core-spreading: same blob at t and t+Δt —
  σ grows (dσ²/dt = 2ν), peak ω drops, **Γ unchanged**.
- Physics guardrails: (a) must not imply Γ changes; (b) must change only σ.
- Supplements: curated experimental images (dye-visualized stretching, diffusing vortex
  core, Van Dyke Album candidates) with URLs + attribution; download only after Ryan picks.

### F3s — experimental supplement candidates (curated 2026-08-29, nothing downloaded)

**Vortex stretching:**
1. Petitjeans, *Europhysics News* 34(1) 2003 — dye filament stretched between co-rotating disks, thins & intensifies. [PDF](https://www.europhysicsnews.org/articles/epn/pdf/2003/01/epn03105.pdf). Cite journal.
2. Douady, Couder & Brachet, *PRL* 67:983 (1991) — cavitation-bubble visualization of intense stretched filaments in turbulence (iconic first observation). [APS](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.67.983). APS figure, talk use w/ citation OK.
3. Leweke & Williamson, *JFM* 360:85 (1998) — secondary filaments stretched around a vortex pair. [JFM](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/cooperative-elliptic-instability-of-a-vortex-pair/875D8F66F9FD4BEA20F5375B44DB039A).
4. Van Dyke *Album of Fluid Motion* scan (pick a plate): [archive.org](https://archive.org/download/ScientificBooks/An_Album_of_Fluid_Motion_text.pdf). Plates copyrighted; lecture use w/ plate+photographer citation customary.

**Viscous diffusion / core decay:**
1. NASA Wallops colored-smoke wingtip vortex (1990) — **public domain**, credit NASA Langley. [Commons](https://commons.wikimedia.org/wiki/File:Airplane_vortex.jpg).
2. NASA B-727 wingtip smoke generators (1974) — aging/diffusing trailing vortices, **public domain**. [DVIDS](https://www.dvidshub.net/image/689372/b-727-flight-during-vortex-study-with-wingtip-smoke-generators).
3. NASA F/A-18 pair smoke vortices — **public domain**. [DVIDS](https://www.dvidshub.net/image/731688/smoke-generators-show-twisting-paths-wingtip-vortices-behind-two-nasa-dryden-f-18-jets-used).
4. Bhagwat & Leishman rotor tip-vortex core radius vs wake age w/ Lamb–Oseen curve — the quantitative √(νt) counterpart. [figure](https://www.researchgate.net/figure/Measured-growth-of-the-vortex-core-radius-as-a-function-of-wake-age_fig6_245430001); open alt: [Ramasamy et al. NASA PDF](https://rotorcraft.arc.nasa.gov/Publications/files/RamasamyPB2011F_Final_3350_875.pdf).

Scout's pick: Petitjeans (or Douady et al.) for stretching; NASA Wallops photo + Bhagwat–Leishman
core-growth plot for diffusion (one qualitative photo + one quantitative measurement each;
NASA images unambiguously public domain). Awaiting Ryan's selection before downloading.

## Style

Shared TikZ preamble (all five .tex files): `\documentclass[tikz,border=5pt]{standalone}`,
`\pagecolor{white}`, colors `vortprimary` RGB(31,119,180) (blue: vorticity, tubes,
particles), `vortaccent` RGB(255,127,14) (orange: σ arrows, u quiver, zoom lines),
`vortgray` RGB(90,90,90) (transition arrows, subcaptions). Labels `\large`/`\Large`,
main strokes ≥0.8 pt. Notation: ω, Γ_p, σ, ζ_σ, u; bold vectors. Strain quiver =
u=(−x/2,−y/2,z) sampled on a cylindrical shell, depth-keyed opacity, drawn behind bodies.
Soft blobs: stacked concentric fills at r_k = σ√(−2 ln(k/N)) (Gaussian accumulated alpha).
Spheres: `ball color` + faint dark rim (0.8 pt, 55% opacity). F1c registration: one macro
(`\corelen`) drives sphere radius, σ arrow, and plot x-unit; verify by pixel measurement
after layout changes.

## Log

- 2026-08-29: item created; plan at `~/.claude/plans/generate-a-new-brainstorm-flickering-waffle.md`.
  F1, F3a, F3b delegated to Opus subagents; web supplement search launched.
- 2026-08-29 (later): round-1 figures built + linked in notebook `# 20260829` entry; Ryan
  feedback → revision round 1. Standing style rulings from the feedback: (1) images carry NO
  caption text or equations — slide text does that; (2) composites split into standalone
  images; (3) all particles drawn as smooth RBFs — opaque center fading to transparent, no
  outline/rim circle, generous radius (they are NOT compact-support); (4) F1a gets a
  pictorial (unnumbered) color bar for ω; (5) F1c's ζ graph sits above the particle with
  y-axis on particle center and σ tick aligned to the core radius; (6) F3a strain shown as
  an orange 3-D ∇u quiver field, not two big arrows; (7) new F3a-p particle-stretching
  subfigure (rVPM: |Γ_p|↑, σ↓); (8) F3b keeps ν∇²ω label and plot labels only.
- 2026-08-29 (round 2 feedback). AMENDED rulings: soft-fade rendering applies ONLY to
  continuous vorticity; DISCRETE particles are compact shaded spheres (ball style).
  "No text" ⇏ no arrow labels — vector labels (ω, u, Γ_p, σ) stay; mechanism terms go
  over the transition arrows: (ω·∇)u on F3a, (Γ_p·∇)u on F3b, σ′=√(σ²+2νΔt) on F3d.
  F1 recombined into one composite (orange b→c zoom lines back); F1a gets the dashed
  region silhouette to match F1b. Set renumbered F3a–d (tube/particle stretching,
  continuous/particle diffusion); files renamed to match. F3c approved as-is by Ryan
  (rename only). Orange quiver = local strain-flow VELOCITY u (label "u" added).
- 2026-08-29 (round 3 feedback). NEW rulings: center annotations are explicit derivative
  equations close above the transition arrow — $\dot\omega=(\omega\cdot\nabla)u$ (3a),
  $\dot\Gamma_p=(\Gamma_p\cdot\nabla)u$ (3b), $\dot\omega=\nu\nabla^2\omega$ (3c),
  $\dot\sigma=\nu/\sigma$ (3d) — and left↔right separation is uniform across all four
  (3c's spacing is the reference). F1 gets panel subcaptions "(a) continuous vorticity
  field / (b) particle discretization of the vorticity field / (c) radial shape function";
  (a) loses the color bar (ω label moves to the cloud); (b)'s circled particle's Γ_p flips
  to up-left to match (c). 3b: σ labels dropped, Γ_p arrowhead inside sphere, sphere
  enlarged toward 3a's tube scale. "u" quiver label appears on BOTH fields (3a and 3b).
- 2026-08-29 (round 3.5, inline edits by main agent): F3d spheres shrunk 20% per Ryan
  (σ 1.10→0.88, σ′ 1.76→1.41; centers stay at 3c's 1.05/7.20). F3d canvas width pinned to
  F3c's exact MediaBox (336.23 bp / 1401 px both) via invisible `\path (-1.4135,0)
  (10.0965,0);` markers in `_d.tex` only — solved from 3c's MediaBox, 3c untouched;
  stack-verified sphere centers register over 3c's y-axes. Notebook `# 20260829` entry
  links all five current PNGs with r3 captions (entry was link-updated in place at Ryan's
  request; append-only policy otherwise). Session ended with context reset; internet-image
  selection handed to the next agent (see RESET BRIEF).
- 2026-08-29 (round 4, inline edits by main agent): Ryan — strip the σ/σ′ arrows AND the
  |ω| vs r profile plots from F3c; those plots belong with F3d, "which DOES need the σ and
  σ′ annotations". F3c is now blobs + arrow only (no orange anywhere): a pure
  continuous-field picture, consistent with σ being a particle-discretization quantity.
  Both profile panels moved verbatim into `_d.tex` below the spheres (shifted to y=−4.15,
  down from 3c's −3.40, so the enlarged σ′ sphere clears the |ω| axis), keeping their own
  σ/σ′ double-headed arrows; profile half-widths stay 0.62→0.99 (ratio 1.60, same as the
  spheres' 0.88→1.41) and the peak still scales by 0.62/0.99 so the area — Γ_p — is
  conserved. Width pin: 3c's natural box shrank to 302.45 bp once the profiles left, so the
  invisible `\path` markers were MIRRORED into `_c.tex` (previously `_d.tex` only, the
  "3c untouched" note in the r3.5 entry no longer applies); rebuilt widths 336.230 (c) /
  336.246 (d) bp. PNGs regenerated at 300 dpi; both renders visually checked. Note this
  supersedes the round-2 "F3c approved as-is" ruling — F3c is back to awaiting approval.
- 2026-08-29 (round 4b, inline): Ryan — F3b gets a SECOND equation below the transition
  arrow, the rVPM core update
  $\dot\sigma=-\tfrac{1}{5}\tfrac{\sigma}{|\bm\Gamma_p|}[(\bm\Gamma_p\cdot\nabla)\mathbf{u}]
  \cdot\hat{\bm\Gamma}_p$, so 3b now shows both halves of stretching (|Γ_p|↑ AND σ↓).
  Set `\large` (vs `\Large` for the Γ̇_p line above) with the standard white fill patch at
  (6.10,−0.22) anchor=north — at `\Large` it ran into the spheres. Written with the
  subscript-p notation used everywhere else in the figure (Ryan's text wrote bare Γ).
  Pattern now: σ equations appear ONLY on the particle figures (3b, 3d), never on the
  continuous ones (3a, 3c) — same rule that drove the r4 F3c strip.
- 2026-08-29 (round 4c, inline): Ryan on F3d — (i) the σ annotations under the profiles
  must be one-sided rays: tail-less at r=0, arrowhead at 1σ (the old double-headed span
  read as 2σ); (ii) the spheres must use "the same diameter σ used to create the particles
  in 3c (so they match the plot below)". Those two clauses pointed at different numbers —
  3c's blobs are built with σ=0.55/0.88, the profiles inherited σ=0.62/0.99 — so ALL of
  3d was unified onto 3c's 0.55/0.88: sphere radii, profile Gaussian widths, and profile
  σ rays. Legal because the ratio is 1.60 either way, so only the scale changed, and the
  area-conserving peak factor 0.626→0.625 (=0.55/0.88) is unchanged to 3 digits. Net
  effect: **one σ per column drives everything**, and since the profile r axis shares the
  spheres' cm scale, each σ ray is exactly as long as its sphere's radius — check that
  invariant after any future rescale. Spheres shrank 0.88/1.41 → 0.55/0.88 (undoing the
  r3.5 20%-shrink baseline), which freed vertical room, so the profiles moved back up to
  y=−3.40 (the pre-r4 spacing). Width pin holds: 336.230 (c) / 336.246 (d) bp.
