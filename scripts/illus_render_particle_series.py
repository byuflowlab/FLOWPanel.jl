#!/usr/bin/env python3
"""Shared particle-series renderer for the sigma/VPM illustration campaign
(plans/sigma_vpm_illustrations_20260827). Renders a FLOWPanel wake .vtp series
to PNG frames + a looping GIF with fixed camera/colormap so blow-up and fixed
runs play side-by-side.

Color = log10(sigma/sigma_shed), fixed range across all acts. Optional HUD from
a wake-health/death-trajectory CSV (step, max_u, min sigma ratio).

Usage:
  python3 scripts/illus_render_particle_series.py \
      --glob 'data/p018_L2_visc_forensics/p018_L2_visc_wake1_particles.*.vtp' \
      --sigma-shed 2.3006e-3 --nt 36 \
      --hud-csv data/p018_L2_visc_forensics/death_trajectory.csv \
      --out ~/Dropbox/research/notebooks/img/20260827_sigma_vpm_illustrations/actII_blowup \
      --title 'Act II blow-up: viscous L2 (sigma/R=0.0193), explicit Euler'
"""
import argparse, csv, glob, math, os, re, subprocess, sys

import numpy as np
import matplotlib
import pyvista as pv

CLIM = (-0.6, 0.6)            # log10(sigma/sigma_shed) color range: symmetric so
                              # white sits at sigma = sigma_shed; narrow enough that
                              # the healthy bulk (log ~ 0.1-0.3) carries visible
                              # tint; extremes saturate. Fixed across ALL renders
CMAP = "RdBu"                 # diverging: red = shrinking sigma (ignition danger),
                              # white = at shed size, blue = overgrown
GCLIM = (-1.5, 4.0)           # gamma mode: log10(|Gamma|/median of first frame);
                              # wide upper end so decade-scale runaway keeps gradation
WINDOW = (1280, 960)
BAR_X, BAR_Y, BAR_H = 0.90, 0.30, 0.45   # scalar-bar viewport placement


def step_of(path):
    m = re.search(r"\.(\d+)\.vtp$", path)
    return int(m.group(1)) if m else -1


def load_hud(path):
    if not path:
        return {}
    hud = {}
    with open(os.path.expanduser(path)) as f:
        for row in csv.DictReader(f):
            step = int(float(row.get("step", "nan")))
            hud[step] = row
    return hud


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--glob", required=True)
    ap.add_argument("--sigma-shed", type=float, required=True)
    ap.add_argument("--out", required=True, help="output prefix (dir created; frames + .gif)")
    ap.add_argument("--title", default="")
    ap.add_argument("--nt", type=int, default=36, help="steps per revolution (for rev HUD)")
    ap.add_argument("--spinup-steps", type=int, default=0, help="steps before rev 0")
    ap.add_argument("--hud-csv", default=None)
    ap.add_argument("--fps", type=float, default=6.0)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--font-size", type=int, default=18,
                    help="HUD/annotation font size")
    ap.add_argument("--bar-font-size", type=int, default=None,
                    help="scalar-bar title/label font size (default font_size+3)")
    ap.add_argument("--bar-label-gap", type=float, default=0.012,
                    help="viewport-fraction gap between the tick labels and the "
                         "left edge of the scalar bar")
    ap.add_argument("--font-file",
                    default="/System/Library/Fonts/Supplemental/Arial Unicode.ttf",
                    help="TTF used for text; VTK's built-in families have no "
                         "Greek glyphs, so sigma silently vanishes without this")
    ap.add_argument("--bar-title-gap", type=int, default=22,
                    help="px between the scalar-bar title and the bar top")
    ap.add_argument("--bar-title", default=None,
                    help="override the scalar-bar title")
    ap.add_argument("--point-size", type=float, default=2.5)
    ap.add_argument("--camera", default="iso", choices=["iso", "side", "front"])
    ap.add_argument("--zoom", type=float, default=1.0)
    ap.add_argument("--yaw", type=float, default=0.0,
                    help="rotate camera about the rotor (screen-vertical) axis, "
                         "degrees; +90 brings the old screen-left side to front")
    ap.add_argument("--steps", default=None, help="min:max step filter, e.g. 980:996")
    ap.add_argument("--milestone", default=None,
                    help="REV:TEXT - from rev REV on, overlay TEXT (dark red)")
    ap.add_argument("--color-by", default="sigma",
                    choices=["sigma", "gamma", "abssigma"],
                    help="sigma: RdBu on log10(sigma/sigma_shed). gamma: sphere "
                         "glyphs sized 2*sigma, color+opacity from "
                         "log10(|Gamma|/first-frame median). abssigma: Gamma "
                         "arrow glyphs (opaque) + 2*sigma sphere glyphs "
                         "(faint), both colored by |sigma|")
    # --- abssigma mode knobs ---
    ap.add_argument("--sigma-cmap", default="viridis",
                    help="abssigma mode: colormap for |sigma|/sigma_shed")
    ap.add_argument("--sigma-clim", default="0:2.6963",
                    help="abssigma mode: LO:HI color range for |sigma|, in the "
                         "units set by --sigma-clim-units")
    ap.add_argument("--sigma-log", action="store_true",
                    help="abssigma mode: map |sigma| on a log10 axis (LO must "
                         "be > 0); scalar bar is labelled at decades")
    ap.add_argument("--sigma-clim-units", default="shed", choices=["shed", "m"],
                    help="abssigma mode: units of --sigma-clim and of the "
                         "scalar bar - 'shed' = |sigma|/sigma_shed (default), "
                         "'m' = physical metres")
    ap.add_argument("--sphere-opacity", type=float, default=0.03,
                    help="abssigma mode: opacity of the 2*sigma sphere glyphs")
    ap.add_argument("--diam-cap", type=float, default=0.0,
                    help="abssigma mode: cap sphere diameter at N*sigma_shed "
                         "(0 = no cap, true 2*sigma)")
    ap.add_argument("--arrow-ref-pct", type=float, default=99.0,
                    help="abssigma mode: |Gamma| percentile (FIRST frame) used "
                         "as the arrow-length anchor")
    ap.add_argument("--arrow-ref-len", type=float, default=3.0,
                    help="abssigma mode: arrow length (in sigma_shed) drawn at "
                         "the --arrow-ref-pct percentile of |Gamma|")
    ap.add_argument("--arrow-cap", type=float, default=8.0,
                    help="abssigma mode: clip arrow length at N*sigma_shed "
                         "(0 = uncapped)")
    ap.add_argument("--arrow-scale", type=float, default=None,
                    help="abssigma mode: explicit world-length per unit |Gamma| "
                         "(overrides the percentile anchor)")
    ap.add_argument("--rotor-radius", type=float, default=0.11938,
                    help="abssigma mode: rotor radius R [m] (DJI 9443 default); "
                         "sets the arrow-opacity fade scale")
    ap.add_argument("--arrow-opacity-full", type=float, default=1.0,
                    help="abssigma mode: arrows stay 100%% opaque up to "
                         "|Gamma| = N*R")
    ap.add_argument("--arrow-opacity-fade", type=float, default=4.0,
                    help="abssigma mode: opacity reaches its floor at "
                         "|Gamma| = N*R and beyond")
    ap.add_argument("--arrow-opacity-min", type=float, default=0.01,
                    help="abssigma mode: floor opacity for runaway arrows")
    ap.add_argument("--arrow-stride", type=int, default=1,
                    help="abssigma mode: draw every Nth particle's arrow")
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.expanduser(args.glob)), key=step_of)
    if args.steps:
        lo, hi = (int(x) for x in args.steps.split(":"))
        files = [f for f in files if lo <= step_of(f) <= hi]
    files = files[:: args.stride]
    if not files:
        sys.exit(f"no files match {args.glob}")

    hud = load_hud(args.hud_csv)
    out = os.path.expanduser(args.out)
    os.makedirs(out, exist_ok=True)

    # Fixed camera from the union of first+last frame bounds so death spikes
    # don't rescale the view mid-animation.
    bounds = list(pv.read(files[0]).bounds)  # first frame only: death debris
    # flying out of a FIXED frame reads better than a debris-inflated view
    center = [(bounds[0] + bounds[1]) / 2, (bounds[2] + bounds[3]) / 2,
              (bounds[4] + bounds[5]) / 2]
    span = max(bounds[1] - bounds[0], bounds[3] - bounds[2], bounds[5] - bounds[4])

    barfs = args.bar_font_size if args.bar_font_size else args.font_size + 3

    frames = []
    gref = None
    ascale = None
    for path in files:
        step = step_of(path)
        mesh = pv.read(path)
        sig = np.asarray(mesh.point_data["sigma"], dtype=float)
        ratio = np.clip(sig / args.sigma_shed, 1e-6, None)
        mesh.point_data["log10_sigma_ratio"] = np.log10(ratio)

        if args.color_by == "gamma":
            # Splats sized by the physical core (diameter 2*sigma, world
            # units); color+opacity from |Gamma| so the blow-up site lights up.
            gmag = np.linalg.norm(np.asarray(mesh.point_data["gamma"]), axis=1)
            if gref is None:
                gref = np.median(gmag[gmag > 0])
            g = np.log10(np.clip(gmag / gref, 1e-12, None))
            lut = matplotlib.colormaps["RdBu_r"]
            # piecewise norm: white pinned at the healthy median (g=0), blue
            # saturates at GCLIM[0], red at GCLIM[1] decades of runaway
            tnorm = (0.5 + 0.5 * np.clip(g / GCLIM[1], 0.0, 1.0)
                         - 0.5 * np.clip(g / GCLIM[0], 0.0, 1.0))
            rgba = lut(tnorm)
            # opacity ~ Gamma runaway: healthy (<= median) nearly transparent,
            # only decades-above-median particles turn opaque
            rgba[:, 3] = 0.02 + 0.91 * np.clip(g / 3.0, 0.0, 1.0) ** 0.794
            mesh.point_data["rgba"] = (rgba * 255).astype(np.uint8)
            # sphere glyphs at physical core size (diameter 2*sigma), capped
            # at 4*sigma_shed so blown-up cores don't swallow the frame
            mesh.point_data["diam"] = np.clip(2.0 * np.abs(sig), 0.0,
                                              4.0 * args.sigma_shed)
            mesh.point_data.active_scalars_name = "rgba"
            mesh = mesh.glyph(geom=pv.Sphere(theta_resolution=7,
                                             phi_resolution=7),
                              scale="diam", orient=False, factor=1.0)
            bar_title, bar_clim, bar_cmap = "log10(|G|/G0)", (-1.0, 1.0), "RdBu_r"
            bar_annot = {-1.0: f"{GCLIM[0]:g}", -0.5: f"{GCLIM[0]/2:g}",
                         0.0: "0", 0.5: f"+{GCLIM[1]/2:g}", 1.0: f"+{GCLIM[1]:g}"}
        elif args.color_by == "abssigma":
            # |sigma| carries the color for BOTH glyph sets; Gamma vectors are
            # drawn as opaque arrows (length + orientation = Gamma), the
            # physical cores as near-transparent spheres of diameter 2*sigma.
            absig = np.abs(sig)
            mesh.point_data["abs_sigma"] = absig
            slo, shi = (float(x) for x in args.sigma_clim.split(":"))
            sval = (absig if args.sigma_clim_units == "m"
                    else absig / args.sigma_shed)
            if args.sigma_log:
                llo, lhi = math.log10(slo), math.log10(shi)
                t = np.clip((np.log10(np.clip(sval, 1e-300, None)) - llo)
                            / (lhi - llo), 0.0, 1.0)
            else:
                t = np.clip((sval - slo) / (shi - slo), 0.0, 1.0)
            lut = matplotlib.colormaps[args.sigma_cmap]
            rgba = (lut(t) * 255).astype(np.uint8)

            gvec = np.asarray(mesh.point_data["gamma"], dtype=float)
            gmag = np.linalg.norm(gvec, axis=1)
            if ascale is None:
                ref = np.percentile(gmag[gmag > 0], args.arrow_ref_pct)
                ascale = (args.arrow_scale if args.arrow_scale is not None
                          else args.arrow_ref_len * args.sigma_shed / ref)
                print(f"arrow scale {ascale:.4g} "
                      f"(p{args.arrow_ref_pct:g}|Gamma|={ref:.4g})")

            sph = mesh.copy()
            diam = 2.0 * absig
            if args.diam_cap > 0:
                diam = np.clip(diam, 0.0, args.diam_cap * args.sigma_shed)
            sph.point_data["diam"] = diam
            srgba = rgba.copy()
            srgba[:, 3] = int(round(255 * args.sphere_opacity))
            sph.point_data["rgba"] = srgba
            sph.point_data.active_scalars_name = "rgba"
            spheres = sph.glyph(geom=pv.Sphere(theta_resolution=7,
                                               phi_resolution=7),
                                scale="diam", orient=False, factor=1.0)

            arr = mesh.copy()
            alen = gmag * ascale
            if args.arrow_cap > 0:
                alen = np.clip(alen, 0.0, args.arrow_cap * args.sigma_shed)
            arr.point_data["alen"] = alen
            arrgba = rgba.copy()
            # Runaway arrows would otherwise blanket the frame: hold full
            # opacity until |Gamma| = R, then fade linearly to the floor at 4R.
            g_full = args.arrow_opacity_full * args.rotor_radius
            g_fade = args.arrow_opacity_fade * args.rotor_radius
            f = np.clip((gmag - g_full) / max(g_fade - g_full, 1e-30), 0.0, 1.0)
            aopa = 1.0 - (1.0 - args.arrow_opacity_min) * f
            arrgba[:, 3] = np.clip(np.rint(255 * aopa), 0, 255).astype(np.uint8)
            arr.point_data["rgba"] = arrgba
            arr.point_data["gamma"] = gvec
            if args.arrow_stride > 1:
                arr = arr.extract_points(
                    np.arange(0, arr.n_points, args.arrow_stride),
                    adjacent_cells=False).extract_surface()
            arr.point_data.active_scalars_name = "rgba"
            arrows = arr.glyph(geom=pv.Arrow(tip_resolution=6,
                                             shaft_resolution=6,
                                             tip_length=0.3, shaft_radius=0.03,
                                             tip_radius=0.09),
                               scale="alen", orient="gamma", factor=1.0)
            mesh = (arrows, spheres)
            bar_title = ("\u03c3  [m]" if args.sigma_clim_units == "m"
                         else "\u03c3/\u03c3_shed")
            bar_cmap = args.sigma_cmap
            if args.sigma_log:
                # bar runs 0-1 in normalized log space; label the decades
                bar_title += "  (log)"
                bar_clim, bar_annot = (0.0, 1.0), {}
                llo, lhi = math.log10(slo), math.log10(shi)
                d = math.ceil(llo)
                while d <= lhi + 1e-9:
                    v = (d - llo) / (lhi - llo)
                    bar_annot[round(v, 6)] = (f"1e{d:g}" if abs(d) > 2
                                              else f"{10.0**d:g}")
                    d += 1
            else:
                bar_clim, bar_annot = (slo, shi), None
        else:
            # Gaussian splats with per-particle opacity: transparent near
            # sigma = sigma_shed, opaque as sigma blows up or shrinks away.
            logr = mesh.point_data["log10_sigma_ratio"]
            lut = matplotlib.colormaps[CMAP]
            tnorm = np.clip((logr - CLIM[0]) / (CLIM[1] - CLIM[0]), 0.0, 1.0)
            rgba = lut(tnorm)
            rgba[:, 3] = 0.35 + 0.65 * np.clip(np.abs(logr) / 0.7, 0.0, 1.0)
            mesh.point_data["rgba"] = (rgba * 255).astype(np.uint8)
            bar_title, bar_clim, bar_cmap = "log10(sig/sig_shed)", CLIM, CMAP
            bar_annot = None

        pl = pv.Plotter(off_screen=True, window_size=list(WINDOW))
        pl.set_background("white")
        if args.color_by == "abssigma":
            arrows, spheres = mesh
            pl.add_mesh(spheres, scalars="rgba", rgba=True)
            pl.add_mesh(arrows, scalars="rgba", rgba=True)
        elif args.color_by == "gamma":
            pl.add_mesh(mesh, scalars="rgba", rgba=True)
        else:
            pl.add_mesh(
                mesh, scalars="rgba", rgba=True, style="points_gaussian",
                emissive=False, point_size=args.point_size,
            )
        if args.bar_title:
            bar_title = args.bar_title
        # invisible 2-point mesh just to host the scalar bar
        bar = pv.PolyData(np.array([center, center], dtype=float))
        bar.point_data["log10_sigma_ratio"] = np.array(bar_clim)
        pl.add_mesh(
            bar, scalars="log10_sigma_ratio", cmap=bar_cmap, clim=list(bar_clim),
            opacity=0.0, point_size=1,
            annotations=None if bar_annot else bar_annot,
            scalar_bar_args=dict(
                title=bar_title, color="black",
                title_font_size=barfs, label_font_size=barfs,
                vertical=True,
                position_x=BAR_X, position_y=BAR_Y, width=0.042, height=BAR_H,
                n_labels=0 if bar_annot else 4,
            ),
        )
        if bar_annot and args.font_file and os.path.exists(args.font_file):
            tp0 = pl.scalar_bars[bar_title].GetTitleTextProperty()
            tp0.SetFontFamily(4)          # VTK_FONT_FILE
            tp0.SetFontFile(args.font_file)
        if bar_annot:
            # VTK parks its title tight on the bar top, where it collides with
            # the topmost decade label; push it up. (The title must stay a VTK
            # scalar-bar title: standalone text actors drop the Greek glyphs.)
            pl.scalar_bars[bar_title].SetVerticalTitleSeparation(
                args.bar_title_gap)
            # VTK insets the ramp inside its box: ~5 px at the bottom and
            # (20 + title separation) px at the top. Measured; label positions
            # must follow it or they drift off their decades.
            y_bot = BAR_Y + 5.0 / WINDOW[1]
            y_top = BAR_Y + BAR_H - (20.0 + args.bar_title_gap) / WINDOW[1]
            for v, lbl in bar_annot.items():
                t = pl.add_text(lbl, viewport=True, color="black",
                                position=(BAR_X - args.bar_label_gap,
                                          y_bot + v * (y_top - y_bot)))
                tp = t.GetTextProperty()
                if args.font_file and os.path.exists(args.font_file):
                    tp.SetFontFamily(4)
                    tp.SetFontFile(args.font_file)
                else:
                    tp.SetFontFamilyToArial()
                tp.SetFontSize(barfs)
                tp.SetJustificationToRight()
                tp.SetVerticalJustificationToCentered()
        # tick labels to the LEFT of the bar with a real gap (default VTK
        # placement overlaps the bar when custom annotations are used)
        d = span * 2.6 / args.zoom
        if args.camera == "iso":
            off = [-0.35 * d, -0.85 * d, 0.4 * d]
        elif args.camera == "side":
            off = [0.0, -d, 0.0]
        else:
            off = [-d, 0.0, 0.0]
        th = math.radians(args.yaw)   # rotate offset about world x (rotor axis)
        oy = off[1] * math.cos(th) - off[2] * math.sin(th)
        oz = off[1] * math.sin(th) + off[2] * math.cos(th)
        pos = [center[0] + off[0], center[1] + oy, center[2] + oz]
        # x = rotor axis; wake convects +x, so up = -x puts the disk on top
        pl.camera_position = [pos, center, [-1, 0, 0]]
        # A single runaway glyph (arrow length = |Gamma| can reach 1e6 m at a
        # death step) makes VTK auto-fit the depth range around it and clip the
        # entire wake away -> blank frame. Pin the range to the fixed view.
        dist = math.sqrt(sum(o * o for o in (off[0], oy, oz)))
        pl.camera.clipping_range = (0.05 * dist, 3.0 * dist)

        rev = (step - args.spinup_steps) / args.nt
        lines = [args.title] if args.title else []
        lines.append(f"step {step}   rev {rev:.2f}")
        row = hud.get(step)
        if row:
            mu = float(row.get("max_u", "nan"))
            lines.append(f"u_max {mu:.0f} m/s")
        else:
            lines.append(f"N {len(sig)}")
        pl.add_text("\n".join(lines), position="upper_left",
                    font_size=args.font_size, color="black")
        if args.milestone:
            mrev, mtext = args.milestone.split(":", 1)
            if rev >= float(mrev):
                ta = pl.add_text(mtext, position=(0.5, 0.90), viewport=True,
                                 font_size=args.font_size + 1,
                                 color=(0.65, 0.05, 0.05))
                ta.GetTextProperty().SetJustificationToCentered()

        png = os.path.join(out, f"frame_{step:05d}.png")
        pl.screenshot(png)
        pl.close()
        frames.append(png)
        print(f"rendered {png}")

    gif = out.rstrip("/") + ".gif"
    # assemble with imagemagick: per-frame -delay is respected everywhere
    # (imageio>=2.28 treats duration as ms and silently writes 0-delay frames);
    # delay is in centiseconds; final (death) frame holds ~2 s extra
    delay_cs = max(1, round(100.0 / args.fps))
    cmd = (["magick", "-delay", str(delay_cs)] + frames[:-1]
           + ["-delay", str(delay_cs + 200), frames[-1], "-loop", "0", gif])
    subprocess.run(cmd, check=True)
    print(f"wrote {gif} ({os.path.getsize(gif)/1e6:.1f} MB, {len(frames)} frames)")


if __name__ == "__main__":
    main()
