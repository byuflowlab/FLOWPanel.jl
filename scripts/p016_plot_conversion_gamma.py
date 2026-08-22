#!/usr/bin/env python3
"""BRAINSTORM/016 §3 case 2: Γ̄(r/R) for the smooth surface-vorticity
conversion versus the legacy edge-jump conversion, on the rotor.

Companion to `p018_plot_gamma.py` (whose loaders this reuses) but NOT a ladder:
the three curves are categorical identities — legacy, smooth `:upstream`,
smooth `:split` — so they get a categorical palette rather than an ordered
colormap, and the two smooth arms are additionally separated by line style
because they very nearly coincide (that coincidence IS the T4 result).

The lower panel is the quantity M2 reduces to ε_Γ: the difference from the
legacy reference, normalized by max|Γ̄_legacy| over the metric band.

Usage
  p016_plot_conversion_gamma.py [--revs LO HI] [--legacy RUN]

`--legacy` selects the reference arm. Default `p018_b0` is the pre-016-wiring
baseline and is CONFOUNDED by a src change (it has no `conversion` field in its
metadata); pass `--legacy p018_conv_legacy` once the same-code-state control
lands to get the unconfounded comparison.
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p018_analyze import (DATA, gamma_profile, load_ct, _interp)  # noqa: E402

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

OUT = os.path.join(DATA, "p016_conversion_gamma")
BAND = (0.3, 0.95)   # M2 metric band
THRESH = 1.0         # ε_Γ pass threshold, percent

# Okabe-Ito subset; validated (lightness band, chroma floor, CVD separation,
# normal-vision floor, contrast) — worst adjacent pair ΔE 11.0 deutan / 25.8
# normal against a light surface.
C_LEGACY = "#0072B2"   # blue
C_UP = "#D55E00"       # vermillion
C_SPLIT = "#009E73"    # bluish green


def span_integral(prof, lo=None, hi=None):
    r = np.array([p[0] for p in prof])
    g = np.array([p[1] for p in prof])
    if lo is not None:
        m = (r >= lo) & (r <= hi)
        r, g = r[m], g[m]
    return float(np.trapezoid(g, r)) if hasattr(np, "trapezoid") else float(np.trapz(g, r))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--revs", nargs=2, type=int, default=[10, 20])
    p.add_argument("--legacy", default="p018_b0")
    a = p.parse_args()
    lo, hi = a.revs

    arms = [
        (f"legacy edge-jump ({a.legacy})", a.legacy, C_LEGACY, "-", 2.4),
        ("smooth surface-vorticity, :upstream", "p016_smooth_up", C_UP, "-", 2.0),
        ("smooth surface-vorticity, :split", "p016_smooth_split", C_SPLIT, "--", 2.0),
    ]

    profs, labels, styles = [], [], []
    for label, run, color, ls, lw in arms:
        try:
            profs.append(gamma_profile(run, load_ct(run), lo, hi))
        except (SystemExit, FileNotFoundError, OSError) as e:
            print(f"  skip {label} ({run}): {e}")
            continue
        labels.append(label)
        styles.append((color, ls, lw))

    if len(profs) < 2:
        raise SystemExit("need at least two arms to compare")

    ref, ref_label = profs[0], labels[0]
    grid = np.array([r for r, _ in ref if BAND[0] <= r <= BAND[1]])
    scale = max(abs(_interp(ref, r)) for r in grid)

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(9, 8), sharex=True,
        gridspec_kw=dict(height_ratios=[2, 1.3], hspace=0.08))

    for i, (prof, label) in enumerate(zip(profs, labels)):
        color, ls, lw = styles[i]
        r = [q[0] for q in prof]
        g = [q[1] for q in prof]
        ax1.plot(r, g, color=color, ls=ls, lw=lw, label=label)
        if i == 0:
            continue
        d = np.array([(_interp(prof, x) - _interp(ref, x)) / scale * 100
                      for x in grid])
        ax2.plot(grid, d, color=color, ls=ls, lw=lw, label=label)
        # Direct label at the right end, so identity is never color-alone.
        # The two smooth arms very nearly coincide, so they are separated by a
        # vertical offset (and by line style) rather than stacked on the curve.
        ax2.annotate(label.split(",")[-1].strip(),
                     xy=(grid[-1], d[-1]), xytext=(6, 10 if i == 1 else -12),
                     textcoords="offset points", color=color, fontsize=8,
                     ha="left", va="center")

    for ax in (ax1, ax2):
        ax.axvspan(BAND[0], BAND[1], color="0.92", zorder=0)
        ax.grid(alpha=0.3)
        ax.set_axisbelow(True)
    ax2.axhline(0, color="0.4", lw=1)
    for s in (+THRESH, -THRESH):
        ax2.axhline(s, color="crimson", ls=":", lw=1.2)
    ax2.text(0.985, -THRESH * 1.02, f"±{THRESH:g}% M2 threshold",
             color="crimson", fontsize=8, va="top", ha="right")
    ax2.set_ylim(-1.25, 1.35)
    ax2.set_xlim(0, 1.12)

    ax1.set_ylabel(r"$\bar\Gamma$  (TE $\mu$ jump)")
    ax1.set_title("016 §3 case 2 — rotor hover: smooth surface-vorticity vs "
                  f"legacy edge-jump\nB0 carrier, revs {lo}–{hi}; shaded = M2 "
                  r"metric band $0.3\leq r/R\leq0.95$", fontsize=11)
    ax1.legend(fontsize=9, loc="lower left")
    ax2.set_ylabel(r"$\Delta\bar\Gamma\,/\,\max_r|\bar\Gamma_{legacy}|$  [%]")
    ax2.set_xlabel(r"$r/R$")

    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, "gamma_conversion.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    # Table view: the figure must never be the only way to read the result.
    print(f"\n=== 016 conversion A/B — Γ̄(r/R), revs {lo}–{hi}, ref = {ref_label}")
    csv_path = os.path.join(OUT, "gamma_conversion.csv")
    band_ref = span_integral(ref, *BAND)
    print(f"  {'arm':38s} {'eps_max':>8s} {'eps_rms':>8s} {'worst r/R':>10s} "
          f"{'band integral':>14s} {'d band':>8s}")
    for i, (prof, label) in enumerate(zip(profs, labels)):
        band = span_integral(prof, *BAND)
        dband = 100 * (band - band_ref) / abs(band_ref)
        if i == 0:
            print(f"  {label:38s} {'—':>8s} {'—':>8s} {'—':>10s} "
                  f"{band:14.6f} {'—':>8s}")
            continue
        sgn = np.array([(_interp(prof, x) - _interp(ref, x)) / scale
                        for x in grid])
        d = np.abs(sgn)
        print(f"  {label:38s} {100*d.max():7.3f}% {100*np.sqrt((d**2).mean()):7.3f}% "
              f"{grid[int(np.argmax(d))]:10.3f} {band:14.6f} {dband:+7.3f}%")

    with open(csv_path, "w") as fh:
        fh.write("r_over_R," + ",".join(
            lbl.replace(",", ";") for lbl in labels) + "\n")
        rs = sorted({r for prof in profs for r, _ in prof})
        for r in rs:
            fh.write(f"{r:.6f}," + ",".join(
                f"{_interp(prof, r):.8f}" for prof in profs) + "\n")

    print(f"\n  figure -> {path}")
    print(f"  table  -> {csv_path}")


if __name__ == "__main__":
    main()
