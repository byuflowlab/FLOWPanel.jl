#!/usr/bin/env python3
"""BRAINSTORM/018: plot Γ̄(r/R) across a convergence ladder.

M2 (`decision_rules.md`) reduces the circulation distribution to a single
scalar ε_Γ. That scalar cannot say *where* a discrepancy lives — root, tip, or
broad — which is what decides whether a failure is a resolution, clearance, or
tip-vortex problem. This plots the distribution itself, plus the normalized
difference that ε_Γ is the max/RMS of.

Loaders are imported from `p018_analyze.py`; this script adds only plotting
and the span-integral diagnostic.

Usage
  p018_plot_gamma.py sigma|das|nrows|dt|all
  p018_plot_gamma.py sigma --add p018_L2:_s2:20:34:"sigma/c 0.151"

A rung is `run[:segment][:lo:hi[:label]]`; `--add` appends to the named ladder
so L2 / L1_ov3 can join the σ figure when they land, without editing this file.
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p018_analyze import (DATA, gamma_profile, load_ct, stitch, _interp,
                          _mean)  # noqa: E402

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

OUT = os.path.join(DATA, "p018_gamma_ladders")
BAND = (0.3, 0.95)          # M2 metric band
THRESH = 1.0                # ε_Γ pass threshold, percent


# Each rung: (label, run, restart_segment_or_None, rev_lo, rev_hi)
# Reference rung is index `ref` in the list.
LADDERS = {
    "sigma": dict(
        title="σ ladder — Γ̄(r/R) at fixed Das=0.41c, N=4, NT=36",
        ref=0,
        rungs=[
            ("σ/c 0.68  (B0)",  "p018_b0", "p018_b0_s2", 15, 29),
            ("σ/c 0.297 (L1)",  "p018_L1", "p018_L1_s2", 17, 31),
        ],
    ),
    "sigma_matched": dict(
        title="σ ladder — matched window (revs 22–31) robustness check",
        ref=0,
        rungs=[
            ("σ/c 0.68  (B0)",  "p018_b0", "p018_b0_s2", 22, 29),
            ("σ/c 0.297 (L1)",  "p018_L1", "p018_L1_s2", 22, 31),
        ],
    ),
    "das": dict(
        title="Das ladder — Γ̄(r/R) at σ/c=0.68, N=4, NT=36 (matched revs 10–19)",
        ref=1,
        rungs=[
            ("Das 0.205c", "p018_das0p5", None, 10, 19),
            ("Das 0.41c (B0)", "p018_b0", None, 10, 19),
            ("Das 0.82c", "p018_das2p0", None, 10, 19),
            ("Das 1.64c", "p018_das4p0", None, 10, 19),
        ],
    ),
    "nrows": dict(
        title="N ladder — Γ̄(r/R) at σ/c=0.68, Das=0.41c (matched revs 10–19)",
        ref=2,
        rungs=[
            ("N=1  (d=0.6σ)", "p018_nrows1", None, 10, 19),
            ("N=2  (d=1.2σ)", "p018_nrows2", None, 10, 19),
            ("N=4  (B0)", "p018_b0", None, 10, 19),
        ],
    ),
    "dt": dict(
        title="dt pair — Γ̄(r/R) at matched Das=0.205c, σ/c=0.68 (revs 10–19)",
        ref=0,
        rungs=[
            ("NT=36", "p018_das0p5", None, 10, 19),
            ("NT=72", "p018_nt72", None, 10, 19),
        ],
    ),
}


def rung_profile(run, seg, lo, hi):
    if seg:
        rows = stitch(run, seg, 719)
        names = [run, seg]
    else:
        rows = load_ct(run)
        names = [run]
    return gamma_profile(names, rows, lo, hi)


def span_integral(prof, lo=None, hi=None):
    """∫Γ̄ dr over r/R, trapezoid. Links Γ̄(r/R) back to thrust."""
    r = np.array([p[0] for p in prof])
    g = np.array([p[1] for p in prof])
    if lo is not None:
        m = (r >= lo) & (r <= hi)
        r, g = r[m], g[m]
    return float(np.trapezoid(g, r)) if hasattr(np, "trapezoid") else float(np.trapz(g, r))


def plot_ladder(key, spec, extra):
    rungs = list(spec["rungs"]) + extra
    profs, labels = [], []
    for label, run, seg, lo, hi in rungs:
        try:
            profs.append(rung_profile(run, seg, lo, hi))
            labels.append(label)
        except (SystemExit, FileNotFoundError, OSError) as e:
            print(f"  skip {label} ({run}): {e}")
    if len(profs) < 2:
        print(f"  {key}: fewer than 2 rungs available — not plotting")
        return

    ref_i = min(spec["ref"], len(profs) - 1)
    ref = profs[ref_i]
    grid = np.array([r for r, _ in ref if BAND[0] <= r <= BAND[1]])
    scale = max(abs(_interp(ref, r)) for r in grid)

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(9, 8), sharex=True,
        gridspec_kw=dict(height_ratios=[2, 1.3], hspace=0.08))

    colors = plt.cm.viridis(np.linspace(0, 0.85, len(profs)))
    for i, (prof, label) in enumerate(zip(profs, labels)):
        r = [p[0] for p in prof]
        g = [p[1] for p in prof]
        marker = "o" if i == ref_i else None
        ax1.plot(r, g, color=colors[i], lw=2, marker=marker, ms=3,
                 label=label + (" [ref]" if i == ref_i else ""))
        if i == ref_i:
            continue
        d = np.array([(_interp(prof, x) - _interp(ref, x)) / scale * 100
                      for x in grid])
        ax2.plot(grid, d, color=colors[i], lw=2, label=label)

    for ax in (ax1, ax2):
        ax.axvspan(BAND[0], BAND[1], color="0.9", zorder=0)
        ax.grid(alpha=0.3)
    ax2.axhline(0, color="0.4", lw=1)
    for s in (+THRESH, -THRESH):
        ax2.axhline(s, color="crimson", ls="--", lw=1)
    ax2.text(BAND[0] + 0.01, THRESH * 1.15, f"±{THRESH:g}% M2 threshold",
             color="crimson", fontsize=8, va="bottom")

    ax1.set_ylabel(r"$\bar\Gamma$  (TE $\mu$ jump)")
    ax1.set_title(spec["title"], fontsize=11)
    ax1.legend(fontsize=9, loc="best")
    ax2.set_ylabel(r"$\Delta\bar\Gamma\,/\,\max_r|\bar\Gamma_{ref}|$  [%]")
    ax2.set_xlabel(r"$r/R$")
    ax2.set_xlim(0, 1.0)

    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, f"gamma_{key}.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    # Numbers that the figure exists to support.
    print(f"\n=== {key}: {spec['title']}")
    tot_ref = span_integral(ref)
    band_ref = span_integral(ref, *BAND)
    for i, (prof, label) in enumerate(zip(profs, labels)):
        d = np.array([abs(_interp(prof, x) - _interp(ref, x)) / scale
                      for x in grid])
        sgn = np.array([(_interp(prof, x) - _interp(ref, x)) / scale
                        for x in grid])
        # where the worst error sits
        r_worst = grid[int(np.argmax(d))]
        tot = span_integral(prof)
        band = span_integral(prof, *BAND)
        tag = " [ref]" if i == ref_i else ""
        print(f"  {label:22s}{tag}")
        if i != ref_i:
            print(f"     eps_max {100*d.max():6.3f}%  eps_rms {100*np.sqrt((d**2).mean()):6.3f}%"
                  f"   worst at r/R={r_worst:.3f}   sign {'+' if sgn.mean()>0 else '-'}"
                  f" (mean {100*sgn.mean():+.3f}%)")
        print(f"     integral: full span {tot:.6f}"
              f"  band {band:.6f}"
              f"  Δfull {100*(tot-tot_ref)/abs(tot_ref):+.3f}%"
              f"  Δband {100*(band-band_ref)/abs(band_ref):+.3f}%")
    print(f"  figure -> {path}")


def parse_extra(items):
    out = []
    for it in items or []:
        parts = it.split(":")
        run = parts[0]
        seg = parts[1] or None if len(parts) > 1 else None
        lo = int(parts[2]) if len(parts) > 2 and parts[2] else 10
        hi = int(parts[3]) if len(parts) > 3 and parts[3] else 19
        label = parts[4] if len(parts) > 4 else run
        out.append((label, run, seg, lo, hi))
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("ladder", choices=list(LADDERS) + ["all"])
    p.add_argument("--add", action="append",
                   help='extra rung "run[:seg][:lo:hi[:label]]"')
    a = p.parse_args()
    extra = parse_extra(a.add)
    keys = list(LADDERS) if a.ladder == "all" else [a.ladder]
    for k in keys:
        plot_ladder(k, LADDERS[k], extra if k == a.ladder else [])


if __name__ == "__main__":
    sys.exit(main())
