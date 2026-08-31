#!/usr/bin/env python3
"""F1 fig_sigma_map + A1 animated cobweb for the sigma/VPM illustration
campaign (plans/sigma_vpm_illustrations_20260827).

The discrete core-size map, composed stretch + viscous diffusion, in
y = sigma^2 normalized by the per-step viscous floor y_floor = 2*nu*dt:

    u_{n+1} = a * u_n + 1,   a = (1 - dtZ)^2

Regime 1 (dtZ < 2): a < 1, fixed point u* = 1/(1-a)  (sigma_eq analogue).
Regime 2 (dtZ > 2): a > 1, geometric divergence; no physical fixed point.

Emits CSVs for the TikZ figure and renders the A1 GIF (matplotlib).
"""
import csv
import math
import os, subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT = os.path.expanduser(
    "~/Dropbox/research/notebooks/img/20260827_sigma_vpm_illustrations")
DATA = os.path.join(OUT, "fig_sigma_map")
os.makedirs(DATA, exist_ok=True)

CASES = {
    # name: (dtZ, u0 list, n iterations, axis max)
    "regime1": (0.3, [8.0], 12, 9.0),
    "regime2": (2.5, [1.5], 5, 40.0),
}


def cobweb(a, u0, n, cap=None):
    """Return cobweb polyline points [(x, y), ...] for u' = a u + 1."""
    pts = [(u0, 0.0)]
    u = u0
    for _ in range(n):
        v = a * u + 1.0
        pts.append((u, v))     # vertical to the map
        pts.append((v, v))     # horizontal to the diagonal
        u = v
        if cap and u > cap:
            break
    return pts


def write_csv(path, rows, header):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)


for name, (dtZ, u0s, n, umax) in CASES.items():
    a = (1.0 - dtZ) ** 2
    # map line and diagonal sampled over the axis range
    xs = np.linspace(0.0, umax, 200)
    write_csv(os.path.join(DATA, f"map_{name}.csv"),
              [(x, a * x + 1.0) for x in xs], ["u", "unext"])
    # cobweb path (single polyline per start)
    rows = []
    for u0 in u0s:
        for x, y in cobweb(a, u0, n, cap=umax * 2):
            rows.append((x, y))
        rows.append(("nan", "nan"))
    write_csv(os.path.join(DATA, f"cobweb_{name}.csv"), rows, ["u", "unext"])
    fp = 1.0 / (1.0 - a) if a < 1 else float("nan")
    print(f"{name}: dtZ={dtZ} a={a:.3f} fixed_point={fp:.3f}")

# multiplier panel: Euler vs exact, sigma and Gamma channels
zs = np.linspace(0.0, 3.0, 301)
write_csv(os.path.join(DATA, "multipliers.csv"),
          [(z, abs(1 - z), abs(1 - 3 * z), math.exp(-z), math.exp(-3 * z))
           for z in zs],
          ["dtZ", "euler_sigma", "euler_gamma", "exact_sigma", "exact_gamma"])

# ---------------------------------------------------------------- A1 GIF ----
import imageio.v2 as imageio

C_MAP = "#0072B2"      # map line (Okabe-Ito blue)
C_COB = "#D55E00"      # cobweb (vermillion)
C_REF = "#707070"

nframes = 13
frames = []
for k in range(nframes):
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.2), dpi=110)
    for ax, (name, (dtZ, u0s, n, umax)) in zip(axes, CASES.items()):
        a = (1.0 - dtZ) ** 2
        xs = np.linspace(0, umax, 100)
        ax.plot(xs, a * xs + 1, color=C_MAP, lw=2,
                label=r"$u' = (1-\Delta tZ)^2\,u + 1$")
        ax.plot(xs, xs, color=C_REF, lw=1, ls="--")
        pts = cobweb(a, u0s[0], min(k, n), cap=umax * 2)
        if len(pts) > 1:
            px, py = zip(*pts)
            ax.plot(px, py, color=C_COB, lw=1.4)
            ax.plot(px[-1], py[-1], "o", color=C_COB, ms=5)
        if a < 1:
            fp = 1 / (1 - a)
            ax.plot(fp, fp, "s", color=C_MAP, ms=6, mfc="white")
            ax.annotate(r"$\sigma_{eq}^2$", (fp, fp), textcoords="offset points",
                        xytext=(8, -14), color=C_MAP)
        ax.set_xlim(0, umax); ax.set_ylim(0, umax)
        ax.set_xlabel(r"$u_n=\sigma_n^2/2\nu\Delta t$")
        ax.set_ylabel(r"$u_{n+1}$")
        reg = "regime 1 (stable)" if a < 1 else "regime 2 (blow-up)"
        ax.set_title(rf"$\Delta tZ={dtZ}$ — {reg}", fontsize=11)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Core-size map  $y_{n+1}=(1-\\Delta tZ)^2 y_n + 2\\nu\\Delta t$"
                 f"   —   step {min(k, n)}", fontsize=11)
    fig.tight_layout()
    path = os.path.join(DATA, f"_a1_frame_{k:02d}.png")
    fig.savefig(path)
    plt.close(fig)
    frames.append(path)

# magick per-frame -delay (centiseconds): imageio>=2.28 treats duration as ms
# and silently writes 0-delay frames. 60 cs/frame, ~2 s extra hold at the end.
gif = os.path.join(OUT, "a1_sigma_map_cobweb.gif")
subprocess.run(["magick", "-delay", "60"] + frames[:-1]
               + ["-delay", "260", frames[-1], "-loop", "0", gif], check=True)
for f in frames:
    os.remove(f)
print("wrote", gif)
