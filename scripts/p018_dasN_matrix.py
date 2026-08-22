#!/usr/bin/env python3
"""BRAINSTORM/018 Phase 5: Das x N cross-matrix collapse test.

Pre-registered discriminator (phase_05 "STAGED" step 2): plot CT-bar and the
signed Gamma-bar deficit of every {N, Das} run against the clearance group
d/sigma = N*Das/sigma. Collapse onto one curve => the Phase-2 non-monotonicity
was clearance physics and Das is (nearly) inert above the 014 floor. No
collapse => Das carries independent attachment-length physics and must be
converged separately at adequate clearance.

All runs share the B0 carrier (NT=36, sigma/c=0.68, OVERLAP 2.0 / PPS 4).
d/sigma = 0.603 * N * eta  (Das/c = 0.41*eta, sigma/c = 0.68).

Outputs data/p018_gamma_ladders/gamma_dasN_matrix.png and a numeric summary.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p018_analyze import DATA, gamma_profile, load_ct, _interp, _mean  # noqa: E402

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

BAND = (0.3, 0.95)
LO, HI = 10, 19  # matched windows, all runs settled by rev ~9

# (run, N, eta) — B0 carrier throughout
RUNS = [
    ("p018_nrows1",        1, 1.0),
    ("p018_nrows2",        2, 1.0),
    ("p018_b0",            4, 1.0),
    ("p018_das0p5",        4, 0.5),
    ("p018_das2p0",        4, 2.0),
    ("p018_das4p0",        4, 4.0),
    ("p018_nrows1_das2p0", 1, 2.0),
    ("p018_nrows2_das2p0", 2, 2.0),
    ("p018_nrows1_das4p0", 1, 4.0),
    ("p018_nrows2_das4p0", 2, 4.0),
]
REF = "p018_b0"  # global reference for the Gamma overlay


def ct_mean(run):
    rows = load_ct(run)
    vals = [ct for s, rev, ct in rows if LO <= int(rev) <= HI]
    return _mean(vals)


def main():
    ref_prof = gamma_profile(REF, load_ct(REF), LO, HI)
    ref_rs = [r for r, _ in ref_prof if BAND[0] <= r <= BAND[1]]
    gmax = max(abs(g) for r, g in ref_prof if BAND[0] <= r <= BAND[1])

    print(f"{'run':22s} {'N':>2s} {'eta':>4s} {'d/sig':>6s} {'CT_bar':>9s} "
          f"{'dCT%vsB0':>9s} {'eps_max%':>9s} {'inb(0.31)':>10s} {'outb(0.84)':>10s}")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    ct_ref = ct_mean(REF)
    cmap = {0.5: "tab:purple", 1.0: "tab:blue", 2.0: "tab:green", 4.0: "tab:red"}
    mstyle = {1: "o", 2: "s", 4: "^"}

    for run, N, eta in RUNS:
        ds = 0.603 * N * eta
        ct = ct_mean(run)
        prof = gamma_profile(run, load_ct(run), LO, HI)
        diff = [( r, ( _interp(prof, r) - _interp(ref_prof, r) ) / gmax * 100.0)
                for r in ref_rs]
        eps_max = max(abs(d) for _, d in diff)
        inb = _interp(dict(diff), 0.31) if False else [d for r, d in diff if abs(r - 0.31) < 0.02]
        inb = min((abs(r - 0.31), d) for r, d in diff)[1]
        outb = min((abs(r - 0.84), d) for r, d in diff)[1]
        print(f"{run:22s} {N:2d} {eta:4.1f} {ds:6.2f} {ct:9.6f} "
              f"{(ct/ct_ref-1)*100:+8.2f}% {eps_max:8.2f}% {inb:+9.2f}% {outb:+9.2f}%")
        ax1.plot([r for r, _ in diff], [d for _, d in diff],
                 color=cmap[eta], marker=mstyle[N], ms=3, lw=1.2,
                 label=f"N={N}, Das={0.41*eta:.2f}c (d/σ={ds:.1f})")
        ax2.scatter(ds, ct, color=cmap[eta], marker=mstyle[N], s=60, zorder=3)
        ax2.annotate(f"N{N}", (ds, ct), textcoords="offset points",
                     xytext=(5, 3), fontsize=7)

    ax1.axhline(0, color="k", lw=0.5)
    ax1.set_xlabel("r/R"); ax1.set_ylabel("ΔΓ̄ / max|Γ̄_B0|  [%]")
    ax1.set_title(f"Γ̄ difference vs B0 (revs {LO}–{HI})")
    ax1.legend(fontsize=6.5, ncol=2)
    ax2.set_xscale("log")
    ax2.set_xlabel("d/σ = N·Das/σ"); ax2.set_ylabel("CT̄")
    ax2.set_title("Collapse test: CT̄ vs d/σ (color=Das, marker=N)")
    ax2.grid(alpha=0.3)
    out = os.path.join(DATA, "p018_gamma_ladders", "gamma_dasN_matrix.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.tight_layout(); fig.savefig(out, dpi=150)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
