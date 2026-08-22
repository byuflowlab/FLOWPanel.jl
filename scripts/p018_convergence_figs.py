#!/usr/bin/env python3
"""BRAINSTORM/018: emit CSV data for the living notebook convergence figures.

One figure per convergence axis, two panels each: CT-bar vs parameter (left)
and Gamma-bar(r/R) per rung (right). This script computes the numbers with the
same machinery as the campaign scoring (`p018_analyze.py`: m1 windows,
gamma_profile, force-monitor fallback) and writes plain CSVs into the notebook
image directory; the standalone TikZ figures there read them.

    FIGDIR/fig_<axis>/ct.csv       series,label,param,ct,eminus,eplus
    FIGDIR/fig_<axis>/flags.csv    series,label,param,status   (ignited/pending/diverged)
    FIGDIR/fig_<axis>/gamma_<slug>.csv   r,gamma

Update flow when a new rung lands: harvest the run locally, add/flip its entry
in AXES below, rerun this script, then `./build.sh` in FIGDIR. The notebook's
image references never change.

Windows are the ledger's matched windows (see 018 ledger.md); anchors this
script must reproduce: b0@15-29 = 0.071701, dasc ladder = 0.071571 / 0.072154 /
0.073455, nt144@10-19 = 0.072023, L1@17-31 = 0.072669.
"""

import csv
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p018_analyze import (DATA, load_ct, gamma_profile, per_rev, _mean,  # noqa: E402
                          moving_block_ci)

FIGDIR = os.path.expanduser(
    "~/Dropbox/research/notebooks/img/20260813_018_convergence")


def chain(runs):
    """Stitch a chain of restart segments on the step column.

    Boundary is each segment's own first non-placeholder step, so no
    restart-step bookkeeping is needed (placeholder rows are already dropped
    by load_ct).
    """
    rows = load_ct(runs[0])
    for seg in runs[1:]:
        b = load_ct(seg)
        first = min(r[0] for r in b)
        rows = [r for r in rows if r[0] < first] + b
    return sorted(rows, key=lambda r: r[0])


def ct_stats(rows, lo, hi):
    bins = per_rev(rows)
    revs = sorted(r for r in bins if lo <= r <= hi)
    if not revs:
        raise SystemExit(f"no samples in revs [{lo},{hi}]")
    rev_means = [_mean(bins[r]) for r in revs]
    ctbar = _mean(rev_means)
    lo_ci, hi_ci = moving_block_ci(rev_means)
    return ctbar, lo_ci, hi_ci


def slug(s):
    return re.sub(r"[^a-z0-9]+", "_", s.lower()).strip("_")


# A CT rung: (series, label, param, runs_chain, lo, hi)
# A flag:    (series, label, param, status)          -- no CT value plotted
# A gamma:   (label, monitor_runs, rows_chain, lo, hi)
#   monitor_runs may differ from the CT chain when a segment lacks the
#   bound-circulation monitor (e.g. p018_nt144 parent).
AXES = {
    "das": dict(
        rungs=[
            ("eta",   "0.205c",        0.205, ["p018_das0p5"],  10, 19),
            ("eta",   "0.41c (B0)",    0.41,  ["p018_b0"],      10, 19),
            ("eta",   "0.82c",         0.82,  ["p018_das2p0"],  10, 19),
            ("eta",   "1.64c",         1.64,  ["p018_das4p0"],  10, 19),
            ("kappa", "k=0.25",        0.25,  ["p018_dasc0p25"], 10, 19),
            ("kappa", "k=0.41",        0.41,  ["p018_dasc0p41"], 10, 19),
            ("kappa", "k=0.82",        0.82,  ["p018_dasc0p82"], 10, 19),
            # ufront law (phase_13 par.7): Das = D*sigma span-uniform, N=1,
            # viscous, OV 2.75/PPS 12 -- a DIFFERENT carrier from B0, and
            # sigma varies down the series (0.035/0.04/0.05R), so along
            # Das/c this is not a pure Das ladder. param = D*sigma / c(0.75R)
            # with c(0.75R) = 0.1277R (validated: L1 sigma/c 0.297 -> 0.0379R
            # = 019's sigma* 0.0381R). Windows are the ledger's settled
            # windows (30-44 for s040 via _s2, 15-29 for s050), not 10-19.
            ("ufront", "3.4s @ s=0.04R", 1.065,
             ["p018_ufront_n1_s040_visc", "p018_ufront_n1_s040_visc_s2"],
             30, 44),
            ("ufront", "3.4s @ s=0.05R", 1.334,
             ["p018_ufront_n1_s050_visc"], 15, 29),
        ],
        flags=[("ufront", "3.0s @ s=0.035R", 0.820, "ignited")],
        gamma=[
            ("Das 0.205c",    ["p018_das0p5"],   10, 19),
            ("Das 0.41c B0",  ["p018_b0"],       10, 19),
            ("Das 1.64c",     ["p018_das4p0"],   10, 19),
            ("kappa 0.25",    ["p018_dasc0p25"], 10, 19),
            ("kappa 0.82",    ["p018_dasc0p82"], 10, 19),
            ("ufront 3.4s s0.04R",
             ["p018_ufront_n1_s040_visc", "p018_ufront_n1_s040_visc_s2"],
             30, 44),
        ],
    ),
    "dt": dict(
        rungs=[
            ("dt", "NT=36",  36,  ["p018_b0"],        10, 19),
            ("dt", "NT=72",  72,  ["p018_nt72_eta2"], 10, 19),
            ("dt", "NT=144", 144, ["p018_nt144", "p018_nt144_s2"], 10, 19),
        ],
        flags=[],
        gamma=[
            ("NT=36",  ["p018_b0"],        10, 19),
            ("NT=72",  ["p018_nt72_eta2"], 10, 19),
            # parent nt144 has no bound-circulation monitor locally; the _s2
            # segment covers revs ~14-19 of the window (partial, disclosed).
            ("NT=144 (revs 14-19)", ["p018_nt144_s2"], 10, 19),
        ],
        gamma_rows={"NT=144 (revs 14-19)": ["p018_nt144", "p018_nt144_s2"]},
    ),
    "sigma": dict(
        rungs=[
            ("ladderA", "0.68 (B0)", 0.68,  ["p018_b0", "p018_b0_s2"], 15, 29),
            ("ladderA", "0.297 (L1)", 0.297, ["p018_L1", "p018_L1_s2"], 17, 31),
            ("ladderB", "0.68 (B0)", 0.68,  ["p018_b0", "p018_b0_s2"], 15, 29),
            ("ladderB", "0.477",     0.477, ["p018_ov1p4"], 15, 29),
            ("ladderB", "0.34",      0.34,  ["p018_ov1p0"], 15, 29),
        ],
        flags=[("ladderA", "L2", 0.151, "diverged")],
        gamma=[
            ("s/c 0.68 B0", ["p018_b0", "p018_b0_s2"], 15, 29),
            ("s/c 0.297 L1", ["p018_L1", "p018_L1_s2"], 17, 31),
            ("s/c 0.477 fixed h", ["p018_ov1p4"], 15, 29),
            ("s/c 0.34 fixed h", ["p018_ov1p0"], 15, 29),
        ],
    ),
    "h": dict(
        rungs=[
            ("h", "0.5 (B0)", 0.5,   ["p018_b0", "p018_b0_s2"], 15, 29),
            ("h", "0.25",     0.25,  ["p018_h0p25"], 15, 29),
            ("h", "0.125",    0.125, ["p018_h0p125", "p018_h0p125_s2"], 15, 29),
        ],
        flags=[],
        gamma=[
            ("h/s 0.5 B0", ["p018_b0", "p018_b0_s2"], 15, 29),
            ("h/s 0.25",   ["p018_h0p25"], 15, 29),
            ("h/s 0.125",  ["p018_h0p125", "p018_h0p125_s2"], 15, 29),
        ],
    ),
    "nrows": dict(
        rungs=[
            ("N", "N=1", 1, ["p018_nrows1"], 10, 19),
            ("N", "N=2", 2, ["p018_nrows2"], 10, 19),
            ("N", "N=4 (B0)", 4, ["p018_b0"], 10, 19),
        ],
        flags=[],
        gamma=[
            ("N=1", ["p018_nrows1"], 10, 19),
            ("N=2", ["p018_nrows2"], 10, 19),
            ("N=4 B0", ["p018_b0"], 10, 19),
        ],
    ),
    "rlxf": dict(
        # Matched revs 45-57: the _s3 chain segment is the stock point; the
        # upward arms warm-started from the same _s2@1619 state (2026-08-14
        # harvest). Below stock both reduced rungs ignited (flags).
        rungs=[
            ("rlxf", "0.3 (stock)", 0.3,
             ["p018_ufront_n1_s040_visc_s3"], 45, 57),
            ("rlxf", "0.45", 0.45, ["p018_rlxf0p45_n1_s040_visc"], 45, 57),
            ("rlxf", "0.6",  0.6,  ["p018_rlxf0p6_n1_s040_visc"],  45, 57),
        ],
        flags=[
            ("rlxf", "0.075", 0.075, "ignited"),
            ("rlxf", "0.15",  0.15,  "ignited"),
        ],
        gamma=[
            ("rlxf 0.3 (stock)", ["p018_ufront_n1_s040_visc_s3"], 45, 57),
            ("rlxf 0.45", ["p018_rlxf0p45_n1_s040_visc"], 45, 57),
            ("rlxf 0.6",  ["p018_rlxf0p6_n1_s040_visc"],  45, 57),
        ],
    ),
}


def main():
    for axis, spec in AXES.items():
        d = os.path.join(FIGDIR, f"fig_{axis}")
        os.makedirs(d, exist_ok=True)

        with open(os.path.join(d, "ct.csv"), "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["series", "label", "param", "ct", "eminus", "eplus"])
            for series, label, param, runs, lo, hi in spec["rungs"]:
                rows = chain(runs)
                ct, lo_ci, hi_ci = ct_stats(rows, lo, hi)
                w.writerow([series, label, param, f"{ct:.6f}",
                            f"{ct - lo_ci:.6f}", f"{hi_ci - ct:.6f}"])
                print(f"[{axis}] {series:8s} {label:14s} param={param:<7g} "
                      f"CT={ct:.6f} [{lo_ci:.6f},{hi_ci:.6f}] revs {lo}-{hi}")

        with open(os.path.join(d, "flags.csv"), "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["series", "label", "param", "status"])
            for series, label, param, status in spec["flags"]:
                w.writerow([series, label, param, status])

        gamma_rows_override = spec.get("gamma_rows", {})
        for label, runs, lo, hi in spec["gamma"]:
            rows = chain(gamma_rows_override.get(label, runs))
            prof = gamma_profile(runs, rows, lo, hi)
            path = os.path.join(d, f"gamma_{slug(label)}.csv")
            with open(path, "w", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(["r", "gamma"])
                for r, g in prof:
                    w.writerow([f"{r:.6f}", f"{g:.8f}"])
            print(f"[{axis}] gamma -> {os.path.basename(path)} "
                  f"({len(prof)} sections)")


if __name__ == "__main__":
    sys.exit(main())
