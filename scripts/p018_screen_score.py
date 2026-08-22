#!/usr/bin/env python3
"""BRAINSTORM/018 phase_14 screen scorer.

Scores SCREEN runs (8-rev, no freestream pulse, `scr_*`) in rung pairs at the
three stage windows, and reports the delta-vs-length stability that decides
whether a screen verdict is trustworthy. Reuses p018_analyze.py machinery.

PRE-REGISTERED EVIDENCE CLASS: screen deltas rank parameters and prune
ladders; they can never satisfy M1/M2 (>=15 settled revs on the staged
startup, HPC). Screen CT is not a CT prediction.

Usage
  p018_screen_score.py pair RUN_A RUN_B      # deltas at all stage windows
  p018_screen_score.py health RUN            # wake-health tripwire summary
  p018_screen_score.py ladder RUN1 RUN2 ...  # successive pairs down a ladder

Stage windows (integer rev blocks, rev k covers [k, k+1)): revs 1-1, 2-3, 4-7.
Rev 0 is spinup+impulse transient and never scored.
"""

import csv
import glob
import importlib.util
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
DATA = os.path.join(REPO, "data")

spec = importlib.util.spec_from_file_location("pa", os.path.join(HERE, "p018_analyze.py"))
pa = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pa)

STAGES = [(1, 1), (2, 3), (4, 7)]
BAND = (0.3, 0.95)


def ct_mean(rows, lo, hi):
    vals = [ct for _, rev, ct in rows if lo <= math.floor(rev + 1e-9) <= hi]
    if not vals:
        return float("nan")
    return sum(vals) / len(vals)


def load_rows(run):
    try:
        return pa.load_ct(run)
    except (FileNotFoundError, SystemExit):
        return pa.load_ct_from_monitor(run)


def gamma_delta(run_a, run_b, rows_a, rows_b, lo, hi):
    try:
        A = dict(pa.gamma_profile([run_a], rows_a, lo, hi))
        B = dict(pa.gamma_profile([run_b], rows_b, lo, hi))
    except SystemExit:
        return None
    gmax = max((abs(v) for r, v in A.items() if BAND[0] <= r <= BAND[1]), default=float("nan"))
    ds = [(r, 100.0 * (B[r] - A[r]) / gmax) for r in sorted(A) if r in B and BAND[0] <= r <= BAND[1]]
    if not ds:
        return None
    eps = [abs(d) for _, d in ds]
    dip = min(ds, key=lambda t: t[1])
    peak = max(ds, key=lambda t: t[1])
    return {
        "eps_max": max(eps),
        "eps_rms": math.sqrt(sum(e * e for e in eps) / len(eps)),
        "dip": dip, "peak": peak,
        "inboard": ds[0], "tip": ds[-1],
    }


def cmd_pair(run_a, run_b):
    ra, rb = load_rows(run_a), load_rows(run_b)
    print(f"--- SCREEN pair {run_a} -> {run_b}  (screen evidence class; not M1/M2)")
    print(f"{'window':>10} {'CT_A':>10} {'CT_B':>10} {'dCT%':>8} {'eps_max%':>9} "
          f"{'eps_rms%':>9}  extremes (r/R: d%)")
    for lo, hi in STAGES:
        ca, cb = ct_mean(ra, lo, hi), ct_mean(rb, lo, hi)
        dct = 100.0 * (cb - ca) / ca if ca and not math.isnan(ca) else float("nan")
        g = gamma_delta(run_a, run_b, ra, rb, lo, hi)
        ext = (f"in {g['inboard'][0]:.2f}:{g['inboard'][1]:+.1f}  "
               f"dip {g['dip'][0]:.2f}:{g['dip'][1]:+.1f}  "
               f"tip {g['tip'][0]:.2f}:{g['tip'][1]:+.1f}") if g else "n/a"
        em = f"{g['eps_max']:9.2f}" if g else f"{'n/a':>9}"
        er = f"{g['eps_rms']:9.2f}" if g else f"{'n/a':>9}"
        print(f"  revs {lo}-{hi:>2} {ca:10.6f} {cb:10.6f} {dct:8.2f} {em} {er}  {ext}")
    print("  stability: trust the pair only if dCT and the Gamma shape are "
          "consistent across the last two windows.")


def cmd_health(run):
    hits = []
    for f in glob.glob(os.path.join(DATA, run, "monitors", "*.csv")):
        with open(f) as fh:
            if "min_sigma_ratio" in fh.readline():
                hits.append(f)
    if not hits:
        print(f"{run}: no wake-health CSV found")
        return
    rows = list(csv.DictReader(open(hits[0])))
    if not rows:
        print(f"{run}: wake-health CSV empty")
        return
    def col(k):
        return [float(r[k]) for r in rows if r.get(k) not in (None, "", "NaN")]
    n = col("n_particles"); mu = col("max_u"); msr = col("min_sigma_ratio")
    last = rows[-1]
    print(f"--- wake health {run} ({len(rows)} steps)")
    print(f"  n_particles: final {int(n[-1])}, max {int(max(n))}")
    print(f"  max_u: final {mu[-1]:.1f}, worst {max(mu):.1f} m/s")
    print(f"  min_sigma_ratio: final {msr[-1]:.4f}, worst {min(msr):.4f}")
    if len(msr) > 100:
        tail = msr[-100:]
        slope = (tail[-1] - tail[0]) / 100.0
        print(f"  contraction rate (last 100 steps): {slope:+.2e} /step "
              f"(negative = still contracting)")
    print(f"  wall_s: last {last.get('wall_s', 'n/a')}")


def cmd_ladder(runs):
    for a, b in zip(runs, runs[1:]):
        cmd_pair(a, b)
        print()


def main():
    if len(sys.argv) < 3:
        print(__doc__); sys.exit(2)
    cmd, args = sys.argv[1], sys.argv[2:]
    if cmd == "pair" and len(args) == 2:
        cmd_pair(*args)
    elif cmd == "health" and len(args) == 1:
        cmd_health(args[0])
    elif cmd == "ladder" and len(args) >= 2:
        cmd_ladder(args)
    else:
        print(__doc__); sys.exit(2)


if __name__ == "__main__":
    main()
