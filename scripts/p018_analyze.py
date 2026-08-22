#!/usr/bin/env python3
"""BRAINSTORM/018 harvest analysis: M1 (settled-mean CT) and M2 (circulation).

Implements the metrics defined in
BRAINSTORM/018_dji9443_hover_convergence_campaign/decision_rules.md.

M1  CT-bar = cycle-mean of CT_bernoulli over a rev window, with 5-rev block-mean
    drift and a moving-block bootstrap 95% CI on per-rev means.
M2  Gamma-bar(r/R) = per-section cycle-mean of circulation_te (TE mu jump),
    blades averaged; eps_Gamma = max/RMS of |dGamma|/max_r|Gamma| on
    0.3 <= r/R <= 0.95.

Usage
  p018_analyze.py m1 RUN [--revs LO HI] [--stitch SEG --restart-step S]
  p018_analyze.py m2 RUN_A RUN_B --revs-a LO HI --revs-b LO HI
  p018_analyze.py windows RUN            # per-rev table, to pick a window

Runs are directory names under data/ (e.g. p018_b0). Rev windows are inclusive
integer revolution indices, matching the `rev_block` convention in the
*_CT_per_rev.csv files (rev index k covers revolution in [k, k+1)).
"""

import argparse
import csv
import math
import os
import random
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(REPO, "data")

# Placeholder rows written by a restart segment before its first real step have
# CT identically +/-0.0; they are not physical samples.
PLACEHOLDER_TOL = 1e-12

# Rotor speed used to convert monitor time to revolutions in the fallback path.
RPM_DEFAULT = 5400.0


def _ct_path(run):
    return os.path.join(DATA, run, f"{run}_CT_vs_rev.csv")


def _force_monitor_path(run):
    return os.path.join(DATA, run, "monitors",
                        f"{run}_monitor02_force_system1.csv")


def load_ct_from_monitor(run, rpm=RPM_DEFAULT):
    """Rebuild the CT series from the incrementally-written force monitor.

    `<run>_CT_vs_rev.csv` is written only after the time loop finishes
    (examples/rotor_hover_pressure_comparison.jl), so a walltime-killed or
    diverged segment leaves none. The ForceMonitor CSV survives, and maps
    exactly onto it (verified bit-for-bit on p018_dasc0p41, 720 rows):

        CT_bernoulli = -CFx ,  step_csv = step_monitor + 1 ,
        revolution   = time * rpm / 60 .

    Only CT_bernoulli is recoverable this way (the Laplace/KJ columns are not
    in the force monitor); m1/windows/m2 use only CT_bernoulli.
    """
    rows = []
    with open(_force_monitor_path(run)) as fh:
        for rec in csv.DictReader(fh):
            try:
                ct = -float(rec["CFx"])
                t = float(rec["time"])
            except (TypeError, ValueError, KeyError):
                continue
            if not math.isfinite(ct) or abs(ct) < PLACEHOLDER_TOL:
                continue
            rows.append((int(rec["step"]) + 1, t * rpm / 60.0, ct))
    return rows


def load_ct(run):
    """Return [(step, revolution, CT_bernoulli)] with placeholder rows dropped.

    Falls back to the force monitor when the end-of-run CT CSV is absent.
    """
    if not os.path.exists(_ct_path(run)) and os.path.exists(_force_monitor_path(run)):
        sys.stderr.write(
            f"note: {run} has no CT_vs_rev.csv (run did not finish); "
            "reconstructing CT_bernoulli from the force monitor\n")
        return load_ct_from_monitor(run)
    rows = []
    with open(_ct_path(run)) as fh:
        for rec in csv.DictReader(fh):
            try:
                ct = float(rec["CT_bernoulli"])
            except (TypeError, ValueError):
                continue
            if not math.isfinite(ct) or abs(ct) < PLACEHOLDER_TOL:
                continue
            rows.append((int(rec["step"]), float(rec["revolution"]), ct))
    return rows


def stitch(parent, segment, restart_step):
    """Concatenate a parent run and its restart segment on the `step` column.

    The segment restarts from VTU step `restart_step`, so parent rows with
    step > restart_step are superseded by the segment's own integration.
    """
    a = [r for r in load_ct(parent) if r[0] <= restart_step]
    b = [r for r in load_ct(segment) if r[0] > restart_step]
    if not b:
        raise SystemExit(f"segment {segment} has no rows past step {restart_step}")
    return sorted(a + b, key=lambda r: r[0])


def per_rev(rows):
    """Group samples into integer-revolution bins -> {rev: [ct, ...]}."""
    bins = {}
    for _, rev, ct in rows:
        bins.setdefault(int(math.floor(rev + 1e-9)), []).append(ct)
    return bins


def _mean(xs):
    return sum(xs) / len(xs)


def _std(xs):
    if len(xs) < 2:
        return 0.0
    m = _mean(xs)
    return math.sqrt(sum((x - m) ** 2 for x in xs) / (len(xs) - 1))


def moving_block_ci(vals, block=3, n_boot=4000, seed=20260803):
    """95% moving-block bootstrap CI on the mean of a serially-correlated series."""
    n = len(vals)
    if n < 2:
        return (float("nan"), float("nan"))
    block = max(1, min(block, n))
    n_blocks = int(math.ceil(n / block))
    starts = list(range(n - block + 1))
    rng = random.Random(seed)
    means = []
    for _ in range(n_boot):
        sample = []
        for _ in range(n_blocks):
            s = rng.choice(starts)
            sample.extend(vals[s:s + block])
        means.append(_mean(sample[:n]))
    means.sort()
    return (means[int(0.025 * n_boot)], means[int(0.975 * n_boot)])


def m1(rows, lo, hi, label=""):
    """Settled-mean CT over revs [lo, hi] with drift and CI diagnostics."""
    bins = per_rev(rows)
    revs = sorted(r for r in bins if lo <= r <= hi)
    if not revs:
        raise SystemExit(f"no samples in revs [{lo},{hi}]")
    rev_means = [_mean(bins[r]) for r in revs]
    ctbar = _mean(rev_means)
    lo_ci, hi_ci = moving_block_ci(rev_means)

    # 5-rev block means: settledness requires drift < ~0.3% and no monotone trend
    blocks = []
    for i in range(0, len(revs) - 4, 5):
        chunk = rev_means[i:i + 5]
        if len(chunk) == 5:
            blocks.append((revs[i], revs[i + 4], _mean(chunk)))
    drift = 0.0
    monotone = False
    if len(blocks) >= 2:
        bm = [b[2] for b in blocks]
        drift = (max(bm) - min(bm)) / ctbar
        monotone = all(y > x for x, y in zip(bm, bm[1:])) or \
                   all(y < x for x, y in zip(bm, bm[1:]))

    print(f"--- M1 {label}")
    print(f"  window          revs {revs[0]}-{revs[-1]}  ({len(revs)} revs, "
          f"{sum(len(bins[r]) for r in revs)} steps)")
    print(f"  CT_bar          {ctbar:.6f}")
    print(f"  95% CI          [{lo_ci:.6f}, {hi_ci:.6f}]")
    print(f"  per-rev std     {_std(rev_means):.6f}")
    print(f"  5-rev blocks    " + ", ".join(f"{a}-{b}:{m:.6f}" for a, b, m in blocks))
    print(f"  block drift     {100*drift:.3f}%   monotone={monotone}")
    ok = (drift < 0.003 and not monotone and len(revs) >= 15)
    print(f"  M1 settled      {'PASS' if ok else 'CHECK'}"
          f"{'' if len(revs) >= 15 else '  (< 15 revs)'}")
    return ctbar, (lo_ci, hi_ci), revs


def _monitor_path(run):
    d = os.path.join(DATA, run, "monitors")
    hits = [f for f in os.listdir(d) if "bound_circulation" in f and f.endswith(".csv")]
    if not hits:
        raise SystemExit(f"no bound_circulation monitor in {d}")
    return os.path.join(d, sorted(hits)[0])


def gamma_profile(runs, rows, lo, hi):
    """Blade-averaged, window-averaged Gamma(r/R) from circulation_te.

    `runs` may be several run names (a parent plus its restart segments); each
    monitor is read and only steps inside the window contribute, so a stitched
    history is covered by passing every segment that holds part of it.
    """
    if isinstance(runs, str):
        runs = [runs]
    steps = {s for s, rev, _ in rows if lo <= math.floor(rev + 1e-9) <= hi}
    acc = {}
    for run in runs:
        with open(_monitor_path(run)) as fh:
            for rec in csv.DictReader(fh):
                if int(rec["step"]) not in steps:
                    continue
                r = abs(float(rec["r_over_R"]))
                g = float(rec["circulation_te"])
                if not math.isfinite(g):
                    continue
                key = round(r, 6)
                a = acc.setdefault(key, [0.0, 0])
                a[0] += g
                a[1] += 1
    if not acc:
        raise SystemExit(f"{runs}: no circulation samples in revs [{lo},{hi}]")
    return sorted((r, v[0] / v[1]) for r, v in acc.items())


def _interp(prof, x):
    rs = [p[0] for p in prof]
    gs = [p[1] for p in prof]
    if x <= rs[0]:
        return gs[0]
    if x >= rs[-1]:
        return gs[-1]
    for i in range(len(rs) - 1):
        if rs[i] <= x <= rs[i + 1]:
            t = (x - rs[i]) / (rs[i + 1] - rs[i])
            return gs[i] + t * (gs[i + 1] - gs[i])
    return gs[-1]


def m2(prof_a, prof_b, label=""):
    """eps_Gamma between two profiles on 0.3 <= r/R <= 0.95, A's grid."""
    grid = [r for r, _ in prof_a if 0.3 <= r <= 0.95]
    if not grid:
        raise SystemExit("no sections in 0.3 <= r/R <= 0.95")
    scale = max(abs(_interp(prof_a, r)) for r in grid)
    diffs = [abs(_interp(prof_a, r) - _interp(prof_b, r)) / scale for r in grid]
    emax, erms = max(diffs), math.sqrt(_mean([d * d for d in diffs]))
    print(f"--- M2 {label}")
    print(f"  sections        {len(grid)} on 0.3<=r/R<=0.95   max|Gamma_A|={scale:.6f}")
    print(f"  eps_max         {100*emax:.3f}%")
    print(f"  eps_rms         {100*erms:.3f}%")
    print(f"  M2 (<=1% max)   {'PASS' if emax <= 0.01 else 'FAIL'}")
    return emax, erms


def main():
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd", required=True)

    q = sub.add_parser("m1")
    q.add_argument("run")
    q.add_argument("--revs", nargs=2, type=int)
    q.add_argument("--stitch")
    q.add_argument("--restart-step", type=int)

    q = sub.add_parser("m2")
    q.add_argument("run_a")
    q.add_argument("run_b")
    q.add_argument("--revs-a", nargs=2, type=int, required=True)
    q.add_argument("--revs-b", nargs=2, type=int, required=True)
    q.add_argument("--stitch-a")
    q.add_argument("--stitch-b")
    q.add_argument("--restart-step", type=int, default=719)

    q = sub.add_parser("windows")
    q.add_argument("run")
    q.add_argument("--stitch")
    q.add_argument("--restart-step", type=int)

    a = p.parse_args()

    if a.cmd in ("m1", "windows"):
        if a.stitch:
            if a.restart_step is None:
                raise SystemExit("--stitch requires --restart-step")
            rows = stitch(a.run, a.stitch, a.restart_step)
            label = f"{a.run}+{a.stitch}"
        else:
            rows = load_ct(a.run)
            label = a.run
        if a.cmd == "windows":
            bins = per_rev(rows)
            print(f"rev  n     CT_mean     ptp        [{label}]")
            for r in sorted(bins):
                v = bins[r]
                print(f"{r:3d}  {len(v):3d}   {_mean(v):.6f}   {max(v)-min(v):.6f}")
            return
        lo, hi = a.revs if a.revs else (min(per_rev(rows)), max(per_rev(rows)))
        m1(rows, lo, hi, label)
        return

    if a.cmd == "m2":
        def arm(run, seg):
            if seg:
                return stitch(run, seg, a.restart_step), [run, seg]
            return load_ct(run), [run]
        ra, na = arm(a.run_a, a.stitch_a)
        rb, nb = arm(a.run_b, a.stitch_b)
        pa = gamma_profile(na, ra, *a.revs_a)
        pb = gamma_profile(nb, rb, *a.revs_b)
        m2(pa, pb, f"{'+'.join(na)} vs {'+'.join(nb)}")
        return


if __name__ == "__main__":
    sys.exit(main())
