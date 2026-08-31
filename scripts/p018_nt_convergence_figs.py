#!/usr/bin/env python3
"""Emit backing CSVs for the 2026-08-25 018 NT-ladder convergence figures.

Three figures are backed:

  fig_ct_history   CT vs revolution, one panel per NT rung, restart seams marked
  fig_ct_rungs     CTbar vs NT and vs lambda, matched and best-available windows
  fig_gamma_span   Gammabar(r/R) and Delta% for the lambda ladder at NT72
  fig_gamma_nt     Gammabar(r/R) and Delta% for the NT ladder at lambda=3.4
  fig_loading_span spanwise loading dCT/dx and Delta, lambda ladder at NT72
  fig_loading_nt   spanwise loading dCT/dx and Delta, NT ladder at lambda=3.4

This is a sibling of scripts/p018_convergence_figs.py, deliberately NOT an
extension of it: that script has no CLI and regenerates all six published
20260813_018_convergence axes on every invocation, and its AXES schema
(one scalar CTbar per rung) cannot express a time series, a seam record, or
two windows for the same rung.

Why the windows are what they are.  The arms have very different histories
(NT36 = 30 rev, NT72 ~ 25-26.5 rev, NT144 ~ 13.3-13.8 rev), so a naive
"last 5 revs" comparison across rungs conflates settling with the timestep
effect.  Every rung is therefore emitted twice:

  matched          revs 8-12, five complete bins reachable by every arm
  best-available   the last five complete rev bins that arm actually has

Five bins in both so the moving-block CI widths are directly comparable.
For NT144 the two windows coincide -- that is not a bug, it is the finding:
NT144 has no settled window to offer.

Usage:
    python3 scripts/p018_nt_convergence_figs.py check
    python3 scripts/p018_nt_convergence_figs.py emit [--figdir DIR] [--only K]
"""

import argparse
import csv
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p018_analyze import (DATA, load_ct, per_rev, _mean, _interp, _monitor_path,
                          gamma_profile)
from p018_convergence_figs import chain, ct_stats, slug

FIGDIR_DEFAULT = os.path.expanduser(
    "~/Dropbox/research/notebooks/img/20260825_p018_nt_convergence")

MATCHED = (8, 12)          # revs, inclusive; reachable by every arm
NWIN = 5                   # bins in the best-available window
SEAM_GATE = 0.05           # percent, ops_reference.md
M2_GATE = 1.0              # percent, decision_rules.md
BAND = (0.3, 0.95)         # r/R band the M2 metric is evaluated on
BLOWUP_CT = 0.25           # |CT| above this is a divergence, not a sample

# label, slug, NT, lambda, N, chain, model definition, colour key
ARMS = [
    # --- NT36 csarc, exact-rate rlxf 0.3 (the base rung) ---
    ("NT36 $\\lambda$2.4",  "nt36_l2p4",  36, 2.4, 1, ["p018_csarc_l2p4"],       "exact", "lamA"),
    ("NT36 $\\lambda$3.4",  "nt36_l3p4",  36, 3.4, 1, ["p018_csarc_l3p4"],       "exact", "lamB"),
    ("NT36 $\\lambda$4.8",  "nt36_l4p8",  36, 4.8, 1, ["p018_csarc_l4p8"],       "exact", "lamC"),
    ("NT36 $\\lambda$4.8 N0", "nt36_n0_l4p8", 36, 4.8, 0, ["p018_csarc_n0_l4p8"], "exact", "lamN0"),
    # --- NT72 csarc, exact-rate rlxf 0.16334 ---
    ("NT72 $\\lambda$2.4",  "nt72_l2p4",  72, 2.4, 1,
     ["p018_csarc_nt72_l2p4", "p018_csarc_nt72_l2p4_s2"], "exact", "lamA"),
    ("NT72 $\\lambda$3.4",  "nt72_l3p4",  72, 3.4, 1,
     ["p018_csarc_nt72_l3p4", "p018_csarc_nt72_l3p4_s2"], "exact", "lamB"),
    ("NT72 $\\lambda$4.8",  "nt72_l4p8",  72, 4.8, 1,
     ["p018_csarc_nt72_l4p8", "p018_csarc_nt72_l4p8_s2"], "exact", "lamC"),
    ("NT72 $\\lambda$4.8 N0", "nt72_n0_l4p8", 72, 4.8, 0,
     ["p018_csarc_n0_nt72_l4p8", "p018_csarc_n0_nt72_l4p8_s2"], "exact", "lamN0"),
    # --- NT144 csarc, exact-rate rlxf 0.08539 (UNSETTLED) ---
    ("NT144 $\\lambda$2.4", "nt144_l2p4", 144, 2.4, 1,
     ["p018_csarc_nt144_l2p4", "p018_csarc_nt144_l2p4_s2"], "exact", "lamA"),
    ("NT144 $\\lambda$3.4", "nt144_l3p4", 144, 3.4, 1,
     ["p018_csarc_nt144_l3p4", "p018_csarc_nt144_l3p4_s2"], "exact", "lamB"),
    ("NT144 $\\lambda$4.8", "nt144_l4p8", 144, 4.8, 1,
     ["p018_csarc_nt144_l4p8", "p018_csarc_nt144_l4p8_s2"], "exact", "lamC"),
    ("NT144 $\\lambda$4.8 N0", "nt144_n0_l4p8", 144, 4.8, 0,
     ["p018_csarc_n0_nt144_l4p8", "p018_csarc_n0_nt144_l4p8_s2"], "exact", "lamN0"),
]

# Prior / reference arms.  mid_l3p4 is an independent NT36 lambda=3.4 run -- it
# matters right now because nt72_l3p4_s2 is the arm whose seam FAILED.
PRIORS = [
    ("NT36 $\\lambda$3.4 (midpoint A/B)", "prior_mid_l3p4", 36, 3.4, 1,
     ["p018_csarc_mid_l3p4"], "exact", "prior"),
    ("B0 carrier", "prior_b0", 36, None, 1, ["p018_b0", "p018_b0_s2"], "exact", "prior"),
]

# Legacy arms: NAIVE linear rate rlxf = 0.3*36/NT, NOT the exact-rate rule.
# Verified from metadata: upin_nt72 rlxf=0.15, upin_nt144/nt144_s3 rlxf=0.075,
# versus the exact-rate 0.16334 / 0.08539 -- 8.9% and 12.2% lower.  That rate
# sits near the ignition boundary and p018_upin_nt72 did in fact blow up.
# These are a DIFFERENT MODEL DEFINITION and must never be joined by a trend
# line to the csarc rungs.
LEGACY = [
    ("NT144 legacy rlxf 0.075", "legacy_nt144", 144, None, 1,
     ["p018_nt144", "p018_nt144_s2", "p018_nt144_s3"], "naive", "legacy"),
    ("NT72 legacy rlxf 0.15",  "legacy_upin_nt72", 72, None, 1,
     ["p018_upin_nt72"], "naive", "legacy"),
    ("NT144 legacy upin rlxf 0.075", "legacy_upin_nt144", 144, None, 1,
     ["p018_upin_nt144"], "naive", "legacy"),
]

RLXF = {("exact", 36): 0.3, ("exact", 72): 0.16334, ("exact", 144): 0.08539,
        ("naive", 72): 0.15, ("naive", 144): 0.075}

# Anchors from the 2026-08-25 harvest (ledger.md).  `check` must reproduce
# these before anything is written.
ANCHOR_CT = {            # chain -> CTbar over revs 20-25
    "nt72_l2p4": 0.072386,
    "nt72_l3p4": 0.072900,
    "nt72_l4p8": 0.074415,
    "nt72_n0_l4p8": 0.075921,
}
# NOTE: the ledger's 2026-08-25 harvest recorded 2.631% and 3.529%.  Those are
# artifacts of a latent bug in `p018_analyze.py m2`: its `--restart-step`
# (:275) is a SINGLE value applied to BOTH runs via --stitch-a/--stitch-b, but
# the four NT72 arms restarted at DIFFERENT steps (1453/1459/1518/1449).
# Passing arm A's boundary for arm B makes stitch() keep B's parent rows only
# up to A's step while B's segment starts later, punching a hole in B's history
# -- 6 steps for the l2p4/l3p4 pair (shared 1453), 59 steps for l3p4/l4p8
# (shared 1459).  Reproduced exactly: shared 1453 -> 2.631%, shared 1459 ->
# 2.654%; shared 1459 -> 3.529%, shared 1518 -> 3.397%.  chain() has no hole
# (it keeps parent rows strictly below the segment's own first step), so these
# are the correct values.  The M2 verdict is unchanged either way: both still
# FAIL the 1% gate by ~3x.
ANCHOR_M2 = {("nt72_l2p4", "nt72_l3p4"): 2.620,
             ("nt72_l3p4", "nt72_l4p8"): 3.429}
ANCHOR_SEAM = {"nt72_l2p4": 0.0062, "nt72_l3p4": 0.1415,
               "nt72_l4p8": 0.0347, "nt72_n0_l4p8": 0.0012}


# --------------------------------------------------------------------------
# loading

def healthy(rows):
    """Truncate a chain at its first divergent sample.

    p018_upin_nt72 ends at CFz ~ -5.3e5; plotting that compresses every real
    curve to a flat line.  Returns (rows, diverged_at_rev or None).
    """
    for i, (_, rev, ct) in enumerate(rows):
        if not math.isfinite(ct) or abs(ct) > BLOWUP_CT:
            return rows[:i], rev
    return rows, None


def complete_bins(rows, nt, frac=0.9):
    """Rev bins holding at least `frac` of a full revolution of samples."""
    bins = per_rev(rows)
    return sorted(r for r, v in bins.items() if len(v) >= frac * nt)


def last_window(rows, nt, n=NWIN):
    good = complete_bins(rows, nt)
    if len(good) < n:
        return (good[0], good[-1]) if good else None
    tail = good[-n:]
    return (tail[0], tail[-1])


def seam_of(arm_chain, rows):
    """CT jump across each restart boundary, per ops_reference.md:110."""
    out = []
    for seg in arm_chain[1:]:
        seg_rows = load_ct(seg)
        if not seg_rows:
            continue
        first = min(r[0] for r in seg_rows)
        before = [r for r in rows if r[0] < first]
        after = [r for r in seg_rows if r[0] == first]
        if not before or not after:
            continue
        ct_b, ct_a = before[-1][2], after[0][2]
        jump = 100.0 * abs(ct_a - ct_b) / abs(ct_b)
        out.append(dict(step=first, rev=after[0][1], ct_before=ct_b,
                        ct_after=ct_a, jump_pct=jump,
                        gate="fail" if jump > SEAM_GATE else "pass"))
    return out


def gamma_col(runs, rows, lo, hi, col):
    """gamma_profile() with a selectable circulation column.

    p018_analyze.gamma_profile hardcodes `circulation_te`; the M2 gate is
    defined on it.  Emitting `circulation_slice` as well separates a
    force-reconstruction artifact (TE only) from a genuine circulation
    difference (both).  Identical binning: abs() the signed r_over_R, average
    the two blades by rounding r to 6 dp.
    """
    if isinstance(runs, str):
        runs = [runs]
    steps = {s for s, rev, _ in rows if lo <= math.floor(rev + 1e-9) <= hi}
    acc = {}
    for run in runs:
        try:
            path = _monitor_path(run)
        except SystemExit:
            continue
        with open(path) as fh:
            for rec in csv.DictReader(fh):
                if int(rec["step"]) not in steps:
                    continue
                try:
                    g = float(rec[col])
                except (KeyError, ValueError):
                    continue
                if not math.isfinite(g):
                    continue
                key = round(abs(float(rec["r_over_R"])), 6)
                a = acc.setdefault(key, [0.0, 0])
                a[0] += g
                a[1] += 1
    if not acc:
        return None
    return sorted((r, v[0] / v[1]) for r, v in acc.items())


def eps_gamma(prof_a, prof_b):
    """eps_Gamma between two profiles on BAND, on A's grid (mirrors m2())."""
    grid = [r for r, _ in prof_a if BAND[0] <= r <= BAND[1]]
    if not grid:
        return None
    scale = max(abs(_interp(prof_a, r)) for r in grid)
    diffs = [(r, (_interp(prof_b, r) - _interp(prof_a, r)) / scale) for r in grid]
    emax = max(abs(d) for _, d in diffs)
    erms = math.sqrt(_mean([d * d for _, d in diffs]))
    worst = max(diffs, key=lambda t: abs(t[1]))[0]
    return dict(eps_max_pct=100 * emax, eps_rms_pct=100 * erms, r_worst=worst,
                curve=[(r, 100 * d) for r, d in diffs],
                verdict="PASS" if 100 * emax <= M2_GATE else "FAIL")


def build():
    """Load every arm once; returns a list of records."""
    recs = []
    tagged = ([(a, "csarc") for a in ARMS] + [(a, "prior") for a in PRIORS] +
              [(a, "legacy") for a in LEGACY])
    for (label, sl, nt, lam, n, runs, mdef, ckey), family in tagged:
        missing = [r for r in runs if not os.path.isdir(os.path.join(DATA, r))]
        if missing:
            print(f"  SKIP {sl}: missing {missing}")
            continue
        try:
            rows = chain(runs)
        except SystemExit as e:
            print(f"  SKIP {sl}: {e}")
            continue
        rows, div = healthy(rows)
        if not rows:
            print(f"  SKIP {sl}: no healthy samples")
            continue
        recs.append(dict(label=label, slug=sl, nt=nt, lam=lam, n=n, runs=runs,
                         family=family,
                         model_def=mdef, ckey=ckey, rows=rows, diverged=div,
                         rlxf=RLXF.get((mdef, nt)),
                         seams=seam_of(runs, rows),
                         rev_max=max(r[1] for r in rows),
                         step_max=max(r[0] for r in rows),
                         last=last_window(rows, nt)))
    return recs


# --------------------------------------------------------------------------
# check

def check():
    """Reproduce every published anchor. Writes nothing. Exit 1 on mismatch."""
    recs = {r["slug"]: r for r in build()}
    bad = []

    print("\n--- CT anchors (revs 20-25, ledger 2026-08-25)")
    for sl, want in ANCHOR_CT.items():
        if sl not in recs:
            bad.append(f"{sl}: arm missing"); continue
        got, _, _ = ct_stats(recs[sl]["rows"], 20, 25)
        # tolerance is 2e-6 because the ledger records CTbar to 6 dp, so a
        # value at the rounding boundary (l3p4 = 0.0729005) can print either way
        ok = abs(got - want) < 2e-6
        print(f"  {sl:16s} want {want:.6f}  got {got:.6f}  {'OK' if ok else 'MISMATCH'}")
        if not ok:
            bad.append(f"{sl}: CT {got:.6f} != {want:.6f}")

    print("\n--- M2 anchors (circulation_te, revs 20-25)")
    profs = {}
    for sl in ANCHOR_CT:
        if sl in recs:
            profs[sl] = gamma_col(recs[sl]["runs"], recs[sl]["rows"], 20, 25,
                                  "circulation_te")
    for (a, b), want in ANCHOR_M2.items():
        if profs.get(a) is None or profs.get(b) is None:
            bad.append(f"{a}->{b}: profile missing"); continue
        e = eps_gamma(profs[a], profs[b])
        ok = abs(e["eps_max_pct"] - want) < 0.005
        print(f"  {a}->{b:14s} want {want:.3f}%  got {e['eps_max_pct']:.3f}%  "
              f"{'OK' if ok else 'MISMATCH'}")
        if not ok:
            bad.append(f"{a}->{b}: eps {e['eps_max_pct']:.3f} != {want}")

    print("\n--- seam anchors (gate 0.05%)")
    for sl, want in ANCHOR_SEAM.items():
        if sl not in recs or not recs[sl]["seams"]:
            bad.append(f"{sl}: no seam"); continue
        got = recs[sl]["seams"][0]["jump_pct"]
        ok = abs(got - want) < 0.0005
        print(f"  {sl:16s} want {want:.4f}%  got {got:.4f}%  "
              f"[{recs[sl]['seams'][0]['gate']}]  {'OK' if ok else 'MISMATCH'}")
        if not ok:
            bad.append(f"{sl}: seam {got:.4f} != {want}")

    print("\n--- coverage (rev_max vs step_max/NT, and vs ledger)")
    for r in sorted(recs.values(), key=lambda x: (x["nt"] or 0, x["slug"])):
        pred = r["step_max"] / r["nt"] if r["nt"] else float("nan")
        flag = "" if not r["nt"] or abs(pred - r["rev_max"]) / r["rev_max"] < 0.01 \
            else "  <-- REV/NT MISMATCH"
        div = f"  DIVERGED@{r['diverged']:.2f}" if r["diverged"] else ""
        print(f"  {r['slug']:20s} NT{str(r['nt']):4s} rev_max {r['rev_max']:6.2f} "
              f"(step/NT {pred:6.2f})  last5 {r['last']}{flag}{div}")
        if flag:
            bad.append(f"{r['slug']}: rev_max {r['rev_max']:.2f} != step/NT {pred:.2f}")

    if bad:
        print("\nFAILED:")
        for b in bad:
            print("  " + b)
        return 1
    print("\nAll anchors reproduced.")
    return 0


# --------------------------------------------------------------------------
# emit

def _w(path, header, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        for r in rows:
            w.writerow(r)
    for r in rows:
        for v in r:
            if isinstance(v, float) and not math.isfinite(v):
                raise SystemExit(f"non-finite value written to {path}")
    return len(rows)


def emit_history(recs, figdir):
    d = os.path.join(figdir, "fig_ct_history")
    n = 0
    for r in recs:
        n += _w(os.path.join(d, f"hist_{r['slug']}.csv"), ["rev", "ct"],
                [(f"{rev:.6f}", f"{ct:.8f}") for _, rev, ct in r["rows"]])
        bins = per_rev(r["rows"])
        good = complete_bins(r["rows"], r["nt"]) if r["nt"] else sorted(bins)
        _w(os.path.join(d, f"perrev_{r['slug']}.csv"),
           ["rev", "ct_mean", "ptp", "nsamp"],
           [(f"{b + 0.5:.3f}", f"{_mean(bins[b]):.8f}",
             f"{max(bins[b]) - min(bins[b]):.8f}", len(bins[b])) for b in good])
    seams = []
    for r in recs:
        for s in r["seams"]:
            seams.append((r["slug"], r["label"], r["nt"], f"{s['rev']:.4f}",
                          s["step"], f"{s['ct_before']:.8f}",
                          f"{s['ct_after']:.8f}", f"{s['jump_pct']:.4f}",
                          s["gate"]))
    _w(os.path.join(d, "seams.csv"),
       ["series", "label", "nt", "rev", "step", "ct_before", "ct_after",
        "jump_pct", "gate"], seams)
    _w(os.path.join(d, "windows.csv"), ["window", "series", "nt", "lo", "hi"],
       [("matched", r["slug"], r["nt"], MATCHED[0], MATCHED[1]) for r in recs] +
       [("last5", r["slug"], r["nt"], r["last"][0], r["last"][1])
        for r in recs if r["last"]])
    _w(os.path.join(d, "coverage.csv"),
       ["series", "label", "nt", "lam", "family", "model_def", "rlxf",
        "nrows", "step_max", "rev_max", "diverged_rev"],
       [(r["slug"], r["label"], r["nt"], r["lam"] if r["lam"] else "",
         r["family"], r["model_def"], r["rlxf"] if r["rlxf"] else "",
         len(r["rows"]),
         r["step_max"], f"{r['rev_max']:.4f}",
         f"{r['diverged']:.4f}" if r["diverged"] else "") for r in recs])
    print(f"  fig_ct_history: {len(recs)} arms, {n} history rows, "
          f"{len(seams)} seams")


def emit_rungs(recs, figdir):
    d = os.path.join(figdir, "fig_ct_rungs")
    rows = []
    for r in recs:
        for wname, win in (("matched", MATCHED), ("last5", r["last"])):
            if not win:
                continue
            good = [b for b in complete_bins(r["rows"], r["nt"] or 36)
                    if win[0] <= b <= win[1]]
            if len(good) < 2:
                continue
            try:
                ct, lo, hi = ct_stats(r["rows"], win[0], win[1])
            except SystemExit:
                continue
            rows.append((r["slug"], r["label"], wname, r["nt"],
                         r["lam"] if r["lam"] else "", r["n"], r["family"],
                         r["model_def"],
                         r["rlxf"] if r["rlxf"] else "", r["nt"],
                         f"{ct:.8f}", f"{ct - lo:.8f}", f"{hi - ct:.8f}",
                         win[0], win[1], len(good)))
    _w(os.path.join(d, "ct.csv"),
       ["series", "label", "window", "nt", "lam", "nrows", "family",
        "model_def", "rlxf", "param", "ct", "eminus", "eplus", "lo", "hi",
        "nrev"], rows)
    _w(os.path.join(d, "flags.csv"), ["series", "label", "nt", "status"],
       [(r["slug"], r["label"], r["nt"],
         "diverged" if r["diverged"] else
         ("legacy" if r["model_def"] == "naive" else
          ("unsettled" if r["nt"] and r["rev_max"] < 20 else "settled")))
        for r in recs])
    print(f"  fig_ct_rungs: {len(rows)} rung/window points")
    return rows


def _emit_gamma_fig(recs_by_slug, figname, members, ref_slug, figdir, windows,
                    note):
    """Emit Gamma(r/R) and Delta% for every (column, window) combination.

    `windows` is [(wname, lo, hi), ...].  Both a settled and a matched window
    are emitted wherever both exist, because revs 8-12 -- the only window every
    arm reaches -- is still deep in the NT36/NT72 transient (CTbar moves 3.8-4.4%
    between the two), so a matched-window Gamma comparison measures settling,
    not the ladder.  Emitting both is the only honest presentation.
    """
    d = os.path.join(figdir, figname)
    m2rows = []
    for wname, lo, hi in windows:
        for col, tag in (("circulation_te", "te"),
                         ("circulation_slice", "slice")):
            profs = {}
            for sl in members:
                r = recs_by_slug.get(sl)
                if not r:
                    continue
                p = gamma_col(r["runs"], r["rows"], lo, hi, col)
                if p is None:
                    print(f"    {figname}/{wname}/{tag}: {sl} has no "
                          f"circulation monitor")
                    continue
                profs[sl] = p
                _w(os.path.join(d, f"gamma_{wname}_{tag}_{sl}.csv"),
                   ["r", "gamma"],
                   [(f"{rr:.6f}", f"{gg:.8f}") for rr, gg in p])
            ref = profs.get(ref_slug)
            if ref is None:
                continue
            for sl, p in profs.items():
                if sl == ref_slug:
                    continue
                e = eps_gamma(ref, p)
                if e is None:
                    continue
                _w(os.path.join(d, f"delta_{wname}_{tag}_{sl}.csv"),
                   ["r", "dpct"],
                   [(f"{rr:.6f}", f"{dd:.6f}") for rr, dd in e["curve"]])
                m2rows.append((wname, tag, f"{ref_slug}->{sl}", ref_slug, sl,
                               lo, hi, f"{e['eps_max_pct']:.4f}",
                               f"{e['eps_rms_pct']:.4f}",
                               f"{e['r_worst']:.4f}", e["verdict"]))
    _w(os.path.join(d, "m2.csv"),
       ["window", "panel", "pair", "ref", "arm", "lo", "hi", "eps_max_pct",
        "eps_rms_pct", "r_worst", "verdict"], m2rows)
    _w(os.path.join(d, "window.csv"), ["window", "lo", "hi", "note"],
       [(w, lo, hi, note) for w, lo, hi in windows])
    print(f"  {figname}: {len(m2rows)} M2 pairs over {len(windows)} window(s)")
    return m2rows


# Settled window for the NT72 lambda ladder: 20-24 is five COMPLETE rev bins
# for all four arms (the ledger's 20-25 includes a partial bin 25 for three of
# them).  NT144 has no settled window at all -- that is the finding.
SETTLED_NT72 = (20, 24)
SETTLED_NT36 = (25, 29)


def emit_gamma(recs, figdir):
    by = {r["slug"]: r for r in recs}
    a = _emit_gamma_fig(
        by, "fig_gamma_span",
        ["nt72_l2p4", "nt72_l3p4", "nt72_l4p8", "nt72_n0_l4p8"], "nt72_l2p4",
        figdir,
        [("settled", SETTLED_NT72[0], SETTLED_NT72[1]),
         ("matched", MATCHED[0], MATCHED[1])],
        "lambda ladder at NT=72; settled=revs 20-24, matched=revs 8-12")
    b = _emit_gamma_fig(
        by, "fig_gamma_nt", ["nt36_l3p4", "nt72_l3p4", "nt144_l3p4"],
        "nt36_l3p4", figdir,
        [("matched", MATCHED[0], MATCHED[1])],
        "NT ladder at lambda=3.4; matched=revs 8-12 only -- NT144 reaches "
        "13.8 rev and has NO settled window")
    return a + b


# --------------------------------------------------------------------------
# spanwise loading, derived from circulation via Kutta-Joukowski
#
# There is NO sectional-force data on disk: SpanwiseLoadingMonitor exists in
# the source (src/FLOWPanel_simulate_monitors.jl:1675) but was not enabled for
# any csarc arm -- only monitors 02/03/04/05 were written. Loading is therefore
# derived, not measured.
#
# In hover the relative velocity at station r is Omega*r (inflow is a few % of
# it), so KJ gives dT/dr = rho * Gamma(r) * Omega * r per blade. With
# CT = T / (rho n^2 D^4), B blades, x = r/R, D = 2R and Omega = 2*pi*n:
#
#   dCT/dx = [ pi*B / (8*n*R^2) ] * Gamma(x) * x
#
# The bracket is KJ_K below. rho cancels entirely.
#
# VALIDATED: integrating this over the sampled stations reproduces the measured
# CTbar to 0.3-0.8% on all four settled NT72 arms (ratio 1.005/1.008/1.006/
# 1.003) using circulation_te. That is a genuine check of both the derivation
# and the identification of circulation_te as the bound circulation.
#
# circulation_slice does NOT share that normalization -- its KJ integral comes
# out at ratio ~0.47, i.e. ~2.11x small. An ABSOLUTE loading from the slice
# column would be wrong by a factor of two, so the slice panels are rescaled by
# the reference arm's te/slice integral ratio and are a SHAPE comparison only.
# The factor is written to loading_scale.csv and stated on the figure.
R_ROTOR = 0.119        # m, metadata.toml
N_BLADES = 2           # monitor03 `blade` column takes values 1,2
RPM = 5400.0
KJ_K = math.pi * N_BLADES / (8.0 * (RPM / 60.0) * R_ROTOR ** 2)


def loading_profile(prof):
    """dCT/dx from a Gamma(x) profile. Returns [(x, dctdx), ...]."""
    return [(x, KJ_K * g * x) for x, g in prof]


def integrate(curve):
    """Trapezoid over the sampled stations."""
    s = 0.0
    for i in range(len(curve) - 1):
        x0, y0 = curve[i]
        x1, y1 = curve[i + 1]
        s += 0.5 * (y0 + y1) * (x1 - x0)
    return s


def _emit_loading_fig(recs_by_slug, figname, members, ref_slug, figdir,
                      wname, lo, hi, note):
    d = os.path.join(figdir, figname)
    scales, rows = [], []
    curves = {}
    for col, tag in (("circulation_te", "te"), ("circulation_slice", "slice")):
        for sl in members:
            r = recs_by_slug.get(sl)
            if not r:
                continue
            p = gamma_col(r["runs"], r["rows"], lo, hi, col)
            if p is None:
                continue
            curves[(tag, sl)] = loading_profile(p)
    # slice rescale, from the reference arm only, so relative differences
    # between arms are preserved inside the panel
    ref_te = curves.get(("te", ref_slug))
    ref_sl = curves.get(("slice", ref_slug))
    scale = 1.0
    if ref_te and ref_sl:
        i_te, i_sl = integrate(ref_te), integrate(ref_sl)
        scale = i_te / i_sl if i_sl else 1.0
    scales.append(("slice", f"{scale:.4f}", f"{integrate(ref_te):.6f}" if ref_te
                   else "", f"{integrate(ref_sl):.6f}" if ref_sl else "",
                   ref_slug))
    _w(os.path.join(d, "loading_scale.csv"),
       ["panel", "scale_applied", "ref_int_te", "ref_int_slice", "ref"], scales)

    for (tag, sl), c in curves.items():
        k = scale if tag == "slice" else 1.0
        cc = [(x, y * k) for x, y in c]
        _w(os.path.join(d, f"load_{wname}_{tag}_{sl}.csv"), ["r", "dctdx"],
           [(f"{x:.6f}", f"{y:.8f}") for x, y in cc])
        r = recs_by_slug[sl]
        try:
            ct_meas, _, _ = ct_stats(r["rows"], lo, hi)
        except SystemExit:
            ct_meas = float("nan")
        ct_kj = integrate(cc)
        rows.append((sl, r["label"], tag, wname, f"{ct_kj:.6f}",
                     f"{ct_meas:.6f}" if math.isfinite(ct_meas) else "",
                     f"{ct_kj / ct_meas:.4f}" if math.isfinite(ct_meas)
                     and ct_meas else "", f"{k:.4f}"))
        # delta vs reference, on the reference's grid; the AREA under this
        # curve is that arm's contribution to Delta-CT, which is the whole
        # point of the panel
        ref = curves.get((tag, ref_slug))
        if ref is None or sl == ref_slug:
            continue
        refc = [(x, y * k) for x, y in ref]
        grid = [x for x, _ in refc]
        dc = [(x, _interp(cc, x) - _interp(refc, x)) for x in grid]
        _w(os.path.join(d, f"dload_{wname}_{tag}_{sl}.csv"), ["r", "ddctdx"],
           [(f"{x:.6f}", f"{y:.8f}") for x, y in dc])
        rows[-1] = rows[-1] + (f"{integrate(dc):+.6f}",)
    rows = [r if len(r) == 9 else r + ("",) for r in rows]
    _w(os.path.join(d, "loading.csv"),
       ["series", "label", "panel", "window", "ct_kj", "ct_measured",
        "ratio", "scale_applied", "dct_vs_ref"], rows)
    _w(os.path.join(d, "window.csv"), ["window", "lo", "hi", "note"],
       [(wname, lo, hi, note)])
    print(f"  {figname}: {len(rows)} arms, slice rescale x{scale:.3f}")
    return rows


def emit_loading(recs, figdir):
    by = {r["slug"]: r for r in recs}
    a = _emit_loading_fig(
        by, "fig_loading_span",
        ["nt72_l2p4", "nt72_l3p4", "nt72_l4p8", "nt72_n0_l4p8"], "nt72_l2p4",
        figdir, "settled", SETTLED_NT72[0], SETTLED_NT72[1],
        "lambda ladder at NT=72, settled revs 20-24, KJ from circulation")
    b = _emit_loading_fig(
        by, "fig_loading_nt", ["nt36_l3p4", "nt72_l3p4", "nt144_l3p4"],
        "nt36_l3p4", figdir, "matched", MATCHED[0], MATCHED[1],
        "NT ladder at lambda=3.4, matched revs 8-12 (TRANSIENT), KJ")
    return a + b


def emit(figdir, only=None):
    recs = build()
    print(f"\nEmitting into {figdir}")
    if only in (None, "history"):
        emit_history(recs, figdir)
    if only in (None, "rungs"):
        rows = emit_rungs(recs, figdir)
        print("\n  CTbar by rung:")
        print(f"  {'series':20s} {'window':9s} {'NT':>4s} {'lam':>5s} "
              f"{'CTbar':>10s} {'revs':>9s}")
        for r in rows:
            if r[6] != "csarc":
                continue
            print(f"  {r[0]:20s} {r[2]:9s} {str(r[3]):>4s} {str(r[4]):>5s} "
                  f"{float(r[10]):10.6f} {str(r[13]) + '-' + str(r[14]):>9s}")
    if only in (None, "loading"):
        ld = emit_loading(recs, figdir)
        print("\n  Spanwise loading (KJ from circulation; dCT/dx = "
              f"{KJ_K:.4f}*Gamma*x):")
        print(f"  {'series':16s} {'panel':6s} {'CT_KJ':>9s} {'CT_meas':>9s} "
              f"{'ratio':>6s} {'dCT vs ref':>11s}")
        for r in ld:
            print(f"  {r[0]:16s} {r[2]:6s} {float(r[4]):9.6f} "
                  f"{r[5]:>9s} {r[6]:>6s} {r[8]:>11s}")
    if only in (None, "gamma"):
        m2 = emit_gamma(recs, figdir)
        print("\n  M2 (gate 1.0%):")
        for m in m2:
            print(f"  {m[0]:8s} {m[1]:6s} {m[2]:30s} revs {m[5]}-{m[6]:<3} "
                  f"eps_max {float(m[7]):6.3f}%  worst r/R {float(m[9]):.3f}  "
                  f"{m[10]}")
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("check", help="reproduce published anchors; writes nothing")
    e = sub.add_parser("emit", help="write the backing CSVs")
    e.add_argument("--figdir", default=FIGDIR_DEFAULT)
    e.add_argument("--only",
                   choices=["history", "rungs", "gamma", "loading"])
    a = ap.parse_args()
    if a.cmd == "check":
        return check()
    return emit(a.figdir, a.only)


if __name__ == "__main__":
    sys.exit(main())
