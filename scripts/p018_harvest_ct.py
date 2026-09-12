#!/usr/bin/env python3
"""Harvest CT windowed scores for BRAINSTORM 018 runs.

Reads either CT_per_rev.csv (if run completed) or reconstructs from force monitor
(if run was killed before end-of-run CSV writing).

Outputs windowed CT mean ± standard error (SE) over the window revs.
"""

import csv
import math
import sys

RPM_DEFAULT = 5400.0

def harvest_from_ct_per_rev(csv_path, run_name, nt, job_id):
    """Read CT_per_rev.csv and compute windowed stats."""
    rows = []
    with open(csv_path) as fh:
        reader = csv.DictReader(fh)
        for rec in reader:
            try:
                in_window = rec['in_convergence_window'].strip().lower() == 'true'
                ct = float(rec['CT_mean'])
                rev = int(rec['rev_block'])
                step_start = int(rec['step_start'])
                step_stop = int(rec['step_stop'])
                rows.append((rev, ct, in_window, step_start, step_stop))
            except (ValueError, KeyError):
                continue

    if not rows:
        return None

    # Filter to convergence window
    windowed = [ct for rev, ct, in_window, _, _ in rows if in_window]
    if not windowed:
        return None

    # Compute mean and std
    mean_ct = sum(windowed) / len(windowed)
    if len(windowed) > 1:
        variance = sum((x - mean_ct) ** 2 for x in windowed) / (len(windowed) - 1)
        std_ct = math.sqrt(variance)
        sem_ct = std_ct / math.sqrt(len(windowed))  # Standard error of the mean
    else:
        std_ct = 0.0
        sem_ct = 0.0

    # Calculate spread (std as % over window revs)
    rel_std = (100.0 * std_ct / mean_ct) if mean_ct != 0 else 0.0

    # Get step coverage
    windowed_rows = [(s, e) for rev, ct, in_window, s, e in rows if in_window]
    if windowed_rows:
        min_step = min(s for s, e in windowed_rows)
        max_step = max(e for s, e in windowed_rows)
        step_coverage = f"{min_step}-{max_step}"
    else:
        step_coverage = "N/A"

    return {
        'run': run_name,
        'job': job_id,
        'nt': nt,
        'ct_mean': mean_ct,
        'ct_std': std_ct,
        'ct_sem': sem_ct,
        'spread_pct': rel_std,
        'spread_sem_pct': 100.0 * sem_ct / mean_ct if mean_ct != 0 else 0.0,
        'n_revs': len(windowed),
        'step_coverage': step_coverage,
        'source': 'CT_per_rev.csv',
        'n_rows': len(rows),
        'complete': True,
    }


def harvest_from_force_monitor(force_csv, run_name, nt, job_id, rpm=RPM_DEFAULT):
    """Reconstruct CT from force monitor and estimate per-rev windows.

    Since the run was killed, we have no pre-computed convergence window flag,
    so we estimate it as the final 10 revolutions of available data.
    """
    rows = []
    with open(force_csv) as fh:
        reader = csv.DictReader(fh)
        for rec in reader:
            try:
                cfx = float(rec['CFx'])
                time = float(rec['time'])
                step = int(rec['step'])
                # CT = -CFx (from BRAINSTORM/018 conventions)
                ct = -cfx
                rev = time * rpm / 60.0
                if math.isfinite(ct) and abs(ct) > 1e-12:  # Skip placeholder rows
                    rows.append((step, rev, ct))
            except (ValueError, KeyError):
                continue

    if not rows:
        return None

    # Group into per-rev bins (integer revolution)
    bins = {}
    for step, rev, ct in rows:
        rev_int = int(math.floor(rev + 1e-9))
        if rev_int not in bins:
            bins[rev_int] = []
        bins[rev_int].append(ct)

    if not bins:
        return None

    # Compute per-rev means
    rev_means = []
    rev_nums = sorted(bins.keys())
    for rev_num in rev_nums:
        rev_means.append((rev_num, sum(bins[rev_num]) / len(bins[rev_num])))

    if not rev_means:
        return None

    # Estimate convergence window as final 10 complete revs
    # (or fewer if run is shorter)
    n_window = min(10, len(rev_means))
    windowed = [ct for _, ct in rev_means[-n_window:]]

    if not windowed:
        return None

    # Compute mean and std
    mean_ct = sum(windowed) / len(windowed)
    if len(windowed) > 1:
        variance = sum((x - mean_ct) ** 2 for x in windowed) / (len(windowed) - 1)
        std_ct = math.sqrt(variance)
        sem_ct = std_ct / math.sqrt(len(windowed))
    else:
        std_ct = 0.0
        sem_ct = 0.0

    rel_std = (100.0 * std_ct / mean_ct) if mean_ct != 0 else 0.0
    rel_sem = (100.0 * sem_ct / mean_ct) if mean_ct != 0 else 0.0

    # Get step coverage for the window
    window_rev_nums = [r for r, _ in rev_means[-n_window:]]
    if window_rev_nums:
        steps_in_window = [r for r, rev, ct in rows
                           if int(math.floor(rev + 1e-9)) in window_rev_nums]
        if steps_in_window:
            min_step = min(steps_in_window)
            max_step = max(steps_in_window)
            step_coverage = f"{min_step}-{max_step}"
        else:
            step_coverage = "N/A"
    else:
        step_coverage = "N/A"

    max_rev = rev_nums[-1] if rev_nums else None

    return {
        'run': run_name,
        'job': job_id,
        'nt': nt,
        'ct_mean': mean_ct,
        'ct_std': std_ct,
        'ct_sem': sem_ct,
        'spread_pct': rel_std,
        'spread_sem_pct': rel_sem,
        'n_revs': len(windowed),
        'step_coverage': step_coverage,
        'source': f'force_monitor (final {n_window} revs)',
        'n_rows': len(rows),
        'complete': False,
        'max_rev': max_rev,
    }


def report(result, ref_ct=0.070775):
    """Print result as a table row."""
    if result is None:
        return None

    # Δ vs reference NT36 CT = 0.070775
    delta = 100.0 * (result['ct_mean'] - ref_ct) / ref_ct

    return {
        'case': result['run'].split('_')[-1] if '_' in result['run'] else result['run'],
        'job': result['job'],
        'nt': result['nt'],
        'ct': result['ct_mean'],
        'ct_sem': result['ct_sem'],
        'spread': result['spread_pct'],
        'delta_vs_ref': delta,
        'coverage': f"{result['n_revs']} revs, {result['step_coverage']}",
        'source': result['source'],
    }


if __name__ == '__main__':
    # Test run 1: has CT_per_rev.csv
    r1 = harvest_from_ct_per_rev(
        '/private/tmp/claude-502/-Users-ryan-Dropbox-research-projects-FLOWPanel-jl/64e502fc-4868-4e85-9e5d-53ef86f2df4a/scratchpad/sfs3nb_CT_per_rev.csv',
        'p018_csarc_l3p0_3r_sfs3nb',
        36,
        13593711
    )

    # Test run 2: killed, needs force monitor reconstruction
    r2 = harvest_from_force_monitor(
        '/private/tmp/claude-502/-Users-ryan-Dropbox-research-projects-FLOWPanel-jl/64e502fc-4868-4e85-9e5d-53ef86f2df4a/scratchpad/sv_h2p0_force_monitor.csv',
        'p018_csarc_l3p0_3r_sv_h2p0',
        36,
        13592759
    )

    print("\n" + "=" * 100)
    print("BRAINSTORM 018 HARVEST RESULTS")
    print("=" * 100)

    if r1:
        rep1 = report(r1)
        if rep1:
            print(f"\nRun 1: _3r_sfs3nb (Job {rep1['job']}, NT{rep1['nt']})")
            print(f"  CT = {rep1['ct']:.6f} ± {rep1['ct_sem']:.6f} (SE)")
            print(f"  Spread over {rep1['coverage']}: {rep1['spread']:.2f}%")
            print(f"  Δ vs ref NT36 (0.070775): {rep1['delta_vs_ref']:+.2f}%")
            print(f"  Source: {rep1['source']}")

    if r2:
        rep2 = report(r2)
        if rep2:
            print(f"\nRun 2: _3r_sv_h2p0 (Job {rep2['job']}, NT{rep2['nt']})")
            print(f"  CT = {rep2['ct']:.6f} ± {rep2['ct_sem']:.6f} (SE)")
            print(f"  Spread over {rep2['coverage']}: {rep2['spread']:.2f}%")
            print(f"  Δ vs ref NT36 (0.070775): {rep2['delta_vs_ref']:+.2f}%")
            print(f"  Source: {rep2['source']}")
            print(f"  Note: Run killed at rev {r2['max_rev']:.2f} (expected 30 revs for full convergence window)")

    print("\n" + "=" * 100)
    if r1 and r2:
        print("Summary Table:")
        print(f"{'Case':<20} {'Job':<10} {'NT':<4} {'Windowed CT':<12} {'SE':<10} {'Spread':<8} {'Δ vs ref':<10} {'Coverage':<25}")
        print("-" * 100)
        rep1 = report(r1)
        rep2 = report(r2)
        print(f"{'_3r_sfs3nb':<20} {rep1['job']:<10} {rep1['nt']:<4} {rep1['ct']:<12.6f} {rep1['ct_sem']:<10.6f} {rep1['spread']:<8.2f}% {rep1['delta_vs_ref']:<+10.2f}% {rep1['coverage']:<25}")
        print(f"{'_3r_sv_h2p0':<20} {rep2['job']:<10} {rep2['nt']:<4} {rep2['ct']:<12.6f} {rep2['ct_sem']:<10.6f} {rep2['spread']:<8.2f}% {rep2['delta_vs_ref']:<+10.2f}% {rep2['coverage']:<25}")
