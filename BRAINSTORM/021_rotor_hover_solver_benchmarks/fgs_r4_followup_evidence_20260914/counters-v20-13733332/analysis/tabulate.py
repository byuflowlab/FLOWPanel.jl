#!/usr/bin/env python3
"""Tabulate perf-counter and stage-activity data for the v20 counters campaign.

Reads (per arm, under rundir/<arm>/):
  perf-stat.csv                       (perf stat -x, CSV)
  results/stage_thread_activity.csv
  results/counter_validation.csv

Writes analysis/counters_summary.md (sibling of rundir/).
"""
import csv
from pathlib import Path
from collections import defaultdict

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RUNDIR = ROOT / "rundir"
ARMS = ["j4-b1", "j64-b1"]

with open(RUNDIR / "clock_ticks_per_second.txt") as f:
    CLK_TCK = float(f.read().strip())


def parse_perf_stat(path):
    """Parse `perf stat -x,` CSV lines into event -> dict."""
    events = {}
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split(",")
            # value,unit,event,run-time,percentage,metric-value,metric-unit
            value_s, unit, event = parts[0], parts[1], parts[2]
            runtime_ns = parts[3] if len(parts) > 3 else ""
            pct = parts[4] if len(parts) > 4 else ""
            try:
                value = float(value_s)
            except ValueError:
                value = None
            events[event] = dict(
                value=value, unit=unit,
                runtime_ns=float(runtime_ns) if runtime_ns else None,
                pct=float(pct) if pct else None,
            )
    return events


def task_a(arm):
    ev = parse_perf_stat(RUNDIR / arm / "perf-stat.csv")
    cycles = ev["cycles:u"]["value"]
    instr = ev["instructions:u"]["value"]
    refs = ev["cache-references:u"]["value"]
    misses = ev["cache-misses:u"]["value"]
    clock_ms = ev["task-clock:u"]["value"]
    ipc = instr / cycles
    miss_ratio = misses / refs
    pcts = {k: v["pct"] for k, v in ev.items()}
    muxed = [k for k, p in pcts.items() if p is not None and p < 100.0]
    return dict(
        arm=arm, cycles=cycles, instructions=instr, ipc=ipc,
        cache_references=refs, cache_misses=misses, miss_ratio=miss_ratio,
        task_clock_ms=clock_ms, pcts=pcts, muxed=muxed,
    )


def task_b(arm):
    path = RUNDIR / arm / "results" / "stage_thread_activity.csv"
    # key: (sequence, stage) -> accumulate
    groups = defaultdict(lambda: dict(span=None, tids=set(), total_ticks=0,
                                       n_excluded=-1, complete=0, total=0))
    excluded_count = defaultdict(int)
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            key = (int(row["sequence"]), row["stage"])
            g = groups[key]
            g["span"] = float(row["span_seconds"])
            g["tids"].add(row["tid"])
            ticks = int(row["cpu_ticks"])
            g["total"] += 1
            if ticks == -1:
                excluded_count[key] += 1
            else:
                g["total_ticks"] += ticks
            if row["complete_endpoints"].strip().lower() == "true":
                g["complete"] += 1
    rows = []
    for key in sorted(groups.keys(), key=lambda k: (k[0], k[1])):
        g = groups[key]
        seq, stage = key
        n_threads = len(g["tids"])
        busy_cpu_seconds = g["total_ticks"] / CLK_TCK
        span = g["span"]
        avg_active = busy_cpu_seconds / span if span else float("nan")
        n_excl = excluded_count[key]
        pct_incomplete = 100.0 * (g["total"] - g["complete"]) / g["total"] if g["total"] else float("nan")
        rows.append(dict(
            sequence=seq, stage=stage, span_seconds=span, n_threads=n_threads,
            total_cpu_ticks=g["total_ticks"], busy_cpu_seconds=busy_cpu_seconds,
            avg_active_threads=avg_active, n_excluded=n_excl,
            pct_incomplete=pct_incomplete,
        ))
    return rows


def task_c(arm, stage_rows):
    path = RUNDIR / arm / "results" / "counter_validation.csv"
    modes = {}
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            modes[row["mode"]] = float(row["diagnostic_seconds"])
    span_sum = sum(r["span_seconds"] for r in stage_rows)
    activity_diag = modes.get("activity")
    counters_diag = modes.get("counters")
    baseline_diag = modes.get("baseline")
    diff = span_sum - activity_diag if activity_diag is not None else float("nan")
    return dict(
        arm=arm, baseline_s=baseline_diag, activity_s=activity_diag,
        counters_s=counters_diag, span_sum_s=span_sum, diff_span_minus_activity=diff,
    )


def fmt(x, nd=4):
    if x is None:
        return "NA"
    if isinstance(x, float):
        return f"{x:.{nd}g}"
    return str(x)


def main():
    lines = []
    lines.append(f"CLK_TCK = {CLK_TCK:g}\n")

    # Task A
    lines.append("## Task A: perf counters (one warmed prepared solve, FIFO-gated, user-space only)\n")
    lines.append("| arm | cycles | instructions | IPC | cache-refs | cache-misses | miss ratio | task-clock (ms) | multiplexed events |")
    lines.append("|---|---|---|---|---|---|---|---|---|")
    a_results = {}
    for arm in ARMS:
        r = task_a(arm)
        a_results[arm] = r
        mux = ", ".join(r["muxed"]) if r["muxed"] else "none (all 100%)"
        lines.append(
            f"| {arm} | {r['cycles']:.0f} | {r['instructions']:.0f} | {r['ipc']:.3f} | "
            f"{r['cache_references']:.0f} | {r['cache_misses']:.0f} | {r['miss_ratio']*100:.3f}% | "
            f"{r['task_clock_ms']:.2f} | {mux} |"
        )
    lines.append("")

    # Task B
    lines.append("## Task B: stage thread activity (per arm, ordered by sequence)\n")
    b_results = {}
    for arm in ARMS:
        rows = task_b(arm)
        b_results[arm] = rows
        lines.append(f"### {arm}\n")
        lines.append("| sequence | stage | span_seconds | n_threads | total_cpu_ticks | busy_cpu_seconds | avg_active_threads | n_excluded(-1) | pct_incomplete |")
        lines.append("|---|---|---|---|---|---|---|---|---|")
        for r in rows:
            lines.append(
                f"| {r['sequence']} | {r['stage']} | {r['span_seconds']:.6f} | {r['n_threads']} | "
                f"{r['total_cpu_ticks']} | {r['busy_cpu_seconds']:.4f} | {r['avg_active_threads']:.3f} | "
                f"{r['n_excluded']} | {r['pct_incomplete']:.2f}% |"
            )
        lines.append("")

    # Task C
    lines.append("## Task C: cross-checks\n")
    lines.append("| arm | baseline diagnostic_s | activity diagnostic_s | counters diagnostic_s | sum(span_seconds) | span_sum - activity_diag |")
    lines.append("|---|---|---|---|---|---|")
    for arm in ARMS:
        c = task_c(arm, b_results[arm])
        lines.append(
            f"| {arm} | {fmt(c['baseline_s'])} | {fmt(c['activity_s'])} | {fmt(c['counters_s'])} | "
            f"{fmt(c['span_sum_s'])} | {fmt(c['diff_span_minus_activity'])} |"
        )
    lines.append("")

    # Anomalies
    lines.append("## Anomalies\n")
    anomalies = []
    for arm in ARMS:
        if a_results[arm]["muxed"]:
            anomalies.append(f"{arm}: multiplexed perf events {a_results[arm]['muxed']}")
        for r in b_results[arm]:
            if r["total_cpu_ticks"] == 0:
                anomalies.append(f"{arm}: stage '{r['stage']}' (seq {r['sequence']}) has zero total cpu_ticks")
            if r["n_excluded"] > 0:
                anomalies.append(f"{arm}: stage '{r['stage']}' (seq {r['sequence']}) excluded {r['n_excluded']} rows with cpu_ticks==-1")
            if r["pct_incomplete"] > 0:
                anomalies.append(f"{arm}: stage '{r['stage']}' (seq {r['sequence']}) has {r['pct_incomplete']:.2f}% incomplete endpoints")
    if not anomalies:
        anomalies.append("none found: all perf events at 100% (no multiplexing), no -1 cpu_ticks, no zero-tick stages, no incomplete endpoints.")
    for a in anomalies:
        lines.append(f"- {a}")

    out_path = HERE / "counters_summary.md"
    out_path.write_text("\n".join(lines) + "\n")
    print(f"wrote {out_path}")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
