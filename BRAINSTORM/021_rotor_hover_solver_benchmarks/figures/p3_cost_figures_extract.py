#!/usr/bin/env python3
"""Emit TikZ backing CSVs for Phase 3 cost-per-step and total-cost figures."""
import csv, os, collections

SRC = os.path.dirname(os.path.abspath(__file__))
FIG = "/Users/ryan/Dropbox/research/projects/FLOWPanel.jl/BRAINSTORM/021_rotor_hover_solver_benchmarks/figures"
MODES = ["cold", "prev", "extrap"]
NT = 36.0  # steps per rev

# (rung, config, warmstart, step) -> row, last-wins dedup
data = {}
for rung in ["R1", "R2", "R3"]:
    with open(f"{SRC}/{rung}.csv") as f:
        for row in csv.DictReader(f):
            if row["n_steps"] == "468":
                continue  # checkpoint rows
            key = (rung, row["config"], row["warmstart"], int(row["step"]))
            data[key] = row

# per-step files
series = collections.defaultdict(list)  # (mode, rung, config) -> [(step, tnet)]
for (rung, cfg, mode, step), row in sorted(data.items()):
    series[(mode, rung, cfg)].append((step, float(row["t_step_net"])))

counts = collections.Counter()
for mode in MODES:
    d = f"{FIG}/p3_cost_per_step_{mode}"
    os.makedirs(d, exist_ok=True)
    for (m, rung, cfg), pts in series.items():
        if m != mode:
            continue
        with open(f"{d}/{rung}_{cfg}.csv", "w") as f:
            f.write("rev,t_step_net\n")
            for step, t in pts:
                f.write(f"{step/NT:.6f},{t:.4f}\n")
        counts[(mode, rung, cfg)] = len(pts)

# totals: only complete arms (all present modes have full step count)
per_arm = collections.defaultdict(dict)  # (rung,cfg) -> mode -> (n, total)
for (mode, rung, cfg), pts in series.items():
    per_arm[(rung, cfg)][mode] = (len(pts), sum(t for _, t in pts))

nmax = max(n for d in per_arm.values() for n, _ in d.values())
os.makedirs(f"{FIG}/p3_total_cost_vs_rung", exist_ok=True)
solvers = sorted({cfg for _, cfg in per_arm})
for cfg in solvers:
    with open(f"{FIG}/p3_total_cost_vs_rung/{cfg}.csv", "w") as f:
        f.write("rung,total_cold,total_prev,total_extrap\n")
        for rung in ["R1", "R2", "R3"]:
            d = per_arm.get((rung, cfg), {})
            vals = []
            complete = all(n == nmax for n, _ in d.values()) and d
            for m in MODES:
                n, tot = d.get(m, (0, 0.0))
                vals.append(f"{tot:.2f}" if (n == nmax) else "nan")
            f.write(f"{rung},{','.join(vals)}\n")

print("steps/rev assumed:", NT, "| full arm segment length:", nmax)
for k in sorted(counts):
    print(k, counts[k])
