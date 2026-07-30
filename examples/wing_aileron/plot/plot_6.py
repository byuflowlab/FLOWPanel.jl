import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# ------------------------------------------------------------
# Load CSV
# ------------------------------------------------------------
df = pd.read_csv("examples/wing_aileron/results/third.csv")

# ------------------------------------------------------------
# Filter panel count
# ------------------------------------------------------------
df = df[df["nps"] == 15152]

coupled = df[df["solver"] == "BackslashCoupled"].sort_values("AOA")
iterative = df[df["solver"] == "BackslashIterative"].sort_values("AOA")

# ------------------------------------------------------------
# Extract data
# ------------------------------------------------------------
aoa = coupled["AOA"].values

build_c = coupled["t_build"].values
solve_c = coupled["t_solve"].values

build_i = iterative["t_build"].values
solve_i = iterative["t_solve"].values

# ------------------------------------------------------------
# Plot setup
# ------------------------------------------------------------
fig, ax = plt.subplots(figsize=(12, 6))

x = np.arange(len(aoa))
width = 0.38

# ------------------------------------------------------------
# BackslashCoupled (stacked)
# ------------------------------------------------------------
ax.bar(
    x - width/2,
    build_c,
    width,
    label="Coupled build",
)

ax.bar(
    x - width/2,
    solve_c,
    width,
    bottom=build_c,
    label="Coupled solve",
)

# ------------------------------------------------------------
# BackslashIterative
# ------------------------------------------------------------
ax.bar(
    x + width/2,
    build_i,
    width,
    label="Iterative build",
)

ax.bar(
    x + width/2,
    solve_i,
    width,
    bottom=build_i,
    label="Iterative solve",
)

# ------------------------------------------------------------
# Formatting
# ------------------------------------------------------------
ax.set_xlabel("Angle of Attack (deg)")
ax.set_ylabel("Time (s)")
# ax.set_title("15152 Panels: Build + Solve Timing")

# Use actual AOA values as ticks
ax.set_xticks(x)
ax.set_xticklabels([f"{v:.1f}" for v in aoa], rotation=45)

ax.legend()
ax.set_yscale("log")
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
# ax.grid(axis="y", alpha=0.3)

plt.tight_layout()
plt.show()