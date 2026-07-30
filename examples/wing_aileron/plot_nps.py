import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# ------------------------------------------------------------
# Load CSV
# ------------------------------------------------------------
df = pd.read_csv("examples/wing_aileron/krylov_solvers.csv")

# ------------------------------------------------------------
# Panel counts present in the data (sorted ascending), one bar group each
# ------------------------------------------------------------
df_nps = np.sort(df["nps"].unique())

m_b_c = []
m_s_c = []
m_b_i = []
m_s_i = []

fig, ax = plt.subplots(figsize=(12, 6))

x = np.arange(len(df_nps))
width = 0.38

for df_n in df_nps:
    df_panel = df[df["nps"] == df_n]
    coupled = df_panel[df_panel["solver"] == "KrylovCoupled"].sort_values("AOA")
    iterative = df_panel[df_panel["solver"] == "KrylovCoupled-FMM"].sort_values("AOA")

    # ------------------------------------------------------------
    # Extract data
    # ------------------------------------------------------------
    build_c = coupled["t_build"].values
    solve_c = coupled["t_solve"].values

    build_i = iterative["t_build"].values
    solve_i = iterative["t_solve"].values

    avg_b_c = np.mean(build_c)
    avg_s_c = np.mean(solve_c)
    avg_b_i = np.mean(build_i)
    avg_s_i = np.mean(solve_i)

    m_b_c.append(avg_b_c)
    m_s_c.append(avg_s_c)
    m_b_i.append(avg_b_i)
    m_s_i.append(avg_s_i)

m_b_c = np.array(m_b_c)
m_s_c = np.array(m_s_c)
m_b_i = np.array(m_b_i)
m_s_i = np.array(m_s_i)

# ------------------------------------------------------------
# KrylovCoupled (stacked)
# ------------------------------------------------------------
bars_b_c = ax.bar(
    x - width/2,
    m_b_c,
    width,
    label="Direct build",
)

bars_s_c = ax.bar(
    x - width/2,
    m_s_c,
    width,
    bottom=m_b_c,
    label="Direct solve",
)

# ------------------------------------------------------------
# KrylovCoupled-FMM
# ------------------------------------------------------------
bars_b_i = ax.bar(
    x + width/2,
    m_b_i,
    width,
    label="FMM build",
)

bars_s_i = ax.bar(
    x + width/2,
    m_s_i,
    width,
    bottom=m_b_i,
    label="FMM solve",
)

# ------------------------------------------------------------
# Total time per solver/panel-count, used for the percentage callout
# ------------------------------------------------------------
total_c = m_b_c + m_s_c
total_i = m_b_i + m_s_i

# ------------------------------------------------------------
# Callouts: % difference between solvers' total time (Iterative vs Coupled)
# ------------------------------------------------------------
pct_diff = 100 * (total_i - total_c) / total_c
for xi, tc, ti, pd in zip(x, total_c, total_i, pct_diff):
    label = f"{abs(pd):.1f}% {'faster' if pd < 0 else 'slower'}"
    y_top = max(tc, ti)
    ax.annotate(
        label,
        xy=(xi, y_top),
        xytext=(0, 14),
        textcoords="offset points",
        ha="center",
        fontsize=9,
        fontweight="bold",
    )

# ------------------------------------------------------------
# Formatting
# ------------------------------------------------------------
ax.set_xlabel("Panel Count")
ax.set_ylabel("Time (s)")

# Use actual panel counts as ticks
ax.set_xticks(x)
ax.set_xticklabels([f"{v}" for v in df_nps], rotation=45)
ax.legend(loc="upper left")
ax.set_yscale("log")
ax.margins(y=0.15)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
# ax.grid(axis="y", alpha=0.3)

plt.tight_layout()
plt.savefig("examples/wing_aileron/wing_aileron_benchmarks_nps_krylov.png")
plt.show()