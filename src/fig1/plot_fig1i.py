"""Figure 1i: quenched-state dwell-time distribution per force condition.

IQR box (no whiskers) plus jittered points of the quenched-state dwell time for
each force, on a log y-axis. Only complete quenched->unquenched (QU) transitions
contribute a dwell time; censored observations (END/BL) are excluded. Under-
sampled 10 and 15 pN conditions are dropped, matching the original
plot_QU_force_dependence.py defaults. Reads data/Fig1/Fig1i.csv and writes
figures/Fig1i.svg and .png.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from style import apply_style, save, DATA_DIR  # noqa: E402

apply_style()

EXCLUDE = {10, 15}
# Okabe-Ito palette. 8 pN is the reddish-purple #CC79A7, sampled from the
# published figure -- paler than the #CC0066 used previously.
COLORS = {3: "#E6AB02", 4: "#000000", 5: "#0072B2",
          6: "#D55E00", 7: "#009E73", 8: "#CC79A7"}


def _rgba(hex_color, alpha):
    r, g, b = plt.matplotlib.colors.to_rgb(hex_color)
    return (r, g, b, alpha)


def main():
    df = pd.read_csv(DATA_DIR / "Fig1" / "Fig1i.csv")
    df = df[~df["force_pN"].isin(EXCLUDE)]
    obs = df[df["observed"]]
    forces = [f for f in sorted(obs["force_pN"].unique())
              if (obs["force_pN"] == f).sum() > 0]
    box_data = [obs.loc[obs["force_pN"] == f, "dwell_time_s"].to_numpy()
                for f in forces]
    positions = list(range(len(forces)))

    fig, ax = plt.subplots(figsize=(1.8, 2.0))
    rng = np.random.default_rng(0)

    bp = ax.boxplot(box_data, positions=positions, widths=0.6,
                    patch_artist=True, showfliers=False, showcaps=False,
                    whis=0, whiskerprops={"linewidth": 0},
                    medianprops={"linewidth": 1.2, "solid_capstyle": "butt"})
    for patch, f in zip(bp["boxes"], forces):
        patch.set_facecolor(_rgba(COLORS[f], 0.20))
        patch.set_edgecolor(_rgba(COLORS[f], 0.55))
        patch.set_linewidth(0.5)
    for med, f in zip(bp["medians"], forces):
        med.set_color(COLORS[f])

    for pos, f, data in zip(positions, forces, box_data):
        jitter = rng.uniform(-0.12, 0.12, size=data.size)
        ax.scatter(pos + jitter, data, s=8, color=COLORS[f],
                   edgecolors="none", alpha=0.85, zorder=3)

    all_vals = np.concatenate(box_data)
    ax.set_yscale("log")
    ax.set_ylim(all_vals.min() / 2.0, all_vals.max() * 2.0)
    ax.set_xticks(positions)
    ax.set_xticklabels([f"{f} pN" for f in forces])
    ax.set_xlabel("Force condition")
    ax.set_ylabel("Quenched state dwell time (s)")
    ax.grid(True, axis="y", which="both", linestyle="-", linewidth=0.3,
            color="0.85", alpha=0.8)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)

    handles = [Line2D([0], [0], color="0.2", lw=1.2, label="Median"),
               Patch(facecolor="0.75", edgecolor="0.45", lw=0.5, alpha=0.6,
                     label="IQR")]
    ax.legend(handles=handles, loc="upper right", frameon=False, fontsize=6)

    fig.tight_layout(pad=0.4)
    save(fig, "Fig1i")


if __name__ == "__main__":
    main()
