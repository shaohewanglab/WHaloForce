"""Figure 2k: per-FOV bright-FA lifetime-beating coupling r, Vin-WF vs Vin-cpH
strip plot. Filled dots = surrogate-nonsignificant, open rings = surrogate
significant (p < 0.05); horizontal bars = construct means. Exact Mann-Whitney U
p-value (two-sided) between constructs is annotated.

Reads data/Fig2/Fig2k.csv and writes figures/Fig2k.svg and .png.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from style import apply_style, save, DATA_DIR, p_text  # noqa: E402

apply_style()

COLOR = {"Vin-WF": "#4C72B0", "Vin-cpH": "#DD8452"}
ORDER = ["Vin-WF", "Vin-cpH"]


def main():
    df = pd.read_csv(DATA_DIR / "Fig2" / "Fig2k.csv")
    grp = {c: df[df["construct"] == c] for c in ORDER}

    U, p_mwu = mannwhitneyu(grp["Vin-WF"]["coupling_r"],
                            grp["Vin-cpH"]["coupling_r"], alternative="two-sided")
    print(f"MWU U={U:.0f} p={p_mwu:.3e}  "
          f"(n Vin-WF={len(grp['Vin-WF'])}, n Vin-cpH={len(grp['Vin-cpH'])})")

    fig, ax = plt.subplots(figsize=(1.8, 1.8))
    rng = np.random.RandomState(0)
    for i, c in enumerate(ORDER):
        d = grp[c]
        r = d["coupling_r"].to_numpy()
        sig = d["significant"].to_numpy().astype(bool)
        x = rng.normal(i, 0.07, len(d))
        ax.scatter(x, r, s=9, color=COLOR[c], alpha=0.85, linewidths=0, zorder=3)
        ax.scatter(x[sig], r[sig], s=34, facecolors="none", edgecolors="0.4",
                   linewidths=0.8, zorder=4)
        ax.plot([i - 0.26, i + 0.26], [r.mean()] * 2, color="0.15", lw=1.2,
                zorder=5)

    ymax = float(df["coupling_r"].max())
    yb = ymax + 0.07
    ax.plot([0, 0, 1, 1], [yb, yb + 0.02, yb + 0.02, yb], color="0.15", lw=0.6,
            zorder=6)
    p_text(ax, 0.5, yb + 0.02, p_mwu)

    ax.set_xticks([0, 1])
    ax.set_xticklabels(ORDER)
    ax.set_xlim(-0.55, 1.55)
    ax.set_ylim(0, yb + 0.16)
    ax.set_ylabel("Lifetime-beating coupling r")
    ax.spines[["top", "right"]].set_visible(False)

    ring = Line2D([], [], marker="o", linestyle="none", markerfacecolor="none",
                  markeredgecolor="0.4", markeredgewidth=0.8, markersize=6)
    ax.legend([ring], ["surrogate\np < 0.05"], loc="upper right", frameon=False,
              fontsize=7, handletextpad=0.3, borderpad=0.1,
              bbox_to_anchor=(1.04, 0.78))

    fig.subplots_adjust(left=0.24, right=0.97, top=0.95, bottom=0.10)
    save(fig, "Fig2k")


if __name__ == "__main__":
    main()
