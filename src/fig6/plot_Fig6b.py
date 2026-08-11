"""
Fig. 6b: spermathecal local lifetime contrast (tau_Spt(norm)) versus adult age,
for three constructs (CeLam-WF, CeLam-cpH, CeLam-C-WF).

Reads Fig6b.csv (relative path) and regenerates the panel: individual worms
as jittered, translucent raw points, mean +/- SD per construct per age drawn
on top, and a per-age n-count table under the x-axis, in the house style.

STATISTICS: none are drawn on this panel -- 3 constructs x 6 ages leaves no
room for brackets. The tests behind the claims made about this panel in the
Results (Mann-Whitney U with Holm-Bonferroni correction) are computed by
stats_Fig6b.py in this directory, which writes Fig6b_stats.csv.

Run:
    python plot_Fig6b.py
    python stats_Fig6b.py      # p-values for this panel

Writes Fig6b.svg and Fig6b.png next to this file.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

# --- house style ---
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from style import apply_style, save, DATA_DIR, TAU  # noqa: E402

apply_style()


VALUE_COL = "tau_Spt_norm_ns"
AGE_ORDER = ["L4", "D1", "D2", "D3", "D4", "D5"]
CONSTRUCT_ORDER = ["CeLam-WF", "CeLam-cpH", "CeLam-C-WF"]
COLORS = {"CeLam-WF": "#1f77b4", "CeLam-cpH": "#ff7f0e", "CeLam-C-WF": "#d62728"}
MARKERS = {"CeLam-WF": "o", "CeLam-cpH": "s", "CeLam-C-WF": "^"}
JITTER = {"CeLam-WF": -0.14, "CeLam-cpH": 0.0, "CeLam-C-WF": 0.14}
YLIM = (-0.80, 0.85)   # top raised: the largest raw point is 0.760 and was clipped
YTICKS = [-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75]

# Layout in inches. Canvas size is fixed to the figure slot; the margins are
# absolute so the n-count table below the axis is actually budgeted for -- it
# sits at fixed point offsets, not axes fractions, so it does not drift when the
# panel is resized.
FIG_W, FIG_H = 2.3, 2.13
LEFT_IN, RIGHT_IN, TOP_IN = 0.48, 0.03, 0.04   # left fits the y label + tick labels
BOT_IN = 0.50                                   # x tick labels + three n rows

# Row labels for the n table are the construct suffixes only ("C-WF:" not
# "CeLam-C-WF:"): at 6 pt the full names are 0.53 in wide, a quarter of the
# canvas. They stay colour-matched to the legend, which carries the full names.
ROW_LABEL = {c: c.split("-", 1)[1] for c in CONSTRUCT_ORDER}


def main():
    df = pd.read_csv(DATA_DIR / "Fig6" / "Fig6b.csv")
    df["age_day"] = pd.Categorical(df["age_day"], categories=AGE_ORDER, ordered=True)

    rng = np.random.default_rng(0)
    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    x_base = np.arange(len(AGE_ORDER))

    n_table = {}
    for construct in CONSTRUCT_ORDER:
        sub = df[df["construct"] == construct]
        grouped = sub.groupby("age_day", observed=True)[VALUE_COL]
        mean = grouped.mean().reindex(AGE_ORDER)
        sd = grouped.std().reindex(AGE_ORDER)
        n = grouped.size().reindex(AGE_ORDER)
        n_table[construct] = n

        # raw per-worm points: jittered around each age's x position
        for i, age in enumerate(AGE_ORDER):
            vals = sub.loc[sub["age_day"] == age, VALUE_COL].to_numpy()
            if vals.size == 0:
                continue
            xj = x_base[i] + JITTER[construct] + rng.uniform(-0.05, 0.05, size=vals.size)
            ax.scatter(xj, vals, marker=MARKERS[construct], s=6,
                       facecolor=COLORS[construct], edgecolor="none",
                       alpha=0.30, zorder=1, linewidths=0)

        # mean +/- SD, drawn on top, black error bars
        ax.errorbar(
            x_base, mean.to_numpy(), yerr=sd.to_numpy(),
            marker=MARKERS[construct], markersize=3.5, linewidth=0.9,
            markerfacecolor=COLORS[construct], markeredgecolor=COLORS[construct],
            color=COLORS[construct], ecolor="black", elinewidth=0.7, capsize=1.8,
            label=construct, zorder=3,
        )
        print(f"{construct}: n = {n.to_dict()}")

    ax.axhline(0, color="0.6", linewidth=0.6, linestyle="--", zorder=0)
    ax.set_xlim(-0.6, len(AGE_ORDER) - 1 + 0.6)
    ax.set_xticks(range(len(AGE_ORDER)))
    ax.set_xticklabels(AGE_ORDER)
    ax.set_ylim(*YLIM)
    ax.set_yticks(YTICKS)
    ax.set_ylabel(f"{TAU}, spermatheca (norm) (ns)")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=6, loc="upper right",
              handlelength=1.1, handletextpad=0.4, labelspacing=0.2,
              borderpad=0.1, borderaxespad=0.15)

    # n-count table under the x-axis, one row per construct. Rows are offset in
    # points from the axis edge, so they hold their spacing at any panel size.
    trans = ax.get_xaxis_transform()
    row_dy = [-15.0, -22.5, -30.0]
    for dy, construct in zip(row_dy, CONSTRUCT_ORDER):
        ax.annotate(f"{ROW_LABEL[construct]}:", xy=(0.0, 0.0),
                    xycoords="axes fraction",
                    xytext=(-2.0, dy), textcoords="offset points", ha="right",
                    va="center", fontsize=6, color=COLORS[construct],
                    annotation_clip=False)
        for i, age in enumerate(AGE_ORDER):
            n_val = n_table[construct].reindex(AGE_ORDER)[age]
            ax.annotate(f"n={n_val}", xy=(x_base[i], 0.0), xycoords=trans,
                        xytext=(0.0, dy), textcoords="offset points",
                        ha="center", va="center", fontsize=6,
                        color=COLORS[construct], annotation_clip=False)

    fig.subplots_adjust(left=LEFT_IN / FIG_W, right=1 - RIGHT_IN / FIG_W,
                        top=1 - TOP_IN / FIG_H, bottom=BOT_IN / FIG_H)
    save(fig, "Fig6b")


if __name__ == "__main__":
    main()
