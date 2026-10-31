"""
Extended Data Fig. 5d: representative single-cell four-region lifetimes (top/bottom/left/
right), for the same stretched aCat-WF cell as Extended Data Fig. 5c.

Reads data/ExtDataFig5/ExtDataFig5d.csv by a relative path and regenerates the bar plot.
No statistics (single cell, illustrative panel).

Run:
    python plot_ExtDataFig5d.py

Writes ExtDataFig5d.svg and ExtDataFig5d.png to ../figures.
"""

from pathlib import Path

import pandas as pd
from matplotlib import pyplot as plt

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from style import apply_style, save, DATA_DIR, TAU  # noqa: E402

apply_style()

# Canvas fixed to the figure slot (the published axes measured 1.913 x 1.208 in,
# but SW asked for a larger panel). Margins are absolute, so the plot region is
# whatever is left: 2.81 x 1.73 in.
FIG_W, FIG_H = 3.2, 2.2
LEFT_IN, RIGHT_IN, TOP_IN, BOT_IN = 0.34, 0.05, 0.17, 0.30

DATA_CSV = DATA_DIR / "ExtDataFig5" / "ExtDataFig5d.csv"

REGION_ORDER = ["T", "R", "B", "L"]
REGION_LABEL = {"T": "top", "R": "right", "B": "bottom", "L": "left"}
COLOR = {"T": "#4C72B0", "R": "#DD8452", "B": "#55A868", "L": "#8172B2"}

TB_LINE_C = "#1f3864"   # T/B average: solid
LR_LINE_C = "#8B4513"   # L/R average: dashed


def main():
    df = pd.read_csv(DATA_CSV).set_index("region").loc[REGION_ORDER].reset_index()
    vals = df.set_index("region")["tau_int_ns"]

    tb_avg = vals.loc[["T", "B"]].mean()
    lr_avg = vals.loc[["L", "R"]].mean()
    dtau_ps = (tb_avg - lr_avg) * 1000

    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    bars = ax.bar(REGION_ORDER, vals.loc[REGION_ORDER],
                  color=[COLOR[r] for r in REGION_ORDER], width=0.6, zorder=2)

    for rect, r in zip(bars, REGION_ORDER):
        ax.text(rect.get_x() + rect.get_width() / 2, rect.get_height() + 0.006,
                f"{vals[r]:.3f}", ha="center", va="bottom", fontsize=6)

    ax.axhline(tb_avg, color=TB_LINE_C, ls="-", lw=0.9, zorder=1)
    ax.axhline(lr_avg, color=LR_LINE_C, ls="--", lw=0.9, zorder=1)

    ax.set_ylim(2.1, 2.7)
    ax.set_xticks(range(len(REGION_ORDER)))
    ax.set_xticklabels([REGION_LABEL[r] for r in REGION_ORDER])
    ax.set_ylabel(f"{TAU} (ns)", labelpad=1.5)
    ax.set_xlabel("Region")
    ax.set_title(f"T/B \u2212 L/R = {dtau_ps:+.1f} ps", fontsize=7, pad=3)
    for sp in ax.spines.values():      # framed, as published
        sp.set_visible(True)
    ax.tick_params(axis="both", labelsize=7, pad=1.5)

    fig.subplots_adjust(left=LEFT_IN / FIG_W, right=1 - RIGHT_IN / FIG_W,
                        top=1 - TOP_IN / FIG_H, bottom=BOT_IN / FIG_H)
    save(fig, "ExtDataFig5d")


if __name__ == "__main__":
    main()
