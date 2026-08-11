"""
Fig. 6d: absolute spermathecal lifetime (tau_Spt) of CeLam-WF, control versus
glp-1(ts), at permissive (15 C) and restrictive (25 C) growth temperature.

Reads Fig6d.csv (relative path) and regenerates the panel in the published
style: gray per-worm dots with the mean drawn as a black bar and SD whiskers on
top, a temperature header rule over each pair, and an n row under the axis. A
two-sided Mann-Whitney U test compares control with glp-1(ts) within each
temperature, Holm-Bonferroni corrected across the two comparisons.

NOTE ON STATISTICS: this panel previously plotted unadjusted two-sided Welch's
t-test p-values (p = 0.66 at 15 C, p = 1.1e-4 at 25 C, verified against
results/glp-1/glp1_20260116_15C25C_DDD1_DDD52_spt_stat_tests.csv, column
welch_t_p / "test_primary") while the legend claimed Mann-Whitney U with
Holm-Bonferroni. The test has now been changed to match the legend and to
unify with Fig. 6b and 6h: Mann-Whitney U, Holm-corrected across the two
temperature comparisons as one family. The conclusion is unchanged -- no
difference at the permissive temperature, control above glp-1(ts) at the
restrictive temperature -- but the printed p-values differ from the earlier
version of the figure.

CONSEQUENCE FOR THE METHODS TEXT: the Methods currently state that "the
permissive and restrictive temperature conditions in Fig. 6d were analyzed as
independent experiments and were not corrected against each other." That is no
longer true -- the two comparisons are now one Holm family (m = 2) -- so that
clause needs to be removed or rewritten.

Run:
    python plot_Fig6d.py

Writes Fig6d.svg and Fig6d.png next to this file.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.transforms import blended_transform_factory

# --- shared Fig. 6 statistics ---
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _stats import mwu_holm  # noqa: E402

# --- house style ---
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from style import (apply_style, save, p_text, swarm_offsets, dot_diameter_in,
                   DOT_SIZE, DOT_COLOR, DOT_ALPHA, DATA_DIR, TAU)  # noqa: E402

apply_style()


VALUE_COL = "tau_Spt_ns"
GROUP_ORDER = [(15, "control"), (15, "glp-1_ts"), (25, "control"), (25, "glp-1_ts")]
GROUP_LABELS = ["Control", "glp-1(ts)", "Control", "glp-1(ts)"]
ITALIC = [False, True, False, True]          # gene name is italic, as published

# x positions: the two temperature blocks are separated by a gap.
POS = [0.0, 1.0, 2.35, 3.35]
BLOCKS = [(0, 1, "15 °C\ngrowth", "black"), (2, 3, "25 °C\ngrowth", "#E8000B")]

# Per-worm dots use the shared Fig. 6 dot style; black mean bar with SD whiskers.
BAR_W = 0.46
XLIM = (POS[0] - 0.75, POS[-1] + 0.75)
MAX_SWARM = 0.40           # half-width of the beeswarm, well inside the group gap

YLIM = (1.58, 2.46)
YTICKS = [1.6, 1.8, 2.0, 2.2]
BRACKET_Y = 2.265          # significance bar height, above the highest dot (2.23)

# Layout in inches, fixed to the figure slot. Margins are absolute so the
# rotated x labels and the n row below them are actually budgeted for.
FIG_W, FIG_H = 1.466, 2.142
LEFT_IN, RIGHT_IN = 0.40, 0.03
TOP_IN = 0.33              # two-line temperature header + its rule
BOT_IN = 0.60              # the n row, then the 45-degree x labels below it


def main():
    df = pd.read_csv(DATA_DIR / "Fig6" / "Fig6d.csv")

    groups = [
        df.loc[(df["temperature_C"] == t) & (df["genotype"] == g), VALUE_COL].to_numpy()
        for t, g in GROUP_ORDER
    ]

    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))

    # marker diameter in data units, for the beeswarm spacing
    d_in = dot_diameter_in()
    dx = d_in / ((FIG_W - LEFT_IN - RIGHT_IN) / (XLIM[1] - XLIM[0]))
    dy = d_in / ((FIG_H - TOP_IN - BOT_IN) / (YLIM[1] - YLIM[0]))

    for i, vals in enumerate(groups):
        x, m = POS[i], vals.mean()
        sd = vals.std(ddof=1)
        # per-worm points, fanned out so large markers do not sit on top of
        # each other (the glp-1(ts) columns are the dense ones)
        ax.scatter(x + swarm_offsets(vals, dx, dy, MAX_SWARM), vals,
                   s=DOT_SIZE, color=DOT_COLOR, alpha=DOT_ALPHA, zorder=2,
                   linewidths=0)
        # mean bar with SD whiskers on top
        ax.errorbar(x, m, yerr=sd, fmt="none", ecolor="black", elinewidth=0.9,
                    capsize=2.2, capthick=0.9, zorder=3)
        ax.plot([x - BAR_W / 2, x + BAR_W / 2], [m, m], color="black", lw=1.2,
                solid_capstyle="butt", zorder=4)

    # Mann-Whitney U within each temperature, Holm-corrected across the two
    # comparisons; the adjusted p is what goes on the panel.
    stats = mwu_holm([(groups[i0], groups[i1]) for i0, i1, _, _ in BLOCKS])
    for (i0, i1, title, _), (U, p_raw, p_adj, eff) in zip(BLOCKS, stats):
        a, b = groups[i0], groups[i1]
        x0, x1 = POS[i0], POS[i1]
        drop = (YLIM[1] - YLIM[0]) * 0.022
        ax.plot([x0, x0, x1, x1],
                [BRACKET_Y - drop, BRACKET_Y, BRACKET_Y, BRACKET_Y - drop],
                color="black", lw=0.7, zorder=5)
        p_text(ax, (x0 + x1) / 2, BRACKET_Y, p_adj, fontsize=7, zorder=6)
        print(f"{title.splitlines()[0]}: MWU U = {U:.1f}, PS = {eff:.3f}, "
              f"p = {p_raw:.4g}, Holm p = {p_adj:.4g}, n = {len(a)}/{len(b)}")

    ax.set_ylim(*YLIM)
    ax.set_yticks(YTICKS)
    ax.set_xlim(*XLIM)
    ax.set_ylabel(f"{TAU}, spermatheca (ns)", labelpad=1.5)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="both", labelsize=7, pad=1.5)
    ax.tick_params(axis="x", pad=15.0)   # leaves room for the n row above

    ax.set_xticks(POS)
    ax.set_xticklabels(GROUP_LABELS, rotation=45, ha="right", fontsize=7)
    for lab, italic, (_, _, _, colr) in zip(ax.get_xticklabels(), ITALIC,
                                            [BLOCKS[0]] * 2 + [BLOCKS[1]] * 2):
        lab.set_color(colr)
        if italic:
            lab.set_fontstyle("italic")

    # temperature header above the axes: label plus a rule spanning its pair
    trans = blended_transform_factory(ax.transData, ax.transAxes)
    for i0, i1, title, colr in BLOCKS:
        mid = (POS[i0] + POS[i1]) / 2
        ax.plot([POS[i0] - 0.45, POS[i1] + 0.45], [1.045, 1.045], transform=trans,
                color=colr, lw=1.0, clip_on=False, zorder=5)
        ax.text(mid, 1.075, title, transform=trans, ha="center", va="bottom",
                fontsize=7, color=colr, linespacing=1.05)

    # n row sits just under the ticks and above the x labels, as published, at a
    # fixed point offset so it holds its place at any panel size
    for i, vals in enumerate(groups):
        ax.annotate(f"{len(vals)}", xy=(POS[i], 0.0), xycoords=trans,
                    xytext=(0.0, -8.5), textcoords="offset points",
                    ha="center", va="center", fontsize=7, color="0.35",
                    annotation_clip=False)
    ax.annotate("n:", xy=(0.0, 0.0), xycoords="axes fraction",
                xytext=(-2.0, -8.5), textcoords="offset points", ha="right",
                va="center", fontsize=7, color="0.35", annotation_clip=False)

    fig.subplots_adjust(left=LEFT_IN / FIG_W, right=1 - RIGHT_IN / FIG_W,
                        top=1 - TOP_IN / FIG_H, bottom=BOT_IN / FIG_H)
    save(fig, "Fig6d")


if __name__ == "__main__":
    main()
