"""
Fig. 6j: spermathecal local lifetime contrast, tau_Spt(norm), of CeLam-WF on
adult day 4 in worms labeled with JFX554-HTL at L4 versus on adult day 3.

tau_Spt(norm) = tau(spermatheca ROI) - mean(tau of the two adjacent-tissue ROIs),
as defined for Fig. 6b (HyD X 2 channel, intensity-weighted mean lifetime from
the LAS X per-ROI fit). One spermatheca per worm. Two independent experiments
(2026-10-19, 2026-10-21) pooled, one symbol. The L4 group includes 17 worms
that were also imaged on day 1 before the day-4 image (`also_imaged_day1`);
only their day-4 images are used. One image with a failed fit (an ROI tau <= 0)
was excluded upstream and is not in the CSV.

Reads Fig6j.csv (relative path), one row per worm; draws per-worm dots with a
black mean bar and SD whiskers (Fig. 6d/6h style) and one two-sided
Mann-Whitney U test (a single comparison, so no multiplicity correction).

Run:
    python plot_Fig6j.py

Writes Fig6j.svg and Fig6j.png.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.transforms import blended_transform_factory

# --- shared Fig. 6 statistics ---
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _stats import mwu  # noqa: E402

# --- house style ---
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from style import (apply_style, save, p_text, swarm_offsets, dot_diameter_in,
                   DOT_COLOR, DOT_ALPHA, DATA_DIR, TAU)  # noqa: E402

apply_style()


VALUE_COL = "tau_Spt_norm_ns"
ORDER = ["L4 label -> D4 FLIM", "D3 label -> D4 FLIM"]
XLABELS = ["L4 label\n→ D4 FLIM", "D3 label\n→ D4 FLIM"]

BAR_W = 0.46
XLIM = (-0.48, 1.48)
MAX_SWARM = 0.30
DOT_S = 12                      # scatter area (pt^2); house DOT_SIZE is 18

# Layout in inches, fixed to the final figure slot.
FIG_W, FIG_H = 1.5, 2.3
LEFT_IN, RIGHT_IN, TOP_IN, BOT_IN = 0.40, 0.03, 0.04, 0.42


def main():
    df = pd.read_csv(DATA_DIR / "Fig6" / "Fig6j.csv")
    vals = [df.loc[df.labeling_schedule == g, VALUE_COL].to_numpy() for g in ORDER]
    U, p, eff = mwu(vals[0], vals[1])
    for g, v in zip(ORDER, vals):
        print(f"{g}: n = {v.size}, mean = {v.mean():.3f} ns, SD = {v.std(ddof=1):.3f} ns")
    print(f"MWU U = {U:.1f}, PS = {eff:.3f}, p = {p:.3g}")

    # y range: 0 to just above the p label
    y_brk = max(v.max() for v in vals) + 0.06
    ylim = (0.0, y_brk + 0.16)
    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    d_in = dot_diameter_in(DOT_S)
    dx = d_in / ((FIG_W - LEFT_IN - RIGHT_IN) / (XLIM[1] - XLIM[0]))
    dy = d_in / ((FIG_H - TOP_IN - BOT_IN) / (ylim[1] - ylim[0]))
    for i, v in enumerate(vals):
        ax.scatter(i + swarm_offsets(v, dx, dy, MAX_SWARM), v, s=DOT_S,
                   color=DOT_COLOR, alpha=DOT_ALPHA, zorder=2, linewidths=0)
        m, sd = v.mean(), v.std(ddof=1)
        ax.errorbar(i, m, yerr=sd, fmt="none", ecolor="black", elinewidth=0.9,
                    capsize=2.2, capthick=0.9, zorder=3)
        ax.plot([i - BAR_W / 2, i + BAR_W / 2], [m, m], color="black", lw=1.2,
                solid_capstyle="butt", zorder=4)
    ax.plot([0, 0, 1, 1], [y_brk - 0.02, y_brk, y_brk, y_brk - 0.02],
            color="black", lw=0.7, zorder=5)
    p_text(ax, 0.5, y_brk, p, fontsize=7, zorder=6)

    ax.set_xlim(*XLIM)
    ax.set_ylim(*ylim)
    ax.set_yticks(np.arange(0, ylim[1] + 1e-9, 0.2))
    ax.set_ylabel(f"{TAU}, spermatheca (norm) (ns)", labelpad=1.5)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="both", labelsize=7, pad=1.5)
    ax.tick_params(axis="x", pad=11.0, length=0)
    ax.set_xticks(range(len(ORDER)))
    ax.set_xticklabels(XLABELS, fontsize=7, linespacing=1.1)
    trans = blended_transform_factory(ax.transData, ax.transAxes)
    for i, v in enumerate(vals):
        ax.annotate(f"n={len(v)}", xy=(i, 0.0), xycoords=trans, xytext=(0.0, -6.0),
                    textcoords="offset points", ha="center", va="center",
                    fontsize=7, color="0.35", annotation_clip=False)
    fig.subplots_adjust(left=LEFT_IN / FIG_W, right=1 - RIGHT_IN / FIG_W,
                        top=1 - TOP_IN / FIG_H, bottom=BOT_IN / FIG_H)
    save(fig, "Fig6j")


if __name__ == "__main__":
    main()
