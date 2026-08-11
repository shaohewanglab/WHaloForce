"""
Fig. 6h: absolute spermathecal lifetime (tau_Spt) of CeLam-WF, control versus
plc-1(RNAi), in young adults with fewer than 5 ovulations per spermatheca and,
separately, in adults with more than 50 ovulations per spermatheca.

NOTE ON THE GROUP THRESHOLD: both panel titles count ovulations per
spermatheca. The high group was originally labelled "more than 100 ovulations",
which counted per animal -- a hermaphrodite has two spermathecae, so the
per-spermatheca equivalent is more than 50, and the two panels were comparing
counts on different units. It now reads "More than 50 ovulations", which
matches the per-spermatheca basis of the "fewer than 5" panel, the unit of
analysis stated in the Methods, and the Results text ("controls that had
undergone more than 50 ovulations").

The `ovulation_group` key in Fig6h.csv was renamed to match, from
`more_than_100_ovulations` to `more_than_50_ovulations` (20 rows). Upstream
code that still emits the old string -- analysis_group() in
scripts-spt-aging/plot_plc1_rnai_embryo_groups.py -- will need the same rename
before it is used to regenerate this CSV.

Reads Fig6h.csv (relative path) and regenerates both sub-panels: box-and-
scatter per group, with a two-sided Mann-Whitney U test between control and
plc-1(RNAi) within each ovulation-count group, Holm-Bonferroni corrected
across the two comparisons.

NOTE ON STATISTICS: this panel previously plotted unadjusted two-sided Welch's
t-test p-values (p = 1.5e-5 for fewer than 5 ovulations, p = 0.93 for more
than 50 ovulations, verified against
scripts-spt-aging/plot_plc1_rnai_embryo_groups.py --
condition_significance_tests(), which stores welch_p as "plotted_p" -- and the
corresponding results/plc-1-RNAi/*_tests.csv) while the legend claimed
Mann-Whitney U. The test has now been changed to match the legend and to unify
with Fig. 6b and 6d: Mann-Whitney U, Holm-corrected across the two ovulation-
group comparisons as one family. The conclusion is unchanged -- plc-1(RNAi)
above control in young adults, indistinguishable after more than 50
ovulations -- but the printed p-values differ from the earlier version of the
figure.

NOTE ON GROUPING: the high group is not a directly measured ovulation count.
In the upstream pipeline it is assigned to any sample acquired on 2025-05-21
(date-based override in analysis_group() in
scripts-spt-aging/plot_plc1_rnai_embryo_groups.py), on the premise that all
worms imaged on that date were confirmed by a separate assay to be well past
100 ovulations per animal, i.e. more than 50 per spermatheca. No per-worm
ovulation count is available for this group in the committed data -- flagged
for the authors to confirm or supply exact counts.

Run:
    python plot_Fig6h.py

Writes Fig6h.svg and Fig6h.png next to this file.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

# --- shared Fig. 6 statistics ---
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _stats import mwu_holm  # noqa: E402

# --- house style ---
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from style import (apply_style, save, p_text, swarm_offsets, dot_diameter_in,
                   DOT_SIZE, DOT_COLOR, DOT_ALPHA, DATA_DIR, TAU)  # noqa: E402

apply_style()


VALUE_COL = "tau_Spt_ns"
GROUP_ORDER = ["fewer_than_5_ovulations", "more_than_50_ovulations"]
GROUP_TITLES = ["Fewer than 5 ovulations", "More than 50 ovulations"]
CONDITION_ORDER = ["control", "plc-1_RNAi"]
CONDITION_LABELS = ["Control", "plc-1 RNAi"]

BAR_W = 0.46
XLIM = (-0.55, 1.55)
MAX_SWARM = 0.30

YLIM = (1.42, 3.20)
YTICKS = [1.5, 2.0, 2.5, 3.0]

# Layout in inches, fixed to the figure slot.
FIG_W, FIG_H = 2.6, 1.4
LEFT_IN, RIGHT_IN = 0.34, 0.02
TOP_IN, BOT_IN = 0.14, 0.26     # panel title / two-line x labels
WSPACE_IN = 0.10                # y ticks are shared, so the gap can be small


def main():
    df = pd.read_csv(DATA_DIR / "Fig6" / "Fig6h.csv")

    fig, axes = plt.subplots(1, 2, figsize=(FIG_W, FIG_H), sharey=True)

    # marker diameter in data units, for the beeswarm spacing
    sub_w = (FIG_W - LEFT_IN - RIGHT_IN - WSPACE_IN) / 2
    d_in = dot_diameter_in()
    dx = d_in / (sub_w / (XLIM[1] - XLIM[0]))
    dy = d_in / ((FIG_H - TOP_IN - BOT_IN) / (YLIM[1] - YLIM[0]))

    # Both groups' samples are needed up front: the Holm family is the two
    # control-vs-plc-1(RNAi) comparisons together, so neither adjusted p is
    # known until both raw p-values are in hand.
    per_group = []
    for group in GROUP_ORDER:
        sub = df[df["ovulation_group"] == group]
        per_group.append([sub.loc[sub["condition"] == c, VALUE_COL].to_numpy()
                          for c in CONDITION_ORDER])
    stats = mwu_holm([(v[0], v[1]) for v in per_group])

    for ax, group, title, vals, (U, p_raw, p_adj, eff) in zip(
            axes, GROUP_ORDER, GROUP_TITLES, per_group, stats):
        # Published style: per-worm dots with the mean drawn as a black bar and
        # SD whiskers on top -- no box.
        for i, v in enumerate(vals):
            ax.scatter(i + swarm_offsets(v, dx, dy, MAX_SWARM), v, s=DOT_SIZE,
                       color=DOT_COLOR, alpha=DOT_ALPHA, zorder=2, linewidths=0)
            m, sd = v.mean(), v.std(ddof=1)
            ax.errorbar(i, m, yerr=sd, fmt="none", ecolor="black", elinewidth=0.9,
                        capsize=2.2, capthick=0.9, zorder=3)
            ax.plot([i - BAR_W / 2, i + BAR_W / 2], [m, m], color="black", lw=1.2,
                    solid_capstyle="butt", zorder=4)

        ymax = max(v.max() for v in vals)
        y_brk = ymax + 0.06
        ax.plot([0, 0, 1, 1], [y_brk - 0.02, y_brk, y_brk, y_brk - 0.02],
                color="black", lw=0.7, zorder=5)
        p_text(ax, 0.5, y_brk, p_adj, fontsize=6.5, zorder=6)
        print(f"{group}: MWU U = {U:.1f}, PS = {eff:.3f}, p = {p_raw:.4g}, "
              f"Holm p = {p_adj:.4g}, n = {len(vals[0])}/{len(vals[1])}")

        ax.set_xticks(range(len(CONDITION_ORDER)))
        ax.set_xticklabels(
            [f"{lab}\n(n={len(v)})" for lab, v in zip(CONDITION_LABELS, vals)]
        )
        ax.set_title(title, fontsize=7, pad=2.5)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_xlim(*XLIM)
        ax.set_ylim(*YLIM)
        ax.set_yticks(YTICKS)
        ax.tick_params(axis="both", labelsize=7, pad=1.5)

    axes[0].set_ylabel(f"{TAU}, spermatheca (ns)", labelpad=1.5)

    plot_w = FIG_W - LEFT_IN - RIGHT_IN
    fig.subplots_adjust(left=LEFT_IN / FIG_W, right=1 - RIGHT_IN / FIG_W,
                        top=1 - TOP_IN / FIG_H, bottom=BOT_IN / FIG_H,
                        wspace=WSPACE_IN / ((plot_w - WSPACE_IN) / 2))
    save(fig, "Fig6h")


if __name__ == "__main__":
    main()
