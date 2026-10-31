"""
Fig. 3l: junctional lifetime anisotropy (dtau) for aCat-WF, aCat-cpH, and
aCat-C-WF, before stretch / after stretch / after LatA.

Reads data/Fig3/Fig3l.csv (one row per field of view) by a relative path and
regenerates the three-panel figure: each FOV is one jittered point, overlaid
with a box spanning the 25th-75th percentile (median line) and Tukey
whiskers (1.5x interquartile range). Statistics are two-sided Mann-Whitney
U for the two displayed comparisons, before vs stretched and stretched vs
stretched + LatA, with Holm-Bonferroni correction applied separately within
each construct across those two comparisons. The plot shows the Holm-adjusted
p values; raw and adjusted values are printed to the console. The undisplayed
before vs LatA comparison is not part of any Holm family; its raw p is printed
for reference only.

Run:
    python plot_fig3l.py

Writes Fig3l.svg and Fig3l.png to ../figures.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import mannwhitneyu

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from style import apply_style, save, DATA_DIR, p_text, TAU, DELTA  # noqa: E402

apply_style()

DATA_CSV = DATA_DIR / "Fig3" / "Fig3l.csv"

VALUE_COL = "dtau_anisotropy_ns"
CONSTRUCT_ORDER = ["aCat-WF", "aCat-cpH", "aCat-C-WF"]
CONDITION_ORDER = ["before", "stretched", "latA"]
CONDITION_LABELS = {"before": "Before", "stretched": "Stretched", "latA": "LatA"}
ALL_PAIRS = [("before", "stretched"), ("stretched", "latA"), ("before", "latA")]
PLOTTED_PAIRS = [("before", "stretched"), ("stretched", "latA")]
COLORS = {"aCat-WF": "#2c7fb8", "aCat-cpH": "#e6842e", "aCat-C-WF": "#d62728"}


def holm(pvals):
    p = np.asarray(pvals, float)
    order = np.argsort(p)
    adj = np.empty(p.size)
    run = 0.0
    for i, idx in enumerate(order):
        run = max(run, (p.size - i) * p[idx])
        adj[idx] = min(run, 1.0)
    return adj


def pctl_box(ax, x, vals, color, width=0.25):
    q1, med, q3 = np.percentile(vals, [25, 50, 75])
    iqr = q3 - q1
    if iqr > 0:
        low, high = q1 - 1.5 * iqr, q3 + 1.5 * iqr
        whisker_low = vals[vals >= low].min()
        whisker_high = vals[vals <= high].max()
    else:
        whisker_low, whisker_high = vals.min(), vals.max()
    ax.add_patch(mpatches.Rectangle((x - width / 2, q1), width, q3 - q1,
                 lw=0.75, edgecolor=color, facecolor=color, alpha=0.10, zorder=3))
    ax.plot([x, x], [whisker_low, q1], color=color, lw=0.75, zorder=3)
    ax.plot([x, x], [q3, whisker_high], color=color, lw=0.75, zorder=3)
    ax.plot([x - width * 0.35, x + width * 0.35], [whisker_low, whisker_low], color=color, lw=0.75, zorder=3)
    ax.plot([x - width * 0.35, x + width * 0.35], [whisker_high, whisker_high], color=color, lw=0.75, zorder=3)
    ax.plot([x - width / 2, x + width / 2], [med, med], color=color, lw=0.75, zorder=4)


def main():
    df = pd.read_csv(DATA_CSV)
    df[VALUE_COL] = df[VALUE_COL] * 1000.0  # ns -> ps, matching the figure's axis units

    # One shared y-axis across the three constructs, so the effect sizes are
    # directly comparable; bracket levels are therefore global, not per panel.
    fig, axes = plt.subplots(1, len(CONSTRUCT_ORDER), figsize=(3.58, 1.5), sharey=True)
    allv = df[VALUE_COL].to_numpy()
    g_min, g_max = float(np.nanmin(allv)), float(np.nanmax(allv))
    g_span = g_max - g_min
    LEVELS = [g_max + g_span * 0.10, g_max + g_span * 0.34]
    rng = np.random.default_rng(42)

    for ax, construct in zip(axes, CONSTRUCT_ORDER):
        sub = df[df["construct"] == construct]
        groups = {c: sub.loc[sub["condition"] == c, VALUE_COL].to_numpy() for c in CONDITION_ORDER}
        color = COLORS[construct]

        raw_p = {}
        for a, b in ALL_PAIRS:
            _, p = mannwhitneyu(groups[a], groups[b], alternative="two-sided")
            raw_p[(a, b)] = p
        # Holm family: the two displayed comparisons of this construct only
        adj_p = dict(zip(PLOTTED_PAIRS, holm([raw_p[pr] for pr in PLOTTED_PAIRS])))
        for a, b in PLOTTED_PAIRS:
            print(f"{construct}: {a} vs {b}: MWU p={raw_p[(a, b)]:.3g}, "
                  f"Holm p={adj_p[(a, b)]:.3g}")
        print(f"{construct}: before vs latA (not displayed, not in the Holm family): "
              f"MWU p={raw_p[('before', 'latA')]:.3g}")

        for i, cond in enumerate(CONDITION_ORDER):
            vals = groups[cond]
            pctl_box(ax, i, vals, color)
            x = i + rng.uniform(-0.10, 0.10, size=len(vals))
            ax.scatter(x, vals, s=6, color=color, alpha=0.62, zorder=5, linewidths=0)

        ax.axhline(0.0, color="0.6", lw=0.5, ls=":", zorder=0)

        for (a, b), y_brk in zip(PLOTTED_PAIRS, LEVELS):
            x0, x1 = CONDITION_ORDER.index(a), CONDITION_ORDER.index(b)
            # shorter bracket legs than before (0.03 -> 0.02 of the span)
            ax.plot([x0, x0, x1, x1],
                    [y_brk - g_span * 0.02, y_brk, y_brk, y_brk - g_span * 0.02],
                    color="#333333", lw=0.6)
            # stacked ("p =" over the value) so the two labels in a panel do not
            # collide at this width
            p_text(ax, (x0 + x1) / 2, y_brk, adj_p[(a, b)], fontsize=7, stacked=True)

        ax.set_xticks(range(len(CONDITION_ORDER)))
        ax.set_xticklabels([CONDITION_LABELS[c] for c in CONDITION_ORDER], rotation=30, ha="right")
        ax.set_title(construct, pad=2.0)
        if construct == CONSTRUCT_ORDER[0]:
            ax.set_ylabel(f"{TAU}(T+B) \u2212 {TAU}(L+R) (ps)")
        ax.spines[["top", "right"]].set_visible(False)

    # two tick labels only, as in the published panel
    axes[0].set_yticks([0, 100])
    # room above the upper bracket so its label clears the panel title
    axes[0].set_ylim(g_min - g_span * 0.08, LEVELS[-1] + g_span * 0.62)

    for a in axes:
        a.tick_params(axis='both', labelsize=7, pad=1.0)
    fig.subplots_adjust(left=0.155, right=0.995, bottom=0.30, top=0.865, wspace=0.18)
    save(fig, "Fig3l")


if __name__ == "__main__":
    main()
