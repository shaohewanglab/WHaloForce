"""
Fig4g: tissue-dependent LAM-2/laminin gamma1 lifetime (pharynx vs body-wall muscle).
Reads Fig4g.csv (tidy, one row per worm per tissue) and regenerates the
paired dot plot with a two-sided paired t-test per construct, matching
Figure 4g of Wu-2026-WHaloForce.

Run:
    python plot_fig4g.py
Writes Fig4g.svg and Fig4g.png next to this file.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.stats import ttest_rel

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from style import apply_style, save, DATA_DIR, p_text, TAU  # noqa: E402

apply_style()


VALUE_COL = "tau_int_ns"
GROUP_COL = "construct"
COND_COL = "tissue"
COND_ORDER = ["pharynx", "body_wall_muscle"]
COND_LABELS = ["Pha", "BWM"]
GROUP_ORDER = ["CeLam-WF", "CeLam-cpH", "CeLam-C-WF"]
ID_COL = "worm_id"

# y ranges and ticks read off the published panel, per construct
YLIM = {"CeLam-WF": (1.2, 4.2), "CeLam-cpH": (2.8, 5.8), "CeLam-C-WF": (1.2, 4.2)}


def paired_panel(ax, sub, title):
    wide = sub.pivot(index=ID_COL, columns=COND_COL, values=VALUE_COL)[COND_ORDER].dropna()
    x = np.arange(len(COND_ORDER))
    for _, row in wide.iterrows():
        ax.plot(x, row.values, color="0.75", linewidth=0.6, zorder=1)
    for xi, c in zip(x, COND_ORDER):
        ax.scatter(np.full(len(wide), xi), wide[c], s=2.5, color="0.35", zorder=2, linewidths=0)
    means = wide.mean()
    sds = wide.std()
    # mean drawn as a horizontal bar (marker "_"), with SD whiskers
    ax.errorbar(x, means[COND_ORDER], yerr=sds[COND_ORDER], fmt="_", color="black",
                markersize=7, markeredgewidth=1.3, capsize=2, elinewidth=0.8,
                zorder=3)
    stat, p = ttest_rel(wide[COND_ORDER[0]], wide[COND_ORDER[1]])
    ax.set_xlim(-0.34, len(COND_ORDER) - 1 + 0.34)
    ax.set_xticks(x)
    ax.set_xticklabels(COND_LABELS, fontsize=6)
    lo, hi = YLIM[title]
    ax.set_ylim(lo, hi)
    ax.set_yticks(np.round(np.linspace(lo, hi, 4), 2))
    # n in the title, p inside the axes, as in the published panel
    ax.set_title(f"{title}\n(n = {len(wide)})", fontsize=7, pad=3)
    p_text(ax, 0.03, 0.97, p, fontsize=7, ha="left", va="top",
           xycoords="axes fraction", pad_pt=0.0, zorder=6, stacked=True)
    for sp in ax.spines.values():        # boxed frame, as in the published panel
        sp.set_visible(True)
    return len(wide), stat, p


def main():
    df = pd.read_csv(DATA_DIR / "Fig4" / "Fig4g.csv")
    fig, axes = plt.subplots(1, len(GROUP_ORDER), figsize=(2.128, 1.446), sharey=False)
    results = []
    for ax, g in zip(axes, GROUP_ORDER):
        sub = df[df[GROUP_COL] == g]
        n, stat, p = paired_panel(ax, sub, g)
        results.append((g, n, stat, p))
    axes[0].set_ylabel(f"Lifetime ({TAU}) (ns)", fontsize=7, labelpad=1.5)
    for a in axes:
        # tick labels close to the axis buys the width that 7 pt needs
        a.tick_params(axis="both", labelsize=7, pad=1.0)
    # Panel width is set by the widest stacked p label (measured at 0.345 in,
    # so 0.392 in of axes); the rest of the 2.128 in canvas goes to the y label
    # on the left, a small right margin for the last title's overhang, and the
    # two gaps between frames.
    fig.subplots_adjust(left=0.145, right=0.952, bottom=0.16, top=0.80, wspace=0.518)
    save(fig, "Fig4g")
    for g, n, stat, p in results:
        print(f"{g}: n={n}, paired t = {stat:.3f}, p = {p:.4g}")


if __name__ == "__main__":
    main()
