"""Figure 2d: per-cell standard deviation of the per-FA lifetime (a spatial-
heterogeneity metric) for the three untreated vinculin constructs, log y-axis,
one dot per cell.

SD = sqrt(per-cell sample variance of per-FA tau, ddof=1). Because sqrt is
monotonic and the test is rank-based, the p-values are identical to those of the
variance version of this panel; only the axis values and units change (ns rather
than ns^2).

Reads data/Fig2/Fig2d.csv and writes figures/Fig2d.svg and .png. Box + jittered
dots on a log axis, colored per construct; each construct compared with Vin-WF by
two-sided Mann-Whitney U with Holm-Bonferroni across the two comparisons
(adjusted p printed).
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from style import (p_text, TAU, apply_style, save, DATA_DIR, VIN_COLORS,  # noqa: E402
                   BOX_ALPHA, DOT_ALPHA, DOT_EDGE)

apply_style()

ORDER = ["Vin-WF", "Vin-cpH", "Vin-N-WF"]
REF = "Vin-WF"
COLORS = [VIN_COLORS[c] for c in ORDER]


def holm(pvals):
    p = np.asarray(pvals, float)
    order = np.argsort(p)
    adj = np.empty(p.size)
    run = 0.0
    for i, idx in enumerate(order):
        run = max(run, (p.size - i) * p[idx])
        adj[idx] = min(run, 1.0)
    return adj


def vs_ref_stats(groups):
    ri = ORDER.index(REF)
    comps, praw = [], []
    for j in range(len(ORDER)):
        if j == ri:
            continue
        _, p = mannwhitneyu(groups[ri], groups[j], alternative="two-sided")
        comps.append((ri, j))
        praw.append(p)
    padj = holm(praw)
    return [(i, j, pr, pa) for (i, j), pr, pa in zip(comps, praw, padj)]


def main():
    df = pd.read_csv(DATA_DIR / "Fig2" / "Fig2d.csv")
    groups = [df.loc[df.construct == c, "tau_sd_ns"].to_numpy()
              for c in ORDER]
    n = [len(g) for g in groups]

    stats = vs_ref_stats(groups)
    for _, j, praw, padj in stats:
        print(f"[Fig2d] {ORDER[j]} vs {REF}: MWU p={praw:.3g}, Holm p={padj:.3g}")

    fig, ax = plt.subplots(figsize=(1.8, 1.8))
    positions = np.arange(1, len(ORDER) + 1)
    bp = ax.boxplot(groups, positions=positions, widths=0.55, patch_artist=True,
                    showfliers=False, zorder=2,
                    medianprops=dict(color="black", linewidth=1.1),
                    whiskerprops=dict(color="black", linewidth=0.8),
                    capprops=dict(color="black", linewidth=0.8))
    for patch, col in zip(bp["boxes"], COLORS):
        patch.set_facecolor(col); patch.set_alpha(BOX_ALPHA)
        patch.set_edgecolor(col); patch.set_linewidth(0.8)

    rng = np.random.default_rng(0)
    for pos, v, col in zip(positions, groups, COLORS):
        x = pos + (rng.random(len(v)) - 0.5) * 0.25
        ax.scatter(x, v, s=8, color=col, alpha=DOT_ALPHA, edgecolors=DOT_EDGE,
                   linewidths=0.3, zorder=3)

    allv = np.concatenate(groups)
    data_top = float(allv.max())
    ax.set_yscale("log")
    ax.set_ylim(float(allv.min()) / 1.4, data_top * 5.0)
    ax.set_xlim(0.4, len(ORDER) + 0.6)

    ax_lo, ax_hi = ax.get_ylim()
    step = 0.11 * (np.log10(ax_hi) - np.log10(ax_lo))
    stats_sorted = sorted(stats, key=lambda s: abs(positions[s[1]] - positions[s[0]]))
    for m, (i, j, praw, padj) in zip([0.5, 2.0], stats_sorted):
        x1, x2 = positions[i], positions[j]
        y = 10 ** (np.log10(data_top) + step * m)
        yt = 10 ** (np.log10(y) - 0.30 * step)
        ax.plot([x1, x1, x2, x2], [yt, y, y, yt], lw=0.8, color="black",
                clip_on=False, zorder=5)
        p_text(ax, (x1 + x2) / 2, y, padj)

    ax.set_xticks(positions)
    ax.set_xticklabels([f"{c}\n(n={k})" for c, k in zip(ORDER, n)], fontsize=7)
    ax.set_ylabel(f"Per-cell SD of per-FA {TAU} (ns)")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="x", length=0)
    fig.tight_layout(pad=0.3)
    save(fig, "Fig2d")


if __name__ == "__main__":
    main()
