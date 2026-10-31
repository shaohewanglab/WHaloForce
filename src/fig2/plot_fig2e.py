"""Figure 2e: per-FOV in-FA lifetime (pooled-NNLS) under DMSO vs blebbistatin
(20 uM), one panel per vinculin construct (Vin-WF, Vin-cpH, Vin-N-WF), one dot
per field of view (FOV); different FOVs contain different cells (earlier
versions labelled the same values per-cell). Pools
two independent experiments per construct, drawn with a single symbol: the
`dataset` column is retained in the source CSV but is no longer distinguished in
the plot (the two experiments are reported in the reporting summary instead).

Reads data/Fig2/Fig2e.csv and writes figures/Fig2e.svg and .png. Within each
construct DMSO vs Bleb is compared by two-sided Mann-Whitney U (pooled across
datasets); p values are Holm-Bonferroni adjusted across the three constructs
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
                   DOT_EDGE)

apply_style()

CONSTRUCTS = ["Vin-WF", "Vin-cpH", "Vin-N-WF"]
CONTROL, TREATMENT = "DMSO", "Bleb-20uM"
ORDER = [CONTROL, TREATMENT]
XLABEL = {"DMSO": "DMSO", "Bleb-20uM": "Bleb"}
# Colored theme of plot_dmso_bleb20_box.py: one hue per construct, with the
# vehicle drawn lighter than the treatment so the pair reads as a comparison.
BOX_ALPHA_BY_GROUP = {"DMSO": 0.25, "Bleb-20uM": 0.45}
DOT_ALPHA_BY_GROUP = {"DMSO": 0.60, "Bleb-20uM": 0.90}


def holm(pvals):
    p = np.asarray(pvals, float)
    order = np.argsort(p)
    adj = np.empty(p.size)
    run = 0.0
    for i, idx in enumerate(order):
        run = max(run, (p.size - i) * p[idx])
        adj[idx] = min(run, 1.0)
    return adj


def main():
    df = pd.read_csv(DATA_DIR / "Fig2" / "Fig2e.csv")
    data, praw = {}, []
    for c in CONSTRUCTS:
        d = df[(df.construct == c) & (df.treatment == CONTROL)]
        y = df[(df.construct == c) & (df.treatment == TREATMENT)]
        data[c] = (d.tau_inFA_ns.to_numpy(), y.tau_inFA_ns.to_numpy())
        _, p = mannwhitneyu(*data[c], alternative="two-sided")
        praw.append(p)
    padj = holm(praw)
    for c, pr, pa in zip(CONSTRUCTS, praw, padj):
        print(f"[Fig2e] {c}: DMSO vs Bleb MWU p={pr:.3g}, Holm p={pa:.3g} "
              f"(n DMSO={len(data[c][0])}, n Bleb={len(data[c][1])})")

    allv = np.concatenate([np.concatenate(data[c]) for c in CONSTRUCTS])
    lo, hi = float(allv.min()), float(allv.max())
    pad = 0.06 * (hi - lo)
    # headroom above the bracket for the p-value label
    ylo, yhi = lo - pad, hi + pad * 5.6

    positions = [1, 2]
    rng = np.random.default_rng(0)
    fig, axes = plt.subplots(1, len(CONSTRUCTS), sharey=True, figsize=(3.3, 1.8))
    for ax, c, pa in zip(axes, CONSTRUCTS, padj):
        d, y = data[c]
        groups = [d, y]
        bp = ax.boxplot(groups, positions=positions, widths=0.55,
                        patch_artist=True, showfliers=False, zorder=2,
                        medianprops=dict(color="black", linewidth=1.1),
                        whiskerprops=dict(color="black", linewidth=0.8),
                        capprops=dict(color="black", linewidth=0.8))
        hue = VIN_COLORS[c]
        for patch, g in zip(bp["boxes"], ORDER):
            patch.set_facecolor(hue); patch.set_alpha(BOX_ALPHA_BY_GROUP[g])
            patch.set_edgecolor(hue); patch.set_linewidth(0.8)
        for pos, v, g in zip(positions, groups, ORDER):
            x = pos + (rng.random(len(v)) - 0.5) * 0.25
            ax.scatter(x, v, s=8, color=hue, alpha=DOT_ALPHA_BY_GROUP[g],
                       edgecolors=DOT_EDGE, linewidths=0.3, zorder=3,
                       marker="o")

        data_top = max(float(d.max()), float(y.max()))
        step = 0.11 * (yhi - ylo)
        ybr = data_top + step
        yt = ybr - 0.30 * step
        ax.plot([1, 1, 2, 2], [yt, ybr, ybr, yt], lw=0.8, color="black",
                clip_on=False, zorder=5)
        p_text(ax, 1.5, ybr, pa)

        ax.set_xticks(positions)
        ax.set_xticklabels([f"{XLABEL[g]}\n(n={len(v)})"
                            for g, v in zip(ORDER, groups)], fontsize=7)
        ax.set_xlim(0.4, 2.6)
        ax.set_title(c, fontsize=7, color=hue)   # title colour-coded per construct
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(axis="x", length=0)

    axes[0].set_ylim(ylo, yhi)
    axes[0].set_ylabel(f"Per-FOV in-FA {TAU} (ns)")

    fig.tight_layout(pad=0.3)
    save(fig, "Fig2e")


if __name__ == "__main__":
    main()
