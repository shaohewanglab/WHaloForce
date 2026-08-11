"""Figure 3f: alpha-catenin per-FOV mean junction lifetime (intensity-weighted
mean fluorescence lifetime, tau_int, ns).

Box + per-FOV points: each point is one FOV (mean junction tau across its
segments); the within-FOV SD across segments is drawn as an error bar. The box
summarizes the per-FOV means (the biological replicate unit). Significance is
tested on the per-FOV means, each non-WF construct vs the aCat-WF reference,
Mann-Whitney U with Holm-Bonferroni correction; adjusted p values are printed.

Reads data/Fig3/Fig3f.csv (one row per junction segment; junctions already
filtered to >=2 um and >=5000 photons) and writes figures/Fig3f.svg and .png.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from style import apply_style, save, DATA_DIR, p_text, TAU  # noqa: E402

apply_style()

ORDER = ["aCat-WF", "aCat-cpH", "aCat-C-WF"]
REF = "aCat-WF"
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


def main():
    df = pd.read_csv(DATA_DIR / "Fig3" / "Fig3f.csv")
    means, sds, ns = [], [], []
    for c in ORDER:
        g = df[df.construct == c].groupby("fov_id")["tau_int_ns"]
        m = g.mean().dropna()
        means.append(m.values)
        sds.append(g.std().reindex(m.index).fillna(0.0).values)
        ns.append(len(m))

    ref_i = ORDER.index(REF)
    comp_idx = [i for i in range(len(ORDER)) if i != ref_i]
    praw = [mannwhitneyu(means[ref_i], means[i], alternative="two-sided")[1]
            for i in comp_idx]
    padj = holm(praw)
    pmap = {i: padj[k] for k, i in enumerate(comp_idx)}
    print("[Fig3f] Holm-adjusted Mann-Whitney U vs aCat-WF (per-FOV means):")
    for k, i in enumerate(comp_idx):
        print(f"    {ORDER[i]} vs {REF}: p_raw={praw[k]:.3g}, p_Holm={pmap[i]:.3g}")

    positions = np.arange(1, len(ORDER) + 1)
    rng = np.random.default_rng(0)
    fig, ax = plt.subplots(figsize=(1.9, 2.0))
    bp = ax.boxplot(means, positions=positions, widths=0.55, patch_artist=True,
                    showfliers=False, zorder=2,
                    medianprops=dict(color="black", linewidth=1.1),
                    whiskerprops=dict(color="black", linewidth=0.8),
                    capprops=dict(color="black", linewidth=0.8),
                    boxprops=dict(linewidth=0.8))
    for patch, c in zip(bp["boxes"], ORDER):
        patch.set_facecolor(COLORS[c])
        patch.set_alpha(0.30)
        patch.set_edgecolor(COLORS[c])
    for pos, m, sd, c in zip(positions, means, sds, ORDER):
        x = pos + (rng.random(len(m)) - 0.5) * 0.30
        ax.errorbar(x, m, yerr=sd, fmt="none", ecolor=COLORS[c],
                    elinewidth=0.6, capsize=1.5, capthick=0.6, alpha=0.8, zorder=3)
        ax.scatter(x, m, s=9, color=COLORS[c], edgecolors="white",
                   linewidths=0.3, zorder=4)

    allv = np.concatenate([m + s for m, s in zip(means, sds)] +
                          [m - s for m, s in zip(means, sds)])
    lo, hi = float(allv.min()), float(allv.max())
    pad = 0.06 * (hi - lo)
    ylo, yhi = lo - pad, hi + pad * 7.0
    step = 0.13 * (yhi - ylo)
    for k, i in enumerate(comp_idx):
        a, b = sorted([ref_i, i])
        ybr = hi + step * (0.8 + 1.8 * k)
        yt = ybr - 0.30 * step
        ax.plot([positions[a], positions[a], positions[b], positions[b]],
                [yt, ybr, ybr, yt], lw=0.8, color="black", clip_on=False, zorder=6)
        p_text(ax, (positions[a] + positions[b]) / 2, ybr, pmap[i], zorder=7)
        yhi = max(yhi, ybr + step)

    ax.set_xticks(positions)
    ax.set_xticklabels([f"{c.split('-', 1)[1]}\n(n={n})"
                        for c, n in zip(ORDER, ns)], fontsize=7)
    ax.text(0.005, -0.085, "aCat:", transform=ax.transAxes, ha="right",
            va="center", fontsize=7, clip_on=False)
    ax.set_xlim(0.4, len(ORDER) + 0.6)
    ax.set_ylim(ylo, yhi)
    # Integer ticks only -- the half-steps crowded a 1.9 in wide panel.
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_ylabel(f"Per-FOV mean junction {TAU} (ns)")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="x", length=0)
    # tight_layout cannot see the "aCat:" label, which is drawn outside the axes
    # in axes coordinates, so it under-reserves the left and bottom margins and
    # the panel gets clipped at the figure edge. Set the margins explicitly.
    fig.subplots_adjust(left=0.235, right=0.985, bottom=0.175, top=0.975)
    save(fig, "Fig3f")


if __name__ == "__main__":
    main()
