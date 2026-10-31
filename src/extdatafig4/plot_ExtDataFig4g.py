"""Extended Data Fig. 4g: E-cadherin per-FOV SD of per-junction lifetime.

Each dot is one FOV: the within-FOV standard deviation (ns) of the per-junction
intensity-weighted mean lifetime (tau_int) across its junction segments
(segments filtered to >=2 um and >=5000 photons). Box per construct, log y.
Significance is tested on the per-FOV SDs, each non-WF construct vs the
Ecad-WF reference, Mann-Whitney U with Holm-Bonferroni correction; adjusted
p values are printed.

Reads data/ExtDataFig4/ExtDataFig4g.csv (one row per FOV) and writes figures/ExtDataFig4g.svg
and .png.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from style import apply_style, save, DATA_DIR, p_text, TAU  # noqa: E402

apply_style()

ORDER = ["Ecad-WF", "Ecad-cpH", "Ecad-C-WF"]
REF = "Ecad-WF"
COLORS = {"Ecad-WF": "#2c7fb8", "Ecad-cpH": "#e6842e", "Ecad-C-WF": "#d62728"}
CSV = "ExtDataFig4g.csv"
PROT = "Ecad"
STEM = "ExtDataFig4g"

# Axes size measured from the published panel; margins are absolute so the
# canvas can be swapped in directly.
AX_W, AX_H = 1.045, 1.292
LEFT_IN, RIGHT_IN, TOP_IN, BOT_IN = 0.46, 0.03, 0.04, 0.32
FIG_W, FIG_H = AX_W + LEFT_IN + RIGHT_IN, AX_H + TOP_IN + BOT_IN


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
    df = pd.read_csv(DATA_DIR / "ExtDataFig4" / CSV)
    vals, ns = [], []
    for c in ORDER:
        v = np.sqrt(df[df.construct == c]["tau_variance_ns2"].to_numpy())
        vals.append(v)
        ns.append(len(v))

    ref_i = ORDER.index(REF)
    comp_idx = [i for i in range(len(ORDER)) if i != ref_i]
    praw = [mannwhitneyu(vals[ref_i], vals[i], alternative="two-sided")[1]
            for i in comp_idx]
    padj = holm(praw)
    pmap = {i: padj[k] for k, i in enumerate(comp_idx)}
    print(f"[{STEM}] Holm-adjusted Mann-Whitney U vs {REF} (per-FOV SD):")
    for k, i in enumerate(comp_idx):
        print(f"    {ORDER[i]} vs {REF}: p_raw={praw[k]:.3g}, p_Holm={pmap[i]:.3g}")

    positions = np.arange(1, len(ORDER) + 1)
    rng = np.random.default_rng(0)
    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    ax.set_yscale("log")
    bp = ax.boxplot(vals, positions=positions, widths=0.55, patch_artist=True,
                    showfliers=False, zorder=2,
                    medianprops=dict(color="black", linewidth=1.1),
                    whiskerprops=dict(color="black", linewidth=0.8),
                    capprops=dict(color="black", linewidth=0.8),
                    boxprops=dict(linewidth=0.8))
    for patch, c in zip(bp["boxes"], ORDER):
        patch.set_facecolor(COLORS[c])
        patch.set_alpha(0.30)
        patch.set_edgecolor(COLORS[c])
    for pos, v, c in zip(positions, vals, ORDER):
        x = pos + (rng.random(len(v)) - 0.5) * 0.30
        ax.scatter(x, v, s=9, color=COLORS[c], edgecolors="white",
                   linewidths=0.3, zorder=4)

    allv = np.concatenate(vals)
    loglo, loghi = np.log10(allv.min()), np.log10(allv.max())
    rng_log = loghi - loglo
    step = 0.13 * rng_log
    ylo = 10 ** (loglo - 0.10 * rng_log)
    base, inc = 1.7, 2.8
    max_lybr = loghi
    for k, i in enumerate(comp_idx):
        a, b = sorted([ref_i, i])
        lybr = loghi + step * (base + inc * k)
        ybr = 10 ** lybr
        yt = 10 ** (lybr - 0.30 * step)
        ax.plot([positions[a], positions[a], positions[b], positions[b]],
                [yt, ybr, ybr, yt], lw=0.8, color="black", clip_on=False, zorder=6)
        p_text(ax, (positions[a] + positions[b]) / 2, ybr, pmap[i], zorder=7)
        max_lybr = max(max_lybr, lybr)
    yhi = 10 ** (max_lybr + 0.30 * rng_log)

    ax.set_xticks(positions)
    ax.set_xticklabels([f"{c.split('-', 1)[1]}\n(n={n})"
                        for c, n in zip(ORDER, ns)], fontsize=7)
    ax.text(0.005, -0.085, f"{PROT}:", transform=ax.transAxes, ha="right",
            va="center", fontsize=7, clip_on=False)
    ax.set_xlim(0.4, len(ORDER) + 0.6)
    # SD spans less than a decade, so the default decade ticks would label only
    # one value. Drop the lower limit to the next decade to keep two labelled
    # decades, as in Fig. 2d.
    if np.floor(np.log10(yhi)) - np.ceil(np.log10(ylo)) < 1:
        ylo = min(ylo, 10 ** np.floor(np.log10(ylo)))
    ax.set_ylim(ylo, yhi)
    ax.set_ylabel(f"Per-FOV SD of junction {TAU} (ns)")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="x", length=0)
    fig.subplots_adjust(left=LEFT_IN / FIG_W, right=1 - RIGHT_IN / FIG_W,
                        top=1 - TOP_IN / FIG_H, bottom=BOT_IN / FIG_H)
    save(fig, STEM)


if __name__ == "__main__":
    main()
