"""Figure 1d: mean fluorescence lifetime (intensity- and amplitude-weighted)
for cpHalo-JFX554, cpWHalo-JFX554, and the JFX554-HaloTag ligand control.

Reads data/Fig1/Fig1d.csv and writes figures/Fig1d.svg and .png. One-way ANOVA
and Tukey HSD p values are printed for the legend.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from style import apply_style, save, DATA_DIR  # noqa: E402

apply_style()

ORDER = ["cpHalo-JFX554", "cpWHalo-JFX554", "JFX554-HTL"]
SHADE = {"intensity": "#444444", "amplitude": "#bbbbbb"}
OFFSET = {"intensity": -0.2, "amplitude": 0.2}
LABEL = {"intensity": "intensity weighted", "amplitude": "amplitude weighted"}


def main():
    df = pd.read_csv(DATA_DIR / "Fig1" / "Fig1d.csv")
    vals = {w: {g: df[(df.group == g) & (df.weighting == w)].tau_ns.to_numpy()
                for g in ORDER}
            for w in ("intensity", "amplitude")}

    for w in ("intensity", "amplitude"):
        groups = [vals[w][g] for g in ORDER]
        f = stats.f_oneway(*groups)
        tuk = stats.tukey_hsd(*groups)
        print(f"[{w}] ANOVA F={f.statistic:.2f}, p={f.pvalue:.2e}")
        for i, a in enumerate(ORDER):
            for j, b in enumerate(ORDER):
                if i < j:
                    print(f"    {a} vs {b}: p={tuk.pvalue[i, j]:.2e}")

    fig, ax = plt.subplots(figsize=(2.7, 2.2))
    rng = np.random.default_rng(0)
    for w in ("intensity", "amplitude"):
        for i, g in enumerate(ORDER):
            y = vals[w][g]
            xc = i + OFFSET[w]
            ax.scatter(xc + rng.uniform(-0.07, 0.07, size=len(y)), y, s=22,
                       facecolors="none", edgecolors=SHADE[w], alpha=0.6,
                       linewidths=0.8, zorder=2)
            m, sd = np.mean(y), np.std(y, ddof=1)
            ax.plot([xc - 0.14, xc + 0.14], [m, m], color="black", lw=1.0, zorder=4)
            ax.errorbar(xc, m, yerr=sd, color="black", lw=1.2, capsize=5,
                        capthick=1.2, zorder=5, fmt="none")

    allv = np.concatenate([vals[w][g] for w in vals for g in ORDER])
    ylo, yhi = np.floor(allv.min() - 0.3), np.ceil(allv.max() + 0.3)
    ax.set_xticks(range(len(ORDER)))
    ax.set_xticklabels(ORDER, rotation=30, ha="right")
    ax.set_ylabel(r"Mean $\tau$ (ns)")
    ax.set_xlim(-0.5, len(ORDER) - 0.5)
    ax.set_ylim(ylo, yhi)
    ax.set_yticks(np.arange(ylo, yhi + 0.5, 1))
    ax.spines[["top", "right"]].set_visible(False)

    handles = [Line2D([0], [0], marker="o", linestyle="none", markerfacecolor="none",
                      markeredgecolor=SHADE[w], markersize=4, markeredgewidth=0.8,
                      label=LABEL[w]) for w in ("intensity", "amplitude")]
    ax.legend(handles=handles, loc="upper right", frameon=False, fontsize=6,
              handletextpad=0.4, borderaxespad=0.2)

    fig.tight_layout(pad=0.3)
    save(fig, "Fig1d")


if __name__ == "__main__":
    main()
