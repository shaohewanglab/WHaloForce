"""
Fig. 3k: junction lifetime profiles for aCat-WF, aCat-cpH, and aCat-C-WF,
before stretch / after stretch / after LatA.

Reads data/Fig3/Fig3k_profiles.csv (one row per position per cell) by a
relative path and regenerates the three-panel line-plus-shaded-SEM-band
figure. Cells are pooled directly across all fields of view within a
construct and condition, and the shaded band is the SEM across those pooled
cells (not a per-FOV average taken first). The observation/replicate unit
for this panel is therefore the cell, pooled irrespective of which FOV it
came from.

All three condition lines are solid and are distinguished only by alpha
(Pre 1.0, Stretched 0.6, LatA 0.3). Dotted vertical guides mark the quadrant
boundaries at 90, 180 and 270 deg.

Run:
    python plot_fig3k.py

Writes Fig3k.svg and Fig3k.png to ../figures.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from style import apply_style, save, DATA_DIR, TAU  # noqa: E402

apply_style()

DATA_CSV = DATA_DIR / "Fig3" / "Fig3k_profiles.csv"

CONSTRUCT_ORDER = ["aCat-WF", "aCat-cpH", "aCat-C-WF"]
CONDITION_ORDER = ["before", "stretched", "latA"]
CONDITION_STYLES = {"before": "-", "stretched": "-", "latA": "-"}
CONDITION_ALPHA = {"before": 1.0, "stretched": 0.6, "latA": 0.3}
COLORS = {"aCat-WF": "#2c7fb8", "aCat-cpH": "#e6842e", "aCat-C-WF": "#d62728"}
YLIM = {"aCat-WF": (2.1, 2.5), "aCat-cpH": (3.6, 4.0), "aCat-C-WF": (2.1, 2.5)}


def sem(a, axis=0):
    n = np.sum(np.isfinite(a), axis=axis)
    return np.nanstd(a, axis=axis, ddof=1) / np.sqrt(np.maximum(n, 1))


def main():
    df = pd.read_csv(DATA_CSV)
    grid = np.sort(df["arc_frac"].unique())

    fig, axes = plt.subplots(1, len(CONSTRUCT_ORDER), figsize=(3.58, 1.5), sharex=True)

    for ax, construct in zip(axes, CONSTRUCT_ORDER):
        color = COLORS[construct]
        for condition in CONDITION_ORDER:
            sub = df[(df["construct"] == construct) & (df["condition"] == condition)]
            n_cells = sub["cell_id"].nunique()
            n_fov = sub["fov_id"].nunique()
            if n_cells == 0:
                continue
            # pool cells directly across FOVs
            pivot = sub.pivot(index="cell_id", columns="arc_frac", values="tau_int_ns").reindex(columns=grid)
            mat = pivot.to_numpy()
            mean = np.nanmean(mat, axis=0)
            err = sem(mat, axis=0)
            ax.plot(
                grid,
                mean,
                color=color,
                lw=0.9,
                ls=CONDITION_STYLES[condition],
                alpha=CONDITION_ALPHA[condition],
            )
            ax.fill_between(grid, mean - err, mean + err, color=color, alpha=0.18, lw=0)
            print(f"{construct} / {condition}: n = {n_cells} cells pooled from {n_fov} FOVs")

        # quadrant boundaries (90, 180, 270 deg), as in the published panel
        for xg in (0.25, 0.50, 0.75):
            ax.axvline(xg, color="0.6", lw=0.5, ls=":", zorder=0)

        ax.set_title(construct, pad=2.0)
        ax.set_xlim(0, 1)
        ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.set_xticklabels(["0", "90", "180", "270", "360"])
        ax.set_xlabel("Arc position")
        ylo, yhi = YLIM[construct]
        ax.set_ylim(ylo, yhi)
        ax.set_yticks(np.round(np.arange(ylo, yhi + 1e-9, 0.1), 2))
        if construct == CONSTRUCT_ORDER[0]:
            ax.set_ylabel(f"{TAU} (ns)")
        ax.spines[["top", "right"]].set_visible(False)
        # No in-figure legend (house style). All conditions are solid lines,
        # distinguished by alpha: Pre 1.0, Stretched 0.6, LatA 0.3 (see
        # figure caption).

    for a in axes:
        a.tick_params(axis='both', labelsize=7, pad=1.0)
    fig.subplots_adjust(left=0.115, right=0.977, bottom=0.22, top=0.90, wspace=0.30)
    save(fig, "Fig3k")


if __name__ == "__main__":
    main()
